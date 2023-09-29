!----------------------------------------------------------------------------------------------------
!> @brief internal microstructure state for all thermal sources and kinematics constitutive models
!----------------------------------------------------------------------------------------------------
submodule(phase) thermal

  type :: tThermalParameters
    real(pREAL) :: C_p = 0.0_pREAL                                                                  !< heat capacity
    real(pREAL), dimension(3,3) :: K = 0.0_pREAL                                                    !< thermal conductivity
    character(len=pSTRLEN), allocatable, dimension(:) :: output
  end type tThermalParameters

  integer, dimension(:), allocatable :: &
    thermal_Nsources

  type(tSourceState),  allocatable, dimension(:) :: &
    thermalState

  type :: tFieldQuantities
    real(pREAL), dimension(:), allocatable :: T, dot_T
  end type tFieldQuantities

  type(tFieldQuantities), dimension(:), allocatable :: current

  type(tThermalParameters), dimension(:), allocatable :: param

  integer :: thermal_source_maxSizeDotState


  interface

    module function source_dissipation_init(maxNsources) result(isMySource)
      integer, intent(in) :: maxNsources
      logical, dimension(:,:), allocatable :: isMySource
    end function source_dissipation_init

    module function source_externalheat_init(maxNsources) result(isMySource)
      integer, intent(in) :: maxNsources
      logical, dimension(:,:), allocatable :: isMySource
    end function source_externalheat_init


    module subroutine source_externalheat_dotState(ph, en)
      integer, intent(in) :: &
        ph, &
        en
    end subroutine source_externalheat_dotState

    module function source_dissipation_f_T(ph,en) result(f_T)
      integer, intent(in) :: &
        ph, &
        en
      real(pREAL) :: f_T
    end function source_dissipation_f_T

    module function source_externalheat_f_T(ph,en)  result(f_T)
      integer, intent(in) :: &
        ph, &
        en
      real(pREAL) :: f_T
    end function source_externalheat_f_T

 end interface

contains

!----------------------------------------------------------------------------------------------
!< @brief Initializes thermal sources and kinematics mechanism.
!----------------------------------------------------------------------------------------------
module subroutine thermal_init(phases)

  type(tDict), pointer :: &
    phases

  type(tDict), pointer :: &
    phase, &
    thermal
  type(tList), pointer :: &
    sources
  character(len=:), allocatable :: &
    refs, &
    extmsg
  integer :: &
    ph, so, &
    Nmembers


  print'(/,1x,a)', '<<<+-  phase:thermal init  -+>>>'

  allocate(current(phases%length))
  allocate(thermalState(phases%length))
  allocate(thermal_Nsources(phases%length),source = 0)
  allocate(param(phases%length))
  extmsg = ''

  do ph = 1, phases%length
    Nmembers = count(material_ID_phase == ph)
    allocate(current(ph)%T(Nmembers),source=T_ROOM)
    allocate(current(ph)%dot_T(Nmembers),source=0.0_pREAL)
    phase => phases%get_dict(ph)
    thermal => phase%get_dict('thermal',defaultVal=emptyDict)

    ! ToDo: temperature dependency of K and C_p
    if (thermal%length > 0) then
      print'(/,1x,a,i0,a)', 'phase ',ph,': '//phases%key(ph)
      refs = config_listReferences(thermal,indent=3)
      if (len(refs) > 0) print'(/,1x,a)', refs
      param(ph)%C_p = thermal%get_asReal('C_p')
      param(ph)%K(1,1) = thermal%get_asReal('K_11')
      if (any(phase_lattice(ph) == ['hP','tI'])) param(ph)%K(3,3) = thermal%get_asReal('K_33')
      param(ph)%K = crystal_symmetrize_33(param(ph)%K,phase_lattice(ph))

      ! sanity checks
      if (    param(ph)%C_p <= 0.0_pREAL )  extmsg = trim(extmsg)//' C_p'
      if (any(param(ph)%K   <  0.0_pREAL))  extmsg = trim(extmsg)//' K'
      if (    phase_rho(ph) <= 0.0_pREAL )  extmsg = trim(extmsg)//' rho'

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
      if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg))

#if defined(__GFORTRAN__)
      param(ph)%output = output_as1dStr(thermal)
#else
      param(ph)%output = thermal%get_as1dStr('output',defaultVal=emptyStrArray)
#endif
      sources => thermal%get_list('source',defaultVal=emptyList)
      thermal_Nsources(ph) = sources%length
    else
      thermal_Nsources(ph) = 0
    end if

    allocate(thermalstate(ph)%p(thermal_Nsources(ph)))

  end do

  allocate(thermal_source_type(maxval(thermal_Nsources),phases%length), source = UNDEFINED)

  if (maxval(thermal_Nsources) /= 0) then
    where(source_dissipation_init (maxval(thermal_Nsources))) thermal_source_type = THERMAL_SOURCE_DISSIPATION
    where(source_externalheat_init(maxval(thermal_Nsources))) thermal_source_type = THERMAL_SOURCE_EXTERNALHEAT
  end if

  thermal_source_maxSizeDotState = 0
  do ph = 1,phases%length

    do so = 1,thermal_Nsources(ph)
      thermalState(ph)%p(so)%state  = thermalState(ph)%p(so)%state0
    end do

    thermal_source_maxSizeDotState  = max(thermal_source_maxSizeDotState, &
                                          maxval(thermalState(ph)%p%sizeDotState))
  end do

end subroutine thermal_init


!----------------------------------------------------------------------------------------------
!< @brief Calculate thermal source (forcing term).
!----------------------------------------------------------------------------------------------
module function phase_f_T(ph,en) result(f)

  integer, intent(in) :: ph, en
  real(pREAL) :: f


  integer :: so


  f = 0.0_pREAL

  do so = 1, thermal_Nsources(ph)
   select case(thermal_source_type(so,ph))

     case (THERMAL_SOURCE_DISSIPATION)
       f = f + source_dissipation_f_T(ph,en)

     case (THERMAL_SOURCE_EXTERNALHEAT)
       f = f + source_externalheat_f_T(ph,en)

   end select

  end do

end function phase_f_T


!--------------------------------------------------------------------------------------------------
!> @brief tbd.
!--------------------------------------------------------------------------------------------------
function phase_thermal_collectDotState(ph,en) result(ok)

  integer, intent(in) :: ph, en
  logical :: ok

  integer :: i


  ok = .true.

  SourceLoop: do i = 1, thermal_Nsources(ph)

    if (thermal_source_type(i,ph) == THERMAL_SOURCE_EXTERNALHEAT) &
      call source_externalheat_dotState(ph,en)

    ok = ok .and. .not. any(IEEE_is_NaN(thermalState(ph)%p(i)%dotState(:,en)))

  end do SourceLoop

end function phase_thermal_collectDotState


!--------------------------------------------------------------------------------------------------
!> @brief Thermal viscosity.
!--------------------------------------------------------------------------------------------------
module function phase_mu_T(co,ce) result(mu)

  integer, intent(in) :: co, ce
  real(pREAL) :: mu


  mu = phase_rho(material_ID_phase(co,ce)) &
     * param(material_ID_phase(co,ce))%C_p

end function phase_mu_T


!--------------------------------------------------------------------------------------------------
!> @brief Thermal conductivity in reference configuration.
!--------------------------------------------------------------------------------------------------
module function phase_K_T(co,ce) result(K)

  integer, intent(in) :: co, ce
  real(pREAL), dimension(3,3) :: K


  K = crystallite_push33ToRef(co,ce,param(material_ID_phase(co,ce))%K)

end function phase_K_T


module function phase_thermal_constitutive(Delta_t,ph,en) result(converged_)

  real(pREAL), intent(in) :: Delta_t
  integer, intent(in) :: ph, en
  logical :: converged_


  converged_ = integrateThermalState(Delta_t,ph,en)

end function phase_thermal_constitutive


!--------------------------------------------------------------------------------------------------
!> @brief Integrate state with 1st order explicit Euler method.
!--------------------------------------------------------------------------------------------------
function integrateThermalState(Delta_t, ph,en) result(converged)

  real(pREAL), intent(in) :: Delta_t
  integer, intent(in) :: ph, en
  logical :: converged

  integer :: &
    so, &
    sizeDotState


  converged = phase_thermal_collectDotState(ph,en)
  if (converged) then

    do so = 1, thermal_Nsources(ph)
      sizeDotState = thermalState(ph)%p(so)%sizeDotState
      thermalState(ph)%p(so)%state(1:sizeDotState,en) = thermalState(ph)%p(so)%state0(1:sizeDotState,en) &
                                                      + thermalState(ph)%p(so)%dotState(1:sizeDotState,en) * Delta_t
    end do

  end if

end function integrateThermalState


module subroutine thermal_restartWrite(groupHandle,ph)

  integer(HID_T), intent(in) :: groupHandle
  integer, intent(in) :: ph

  integer :: so


  do so = 1,thermal_Nsources(ph)
    call HDF5_write(thermalState(ph)%p(so)%state,groupHandle,'omega_thermal')
  end do

end subroutine thermal_restartWrite


module subroutine thermal_restartRead(groupHandle,ph)

  integer(HID_T), intent(in) :: groupHandle
  integer, intent(in) :: ph

  integer :: so


  do so = 1,thermal_Nsources(ph)
    call HDF5_read(thermalState(ph)%p(so)%state0,groupHandle,'omega_thermal')
  end do

end subroutine thermal_restartRead


module subroutine thermal_forward()

  integer :: ph, so


  do ph = 1, size(thermalState)
    do so = 1, size(thermalState(ph)%p)
      thermalState(ph)%p(so)%state0 = thermalState(ph)%p(so)%state
    end do
  end do

end subroutine thermal_forward


!----------------------------------------------------------------------------------------------
!< @brief Get temperature (for use by non-thermal physics).
!----------------------------------------------------------------------------------------------
pure module function thermal_T(ph,en) result(T)

  integer, intent(in) :: ph, en
  real(pREAL) :: T


  T = current(ph)%T(en)

end function thermal_T


!----------------------------------------------------------------------------------------------
!< @brief Get rate of temperature (for use by non-thermal physics).
!----------------------------------------------------------------------------------------------
module function thermal_dot_T(ph,en) result(dot_T)

  integer, intent(in) :: ph, en
  real(pREAL) :: dot_T


  dot_T = current(ph)%dot_T(en)

end function thermal_dot_T


!----------------------------------------------------------------------------------------------
!< @brief Set temperature
!----------------------------------------------------------------------------------------------
module subroutine phase_thermal_setField(T,dot_T, co,ce)

  real(pREAL), intent(in) :: T, dot_T
  integer, intent(in) :: ce, co


  current(material_ID_phase(co,ce))%T(material_entry_phase(co,ce)) = T
  current(material_ID_phase(co,ce))%dot_T(material_entry_phase(co,ce)) = dot_T

end subroutine phase_thermal_setField



!--------------------------------------------------------------------------------------------------
!> @brief checks if a source mechanism is active or not
!--------------------------------------------------------------------------------------------------
function thermal_active(source_label,src_length)  result(active_source)

  character(len=*), intent(in)         :: source_label                                              !< name of source mechanism
  integer,          intent(in)         :: src_length                                                !< max. number of sources in system
  logical, dimension(:,:), allocatable :: active_source

  type(tDict), pointer :: &
    phases, &
    phase, &
    thermal, &
    src
  type(tList), pointer :: &
    sources
  integer :: p,s

  phases => config_material%get_dict('phase')
  allocate(active_source(src_length,phases%length), source = .false. )
  do p = 1, phases%length
    phase => phases%get_dict(p)
    thermal => phase%get_dict('thermal',defaultVal=emptyDict)
    sources => thermal%get_list('source',defaultVal=emptyList)
    do s = 1, sources%length
      src => sources%get_dict(s)
      active_source(s,p) = src%get_asStr('type') == source_label
    end do
  end do


end function thermal_active


!----------------------------------------------------------------------------------------------
!< @brief Write thermal sources results to HDF5 output file.
!----------------------------------------------------------------------------------------------
module subroutine thermal_result(group,ph)

  character(len=*), intent(in) :: group
  integer,          intent(in) :: ph


  integer :: ou

  if (.not. allocated(param(ph)%output)) return

  call result_closeGroup(result_addGroup(group//'thermal'))

  do ou = 1, size(param(ph)%output)

    select case(trim(param(ph)%output(ou)))

      case ('T')
        call result_writeDataset(current(ph)%T,group//'thermal','T', 'temperature','K')

    end select

  end do

end subroutine thermal_result


end submodule thermal
