!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic brittle damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule (phase:damagee) anisobrittle

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pReal) :: &
      dot_o, &                                                                                      !< opening rate of cleavage planes
      q                                                                                             !< damage rate sensitivity
    real(pReal), dimension(:), allocatable :: &
      s_crit, &                                                                                     !< critical displacement
      g_crit                                                                                        !< critical load
    real(pReal), dimension(:,:,:,:), allocatable :: &
      cleavage_systems
    integer :: &
      sum_N_cl                                                                                      !< total number of cleavage planes
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstances)


contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function anisobrittle_init(source_length) result(mySources)

  integer, intent(in)                  :: source_length
  logical, dimension(:,:), allocatable :: mySources

  class(tNode), pointer :: &
    phases, &
    phase, &
    sources, &
    src
  integer :: Ninstances,sourceOffset,Nconstituents,p
  integer, dimension(:), allocatable :: N_cl
  character(len=pStringLen) :: extmsg = ''

  print'(/,a)', ' <<<+-  phase:damage:anisobrittle init  -+>>>'

  mySources = source_active('anisobrittle',source_length)
  Ninstances = count(mySources)
  print'(a,i2)', ' # instances: ',Ninstances; flush(IO_STDOUT)
  if(Ninstances == 0) return

  phases => config_material%get('phase')
  allocate(param(phases%length))


  do p = 1, phases%length
    phase => phases%get(p)
    if(count(mySources(:,p)) == 0) cycle
    sources => phase%get('damage')
    do sourceOffset = 1, sources%length
      if(mySources(sourceOffset,p)) then
        associate(prm  => param(p))
        src => sources%get(sourceOffset)

        N_cl = src%get_asInts('N_cl',defaultVal=emptyIntArray)
        prm%sum_N_cl = sum(abs(N_cl))

        prm%q       = src%get_asFloat('q')
        prm%dot_o   = src%get_asFloat('dot_o')

        prm%s_crit  = src%get_asFloats('s_crit',  requiredSize=size(N_cl))
        prm%g_crit  = src%get_asFloats('g_crit',  requiredSize=size(N_cl))

        prm%cleavage_systems = lattice_SchmidMatrix_cleavage(N_cl,phase%get_asString('lattice'),&
                                                             phase%get_asFloat('c/a',defaultVal=0.0_pReal))

        ! expand: family => system
        prm%s_crit = math_expand(prm%s_crit,N_cl)
        prm%g_crit = math_expand(prm%g_crit,N_cl)

#if defined (__GFORTRAN__)
        prm%output = output_asStrings(src)
#else
        prm%output = src%get_asStrings('output',defaultVal=emptyStringArray)
#endif

          ! sanity checks
        if (prm%q          <= 0.0_pReal)  extmsg = trim(extmsg)//' q'
        if (prm%dot_o      <= 0.0_pReal)  extmsg = trim(extmsg)//' dot_o'
        if (any(prm%g_crit <  0.0_pReal)) extmsg = trim(extmsg)//' g_crit'
        if (any(prm%s_crit <  0.0_pReal)) extmsg = trim(extmsg)//' s_crit'

        Nconstituents = count(material_phaseAt==p) * discretization_nIPs
        call phase_allocateState(damageState(p),Nconstituents,1,1,0)
        damageState(p)%atol = src%get_asFloat('anisobrittle_atol',defaultVal=1.0e-3_pReal)
        if(any(damageState(p)%atol < 0.0_pReal)) extmsg = trim(extmsg)//' anisobrittle_atol'

        end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
        if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(damage_anisoBrittle)')
      endif
    enddo
  enddo

end function anisobrittle_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
module subroutine anisobrittle_dotState(S, ph,me)

  integer, intent(in) :: &
    ph,me
  real(pReal),  intent(in), dimension(3,3) :: &
    S

  integer :: &
    sourceOffset, &
    damageOffset, &
    homog, &
    i
  real(pReal) :: &
    traction_d, traction_t, traction_n, traction_crit


  associate(prm => param(ph))
    damageState(ph)%dotState(1,me) = 0.0_pReal
    do i = 1, prm%sum_N_cl
      traction_d = math_tensordot(S,prm%cleavage_systems(1:3,1:3,1,i))
      traction_t = math_tensordot(S,prm%cleavage_systems(1:3,1:3,2,i))
      traction_n = math_tensordot(S,prm%cleavage_systems(1:3,1:3,3,i))

      traction_crit = prm%g_crit(i)*damage_phi(ph,me)**2.0_pReal

      damageState(ph)%dotState(1,me) = damageState(ph)%dotState(1,me) &
          + prm%dot_o / prm%s_crit(i) &
            * ((max(0.0_pReal, abs(traction_d) - traction_crit)/traction_crit)**prm%q + &
               (max(0.0_pReal, abs(traction_t) - traction_crit)/traction_crit)**prm%q + &
               (max(0.0_pReal, abs(traction_n) - traction_crit)/traction_crit)**prm%q)
    enddo
  end associate

end subroutine anisobrittle_dotState


!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
module subroutine anisobrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, ph, me)

  integer, intent(in) :: &
    ph, &
    me
  real(pReal),  intent(in) :: &
    phi
  real(pReal),  intent(out) :: &
    localphiDot, &
    dLocalphiDot_dPhi


  dLocalphiDot_dPhi = -damageState(ph)%state(1,me)

  localphiDot = 1.0_pReal &
              + dLocalphiDot_dPhi*phi

end subroutine anisobrittle_getRateAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
module subroutine anisobrittle_results(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o


  associate(prm => param(phase), stt => damageState(phase)%state)
    outputsLoop: do o = 1,size(prm%output)
      select case(trim(prm%output(o)))
        case ('f_phi')
          call results_writeDataset(group,stt,trim(prm%output(o)),'driving force','J/m³')
      end select
    enddo outputsLoop
  end associate

end subroutine anisobrittle_results

end submodule anisobrittle
