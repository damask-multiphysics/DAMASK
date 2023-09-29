!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic brittle damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule (phase:damage) anisobrittle

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pREAL) :: &
      dot_o_0, &                                                                                    !< opening rate of cleavage planes
      p                                                                                             !< damage rate sensitivity
    real(pREAL), dimension(:), allocatable :: &
      s_crit, &                                                                                     !< critical displacement
      g_crit                                                                                        !< critical load
    real(pREAL), dimension(:,:,:,:), allocatable :: &
      cleavage_systems
    integer :: &
      sum_N_cl                                                                                      !< total number of cleavage planes
    character(len=pSTRLEN), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstances)


contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function anisobrittle_init() result(mySources)

  logical, dimension(:), allocatable :: mySources

  type(tDict), pointer :: &
    phases, &
    phase, &
    src
  integer :: Nmembers,ph
  integer, dimension(:), allocatable :: N_cl
  character(len=:), allocatable :: &
    refs, &
    extmsg


  mySources = source_active('anisobrittle')
  if (count(mySources) == 0) return

  print'(/,1x,a)', '<<<+-  phase:damage:anisobrittle init  -+>>>'
  print'(/,1x,a,1x,i0)', '# phases:',count(mySources); flush(IO_STDOUT)


  phases => config_material%get_dict('phase')
  allocate(param(phases%length))
  extmsg = ''

  do ph = 1, phases%length
    if (mySources(ph)) then
      phase => phases%get_dict(ph)
      src => phase%get_dict('damage')

      associate(prm  => param(ph))

        print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)
        refs = config_listReferences(src,indent=3)
        if (len(refs) > 0) print'(/,1x,a)', refs

        N_cl = src%get_as1dInt('N_cl',defaultVal=emptyIntArray)
        prm%sum_N_cl = sum(abs(N_cl))

        prm%p       = src%get_asReal('p')
        prm%dot_o_0 = src%get_asReal('dot_o_0')

        prm%s_crit  = src%get_as1dReal('s_crit',requiredSize=size(N_cl))
        prm%g_crit  = src%get_as1dReal('g_crit',requiredSize=size(N_cl))

        prm%cleavage_systems = crystal_SchmidMatrix_cleavage(N_cl,phase_lattice(ph),phase_cOverA(ph))

        ! expand: family => system
        prm%s_crit = math_expand(prm%s_crit,N_cl)
        prm%g_crit = math_expand(prm%g_crit,N_cl)

#if defined (__GFORTRAN__)
        prm%output = output_as1dStr(src)
#else
        prm%output = src%get_as1dStr('output',defaultVal=emptyStrArray)
#endif

        ! sanity checks
        if (prm%p          <= 0.0_pREAL)  extmsg = trim(extmsg)//' p'
        if (prm%dot_o_0    <= 0.0_pREAL)  extmsg = trim(extmsg)//' dot_o_0'
        if (any(prm%g_crit <  0.0_pREAL)) extmsg = trim(extmsg)//' g_crit'
        if (any(prm%s_crit <  0.0_pREAL)) extmsg = trim(extmsg)//' s_crit'

        Nmembers = count(material_ID_phase==ph)
        call phase_allocateState(damageState(ph),Nmembers,1,1,0)
        damageState(ph)%atol = src%get_asReal('atol_phi',defaultVal=1.0e-9_pREAL)
        if (any(damageState(ph)%atol < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_phi'

      end associate

      if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(damage_anisoBrittle)')
    end if

  end do

end function anisobrittle_init


!--------------------------------------------------------------------------------------------------
!> @brief
!--------------------------------------------------------------------------------------------------
module subroutine anisobrittle_dotState(M_i, ph,en)

  integer, intent(in) :: &
    ph,en
  real(pREAL),  intent(in), dimension(3,3) :: &
    M_i

  integer :: &
    a, i
  real(pREAL) :: &
    traction, traction_crit


  associate(prm => param(ph))
    damageState(ph)%dotState(1,en) = 0.0_pREAL
    do a = 1, prm%sum_N_cl
      traction_crit = damage_phi(ph,en)**2 * prm%g_crit(a)
      do i = 1,3
        traction = math_tensordot(M_i,prm%cleavage_systems(1:3,1:3,i,a))

        damageState(ph)%dotState(1,en) = damageState(ph)%dotState(1,en) &
          + prm%dot_o_0 / prm%s_crit(a) &
            * (max(0.0_pREAL, abs(traction) - traction_crit)/traction_crit)**prm%p
      end do
    end do
  end associate

end subroutine anisobrittle_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine anisobrittle_result(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o


  associate(prm => param(phase), stt => damageState(phase)%state)
    outputsLoop: do o = 1,size(prm%output)
      select case(trim(prm%output(o)))
        case ('Psi_D')
          call result_writeDataset(stt,group,trim(prm%output(o)),'damage energy density','J/m³')
      end select
    end do outputsLoop
  end associate

end subroutine anisobrittle_result


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
!--------------------------------------------------------------------------------------------------
module subroutine damage_anisobrittle_LiAndItsTangent(L_i, dL_i_dM_i, M_i, ph,en)

  integer, intent(in) :: &
    ph,en
  real(pREAL),   intent(in),  dimension(3,3) :: &
    M_i
  real(pREAL),   intent(out), dimension(3,3) :: &
    L_i                                                                                             !< damage velocity gradient
  real(pREAL),   intent(out), dimension(3,3,3,3) :: &
    dL_i_dM_i                                                                                       !< derivative of L_i with respect to M_i

  integer :: &
    a, k, l, m, n, i
  real(pREAL) :: &
    traction, traction_crit, &
    udot, dudot_dt


  L_i = 0.0_pREAL
  dL_i_dM_i = 0.0_pREAL
  associate(prm => param(ph))
    do a = 1,prm%sum_N_cl
      traction_crit = damage_phi(ph,en)**2 * prm%g_crit(a)

      do i = 1, 3
        traction = math_tensordot(M_i,prm%cleavage_systems(1:3,1:3,i,a))
        if (abs(traction) > traction_crit + tol_math_check) then
          udot = sign(1.0_pREAL,traction)* prm%dot_o_0 * ((abs(traction) - traction_crit)/traction_crit)**prm%p
          L_i = L_i + udot*prm%cleavage_systems(1:3,1:3,i,a)
          dudot_dt = sign(1.0_pREAL,traction)*udot*prm%p / (abs(traction) - traction_crit)
          forall (k=1:3,l=1:3,m=1:3,n=1:3) &
            dL_i_dM_i(k,l,m,n) = dL_i_dM_i(k,l,m,n) &
                               + dudot_dt*prm%cleavage_systems(k,l,i,a) * prm%cleavage_systems(m,n,i,a)
        end if
      end do
    end do
  end associate

end subroutine damage_anisobrittle_LiAndItsTangent

end submodule anisobrittle
