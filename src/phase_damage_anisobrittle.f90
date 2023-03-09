!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic brittle damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule (phase:damage) anisobrittle

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
  print'(/,a,i0)', ' # phases: ',count(mySources); flush(IO_STDOUT)


  phases => config_material%get_dict('phase')
  allocate(param(phases%length))
  extmsg = ''

  do ph = 1, phases%length
    if (mySources(ph)) then
      phase => phases%get_dict(ph)
      src => phase%get_dict('damage')

      associate(prm  => param(ph))

        print'(/,1x,a,i0,a)', 'phase ',ph,': '//phases%key(ph)
        refs = config_listReferences(src,indent=3)
        if (len(refs) > 0) print'(/,1x,a)', refs

        N_cl = src%get_as1dInt('N_cl',defaultVal=emptyIntArray)
        prm%sum_N_cl = sum(abs(N_cl))

        prm%q       = src%get_asFloat('q')
        prm%dot_o   = src%get_asFloat('dot_o')

        prm%s_crit  = src%get_as1dFloat('s_crit',  requiredSize=size(N_cl))
        prm%g_crit  = src%get_as1dFloat('g_crit',  requiredSize=size(N_cl))

        prm%cleavage_systems = lattice_SchmidMatrix_cleavage(N_cl,phase_lattice(ph),phase_cOverA(ph))

        ! expand: family => system
        prm%s_crit = math_expand(prm%s_crit,N_cl)
        prm%g_crit = math_expand(prm%g_crit,N_cl)

#if defined (__GFORTRAN__)
        prm%output = output_as1dString(src)
#else
        prm%output = src%get_as1dString('output',defaultVal=emptyStringArray)
#endif

        ! sanity checks
        if (prm%q          <= 0.0_pReal)  extmsg = trim(extmsg)//' q'
        if (prm%dot_o      <= 0.0_pReal)  extmsg = trim(extmsg)//' dot_o'
        if (any(prm%g_crit <  0.0_pReal)) extmsg = trim(extmsg)//' g_crit'
        if (any(prm%s_crit <  0.0_pReal)) extmsg = trim(extmsg)//' s_crit'

        Nmembers = count(material_ID_phase==ph)
        call phase_allocateState(damageState(ph),Nmembers,1,1,0)
        damageState(ph)%atol = src%get_asFloat('atol_phi',defaultVal=1.0e-9_pReal)
        if (any(damageState(ph)%atol < 0.0_pReal)) extmsg = trim(extmsg)//' atol_phi'

      end associate

      if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(damage_anisoBrittle)')
    end if

  end do

end function anisobrittle_init


!--------------------------------------------------------------------------------------------------
!> @brief
!--------------------------------------------------------------------------------------------------
module subroutine anisobrittle_dotState(S, ph,en)

  integer, intent(in) :: &
    ph,en
  real(pReal),  intent(in), dimension(3,3) :: &
    S

  integer :: &
    a
  real(pReal) :: &
    traction_d, traction_t, traction_n, traction_crit


  associate(prm => param(ph))
    damageState(ph)%dotState(1,en) = 0.0_pReal
    do a = 1, prm%sum_N_cl
      traction_d = math_tensordot(S,prm%cleavage_systems(1:3,1:3,1,a))
      traction_t = math_tensordot(S,prm%cleavage_systems(1:3,1:3,2,a))
      traction_n = math_tensordot(S,prm%cleavage_systems(1:3,1:3,3,a))

      traction_crit = prm%g_crit(a)*damage_phi(ph,en)**2

      damageState(ph)%dotState(1,en) = damageState(ph)%dotState(1,en) &
          + prm%dot_o / prm%s_crit(a) &
            * ((max(0.0_pReal, abs(traction_d) - traction_crit)/traction_crit)**prm%q + &
               (max(0.0_pReal, abs(traction_t) - traction_crit)/traction_crit)**prm%q + &
               (max(0.0_pReal, abs(traction_n) - traction_crit)/traction_crit)**prm%q)
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
        case ('f_phi')
          call result_writeDataset(stt,group,trim(prm%output(o)),'driving force','-')
      end select
    end do outputsLoop
  end associate

end subroutine anisobrittle_result


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
!--------------------------------------------------------------------------------------------------
module subroutine damage_anisobrittle_LiAndItsTangent(Ld, dLd_dTstar, S, ph,en)

  integer, intent(in) :: &
    ph,en
  real(pReal),   intent(in),  dimension(3,3) :: &
    S
  real(pReal),   intent(out), dimension(3,3) :: &
    Ld                                                                                              !< damage velocity gradient
  real(pReal),   intent(out), dimension(3,3,3,3) :: &
    dLd_dTstar                                                                                      !< derivative of Ld with respect to Tstar (4th-order tensor)

  integer :: &
    a, k, l, m, n
  real(pReal) :: &
    traction, traction_crit, &
    udot, dudot_dt


  Ld = 0.0_pReal
  dLd_dTstar = 0.0_pReal
  associate(prm => param(ph))
    do a = 1,prm%sum_N_cl
      traction_crit = prm%g_crit(a)*damage_phi(ph,en)**2

      traction = math_tensordot(S,prm%cleavage_systems(1:3,1:3,1,a))
      if (abs(traction) > traction_crit + tol_math_check) then
        udot = sign(1.0_pReal,traction)* prm%dot_o * ((abs(traction) - traction_crit)/traction_crit)**prm%q
        Ld = Ld + udot*prm%cleavage_systems(1:3,1:3,1,a)
        dudot_dt = sign(1.0_pReal,traction)*udot*prm%q / (abs(traction) - traction_crit)
        forall (k=1:3,l=1:3,m=1:3,n=1:3) &
          dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) &
                              + dudot_dt*prm%cleavage_systems(k,l,1,a) * prm%cleavage_systems(m,n,1,a)
      end if

      traction = math_tensordot(S,prm%cleavage_systems(1:3,1:3,2,a))
      if (abs(traction) > traction_crit + tol_math_check) then
        udot = sign(1.0_pReal,traction)* prm%dot_o * ((abs(traction) - traction_crit)/traction_crit)**prm%q
        Ld = Ld + udot*prm%cleavage_systems(1:3,1:3,2,a)
        dudot_dt = sign(1.0_pReal,traction)*udot*prm%q / (abs(traction) - traction_crit)
        forall (k=1:3,l=1:3,m=1:3,n=1:3) &
          dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) &
                              + dudot_dt*prm%cleavage_systems(k,l,2,a) * prm%cleavage_systems(m,n,2,a)
      end if

      traction = math_tensordot(S,prm%cleavage_systems(1:3,1:3,3,a))
      if (abs(traction) > traction_crit + tol_math_check) then
        udot = sign(1.0_pReal,traction)* prm%dot_o * ((abs(traction) - traction_crit)/traction_crit)**prm%q
        Ld = Ld + udot*prm%cleavage_systems(1:3,1:3,3,a)
        dudot_dt = sign(1.0_pReal,traction)*udot*prm%q / (abs(traction) - traction_crit)
        forall (k=1:3,l=1:3,m=1:3,n=1:3) &
          dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) &
                              + dudot_dt*prm%cleavage_systems(k,l,3,a) * prm%cleavage_systems(m,n,3,a)
      end if
    end do
  end associate

end subroutine damage_anisobrittle_LiAndItsTangent

end submodule anisobrittle
