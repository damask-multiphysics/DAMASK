!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from opening of slip planes
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:eigen) slipplaneopening

  integer, dimension(:), allocatable :: damage_isoductile_instance

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    integer :: &
      sum_N_sl                                                                                      !< total number of cleavage planes
    real(pReal) :: &
      dot_o, &                                                                                      !< opening rate of cleavage planes
      q                                                                                             !< damage rate sensitivity
    real(pReal), dimension(:),   allocatable :: &
      g_crit
    real(pReal), dimension(:,:,:), allocatable     :: &
      P_d, &
      P_t, &
      P_n
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstances)


contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function damage_isoductile_init() result(myKinematics)

  logical, dimension(:), allocatable :: myKinematics

  integer :: p,i
  character(len=pStringLen) :: extmsg = ''
  integer,     dimension(:),   allocatable :: N_sl
  real(pReal), dimension(:,:), allocatable :: d,n,t
  class(tNode), pointer :: &
    phases, &
    phase, &
    mech, &
    pl, &
    kinematics, &
    kinematic_type


  myKinematics = kinematics_active2('isoductile')
  if(count(myKinematics) == 0) return
  print'(/,a)', ' <<<+-  phase:mechanical:eigen:slipplaneopening init  -+>>>'
  print'(a,i2)', ' # phases: ',count(myKinematics); flush(IO_STDOUT)


  phases => config_material%get('phase')
  allocate(param(phases%length))

  do p = 1, phases%length
    if(myKinematics(p)) then
      phase => phases%get(p)
      mech  => phase%get('mechanical')
      pl    => mech%get('plastic')

      kinematics => phase%get('damage')

      associate(prm  => param(p))
        kinematic_type => kinematics%get(1)

        prm%dot_o    = kinematic_type%get_asFloat('dot_o')
        prm%q        = kinematic_type%get_asFloat('q')
        N_sl         = pl%get_as1dInt('N_sl')
        prm%sum_N_sl = sum(abs(N_sl))

        d = lattice_slip_direction (N_sl,phase%get_asString('lattice'),&
                                    phase%get_asFloat('c/a',defaultVal=0.0_pReal))
        t = lattice_slip_transverse(N_sl,phase%get_asString('lattice'),&
                                    phase%get_asFloat('c/a',defaultVal=0.0_pReal))
        n = lattice_slip_normal    (N_sl,phase%get_asString('lattice'),&
                                    phase%get_asFloat('c/a',defaultVal=0.0_pReal))
        allocate(prm%P_d(3,3,size(d,2)),prm%P_t(3,3,size(t,2)),prm%P_n(3,3,size(n,2)))

        do i=1, size(n,2)
          prm%P_d(1:3,1:3,i) = math_outer(d(1:3,i), n(1:3,i))
          prm%P_t(1:3,1:3,i) = math_outer(t(1:3,i), n(1:3,i))
          prm%P_n(1:3,1:3,i) = math_outer(n(1:3,i), n(1:3,i))
        enddo

        prm%g_crit = kinematic_type%get_as1dFloat('g_crit',requiredSize=size(N_sl))

        ! expand: family => system
        prm%g_crit = math_expand(prm%g_crit,N_sl)

        ! sanity checks
        if (prm%q          <= 0.0_pReal)  extmsg = trim(extmsg)//' anisoDuctile_n'
        if (prm%dot_o      <= 0.0_pReal)  extmsg = trim(extmsg)//' anisoDuctile_sdot0'
        if (any(prm%g_crit <  0.0_pReal)) extmsg = trim(extmsg)//' anisoDuctile_critLoad'

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
        if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(slipplane_opening)')

      end associate
    endif
  enddo


end function damage_isoductile_init


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
!--------------------------------------------------------------------------------------------------
module subroutine damage_isoductile_LiAndItsTangent(Ld, dLd_dTstar, S, ph,me)

  integer, intent(in) :: &
    ph, me
  real(pReal),   intent(in),  dimension(3,3) :: &
    S
  real(pReal),   intent(out), dimension(3,3) :: &
    Ld                                                                                              !< damage velocity gradient
  real(pReal),   intent(out), dimension(3,3,3,3) :: &
    dLd_dTstar                                                                                      !< derivative of Ld with respect to Tstar (4th-order tensor)

  integer :: &
    i, k, l, m, n
  real(pReal) :: &
    traction_d, traction_t, traction_n, traction_crit, &
    udotd, dudotd_dt, udott, dudott_dt, udotn, dudotn_dt


  associate(prm => param(ph))
  Ld = 0.0_pReal
  dLd_dTstar = 0.0_pReal
  do i = 1, prm%sum_N_sl

    traction_d = math_tensordot(S,prm%P_d(1:3,1:3,i))
    traction_t = math_tensordot(S,prm%P_t(1:3,1:3,i))
    traction_n = math_tensordot(S,prm%P_n(1:3,1:3,i))

    traction_crit = prm%g_crit(i)* damage_phi(ph,me)

    udotd = sign(1.0_pReal,traction_d)* prm%dot_o* (  abs(traction_d)/traction_crit &
                                                    - abs(traction_d)/prm%g_crit(i))**prm%q
    udott = sign(1.0_pReal,traction_t)* prm%dot_o* (  abs(traction_t)/traction_crit &
                                                    - abs(traction_t)/prm%g_crit(i))**prm%q
    udotn = prm%dot_o* (  max(0.0_pReal,traction_n)/traction_crit &
                        - max(0.0_pReal,traction_n)/prm%g_crit(i))**prm%q

    if (dNeq0(traction_d)) then
      dudotd_dt = udotd*prm%q/traction_d
    else
      dudotd_dt = 0.0_pReal
    endif
    if (dNeq0(traction_t)) then
      dudott_dt = udott*prm%q/traction_t
    else
      dudott_dt = 0.0_pReal
    endif
    if (dNeq0(traction_n)) then
      dudotn_dt = udotn*prm%q/traction_n
    else
      dudotn_dt = 0.0_pReal
    endif

    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) &
                          + dudotd_dt*prm%P_d(k,l,i)*prm%P_d(m,n,i) &
                          + dudott_dt*prm%P_t(k,l,i)*prm%P_t(m,n,i) &
                          + dudotn_dt*prm%P_n(k,l,i)*prm%P_n(m,n,i)

    Ld = Ld &
       + udotd*prm%P_d(1:3,1:3,i) &
       + udott*prm%P_t(1:3,1:3,i) &
       + udotn*prm%P_n(1:3,1:3,i)
  enddo

  end associate

end subroutine damage_isoductile_LiAndItsTangent

end submodule slipplaneopening
