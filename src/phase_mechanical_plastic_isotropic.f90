!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for isotropic plasticity
!> @details Isotropic Plasticity which resembles the phenopowerlaw plasticity without
!! resolving the stress on the slip systems. Will give the response of phenopowerlaw for an
!! untextured polycrystal
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) isotropic

  type :: tParameters
    real(pReal) :: &
      M, &                                                                                          !< Taylor factor
      dot_gamma_0, &                                                                                !< reference strain rate
      n, &                                                                                          !< stress exponent
      h_0, &
      h_ln, &
      xi_inf, &                                                                                     !< maximum critical stress
      a, &
      c_1, &
      c_4, &
      c_3, &
      c_2
    logical :: &
      dilatation
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type :: tIsotropicState
    real(pReal), pointer, dimension(:) :: &
      xi, &
      gamma
  end type tIsotropicState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
  type(tParameters),     allocatable, dimension(:) :: param
  type(tIsotropicState), allocatable, dimension(:) :: &
    dotState, &
    state

contains

!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_isotropic_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    ph, &
    Nmembers, &
    sizeState, sizeDotState
  real(pReal) :: &
    xi_0                                                                                            !< initial critical stress
  character(len=pStringLen) :: &
    extmsg = ''
  class(tNode), pointer :: &
    phases, &
    phase, &
    mech, &
    pl


  myPlasticity = plastic_active('isotropic')
  if(count(myPlasticity) == 0) return

  print'(/,a)', ' <<<+-  phase:mechanical:plastic:isotropic init  -+>>>'
  print'(a,i0)', ' # phases: ',count(myPlasticity); flush(IO_STDOUT)

  print*, 'T. Maiti and P. Eisenlohr, Scripta Materialia 145:37–40, 2018'
  print*, 'https://doi.org/10.1016/j.scriptamat.2017.09.047'

  phases => config_material%get('phase')
  allocate(param(phases%length))
  allocate(state(phases%length))
  allocate(dotState(phases%length))

  do ph = 1, phases%length
    if(.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), dot => dotState(ph), stt => state(ph))

    phase => phases%get(ph)
    mech  => phase%get('mechanical')
    pl  => mech%get('plastic')

#if defined (__GFORTRAN__)
    prm%output = output_as1dString(pl)
#else
    prm%output = pl%get_as1dString('output',defaultVal=emptyStringArray)
#endif

    xi_0            = pl%get_asFloat('xi_0')
    prm%xi_inf      = pl%get_asFloat('xi_inf')
    prm%dot_gamma_0 = pl%get_asFloat('dot_gamma_0')
    prm%n           = pl%get_asFloat('n')
    prm%h_0         = pl%get_asFloat('h_0')
    prm%M           = pl%get_asFloat('M')
    prm%h_ln        = pl%get_asFloat('h_ln', defaultVal=0.0_pReal)
    prm%c_1         = pl%get_asFloat('c_1',  defaultVal=0.0_pReal)
    prm%c_4         = pl%get_asFloat('c_4',  defaultVal=0.0_pReal)
    prm%c_3         = pl%get_asFloat('c_3',  defaultVal=0.0_pReal)
    prm%c_2         = pl%get_asFloat('c_2',  defaultVal=0.0_pReal)
    prm%a           = pl%get_asFloat('a')

    prm%dilatation  = pl%get_AsBool('dilatation',defaultVal = .false.)

!--------------------------------------------------------------------------------------------------
!  sanity checks
    if (xi_0            <  0.0_pReal) extmsg = trim(extmsg)//' xi_0'
    if (prm%dot_gamma_0 <= 0.0_pReal) extmsg = trim(extmsg)//' dot_gamma_0'
    if (prm%n           <= 0.0_pReal) extmsg = trim(extmsg)//' n'
    if (prm%a           <= 0.0_pReal) extmsg = trim(extmsg)//' a'
    if (prm%M           <= 0.0_pReal) extmsg = trim(extmsg)//' M'

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers = count(material_phaseID == ph)
    sizeDotState = size(['xi   ','gamma'])
    sizeState = sizeDotState

    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,0)

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    stt%xi  => plasticState(ph)%state   (1,:)
    stt%xi  =  xi_0
    dot%xi  => plasticState(ph)%dotState(1,:)
    plasticState(ph)%atol(1) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)
    if (plasticState(ph)%atol(1) < 0.0_pReal) extmsg = trim(extmsg)//' atol_xi'

    stt%gamma  => plasticState(ph)%state   (2,:)
    dot%gamma  => plasticState(ph)%dotState(2,:)
    plasticState(ph)%atol(2) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)
    if (plasticState(ph)%atol(2) < 0.0_pReal) extmsg = trim(extmsg)//' atol_gamma'
    ! global alias
    plasticState(ph)%slipRate => plasticState(ph)%dotState(2:2,:)

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(isotropic)')

  enddo

end function plastic_isotropic_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!--------------------------------------------------------------------------------------------------
module subroutine isotropic_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

  real(pReal), dimension(3,3),     intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pReal), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to the Mandel stress

  real(pReal), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                     intent(in) :: &
    ph, &
    en

  real(pReal), dimension(3,3) :: &
    Mp_dev                                                                                          !< deviatoric part of the Mandel stress
  real(pReal) :: &
    dot_gamma, &                                                                                    !< strainrate
    norm_Mp_dev, &                                                                                  !< norm of the deviatoric part of the Mandel stress
    squarenorm_Mp_dev                                                                               !< square of the norm of the deviatoric part of the Mandel stress
  integer :: &
    k, l, m, n

  associate(prm => param(ph), stt => state(ph))

  Mp_dev = math_deviatoric33(Mp)
  squarenorm_Mp_dev = math_tensordot(Mp_dev,Mp_dev)
  norm_Mp_dev = sqrt(squarenorm_Mp_dev)

  if (norm_Mp_dev > 0.0_pReal) then
    dot_gamma = prm%dot_gamma_0 * (sqrt(1.5_pReal) * norm_Mp_dev/(prm%M*stt%xi(en))) **prm%n

    Lp = dot_gamma/prm%M * Mp_dev/norm_Mp_dev
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = (prm%n-1.0_pReal) * Mp_dev(k,l)*Mp_dev(m,n) / squarenorm_Mp_dev
    forall (k=1:3,l=1:3) &
      dLp_dMp(k,l,k,l) = dLp_dMp(k,l,k,l) + 1.0_pReal
    forall (k=1:3,m=1:3) &
      dLp_dMp(k,k,m,m) = dLp_dMp(k,k,m,m) - 1.0_pReal/3.0_pReal
    dLp_dMp = dot_gamma / prm%M * dLp_dMp / norm_Mp_dev
  else
    Lp = 0.0_pReal
    dLp_dMp = 0.0_pReal
  end if

  end associate

end subroutine isotropic_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate inelastic velocity gradient and its tangent.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_isotropic_LiAndItsTangent(Li,dLi_dMi,Mi,ph,en)

  real(pReal), dimension(3,3),     intent(out) :: &
    Li                                                                                              !< inleastic velocity gradient
  real(pReal), dimension(3,3,3,3), intent(out)  :: &
    dLi_dMi                                                                                         !< derivative of Li with respect to Mandel stress

  real(pReal), dimension(3,3), intent(in) :: &
    Mi                                                                                              !< Mandel stress
  integer,                     intent(in) :: &
    ph, &
    en

  real(pReal) :: &
    tr                                                                                              !< trace of spherical part of Mandel stress (= 3 x pressure)
  integer :: &
    k, l, m, n

  associate(prm => param(ph), stt => state(ph))

  tr=math_trace33(math_spherical33(Mi))

  if (prm%dilatation .and. abs(tr) > 0.0_pReal) then                                                ! no stress or J2 plasticity --> Li and its derivative are zero
    Li = math_I3 &
       * prm%dot_gamma_0/prm%M * (3.0_pReal*prm%M*stt%xi(en))**(-prm%n) &
       * tr * abs(tr)**(prm%n-1.0_pReal)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) dLi_dMi(k,l,m,n) = prm%n / tr * Li(k,l) * math_I3(m,n)
  else
    Li      = 0.0_pReal
    dLi_dMi = 0.0_pReal
  endif

  end associate

 end subroutine plastic_isotropic_LiAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module subroutine isotropic_dotState(Mp,ph,en)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en

  real(pReal) :: &
    dot_gamma, &                                                                                    !< strainrate
    xi_inf_star, &                                                                                  !< saturation xi
    norm_Mp                                                                                         !< norm of the (deviatoric) Mandel stress

  associate(prm => param(ph), stt => state(ph), &
   dot => dotState(ph))

  if (prm%dilatation) then
    norm_Mp = sqrt(math_tensordot(Mp,Mp))
  else
    norm_Mp = sqrt(math_tensordot(math_deviatoric33(Mp),math_deviatoric33(Mp)))
  endif

  dot_gamma = prm%dot_gamma_0 * (sqrt(1.5_pReal) * norm_Mp /(prm%M*stt%xi(en))) **prm%n

  if (dot_gamma > 1e-12_pReal) then
    if (dEq0(prm%c_1)) then
      xi_inf_star = prm%xi_inf
    else
      xi_inf_star = prm%xi_inf &
                  + asinh( (dot_gamma / prm%c_1)**(1.0_pReal / prm%c_2))**(1.0_pReal / prm%c_3) &
                  / prm%c_4 * (dot_gamma / prm%dot_gamma_0)**(1.0_pReal / prm%n)
    endif
    dot%xi(en) = dot_gamma &
               * ( prm%h_0 + prm%h_ln * log(dot_gamma) ) &
               * abs( 1.0_pReal - stt%xi(en)/xi_inf_star )**prm%a &
               * sign(1.0_pReal, 1.0_pReal - stt%xi(en)/xi_inf_star)
  else
    dot%xi(en) = 0.0_pReal
  endif

  dot%gamma(en) = dot_gamma                                                                         ! ToDo: not really used

  end associate

end subroutine isotropic_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_isotropic_results(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(ph), stt => state(ph))
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case ('xi')
        call results_writeDataset(stt%xi,group,trim(prm%output(o)), &
                                    'resistance against plastic flow','Pa')
    end select
  enddo outputsLoop
  end associate

end subroutine plastic_isotropic_results


end submodule isotropic
