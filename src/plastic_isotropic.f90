!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for isotropic plasticity
!> @details Isotropic Plasticity which resembles the phenopowerlaw plasticity without
!! resolving the stress on the slip systems. Will give the response of phenopowerlaw for an
!! untextured polycrystal
!--------------------------------------------------------------------------------------------------
module plastic_isotropic
  use prec, only: &
    pReal
 
  implicit none
  private
  integer,           dimension(:,:),   allocatable, target, public :: &
    plastic_isotropic_sizePostResult                                                                !< size of each post result output
  character(len=64), dimension(:,:),   allocatable, target, public :: &
    plastic_isotropic_output                                                                        !< name of each post result output
 
  enum, bind(c)
    enumerator :: &
      undefined_ID, &
      xi_ID, &
      dot_gamma_ID
  end enum
 
  type, private :: tParameters
    real(pReal) :: &
      M, &                                                                                          !< Taylor factor
      xi_0, &                                                                                       !< initial critical stress
      dot_gamma_0, &                                                                                !< reference strain rate
      n, &                                                                                          !< stress exponent
      h0, &
      h_ln, &
      xi_inf, &                                                                                     !< maximum critical stress
      a, &
      c_1, &
      c_4, &
      c_3, &
      c_2, &
      aTol_xi, &
      aTol_gamma
    integer :: &
      of_debug = 0
    integer(kind(undefined_ID)), allocatable, dimension(:) :: &
      outputID
    logical :: &
      dilatation
  end type tParameters
 
  type, private :: tIsotropicState
    real(pReal), pointer, dimension(:) :: &
      xi, &
      gamma
  end type tIsotropicState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
 type(tParameters),     allocatable, dimension(:), private :: param
 type(tIsotropicState), allocatable, dimension(:), private :: &
   dotState, &
   state

 public :: &
   plastic_isotropic_init, &
   plastic_isotropic_LpAndItsTangent, &
   plastic_isotropic_LiAndItsTangent, &
   plastic_isotropic_dotState, &
   plastic_isotropic_postResults, &
   plastic_isotropic_results

contains

!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_init
  use prec, only: &
    pStringLen
  use debug, only: &
#ifdef DEBUG
    debug_e, &
    debug_i, &
    debug_g, &
    debug_levelExtensive, &
#endif
    debug_level, &
    debug_constitutive, &
    debug_levelBasic
  use IO, only: &
    IO_error
  use material, only: &
#ifdef DEBUG
    phasememberAt, &
#endif
    phase_plasticity, &
    phase_plasticityInstance, &
    phase_Noutput, &
    material_allocatePlasticState, &
    PLASTICITY_ISOTROPIC_label, &
    PLASTICITY_ISOTROPIC_ID, &
    material_phase, &
    plasticState
  use config, only: &
    config_phase
  use lattice
 
  implicit none
  integer :: &
    Ninstance, &
    p, i, &
    NipcMyPhase, &
    sizeState, sizeDotState
 
  character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
 
  integer(kind(undefined_ID)) :: &
    outputID
 
  character(len=pStringLen) :: &
    extmsg = ''
  character(len=65536), dimension(:), allocatable :: &
    outputs
 
  write(6,'(/,a)')   ' <<<+-  plastic_'//PLASTICITY_ISOTROPIC_label//' init  -+>>>'
 
  write(6,'(/,a)')   ' Maiti and Eisenlohr, Scripta Materialia 145:37–40, 2018'
  write(6,'(a)')     ' https://doi.org/10.1016/j.scriptamat.2017.09.047'
 
  Ninstance = count(phase_plasticity == PLASTICITY_ISOTROPIC_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance
 
  allocate(plastic_isotropic_sizePostResult(maxval(phase_Noutput),Ninstance),source=0)
  allocate(plastic_isotropic_output(maxval(phase_Noutput),Ninstance))
           plastic_isotropic_output = ''
 
  allocate(param(Ninstance))
  allocate(state(Ninstance))
  allocate(dotState(Ninstance))
 
  do p = 1, size(phase_plasticity)
    if (phase_plasticity(p) /= PLASTICITY_ISOTROPIC_ID) cycle
    associate(prm => param(phase_plasticityInstance(p)), &
              dot => dotState(phase_plasticityInstance(p)), &
              stt => state(phase_plasticityInstance(p)), &
              config => config_phase(p))
 
#ifdef DEBUG
    if  (p==material_phase(debug_g,debug_i,debug_e)) &
      prm%of_debug = phasememberAt(debug_g,debug_i,debug_e)
#endif
 
    prm%xi_0            = config%getFloat('tau0')
    prm%xi_inf          = config%getFloat('tausat')
    prm%dot_gamma_0     = config%getFloat('gdot0')
    prm%n               = config%getFloat('n')
    prm%h0              = config%getFloat('h0')
    prm%M               = config%getFloat('m')
    prm%h_ln            = config%getFloat('h0_slopelnrate', defaultVal=0.0_pReal)
    prm%c_1             = config%getFloat('tausat_sinhfita',defaultVal=0.0_pReal)
    prm%c_4             = config%getFloat('tausat_sinhfitb',defaultVal=0.0_pReal)
    prm%c_3             = config%getFloat('tausat_sinhfitc',defaultVal=0.0_pReal)
    prm%c_2             = config%getFloat('tausat_sinhfitd',defaultVal=0.0_pReal)
    prm%a               = config%getFloat('a')
    prm%aTol_xi         = config%getFloat('atol_flowstress',defaultVal=1.0_pReal)
    prm%aTol_gamma      = config%getFloat('atol_shear',     defaultVal=1.0e-6_pReal)
 
    prm%dilatation      = config%keyExists('/dilatation/')
 
!--------------------------------------------------------------------------------------------------
!  sanity checks
    extmsg = ''
    if (prm%aTol_gamma     <= 0.0_pReal) extmsg = trim(extmsg)//' aTol_gamma'
    if (prm%xi_0           < 0.0_pReal) extmsg = trim(extmsg)//' xi_0'
    if (prm%dot_gamma_0    <= 0.0_pReal) extmsg = trim(extmsg)//' dot_gamma_0'
    if (prm%n              <= 0.0_pReal) extmsg = trim(extmsg)//' n'
    if (prm%a              <= 0.0_pReal) extmsg = trim(extmsg)//' a'
    if (prm%M              <= 0.0_pReal) extmsg = trim(extmsg)//' m'
    if (prm%aTol_xi        <= 0.0_pReal) extmsg = trim(extmsg)//' atol_xi'
    if (prm%aTol_gamma     <= 0.0_pReal) extmsg = trim(extmsg)//' atol_shear'
 
!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') &
      call IO_error(211,ext_msg=trim(extmsg)//'('//PLASTICITY_ISOTROPIC_label//')')
 
!--------------------------------------------------------------------------------------------------
!  output pararameters
    outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
    allocate(prm%outputID(0))
    do i=1, size(outputs)
      outputID = undefined_ID
      select case(outputs(i))
 
        case ('flowstress')
          outputID = xi_ID
        case ('strainrate')
          outputID = dot_gamma_ID
 
      end select
 
      if (outputID /= undefined_ID) then
        plastic_isotropic_output(i,phase_plasticityInstance(p)) = outputs(i)
        plastic_isotropic_sizePostResult(i,phase_plasticityInstance(p)) = 1
        prm%outputID = [prm%outputID, outputID]
     endif
 
    enddo
 
!--------------------------------------------------------------------------------------------------
! allocate state arrays
    NipcMyPhase = count(material_phase == p)
    sizeDotState = size(['xi               ','accumulated_shear'])
    sizeState = sizeDotState
 
    call material_allocatePlasticState(p,NipcMyPhase,sizeState,sizeDotState,0, &
                                       1,0,0)
    plasticState(p)%sizePostResults = sum(plastic_isotropic_sizePostResult(:,phase_plasticityInstance(p)))
 
!--------------------------------------------------------------------------------------------------
! locally defined state aliases and initialization of state0 and aTolState
    stt%xi  => plasticState(p)%state   (1,:)
    stt%xi  = prm%xi_0
    dot%xi  => plasticState(p)%dotState(1,:)
    plasticState(p)%aTolState(1) = prm%aTol_xi
 
    stt%gamma  => plasticState(p)%state   (2,:)
    dot%gamma  => plasticState(p)%dotState(2,:)
    plasticState(p)%aTolState(2) = prm%aTol_gamma
    ! global alias
    plasticState(p)%slipRate        => plasticState(p)%dotState(2:2,:)
    plasticState(p)%accumulatedSlip => plasticState(p)%state   (2:2,:)
 
    plasticState(p)%state0 = plasticState(p)%state                                                  ! ToDo: this could be done centrally
 
    end associate
 
  enddo

end subroutine plastic_isotropic_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,of)
#ifdef DEBUG
  use debug, only: &
    debug_level, &
    debug_constitutive,&
    debug_levelExtensive, &
    debug_levelSelective
#endif
  use math, only: &
    math_deviatoric33, &
    math_mul33xx33
 
  implicit none
  real(pReal), dimension(3,3),     intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pReal), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to the Mandel stress
 
  real(pReal), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                     intent(in) :: &
    instance, &
    of
 
  real(pReal), dimension(3,3) :: &
    Mp_dev                                                                                          !< deviatoric part of the Mandel stress
  real(pReal) :: &
    dot_gamma, &                                                                                    !< strainrate
    norm_Mp_dev, &                                                                                  !< norm of the deviatoric part of the Mandel stress
    squarenorm_Mp_dev                                                                               !< square of the norm of the deviatoric part of the Mandel stress
  integer :: &
    k, l, m, n
 
  associate(prm => param(instance), stt => state(instance))
 
  Mp_dev = math_deviatoric33(Mp)
  squarenorm_Mp_dev = math_mul33xx33(Mp_dev,Mp_dev)
  norm_Mp_dev = sqrt(squarenorm_Mp_dev)
 
  if (norm_Mp_dev > 0.0_pReal) then
    dot_gamma = prm%dot_gamma_0 * (sqrt(1.5_pReal) * norm_Mp_dev/(prm%M*stt%xi(of))) **prm%n
 
    Lp = dot_gamma/prm%M * Mp_dev/norm_Mp_dev
#ifdef DEBUG
    if (iand(debug_level(debug_constitutive), debug_levelExtensive) /= 0 &
        .and. (of == prm%of_debug .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0)) then
      write(6,'(/,a,/,3(12x,3(f12.4,1x)/))') '<< CONST isotropic >> Tstar (dev) / MPa', &
                                       transpose(Mp_dev)*1.0e-6_pReal
      write(6,'(/,a,/,f12.5)') '<< CONST isotropic >> norm Tstar / MPa', norm_Mp_dev*1.0e-6_pReal
      write(6,'(/,a,/,f12.5)') '<< CONST isotropic >> gdot', dot_gamma
    end if
#endif
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

end subroutine plastic_isotropic_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
! ToDo: Rename Tstar to Mi?
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_LiAndItsTangent(Li,dLi_dTstar,Tstar,instance,of)
  use math, only: &
    math_I3, &
    math_spherical33, &
    math_mul33xx33
 
  implicit none
  real(pReal), dimension(3,3), intent(out) :: &
    Li                                                                                              !< inleastic velocity gradient
  real(pReal), dimension(3,3,3,3), intent(out)  :: &
    dLi_dTstar                                                                                      !< derivative of Li with respect to the Mandel stress
 
  real(pReal), dimension(3,3),   intent(in) :: &
    Tstar                                                                                           !< Mandel stress ToDo: Mi?
  integer,                       intent(in) :: &
    instance, &
    of
 
  real(pReal), dimension(3,3) :: &
    Tstar_sph                                                                                       !< sphiatoric part of the Mandel stress
  real(pReal) :: &
    dot_gamma, &                                                                                    !< strainrate
    norm_Tstar_sph, &                                                                               !< euclidean norm of Tstar_sph
    squarenorm_Tstar_sph                                                                            !< square of the euclidean norm of Tstar_sph
  integer :: &
    k, l, m, n
 
  associate(prm => param(instance), stt => state(instance))
 
  Tstar_sph = math_spherical33(Tstar)
  squarenorm_Tstar_sph = math_mul33xx33(Tstar_sph,Tstar_sph)
  norm_Tstar_sph = sqrt(squarenorm_Tstar_sph)
 
  if (prm%dilatation .and. norm_Tstar_sph > 0.0_pReal) then                                         ! no stress or J2 plastitiy --> Li and its derivative are zero
    dot_gamma = prm%dot_gamma_0 * (sqrt(1.5_pReal) * norm_Tstar_sph /(prm%M*stt%xi(of))) **prm%n
 
    Li = math_I3/sqrt(3.0_pReal) * dot_gamma/prm%M
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLi_dTstar(k,l,m,n) = (prm%n-1.0_pReal) * Tstar_sph(k,l)*Tstar_sph(m,n) / squarenorm_Tstar_sph
    forall (k=1:3,l=1:3) &
      dLi_dTstar(k,l,k,l) = dLi_dTstar(k,l,k,l) + 1.0_pReal
 
    dLi_dTstar = dot_gamma / prm%M * dLi_dTstar / norm_Tstar_sph
  else
    Li = 0.0_pReal
    dLi_dTstar = 0.0_pReal
  endif
 
  end associate

 end subroutine plastic_isotropic_LiAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_dotState(Mp,instance,of)
  use prec, only: &
    dEq0
  use math, only: &
    math_mul33xx33, &
    math_deviatoric33
 
  implicit none
  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    instance, &
    of
 
  real(pReal) :: &
    dot_gamma, &                                                                                    !< strainrate
    xi_inf_star, &                                                                                  !< saturation xi
    norm_Mp                                                                                         !< norm of the (deviatoric) Mandel stress
 
  associate(prm => param(instance), stt => state(instance), dot => dotState(instance))
 
  if (prm%dilatation) then
    norm_Mp = sqrt(math_mul33xx33(Mp,Mp))
  else
    norm_Mp = sqrt(math_mul33xx33(math_deviatoric33(Mp),math_deviatoric33(Mp)))
  endif
 
  dot_gamma = prm%dot_gamma_0 * (sqrt(1.5_pReal) * norm_Mp /(prm%M*stt%xi(of))) **prm%n
 
  if (dot_gamma > 1e-12_pReal) then
    if (dEq0(prm%c_1)) then
      xi_inf_star = prm%xi_inf
    else
      xi_inf_star = prm%xi_inf &
                  + asinh( (dot_gamma / prm%c_1)**(1.0_pReal / prm%c_2))**(1.0_pReal / prm%c_3) &
                     / prm%c_4 * (dot_gamma / prm%dot_gamma_0)**(1.0_pReal / prm%n)
    endif
    dot%xi(of) = dot_gamma &
               * ( prm%h0 + prm%h_ln * log(dot_gamma) ) &
               * abs( 1.0_pReal - stt%xi(of)/xi_inf_star )**prm%a &
               * sign(1.0_pReal, 1.0_pReal - stt%xi(of)/xi_inf_star)
  else
    dot%xi(of) = 0.0_pReal
  endif
 
  dot%gamma(of) = dot_gamma                                                                         ! ToDo: not really used
 
  end associate

end subroutine plastic_isotropic_dotState


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_isotropic_postResults(Mp,instance,of) result(postResults)
  use math, only: &
    math_mul33xx33, &
    math_deviatoric33
 
  implicit none
  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    instance, &
    of
 
  real(pReal), dimension(sum(plastic_isotropic_sizePostResult(:,instance))) :: &
    postResults
 
  real(pReal) :: &
    norm_Mp                                                                                         !< norm of the Mandel stress
  integer :: &
    o,c
 
  associate(prm => param(instance), stt => state(instance))
 
  if (prm%dilatation) then
    norm_Mp = sqrt(math_mul33xx33(Mp,Mp))
  else
    norm_Mp = sqrt(math_mul33xx33(math_deviatoric33(Mp),math_deviatoric33(Mp)))
  endif
 
  c = 0
 
  outputsLoop: do o = 1,size(prm%outputID)
    select case(prm%outputID(o))
 
      case (xi_ID)
        postResults(c+1) = stt%xi(of)
        c = c + 1
      case (dot_gamma_ID)
        postResults(c+1) = prm%dot_gamma_0 &
                         * (sqrt(1.5_pReal) * norm_Mp /(prm%M * stt%xi(of)))**prm%n
        c = c + 1
 
    end select
  enddo outputsLoop
 
  end associate

end function plastic_isotropic_postResults


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_results(instance,group)
#if defined(PETSc) || defined(DAMASK_HDF5)
  use results

  implicit none
  integer,          intent(in) :: instance
  character(len=*), intent(in) :: group
  
  integer :: o

  associate(prm => param(instance), stt => state(instance))
  outputsLoop: do o = 1,size(prm%outputID)
    select case(prm%outputID(o))
      case (xi_ID)
        call results_writeDataset(group,stt%xi,'xi','resistance against plastic flow','Pa')
    end select
  enddo outputsLoop
  end associate
  
#else
  integer,          intent(in) :: instance
  character(len=*), intent(in) :: group
#endif

end subroutine plastic_isotropic_results


end module plastic_isotropic
