!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Chen Zhang, Michigan State University
!> @brief crystallite state integration functions and reporting of results
!--------------------------------------------------------------------------------------------------

module crystallite
  use prec
  use IO
  use config
  use debug
  use numerics
  use rotations
  use math
  use FEsolving
  use material
  use constitutive
  use discretization
  use lattice
  use plastic_nonlocal
  use results
 
  implicit none 
  private
 
  real(pReal),               dimension(:,:,:),        allocatable, public :: &
    crystallite_dt                                                                                  !< requested time increment of each grain
  real(pReal),               dimension(:,:,:),        allocatable :: &
    crystallite_subdt, &                                                                            !< substepped time increment of each grain
    crystallite_subFrac, &                                                                          !< already calculated fraction of increment
    crystallite_subStep                                                                             !< size of next integration step
  type(rotation),            dimension(:,:,:),        allocatable :: &
    crystallite_orientation                                                                         !< current orientation 
  real(pReal),               dimension(:,:,:,:,:),    allocatable, public, protected :: &
    crystallite_Fe, &                                                                               !< current "elastic" def grad (end of converged time step)
    crystallite_P                                                                                   !< 1st Piola-Kirchhoff stress per grain
  real(pReal),               dimension(:,:,:,:,:),    allocatable, public :: &
    crystallite_S, &                                                                                !< current 2nd Piola-Kirchhoff stress vector (end of converged time step)
    crystallite_S0, &                                                                               !< 2nd Piola-Kirchhoff stress vector at start of FE inc
    crystallite_partionedS0, &                                                                      !< 2nd Piola-Kirchhoff stress vector at start of homog inc
    crystallite_Fp, &                                                                               !< current plastic def grad (end of converged time step)
    crystallite_Fp0, &                                                                              !< plastic def grad at start of FE inc
    crystallite_partionedFp0,&                                                                      !< plastic def grad at start of homog inc
    crystallite_Fi, &                                                                               !< current intermediate def grad (end of converged time step)
    crystallite_Fi0, &                                                                              !< intermediate def grad at start of FE inc
    crystallite_partionedFi0,&                                                                      !< intermediate def grad at start of homog inc
    crystallite_F0, &                                                                               !< def grad at start of FE inc
    crystallite_partionedF,  &                                                                      !< def grad to be reached at end of homog inc
    crystallite_partionedF0, &                                                                      !< def grad at start of homog inc
    crystallite_Lp, &                                                                               !< current plastic velocitiy grad (end of converged time step)
    crystallite_Lp0, &                                                                              !< plastic velocitiy grad at start of FE inc
    crystallite_partionedLp0, &                                                                     !< plastic velocity grad at start of homog inc
    crystallite_Li, &                                                                               !< current intermediate velocitiy grad (end of converged time step)
    crystallite_Li0, &                                                                              !< intermediate velocitiy grad at start of FE inc
    crystallite_partionedLi0                                                                        !< intermediate velocity grad at start of homog inc
  real(pReal),                dimension(:,:,:,:,:),    allocatable :: &
    crystallite_subS0, &                                                                            !< 2nd Piola-Kirchhoff stress vector at start of crystallite inc
    crystallite_invFp, &                                                                            !< inverse of current plastic def grad (end of converged time step)
    crystallite_subFp0,&                                                                            !< plastic def grad at start of crystallite inc
    crystallite_invFi, &                                                                            !< inverse of current intermediate def grad (end of converged time step)
    crystallite_subFi0,&                                                                            !< intermediate def grad at start of crystallite inc
    crystallite_subF,  &                                                                            !< def grad to be reached at end of crystallite inc
    crystallite_subF0, &                                                                            !< def grad at start of crystallite inc
    crystallite_subLp0,&                                                                            !< plastic velocity grad at start of crystallite inc
    crystallite_subLi0                                                                              !< intermediate velocity grad at start of crystallite inc
  real(pReal),                dimension(:,:,:,:,:,:,:), allocatable, public :: &
    crystallite_dPdF                                                                                !< current individual dPdF per grain (end of converged time step)
  logical,                    dimension(:,:,:),         allocatable, public :: &
    crystallite_requested                                                                           !< used by upper level (homogenization) to request crystallite calculation
  logical,                    dimension(:,:,:),         allocatable :: &
    crystallite_converged, &                                                                        !< convergence flag
    crystallite_todo, &                                                                             !< flag to indicate need for further computation
    crystallite_localPlasticity                                                                     !< indicates this grain to have purely local constitutive law
 
  type :: tOutput                                                                                   !< new requested output (per phase)
    character(len=65536), allocatable, dimension(:) :: &
      label
  end type tOutput
  type(tOutput), allocatable, dimension(:) :: output_constituent
  
  type :: tNumerics
    integer :: &
      iJacoLpresiduum, &                                                                            !< frequency of Jacobian update of residuum in Lp
      nState, &                                                                                     !< state loop limit
      nStress                                                                                       !< stress loop limit                                                                   
    real(pReal) :: &
      subStepMinCryst, &                                                                            !< minimum (relative) size of sub-step allowed during cutback
      subStepSizeCryst, &                                                                           !< size of first substep when cutback
      subStepSizeLp, &                                                                              !< size of first substep when cutback in Lp calculation
      subStepSizeLi, &                                                                              !< size of first substep when cutback in Li calculation
      stepIncreaseCryst, &                                                                          !< increase of next substep size when previous substep converged
      rTol_crystalliteState, &                                                                      !< relative tolerance in state loop 
      rTol_crystalliteStress, &                                                                     !< relative tolerance in stress loop
      aTol_crystalliteStress                                                                        !< absolute tolerance in stress loop
  end type tNumerics
  
  type(tNumerics) :: num                                                                            ! numerics parameters. Better name?
  
  procedure(), pointer :: integrateState
 
  public :: &
    crystallite_init, &
    crystallite_stress, &
    crystallite_stressTangent, &
    crystallite_orientations, &
    crystallite_push33ToRef, &
    crystallite_results

contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates and initialize per grain variables
!--------------------------------------------------------------------------------------------------
subroutine crystallite_init
 
  logical, dimension(:,:), allocatable :: devNull
  integer :: &
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e, &                                                                                            !< counter in element loop
    cMax, &                                                                                         !< maximum number of  integration point components
    iMax, &                                                                                         !< maximum number of integration points
    eMax, &                                                                                         !< maximum number of elements
    myNcomponents                                                                                   !< number of components at current IP
 
  write(6,'(/,a)')   ' <<<+-  crystallite init  -+>>>'
 
  cMax = homogenization_maxNgrains
  iMax = discretization_nIP
  eMax = discretization_nElem
 
  allocate(crystallite_S0(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_partionedS0(3,3,cMax,iMax,eMax),       source=0.0_pReal)
  allocate(crystallite_S(3,3,cMax,iMax,eMax),                 source=0.0_pReal)
  allocate(crystallite_subS0(3,3,cMax,iMax,eMax),             source=0.0_pReal)
  allocate(crystallite_P(3,3,cMax,iMax,eMax),                 source=0.0_pReal)
  allocate(crystallite_F0(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_partionedF0(3,3,cMax,iMax,eMax),       source=0.0_pReal)
  allocate(crystallite_partionedF(3,3,cMax,iMax,eMax),        source=0.0_pReal)
  allocate(crystallite_subF0(3,3,cMax,iMax,eMax),             source=0.0_pReal)
  allocate(crystallite_subF(3,3,cMax,iMax,eMax),              source=0.0_pReal)
  allocate(crystallite_Fp0(3,3,cMax,iMax,eMax),               source=0.0_pReal)
  allocate(crystallite_partionedFp0(3,3,cMax,iMax,eMax),      source=0.0_pReal)
  allocate(crystallite_subFp0(3,3,cMax,iMax,eMax),            source=0.0_pReal)
  allocate(crystallite_Fp(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_invFp(3,3,cMax,iMax,eMax),             source=0.0_pReal)
  allocate(crystallite_Fi0(3,3,cMax,iMax,eMax),               source=0.0_pReal)
  allocate(crystallite_partionedFi0(3,3,cMax,iMax,eMax),      source=0.0_pReal)
  allocate(crystallite_subFi0(3,3,cMax,iMax,eMax),            source=0.0_pReal)
  allocate(crystallite_Fi(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_invFi(3,3,cMax,iMax,eMax),             source=0.0_pReal)
  allocate(crystallite_Fe(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_Lp0(3,3,cMax,iMax,eMax),               source=0.0_pReal)
  allocate(crystallite_partionedLp0(3,3,cMax,iMax,eMax),      source=0.0_pReal)
  allocate(crystallite_subLp0(3,3,cMax,iMax,eMax),            source=0.0_pReal)
  allocate(crystallite_Lp(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_Li0(3,3,cMax,iMax,eMax),               source=0.0_pReal)
  allocate(crystallite_partionedLi0(3,3,cMax,iMax,eMax),      source=0.0_pReal)
  allocate(crystallite_subLi0(3,3,cMax,iMax,eMax),            source=0.0_pReal)
  allocate(crystallite_Li(3,3,cMax,iMax,eMax),                source=0.0_pReal)
  allocate(crystallite_dPdF(3,3,3,3,cMax,iMax,eMax),          source=0.0_pReal)
  allocate(crystallite_dt(cMax,iMax,eMax),                    source=0.0_pReal)
  allocate(crystallite_subdt(cMax,iMax,eMax),                 source=0.0_pReal)
  allocate(crystallite_subFrac(cMax,iMax,eMax),               source=0.0_pReal)
  allocate(crystallite_subStep(cMax,iMax,eMax),               source=0.0_pReal)
  allocate(crystallite_orientation(cMax,iMax,eMax))
  allocate(crystallite_localPlasticity(cMax,iMax,eMax),       source=.true.)
  allocate(crystallite_requested(cMax,iMax,eMax),             source=.false.)
  allocate(crystallite_todo(cMax,iMax,eMax),                  source=.false.)
  allocate(crystallite_converged(cMax,iMax,eMax),             source=.true.)
                                      
  num%subStepMinCryst        = config_numerics%getFloat('substepmincryst',       defaultVal=1.0e-3_pReal)
  num%subStepSizeCryst       = config_numerics%getFloat('substepsizecryst',      defaultVal=0.25_pReal)
  num%stepIncreaseCryst      = config_numerics%getFloat('stepincreasecryst',     defaultVal=1.5_pReal)
  
  num%subStepSizeLp          = config_numerics%getFloat('substepsizelp',         defaultVal=0.5_pReal)
  num%subStepSizeLi          = config_numerics%getFloat('substepsizeli',         defaultVal=0.5_pReal)

  num%rTol_crystalliteState  = config_numerics%getFloat('rtol_crystallitestate', defaultVal=1.0e-6_pReal)
  num%rTol_crystalliteStress = config_numerics%getFloat('rtol_crystallitestress',defaultVal=1.0e-6_pReal)
  num%aTol_crystalliteStress = config_numerics%getFloat('atol_crystallitestress',defaultVal=1.0e-8_pReal)
  
  num%iJacoLpresiduum        = config_numerics%getInt  ('ijacolpresiduum',       defaultVal=1)
  
  num%nState                 = config_numerics%getInt  ('nstate',                defaultVal=20)
  num%nStress                = config_numerics%getInt  ('nstress',               defaultVal=40)
  
  if(num%subStepMinCryst   <= 0.0_pReal)      call IO_error(301,ext_msg='subStepMinCryst')
  if(num%subStepSizeCryst  <= 0.0_pReal)      call IO_error(301,ext_msg='subStepSizeCryst')
  if(num%stepIncreaseCryst <= 0.0_pReal)      call IO_error(301,ext_msg='stepIncreaseCryst')

  if(num%subStepSizeLp <= 0.0_pReal)          call IO_error(301,ext_msg='subStepSizeLp')
  if(num%subStepSizeLi <= 0.0_pReal)          call IO_error(301,ext_msg='subStepSizeLi')

  if(num%rTol_crystalliteState  <= 0.0_pReal) call IO_error(301,ext_msg='rTol_crystalliteState')
  if(num%rTol_crystalliteStress <= 0.0_pReal) call IO_error(301,ext_msg='rTol_crystalliteStress')
  if(num%aTol_crystalliteStress <= 0.0_pReal) call IO_error(301,ext_msg='aTol_crystalliteStress')
  
  if(num%iJacoLpresiduum < 1)                 call IO_error(301,ext_msg='iJacoLpresiduum')
  
  if(num%nState < 1)                          call IO_error(301,ext_msg='nState')
  if(num%nStress< 1)                          call IO_error(301,ext_msg='nStress')
 
  select case(numerics_integrator)
    case(1)
      integrateState => integrateStateFPI
    case(2)
      integrateState => integrateStateEuler
    case(3)
      integrateState => integrateStateAdaptiveEuler
    case(4)
      integrateState => integrateStateRK4
    case(5)
      integrateState => integrateStateRKCK45
  end select
 
  allocate(output_constituent(size(config_phase)))
  do c = 1, size(config_phase)
#if defined(__GFORTRAN__)
    allocate(output_constituent(c)%label(1)) 
    output_constituent(c)%label(1)= 'GfortranBug86277'
    output_constituent(c)%label  = config_phase(c)%getStrings('(output)',defaultVal=output_constituent(c)%label )
    if (output_constituent(c)%label (1) == 'GfortranBug86277') output_constituent(c)%label  = [character(len=pStringLen)::]
#else
    output_constituent(c)%label  = config_phase(c)%getStrings('(output)',defaultVal=[character(len=pStringLen)::])
#endif
  enddo

  call config_deallocate('material.config/phase')

!--------------------------------------------------------------------------------------------------
! initialize
 !$OMP PARALLEL DO PRIVATE(myNcomponents,i,c)
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    myNcomponents = homogenization_Ngrains(material_homogenizationAt(e))
    do i = FEsolving_execIP(1,e), FEsolving_execIP(2,e); do c = 1, myNcomponents
      crystallite_Fp0(1:3,1:3,c,i,e) = math_EulerToR(material_Eulers(1:3,c,i,e))                    ! plastic def gradient reflects init orientation
      crystallite_Fi0(1:3,1:3,c,i,e) = constitutive_initialFi(c,i,e)
      crystallite_F0(1:3,1:3,c,i,e)  = math_I3
      crystallite_localPlasticity(c,i,e) = phase_localPlasticity(material_phaseAt(c,e))
      crystallite_Fe(1:3,1:3,c,i,e)  = math_inv33(matmul(crystallite_Fi0(1:3,1:3,c,i,e), &
                                                         crystallite_Fp0(1:3,1:3,c,i,e)))           ! assuming that euler angles are given in internal strain free configuration
      crystallite_Fp(1:3,1:3,c,i,e)  = crystallite_Fp0(1:3,1:3,c,i,e)
      crystallite_Fi(1:3,1:3,c,i,e)  = crystallite_Fi0(1:3,1:3,c,i,e)
      crystallite_requested(c,i,e) = .true.
    enddo; enddo
  enddo
  !$OMP END PARALLEL DO
 
  if(any(.not. crystallite_localPlasticity) .and. .not. usePingPong) call IO_error(601)             ! exit if nonlocal but no ping-pong ToDo: Why not check earlier? or in nonlocal?
 
  crystallite_partionedFp0 = crystallite_Fp0
  crystallite_partionedFi0 = crystallite_Fi0
  crystallite_partionedF0  = crystallite_F0
  crystallite_partionedF   = crystallite_F0
 
  call crystallite_orientations()
 
  !$OMP PARALLEL DO
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      do c = 1,homogenization_Ngrains(material_homogenizationAt(e))
        call constitutive_microstructure(crystallite_Fe(1:3,1:3,c,i,e), &
                                         crystallite_Fp(1:3,1:3,c,i,e), &
                                         c,i,e)                                                     ! update dependent state variables to be consistent with basic states
     enddo
    enddo
  enddo
  !$OMP END PARALLEL DO
 
  devNull = crystallite_stress()
  call crystallite_stressTangent

#ifdef DEBUG
  if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) then
    write(6,'(a42,1x,i10)') '    # of elements:                       ', eMax
    write(6,'(a42,1x,i10)') 'max # of integration points/element:     ', iMax
    write(6,'(a42,1x,i10)') 'max # of constituents/integration point: ', cMax
    write(6,'(a42,1x,i10)') '    # of nonlocal constituents:          ',count(.not. crystallite_localPlasticity)
    flush(6)
  endif
 
  call debug_info
  call debug_reset
#endif

end subroutine crystallite_init


!--------------------------------------------------------------------------------------------------
!> @brief calculate stress (P)
!--------------------------------------------------------------------------------------------------
function crystallite_stress(dummyArgumentToPreventInternalCompilerErrorWithGCC)
 
  logical, dimension(discretization_nIP,discretization_nElem) :: crystallite_stress
  real(pReal), intent(in), optional :: &
    dummyArgumentToPreventInternalCompilerErrorWithGCC
  real(pReal) :: &
    formerSubStep
  integer :: &
    NiterationCrystallite, &                                                                        ! number of iterations in crystallite loop
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e, &                                                                                            !< counter in element loop
    startIP, endIP, &
    s
 
#ifdef DEBUG
  if (iand(debug_level(debug_crystallite),debug_levelSelective) /= 0 &
      .and. FEsolving_execElem(1) <= debug_e &
      .and.                          debug_e <= FEsolving_execElem(2)) then
      write(6,'(/,a,i8,1x,i2,1x,i3)')    '<< CRYST stress >> boundary and initial values at el ip ipc ', &
        debug_e,debug_i, debug_g
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST stress >> F  ', &
                                          transpose(crystallite_partionedF(1:3,1:3,debug_g,debug_i,debug_e))
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST stress >> F0 ', &
                                          transpose(crystallite_partionedF0(1:3,1:3,debug_g,debug_i,debug_e))
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST stress >> Fp0', &
                                          transpose(crystallite_partionedFp0(1:3,1:3,debug_g,debug_i,debug_e))
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST stress >> Fi0', &
                                          transpose(crystallite_partionedFi0(1:3,1:3,debug_g,debug_i,debug_e))
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST stress >> Lp0', &
                                          transpose(crystallite_partionedLp0(1:3,1:3,debug_g,debug_i,debug_e))
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST stress >> Li0', &
                                          transpose(crystallite_partionedLi0(1:3,1:3,debug_g,debug_i,debug_e))
  endif
#endif

!--------------------------------------------------------------------------------------------------
! initialize to starting condition
  crystallite_subStep = 0.0_pReal
  !$OMP PARALLEL DO
  elementLooping1: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e); do c = 1,homogenization_Ngrains(material_homogenizationAt(e))
      homogenizationRequestsCalculation: if (crystallite_requested(c,i,e)) then
        plasticState    (material_phaseAt(c,e))%subState0(      :,material_phaseMemberAt(c,i,e)) = &
        plasticState    (material_phaseAt(c,e))%partionedState0(:,material_phaseMemberAt(c,i,e))

        do s = 1, phase_Nsources(material_phaseAt(c,e))
          sourceState(material_phaseAt(c,e))%p(s)%subState0(      :,material_phaseMemberAt(c,i,e)) = &
          sourceState(material_phaseAt(c,e))%p(s)%partionedState0(:,material_phaseMemberAt(c,i,e))
        enddo
        crystallite_subFp0(1:3,1:3,c,i,e) = crystallite_partionedFp0(1:3,1:3,c,i,e)
        crystallite_subLp0(1:3,1:3,c,i,e) = crystallite_partionedLp0(1:3,1:3,c,i,e)
        crystallite_subFi0(1:3,1:3,c,i,e) = crystallite_partionedFi0(1:3,1:3,c,i,e)
        crystallite_subLi0(1:3,1:3,c,i,e) = crystallite_partionedLi0(1:3,1:3,c,i,e)
        crystallite_subF0(1:3,1:3,c,i,e)  = crystallite_partionedF0(1:3,1:3,c,i,e)
        crystallite_subS0(1:3,1:3,c,i,e)  = crystallite_partionedS0(1:3,1:3,c,i,e)
        crystallite_subFrac(c,i,e) = 0.0_pReal
        crystallite_subStep(c,i,e) = 1.0_pReal/num%subStepSizeCryst
        crystallite_todo(c,i,e) = .true.
        crystallite_converged(c,i,e) = .false.                                                      ! pretend failed step of 1/subStepSizeCryst
      endif homogenizationRequestsCalculation
    enddo; enddo
  enddo elementLooping1
  !$OMP END PARALLEL DO

  singleRun: if (FEsolving_execELem(1) == FEsolving_execElem(2) .and. &
      FEsolving_execIP(1,FEsolving_execELem(1))==FEsolving_execIP(2,FEsolving_execELem(1))) then
    startIP = FEsolving_execIP(1,FEsolving_execELem(1))
    endIP   = startIP
  else singleRun
    startIP = 1
    endIP = discretization_nIP
  endif singleRun

  NiterationCrystallite = 0
  cutbackLooping: do while (any(crystallite_todo(:,startIP:endIP,FEsolving_execELem(1):FEsolving_execElem(2))))
    NiterationCrystallite = NiterationCrystallite + 1

#ifdef DEBUG
    if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0) &
      write(6,'(a,i6)') '<< CRYST stress >> crystallite iteration ',NiterationCrystallite
#endif
    !$OMP PARALLEL DO PRIVATE(formerSubStep)
    elementLooping3: do e = FEsolving_execElem(1),FEsolving_execElem(2)
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
        do c = 1,homogenization_Ngrains(material_homogenizationAt(e))
!--------------------------------------------------------------------------------------------------
!  wind forward
          if (crystallite_converged(c,i,e)) then
            formerSubStep = crystallite_subStep(c,i,e)
            crystallite_subFrac(c,i,e) = crystallite_subFrac(c,i,e) + crystallite_subStep(c,i,e)
            crystallite_subStep(c,i,e) = min(1.0_pReal - crystallite_subFrac(c,i,e), &
                                             num%stepIncreaseCryst * crystallite_subStep(c,i,e))

            crystallite_todo(c,i,e) = crystallite_subStep(c,i,e) > 0.0_pReal                        ! still time left to integrate on?
            if (crystallite_todo(c,i,e)) then
              crystallite_subF0 (1:3,1:3,c,i,e) = crystallite_subF(1:3,1:3,c,i,e)
              crystallite_subLp0(1:3,1:3,c,i,e) = crystallite_Lp  (1:3,1:3,c,i,e)
              crystallite_subLi0(1:3,1:3,c,i,e) = crystallite_Li  (1:3,1:3,c,i,e)
              crystallite_subFp0(1:3,1:3,c,i,e) = crystallite_Fp  (1:3,1:3,c,i,e)
              crystallite_subFi0(1:3,1:3,c,i,e) = crystallite_Fi  (1:3,1:3,c,i,e)
              crystallite_subS0 (1:3,1:3,c,i,e) = crystallite_S   (1:3,1:3,c,i,e)
              !if abbrevation, make c and p private in omp
              plasticState(    material_phaseAt(c,e))%subState0(:,material_phaseMemberAt(c,i,e)) &
                = plasticState(material_phaseAt(c,e))%state(    :,material_phaseMemberAt(c,i,e))
              do s = 1, phase_Nsources(material_phaseAt(c,e))
                sourceState(    material_phaseAt(c,e))%p(s)%subState0(:,material_phaseMemberAt(c,i,e)) &
                  = sourceState(material_phaseAt(c,e))%p(s)%state(    :,material_phaseMemberAt(c,i,e))
              enddo
#ifdef DEBUG
              if (iand(debug_level(debug_crystallite),debug_levelBasic) /= 0 &
                  .and. ((e == debug_e .and. i == debug_i .and. c == debug_g) &
                         .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) &
                write(6,'(a,f12.8,a,f12.8,a,i8,1x,i2,1x,i3,/)') '<< CRYST stress >> winding forward from ', &
                  crystallite_subFrac(c,i,e)-formerSubStep,' to current crystallite_subfrac ', &
                  crystallite_subFrac(c,i,e),' in crystallite_stress at el ip ipc ',e,i,c
#endif
            endif

!--------------------------------------------------------------------------------------------------
!  cut back (reduced time and restore)
          else
            crystallite_subStep(c,i,e)       = num%subStepSizeCryst * crystallite_subStep(c,i,e)
            crystallite_Fp   (1:3,1:3,c,i,e) =            crystallite_subFp0(1:3,1:3,c,i,e)
            crystallite_invFp(1:3,1:3,c,i,e) = math_inv33(crystallite_Fp    (1:3,1:3,c,i,e))
            crystallite_Fi   (1:3,1:3,c,i,e) =            crystallite_subFi0(1:3,1:3,c,i,e)
            crystallite_invFi(1:3,1:3,c,i,e) = math_inv33(crystallite_Fi    (1:3,1:3,c,i,e))
            crystallite_S    (1:3,1:3,c,i,e) =            crystallite_S0    (1:3,1:3,c,i,e)
            if (crystallite_subStep(c,i,e) < 1.0_pReal) then                                        ! actual (not initial) cutback
              crystallite_Lp (1:3,1:3,c,i,e) =            crystallite_subLp0(1:3,1:3,c,i,e)
              crystallite_Li (1:3,1:3,c,i,e) =            crystallite_subLi0(1:3,1:3,c,i,e)
            endif
            plasticState    (material_phaseAt(c,e))%state(    :,material_phaseMemberAt(c,i,e)) &
              = plasticState(material_phaseAt(c,e))%subState0(:,material_phaseMemberAt(c,i,e))
            do s = 1, phase_Nsources(material_phaseAt(c,e))
              sourceState(    material_phaseAt(c,e))%p(s)%state(    :,material_phaseMemberAt(c,i,e)) &
                = sourceState(material_phaseAt(c,e))%p(s)%subState0(:,material_phaseMemberAt(c,i,e))
            enddo

                                                                                                    ! cant restore dotState here, since not yet calculated in first cutback after initialization
            crystallite_todo(c,i,e) = crystallite_subStep(c,i,e) > num%subStepMinCryst              ! still on track or already done (beyond repair)
#ifdef DEBUG
            if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
               .and. ((e == debug_e .and. i == debug_i .and. c == debug_g) &
                      .or. .not. iand(debug_level(debug_crystallite),debug_levelSelective) /= 0)) then
              if (crystallite_todo(c,i,e)) then
                write(6,'(a,f12.8,a,i8,1x,i2,1x,i3,/)') '<< CRYST stress >> cutback with new crystallite_subStep: ', &
                                                        crystallite_subStep(c,i,e),' at el ip ipc ',e,i,c
              else
                write(6,'(a,i8,1x,i2,1x,i3,/)') '<< CRYST stress >> reached minimum step size at el ip ipc ',e,i,c
              endif
            endif
#endif
          endif

!--------------------------------------------------------------------------------------------------
!  prepare for integration
          if (crystallite_todo(c,i,e)) then
            crystallite_subF(1:3,1:3,c,i,e) = crystallite_subF0(1:3,1:3,c,i,e) &
                                            + crystallite_subStep(c,i,e) * (crystallite_partionedF (1:3,1:3,c,i,e) &
                                                                          - crystallite_partionedF0(1:3,1:3,c,i,e))
            crystallite_Fe(1:3,1:3,c,i,e) = matmul(matmul(crystallite_subF (1:3,1:3,c,i,e), &
                                                          crystallite_invFp(1:3,1:3,c,i,e)), &
                                                   crystallite_invFi(1:3,1:3,c,i,e))
            crystallite_subdt(c,i,e) = crystallite_subStep(c,i,e) * crystallite_dt(c,i,e)
            crystallite_converged(c,i,e) = .false.
          endif

        enddo
      enddo
    enddo elementLooping3
    !$OMP END PARALLEL DO

#ifdef DEBUG
    if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0) then
      write(6,'(/,a,f8.5,a,f8.5,/)') '<< CRYST stress >> ',minval(crystallite_subStep), &
                                           ' ≤ subStep ≤ ',maxval(crystallite_subStep)
      write(6,'(/,a,f8.5,a,f8.5,/)') '<< CRYST stress >> ',minval(crystallite_subFrac), &
                                           ' ≤ subFrac ≤ ',maxval(crystallite_subFrac)
      flush(6)
      if (iand(debug_level(debug_crystallite),debug_levelSelective) /= 0) then
        write(6,'(/,a,f8.5,1x,a,1x,f8.5,1x,a)') '<< CRYST stress >> subFrac + subStep = ',&
           crystallite_subFrac(debug_g,debug_i,debug_e),'+',crystallite_subStep(debug_g,debug_i,debug_e),'@selective'
        flush(6)
      endif
    endif
#endif
!--------------------------------------------------------------------------------------------------
!  integrate --- requires fully defined state array (basic + dependent state)
    if (any(crystallite_todo)) call integrateState                                                  ! TODO: unroll into proper elementloop to avoid N^2 for single point evaluation
    where(.not. crystallite_converged .and. crystallite_subStep > num%subStepMinCryst) &            ! do not try non-converged but fully cutbacked any further
      crystallite_todo = .true.                                                                     ! TODO: again unroll this into proper elementloop to avoid N^2 for single point evaluation


  enddo cutbackLooping

! return whether converged or not
  crystallite_stress = .false.
  elementLooping5: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      crystallite_stress(i,e) = all(crystallite_converged(:,i,e)) 
    enddo
  enddo elementLooping5

#ifdef DEBUG
  elementLooping6: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      do c = 1,homogenization_Ngrains(material_homogenizationAt(e))
        if (.not. crystallite_converged(c,i,e)) then
          if(iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
            write(6,'(a,i8,1x,i2,1x,i3,/)') '<< CRYST stress >> no convergence at el ip ipc ', &
              e,i,c
        endif
        if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
            .and. ((e == debug_e .and. i == debug_i .and. c == debug_g) &
                   .or. .not. iand(debug_level(debug_crystallite),debug_levelSelective) /= 0)) then
          write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST stress >> solution at el ip ipc ',e,i,c
          write(6,'(/,a,/,3(12x,3(f12.4,1x)/))') '<< CRYST stress >> P / MPa', &
                                           transpose(crystallite_P(1:3,1:3,c,i,e))*1.0e-6_pReal
          write(6,'(a,/,3(12x,3(f14.9,1x)/))')   '<< CRYST stress >> Fp', &
                                           transpose(crystallite_Fp(1:3,1:3,c,i,e))
          write(6,'(a,/,3(12x,3(f14.9,1x)/))')   '<< CRYST stress >> Fi', &
                                           transpose(crystallite_Fi(1:3,1:3,c,i,e))
          write(6,'(a,/,3(12x,3(f14.9,1x)/),/)') '<< CRYST stress >> Lp', &
                                           transpose(crystallite_Lp(1:3,1:3,c,i,e))
          write(6,'(a,/,3(12x,3(f14.9,1x)/),/)') '<< CRYST stress >> Li', &
                                           transpose(crystallite_Li(1:3,1:3,c,i,e))
          flush(6)
        endif
      enddo
    enddo
  enddo elementLooping6
#endif

end function crystallite_stress


!--------------------------------------------------------------------------------------------------
!> @brief calculate tangent (dPdF)
!--------------------------------------------------------------------------------------------------
subroutine crystallite_stressTangent

  integer :: &
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e, &                                                                                            !< counter in element loop
    o, &
    p

  real(pReal), dimension(3,3)     ::   temp_33_1, devNull,invSubFi0, temp_33_2, temp_33_3, temp_33_4
  real(pReal), dimension(3,3,3,3) ::   dSdFe, &
                                       dSdF, &
                                       dSdFi, &
                                       dLidS, &
                                       dLidFi, &
                                       dLpdS, &
                                       dLpdFi, &
                                       dFidS, &
                                       dFpinvdF, &
                                       rhs_3333, &
                                       lhs_3333, &
                                       temp_3333
  real(pReal), dimension(9,9)::        temp_99
  logical :: error

  !$OMP PARALLEL DO PRIVATE(dSdF,dSdFe,dSdFi,dLpdS,dLpdFi,dFpinvdF,dLidS,dLidFi,dFidS,invSubFi0,o,p, &
  !$OMP                     rhs_3333,lhs_3333,temp_99,temp_33_1,temp_33_2,temp_33_3,temp_33_4,temp_3333,error)
  elementLooping: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      do c = 1,homogenization_Ngrains(material_homogenizationAt(e))

        call constitutive_SandItsTangents(devNull,dSdFe,dSdFi, &
                                         crystallite_Fe(1:3,1:3,c,i,e), &
                                         crystallite_Fi(1:3,1:3,c,i,e),c,i,e)                       ! call constitutive law to calculate elastic stress tangent
        call constitutive_LiAndItsTangents(devNull,dLidS,dLidFi, &
                                           crystallite_S (1:3,1:3,c,i,e), &
                                           crystallite_Fi(1:3,1:3,c,i,e), &
                                           c,i,e)                                                   ! call constitutive law to calculate Li tangent in lattice configuration

        if (sum(abs(dLidS)) < tol_math_check) then
          dFidS = 0.0_pReal
        else
          invSubFi0 = math_inv33(crystallite_subFi0(1:3,1:3,c,i,e))
          lhs_3333 = 0.0_pReal; rhs_3333 = 0.0_pReal
          do o=1,3; do p=1,3
            lhs_3333(1:3,1:3,o,p) = lhs_3333(1:3,1:3,o,p) &
                                  + crystallite_subdt(c,i,e)*matmul(invSubFi0,dLidFi(1:3,1:3,o,p))
            lhs_3333(1:3,o,1:3,p) = lhs_3333(1:3,o,1:3,p) &
                                  + crystallite_invFi(1:3,1:3,c,i,e)*crystallite_invFi(p,o,c,i,e)
            rhs_3333(1:3,1:3,o,p) = rhs_3333(1:3,1:3,o,p) &
                                  - crystallite_subdt(c,i,e)*matmul(invSubFi0,dLidS(1:3,1:3,o,p))
          enddo; enddo
          call math_invert(temp_99,error,math_3333to99(lhs_3333))
          if (error) then
            call IO_warning(warning_ID=600,el=e,ip=i,g=c, &
                            ext_msg='inversion error in analytic tangent calculation')
            dFidS = 0.0_pReal
          else
            dFidS = math_mul3333xx3333(math_99to3333(temp_99),rhs_3333)
          endif
          dLidS = math_mul3333xx3333(dLidFi,dFidS) + dLidS
        endif

        call constitutive_LpAndItsTangents(devNull,dLpdS,dLpdFi, &
                                           crystallite_S (1:3,1:3,c,i,e), &
                                           crystallite_Fi(1:3,1:3,c,i,e),c,i,e)                     ! call constitutive law to calculate Lp tangent in lattice configuration
        dLpdS = math_mul3333xx3333(dLpdFi,dFidS) + dLpdS

!--------------------------------------------------------------------------------------------------
! calculate dSdF
        temp_33_1 = transpose(matmul(crystallite_invFp(1:3,1:3,c,i,e), &
                                            crystallite_invFi(1:3,1:3,c,i,e)))
        temp_33_2 = matmul(           crystallite_subF  (1:3,1:3,c,i,e), &
                                  math_inv33(crystallite_subFp0(1:3,1:3,c,i,e)))
        temp_33_3 = matmul(matmul(crystallite_subF  (1:3,1:3,c,i,e), &
                                                crystallite_invFp (1:3,1:3,c,i,e)), &
                                     math_inv33(crystallite_subFi0(1:3,1:3,c,i,e)))

        do o=1,3; do p=1,3
          rhs_3333(p,o,1:3,1:3)  = matmul(dSdFe(p,o,1:3,1:3),temp_33_1)
          temp_3333(1:3,1:3,p,o) = matmul(matmul(temp_33_2,dLpdS(1:3,1:3,p,o)), &
                                                 crystallite_invFi(1:3,1:3,c,i,e)) &
                                 + matmul(temp_33_3,dLidS(1:3,1:3,p,o))
        enddo; enddo
        lhs_3333 = crystallite_subdt(c,i,e)*math_mul3333xx3333(dSdFe,temp_3333) &
                 + math_mul3333xx3333(dSdFi,dFidS)

        call math_invert(temp_99,error,math_identity2nd(9)+math_3333to99(lhs_3333))
        if (error) then
          call IO_warning(warning_ID=600,el=e,ip=i,g=c, &
                          ext_msg='inversion error in analytic tangent calculation')
          dSdF = rhs_3333
        else
          dSdF = math_mul3333xx3333(math_99to3333(temp_99),rhs_3333)
        endif

!--------------------------------------------------------------------------------------------------
! calculate dFpinvdF
        temp_3333 = math_mul3333xx3333(dLpdS,dSdF)
        do o=1,3; do p=1,3
          dFpinvdF(1:3,1:3,p,o) &
            = -crystallite_subdt(c,i,e) &
            * matmul(math_inv33(crystallite_subFp0(1:3,1:3,c,i,e)), &
                            matmul(temp_3333(1:3,1:3,p,o),crystallite_invFi(1:3,1:3,c,i,e)))
        enddo; enddo

!--------------------------------------------------------------------------------------------------
! assemble dPdF
        temp_33_1 = matmul(crystallite_invFp(1:3,1:3,c,i,e), &
                                  matmul(crystallite_S(1:3,1:3,c,i,e), &
                                                transpose(crystallite_invFp(1:3,1:3,c,i,e))))
        temp_33_2 = matmul(crystallite_S(1:3,1:3,c,i,e), &
                                  transpose(crystallite_invFp(1:3,1:3,c,i,e)))
        temp_33_3 = matmul(crystallite_subF(1:3,1:3,c,i,e), &
                                  crystallite_invFp(1:3,1:3,c,i,e))
        temp_33_4 = matmul(matmul(crystallite_subF(1:3,1:3,c,i,e), &
                                                crystallite_invFp(1:3,1:3,c,i,e)), &
                                               crystallite_S(1:3,1:3,c,i,e))

        crystallite_dPdF(1:3,1:3,1:3,1:3,c,i,e) = 0.0_pReal
        do p=1,3
          crystallite_dPdF(p,1:3,p,1:3,c,i,e) = transpose(temp_33_1)
        enddo
        do o=1,3; do p=1,3
          crystallite_dPdF(1:3,1:3,p,o,c,i,e) = crystallite_dPdF(1:3,1:3,p,o,c,i,e) + &
            matmul(matmul(crystallite_subF(1:3,1:3,c,i,e),dFpinvdF(1:3,1:3,p,o)),temp_33_2) + &
            matmul(matmul(temp_33_3,dSdF(1:3,1:3,p,o)),transpose(crystallite_invFp(1:3,1:3,c,i,e))) + &
            matmul(temp_33_4,transpose(dFpinvdF(1:3,1:3,p,o)))
        enddo; enddo

    enddo; enddo
  enddo elementLooping
  !$OMP END PARALLEL DO

end subroutine crystallite_stressTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates orientations
!--------------------------------------------------------------------------------------------------
subroutine crystallite_orientations
 
  integer &
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e                                                                                               !< counter in element loop
 
  !$OMP PARALLEL DO
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      do c = 1,homogenization_Ngrains(material_homogenizationAt(e))
        call crystallite_orientation(c,i,e)%fromMatrix(transpose(math_rotationalPart33(crystallite_Fe(1:3,1:3,c,i,e))))
  enddo; enddo; enddo
  !$OMP END PARALLEL DO
  
  nonlocalPresent: if (any(plasticState%nonLocal)) then
    !$OMP PARALLEL DO
    do e = FEsolving_execElem(1),FEsolving_execElem(2)
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
        if (plasticState(material_phaseAt(1,e))%nonLocal) &                                         ! if nonlocal model
          call plastic_nonlocal_updateCompatibility(crystallite_orientation,i,e)
    enddo; enddo
    !$OMP END PARALLEL DO
  endif nonlocalPresent

end subroutine crystallite_orientations


!--------------------------------------------------------------------------------------------------
!> @brief Map 2nd order tensor to reference config
!--------------------------------------------------------------------------------------------------
function crystallite_push33ToRef(ipc,ip,el, tensor33)
 
  real(pReal), dimension(3,3) :: crystallite_push33ToRef
  real(pReal), dimension(3,3), intent(in) :: tensor33
  real(pReal), dimension(3,3)             :: T
  integer, intent(in):: &
    el, &
    ip, &
    ipc
 
  T = matmul(material_orientation0(ipc,ip,el)%asMatrix(), &                                         ! ToDo: initial orientation correct?
             transpose(math_inv33(crystallite_subF(1:3,1:3,ipc,ip,el))))
  crystallite_push33ToRef = matmul(transpose(T),matmul(tensor33,T))

end function crystallite_push33ToRef


!--------------------------------------------------------------------------------------------------
!> @brief writes crystallite results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine crystallite_results
#if defined(PETSc) || defined(DAMASK_HDF5)
  integer :: p,o
  real(pReal),    allocatable, dimension(:,:,:) :: selected_tensors
  type(rotation), allocatable, dimension(:)     :: selected_rotations
  character(len=256) :: group,lattice_label
                                             
  do p=1,size(config_name_phase)
    group = trim('current/constituent')//'/'//trim(config_name_phase(p))//'/generic'
    
    call results_closeGroup(results_addGroup(group))  

    do o = 1, size(output_constituent(p)%label)
      select case (output_constituent(p)%label(o))
        case('f')
          selected_tensors = select_tensors(crystallite_partionedF,p)
          call results_writeDataset(group,selected_tensors,'F',&
                                   'deformation gradient','1')
        case('fe')
          selected_tensors = select_tensors(crystallite_Fe,p)
          call results_writeDataset(group,selected_tensors,'Fe',&
                                   'elastic deformation gradient','1')
        case('fp')
          selected_tensors = select_tensors(crystallite_Fp,p)
          call results_writeDataset(group,selected_tensors,'Fp',&
                                   'plastic deformation gradient','1')
        case('fi')
          selected_tensors = select_tensors(crystallite_Fi,p)
          call results_writeDataset(group,selected_tensors,'Fi',&
                                   'inelastic deformation gradient','1')
        case('lp')
          selected_tensors = select_tensors(crystallite_Lp,p)
          call results_writeDataset(group,selected_tensors,'Lp',&
                                   'plastic velocity gradient','1/s')
        case('li')
          selected_tensors = select_tensors(crystallite_Li,p)
          call results_writeDataset(group,selected_tensors,'Li',&
                                   'inelastic velocity gradient','1/s')
        case('p')
          selected_tensors = select_tensors(crystallite_P,p)
          call results_writeDataset(group,selected_tensors,'P',&
                                   '1st Piola-Kirchoff stress','Pa')
        case('s')
          selected_tensors = select_tensors(crystallite_S,p)
          call results_writeDataset(group,selected_tensors,'S',&
                                   '2nd Piola-Kirchoff stress','Pa')
        case('orientation')
          select case(lattice_structure(p))
            case(LATTICE_iso_ID)
              lattice_label = 'iso'
            case(LATTICE_fcc_ID)
              lattice_label = 'fcc'
            case(LATTICE_bcc_ID)
              lattice_label = 'bcc'
            case(LATTICE_bct_ID)
              lattice_label = 'bct'
            case(LATTICE_hex_ID)
              lattice_label = 'hex'
            case(LATTICE_ort_ID)
              lattice_label = 'ort'
          end select
          selected_rotations = select_rotations(crystallite_orientation,p)
          call results_writeDataset(group,selected_rotations,'orientation',&
                                   'crystal orientation as quaternion',lattice_label)
      end select
    enddo
  enddo
 
  contains

  !------------------------------------------------------------------------------------------------
  !> @brief select tensors for output
  !------------------------------------------------------------------------------------------------
  function select_tensors(dataset,instance)
 
    integer, intent(in) :: instance
    real(pReal), dimension(:,:,:,:,:), intent(in) :: dataset
    real(pReal), allocatable, dimension(:,:,:) :: select_tensors
    integer :: e,i,c,j
     
    allocate(select_tensors(3,3,count(material_phaseAt==instance)*homogenization_maxNgrains*discretization_nIP))

    j=0
    do e = 1, size(material_phaseAt,2)
      do i = 1, discretization_nIP
        do c = 1, size(material_phaseAt,1)                                                          !ToDo: this needs to be changed for varying Ngrains
          if (material_phaseAt(c,e) == instance) then
            j = j + 1
            select_tensors(1:3,1:3,j) = dataset(1:3,1:3,c,i,e)
          endif
        enddo
      enddo
    enddo
   
  end function select_tensors
 
 
!--------------------------------------------------------------------------------------------------
!> @brief select rotations for output
!-------------------------------------------------------------------------------------------------- 
  function select_rotations(dataset,instance)
 
    integer, intent(in) :: instance
    type(rotation), dimension(:,:,:), intent(in) :: dataset
    type(rotation), allocatable, dimension(:) :: select_rotations
    integer :: e,i,c,j
     
    allocate(select_rotations(count(material_phaseAt==instance)*homogenization_maxNgrains*discretization_nIP))

    j=0
    do e = 1, size(material_phaseAt,2)
      do i = 1, discretization_nIP
        do c = 1, size(material_phaseAt,1)                                                          !ToDo: this needs to be changed for varying Ngrains
           if (material_phaseAt(c,e) == instance) then
             j = j + 1
             select_rotations(j) = dataset(c,i,e)
           endif
        enddo
      enddo
   enddo
   
 end function select_rotations
#endif
end subroutine crystallite_results


!--------------------------------------------------------------------------------------------------
!> @brief calculation of stress (P) with time integration based on a residuum in Lp and
!> intermediate acceleration of the Newton-Raphson correction
!--------------------------------------------------------------------------------------------------
logical function integrateStress(ipc,ip,el,timeFraction)
 
  integer, intent(in)::         el, &                                                               ! element index
                                      ip, &                                                         ! integration point index
                                      ipc                                                           ! grain index
  real(pReal), optional, intent(in) :: timeFraction                                                 ! fraction of timestep
 
  real(pReal), dimension(3,3)::       Fg_new, &                                                     ! deformation gradient at end of timestep
                                      Fp_new, &                                                     ! plastic deformation gradient at end of timestep
                                      Fe_new, &                                                     ! elastic deformation gradient at end of timestep
                                      invFp_new, &                                                  ! inverse of Fp_new
                                      Fi_new, &                                                     ! gradient of intermediate deformation stages
                                      invFi_new, &
                                      invFp_current, &                                              ! inverse of Fp_current
                                      invFi_current, &                                              ! inverse of Fp_current
                                      Lpguess, &                                                    ! current guess for plastic velocity gradient
                                      Lpguess_old, &                                                ! known last good guess for plastic velocity gradient
                                      Lp_constitutive, &                                            ! plastic velocity gradient resulting from constitutive law
                                      residuumLp, &                                                 ! current residuum of plastic velocity gradient
                                      residuumLp_old, &                                             ! last residuum of plastic velocity gradient
                                      deltaLp, &                                                    ! direction of next guess
                                      Liguess, &                                                    ! current guess for intermediate velocity gradient
                                      Liguess_old, &                                                ! known last good guess for intermediate velocity gradient
                                      Li_constitutive, &                                            ! intermediate velocity gradient resulting from constitutive law
                                      residuumLi, &                                                 ! current residuum of intermediate velocity gradient
                                      residuumLi_old, &                                             ! last residuum of intermediate velocity gradient
                                      deltaLi, &                                                    ! direction of next guess
                                      S, &                                                          ! 2nd Piola-Kirchhoff Stress in plastic (lattice) configuration
                                      A, &
                                      B, &
                                      Fe, &                                                         ! elastic deformation gradient
                                      temp_33
  real(pReal), dimension(9) ::        work                                                          ! needed for matrix inversion by LAPACK
  integer,    dimension(9) ::         devNull                                                       ! needed for matrix inversion by LAPACK
  real(pReal), dimension(9,9) ::      dRLp_dLp, &                                                   ! partial derivative of residuum (Jacobian for Newton-Raphson scheme)
                                      dRLp_dLp2, &                                                  ! working copy of dRdLp
                                      dRLi_dLi                                                      ! partial derivative of residuumI (Jacobian for Newton-Raphson scheme)
  real(pReal), dimension(3,3,3,3)::   dS_dFe, &                                                     ! partial derivative of 2nd Piola-Kirchhoff stress
                                      dS_dFi, &
                                      dFe_dLp, &                                                    ! partial derivative of elastic deformation gradient
                                      dFe_dLi, &
                                      dFi_dLi, &
                                      dLp_dFi, &
                                      dLi_dFi, &
                                      dLp_dS, &
                                      dLi_dS
  real(pReal)                         detInvFi, &                                                   ! determinant of InvFi
                                      steplengthLp, &
                                      steplengthLi, &
                                      dt, &                                                         ! time increment
                                      aTolLp, &
                                      aTolLi
  integer                             NiterationStressLp, &                                         ! number of stress integrations
                                      NiterationStressLi, &                                         ! number of inner stress integrations
                                      ierr, &                                                       ! error indicator for LAPACK
                                      o, &
                                      p, &
                                      jacoCounterLp, &
                                      jacoCounterLi                                                 ! counters to check for Jacobian update
  external :: &
    dgesv
 
  !* be pessimistic
  integrateStress = .false.
#ifdef DEBUG
  if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
      .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
             .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) &
  write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST integrateStress >> at el ip ipc ',el,ip,ipc
#endif

  if (present(timeFraction)) then
    dt = crystallite_subdt(ipc,ip,el) * timeFraction
    Fg_new = crystallite_subF0(1:3,1:3,ipc,ip,el) &
           + (crystallite_subF(1:3,1:3,ipc,ip,el) - crystallite_subF0(1:3,1:3,ipc,ip,el)) * timeFraction
  else
    dt = crystallite_subdt(ipc,ip,el)
    Fg_new = crystallite_subF(1:3,1:3,ipc,ip,el)
  endif

  Lpguess     =   crystallite_Lp(1:3,1:3,ipc,ip,el)                                                 ! take as first guess
  Liguess     =   crystallite_Li(1:3,1:3,ipc,ip,el)                                                 ! take as first guess
  Liguess_old =   Liguess
 
  invFp_current = math_inv33(crystallite_subFp0(1:3,1:3,ipc,ip,el))
  failedInversionFp: if (all(dEq0(invFp_current))) then
#ifdef DEBUG
    if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
      write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST integrateStress >> failed on inversion of current Fp at el ip ipc ',&
                                     el,ip,ipc
    if (iand(debug_level(debug_crystallite), debug_levelExtensive) > 0) &
      write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> current Fp ',transpose(crystallite_subFp0(1:3,1:3,ipc,ip,el))
#endif
    return
  endif failedInversionFp
  A = matmul(Fg_new,invFp_current)                                                            ! intermediate tensor needed later to calculate dFe_dLp

  invFi_current = math_inv33(crystallite_subFi0(1:3,1:3,ipc,ip,el))
  failedInversionFi: if (all(dEq0(invFi_current))) then
#ifdef DEBUG
    if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
      write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST integrateStress >> failed on inversion of current Fi at el ip ipc ',&
                                     el,ip,ipc
    if (iand(debug_level(debug_crystallite), debug_levelExtensive) > 0) &
      write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> current Fi ', &
                                             transpose(crystallite_subFi0(1:3,1:3,ipc,ip,el))
#endif
    return
  endif failedInversionFi
 
  !* start Li loop with normal step length
  NiterationStressLi = 0
  jacoCounterLi      = 0
  steplengthLi       = 1.0_pReal
  residuumLi_old     = 0.0_pReal
 
  LiLoop: do
    NiterationStressLi = NiterationStressLi + 1
    LiLoopLimit: if (NiterationStressLi > num%nStress) then
#ifdef DEBUG
      if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
        write(6,'(a,i3,a,i8,1x,i2,1x,i3,/)') '<< CRYST integrateStress >> reached Li loop limit',num%nStress, &
                                            ' at el ip ipc ', el,ip,ipc
#endif
      return
    endif LiLoopLimit
 
    invFi_new = matmul(invFi_current,math_I3 - dt*Liguess)
    Fi_new    = math_inv33(invFi_new)
    detInvFi  = math_det33(invFi_new)
 
    !* start Lp loop with normal step length
    NiterationStressLp = 0
    jacoCounterLp      = 0
    steplengthLp       = 1.0_pReal
    residuumLp_old     = 0.0_pReal
    Lpguess_old        = Lpguess
 
    LpLoop: do
      NiterationStressLp = NiterationStressLp + 1
      LpLoopLimit: if (NiterationStressLp > num%nStress) then
#ifdef DEBUG
        if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
           write(6,'(a,i3,a,i8,1x,i2,1x,i3,/)') '<< CRYST integrateStress >> reached Lp loop limit',num%nStress, &
                                                ' at el ip ipc ', el,ip,ipc
#endif
        return
      endif LpLoopLimit
 
      !* calculate (elastic) 2nd Piola--Kirchhoff stress tensor and its tangent from constitutive law
 
      B  = math_I3 - dt*Lpguess
      Fe = matmul(matmul(A,B), invFi_new)
      call constitutive_SandItsTangents(S, dS_dFe, dS_dFi, &
                                        Fe, Fi_new, ipc, ip, el)                                     ! call constitutive law to calculate 2nd Piola-Kirchhoff stress and its derivative in unloaded configuration
 
      !* calculate plastic velocity gradient and its tangent from constitutive law
      call constitutive_LpAndItsTangents(Lp_constitutive, dLp_dS, dLp_dFi, &
                                         S, Fi_new, ipc, ip, el)
 
#ifdef DEBUG
      if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
          .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                 .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
        write(6,'(a,i3,/)')                   '<< CRYST integrateStress >> Lp iteration ', NiterationStressLp
        write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST integrateStress >> Lpguess', transpose(Lpguess)
        write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST integrateStress >> Lp_constitutive', transpose(Lp_constitutive)
        write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST integrateStress >> Fi', transpose(Fi_new)
        write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST integrateStress >> Fe', transpose(Fe)
        write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST integrateStress >> S', transpose(S)
      endif
#endif

      !* update current residuum and check for convergence of loop
      aTolLp = max(num%rTol_crystalliteStress * max(norm2(Lpguess),norm2(Lp_constitutive)), &       ! absolute tolerance from largest acceptable relative error
                   num%aTol_crystalliteStress)                                                      ! minimum lower cutoff
      residuumLp = Lpguess - Lp_constitutive
 
      if (any(IEEE_is_NaN(residuumLp))) then
#ifdef DEBUG
        if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
          write(6,'(a,i8,1x,i2,1x,i3,a,i3,a)') '<< CRYST integrateStress >> encountered NaN for Lp-residuum at el ip ipc ', &
                                               el,ip,ipc, &
                                               ' ; iteration ', NiterationStressLp,&
                                               ' >> returning..!'
#endif
        return                                                                                      ! ...me = .false. to inform integrator about problem
      elseif (norm2(residuumLp) < aTolLp) then                                                      ! converged if below absolute tolerance
        exit LpLoop                                                                                 ! ...leave iteration loop
      elseif (     NiterationStressLp == 1 &
              .or. norm2(residuumLp) < norm2(residuumLp_old)) then                                  ! not converged, but improved norm of residuum (always proceed in first iteration)...
        residuumLp_old = residuumLp                                                                 ! ...remember old values and...
        Lpguess_old    = Lpguess
        steplengthLp   = 1.0_pReal                                                                  ! ...proceed with normal step length (calculate new search direction)
      else                                                                                          ! not converged and residuum not improved...
        steplengthLp = num%subStepSizeLp * steplengthLp                                             ! ...try with smaller step length in same direction
        Lpguess      = Lpguess_old + steplengthLp * deltaLp
#ifdef DEBUG
        if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
            .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                   .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
          write(6,'(a,1x,f7.4)') '<< CRYST integrateStress >> linear search for Lpguess with step', steplengthLp
        endif
#endif
        cycle LpLoop
      endif


     !* calculate Jacobian for correction term
      if (mod(jacoCounterLp, num%iJacoLpresiduum) == 0) then
        do o=1,3; do p=1,3
          dFe_dLp(o,1:3,p,1:3) = A(o,p)*transpose(invFi_new)                                        ! dFe_dLp(i,j,k,l) = -dt * A(i,k) invFi(l,j)
        enddo; enddo
        dFe_dLp  = - dt * dFe_dLp
        dRLp_dLp = math_identity2nd(9) &
                 - math_3333to99(math_mul3333xx3333(math_mul3333xx3333(dLp_dS,dS_dFe),dFe_dLp))
#ifdef DEBUG
        if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
            .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                   .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
          write(6,'(a,/,9(12x,9(e12.4,1x)/))') '<< CRYST integrateStress >> dLp_dS', math_3333to99(dLp_dS)
          write(6,'(a,1x,e20.10)')             '<< CRYST integrateStress >> dLp_dS norm', norm2(math_3333to99(dLp_dS))
          write(6,'(a,/,9(12x,9(e12.4,1x)/))') '<< CRYST integrateStress >> dRLp_dLp', dRLp_dLp-math_identity2nd(9)
          write(6,'(a,1x,e20.10)')             '<< CRYST integrateStress >> dRLp_dLp norm', norm2(dRLp_dLp-math_identity2nd(9))
        endif
#endif
        dRLp_dLp2 = dRLp_dLp                                                                        ! will be overwritten in first call to LAPACK routine
        work = math_33to9(residuumLp)
        call dgesv(9,1,dRLp_dLp2,9,devNull,work,9,ierr)                                             ! solve dRLp/dLp * delta Lp = -res for delta Lp
        if (ierr /= 0) then
#ifdef DEBUG
          if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) then
           write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST integrateStress >> failed on dR/dLp inversion at el ip ipc ', &
                                         el,ip,ipc
            if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
                .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g)&
                       .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
              write(6,*)
              write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dR_dLp',transpose(dRLp_dLp)
              write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dFe_dLp',transpose(math_3333to99(dFe_dLp))
              write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dS_dFe (cnst)',transpose(math_3333to99(dS_dFe))
              write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dLp_dS (cnst)',transpose(math_3333to99(dLp_dS))
              write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> A',transpose(A)
              write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> B',transpose(B)
              write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> Lp_constitutive',transpose(Lp_constitutive)
              write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> Lpguess',transpose(Lpguess)
            endif
          endif
#endif
          return
        endif
        deltaLp = - math_9to33(work)
      endif
      jacoCounterLp = jacoCounterLp + 1

      Lpguess = Lpguess + steplengthLp * deltaLp

    enddo LpLoop

   !* calculate intermediate velocity gradient and its tangent from constitutive law
    call constitutive_LiAndItsTangents(Li_constitutive, dLi_dS, dLi_dFi, &
                                       S, Fi_new, ipc, ip, el)

#ifdef DEBUG
      if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
          .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                 .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
        write(6,'(a,i3,/)')                  '<< CRYST integrateStress >> Li iteration ', NiterationStressLi
        write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> Li_constitutive', transpose(Li_constitutive)
        write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> Liguess', transpose(Liguess)
      endif
#endif

    !* update current residuum and check for convergence of loop
    aTolLi = max(num%rTol_crystalliteStress * max(norm2(Liguess),norm2(Li_constitutive)), &         ! absolute tolerance from largest acceptable relative error
                 num%aTol_crystalliteStress)                                                        ! minimum lower cutoff
    residuumLi = Liguess - Li_constitutive
    if (any(IEEE_is_NaN(residuumLi))) then                                                          ! NaN in residuum...
#ifdef DEBUG
      if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) &
        write(6,'(a,i8,1x,i2,1x,i3,a,i3,a)') '<< CRYST integrateStress >> encountered NaN for Li-residuum at el ip ipc ', &
                                             el,ip,ipc, &
                                             ' ; iteration ', NiterationStressLi,&
                                             ' >> returning..!'
#endif
      return                                                                                        ! ...me = .false. to inform integrator about problem
    elseif (norm2(residuumLi) < aTolLi) then                                                        ! converged if below absolute tolerance
      exit LiLoop                                                                                   ! ...leave iteration loop
    elseif (     NiterationStressLi == 1 &
            .or. norm2(residuumLi) < norm2(residuumLi_old)) then                                    ! not converged, but improved norm of residuum (always proceed in first iteration)...
      residuumLi_old = residuumLi                                                                   ! ...remember old values and...
      Liguess_old    = Liguess
      steplengthLi   = 1.0_pReal                                                                    ! ...proceed with normal step length (calculate new search direction)
    else                                                                                            ! not converged and residuum not improved...
      steplengthLi   = num%subStepSizeLi * steplengthLi                                             ! ...try with smaller step length in same direction
      Liguess        = Liguess_old + steplengthLi * deltaLi
#ifdef DEBUG
        if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
            .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                   .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
          write(6,'(a,1x,f7.4)') '<< CRYST integrateStress >> linear search for Liguess with step', steplengthLi
        endif
#endif
      cycle LiLoop
    endif
 
    !* calculate Jacobian for correction term
    if (mod(jacoCounterLi, num%iJacoLpresiduum) == 0) then
      temp_33     = matmul(matmul(A,B),invFi_current)
      do o=1,3; do p=1,3
        dFe_dLi(1:3,o,1:3,p) = -dt*math_I3(o,p)*temp_33                                             ! dFe_dLp(i,j,k,l) = -dt * A(i,k) invFi(l,j)
        dFi_dLi(1:3,o,1:3,p) = -dt*math_I3(o,p)*invFi_current
      enddo; enddo
      do o=1,3; do p=1,3
        dFi_dLi(1:3,1:3,o,p) = matmul(matmul(Fi_new,dFi_dLi(1:3,1:3,o,p)),Fi_new)
      enddo; enddo
      dRLi_dLi  = math_identity2nd(9) &
                - math_3333to99(math_mul3333xx3333(dLi_dS, math_mul3333xx3333(dS_dFe, dFe_dLi) + &
                                                  math_mul3333xx3333(dS_dFi, dFi_dLi)))  &
                - math_3333to99(math_mul3333xx3333(dLi_dFi, dFi_dLi))
      work = math_33to9(residuumLi)
      call dgesv(9,1,dRLi_dLi,9,devNull,work,9,ierr)                                                ! solve dRLi/dLp * delta Li = -res for delta Li
      if (ierr /= 0) then
#ifdef DEBUG
        if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) then
          write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST integrateStress >> failed on dR/dLi inversion at el ip ipc ', &
                                        el,ip,ipc
          if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
              .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g)&
                     .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
            write(6,*)
            write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dR_dLi',transpose(dRLi_dLi)
            write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dFe_dLi',transpose(math_3333to99(dFe_dLi))
            write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dS_dFi (cnst)',transpose(math_3333to99(dS_dFi))
            write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST integrateStress >> dLi_dS (cnst)',transpose(math_3333to99(dLi_dS))
            write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> Li_constitutive',transpose(Li_constitutive)
            write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> Liguess',transpose(Liguess)
          endif
        endif
#endif
        return
      endif
 
      deltaLi = - math_9to33(work)
    endif
    jacoCounterLi = jacoCounterLi + 1
 
    Liguess = Liguess + steplengthLi * deltaLi
#ifdef DEBUG
        if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
            .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                   .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
          write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST integrateStress >> corrected Liguess by', transpose(deltaLi)
        endif
#endif
  enddo LiLoop
 
  !* calculate new plastic and elastic deformation gradient
  invFp_new = matmul(invFp_current,B)
  invFp_new = invFp_new / math_det33(invFp_new)**(1.0_pReal/3.0_pReal)                              ! regularize
  Fp_new = math_inv33(invFp_new)
  failedInversionInvFp: if (all(dEq0(Fp_new))) then
#ifdef DEBUG
    if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0) then
     write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST integrateStress >> failed on invFp_new inversion at el ip ipc ', &
                                   el,ip,ipc
      if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
          .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                 .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) &
        write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> invFp_new',transpose(invFp_new)
    endif
#endif
    return
  endif failedInversionInvFp
  Fe_new = matmul(matmul(Fg_new,invFp_new),invFi_new)

!--------------------------------------------------------------------------------------------------
! stress integration was successful
  integrateStress = .true.
  crystallite_P    (1:3,1:3,ipc,ip,el) = matmul(matmul(Fg_new,invFp_new),matmul(S,transpose(invFp_new)))
  crystallite_S    (1:3,1:3,ipc,ip,el) = S
  crystallite_Lp   (1:3,1:3,ipc,ip,el) = Lpguess
  crystallite_Li   (1:3,1:3,ipc,ip,el) = Liguess
  crystallite_Fp   (1:3,1:3,ipc,ip,el) = Fp_new
  crystallite_Fi   (1:3,1:3,ipc,ip,el) = Fi_new
  crystallite_Fe   (1:3,1:3,ipc,ip,el) = Fe_new
  crystallite_invFp(1:3,1:3,ipc,ip,el) = invFp_new
  crystallite_invFi(1:3,1:3,ipc,ip,el) = invFi_new

#ifdef DEBUG
  if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0 &
      .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
              .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
    write(6,'(a,/)')                     '<< CRYST integrateStress >> successful integration'
    write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> P / MPa', &
               transpose(crystallite_P(1:3,1:3,ipc,ip,el))*1.0e-6_pReal
    write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> Cauchy / MPa', &
               matmul(crystallite_P(1:3,1:3,ipc,ip,el), transpose(Fg_new)) * 1.0e-6_pReal / math_det33(Fg_new)
    write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> Fe Lp Fe^-1', &
               transpose(matmul(Fe_new, matmul(crystallite_Lp(1:3,1:3,ipc,ip,el), math_inv33(Fe_new))))
    write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> Fp',transpose(crystallite_Fp(1:3,1:3,ipc,ip,el))
    write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST integrateStress >> Fi',transpose(crystallite_Fi(1:3,1:3,ipc,ip,el))
  endif
#endif

end function integrateStress


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with adaptive 1st order explicit Euler method
!> using Fixed Point Iteration to adapt the stepsize
!--------------------------------------------------------------------------------------------------
subroutine integrateStateFPI

 integer :: &
   NiterationState, &                                                                               !< number of iterations in state loop
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g, &                                                                                             !< grain index in grain loop
   p, &
   c, &
   s, &
   sizeDotState
 real(pReal) :: &
   zeta
 real(pReal), dimension(constitutive_plasticity_maxSizeDotState) :: &
   residuum_plastic                                                                                 ! residuum for plastic state
 real(pReal), dimension(constitutive_source_maxSizeDotState) :: &
   residuum_source                                                                                  ! residuum for source state
 logical :: &
   doneWithIntegration

 ! --+>> PREGUESS FOR STATE <<+--
 call update_dotState(1.0_pReal)
 call update_state(1.0_pReal)

 NiterationState = 0
 doneWithIntegration = .false.
 crystalliteLooping: do while (.not. doneWithIntegration .and. NiterationState < num%nState)
   NiterationState = NiterationState + 1

#ifdef DEBUG
   if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0) &
     write(6,'(a,i6)') '<< CRYST stateFPI >> state iteration ',NiterationState
#endif

   ! store previousDotState and previousDotState2
   
   !$OMP PARALLEL DO PRIVATE(p,c)
     do e = FEsolving_execElem(1),FEsolving_execElem(2)
       do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
         do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
           if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
             p = material_phaseAt(g,e); c = material_phaseMemberAt(g,i,e)

             plasticState(p)%previousDotState2(:,c) = merge(plasticState(p)%previousDotState(:,c),&
                                                            0.0_pReal,&
                                                            NiterationState > 1)
             plasticState(p)%previousDotState (:,c) = plasticState(p)%dotState(:,c)
             do s = 1, phase_Nsources(p)
               sourceState(p)%p(s)%previousDotState2(:,c) = merge(sourceState(p)%p(s)%previousDotState(:,c),&
                                                                  0.0_pReal, &
                                                                  NiterationState > 1)
               sourceState(p)%p(s)%previousDotState (:,c) = sourceState(p)%p(s)%dotState(:,c)
             enddo
           endif
       enddo
     enddo
   enddo
   !$OMP END PARALLEL DO

   call update_dependentState
   call update_stress(1.0_pReal)
   call update_dotState(1.0_pReal)
   
   !$OMP PARALLEL
   !$OMP DO PRIVATE(sizeDotState,residuum_plastic,residuum_source,zeta,p,c)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
         if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
           p = material_phaseAt(g,e); c = material_phaseMemberAt(g,i,e)
           sizeDotState = plasticState(p)%sizeDotState

           zeta = damper(plasticState(p)%dotState         (:,c), &
                         plasticState(p)%previousDotState (:,c), &
                         plasticState(p)%previousDotState2(:,c))
          
           residuum_plastic(1:SizeDotState) = plasticState(p)%state    (1:sizeDotState,c) &
                                            - plasticState(p)%subState0(1:sizeDotState,c)  &
                                            - (  plasticState(p)%dotState        (:,c) * zeta &
                                               + plasticState(p)%previousDotState(:,c) * (1.0_pReal-zeta) &
                                              ) * crystallite_subdt(g,i,e)

           plasticState(p)%state(1:sizeDotState,c) = plasticState(p)%state(1:sizeDotState,c) &
                                                   - residuum_plastic(1:sizeDotState) 
           plasticState(p)%dotState(:,c) = plasticState(p)%dotState(:,c) * zeta &
                                         + plasticState(p)%previousDotState(:,c) * (1.0_pReal - zeta)
           
           crystallite_converged(g,i,e) = converged(residuum_plastic(1:sizeDotState), &
                                                    plasticState(p)%state(1:sizeDotState,c), &
                                                    plasticState(p)%aTolState(1:sizeDotState))
                             

           do s = 1, phase_Nsources(p)
             sizeDotState  = sourceState(p)%p(s)%sizeDotState
             
             zeta = damper(sourceState(p)%p(s)%dotState         (:,c), &
                           sourceState(p)%p(s)%previousDotState (:,c), &
                           sourceState(p)%p(s)%previousDotState2(:,c))

             residuum_source(1:sizeDotState) = sourceState(p)%p(s)%state    (1:sizeDotState,c)  &
                                             - sourceState(p)%p(s)%subState0(1:sizeDotState,c)  &
                                             - (  sourceState(p)%p(s)%dotState         (:,c) * zeta &
                                                 + sourceState(p)%p(s)%previousDotState(:,c) * (1.0_pReal - zeta) &
                                               ) * crystallite_subdt(g,i,e)

             sourceState(p)%p(s)%state(1:sizeDotState,c) = sourceState(p)%p(s)%state(1:sizeDotState,c) &
                                                         - residuum_source(1:sizeDotState)
             sourceState(p)%p(s)%dotState(:,c) = sourceState(p)%p(s)%dotState(:,c) * zeta &
                                               + sourceState(p)%p(s)%previousDotState(:,c)* (1.0_pReal - zeta)

             crystallite_converged(g,i,e) = &
             crystallite_converged(g,i,e) .and. converged(residuum_source(1:sizeDotState), &
                                                          sourceState(p)%p(s)%state(1:sizeDotState,c), &
                                                          sourceState(p)%p(s)%aTolState(1:sizeDotState))
           enddo
         endif
   enddo; enddo; enddo
   !$OMP ENDDO

   !$OMP DO
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. crystallite_converged(g,i,e)) then                                  ! converged and still alive...
         crystallite_todo(g,i,e) = stateJump(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e)) then                                                             ! if state jump fails, then convergence is broken
           crystallite_converged(g,i,e) = .false.
           if (.not. crystallite_localPlasticity(g,i,e)) then                                                ! if broken non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP END PARALLEL


   if (any(plasticState(:)%nonlocal)) call nonlocalConvergenceCheck


   ! --- CHECK IF DONE WITH INTEGRATION ---
   doneWithIntegration = .true.
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         doneWithIntegration = .false.
         exit
       endif
     enddo; enddo
   enddo

 enddo crystalliteLooping


 contains

 !--------------------------------------------------------------------------------------------------
 !> @brief calculate the damping for correction of state and dot state
 !--------------------------------------------------------------------------------------------------
 real(pReal) pure function damper(current,previous,previous2)
 
 real(pReal), dimension(:), intent(in) ::&
   current, previous, previous2
 
 real(pReal) :: dot_prod12, dot_prod22
   
 dot_prod12 = dot_product(current  - previous,  previous - previous2)
 dot_prod22 = dot_product(previous - previous2, previous - previous2)
 if ((dot_product(current,previous) < 0.0_pReal .or. dot_prod12 < 0.0_pReal) .and. dot_prod22 > 0.0_pReal) then
   damper = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)
 else
   damper = 1.0_pReal
 endif
   
 end function damper

end subroutine integrateStateFPI


!--------------------------------------------------------------------------------------------------
!> @brief integrate state with 1st order explicit Euler method
!--------------------------------------------------------------------------------------------------
subroutine integrateStateEuler

 call update_dotState(1.0_pReal)
 call update_state(1.0_pReal)
 call update_deltaState
 call update_dependentState
 call update_stress(1.0_pReal)
 call setConvergenceFlag
 if (any(plasticState(:)%nonlocal)) call nonlocalConvergenceCheck

end subroutine integrateStateEuler


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with 1st order Euler method with adaptive step size
!--------------------------------------------------------------------------------------------------
subroutine integrateStateAdaptiveEuler

 integer :: &
   e, &                                                                                             ! element index in element loop
   i, &                                                                                             ! integration point index in ip loop
   g, &                                                                                             ! grain index in grain loop
   p, &
   c, &
   s, &
   sizeDotState
   
  ! ToDo: MD: once all constitutives use allocate state, attach residuum arrays to the state in case of adaptive Euler
 real(pReal), dimension(constitutive_plasticity_maxSizeDotState,            &
                        homogenization_maxNgrains,discretization_nIP,discretization_nElem) :: &
   residuum_plastic
 real(pReal), dimension(constitutive_source_maxSizeDotState,&
                        maxval(phase_Nsources), &
                        homogenization_maxNgrains,discretization_nIP,discretization_nElem) :: &
   residuum_source

!--------------------------------------------------------------------------------------------------
! contribution to state and relative residui and from Euler integration
 call update_dotState(1.0_pReal)

 !$OMP PARALLEL DO PRIVATE(sizeDotState,p,c)
 do e = FEsolving_execElem(1),FEsolving_execElem(2)
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
       if (crystallite_todo(g,i,e)) then
         p = material_phaseAt(g,e); c = material_phaseMemberAt(g,i,e)
         sizeDotState = plasticState(p)%sizeDotState
         
         residuum_plastic(1:sizeDotState,g,i,e) = plasticState(p)%dotstate(1:sizeDotState,c) &
                                                * (- 0.5_pReal * crystallite_subdt(g,i,e))
         plasticState(p)%state(1:sizeDotState,c) = &
         plasticState(p)%state(1:sizeDotState,c) + plasticState(p)%dotstate(1:sizeDotState,c) * crystallite_subdt(g,i,e) !ToDo: state, partitioned state?
         do s = 1, phase_Nsources(p)
           sizeDotState = sourceState(p)%p(s)%sizeDotState
           
           residuum_source(1:sizeDotState,s,g,i,e) = sourceState(p)%p(s)%dotstate(1:sizeDotState,c) &
                                                   * (- 0.5_pReal * crystallite_subdt(g,i,e))
           sourceState(p)%p(s)%state(1:sizeDotState,c) = &
           sourceState(p)%p(s)%state(1:sizeDotState,c) + sourceState(p)%p(s)%dotstate(1:sizeDotState,c) * crystallite_subdt(g,i,e) !ToDo: state, partitioned state?
         enddo
       endif
     enddo; enddo; enddo
 !$OMP END PARALLEL DO

 call update_deltaState
 call update_dependentState
 call update_stress(1.0_pReal)
 call update_dotState(1.0_pReal)

 !$OMP PARALLEL DO PRIVATE(sizeDotState,p,c)
 do e = FEsolving_execElem(1),FEsolving_execElem(2)
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
       if (crystallite_todo(g,i,e)) then
         p = material_phaseAt(g,e); c = material_phaseMemberAt(g,i,e)
         sizeDotState = plasticState(p)%sizeDotState
         
         residuum_plastic(1:sizeDotState,g,i,e) = residuum_plastic(1:sizeDotState,g,i,e) &
                                                + 0.5_pReal * plasticState(p)%dotState(:,c) * crystallite_subdt(g,i,e)
              
         crystallite_converged(g,i,e) = converged(residuum_plastic(1:sizeDotState,g,i,e), &
                                                  plasticState(p)%state(1:sizeDotState,c), &
                                                  plasticState(p)%aTolState(1:sizeDotState))

         do s = 1, phase_Nsources(p)
           sizeDotState = sourceState(p)%p(s)%sizeDotState
           
           residuum_source(1:sizeDotState,s,g,i,e) = &
           residuum_source(1:sizeDotState,s,g,i,e) + 0.5_pReal * sourceState(p)%p(s)%dotState(:,c) * crystallite_subdt(g,i,e)

           crystallite_converged(g,i,e) = &
           crystallite_converged(g,i,e) .and. converged(residuum_source(1:sizeDotState,s,g,i,e), &
                                                        sourceState(p)%p(s)%state(1:sizeDotState,c), &
                                                        sourceState(p)%p(s)%aTolState(1:sizeDotState))
          enddo
          
       endif
 enddo; enddo; enddo
 !$OMP END PARALLEL DO

 if (any(plasticState(:)%nonlocal)) call nonlocalConvergenceCheck
 
end subroutine integrateStateAdaptiveEuler


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with 4th order explicit Runge Kutta method
! ToDo: This is totally BROKEN: RK4dotState is never used!!!
!--------------------------------------------------------------------------------------------------
subroutine integrateStateRK4

 real(pReal), dimension(4), parameter :: &
   TIMESTEPFRACTION = [0.5_pReal, 0.5_pReal, 1.0_pReal, 1.0_pReal]                                   ! factor giving the fraction of the original timestep used for Runge Kutta Integration
 real(pReal), dimension(4), parameter :: &
   WEIGHT = [1.0_pReal, 2.0_pReal, 2.0_pReal, 1.0_pReal/6.0_pReal]                                   ! weight of slope used for Runge Kutta integration (final weight divided by 6)

 integer ::                                    e, &                                                  ! element index in element loop
                                               i, &                                                  ! integration point index in ip loop
                                               g, &                                                  ! grain index in grain loop
                                               p, &                                                  ! phase loop
                                               c, &
                                               n, &
                                               s

 call update_dotState(1.0_pReal)


 do n = 1,4

   !$OMP PARALLEL DO PRIVATE(p,c)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
       if (crystallite_todo(g,i,e)) then
         p = material_phaseAt(g,e); c = material_phaseMemberAt(g,i,e)

         plasticState(p)%RK4dotState(:,c) = WEIGHT(n)*plasticState(p)%dotState(:,c) &
                                          + merge(plasticState(p)%RK4dotState(:,c),0.0_pReal,n>1)
         do s = 1, phase_Nsources(p)
           sourceState(p)%p(s)%RK4dotState(:,c) = WEIGHT(n)*sourceState(p)%p(s)%dotState(:,c) &
                                                + merge(sourceState(p)%p(s)%RK4dotState(:,c),0.0_pReal,n>1)
         enddo
       endif
     enddo; enddo; enddo
   !$OMP END PARALLEL DO

   call update_state(TIMESTEPFRACTION(n))
   call update_deltaState
   call update_dependentState
   call update_stress(TIMESTEPFRACTION(n))
   ! --- dot state and RK dot state---

   first3steps: if (n < 4) then
     call update_dotState(TIMESTEPFRACTION(n))
   endif first3steps

 enddo

 call setConvergenceFlag
 if (any(plasticState(:)%nonlocal)) call nonlocalConvergenceCheck

end subroutine integrateStateRK4


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with 5th order Runge-Kutta Cash-Karp method with
!> adaptive step size  (use 5th order solution to advance = "local extrapolation")
!--------------------------------------------------------------------------------------------------
subroutine integrateStateRKCK45

 real(pReal), dimension(5,5), parameter :: &
   A = reshape([&
     .2_pReal, .075_pReal,   .3_pReal, -11.0_pReal/54.0_pReal,  1631.0_pReal/55296.0_pReal, &
     .0_pReal, .225_pReal,  -.9_pReal,   2.5_pReal,              175.0_pReal/512.0_pReal, &
     .0_pReal,   .0_pReal,  1.2_pReal, -70.0_pReal/27.0_pReal,   575.0_pReal/13824.0_pReal, &
     .0_pReal,   .0_pReal,   .0_pReal,  35.0_pReal/27.0_pReal, 44275.0_pReal/110592.0_pReal, &
     .0_pReal,   .0_pReal,   .0_pReal,    .0_pReal,              253.0_pReal/4096.0_pReal], &
     [5,5], order=[2,1])                                                                            !< coefficients in Butcher tableau (used for preliminary integration in stages 2 to 6)

 real(pReal), dimension(6), parameter :: &
   B = &
     [37.0_pReal/378.0_pReal, .0_pReal, 250.0_pReal/621.0_pReal, &
     125.0_pReal/594.0_pReal, .0_pReal, 512.0_pReal/1771.0_pReal], &                                !< coefficients in Butcher tableau (used for final integration and error estimate)
   DB = B - &
     [2825.0_pReal/27648.0_pReal,    .0_pReal,                18575.0_pReal/48384.0_pReal,&
     13525.0_pReal/55296.0_pReal, 277.0_pReal/14336.0_pReal,      0.25_pReal]                       !< coefficients in Butcher tableau (used for final integration and error estimate)

 real(pReal), dimension(5), parameter :: &
   C = [0.2_pReal, 0.3_pReal, 0.6_pReal, 1.0_pReal, 0.875_pReal]                                    !< coefficients in Butcher tableau (fractions of original time step in stages 2 to 6)

 integer :: &
   e, &                                                                                             ! element index in element loop
   i, &                                                                                             ! integration point index in ip loop
   g, &                                                                                             ! grain index in grain loop
   stage, &                                                                                         ! stage index in integration stage loop
   n, &
   p, &
   cc, &
   s, &
   sizeDotState

   ! ToDo: MD: once all constitutives use allocate state, attach residuum arrays to the state in case of RKCK45

 real(pReal), dimension(constitutive_plasticity_maxSizeDotState,            &
                        homogenization_maxNgrains,discretization_nIP,discretization_nElem) :: &
   residuum_plastic                                                                         ! relative residuum from evolution in microstructure
 real(pReal), dimension(constitutive_source_maxSizeDotState, &
                        maxval(phase_Nsources), &
                        homogenization_maxNgrains,discretization_nIP,discretization_nElem) :: &
   residuum_source                                                                   ! relative residuum from evolution in microstructure


 call update_dotState(1.0_pReal)

 ! --- SECOND TO SIXTH RUNGE KUTTA STEP ---

 do stage = 1,5

   ! --- state update ---

   !$OMP PARALLEL DO PRIVATE(p,cc)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
       if (crystallite_todo(g,i,e)) then
         p = material_phaseAt(g,e); cc = material_phaseMemberAt(g,i,e)

         plasticState(p)%RKCK45dotState(stage,:,cc) = plasticState(p)%dotState(:,cc)
         plasticState(p)%dotState(:,cc) = A(1,stage) * plasticState(p)%RKCK45dotState(1,:,cc)

         do s = 1, phase_Nsources(p)
           sourceState(p)%p(s)%RKCK45dotState(stage,:,cc) = sourceState(p)%p(s)%dotState(:,cc)
           sourceState(p)%p(s)%dotState(:,cc) = A(1,stage) * sourceState(p)%p(s)%RKCK45dotState(1,:,cc)
         enddo

         do n = 2, stage
           plasticState(p)%dotState(:,cc) = plasticState(p)%dotState(:,cc) &
                                          + A(n,stage) * plasticState(p)%RKCK45dotState(n,:,cc)
           do s = 1, phase_Nsources(p)
             sourceState(p)%p(s)%dotState(:,cc) = sourceState(p)%p(s)%dotState(:,cc) &
                                                + A(n,stage) * sourceState(p)%p(s)%RKCK45dotState(n,:,cc)
           enddo
         enddo

       endif
     enddo; enddo; enddo
   !$OMP END PARALLEL DO

    call update_state(1.0_pReal) !MD: 1.0 correct?
    call update_deltaState
    call update_dependentState
    call update_stress(C(stage))
    call update_dotState(C(stage))

 enddo


!--------------------------------------------------------------------------------------------------
! --- STATE UPDATE WITH ERROR ESTIMATE FOR STATE ---

 !$OMP PARALLEL DO PRIVATE(sizeDotState,p,cc)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
     if (crystallite_todo(g,i,e)) then
       p = material_phaseAt(g,e); cc = material_phaseMemberAt(g,i,e)
       
       sizeDotState = plasticState(p)%sizeDotState
              
       plasticState(p)%RKCK45dotState(6,:,cc) = plasticState (p)%dotState(:,cc)
       
       residuum_plastic(1:sizeDotState,g,i,e) = &
         matmul(transpose(plasticState(p)%RKCK45dotState(1:6,1:sizeDotState,cc)),DB) &              ! why transpose? Better to transpose constant DB
       * crystallite_subdt(g,i,e)
       
        plasticState(p)%dotState(:,cc) =  &
         matmul(transpose(plasticState(p)%RKCK45dotState(1:6,1:sizeDotState,cc)), B)                ! why transpose? Better to transpose constant B
         
       do s = 1, phase_Nsources(p)
         sizeDotState = sourceState(p)%p(s)%sizeDotState
       
         sourceState(p)%p(s)%RKCK45dotState(6,:,cc) = sourceState(p)%p(s)%dotState(:,cc)
         
         residuum_source(1:sizeDotState,s,g,i,e) = &
           matmul(transpose(sourceState(p)%p(s)%RKCK45dotState(1:6,1:sizeDotState,cc)),DB) &
         * crystallite_subdt(g,i,e)

         sourceState(p)%p(s)%dotState(:,cc)  = &
           matmul(transpose(sourceState(p)%p(s)%RKCK45dotState(1:6,1:sizeDotState,cc)),B)
       enddo
       
     endif
   enddo; enddo; enddo
 !$OMP END PARALLEL DO

 call update_state(1.0_pReal)
 
 ! --- relative residui and state convergence ---

 !$OMP PARALLEL DO PRIVATE(sizeDotState,p,cc)
 do e = FEsolving_execElem(1),FEsolving_execElem(2)
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
       if (crystallite_todo(g,i,e)) then
         p  = material_phaseAt(g,e); cc = material_phaseMemberAt(g,i,e)
       
         sizeDotState = plasticState(p)%sizeDotState
         
         crystallite_todo(g,i,e) = converged(residuum_plastic(1:sizeDotState,g,i,e), &
                                             plasticState(p)%state(1:sizeDotState,cc), &
                                             plasticState(p)%aTolState(1:sizeDotState))

         do s = 1, phase_Nsources(p)
           sizeDotState = sourceState(p)%p(s)%sizeDotState
         
                  crystallite_todo(g,i,e) = &
                  crystallite_todo(g,i,e) .and. converged(residuum_source(1:sizeDotState,s,g,i,e), &
                                                          sourceState(p)%p(s)%state(1:sizeDotState,cc), &
                                                          sourceState(p)%p(s)%aTolState(1:sizeDotState))
         enddo
     endif
   enddo; enddo; enddo
 !$OMP END PARALLEL DO

 call update_deltaState
 call update_dependentState
 call update_stress(1.0_pReal)
 call setConvergenceFlag
 if (any(plasticState(:)%nonlocal)) call nonlocalConvergenceCheck
 
end subroutine integrateStateRKCK45


!--------------------------------------------------------------------------------------------------
!> @brief sets convergence flag for nonlocal calculations
!> @detail one non-converged nonlocal sets all other nonlocals to non-converged to trigger cut back
!--------------------------------------------------------------------------------------------------
subroutine nonlocalConvergenceCheck
 
 if (any(.not. crystallite_converged .and. .not. crystallite_localPlasticity)) &                    ! any non-local not yet converged (or broken)...
   where( .not. crystallite_localPlasticity) crystallite_converged = .false.

end subroutine nonlocalConvergenceCheck


!--------------------------------------------------------------------------------------------------
!> @brief Sets convergence flag based on "todo": every point that survived the integration (todo is
! still .true. is considered as converged
!> @details: For explicitEuler, RK4 and RKCK45, adaptive Euler and FPI have their on criteria
!--------------------------------------------------------------------------------------------------
subroutine setConvergenceFlag

 integer :: &
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g                                                                                                !< grain index in grain loop
 
 !OMP DO PARALLEL PRIVATE
 do e = FEsolving_execElem(1),FEsolving_execElem(2)
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
       crystallite_converged(g,i,e) = crystallite_todo(g,i,e) .or. crystallite_converged(g,i,e)     ! if still "to do" then converged per definition
 enddo; enddo; enddo
 !OMP END DO PARALLEL

end subroutine setConvergenceFlag


 !--------------------------------------------------------------------------------------------------
 !> @brief determines whether a point is converged
 !--------------------------------------------------------------------------------------------------
 logical pure function converged(residuum,state,aTol)
    
  real(pReal), intent(in), dimension(:) ::&
    residuum, state, aTol
  real(pReal) :: &
    rTol

  rTol = num%rTol_crystalliteState

  converged = all(abs(residuum) <= max(aTol, rTol*abs(state)))

 end function converged


!--------------------------------------------------------------------------------------------------
!> @brief Standard forwarding of state as state = state0 + dotState * (delta t)                        comment seems wrong!
!--------------------------------------------------------------------------------------------------
subroutine update_stress(timeFraction)

 real(pReal), intent(in) :: &
   timeFraction
 integer :: &
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g

 !$OMP PARALLEL DO
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
       crystallite_todo(g,i,e) = integrateStress(g,i,e,timeFraction)
       !$OMP FLUSH(crystallite_todo)
       if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then                  ! if broken non-local...
         !$OMP CRITICAL (checkTodo)
           crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                           ! ...all non-locals skipped
         !$OMP END CRITICAL (checkTodo)
       endif
     endif
   enddo; enddo; enddo
 !$OMP END PARALLEL DO

end subroutine update_stress

!--------------------------------------------------------------------------------------------------
!> @brief tbd
!--------------------------------------------------------------------------------------------------
subroutine update_dependentState
 use constitutive, only: &
   constitutive_dependentState => constitutive_microstructure

 integer ::                                    e, &                                                  ! element index in element loop
                                               i, &                                                  ! integration point index in ip loop
                                               g                                                     ! grain index in grain loop

 !$OMP PARALLEL DO
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
         if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) &
         call constitutive_dependentState(crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), &
                                          g, i, e)
   enddo; enddo; enddo
 !$OMP END PARALLEL DO

end subroutine update_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief Standard forwarding of state as state = state0 + dotState * (delta t)
!--------------------------------------------------------------------------------------------------
subroutine update_state(timeFraction)

 real(pReal), intent(in) :: &
   timeFraction
 integer :: &
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g, &                                                                                             !< grain index in grain loop
   p, &
   c, &
   s, &
   mySize

 !$OMP PARALLEL DO PRIVATE(mySize,p,c)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
         if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
       p = material_phaseAt(g,e); c = material_phaseMemberAt(g,i,e)

       mySize = plasticState(p)%sizeDotState
       plasticState(p)%state(1:mySize,c) = plasticState(p)%subState0(1:mySize,c) &
                                         + plasticState(p)%dotState (1:mySize,c) &
                                         * crystallite_subdt(g,i,e) * timeFraction
       do s = 1, phase_Nsources(p)
         mySize = sourceState(p)%p(s)%sizeDotState
         sourceState(p)%p(s)%state(1:mySize,c) = sourceState(p)%p(s)%subState0(1:mySize,c) &
                                               + sourceState(p)%p(s)%dotState (1:mySize,c) &
                                               * crystallite_subdt(g,i,e) * timeFraction
       enddo
     endif
 enddo; enddo; enddo
 !$OMP END PARALLEL DO

end subroutine update_state


!--------------------------------------------------------------------------------------------------
!> @brief triggers calculation of all new rates
!> if NaN occurs, crystallite_todo is set to FALSE. Any NaN in a nonlocal propagates to all others
!--------------------------------------------------------------------------------------------------
subroutine update_dotState(timeFraction)

 real(pReal), intent(in) :: &
   timeFraction
 integer :: &
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g, &                                                                                             !< grain index in grain loop
   p, &
   c, &
   s 
 logical :: &
   NaN, &
   nonlocalStop
   
   nonlocalStop = .false.

   !$OMP PARALLEL DO PRIVATE (p,c,NaN)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
         !$OMP FLUSH(nonlocalStop)
        if ((crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) .and. .not. nonlocalStop) then
           call constitutive_collectDotState(crystallite_S(1:3,1:3,g,i,e), &
                                             crystallite_Fe, &
                                             crystallite_Fi(1:3,1:3,g,i,e), &
                                             crystallite_Fp, &
                                             crystallite_subdt(g,i,e)*timeFraction, g,i,e)
           p = material_phaseAt(g,e); c = material_phaseMemberAt(g,i,e)
           NaN = any(IEEE_is_NaN(plasticState(p)%dotState(:,c)))
           do s = 1, phase_Nsources(p)
             NaN = NaN .or. any(IEEE_is_NaN(sourceState(p)%p(s)%dotState(:,c)))
           enddo
           if (NaN) then
             crystallite_todo(g,i,e) = .false.                                                      ! this one done (and broken)
             if (.not. crystallite_localPlasticity(g,i,e)) nonlocalStop = .True.
           endif
         endif
   enddo; enddo; enddo
   !$OMP END PARALLEL DO

 if (nonlocalStop) crystallite_todo = crystallite_todo .and. crystallite_localPlasticity 

end subroutine update_DotState


subroutine update_deltaState

 integer :: &
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g, &                                                                                             !< grain index in grain loop
   p, &
   mySize, &
   myOffset, &
   c, &
   s 
 logical :: &
   NaN, &
   nonlocalStop
   
   nonlocalStop = .false.

   !$OMP PARALLEL DO PRIVATE(p,c,myOffset,mySize,NaN)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do g = 1,homogenization_Ngrains(material_homogenizationAt(e))
         !$OMP FLUSH(nonlocalStop)
         if ((crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) .and. .not. nonlocalStop) then
        call constitutive_collectDeltaState(crystallite_S(1:3,1:3,g,i,e), &
                                            crystallite_Fe(1:3,1:3,g,i,e), &
                                            crystallite_Fi(1:3,1:3,g,i,e), &
                                            g,i,e)
         p = material_phaseAt(g,e); c = material_phaseMemberAt(g,i,e)
         myOffset = plasticState(p)%offsetDeltaState
         mySize   = plasticState(p)%sizeDeltaState
         NaN = any(IEEE_is_NaN(plasticState(p)%deltaState(1:mySize,c)))
         
         if (.not. NaN) then
         
           plasticState(p)%state(myOffset + 1: myOffset + mySize,c) = &
           plasticState(p)%state(myOffset + 1: myOffset + mySize,c) + plasticState(p)%deltaState(1:mySize,c)
           do s = 1, phase_Nsources(p)
             myOffset = sourceState(p)%p(s)%offsetDeltaState
             mySize   = sourceState(p)%p(s)%sizeDeltaState
             NaN = NaN .or. any(IEEE_is_NaN(sourceState(p)%p(s)%deltaState(1:mySize,c)))
             
             if (.not. NaN) then
               sourceState(p)%p(s)%state(myOffset + 1:myOffset + mySize,c) = &
               sourceState(p)%p(s)%state(myOffset + 1:myOffset + mySize,c) + sourceState(p)%p(s)%deltaState(1:mySize,c)
             endif
          enddo
        endif
         
         crystallite_todo(g,i,e) = .not. NaN
         if (.not. crystallite_todo(g,i,e)) then                                                             ! if state jump fails, then convergence is broken
           crystallite_converged(g,i,e) = .false.
           if (.not. crystallite_localPlasticity(g,i,e)) nonlocalStop = .true.
         endif
       endif
     enddo; enddo; enddo
   !$OMP END PARALLEL DO
 if (nonlocalStop) crystallite_todo = crystallite_todo .and. crystallite_localPlasticity
 
end subroutine update_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief calculates a jump in the state according to the current state and the current stress
!> returns true, if state jump was successfull or not needed. false indicates NaN in delta state
!--------------------------------------------------------------------------------------------------
logical function stateJump(ipc,ip,el)

 integer, intent(in):: &
   el, &                       ! element index
   ip, &                       ! integration point index
   ipc                         ! grain index

 integer :: &
   c, &
   p, &
   mySource, &
   myOffset, &
   mySize

 c = material_phaseMemberAt(ipc,ip,el)
 p = material_phaseAt(ipc,el)

 call constitutive_collectDeltaState(crystallite_S(1:3,1:3,ipc,ip,el), &
                                     crystallite_Fe(1:3,1:3,ipc,ip,el), &
                                     crystallite_Fi(1:3,1:3,ipc,ip,el), &
                                     ipc,ip,el)

 myOffset = plasticState(p)%offsetDeltaState
 mySize   = plasticState(p)%sizeDeltaState

 if( any(IEEE_is_NaN(plasticState(p)%deltaState(1:mySize,c)))) then                                       ! NaN occured in deltaState
   stateJump = .false.
   return
 endif

 plasticState(p)%state(myOffset + 1:myOffset + mySize,c) = &
 plasticState(p)%state(myOffset + 1:myOffset + mySize,c) + plasticState(p)%deltaState(1:mySize,c)

 do mySource = 1, phase_Nsources(p)
   myOffset = sourceState(p)%p(mySource)%offsetDeltaState
   mySize   = sourceState(p)%p(mySource)%sizeDeltaState
   if (any(IEEE_is_NaN(sourceState(p)%p(mySource)%deltaState(1:mySize,c)))) then   ! NaN occured in deltaState
     stateJump = .false.
     return
   endif
   sourceState(p)%p(mySource)%state(myOffset + 1: myOffset + mySize,c) = &
   sourceState(p)%p(mySource)%state(myOffset + 1: myOffset + mySize,c) + sourceState(p)%p(mySource)%deltaState(1:mySize,c)
 enddo

#ifdef DEBUG
 if (any(dNeq0(plasticState(p)%deltaState(1:mySize,c))) &
     .and. iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0 &
     .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
             .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0)) then
   write(6,'(a,i8,1x,i2,1x,i3, /)') '<< CRYST >> update state at el ip ipc ',el,ip,ipc
   write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> deltaState', plasticState(p)%deltaState(1:mySize,c)
   write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> new state', &
     plasticState(p)%state(myOffset + 1                : &
                           myOffset + mySize,c)
 endif
#endif

 stateJump = .true.

end function stateJump

end module crystallite
