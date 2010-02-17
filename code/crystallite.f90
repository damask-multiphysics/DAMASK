!* $Id$
!***************************************
!*      Module: CRYSTALLITE            *
!***************************************
!* contains:                           *
!* - _init                             *
!* - materialpoint_stressAndItsTangent *
!* - _partitionDeformation             *
!* - _updateState                      *
!* - _stressAndItsTangent       *
!* - _postResults                      *
!***************************************

MODULE crystallite

use prec, only: pReal, pInt
implicit none
!
! ****************************************************************
! *** General variables for the crystallite calculation        ***
! ****************************************************************
integer(pInt), parameter :: crystallite_Nresults = 14_pInt                                   ! phaseID, volume, Euler angles, def gradient

real(pReal), dimension (:,:,:), allocatable ::          crystallite_dt, &                    ! requested time increment of each grain
                                                        crystallite_subdt, &                 ! substepped time increment of each grain
                                                        crystallite_subFrac, &               ! already calculated fraction of increment
                                                        crystallite_subStep, &               ! size of next integration step
                                                        crystallite_Temperature, &           ! Temp of each grain 
                                                        crystallite_partionedTemperature0, & ! Temp of each grain at start of homog inc
                                                        crystallite_subTemperature0          ! Temp of each grain at start of crystallite inc
real(pReal), dimension (:,:,:,:), allocatable ::        crystallite_Tstar_v, &               ! current 2nd Piola-Kirchhoff stress vector (end of converged time step)
                                                        crystallite_Tstar0_v, &              ! 2nd Piola-Kirchhoff stress vector at start of FE inc
                                                        crystallite_partionedTstar0_v, &     ! 2nd Piola-Kirchhoff stress vector at start of homog inc
                                                        crystallite_subTstar0_v, &           ! 2nd Piola-Kirchhoff stress vector at start of crystallite inc
                                                        crystallite_eulerangles              ! euler angles phi1 Phi phi2
real(pReal), dimension (:,:,:,:,:), allocatable ::      crystallite_Fe, &                    ! current "elastic" def grad (end of converged time step)
                                                        crystallite_Fp, &                    ! current plastic def grad (end of converged time step)
                                                        crystallite_invFp, &                 ! inverse of current plastic def grad (end of converged time step)
                                                        crystallite_Fp0, &                   ! plastic def grad at start of FE inc
                                                        crystallite_partionedFp0,&           ! plastic def grad at start of homog inc
                                                        crystallite_subFp0,&                 ! plastic def grad at start of crystallite inc
                                                        crystallite_F0, &                    ! def grad at start of FE inc
                                                        crystallite_partionedF,  &           ! def grad to be reached at end of homog inc
                                                        crystallite_partionedF0, &           ! def grad at start of homog inc
                                                        crystallite_subF,  &                 ! def grad to be reached at end of crystallite inc
                                                        crystallite_subF0, &                 ! def grad at start of crystallite inc
                                                        crystallite_Lp, &                    ! current plastic velocitiy grad (end of converged time step)
                                                        crystallite_Lp0, &                   ! plastic velocitiy grad at start of FE inc
                                                        crystallite_partionedLp0,&           ! plastic velocity grad at start of homog inc
                                                        crystallite_subLp0,&                 ! plastic velocity grad at start of crystallite inc
                                                        crystallite_P, &                     ! 1st Piola-Kirchhoff stress per grain
                                                        crystallite_R, &                     ! crystal orientation (rotation matrix current -> lattice conf)
                                                        crystallite_misorientation           ! misorientation between two neighboring ips (only calculated for single grain IPs)
real(pReal), dimension (:,:,:,:,:,:,:), allocatable ::  crystallite_dPdF, &                  ! individual dPdF per grain
                                                        crystallite_fallbackdPdF             ! dPdF fallback for non-converged grains (elastic prediction)
real(pReal)                                             crystallite_statedamper              ! damping for state update

logical, dimension (:,:,:), allocatable ::              crystallite_localConstitution, &     ! indicates this grain to have purely local constitutive law
                                                        crystallite_requested, &             ! flag to request crystallite calculation
                                                        crystallite_onTrack, &               ! flag to indicate ongoing calculation
                                                        crystallite_converged, &             ! convergence flag
                                                        crystallite_stateConverged, &        ! flag indicating convergence of state
                                                        crystallite_temperatureConverged, &  ! flag indicating convergence of temperature
                                                        crystallite_todo                     ! requested and ontrack but not converged

CONTAINS

!********************************************************************
! allocate and initialize per grain variables
!********************************************************************
subroutine crystallite_init(Temperature)

  !*** variables and functions from other modules ***!
  use prec, only:             pInt, &
                              pReal
  use debug, only:            debug_info, &
                              debug_reset
  use math, only:             math_I3, &
                              math_EulerToR
  use FEsolving, only:        FEsolving_execElem, &
                              FEsolving_execIP
  use mesh, only:             mesh_element, &
                              mesh_NcpElems, &
                              mesh_maxNips, &
                              mesh_maxNipNeighbors
  use material, only:         homogenization_Ngrains, &
                              homogenization_maxNgrains, &
                              material_EulerAngles, &
                              material_phase, &
                              phase_localConstitution
  implicit none

  !*** input variables ***!
  real(pReal) Temperature
 
  !*** output variables ***!
 
  !*** local variables ***!
  integer(pInt)               g, &                          ! grain number
                              i, &                          ! integration point number
                              e, &                          ! element number
                              gMax, &                       ! maximum number of grains
                              iMax, &                       ! maximum number of integration points
                              eMax, &                       ! maximum  number of elements
                              nMax, &                       ! maximum number of ip neighbors
                              myNgrains

  gMax = homogenization_maxNgrains
  iMax = mesh_maxNips
  eMax = mesh_NcpElems
  nMax = mesh_maxNipNeighbors

  allocate(crystallite_Temperature(gMax,iMax,eMax));                     crystallite_Temperature = Temperature
  allocate(crystallite_P(3,3,gMax,iMax,eMax));                                     crystallite_P = 0.0_pReal
  allocate(crystallite_Fe(3,3,gMax,iMax,eMax));                                   crystallite_Fe = 0.0_pReal
  allocate(crystallite_Fp(3,3,gMax,iMax,eMax));                                   crystallite_Fp = 0.0_pReal
  allocate(crystallite_invFp(3,3,gMax,iMax,eMax));                             crystallite_invFp = 0.0_pReal
  allocate(crystallite_Lp(3,3,gMax,iMax,eMax));                                   crystallite_Lp = 0.0_pReal
  allocate(crystallite_Tstar_v(6,gMax,iMax,eMax));                           crystallite_Tstar_v = 0.0_pReal
  allocate(crystallite_F0(3,3,gMax,iMax,eMax));                                   crystallite_F0 = 0.0_pReal
  allocate(crystallite_Fp0(3,3,gMax,iMax,eMax));                                 crystallite_Fp0 = 0.0_pReal
  allocate(crystallite_Lp0(3,3,gMax,iMax,eMax));                                 crystallite_Lp0 = 0.0_pReal
  allocate(crystallite_Tstar0_v(6,gMax,iMax,eMax));                         crystallite_Tstar0_v = 0.0_pReal
  allocate(crystallite_partionedTemperature0(gMax,iMax,eMax)); crystallite_partionedTemperature0 = 0.0_pReal
  allocate(crystallite_partionedF(3,3,gMax,iMax,eMax));                   crystallite_partionedF = 0.0_pReal
  allocate(crystallite_partionedF0(3,3,gMax,iMax,eMax));                 crystallite_partionedF0 = 0.0_pReal
  allocate(crystallite_partionedFp0(3,3,gMax,iMax,eMax));               crystallite_partionedFp0 = 0.0_pReal
  allocate(crystallite_partionedLp0(3,3,gMax,iMax,eMax));               crystallite_partionedLp0 = 0.0_pReal
  allocate(crystallite_partionedTstar0_v(6,gMax,iMax,eMax));       crystallite_partionedTstar0_v = 0.0_pReal
  allocate(crystallite_subTemperature0(gMax,iMax,eMax));             crystallite_subTemperature0 = 0.0_pReal
  allocate(crystallite_subF(3,3,gMax,iMax,eMax));                               crystallite_subF = 0.0_pReal
  allocate(crystallite_subF0(3,3,gMax,iMax,eMax));                             crystallite_subF0 = 0.0_pReal
  allocate(crystallite_subFp0(3,3,gMax,iMax,eMax));                           crystallite_subFp0 = 0.0_pReal
  allocate(crystallite_subLp0(3,3,gMax,iMax,eMax));                           crystallite_subLp0 = 0.0_pReal
  allocate(crystallite_R(3,3,gMax,iMax,eMax));                                     crystallite_R = 0.0_pReal
  allocate(crystallite_eulerangles(3,gMax,iMax,eMax));                   crystallite_eulerangles = 0.0_pReal
  allocate(crystallite_misorientation(4,nMax,gMax,iMax,eMax));        crystallite_misorientation = 0.0_pReal
  allocate(crystallite_subTstar0_v(6,gMax,iMax,eMax));                   crystallite_subTstar0_v = 0.0_pReal
  allocate(crystallite_dPdF(3,3,3,3,gMax,iMax,eMax));                           crystallite_dPdF = 0.0_pReal
  allocate(crystallite_fallbackdPdF(3,3,3,3,gMax,iMax,eMax));           crystallite_fallbackdPdF = 0.0_pReal
  allocate(crystallite_dt(gMax,iMax,eMax));                                       crystallite_dt = 0.0_pReal
  allocate(crystallite_subdt(gMax,iMax,eMax));                                 crystallite_subdt = 0.0_pReal
  allocate(crystallite_subFrac(gMax,iMax,eMax));                             crystallite_subFrac = 0.0_pReal
  allocate(crystallite_subStep(gMax,iMax,eMax));                             crystallite_subStep = 0.0_pReal
  allocate(crystallite_localConstitution(gMax,iMax,eMax));         crystallite_localConstitution = .true.
  allocate(crystallite_requested(gMax,iMax,eMax));                         crystallite_requested = .false.
  allocate(crystallite_onTrack(gMax,iMax,eMax));                             crystallite_onTrack = .true.
  allocate(crystallite_converged(gMax,iMax,eMax));                         crystallite_converged = .true.
  allocate(crystallite_stateConverged(gMax,iMax,eMax));               crystallite_stateConverged = .false.
  allocate(crystallite_temperatureConverged(gMax,iMax,eMax));   crystallite_temperatureConverged = .false.
  allocate(crystallite_todo(gMax,iMax,eMax));                     crystallite_todo = .true.

  !$OMP PARALLEL DO
    do e = FEsolving_execElem(1),FEsolving_execElem(2)           ! iterate over all cp elements
      myNgrains = homogenization_Ngrains(mesh_element(3,e))
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)         ! iterate over IPs of this element
        do g = 1,myNgrains
          crystallite_partionedTemperature0(g,i,e) = Temperature                       ! isothermal assumption
          crystallite_Fp0(:,:,g,i,e) = math_EulerToR(material_EulerAngles(:,g,i,e))    ! plastic def gradient reflects init orientation
          crystallite_Fe(:,:,g,i,e)  = transpose(crystallite_Fp0(:,:,g,i,e))
          crystallite_F0(:,:,g,i,e)  = math_I3
          crystallite_partionedFp0(:,:,g,i,e) = crystallite_Fp0(:,:,g,i,e)
          crystallite_partionedF0(:,:,g,i,e)  = crystallite_F0(:,:,g,i,e)
          crystallite_partionedF(:,:,g,i,e)   = crystallite_F0(:,:,g,i,e)
          crystallite_requested(g,i,e)        = .true.
          crystallite_localConstitution(g,i,e) = phase_localConstitution(material_phase(g,i,e))
        enddo
      enddo
    enddo
  !$OMPEND PARALLEL DO

  call crystallite_orientations()
  call crystallite_stressAndItsTangent(.true.)                 ! request elastic answers
  crystallite_fallbackdPdF = crystallite_dPdF                  ! use initial elastic stiffness as fallback
 
  !    *** Output to MARC output file ***
  !$OMP CRITICAL (write2out)
    write(6,*)
    write(6,*) '<<<+-  crystallite init  -+>>>'
    write(6,*) '$Id$'
    write(6,*)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Nresults:              ', crystallite_Nresults
    write(6,*)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Temperature:           ', shape(crystallite_Temperature)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Fe:                    ', shape(crystallite_Fe)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Fp:                    ', shape(crystallite_Fp)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Lp:                    ', shape(crystallite_Lp)
    write(6,'(a35,x,7(i5,x))') 'crystallite_F0:                    ', shape(crystallite_F0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Fp0:                   ', shape(crystallite_Fp0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Lp0:                   ', shape(crystallite_Lp0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_partionedF:            ', shape(crystallite_partionedF)
    write(6,'(a35,x,7(i5,x))') 'crystallite_partionedTemp0:        ', shape(crystallite_partionedTemperature0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_partionedF0:           ', shape(crystallite_partionedF0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_partionedFp0:          ', shape(crystallite_partionedFp0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_partionedLp0:          ', shape(crystallite_partionedLp0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subF:                  ', shape(crystallite_subF)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subTemperature0:       ', shape(crystallite_subTemperature0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subF0:                 ', shape(crystallite_subF0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subFp0:                ', shape(crystallite_subFp0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subLp0:                ', shape(crystallite_subLp0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_P:                     ', shape(crystallite_P)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Tstar_v:               ', shape(crystallite_Tstar_v)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Tstar0_v:              ', shape(crystallite_Tstar0_v)
    write(6,'(a35,x,7(i5,x))') 'crystallite_partionedTstar0_v:     ', shape(crystallite_partionedTstar0_v)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subTstar0_v:           ', shape(crystallite_subTstar0_v)
    write(6,'(a35,x,7(i5,x))') 'crystallite_dPdF:                  ', shape(crystallite_dPdF)
    write(6,'(a35,x,7(i5,x))') 'crystallite_fallbackdPdF:          ', shape(crystallite_fallbackdPdF)
    write(6,'(a35,x,7(i5,x))') 'crystallite_R:                     ', shape(crystallite_R)
    write(6,'(a35,x,7(i5,x))') 'crystallite_eulerangles:           ', shape(crystallite_eulerangles)
    write(6,'(a35,x,7(i5,x))') 'crystallite_misorientation:        ', shape(crystallite_misorientation)
    write(6,'(a35,x,7(i5,x))') 'crystallite_dt:                    ', shape(crystallite_dt)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subdt:                 ', shape(crystallite_subdt)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subFrac:               ', shape(crystallite_subFrac)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subStep:               ', shape(crystallite_subStep)
    write(6,'(a35,x,7(i5,x))') 'crystallite_localConstitution:     ', shape(crystallite_localConstitution)
    write(6,'(a35,x,7(i5,x))') 'crystallite_requested:             ', shape(crystallite_requested)
    write(6,'(a35,x,7(i5,x))') 'crystallite_onTrack:               ', shape(crystallite_onTrack)
    write(6,'(a35,x,7(i5,x))') 'crystallite_converged:             ', shape(crystallite_converged)
    write(6,'(a35,x,7(i5,x))') 'crystallite_stateConverged:        ', shape(crystallite_stateConverged)
    write(6,'(a35,x,7(i5,x))') 'crystallite_temperatureConverged:  ', shape(crystallite_temperatureConverged)
    write(6,'(a35,x,7(i5,x))') 'crystallite_todo:                  ', shape(crystallite_todo)
    write(6,*)
    write(6,*) 'Number of nonlocal grains: ',count(.not. crystallite_localConstitution)
    call flush(6)
  !$OMPEND CRITICAL (write2out)

  call debug_info()
  call debug_reset()

  return

endsubroutine


 
!********************************************************************
! calculate stress (P) and tangent (dPdF) for crystallites
!********************************************************************
subroutine crystallite_stressAndItsTangent(updateJaco)

  !*** variables and functions from other modules ***!
  use prec, only:                                       pInt, &
                                                        pReal
  use numerics, only:                                   subStepMinCryst, &
                                                        subStepSizeCryst, &
                                                        stepIncreaseCryst, &
                                                        pert_Fg, &
                                                        pert_method, &
                                                        nState, &
                                                        nCryst
  use debug, only:                                      debugger, &
                                                        selectiveDebugger, &
                                                        debug_e, &
                                                        debug_i, &
                                                        debug_g, &
                                                        debug_CrystalliteLoopDistribution, &
                                                        debug_CrystalliteStateLoopDistribution, &
                                                        debug_StiffnessStateLoopDistribution
  use IO, only:                                         IO_warning
  use math, only:                                       math_inv3x3, &
                                                        math_mul33x33, &
                                                        math_mul66x6, &
                                                        math_Mandel6to33, &
                                                        math_Mandel33to6, &
                                                        math_I3, &
                                                        math_Plain3333to99
  use FEsolving, only:                                  FEsolving_execElem, & 
                                                        FEsolving_execIP, &
                                                        theInc, &
                                                        cycleCounter
  use mesh, only:                                       mesh_element, &
                                                        mesh_NcpElems, &
                                                        mesh_maxNips
  use material, only:                                   homogenization_Ngrains, &
                                                        homogenization_maxNgrains
  use constitutive, only:                               constitutive_maxSizeState, &
                                                        constitutive_maxSizeDotState, &
                                                        constitutive_sizeState, &
                                                        constitutive_sizeDotState, &
                                                        constitutive_state, &
                                                        constitutive_subState0, &
                                                        constitutive_partionedState0, &
                                                        constitutive_homogenizedC, &
                                                        constitutive_dotState, &
                                                        constitutive_previousDotState, &
                                                        constitutive_previousDotState2, &
                                                        constitutive_collectDotState, &
                                                        constitutive_dotTemperature, &
                                                        constitutive_microstructure

  implicit none

  !*** input variables ***!
  logical, intent(in) ::                                updateJaco                    ! flag indicating wehther we want to update the Jacobian (stiffness) or not

  !*** output variables ***!

  !*** local variables ***!
  real(pReal)                                           myTemperature, &              ! local copy of the temperature
                                                        myPert                        ! perturbation with correct sign
  real(pReal), dimension(3,3) ::                        invFp, &                      ! inverse of the plastic deformation gradient
                                                        Fe_guess, &                   ! guess for elastic deformation gradient
                                                        Tstar                         ! 2nd Piola-Kirchhoff stress tensor
  integer(pInt)                                         NiterationCrystallite, &      ! number of iterations in crystallite loop
                                                        NiterationState               ! number of iterations in state loop
  integer(pInt)                                         e, ee, &                      ! element index
                                                        i, ii, &                      ! integration point index
                                                        g, gg, &                      ! grain index
                                                        k, &
                                                        l, &
                                                        perturbation , &              ! loop counter for forward,backward perturbation mode
                                                        comp, &
                                                        myNgrains, &
                                                        mySizeState, &
                                                        mySizeDotState
  integer(pInt), dimension(2,9) ::                      kl
  logical                                               onTrack, &                    ! flag indicating whether we are still on track
                                                        temperatureConverged, &       ! flag indicating if temperature converged
                                                        stateConverged, &             ! flag indicating if state converged
                                                        converged                     ! flag indicating if iteration converged
  real(pReal), dimension(9,9) ::                        dPdF99
  real(pReal), dimension(3,3,3,3,2) ::                  dPdF_perturbation
  real(pReal), dimension(constitutive_maxSizeDotState) :: delta_dotState1, &          ! difference between current and previous dotstate
                                                        delta_dotState2               ! difference between previousDotState and previousDotState2
  real(pReal)                                           dot_prod12, &
                                                        dot_prod22  
  real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                                        storedF, &
                                                        storedFp, &
                                                        storedInvFp, &
                                                        storedFe, &
                                                        storedLp, &
                                                        storedP
  real(pReal), dimension(6,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                                        storedTstar_v
  real(pReal), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                                        storedTemperature
  real(pReal), dimension(constitutive_maxSizeState,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                                        storedState
  real(pReal), dimension(constitutive_maxSizeDotState,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                                        storedDotState
  logical, dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                                        storedConvergenceFlag
  logical, dimension(3,3) ::                            mask
  logical																								forceLocalStiffnessCalculation = .true. ! flag indicating that stiffness calculation is always done locally 
                                                        
  ! ------ initialize to starting condition ------
  
!$OMP CRITICAL (write2out)
!  write (6,*)
!  write (6,*) 'Crystallite request from Materialpoint'
!  write (6,'(a,/,(f12.7,x))')      'crystallite_partionedTemp0 of 1 1 1' ,crystallite_partionedTemperature0(1,1,1)
!  write (6,'(a,/,3(3(f12.7,x)/))') 'crystallite_partionedF0  of 1 1 1'   ,crystallite_partionedF0(1:3,:,1,1,1)
!  write (6,'(a,/,3(3(f12.7,x)/))') 'crystallite_partionedFp0 of 1 1 1'   ,crystallite_partionedFp0(1:3,:,1,1,1)
!  write (6,'(a,/,3(3(f12.7,x)/))') 'crystallite_partionedF   of 1 1 1'   ,crystallite_partionedF(1:3,:,1,1,1)
!  write (6,'(a,/,3(3(f12.7,x)/))') 'crystallite_partionedLp0 of 1 1 1'   ,crystallite_partionedLp0(1:3,:,1,1,1)
!$OMPEND CRITICAL (write2out)

  crystallite_subStep = 0.0_pReal

  !$OMP PARALLEL DO
    do e = FEsolving_execElem(1),FEsolving_execElem(2)                 ! iterate over elements to be processed
      myNgrains = homogenization_Ngrains(mesh_element(3,e))
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)               ! iterate over IPs of this element to be processed
        do g = 1,myNgrains
          if (crystallite_requested(g,i,e)) then                       ! initialize restoration point of ...
            crystallite_subTemperature0(g,i,e) = crystallite_partionedTemperature0(g,i,e) ! ...temperature
            constitutive_subState0(g,i,e)%p   = constitutive_partionedState0(g,i,e)%p     ! ...microstructure
            crystallite_subFp0(:,:,g,i,e)     = crystallite_partionedFp0(:,:,g,i,e)       ! ...plastic def grad
            crystallite_subLp0(:,:,g,i,e)     = crystallite_partionedLp0(:,:,g,i,e)       ! ...plastic velocity grad
            crystallite_subF0(:,:,g,i,e)      = crystallite_partionedF0(:,:,g,i,e)        ! ...def grad
            crystallite_subTstar0_v(:,g,i,e)  = crystallite_partionedTstar0_v(:,g,i,e)    ! ...2nd PK stress

            crystallite_subFrac(g,i,e) = 0.0_pReal
            crystallite_subStep(g,i,e) = 1.0_pReal/subStepSizeCryst    ! <<added flexibility in cutback size>>
            crystallite_onTrack(g,i,e) = .true.
            crystallite_converged(g,i,e) = .false.                     ! pretend failed step of twice the required size
          endif
        enddo
      enddo
    enddo
  !$OMPEND PARALLEL DO
 

  ! --+>> crystallite loop <<+--

  NiterationCrystallite = 0_pInt

  do while (any(crystallite_subStep(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) > subStepMinCryst))      ! cutback loop for crystallites

    !$OMP PARALLEL DO
      do e = FEsolving_execElem(1),FEsolving_execElem(2)               ! iterate over elements to be processed
        myNgrains = homogenization_Ngrains(mesh_element(3,e))
        do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)             ! iterate over IPs of this element to be processed
          do g = 1,myNgrains
            selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g) 
            if (crystallite_converged(g,i,e)) then
              if (selectiveDebugger) then
                !$OMP CRITICAL (write2out)
                  write(6,'(a21,f10.8,a32,f10.8,a35)') 'winding forward from ', &
                    crystallite_subFrac(g,i,e),' to current crystallite_subfrac ', &
                    crystallite_subFrac(g,i,e)+crystallite_subStep(g,i,e),' in crystallite_stressAndItsTangent'
                  write(6,*)
                !$OMPEND CRITICAL (write2out)
              endif
              crystallite_subFrac(g,i,e) = crystallite_subFrac(g,i,e) + crystallite_subStep(g,i,e)
              crystallite_subStep(g,i,e) = min(1.0_pReal-crystallite_subFrac(g,i,e), &
                                               stepIncreaseCryst*crystallite_subStep(g,i,e))   ! <<introduce possibility for acceleration>>
              if (crystallite_subStep(g,i,e) > subStepMinCryst) then
                crystallite_subTemperature0(g,i,e) = crystallite_Temperature(g,i,e)     ! wind forward...
                crystallite_subF0(:,:,g,i,e)      = crystallite_subF(:,:,g,i,e)         ! ...def grad
                crystallite_subFp0(:,:,g,i,e)     = crystallite_Fp(:,:,g,i,e)           ! ...plastic def grad
                crystallite_subLp0(:,:,g,i,e)     = crystallite_Lp(:,:,g,i,e)           ! ...plastic velocity gradient
                constitutive_subState0(g,i,e)%p   = constitutive_state(g,i,e)%p         ! ...microstructure
                crystallite_subTstar0_v(:,g,i,e)  = crystallite_Tstar_v(:,g,i,e)        ! ...2nd PK stress
              elseif (crystallite_onTrack(g,i,e)) then                                       ! this crystallite just converged
               !$OMP CRITICAL (distributionCrystallite)
                 debug_CrystalliteLoopDistribution(min(nCryst+1,NiterationCrystallite)) = &
                   debug_CrystalliteLoopDistribution(min(nCryst+1,NiterationCrystallite)) + 1
               !$OMPEND CRITICAL (distributionCrystallite)
              endif
            else
              crystallite_subStep(g,i,e)    = subStepSizeCryst*crystallite_subStep(g,i,e)   ! cut step in half and restore...
              crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e)       ! ...temperature
              crystallite_Fp(:,:,g,i,e)     = crystallite_subFp0(:,:,g,i,e)             ! ...plastic def grad
              crystallite_invFp(:,:,g,i,e)  = math_inv3x3(crystallite_Fp(:,:,g,i,e))
              crystallite_Lp(:,:,g,i,e)     = crystallite_subLp0(:,:,g,i,e)             ! ...plastic velocity grad
              constitutive_state(g,i,e)%p   = constitutive_subState0(g,i,e)%p           ! ...microstructure
              crystallite_Tstar_v(:,g,i,e)  = crystallite_subTstar0_v(:,g,i,e)          ! ...2nd PK stress
              if (selectiveDebugger) then
                !$OMP CRITICAL (write2out)
                  write(6,'(a78,f10.8)') 'cutback step in crystallite_stressAndItsTangent with new crystallite_subStep: ',&
                                         crystallite_subStep(g,i,e)
                  write(6,*)
                !$OMPEND CRITICAL (write2out)
              endif
            endif

            crystallite_onTrack(g,i,e) = crystallite_subStep(g,i,e) > subStepMinCryst   ! still on track or already done (beyond repair)
            if (crystallite_onTrack(g,i,e)) then                                        ! specify task (according to substep)
              crystallite_subF(:,:,g,i,e)  = crystallite_subF0(:,:,g,i,e) + &
                                             crystallite_subStep(g,i,e) * &
                                             (crystallite_partionedF(:,:,g,i,e) - crystallite_partionedF0(:,:,g,i,e))
              crystallite_Fe(:,:,g,i,e)    = math_mul33x33(crystallite_subF(:,:,g,i,e),crystallite_invFp(:,:,g,i,e))
              crystallite_subdt(g,i,e)     = crystallite_subStep(g,i,e) * crystallite_dt(g,i,e)
              crystallite_converged(g,i,e) = .false.                                  ! start out non-converged
            endif
          enddo
        enddo
      enddo
    !$OMPEND PARALLEL DO

    crystallite_todo = (            crystallite_requested & 
                        .and.       crystallite_onTrack &
                        .and. .not. crystallite_converged)

    ! --+>> preguess for state <<+--
    !
    ! incrementing by crystallite_subdt
    ! based on constitutive_subState0
    ! results in constitutive_state
    ! first loop for collection of state evolution based on old state
    ! second loop for updating to new state
   
    !$OMP PARALLEL DO
      do e = FEsolving_execElem(1),FEsolving_execElem(2)                              ! iterate over elements to be processed
        myNgrains = homogenization_Ngrains(mesh_element(3,e))
        do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                            ! iterate over IPs of this element to be processed
          do g = 1,myNgrains
            if (crystallite_todo(g,i,e)) then                                         ! all undone crystallites
              call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Tstar_v(:,g,i,e), crystallite_Fe, &
                                               crystallite_Fp, g, i, e)               ! update dependent state variables to be consistent with basic states
              constitutive_previousDotState2(g,i,e)%p = 0.0_pReal
              constitutive_previousDotState(g,i,e)%p = 0.0_pReal
              constitutive_dotState(g,i,e)%p = 0.0_pReal                              ! zero out dotStates
            endif
     enddo; enddo; enddo
    !$OMPEND PARALLEL DO
    !$OMP PARALLEL DO
      do e = FEsolving_execElem(1),FEsolving_execElem(2)                              ! iterate over elements to be processed
        myNgrains = homogenization_Ngrains(mesh_element(3,e))
        do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                            ! iterate over IPs of this element to be processed
          do g = 1,myNgrains
!            selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g)
            if (crystallite_todo(g,i,e)) then                                         ! all undone crystallites
              call constitutive_collectDotState(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), &
                                                crystallite_Fe, crystallite_Fp, crystallite_Temperature(g,i,e), & 
                                                crystallite_misorientation(:,:,g,i,e), crystallite_subdt(g,i,e), g, i, e)
            endif
      enddo; enddo; enddo
    !$OMPEND PARALLEL DO

    crystallite_statedamper = 1.0_pReal

    !$OMP PARALLEL DO
      do e = FEsolving_execElem(1),FEsolving_execElem(2)                                ! iterate over elements to be processed
        myNgrains = homogenization_Ngrains(mesh_element(3,e))
        do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                              ! iterate over IPs of this element to be processed
          do g = 1,myNgrains
!            selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g)
            if (crystallite_todo(g,i,e)) then                                           ! all undone crystallites
              crystallite_stateConverged(g,i,e) = crystallite_updateState(g,i,e)        ! update state
              crystallite_temperatureConverged(g,i,e) = crystallite_updateTemperature(g,i,e) ! update temperature
              crystallite_converged(g,i,e) = .false.                                    ! force at least one iteration step even if state already converged
            endif
      enddo; enddo; enddo
    !$OMPEND PARALLEL DO
   
    ! --+>> state loop <<+--
   
    NiterationState = 0_pInt
    
    do while ( any(crystallite_todo(:,:,FEsolving_execELem(1):FEsolving_execElem(2))) &
             .and. NiterationState < nState)                                          ! convergence loop for crystallite
     
      NiterationState = NiterationState + 1_pInt
      ! --+>> stress integration <<+--
      !
      ! incrementing by crystallite_subdt
      ! based on crystallite_subF0,.._subFp0,.._subLp0
      !          constitutive_state is internally interpolated with .._subState0
      !          to account for substepping within _integrateStress
      ! results in crystallite_Fp,.._Lp
      !$OMP PARALLEL DO
      do e = FEsolving_execElem(1),FEsolving_execElem(2)                            ! iterate over elements to be processed
        myNgrains = homogenization_Ngrains(mesh_element(3,e))
        do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                          ! iterate over IPs of this element to be processed
          do g = 1,myNgrains
!            selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g)
            if (crystallite_todo(g,i,e)) &                                          ! all undone crystallites
              crystallite_onTrack(g,i,e) = crystallite_integrateStress(g,i,e)
      enddo; enddo; enddo
      !$OMPEND PARALLEL DO

      if (debugger) then
        !$OMP CRITICAL (write2out)
          write(6,*) count(crystallite_onTrack(:,:,:)),'grains onTrack after stress integration'
        !$OMPEND CRITICAL (write2out)
      endif

      crystallite_todo = crystallite_todo .and. crystallite_onTrack                   ! continue with non-broken grains

      if (any(.not. crystallite_onTrack .and. .not. crystallite_localConstitution)) & ! any non-local is broken?
          crystallite_todo = crystallite_todo .and. crystallite_localConstitution     ! all nonlocal crystallites can be skipped

      if (debugger) then
        !$OMP CRITICAL (write2out)
          write(6,*) count(crystallite_todo(:,:,:)),'grains todo after stress integration'
        !$OMPEND CRITICAL (write2out)
      endif

      ! --+>> state integration <<+--
      !
      ! incrementing by crystallite_subdt
      ! based on constitutive_subState0
      ! results in constitutive_state
      ! first loop for collection of state evolution based on old state
      ! second loop for updating to new state

      !$OMP PARALLEL DO
      do e = FEsolving_execElem(1),FEsolving_execElem(2)                              ! iterate over elements to be processed
        myNgrains = homogenization_Ngrains(mesh_element(3,e))
        do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                            ! iterate over IPs of this element to be processed
          do g = 1,myNgrains
            if (crystallite_todo(g,i,e)) then                                         ! all undone crystallites
              constitutive_previousDotState2(g,i,e)%p = constitutive_previousDotState(g,i,e)%p
              constitutive_previousDotState(g,i,e)%p = constitutive_dotState(g,i,e)%p
              constitutive_dotState(g,i,e)%p = 0.0_pReal                              ! zero out dotState
            endif
      enddo; enddo; enddo
      !$OMPEND PARALLEL DO

      crystallite_statedamper = 1.0_pReal

      !$OMP PARALLEL DO
      do e = FEsolving_execElem(1),FEsolving_execElem(2)                              ! iterate over elements to be processed
        myNgrains = homogenization_Ngrains(mesh_element(3,e))
        do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                            ! iterate over IPs of this element to be processed
          do g = 1,myNgrains
!            selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g)
            if (crystallite_todo(g,i,e)) then                                         ! all undone crystallites
              call constitutive_collectDotState(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), &
                                                crystallite_Fe, crystallite_Fp, crystallite_Temperature(g,i,e), & 
                                                crystallite_misorientation(:,:,g,i,e), crystallite_subdt(g,i,e), g, i, e)                
              delta_dotState1 = constitutive_dotState(g,i,e)%p - constitutive_previousDotState(g,i,e)%p
              delta_dotState2 = constitutive_previousDotState(g,i,e)%p - constitutive_previousDotState2(g,i,e)%p
              dot_prod12 = dot_product(delta_dotState1, delta_dotState2)
              dot_prod22 = dot_product(delta_dotState2, delta_dotState2)
              if (      dot_prod22 > 0.0_pReal &
                  .and. (     dot_prod12 < 0.0_pReal &
                         .or. dot_product(constitutive_dotState(g,i,e)%p, constitutive_previousDotState(g,i,e)%p) < 0.0_pReal) ) &
                crystallite_statedamper = min(crystallite_statedamper, &
                                              0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22) )
              endif
      enddo; enddo; enddo
      !$OMPEND PARALLEL DO

      !$OMP PARALLEL DO
        do e = FEsolving_execElem(1),FEsolving_execElem(2)                            ! iterate over elements to be processed
          myNgrains = homogenization_Ngrains(mesh_element(3,e))
          do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                          ! iterate over IPs of this element to be processed
            do g = 1,myNgrains
!              selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g)
              if (crystallite_todo(g,i,e)) then                                       ! all undone crystallites
                crystallite_stateConverged(g,i,e) = crystallite_updateState(g,i,e)    ! update state
                crystallite_temperatureConverged(g,i,e) = crystallite_updateTemperature(g,i,e)  ! update temperature
                crystallite_converged(g,i,e) = crystallite_stateConverged(g,i,e) .and. crystallite_temperatureConverged(g,i,e)
                if (crystallite_converged(g,i,e)) then
                  !$OMP CRITICAL (distributionState)
                    debug_CrystalliteStateLoopDistribution(NiterationState) = &
                      debug_CrystalliteStateLoopDistribution(NiterationState) + 1
                  !$OMPEND CRITICAL (distributionState)
                endif
              endif
            enddo
          enddo
        enddo
      !$OMPEND PARALLEL DO

      if (debugger) then
        !$OMP CRITICAL (write2out)
          write(6,*) count(crystallite_converged(:,:,:)),'grains converged after state integration no.', NiterationState
          write(6,*)
!          write(6,'(8(L,x))') crystallite_converged(:,:,:)
!          do e = FEsolving_execElem(1),FEsolving_execElem(2)
!            if (any(.not. crystallite_converged(:,:,e))) &
!              write(6,'(i4,8(x,L))') e, crystallite_converged(:,:,e)
!          enddo
        !$OMPEND CRITICAL (write2out)
      endif
      if (any(.not. crystallite_converged .and. .not. crystallite_localConstitution)) & ! any non-local not yet converged?
          crystallite_converged = crystallite_converged .and. crystallite_localConstitution       ! all non-local not converged
         
      crystallite_todo = crystallite_todo .and. .not. crystallite_converged             ! skip all converged
      
      if (debugger) then
        !$OMP CRITICAL (write2out)
          write(6,*) count(crystallite_converged(:,:,:)),'grains converged after non-local check'
          write(6,*) count(crystallite_todo(:,:,:)),'grains todo after state integration no.', NiterationState
          write(6,*)
        !$OMPEND CRITICAL (write2out)
      endif
      
    enddo                                                                             ! crystallite convergence loop  
    NiterationCrystallite = NiterationCrystallite + 1
    
  enddo                                                                               ! cutback loop
  
  ! ------ check for non-converged crystallites ------
  !$OMP PARALLEL DO
    do e = FEsolving_execElem(1),FEsolving_execElem(2)                                ! iterate over elements to be processed
      myNgrains = homogenization_Ngrains(mesh_element(3,e))
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                              ! iterate over IPs of this element to be processed
        do g = 1,myNgrains
          if (.not. crystallite_converged(g,i,e)) then                                ! respond fully elastically (might be not required due to becoming terminally ill anyway)
!            call IO_warning(600,e,i,g)
            invFp = math_inv3x3(crystallite_partionedFp0(:,:,g,i,e))
            Fe_guess = math_mul33x33(crystallite_partionedF(:,:,g,i,e),invFp)
            Tstar = math_Mandel6to33( &
                                      math_mul66x6( 0.5_pReal*constitutive_homogenizedC(g,i,e), &
                                                    math_Mandel33to6( math_mul33x33(transpose(Fe_guess),Fe_guess) - math_I3 ) &
                                                  ) &
                                    )
            crystallite_P(:,:,g,i,e) = math_mul33x33(Fe_guess,math_mul33x33(Tstar,transpose(invFp)))
          endif
        enddo
      enddo
    enddo
  !$OMPEND PARALLEL DO


  ! --+>> stiffness calculation <<+--
  
  if(updateJaco) then                                                                                   ! Jacobian required
    
    crystallite_statedamper = 1.0_pReal
    
    do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                  ! iterate over elements to be processed
      myNgrains = homogenization_Ngrains(mesh_element(3,e))
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                ! iterate over IPs of this element to be processed
        do g = 1,myNgrains
          mySizeState = constitutive_sizeState(g,i,e)                                                   ! number of state variables for this grain
          mySizeDotState = constitutive_sizeDotState(g,i,e)                                             ! number of dotStates for this grain
          storedState(1:mySizeState,g,i,e) = constitutive_state(g,i,e)%p                                ! remember unperturbed, converged state, ...
          storedDotState(1:mySizeDotState,g,i,e) = constitutive_dotState(g,i,e)%p                       ! ... dotStates, ...
    enddo; enddo; enddo
		storedTemperature = crystallite_Temperature                                         								! ... Temperature, ...
		storedF = crystallite_subF                                            															! ... and kinematics
		storedFp = crystallite_Fp
		storedInvFp = crystallite_invFp
		storedFe = crystallite_Fe
		storedLp = crystallite_Lp
		storedTstar_v = crystallite_Tstar_v
		storedP = crystallite_P
		storedConvergenceFlag = crystallite_converged

    if (all(crystallite_localConstitution) .or. theInc < 2 .or. forceLocalStiffnessCalculation) then      ! all grains have local constitution, so local convergence of perturbed grain is sufficient
    
      !$OMP PARALLEL DO
        do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                ! iterate over elements to be processed
          myNgrains = homogenization_Ngrains(mesh_element(3,e))
          do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                              ! iterate over IPs of this element to be processed
            do g = 1,myNgrains
!               selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g)
              if (crystallite_requested(g,i,e)) then                                                      ! first check whether is requested at all!
                if (crystallite_converged(g,i,e)) then                                                    ! grain converged in above iteration
                  if (selectiveDebugger) then
                    !$OMP CRITICAL (write2out)
                    write (6,*) '#############'
                    write (6,*) 'central solution of cryst_StressAndTangent'
                    write (6,*) '#############'
                    write (6,'(a8,3(x,i4),/,3(3(f12.4,x)/))') '    P of', g, i, e, storedP(1:3,:,g,i,e)/1e6
                    write (6,'(a8,3(x,i4),/,3(3(f12.8,x)/))') '   Fp of', g, i, e, storedFp(1:3,:,g,i,e)
                    write (6,'(a8,3(x,i4),/,3(3(f12.8,x)/))') '   Lp of', g, i, e, storedLp(1:3,:,g,i,e)
                    !$OMPEND CRITICAL (write2out)
                  endif
                  
                  do perturbation = 1,2
                    if (iand(pert_method,perturbation) > 0) then
                      myPert = -pert_Fg * (-1.0_pReal)**perturbation                                        ! forward or backward perturbation
                      do k = 1,3                                                                            ! perturbation...
                      do l = 1,3                                                                            ! ...components to the positive direction
                      crystallite_subF(k,l,g,i,e) = crystallite_subF(k,l,g,i,e) + myPert                    ! perturb single component (either forward or backward)
                      if (selectiveDebugger) then
                        !$OMP CRITICAL (write2out)
                        write (6,'(i1,x,i1)') k,l
                        write (6,'(a8,3(x,i4),/,3(3(f12.6,x)/))') 'pertF of', g, i, e, crystallite_subF(1:3,:,g,i,e)
                        !$OMPEND CRITICAL (write2out)
                      endif
                      onTrack = .true.
                      converged = .false.
                      NiterationState = 0_pInt
                      do while(.not. converged .and. onTrack .and. NiterationState < nState)              ! keep cycling until done (potentially non-converged)
                        NiterationState = NiterationState + 1_pInt
                        onTrack = crystallite_integrateStress(g,i,e)                                      ! stress of perturbed situation (overwrites _P,_Tstar_v,_Fp,_Lp,_Fe)
                        if (onTrack) then 
                          constitutive_dotState(g,i,e)%p = 0.0_pReal
                          call constitutive_collectDotState(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), &
                                                            crystallite_Fe, crystallite_Fp, crystallite_Temperature(g,i,e), &
                                                            crystallite_misorientation(:,:,g,i,e), crystallite_subdt(g,i,e), &
                                                            g,i,e)
                          stateConverged = crystallite_updateState(g,i,e)                                 ! update state
                          temperatureConverged = crystallite_updateTemperature(g,i,e)                     ! update temperature
                          converged = stateConverged .and. temperatureConverged
                        endif
                        if (selectiveDebugger) then
                          !$OMP CRITICAL (write2out)
                          write (6,*) '-------------'
                          write (6,'(a,x,l,x,l)') 'ontrack + converged:',onTrack,converged
                          write (6,'(a12,3(x,i4),/,3(3(f12.4,x)/))') 'pertP/MPa of', g, i, e, crystallite_P(1:3,:,g,i,e)/1e6
                          write (6,'(a12,3(x,i4),/,3(3(f12.4,x)/))') 'DP/MPa    of', g, i, e, &
                                                                               (crystallite_P(1:3,:,g,i,e)-storedP(1:3,:,g,i,e))/1e6
                          !$OMPEND CRITICAL (write2out)
                        endif
                      enddo
                      if (converged) &                                                                    ! converged state warrants stiffness update
                        dPdF_perturbation(:,:,k,l,perturbation) = (crystallite_P(:,:,g,i,e) - storedP(:,:,g,i,e))/myPert       ! tangent dP_ij/dFg_kl

                      mySizeState = constitutive_sizeState(g,i,e)                                                   ! number of state variables for this grain
                      mySizeDotState = constitutive_sizeDotState(g,i,e)                                             ! number of dotStates for this grain
                      constitutive_state(g,i,e)%p = storedState(1:mySizeState,g,i,e)
                      constitutive_dotState(g,i,e)%p = storedDotState(1:mySizeDotState,g,i,e)
                      crystallite_Temperature(g,i,e) = storedTemperature(g,i,e)
                      crystallite_subF(:,:,g,i,e) = storedF(:,:,g,i,e)
                      crystallite_Fp(:,:,g,i,e) = storedFp(:,:,g,i,e) 
                      crystallite_invFp(:,:,g,i,e) = storedInvFp(:,:,g,i,e)
                      crystallite_Fe(:,:,g,i,e) = storedFe(:,:,g,i,e)
                      crystallite_Lp(:,:,g,i,e) = storedLp(:,:,g,i,e)
                      crystallite_Tstar_v(:,g,i,e) = storedTstar_v(:,g,i,e)
                      crystallite_P(:,:,g,i,e) = storedP(:,:,g,i,e)
                      !$OMP CRITICAL (out)
                        debug_StiffnessStateLoopDistribution(NiterationState) = &
                          debug_StiffnessstateLoopDistribution(NiterationState) + 1
                      !$OMPEND CRITICAL (out)
                      enddo; enddo
                    endif
                  enddo                                                                                    ! perturbation direction
                  select case(pert_method)
                    case (1)
                       crystallite_dPdF(:,:,:,:,g,i,e) = dPdF_perturbation(:,:,:,:,1)
                    case (2)
                       crystallite_dPdF(:,:,:,:,g,i,e) = dPdF_perturbation(:,:,:,:,2)
                    case (3)
                       crystallite_dPdF(:,:,:,:,g,i,e) = 0.5_pReal*(dPdF_perturbation(:,:,:,:,1)+dPdF_perturbation(:,:,:,:,2))
                  end select
                else                                                                                    ! grain did not converge
                  crystallite_dPdF(:,:,:,:,g,i,e) = crystallite_fallbackdPdF(:,:,:,:,g,i,e)             ! use (elastic) fallback
                endif               ! grain convergence
              endif                 ! grain request
            enddo                   ! grain   loop
          enddo                     ! ip      loop
        enddo                       ! element loop
      !$OMPEND PARALLEL DO
      
    elseif (any(.not. crystallite_localConstitution)) then                                              ! if any nonlocal grain present, we have to do a full loop over all grains after each perturbance
      
      do e = FEsolving_execElem(1),FEsolving_execElem(2)
        myNgrains = homogenization_Ngrains(mesh_element(3,e))
        do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
          do g = 1,myNgrains
            
            ! perturb components in the order of biggest change in F (-> component with biggest change in F is perturbed in first cycle, component with second biggest change in next cycle, ...)
            mask = .true.
            do comp = 1,9
              kl(:,comp) = maxloc(abs(crystallite_subF(:,:,g,i,e)-crystallite_F0(:,:,g,i,e)), mask)
              mask(kl(1,comp),kl(2,comp)) = .false.
            enddo
            k = kl(1,mod((cycleCounter-1)/2+1,9))
            l = kl(2,mod((cycleCounter-1)/2+1,9))
            
            crystallite_subF(k,l,g,i,e) = crystallite_subF(k,l,g,i,e) + pert_Fg                         ! perturb single component
            
            NiterationState = 0_pInt
            crystallite_todo = .true.
            do while ( any(crystallite_todo(:,:,FEsolving_execELem(1):FEsolving_execElem(2))) &
                     .and. NiterationState < nState)
              NiterationState = NiterationState + 1_pInt

              do ee = FEsolving_execElem(1),FEsolving_execElem(2)
                myNgrains = homogenization_Ngrains(mesh_element(3,ee))
                do ii = FEsolving_execIP(1,ee),FEsolving_execIP(2,ee)
                  do gg = 1,myNgrains
                    if (crystallite_todo(gg,ii,ee)) &
                      crystallite_onTrack(gg,ii,ee) = crystallite_integrateStress(gg,ii,ee)             ! stress integration
               enddo; enddo; enddo
                    
              crystallite_todo = crystallite_todo .and. crystallite_onTrack                             ! continue with non-broken grains
              if (any(.not. crystallite_onTrack .and. .not. crystallite_localConstitution)) &           ! any non-local is broken?
                  crystallite_todo = crystallite_todo .and. crystallite_localConstitution               ! all nonlocal crystallites can be skipped
        
              do ee = FEsolving_execElem(1),FEsolving_execElem(2)
                myNgrains = homogenization_Ngrains(mesh_element(3,ee))
                do ii = FEsolving_execIP(1,ee),FEsolving_execIP(2,ee)
                  do gg = 1,myNgrains
                    if (crystallite_todo(gg,ii,ee)) &
                      constitutive_dotState(gg,ii,ee)%p = 0.0_pReal                                        ! zero out dotState
              enddo; enddo; enddo

              do ee = FEsolving_execElem(1),FEsolving_execElem(2)
                myNgrains = homogenization_Ngrains(mesh_element(3,ee))
                do ii = FEsolving_execIP(1,ee),FEsolving_execIP(2,ee)
                  do gg = 1,myNgrains
                    if (crystallite_todo(gg,ii,ee)) &
                      call constitutive_collectDotState(crystallite_Tstar_v(:,gg,ii,ee), crystallite_subTstar0_v(:,gg,ii,ee), &
                                                        crystallite_Fe, crystallite_Fp, crystallite_Temperature(gg,ii,ee), & 
                                                        crystallite_misorientation(:,:,g,i,e), crystallite_subdt(gg,ii,ee), &
                                                        gg, ii, ee)                                        ! collect dot state                
              enddo; enddo; enddo

              do ee = FEsolving_execElem(1),FEsolving_execElem(2)
                myNgrains = homogenization_Ngrains(mesh_element(3,ee))
                do ii = FEsolving_execIP(1,ee),FEsolving_execIP(2,ee)
                  do gg = 1,myNgrains
                    if (crystallite_todo(gg,ii,ee)) then
                      crystallite_stateConverged(gg,ii,ee) = crystallite_updateState(gg,ii,ee)             ! update state
                      crystallite_temperatureConverged(gg,ii,ee) = crystallite_updateTemperature(gg,ii,ee) ! update temperature
                      crystallite_converged(gg,ii,ee) =       crystallite_stateConverged(gg,ii,ee) &
                                                        .and. crystallite_temperatureConverged(gg,ii,ee)
                    endif
                  enddo
                enddo
              enddo
        
              if (any(.not. crystallite_converged .and. .not. crystallite_localConstitution)) &            ! any non-local not yet converged?
                  crystallite_converged = crystallite_converged .and. crystallite_localConstitution        ! all non-local not converged
                 
              crystallite_todo = crystallite_todo .and. .not. crystallite_converged                        ! skip all converged
                              
            enddo       ! state loop
            
            if (all(crystallite_converged(:,:,FEsolving_execELem(1):FEsolving_execElem(2)))) then
              crystallite_dPdF(:,:,k,l,g,i,e) = (crystallite_P(:,:,g,i,e) - storedP(:,:,g,i,e))/pert_Fg    ! tangent dP_ij/dFg_kl
            else                                                                                           ! grain did not converge
              crystallite_dPdF(:,:,k,l,g,i,e) = crystallite_fallbackdPdF(:,:,k,l,g,i,e)                    ! use (elastic) fallback
            endif 
            
            do ee = FEsolving_execElem(1),FEsolving_execElem(2)
              myNgrains = homogenization_Ngrains(mesh_element(3,ee))
              do ii = FEsolving_execIP(1,ee),FEsolving_execIP(2,ee)
                do gg = 1,myNgrains
                  mySizeState = constitutive_sizeState(gg,ii,ee)
                  mySizeDotState = constitutive_sizeDotState(gg,ii,ee)
                  constitutive_state(gg,ii,ee)%p = storedState(1:mySizeState,gg,ii,ee)
                  constitutive_dotState(gg,ii,ee)%p = storedDotState(1:mySizeDotState,gg,ii,ee)
            enddo; enddo; enddo
						crystallite_Temperature = storedTemperature
						crystallite_subF = storedF
						crystallite_Fp = storedFp 
						crystallite_invFp = storedInvFp
						crystallite_Fe = storedFe
						crystallite_Lp = storedLp
						crystallite_Tstar_v = storedTstar_v
						crystallite_P = storedP
            
            !$OMP CRITICAL (out)
              debug_StiffnessStateLoopDistribution(NiterationState) = debug_StiffnessstateLoopDistribution(NiterationState) + 1
            !$OMPEND CRITICAL (out)
              
      enddo; enddo; enddo         ! element,ip,grain loop (e,i,g)
      crystallite_converged = storedConvergenceFlag
      
    endif

  endif                           ! jacobian calculation
 
endsubroutine



!********************************************************************
! update the internal state of the constitutive law
! and tell whether state has converged
!********************************************************************
 function crystallite_updateState(&
   g,&              ! grain number
   i,&              ! integration point number
   e &              ! element number
 )
 
 !*** variables and functions from other modules ***!
 use prec, only:                      pReal, &
                                      pInt, &
                                      pLongInt
 use numerics, only:                  rTol_crystalliteState
 use constitutive, only:              constitutive_dotState, &
                                      constitutive_previousDotState, &
                                      constitutive_sizeDotState, &
                                      constitutive_subState0, &
                                      constitutive_state, &
                                      constitutive_relevantState, &
                                      constitutive_microstructure
 use debug, only:                     debugger, &
                                      selectiveDebugger
 use FEsolving, only:                 cycleCounter, theInc
 
 !*** input variables ***!
 integer(pInt), intent(in)::          e, &                          ! element index
                                      i, &                          ! integration point index
                                      g                             ! grain index
 
 !*** output variables ***!
 logical                              crystallite_updateState       ! flag indicating if integration suceeded

 !*** local variables ***!
 real(pReal), dimension(constitutive_sizeDotState(g,i,e)) :: residuum ! residuum from evolution of microstructure
 integer(pInt)                        mySize

 
 mySize = constitutive_sizeDotState(g,i,e)

 ! correct my dotState 
 constitutive_dotState(g,i,e)%p(1:mySize) = constitutive_dotState(g,i,e)%p(1:mySize) * crystallite_statedamper &
                                          + constitutive_previousDotState(g,i,e)%p(1:mySize) * (1.0_pReal-crystallite_statedamper)
 ! calculate the residuum
 residuum = constitutive_state(g,i,e)%p(1:mySize) - constitutive_subState0(g,i,e)%p(1:mySize) &
                                                  - constitutive_dotState(g,i,e)%p(1:mySize) * crystallite_subdt(g,i,e)
 
 ! if NaN occured then return without changing the state
 if (any(residuum/=residuum)) then
   crystallite_updateState = .false.                                  ! indicate state update failed
   crystallite_onTrack(g,i,e) = .false.                               ! no need to calculate any further
   if (debugger) then
     !$OMP CRITICAL (write2out)
       write(6,*) '::: updateState encountered NaN',g,i,e
     !$OMPEND CRITICAL (write2out)
   endif   
   return
 endif
 
 ! update the microstructure
 constitutive_state(g,i,e)%p(1:mySize) = constitutive_state(g,i,e)%p(1:mySize) - residuum
 call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Tstar_v(:,g,i,e), crystallite_Fe, crystallite_Fp, &
                                  g, i, e)
 
 
 ! setting flag to true if state is below relative tolerance, otherwise set it to false
 crystallite_updateState = all(    constitutive_state(g,i,e)%p(1:mySize) < constitutive_relevantState(g,i,e)%p(1:mySize) &
                              .or. abs(residuum) < rTol_crystalliteState*abs(constitutive_state(g,i,e)%p(1:mySize)))
 if (selectiveDebugger) then
   !$OMP CRITICAL (write2out)
     if (crystallite_updateState) then
       write(6,*) '::: updateState converged',g,i,e
     else
       write(6,*) '::: updateState did not converge',g,i,e
     endif
     write(6,*)
     write(6,'(a,f6.1)') 'crystallite_statedamper',crystallite_statedamper
     write(6,*)
     write(6,'(a,/,12(e12.5,x))') 'dotState',constitutive_dotState(g,i,e)%p(1:mySize)
     write(6,*)
     write(6,'(a,/,12(e12.5,x))') 'new state',constitutive_state(g,i,e)%p(1:mySize)
     write(6,*)
     write(6,'(a,/,12(f12.1,x))') 'resid tolerance',abs(residuum/rTol_crystalliteState/constitutive_state(g,i,e)%p(1:mySize))
     write(6,*)
   !$OMPEND CRITICAL (write2out)
 endif
 return

 endfunction


!********************************************************************
! update the temperature of the grain
! and tell whether it has converged
!********************************************************************
 function crystallite_updateTemperature(&
   g,&              ! grain number
   i,&              ! integration point number
   e &              ! element number
 )
 
 !*** variables and functions from other modules ***!
 use prec, only:                      pReal, &
                                      pInt, &
                                      pLongInt
 use numerics, only:                  rTol_crystalliteTemperature
 use constitutive, only:              constitutive_dotTemperature
 use debug, only:                     debugger, &
                                      debug_cumDotTemperatureCalls, &
                                      debug_cumDotTemperatureTicks
 
 !*** input variables ***!
 integer(pInt), intent(in)::          e, &                          ! element index
                                      i, &                          ! integration point index
                                      g                             ! grain index
 
 !*** output variables ***!
 logical                              crystallite_updateTemperature ! flag indicating if integration suceeded

 !*** local variables ***!
 real(pReal) residuum                                               ! residuum from evolution of temperature
 integer(pLongInt)                    tick, &
                                      tock, &
                                      tickrate, &
                                      maxticks
 
 ! calculate the residuum 
 call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
 residuum = crystallite_Temperature(g,i,e) - crystallite_subTemperature0(g,i,e) - &
            crystallite_subdt(g,i,e) * &
            constitutive_dotTemperature(crystallite_Tstar_v(:,g,i,e),crystallite_Temperature(g,i,e),g,i,e)
 call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
 debug_cumDotTemperatureCalls = debug_cumDotTemperatureCalls + 1_pInt
 debug_cumDotTemperatureTicks = debug_cumDotTemperatureTicks + tock-tick
 if (tock < tick) debug_cumDotTemperatureTicks  = debug_cumDotTemperatureTicks + maxticks
 
 ! if NaN occured then return without changing the state
 if (residuum/=residuum) then
   crystallite_updateTemperature = .false.                                  ! indicate update failed
   !$OMP CRITICAL (write2out)
   write(6,*) '::: updateTemperature encountered NaN',g,i,e
   !$OMPEND CRITICAL (write2out)
   return
 endif
 
 ! update the microstructure
 crystallite_Temperature(g,i,e) = crystallite_Temperature(g,i,e) - residuum
 
 ! setting flag to true if residuum is below relative tolerance (or zero Kelvin), otherwise set it to false
 crystallite_updateTemperature = crystallite_Temperature(g,i,e) == 0.0_pReal .or. &
                                  abs(residuum) < rTol_crystalliteTemperature*crystallite_Temperature(g,i,e)
 
 return

 endfunction



!***********************************************************************
!***     calculation of stress (P) with time integration             ***
!***     based on a residuum in Lp and intermediate                  ***
!***     acceleration of the Newton-Raphson correction               ***
!***********************************************************************
 function crystallite_integrateStress(&
     g,&          ! grain number
     i,&          ! integration point number
     e)           ! element number

 !*** variables and functions from other modules ***!
 use prec, only:                      pReal, &
                                      pInt, &
                                      pLongInt
 use numerics, only:                  nStress, &
                                      aTol_crystalliteStress, &
                                      rTol_crystalliteStress, &
                                      iJacoLpresiduum, &
                                      relevantStrain
 use debug, only:                     debugger, &
                                      selectiveDebugger, &
                                      debug_cumLpCalls, &
                                      debug_cumLpTicks, &
                                      debug_StressLoopDistribution
 use constitutive, only:              constitutive_homogenizedC, &
                                      constitutive_LpAndItsTangent
 use math, only:                      math_mul33x33, &
                                      math_mul66x6, &
                                      math_mul99x99, &
                                      math_inv3x3, &
                                      math_invert3x3, &
                                      math_invert, &
                                      math_det3x3, &
                                      math_I3, &
                                      math_identity2nd, &
                                      math_Mandel66to3333, &
                                      math_Mandel6to33, &
                                      math_mandel33to6

 implicit none

 !*** input variables ***!
 integer(pInt), intent(in)::          e, &                          ! element index
                                      i, &                          ! integration point index
                                      g                             ! grain index
 
 !*** output variables ***!
 logical                              crystallite_integrateStress   ! flag indicating if integration suceeded
 
 !*** local variables ***!
 real(pReal), dimension(3,3)::        Fg_new, &                     ! deformation gradient at end of timestep
                                      Fp_current, &                 ! plastic deformation gradient at start of timestep
                                      Fp_new, &                     ! plastic deformation gradient at end of timestep
                                      Fe_new, &                     ! elastic deformation gradient at end of timestep
                                      invFp_new, &                  ! inverse of Fp_new
                                      invFp_current, &              ! inverse of Fp_current
                                      Lpguess, &                    ! current guess for plastic velocity gradient
                                      Lpguess_old, &                ! known last good guess for plastic velocity gradient
                                      Lp_constitutive, &            ! plastic velocity gradient resulting from constitutive law
                                      residuum, &                   ! current residuum of plastic velocity gradient
                                      residuum_old, &               ! last residuum of plastic velocity gradient
                                      A, &
                                      B, &
                                      BT, &
                                      AB, &
                                      BTA
 real(pReal), dimension(6)::          Tstar_v                       ! 2nd Piola-Kirchhoff Stress in Mandel-Notation
 real(pReal), dimension(9,9)::        dLpdT_constitutive, &         ! partial derivative of plastic velocity gradient calculated by constitutive law
                                      dTdLp, &                      ! partial derivative of 2nd Piola-Kirchhoff stress
                                      dRdLp, &                      ! partial derivative of residuum (Jacobian for NEwton-Raphson scheme)
                                      invdRdLp                      ! inverse of dRdLp
 real(pReal), dimension(3,3,3,3)::    C                             ! 4th rank elasticity tensor
 real(pReal), dimension(6,6)::        C_66                          ! simplified 2nd rank elasticity tensor 
 real(pReal)                          p_hydro, &                    ! volumetric part of 2nd Piola-Kirchhoff Stress
                                      det, &                        ! determinant
                                      leapfrog, &                   ! acceleration factor for Newton-Raphson scheme
                                      maxleap                       ! maximum acceleration factor
 logical                              error                         ! flag indicating an error
 integer(pInt)                        NiterationStress, &           ! number of stress integrations
                                      dummy, &
                                      h, &
                                      j, &
                                      k, &
                                      l, &
                                      m, &
                                      n, &
                                      jacoCounter                   ! counter to check for Jacobian update
 integer(pLongInt)                    tick, &
                                      tock, &
                                      tickrate, &
                                      maxticks
 
 ! be pessimistic
 crystallite_integrateStress = .false.

 ! feed local variables
 Fg_new =       crystallite_subF(:,:,g,i,e)
 Fp_current =   crystallite_subFp0(:,:,g,i,e)
 Tstar_v =      crystallite_Tstar_v(:,g,i,e)
 Lpguess_old =  crystallite_Lp(:,:,g,i,e)                           ! consider present Lp good (i.e. worth remembering) ...
 Lpguess =      crystallite_Lp(:,:,g,i,e)                           ! ... and take it as first guess

 
 ! inversion of Fp_current...
 invFp_current = math_inv3x3(Fp_current)                            
 if (all(invFp_current == 0.0_pReal)) then                          ! ... failed?
   if (debugger) then 
     !$OMP CRITICAL (write2out)
       write(6,*) '::: integrateStress failed on invFp_current inversion',g,i,e
       write(6,*)
       write(6,'(a11,3(i3,x),/,3(3(f12.7,x)/))') 'invFp_new at ',g,i,e,invFp_new
     !$OMPEND CRITICAL (write2out)
   endif
   return
 endif
 
 A = math_mul33x33(transpose(invFp_current), math_mul33x33(transpose(Fg_new),math_mul33x33(Fg_new,invFp_current)))
 
 ! get elasticity tensor
 C_66 = constitutive_homogenizedC(g,i,e)
! if (debugger) write(6,'(a,/,6(6(f10.4,x)/))') 'elasticity',C_66(1:6,:)/1e9
 C = math_Mandel66to3333(C_66)
 
 ! start LpLoop with no acceleration
 NiterationStress = 0_pInt
 leapfrog = 1.0_pReal
 maxleap = 1024.0_pReal
 jacoCounter = 0_pInt

LpLoop: do
   
   ! increase loop counter
   NiterationStress = NiterationStress + 1
   
   ! too many loops required ?
   if (NiterationStress > nStress) then
     if (debugger) then 
       !$OMP CRITICAL (write2out)
         write(6,*) '::: integrateStress reached loop limit at ',g,i,e
         write(6,*)
       !$OMPEND CRITICAL (write2out)
     endif
     return
   endif
   
   B = math_I3 - crystallite_subdt(g,i,e)*Lpguess
   BT = transpose(B)
   AB = math_mul33x33(A,B)
   BTA = math_mul33x33(BT,A)
   
   ! calculate 2nd Piola-Kirchhoff stress tensor
   Tstar_v = 0.5_pReal*math_mul66x6(C_66,math_mandel33to6(math_mul33x33(BT,AB)-math_I3))
   p_hydro = sum(Tstar_v(1:3))/3.0_pReal
   forall(n=1:3) Tstar_v(n) = Tstar_v(n) - p_hydro                  ! get deviatoric stress tensor
   
   ! calculate plastic velocity gradient and its tangent according to constitutive law
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
   call constitutive_LpAndItsTangent(Lp_constitutive, dLpdT_constitutive, Tstar_v, crystallite_Temperature(g,i,e), g, i, e)
   call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
   debug_cumLpCalls = debug_cumLpCalls + 1_pInt
   debug_cumLpTicks  = debug_cumLpTicks + tock-tick
   if (tock < tick) debug_cumLpTicks = debug_cumLpTicks + maxticks
   if (selectiveDebugger) then
     !$OMP CRITICAL (write2out)
       write(6,*) '::: integrateStress at ' ,g,i,e, ' ; iteration ', NiterationStress
       write(6,*)
       write(6,'(a,/,3(3(f20.7,x)/))') 'Lp_constitutive', Lp_constitutive
       write(6,'(a,/,3(3(f20.7,x)/))') 'Lpguess', Lpguess
     !$OMPEND CRITICAL (write2out)
   endif

   ! update current residuum
   residuum = Lpguess - Lp_constitutive

   ! Check for convergence of loop
   if (.not.(any(residuum/=residuum)) .and. &                       ! exclude any NaN in residuum
       ( maxval(abs(residuum)) < aTol_crystalliteStress .or. &      ! below absolute tolerance .or.
         ( any(abs(crystallite_subdt(g,i,e)*Lpguess) > relevantStrain) .and. & ! worth checking? .and.
             maxval(abs(residuum/Lpguess), abs(crystallite_subdt(g,i,e)*Lpguess) > relevantStrain) < rTol_crystalliteStress & ! below relative tolerance
         ) &
       ) &
      ) &
     exit LpLoop
   
   ! NaN occured at regular speed?
   if (any(residuum/=residuum) .and. leapfrog == 1.0) then
     if (debugger) then 
       !$OMP CRITICAL (write2out)
         write(6,*) '::: integrateStress encountered NaN at ',g,i,e,' ; iteration ', NiterationStress
       !$OMPEND CRITICAL (write2out)
     endif
     return

   ! something went wrong at accelerated speed?
   elseif (leapfrog > 1.0_pReal .and. &                             ! at fast pace .and.
            ( sum(residuum*residuum) > sum(residuum_old*residuum_old) .or. &  ! worse residuum .or.
              sum(residuum*residuum_old) < 0.0_pReal .or. &         ! residuum changed sign (overshoot) .or.
              any(residuum/=residuum) &                             ! NaN occured
            ) &
          ) then
     maxleap = 0.5_pReal * leapfrog                                 ! limit next acceleration
     leapfrog = 1.0_pReal                                           ! grinding halt
     jacoCounter = 0_pInt                                           ! reset counter for Jacobian update (we want to do an update next time!)
     
     ! restore old residuum and Lp
     Lpguess = Lpguess_old                                       
     residuum  = residuum_old

   ! residuum got better
   else
     ! calculate Jacobian for correction term 
     if (mod(jacoCounter, iJacoLpresiduum) == 0_pInt) then
       dTdLp = 0.0_pReal
       forall (h=1:3,j=1:3,k=1:3,l=1:3,m=1:3) &
         dTdLp(3*(h-1)+j,3*(k-1)+l) = dTdLp(3*(h-1)+j,3*(k-1)+l) + C(h,j,l,m)*AB(k,m)+C(h,j,m,l)*BTA(m,k)
       dTdLp = -0.5_pReal*crystallite_subdt(g,i,e)*dTdLp
       dRdLp = math_identity2nd(9) - math_mul99x99(dLpdT_constitutive,dTdLp)
       invdRdLp = 0.0_pReal
       call math_invert(9,dRdLp,invdRdLp,dummy,error)               ! invert dR/dLp --> dLp/dR
       if (error) then
         if (debugger) then
           !$OMP CRITICAL (write2out)
             write(6,*) '::: integrateStress failed on dR/dLp inversion at ',g,i,e,' ; iteration ', NiterationStress
             write(6,*)
             write(6,'(a,/,9(9(f15.3,x)/))') 'dRdLp',dRdLp
             write(6,'(a,/,9(9(f15.3,x)/))') 'dLpdT_constitutive',dLpdT_constitutive
             write(6,'(a,/,3(3(f20.7,x)/))') 'Lp_constitutive',Lp_constitutive
             write(6,'(a,/,3(3(f20.7,x)/))') 'Lpguess',Lpguess
           !$OMPEND CRITICAL (write2out)
         endif
         return
       endif
     endif
     jacoCounter = jacoCounter + 1_pInt                             ! increase counter for jaco update
     
     ! remember current residuum and Lpguess
     residuum_old = residuum
     Lpguess_old = Lpguess 
     
     ! accelerate?
     if (NiterationStress > 1 .and. leapfrog < maxleap) leapfrog = 2.0_pReal * leapfrog
   endif

   ! leapfrog to updated Lp
   forall (k=1:3,l=1:3,m=1:3,n=1:3) & 
     Lpguess(k,l) = Lpguess(k,l) - leapfrog*invdRdLp(3*(k-1)+l,3*(m-1)+n)*residuum(m,n)
 enddo LpLoop

 ! calculate new plastic and elastic deformation gradient
 invFp_new = math_mul33x33(invFp_current,B)
 invFp_new = invFp_new/math_det3x3(invFp_new)**(1.0_pReal/3.0_pReal)  ! regularize by det
 call math_invert3x3(invFp_new,Fp_new,det,error)
 if (error) then
   if (debugger) then
     !$OMP CRITICAL (write2out)
       write(6,*) '::: integrateStress failed on invFp_new inversion at ',g,i,e,' ; iteration ', NiterationStress
       write(6,*)
       write(6,'(a11,3(i3,x),/,3(3(f12.7,x)/))') 'invFp_new at ',g,i,e,invFp_new
     !$OMPEND CRITICAL (write2out)
   endif
   return
 endif
 Fe_new = math_mul33x33(Fg_new,invFp_new)                             ! calc resulting Fe

 ! add volumetric component to 2nd Piola-Kirchhoff stress
 forall (n=1:3) Tstar_v(n) = Tstar_v(n) + p_hydro
 
 ! calculate 1st Piola-Kirchhoff stress
 crystallite_P(:,:,g,i,e) = math_mul33x33(Fe_new,math_mul33x33(math_Mandel6to33(Tstar_v),transpose(invFp_new)))
 
 ! store local values in global variables
 crystallite_Lp(:,:,g,i,e) = Lpguess
 crystallite_Tstar_v(:,g,i,e) = Tstar_v
 crystallite_Fp(:,:,g,i,e) = Fp_new
 crystallite_Fe(:,:,g,i,e) = Fe_new
 crystallite_invFp(:,:,g,i,e) = invFp_new

 ! set return flag to true
 crystallite_integrateStress = .true.
 if (selectiveDebugger) then 
   !$OMP CRITICAL (write2out)
   write(6,*) '::: integrateStress converged at ',g,i,e,' ; iteration ', NiterationStress
   write(6,*)
   write(6,'(a,/,3(3(f12.7,x)/))') 'P / MPa',crystallite_P(:,:,g,i,e)/1e6
   write(6,'(a,/,3(3(f12.7,x)/))') 'Lp',crystallite_Lp(:,:,g,i,e)
   write(6,'(a,/,3(3(f12.7,x)/))') 'Fp',crystallite_Fp(:,:,g,i,e)
   !$OMP CRITICAL (write2out)
 endif

 !$OMP CRITICAL (distributionStress)
 debug_StressLoopDistribution(NiterationStress) = debug_StressLoopDistribution(NiterationStress) + 1
 !$OMPEND CRITICAL (distributionStress)

 return

 endfunction
 
 
 
!********************************************************************
! calculates orientations and misorientations (in case of single grain ips)
!******************************************************************** 
subroutine crystallite_orientations()
  
!*** variables and functions from other modules ***!
use prec, only:												pInt, &
																			pReal
use math, only:               				math_pDecomposition, &
																			math_RtoEuler, &
																			math_misorientation, &
																			inDeg
use FEsolving, only:          				FEsolving_execElem, & 
														  				FEsolving_execIP
use IO, only:                 				IO_warning
use material, only:           				material_phase, &
																			homogenization_Ngrains, &
																			phase_constitution, &
                              				phase_constitutionInstance
use mesh, only:               				mesh_element, &
                              				mesh_ipNeighborhood, &
														  				FE_NipNeighbors
use constitutive_phenopowerlaw, only: constitutive_phenopowerlaw_label, &
																			constitutive_phenopowerlaw_structure
use constitutive_dislotwin, only: 		constitutive_dislotwin_label, &
																			constitutive_dislotwin_structure
use constitutive_nonlocal, only: 			constitutive_nonlocal_label, &
																			constitutive_nonlocal_structure

implicit none

!*** input variables ***!

!*** output variables ***!

!*** local variables ***!
integer(pInt)                 	e, &                          ! element index
																i, &                          ! integration point index
																g, &                          ! grain index
																n, &													! neighbor index 
																myPhase, &									  ! phase
																myStructure, &                ! lattice structure
																neighboring_e, &              ! element index of my neighbor
																neighboring_i, &              ! integration point index of my neighbor
																neighboringPhase, &           ! phase of my neighbor
																neighboringStructure, &       ! lattice structure of my neighbor
																symmetryType                  ! type of crystal symmetry
real(pReal), dimension(3,3) ::  U, R, &											  ! polar decomposition of Fe
																netRotation                   ! net rotation between two orientations
logical error


!$OMP PARALLEL DO
	do e = FEsolving_execElem(1),FEsolving_execElem(2)
		do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
			do g = 1,homogenization_Ngrains(mesh_element(3,e))
			  
			  ! calculate orientation in terms of rotation matrix and euler angles
				call math_pDecomposition(crystallite_Fe(:,:,g,i,e), U, R, error)          		! polar decomposition of Fe
				if (error) then
					call IO_warning(650, e, i, g)
					crystallite_R(:,:,g,i,e) = 0.0_pReal
					crystallite_eulerangles(:,g,i,e) = (/400.0, 400.0, 400.0/)                 		! fake orientation
				else
					crystallite_R(:,:,g,i,e) = transpose(R)
					crystallite_eulerangles(:,g,i,e) = math_RtoEuler(crystallite_R(:,:,g,i,e)) * inDeg
				endif
				
			enddo
		enddo
	enddo
!$OMPEND PARALLEL DO

!$OMP PARALLEL DO
	do e = FEsolving_execElem(1),FEsolving_execElem(2)
		do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
		  if (homogenization_Ngrains(mesh_element(3,e)) == 1_pInt) then                   ! if single grain ip

        myPhase = material_phase(1,i,e)                																! get my crystal structure
				select case (phase_constitution(myPhase))
					case (constitutive_phenopowerlaw_label)
						myStructure = constitutive_phenopowerlaw_structure(phase_constitutionInstance(myPhase))
					case (constitutive_dislotwin_label)
						myStructure = constitutive_dislotwin_structure(phase_constitutionInstance(myPhase))
					case (constitutive_nonlocal_label)
						myStructure = constitutive_nonlocal_structure(phase_constitutionInstance(myPhase))
					case default
						myStructure = ''
				end select   
			  
			  do n = 1,FE_NipNeighbors(mesh_element(2,e))                                   ! loop through my neighbors
					
					neighboring_e = mesh_ipNeighborhood(1,n,i,e)
					neighboring_i = mesh_ipNeighborhood(2,n,i,e)

					if ((neighboring_e > 0) .and. (neighboring_i > 0)) then										  ! if neighbor exists

						neighboringPhase = material_phase(1,neighboring_i,neighboring_e)        	! get my neighbor's crystal structure								
						if (myPhase == neighboringPhase) then                           	        ! if my neighbor has same phase like me
						
							select case (myStructure)                                               ! get type of symmetry
								case (1_pInt, 2_pInt)                                                 ! fcc and bcc:
									symmetryType = 1_pInt																								! -> cubic symmetry
								case (3_pInt)																													! hex:
									symmetryType = 2_pInt																								! -> hexagonal symmetry
								case default
									symmetryType = 0_pInt
							end select
							
							call math_misorientation( crystallite_misorientation(1:3,n,1,i,e), &
																				crystallite_misorientation(4,n,1,i,e), &
																				netRotation, &
																				crystallite_R(:,:,1,i,e), &
																				crystallite_R(:,:,1,neighboring_i,neighboring_e), & 
																				symmetryType)   																! calculate misorientation
				    
				    else                                                                      ! for neighbor with different phase
				      crystallite_misorientation(4,n,1,i,e) = 400.0_pReal                     ! set misorientation angle to 400
				      
				    endif
				  else                                                                        ! no existing neighbor
				    crystallite_misorientation(4,n,1,i,e) = 0.0_pReal                         ! set misorientation angle to zero
				  endif
				enddo
			endif
		enddo
	enddo
!$OMPEND PARALLEL DO

endsubroutine


 
!********************************************************************
! return results of particular grain
!********************************************************************
function crystallite_postResults(&
   dt,&             ! time increment
   g,&              ! grain number
   i,&              ! integration point number
   e &              ! element number
 )

 !*** variables and functions from other modules ***!
 use prec, only:                      pInt, &
                                      pReal
 use material, only:                  material_phase, &
                                      material_volume
 use constitutive, only:              constitutive_sizePostResults, &
                                      constitutive_postResults
 
 implicit none

 !*** input variables ***!
 integer(pInt), intent(in)::          e, &                          ! element index
                                      i, &                          ! integration point index
                                      g                             ! grain index
 real(pReal), intent(in)::            dt                            ! time increment

 !*** output variables ***!
 real(pReal), dimension(1+crystallite_Nresults + 1+constitutive_sizePostResults(g,i,e)) :: crystallite_postResults
 
 !*** local variables ***!
 real(pReal), dimension(3,3) ::       U, R
 integer(pInt)              k,l,c
 logical error

 c = 0_pInt 
 crystallite_postResults(c+1) = crystallite_Nresults; c = c+1_pInt         ! size of (hardwired) results
 if (crystallite_Nresults >= 2) then
   crystallite_postResults(c+1) = material_phase(g,i,e)
   crystallite_postResults(c+2) = material_volume(g,i,e)
   c = c+2_pInt
 endif
 if (crystallite_Nresults >= 5) then
   crystallite_postResults(c+1:c+3) = crystallite_eulerangles(:,i,e,g)                 ! fake orientation
   c = c+3_pInt
 endif
 if (crystallite_Nresults >= 14) then                                     ! deformation gradient
   forall (k=0:2,l=0:2) crystallite_postResults(c+1+k*3+l) = crystallite_partionedF(k+1,l+1,g,i,e)
   c = c+9_pInt
 endif
 
 crystallite_postResults(c+1) = constitutive_sizePostResults(g,i,e); c = c+1_pInt  ! size of constitutive results
 crystallite_postResults(c+1:c+constitutive_sizePostResults(g,i,e)) = &
         constitutive_postResults(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), crystallite_Fe, crystallite_Fp, &
                                  crystallite_Temperature(g,i,e), crystallite_misorientation(:,:,g,i,e), dt, &
                                  crystallite_subdt(g,i,e), g, i, e)
 c = c + constitutive_sizePostResults(g,i,e)
 
 return 
 
endfunction


END MODULE
!##############################################################

