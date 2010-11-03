!* $Id$
!************************************ 
!*  Module: CONSTITUTIVE_NONLOCAL   *
!************************************
!* contains:                        *
!* - constitutive equations         *
!* - parameters definition          *
!************************************


MODULE constitutive_nonlocal

!* Include other modules
use prec, only: pReal,pInt
implicit none


!* Definition of parameters
character (len=*), parameter :: constitutive_nonlocal_label = 'nonlocal'
character(len=22), dimension(10), parameter ::    constitutive_nonlocal_listBasicStates = (/'rhoSglEdgePosMobile   ', &
                                                                                            'rhoSglEdgeNegMobile   ', &
                                                                                            'rhoSglScrewPosMobile  ', &
                                                                                            'rhoSglScrewNegMobile  ', &
                                                                                            'rhoSglEdgePosImmobile ', &
                                                                                            'rhoSglEdgeNegImmobile ', &
                                                                                            'rhoSglScrewPosImmobile', &
                                                                                            'rhoSglScrewNegImmobile', &
                                                                                            'rhoDipEdge            ', &
                                                                                            'rhoDipScrew           ' /) ! list of "basic" microstructural state variables that are independent from other state variables
character(len=15), dimension(3), parameter :: constitutive_nonlocal_listDependentStates = (/'rhoForest      ', &
                                                                                            'tauThreshold   ', &
                                                                                            'Tdislocation_v ' /) ! list of microstructural state variables that depend on other state variables
real(pReal), parameter :: kB = 1.38e-23_pReal                                                                   ! Physical parameter, Boltzmann constant in J/Kelvin

!* Definition of global variables
integer(pInt), dimension(:), allocatable ::               constitutive_nonlocal_sizeDotState, &                 ! number of dotStates
                                                          constitutive_nonlocal_sizeState, &                    ! total number of microstructural state variables
                                                          constitutive_nonlocal_sizePostResults                 ! cumulative size of post results
integer(pInt), dimension(:,:), allocatable, target ::     constitutive_nonlocal_sizePostResult                  ! size of each post result output
character(len=64), dimension(:,:), allocatable, target :: constitutive_nonlocal_output                          ! name of each post result output

character(len=32), dimension(:), allocatable ::           constitutive_nonlocal_structureName                   ! name of the lattice structure
integer(pInt), dimension(:), allocatable ::               constitutive_nonlocal_structure, &                    ! number representing the kind of lattice structure
                                                          constitutive_nonlocal_totalNslip                      ! total number of active slip systems for each instance
integer(pInt), dimension(:,:), allocatable ::             constitutive_nonlocal_Nslip, &                        ! number of active slip systems for each family and instance
                                                          constitutive_nonlocal_slipFamily, &                   ! lookup table relating active slip system to slip family for each instance
                                                          constitutive_nonlocal_slipSystemLattice               ! lookup table relating active slip system index to lattice slip system index for each instance

real(pReal), dimension(:), allocatable ::                 constitutive_nonlocal_CoverA, &                       ! c/a ratio for hex type lattice
                                                          constitutive_nonlocal_C11, &                          ! C11 element in elasticity matrix
                                                          constitutive_nonlocal_C12, &                          ! C12 element in elasticity matrix
                                                          constitutive_nonlocal_C13, &                          ! C13 element in elasticity matrix
                                                          constitutive_nonlocal_C33, &                          ! C33 element in elasticity matrix
                                                          constitutive_nonlocal_C44, &                          ! C44 element in elasticity matrix
                                                          constitutive_nonlocal_Gmod, &                         ! shear modulus
                                                          constitutive_nonlocal_nu, &                           ! poisson's ratio
                                                          constitutive_nonlocal_Q0, &                           ! activation energy for dislocation glide
                                                          constitutive_nonlocal_atomicVolume, &                 ! atomic volume
                                                          constitutive_nonlocal_D0, &                           ! prefactor for self-diffusion coefficient
                                                          constitutive_nonlocal_Qsd, &                          ! activation enthalpy for diffusion
                                                          constitutive_nonlocal_aTolRho, &                      ! absolute tolerance for dislocation density in state integration
                                                          constitutive_nonlocal_R                               ! cutoff radius for dislocation stress
real(pReal), dimension(:,:,:), allocatable ::             constitutive_nonlocal_Cslip_66                        ! elasticity matrix in Mandel notation for each instance
real(pReal), dimension(:,:,:,:,:), allocatable ::         constitutive_nonlocal_Cslip_3333                      ! elasticity matrix for each instance
real(pReal), dimension(:,:), allocatable ::               constitutive_nonlocal_rhoSglEdgePos0, &               ! initial edge_pos dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoSglEdgeNeg0, &               ! initial edge_neg dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoSglScrewPos0, &              ! initial screw_pos dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoSglScrewNeg0, &              ! initial screw_neg dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoDipEdge0, &                  ! initial edge dipole dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoDipScrew0, &                 ! initial screw dipole dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_v0PerSlipFamily, &              ! dislocation velocity prefactor [m/s] for each family and instance
                                                          constitutive_nonlocal_v0PerSlipSystem, &              ! dislocation velocity prefactor [m/s] for each slip system and instance
                                                          constitutive_nonlocal_lambda0PerSlipFamily, &         ! mean free path prefactor for each family and instance
                                                          constitutive_nonlocal_lambda0PerSlipSystem, &         ! mean free path prefactor for each slip system and instance
                                                          constitutive_nonlocal_burgersPerSlipFamily, &         ! absolute length of burgers vector [m] for each family and instance
                                                          constitutive_nonlocal_burgersPerSlipSystem, &         ! absolute length of burgers vector [m] for each slip system and instance
                                                          constitutive_nonlocal_dLowerEdgePerSlipFamily, &      ! minimum stable edge dipole height for each family and instance
                                                          constitutive_nonlocal_dLowerEdgePerSlipSystem, &      ! minimum stable edge dipole height for each slip system and instance
                                                          constitutive_nonlocal_dLowerScrewPerSlipFamily, &     ! minimum stable screw dipole height for each family and instance
                                                          constitutive_nonlocal_dLowerScrewPerSlipSystem, &     ! minimum stable screw dipole height for each slip system and instance
                                                          constitutive_nonlocal_interactionSlipSlip             ! coefficients for slip-slip interaction for each interaction type and instance
real(pReal), dimension(:,:,:,:,:), allocatable ::         constitutive_nonlocal_v, &                            ! dislocation velocity
                                                          constitutive_nonlocal_rhoDotFlux                      ! dislocation convection term
real(pReal), dimension(:,:,:,:,:,:), allocatable ::       constitutive_nonlocal_compatibility                   ! slip system compatibility between me and my neighbors
real(pReal), dimension(:,:,:), allocatable ::             constitutive_nonlocal_forestProjectionEdge, &         ! matrix of forest projections of edge dislocations for each instance
                                                          constitutive_nonlocal_forestProjectionScrew, &        ! matrix of forest projections of screw dislocations for each instance
                                                          constitutive_nonlocal_interactionMatrixSlipSlip       ! interaction matrix of the different slip systems for each instance

CONTAINS
!****************************************
!* - constitutive_init
!* - constitutive_stateInit
!* - constitutive_homogenizedC
!* - constitutive_microstructure
!* - constitutive_LpAndItsTangent
!* - constitutive_dotState
!* - constitutive_dotTemperature
!* - constitutive_postResults
!****************************************


!**************************************
!*      Module initialization         *
!**************************************
subroutine constitutive_nonlocal_init(file)

use prec,     only: pInt, pReal
use math,     only: math_Mandel3333to66, & 
                    math_Voigt66to3333, & 
                    math_mul3x3
use IO,       only: IO_lc, &
                    IO_getTag, &
                    IO_isBlank, &
                    IO_stringPos, &
                    IO_stringValue, &
                    IO_floatValue, &
                    IO_intValue, &
                    IO_error
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    FE_maxNipNeighbors
use material, only: homogenization_maxNgrains, &
                    phase_constitution, &
                    phase_constitutionInstance, &
                    phase_Noutput
use lattice,  only: lattice_maxNslipFamily, &
                    lattice_maxNtwinFamily, &
                    lattice_maxNslip, &
                    lattice_maxNtwin, &
                    lattice_maxNinteraction, &
                    lattice_NslipSystem, &
                    lattice_NtwinSystem, &
                    lattice_initializeStructure, &
                    lattice_Qtwin, &
                    lattice_sd, &
                    lattice_sn, &
                    lattice_st, &
                    lattice_interactionSlipSlip

!*** output variables

!*** input variables
integer(pInt), intent(in) ::                file

!*** local variables
integer(pInt), parameter ::                 maxNchunks = 21
integer(pInt), dimension(1+2*maxNchunks) :: positions
integer(pInt)                               section, &
                                            maxNinstance, &
                                            maxTotalNslip, &
                                            myStructure, &
                                            f, &                ! index of my slip family
                                            i, &                ! index of my instance of this constitution
                                            j, &
                                            k, &
                                            l, &
                                            ns, &               ! short notation for total number of active slip systems for the current instance
                                            o, &                ! index of my output
                                            s, &                ! index of my slip system
                                            s1, &               ! index of my slip system
                                            s2, &               ! index of my slip system
                                            it, &               ! index of my interaction type
                                            output, &
                                            mySize
character(len=64)                           tag
character(len=1024)                         line


write(6,*)
write(6,'(a20,a20,a12)') '<<<+-  constitutive_',constitutive_nonlocal_label,' init  -+>>>'
write(6,*) '$Id$'
write(6,*)

maxNinstance = count(phase_constitution == constitutive_nonlocal_label)
if (maxNinstance == 0) return                                                                                                       ! we don't have to do anything if there's no instance for this constitutive law


!*** space allocation for global variables

allocate(constitutive_nonlocal_sizeDotState(maxNinstance))
allocate(constitutive_nonlocal_sizeState(maxNinstance))
allocate(constitutive_nonlocal_sizePostResults(maxNinstance))
allocate(constitutive_nonlocal_sizePostResult(maxval(phase_Noutput), maxNinstance))
allocate(constitutive_nonlocal_output(maxval(phase_Noutput), maxNinstance))
constitutive_nonlocal_sizeDotState = 0_pInt
constitutive_nonlocal_sizeState = 0_pInt
constitutive_nonlocal_sizePostResults = 0_pInt
constitutive_nonlocal_sizePostResult = 0_pInt
constitutive_nonlocal_output = ''

allocate(constitutive_nonlocal_structureName(maxNinstance))
allocate(constitutive_nonlocal_structure(maxNinstance))
allocate(constitutive_nonlocal_Nslip(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_slipFamily(lattice_maxNslip, maxNinstance))
allocate(constitutive_nonlocal_slipSystemLattice(lattice_maxNslip, maxNinstance))
allocate(constitutive_nonlocal_totalNslip(maxNinstance))
constitutive_nonlocal_structureName = ''
constitutive_nonlocal_structure = 0_pInt
constitutive_nonlocal_Nslip = 0_pInt
constitutive_nonlocal_slipFamily = 0_pInt
constitutive_nonlocal_slipSystemLattice = 0_pInt
constitutive_nonlocal_totalNslip = 0_pInt

allocate(constitutive_nonlocal_CoverA(maxNinstance))
allocate(constitutive_nonlocal_C11(maxNinstance))
allocate(constitutive_nonlocal_C12(maxNinstance))
allocate(constitutive_nonlocal_C13(maxNinstance))
allocate(constitutive_nonlocal_C33(maxNinstance))
allocate(constitutive_nonlocal_C44(maxNinstance))
allocate(constitutive_nonlocal_Gmod(maxNinstance))
allocate(constitutive_nonlocal_nu(maxNinstance))
allocate(constitutive_nonlocal_Q0(maxNinstance))
allocate(constitutive_nonlocal_atomicVolume(maxNinstance))
allocate(constitutive_nonlocal_D0(maxNinstance))
allocate(constitutive_nonlocal_Qsd(maxNinstance))
allocate(constitutive_nonlocal_aTolRho(maxNinstance))
allocate(constitutive_nonlocal_Cslip_66(6,6,maxNinstance))
allocate(constitutive_nonlocal_Cslip_3333(3,3,3,3,maxNinstance))
allocate(constitutive_nonlocal_R(maxNinstance))
constitutive_nonlocal_CoverA = 0.0_pReal 
constitutive_nonlocal_C11 = 0.0_pReal
constitutive_nonlocal_C12 = 0.0_pReal
constitutive_nonlocal_C13 = 0.0_pReal
constitutive_nonlocal_C33 = 0.0_pReal
constitutive_nonlocal_C44 = 0.0_pReal
constitutive_nonlocal_Gmod = 0.0_pReal
constitutive_nonlocal_Q0 = 0.0_pReal
constitutive_nonlocal_atomicVolume = 0.0_pReal
constitutive_nonlocal_D0 = 0.0_pReal
constitutive_nonlocal_Qsd = 0.0_pReal
constitutive_nonlocal_aTolRho = 0.0_pReal
constitutive_nonlocal_nu = 0.0_pReal
constitutive_nonlocal_Cslip_66 = 0.0_pReal
constitutive_nonlocal_Cslip_3333 = 0.0_pReal
constitutive_nonlocal_R = 0.0_pReal

allocate(constitutive_nonlocal_rhoSglEdgePos0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoSglEdgeNeg0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoSglScrewPos0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoSglScrewNeg0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoDipEdge0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoDipScrew0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_v0PerSlipFamily(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_burgersPerSlipFamily(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_Lambda0PerSlipFamily(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_interactionSlipSlip(lattice_maxNinteraction, maxNinstance))
allocate(constitutive_nonlocal_dLowerEdgePerSlipFamily(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_dLowerScrewPerSlipFamily(lattice_maxNslipFamily, maxNinstance))
constitutive_nonlocal_rhoSglEdgePos0 = -1.0_pReal
constitutive_nonlocal_rhoSglEdgeNeg0 = -1.0_pReal
constitutive_nonlocal_rhoSglScrewPos0 = -1.0_pReal
constitutive_nonlocal_rhoSglScrewNeg0 = -1.0_pReal
constitutive_nonlocal_rhoDipEdge0 = -1.0_pReal
constitutive_nonlocal_rhoDipScrew0 = -1.0_pReal
constitutive_nonlocal_v0PerSlipFamily = 0.0_pReal
constitutive_nonlocal_burgersPerSlipFamily = 0.0_pReal
constitutive_nonlocal_lambda0PerSlipFamily = 0.0_pReal
constitutive_nonlocal_interactionSlipSlip = 0.0_pReal
constitutive_nonlocal_dLowerEdgePerSlipFamily = 0.0_pReal
constitutive_nonlocal_dLowerScrewPerSlipFamily = 0.0_pReal

!*** readout data from material.config file

rewind(file)
line = ''
section = 0

do while (IO_lc(IO_getTag(line,'<','>')) /= 'phase')                                                                                ! wind forward to <phase>
  read(file,'(a1024)',END=100) line
enddo

do                                                                                                                                  ! read thru sections of phase part
  read(file,'(a1024)',END=100) line
  if (IO_isBlank(line)) cycle                                                                                                       ! skip empty lines
  if (IO_getTag(line,'<','>') /= '') exit                                                                                           ! stop at next part
  if (IO_getTag(line,'[',']') /= '') then                                                                                           ! next section
    section = section + 1
    output = 0                                                                                                                      ! reset output counter
  endif
  if (section > 0 .and. phase_constitution(section) == constitutive_nonlocal_label) then                                            ! one of my sections
    i = phase_constitutionInstance(section)                                                                                         ! which instance of my constitution is present phase
    positions = IO_stringPos(line,maxNchunks)
    tag = IO_lc(IO_stringValue(line,positions,1))                                                                                   ! extract key
    select case(tag)
      case ('(output)')
        output = output + 1
        constitutive_nonlocal_output(output,i) = IO_lc(IO_stringValue(line,positions,2))
      case ('lattice_structure')
        constitutive_nonlocal_structureName(i) = IO_lc(IO_stringValue(line,positions,2))
      case ('covera_ratio')
        constitutive_nonlocal_CoverA(i) = IO_floatValue(line,positions,2)
      case ('c11')
        constitutive_nonlocal_C11(i) = IO_floatValue(line,positions,2)
      case ('c12')
        constitutive_nonlocal_C12(i) = IO_floatValue(line,positions,2)
      case ('c13')
        constitutive_nonlocal_C13(i) = IO_floatValue(line,positions,2)
      case ('c33')
        constitutive_nonlocal_C33(i) = IO_floatValue(line,positions,2)
      case ('c44')
        constitutive_nonlocal_C44(i) = IO_floatValue(line,positions,2)
      case ('nslip')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_Nslip(f,i) = IO_intValue(line,positions,1+f)
      case ('rhosgledgepos0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoSglEdgePos0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhosgledgeneg0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoSglEdgeNeg0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhosglscrewpos0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoSglScrewPos0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhosglscrewneg0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoSglScrewNeg0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhodipedge0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoDipEdge0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhodipscrew0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoDipScrew0(f,i) = IO_floatValue(line,positions,1+f)
      case ('v0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_v0PerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case ('lambda0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_lambda0PerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case ('burgers')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_burgersPerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case('r')
        constitutive_nonlocal_R(i) = IO_floatValue(line,positions,2)
      case('ddipminedge')
        forall (f = 1:lattice_maxNslipFamily) & 
            constitutive_nonlocal_dLowerEdgePerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case('ddipminscrew')
        forall (f = 1:lattice_maxNslipFamily) & 
            constitutive_nonlocal_dLowerScrewPerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case('q0')
        constitutive_nonlocal_Q0(i) = IO_floatValue(line,positions,2)
      case('atomicvolume')
        constitutive_nonlocal_atomicVolume(i) = IO_floatValue(line,positions,2)
      case('d0')
        constitutive_nonlocal_D0(i) = IO_floatValue(line,positions,2)
      case('qsd')
        constitutive_nonlocal_Qsd(i) = IO_floatValue(line,positions,2)
      case('atol_rho')
        constitutive_nonlocal_aTolRho(i) = IO_floatValue(line,positions,2)
      case ('interaction_slipslip')
        forall (it = 1:lattice_maxNinteraction) constitutive_nonlocal_interactionSlipSlip(it,i) = IO_floatValue(line,positions,1+it)
    end select
  endif
enddo


100 do i = 1,maxNinstance

  constitutive_nonlocal_structure(i) = &
    lattice_initializeStructure(constitutive_nonlocal_structureName(i), constitutive_nonlocal_CoverA(i))                            ! our lattice structure is defined in the material.config file by the structureName (and the c/a ratio)
  myStructure = constitutive_nonlocal_structure(i)
  
!*** sanity checks
  
  if (myStructure < 1 .or. myStructure > 3)                                                 call IO_error(205)
  if (sum(constitutive_nonlocal_Nslip(:,i)) <= 0_pInt)                                      call IO_error(225)
  do o = 1,maxval(phase_Noutput)
    if(len(constitutive_nonlocal_output(o,i)) > 64)                                         call IO_error(666)
  enddo
  do f = 1,lattice_maxNslipFamily
    if (constitutive_nonlocal_Nslip(f,i) > 0_pInt) then
      if (constitutive_nonlocal_rhoSglEdgePos0(f,i) < 0.0_pReal)                            call IO_error(220)
      if (constitutive_nonlocal_rhoSglEdgeNeg0(f,i) < 0.0_pReal)                            call IO_error(220)
      if (constitutive_nonlocal_rhoSglScrewPos0(f,i) < 0.0_pReal)                           call IO_error(220)
      if (constitutive_nonlocal_rhoSglScrewNeg0(f,i) < 0.0_pReal)                           call IO_error(220)
      if (constitutive_nonlocal_rhoDipEdge0(f,i) < 0.0_pReal)                               call IO_error(220)
      if (constitutive_nonlocal_rhoDipScrew0(f,i) < 0.0_pReal)                              call IO_error(220)
      if (constitutive_nonlocal_burgersPerSlipFamily(f,i) <= 0.0_pReal)                     call IO_error(221)
      if (constitutive_nonlocal_v0PerSlipFamily(f,i) <= 0.0_pReal)                          call IO_error(226)
      if (constitutive_nonlocal_lambda0PerSlipFamily(f,i) <= 0.0_pReal)                     call IO_error(227)
      if (constitutive_nonlocal_dLowerEdgePerSlipFamily(f,i) <= 0.0_pReal)                  call IO_error(228)
      if (constitutive_nonlocal_dLowerScrewPerSlipFamily(f,i) <= 0.0_pReal)                 call IO_error(228)
    endif
  enddo
  if (any(constitutive_nonlocal_interactionSlipSlip(1:maxval(lattice_interactionSlipSlip(:,:,myStructure)),i) < 0.0_pReal)) &
                                                                                            call IO_error(229)
  if (constitutive_nonlocal_Q0(i) <= 0.0_pReal)                                             call IO_error(-1)
  if (constitutive_nonlocal_R(i) <= 0.0_pReal)                                              call IO_error(-1)
  if (constitutive_nonlocal_atomicVolume(i) <= 0.0_pReal)                                   call IO_error(230)
  if (constitutive_nonlocal_D0(i) <= 0.0_pReal)                                             call IO_error(231)
  if (constitutive_nonlocal_Qsd(i) <= 0.0_pReal)                                            call IO_error(232)
  if (constitutive_nonlocal_aTolRho(i) <= 0.0_pReal)                                        call IO_error(233)
    
  
!*** determine total number of active slip systems
  
  constitutive_nonlocal_Nslip(:,i) = min( lattice_NslipSystem(:, myStructure), constitutive_nonlocal_Nslip(:,i) )                   ! we can't use more slip systems per family than specified in lattice 
  constitutive_nonlocal_totalNslip(i) = sum(constitutive_nonlocal_Nslip(:,i))

enddo


!*** allocation of variables whose size depends on the total number of active slip systems

maxTotalNslip = maxval(constitutive_nonlocal_totalNslip)

allocate(constitutive_nonlocal_burgersPerSlipSystem(maxTotalNslip, maxNinstance))
constitutive_nonlocal_burgersPerSlipSystem = 0.0_pReal

allocate(constitutive_nonlocal_v0PerSlipSystem(maxTotalNslip, maxNinstance))
constitutive_nonlocal_v0PerSlipSystem = 0.0_pReal

allocate(constitutive_nonlocal_lambda0PerSlipSystem(maxTotalNslip, maxNinstance))
constitutive_nonlocal_lambda0PerSlipSystem = 0.0_pReal

allocate(constitutive_nonlocal_dLowerEdgePerSlipSystem(maxTotalNslip, maxNinstance))
constitutive_nonlocal_dLowerEdgePerSlipSystem = 0.0_pReal

allocate(constitutive_nonlocal_dLowerScrewPerSlipSystem(maxTotalNslip, maxNinstance))
constitutive_nonlocal_dLowerScrewPerSlipSystem = 0.0_pReal

allocate(constitutive_nonlocal_forestProjectionEdge(maxTotalNslip, maxTotalNslip, maxNinstance))
constitutive_nonlocal_forestProjectionEdge = 0.0_pReal

allocate(constitutive_nonlocal_forestProjectionScrew(maxTotalNslip, maxTotalNslip, maxNinstance))
constitutive_nonlocal_forestProjectionScrew = 0.0_pReal

allocate(constitutive_nonlocal_interactionMatrixSlipSlip(maxTotalNslip, maxTotalNslip, maxNinstance))
constitutive_nonlocal_interactionMatrixSlipSlip = 0.0_pReal

allocate(constitutive_nonlocal_v(maxTotalNslip, 4, homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems))
constitutive_nonlocal_v = 0.0_pReal

allocate(constitutive_nonlocal_rhoDotFlux(maxTotalNslip, 8, homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems))
constitutive_nonlocal_rhoDotFlux = 0.0_pReal

allocate(constitutive_nonlocal_compatibility(2,maxTotalNslip, maxTotalNslip, FE_maxNipNeighbors, mesh_maxNips, mesh_NcpElems))
constitutive_nonlocal_compatibility = 0.0_pReal

do i = 1,maxNinstance
  
  myStructure = constitutive_nonlocal_structure(i)                                                                                  ! lattice structure of this instance
    
!*** Inverse lookup of my slip system family and the slip system in lattice
  
  l = 0_pInt
  do f = 1,lattice_maxNslipFamily
    do s = 1,constitutive_nonlocal_Nslip(f,i)
      l = l + 1
      constitutive_nonlocal_slipFamily(l,i) = f
      constitutive_nonlocal_slipSystemLattice(l,i) = sum(lattice_NslipSystem(1:f-1, myStructure)) + s
  enddo; enddo
  
  
!*** determine size of state array
  
  ns = constitutive_nonlocal_totalNslip(i)
  constitutive_nonlocal_sizeState(i) = size(constitutive_nonlocal_listBasicStates) * ns & 
                                       + ( size(constitutive_nonlocal_listDependentStates) - 1_pInt ) * ns + 6_pInt
  constitutive_nonlocal_sizeDotState(i) = size(constitutive_nonlocal_listBasicStates) * ns
  
  
!*** determine size of postResults array
  
  do o = 1,maxval(phase_Noutput)
    select case(constitutive_nonlocal_output(o,i))
      case( 'rho', &
            'delta', &
            'rho_edge', &
            'rho_screw', &
            'rho_sgl', &
            'delta_sgl', &
            'rho_sgl_edge', &
            'rho_sgl_edge_pos', &
            'rho_sgl_edge_neg', &
            'rho_sgl_screw', &
            'rho_sgl_screw_pos', &
            'rho_sgl_screw_neg', &
            'rho_sgl_mobile', &
            'rho_sgl_edge_mobile', &
            'rho_sgl_edge_pos_mobile', &
            'rho_sgl_edge_neg_mobile', &
            'rho_sgl_screw_mobile', &
            'rho_sgl_screw_pos_mobile', &
            'rho_sgl_screw_neg_mobile', &
            'rho_sgl_immobile', &
            'rho_sgl_edge_immobile', &
            'rho_sgl_edge_pos_immobile', &
            'rho_sgl_edge_neg_immobile', &
            'rho_sgl_screw_immobile', &
            'rho_sgl_screw_pos_immobile', &
            'rho_sgl_screw_neg_immobile', &
            'rho_dip', &
            'delta_dip', &
            'rho_dip_edge', &
            'rho_dip_screw', &
            'excess_rho', &
            'excess_rho_edge', &
            'excess_rho_screw', &
            'rho_forest', &
            'shearrate', &
            'resolvedstress', &
            'resolvedstress_internal', &
            'resolvedstress_external', &
            'resistance', &
            'rho_dot', &
            'rho_dot_sgl', &
            'rho_dot_dip', &
            'rho_dot_gen', &
            'rho_dot_gen_edge', &
            'rho_dot_gen_screw', &
            'rho_dot_sgl2dip', &
            'rho_dot_dip2sgl', &
            'rho_dot_ann_ath', &
            'rho_dot_ann_the', &
            'rho_dot_flux', &
            'rho_dot_flux_edge', &
            'rho_dot_flux_screw', &
            'dislocationvelocity', &
            'fluxdensity_edge_pos_x', &
            'fluxdensity_edge_pos_y', &
            'fluxdensity_edge_pos_z', &
            'fluxdensity_edge_neg_x', &
            'fluxdensity_edge_neg_y', &
            'fluxdensity_edge_neg_z', &
            'fluxdensity_screw_pos_x', &
            'fluxdensity_screw_pos_y', &
            'fluxdensity_screw_pos_z', &
            'fluxdensity_screw_neg_x', &
            'fluxdensity_screw_neg_y', &
            'fluxdensity_screw_neg_z', &
            'd_upper_edge', &
            'd_upper_screw', &
            'd_upper_dot_edge', &
            'd_upper_dot_screw' )
        mySize = constitutive_nonlocal_totalNslip(i)
      case default
        mySize = 0_pInt
    end select

    if (mySize > 0_pInt) then                                                                                                       ! any meaningful output found                               
      constitutive_nonlocal_sizePostResult(o,i) = mySize
      constitutive_nonlocal_sizePostResults(i)  = constitutive_nonlocal_sizePostResults(i) + mySize
    endif
  enddo
  
  
!*** elasticity matrix and shear modulus according to material.config
  
  select case (myStructure)
    case(1:2)                                                                                                                       ! cubic(s)
      forall(k=1:3)
        forall(j=1:3) constitutive_nonlocal_Cslip_66(k,j,i) = constitutive_nonlocal_C12(i)
        constitutive_nonlocal_Cslip_66(k,k,i) = constitutive_nonlocal_C11(i)
        constitutive_nonlocal_Cslip_66(k+3,k+3,i) = constitutive_nonlocal_C44(i)
      end forall
    case(3:)                                                                                                                        ! all hex
      constitutive_nonlocal_Cslip_66(1,1,i) = constitutive_nonlocal_C11(i)
      constitutive_nonlocal_Cslip_66(2,2,i) = constitutive_nonlocal_C11(i)
      constitutive_nonlocal_Cslip_66(3,3,i) = constitutive_nonlocal_C33(i)
      constitutive_nonlocal_Cslip_66(1,2,i) = constitutive_nonlocal_C12(i)
      constitutive_nonlocal_Cslip_66(2,1,i) = constitutive_nonlocal_C12(i)
      constitutive_nonlocal_Cslip_66(1,3,i) = constitutive_nonlocal_C13(i)
      constitutive_nonlocal_Cslip_66(3,1,i) = constitutive_nonlocal_C13(i)
      constitutive_nonlocal_Cslip_66(2,3,i) = constitutive_nonlocal_C13(i)
      constitutive_nonlocal_Cslip_66(3,2,i) = constitutive_nonlocal_C13(i)
      constitutive_nonlocal_Cslip_66(4,4,i) = constitutive_nonlocal_C44(i)
      constitutive_nonlocal_Cslip_66(5,5,i) = constitutive_nonlocal_C44(i)
      constitutive_nonlocal_Cslip_66(6,6,i) = 0.5_pReal*(constitutive_nonlocal_C11(i)- constitutive_nonlocal_C12(i))
  end select
  constitutive_nonlocal_Cslip_66(:,:,i) = math_Mandel3333to66(math_Voigt66to3333(constitutive_nonlocal_Cslip_66(:,:,i)))
  constitutive_nonlocal_Cslip_3333(:,:,:,:,i) = math_Voigt66to3333(constitutive_nonlocal_Cslip_66(:,:,i))

  constitutive_nonlocal_Gmod(i) = 0.2_pReal * ( constitutive_nonlocal_C11(i) - constitutive_nonlocal_C12(i) &
                                                + 3.0_pReal*constitutive_nonlocal_C44(i) )                                          ! (C11iso-C12iso)/2 with C11iso=(3*C11+2*C12+4*C44)/5 and C12iso=(C11+4*C12-2*C44)/5
  constitutive_nonlocal_nu(i) =   ( constitutive_nonlocal_C11(i) + 4.0_pReal*constitutive_nonlocal_C12(i) &
                                    - 2.0_pReal*constitutive_nonlocal_C44(i) ) &
                                / ( 4.0_pReal*constitutive_nonlocal_C11(i) + 6.0_pReal*constitutive_nonlocal_C12(i) &
                                    + 2.0_pReal*constitutive_nonlocal_C44(i) )                                                      ! C12iso/(C11iso+C12iso) with C11iso=(3*C11+2*C12+4*C44)/5 and C12iso=(C11+4*C12-2*C44)/5
  
  
  do s1 = 1,ns
    
    f = constitutive_nonlocal_slipFamily(s1,i)
    
!*** burgers vector, dislocation velocity prefactor, mean free path prefactor and minimum dipole distance for each slip system
  
    constitutive_nonlocal_burgersPerSlipSystem(s1,i) = constitutive_nonlocal_burgersPerSlipFamily(f,i)
    constitutive_nonlocal_v0PerSlipSystem(s1,i) = constitutive_nonlocal_v0PerSlipFamily(f,i)
    constitutive_nonlocal_lambda0PerSlipSystem(s1,i) = constitutive_nonlocal_lambda0PerSlipFamily(f,i)
    constitutive_nonlocal_dLowerEdgePerSlipSystem(s1,i) = constitutive_nonlocal_dLowerEdgePerSlipFamily(f,i)
    constitutive_nonlocal_dLowerScrewPerSlipSystem(s1,i) = constitutive_nonlocal_dLowerScrewPerSlipFamily(f,i)

    do s2 = 1,ns
      
!*** calculation of forest projections for edge and screw dislocations. s2 acts as forest to s1

      constitutive_nonlocal_forestProjectionEdge(s1, s2, i) &
          = abs(math_mul3x3(lattice_sn(:, constitutive_nonlocal_slipSystemLattice(s1,i), myStructure), &
                            lattice_st(:, constitutive_nonlocal_slipSystemLattice(s2,i), myStructure)))                             ! forest projection of edge dislocations is the projection of (t = b x n) onto the slip normal of the respective slip plane
      
      constitutive_nonlocal_forestProjectionScrew(s1, s2, i) &
          = abs(math_mul3x3(lattice_sn(:, constitutive_nonlocal_slipSystemLattice(s1,i), myStructure), &
                            lattice_sd(:, constitutive_nonlocal_slipSystemLattice(s2,i), myStructure)))                             ! forest projection of screw dislocations is the projection of b onto the slip normal of the respective splip plane
  
!*** calculation of interaction matrices

      constitutive_nonlocal_interactionMatrixSlipSlip(s1, s2, i) &
          = constitutive_nonlocal_interactionSlipSlip( lattice_interactionSlipSlip(constitutive_nonlocal_slipSystemLattice(s1,i), &
                                                                                   constitutive_nonlocal_slipSystemLattice(s2,i), &
                                                                                   myStructure), &
                                                       i )
  
  enddo; enddo

enddo

endsubroutine



!*********************************************************************
!* initial microstructural state (just the "basic" states)           *
!*********************************************************************
pure function constitutive_nonlocal_stateInit(myInstance)

use prec,     only: pReal, &
                    pInt
use lattice,  only: lattice_maxNslipFamily
implicit none

!*** input variables
integer(pInt), intent(in) ::  myInstance                      ! number specifying the current instance of the constitution

!*** output variables
real(pReal), dimension(constitutive_nonlocal_sizeState(myInstance)) :: &
                              constitutive_nonlocal_stateInit

!*** local variables
real(pReal), dimension(constitutive_nonlocal_totalNslip(myInstance)) :: &              
                              rhoSglEdgePos, &                ! positive edge dislocation density
                              rhoSglEdgeNeg, &                ! negative edge dislocation density
                              rhoSglScrewPos, &               ! positive screw dislocation density
                              rhoSglScrewNeg, &               ! negative screw dislocation density
                              rhoSglEdgePosUsed, &            ! used positive edge dislocation density
                              rhoSglEdgeNegUsed, &            ! used negative edge dislocation density
                              rhoSglScrewPosUsed, &           ! used positive screw dislocation density
                              rhoSglScrewNegUsed, &           ! used negative screw dislocation density
                              rhoDipEdge, &                   ! edge dipole dislocation density
                              rhoDipScrew, &                  ! screw dipole dislocation density
                              rhoForest, &                    ! forest dislocation density
                              tauSlipThreshold                ! threshold shear stress for slip
integer(pInt)                 ns, &                           ! short notation for total number of active slip systems 
                              f, &                            ! index of lattice family
                              from, &
                              upto, &
                              s                               ! index of slip system


constitutive_nonlocal_stateInit = 0.0_pReal
ns = constitutive_nonlocal_totalNslip(myInstance)

!*** set the basic state variables

do f = 1,lattice_maxNslipFamily
  from = 1+sum(constitutive_nonlocal_Nslip(1:f-1,myInstance))
  upto = sum(constitutive_nonlocal_Nslip(1:f,myInstance))
  rhoSglEdgePos(from:upto)  = constitutive_nonlocal_rhoSglEdgePos0(f, myInstance)
  rhoSglEdgeNeg(from:upto)  = constitutive_nonlocal_rhoSglEdgeNeg0(f, myInstance)
  rhoSglScrewPos(from:upto) = constitutive_nonlocal_rhoSglScrewPos0(f, myInstance)
  rhoSglScrewNeg(from:upto) = constitutive_nonlocal_rhoSglScrewNeg0(f, myInstance)
  rhoSglEdgePosUsed(from:upto)  = 0.0_pReal
  rhoSglEdgeNegUsed(from:upto)  = 0.0_pReal
  rhoSglScrewPosUsed(from:upto) = 0.0_pReal
  rhoSglScrewNegUsed(from:upto) = 0.0_pReal
  rhoDipEdge(from:upto)  = constitutive_nonlocal_rhoDipEdge0(f, myInstance)
  rhoDipScrew(from:upto) = constitutive_nonlocal_rhoDipScrew0(f, myInstance)
enddo


!*** put everything together and in right order

constitutive_nonlocal_stateInit(      1:   ns) = rhoSglEdgePos
constitutive_nonlocal_stateInit(   ns+1: 2*ns) = rhoSglEdgeNeg
constitutive_nonlocal_stateInit( 2*ns+1: 3*ns) = rhoSglScrewPos
constitutive_nonlocal_stateInit( 3*ns+1: 4*ns) = rhoSglScrewNeg
constitutive_nonlocal_stateInit( 4*ns+1: 5*ns) = rhoSglEdgePosUsed
constitutive_nonlocal_stateInit( 5*ns+1: 6*ns) = rhoSglEdgeNegUsed
constitutive_nonlocal_stateInit( 6*ns+1: 7*ns) = rhoSglScrewPosUsed
constitutive_nonlocal_stateInit( 7*ns+1: 8*ns) = rhoSglScrewNegUsed
constitutive_nonlocal_stateInit( 8*ns+1: 9*ns) = rhoDipEdge
constitutive_nonlocal_stateInit( 9*ns+1:10*ns) = rhoDipScrew

endfunction



!*********************************************************************
!* absolute state tolerance                                          *
!*********************************************************************
pure function constitutive_nonlocal_aTolState(myInstance)

use prec,     only: pReal, &
                    pInt
implicit none

!*** input variables
integer(pInt), intent(in) ::  myInstance                      ! number specifying the current instance of the constitution

!*** output variables
real(pReal), dimension(constitutive_nonlocal_sizeState(myInstance)) :: &
                              constitutive_nonlocal_aTolState ! absolute state tolerance for the current instance of this constitution

!*** local variables

constitutive_nonlocal_aTolState = constitutive_nonlocal_aTolRho(myInstance)

endfunction



!*********************************************************************
!* calculates homogenized elacticity matrix                          *
!*********************************************************************
pure function constitutive_nonlocal_homogenizedC(state,g,ip,el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_constitutionInstance
implicit none

!*** input variables
integer(pInt), intent(in) ::    g, &                                ! current grain ID
                                ip, &                               ! current integration point
                                el                                  ! current element
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state ! microstructural state

!*** output variables
real(pReal), dimension(6,6) ::  constitutive_nonlocal_homogenizedC  ! homogenized elasticity matrix

!*** local variables
integer(pInt)                   myInstance                          ! current instance of this constitution

myInstance = phase_constitutionInstance(material_phase(g,ip,el))

constitutive_nonlocal_homogenizedC = constitutive_nonlocal_Cslip_66(:,:,myInstance)
 
endfunction



!*********************************************************************
!* calculates quantities characterizing the microstructure           *
!*********************************************************************
subroutine constitutive_nonlocal_microstructure(state, Temperature, Tstar_v, Fe, Fp, g, ip, el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use math,     only: math_Plain3333to99, &
                    math_Mandel33to6, &
                    math_Mandel6to33, &
                    math_mul33x33, &
                    math_mul3x3, &
                    math_mul33x3, &
                    math_inv3x3, &
                    math_det3x3, &
                    pi
use debug,    only: debugger, &
                    verboseDebugger
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_maxNipNeighbors, &
                    mesh_element, &
                    FE_NipNeighbors, &
                    mesh_ipNeighborhood, &
                    mesh_ipVolume, &
                    mesh_ipCenterOfGravity
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_localConstitution, &
                    phase_constitutionInstance
use lattice,  only: lattice_Sslip, &
                    lattice_Sslip_v, &
                    lattice_maxNslipFamily, &
                    lattice_NslipSystem, &
                    lattice_maxNslip, &
                    lattice_sd, &
                    lattice_sn, &
                    lattice_st

implicit none

!*** input variables
integer(pInt), intent(in) ::      g, &                        ! current grain ID
                                  ip, &                       ! current integration point
                                  el                          ! current element
real(pReal), intent(in) ::        Temperature                 ! temperature
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                  Fe, &                       ! elastic deformation gradient
                                  Fp                          ! plastic deformation gradient
real(pReal), dimension(6), intent(in) :: &
                                  Tstar_v                     ! 2nd Piola-Kirchhoff stress in Mandel notation


!*** input/output variables
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(inout) :: &
                                  state                       ! microstructural state

!*** output variables

!*** local variables
integer(pInt)                     myInstance, &               ! current instance of this constitution
                                  myStructure, &              ! current lattice structure
                                  ns, &                       ! short notation for the total number of active slip systems
                                  neighboring_el, &           ! element number of my neighbor
                                  neighboring_ip, &           ! integration point of my neighbor
                                  c, &                        ! index of dilsocation character (edge, screw)
                                  n, &                        ! index of my current neighbor
                                  s, &                        ! index of my current slip system
                                  t, &                        ! index of dilsocation type (e+, e-, s+, s-, used e+, used e-, used s+, used s-)
                                  sLattice, &                 ! index of my current slip system according to lattice order
                                  i, &
                                  j
real(pReal)                       nu                          ! poisson's ratio
real(pReal), dimension(3,2) ::    rhoExcessDifference, &      ! finite differences of excess density (in 3 directions for edge and screw)
                                  disloGradients              ! spatial gradient in excess dislocation density (in 3 directions for edge and screw)
real(pReal), dimension(3,3) ::    sigma, &                    ! dislocation stress for one slip system in its slip system frame
                                  lattice2slip, &             ! orthogonal transformation matrix from lattice coordinate system to slip coordinate system with e1=bxn, e2=b, e3=n (passive rotation!!!)
                                  F, &                        ! total deformation gradient
                                  neighboring_F, &            ! total deformation gradient of neighbor
                                  invFe, &                    ! inverse elastic deformation gradient
                                  invPositionDifference       ! inverse of a 3x3 matrix containing finite differences of pairs of position vectors
real(pReal), dimension(6) ::      Tdislocation_v              ! dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress in Mandel notation
real(pReal), dimension(2,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                  rhoExcess                   ! central excess density
real(pReal), dimension(6,2,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                  neighboring_rhoExcess       ! excess density for each neighbor, dislo character and slip system
real(pReal), dimension(6,3) ::    neighboring_position        ! position vector of each neighbor when seen from the centreal material point's lattice frame
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),8) :: &
                                  rhoSgl                      ! single dislocation density (edge+, edge-, screw+, screw-, used edge+, used edge-, used screw+, used screw-)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                  rhoDip                      ! dipole dislocation density (edge, screw)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                  transmissivity, &           ! transmissivity
                                  rhoForest, &                ! forest dislocation density
                                  tauThreshold, &             ! threshold shear stress
                                  tau                         ! resolved shear stress

myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance)
ns = constitutive_nonlocal_totalNslip(myInstance)


!**********************************************************************
!*** set fluxes to zero

constitutive_nonlocal_rhoDotFlux(:,:,g,ip,el) = 0.0_pReal


!**********************************************************************
!*** get basic states

forall (t = 1:4) rhoSgl(:,t) = max(state(g,ip,el)%p((t-1)*ns+1:t*ns), 0.0_pReal)                          ! ensure positive single mobile densities
forall (t = 5:8) rhoSgl(:,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (c = 1:2) rhoDip(:,c) = max(state(g,ip,el)%p((c+7)*ns+1:(c+8)*ns), 0.0_pReal)                      ! ensure positive dipole densities
where(rhoSgl(:,1:4) < min(0.1, 0.01*constitutive_nonlocal_aTolRho(myInstance))) rhoSgl(:,1:4) = 0.0_pReal ! delete non-significant single density


!**********************************************************************
!*** calculate dependent states

!*** calculate the forest dislocation density

forall (s = 1:ns) &
  rhoForest(s) =  dot_product( ( sum(abs(rhoSgl(:,(/1,2,5,6/))),2) + rhoDip(:,1) ), &
                               constitutive_nonlocal_forestProjectionEdge(s, 1:ns, myInstance) ) & 
                + dot_product( ( sum(abs(rhoSgl(:,(/3,4,7,8/))),2) + rhoDip(:,2) ), &
                               constitutive_nonlocal_forestProjectionScrew(s, 1:ns, myInstance) )         ! calculation of forest dislocation density as projection of screw and edge dislocations
! if (debugger) write(6,'(a30,3(i3,x),/,12(e10.3,x),/)') 'forest dislocation density at ',g,ip,el, rhoForest


!*** calculate the threshold shear stress for dislocation slip 

forall (s = 1:ns) &
  tauThreshold(s) =   constitutive_nonlocal_Gmod(myInstance) & 
                    * constitutive_nonlocal_burgersPerSlipSystem(s, myInstance) &
                    * sqrt( dot_product( ( sum(abs(rhoSgl),2) + sum(abs(rhoDip),2) ), &
                                         constitutive_nonlocal_interactionMatrixSlipSlip(s, 1:ns, myInstance) ) )
! if (debugger) write(6,'(a22,3(i3,x),/,12(f10.5,x),/)') 'tauThreshold / MPa at ',g,ip,el, tauThreshold/1e6


!*** calculate the dislocation stress of the neighboring excess dislocation densities

Tdislocation_v = 0.0_pReal
F = math_mul33x33(Fe(:,:,g,ip,el), Fp(:,:,g,ip,el))
invFe = math_inv3x3(Fe(:,:,g,ip,el))
nu = constitutive_nonlocal_nu(myInstance)

forall (s = 1:ns, c = 1:2) &
  rhoExcess(c,s) =  state(g,ip,el)%p((2*c-2)*ns+s) + abs(state(g,ip,el)%p((2*c+2)*ns+s)) &
                  - state(g,ip,el)%p((2*c-1)*ns+s) - abs(state(g,ip,el)%p((2*c+3)*ns+s))
do n = 1,6
  neighboring_el = mesh_ipNeighborhood(1,n,ip,el)
  neighboring_ip = mesh_ipNeighborhood(2,n,ip,el)  
  if ( neighboring_ip == 0 .or. neighboring_el == 0 ) then                                                ! at free surfaces ...
    neighboring_el = el                                                                                   ! ... use central values instead of neighboring values
    neighboring_ip = ip
    neighboring_position(n,:) = 0.0_pReal
    neighboring_rhoExcess(n,:,:) = rhoExcess
  elseif (phase_localConstitution(material_phase(1,neighboring_ip,neighboring_el))) then                  ! for neighbors with local constitution
    neighboring_el = el                                                                                   ! ... use central values instead of neighboring values
    neighboring_ip = ip
    neighboring_position(n,:) = 0.0_pReal
    neighboring_rhoExcess(n,:,:) = rhoExcess
  elseif (myStructure /= &
          constitutive_nonlocal_structure(phase_constitutionInstance(material_phase(1,neighboring_ip,neighboring_el)))) then ! for neighbors with different crystal structure
    neighboring_el = el                                                                                   ! ... use central values instead of neighboring values
    neighboring_ip = ip
    neighboring_position(n,:) = 0.0_pReal
    neighboring_rhoExcess(n,:,:) = rhoExcess
  else
    forall (s = 1:ns, c = 1:2) &
      neighboring_rhoExcess(n,c,s) =      state(g,neighboring_ip,neighboring_el)%p((2*c-2)*ns+s) &
                                    + abs(state(g,neighboring_ip,neighboring_el)%p((2*c+2)*ns+s)) &
                                    -     state(g,neighboring_ip,neighboring_el)%p((2*c-1)*ns+s) &
                                    - abs(state(g,neighboring_ip,neighboring_el)%p((2*c+3)*ns+s))
    transmissivity = sum(constitutive_nonlocal_compatibility(2,:,:,n,ip,el)**2.0_pReal, 1)
    if ( any(transmissivity < 0.99_pReal) ) then                                                          ! at grain boundary (=significantly decreased transmissivity) ...
      neighboring_el = el                                                                                 ! ... use central values instead of neighboring values
      neighboring_ip = ip
      neighboring_position(n,:) = 0.0_pReal
      neighboring_rhoExcess(n,:,:) = rhoExcess
    else 
      neighboring_F = math_mul33x33(Fe(:,:,g,neighboring_ip,neighboring_el), Fp(:,:,g,neighboring_ip,neighboring_el))
      neighboring_position(n,:) = &
            0.5_pReal * math_mul33x3( math_mul33x33(invFe,neighboring_F) + Fp(:,:,g,ip,el), &
                                      mesh_ipCenterOfGravity(:,neighboring_ip,neighboring_el) - mesh_ipCenterOfGravity(:,ip,el) )
    endif
  endif  
enddo

invPositionDifference = math_inv3x3(neighboring_position((/1,3,5/),:) - neighboring_position((/2,4,6/),:))

do s = 1,ns

  lattice2slip = transpose( reshape( (/ lattice_st(:, constitutive_nonlocal_slipSystemLattice(s,myInstance), myStructure), &
                                        lattice_sd(:, constitutive_nonlocal_slipSystemLattice(s,myInstance), myStructure), &
                                        lattice_sn(:, constitutive_nonlocal_slipSystemLattice(s,myInstance), myStructure) /), &
                                     (/ 3,3 /) ) )

  rhoExcessDifference = neighboring_rhoExcess((/1,3,5/),:,s) - neighboring_rhoExcess((/2,4,6/),:,s)
  forall (c = 1:2) &    
    disloGradients(:,c) = math_mul33x3( lattice2slip, math_mul33x3(invPositionDifference, rhoExcessDifference(:,c)) )
  
  sigma = 0.0_pReal
  sigma(1,1) = + (-0.06066_pReal + nu*0.41421_pReal) / (1.0_pReal-nu) * disloGradients(3,1)
  sigma(2,2) = + 0.32583_pReal / (1.0_pReal-nu) * disloGradients(3,1)
  sigma(3,3) = + 0.14905_pReal / (1.0_pReal-nu) * disloGradients(3,1)
  sigma(1,2) = + 0.20711_pReal * disloGradients(3,2)
  sigma(2,3) = - 0.08839_pReal / (1.0_pReal-nu) * disloGradients(2,1) - 0.20711_pReal * disloGradients(1,2)
  sigma(2,1) = sigma(1,2)
  sigma(3,2) = sigma(2,3)
  
  forall (i=1:3, j=1:3) &
    sigma(i,j) = sigma(i,j) * constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgersPerSlipSystem(s,myInstance) &
                            * constitutive_nonlocal_R(myInstance)**2.0_pReal

  Tdislocation_v = Tdislocation_v + math_Mandel33to6( math_mul33x33(transpose(lattice2slip), math_mul33x33(sigma, lattice2slip) ) )
    
enddo

!**********************************************************************
!*** set states

state(g,ip,el)%p(1:8*ns) = reshape(rhoSgl,(/8*ns/))                                                       ! ensure positive single mobile densities
state(g,ip,el)%p(8*ns+1:10*ns) = reshape(rhoDip,(/2*ns/))                                                 ! ensure positive dipole densities
state(g,ip,el)%p(10*ns+1:11*ns) = rhoForest
state(g,ip,el)%p(11*ns+1:12*ns) = tauThreshold
state(g,ip,el)%p(12*ns+1:12*ns+6) = Tdislocation_v

endsubroutine



!*********************************************************************
!* calculates kinetics                                               *
!*********************************************************************
subroutine constitutive_nonlocal_kinetics(Tstar_v, Temperature, state, g, ip, el, dv_dtau)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use math,     only: math_mul6x6, &
                    math_Mandel6to33
use debug,    only: debugger, &
                    verboseDebugger, &
                    debug_g, &
                    debug_i, &
                    debug_e
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_constitutionInstance
use lattice,  only: lattice_Sslip, &
                    lattice_Sslip_v

implicit none

!*** input variables
integer(pInt), intent(in) ::                g, &                        ! current grain number
                                            ip, &                       ! current integration point
                                            el                          ! current element number
real(pReal), intent(in) ::                  Temperature                 ! temperature
type(p_vec), intent(in) ::                  state                       ! microstructural state
real(pReal), dimension(6), intent(in) ::    Tstar_v                     ! 2nd Piola-Kirchhoff stress in Mandel notation

!*** output variables
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))), &
                   intent(out), optional :: dv_dtau                     ! velocity derivative with respect to resolved shear stress

!*** local variables
integer(pInt)                               myInstance, &               ! current instance of this constitution
                                            myStructure, &              ! current lattice structure
                                            ns, &                       ! short notation for the total number of active slip systems
                                            t, &                        ! dislocation type
                                            s                           ! index of my current slip system
real(pReal), dimension(6) ::                Tdislocation_v              ! dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            tauThreshold, &             ! threshold shear stress
                                            tau, &                      ! resolved shear stress
                                            rhoForest                   ! forest dislocation density
real(pReal)                                 boltzmannProbability


myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)

rhoForest = state%p(10*ns+1:11*ns)
tauThreshold = state%p(11*ns+1:12*ns)
Tdislocation_v = state%p(12*ns+1:12*ns+6)

tau = 0.0_pReal
constitutive_nonlocal_v(:,:,g,ip,el) = 0.0_pReal
if (present(dv_dtau)) dv_dtau = 0.0_pReal

if ( Temperature > 0.0_pReal ) then
  do s = 1,ns
  
    tau(s) = math_mul6x6( Tstar_v + Tdislocation_v, &
                          lattice_Sslip_v(:,constitutive_nonlocal_slipSystemLattice(s,myInstance),myStructure) )
      
    if ( abs(tau(s)) > 0.0_pReal ) then
      boltzmannProbability = dexp( -constitutive_nonlocal_Q0(myInstance) * dsqrt(rhoForest(s)) / ( abs(tau(s)) * kB * Temperature) )
      constitutive_nonlocal_v(s,:,g,ip,el) = sign(constitutive_nonlocal_v0PerSlipSystem(s,myInstance), tau(s)) &
          / ( boltzmannProbability + constitutive_nonlocal_burgersPerSlipSystem(s,myInstance) * dsqrt(rhoForest(s)) ) &
          * boltzmannProbability

      
      if (present(dv_dtau)) &
        dv_dtau(s) = abs(constitutive_nonlocal_v(s,1,g,ip,el)) * constitutive_nonlocal_Q0(myInstance) &
          * constitutive_nonlocal_burgersPerSlipSystem(s,myInstance) * rhoForest(s) / ( tau(s)**2.0_pReal * kB * Temperature ) &
          / ( boltzmannProbability + constitutive_nonlocal_burgersPerSlipSystem(s,myInstance) * dsqrt(rhoForest(s)) )
    endif
  enddo
endif

!if (verboseDebugger .and. s) then 
!  !$OMP CRITICAL (write2out)
!    write(6,*) '::: kinetics',g,ip,el
!    write(6,*)
!    write(6,'(a,/,3(3(f12.3,x)/))') 'Tdislocation / MPa', math_Mandel6to33(Tdislocation_v/1e6)
!    write(6,'(a,/,3(3(f12.3,x)/))') 'Tstar / MPa', math_Mandel6to33(Tstar_v/1e6)
!    write(6,'(a,/,12(f12.5,x),/)') 'tau / MPa', tau/1e6_pReal
!    write(6,'(a,/,12(e12.5,x),/)') 'rhoForest / 1/m**2', rhoForest
!    write(6,'(a,/,4(12(f12.5,x),/))') 'v / 1e-3m/s', constitutive_nonlocal_v(:,:,g,ip,el)*1e3
!  !$OMP END CRITICAL (write2out)
!endif

endsubroutine



!*********************************************************************
!* calculates plastic velocity gradient and its tangent              *
!*********************************************************************
subroutine constitutive_nonlocal_LpAndItsTangent(Lp, dLp_dTstar99, Tstar_v, Temperature, state, g, ip, el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use math,     only: math_Plain3333to99, &
                    math_mul6x6, &
                    math_Mandel6to33
use debug,    only: debugger, &
                    verboseDebugger, &
                    debug_g, &
                    debug_i, &
                    debug_e
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_constitutionInstance
use lattice,  only: lattice_Sslip, &
                    lattice_Sslip_v

implicit none

!*** input variables
integer(pInt), intent(in) ::                g, &                        ! current grain number
                                            ip, &                       ! current integration point
                                            el                          ! current element number
real(pReal), intent(in) ::                  Temperature                 ! temperature
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            state                       ! microstructural state
real(pReal), dimension(6), intent(in) ::    Tstar_v                     ! 2nd Piola-Kirchhoff stress in Mandel notation

!*** output variables
real(pReal), dimension(3,3), intent(out) :: Lp                          ! plastic velocity gradient
real(pReal), dimension(9,9), intent(out) :: dLp_dTstar99                ! derivative of Lp with respect to Tstar (9x9 matrix)

!*** local variables
integer(pInt)                               myInstance, &               ! current instance of this constitution
                                            myStructure, &              ! current lattice structure
                                            ns, &                       ! short notation for the total number of active slip systems
                                            i, &
                                            j, &
                                            k, &
                                            l, &
                                            t, &                        ! dislocation type
                                            s, &                        ! index of my current slip system
                                            sLattice                    ! index of my current slip system according to lattice order
real(pReal), dimension(3,3,3,3) ::          dLp_dTstar3333              ! derivative of Lp with respect to Tstar (3x3x3x3 matrix)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),8) :: &
                                            rhoSgl                      ! single dislocation densities (including used) 
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            gdot                        ! shear rate per dislocation type 
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            tauThreshold, &             ! threshold shear stress
                                            gdotTotal, &                ! shear rate
                                            dv_dtau, &                  ! velocity derivative with respect to the shear stress
                                            dgdotTotal_dtau, &          ! derivative of the shear rate with respect to the shear stress
                                            rhoForest                   ! forest dislocation density


!*** initialize local variables

gdot = 0.0_pReal
Lp = 0.0_pReal
dLp_dTstar3333 = 0.0_pReal

myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)

!*** shortcut to state variables 

forall (t = 1:8) &
  rhoSgl(:,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (s = 1:ns, t = 5:8, rhoSgl(s,t) * constitutive_nonlocal_v(s,t-4,g,ip,el) < 0.0_pReal) &                                      ! contribution of used rho for changing sign of v
  rhoSgl(s,t-4) = rhoSgl(s,t-4) + abs(rhoSgl(s,t))

rhoForest = state(g,ip,el)%p(10*ns+1:11*ns)
tauThreshold = state(g,ip,el)%p(11*ns+1:12*ns)

call constitutive_nonlocal_kinetics(Tstar_v, Temperature, state(g,ip,el), g, ip, el, dv_dtau)                                       ! update dislocation velocity

!*** Calculation of gdot and its tangent

forall (t = 1:4) &
  gdot(:,t) =  rhoSgl(:,t) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) * constitutive_nonlocal_v(:,t,g,ip,el)
gdotTotal = sum(gdot,2)

dgdotTotal_dtau = sum(rhoSgl,2) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) * dv_dtau

!*** Calculation of Lp and its tangent

do s = 1,ns
  sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
  
  Lp = Lp + gdotTotal(s) * lattice_Sslip(:,:,sLattice,myStructure)
  
  forall (i=1:3,j=1:3,k=1:3,l=1:3) &
    dLp_dTstar3333(i,j,k,l) = dLp_dTstar3333(i,j,k,l) + dgdotTotal_dtau(s) * lattice_Sslip(i,j, sLattice,myStructure) &
                                                                           * lattice_Sslip(k,l, sLattice,myStructure) 
enddo

dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)

!if (verboseDebugger .and. (debug_g==g .and. debug_i==i .and. debug_e==e)) then 
!  !$OMP CRITICAL (write2out)
!    write(6,*) '::: LpandItsTangent',g,ip,el
!    write(6,*)
!    write(6,'(a,/,12(f12.5,x),/)') 'v / 1e-3m/s', constitutive_nonlocal_v(:,:,g,ip,el)*1e3
!    write(6,'(a,/,12(f12.5,x),/)') 'gdot / 1e-3',gdot*1e3_pReal
!    write(6,'(a,/,12(f12.5,x),/)') 'gdot total / 1e-3',gdotTotal*1e3_pReal
!    write(6,'(a,/,3(3(f12.7,x)/))') 'Lp',Lp
!    ! call flush(6)
!  !$OMP END CRITICAL (write2out)
!endif

endsubroutine



!*********************************************************************
!* rate of change of microstructure                                  *
!*********************************************************************
subroutine constitutive_nonlocal_dotState(dotState, Tstar_v, previousTstar_v, Fe, Fp, Temperature, dt_previous, &
                                          state, previousState, aTolState, timestep, orientation, g,ip,el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use IO,       only: IO_error
use debug,    only: debugger, &
                    debug_g, &
                    debug_i, &
                    debug_e, &
                    verboseDebugger
use math,     only: math_norm3, &
                    math_mul6x6, &
                    math_mul3x3, &
                    math_mul33x3, &
                    math_mul33x33, &
                    math_inv3x3, &
                    math_det3x3, &
                    math_Mandel6to33, &
                    math_QuaternionDisorientation, &
                    math_qRot, &
                    pi, &
                    NaN
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_maxNipNeighbors, &
                    mesh_element, &
                    FE_NipNeighbors, &
                    mesh_ipNeighborhood, &
                    mesh_ipVolume, &
                    mesh_ipArea, &
                    mesh_ipAreaNormal, &
                    mesh_ipCenterOfGravity
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_constitutionInstance, &
                    phase_localConstitution
use lattice,  only: lattice_Sslip, &
                    lattice_Sslip_v, &
                    lattice_sd, &
                    lattice_sn, &
                    lattice_st, &
                    lattice_maxNslipFamily, &
                    lattice_NslipSystem 
use FEsolving, only:theInc, &
                    FEsolving_execElem, & 
                    FEsolving_execIP

implicit none

!*** input variables
integer(pInt), intent(in) ::                g, &                      ! current grain number
                                            ip, &                     ! current integration point
                                            el                        ! current element number
real(pReal), intent(in) ::                  Temperature, &            ! temperature
                                            timestep, &               ! substepped crystallite time increment
                                            dt_previous               ! time increment between previous and current state
real(pReal), dimension(6), intent(in) ::    Tstar_v, &                ! current 2nd Piola-Kirchhoff stress in Mandel notation
                                            previousTstar_v           ! previous 2nd Piola-Kirchhoff stress in Mandel notation
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            Fe, &                     ! elastic deformation gradient
                                            Fp                        ! plastic deformation gradient
real(pReal), dimension(4,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            orientation               ! crystal lattice orientation
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            state, &                  ! current microstructural state
                                            previousState, &          ! previous microstructural state
                                            aTolState                 ! absolute state tolerance

!*** input/output variables
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(inout) :: &
                                            dotState                  ! evolution of state variables / microstructure
 
!*** output variables
 
!*** local variables
integer(pInt)                               myInstance, &             ! current instance of this constitution
                                            myStructure, &            ! current lattice structure
                                            ns, &                     ! short notation for the total number of active slip systems
                                            c, &                      ! character of dislocation
                                            n, &                      ! index of my current neighbor
                                            neighboring_el, &         ! element number of my neighbor
                                            neighboring_ip, &         ! integration point of my neighbor
                                            opposite_n, &             ! index of my opposite neighbor
                                            opposite_ip, &            ! ip of my opposite neighbor
                                            opposite_el, &            ! element index of my opposite neighbor
                                            t, &                      ! type of dislocation
                                            topp, &                   ! type of dislocation with opposite sign to t
                                            s, &                      ! index of my current slip system
                                            sLattice, &               ! index of my current slip system according to lattice order
                                            i
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),10) :: &
                                            rhoDot, &                 ! density evolution
                                            rhoDotRemobilization, &   ! density evolution by remobilization
                                            rhoDotMultiplication, &   ! density evolution by multiplication
                                            rhoDotFlux, &             ! density evolution by flux
                                            neighboring_rhoDotFlux, & ! density evolution by flux at neighbor
                                            rhoDotSingle2DipoleGlide, & ! density evolution by dipole formation (by glide)
                                            rhoDotAthermalAnnihilation, & ! density evolution by athermal annihilation
                                            rhoDotThermalAnnihilation, & ! density evolution by thermal annihilation
                                            rhoDotDipole2SingleStressChange, & ! density evolution by dipole dissociation (by stress increase)
                                            rhoDotSingle2DipoleStressChange ! density evolution by dipole formation (by stress decrease)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),8) :: &
                                            rhoSgl, &                 ! current single dislocation densities (positive/negative screw and edge without dipoles)
                                            previousRhoSgl            ! previous single dislocation densities (positive/negative screw and edge without dipoles)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            fluxdensity, &            ! flux density at central material point
                                            gdot                      ! shear rates
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            rhoForest, &              ! forest dislocation density
                                            tauThreshold, &           ! threshold shear stress
                                            tau, &                    ! current resolved shear stress
                                            previousTau, &            ! previous resolved shear stress
                                            invLambda, &              ! inverse of mean free path for dislocations
                                            vClimb                    ! climb velocity of edge dipoles
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                            rhoDip, &                 ! current dipole dislocation densities (screw and edge dipoles)
                                            previousRhoDip, &         ! previous dipole dislocation densities (screw and edge dipoles)
                                            dLower, &                 ! minimum stable dipole distance for edges and screws
                                            dUpper, &                 ! current maximum stable dipole distance for edges and screws
                                            previousDUpper, &         ! previous maximum stable dipole distance for edges and screws
                                            dUpperDot                 ! rate of change of the maximum stable dipole distance for edges and screws
real(pReal), dimension(3,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            m                         ! direction of dislocation motion
real(pReal), dimension(3,3) ::              F, &                      ! total deformation gradient
                                            neighboring_F, &          ! total deformation gradient of my neighbor
                                            Favg                      ! average total deformation gradient of me and my neighbor
real(pReal), dimension(6) ::                Tdislocation_v, &         ! current dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress
                                            previousTdislocation_v    ! previous dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress
real(pReal), dimension(3) ::                surfaceNormal, &          ! surface normal in lattice configuration
                                            surfaceNormal_currentconf ! surface normal in current configuration
real(pReal)                                 area, &                   ! area of the current interface
                                            detFe, &                  ! determinant of elastic defornmation gradient
                                            transmissivity, &         ! overall transmissivity of dislocation flux to neighboring material point
                                            lineLength, &             ! dislocation line length leaving the current interface
                                            D, &                      ! self diffusion
                                            correction
logical, dimension(3) ::                    periodicSurfaceFlux       ! flag indicating periodic fluxes at surfaces when surface normal points mainly in x, y and z direction respectively (in reference configuration)

if (verboseDebugger .and. (debug_g==g .and. debug_i==ip .and. debug_e==el)) then 
  !$OMP CRITICAL (write2out)
    write(6,*) '::: constitutive_nonlocal_dotState at ',g,ip,el
    write(6,*)
  !$OMP END CRITICAL (write2out)
endif

select case(mesh_element(2,el))
  case (1,6,7,8,9)
    ! all fine
  case default
    call IO_error(-1,el,ip,g,'element type not supported for nonlocal constitution')
  

end select

myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)

tau = 0.0_pReal
previousTau = 0.0_pReal
gdot = 0.0_pReal
dLower = 0.0_pReal
dUpper = 0.0_pReal
previousDUpper = 0.0_pReal
dUpperDot = 0.0_pReal

!*** shortcut to state variables 

forall (t = 1:8) rhoSgl(:,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (t = 1:8) previousRhoSgl(:,t) = previousState(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (c = 1:2) rhoDip(:,c) = state(g,ip,el)%p((7+c)*ns+1:(8+c)*ns)
forall (c = 1:2) previousRhoDip(:,c) = previousState(g,ip,el)%p((7+c)*ns+1:(8+c)*ns)
rhoForest = state(g,ip,el)%p(10*ns+1:11*ns)
tauThreshold = state(g,ip,el)%p(11*ns+1:12*ns)
Tdislocation_v = state(g,ip,el)%p(12*ns+1:12*ns+6)
previousTdislocation_v = previousState(g,ip,el)%p(12*ns+1:12*ns+6)

!*** sanity check for timestep

if (timestep <= 0.0_pReal) then                                                                                                     ! if illegal timestep...
  dotState(g,ip,el)%p = 0.0_pReal                                                                                                   ! ...return without doing anything (-> zero dotState)
  return
endif


!****************************************************************************
!*** Calculate shear rate

call constitutive_nonlocal_kinetics(Tstar_v, Temperature, state(g,ip,el), g, ip, el)                                                ! get velocities

forall (t = 1:4) &
  gdot(:,t) = rhoSgl(:,t) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) * constitutive_nonlocal_v(:,t,g,ip,el)
forall (s = 1:ns, t = 1:4, rhoSgl(s,t+4) * constitutive_nonlocal_v(s,t,g,ip,el) < 0.0_pReal) &                                      ! contribution of used rho for changing sign of v
  gdot(s,t) = gdot(s,t) + abs(rhoSgl(s,t+4)) * constitutive_nonlocal_burgersPerSlipSystem(s,myInstance) &
                                             * constitutive_nonlocal_v(s,t,g,ip,el)

if (verboseDebugger .and. (debug_g==g .and. debug_i==ip .and. debug_e==el)) then 
  !$OMP CRITICAL (write2out)
    write(6,'(a,/,10(12(e12.5,x),/))') 'rho / 1/m^2', rhoSgl, rhoDip
    write(6,'(a,/,4(12(e12.5,x),/))') 'v / m/s', constitutive_nonlocal_v(:,:,g,ip,el)
    write(6,'(a,/,4(12(e12.5,x),/))') 'gdot / 1/s',gdot
  !$OMP END CRITICAL (write2out)
endif


!****************************************************************************
!*** calculate limits for stable dipole height and its rate of change

do s = 1,ns   ! loop over slip systems
  sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
  
  tau(s) = math_mul6x6( Tstar_v + Tdislocation_v, lattice_Sslip_v(:,sLattice,myStructure) )
  previousTau(s) = math_mul6x6( previousTstar_v + previousTdislocation_v, lattice_Sslip_v(:,sLattice,myStructure) )
  
enddo

dLower(:,1) = constitutive_nonlocal_dLowerEdgePerSlipSystem(:,myInstance)
dLower(:,2) = constitutive_nonlocal_dLowerScrewPerSlipSystem(:,myInstance)
dUpper(:,2) = min( 1.0_pReal / sqrt( sum(abs(rhoSgl),2)+sum(rhoDip,2) ), &
                   constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                                                          / ( 8.0_pReal * pi * abs(tau) ) )
dUpper(:,1) = dUpper(:,2) / ( 1.0_pReal - constitutive_nonlocal_nu(myInstance) )
previousDUpper(:,2) = min( 1.0_pReal / sqrt(  sum(abs(previousRhoSgl),2) + sum(previousRhoDip,2) ), &
                           constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                                                                  / ( 8.0_pReal * pi * abs(previousTau) ) )
previousDUpper(:,1) = previousDUpper(:,2) / ( 1.0_pReal - constitutive_nonlocal_nu(myInstance) )

if (dt_previous > 0.0_pReal) dUpperDot = (dUpper - previousDUpper) / dt_previous


!****************************************************************************
!*** dislocation remobilization (bauschinger effect)

rhoDotRemobilization = 0.0_pReal
if (timestep > 0.0_pReal) then
  do t = 1,4
    do s = 1,ns
      if (rhoSgl(s,t+4) * constitutive_nonlocal_v(s,t,g,ip,el) < 0.0_pReal) then
        rhoDotRemobilization(s,t) = abs(rhoSgl(s,t+4)) / timestep
        rhoSgl(s,t) = rhoSgl(s,t) + abs(rhoSgl(s,t+4))
        rhoDotRemobilization(s,t+4) = - rhoSgl(s,t+4) / timestep
        rhoSgl(s,t+4) = 0.0_pReal
      endif
    enddo
  enddo
endif


!****************************************************************************
!*** calculate dislocation multiplication

rhoDotMultiplication = 0.0_pReal
where (rhoSgl(:,3:4) > 0.0_pReal) &
  rhoDotMultiplication(:,1:2) = spread(0.5_pReal * sum(abs(gdot(:,3:4)),2) * sqrt(rhoForest)  &
                                                 / constitutive_nonlocal_lambda0PerSlipSystem(:,myInstance) &
                                                 / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance), 2, 2)
where (rhoSgl(:,1:2) > 0.0_pReal) &
  rhoDotMultiplication(:,3:4) = spread(0.5_pReal * sum(abs(gdot(:,1:2)),2) * sqrt(rhoForest)  &
                                                 / constitutive_nonlocal_lambda0PerSlipSystem(:,myInstance) &
                                                 / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance), 2, 2)


!****************************************************************************
!*** calculate dislocation fluxes

rhoDotFlux = 0.0_pReal
periodicSurfaceFlux = (/.false.,.true.,.false./)

m(:,:,1) =  lattice_sd(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)
m(:,:,2) = -lattice_sd(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)
m(:,:,3) =  lattice_st(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)
m(:,:,4) = -lattice_st(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)

F = math_mul33x33(Fe(:,:,g,ip,el), Fp(:,:,g,ip,el))
detFe = math_det3x3(Fe(:,:,g,ip,el)) 

fluxdensity = rhoSgl(:,1:4) * constitutive_nonlocal_v(:,:,g,ip,el)
  
!if ((debug_g==g .and. debug_i==ip .and. debug_e==el)) write(6,*) '--> dislocation flux <---'
do n = 1,FE_NipNeighbors(mesh_element(2,el))                                                                                        ! loop through my neighbors

  neighboring_el = mesh_ipNeighborhood(1,n,ip,el)
  neighboring_ip = mesh_ipNeighborhood(2,n,ip,el)
  
  opposite_n = n + mod(n,2) - mod(n+1,2)
  opposite_el = mesh_ipNeighborhood(1,opposite_n,ip,el)
  opposite_ip = mesh_ipNeighborhood(2,opposite_n,ip,el)

  if ( neighboring_el > 0_pInt .and. neighboring_ip > 0_pInt ) then                                                                 ! if neighbor exists, average deformation gradient
    neighboring_F = math_mul33x33(Fe(:,:,g,neighboring_ip,neighboring_el), Fp(:,:,g,neighboring_ip,neighboring_el))
    Favg = 0.5_pReal * (F + neighboring_F)
  else                                                                                                                              ! if no neighbor, take my value as average
    Favg = F
  endif
  
  surfaceNormal_currentconf = math_det3x3(Favg) * math_mul33x3(math_inv3x3(transpose(Favg)), mesh_ipAreaNormal(:,n,ip,el))          ! calculate the normal of the interface in current ...
  surfaceNormal = math_mul33x3(transpose(Fe(:,:,g,ip,el)), surfaceNormal_currentconf) / detFe                                       ! ... and lattice configuration
  area = mesh_ipArea(n,ip,el) * math_norm3(surfaceNormal)
  surfaceNormal = surfaceNormal / math_norm3(surfaceNormal)                                                                         ! normalize the surface normal to unit length
    
  neighboring_rhoDotFlux = 0.0_pReal
!  if ((debug_g==g .and. debug_i==ip .and. debug_e==el)) write(6,'(a,x,i2)') 'neighbor',n
  do s = 1,ns
!    if ((debug_g==g .and. debug_i==ip .and. debug_e==el)) write(6,'(a,x,i2)') '  system',s
    do t = 1,4
!      if ((debug_g==g .and. debug_i==ip .and. debug_e==el)) write(6,'(a,x,i2)') '    type',t
      c = (t + 1) / 2
      topp = t + mod(t,2) - mod(t+1,2)
      
      if ( abs(math_mul3x3(m(:,s,t),surfaceNormal)) > 0.01_pReal &
          .and. fluxdensity(s,t) * math_mul3x3(m(:,s,t),surfaceNormal) > 0.0_pReal ) then                                           ! outgoing flux        
        
        lineLength = fluxdensity(s,t) * math_mul3x3(m(:,s,t),surfaceNormal) * area                                                  ! line length that wants to leave thrugh this interface
        
        if ( (opposite_el > 0 .and. opposite_ip > 0) &
             .or. .not. all(periodicSurfaceFlux(maxloc(abs(mesh_ipAreaNormal(:,opposite_n,ip,el))))) ) then
          rhoDotFlux(s,t) = rhoDotFlux(s,t) - lineLength / mesh_ipVolume(ip,el)                                                     ! subtract dislocation flux from cuurent mobile type
!          if ((debug_g==g .and. debug_i==ip .and. debug_e==el)) write(6,'(a,x,e12.5)') '      outgoing flux:', lineLength / mesh_ipVolume(ip,el)
        endif
        rhoDotFlux(s,t+4) = rhoDotFlux(s,t+4) + lineLength / mesh_ipVolume(ip,el) &
                          * (1.0_pReal - sum(constitutive_nonlocal_compatibility(c,:,s,n,ip,el)**2.0_pReal)) &
                          * sign(1.0_pReal, fluxdensity(s,t))                                                                       ! dislocation flux that is not able to leave through interface (because of low transmissivity) will remain as immobile single density at the material point        
        
        if (neighboring_el > 0 .and. neighboring_ip > 0) then                                                                       ! neighbor present
          where (constitutive_nonlocal_compatibility(c,:,s,n,ip,el) > 0.0_pReal) &                                                  ! ..positive compatibility
            neighboring_rhoDotFlux(:,t) = neighboring_rhoDotFlux(:,t) &                                                             ! ....transferring to equally signed dislocation type at neighbor
                                        + lineLength / mesh_ipVolume(neighboring_ip,neighboring_el) &
                                                     * constitutive_nonlocal_compatibility(c,:,s,n,ip,el) ** 2.0_pReal
          where (constitutive_nonlocal_compatibility(c,:,s,n,ip,el) < 0.0_pReal) &                                                  ! ..negative compatibility
            neighboring_rhoDotFlux(:,topp) = neighboring_rhoDotFlux(:,topp) &                                                       ! ....transferring to opposite signed dislocation type at neighbor
                                           + lineLength / mesh_ipVolume(neighboring_ip,neighboring_el) &
                                                        * constitutive_nonlocal_compatibility(c,:,s,n,ip,el) ** 2.0_pReal
!          if ((debug_g==g .and. debug_i==ip .and. debug_e==el)) write(6,'(a,x,e12.5)') '      entering flux at neighbor:', lineLength / mesh_ipVolume(ip,el) &
!                * sum(constitutive_nonlocal_compatibility(c,:,s,n,ip,el) ** 2.0_pReal)
        endif
        
      endif

    enddo ! dislocation type loop
  enddo ! slip system loop

  if (any(abs(neighboring_rhoDotFlux) > 10.0_pReal)) then                                                                            ! only significant density change in neighbr is considered
    !$OMP CRITICAL (fluxes)
      constitutive_nonlocal_rhoDotFlux(:,:,g,neighboring_ip,neighboring_el) = &
        constitutive_nonlocal_rhoDotFlux(:,:,g,neighboring_ip,neighboring_el) + neighboring_rhoDotFlux
      dotState(g,neighboring_ip,neighboring_el)%p(1:10*ns) = &
        dotState(g,neighboring_ip,neighboring_el)%p(1:10*ns) + reshape(neighboring_rhoDotFlux,(/10*ns/))
    !$OMP END CRITICAL (fluxes)
  else
    neighboring_rhoDotFlux = 0.0_pReal
  endif

enddo ! neighbor loop

if (any(abs(rhoDotFlux) > 0.0_pReal)) then
  !$OMP CRITICAL (fluxes)
    constitutive_nonlocal_rhoDotFlux(:,:,g,ip,el) = constitutive_nonlocal_rhoDotFlux(:,:,g,ip,el) + rhoDotFlux
  !$OMP END CRITICAL (fluxes)
endif


!****************************************************************************
!*** calculate dipole formation and annihilation

!*** formation by glide

do c = 1,2

  rhoDotSingle2DipoleGlide(:,2*c-1) = - 2.0_pReal * dUpper(:,c) / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                                                  * (rhoSgl(:,2*c-1) * abs(gdot(:,2*c)) + rhoSgl(:,2*c) * abs(gdot(:,2*c-1)) &      ! negative mobile <-> positive mobile
                                                                                        + abs(rhoSgl(:,2*c+4)) * abs(gdot(:,2*c-1)))! negative immobile <-> positive mobile

  rhoDotSingle2DipoleGlide(:,2*c) = - 2.0_pReal * dUpper(:,c) / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                                                * (rhoSgl(:,2*c-1) * abs(gdot(:,2*c)) + rhoSgl(:,2*c) * abs(gdot(:,2*c-1)) &        ! negative mobile <-> positive mobile
                                                                                      + abs(rhoSgl(:,2*c+3)) * abs(gdot(:,2*c)))    ! negative mobile <-> positive immobile

  rhoDotSingle2DipoleGlide(:,2*c+3) = - 2.0_pReal * dUpper(:,c) / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &        ! negative mobile <-> positive immobile
                                                  * rhoSgl(:,2*c+3) * abs(gdot(:,2*c))

  rhoDotSingle2DipoleGlide(:,2*c+4) = - 2.0_pReal * dUpper(:,c) / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                                                  * rhoSgl(:,2*c+4) * abs(gdot(:,2*c-1))                                            ! negative immobile <-> positive mobile

  rhoDotSingle2DipoleGlide(:,c+8) = - rhoDotSingle2DipoleGlide(:,2*c-1) - rhoDotSingle2DipoleGlide(:,2*c) &
                                    + abs(rhoDotSingle2DipoleGlide(:,2*c+3)) + abs(rhoDotSingle2DipoleGlide(:,2*c+4))
enddo


!*** athermal annihilation

rhoDotAthermalAnnihilation = 0.0_pReal

forall (c=1:2) &  
  rhoDotAthermalAnnihilation(:,c+8) = - 2.0_pReal * dLower(:,c) / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                         * (  2.0_pReal * ( rhoSgl(:,2*c-1) * abs(gdot(:,2*c)) + rhoSgl(:,2*c) * abs(gdot(:,2*c-1)) ) &             ! was single hitting single
                            + 2.0_pReal * ( abs(rhoSgl(:,2*c+3)) * abs(gdot(:,2*c)) + abs(rhoSgl(:,2*c+4)) * abs(gdot(:,2*c-1)) ) & ! was single hitting immobile single or was immobile single hit by single
                            + rhoDip(:,c) * ( abs(gdot(:,2*c-1)) + abs(gdot(:,2*c)) ) )                                             ! single knocks dipole constituent
  
  
!*** thermally activated annihilation of dipoles

rhoDotThermalAnnihilation = 0.0_pReal

D = constitutive_nonlocal_D0(myInstance) * exp(-constitutive_nonlocal_Qsd(myInstance) / (kB * Temperature))

vClimb =  constitutive_nonlocal_atomicVolume(myInstance) * D / ( kB * Temperature ) &
          * constitutive_nonlocal_Gmod(myInstance) / ( 2.0_pReal * pi * (1.0_pReal-constitutive_nonlocal_nu(myInstance)) ) &
          * 2.0_pReal / ( dUpper(:,1) + dLower(:,1) )
          
rhoDotThermalAnnihilation(:,9) = - 4.0_pReal * rhoDip(:,1) * vClimb / ( dUpper(:,1) - dLower(:,1) )                                 ! edge climb
rhoDotThermalAnnihilation(:,10) = 0.0_pReal                                                                                         !!! cross slipping still has to be implemented !!!


!*** formation/dissociation by stress change = alteration in dUpper

rhoDotDipole2SingleStressChange = 0.0_pReal
forall (c=1:2, s=1:ns, dUpperDot(s,c) < 0.0_pReal) &                                                                                ! increased stress => dipole dissociation
  rhoDotDipole2SingleStressChange(s,8+c) = rhoDip(s,c) * dUpperDot(s,c) / (previousDUpper(s,c) - dLower(s,c))

forall (t=1:4) &
  rhoDotDipole2SingleStressChange(:,t) = -0.5_pReal * rhoDotDipole2SingleStressChange(:,(t-1)/2+9)


rhoDotSingle2DipoleStressChange = 0.0_pReal
do c = 1,2
  do s = 1,ns
    if (dUpperDot(s,c) > 0.0_pReal) then                                                                                            ! stress decrease => dipole formation
      rhoDotSingle2DipoleStressChange(s,2*(c-1)+1) = -4.0_pReal * dUpperDot(s,c) * previousDUpper(s,c) * rhoSgl(s,2*(c-1)+1) &
                                                                * ( rhoSgl(s,2*(c-1)+2) + abs(rhoSgl(s,2*(c-1)+6)) )
      rhoDotSingle2DipoleStressChange(s,2*(c-1)+2) = -4.0_pReal * dUpperDot(s,c) * previousDUpper(s,c) * rhoSgl(s,2*(c-1)+2) &
                                                                * ( rhoSgl(s,2*(c-1)+1) + abs(rhoSgl(s,2*(c-1)+5)) )
      rhoDotSingle2DipoleStressChange(s,2*(c-1)+5) = -4.0_pReal * dUpperDot(s,c) * previousDUpper(s,c) * rhoSgl(s,2*(c-1)+5) &
                                                                * ( rhoSgl(s,2*(c-1)+2) + abs(rhoSgl(s,2*(c-1)+6)) )
      rhoDotSingle2DipoleStressChange(s,2*(c-1)+6) = -4.0_pReal * dUpperDot(s,c) * previousDUpper(s,c) * rhoSgl(s,2*(c-1)+6) &
                                                                * ( rhoSgl(s,2*(c-1)+1) + abs(rhoSgl(s,2*(c-1)+5)) )
    endif
  enddo
enddo

forall (c = 1:2) &
    rhoDotSingle2DipoleStressChange(:,8+c) =   abs(rhoDotSingle2DipoleStressChange(:,2*(c-1)+1)) &
                                             + abs(rhoDotSingle2DipoleStressChange(:,2*(c-1)+2)) &
                                             + abs(rhoDotSingle2DipoleStressChange(:,2*(c-1)+5)) &
                                             + abs(rhoDotSingle2DipoleStressChange(:,2*(c-1)+6))


!****************************************************************************
!*** assign the rates of dislocation densities to my dotState

rhoDot = 0.0_pReal
forall (t = 1:10) &
  rhoDot(:,t) = rhoDotFlux(:,t) &
              + rhoDotMultiplication(:,t) &
              + rhoDotRemobilization(:,t) &
              + rhoDotSingle2DipoleGlide(:,t) &
              + rhoDotAthermalAnnihilation(:,t) &
              + rhoDotThermalAnnihilation(:,t) 
!              + rhoDotDipole2SingleStressChange(:,t) 
!              + rhoDotSingle2DipoleStressChange(:,t)

if (verboseDebugger .and. (debug_g==g .and. debug_i==ip .and. debug_e==el)) then
  !$OMP CRITICAL (write2out)
    write(6,'(a,/,8(12(e12.5,x),/))') 'dislocation remobilization', rhoDotRemobilization(:,1:8) * timestep
    write(6,'(a,/,4(12(e12.5,x),/))') 'dislocation multiplication', rhoDotMultiplication(:,1:4) * timestep
    write(6,'(a,/,8(12(e12.5,x),/))') 'dislocation flux (outgoing)', rhoDotFlux(:,1:8) * timestep
    write(6,'(a,/,10(12(e12.5,x),/))') 'dipole formation by glide', rhoDotSingle2DipoleGlide * timestep
    write(6,'(a,/,2(12(e12.5,x),/))') 'athermal dipole annihilation', rhoDotAthermalAnnihilation(:,1:2) * timestep
    write(6,'(a,/,2(12(e12.5,x),/))') 'thermally activated dipole annihilation', rhoDotThermalAnnihilation(:,9:10) * timestep
!    write(6,'(a,/,10(12(e12.5,x),/))') 'dipole dissociation by stress increase', rhoDotDipole2SingleStressChange * timestep
!    write(6,'(a,/,10(12(e12.5,x),/))') 'dipole formation by stress decrease', rhoDotSingle2DipoleStressChange * timestep
    write(6,'(a,/,10(12(e12.5,x),/))') 'total density change', rhoDot * timestep
    write(6,'(a,/,10(12(f12.7,x),/))') 'relative density change', rhoDot(:,1:8) * timestep / (abs(rhoSgl)+1.0e-10), &
                                                                  rhoDot(:,9:10) * timestep / (rhoDip+1.0e-10)
    write(6,*)
  !$OMP END CRITICAL (write2out)
endif

!$OMP CRITICAL (copy2dotState)
  dotState(g,ip,el)%p(1:10*ns) = dotState(g,ip,el)%p(1:10*ns) + reshape(rhoDot,(/10*ns/))
!$OMP END CRITICAL (copy2dotState)

endsubroutine



!*********************************************************************
!* COMPATIBILITY UPDATE                                              *
!* Compatibility is defined as normalized product of signed cosine   *
!* of the angle between the slip plane normals and signed cosine of  *
!* the angle between the slip directions. Only the largest values    *
!* that sum up to a total of 1 are considered, all others are set to *
!* zero.                                                             *
!*********************************************************************
subroutine constitutive_nonlocal_updateCompatibility(orientation,i,e)

use prec,     only:   pReal, &
                      pInt
use math, only:       math_QuaternionDisorientation, &
                      math_mul3x3, &
                      math_qRot
use material, only:   material_phase, &
                      phase_constitution, &
                      phase_localConstitution, &
                      phase_constitutionInstance, &
                      homogenization_maxNgrains
use mesh, only:       mesh_element, &
                      mesh_ipNeighborhood, &
                      FE_NipNeighbors, &
                      mesh_maxNips, &
                      mesh_NcpElems
use lattice, only:    lattice_sn, &
                      lattice_sd, &
                      lattice_st
use debug, only:      debugger, &
                      debug_e, debug_i, debug_g, &
                      verboseDebugger

implicit none

!* input variables
integer(pInt), intent(in) ::                    i, &                          ! ip index
                                                e                             ! element index
real(pReal), dimension(4,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                                orientation                   ! crystal orientation in quaternions
                                            
!* output variables

!* local variables
integer(pInt)                                   n, &                          ! neighbor index 
                                                neighboring_e, &              ! element index of my neighbor
                                                neighboring_i, &              ! integration point index of my neighbor
                                                myPhase, &                    ! phase 
                                                neighboringPhase, &
                                                myInstance, &                 ! instance of constitution
                                                neighboringInstance, &
                                                myStructure, &                ! lattice structure
                                                neighboringStructure, &
                                                myNSlipSystems, &             ! number of active slip systems
                                                neighboringNSlipSystems, &
                                                s1, &                         ! slip system index (me)
                                                s2                            ! slip system index (my neighbor)
integer(pInt), dimension(maxval(constitutive_nonlocal_totalNslip)) :: &
                                                mySlipSystems, &              ! slip system numbering according to lattice 
                                                neighboringSlipSystems
real(pReal), dimension(4) ::                    absoluteMisorientation        ! absolute misorientation (without symmetry) between me and my neighbor
real(pReal), dimension(3,maxval(constitutive_nonlocal_totalNslip)) :: &  
                                                myNormals, &                  ! slip plane normals
                                                neighboringNormals, &
                                                mySlipDirections, &           ! slip directions
                                                neighboringSlipDirections
real(pReal)                                     compatibilitySum, &
                                                compatibilityMax, &
                                                compatibilityMaxCount
logical, dimension(maxval(constitutive_nonlocal_totalNslip)) :: & 
                                                compatibilityMask


myPhase = material_phase(1,i,e)
myInstance = phase_constitutionInstance(myPhase)
myStructure = constitutive_nonlocal_structure(myInstance)
myNSlipSystems = constitutive_nonlocal_totalNslip(myInstance)
mySlipSystems(1:myNSlipSystems) = constitutive_nonlocal_slipSystemLattice(1:myNSlipSystems,myInstance)
myNormals = lattice_sn(:, mySlipSystems, myStructure)
mySlipDirections = lattice_sd(:, mySlipSystems, myStructure)

do n = 1,FE_NipNeighbors(mesh_element(2,e))                                                               ! loop through my neighbors          
  neighboring_e = mesh_ipNeighborhood(1,n,i,e)
  neighboring_i = mesh_ipNeighborhood(2,n,i,e)    
  if ((neighboring_e > 0) .and. (neighboring_i > 0)) then                                                 ! if neighbor exists
    neighboringPhase = material_phase(1,neighboring_i,neighboring_e)
    
    if (.not. phase_localConstitution(neighboringPhase)) then                                             ! neighbor got also nonlocal constitution
      neighboringInstance = phase_constitutionInstance(neighboringPhase)        
      neighboringStructure = constitutive_nonlocal_structure(neighboringInstance)
      neighboringNSlipSystems = constitutive_nonlocal_totalNslip(neighboringInstance)
      neighboringSlipSystems(1:neighboringNSlipSystems) = constitutive_nonlocal_slipSystemLattice(1:neighboringNSlipSystems,&
                                                                                                  neighboringInstance)
      neighboringNormals = lattice_sn(:, neighboringSlipSystems, neighboringStructure)
      neighboringSlipDirections = lattice_sd(:, neighboringSlipSystems, neighboringStructure)
  
      if (myStructure == neighboringStructure) then                                                       ! if my neighbor has same crystal structure like me
        absoluteMisorientation = math_QuaternionDisorientation( orientation(:,1,i,e), &
                                                                orientation(:,1,neighboring_i,neighboring_e), 0_pInt)
  
        do s1 = 1,myNSlipSystems                                                                          ! loop through my slip systems
          do s2 = 1,neighboringNSlipSystems                                                               ! loop through my neighbors' slip systems
            constitutive_nonlocal_compatibility(1,s2,s1,n,i,e) = &
                  math_mul3x3(myNormals(:,s1), math_qRot(absoluteMisorientation, neighboringNormals(:,s2))) &
                * abs(math_mul3x3(mySlipDirections(:,s1), math_qRot(absoluteMisorientation, neighboringSlipDirections(:,s2))))
            constitutive_nonlocal_compatibility(2,s2,s1,n,i,e) = &
                  abs(math_mul3x3(myNormals(:,s1), math_qRot(absoluteMisorientation, neighboringNormals(:,s2)))) &
                * abs(math_mul3x3(mySlipDirections(:,s1), math_qRot(absoluteMisorientation, neighboringSlipDirections(:,s2))))
          enddo
          compatibilitySum = 0.0_pReal
          compatibilityMask = .true.
          do while ( (1.0_pReal - compatibilitySum > 0.0_pReal) .and. any(compatibilityMask) )            ! only those largest values that sum up to 1 are considered (round off of the smallest considered values to ensure sum to be exactly 1)
            compatibilityMax = maxval(constitutive_nonlocal_compatibility(2,:,s1,n,i,e), compatibilityMask) ! screws always positive
            compatibilityMaxCount = dble(count(constitutive_nonlocal_compatibility(2,:,s1,n,i,e) == compatibilityMax))
            where (constitutive_nonlocal_compatibility(2,:,s1,n,i,e) >= compatibilityMax) compatibilityMask = .false.
            if (compatibilitySum + compatibilityMax * compatibilityMaxCount > 1.0_pReal) &                ! if compatibility sum exceeds 1...
              where (abs(constitutive_nonlocal_compatibility(:,:,s1,n,i,e)) == compatibilityMax) &        ! ... equally distribute what is left
                constitutive_nonlocal_compatibility(:,:,s1,n,i,e) = sign((1.0_pReal - compatibilitySum) / compatibilityMaxCount, &
                                                                         constitutive_nonlocal_compatibility(:,:,s1,n,i,e))
            compatibilitySum = compatibilitySum + compatibilityMaxCount * compatibilityMax
          enddo
          where (compatibilityMask) constitutive_nonlocal_compatibility(1,:,s1,n,i,e) = 0.0_pReal
          where (compatibilityMask) constitutive_nonlocal_compatibility(2,:,s1,n,i,e) = 0.0_pReal
        enddo          
      else                                                                                                ! neighbor has different crystal structure
        constitutive_nonlocal_compatibility(:,:,:,n,i,e) = 0.0_pReal                                      ! no compatibility              
      endif
    else                                                                                                  ! neighbor has local constitution
      constitutive_nonlocal_compatibility(:,:,:,n,i,e) = 0.0_pReal
      forall(s1 = 1:maxval(constitutive_nonlocal_totalNslip)) &
        constitutive_nonlocal_compatibility(:,s1,s1,n,i,e) = 1.0_pReal                                    ! assume perfect compatibility for equal slip system index
    endif
  else                                                                                                    ! no neighbor present
    constitutive_nonlocal_compatibility(:,:,:,n,i,e) = 0.0_pReal
    forall(s1 = 1:maxval(constitutive_nonlocal_totalNslip)) &
      constitutive_nonlocal_compatibility(:,s1,s1,n,i,e) = 1.0_pReal                                      ! perfect compatibility for equal slip system index
  endif

enddo
  
endsubroutine 



!*********************************************************************
!* rate of change of temperature                                     *
!*********************************************************************
pure function constitutive_nonlocal_dotTemperature(Tstar_v,Temperature,state,g,ip,el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips
use material, only: homogenization_maxNgrains
implicit none

!* input variables
integer(pInt), intent(in) ::              g, &              ! current grain ID
                                          ip, &             ! current integration point
                                          el                ! current element
real(pReal), intent(in) ::                Temperature       ! temperature
real(pReal), dimension(6), intent(in) ::  Tstar_v           ! 2nd Piola-Kirchhoff stress in Mandel notation
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: & 
                                          state             ! microstructural state

!* output variables
real(pReal) constitutive_nonlocal_dotTemperature        ! evolution of Temperature

!* local variables
   
constitutive_nonlocal_dotTemperature = 0.0_pReal

endfunction



!*********************************************************************
!* return array of constitutive results                              *
!*********************************************************************
function constitutive_nonlocal_postResults(Tstar_v, previousTstar_v, Fe, Fp, Temperature, disorientation, dt, dt_previous, &
                                                state, previousState, dotState, g,ip,el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use math,     only: math_norm3, &
                    math_mul6x6, &
                    math_mul3x3, &
                    math_mul33x3, &
                    math_mul33x33, &
                    math_inv3x3, &
                    math_det3x3, &
                    math_Mandel6to33, &
                    pi
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_maxNipNeighbors, &
                    mesh_element, &
                    FE_NipNeighbors, &
                    mesh_ipNeighborhood, &
                    mesh_ipVolume, &
                    mesh_ipArea, &
                    mesh_ipAreaNormal
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_constitutionInstance, &
                    phase_Noutput
use lattice,  only: lattice_Sslip, &
                    lattice_Sslip_v, &
                    lattice_sd, &
                    lattice_sn, &
                    lattice_st, &
                    lattice_maxNslipFamily, &
                    lattice_NslipSystem
implicit none

!*** input variables
integer(pInt), intent(in) ::                g, &                      ! current grain number
                                            ip, &                     ! current integration point
                                            el                        ! current element number
real(pReal), intent(in) ::                  Temperature, &            ! temperature
                                            dt, &                     ! time increment
                                            dt_previous               ! time increment between previous and current state
real(pReal), dimension(6), intent(in) ::    Tstar_v, &                ! current 2nd Piola-Kirchhoff stress in Mandel notation
                                            previousTstar_v           ! previous 2nd Piola-Kirchhoff stress in Mandel notation
real(pReal), dimension(4,mesh_maxNipNeighbors), intent(in) :: &
                                            disorientation            ! crystal disorientation between me and my neighbor (axis, angle pair)
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            Fe, &                     ! elastic deformation gradient
                                            Fp                        ! plastic deformation gradient
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            state, &                  ! current microstructural state
                                            previousState, &          ! previous microstructural state
                                            dotState                  ! evolution rate of microstructural state

!*** output variables
real(pReal), dimension(constitutive_nonlocal_sizePostResults(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            constitutive_nonlocal_postResults

!*** local variables
integer(pInt)                               myInstance, &             ! current instance of this constitution
                                            myStructure, &            ! current lattice structure
                                            ns, &                     ! short notation for the total number of active slip systems
                                            neighboring_el, &         ! element number of my neighbor
                                            neighboring_ip, &         ! integration point of my neighbor
                                            c, &                      ! character of dislocation
                                            cs, &                     ! constitutive result index
                                            n, &                      ! index of my current neighbor
                                            o, &                      ! index of current output
                                            t, &                      ! type of dislocation
                                            s, &                      ! index of my current slip system
                                            sLattice                  ! index of my current slip system according to lattice order
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),6,4) :: &
                                            fluxes                    ! outgoing fluxes per slipsystem, neighbor and dislocation type
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),8) :: &
                                            rhoSgl, &                 ! current single dislocation densities (positive/negative screw and edge without dipoles)
                                            previousRhoSgl, &         ! previous single dislocation densities (positive/negative screw and edge without dipoles)
                                            rhoDotSgl                 ! evolution rate of single dislocation densities (positive/negative screw and edge without dipoles)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            gdot, &                   ! shear rates
                                            lineLength                ! dislocation line length leaving the current interface
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            rhoForest, &              ! forest dislocation density
                                            tauThreshold, &           ! threshold shear stress
                                            tau, &                    ! current resolved shear stress
                                            previousTau, &            ! previous resolved shear stress
                                            invLambda, &              ! inverse of mean free path for dislocations
                                            vClimb                    ! climb velocity of edge dipoles
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                            rhoDip, &                 ! current dipole dislocation densities (screw and edge dipoles)
                                            previousRhoDip, &         ! previous dipole dislocation densities (screw and edge dipoles)
                                            rhoDotDip, &              ! evolution rate of dipole dislocation densities (screw and edge dipoles)
                                            dLower, &                 ! minimum stable dipole distance for edges and screws
                                            dUpper, &                 ! current maximum stable dipole distance for edges and screws
                                            previousDUpper, &         ! previous maximum stable dipole distance for edges and screws
                                            dUpperDot                 ! rate of change of the maximum stable dipole distance for edges and screws
real(pReal), dimension(3,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                            m, &                      ! direction of dislocation motion for edge and screw (unit vector)
                                            m_currentconf             ! direction of dislocation motion for edge and screw (unit vector) in current configuration
real(pReal), dimension(3,3) ::              F, &                      ! total deformation gradient
                                            neighboring_F, &          ! total deformation gradient of my neighbor
                                            Favg                      ! average total deformation gradient of me and my neighbor
real(pReal), dimension(6) ::                Tdislocation_v, &         ! current dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress
                                            previousTdislocation_v    ! previous dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress
real(pReal), dimension(3) ::                surfaceNormal, &          ! surface normal in lattice configuration
                                            surfaceNormal_currentconf ! surface normal in current configuration
real(pReal)                                 area, &                   ! area of the current interface
                                            detFe, &                  ! determinant of elastic defornmation gradient
                                            D                         ! self diffusion


myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)

cs = 0_pInt
constitutive_nonlocal_postResults = 0.0_pReal


!* short hand notations for state variables

forall (t = 1:8) rhoSgl(:,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (t = 1:8) previousRhoSgl(:,t) = previousState(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (c = 1:2) rhoDip(:,c) = state(g,ip,el)%p((7+c)*ns+1:(8+c)*ns)
forall (c = 1:2) previousRhoDip(:,c) = previousState(g,ip,el)%p((7+c)*ns+1:(8+c)*ns)
rhoForest = state(g,ip,el)%p(10*ns+1:11*ns)
tauThreshold = state(g,ip,el)%p(11*ns+1:12*ns)
Tdislocation_v = state(g,ip,el)%p(12*ns+1:12*ns+6)
previousTdislocation_v = previousState(g,ip,el)%p(12*ns+1:12*ns+6)
forall (t = 1:8) rhoDotSgl(:,t) = dotState(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (c = 1:2) rhoDotDip(:,c) = dotState(g,ip,el)%p((7+c)*ns+1:(8+c)*ns)


!* Calculate shear rate

do t = 1,4
  do s = 1,ns
    if (rhoSgl(s,t+4) * constitutive_nonlocal_v(s,t,g,ip,el) < 0.0_pReal) then
      rhoSgl(s,t) = rhoSgl(s,t) + abs(rhoSgl(s,t+4))                                                                  ! remobilization of immobile singles for changing sign of v (bauschinger effect)
      rhoSgl(s,t+4) = 0.0_pReal                                                                                       ! remobilization of immobile singles for changing sign of v (bauschinger effect)
    endif
  enddo
enddo

forall (t = 1:4) &
  gdot(:,t) =  rhoSgl(:,t) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) * constitutive_nonlocal_v(:,t,g,ip,el)
  
!* calculate limits for stable dipole height and its rate of change

do s = 1,ns
  sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
  tau(s) = math_mul6x6( Tstar_v + Tdislocation_v, lattice_Sslip_v(:,sLattice,myStructure) )
  previousTau(s) = math_mul6x6( previousTstar_v + previousTdislocation_v, lattice_Sslip_v(:,sLattice,myStructure) )
enddo

dLower(:,1) = constitutive_nonlocal_dLowerEdgePerSlipSystem(:,myInstance)
dLower(:,2) = constitutive_nonlocal_dLowerScrewPerSlipSystem(:,myInstance)
dUpper(:,2) = min( constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                                                           / ( 8.0_pReal * pi * abs(tau) ), &
                   1.0_pReal / sqrt( sum(abs(rhoSgl),2)+sum(rhoDip,2) ) )
dUpper(:,1) = dUpper(:,2) / ( 1.0_pReal - constitutive_nonlocal_nu(myInstance) )
previousDUpper(:,2) = min( constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                                                            / ( 8.0_pReal * pi * abs(previousTau) ), &
                    1.0_pReal / sqrt(  sum(abs(previousRhoSgl),2) + sum(previousRhoDip,2) ) )
previousDUpper(:,1) = previousDUpper(:,2) / ( 1.0_pReal - constitutive_nonlocal_nu(myInstance) )

if (dt_previous > 0.0_pReal) then
  dUpperDot = (dUpper - previousDUpper) / dt_previous
else 
  dUpperDot = 0.0_pReal
endif


!*** dislocation motion

m(:,:,1) = lattice_sd(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)
m(:,:,2) = lattice_st(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)
forall (c = 1:2, s = 1:ns) &
  m_currentconf(:,s,c) = math_mul33x3(Fe(:,:,g,ip,el), m(:,s,c))


do o = 1,phase_Noutput(material_phase(g,ip,el))
  select case(constitutive_nonlocal_output(o,myInstance))
    
    case ('rho')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl),2) + sum(rhoDip,2)
      cs = cs + ns
      
    case ('rho_sgl')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl),2)
      cs = cs + ns
      
    case ('rho_sgl_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(:,1:4)),2)
      cs = cs + ns
      
    case ('rho_sgl_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(:,5:8)),2)
      cs = cs + ns
      
    case ('rho_dip')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoDip,2)
      cs = cs + ns
      
    case ('rho_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(:,(/1,2,5,6/))),2) + rhoDip(:,1)
      cs = cs + ns
      
    case ('rho_sgl_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(:,(/1,2,5,6/))),2)
      cs = cs + ns
      
    case ('rho_sgl_edge_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoSgl(:,1:2),2)
      cs = cs + ns
      
    case ('rho_sgl_edge_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(:,5:6)),2)
      cs = cs + ns
      
    case ('rho_sgl_edge_pos')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,1) + abs(rhoSgl(:,5))
      cs = cs + ns
      
    case ('rho_sgl_edge_pos_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,1)
      cs = cs + ns
      
    case ('rho_sgl_edge_pos_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = abs(rhoSgl(:,5))
      cs = cs + ns
      
    case ('rho_sgl_edge_neg')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,2) + abs(rhoSgl(:,6))
      cs = cs + ns
      
    case ('rho_sgl_edge_neg_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,2)
      cs = cs + ns
      
    case ('rho_sgl_edge_neg_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = abs(rhoSgl(:,6))
      cs = cs + ns
      
    case ('rho_dip_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoDip(:,1)
      cs = cs + ns
      
    case ('rho_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(:,(/3,4,7,8/))),2) + rhoDip(:,2)
      cs = cs + ns
      
    case ('rho_sgl_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(:,(/3,4,7,8/))),2)
      cs = cs + ns
            
    case ('rho_sgl_screw_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoSgl(:,3:4),2)
      cs = cs + ns
      
    case ('rho_sgl_screw_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(:,7:8)),2)
      cs = cs + ns
      
    case ('rho_sgl_screw_pos')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,3) + abs(rhoSgl(:,7))
      cs = cs + ns
      
    case ('rho_sgl_screw_pos_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,3)
      cs = cs + ns
      
    case ('rho_sgl_screw_pos_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = abs(rhoSgl(:,7))
      cs = cs + ns
      
    case ('rho_sgl_screw_neg')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,4) + abs(rhoSgl(:,8))
      cs = cs + ns

    case ('rho_sgl_screw_neg_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,4)
      cs = cs + ns

    case ('rho_sgl_screw_neg_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = abs(rhoSgl(:,8))
      cs = cs + ns

    case ('rho_dip_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoDip(:,2)
      cs = cs + ns
      
    case ('excess_rho')
      constitutive_nonlocal_postResults(cs+1:cs+ns) =   (rhoSgl(:,1) + abs(rhoSgl(:,5))) - (rhoSgl(:,2) + abs(rhoSgl(:,6))) &
                                                      + (rhoSgl(:,3) + abs(rhoSgl(:,7))) - (rhoSgl(:,4) + abs(rhoSgl(:,8)))
      cs = cs + ns
      
    case ('excess_rho_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = (rhoSgl(:,1) + abs(rhoSgl(:,5))) - (rhoSgl(:,2) + abs(rhoSgl(:,6)))
      cs = cs + ns
      
    case ('excess_rho_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = (rhoSgl(:,3) + abs(rhoSgl(:,7))) - (rhoSgl(:,4) + abs(rhoSgl(:,8)))
      cs = cs + ns
      
    case ('rho_forest')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoForest
      cs = cs + ns
    
    case ('delta')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = 1.0_pReal / sqrt( sum(abs(rhoSgl),2) + sum(rhoDip,2) )
      cs = cs + ns
      
    case ('delta_sgl')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = 1.0_pReal / sqrt( sum(abs(rhoSgl),2))
      cs = cs + ns
      
    case ('delta_dip')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = 1.0_pReal / sqrt( sum(rhoDip,2) )
      cs = cs + ns
      
    case ('shearrate')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(gdot),2)
      cs = cs + ns
      
    case ('resolvedstress')
      do s = 1,ns  
        sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
        constitutive_nonlocal_postResults(cs+s) = math_mul6x6( Tstar_v + Tdislocation_v, lattice_Sslip_v(:,sLattice,myStructure) )
      enddo
      cs = cs + ns
      
    case ('resolvedstress_internal')
      do s = 1,ns  
        sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
        constitutive_nonlocal_postResults(cs+s) = math_mul6x6(Tdislocation_v, lattice_Sslip_v(:,sLattice,myStructure) )
      enddo
      cs = cs + ns
      
    case ('resolvedstress_external')
      do s = 1,ns  
        sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
        constitutive_nonlocal_postResults(cs+s) = math_mul6x6(Tstar_v, lattice_Sslip_v(:,sLattice,myStructure) )
      enddo
      cs = cs + ns
      
    case ('resistance')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = tauThreshold
      cs = cs + ns
    
    case ('rho_dot')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoDotSgl,2) + sum(rhoDotDip,2)
      cs = cs + ns
      
    case ('rho_dot_sgl')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoDotSgl,2)
      cs = cs + ns
      
    case ('rho_dot_dip')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoDotDip,2)
      cs = cs + ns
    
    case ('rho_dot_gen')
      constitutive_nonlocal_postResults(cs+1:cs+ns) =   sum(abs(gdot),2) * sqrt(rhoForest)  &
                                                      / constitutive_nonlocal_lambda0PerSlipSystem(:,myInstance) &
                                                      / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance)
      cs = cs + ns

    case ('rho_dot_gen_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) =   sum(abs(gdot(:,3:4)),2) * sqrt(rhoForest)  &
                                                      / constitutive_nonlocal_lambda0PerSlipSystem(:,myInstance) &
                                                      / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance)
      cs = cs + ns

    case ('rho_dot_gen_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) =   sum(abs(gdot(:,1:2)),2) * sqrt(rhoForest)  &
                                                      / constitutive_nonlocal_lambda0PerSlipSystem(:,myInstance) &
                                                      / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance)
      cs = cs + ns
      
    case ('rho_dot_sgl2dip')
      do c=1,2                                 ! dipole formation by glide
        constitutive_nonlocal_postResults(cs+1:cs+ns) = constitutive_nonlocal_postResults(cs+1:cs+ns) + &
            2.0_pReal * dUpper(:,c) / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                      * (  2.0_pReal * ( rhoSgl(:,2*c-1) * abs(gdot(:,2*c)) + rhoSgl(:,2*c) * abs(gdot(:,2*c-1)) ) &                ! was single hitting single
                         + 2.0_pReal * ( abs(rhoSgl(:,2*c+3)) * abs(gdot(:,2*c)) + abs(rhoSgl(:,2*c+4)) * abs(gdot(:,2*c-1)) ) )    ! was single hitting immobile/used single
      enddo
!      do c=1,2
!        forall (s=1:ns, dUpperDot(s,c) > 0.0_pReal) &    ! dipole formation by stress decrease
!          constitutive_nonlocal_postResults(cs+s) = constitutive_nonlocal_postResults(cs+s) + &
!                                                8.0_pReal * rhoSgl(s,2*c-1) * rhoSgl(s,2*c) * previousDUpper(s,c) * dUpperDot(s,c)
!      enddo
      cs = cs + ns
    
    case ('rho_dot_dip2sgl')
      do c=1,2
        forall (s=1:ns, dUpperDot(s,c) < 0.0_pReal) &
          constitutive_nonlocal_postResults(cs+s) = constitutive_nonlocal_postResults(cs+s) - &
                                                    rhoDip(s,c) * dUpperDot(s,c) / (previousDUpper(s,c) - dLower(s,c))
      enddo
      cs = cs + ns
    
    case ('rho_dot_ann_ath')
      do c=1,2
        constitutive_nonlocal_postResults(cs+1:cs+ns) = constitutive_nonlocal_postResults(cs+1:cs+ns) + &
            2.0_pReal * dLower(:,c) / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                         * (  2.0_pReal * ( rhoSgl(:,2*c-1) * abs(gdot(:,2*c)) + rhoSgl(:,2*c) * abs(gdot(:,2*c-1)) ) &             ! was single hitting single
                            + 2.0_pReal * ( abs(rhoSgl(:,2*c+3)) * abs(gdot(:,2*c)) + abs(rhoSgl(:,2*c+4)) * abs(gdot(:,2*c-1)) ) & ! was single hitting immobile/used single
                            + rhoDip(:,c) * ( abs(gdot(:,2*c-1)) + abs(gdot(:,2*c)) ) )                                             ! single knocks dipole constituent
      enddo
      cs = cs + ns
      
    case ('rho_dot_ann_the') 
      D = constitutive_nonlocal_D0(myInstance) * exp(-constitutive_nonlocal_Qsd(myInstance) / (kB * Temperature))

      vClimb =  constitutive_nonlocal_atomicVolume(myInstance) * D / ( kB * Temperature ) &
          * constitutive_nonlocal_Gmod(myInstance) / ( 2.0_pReal * pi * (1.0_pReal-constitutive_nonlocal_nu(myInstance)) ) &
          * 2.0_pReal / ( dUpper(:,1) + dLower(:,1) )
          
      constitutive_nonlocal_postResults(cs+1:cs+ns) = 4.0_pReal * rhoDip(:,1) * vClimb / ( dUpper(:,1) - dLower(:,1) )
      ! !!! cross-slip of screws missing !!!
      cs = cs + ns

    case ('rho_dot_flux')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(constitutive_nonlocal_rhoDotFlux(:,1:4,g,ip,el),2) &
                                                      + sum(abs(constitutive_nonlocal_rhoDotFlux(:,5:8,g,ip,el)),2)
      cs = cs + ns
    
    case ('rho_dot_flux_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(constitutive_nonlocal_rhoDotFlux(:,1:2,g,ip,el),2) &
                                                      + sum(abs(constitutive_nonlocal_rhoDotFlux(:,5:6,g,ip,el)),2)
      cs = cs + ns
      
    case ('rho_dot_flux_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(constitutive_nonlocal_rhoDotFlux(:,3:4,g,ip,el),2) &
                                                      + sum(abs(constitutive_nonlocal_rhoDotFlux(:,7:8,g,ip,el)),2)
      cs = cs + ns
            
    case ('dislocationvelocity')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = constitutive_nonlocal_v(:,1,g,ip,el)
      cs = cs + ns
    
    case ('fluxdensity_edge_pos_x')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,1) * constitutive_nonlocal_v(:,1,g,ip,el) * m_currentconf(1,:,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_pos_y')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,1) * constitutive_nonlocal_v(:,1,g,ip,el) * m_currentconf(2,:,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_pos_z')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,1) * constitutive_nonlocal_v(:,1,g,ip,el) * m_currentconf(3,:,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_neg_x')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,2) * constitutive_nonlocal_v(:,2,g,ip,el) * m_currentconf(1,:,2)
      cs = cs + ns
    
    case ('fluxdensity_edge_neg_y')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,2) * constitutive_nonlocal_v(:,2,g,ip,el) * m_currentconf(2,:,2)
      cs = cs + ns
    
    case ('fluxdensity_edge_neg_z')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,2) * constitutive_nonlocal_v(:,2,g,ip,el) * m_currentconf(3,:,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_pos_x')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,3) * constitutive_nonlocal_v(:,3,g,ip,el) * m_currentconf(1,:,3)
      cs = cs + ns
    
    case ('fluxdensity_screw_pos_y')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,3) * constitutive_nonlocal_v(:,3,g,ip,el) * m_currentconf(2,:,3)
      cs = cs + ns
    
    case ('fluxdensity_screw_pos_z')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,3) * constitutive_nonlocal_v(:,3,g,ip,el) * m_currentconf(3,:,3)
      cs = cs + ns
    
    case ('fluxdensity_screw_neg_x')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,4) * constitutive_nonlocal_v(:,4,g,ip,el) * m_currentconf(1,:,4)
      cs = cs + ns
    
    case ('fluxdensity_screw_neg_y')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,4) * constitutive_nonlocal_v(:,4,g,ip,el) * m_currentconf(2,:,4)
      cs = cs + ns
    
    case ('fluxdensity_screw_neg_z')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(:,4) * constitutive_nonlocal_v(:,4,g,ip,el) * m_currentconf(3,:,4)
      cs = cs + ns
    
    case ('d_upper_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = dUpper(:,1)
      cs = cs + ns
      
    case ('d_upper_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = dUpper(:,2)
      cs = cs + ns
      
    case ('d_upper_dot_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = dUpperDot(:,1)
      cs = cs + ns
      
    case ('d_upper_dot_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = dUpperDot(:,2)
      cs = cs + ns
      
 end select
enddo

endfunction
END MODULE
