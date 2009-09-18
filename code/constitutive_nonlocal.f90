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
character(len=16), dimension(6), parameter ::     constitutive_nonlocal_listBasicStates = (/'rhoEdgePos      ', &
                                                                                            'rhoEdgeNeg      ', &
                                                                                            'rhoScrewPos     ', &
                                                                                            'rhoScrewNeg     ', &
                                                                                            'rhoEdgeDip      ', &
                                                                                            'rhoScrewDip     ' /) ! list of "basic" microstructural state variables that are independent from other state variables
character(len=16), dimension(3), parameter :: constitutive_nonlocal_listDependentStates = (/'rhoForest       ', &
                                                                                            'tauSlipThreshold', &
                                                                                            'Tdislocation_v  ' /) ! list of microstructural state variables that depend on other state variables
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
                                                          constitutive_nonlocal_atomicVolume, &                 ! atomic volume
                                                          constitutive_nonlocal_D0, &                           ! prefactor for self-diffusion coefficient
                                                          constitutive_nonlocal_Qsd, &                          ! activation energy for dislocation climb
                                                          constitutive_nonlocal_relevantRho                     ! dislocation density considered relevant
real(pReal), dimension(:,:,:), allocatable ::             constitutive_nonlocal_Cslip_66                        ! elasticity matrix in Mandel notation for each instance
real(pReal), dimension(:,:,:,:,:), allocatable ::         constitutive_nonlocal_Cslip_3333                      ! elasticity matrix for each instance
real(pReal), dimension(:,:), allocatable ::               constitutive_nonlocal_rhoEdgePos0, &                  ! initial edge_pos dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoEdgeNeg0, &                  ! initial edge_neg dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoScrewPos0, &                 ! initial screw_pos dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoScrewNeg0, &                 ! initial screw_neg dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoEdgeDip0, &                  ! initial edge dipole dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoScrewDip0, &                 ! initial screw dipole dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_v0PerSlipFamily, &              ! dislocation velocity prefactor [m/s] for each family and instance
                                                          constitutive_nonlocal_v0PerSlipSystem, &              ! dislocation velocity prefactor [m/s] for each slip system and instance
                                                          constitutive_nonlocal_lambda0PerSlipFamily, &         ! mean free path prefactor for each family and instance
                                                          constitutive_nonlocal_lambda0PerSlipSystem, &         ! mean free path prefactor for each slip system and instance
                                                          constitutive_nonlocal_burgersPerSlipFamily, &         ! absolute length of burgers vector [m] for each family and instance
                                                          constitutive_nonlocal_burgersPerSlipSystem, &         ! absolute length of burgers vector [m] for each slip system and instance
                                                          constitutive_nonlocal_dDipMinEdgePerSlipFamily, &     ! minimum stable edge dipole height for each family and instance
                                                          constitutive_nonlocal_dDipMinEdgePerSlipSystem, &     ! minimum stable edge dipole height for each slip system and instance
                                                          constitutive_nonlocal_dDipMinScrewPerSlipFamily, &    ! minimum stable screw dipole height for each family and instance
                                                          constitutive_nonlocal_dDipMinScrewPerSlipSystem, &    ! minimum stable screw dipole height for each slip system and instance
                                                          constitutive_nonlocal_interactionSlipSlip             ! coefficients for slip-slip interaction for each interaction type and instance

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
use material, only: phase_constitution, &
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
allocate(constitutive_nonlocal_atomicVolume(maxNinstance))
allocate(constitutive_nonlocal_D0(maxNinstance))
allocate(constitutive_nonlocal_Qsd(maxNinstance))
allocate(constitutive_nonlocal_relevantRho(maxNinstance))
allocate(constitutive_nonlocal_Cslip_66(6,6,maxNinstance))
allocate(constitutive_nonlocal_Cslip_3333(3,3,3,3,maxNinstance))
constitutive_nonlocal_CoverA = 0.0_pReal 
constitutive_nonlocal_C11 = 0.0_pReal
constitutive_nonlocal_C12 = 0.0_pReal
constitutive_nonlocal_C13 = 0.0_pReal
constitutive_nonlocal_C33 = 0.0_pReal
constitutive_nonlocal_C44 = 0.0_pReal
constitutive_nonlocal_Gmod = 0.0_pReal
constitutive_nonlocal_atomicVolume = 0.0_pReal
constitutive_nonlocal_D0 = 0.0_pReal
constitutive_nonlocal_Qsd = 0.0_pReal
constitutive_nonlocal_relevantRho = 0.0_pReal
constitutive_nonlocal_nu = 0.0_pReal
constitutive_nonlocal_Cslip_66 = 0.0_pReal
constitutive_nonlocal_Cslip_3333 = 0.0_pReal

allocate(constitutive_nonlocal_rhoEdgePos0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoEdgeNeg0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoScrewPos0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoScrewNeg0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoEdgeDip0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoScrewDip0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_v0PerSlipFamily(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_burgersPerSlipFamily(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_Lambda0PerSlipFamily(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_interactionSlipSlip(lattice_maxNinteraction, maxNinstance))
allocate(constitutive_nonlocal_dDipMinEdgePerSlipFamily(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_dDipMinScrewPerSlipFamily(lattice_maxNslipFamily, maxNinstance))
constitutive_nonlocal_rhoEdgePos0 = 0.0_pReal
constitutive_nonlocal_rhoEdgeNeg0 = 0.0_pReal
constitutive_nonlocal_rhoScrewPos0 = 0.0_pReal
constitutive_nonlocal_rhoScrewNeg0 = 0.0_pReal
constitutive_nonlocal_rhoEdgeDip0 = 0.0_pReal
constitutive_nonlocal_rhoScrewDip0 = 0.0_pReal
constitutive_nonlocal_v0PerSlipFamily = 0.0_pReal
constitutive_nonlocal_burgersPerSlipFamily = 0.0_pReal
constitutive_nonlocal_lambda0PerSlipFamily = 0.0_pReal
constitutive_nonlocal_interactionSlipSlip = 0.0_pReal
constitutive_nonlocal_dDipMinEdgePerSlipFamily = 0.0_pReal
constitutive_nonlocal_dDipMinScrewPerSlipFamily = 0.0_pReal


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
      case ('rhoedgepos0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoEdgePos0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhoedgeneg0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoEdgeNeg0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhoscrewpos0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoScrewPos0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhoscrewneg0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoScrewNeg0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhoedgedip0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoEdgeDip0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhoscrewdip0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoScrewDip0(f,i) = IO_floatValue(line,positions,1+f)
      case ('v0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_v0PerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case ('lambda0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_lambda0PerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case ('burgers')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_burgersPerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case('ddipminedge')
        forall (f = 1:lattice_maxNslipFamily) & 
            constitutive_nonlocal_dDipMinEdgePerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case('ddipminscrew')
        forall (f = 1:lattice_maxNslipFamily) & 
            constitutive_nonlocal_dDipMinScrewPerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case('atomicvolume')
        constitutive_nonlocal_atomicVolume(i) = IO_floatValue(line,positions,2)
      case('d0')
        constitutive_nonlocal_D0(i) = IO_floatValue(line,positions,2)
      case('qsd')
        constitutive_nonlocal_Qsd(i) = IO_floatValue(line,positions,2)
      case('relevantrho')
        constitutive_nonlocal_relevantRho(i) = IO_floatValue(line,positions,2)
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
  do f = 1,lattice_maxNslipFamily
    if (constitutive_nonlocal_Nslip(f,i) > 0_pInt) then
      if (constitutive_nonlocal_rhoEdgePos0(f,i) < 0.0_pReal)                               call IO_error(220)
      if (constitutive_nonlocal_rhoEdgeNeg0(f,i) < 0.0_pReal)                               call IO_error(220)
      if (constitutive_nonlocal_rhoScrewPos0(f,i) < 0.0_pReal)                              call IO_error(220)
      if (constitutive_nonlocal_rhoScrewNeg0(f,i) < 0.0_pReal)                              call IO_error(220)
      if (constitutive_nonlocal_rhoEdgeDip0(f,i) < 0.0_pReal)                               call IO_error(220)
      if (constitutive_nonlocal_rhoScrewDip0(f,i) < 0.0_pReal)                              call IO_error(220)
      if (constitutive_nonlocal_burgersPerSlipFamily(f,i) <= 0.0_pReal)                     call IO_error(221)
      if (constitutive_nonlocal_v0PerSlipFamily(f,i) <= 0.0_pReal)                          call IO_error(226)
      if (constitutive_nonlocal_lambda0PerSlipFamily(f,i) <= 0.0_pReal)                     call IO_error(227)
      if (constitutive_nonlocal_dDipMinEdgePerSlipFamily(f,i) <= 0.0_pReal)                 call IO_error(228)
      if (constitutive_nonlocal_dDipMinScrewPerSlipFamily(f,i) <= 0.0_pReal)                call IO_error(228)
    endif
  enddo
  if (any(constitutive_nonlocal_interactionSlipSlip(1:maxval(lattice_interactionSlipSlip(:,:,myStructure)),i) < 1.0_pReal)) &
                                                                                            call IO_error(229)
  if (constitutive_nonlocal_atomicVolume(i) <= 0.0_pReal)                                   call IO_error(230)
  if (constitutive_nonlocal_D0(i) <= 0.0_pReal)                                             call IO_error(231)
  if (constitutive_nonlocal_Qsd(i) <= 0.0_pReal)                                            call IO_error(232)
  if (constitutive_nonlocal_relevantRho(i) <= 0.0_pReal)                                    call IO_error(233)
    
  
!*** determine total number of active slip systems
  
  constitutive_nonlocal_Nslip(:,i) = min( lattice_NslipSystem(:, myStructure), constitutive_nonlocal_Nslip(:,i) )                         ! we can't use more slip systems per family than specified in lattice 
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

allocate(constitutive_nonlocal_dDipMinEdgePerSlipSystem(maxTotalNslip, maxNinstance))
constitutive_nonlocal_dDipMinEdgePerSlipSystem = 0.0_pReal

allocate(constitutive_nonlocal_dDipMinScrewPerSlipSystem(maxTotalNslip, maxNinstance))
constitutive_nonlocal_dDipMinScrewPerSlipSystem = 0.0_pReal

allocate(constitutive_nonlocal_forestProjectionEdge(maxTotalNslip, maxTotalNslip, maxNinstance))
constitutive_nonlocal_forestProjectionEdge = 0.0_pReal

allocate(constitutive_nonlocal_forestProjectionScrew(maxTotalNslip, maxTotalNslip, maxNinstance))
constitutive_nonlocal_forestProjectionScrew = 0.0_pReal

allocate(constitutive_nonlocal_interactionMatrixSlipSlip(maxTotalNslip, maxTotalNslip, maxNinstance))
constitutive_nonlocal_interactionMatrixSlipSlip = 0.0_pReal


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
            'rho_edge', &
            'rho_screw', &
            'excess_rho', &
            'excess_rho_edge', &
            'excess_rho_screw', &
            'rho_forest', &
            'rho_dip', &
            'rho_edge_dip', &
            'rho_screw_dip', &
            'shearrate', &
            'resolvedstress', &
            'resistance')
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

  constitutive_nonlocal_Gmod(i) = constitutive_nonlocal_C44(i)
  constitutive_nonlocal_nu(i) = constitutive_nonlocal_C12(i) / constitutive_nonlocal_C11(i)
  
  
!*** burgers vector, dislocation velocity prefactor, mean free path prefactor and minimum dipole distance for each slip system
  
  do s = 1,constitutive_nonlocal_totalNslip(i)
    
    f = constitutive_nonlocal_slipFamily(s,i)
    
    constitutive_nonlocal_burgersPerSlipSystem(s,i) = constitutive_nonlocal_burgersPerSlipFamily(f,i)
    constitutive_nonlocal_v0PerSlipSystem(s,i) = constitutive_nonlocal_v0PerSlipFamily(f,i)
    constitutive_nonlocal_lambda0PerSlipSystem(s,i) = constitutive_nonlocal_lambda0PerSlipFamily(f,i)
    constitutive_nonlocal_dDipMinEdgePerSlipSystem(s,i) = constitutive_nonlocal_dDipMinEdgePerSlipFamily(f,i)
    constitutive_nonlocal_dDipMinScrewPerSlipSystem(s,i) = constitutive_nonlocal_dDipMinScrewPerSlipFamily(f,i)
  
  enddo
  
  
!*** calculation of forest projections for edge and screw dislocations
  
  do s1 = 1,constitutive_nonlocal_totalNslip(i)
    do s2 = 1,constitutive_nonlocal_totalNslip(i)
      
      constitutive_nonlocal_forestProjectionEdge(s1, s2, i) &
          = abs(math_mul3x3(lattice_sn(:, constitutive_nonlocal_slipSystemLattice(s1,i), myStructure), &
                            lattice_st(:, constitutive_nonlocal_slipSystemLattice(s2,i), myStructure)))                              ! forest projection of edge dislocations is the projection of (b x n) onto the slip normal of the respective splip plane
  
      
      constitutive_nonlocal_forestProjectionScrew(s1, s2, i) &
          = abs(math_mul3x3(lattice_sn(:, constitutive_nonlocal_slipSystemLattice(s1,i), myStructure), &
                            lattice_sd(:, constitutive_nonlocal_slipSystemLattice(s2,i), myStructure)))                              ! forest projection of screw dislocations is the projection of b onto the slip normal of the respective splip plane
  
  enddo; enddo
  
!*** calculation of interaction matrices
  
  do s1 = 1,constitutive_nonlocal_totalNslip(i)
    do s2 = 1,constitutive_nonlocal_totalNslip(i)
      
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
                              rhoEdgePos, &                   ! positive edge dislocation density
                              rhoEdgeNeg, &                   ! negative edge dislocation density
                              rhoScrewPos, &                  ! positive screw dislocation density
                              rhoScrewNeg, &                  ! negative screw dislocation density
                              rhoEdgeDip, &                   ! edge dipole dislocation density
                              rhoScrewDip, &                  ! screw dipole dislocation density
                              rhoForest, &                    ! forest dislocation density
                              tauSlipThreshold                ! threshold shear stress for slip
integer(pInt)                 ns, &                           ! short notation for total number of active slip systems 
                              f, &                            ! index of lattice family
                              s0, &
                              s1, &
                              s                               ! index of slip system


constitutive_nonlocal_stateInit = 0.0_pReal
ns = constitutive_nonlocal_totalNslip(myInstance)


!*** set the basic state variables

s1 = 0_pInt
do f = 1,lattice_maxNslipFamily

  s0 = s1 + 1_pInt
  s1 = s0 + constitutive_nonlocal_Nslip(f,myInstance) - 1_pInt
  
  do s = s0,s1
    rhoEdgePos(s) = constitutive_nonlocal_rhoEdgePos0(f, myInstance)
    rhoEdgeNeg(s) = constitutive_nonlocal_rhoEdgeNeg0(f, myInstance)
    rhoScrewPos(s) = constitutive_nonlocal_rhoScrewPos0(f, myInstance)
    rhoScrewNeg(s) = constitutive_nonlocal_rhoScrewNeg0(f, myInstance)
    rhoEdgeDip(s) = constitutive_nonlocal_rhoEdgeDip0(f, myInstance)
    rhoScrewDip(s) = constitutive_nonlocal_rhoScrewDip0(f, myInstance)
  enddo 
enddo


!*** calculate the dependent state variables

! forest dislocation density
forall (s = 1:ns) &
  rhoForest(s) &
      = dot_product( (rhoEdgePos + rhoEdgeNeg + rhoEdgeDip), constitutive_nonlocal_forestProjectionEdge(1:ns, s, myInstance) ) & 
      + dot_product( (rhoScrewPos + rhoScrewNeg + rhoScrewDip), constitutive_nonlocal_forestProjectionScrew(1:ns, s, myInstance) )      ! calculation of forest dislocation density as projection of screw and edge dislocations


! threshold shear stress for dislocation slip 
forall (s = 1:ns) &
  tauSlipThreshold(s) =   constitutive_nonlocal_Gmod(myInstance) & 
                        * constitutive_nonlocal_burgersPerSlipSystem(s, myInstance) &
                        * sqrt( dot_product( (rhoEdgePos + rhoEdgeNeg + rhoScrewPos + rhoScrewNeg), &
                                             constitutive_nonlocal_interactionMatrixSlipSlip(1:ns, s, myInstance) ) )


!*** put everything together and in right order

constitutive_nonlocal_stateInit(     1:  ns) = rhoEdgePos
constitutive_nonlocal_stateInit(  ns+1:2*ns) = rhoEdgeNeg
constitutive_nonlocal_stateInit(2*ns+1:3*ns) = rhoScrewPos
constitutive_nonlocal_stateInit(3*ns+1:4*ns) = rhoScrewNeg
constitutive_nonlocal_stateInit(4*ns+1:5*ns) = rhoEdgeDip
constitutive_nonlocal_stateInit(5*ns+1:6*ns) = rhoScrewDip
constitutive_nonlocal_stateInit(6*ns+1:7*ns) = rhoForest
constitutive_nonlocal_stateInit(7*ns+1:8*ns) = tauSlipThreshold

endfunction



!*********************************************************************
!* relevant microstructural state                                    *
!*********************************************************************
pure function constitutive_nonlocal_relevantState(myInstance)

use prec,     only: pReal, &
                    pInt
implicit none

!*** input variables
integer(pInt), intent(in) ::  myInstance                      ! number specifying the current instance of the constitution

!*** output variables
real(pReal), dimension(constitutive_nonlocal_sizeState(myInstance)) :: &
                              constitutive_nonlocal_relevantState ! relevant state values for the current instance of this constitution

!*** local variables

constitutive_nonlocal_relevantState = constitutive_nonlocal_relevantRho(myInstance)

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
subroutine constitutive_nonlocal_microstructure(Temperature, Fp, state, g, ip, el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use math,     only: math_Plain3333to99, &
                    math_Mandel33to6, &
                    math_Mandel6to33, &
                    math_mul33x33, &
                    math_mul3x3, &
                    math_mul33x3, &
                    pi
use debug,    only: debugger
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_element, &
                    FE_NipNeighbors, &
                    mesh_ipNeighborhood, &
                    mesh_ipVolume, &
                    mesh_ipCenterOfGravity
use material, only: homogenization_maxNgrains, &
                    material_phase, &
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
                                  Fp                          ! plastic deformation gradient

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
                                  n, &                        ! index of my current neighbor
                                  s, &                        ! index of my current slip system
                                  sLattice                    ! index of my current slip system as specified by lattice
real(pReal)                       gb, &                       ! short notation for G*b/2/pi
                                  x, &                        ! coordinate in direction of lvec
                                  y, &                        ! coordinate in direction of bvec
                                  z                           ! coordinate in direction of nvec
real(pReal), dimension(3) ::      connectingVector            ! vector connecting the centers of gravity of me and my neigbor
real(pReal), dimension(6) ::      Tdislocation_v              ! dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress in Mandel notation
real(pReal), dimension(3,3) ::    transform, &                ! orthogonal transformation matrix from slip coordinate system with e1=bxn, e2=b, e3=n to lattice coordinate system 
                                  sigma                       ! Tdislocation resulting from the excess dislocation density of a single slip system and a single neighbor calculated in the coordinate system of the slip system
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                  rhoEdgePos, &               ! positive edge dislocation density
                                  rhoEdgeNeg, &               ! negative edge dislocation density
                                  rhoScrewPos, &              ! positive screw dislocation density
                                  rhoScrewNeg, &              ! negative screw dislocation density
                                  rhoEdgeDip, &               ! edge dipole dislocation density
                                  rhoScrewDip, &              ! screw dipole dislocation density
                                  rhoForest, &                ! forest dislocation density
                                  tauSlipThreshold, &         ! threshold shear stress
                                  neighboring_rhoEdgePos, &   ! positive edge dislocation density of my neighbor
                                  neighboring_rhoEdgeNeg, &   ! negative edge dislocation density of my neighbor
                                  neighboring_rhoScrewPos, &  ! positive screw dislocation density of my neighbor
                                  neighboring_rhoScrewNeg, &  ! negative screw dislocation density of my neighbor
                                  neighboring_rhoEdgeExcess, &! edge excess dislocation density of my neighbor
                                  neighboring_rhoScrewExcess,&! screw excess dislocation density of my neighbor
                                  neighboring_Nedge, &        ! total number of edge excess dislocations in my neighbor
                                  neighboring_Nscrew


myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance)
ns = constitutive_nonlocal_totalNslip(myInstance)


!**********************************************************************
!*** get basic states

rhoEdgePos =  state(g,ip,el)%p(     1:  ns)
rhoEdgeNeg =  state(g,ip,el)%p(  ns+1:2*ns)
rhoScrewPos = state(g,ip,el)%p(2*ns+1:3*ns)
rhoScrewNeg = state(g,ip,el)%p(3*ns+1:4*ns)
rhoEdgeDip =  state(g,ip,el)%p(4*ns+1:5*ns)
rhoScrewDip = state(g,ip,el)%p(5*ns+1:6*ns)


!**********************************************************************
!*** calculate dependent states

!*** calculate the forest dislocation density

forall (s = 1:ns) &
  rhoForest(s) &
      = dot_product( (rhoEdgePos + rhoEdgeNeg + rhoEdgeDip), constitutive_nonlocal_forestProjectionEdge(1:ns, s, myInstance) ) & 
      + dot_product( (rhoScrewPos + rhoScrewNeg + rhoScrewDip), constitutive_nonlocal_forestProjectionScrew(1:ns, s, myInstance) )  ! calculation of forest dislocation density as projection of screw and edge dislocations
! if (debugger) write(6,'(a23,3(i3,x),/,12(e10.3,x),/)') 'forest dislocation density at ',g,ip,el, rhoForest

!*** calculate the threshold shear stress for dislocation slip 

forall (s = 1:ns) &
  tauSlipThreshold(s) =   constitutive_nonlocal_Gmod(myInstance) & 
                        * constitutive_nonlocal_burgersPerSlipSystem(s, myInstance) &
                        * sqrt( dot_product( (rhoEdgePos + rhoEdgeNeg + rhoScrewPos + rhoScrewNeg), &
                                             constitutive_nonlocal_interactionMatrixSlipSlip(1:ns, s, myInstance) ) )
! if (debugger) write(6,'(a26,3(i3,x),/,12(f10.5,x),/)') 'tauSlipThreshold / MPa at ',g,ip,el, tauSlipThreshold/1e6


!*** calculate the dislocation stress of the neighboring excess dislocation densities

Tdislocation_v = 0.0_pReal


! loop through my neighbors (if existent!)

do n = 1,FE_NipNeighbors(mesh_element(2,el))

  neighboring_el = mesh_ipNeighborhood(1,n,ip,el)
  neighboring_ip = mesh_ipNeighborhood(2,n,ip,el)
  
  if ( neighboring_el == 0 .or. neighboring_ip == 0 ) cycle

  
  ! calculate the connecting vector between me and my neighbor and his excess dislocation density
  
  connectingVector = math_mul33x3( Fp(:,:,g,neighboring_ip,neighboring_el), &
                                  (mesh_ipCenterOfGravity(:,ip,el) - mesh_ipCenterOfGravity(:,neighboring_ip,neighboring_el)) )
  
  neighboring_rhoEdgePos  = state(1, neighboring_ip, neighboring_el)%p(     1:  ns)
  neighboring_rhoEdgeNeg  = state(1, neighboring_ip, neighboring_el)%p(  ns+1:2*ns)
  neighboring_rhoScrewPos = state(1, neighboring_ip, neighboring_el)%p(2*ns+1:3*ns)
  neighboring_rhoScrewNeg = state(1, neighboring_ip, neighboring_el)%p(3*ns+1:4*ns)
  
  neighboring_rhoEdgeExcess = neighboring_rhoEdgePos - neighboring_rhoEdgeNeg
  neighboring_rhoScrewExcess = neighboring_rhoScrewPos - neighboring_rhoScrewNeg
  
  neighboring_Nedge = neighboring_rhoEdgeExcess * mesh_ipVolume(neighboring_ip, neighboring_el) ** (2.0_pReal/3.0_pReal)
  neighboring_Nscrew = neighboring_rhoScrewExcess * mesh_ipVolume(neighboring_ip, neighboring_el) ** (2.0_pReal/3.0_pReal)
  
  ! loop over slip systems and get their slip directions, slip normals, and sd x sn
  
  do s = 1,ns

    transform = reshape( (/lattice_st(:, constitutive_nonlocal_slipSystemLattice(s,myInstance), myStructure), &
                           lattice_sd(:, constitutive_nonlocal_slipSystemLattice(s,myInstance), myStructure), &
                           lattice_sn(:, constitutive_nonlocal_slipSystemLattice(s,myInstance), myStructure)/), (/3,3/) )
                           
    
    ! coordinate transformation of p from the lattice coordinate system to the slip coordinate system
    x = math_mul3x3(connectingVector, transform(:,1))
    y = math_mul3x3(connectingVector, transform(:,2))
    z = math_mul3x3(connectingVector, transform(:,3))
    
    
    ! calculate the back stress in the slip coordinate system for this slip system
    gb = constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgersPerSlipSystem(s,myInstance) / (2.0_pReal*pi)
    
    sigma(2,2) = - gb * neighboring_Nedge(s) / (1.0_pReal-constitutive_nonlocal_nu(myInstance)) &
                      * z * (3.0_pReal*y**2.0_pReal + z**2.0_pReal) / (y**2.0_pReal + z**2.0_pReal)**2.0_pReal
    
    sigma(3,3) =   gb * neighboring_Nedge(s) / (1.0_pReal-constitutive_nonlocal_nu(myInstance)) &
                      * z * (y**2.0_pReal - z**2.0_pReal) / (y**2.0_pReal + z**2.0_pReal)**2.0_pReal
    
    sigma(1,1) = constitutive_nonlocal_nu(myInstance) * (sigma(2,2) + sigma(3,3))
    
    sigma(1,2) = gb * neighboring_Nscrew(s) * z / (x**2.0_pReal + z**2.0_pReal)
    
    sigma(2,3) = gb * (   neighboring_Nedge(s) / (1.0_pReal-constitutive_nonlocal_nu(myInstance)) &
                             * (y**2.0_pReal - z**2.0_pReal) / (y**2.0_pReal + z**2.0_pReal) &
                        - neighboring_Nscrew(s) * x / (x**2.0_pReal + z**2.0_pReal) )
    
    sigma(2,1) = sigma(1,2)
    sigma(3,2) = sigma(2,3)
    sigma(1,3) = 0.0_pReal
    sigma(3,1) = 0.0_pReal
    
    ! coordinate transformation from the slip coordinate system to the lattice coordinate system
    Tdislocation_v = Tdislocation_v + math_Mandel33to6( math_mul33x33(transform, math_mul33x33(sigma, transpose(transform)) ) )
    ! if (debugger) write(6,'(a15,3(i3,x),/,3(3(f12.3,x)/))') 'sigma / MPa at ',g,ip,el, sigma/1e6
    ! if (debugger) write(6,'(a15,3(i3,x),/,3(3(f12.3,x)/))') 'Tdislocation / MPa at ',g,ip,el, math_Mandel6to33(Tdislocation_v/1e6)
  enddo
enddo


!**********************************************************************
!*** set dependent states

state(g,ip,el)%p(6*ns+1:7*ns) = rhoForest
state(g,ip,el)%p(7*ns+1:8*ns) = tauSlipThreshold
state(g,ip,el)%p(8*ns+1:8*ns+6) = Tdislocation_v

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
use debug,    only: debugger
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
                                            sLattice                    ! index of my current slip system as specified by lattice
real(pReal), dimension(6) ::                Tdislocation_v              ! dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress
real(pReal), dimension(3,3,3,3) ::          dLp_dTstar3333              ! derivative of Lp with respect to Tstar (3x3x3x3 matrix)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            rho                         ! dislocation densities
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            rhoForest, &                ! forest dislocation density
                                            tauSlipThreshold, &         ! threshold shear stress
                                            tauSlip, &                  ! resolved shear stress
                                            gdotSlip, &                 ! shear rate
                                            dgdot_dtauSlip, &           ! derivative of the shear rate with respect to the shear stress
                                            v                           ! dislocation velocity


!*** initialize local variables

v = 0.0_pReal
tauSlip = 0.0_pReal
gdotSlip = 0.0_pReal
Lp = 0.0_pReal
dLp_dTstar3333 = 0.0_pReal

myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)

!*** shortcut to state variables 

forall (t = 1:4) rho(:,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
rhoForest        = state(g,ip,el)%p(6*ns+1:7*ns)
tauSlipThreshold = state(g,ip,el)%p(7*ns+1:8*ns)
Tdislocation_v   = state(g,ip,el)%p(8*ns+1:8*ns+6)
! if (debugger) write(6,'(a20,3(i3,x),/,3(3(f12.3,x)/))') 'Tdislocation / MPa at ', g,ip,el, math_Mandel6to33(Tdislocation_v/1e6)
! if (debugger) write(6,'(a15,3(i3,x),/,3(3(f12.3,x)/))') 'Tstar / MPa at ',g,ip,el, math_Mandel6to33(Tstar_v/1e6)


!*** calculation of resolved stress

forall (s =1:ns) & 
  tauSlip(s) = math_mul6x6( Tstar_v + Tdislocation_v, &
                            lattice_Sslip_v(:,constitutive_nonlocal_slipSystemLattice(s,myInstance),myStructure) )


!*** Calculation of gdot and its tangent

v = constitutive_nonlocal_v0PerSlipSystem(:,myInstance) &
    * exp( - ( tauSlipThreshold - abs(tauSlip) ) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance)**2.0_pReal &
             / ( kB * Temperature * sqrt(rhoForest) ) ) &
    * sign(1.0_pReal,tauSlip)

gdotSlip =  sum(rho,2) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) * v    

dgdot_dtauSlip = gdotSlip * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance)**2.0_pReal &
                          / ( kB * Temperature * sqrt(rhoForest) )


!*** Calculation of Lp and its tangent

do s = 1,ns

  sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
  
  Lp = Lp + gdotSlip(s) * lattice_Sslip(:,:,sLattice,myStructure)
  ! if (debugger) write(6,'(a4,i2,a3,/,3(3(f15.7)/))') 'dLp(',s,'): ',gdotSlip(s) * lattice_Sslip(:,:,sLattice,myStructure)
  
  forall (i=1:3,j=1:3,k=1:3,l=1:3) &
    dLp_dTstar3333(i,j,k,l) = dLp_dTstar3333(i,j,k,l) + dgdot_dtauSlip(s) * lattice_Sslip(i,j, sLattice,myStructure) &
                                                                          * lattice_Sslip(k,l, sLattice,myStructure) 
enddo

dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)

! if (debugger) then 
 ! !$OMP CRITICAL (write2out)
   ! write(6,*) '::: LpandItsTangent',g,ip,el
   ! write(6,*)
   ! write(6,'(a,/,12(f12.5,x))') 'gdot/1e-3',gdotSlip*1e3_pReal
   ! write(6,*)
   ! write(6,'(a,/,3(3(f12.7,x)/))') 'Lp',Lp
   ! write(6,*)
 ! !$OMPEND CRITICAL (write2out)
! endif

endsubroutine



!*********************************************************************
!* rate of change of microstructure                                  *
!*********************************************************************
subroutine constitutive_nonlocal_dotState(dotState, Tstar_v, subTstar0_v, Fp, invFp, Temperature, subdt, state, subState0, g,ip,el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use debug,    only: debugger
use math,     only: math_norm3, &
                    math_mul6x6, &
                    math_mul3x3, &
                    math_mul33x3, &
                    math_transpose3x3, &
                    pi
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_element, &
                    FE_NipNeighbors, &
                    mesh_ipNeighborhood, &
                    mesh_ipVolume, &
                    mesh_ipArea, &
                    mesh_ipAreaNormal
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_constitutionInstance
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
                                            subdt                     ! substepped crystallite time increment
real(pReal), dimension(6), intent(in) ::    Tstar_v, &                ! current 2nd Piola-Kirchhoff stress in Mandel notation
                                            subTstar0_v               ! 2nd Piola-Kirchhoff stress in Mandel notation at start of crystallite increment
real(pReal), dimension(3,3), intent(in) ::  Fp, &                     ! plastic deformation gradient
                                            invFp                     ! inverse of plastic deformation gradient
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            state, &                  ! current microstructural state
                                            subState0                 ! microstructural state at start of crystallite increment

!*** input/output variables
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(inout) :: &
                                            dotState                  ! evolution of state variables / microstructure
 
!*** output variables
 
!*** local variables
integer(pInt)                               myInstance, &             ! current instance of this constitution
                                            myStructure, &            ! current lattice structure
                                            ns, &                     ! short notation for the total number of active slip systems
                                            neighboring_el, &         ! element number of my neighbor
                                            neighboring_ip, &         ! integration point of my neighbor
                                            c, &                      ! character of dislocation
                                            n, &                      ! index of my current neighbor
                                            t, &                      ! type of dislocation
                                            s                         ! index of my current slip system
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            rho, &                    ! dislocation densities (positive/negative screw and edge without dipoles)
                                            rhoDot, &                 ! rate of change of dislocation densities
                                            gdot, &                   ! shear rates
                                            lineLength                ! dislocation line length leaving the current interface
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            rhoForest, &              ! forest dislocation density
                                            tauSlipThreshold, &       ! threshold shear stress
                                            tauSlip, &                ! current resolved shear stress
                                            subTauSlip0, &            ! resolved shear stress at start of crystallite increment
                                            v, &                      ! dislocation velocity
                                            invLambda, &              ! inverse of mean free path for dislocations
                                            vClimb                    ! climb velocity of edge dipoles
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                            rhoDip, &                 ! dipole dislocation densities (screw and edge dipoles)
                                            rhoDipDot, &              ! rate of change of dipole dislocation densities
                                            rhoDotTransfer, &         ! dislocation density rate that is transferred from single dislocation to dipole dislocation
                                            dDipMin, &                ! minimum stable dipole distance for edges and screws
                                            dDipMax, &                ! current maximum stable dipole distance for edges and screws
                                            dDipMax0, &               ! maximum stable dipole distance for edges and screws at start of crystallite increment
                                            dDipMaxDot                ! rate of change of the maximum stable dipole distance for edges and screws
real(pReal), dimension(3,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            m                         ! direction of dislocation motion
real(pReal), dimension(6) ::                Tdislocation_v, &         ! current dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress
                                            subTdislocation0_v        ! dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress at start of crystallite increment
real(pReal), dimension(3) ::                surfaceNormal             ! surface normal of the current interface
real(pReal)                                 norm_surfaceNormal, &     ! euclidic norm of the surface normal
                                            area, &                   ! area of the current interface
                                            D                         ! self diffusion
                                            

myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)

tauSlip = 0.0_pReal
subTauSlip0 = 0.0_pReal
v = 0.0_pReal
gdot = 0.0_pReal
dDipMin = 0.0_pReal
dDipMax = 0.0_pReal
dDipMax0 = 0.0_pReal
dDipMaxDot = 0.0_pReal
rhoDot = 0.0_pReal
rhoDipDot = 0.0_pReal
rhoDotTransfer = 0.0_pReal

!*** shortcut to state variables 

forall (t = 1:4) rho(:,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (c = 1:2) rhoDip(:,c) = state(g,ip,el)%p((3+c)*ns+1:(4+c)*ns)
rhoForest = state(g,ip,el)%p(6*ns+1:7*ns)
tauSlipThreshold = state(g,ip,el)%p(7*ns+1:8*ns)
Tdislocation_v = state(g,ip,el)%p(8*ns+1:8*ns+6)
subTdislocation0_v = subState0(g,ip,el)%p(8*ns+1:8*ns+6)


!****************************************************************************
!*** Calculate shear rate

do s = 1,ns   ! loop over slip systems

  tauSlip(s) = math_mul6x6( Tstar_v + Tdislocation_v, &
                         lattice_Sslip_v(:,constitutive_nonlocal_slipSystemLattice(s,myInstance),myStructure) )
  
  subTauSlip0(s) = math_mul6x6( subTstar0_v + subTdislocation0_v, &
                         lattice_Sslip_v(:,constitutive_nonlocal_slipSystemLattice(s,myInstance),myStructure) )
enddo

v = constitutive_nonlocal_v0PerSlipSystem(:,myInstance) &
    * exp( - ( tauSlipThreshold - abs(tauSlip) ) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance)**2.0_pReal &
             / ( kB * Temperature * sqrt(rhoForest) ) ) &
    * sign(1.0_pReal,tauSlip)
    
forall (t = 1:4) &
  gdot(:,t) = rho(:,t) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) * v


!****************************************************************************
!*** calculate dislocation multiplication

invLambda = sqrt(rhoForest) / constitutive_nonlocal_lambda0PerSlipSystem(:,myInstance)

rhoDot = rhoDot + spread(0.25_pReal * sum(abs(gdot),2) * invLambda / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance), 2, 4)
if (debugger) then 
  write(6,*) '::: constitutive_nonlocal_dotState at ',g,ip,el
  write(6,*)
  write(6,'(a,/,12(f12.5,x),/)') 'tauSlip / MPa', tauSlip/1e6_pReal
  write(6,'(a,/,12(f12.5,x),/)') 'tauSlipThreshold / MPa', tauSlipThreshold/1e6_pReal
  write(6,'(a,/,12(e10.3,x),/)') 'v', v
  write(6,'(a,/,4(12(f12.5,x),/))') 'gdot / 1e-3', gdot*1e3_pReal
  write(6,'(a,/,(12(f12.5,x),/))') 'gdot total/ 1e-3', sum(gdot,2)*1e3_pReal
  write(6,'(a,/,6(12(e12.5,x),/))') 'dislocation multiplication', &
        spread(0.25_pReal * sum(abs(gdot),2) * invLambda / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance), 2, 4)*subdt, &
        0.0_pReal*rhoDotTransfer
endif


!****************************************************************************
!*** calculate dislocation fluxes

m(:,:,1) =  lattice_sd(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)
m(:,:,2) = -lattice_sd(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)
m(:,:,3) =  lattice_st(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)
m(:,:,4) = -lattice_st(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)

do n = 1,FE_NipNeighbors(mesh_element(2,el))                                                        ! loop through my neighbors

  neighboring_el = mesh_ipNeighborhood(1,n,ip,el)
  neighboring_ip = mesh_ipNeighborhood(2,n,ip,el)
   
  ! calculate the area and the surface normal of the interface
  surfaceNormal = math_mul33x3(math_transpose3x3(invFp), mesh_ipAreaNormal(:,n,ip,el))
  norm_surfaceNormal = math_norm3(surfaceNormal)
  surfaceNormal = surfaceNormal / norm_surfaceNormal
  area = mesh_ipArea(n,ip,el) / norm_surfaceNormal

  lineLength = 0.0_pReal
  
  do s = 1,ns                                                                                       ! loop over slip systems
  
    do t = 1,4                                                                                      ! loop over dislocation types
      
      if ( sign(1.0_pReal,math_mul3x3(m(:,s,t),surfaceNormal)) == sign(1.0_pReal,gdot(s,t)) ) then
        
        lineLength(s,t) = gdot(s,t) / constitutive_nonlocal_burgersPerSlipSystem(s,myInstance) &
                                    * math_mul3x3(m(:,s,t),surfaceNormal) * area                    ! dislocation line length that leaves this interface per second
        
        rhoDot(s,t) = rhoDot(s,t) - lineLength(s,t) / mesh_ipVolume(ip,el)                          ! subtract dislocation density rate (= line length over volume) that leaves through an interface from my dotState ...
        
        if ( neighboring_el > 0 .and. neighboring_ip > 0 ) then
!*****************************************************************************************************
!***   OMP locking for this neighbor missing
!*****************************************************************************************************
          dotState(1,neighboring_ip,neighboring_el)%p((t-1)*ns+s) = dotState(1,neighboring_ip,neighboring_el)%p((t-1)*ns+s) &
                                                                    + lineLength(s,t) / mesh_ipVolume(neighboring_ip,neighboring_el)  ! ... and add it to the neighboring dotState (if neighbor exists)
        endif
      endif
    enddo
  enddo
enddo
if (debugger) write(6,'(a,/,6(12(e12.5,x),/))') 'dislocation flux', lineLength/mesh_ipVolume(ip,el)*subdt, 0.0_pReal*rhoDotTransfer


!****************************************************************************
!*** calculate dipole formation and annihilation

!*** limits for stable dipole height and its tate of change
  
dDipMin(:,1) = constitutive_nonlocal_dDipMinEdgePerSlipSystem(:,myInstance)
dDipMin(:,2) = constitutive_nonlocal_dDipMinScrewPerSlipSystem(:,myInstance)
dDipMax(:,2) = constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                  / ( 8.0_pReal * pi * abs(tauSlip) )
dDipMax(:,1) = dDipMax(:,2) / ( 1.0_pReal - constitutive_nonlocal_nu(myInstance) )
dDipMax0(:,2) = constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                  / ( 8.0_pReal * pi * abs(subTauSlip0) )
dDipMax0(:,1) = dDipMax0(:,2) / ( 1.0_pReal - constitutive_nonlocal_nu(myInstance) )
  
dDipMaxDot(:,1) = (dDipMax(:,1) - dDipMax0(:,1)) / subdt
dDipMaxDot(:,2) = (dDipMax(:,2) - dDipMax0(:,2)) / subdt
! if (debugger) write(6,'(a,/,2(12(e12.5,x),/))') 'dDipMax:',dDipMax
! if (debugger) write(6,'(a,/,2(12(e12.5,x),/))') 'dDipMaxDot:',dDipMaxDot
  
  
!*** formation by glide

forall (c=1:2) &  
  rhoDotTransfer(:,c) = 2.0_pReal * dDipMax(:,c) / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                                  * ( rho(:,2*c-1)*abs(gdot(:,2*c)) + rho(:,2*c)*abs(gdot(:,2*c-1)) )
if (debugger) write(6,'(a,/,6(12(e12.5,x),/))') 'dipole formation by glide', &
                                                        -rhoDotTransfer*subdt,-rhoDotTransfer*subdt,2.0_pReal*rhoDotTransfer*subdt
  
rhoDot(:,(/1,3/)) = rhoDot(:,(/1,3/)) - rhoDotTransfer                                            ! subtract from positive single dislocation density of this character
rhoDot(:,(/2,4/)) = rhoDot(:,(/2,4/)) - rhoDotTransfer                                            ! subtract from negative single dislocation density of this character
rhoDipDot = rhoDipDot + 2.0_pReal * rhoDotTransfer                                                ! add twice to dipole dislocation density of this character
  

!*** athermal annihilation

forall (c=1:2) &  
  rhoDotTransfer(:,c) = - 2.0_pReal * dDipMin(:,c) / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) &
                                    * ( rho(:,2*c-1)*abs(gdot(:,2*c)) + rho(:,2*c)*abs(gdot(:,2*c-1)) )
if (debugger) write(6,'(a,/,6(12(e12.5,x),/))') 'athermal dipole annihilation', &
                                                 0.0_pReal*rhoDotTransfer,0.0_pReal*rhoDotTransfer,2.0_pReal*rhoDotTransfer*subdt

rhoDipDot = rhoDipDot + 2.0_pReal * rhoDotTransfer                                                ! add twice to dipole dislocation density of this character
  
  
!*** thermally activated annihilation

D = constitutive_nonlocal_D0(myInstance) * exp(-constitutive_nonlocal_Qsd(myInstance) / (kB * Temperature))

vClimb =  constitutive_nonlocal_atomicVolume(myInstance) * D / ( kB * Temperature ) &
          * constitutive_nonlocal_Gmod(myInstance) / ( 2.0_pReal * pi * (1.0_pReal-constitutive_nonlocal_nu(myInstance)) ) &
          * 2.0_pReal / ( dDipMax(:,1) + dDipMin(:,1) )
          
rhoDipDot(:,1) = rhoDipDot(:,1) - 4.0_pReal * rho(:,1) * vClimb / ( dDipMax(:,1) - dDipMin(:,1) )
if (debugger) write(6,'(a,/,6(12(e12.5,x),/))') 'thermally activated dipole annihilation', &
      0.0_pReal*rhoDotTransfer, 0.0_pReal*rhoDotTransfer, - 4.0_pReal * rho(:,1) * vClimb / ( dDipMax(:,1) - dDipMin(:,1) )*subdt, &
      0.0_pReal*vClimb


! !*** formation by stress decrease = increase in dDipMax

! forall (s=1:ns, dDipMaxDot(s,1) > 0.0_pReal) &  
  ! rhoDotTransfer(s,:) = 4.0_pReal * rho(s,(/1,3/)) * rho(s,(/2,4/)) * dDipMax0(s,:) * dDipMaxDot(s,:)

! if (debugger) write(6,'(a,/,6(12(e12.5,x),/))') 'dipole formation by stress decrease',& 
                                                        ! -rhoDotTransfer*subdt,-rhoDotTransfer*subdt,2.0_pReal*rhoDotTransfer*subdt

! rhoDot(:,(/1,3/)) = rhoDot(:,(/1,3/)) - rhoDotTransfer                                            ! subtract from positive single dislocation density of this character
! rhoDot(:,(/2,4/)) = rhoDot(:,(/2,4/)) - rhoDotTransfer                                            ! subtract from negative single dislocation density of this character
! rhoDipDot = rhoDipDot + 2.0_pReal * rhoDotTransfer                                                ! add twice to dipole dislocation density of this character


! !*** dipole dissociation by increased stress = decrease in dDipMax

! forall (s=1:ns, dDipMaxDot(s,1) < 0.0_pReal) &  
  ! rhoDotTransfer(s,:) = 0.5_pReal * rhoDip(s,:) * dDipMaxDot(s,:) / (dDipMax0(s,:) - dDipMin(s,:))

! if (debugger) write(6,'(a,/,6(12(e12.5,x),/))') 'dipole formation by stress decrease',& 
                                                        ! -rhoDotTransfer*subdt,-rhoDotTransfer*subdt,2.0_pReal*rhoDotTransfer*subdt
  
! rhoDot(:,(/1,3/)) = rhoDot(:,(/1,3/)) - rhoDotTransfer                                            ! subtract from positive single dislocation density of this character
! rhoDot(:,(/2,4/)) = rhoDot(:,(/2,4/)) - rhoDotTransfer                                            ! subtract from negative single dislocation density of this character
! rhoDipDot = rhoDipDot + 2.0_pReal * rhoDotTransfer                                                ! add twice to dipole dislocation density of this character


!****************************************************************************
!*** assign the rates of dislocation densities to my dotState

dotState(1,ip,el)%p(1:4*ns) = reshape(rhoDot,(/4*ns/))
dotState(1,ip,el)%p(4*ns+1:6*ns) = reshape(rhoDipDot,(/2*ns/))

if (debugger) write(6,'(a,/,4(12(e12.5,x),/))') 'deltaRho:',rhoDot*subdt
if (debugger) write(6,'(a,/,2(12(e12.5,x),/))') 'deltaRhoDip:',rhoDipDot*subdt

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
pure function constitutive_nonlocal_postResults(Tstar_v, Temperature, dt, state, g, ip, el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use math,     only: math_mul6x6
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_constitutionInstance, &
                    phase_Noutput
use lattice,  only: lattice_Sslip_v, &
                    lattice_NslipSystem
implicit none

!*** input variables
integer(pInt), intent(in) ::              g, &                ! current grain number
                                          ip, &               ! current integration point
                                          el                  ! current element number
real(pReal), intent(in) ::                dt, &               ! time increment
                                          Temperature         ! temperature
real(pReal), dimension(6), intent(in) ::  Tstar_v             ! 2nd Piola-Kirchhoff stress in Mandel notation
type(p_vec), dimension(homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems), intent(in) :: &
                                          state               ! microstructural state

!*** output variables
real(pReal), dimension(constitutive_nonlocal_sizePostResults(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                          constitutive_nonlocal_postResults

!*** local variables
integer(pInt)                             myInstance, &       ! current instance of this constitution
                                          myStructure, &      ! current lattice structure
                                          ns, &               ! short notation for the total number of active slip systems
                                          o, &                ! index of current output
                                          s, &                ! index of current slip system
                                          sLattice, &         ! index of current slip system as specified by lattice
                                          c
real(pReal)                               tau, &              ! resolved shear stress on current slip system
                                          v                   ! dislocation velocity on current slip system


myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)


c = 0_pInt
constitutive_nonlocal_postResults = 0.0_pReal

do o = 1,phase_Noutput(material_phase(g,ip,el))
  select case(constitutive_nonlocal_output(o,myInstance))
    
    case ('rho')
      constitutive_nonlocal_postResults(c+1:c+ns) =   state(g,ip,el)%p(1:ns) + state(g,ip,el)%p(ns+1:2*ns) &
                                                    + state(g,ip,el)%p(2*ns+1:3*ns) + state(g,ip,el)%p(3*ns+1:4*ns)
      c = c + ns
      
    case ('rho_edge')
      constitutive_nonlocal_postResults(c+1:c+ns) = state(g,ip,el)%p(1:ns) + state(g,ip,el)%p(ns+1:2*ns)
      c = c + ns
      
    case ('rho_screw')
      constitutive_nonlocal_postResults(c+1:c+ns) = state(g,ip,el)%p(2*ns+1:3*ns) + state(g,ip,el)%p(3*ns+1:4*ns)
      c = c + ns
      
    case ('excess_rho')
      constitutive_nonlocal_postResults(c+1:c+ns) =   state(g,ip,el)%p(1:ns) - state(g,ip,el)%p(ns+1:2*ns) &
                                                    + state(g,ip,el)%p(2*ns+1:3*ns) - state(g,ip,el)%p(3*ns+1:4*ns)
      c = c + ns
      
    case ('excess_rho_edge')
      constitutive_nonlocal_postResults(c+1:c+ns) = state(g,ip,el)%p(1:ns) - state(g,ip,el)%p(ns+1:2*ns)
      c = c + ns
      
    case ('excess_rho_screw')
      constitutive_nonlocal_postResults(c+1:c+ns) = state(g,ip,el)%p(2*ns+1:3*ns) - state(g,ip,el)%p(3*ns+1:4*ns)
      c = c + ns
      
    case ('rho_forest')
      constitutive_nonlocal_postResults(c+1:c+ns) = state(g,ip,el)%p(6*ns+1:7*ns)
      c = c + ns
    
    case ('rho_dip')
      constitutive_nonlocal_postResults(c+1:c+ns) =   state(g,ip,el)%p(4*ns+1:5*ns) + state(g,ip,el)%p(5*ns+1:6*ns)
      c = c + ns
      
    case ('rho_edge_dip')
      constitutive_nonlocal_postResults(c+1:c+ns) = state(g,ip,el)%p(4*ns+1:5*ns)
      c = c + ns
      
    case ('rho_screw_dip')
      constitutive_nonlocal_postResults(c+1:c+ns) = state(g,ip,el)%p(5*ns+1:6*ns)
      c = c + ns
    
    case ('shearrate')
      do s = 1,ns
        sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
        tau = math_mul6x6( Tstar_v + state(g,ip,el)%p(8*ns+1:8*ns+6), lattice_Sslip_v(:,sLattice,myStructure) )
        
        if (state(g,ip,el)%p(4*ns+s) > 0.0_pReal) then
          v =  constitutive_nonlocal_v0PerSlipSystem(s,myInstance) &
             * exp( - ( state(g,ip,el)%p(7*ns+s) - abs(tau) ) * constitutive_nonlocal_burgersPerSlipSystem(s,myInstance)**2.0_pReal &
                      / ( kB * Temperature * sqrt(state(g,ip,el)%p(4*ns+s)) ) ) &
             * sign(1.0_pReal,tau)
        else
          v = 0.0_pReal
        endif
        
        constitutive_nonlocal_postResults(c+s) =  (   state(g,ip,el)%p(s) + state(g,ip,el)%p(ns+s) &
                                                    + state(g,ip,el)%p(2*ns+s) + state(g,ip,el)%p(3*ns+s) ) &
                                                  * constitutive_nonlocal_burgersPerSlipSystem(s,myInstance) * v
      enddo
      c = c + ns
      
    case ('resolvedstress')
      do s = 1,ns  
        sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
        constitutive_nonlocal_postResults(c+s) = math_mul6x6( Tstar_v + state(g,ip,el)%p(8*ns+1:8*ns+6), &
                                                              lattice_Sslip_v(:,sLattice,myStructure) )
      enddo
      c = c + ns
      
    case ('resistance')
      constitutive_nonlocal_postResults(c+1:c+ns) = state(g,ip,el)%p(7*ns+1:8*ns)
      c = c + ns

 end select
enddo

endfunction
END MODULE
