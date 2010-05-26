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
                                                          constitutive_nonlocal_relevantRho                     ! dislocation density considered relevant
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
                    mesh_maxNips
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
constitutive_nonlocal_Q0 = 0.0_pReal
constitutive_nonlocal_atomicVolume = 0.0_pReal
constitutive_nonlocal_D0 = 0.0_pReal
constitutive_nonlocal_Qsd = 0.0_pReal
constitutive_nonlocal_relevantRho = 0.0_pReal
constitutive_nonlocal_nu = 0.0_pReal
constitutive_nonlocal_Cslip_66 = 0.0_pReal
constitutive_nonlocal_Cslip_3333 = 0.0_pReal

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
constitutive_nonlocal_rhoSglEdgePos0 = 0.0_pReal
constitutive_nonlocal_rhoSglEdgeNeg0 = 0.0_pReal
constitutive_nonlocal_rhoSglScrewPos0 = 0.0_pReal
constitutive_nonlocal_rhoSglScrewNeg0 = 0.0_pReal
constitutive_nonlocal_rhoDipEdge0 = 0.0_pReal
constitutive_nonlocal_rhoDipScrew0 = 0.0_pReal
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
  if (constitutive_nonlocal_atomicVolume(i) <= 0.0_pReal)                                   call IO_error(230)
  if (constitutive_nonlocal_D0(i) <= 0.0_pReal)                                             call IO_error(231)
  if (constitutive_nonlocal_Qsd(i) <= 0.0_pReal)                                            call IO_error(232)
  if (constitutive_nonlocal_relevantRho(i) <= 0.0_pReal)                                    call IO_error(233)
    
  
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

  constitutive_nonlocal_Gmod(i) = constitutive_nonlocal_C44(i)
  constitutive_nonlocal_nu(i) = 0.5_pReal * constitutive_nonlocal_C12(i) &
                                          / (constitutive_nonlocal_C12(i) + constitutive_nonlocal_C44(i))
  
  
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
                            lattice_st(:, constitutive_nonlocal_slipSystemLattice(s2,i), myStructure)))                              ! forest projection of edge dislocations is the projection of (t = b x n) onto the slip normal of the respective slip plane
      
      constitutive_nonlocal_forestProjectionScrew(s1, s2, i) &
          = abs(math_mul3x3(lattice_sn(:, constitutive_nonlocal_slipSystemLattice(s1,i), myStructure), &
                            lattice_sd(:, constitutive_nonlocal_slipSystemLattice(s2,i), myStructure)))                              ! forest projection of screw dislocations is the projection of b onto the slip normal of the respective splip plane
  
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
                    verboseDebugger, &
                    selectiveDebugger
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
                                  opposite_n, &               ! index of my opposite neighbor
                                  opposite_ip, &              ! ip of my opposite neighbor
                                  opposite_el, &              ! element index of my opposite neighbor
                                  s, &                        ! index of my current slip system
                                  t, &                        ! index of dilsocation type (e+, e-, s+, s-, used e+, used e-, used s+, used s-)
                                  sLattice                    ! index of my current slip system according to lattice order
real(pReal)                       gb, &                       ! short notation for G*b/2/pi
                                  x, &                        ! coordinate in direction of lvec
                                  y, &                        ! coordinate in direction of bvec
                                  z, &                        ! coordinate in direction of nvec
                                  detFe, &                    ! determinant of elastic deformation gradient
                                  neighboring_detFe           ! determinant of my neighboring elastic deformation gradient
real(pReal), dimension(3,6) ::    connectingVector            ! vector connecting the centers of gravity of me and my neigbor (for each neighbor)
real(pReal), dimension(6) ::      Tdislocation_v              ! dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress in Mandel notation
real(pReal), dimension(3,3) ::    lattice2slip, &             ! orthogonal transformation matrix from lattice coordinate system to slip coordinate system with e1=bxn, e2=b, e3=n
                                  neighboringSlip2myLattice, &! mapping from my neighbors slip coordinate system to my lattice coordinate system
                                  sigma, &                    ! Tdislocation resulting from the excess dislocation density of a single slip system and a single neighbor calculated in the coordinate system of the slip system
                                  F, &                        ! total deformation gradient
                                  neighboring_F, &            ! total deformation gradient of my neighbor
                                  Favg, &                     ! average total deformation gradient of me and my neighbor
                                  invFe, &                    ! inverse of elastic deformation gradient
                                  neighboring_invFe           ! inverse of my neighboring elastic deformation gradient 
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),8) :: &
                                  rhoSgl, &                   ! single dislocation density (edge+, edge-, screw+, screw-, used edge+, used edge-, used screw+, used screw-)
                                  neighboring_rhoSgl          ! single dislocation density of my neighbor (edge+, edge-, screw+, screw-, used edge+, used edge-, used screw+, used screw-)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                  rhoDip                      ! dipole dislocation density (edge, screw)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                  rhoForest, &                ! forest dislocation density
                                  tauThreshold, &             ! threshold shear stress
                                  tau, &                      ! resolved shear stress
                                  neighboring_rhoEdgeExcess, &! edge excess dislocation density of my neighbor
                                  neighboring_rhoScrewExcess,&! screw excess dislocation density of my neighbor
                                  neighboring_Nedge, &        ! total number of edge excess dislocations in my neighbor
                                  neighboring_Nscrew          ! total number of screw excess dislocations in my neighbor


myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance)
ns = constitutive_nonlocal_totalNslip(myInstance)


!**********************************************************************
!*** get basic states

forall (t = 1:8) rhoSgl(:,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (c = 1:2) rhoDip(:,c) = state(g,ip,el)%p((c+7)*ns+1:(c+8)*ns)


!**********************************************************************
!*** calculate dependent states

!*** calculate the forest dislocation density

forall (s = 1:ns) &
  rhoForest(s) =  dot_product( ( sum(abs(rhoSgl(:,(/1,2,5,6/))),2) + rhoDip(:,1) ), &
                               constitutive_nonlocal_forestProjectionEdge(s, 1:ns, myInstance) ) & 
                + dot_product( ( sum(abs(rhoSgl(:,(/3,4,7,8/))),2) + rhoDip(:,2) ), &
                               constitutive_nonlocal_forestProjectionScrew(s, 1:ns, myInstance) )      ! calculation of forest dislocation density as projection of screw and edge dislocations
! if (debugger) write(6,'(a30,3(i3,x),/,12(e10.3,x),/)') 'forest dislocation density at ',g,ip,el, rhoForest


!*** calculate the threshold shear stress for dislocation slip 

forall (s = 1:ns) &
  tauThreshold(s) =   constitutive_nonlocal_Gmod(myInstance) & 
                    * constitutive_nonlocal_burgersPerSlipSystem(s, myInstance) &
                    * sqrt( dot_product( ( sum(abs(rhoSgl),2) + sum(abs(rhoDip),2) ), &
                                         constitutive_nonlocal_interactionMatrixSlipSlip(s, 1:ns, myInstance) ) )
! if (debugger) write(6,'(a26,3(i3,x),/,12(f10.5,x),/)') 'tauSlipThreshold / MPa at ',g,ip,el, tauSlipThreshold/1e6


!*** calculate the dislocation stress of the neighboring excess dislocation densities

Tdislocation_v = 0.0_pReal
connectingVector = 0.0_pReal
F = math_mul33x33(Fe(:,:,g,ip,el), Fp(:,:,g,ip,el))
detFe = math_det3x3(Fe)
invFe = math_inv3x3(Fe)

do n = 1,FE_NipNeighbors(mesh_element(2,el))
  
  neighboring_el = mesh_ipNeighborhood(1,n,ip,el)
  neighboring_ip = mesh_ipNeighborhood(2,n,ip,el)  
  if ( neighboring_ip == 0 ) &
    cycle
  
  neighboring_F = math_mul33x33(Fe(:,:,g,neighboring_ip,neighboring_el), Fp(:,:,g,neighboring_ip,neighboring_el))
  Favg = 0.5_pReal * (F + neighboring_F)
  neighboring_detFe = math_det3x3(Fe(:,:,g,neighboring_ip,neighboring_el))
  neighboring_invFe = math_inv3x3(Fe(:,:,g,neighboring_ip,neighboring_el))
  
  connectingVector(:,n) = math_mul33x3(neighboring_invFe, math_mul33x3(Favg, &
                                  (mesh_ipCenterOfGravity(:,neighboring_ip,neighboring_el) - mesh_ipCenterOfGravity(:,ip,el)) ) )   ! calculate connection vector between me and my neighbor in its lattice configuration
  
  opposite_n = n - 1_pInt + 2_pInt*mod(n,2_pInt)
  opposite_el = mesh_ipNeighborhood(1,opposite_n,ip,el)
  opposite_ip = mesh_ipNeighborhood(2,opposite_n,ip,el)
  if ( opposite_ip == 0 ) &                                                                                                         ! if no opposite neighbor present...
    connectingVector(:,opposite_n) = -connectingVector(:,n)                                                                         ! ... assume "ghost" IP at mirrored position of this neighbor 
  
enddo


do n = 1,FE_NipNeighbors(mesh_element(2,el))                                                                                        ! loop through my neighbors

  neighboring_el = mesh_ipNeighborhood(1,n,ip,el)
  neighboring_ip = mesh_ipNeighborhood(2,n,ip,el)
  
  if ( neighboring_ip == 0 ) then                                                                                                   ! if no neighbor present...
    opposite_n = n - 1_pInt + 2_pInt*mod(n,2_pInt)
    opposite_el = mesh_ipNeighborhood(1,opposite_n,ip,el)
    opposite_ip = mesh_ipNeighborhood(2,opposite_n,ip,el)
    if ( opposite_ip == 0 ) &                                                                                                       ! if both neighbors not present...
      cycle    neighboring_el = opposite_el                                                                                         ! ... skip this element
    neighboring_ip = opposite_ip
    forall (t = 1:8)  neighboring_rhoSgl(:,t) = max(0.0_pReal, 2.0_pReal * state(g,ip,el)%p((t-1)*ns+1:t*ns) &
                                                                - state(g,opposite_ip,opposite_el)%p((t-1)*ns+1:t*ns) )             ! ... extrapolate density from opposite neighbor (but assure positive value for density)
  else
    forall (t = 1:8)  neighboring_rhoSgl(:,t) = state(g,neighboring_ip,neighboring_el)%p((t-1)*ns+1:t*ns)
  endif
  
  neighboring_rhoEdgeExcess  = sum(abs(neighboring_rhoSgl(:,(/1,5/))),2) - sum(abs(neighboring_rhoSgl(:,(/2,6/))),2)
  neighboring_rhoScrewExcess = sum(abs(neighboring_rhoSgl(:,(/3,7/))),2) - sum(abs(neighboring_rhoSgl(:,(/4,8/))),2)  
  
  neighboring_Nedge = neighboring_rhoEdgeExcess * mesh_ipVolume(neighboring_ip,neighboring_el) ** (2.0_pReal/3.0_pReal)
  neighboring_Nscrew = neighboring_rhoScrewExcess * mesh_ipVolume(neighboring_ip,neighboring_el) ** (2.0_pReal/3.0_pReal)
    
  do s = 1,ns

    lattice2slip = reshape( (/lattice_st(:, constitutive_nonlocal_slipSystemLattice(s,myInstance), myStructure), &
                              lattice_sd(:, constitutive_nonlocal_slipSystemLattice(s,myInstance), myStructure), &
                              lattice_sn(:, constitutive_nonlocal_slipSystemLattice(s,myInstance), myStructure)/), (/3,3/) )
                           
    x = math_mul3x3(lattice2slip(1,:), -connectingVector(:,n))                                                                      ! coordinate transformation of connecting vector from the lattice coordinate system to the slip coordinate system
    y = math_mul3x3(lattice2slip(2,:), -connectingVector(:,n))
    z = math_mul3x3(lattice2slip(3,:), -connectingVector(:,n))
        
    gb = constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgersPerSlipSystem(s,myInstance) / (2.0_pReal*pi)
    
    sigma(2,2) = - gb * neighboring_Nedge(s) / (1.0_pReal-constitutive_nonlocal_nu(myInstance)) &
                      * z * (3.0_pReal*y**2.0_pReal + z**2.0_pReal) / (y**2.0_pReal + z**2.0_pReal)**2.0_pReal
    
    sigma(3,3) =   gb * neighboring_Nedge(s) / (1.0_pReal-constitutive_nonlocal_nu(myInstance)) &
                      * z * (y**2.0_pReal - z**2.0_pReal) / (y**2.0_pReal + z**2.0_pReal)**2.0_pReal
    
    sigma(1,1) = constitutive_nonlocal_nu(myInstance) * (sigma(2,2) + sigma(3,3))
    
    sigma(1,2) = - gb * neighboring_Nscrew(s) * z / (x**2.0_pReal + z**2.0_pReal)
    
    sigma(2,3) = gb * (   neighboring_Nedge(s) / (1.0_pReal-constitutive_nonlocal_nu(myInstance)) &
                             * y * (y**2.0_pReal - z**2.0_pReal) / (y**2.0_pReal + z**2.0_pReal)**2.0_pReal &
                        + neighboring_Nscrew(s) * x / (x**2.0_pReal + z**2.0_pReal) )
    
    sigma(2,1) = sigma(1,2)
    sigma(3,2) = sigma(2,3)
    sigma(1,3) = 0.0_pReal
    sigma(3,1) = 0.0_pReal
    
    neighboringSlip2myLattice = math_mul33x33(invFe,math_mul33x33(Fe(:,:,g,neighboring_ip,neighboring_el),transpose(lattice2slip))) ! coordinate transformation from the slip coordinate system to the lattice coordinate system
    Tdislocation_v = Tdislocation_v + math_Mandel33to6( detFe / math_det3x3(Fe(:,:,g,neighboring_ip,neighboring_el)) &
                           * math_mul33x33(neighboringSlip2myLattice, math_mul33x33(sigma, transpose(neighboringSlip2myLattice))) )
  enddo
enddo

!**********************************************************************
!*** set dependent states

state(g,ip,el)%p(10*ns+1:11*ns) = rhoForest
state(g,ip,el)%p(11*ns+1:12*ns) = tauThreshold
state(g,ip,el)%p(12*ns+1:12*ns+6) = Tdislocation_v


!**********************************************************************
!*** calculate the dislocation velocity

call constitutive_nonlocal_kinetics(Tstar_v, Temperature, state, g, ip, el)

endsubroutine



!*********************************************************************
!* calculates kinetics                                               *
!*********************************************************************
subroutine constitutive_nonlocal_kinetics(Tstar_v, Temperature, state, g, ip, el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use math,     only: math_mul6x6, &
                    math_Mandel6to33
use debug,    only: debugger, &
                    selectiveDebugger, &
                    verboseDebugger
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

!*** local variables
integer(pInt)                               myInstance, &               ! current instance of this constitution
                                            myStructure, &              ! current lattice structure
                                            ns, &                       ! short notation for the total number of active slip systems
                                            t, &                        ! dislocation type
                                            s                           ! index of my current slip system
real(pReal), dimension(6) ::                Tdislocation_v              ! dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            tauThreshold, &             ! threshold shear stress
                                            tau                         ! resolved shear stress


myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)

tauThreshold = state(g,ip,el)%p(11*ns+1:12*ns)
Tdislocation_v = state(g,ip,el)%p(12*ns+1:12*ns+6)

tau = 0.0_pReal
constitutive_nonlocal_v(:,:,g,ip,el) = 0.0_pReal

do s = 1,ns
  if ((tauThreshold(s) > 0.0_pReal) .and. (Temperature > 0.0_pReal)) then
  
    tau(s) = math_mul6x6( Tstar_v + Tdislocation_v, &
                        lattice_Sslip_v(:,constitutive_nonlocal_slipSystemLattice(s,myInstance),myStructure) )

    constitutive_nonlocal_v(s,:,g,ip,el) = constitutive_nonlocal_v0PerSlipSystem(s,myInstance) &
            * exp( - constitutive_nonlocal_Q0(myInstance) / ( kB * Temperature) * (1.0_pReal - (abs(tau(s))/tauThreshold(s)) ) ) &
            * sign(1.0_pReal,tau(s))
  
  endif
enddo

if (verboseDebugger .and. selectiveDebugger) then 
  !$OMP CRITICAL (write2out)
    write(6,*) '::: kinetics',g,ip,el
    write(6,*)
!    write(6,'(a,/,3(3(f12.3,x)/))') 'Tdislocation / MPa', math_Mandel6to33(Tdislocation_v/1e6)
!    write(6,'(a,/,3(3(f12.3,x)/))') 'Tstar / MPa', math_Mandel6to33(Tstar_v/1e6)
    write(6,'(a,/,12(f12.5,x),/)') 'tau / MPa', tau/1e6_pReal
    write(6,'(a,/,12(f12.5,x),/)') 'tauThreshold / MPa', tauThreshold/1e6_pReal
    write(6,'(a,/,4(12(f12.5,x),/))') 'v / 1e-3m/s', constitutive_nonlocal_v(:,:,g,ip,el)*1e3
  !$OMPEND CRITICAL (write2out)
endif

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
                    selectiveDebugger, &
                    verboseDebugger
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
                                            dgdotTotal_dtau             ! derivative of the shear rate with respect to the shear stress


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
forall (s = 1:ns, t = 5:8, rhoSgl(s,t) * constitutive_nonlocal_v(s,t-4,g,ip,el) < 0.0_pReal) &                                        ! contribution of used rho for changing sign of v
  rhoSgl(s,t-4) = rhoSgl(s,t-4) + abs(rhoSgl(s,t))

tauThreshold = state(g,ip,el)%p(11*ns+1:12*ns)

call constitutive_nonlocal_kinetics(Tstar_v, Temperature, state, g, ip, el)                                                         ! update dislocation velocity

!*** Calculation of gdot and its tangent

forall (t = 1:4 ) &
  gdot(:,t) =  rhoSgl(:,t) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) * constitutive_nonlocal_v(:,t,g,ip,el)
gdotTotal = sum(gdot,2)

dgdotTotal_dtau = abs(gdotTotal) * constitutive_nonlocal_Q0(myInstance) / ( kB * Temperature * tauThreshold )

!*** Calculation of Lp and its tangent

do s = 1,ns
  sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
  
  Lp = Lp + gdotTotal(s) * lattice_Sslip(:,:,sLattice,myStructure)
  
  forall (i=1:3,j=1:3,k=1:3,l=1:3) &
    dLp_dTstar3333(i,j,k,l) = dLp_dTstar3333(i,j,k,l) + dgdotTotal_dtau(s) * lattice_Sslip(i,j, sLattice,myStructure) &
                                                                           * lattice_Sslip(k,l, sLattice,myStructure) 
enddo

dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)

!if (verboseDebugger .and. selectiveDebugger) then 
!  !$OMP CRITICAL (write2out)
!    write(6,*) '::: LpandItsTangent',g,ip,el
!    write(6,*)
!    ! write(6,'(a,/,12(f12.5,x),/)') 'v', constitutive_nonlocal_v(:,t,g,ip,el)
!    write(6,'(a,/,12(f12.5,x),/)') 'gdot /1e-3',gdot*1e3_pReal
!    write(6,'(a,/,12(f12.5,x),/)') 'gdot total /1e-3',gdotTotal*1e3_pReal
!    write(6,'(a,/,3(3(f12.7,x)/))') 'Lp',Lp
!    ! call flush(6)
!  !$OMPEND CRITICAL (write2out)
!endif

endsubroutine



!*********************************************************************
!* rate of change of microstructure                                  *
!*********************************************************************
subroutine constitutive_nonlocal_dotState(dotState, Tstar_v, previousTstar_v, Fe, Fp, Temperature, disorientation, dt_previous, &
                                          state, previousState, relevantState, timestep, g,ip,el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use IO,       only: IO_error
use debug,    only: debugger, &
                    selectiveDebugger, &
                    verboseDebugger
use math,     only: math_norm3, &
                    math_mul6x6, &
                    math_mul3x3, &
                    math_mul33x3, &
                    math_mul33x33, &
                    math_inv3x3, &
                    math_det3x3, &
                    math_Mandel6to33, &
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
real(pReal), dimension(4,mesh_maxNipNeighbors), intent(in) :: &
                                            disorientation            ! crystal disorientation between me and my neighbor (axis, angle pair)
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            Fe, &                     ! elastic deformation gradient
                                            Fp                        ! plastic deformation gradient
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            state, &                  ! current microstructural state
                                            previousState, &          ! previous microstructural state
                                            relevantState             ! relevant microstructural state

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
                                            opposite_n, &             ! index of my opposite neighbor
                                            opposite_ip, &            ! ip of my opposite neighbor
                                            opposite_el, &            ! element index of my opposite neighbor
                                            t, &                      ! type of dislocation
                                            s, &                      ! index of my current slip system
                                            sLattice, &               ! index of my current slip system according to lattice order
                                            i
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),10) :: &
                                            rhoDot, &                 ! density evolution
                                            rhoDotRemobilization, &   ! density evolution by remobilization
                                            rhoDotMultiplication, &   ! density evolution by multiplication
                                            rhoDotFlux, &             ! density evolution by flux
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
                                            neighboring_fluxdensity, &! flux density at neighbroing material point
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
                                            transmissivity, &         ! transmissivity of interfaces for dislocation flux
                                            average_fluxdensity, &    ! average flux density at interface
                                            maximum_fluxdensity, &    ! upper bound for flux density at interface
                                            weight, &                 ! weight for interpolation of flux density
                                            lineLength, &             ! dislocation line length leaving the current interface
                                            D                         ! self diffusion
logical                                     highOrderScheme           ! flag indicating whether we use a high order interpolation scheme or not

if (verboseDebugger .and. selectiveDebugger) then 
  !$OMP CRITICAL (write2out)
    write(6,*) '::: constitutive_nonlocal_dotState at ',g,ip,el
    write(6,*)
  !$OMPEND CRITICAL (write2out)
endif

if (.not. (      mesh_element(2,el) == 1_pInt &
            .or. mesh_element(2,el) == 6_pInt &
            .or. mesh_element(2,el) == 7_pInt &
            .or. mesh_element(2,el) == 8_pInt )) &
  call IO_error(-1,el,ip,g,'element type not supported for nonlocal constitution')

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
  dotState(1,ip,el)%p(1:10*ns) = 0.0_pReal                                                                                          ! ...return without doing anything (-> zero dotState)
  return
endif

if (any(constitutive_nonlocal_v(:,:,g,ip,el)*timestep > mesh_ipVolume(ip,el)**(1.0_pReal/3.0_pReal))) then                          ! if timestep is too large,...
  dotState(1,ip,el)%p(1:10*ns) = NaN                                                                                                ! ...assign NaN and enforce a cutback
  if (verboseDebugger) then 
    !$OMP CRITICAL (write2out)
      write(6,*) 'exceeded maximum allowed dislocation velocity at ',g,ip,el
      write(6,*)
    !$OMPEND CRITICAL (write2out)
  endif
  return
endif


!****************************************************************************
!*** Calculate shear rate

forall (t = 1:4) &
  gdot(:,t) = rhoSgl(:,t) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) * constitutive_nonlocal_v(:,t,g,ip,el)
forall (s = 1:ns, t = 1:4, rhoSgl(s,t+4) * constitutive_nonlocal_v(s,t,g,ip,el) < 0.0_pReal) &                                      ! contribution of used rho for changing sign of v
  gdot(s,t) = gdot(s,t) + abs(rhoSgl(s,t+4)) * constitutive_nonlocal_burgersPerSlipSystem(s,myInstance) &
                                             * constitutive_nonlocal_v(s,t,g,ip,el)

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

rhoDotMultiplication(:,1:2) = spread(0.5_pReal * sum(abs(gdot(:,3:4)),2) * sqrt(rhoForest)  &
                                               / constitutive_nonlocal_lambda0PerSlipSystem(:,myInstance) &
                                               / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance), 2, 2)
rhoDotMultiplication(:,3:4) = spread(0.5_pReal * sum(abs(gdot(:,1:2)),2) * sqrt(rhoForest)  &
                                               / constitutive_nonlocal_lambda0PerSlipSystem(:,myInstance) &
                                               / constitutive_nonlocal_burgersPerSlipSystem(:,myInstance), 2, 2)
rhoDotMultiplication(:,5:10) = 0.0_pReal      ! used dislocation densities and dipoles don't multiplicate


!****************************************************************************
!*** calculate dislocation fluxes

rhoDotFlux = 0.0_pReal

m(:,:,1) =  lattice_sd(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)
m(:,:,2) = -lattice_sd(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)
m(:,:,3) =  lattice_st(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)
m(:,:,4) = -lattice_st(:, constitutive_nonlocal_slipSystemLattice(:,myInstance), myStructure)

F = math_mul33x33(Fe(:,:,g,ip,el), Fp(:,:,g,ip,el))
detFe = math_det3x3(Fe(:,:,g,ip,el)) 

fluxdensity = rhoSgl(:,1:4) * constitutive_nonlocal_v(:,:,g,ip,el)

do n = 1,FE_NipNeighbors(mesh_element(2,el))                                                                                        ! loop through my neighbors
  opposite_n = n - 1_pInt + 2_pInt*mod(n,2_pInt)
  
  neighboring_el = mesh_ipNeighborhood(1,n,ip,el)
  neighboring_ip = mesh_ipNeighborhood(2,n,ip,el)
  
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
  
  transmissivity = constitutive_nonlocal_transmissivity(disorientation(:,n))
  
  highOrderScheme = .false.
  if ( neighboring_el > 0 .and. neighboring_ip > 0 ) then                                                                           ! if neighbor exists...
    if ( .not. phase_localConstitution(material_phase(1,neighboring_ip,neighboring_el))) then                                       !   ... and is of nonlocal constitution...
      forall (t = 1:4) &                                                                                                            !     ... then calculate neighboring flux density
        neighboring_fluxdensity(:,t) = state(g,neighboring_ip,neighboring_el)%p((t-1)*ns+1:t*ns) &
                                     * constitutive_nonlocal_v(:,t,g,neighboring_ip,neighboring_el)
      if (transmissivity == 1.0_pReal) then                                                                                         !     ... if additionally interface's transmission is perfect...
        highOrderScheme = .true.                                                                                                    !       ... then use high order interpolation scheme
        weight = 0.5_pReal * mesh_ipVolume(ip,el) / area &
                           / math_norm3(math_mul33x3(Favg,(  mesh_ipCenterOfGravity(:,neighboring_ip,neighboring_el) &
                                                           - mesh_ipCenterOfGravity(:,ip,el))))
      endif
    else                                                                                                                            !   ... and is of local constitution...
      neighboring_fluxdensity = fluxdensity                                                                                         !     ... then copy flux density to neighbor
    endif
  else                                                                                                                              ! if no neighbor existent...
    neighboring_fluxdensity = 0.0_pReal                                                                                             !   ... assume zero density
  endif
    
  do s = 1,ns
    do t = 1,4
      
      if ( fluxdensity(s,t) * math_mul3x3(m(:,s,t), surfaceNormal) > 0.0_pReal ) then                                               ! outgoing flux        
        if ( highOrderScheme ) then
          average_fluxdensity = fluxdensity(s,t) + weight * (neighboring_fluxdensity(s,t) - fluxdensity(s,t))
          maximum_fluxdensity = rhoSgl(s,t) * mesh_ipVolume(ip,el)**(1.0_pReal/3.0_pReal) / timestep
          average_fluxdensity = min(abs(average_fluxdensity),maximum_fluxdensity) * sign(1.0_pReal,average_fluxdensity)
        else
          average_fluxdensity = fluxdensity(s,t)
        endif

        lineLength = average_fluxdensity * math_mul3x3(m(:,s,t), surfaceNormal) * area                                              !   line length that wants to leave this interface
        rhoDotFlux(s,t) = rhoDotFlux(s,t) - lineLength / mesh_ipVolume(ip,el)                                                       !   subtract positive dislocation flux that leaves the material point
        rhoDotFlux(s,t+4) = rhoDotFlux(s,t+4) + lineLength / mesh_ipVolume(ip,el) * (1.0_pReal - transmissivity) &
                                                                                  * sign(1.0_pReal, average_fluxdensity)            !   dislocation flux that is not able to leave through interface (because of low transmissivity) will remain as immobile single density at the material point        
!        if (selectiveDebugger .and. s==1) & 
!          write(6,'(a22,i2,a15,i2,a3,4(e20.10,x))') 'outgoing flux of type ',t,'   to neighbor ',n,' : ', &
!            -lineLength / mesh_ipVolume(ip,el), average_fluxdensity, maximum_fluxdensity, &
!            average_fluxdensity / maximum_fluxdensity
      else                                                                                                                          ! incoming flux
        if ( highOrderScheme ) then
          average_fluxdensity = fluxdensity(s,t) + weight * (neighboring_fluxdensity(s,t) - fluxdensity(s,t))
          maximum_fluxdensity = state(g,neighboring_ip,neighboring_el)%p((t-1)*ns+s) &
                              * mesh_ipVolume(neighboring_ip,neighboring_el)**(1.0_pReal/3.0_pReal) / timestep
          average_fluxdensity = min(abs(average_fluxdensity),maximum_fluxdensity) * sign(1.0_pReal,average_fluxdensity)
        else
          average_fluxdensity = neighboring_fluxdensity(s,t)
        endif

        lineLength = average_fluxdensity * math_mul3x3(m(:,s,t), surfaceNormal) * area                                              !   line length that wants to leave this interface
        rhoDotFlux(s,t) = rhoDotFlux(s,t) - lineLength / mesh_ipVolume(ip,el) * transmissivity                                      !   subtract negative dislocation flux that enters the material point
!        if (selectiveDebugger .and. s==1) & 
!          write(6,'(a22,i2,a15,i2,a3,4(e20.10,x))') 'incoming flux of type ',t,' from neighbor ',n,' : ', &
!            -lineLength / mesh_ipVolume(ip,el) * transmissivity, average_fluxdensity, maximum_fluxdensity, &
!            average_fluxdensity / maximum_fluxdensity
      endif
    
    enddo
  enddo
  
enddo

constitutive_nonlocal_rhoDotFlux(:,:,g,ip,el) = rhoDotFlux


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
              + rhoDotSingle2DipoleGlide(:,t) &
              + rhoDotAthermalAnnihilation(:,t) &
              + rhoDotRemobilization(:,t) &
              + rhoDotMultiplication(:,t) &
              + rhoDotThermalAnnihilation(:,t)
!              + rhoDotDipole2SingleStressChange(:,t) &
!              + rhoDotSingle2DipoleStressChange(:,t)

dotState(g,ip,el)%p(1:10*ns) = dotState(g,ip,el)%p(1:10*ns) + reshape(rhoDot,(/10*ns/))

do i = 1,4*ns
  if (previousState(g,ip,el)%p(i) + dotState(g,ip,el)%p(i)*timestep < 0.0_pReal) then                                               ! if single mobile densities become negative...
    if (previousState(g,ip,el)%p(i) < relevantState(g,ip,el)%p(i)) then                                                             !   ... and density is already below relevance...
      dotState(g,ip,el)%p(i) = 0.0_pReal                                                                                            !     ... set dotState to zero
    else                                                                                                                            !   ... otherwise...
      if (verboseDebugger) then 
        !$OMP CRITICAL (write2out)
          write(6,*) 'negative dislocation density at ',g,ip,el
          write(6,*)
        !$OMPEND CRITICAL (write2out)
      endif
      dotState(g,ip,el)%p(i) = NaN                                                                                                  !     ... assign NaN and enforce a cutback
    endif
  endif
enddo

if (verboseDebugger .and. selectiveDebugger) then
  !$OMP CRITICAL (write2out)
    write(6,'(a,/,8(12(e12.5,x),/))') 'dislocation remobilization', rhoDotRemobilization(:,1:8) * timestep
    write(6,'(a,/,4(12(e12.5,x),/))') 'dislocation multiplication', rhoDotMultiplication(:,1:4) * timestep
    write(6,'(a,/,8(12(e12.5,x),/))') 'dislocation flux', rhoDotFlux(:,1:8) * timestep
    write(6,'(a,/,10(12(e12.5,x),/))') 'dipole formation by glide', rhoDotSingle2DipoleGlide * timestep
    write(6,'(a,/,2(12(e12.5,x),/))') 'athermal dipole annihilation', rhoDotAthermalAnnihilation(:,1:2) * timestep
    write(6,'(a,/,2(12(e12.5,x),/))') 'thermally activated dipole annihilation', rhoDotThermalAnnihilation(:,9:10) * timestep
!    write(6,'(a,/,10(12(e12.5,x),/))') 'dipole dissociation by stress increase', rhoDotDipole2SingleStressChange * timestep
!    write(6,'(a,/,10(12(e12.5,x),/))') 'dipole formation by stress decrease', rhoDotSingle2DipoleStressChange * timestep
    write(6,'(a,/,10(12(e12.5,x),/))') 'total density change', rhoDot * timestep
  !$OMPEND CRITICAL (write2out)
endif

endsubroutine



!*********************************************************************
!* transmissivity of IP interface                                    *
!*********************************************************************
function constitutive_nonlocal_transmissivity(disorientation)

use prec,     only: pReal, &
                    pInt
use math,     only: math_QuaternionToAxisAngle

implicit none

!* input variables
real(pReal), dimension(4), intent(in) ::    disorientation        ! disorientation as quaternion
                                            
!* output variables
real(pReal) constitutive_nonlocal_transmissivity                  ! transmissivity of an IP interface for dislocations

!* local variables
real(pReal)                                 disorientationAngle
real(pReal), dimension(3) ::                disorientationAxis
real(pReal), dimension(4) ::                disorientationAxisAngle


disorientationAxisAngle = math_QuaternionToAxisAngle(disorientation)

disorientationAxis = disorientationAxisAngle(1:3)
disorientationAngle = disorientationAxisAngle(4)
 
if (disorientationAngle < 3.0_pReal) then
  constitutive_nonlocal_transmissivity = 1.0_pReal
else
  constitutive_nonlocal_transmissivity = 0.5_pReal
endif     
  
endfunction 



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
  gdot(:,t) = rhoSgl(:,t) * constitutive_nonlocal_burgersPerSlipSystem(:,myInstance) * constitutive_nonlocal_v(:,t,g,ip,el)
  
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
