
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
character(len=16), dimension(7), parameter :: constitutive_nonlocal_stateList = (/ 'rhoEdgePos      ', &
                                                                                   'rhoEdgeNeg      ', &
                                                                                   'rhoScrewPos     ', &
                                                                                   'rhoScrewNeg     ', &
                                                                                   'rhoForest       ', &
                                                                                   'tauSlipThreshold', &
                                                                                   'backStress_v    ' /)                            ! list of microstructural state variables
character(len=16), dimension(4), parameter :: constitutive_nonlocal_stateListBasic = constitutive_nonlocal_stateList(1:4)           ! list of "basic" microstructural state variables that are independent from other state variables
character(len=16), dimension(3), parameter :: constitutive_nonlocal_stateListDependent = constitutive_nonlocal_stateList(5:7)       ! list of microstructural state variables that depend on other state variables
real(pReal), parameter :: kB = 1.38e-23_pReal                                                                                       ! Physical parameter, Boltzmann constant in J/Kelvin


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
                                                          constitutive_nonlocal_nu                              ! poisson's ratio
real(pReal), dimension(:,:,:), allocatable ::             constitutive_nonlocal_Cslip_66                        ! elasticity matrix in Mandel notation for each instance
real(pReal), dimension(:,:,:,:,:), allocatable ::         constitutive_nonlocal_Cslip_3333                      ! elasticity matrix for each instance
real(pReal), dimension(:,:), allocatable ::               constitutive_nonlocal_rhoEdgePos0, &                  ! initial edge_pos dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoEdgeNeg0, &                  ! initial edge_neg dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoScrewPos0, &                 ! initial screw_pos dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoScrewNeg0, &                 ! initial screw_neg dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_v0BySlipFamily, &               ! initial dislocation velocity [m/s] for each family and instance
                                                          constitutive_nonlocal_v0BySlipSystem, &               ! initial dislocation velocity [m/s] for each slip system and instance
                                                          constitutive_nonlocal_burgersBySlipFamily, &          ! absolute length of burgers vector [m] for each family and instance
                                                          constitutive_nonlocal_burgersBySlipSystem, &          ! absolute length of burgers vector [m] for each slip system and instance
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
                                            o, &                ! index of my output
                                            s, &                ! index of my slip system
                                            s1, &               ! index of my slip system
                                            s2, &               ! index of my slip system
                                            output, &
                                            mySize
character(len=64)                           tag
character(len=1024)                         line


write(6,*)
write(6,'(a20,a20,a12)') '<<<+-  constitutive_',constitutive_nonlocal_label,' init  -+>>>'
write(6,*)

maxNinstance = count(phase_constitution == constitutive_nonlocal_label)
if (maxNinstance == 0) return                                                                                                       ! we don't have to do anything if there's no instance for this constitutive law


!*** space allocation for global variables

allocate(constitutive_nonlocal_sizeDotState(maxNinstance));                         constitutive_nonlocal_sizeDotState = 0_pInt
allocate(constitutive_nonlocal_sizeState(maxNinstance));                            constitutive_nonlocal_sizeState = 0_pInt
allocate(constitutive_nonlocal_sizePostResults(maxNinstance));                      constitutive_nonlocal_sizePostResults = 0_pInt
allocate(constitutive_nonlocal_sizePostResult(maxval(phase_Noutput), maxNinstance));constitutive_nonlocal_sizePostResult = 0_pInt
allocate(constitutive_nonlocal_output(maxval(phase_Noutput), maxNinstance));        constitutive_nonlocal_output = ''

allocate(constitutive_nonlocal_structureName(maxNinstance));                        constitutive_nonlocal_structureName = ''
allocate(constitutive_nonlocal_structure(maxNinstance));                            constitutive_nonlocal_structure = 0_pInt
allocate(constitutive_nonlocal_Nslip(lattice_maxNslipFamily, maxNinstance));        constitutive_nonlocal_Nslip = 0_pInt
allocate(constitutive_nonlocal_slipFamily(lattice_maxNslip, maxNinstance));         constitutive_nonlocal_slipFamily = 0_pInt
allocate(constitutive_nonlocal_slipSystemLattice(lattice_maxNslip, maxNinstance));  constitutive_nonlocal_slipSystemLattice = 0_pInt
allocate(constitutive_nonlocal_totalNslip(maxNinstance));                           constitutive_nonlocal_totalNslip = 0_pInt

allocate(constitutive_nonlocal_CoverA(maxNinstance));                               constitutive_nonlocal_CoverA = 0.0_pReal 
allocate(constitutive_nonlocal_C11(maxNinstance));                                  constitutive_nonlocal_C11 = 0.0_pReal
allocate(constitutive_nonlocal_C12(maxNinstance));                                  constitutive_nonlocal_C12 = 0.0_pReal
allocate(constitutive_nonlocal_C13(maxNinstance));                                  constitutive_nonlocal_C13 = 0.0_pReal
allocate(constitutive_nonlocal_C33(maxNinstance));                                  constitutive_nonlocal_C33 = 0.0_pReal
allocate(constitutive_nonlocal_C44(maxNinstance));                                  constitutive_nonlocal_C44 = 0.0_pReal
allocate(constitutive_nonlocal_Gmod(maxNinstance));                                 constitutive_nonlocal_Gmod = 0.0_pReal
allocate(constitutive_nonlocal_nu(maxNinstance));                                   constitutive_nonlocal_nu = 0.0_pReal
allocate(constitutive_nonlocal_Cslip_66(6,6,maxNinstance));                         constitutive_nonlocal_Cslip_66 = 0.0_pReal
allocate(constitutive_nonlocal_Cslip_3333(3,3,3,3,maxNinstance));                   constitutive_nonlocal_Cslip_3333 = 0.0_pReal

allocate(constitutive_nonlocal_rhoEdgePos0(lattice_maxNslipFamily, maxNinstance));  constitutive_nonlocal_rhoEdgePos0 = 0.0_pReal
allocate(constitutive_nonlocal_rhoEdgeNeg0(lattice_maxNslipFamily, maxNinstance));  constitutive_nonlocal_rhoEdgeNeg0 = 0.0_pReal
allocate(constitutive_nonlocal_rhoScrewPos0(lattice_maxNslipFamily, maxNinstance)); constitutive_nonlocal_rhoScrewPos0 = 0.0_pReal
allocate(constitutive_nonlocal_rhoScrewNeg0(lattice_maxNslipFamily, maxNinstance)); constitutive_nonlocal_rhoScrewNeg0 = 0.0_pReal
allocate(constitutive_nonlocal_v0BySlipFamily(lattice_maxNslipFamily, maxNinstance));
                                                                              constitutive_nonlocal_v0BySlipFamily = 0.0_pReal
allocate(constitutive_nonlocal_burgersBySlipFamily(lattice_maxNslipFamily, maxNinstance));
                                                                              constitutive_nonlocal_burgersBySlipFamily = 0.0_pReal

allocate(constitutive_nonlocal_interactionSlipSlip(lattice_maxNinteraction, maxNinstance))
                                                                              constitutive_nonlocal_interactionSlipSlip = 0.0_pReal


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
      case ('v0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_v0BySlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case ('burgers')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_burgersBySlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case ('interaction_slipslip')
        forall (f = 1:lattice_maxNinteraction) constitutive_nonlocal_interactionSlipSlip(f,i) = IO_floatValue(line,positions,1+f)
    end select
  endif
enddo


100 do i = 1,maxNinstance

  constitutive_nonlocal_structure(i) = &
    lattice_initializeStructure(constitutive_nonlocal_structureName(i), constitutive_nonlocal_CoverA(i))                            ! our lattice structure is defined in the material.config file by the structureName (and the c/a ratio)
  
  
!*** sanity checks
  
  if (constitutive_nonlocal_structure(i) < 1 .or. constitutive_nonlocal_structure(i) > 3) call IO_error(205)
  
  if (sum(constitutive_nonlocal_Nslip(:,i)) <= 0)                                         call IO_error(225)
  
  if (any(      constitutive_nonlocal_rhoEdgePos0(:,i) < 0.0_pReal &
          .and. constitutive_nonlocal_Nslip(:,i) > 0))                                    call IO_error(220)
  
  if (any(      constitutive_nonlocal_rhoEdgeNeg0(:,i) < 0.0_pReal &
          .and. constitutive_nonlocal_Nslip(:,i) > 0))                                    call IO_error(220)
  
  if (any(      constitutive_nonlocal_rhoScrewPos0(:,i) < 0.0_pReal &
          .and. constitutive_nonlocal_Nslip(:,i) > 0))                                    call IO_error(220)
  
  if (any(      constitutive_nonlocal_rhoScrewNeg0(:,i) < 0.0_pReal &
          .and. constitutive_nonlocal_Nslip(:,i) > 0))                                    call IO_error(220)
  
  if (any(      constitutive_nonlocal_burgersBySlipFamily(:,i) <= 0.0_pReal &
          .and. constitutive_nonlocal_Nslip(:,i) > 0))                                    call IO_error(221)
  
  if (any(      constitutive_nonlocal_v0BySlipFamily(:,i) <= 0.0_pReal &
          .and. constitutive_nonlocal_Nslip(:,i) > 0))                                    call IO_error(-1)
  
  
!*** determine total number of active slip systems
  
  constitutive_nonlocal_Nslip(:,i)  & 
      = min( lattice_NslipSystem(:, constitutive_nonlocal_structure(i)), constitutive_nonlocal_Nslip(:,i) )                         ! we can't use more slip systems per family than specified in lattice 
  constitutive_nonlocal_totalNslip(i) = sum(constitutive_nonlocal_Nslip(:,i))

enddo


!*** allocation of variables whose size depends on the total number of active slip systems

maxTotalNslip = maxval(constitutive_nonlocal_totalNslip)
allocate(constitutive_nonlocal_burgersBySlipSystem(maxTotalNslip, maxNinstance))
                                                                               constitutive_nonlocal_burgersBySlipSystem = 0.0_pReal
allocate(constitutive_nonlocal_v0BySlipSystem(maxTotalNslip, maxNinstance))
                                                                                    constitutive_nonlocal_v0BySlipSystem = 0.0_pReal
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
  
  constitutive_nonlocal_sizeState(i) = (size(constitutive_nonlocal_stateList)-1) * constitutive_nonlocal_totalNslip(i) + 6_pInt     ! the size of the list of states times the number of active slip systems gives the required size for the state array
  constitutive_nonlocal_sizeDotState(i) = size(constitutive_nonlocal_stateListBasic) * constitutive_nonlocal_totalNslip(i)          ! the size of the list of basic states times the number of active slip systems gives the required size for the dotState array
  
  
!*** determine size of postResults array
  
    do o = 1,maxval(phase_Noutput)
    select case(constitutive_nonlocal_output(o,i))
      case( 'dislocationdensity', &
            'shearrate_slip', &
            'resolvedstress_slip', &
            'resistance_slip')
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

  constitutive_nonlocal_Gmod(i) = constitutive_nonlocal_C44(i)                                                                      ! shear modulus is given by elastic constant C44
  constitutive_nonlocal_nu(i) = constitutive_nonlocal_C12(i) / constitutive_nonlocal_C11(i)
  
  
!*** burgers vector and initial dislocation velocity for each slip system
  
  do s = 1,constitutive_nonlocal_totalNslip(i)
    constitutive_nonlocal_burgersBySlipSystem(s,i) &
        = constitutive_nonlocal_burgersBySlipFamily( constitutive_nonlocal_slipFamily(s,i), i )
    constitutive_nonlocal_v0BySlipSystem(s,i) = constitutive_nonlocal_v0BySlipFamily(constitutive_nonlocal_slipFamily(s,i),i)
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
use IO,       only: IO_error
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
    
enddo; enddo


!*** calculate the dependent state variables

! forest dislocation density
forall (s = 1:ns) &
  rhoForest(s) =  dot_product( (rhoEdgePos + rhoEdgeNeg), constitutive_nonlocal_forestProjectionEdge(1:ns, s, myInstance) ) & 
                + dot_product( (rhoScrewPos + rhoScrewNeg), constitutive_nonlocal_forestProjectionScrew(1:ns, s, myInstance) )      ! calculation of forest dislocation density as projection of screw and edge dislocations


! threshold shear stress for dislocation slip 
forall (s = 1:ns) &
  tauSlipThreshold(s) =   constitutive_nonlocal_Gmod(myInstance) & 
                        * constitutive_nonlocal_burgersBySlipSystem(s, myInstance) &
                        * sqrt( dot_product( (rhoEdgePos + rhoEdgeNeg + rhoScrewPos + rhoScrewNeg), &
                                             constitutive_nonlocal_interactionMatrixSlipSlip(1:ns, s, myInstance) ) )


!*** put everything together and in right order

constitutive_nonlocal_stateInit(     1:  ns) = rhoEdgePos
constitutive_nonlocal_stateInit(  ns+1:2*ns) = rhoEdgeNeg
constitutive_nonlocal_stateInit(2*ns+1:3*ns) = rhoScrewPos
constitutive_nonlocal_stateInit(3*ns+1:4*ns) = rhoScrewNeg
constitutive_nonlocal_stateInit(4*ns+1:5*ns) = rhoForest
constitutive_nonlocal_stateInit(5*ns+1:6*ns) = tauSlipThreshold

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
                    math_transpose3x3, &
                    pi
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
real(pReal), dimension(6) ::      backStress_v                ! backstress resulting from the neighboring excess dislocation densities as 2nd Piola-Kirchhoff stress in Mandel notation
real(pReal), dimension(3,3) ::    transform, &                ! orthogonal transformation matrix from slip coordinate system with e1=bxn, e2=b, e3=n to lattice coordinate system 
                                  sigma                       ! backstress resulting from the excess dislocation density of a single slip system and a single neighbor calculated in the coordinate system of the slip system
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                  rhoEdgePos, &               ! positive edge dislocation density
                                  rhoEdgeNeg, &               ! negative edge dislocation density
                                  rhoScrewPos, &              ! positive screw dislocation density
                                  rhoScrewNeg, &              ! negative screw dislocation density
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

rhoEdgePos = state(g,ip,el)%p(      1:  ns)
rhoEdgeNeg = state(g,ip,el)%p(   ns+1:2*ns)
rhoScrewPos = state(g,ip,el)%p(2*ns+1:3*ns)
rhoScrewNeg = state(g,ip,el)%p(3*ns+1:4*ns)


!**********************************************************************
!*** calculate dependent states

!*** calculate the forest dislocation density

forall (s = 1:ns) &
  rhoForest(s) =  dot_product( (rhoEdgePos + rhoEdgeNeg), constitutive_nonlocal_forestProjectionEdge(1:ns, s, myInstance) ) & 
                + dot_product( (rhoScrewPos + rhoScrewNeg), constitutive_nonlocal_forestProjectionScrew(1:ns, s, myInstance) )      ! calculation of forest dislocation density as projection of screw and edge dislocations


!*** calculate the threshold shear stress for dislocation slip 

forall (s = 1:ns) &
  tauSlipThreshold(s) =   constitutive_nonlocal_Gmod(myInstance) & 
                        * constitutive_nonlocal_burgersBySlipSystem(s, myInstance) &
                        * sqrt( dot_product( (rhoEdgePos + rhoEdgeNeg + rhoScrewPos + rhoScrewNeg), &
                                             constitutive_nonlocal_interactionMatrixSlipSlip(1:ns, s, myInstance) ) )
  

!*** calculate the backstress of the neighboring excess dislocation densities

backStress_v = 0.0_pReal


! loop through my neighbors (if it exists!)

do n = 1,FE_NipNeighbors(mesh_element(2,el))

  neighboring_el = mesh_ipNeighborhood(1,n,ip,el)
  neighboring_ip = mesh_ipNeighborhood(2,n,ip,el)
  
  if ( neighboring_el == 0 .or. neighboring_ip == 0 ) cycle

  
  ! calculate the connecting vector between me and my neighbor and his excess dislocation density
  
  connectingVector = math_mul33x3( Fp(:,:,g,ip,el), &
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
    gb = constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgersBySlipSystem(s,myInstance) / (2.0_pReal*pi)
    
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
    backStress_v = backStress_v + math_Mandel33to6( math_mul33x33(math_transpose3x3(transform), math_mul33x33(sigma, transform) ) )
    
  enddo
enddo


!**********************************************************************
!*** set dependent states

state(g,ip,el)%p(4*ns+1:5*ns) = rhoForest
state(g,ip,el)%p(5*ns+1:6*ns) = tauSlipThreshold
state(g,ip,el)%p(6*ns+1:6*ns+6) = backstress_v

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
real(pReal), dimension(6) ::                backStress_v                ! backstress resulting from the neighboring excess dislocation densities as 2nd Piola-Kirchhoff stress
real(pReal), dimension(3,3,3,3) ::          dLp_dTstar3333              ! derivative of Lp with respect to Tstar (3x3x3x3 matrix)
real(pReal), dimension(4,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
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

forall (t = 1:4) rho(t,:) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
rhoForest        = state(g,ip,el)%p(4*ns+1:5*ns)
tauSlipThreshold = state(g,ip,el)%p(5*ns+1:6*ns)
backStress_v     = state(g,ip,el)%p(6*ns+1:6*ns+6)
! if (debugger) write(6,'(a20,3(i3,x),/,3(3(f12.3,x)/))') 'backstress / MPa at ', g,ip,el, math_Mandel6to33(backStress_v/1e6)
! if (debugger) write(6,'(a15,3(i3,x),/,3(3(f12.3,x)/))') 'Tstar / MPa at ',g,ip,el, math_Mandel6to33(Tstar_v/1e6)

!*** loop over slip systems
  
do s = 1,constitutive_nonlocal_totalNslip(myInstance)

  sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)

  
  !*** Calculation of Lp
  
  tauSlip(s) = math_mul6x6( Tstar_v + backStress_v, lattice_Sslip_v(:,sLattice,myStructure) )
  
  if (rhoForest(s) > 0.0_pReal) &
    v(s) =  constitutive_nonlocal_v0BySlipSystem(s,myInstance) &
          * exp( - ( tauSlipThreshold(s) - abs(tauSlip(s)) ) * constitutive_nonlocal_burgersBySlipSystem(s,myInstance)**2.0_pReal &
                   / ( kB * Temperature * sqrt(rhoForest(s)) ) ) &
          * sign(1.0_pReal,tauSlip(s))
  
  gdotSlip(s) =  sum(rho(:,s)) * constitutive_nonlocal_burgersBySlipSystem(s,myInstance) * v(s)
  
  Lp = Lp + gdotSlip(s) * lattice_Sslip(:,:,sLattice,myStructure)
  
  
  !*** Calculation of the tangent of Lp
  
  dgdot_dtauSlip(s) = gdotSlip(s) * constitutive_nonlocal_burgersBySlipSystem(s,myInstance)**2.0_pReal &
                                  / ( kB * Temperature * sqrt(rhoForest(s)) )
  
  forall (i=1:3,j=1:3,k=1:3,l=1:3) &
    dLp_dTstar3333(i,j,k,l) = dLp_dTstar3333(i,j,k,l) + dgdot_dtauSlip(s) * lattice_Sslip(i,j, sLattice,myStructure) &
                                                                          * lattice_Sslip(k,l, sLattice,myStructure) 
enddo
!if (debugger) write(6,'(a26,3(i3,x),/,12(f10.3,x),/)') 'tauSlipThreshold / MPa at ',g,ip,el, tauSlipThreshold/1e6
!if (debugger) write(6,'(a15,3(i3,x),/,12(f10.3,x),/)') 'tauSlip / MPa at ',g,ip,el, tauSlip/1e6
!if (debugger) write(6,'(a15,3(i3,x),/,12(e10.3,x),/)') 'gdotSlip at ',g,ip,el, gdotSlip

dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)

endsubroutine



!*********************************************************************
!* rate of change of microstructure                                  *
!*********************************************************************
subroutine constitutive_nonlocal_dotState(dotState, Tstar_v, Fp, invFp, Temperature, state, g, ip, el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use debug,    only: debugger
use math,     only: math_norm3, &
                    math_mul6x6, &
                    math_mul3x3, &
                    math_mul33x3, &
                    math_transpose3x3
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_element, &
                    FE_NipNeighbors, &
                    mesh_ipNeighborhood, &
                    mesh_ipVolume, &
                    mesh_ipCenterOfGravity, &
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
real(pReal), intent(in) ::                  Temperature               ! temperature
real(pReal), dimension(6), intent(in) ::    Tstar_v                   ! 2nd Piola-Kirchhoff stress in Mandel notation
real(pReal), dimension(3,3), intent(in) ::  Fp, &                     ! plastic deformation gradient
                                            invFp                     ! inverse of plastic deformation gradient
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            state                     ! microstructural state

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
                                            n, &                      ! index of my current neighbor
                                            t, &                      ! type of dislocation
                                            s, &                      ! index of my current slip system
                                            sLattice                  ! index of my current slip system as specified by lattice
real(pReal), dimension(4,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            rho, &                    ! dislocation densities
                                            gdot                      ! shear rates
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            rhoForest, &              ! forest dislocation density
                                            tauSlipThreshold, &       ! threshold shear stress
                                            tauSlip, &                ! resolved shear stress
                                            v                         ! dislocation velocity
real(pReal), dimension(3,4,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            m                         ! direction of dislocation motion
real(pReal), dimension(6) ::                backStress_v              ! backstress resulting from the neighboring excess dislocation densities as 2nd Piola-Kirchhoff stress
real(pReal), dimension(4) ::                leavingRho                ! dislocation densities leaving the current interface
real(pReal), dimension(3) ::                surfaceNormal             ! surface normal of the current interface
real(pReal)                                 norm_surfaceNormal, &     ! euclidic norm of the surface normal
                                            area                      ! area of the current interface

 

myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)


!*** shortcut to state variables 

forall (t = 1:4) rho(t,:) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
rhoForest        = state(g,ip,el)%p(4*ns+1:5*ns)
tauSlipThreshold = state(g,ip,el)%p(5*ns+1:6*ns)
backStress_v     = state(g,ip,el)%p(6*ns+1:6*ns+6)


do s = 1,ns

!*** Calculate gdot for each type of dislocation
  
  sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
  tauSlip(s) = math_mul6x6( Tstar_v + backStress_v, lattice_Sslip_v(:,sLattice,myStructure) )
  
  v(s) =  constitutive_nonlocal_v0BySlipSystem(s,myInstance) &
        * exp( - ( tauSlipThreshold(s) - abs(tauSlip(s)) ) &
               * constitutive_nonlocal_burgersBySlipSystem(s,myInstance)**2.0_pReal &
               / ( kB * Temperature * sqrt(rhoForest(s)) ) ) &
        * sign(1.0_pReal,tauSlip(s))
  
  forall (t = 1:4) gdot(t,s) = rho(t,s) * constitutive_nonlocal_burgersBySlipSystem(s,myInstance) * v(s)

  
!*** Direction of dislocation motion

  m(:,1,s) =  lattice_sd(:, sLattice, myStructure)
  m(:,2,s) = -lattice_sd(:, sLattice, myStructure)
  m(:,3,s) =  lattice_st(:, sLattice, myStructure)
  m(:,4,s) = -lattice_st(:, sLattice, myStructure)

enddo


!*** loop through my neighbors if existent and calculate the area and the surface normal of the interface

do n = 1,FE_NipNeighbors(mesh_element(2,el))

  neighboring_el = mesh_ipNeighborhood(1,n,ip,el)
  neighboring_ip = mesh_ipNeighborhood(2,n,ip,el)
   
  surfaceNormal = math_mul33x3(math_transpose3x3(invFp), mesh_ipAreaNormal(:,n,ip,el))
  norm_surfaceNormal = math_norm3(surfaceNormal)
  surfaceNormal = surfaceNormal / norm_surfaceNormal
  area = mesh_ipArea(n,ip,el) / norm_surfaceNormal

  
!*** loop through my interfaces
!*** subtract dislocation densities that leave through an interface from my dotState and add them to the neighboring dotState
  
  do s = 1,ns
  
    leavingRho = 0.0_pReal
    do t = 1,4
      if ( sign(1.0_pReal,math_mul3x3(m(:,t,s),surfaceNormal)) == sign(1.0_pReal,gdot(t,s)) ) then
        
        leavingRho(t) = gdot(t,s) / constitutive_nonlocal_burgersBySlipSystem(s,myInstance) &
                                  * math_mul3x3(m(:,t,s),surfaceNormal) * area
        
        dotState(1,ip,el)%p((t-1)*ns+s) = dotState(1,ip,el)%p((t-1)*ns+s) - leavingRho(t)
        
        if ( neighboring_el > 0 .and. neighboring_ip > 0 ) then
!*****************************************************************************************************
!***   OMP locking for this neighbor missing
!*****************************************************************************************************
          dotState(1,neighboring_ip,neighboring_el)%p((t-1)*ns+s) = &
                dotState(1,neighboring_ip,neighboring_el)%p((t-1)*ns+s) + leavingRho(t)
        endif
      endif
    enddo
    
  enddo

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
!* INPUT:                                                            *
!*  - Temperature     : temperature                                  *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - dt              : current time increment                       *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
pure function constitutive_nonlocal_postResults(Tstar_v,Temperature,dt,state,ipc,ip,el)

use prec,     only: pReal,pInt,p_vec
use mesh,     only: mesh_NcpElems,mesh_maxNips
use material, only: homogenization_maxNgrains,material_phase,phase_constitutionInstance,phase_Noutput
use lattice,  only: lattice_Sslip_v,lattice_Stwin_v,lattice_maxNslipFamily,lattice_maxNtwinFamily, &
         lattice_NslipSystem,lattice_NtwinSystem,lattice_shearTwin  
implicit none

!* Definition of variables
integer(pInt), intent(in) :: ipc,ip,el
real(pReal), intent(in) :: dt,Temperature
real(pReal), dimension(6), intent(in) :: Tstar_v
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
integer(pInt) matID,structID,ns,nt,f,o,i,c,j,index_myFamily
real(pReal) sumf,tau
real(pReal), dimension(constitutive_nonlocal_sizePostResults(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
 constitutive_nonlocal_postResults

!* Shortened notation
matID    = phase_constitutionInstance(material_phase(ipc,ip,el))
structID = constitutive_nonlocal_structure(matID) 
ns = constitutive_nonlocal_totalNslip(matID)


!* Required output 
c = 0_pInt
constitutive_nonlocal_postResults = 0.0_pReal

do o = 1,phase_Noutput(material_phase(ipc,ip,el))
 select case(constitutive_nonlocal_output(o,matID))

   case ('dislocationdensity')
     constitutive_nonlocal_postResults(c+1:c+ns) = state(ipc,ip,el)%p(1:ns)
     c = c + ns

 end select
enddo

endfunction

END MODULE