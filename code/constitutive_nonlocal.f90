! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!##############################################################
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
character(len=16), dimension(3), parameter :: constitutive_nonlocal_listDependentStates = (/'rhoForest       ', &
                                                                                            'tauThreshold    ', &
                                                                                            'tauBack         ' /) ! list of microstructural state variables that depend on other state variables
character(len=16), dimension(4), parameter ::     constitutive_nonlocal_listOtherStates = (/'velocityEdgePos ', &
                                                                                            'velocityEdgeNeg ', &
                                                                                            'velocityScrewPos', &
                                                                                            'velocityScrewNeg' /) ! list of other dependent state variables that are not updated by microstructure
real(pReal), parameter :: kB = 1.38e-23_pReal                                                                   ! Physical parameter, Boltzmann constant in J/Kelvin

!* Definition of global variables
integer(pInt), dimension(:), allocatable ::               constitutive_nonlocal_sizeDotState, &                 ! number of dotStates = number of basic state variables
                                                          constitutive_nonlocal_sizeDependentState, &           ! number of dependent state variables
                                                          constitutive_nonlocal_sizeState, &                    ! total number of state variables
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
                                                          constitutive_nonlocal_Dsd0, &                         ! prefactor for self-diffusion coefficient
                                                          constitutive_nonlocal_Qsd, &                          ! activation enthalpy for diffusion
                                                          constitutive_nonlocal_aTolRho, &                      ! absolute tolerance for dislocation density in state integration
                                                          constitutive_nonlocal_R, &                            ! cutoff radius for dislocation stress
                                                          constitutive_nonlocal_solidSolutionStrength, &        ! solid solution obstacle strength in Pa
                                                          constitutive_nonlocal_solidSolutionEnergy, &          ! solid solution obstacle strength in Pa
                                                          constitutive_nonlocal_p, &                            ! parameter for kinetic law (Kocks,Argon,Ashby)
                                                          constitutive_nonlocal_q, &                            ! parameter for kinetic law (Kocks,Argon,Ashby)
                                                          constitutive_nonlocal_viscosity, &                    ! viscosity for dislocation glide in Pa s
                                                          constitutive_nonlocal_fattack, &                      ! attack frequency in Hz
                                                          constitutive_nonlocal_rhoSglScatter, &                ! standard deviation of scatter in initial dislocation density
                                                          constitutive_nonlocal_surfaceTransmissivity           ! transmissivity at free surface
real(pReal), dimension(:,:,:), allocatable ::             constitutive_nonlocal_Cslip_66                        ! elasticity matrix in Mandel notation for each instance
real(pReal), dimension(:,:,:,:,:), allocatable ::         constitutive_nonlocal_Cslip_3333                      ! elasticity matrix for each instance
real(pReal), dimension(:,:), allocatable ::               constitutive_nonlocal_rhoSglEdgePos0, &               ! initial edge_pos dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoSglEdgeNeg0, &               ! initial edge_neg dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoSglScrewPos0, &              ! initial screw_pos dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoSglScrewNeg0, &              ! initial screw_neg dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoDipEdge0, &                  ! initial edge dipole dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoDipScrew0, &                 ! initial screw dipole dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_lambda0PerSlipFamily, &         ! mean free path prefactor for each family and instance
                                                          constitutive_nonlocal_lambda0, &                      ! mean free path prefactor for each slip system and instance
                                                          constitutive_nonlocal_burgersPerSlipFamily, &         ! absolute length of burgers vector [m] for each family and instance
                                                          constitutive_nonlocal_burgers, &                      ! absolute length of burgers vector [m] for each slip system and instance
                                                          constitutive_nonlocal_interactionSlipSlip             ! coefficients for slip-slip interaction for each interaction type and instance
real(pReal), dimension(:,:,:), allocatable ::             constitutive_nonlocal_minimumDipoleHeightPerSlipFamily, & ! minimum stable edge/screw dipole height for each family and instance
                                                          constitutive_nonlocal_minimumDipoleHeight, &          ! minimum stable edge/screw dipole height for each slip system and instance
                                                          constitutive_nonlocal_peierlsStressPerSlipFamily, &   ! Peierls stress (edge and screw) 
                                                          constitutive_nonlocal_peierlsStress, &                ! Peierls stress (edge and screw) 
                                                          constitutive_nonlocal_peierlsEnergyPerSlipFamily, &   ! activation energy of peierls barrier (edge and screw)
                                                          constitutive_nonlocal_peierlsEnergy                   ! activation energy of peierls barrier (edge and screw)
real(pReal), dimension(:,:,:,:,:), allocatable ::         constitutive_nonlocal_rhoDotFlux                      ! dislocation convection term
real(pReal), dimension(:,:,:,:,:,:), allocatable ::       constitutive_nonlocal_compatibility                   ! slip system compatibility between me and my neighbors
real(pReal), dimension(:,:,:), allocatable ::             constitutive_nonlocal_forestProjectionEdge, &         ! matrix of forest projections of edge dislocations for each instance
                                                          constitutive_nonlocal_forestProjectionScrew, &        ! matrix of forest projections of screw dislocations for each instance
                                                          constitutive_nonlocal_interactionMatrixSlipSlip       ! interaction matrix of the different slip systems for each instance
real(pReal), dimension(:,:,:,:), allocatable ::           constitutive_nonlocal_lattice2slip, &                 ! orthogonal transformation matrix from lattice coordinate system to slip coordinate system (passive rotation !!!)
                                                          constitutive_nonlocal_accumulatedShear                ! accumulated shear per slip system up to the start of the FE increment


CONTAINS
!****************************************
!* - constitutive_nonlocal_init
!* - constitutive_nonlocal_stateInit
!* - constitutive_nonlocal_aTolState
!* - constitutive_nonlocal_homogenizedC
!* - constitutive_nonlocal_microstructure
!* - constitutive_nonlocal_kinetics
!* - constitutive_nonlocal_LpAndItsTangent
!* - constitutive_nonlocal_dotState
!* - constitutive_nonlocal_dotTemperature
!* - constitutive_nonlocal_updateCompatibility
!* - constitutive_nonlocal_postResults
!****************************************


!**************************************
!*      Module initialization         *
!**************************************
subroutine constitutive_nonlocal_init(file)

use prec,     only: pInt, pReal
use math,     only: math_Mandel3333to66, & 
                    math_Voigt66to3333, & 
                    math_mul3x3, &
                    math_transpose3x3
use IO,       only: IO_lc, &
                    IO_getTag, &
                    IO_isBlank, &
                    IO_stringPos, &
                    IO_stringValue, &
                    IO_floatValue, &
                    IO_intValue, &
                    IO_error
use debug,    only: debug_verbosity
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


!$OMP CRITICAL (write2out)
  write(6,*)
  write(6,*) '<<<+-  constitutive_',trim(constitutive_nonlocal_label),' init  -+>>>'
  write(6,*) '$Id$'
  write(6,*)
!$OMP END CRITICAL (write2out)

maxNinstance = count(phase_constitution == constitutive_nonlocal_label)
if (maxNinstance == 0) return                                                                                                       ! we don't have to do anything if there's no instance for this constitutive law

if (debug_verbosity > 0) then
  !$OMP CRITICAL (write2out)
    write(6,'(a16,x,i5)') '# instances:',maxNinstance
  !$OMP END CRITICAL (write2out)
endif


!*** space allocation for global variables

allocate(constitutive_nonlocal_sizeDotState(maxNinstance))
allocate(constitutive_nonlocal_sizeDependentState(maxNinstance))
allocate(constitutive_nonlocal_sizeState(maxNinstance))
allocate(constitutive_nonlocal_sizePostResults(maxNinstance))
allocate(constitutive_nonlocal_sizePostResult(maxval(phase_Noutput), maxNinstance))
allocate(constitutive_nonlocal_output(maxval(phase_Noutput), maxNinstance))
constitutive_nonlocal_sizeDotState = 0_pInt
constitutive_nonlocal_sizeDependentState = 0_pInt
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
allocate(constitutive_nonlocal_Dsd0(maxNinstance))
allocate(constitutive_nonlocal_Qsd(maxNinstance))
allocate(constitutive_nonlocal_aTolRho(maxNinstance))
allocate(constitutive_nonlocal_Cslip_66(6,6,maxNinstance))
allocate(constitutive_nonlocal_Cslip_3333(3,3,3,3,maxNinstance))
allocate(constitutive_nonlocal_R(maxNinstance))
allocate(constitutive_nonlocal_solidSolutionStrength(maxNinstance))
allocate(constitutive_nonlocal_solidSolutionEnergy(maxNinstance))
allocate(constitutive_nonlocal_p(maxNinstance))
allocate(constitutive_nonlocal_q(maxNinstance))
allocate(constitutive_nonlocal_viscosity(maxNinstance))
allocate(constitutive_nonlocal_fattack(maxNinstance))
allocate(constitutive_nonlocal_rhoSglScatter(maxNinstance))
allocate(constitutive_nonlocal_surfaceTransmissivity(maxNinstance))
constitutive_nonlocal_CoverA = 0.0_pReal 
constitutive_nonlocal_C11 = 0.0_pReal
constitutive_nonlocal_C12 = 0.0_pReal
constitutive_nonlocal_C13 = 0.0_pReal
constitutive_nonlocal_C33 = 0.0_pReal
constitutive_nonlocal_C44 = 0.0_pReal
constitutive_nonlocal_Gmod = 0.0_pReal
constitutive_nonlocal_atomicVolume = 0.0_pReal
constitutive_nonlocal_Dsd0 = 0.0_pReal
constitutive_nonlocal_Qsd = 0.0_pReal
constitutive_nonlocal_aTolRho = 0.0_pReal
constitutive_nonlocal_nu = 0.0_pReal
constitutive_nonlocal_Cslip_66 = 0.0_pReal
constitutive_nonlocal_Cslip_3333 = 0.0_pReal
constitutive_nonlocal_R = -1.0_pReal
constitutive_nonlocal_solidSolutionStrength = 0.0_pReal
constitutive_nonlocal_solidSolutionEnergy = 0.0_pReal
constitutive_nonlocal_p = 1.0_pReal
constitutive_nonlocal_q = 1.0_pReal
constitutive_nonlocal_viscosity = 0.0_pReal
constitutive_nonlocal_fattack = 0.0_pReal
constitutive_nonlocal_rhoSglScatter = 0.0_pReal
constitutive_nonlocal_surfaceTransmissivity = 1.0_pReal

allocate(constitutive_nonlocal_rhoSglEdgePos0(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_nonlocal_rhoSglEdgeNeg0(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_nonlocal_rhoSglScrewPos0(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_nonlocal_rhoSglScrewNeg0(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_nonlocal_rhoDipEdge0(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_nonlocal_rhoDipScrew0(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_nonlocal_burgersPerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_nonlocal_Lambda0PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_nonlocal_interactionSlipSlip(lattice_maxNinteraction,maxNinstance))
constitutive_nonlocal_rhoSglEdgePos0 = -1.0_pReal
constitutive_nonlocal_rhoSglEdgeNeg0 = -1.0_pReal
constitutive_nonlocal_rhoSglScrewPos0 = -1.0_pReal
constitutive_nonlocal_rhoSglScrewNeg0 = -1.0_pReal
constitutive_nonlocal_rhoDipEdge0 = -1.0_pReal
constitutive_nonlocal_rhoDipScrew0 = -1.0_pReal
constitutive_nonlocal_burgersPerSlipFamily = 0.0_pReal
constitutive_nonlocal_lambda0PerSlipFamily = 0.0_pReal
constitutive_nonlocal_interactionSlipSlip = 0.0_pReal

allocate(constitutive_nonlocal_minimumDipoleHeightPerSlipFamily(lattice_maxNslipFamily,2,maxNinstance))
allocate(constitutive_nonlocal_peierlsEnergyPerSlipFamily(lattice_maxNslipFamily,2,maxNinstance))
allocate(constitutive_nonlocal_peierlsStressPerSlipFamily(lattice_maxNslipFamily,2,maxNinstance))
constitutive_nonlocal_minimumDipoleHeightPerSlipFamily = 0.0_pReal
constitutive_nonlocal_peierlsEnergyPerSlipFamily = 0.0_pReal
constitutive_nonlocal_peierlsStressPerSlipFamily = 0.0_pReal

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
    cycle
  endif
  if (section > 0 .and. phase_constitution(section) == constitutive_nonlocal_label) then                                            ! one of my sections
    i = phase_constitutionInstance(section)                                                                                         ! which instance of my constitution is present phase
    positions = IO_stringPos(line,maxNchunks)
    tag = IO_lc(IO_stringValue(line,positions,1))                                                                                   ! extract key
    select case(tag)
      case('constitution','/nonlocal/')
        cycle
      case ('(output)')
        output = output + 1
        constitutive_nonlocal_output(output,i) = IO_lc(IO_stringValue(line,positions,2))
      case ('lattice_structure')
        constitutive_nonlocal_structureName(i) = IO_lc(IO_stringValue(line,positions,2))
      case ('c/a_ratio','covera_ratio')
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
      case ('lambda0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_lambda0PerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case ('burgers')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_burgersPerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case('cutoffradius','r')
        constitutive_nonlocal_R(i) = IO_floatValue(line,positions,2)
      case('minimumdipoleheightedge','ddipminedge')
        forall (f = 1:lattice_maxNslipFamily) & 
          constitutive_nonlocal_minimumDipoleHeightPerSlipFamily(f,1,i) = IO_floatValue(line,positions,1+f)
      case('minimumdipoleheightscrew','ddipminscrew')
        forall (f = 1:lattice_maxNslipFamily) & 
          constitutive_nonlocal_minimumDipoleHeightPerSlipFamily(f,2,i) = IO_floatValue(line,positions,1+f)
      case('atomicvolume')
        constitutive_nonlocal_atomicVolume(i) = IO_floatValue(line,positions,2)
      case('selfdiffusionprefactor','dsd0')
        constitutive_nonlocal_Dsd0(i) = IO_floatValue(line,positions,2)
      case('selfdiffusionenergy','qsd')
        constitutive_nonlocal_Qsd(i) = IO_floatValue(line,positions,2)
      case('atol_rho')
        constitutive_nonlocal_aTolRho(i) = IO_floatValue(line,positions,2)
      case ('interaction_slipslip')
        forall (it = 1:lattice_maxNinteraction) constitutive_nonlocal_interactionSlipSlip(it,i) = IO_floatValue(line,positions,1+it)
      case('peierlsstressedge')
        forall (f = 1:lattice_maxNslipFamily) &
          constitutive_nonlocal_peierlsStressPerSlipFamily(f,1,i) = IO_floatValue(line,positions,1+f)
      case('peierlsstressscrew')
        forall (f = 1:lattice_maxNslipFamily) &
          constitutive_nonlocal_peierlsStressPerSlipFamily(f,2,i) = IO_floatValue(line,positions,1+f)
      case('peierlsenergyedge')
        forall (f = 1:lattice_maxNslipFamily) &
          constitutive_nonlocal_peierlsEnergyPerSlipFamily(f,1,i) = IO_floatValue(line,positions,1+f)
      case('peierlsenergyscrew')
        forall (f = 1:lattice_maxNslipFamily) &
          constitutive_nonlocal_peierlsEnergyPerSlipFamily(f,2,i) = IO_floatValue(line,positions,1+f)
      case('solidsolutionstrength')
        constitutive_nonlocal_solidSolutionStrength(i) = IO_floatValue(line,positions,2)
      case('solidsolutionenergy')
        constitutive_nonlocal_solidSolutionEnergy(i) = IO_floatValue(line,positions,2)
      case('p')
        constitutive_nonlocal_p(i) = IO_floatValue(line,positions,2)
      case('q')
        constitutive_nonlocal_q(i) = IO_floatValue(line,positions,2)
      case('viscosity','glideviscosity')
        constitutive_nonlocal_viscosity(i) = IO_floatValue(line,positions,2)
      case('attackfrequency','fattack')
        constitutive_nonlocal_fattack(i) = IO_floatValue(line,positions,2)
      case('rhosglscatter')
        constitutive_nonlocal_rhoSglScatter(i) = IO_floatValue(line,positions,2)
      case('surfacetransmissivity')
        constitutive_nonlocal_surfaceTransmissivity(i) = IO_floatValue(line,positions,2)
      case default
        call IO_error(236,ext_msg=tag)
    end select
  endif
enddo


100 do i = 1,maxNinstance

  constitutive_nonlocal_structure(i) = &
    lattice_initializeStructure(constitutive_nonlocal_structureName(i), constitutive_nonlocal_CoverA(i))                            ! our lattice structure is defined in the material.config file by the structureName (and the c/a ratio)
  myStructure = constitutive_nonlocal_structure(i)
  
  
  !*** sanity checks
  
  if (myStructure < 1 .or. myStructure > 3)                                   call IO_error(205)
  if (sum(constitutive_nonlocal_Nslip(:,i)) <= 0_pInt)                        call IO_error(235,ext_msg='Nslip')
  do o = 1,maxval(phase_Noutput)
    if(len(constitutive_nonlocal_output(o,i)) > 64)                           call IO_error(666)
  enddo
  do f = 1,lattice_maxNslipFamily
    if (constitutive_nonlocal_Nslip(f,i) > 0_pInt) then
      if (constitutive_nonlocal_rhoSglEdgePos0(f,i) < 0.0_pReal)                call IO_error(235,ext_msg='rhoSglEdgePos0')
      if (constitutive_nonlocal_rhoSglEdgeNeg0(f,i) < 0.0_pReal)                call IO_error(235,ext_msg='rhoSglEdgeNeg0')
      if (constitutive_nonlocal_rhoSglScrewPos0(f,i) < 0.0_pReal)               call IO_error(235,ext_msg='rhoSglScrewPos0')
      if (constitutive_nonlocal_rhoSglScrewNeg0(f,i) < 0.0_pReal)               call IO_error(235,ext_msg='rhoSglScrewNeg0')
      if (constitutive_nonlocal_rhoDipEdge0(f,i) < 0.0_pReal)                   call IO_error(235,ext_msg='rhoDipEdge0')
      if (constitutive_nonlocal_rhoDipScrew0(f,i) < 0.0_pReal)                  call IO_error(235,ext_msg='rhoDipScrew0')
      if (constitutive_nonlocal_burgersPerSlipFamily(f,i) <= 0.0_pReal)         call IO_error(235,ext_msg='burgers')
      if (constitutive_nonlocal_lambda0PerSlipFamily(f,i) <= 0.0_pReal)         call IO_error(235,ext_msg='lambda0')
      if (constitutive_nonlocal_minimumDipoleHeightPerSlipFamily(f,1,i) <= 0.0_pReal) &
                                                                              call IO_error(235,ext_msg='minimumDipoleHeightEdge')
      if (constitutive_nonlocal_minimumDipoleHeightPerSlipFamily(f,2,i) <= 0.0_pReal) &
                                                                              call IO_error(235,ext_msg='minimumDipoleHeightScrew')
      if (constitutive_nonlocal_peierlsStressPerSlipFamily(f,1,i) <= 0.0_pReal) call IO_error(235,ext_msg='peierlsStressEdge')
      if (constitutive_nonlocal_peierlsStressPerSlipFamily(f,2,i) <= 0.0_pReal) call IO_error(235,ext_msg='peierlsStressScrew')
      if (constitutive_nonlocal_peierlsEnergyPerSlipFamily(f,1,i) <= 0.0_pReal) call IO_error(235,ext_msg='peierlsEnergyEdge')
      if (constitutive_nonlocal_peierlsEnergyPerSlipFamily(f,2,i) <= 0.0_pReal) call IO_error(235,ext_msg='peierlsEnergyScrew')
    endif
  enddo
  if (any(constitutive_nonlocal_interactionSlipSlip(1:maxval(lattice_interactionSlipSlip(:,:,myStructure)),i) < 0.0_pReal)) &
                                                                                call IO_error(235,ext_msg='interaction_SlipSlip')
  if (constitutive_nonlocal_R(i) < 0.0_pReal)                                   call IO_error(235,ext_msg='r')
  if (constitutive_nonlocal_atomicVolume(i) <= 0.0_pReal)                       call IO_error(235,ext_msg='atomicVolume')
  if (constitutive_nonlocal_Dsd0(i) <= 0.0_pReal)                               call IO_error(235,ext_msg='selfDiffusionPrefactor')
  if (constitutive_nonlocal_Qsd(i) <= 0.0_pReal)                                call IO_error(235,ext_msg='selfDiffusionEnergy')
  if (constitutive_nonlocal_aTolRho(i) <= 0.0_pReal)                            call IO_error(235,ext_msg='aTol_rho')
  if (constitutive_nonlocal_solidSolutionStrength(i) <= 0.0_pReal)              call IO_error(235,ext_msg='solidSolutionStrength')
  if (constitutive_nonlocal_solidSolutionEnergy(i) <= 0.0_pReal)                call IO_error(235,ext_msg='solidSolutionEnergy')
  if (constitutive_nonlocal_p(i) <= 0.0_pReal .or. constitutive_nonlocal_p(i) > 1.0_pReal) call IO_error(235,ext_msg='p')
  if (constitutive_nonlocal_q(i) < 1.0_pReal .or. constitutive_nonlocal_q(i) > 2.0_pReal) call IO_error(235,ext_msg='q')
  if (constitutive_nonlocal_viscosity(i) <= 0.0_pReal)                          call IO_error(235,ext_msg='viscosity')
  if (constitutive_nonlocal_fattack(i) <= 0.0_pReal)                            call IO_error(235,ext_msg='attackFrequency')
  if (constitutive_nonlocal_rhoSglScatter(i) < 0.0_pReal)                       call IO_error(235,ext_msg='rhoSglScatter')
  if (constitutive_nonlocal_surfaceTransmissivity(i) < 0.0_pReal &
      .or. constitutive_nonlocal_surfaceTransmissivity(i) > 1.0_pReal)          call IO_error(235,ext_msg='surfaceTransmissivity')
  
  
  !*** determine total number of active slip systems
  
  constitutive_nonlocal_Nslip(1:lattice_maxNslipFamily,i) = min( lattice_NslipSystem(1:lattice_maxNslipFamily, myStructure), &
                                                                constitutive_nonlocal_Nslip(1:lattice_maxNslipFamily,i) )           ! we can't use more slip systems per family than specified in lattice 
  constitutive_nonlocal_totalNslip(i) = sum(constitutive_nonlocal_Nslip(1:lattice_maxNslipFamily,i))

enddo


!*** allocation of variables whose size depends on the total number of active slip systems

maxTotalNslip = maxval(constitutive_nonlocal_totalNslip)

allocate(constitutive_nonlocal_burgers(maxTotalNslip, maxNinstance))
constitutive_nonlocal_burgers = 0.0_pReal

allocate(constitutive_nonlocal_lambda0(maxTotalNslip, maxNinstance))
constitutive_nonlocal_lambda0 = 0.0_pReal

allocate(constitutive_nonlocal_minimumDipoleHeight(maxTotalNslip,2,maxNinstance))
constitutive_nonlocal_minimumDipoleHeight = 0.0_pReal

allocate(constitutive_nonlocal_forestProjectionEdge(maxTotalNslip, maxTotalNslip, maxNinstance))
constitutive_nonlocal_forestProjectionEdge = 0.0_pReal

allocate(constitutive_nonlocal_forestProjectionScrew(maxTotalNslip, maxTotalNslip, maxNinstance))
constitutive_nonlocal_forestProjectionScrew = 0.0_pReal

allocate(constitutive_nonlocal_interactionMatrixSlipSlip(maxTotalNslip, maxTotalNslip, maxNinstance))
constitutive_nonlocal_interactionMatrixSlipSlip = 0.0_pReal

allocate(constitutive_nonlocal_lattice2slip(1:3, 1:3, maxTotalNslip, maxNinstance))
constitutive_nonlocal_lattice2slip = 0.0_pReal

allocate(constitutive_nonlocal_accumulatedShear(maxTotalNslip, homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems))
constitutive_nonlocal_accumulatedShear = 0.0_pReal

allocate(constitutive_nonlocal_rhoDotFlux(maxTotalNslip, 10, homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems))
constitutive_nonlocal_rhoDotFlux = 0.0_pReal

allocate(constitutive_nonlocal_compatibility(2,maxTotalNslip, maxTotalNslip, FE_maxNipNeighbors, mesh_maxNips, mesh_NcpElems))
constitutive_nonlocal_compatibility = 0.0_pReal

allocate(constitutive_nonlocal_peierlsEnergy(maxTotalNslip,2,maxNinstance))
constitutive_nonlocal_peierlsEnergy = 0.0_pReal

allocate(constitutive_nonlocal_peierlsStress(maxTotalNslip,2,maxNinstance))
constitutive_nonlocal_peierlsStress = 0.0_pReal

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
  constitutive_nonlocal_sizeDotState(i) = size(constitutive_nonlocal_listBasicStates) * ns
  constitutive_nonlocal_sizeDependentState(i) = size(constitutive_nonlocal_listDependentStates) * ns
  constitutive_nonlocal_sizeState(i) = constitutive_nonlocal_sizeDotState(i) &
                                     + constitutive_nonlocal_sizeDependentState(i) &
                                     + size(constitutive_nonlocal_listOtherStates) * ns

  
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
            'resolvedstress_external', &
            'resolvedstress_back', &
            'resistance', &
            'rho_dot', &
            'rho_dot_sgl', &
            'rho_dot_dip', &
            'rho_dot_gen', &
            'rho_dot_gen_edge', &
            'rho_dot_gen_screw', &
            'rho_dot_sgl2dip', &
            'rho_dot_ann_ath', &
            'rho_dot_ann_the', &
            'rho_dot_flux', &
            'rho_dot_flux_edge', &
            'rho_dot_flux_screw', &
            'velocity_edge_pos', &
            'velocity_edge_neg', &
            'velocity_screw_pos', &
            'velocity_screw_neg', &
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
            'maximumdipoleheight_edge', &
            'maximumdipoleheight_screw', &
            'accumulatedshear' )
        mySize = constitutive_nonlocal_totalNslip(i)
      case('dislocationstress')
        mySize = 6_pInt
      case default
        call IO_error(237,ext_msg=constitutive_nonlocal_output(o,i))
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
  constitutive_nonlocal_Cslip_66(1:6,1:6,i) = math_Mandel3333to66(math_Voigt66to3333(constitutive_nonlocal_Cslip_66(1:6,1:6,i)))
  constitutive_nonlocal_Cslip_3333(1:3,1:3,1:3,1:3,i) = math_Voigt66to3333(constitutive_nonlocal_Cslip_66(1:6,1:6,i))

  constitutive_nonlocal_Gmod(i) = 0.2_pReal * ( constitutive_nonlocal_C11(i) - constitutive_nonlocal_C12(i) &
                                                + 3.0_pReal*constitutive_nonlocal_C44(i) )                                          ! (C11iso-C12iso)/2 with C11iso=(3*C11+2*C12+4*C44)/5 and C12iso=(C11+4*C12-2*C44)/5
  constitutive_nonlocal_nu(i) =   ( constitutive_nonlocal_C11(i) + 4.0_pReal*constitutive_nonlocal_C12(i) &
                                    - 2.0_pReal*constitutive_nonlocal_C44(i) ) &
                                / ( 4.0_pReal*constitutive_nonlocal_C11(i) + 6.0_pReal*constitutive_nonlocal_C12(i) &
                                    + 2.0_pReal*constitutive_nonlocal_C44(i) )                                                      ! C12iso/(C11iso+C12iso) with C11iso=(3*C11+2*C12+4*C44)/5 and C12iso=(C11+4*C12-2*C44)/5
  
  do s1 = 1,ns      
    f = constitutive_nonlocal_slipFamily(s1,i)
    
    !*** burgers vector, mean free path prefactor and minimum dipole distance for each slip system
  
    constitutive_nonlocal_burgers(s1,i) = constitutive_nonlocal_burgersPerSlipFamily(f,i)
    constitutive_nonlocal_lambda0(s1,i) = constitutive_nonlocal_lambda0PerSlipFamily(f,i)
    constitutive_nonlocal_minimumDipoleHeight(s1,1:2,i) = constitutive_nonlocal_minimumDipoleHeightPerSlipFamily(f,1:2,i)
    constitutive_nonlocal_peierlsStress(s1,1:2,i) = constitutive_nonlocal_peierlsStressPerSlipFamily(f,1:2,i)
    constitutive_nonlocal_peierlsEnergy(s1,1:2,i) = constitutive_nonlocal_peierlsEnergyPerSlipFamily(f,1:2,i)

    do s2 = 1,ns
      
      !*** calculation of forest projections for edge and screw dislocations. s2 acts as forest for s1

      constitutive_nonlocal_forestProjectionEdge(s1,s2,i) &
          = abs(math_mul3x3(lattice_sn(1:3,constitutive_nonlocal_slipSystemLattice(s1,i),myStructure), &
                            lattice_st(1:3,constitutive_nonlocal_slipSystemLattice(s2,i),myStructure)))                             ! forest projection of edge dislocations is the projection of (t = b x n) onto the slip normal of the respective slip plane
      
      constitutive_nonlocal_forestProjectionScrew(s1,s2,i) &
          = abs(math_mul3x3(lattice_sn(1:3,constitutive_nonlocal_slipSystemLattice(s1,i),myStructure), &
                            lattice_sd(1:3,constitutive_nonlocal_slipSystemLattice(s2,i),myStructure)))                             ! forest projection of screw dislocations is the projection of b onto the slip normal of the respective splip plane
  
      !*** calculation of interaction matrices

      constitutive_nonlocal_interactionMatrixSlipSlip(s1,s2,i) &
          = constitutive_nonlocal_interactionSlipSlip(lattice_interactionSlipSlip(constitutive_nonlocal_slipSystemLattice(s1,i), &
                                                                                  constitutive_nonlocal_slipSystemLattice(s2,i), &
                                                                                  myStructure), &
                                                      i)
  
    enddo

    !*** rotation matrix from lattice configuration to slip system

    constitutive_nonlocal_lattice2slip(1:3,1:3,s1,i) &
        = math_transpose3x3( reshape((/ lattice_sd(1:3, constitutive_nonlocal_slipSystemLattice(s1,i), myStructure), &
                                       -lattice_st(1:3, constitutive_nonlocal_slipSystemLattice(s1,i), myStructure), &
                                        lattice_sn(1:3, constitutive_nonlocal_slipSystemLattice(s1,i), myStructure)/), (/3,3/)))
  enddo
  
enddo

endsubroutine



!*********************************************************************
!* initial microstructural state (just the "basic" states)           *
!*********************************************************************
function constitutive_nonlocal_stateInit(myInstance)

use prec,     only: pReal, &
                    pInt
use lattice,  only: lattice_maxNslipFamily
use math,     only: math_sampleGaussVar

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
                              rhoDipScrew                     ! screw dipole dislocation density
integer(pInt)                 ns, &                           ! short notation for total number of active slip systems 
                              f, &                            ! index of lattice family
                              from, &
                              upto, &
                              s, &                            ! index of slip system
                              i
real(pReal), dimension(2) ::  noise

constitutive_nonlocal_stateInit = 0.0_pReal
ns = constitutive_nonlocal_totalNslip(myInstance)

!*** set the basic state variables

do f = 1,lattice_maxNslipFamily
  from = 1 + sum(constitutive_nonlocal_Nslip(1:f-1,myInstance))
  upto = sum(constitutive_nonlocal_Nslip(1:f,myInstance))
  do s = from,upto
    do i = 1,2
      noise(i) = math_sampleGaussVar(0.0_pReal, constitutive_nonlocal_rhoSglScatter(myInstance))
    enddo
    rhoSglEdgePos(s) = constitutive_nonlocal_rhoSglEdgePos0(f, myInstance) + noise(1)
    rhoSglEdgeNeg(s) = constitutive_nonlocal_rhoSglEdgeNeg0(f, myInstance) + noise(1)
    rhoSglScrewPos(s) = constitutive_nonlocal_rhoSglScrewPos0(f, myInstance) + noise(2)
    rhoSglScrewNeg(s) = constitutive_nonlocal_rhoSglScrewNeg0(f, myInstance) + noise(2)
  enddo 
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

constitutive_nonlocal_homogenizedC = constitutive_nonlocal_Cslip_66(1:6,1:6,myInstance)
 
endfunction



!*********************************************************************
!* calculates quantities characterizing the microstructure           *
!*********************************************************************
subroutine constitutive_nonlocal_microstructure(state, Temperature, Fe, Fp, g, ip, el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use IO,       only: IO_error
use math,     only: math_Mandel33to6, &
                    math_mul33x33, &
                    math_mul33x3, &
                    math_mul3x3, &
                    math_norm3, &
                    math_inv3x3, &
                    math_invert3x3, &
                    math_transpose3x3, &
                    pi
use debug,    only: debug_verbosity, &
                    debug_selectiveDebugger, &
                    debug_g, &
                    debug_i, &
                    debug_e
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_element, &
                    FE_NipNeighbors, &
                    FE_maxNipNeighbors, &
                    mesh_ipNeighborhood, &
                    mesh_ipCenterOfGravity, &
                    mesh_ipVolume, &
                    mesh_ipAreaNormal
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_localConstitution, &
                    phase_constitutionInstance
use lattice,  only: lattice_sd, &
                    lattice_st

implicit none

!*** input variables
integer(pInt), intent(in) ::    g, &                          ! current grain ID
                                ip, &                         ! current integration point
                                el                            ! current element
real(pReal), intent(in) ::      Temperature                   ! temperature
real(pReal), dimension(3,3), intent(in) :: &
                                Fe, &                         ! elastic deformation gradient
                                Fp                            ! elastic deformation gradient

!*** input/output variables
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(inout) :: &
                                state                         ! microstructural state

!*** output variables

!*** local variables
integer(pInt)                   neighboring_el, &             ! element number of neighboring material point
                                neighboring_ip, &             ! integration point of neighboring material point
                                instance, &                   ! my instance of this constitution
                                neighboring_instance, &       ! instance of this constitution of neighboring material point
                                latticeStruct, &              ! my lattice structure
                                neighboring_latticeStruct, &  ! lattice structure of neighboring material point
                                phase, &
                                neighboring_phase, &
                                ns, &                         ! total number of active slip systems at my material point
                                neighboring_ns, &             ! total number of active slip systems at neighboring material point
                                c, &                          ! index of dilsocation character (edge, screw)
                                s, &                          ! slip system index
                                t, &                          ! index of dilsocation type (e+, e-, s+, s-, used e+, used e-, used s+, used s-)
                                dir, &
                                side, &
                                n
integer(pInt), dimension(2) ::  neighbor
real(pReal)                     nu, &                         ! poisson's ratio
                                mu, &
                                b, &
                                detFe, &
                                detFp, &
                                FVsize, &
                                rhoExcessGradient
real(pReal), dimension(2) ::    rhoExcessGradient_over_rho, &
                                gradient, &
                                gradientDeads, &
                                gradientInter, &
                                gradientDistance, &
                                gradientDistanceDeads, &
                                gradientDistanceInter, &
                                rhoExcessAtSampledPoint
real(pReal), dimension(3) ::    ipCoords, &
                                neighboring_ipCoords
real(pReal), dimension(FE_maxNipNeighbors) :: &
                                distance                      ! length of connection vector
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                rhoForest, &                  ! forest dislocation density
                                tauBack, &                    ! back stress from pileup on same slip system
                                tauThreshold                  ! threshold shear stress
real(pReal), dimension(3,2) ::  rhoExcessDifferences, &
                                sampledPoint
real(pReal), dimension(3,3) ::  invFe, &                      ! inverse of elastic deformation gradient
                                invFp                         ! inverse of plastic deformation gradient
real(pReal), dimension(3,FE_maxNipNeighbors) :: &
                                connection_latticeConf, &
                                areaNormal_latticeConf
real(pReal), dimension(2,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                rhoExcess
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                rhoDip                        ! dipole dislocation density (edge, screw)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),8) :: &
                                rhoSgl                        ! single dislocation density (edge+, edge-, screw+, screw-, used edge+, used edge-, used screw+, used screw-)
real(pReal), dimension(3,3,2) :: connections
real(pReal), dimension(2,maxval(constitutive_nonlocal_totalNslip),FE_maxNipNeighbors) :: &
                                neighboring_rhoExcess         ! excess density at neighboring material point
real(pReal), dimension(3,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                m                             ! direction of dislocation motion
logical                         inversionError


phase = material_phase(g,ip,el)
instance = phase_constitutionInstance(phase)
latticeStruct = constitutive_nonlocal_structure(instance)
ns = constitutive_nonlocal_totalNslip(instance)



!*** get basic states

forall (s = 1:ns, t = 1:4) &
  rhoSgl(s,t) = max(state(g,ip,el)%p((t-1)*ns+s), 0.0_pReal)                                                        ! ensure positive single mobile densities
forall (t = 5:8) & 
  rhoSgl(1:ns,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (s = 1:ns, c = 1:2) &
  rhoDip(s,c) = max(state(g,ip,el)%p((7+c)*ns+s), 0.0_pReal)                                                        ! ensure positive dipole densities



!*** calculate the forest dislocation density
!*** (= projection of screw and edge dislocations)

forall (s = 1:ns) &
  rhoForest(s) = dot_product((sum(abs(rhoSgl(1:ns,(/1,2,5,6/))),2) + rhoDip(1:ns,1)), &
                              constitutive_nonlocal_forestProjectionEdge(s,1:ns,instance)) & 
                + dot_product((sum(abs(rhoSgl(1:ns,(/3,4,7,8/))),2) + rhoDip(1:ns,2)), &
                              constitutive_nonlocal_forestProjectionScrew(s,1:ns,instance))



!*** calculate the threshold shear stress for dislocation slip 

forall (s = 1:ns) &
  tauThreshold(s) =   constitutive_nonlocal_Gmod(instance) * constitutive_nonlocal_burgers(s,instance) &
                    * sqrt(dot_product((sum(abs(rhoSgl),2) + sum(abs(rhoDip),2)), &
                                        constitutive_nonlocal_interactionMatrixSlipSlip(s,1:ns,instance)))



!*** calculate the dislocation stress of the neighboring excess dislocation densities
!*** zero for material points of local constitution

tauBack = 0.0_pReal

if (.not. phase_localConstitution(phase)) then
  call math_invert3x3(Fe, invFe, detFe, inversionError)
  call math_invert3x3(Fp, invFp, detFp, inversionError)
  ipCoords = mesh_ipCenterOfGravity(1:3,ip,el)
  rhoExcess(1,1:ns) = rhoSgl(1:ns,1) - rhoSgl(1:ns,2)
  rhoExcess(2,1:ns) = rhoSgl(1:ns,3) - rhoSgl(1:ns,4)
  FVsize = mesh_ipVolume(ip,el) ** (1.0_pReal/3.0_pReal)
  nu = constitutive_nonlocal_nu(instance)
  mu = constitutive_nonlocal_Gmod(instance)
  
  !* loop through my neighborhood and get the connection vectors (in lattice frame) and the excess densities
  
  do n = 1,FE_NipNeighbors(mesh_element(2,el))
    neighboring_el = mesh_ipNeighborhood(1,n,ip,el)
    neighboring_ip = mesh_ipNeighborhood(2,n,ip,el)
    areaNormal_latticeConf(1:3,n) = detFp * math_mul33x3(math_transpose3x3(invFp), mesh_ipAreaNormal(1:3,n,ip,el))  ! calculate the normal of the interface in lattice configuration
    areaNormal_latticeConf(1:3,n) = areaNormal_latticeConf(1:3,n) / math_norm3(areaNormal_latticeConf(1:3,n))       ! normalize the surface normal to unit length
    if (neighboring_el > 0 .and. neighboring_ip > 0) then
      neighboring_phase = material_phase(g,neighboring_ip,neighboring_el)
      neighboring_instance = phase_constitutionInstance(neighboring_phase)
      neighboring_latticeStruct = constitutive_nonlocal_structure(neighboring_instance)
      neighboring_ns = constitutive_nonlocal_totalNslip(neighboring_instance)
      neighboring_ipCoords = mesh_ipCenterOfGravity(1:3,neighboring_ip,neighboring_el)
      if (.not. phase_localConstitution(neighboring_phase) &
          .and. neighboring_latticeStruct == latticeStruct & 
          .and. neighboring_instance == instance) then
        if (neighboring_ns == ns) then
          if (neighboring_el /= el .or. neighboring_ip /= ip) then
            connection_latticeConf(1:3,n) = math_mul33x3(invFe, neighboring_ipCoords - ipCoords)
            forall (s = 1:ns, c = 1:2) &
              neighboring_rhoExcess(c,s,n) = state(g,neighboring_ip,neighboring_el)%p((2*c-2)*ns+s) &  ! positive mobiles
                                           - state(g,neighboring_ip,neighboring_el)%p((2*c-1)*ns+s)    ! negative mobiles
          else
            ! thats myself! probably using periodic images
            connection_latticeConf(1:3,n) = areaNormal_latticeConf(1:3,n) * FVsize
            neighboring_rhoExcess(1:2,1:ns,n) = rhoExcess
          endif
        else
          ! different number of active slip systems
          call IO_error(-1,ext_msg='different number of active slip systems in neighboring IPs of same crystal structure')
        endif
      else
        ! local neighbor or different lattice structure or different constitution instance
        connection_latticeConf(1:3,n) = math_mul33x3(invFe, neighboring_ipCoords - ipCoords)
        neighboring_rhoExcess(1:2,1:ns,n) = rhoExcess
      endif
    else
      ! free surface
      connection_latticeConf(1:3,n) = areaNormal_latticeConf(1:3,n) * FVsize
      neighboring_rhoExcess(1:2,1:ns,n) = rhoExcess
    endif
    distance(n) = math_norm3(connection_latticeConf(1:3,n))
  enddo
  
  !* loop through the slip systems
  !* calculate the dislocation gradient in both directions of m with two different methods:
  !* 1. gradient between central excess density and dead dislocations in central ip
  !* 2. interpolate gradient from excess density in three neighboring ips
  !* take the heigher gradient in both directions and do a weighted sum with weights according to the distance
  
  m(1:3,1:ns,1) =  lattice_sd(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,instance), latticeStruct)
  m(1:3,1:ns,2) = -lattice_st(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,instance), latticeStruct)

  do s = 1,ns
    rhoExcessGradient_over_rho = 0.0_pReal
    do c = 1,2
      if (rhoSgl(s,2*c-1) + rhoSgl(s,2*c) < 1.0_pReal) then
        cycle       ! no siginificant density
      endif
    
      !* gradient from dead dislocations
      
      gradientDeads = 0.0_pReal
      if (rhoSgl(s,2*c+3) > 0.0_pReal) then                                         ! positive deads
        gradientDeads(1) = + 2.0_pReal * rhoSgl(s,2*c+3) / FVsize                   ! on positive side
      else
        gradientDeads(2) = - 2.0_pReal * rhoSgl(s,2*c+3) / FVsize                   ! on negative side
      endif
      if (rhoSgl(s,2*c+4) > 0.0_pReal) then                                         ! negative deads
        gradientDeads(2) = gradientDeads(2) + 2.0_pReal * rhoSgl(s,2*c+4) / FVsize  ! on negative side
      else
        gradientDeads(1) = gradientDeads(1) - 2.0_pReal * rhoSgl(s,2*c+4) / FVsize  ! on positive side
      endif
      gradientDistanceDeads(1:2) = 0.5_pReal * FVsize
      
      !* gradient from interpolation

      gradientInter = 0.0_pReal
      rhoExcessDifferences = 0.0_pReal
      connections = 0.0_pReal
      gradientDistanceInter = 0.0_pReal
      do dir = 1,3
        if (math_mul3x3(areaNormal_latticeConf(1:3,2*dir-1),m(1:3,s,c)) > 0.0_pReal) then ! on positive side
          neighbor(1) = 2 * dir - 1
          neighbor(2) = 2 * dir
        else                                                                              ! on negative side
          neighbor(1) = 2 * dir
          neighbor(2) = 2 * dir - 1
        endif
        do side = 1,2
          n = neighbor(side)
          rhoExcessDifferences(dir,side) = neighboring_rhoExcess(c,s,n) - rhoExcess(c,s)
          connections(dir,1:3,side) = connection_latticeConf(1:3,n)
          gradientDistanceInter(side) = gradientDistanceInter(side) &
                                      + (math_mul3x3(connection_latticeConf(1:3,n),m(1:3,s,c))) ** 2.0_pReal / distance(n)
        enddo
      enddo
      sampledPoint(1:3,1) = + gradientDistanceInter(1) * m(1:3,s,c)
      sampledPoint(1:3,2) = - gradientDistanceInter(2) * m(1:3,s,c)
      do side = 1,2
        rhoExcessAtSampledPoint(side) = math_mul3x3(math_mul33x3(math_inv3x3(connections(1:3,1:3,side)), &
                                                                 rhoExcessDifferences(1:3,side)), &
                                                    sampledPoint(1:3,side)) &
                                      + rhoExcess(c,s)
      enddo
      gradientInter(1) = (rhoExcessAtSampledPoint(1) - rhoExcess(c,s)) / gradientDistanceInter(1)
      gradientInter(2) = (rhoExcess(c,s) - rhoExcessAtSampledPoint(2)) / gradientDistanceInter(2)
       
      !* take maximum of both gradients and mix contributions from both sides according to weighted distances

      do dir = 1,2
        if (abs(gradientDeads(dir)) > abs(gradientInter(dir))) then
          gradient(dir) = gradientDeads(dir)
          gradientDistance(dir) = gradientDistanceDeads(dir)
        else
          gradient(dir) = gradientInter(dir)
          gradientDistance(dir) = gradientDistanceInter(dir)
        endif
      enddo 
      rhoExcessGradient = (gradient(1) * gradientDistance(2) + gradient(2) * gradientDistance(1)) &
                        / (gradientDistance(1) + gradientDistance(2))
      
      !* excess gradient over density: in case of vanishing central total density we take the distance squared instead!!!

      rhoExcessGradient_over_rho(c) = rhoExcessGradient / (rhoSgl(s,2*c-1) + rhoSgl(s,2*c))
    enddo
    
    b = constitutive_nonlocal_burgers(s,instance)
    tauBack(s) = - mu * b / (2.0_pReal * pi) * (rhoExcessGradient_over_rho(1) / (1.0_pReal - nu) + rhoExcessGradient_over_rho(2))

  enddo
endif


!*** set dependent states

state(g,ip,el)%p(10*ns+1:11*ns) = rhoForest
state(g,ip,el)%p(11*ns+1:12*ns) = tauThreshold
state(g,ip,el)%p(12*ns+1:13*ns) = tauBack


#ifndef _OPENMP
  if (debug_verbosity > 6 .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g) .or. .not. debug_selectiveDebugger)) then
    write(6,*)
    write(6,'(a,i8,x,i2,x,i1)') '<< CONST >> nonlocal_microstructure at el ip g',el,ip,g
    write(6,*)
    write(6,'(a,/,12(x),12(e10.3,x))') '<< CONST >> rhoForest', rhoForest
    write(6,'(a,/,12(x),12(f10.5,x))') '<< CONST >> tauThreshold / MPa', tauThreshold/1e6
    write(6,'(a,/,12(x),12(f10.5,x))') '<< CONST >> tauBack / MPa', tauBack/1e6
  endif
#endif

endsubroutine



!*********************************************************************
!* calculates kinetics                                               *
!*********************************************************************
subroutine constitutive_nonlocal_kinetics(v, tau, c, Temperature, state, g, ip, el, dv_dtau)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use math,     only: math_mul6x6, &
                    math_Mandel6to33
use debug,    only: debug_verbosity, &
                    debug_selectiveDebugger, &
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
                                            el, &                       ! current element number
                                            c                           ! dislocation character (1:edge, 2:screw)
real(pReal), intent(in) ::                  Temperature                 ! temperature
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))), &
             intent(in) ::                  tau                         ! resolved external shear stress (for bcc this already contains non Schmid effects)
type(p_vec), intent(in) ::                  state                       ! microstructural state

!*** input/output variables

!*** output variables
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))), &
                            intent(out) ::  v                           ! velocity
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))), &
                   intent(out), optional :: dv_dtau                     ! velocity derivative with respect to resolved shear stress

!*** local variables
integer(pInt)                               myInstance, &               ! current instance of this constitution
                                            myStructure, &              ! current lattice structure
                                            ns, &                       ! short notation for the total number of active slip systems
                                            s                           ! index of my current slip system
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            tauThreshold, &             ! threshold shear stress
                                            rhoForest, &                ! forest dislocation density
                                            meanfreepath, &             ! mean free travel distance for dislocations between two strong obstacles
                                            sweaptArea, &               ! area that is swept when one strong obstacle is surmounted
                                            b                           ! shortcut for burgers vector length
real(pReal)                                 tauRelPeierls, & 
                                            tauRelSS, &
                                            tPeierls, &                 ! waiting time in front of a peierls barriers
                                            tSS, &                      ! waiting time in front of a solid solution obstacle
                                            tViscous, &                 ! travel time for mean freepath in case of viscous glide
                                            dtPeierls_dtau, &           ! derivative with respect to resolved shear stress
                                            dtSS_dtau, &                ! derivative with respect to resolved shear stress
                                            dtViscous_dtau, &           ! derivative with respect to resolved shear stress
                                            p, &                        ! shortcut to Kocks,Argon,Ashby parameter p
                                            q                           ! shortcut to Kocks,Argon,Ashby parameter q


myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)

rhoForest = state%p(10*ns+1:11*ns)
tauThreshold = state%p(11*ns+1:12*ns)
meanfreepath = 1.0_pReal / sqrt(rhoForest)
sweaptArea = 1.0_pReal / rhoForest

p = constitutive_nonlocal_p(myInstance)
q = constitutive_nonlocal_q(myInstance)
b = constitutive_nonlocal_burgers(1:ns,myInstance)

v = 0.0_pReal
if (present(dv_dtau)) dv_dtau = 0.0_pReal


if (Temperature > 0.0_pReal) then
  do s = 1,ns
    if (abs(tau(s)) > tauThreshold(s)) then
      
      !* Peierls contribution
      !* The derivative only gives absolute values; the correct sign is taken care of in the formula for the derivative of the velocity

      tauRelPeierls = (abs(tau(s)) - tauThreshold(s)) / constitutive_nonlocal_peierlsStress(s,c,myInstance)
      tPeierls = constitutive_nonlocal_fattack(myInstance) &
               * exp(constitutive_nonlocal_peierlsEnergy(s,c,myInstance) / (kB * Temperature) * (1.0_pReal - tauRelPeierls**p)**q )
      if (present(dv_dtau)) then
        dtPeierls_dtau = tPeierls * p * q * constitutive_nonlocal_peierlsEnergy(s,c,myInstance) & 
                       / (kB * Temperature * constitutive_nonlocal_peierlsStress(s,c,myInstance)) &
                       * (1.0_pReal - tauRelPeierls**p)**(q-1.0_pReal) * tauRelPeierls**(p-1.0_pReal) 
      endif


      !* Contribution from solid solution strengthening
      !* The derivative only gives absolute values; the correct sign is taken care of in the formula for the derivative of the velocity

      tauRelSS = (abs(tau(s)) - tauThreshold(s)) / constitutive_nonlocal_solidSolutionStrength(myInstance)
      tSS = constitutive_nonlocal_fattack(myInstance) &
          * exp(constitutive_nonlocal_solidSolutionEnergy(myInstance) / (kB * Temperature) * (1.0_pReal - tauRelSS**p)**q )
      if (present(dv_dtau)) then
        dtSS_dtau = tSS * p * q * constitutive_nonlocal_solidSolutionEnergy(myInstance) & 
                       / (kB * Temperature * constitutive_nonlocal_solidsolutionStrength(myInstance)) &
                       * (1.0_pReal - tauRelSS**p)**(q-1.0_pReal) * tauRelSS**(p-1.0_pReal) 
      endif


      !* Contribution from viscous glide
      !* The derivative only gives absolute values; the correct sign is taken care of in the formula for the derivative of the velocity
      
      tViscous = meanfreepath(s) * constitutive_nonlocal_viscosity(myInstance) / (b(s) * abs(tau(s)))
      dtViscous_dtau = tViscous / abs(tau(s))


      !* velocity = travel distance over travel time times correction term for backward jumps
      
      v(s) = meanfreepath(s) / (tPeierls + tSS + tViscous) &
           * (1.0_pReal - exp(-(abs(tau(s)) - tauThreshold(s)) * sweaptArea(s) * b(s) / (kB * Temperature)))
      v(s) = sign(v(s),tau(s))
      if (present(dv_dtau)) then
        dv_dtau(s) = abs(v(s)) * (dtPeierls_dtau + dtSS_dtau + dtViscous_dtau) / (tPeierls + tSS + tViscous) &
                   + sweaptArea(s) * b(s) / (kB * Temperature) * meanfreepath(s) / (tPeierls + tSS + tViscous) &
                                * exp(-(abs(tau(s)) - tauThreshold(s)) * sweaptArea(s) * b(s) / (kB * Temperature))
      endif

    endif
  enddo
endif
    

#ifndef _OPENMP
  if (debug_verbosity > 6 .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g) .or. .not. debug_selectiveDebugger)) then
    write(6,*)
    write(6,'(a,i8,x,i2,x,i1)') '<< CONST >> nonlocal_kinetics at el ip g',el,ip,g
    write(6,*)
    write(6,'(a,/,12(x),12(f12.5,x))') '<< CONST >> tau / MPa', tau / 1e6_pReal
    write(6,'(a,/,4(12(x),12(f12.5,x),/))') '<< CONST >> v / 1e-3m/s', v * 1e3
  endif
#endif

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
use debug,    only: debug_verbosity, &
                    debug_selectiveDebugger, &
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
real(pReal), dimension(6), intent(in) ::    Tstar_v                     ! 2nd Piola-Kirchhoff stress in Mandel notation

!*** input/output variables
type(p_vec), intent(inout) ::               state                       ! microstructural state

!*** output variables
real(pReal), dimension(3,3), intent(out) :: Lp                          ! plastic velocity gradient
real(pReal), dimension(9,9), intent(out) :: dLp_dTstar99                ! derivative of Lp with respect to Tstar (9x9 matrix)

!*** local variables
integer(pInt)                               myInstance, &               ! current instance of this constitution
                                            myStructure, &              ! current lattice structure
                                            ns, &                       ! short notation for the total number of active slip systems
                                            c, &
                                            i, &
                                            j, &
                                            k, &
                                            l, &
                                            t, &                        ! dislocation type
                                            s, &                        ! index of my current slip system
                                            sLattice                    ! index of my current slip system according to lattice order
real(pReal), dimension(3,3,3,3) ::          dLp_dTstar3333              ! derivative of Lp with respect to Tstar (3x3x3x3 matrix)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            rhoSgl, &                   ! single dislocation densities (including used) 
                                            v, &                        ! velocity
                                            dv_dtau                     ! velocity derivative with respect to the shear stress
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            tau, &                      ! resolved shear stress including non Schmid and backstress terms
                                            gdotTotal, &                ! shear rate
                                            dgdotTotal_dtau, &          ! derivative of the shear rate with respect to the shear stress
                                            tauBack                     ! back stress from dislocation gradients on same slip system


!*** initialize local variables

Lp = 0.0_pReal
dLp_dTstar3333 = 0.0_pReal

myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)


!*** shortcut to state variables 

forall (s = 1:ns, t = 1:4) &
  rhoSgl(s,t) = max(state%p((t-1)*ns+s), 0.0_pReal)
tauBack = state%p(12*ns+1:13*ns)


!*** get effective resolved shear stress

do s = 1,ns
  tau(s) = math_mul6x6(Tstar_v, lattice_Sslip_v(:,constitutive_nonlocal_slipSystemLattice(s,myInstance),myStructure)) &
         + tauBack(s)
enddo


!*** get dislocation velocity and its tangent and store the velocity in the state array

if (myStructure == 1_pInt) then   ! for fcc all velcities are equal
  call constitutive_nonlocal_kinetics(v(1:ns,1), tau, 1, Temperature, state, g, ip, el, dv_dtau(1:ns,1))
  do t = 1,4
    v(1:ns,t) = v(1:ns,1)
    dv_dtau(1:ns,t) = dv_dtau(1:ns,1)
    state%p((12+t)*ns+1:(13+t)*ns) = v(1:ns,1)
  enddo
else                              ! for all other lattice structures the velcities may vary with character and sign
  do t = 1,4
    c = (t-1)/2+1
    call constitutive_nonlocal_kinetics(v(1:ns,t), tau, c, Temperature, state, g, ip, el, dv_dtau(1:ns,t))
    state%p((12+t)*ns+1:(13+t)*ns) = v(1:ns,t)
  enddo
endif


!*** Bauschinger effect

forall (s = 1:ns, t = 5:8, state%p((t-1)*ns+s) * v(s,t-4) < 0.0_pReal) &
  rhoSgl(s,t-4) = rhoSgl(s,t-4) + abs(state%p((t-1)*ns+s))


!*** Calculation of gdot and its tangent

gdotTotal = sum(rhoSgl * v, 2) * constitutive_nonlocal_burgers(1:ns,myInstance)
dgdotTotal_dtau = sum(rhoSgl * dv_dtau, 2) * constitutive_nonlocal_burgers(1:ns,myInstance) 


!*** Calculation of Lp and its tangent

do s = 1,ns
  sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)  
  Lp = Lp + gdotTotal(s) * lattice_Sslip(1:3,1:3,sLattice,myStructure)
  forall (i=1:3,j=1:3,k=1:3,l=1:3) &
    dLp_dTstar3333(i,j,k,l) = dLp_dTstar3333(i,j,k,l) + dgdotTotal_dtau(s) * lattice_Sslip(i,j, sLattice,myStructure) &
                                                                           * lattice_Sslip(k,l, sLattice,myStructure) 
enddo
dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)


#ifndef _OPENMP
  if (debug_verbosity > 6 .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g) .or. .not. debug_selectiveDebugger)) then
    write(6,*)
    write(6,'(a,i8,x,i2,x,i1)') '<< CONST >> nonlocal_LpandItsTangent at el ip g ',el,ip,g
    write(6,*)
    write(6,'(a,/,12(x),12(f12.5,x))') '<< CONST >> gdot total / 1e-3',gdotTotal*1e3_pReal
    write(6,'(a,/,3(12(x),3(f12.7,x),/))') '<< CONST >> Lp',Lp
  endif
#endif

endsubroutine



!*********************************************************************
!* rate of change of microstructure                                  *
!*********************************************************************
subroutine constitutive_nonlocal_dotState(dotState, Tstar_v, Fe, Fp, Temperature, state, aTolState, timestep, orientation, g,ip,el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec, &
                    DAMASK_NaN
use numerics, only: numerics_integrationMode
use IO,       only: IO_error
use debug,    only: debug_verbosity, &
                    debug_selectiveDebugger, &
                    debug_g, &
                    debug_i, &
                    debug_e
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
                    pi                
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
                    phase_localConstitution, &
                    phase_constitution
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
                                            timestep                  ! substepped crystallite time increment
real(pReal), dimension(6), intent(in) ::    Tstar_v                   ! current 2nd Piola-Kirchhoff stress in Mandel notation
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            Fe, &                     ! elastic deformation gradient
                                            Fp                        ! plastic deformation gradient
real(pReal), dimension(4,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            orientation               ! crystal lattice orientation
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            state, &                  ! current microstructural state
                                            aTolState                 ! absolute state tolerance

!*** input/output variables
type(p_vec), intent(inout) ::               dotState                  ! evolution of state variables / microstructure
 
!*** output variables
 
!*** local variables
integer(pInt)                               myInstance, &             ! current instance of this constitution
                                            myStructure, &            ! current lattice structure
                                            ns, &                     ! short notation for the total number of active slip systems
                                            c, &                      ! character of dislocation
                                            n, &                      ! index of my current neighbor
                                            neighboring_el, &         ! element number of my neighbor
                                            neighboring_ip, &         ! integration point of my neighbor
                                            neighboring_n, &          ! neighbor index pointing to me when looking from my neighbor
                                            opposite_n, &             ! index of my opposite neighbor
                                            opposite_ip, &            ! ip of my opposite neighbor
                                            opposite_el, &            ! element index of my opposite neighbor
                                            t, &                      ! type of dislocation
                                            topp, &                   ! type of dislocation with opposite sign to t
                                            s, &                      ! index of my current slip system
                                            sLattice                  ! index of my current slip system according to lattice order
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),10) :: &
                                            rhoDot, &                     ! density evolution
                                            rhoDotRemobilization, &       ! density evolution by remobilization
                                            rhoDotMultiplication, &       ! density evolution by multiplication
                                            rhoDotFlux, &                 ! density evolution by flux
                                            rhoDotSingle2DipoleGlide, &   ! density evolution by dipole formation (by glide)
                                            rhoDotAthermalAnnihilation, & ! density evolution by athermal annihilation
                                            rhoDotThermalAnnihilation     ! density evolution by thermal annihilation
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),8) :: &
                                            rhoSgl                        ! current single dislocation densities (positive/negative screw and edge without dipoles)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            v, &                          ! dislocation glide velocity
                                            fluxdensity, &                ! flux density at central material point
                                            neighboring_fluxdensity, &    ! flux density at neighboring material point
                                            gdot                          ! shear rates
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            rhoForest, &                  ! forest dislocation density
                                            tauThreshold, &               ! threshold shear stress
                                            tau, &                        ! current resolved shear stress
                                            tauBack, &                    ! current back stress from pileups on same slip system
                                            vClimb                        ! climb velocity of edge dipoles
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                            rhoDip, &                     ! current dipole dislocation densities (screw and edge dipoles)
                                            dLower, &                     ! minimum stable dipole distance for edges and screws
                                            dUpper                        ! current maximum stable dipole distance for edges and screws
real(pReal), dimension(3,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            m                             ! direction of dislocation motion
real(pReal), dimension(3,3) ::              my_F, &                       ! my total deformation gradient
                                            neighboring_F, &              ! total deformation gradient of my neighbor
                                            my_Fe, &                      ! my elastic deformation gradient
                                            neighboring_Fe, &             ! elastic deformation gradient of my neighbor
                                            Favg                          ! average total deformation gradient of me and my neighbor
real(pReal), dimension(3) ::                normal_neighbor2me, &         ! interface normal pointing from my neighbor to me in neighbor's lattice configuration
                                            normal_neighbor2me_defConf, & ! interface normal pointing from my neighbor to me in shared deformed configuration
                                            normal_me2neighbor, &         ! interface normal pointing from me to my neighbor in my lattice configuration
                                            normal_me2neighbor_defConf    ! interface normal pointing from me to my neighbor in shared deformed configuration
real(pReal)                                 area, &                       ! area of the current interface
                                            transmissivity, &             ! overall transmissivity of dislocation flux to neighboring material point
                                            lineLength, &                 ! dislocation line length leaving the current interface
                                            D                             ! self diffusion
logical                                     considerEnteringFlux, &
                                            considerLeavingFlux

#ifndef _OPENMP
  if (debug_verbosity > 6 .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g) .or. .not. debug_selectiveDebugger)) then
    write(6,*)
    write(6,'(a,i8,x,i2,x,i1)') '<< CONST >> nonlocal_dotState at el ip g ',el,ip,g
    write(6,*)
  endif
#endif

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
gdot = 0.0_pReal
dLower = 0.0_pReal
dUpper = 0.0_pReal


!*** shortcut to state variables 

forall (s = 1:ns, t = 1:4) &
  rhoSgl(s,t) = max(state(g,ip,el)%p((t-1)*ns+s), 0.0_pReal)
forall (s = 1:ns, t = 5:8) &
  rhoSgl(s,t) = state(g,ip,el)%p((t-1)*ns+s)
forall (s = 1:ns, c = 1:2) &
  rhoDip(s,c) = max(state(g,ip,el)%p((7+c)*ns+s), 0.0_pReal)
rhoForest = state(g,ip,el)%p(10*ns+1:11*ns)
tauThreshold = state(g,ip,el)%p(11*ns+1:12*ns)
tauBack = state(g,ip,el)%p(12*ns+1:13*ns)
forall (t = 1:4) &
  v(1:ns,t) = state(g,ip,el)%p((12+t)*ns+1:(13+t)*ns)


!*** sanity check for timestep

if (timestep <= 0.0_pReal) then                                                                                                     ! if illegal timestep...
  dotState%p = 0.0_pReal                                                                                                            ! ...return without doing anything (-> zero dotState)
  return
endif



!****************************************************************************
!*** Calculate shear rate

forall (t = 1:4) &
  gdot(1:ns,t) = rhoSgl(1:ns,t) * constitutive_nonlocal_burgers(1:ns,myInstance) * v(1:ns,t)
forall (s = 1:ns, t = 1:4, rhoSgl(s,t+4) * v(s,t) < 0.0_pReal) &                                                                    ! contribution of used rho for changing sign of v
  gdot(s,t) = gdot(s,t) + abs(rhoSgl(s,t+4)) * constitutive_nonlocal_burgers(s,myInstance) * v(s,t)

#ifndef _OPENMP
  if (debug_verbosity > 6 .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g) .or. .not. debug_selectiveDebugger)) then
    write(6,'(a,/,10(12(x),12(e12.5,x),/))') '<< CONST >> rho / 1/m^2', rhoSgl, rhoDip
    write(6,'(a,/,4(12(x),12(e12.5,x),/))') '<< CONST >> gdot / 1/s',gdot
  endif
#endif



!****************************************************************************
!*** check LFC condition for flux

if (any(abs(gdot) > 0.0_pReal .and. 2.0_pReal * v * timestep > mesh_ipVolume(ip,el) / maxval(mesh_ipArea(:,ip,el)))) then           ! safety factor 2.0 (we use the reference volume and are for simplicity here)
#ifndef _OPENMP
  if (debug_verbosity > 6) then
    write(6,'(a,i5,a,i2)') '<< CONST >> LFC condition not fullfilled at el ',el,' ip ',ip
  endif
#endif
  dotState%p = DAMASK_NaN
  return
endif



!****************************************************************************
!*** calculate limits for stable dipole height

do s = 1,ns   ! loop over slip systems
  sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)  
  tau(s) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,sLattice,myStructure)) + tauBack(s)
enddo

dLower = constitutive_nonlocal_minimumDipoleHeight(1:ns,1:2,myInstance)
dUpper(1:ns,2) = min( 1.0_pReal / sqrt( sum(abs(rhoSgl),2)+sum(rhoDip,2) ), &
                      constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgers(1:ns,myInstance) &
                                                             / ( 8.0_pReal * pi * abs(tau) ) )
dUpper(1:ns,1) = dUpper(1:ns,2) / ( 1.0_pReal - constitutive_nonlocal_nu(myInstance) )



!****************************************************************************
!*** dislocation remobilization (bauschinger effect)

rhoDotRemobilization = 0.0_pReal
if (timestep > 0.0_pReal) then
  do t = 1,4
    do s = 1,ns
      if (rhoSgl(s,t+4) * v(s,t) < 0.0_pReal) then
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
where (rhoSgl(1:ns,3:4) > 0.0_pReal) &
  rhoDotMultiplication(1:ns,1:2) = spread(0.5_pReal * sum(abs(gdot(1:ns,3:4)),2) * sqrt(rhoForest)  &
                                                    / constitutive_nonlocal_lambda0(1:ns,myInstance) &
                                                    / constitutive_nonlocal_burgers(1:ns,myInstance), 2, 2)
where (rhoSgl(1:ns,1:2) > 0.0_pReal) &
  rhoDotMultiplication(1:ns,3:4) = spread(0.5_pReal * sum(abs(gdot(1:ns,1:2)),2) * sqrt(rhoForest)  &
                                                    / constitutive_nonlocal_lambda0(1:ns,myInstance) &
                                                    / constitutive_nonlocal_burgers(1:ns,myInstance), 2, 2)



!****************************************************************************
!*** calculate dislocation fluxes (only for nonlocal constitution)

rhoDotFlux = 0.0_pReal

if (.not. phase_localConstitution(material_phase(g,ip,el))) then                                                                    ! only for nonlocal constitution
  
  !*** take care of the definition of lattice_st = lattice_sd x lattice_sn !!!
  !*** opposite sign to our p vector in the (s,p,n) triplet !!!
  
  m(1:3,1:ns,1) =  lattice_sd(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,myInstance), myStructure)
  m(1:3,1:ns,2) = -lattice_sd(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,myInstance), myStructure)
  m(1:3,1:ns,3) = -lattice_st(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,myInstance), myStructure)
  m(1:3,1:ns,4) =  lattice_st(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,myInstance), myStructure)
  
  my_Fe = Fe(1:3,1:3,g,ip,el)
  my_F = math_mul33x33(my_Fe, Fp(1:3,1:3,g,ip,el))
  
  fluxdensity = rhoSgl(1:ns,1:4) * v
    
  do n = 1,FE_NipNeighbors(mesh_element(2,el))                                                                                      ! loop through my neighbors
    neighboring_el = mesh_ipNeighborhood(1,n,ip,el)
    neighboring_ip = mesh_ipNeighborhood(2,n,ip,el)
    if (neighboring_el > 0_pInt .and. neighboring_ip > 0_pInt) then                                                                 ! if neighbor exists ...
      do neighboring_n = 1,FE_NipNeighbors(mesh_element(2,neighboring_el))                                                          ! find neighboring index that points from my neighbor to myself
        if (      el == mesh_ipNeighborhood(1,neighboring_n,neighboring_ip,neighboring_el) &
            .and. ip == mesh_ipNeighborhood(2,neighboring_n,neighboring_ip,neighboring_el)) then                                    ! possible candidate
          if (math_mul3x3(mesh_ipAreaNormal(1:3,n,ip,el),&
                          mesh_ipAreaNormal(1:3,neighboring_n,neighboring_ip,neighboring_el)) < 0.0_pReal) then                     ! area normals have opposite orientation (we have to check that because of special case for single element with two ips and periodicity. In this case the neighbor is identical in two different directions.)
            exit
          endif
        endif
      enddo
    endif
  
    opposite_n = n + mod(n,2) - mod(n+1,2)
    opposite_el = mesh_ipNeighborhood(1,opposite_n,ip,el)
    opposite_ip = mesh_ipNeighborhood(2,opposite_n,ip,el)
  
    if (neighboring_el > 0_pInt .and. neighboring_ip > 0_pInt) then                                                                 ! if neighbor exists, average deformation gradient
      neighboring_Fe = Fe(1:3,1:3,g,neighboring_ip,neighboring_el)
      neighboring_F = math_mul33x33(neighboring_Fe, Fp(1:3,1:3,g,neighboring_ip,neighboring_el))
      Favg = 0.5_pReal * (my_F + neighboring_F)
    else                                                                                                                            ! if no neighbor, take my value as average
      Favg = my_F
    endif
    

    !* FLUX FROM MY NEIGHBOR TO ME
    !* This is only considered, if I have a neighbor of nonlocal constitution (also nonlocal constitutive law with local properties) that is at least a little bit compatible.
    !* If it's not at all compatible, no flux is arriving, because everything is dammed in front of my neighbor's interface.
    !* The entering flux from my neighbor will be distributed on my slip systems according to the compatibility
    
    considerEnteringFlux = .false.
    neighboring_fluxdensity = 0.0_pReal   ! needed for check of sign change in flux density below 
    if (neighboring_el > 0_pInt .or. neighboring_ip > 0_pInt) then
      if (phase_constitution(material_phase(1,neighboring_ip,neighboring_el)) == constitutive_nonlocal_label &
          .and. any(constitutive_nonlocal_compatibility(:,:,:,n,ip,el) > 0.0_pReal)) &
        considerEnteringFlux = .true.
    endif
    
    if (considerEnteringFlux) then
      forall (t = 1:4) &
        neighboring_fluxdensity(1:ns,t) = state(g,neighboring_ip,neighboring_el)%p((t-1)*ns+1:t*ns) &
                                        * state(g,neighboring_ip,neighboring_el)%p((12+t)*ns+1:(13+t)*ns)
      normal_neighbor2me_defConf = math_det3x3(Favg) &
                  * math_mul33x3(math_inv3x3(transpose(Favg)), mesh_ipAreaNormal(1:3,neighboring_n,neighboring_ip,neighboring_el))  ! calculate the normal of the interface in (average) deformed configuration (now pointing from my neighbor to me!!!)
      normal_neighbor2me = math_mul33x3(transpose(neighboring_Fe), normal_neighbor2me_defConf) / math_det3x3(neighboring_Fe)        ! interface normal in the lattice configuration of my neighbor
      area = mesh_ipArea(neighboring_n,neighboring_ip,neighboring_el) * math_norm3(normal_neighbor2me)
      normal_neighbor2me = normal_neighbor2me / math_norm3(normal_neighbor2me)                                                      ! normalize the surface normal to unit length
      do s = 1,ns
        do t = 1,4
          c = (t + 1) / 2
          topp = t + mod(t,2) - mod(t+1,2)
          if (neighboring_fluxdensity(s,t) * math_mul3x3(m(1:3,s,t), normal_neighbor2me) > 0.0_pReal &                              ! flux from my neighbor to me == entering flux for me
              .and. fluxdensity(s,t) * neighboring_fluxdensity(s,t) >= 0.0_pReal ) then                                             ! ... only if no sign change in flux density  
            lineLength = neighboring_fluxdensity(s,t) * math_mul3x3(m(1:3,s,t), normal_neighbor2me) * area                          ! positive line length that wants to enter through this interface
            where (constitutive_nonlocal_compatibility(c,1:ns,s,n,ip,el) > 0.0_pReal) &                                             ! positive compatibility...
              rhoDotFlux(1:ns,t) = rhoDotFlux(1:ns,t) + lineLength / mesh_ipVolume(ip,el) &                                         ! ... transferring to equally signed dislocation type
                                                      * constitutive_nonlocal_compatibility(c,1:ns,s,n,ip,el) ** 2.0_pReal
            where (constitutive_nonlocal_compatibility(c,1:ns,s,n,ip,el) < 0.0_pReal) &                                             ! ..negative compatibility...
              rhoDotFlux(1:ns,topp) = rhoDotFlux(1:ns,topp) + lineLength / mesh_ipVolume(ip,el) &                                   ! ... transferring to opposite signed dislocation type
                                                      * constitutive_nonlocal_compatibility(c,1:ns,s,n,ip,el) ** 2.0_pReal
          endif
        enddo
      enddo
    endif
   
 
    !* FLUX FROM ME TO MY NEIGHBOR
    !* This is not considered, if my opposite neighbor has a different constitutive law than nonlocal (still considered for nonlocal law with lcal properties). 
    !* Then, we assume, that the opposite(!) neighbor sends an equal amount of dislocations to me.
    !* So the net flux in the direction of my neighbor is equal to zero:
    !*    leaving flux to neighbor == entering flux from opposite neighbor
    !* In case of reduced transmissivity, part of the leaving flux is stored as dead dislocation density.
    !* That means for an interface of zero transmissivity the leaving flux is fully converted to dead dislocations.
    
    considerLeavingFlux = .true.
    if (opposite_el > 0 .and. opposite_ip > 0) then
      if (phase_constitution(material_phase(1,opposite_ip,opposite_el)) /= constitutive_nonlocal_label) &
        considerLeavingFlux = .false.
    endif

    if (considerLeavingFlux) then
      normal_me2neighbor_defConf = math_det3x3(Favg) * math_mul33x3(math_inv3x3(transpose(Favg)), mesh_ipAreaNormal(1:3,n,ip,el))   ! calculate the normal of the interface in (average) deformed configuration (pointing from me to my neighbor!!!)
      normal_me2neighbor = math_mul33x3(transpose(my_Fe), normal_me2neighbor_defConf) / math_det3x3(my_Fe)                          ! interface normal in my lattice configuration
      area = mesh_ipArea(n,ip,el) * math_norm3(normal_me2neighbor)
      normal_me2neighbor = normal_me2neighbor / math_norm3(normal_me2neighbor)                                                      ! normalize the surface normal to unit length    
      do s = 1,ns
        do t = 1,4
          c = (t + 1) / 2        
          if (fluxdensity(s,t) * math_mul3x3(m(1:3,s,t), normal_me2neighbor) > 0.0_pReal ) then                                     ! flux from me to my neighbor == leaving flux for me (might also be a pure flux from my mobile density to dead density if interface not at all transmissive)
            lineLength = fluxdensity(s,t) * math_mul3x3(m(1:3,s,t), normal_me2neighbor) * area                                      ! positive line length that wants to leave through this interface
            if (fluxdensity(s,t) * neighboring_fluxdensity(s,t) >= 0.0_pReal) then                                                  ! no sign change in flux density
              transmissivity = sum(constitutive_nonlocal_compatibility(c,1:ns,s,n,ip,el)**2.0_pReal)                                ! overall transmissivity from this slip system to my neighbor
            else                                                                                                                    ! sign change in flux density means sign change in stress which does not allow for dislocations to arive at the neighbor
              transmissivity = 0.0_pReal
            endif
            rhoDotFlux(s,t) = rhoDotFlux(s,t) - lineLength / mesh_ipVolume(ip,el)                                                   ! subtract dislocation flux from current mobile type
            rhoDotFlux(s,t+4) = rhoDotFlux(s,t+4) + lineLength / mesh_ipVolume(ip,el) * (1.0_pReal - transmissivity) &
                                                               * sign(1.0_pReal, fluxdensity(s,t))                                  ! dislocation flux that is not able to leave through interface (because of low transmissivity) will remain as immobile single density at the material point
          endif
        enddo
      enddo
    endif    
    
  enddo ! neighbor loop  
endif

if (numerics_integrationMode == 1_pInt) then
  constitutive_nonlocal_rhoDotFlux(1:ns,1:10,g,ip,el) = rhoDotFlux(1:ns,1:10)                                                       ! save flux calculation for output (if in central integration mode)
endif



!****************************************************************************
!*** calculate dipole formation and annihilation

!*** formation by glide

do c = 1,2

  rhoDotSingle2DipoleGlide(1:ns,2*c-1) = -2.0_pReal * dUpper(1:ns,c) / constitutive_nonlocal_burgers(1:ns,myInstance) &
                                                    * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) &                                   ! negative mobile --> positive mobile
                                                       + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1)) &                                 ! positive mobile --> negative mobile
                                                       + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1)))                           ! positive mobile --> negative immobile

  rhoDotSingle2DipoleGlide(1:ns,2*c) = -2.0_pReal * dUpper(1:ns,c) / constitutive_nonlocal_burgers(1:ns,myInstance) &
                                                  * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) &                                     ! negative mobile --> positive mobile
                                                     + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1)) &                                   ! positive mobile --> negative mobile
                                                     + abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c)))                               ! negative mobile --> positive immobile

  rhoDotSingle2DipoleGlide(1:ns,2*c+3) = -2.0_pReal * dUpper(1:ns,c) / constitutive_nonlocal_burgers(1:ns,myInstance) &
                                                    * rhoSgl(1:ns,2*c+3) * abs(gdot(1:ns,2*c))                                      ! negative mobile --> positive immobile

  rhoDotSingle2DipoleGlide(1:ns,2*c+4) = -2.0_pReal * dUpper(1:ns,c) / constitutive_nonlocal_burgers(1:ns,myInstance) &
                                                    * rhoSgl(1:ns,2*c+4) * abs(gdot(1:ns,2*c-1))                                    ! positive mobile --> negative immobile

  rhoDotSingle2DipoleGlide(1:ns,c+8) = - rhoDotSingle2DipoleGlide(1:ns,2*c-1) - rhoDotSingle2DipoleGlide(1:ns,2*c) &
                                       + abs(rhoDotSingle2DipoleGlide(1:ns,2*c+3)) + abs(rhoDotSingle2DipoleGlide(1:ns,2*c+4))
enddo


!*** athermal annihilation

rhoDotAthermalAnnihilation = 0.0_pReal

forall (c=1:2) &  
  rhoDotAthermalAnnihilation(1:ns,c+8) = -2.0_pReal * dLower(1:ns,c) / constitutive_nonlocal_burgers(1:ns,myInstance) &
               * (  2.0_pReal * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1))) &             ! was single hitting single
                  + 2.0_pReal * (abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c)) + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1))) & ! was single hitting immobile single or was immobile single hit by single
                  + rhoDip(1:ns,c) * (abs(gdot(1:ns,2*c-1)) + abs(gdot(1:ns,2*c))))                                                 ! single knocks dipole constituent
  
  
!*** thermally activated annihilation of dipoles

rhoDotThermalAnnihilation = 0.0_pReal

D = constitutive_nonlocal_Dsd0(myInstance) * exp(-constitutive_nonlocal_Qsd(myInstance) / (kB * Temperature))

vClimb =  constitutive_nonlocal_atomicVolume(myInstance) * D / ( kB * Temperature ) &
          * constitutive_nonlocal_Gmod(myInstance) / ( 2.0_pReal * pi * (1.0_pReal-constitutive_nonlocal_nu(myInstance)) ) &
          * 2.0_pReal / ( dUpper(1:ns,1) + dLower(1:ns,1) )
          
rhoDotThermalAnnihilation(1:ns,9) = - 4.0_pReal * rhoDip(1:ns,1) * vClimb / (dUpper(1:ns,1) - dLower(1:ns,1))                       ! edge climb
rhoDotThermalAnnihilation(1:ns,10) = 0.0_pReal                                                                                      !!! cross slipping still has to be implemented !!!


!****************************************************************************
!*** assign the rates of dislocation densities to my dotState

rhoDot = 0.0_pReal
rhoDot = rhoDotFlux &
       + rhoDotMultiplication &
       + rhoDotRemobilization &
       + rhoDotSingle2DipoleGlide &
       + rhoDotAthermalAnnihilation &
       + rhoDotThermalAnnihilation 

dotState%p(1:10*ns) = dotState%p(1:10*ns) + reshape(rhoDot,(/10*ns/))

#ifndef _OPENMP
  if (debug_verbosity > 6 .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g) .or. .not. debug_selectiveDebugger)) then
    write(6,'(a,/,8(12(x),12(e12.5,x),/))') '<< CONST >> dislocation remobilization', rhoDotRemobilization(1:ns,1:8) * timestep
    write(6,'(a,/,4(12(x),12(e12.5,x),/))') '<< CONST >> dislocation multiplication', rhoDotMultiplication(1:ns,1:4) * timestep
    write(6,'(a,/,8(12(x),12(e12.5,x),/))') '<< CONST >> dislocation flux', rhoDotFlux(1:ns,1:8) * timestep
    write(6,'(a,/,10(12(x),12(e12.5,x),/))') '<< CONST >> dipole formation by glide', rhoDotSingle2DipoleGlide * timestep
    write(6,'(a,/,2(12(x),12(e12.5,x),/))') '<< CONST >> athermal dipole annihilation', &
                                            rhoDotAthermalAnnihilation(1:ns,1:2) * timestep
    write(6,'(a,/,2(12(x),12(e12.5,x),/))') '<< CONST >> thermally activated dipole annihilation', &
                                            rhoDotThermalAnnihilation(1:ns,9:10) * timestep
    write(6,'(a,/,10(12(x),12(e12.5,x),/))') '<< CONST >> total density change', rhoDot * timestep
    write(6,'(a,/,10(12(x),12(f12.7,x),/))') '<< CONST >> relative density change', &
                                            rhoDot(1:ns,1:8) * timestep / (abs(rhoSgl)+1.0e-10), &
                                            rhoDot(1:ns,9:10) * timestep / (rhoDip+1.0e-10)
    write(6,*)
  endif
#endif

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
                      FE_maxNipNeighbors, &
                      mesh_maxNips, &
                      mesh_NcpElems
use lattice, only:    lattice_sn, &
                      lattice_sd, &
                      lattice_st

implicit none

!* input variables
integer(pInt), intent(in) ::                    i, &                          ! ip index
                                                e                             ! element index
real(pReal), dimension(4,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                                orientation                   ! crystal orientation in quaternions
                                            
!* output variables

!* local variables
integer(pInt)                                   Nneighbors, &                 ! number of neighbors
                                                n, &                          ! neighbor index 
                                                neighboring_e, &              ! element index of my neighbor
                                                neighboring_i, &              ! integration point index of my neighbor
                                                my_phase, &
                                                neighboring_phase, &
                                                my_structure, &               ! lattice structure
                                                my_instance, &                ! instance of constitution
                                                ns, &                         ! number of active slip systems
                                                s1, &                         ! slip system index (me)
                                                s2                            ! slip system index (my neighbor)
real(pReal), dimension(4) ::                    absoluteMisorientation        ! absolute misorientation (without symmetry) between me and my neighbor
real(pReal), dimension(2,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(1,i,e))),&
                         constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(1,i,e))),&
                         FE_NipNeighbors(mesh_element(2,e))) :: &  
                                                compatibility                 ! compatibility for current element and ip
real(pReal), dimension(3,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(1,i,e)))) :: &  
                                                slipNormal, &
                                                slipDirection
real(pReal)                                     compatibilitySum, &
                                                thresholdValue, &
                                                nThresholdValues
logical, dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(1,i,e)))) :: & 
                                                belowThreshold


Nneighbors = FE_NipNeighbors(mesh_element(2,e))
my_phase = material_phase(1,i,e)
my_instance = phase_constitutionInstance(my_phase)
my_structure = constitutive_nonlocal_structure(my_instance)
ns = constitutive_nonlocal_totalNslip(my_instance)
slipNormal(1:3,1:ns) =    lattice_sn(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,my_instance), my_structure)
slipDirection(1:3,1:ns) = lattice_sd(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,my_instance), my_structure)


!*** start out fully compatible

compatibility = 0.0_pReal
forall(s1 = 1:ns) &
  compatibility(1:2,s1,s1,1:Nneighbors) = 1.0_pReal


!*** Loop thrugh neighbors and check whether there is any compatibility.

do n = 1,Nneighbors
  neighboring_e = mesh_ipNeighborhood(1,n,i,e)
  neighboring_i = mesh_ipNeighborhood(2,n,i,e)
  
  
  !* FREE SURFACE
  !* Set surface transmissivity to the value specified in the material.config
  
  if (neighboring_e <= 0 .or. neighboring_i <= 0) then
    forall(s1 = 1:ns) &
      compatibility(1:2,s1,s1,n) = sqrt(constitutive_nonlocal_surfaceTransmissivity(my_instance))
    cycle
  endif
  
  
  !* PHASE BOUNDARY
  !* If we encounter a different nonlocal "cpfem" phase at the neighbor, 
  !* we consider this to be a real "physical" phase boundary, so completely incompatible.
  !* If the neighboring "cpfem" phase has a local constitution, 
  !* we do not consider this to be a phase boundary, so completely compatible.
  
  neighboring_phase = material_phase(1,neighboring_i,neighboring_e)
  if (neighboring_phase /= my_phase) then
    if (.not. phase_localConstitution(neighboring_phase)) then
      forall(s1 = 1:ns) &
        compatibility(1:2,s1,s1,n) = 0.0_pReal ! = sqrt(0.0)
    endif
    cycle
  endif

    
  !* GRAIN BOUNDARY ?
  !* The compatibility value is defined as the product of the slip normal projection and the slip direction projection.
  !* Its sign is always positive for screws, for edges it has the same sign as the slip normal projection. 
  !* Since the sum for each slip system can easily exceed one (which would result in a transmissivity larger than one), 
  !* only values above or equal to a certain threshold value are considered. This threshold value is chosen, such that
  !* the number of compatible slip systems is minimized with the sum of the original compatibility values exceeding one. 
  !* Finally the smallest compatibility value is decreased until the sum is exactly equal to one. 
  !* All values below the threshold are set to zero.
  
  absoluteMisorientation = math_QuaternionDisorientation(orientation(1:4,1,i,e), &
                                                         orientation(1:4,1,neighboring_i,neighboring_e), &
                                                         0_pInt)      ! no symmetry
                                                         
  do s1 = 1,ns    ! my slip systems
    do s2 = 1,ns  ! my neighbor's slip systems
      compatibility(1,s2,s1,n) =     math_mul3x3(slipNormal(1:3,s1), math_qRot(absoluteMisorientation, slipNormal(1:3,s2))) &
                               * abs(math_mul3x3(slipDirection(1:3,s1), math_qRot(absoluteMisorientation, slipDirection(1:3,s2))))
      compatibility(2,s2,s1,n) = abs(math_mul3x3(slipNormal(1:3,s1), math_qRot(absoluteMisorientation, slipNormal(1:3,s2)))) &
                               * abs(math_mul3x3(slipDirection(1:3,s1), math_qRot(absoluteMisorientation, slipDirection(1:3,s2))))
    enddo
    
    compatibilitySum = 0.0_pReal
    belowThreshold = .true.
    do while (compatibilitySum < 1.0_pReal .and. any(belowThreshold(1:ns)))
      thresholdValue = maxval(compatibility(2,1:ns,s1,n), belowThreshold(1:ns))              ! screws always positive
      nThresholdValues = real(count(compatibility(2,1:ns,s1,n) == thresholdValue),pReal)
      where (compatibility(2,1:ns,s1,n) >= thresholdValue) &
        belowThreshold(1:ns) = .false.
      if (compatibilitySum + thresholdValue * nThresholdValues > 1.0_pReal) &
        where (abs(compatibility(1:2,1:ns,s1,n)) == thresholdValue) &
          compatibility(1:2,1:ns,s1,n) = sign((1.0_pReal - compatibilitySum) / nThresholdValues, compatibility(1:2,1:ns,s1,n))
      compatibilitySum = compatibilitySum + nThresholdValues * thresholdValue
    enddo
    where (belowThreshold(1:ns)) compatibility(1,1:ns,s1,n) = 0.0_pReal
    where (belowThreshold(1:ns)) compatibility(2,1:ns,s1,n) = 0.0_pReal
  enddo ! my slip systems cycle
enddo   ! neighbor cycle

constitutive_nonlocal_compatibility(1:2,1:ns,1:ns,1:Nneighbors,i,e) = compatibility

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
real(pReal) constitutive_nonlocal_dotTemperature            ! evolution of Temperature

!* local variables
   
constitutive_nonlocal_dotTemperature = 0.0_pReal

endfunction




!*********************************************************************
!* calculates quantities characterizing the microstructure           *
!*********************************************************************
function constitutive_nonlocal_dislocationstress(state, Fe, g, ip, el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use math,     only: math_mul33x33, &
                    math_mul33x3, &
                    math_invert3x3, &
                    math_transpose3x3, &
                    pi
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_element, &
                    mesh_node0, &
                    FE_Nips, &
                    mesh_ipCenterOfGravity, &
                    mesh_ipVolume, &
                    mesh_periodicSurface
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_localConstitution, &
                    phase_constitutionInstance

implicit none


!*** input variables
integer(pInt), intent(in) ::    g, &                          ! current grain ID
                                ip, &                         ! current integration point
                                el                            ! current element
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                Fe                            ! elastic deformation gradient
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                state                         ! microstructural state

!*** input/output variables

!*** output variables
real(pReal), dimension(3,3) ::  constitutive_nonlocal_dislocationstress

!*** local variables
integer(pInt)                   neighboring_el, &             ! element number of neighboring material point
                                neighboring_ip, &             ! integration point of neighboring material point
                                instance, &                   ! my instance of this constitution
                                neighboring_instance, &       ! instance of this constitution of neighboring material point
                                latticeStruct, &              ! my lattice structure
                                neighboring_latticeStruct, &  ! lattice structure of neighboring material point
                                phase, &
                                neighboring_phase, &
                                ns, &                         ! total number of active slip systems at my material point
                                neighboring_ns, &             ! total number of active slip systems at neighboring material point
                                c, &                          ! index of dilsocation character (edge, screw)
                                s, &                          ! slip system index
                                t, &                          ! index of dilsocation type (e+, e-, s+, s-, used e+, used e-, used s+, used s-)
                                dir, &
                                deltaX, deltaY, deltaZ, &
                                side, &
                                j
integer(pInt), dimension(2,3) :: periodicImages
real(pReal)                     nu, &                         ! poisson's ratio
                                x, y, z, &                    ! coordinates of connection vector in neighboring lattice frame
                                xsquare, ysquare, zsquare, &  ! squares of respective coordinates
                                distance, &                   ! length of connection vector
                                segmentLength, &              ! segment length of dislocations
                                lambda, &
                                R, Rsquare, Rcube, &
                                denominator, &
                                flipSign, &
                                neighboring_ipVolumeSideLength, &
                                detFe
real(pReal), dimension(3) ::    connection, &                 ! connection vector between me and my neighbor in the deformed configuration
                                connection_neighboringLattice, & ! connection vector between me and my neighbor in the lattice configuration of my neighbor
                                connection_neighboringSlip, & ! connection vector between me and my neighbor in the slip system frame of my neighbor
                                maxCoord, minCoord, &
                                meshSize, &
                                ipCoords, &
                                neighboring_ipCoords
real(pReal), dimension(3,3) ::  sigma, &                      ! dislocation stress for one slip system in neighboring material point's slip system frame
                                Tdislo_neighboringLattice, &  ! dislocation stress as 2nd Piola-Kirchhoff stress at neighboring material point
                                Tdislo, &                     ! dislocation stress as 2nd Piola-Kirchhoff stress at my material point
                                invFe, &                      ! inverse of my elastic deformation gradient
                                neighboring_invFe, &
                                neighboringLattice2myLattice  ! mapping from neighboring MPs lattice configuration to my lattice configuration
real(pReal), dimension(2,2,maxval(constitutive_nonlocal_totalNslip)) :: &
                                neighboring_rhoExcess         ! excess density at neighboring material point (edge/screw,mobile/dead,slipsystem)
real(pReal), dimension(2,maxval(constitutive_nonlocal_totalNslip)) :: &
                                rhoExcessDead
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),8) :: &
                                rhoSgl                        ! single dislocation density (edge+, edge-, screw+, screw-, used edge+, used edge-, used screw+, used screw-)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                rhoForest, &                  ! forest dislocation density
                                tauThreshold                  ! threshold shear stress
logical                         inversionError

phase = material_phase(g,ip,el)
instance = phase_constitutionInstance(phase)
latticeStruct = constitutive_nonlocal_structure(instance)
ns = constitutive_nonlocal_totalNslip(instance)



!*** get basic states

forall (s = 1:ns, t = 1:4) &
  rhoSgl(s,t) = max(state(g,ip,el)%p((t-1)*ns+s), 0.0_pReal)                                                        ! ensure positive single mobile densities
forall (t = 5:8) & 
  rhoSgl(1:ns,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)



!*** calculate the dislocation stress of the neighboring excess dislocation densities
!*** zero for material points of local constitution

constitutive_nonlocal_dislocationstress = 0.0_pReal

if (.not. phase_localConstitution(phase)) then
  call math_invert3x3(Fe(1:3,1:3,1,ip,el), invFe, detFe, inversionError)
!  if (inversionError) then
!    return
!  endif

  !* in case of periodic surfaces we have to find out how many periodic images in each direction we need
  
  do dir = 1,3
    maxCoord(dir) = maxval(mesh_node0(dir,:))
    minCoord(dir) = minval(mesh_node0(dir,:))
  enddo
  meshSize = maxCoord - minCoord
  ipCoords = mesh_ipCenterOfGravity(1:3,ip,el)
  periodicImages = 0_pInt
  do dir = 1,3
    if (mesh_periodicSurface(dir)) then
      periodicImages(1,dir) =   floor((ipCoords(dir) - constitutive_nonlocal_R(instance) - minCoord(dir)) / meshSize(dir), pInt)
      periodicImages(2,dir) = ceiling((ipCoords(dir) + constitutive_nonlocal_R(instance) - maxCoord(dir)) / meshSize(dir), pInt)
    endif
  enddo

      
  !* loop through all material points (also through their periodic images if present), 
  !* but only consider nonlocal neighbors within a certain cutoff radius R
  
  do neighboring_el = 1,mesh_NcpElems
ipLoop: do neighboring_ip = 1,FE_Nips(mesh_element(2,neighboring_el))
      neighboring_phase = material_phase(g,neighboring_ip,neighboring_el)
      if (phase_localConstitution(neighboring_phase)) then
        cycle
      endif
      neighboring_instance = phase_constitutionInstance(neighboring_phase)
      neighboring_latticeStruct = constitutive_nonlocal_structure(neighboring_instance)
      neighboring_ns = constitutive_nonlocal_totalNslip(neighboring_instance)
      call math_invert3x3(Fe(1:3,1:3,1,neighboring_ip,neighboring_el), neighboring_invFe, detFe, inversionError)
!      if (inversionError) then
!        return
!      endif
      neighboring_ipVolumeSideLength = mesh_ipVolume(neighboring_ip,neighboring_el) ** (1.0_pReal/3.0_pReal) ! reference volume used here
      forall (s = 1:neighboring_ns, c = 1:2) &
        neighboring_rhoExcess(c,1,s) = state(g,neighboring_ip,neighboring_el)%p((2*c-2)*neighboring_ns+s) &  ! positive mobiles
                                     - state(g,neighboring_ip,neighboring_el)%p((2*c-1)*neighboring_ns+s)    ! negative mobiles
      forall (s = 1:neighboring_ns, c = 1:2) &
        neighboring_rhoExcess(c,2,s) = abs(state(g,neighboring_ip,neighboring_el)%p((2*c+2)*neighboring_ns+s)) & ! positive deads
                                     - abs(state(g,neighboring_ip,neighboring_el)%p((2*c+3)*neighboring_ns+s))   ! negative deads
      nu = constitutive_nonlocal_nu(neighboring_instance)
      Tdislo_neighboringLattice = 0.0_pReal
      do deltaX = periodicImages(1,1),periodicImages(2,1)
        do deltaY = periodicImages(1,2),periodicImages(2,2)
          do deltaZ = periodicImages(1,3),periodicImages(2,3)
            
            
            !* regular case
            
            if (neighboring_el /= el .or. neighboring_ip /= ip &
                .or. deltaX /= 0_pInt .or. deltaY /= 0_pInt .or. deltaZ /= 0_pInt) then
            
              neighboring_ipCoords = mesh_ipCenterOfGravity(1:3,neighboring_ip,neighboring_el) &
                                   + (/real(deltaX,pReal), real(deltaY,pReal), real(deltaZ,pReal)/) * meshSize
              connection = neighboring_ipCoords - ipCoords
              distance = sqrt(sum(connection * connection))
              if (distance > constitutive_nonlocal_R(instance)) then
                cycle
              endif
                

              !* the segment length is the minimum of the third root of the control volume and the ip distance
              !* this ensures, that the central MP never sits on a neighboring dislocation segment
              
              connection_neighboringLattice = math_mul33x3(neighboring_invFe, connection)
              segmentLength = min(neighboring_ipVolumeSideLength, distance)
      

              !* loop through all slip systems of the neighboring material point
              !* and add up the stress contributions from egde and screw excess on these slip systems (if significant)
      
              do s = 1,neighboring_ns
                if (all(abs(neighboring_rhoExcess(:,:,s)) < constitutive_nonlocal_aTolRho(instance))) then
                  cycle                                                                             ! not significant
                endif
                
                
                !* map the connection vector from the lattice into the slip system frame
                
                connection_neighboringSlip = math_mul33x3(constitutive_nonlocal_lattice2slip(1:3,1:3,s,neighboring_instance), &
                                                          connection_neighboringLattice)
                
                
                !* edge contribution to stress
                sigma = 0.0_pReal
                
                x = connection_neighboringSlip(1)
                y = connection_neighboringSlip(2)
                z = connection_neighboringSlip(3)
                xsquare = x * x
                ysquare = y * y
                zsquare = z * z
                do j = 1,2
                  if (abs(neighboring_rhoExcess(1,j,s)) < constitutive_nonlocal_aTolRho(instance)) then
                    cycle 
                  elseif (j > 1) then
                    x = connection_neighboringSlip(1) + sign(0.5_pReal * segmentLength, &
                                                               state(g,neighboring_ip,neighboring_el)%p(4*neighboring_ns+s) &
                                                             - state(g,neighboring_ip,neighboring_el)%p(5*neighboring_ns+s))
                    xsquare = x * x
                  endif
                   
                  flipSign = sign(1.0_pReal, -y)
                  do side = 1,-1,-2
                    lambda = real(side,pReal) * 0.5_pReal * segmentLength - y
                    R = sqrt(xsquare + zsquare + lambda * lambda)
                    Rsquare = R * R
                    Rcube = Rsquare * R 
                    denominator = R * (R + flipSign * lambda)
                    if (denominator == 0.0_pReal) then
                      exit ipLoop
                    endif
                      
                    sigma(1,1) = sigma(1,1) - real(side,pReal) * flipSign * z / denominator &
                                                               * (1.0_pReal + xsquare / Rsquare + xsquare / denominator) &
                                                               * neighboring_rhoExcess(1,j,s)
                    sigma(2,2) = sigma(2,2) - real(side,pReal) * (flipSign * 2.0_pReal * nu * z / denominator + z * lambda / Rcube)&
                                                               * neighboring_rhoExcess(1,j,s)
                    sigma(3,3) = sigma(3,3) + real(side,pReal) * flipSign * z / denominator &
                                                               * (1.0_pReal - zsquare / Rsquare - zsquare / denominator) &
                                                               * neighboring_rhoExcess(1,j,s)
                    sigma(1,2) = sigma(1,2) + real(side,pReal) * x * z / Rcube * neighboring_rhoExcess(1,j,s)
                    sigma(1,3) = sigma(1,3) + real(side,pReal) * flipSign * x / denominator &
                                                               * (1.0_pReal - zsquare / Rsquare - zsquare / denominator) &
                                                               * neighboring_rhoExcess(1,j,s)
                    sigma(2,3) = sigma(2,3) - real(side,pReal) * (nu / R - zsquare / Rcube) * neighboring_rhoExcess(1,j,s)
                  enddo
                enddo 
                
                !* screw contribution to stress
                
                x = connection_neighboringSlip(1)   ! have to restore this value, because position might have been adapted for edge deads before
                do j = 1,2
                  if (abs(neighboring_rhoExcess(2,j,s)) < constitutive_nonlocal_aTolRho(instance)) then
                    cycle 
                  elseif (j > 1) then
                    y = connection_neighboringSlip(2) + sign(0.5_pReal * segmentLength, &
                                                               state(g,neighboring_ip,neighboring_el)%p(6*neighboring_ns+s) &
                                                             - state(g,neighboring_ip,neighboring_el)%p(7*neighboring_ns+s))
                    ysquare = y * y
                  endif

                  flipSign = sign(1.0_pReal, x)
                  do side = 1,-1,-2
                    lambda = x + real(side,pReal) * 0.5_pReal * segmentLength
                    R = sqrt(ysquare + zsquare + lambda * lambda)
                    Rsquare = R * R
                    Rcube = Rsquare * R 
                    denominator = R * (R + flipSign * lambda)
                    if (denominator == 0.0_pReal) then
                      exit ipLoop
                    endif
                    
                    sigma(1,2) = sigma(1,2) - real(side,pReal) * flipSign * z * (1.0_pReal - nu) / denominator &
                                                                              * neighboring_rhoExcess(2,j,s)
                    sigma(1,3) = sigma(1,3) + real(side,pReal) * flipSign * y * (1.0_pReal - nu) / denominator &
                                                                              * neighboring_rhoExcess(2,j,s)
                  enddo
                enddo
               
                if (all(abs(sigma) < 1.0e-10_pReal)) then ! SIGMA IS NOT A REAL STRESS, THATS WHY WE NEED A REALLY SMALL VALUE HERE
                  cycle
                endif

                !* copy symmetric parts
                
                sigma(2,1) = sigma(1,2)
                sigma(3,1) = sigma(1,3)
                sigma(3,2) = sigma(2,3)

                
                !* scale stresses and map them into the neighboring material point's lattice configuration
                
                sigma = sigma * constitutive_nonlocal_Gmod(neighboring_instance) &
                              * constitutive_nonlocal_burgers(s,neighboring_instance) &
                              / (4.0_pReal * pi * (1.0_pReal - nu)) &
                              * mesh_ipVolume(neighboring_ip,neighboring_el) / segmentLength      ! reference volume is used here (according to the segment length calculation)
                Tdislo_neighboringLattice = Tdislo_neighboringLattice &
                      + math_mul33x33(math_transpose3x3(constitutive_nonlocal_lattice2slip(1:3,1:3,s,neighboring_instance)), &
                        math_mul33x33(sigma, constitutive_nonlocal_lattice2slip(1:3,1:3,s,neighboring_instance)))
                                            
              enddo ! slip system loop


            !* special case of central ip volume
            !* only consider dead dislocations
            !* we assume that they all sit at a distance equal to half the third root of V
            !* in direction of the according slip direction
            
            else
              
              forall (s = 1:ns, c = 1:2) &
                rhoExcessDead(c,s) = state(g,ip,el)%p((2*c+2)*ns+s) &  ! positive deads (here we use symmetry: if this has negative sign it is treated as negative density at positive position instead of positive density at negative position)
                                   + state(g,ip,el)%p((2*c+3)*ns+s)    ! negative deads (here we use symmetry: if this has negative sign it is treated as positive density at positive position instead of negative density at negative position)
              do s = 1,ns
                if (all(abs(rhoExcessDead(:,s)) < constitutive_nonlocal_aTolRho(instance))) then
                  cycle                                                                             ! not significant
                endif
                sigma = 0.0_pReal                                                                   ! all components except for sigma13 are zero
                sigma(1,3) = - (rhoExcessDead(1,s) + rhoExcessDead(2,s) * (1.0_pReal - nu)) * neighboring_ipVolumeSideLength &
                             * constitutive_nonlocal_Gmod(instance) * constitutive_nonlocal_burgers(s,instance) &
                             / (sqrt(2.0_pReal) * pi * (1.0_pReal - nu))
                sigma(3,1) = sigma(1,3)
                
                Tdislo_neighboringLattice = Tdislo_neighboringLattice &
                                      + math_mul33x33(math_transpose3x3(constitutive_nonlocal_lattice2slip(1:3,1:3,s,instance)), &
                                                      math_mul33x33(sigma, constitutive_nonlocal_lattice2slip(1:3,1:3,s,instance)))
                                            
              enddo ! slip system loop

            endif

          enddo ! deltaZ loop
        enddo ! deltaY loop
      enddo ! deltaX loop


      !* map the stress from the neighboring MP's lattice configuration into the deformed configuration 
      !* and back into my lattice configuration

      neighboringLattice2myLattice = math_mul33x33(invFe, Fe(1:3,1:3,1,neighboring_ip,neighboring_el))
      constitutive_nonlocal_dislocationstress = constitutive_nonlocal_dislocationstress &
                                              + math_mul33x33(neighboringLattice2myLattice, &
                                                math_mul33x33(Tdislo_neighboringLattice, &
                                                math_transpose3x3(neighboringLattice2myLattice)))
                        
    enddo ipLoop
  enddo ! element loop
    
endif

endfunction


!*********************************************************************
!* return array of constitutive results                              *
!*********************************************************************
function constitutive_nonlocal_postResults(Tstar_v, Fe, Temperature, dt, state, dotState, g,ip,el)

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
                    math_transpose3x3, &
                    pi
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_element
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_constitutionInstance, &
                    phase_Noutput
use lattice,  only: lattice_Sslip, &
                    lattice_Sslip_v, &
                    lattice_sd, &
                    lattice_st, &
                    lattice_maxNslipFamily

implicit none

!*** input variables
integer(pInt), intent(in) ::                g, &                      ! current grain number
                                            ip, &                     ! current integration point
                                            el                        ! current element number
real(pReal), intent(in) ::                  Temperature, &            ! temperature
                                            dt                        ! time increment
real(pReal), dimension(6), intent(in) ::    Tstar_v                   ! current 2nd Piola-Kirchhoff stress in Mandel notation
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) ::  &
                                            Fe                        ! elastic deformation gradient
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            state                     ! current microstructural state
type(p_vec), intent(in) ::                  dotState                  ! evolution rate of microstructural state

!*** output variables
real(pReal), dimension(constitutive_nonlocal_sizePostResults(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            constitutive_nonlocal_postResults

!*** local variables
integer(pInt)                               myInstance, &             ! current instance of this constitution
                                            myStructure, &            ! current lattice structure
                                            ns, &                     ! short notation for the total number of active slip systems
                                            c, &                      ! character of dislocation
                                            cs, &                     ! constitutive result index
                                            o, &                      ! index of current output
                                            t, &                      ! type of dislocation
                                            s, &                      ! index of my current slip system
                                            sLattice                  ! index of my current slip system according to lattice order
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),8) :: &
                                            rhoSgl, &                 ! current single dislocation densities (positive/negative screw and edge without dipoles)
                                            rhoDotSgl                 ! evolution rate of single dislocation densities (positive/negative screw and edge without dipoles)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            gdot, &                   ! shear rates
                                            v                         ! velocities
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            rhoForest, &              ! forest dislocation density
                                            tauThreshold, &           ! threshold shear stress
                                            tau, &                    ! current resolved shear stress
                                            tauBack, &                ! back stress from pileups on same slip system
                                            vClimb                    ! climb velocity of edge dipoles
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                            rhoDip, &                 ! current dipole dislocation densities (screw and edge dipoles)
                                            rhoDotDip, &              ! evolution rate of dipole dislocation densities (screw and edge dipoles)
                                            dLower, &                 ! minimum stable dipole distance for edges and screws
                                            dUpper                    ! current maximum stable dipole distance for edges and screws
real(pReal), dimension(3,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                            m, &                      ! direction of dislocation motion for edge and screw (unit vector)
                                            m_currentconf             ! direction of dislocation motion for edge and screw (unit vector) in current configuration
real(pReal)                                 D                         ! self diffusion
real(pReal), dimension(3,3) ::              sigma

myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)

cs = 0_pInt
constitutive_nonlocal_postResults = 0.0_pReal


!* short hand notations for state variables

forall (t = 1:8) rhoSgl(1:ns,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (c = 1:2) rhoDip(1:ns,c) = state(g,ip,el)%p((7+c)*ns+1:(8+c)*ns)
rhoForest = state(g,ip,el)%p(10*ns+1:11*ns)
tauThreshold = state(g,ip,el)%p(11*ns+1:12*ns)
tauBack = state(g,ip,el)%p(12*ns+1:13*ns)
forall (t = 1:8) rhoDotSgl(1:ns,t) = dotState%p((t-1)*ns+1:t*ns)
forall (c = 1:2) rhoDotDip(1:ns,c) = dotState%p((7+c)*ns+1:(8+c)*ns)
forall (t = 1:4) v(1:ns,t) = state(g,ip,el)%p((12+t)*ns+1:(13+t)*ns)


!* Calculate shear rate

do t = 1,4
  do s = 1,ns
    if (rhoSgl(s,t+4) * v(s,t) < 0.0_pReal) then
      rhoSgl(s,t) = rhoSgl(s,t) + abs(rhoSgl(s,t+4))                                                                  ! remobilization of immobile singles for changing sign of v (bauschinger effect)
      rhoSgl(s,t+4) = 0.0_pReal                                                                                       ! remobilization of immobile singles for changing sign of v (bauschinger effect)
    endif
  enddo
enddo

forall (t = 1:4) &
  gdot(1:ns,t) = rhoSgl(1:ns,t) * constitutive_nonlocal_burgers(1:ns,myInstance) * v(1:ns,t)
  

!* calculate limits for stable dipole height

do s = 1,ns
  sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
  tau(s) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,sLattice,myStructure)) + tauBack(s)
enddo

dLower = constitutive_nonlocal_minimumDipoleHeight(1:ns,1:2,myInstance)
dUpper(1:ns,2) = min( constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgers(1:ns,myInstance) &
                                                             / (8.0_pReal * pi * abs(tau)), &
                      1.0_pReal / sqrt(sum(abs(rhoSgl),2)+sum(rhoDip,2)) )
dUpper(1:ns,1) = dUpper(1:ns,2) / (1.0_pReal - constitutive_nonlocal_nu(myInstance))


!*** dislocation motion

m(1:3,1:ns,1) = lattice_sd(1:3,constitutive_nonlocal_slipSystemLattice(1:ns,myInstance),myStructure)
m(1:3,1:ns,2) = -lattice_st(1:3,constitutive_nonlocal_slipSystemLattice(1:ns,myInstance),myStructure)
forall (c = 1:2, s = 1:ns) &
  m_currentconf(1:3,s,c) = math_mul33x3(Fe, m(1:3,s,c))


do o = 1,phase_Noutput(material_phase(g,ip,el))
  select case(constitutive_nonlocal_output(o,myInstance))
    
    case ('rho')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl),2) + sum(rhoDip,2)
      cs = cs + ns
      
    case ('rho_sgl')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl),2)
      cs = cs + ns
      
    case ('rho_sgl_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(1:ns,1:4)),2)
      cs = cs + ns
      
    case ('rho_sgl_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoSgl(1:ns,5:8),2)
      cs = cs + ns
      
    case ('rho_dip')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoDip,2)
      cs = cs + ns
      
    case ('rho_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(1:ns,(/1,2,5,6/))),2) + rhoDip(1:ns,1)
      cs = cs + ns
      
    case ('rho_sgl_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(1:ns,(/1,2,5,6/))),2)
      cs = cs + ns
      
    case ('rho_sgl_edge_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoSgl(1:ns,1:2),2)
      cs = cs + ns
      
    case ('rho_sgl_edge_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoSgl(1:ns,5:6),2)
      cs = cs + ns
      
    case ('rho_sgl_edge_pos')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,1) + abs(rhoSgl(1:ns,5))
      cs = cs + ns
      
    case ('rho_sgl_edge_pos_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,1)
      cs = cs + ns
      
    case ('rho_sgl_edge_pos_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,5)
      cs = cs + ns
      
    case ('rho_sgl_edge_neg')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,2) + abs(rhoSgl(1:ns,6))
      cs = cs + ns
      
    case ('rho_sgl_edge_neg_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,2)
      cs = cs + ns
      
    case ('rho_sgl_edge_neg_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,6)
      cs = cs + ns
      
    case ('rho_dip_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoDip(1:ns,1)
      cs = cs + ns
      
    case ('rho_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(1:ns,(/3,4,7,8/))),2) + rhoDip(1:ns,2)
      cs = cs + ns
      
    case ('rho_sgl_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(1:ns,(/3,4,7,8/))),2)
      cs = cs + ns
            
    case ('rho_sgl_screw_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoSgl(1:ns,3:4),2)
      cs = cs + ns
      
    case ('rho_sgl_screw_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoSgl(1:ns,7:8),2)
      cs = cs + ns
      
    case ('rho_sgl_screw_pos')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,3) + abs(rhoSgl(1:ns,7))
      cs = cs + ns
      
    case ('rho_sgl_screw_pos_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,3)
      cs = cs + ns
      
    case ('rho_sgl_screw_pos_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,7)
      cs = cs + ns
      
    case ('rho_sgl_screw_neg')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,4) + abs(rhoSgl(1:ns,8))
      cs = cs + ns

    case ('rho_sgl_screw_neg_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,4)
      cs = cs + ns

    case ('rho_sgl_screw_neg_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,8)
      cs = cs + ns

    case ('rho_dip_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoDip(1:ns,2)
      cs = cs + ns
      
    case ('excess_rho')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = (rhoSgl(1:ns,1) + abs(rhoSgl(1:ns,5))) &
                                                    - (rhoSgl(1:ns,2) + abs(rhoSgl(1:ns,6))) &
                                                    + (rhoSgl(1:ns,3) + abs(rhoSgl(1:ns,7))) &
                                                    - (rhoSgl(1:ns,4) + abs(rhoSgl(1:ns,8)))
      cs = cs + ns
      
    case ('excess_rho_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = (rhoSgl(1:ns,1) + abs(rhoSgl(1:ns,5))) &
                                                    - (rhoSgl(1:ns,2) + abs(rhoSgl(1:ns,6)))
      cs = cs + ns
      
    case ('excess_rho_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = (rhoSgl(1:ns,3) + abs(rhoSgl(1:ns,7))) &
                                                    - (rhoSgl(1:ns,4) + abs(rhoSgl(1:ns,8)))
      cs = cs + ns
      
    case ('rho_forest')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoForest
      cs = cs + ns
    
    case ('delta')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = 1.0_pReal / sqrt(sum(abs(rhoSgl),2) + sum(rhoDip,2))
      cs = cs + ns
      
    case ('delta_sgl')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = 1.0_pReal / sqrt(sum(abs(rhoSgl),2))
      cs = cs + ns
      
    case ('delta_dip')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = 1.0_pReal / sqrt(sum(rhoDip,2))
      cs = cs + ns
      
    case ('shearrate')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(gdot,2)
      cs = cs + ns
      
    case ('resolvedstress')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = tau
      cs = cs + ns
      
    case ('resolvedstress_back')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = tauBack
      cs = cs + ns
      
    case ('resolvedstress_external')
      do s = 1,ns  
        sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
        constitutive_nonlocal_postResults(cs+s) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,sLattice,myStructure))
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
                                                      / constitutive_nonlocal_lambda0(1:ns,myInstance) &
                                                      / constitutive_nonlocal_burgers(1:ns,myInstance)
      cs = cs + ns

    case ('rho_dot_gen_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) =   sum(abs(gdot(1:ns,3:4)),2) * sqrt(rhoForest)  &
                                                      / constitutive_nonlocal_lambda0(1:ns,myInstance) &
                                                      / constitutive_nonlocal_burgers(1:ns,myInstance)
      cs = cs + ns

    case ('rho_dot_gen_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) =   sum(abs(gdot(1:ns,1:2)),2) * sqrt(rhoForest)  &
                                                      / constitutive_nonlocal_lambda0(1:ns,myInstance) &
                                                      / constitutive_nonlocal_burgers(1:ns,myInstance)
      cs = cs + ns
      
    case ('rho_dot_sgl2dip')
      do c=1,2                                                                                                                      ! dipole formation by glide
        constitutive_nonlocal_postResults(cs+1:cs+ns) = constitutive_nonlocal_postResults(cs+1:cs+ns) + &
            2.0_pReal * dUpper(1:ns,c) / constitutive_nonlocal_burgers(1:ns,myInstance) &
                      * (  2.0_pReal * (  rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) &
                                        + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1))) &                                               ! was single hitting single
                         + 2.0_pReal * (  abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c)) &
                                        + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1))))                                         ! was single hitting immobile/used single
      enddo
      cs = cs + ns
    
    case ('rho_dot_ann_ath')
      do c=1,2
        constitutive_nonlocal_postResults(cs+1:cs+ns) = constitutive_nonlocal_postResults(cs+1:cs+ns) + &
            2.0_pReal * dLower(1:ns,c) / constitutive_nonlocal_burgers(1:ns,myInstance) &
                      * (  2.0_pReal * (  rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) &
                                        + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1))) &                                               ! was single hitting single
                         + 2.0_pReal * (  abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c)) &
                                        + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1))) &                                        ! was single hitting immobile/used single
                         + rhoDip(1:ns,c) * (abs(gdot(1:ns,2*c-1)) + abs(gdot(1:ns,2*c))))                                          ! single knocks dipole constituent
      enddo
      cs = cs + ns
      
    case ('rho_dot_ann_the') 
      D = constitutive_nonlocal_Dsd0(myInstance) * exp(-constitutive_nonlocal_Qsd(myInstance) / (kB * Temperature))

      vClimb =  constitutive_nonlocal_atomicVolume(myInstance) * D / (kB * Temperature) &
          * constitutive_nonlocal_Gmod(myInstance) / (2.0_pReal * pi * (1.0_pReal-constitutive_nonlocal_nu(myInstance))) &
          * 2.0_pReal / (dUpper(1:ns,1) + dLower(1:ns,1))
          
      constitutive_nonlocal_postResults(cs+1:cs+ns) = 4.0_pReal * rhoDip(1:ns,1) * vClimb / (dUpper(1:ns,1) - dLower(1:ns,1))
      ! !!! cross-slip of screws missing !!!
      cs = cs + ns

    case ('rho_dot_flux')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(constitutive_nonlocal_rhoDotFlux(1:ns,1:4,g,ip,el),2) &
                                                      + sum(abs(constitutive_nonlocal_rhoDotFlux(1:ns,5:8,g,ip,el)),2)
      cs = cs + ns
    
    case ('rho_dot_flux_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(constitutive_nonlocal_rhoDotFlux(1:ns,1:2,g,ip,el),2) &
                                                      + sum(abs(constitutive_nonlocal_rhoDotFlux(1:ns,5:6,g,ip,el)),2)
      cs = cs + ns
      
    case ('rho_dot_flux_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(constitutive_nonlocal_rhoDotFlux(1:ns,3:4,g,ip,el),2) &
                                                      + sum(abs(constitutive_nonlocal_rhoDotFlux(1:ns,7:8,g,ip,el)),2)
      cs = cs + ns
            
    case ('velocity_edge_pos')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = v(1:ns,1)
      cs = cs + ns
    
    case ('velocity_edge_neg')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = v(1:ns,2)
      cs = cs + ns
    
    case ('velocity_screw_pos')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = v(1:ns,3)
      cs = cs + ns
    
    case ('velocity_screw_neg')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = v(1:ns,4)
      cs = cs + ns
    
    case ('fluxdensity_edge_pos_x')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,1) * v(1:ns,1) * m_currentconf(1,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_pos_y')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,1) * v(1:ns,1) * m_currentconf(2,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_pos_z')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,1) * v(1:ns,1) * m_currentconf(3,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_neg_x')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = - rhoSgl(1:ns,2) * v(1:ns,2) * m_currentconf(1,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_neg_y')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = - rhoSgl(1:ns,2) * v(1:ns,2) * m_currentconf(2,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_neg_z')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = - rhoSgl(1:ns,2) * v(1:ns,2) * m_currentconf(3,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_screw_pos_x')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,3) * v(1:ns,3) * m_currentconf(1,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_pos_y')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,3) * v(1:ns,3) * m_currentconf(2,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_pos_z')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,3) * v(1:ns,3) * m_currentconf(3,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_neg_x')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = - rhoSgl(1:ns,4) * v(1:ns,4) * m_currentconf(1,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_neg_y')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = - rhoSgl(1:ns,4) * v(1:ns,4) * m_currentconf(2,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_neg_z')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = - rhoSgl(1:ns,4) * v(1:ns,4) * m_currentconf(3,1:ns,2)
      cs = cs + ns
    
    case ('maximumdipoleheight_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = dUpper(1:ns,1)
      cs = cs + ns
      
    case ('maximumdipoleheight_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = dUpper(1:ns,2)
      cs = cs + ns
    
    case('dislocationstress')
      sigma = constitutive_nonlocal_dislocationstress(state, Fe, g, ip, el)
      constitutive_nonlocal_postResults(cs+1) = sigma(1,1)
      constitutive_nonlocal_postResults(cs+2) = sigma(2,2)
      constitutive_nonlocal_postResults(cs+3) = sigma(3,3)
      constitutive_nonlocal_postResults(cs+4) = sigma(1,2)
      constitutive_nonlocal_postResults(cs+5) = sigma(2,3)
      constitutive_nonlocal_postResults(cs+6) = sigma(3,1)
      cs = cs + 6_pInt
    
    case('accumulatedshear')
      constitutive_nonlocal_accumulatedShear(1:ns,g,ip,el) = constitutive_nonlocal_accumulatedShear(1:ns,g,ip,el) + sum(gdot,2)*dt
      constitutive_nonlocal_postResults(cs+1:cs+ns) = constitutive_nonlocal_accumulatedShear(1:ns,g,ip,el)
      cs = cs + ns

 end select
enddo

endfunction

END MODULE
