! Copyright 2011-13 Max-Planck-Institut für Eisenforschung GmbH
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
use prec, only: &
  pReal, &
  pInt, &
  p_vec

implicit none
private


!* Definition of parameters

character (len=*), parameter, public :: &
CONSTITUTIVE_NONLOCAL_LABEL = 'nonlocal'

character(len=22), dimension(11), parameter, private :: &
BASICSTATES = (/'rhoSglEdgePosMobile   ', &
                'rhoSglEdgeNegMobile   ', &
                'rhoSglScrewPosMobile  ', &
                'rhoSglScrewNegMobile  ', &
                'rhoSglEdgePosImmobile ', &
                'rhoSglEdgeNegImmobile ', &
                'rhoSglScrewPosImmobile', &
                'rhoSglScrewNegImmobile', &
                'rhoDipEdge            ', &
                'rhoDipScrew           ', &
                'accumulatedshear      ' /)                          !< list of "basic" microstructural state variables that are independent from other state variables

character(len=16), dimension(3), parameter, private :: &
DEPENDENTSTATES = (/'rhoForest       ', &
                    'tauThreshold    ', &
                    'tauBack         ' /)                            !< list of microstructural state variables that depend on other state variables

character(len=20), dimension(6), parameter, private :: &
OTHERSTATES = (/'velocityEdgePos     ', &
                'velocityEdgeNeg     ', &
                'velocityScrewPos    ', &
                'velocityScrewNeg    ', &
                'maxDipoleHeightEdge ', &
                'maxDipoleHeightScrew'  /)                           !< list of other dependent state variables that are not updated by microstructure

real(pReal), parameter, private :: &
KB = 1.38e-23_pReal                                                  !< Physical parameter, Boltzmann constant in J/Kelvin


!* Definition of global variables

integer(pInt), dimension(:), allocatable, public :: &
constitutive_nonlocal_sizeDotState, &                                !< number of dotStates = number of basic state variables
constitutive_nonlocal_sizeDependentState, &                          !< number of dependent state variables
constitutive_nonlocal_sizeState, &                                   !< total number of state variables
constitutive_nonlocal_sizePostResults                                !< cumulative size of post results

integer(pInt), dimension(:,:), allocatable, target, public :: &
constitutive_nonlocal_sizePostResult                                 !< size of each post result output

character(len=64), dimension(:,:), allocatable, target, public :: &
constitutive_nonlocal_output                                         !< name of each post result output

integer(pInt), dimension(:), allocatable, private :: &
Noutput                                                              !< number of outputs per instance of this plasticity 

integer(pInt), dimension(:,:), allocatable, private :: &
iGamma, &                                                            !< state indices for accumulated shear
iRhoF, &                                                             !< state indices for forest density
iTauF, &                                                             !< state indices for critical resolved shear stress
iTauB                                                                !< state indices for backstress
integer(pInt), dimension(:,:,:), allocatable, private :: &
iRhoU, &                                                             !< state indices for unblocked density
iRhoB, &                                                             !< state indices for blocked density
iRhoD, &                                                             !< state indices for dipole density
iV, &                                                                !< state indices for dislcation velocities
iD                                                                   !< state indices for stable dipole height


character(len=32), dimension(:), allocatable, public :: &
constitutive_nonlocal_structureName                                  !< name of the lattice structure

integer(pInt), dimension(:), allocatable, public :: &
constitutive_nonlocal_structure                                      !< number representing the kind of lattice structure

integer(pInt), dimension(:), allocatable, private :: &
totalNslip                                                           !< total number of active slip systems for each instance

integer(pInt), dimension(:,:), allocatable, private :: &
Nslip, &                                                             !< number of active slip systems for each family and instance
slipFamily, &                                                        !< lookup table relating active slip system to slip family for each instance
slipSystemLattice, &                                                 !< lookup table relating active slip system index to lattice slip system index for each instance
colinearSystem                                                       !< colinear system to the active slip system (only valid for fcc!)

real(pReal), dimension(:), allocatable, private :: &
CoverA, &                                                            !< c/a ratio for hex type lattice
mu, &                                                                !< shear modulus
nu, &                                                                !< poisson's ratio
atomicVolume, &                                                      !< atomic volume
Dsd0, &                                                              !< prefactor for self-diffusion coefficient
selfDiffusionEnergy, &                                               !< activation enthalpy for diffusion
aTolRho, &                                                           !< absolute tolerance for dislocation density in state integration
aTolShear, &                                                         !< absolute tolerance for accumulated shear in state integration
significantRho, &                                                    !< density considered significant
significantN, &                                                      !< number of dislocations considered significant
cutoffRadius, &                                                      !< cutoff radius for dislocation stress
doublekinkwidth, &                                                   !< width of a doubkle kink in multiples of the burgers vector length b
solidSolutionEnergy, &                                               !< activation energy for solid solution in J
solidSolutionSize, &                                                 !< solid solution obstacle size in multiples of the burgers vector length
solidSolutionConcentration, &                                        !< concentration of solid solution in atomic parts
pParam, &                                                            !< parameter for kinetic law (Kocks,Argon,Ashby)
qParam, &                                                            !< parameter for kinetic law (Kocks,Argon,Ashby)
viscosity, &                                                         !< viscosity for dislocation glide in Pa s
fattack, &                                                           !< attack frequency in Hz
rhoSglScatter, &                                                     !< standard deviation of scatter in initial dislocation density
surfaceTransmissivity, &                                             !< transmissivity at free surface
grainboundaryTransmissivity, &                                       !< transmissivity at grain boundary (identified by different texture)
CFLfactor, &                                                         !< safety factor for CFL flux condition
fEdgeMultiplication, &                                               !< factor that determines how much edge dislocations contribute to multiplication (0...1)
rhoSglRandom, &
rhoSglRandomBinning, &
linetensionEffect, &
edgeJogFactor

real(pReal), dimension(:,:), allocatable, private :: &
rhoSglEdgePos0, &                                                    !< initial edge_pos dislocation density per slip system for each family and instance
rhoSglEdgeNeg0, &                                                    !< initial edge_neg dislocation density per slip system for each family and instance
rhoSglScrewPos0, &                                                   !< initial screw_pos dislocation density per slip system for each family and instance
rhoSglScrewNeg0, &                                                   !< initial screw_neg dislocation density per slip system for each family and instance
rhoDipEdge0, &                                                       !< initial edge dipole dislocation density per slip system for each family and instance
rhoDipScrew0, &                                                      !< initial screw dipole dislocation density per slip system for each family and instance
lambda0PerSlipFamily, &                                              !< mean free path prefactor for each family and instance
lambda0, &                                                           !< mean free path prefactor for each slip system and instance
burgersPerSlipFamily, &                                              !< absolute length of burgers vector [m] for each family and instance
burgers, &                                                           !< absolute length of burgers vector [m] for each slip system and instance
interactionSlipSlip                                                  !< coefficients for slip-slip interaction for each interaction type and instance

real(pReal), dimension(:,:,:), allocatable, private :: &
Cslip66, &                                                           !< elasticity matrix in Mandel notation for each instance
minDipoleHeightPerSlipFamily, &                                      !< minimum stable edge/screw dipole height for each family and instance
minDipoleHeight, &                                                   !< minimum stable edge/screw dipole height for each slip system and instance
peierlsStressPerSlipFamily, &                                        !< Peierls stress (edge and screw) 
peierlsStress, &                                                     !< Peierls stress (edge and screw) 
forestProjectionEdge, &                                              !< matrix of forest projections of edge dislocations for each instance
forestProjectionScrew, &                                             !< matrix of forest projections of screw dislocations for each instance
interactionMatrixSlipSlip                                            !< interaction matrix of the different slip systems for each instance

real(pReal), dimension(:,:,:,:), allocatable, private :: &
lattice2slip, &                                                      !< orthogonal transformation matrix from lattice coordinate system to slip coordinate system (passive rotation !!!)
rhoDotEdgeJogsOutput, &
sourceProbability

real(pReal), dimension(:,:,:,:,:), allocatable, private :: &
Cslip3333, &                                                         !< elasticity matrix for each instance
rhoDotFluxOutput, & 
rhoDotMultiplicationOutput, &
rhoDotSingle2DipoleGlideOutput, &
rhoDotAthermalAnnihilationOutput, &
rhoDotThermalAnnihilationOutput

real(pReal), dimension(:,:,:,:,:,:), allocatable, private :: &
compatibility                                                        !< slip system compatibility between me and my neighbors

real(pReal), dimension(:,:), allocatable, private :: &
nonSchmidCoeff

logical, dimension(:), allocatable, private :: &
shortRangeStressCorrection, &                                        !< flag indicating the use of the short range stress correction by a excess density gradient term
deadZoneScaling, &
probabilisticMultiplication

public :: &
constitutive_nonlocal_init, &
constitutive_nonlocal_stateInit, &
constitutive_nonlocal_aTolState, &
constitutive_nonlocal_homogenizedC, &
constitutive_nonlocal_microstructure, &
constitutive_nonlocal_LpAndItsTangent, &
constitutive_nonlocal_dotState, &
constitutive_nonlocal_deltaState, &
constitutive_nonlocal_dotTemperature, &
constitutive_nonlocal_updateCompatibility, &
constitutive_nonlocal_postResults

private :: &
constitutive_nonlocal_kinetics, &
constitutive_nonlocal_dislocationstress


CONTAINS

!**************************************
!*      Module initialization         *
!**************************************
subroutine constitutive_nonlocal_init(myFile)

use, intrinsic :: iso_fortran_env                                          ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
use math,     only: math_Mandel3333to66, & 
                    math_Voigt66to3333, & 
                    math_mul3x3, &
                    math_transpose33
use IO,       only: IO_read, &
                    IO_lc, &
                    IO_getTag, &
                    IO_isBlank, &
                    IO_stringPos, &
                    IO_stringValue, &
                    IO_floatValue, &
                    IO_intValue, &
                    IO_error, &
                    IO_timeStamp
use debug,    only: debug_level, &
                    debug_constitutive, &
                    debug_levelBasic
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_maxNipNeighbors
use material, only: homogenization_maxNgrains, &
                    phase_plasticity, &
                    phase_plasticityInstance, &
                    phase_Noutput
use lattice

!*** output variables

!*** input variables
integer(pInt), intent(in) ::                myFile

!*** local variables
integer(pInt), parameter ::                 maxNchunks = 21_pInt
integer(pInt), &
    dimension(1_pInt+2_pInt*maxNchunks) ::  positions
integer(pInt), dimension(6) ::              configNchunks
integer(pInt)          ::                   section = 0_pInt, &
                                            maxNinstance, &
                                            maxTotalNslip, &
                                            myStructure, &
                                            f, &                ! index of my slip family
                                            i, &                ! index of my instance of this plasticity
                                            l, &
                                            ns, &               ! short notation for total number of active slip systems for the current instance
                                            o, &                ! index of my output
                                            s, &                ! index of my slip system
                                            s1, &               ! index of my slip system
                                            s2, &               ! index of my slip system
                                            it, &               ! index of my interaction type
                                            t, &                ! index of dislocation type
                                            c, &                ! index of dislocation character
                                            Nchunks_SlipSlip = 0_pInt, &
                                            Nchunks_SlipFamilies = 0_pInt, &
                                            mySize = 0_pInt     ! to suppress warnings, safe as init is called only once
character(len=65536)                        tag
character(len=65536) ::                     line = ''                                ! to start initialized
 
 write(6,*)
 write(6,*) '<<<+-  constitutive_',trim(CONSTITUTIVE_NONLOCAL_LABEL),' init  -+>>>'
 write(6,*) '$Id$'
 write(6,'(a16,a)')   ' Current time : ',IO_timeStamp()
#include "compilation_info.f90"

maxNinstance = int(count(phase_plasticity == CONSTITUTIVE_NONLOCAL_LABEL),pInt)
if (maxNinstance == 0) return                                                                                                       ! we don't have to do anything if there's no instance for this constitutive law

if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) then
  write(6,'(a16,1x,i5)') '# instances:',maxNinstance
  write(6,*)
endif

!*** memory allocation for global variables

allocate(constitutive_nonlocal_sizeDotState(maxNinstance))
allocate(constitutive_nonlocal_sizeDependentState(maxNinstance))
allocate(constitutive_nonlocal_sizeState(maxNinstance))
allocate(constitutive_nonlocal_sizePostResults(maxNinstance))
allocate(constitutive_nonlocal_sizePostResult(maxval(phase_Noutput), maxNinstance))
allocate(constitutive_nonlocal_output(maxval(phase_Noutput), maxNinstance))
allocate(Noutput(maxNinstance))
constitutive_nonlocal_sizeDotState = 0_pInt
constitutive_nonlocal_sizeDependentState = 0_pInt
constitutive_nonlocal_sizeState = 0_pInt
constitutive_nonlocal_sizePostResults = 0_pInt
constitutive_nonlocal_sizePostResult = 0_pInt
constitutive_nonlocal_output = ''
Noutput = 0_pInt

allocate(constitutive_nonlocal_structureName(maxNinstance))
allocate(constitutive_nonlocal_structure(maxNinstance))
allocate(Nslip(lattice_maxNslipFamily, maxNinstance))
allocate(slipFamily(lattice_maxNslip, maxNinstance))
allocate(slipSystemLattice(lattice_maxNslip, maxNinstance))
allocate(totalNslip(maxNinstance))
constitutive_nonlocal_structureName = ''
constitutive_nonlocal_structure = 0_pInt
Nslip = 0_pInt
slipFamily = 0_pInt
slipSystemLattice = 0_pInt
totalNslip = 0_pInt

allocate(CoverA(maxNinstance))
allocate(mu(maxNinstance))
allocate(nu(maxNinstance))
allocate(atomicVolume(maxNinstance))
allocate(Dsd0(maxNinstance))
allocate(selfDiffusionEnergy(maxNinstance))
allocate(aTolRho(maxNinstance))
allocate(aTolShear(maxNinstance))
allocate(significantRho(maxNinstance))
allocate(significantN(maxNinstance))
allocate(Cslip66(6,6,maxNinstance))
allocate(Cslip3333(3,3,3,3,maxNinstance))
allocate(cutoffRadius(maxNinstance))
allocate(doublekinkwidth(maxNinstance))
allocate(solidSolutionEnergy(maxNinstance))
allocate(solidSolutionSize(maxNinstance))
allocate(solidSolutionConcentration(maxNinstance))
allocate(pParam(maxNinstance))
allocate(qParam(maxNinstance))
allocate(viscosity(maxNinstance))
allocate(fattack(maxNinstance))
allocate(rhoSglScatter(maxNinstance))
allocate(rhoSglRandom(maxNinstance))
allocate(rhoSglRandomBinning(maxNinstance))
allocate(surfaceTransmissivity(maxNinstance))
allocate(grainboundaryTransmissivity(maxNinstance))
allocate(shortRangeStressCorrection(maxNinstance))
allocate(deadZoneScaling(maxNinstance))
allocate(probabilisticMultiplication(maxNinstance))
allocate(CFLfactor(maxNinstance))
allocate(fEdgeMultiplication(maxNinstance))
allocate(linetensionEffect(maxNinstance))
allocate(edgeJogFactor(maxNinstance))
CoverA = 0.0_pReal 
mu = 0.0_pReal
atomicVolume = 0.0_pReal
Dsd0 = -1.0_pReal
selfDiffusionEnergy = 0.0_pReal
aTolRho = 0.0_pReal
aTolShear = 0.0_pReal
significantRho = 0.0_pReal
significantN = 0.0_pReal
nu = 0.0_pReal
Cslip66 = 0.0_pReal
Cslip3333 = 0.0_pReal
cutoffRadius = -1.0_pReal
doublekinkwidth = 0.0_pReal
solidSolutionEnergy = 0.0_pReal
solidSolutionSize = 0.0_pReal
solidSolutionConcentration = 0.0_pReal
pParam = 1.0_pReal
qParam = 1.0_pReal
viscosity = 0.0_pReal
fattack = 0.0_pReal
rhoSglScatter = 0.0_pReal
rhoSglRandom = 0.0_pReal
rhoSglRandomBinning = 1.0_pReal
surfaceTransmissivity = 1.0_pReal
grainboundaryTransmissivity = -1.0_pReal
CFLfactor = 2.0_pReal
fEdgeMultiplication = 0.0_pReal
linetensionEffect = 0.0_pReal
edgeJogFactor = 1.0_pReal
shortRangeStressCorrection = .false.
deadZoneScaling = .false.
probabilisticMultiplication = .false.

allocate(rhoSglEdgePos0(lattice_maxNslipFamily,maxNinstance))
allocate(rhoSglEdgeNeg0(lattice_maxNslipFamily,maxNinstance))
allocate(rhoSglScrewPos0(lattice_maxNslipFamily,maxNinstance))
allocate(rhoSglScrewNeg0(lattice_maxNslipFamily,maxNinstance))
allocate(rhoDipEdge0(lattice_maxNslipFamily,maxNinstance))
allocate(rhoDipScrew0(lattice_maxNslipFamily,maxNinstance))
allocate(burgersPerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(lambda0PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(interactionSlipSlip(lattice_maxNinteraction,maxNinstance))
rhoSglEdgePos0 = -1.0_pReal
rhoSglEdgeNeg0 = -1.0_pReal
rhoSglScrewPos0 = -1.0_pReal
rhoSglScrewNeg0 = -1.0_pReal
rhoDipEdge0 = -1.0_pReal
rhoDipScrew0 = -1.0_pReal
burgersPerSlipFamily = 0.0_pReal
lambda0PerSlipFamily = 0.0_pReal
interactionSlipSlip = 0.0_pReal

allocate(minDipoleHeightPerSlipFamily(lattice_maxNslipFamily,2,maxNinstance))
allocate(peierlsStressPerSlipFamily(lattice_maxNslipFamily,2,maxNinstance))
minDipoleHeightPerSlipFamily = -1.0_pReal
peierlsStressPerSlipFamily = 0.0_pReal

allocate(nonSchmidCoeff(lattice_maxNonSchmid,maxNinstance))
nonSchmidCoeff = 0.0_pReal

!*** readout data from material.config file

rewind(myFile)
do while (trim(line) /= '#EOF#' .and. IO_lc(IO_getTag(line,'<','>')) /= 'phase')                                                  ! wind forward to <phase>
  line = IO_read(myFile)
enddo
 
do while (trim(line) /= '#EOF#')                                                                                                                                ! read thru sections of phase part
  line = IO_read(myFile)
  if (IO_isBlank(line)) cycle                                                                                                      ! skip empty lines
  if (IO_getTag(line,'<','>') /= '') exit                                                                                          ! stop at next part
  if (IO_getTag(line,'[',']') /= '') then                                                                                          ! next section
    section = section + 1_pInt                                                                                                     ! advance section counter
    cycle
  endif
  if (section > 0_pInt ) then                                                                       ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
    if (phase_plasticity(section) == CONSTITUTIVE_NONLOCAL_LABEL) then                              ! one of my sections
      i = phase_plasticityInstance(section)                                                                                          ! which instance of my plasticity is present phase
      positions = IO_stringPos(line,maxNchunks)
      tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                                                             ! extract key
      select case(tag)
        case('plasticity','elasticity','/nonlocal/')
          cycle
        case ('(output)')
          Noutput(i) = Noutput(i) + 1_pInt
          constitutive_nonlocal_output(Noutput(i),i) = IO_lc(IO_stringValue(line,positions,2_pInt))
        case ('lattice_structure')
          constitutive_nonlocal_structureName(i) = IO_lc(IO_stringValue(line,positions,2_pInt))
          configNchunks = lattice_configNchunks(constitutive_nonlocal_structureName(i))
          Nchunks_SlipFamilies = configNchunks(1)
          Nchunks_SlipSlip =     configNchunks(3)
        case ('c/a_ratio','covera_ratio')
          CoverA(i) = IO_floatValue(line,positions,2_pInt)
         case ('c11')
           Cslip66(1,1,i) = IO_floatValue(line,positions,2_pInt)
         case ('c12')
           Cslip66(1,2,i) = IO_floatValue(line,positions,2_pInt)
         case ('c13')
           Cslip66(1,3,i) = IO_floatValue(line,positions,2_pInt)
         case ('c22')
           Cslip66(2,2,i) = IO_floatValue(line,positions,2_pInt)
         case ('c23')
           Cslip66(2,3,i) = IO_floatValue(line,positions,2_pInt)
         case ('c33')
           Cslip66(3,3,i) = IO_floatValue(line,positions,2_pInt)
         case ('c44')
           Cslip66(4,4,i) = IO_floatValue(line,positions,2_pInt)
         case ('c55')
           Cslip66(5,5,i) = IO_floatValue(line,positions,2_pInt)
         case ('c66')
           Cslip66(6,6,i) = IO_floatValue(line,positions,2_pInt)
        case ('nslip')
          do f = 1_pInt, Nchunks_SlipFamilies
            Nslip(f,i) = IO_intValue(line,positions,1_pInt+f)
          enddo
        case ('rhosgledgepos0')
          do f = 1_pInt, Nchunks_SlipFamilies
            rhoSglEdgePos0(f,i) = IO_floatValue(line,positions,1_pInt+f)
          enddo
        case ('rhosgledgeneg0')
          do f = 1_pInt, Nchunks_SlipFamilies
            rhoSglEdgeNeg0(f,i) = IO_floatValue(line,positions,1_pInt+f)
          enddo
        case ('rhosglscrewpos0')
          do f = 1_pInt, Nchunks_SlipFamilies
            rhoSglScrewPos0(f,i) = IO_floatValue(line,positions,1_pInt+f)
          enddo
        case ('rhosglscrewneg0')
          do f = 1_pInt, Nchunks_SlipFamilies
            rhoSglScrewNeg0(f,i) = IO_floatValue(line,positions,1_pInt+f)
          enddo
        case ('rhodipedge0')
          do f = 1_pInt, Nchunks_SlipFamilies
            rhoDipEdge0(f,i) = IO_floatValue(line,positions,1_pInt+f)
          enddo
        case ('rhodipscrew0')
          do f = 1_pInt, Nchunks_SlipFamilies
            rhoDipScrew0(f,i) = IO_floatValue(line,positions,1_pInt+f)
          enddo
        case ('lambda0')
          do f = 1_pInt, Nchunks_SlipFamilies
            lambda0PerSlipFamily(f,i) = IO_floatValue(line,positions,1_pInt+f)
          enddo
        case ('burgers')
          do f = 1_pInt, Nchunks_SlipFamilies
            burgersPerSlipFamily(f,i) = IO_floatValue(line,positions,1_pInt+f)
          enddo
        case('cutoffradius','r')
          cutoffRadius(i) = IO_floatValue(line,positions,2_pInt)
        case('minimumdipoleheightedge','ddipminedge')
          do f = 1_pInt, Nchunks_SlipFamilies
            minDipoleHeightPerSlipFamily(f,1_pInt,i) = IO_floatValue(line,positions,1_pInt+f)
          enddo
        case('minimumdipoleheightscrew','ddipminscrew')
          do f = 1_pInt, Nchunks_SlipFamilies
            minDipoleHeightPerSlipFamily(f,2_pInt,i) = IO_floatValue(line,positions,1_pInt+f)
          enddo
        case('atomicvolume')
          atomicVolume(i) = IO_floatValue(line,positions,2_pInt)
        case('selfdiffusionprefactor','dsd0')
          Dsd0(i) = IO_floatValue(line,positions,2_pInt)
        case('selfdiffusionenergy','qsd')
          selfDiffusionEnergy(i) = IO_floatValue(line,positions,2_pInt)
        case('atol_rho','atol_density','absolutetolerancedensity','absolutetolerance_density')
          aTolRho(i) = IO_floatValue(line,positions,2_pInt)
        case('atol_shear','atol_plasticshear','atol_accumulatedshear','absolutetoleranceshear','absolutetolerance_shear')
          aTolShear(i) = IO_floatValue(line,positions,2_pInt)
        case('significantrho','significant_rho','significantdensity','significant_density')
          significantRho(i) = IO_floatValue(line,positions,2_pInt)
        case('significantn','significant_n','significantdislocations','significant_dislcations')
          significantN(i) = IO_floatValue(line,positions,2_pInt)
        case ('interaction_slipslip')
          do it = 1_pInt, Nchunks_SlipSlip
            interactionSlipSlip(it,i) = IO_floatValue(line,positions,1_pInt+it)
          enddo
        case('linetension','linetensioneffect','linetension_effect')
          linetensionEffect(i) = IO_floatValue(line,positions,2_pInt)
        case('edgejog','edgejogs','edgejogeffect','edgejog_effect')
          edgeJogFactor(i) = IO_floatValue(line,positions,2_pInt)
        case('peierlsstressedge','peierlsstress_edge')
          do f = 1_pInt, Nchunks_SlipFamilies
            peierlsStressPerSlipFamily(f,1_pInt,i) = IO_floatValue(line,positions,1_pInt+f)
          enddo
        case('peierlsstressscrew','peierlsstress_screw')
          do f = 1_pInt, Nchunks_SlipFamilies
            peierlsStressPerSlipFamily(f,2_pInt,i) = IO_floatValue(line,positions,1_pInt+f)
          enddo
        case('doublekinkwidth')
          doublekinkwidth(i) = IO_floatValue(line,positions,2_pInt)
        case('solidsolutionenergy')
          solidSolutionEnergy(i) = IO_floatValue(line,positions,2_pInt)
        case('solidsolutionsize')
          solidSolutionSize(i) = IO_floatValue(line,positions,2_pInt)
        case('solidsolutionconcentration')
          solidSolutionConcentration(i) = IO_floatValue(line,positions,2_pInt)
        case('p')
          pParam(i) = IO_floatValue(line,positions,2_pInt)
        case('q')
          qParam(i) = IO_floatValue(line,positions,2_pInt)
        case('viscosity','glideviscosity')
          viscosity(i) = IO_floatValue(line,positions,2_pInt)
        case('attackfrequency','fattack')
          fattack(i) = IO_floatValue(line,positions,2_pInt)
        case('rhosglscatter')
          rhoSglScatter(i) = IO_floatValue(line,positions,2_pInt)
        case('rhosglrandom')
          rhoSglRandom(i) = IO_floatValue(line,positions,2_pInt)
        case('rhosglrandombinning')
          rhoSglRandomBinning(i) = IO_floatValue(line,positions,2_pInt)
        case('surfacetransmissivity')
          surfaceTransmissivity(i) = IO_floatValue(line,positions,2_pInt)
        case('grainboundarytransmissivity')
          grainboundaryTransmissivity(i) = IO_floatValue(line,positions,2_pInt)
        case('cflfactor')
          CFLfactor(i) = IO_floatValue(line,positions,2_pInt)
        case('fedgemultiplication','edgemultiplicationfactor','edgemultiplication')
          fEdgeMultiplication(i) = IO_floatValue(line,positions,2_pInt)
        case('shortrangestresscorrection')
          shortRangeStressCorrection(i) = IO_floatValue(line,positions,2_pInt) > 0.0_pReal
        case ('nonschmid_coefficients')
          do f = 1_pInt, lattice_maxNonSchmid
            nonSchmidCoeff(f,i) = IO_floatValue(line,positions,1_pInt+f)
          enddo
        case('deadzonescaling','deadzone','deadscaling')
          deadZoneScaling(i) = IO_floatValue(line,positions,2_pInt) > 0.0_pReal
        case('probabilisticmultiplication','randomsources','randommultiplication','discretesources')
          probabilisticMultiplication(i) = IO_floatValue(line,positions,2_pInt) > 0.0_pReal
        case default
          call IO_error(210_pInt,ext_msg=trim(tag)//' ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
      end select
    endif
  endif
enddo


do i = 1_pInt,maxNinstance

  constitutive_nonlocal_structure(i) = &
    lattice_initializeStructure(constitutive_nonlocal_structureName(i), CoverA(i))                            ! our lattice structure is defined in the material.config file by the structureName (and the c/a ratio)
  myStructure = constitutive_nonlocal_structure(i)
  
  
  !*** sanity checks
  
  if (myStructure < 1_pInt) &
    call IO_error(205_pInt,e=i)
  if (sum(Nslip(:,i)) <= 0_pInt) &
    call IO_error(211_pInt,ext_msg='Nslip ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  do o = 1_pInt,maxval(phase_Noutput)
    if(len(constitutive_nonlocal_output(o,i)) > 64_pInt) &
      call IO_error(666_pInt)
  enddo
  do f = 1_pInt,lattice_maxNslipFamily
    if (Nslip(f,i) > 0_pInt) then
      if (rhoSglEdgePos0(f,i) < 0.0_pReal) &
        call IO_error(211_pInt,ext_msg='rhoSglEdgePos0 ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
      if (rhoSglEdgeNeg0(f,i) < 0.0_pReal) &
        call IO_error(211_pInt,ext_msg='rhoSglEdgeNeg0 ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
      if (rhoSglScrewPos0(f,i) < 0.0_pReal) &
        call IO_error(211_pInt,ext_msg='rhoSglScrewPos0 ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
      if (rhoSglScrewNeg0(f,i) < 0.0_pReal) &
        call IO_error(211_pInt,ext_msg='rhoSglScrewNeg0 ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
      if (rhoDipEdge0(f,i) < 0.0_pReal) &
        call IO_error(211_pInt,ext_msg='rhoDipEdge0 ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
      if (rhoDipScrew0(f,i) < 0.0_pReal) &
        call IO_error(211_pInt,ext_msg='rhoDipScrew0 ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
      if (burgersPerSlipFamily(f,i) <= 0.0_pReal) &
        call IO_error(211_pInt,ext_msg='Burgers ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
      if (lambda0PerSlipFamily(f,i) <= 0.0_pReal) &
        call IO_error(211_pInt,ext_msg='lambda0 ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
      if (minDipoleHeightPerSlipFamily(f,1,i) < 0.0_pReal) &
        call IO_error(211_pInt,ext_msg='minimumDipoleHeightEdge ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
      if (minDipoleHeightPerSlipFamily(f,2,i) < 0.0_pReal) &
        call IO_error(211_pInt,ext_msg='minimumDipoleHeightScrew ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
      if (peierlsStressPerSlipFamily(f,1,i) <= 0.0_pReal) &
        call IO_error(211_pInt,ext_msg='peierlsStressEdge ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
      if (peierlsStressPerSlipFamily(f,2,i) <= 0.0_pReal) &
        call IO_error(211_pInt,ext_msg='peierlsStressScrew ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
    endif
  enddo
  if (any(interactionSlipSlip(1:maxval(lattice_interactionSlipSlip(:,:,myStructure)),i) < 0.0_pReal)) &
    call IO_error(211_pInt,ext_msg='interaction_SlipSlip ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (linetensionEffect(i) < 0.0_pReal .or. linetensionEffect(i) > 1.0_pReal) &
    call IO_error(211_pInt,ext_msg='linetension ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (edgeJogFactor(i) < 0.0_pReal .or. edgeJogFactor(i) > 1.0_pReal) &
    call IO_error(211_pInt,ext_msg='edgejog ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (cutoffRadius(i) < 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='r ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (atomicVolume(i) <= 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='atomicVolume ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (Dsd0(i) < 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='selfDiffusionPrefactor ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (selfDiffusionEnergy(i) <= 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='selfDiffusionEnergy ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (aTolRho(i) <= 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='aTol_rho ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (aTolShear(i) <= 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='aTol_shear ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (significantRho(i) < 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='significantRho ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (significantN(i) < 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='significantN ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (doublekinkwidth(i) <= 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='doublekinkwidth ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (solidSolutionEnergy(i) <= 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='solidSolutionEnergy ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (solidSolutionSize(i) <= 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='solidSolutionSize ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (solidSolutionConcentration(i) <= 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='solidSolutionConcentration ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (pParam(i) <= 0.0_pReal .or. pParam(i) > 1.0_pReal) &
    call IO_error(211_pInt,ext_msg='p ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (qParam(i) < 1.0_pReal .or. qParam(i) > 2.0_pReal) &
    call IO_error(211_pInt,ext_msg='q ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (viscosity(i) <= 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='viscosity ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (fattack(i) <= 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='attackFrequency ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (rhoSglScatter(i) < 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='rhoSglScatter ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (rhoSglRandom(i) < 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='rhoSglRandom ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (rhoSglRandomBinning(i) <= 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='rhoSglRandomBinning ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (surfaceTransmissivity(i) < 0.0_pReal .or. surfaceTransmissivity(i) > 1.0_pReal) &
    call IO_error(211_pInt,ext_msg='surfaceTransmissivity ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (grainboundaryTransmissivity(i) > 1.0_pReal) &
    call IO_error(211_pInt,ext_msg='grainboundaryTransmissivity ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (CFLfactor(i) < 0.0_pReal) &
    call IO_error(211_pInt,ext_msg='CFLfactor ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  if (fEdgeMultiplication(i) < 0.0_pReal .or. fEdgeMultiplication(i) > 1.0_pReal) &
    call IO_error(211_pInt,ext_msg='edgemultiplicationfactor ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  
  
  !*** determine total number of active slip systems
  
  Nslip(1:lattice_maxNslipFamily,i) = min(lattice_NslipSystem(1:lattice_maxNslipFamily,myStructure), &
                                          Nslip(1:lattice_maxNslipFamily,i) )                         ! we can't use more slip systems per family than specified in lattice 
  totalNslip(i) = sum(Nslip(1:lattice_maxNslipFamily,i))

enddo


!*** allocation of variables whose size depends on the total number of active slip systems

maxTotalNslip = maxval(totalNslip)

allocate(iRhoU(maxTotalNslip,4,maxNinstance))
allocate(iRhoB(maxTotalNslip,4,maxNinstance))
allocate(iRhoD(maxTotalNslip,2,maxNinstance))
allocate(iV(maxTotalNslip,4,maxNinstance))
allocate(iD(maxTotalNslip,2,maxNinstance))
allocate(iGamma(maxTotalNslip,maxNinstance))
allocate(iRhoF(maxTotalNslip,maxNinstance))
allocate(iTauF(maxTotalNslip,maxNinstance))
allocate(iTauB(maxTotalNslip,maxNinstance))
iRhoU = 0_pInt
iRhoB = 0_pInt
iRhoD = 0_pInt
iV = 0_pInt
iD = 0_pInt
iGamma = 0_pInt
iRhoF = 0_pInt
iTauF = 0_pInt
iTauB = 0_pInt

allocate(burgers(maxTotalNslip, maxNinstance))
burgers = 0.0_pReal

allocate(lambda0(maxTotalNslip, maxNinstance))
lambda0 = 0.0_pReal

allocate(minDipoleHeight(maxTotalNslip,2,maxNinstance))
minDipoleHeight = -1.0_pReal

allocate(forestProjectionEdge(maxTotalNslip, maxTotalNslip, maxNinstance))
forestProjectionEdge = 0.0_pReal

allocate(forestProjectionScrew(maxTotalNslip, maxTotalNslip, maxNinstance))
forestProjectionScrew = 0.0_pReal

allocate(interactionMatrixSlipSlip(maxTotalNslip, maxTotalNslip, maxNinstance))
interactionMatrixSlipSlip = 0.0_pReal

allocate(lattice2slip(1:3, 1:3, maxTotalNslip, maxNinstance))
lattice2slip = 0.0_pReal

allocate(sourceProbability(maxTotalNslip, homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems))
sourceProbability = 2.0_pReal

allocate(rhoDotFluxOutput(maxTotalNslip, 8, homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems))
allocate(rhoDotMultiplicationOutput(maxTotalNslip, 2, homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems))
allocate(rhoDotSingle2DipoleGlideOutput(maxTotalNslip, 2, homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems))
allocate(rhoDotAthermalAnnihilationOutput(maxTotalNslip, 2, homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems))
allocate(rhoDotThermalAnnihilationOutput(maxTotalNslip, 2, homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems))
allocate(rhoDotEdgeJogsOutput(maxTotalNslip, homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems))
rhoDotFluxOutput = 0.0_pReal
rhoDotMultiplicationOutput = 0.0_pReal
rhoDotSingle2DipoleGlideOutput = 0.0_pReal
rhoDotAthermalAnnihilationOutput = 0.0_pReal
rhoDotThermalAnnihilationOutput = 0.0_pReal
rhoDotEdgeJogsOutput = 0.0_pReal

allocate(compatibility(2,maxTotalNslip, maxTotalNslip, mesh_maxNipNeighbors, mesh_maxNips, mesh_NcpElems))
compatibility = 0.0_pReal

allocate(peierlsStress(maxTotalNslip,2,maxNinstance))
peierlsStress = 0.0_pReal

allocate(colinearSystem(maxTotalNslip,maxNinstance))
colinearSystem = 0_pInt

do i = 1,maxNinstance
  
  myStructure = constitutive_nonlocal_structure(i)                                                                                  ! lattice structure of this instance
    

  !*** Inverse lookup of my slip system family and the slip system in lattice
  
  l = 0_pInt
  do f = 1_pInt,lattice_maxNslipFamily
    do s = 1_pInt,Nslip(f,i)
      l = l + 1_pInt
      slipFamily(l,i) = f
      slipSystemLattice(l,i) = sum(lattice_NslipSystem(1:f-1_pInt, myStructure)) + s
  enddo; enddo
  
  
  !*** determine size of state array
  
  ns = totalNslip(i)
  constitutive_nonlocal_sizeDotState(i) = int(size(BASICSTATES),pInt) * ns
  constitutive_nonlocal_sizeDependentState(i) = int(size(DEPENDENTSTATES),pInt) * ns
  constitutive_nonlocal_sizeState(i) = constitutive_nonlocal_sizeDotState(i) &
                                     + constitutive_nonlocal_sizeDependentState(i) &
                                     + int(size(OTHERSTATES),pInt) * ns

  !*** determine indices to state array

  l = 0_pInt
  do t = 1_pInt,4_pInt
    do s = 1_pInt,ns
      l = l + 1_pInt
      iRhoU(s,t,i) = l
    enddo
  enddo
  do t = 1_pInt,4_pInt
    do s = 1_pInt,ns
      l = l + 1_pInt
      iRhoB(s,t,i) = l
    enddo
  enddo
  do c = 1_pInt,2_pInt
    do s = 1_pInt,ns
      l = l + 1_pInt
      iRhoD(s,c,i) = l
    enddo
  enddo
  do s = 1_pInt,ns
    l = l + 1_pInt
    iGamma(s,i) = l
  enddo
  do s = 1_pInt,ns
    l = l + 1_pInt
    iRhoF(s,i) = l
  enddo
  do s = 1_pInt,ns
    l = l + 1_pInt
    iTauF(s,i) = l
  enddo
  do s = 1_pInt,ns
    l = l + 1_pInt
    iTauB(s,i) = l
  enddo
  do t = 1_pInt,4_pInt
    do s = 1_pInt,ns
      l = l + 1_pInt
      iV(s,t,i) = l
    enddo
  enddo
  do c = 1_pInt,2_pInt
    do s = 1_pInt,ns
      l = l + 1_pInt
      iD(s,c,i) = l
    enddo
  enddo
  if (iD(ns,2,i) /= constitutive_nonlocal_sizeState(i)) &  ! check if last index is equal to size of state
    call IO_error(0_pInt, ext_msg = 'state indices not properly set ('//CONSTITUTIVE_NONLOCAL_LABEL//')')
  

  !*** determine size of postResults array
  
  do o = 1_pInt,Noutput(i)
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
            'rho_dot_sgl2dip_edge', &
            'rho_dot_sgl2dip_screw', &
            'rho_dot_ann_ath', &
            'rho_dot_ann_the', &
            'rho_dot_ann_the_edge', &
            'rho_dot_ann_the_screw', &
            'rho_dot_edgejogs', &
            'rho_dot_flux', &
            'rho_dot_flux_edge', &
            'rho_dot_flux_screw', &
            'velocity_edge_pos', &
            'velocity_edge_neg', &
            'velocity_screw_pos', &
            'velocity_screw_neg', &
            'slipdirection.x', &
            'slipdirection.y', &
            'slipdirection.z', &
            'slipnormal.x', &
            'slipnormal.y', &
            'slipnormal.z', &
            'fluxdensity_edge_pos.x', &
            'fluxdensity_edge_pos.y', &
            'fluxdensity_edge_pos.z', &
            'fluxdensity_edge_neg.x', &
            'fluxdensity_edge_neg.y', &
            'fluxdensity_edge_neg.z', &
            'fluxdensity_screw_pos.x', &
            'fluxdensity_screw_pos.y', &
            'fluxdensity_screw_pos.z', &
            'fluxdensity_screw_neg.x', &
            'fluxdensity_screw_neg.y', &
            'fluxdensity_screw_neg.z', &
            'maximumdipoleheight_edge', &
            'maximumdipoleheight_screw', &
            'accumulatedshear', &
            'boundarylayer' )
        mySize = totalNslip(i)
      case('dislocationstress')
        mySize = 6_pInt
      case default
        call IO_error(212_pInt,ext_msg=constitutive_nonlocal_output(o,i)//&
                                       '('//CONSTITUTIVE_NONLOCAL_LABEL//')')
    end select

    if (mySize > 0_pInt) then                                                                       ! any meaningful output found                               
      constitutive_nonlocal_sizePostResult(o,i) = mySize
      constitutive_nonlocal_sizePostResults(i)  = constitutive_nonlocal_sizePostResults(i) + mySize
    endif
  enddo
  
  
  !*** elasticity matrix and shear modulus according to material.config
  
  Cslip66(:,:,i) = lattice_symmetrizeC66(constitutive_nonlocal_structureName(i), Cslip66(:,:,i)) 
  mu(i) = 0.2_pReal * ( Cslip66(1,1,i) - Cslip66(1,2,i) + 3.0_pReal*Cslip66(4,4,i))                 ! (C11iso-C12iso)/2 with C11iso=(3*C11+2*C12+4*C44)/5 and C12iso=(C11+4*C12-2*C44)/5
  nu(i) = (Cslip66(1,1,i) + 4.0_pReal*Cslip66(1,2,i) - 2.0_pReal*Cslip66(4,4,i)) &
        / (4.0_pReal*Cslip66(1,1,i) + 6.0_pReal*Cslip66(1,2,i) + 2.0_pReal*Cslip66(4,4,i))          ! C12iso/(C11iso+C12iso) with C11iso=(3*C11+2*C12+4*C44)/5 and C12iso=(C11+4*C12-2*C44)/5
  Cslip66(1:6,1:6,i) = math_Mandel3333to66(math_Voigt66to3333(Cslip66(1:6,1:6,i)))
  Cslip3333(1:3,1:3,1:3,1:3,i) = math_Voigt66to3333(Cslip66(1:6,1:6,i))
  
  do s1 = 1_pInt,ns 
    f = slipFamily(s1,i)
    
    !*** burgers vector, mean free path prefactor and minimum dipole distance for each slip system
  
    burgers(s1,i) = burgersPerSlipFamily(f,i)
    lambda0(s1,i) = lambda0PerSlipFamily(f,i)
    minDipoleHeight(s1,1:2,i) = minDipoleHeightPerSlipFamily(f,1:2,i)
    peierlsStress(s1,1:2,i) = peierlsStressPerSlipFamily(f,1:2,i)

    do s2 = 1_pInt,ns
      
      !*** calculation of forest projections for edge and screw dislocations. s2 acts as forest for s1

      forestProjectionEdge(s1,s2,i) &
          = abs(math_mul3x3(lattice_sn(1:3,slipSystemLattice(s1,i),myStructure), &
                            lattice_st(1:3,slipSystemLattice(s2,i),myStructure)))                   ! forest projection of edge dislocations is the projection of (t = b x n) onto the slip normal of the respective slip plane
      
      forestProjectionScrew(s1,s2,i) &
          = abs(math_mul3x3(lattice_sn(1:3,slipSystemLattice(s1,i),myStructure), &
                            lattice_sd(1:3,slipSystemLattice(s2,i),myStructure)))                   ! forest projection of screw dislocations is the projection of b onto the slip normal of the respective splip plane
  
      !*** calculation of interaction matrices

      interactionMatrixSlipSlip(s1,s2,i) &
          = interactionSlipSlip(lattice_interactionSlipSlip(slipSystemLattice(s1,i), &
                                                            slipSystemLattice(s2,i), &
                                                            myStructure), i)
      
      !*** colinear slip system (only makes sense for fcc like it is defined here)
      
      if (lattice_interactionSlipSlip(slipSystemLattice(s1,i), &
                                      slipSystemLattice(s2,i), &
                                      myStructure) == 3_pInt) then
        colinearSystem(s1,i) = s2
      endif
  
    enddo

    !*** rotation matrix from lattice configuration to slip system

    lattice2slip(1:3,1:3,s1,i) &
        = math_transpose33( reshape([ lattice_sd(1:3, slipSystemLattice(s1,i), myStructure), &
                                     -lattice_st(1:3, slipSystemLattice(s1,i), myStructure), &
                                      lattice_sn(1:3, slipSystemLattice(s1,i), myStructure)], [3,3]))
  enddo
  
enddo

endsubroutine



!*********************************************************************
!* initial microstructural state (just the "basic" states)           *
!*********************************************************************
subroutine constitutive_nonlocal_stateInit(state)

use IO,       only: IO_error
use lattice,  only: lattice_maxNslipFamily
use math,     only: math_sampleGaussVar
use mesh,     only: mesh_ipVolume, &
                    mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_element, &
                    FE_Nips, &
                    FE_geomtype
use material, only: material_phase, &
                    phase_plasticityInstance, &
                    phase_plasticity, &
                    homogenization_Ngrains

implicit none

!*** input/output variables
type(p_vec), dimension(1,mesh_maxNips,mesh_NcpElems), intent(inout) :: &
                              state                           ! microstructural state

!*** local variables
integer(pInt)                 el, &
                              ip, &
                              e, &
                              i, &
                              idx, &
                              ns, &                           ! short notation for total number of active slip systems 
                              f, &                            ! index of lattice family
                              from, &
                              upto, &
                              s, &                            ! index of slip system
                              t, &
                              j, &
                              myInstance, &
                              maxNinstance
real(pReal), dimension(2) ::  noise
real(pReal), dimension(4) ::  rnd   
real(pReal)                   meanDensity, &
                              totalVolume, &
                              densityBinning, &
                              minimumIpVolume


maxNinstance = int(count(phase_plasticity == CONSTITUTIVE_NONLOCAL_LABEL),pInt)


! ititalize all states to zero

do e = 1_pInt,mesh_NcpElems
  do i = 1_pInt,FE_Nips(FE_geomtype(mesh_element(2,e)))
    if (CONSTITUTIVE_NONLOCAL_LABEL == phase_plasticity(material_phase(1,i,e))) &
      state(1,i,e)%p = 0.0_pReal
  enddo
enddo


do myInstance = 1_pInt,maxNinstance
  ns = totalNslip(myInstance)

  ! randomly distribute dislocation segments on random slip system and of random type in the volume 
  if (rhoSglRandom(myInstance) > 0.0_pReal) then

    ! get the total volume of the instance

    minimumIpVolume = 1e99_pReal
    totalVolume = 0.0_pReal
    do e = 1_pInt,mesh_NcpElems
      do i = 1_pInt,FE_Nips(FE_geomtype(mesh_element(2,e)))
        if (CONSTITUTIVE_NONLOCAL_LABEL == phase_plasticity(material_phase(1,i,e)) &
            .and. myInstance == phase_plasticityInstance(material_phase(1,i,e))) then
          totalVolume = totalVolume + mesh_ipVolume(i,e)
          minimumIpVolume = min(minimumIpVolume, mesh_ipVolume(i,e))
        endif
      enddo
    enddo
    densityBinning = rhoSglRandomBinning(myInstance) / minimumIpVolume ** (2.0_pReal / 3.0_pReal)

    ! subsequently fill random ips with dislocation segments until we reach the desired overall density

    meanDensity = 0.0_pReal
    do while(meanDensity < rhoSglRandom(myInstance))
      call random_number(rnd)
      el = nint(rnd(1)*real(mesh_NcpElems,pReal)+0.5_pReal,pInt)
      ip = nint(rnd(2)*real(FE_Nips(FE_geomtype(mesh_element(2,el))),pReal)+0.5_pReal,pInt)
      if (CONSTITUTIVE_NONLOCAL_LABEL == phase_plasticity(material_phase(1,ip,el)) &
          .and. myInstance == phase_plasticityInstance(material_phase(1,ip,el))) then
        s = nint(rnd(3)*real(ns,pReal)+0.5_pReal,pInt)
        t = nint(rnd(4)*4.0_pReal+0.5_pReal,pInt)
        meanDensity = meanDensity + densityBinning * mesh_ipVolume(ip,el) / totalVolume
        state(1,ip,el)%p(iRhoU(s,t,myInstance)) = state(1,ip,el)%p(iRhoU(s,t,myInstance)) + densityBinning
      endif
    enddo

  ! homogeneous distribution of density with some noise
  else
    do e = 1_pInt,mesh_NcpElems
      do i = 1_pInt,FE_Nips(FE_geomtype(mesh_element(2,e)))
        if (CONSTITUTIVE_NONLOCAL_LABEL == phase_plasticity(material_phase(1,i,e)) &
            .and. myInstance == phase_plasticityInstance(material_phase(1,i,e))) then
          do f = 1_pInt,lattice_maxNslipFamily
            from = 1_pInt + sum(Nslip(1:f-1_pInt,myInstance))
            upto = sum(Nslip(1:f,myInstance))
            do s = from,upto
              do j = 1_pInt,2_pInt
                noise(j) = math_sampleGaussVar(0.0_pReal, rhoSglScatter(myInstance))
              enddo
              state(1,i,e)%p(iRhoU(s,1,myInstance)) = rhoSglEdgePos0(f, myInstance) + noise(1)
              state(1,i,e)%p(iRhoU(s,2,myInstance)) = rhoSglEdgeNeg0(f, myInstance) + noise(1)
              state(1,i,e)%p(iRhoU(s,3,myInstance)) = rhoSglScrewPos0(f, myInstance) + noise(2)
              state(1,i,e)%p(iRhoU(s,4,myInstance)) = rhoSglScrewNeg0(f, myInstance) + noise(2)
            enddo 
            state(1,i,e)%p(iRhoD(from:upto,1,myInstance)) = rhoDipEdge0(f, myInstance)
            state(1,i,e)%p(iRhoD(from:upto,2,myInstance)) = rhoDipScrew0(f, myInstance)
          enddo
        endif
      enddo
    enddo
  endif
enddo

endsubroutine



!*********************************************************************
!* absolute state tolerance                                          *
!*********************************************************************
pure function constitutive_nonlocal_aTolState(myInstance)

implicit none

!*** input variables
integer(pInt), intent(in) :: myInstance                       ! number specifying the current instance of the plasticity

!*** output variables
real(pReal), dimension(constitutive_nonlocal_sizeState(myInstance)) :: &
                              constitutive_nonlocal_aTolState ! absolute state tolerance for the current instance of this plasticity

!*** local variables
integer(pInt)             :: ns, t, c

ns = totalNslip(myInstance)
constitutive_nonlocal_aTolState = 0.0_pReal
forall (t = 1_pInt:4_pInt)
  constitutive_nonlocal_aTolState(iRhoU(1:ns,t,myInstance)) = aTolRho(myInstance)
  constitutive_nonlocal_aTolState(iRhoB(1:ns,t,myInstance)) = aTolRho(myInstance)
endforall
forall (c = 1_pInt:2_pInt) &
  constitutive_nonlocal_aTolState(iRhoD(1:ns,c,myInstance)) = aTolRho(myInstance)
constitutive_nonlocal_aTolState(iGamma(1:ns,myInstance)) = aTolShear(myInstance)

endfunction



!*********************************************************************
!* calculates homogenized elacticity matrix                          *
!*********************************************************************
pure function constitutive_nonlocal_homogenizedC(state,g,ip,el)

use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_plasticityInstance
implicit none

!*** input variables
integer(pInt), intent(in) ::    g, &                                ! current grain ID
                                ip, &                               ! current integration point
                                el                                  ! current element
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state ! microstructural state

!*** output variables
real(pReal), dimension(6,6) ::  constitutive_nonlocal_homogenizedC  ! homogenized elasticity matrix

!*** local variables
integer(pInt)                   myInstance                          ! current instance of this plasticity

myInstance = phase_plasticityInstance(material_phase(g,ip,el))

constitutive_nonlocal_homogenizedC = Cslip66(1:6,1:6,myInstance)
 
endfunction



!*********************************************************************
!* calculates quantities characterizing the microstructure           *
!*********************************************************************
subroutine constitutive_nonlocal_microstructure(state, Temperature, Fe, Fp, gr, ip, el)

use IO, only: &
  IO_error
use math, only: &
  pi, &
  math_mul33x3, &
  math_mul3x3, &
  math_norm3, &
  math_invert33, &
  math_transpose33
use debug, only: &
  debug_level, &
  debug_constitutive, &
  debug_levelBasic, &
  debug_levelExtensive, &
  debug_levelSelective, &
  debug_g, &
  debug_i, &
  debug_e
use mesh, only: &
  mesh_NcpElems, &
  mesh_maxNips, &
  mesh_element, &
  mesh_ipNeighborhood, &
  mesh_ipCoordinates, &
  mesh_ipVolume, &
  mesh_ipAreaNormal, &
  mesh_ipArea, &
  FE_NipNeighbors, &
  mesh_maxNipNeighbors, &
  FE_geomtype, &
  FE_celltype
use material, only: &
  homogenization_maxNgrains, &
  material_phase, &
  phase_localPlasticity, &
  phase_plasticityInstance
use lattice, only: &
  lattice_sd, &
  lattice_st, &
  lattice_interactionSlipSlip

implicit none

!*** input variables
integer(pInt), intent(in) ::    gr, &                          ! current grain ID
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
                                instance, &                   ! my instance of this plasticity
                                neighboring_instance, &       ! instance of this plasticity of neighboring material point
                                latticeStruct, &              ! my lattice structure
                                neighboring_latticeStruct, &  ! lattice structure of neighboring material point
                                phase, &
                                neighboring_phase, &
                                ns, &                         ! total number of active slip systems at my material point
                                neighboring_ns, &             ! total number of active slip systems at neighboring material point
                                c, &                          ! index of dilsocation character (edge, screw)
                                s, &                          ! slip system index
                                s2, &                         ! slip system index
                                t, &                          ! index of dilsocation type (e+, e-, s+, s-, used e+, used e-, used s+, used s-)
                                dir, &
                                n, &
                                nRealNeighbors, &             ! number of really existing neighbors
                                interactionCoefficient
integer(pInt), dimension(2) ::  neighbor
real(pReal)                     detFe, &
                                detFp, &
                                FVsize, &
                                temp, &
                                correction, &
                                myRhoForest
real(pReal), dimension(2) ::    rhoExcessGradient, &
                                rhoExcessGradient_over_rho, &
                                rhoTotal
real(pReal), dimension(3) ::    rhoExcessDifferences, &
                                normal_latticeConf
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(gr,ip,el)))) :: &
                                rhoForest, &                  ! forest dislocation density
                                tauBack, &                    ! back stress from pileup on same slip system
                                tauThreshold                  ! threshold shear stress
real(pReal), dimension(3,3) ::  invFe, &                      ! inverse of elastic deformation gradient
                                invFp, &                      ! inverse of plastic deformation gradient
                                connections, &
                                invConnections
real(pReal), dimension(3,mesh_maxNipNeighbors) :: &
                                connection_latticeConf
real(pReal), dimension(2,totalNslip(phase_plasticityInstance(material_phase(gr,ip,el)))) :: &
                                rhoExcess
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(gr,ip,el))),2) :: &
                                rhoDip                        ! dipole dislocation density (edge, screw)
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(gr,ip,el))),8) :: &
                                rhoSgl                        ! single dislocation density (edge+, edge-, screw+, screw-, used edge+, used edge-, used screw+, used screw-)
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(gr,ip,el))), &
                       totalNslip(phase_plasticityInstance(material_phase(gr,ip,el)))) :: &
                                myInteractionMatrix           ! corrected slip interaction matrix
real(pReal), dimension(2,maxval(totalNslip),mesh_maxNipNeighbors) :: &
                                neighboring_rhoExcess, &      ! excess density at neighboring material point
                                neighboring_rhoTotal          ! total density at neighboring material point
real(pReal), dimension(3,totalNslip(phase_plasticityInstance(material_phase(gr,ip,el))),2) :: &
                                m                             ! direction of dislocation motion
logical                         inversionError


phase = material_phase(gr,ip,el)
instance = phase_plasticityInstance(phase)
latticeStruct = constitutive_nonlocal_structure(instance)
ns = totalNslip(instance)


!*** get basic states

forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)
  rhoSgl(s,t) = max(state(gr,ip,el)%p(iRhoU(s,t,instance)), 0.0_pReal)                              ! ensure positive single mobile densities
  rhoSgl(s,t+4_pInt) = state(gr,ip,el)%p(iRhoB(s,t,instance))
endforall
forall (s = 1_pInt:ns, c = 1_pInt:2_pInt) &
  rhoDip(s,c) = max(state(gr,ip,el)%p(iRhoD(s,c,instance)), 0.0_pReal)                              ! ensure positive dipole densities
where (abs(rhoSgl) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(instance) &
  .or. abs(rhoSgl) < significantRho(instance)) &
  rhoSgl = 0.0_pReal
where (abs(rhoDip) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(instance) &
  .or. abs(rhoDip) < significantRho(instance)) &
  rhoDip = 0.0_pReal


!*** calculate the forest dislocation density
!*** (= projection of screw and edge dislocations)

forall (s = 1_pInt:ns) &
  rhoForest(s) = dot_product((sum(abs(rhoSgl(1:ns,[1,2,5,6])),2) + rhoDip(1:ns,1)), &
                              forestProjectionEdge(s,1:ns,instance)) & 
               + dot_product((sum(abs(rhoSgl(1:ns,[3,4,7,8])),2) + rhoDip(1:ns,2)), &
                              forestProjectionScrew(s,1:ns,instance))



!*** calculate the threshold shear stress for dislocation slip 

myInteractionMatrix = 0.0_pReal
myInteractionMatrix(1:ns,1:ns) = interactionMatrixSlipSlip(1:ns,1:ns,instance)
if (latticeStruct == 1_pInt) then                                                                   ! in case of fcc: coefficients are corrected for the line tension effect (see Kubin,Devincre,Hoc; 2008; Modeling dislocation storage rates and mean free paths in face-centered cubic crystals)
  do s = 1_pInt,ns 
    myRhoForest = max(rhoForest(s),significantRho(instance))
    correction = (  1.0_pReal - linetensionEffect(instance) &
                  + linetensionEffect(instance) &
                  * log(0.35_pReal * burgers(s,instance) * sqrt(myRhoForest)) &
                  / log(0.35_pReal * burgers(s,instance) * 1e6_pReal)) ** 2.0_pReal
    do s2 = 1_pInt,ns
      interactionCoefficient = &
        lattice_interactionSlipSlip(slipSystemLattice(s,instance), &
                                    slipSystemLattice(s2,instance), &
                                    latticeStruct)
      select case(interactionCoefficient)
        case(4_pInt,5_pInt,6_pInt)                                                                  ! only correct junction forming interactions (4,5,6)
          myInteractionMatrix(s,s2) = correction * myInteractionMatrix(s,s2) 
      endselect
    enddo
  enddo
endif
forall (s = 1_pInt:ns) &
  tauThreshold(s) = mu(instance) * burgers(s,instance) &
                  * sqrt(dot_product((sum(abs(rhoSgl),2) + sum(abs(rhoDip),2)), myInteractionMatrix(s,1:ns)))



!*** calculate the dislocation stress of the neighboring excess dislocation densities
!*** zero for material points of local plasticity

tauBack = 0.0_pReal

if (.not. phase_localPlasticity(phase) .and. shortRangeStressCorrection(instance)) then
  call math_invert33(Fe, invFe, detFe, inversionError)
  call math_invert33(Fp, invFp, detFp, inversionError)
  rhoExcess(1,1:ns) = rhoSgl(1:ns,1) - rhoSgl(1:ns,2)
  rhoExcess(2,1:ns) = rhoSgl(1:ns,3) - rhoSgl(1:ns,4)
  FVsize = mesh_ipVolume(ip,el) ** (1.0_pReal/3.0_pReal)
  
  !* loop through my neighborhood and get the connection vectors (in lattice frame) and the excess densities
  
  nRealNeighbors = 0_pInt
  neighboring_rhoTotal = 0.0_pReal
  do n = 1_pInt,FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,el))))
    neighboring_el = mesh_ipNeighborhood(1,n,ip,el)
    neighboring_ip = mesh_ipNeighborhood(2,n,ip,el)
    if (neighboring_el > 0 .and. neighboring_ip > 0) then
      neighboring_phase = material_phase(gr,neighboring_ip,neighboring_el)
      neighboring_instance = phase_plasticityInstance(neighboring_phase)
      neighboring_latticeStruct = constitutive_nonlocal_structure(neighboring_instance)
      neighboring_ns = totalNslip(neighboring_instance)
      if (.not. phase_localPlasticity(neighboring_phase) &
          .and. neighboring_latticeStruct == latticeStruct & 
          .and. neighboring_instance == instance) then
        if (neighboring_ns == ns) then
          nRealNeighbors = nRealNeighbors + 1_pInt
          forall (s = 1_pInt:ns, c = 1_pInt:2_pInt)
            neighboring_rhoExcess(c,s,n) = &
                max(state(gr,neighboring_ip,neighboring_el)%p(iRhoU(s,2*c-1,neighboring_instance)), 0.0_pReal) &! positive mobiles
              - max(state(gr,neighboring_ip,neighboring_el)%p(iRhoU(s,2*c,neighboring_instance)), 0.0_pReal)    ! negative mobiles
            neighboring_rhoTotal(c,s,n) = &
                max(state(gr,neighboring_ip,neighboring_el)%p(iRhoU(s,2*c-1,neighboring_instance)), 0.0_pReal) &! positive mobiles
              + max(state(gr,neighboring_ip,neighboring_el)%p(iRhoU(s,2*c,neighboring_instance)), 0.0_pReal) &  ! negative mobiles
              + abs(state(gr,neighboring_ip,neighboring_el)%p(iRhoB(s,2*c-1,neighboring_instance))) &           ! positive deads
              + abs(state(gr,neighboring_ip,neighboring_el)%p(iRhoB(s,2*c,neighboring_instance))) &             ! negative deads
              + max(state(gr,neighboring_ip,neighboring_el)%p(iRhoD(s,c,neighboring_instance)), 0.0_pReal)      ! dipoles
          endforall
          connection_latticeConf(1:3,n) = &
            math_mul33x3(invFe, mesh_ipCoordinates(1:3,neighboring_ip,neighboring_el) &
                                - mesh_ipCoordinates(1:3,ip,el))
          normal_latticeConf = math_mul33x3(math_transpose33(invFp), mesh_ipAreaNormal(1:3,n,ip,el))
          if (math_mul3x3(normal_latticeConf,connection_latticeConf(1:3,n)) < 0.0_pReal) then       ! neighbor connection points in opposite direction to face normal: must be periodic image
            connection_latticeConf(1:3,n) = normal_latticeConf * mesh_ipVolume(ip,el) &
                                                               / mesh_ipArea(n,ip,el)               ! instead take the surface normal scaled with the diameter of the cell
          endif
        else
          ! different number of active slip systems
          call IO_error(-1_pInt,ext_msg='different number of active slip systems in neighboring IPs of same crystal structure')
        endif
      else
        ! local neighbor or different lattice structure or different constitution instance -> use central values instead
        connection_latticeConf(1:3,n) = 0.0_pReal
        neighboring_rhoExcess(1:2,1:ns,n) = rhoExcess
      endif
    else
      ! free surface -> use central values instead
      connection_latticeConf(1:3,n) = 0.0_pReal
      neighboring_rhoExcess(1:2,1:ns,n) = rhoExcess
    endif
  enddo
  

  !* loop through the slip systems and calculate the dislocation gradient by
  !* 1. interpolation of the excess density in the neighorhood
  !* 2. interpolation of the dead dislocation density in the central volume
  
  m(1:3,1:ns,1) =  lattice_sd(1:3,slipSystemLattice(1:ns,instance),latticeStruct)
  m(1:3,1:ns,2) = -lattice_st(1:3,slipSystemLattice(1:ns,instance),latticeStruct)

  do s = 1_pInt,ns
    
    !* gradient from interpolation of neighboring excess density

    do c = 1_pInt,2_pInt
      do dir = 1_pInt,3_pInt
        neighbor(1) = 2_pInt * dir - 1_pInt
        neighbor(2) = 2_pInt * dir
        connections(dir,1:3) = connection_latticeConf(1:3,neighbor(1)) &
                             - connection_latticeConf(1:3,neighbor(2))
        rhoExcessDifferences(dir) = neighboring_rhoExcess(c,s,neighbor(1)) &
                                  - neighboring_rhoExcess(c,s,neighbor(2))
      enddo
      call math_invert33(connections,invConnections,temp,inversionError)
      if (inversionError) then
        call IO_error(-1_pInt,ext_msg='back stress calculation: inversion error')
      endif
      rhoExcessGradient(c) = math_mul3x3(m(1:3,s,c), &
                                         math_mul33x3(invConnections,rhoExcessDifferences))
    enddo
      
    !* plus gradient from deads
    
    do t = 1_pInt,4_pInt
      c = (t - 1_pInt) / 2_pInt + 1_pInt
      rhoExcessGradient(c) = rhoExcessGradient(c) + rhoSgl(s,t+4_pInt) / FVsize
    enddo

    !* normalized with the total density
    
    rhoExcessGradient_over_rho = 0.0_pReal
    forall (c = 1_pInt:2_pInt) &
      rhoTotal(c) = (sum(abs(rhoSgl(s,[2*c-1,2*c,2*c+3,2*c+4]))) + rhoDip(s,c) + sum(neighboring_rhoTotal(c,s,:))) &
                  / real(1_pInt + nRealNeighbors,pReal)
    forall (c = 1_pInt:2_pInt, rhoTotal(c) > 0.0_pReal) &
      rhoExcessGradient_over_rho(c) = rhoExcessGradient(c) / rhoTotal(c)
    
    !* gives the local stress correction when multiplied with a factor

    tauBack(s) = - mu(instance) * burgers(s,instance) / (2.0_pReal * pi) &
               * (rhoExcessGradient_over_rho(1) / (1.0_pReal - nu(instance)) + rhoExcessGradient_over_rho(2))

  enddo
endif


!*** set dependent states

state(gr,ip,el)%p(iRhoF(1:ns,instance)) = rhoForest
state(gr,ip,el)%p(iTauF(1:ns,instance)) = tauThreshold
state(gr,ip,el)%p(iTauB(1:ns,instance)) = tauBack


#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip .and. debug_g == gr)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
    write(6,*)
    write(6,'(a,i8,1x,i2,1x,i1)') '<< CONST >> nonlocal_microstructure at el ip g',el,ip,gr
    write(6,*)
    write(6,'(a,/,12x,12(e10.3,1x))') '<< CONST >> rhoForest', rhoForest
    write(6,'(a,/,12x,12(f10.5,1x))') '<< CONST >> tauThreshold / MPa', tauThreshold/1e6
    write(6,'(a,/,12x,12(f10.5,1x))') '<< CONST >> tauBack / MPa', tauBack/1e6
    write(6,*)
  endif
#endif

endsubroutine



!*********************************************************************
!* calculates kinetics                                               *
!*********************************************************************
subroutine constitutive_nonlocal_kinetics(v, tau, c, Temperature, state, g, ip, el, dv_dtau)

use debug,    only: debug_level, &
                    debug_constitutive, &
                    debug_levelBasic, &
                    debug_levelExtensive, &
                    debug_levelSelective, &
                    debug_g, &
                    debug_i, &
                    debug_e
use material, only: material_phase, &
                    phase_plasticityInstance

implicit none

!*** input variables
integer(pInt), intent(in) ::                g, &                        ! current grain number
                                            ip, &                       ! current integration point
                                            el, &                       ! current element number
                                            c                           ! dislocation character (1:edge, 2:screw)
real(pReal), intent(in) ::                  Temperature                 ! temperature
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el)))), &
             intent(in) ::                  tau                         ! resolved external shear stress (for bcc this already contains non Schmid effects)
type(p_vec), intent(in) ::                  state                       ! microstructural state

!*** input/output variables

!*** output variables
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el)))), &
                            intent(out) ::  v                           ! velocity
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el)))), &
                   intent(out), optional :: dv_dtau                     ! velocity derivative with respect to resolved shear stress

!*** local variables
integer(pInt)    ::                         instance, &                 ! current instance of this plasticity
                                            ns, &                       ! short notation for the total number of active slip systems
                                            s                           ! index of my current slip system
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el)))) :: &
                                            tauThreshold, &             ! threshold shear stress
                                            tauEff                      ! effective shear stress
real(pReal)                                 tauRel_P, & 
                                            tauRel_S, &
                                            tPeierls, &                 ! waiting time in front of a peierls barriers
                                            tSolidSolution, &           ! waiting time in front of a solid solution obstacle
                                            vViscous, &                 ! viscous glide velocity
                                            dtPeierls_dtau, &           ! derivative with respect to resolved shear stress
                                            dtSolidSolution_dtau, &     ! derivative with respect to resolved shear stress
                                            meanfreepath_S, &           ! mean free travel distance for dislocations between two solid solution obstacles
                                            meanfreepath_P, &           ! mean free travel distance for dislocations between two Peierls barriers
                                            jumpWidth_P, &              ! depth of activated area
                                            jumpWidth_S, &              ! depth of activated area
                                            activationLength_P, &       ! length of activated dislocation line
                                            activationLength_S, &       ! length of activated dislocation line
                                            activationVolume_P, &       ! volume that needs to be activated to overcome barrier
                                            activationVolume_S, &       ! volume that needs to be activated to overcome barrier
                                            activationEnergy_P, &       ! energy that is needed to overcome barrier
                                            activationEnergy_S, &       ! energy that is needed to overcome barrier
                                            criticalStress_P, &         ! maximum obstacle strength
                                            criticalStress_S, &         ! maximum obstacle strength
                                            mobility                    ! dislocation mobility


instance = phase_plasticityInstance(material_phase(g,ip,el))
ns = totalNslip(instance)

tauThreshold = state%p(iTauF(1:ns,instance))
tauEff = abs(tau) - tauThreshold

v = 0.0_pReal
if (present(dv_dtau)) dv_dtau = 0.0_pReal


if (Temperature > 0.0_pReal) then
  do s = 1_pInt,ns
    if (tauEff(s) > 0.0_pReal) then
      
      !* Peierls contribution
      !* The derivative only gives absolute values; the correct sign is taken care of in the formula for the derivative of the velocity
      
      meanfreepath_P = burgers(s,instance)
      jumpWidth_P = burgers(s,instance)
      activationLength_P = doublekinkwidth(instance) * burgers(s,instance)
      activationVolume_P = activationLength_P * jumpWidth_P * burgers(s,instance)
      criticalStress_P = peierlsStress(s,c,instance)
      activationEnergy_P = criticalStress_P * activationVolume_P
      tauRel_P = min(1.0_pReal, tauEff(s) / criticalStress_P)       ! ensure that the activation probability cannot become greater than one
      tPeierls = 1.0_pReal / fattack(instance) &
               * exp(activationEnergy_P / (KB * Temperature) &
                     * (1.0_pReal - tauRel_P**pParam(instance))**qParam(instance))
      if (present(dv_dtau)) then
        if (tauEff(s) < criticalStress_P) then
          dtPeierls_dtau = tPeierls * pParam(instance) * qParam(instance) * activationVolume_P / (KB * Temperature) &
                         * (1.0_pReal - tauRel_P**pParam(instance))**(qParam(instance)-1.0_pReal) &
                                      * tauRel_P**(pParam(instance)-1.0_pReal) 
        else
          dtPeierls_dtau = 0.0_pReal
        endif
      endif


      !* Contribution from solid solution strengthening
      !* The derivative only gives absolute values; the correct sign is taken care of in the formula for the derivative of the velocity

      meanfreepath_S = burgers(s,instance) / sqrt(solidSolutionConcentration(instance))
      jumpWidth_S = solidSolutionSize(instance) * burgers(s,instance)
      activationLength_S = burgers(s,instance) / sqrt(solidSolutionConcentration(instance))
      activationVolume_S = activationLength_S * jumpWidth_S * burgers(s,instance)
      activationEnergy_S = solidSolutionEnergy(instance)
      criticalStress_S = activationEnergy_S / activationVolume_S
      tauRel_S = min(1.0_pReal, tauEff(s) / criticalStress_S)       ! ensure that the activation probability cannot become greater than one
      tSolidSolution = 1.0_pReal / fattack(instance) &
                     * exp(activationEnergy_S / (KB * Temperature) &
                           * (1.0_pReal - tauRel_S**pParam(instance))**qParam(instance))
      if (present(dv_dtau)) then
        if (tauEff(s) < criticalStress_S) then
          dtSolidSolution_dtau = tSolidSolution * pParam(instance) * qParam(instance) &
                               * activationVolume_S / (KB * Temperature) &
                               * (1.0_pReal - tauRel_S**pParam(instance))**(qParam(instance)-1.0_pReal) &
                                              * tauRel_S**(pParam(instance)-1.0_pReal) 
        else
          dtSolidSolution_dtau = 0.0_pReal
        endif
      endif


      !* viscous glide velocity
      
      mobility = burgers(s,instance) / viscosity(instance)
      vViscous = mobility * tauEff(s)


      !* Mean velocity results from waiting time at peierls barriers and solid solution obstacles with respective meanfreepath of 
      !* free flight at glide velocity in between.       
      !* adopt sign from resolved stress

      v(s) = sign(1.0_pReal,tau(s)) &
           / (tPeierls / meanfreepath_P + tSolidSolution / meanfreepath_S + 1.0_pReal / vViscous)
      if (present(dv_dtau)) then
        dv_dtau(s) = v(s) * v(s) &
                   * (dtPeierls_dtau / meanfreepath_P &
                      + dtSolidSolution_dtau / meanfreepath_S &
                      + 1.0_pReal / (mobility * tauEff(s)*tauEff(s)))
      endif


    endif
  enddo
endif
    

#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
    write(6,*)
    write(6,'(a,i8,1x,i2,1x,i1)') '<< CONST >> nonlocal_kinetics at el ip g',el,ip,g
    write(6,*)
    write(6,'(a,/,12x,12(f12.5,1x))') '<< CONST >> tau / MPa', tau / 1e6_pReal
    write(6,'(a,/,12x,12(f12.5,1x))') '<< CONST >> tauEff / MPa', tauEff / 1e6_pReal
    write(6,'(a,/,12x,12(f12.5,1x))') '<< CONST >> v / 1e-3m/s', v * 1e3
    if (present(dv_dtau)) then
      write(6,'(a,/,12x,12(e12.5,1x))') '<< CONST >> dv_dtau', dv_dtau
    endif
  endif
#endif

endsubroutine



!*********************************************************************
!* calculates plastic velocity gradient and its tangent              *
!*********************************************************************
subroutine constitutive_nonlocal_LpAndItsTangent(Lp, dLp_dTstar99, Tstar_v, Temperature, state, g, ip, el)

use math,     only: math_Plain3333to99, &
                    math_mul6x6, &
                    math_Mandel6to33
use debug,    only: debug_level, &
                    debug_constitutive, &
                    debug_levelBasic, &
                    debug_levelExtensive, &
                    debug_levelSelective, &
                    debug_g, &
                    debug_i, &
                    debug_e
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_plasticityInstance
use lattice,  only: lattice_Sslip, &
                    lattice_Sslip_v, &
                    NnonSchmid
use mesh,     only: mesh_ipVolume

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
integer(pInt)                               myInstance, &               ! current instance of this plasticity
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
real(pReal), dimension(3,3,2,totalNslip(phase_plasticityInstance(material_phase(g,ip,el)))) :: & 
                                            nonSchmidTensor
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),8) :: &
                                            rhoSgl                      ! single dislocation densities (including blocked) 
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),4) :: &
                                            v, &                        ! velocity
                                            tau, &                      ! resolved shear stress including non Schmid and backstress terms
                                            dgdot_dtau, &               ! derivative of the shear rate with respect to the shear stress
                                            dv_dtau                     ! velocity derivative with respect to the shear stress
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el)))) :: &
                                            gdotTotal, &                ! shear rate
                                            tauBack, &                  ! back stress from dislocation gradients on same slip system
                                            deadZoneSize


!*** initialize local variables

Lp = 0.0_pReal
dLp_dTstar3333 = 0.0_pReal
nonSchmidTensor = 0.0_pReal

myInstance = phase_plasticityInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = totalNslip(myInstance)


!*** shortcut to state variables 


forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)
  rhoSgl(s,t) = max(state%p(iRhoU(s,t,myInstance)), 0.0_pReal)                              ! ensure positive single mobile densities
  rhoSgl(s,t+4_pInt) = state%p(iRhoB(s,t,myInstance))
endforall
where (abs(rhoSgl) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(myInstance) &
  .or. abs(rhoSgl) < significantRho(myInstance)) &
  rhoSgl = 0.0_pReal

tauBack = state%p(iTauB(1:ns,myInstance))


!*** get effective resolved shear stress
!*** add non schmid contributions to ONLY screw components if present (i.e.  if NnonSchmid(myStructure) > 0)

do s = 1_pInt,ns
  sLattice = slipSystemLattice(s,myInstance)  
  tau(s,1:4) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,1,sLattice,myStructure)) + tauBack(s)
  nonSchmidTensor(1:3,1:3,1,s) = lattice_Sslip(1:3,1:3,sLattice,myStructure)
  nonSchmidTensor(1:3,1:3,2,s) = nonSchmidTensor(1:3,1:3,1,s)
  do k = 1_pInt, NnonSchmid(myStructure)
    tau(s,3) = tau(s,3) + nonSchmidCoeff(k,myInstance) &
                        * math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,2*k,sLattice,myStructure))
    tau(s,4) = tau(s,4) + nonSchmidCoeff(k,myInstance) &
                        * math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,2*k+1,sLattice,myStructure))
    nonSchmidTensor(1:3,1:3,1,s) = nonSchmidTensor(1:3,1:3,1,s) &
                                 + nonSchmidCoeff(k,myInstance) &
                                 * math_Mandel6to33(lattice_Sslip_v(1:6,2*k,sLattice,myStructure))
    nonSchmidTensor(1:3,1:3,2,s) = nonSchmidTensor(1:3,1:3,2,s) &
                                 + nonSchmidCoeff(k,myInstance) &
                                 * math_Mandel6to33(lattice_Sslip_v(1:6,2*k+1,sLattice,myStructure))
  enddo
enddo


!*** get dislocation velocity and its tangent and store the velocity in the state array

if (myStructure == 1_pInt .and. NnonSchmid(myStructure) == 0_pInt) then                            ! for fcc all velcities are equal
  call constitutive_nonlocal_kinetics(v(1:ns,1), tau(1:ns,1), 1_pInt, Temperature, state, &
                                      g, ip, el, dv_dtau(1:ns,1))
  do t = 1_pInt,4_pInt
    v(1:ns,t) = v(1:ns,1)
    dv_dtau(1:ns,t) = dv_dtau(1:ns,1)
    state%p(iV(1:ns,t,myInstance)) = v(1:ns,1)
  enddo
else                                                                                               ! for all other lattice structures the velocities may vary with character and sign
  do t = 1_pInt,4_pInt
    c = (t-1_pInt)/2_pInt+1_pInt
    call constitutive_nonlocal_kinetics(v(1:ns,t), tau(1:ns,t), c, Temperature, state, &
                                        g, ip, el, dv_dtau(1:ns,t))
    state%p(iV(1:ns,t,myInstance)) = v(1:ns,t)
  enddo
endif


!*** Bauschinger effect

forall (s = 1_pInt:ns, t = 5_pInt:8_pInt, rhoSgl(s,t) * v(s,t-4_pInt) < 0.0_pReal) &
  rhoSgl(s,t-4_pInt) = rhoSgl(s,t-4_pInt) + abs(rhoSgl(s,t))


!*** Calculation of gdot and its tangent

deadZoneSize = 0.0_pReal
if (deadZoneScaling(myInstance)) then
  forall(s = 1_pInt:ns, sum(abs(rhoSgl(s,1:8))) > 0.0_pReal) &
    deadZoneSize(s) = maxval(abs(rhoSgl(s,5:8)) / (rhoSgl(s,1:4) + abs(rhoSgl(s,5:8))))
endif
gdotTotal = sum(rhoSgl(1:ns,1:4) * v, 2) * burgers(1:ns,myInstance) * (1.0_pReal - deadZoneSize)
do t = 1_pInt,4_pInt
  dgdot_dtau(:,t) = rhoSgl(1:ns,t) * dv_dtau(1:ns,t) * burgers(1:ns,myInstance) * (1.0_pReal - deadZoneSize)
enddo


!*** Calculation of Lp and its tangent

do s = 1_pInt,ns
  sLattice = slipSystemLattice(s,myInstance)  
  Lp = Lp + gdotTotal(s) * lattice_Sslip(1:3,1:3,sLattice,myStructure)
  forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt,k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
    dLp_dTstar3333(i,j,k,l) = dLp_dTstar3333(i,j,k,l) &
      + dgdot_dtau(s,1) * lattice_Sslip(i,j,sLattice,myStructure) * lattice_Sslip(k,l,sLattice,myStructure) &
      + dgdot_dtau(s,2) * lattice_Sslip(i,j,sLattice,myStructure) * lattice_Sslip(k,l,sLattice,myStructure) & 
      + dgdot_dtau(s,3) * lattice_Sslip(i,j,sLattice,myStructure) * nonSchmidTensor(k,l,1,s) &
      + dgdot_dtau(s,4) * lattice_Sslip(i,j,sLattice,myStructure) * nonSchmidTensor(k,l,2,s)
enddo
dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)


#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt )) then
    write(6,*)
    write(6,'(a,i8,1x,i2,1x,i1)') '<< CONST >> nonlocal_LpandItsTangent at el ip g ',el,ip,g
    write(6,*)
    write(6,'(a,/,12x,12(f12.5,1x))') '<< CONST >> gdot total / 1e-3',gdotTotal*1e3_pReal
    write(6,'(a,/,3(12x,3(f12.7,1x),/))') '<< CONST >> Lp',transpose(Lp)
  endif
#endif

endsubroutine



!*********************************************************************
!* incremental change of microstructure                              *
!*********************************************************************
subroutine constitutive_nonlocal_deltaState(deltaState, state, Tstar_v, Temperature, g,ip,el)

use debug,    only: debug_level, &
                    debug_constitutive, &
                    debug_levelBasic, &
                    debug_levelExtensive, &
                    debug_levelSelective, &
                    debug_g, &
                    debug_i, &
                    debug_e
use math,     only: pi, &
                    math_mul6x6
use lattice,  only: lattice_Sslip_v
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_ipVolume
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_plasticityInstance

implicit none

!*** input variables
integer(pInt), intent(in) ::                g, &                      ! current grain number
                                            ip, &                     ! current integration point
                                            el                        ! current element number
real(pReal), intent(in) ::                  Temperature               ! temperature
real(pReal), dimension(6), intent(in) ::    Tstar_v                   ! current 2nd Piola-Kirchhoff stress in Mandel notation

!*** input/output variables
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(inout) :: &
                                            state                     ! current microstructural state

!*** output variables
type(p_vec), intent(out) ::                 deltaState                ! change of state variables / microstructure
 
!*** local variables
integer(pInt)                               myInstance, &             ! current instance of this plasticity
                                            myStructure, &            ! current lattice structure
                                            ns, &                     ! short notation for the total number of active slip systems
                                            c, &                      ! character of dislocation
                                            t, &                      ! type of dislocation
                                            s, &                      ! index of my current slip system
                                            sLattice                  ! index of my current slip system according to lattice order
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),10) :: &
                                            deltaRho, &                     ! density increment
                                            deltaRhoRemobilization, &       ! density increment by remobilization
                                            deltaRhoDipole2SingleStress     ! density increment by dipole dissociation (by stress change)
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),8) :: &
                                            rhoSgl                        ! current single dislocation densities (positive/negative screw and edge without dipoles)
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),4) :: &
                                            v                             ! dislocation glide velocity
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el)))) :: &
                                            tau, &                        ! current resolved shear stress
                                            tauBack                       ! current back stress from pileups on same slip system
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),2) :: &
                                            rhoDip, &                     ! current dipole dislocation densities (screw and edge dipoles)
                                            dLower, &                     ! minimum stable dipole distance for edges and screws
                                            dUpper, &                     ! current maximum stable dipole distance for edges and screws
                                            dUpperOld, &                  ! old maximum stable dipole distance for edges and screws
                                            deltaDUpper                   ! change in maximum stable dipole distance for edges and screws


#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
    write(6,*)
    write(6,'(a,i8,1x,i2,1x,i1)') '<< CONST >> nonlocal_deltaState at el ip g ',el,ip,g
    write(6,*)
  endif
#endif

myInstance = phase_plasticityInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = totalNslip(myInstance)


!*** shortcut to state variables 


forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)
  rhoSgl(s,t) = max(state(g,ip,el)%p(iRhoU(s,t,myInstance)), 0.0_pReal)                              ! ensure positive single mobile densities
  rhoSgl(s,t+4_pInt) = state(g,ip,el)%p(iRhoB(s,t,myInstance))
  v(s,t) = state(g,ip,el)%p(iV(s,t,myInstance))
endforall
forall (s = 1_pInt:ns, c = 1_pInt:2_pInt)
  rhoDip(s,c) = max(state(g,ip,el)%p(iRhoD(s,c,myInstance)), 0.0_pReal)                               ! ensure positive dipole densities
  dUpperOld(s,c) = state(g,ip,el)%p(iD(s,c,myInstance))
endforall
tauBack = state(g,ip,el)%p(iTauB(1:ns,myInstance))

where (abs(rhoSgl) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(myInstance) &
  .or. abs(rhoSgl) < significantRho(myInstance)) &
  rhoSgl = 0.0_pReal
where (abs(rhoDip) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(myInstance) &
  .or. abs(rhoDip) < significantRho(myInstance)) &
  rhoDip = 0.0_pReal




!****************************************************************************
!*** dislocation remobilization (bauschinger effect)

deltaRhoRemobilization = 0.0_pReal
do t = 1_pInt,4_pInt
  do s = 1_pInt,ns
    if (rhoSgl(s,t+4_pInt) * v(s,t) < 0.0_pReal) then
      deltaRhoRemobilization(s,t) = abs(rhoSgl(s,t+4_pInt))
      rhoSgl(s,t) = rhoSgl(s,t) + abs(rhoSgl(s,t+4_pInt))
      deltaRhoRemobilization(s,t+4_pInt) = - rhoSgl(s,t+4_pInt)
      rhoSgl(s,t+4_pInt) = 0.0_pReal
    endif
  enddo
enddo



!****************************************************************************
!*** calculate dipole formation and dissociation by stress change

!*** calculate limits for stable dipole height

do s = 1_pInt,ns
  sLattice = slipSystemLattice(s,myInstance)  
  tau(s) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,1,sLattice,myStructure)) + tauBack(s)
  if (abs(tau(s)) < 1.0e-15_pReal) tau(s) = 1.0e-15_pReal
enddo
dLower = minDipoleHeight(1:ns,1:2,myInstance)
dUpper(1:ns,1) = mu(myInstance) * burgers(1:ns,myInstance) &
               / (8.0_pReal * pi * (1.0_pReal - nu(myInstance)) * abs(tau))
dUpper(1:ns,2) = mu(myInstance) * burgers(1:ns,myInstance) / (4.0_pReal * pi * abs(tau))
forall (c = 1_pInt:2_pInt) &
  dUpper(1:ns,c) = min(1.0_pReal / sqrt(rhoSgl(1:ns,2*c-1) + rhoSgl(1:ns,2*c) & 
                       + abs(rhoSgl(1:ns,2*c+3)) + abs(rhoSgl(1:ns,2*c+4)) + rhoDip(1:ns,c)), &
                       dUpper(1:ns,c))
dUpper = max(dUpper,dLower)
deltaDUpper = dUpper - dUpperOld


!*** dissociation by stress increase

deltaRhoDipole2SingleStress = 0.0_pReal
forall (c=1_pInt:2_pInt, s=1_pInt:ns, deltaDUpper(s,c) < 0.0_pReal) &
  deltaRhoDipole2SingleStress(s,8_pInt+c) = rhoDip(s,c) * deltaDUpper(s,c) / (dUpperOld(s,c) - dLower(s,c))

forall (t=1_pInt:4_pInt) &
  deltaRhoDipole2SingleStress(1_pInt:ns,t) = -0.5_pReal * deltaRhoDipole2SingleStress(1_pInt:ns,(t-1_pInt)/2_pInt+9_pInt)

 

!*** store new maximum dipole height in state

forall (s = 1_pInt:ns, c = 1_pInt:2_pInt) &
  state(g,ip,el)%p(iD(s,c,myInstance)) = dUpper(s,c) 



!****************************************************************************
!*** assign the changes in the dislocation densities to deltaState

deltaRho = deltaRhoRemobilization &
         + deltaRhoDipole2SingleStress

deltaState%p = 0.0_pReal
forall (s = 1:ns, t = 1_pInt:4_pInt) 
  deltaState%p(iRhoU(s,t,myInstance)) = deltaRho(s,t)
  deltaState%p(iRhoB(s,t,myInstance)) = deltaRho(s,t+4_pInt)
endforall
forall (s = 1:ns, c = 1_pInt:2_pInt) &
  deltaState%p(iRhoD(s,c,myInstance)) = deltaRho(s,c+8_pInt)


#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt )) then
    write(6,'(a,/,8(12x,12(e12.5,1x),/))') '<< CONST >> dislocation remobilization', deltaRhoRemobilization(1:ns,1:8)
    write(6,'(a,/,10(12x,12(e12.5,1x),/))') '<< CONST >> dipole dissociation by stress increase', deltaRhoDipole2SingleStress
    write(6,*)
  endif
#endif

endsubroutine



!*********************************************************************
!* rate of change of microstructure                                  *
!*********************************************************************
function constitutive_nonlocal_dotState(Tstar_v, Fe, Fp, Temperature, state, state0, timestep, subfrac, g,ip,el)

use prec,     only: DAMASK_NaN
use numerics, only: numerics_integrationMode, &
                    numerics_timeSyncing
use IO,       only: IO_error
use debug,    only: debug_level, &
                    debug_constitutive, &
                    debug_levelBasic, &
                    debug_levelExtensive, &
                    debug_levelSelective, &
                    debug_g, &
                    debug_i, &
                    debug_e
use math,     only: math_norm3, &
                    math_mul6x6, &
                    math_mul3x3, &
                    math_mul33x3, &
                    math_mul33x33, &
                    math_inv33, &
                    math_det33, &
                    math_transpose33, &  
                    pi
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_element, &
                    mesh_ipNeighborhood, &
                    mesh_ipVolume, &
                    mesh_ipArea, &
                    mesh_ipAreaNormal, &
                    FE_NipNeighbors, &
                    FE_geomtype, &
                    FE_celltype
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_plasticityInstance, &
                    phase_localPlasticity, &
                    phase_plasticity
use lattice,  only: lattice_Sslip_v, &
                    lattice_sd, &
                    lattice_st

implicit none

!*** input variables
integer(pInt), intent(in) ::                g, &                      !< current grain number
                                            ip, &                     !< current integration point
                                            el                        !< current element number
real(pReal), intent(in) ::                  Temperature, &            !< temperature
                                            timestep                  !< substepped crystallite time increment
real(pReal), dimension(6), intent(in) ::    Tstar_v                   !< current 2nd Piola-Kirchhoff stress in Mandel notation
real(pReal), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            subfrac                   !< fraction of timestep at the beginning of the substepped crystallite time increment 
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            Fe, &                     !< elastic deformation gradient
                                            Fp                        !< plastic deformation gradient
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            state, &                  !< current microstructural state
                                            state0                    !< microstructural state at beginning of crystallite increment

!*** input/output variables
 
!*** output variables
real(pReal), dimension(constitutive_nonlocal_sizeDotState(phase_plasticityInstance(material_phase(g,ip,el)))) :: &
                                            constitutive_nonlocal_dotState !< evolution of state variables / microstructure
 
!*** local variables
integer(pInt)                               myInstance, &             !< current instance of this plasticity
                                            neighbor_instance, &      !< instance of my neighbor's plasticity
                                            myStructure, &            !< current lattice structure
                                            ns, &                     !< short notation for the total number of active slip systems
                                            c, &                      !< character of dislocation
                                            n, &                      !< index of my current neighbor
                                            neighbor_el, &            !< element number of my neighbor
                                            neighbor_ip, &            !< integration point of my neighbor
                                            neighbor_n, &             !< neighbor index pointing to me when looking from my neighbor
                                            opposite_neighbor, &      !< index of my opposite neighbor
                                            opposite_ip, &            !< ip of my opposite neighbor
                                            opposite_el, &            !< element index of my opposite neighbor
                                            opposite_n, &             !< neighbor index pointing to me when looking from my opposite neighbor
                                            t, &                      !< type of dislocation
                                            topp, &                   !< type of dislocation with opposite sign to t
                                            s, &                      !< index of my current slip system
                                            sLattice, &               !< index of my current slip system according to lattice order
                                            deads
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),10) :: &
                                            rhoDot, &                     !< density evolution
                                            rhoDotMultiplication, &       !< density evolution by multiplication
                                            rhoDotFlux, &                 !< density evolution by flux
                                            rhoDotSingle2DipoleGlide, &   !< density evolution by dipole formation (by glide)
                                            rhoDotAthermalAnnihilation, & !< density evolution by athermal annihilation
                                            rhoDotThermalAnnihilation     !< density evolution by thermal annihilation
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),8) :: &
                                            rhoSgl, &                     !< current single dislocation densities (positive/negative screw and edge without dipoles)
                                            rhoSglOriginal, &
                                            neighbor_rhoSgl, &            !< current single dislocation densities of neighboring ip (positive/negative screw and edge without dipoles)
                                            rhoSgl0, &                    !< single dislocation densities at start of cryst inc (positive/negative screw and edge without dipoles)
                                            my_rhoSgl                     !< single dislocation densities of central ip (positive/negative screw and edge without dipoles)
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),4) :: &
                                            v, &                          !< current dislocation glide velocity
                                            v0, &                         !< dislocation glide velocity at start of cryst inc
                                            my_v, &                       !< dislocation glide velocity of central ip
                                            neighbor_v, &                 !< dislocation glide velocity of enighboring ip
                                            gdot                          !< shear rates
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el)))) :: &
                                            rhoForest, &                  !< forest dislocation density
                                            tauThreshold, &               !< threshold shear stress
                                            tau, &                        !< current resolved shear stress
                                            tauBack, &                    !< current back stress from pileups on same slip system
                                            vClimb, &                     !< climb velocity of edge dipoles
                                            nSources
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),2) :: &
                                            rhoDip, &                     !< current dipole dislocation densities (screw and edge dipoles)
                                            rhoDipOriginal, & 
                                            dLower, &                     !< minimum stable dipole distance for edges and screws
                                            dUpper                        !< current maximum stable dipole distance for edges and screws
real(pReal), dimension(3,totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),4) :: &
                                            m                             !< direction of dislocation motion
real(pReal), dimension(3,3) ::              my_F, &                       !< my total deformation gradient
                                            neighbor_F, &                 !< total deformation gradient of my neighbor
                                            my_Fe, &                      !< my elastic deformation gradient
                                            neighbor_Fe, &                !< elastic deformation gradient of my neighbor
                                            Favg                          !< average total deformation gradient of me and my neighbor
real(pReal), dimension(3) ::                normal_neighbor2me, &         !< interface normal pointing from my neighbor to me in neighbor's lattice configuration
                                            normal_neighbor2me_defConf, & !< interface normal pointing from my neighbor to me in shared deformed configuration
                                            normal_me2neighbor, &         !< interface normal pointing from me to my neighbor in my lattice configuration
                                            normal_me2neighbor_defConf    !< interface normal pointing from me to my neighbor in shared deformed configuration
real(pReal)                                 area, &                       !< area of the current interface
                                            transmissivity, &             !< overall transmissivity of dislocation flux to neighboring material point
                                            lineLength, &                 !< dislocation line length leaving the current interface
                                            selfDiffusion, &              !< self diffusion
                                            rnd, & 
                                            meshlength
logical                                     considerEnteringFlux, &
                                            considerLeavingFlux

#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
    write(6,*)
    write(6,'(a,i8,1x,i2,1x,i1)') '<< CONST >> nonlocal_dotState at el ip g ',el,ip,g
    write(6,*)
  endif
#endif


myInstance = phase_plasticityInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = totalNslip(myInstance)

tau = 0.0_pReal
gdot = 0.0_pReal


!*** shortcut to state variables 


forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)
  rhoSgl(s,t) = max(state(g,ip,el)%p(iRhoU(s,t,myInstance)), 0.0_pReal)                             ! ensure positive single mobile densities
  rhoSgl(s,t+4_pInt) = state(g,ip,el)%p(iRhoB(s,t,myInstance))
  v(s,t) = state(g,ip,el)%p(iV(s,t,myInstance))
endforall
forall (s = 1_pInt:ns, c = 1_pInt:2_pInt)
  rhoDip(s,c) = max(state(g,ip,el)%p(iRhoD(s,c,myInstance)), 0.0_pReal)                             ! ensure positive dipole densities
endforall
rhoForest = state(g,ip,el)%p(iRhoF(1:ns,myInstance))
tauThreshold = state(g,ip,el)%p(iTauF(1:ns,myInstance))
tauBack = state(g,ip,el)%p(iTauB(1:ns,myInstance))

rhoSglOriginal = rhoSgl
rhoDipOriginal = rhoDip
where (abs(rhoSgl) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(myInstance) &
  .or. abs(rhoSgl) < significantRho(myInstance)) &
  rhoSgl = 0.0_pReal
where (abs(rhoDip) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(myInstance) &
  .or. abs(rhoDip) < significantRho(myInstance)) &
  rhoDip = 0.0_pReal

if (numerics_timeSyncing) then
  forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)
    rhoSgl0(s,t) = max(state0(g,ip,el)%p(iRhoU(s,t,myInstance)), 0.0_pReal)
    rhoSgl0(s,t+4_pInt) = state0(g,ip,el)%p(iRhoB(s,t,myInstance))
    v0(s,t) = state0(g,ip,el)%p(iV(s,t,myInstance))
  endforall
  where (abs(rhoSgl0) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(myInstance) &
    .or. abs(rhoSgl0) < significantRho(myInstance)) &
    rhoSgl0 = 0.0_pReal
endif
  


!*** sanity check for timestep

if (timestep <= 0.0_pReal) then                                                                     ! if illegal timestep...
  constitutive_nonlocal_dotState = 0.0_pReal                                                        ! ...return without doing anything (-> zero dotState)
  return
endif



!****************************************************************************
!*** Calculate shear rate

forall (t = 1_pInt:4_pInt) &
  gdot(1_pInt:ns,t) = rhoSgl(1_pInt:ns,t) * burgers(1:ns,myInstance) * v(1:ns,t)

#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt )) then
    write(6,'(a,/,10(12x,12(e12.5,1x),/))') '<< CONST >> rho / 1/m^2', rhoSgl, rhoDip
    write(6,'(a,/,4(12x,12(e12.5,1x),/))') '<< CONST >> gdot / 1/s',gdot
  endif
#endif



!****************************************************************************
!*** calculate limits for stable dipole height

do s = 1_pInt,ns   ! loop over slip systems
  sLattice = slipSystemLattice(s,myInstance)  
  tau(s) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,1,sLattice,myStructure)) + tauBack(s)
  if (abs(tau(s)) < 1.0e-15_pReal) tau(s) = 1.0e-15_pReal
enddo

dLower = minDipoleHeight(1:ns,1:2,myInstance)
dUpper(1:ns,1) = mu(myInstance) * burgers(1:ns,myInstance) &
               / (8.0_pReal * pi * (1.0_pReal - nu(myInstance)) * abs(tau))
dUpper(1:ns,2) = mu(myInstance) * burgers(1:ns,myInstance) &
               / (4.0_pReal * pi * abs(tau))
forall (c = 1_pInt:2_pInt) &
  dUpper(1:ns,c) = min(1.0_pReal / sqrt(rhoSgl(1:ns,2*c-1) + rhoSgl(1:ns,2*c) & 
                       + abs(rhoSgl(1:ns,2*c+3)) + abs(rhoSgl(1:ns,2*c+4)) + rhoDip(1:ns,c)), &
                       dUpper(1:ns,c))
dUpper = max(dUpper,dLower)



!****************************************************************************
!*** calculate dislocation multiplication

rhoDotMultiplication = 0.0_pReal
if (probabilisticMultiplication(myInstance)) then
  meshlength = mesh_ipVolume(ip,el)**0.333_pReal
  where(sum(rhoSgl(1:ns,1:4),2) > 0.0_pReal)
    nSources = (sum(rhoSgl(1:ns,1:2),2) * fEdgeMultiplication(myInstance) + sum(rhoSgl(1:ns,3:4),2)) &
             / sum(rhoSgl(1:ns,1:4),2) * meshlength / lambda0(1:ns,myInstance) * sqrt(rhoForest(1:ns))
  elsewhere
    nSources = meshlength / lambda0(1:ns,myInstance) * sqrt(rhoForest(1:ns))
  endwhere
  do s = 1_pInt,ns
    if (nSources(s) < 1.0_pReal) then
      if (sourceProbability(s,g,ip,el) > 1.0_pReal) then
        call random_number(rnd)
        sourceProbability(s,g,ip,el) = rnd
        !$OMP FLUSH(sourceProbability)
      endif
      if (sourceProbability(s,g,ip,el) > 1.0_pReal - nSources(s)) then
        rhoDotMultiplication(s,1:4) = sum(rhoSglOriginal(s,1:4) * abs(v(s,1:4))) / meshlength
      endif
    else
      sourceProbability(s,g,ip,el) = 2.0_pReal
      rhoDotMultiplication(s,1:4) = &
        (sum(abs(gdot(s,1:2))) * fEdgeMultiplication(myInstance) + sum(abs(gdot(s,3:4)))) &
        / burgers(s,myInstance) * sqrt(rhoForest(s)) / lambda0(s,myInstance)
    endif
  enddo
#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt )) then
    write(6,'(a,/,4(12x,12(f12.5,1x),/))') '<< CONST >> sources', nSources
    write(6,*)
  endif
#endif
else
  rhoDotMultiplication(1:ns,1:4) = spread( &
      (sum(abs(gdot(1:ns,1:2)),2) * fEdgeMultiplication(myInstance) + sum(abs(gdot(1:ns,3:4)),2)) &
    * sqrt(rhoForest(1:ns)) / lambda0(1:ns,myInstance) / burgers(1:ns,myInstance), 2, 4)
endif



!****************************************************************************
!*** calculate dislocation fluxes (only for nonlocal plasticity)

rhoDotFlux = 0.0_pReal

if (.not. phase_localPlasticity(material_phase(g,ip,el))) then                                      ! only for nonlocal plasticity
  

  !*** check CFL (Courant-Friedrichs-Lewy) condition for flux

  if (any( abs(gdot) > 0.0_pReal &                                                                  ! any active slip system ...
          .and. CFLfactor(myInstance) * abs(v) * timestep &
              > mesh_ipVolume(ip,el) / maxval(mesh_ipArea(:,ip,el)))) then                          ! ...with velocity above critical value (we use the reference volume and area for simplicity here)
#ifndef _OPENMP
    if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt) then 
      write(6,'(a,i5,a,i2)') '<< CONST >> CFL condition not fullfilled at el ',el,' ip ',ip
      write(6,'(a,e10.3,a,e10.3)') '<< CONST >> velocity is at  ', &
        maxval(abs(v), abs(gdot) > 0.0_pReal &
                       .and. CFLfactor(myInstance) * abs(v) * timestep &
                             > mesh_ipVolume(ip,el) / maxval(mesh_ipArea(:,ip,el))), &
        ' at a timestep of ',timestep
      write(6,'(a)') '<< CONST >> enforcing cutback !!!'
    endif
#endif
    constitutive_nonlocal_dotState = DAMASK_NaN                                                     ! -> return NaN and, hence, enforce cutback
    return
  endif


  !*** be aware of the definition of lattice_st = lattice_sd x lattice_sn !!!
  !*** opposite sign to our p vector in the (s,p,n) triplet !!!
  
  m(1:3,1:ns,1) =  lattice_sd(1:3, slipSystemLattice(1:ns,myInstance), myStructure)
  m(1:3,1:ns,2) = -lattice_sd(1:3, slipSystemLattice(1:ns,myInstance), myStructure)
  m(1:3,1:ns,3) = -lattice_st(1:3, slipSystemLattice(1:ns,myInstance), myStructure)
  m(1:3,1:ns,4) =  lattice_st(1:3, slipSystemLattice(1:ns,myInstance), myStructure)
  
  my_Fe = Fe(1:3,1:3,g,ip,el)
  my_F = math_mul33x33(my_Fe, Fp(1:3,1:3,g,ip,el))
  
  do n = 1_pInt,FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,el))))                       ! loop through my neighbors
    neighbor_el = mesh_ipNeighborhood(1,n,ip,el)
    neighbor_ip = mesh_ipNeighborhood(2,n,ip,el)
    neighbor_n  = mesh_ipNeighborhood(3,n,ip,el)
  
    opposite_neighbor = n + mod(n,2_pInt) - mod(n+1_pInt,2_pInt)
    opposite_el = mesh_ipNeighborhood(1,opposite_neighbor,ip,el)
    opposite_ip = mesh_ipNeighborhood(2,opposite_neighbor,ip,el)
    opposite_n  = mesh_ipNeighborhood(3,opposite_neighbor,ip,el)
  
    if (neighbor_n > 0_pInt) then                                                                   ! if neighbor exists, average deformation gradient
      neighbor_instance = phase_plasticityInstance(material_phase(g,neighbor_ip,neighbor_el))
      neighbor_Fe = Fe(1:3,1:3,g,neighbor_ip,neighbor_el)
      neighbor_F = math_mul33x33(neighbor_Fe, Fp(1:3,1:3,g,neighbor_ip,neighbor_el))
      Favg = 0.5_pReal * (my_F + neighbor_F)
    else                                                                                            ! if no neighbor, take my value as average
      Favg = my_F
    endif
    

    !* FLUX FROM MY NEIGHBOR TO ME
    !* This is only considered, if I have a neighbor of nonlocal plasticity (also nonlocal constitutive law with local properties) that is at least a little bit compatible.
    !* If it's not at all compatible, no flux is arriving, because everything is dammed in front of my neighbor's interface.
    !* The entering flux from my neighbor will be distributed on my slip systems according to the compatibility
    
    considerEnteringFlux = .false.
    neighbor_v = 0.0_pReal        ! needed for check of sign change in flux density below 
    neighbor_rhoSgl = 0.0_pReal
    if (neighbor_n > 0_pInt) then
      if (phase_plasticity(material_phase(1,neighbor_ip,neighbor_el)) == CONSTITUTIVE_NONLOCAL_LABEL &
          .and. any(compatibility(:,:,:,n,ip,el) > 0.0_pReal)) &
        considerEnteringFlux = .true.
    endif
    
    if (considerEnteringFlux) then
      if(numerics_timeSyncing .and. (subfrac(g,neighbor_ip,neighbor_el) /= subfrac(g,ip,el))) then  ! for timesyncing: in case of a timestep at the interface we have to use "state0" to make sure that fluxes n both sides are equal
        forall (s = 1:ns, t = 1_pInt:4_pInt)
          neighbor_v(s,t) = state0(g,neighbor_ip,neighbor_el)%p(iV(s,t,neighbor_instance))
          neighbor_rhoSgl(s,t) = max(state0(g,neighbor_ip,neighbor_el)%p(iRhoU(s,t,neighbor_instance)), 0.0_pReal)
          neighbor_rhoSgl(s,t+4_pInt) = state0(g,neighbor_ip,neighbor_el)%p(iRhoB(s,t,neighbor_instance))
        endforall
      else
        forall (s = 1:ns, t = 1_pInt:4_pInt)
          neighbor_v(s,t) = state(g,neighbor_ip,neighbor_el)%p(iV(s,t,neighbor_instance))
          neighbor_rhoSgl(s,t) = max(state(g,neighbor_ip,neighbor_el)%p(iRhoU(s,t,neighbor_instance)), 0.0_pReal)
          neighbor_rhoSgl(s,t+4_pInt) = state(g,neighbor_ip,neighbor_el)%p(iRhoB(s,t,neighbor_instance))
        endforall
      endif
      where (abs(neighbor_rhoSgl) * mesh_ipVolume(neighbor_ip,neighbor_el) ** 0.667_pReal &
             < significantN(myInstance) &
        .or. abs(neighbor_rhoSgl) < significantRho(myInstance)) &
        neighbor_rhoSgl = 0.0_pReal
      normal_neighbor2me_defConf = math_det33(Favg) * math_mul33x3(math_inv33(transpose(Favg)), &
                                   mesh_ipAreaNormal(1:3,neighbor_n,neighbor_ip,neighbor_el))       ! calculate the normal of the interface in (average) deformed configuration (now pointing from my neighbor to me!!!)
      normal_neighbor2me = math_mul33x3(transpose(neighbor_Fe), normal_neighbor2me_defConf) &
                         / math_det33(neighbor_Fe)                                                  ! interface normal in the lattice configuration of my neighbor
      area = mesh_ipArea(neighbor_n,neighbor_ip,neighbor_el) * math_norm3(normal_neighbor2me)
      normal_neighbor2me = normal_neighbor2me / math_norm3(normal_neighbor2me)                      ! normalize the surface normal to unit length
      do s = 1_pInt,ns
        do t = 1_pInt,4_pInt
          c = (t + 1_pInt) / 2
          topp = t + mod(t,2_pInt) - mod(t+1_pInt,2_pInt)
          if (neighbor_v(s,t) * math_mul3x3(m(1:3,s,t), normal_neighbor2me) > 0.0_pReal &           ! flux from my neighbor to me == entering flux for me
              .and. v(s,t) * neighbor_v(s,t) > 0.0_pReal ) then                                     ! ... only if no sign change in flux density  
            do deads = 0_pInt,4_pInt,4_pInt
              lineLength = abs(neighbor_rhoSgl(s,t+deads)) * neighbor_v(s,t) &
                         * math_mul3x3(m(1:3,s,t), normal_neighbor2me) * area                       ! positive line length that wants to enter through this interface
              where (compatibility(c,1_pInt:ns,s,n,ip,el) > 0.0_pReal) &                            ! positive compatibility...
                rhoDotFlux(1_pInt:ns,t) = rhoDotFlux(1_pInt:ns,t) &
                                        + lineLength / mesh_ipVolume(ip,el) &                       ! ... transferring to equally signed mobile dislocation type
                                        * compatibility(c,1_pInt:ns,s,n,ip,el) ** 2.0_pReal
              where (compatibility(c,1_pInt:ns,s,n,ip,el) < 0.0_pReal) &                            ! ..negative compatibility...
                rhoDotFlux(1_pInt:ns,topp) = rhoDotFlux(1_pInt:ns,topp) &
                                           + lineLength / mesh_ipVolume(ip,el) &                    ! ... transferring to opposite signed mobile dislocation type
                                           * compatibility(c,1_pInt:ns,s,n,ip,el) ** 2.0_pReal
            enddo
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
    if (opposite_n > 0_pInt) then
      if (phase_plasticity(material_phase(1,opposite_ip,opposite_el)) /= CONSTITUTIVE_NONLOCAL_LABEL) &
        considerLeavingFlux = .false.
    endif

    if (considerLeavingFlux) then
      
      !* timeSyncing mode: If the central ip has zero subfraction, always use "state0". This is needed in case of 
      !* a synchronization step for the central ip, because then "state" contains the values at the end of the
      !* previously converged full time step. Also, if either me or my neighbor has zero subfraction, we have to 
      !* use "state0" to make sure that fluxes on both sides of the (potential) timestep are equal.
      my_rhoSgl = rhoSgl
      my_v = v
      if(numerics_timeSyncing) then
        if (subfrac(g,ip,el) == 0.0_pReal) then
          my_rhoSgl = rhoSgl0
          my_v = v0
        elseif (neighbor_n > 0_pInt) then
          if (subfrac(g,neighbor_ip,neighbor_el) == 0.0_pReal) then
            my_rhoSgl = rhoSgl0
            my_v = v0
          endif
        endif
      endif
      
      normal_me2neighbor_defConf = math_det33(Favg) &
                                 * math_mul33x3(math_inv33(math_transpose33(Favg)), & 
                                                           mesh_ipAreaNormal(1:3,n,ip,el))          ! calculate the normal of the interface in (average) deformed configuration (pointing from me to my neighbor!!!)
      normal_me2neighbor = math_mul33x3(math_transpose33(my_Fe), normal_me2neighbor_defConf) &
                         / math_det33(my_Fe)                                                        ! interface normal in my lattice configuration
      area = mesh_ipArea(n,ip,el) * math_norm3(normal_me2neighbor)
      normal_me2neighbor = normal_me2neighbor / math_norm3(normal_me2neighbor)                      ! normalize the surface normal to unit length    
      do s = 1_pInt,ns
        do t = 1_pInt,4_pInt
          c = (t + 1_pInt) / 2_pInt
          if (my_v(s,t) * math_mul3x3(m(1:3,s,t), normal_me2neighbor) > 0.0_pReal ) then            ! flux from me to my neighbor == leaving flux for me (might also be a pure flux from my mobile density to dead density if interface not at all transmissive)
            if (my_v(s,t) * neighbor_v(s,t) > 0.0_pReal) then                                       ! no sign change in flux density
              transmissivity = sum(compatibility(c,1_pInt:ns,s,n,ip,el)**2.0_pReal)                 ! overall transmissivity from this slip system to my neighbor
            else                                                                                    ! sign change in flux density means sign change in stress which does not allow for dislocations to arive at the neighbor
              transmissivity = 0.0_pReal
            endif
            lineLength = my_rhoSgl(s,t) * my_v(s,t) &
                       * math_mul3x3(m(1:3,s,t), normal_me2neighbor) * area                         ! positive line length of mobiles that wants to leave through this interface
            rhoDotFlux(s,t) = rhoDotFlux(s,t) - lineLength / mesh_ipVolume(ip,el)                   ! subtract dislocation flux from current type
            rhoDotFlux(s,t+4_pInt) = rhoDotFlux(s,t+4_pInt) &
                                   + lineLength / mesh_ipVolume(ip,el) * (1.0_pReal - transmissivity) &
                                   * sign(1.0_pReal, my_v(s,t))                                     ! dislocation flux that is not able to leave through interface (because of low transmissivity) will remain as immobile single density at the material point
            lineLength = my_rhoSgl(s,t+4_pInt) * my_v(s,t) &
                       * math_mul3x3(m(1:3,s,t), normal_me2neighbor) * area                         ! positive line length of deads that wants to leave through this interface
            rhoDotFlux(s,t+4_pInt) = rhoDotFlux(s,t+4_pInt) &
                                   - lineLength / mesh_ipVolume(ip,el) * transmissivity             ! dead dislocations leaving through this interface
          endif
        enddo
      enddo
    endif    
    
  enddo ! neighbor loop  
endif



!****************************************************************************
!*** calculate dipole formation and annihilation

!*** formation by glide

do c = 1_pInt,2_pInt

  rhoDotSingle2DipoleGlide(1:ns,2*c-1) = -2.0_pReal * dUpper(1:ns,c) / burgers(1:ns,myInstance) &
                                                    * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) &                                   ! negative mobile --> positive mobile
                                                       + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1)) &                                 ! positive mobile --> negative mobile
                                                       + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1)))                           ! positive mobile --> negative immobile

  rhoDotSingle2DipoleGlide(1:ns,2*c) = -2.0_pReal * dUpper(1:ns,c) / burgers(1:ns,myInstance) &
                                                  * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) &                                     ! negative mobile --> positive mobile
                                                     + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1)) &                                   ! positive mobile --> negative mobile
                                                     + abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c)))                               ! negative mobile --> positive immobile

  rhoDotSingle2DipoleGlide(1:ns,2*c+3) = -2.0_pReal * dUpper(1:ns,c) / burgers(1:ns,myInstance) &
                                                    * rhoSgl(1:ns,2*c+3) * abs(gdot(1:ns,2*c))                                      ! negative mobile --> positive immobile

  rhoDotSingle2DipoleGlide(1:ns,2*c+4) = -2.0_pReal * dUpper(1:ns,c) / burgers(1:ns,myInstance) &
                                                    * rhoSgl(1:ns,2*c+4) * abs(gdot(1:ns,2*c-1))                                    ! positive mobile --> negative immobile

  rhoDotSingle2DipoleGlide(1:ns,c+8) = - rhoDotSingle2DipoleGlide(1:ns,2*c-1) - rhoDotSingle2DipoleGlide(1:ns,2*c) &
                                       + abs(rhoDotSingle2DipoleGlide(1:ns,2*c+3)) + abs(rhoDotSingle2DipoleGlide(1:ns,2*c+4))
enddo


!*** athermal annihilation

rhoDotAthermalAnnihilation = 0.0_pReal

forall (c=1_pInt:2_pInt) &  
  rhoDotAthermalAnnihilation(1:ns,c+8_pInt) = -2.0_pReal * dLower(1:ns,c) / burgers(1:ns,myInstance) &
               * (  2.0_pReal * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1))) &             ! was single hitting single
                  + 2.0_pReal * (abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c)) + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1))) & ! was single hitting immobile single or was immobile single hit by single
                  + rhoDip(1:ns,c) * (abs(gdot(1:ns,2*c-1)) + abs(gdot(1:ns,2*c))))                                                 ! single knocks dipole constituent
! annihilated screw dipoles leave edge jogs behind on the colinear system
if (myStructure == 1_pInt) then ! only fcc
  forall (s = 1:ns, colinearSystem(s,myInstance) > 0_pInt) &
    rhoDotAthermalAnnihilation(colinearSystem(s,myInstance),1:2) = - rhoDotAthermalAnnihilation(s,10) &
      * 0.25_pReal * sqrt(rhoForest(s)) * (dUpper(s,2) + dLower(s,2)) * edgeJogFactor(myInstance)
endif
  
  
!*** thermally activated annihilation of edge dipoles by climb

rhoDotThermalAnnihilation = 0.0_pReal
selfDiffusion = Dsd0(myInstance) * exp(-selfDiffusionEnergy(myInstance) / (KB * Temperature))
vClimb =  atomicVolume(myInstance) * selfDiffusion / ( KB * Temperature ) &
          * mu(myInstance) / ( 2.0_pReal * PI * (1.0_pReal-nu(myInstance)) ) &
          * 2.0_pReal / ( dUpper(1:ns,1) + dLower(1:ns,1) )
forall (s = 1_pInt:ns, dUpper(s,1) > dLower(s,1)) &
  rhoDotThermalAnnihilation(s,9) = max(- 4.0_pReal * rhoDip(s,1) * vClimb(s) / (dUpper(s,1) - dLower(s,1)), &
                                       - rhoDip(s,1) / timestep - rhoDotAthermalAnnihilation(s,9) - rhoDotSingle2DipoleGlide(s,9))  ! make sure that we do not annihilate more dipoles than we have



!****************************************************************************
!*** assign the rates of dislocation densities to my dotState
!*** if evolution rates lead to negative densities, a cutback is enforced

rhoDot = 0.0_pReal
rhoDot = rhoDotFlux &
       + rhoDotMultiplication &
       + rhoDotSingle2DipoleGlide &
       + rhoDotAthermalAnnihilation &
       + rhoDotThermalAnnihilation 

if (numerics_integrationMode == 1_pInt) then                                       ! save rates for output if in central integration mode
  rhoDotFluxOutput(1:ns,1:8,g,ip,el) = rhoDotFlux(1:ns,1:8)
  rhoDotMultiplicationOutput(1:ns,1:2,g,ip,el) = rhoDotMultiplication(1:ns,[1,3])
  rhoDotSingle2DipoleGlideOutput(1:ns,1:2,g,ip,el) = rhoDotSingle2DipoleGlide(1:ns,9:10)
  rhoDotAthermalAnnihilationOutput(1:ns,1:2,g,ip,el) = rhoDotAthermalAnnihilation(1:ns,9:10)
  rhoDotThermalAnnihilationOutput(1:ns,1:2,g,ip,el) = rhoDotThermalAnnihilation(1:ns,9:10)
  rhoDotEdgeJogsOutput(1:ns,g,ip,el) = 2.0_pReal * rhoDotThermalAnnihilation(1:ns,1)
endif


#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt )) then
    write(6,'(a,/,4(12x,12(e12.5,1x),/))') '<< CONST >> dislocation multiplication', rhoDotMultiplication(1:ns,1:4) * timestep
    write(6,'(a,/,8(12x,12(e12.5,1x),/))') '<< CONST >> dislocation flux', rhoDotFlux(1:ns,1:8) * timestep
    write(6,'(a,/,10(12x,12(e12.5,1x),/))') '<< CONST >> dipole formation by glide', rhoDotSingle2DipoleGlide * timestep
    write(6,'(a,/,10(12x,12(e12.5,1x),/))') '<< CONST >> athermal dipole annihilation', &
                                            rhoDotAthermalAnnihilation * timestep
    write(6,'(a,/,2(12x,12(e12.5,1x),/))') '<< CONST >> thermally activated dipole annihilation', &
                                            rhoDotThermalAnnihilation(1:ns,9:10) * timestep
    write(6,'(a,/,10(12x,12(e12.5,1x),/))') '<< CONST >> total density change', rhoDot * timestep
    write(6,'(a,/,10(12x,12(f12.5,1x),/))') '<< CONST >> relative density change', &
                                            rhoDot(1:ns,1:8) * timestep / (abs(rhoSglOriginal)+1.0e-10), &
                                            rhoDot(1:ns,9:10) * timestep / (rhoDipOriginal+1.0e-10)
    write(6,*)
  endif
#endif


if (    any(rhoSglOriginal(1:ns,1:4) + rhoDot(1:ns,1:4) * timestep < -aTolRho(myInstance)) &
   .or. any(rhoDipOriginal(1:ns,1:2) + rhoDot(1:ns,9:10) * timestep < -aTolRho(myInstance))) then
#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt) then 
    write(6,'(a,i5,a,i2)') '<< CONST >> evolution rate leads to negative density at el ',el,' ip ',ip
    write(6,'(a)') '<< CONST >> enforcing cutback !!!'
  endif
#endif
  constitutive_nonlocal_dotState = DAMASK_NaN
  return
else
  forall (s = 1:ns, t = 1_pInt:4_pInt) 
    constitutive_nonlocal_dotState(iRhoU(s,t,myInstance)) = rhoDot(s,t)
    constitutive_nonlocal_dotState(iRhoB(s,t,myInstance)) = rhoDot(s,t+4_pInt)
  endforall
  forall (s = 1:ns, c = 1_pInt:2_pInt) &
    constitutive_nonlocal_dotState(iRhoD(s,c,myInstance)) = rhoDot(s,c+8_pInt)
  forall (s = 1:ns) &
    constitutive_nonlocal_dotState(iGamma(s,myInstance)) = sum(gdot(s,1:4))
endif

endfunction



!*********************************************************************
!* COMPATIBILITY UPDATE                                              *
!* Compatibility is defined as normalized product of signed cosine   *
!* of the angle between the slip plane normals and signed cosine of  *
!* the angle between the slip directions. Only the largest values    *
!* that sum up to a total of 1 are considered, all others are set to *
!* zero.                                                             *
!*********************************************************************
subroutine constitutive_nonlocal_updateCompatibility(orientation,i,e)

use math, only:       math_qDisorientation, &
                      math_mul3x3, &
                      math_qRot
use material, only:   material_phase, &
                      material_texture, &
                      phase_localPlasticity, &
                      phase_plasticityInstance, &
                      homogenization_maxNgrains
use mesh, only:       mesh_element, &
                      mesh_ipNeighborhood, &
                      mesh_maxNips, &
                      mesh_NcpElems, &
                      FE_NipNeighbors, &
                      FE_geomtype, &
                      FE_celltype
use lattice, only:    lattice_sn, &
                      lattice_sd

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
                                                neighbor_e, &                 ! element index of my neighbor
                                                neighbor_i, &                 ! integration point index of my neighbor
                                                my_phase, &
                                                neighbor_phase, &
                                                my_texture, &
                                                neighbor_texture, &
                                                my_structure, &               ! lattice structure
                                                my_instance, &                ! instance of plasticity
                                                ns, &                         ! number of active slip systems
                                                s1, &                         ! slip system index (me)
                                                s2                            ! slip system index (my neighbor)
real(pReal), dimension(4) ::                    absoluteMisorientation        ! absolute misorientation (without symmetry) between me and my neighbor
real(pReal), dimension(2,totalNslip(phase_plasticityInstance(material_phase(1,i,e))),&
                         totalNslip(phase_plasticityInstance(material_phase(1,i,e))),&
                         FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,e))))) :: &  
                                                my_compatibility              ! my_compatibility for current element and ip
real(pReal), dimension(3,totalNslip(phase_plasticityInstance(material_phase(1,i,e)))) :: &  
                                                slipNormal, &
                                                slipDirection
real(pReal)                                     my_compatibilitySum, &
                                                thresholdValue, &
                                                nThresholdValues
logical, dimension(totalNslip(phase_plasticityInstance(material_phase(1,i,e)))) :: & 
                                                belowThreshold


Nneighbors = FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,e))))
my_phase = material_phase(1,i,e)
my_texture = material_texture(1,i,e)
my_instance = phase_plasticityInstance(my_phase)
my_structure = constitutive_nonlocal_structure(my_instance)
ns = totalNslip(my_instance)
slipNormal(1:3,1:ns) = lattice_sn(1:3, slipSystemLattice(1:ns,my_instance), my_structure)
slipDirection(1:3,1:ns) = lattice_sd(1:3, slipSystemLattice(1:ns,my_instance), my_structure)


!*** start out fully compatible

my_compatibility = 0.0_pReal
forall(s1 = 1_pInt:ns) &
  my_compatibility(1:2,s1,s1,1:Nneighbors) = 1.0_pReal


!*** Loop thrugh neighbors and check whether there is any my_compatibility.

do n = 1_pInt,Nneighbors
  neighbor_e = mesh_ipNeighborhood(1,n,i,e)
  neighbor_i = mesh_ipNeighborhood(2,n,i,e)
  
  
  !* FREE SURFACE
  !* Set surface transmissivity to the value specified in the material.config
  
  if (neighbor_e <= 0_pInt .or. neighbor_i <= 0_pInt) then
    forall(s1 = 1_pInt:ns) &
      my_compatibility(1:2,s1,s1,n) = sqrt(surfaceTransmissivity(my_instance))
    cycle
  endif
  
  
  !* PHASE BOUNDARY
  !* If we encounter a different nonlocal "cpfem" phase at the neighbor, 
  !* we consider this to be a real "physical" phase boundary, so completely incompatible.
  !* If one of the two "CPFEM" phases has a local plasticity law, 
  !* we do not consider this to be a phase boundary, so completely compatible.
  
  neighbor_phase = material_phase(1,neighbor_i,neighbor_e)
  if (neighbor_phase /= my_phase) then
    if (.not. phase_localPlasticity(neighbor_phase) .and. .not. phase_localPlasticity(my_phase)) then
      forall(s1 = 1_pInt:ns) &
        my_compatibility(1:2,s1,s1,n) = 0.0_pReal ! = sqrt(0.0)
    endif
    cycle
  endif

    
  !* GRAIN BOUNDARY !
  !* fixed transmissivity for adjacent ips with different texture (only if explicitly given in material.config)

  if (grainboundaryTransmissivity(my_instance) >= 0.0_pReal) then
    neighbor_texture = material_texture(1,neighbor_i,neighbor_e)
    if (neighbor_texture /= my_texture) then
      if (.not. phase_localPlasticity(neighbor_phase)) then
        forall(s1 = 1_pInt:ns) &
          my_compatibility(1:2,s1,s1,n) = sqrt(grainboundaryTransmissivity(my_instance))
      endif
      cycle
    endif
  

  !* GRAIN BOUNDARY ?
  !* Compatibility defined by relative orientation of slip systems:
  !* The my_compatibility value is defined as the product of the slip normal projection and the slip direction projection.
  !* Its sign is always positive for screws, for edges it has the same sign as the slip normal projection. 
  !* Since the sum for each slip system can easily exceed one (which would result in a transmissivity larger than one), 
  !* only values above or equal to a certain threshold value are considered. This threshold value is chosen, such that
  !* the number of compatible slip systems is minimized with the sum of the original my_compatibility values exceeding one. 
  !* Finally the smallest my_compatibility value is decreased until the sum is exactly equal to one. 
  !* All values below the threshold are set to zero.
  else
    absoluteMisorientation = math_qDisorientation(orientation(1:4,1,i,e), &
                                                  orientation(1:4,1,neighbor_i,neighbor_e), &
                                                  0_pInt)      ! no symmetry
    do s1 = 1_pInt,ns    ! my slip systems
      do s2 = 1_pInt,ns  ! my neighbor's slip systems
        my_compatibility(1,s2,s1,n) =  math_mul3x3(slipNormal(1:3,s1), math_qRot(absoluteMisorientation, slipNormal(1:3,s2))) &
                                * abs(math_mul3x3(slipDirection(1:3,s1), math_qRot(absoluteMisorientation, slipDirection(1:3,s2))))
        my_compatibility(2,s2,s1,n) = abs(math_mul3x3(slipNormal(1:3,s1), math_qRot(absoluteMisorientation, slipNormal(1:3,s2)))) &
                                * abs(math_mul3x3(slipDirection(1:3,s1), math_qRot(absoluteMisorientation, slipDirection(1:3,s2))))
      enddo
      
      my_compatibilitySum = 0.0_pReal
      belowThreshold = .true.
      do while (my_compatibilitySum < 1.0_pReal .and. any(belowThreshold(1:ns)))
        thresholdValue = maxval(my_compatibility(2,1:ns,s1,n), belowThreshold(1:ns))              ! screws always positive
        nThresholdValues = real(count(my_compatibility(2,1:ns,s1,n) == thresholdValue),pReal)
        where (my_compatibility(2,1:ns,s1,n) >= thresholdValue) &
          belowThreshold(1:ns) = .false.
        if (my_compatibilitySum + thresholdValue * nThresholdValues > 1.0_pReal) &
          where (abs(my_compatibility(1:2,1:ns,s1,n)) == thresholdValue) &
            my_compatibility(1:2,1:ns,s1,n) = sign((1.0_pReal - my_compatibilitySum) &
                                                 / nThresholdValues, my_compatibility(1:2,1:ns,s1,n))
        my_compatibilitySum = my_compatibilitySum + nThresholdValues * thresholdValue
      enddo
      where (belowThreshold(1:ns)) my_compatibility(1,1:ns,s1,n) = 0.0_pReal
      where (belowThreshold(1:ns)) my_compatibility(2,1:ns,s1,n) = 0.0_pReal
    enddo ! my slip systems cycle
  endif

enddo   ! neighbor cycle

compatibility(1:2,1:ns,1:ns,1:Nneighbors,i,e) = my_compatibility

endsubroutine 



!*********************************************************************
!* rate of change of temperature                                     *
!*********************************************************************
pure function constitutive_nonlocal_dotTemperature(Tstar_v,Temperature,state,g,ip,el)

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

use math,     only: math_mul33x33, &
                    math_mul33x3, &
                    math_invert33, &
                    math_transpose33, &
                    pi
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_element, &
                    mesh_node0, &
                    mesh_cellCenterCoordinates, &
                    mesh_ipVolume, &
                    mesh_periodicSurface, &
                    FE_Nips, &
                    FE_geomtype
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_localPlasticity, &
                    phase_plasticityInstance

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
integer(pInt)                   neighbor_el, &                ! element number of neighbor material point
                                neighbor_ip, &                ! integration point of neighbor material point
                                instance, &                   ! my instance of this plasticity
                                neighbor_instance, &          ! instance of this plasticity of neighbor material point
                                latticeStruct, &              ! my lattice structure
                                neighbor_latticeStruct, &     ! lattice structure of neighbor material point
                                phase, &
                                neighbor_phase, &
                                ns, &                         ! total number of active slip systems at my material point
                                neighbor_ns, &                ! total number of active slip systems at neighbor material point
                                c, &                          ! index of dilsocation character (edge, screw)
                                s, &                          ! slip system index
                                t, &                          ! index of dilsocation type (e+, e-, s+, s-, used e+, used e-, used s+, used s-)
                                dir, &
                                deltaX, deltaY, deltaZ, &
                                side, &
                                j
integer(pInt), dimension(2,3) :: periodicImages
real(pReal)                     x, y, z, &                    ! coordinates of connection vector in neighbor lattice frame
                                xsquare, ysquare, zsquare, &  ! squares of respective coordinates
                                distance, &                   ! length of connection vector
                                segmentLength, &              ! segment length of dislocations
                                lambda, &
                                R, Rsquare, Rcube, &
                                denominator, &
                                flipSign, &
                                neighbor_ipVolumeSideLength, &
                                detFe
real(pReal), dimension(3) ::    connection, &                 ! connection vector between me and my neighbor in the deformed configuration
                                connection_neighborLattice, & ! connection vector between me and my neighbor in the lattice configuration of my neighbor
                                connection_neighborSlip, &    ! connection vector between me and my neighbor in the slip system frame of my neighbor
                                maxCoord, minCoord, &
                                meshSize, &
                                coords, &                     ! x,y,z coordinates of cell center of ip volume
                                neighbor_coords               ! x,y,z coordinates of cell center of neighbor ip volume
real(pReal), dimension(3,3) ::  sigma, &                      ! dislocation stress for one slip system in neighbor material point's slip system frame
                                Tdislo_neighborLattice, &     ! dislocation stress as 2nd Piola-Kirchhoff stress at neighbor material point
                                invFe, &                      ! inverse of my elastic deformation gradient
                                neighbor_invFe, &
                                neighborLattice2myLattice     ! mapping from neighbor MPs lattice configuration to my lattice configuration
real(pReal), dimension(2,2,maxval(totalNslip)) :: &
                                neighbor_rhoExcess            ! excess density at neighbor material point (edge/screw,mobile/dead,slipsystem)
real(pReal), dimension(2,maxval(totalNslip)) :: &
                                rhoExcessDead
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),8) :: &
                                rhoSgl                        ! single dislocation density (edge+, edge-, screw+, screw-, used edge+, used edge-, used screw+, used screw-)
logical                         inversionError

phase = material_phase(g,ip,el)
instance = phase_plasticityInstance(phase)
latticeStruct = constitutive_nonlocal_structure(instance)
ns = totalNslip(instance)



!*** get basic states

forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)
  rhoSgl(s,t) = max(state(g,ip,el)%p(iRhoU(s,t,instance)), 0.0_pReal)                               ! ensure positive single mobile densities
  rhoSgl(s,t+4_pInt) = state(g,ip,el)%p(iRhoB(s,t,instance))
endforall



!*** calculate the dislocation stress of the neighboring excess dislocation densities
!*** zero for material points of local plasticity

constitutive_nonlocal_dislocationstress = 0.0_pReal

if (.not. phase_localPlasticity(phase)) then
  call math_invert33(Fe(1:3,1:3,g,ip,el), invFe, detFe, inversionError)

  !* in case of periodic surfaces we have to find out how many periodic images in each direction we need
  
  do dir = 1_pInt,3_pInt
    maxCoord(dir) = maxval(mesh_node0(dir,:))
    minCoord(dir) = minval(mesh_node0(dir,:))
  enddo
  meshSize = maxCoord - minCoord
  coords = mesh_cellCenterCoordinates(ip,el)
  periodicImages = 0_pInt
  do dir = 1_pInt,3_pInt
    if (mesh_periodicSurface(dir)) then
      periodicImages(1,dir) =   floor((coords(dir) - cutoffRadius(instance) - minCoord(dir)) / meshSize(dir), pInt)
      periodicImages(2,dir) = ceiling((coords(dir) + cutoffRadius(instance) - maxCoord(dir)) / meshSize(dir), pInt)
    endif
  enddo

      
  !* loop through all material points (also through their periodic images if present), 
  !* but only consider nonlocal neighbors within a certain cutoff radius R
  
  do neighbor_el = 1_pInt,mesh_NcpElems
ipLoop: do neighbor_ip = 1_pInt,FE_Nips(FE_geomtype(mesh_element(2,neighbor_el)))
      neighbor_phase = material_phase(g,neighbor_ip,neighbor_el)
      if (phase_localPlasticity(neighbor_phase)) then
        cycle
      endif
      neighbor_instance = phase_plasticityInstance(neighbor_phase)
      neighbor_latticeStruct = constitutive_nonlocal_structure(neighbor_instance)
      neighbor_ns = totalNslip(neighbor_instance)
      call math_invert33(Fe(1:3,1:3,1,neighbor_ip,neighbor_el), neighbor_invFe, detFe, inversionError)
      neighbor_ipVolumeSideLength = mesh_ipVolume(neighbor_ip,neighbor_el) ** (1.0_pReal/3.0_pReal) ! reference volume used here
      forall (s = 1_pInt:neighbor_ns, c = 1_pInt:2_pInt)
        neighbor_rhoExcess(c,1,s) = state(g,neighbor_ip,neighbor_el)%p(iRhoU(s,2*c-1,neighbor_instance)) &  ! positive mobiles
                                  - state(g,neighbor_ip,neighbor_el)%p(iRhoU(s,2*c,neighbor_instance))      ! negative mobiles
        neighbor_rhoExcess(c,2,s) = abs(state(g,neighbor_ip,neighbor_el)%p(iRhoB(s,2*c-1,neighbor_instance))) & ! positive deads
                                  - abs(state(g,neighbor_ip,neighbor_el)%p(iRhoB(s,2*c,neighbor_instance)))     ! negative deads
      endforall
      Tdislo_neighborLattice = 0.0_pReal
      do deltaX = periodicImages(1,1),periodicImages(2,1)
        do deltaY = periodicImages(1,2),periodicImages(2,2)
          do deltaZ = periodicImages(1,3),periodicImages(2,3)
            
            
            !* regular case
            
            if (neighbor_el /= el .or. neighbor_ip /= ip &
                .or. deltaX /= 0_pInt .or. deltaY /= 0_pInt .or. deltaZ /= 0_pInt) then
            
              neighbor_coords = mesh_cellCenterCoordinates(neighbor_ip,neighbor_el) &
                                 + (/real(deltaX,pReal), real(deltaY,pReal), real(deltaZ,pReal)/) * meshSize
              connection = neighbor_coords - coords
              distance = sqrt(sum(connection * connection))
              if (distance > cutoffRadius(instance)) then
                cycle
              endif
                

              !* the segment length is the minimum of the third root of the control volume and the ip distance
              !* this ensures, that the central MP never sits on a neighbor dislocation segment
              
              connection_neighborLattice = math_mul33x3(neighbor_invFe, connection)
              segmentLength = min(neighbor_ipVolumeSideLength, distance)
      

              !* loop through all slip systems of the neighbor material point
              !* and add up the stress contributions from egde and screw excess on these slip systems (if significant)
      
              do s = 1_pInt,neighbor_ns
                if (all(abs(neighbor_rhoExcess(:,:,s)) < significantRho(instance))) then
                  cycle                                                                             ! not significant
                endif
                
                
                !* map the connection vector from the lattice into the slip system frame
                
                connection_neighborSlip = math_mul33x3(lattice2slip(1:3,1:3,s,neighbor_instance), &
                                                          connection_neighborLattice)
                
                
                !* edge contribution to stress
                sigma = 0.0_pReal
                
                x = connection_neighborSlip(1)
                y = connection_neighborSlip(2)
                z = connection_neighborSlip(3)
                xsquare = x * x
                ysquare = y * y
                zsquare = z * z

                do j = 1_pInt,2_pInt
                  if (abs(neighbor_rhoExcess(1,j,s)) < significantRho(instance)) then
                    cycle 
                  elseif (j > 1_pInt) then
                    x = connection_neighborSlip(1) + sign(0.5_pReal * segmentLength, &
                                                          state(g,neighbor_ip,neighbor_el)%p(iRhoB(s,1,neighbor_instance)) &
                                                        - state(g,neighbor_ip,neighbor_el)%p(iRhoB(s,2,neighbor_instance)))
                    xsquare = x * x
                  endif
                   
                  flipSign = sign(1.0_pReal, -y)
                  do side = 1_pInt,-1_pInt,-2_pInt
                    lambda = real(side,pReal) * 0.5_pReal * segmentLength - y
                    R = sqrt(xsquare + zsquare + lambda * lambda)
                    Rsquare = R * R
                    Rcube = Rsquare * R 
                    denominator = R * (R + flipSign * lambda)
                    if (denominator == 0.0_pReal) then
                      exit ipLoop
                    endif
                      
                    sigma(1,1) = sigma(1,1) - real(side,pReal) &
                                            * flipSign * z / denominator &
                                            * (1.0_pReal + xsquare / Rsquare + xsquare / denominator) &
                                            * neighbor_rhoExcess(1,j,s)
                    sigma(2,2) = sigma(2,2) - real(side,pReal) &
                                            * (flipSign * 2.0_pReal * nu(instance) * z / denominator + z * lambda / Rcube) &
                                            * neighbor_rhoExcess(1,j,s)
                    sigma(3,3) = sigma(3,3) + real(side,pReal) &
                                            * flipSign * z / denominator &
                                            * (1.0_pReal - zsquare / Rsquare - zsquare / denominator) &
                                            * neighbor_rhoExcess(1,j,s)
                    sigma(1,2) = sigma(1,2) + real(side,pReal) &
                                            * x * z / Rcube * neighbor_rhoExcess(1,j,s)
                    sigma(1,3) = sigma(1,3) + real(side,pReal) &
                                            * flipSign * x / denominator &
                                            * (1.0_pReal - zsquare / Rsquare - zsquare / denominator) &
                                            * neighbor_rhoExcess(1,j,s)
                    sigma(2,3) = sigma(2,3) - real(side,pReal) &
                                            * (nu(instance) / R - zsquare / Rcube) * neighbor_rhoExcess(1,j,s)
                  enddo
                enddo 
                
                !* screw contribution to stress
                
                x = connection_neighborSlip(1)   ! have to restore this value, because position might have been adapted for edge deads before
                do j = 1_pInt,2_pInt
                  if (abs(neighbor_rhoExcess(2,j,s)) < significantRho(instance)) then
                    cycle 
                  elseif (j > 1_pInt) then
                    y = connection_neighborSlip(2) + sign(0.5_pReal * segmentLength, &
                                                          state(g,neighbor_ip,neighbor_el)%p(iRhoB(s,3,neighbor_instance)) &
                                                        - state(g,neighbor_ip,neighbor_el)%p(iRhoB(s,4,neighbor_instance)))
                    ysquare = y * y
                  endif

                  flipSign = sign(1.0_pReal, x)
                  do side = 1_pInt,-1_pInt,-2_pInt
                    lambda = x + real(side,pReal) * 0.5_pReal * segmentLength
                    R = sqrt(ysquare + zsquare + lambda * lambda)
                    Rsquare = R * R
                    Rcube = Rsquare * R 
                    denominator = R * (R + flipSign * lambda)
                    if (denominator == 0.0_pReal) then
                      exit ipLoop
                    endif
                    
                    sigma(1,2) = sigma(1,2) - real(side,pReal) * flipSign * z * (1.0_pReal - nu(instance)) / denominator &
                                                                              * neighbor_rhoExcess(2,j,s)
                    sigma(1,3) = sigma(1,3) + real(side,pReal) * flipSign * y * (1.0_pReal - nu(instance)) / denominator &
                                                                              * neighbor_rhoExcess(2,j,s)
                  enddo
                enddo
               
                if (all(abs(sigma) < 1.0e-10_pReal)) then ! SIGMA IS NOT A REAL STRESS, THATS WHY WE NEED A REALLY SMALL VALUE HERE
                  cycle
                endif

                !* copy symmetric parts
                
                sigma(2,1) = sigma(1,2)
                sigma(3,1) = sigma(1,3)
                sigma(3,2) = sigma(2,3)

                
                !* scale stresses and map them into the neighbor material point's lattice configuration
                
                sigma = sigma * mu(neighbor_instance) * burgers(s,neighbor_instance) &
                              / (4.0_pReal * pi * (1.0_pReal - nu(neighbor_instance))) &
                              * mesh_ipVolume(neighbor_ip,neighbor_el) / segmentLength      ! reference volume is used here (according to the segment length calculation)
                Tdislo_neighborLattice = Tdislo_neighborLattice &
                      + math_mul33x33(math_transpose33(lattice2slip(1:3,1:3,s,neighbor_instance)), &
                        math_mul33x33(sigma, lattice2slip(1:3,1:3,s,neighbor_instance)))
                                            
              enddo ! slip system loop


            !* special case of central ip volume
            !* only consider dead dislocations
            !* we assume that they all sit at a distance equal to half the third root of V
            !* in direction of the according slip direction
            
            else
              
              forall (s = 1_pInt:ns, c = 1_pInt:2_pInt) &
                rhoExcessDead(c,s) = state(g,ip,el)%p(iRhoB(s,2*c-1,instance)) &  ! positive deads (here we use symmetry: if this has negative sign it is treated as negative density at positive position instead of positive density at negative position)
                                   + state(g,ip,el)%p(iRhoB(s,2*c,instance))      ! negative deads (here we use symmetry: if this has negative sign it is treated as positive density at positive position instead of negative density at negative position)

              do s = 1_pInt,ns
                if (all(abs(rhoExcessDead(:,s)) < significantRho(instance))) then
                  cycle                                                                             ! not significant
                endif
                sigma = 0.0_pReal                                                                   ! all components except for sigma13 are zero
                sigma(1,3) = - (rhoExcessDead(1,s) + rhoExcessDead(2,s) * (1.0_pReal - nu(instance))) &
                           * neighbor_ipVolumeSideLength * mu(instance) * burgers(s,instance) &
                           / (sqrt(2.0_pReal) * pi * (1.0_pReal - nu(instance)))
                sigma(3,1) = sigma(1,3)
                
                Tdislo_neighborLattice = Tdislo_neighborLattice &
                                      + math_mul33x33(math_transpose33(lattice2slip(1:3,1:3,s,instance)), &
                                                      math_mul33x33(sigma, lattice2slip(1:3,1:3,s,instance)))
                                            
              enddo ! slip system loop

            endif

          enddo ! deltaZ loop
        enddo ! deltaY loop
      enddo ! deltaX loop


      !* map the stress from the neighbor MP's lattice configuration into the deformed configuration 
      !* and back into my lattice configuration

      neighborLattice2myLattice = math_mul33x33(invFe, Fe(1:3,1:3,1,neighbor_ip,neighbor_el))
      constitutive_nonlocal_dislocationstress = constitutive_nonlocal_dislocationstress &
                                              + math_mul33x33(neighborLattice2myLattice, &
                                                math_mul33x33(Tdislo_neighborLattice, &
                                                math_transpose33(neighborLattice2myLattice)))
                        
    enddo ipLoop
  enddo ! element loop
    
endif

endfunction


!*********************************************************************
!* return array of constitutive results                              *
!*********************************************************************
function constitutive_nonlocal_postResults(Tstar_v, Fe, Temperature, dt, state, dotState, g,ip,el)

use math,     only: math_mul6x6, &
                    math_mul33x3, &
                    math_mul33x33, &
                    pi
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_ipVolume
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_plasticityInstance, &
                    phase_Noutput
use lattice,  only: lattice_Sslip_v, &
                    lattice_sd, &
                    lattice_st, &
                    lattice_sn

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
real(pReal), dimension(constitutive_nonlocal_sizePostResults(phase_plasticityInstance(material_phase(g,ip,el)))) :: &
                                            constitutive_nonlocal_postResults

!*** local variables
integer(pInt)                               myInstance, &             ! current instance of this plasticity
                                            myStructure, &            ! current lattice structure
                                            ns, &                     ! short notation for the total number of active slip systems
                                            c, &                      ! character of dislocation
                                            cs, &                     ! constitutive result index
                                            o, &                      ! index of current output
                                            t, &                      ! type of dislocation
                                            s, &                      ! index of my current slip system
                                            sLattice                  ! index of my current slip system according to lattice order
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),8) :: &
                                            rhoSgl, &                 ! current single dislocation densities (positive/negative screw and edge without dipoles)
                                            rhoDotSgl                 ! evolution rate of single dislocation densities (positive/negative screw and edge without dipoles)
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),4) :: &
                                            gdot, &                   ! shear rates
                                            v                         ! velocities
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el)))) :: &
                                            rhoForest, &              ! forest dislocation density
                                            tauThreshold, &           ! threshold shear stress
                                            tau, &                    ! current resolved shear stress
                                            tauBack                   ! back stress from pileups on same slip system
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),2) :: &
                                            rhoDip, &                 ! current dipole dislocation densities (screw and edge dipoles)
                                            rhoDotDip, &              ! evolution rate of dipole dislocation densities (screw and edge dipoles)
                                            dLower, &                 ! minimum stable dipole distance for edges and screws
                                            dUpper                    ! current maximum stable dipole distance for edges and screws
real(pReal), dimension(3,totalNslip(phase_plasticityInstance(material_phase(g,ip,el))),2) :: &
                                            m, &                      ! direction of dislocation motion for edge and screw (unit vector)
                                            m_currentconf             ! direction of dislocation motion for edge and screw (unit vector) in current configuration
real(pReal), dimension(3,totalNslip(phase_plasticityInstance(material_phase(g,ip,el)))) :: &
                                            n_currentconf             ! slip system normal (unit vector) in current configuration
real(pReal), dimension(3,3) ::              sigma

myInstance = phase_plasticityInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = totalNslip(myInstance)

cs = 0_pInt
constitutive_nonlocal_postResults = 0.0_pReal


!* short hand notations for state variables

forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)
  rhoSgl(s,t) = state(g,ip,el)%p(iRhoU(s,t,myInstance))
  rhoSgl(s,t+4_pInt) = state(g,ip,el)%p(iRhoB(s,t,myInstance))
  v(s,t) = state(g,ip,el)%p(iV(s,t,myInstance))
  rhoDotSgl(s,t) = dotState%p(iRhoU(s,t,myInstance))
  rhoDotSgl(s,t+4_pInt) = dotState%p(iRhoB(s,t,myInstance))
endforall
forall (s = 1_pInt:ns, c = 1_pInt:2_pInt)
  rhoDip(s,c) = state(g,ip,el)%p(iRhoD(s,c,myInstance))
  rhoDotDip(s,c) = dotState%p(iRhoD(s,c,myInstance))
endforall
rhoForest = state(g,ip,el)%p(iRhoF(1:ns,myInstance))
tauThreshold = state(g,ip,el)%p(iTauF(1:ns,myInstance))
tauBack = state(g,ip,el)%p(iTauB(1:ns,myInstance))



!* Calculate shear rate

forall (t = 1_pInt:4_pInt) &
  gdot(1:ns,t) = rhoSgl(1:ns,t) * burgers(1:ns,myInstance) * v(1:ns,t)
  

!* calculate limits for stable dipole height

do s = 1_pInt,ns
  sLattice = slipSystemLattice(s,myInstance)
  tau(s) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,1,sLattice,myStructure)) + tauBack(s)
  if (abs(tau(s)) < 1.0e-15_pReal) tau(s) = 1.0e-15_pReal
enddo

dLower = minDipoleHeight(1:ns,1:2,myInstance)
dUpper(1:ns,1) = mu(myInstance) * burgers(1:ns,myInstance) &
               / (8.0_pReal * pi * (1.0_pReal - nu(myInstance)) * abs(tau))
dUpper(1:ns,2) = mu(myInstance) * burgers(1:ns,myInstance) &
               / (4.0_pReal * pi * abs(tau))
forall (c = 1_pInt:2_pInt) &
  dUpper(1:ns,c) = min(1.0_pReal / sqrt(rhoSgl(1:ns,2*c-1) + rhoSgl(1:ns,2*c) & 
                       + abs(rhoSgl(1:ns,2*c+3)) + abs(rhoSgl(1:ns,2*c+4)) + rhoDip(1:ns,c)), &
                       dUpper(1:ns,c))
dUpper = max(dUpper,dLower)


!*** dislocation motion

m(1:3,1:ns,1) = lattice_sd(1:3,slipSystemLattice(1:ns,myInstance),myStructure)
m(1:3,1:ns,2) = -lattice_st(1:3,slipSystemLattice(1:ns,myInstance),myStructure)
forall (c = 1_pInt:2_pInt, s = 1_pInt:ns) &
  m_currentconf(1:3,s,c) = math_mul33x3(Fe(1:3,1:3,g,ip,el), m(1:3,s,c))
forall (s = 1_pInt:ns) &
  n_currentconf(1:3,s) = math_mul33x3(Fe(1:3,1:3,g,ip,el), &
                                      lattice_sn(1:3,slipSystemLattice(s,myInstance),myStructure))


do o = 1_pInt,phase_Noutput(material_phase(g,ip,el))
  select case(constitutive_nonlocal_output(o,myInstance))
    
    case ('rho')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl),2) + sum(rhoDip,2)
      cs = cs + ns
      
    case ('rho_sgl')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl),2)
      cs = cs + ns
      
    case ('rho_sgl_mobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl(1:ns,1:4)),2)
      cs = cs + ns
      
    case ('rho_sgl_immobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoSgl(1:ns,5:8),2)
      cs = cs + ns
      
    case ('rho_dip')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDip,2)
      cs = cs + ns
      
    case ('rho_edge')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl(1:ns,(/1,2,5,6/))),2) + rhoDip(1:ns,1)
      cs = cs + ns
      
    case ('rho_sgl_edge')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl(1:ns,(/1,2,5,6/))),2)
      cs = cs + ns
      
    case ('rho_sgl_edge_mobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoSgl(1:ns,1:2),2)
      cs = cs + ns
      
    case ('rho_sgl_edge_immobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoSgl(1:ns,5:6),2)
      cs = cs + ns
      
    case ('rho_sgl_edge_pos')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,1) + abs(rhoSgl(1:ns,5))
      cs = cs + ns
      
    case ('rho_sgl_edge_pos_mobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,1)
      cs = cs + ns
      
    case ('rho_sgl_edge_pos_immobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,5)
      cs = cs + ns
      
    case ('rho_sgl_edge_neg')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,2) + abs(rhoSgl(1:ns,6))
      cs = cs + ns
      
    case ('rho_sgl_edge_neg_mobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,2)
      cs = cs + ns
      
    case ('rho_sgl_edge_neg_immobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,6)
      cs = cs + ns
      
    case ('rho_dip_edge')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDip(1:ns,1)
      cs = cs + ns
      
    case ('rho_screw')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl(1:ns,(/3,4,7,8/))),2) + rhoDip(1:ns,2)
      cs = cs + ns
      
    case ('rho_sgl_screw')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl(1:ns,(/3,4,7,8/))),2)
      cs = cs + ns
            
    case ('rho_sgl_screw_mobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoSgl(1:ns,3:4),2)
      cs = cs + ns
      
    case ('rho_sgl_screw_immobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoSgl(1:ns,7:8),2)
      cs = cs + ns
      
    case ('rho_sgl_screw_pos')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,3) + abs(rhoSgl(1:ns,7))
      cs = cs + ns
      
    case ('rho_sgl_screw_pos_mobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,3)
      cs = cs + ns
      
    case ('rho_sgl_screw_pos_immobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,7)
      cs = cs + ns
      
    case ('rho_sgl_screw_neg')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,4) + abs(rhoSgl(1:ns,8))
      cs = cs + ns

    case ('rho_sgl_screw_neg_mobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,4)
      cs = cs + ns

    case ('rho_sgl_screw_neg_immobile')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,8)
      cs = cs + ns

    case ('rho_dip_screw')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDip(1:ns,2)
      cs = cs + ns
      
    case ('excess_rho')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = (rhoSgl(1:ns,1) + abs(rhoSgl(1:ns,5))) &
                                                         - (rhoSgl(1:ns,2) + abs(rhoSgl(1:ns,6))) &
                                                         + (rhoSgl(1:ns,3) + abs(rhoSgl(1:ns,7))) &
                                                         - (rhoSgl(1:ns,4) + abs(rhoSgl(1:ns,8)))
      cs = cs + ns
      
    case ('excess_rho_edge')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = (rhoSgl(1:ns,1) + abs(rhoSgl(1:ns,5))) &
                                                         - (rhoSgl(1:ns,2) + abs(rhoSgl(1:ns,6)))
      cs = cs + ns
      
    case ('excess_rho_screw')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = (rhoSgl(1:ns,3) + abs(rhoSgl(1:ns,7))) &
                                                         - (rhoSgl(1:ns,4) + abs(rhoSgl(1:ns,8)))
      cs = cs + ns
      
    case ('rho_forest')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoForest
      cs = cs + ns
    
    case ('delta')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = 1.0_pReal / sqrt(sum(abs(rhoSgl),2) + sum(rhoDip,2))
      cs = cs + ns
      
    case ('delta_sgl')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = 1.0_pReal / sqrt(sum(abs(rhoSgl),2))
      cs = cs + ns
      
    case ('delta_dip')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = 1.0_pReal / sqrt(sum(rhoDip,2))
      cs = cs + ns
      
    case ('shearrate')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(gdot,2)
      cs = cs + ns
      
    case ('resolvedstress')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = tau
      cs = cs + ns
      
    case ('resolvedstress_back')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = tauBack
      cs = cs + ns
      
    case ('resolvedstress_external')
      do s = 1_pInt,ns  
        sLattice = slipSystemLattice(s,myInstance)
        constitutive_nonlocal_postResults(cs+s) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,1,sLattice,myStructure))
      enddo
      cs = cs + ns
      
    case ('resistance')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = tauThreshold
      cs = cs + ns
    
    case ('rho_dot')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotSgl,2) + sum(rhoDotDip,2)
      cs = cs + ns
      
    case ('rho_dot_sgl')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotSgl,2)
      cs = cs + ns
      
    case ('rho_dot_dip')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotDip,2)
      cs = cs + ns
    
    case ('rho_dot_gen')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotMultiplicationOutput(1:ns,1,g,ip,el) &
                                                         + rhoDotMultiplicationOutput(1:ns,2,g,ip,el)
      cs = cs + ns

    case ('rho_dot_gen_edge')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotMultiplicationOutput(1:ns,1,g,ip,el)
      cs = cs + ns

    case ('rho_dot_gen_screw')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotMultiplicationOutput(1:ns,2,g,ip,el)
      cs = cs + ns
      
    case ('rho_dot_sgl2dip')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotSingle2DipoleGlideOutput(1:ns,1,g,ip,el) &
                                                         + rhoDotSingle2DipoleGlideOutput(1:ns,2,g,ip,el)
      cs = cs + ns
    
    case ('rho_dot_sgl2dip_edge')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotSingle2DipoleGlideOutput(1:ns,1,g,ip,el)
      cs = cs + ns
    
    case ('rho_dot_sgl2dip_screw')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotSingle2DipoleGlideOutput(1:ns,2,g,ip,el)
      cs = cs + ns
    
    case ('rho_dot_ann_ath')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotAthermalAnnihilationOutput(1:ns,1,g,ip,el) & 
                                                         + rhoDotAthermalAnnihilationOutput(1:ns,2,g,ip,el)
      cs = cs + ns
      
    case ('rho_dot_ann_the') 
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotThermalAnnihilationOutput(1:ns,1,g,ip,el) & 
                                                         + rhoDotThermalAnnihilationOutput(1:ns,2,g,ip,el)
      cs = cs + ns

    case ('rho_dot_ann_the_edge') 
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotThermalAnnihilationOutput(1:ns,1,g,ip,el) 
      cs = cs + ns

    case ('rho_dot_ann_the_screw') 
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotThermalAnnihilationOutput(1:ns,2,g,ip,el)
      cs = cs + ns

    case ('rho_dot_edgejogs') 
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotEdgeJogsOutput(1:ns,g,ip,el)
      cs = cs + ns

    case ('rho_dot_flux')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotFluxOutput(1:ns,1:4,g,ip,el),2) &
                                                      + sum(abs(rhoDotFluxOutput(1:ns,5:8,g,ip,el)),2)
      cs = cs + ns
    
    case ('rho_dot_flux_edge')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotFluxOutput(1:ns,1:2,g,ip,el),2) &
                                                      + sum(abs(rhoDotFluxOutput(1:ns,5:6,g,ip,el)),2)
      cs = cs + ns
      
    case ('rho_dot_flux_screw')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotFluxOutput(1:ns,3:4,g,ip,el),2) &
                                                      + sum(abs(rhoDotFluxOutput(1:ns,7:8,g,ip,el)),2)
      cs = cs + ns
            
    case ('velocity_edge_pos')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = v(1:ns,1)
      cs = cs + ns
    
    case ('velocity_edge_neg')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = v(1:ns,2)
      cs = cs + ns
    
    case ('velocity_screw_pos')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = v(1:ns,3)
      cs = cs + ns
    
    case ('velocity_screw_neg')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = v(1:ns,4)
      cs = cs + ns
    
    case ('slipdirection.x')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = m_currentconf(1,1:ns,1)
      cs = cs + ns
    
    case ('slipdirection.y')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = m_currentconf(2,1:ns,1)
      cs = cs + ns
    
    case ('slipdirection.z')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = m_currentconf(3,1:ns,1)
      cs = cs + ns
    
    case ('slipnormal.x')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = n_currentconf(1,1:ns)
      cs = cs + ns
    
    case ('slipnormal.y')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = n_currentconf(2,1:ns)
      cs = cs + ns
    
    case ('slipnormal.z')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = n_currentconf(3,1:ns)
      cs = cs + ns
    
    case ('fluxdensity_edge_pos.x')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,1) * v(1:ns,1) * m_currentconf(1,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_pos.y')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,1) * v(1:ns,1) * m_currentconf(2,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_pos.z')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,1) * v(1:ns,1) * m_currentconf(3,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_neg.x')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = - rhoSgl(1:ns,2) * v(1:ns,2) * m_currentconf(1,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_neg.y')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = - rhoSgl(1:ns,2) * v(1:ns,2) * m_currentconf(2,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_neg.z')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = - rhoSgl(1:ns,2) * v(1:ns,2) * m_currentconf(3,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_screw_pos.x')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,3) * v(1:ns,3) * m_currentconf(1,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_pos.y')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,3) * v(1:ns,3) * m_currentconf(2,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_pos.z')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,3) * v(1:ns,3) * m_currentconf(3,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_neg.x')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = - rhoSgl(1:ns,4) * v(1:ns,4) * m_currentconf(1,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_neg.y')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = - rhoSgl(1:ns,4) * v(1:ns,4) * m_currentconf(2,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_neg.z')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = - rhoSgl(1:ns,4) * v(1:ns,4) * m_currentconf(3,1:ns,2)
      cs = cs + ns
    
    case ('maximumdipoleheight_edge')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = dUpper(1:ns,1)
      cs = cs + ns
      
    case ('maximumdipoleheight_screw')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = dUpper(1:ns,2)
      cs = cs + ns
    
    case('dislocationstress')
      sigma = constitutive_nonlocal_dislocationstress(state, Fe, g, ip, el)
      constitutive_nonlocal_postResults(cs+1_pInt) = sigma(1,1)
      constitutive_nonlocal_postResults(cs+2_pInt) = sigma(2,2)
      constitutive_nonlocal_postResults(cs+3_pInt) = sigma(3,3)
      constitutive_nonlocal_postResults(cs+4_pInt) = sigma(1,2)
      constitutive_nonlocal_postResults(cs+5_pInt) = sigma(2,3)
      constitutive_nonlocal_postResults(cs+6_pInt) = sigma(3,1)
      cs = cs + 6_pInt
    
    case('accumulatedshear')
      constitutive_nonlocal_postResults(cs+1_pInt:cs+ns) = state(g,ip,el)%p(iGamma(1:ns,myInstance))
      cs = cs + ns
    
    case('boundarylayer')
      do s = 1_pInt,ns
        if (sum(abs(rhoSgl(s,1:8))) > 0.0_pReal) then
          constitutive_nonlocal_postResults(cs+s) = maxval(abs(rhoSgl(s,5:8))/(rhoSgl(s,1:4)+abs(rhoSgl(s,5:8))))
        else
          constitutive_nonlocal_postResults(cs+s) = 0.0_pReal
        endif
      enddo
      cs = cs + ns

 end select
enddo

endfunction

END MODULE
