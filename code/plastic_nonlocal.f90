!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for plasticity including dislocation flux
!--------------------------------------------------------------------------------------------------
module plastic_nonlocal
 use prec, only: &
   pReal, &
   pInt
 
 implicit none
 private
 character(len=22), dimension(11), parameter, private :: &
   BASICSTATES = ['rhoSglEdgePosMobile   ', &
                  'rhoSglEdgeNegMobile   ', &
                  'rhoSglScrewPosMobile  ', &
                  'rhoSglScrewNegMobile  ', &
                  'rhoSglEdgePosImmobile ', &
                  'rhoSglEdgeNegImmobile ', &
                  'rhoSglScrewPosImmobile', &
                  'rhoSglScrewNegImmobile', &
                  'rhoDipEdge            ', &
                  'rhoDipScrew           ', &
                  'accumulatedshear      ' ]                                                        !< list of "basic" microstructural state variables that are independent from other state variables
   
 character(len=16), dimension(3), parameter, private :: &
   DEPENDENTSTATES = ['rhoForest       ', & 
                      'tauThreshold    ', & 
                      'tauBack         ' ]                                                          !< list of microstructural state variables that depend on other state variables
   
 character(len=20), dimension(6), parameter, private :: &
 OTHERSTATES = ['velocityEdgePos     ', & 
                'velocityEdgeNeg     ', & 
                'velocityScrewPos    ', &
                'velocityScrewNeg    ', &
                'maxDipoleHeightEdge ', &
                'maxDipoleHeightScrew'  ]                                                           !< list of other dependent state variables that are not updated by microstructure
 
 real(pReal), parameter, private :: &
   KB = 1.38e-23_pReal                                                                              !< Physical parameter, Boltzmann constant in J/Kelvin

 integer(pInt), dimension(:), allocatable, public, protected :: &
   plastic_nonlocal_sizeDotState, &                                                            !< number of dotStates = number of basic state variables
   plastic_nonlocal_sizeDependentState, &                                                      !< number of dependent state variables
   plastic_nonlocal_sizeState, &                                                               !< total number of state variables
   plastic_nonlocal_sizePostResults                                                            !< cumulative size of post results

 integer(pInt), dimension(:,:), allocatable, target, public :: &
   plastic_nonlocal_sizePostResult                                                             !< size of each post result output
 
 character(len=64), dimension(:,:), allocatable, target, public :: &
   plastic_nonlocal_output                                                                     !< name of each post result output
 
 integer(pInt), dimension(:), allocatable, target, public :: &
   plastic_nonlocal_Noutput                                                                    !< number of outputs per instance of this plasticity 
 
 integer(pInt), dimension(:,:), allocatable, private :: &
   iGamma, &                                                                                        !< state indices for accumulated shear
   iRhoF, &                                                                                         !< state indices for forest density
   iTauF, &                                                                                         !< state indices for critical resolved shear stress
   iTauB                                                                                            !< state indices for backstress
 integer(pInt), dimension(:,:,:), allocatable, private :: &
   iRhoU, &                                                                                         !< state indices for unblocked density
   iRhoB, &                                                                                         !< state indices for blocked density
   iRhoD, &                                                                                         !< state indices for dipole density
   iV, &                                                                                            !< state indices for dislcation velocities
   iD                                                                                               !< state indices for stable dipole height
 
 integer(pInt), dimension(:), allocatable, public, protected :: &
   totalNslip                                                                                       !< total number of active slip systems for each instance
 
 integer(pInt), dimension(:,:), allocatable, private :: &
   Nslip, &                                                                                         !< number of active slip systems for each family and instance
   slipFamily, &                                                                                    !< lookup table relating active slip system to slip family for each instance
   slipSystemLattice, &                                                                             !< lookup table relating active slip system index to lattice slip system index for each instance
   colinearSystem                                                                                   !< colinear system to the active slip system (only valid for fcc!)
 
 real(pReal), dimension(:), allocatable, private :: &
   atomicVolume, &                                                                                  !< atomic volume
   Dsd0, &                                                                                          !< prefactor for self-diffusion coefficient
   selfDiffusionEnergy, &                                                                           !< activation enthalpy for diffusion
   aTolRho, &                                                                                       !< absolute tolerance for dislocation density in state integration
   aTolShear, &                                                                                     !< absolute tolerance for accumulated shear in state integration
   significantRho, &                                                                                !< density considered significant
   significantN, &                                                                                  !< number of dislocations considered significant
   cutoffRadius, &                                                                                  !< cutoff radius for dislocation stress
   doublekinkwidth, &                                                                               !< width of a doubkle kink in multiples of the burgers vector length b
   solidSolutionEnergy, &                                                                           !< activation energy for solid solution in J
   solidSolutionSize, &                                                                             !< solid solution obstacle size in multiples of the burgers vector length
   solidSolutionConcentration, &                                                                    !< concentration of solid solution in atomic parts
   pParam, &                                                                                        !< parameter for kinetic law (Kocks,Argon,Ashby)
   qParam, &                                                                                        !< parameter for kinetic law (Kocks,Argon,Ashby)
   viscosity, &                                                                                     !< viscosity for dislocation glide in Pa s
   fattack, &                                                                                       !< attack frequency in Hz
   rhoSglScatter, &                                                                                 !< standard deviation of scatter in initial dislocation density
   surfaceTransmissivity, &                                                                         !< transmissivity at free surface
   grainboundaryTransmissivity, &                                                                   !< transmissivity at grain boundary (identified by different texture)
   CFLfactor, &                                                                                     !< safety factor for CFL flux condition
   fEdgeMultiplication, &                                                                           !< factor that determines how much edge dislocations contribute to multiplication (0...1)
   rhoSglRandom, &
   rhoSglRandomBinning, &
   linetensionEffect, &
   edgeJogFactor
 
 real(pReal), dimension(:,:), allocatable, private :: &
   rhoSglEdgePos0, &                                                                                !< initial edge_pos dislocation density per slip system for each family and instance
   rhoSglEdgeNeg0, &                                                                                !< initial edge_neg dislocation density per slip system for each family and instance
   rhoSglScrewPos0, &                                                                               !< initial screw_pos dislocation density per slip system for each family and instance
   rhoSglScrewNeg0, &                                                                               !< initial screw_neg dislocation density per slip system for each family and instance
   rhoDipEdge0, &                                                                                   !< initial edge dipole dislocation density per slip system for each family and instance
   rhoDipScrew0, &                                                                                  !< initial screw dipole dislocation density per slip system for each family and instance
   lambda0PerSlipFamily, &                                                                          !< mean free path prefactor for each family and instance
   lambda0, &                                                                                       !< mean free path prefactor for each slip system and instance
   burgersPerSlipFamily, &                                                                          !< absolute length of burgers vector [m] for each family and instance
   burgers, &                                                                                       !< absolute length of burgers vector [m] for each slip system and instance
   interactionSlipSlip                                                                              !< coefficients for slip-slip interaction for each interaction type and instance
 
 real(pReal), dimension(:,:,:), allocatable, private :: &
   minDipoleHeightPerSlipFamily, &                                                                  !< minimum stable edge/screw dipole height for each family and instance
   minDipoleHeight, &                                                                               !< minimum stable edge/screw dipole height for each slip system and instance
   peierlsStressPerSlipFamily, &                                                                    !< Peierls stress (edge and screw) 
   peierlsStress, &                                                                                 !< Peierls stress (edge and screw) 
   forestProjectionEdge, &                                                                          !< matrix of forest projections of edge dislocations for each instance
   forestProjectionScrew, &                                                                         !< matrix of forest projections of screw dislocations for each instance
   interactionMatrixSlipSlip                                                                        !< interaction matrix of the different slip systems for each instance
 
 real(pReal), dimension(:,:,:,:), allocatable, private :: &
   lattice2slip, &                                                                                  !< orthogonal transformation matrix from lattice coordinate system to slip coordinate system (passive rotation !!!)
   rhoDotEdgeJogsOutput, &
   sourceProbability
 
 real(pReal), dimension(:,:,:,:,:), allocatable, private :: &
   rhoDotFluxOutput, & 
   rhoDotMultiplicationOutput, &
   rhoDotSingle2DipoleGlideOutput, &
   rhoDotAthermalAnnihilationOutput, &
   rhoDotThermalAnnihilationOutput, &
   nonSchmidProjection                                                                              !< combined projection of Schmid and non-Schmid contributions to the resolved shear stress (only for screws)
 
 real(pReal), dimension(:,:,:,:,:,:), allocatable, private :: &
   compatibility                                                                                    !< slip system compatibility between me and my neighbors
 
 real(pReal), dimension(:,:), allocatable, private :: &
   nonSchmidCoeff
 
 logical, dimension(:), allocatable, private :: &
   shortRangeStressCorrection, &                                                                    !< flag indicating the use of the short range stress correction by a excess density gradient term
   probabilisticMultiplication
 
 enum, bind(c) 
   enumerator :: undefined_ID, &
                 rho_ID, &
                 delta_ID, &
                 rho_edge_ID, &
                 rho_screw_ID, &
                 rho_sgl_ID, &
                 delta_sgl_ID, &
                 rho_sgl_edge_ID, &
                 rho_sgl_edge_pos_ID, &
                 rho_sgl_edge_neg_ID, &
                 rho_sgl_screw_ID, &
                 rho_sgl_screw_pos_ID, &
                 rho_sgl_screw_neg_ID, &
                 rho_sgl_mobile_ID, &
                 rho_sgl_edge_mobile_ID, &
                 rho_sgl_edge_pos_mobile_ID, &
                 rho_sgl_edge_neg_mobile_ID, &
                 rho_sgl_screw_mobile_ID, &
                 rho_sgl_screw_pos_mobile_ID, &
                 rho_sgl_screw_neg_mobile_ID, &
                 rho_sgl_immobile_ID, &
                 rho_sgl_edge_immobile_ID, &
                 rho_sgl_edge_pos_immobile_ID, &
                 rho_sgl_edge_neg_immobile_ID, &
                 rho_sgl_screw_immobile_ID, &
                 rho_sgl_screw_pos_immobile_ID, &
                 rho_sgl_screw_neg_immobile_ID, &
                 rho_dip_ID, &
                 delta_dip_ID, &
                 rho_dip_edge_ID, &
                 rho_dip_screw_ID, &
                 excess_rho_ID, &
                 excess_rho_edge_ID, &
                 excess_rho_screw_ID, &
                 rho_forest_ID, &
                 shearrate_ID, &
                 resolvedstress_ID, &
                 resolvedstress_external_ID, &
                 resolvedstress_back_ID, &
                 resistance_ID, &
                 rho_dot_ID, &
                 rho_dot_sgl_ID, &
                 rho_dot_sgl_mobile_ID, &
                 rho_dot_dip_ID, &
                 rho_dot_gen_ID, &
                 rho_dot_gen_edge_ID, &
                 rho_dot_gen_screw_ID, &
                 rho_dot_sgl2dip_ID, &
                 rho_dot_sgl2dip_edge_ID, &
                 rho_dot_sgl2dip_screw_ID, &
                 rho_dot_ann_ath_ID, &
                 rho_dot_ann_the_ID, &
                 rho_dot_ann_the_edge_ID, &
                 rho_dot_ann_the_screw_ID, &
                 rho_dot_edgejogs_ID, &
                 rho_dot_flux_ID, &
                 rho_dot_flux_mobile_ID, &
                 rho_dot_flux_edge_ID, &
                 rho_dot_flux_screw_ID, &
                 velocity_edge_pos_ID, &
                 velocity_edge_neg_ID, &
                 velocity_screw_pos_ID, &
                 velocity_screw_neg_ID, &
                 slipdirectionx_ID, &
                 slipdirectiony_ID, &
                 slipdirectionz_ID, &
                 slipnormalx_ID, &
                 slipnormaly_ID, &
                 slipnormalz_ID, &
                 fluxdensity_edge_posx_ID, &
                 fluxdensity_edge_posy_ID, &
                 fluxdensity_edge_posz_ID, &
                 fluxdensity_edge_negx_ID, &
                 fluxdensity_edge_negy_ID, &
                 fluxdensity_edge_negz_ID, &
                 fluxdensity_screw_posx_ID, &
                 fluxdensity_screw_posy_ID, &
                 fluxdensity_screw_posz_ID, &
                 fluxdensity_screw_negx_ID, &
                 fluxdensity_screw_negy_ID, &
                 fluxdensity_screw_negz_ID, &
                 maximumdipoleheight_edge_ID, &
                 maximumdipoleheight_screw_ID, &
                 accumulatedshear_ID, &
                 dislocationstress_ID
 end enum
 integer(kind(undefined_ID)), dimension(:,:),   allocatable, private :: & 
   plastic_nonlocal_outputID                                                              !< ID of each post result output
 
 public :: &
 plastic_nonlocal_init, &
 plastic_nonlocal_stateInit, &
 plastic_nonlocal_aTolState, &
 plastic_nonlocal_microstructure, &
 plastic_nonlocal_LpAndItsTangent, &
 plastic_nonlocal_dotState, &
 plastic_nonlocal_deltaState, &
 plastic_nonlocal_updateCompatibility, &
 plastic_nonlocal_postResults
 
 private :: &
 plastic_nonlocal_kinetics, &
 plastic_nonlocal_dislocationstress
 

contains

!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_init(fileUnit)
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
                    IO_warning, &
                    IO_error, &
                    IO_timeStamp, &
                    IO_EOF
use debug,    only: debug_level, &
                    debug_constitutive, &
                    debug_levelBasic
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_maxNipNeighbors
use material, only: homogenization_maxNgrains, &
                    phase_plasticity, &
                    phase_plasticityInstance, &
                    phase_Noutput, &
                    PLASTICITY_NONLOCAL_label, &
                    PLASTICITY_NONLOCAL_ID, &
                    plasticState, &
!                    material_phase, &
                    material_Nphase, &
                    MATERIAL_partPhase ,&
                    material_phase
use lattice
use numerics,only: &
   worldrank, &
  numerics_integrator


implicit none
integer(pInt), intent(in) ::                fileUnit

!*** local variables
 integer(pInt), parameter :: MAXNCHUNKS = LATTICE_maxNinteraction + 1_pInt
integer(pInt), &
    dimension(1_pInt+2_pInt*MAXNCHUNKS) ::  positions
integer(pInt)          ::                   phase, &
                                            maxNinstances, &
                                            maxTotalNslip, &
                                            f, &                ! index of my slip family
                                            instance, &                ! index of my instance of this plasticity
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
                                            Nchunks_nonSchmid = 0_pInt, &
                                            mySize = 0_pInt     ! to suppress warnings, safe as init is called only once
 character(len=65536) :: &
   tag  = '', &
   line = ''

 integer(pInt) :: sizeState, sizeDotState,sizeDependentState


 integer(pInt) :: NofMyPhase 
 
 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_NONLOCAL_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 maxNinstances = int(count(phase_plasticity == PLASTICITY_NONLOCAL_ID),pInt)
 if (maxNinstances == 0) return                                              ! we don't have to do anything if there's no instance for this constitutive law

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstances

!*** memory allocation for global variables

allocate(plastic_nonlocal_sizeDotState(maxNinstances),                          source=0_pInt)
allocate(plastic_nonlocal_sizeDependentState(maxNinstances),                    source=0_pInt)
allocate(plastic_nonlocal_sizeState(maxNinstances),                             source=0_pInt)
allocate(plastic_nonlocal_sizePostResults(maxNinstances),                       source=0_pInt)
allocate(plastic_nonlocal_sizePostResult(maxval(phase_Noutput), maxNinstances), source=0_pInt)
allocate(plastic_nonlocal_Noutput(maxNinstances),                               source=0_pInt)
allocate(plastic_nonlocal_output(maxval(phase_Noutput), maxNinstances))
         plastic_nonlocal_output = ''
allocate(plastic_nonlocal_outputID(maxval(phase_Noutput), maxNinstances),       source=undefined_ID)
allocate(Nslip(lattice_maxNslipFamily,maxNinstances),       source=0_pInt)
allocate(slipFamily(lattice_maxNslip,maxNinstances),        source=0_pInt)
allocate(slipSystemLattice(lattice_maxNslip,maxNinstances), source=0_pInt)
allocate(totalNslip(maxNinstances),                         source=0_pInt)
allocate(atomicVolume(maxNinstances),                       source=0.0_pReal)
allocate(Dsd0(maxNinstances),                               source=-1.0_pReal)
allocate(selfDiffusionEnergy(maxNinstances),                source=0.0_pReal)
allocate(aTolRho(maxNinstances),                            source=0.0_pReal)
allocate(aTolShear(maxNinstances),                          source=0.0_pReal)
allocate(significantRho(maxNinstances),                     source=0.0_pReal)
allocate(significantN(maxNinstances),                       source=0.0_pReal)
allocate(cutoffRadius(maxNinstances),                       source=-1.0_pReal)
allocate(doublekinkwidth(maxNinstances),                    source=0.0_pReal)
allocate(solidSolutionEnergy(maxNinstances),                source=0.0_pReal)
allocate(solidSolutionSize(maxNinstances),                  source=0.0_pReal)
allocate(solidSolutionConcentration(maxNinstances),         source=0.0_pReal)
allocate(pParam(maxNinstances),                             source=1.0_pReal)
allocate(qParam(maxNinstances),                             source=1.0_pReal)
allocate(viscosity(maxNinstances),                          source=0.0_pReal)
allocate(fattack(maxNinstances),                            source=0.0_pReal)
allocate(rhoSglScatter(maxNinstances),                      source=0.0_pReal)
allocate(rhoSglRandom(maxNinstances),                       source=0.0_pReal)
allocate(rhoSglRandomBinning(maxNinstances),                source=1.0_pReal)
allocate(surfaceTransmissivity(maxNinstances),              source=1.0_pReal)
allocate(grainboundaryTransmissivity(maxNinstances),        source=-1.0_pReal)
allocate(CFLfactor(maxNinstances),                          source=2.0_pReal)
allocate(fEdgeMultiplication(maxNinstances),                source=0.0_pReal)
allocate(linetensionEffect(maxNinstances),                  source=0.0_pReal)
allocate(edgeJogFactor(maxNinstances),                      source=1.0_pReal)
allocate(shortRangeStressCorrection(maxNinstances),         source=.false.)
allocate(probabilisticMultiplication(maxNinstances),        source=.false.)

allocate(rhoSglEdgePos0(lattice_maxNslipFamily,maxNinstances),                 source=-1.0_pReal)
allocate(rhoSglEdgeNeg0(lattice_maxNslipFamily,maxNinstances),                 source=-1.0_pReal)
allocate(rhoSglScrewPos0(lattice_maxNslipFamily,maxNinstances),                source=-1.0_pReal)
allocate(rhoSglScrewNeg0(lattice_maxNslipFamily,maxNinstances),                source=-1.0_pReal)
allocate(rhoDipEdge0(lattice_maxNslipFamily,maxNinstances),                    source=-1.0_pReal)
allocate(rhoDipScrew0(lattice_maxNslipFamily,maxNinstances),                   source=-1.0_pReal)
allocate(burgersPerSlipFamily(lattice_maxNslipFamily,maxNinstances),           source=0.0_pReal)
allocate(lambda0PerSlipFamily(lattice_maxNslipFamily,maxNinstances),           source=0.0_pReal)
allocate(interactionSlipSlip(lattice_maxNinteraction,maxNinstances),           source=0.0_pReal)
allocate(minDipoleHeightPerSlipFamily(lattice_maxNslipFamily,2,maxNinstances), source=-1.0_pReal)
allocate(peierlsStressPerSlipFamily(lattice_maxNslipFamily,2,maxNinstances),   source=0.0_pReal)
allocate(nonSchmidCoeff(lattice_maxNnonSchmid,maxNinstances),                  source=0.0_pReal)


 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= MATERIAL_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo
  
 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through phases of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '')  then                                                         ! stop at next part 
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit
   endif
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next phase
     phase = phase + 1_pInt                                                                         ! advance phase section counter
     if (phase_plasticity(phase) == PLASTICITY_NONLOCAL_ID) then
       Nchunks_SlipFamilies = count(lattice_NslipSystem(:,phase) > 0_pInt)
       Nchunks_SlipSlip     = maxval(lattice_InteractionSlipSlip(:,:,phase))
       Nchunks_nonSchmid    = lattice_NnonSchmid(phase)
     endif
     cycle
   endif
   if (phase > 0_pInt ) then; if (phase_plasticity(phase) == PLASTICITY_NONLOCAL_ID) then           ! one of my phases. do not short-circuit here (.and. with next if statement). It's not safe in Fortran
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my plasticity is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('rho')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('delta')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = delta_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_edge')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_edge_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_screw')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_screw_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('delta_sgl')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = delta_sgl_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_edge')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_edge_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_edge_pos')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_edge_pos_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_edge_neg')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_edge_neg_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_screw')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_screw_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_screw_pos')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_screw_pos_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_screw_neg')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_screw_neg_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_mobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_mobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_edge_mobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_edge_mobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_edge_pos_mobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_edge_pos_mobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_edge_neg_mobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_edge_neg_mobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_screw_mobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_screw_mobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_screw_pos_mobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_screw_pos_mobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_screw_neg_mobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_screw_neg_mobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_immobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_immobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_edge_immobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_edge_immobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_edge_pos_immobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_edge_pos_immobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_edge_neg_immobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_edge_neg_immobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_screw_immobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_screw_immobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_screw_pos_immobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_screw_pos_immobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_sgl_screw_neg_immobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_sgl_screw_neg_immobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dip')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dip_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('delta_dip')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = delta_dip_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dip_edge')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dip_edge_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dip_screw')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dip_screw_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('excess_rho')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = excess_rho_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('excess_rho_edge')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = excess_rho_edge_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('excess_rho_screw')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = excess_rho_screw_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_forest')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_forest_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('shearrate')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = shearrate_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('resolvedstress')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = resolvedstress_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('resolvedstress_external')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = resolvedstress_external_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('resolvedstress_back')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = resolvedstress_back_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('resistance')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = resistance_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_sgl')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_sgl_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_sgl_mobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_sgl_mobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_dip')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_dip_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_gen')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_gen_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_gen_edge')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_gen_edge_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_gen_screw')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_gen_screw_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_sgl2dip')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_sgl2dip_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_sgl2dip_edge')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_sgl2dip_edge_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_sgl2dip_screw')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_sgl2dip_screw_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_ann_ath')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_ann_ath_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_ann_the')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_ann_the_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_ann_the_edge')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_ann_the_edge_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_ann_the_screw')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_ann_the_screw_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_edgejogs')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_edgejogs_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_flux')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_flux_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_flux_mobile')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_flux_mobile_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_flux_edge')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_flux_edge_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('rho_dot_flux_screw')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = rho_dot_flux_screw_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('velocity_edge_pos')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = velocity_edge_pos_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('velocity_edge_neg')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = velocity_edge_neg_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('velocity_screw_pos')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = velocity_screw_pos_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('velocity_screw_neg')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = velocity_screw_neg_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('slipdirection.x')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = slipdirectionx_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('slipdirection.y')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = slipdirectiony_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('slipdirection.z')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = slipdirectionz_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('slipnormal.x')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = slipnormalx_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('slipnormal.y')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = slipnormaly_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('slipnormal.z')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = slipnormalz_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('fluxdensity_edge_pos.x')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = fluxdensity_edge_posx_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('fluxdensity_edge_pos.y')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = fluxdensity_edge_posy_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('fluxdensity_edge_pos.z')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = fluxdensity_edge_posz_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('fluxdensity_edge_neg.x')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = fluxdensity_edge_negx_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('fluxdensity_edge_neg.y')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = fluxdensity_edge_negy_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('fluxdensity_edge_neg.z')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = fluxdensity_edge_negz_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('fluxdensity_screw_pos.x')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = fluxdensity_screw_posx_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('fluxdensity_screw_pos.y')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = fluxdensity_screw_posy_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('fluxdensity_screw_pos.z')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = fluxdensity_screw_posz_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('fluxdensity_screw_neg.x')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = fluxdensity_screw_negx_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('fluxdensity_screw_neg.y')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = fluxdensity_screw_negy_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('fluxdensity_screw_neg.z')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = fluxdensity_screw_negz_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('maximumdipoleheight_edge')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = maximumdipoleheight_edge_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('maximumdipoleheight_screw')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = maximumdipoleheight_screw_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('accumulatedshear','accumulated_shear')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = accumulatedshear_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('dislocationstress')
             plastic_nonlocal_Noutput(instance) = plastic_nonlocal_Noutput(instance) + 1_pInt
             plastic_nonlocal_outputID(plastic_nonlocal_Noutput(instance),instance) = dislocationstress_ID
             plastic_nonlocal_output(plastic_nonlocal_Noutput(instance),instance) = &
               IO_lc(IO_stringValue(line,positions,2_pInt))
         end select
       case ('nslip')
         if (positions(1) < 1_pInt + Nchunks_SlipFamilies) &
           call IO_warning(50_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_NONLOCAL_LABEL//')')
         Nchunks_SlipFamilies = positions(1) - 1_pInt
         do f = 1_pInt, Nchunks_SlipFamilies
           Nslip(f,instance) = IO_intValue(line,positions,1_pInt+f)
         enddo
       case ('rhosgledgepos0')
         do f = 1_pInt, Nchunks_SlipFamilies
           rhoSglEdgePos0(f,instance) = IO_floatValue(line,positions,1_pInt+f)
         enddo
       case ('rhosgledgeneg0')
         do f = 1_pInt, Nchunks_SlipFamilies
           rhoSglEdgeNeg0(f,instance) = IO_floatValue(line,positions,1_pInt+f)
         enddo
       case ('rhosglscrewpos0')
         do f = 1_pInt, Nchunks_SlipFamilies
           rhoSglScrewPos0(f,instance) = IO_floatValue(line,positions,1_pInt+f)
         enddo
       case ('rhosglscrewneg0')
         do f = 1_pInt, Nchunks_SlipFamilies
           rhoSglScrewNeg0(f,instance) = IO_floatValue(line,positions,1_pInt+f)
         enddo
       case ('rhodipedge0')
         do f = 1_pInt, Nchunks_SlipFamilies
           rhoDipEdge0(f,instance) = IO_floatValue(line,positions,1_pInt+f)
         enddo
       case ('rhodipscrew0')
         do f = 1_pInt, Nchunks_SlipFamilies
           rhoDipScrew0(f,instance) = IO_floatValue(line,positions,1_pInt+f)
         enddo
       case ('lambda0')
         do f = 1_pInt, Nchunks_SlipFamilies
           lambda0PerSlipFamily(f,instance) = IO_floatValue(line,positions,1_pInt+f)
         enddo
       case ('burgers')
         do f = 1_pInt, Nchunks_SlipFamilies
           burgersPerSlipFamily(f,instance) = IO_floatValue(line,positions,1_pInt+f)
         enddo
       case('cutoffradius','r')
         cutoffRadius(instance) = IO_floatValue(line,positions,2_pInt)
       case('minimumdipoleheightedge','ddipminedge')
         do f = 1_pInt, Nchunks_SlipFamilies
           minDipoleHeightPerSlipFamily(f,1_pInt,instance) = IO_floatValue(line,positions,1_pInt+f)
         enddo
       case('minimumdipoleheightscrew','ddipminscrew')
         do f = 1_pInt, Nchunks_SlipFamilies
           minDipoleHeightPerSlipFamily(f,2_pInt,instance) = IO_floatValue(line,positions,1_pInt+f)
         enddo
       case('atomicvolume')
         atomicVolume(instance) = IO_floatValue(line,positions,2_pInt)
       case('selfdiffusionprefactor','dsd0')
         Dsd0(instance) = IO_floatValue(line,positions,2_pInt)
       case('selfdiffusionenergy','qsd')
         selfDiffusionEnergy(instance) = IO_floatValue(line,positions,2_pInt)
       case('atol_rho','atol_density','absolutetolerancedensity','absolutetolerance_density')
         aTolRho(instance) = IO_floatValue(line,positions,2_pInt)
       case('atol_shear','atol_plasticshear','atol_accumulatedshear','absolutetoleranceshear','absolutetolerance_shear')
         aTolShear(instance) = IO_floatValue(line,positions,2_pInt)
       case('significantrho','significant_rho','significantdensity','significant_density')
         significantRho(instance) = IO_floatValue(line,positions,2_pInt)
       case('significantn','significant_n','significantdislocations','significant_dislcations')
         significantN(instance) = IO_floatValue(line,positions,2_pInt)
       case ('interaction_slipslip')
          if (positions(1) < 1_pInt + Nchunks_SlipSlip) &
            call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_NONLOCAL_LABEL//')')
         do it = 1_pInt,Nchunks_SlipSlip
           interactionSlipSlip(it,instance) = IO_floatValue(line,positions,1_pInt+it)
         enddo
       case('linetension','linetensioneffect','linetension_effect')
         linetensionEffect(instance) = IO_floatValue(line,positions,2_pInt)
       case('edgejog','edgejogs','edgejogeffect','edgejog_effect')
         edgeJogFactor(instance) = IO_floatValue(line,positions,2_pInt)
       case('peierlsstressedge','peierlsstress_edge')
         do f = 1_pInt, Nchunks_SlipFamilies
           peierlsStressPerSlipFamily(f,1_pInt,instance) = IO_floatValue(line,positions,1_pInt+f)
         enddo
       case('peierlsstressscrew','peierlsstress_screw')
         do f = 1_pInt, Nchunks_SlipFamilies
           peierlsStressPerSlipFamily(f,2_pInt,instance) = IO_floatValue(line,positions,1_pInt+f)
         enddo
       case('doublekinkwidth')
         doublekinkwidth(instance) = IO_floatValue(line,positions,2_pInt)
       case('solidsolutionenergy')
         solidSolutionEnergy(instance) = IO_floatValue(line,positions,2_pInt)
       case('solidsolutionsize')
         solidSolutionSize(instance) = IO_floatValue(line,positions,2_pInt)
       case('solidsolutionconcentration')
         solidSolutionConcentration(instance) = IO_floatValue(line,positions,2_pInt)
       case('p')
         pParam(instance) = IO_floatValue(line,positions,2_pInt)
       case('q')
         qParam(instance) = IO_floatValue(line,positions,2_pInt)
       case('viscosity','glideviscosity')
         viscosity(instance) = IO_floatValue(line,positions,2_pInt)
       case('attackfrequency','fattack')
         fattack(instance) = IO_floatValue(line,positions,2_pInt)
       case('rhosglscatter')
         rhoSglScatter(instance) = IO_floatValue(line,positions,2_pInt)
       case('rhosglrandom')
         rhoSglRandom(instance) = IO_floatValue(line,positions,2_pInt)
       case('rhosglrandombinning')
         rhoSglRandomBinning(instance) = IO_floatValue(line,positions,2_pInt)
       case('surfacetransmissivity')
         surfaceTransmissivity(instance) = IO_floatValue(line,positions,2_pInt)
       case('grainboundarytransmissivity')
         grainboundaryTransmissivity(instance) = IO_floatValue(line,positions,2_pInt)
       case('cflfactor')
         CFLfactor(instance) = IO_floatValue(line,positions,2_pInt)
       case('fedgemultiplication','edgemultiplicationfactor','edgemultiplication')
         fEdgeMultiplication(instance) = IO_floatValue(line,positions,2_pInt)
       case('shortrangestresscorrection')
         shortRangeStressCorrection(instance) = IO_floatValue(line,positions,2_pInt) > 0.0_pReal
       case ('nonschmid_coefficients')
         if (positions(1) < 1_pInt + Nchunks_nonSchmid) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_NONLOCAL_label//')')
         do f = 1_pInt,Nchunks_nonSchmid
           nonSchmidCoeff(f,instance) = IO_floatValue(line,positions,1_pInt+f)
         enddo
       case('probabilisticmultiplication','randomsources','randommultiplication','discretesources')
         probabilisticMultiplication(instance) = IO_floatValue(line,positions,2_pInt) > 0.0_pReal
     end select
   endif; endif
 enddo parsingFile

 sanityChecks: do phase = 1_pInt, size(phase_plasticity)
   myPhase: if (phase_plasticity(phase) == PLASTICITY_NONLOCAL_ID) then
    instance = phase_plasticityInstance(phase) 
    if (sum(Nslip(:,instance)) <= 0_pInt) &
      call IO_error(211_pInt,ext_msg='Nslip ('//PLASTICITY_NONLOCAL_label//')')
    do o = 1_pInt,maxval(phase_Noutput)
      if(len(plastic_nonlocal_output(o,instance)) > 64_pInt) &
        call IO_error(666_pInt)
    enddo
    do f = 1_pInt,lattice_maxNslipFamily
      if (Nslip(f,instance) > 0_pInt) then
        if (rhoSglEdgePos0(f,instance) < 0.0_pReal) &
          call IO_error(211_pInt,ext_msg='rhoSglEdgePos0 ('//PLASTICITY_NONLOCAL_label//')')
        if (rhoSglEdgeNeg0(f,instance) < 0.0_pReal) &
          call IO_error(211_pInt,ext_msg='rhoSglEdgeNeg0 ('//PLASTICITY_NONLOCAL_label//')')
        if (rhoSglScrewPos0(f,instance) < 0.0_pReal) &
          call IO_error(211_pInt,ext_msg='rhoSglScrewPos0 ('//PLASTICITY_NONLOCAL_label//')')
        if (rhoSglScrewNeg0(f,instance) < 0.0_pReal) &
          call IO_error(211_pInt,ext_msg='rhoSglScrewNeg0 ('//PLASTICITY_NONLOCAL_label//')')
        if (rhoDipEdge0(f,instance) < 0.0_pReal) &
          call IO_error(211_pInt,ext_msg='rhoDipEdge0 ('//PLASTICITY_NONLOCAL_label//')')
        if (rhoDipScrew0(f,instance) < 0.0_pReal) &
          call IO_error(211_pInt,ext_msg='rhoDipScrew0 ('//PLASTICITY_NONLOCAL_label//')')
        if (burgersPerSlipFamily(f,instance) <= 0.0_pReal) &
          call IO_error(211_pInt,ext_msg='Burgers ('//PLASTICITY_NONLOCAL_label//')')
        if (lambda0PerSlipFamily(f,instance) <= 0.0_pReal) &
          call IO_error(211_pInt,ext_msg='lambda0 ('//PLASTICITY_NONLOCAL_label//')')
        if (minDipoleHeightPerSlipFamily(f,1,instance) < 0.0_pReal) &
          call IO_error(211_pInt,ext_msg='minimumDipoleHeightEdge ('//PLASTICITY_NONLOCAL_label//')')
        if (minDipoleHeightPerSlipFamily(f,2,instance) < 0.0_pReal) &
          call IO_error(211_pInt,ext_msg='minimumDipoleHeightScrew ('//PLASTICITY_NONLOCAL_label//')')
        if (peierlsStressPerSlipFamily(f,1,instance) <= 0.0_pReal) &
          call IO_error(211_pInt,ext_msg='peierlsStressEdge ('//PLASTICITY_NONLOCAL_label//')')
        if (peierlsStressPerSlipFamily(f,2,instance) <= 0.0_pReal) &
          call IO_error(211_pInt,ext_msg='peierlsStressScrew ('//PLASTICITY_NONLOCAL_label//')')
      endif
    enddo
    if (any(interactionSlipSlip(1:maxval(lattice_interactionSlipSlip(:,:,phase)),instance) < 0.0_pReal)) &
      call IO_error(211_pInt,ext_msg='interaction_SlipSlip ('//PLASTICITY_NONLOCAL_label//')')
    if (linetensionEffect(instance) < 0.0_pReal .or. linetensionEffect(instance) > 1.0_pReal) &
      call IO_error(211_pInt,ext_msg='linetension ('//PLASTICITY_NONLOCAL_label//')')
    if (edgeJogFactor(instance) < 0.0_pReal .or. edgeJogFactor(instance) > 1.0_pReal) &
      call IO_error(211_pInt,ext_msg='edgejog ('//PLASTICITY_NONLOCAL_label//')')
    if (cutoffRadius(instance) < 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='r ('//PLASTICITY_NONLOCAL_label//')')
    if (atomicVolume(instance) <= 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='atomicVolume ('//PLASTICITY_NONLOCAL_label//')')
    if (Dsd0(instance) < 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='selfDiffusionPrefactor ('//PLASTICITY_NONLOCAL_label//')')
    if (selfDiffusionEnergy(instance) <= 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='selfDiffusionEnergy ('//PLASTICITY_NONLOCAL_label//')')
    if (aTolRho(instance) <= 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='aTol_rho ('//PLASTICITY_NONLOCAL_label//')')
    if (aTolShear(instance) <= 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='aTol_shear ('//PLASTICITY_NONLOCAL_label//')')
    if (significantRho(instance) < 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='significantRho ('//PLASTICITY_NONLOCAL_label//')')
    if (significantN(instance) < 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='significantN ('//PLASTICITY_NONLOCAL_label//')')
    if (doublekinkwidth(instance) <= 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='doublekinkwidth ('//PLASTICITY_NONLOCAL_label//')')
    if (solidSolutionEnergy(instance) <= 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='solidSolutionEnergy ('//PLASTICITY_NONLOCAL_label//')')
    if (solidSolutionSize(instance) <= 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='solidSolutionSize ('//PLASTICITY_NONLOCAL_label//')')
    if (solidSolutionConcentration(instance) <= 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='solidSolutionConcentration ('//PLASTICITY_NONLOCAL_label//')')
    if (pParam(instance) <= 0.0_pReal .or. pParam(instance) > 1.0_pReal) &
      call IO_error(211_pInt,ext_msg='p ('//PLASTICITY_NONLOCAL_label//')')
    if (qParam(instance) < 1.0_pReal .or. qParam(instance) > 2.0_pReal) &
      call IO_error(211_pInt,ext_msg='q ('//PLASTICITY_NONLOCAL_label//')')
    if (viscosity(instance) <= 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='viscosity ('//PLASTICITY_NONLOCAL_label//')')
    if (fattack(instance) <= 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='attackFrequency ('//PLASTICITY_NONLOCAL_label//')')
    if (rhoSglScatter(instance) < 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='rhoSglScatter ('//PLASTICITY_NONLOCAL_label//')')
    if (rhoSglRandom(instance) < 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='rhoSglRandom ('//PLASTICITY_NONLOCAL_label//')')
    if (rhoSglRandomBinning(instance) <= 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='rhoSglRandomBinning ('//PLASTICITY_NONLOCAL_label//')')
    if (surfaceTransmissivity(instance) < 0.0_pReal .or. surfaceTransmissivity(instance) > 1.0_pReal) &
      call IO_error(211_pInt,ext_msg='surfaceTransmissivity ('//PLASTICITY_NONLOCAL_label//')')
    if (grainboundaryTransmissivity(instance) > 1.0_pReal) &
      call IO_error(211_pInt,ext_msg='grainboundaryTransmissivity ('//PLASTICITY_NONLOCAL_label//')')
    if (CFLfactor(instance) < 0.0_pReal) &
      call IO_error(211_pInt,ext_msg='CFLfactor ('//PLASTICITY_NONLOCAL_label//')')
    if (fEdgeMultiplication(instance) < 0.0_pReal .or. fEdgeMultiplication(instance) > 1.0_pReal) &
      call IO_error(211_pInt,ext_msg='edgemultiplicationfactor ('//PLASTICITY_NONLOCAL_label//')')
    
    
    !*** determine total number of active slip systems
    Nslip(1:lattice_maxNslipFamily,instance) = min(lattice_NslipSystem(1:lattice_maxNslipFamily,phase), &
                                            Nslip(1:lattice_maxNslipFamily,instance) )              ! we can't use more slip systems per family than specified in lattice
    totalNslip(instance) = sum(Nslip(1:lattice_maxNslipFamily,instance))
  endif myPhase 
enddo sanityChecks


!*** allocation of variables whose size depends on the total number of active slip systems

maxTotalNslip = maxval(totalNslip)

allocate(iRhoU(maxTotalNslip,4,maxNinstances), source=0_pInt)
allocate(iRhoB(maxTotalNslip,4,maxNinstances), source=0_pInt)
allocate(iRhoD(maxTotalNslip,2,maxNinstances), source=0_pInt)
allocate(iV(maxTotalNslip,4,maxNinstances),    source=0_pInt)
allocate(iD(maxTotalNslip,2,maxNinstances),    source=0_pInt)
allocate(iGamma(maxTotalNslip,maxNinstances),  source=0_pInt)
allocate(iRhoF(maxTotalNslip,maxNinstances),   source=0_pInt)
allocate(iTauF(maxTotalNslip,maxNinstances),   source=0_pInt)
allocate(iTauB(maxTotalNslip,maxNinstances),   source=0_pInt)
allocate(burgers(maxTotalNslip,maxNinstances),                                        source=0.0_pReal)
allocate(lambda0(maxTotalNslip,maxNinstances),                                        source=0.0_pReal)
allocate(minDipoleHeight(maxTotalNslip,2,maxNinstances),                              source=-1.0_pReal)
allocate(forestProjectionEdge(maxTotalNslip,maxTotalNslip,maxNinstances),             source=0.0_pReal)
allocate(forestProjectionScrew(maxTotalNslip,maxTotalNslip,maxNinstances),            source=0.0_pReal)
allocate(interactionMatrixSlipSlip(maxTotalNslip,maxTotalNslip,maxNinstances),        source=0.0_pReal)
allocate(lattice2slip(1:3, 1:3, maxTotalNslip,maxNinstances),                         source=0.0_pReal)
allocate(sourceProbability(maxTotalNslip,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), &
                                                                                   source=2.0_pReal)

allocate(rhoDotFluxOutput(maxTotalNslip,8,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), &
                                                                                   source=0.0_pReal)
allocate(rhoDotMultiplicationOutput(maxTotalNslip,2,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), &
                                                                                   source=0.0_pReal)
allocate(rhoDotSingle2DipoleGlideOutput(maxTotalNslip,2,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), &
                                                                                   source=0.0_pReal)
allocate(rhoDotAthermalAnnihilationOutput(maxTotalNslip,2,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), &
                                                                                   source=0.0_pReal)
allocate(rhoDotThermalAnnihilationOutput(maxTotalNslip,2,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), &
                                                                                   source=0.0_pReal)
allocate(rhoDotEdgeJogsOutput(maxTotalNslip,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), &
                                                                                   source=0.0_pReal)

allocate(compatibility(2,maxTotalNslip,maxTotalNslip,mesh_maxNipNeighbors,mesh_maxNips,mesh_NcpElems), &
                                                                                   source=0.0_pReal)
allocate(peierlsStress(maxTotalNslip,2,maxNinstances),                                source=0.0_pReal)
allocate(colinearSystem(maxTotalNslip,maxNinstances),                                 source=0_pInt)
allocate(nonSchmidProjection(3,3,4,maxTotalNslip,maxNinstances),                      source=0.0_pReal)
                                            
 initializeInstances: do phase = 1_pInt, size(phase_plasticity)
   NofMyPhase=count(material_phase==phase)
   if (phase_plasticity(phase) == PLASTICITY_NONLOCAL_ID .and. NofMyPhase/=0) then
     instance = phase_plasticityInstance(phase)
     !*** Inverse lookup of my slip system family and the slip system in lattice
     
     l = 0_pInt
     do f = 1_pInt,lattice_maxNslipFamily
       do s = 1_pInt,Nslip(f,instance)
         l = l + 1_pInt
         slipFamily(l,instance) = f
         slipSystemLattice(l,instance) = sum(lattice_NslipSystem(1:f-1_pInt, phase)) + s
     enddo; enddo
     
     
     !*** determine size of state array
     
     ns = totalNslip(instance)

     sizeDotState              = int(size(BASICSTATES),pInt) * ns
     sizeDependentState        = int(size(DEPENDENTSTATES),pInt) * ns
     sizeState                 = sizeDotState + sizeDependentState &
                               + int(size(OTHERSTATES),pInt) * ns

     !*** determine indices to state array

     l = 0_pInt
     do t = 1_pInt,4_pInt
       do s = 1_pInt,ns
         l = l + 1_pInt
         iRhoU(s,t,instance) = l
       enddo
     enddo
     do t = 1_pInt,4_pInt
       do s = 1_pInt,ns
         l = l + 1_pInt
         iRhoB(s,t,instance) = l
       enddo
     enddo
     do c = 1_pInt,2_pInt
       do s = 1_pInt,ns
         l = l + 1_pInt
         iRhoD(s,c,instance) = l
       enddo
     enddo
     do s = 1_pInt,ns
       l = l + 1_pInt
       iGamma(s,instance) = l
     enddo
     do s = 1_pInt,ns
       l = l + 1_pInt
       iRhoF(s,instance) = l
     enddo
     do s = 1_pInt,ns
       l = l + 1_pInt
       iTauF(s,instance) = l
     enddo
     do s = 1_pInt,ns
       l = l + 1_pInt
       iTauB(s,instance) = l
     enddo
     do t = 1_pInt,4_pInt
       do s = 1_pInt,ns
         l = l + 1_pInt
         iV(s,t,instance) = l
       enddo
     enddo
     do c = 1_pInt,2_pInt
       do s = 1_pInt,ns
         l = l + 1_pInt
         iD(s,c,instance) = l
       enddo
     enddo
     if (iD(ns,2,instance) /= sizeState) &  ! check if last index is equal to size of state
       call IO_error(0_pInt, ext_msg = 'state indices not properly set ('//PLASTICITY_NONLOCAL_label//')')
     
   
     !*** determine size of postResults array
     
     outputsLoop: do o = 1_pInt,plastic_nonlocal_Noutput(instance)
       select case(plastic_nonlocal_outputID(o,instance))
         case( rho_ID, &
               delta_ID, &
               rho_edge_ID, &
               rho_screw_ID, &
               rho_sgl_ID, &
               delta_sgl_ID, &
               rho_sgl_edge_ID, &
               rho_sgl_edge_pos_ID, &
               rho_sgl_edge_neg_ID, &
               rho_sgl_screw_ID, &
               rho_sgl_screw_pos_ID, &
               rho_sgl_screw_neg_ID, &
               rho_sgl_mobile_ID, &
               rho_sgl_edge_mobile_ID, &
               rho_sgl_edge_pos_mobile_ID, &
               rho_sgl_edge_neg_mobile_ID, &
               rho_sgl_screw_mobile_ID, &
               rho_sgl_screw_pos_mobile_ID, &
               rho_sgl_screw_neg_mobile_ID, &
               rho_sgl_immobile_ID, &
               rho_sgl_edge_immobile_ID, &
               rho_sgl_edge_pos_immobile_ID, &
               rho_sgl_edge_neg_immobile_ID, &
               rho_sgl_screw_immobile_ID, &
               rho_sgl_screw_pos_immobile_ID, &
               rho_sgl_screw_neg_immobile_ID, &
               rho_dip_ID, &
               delta_dip_ID, &
               rho_dip_edge_ID, &
               rho_dip_screw_ID, &
               excess_rho_ID, &
               excess_rho_edge_ID, &
               excess_rho_screw_ID, &
               rho_forest_ID, &
               shearrate_ID, &
               resolvedstress_ID, &
               resolvedstress_external_ID, &
               resolvedstress_back_ID, &
               resistance_ID, &
               rho_dot_ID, &
               rho_dot_sgl_ID, &
               rho_dot_sgl_mobile_ID, &
               rho_dot_dip_ID, &
               rho_dot_gen_ID, &
               rho_dot_gen_edge_ID, &
               rho_dot_gen_screw_ID, &
               rho_dot_sgl2dip_ID, &
               rho_dot_sgl2dip_edge_ID, &
               rho_dot_sgl2dip_screw_ID, &
               rho_dot_ann_ath_ID, &
               rho_dot_ann_the_ID, &
               rho_dot_ann_the_edge_ID, &
               rho_dot_ann_the_screw_ID, &
               rho_dot_edgejogs_ID, &
               rho_dot_flux_ID, &
               rho_dot_flux_mobile_ID, &
               rho_dot_flux_edge_ID, &
               rho_dot_flux_screw_ID, &
               velocity_edge_pos_ID, &
               velocity_edge_neg_ID, &
               velocity_screw_pos_ID, &
               velocity_screw_neg_ID, &
               slipdirectionx_ID, &
               slipdirectiony_ID, &
               slipdirectionz_ID, &
               slipnormalx_ID, &
               slipnormaly_ID, &
               slipnormalz_ID, &
               fluxdensity_edge_posx_ID, &
               fluxdensity_edge_posy_ID, &
               fluxdensity_edge_posz_ID, &
               fluxdensity_edge_negx_ID, &
               fluxdensity_edge_negy_ID, &
               fluxdensity_edge_negz_ID, &
               fluxdensity_screw_posx_ID, &
               fluxdensity_screw_posy_ID, &
               fluxdensity_screw_posz_ID, &
               fluxdensity_screw_negx_ID, &
               fluxdensity_screw_negy_ID, &
               fluxdensity_screw_negz_ID, &
               maximumdipoleheight_edge_ID, &
               maximumdipoleheight_screw_ID, &
               accumulatedshear_ID )
           mySize = totalNslip(instance)
         case(dislocationstress_ID)
           mySize = 6_pInt
         case default
       end select
   
       if (mySize > 0_pInt) then                                                                       ! any meaningful output found                               
         plastic_nonlocal_sizePostResult(o,instance) = mySize
         plastic_nonlocal_sizePostResults(instance)  = plastic_nonlocal_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop
                
     plasticState(phase)%sizeState    = sizeState
     plasticState(phase)%sizeDotState = sizeDotState
     plasticState(phase)%sizePostResults = plastic_nonlocal_sizePostResults(instance)
     plasticState(phase)%nonlocal = .true.
     plasticState(phase)%nSlip = totalNslip(instance)
     plasticState(phase)%nTwin = 0_pInt
     plasticState(phase)%nTrans= 0_pInt
     allocate(plasticState(phase)%aTolState           (sizeState),                source=0.0_pReal)
     allocate(plasticState(phase)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%state               (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%state_backup        (sizeState,NofMyPhase),     source=0.0_pReal)

     allocate(plasticState(phase)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(plasticState(phase)%deltaState          (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(plasticState(phase)%dotState_backup     (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(phase)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
       allocate(plasticState(phase)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(phase)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(phase)%RKCK45dotState    (6,sizeDotState,NofMyPhase),source=0.0_pReal)
     plasticState(phase)%slipRate => &
       plasticState(phase)%dotState(iGamma(1,instance):iGamma(ns,instance),1:NofMyPhase)
     plasticState(phase)%accumulatedSlip => &
       plasticState(phase)%state   (iGamma(1,instance):iGamma(ns,instance),1:NofMyPhase)
       
     do s1 = 1_pInt,ns 
       f = slipFamily(s1,instance)
       
       !*** burgers vector, mean free path prefactor and minimum dipole distance for each slip system
     
       burgers(s1,instance) = burgersPerSlipFamily(f,instance)
       lambda0(s1,instance) = lambda0PerSlipFamily(f,instance)
       minDipoleHeight(s1,1:2,instance) = minDipoleHeightPerSlipFamily(f,1:2,instance)
       peierlsStress(s1,1:2,instance) = peierlsStressPerSlipFamily(f,1:2,instance)
   
       do s2 = 1_pInt,ns
         
         !*** calculation of forest projections for edge and screw dislocations. s2 acts as forest for s1
   
         forestProjectionEdge(s1,s2,instance) &
             = abs(math_mul3x3(lattice_sn(1:3,slipSystemLattice(s1,instance),phase), &
                               lattice_st(1:3,slipSystemLattice(s2,instance),phase)))                   ! forest projection of edge dislocations is the projection of (t = b x n) onto the slip normal of the respective slip plane
         
         forestProjectionScrew(s1,s2,instance) &
             = abs(math_mul3x3(lattice_sn(1:3,slipSystemLattice(s1,instance),phase), &
                               lattice_sd(1:3,slipSystemLattice(s2,instance),phase)))                   ! forest projection of screw dislocations is the projection of b onto the slip normal of the respective splip plane
     
         !*** calculation of interaction matrices
   
         interactionMatrixSlipSlip(s1,s2,instance) &
             = interactionSlipSlip(lattice_interactionSlipSlip(slipSystemLattice(s1,instance), &
                                                               slipSystemLattice(s2,instance), &
                                                               phase), instance)
         
         !*** colinear slip system (only makes sense for fcc like it is defined here)
         
         if (lattice_interactionSlipSlip(slipSystemLattice(s1,instance), &
                                         slipSystemLattice(s2,instance), &
                                         phase) == 3_pInt) then
           colinearSystem(s1,instance) = s2
         endif
     
       enddo
   
       !*** rotation matrix from lattice configuration to slip system
   
       lattice2slip(1:3,1:3,s1,instance) &
           = math_transpose33( reshape([ lattice_sd(1:3, slipSystemLattice(s1,instance), phase), &
                                        -lattice_st(1:3, slipSystemLattice(s1,instance), phase), &
                                         lattice_sn(1:3, slipSystemLattice(s1,instance), phase)], [3,3]))
     enddo
   
     
     !*** combined projection of Schmid and non-Schmid contributions to the resolved shear stress (only for screws)
     !* four types t: 
     !*   1) positive screw at positive resolved stress
     !*   2) positive screw at negative resolved stress
     !*   3) negative screw at positive resolved stress
     !*   4) negative screw at negative resolved stress
     
     do s = 1_pInt,ns 
       do l = 1_pInt,lattice_NnonSchmid(phase)
         nonSchmidProjection(1:3,1:3,1,s,instance) = nonSchmidProjection(1:3,1:3,1,s,instance) &
             + nonSchmidCoeff(l,instance) * lattice_Sslip(1:3,1:3,2*l,slipSystemLattice(s,instance),phase)
         nonSchmidProjection(1:3,1:3,2,s,instance) = nonSchmidProjection(1:3,1:3,2,s,instance) &
             + nonSchmidCoeff(l,instance) * lattice_Sslip(1:3,1:3,2*l+1,slipSystemLattice(s,instance),phase)
       enddo
    nonSchmidProjection(1:3,1:3,3,s,instance) = -nonSchmidProjection(1:3,1:3,2,s,instance)
    nonSchmidProjection(1:3,1:3,4,s,instance) = -nonSchmidProjection(1:3,1:3,1,s,instance)
    forall (t = 1:4) &
      nonSchmidProjection(1:3,1:3,t,s,instance) = nonSchmidProjection(1:3,1:3,t,s,instance) &
             + lattice_Sslip(1:3,1:3,1,slipSystemLattice(s,instance),phase)
     enddo
   endif
   call plastic_nonlocal_aTolState(phase,instance)
 enddo initializeInstances

end subroutine plastic_nonlocal_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the initial microstructural state for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------

subroutine plastic_nonlocal_stateInit()
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
                    plasticState, &
                    mappingConstitutive, &
                    material_Nphase, &
                    phase_plasticity ,&
                    PLASTICITY_NONLOCAL_ID
 use numerics,only: &
                    numerics_integrator

implicit none

integer(pInt)        ::       e, &
                              i, &
                              ns, &                           ! short notation for total number of active slip systems 
                              f, &                            ! index of lattice family
                              from, &
                              upto, &
                              s, &                            ! index of slip system
                              t, &
                              j, &
                              instance, &
                              maxNinstances
real(pReal), dimension(2) ::  noise
real(pReal), dimension(4) ::  rnd   
real(pReal)                   meanDensity, &
                              totalVolume, &
                              densityBinning, &
                              minimumIpVolume

maxNinstances = int(count(phase_plasticity == PLASTICITY_NONLOCAL_ID),pInt)

do instance = 1_pInt,maxNinstances
  ns = totalNslip(instance)

  ! randomly distribute dislocation segments on random slip system and of random type in the volume 
  if (rhoSglRandom(instance) > 0.0_pReal) then

    ! get the total volume of the instance

    minimumIpVolume = huge(1.0_pReal)
    totalVolume = 0.0_pReal
    do e = 1_pInt,mesh_NcpElems
      do i = 1_pInt,FE_Nips(FE_geomtype(mesh_element(2,e)))
        if (PLASTICITY_NONLOCAL_ID == phase_plasticity(material_phase(1,i,e)) &
            .and. instance == phase_plasticityInstance(material_phase(1,i,e))) then
          totalVolume = totalVolume + mesh_ipVolume(i,e)
          minimumIpVolume = min(minimumIpVolume, mesh_ipVolume(i,e))
        endif
      enddo
    enddo
    densityBinning = rhoSglRandomBinning(instance) / minimumIpVolume ** (2.0_pReal / 3.0_pReal)

    ! subsequently fill random ips with dislocation segments until we reach the desired overall density

    meanDensity = 0.0_pReal
    do while(meanDensity < rhoSglRandom(instance))
      call random_number(rnd)
      e = nint(rnd(1)*real(mesh_NcpElems,pReal)+0.5_pReal,pInt)
      i = nint(rnd(2)*real(FE_Nips(FE_geomtype(mesh_element(2,e))),pReal)+0.5_pReal,pInt)
      if (PLASTICITY_NONLOCAL_ID == phase_plasticity(material_phase(1,i,e)) &
          .and. instance == phase_plasticityInstance(material_phase(1,i,e))) then
        s = nint(rnd(3)*real(ns,pReal)+0.5_pReal,pInt)
        t = nint(rnd(4)*4.0_pReal+0.5_pReal,pInt)
        meanDensity = meanDensity + densityBinning * mesh_ipVolume(i,e) / totalVolume
        plasticState(mappingConstitutive(2,1,i,e))%state0(iRhoU(s,t,instance),mappingConstitutive(2,1,i,e)) = &
        plasticState(mappingConstitutive(2,1,i,e))%state0(iRhoU(s,t,instance),mappingConstitutive(2,1,i,e)) &
        + densityBinning
      endif
    enddo
  ! homogeneous distribution of density with some noise
  else
    do e = 1_pInt,mesh_NcpElems
      do i = 1_pInt,FE_Nips(FE_geomtype(mesh_element(2,e)))
        if (PLASTICITY_NONLOCAL_ID == phase_plasticity(material_phase(1,i,e)) &
            .and. instance == phase_plasticityInstance(material_phase(1,i,e))) then
          do f = 1_pInt,lattice_maxNslipFamily
            from = 1_pInt + sum(Nslip(1:f-1_pInt,instance))
            upto = sum(Nslip(1:f,instance))
            do s = from,upto
              do j = 1_pInt,2_pInt
                noise(j) = math_sampleGaussVar(0.0_pReal, rhoSglScatter(instance))
              enddo
              plasticState(mappingConstitutive(2,1,i,e))%state0(iRhoU(s,1,instance),mappingConstitutive(1,1,i,e)) = &
                       rhoSglEdgePos0(f,instance) + noise(1)
              plasticState(mappingConstitutive(2,1,i,e))%state0(iRhoU(s,2,instance),mappingConstitutive(1,1,i,e)) = &
                       rhoSglEdgeNeg0(f,instance) + noise(1)
              plasticState(mappingConstitutive(2,1,i,e))%state0(iRhoU(s,3,instance),mappingConstitutive(1,1,i,e)) = &
                       rhoSglScrewPos0(f,instance) + noise(2)
              plasticState(mappingConstitutive(2,1,i,e))%state0(iRhoU(s,4,instance),mappingConstitutive(1,1,i,e)) = &
                       rhoSglScrewNeg0(f,instance) + noise(2)
            enddo
            plasticState(mappingConstitutive(2,1,i,e))%state0(iRhoD(from:upto,1,instance),mappingConstitutive(1,1,i,e)) = &
              rhoDipEdge0(f,instance)
            plasticState(mappingConstitutive(2,1,i,e))%state0(iRhoD(from:upto,2,instance),mappingConstitutive(1,1,i,e)) = &
              rhoDipScrew0(f,instance)
          enddo
        endif
      enddo
    enddo
  endif
enddo

end subroutine plastic_nonlocal_stateInit


!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_aTolState(ph,instance)
 use material, only: &
   plasticState

 implicit none
 integer(pInt), intent(in) :: & 
   instance, &                                                              !< number specifying the instance of the plasticity
   ph
 integer(pInt) :: &
   ns, &
   t, c

 ns = totalNslip(instance)
 forall (t = 1_pInt:4_pInt)
   plasticState(ph)%aTolState(iRhoU(1:ns,t,instance)) = aTolRho(instance)
   plasticState(ph)%aTolState(iRhoB(1:ns,t,instance)) = aTolRho(instance)
 end forall
 forall (c = 1_pInt:2_pInt) &
   plasticState(ph)%aTolState(iRhoD(1:ns,c,instance)) = aTolRho(instance)
 
 plasticState(ph)%aTolState(iGamma(1:ns,instance))  = aTolShear(instance)

end subroutine plastic_nonlocal_aTolState

!--------------------------------------------------------------------------------------------------
!> @brief calculates quantities characterizing the microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_microstructure(Fe, Fp, ip, el)
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
  plasticState, &
  mappingConstitutive, &
  phase_plasticityInstance
use lattice, only: &
  lattice_sd, &
  lattice_st, &
  lattice_mu, &
  lattice_nu, &
  lattice_structure, &
  LATTICE_bcc_ID, &
  LATTICE_fcc_ID

implicit none

integer(pInt), intent(in) ::    ip, &                         ! current integration point
                                el                            ! current element
real(pReal), dimension(3,3), intent(in) :: &
                                Fe, &                         ! elastic deformation gradient
                                Fp                            ! elastic deformation gradient

 integer(pInt) :: &
   ph, &                                                                                             !< phase
   of, &                                                                                             !< offset
   np, &                                                                                            !< neighbor phase
   no                                                                                               !< nieghbor offset

integer(pInt)                   neighbor_el, &                ! element number of neighboring material point
                                neighbor_ip, &                ! integration point of neighboring material point
                                instance, &                      ! my instance of this plasticity
                                neighbor_instance, &             ! instance of this plasticity of neighboring material point
                                neighbor_phase, &
                                ns, &                         ! total number of active slip systems at my material point
                                neighbor_ns, &                ! total number of active slip systems at neighboring material point
                                c, &                          ! index of dilsocation character (edge, screw)
                                s, &                          ! slip system index
                                t, &                          ! index of dilsocation type (e+, e-, s+, s-, used e+, used e-, used s+, used s-)
                                dir, &
                                n, &
                                nRealNeighbors                ! number of really existing neighbors
integer(pInt), dimension(2) ::  neighbors
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
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el)))) :: &
                                rhoForest, &                  ! forest dislocation density
                                tauBack, &                    ! back stress from pileup on same slip system
                                tauThreshold                  ! threshold shear stress
real(pReal), dimension(3,3) ::  invFe, &                      ! inverse of elastic deformation gradient
                                invFp, &                      ! inverse of plastic deformation gradient
                                connections, &
                                invConnections
real(pReal), dimension(3,mesh_maxNipNeighbors) :: &
                                connection_latticeConf
real(pReal), dimension(2,totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el)))) :: &
                                rhoExcess
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),2) :: &
                                rhoDip                        ! dipole dislocation density (edge, screw)
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),8) :: &
                                rhoSgl                        ! single dislocation density (edge+, edge-, screw+, screw-, used edge+, used edge-, used screw+, used screw-)
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))), &
                       totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el)))) :: &
                                myInteractionMatrix           ! corrected slip interaction matrix
real(pReal), dimension(2,maxval(totalNslip),mesh_maxNipNeighbors) :: &
                                neighbor_rhoExcess, &      ! excess density at neighboring material point
                                neighbor_rhoTotal          ! total density at neighboring material point
real(pReal), dimension(3,totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),2) :: &
                                m                             ! direction of dislocation motion
logical                         inversionError

ph = mappingConstitutive(2,1,ip,el)
of = mappingConstitutive(1,1,ip,el)
instance = phase_plasticityInstance(ph)
ns = totalNslip(instance)

!*** get basic states


forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)
  rhoSgl(s,t) = max(plasticState(ph)%state(iRhoU(s,t,instance),of), 0.0_pReal)                              ! ensure positive single mobile densities
  rhoSgl(s,t+4_pInt) = plasticState(ph)%state(iRhoB(s,t,instance),of)
endforall
forall (s = 1_pInt:ns, c = 1_pInt:2_pInt) &
  rhoDip(s,c) = max(plasticState(ph)%state(iRhoD(s,c,instance),of), 0.0_pReal)                              ! ensure positive dipole densities

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
!*** coefficients are corrected for the line tension effect 
!*** (see Kubin,Devincre,Hoc; 2008; Modeling dislocation storage rates and mean free paths in face-centered cubic crystals)

myInteractionMatrix = 0.0_pReal
myInteractionMatrix(1:ns,1:ns) = interactionMatrixSlipSlip(1:ns,1:ns,instance)
if (lattice_structure(ph) ==  LATTICE_bcc_ID .or. lattice_structure(ph) == LATTICE_fcc_ID) then     ! only fcc and bcc
  do s = 1_pInt,ns 
    myRhoForest = max(rhoForest(s),significantRho(instance))
    correction = (  1.0_pReal - linetensionEffect(instance) &
                  + linetensionEffect(instance) &
                  * log(0.35_pReal * burgers(s,instance) * sqrt(myRhoForest)) &
                  / log(0.35_pReal * burgers(s,instance) * 1e6_pReal)) ** 2.0_pReal
    myInteractionMatrix(s,1:ns) = correction * myInteractionMatrix(s,1:ns) 
  enddo
endif
forall (s = 1_pInt:ns) &
  tauThreshold(s) = lattice_mu(ph) * burgers(s,instance) &
                  * sqrt(dot_product((sum(abs(rhoSgl),2) + sum(abs(rhoDip),2)), myInteractionMatrix(s,1:ns)))


!*** calculate the dislocation stress of the neighboring excess dislocation densities
!*** zero for material points of local plasticity

tauBack = 0.0_pReal

if (.not. phase_localPlasticity(ph) .and. shortRangeStressCorrection(instance)) then
  call math_invert33(Fe, invFe, detFe, inversionError)
  call math_invert33(Fp, invFp, detFp, inversionError)
  rhoExcess(1,1:ns) = rhoSgl(1:ns,1) - rhoSgl(1:ns,2)
  rhoExcess(2,1:ns) = rhoSgl(1:ns,3) - rhoSgl(1:ns,4)
  FVsize = mesh_ipVolume(ip,el) ** (1.0_pReal/3.0_pReal)
  
  !* loop through my neighborhood and get the connection vectors (in lattice frame) and the excess densities
  
  nRealNeighbors = 0_pInt
  neighbor_rhoTotal = 0.0_pReal
  do n = 1_pInt,FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,el))))
    neighbor_el = mesh_ipNeighborhood(1,n,ip,el)
    neighbor_ip = mesh_ipNeighborhood(2,n,ip,el)
    np = mappingConstitutive(2,1,neighbor_ip,neighbor_el)
    no = mappingConstitutive(1,1,neighbor_ip,neighbor_el)
    if (neighbor_el > 0 .and. neighbor_ip > 0) then
      neighbor_phase = material_phase(1,neighbor_ip,neighbor_el)
      neighbor_instance = phase_plasticityInstance(neighbor_phase)
      neighbor_ns = totalNslip(neighbor_instance)
      if (.not. phase_localPlasticity(neighbor_phase) &
          .and. neighbor_instance == instance) then                                                    ! same instance should be same structure
        if (neighbor_ns == ns) then
          nRealNeighbors = nRealNeighbors + 1_pInt
          forall (s = 1_pInt:ns, c = 1_pInt:2_pInt)

            neighbor_rhoExcess(c,s,n) = &
                max(plasticState(np)%state(iRhoU(s,2*c-1,neighbor_instance),no), 0.0_pReal) &       ! positive mobiles
              - max(plasticState(np)%state(iRhoU(s,2*c,neighbor_instance),  no), 0.0_pReal)         ! negative mobiles
            neighbor_rhoTotal(c,s,n) = &
                max(plasticState(np)%state(iRhoU(s,2*c-1,neighbor_instance),no), 0.0_pReal) &       ! positive mobiles
              + max(plasticState(np)%state(iRhoU(s,2*c,neighbor_instance),  no), 0.0_pReal) &       ! negative mobiles
              + abs(plasticState(np)%state(iRhoB(s,2*c-1,neighbor_instance),no)) &                  ! positive deads
              + abs(plasticState(np)%state(iRhoB(s,2*c,neighbor_instance),  no)) &                  ! negative deads
              + max(plasticState(np)%state(iRhoD(s,c,neighbor_instance),    no), 0.0_pReal)         ! dipoles

          endforall
          connection_latticeConf(1:3,n) = &
            math_mul33x3(invFe, mesh_ipCoordinates(1:3,neighbor_ip,neighbor_el) &
                                - mesh_ipCoordinates(1:3,ip,el))
          normal_latticeConf = math_mul33x3(math_transpose33(invFp), mesh_ipAreaNormal(1:3,n,ip,el))
          if (math_mul3x3(normal_latticeConf,connection_latticeConf(1:3,n)) < 0.0_pReal) then       ! neighboring connection points in opposite direction to face normal: must be periodic image
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
        neighbor_rhoExcess(1:2,1:ns,n) = rhoExcess
      endif
    else
      ! free surface -> use central values instead
      connection_latticeConf(1:3,n) = 0.0_pReal
      neighbor_rhoExcess(1:2,1:ns,n) = rhoExcess
    endif
  enddo
  

  !* loop through the slip systems and calculate the dislocation gradient by
  !* 1. interpolation of the excess density in the neighorhood
  !* 2. interpolation of the dead dislocation density in the central volume
  
  m(1:3,1:ns,1) =  lattice_sd(1:3,slipSystemLattice(1:ns,instance),ph)
  m(1:3,1:ns,2) = -lattice_st(1:3,slipSystemLattice(1:ns,instance),ph)

  do s = 1_pInt,ns
    
    !* gradient from interpolation of neighboring excess density

    do c = 1_pInt,2_pInt
      do dir = 1_pInt,3_pInt
        neighbors(1) = 2_pInt * dir - 1_pInt
        neighbors(2) = 2_pInt * dir
        connections(dir,1:3) = connection_latticeConf(1:3,neighbors(1)) &
                             - connection_latticeConf(1:3,neighbors(2))
        rhoExcessDifferences(dir) = neighbor_rhoExcess(c,s,neighbors(1)) &
                                  - neighbor_rhoExcess(c,s,neighbors(2))
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
      rhoTotal(c) = (sum(abs(rhoSgl(s,[2*c-1,2*c,2*c+3,2*c+4]))) + rhoDip(s,c) &
                            + sum(neighbor_rhoTotal(c,s,:)))  / real(1_pInt + nRealNeighbors,pReal)
    forall (c = 1_pInt:2_pInt, rhoTotal(c) > 0.0_pReal) &
      rhoExcessGradient_over_rho(c) = rhoExcessGradient(c) / rhoTotal(c)
    
    !* gives the local stress correction when multiplied with a factor

    tauBack(s) = - lattice_mu(ph) * burgers(s,instance) / (2.0_pReal * pi) &
               * (rhoExcessGradient_over_rho(1) / (1.0_pReal - lattice_nu(ph)) &
                  + rhoExcessGradient_over_rho(2))

  enddo
endif


!*** set dependent states
plasticState(ph)%state(iRhoF(1:ns,instance),of) = rhoForest
plasticState(ph)%state(iTauF(1:ns,instance),of) = tauThreshold
plasticState(ph)%state(iTauB(1:ns,instance),of) = tauBack

#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
    write(6,'(/,a,i8,1x,i2,1x,i1,/)') '<< CONST >> nonlocal_microstructure at el ip ',el,ip
    write(6,'(a,/,12x,12(e10.3,1x))') '<< CONST >> rhoForest', rhoForest
    write(6,'(a,/,12x,12(f10.5,1x))') '<< CONST >> tauThreshold / MPa', tauThreshold/1e6
    write(6,'(a,/,12x,12(f10.5,1x),/)') '<< CONST >> tauBack / MPa', tauBack/1e6
  endif
#endif

end subroutine plastic_nonlocal_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief calculates kinetics 
!--------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_kinetics(v, dv_dtau, dv_dtauNS, tau, tauNS, &
                                          tauThreshold, c, Temperature, ip, el)

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
integer(pInt), intent(in) ::                ip, &                       !< current integration point
                                            el, &                       !< current element number
                                            c                           !< dislocation character (1:edge, 2:screw)
real(pReal), intent(in) ::                  Temperature                 !< temperature
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el)))), &
             intent(in) ::                  tau, &                      !< resolved external shear stress (without non Schmid effects)
                                            tauNS, &                    !< resolved external shear stress (including non Schmid effects)
                                            tauThreshold                !< threshold shear stress

!*** output variables
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el)))), &
                            intent(out) ::  v, &                        !< velocity
                                            dv_dtau, &                  !< velocity derivative with respect to resolved shear stress (without non Schmid contributions)
                                            dv_dtauNS                   !< velocity derivative with respect to resolved shear stress (including non Schmid contributions)

!*** local variables
integer(pInt)    ::                         instance, &                    !< current instance of this plasticity
                                            ns, &                       !< short notation for the total number of active slip systems
                                            s                           !< index of my current slip system
real(pReal)                                 tauRel_P, & 
                                            tauRel_S, &
                                            tauEff, &                   !< effective shear stress
                                            tPeierls, &                 !< waiting time in front of a peierls barriers
                                            tSolidSolution, &           !< waiting time in front of a solid solution obstacle
                                            vViscous, &                 !< viscous glide velocity
                                            dtPeierls_dtau, &           !< derivative with respect to resolved shear stress
                                            dtSolidSolution_dtau, &     !< derivative with respect to resolved shear stress
                                            meanfreepath_S, &           !< mean free travel distance for dislocations between two solid solution obstacles
                                            meanfreepath_P, &           !< mean free travel distance for dislocations between two Peierls barriers
                                            jumpWidth_P, &              !< depth of activated area
                                            jumpWidth_S, &              !< depth of activated area
                                            activationLength_P, &       !< length of activated dislocation line
                                            activationLength_S, &       !< length of activated dislocation line
                                            activationVolume_P, &       !< volume that needs to be activated to overcome barrier
                                            activationVolume_S, &       !< volume that needs to be activated to overcome barrier
                                            activationEnergy_P, &       !< energy that is needed to overcome barrier
                                            activationEnergy_S, &       !< energy that is needed to overcome barrier
                                            criticalStress_P, &         !< maximum obstacle strength
                                            criticalStress_S, &         !< maximum obstacle strength
                                            mobility                    !< dislocation mobility


instance = phase_plasticityInstance(material_phase(1_pInt,ip,el))
ns = totalNslip(instance)

v = 0.0_pReal
dv_dtau = 0.0_pReal
dv_dtauNS = 0.0_pReal


if (Temperature > 0.0_pReal) then
  do s = 1_pInt,ns
    if (abs(tau(s)) > tauThreshold(s)) then

      !* Peierls contribution
      !* Effective stress includes non Schmid constributions
      !* The derivative only gives absolute values; the correct sign is taken care of in the formula for the derivative of the velocity
      
      tauEff = max(0.0_pReal, abs(tauNS(s)) - tauThreshold(s))                                      ! ensure that the effective stress is positive
      meanfreepath_P = burgers(s,instance)
      jumpWidth_P = burgers(s,instance)
      activationLength_P = doublekinkwidth(instance) * burgers(s,instance)
      activationVolume_P = activationLength_P * jumpWidth_P * burgers(s,instance)
      criticalStress_P = peierlsStress(s,c,instance)
      activationEnergy_P = criticalStress_P * activationVolume_P
      tauRel_P = min(1.0_pReal, tauEff / criticalStress_P)                                          ! ensure that the activation probability cannot become greater than one
      tPeierls = 1.0_pReal / fattack(instance) &
               * exp(activationEnergy_P / (KB * Temperature) &
                     * (1.0_pReal - tauRel_P**pParam(instance))**qParam(instance))
      if (tauEff < criticalStress_P) then
        dtPeierls_dtau = tPeierls * pParam(instance) * qParam(instance) * activationVolume_P / (KB * Temperature) &
                       * (1.0_pReal - tauRel_P**pParam(instance))**(qParam(instance)-1.0_pReal) &
                                    * tauRel_P**(pParam(instance)-1.0_pReal) 
      else
        dtPeierls_dtau = 0.0_pReal
      endif


      !* Contribution from solid solution strengthening
      !* The derivative only gives absolute values; the correct sign is taken care of in the formula for the derivative of the velocity

      tauEff = abs(tau(s)) - tauThreshold(s)
      meanfreepath_S = burgers(s,instance) / sqrt(solidSolutionConcentration(instance))
      jumpWidth_S = solidSolutionSize(instance) * burgers(s,instance)
      activationLength_S = burgers(s,instance) / sqrt(solidSolutionConcentration(instance))
      activationVolume_S = activationLength_S * jumpWidth_S * burgers(s,instance)
      activationEnergy_S = solidSolutionEnergy(instance)
      criticalStress_S = activationEnergy_S / activationVolume_S
      tauRel_S = min(1.0_pReal, tauEff / criticalStress_S)                                          ! ensure that the activation probability cannot become greater than one
      tSolidSolution = 1.0_pReal / fattack(instance) &
                     * exp(activationEnergy_S / (KB * Temperature) &
                           * (1.0_pReal - tauRel_S**pParam(instance))**qParam(instance))
      if (tauEff < criticalStress_S) then
        dtSolidSolution_dtau = tSolidSolution * pParam(instance) * qParam(instance) &
                             * activationVolume_S / (KB * Temperature) &
                             * (1.0_pReal - tauRel_S**pParam(instance))**(qParam(instance)-1.0_pReal) &
                                            * tauRel_S**(pParam(instance)-1.0_pReal) 
      else
        dtSolidSolution_dtau = 0.0_pReal
      endif


      !* viscous glide velocity
      
      tauEff = abs(tau(s)) - tauThreshold(s)
      mobility = burgers(s,instance) / viscosity(instance)
      vViscous = mobility * tauEff


      !* Mean velocity results from waiting time at peierls barriers and solid solution obstacles with respective meanfreepath of 
      !* free flight at glide velocity in between.       
      !* adopt sign from resolved stress

      v(s) = sign(1.0_pReal,tau(s)) &
           / (tPeierls / meanfreepath_P + tSolidSolution / meanfreepath_S + 1.0_pReal / vViscous)
      dv_dtau(s) = v(s) * v(s) * (dtSolidSolution_dtau / meanfreepath_S &
                                 + mobility / (vViscous * vViscous))
      dv_dtauNS(s) = v(s) * v(s) * dtPeierls_dtau / meanfreepath_P 
    endif
  enddo
endif
    

#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
    write(6,'(/,a,i8,1x,i2,1x,i1,/)') '<< CONST >> nonlocal_kinetics at el ip',el,ip
    write(6,'(a,/,12x,12(f12.5,1x))') '<< CONST >> tauThreshold / MPa', tauThreshold / 1e6_pReal
    write(6,'(a,/,12x,12(f12.5,1x))') '<< CONST >> tau / MPa', tau / 1e6_pReal
    write(6,'(a,/,12x,12(f12.5,1x))') '<< CONST >> tauNS / MPa', tauNS / 1e6_pReal
    write(6,'(a,/,12x,12(f12.5,1x))') '<< CONST >> v / 1e-3m/s', v * 1e3
    write(6,'(a,/,12x,12(e12.5,1x))') '<< CONST >> dv_dtau', dv_dtau
    write(6,'(a,/,12x,12(e12.5,1x))') '<< CONST >> dv_dtauNS', dv_dtauNS
  endif
#endif

end subroutine plastic_nonlocal_kinetics

!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_LpAndItsTangent(Lp, dLp_dTstar99, Tstar_v, Temperature, ipc, ip, el)

use math,     only: math_Plain3333to99, &
                    math_mul6x6, &
                    math_mul33xx33, &
                    math_Mandel6to33
use debug,    only: debug_level, &
                    debug_constitutive, &
                    debug_levelBasic, &
                    debug_levelExtensive, &
                    debug_levelSelective, &
                    debug_g, &
                    debug_i, &
                    debug_e
use material, only: material_phase, &
                    plasticState, &
                    mappingConstitutive,&
                    phase_plasticityInstance
use lattice,  only: lattice_Sslip, &
                    lattice_Sslip_v, &
                    lattice_NnonSchmid
use mesh,     only: mesh_ipVolume

implicit none

!*** input variables
integer(pInt), intent(in) ::                ipc, &
                                            ip, &                                                   !< current integration point
                                            el                                                      !< current element number
real(pReal), intent(in) ::                  Temperature                                             !< temperature
real(pReal), dimension(6), intent(in) ::    Tstar_v                                                 !< 2nd Piola-Kirchhoff stress in Mandel notation


!*** output variables
real(pReal), dimension(3,3), intent(out) :: Lp                                                      !< plastic velocity gradient
real(pReal), dimension(9,9), intent(out) :: dLp_dTstar99                                            !< derivative of Lp with respect to Tstar (9x9 matrix)

!*** local variables
integer(pInt)                               instance, &                                             !< current instance of this plasticity
                                            ns, &                                                   !< short notation for the total number of active slip systems
                                            i, &
                                            j, &
                                            k, &
                                            l, &
                                            ph, &                                                    !phase number
                                            of, &                                                    !offset
                                            t, &                                                    !< dislocation type
                                            s, &                                                    !< index of my current slip system
                                            sLattice                                                !< index of my current slip system according to lattice order
real(pReal), dimension(3,3,3,3) ::          dLp_dTstar3333                                          !< derivative of Lp with respect to Tstar (3x3x3x3 matrix)
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),8) :: &
                                            rhoSgl                                                  !< single dislocation densities (including blocked) 
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),4) :: &
                                            v, &                                                    !< velocity
                                            tauNS, &                                                !< resolved shear stress including non Schmid and backstress terms
                                            dv_dtau, &                                              !< velocity derivative with respect to the shear stress
                                            dv_dtauNS                                               !< velocity derivative with respect to the shear stress
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el)))) :: &
                                            tau, &                                                  !< resolved shear stress including backstress terms
                                            gdotTotal, &                                            !< shear rate
                                            tauBack, &                                              !< back stress from dislocation gradients on same slip system
                                            tauThreshold                                            !< threshold shear stress
!*** shortcut for mapping
ph = mappingConstitutive(2,1_pInt,ip,el)
of = mappingConstitutive(1,1_pInt,ip,el)

!*** initialize local variables

Lp = 0.0_pReal
dLp_dTstar3333 = 0.0_pReal

instance = phase_plasticityInstance(ph)
ns = totalNslip(instance)


!*** shortcut to state variables 


forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)

  rhoSgl(s,t) = max(plasticState(ph)%state(iRhoU(s,t,instance),of), 0.0_pReal)                         ! ensure positive single mobile densities
  rhoSgl(s,t+4_pInt) = plasticState(ph)%state(iRhoB(s,t,instance),of)
endforall
where (abs(rhoSgl) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(instance) &
  .or. abs(rhoSgl) < significantRho(instance)) &
  rhoSgl = 0.0_pReal

tauBack = plasticState(ph)%state(iTauB(1:ns,instance),of)
tauThreshold = plasticState(ph)%state(iTauF(1:ns,instance),of)


!*** get resolved shear stress
!*** for screws possible non-schmid contributions are also taken into account

do s = 1_pInt,ns
  sLattice = slipSystemLattice(s,instance)
  tau(s) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,1,sLattice,ph))
  tauNS(s,1) = tau(s)
  tauNS(s,2) = tau(s)
  if (tau(s) > 0.0_pReal) then
    tauNS(s,3) = math_mul33xx33(math_Mandel6to33(Tstar_v), nonSchmidProjection(1:3,1:3,1,s,instance))
    tauNS(s,4) = math_mul33xx33(math_Mandel6to33(Tstar_v), nonSchmidProjection(1:3,1:3,3,s,instance))
  else
    tauNS(s,3) = math_mul33xx33(math_Mandel6to33(Tstar_v), nonSchmidProjection(1:3,1:3,2,s,instance))
    tauNS(s,4) = math_mul33xx33(math_Mandel6to33(Tstar_v), nonSchmidProjection(1:3,1:3,4,s,instance))
  endif
enddo
forall (t = 1_pInt:4_pInt) &
  tauNS(1:ns,t) = tauNS(1:ns,t) + tauBack                                                           ! add backstress
tau = tau + tauBack                                                                                 ! add backstress


!*** get dislocation velocity and its tangent and store the velocity in the state array

! edges 
call plastic_nonlocal_kinetics(v(1:ns,1), dv_dtau(1:ns,1), dv_dtauNS(1:ns,1), &
                                    tau(1:ns), tauNS(1:ns,1), tauThreshold(1:ns), &
                                    1_pInt, Temperature, ip, el)
v(1:ns,2) = v(1:ns,1)
dv_dtau(1:ns,2) = dv_dtau(1:ns,1)
dv_dtauNS(1:ns,2) = dv_dtauNS(1:ns,1)

!screws
if (lattice_NnonSchmid(ph) == 0_pInt) then                                                 ! no non-Schmid contributions
  forall(t = 3_pInt:4_pInt)
    v(1:ns,t) = v(1:ns,1)
    dv_dtau(1:ns,t) = dv_dtau(1:ns,1)
    dv_dtauNS(1:ns,t) = dv_dtauNS(1:ns,1)
  endforall
else                                                                                                ! take non-Schmid contributions into account
  do t = 3_pInt,4_pInt
    call plastic_nonlocal_kinetics(v(1:ns,t), dv_dtau(1:ns,t), dv_dtauNS(1:ns,t), &
                                        tau(1:ns), tauNS(1:ns,t), tauThreshold(1:ns), &
                                        2_pInt , Temperature, ip, el)
  enddo
endif


!*** store velocity in state

forall (t = 1_pInt:4_pInt) &
  plasticState(ph)%state(iV(1:ns,t,instance),of) = v(1:ns,t)
!*** Bauschinger effect

forall (s = 1_pInt:ns, t = 5_pInt:8_pInt, rhoSgl(s,t) * v(s,t-4_pInt) < 0.0_pReal) &
  rhoSgl(s,t-4_pInt) = rhoSgl(s,t-4_pInt) + abs(rhoSgl(s,t))


!*** Calculation of Lp and its tangent

gdotTotal = sum(rhoSgl(1:ns,1:4) * v, 2) * burgers(1:ns,instance)

do s = 1_pInt,ns
  sLattice = slipSystemLattice(s,instance)  
  Lp = Lp + gdotTotal(s) * lattice_Sslip(1:3,1:3,1,sLattice,ph)

  ! Schmid contributions to tangent
  forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt,k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
    dLp_dTstar3333(i,j,k,l) = dLp_dTstar3333(i,j,k,l) &
        + lattice_Sslip(i,j,1,sLattice,ph) * lattice_Sslip(k,l,1,sLattice,ph) &
        * sum(rhoSgl(s,1:4) * dv_dtau(s,1:4)) * burgers(s,instance)

  ! non Schmid contributions to tangent
  if (tau(s) > 0.0_pReal) then
    forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt,k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
      dLp_dTstar3333(i,j,k,l) = dLp_dTstar3333(i,j,k,l) &
          + lattice_Sslip(i,j,1,sLattice,ph) &
          * ( nonSchmidProjection(k,l,1,s,instance) * rhoSgl(s,3) * dv_dtauNS(s,3) & 
            + nonSchmidProjection(k,l,3,s,instance) * rhoSgl(s,4) * dv_dtauNS(s,4) ) &
          * burgers(s,instance)
  else
    forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt,k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
      dLp_dTstar3333(i,j,k,l) = dLp_dTstar3333(i,j,k,l) &
          + lattice_Sslip(i,j,1,sLattice,ph) &
          * ( nonSchmidProjection(k,l,2,s,instance) * rhoSgl(s,3) * dv_dtauNS(s,3) & 
            + nonSchmidProjection(k,l,4,s,instance) * rhoSgl(s,4) * dv_dtauNS(s,4) ) &
          * burgers(s,instance)
  endif
enddo
dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)


#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt )) then
    write(6,'(/,a,i8,1x,i2,1x,i1,/)') '<< CONST >> nonlocal_LpandItsTangent at el ip',el,ip
    write(6,'(a,/,12x,12(f12.5,1x))') '<< CONST >> gdot total / 1e-3',gdotTotal*1e3_pReal
    write(6,'(a,/,3(12x,3(f12.7,1x),/))') '<< CONST >> Lp',transpose(Lp)
  endif
#endif

end subroutine plastic_nonlocal_LpAndItsTangent



!--------------------------------------------------------------------------------------------------
!> @brief (instantaneous) incremental change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_deltaState(Tstar_v,ip,el)
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
use lattice,  only: lattice_Sslip_v ,&
                    lattice_mu, &
                    lattice_nu
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_ipVolume
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    plasticState, &
                    mappingConstitutive, &
                    phase_plasticityInstance

implicit none
integer(pInt), intent(in) ::                ip, &                    ! current grain number
                                            el                        ! current element number
real(pReal), dimension(6), intent(in) ::    Tstar_v                   ! current 2nd Piola-Kirchhoff stress in Mandel notation


 integer(pInt) :: &
   ph, &                                                                   !< phase
   of                                                                                        !< offset

integer(pInt) ::instance, &               ! current instance of this plasticity
                                            ns, &                     ! short notation for the total number of active slip systems
                                            c, &                      ! character of dislocation
                                            t, &                      ! type of dislocation
                                            s, &                      ! index of my current slip system
                                            sLattice                  ! index of my current slip system according to lattice order
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1,ip,el))),10) :: &
                                            deltaRho, &                     ! density increment
                                            deltaRhoRemobilization, &       ! density increment by remobilization
                                            deltaRhoDipole2SingleStress     ! density increment by dipole dissociation (by stress change)
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1,ip,el))),8) :: &
                                            rhoSgl                        ! current single dislocation densities (positive/negative screw and edge without dipoles)
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1,ip,el))),4) :: &
                                            v                             ! dislocation glide velocity
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1,ip,el)))) :: &
                                            tau, &                        ! current resolved shear stress
                                            tauBack                       ! current back stress from pileups on same slip system
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1,ip,el))),2) :: &
                                            rhoDip, &                     ! current dipole dislocation densities (screw and edge dipoles)
                                            dLower, &                     ! minimum stable dipole distance for edges and screws
                                            dUpper, &                     ! current maximum stable dipole distance for edges and screws
                                            dUpperOld, &                  ! old maximum stable dipole distance for edges and screws
                                            deltaDUpper                   ! change in maximum stable dipole distance for edges and screws

#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) &
    write(6,'(/,a,i8,1x,i2,1x,i1,/)') '<< CONST >> nonlocal_deltaState at el ip ',el,ip
#endif

 ph = mappingConstitutive(2,1,ip,el)
 of = mappingConstitutive(1,1,ip,el)
 instance = phase_plasticityInstance(ph)
 ns = totalNslip(instance)


!*** shortcut to state variables 

 forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)
  rhoSgl(s,t) = max(plasticState(ph)%state(iRhoU(s,t,instance),of), 0.0_pReal)                              ! ensure positive single mobile densities
  rhoSgl(s,t+4_pInt) = plasticState(ph)%state(iRhoB(s,t,instance),of)
  v(s,t) = plasticState(ph)%state(iV(s,t,instance),of)
endforall
forall (s = 1_pInt:ns, c = 1_pInt:2_pInt)
  rhoDip(s,c) = max(plasticState(ph)%state(iRhoD(s,c,instance),of), 0.0_pReal)                              ! ensure positive dipole densities
  dUpperOld(s,c) = plasticState(ph)%state(iD(s,c,instance),of)
endforall
  tauBack = plasticState(ph)%state(iTauB(1:ns,instance),of)

where (abs(rhoSgl) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(instance) &
  .or. abs(rhoSgl) < significantRho(instance)) &
  rhoSgl = 0.0_pReal
where (abs(rhoDip) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(instance) &
  .or. abs(rhoDip) < significantRho(instance)) &
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
  sLattice = slipSystemLattice(s,instance)  
  tau(s) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,1,sLattice,ph)) + tauBack(s)
  if (abs(tau(s)) < 1.0e-15_pReal) tau(s) = 1.0e-15_pReal
enddo
dLower = minDipoleHeight(1:ns,1:2,instance)
dUpper(1:ns,1) = lattice_mu(ph) * burgers(1:ns,instance) &
               / (8.0_pReal * pi * (1.0_pReal - lattice_nu(ph)) * abs(tau))
dUpper(1:ns,2) = lattice_mu(ph) * burgers(1:ns,instance) / (4.0_pReal * pi * abs(tau))


forall (c = 1_pInt:2_pInt)
  where(sqrt(rhoSgl(1:ns,2*c-1)+rhoSgl(1:ns,2*c)+&
        abs(rhoSgl(1:ns,2*c+3))+abs(rhoSgl(1:ns,2*c+4))+rhoDip(1:ns,c)) >= tiny(0.0_pReal)) &
    dUpper(1:ns,c) = min(1.0_pReal / sqrt(rhoSgl(1:ns,2*c-1) + rhoSgl(1:ns,2*c) & 
                       + abs(rhoSgl(1:ns,2*c+3)) + abs(rhoSgl(1:ns,2*c+4)) + rhoDip(1:ns,c)), &
                       dUpper(1:ns,c))
end forall
dUpper = max(dUpper,dLower)
deltaDUpper = dUpper - dUpperOld


!*** dissociation by stress increase
deltaRhoDipole2SingleStress = 0.0_pReal
forall (c=1_pInt:2_pInt, s=1_pInt:ns, deltaDUpper(s,c) < 0.0_pReal .and. &
                                           (dUpperOld(s,c) - dLower(s,c))/= 0.0_pReal) &
  deltaRhoDipole2SingleStress(s,8_pInt+c) = rhoDip(s,c) * deltaDUpper(s,c) &
                                           / (dUpperOld(s,c) - dLower(s,c))

forall (t=1_pInt:4_pInt) &
  deltaRhoDipole2SingleStress(1_pInt:ns,t) = -0.5_pReal &
                                  * deltaRhoDipole2SingleStress(1_pInt:ns,(t-1_pInt)/2_pInt+9_pInt)


!*** store new maximum dipole height in state

forall (s = 1_pInt:ns, c = 1_pInt:2_pInt) &
  plasticState(ph)%state(iD(s,c,instance),of) = dUpper(s,c) 



!****************************************************************************
!*** assign the changes in the dislocation densities to deltaState

deltaRho = deltaRhoRemobilization &
         + deltaRhoDipole2SingleStress
plasticState(ph)%deltaState(:,of) = 0.0_pReal
forall (s = 1:ns, t = 1_pInt:4_pInt)
  plasticState(ph)%deltaState(iRhoU(s,t,instance),of)= deltaRho(s,t)
  plasticState(ph)%deltaState(iRhoB(s,t,instance),of) = deltaRho(s,t+4_pInt)
endforall
forall (s = 1:ns, c = 1_pInt:2_pInt) &
  plasticState(ph)%deltaState(iRhoD(s,c,instance),of) = deltaRho(s,c+8_pInt)


#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt )) then
    write(6,'(a,/,8(12x,12(e12.5,1x),/))') '<< CONST >> dislocation remobilization', deltaRhoRemobilization(1:ns,1:8)
    write(6,'(a,/,10(12x,12(e12.5,1x),/),/)') '<< CONST >> dipole dissociation by stress increase', deltaRhoDipole2SingleStress
  endif
#endif

end subroutine plastic_nonlocal_deltaState

!---------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!---------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_dotState(Tstar_v, Fe, Fp, Temperature, &
                                     timestep,subfrac, ip,el)

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
                    plasticState, &
                    mappingConstitutive, &
                    phase_plasticity ,&
                    PLASTICITY_NONLOCAL_ID
use lattice,  only: lattice_Sslip_v, &
                    lattice_sd, &
                    lattice_st ,&
                    lattice_mu, &
                    lattice_nu, &
                    lattice_structure, &
                    LATTICE_bcc_ID, & 
                    LATTICE_fcc_ID

implicit none

!*** input variables
integer(pInt), intent(in) ::                ip, &                                                   !< current integration point
                                            el                                                      !< current element number
real(pReal), intent(in) ::                  Temperature, &                                          !< temperature
                                            timestep                                                !< substepped crystallite time increment
real(pReal), dimension(6), intent(in) ::    Tstar_v                                                 !< current 2nd Piola-Kirchhoff stress in Mandel notation
real(pReal), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            subfrac                                                 !< fraction of timestep at the beginning of the substepped crystallite time increment 
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            Fe, &                                                   !< elastic deformation gradient
                                            Fp                                                      !< plastic deformation gradient

 
!*** local variables
integer(pInt) ::                            ph, &                              
                                            instance, &                                             !< current instance of this plasticity
                                            neighbor_instance, &                                    !< instance of my neighbor's plasticity
                                            ns, &                                                   !< short notation for the total number of active slip systems
                                            c, &                                                    !< character of dislocation
                                            n, &                                                    !< index of my current neighbor
                                            neighbor_el, &                                          !< element number of my neighbor
                                            neighbor_ip, &                                          !< integration point of my neighbor
                                            neighbor_n, &                                           !< neighbor index pointing to me when looking from my neighbor
                                            opposite_neighbor, &                                    !< index of my opposite neighbor
                                            opposite_ip, &                                          !< ip of my opposite neighbor
                                            opposite_el, &                                          !< element index of my opposite neighbor
                                            opposite_n, &                                           !< neighbor index pointing to me when looking from my opposite neighbor
                                            t, &                                                    !< type of dislocation
                                            o,&                                                     !< offset shortcut
                                            no,&                                                    !< neighbour offset shortcut
                                            p,&                                                     !< phase shortcut
                                            np,&                                                    !< neighbour phase shortcut
                                            topp, &                                                 !< type of dislocation with opposite sign to t
                                            s, &                                                    !< index of my current slip system
                                            sLattice                                                !< index of my current slip system according to lattice order
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),10) :: &
                                            rhoDot, &                                               !< density evolution
                                            rhoDotMultiplication, &                                 !< density evolution by multiplication
                                            rhoDotFlux, &                                           !< density evolution by flux
                                            rhoDotSingle2DipoleGlide, &                             !< density evolution by dipole formation (by glide)
                                            rhoDotAthermalAnnihilation, &                           !< density evolution by athermal annihilation
                                            rhoDotThermalAnnihilation                               !< density evolution by thermal annihilation
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),8) :: &
                                            rhoSgl, &                                               !< current single dislocation densities (positive/negative screw and edge without dipoles)
                                            rhoSglOriginal, &
                                            neighbor_rhoSgl, &                                      !< current single dislocation densities of neighboring ip (positive/negative screw and edge without dipoles)
                                            rhoSgl0, &                                              !< single dislocation densities at start of cryst inc (positive/negative screw and edge without dipoles)
                                            my_rhoSgl                                               !< single dislocation densities of central ip (positive/negative screw and edge without dipoles)
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),4) :: &
                                            v, &                                                    !< current dislocation glide velocity
                                            v0, &                                                   !< dislocation glide velocity at start of cryst inc
                                            my_v, &                                                 !< dislocation glide velocity of central ip
                                            neighbor_v, &                                           !< dislocation glide velocity of enighboring ip
                                            gdot                                                    !< shear rates
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el)))) :: &
                                            rhoForest, &                                            !< forest dislocation density
                                            tauThreshold, &                                         !< threshold shear stress
                                            tau, &                                                  !< current resolved shear stress
                                            tauBack, &                                              !< current back stress from pileups on same slip system
                                            vClimb, &                                               !< climb velocity of edge dipoles
                                            nSources
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),2) :: &
                                            rhoDip, &                                               !< current dipole dislocation densities (screw and edge dipoles)
                                            rhoDipOriginal, & 
                                            dLower, &                                               !< minimum stable dipole distance for edges and screws
                                            dUpper                                                  !< current maximum stable dipole distance for edges and screws
real(pReal), dimension(3,totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),4) :: &
                                            m                                                       !< direction of dislocation motion
real(pReal), dimension(3,3) ::              my_F, &                                                 !< my total deformation gradient
                                            neighbor_F, &                                           !< total deformation gradient of my neighbor
                                            my_Fe, &                                                !< my elastic deformation gradient
                                            neighbor_Fe, &                                          !< elastic deformation gradient of my neighbor
                                            Favg                                                    !< average total deformation gradient of me and my neighbor
real(pReal), dimension(3) ::                normal_neighbor2me, &                                   !< interface normal pointing from my neighbor to me in neighbor's lattice configuration
                                            normal_neighbor2me_defConf, &                           !< interface normal pointing from my neighbor to me in shared deformed configuration
                                            normal_me2neighbor, &                                   !< interface normal pointing from me to my neighbor in my lattice configuration
                                            normal_me2neighbor_defConf                              !< interface normal pointing from me to my neighbor in shared deformed configuration
real(pReal)                                 area, &                                                 !< area of the current interface
                                            transmissivity, &                                       !< overall transmissivity of dislocation flux to neighboring material point
                                            lineLength, &                                           !< dislocation line length leaving the current interface
                                            selfDiffusion, &                                        !< self diffusion
                                            rnd, & 
                                            meshlength
logical                                     considerEnteringFlux, &
                                            considerLeavingFlux
                                            

 p = mappingConstitutive(2,1,ip,el)
 o = mappingConstitutive(1,1,ip,el)



#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) &
    write(6,'(/,a,i8,1x,i2,1x,i1,/)') '<< CONST >> nonlocal_dotState at el ip ',el,ip
#endif

ph = material_phase(1_pInt,ip,el)
instance = phase_plasticityInstance(ph)
ns = totalNslip(instance)

tau = 0.0_pReal
gdot = 0.0_pReal


!*** shortcut to state variables 


forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)
  rhoSgl(s,t) =     max(plasticState(p)%state(iRhoU(s,t,instance),o), 0.0_pReal)                      ! ensure positive single mobile densities
  rhoSgl(s,t+4_pInt) =  plasticState(p)%state(iRhoB(s,t,instance),o)
  v(s,t) =              plasticState(p)%state(iV   (s,t,instance),o)
endforall
forall (s = 1_pInt:ns, c = 1_pInt:2_pInt)
  rhoDip(s,c) =     max(plasticState(p)%state(iRhoD(s,c,instance),o), 0.0_pReal)                      ! ensure positive dipole densities
endforall
rhoForest =              plasticState(p)%state(iRhoF(1:ns,instance),o)
tauThreshold =           plasticState(p)%state(iTauF(1:ns,instance),o)
tauBack =                plasticState(p)%state(iTauB(1:ns,instance),o)

rhoSglOriginal = rhoSgl
rhoDipOriginal = rhoDip
where (abs(rhoSgl) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(instance) &
  .or. abs(rhoSgl) < significantRho(instance)) &
  rhoSgl = 0.0_pReal
where (abs(rhoDip) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(instance) &
  .or. abs(rhoDip) < significantRho(instance)) &
  rhoDip = 0.0_pReal

if (numerics_timeSyncing) then
  forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)
    rhoSgl0(s,t) =    max(plasticState(p)%state0(iRhoU(s,t,instance),o), 0.0_pReal)
    rhoSgl0(s,t+4_pInt) = plasticState(p)%state0(iRhoB(s,t,instance),o)
    v0(s,t) =             plasticState(p)%state0(iV   (s,t,instance),o)
  endforall
  where (abs(rhoSgl0) * mesh_ipVolume(ip,el) ** 0.667_pReal < significantN(instance) &
    .or. abs(rhoSgl0) < significantRho(instance)) &
    rhoSgl0 = 0.0_pReal
endif
  


!*** sanity check for timestep

if (timestep <= 0.0_pReal) then                                                                      ! if illegal timestep... Why here and not on function entry??
  plasticState(p)%dotState = 0.0_pReal                                                         ! ...return without doing anything (-> zero dotState)
  return
endif



!****************************************************************************
!*** Calculate shear rate

forall (t = 1_pInt:4_pInt) &
  gdot(1_pInt:ns,t) = rhoSgl(1_pInt:ns,t) * burgers(1:ns,instance) * v(1:ns,t)

#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt )) then
    write(6,'(a,/,10(12x,12(e12.5,1x),/))') '<< CONST >> rho / 1/m^2', rhoSgl, rhoDip
    write(6,'(a,/,4(12x,12(e12.5,1x),/))') '<< CONST >> gdot / 1/s',gdot
  endif
#endif



!****************************************************************************
!*** calculate limits for stable dipole height

do s = 1_pInt,ns   ! loop over slip systems
  sLattice = slipSystemLattice(s,instance)  
  tau(s) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,1,sLattice,ph)) + tauBack(s)
  if (abs(tau(s)) < 1.0e-15_pReal) tau(s) = 1.0e-15_pReal
enddo

dLower = minDipoleHeight(1:ns,1:2,instance)
dUpper(1:ns,1) = lattice_mu(ph) * burgers(1:ns,instance) &
               / (8.0_pReal * pi * (1.0_pReal - lattice_nu(ph)) * abs(tau))
dUpper(1:ns,2) = lattice_mu(ph) * burgers(1:ns,instance) &
               / (4.0_pReal * pi * abs(tau))
forall (c = 1_pInt:2_pInt)
  where(sqrt(rhoSgl(1:ns,2*c-1)+rhoSgl(1:ns,2*c)+&
        abs(rhoSgl(1:ns,2*c+3))+abs(rhoSgl(1:ns,2*c+4))+rhoDip(1:ns,c)) >= tiny(0.0_pReal)) &
    dUpper(1:ns,c) = min(1.0_pReal / sqrt(rhoSgl(1:ns,2*c-1) + rhoSgl(1:ns,2*c) & 
                       + abs(rhoSgl(1:ns,2*c+3)) + abs(rhoSgl(1:ns,2*c+4)) + rhoDip(1:ns,c)), &
                       dUpper(1:ns,c))
end forall
dUpper = max(dUpper,dLower)

!****************************************************************************
!*** calculate dislocation multiplication

rhoDotMultiplication = 0.0_pReal
if (lattice_structure(ph) == LATTICE_bcc_ID) then                                                 ! BCC
  forall (s = 1:ns, sum(abs(v(s,1:4))) > 0.0_pReal)
    rhoDotMultiplication(s,1:2) = sum(abs(gdot(s,3:4))) / burgers(s,instance) &                        ! assuming double-cross-slip of screws to be decisive for multiplication
                                * sqrt(rhoForest(s)) / lambda0(s,instance) ! &                         ! mean free path
                                ! * 2.0_pReal * sum(abs(v(s,3:4))) / sum(abs(v(s,1:4)))             ! ratio of screw to overall velocity determines edge generation
    rhoDotMultiplication(s,3:4) = sum(abs(gdot(s,3:4))) / burgers(s,instance) &                        ! assuming double-cross-slip of screws to be decisive for multiplication
                                * sqrt(rhoForest(s)) / lambda0(s,instance) ! &                         ! mean free path
                                ! * 2.0_pReal * sum(abs(v(s,1:2))) / sum(abs(v(s,1:4)))             ! ratio of edge to overall velocity determines screw generation
  endforall

else                                                                                                ! ALL OTHER STRUCTURES
  if (probabilisticMultiplication(instance)) then
    meshlength = mesh_ipVolume(ip,el)**0.333_pReal
    where(sum(rhoSgl(1:ns,1:4),2) > 0.0_pReal)
      nSources = (sum(rhoSgl(1:ns,1:2),2) * fEdgeMultiplication(instance) + sum(rhoSgl(1:ns,3:4),2)) &
               / sum(rhoSgl(1:ns,1:4),2) * meshlength / lambda0(1:ns,instance)*sqrt(rhoForest(1:ns))
    elsewhere
      nSources = meshlength / lambda0(1:ns,instance) * sqrt(rhoForest(1:ns))
    endwhere
    do s = 1_pInt,ns
      if (nSources(s) < 1.0_pReal) then
        if (sourceProbability(s,1_pInt,ip,el) > 1.0_pReal) then
          call random_number(rnd)
          sourceProbability(s,1_pInt,ip,el) = rnd
          !$OMP FLUSH(sourceProbability)
        endif
        if (sourceProbability(s,1_pInt,ip,el) > 1.0_pReal - nSources(s)) then
          rhoDotMultiplication(s,1:4) = sum(rhoSglOriginal(s,1:4) * abs(v(s,1:4))) / meshlength
        endif
      else
        sourceProbability(s,1_pInt,ip,el) = 2.0_pReal
        rhoDotMultiplication(s,1:4) = &
          (sum(abs(gdot(s,1:2))) * fEdgeMultiplication(instance) + sum(abs(gdot(s,3:4)))) &
          / burgers(s,instance) * sqrt(rhoForest(s)) / lambda0(s,instance)
      endif
    enddo
#ifndef _OPENMP
    if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt &
        .and. ((debug_e == el .and. debug_i == ip)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt )) &
      write(6,'(a,/,4(12x,12(f12.5,1x),/,/))') '<< CONST >> sources', nSources
#endif
  else
    rhoDotMultiplication(1:ns,1:4) = spread( &
        (sum(abs(gdot(1:ns,1:2)),2) * fEdgeMultiplication(instance) + sum(abs(gdot(1:ns,3:4)),2)) &
      * sqrt(rhoForest(1:ns)) / lambda0(1:ns,instance) / burgers(1:ns,instance), 2, 4)
  endif
endif



!****************************************************************************
!*** calculate dislocation fluxes (only for nonlocal plasticity)

rhoDotFlux = 0.0_pReal
!? why needed here
if (.not. phase_localPlasticity(material_phase(1_pInt,ip,el))) then                                    ! only for nonlocal plasticity

  !*** check CFL (Courant-Friedrichs-Lewy) condition for flux

  if (any( abs(gdot) > 0.0_pReal &                                                                  ! any active slip system ...
          .and. CFLfactor(instance) * abs(v) * timestep &
              > mesh_ipVolume(ip,el) / maxval(mesh_ipArea(:,ip,el)))) then                          ! ...with velocity above critical value (we use the reference volume and area for simplicity here)
#ifndef _OPENMP
    if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt) then 
      write(6,'(a,i5,a,i2)') '<< CONST >> CFL condition not fullfilled at el ',el,' ip ',ip
      write(6,'(a,e10.3,a,e10.3)') '<< CONST >> velocity is at  ', &
        maxval(abs(v), abs(gdot) > 0.0_pReal &
                       .and. CFLfactor(instance) * abs(v) * timestep &
                             > mesh_ipVolume(ip,el) / maxval(mesh_ipArea(:,ip,el))), &
        ' at a timestep of ',timestep
      write(6,'(a)') '<< CONST >> enforcing cutback !!!'
    endif
#endif
    plasticState(p)%dotState =  DAMASK_NaN                                                     ! -> return NaN and, hence, enforce cutback
    return
  endif


  !*** be aware of the definition of lattice_st = lattice_sd x lattice_sn !!!
  !*** opposite sign to our p vector in the (s,p,n) triplet !!!
  
  m(1:3,1:ns,1) =  lattice_sd(1:3, slipSystemLattice(1:ns,instance), ph)
  m(1:3,1:ns,2) = -lattice_sd(1:3, slipSystemLattice(1:ns,instance), ph)
  m(1:3,1:ns,3) = -lattice_st(1:3, slipSystemLattice(1:ns,instance), ph)
  m(1:3,1:ns,4) =  lattice_st(1:3, slipSystemLattice(1:ns,instance), ph)
  
  my_Fe = Fe(1:3,1:3,1_pInt,ip,el)
  my_F = math_mul33x33(my_Fe, Fp(1:3,1:3,1_pInt,ip,el))
  
  do n = 1_pInt,FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,el))))                       ! loop through my neighbors
!       write(6,*) 'c'
    neighbor_el = mesh_ipNeighborhood(1,n,ip,el)
    neighbor_ip = mesh_ipNeighborhood(2,n,ip,el)
    neighbor_n  = mesh_ipNeighborhood(3,n,ip,el)
    np = mappingConstitutive(2,1,neighbor_ip,neighbor_el)
    no = mappingConstitutive(1,1,neighbor_ip,neighbor_el)

    opposite_neighbor = n + mod(n,2_pInt) - mod(n+1_pInt,2_pInt)
    opposite_el = mesh_ipNeighborhood(1,opposite_neighbor,ip,el)
    opposite_ip = mesh_ipNeighborhood(2,opposite_neighbor,ip,el)
    opposite_n  = mesh_ipNeighborhood(3,opposite_neighbor,ip,el)
  
    if (neighbor_n > 0_pInt) then                                                                    ! if neighbor exists, average deformation gradient
      neighbor_instance = phase_plasticityInstance(material_phase(1_pInt,neighbor_ip,neighbor_el))
      neighbor_Fe = Fe(1:3,1:3,1_pInt,neighbor_ip,neighbor_el)
      neighbor_F = math_mul33x33(neighbor_Fe, Fp(1:3,1:3,1_pInt,neighbor_ip,neighbor_el))
      Favg = 0.5_pReal * (my_F + neighbor_F)
    else                                                                                            ! if no neighbor, take my value as average
      Favg = my_F
    endif
    

    !* FLUX FROM MY NEIGHBOR TO ME
    !* This is only considered, if I have a neighbor of nonlocal plasticity 
    !* (also nonlocal constitutive law with local properties) that is at least a little bit 
    !* compatible.
    !* If it's not at all compatible, no flux is arriving, because everything is dammed in front of
    !* my neighbor's interface.
    !* The entering flux from my neighbor will be distributed on my slip systems according to the 
    !*compatibility
    
    considerEnteringFlux = .false.
    neighbor_v = 0.0_pReal        ! needed for check of sign change in flux density below 
    neighbor_rhoSgl = 0.0_pReal
    if (neighbor_n > 0_pInt) then
      if (phase_plasticity(material_phase(1,neighbor_ip,neighbor_el)) == PLASTICITY_NONLOCAL_ID &
          .and. any(compatibility(:,:,:,n,ip,el) > 0.0_pReal)) &
        considerEnteringFlux = .true.
    endif
    
    if (considerEnteringFlux) then
      if(numerics_timeSyncing .and. (subfrac(1_pInt,neighbor_ip,neighbor_el) /= subfrac(1_pInt,ip,el))) &
                                                                                              then  ! for timesyncing: in case of a timestep at the interface we have to use "state0" to make sure that fluxes n both sides are equal
        forall (s = 1:ns, t = 1_pInt:4_pInt)

          neighbor_v(s,t)  =         plasticState(np)%state0(iV   (s,t,neighbor_instance),no)
          neighbor_rhoSgl(s,t) = max(plasticState(np)%state0(iRhoU(s,t,neighbor_instance),no),0.0_pReal)

        endforall
      else
        forall (s = 1:ns, t = 1_pInt:4_pInt)
          neighbor_v(s,t) =          plasticState(np)%state(iV   (s,t,neighbor_instance),no)
          neighbor_rhoSgl(s,t) = max(plasticState(np)%state(iRhoU(s,t,neighbor_instance),no), &
                                                                                          0.0_pReal)
        endforall
      endif

      where (neighbor_rhoSgl * mesh_ipVolume(neighbor_ip,neighbor_el) ** 0.667_pReal < significantN(instance) &
        .or. neighbor_rhoSgl < significantRho(instance)) &
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
              .and. v(s,t) * neighbor_v(s,t) >= 0.0_pReal ) then                                    ! ... only if no sign change in flux density  
            lineLength = neighbor_rhoSgl(s,t) * neighbor_v(s,t) &
                       * math_mul3x3(m(1:3,s,t), normal_neighbor2me) * area                         ! positive line length that wants to enter through this interface
            where (compatibility(c,1_pInt:ns,s,n,ip,el) > 0.0_pReal) &                              ! positive compatibility...
              rhoDotFlux(1_pInt:ns,t) = rhoDotFlux(1_pInt:ns,t) &
                                      + lineLength / mesh_ipVolume(ip,el) &                         ! ... transferring to equally signed mobile dislocation type
                                      * compatibility(c,1_pInt:ns,s,n,ip,el) ** 2.0_pReal
            where (compatibility(c,1_pInt:ns,s,n,ip,el) < 0.0_pReal) &                              ! ..negative compatibility...
              rhoDotFlux(1_pInt:ns,topp) = rhoDotFlux(1_pInt:ns,topp) &
                                         + lineLength / mesh_ipVolume(ip,el) &                      ! ... transferring to opposite signed mobile dislocation type
                                         * compatibility(c,1_pInt:ns,s,n,ip,el) ** 2.0_pReal
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
      if (phase_plasticity(material_phase(1,opposite_ip,opposite_el)) /= PLASTICITY_NONLOCAL_ID) &
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
        if (subfrac(1_pInt,ip,el) == 0.0_pReal) then
          my_rhoSgl = rhoSgl0
          my_v = v0
        elseif (neighbor_n > 0_pInt) then
          if (subfrac(1_pInt,neighbor_ip,neighbor_el) == 0.0_pReal) then
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
            if (my_v(s,t) * neighbor_v(s,t) >= 0.0_pReal) then                                      ! no sign change in flux density
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
  rhoDotSingle2DipoleGlide(1:ns,2*c-1) = -2.0_pReal * dUpper(1:ns,c) / burgers(1:ns,instance) &
                                                    * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) &         ! negative mobile --> positive mobile
                                                       + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1)) &       ! positive mobile --> negative mobile
                                                       + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1))) ! positive mobile --> negative immobile

  rhoDotSingle2DipoleGlide(1:ns,2*c) = -2.0_pReal * dUpper(1:ns,c) / burgers(1:ns,instance) &
                                                  * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) &           ! negative mobile --> positive mobile
                                                     + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1)) &         ! positive mobile --> negative mobile
                                                     + abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c)))     ! negative mobile --> positive immobile

  rhoDotSingle2DipoleGlide(1:ns,2*c+3) = -2.0_pReal * dUpper(1:ns,c) / burgers(1:ns,instance) &
                                                    * rhoSgl(1:ns,2*c+3) * abs(gdot(1:ns,2*c))            ! negative mobile --> positive immobile

  rhoDotSingle2DipoleGlide(1:ns,2*c+4) = -2.0_pReal * dUpper(1:ns,c) / burgers(1:ns,instance) &
                                                    * rhoSgl(1:ns,2*c+4) * abs(gdot(1:ns,2*c-1))          ! positive mobile --> negative immobile

  rhoDotSingle2DipoleGlide(1:ns,c+8) = - rhoDotSingle2DipoleGlide(1:ns,2*c-1) &
                                       - rhoDotSingle2DipoleGlide(1:ns,2*c) &
                                       + abs(rhoDotSingle2DipoleGlide(1:ns,2*c+3)) &
                                       + abs(rhoDotSingle2DipoleGlide(1:ns,2*c+4))
enddo


!*** athermal annihilation

rhoDotAthermalAnnihilation = 0.0_pReal

forall (c=1_pInt:2_pInt) &  
  rhoDotAthermalAnnihilation(1:ns,c+8_pInt) = -2.0_pReal * dLower(1:ns,c) / burgers(1:ns,instance) &
               * (  2.0_pReal * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1))) &             ! was single hitting single
                  + 2.0_pReal * (abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c)) + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1))) & ! was single hitting immobile single or was immobile single hit by single
                  + rhoDip(1:ns,c) * (abs(gdot(1:ns,2*c-1)) + abs(gdot(1:ns,2*c))))                                                 ! single knocks dipole constituent
! annihilated screw dipoles leave edge jogs behind on the colinear system
if (lattice_structure(ph) == LATTICE_fcc_ID) & ! only fcc
  forall (s = 1:ns, colinearSystem(s,instance) > 0_pInt) &
    rhoDotAthermalAnnihilation(colinearSystem(s,instance),1:2) = - rhoDotAthermalAnnihilation(s,10) &
      * 0.25_pReal * sqrt(rhoForest(s)) * (dUpper(s,2) + dLower(s,2)) * edgeJogFactor(instance)

  
  
!*** thermally activated annihilation of edge dipoles by climb

rhoDotThermalAnnihilation = 0.0_pReal
selfDiffusion = Dsd0(instance) * exp(-selfDiffusionEnergy(instance) / (KB * Temperature))
vClimb =  atomicVolume(instance) * selfDiffusion / ( KB * Temperature ) &
          * lattice_mu(ph) / ( 2.0_pReal * PI * (1.0_pReal-lattice_nu(ph)) ) &
          * 2.0_pReal / ( dUpper(1:ns,1) + dLower(1:ns,1) )
forall (s = 1_pInt:ns, dUpper(s,1) > dLower(s,1)) &
  rhoDotThermalAnnihilation(s,9) = max(- 4.0_pReal * rhoDip(s,1) * vClimb(s) / (dUpper(s,1) - dLower(s,1)), &
                                       - rhoDip(s,1) / timestep - rhoDotAthermalAnnihilation(s,9) &
                                                                - rhoDotSingle2DipoleGlide(s,9))    ! make sure that we do not annihilate more dipoles than we have



!****************************************************************************
!*** assign the rates of dislocation densities to my dotState
!*** if evolution rates lead to negative densities, a cutback is enforced

rhoDot = 0.0_pReal
rhoDot = rhoDotFlux &
       + rhoDotMultiplication &
       + rhoDotSingle2DipoleGlide &
       + rhoDotAthermalAnnihilation &
       + rhoDotThermalAnnihilation 

if (numerics_integrationMode == 1_pInt) then                                                        ! save rates for output if in central integration mode
  rhoDotFluxOutput(1:ns,1:8,1_pInt,ip,el) = rhoDotFlux(1:ns,1:8)
  rhoDotMultiplicationOutput(1:ns,1:2,1_pInt,ip,el) = rhoDotMultiplication(1:ns,[1,3])
  rhoDotSingle2DipoleGlideOutput(1:ns,1:2,1_pInt,ip,el) = rhoDotSingle2DipoleGlide(1:ns,9:10)
  rhoDotAthermalAnnihilationOutput(1:ns,1:2,1_pInt,ip,el) = rhoDotAthermalAnnihilation(1:ns,9:10)
  rhoDotThermalAnnihilationOutput(1:ns,1:2,1_pInt,ip,el) = rhoDotThermalAnnihilation(1:ns,9:10)
  rhoDotEdgeJogsOutput(1:ns,1_pInt,ip,el) = 2.0_pReal * rhoDotThermalAnnihilation(1:ns,1)
endif


#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt &
      .and. ((debug_e == el .and. debug_i == ip .and. debug_g == 1_pInt)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt )) then
    write(6,'(a,/,4(12x,12(e12.5,1x),/))')  '<< CONST >> dislocation multiplication', &
                                            rhoDotMultiplication(1:ns,1:4) * timestep
    write(6,'(a,/,8(12x,12(e12.5,1x),/))')  '<< CONST >> dislocation flux', &
                                            rhoDotFlux(1:ns,1:8) * timestep
    write(6,'(a,/,10(12x,12(e12.5,1x),/))') '<< CONST >> dipole formation by glide', &
                                            rhoDotSingle2DipoleGlide * timestep
    write(6,'(a,/,10(12x,12(e12.5,1x),/))') '<< CONST >> athermal dipole annihilation', &
                                            rhoDotAthermalAnnihilation * timestep
    write(6,'(a,/,2(12x,12(e12.5,1x),/))')  '<< CONST >> thermally activated dipole annihilation', &
                                            rhoDotThermalAnnihilation(1:ns,9:10) * timestep
    write(6,'(a,/,10(12x,12(e12.5,1x),/))') '<< CONST >> total density change', &
                                            rhoDot * timestep
    write(6,'(a,/,10(12x,12(f12.5,1x),/))') '<< CONST >> relative density change', &
                                            rhoDot(1:ns,1:8) * timestep / (abs(rhoSglOriginal)+1.0e-10), &
                                            rhoDot(1:ns,9:10) * timestep / (rhoDipOriginal+1.0e-10)
    write(6,*)
  endif
#endif


if (    any(rhoSglOriginal(1:ns,1:4) + rhoDot(1:ns,1:4) * timestep < -aTolRho(instance)) &
   .or. any(rhoDipOriginal(1:ns,1:2) + rhoDot(1:ns,9:10) * timestep < -aTolRho(instance))) then
#ifndef _OPENMP
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0_pInt) then 
    write(6,'(a,i5,a,i2)') '<< CONST >> evolution rate leads to negative density at el ',el,' ip ',ip
    write(6,'(a)') '<< CONST >> enforcing cutback !!!'
  endif
#endif
  plasticState(p)%dotState = DAMASK_NaN
  return
else
  forall (s = 1:ns, t = 1_pInt:4_pInt) 
    plasticState(p)%dotState(iRhoU(s,t,instance),o) = rhoDot(s,t)
    plasticState(p)%dotState(iRhoB(s,t,instance),o) = rhoDot(s,t+4_pInt)
  endforall
  forall (s = 1:ns, c = 1_pInt:2_pInt) &
    plasticState(p)%dotState(iRhoD(s,c,instance),o) = rhoDot(s,c+8_pInt)
  forall (s = 1:ns) &
    plasticState(p)%dotState(iGamma(s,instance),o) = sum(gdot(s,1:4))
endif

end subroutine plastic_nonlocal_dotState


!*********************************************************************
!* COMPATIBILITY UPDATE                                              *
!* Compatibility is defined as normalized product of signed cosine   *
!* of the angle between the slip plane normals and signed cosine of  *
!* the angle between the slip directions. Only the largest values    *
!* that sum up to a total of 1 are considered, all others are set to *
!* zero.                                                             *
!*********************************************************************
subroutine plastic_nonlocal_updateCompatibility(orientation,i,e)

use math, only:       math_mul3x3, &
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
                      lattice_sd, &
                      lattice_qDisorientation

implicit none

!* input variables
integer(pInt), intent(in) ::                    i, &                                                ! ip index
                                                e                                                   ! element index
real(pReal), dimension(4,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                                orientation                                         ! crystal orientation in quaternions
                                            
!* local variables
integer(pInt)                                   Nneighbors, &                                       ! number of neighbors
                                                n, &                                                ! neighbor index 
                                                neighbor_e, &                                       ! element index of my neighbor
                                                neighbor_i, &                                       ! integration point index of my neighbor
                                                ph, &
                                                neighbor_phase, &
                                                textureID, &
                                                neighbor_textureID, &
                                                instance, &                                         ! instance of plasticity
                                                ns, &                                               ! number of active slip systems
                                                s1, &                                               ! slip system index (me)
                                                s2                                                  ! slip system index (my neighbor)
real(pReal), dimension(4) ::                    absoluteMisorientation                              ! absolute misorientation (without symmetry) between me and my neighbor
real(pReal), dimension(2,totalNslip(phase_plasticityInstance(material_phase(1,i,e))),&
                         totalNslip(phase_plasticityInstance(material_phase(1,i,e))),&
                         FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,e))))) :: &  
                                                my_compatibility                                    ! my_compatibility for current element and ip
real(pReal), dimension(3,totalNslip(phase_plasticityInstance(material_phase(1,i,e)))) :: &  
                                                slipNormal, &
                                                slipDirection
real(pReal)                                     my_compatibilitySum, &
                                                thresholdValue, &
                                                nThresholdValues
logical, dimension(totalNslip(phase_plasticityInstance(material_phase(1,i,e)))) :: & 
                                                belowThreshold


Nneighbors = FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,e))))
ph = material_phase(1,i,e)
textureID = material_texture(1,i,e)
instance = phase_plasticityInstance(ph)
ns = totalNslip(instance)
slipNormal(1:3,1:ns) = lattice_sn(1:3, slipSystemLattice(1:ns,instance), ph)
slipDirection(1:3,1:ns) = lattice_sd(1:3, slipSystemLattice(1:ns,instance), ph)


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
      my_compatibility(1:2,s1,s1,n) = sqrt(surfaceTransmissivity(instance))
    cycle
  endif
  
  
  !* PHASE BOUNDARY
  !* If we encounter a different nonlocal "cpfem" phase at the neighbor, 
  !* we consider this to be a real "physical" phase boundary, so completely incompatible.
  !* If one of the two "CPFEM" phases has a local plasticity law, 
  !* we do not consider this to be a phase boundary, so completely compatible.
  
  neighbor_phase = material_phase(1,neighbor_i,neighbor_e)
  if (neighbor_phase /= ph) then
    if (.not. phase_localPlasticity(neighbor_phase) .and. .not. phase_localPlasticity(ph)) then
      forall(s1 = 1_pInt:ns) &
        my_compatibility(1:2,s1,s1,n) = 0.0_pReal ! = sqrt(0.0)
    endif
    cycle
  endif

    
  !* GRAIN BOUNDARY !
  !* fixed transmissivity for adjacent ips with different texture (only if explicitly given in material.config)

  if (grainboundaryTransmissivity(instance) >= 0.0_pReal) then
    neighbor_textureID = material_texture(1,neighbor_i,neighbor_e)
    if (neighbor_textureID /= textureID) then
      if (.not. phase_localPlasticity(neighbor_phase)) then
        forall(s1 = 1_pInt:ns) &
          my_compatibility(1:2,s1,s1,n) = sqrt(grainboundaryTransmissivity(instance))
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
    absoluteMisorientation = lattice_qDisorientation(orientation(1:4,1,i,e), &
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

end subroutine plastic_nonlocal_updateCompatibility

!*********************************************************************
!* calculates quantities characterizing the microstructure           *
!*********************************************************************
function plastic_nonlocal_dislocationstress(Fe, ip, el)
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
                    plasticState, &
                    mappingConstitutive,&
                    phase_localPlasticity, &
                    phase_plasticityInstance
use lattice, only:  lattice_mu, &
                    lattice_nu

implicit none

!*** input variables
integer(pInt), intent(in) ::    ip, &                                                               !< current integration point
                                el                                                                  !< current element
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                Fe                                                                  !< elastic deformation gradient

!*** output variables
real(pReal), dimension(3,3) ::  plastic_nonlocal_dislocationstress

!*** local variables
integer(pInt)                   neighbor_el, &                                                      !< element number of neighbor material point
                                neighbor_ip, &                                                      !< integration point of neighbor material point
                                instance, &                                                         !< my instance of this plasticity
                                neighbor_instance, &                                                !< instance of this plasticity of neighbor material point
                                ph, &
                                neighbor_phase, &
                                ns, &                                                               !< total number of active slip systems at my material point
                                neighbor_ns, &                                                      !< total number of active slip systems at neighbor material point
                                c, &                                                                !< index of dilsocation character (edge, screw)
                                s, &                                                                !< slip system index
                                o,&                                                                 !< offset shortcut
                                no,&                                                                !< neighbour offset shortcut
                                p,&                                                                 !< phase shortcut
                                np,&                                                                !< neighbour phase shortcut
                                t, &                                                                !< index of dilsocation type (e+, e-, s+, s-, used e+, used e-, used s+, used s-)
                                dir, &
                                deltaX, deltaY, deltaZ, &
                                side, &
                                j
integer(pInt), dimension(2,3) :: periodicImages
real(pReal)                     x, y, z, &                                                          !< coordinates of connection vector in neighbor lattice frame
                                xsquare, ysquare, zsquare, &                                        !< squares of respective coordinates
                                distance, &                                                         !< length of connection vector
                                segmentLength, &                                                    !< segment length of dislocations
                                lambda, &
                                R, Rsquare, Rcube, &
                                denominator, &
                                flipSign, &
                                neighbor_ipVolumeSideLength, &                                      
                                detFe
real(pReal), dimension(3) ::    connection, &                                                       !< connection vector between me and my neighbor in the deformed configuration
                                connection_neighborLattice, &                                       !< connection vector between me and my neighbor in the lattice configuration of my neighbor
                                connection_neighborSlip, &                                          !< connection vector between me and my neighbor in the slip system frame of my neighbor
                                maxCoord, minCoord, &
                                meshSize, &
                                coords, &                                                           !< x,y,z coordinates of cell center of ip volume
                                neighbor_coords                                                     !< x,y,z coordinates of cell center of neighbor ip volume
real(pReal), dimension(3,3) ::  sigma, &                                                            !< dislocation stress for one slip system in neighbor material point's slip system frame
                                Tdislo_neighborLattice, &                                           !< dislocation stress as 2nd Piola-Kirchhoff stress at neighbor material point
                                invFe, &                                                            !< inverse of my elastic deformation gradient
                                neighbor_invFe, &
                                neighborLattice2myLattice                                           !< mapping from neighbor MPs lattice configuration to my lattice configuration
real(pReal), dimension(2,2,maxval(totalNslip)) :: &
                                neighbor_rhoExcess                                                  !< excess density at neighbor material point (edge/screw,mobile/dead,slipsystem)
real(pReal), dimension(2,maxval(totalNslip)) :: &
                                rhoExcessDead
real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),8) :: &
                                rhoSgl                                                              ! single dislocation density (edge+, edge-, screw+, screw-, used edge+, used edge-, used screw+, used screw-)
logical                         inversionError

ph = material_phase(1_pInt,ip,el)
instance = phase_plasticityInstance(ph)
ns = totalNslip(instance)
p = mappingConstitutive(2,1,ip,el)
o = mappingConstitutive(1,1,ip,el)



!*** get basic states

forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)
  rhoSgl(s,t) = max(plasticState(p)%state(iRhoU(s,t,instance),o), 0.0_pReal)                              ! ensure positive single mobile densities
  rhoSgl(s,t+4_pInt) = plasticState(p)%state(iRhoB(s,t,instance),o)
endforall



!*** calculate the dislocation stress of the neighboring excess dislocation densities
!*** zero for material points of local plasticity

plastic_nonlocal_dislocationstress = 0.0_pReal

if (.not. phase_localPlasticity(ph)) then
  call math_invert33(Fe(1:3,1:3,1_pInt,ip,el), invFe, detFe, inversionError)

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
      neighbor_phase = material_phase(1_pInt,neighbor_ip,neighbor_el)
      np = mappingConstitutive(2,1,neighbor_ip,neighbor_el)
      no = mappingConstitutive(1,1,neighbor_ip,neighbor_el)

      if (phase_localPlasticity(neighbor_phase)) then
        cycle
      endif
      neighbor_instance = phase_plasticityInstance(neighbor_phase)
      neighbor_ns = totalNslip(neighbor_instance)
      call math_invert33(Fe(1:3,1:3,1,neighbor_ip,neighbor_el), neighbor_invFe, detFe, inversionError)
      neighbor_ipVolumeSideLength = mesh_ipVolume(neighbor_ip,neighbor_el) ** (1.0_pReal/3.0_pReal) ! reference volume used here
      
      forall (s = 1_pInt:neighbor_ns, c = 1_pInt:2_pInt)
        neighbor_rhoExcess(c,1,s) = plasticState(np)%state(iRhoU(s,2*c-1,neighbor_instance),no) &      ! positive mobiles
                                  - plasticState(np)%state(iRhoU(s,2*c,neighbor_instance),no)          ! negative mobiles
        neighbor_rhoExcess(c,2,s) = abs(plasticState(np)%state(iRhoB(s,2*c-1,neighbor_instance),no)) & ! positive deads
                                  - abs(plasticState(np)%state(iRhoB(s,2*c,neighbor_instance),no))     ! negative deads

      endforall
      Tdislo_neighborLattice = 0.0_pReal
      do deltaX = periodicImages(1,1),periodicImages(2,1)
        do deltaY = periodicImages(1,2),periodicImages(2,2)
          do deltaZ = periodicImages(1,3),periodicImages(2,3)
            
            
            !* regular case
            
            if (neighbor_el /= el .or. neighbor_ip /= ip &
                .or. deltaX /= 0_pInt .or. deltaY /= 0_pInt .or. deltaZ /= 0_pInt) then
            
              neighbor_coords = mesh_cellCenterCoordinates(neighbor_ip,neighbor_el) &
                                 + [real(deltaX,pReal), real(deltaY,pReal), real(deltaZ,pReal)] * meshSize
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
                if (all(abs(neighbor_rhoExcess(:,:,s)) < significantRho(instance))) cycle           ! not significant
                
                
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
                    x = connection_neighborSlip(1) &
                      + sign(0.5_pReal * segmentLength, &
                              plasticState(np)%state(iRhoB(s,1,neighbor_instance),no) &
                              - plasticState(np)%state(iRhoB(s,2,neighbor_instance),no))

                    xsquare = x * x
                  endif
                   
                  flipSign = sign(1.0_pReal, -y)
                  do side = 1_pInt,-1_pInt,-2_pInt
                    lambda = real(side,pReal) * 0.5_pReal * segmentLength - y
                    R = sqrt(xsquare + zsquare + lambda * lambda)
                    Rsquare = R * R
                    Rcube = Rsquare * R 
                    denominator = R * (R + flipSign * lambda)
                    if (denominator == 0.0_pReal) exit ipLoop
                      
                    sigma(1,1) = sigma(1,1) - real(side,pReal) &
                                            * flipSign * z / denominator &
                                            * (1.0_pReal + xsquare / Rsquare + xsquare / denominator) &
                                            * neighbor_rhoExcess(1,j,s)
                    sigma(2,2) = sigma(2,2) - real(side,pReal) &
                                            * (flipSign * 2.0_pReal * lattice_nu(ph) * z / denominator + z * lambda / Rcube) &
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
                                            * (lattice_nu(ph) / R - zsquare / Rcube) * neighbor_rhoExcess(1,j,s)
                  enddo
                enddo 
                
                !* screw contribution to stress
                
                x = connection_neighborSlip(1)   ! have to restore this value, because position might have been adapted for edge deads before
                do j = 1_pInt,2_pInt
                  if (abs(neighbor_rhoExcess(2,j,s)) < significantRho(instance)) then
                    cycle 
                  elseif (j > 1_pInt) then
                    y = connection_neighborSlip(2) &
                      + sign(0.5_pReal * segmentLength, &
                             plasticState(np)%state(iRhoB(s,3,neighbor_instance),no) &
                             - plasticState(np)%state(iRhoB(s,4,neighbor_instance),no))
                    ysquare = y * y
                  endif

                  flipSign = sign(1.0_pReal, x)
                  do side = 1_pInt,-1_pInt,-2_pInt
                    lambda = x + real(side,pReal) * 0.5_pReal * segmentLength
                    R = sqrt(ysquare + zsquare + lambda * lambda)
                    Rsquare = R * R
                    Rcube = Rsquare * R 
                    denominator = R * (R + flipSign * lambda)
                    if (denominator == 0.0_pReal) exit ipLoop
                    
                    sigma(1,2) = sigma(1,2) - real(side,pReal) * flipSign * z &
                                                               * (1.0_pReal - lattice_nu(ph)) / denominator &
                                                               * neighbor_rhoExcess(2,j,s)
                    sigma(1,3) = sigma(1,3) + real(side,pReal) * flipSign * y &
                                                               * (1.0_pReal - lattice_nu(ph)) / denominator &
                                                               * neighbor_rhoExcess(2,j,s)
                  enddo
                enddo
               
                if (all(abs(sigma) < 1.0e-10_pReal)) cycle ! SIGMA IS NOT A REAL STRESS, THATS WHY WE NEED A REALLY SMALL VALUE HERE

                !* copy symmetric parts
                
                sigma(2,1) = sigma(1,2)
                sigma(3,1) = sigma(1,3)
                sigma(3,2) = sigma(2,3)

                
                !* scale stresses and map them into the neighbor material point's lattice configuration
                
                sigma = sigma * lattice_mu(neighbor_phase) * burgers(s,neighbor_instance) &
                              / (4.0_pReal * pi * (1.0_pReal - lattice_nu(neighbor_phase))) &
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

                rhoExcessDead(c,s) = plasticState(p)%state(iRhoB(s,2*c-1,instance),o) &             ! positive deads (here we use symmetry: if this has negative sign it is 
                                                                                                    !treated as negative density at positive position instead of positive 
                                                                                                    !density at negative position)
                                   + plasticState(p)%state(iRhoB(s,2*c,instance),o)                 ! negative deads (here we use symmetry: if this has negative sign it is 
                                                                                                    !treated as positive density at positive position instead of negative 
                                                                                                    !density at negative position)
              do s = 1_pInt,ns
                if (all(abs(rhoExcessDead(:,s)) < significantRho(instance))) cycle                  ! not significant
                sigma = 0.0_pReal                                                                   ! all components except for sigma13 are zero
                sigma(1,3) = - (rhoExcessDead(1,s) + rhoExcessDead(2,s) * (1.0_pReal - lattice_nu(ph))) &
                           * neighbor_ipVolumeSideLength * lattice_mu(ph) * burgers(s,instance) &
                           / (sqrt(2.0_pReal) * pi * (1.0_pReal - lattice_nu(ph)))
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
      plastic_nonlocal_dislocationstress = plastic_nonlocal_dislocationstress &
                                              + math_mul33x33(neighborLattice2myLattice, &
                                                math_mul33x33(Tdislo_neighborLattice, &
                                                math_transpose33(neighborLattice2myLattice)))
                        
    enddo ipLoop
  enddo ! element loop
    
endif

end function plastic_nonlocal_dislocationstress

 
!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_nonlocal_postResults(Tstar_v,Fe,ip,el)
 use math, only: &
   math_mul6x6, &
   math_mul33x3, &
   math_mul33x33, &
   pi
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   mappingConstitutive, &
   plasticState, &
   phase_plasticityInstance, &
   phase_Noutput
 use lattice, only: &
   lattice_Sslip_v, &
   lattice_sd, &
   lattice_st, &
   lattice_sn, &
   lattice_mu, &
   lattice_nu

 implicit none
 real(pReal),   dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),   dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
   Fe                                                                                               !< elastic deformation gradient
 integer(pInt),                intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 
 real(pReal),   dimension(plastic_nonlocal_sizePostResults(&
                                         phase_plasticityInstance(material_phase(1_pInt,ip,el)))) :: &
   plastic_nonlocal_postResults
 
 integer(pInt) :: &
   ph, &
   instance, &                                                                                      !< current instance of this plasticity
   ns, &                                                                                            !< short notation for the total number of active slip systems
   c, &                                                                                             !< character of dislocation
   cs, &                                                                                            !< constitutive result index
   o, &                                                                                             !< index of current output
   of,&                                                                                             !< offset shortcut
   t, &                                                                                             !< type of dislocation
   s, &                                                                                             !< index of my current slip system
   sLattice                                                                                         !< index of my current slip system according to lattice order
 real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),8) :: &
   rhoSgl, &                                                                                        !< current single dislocation densities (positive/negative screw and edge without dipoles)
   rhoDotSgl                                                                                        !< evolution rate of single dislocation densities (positive/negative screw and edge without dipoles)
 real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),4) :: &
   gdot, &                                                                                          !< shear rates
   v                                                                                                !< velocities
 real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el)))) :: &
   rhoForest, &                                                                                     !< forest dislocation density
   tauThreshold, &                                                                                  !< threshold shear stress
   tau, &                                                                                           !< current resolved shear stress
   tauBack                                                                                          !< back stress from pileups on same slip system
 real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),2) :: &
   rhoDip, &                                                                                        !< current dipole dislocation densities (screw and edge dipoles)
   rhoDotDip, &                                                                                     !< evolution rate of dipole dislocation densities (screw and edge dipoles)
   dLower, &                                                                                        !< minimum stable dipole distance for edges and screws
   dUpper                                                                                           !< current maximum stable dipole distance for edges and screws
 real(pReal), dimension(3,totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el))),2) :: &
   m, &                                                                                             !< direction of dislocation motion for edge and screw (unit vector)
   m_currentconf                                                                                    !< direction of dislocation motion for edge and screw (unit vector) in current configuration
 real(pReal), dimension(3,totalNslip(phase_plasticityInstance(material_phase(1_pInt,ip,el)))) :: &
   n_currentconf                                                                                    !< slip system normal (unit vector) in current configuration
 real(pReal), dimension(3,3) :: &
   sigma

ph  = mappingConstitutive(2,1,ip,el)
of = mappingConstitutive(1,1,ip,el)
instance = phase_plasticityInstance(ph)
ns = totalNslip(instance)

cs = 0_pInt
plastic_nonlocal_postResults = 0.0_pReal


!* short hand notations for state variables

forall (s = 1_pInt:ns, t = 1_pInt:4_pInt)
  rhoSgl(s,t)           = plasticState(ph)%State(iRhoU(s,t,instance),of)
  rhoSgl(s,t+4_pInt)    = plasticState(ph)%State(iRhoB(s,t,instance),of)
  v(s,t)                = plasticState(ph)%State(iV(s,t,instance),of)
  rhoDotSgl(s,t)        = plasticState(ph)%dotState(iRhoU(s,t,instance),of)
  rhoDotSgl(s,t+4_pInt) = plasticState(ph)%dotState(iRhoB(s,t,instance),of)
endforall
forall (s = 1_pInt:ns, c = 1_pInt:2_pInt)
  rhoDip(s,c)    = plasticState(ph)%State(iRhoD(s,c,instance),of)
  rhoDotDip(s,c) = plasticState(ph)%dotState(iRhoD(s,c,instance),of)
endforall
rhoForest    = plasticState(ph)%State(iRhoF(1:ns,instance),of)
tauThreshold = plasticState(ph)%State(iTauF(1:ns,instance),of)
tauBack      = plasticState(ph)%State(iTauB(1:ns,instance),of)

!* Calculate shear rate

forall (t = 1_pInt:4_pInt) &
  gdot(1:ns,t) = rhoSgl(1:ns,t) * burgers(1:ns,instance) * v(1:ns,t)
  

!* calculate limits for stable dipole height

do s = 1_pInt,ns
  sLattice = slipSystemLattice(s,instance)
  tau(s) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,1,sLattice,ph)) + tauBack(s)
  if (abs(tau(s)) < 1.0e-15_pReal) tau(s) = 1.0e-15_pReal
enddo

dLower = minDipoleHeight(1:ns,1:2,instance)
dUpper(1:ns,1) = lattice_mu(ph) * burgers(1:ns,instance) &
               / (8.0_pReal * pi * (1.0_pReal - lattice_nu(ph)) * abs(tau))
dUpper(1:ns,2) = lattice_mu(ph) * burgers(1:ns,instance) &
               / (4.0_pReal * pi * abs(tau))
forall (c = 1_pInt:2_pInt)
  where(sqrt(rhoSgl(1:ns,2*c-1)+rhoSgl(1:ns,2*c)+&
        abs(rhoSgl(1:ns,2*c+3))+abs(rhoSgl(1:ns,2*c+4))+rhoDip(1:ns,c)) >= tiny(0.0_pReal)) &
    dUpper(1:ns,c) = min(1.0_pReal / sqrt(rhoSgl(1:ns,2*c-1) + rhoSgl(1:ns,2*c) & 
                       + abs(rhoSgl(1:ns,2*c+3)) + abs(rhoSgl(1:ns,2*c+4)) + rhoDip(1:ns,c)), &
                       dUpper(1:ns,c))
end forall
dUpper = max(dUpper,dLower)


!*** dislocation motion

m(1:3,1:ns,1) = lattice_sd(1:3,slipSystemLattice(1:ns,instance),ph)
m(1:3,1:ns,2) = -lattice_st(1:3,slipSystemLattice(1:ns,instance),ph)
forall (c = 1_pInt:2_pInt, s = 1_pInt:ns) &
  m_currentconf(1:3,s,c) = math_mul33x3(Fe(1:3,1:3,1_pInt,ip,el), m(1:3,s,c))
forall (s = 1_pInt:ns) &
  n_currentconf(1:3,s) = math_mul33x3(Fe(1:3,1:3,1_pInt,ip,el), &
                                      lattice_sn(1:3,slipSystemLattice(s,instance),ph))


outputsLoop: do o = 1_pInt,plastic_nonlocal_Noutput(instance)
  select case(plastic_nonlocal_outputID(o,instance))
    case (rho_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl),2) + sum(rhoDip,2)
      cs = cs + ns
      
    case (rho_sgl_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl),2)
      cs = cs + ns
      
    case (rho_sgl_mobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl(1:ns,1:4)),2)
      cs = cs + ns
      
    case (rho_sgl_immobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoSgl(1:ns,5:8),2)
      cs = cs + ns
      
    case (rho_dip_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDip,2)
      cs = cs + ns
      
    case (rho_edge_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl(1:ns,[1,2,5,6])),2) + rhoDip(1:ns,1)
      cs = cs + ns
      
    case (rho_sgl_edge_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl(1:ns,[1,2,5,6])),2)
      cs = cs + ns
      
    case (rho_sgl_edge_mobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoSgl(1:ns,1:2),2)
      cs = cs + ns
      
    case (rho_sgl_edge_immobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoSgl(1:ns,5:6),2)
      cs = cs + ns
      
    case (rho_sgl_edge_pos_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,1) + abs(rhoSgl(1:ns,5))
      cs = cs + ns
      
    case (rho_sgl_edge_pos_mobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,1)
      cs = cs + ns
      
    case (rho_sgl_edge_pos_immobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,5)
      cs = cs + ns
      
    case (rho_sgl_edge_neg_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,2) + abs(rhoSgl(1:ns,6))
      cs = cs + ns
      
    case (rho_sgl_edge_neg_mobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,2)
      cs = cs + ns
      
    case (rho_sgl_edge_neg_immobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,6)
      cs = cs + ns
      
    case (rho_dip_edge_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDip(1:ns,1)
      cs = cs + ns
      
    case (rho_screw_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl(1:ns,[3,4,7,8])),2) + rhoDip(1:ns,2)
      cs = cs + ns
      
    case (rho_sgl_screw_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(abs(rhoSgl(1:ns,[3,4,7,8])),2)
      cs = cs + ns
            
    case (rho_sgl_screw_mobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoSgl(1:ns,3:4),2)
      cs = cs + ns
      
    case (rho_sgl_screw_immobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoSgl(1:ns,7:8),2)
      cs = cs + ns
      
    case (rho_sgl_screw_pos_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,3) + abs(rhoSgl(1:ns,7))
      cs = cs + ns
      
    case (rho_sgl_screw_pos_mobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,3)
      cs = cs + ns
      
    case (rho_sgl_screw_pos_immobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,7)
      cs = cs + ns
      
    case (rho_sgl_screw_neg_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,4) + abs(rhoSgl(1:ns,8))
      cs = cs + ns

    case (rho_sgl_screw_neg_mobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,4)
      cs = cs + ns

    case (rho_sgl_screw_neg_immobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,8)
      cs = cs + ns

    case (rho_dip_screw_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDip(1:ns,2)
      cs = cs + ns
      
    case (excess_rho_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = (rhoSgl(1:ns,1) + abs(rhoSgl(1:ns,5))) &
                                                         - (rhoSgl(1:ns,2) + abs(rhoSgl(1:ns,6))) &
                                                         + (rhoSgl(1:ns,3) + abs(rhoSgl(1:ns,7))) &
                                                         - (rhoSgl(1:ns,4) + abs(rhoSgl(1:ns,8)))
      cs = cs + ns
      
    case (excess_rho_edge_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = (rhoSgl(1:ns,1) + abs(rhoSgl(1:ns,5))) &
                                                         - (rhoSgl(1:ns,2) + abs(rhoSgl(1:ns,6)))
      cs = cs + ns
      
    case (excess_rho_screw_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = (rhoSgl(1:ns,3) + abs(rhoSgl(1:ns,7))) &
                                                         - (rhoSgl(1:ns,4) + abs(rhoSgl(1:ns,8)))
      cs = cs + ns
      
    case (rho_forest_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoForest
      cs = cs + ns
    
    case (delta_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = 1.0_pReal / sqrt(sum(abs(rhoSgl),2) + sum(rhoDip,2))
      cs = cs + ns
      
    case (delta_sgl_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = 1.0_pReal / sqrt(sum(abs(rhoSgl),2))
      cs = cs + ns
      
    case (delta_dip_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = 1.0_pReal / sqrt(sum(rhoDip,2))
      cs = cs + ns
      
    case (shearrate_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(gdot,2)
      cs = cs + ns
      
    case (resolvedstress_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = tau
      cs = cs + ns
      
    case (resolvedstress_back_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = tauBack
      cs = cs + ns
      
    case (resolvedstress_external_ID)
      do s = 1_pInt,ns  
        sLattice = slipSystemLattice(s,instance)
        plastic_nonlocal_postResults(cs+s) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,1,sLattice,ph))
      enddo
      cs = cs + ns
      
    case (resistance_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = tauThreshold
      cs = cs + ns
    
    case (rho_dot_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotSgl(1:ns,1:4),2) &
                                                         + sum(rhoDotSgl(1:ns,5:8)*sign(1.0_pReal,rhoSgl(1:ns,5:8)),2) &
                                                         + sum(rhoDotDip,2)
      cs = cs + ns
      
    case (rho_dot_sgl_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotSgl(1:ns,1:4),2) &
                                                         + sum(rhoDotSgl(1:ns,5:8)*sign(1.0_pReal,rhoSgl(1:ns,5:8)),2)
      cs = cs + ns
      
    case (rho_dot_sgl_mobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotSgl(1:ns,1:4),2)
      cs = cs + ns
      
    case (rho_dot_dip_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotDip,2)
      cs = cs + ns
    
    case (rho_dot_gen_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotMultiplicationOutput(1:ns,1,1_pInt,ip,el) &
                                                         + rhoDotMultiplicationOutput(1:ns,2,1_pInt,ip,el)
      cs = cs + ns

    case (rho_dot_gen_edge_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotMultiplicationOutput(1:ns,1,1_pInt,ip,el)
      cs = cs + ns

    case (rho_dot_gen_screw_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotMultiplicationOutput(1:ns,2,1_pInt,ip,el)
      cs = cs + ns
      
    case (rho_dot_sgl2dip_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotSingle2DipoleGlideOutput(1:ns,1,1_pInt,ip,el) &
                                                         + rhoDotSingle2DipoleGlideOutput(1:ns,2,1_pInt,ip,el)
      cs = cs + ns
    
    case (rho_dot_sgl2dip_edge_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotSingle2DipoleGlideOutput(1:ns,1,1_pInt,ip,el)
      cs = cs + ns
    
    case (rho_dot_sgl2dip_screw_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotSingle2DipoleGlideOutput(1:ns,2,1_pInt,ip,el)
      cs = cs + ns
    
    case (rho_dot_ann_ath_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotAthermalAnnihilationOutput(1:ns,1,1_pInt,ip,el) & 
                                                         + rhoDotAthermalAnnihilationOutput(1:ns,2,1_pInt,ip,el)
      cs = cs + ns
      
    case (rho_dot_ann_the_ID) 
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotThermalAnnihilationOutput(1:ns,1,1_pInt,ip,el) & 
                                                         + rhoDotThermalAnnihilationOutput(1:ns,2,1_pInt,ip,el)
      cs = cs + ns

    case (rho_dot_ann_the_edge_ID) 
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotThermalAnnihilationOutput(1:ns,1,1_pInt,ip,el) 
      cs = cs + ns

    case (rho_dot_ann_the_screw_ID) 
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotThermalAnnihilationOutput(1:ns,2,1_pInt,ip,el)
      cs = cs + ns

    case (rho_dot_edgejogs_ID) 
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoDotEdgeJogsOutput(1:ns,1_pInt,ip,el)
      cs = cs + ns

    case (rho_dot_flux_mobile_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotFluxOutput(1:ns,1:4,1_pInt,ip,el),2)
      cs = cs + ns
    
    case (rho_dot_flux_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotFluxOutput(1:ns,1:4,1_pInt,ip,el),2) &
                          + sum(rhoDotFluxOutput(1:ns,5:8,1_pInt,ip,el)*sign(1.0_pReal,rhoSgl(1:ns,5:8)),2)
      cs = cs + ns
    
    case (rho_dot_flux_edge_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotFluxOutput(1:ns,1:2,1_pInt,ip,el),2) &
                          + sum(rhoDotFluxOutput(1:ns,5:6,1_pInt,ip,el)*sign(1.0_pReal,rhoSgl(1:ns,5:6)),2)
      cs = cs + ns
      
    case (rho_dot_flux_screw_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = sum(rhoDotFluxOutput(1:ns,3:4,1_pInt,ip,el),2) &
                          + sum(rhoDotFluxOutput(1:ns,7:8,1_pInt,ip,el)*sign(1.0_pReal,rhoSgl(1:ns,7:8)),2)
      cs = cs + ns
            
    case (velocity_edge_pos_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = v(1:ns,1)
      cs = cs + ns
    
    case (velocity_edge_neg_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = v(1:ns,2)
      cs = cs + ns
    
    case (velocity_screw_pos_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = v(1:ns,3)
      cs = cs + ns
    
    case (velocity_screw_neg_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = v(1:ns,4)
      cs = cs + ns
    
    case (slipdirectionx_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = m_currentconf(1,1:ns,1)
      cs = cs + ns
    
    case (slipdirectiony_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = m_currentconf(2,1:ns,1)
      cs = cs + ns
    
    case (slipdirectionz_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = m_currentconf(3,1:ns,1)
      cs = cs + ns
    
    case (slipnormalx_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = n_currentconf(1,1:ns)
      cs = cs + ns
    
    case (slipnormaly_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = n_currentconf(2,1:ns)
      cs = cs + ns
    
    case (slipnormalz_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = n_currentconf(3,1:ns)
      cs = cs + ns
    
    case (fluxdensity_edge_posx_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,1) * v(1:ns,1) * m_currentconf(1,1:ns,1)
      cs = cs + ns
    
    case (fluxdensity_edge_posy_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,1) * v(1:ns,1) * m_currentconf(2,1:ns,1)
      cs = cs + ns
    
    case (fluxdensity_edge_posz_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,1) * v(1:ns,1) * m_currentconf(3,1:ns,1)
      cs = cs + ns
    
    case (fluxdensity_edge_negx_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = - rhoSgl(1:ns,2) * v(1:ns,2) * m_currentconf(1,1:ns,1)
      cs = cs + ns
    
    case (fluxdensity_edge_negy_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = - rhoSgl(1:ns,2) * v(1:ns,2) * m_currentconf(2,1:ns,1)
      cs = cs + ns
    
    case (fluxdensity_edge_negz_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = - rhoSgl(1:ns,2) * v(1:ns,2) * m_currentconf(3,1:ns,1)
      cs = cs + ns
    
    case (fluxdensity_screw_posx_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,3) * v(1:ns,3) * m_currentconf(1,1:ns,2)
      cs = cs + ns
    
    case (fluxdensity_screw_posy_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,3) * v(1:ns,3) * m_currentconf(2,1:ns,2)
      cs = cs + ns
    
    case (fluxdensity_screw_posz_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = rhoSgl(1:ns,3) * v(1:ns,3) * m_currentconf(3,1:ns,2)
      cs = cs + ns
    
    case (fluxdensity_screw_negx_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = - rhoSgl(1:ns,4) * v(1:ns,4) * m_currentconf(1,1:ns,2)
      cs = cs + ns
    
    case (fluxdensity_screw_negy_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = - rhoSgl(1:ns,4) * v(1:ns,4) * m_currentconf(2,1:ns,2)
      cs = cs + ns
    
    case (fluxdensity_screw_negz_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = - rhoSgl(1:ns,4) * v(1:ns,4) * m_currentconf(3,1:ns,2)
      cs = cs + ns
    
    case (maximumdipoleheight_edge_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = dUpper(1:ns,1)
      cs = cs + ns
      
    case (maximumdipoleheight_screw_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = dUpper(1:ns,2)
      cs = cs + ns
    
    case(dislocationstress_ID)
      sigma = plastic_nonlocal_dislocationstress(Fe, ip, el)
      plastic_nonlocal_postResults(cs+1_pInt) = sigma(1,1)
      plastic_nonlocal_postResults(cs+2_pInt) = sigma(2,2)
      plastic_nonlocal_postResults(cs+3_pInt) = sigma(3,3)
      plastic_nonlocal_postResults(cs+4_pInt) = sigma(1,2)
      plastic_nonlocal_postResults(cs+5_pInt) = sigma(2,3)
      plastic_nonlocal_postResults(cs+6_pInt) = sigma(3,1)
      cs = cs + 6_pInt
    
    case(accumulatedshear_ID)
      plastic_nonlocal_postResults(cs+1_pInt:cs+ns) = plasticState(ph)%state(iGamma(1:ns,instance),of)
      cs = cs + ns
    
  end select
enddo outputsLoop

end function plastic_nonlocal_postResults

end module plastic_nonlocal
