!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for phenomenological crystal plasticity formulation using a powerlaw 
!! fitting
!--------------------------------------------------------------------------------------------------
module plastic_phenopowerlaw
 use prec, only: &
   pReal,&
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,         public, protected :: &
   plastic_phenopowerlaw_sizePostResults                                                       !< cumulative size of post results

 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   plastic_phenopowerlaw_sizePostResult                                                        !< size of each post result output

 character(len=64),                   dimension(:,:),   allocatable, target, public :: & 
   plastic_phenopowerlaw_output                                                                !< name of each post result output

 integer(pInt),                       dimension(:),     allocatable, target, public :: &
   plastic_phenopowerlaw_Noutput                                                               !< number of outputs per instance of this constitution 

 integer(pInt),                       dimension(:),     allocatable,         public, protected :: &
   plastic_phenopowerlaw_totalNslip, &                                                         !< no. of slip system used in simulation
   plastic_phenopowerlaw_totalNtwin, &                                                         !< no. of twin system used in simulation
   plastic_phenopowerlaw_totalNtrans                                                           !< no. of trans system used in simulation

 integer(pInt),                       dimension(:,:),   allocatable,         private :: &
   plastic_phenopowerlaw_Nslip, &                                                              !< active number of slip systems per family (input parameter, per family)
   plastic_phenopowerlaw_Ntwin, &                                                              !< active number of twin systems per family (input parameter, per family)
   plastic_phenopowerlaw_Ntrans                                                                !< active number of trans systems per family (input parameter, per family)

 real(pReal),                         dimension(:),     allocatable,         private :: &
   plastic_phenopowerlaw_gdot0_slip, &                                                         !< reference shear strain rate for slip (input parameter)
   plastic_phenopowerlaw_gdot0_twin, &                                                         !< reference shear strain rate for twin (input parameter)
   plastic_phenopowerlaw_n_slip, &                                                             !< stress exponent for slip (input parameter)
   plastic_phenopowerlaw_n_twin, &                                                             !< stress exponent for twin (input parameter)
   plastic_phenopowerlaw_spr, &                                                                !< push-up factor for slip saturation due to twinning
   plastic_phenopowerlaw_twinB, &
   plastic_phenopowerlaw_twinC, &
   plastic_phenopowerlaw_twinD, &
   plastic_phenopowerlaw_twinE, &
   plastic_phenopowerlaw_h0_SlipSlip, &                                                        !< reference hardening slip - slip (input parameter)
   plastic_phenopowerlaw_h0_TwinSlip, &                                                        !< reference hardening twin - slip (input parameter)
   plastic_phenopowerlaw_h0_TwinTwin, &                                                        !< reference hardening twin - twin (input parameter)
   plastic_phenopowerlaw_a_slip, &
   plastic_phenopowerlaw_aTolResistance, &
   plastic_phenopowerlaw_aTolShear, &
   plastic_phenopowerlaw_aTolTwinfrac, &
   plastic_phenopowerlaw_aTolTransfrac, &
   plastic_phenopowerlaw_Cnuc, &                                                               !< coefficient for strain-induced martensite nucleation  
   plastic_phenopowerlaw_Cdwp, &                                                               !< coefficient for double well potential 
   plastic_phenopowerlaw_Cgro, &                                                               !< coefficient for stress-assisted martensite growth  
   plastic_phenopowerlaw_deltaG                                                                !< free energy difference between austensite and martensite [MPa]

 real(pReal),                         dimension(:,:),   allocatable,          private :: &
   plastic_phenopowerlaw_tau0_slip, &                                                          !< initial critical shear stress for slip (input parameter, per family)
   plastic_phenopowerlaw_tau0_twin, &                                                          !< initial critical shear stress for twin (input parameter, per family)
   plastic_phenopowerlaw_tausat_slip, &                                                        !< maximum critical shear stress for slip (input parameter, per family)
   plastic_phenopowerlaw_nonSchmidCoeff, &

   plastic_phenopowerlaw_interaction_SlipSlip, &                                               !< interaction factors slip - slip (input parameter)
   plastic_phenopowerlaw_interaction_SlipTwin, &                                               !< interaction factors slip - twin (input parameter)
   plastic_phenopowerlaw_interaction_TwinSlip, &                                               !< interaction factors twin - slip (input parameter)
   plastic_phenopowerlaw_interaction_TwinTwin                                                  !< interaction factors twin - twin (input parameter)

 real(pReal),                         dimension(:,:,:), allocatable,          private :: &
   plastic_phenopowerlaw_hardeningMatrix_SlipSlip, &
   plastic_phenopowerlaw_hardeningMatrix_SlipTwin, &
   plastic_phenopowerlaw_hardeningMatrix_TwinSlip, &
   plastic_phenopowerlaw_hardeningMatrix_TwinTwin

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 resistance_slip_ID, &
                 accumulatedshear_slip_ID, &
                 shearrate_slip_ID, &
                 resolvedstress_slip_ID, &
                 totalshear_ID, &
                 resistance_twin_ID, &
                 accumulatedshear_twin_ID, &
                 shearrate_twin_ID, &
                 resolvedstress_twin_ID, &
                 totalvolfrac_twin_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),   allocatable,          private :: & 
   plastic_phenopowerlaw_outputID                                                              !< ID of each post result output
 
 public :: &
   plastic_phenopowerlaw_init, &
   plastic_phenopowerlaw_LpAndItsTangent, &
   plastic_phenopowerlaw_dotState, &
   plastic_phenopowerlaw_postResults
 private :: &
   plastic_phenopowerlaw_aTolState, &
   plastic_phenopowerlaw_stateInit


contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_phenopowerlaw_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level, &
   debug_constitutive,&
   debug_levelBasic
 use math, only: &
   math_Mandel3333to66, &
   math_Voigt66to3333
 use IO, only: &
   IO_read, &
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
 use material, only: &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_PHENOPOWERLAW_label, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   material_phase, &
   plasticState, &
   MATERIAL_partPhase
 use lattice
 use numerics,only: &
   worldrank, &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: &
   maxNinstance, &
   instance,phase,j,k, f,o, &
   Nchunks_SlipSlip = 0_pInt, Nchunks_SlipTwin = 0_pInt, &
   Nchunks_TwinSlip = 0_pInt, Nchunks_TwinTwin = 0_pInt, &
   Nchunks_SlipFamilies = 0_pInt, Nchunks_TwinFamilies = 0_pInt, &
   Nchunks_TransFamilies = 0_pInt, Nchunks_nonSchmid = 0_pInt, &
   NipcMyPhase, &
   offset_slip, index_myFamily, index_otherFamily, &
   mySize=0_pInt,sizeState,sizeDotState, sizeDeltaState
 character(len=65536) :: &
   tag  = '', &
   line = ''
 real(pReal), dimension(:), allocatable :: tempPerSlip
 
 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_PHENOPOWERLAW_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_PHENOPOWERLAW_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 allocate(plastic_phenopowerlaw_sizePostResults(maxNinstance),               source=0_pInt)
 allocate(plastic_phenopowerlaw_sizePostResult(maxval(phase_Noutput),maxNinstance), &
                                                                                  source=0_pInt)
 allocate(plastic_phenopowerlaw_output(maxval(phase_Noutput),maxNinstance))
          plastic_phenopowerlaw_output               = ''
 allocate(plastic_phenopowerlaw_outputID(maxval(phase_Noutput),maxNinstance),source=undefined_ID)
 allocate(plastic_phenopowerlaw_Noutput(maxNinstance),                       source=0_pInt)
 allocate(plastic_phenopowerlaw_Nslip(lattice_maxNslipFamily,maxNinstance),  source=0_pInt)
 allocate(plastic_phenopowerlaw_Ntwin(lattice_maxNtwinFamily,maxNinstance),  source=0_pInt)
 allocate(plastic_phenopowerlaw_Ntrans(lattice_maxNtransFamily,maxNinstance),source=0_pInt)
 allocate(plastic_phenopowerlaw_totalNslip(maxNinstance),                    source=0_pInt)
 allocate(plastic_phenopowerlaw_totalNtwin(maxNinstance),                    source=0_pInt)
 allocate(plastic_phenopowerlaw_totalNtrans(maxNinstance),                   source=0_pInt)
 allocate(plastic_phenopowerlaw_gdot0_slip(maxNinstance),                    source=0.0_pReal)
 allocate(plastic_phenopowerlaw_n_slip(maxNinstance),                        source=0.0_pReal)
 allocate(plastic_phenopowerlaw_tau0_slip(lattice_maxNslipFamily,maxNinstance),  &
                                                                                  source=0.0_pReal)
 allocate(plastic_phenopowerlaw_tausat_slip(lattice_maxNslipFamily,maxNinstance), &
                                                                                  source=0.0_pReal)
 allocate(plastic_phenopowerlaw_gdot0_twin(maxNinstance),                    source=0.0_pReal)
 allocate(plastic_phenopowerlaw_n_twin(maxNinstance),                        source=0.0_pReal)
 allocate(plastic_phenopowerlaw_tau0_twin(lattice_maxNtwinFamily,maxNinstance), &
                                                                                  source=0.0_pReal)
 allocate(plastic_phenopowerlaw_spr(maxNinstance),                           source=0.0_pReal)
 allocate(plastic_phenopowerlaw_twinB(maxNinstance),                         source=0.0_pReal)
 allocate(plastic_phenopowerlaw_twinC(maxNinstance),                         source=0.0_pReal)
 allocate(plastic_phenopowerlaw_twinD(maxNinstance),                         source=0.0_pReal)
 allocate(plastic_phenopowerlaw_twinE(maxNinstance),                         source=0.0_pReal)
 allocate(plastic_phenopowerlaw_h0_SlipSlip(maxNinstance),                   source=0.0_pReal)
 allocate(plastic_phenopowerlaw_h0_TwinSlip(maxNinstance),                   source=0.0_pReal)
 allocate(plastic_phenopowerlaw_h0_TwinTwin(maxNinstance),                   source=0.0_pReal)
 allocate(plastic_phenopowerlaw_interaction_SlipSlip(lattice_maxNinteraction,maxNinstance), &
                                                                                  source=0.0_pReal)
 allocate(plastic_phenopowerlaw_interaction_SlipTwin(lattice_maxNinteraction,maxNinstance), &
                                                                                  source=0.0_pReal)
 allocate(plastic_phenopowerlaw_interaction_TwinSlip(lattice_maxNinteraction,maxNinstance), &
                                                                                  source=0.0_pReal)
 allocate(plastic_phenopowerlaw_interaction_TwinTwin(lattice_maxNinteraction,maxNinstance), &
                                                                                  source=0.0_pReal)
 allocate(plastic_phenopowerlaw_a_slip(maxNinstance),                        source=0.0_pReal)
 allocate(plastic_phenopowerlaw_aTolResistance(maxNinstance),                source=0.0_pReal)
 allocate(plastic_phenopowerlaw_aTolShear(maxNinstance),                     source=0.0_pReal)
 allocate(plastic_phenopowerlaw_aTolTwinfrac(maxNinstance),                  source=0.0_pReal)
 allocate(plastic_phenopowerlaw_aTolTransfrac(maxNinstance),                 source=0.0_pReal)
 allocate(plastic_phenopowerlaw_nonSchmidCoeff(lattice_maxNnonSchmid,maxNinstance), &
                                                                                  source=0.0_pReal)
 allocate(plastic_phenopowerlaw_Cnuc(maxNinstance),                          source=0.0_pReal)
 allocate(plastic_phenopowerlaw_Cdwp(maxNinstance),                          source=0.0_pReal)
 allocate(plastic_phenopowerlaw_Cgro(maxNinstance),                          source=0.0_pReal)
 allocate(plastic_phenopowerlaw_deltaG(maxNinstance),                        source=0.0_pReal)

 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo

 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next phase
     phase = phase + 1_pInt                                                                         ! advance phase section counter
     if (phase_plasticity(phase) == PLASTICITY_PHENOPOWERLAW_ID) then
       Nchunks_SlipFamilies  = count(lattice_NslipSystem(:,phase) > 0_pInt)                         ! maximum number of slip families according to lattice type of current phase
       Nchunks_TwinFamilies  = count(lattice_NtwinSystem(:,phase) > 0_pInt)                         ! maximum number of twin families according to lattice type of current phase
       Nchunks_TransFamilies = count(lattice_NtransSystem(:,phase) > 0_pInt)                        ! maximum number of trans families according to lattice type of current phase
       Nchunks_SlipSlip =     maxval(lattice_interactionSlipSlip(:,:,phase))
       Nchunks_SlipTwin =     maxval(lattice_interactionSlipTwin(:,:,phase))
       Nchunks_TwinSlip =     maxval(lattice_interactionTwinSlip(:,:,phase))
       Nchunks_TwinTwin =     maxval(lattice_interactionTwinTwin(:,:,phase))
       Nchunks_nonSchmid =    lattice_NnonSchmid(phase)
       if(allocated(tempPerSlip)) deallocate(tempPerSlip)
       allocate(tempPerSlip(Nchunks_SlipFamilies))
     endif
     cycle                                                                                          ! skip to next line
   endif
   if (phase > 0_pInt ) then; if (phase_plasticity(phase) == PLASTICITY_PHENOPOWERLAW_ID) then      ! one of my phases. Do not short-circuit here (.and. between if-statements), it's not safe in Fortran
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my plasticity is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case ('resistance_slip')
             plastic_phenopowerlaw_Noutput(instance) = plastic_phenopowerlaw_Noutput(instance) + 1_pInt
             plastic_phenopowerlaw_outputID(plastic_phenopowerlaw_Noutput(instance),instance) = resistance_slip_ID
             plastic_phenopowerlaw_output(plastic_phenopowerlaw_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('accumulatedshear_slip','accumulated_shear_slip')
             plastic_phenopowerlaw_Noutput(instance) = plastic_phenopowerlaw_Noutput(instance) + 1_pInt
             plastic_phenopowerlaw_outputID(plastic_phenopowerlaw_Noutput(instance),instance) = accumulatedshear_slip_ID
             plastic_phenopowerlaw_output(plastic_phenopowerlaw_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('shearrate_slip')
             plastic_phenopowerlaw_Noutput(instance) = plastic_phenopowerlaw_Noutput(instance) + 1_pInt
             plastic_phenopowerlaw_outputID(plastic_phenopowerlaw_Noutput(instance),instance) = shearrate_slip_ID
             plastic_phenopowerlaw_output(plastic_phenopowerlaw_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('resolvedstress_slip')
             plastic_phenopowerlaw_Noutput(instance) = plastic_phenopowerlaw_Noutput(instance) + 1_pInt
             plastic_phenopowerlaw_outputID(plastic_phenopowerlaw_Noutput(instance),instance) = resolvedstress_slip_ID
             plastic_phenopowerlaw_output(plastic_phenopowerlaw_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('totalshear')
             plastic_phenopowerlaw_Noutput(instance) = plastic_phenopowerlaw_Noutput(instance) + 1_pInt
             plastic_phenopowerlaw_outputID(plastic_phenopowerlaw_Noutput(instance),instance) = totalshear_ID
             plastic_phenopowerlaw_output(plastic_phenopowerlaw_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('resistance_twin')
             plastic_phenopowerlaw_Noutput(instance) = plastic_phenopowerlaw_Noutput(instance) + 1_pInt
             plastic_phenopowerlaw_outputID(plastic_phenopowerlaw_Noutput(instance),instance) = resistance_twin_ID
             plastic_phenopowerlaw_output(plastic_phenopowerlaw_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('accumulatedshear_twin','accumulated_shear_twin')
             plastic_phenopowerlaw_Noutput(instance) = plastic_phenopowerlaw_Noutput(instance) + 1_pInt
             plastic_phenopowerlaw_outputID(plastic_phenopowerlaw_Noutput(instance),instance) = accumulatedshear_twin_ID
             plastic_phenopowerlaw_output(plastic_phenopowerlaw_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('shearrate_twin')
             plastic_phenopowerlaw_Noutput(instance) = plastic_phenopowerlaw_Noutput(instance) + 1_pInt
             plastic_phenopowerlaw_outputID(plastic_phenopowerlaw_Noutput(instance),instance) = shearrate_twin_ID
             plastic_phenopowerlaw_output(plastic_phenopowerlaw_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('resolvedstress_twin')
             plastic_phenopowerlaw_Noutput(instance) = plastic_phenopowerlaw_Noutput(instance) + 1_pInt
             plastic_phenopowerlaw_outputID(plastic_phenopowerlaw_Noutput(instance),instance) = resolvedstress_twin_ID
             plastic_phenopowerlaw_output(plastic_phenopowerlaw_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('totalvolfrac_twin')
             plastic_phenopowerlaw_Noutput(instance) = plastic_phenopowerlaw_Noutput(instance) + 1_pInt
             plastic_phenopowerlaw_outputID(plastic_phenopowerlaw_Noutput(instance),instance) = totalvolfrac_twin_ID
             plastic_phenopowerlaw_output(plastic_phenopowerlaw_Noutput(instance),instance) = &
                                                           IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case default

         end select
!--------------------------------------------------------------------------------------------------
! parameters depending on number of slip families
       case ('nslip')
         if (chunkPos(1) < Nchunks_SlipFamilies + 1_pInt) &
           call IO_warning(50_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         if (chunkPos(1) > Nchunks_SlipFamilies + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         Nchunks_SlipFamilies = chunkPos(1) - 1_pInt                                                 ! user specified number of (possibly) active slip families (e.g. 6 0 6 --> 3)
         do j = 1_pInt, Nchunks_SlipFamilies
           plastic_phenopowerlaw_Nslip(j,instance) = IO_intValue(line,chunkPos,1_pInt+j)
         enddo
       case ('tausat_slip','tau0_slip')
         tempPerSlip = 0.0_pReal
         do j = 1_pInt, Nchunks_SlipFamilies
           if (plastic_phenopowerlaw_Nslip(j,instance) > 0_pInt) &
             tempPerSlip(j) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         select case(tag)
           case ('tausat_slip')
             plastic_phenopowerlaw_tausat_slip(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('tau0_slip')
             plastic_phenopowerlaw_tau0_slip(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
         end select
!--------------------------------------------------------------------------------------------------
! parameters depending on number of twin families
       case ('ntwin')
         if (chunkPos(1) < Nchunks_TwinFamilies + 1_pInt) &
           call IO_warning(51_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         if (chunkPos(1) > Nchunks_TwinFamilies + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         Nchunks_TwinFamilies = chunkPos(1) - 1_pInt
         do j = 1_pInt, Nchunks_TwinFamilies
             plastic_phenopowerlaw_Ntwin(j,instance) = IO_intValue(line,chunkPos,1_pInt+j)
         enddo
       case ('tau0_twin')
         do j = 1_pInt, Nchunks_TwinFamilies
           if (plastic_phenopowerlaw_Ntwin(j,instance) > 0_pInt) &
             plastic_phenopowerlaw_tau0_twin(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
!--------------------------------------------------------------------------------------------------
! parameters depending on number of transformation families
       case ('ntrans')
         if (chunkPos(1) < Nchunks_TransFamilies + 1_pInt) &
           call IO_warning(51_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         if (chunkPos(1) > Nchunks_TransFamilies + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         Nchunks_TransFamilies = chunkPos(1) - 1_pInt
         do j = 1_pInt, Nchunks_TransFamilies
             plastic_phenopowerlaw_Ntrans(j,instance) = IO_intValue(line,chunkPos,1_pInt+j)
         enddo
!--------------------------------------------------------------------------------------------------
! parameters depending on number of interactions
       case ('interaction_sliptwin')
         if (chunkPos(1) < 1_pInt + Nchunks_SlipTwin) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         do j = 1_pInt, Nchunks_SlipTwin
           plastic_phenopowerlaw_interaction_SlipTwin(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
       case ('interaction_twinslip')
         if (chunkPos(1) < 1_pInt + Nchunks_TwinSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         do j = 1_pInt, Nchunks_TwinSlip
           plastic_phenopowerlaw_interaction_TwinSlip(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
       case ('interaction_twintwin')
         if (chunkPos(1) < 1_pInt + Nchunks_TwinTwin) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         do j = 1_pInt, Nchunks_TwinTwin
           plastic_phenopowerlaw_interaction_TwinTwin(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
       case ('nonschmid_coefficients')
         if (chunkPos(1) < 1_pInt + Nchunks_nonSchmid) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         do j = 1_pInt,Nchunks_nonSchmid
           plastic_phenopowerlaw_nonSchmidCoeff(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
!--------------------------------------------------------------------------------------------------
! parameters independent of number of slip/twin systems
       case ('gdot0_slip')
         plastic_phenopowerlaw_gdot0_slip(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('n_slip')
         plastic_phenopowerlaw_n_slip(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('a_slip', 'w0_slip')
         plastic_phenopowerlaw_a_slip(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('gdot0_twin')
         plastic_phenopowerlaw_gdot0_twin(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('n_twin')
         plastic_phenopowerlaw_n_twin(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('s_pr')
         plastic_phenopowerlaw_spr(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('twin_b')
         plastic_phenopowerlaw_twinB(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('twin_c')
         plastic_phenopowerlaw_twinC(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('twin_d')
         plastic_phenopowerlaw_twinD(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('twin_e')
         plastic_phenopowerlaw_twinE(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('h0_slipslip')
         plastic_phenopowerlaw_h0_SlipSlip(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('h0_twinslip')
         plastic_phenopowerlaw_h0_TwinSlip(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('h0_twintwin')
         plastic_phenopowerlaw_h0_TwinTwin(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('atol_resistance')
         plastic_phenopowerlaw_aTolResistance(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('atol_shear')
         plastic_phenopowerlaw_aTolShear(instance)      = IO_floatValue(line,chunkPos,2_pInt)
       case ('atol_twinfrac')
         plastic_phenopowerlaw_aTolTwinfrac(instance)   = IO_floatValue(line,chunkPos,2_pInt)
       case ('atol_transfrac')
         plastic_phenopowerlaw_aTolTransfrac(instance)  = IO_floatValue(line,chunkPos,2_pInt)
       case ('cnuc')
         plastic_phenopowerlaw_Cnuc(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('cdwp')
         plastic_phenopowerlaw_Cdwp(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('cgro')
         plastic_phenopowerlaw_Cgro(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('deltag')
         plastic_phenopowerlaw_deltaG(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('interaction_slipslip')
         if (chunkPos(1) < 1_pInt + Nchunks_SlipSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         do j = 1_pInt, Nchunks_SlipSlip
           plastic_phenopowerlaw_interaction_SlipSlip(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
       case default

     end select
   endif; endif
 enddo parsingFile

 sanityChecks: do phase = 1_pInt, size(phase_plasticity)
   myPhase: if (phase_plasticity(phase) == PLASTICITY_phenopowerlaw_ID) then
     instance = phase_plasticityInstance(phase)           
     plastic_phenopowerlaw_Nslip(1:lattice_maxNslipFamily,instance) = &
       min(lattice_NslipSystem(1:lattice_maxNslipFamily,phase),&                                    ! limit active slip systems per family to min of available and requested
           plastic_phenopowerlaw_Nslip(1:lattice_maxNslipFamily,instance))
     plastic_phenopowerlaw_Ntwin(1:lattice_maxNtwinFamily,instance) = &
       min(lattice_NtwinSystem(1:lattice_maxNtwinFamily,phase),&                                    ! limit active twin systems per family to min of available and requested
           plastic_phenopowerlaw_Ntwin(:,instance))
     plastic_phenopowerlaw_totalNslip(instance)  = sum(plastic_phenopowerlaw_Nslip(:,instance))           ! how many slip systems altogether
     plastic_phenopowerlaw_totalNtwin(instance)  = sum(plastic_phenopowerlaw_Ntwin(:,instance))           ! how many twin systems altogether
     plastic_phenopowerlaw_totalNtrans(instance) = sum(plastic_phenopowerlaw_Ntrans(:,instance))          ! how many trans systems altogether

     if (any(plastic_phenopowerlaw_tau0_slip(:,instance) < 0.0_pReal .and. &
             plastic_phenopowerlaw_Nslip(:,instance) > 0)) &
       call IO_error(211_pInt,el=instance,ext_msg='tau0_slip ('//PLASTICITY_PHENOPOWERLAW_label//')')
     if (plastic_phenopowerlaw_gdot0_slip(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='gdot0_slip ('//PLASTICITY_PHENOPOWERLAW_label//')')
     if (plastic_phenopowerlaw_n_slip(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='n_slip ('//PLASTICITY_PHENOPOWERLAW_label//')')
     if (any(plastic_phenopowerlaw_tausat_slip(:,instance) <= 0.0_pReal .and. &
             plastic_phenopowerlaw_Nslip(:,instance) > 0)) &
       call IO_error(211_pInt,el=instance,ext_msg='tausat_slip ('//PLASTICITY_PHENOPOWERLAW_label//')')
     if (any(abs(plastic_phenopowerlaw_a_slip(instance)) <= tiny(0.0_pReal) .and. &
             plastic_phenopowerlaw_Nslip(:,instance) > 0)) &
       call IO_error(211_pInt,el=instance,ext_msg='a_slip ('//PLASTICITY_PHENOPOWERLAW_label//')')
     if (any(plastic_phenopowerlaw_tau0_twin(:,instance) < 0.0_pReal .and. &
             plastic_phenopowerlaw_Ntwin(:,instance) > 0)) &
       call IO_error(211_pInt,el=instance,ext_msg='tau0_twin ('//PLASTICITY_PHENOPOWERLAW_label//')')
     if (    plastic_phenopowerlaw_gdot0_twin(instance) <= 0.0_pReal .and. &
         any(plastic_phenopowerlaw_Ntwin(:,instance) > 0)) &
       call IO_error(211_pInt,el=instance,ext_msg='gdot0_twin ('//PLASTICITY_PHENOPOWERLAW_label//')')
     if (    plastic_phenopowerlaw_n_twin(instance) <= 0.0_pReal .and. &
        any(plastic_phenopowerlaw_Ntwin(:,instance) > 0)) &
       call IO_error(211_pInt,el=instance,ext_msg='n_twin ('//PLASTICITY_PHENOPOWERLAW_label//')')
     if (plastic_phenopowerlaw_aTolResistance(instance) <= 0.0_pReal) &
       plastic_phenopowerlaw_aTolResistance(instance) = 1.0_pReal                                       ! default absolute tolerance 1 Pa
     if (plastic_phenopowerlaw_aTolShear(instance) <= 0.0_pReal) &
       plastic_phenopowerlaw_aTolShear(instance) = 1.0e-6_pReal                                         ! default absolute tolerance 1e-6
     if (plastic_phenopowerlaw_aTolTwinfrac(instance) <= 0.0_pReal) &
       plastic_phenopowerlaw_aTolTwinfrac(instance) = 1.0e-6_pReal                                      ! default absolute tolerance 1e-6
     if (plastic_phenopowerlaw_aTolTransfrac(instance) <= 0.0_pReal) &
       plastic_phenopowerlaw_aTolTransfrac(instance) = 1.0e-6_pReal                                     ! default absolute tolerance 1e-6
   endif myPhase
 enddo sanityChecks

!--------------------------------------------------------------------------------------------------
! allocation of variables whose size depends on the total number of active slip systems
 allocate(plastic_phenopowerlaw_hardeningMatrix_SlipSlip(maxval(plastic_phenopowerlaw_totalNslip),&   ! slip resistance from slip activity
                                                              maxval(plastic_phenopowerlaw_totalNslip),&
                                                              maxNinstance), source=0.0_pReal)
 allocate(plastic_phenopowerlaw_hardeningMatrix_SlipTwin(maxval(plastic_phenopowerlaw_totalNslip),&   ! slip resistance from twin activity
                                                              maxval(plastic_phenopowerlaw_totalNtwin),&
                                                              maxNinstance), source=0.0_pReal)
 allocate(plastic_phenopowerlaw_hardeningMatrix_TwinSlip(maxval(plastic_phenopowerlaw_totalNtwin),&   ! twin resistance from slip activity
                                                              maxval(plastic_phenopowerlaw_totalNslip),&
                                                              maxNinstance), source=0.0_pReal)
 allocate(plastic_phenopowerlaw_hardeningMatrix_TwinTwin(maxval(plastic_phenopowerlaw_totalNtwin),&   ! twin resistance from twin activity
                                                              maxval(plastic_phenopowerlaw_totalNtwin),&
                                                              maxNinstance), source=0.0_pReal)

 initializeInstances: do phase = 1_pInt, size(phase_plasticity)                                     ! loop through all phases in material.config
   myPhase2: if (phase_plasticity(phase) == PLASTICITY_phenopowerlaw_ID) then                       ! only consider my phase
     NipcMyPhase = count(material_phase == phase)                                                   ! number of IPCs containing my phase
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my phase

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,plastic_phenopowerlaw_Noutput(instance)
       select case(plastic_phenopowerlaw_outputID(o,instance))
         case(resistance_slip_ID, &
              shearrate_slip_ID, &
              accumulatedshear_slip_ID, &
              resolvedstress_slip_ID &
              )
           mySize = plastic_phenopowerlaw_totalNslip(instance)
         case(resistance_twin_ID, &
              shearrate_twin_ID, &
              accumulatedshear_twin_ID, &
              resolvedstress_twin_ID &
              )
           mySize = plastic_phenopowerlaw_totalNtwin(instance)
         case(totalshear_ID, &
              totalvolfrac_twin_ID &
              )
           mySize = 1_pInt
         case default
       end select
  
       outputFound: if (mySize > 0_pInt) then
         plastic_phenopowerlaw_sizePostResult(o,instance) = mySize
         plastic_phenopowerlaw_sizePostResults(instance)  = plastic_phenopowerlaw_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop 
!--------------------------------------------------------------------------------------------------
! allocate state arrays
     sizeState = plastic_phenopowerlaw_totalNslip(instance) &                         ! s_slip
               + plastic_phenopowerlaw_totalNtwin(instance) &                         ! s_twin
               + 2_pInt &                                                             ! sum(gamma) + sum(f)
               + plastic_phenopowerlaw_totalNslip(instance) &                         ! accshear_slip
               + plastic_phenopowerlaw_totalNtwin(instance)                           ! accshear_twin
               
     sizeDotState = sizeState
     sizeDeltaState = 0_pInt
     plasticState(phase)%sizeState = sizeState
     plasticState(phase)%sizeDotState = sizeDotState
     plasticState(phase)%sizeDeltaState = sizeDeltaState
     plasticState(phase)%sizePostResults = plastic_phenopowerlaw_sizePostResults(instance)
     plasticState(phase)%nSlip =plastic_phenopowerlaw_totalNslip(instance)
     plasticState(phase)%nTwin =plastic_phenopowerlaw_totalNtwin(instance)
     plasticState(phase)%nTrans=plastic_phenopowerlaw_totalNtrans(instance)
     allocate(plasticState(phase)%aTolState          (   sizeState),             source=0.0_pReal)
     allocate(plasticState(phase)%state0             (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%partionedState0    (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%subState0          (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%state              (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%state_backup       (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%dotState           (sizeDotState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%deltaState       (sizeDeltaState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%dotState_backup    (sizeDotState,NipcMyPhase), source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(phase)%previousDotState (sizeDotState,NipcMyPhase),source=0.0_pReal)
       allocate(plasticState(phase)%previousDotState2(sizeDotState,NipcMyPhase),source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(phase)%RK4dotState      (sizeDotState,NipcMyPhase), source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(phase)%RKCK45dotState (6,sizeDotState,NipcMyPhase), source=0.0_pReal)
     
     offset_slip = plasticState(phase)%nSlip+plasticState(phase)%nTwin+2_pInt
     plasticState(phase)%slipRate => &
       plasticState(phase)%dotState(offset_slip+1:offset_slip+plasticState(phase)%nSlip,1:NipcMyPhase)
     plasticState(phase)%accumulatedSlip => &
       plasticState(phase)%state(offset_slip+1:offset_slip+plasticState(phase)%nSlip,1:NipcMyPhase)
       
     do f = 1_pInt,lattice_maxNslipFamily                                                                    ! >>> interaction slip -- X
       index_myFamily = sum(plastic_phenopowerlaw_Nslip(1:f-1_pInt,instance))
       do j = 1_pInt,plastic_phenopowerlaw_Nslip(f,instance)                                            ! loop over (active) systems in my family (slip)
         do o = 1_pInt,lattice_maxNslipFamily
           index_otherFamily = sum(plastic_phenopowerlaw_Nslip(1:o-1_pInt,instance))
           do k = 1_pInt,plastic_phenopowerlaw_Nslip(o,instance)                                        ! loop over (active) systems in other family (slip)
             plastic_phenopowerlaw_hardeningMatrix_SlipSlip(index_myFamily+j,index_otherFamily+k,instance) = &
                 plastic_phenopowerlaw_interaction_SlipSlip(lattice_interactionSlipSlip( &
                                                                   sum(lattice_NslipSystem(1:f-1,phase))+j, &
                                                                   sum(lattice_NslipSystem(1:o-1,phase))+k, &
                                                                   phase), instance )
         enddo; enddo
  
         do o = 1_pInt,lattice_maxNtwinFamily
           index_otherFamily = sum(plastic_phenopowerlaw_Ntwin(1:o-1_pInt,instance))
           do k = 1_pInt,plastic_phenopowerlaw_Ntwin(o,instance)                                        ! loop over (active) systems in other family (twin)
             plastic_phenopowerlaw_hardeningMatrix_SlipTwin(index_myFamily+j,index_otherFamily+k,instance) = &
                 plastic_phenopowerlaw_interaction_SlipTwin(lattice_interactionSlipTwin( &
                                                                   sum(lattice_NslipSystem(1:f-1_pInt,phase))+j, &
                                                                   sum(lattice_NtwinSystem(1:o-1_pInt,phase))+k, &
                                                                   phase), instance )
         enddo; enddo
  
     enddo; enddo
  
     do f = 1_pInt,lattice_maxNtwinFamily                                                                    ! >>> interaction twin -- X
       index_myFamily = sum(plastic_phenopowerlaw_Ntwin(1:f-1_pInt,instance))
       do j = 1_pInt,plastic_phenopowerlaw_Ntwin(f,instance)                                            ! loop over (active) systems in my family (twin)
  
         do o = 1_pInt,lattice_maxNslipFamily
           index_otherFamily = sum(plastic_phenopowerlaw_Nslip(1:o-1_pInt,instance))
           do k = 1_pInt,plastic_phenopowerlaw_Nslip(o,instance)                                        ! loop over (active) systems in other family (slip)
             plastic_phenopowerlaw_hardeningMatrix_TwinSlip(index_myFamily+j,index_otherFamily+k,instance) = &
                 plastic_phenopowerlaw_interaction_TwinSlip(lattice_interactionTwinSlip( &
                                                                   sum(lattice_NtwinSystem(1:f-1_pInt,phase))+j, &
                                                                   sum(lattice_NslipSystem(1:o-1_pInt,phase))+k, &
                                                                   phase), instance )
         enddo; enddo
  
         do o = 1_pInt,lattice_maxNtwinFamily
           index_otherFamily = sum(plastic_phenopowerlaw_Ntwin(1:o-1_pInt,instance))
           do k = 1_pInt,plastic_phenopowerlaw_Ntwin(o,instance)                                        ! loop over (active) systems in other family (twin)
             plastic_phenopowerlaw_hardeningMatrix_TwinTwin(index_myFamily+j,index_otherFamily+k,instance) = &
                 plastic_phenopowerlaw_interaction_TwinTwin(lattice_interactionTwinTwin( &
                                                                   sum(lattice_NtwinSystem(1:f-1_pInt,phase))+j, &
                                                                   sum(lattice_NtwinSystem(1:o-1_pInt,phase))+k, &
                                                                   phase), instance )
         enddo; enddo
  
     enddo; enddo

     call plastic_phenopowerlaw_stateInit(phase,instance)
     call plastic_phenopowerlaw_aTolState(phase,instance)
   endif myPhase2 
 enddo initializeInstances

end subroutine plastic_phenopowerlaw_init


!--------------------------------------------------------------------------------------------------
!> @brief sets the initial microstructural state for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
subroutine plastic_phenopowerlaw_stateInit(ph,instance)
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily
 use material, only: &
   plasticState
 
 implicit none
 integer(pInt), intent(in) :: &
   instance, &                                                                                        !< number specifying the instance of the plasticity
   ph 
 integer(pInt) :: &
   i
 real(pReal),   dimension(plasticState(ph)%sizeState) :: &
   tempState

 tempState = 0.0_pReal
 do i = 1_pInt,lattice_maxNslipFamily
   tempState(1+sum(plastic_phenopowerlaw_Nslip(1:i-1,instance)) : &
               sum(plastic_phenopowerlaw_Nslip(1:i  ,instance))) = &
     plastic_phenopowerlaw_tau0_slip(i,instance)
 enddo

 do i = 1_pInt,lattice_maxNtwinFamily
   tempState(1+sum(plastic_phenopowerlaw_Nslip(:,instance))+&
               sum(plastic_phenopowerlaw_Ntwin(1:i-1,instance)) : &
               sum(plastic_phenopowerlaw_Nslip(:,instance))+&
               sum(plastic_phenopowerlaw_Ntwin(1:i  ,instance))) = &
     plastic_phenopowerlaw_tau0_twin(i,instance)
 enddo

 plasticState(ph)%state0(:,:) = spread(tempState, &                                                 ! spread single tempstate array
                                       2, &                                                         ! along dimension 2
                                       size(plasticState(ph)%state0(1,:)))                          ! number of copies (number of IPCs)

end subroutine plastic_phenopowerlaw_stateInit


!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
subroutine plastic_phenopowerlaw_aTolState(ph,instance)
  use material, only: &
   plasticState

 implicit none
 integer(pInt), intent(in) :: & 
   instance, &                                                              !< number specifying the instance of the plasticity
   ph

 plasticState(ph)%aTolState(1:plastic_phenopowerlaw_totalNslip(instance)+ &
                              plastic_phenopowerlaw_totalNtwin(instance)) = &
                                              plastic_phenopowerlaw_aTolResistance(instance)
 plasticState(ph)%aTolState(1+plastic_phenopowerlaw_totalNslip(instance)+ &
                              plastic_phenopowerlaw_totalNtwin(instance)) = &
                                              plastic_phenopowerlaw_aTolShear(instance)
 plasticState(ph)%aTolState(2+plastic_phenopowerlaw_totalNslip(instance)+ &
                              plastic_phenopowerlaw_totalNtwin(instance)) = &
                                             plastic_phenopowerlaw_aTolTwinFrac(instance)
 plasticState(ph)%aTolState(3+plastic_phenopowerlaw_totalNslip(instance)+ &
                              plastic_phenopowerlaw_totalNtwin(instance): &
                            2+2*(plastic_phenopowerlaw_totalNslip(instance)+ &
                                 plastic_phenopowerlaw_totalNtwin(instance))) = &
                                             plastic_phenopowerlaw_aTolShear(instance)

end subroutine plastic_phenopowerlaw_aTolState


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_phenopowerlaw_LpAndItsTangent(Lp,dLp_dTstar99,Tstar_v,ipc,ip,el)
 use math, only: &
   math_Plain3333to99, &
   math_Mandel6to33
 use lattice, only: &
   lattice_Sslip, &
   lattice_Sslip_v, &
   lattice_Stwin, &
   lattice_Stwin_v, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_NnonSchmid
 use material, only: &
   plasticState, &
   mappingConstitutive, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3), intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(9,9), intent(out) :: &
   dLp_dTstar99                                                                                     !< derivative of Lp with respect to 2nd Piola Kirchhoff stress

 integer(pInt),               intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(6),   intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation

 integer(pInt) :: &
   instance, & 
   nSlip, &
   nTwin,index_Gamma,index_F,index_myFamily, &
   f,i,j,k,l,m,n, &
   of, &
   ph
 real(pReal) :: &
   tau_slip_pos,tau_slip_neg, &
   gdot_slip_pos,gdot_slip_neg, &
   dgdot_dtauslip_pos,dgdot_dtauslip_neg, &
   gdot_twin,dgdot_dtautwin,tau_twin
 real(pReal), dimension(3,3,3,3) :: &
   dLp_dTstar3333                                                                                   !< derivative of Lp with respect to Tstar as 4th order tensor
 real(pReal), dimension(3,3,2) :: &
   nonSchmid_tensor
 
 of = mappingConstitutive(1,ipc,ip,el)
 ph = mappingConstitutive(2,ipc,ip,el)
 instance = phase_plasticityInstance(ph)
 nSlip = plastic_phenopowerlaw_totalNslip(instance)
 nTwin = plastic_phenopowerlaw_totalNtwin(instance)
 index_Gamma = nSlip + nTwin + 1_pInt
 index_F     = nSlip + nTwin + 2_pInt

 Lp = 0.0_pReal
 dLp_dTstar3333 = 0.0_pReal
 dLp_dTstar99 = 0.0_pReal

!--------------------------------------------------------------------------------------------------
! Slip part
 j = 0_pInt
 slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                          ! at which index starts my family
   slipSystems: do i = 1_pInt,plastic_phenopowerlaw_Nslip(f,instance)
     j = j+1_pInt
     
     ! Calculation of Lp
     tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
     tau_slip_neg  = tau_slip_pos
     nonSchmid_tensor(1:3,1:3,1) = lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph)
     nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,1)
     do k = 1,lattice_NnonSchmid(ph) 
       tau_slip_pos = tau_slip_pos + plastic_phenopowerlaw_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,ph))
       tau_slip_neg = tau_slip_neg + plastic_phenopowerlaw_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
       nonSchmid_tensor(1:3,1:3,1) = nonSchmid_tensor(1:3,1:3,1) + plastic_phenopowerlaw_nonSchmidCoeff(k,instance)*&
                                           lattice_Sslip(1:3,1:3,2*k,index_myFamily+i,ph)
       nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,2) + plastic_phenopowerlaw_nonSchmidCoeff(k,instance)*&
                                           lattice_Sslip(1:3,1:3,2*k+1,index_myFamily+i,ph)
     enddo
     gdot_slip_pos = 0.5_pReal*plastic_phenopowerlaw_gdot0_slip(instance)* &
                    ((abs(tau_slip_pos)/(plasticState(ph)%state(j,of))) &
                    **plastic_phenopowerlaw_n_slip(instance))*sign(1.0_pReal,tau_slip_pos)

     gdot_slip_neg = 0.5_pReal*plastic_phenopowerlaw_gdot0_slip(instance)* &
                    ((abs(tau_slip_neg)/(plasticState(ph)%state(j,of))) &
                    **plastic_phenopowerlaw_n_slip(instance))*sign(1.0_pReal,tau_slip_neg)
                                                                    
     Lp = Lp + (1.0_pReal-plasticState(ph)%state(index_F,of))*&                                  ! 1-F
               (gdot_slip_pos+gdot_slip_neg)*lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph)

     ! Calculation of the tangent of Lp
     if (abs(gdot_slip_pos) > tiny(0.0_pReal)) then
       dgdot_dtauslip_pos = gdot_slip_pos*plastic_phenopowerlaw_n_slip(instance)/tau_slip_pos
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtauslip_pos*lattice_Sslip(k,l,1,index_myFamily+i,ph)* &
                                                     nonSchmid_tensor(m,n,1)
     endif
     
     if (abs(gdot_slip_neg) > tiny(0.0_pReal)) then
       dgdot_dtauslip_neg = gdot_slip_neg*plastic_phenopowerlaw_n_slip(instance)/tau_slip_neg
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtauslip_neg*lattice_Sslip(k,l,1,index_myFamily+i,ph)* &
                                                     nonSchmid_tensor(m,n,2)
     endif
   enddo slipSystems
 enddo slipFamilies

!--------------------------------------------------------------------------------------------------
! Twinning part
 j = 0_pInt
 twinFamilies: do f = 1_pInt,lattice_maxNtwinFamily
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                      ! at which index starts my family
   twinSystems: do i = 1_pInt,plastic_phenopowerlaw_Ntwin(f,instance)
     j = j+1_pInt

     ! Calculation of Lp
     tau_twin  = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,ph)) 
     gdot_twin = (1.0_pReal-plasticState(ph)%state(index_F,of))*&                                                  ! 1-F
                    plastic_phenopowerlaw_gdot0_twin(instance)*&
                    (abs(tau_twin)/plasticState(ph)%state(nSlip+j,of))**&
                    plastic_phenopowerlaw_n_twin(instance)*max(0.0_pReal,sign(1.0_pReal,tau_twin))               
     Lp = Lp + gdot_twin*lattice_Stwin(1:3,1:3,index_myFamily+i,ph)

     ! Calculation of the tangent of Lp
     if (abs(gdot_twin) > tiny(0.0_pReal)) then
       dgdot_dtautwin = gdot_twin*plastic_phenopowerlaw_n_twin(instance)/tau_twin
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtautwin*lattice_Stwin(k,l,index_myFamily+i,ph)* &
                                                  lattice_Stwin(m,n,index_myFamily+i,ph)
     endif
   enddo twinSystems
 enddo twinFamilies

 dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)


end subroutine plastic_phenopowerlaw_LpAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_phenopowerlaw_dotState(Tstar_v,ipc,ip,el)
 use lattice, only: &
   lattice_Sslip_v, &
   lattice_Stwin_v, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_shearTwin, &
   lattice_NnonSchmid   
 use material, only: &
   material_phase, &
   mappingConstitutive, &
   plasticState, &
   phase_plasticityInstance
 
 implicit none
 real(pReal), dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element                                                                                    !< microstructure state

 integer(pInt) :: &
   instance,ph, &
   nSlip,nTwin, &
   f,i,j,k, &
   index_Gamma,index_F,index_myFamily, &
   offset_accshear_slip,offset_accshear_twin, &
   of
 real(pReal) :: &
   c_SlipSlip,c_TwinSlip,c_TwinTwin, &
   ssat_offset, &
   tau_slip_pos,tau_slip_neg,tau_twin

 real(pReal), dimension(plastic_phenopowerlaw_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip,left_SlipSlip,left_SlipTwin,right_SlipSlip,right_TwinSlip
 real(pReal), dimension(plastic_phenopowerlaw_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_twin,left_TwinSlip,left_TwinTwin,right_SlipTwin,right_TwinTwin
 
 of = mappingConstitutive(1,ipc,ip,el)
 ph = mappingConstitutive(2,ipc,ip,el)
 instance = phase_plasticityInstance(ph)
 
 nSlip = plastic_phenopowerlaw_totalNslip(instance)
 nTwin = plastic_phenopowerlaw_totalNtwin(instance)

 index_Gamma = nSlip + nTwin + 1_pInt
 index_F     = nSlip + nTwin + 2_pInt
 offset_accshear_slip = nSlip + nTwin + 2_pInt
 offset_accshear_twin = nSlip + nTwin + 2_pInt + nSlip
 plasticState(ph)%dotState(:,of) = 0.0_pReal
 

!--------------------------------------------------------------------------------------------------
! system-independent (nonlinear) prefactors to M_Xx (X influenced by x) matrices
 c_SlipSlip = plastic_phenopowerlaw_h0_SlipSlip(instance)*&
              (1.0_pReal + plastic_phenopowerlaw_twinC(instance)*plasticState(ph)%state(index_F,of)**&
                                                           plastic_phenopowerlaw_twinB(instance))
 c_TwinSlip = plastic_phenopowerlaw_h0_TwinSlip(instance)*&
              plasticState(ph)%state(index_Gamma,of)**plastic_phenopowerlaw_twinE(instance)
 c_TwinTwin = plastic_phenopowerlaw_h0_TwinTwin(instance)*&
              plasticState(ph)%state(index_F,of)**plastic_phenopowerlaw_twinD(instance)

!--------------------------------------------------------------------------------------------------
!  calculate left and right vectors and calculate dot gammas
 ssat_offset = plastic_phenopowerlaw_spr(instance)*sqrt(plasticState(ph)%state(index_F,of))
 j = 0_pInt
 slipFamilies1: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                         ! at which index starts my family
   slipSystems1: do i = 1_pInt,plastic_phenopowerlaw_Nslip(f,instance)
     j = j+1_pInt
     left_SlipSlip(j) = 1.0_pReal                                                                   ! no system-dependent left part
     left_SlipTwin(j) = 1.0_pReal                                                                   ! no system-dependent left part
     right_SlipSlip(j) = abs(1.0_pReal-plasticState(ph)%state(j,of) / &
                                    (plastic_phenopowerlaw_tausat_slip(f,instance)+ssat_offset)) &
                         **plastic_phenopowerlaw_a_slip(instance)&
                         *sign(1.0_pReal,1.0_pReal-plasticState(ph)%state(j,of) / &
                                    (plastic_phenopowerlaw_tausat_slip(f,instance)+ssat_offset))
     right_TwinSlip(j) = 1.0_pReal                                                                  ! no system-dependent part
     
!--------------------------------------------------------------------------------------------------
! Calculation of dot gamma 
     tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
     tau_slip_neg  = tau_slip_pos
     nonSchmidSystems: do k = 1,lattice_NnonSchmid(ph) 
       tau_slip_pos = tau_slip_pos + plastic_phenopowerlaw_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,  index_myFamily+i,ph))
       tau_slip_neg = tau_slip_neg + plastic_phenopowerlaw_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
     enddo nonSchmidSystems
     gdot_slip(j) = plastic_phenopowerlaw_gdot0_slip(instance)*0.5_pReal* &
                  ((abs(tau_slip_pos)/(plasticState(ph)%state(j,of)))**plastic_phenopowerlaw_n_slip(instance) &
                  +(abs(tau_slip_neg)/(plasticState(ph)%state(j,of)))**plastic_phenopowerlaw_n_slip(instance))&
                  *sign(1.0_pReal,tau_slip_pos) 
   enddo slipSystems1
 enddo slipFamilies1


 j = 0_pInt
 twinFamilies1: do f = 1_pInt,lattice_maxNtwinFamily
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                         ! at which index starts my family
   twinSystems1: do i = 1_pInt,plastic_phenopowerlaw_Ntwin(f,instance)
     j = j+1_pInt
     left_TwinSlip(j)  = 1.0_pReal                                                                  ! no system-dependent right part
     left_TwinTwin(j)  = 1.0_pReal                                                                  ! no system-dependent right part
     right_SlipTwin(j) = 1.0_pReal                                                                  ! no system-dependent right part
     right_TwinTwin(j) = 1.0_pReal                                                                  ! no system-dependent right part

!--------------------------------------------------------------------------------------------------
! Calculation of dot vol frac
     tau_twin  = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,ph)) 
     gdot_twin(j) = (1.0_pReal-plasticState(ph)%state(index_F,of))*&                                       ! 1-F
                    plastic_phenopowerlaw_gdot0_twin(instance)*&
                    (abs(tau_twin)/plasticState(ph)%state(nslip+j,of))**&
                    plastic_phenopowerlaw_n_twin(instance)*max(0.0_pReal,sign(1.0_pReal,tau_twin))
    enddo twinSystems1
  enddo twinFamilies1

!--------------------------------------------------------------------------------------------------
! calculate the overall hardening based on above
 j = 0_pInt
 slipFamilies2: do f = 1_pInt,lattice_maxNslipFamily
   slipSystems2: do i = 1_pInt,plastic_phenopowerlaw_Nslip(f,instance)
     j = j+1_pInt
     plasticState(ph)%dotState(j,of) = &                                                            ! evolution of slip resistance j
       c_SlipSlip * left_SlipSlip(j) * &
       dot_product(plastic_phenopowerlaw_hardeningMatrix_SlipSlip(j,1:nSlip,instance), &
                   right_SlipSlip*abs(gdot_slip)) + &                                               ! dot gamma_slip modulated by right-side slip factor
       dot_product(plastic_phenopowerlaw_hardeningMatrix_SlipTwin(j,1:nTwin,instance), &
                   right_SlipTwin*gdot_twin)                                                        ! dot gamma_twin modulated by right-side twin factor
     plasticState(ph)%dotState(index_Gamma,of) = plasticState(ph)%dotState(index_Gamma,of) + &
                                                        abs(gdot_slip(j))
     plasticState(ph)%dotState(offset_accshear_slip+j,of) = abs(gdot_slip(j))
   enddo slipSystems2
 enddo slipFamilies2

 j = 0_pInt
 twinFamilies2: do f = 1_pInt,lattice_maxNtwinFamily
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                         ! at which index starts my family
   twinSystems2: do i = 1_pInt,plastic_phenopowerlaw_Ntwin(f,instance)
     j = j+1_pInt
     plasticState(ph)%dotState(j+nSlip,of) = &                                                      ! evolution of twin resistance j
       c_TwinSlip * left_TwinSlip(j) * &
       dot_product(plastic_phenopowerlaw_hardeningMatrix_TwinSlip(j,1:nSlip,instance), &
                   right_TwinSlip*abs(gdot_slip)) + &                                               ! dot gamma_slip modulated by right-side slip factor
       c_TwinTwin * left_TwinTwin(j) * &
       dot_product(plastic_phenopowerlaw_hardeningMatrix_TwinTwin(j,1:nTwin,instance), &
                   right_TwinTwin*gdot_twin)                                                        ! dot gamma_twin modulated by right-side twin factor
     if (plasticState(ph)%state(index_F,of) < 0.98_pReal) &                                         ! ensure twin volume fractions stays below 1.0
       plasticState(ph)%dotState(index_F,of) = plasticState(ph)%dotState(index_F,of) + &
                                                      gdot_twin(j)/lattice_shearTwin(index_myFamily+i,ph)
     plasticState(ph)%dotState(offset_accshear_twin+j,of) = abs(gdot_twin(j))
   enddo twinSystems2
 enddo twinFamilies2

 
end subroutine plastic_phenopowerlaw_dotState

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_phenopowerlaw_postResults(Tstar_v,ipc,ip,el)
 use material, only: &
   material_phase, &
   plasticState, &
   mappingConstitutive, &
   phase_plasticityInstance
 use lattice, only: &
   lattice_Sslip_v, &
   lattice_Stwin_v, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_NnonSchmid 

 implicit none
 real(pReal), dimension(6), intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),             intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element                                                                                        !< microstructure state

 real(pReal), dimension(plastic_phenopowerlaw_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   plastic_phenopowerlaw_postResults

 integer(pInt) :: &
   instance,ph, of, &
   nSlip,nTwin, &
   o,f,i,c,j,k, &
   index_Gamma,index_F,index_accshear_slip,index_accshear_twin,index_myFamily 
 real(pReal) :: &
   tau_slip_pos,tau_slip_neg,tau

 of = mappingConstitutive(1,ipc,ip,el)
 ph = mappingConstitutive(2,ipc,ip,el)
 instance = phase_plasticityInstance(ph)

 nSlip = plastic_phenopowerlaw_totalNslip(instance)
 nTwin = plastic_phenopowerlaw_totalNtwin(instance)

 index_Gamma = nSlip + nTwin + 1_pInt
 index_F     = nSlip + nTwin + 2_pInt
 index_accshear_slip = nSlip + nTwin + 3_pInt
 index_accshear_twin = nSlip + nTwin + 3_pInt + nSlip

 plastic_phenopowerlaw_postResults = 0.0_pReal
 c = 0_pInt

 outputsLoop: do o = 1_pInt,plastic_phenopowerlaw_Noutput(instance)
   select case(plastic_phenopowerlaw_outputID(o,instance))
     case (resistance_slip_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt:c+nSlip) = plasticState(ph)%state(1:nSlip,of)
       c = c + nSlip

     case (accumulatedshear_slip_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt:c+nSlip) = plasticState(ph)%state(index_accshear_slip:&
                                                                        index_accshear_slip+nSlip-1_pInt,of)
       c = c + nSlip

     case (shearrate_slip_ID)
       j = 0_pInt
       slipFamilies1: do f = 1_pInt,lattice_maxNslipFamily
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         slipSystems1: do i = 1_pInt,plastic_phenopowerlaw_Nslip(f,instance)
           j = j + 1_pInt
           tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
           tau_slip_neg  = tau_slip_pos
           do k = 1,lattice_NnonSchmid(ph) 
             tau_slip_pos = tau_slip_pos + plastic_phenopowerlaw_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,ph))
             tau_slip_neg = tau_slip_neg + plastic_phenopowerlaw_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
           enddo
           plastic_phenopowerlaw_postResults(c+j) = plastic_phenopowerlaw_gdot0_slip(instance)*0.5_pReal* &
                    ((abs(tau_slip_pos)/plasticState(ph)%state(j,of))**plastic_phenopowerlaw_n_slip(instance) &
                    +(abs(tau_slip_neg)/plasticState(ph)%state(j,of))**plastic_phenopowerlaw_n_slip(instance))&
                    *sign(1.0_pReal,tau_slip_pos)

         enddo slipSystems1
       enddo slipFamilies1
       c = c + nSlip

     case (resolvedstress_slip_ID)
       j = 0_pInt
       slipFamilies2: do f = 1_pInt,lattice_maxNslipFamily
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         slipSystems2: do i = 1_pInt,plastic_phenopowerlaw_Nslip(f,instance)
           j = j + 1_pInt
           plastic_phenopowerlaw_postResults(c+j) = &
                             dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
         enddo slipSystems2
       enddo slipFamilies2
       c = c + nSlip

     case (totalshear_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt) = &
                             plasticState(ph)%state(index_Gamma,of)
       c = c + 1_pInt

     case (resistance_twin_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt:c+nTwin) = &
                             plasticState(ph)%state(1_pInt+nSlip:nTwin+nSlip-1_pInt,of)
       c = c + nTwin

     case (accumulatedshear_twin_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt:c+nTwin) = &
                             plasticState(ph)%state(index_accshear_twin:index_accshear_twin+nTwin-1_pInt,of)
       c = c + nTwin
     case (shearrate_twin_ID)
       j = 0_pInt
       twinFamilies1: do f = 1_pInt,lattice_maxNtwinFamily
         index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         twinSystems1: do i = 1_pInt,plastic_phenopowerlaw_Ntwin(f,instance)
           j = j + 1_pInt
           tau = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,ph))
           plastic_phenopowerlaw_postResults(c+j) = (1.0_pReal-plasticState(ph)%state(index_F,of))*&  ! 1-F
                                                         plastic_phenopowerlaw_gdot0_twin(instance)*&
                                                         (abs(tau)/plasticState(ph)%state(j+nSlip,of))**&
                                           plastic_phenopowerlaw_n_twin(instance)*max(0.0_pReal,sign(1.0_pReal,tau))
         enddo twinSystems1
       enddo twinFamilies1
       c = c + nTwin

     case (resolvedstress_twin_ID)
       j = 0_pInt
       twinFamilies2: do f = 1_pInt,lattice_maxNtwinFamily
         index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         twinSystems2: do i = 1_pInt,plastic_phenopowerlaw_Ntwin(f,instance)
           j = j + 1_pInt
           plastic_phenopowerlaw_postResults(c+j) = &
                             dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,ph))
         enddo twinSystems2
       enddo twinFamilies2
       c = c + nTwin

     case (totalvolfrac_twin_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt) = plasticState(ph)%state(index_F,of)
       c = c + 1_pInt

   end select
 enddo outputsLoop

end function plastic_phenopowerlaw_postResults

end module plastic_phenopowerlaw
