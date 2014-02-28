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
!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for phenomenological crystal plasticity formulation using a powerlaw 
!! fitting
!--------------------------------------------------------------------------------------------------
module constitutive_phenopowerlaw
 use prec, only: &
   pReal,&
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,         public, protected :: &
   constitutive_phenopowerlaw_sizeDotState, &
   constitutive_phenopowerlaw_sizeState, &
   constitutive_phenopowerlaw_sizePostResults, &                                                    !< cumulative size of post results
   constitutive_phenopowerlaw_structure

 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   constitutive_phenopowerlaw_sizePostResult                                                        !< size of each post result output

 character(len=64),                   dimension(:,:),   allocatable, target, public :: & 
   constitutive_phenopowerlaw_output                                                                !< name of each post result output

 integer(pInt),                       dimension(:),     allocatable,         private :: &
   constitutive_phenopowerlaw_Noutput, &                                                            !< number of outputs per instance of this constitution 
   constitutive_phenopowerlaw_totalNslip, &                                                         !< no. of slip system used in simulation
   constitutive_phenopowerlaw_totalNtwin                                                            !< no. of twin system used in simulation

 integer(pInt),                       dimension(:,:),   allocatable,         private :: &
   constitutive_phenopowerlaw_Nslip, &                                                              !< active number of slip systems per family (input parameter, per family)
   constitutive_phenopowerlaw_Ntwin                                                                 !< active number of twin systems per family (input parameter, per family)

 real(pReal),                         dimension(:),     allocatable,         private :: &
   constitutive_phenopowerlaw_gdot0_slip, &                                                         !< reference shear strain rate for slip (input parameter)
   constitutive_phenopowerlaw_gdot0_twin, &                                                         !< reference shear strain rate for twin (input parameter)
   constitutive_phenopowerlaw_n_slip, &                                                             !< stress exponent for slip (input parameter)
   constitutive_phenopowerlaw_n_twin, &                                                             !< stress exponent for twin (input parameter)

   constitutive_phenopowerlaw_spr, &                                                                !< push-up factor for slip saturation due to twinning
   constitutive_phenopowerlaw_twinB, &
   constitutive_phenopowerlaw_twinC, &
   constitutive_phenopowerlaw_twinD, &
   constitutive_phenopowerlaw_twinE, &
   constitutive_phenopowerlaw_h0_SlipSlip, &                                                        !< reference hardening slip - slip (input parameter)
   constitutive_phenopowerlaw_h0_SlipTwin, &                                                        !< reference hardening slip - twin (input parameter, no effect at the moment)
   constitutive_phenopowerlaw_h0_TwinSlip, &                                                        !< reference hardening twin - slip (input parameter)
   constitutive_phenopowerlaw_h0_TwinTwin, &                                                        !< reference hardening twin - twin (input parameter)
   constitutive_phenopowerlaw_a_slip, &
   constitutive_phenopowerlaw_aTolResistance, &
   constitutive_phenopowerlaw_aTolShear, &
   constitutive_phenopowerlaw_aTolTwinfrac

 real(pReal),                         dimension(:,:),   allocatable,          private :: &
   constitutive_phenopowerlaw_tau0_slip, &                                                          !< initial critical shear stress for slip (input parameter, per family)
   constitutive_phenopowerlaw_tau0_twin, &                                                          !< initial critical shear stress for twin (input parameter, per family)
   constitutive_phenopowerlaw_tausat_slip, &                                                        !< maximum critical shear stress for slip (input parameter, per family)
   constitutive_phenopowerlaw_nonSchmidCoeff, &

   constitutive_phenopowerlaw_interaction_SlipSlip, &                                               !< interaction factors slip - slip (input parameter)
   constitutive_phenopowerlaw_interaction_SlipTwin, &                                               !< interaction factors slip - twin (input parameter)
   constitutive_phenopowerlaw_interaction_TwinSlip, &                                               !< interaction factors twin - slip (input parameter)
   constitutive_phenopowerlaw_interaction_TwinTwin                                                  !< interaction factors twin - twin (input parameter)

 real(pReal),                         dimension(:,:,:), allocatable,          private :: &
   constitutive_phenopowerlaw_hardeningMatrix_SlipSlip, &
   constitutive_phenopowerlaw_hardeningMatrix_SlipTwin, &
   constitutive_phenopowerlaw_hardeningMatrix_TwinSlip, &
   constitutive_phenopowerlaw_hardeningMatrix_TwinTwin, &
   constitutive_phenopowerlaw_Cslip_66
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
                 totalvolfrac_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),   allocatable,          private :: & 
   constitutive_phenopowerlaw_outputID                                                              !< ID of each post result output
 
 public :: &
   constitutive_phenopowerlaw_init, &
   constitutive_phenopowerlaw_stateInit, &
   constitutive_phenopowerlaw_aTolState, &
   constitutive_phenopowerlaw_homogenizedC, &
   constitutive_phenopowerlaw_LpAndItsTangent, &
   constitutive_phenopowerlaw_dotState, &
   constitutive_phenopowerlaw_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine constitutive_phenopowerlaw_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use prec, only: &
   tol_math_check
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
   homogenization_maxNgrains, &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_PHENOPOWERLAW_label, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   MATERIAL_partPhase
 use lattice

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), parameter :: MAXNCHUNKS = LATTICE_maxNinteraction + 1_pInt
 integer(pInt), dimension(1_pInt+2_pInt*MAXNCHUNKS) :: positions
 integer(pInt), dimension(7) :: configNchunks
 integer(pInt) :: &
   maxNinstance, &
   instance,j,k, f,o, &
   Nchunks_SlipSlip, Nchunks_SlipTwin, Nchunks_TwinSlip, Nchunks_TwinTwin, &
   Nchunks_SlipFamilies, Nchunks_TwinFamilies, Nchunks_nonSchmid, &
   structID, index_myFamily, index_otherFamily, &
   mySize=0_pInt, section = 0_pInt
 character(len=32) :: &
   structure  = ''
 character(len=65536) :: &
   tag  = '', &
   line = ''  
 
 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_PHENOPOWERLAW_label//' init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_PHENOPOWERLAW_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 Nchunks_SlipFamilies = lattice_maxNslipFamily
 Nchunks_TwinFamilies = lattice_maxNtwinFamily
 Nchunks_SlipSlip = lattice_maxNinteraction
 Nchunks_SlipTwin = lattice_maxNinteraction
 Nchunks_TwinSlip = lattice_maxNinteraction
 Nchunks_TwinTwin = lattice_maxNinteraction
 Nchunks_nonSchmid = lattice_maxNnonSchmid
 
 allocate(constitutive_phenopowerlaw_sizeDotState(maxNinstance),                  source=0_pInt)
 allocate(constitutive_phenopowerlaw_sizeState(maxNinstance),                     source=0_pInt)
 allocate(constitutive_phenopowerlaw_sizePostResults(maxNinstance),               source=0_pInt)
 allocate(constitutive_phenopowerlaw_sizePostResult(maxval(phase_Noutput),maxNinstance), &
                                                                                  source=0_pInt)
 allocate(constitutive_phenopowerlaw_output(maxval(phase_Noutput),maxNinstance))
          constitutive_phenopowerlaw_output               = ''
 allocate(constitutive_phenopowerlaw_outputID(maxval(phase_Noutput),maxNinstance),source=undefined_ID)
 allocate(constitutive_phenopowerlaw_Noutput(maxNinstance),                       source=0_pInt)
 allocate(constitutive_phenopowerlaw_structure(maxNinstance),                     source=0_pInt)
 allocate(constitutive_phenopowerlaw_Nslip(lattice_maxNslipFamily,maxNinstance),  source=0_pInt)
 allocate(constitutive_phenopowerlaw_Ntwin(lattice_maxNtwinFamily,maxNinstance),  source=0_pInt)
 allocate(constitutive_phenopowerlaw_totalNslip(maxNinstance),                    source=0_pInt)
 allocate(constitutive_phenopowerlaw_totalNtwin(maxNinstance),                    source=0_pInt)
 allocate(constitutive_phenopowerlaw_Cslip_66(6,6,maxNinstance),                  source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_gdot0_slip(maxNinstance),                    source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_n_slip(maxNinstance),                        source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_tau0_slip(lattice_maxNslipFamily,maxNinstance),  &
                                                                                  source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_tausat_slip(lattice_maxNslipFamily,maxNinstance), &
                                                                                  source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_gdot0_twin(maxNinstance),                    source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_n_twin(maxNinstance),                        source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_tau0_twin(lattice_maxNtwinFamily,maxNinstance), &
                                                                                  source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_spr(maxNinstance),                           source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_twinB(maxNinstance),                         source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_twinC(maxNinstance),                         source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_twinD(maxNinstance),                         source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_twinE(maxNinstance),                         source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_h0_SlipSlip(maxNinstance),                   source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_h0_SlipTwin(maxNinstance),                   source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_h0_TwinSlip(maxNinstance),                   source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_h0_TwinTwin(maxNinstance),                   source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_interaction_SlipSlip(lattice_maxNinteraction,maxNinstance), &
                                                                                  source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_interaction_SlipTwin(lattice_maxNinteraction,maxNinstance), &
                                                                                  source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_interaction_TwinSlip(lattice_maxNinteraction,maxNinstance), &
                                                                                  source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_interaction_TwinTwin(lattice_maxNinteraction,maxNinstance), &
                                                                                  source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_a_slip(maxNinstance),                        source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_aTolResistance(maxNinstance),                source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_aTolShear(maxNinstance),                     source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_aTolTwinfrac(maxNinstance),                  source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_nonSchmidCoeff(lattice_maxNnonSchmid,maxNinstance), &
                                                                                  source=0.0_pReal)

 rewind(fileUnit)
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo

 do while (trim(line) /= IO_EOF)                                                                    ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt                                                                     ! advance section counter
     if (phase_plasticity(section) == PLASTICITY_PHENOPOWERLAW_ID) then
       instance = phase_plasticityInstance(section)
       constitutive_phenopowerlaw_Cslip_66(1:6,1:6,instance)  = lattice_Cslip_66(1:6,1:6,section)
       constitutive_phenopowerlaw_structure(instance)  = lattice_structure(section)
       configNchunks = lattice_configNchunks(lattice_structureID(section))
       Nchunks_SlipFamilies = configNchunks(1)
       Nchunks_TwinFamilies = configNchunks(2)
       Nchunks_SlipSlip =     configNchunks(3)
       Nchunks_SlipTwin =     configNchunks(4)
       Nchunks_TwinSlip =     configNchunks(5)
       Nchunks_TwinTwin =     configNchunks(6)
       Nchunks_nonSchmid =    configNchunks(7)
     endif
     cycle                                                                                          ! skip to next line
   endif
   if (section > 0_pInt ) then; if (phase_plasticity(section) == PLASTICITY_PHENOPOWERLAW_ID) then  ! one of my sections. Do not short-circuit here (.and. between if-statements), it's not safe in Fortran
     instance = phase_plasticityInstance(section)                                                   ! which instance of my plasticity is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('plasticity','elasticity','lattice_structure','covera_ratio',&
             'c11','c12','c13','c22','c23','c33','c44','c55','c66')
       case ('(output)')
         constitutive_phenopowerlaw_Noutput(instance) = constitutive_phenopowerlaw_Noutput(instance) + 1_pInt
         constitutive_phenopowerlaw_output(constitutive_phenopowerlaw_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('resistance_slip')
             constitutive_phenopowerlaw_outputID(constitutive_phenopowerlaw_Noutput(instance),instance) = resistance_slip_ID
           case ('accumulatedshear_slip')
             constitutive_phenopowerlaw_outputID(constitutive_phenopowerlaw_Noutput(instance),instance) = accumulatedshear_slip_ID
           case ('shearrate_slip')
             constitutive_phenopowerlaw_outputID(constitutive_phenopowerlaw_Noutput(instance),instance) = shearrate_slip_ID
           case ('resolvedstress_slip')
             constitutive_phenopowerlaw_outputID(constitutive_phenopowerlaw_Noutput(instance),instance) = resolvedstress_slip_ID
           case ('totalshear')
             constitutive_phenopowerlaw_outputID(constitutive_phenopowerlaw_Noutput(instance),instance) = totalshear_ID
           case ('resistance_twin')
             constitutive_phenopowerlaw_outputID(constitutive_phenopowerlaw_Noutput(instance),instance) = resistance_twin_ID
           case ('accumulatedshear_twin')
             constitutive_phenopowerlaw_outputID(constitutive_phenopowerlaw_Noutput(instance),instance) = accumulatedshear_twin_ID
           case ('shearrate_twin')
             constitutive_phenopowerlaw_outputID(constitutive_phenopowerlaw_Noutput(instance),instance) = shearrate_twin_ID
           case ('resolvedstress_twin')
             constitutive_phenopowerlaw_outputID(constitutive_phenopowerlaw_Noutput(instance),instance) = resolvedstress_twin_ID
           case ('totalvolfrac')
             constitutive_phenopowerlaw_outputID(constitutive_phenopowerlaw_Noutput(instance),instance) = totalvolfrac_ID
           case default
             call IO_error(105_pInt,ext_msg=IO_stringValue(line,positions,2_pInt)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         end select
       case ('nslip')
         if (positions(1) < 1_pInt + Nchunks_SlipFamilies) &
           call IO_warning(50_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         Nchunks_SlipFamilies = positions(1) - 1_pInt
         do j = 1_pInt, Nchunks_SlipFamilies
            constitutive_phenopowerlaw_Nslip(j,instance) = IO_intValue(line,positions,1_pInt+j)
         enddo
       case ('gdot0_slip')
         constitutive_phenopowerlaw_gdot0_slip(instance) = IO_floatValue(line,positions,2_pInt)
       case ('n_slip')
         constitutive_phenopowerlaw_n_slip(instance) = IO_floatValue(line,positions,2_pInt)
       case ('tau0_slip')
         do j = 1_pInt,Nchunks_SlipFamilies
           constitutive_phenopowerlaw_tau0_slip(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('tausat_slip')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_phenopowerlaw_tausat_slip(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('a_slip', 'w0_slip')
         constitutive_phenopowerlaw_a_slip(instance) = IO_floatValue(line,positions,2_pInt)
       case ('ntwin')
         if (positions(1) < 1_pInt + Nchunks_TwinFamilies) &
           call IO_warning(51_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         Nchunks_TwinFamilies = positions(1) - 1_pInt
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_phenopowerlaw_Ntwin(j,instance) = IO_intValue(line,positions,1_pInt+j)
         enddo
       case ('gdot0_twin')
         constitutive_phenopowerlaw_gdot0_twin(instance) = IO_floatValue(line,positions,2_pInt)
       case ('n_twin')
         constitutive_phenopowerlaw_n_twin(instance) = IO_floatValue(line,positions,2_pInt)
       case ('tau0_twin')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_phenopowerlaw_tau0_twin(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('s_pr')
         constitutive_phenopowerlaw_spr(instance) = IO_floatValue(line,positions,2_pInt)
       case ('twin_b')
         constitutive_phenopowerlaw_twinB(instance) = IO_floatValue(line,positions,2_pInt)
       case ('twin_c')
         constitutive_phenopowerlaw_twinC(instance) = IO_floatValue(line,positions,2_pInt)
       case ('twin_d')
         constitutive_phenopowerlaw_twinD(instance) = IO_floatValue(line,positions,2_pInt)
       case ('twin_e')
         constitutive_phenopowerlaw_twinE(instance) = IO_floatValue(line,positions,2_pInt)
       case ('h0_slipslip')
         constitutive_phenopowerlaw_h0_SlipSlip(instance) = IO_floatValue(line,positions,2_pInt)
       case ('h0_sliptwin')
         constitutive_phenopowerlaw_h0_SlipTwin(instance) = IO_floatValue(line,positions,2_pInt)
         call IO_warning(42_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
       case ('h0_twinslip')
         constitutive_phenopowerlaw_h0_TwinSlip(instance) = IO_floatValue(line,positions,2_pInt)
       case ('h0_twintwin')
         constitutive_phenopowerlaw_h0_TwinTwin(instance) = IO_floatValue(line,positions,2_pInt)
       case ('atol_resistance')
         constitutive_phenopowerlaw_aTolResistance(instance) = IO_floatValue(line,positions,2_pInt)
       case ('atol_shear')
         constitutive_phenopowerlaw_aTolShear(instance)      = IO_floatValue(line,positions,2_pInt)
       case ('atol_twinfrac')
         constitutive_phenopowerlaw_aTolTwinfrac(instance)   = IO_floatValue(line,positions,2_pInt)
       case ('interaction_slipslip')
         if (positions(1) < 1_pInt + Nchunks_SlipSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         do j = 1_pInt, Nchunks_SlipSlip
           constitutive_phenopowerlaw_interaction_SlipSlip(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_sliptwin')
         if (positions(1) < 1_pInt + Nchunks_SlipTwin) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         do j = 1_pInt, Nchunks_SlipTwin
           constitutive_phenopowerlaw_interaction_SlipTwin(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_twinslip')
         if (positions(1) < 1_pInt + Nchunks_TwinSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         do j = 1_pInt, Nchunks_TwinSlip
           constitutive_phenopowerlaw_interaction_TwinSlip(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_twintwin')
         if (positions(1) < 1_pInt + Nchunks_TwinTwin) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         do j = 1_pInt, Nchunks_TwinTwin
           constitutive_phenopowerlaw_interaction_TwinTwin(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('nonschmid_coefficients')
         if (positions(1) < 1_pInt + Nchunks_nonSchmid) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
         do j = 1_pInt,Nchunks_nonSchmid
           constitutive_phenopowerlaw_nonSchmidCoeff(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case default
         call IO_error(210_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_PHENOPOWERLAW_label//')')
     end select
   endif; endif
 enddo

 sanityChecks: do instance = 1_pInt,maxNinstance           
   constitutive_phenopowerlaw_Nslip(1:lattice_maxNslipFamily,instance) = &
     min(lattice_NslipSystem(1:lattice_maxNslipFamily,constitutive_phenopowerlaw_structure(instance)),& ! limit active slip systems per family to min of available and requested
                                  constitutive_phenopowerlaw_Nslip(1:lattice_maxNslipFamily,instance))
   constitutive_phenopowerlaw_Ntwin(1:lattice_maxNtwinFamily,instance) = &
     min(lattice_NtwinSystem(1:lattice_maxNtwinFamily,constitutive_phenopowerlaw_structure(instance)),& ! limit active twin systems per family to min of available and requested
                                  constitutive_phenopowerlaw_Ntwin(:,instance))
   constitutive_phenopowerlaw_totalNslip(instance) = sum(constitutive_phenopowerlaw_Nslip(:,instance))            ! how many slip systems altogether
   constitutive_phenopowerlaw_totalNtwin(instance) = sum(constitutive_phenopowerlaw_Ntwin(:,instance))            ! how many twin systems altogether

   if (any(constitutive_phenopowerlaw_tau0_slip(:,instance) < 0.0_pReal .and. &
           constitutive_phenopowerlaw_Nslip(:,instance) > 0))        call IO_error(211_pInt,el=instance,ext_msg='tau0_slip (' &
                                                                 //PLASTICITY_PHENOPOWERLAW_label//')')
   if (constitutive_phenopowerlaw_gdot0_slip(instance) <= 0.0_pReal) call IO_error(211_pInt,el=instance,ext_msg='gdot0_slip (' &
                                                                 //PLASTICITY_PHENOPOWERLAW_label//')')
   if (constitutive_phenopowerlaw_n_slip(instance) <= 0.0_pReal)     call IO_error(211_pInt,el=instance,ext_msg='n_slip (' &
                                                                 //PLASTICITY_PHENOPOWERLAW_label//')')
   if (any(constitutive_phenopowerlaw_tausat_slip(:,instance) <= 0.0_pReal .and. &
           constitutive_phenopowerlaw_Nslip(:,instance) > 0))        call IO_error(211_pInt,el=instance,ext_msg='tausat_slip (' &
                                                                 //PLASTICITY_PHENOPOWERLAW_label//')')
   if (any(constitutive_phenopowerlaw_a_slip(instance) == 0.0_pReal .and. &
           constitutive_phenopowerlaw_Nslip(:,instance) > 0))        call IO_error(211_pInt,el=instance,ext_msg='a_slip (' &
                                                                 //PLASTICITY_PHENOPOWERLAW_label//')')
   if (any(constitutive_phenopowerlaw_tau0_twin(:,instance) < 0.0_pReal .and. &
           constitutive_phenopowerlaw_Ntwin(:,instance) > 0))        call IO_error(211_pInt,el=instance,ext_msg='tau0_twin (' &
                                                                 //PLASTICITY_PHENOPOWERLAW_label//')')
   if (    constitutive_phenopowerlaw_gdot0_twin(instance) <= 0.0_pReal .and. &
       any(constitutive_phenopowerlaw_Ntwin(:,instance) > 0))        call IO_error(211_pInt,el=instance,ext_msg='gdot0_twin (' &
                                                                 //PLASTICITY_PHENOPOWERLAW_label//')')
   if (    constitutive_phenopowerlaw_n_twin(instance) <= 0.0_pReal .and. &
       any(constitutive_phenopowerlaw_Ntwin(:,instance) > 0))        call IO_error(211_pInt,el=instance,ext_msg='n_twin (' &
                                                                 //PLASTICITY_PHENOPOWERLAW_label//')')
   if (constitutive_phenopowerlaw_aTolResistance(instance) <= 0.0_pReal) &
     constitutive_phenopowerlaw_aTolResistance(instance) = 1.0_pReal                                       ! default absolute tolerance 1 Pa
   if (constitutive_phenopowerlaw_aTolShear(instance) <= 0.0_pReal) &
     constitutive_phenopowerlaw_aTolShear(instance) = 1.0e-6_pReal                                         ! default absolute tolerance 1e-6
   if (constitutive_phenopowerlaw_aTolTwinfrac(instance) <= 0.0_pReal) &
     constitutive_phenopowerlaw_aTolTwinfrac(instance) = 1.0e-6_pReal                                      ! default absolute tolerance 1e-6

 enddo sanityChecks

!--------------------------------------------------------------------------------------------------
! allocation of variables whose size depends on the total number of active slip systems
 allocate(constitutive_phenopowerlaw_hardeningMatrix_SlipSlip(maxval(constitutive_phenopowerlaw_totalNslip),&   ! slip resistance from slip activity
                                                              maxval(constitutive_phenopowerlaw_totalNslip),&
                                                              maxNinstance), source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_hardeningMatrix_SlipTwin(maxval(constitutive_phenopowerlaw_totalNslip),&   ! slip resistance from twin activity
                                                              maxval(constitutive_phenopowerlaw_totalNtwin),&
                                                              maxNinstance), source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_hardeningMatrix_TwinSlip(maxval(constitutive_phenopowerlaw_totalNtwin),&   ! twin resistance from slip activity
                                                              maxval(constitutive_phenopowerlaw_totalNslip),&
                                                              maxNinstance), source=0.0_pReal)
 allocate(constitutive_phenopowerlaw_hardeningMatrix_TwinTwin(maxval(constitutive_phenopowerlaw_totalNtwin),&   ! twin resistance from twin activity
                                                              maxval(constitutive_phenopowerlaw_totalNtwin),&
                                                              maxNinstance), source=0.0_pReal)

 instancesLoop: do instance = 1_pInt,maxNinstance
   outputsLoop: do o = 1_pInt,constitutive_phenopowerlaw_Noutput(instance)
    select case(constitutive_phenopowerlaw_outputID(o,instance))
       case(resistance_slip_ID, &
            shearrate_slip_ID, &
            accumulatedshear_slip_ID, &
            resolvedstress_slip_ID &
            )
         mySize = constitutive_phenopowerlaw_totalNslip(instance)
       case(resistance_twin_ID, &
            shearrate_twin_ID, &
            accumulatedshear_twin_ID, &
            resolvedstress_twin_ID &
            )
         mySize = constitutive_phenopowerlaw_totalNtwin(instance)
       case(totalshear_ID, &
            totalvolfrac_ID &
            )
         mySize = 1_pInt
       case default
     end select

     outputFound: if (mySize > 0_pInt) then
       constitutive_phenopowerlaw_sizePostResult(o,instance) = mySize
       constitutive_phenopowerlaw_sizePostResults(instance)  = constitutive_phenopowerlaw_sizePostResults(instance) + mySize
     endif outputFound
   enddo outputsLoop 

   constitutive_phenopowerlaw_sizeDotState(instance) = constitutive_phenopowerlaw_totalNslip(instance)+ &
                                                constitutive_phenopowerlaw_totalNtwin(instance)+ &
                                                2_pInt + &
                                                constitutive_phenopowerlaw_totalNslip(instance)+ &
                                                constitutive_phenopowerlaw_totalNtwin(instance)            ! s_slip, s_twin, sum(gamma), sum(f), accshear_slip, accshear_twin
   constitutive_phenopowerlaw_sizeState(instance)    = constitutive_phenopowerlaw_sizeDotState(instance)

   structID = constitutive_phenopowerlaw_structure(instance)

   do f = 1_pInt,lattice_maxNslipFamily                                                             ! >>> interaction slip -- X
     index_myFamily = sum(constitutive_phenopowerlaw_Nslip(1:f-1_pInt,instance))
     do j = 1_pInt,constitutive_phenopowerlaw_Nslip(f,instance)                                            ! loop over (active) systems in my family (slip)
       do o = 1_pInt,lattice_maxNslipFamily
         index_otherFamily = sum(constitutive_phenopowerlaw_Nslip(1:o-1_pInt,instance))
         do k = 1_pInt,constitutive_phenopowerlaw_Nslip(o,instance)                                        ! loop over (active) systems in other family (slip)
           constitutive_phenopowerlaw_hardeningMatrix_SlipSlip(index_myFamily+j,index_otherFamily+k,instance) = &
               constitutive_phenopowerlaw_interaction_SlipSlip(lattice_interactionSlipSlip( &
                                                                 sum(lattice_NslipSystem(1:f-1,structID))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1,structID))+k, &
                                                                 structID), instance )
       enddo; enddo

       do o = 1_pInt,lattice_maxNtwinFamily
         index_otherFamily = sum(constitutive_phenopowerlaw_Ntwin(1:o-1_pInt,instance))
         do k = 1_pInt,constitutive_phenopowerlaw_Ntwin(o,instance)                                        ! loop over (active) systems in other family (twin)
           constitutive_phenopowerlaw_hardeningMatrix_SlipTwin(index_myFamily+j,index_otherFamily+k,instance) = &
               constitutive_phenopowerlaw_interaction_SlipTwin(lattice_interactionSlipTwin( &
                                                                 sum(lattice_NslipSystem(1:f-1_pInt,structID))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,structID))+k, &
                                                                 structID), instance )
       enddo; enddo

   enddo; enddo

   do f = 1_pInt,lattice_maxNtwinFamily                                                             ! >>> interaction twin -- X
     index_myFamily = sum(constitutive_phenopowerlaw_Ntwin(1:f-1_pInt,instance))
     do j = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,instance)                                            ! loop over (active) systems in my family (twin)

       do o = 1_pInt,lattice_maxNslipFamily
         index_otherFamily = sum(constitutive_phenopowerlaw_Nslip(1:o-1_pInt,instance))
         do k = 1_pInt,constitutive_phenopowerlaw_Nslip(o,instance)                                        ! loop over (active) systems in other family (slip)
           constitutive_phenopowerlaw_hardeningMatrix_TwinSlip(index_myFamily+j,index_otherFamily+k,instance) = &
               constitutive_phenopowerlaw_interaction_TwinSlip(lattice_interactionTwinSlip( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,structID))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1_pInt,structID))+k, &
                                                                 structID), instance )
       enddo; enddo

       do o = 1_pInt,lattice_maxNtwinFamily
         index_otherFamily = sum(constitutive_phenopowerlaw_Ntwin(1:o-1_pInt,instance))
         do k = 1_pInt,constitutive_phenopowerlaw_Ntwin(o,instance)                                        ! loop over (active) systems in other family (twin)
           constitutive_phenopowerlaw_hardeningMatrix_TwinTwin(index_myFamily+j,index_otherFamily+k,instance) = &
               constitutive_phenopowerlaw_interaction_TwinTwin(lattice_interactionTwinTwin( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,structID))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,structID))+k, &
                                                                 structID), instance )
       enddo; enddo

   enddo; enddo

 enddo instancesLoop

end subroutine constitutive_phenopowerlaw_init


!--------------------------------------------------------------------------------------------------
!> @brief sets the initial microstructural state for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
pure function constitutive_phenopowerlaw_stateInit(instance)
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily
 
 implicit none
 integer(pInt), intent(in) :: &
    instance                                                                                  !< number specifying the instance of the plasticity
 real(pReal), dimension(constitutive_phenopowerlaw_sizeDotState(instance)) :: &
   constitutive_phenopowerlaw_stateInit
 integer(pInt) :: &
   i

 constitutive_phenopowerlaw_stateInit = 0.0_pReal
 
 do i = 1_pInt,lattice_maxNslipFamily
   constitutive_phenopowerlaw_stateInit(1+&
                                        sum(constitutive_phenopowerlaw_Nslip(1:i-1,instance)) : &
                                        sum(constitutive_phenopowerlaw_Nslip(1:i  ,instance))) = &
     constitutive_phenopowerlaw_tau0_slip(i,instance)
 enddo

 do i = 1_pInt,lattice_maxNtwinFamily
   constitutive_phenopowerlaw_stateInit(1+sum(constitutive_phenopowerlaw_Nslip(:,instance))+&
                                        sum(constitutive_phenopowerlaw_Ntwin(1:i-1,instance)) : &
                                          sum(constitutive_phenopowerlaw_Nslip(:,instance))+&
                                        sum(constitutive_phenopowerlaw_Ntwin(1:i  ,instance))) = &
     constitutive_phenopowerlaw_tau0_twin(i,instance)
 enddo

end function constitutive_phenopowerlaw_stateInit


!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
pure function constitutive_phenopowerlaw_aTolState(instance)
 
 implicit none
 integer(pInt), intent(in) :: instance                                                         !< number specifying the instance of the plasticity
 
real(pReal), dimension(constitutive_phenopowerlaw_sizeState(instance)) :: &
   constitutive_phenopowerlaw_aTolState

 constitutive_phenopowerlaw_aTolState(1:constitutive_phenopowerlaw_totalNslip(instance)+ &
                                       constitutive_phenopowerlaw_totalNtwin(instance)) = &
          constitutive_phenopowerlaw_aTolResistance(instance)
 constitutive_phenopowerlaw_aTolState(1+constitutive_phenopowerlaw_totalNslip(instance)+ &
                                       constitutive_phenopowerlaw_totalNtwin(instance)) = &
          constitutive_phenopowerlaw_aTolShear(instance)
 constitutive_phenopowerlaw_aTolState(2+constitutive_phenopowerlaw_totalNslip(instance)+ &
                                       constitutive_phenopowerlaw_totalNtwin(instance)) = &
          constitutive_phenopowerlaw_aTolTwinFrac(instance)
 constitutive_phenopowerlaw_aTolState(3+constitutive_phenopowerlaw_totalNslip(instance)+ &
                                       constitutive_phenopowerlaw_totalNtwin(instance): &
                                     2+2*(constitutive_phenopowerlaw_totalNslip(instance)+ &
                                          constitutive_phenopowerlaw_totalNtwin(instance))) = &
          constitutive_phenopowerlaw_aTolShear(instance)

end function constitutive_phenopowerlaw_aTolState


!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenized elasticity matrix
!--------------------------------------------------------------------------------------------------
pure function constitutive_phenopowerlaw_homogenizedC(ipc,ip,el)
 use prec, only: &
   p_vec
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance
 
 implicit none
 real(pReal), dimension(6,6) :: &
   constitutive_phenopowerlaw_homogenizedC
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 constitutive_phenopowerlaw_homogenizedC = constitutive_phenopowerlaw_Cslip_66(1:6,1:6,&
                                              phase_plasticityInstance(material_phase(ipc,ip,el)))

end function constitutive_phenopowerlaw_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
pure subroutine constitutive_phenopowerlaw_LpAndItsTangent(Lp,dLp_dTstar99,Tstar_v,state,ipc,ip,el)
 use prec, only: &
   p_vec
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
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3),                                                  intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(9,9),                                                  intent(out) :: &
   dLp_dTstar99                                                                                     !< derivative of Lp with respect to 2nd Piola Kirchhoff stress

 real(pReal), dimension(6),                                                    intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),                                                                intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
   state                                                                                            !< microstructure state

 integer(pInt) :: &
   instance, & 
   nSlip, &
   nTwin,structID,index_Gamma,index_F,index_myFamily, &
   f,i,j,k,l,m,n
 real(pReal), dimension(3,3,3,3) :: &
   dLp_dTstar3333                                                                                   !< derivative of Lp with respect to Tstar as 4th order tensor
 real(pReal), dimension(3,3,2) :: &
   nonSchmid_tensor
 real(pReal), dimension(constitutive_phenopowerlaw_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip_pos,gdot_slip_neg,dgdot_dtauslip_pos,dgdot_dtauslip_neg,tau_slip_pos,tau_slip_neg
 real(pReal), dimension(constitutive_phenopowerlaw_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_twin,dgdot_dtautwin,tau_twin

 instance    = phase_plasticityInstance(material_phase(ipc,ip,el))
 structID = constitutive_phenopowerlaw_structure(instance)

 nSlip = constitutive_phenopowerlaw_totalNslip(instance)
 nTwin = constitutive_phenopowerlaw_totalNtwin(instance)
 
 index_Gamma = nSlip + nTwin + 1_pInt
 index_F     = nSlip + nTwin + 2_pInt

 Lp = 0.0_pReal
 dLp_dTstar3333 = 0.0_pReal
 dLp_dTstar99 = 0.0_pReal

 j = 0_pInt
 slipFamiliesLoop: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,structID))                                   ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,instance)                                          ! process each (active) slip system in family
     j = j+1_pInt
     
!--------------------------------------------------------------------------------------------------
! Calculation of Lp
     tau_slip_pos(j)  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,structID))
     tau_slip_neg(j)  = tau_slip_pos(j)
     nonSchmid_tensor(1:3,1:3,1) = lattice_Sslip(1:3,1:3,1,index_myFamily+i,structID)
     nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,1)
     do k = 1,lattice_NnonSchmid(structID) 
       tau_slip_pos(j) = tau_slip_pos(j) + constitutive_phenopowerlaw_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,structID))
       tau_slip_neg(j) = tau_slip_neg(j) + constitutive_phenopowerlaw_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,structID))
       nonSchmid_tensor(1:3,1:3,1) = nonSchmid_tensor(1:3,1:3,1) + constitutive_phenopowerlaw_nonSchmidCoeff(k,instance)*&
                                           lattice_Sslip(1:3,1:3,2*k,index_myFamily+i,structID)
       nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,2) + constitutive_phenopowerlaw_nonSchmidCoeff(k,instance)*&
                                           lattice_Sslip(1:3,1:3,2*k+1,index_myFamily+i,structID)
     enddo
     gdot_slip_pos(j) = 0.5_pReal*constitutive_phenopowerlaw_gdot0_slip(instance)* &
                    ((abs(tau_slip_pos(j))/state(ipc,ip,el)%p(j))**constitutive_phenopowerlaw_n_slip(instance))*&
                                                                    sign(1.0_pReal,tau_slip_pos(j))
     gdot_slip_neg(j) = 0.5_pReal*constitutive_phenopowerlaw_gdot0_slip(instance)* &
                    ((abs(tau_slip_neg(j))/state(ipc,ip,el)%p(j))**constitutive_phenopowerlaw_n_slip(instance))*&
                                                                    sign(1.0_pReal,tau_slip_neg(j))
     Lp = Lp + (1.0_pReal-state(ipc,ip,el)%p(index_F))*&                     ! 1-F
               (gdot_slip_pos(j)+gdot_slip_neg(j))*lattice_Sslip(1:3,1:3,1,index_myFamily+i,structID)

!--------------------------------------------------------------------------------------------------
! Calculation of the tangent of Lp
     if (gdot_slip_pos(j) /= 0.0_pReal) then
       dgdot_dtauslip_pos(j) = gdot_slip_pos(j)*constitutive_phenopowerlaw_n_slip(instance)/tau_slip_pos(j)
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtauslip_pos(j)*lattice_Sslip(k,l,1,index_myFamily+i,structID)* &
                                                     nonSchmid_tensor(m,n,1)
     endif
     
     if (gdot_slip_neg(j) /= 0.0_pReal) then
       dgdot_dtauslip_neg(j) = gdot_slip_neg(j)*constitutive_phenopowerlaw_n_slip(instance)/tau_slip_neg(j)
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtauslip_neg(j)*lattice_Sslip(k,l,1,index_myFamily+i,structID)* &
                                                     nonSchmid_tensor(m,n,2)
     endif
   enddo
 enddo slipFamiliesLoop

 j = 0_pInt
 twinFamiliesLoop: do f = 1_pInt,lattice_maxNtwinFamily
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                                   ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,instance)                                          ! process each (active) twin system in family
     j = j+1_pInt

!--------------------------------------------------------------------------------------------------
! Calculation of Lp
     tau_twin(j)  = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,structID)) 
     gdot_twin(j) = (1.0_pReal-state(ipc,ip,el)%p(index_F))*&                                       ! 1-F
                    constitutive_phenopowerlaw_gdot0_twin(instance)*&
                    (abs(tau_twin(j))/state(ipc,ip,el)%p(nSlip+j))**&
                    constitutive_phenopowerlaw_n_twin(instance)*max(0.0_pReal,sign(1.0_pReal,tau_twin(j)))
     Lp = Lp + gdot_twin(j)*lattice_Stwin(1:3,1:3,index_myFamily+i,structID)

!--------------------------------------------------------------------------------------------------
! Calculation of the tangent of Lp
     if (gdot_twin(j) /= 0.0_pReal) then
       dgdot_dtautwin(j) = gdot_twin(j)*constitutive_phenopowerlaw_n_twin(instance)/tau_twin(j)
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtautwin(j)*lattice_Stwin(k,l,index_myFamily+i,structID)* &
                                                     lattice_Stwin(m,n,index_myFamily+i,structID)
     endif
   enddo
 enddo twinFamiliesLoop

 dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)

end subroutine constitutive_phenopowerlaw_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
function constitutive_phenopowerlaw_dotState(Tstar_v,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use lattice, only: &
   lattice_Sslip_v, &
   lattice_Stwin_v, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_shearTwin, &
   lattice_NnonSchmid   
 use mesh, only: &
   mesh_NcpElems,&
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance
 
 implicit none
 real(pReal), dimension(6),                                                    intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),                                                                intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
   state                                                                                            !< microstructure state

 real(pReal), dimension(constitutive_phenopowerlaw_sizeDotState(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_phenopowerlaw_dotState

 integer(pInt) :: &
   instance,structID, &
   nSlip,nTwin, &
   f,i,j,k, &
   index_Gamma,index_F,index_myFamily, &
   offset_accshear_slip,offset_accshear_twin
 real(pReal) :: &
   c_SlipSlip,c_SlipTwin,c_TwinSlip,c_TwinTwin, &
   ssat_offset

 real(pReal), dimension(constitutive_phenopowerlaw_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip,tau_slip_pos,tau_slip_neg,left_SlipSlip,left_SlipTwin,right_SlipSlip,right_TwinSlip
 real(pReal), dimension(constitutive_phenopowerlaw_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_twin,tau_twin,left_TwinSlip,left_TwinTwin,right_SlipTwin,right_TwinTwin


 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 structID = constitutive_phenopowerlaw_structure(instance)
 
 nSlip = constitutive_phenopowerlaw_totalNslip(instance)
 nTwin = constitutive_phenopowerlaw_totalNtwin(instance)

 index_Gamma = nSlip + nTwin + 1_pInt
 index_F     = nSlip + nTwin + 2_pInt
 offset_accshear_slip = nSlip + nTwin + 2_pInt
 offset_accshear_twin = nSlip + nTwin + 2_pInt + nSlip

 constitutive_phenopowerlaw_dotState = 0.0_pReal
 
!--------------------------------------------------------------------------------------------------
! system-independent (nonlinear) prefactors to M_Xx (X influenced by x) matrices
 c_SlipSlip = constitutive_phenopowerlaw_h0_SlipSlip(instance)*&
              (1.0_pReal + constitutive_phenopowerlaw_twinC(instance)*state(ipc,ip,el)%p(index_F)**&
                                                           constitutive_phenopowerlaw_twinB(instance))
 c_SlipTwin = 0.0_pReal
 c_TwinSlip = constitutive_phenopowerlaw_h0_TwinSlip(instance)*&
              state(ipc,ip,el)%p(index_Gamma)**constitutive_phenopowerlaw_twinE(instance)
 c_TwinTwin = constitutive_phenopowerlaw_h0_TwinTwin(instance)*&
              state(ipc,ip,el)%p(index_F)**constitutive_phenopowerlaw_twinD(instance)

!--------------------------------------------------------------------------------------------------
!  calculate left and right vectors and calculate dot gammas
 ssat_offset = constitutive_phenopowerlaw_spr(instance)*sqrt(state(ipc,ip,el)%p(index_F))
 j = 0_pInt
 slipFamiliesLoop1: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,structID))                                   ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,instance)                                          ! process each (active) slip system in family
     j = j+1_pInt
     left_SlipSlip(j) = 1.0_pReal                                                                   ! no system-dependent left part
     left_SlipTwin(j) = 1.0_pReal                                                                   ! no system-dependent left part
     right_SlipSlip(j) = abs(1.0_pReal-state(ipc,ip,el)%p(j) / &
                                    (constitutive_phenopowerlaw_tausat_slip(f,instance)+ssat_offset)) &
                         **constitutive_phenopowerlaw_a_slip(instance)&
                         *sign(1.0_pReal,1.0_pReal-state(ipc,ip,el)%p(j) / &
                                    (constitutive_phenopowerlaw_tausat_slip(f,instance)+ssat_offset))
     right_TwinSlip(j) = 1.0_pReal                                                                  ! no system-dependent part
     
!--------------------------------------------------------------------------------------------------
! Calculation of dot gamma 
     tau_slip_pos(j)  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,structID))
     tau_slip_neg(j)  = tau_slip_pos(j)
     do k = 1,lattice_NnonSchmid(structID) 
       tau_slip_pos(j) = tau_slip_pos(j) + constitutive_phenopowerlaw_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,structID))
       tau_slip_neg(j) = tau_slip_neg(j) + constitutive_phenopowerlaw_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,structID))
     enddo
     gdot_slip(j) = constitutive_phenopowerlaw_gdot0_slip(instance)*0.5_pReal* &
                  ((abs(tau_slip_pos(j))/state(ipc,ip,el)%p(j))**constitutive_phenopowerlaw_n_slip(instance) &
                  +(abs(tau_slip_neg(j))/state(ipc,ip,el)%p(j))**constitutive_phenopowerlaw_n_slip(instance))&
                  *sign(1.0_pReal,tau_slip_pos(j)) 
   enddo
 enddo slipFamiliesLoop1

 j = 0_pInt
 twinFamiliesLoop1: do f = 1_pInt,lattice_maxNtwinFamily
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                                   ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,instance)                                          ! process each (active) twin system in family
     j = j+1_pInt
     left_TwinSlip(j)  = 1.0_pReal                                                                  ! no system-dependent right part
     left_TwinTwin(j)  = 1.0_pReal                                                                  ! no system-dependent right part
     right_SlipTwin(j) = 1.0_pReal                                                                  ! no system-dependent right part
     right_TwinTwin(j) = 1.0_pReal                                                                  ! no system-dependent right part

!--------------------------------------------------------------------------------------------------
! Calculation of dot vol frac
     tau_twin(j)  = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,structID)) 
     gdot_twin(j) = (1.0_pReal-state(ipc,ip,el)%p(index_F))*&                                       ! 1-F
                    constitutive_phenopowerlaw_gdot0_twin(instance)*&
                    (abs(tau_twin(j))/state(ipc,ip,el)%p(nSlip+j))**&
                    constitutive_phenopowerlaw_n_twin(instance)*max(0.0_pReal,sign(1.0_pReal,tau_twin(j)))
    enddo
  enddo twinFamiliesLoop1

!--------------------------------------------------------------------------------------------------
! calculate the overall hardening based on above
 j = 0_pInt
 slipFamiliesLoop2: do f = 1_pInt,lattice_maxNslipFamily
   do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,instance)                                          ! process each (active) slip system in family
     j = j+1_pInt
     constitutive_phenopowerlaw_dotState(j) = &                                                     ! evolution of slip resistance j
       c_SlipSlip * left_SlipSlip(j) * &
       dot_product(constitutive_phenopowerlaw_hardeningMatrix_SlipSlip(j,1:nSlip,instance), &
                   right_SlipSlip*abs(gdot_slip)) + &                                               ! dot gamma_slip modulated by right-side slip factor
       c_SlipTwin * left_SlipTwin(j) * &
       dot_product(constitutive_phenopowerlaw_hardeningMatrix_SlipTwin(j,1:nTwin,instance), &
                   right_SlipTwin*gdot_twin)                                                        ! dot gamma_twin modulated by right-side twin factor
     constitutive_phenopowerlaw_dotState(index_Gamma) = constitutive_phenopowerlaw_dotState(index_Gamma) + &
                                                        abs(gdot_slip(j))
     constitutive_phenopowerlaw_dotState(offset_accshear_slip+j) = abs(gdot_slip(j))
   enddo
 enddo slipFamiliesLoop2
 
 j = 0_pInt
 twinFamiliesLoop2: do f = 1_pInt,lattice_maxNtwinFamily
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                                   ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,instance)                                          ! process each (active) twin system in family
     j = j+1_pInt
     constitutive_phenopowerlaw_dotState(j+nSlip) = &                                               ! evolution of twin resistance j
       c_TwinSlip * left_TwinSlip(j) * &
       dot_product(constitutive_phenopowerlaw_hardeningMatrix_TwinSlip(j,1:nSlip,instance), &
                   right_TwinSlip*abs(gdot_slip)) + &                                               ! dot gamma_slip modulated by right-side slip factor
       c_TwinTwin * left_TwinTwin(j) * &
       dot_product(constitutive_phenopowerlaw_hardeningMatrix_TwinTwin(j,1:nTwin,instance), &
                   right_TwinTwin*gdot_twin)                                                        ! dot gamma_twin modulated by right-side twin factor
     if (state(ipc,ip,el)%p(index_F) < 0.98_pReal) &                                                ! ensure twin volume fractions stays below 1.0
       constitutive_phenopowerlaw_dotState(index_F) = constitutive_phenopowerlaw_dotState(index_F) + &
                                                      gdot_twin(j)/lattice_shearTwin(index_myFamily+i,structID)
     constitutive_phenopowerlaw_dotState(offset_accshear_twin+j) = abs(gdot_twin(j))
   enddo
 enddo twinFamiliesLoop2

end function constitutive_phenopowerlaw_dotState


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
pure function constitutive_phenopowerlaw_postResults(Tstar_v,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance, &
   phase_Noutput
 use lattice, only: &
   lattice_Sslip_v, &
   lattice_Stwin_v, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_NnonSchmid 
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips

 implicit none
 real(pReal), dimension(6),                                                    intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),                                                                intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
   state                                                                                            !< microstructure state

 real(pReal), dimension(constitutive_phenopowerlaw_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_phenopowerlaw_postResults

 integer(pInt) :: &
   instance,structID, &
   nSlip,nTwin, &
   o,f,i,c,j,k, &
   index_Gamma,index_F,index_accshear_slip,index_accshear_twin,index_myFamily 
 real(pReal) :: &
   tau_slip_pos,tau_slip_neg,tau


 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 structID = constitutive_phenopowerlaw_structure(instance)

 nSlip = constitutive_phenopowerlaw_totalNslip(instance)
 nTwin = constitutive_phenopowerlaw_totalNtwin(instance)

 index_Gamma = nSlip + nTwin + 1_pInt
 index_F     = nSlip + nTwin + 2_pInt
 index_accshear_slip = nSlip + nTwin + 3_pInt
 index_accshear_twin = nSlip + nTwin + 3_pInt + nSlip

 constitutive_phenopowerlaw_postResults = 0.0_pReal
 c = 0_pInt

 outputsLoop: do o = 1_pInt,phase_Noutput(material_phase(ipc,ip,el))
   select case(constitutive_phenopowerlaw_outputID(o,instance))
     case (resistance_slip_ID)
       constitutive_phenopowerlaw_postResults(c+1_pInt:c+nSlip) = state(ipc,ip,el)%p(1:nSlip)
       c = c + nSlip

     case (accumulatedshear_slip_ID)
       constitutive_phenopowerlaw_postResults(c+1_pInt:c+nSlip) = state(ipc,ip,el)%p(index_accshear_slip:&
                                                                                     index_accshear_slip+nSlip)
       c = c + nSlip

     case (shearrate_slip_ID)
       j = 0_pInt
       slipFamiliesLoop1: do f = 1_pInt,lattice_maxNslipFamily
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,structID))                             ! at which index starts my family
         do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,instance)                                    ! process each (active) slip system in family
           j = j + 1_pInt
           tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,structID))
           tau_slip_neg  = tau_slip_pos
           do k = 1,lattice_NnonSchmid(structID) 
             tau_slip_pos = tau_slip_pos + constitutive_phenopowerlaw_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,structID))
             tau_slip_neg = tau_slip_neg + constitutive_phenopowerlaw_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,structID))
           enddo
           constitutive_phenopowerlaw_postResults(c+j) = constitutive_phenopowerlaw_gdot0_slip(instance)*0.5_pReal* &
                    ((abs(tau_slip_pos)/state(ipc,ip,el)%p(j))**constitutive_phenopowerlaw_n_slip(instance) &
                    +(abs(tau_slip_neg)/state(ipc,ip,el)%p(j))**constitutive_phenopowerlaw_n_slip(instance))&
                    *sign(1.0_pReal,tau_slip_pos)
         enddo
       enddo slipFamiliesLoop1
       c = c + nSlip

     case (resolvedstress_slip_ID)
       j = 0_pInt
       slipFamiliesLoop2: do f = 1_pInt,lattice_maxNslipFamily
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,structID))                             ! at which index starts my family
         do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,instance)                                    ! process each (active) slip system in family
           j = j + 1_pInt
           constitutive_phenopowerlaw_postResults(c+j) = &
                             dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,structID))
         enddo
       enddo slipFamiliesLoop2
       c = c + nSlip

     case (totalshear_ID)
       constitutive_phenopowerlaw_postResults(c+1_pInt) = &
                             state(ipc,ip,el)%p(index_Gamma)
       c = c + 1_pInt

     case (resistance_twin_ID)
       constitutive_phenopowerlaw_postResults(c+1_pInt:c+nTwin) = &
                             state(ipc,ip,el)%p(1_pInt+nSlip:nTwin+nSlip)
       c = c + nTwin

     case (accumulatedshear_twin_ID)
       constitutive_phenopowerlaw_postResults(c+1_pInt:c+nTwin) = &
                             state(ipc,ip,el)%p(index_accshear_twin:index_accshear_twin+nTwin)
       c = c + nTwin

     case (shearrate_twin_ID)
       j = 0_pInt
       twinFamiliesLoop1: do f = 1_pInt,lattice_maxNtwinFamily
         index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                             ! at which index starts my family
         do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,instance)                                    ! process each (active) twin system in family
           j = j + 1_pInt
           tau = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,structID))
           constitutive_phenopowerlaw_postResults(c+j) = (1.0_pReal-state(ipc,ip,el)%p(index_F))*&          ! 1-F
                                                         constitutive_phenopowerlaw_gdot0_twin(instance)*&
                                                         (abs(tau)/state(ipc,ip,el)%p(j+nSlip))**&
                                           constitutive_phenopowerlaw_n_twin(instance)*max(0.0_pReal,sign(1.0_pReal,tau))
         enddo
       enddo twinFamiliesLoop1
       c = c + nTwin

     case (resolvedstress_twin_ID)
       j = 0_pInt
       twinFamiliesLoop2: do f = 1_pInt,lattice_maxNtwinFamily
         index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                             ! at which index starts my family
         do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,instance)                                    ! process each (active) twin system in family
           j = j + 1_pInt
           constitutive_phenopowerlaw_postResults(c+j) = &
                             dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,structID))
         enddo
       enddo twinFamiliesLoop2
       c = c + nTwin

     case (totalvolfrac_ID)
       constitutive_phenopowerlaw_postResults(c+1_pInt) = state(ipc,ip,el)%p(index_F)
       c = c + 1_pInt

   end select
 enddo outputsLoop

end function constitutive_phenopowerlaw_postResults

end module constitutive_phenopowerlaw
