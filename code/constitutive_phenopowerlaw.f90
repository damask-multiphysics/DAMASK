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
!*****************************************************
!*      Module: CONSTITUTIVE_PHENOPOWERLAW           *
!*****************************************************
!* contains:                                         *
!* - constitutive equations                          *
!* - parameters definition                           *
!*****************************************************

![Alu]
!plasticity              phenopowerlaw
!(output)                resistance_slip
!(output)                shearrate_slip
!(output)                resolvedstress_slip
!(output)                totalshear
!(output)                resistance_twin
!(output)                shearrate_twin
!(output)                resolvedstress_twin
!(output)                totalvolfrac
!lattice_structure       hex
!covera_ratio            1.587
!Nslip                   3 3 6 12               # per family
!Ntwin                   6 6 6 6                # per family
!
!c11                     162.2e9
!c12                     91.8e9
!c13                     68.8e9
!c33                     180.5e9
!c44                     46.7e9
!
!gdot0_slip              0.001
!n_slip                  50
!tau0_slip               65e6 22e6 52e6 50e6               # per family
!tausat_slip             80e6 180e6 140e6 140e6            # per family
!a_slip                  1
!gdot0_twin              0.001
!n_twin                  50
!tau0_twin               52e6 52e6 52e6 52e6              # per family
!s_pr                    50e6                             # push-up stress for slip saturation due to twinning
!twin_b                  2
!twin_C                  25
!twin_d                  0.1
!twin_e                  0.1
!h0_slipslip             10e6
!h0_sliptwin             0
!h0_twinslip             625e6
!h0_twintwin             400e6
!interaction_slipslip    5.5 5.5 1.0 52.0 5.5 5.5 1.0 52.0 27.5 0.2 72.8 1.0 72.8 72.8 27.5 1.1 1.4 5.5 7.7 7.7
!interaction_sliptwin    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
!interaction_twinslip    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
!interaction_twintwin    1 1 1 1 1 1 1 1 10 10 10 10 10 10 10 10 10 10 10 10 
!relevantResistance      1

module constitutive_phenopowerlaw

 use prec, only: pReal,pInt

 implicit none
 private
 character (len=*), parameter, public :: &
   constitutive_phenopowerlaw_label = 'phenopowerlaw'
    
 integer(pInt), dimension(:), allocatable, public :: &
   constitutive_phenopowerlaw_sizeDotState, &
   constitutive_phenopowerlaw_sizeState, &
   constitutive_phenopowerlaw_sizePostResults, &                                                    ! cumulative size of post results
   constitutive_phenopowerlaw_structure

 integer(pInt), dimension(:), allocatable, private :: &
   constitutive_phenopowerlaw_Noutput, &                                                            ! number of outputs per instance of this constitution 
   constitutive_phenopowerlaw_totalNslip, &                                                         ! no. of slip system used in simulation
   constitutive_phenopowerlaw_totalNtwin                                                            ! no. of twin system used in simulation
 
 integer(pInt), dimension(:,:), allocatable, target, public :: &
   constitutive_phenopowerlaw_sizePostResult                                                        ! size of each post result output

 integer(pInt), dimension(:,:), allocatable, private :: &
   constitutive_phenopowerlaw_Nslip, &                                                              ! active number of slip systems per family
   constitutive_phenopowerlaw_Ntwin                                                                 ! active number of twin systems per family

 character(len=64), dimension(:,:), allocatable, target, public :: & 
   constitutive_phenopowerlaw_output                                                                ! name of each post result output

 character(len=32), dimension(:), allocatable, private :: &
   constitutive_phenopowerlaw_structureName

 real(pReal), dimension(:), allocatable, private :: &
   constitutive_phenopowerlaw_CoverA, &
   constitutive_phenopowerlaw_C11, &
   constitutive_phenopowerlaw_C12, &
   constitutive_phenopowerlaw_C13, &
   constitutive_phenopowerlaw_C33, &
   constitutive_phenopowerlaw_C44, &
   constitutive_phenopowerlaw_gdot0_slip, &
   constitutive_phenopowerlaw_n_slip, &
   constitutive_phenopowerlaw_n_twin, &
   constitutive_phenopowerlaw_gdot0_twin

 real(pReal), dimension(:,:), allocatable, private :: &
   constitutive_phenopowerlaw_tau0_slip, &
   constitutive_phenopowerlaw_tausat_slip, &
   constitutive_phenopowerlaw_tau0_twin

 real(pReal), dimension(:), allocatable, private :: &
   constitutive_phenopowerlaw_spr, &
   constitutive_phenopowerlaw_twinB, &
   constitutive_phenopowerlaw_twinC, &
   constitutive_phenopowerlaw_twinD, &
   constitutive_phenopowerlaw_twinE, &
   constitutive_phenopowerlaw_h0_slipslip, &
   constitutive_phenopowerlaw_h0_sliptwin, &
   constitutive_phenopowerlaw_h0_twinslip, &
   constitutive_phenopowerlaw_h0_twintwin, &
   constitutive_phenopowerlaw_a_slip, &
   constitutive_phenopowerlaw_aTolResistance

 real(pReal), dimension(:,:), allocatable, private :: &
   constitutive_phenopowerlaw_interaction_slipslip, &
   constitutive_phenopowerlaw_interaction_sliptwin, &
   constitutive_phenopowerlaw_interaction_twinslip, &
   constitutive_phenopowerlaw_interaction_twintwin

 real(pReal), dimension(:,:,:), allocatable, private :: &
   constitutive_phenopowerlaw_hardeningMatrix_slipslip, &
   constitutive_phenopowerlaw_hardeningMatrix_sliptwin, &
   constitutive_phenopowerlaw_hardeningMatrix_twinslip, &
   constitutive_phenopowerlaw_hardeningMatrix_twintwin, &
   constitutive_phenopowerlaw_Cslip_66

 public :: &
   constitutive_phenopowerlaw_init, &
   constitutive_phenopowerlaw_homogenizedC, &
   constitutive_phenopowerlaw_aTolState, &
   constitutive_phenopowerlaw_dotState, &
   constitutive_phenopowerlaw_dotTemperature, &
   constitutive_phenopowerlaw_microstructure, &
   constitutive_phenopowerlaw_LpAndItsTangent, &
   constitutive_phenopowerlaw_postResults, &
   constitutive_phenopowerlaw_stateInit

contains

subroutine constitutive_phenopowerlaw_init(myFile)
!**************************************
!*      Module initialization         *
!**************************************
 use, intrinsic :: iso_fortran_env                                ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use math,    only: math_Mandel3333to66, &
                    math_Voigt66to3333
 use IO
 use material
 use debug,   only: debug_what,&
                    debug_constitutive,&
                    debug_levelBasic
 use lattice, only: lattice_initializeStructure, lattice_symmetryType, &
                    lattice_maxNslipFamily, lattice_maxNtwinFamily, &
                    lattice_maxNinteraction, lattice_NslipSystem, lattice_NtwinSystem, &
                    lattice_interactionSlipSlip, &
                    lattice_interactionSlipTwin, &
                    lattice_interactionTwinSlip, &
                    lattice_interactionTwinTwin

 implicit none
 integer(pInt), intent(in) :: myFile
 integer(pInt), parameter :: maxNchunks = lattice_maxNinteraction + 1_pInt
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) section, maxNinstance, i,j,k, f,o, &
               mySize, myStructure, index_myFamily, index_otherFamily
 character(len=64)   :: tag
 character(len=1024) :: line

 !$OMP CRITICAL (write2out)
   write(6,*)
   write(6,*) '<<<+-  constitutive_',trim(constitutive_phenopowerlaw_label),' init  -+>>>'
   write(6,*) '$Id$'
#include "compilation_info.f90"
 !$OMP END CRITICAL (write2out)
 
 maxNinstance = int(count(phase_plasticity == constitutive_phenopowerlaw_label),pInt)
 if (maxNinstance == 0) return

 if (iand(debug_what(debug_constitutive),debug_levelBasic) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
     write(6,'(a16,1x,i5)') '# instances:',maxNinstance
     write(6,*)
   !$OMP END CRITICAL (write2out)
 endif
 allocate(constitutive_phenopowerlaw_sizeDotState(maxNinstance))
          constitutive_phenopowerlaw_sizeDotState         = 0_pInt
 allocate(constitutive_phenopowerlaw_sizeState(maxNinstance))
          constitutive_phenopowerlaw_sizeState            = 0_pInt
 allocate(constitutive_phenopowerlaw_sizePostResults(maxNinstance))
          constitutive_phenopowerlaw_sizePostResults      = 0_pInt
 allocate(constitutive_phenopowerlaw_sizePostResult(maxval(phase_Noutput),maxNinstance))
          constitutive_phenopowerlaw_sizePostResult       = 0_pInt
 allocate(constitutive_phenopowerlaw_output(maxval(phase_Noutput),maxNinstance))
          constitutive_phenopowerlaw_output               = ''
 allocate(constitutive_phenopowerlaw_Noutput(maxNinstance))
          constitutive_phenopowerlaw_Noutput              = 0_pInt
 allocate(constitutive_phenopowerlaw_structureName(maxNinstance))
          constitutive_phenopowerlaw_structureName        = ''
 allocate(constitutive_phenopowerlaw_structure(maxNinstance))
          constitutive_phenopowerlaw_structure            = 0_pInt
 allocate(constitutive_phenopowerlaw_Nslip(lattice_maxNslipFamily,maxNinstance))
          constitutive_phenopowerlaw_Nslip                = 0_pInt           
 allocate(constitutive_phenopowerlaw_Ntwin(lattice_maxNtwinFamily,maxNinstance)) 
          constitutive_phenopowerlaw_Ntwin                = 0_pInt
 allocate(constitutive_phenopowerlaw_totalNslip(maxNinstance))
          constitutive_phenopowerlaw_totalNslip           = 0_pInt
 allocate(constitutive_phenopowerlaw_totalNtwin(maxNinstance))
          constitutive_phenopowerlaw_totalNtwin           = 0_pInt
 allocate(constitutive_phenopowerlaw_CoverA(maxNinstance)) 
          constitutive_phenopowerlaw_CoverA               = 0.0_pReal
 allocate(constitutive_phenopowerlaw_C11(maxNinstance))
          constitutive_phenopowerlaw_C11                  = 0.0_pReal
 allocate(constitutive_phenopowerlaw_C12(maxNinstance))
          constitutive_phenopowerlaw_C12                  = 0.0_pReal
 allocate(constitutive_phenopowerlaw_C13(maxNinstance))
          constitutive_phenopowerlaw_C13                  = 0.0_pReal
 allocate(constitutive_phenopowerlaw_C33(maxNinstance))
          constitutive_phenopowerlaw_C33                  = 0.0_pReal
 allocate(constitutive_phenopowerlaw_C44(maxNinstance))
          constitutive_phenopowerlaw_C44                  = 0.0_pReal
 allocate(constitutive_phenopowerlaw_Cslip_66(6,6,maxNinstance))
          constitutive_phenopowerlaw_Cslip_66             = 0.0_pReal
 allocate(constitutive_phenopowerlaw_gdot0_slip(maxNinstance))
          constitutive_phenopowerlaw_gdot0_slip           = 0.0_pReal
 allocate(constitutive_phenopowerlaw_n_slip(maxNinstance))
          constitutive_phenopowerlaw_n_slip               = 0.0_pReal
 allocate(constitutive_phenopowerlaw_tau0_slip(lattice_maxNslipFamily,maxNinstance))
          constitutive_phenopowerlaw_tau0_slip            = 0.0_pReal
 allocate(constitutive_phenopowerlaw_tausat_slip(lattice_maxNslipFamily,maxNinstance))
          constitutive_phenopowerlaw_tausat_slip          = 0.0_pReal
 allocate(constitutive_phenopowerlaw_gdot0_twin(maxNinstance))
          constitutive_phenopowerlaw_gdot0_twin           = 0.0_pReal
 allocate(constitutive_phenopowerlaw_n_twin(maxNinstance))
          constitutive_phenopowerlaw_n_twin               = 0.0_pReal
 allocate(constitutive_phenopowerlaw_tau0_twin(lattice_maxNtwinFamily,maxNinstance))
          constitutive_phenopowerlaw_tau0_twin            = 0.0_pReal
 allocate(constitutive_phenopowerlaw_spr(maxNinstance))
          constitutive_phenopowerlaw_spr                  = 0.0_pReal
 allocate(constitutive_phenopowerlaw_twinB(maxNinstance))
          constitutive_phenopowerlaw_twinB                = 0.0_pReal
 allocate(constitutive_phenopowerlaw_twinC(maxNinstance))
          constitutive_phenopowerlaw_twinC                = 0.0_pReal
 allocate(constitutive_phenopowerlaw_twinD(maxNinstance))
          constitutive_phenopowerlaw_twinD                = 0.0_pReal
 allocate(constitutive_phenopowerlaw_twinE(maxNinstance))
          constitutive_phenopowerlaw_twinE                = 0.0_pReal
 allocate(constitutive_phenopowerlaw_h0_slipslip(maxNinstance))
          constitutive_phenopowerlaw_h0_slipslip          = 0.0_pReal
 allocate(constitutive_phenopowerlaw_h0_sliptwin(maxNinstance))
          constitutive_phenopowerlaw_h0_sliptwin          = 0.0_pReal
 allocate(constitutive_phenopowerlaw_h0_twinslip(maxNinstance))  
          constitutive_phenopowerlaw_h0_twinslip          = 0.0_pReal
 allocate(constitutive_phenopowerlaw_h0_twintwin(maxNinstance))  
          constitutive_phenopowerlaw_h0_twintwin          = 0.0_pReal
 allocate(constitutive_phenopowerlaw_interaction_slipslip(lattice_maxNinteraction,maxNinstance))
          constitutive_phenopowerlaw_interaction_slipslip = 0.0_pReal
 allocate(constitutive_phenopowerlaw_interaction_sliptwin(lattice_maxNinteraction,maxNinstance))
          constitutive_phenopowerlaw_interaction_sliptwin = 0.0_pReal
 allocate(constitutive_phenopowerlaw_interaction_twinslip(lattice_maxNinteraction,maxNinstance))
          constitutive_phenopowerlaw_interaction_twinslip = 0.0_pReal
 allocate(constitutive_phenopowerlaw_interaction_twintwin(lattice_maxNinteraction,maxNinstance))
          constitutive_phenopowerlaw_interaction_twintwin = 0.0_pReal
 allocate(constitutive_phenopowerlaw_a_slip(maxNinstance))
          constitutive_phenopowerlaw_a_slip               = 0.0_pReal
 allocate(constitutive_phenopowerlaw_aTolResistance(maxNinstance))
          constitutive_phenopowerlaw_aTolResistance       = 0.0_pReal

 rewind(myFile)
 section = 0_pInt
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= 'phase')     ! wind forward to <phase>
   read(myFile,'(a1024)',END=100) line
 enddo

 do                                                       ! read thru sections of phase part
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                ! next section
     section = section + 1_pInt                           ! advance section counter
     cycle                                                ! skip to next line
   endif
   if (section > 0_pInt .and. phase_plasticity(section) == constitutive_phenopowerlaw_label) then  ! one of my sections
     i = phase_plasticityInstance(section)                ! which instance of my plasticity is present phase
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))   ! extract key
     select case(tag)
       case ('plasticity','elasticity')
         cycle
       case ('(output)')
         constitutive_phenopowerlaw_Noutput(i) = constitutive_phenopowerlaw_Noutput(i) + 1_pInt
         constitutive_phenopowerlaw_output(constitutive_phenopowerlaw_Noutput(i),i) = IO_lc(IO_stringValue(line,positions,2_pInt))
       case ('lattice_structure')
              constitutive_phenopowerlaw_structureName(i) = IO_lc(IO_stringValue(line,positions,2_pInt))
       case ('covera_ratio')
              constitutive_phenopowerlaw_CoverA(i) = IO_floatValue(line,positions,2_pInt)
       case ('c11')
              constitutive_phenopowerlaw_C11(i) = IO_floatValue(line,positions,2_pInt)
       case ('c12')
              constitutive_phenopowerlaw_C12(i) = IO_floatValue(line,positions,2_pInt)
       case ('c13')
              constitutive_phenopowerlaw_C13(i) = IO_floatValue(line,positions,2_pInt)
       case ('c33')
              constitutive_phenopowerlaw_C33(i) = IO_floatValue(line,positions,2_pInt)
       case ('c44')
              constitutive_phenopowerlaw_C44(i) = IO_floatValue(line,positions,2_pInt)
       case ('nslip')
              forall (j = 1_pInt:lattice_maxNslipFamily)&
                 constitutive_phenopowerlaw_Nslip(j,i) = IO_intValue(line,positions,1_pInt+j)
       case ('gdot0_slip')
              constitutive_phenopowerlaw_gdot0_slip(i) = IO_floatValue(line,positions,2_pInt)
       case ('n_slip')
              constitutive_phenopowerlaw_n_slip(i) = IO_floatValue(line,positions,2_pInt)
       case ('tau0_slip')
              forall (j = 1_pInt:lattice_maxNslipFamily)&
                 constitutive_phenopowerlaw_tau0_slip(j,i) = IO_floatValue(line,positions,1_pInt+j)
       case ('tausat_slip')
              forall (j = 1_pInt:lattice_maxNslipFamily)&
                 constitutive_phenopowerlaw_tausat_slip(j,i) = IO_floatValue(line,positions,1_pInt+j)
       case ('a_slip', 'w0_slip')
              constitutive_phenopowerlaw_a_slip(i) = IO_floatValue(line,positions,2_pInt)
       case ('ntwin')
              forall (j = 1_pInt:lattice_maxNtwinFamily)& 
                 constitutive_phenopowerlaw_Ntwin(j,i) = IO_intValue(line,positions,1_pInt+j)
       case ('gdot0_twin')
              constitutive_phenopowerlaw_gdot0_twin(i) = IO_floatValue(line,positions,2_pInt)
       case ('n_twin')
              constitutive_phenopowerlaw_n_twin(i) = IO_floatValue(line,positions,2_pInt)
       case ('tau0_twin')
              forall (j = 1_pInt:lattice_maxNtwinFamily)&
                 constitutive_phenopowerlaw_tau0_twin(j,i) = IO_floatValue(line,positions,1_pInt+j)
       case ('s_pr')
              constitutive_phenopowerlaw_spr(i) = IO_floatValue(line,positions,2_pInt)
       case ('twin_b')
              constitutive_phenopowerlaw_twinB(i) = IO_floatValue(line,positions,2_pInt)
       case ('twin_c')
              constitutive_phenopowerlaw_twinC(i) = IO_floatValue(line,positions,2_pInt)
       case ('twin_d')
              constitutive_phenopowerlaw_twinD(i) = IO_floatValue(line,positions,2_pInt)
       case ('twin_e')
              constitutive_phenopowerlaw_twinE(i) = IO_floatValue(line,positions,2_pInt)
       case ('h0_slipslip')
              constitutive_phenopowerlaw_h0_slipslip(i) = IO_floatValue(line,positions,2_pInt)
       case ('h0_sliptwin')
              constitutive_phenopowerlaw_h0_sliptwin(i) = IO_floatValue(line,positions,2_pInt)
       case ('h0_twinslip')
              constitutive_phenopowerlaw_h0_twinslip(i) = IO_floatValue(line,positions,2_pInt)
       case ('h0_twintwin')
              constitutive_phenopowerlaw_h0_twintwin(i) = IO_floatValue(line,positions,2_pInt)
       case ('atol_resistance')
              constitutive_phenopowerlaw_aTolResistance(i) = IO_floatValue(line,positions,2_pInt)
       case ('interaction_slipslip')
              forall (j = 1_pInt:lattice_maxNinteraction) &
                constitutive_phenopowerlaw_interaction_slipslip(j,i) = IO_floatValue(line,positions,1_pInt+j)
       case ('interaction_sliptwin')
              forall (j = 1_pInt:lattice_maxNinteraction) &
                constitutive_phenopowerlaw_interaction_sliptwin(j,i) = IO_floatValue(line,positions,1_pInt+j)
       case ('interaction_twinslip')
              forall (j = 1_pInt:lattice_maxNinteraction) &
                constitutive_phenopowerlaw_interaction_twinslip(j,i) = IO_floatValue(line,positions,1_pInt+j)
       case ('interaction_twintwin')
              forall (j = 1_pInt:lattice_maxNinteraction) &
                constitutive_phenopowerlaw_interaction_twintwin(j,i) = IO_floatValue(line,positions,1_pInt+j)
       case default
              call IO_error(220_pInt,ext_msg=tag)
     end select
   endif
 enddo

100 do i = 1_pInt,maxNinstance

   constitutive_phenopowerlaw_structure(i) = lattice_initializeStructure(constitutive_phenopowerlaw_structureName(i), &    ! get structure
                                                                         constitutive_phenopowerlaw_CoverA(i))
   constitutive_phenopowerlaw_Nslip(1:lattice_maxNslipFamily,i) = &
            min(lattice_NslipSystem(1:lattice_maxNslipFamily,constitutive_phenopowerlaw_structure(i)),&       ! limit active slip systems per family to min of available and requested
                constitutive_phenopowerlaw_Nslip(1:lattice_maxNslipFamily,i))
   constitutive_phenopowerlaw_Ntwin(1:lattice_maxNtwinFamily,i) = &
            min(lattice_NtwinSystem(1:lattice_maxNtwinFamily,constitutive_phenopowerlaw_structure(i)),&       ! limit active twin systems per family to min of available and requested
                constitutive_phenopowerlaw_Ntwin(:,i))
   constitutive_phenopowerlaw_totalNslip(i) = sum(constitutive_phenopowerlaw_Nslip(:,i))      ! how many slip systems altogether
   constitutive_phenopowerlaw_totalNtwin(i) = sum(constitutive_phenopowerlaw_Ntwin(:,i))      ! how many twin systems altogether

   if (constitutive_phenopowerlaw_structure(i) < 1 )                  call IO_error(205_pInt,e=i)
   if (any(constitutive_phenopowerlaw_tau0_slip(:,i) < 0.0_pReal .and. &
           constitutive_phenopowerlaw_Nslip(:,i) > 0))                call IO_error(221_pInt,e=i,ext_msg='tau0_slip')
   if (constitutive_phenopowerlaw_gdot0_slip(i) <= 0.0_pReal)         call IO_error(221_pInt,e=i,ext_msg='gdot0_slip')
   if (constitutive_phenopowerlaw_n_slip(i) <= 0.0_pReal)             call IO_error(221_pInt,e=i,ext_msg='n_slip')
   if (any(constitutive_phenopowerlaw_tausat_slip(:,i) <= 0.0_pReal .and. &
           constitutive_phenopowerlaw_Nslip(:,i) > 0))                call IO_error(221_pInt,e=i,ext_msg='tausat_slip')
   if (any(constitutive_phenopowerlaw_a_slip(i) == 0.0_pReal .and. &
           constitutive_phenopowerlaw_Nslip(:,i) > 0))                call IO_error(221_pInt,e=i,ext_msg='a_slip')
   if (any(constitutive_phenopowerlaw_tau0_twin(:,i) < 0.0_pReal .and. &
           constitutive_phenopowerlaw_Ntwin(:,i) > 0))                call IO_error(221_pInt,e=i,ext_msg='tau0_twin')
   if (    constitutive_phenopowerlaw_gdot0_twin(i) <= 0.0_pReal .and. &
       any(constitutive_phenopowerlaw_Ntwin(:,i) > 0))                call IO_error(221_pInt,e=i,ext_msg='gdot0_twin')
   if (    constitutive_phenopowerlaw_n_twin(i) <= 0.0_pReal .and. &
       any(constitutive_phenopowerlaw_Ntwin(:,i) > 0))                call IO_error(221_pInt,e=i,ext_msg='n_twin')
   if (constitutive_phenopowerlaw_aTolResistance(i) <= 0.0_pReal) &
     constitutive_phenopowerlaw_aTolResistance(i) = 1.0_pReal              ! default absolute tolerance 1 Pa

 enddo

 allocate(constitutive_phenopowerlaw_hardeningMatrix_slipslip(maxval(constitutive_phenopowerlaw_totalNslip),&   ! slip resistance from slip activity
                                                              maxval(constitutive_phenopowerlaw_totalNslip),&
                                                              maxNinstance))
 allocate(constitutive_phenopowerlaw_hardeningMatrix_sliptwin(maxval(constitutive_phenopowerlaw_totalNtwin),&   ! slip resistance from twin activity
                                                              maxval(constitutive_phenopowerlaw_totalNslip),&
                                                              maxNinstance))
 allocate(constitutive_phenopowerlaw_hardeningMatrix_twinslip(maxval(constitutive_phenopowerlaw_totalNslip),&   ! twin resistance from slip activity
                                                              maxval(constitutive_phenopowerlaw_totalNtwin),&
                                                              maxNinstance))
 allocate(constitutive_phenopowerlaw_hardeningMatrix_twintwin(maxval(constitutive_phenopowerlaw_totalNtwin),&   ! twin resistance from twin activity
                                                              maxval(constitutive_phenopowerlaw_totalNtwin),&
                                                              maxNinstance))
 constitutive_phenopowerlaw_hardeningMatrix_slipslip = 0.0_pReal
 constitutive_phenopowerlaw_hardeningMatrix_sliptwin = 0.0_pReal
 constitutive_phenopowerlaw_hardeningMatrix_twinslip = 0.0_pReal
 constitutive_phenopowerlaw_hardeningMatrix_twintwin = 0.0_pReal

 
 do i = 1_pInt,maxNinstance
   do j = 1_pInt,constitutive_phenopowerlaw_Noutput(i)
     select case(constitutive_phenopowerlaw_output(j,i))
       case('resistance_slip', &
            'shearrate_slip', &
            'resolvedstress_slip' &
            )
         mySize = constitutive_phenopowerlaw_totalNslip(i)
       case('resistance_twin', &
            'shearrate_twin', &
            'resolvedstress_twin' &
            )
         mySize = constitutive_phenopowerlaw_totalNtwin(i)
       case('totalshear', &
            'totalvolfrac' &
            )
         mySize = 1_pInt
       case default
         call IO_error(222_pInt,ext_msg=constitutive_phenopowerlaw_output(j,i))
     end select

     if (mySize > 0_pInt) then                               ! any meaningful output found
       constitutive_phenopowerlaw_sizePostResult(j,i) = mySize
       constitutive_phenopowerlaw_sizePostResults(i) = &
       constitutive_phenopowerlaw_sizePostResults(i) + mySize
     endif
   enddo

   constitutive_phenopowerlaw_sizeDotState(i) = constitutive_phenopowerlaw_totalNslip(i)+ &
                                                constitutive_phenopowerlaw_totalNtwin(i)+ 2_pInt    ! s_slip, s_twin, sum(gamma), sum(f)
   constitutive_phenopowerlaw_sizeState(i)    = constitutive_phenopowerlaw_totalNslip(i)+ &
                                                constitutive_phenopowerlaw_totalNtwin(i)+ 2_pInt    ! s_slip, s_twin, sum(gamma), sum(f)

   myStructure = constitutive_phenopowerlaw_structure(i)

   select case (lattice_symmetryType(myStructure))                                             ! assign elasticity tensor
     case(1_pInt)                                                                                   ! cubic(s)
       forall(k=1_pInt:3_pInt)
         forall(j=1_pInt:3_pInt) &
           constitutive_phenopowerlaw_Cslip_66(k,j,i) =   constitutive_phenopowerlaw_C12(i)
         constitutive_phenopowerlaw_Cslip_66(k,k,i) =     constitutive_phenopowerlaw_C11(i)
         constitutive_phenopowerlaw_Cslip_66(k+3_pInt,k+3_pInt,i) = constitutive_phenopowerlaw_C44(i)
       end forall
     case(2_pInt)                                                                                   ! hex
       constitutive_phenopowerlaw_Cslip_66(1,1,i) = constitutive_phenopowerlaw_C11(i)
       constitutive_phenopowerlaw_Cslip_66(2,2,i) = constitutive_phenopowerlaw_C11(i)
       constitutive_phenopowerlaw_Cslip_66(3,3,i) = constitutive_phenopowerlaw_C33(i)
       constitutive_phenopowerlaw_Cslip_66(1,2,i) = constitutive_phenopowerlaw_C12(i)
       constitutive_phenopowerlaw_Cslip_66(2,1,i) = constitutive_phenopowerlaw_C12(i)
       constitutive_phenopowerlaw_Cslip_66(1,3,i) = constitutive_phenopowerlaw_C13(i)
       constitutive_phenopowerlaw_Cslip_66(3,1,i) = constitutive_phenopowerlaw_C13(i)
       constitutive_phenopowerlaw_Cslip_66(2,3,i) = constitutive_phenopowerlaw_C13(i)
       constitutive_phenopowerlaw_Cslip_66(3,2,i) = constitutive_phenopowerlaw_C13(i)
       constitutive_phenopowerlaw_Cslip_66(4,4,i) = constitutive_phenopowerlaw_C44(i)
       constitutive_phenopowerlaw_Cslip_66(5,5,i) = constitutive_phenopowerlaw_C44(i)
       constitutive_phenopowerlaw_Cslip_66(6,6,i) = 0.5_pReal*(constitutive_phenopowerlaw_C11(i)- &
                                                                  constitutive_phenopowerlaw_C12(i))
   end select
   constitutive_phenopowerlaw_Cslip_66(:,:,i) = &
     math_Mandel3333to66(math_Voigt66to3333(constitutive_phenopowerlaw_Cslip_66(:,:,i)))

   do f = 1_pInt,lattice_maxNslipFamily                                          ! >>> interaction slip -- X
     index_myFamily = sum(constitutive_phenopowerlaw_Nslip(1:f-1_pInt,i))
     do j = 1_pInt,constitutive_phenopowerlaw_Nslip(f,i)                         ! loop over (active) systems in my family (slip)
       do o = 1_pInt,lattice_maxNslipFamily
         index_otherFamily = sum(constitutive_phenopowerlaw_Nslip(1:o-1_pInt,i))
         do k = 1_pInt,constitutive_phenopowerlaw_Nslip(o,i)                     ! loop over (active) systems in other family (slip)
           constitutive_phenopowerlaw_hardeningMatrix_slipslip(index_otherFamily+k,index_myFamily+j,i) = &
               constitutive_phenopowerlaw_interaction_slipslip(lattice_interactionSlipSlip( &
                                                        sum(lattice_NslipSystem(1:o-1,myStructure))+k, &
                                                        sum(lattice_NslipSystem(1:f-1,myStructure))+j, &
                                                        myStructure), i )
       enddo; enddo

       do o = 1_pInt,lattice_maxNtwinFamily
         index_otherFamily = sum(constitutive_phenopowerlaw_Ntwin(1:o-1_pInt,i))
         do k = 1_pInt,constitutive_phenopowerlaw_Ntwin(o,i)                     ! loop over (active) systems in other family (twin)
           constitutive_phenopowerlaw_hardeningMatrix_sliptwin(index_otherFamily+k,index_myFamily+j,i) = &
               constitutive_phenopowerlaw_interaction_sliptwin(lattice_interactionSlipTwin( &
                                                        sum(lattice_NtwinSystem(1:o-1_pInt,myStructure))+k, &
                                                        sum(lattice_NslipSystem(1:f-1_pInt,myStructure))+j, &
                                                        myStructure), i )
       enddo; enddo

   enddo; enddo

   do f = 1_pInt,lattice_maxNtwinFamily                                          ! >>> interaction twin -- X
     index_myFamily = sum(constitutive_phenopowerlaw_Ntwin(1:f-1_pInt,i))
     do j = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,i)                         ! loop over (active) systems in my family (twin)

       do o = 1_pInt,lattice_maxNslipFamily
         index_otherFamily = sum(constitutive_phenopowerlaw_Nslip(1:o-1_pInt,i))
         do k = 1_pInt,constitutive_phenopowerlaw_Nslip(o,i)                     ! loop over (active) systems in other family (slip)
           constitutive_phenopowerlaw_hardeningMatrix_twinslip(index_otherFamily+k,index_myFamily+j,i) = &
               constitutive_phenopowerlaw_interaction_twinslip(lattice_interactionTwinSlip( &
                                                        sum(lattice_NslipSystem(1:o-1_pInt,myStructure))+k, &
                                                        sum(lattice_NtwinSystem(1:f-1_pInt,myStructure))+j, &
                                                        myStructure), i )
       enddo; enddo

       do o = 1_pInt,lattice_maxNtwinFamily
         index_otherFamily = sum(constitutive_phenopowerlaw_Ntwin(1:o-1_pInt,i))
         do k = 1_pInt,constitutive_phenopowerlaw_Ntwin(o,i)                     ! loop over (active) systems in other family (twin)
           constitutive_phenopowerlaw_hardeningMatrix_twintwin(index_otherFamily+k,index_myFamily+j,i) = &
               constitutive_phenopowerlaw_interaction_twintwin(lattice_interactionTwinTwin( &
                                                        sum(lattice_NtwinSystem(1:o-1_pInt,myStructure))+k, &
                                                        sum(lattice_NtwinSystem(1:f-1_pInt,myStructure))+j, &
                                                        myStructure), i )
       enddo; enddo

   enddo; enddo

! report to out file...

 enddo

 return

end subroutine constitutive_phenopowerlaw_init


function constitutive_phenopowerlaw_stateInit(myInstance)
!*********************************************************************
!* initial microstructural state                                     *
!*********************************************************************
 use lattice, only: lattice_maxNslipFamily, lattice_maxNtwinFamily
 
 implicit none
 integer(pInt), intent(in) :: myInstance
 integer(pInt) :: i
 real(pReal), dimension(constitutive_phenopowerlaw_sizeDotState(myInstance)) :: constitutive_phenopowerlaw_stateInit

 constitutive_phenopowerlaw_stateInit = 0.0_pReal
 
 do i = 1_pInt,lattice_maxNslipFamily
   constitutive_phenopowerlaw_stateInit(1+&
                                        sum(constitutive_phenopowerlaw_Nslip(1:i-1,myInstance)) : &
                                        sum(constitutive_phenopowerlaw_Nslip(1:i  ,myInstance))) = &
     constitutive_phenopowerlaw_tau0_slip(i,myInstance)
 enddo

 do i = 1_pInt,lattice_maxNtwinFamily
   constitutive_phenopowerlaw_stateInit(1+sum(constitutive_phenopowerlaw_Nslip(:,myInstance))+&
                                        sum(constitutive_phenopowerlaw_Ntwin(1:i-1,myInstance)) : &
                                          sum(constitutive_phenopowerlaw_Nslip(:,myInstance))+&
                                        sum(constitutive_phenopowerlaw_Ntwin(1:i  ,myInstance))) = &
     constitutive_phenopowerlaw_tau0_twin(i,myInstance)
 enddo
 return

end function constitutive_phenopowerlaw_stateInit


!*********************************************************************
!* absolute state tolerance                                          *
!*********************************************************************
pure function constitutive_phenopowerlaw_aTolState(myInstance)

implicit none
!*** input variables
integer(pInt), intent(in) ::  myInstance                      ! number specifying the current instance of the plasticity

!*** output variables
real(pReal), dimension(constitutive_phenopowerlaw_sizeState(myInstance)) :: &
                              constitutive_phenopowerlaw_aTolState ! relevant state values for the current instance of this plasticity

!*** local variables

constitutive_phenopowerlaw_aTolState = constitutive_phenopowerlaw_aTolResistance(myInstance)

end function constitutive_phenopowerlaw_aTolState


function constitutive_phenopowerlaw_homogenizedC(state,ipc,ip,el)
!*********************************************************************
!* homogenized elacticity matrix                                     *
!* INPUT:                                                            *
!*  - state           : state variables                              *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec, only: p_vec
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_plasticityInstance
 
 implicit none
 integer(pInt), intent(in) :: ipc,ip,el
 integer(pInt) matID
 real(pReal), dimension(6,6) :: constitutive_phenopowerlaw_homogenizedC
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 
 matID = phase_plasticityInstance(material_phase(ipc,ip,el))
 constitutive_phenopowerlaw_homogenizedC = constitutive_phenopowerlaw_Cslip_66(:,:,matID)

 return

end function constitutive_phenopowerlaw_homogenizedC


subroutine constitutive_phenopowerlaw_microstructure(Temperature,state,ipc,ip,el)
!*********************************************************************
!* calculate derived quantities from state (not used here)           *
!* INPUT:                                                            *
!*  - Tp              : temperature                                  *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_plasticityInstance
 
 implicit none
 integer(pInt) ipc,ip,el, matID
 real(pReal) Temperature
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state

 matID = phase_plasticityInstance(material_phase(ipc,ip,el))
  
end subroutine constitutive_phenopowerlaw_microstructure


subroutine constitutive_phenopowerlaw_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,state,ipc,ip,el)
!*********************************************************************
!* plastic velocity gradient and its tangent                         *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - Lp              : plastic velocity gradient                    *
!*  - dLp_dTstar      : derivative of Lp (4th-rank tensor)           *
!*********************************************************************
 use prec, only: p_vec
 use math, only: math_Plain3333to99
 use lattice, only: lattice_Sslip,lattice_Sslip_v,lattice_Stwin,lattice_Stwin_v, lattice_maxNslipFamily, lattice_maxNtwinFamily, &
                    lattice_NslipSystem,lattice_NtwinSystem
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_plasticityInstance

 implicit none
 integer(pInt) ipc,ip,el
 integer(pInt) matID,nSlip,nTwin,f,i,j,k,l,m,n, structID,index_Gamma,index_F,index_myFamily
 real(pReal) Temperature
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(3,3) :: Lp
 real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
 real(pReal), dimension(9,9) :: dLp_dTstar
 real(pReal), dimension(constitutive_phenopowerlaw_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip,dgdot_dtauslip,tau_slip
 real(pReal), dimension(constitutive_phenopowerlaw_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_twin,dgdot_dtautwin,tau_twin

 matID    = phase_plasticityInstance(material_phase(ipc,ip,el))
 structID = constitutive_phenopowerlaw_structure(matID)

 nSlip = constitutive_phenopowerlaw_totalNslip(matID)
 nTwin = constitutive_phenopowerlaw_totalNtwin(matID)
 
 index_Gamma = nSlip + nTwin + 1_pInt
 index_F     = nSlip + nTwin + 2_pInt

 Lp = 0.0_pReal
 dLp_dTstar3333 = 0.0_pReal
 dLp_dTstar = 0.0_pReal

 j = 0_pInt
 do f = 1_pInt,lattice_maxNslipFamily                                             ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,structID))                 ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,matID)                        ! process each (active) slip system in family
     j = j+1_pInt

!* Calculation of Lp

     tau_slip(j)  = dot_product(Tstar_v,lattice_Sslip_v(1:6,index_myFamily+i,structID)) 
     gdot_slip(j) = constitutive_phenopowerlaw_gdot0_slip(matID)*(abs(tau_slip(j))/state(ipc,ip,el)%p(j))**&
                    constitutive_phenopowerlaw_n_slip(matID)*sign(1.0_pReal,tau_slip(j))
     Lp = Lp + (1.0_pReal-state(ipc,ip,el)%p(index_F))*&                     ! 1-F
               gdot_slip(j)*lattice_Sslip(1:3,1:3,index_myFamily+i,structID)

!* Calculation of the tangent of Lp

     if (gdot_slip(j) /= 0.0_pReal) then
       dgdot_dtauslip(j) = gdot_slip(j)*constitutive_phenopowerlaw_n_slip(matID)/tau_slip(j)
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtauslip(j)*lattice_Sslip(k,l,index_myFamily+i,structID)* &
                                                     lattice_Sslip(m,n,index_myFamily+i,structID)
     endif
   enddo
 enddo

 j = 0_pInt
 do f = 1_pInt,lattice_maxNtwinFamily                                             ! loop over all twin families
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                 ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,matID)                        ! process each (active) twin system in family
     j = j+1_pInt

!* Calculation of Lp

     tau_twin(j)  = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,structID)) 
     gdot_twin(j) = (1.0_pReal-state(ipc,ip,el)%p(index_F))*&                ! 1-F
                    constitutive_phenopowerlaw_gdot0_twin(matID)*&
                    (abs(tau_twin(j))/state(ipc,ip,el)%p(nSlip+j))**&
                    constitutive_phenopowerlaw_n_twin(matID)*max(0.0_pReal,sign(1.0_pReal,tau_twin(j)))
     Lp = Lp + gdot_twin(j)*lattice_Stwin(1:3,1:3,index_myFamily+i,structID)

!* Calculation of the tangent of Lp

     if (gdot_twin(j) /= 0.0_pReal) then
       dgdot_dtautwin(j) = gdot_twin(j)*constitutive_phenopowerlaw_n_twin(matID)/tau_twin(j)
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtautwin(j)*lattice_Stwin(k,l,index_myFamily+i,structID)* &
                                                     lattice_Stwin(m,n,index_myFamily+i,structID)
     endif
   enddo
 enddo

 dLp_dTstar = math_Plain3333to99(dLp_dTstar3333)

 return
end subroutine constitutive_phenopowerlaw_LpAndItsTangent


function constitutive_phenopowerlaw_dotState(Tstar_v,Temperature,state,ipc,ip,el)
!*********************************************************************
!* rate of change of microstructure                                  *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - constitutive_dotState : evolution of state variable            *
!*********************************************************************
 use prec,     only: p_vec
 use lattice,  only: lattice_Sslip_v, lattice_Stwin_v, lattice_maxNslipFamily, lattice_maxNtwinFamily, &
                     lattice_NslipSystem,lattice_NtwinSystem,lattice_shearTwin   
 use mesh,     only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_plasticityInstance
 
 implicit none
 integer(pInt) ipc,ip,el
 integer(pInt) matID,nSlip,nTwin,f,i,j, structID,index_Gamma,index_F,index_myFamily 
 real(pReal) Temperature,c_slipslip,c_sliptwin,c_twinslip,c_twintwin, ssat_offset
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(constitutive_phenopowerlaw_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip,tau_slip,h_slipslip,h_sliptwin
 real(pReal), dimension(constitutive_phenopowerlaw_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_twin,tau_twin,h_twinslip,h_twintwin
 real(pReal), dimension(constitutive_phenopowerlaw_sizeDotState(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_phenopowerlaw_dotState

 matID = phase_plasticityInstance(material_phase(ipc,ip,el))
 structID = constitutive_phenopowerlaw_structure(matID)
 
 nSlip = constitutive_phenopowerlaw_totalNslip(matID)
 nTwin = constitutive_phenopowerlaw_totalNtwin(matID)

 index_Gamma = nSlip + nTwin + 1_pInt
 index_F     = nSlip + nTwin + 2_pInt

 constitutive_phenopowerlaw_dotState = 0.0_pReal

!-- system-independent (nonlinear) prefactors to M_xx matrices

 c_slipslip = constitutive_phenopowerlaw_h0_slipslip(matID)*&
              (1.0_pReal + &
               constitutive_phenopowerlaw_twinC(matID)*state(ipc,ip,el)%p(index_F)**constitutive_phenopowerlaw_twinB(matID))
 c_sliptwin = 0.0_pReal
 c_twinslip = constitutive_phenopowerlaw_h0_twinslip(matID)*&
              state(ipc,ip,el)%p(index_Gamma)**constitutive_phenopowerlaw_twinE(matID)
 c_twintwin = constitutive_phenopowerlaw_h0_twintwin(matID)*&
              state(ipc,ip,el)%p(index_F)**constitutive_phenopowerlaw_twinD(matID)

!-- add system-dependent part and calculate dot gammas

 ssat_offset = constitutive_phenopowerlaw_spr(matID)*sqrt(state(ipc,ip,el)%p(index_F))
 j = 0_pInt
 do f = 1_pInt,lattice_maxNslipFamily                                             ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,structID))                 ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,matID)                        ! process each (active) slip system in family
     j = j+1_pInt
     h_slipslip(j) = c_slipslip*(1.0_pReal-state(ipc,ip,el)%p(j) / &         ! system-dependent prefactor for slip--slip interaction
                                 (constitutive_phenopowerlaw_tausat_slip(f,matID)+ssat_offset))** &
                                constitutive_phenopowerlaw_a_slip(matID)
                     
     h_sliptwin(j) = c_sliptwin                                              ! no system-dependent part
     
!* Calculation of dot gamma 

     tau_slip(j)  = dot_product(Tstar_v,lattice_Sslip_v(1:6,index_myFamily+i,structID)) 
     gdot_slip(j) = constitutive_phenopowerlaw_gdot0_slip(matID)*(abs(tau_slip(j))/state(ipc,ip,el)%p(j))**&
                    constitutive_phenopowerlaw_n_slip(matID)*sign(1.0_pReal,tau_slip(j)) 
    enddo
  enddo

 j = 0_pInt
 do f = 1_pInt,lattice_maxNtwinFamily                                             ! loop over all twin families
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                 ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,matID)                        ! process each (active) twin system in family
     j = j+1_pInt
     h_twinslip(j) = c_twinslip                                              ! no system-dependent parts
     h_twintwin(j) = c_twintwin

!* Calculation of dot vol frac

     tau_twin(j)  = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,structID)) 
     gdot_twin(j) = (1.0_pReal-state(ipc,ip,el)%p(index_F))*&                ! 1-F
                    constitutive_phenopowerlaw_gdot0_twin(matID)*&
                    (abs(tau_twin(j))/state(ipc,ip,el)%p(nSlip+j))**&
                    constitutive_phenopowerlaw_n_twin(matID)*max(0.0_pReal,sign(1.0_pReal,tau_twin(j)))
    enddo
  enddo

!-- calculate the overall hardening based on above

 j = 0_pInt
 do f = 1_pInt,lattice_maxNslipFamily                                             ! loop over all slip families
   do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,matID)                        ! process each (active) slip system in family
     j = j+1_pInt
     constitutive_phenopowerlaw_dotState(j) = &                                                                             ! evolution of slip resistance j
       h_slipslip(j) * dot_product(constitutive_phenopowerlaw_hardeningMatrix_slipslip(1:nSlip,j,matID),abs(gdot_slip)) + & ! dot gamma_slip
       h_sliptwin(j) * dot_product(constitutive_phenopowerlaw_hardeningMatrix_sliptwin(1:nTwin,j,matID),gdot_twin)          ! dot gamma_twin
     constitutive_phenopowerlaw_dotState(index_Gamma) = constitutive_phenopowerlaw_dotState(index_Gamma) + &
                                                        abs(gdot_slip(j))
   enddo
 enddo
 
 j = 0_pInt
 do f = 1_pInt,lattice_maxNtwinFamily                                             ! loop over all twin families
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                 ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,matID)                        ! process each (active) twin system in family
     j = j+1_pInt
     constitutive_phenopowerlaw_dotState(j+nSlip) = &                                                                       ! evolution of twin resistance j
       h_twinslip(j) * dot_product(constitutive_phenopowerlaw_hardeningMatrix_twinslip(1:nSlip,j,matID),abs(gdot_slip)) + & ! dot gamma_slip
       h_twintwin(j) * dot_product(constitutive_phenopowerlaw_hardeningMatrix_twintwin(1:nTwin,j,matID),gdot_twin)          ! dot gamma_twin
     constitutive_phenopowerlaw_dotState(index_F) = constitutive_phenopowerlaw_dotState(index_F) + &
                                                    gdot_twin(j)/lattice_shearTwin(index_myFamily+i,structID)
   enddo
 enddo

end function constitutive_phenopowerlaw_dotState


!****************************************************************
!* calculates the rate of change of temperature                 *
!****************************************************************
pure function constitutive_phenopowerlaw_dotTemperature(Tstar_v,Temperature,state,ipc,ip,el)

  !*** variables and functions from other modules ***!
  use prec,     only: pReal,pInt,p_vec
  use mesh,     only: mesh_NcpElems, mesh_maxNips
  use material, only: homogenization_maxNgrains
  
  implicit none
  !*** input variables ***!
  real(pReal), dimension(6), intent(in) ::  Tstar_v                   ! 2nd Piola Kirchhoff stress tensor in Mandel notation
  real(pReal), intent(in) ::                Temperature
  integer(pInt), intent(in)::               ipc, &                    ! grain number
                                            ip, &                     ! integration point number
                                            el                        ! element number
  type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state ! state of the current microstructure
  
  !*** output variables ***!
  real(pReal) constitutive_phenopowerlaw_dotTemperature               ! rate of change of temparature
  
  ! calculate dotTemperature
  constitutive_phenopowerlaw_dotTemperature = 0.0_pReal

end function constitutive_phenopowerlaw_dotTemperature



pure function constitutive_phenopowerlaw_postResults(Tstar_v,Temperature,dt,state,ipc,ip,el)
!*********************************************************************
!* return array of constitutive results                              *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - dt              : current time increment                       *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec, only: pReal,pInt,p_vec
 use lattice, only: lattice_Sslip_v,lattice_Stwin_v, lattice_maxNslipFamily, lattice_maxNtwinFamily, &
                    lattice_NslipSystem,lattice_NtwinSystem   
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase,phase_plasticityInstance,phase_Noutput
 
 implicit none
 integer(pInt), intent(in) :: ipc,ip,el
 real(pReal), intent(in) :: dt,Temperature
 real(pReal), dimension(6), intent(in) :: Tstar_v
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
 integer(pInt) matID,o,f,i,c,nSlip,nTwin,j, structID,index_Gamma,index_F,index_myFamily 
 real(pReal) tau
 real(pReal), dimension(constitutive_phenopowerlaw_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_phenopowerlaw_postResults

 matID = phase_plasticityInstance(material_phase(ipc,ip,el))
 structID = constitutive_phenopowerlaw_structure(matID)

 nSlip = constitutive_phenopowerlaw_totalNslip(matID)
 nTwin = constitutive_phenopowerlaw_totalNtwin(matID)

 index_Gamma = nSlip + nTwin + 1_pInt
 index_F     = nSlip + nTwin + 2_pInt

 constitutive_phenopowerlaw_postResults = 0.0_pReal
 c = 0_pInt

 do o = 1_pInt,phase_Noutput(material_phase(ipc,ip,el))
   select case(constitutive_phenopowerlaw_output(o,matID))
     case ('resistance_slip')
       constitutive_phenopowerlaw_postResults(c+1_pInt:c+nSlip) = state(ipc,ip,el)%p(1:nSlip)
       c = c + nSlip

     case ('shearrate_slip')
       j = 0_pInt
       do f = 1_pInt,lattice_maxNslipFamily                                             ! loop over all slip families
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,structID))                 ! at which index starts my family
         do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,matID)                        ! process each (active) slip system in family
           j = j + 1_pInt
           tau = dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,structID))
           constitutive_phenopowerlaw_postResults(c+j) = constitutive_phenopowerlaw_gdot0_slip(matID)*&
                                                         (abs(tau)/state(ipc,ip,el)%p(j))**&
                                                         constitutive_phenopowerlaw_n_slip(matID)*sign(1.0_pReal,tau)
       enddo; enddo
       c = c + nSlip

     case ('resolvedstress_slip')
       j = 0_pInt
       do f = 1_pInt,lattice_maxNslipFamily                                             ! loop over all slip families
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,structID))                 ! at which index starts my family
         do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,matID)                        ! process each (active) slip system in family
           j = j + 1_pInt
           constitutive_phenopowerlaw_postResults(c+j) = dot_product(Tstar_v,lattice_Sslip_v(1:6,index_myFamily+i,structID))
       enddo; enddo
       c = c + nSlip

     case ('totalshear')
       constitutive_phenopowerlaw_postResults(c+1_pInt) = state(ipc,ip,el)%p(index_Gamma)
       c = c + 1_pInt

     case ('resistance_twin')
       constitutive_phenopowerlaw_postResults(c+1_pInt:c+nTwin) = state(ipc,ip,el)%p(1_pInt+nSlip:nTwin+nSlip)
       c = c + nTwin

     case ('shearrate_twin')
       j = 0_pInt
       do f = 1_pInt,lattice_maxNtwinFamily                                             ! loop over all twin families
         index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                 ! at which index starts my family
         do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,matID)                        ! process each (active) twin system in family
           j = j + 1_pInt
           tau = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,structID))
           constitutive_phenopowerlaw_postResults(c+j) = (1.0_pReal-state(ipc,ip,el)%p(index_F))*&          ! 1-F
                                                         constitutive_phenopowerlaw_gdot0_twin(matID)*&
                                                         (abs(tau)/state(ipc,ip,el)%p(j+nSlip))**&
                                                         constitutive_phenopowerlaw_n_twin(matID)*max(0.0_pReal,sign(1.0_pReal,tau))
       enddo; enddo
       c = c + nTwin

     case ('resolvedstress_twin')
       j = 0_pInt
       do f = 1_pInt,lattice_maxNtwinFamily                                             ! loop over all twin families
         index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                 ! at which index starts my family
         do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,matID)                        ! process each (active) twin system in family
           j = j + 1_pInt
           constitutive_phenopowerlaw_postResults(c+j) = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,structID))
       enddo; enddo
       c = c + nTwin

     case ('totalvolfrac')
       constitutive_phenopowerlaw_postResults(c+1_pInt) = state(ipc,ip,el)%p(index_F)
       c = c + 1_pInt

   end select
 enddo

end function constitutive_phenopowerlaw_postResults

end module constitutive_phenopowerlaw
