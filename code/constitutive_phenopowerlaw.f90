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
!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief phenomenological crystal plasticity formulation using a powerlaw fitting
!--------------------------------------------------------------------------------------------------
module constitutive_phenopowerlaw

 use prec, only: pReal,pInt

 implicit none
 private
 character (len=*), parameter, public :: &
   constitutive_phenopowerlaw_label = 'phenopowerlaw'
    
 integer(pInt), dimension(:), allocatable, public :: &
   constitutive_phenopowerlaw_sizeDotState, &
   constitutive_phenopowerlaw_sizeState, &
   constitutive_phenopowerlaw_sizePostResults, &                                                    !< cumulative size of post results
   constitutive_phenopowerlaw_structure

 integer(pInt), dimension(:), allocatable, private :: &
   constitutive_phenopowerlaw_Noutput, &                                                            !< number of outputs per instance of this constitution 
   constitutive_phenopowerlaw_totalNslip, &                                                         !< no. of slip system used in simulation
   constitutive_phenopowerlaw_totalNtwin                                                            !< no. of twin system used in simulation
 
 integer(pInt), dimension(:,:), allocatable, target, public :: &
   constitutive_phenopowerlaw_sizePostResult                                                        !< size of each post result output

 integer(pInt), dimension(:,:), allocatable, private :: &
   constitutive_phenopowerlaw_Nslip, &                                                              !< active number of slip systems per family (input parameter, per family)
   constitutive_phenopowerlaw_Ntwin                                                                 !< active number of twin systems per family (input parameter, per family)

 character(len=64), dimension(:,:), allocatable, target, public :: & 
   constitutive_phenopowerlaw_output                                                                !< name of each post result output

 character(len=32), dimension(:), allocatable, public :: &
   constitutive_phenopowerlaw_structureName

 real(pReal), dimension(:), allocatable, private :: &
   constitutive_phenopowerlaw_CoverA, &                                                             !< c/a of the crystal (input parameter)
   constitutive_phenopowerlaw_gdot0_slip, &                                                         !< reference shear strain rate for slip (input parameter)
   constitutive_phenopowerlaw_gdot0_twin, &                                                         !< reference shear strain rate for twin (input parameter)
   constitutive_phenopowerlaw_n_slip, &                                                             !< stress exponent for slip (input parameter)
   constitutive_phenopowerlaw_n_twin                                                                !< stress exponent for twin (input parameter)

 real(pReal), dimension(:,:), allocatable, private :: &
   constitutive_phenopowerlaw_tau0_slip, &                                                          !< initial critical shear stress for slip (input parameter, per family)
   constitutive_phenopowerlaw_tau0_twin, &                                                          !< initial critical shear stress for twin (input parameter, per family)
   constitutive_phenopowerlaw_tausat_slip                                                           !< maximum critical shear stress for slip (input parameter, per family)

 real(pReal), dimension(:), allocatable, private :: &
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

 real(pReal), dimension(:,:), allocatable, private :: &
   constitutive_phenopowerlaw_interaction_SlipSlip, &                                               !< interaction factors slip - slip (input parameter)
   constitutive_phenopowerlaw_interaction_SlipTwin, &                                               !< interaction factors slip - twin (input parameter)
   constitutive_phenopowerlaw_interaction_TwinSlip, &                                               !< interaction factors twin - slip (input parameter)
   constitutive_phenopowerlaw_interaction_TwinTwin                                                  !< interaction factors twin - twin (input parameter)

 real(pReal), dimension(:,:,:), allocatable, private :: &
   constitutive_phenopowerlaw_hardeningMatrix_SlipSlip, &
   constitutive_phenopowerlaw_hardeningMatrix_SlipTwin, &
   constitutive_phenopowerlaw_hardeningMatrix_TwinSlip, &
   constitutive_phenopowerlaw_hardeningMatrix_TwinTwin, &
   constitutive_phenopowerlaw_Cslip_66

 real(pReal), dimension(:,:), allocatable, private :: &
   constitutive_phenopowerlaw_nonSchmidCoeff

 public :: &
   constitutive_phenopowerlaw_init, &
   constitutive_phenopowerlaw_homogenizedC, &
   constitutive_phenopowerlaw_aTolState, &
   constitutive_phenopowerlaw_dotState, &
   constitutive_phenopowerlaw_deltaState, &
   constitutive_phenopowerlaw_dotTemperature, &
   constitutive_phenopowerlaw_microstructure, &
   constitutive_phenopowerlaw_LpAndItsTangent, &
   constitutive_phenopowerlaw_postResults, &
   constitutive_phenopowerlaw_stateInit

contains

!--------------------------------------------------------------------------------------------------
!> @brief reading in parameters from material config and doing consistency checks
!--------------------------------------------------------------------------------------------------
subroutine constitutive_phenopowerlaw_init(myFile)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use math,    only: math_Mandel3333to66, &
                    math_Voigt66to3333
 use IO
 use material
 use debug,   only: debug_level,&
                    debug_constitutive,&
                    debug_levelBasic
 use lattice

 implicit none
 integer(pInt), intent(in) :: myFile
 integer(pInt), parameter :: maxNchunks = lattice_maxNinteraction + 1_pInt
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt), dimension(6) :: configNchunks
 integer(pInt) :: section, maxNinstance, i,j,k, f,o, &
                  Nchunks_SlipSlip, Nchunks_SlipTwin, Nchunks_TwinSlip, Nchunks_TwinTwin, &
                  Nchunks_SlipFamilies, Nchunks_TwinFamilies, &
                  mySize=0_pInt, myStructure, index_myFamily, index_otherFamily
 character(len=64)   :: tag
 character(len=1024) :: line = ''                                                                   ! to start initialized
 
 write(6,*)
 write(6,*) '<<<+-  constitutive_',trim(constitutive_phenopowerlaw_label),' init  -+>>>'
 write(6,*) '$Id$'
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_plasticity == constitutive_phenopowerlaw_label),pInt)
 if (maxNinstance == 0) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) then
   write(6,'(a16,1x,i5)') '# instances:',maxNinstance
   write(6,*)
 endif

 Nchunks_SlipFamilies = lattice_maxNslipFamily
 Nchunks_TwinFamilies = lattice_maxNtwinFamily
 Nchunks_SlipSlip = lattice_maxNinteraction
 Nchunks_SlipTwin = lattice_maxNinteraction
 Nchunks_TwinSlip = lattice_maxNinteraction
 Nchunks_TwinTwin = lattice_maxNinteraction
 
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
 allocate(constitutive_phenopowerlaw_h0_SlipSlip(maxNinstance))
          constitutive_phenopowerlaw_h0_SlipSlip          = 0.0_pReal
 allocate(constitutive_phenopowerlaw_h0_SlipTwin(maxNinstance))
          constitutive_phenopowerlaw_h0_SlipTwin          = 0.0_pReal
 allocate(constitutive_phenopowerlaw_h0_TwinSlip(maxNinstance))  
          constitutive_phenopowerlaw_h0_TwinSlip          = 0.0_pReal
 allocate(constitutive_phenopowerlaw_h0_TwinTwin(maxNinstance))  
          constitutive_phenopowerlaw_h0_TwinTwin          = 0.0_pReal
 allocate(constitutive_phenopowerlaw_interaction_SlipSlip(lattice_maxNinteraction,maxNinstance))
          constitutive_phenopowerlaw_interaction_SlipSlip = 0.0_pReal
 allocate(constitutive_phenopowerlaw_interaction_SlipTwin(lattice_maxNinteraction,maxNinstance))
          constitutive_phenopowerlaw_interaction_SlipTwin = 0.0_pReal
 allocate(constitutive_phenopowerlaw_interaction_TwinSlip(lattice_maxNinteraction,maxNinstance))
          constitutive_phenopowerlaw_interaction_TwinSlip = 0.0_pReal
 allocate(constitutive_phenopowerlaw_interaction_TwinTwin(lattice_maxNinteraction,maxNinstance))
          constitutive_phenopowerlaw_interaction_TwinTwin = 0.0_pReal
 allocate(constitutive_phenopowerlaw_a_slip(maxNinstance))
          constitutive_phenopowerlaw_a_slip               = 0.0_pReal
 allocate(constitutive_phenopowerlaw_aTolResistance(maxNinstance))
          constitutive_phenopowerlaw_aTolResistance       = 0.0_pReal
 allocate(constitutive_phenopowerlaw_aTolShear(maxNinstance))
          constitutive_phenopowerlaw_aTolShear            = 0.0_pReal
 allocate(constitutive_phenopowerlaw_aTolTwinfrac(maxNinstance))
          constitutive_phenopowerlaw_aTolTwinfrac         = 0.0_pReal
 allocate(constitutive_phenopowerlaw_nonSchmidCoeff(lattice_maxNonSchmid,maxNinstance))
          constitutive_phenopowerlaw_nonSchmidCoeff = 0.0_pReal

 rewind(myFile)
 section = 0_pInt
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= 'phase')                                               ! wind forward to <phase>
   read(myFile,'(a1024)',END=100) line
 enddo

 do                                                                                                 ! read thru sections of phase part
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt                                                                     ! advance section counter
     cycle                                                                                          ! skip to next line
   endif
   if (section > 0_pInt .and. phase_plasticity(section) == constitutive_phenopowerlaw_label) then   ! one of my sections
     i = phase_plasticityInstance(section)                                                          ! which instance of my plasticity is present phase
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('plasticity','elasticity')
         cycle
       case ('(output)')
         constitutive_phenopowerlaw_Noutput(i) = constitutive_phenopowerlaw_Noutput(i) + 1_pInt
         constitutive_phenopowerlaw_output(constitutive_phenopowerlaw_Noutput(i),i) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
       case ('lattice_structure')
         constitutive_phenopowerlaw_structureName(i) = IO_lc(IO_stringValue(line,positions,2_pInt))
         configNchunks = lattice_configNchunks(constitutive_phenopowerlaw_structureName(i))
         Nchunks_SlipFamilies = configNchunks(1)
         Nchunks_TwinFamilies = configNchunks(2)
         Nchunks_SlipSlip =     configNchunks(3)
         Nchunks_SlipTwin =     configNchunks(4)
         Nchunks_TwinSlip =     configNchunks(5)
         Nchunks_TwinTwin =     configNchunks(6)
       case ('covera_ratio')
         constitutive_phenopowerlaw_CoverA(i) = IO_floatValue(line,positions,2_pInt)
       case ('c11')
         constitutive_phenopowerlaw_Cslip_66(1,1,i) = IO_floatValue(line,positions,2_pInt)
       case ('c12')
         constitutive_phenopowerlaw_Cslip_66(1,2,i) = IO_floatValue(line,positions,2_pInt)
       case ('c13')
         constitutive_phenopowerlaw_Cslip_66(1,3,i) = IO_floatValue(line,positions,2_pInt)
       case ('c22')
         constitutive_phenopowerlaw_Cslip_66(2,2,i) = IO_floatValue(line,positions,2_pInt)
       case ('c23')
         constitutive_phenopowerlaw_Cslip_66(2,3,i) = IO_floatValue(line,positions,2_pInt)
       case ('c33')
         constitutive_phenopowerlaw_Cslip_66(3,3,i) = IO_floatValue(line,positions,2_pInt)
       case ('c44')
         constitutive_phenopowerlaw_Cslip_66(4,4,i) = IO_floatValue(line,positions,2_pInt)
       case ('c55')
         constitutive_phenopowerlaw_Cslip_66(5,5,i) = IO_floatValue(line,positions,2_pInt)
       case ('c66')
         constitutive_phenopowerlaw_Cslip_66(6,6,i) = IO_floatValue(line,positions,2_pInt)
       case ('nslip')
         do j = 1_pInt, Nchunks_SlipFamilies
            constitutive_phenopowerlaw_Nslip(j,i) = IO_intValue(line,positions,1_pInt+j)
          enddo
       case ('gdot0_slip')
         constitutive_phenopowerlaw_gdot0_slip(i) = IO_floatValue(line,positions,2_pInt)
       case ('n_slip')
         constitutive_phenopowerlaw_n_slip(i) = IO_floatValue(line,positions,2_pInt)
       case ('tau0_slip')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_phenopowerlaw_tau0_slip(j,i) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('tausat_slip')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_phenopowerlaw_tausat_slip(j,i) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('a_slip', 'w0_slip')
         constitutive_phenopowerlaw_a_slip(i) = IO_floatValue(line,positions,2_pInt)
       case ('ntwin')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_phenopowerlaw_Ntwin(j,i) = IO_intValue(line,positions,1_pInt+j)
         enddo
       case ('gdot0_twin')
         constitutive_phenopowerlaw_gdot0_twin(i) = IO_floatValue(line,positions,2_pInt)
       case ('n_twin')
         constitutive_phenopowerlaw_n_twin(i) = IO_floatValue(line,positions,2_pInt)
       case ('tau0_twin')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_phenopowerlaw_tau0_twin(j,i) = IO_floatValue(line,positions,1_pInt+j)
         enddo
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
         constitutive_phenopowerlaw_h0_SlipSlip(i) = IO_floatValue(line,positions,2_pInt)
       case ('h0_sliptwin')
         constitutive_phenopowerlaw_h0_SlipTwin(i) = IO_floatValue(line,positions,2_pInt)
         call IO_warning(42_pInt,ext_msg=trim(tag)//' ('//constitutive_phenopowerlaw_label//')')
       case ('h0_twinslip')
         constitutive_phenopowerlaw_h0_TwinSlip(i) = IO_floatValue(line,positions,2_pInt)
       case ('h0_twintwin')
         constitutive_phenopowerlaw_h0_TwinTwin(i) = IO_floatValue(line,positions,2_pInt)
       case ('atol_resistance')
         constitutive_phenopowerlaw_aTolResistance(i) = IO_floatValue(line,positions,2_pInt)
       case ('atol_shear')
         constitutive_phenopowerlaw_aTolShear(i)      = IO_floatValue(line,positions,2_pInt)
       case ('atol_twinfrac')
         constitutive_phenopowerlaw_aTolTwinfrac(i)   = IO_floatValue(line,positions,2_pInt)
       case ('interaction_slipslip')
         do j = 1_pInt, Nchunks_SlipSlip
           constitutive_phenopowerlaw_interaction_SlipSlip(j,i) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_sliptwin')
         do j = 1_pInt, Nchunks_SlipTwin
           constitutive_phenopowerlaw_interaction_SlipTwin(j,i) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_twinslip')
         do j = 1_pInt, Nchunks_TwinSlip
           constitutive_phenopowerlaw_interaction_TwinSlip(j,i) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_twintwin')
         do j = 1_pInt, Nchunks_TwinTwin
           constitutive_phenopowerlaw_interaction_TwinTwin(j,i) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('nonschmid_coefficients')
         do j = 1_pInt, lattice_maxNonSchmid
           constitutive_phenopowerlaw_nonSchmidCoeff(j,i) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case default
         call IO_error(210_pInt,ext_msg=tag//' ('//constitutive_phenopowerlaw_label//')')
     end select
   endif
 enddo

100 do i = 1_pInt,maxNinstance

   constitutive_phenopowerlaw_structure(i) = lattice_initializeStructure(constitutive_phenopowerlaw_structureName(i), &    ! get structure
                                                                         constitutive_phenopowerlaw_CoverA(i))
   constitutive_phenopowerlaw_Nslip(1:lattice_maxNslipFamily,i) = &
            min(lattice_NslipSystem(1:lattice_maxNslipFamily,constitutive_phenopowerlaw_structure(i)),& ! limit active slip systems per family to min of available and requested
                constitutive_phenopowerlaw_Nslip(1:lattice_maxNslipFamily,i))
   constitutive_phenopowerlaw_Ntwin(1:lattice_maxNtwinFamily,i) = &
            min(lattice_NtwinSystem(1:lattice_maxNtwinFamily,constitutive_phenopowerlaw_structure(i)),& ! limit active twin systems per family to min of available and requested
                constitutive_phenopowerlaw_Ntwin(:,i))
   constitutive_phenopowerlaw_totalNslip(i) = sum(constitutive_phenopowerlaw_Nslip(:,i))                ! how many slip systems altogether
   constitutive_phenopowerlaw_totalNtwin(i) = sum(constitutive_phenopowerlaw_Ntwin(:,i))                ! how many twin systems altogether

   if (constitutive_phenopowerlaw_structure(i) < 1 )                  call IO_error(205_pInt,e=i)
   if (any(constitutive_phenopowerlaw_tau0_slip(:,i) < 0.0_pReal .and. &
           constitutive_phenopowerlaw_Nslip(:,i) > 0))                call IO_error(211_pInt,e=i,ext_msg='tau0_slip (' &
                                                                           //constitutive_phenopowerlaw_label//')')
   if (constitutive_phenopowerlaw_gdot0_slip(i) <= 0.0_pReal)         call IO_error(211_pInt,e=i,ext_msg='gdot0_slip (' &
                                                                           //constitutive_phenopowerlaw_label//')')
   if (constitutive_phenopowerlaw_n_slip(i) <= 0.0_pReal)             call IO_error(211_pInt,e=i,ext_msg='n_slip (' &
                                                                           //constitutive_phenopowerlaw_label//')')
   if (any(constitutive_phenopowerlaw_tausat_slip(:,i) <= 0.0_pReal .and. &
           constitutive_phenopowerlaw_Nslip(:,i) > 0))                call IO_error(211_pInt,e=i,ext_msg='tausat_slip (' &
                                                                           //constitutive_phenopowerlaw_label//')')
   if (any(constitutive_phenopowerlaw_a_slip(i) == 0.0_pReal .and. &
           constitutive_phenopowerlaw_Nslip(:,i) > 0))                call IO_error(211_pInt,e=i,ext_msg='a_slip (' &
                                                                           //constitutive_phenopowerlaw_label//')')
   if (any(constitutive_phenopowerlaw_tau0_twin(:,i) < 0.0_pReal .and. &
           constitutive_phenopowerlaw_Ntwin(:,i) > 0))                call IO_error(211_pInt,e=i,ext_msg='tau0_twin (' &
                                                                           //constitutive_phenopowerlaw_label//')')
   if (    constitutive_phenopowerlaw_gdot0_twin(i) <= 0.0_pReal .and. &
       any(constitutive_phenopowerlaw_Ntwin(:,i) > 0))                call IO_error(211_pInt,e=i,ext_msg='gdot0_twin (' &
                                                                           //constitutive_phenopowerlaw_label//')')
   if (    constitutive_phenopowerlaw_n_twin(i) <= 0.0_pReal .and. &
       any(constitutive_phenopowerlaw_Ntwin(:,i) > 0))                call IO_error(211_pInt,e=i,ext_msg='n_twin (' &
                                                                           //constitutive_phenopowerlaw_label//')')
   if (constitutive_phenopowerlaw_aTolResistance(i) <= 0.0_pReal) &
     constitutive_phenopowerlaw_aTolResistance(i) = 1.0_pReal              ! default absolute tolerance 1 Pa
   if (constitutive_phenopowerlaw_aTolShear(i) <= 0.0_pReal) &
     constitutive_phenopowerlaw_aTolShear(i) = 1.0e-6_pReal                ! default absolute tolerance 1e-6
   if (constitutive_phenopowerlaw_aTolTwinfrac(i) <= 0.0_pReal) &
     constitutive_phenopowerlaw_aTolTwinfrac(i) = 1.0e-6_pReal             ! default absolute tolerance 1e-6

 enddo

 allocate(constitutive_phenopowerlaw_hardeningMatrix_SlipSlip(maxval(constitutive_phenopowerlaw_totalNslip),&   ! slip resistance from slip activity
                                                              maxval(constitutive_phenopowerlaw_totalNslip),&
                                                              maxNinstance))
 allocate(constitutive_phenopowerlaw_hardeningMatrix_SlipTwin(maxval(constitutive_phenopowerlaw_totalNslip),&   ! slip resistance from twin activity
                                                              maxval(constitutive_phenopowerlaw_totalNtwin),&
                                                              maxNinstance))
 allocate(constitutive_phenopowerlaw_hardeningMatrix_TwinSlip(maxval(constitutive_phenopowerlaw_totalNtwin),&   ! twin resistance from slip activity
                                                              maxval(constitutive_phenopowerlaw_totalNslip),&
                                                              maxNinstance))
 allocate(constitutive_phenopowerlaw_hardeningMatrix_TwinTwin(maxval(constitutive_phenopowerlaw_totalNtwin),&   ! twin resistance from twin activity
                                                              maxval(constitutive_phenopowerlaw_totalNtwin),&
                                                              maxNinstance))
 constitutive_phenopowerlaw_hardeningMatrix_SlipSlip = 0.0_pReal
 constitutive_phenopowerlaw_hardeningMatrix_SlipTwin = 0.0_pReal
 constitutive_phenopowerlaw_hardeningMatrix_TwinSlip = 0.0_pReal
 constitutive_phenopowerlaw_hardeningMatrix_TwinTwin = 0.0_pReal

 do i = 1_pInt,maxNinstance
   do o = 1_pInt,constitutive_phenopowerlaw_Noutput(i)
     select case(constitutive_phenopowerlaw_output(o,i))
       case('resistance_slip', &
            'shearrate_slip', &
            'accumulatedshear_slip', &
            'resolvedstress_slip' &
            )
         mySize = constitutive_phenopowerlaw_totalNslip(i)
       case('resistance_twin', &
            'shearrate_twin', &
            'accumulatedshear_twin', &
            'resolvedstress_twin' &
            )
         mySize = constitutive_phenopowerlaw_totalNtwin(i)
       case('totalshear', &
            'totalvolfrac' &
            )
         mySize = 1_pInt
       case default
         call IO_error(212_pInt,ext_msg=constitutive_phenopowerlaw_output(o,i)//' ('//constitutive_phenopowerlaw_label//')')
     end select

     if (mySize > 0_pInt) then                                                                      ! any meaningful output found
       constitutive_phenopowerlaw_sizePostResult(o,i) = mySize
       constitutive_phenopowerlaw_sizePostResults(i) = &
       constitutive_phenopowerlaw_sizePostResults(i) + mySize
     endif
   enddo    ! outputs

   constitutive_phenopowerlaw_sizeDotState(i) = constitutive_phenopowerlaw_totalNslip(i)+ &
                                                constitutive_phenopowerlaw_totalNtwin(i)+ &
                                                2_pInt + &
                                                constitutive_phenopowerlaw_totalNslip(i)+ &
                                                constitutive_phenopowerlaw_totalNtwin(i)            ! s_slip, s_twin, sum(gamma), sum(f), accshear_slip, accshear_twin
   constitutive_phenopowerlaw_sizeState(i)    = constitutive_phenopowerlaw_sizeDotState(i)

   myStructure = constitutive_phenopowerlaw_structure(i)

   constitutive_phenopowerlaw_Cslip_66(:,:,i) = lattice_symmetrizeC66(constitutive_phenopowerlaw_structureName(i),&
                                                                      constitutive_phenopowerlaw_Cslip_66(:,:,i))   
                                                                                                    ! assign elasticity tensor
   constitutive_phenopowerlaw_Cslip_66(:,:,i) = &
     math_Mandel3333to66(math_Voigt66to3333(constitutive_phenopowerlaw_Cslip_66(:,:,i)))

   do f = 1_pInt,lattice_maxNslipFamily                                                             ! >>> interaction slip -- X
     index_myFamily = sum(constitutive_phenopowerlaw_Nslip(1:f-1_pInt,i))
     do j = 1_pInt,constitutive_phenopowerlaw_Nslip(f,i)                                            ! loop over (active) systems in my family (slip)
       do o = 1_pInt,lattice_maxNslipFamily
         index_otherFamily = sum(constitutive_phenopowerlaw_Nslip(1:o-1_pInt,i))
         do k = 1_pInt,constitutive_phenopowerlaw_Nslip(o,i)                                        ! loop over (active) systems in other family (slip)
           constitutive_phenopowerlaw_hardeningMatrix_SlipSlip(index_myFamily+j,index_otherFamily+k,i) = &
               constitutive_phenopowerlaw_interaction_SlipSlip(lattice_interactionSlipSlip( &
                                                                 sum(lattice_NslipSystem(1:f-1,myStructure))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1,myStructure))+k, &
                                                                 myStructure), i )
       enddo; enddo

       do o = 1_pInt,lattice_maxNtwinFamily
         index_otherFamily = sum(constitutive_phenopowerlaw_Ntwin(1:o-1_pInt,i))
         do k = 1_pInt,constitutive_phenopowerlaw_Ntwin(o,i)                                        ! loop over (active) systems in other family (twin)
           constitutive_phenopowerlaw_hardeningMatrix_SlipTwin(index_myFamily+j,index_otherFamily+k,i) = &
               constitutive_phenopowerlaw_interaction_SlipTwin(lattice_interactionSlipTwin( &
                                                                 sum(lattice_NslipSystem(1:f-1_pInt,myStructure))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,myStructure))+k, &
                                                                 myStructure), i )
       enddo; enddo

   enddo; enddo

   do f = 1_pInt,lattice_maxNtwinFamily                                                             ! >>> interaction twin -- X
     index_myFamily = sum(constitutive_phenopowerlaw_Ntwin(1:f-1_pInt,i))
     do j = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,i)                                            ! loop over (active) systems in my family (twin)

       do o = 1_pInt,lattice_maxNslipFamily
         index_otherFamily = sum(constitutive_phenopowerlaw_Nslip(1:o-1_pInt,i))
         do k = 1_pInt,constitutive_phenopowerlaw_Nslip(o,i)                                        ! loop over (active) systems in other family (slip)
           constitutive_phenopowerlaw_hardeningMatrix_TwinSlip(index_myFamily+j,index_otherFamily+k,i) = &
               constitutive_phenopowerlaw_interaction_TwinSlip(lattice_interactionTwinSlip( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,myStructure))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1_pInt,myStructure))+k, &
                                                                 myStructure), i )
       enddo; enddo

       do o = 1_pInt,lattice_maxNtwinFamily
         index_otherFamily = sum(constitutive_phenopowerlaw_Ntwin(1:o-1_pInt,i))
         do k = 1_pInt,constitutive_phenopowerlaw_Ntwin(o,i)                                        ! loop over (active) systems in other family (twin)
           constitutive_phenopowerlaw_hardeningMatrix_TwinTwin(index_myFamily+j,index_otherFamily+k,i) = &
               constitutive_phenopowerlaw_interaction_TwinTwin(lattice_interactionTwinTwin( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,myStructure))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,myStructure))+k, &
                                                                 myStructure), i )
       enddo; enddo

   enddo; enddo

! report to out file...

 enddo

end subroutine constitutive_phenopowerlaw_init


!--------------------------------------------------------------------------------------------------
!> @brief initial microstructural state
!--------------------------------------------------------------------------------------------------
function constitutive_phenopowerlaw_stateInit(myInstance)
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

end function constitutive_phenopowerlaw_stateInit


!--------------------------------------------------------------------------------------------------
!> @brief absolute state tolerance
!--------------------------------------------------------------------------------------------------
pure function constitutive_phenopowerlaw_aTolState(myInstance)
 implicit none
 integer(pInt), intent(in) ::  myInstance                                                           ! number specifying the current instance of the plasticity
 real(pReal), dimension(constitutive_phenopowerlaw_sizeState(myInstance)) :: &
                              constitutive_phenopowerlaw_aTolState                                  ! relevant state values for the current instance of this plasticity


constitutive_phenopowerlaw_aTolState(1:constitutive_phenopowerlaw_totalNslip(myInstance)+ &
                                       constitutive_phenopowerlaw_totalNtwin(myInstance)) = &
  constitutive_phenopowerlaw_aTolResistance(myInstance)
constitutive_phenopowerlaw_aTolState(1+constitutive_phenopowerlaw_totalNslip(myInstance)+ &
                                       constitutive_phenopowerlaw_totalNtwin(myInstance)) = &
  constitutive_phenopowerlaw_aTolShear(myInstance)
constitutive_phenopowerlaw_aTolState(2+constitutive_phenopowerlaw_totalNslip(myInstance)+ &
                                       constitutive_phenopowerlaw_totalNtwin(myInstance)) = &
  constitutive_phenopowerlaw_aTolTwinFrac(myInstance)
constitutive_phenopowerlaw_aTolState(3+constitutive_phenopowerlaw_totalNslip(myInstance)+ &
                                       constitutive_phenopowerlaw_totalNtwin(myInstance): &
                                     2+2*(constitutive_phenopowerlaw_totalNslip(myInstance)+ &
                                          constitutive_phenopowerlaw_totalNtwin(myInstance))) = &
  constitutive_phenopowerlaw_aTolShear(myInstance)

end function constitutive_phenopowerlaw_aTolState


!--------------------------------------------------------------------------------------------------
!> @brief homogenized elacticity matrix
!--------------------------------------------------------------------------------------------------
pure function constitutive_phenopowerlaw_homogenizedC(state,ipc,ip,el)
 use prec, only: p_vec
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_plasticityInstance
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !component-ID of current integration point
   ip, &                                                                                            !current integration point
   el                                                                                               !current element
 integer(pInt) matID
 real(pReal), dimension(6,6) :: constitutive_phenopowerlaw_homogenizedC
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
   state                                                                                            ! state variables
 
 matID = phase_plasticityInstance(material_phase(ipc,ip,el))
 constitutive_phenopowerlaw_homogenizedC = constitutive_phenopowerlaw_Cslip_66(1:6,1:6,matID)

end function constitutive_phenopowerlaw_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief calculate derived quantities from state (dummy subroutine, not used here) 
!--------------------------------------------------------------------------------------------------
pure subroutine constitutive_phenopowerlaw_microstructure(Temperature,state,ipc,ip,el)
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
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !component-ID of current integration point
   ip, &                                                                                            !current integration point
   el                                                                                               !current element
 integer(pInt) :: matID
 real(pReal), intent(in) :: Temperature                                                             ! temperature
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state

 matID = phase_plasticityInstance(material_phase(ipc,ip,el))
  
end subroutine constitutive_phenopowerlaw_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine constitutive_phenopowerlaw_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,state,ipc,ip,el)
 use prec, only: p_vec
 use math, only: math_Plain3333to99,math_Mandel6to33
 use lattice, only: lattice_Sslip,lattice_Sslip_v,lattice_Stwin,lattice_Stwin_v, lattice_maxNslipFamily, lattice_maxNtwinFamily, &
                    lattice_NslipSystem,lattice_NtwinSystem,NnonSchmid
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_plasticityInstance

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           ! component-ID at current integration point
   ip, &                                                                                            ! current integration point
   el                                                                                               ! current element
 integer(pInt) matID,nSlip,nTwin,f,i,j,k,l,m,n, structID,index_Gamma,index_F,index_myFamily
 real(pReal) Temperature
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6),    intent(in) :: Tstar_v                                                ! 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal), dimension(3,3), intent(out) :: Lp                                                     ! plastic velocity gradient
 real(pReal), dimension(3,3,3,3)          :: dLp_dTstar3333                                         ! derivative of Lp (4th-rank tensor)
 real(pReal), dimension(9,9), intent(out) :: dLp_dTstar
 real(pReal), dimension(3,3,2) :: nonSchmid_tensor
 real(pReal), dimension(constitutive_phenopowerlaw_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip_pos,gdot_slip_neg,dgdot_dtauslip_pos,dgdot_dtauslip_neg,tau_slip_pos,tau_slip_neg
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
 do f = 1_pInt,lattice_maxNslipFamily                                                               ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,structID))                                   ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,matID)                                          ! process each (active) slip system in family
     j = j+1_pInt
     
!--------------------------------------------------------------------------------------------------
! Calculation of Lp
     tau_slip_pos(j)  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,structID))
     tau_slip_neg(j)  = tau_slip_pos(j)
     nonSchmid_tensor(1:3,1:3,1)  = math_Mandel6to33(lattice_Sslip_v(1:6,1,index_myFamily+i,structID))
     nonSchmid_tensor(1:3,1:3,2)  = nonSchmid_tensor(1:3,1:3,1)
     do k = 1, NnonSchmid(structID) 
       tau_slip_pos(j) = tau_slip_pos(j) + constitutive_phenopowerlaw_nonSchmidCoeff(k,matID)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,structID))
       tau_slip_neg(j) = tau_slip_neg(j) + constitutive_phenopowerlaw_nonSchmidCoeff(k,matID)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,structID))
       nonSchmid_tensor(1:3,1:3,1) = nonSchmid_tensor(1:3,1:3,1) + constitutive_phenopowerlaw_nonSchmidCoeff(k,matID)*&
                                           math_Mandel6to33(lattice_Sslip_v(1:6,2*k,index_myFamily+i,structID))
       nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,2) + constitutive_phenopowerlaw_nonSchmidCoeff(k,matID)*&
                                           math_Mandel6to33(lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,structID))
     enddo
     gdot_slip_pos(j) = 0.5_pReal*constitutive_phenopowerlaw_gdot0_slip(matID)* &
                    ((abs(tau_slip_pos(j))/state(ipc,ip,el)%p(j))**constitutive_phenopowerlaw_n_slip(matID))*&
                                                                    sign(1.0_pReal,tau_slip_pos(j))
     gdot_slip_neg(j) = 0.5_pReal*constitutive_phenopowerlaw_gdot0_slip(matID)* &
                    ((abs(tau_slip_neg(j))/state(ipc,ip,el)%p(j))**constitutive_phenopowerlaw_n_slip(matID))*&
                                                                    sign(1.0_pReal,tau_slip_neg(j))
     Lp = Lp + (1.0_pReal-state(ipc,ip,el)%p(index_F))*&                     ! 1-F
               (gdot_slip_pos(j)+gdot_slip_neg(j))*lattice_Sslip(1:3,1:3,index_myFamily+i,structID)

!--------------------------------------------------------------------------------------------------
! Calculation of the tangent of Lp
     if (gdot_slip_pos(j) /= 0.0_pReal) then
       dgdot_dtauslip_pos(j) = gdot_slip_pos(j)*constitutive_phenopowerlaw_n_slip(matID)/tau_slip_pos(j)
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtauslip_pos(j)*lattice_Sslip(k,l,index_myFamily+i,structID)* &
                                                     nonSchmid_tensor(m,n,1)
     endif
     
     if (gdot_slip_neg(j) /= 0.0_pReal) then
       dgdot_dtauslip_neg(j) = gdot_slip_neg(j)*constitutive_phenopowerlaw_n_slip(matID)/tau_slip_neg(j)
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtauslip_neg(j)*lattice_Sslip(k,l,index_myFamily+i,structID)* &
                                                     nonSchmid_tensor(m,n,2)
     endif
   enddo
 enddo

 j = 0_pInt
 do f = 1_pInt,lattice_maxNtwinFamily                                                               ! loop over all twin families
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                                   ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,matID)                                          ! process each (active) twin system in family
     j = j+1_pInt

!--------------------------------------------------------------------------------------------------
! Calculation of Lp
     tau_twin(j)  = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,structID)) 
     gdot_twin(j) = (1.0_pReal-state(ipc,ip,el)%p(index_F))*&                                       ! 1-F
                    constitutive_phenopowerlaw_gdot0_twin(matID)*&
                    (abs(tau_twin(j))/state(ipc,ip,el)%p(nSlip+j))**&
                    constitutive_phenopowerlaw_n_twin(matID)*max(0.0_pReal,sign(1.0_pReal,tau_twin(j)))
     Lp = Lp + gdot_twin(j)*lattice_Stwin(1:3,1:3,index_myFamily+i,structID)

!--------------------------------------------------------------------------------------------------
! Calculation of the tangent of Lp
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

end subroutine constitutive_phenopowerlaw_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief of change of microstructure, evolution of state variable
!--------------------------------------------------------------------------------------------------
function constitutive_phenopowerlaw_dotState(Tstar_v,Temperature,state,ipc,ip,el)
 use prec,     only: p_vec
 use lattice,  only: lattice_Sslip_v, lattice_Stwin_v, lattice_maxNslipFamily, lattice_maxNtwinFamily, &
                     lattice_NslipSystem,lattice_NtwinSystem,lattice_shearTwin,NnonSchmid   
 use mesh,     only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_plasticityInstance
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID at current integration point
   ip, &                                                                                            !< current integration point
   el                                                                                               !< current element
 integer(pInt) matID,nSlip,nTwin,f,i,j,k,structID, &
               index_Gamma,index_F,offset_accshear_slip,offset_accshear_twin,index_myFamily 
 real(pReal) Temperature,c_SlipSlip,c_SlipTwin,c_TwinSlip,c_TwinTwin, ssat_offset
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
 real(pReal), dimension(6), intent(in) :: Tstar_v                                                   !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal), dimension(constitutive_phenopowerlaw_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip,tau_slip_pos,tau_slip_neg,left_SlipSlip,left_SlipTwin,right_SlipSlip,right_TwinSlip
 real(pReal), dimension(constitutive_phenopowerlaw_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_twin,tau_twin,left_TwinSlip,left_TwinTwin,right_SlipTwin,right_TwinTwin
 real(pReal), dimension(constitutive_phenopowerlaw_sizeDotState(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_phenopowerlaw_dotState

 matID = phase_plasticityInstance(material_phase(ipc,ip,el))
 structID = constitutive_phenopowerlaw_structure(matID)
 
 nSlip = constitutive_phenopowerlaw_totalNslip(matID)
 nTwin = constitutive_phenopowerlaw_totalNtwin(matID)

 index_Gamma = nSlip + nTwin + 1_pInt
 index_F     = nSlip + nTwin + 2_pInt
 offset_accshear_slip = nSlip + nTwin + 2_pInt
 offset_accshear_twin = nSlip + nTwin + 2_pInt + nSlip

 constitutive_phenopowerlaw_dotState = 0.0_pReal
 
!--------------------------------------------------------------------------------------------------
! system-independent (nonlinear) prefactors to M_Xx (X influenced by x) matrices
 c_SlipSlip = constitutive_phenopowerlaw_h0_SlipSlip(matID)*&
              (1.0_pReal + &
               constitutive_phenopowerlaw_twinC(matID)*state(ipc,ip,el)%p(index_F)**constitutive_phenopowerlaw_twinB(matID))
 c_SlipTwin = 0.0_pReal
 c_TwinSlip = constitutive_phenopowerlaw_h0_TwinSlip(matID)*&
              state(ipc,ip,el)%p(index_Gamma)**constitutive_phenopowerlaw_twinE(matID)
 c_TwinTwin = constitutive_phenopowerlaw_h0_TwinTwin(matID)*&
              state(ipc,ip,el)%p(index_F)**constitutive_phenopowerlaw_twinD(matID)

!-- calculate left and right vectors and calculate dot gammas

 ssat_offset = constitutive_phenopowerlaw_spr(matID)*sqrt(state(ipc,ip,el)%p(index_F))
 j = 0_pInt
 do f = 1_pInt,lattice_maxNslipFamily                                                               ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,structID))                                   ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,matID)                                          ! process each (active) slip system in family
     j = j+1_pInt
     left_SlipSlip(j) = 1.0_pReal                                                 ! no system-dependent left part
     left_SlipTwin(j) = 1.0_pReal                                                 ! no system-dependent left part
     right_SlipSlip(j) = abs(1.0_pReal-state(ipc,ip,el)%p(j) / &
                                    (constitutive_phenopowerlaw_tausat_slip(f,matID)+ssat_offset)) &
                         **constitutive_phenopowerlaw_a_slip(matID)&
                         *sign(1.0_pReal,1.0_pReal-state(ipc,ip,el)%p(j) / &
                                    (constitutive_phenopowerlaw_tausat_slip(f,matID)+ssat_offset))
     right_TwinSlip(j) = 1.0_pReal                                                ! no system-dependent part
     
!--------------------------------------------------------------------------------------------------
! Calculation of dot gamma 
     tau_slip_pos(j)  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,structID))
     tau_slip_neg(j)  = tau_slip_pos(j)
     do k = 1, NnonSchmid(structID) 
       tau_slip_pos(j) = tau_slip_pos(j) + constitutive_phenopowerlaw_nonSchmidCoeff(k,matID)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,structID))
       tau_slip_neg(j) = tau_slip_neg(j) + constitutive_phenopowerlaw_nonSchmidCoeff(k,matID)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,structID))
     enddo
     gdot_slip(j) = constitutive_phenopowerlaw_gdot0_slip(matID)*0.5_pReal* &
                    ((abs(tau_slip_pos(j))/state(ipc,ip,el)%p(j))**constitutive_phenopowerlaw_n_slip(matID) &
                    +(abs(tau_slip_neg(j))/state(ipc,ip,el)%p(j))**constitutive_phenopowerlaw_n_slip(matID))&
                    *sign(1.0_pReal,tau_slip_pos(j)) 
    enddo
  enddo

 j = 0_pInt
 do f = 1_pInt,lattice_maxNtwinFamily                                                               ! loop over all twin families
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                                   ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,matID)                                          ! process each (active) twin system in family
     j = j+1_pInt
     left_TwinSlip(j)  = 1.0_pReal                                                                  ! no system-dependent right part
     left_TwinTwin(j)  = 1.0_pReal                                                                  ! no system-dependent right part
     right_SlipTwin(j) = 1.0_pReal                                                                  ! no system-dependent right part
     right_TwinTwin(j) = 1.0_pReal                                                                  ! no system-dependent right part

!* Calculation of dot vol frac

     tau_twin(j)  = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,structID)) 
     gdot_twin(j) = (1.0_pReal-state(ipc,ip,el)%p(index_F))*&                                       ! 1-F
                    constitutive_phenopowerlaw_gdot0_twin(matID)*&
                    (abs(tau_twin(j))/state(ipc,ip,el)%p(nSlip+j))**&
                    constitutive_phenopowerlaw_n_twin(matID)*max(0.0_pReal,sign(1.0_pReal,tau_twin(j)))
    enddo
  enddo

!--------------------------------------------------------------------------------------------------
! calculate the overall hardening based on above
 j = 0_pInt
 do f = 1_pInt,lattice_maxNslipFamily                                                               ! loop over all slip families
   do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,matID)                                          ! process each (active) slip system in family
     j = j+1_pInt
     constitutive_phenopowerlaw_dotState(j) = &                                                     ! evolution of slip resistance j
       c_SlipSlip * left_SlipSlip(j) * &
       dot_product(constitutive_phenopowerlaw_hardeningMatrix_SlipSlip(j,1:nSlip,matID), &
                   right_SlipSlip*abs(gdot_slip)) + &                                               ! dot gamma_slip modulated by right-side slip factor
       c_SlipTwin * left_SlipTwin(j) * &
       dot_product(constitutive_phenopowerlaw_hardeningMatrix_SlipTwin(j,1:nTwin,matID), &
                   right_SlipTwin*gdot_twin)                                                        ! dot gamma_twin modulated by right-side twin factor
     constitutive_phenopowerlaw_dotState(index_Gamma) = constitutive_phenopowerlaw_dotState(index_Gamma) + &
                                                        abs(gdot_slip(j))
     constitutive_phenopowerlaw_dotState(offset_accshear_slip+j) = abs(gdot_slip(j))
   enddo
 enddo
 
 j = 0_pInt
 do f = 1_pInt,lattice_maxNtwinFamily                                                               ! loop over all twin families
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                                   ! at which index starts my family
   do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,matID)                                          ! process each (active) twin system in family
     j = j+1_pInt
     constitutive_phenopowerlaw_dotState(j+nSlip) = &                                               ! evolution of twin resistance j
       c_TwinSlip * left_TwinSlip(j) * &
       dot_product(constitutive_phenopowerlaw_hardeningMatrix_TwinSlip(j,1:nSlip,matID), &
                   right_TwinSlip*abs(gdot_slip)) + &                                               ! dot gamma_slip modulated by right-side slip factor
       c_TwinTwin * left_TwinTwin(j) * &
       dot_product(constitutive_phenopowerlaw_hardeningMatrix_TwinTwin(j,1:nTwin,matID), &
                   right_TwinTwin*gdot_twin)                                                        ! dot gamma_twin modulated by right-side twin factor
     constitutive_phenopowerlaw_dotState(index_F) = constitutive_phenopowerlaw_dotState(index_F) + &
                                                    gdot_twin(j)/lattice_shearTwin(index_myFamily+i,structID)
     constitutive_phenopowerlaw_dotState(offset_accshear_twin+j) = abs(gdot_twin(j))
   enddo
 enddo

end function constitutive_phenopowerlaw_dotState


!--------------------------------------------------------------------------------------------------
!> @brief (instantaneous) incremental change of microstructure
!--------------------------------------------------------------------------------------------------
function constitutive_phenopowerlaw_deltaState(Tstar_v, Temperature, state, g,ip,el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_plasticityInstance

implicit none

integer(pInt), intent(in) ::                g, &                                                    ! current grain number
                                            ip, &                                                   ! current integration point
                                            el                                                      ! current element number
real(pReal), intent(in) ::                  Temperature                                             ! temperature
real(pReal), dimension(6), intent(in) ::    Tstar_v                                                 ! current 2nd Piola-Kirchhoff stress in Mandel notation
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            state                                                   ! current microstructural state

real(pReal), dimension(constitutive_phenopowerlaw_sizeDotState(phase_plasticityInstance(material_phase(g,ip,el)))) :: &
                                            constitutive_phenopowerlaw_deltaState                   ! change of state variables / microstructure


constitutive_phenopowerlaw_deltaState = 0.0_pReal

end function constitutive_phenopowerlaw_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of temperature (dummy function)
!--------------------------------------------------------------------------------------------------
pure function constitutive_phenopowerlaw_dotTemperature(Tstar_v,Temperature,state,ipc,ip,el)
  use prec,     only: pReal,pInt,p_vec
  use mesh,     only: mesh_NcpElems, mesh_maxNips
  use material, only: homogenization_maxNgrains
  
  implicit none
  real(pReal), dimension(6), intent(in) ::  Tstar_v                                                 ! 2nd Piola Kirchhoff stress tensor in Mandel notation
  real(pReal), intent(in) ::                Temperature
  integer(pInt), intent(in)::               ipc, &                                                  ! grain number
                                            ip, &                                                   ! integration point number
                                            el                                                      ! element number
  type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state ! state of the current microstructure
  
  real(pReal) constitutive_phenopowerlaw_dotTemperature                                             ! rate of change of temparature
  
  constitutive_phenopowerlaw_dotTemperature = 0.0_pReal

end function constitutive_phenopowerlaw_dotTemperature


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
pure function constitutive_phenopowerlaw_postResults(Tstar_v,Temperature,dt,state,ipc,ip,el)
 use prec, only: pReal,pInt,p_vec
 use lattice, only: lattice_Sslip_v,lattice_Stwin_v, lattice_maxNslipFamily, lattice_maxNtwinFamily, &
                    lattice_NslipSystem,lattice_NtwinSystem,NnonSchmid 
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase,phase_plasticityInstance,phase_Noutput
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !component-ID at current integration point
   ip, &                                                                                            !current integration point
   el                                                                                               !current element
 real(pReal), intent(in) :: &
   dt, &                                                                                            !current time increment
   Temperature
 real(pReal), dimension(6), intent(in) :: Tstar_v                                                   ! 2nd Piola Kirchhoff stress tensor (Mandel)
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
 integer(pInt) matID,o,f,i,c,nSlip,nTwin,j,k,structID, &
               index_Gamma,index_F,index_accshear_slip,index_accshear_twin,index_myFamily 
 real(pReal) tau_slip_pos,tau_slip_neg,tau
 real(pReal), dimension(constitutive_phenopowerlaw_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_phenopowerlaw_postResults

 matID = phase_plasticityInstance(material_phase(ipc,ip,el))
 structID = constitutive_phenopowerlaw_structure(matID)

 nSlip = constitutive_phenopowerlaw_totalNslip(matID)
 nTwin = constitutive_phenopowerlaw_totalNtwin(matID)

 index_Gamma = nSlip + nTwin + 1_pInt
 index_F     = nSlip + nTwin + 2_pInt
 index_accshear_slip = nSlip + nTwin + 3_pInt
 index_accshear_twin = nSlip + nTwin + 3_pInt + nSlip

 constitutive_phenopowerlaw_postResults = 0.0_pReal
 c = 0_pInt

 do o = 1_pInt,phase_Noutput(material_phase(ipc,ip,el))
   select case(constitutive_phenopowerlaw_output(o,matID))
     case ('resistance_slip')
       constitutive_phenopowerlaw_postResults(c+1_pInt:c+nSlip) = state(ipc,ip,el)%p(1:nSlip)
       c = c + nSlip

     case ('accumulatedshear_slip')
       constitutive_phenopowerlaw_postResults(c+1_pInt:c+nSlip) = state(ipc,ip,el)%p(index_accshear_slip:&
                                                                                     index_accshear_slip+nSlip)
       c = c + nSlip

     case ('shearrate_slip')
       j = 0_pInt
       do f = 1_pInt,lattice_maxNslipFamily                                                         ! loop over all slip families
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,structID))                             ! at which index starts my family
         do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,matID)                                    ! process each (active) slip system in family
           j = j + 1_pInt
           tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,structID))
           tau_slip_neg  = tau_slip_pos
           do k = 1, NnonSchmid(structID) 
             tau_slip_pos = tau_slip_pos + constitutive_phenopowerlaw_nonSchmidCoeff(k,matID)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,structID))
             tau_slip_neg = tau_slip_neg + constitutive_phenopowerlaw_nonSchmidCoeff(k,matID)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,structID))
           enddo
           constitutive_phenopowerlaw_postResults(c+j) = constitutive_phenopowerlaw_gdot0_slip(matID)*0.5_pReal* &
                    ((abs(tau_slip_pos)/state(ipc,ip,el)%p(j))**constitutive_phenopowerlaw_n_slip(matID) &
                    +(abs(tau_slip_neg)/state(ipc,ip,el)%p(j))**constitutive_phenopowerlaw_n_slip(matID))&
                    *sign(1.0_pReal,tau_slip_pos)
       enddo; enddo
       c = c + nSlip

     case ('resolvedstress_slip')
       j = 0_pInt
       do f = 1_pInt,lattice_maxNslipFamily                                                         ! loop over all slip families
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,structID))                             ! at which index starts my family
         do i = 1_pInt,constitutive_phenopowerlaw_Nslip(f,matID)                                    ! process each (active) slip system in family
           j = j + 1_pInt
           constitutive_phenopowerlaw_postResults(c+j) = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,structID))
       enddo; enddo
       c = c + nSlip

     case ('totalshear')
       constitutive_phenopowerlaw_postResults(c+1_pInt) = state(ipc,ip,el)%p(index_Gamma)
       c = c + 1_pInt

     case ('resistance_twin')
       constitutive_phenopowerlaw_postResults(c+1_pInt:c+nTwin) = state(ipc,ip,el)%p(1_pInt+nSlip:nTwin+nSlip)
       c = c + nTwin

     case ('accumulatedshear_twin')
       constitutive_phenopowerlaw_postResults(c+1_pInt:c+nTwin) = state(ipc,ip,el)%p(index_accshear_twin:&
                                                                                     index_accshear_twin+nTwin)
       c = c + nTwin

     case ('shearrate_twin')
       j = 0_pInt
       do f = 1_pInt,lattice_maxNtwinFamily                                                         ! loop over all twin families
         index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                             ! at which index starts my family
         do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,matID)                                    ! process each (active) twin system in family
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
       do f = 1_pInt,lattice_maxNtwinFamily                                                         ! loop over all twin families
         index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,structID))                             ! at which index starts my family
         do i = 1_pInt,constitutive_phenopowerlaw_Ntwin(f,matID)                                    ! process each (active) twin system in family
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
