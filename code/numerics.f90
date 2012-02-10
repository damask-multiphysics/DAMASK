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
!##############################################################
MODULE numerics
!##############################################################

use prec, only: pInt, pReal
use IO, only: IO_warning
implicit none

character(len=64), parameter :: numerics_configFile        = 'numerics.config'               ! name of configuration file
integer(pInt) ::                iJacoStiffness             =  1_pInt, &                      ! frequency of stiffness update
                                iJacoLpresiduum            =  1_pInt, &                      ! frequency of Jacobian update of residuum in Lp
                                nHomog                     = 20_pInt, &                      ! homogenization loop limit (only for debugging info, loop limit is determined by "subStepMinHomog")
                                nMPstate                   = 10_pInt, &                      ! materialpoint state loop limit
                                nCryst                     = 20_pInt, &                      ! crystallite loop limit (only for debugging info, loop limit is determined by "subStepMinCryst")
                                nState                     = 10_pInt, &                      ! state loop limit
                                nStress                    = 40_pInt, &                      ! stress loop limit
                                pert_method                =  1_pInt, &                      ! method used in perturbation technique for tangent
                                numerics_integrationMode   =  0_pInt                         ! integrationMode 1 = central solution ; integrationMode 2 = perturbation, Default 0: undefined, is not read from file
integer(pInt), dimension(2) ::  numerics_integrator        =  1_pInt                         ! method used for state integration (central & perturbed state), Default 1: fix-point iteration for both states
real(pReal) ::                  relevantStrain             =  1.0e-7_pReal, &                ! strain increment considered significant (used by crystallite to determine whether strain inc is considered significant)
                                defgradTolerance           =  1.0e-7_pReal, &                ! deviation of deformation gradient that is still allowed (used by CPFEM to determine outdated ffn1)
                                pert_Fg                    =  1.0e-7_pReal, &                ! strain perturbation for FEM Jacobi
                                subStepMinCryst            =  1.0e-3_pReal, &                ! minimum (relative) size of sub-step allowed during cutback in crystallite
                                subStepMinHomog            =  1.0e-3_pReal, &                ! minimum (relative) size of sub-step allowed during cutback in homogenization
                                subStepSizeCryst           =  0.25_pReal, &                  ! size of first substep when cutback in crystallite
                                subStepSizeHomog           =  0.25_pReal, &                  ! size of first substep when cutback in homogenization
                                stepIncreaseCryst          =  1.5_pReal, &                   ! increase of next substep size when previous substep converged in crystallite
                                stepIncreaseHomog          =  1.5_pReal, &                   ! increase of next substep size when previous substep converged in homogenization
                                rTol_crystalliteState      =  1.0e-6_pReal, &                ! relative tolerance in crystallite state loop 
                                rTol_crystalliteTemperature=  1.0e-6_pReal, &                ! relative tolerance in crystallite temperature loop 
                                rTol_crystalliteStress     =  1.0e-6_pReal, &                ! relative tolerance in crystallite stress loop
                                aTol_crystalliteStress     =  1.0e-8_pReal, &                ! absolute tolerance in crystallite stress loop, Default 1.0e-8: residuum is in Lp and hence strain is on this order

                                absTol_RGC                 =  1.0e+4_pReal, &                ! absolute tolerance of RGC residuum
                                relTol_RGC                 =  1.0e-3_pReal, &                ! relative tolerance of RGC residuum
                                absMax_RGC                 =  1.0e+10_pReal, &               ! absolute maximum of RGC residuum
                                relMax_RGC                 =  1.0e+2_pReal, &                ! relative maximum of RGC residuum
                                pPert_RGC                  =  1.0e-7_pReal, &                ! perturbation for computing RGC penalty tangent
                                xSmoo_RGC                  =  1.0e-5_pReal, &                ! RGC penalty smoothing parameter (hyperbolic tangent)
                                viscPower_RGC              =  1.0e+0_pReal, &                ! power (sensitivity rate) of numerical viscosity in RGC scheme, Default 1.0e0: Newton viscosity (linear model)
                                viscModus_RGC              =  0.0e+0_pReal, &                ! stress modulus of RGC numerical viscosity, Default 0.0e0: No viscosity is applied
                                refRelaxRate_RGC           =  1.0e-3_pReal, &                ! reference relaxation rate in RGC viscosity
                                maxdRelax_RGC              =  1.0e+0_pReal, &                ! threshold of maximum relaxation vector increment (if exceed this then cutback)
                                maxVolDiscr_RGC            =  1.0e-5_pReal, &                ! threshold of maximum volume discrepancy allowed
                                volDiscrMod_RGC            =  1.0e+12_pReal, &               ! stiffness of RGC volume discrepancy (zero = without volume discrepancy constraint)
                                volDiscrPow_RGC            =  5.0_pReal, &                   ! powerlaw penalty for volume discrepancy
!* spectral parameters:
                                err_div_tol                =  1.0e-4_pReal, &                ! error of divergence in fourier space, Default 1.0e-4: proposed by Suquet
                                err_stress_tolrel          =  0.01_pReal , &                 ! relative tolerance for fullfillment of stress BC, Default: 0.01 allowing deviation of 1% of maximum stress 
                                fftw_timelimit             = -1.0_pReal, &                   ! sets the timelimit of plan creation for FFTW, see manual on www.fftw.org, Default -1.0: disable timelimit
                                rotation_tol               =  1.0e-12_pReal                  ! tolerance of rotation specified in loadcase, Default 1.0e-12: first guess
character(len=64) ::            fftw_planner_string        = 'FFTW_PATIENT'                  ! reads the planing-rigor flag, see manual on www.fftw.org, Default FFTW_PATIENT: use patiant planner flag
integer(pInt) ::                fftw_planner_flag          =  -1_pInt                         ! conversion of fftw_planner_string to integer, basically what is usually done in the include file of fftw
logical ::                      memory_efficient           = .true. ,&                       ! for fast execution (pre calculation of gamma_hat), Default .true.: do not precalculate
                                divergence_correction      = .false.     ,&                  ! correct divergence calculation in fourier space, Default .false.: no correction
                                update_gamma               = .false.,&                       ! update gamma operator with current stiffness, Default .false.: use initial stiffness 
                                simplified_algorithm       = .true.                          ! use short algorithm without fluctuation field, Default .true.: use simplified algorithm
real(pReal) ::                  cut_off_value              =  0.0_pReal                      ! percentage of frequencies to cut away, Default 0.0: use all frequencies
integer(pInt) ::                itmax                      = 20_pInt , &                     ! maximum number of iterations


!* Random seeding parameters
                                fixedSeed                  = 0_pInt                          ! fixed seeding for pseudo-random number generator, Default 0: use random seed
!* OpenMP variable
integer(pInt) ::                DAMASK_NumThreadsInt       = 0_pInt                          ! value stored in environment variable DAMASK_NUM_THREADS, set to zero if no OpenMP directive


CONTAINS
 
!*******************************************
!    initialization subroutine
!*******************************************
subroutine numerics_init()
  
use, intrinsic :: iso_fortran_env      
  !*** variables and functions from other modules ***!
  use prec, only:                             pInt, & 
                                              pReal  
  use IO, only:                               IO_error, &
                                              IO_open_file, &
                                              IO_isBlank, &
                                              IO_stringPos, &
                                              IO_stringValue, &
                                              IO_lc, &
                                              IO_floatValue, &
                                              IO_intValue
!$ use OMP_LIB                                                                ! the openMP function library
  
  implicit none
  
  !*** local variables ***!
  integer(pInt), parameter ::                 fileunit = 300_pInt
  integer(pInt), parameter ::                 maxNchunks = 2_pInt
  integer(pInt) ::                            gotDAMASK_NUM_THREADS = 1_pInt
  integer(pInt), dimension(1+2*maxNchunks) :: positions
  character(len=64) ::                        tag
  character(len=1024) ::                      line
  
! OpenMP variable
!$ character(len=6) DAMASK_NumThreadsString                               !environment variable DAMASK_NUM_THREADS

!$OMP CRITICAL (write2out)
  write(6,*)
  write(6,*) '<<<+-  numerics init  -+>>>'
  write(6,*) '$Id$'
#include "compilation_info.f90"
!$OMP END CRITICAL (write2out)

!$ call GET_ENVIRONMENT_VARIABLE(NAME='DAMASK_NUM_THREADS',VALUE=DAMASK_NumThreadsString,STATUS=gotDAMASK_NUM_THREADS)   ! get environment variable DAMASK_NUM_THREADS...
!$ if(gotDAMASK_NUM_THREADS /= 0_pInt) call IO_warning(47_pInt,ext_msg=DAMASK_NumThreadsString)
!$ read(DAMASK_NumThreadsString,'(i6)') DAMASK_NumThreadsInt                                        ! ...convert it to integer...
!$ if (DAMASK_NumThreadsInt < 1) DAMASK_NumThreadsInt = 1                                           ! ...ensure that its at least one...
!$ call omp_set_num_threads(DAMASK_NumThreadsInt)                                                   ! ...and use it as number of threads for parallel execution

  ! try to open the config file
  if(IO_open_file(fileunit,numerics_configFile)) then 
  
    !$OMP CRITICAL (write2out)
      write(6,*) '   ... using values from config file'
      write(6,*)
    !$OMP END CRITICAL (write2out)
    
    line = ''
    ! read variables from config file and overwrite parameters
    do
      read(fileunit,'(a1024)',END=100) line
      if (IO_isBlank(line)) cycle                           ! skip empty lines
      positions = IO_stringPos(line,maxNchunks)
      tag = IO_lc(IO_stringValue(line,positions,1_pInt))         ! extract key
      select case(tag)
        case ('relevantstrain')
              relevantStrain = IO_floatValue(line,positions,2_pInt)
        case ('defgradtolerance')
              defgradTolerance = IO_floatValue(line,positions,2_pInt)
        case ('ijacostiffness')
              iJacoStiffness = IO_intValue(line,positions,2_pInt)
        case ('ijacolpresiduum')
              iJacoLpresiduum = IO_intValue(line,positions,2_pInt)
        case ('pert_fg')
              pert_Fg = IO_floatValue(line,positions,2_pInt)
        case ('pert_method')
              pert_method = IO_intValue(line,positions,2_pInt)
        case ('nhomog')
              nHomog = IO_intValue(line,positions,2_pInt)
        case ('nmpstate')
              nMPstate = IO_intValue(line,positions,2_pInt)
        case ('ncryst')
              nCryst = IO_intValue(line,positions,2_pInt)
        case ('nstate')
              nState = IO_intValue(line,positions,2_pInt)
        case ('nstress')
              nStress = IO_intValue(line,positions,2_pInt)
        case ('substepmincryst')
              subStepMinCryst = IO_floatValue(line,positions,2_pInt)
        case ('substepsizecryst')
              subStepSizeCryst = IO_floatValue(line,positions,2_pInt)
        case ('stepincreasecryst')
              stepIncreaseCryst = IO_floatValue(line,positions,2_pInt)
        case ('substepminhomog')
              subStepMinHomog = IO_floatValue(line,positions,2_pInt)
        case ('substepsizehomog')
              subStepSizeHomog = IO_floatValue(line,positions,2_pInt)
        case ('stepincreasehomog')
              stepIncreaseHomog = IO_floatValue(line,positions,2_pInt)
        case ('rtol_crystallitestate')
              rTol_crystalliteState = IO_floatValue(line,positions,2_pInt)
        case ('rtol_crystallitetemperature')
              rTol_crystalliteTemperature = IO_floatValue(line,positions,2_pInt)
        case ('rtol_crystallitestress')
              rTol_crystalliteStress = IO_floatValue(line,positions,2_pInt)
        case ('atol_crystallitestress')
              aTol_crystalliteStress = IO_floatValue(line,positions,2_pInt)
        case ('integrator')
              numerics_integrator(1) = IO_intValue(line,positions,2_pInt)
        case ('integratorstiffness')
              numerics_integrator(2) = IO_intValue(line,positions,2_pInt)

!* RGC parameters: 
        case ('atol_rgc')
              absTol_RGC = IO_floatValue(line,positions,2_pInt)
        case ('rtol_rgc')
              relTol_RGC = IO_floatValue(line,positions,2_pInt)
        case ('amax_rgc')
              absMax_RGC = IO_floatValue(line,positions,2_pInt)
        case ('rmax_rgc')
              relMax_RGC = IO_floatValue(line,positions,2_pInt)
        case ('perturbpenalty_rgc')
              pPert_RGC = IO_floatValue(line,positions,2_pInt)
        case ('relevantmismatch_rgc')
              xSmoo_RGC = IO_floatValue(line,positions,2_pInt)
        case ('viscositypower_rgc')
              viscPower_RGC = IO_floatValue(line,positions,2_pInt)
        case ('viscositymodulus_rgc')
              viscModus_RGC = IO_floatValue(line,positions,2_pInt)
        case ('refrelaxationrate_rgc')
              refRelaxRate_RGC = IO_floatValue(line,positions,2_pInt)
        case ('maxrelaxation_rgc')
              maxdRelax_RGC = IO_floatValue(line,positions,2_pInt)
        case ('maxvoldiscrepancy_rgc')
              maxVolDiscr_RGC = IO_floatValue(line,positions,2_pInt)
        case ('voldiscrepancymod_rgc')
              volDiscrMod_RGC = IO_floatValue(line,positions,2_pInt)
        case ('discrepancypower_rgc')
              volDiscrPow_RGC = IO_floatValue(line,positions,2_pInt)

!* spectral parameters
        case ('err_div_tol')
              err_div_tol = IO_floatValue(line,positions,2_pInt)
        case ('err_stress_tolrel')
              err_stress_tolrel = IO_floatValue(line,positions,2_pInt)
        case ('itmax')
              itmax = IO_intValue(line,positions,2_pInt)
        case ('memory_efficient')
              memory_efficient = IO_intValue(line,positions,2_pInt)  > 0_pInt
        case ('fftw_timelimit')
              fftw_timelimit = IO_floatValue(line,positions,2_pInt)
        case ('fftw_planner_string')
              fftw_planner_string = IO_stringValue(line,positions,2_pInt)
        case ('rotation_tol')
              rotation_tol = IO_floatValue(line,positions,2_pInt)
        case ('divergence_correction')
              divergence_correction = IO_intValue(line,positions,2_pInt)  > 0_pInt
        case ('update_gamma')
              update_gamma = IO_intValue(line,positions,2_pInt)  > 0_pInt
        case ('simplified_algorithm')
              simplified_algorithm = IO_intValue(line,positions,2_pInt)  > 0_pInt
        case ('cut_off_value')
              cut_off_value = IO_floatValue(line,positions,2_pInt)

!* Random seeding parameters
        case ('fixed_seed')
              fixedSeed = IO_intValue(line,positions,2_pInt)
      endselect
    enddo
    100 close(fileunit)
  
  ! no config file, so we use standard values
  else 
  
    !$OMP CRITICAL (write2out)
      write(6,*) '   ... using standard values'
      write(6,*)
    !$OMP END CRITICAL (write2out)
    
  endif
  select case(IO_lc(fftw_planner_string))                        ! setting parameters for the plan creation of FFTW. Basically a translation from fftw3.f
    case('estimate','fftw_estimate')                             ! ordered from slow execution (but fast plan creation) to fast execution
       fftw_planner_flag = 64_pInt
     case('measure','fftw_measure')
       fftw_planner_flag = 0_pInt
     case('patient','fftw_patient')
       fftw_planner_flag= 32_pInt
     case('exhaustive','fftw_exhaustive')
       fftw_planner_flag = 8_pInt 
     case default
       call IO_warning(warning_ID=47_pInt,ext_msg=trim(IO_lc(fftw_planner_string)))
       fftw_planner_flag = 32_pInt
  end select

  ! writing parameters to output file
  !$OMP CRITICAL (write2out)
    write(6,'(a24,1x,e8.1)') ' relevantStrain:         ',relevantStrain
    write(6,'(a24,1x,e8.1)') ' defgradTolerance:       ',defgradTolerance
    write(6,'(a24,1x,i8)')   ' iJacoStiffness:         ',iJacoStiffness
    write(6,'(a24,1x,i8)')   ' iJacoLpresiduum:        ',iJacoLpresiduum
    write(6,'(a24,1x,e8.1)') ' pert_Fg:                ',pert_Fg
    write(6,'(a24,1x,i8)')   ' pert_method:            ',pert_method
    write(6,'(a24,1x,i8)')   ' nCryst:                 ',nCryst
    write(6,'(a24,1x,e8.1)') ' subStepMinCryst:        ',subStepMinCryst
    write(6,'(a24,1x,e8.1)') ' subStepSizeCryst:       ',subStepSizeCryst
    write(6,'(a24,1x,e8.1)') ' stepIncreaseCryst:      ',stepIncreaseCryst
    write(6,'(a24,1x,i8)')   ' nState:                 ',nState
    write(6,'(a24,1x,i8)')   ' nStress:                ',nStress
    write(6,'(a24,1x,e8.1)') ' rTol_crystalliteState:  ',rTol_crystalliteState
    write(6,'(a24,1x,e8.1)') ' rTol_crystalliteTemp:   ',rTol_crystalliteTemperature
    write(6,'(a24,1x,e8.1)') ' rTol_crystalliteStress: ',rTol_crystalliteStress
    write(6,'(a24,1x,e8.1)') ' aTol_crystalliteStress: ',aTol_crystalliteStress
    write(6,'(a24,2(1x,i8),/)')' integrator:             ',numerics_integrator
  
    write(6,'(a24,1x,i8)')   ' nHomog:                 ',nHomog
    write(6,'(a24,1x,e8.1)') ' subStepMinHomog:        ',subStepMinHomog
    write(6,'(a24,1x,e8.1)') ' subStepSizeHomog:       ',subStepSizeHomog
    write(6,'(a24,1x,e8.1)') ' stepIncreaseHomog:      ',stepIncreaseHomog
    write(6,'(a24,1x,i8,/)') ' nMPstate:               ',nMPstate

!* RGC parameters
    write(6,'(a24,1x,e8.1)') ' aTol_RGC:               ',absTol_RGC
    write(6,'(a24,1x,e8.1)') ' rTol_RGC:               ',relTol_RGC
    write(6,'(a24,1x,e8.1)') ' aMax_RGC:               ',absMax_RGC
    write(6,'(a24,1x,e8.1)') ' rMax_RGC:               ',relMax_RGC
    write(6,'(a24,1x,e8.1)') ' perturbPenalty_RGC:     ',pPert_RGC
    write(6,'(a24,1x,e8.1)') ' relevantMismatch_RGC:   ',xSmoo_RGC
    write(6,'(a24,1x,e8.1)') ' viscosityrate_RGC:      ',viscPower_RGC
    write(6,'(a24,1x,e8.1)') ' viscositymodulus_RGC:   ',viscModus_RGC
    write(6,'(a24,1x,e8.1)') ' maxrelaxation_RGC:      ',maxdRelax_RGC
    write(6,'(a24,1x,e8.1)') ' maxVolDiscrepancy_RGC:  ',maxVolDiscr_RGC
    write(6,'(a24,1x,e8.1)') ' volDiscrepancyMod_RGC:  ',volDiscrMod_RGC
    write(6,'(a24,1x,e8.1,/)') ' discrepancyPower_RGC:   ',volDiscrPow_RGC

!* spectral parameters
    write(6,'(a24,1x,e8.1)')   ' err_div_tol:            ',err_div_tol
    write(6,'(a24,1x,e8.1)')   ' err_stress_tolrel:      ',err_stress_tolrel
    write(6,'(a24,1x,i8)')     ' itmax:                  ',itmax
    write(6,'(a24,1x,L8)')     ' memory_efficient:       ',memory_efficient
    if(fftw_timelimit<0.0_pReal) then
      write(6,'(a24,1x,L8)')   ' fftw_timelimit:         ',.false.
    else    
      write(6,'(a24,1x,e8.1)') ' fftw_timelimit:         ',fftw_timelimit
    endif
    write(6,'(a24,1x,a)')      ' fftw_planner_string:    ',trim(fftw_planner_string)
    write(6,'(a24,1x,i8)')     ' fftw_planner_flag:      ',fftw_planner_flag
    write(6,'(a24,1x,e8.1)')   ' rotation_tol:           ',rotation_tol
    write(6,'(a24,1x,L8,/)')   ' divergence_correction:  ',divergence_correction
    write(6,'(a24,1x,L8,/)')   ' update_gamma:           ',update_gamma
    write(6,'(a24,1x,L8,/)')   ' simplified_algorithm:   ',simplified_algorithm
    write(6,'(a24,1x,e8.1)')   ' cut_off_value:          ',cut_off_value
!* Random seeding parameters
    write(6,'(a24,1x,i16,/)')   ' fixed_seed:             ',fixedSeed
  !$OMP END CRITICAL (write2out)

!* openMP parameter
!$  write(6,'(a24,1x,i8,/)')   ' number of threads:      ',DAMASK_NumThreadsInt
  
  ! sanity check  
  if (relevantStrain <= 0.0_pReal)          call IO_error(260_pInt)
  if (defgradTolerance <= 0.0_pReal)        call IO_error(294_pInt)
  if (iJacoStiffness < 1_pInt)              call IO_error(261_pInt)
  if (iJacoLpresiduum < 1_pInt)             call IO_error(262_pInt)
  if (pert_Fg <= 0.0_pReal)                 call IO_error(263_pInt)
  if (pert_method <= 0_pInt .or. pert_method >= 4_pInt) &
                                            call IO_error(299_pInt)
  if (nHomog < 1_pInt)                      call IO_error(264_pInt)
  if (nMPstate < 1_pInt)                    call IO_error(279_pInt)  !! missing in IO !!
  if (nCryst < 1_pInt)                      call IO_error(265_pInt)
  if (nState < 1_pInt)                      call IO_error(266_pInt)
  if (nStress < 1_pInt)                     call IO_error(267_pInt)
  if (subStepMinCryst <= 0.0_pReal)         call IO_error(268_pInt)
  if (subStepSizeCryst <= 0.0_pReal)        call IO_error(268_pInt)
  if (stepIncreaseCryst <= 0.0_pReal)       call IO_error(268_pInt)
  if (subStepMinHomog <= 0.0_pReal)         call IO_error(268_pInt)
  if (subStepSizeHomog <= 0.0_pReal)        call IO_error(268_pInt)
  if (stepIncreaseHomog <= 0.0_pReal)       call IO_error(268_pInt)
  if (rTol_crystalliteState <= 0.0_pReal)   call IO_error(269_pInt)
  if (rTol_crystalliteTemperature <= 0.0_pReal) call IO_error(276_pInt) !! oops !!
  if (rTol_crystalliteStress <= 0.0_pReal)  call IO_error(270_pInt)
  if (aTol_crystalliteStress <= 0.0_pReal)  call IO_error(271_pInt)
  if (any(numerics_integrator <= 0_pInt) .or. any(numerics_integrator >= 6_pInt)) &
                                            call IO_error(298_pInt)


  if (absTol_RGC <= 0.0_pReal)              call IO_error(272_pInt)
  if (relTol_RGC <= 0.0_pReal)              call IO_error(273_pInt)
  if (absMax_RGC <= 0.0_pReal)              call IO_error(274_pInt)
  if (relMax_RGC <= 0.0_pReal)              call IO_error(275_pInt)
  if (pPert_RGC <= 0.0_pReal)               call IO_error(276_pInt)   !! oops !!
  if (xSmoo_RGC <= 0.0_pReal)               call IO_error(277_pInt)
  if (viscPower_RGC < 0.0_pReal)            call IO_error(278_pInt)
  if (viscModus_RGC < 0.0_pReal)            call IO_error(278_pInt)
  if (refRelaxRate_RGC <= 0.0_pReal)        call IO_error(278_pInt)
  if (maxdRelax_RGC <= 0.0_pReal)           call IO_error(288_pInt)
  if (maxVolDiscr_RGC <= 0.0_pReal)         call IO_error(289_pInt)
  if (volDiscrMod_RGC < 0.0_pReal)          call IO_error(289_pInt)
  if (volDiscrPow_RGC <= 0.0_pReal)         call IO_error(289_pInt)

!* spectral parameters
  if (err_div_tol <= 0.0_pReal)             call IO_error(49_pInt)
  if (err_stress_tolrel <= 0.0_pReal)       call IO_error(49_pInt)
  if (itmax <= 1.0_pInt)                    call IO_error(49_pInt)
  
  if (fixedSeed <= 0_pInt) then
    !$OMP CRITICAL (write2out)
      write(6,'(a)') 'Random is random!'
    !$OMP END CRITICAL (write2out)
  endif
endsubroutine

END MODULE numerics
