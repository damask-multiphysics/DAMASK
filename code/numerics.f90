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
module numerics
!##############################################################

use prec, only: pInt, pReal

implicit none
character(len=64), parameter, private ::&
  numerics_configFile        = 'numerics.config'                                                    !< name of configuration file

integer(pInt) ::                iJacoStiffness             =  1_pInt, &                             !< frequency of stiffness update
                                iJacoLpresiduum            =  1_pInt, &                             !< frequency of Jacobian update of residuum in Lp
                                nHomog                     = 20_pInt, &                             !< homogenization loop limit (only for debugging info, loop limit is determined by "subStepMinHomog")
                                nMPstate                   = 10_pInt, &                             !< materialpoint state loop limit
                                nCryst                     = 20_pInt, &                             !< crystallite loop limit (only for debugging info, loop limit is determined by "subStepMinCryst")
                                nState                     = 10_pInt, &                             !< state loop limit
                                nStress                    = 40_pInt, &                             !< stress loop limit
                                pert_method                =  1_pInt, &                             !< method used in perturbation technique for tangent
                                numerics_integrationMode   =  0_pInt                                !< integrationMode 1 = central solution ; integrationMode 2 = perturbation, Default 0: undefined, is not read from file
integer(pInt), dimension(2) ::  numerics_integrator        =  1_pInt                                !< method used for state integration (central & perturbed state), Default 1: fix-point iteration for both states
real(pReal) ::                  relevantStrain             =  1.0e-7_pReal, &                       !< strain increment considered significant (used by crystallite to determine whether strain inc is considered significant)
                                defgradTolerance           =  1.0e-7_pReal, &                       !< deviation of deformation gradient that is still allowed (used by CPFEM to determine outdated ffn1)
                                pert_Fg                    =  1.0e-7_pReal, &                       !< strain perturbation for FEM Jacobi
                                subStepMinCryst            =  1.0e-3_pReal, &                       !< minimum (relative) size of sub-step allowed during cutback in crystallite
                                subStepMinHomog            =  1.0e-3_pReal, &                       !< minimum (relative) size of sub-step allowed during cutback in homogenization
                                subStepSizeCryst           =  0.25_pReal, &                         !< size of first substep when cutback in crystallite
                                subStepSizeHomog           =  0.25_pReal, &                         !< size of first substep when cutback in homogenization
                                stepIncreaseCryst          =  1.5_pReal, &                          !< increase of next substep size when previous substep converged in crystallite
                                stepIncreaseHomog          =  1.5_pReal, &                          !< increase of next substep size when previous substep converged in homogenization
                                rTol_crystalliteState      =  1.0e-6_pReal, &                       !< relative tolerance in crystallite state loop 
                                rTol_crystalliteTemperature=  1.0e-6_pReal, &                       !< relative tolerance in crystallite temperature loop 
                                rTol_crystalliteStress     =  1.0e-6_pReal, &                       !< relative tolerance in crystallite stress loop
                                aTol_crystalliteStress     =  1.0e-8_pReal, &                       !< absolute tolerance in crystallite stress loop, Default 1.0e-8: residuum is in Lp and hence strain is on this order
                                
                                absTol_RGC                 =  1.0e+4_pReal, &                       !< absolute tolerance of RGC residuum
                                relTol_RGC                 =  1.0e-3_pReal, &                       !< relative tolerance of RGC residuum
                                absMax_RGC                 =  1.0e+10_pReal, &                      !< absolute maximum of RGC residuum
                                relMax_RGC                 =  1.0e+2_pReal, &                       !< relative maximum of RGC residuum
                                pPert_RGC                  =  1.0e-7_pReal, &                       !< perturbation for computing RGC penalty tangent
                                xSmoo_RGC                  =  1.0e-5_pReal, &                       !< RGC penalty smoothing parameter (hyperbolic tangent)
                                viscPower_RGC              =  1.0e+0_pReal, &                       !< power (sensitivity rate) of numerical viscosity in RGC scheme, Default 1.0e0: Newton viscosity (linear model)
                                viscModus_RGC              =  0.0e+0_pReal, &                       !< stress modulus of RGC numerical viscosity, Default 0.0e0: No viscosity is applied
                                refRelaxRate_RGC           =  1.0e-3_pReal, &                       !< reference relaxation rate in RGC viscosity
                                maxdRelax_RGC              =  1.0e+0_pReal, &                       !< threshold of maximum relaxation vector increment (if exceed this then cutback)
                                maxVolDiscr_RGC            =  1.0e-5_pReal, &                       !< threshold of maximum volume discrepancy allowed
                                volDiscrMod_RGC            =  1.0e+12_pReal, &                      !< stiffness of RGC volume discrepancy (zero = without volume discrepancy constraint)
                                volDiscrPow_RGC            =  5.0_pReal                             !< powerlaw penalty for volume discrepancy
logical ::                      analyticJaco               = .false.                                !< use analytic Jacobian or perturbation, Default .false.: calculate Jacobian using perturbations
!* Random seeding parameters
integer(pInt) ::                fixedSeed                  = 0_pInt                                 !< fixed seeding for pseudo-random number generator, Default 0: use random seed
!* OpenMP variable
integer(pInt) ::                DAMASK_NumThreadsInt       = 0_pInt                                 !< value stored in environment variable DAMASK_NUM_THREADS, set to zero if no OpenMP directive


!* spectral parameters:
#ifdef Spectral
real(pReal) ::                  err_div_tol                =  0.1_pReal, &                          !< Div(P)/avg(P)*meter
                                err_stress_tolrel          =  0.01_pReal, &                         !< relative tolerance for fullfillment of stress BC, Default: 0.01 allowing deviation of 1% of maximum stress 
                                err_stress_tolabs          =  huge(1.0_pReal),  &                   !< absolute tolerance for fullfillment of stress BC, Default: 0.01 allowing deviation of 1% of maximum stress 
                                err_f_tol                  =  1e-6_pReal,  &
                                err_p_tol                  =  1e-5_pReal,  &
                                fftw_timelimit             = -1.0_pReal, &                          !< sets the timelimit of plan creation for FFTW, see manual on www.fftw.org, Default -1.0: disable timelimit
                                rotation_tol               =  1.0e-12_pReal                         !< tolerance of rotation specified in loadcase, Default 1.0e-12: first guess
character(len=64) ::            fftw_plan_mode             = 'FFTW_PATIENT', &                      !< reads the planing-rigor flag, see manual on www.fftw.org, Default FFTW_PATIENT: use patient planner flag
                                myspectralsolver           = 'basic'  , &                           !< spectral solution method 
                                myfilter                   = 'none'                                 !< spectral filtering method
character(len=1024) ::          petsc_options              = '-snes_type ngmres -snes_ngmres_anderson -snes_view'
integer(pInt) ::                fftw_planner_flag          =  32_pInt, &                            !< conversion of fftw_plan_mode to integer, basically what is usually done in the include file of fftw
                                itmax                      =  20_pInt, &                            !< maximum number of iterations
                                itmin                      =  2_pInt, &                             !< minimum number of iterations
                                maxCutBack                 =  3_pInt                                !< max number of cut backs
logical ::                      memory_efficient           = .true., &                              !< for fast execution (pre calculation of gamma_hat), Default .true.: do not precalculate
                                divergence_correction      = .false., &                             !< correct divergence calculation in fourier space, Default .false.: no correction
                                update_gamma               = .false.                                !< update gamma operator with current stiffness, Default .false.: use initial stiffness 
#endif



CONTAINS
 
!*******************************************
!    initialization subroutine
!*******************************************
subroutine numerics_init
  
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use IO, only:                               IO_error, &
                                             IO_open_file_stat, &
                                             IO_isBlank, &
                                             IO_stringPos, &
                                             IO_stringValue, &
                                             IO_lc, &
                                             IO_floatValue, &
                                             IO_intValue, &
                                             IO_warning
#ifndef Marc                                                                                      ! Use the standard conforming module file for omp if using the spectral solver
!$ use OMP_LIB, only: omp_set_num_threads
#endif
 implicit none
#ifdef Marc                                                                                       ! use the non F90 standard include file because some versions of Marc and Abaqus crash when using the module
!$ include "omp_lib.h"
#endif
 integer(pInt), parameter ::                 fileunit = 300_pInt ,&
                                             maxNchunks = 2_pInt
!$ integer ::                                gotDAMASK_NUM_THREADS = 1
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 character(len=64) ::                        tag
 character(len=1024) ::                      line
!$ character(len=6) DAMASK_NumThreadsString                                                         ! environment variable DAMASK_NUM_THREADS

 write(6,*)
 write(6,*) '<<<+-  numerics init  -+>>>'
 write(6,*) '$Id$'
#include "compilation_info.f90"


!$ call GET_ENVIRONMENT_VARIABLE(NAME='DAMASK_NUM_THREADS',VALUE=DAMASK_NumThreadsString,STATUS=gotDAMASK_NUM_THREADS)   ! get environment variable DAMASK_NUM_THREADS...
!$ if(gotDAMASK_NUM_THREADS /= 0) call IO_warning(47_pInt,ext_msg=DAMASK_NumThreadsString)
!$ read(DAMASK_NumThreadsString,'(i6)') DAMASK_NumThreadsInt                                        ! ...convert it to integer...
!$ if (DAMASK_NumThreadsInt < 1_pInt) DAMASK_NumThreadsInt = 1_pInt                                 ! ...ensure that its at least one...
!$ call omp_set_num_threads(DAMASK_NumThreadsInt)                                                   ! ...and use it as number of threads for parallel execution

  ! try to open the config file
 if(IO_open_file_stat(fileunit,numerics_configFile)) then 

   write(6,*) '   ... using values from config file'
   write(6,*)
    
   !* read variables from config file and overwrite parameters
   
   line = ''
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
       case ('analyticjaco')
             analyticJaco = IO_intValue(line,positions,2_pInt) > 0_pInt

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
       !* Random seeding parameters
       case ('fixed_seed')
             fixedSeed = IO_intValue(line,positions,2_pInt)
       !* spectral parameters
#ifdef Spectral
       case ('err_div_tol')
             err_div_tol = IO_floatValue(line,positions,2_pInt)
       case ('err_stress_tolrel')
             err_stress_tolrel = IO_floatValue(line,positions,2_pInt)
       case ('err_stress_tolabs')
             err_stress_tolabs = IO_floatValue(line,positions,2_pInt)
       case ('itmax')
             itmax = IO_intValue(line,positions,2_pInt)
       case ('itmin')
             itmin = IO_intValue(line,positions,2_pInt)
       case ('maxcutback')
             maxCutBack = IO_intValue(line,positions,2_pInt)
       case ('memory_efficient')
             memory_efficient = IO_intValue(line,positions,2_pInt)  > 0_pInt
       case ('fftw_timelimit')
             fftw_timelimit = IO_floatValue(line,positions,2_pInt)
       case ('fftw_plan_mode')
             fftw_plan_mode = IO_stringValue(line,positions,2_pInt)
       case ('myfilter')
             myfilter = IO_stringValue(line,positions,2_pInt)
       case ('rotation_tol')
             rotation_tol = IO_floatValue(line,positions,2_pInt)
       case ('divergence_correction')
             divergence_correction = IO_intValue(line,positions,2_pInt)  > 0_pInt
       case ('update_gamma')
             update_gamma = IO_intValue(line,positions,2_pInt)  > 0_pInt
#ifdef PETSc
       case ('petsc_options')
             petsc_options = trim(line(positions(4):))
       case ('myspectralsolver')
             myspectralsolver = IO_stringValue(line,positions,2_pInt)
       case ('err_f_tol')
             err_f_tol = IO_floatValue(line,positions,2_pInt)
       case ('err_p_tol')
             err_p_tol = IO_floatValue(line,positions,2_pInt)
#endif 
#ifndef PETSc
       case ('myspectralsolver', 'petsc_options','err_f_tol', 'err_p_tol')
             call IO_warning(41_pInt,ext_msg=tag)
#endif
#endif
#ifndef Spectral
      case ('err_div_tol','err_stress_tolrel','err_stress_tolabs',&
            'itmax', 'itmin','memory_efficient','fftw_timelimit','fftw_plan_mode','myspectralsolver', &
            'rotation_tol','divergence_correction','update_gamma','petsc_options','myfilter', &
            'err_f_tol', 'err_p_tol', 'maxcutback')
             call IO_warning(40_pInt,ext_msg=tag)
#endif
       case default 
             call IO_error(300_pInt,ext_msg=tag)
     endselect
   enddo
   100 close(fileunit)
  
  ! no config file, so we use standard values
 else 
   write(6,*) '   ... using standard values'
   write(6,*)
 endif
#ifdef Spectral
 select case(IO_lc(fftw_plan_mode))                                                                ! setting parameters for the plan creation of FFTW. Basically a translation from fftw3.f
    case('estimate','fftw_estimate')                                                               ! ordered from slow execution (but fast plan creation) to fast execution
      fftw_planner_flag = 64_pInt
    case('measure','fftw_measure')
      fftw_planner_flag = 0_pInt
    case('patient','fftw_patient')
      fftw_planner_flag= 32_pInt
    case('exhaustive','fftw_exhaustive')
      fftw_planner_flag = 8_pInt 
    case default
      call IO_warning(warning_ID=47_pInt,ext_msg=trim(IO_lc(fftw_plan_mode)))
      fftw_planner_flag = 32_pInt
 end select
#endif
 
 !* writing parameters to output file
   
   write(6,'(a24,1x,es8.1)')  ' relevantStrain:         ',relevantStrain
   write(6,'(a24,1x,es8.1)')  ' defgradTolerance:       ',defgradTolerance
   write(6,'(a24,1x,i8)')     ' iJacoStiffness:         ',iJacoStiffness
   write(6,'(a24,1x,i8)')     ' iJacoLpresiduum:        ',iJacoLpresiduum
   write(6,'(a24,1x,es8.1)')  ' pert_Fg:                ',pert_Fg
   write(6,'(a24,1x,i8)')     ' pert_method:            ',pert_method
   write(6,'(a24,1x,i8)')     ' nCryst:                 ',nCryst
   write(6,'(a24,1x,es8.1)')  ' subStepMinCryst:        ',subStepMinCryst
   write(6,'(a24,1x,es8.1)')  ' subStepSizeCryst:       ',subStepSizeCryst
   write(6,'(a24,1x,es8.1)')  ' stepIncreaseCryst:      ',stepIncreaseCryst
   write(6,'(a24,1x,i8)')     ' nState:                 ',nState
   write(6,'(a24,1x,i8)')     ' nStress:                ',nStress
   write(6,'(a24,1x,es8.1)')  ' rTol_crystalliteState:  ',rTol_crystalliteState
   write(6,'(a24,1x,es8.1)')  ' rTol_crystalliteTemp:   ',rTol_crystalliteTemperature
   write(6,'(a24,1x,es8.1)')  ' rTol_crystalliteStress: ',rTol_crystalliteStress
   write(6,'(a24,1x,es8.1)')  ' aTol_crystalliteStress: ',aTol_crystalliteStress
   write(6,'(a24,2(1x,i8))')  ' integrator:             ',numerics_integrator
   write(6,'(a24,1x,L8,/)')   ' analytic Jacobian:      ',analyticJaco
 
   write(6,'(a24,1x,i8)')     ' nHomog:                 ',nHomog
   write(6,'(a24,1x,es8.1)')  ' subStepMinHomog:        ',subStepMinHomog
   write(6,'(a24,1x,es8.1)')  ' subStepSizeHomog:       ',subStepSizeHomog
   write(6,'(a24,1x,es8.1)')  ' stepIncreaseHomog:      ',stepIncreaseHomog
   write(6,'(a24,1x,i8,/)')   ' nMPstate:               ',nMPstate

   !* RGC parameters

   write(6,'(a24,1x,es8.1)')   ' aTol_RGC:               ',absTol_RGC
   write(6,'(a24,1x,es8.1)')   ' rTol_RGC:               ',relTol_RGC
   write(6,'(a24,1x,es8.1)')   ' aMax_RGC:               ',absMax_RGC
   write(6,'(a24,1x,es8.1)')   ' rMax_RGC:               ',relMax_RGC
   write(6,'(a24,1x,es8.1)')   ' perturbPenalty_RGC:     ',pPert_RGC
   write(6,'(a24,1x,es8.1)')   ' relevantMismatch_RGC:   ',xSmoo_RGC
   write(6,'(a24,1x,es8.1)')   ' viscosityrate_RGC:      ',viscPower_RGC
   write(6,'(a24,1x,es8.1)')   ' viscositymodulus_RGC:   ',viscModus_RGC
   write(6,'(a24,1x,es8.1)')   ' maxrelaxation_RGC:      ',maxdRelax_RGC
   write(6,'(a24,1x,es8.1)')   ' maxVolDiscrepancy_RGC:  ',maxVolDiscr_RGC
   write(6,'(a24,1x,es8.1)')   ' volDiscrepancyMod_RGC:  ',volDiscrMod_RGC
   write(6,'(a24,1x,es8.1,/)') ' discrepancyPower_RGC:   ',volDiscrPow_RGC
   !* Random seeding parameters
   
   write(6,'(a24,1x,i16,/)')    ' fixed_seed:             ',fixedSeed
  !* openMP parameter
  !$  write(6,'(a24,1x,i8,/)')   ' number of threads:      ',DAMASK_NumThreadsInt

  !* spectral parameters
#ifdef Spectral
   write(6,'(a24,1x,es8.1)')   ' err_div_tol:            ',err_div_tol
   write(6,'(a24,1x,es8.1)')   ' err_stress_tolrel:      ',err_stress_tolrel
   write(6,'(a24,1x,es8.1)')   ' err_stress_tolabs:      ',err_stress_tolabs

   write(6,'(a24,1x,i8)')      ' itmax:                  ',itmax
   write(6,'(a24,1x,i8)')      ' itmin:                  ',itmin
   write(6,'(a24,1x,i8)')      ' maxCutBack:             ',maxCutBack
   write(6,'(a24,1x,L8)')      ' memory_efficient:       ',memory_efficient
   if(fftw_timelimit<0.0_pReal) then
     write(6,'(a24,1x,L8)')    ' fftw_timelimit:         ',.false.
   else    
     write(6,'(a24,1x,es8.1)') ' fftw_timelimit:         ',fftw_timelimit
   endif
   write(6,'(a24,1x,a)')       ' fftw_plan_mode:         ',trim(fftw_plan_mode)

   write(6,'(a24,1x,a)')       ' myfilter:               ',trim(myfilter)
   write(6,'(a24,1x,i8)')      ' fftw_planner_flag:      ',fftw_planner_flag
   write(6,'(a24,1x,es8.1)')   ' rotation_tol:           ',rotation_tol
   write(6,'(a24,1x,L8,/)')    ' divergence_correction:  ',divergence_correction
   write(6,'(a24,1x,L8,/)')    ' update_gamma:           ',update_gamma
#ifdef PETSc
   write(6,'(a24,1x,es8.1)')   ' err_f_tol:              ',err_f_tol
   write(6,'(a24,1x,es8.1)')   ' err_p_tol:              ',err_p_tol
   write(6,'(a24,1x,a)')       ' myspectralsolver:       ',trim(myspectralsolver)
   write(6,'(a24,1x,a)')       ' PETSc_options:          ',trim(petsc_options)
#endif
#endif

 !*  sanity check

 if (relevantStrain <= 0.0_pReal)          call IO_error(301_pInt,ext_msg='relevantStrain')
 if (defgradTolerance <= 0.0_pReal)        call IO_error(301_pInt,ext_msg='defgradTolerance')
 if (iJacoStiffness < 1_pInt)              call IO_error(301_pInt,ext_msg='iJacoStiffness')
 if (iJacoLpresiduum < 1_pInt)             call IO_error(301_pInt,ext_msg='iJacoLpresiduum')
 if (pert_Fg <= 0.0_pReal)                 call IO_error(301_pInt,ext_msg='pert_Fg')
 if (pert_method <= 0_pInt .or. pert_method >= 4_pInt) &
                                           call IO_error(301_pInt,ext_msg='pert_method')
 if (nHomog < 1_pInt)                      call IO_error(301_pInt,ext_msg='nHomog')
 if (nMPstate < 1_pInt)                    call IO_error(301_pInt,ext_msg='nMPstate')
 if (nCryst < 1_pInt)                      call IO_error(301_pInt,ext_msg='nCryst')
 if (nState < 1_pInt)                      call IO_error(301_pInt,ext_msg='nState')
 if (nStress < 1_pInt)                     call IO_error(301_pInt,ext_msg='nStress')
 if (subStepMinCryst <= 0.0_pReal)         call IO_error(301_pInt,ext_msg='subStepMinCryst')
 if (subStepSizeCryst <= 0.0_pReal)        call IO_error(301_pInt,ext_msg='subStepSizeCryst')
 if (stepIncreaseCryst <= 0.0_pReal)       call IO_error(301_pInt,ext_msg='stepIncreaseCryst')
 if (subStepMinHomog <= 0.0_pReal)         call IO_error(301_pInt,ext_msg='subStepMinHomog')
 if (subStepSizeHomog <= 0.0_pReal)        call IO_error(301_pInt,ext_msg='subStepSizeHomog')
 if (stepIncreaseHomog <= 0.0_pReal)       call IO_error(301_pInt,ext_msg='stepIncreaseHomog')
 if (rTol_crystalliteState <= 0.0_pReal)   call IO_error(301_pInt,ext_msg='rTol_crystalliteState')
 if (rTol_crystalliteTemperature <= 0.0_pReal) call IO_error(301_pInt,ext_msg='rTol_crystalliteTemperature')
 if (rTol_crystalliteStress <= 0.0_pReal)  call IO_error(301_pInt,ext_msg='rTol_crystalliteStress')
 if (aTol_crystalliteStress <= 0.0_pReal)  call IO_error(301_pInt,ext_msg='aTol_crystalliteStress')
 if (any(numerics_integrator <= 0_pInt) .or. any(numerics_integrator >= 6_pInt)) &
                                           call IO_error(301_pInt,ext_msg='integrator')

 !* RGC parameters
 if (absTol_RGC <= 0.0_pReal)              call IO_error(301_pInt,ext_msg='absTol_RGC')
 if (relTol_RGC <= 0.0_pReal)              call IO_error(301_pInt,ext_msg='relTol_RGC')
 if (absMax_RGC <= 0.0_pReal)              call IO_error(301_pInt,ext_msg='absMax_RGC')
 if (relMax_RGC <= 0.0_pReal)              call IO_error(301_pInt,ext_msg='relMax_RGC')
 if (pPert_RGC <= 0.0_pReal)               call IO_error(301_pInt,ext_msg='pPert_RGC')
 if (xSmoo_RGC <= 0.0_pReal)               call IO_error(301_pInt,ext_msg='xSmoo_RGC')
 if (viscPower_RGC < 0.0_pReal)            call IO_error(301_pInt,ext_msg='viscPower_RGC')
 if (viscModus_RGC < 0.0_pReal)            call IO_error(301_pInt,ext_msg='viscModus_RGC')
 if (refRelaxRate_RGC <= 0.0_pReal)        call IO_error(301_pInt,ext_msg='refRelaxRate_RGC')
 if (maxdRelax_RGC <= 0.0_pReal)           call IO_error(301_pInt,ext_msg='maxdRelax_RGC')
 if (maxVolDiscr_RGC <= 0.0_pReal)         call IO_error(301_pInt,ext_msg='maxVolDiscr_RGC')
 if (volDiscrMod_RGC < 0.0_pReal)          call IO_error(301_pInt,ext_msg='volDiscrMod_RGC')
 if (volDiscrPow_RGC <= 0.0_pReal)         call IO_error(301_pInt,ext_msg='volDiscrPw_RGC')

 !* spectral parameters
#ifdef Spectral
 if (err_div_tol <= 0.0_pReal)             call IO_error(301_pInt,ext_msg='err_div_tol')
 if (err_stress_tolrel <= 0.0_pReal)       call IO_error(301_pInt,ext_msg='err_stress_tolrel')
 if (err_stress_tolabs <= 0.0_pReal)       call IO_error(301_pInt,ext_msg='err_stress_tolabs')
 if (itmax <= 1.0_pInt)                    call IO_error(301_pInt,ext_msg='itmax')
 if (itmin > itmax .or. itmin < 1_pInt)    call IO_error(301_pInt,ext_msg='itmin')
 if (maxCutBack <= 1.0_pInt)               call IO_error(301_pInt,ext_msg='maxCutBack')
 if (update_gamma .and. &
                   .not. memory_efficient) call IO_error(error_ID = 847_pInt)
#ifdef PETSc
 if (err_f_tol <= 0.0_pReal)               call IO_error(301_pInt,ext_msg='err_f_tol')
 if (err_p_tol <= 0.0_pReal)               call IO_error(301_pInt,ext_msg='err_p_tol')
#endif
#endif
 if (fixedSeed <= 0_pInt) then
   write(6,'(a,/)') ' Random is random!'
 endif

end subroutine numerics_init

end module numerics
