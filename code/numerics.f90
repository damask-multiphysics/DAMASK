!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Managing of parameters related to numerics
!--------------------------------------------------------------------------------------------------
module numerics
 use prec, only: &
   pInt, &
   pReal

 implicit none
 private
#ifdef PETSc
#include <finclude/petsc.h90>
#endif
 character(len=64), parameter, private :: &
   numerics_CONFIGFILE        = 'numerics.config'                                                   !< name of configuration file

 integer(pInt), protected, public :: &
   iJacoStiffness             =  1_pInt, &                                                          !< frequency of stiffness update
   iJacoLpresiduum            =  1_pInt, &                                                          !< frequency of Jacobian update of residuum in Lp
   nHomog                     = 20_pInt, &                                                          !< homogenization loop limit (only for debugging info, loop limit is determined by "subStepMinHomog")
   nMPstate                   = 10_pInt, &                                                          !< materialpoint state loop limit
   nCryst                     = 20_pInt, &                                                          !< crystallite loop limit (only for debugging info, loop limit is determined by "subStepMinCryst")
   nState                     = 10_pInt, &                                                          !< state loop limit
   nStress                    = 40_pInt, &                                                          !< stress loop limit
   pert_method                =  1_pInt, &                                                          !< method used in perturbation technique for tangent
   fixedSeed                  =  0_pInt, &                                                          !< fixed seeding for pseudo-random number generator, Default 0: use random seed
   worldrank                  =  0_pInt, &                                                          !< MPI worldrank (/=0 for MPI simulations only)
   worldsize                  =  0_pInt                                                             !< MPI worldsize (/=0 for MPI simulations only)
 integer, protected, public :: &
   DAMASK_NumThreadsInt       =  0                                                                  !< value stored in environment variable DAMASK_NUM_THREADS, set to zero if no OpenMP directive
 integer(pInt), public :: &
   numerics_integrationMode   =  0_pInt                                                             !< integrationMode 1 = central solution ; integrationMode 2 = perturbation, Default 0: undefined, is not read from file
 integer(pInt), dimension(2) , protected, public :: &
   numerics_integrator        =  1_pInt                                                             !< method used for state integration (central & perturbed state), Default 1: fix-point iteration for both states
 real(pReal), protected, public :: &
   relevantStrain             =  1.0e-7_pReal, &                                                    !< strain increment considered significant (used by crystallite to determine whether strain inc is considered significant)
   defgradTolerance           =  1.0e-7_pReal, &                                                    !< deviation of deformation gradient that is still allowed (used by CPFEM to determine outdated ffn1)
   pert_Fg                    =  1.0e-7_pReal, &                                                    !< strain perturbation for FEM Jacobi
   subStepMinCryst            =  1.0e-3_pReal, &                                                    !< minimum (relative) size of sub-step allowed during cutback in crystallite
   subStepMinHomog            =  1.0e-3_pReal, &                                                    !< minimum (relative) size of sub-step allowed during cutback in homogenization
   subStepSizeCryst           =  0.25_pReal, &                                                      !< size of first substep when cutback in crystallite
   subStepSizeHomog           =  0.25_pReal, &                                                      !< size of first substep when cutback in homogenization
   stepIncreaseCryst          =  1.5_pReal, &                                                       !< increase of next substep size when previous substep converged in crystallite
   stepIncreaseHomog          =  1.5_pReal, &                                                       !< increase of next substep size when previous substep converged in homogenization
   rTol_crystalliteState      =  1.0e-6_pReal, &                                                    !< relative tolerance in crystallite state loop 
   rTol_crystalliteStress     =  1.0e-6_pReal, &                                                    !< relative tolerance in crystallite stress loop
   aTol_crystalliteStress     =  1.0e-8_pReal, &                                                    !< absolute tolerance in crystallite stress loop, Default 1.0e-8: residuum is in Lp and hence strain is on this order
   numerics_unitlength        =  1.0_pReal, &                                                       !< determines the physical length of one computational length unit
   absTol_RGC                 =  1.0e+4_pReal, &                                                    !< absolute tolerance of RGC residuum
   relTol_RGC                 =  1.0e-3_pReal, &                                                    !< relative tolerance of RGC residuum
   absMax_RGC                 =  1.0e+10_pReal, &                                                   !< absolute maximum of RGC residuum
   relMax_RGC                 =  1.0e+2_pReal, &                                                    !< relative maximum of RGC residuum
   pPert_RGC                  =  1.0e-7_pReal, &                                                    !< perturbation for computing RGC penalty tangent
   xSmoo_RGC                  =  1.0e-5_pReal, &                                                    !< RGC penalty smoothing parameter (hyperbolic tangent)
   viscPower_RGC              =  1.0e+0_pReal, &                                                    !< power (sensitivity rate) of numerical viscosity in RGC scheme, Default 1.0e0: Newton viscosity (linear model)
   viscModus_RGC              =  0.0e+0_pReal, &                                                    !< stress modulus of RGC numerical viscosity, Default 0.0e0: No viscosity is applied
   refRelaxRate_RGC           =  1.0e-3_pReal, &                                                    !< reference relaxation rate in RGC viscosity
   maxdRelax_RGC              =  1.0e+0_pReal, &                                                    !< threshold of maximum relaxation vector increment (if exceed this then cutback)
   maxVolDiscr_RGC            =  1.0e-5_pReal, &                                                    !< threshold of maximum volume discrepancy allowed
   volDiscrMod_RGC            =  1.0e+12_pReal, &                                                   !< stiffness of RGC volume discrepancy (zero = without volume discrepancy constraint)
   volDiscrPow_RGC            =  5.0_pReal, &                                                       !< powerlaw penalty for volume discrepancy
   charLength                 =  1.0_pReal, &                                                       !< characteristic length scale for gradient problems
   residualStiffness          =  1.0e-6_pReal                                                       !< non-zero residual damage   
 logical, protected, public :: &                                                   
#if defined(Spectral) || defined(FEM)
   analyticJaco               = .true.,  &                                                          !< use analytic Jacobian or perturbation, Default for Spectral solver .true.:
#else
   analyticJaco               = .false., &                                                          !< use analytic Jacobian or perturbation, Default .false.: calculate Jacobian using perturbations
#endif
   usePingPong                = .true., & 
   numerics_timeSyncing       = .false.                                                             !< flag indicating if time synchronization in crystallite is used for nonlocal plasticity

!--------------------------------------------------------------------------------------------------
! spectral parameters:
#ifdef Spectral
 real(pReal), protected, public :: &
   err_div_tolAbs             =  1.0e-10_pReal, &                                                   !< absolute tolerance for equilibrium
   err_div_tolRel             =  5.0e-4_pReal, &                                                    !< relative tolerance for equilibrium
   err_curl_tolAbs            =  1.0e-10_pReal, &                                                   !< absolute tolerance for compatibility
   err_curl_tolRel            =  5.0e-4_pReal, &                                                    !< relative tolerance for compatibility
   err_stress_tolAbs          =  1.0e3_pReal,  &                                                    !< absolute tolerance for fullfillment of stress BC
   err_stress_tolRel          =  0.01_pReal, &                                                      !< relative tolerance for fullfillment of stress BC
   fftw_timelimit             = -1.0_pReal, &                                                       !< sets the timelimit of plan creation for FFTW, see manual on www.fftw.org, Default -1.0: disable timelimit
   rotation_tol               =  1.0e-12_pReal, &                                                   !< tolerance of rotation specified in loadcase, Default 1.0e-12: first guess
   polarAlpha                 =  1.0_pReal, &                                                       !< polarization scheme parameter 0.0 < alpha < 2.0. alpha = 1.0 ==> AL scheme, alpha = 2.0 ==> accelerated scheme 
   polarBeta                  =  1.0_pReal                                                          !< polarization scheme parameter 0.0 < beta < 2.0. beta = 1.0 ==> AL scheme, beta = 2.0 ==> accelerated scheme 
 character(len=64), private :: &
   fftw_plan_mode             = 'FFTW_PATIENT'                                                      !< reads the planing-rigor flag, see manual on www.fftw.org, Default FFTW_PATIENT: use patient planner flag
 character(len=64), protected, public :: & 
   spectral_solver            = 'basic'  , &                                                        !< spectral solution method 
   spectral_filter            = 'none'                                                              !< spectral filtering method
 character(len=1024), protected, public :: &
   petsc_options              = '-snes_type ngmres &
                                &-snes_ngmres_anderson '
 integer(pInt), protected, public :: &
   fftw_planner_flag          =  32_pInt, &                                                         !< conversion of fftw_plan_mode to integer, basically what is usually done in the include file of fftw
   itmax                      =  250_pInt, &                                                        !< maximum number of iterations
   itmin                      =  2_pInt, &                                                          !< minimum number of iterations
   maxCutBack                 =  3_pInt, &                                                          !< max number of cut backs
   continueCalculation        =  0_pInt, &                                                          !< 0: exit if BVP solver does not converge, 1: continue calculation if BVP solver does not converge
   divergence_correction      =  2_pInt                                                             !< correct divergence calculation in fourier space 0: no correction, 1: size scaled to 1, 2: size scaled to Npoints
 logical, protected, public :: &
   memory_efficient           = .true., &                                                           !< for fast execution (pre calculation of gamma_hat), Default .true.: do not precalculate
   update_gamma               = .false.                                                             !< update gamma operator with current stiffness, Default .false.: use initial stiffness 
#endif

!--------------------------------------------------------------------------------------------------
! FEM parameters:
#ifdef FEM
 real(pReal), protected, public :: &
   err_struct_tolAbs          =  1.0e-10_pReal, &                                                   !< absolute tolerance for equilibrium
   err_struct_tolRel          =  1.0e-4_pReal, &                                                    !< relative tolerance for equilibrium
   err_thermal_tol            =  1.0e-1_pReal, &
   err_damage_tol             =  1.0e-2_pReal, &
   err_vacancydiffusion_tol   =  1.0e-8_pReal
 character(len=4096), protected, public :: &
   petsc_optionsFEM        = '-mech_snes_type newtonls &
                             &-mech_snes_linesearch_type cp &
                             &-mech_snes_ksp_ew &
                             &-mech_snes_ksp_ew_rtol0 0.01 &
                             &-mech_snes_ksp_ew_rtolmax 0.01 &
                             &-mech_ksp_type fgmres &
                             &-mech_ksp_max_it 25 &
                             &-mech_pc_type ml &
                             &-mech_pc_ml_maxNlevels 2 &
                             &-mech_mg_coarse_ksp_type preonly &
                             &-mech_mg_coarse_pc_type lu &
                             &-mech_mg_coarse_pc_factor_mat_solver_package superlu_dist &
                             &-mech_mg_levels_ksp_type chebyshev &
                             &-mech_mg_levels_ksp_chebyshev_estimate_eigenvalues 0,0.1,0,1.1 &
                             &-mech_mg_levels_pc_type sor &
                             &-mech_pc_ml_nullspace user &
                             &-damage_snes_type newtonls &
                             &-damage_snes_linesearch_type cp &
                             &-damage_ksp_type fgmres &
                             &-damage_snes_atol 1e-8 &
                             &-damage_pc_type ml &
                             &-damage_mg_levels_ksp_type chebyshev &
                             &-damage_mg_levels_ksp_chebyshev_estimate_eigenvalues 0,0.1,0,1.1 &
                             &-damage_mg_levels_pc_type sor &
                             &-thermal_snes_type newtonls &
                             &-thermal_snes_linesearch_type cp &
                             &-thermal_ksp_type fgmres &
                             &-thermal_snes_atol 1e-1 &
                             &-thermal_pc_type ml &
                             &-thermal_mg_levels_ksp_type chebyshev &
                             &-thermal_mg_levels_ksp_chebyshev_estimate_eigenvalues 0,0.1,0,1.1 &
                             &-thermal_mg_levels_pc_type sor &
                             &-vacancyDiffusion_snes_type newtonls &
                             &-vacancyDiffusion_snes_linesearch_type cp &
                             &-vacancyDiffusion_ksp_type fgmres &
                             &-vacancyDiffusion_snes_atol 1e-6 &
                             &-vacancyDiffusion_pc_type ml &
                             &-vacancyDiffusion_mg_levels_ksp_type chebyshev &
                             &-vacancyDiffusion_mg_levels_ksp_chebyshev_estimate_eigenvalues 0,0.1,0,1.1 &
                             &-vacancyDiffusion_mg_levels_pc_type sor '
 integer(pInt), protected, public :: &
   itmaxFEM                   =  25_pInt, &                                                         !< maximum number of iterations
   itminFEM                   =  1_pInt, &                                                          !< minimum number of iterations
   stagItMax                  =  10_pInt, &                                                         !< max number of field level staggered iterations
   maxCutBackFEM              =  3_pInt, &                                                          !< max number of cut backs
   integrationOrder           =  2_pInt, &
   structOrder                =  2_pInt, &
   thermalOrder               =  2_pInt, &
   damageOrder                =  2_pInt, &
   vacancyDiffusionOrder      =  2_pInt
 logical, protected, public :: & 
   BBarStabilisation          = .false.                                                  
#endif

 public :: numerics_init
  
contains


!--------------------------------------------------------------------------------------------------
!> @brief reads in parameters from numerics.config and sets openMP related parameters. Also does
! a sanity check
!--------------------------------------------------------------------------------------------------
subroutine numerics_init
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use IO, only: &
   IO_read, &
   IO_error, &
   IO_open_file_stat, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_lc, &
   IO_floatValue, &
   IO_intValue, &
   IO_warning, &
   IO_timeStamp, &
   IO_EOF

#if defined(Spectral) || defined(FEM)
!$ use OMP_LIB, only: omp_set_num_threads                                                           ! Use the standard conforming module file for omp if using the spectral solver
 implicit none
#else  
 implicit none
!$ include "omp_lib.h"                                                                              ! use the not F90 standard conforming include file to prevent crashes with some versions of MSC.Marc
#endif
 integer(pInt), parameter ::                 FILEUNIT = 300_pInt ,&
                                             maxNchunks = 2_pInt
!$ integer ::                                gotDAMASK_NUM_THREADS = 1
 integer :: i, ierr                                                                                 ! no pInt
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 character(len=65536) :: &
   tag ,&
   line
!$ character(len=6) DAMASK_NumThreadsString                                                         ! environment variable DAMASK_NUM_THREADS

#ifdef PETSc
 call MPI_Comm_rank(PETSC_COMM_WORLD,worldrank,ierr);CHKERRQ(ierr)
 call MPI_Comm_size(PETSC_COMM_WORLD,worldsize,ierr);CHKERRQ(ierr)
#endif
 mainProcess: if (worldrank == 0) then
   write(6,'(/,a)') ' <<<+-  numerics init  -+>>>'
   write(6,'(a)')   ' $Id$'
   write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
!$ call GET_ENVIRONMENT_VARIABLE(NAME='DAMASK_NUM_THREADS',VALUE=DAMASK_NumThreadsString,STATUS=gotDAMASK_NUM_THREADS)   ! get environment variable DAMASK_NUM_THREADS...
!$ if(gotDAMASK_NUM_THREADS /= 0) then                                                              ! could not get number of threads, set it to 1
!$   call IO_warning(35_pInt,ext_msg='BEGIN:'//DAMASK_NumThreadsString//':END')
!$   DAMASK_NumThreadsInt = 1
!$ else
!$   read(DAMASK_NumThreadsString,'(i6)') DAMASK_NumThreadsInt                                      ! read as integer
!$   if (DAMASK_NumThreadsInt < 1) DAMASK_NumThreadsInt = 1                                         ! in case of string conversion fails, set it to one
!$ endif
!$ call omp_set_num_threads(DAMASK_NumThreadsInt)                                                   ! set number of threads for parallel execution

!--------------------------------------------------------------------------------------------------
! try to open the config file
 fileExists: if(IO_open_file_stat(FILEUNIT,numerics_configFile)) then 
   mainProcess2: if (worldrank == 0) then
     write(6,'(a,/)') ' using values from config file'
     flush(6)
   endif mainProcess2
    
!--------------------------------------------------------------------------------------------------
! read variables from config file and overwrite default parameters if keyword is present
   line = ''
   do while (trim(line) /= IO_EOF)                                                                  ! read thru sections of phase part
     line = IO_read(FILEUNIT)
     do i=1,len(line)
       if(line(i:i) == '=') line(i:i) = ' '                                                         ! also allow keyword = value version
     enddo
     if (IO_isBlank(line)) cycle                                                                    ! skip empty lines
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
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
       case ('usepingpong')
         usepingpong = IO_intValue(line,positions,2_pInt) > 0_pInt
       case ('timesyncing')
         numerics_timeSyncing = IO_intValue(line,positions,2_pInt) > 0_pInt
       case ('unitlength')
         numerics_unitlength = IO_floatValue(line,positions,2_pInt)

!--------------------------------------------------------------------------------------------------
! RGC parameters
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

!--------------------------------------------------------------------------------------------------
! random seeding parameter
       case ('fixed_seed')
         fixedSeed = IO_intValue(line,positions,2_pInt)

!--------------------------------------------------------------------------------------------------
! gradient parameter
       case ('charlength')
         charLength = IO_floatValue(line,positions,2_pInt)
       case ('residualstiffness')
         residualStiffness = IO_floatValue(line,positions,2_pInt)

!--------------------------------------------------------------------------------------------------
! spectral parameters
#ifdef Spectral
       case ('err_div_tolabs')
         err_div_tolAbs = IO_floatValue(line,positions,2_pInt)
       case ('err_div_tolrel')
         err_div_tolRel = IO_floatValue(line,positions,2_pInt)
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
       case ('continuecalculation')
         continueCalculation = IO_intValue(line,positions,2_pInt)
       case ('memory_efficient')
         memory_efficient = IO_intValue(line,positions,2_pInt)  > 0_pInt
       case ('fftw_timelimit')
         fftw_timelimit = IO_floatValue(line,positions,2_pInt)
       case ('fftw_plan_mode')
         fftw_plan_mode = IO_lc(IO_stringValue(line,positions,2_pInt))
       case ('spectralfilter','myfilter')
         spectral_filter = IO_lc(IO_stringValue(line,positions,2_pInt))
       case ('divergence_correction')
         divergence_correction = IO_intValue(line,positions,2_pInt)
       case ('update_gamma')
         update_gamma = IO_intValue(line,positions,2_pInt)  > 0_pInt
       case ('petsc_options')
         petsc_options = trim(line(positions(4):))
       case ('spectralsolver','myspectralsolver')
         spectral_solver = IO_lc(IO_stringValue(line,positions,2_pInt))
       case ('err_curl_tolabs')
         err_curl_tolAbs = IO_floatValue(line,positions,2_pInt)
       case ('err_curl_tolrel')
         err_curl_tolRel = IO_floatValue(line,positions,2_pInt)
       case ('polaralpha')
         polarAlpha = IO_floatValue(line,positions,2_pInt)
       case ('polarbeta')
         polarBeta = IO_floatValue(line,positions,2_pInt)
#else
      case ('err_div_tolabs','err_div_tolrel','err_stress_tolrel','err_stress_tolabs',&             ! found spectral parameter for FEM build
            'itmax', 'itmin','memory_efficient','fftw_timelimit','fftw_plan_mode', &
            'divergence_correction','update_gamma','spectralfilter','myfilter', &
            'err_curl_tolabs','err_curl_tolrel', &
            'maxcutback','polaralpha','polarbeta')
         call IO_warning(40_pInt,ext_msg=tag)
#endif

!--------------------------------------------------------------------------------------------------
! FEM parameters
#ifdef FEM
       case ('err_struct_tolabs')
         err_struct_tolAbs = IO_floatValue(line,positions,2_pInt)
       case ('err_struct_tolrel')
         err_struct_tolRel = IO_floatValue(line,positions,2_pInt)
       case ('err_thermal_tol')
         err_thermal_tol = IO_floatValue(line,positions,2_pInt)
       case ('err_damage_tol')
         err_damage_tol = IO_floatValue(line,positions,2_pInt)
       case ('err_vacancydiffusion_tol')
         err_vacancyDiffusion_tol = IO_floatValue(line,positions,2_pInt)
       case ('itmaxfem')
         itmaxFEM = IO_intValue(line,positions,2_pInt)
       case ('itminfem')
         itminFEM = IO_intValue(line,positions,2_pInt)
       case ('maxcutbackfem')
         maxCutBackFEM = IO_intValue(line,positions,2_pInt)
       case ('maxstaggerediter')
         stagItMax = IO_intValue(line,positions,2_pInt)
       case ('integrationorder')
         integrationorder = IO_intValue(line,positions,2_pInt)
       case ('structorder')
         structorder = IO_intValue(line,positions,2_pInt)
       case ('thermalorder')
         thermalorder = IO_intValue(line,positions,2_pInt)
       case ('damageorder')
         damageorder = IO_intValue(line,positions,2_pInt)
       case ('vacancydiffusionorder')
         vacancyDiffusionOrder = IO_intValue(line,positions,2_pInt)
       case ('petsc_optionsfem')
         petsc_optionsFEM = trim(line(positions(4):))
       case ('bbarstabilisation')
         BBarStabilisation = IO_intValue(line,positions,2_pInt) > 0_pInt
#else
      case ('err_struct_tolabs','err_struct_tolrel','err_thermal_tol','err_damage_tol','err_vacancydiffusion_tol',&         ! found FEM parameter for spectral/Abaqus/Marc build
            'itmaxfem', 'itminfem','maxcutbackfem','maxstaggerediter','integrationorder',&
            'structorder','thermalorder', 'damageorder','petsc_optionsfem','bbarstabilisation')
         call IO_warning(40_pInt,ext_msg=tag)
#endif
       case default                                                                                ! found unknown keyword
         call IO_error(300_pInt,ext_msg=tag)
     endselect
   enddo
   close(FILEUNIT)

 else fileExists
#ifdef FEM
   if (worldrank == 0) then
#endif  
   write(6,'(a,/)') ' using standard values'
   flush(6)
#ifdef FEM
   endif
#endif  
 endif fileExists

#ifdef Spectral
 select case(IO_lc(fftw_plan_mode))                                                                ! setting parameters for the plan creation of FFTW. Basically a translation from fftw3.f
   case('estimate','fftw_estimate')                                                                ! ordered from slow execution (but fast plan creation) to fast execution
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

 numerics_timeSyncing = numerics_timeSyncing .and. all(numerics_integrator==2_pInt)                 ! timeSyncing only allowed for explicit Euler integrator

!--------------------------------------------------------------------------------------------------
! writing parameters to output
 mainProcess3: if (worldrank == 0) then
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
   write(6,'(a24,1x,es8.1)')  ' rTol_crystalliteStress: ',rTol_crystalliteStress
   write(6,'(a24,1x,es8.1)')  ' aTol_crystalliteStress: ',aTol_crystalliteStress
   write(6,'(a24,2(1x,i8))')  ' integrator:             ',numerics_integrator
   write(6,'(a24,1x,L8)')     ' timeSyncing:            ',numerics_timeSyncing
   write(6,'(a24,1x,L8)')     ' analytic Jacobian:      ',analyticJaco
   write(6,'(a24,1x,L8)')     ' use ping pong scheme:   ',usepingpong
   write(6,'(a24,1x,es8.1,/)')' unitlength:             ',numerics_unitlength

   write(6,'(a24,1x,i8)')     ' nHomog:                 ',nHomog
   write(6,'(a24,1x,es8.1)')  ' subStepMinHomog:        ',subStepMinHomog
   write(6,'(a24,1x,es8.1)')  ' subStepSizeHomog:       ',subStepSizeHomog
   write(6,'(a24,1x,es8.1)')  ' stepIncreaseHomog:      ',stepIncreaseHomog
   write(6,'(a24,1x,i8,/)')   ' nMPstate:               ',nMPstate

!--------------------------------------------------------------------------------------------------
! RGC parameters
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

!--------------------------------------------------------------------------------------------------
! Random seeding parameter
   write(6,'(a24,1x,i16,/)')    ' fixed_seed:            ',fixedSeed
   if (fixedSeed <= 0_pInt) &
     write(6,'(a,/)') ' No fixed Seed: Random is random!'

!--------------------------------------------------------------------------------------------------
! gradient parameter
   write(6,'(a24,1x,es8.1)')   ' charLength:             ',charLength
   write(6,'(a24,1x,es8.1)')   ' residualStiffness:      ',residualStiffness

!--------------------------------------------------------------------------------------------------
! openMP parameter
  !$  write(6,'(a24,1x,i8,/)')   ' number of threads:      ',DAMASK_NumThreadsInt

!--------------------------------------------------------------------------------------------------
! spectral parameters
#ifdef Spectral
   write(6,'(a24,1x,i8)')      ' itmax:                  ',itmax
   write(6,'(a24,1x,i8)')      ' itmin:                  ',itmin
   write(6,'(a24,1x,i8)')      ' maxCutBack:             ',maxCutBack
   write(6,'(a24,1x,i8)')      ' continueCalculation:    ',continueCalculation
   write(6,'(a24,1x,L8)')      ' memory_efficient:       ',memory_efficient
   write(6,'(a24,1x,i8)')      ' divergence_correction:  ',divergence_correction
   write(6,'(a24,1x,a)')       ' spectral filter:        ',trim(spectral_filter)
   if(fftw_timelimit<0.0_pReal) then
     write(6,'(a24,1x,L8)')    ' fftw_timelimit:         ',.false.
   else    
     write(6,'(a24,1x,es8.1)') ' fftw_timelimit:         ',fftw_timelimit
   endif
   write(6,'(a24,1x,a)')       ' fftw_plan_mode:         ',trim(fftw_plan_mode)
   write(6,'(a24,1x,i8)')      ' fftw_planner_flag:      ',fftw_planner_flag
   write(6,'(a24,1x,L8,/)')    ' update_gamma:           ',update_gamma
   write(6,'(a24,1x,es8.1)')   ' err_stress_tolAbs:      ',err_stress_tolAbs
   write(6,'(a24,1x,es8.1)')   ' err_stress_tolRel:      ',err_stress_tolRel
   write(6,'(a24,1x,es8.1)')   ' err_div_tolAbs:         ',err_div_tolAbs
   write(6,'(a24,1x,es8.1)')   ' err_div_tolRel:         ',err_div_tolRel
   write(6,'(a24,1x,es8.1)')   ' err_curl_tolAbs:        ',err_curl_tolAbs
   write(6,'(a24,1x,es8.1)')   ' err_curl_tolRel:        ',err_curl_tolRel
   write(6,'(a24,1x,es8.1)')   ' polarAlpha:             ',polarAlpha
   write(6,'(a24,1x,es8.1)')   ' polarBeta:              ',polarBeta
   write(6,'(a24,1x,a)')       ' spectral solver:        ',trim(spectral_solver)
   write(6,'(a24,1x,a)')       ' PETSc_options:          ',trim(petsc_options)
#endif

!--------------------------------------------------------------------------------------------------
! spectral parameters
#ifdef FEM
   write(6,'(a24,1x,i8)')      ' itmaxFEM:               ',itmaxFEM
   write(6,'(a24,1x,i8)')      ' itminFEM:               ',itminFEM
   write(6,'(a24,1x,i8)')      ' maxCutBackFEM:          ',maxCutBackFEM
   write(6,'(a24,1x,i8)')      ' maxStaggeredIter:       ',stagItMax
   write(6,'(a24,1x,i8)')      ' integrationOrder:       ',integrationOrder
   write(6,'(a24,1x,i8)')      ' structOrder:            ',structOrder
   write(6,'(a24,1x,i8)')      ' thermalOrder:           ',thermalOrder
   write(6,'(a24,1x,i8)')      ' damageOrder:            ',damageOrder
   write(6,'(a24,1x,i8)')      ' vacancyDiffusionOrder:  ',vacancyDiffusionOrder
   write(6,'(a24,1x,es8.1)')   ' err_struct_tolAbs:      ',err_struct_tolAbs
   write(6,'(a24,1x,es8.1)')   ' err_struct_tolRel:      ',err_struct_tolRel
   write(6,'(a24,1x,es8.1)')   ' err_thermal_tol:        ',err_thermal_tol
   write(6,'(a24,1x,es8.1)')   ' err_damage_tol:         ',err_damage_tol
   write(6,'(a24,1x,es8.1)')   ' err_vacancyDiff_tol:    ',err_vacancyDiffusion_tol
   write(6,'(a24,1x,a)')       ' PETSc_optionsFEM:       ',trim(petsc_optionsFEM)
   write(6,'(a24,1x,L8)')      ' B-Bar stabilisation:    ',BBarStabilisation
#endif
 endif mainProcess3
 

!--------------------------------------------------------------------------------------------------
! sanity checks
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
 if (rTol_crystalliteStress <= 0.0_pReal)  call IO_error(301_pInt,ext_msg='rTol_crystalliteStress')
 if (aTol_crystalliteStress <= 0.0_pReal)  call IO_error(301_pInt,ext_msg='aTol_crystalliteStress')
 if (any(numerics_integrator <= 0_pInt) .or. any(numerics_integrator >= 6_pInt)) &
                                           call IO_error(301_pInt,ext_msg='integrator')
 if (numerics_unitlength <= 0.0_pReal)     call IO_error(301_pInt,ext_msg='unitlength')
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
 if (residualStiffness <= 0.0_pReal)       call IO_error(301_pInt,ext_msg='residualStiffness')
#ifdef Spectral
 if (itmax <= 1_pInt)                      call IO_error(301_pInt,ext_msg='itmax')
 if (itmin > itmax .or. itmin < 1_pInt)    call IO_error(301_pInt,ext_msg='itmin')
 if (continueCalculation /= 0_pInt .and. &
     continueCalculation /= 1_pInt)        call IO_error(301_pInt,ext_msg='continueCalculation')
 if (divergence_correction < 0_pInt .or. &
     divergence_correction > 2_pInt)       call IO_error(301_pInt,ext_msg='divergence_correction')
 if (maxCutBack < 0_pInt)                  call IO_error(301_pInt,ext_msg='maxCutBack')
 if (update_gamma .and. &
                   .not. memory_efficient) call IO_error(error_ID = 847_pInt)
 if (err_stress_tolrel <= 0.0_pReal)       call IO_error(301_pInt,ext_msg='err_stress_tolRel')
 if (err_stress_tolabs <= 0.0_pReal)       call IO_error(301_pInt,ext_msg='err_stress_tolAbs')
 if (err_div_tolRel <= 0.0_pReal)          call IO_error(301_pInt,ext_msg='err_div_tolRel')
 if (err_div_tolAbs <= 0.0_pReal)          call IO_error(301_pInt,ext_msg='err_div_tolAbs')
 if (err_curl_tolRel <= 0.0_pReal)         call IO_error(301_pInt,ext_msg='err_curl_tolRel')
 if (err_curl_tolAbs <= 0.0_pReal)         call IO_error(301_pInt,ext_msg='err_curl_tolAbs')
 if (polarAlpha <= 0.0_pReal .or. &
     polarAlpha >  2.0_pReal)              call IO_error(301_pInt,ext_msg='polarAlpha')
 if (polarBeta < 0.0_pReal .or. &
     polarBeta >  2.0_pReal)               call IO_error(301_pInt,ext_msg='polarBeta')
#endif
#ifdef FEM
 if (itmaxFEM <= 1_pInt)                   call IO_error(301_pInt,ext_msg='itmaxFEM')
 if (itminFEM > itmaxFEM .or. &
     itminFEM < 0_pInt)                    call IO_error(301_pInt,ext_msg='itminFEM')
 if (maxCutBackFEM < 0_pInt)               call IO_error(301_pInt,ext_msg='maxCutBackFEM')
 if (stagItMax < 0_pInt)                   call IO_error(301_pInt,ext_msg='maxStaggeredIter')
 if (err_struct_tolRel <= 0.0_pReal)       call IO_error(301_pInt,ext_msg='err_struct_tolRel')
 if (err_struct_tolAbs <= 0.0_pReal)       call IO_error(301_pInt,ext_msg='err_struct_tolAbs')
 if (err_thermal_tol <= 0.0_pReal)         call IO_error(301_pInt,ext_msg='err_thermal_tol')
 if (err_damage_tol <= 0.0_pReal)          call IO_error(301_pInt,ext_msg='err_damage_tol')
 if (err_vacancyDiffusion_tol <= 0.0_pReal)call IO_error(301_pInt,ext_msg='err_vacancydiffusion_tol')
#endif

end subroutine numerics_init

end module numerics
