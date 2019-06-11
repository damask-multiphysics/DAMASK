!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Managing of parameters related to numerics
!--------------------------------------------------------------------------------------------------
module numerics
 use prec

 implicit none
 private

 integer, protected, public :: &
   iJacoStiffness             =  1, &                                                              !< frequency of stiffness update
   nMPstate                   = 10, &                                                              !< materialpoint state loop limit
   randomSeed                 =  0, &                                                              !< fixed seeding for pseudo-random number generator, Default 0: use random seed
   worldrank                  =  0, &                                                              !< MPI worldrank (/=0 for MPI simulations only)
   worldsize                  =  1, &                                                              !< MPI worldsize (/=1 for MPI simulations only)
   numerics_integrator        =  1                                                                 !< method used for state integration Default 1: fix-point iteration
 integer(4), protected, public :: &
   DAMASK_NumThreadsInt       =  0                                                                  !< value stored in environment variable DAMASK_NUM_THREADS, set to zero if no OpenMP directive
 real(pReal), protected, public :: &
   defgradTolerance           =  1.0e-7_pReal, &                                                    !< deviation of deformation gradient that is still allowed (used by CPFEM to determine outdated ffn1)
   subStepMinHomog            =  1.0e-3_pReal, &                                                    !< minimum (relative) size of sub-step allowed during cutback in homogenization
   subStepSizeHomog           =  0.25_pReal, &                                                      !< size of first substep when cutback in homogenization
   stepIncreaseHomog          =  1.5_pReal, &                                                       !< increase of next substep size when previous substep converged in homogenization
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
   usePingPong                = .true.

!--------------------------------------------------------------------------------------------------
! field parameters:
 real(pReal), protected, public :: &
   err_struct_tolAbs          =  1.0e-10_pReal, &                                                   !< absolute tolerance for mechanical equilibrium
   err_struct_tolRel          =  1.0e-4_pReal, &                                                    !< relative tolerance for mechanical equilibrium
   err_thermal_tolAbs         =  1.0e-2_pReal, &                                                    !< absolute tolerance for thermal equilibrium
   err_thermal_tolRel         =  1.0e-6_pReal, &                                                    !< relative tolerance for thermal equilibrium
   err_damage_tolAbs          =  1.0e-2_pReal, &                                                    !< absolute tolerance for damage evolution
   err_damage_tolRel          =  1.0e-6_pReal                                                       !< relative tolerance for damage evolution
 integer, protected, public :: &
   itmax                      =  250, &                                                             !< maximum number of iterations
   itmin                      =  1, &                                                               !< minimum number of iterations
   stagItMax                  =  10, &                                                              !< max number of field level staggered iterations
   maxCutBack                 =  3                                                                  !< max number of cut backs

!--------------------------------------------------------------------------------------------------
! spectral parameters:
#ifdef Grid
 real(pReal), protected, public :: &
   err_div_tolAbs             =  1.0e-4_pReal, &                                                    !< absolute tolerance for equilibrium
   err_div_tolRel             =  5.0e-4_pReal, &                                                    !< relative tolerance for equilibrium
   err_curl_tolAbs            =  1.0e-10_pReal, &                                                   !< absolute tolerance for compatibility
   err_curl_tolRel            =  5.0e-4_pReal, &                                                    !< relative tolerance for compatibility
   err_stress_tolAbs          =  1.0e3_pReal,  &                                                    !< absolute tolerance for fullfillment of stress BC
   err_stress_tolRel          =  0.01_pReal, &                                                      !< relative tolerance for fullfillment of stress BC
   rotation_tol               =  1.0e-12_pReal, &                                                   !< tolerance of rotation specified in loadcase, Default 1.0e-12: first guess
   polarAlpha                 =  1.0_pReal, &                                                       !< polarization scheme parameter 0.0 < alpha < 2.0. alpha = 1.0 ==> AL scheme, alpha = 2.0 ==> accelerated scheme 
   polarBeta                  =  1.0_pReal                                                          !< polarization scheme parameter 0.0 < beta < 2.0. beta = 1.0 ==> AL scheme, beta = 2.0 ==> accelerated scheme 
 character(len=1024), protected, public :: &
   petsc_defaultOptions       = '-mech_snes_type ngmres &
                                &-damage_snes_type ngmres &
                                &-thermal_snes_type ngmres ', &
   petsc_options              = ''
 logical, protected, public :: &
   continueCalculation        = .false.                                                          !< false:exit if BVP solver does not converge, true: continue calculation despite BVP solver not converging

#endif

!--------------------------------------------------------------------------------------------------
! FEM parameters:
#ifdef FEM
 integer, protected, public :: &
   integrationOrder           =  2, &                                                              !< order of quadrature rule required
   structOrder                =  2                                                                 !< order of displacement shape functions
 logical, protected, public :: & 
   BBarStabilisation          = .false.                                                  
 character(len=4096), protected, public :: &
   petsc_defaultOptions    = '-mech_snes_type newtonls &
                             &-mech_snes_linesearch_type cp &
                             &-mech_snes_ksp_ew &
                             &-mech_snes_ksp_ew_rtol0 0.01 &
                             &-mech_snes_ksp_ew_rtolmax 0.01 &
                             &-mech_ksp_type fgmres &
                             &-mech_ksp_max_it 25 &
                             &-mech_pc_type ml &
                             &-mech_mg_levels_ksp_type chebyshev &
                             &-mech_mg_levels_pc_type sor &
                             &-mech_pc_ml_nullspace user ',&
   petsc_options           = ''
#endif

 public :: numerics_init
  
contains


!--------------------------------------------------------------------------------------------------
!> @brief reads in parameters from numerics.config and sets openMP related parameters. Also does
! a sanity check
!--------------------------------------------------------------------------------------------------
subroutine numerics_init
 use IO, only: &
   IO_read_ASCII, &
   IO_error, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_lc, &
   IO_floatValue, &
   IO_intValue, &
   IO_warning
#ifdef PETSc
#include <petsc/finclude/petscsys.h>
   use petscsys
#endif
!$ use OMP_LIB, only: omp_set_num_threads
!$ integer ::                                gotDAMASK_NUM_THREADS = 1
 integer :: i,j, ierr
 integer, allocatable, dimension(:) :: chunkPos
 character(len=pStringLen), dimension(:), allocatable :: fileContent
 character(len=pStringLen) :: &
   tag ,&
   line
 logical :: fexist
!$ character(len=6) DAMASK_NumThreadsString                                                         ! environment variable DAMASK_NUM_THREADS

#ifdef PETSc
 call MPI_Comm_rank(PETSC_COMM_WORLD,worldrank,ierr);CHKERRQ(ierr)
 call MPI_Comm_size(PETSC_COMM_WORLD,worldsize,ierr);CHKERRQ(ierr)
#endif
 write(6,'(/,a)') ' <<<+-  numerics init  -+>>>'
 
!$ call GET_ENVIRONMENT_VARIABLE(NAME='DAMASK_NUM_THREADS',VALUE=DAMASK_NumThreadsString,STATUS=gotDAMASK_NUM_THREADS)   ! get environment variable DAMASK_NUM_THREADS...
!$ if(gotDAMASK_NUM_THREADS /= 0) then                                                              ! could not get number of threads, set it to 1
!$   call IO_warning(35,ext_msg='BEGIN:'//DAMASK_NumThreadsString//':END')
!$   DAMASK_NumThreadsInt = 1_4
!$ else
!$   read(DAMASK_NumThreadsString,'(i6)') DAMASK_NumThreadsInt                                      ! read as integer
!$   if (DAMASK_NumThreadsInt < 1_4) DAMASK_NumThreadsInt = 1_4                                     ! in case of string conversion fails, set it to one
!$ endif
!$ call omp_set_num_threads(DAMASK_NumThreadsInt)                                                   ! set number of threads for parallel execution
 
 inquire(file='numerics.config', exist=fexist)
 
 fileExists: if (fexist) then
   write(6,'(a,/)') ' using values from config file'
   flush(6)
   fileContent = IO_read_ASCII('numerics.config') 
   do j=1, size(fileContent)
   
!--------------------------------------------------------------------------------------------------
! read variables from config file and overwrite default parameters if keyword is present
     line = fileContent(j)
     do i=1,len(line)
       if(line(i:i) == '=') line(i:i) = ' '                                                         ! also allow keyword = value version
     enddo
     if (IO_isBlank(line)) cycle                                                                    ! skip empty lines
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1))                                                   ! extract key

     select case(tag)
       case ('defgradtolerance')
         defgradTolerance = IO_floatValue(line,chunkPos,2)
       case ('ijacostiffness')
         iJacoStiffness = IO_intValue(line,chunkPos,2)
       case ('nmpstate')
         nMPstate = IO_intValue(line,chunkPos,2)
       case ('substepminhomog')
         subStepMinHomog = IO_floatValue(line,chunkPos,2)
       case ('substepsizehomog')
         subStepSizeHomog = IO_floatValue(line,chunkPos,2)
       case ('stepincreasehomog')
         stepIncreaseHomog = IO_floatValue(line,chunkPos,2)
       case ('integrator')
         numerics_integrator = IO_intValue(line,chunkPos,2)
       case ('usepingpong')
         usepingpong = IO_intValue(line,chunkPos,2) > 0
       case ('unitlength')
         numerics_unitlength = IO_floatValue(line,chunkPos,2)

!--------------------------------------------------------------------------------------------------
! RGC parameters
       case ('atol_rgc')
         absTol_RGC = IO_floatValue(line,chunkPos,2)
       case ('rtol_rgc')
         relTol_RGC = IO_floatValue(line,chunkPos,2)
       case ('amax_rgc')
         absMax_RGC = IO_floatValue(line,chunkPos,2)
       case ('rmax_rgc')
         relMax_RGC = IO_floatValue(line,chunkPos,2)
       case ('perturbpenalty_rgc')
         pPert_RGC = IO_floatValue(line,chunkPos,2)
       case ('relevantmismatch_rgc')
         xSmoo_RGC = IO_floatValue(line,chunkPos,2)
       case ('viscositypower_rgc')
         viscPower_RGC = IO_floatValue(line,chunkPos,2)
       case ('viscositymodulus_rgc')
         viscModus_RGC = IO_floatValue(line,chunkPos,2)
       case ('refrelaxationrate_rgc')
         refRelaxRate_RGC = IO_floatValue(line,chunkPos,2)
       case ('maxrelaxation_rgc')
         maxdRelax_RGC = IO_floatValue(line,chunkPos,2)
       case ('maxvoldiscrepancy_rgc')
         maxVolDiscr_RGC = IO_floatValue(line,chunkPos,2)
       case ('voldiscrepancymod_rgc')
         volDiscrMod_RGC = IO_floatValue(line,chunkPos,2)
       case ('discrepancypower_rgc')
         volDiscrPow_RGC = IO_floatValue(line,chunkPos,2)

!--------------------------------------------------------------------------------------------------
! random seeding parameter
       case ('random_seed','fixed_seed')
         randomSeed = IO_intValue(line,chunkPos,2)

!--------------------------------------------------------------------------------------------------
! gradient parameter
       case ('charlength')
         charLength = IO_floatValue(line,chunkPos,2)
       case ('residualstiffness')
         residualStiffness = IO_floatValue(line,chunkPos,2)

!--------------------------------------------------------------------------------------------------
! field parameters
       case ('err_struct_tolabs')
         err_struct_tolAbs = IO_floatValue(line,chunkPos,2)
       case ('err_struct_tolrel')
         err_struct_tolRel = IO_floatValue(line,chunkPos,2)
       case ('err_thermal_tolabs')
         err_thermal_tolabs = IO_floatValue(line,chunkPos,2)
       case ('err_thermal_tolrel')
         err_thermal_tolrel = IO_floatValue(line,chunkPos,2)
       case ('err_damage_tolabs')
         err_damage_tolabs = IO_floatValue(line,chunkPos,2)
       case ('err_damage_tolrel')
         err_damage_tolrel = IO_floatValue(line,chunkPos,2)
       case ('itmax')
         itmax = IO_intValue(line,chunkPos,2)
       case ('itmin')
         itmin = IO_intValue(line,chunkPos,2)
       case ('maxcutback')
         maxCutBack = IO_intValue(line,chunkPos,2)
       case ('maxstaggerediter')
         stagItMax = IO_intValue(line,chunkPos,2)

!--------------------------------------------------------------------------------------------------
! spectral parameters
#ifdef Grid
       case ('err_div_tolabs')
         err_div_tolAbs = IO_floatValue(line,chunkPos,2)
       case ('err_div_tolrel')
         err_div_tolRel = IO_floatValue(line,chunkPos,2)
       case ('err_stress_tolrel')
         err_stress_tolrel = IO_floatValue(line,chunkPos,2)
       case ('err_stress_tolabs')
         err_stress_tolabs = IO_floatValue(line,chunkPos,2)
       case ('continuecalculation')
         continueCalculation = IO_intValue(line,chunkPos,2) > 0
       case ('petsc_options')
         petsc_options = trim(line(chunkPos(4):))
       case ('err_curl_tolabs')
         err_curl_tolAbs = IO_floatValue(line,chunkPos,2)
       case ('err_curl_tolrel')
         err_curl_tolRel = IO_floatValue(line,chunkPos,2)
       case ('polaralpha')
         polarAlpha = IO_floatValue(line,chunkPos,2)
       case ('polarbeta')
         polarBeta = IO_floatValue(line,chunkPos,2)
#endif

!--------------------------------------------------------------------------------------------------
! FEM parameters
#ifdef FEM
       case ('integrationorder')
         integrationorder = IO_intValue(line,chunkPos,2)
       case ('structorder')
         structorder = IO_intValue(line,chunkPos,2)
       case ('petsc_options')
         petsc_options = trim(line(chunkPos(4):))
       case ('bbarstabilisation')
         BBarStabilisation = IO_intValue(line,chunkPos,2) > 0
#endif
     end select
   enddo


 else fileExists
   write(6,'(a,/)') ' using standard values'
   flush(6)
 endif fileExists


!--------------------------------------------------------------------------------------------------
! writing parameters to output
 write(6,'(a24,1x,es8.1)')  ' defgradTolerance:       ',defgradTolerance
 write(6,'(a24,1x,i8)')     ' iJacoStiffness:         ',iJacoStiffness
 write(6,'(a24,1x,i8)')     ' integrator:             ',numerics_integrator
 write(6,'(a24,1x,L8)')     ' use ping pong scheme:   ',usepingpong
 write(6,'(a24,1x,es8.1,/)')' unitlength:             ',numerics_unitlength

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
 write(6,'(a16,1x,i16,/)')    ' random_seed:    ',randomSeed
 if (randomSeed <= 0) &
   write(6,'(a,/)')           ' random seed will be generated!'

!--------------------------------------------------------------------------------------------------
! gradient parameter
 write(6,'(a24,1x,es8.1)')   ' charLength:             ',charLength
 write(6,'(a24,1x,es8.1)')   ' residualStiffness:      ',residualStiffness

!--------------------------------------------------------------------------------------------------
! openMP parameter
 !$  write(6,'(a24,1x,i8,/)')   ' number of threads:      ',DAMASK_NumThreadsInt

!--------------------------------------------------------------------------------------------------
! field parameters
 write(6,'(a24,1x,i8)')      ' itmax:                  ',itmax
 write(6,'(a24,1x,i8)')      ' itmin:                  ',itmin
 write(6,'(a24,1x,i8)')      ' maxCutBack:             ',maxCutBack
 write(6,'(a24,1x,i8)')      ' maxStaggeredIter:       ',stagItMax
 write(6,'(a24,1x,es8.1)')   ' err_struct_tolAbs:      ',err_struct_tolAbs
 write(6,'(a24,1x,es8.1)')   ' err_struct_tolRel:      ',err_struct_tolRel
 write(6,'(a24,1x,es8.1)')   ' err_thermal_tolabs:     ',err_thermal_tolabs
 write(6,'(a24,1x,es8.1)')   ' err_thermal_tolrel:     ',err_thermal_tolrel
 write(6,'(a24,1x,es8.1)')   ' err_damage_tolabs:      ',err_damage_tolabs
 write(6,'(a24,1x,es8.1)')   ' err_damage_tolrel:      ',err_damage_tolrel

!--------------------------------------------------------------------------------------------------
! spectral parameters
#ifdef Grid
 write(6,'(a24,1x,L8)')      ' continueCalculation:    ',continueCalculation
 write(6,'(a24,1x,es8.1)')   ' err_stress_tolAbs:      ',err_stress_tolAbs
 write(6,'(a24,1x,es8.1)')   ' err_stress_tolRel:      ',err_stress_tolRel
 write(6,'(a24,1x,es8.1)')   ' err_div_tolAbs:         ',err_div_tolAbs
 write(6,'(a24,1x,es8.1)')   ' err_div_tolRel:         ',err_div_tolRel
 write(6,'(a24,1x,es8.1)')   ' err_curl_tolAbs:        ',err_curl_tolAbs
 write(6,'(a24,1x,es8.1)')   ' err_curl_tolRel:        ',err_curl_tolRel
 write(6,'(a24,1x,es8.1)')   ' polarAlpha:             ',polarAlpha
 write(6,'(a24,1x,es8.1)')   ' polarBeta:              ',polarBeta
 write(6,'(a24,1x,a)')       ' PETSc_options:          ',trim(petsc_defaultOptions)//' '//trim(petsc_options)
#endif

!--------------------------------------------------------------------------------------------------
! spectral parameters
#ifdef FEM
 write(6,'(a24,1x,i8)')      ' integrationOrder:       ',integrationOrder
 write(6,'(a24,1x,i8)')      ' structOrder:            ',structOrder
 write(6,'(a24,1x,a)')       ' PETSc_options:          ',trim(petsc_defaultOptions)//' '//trim(petsc_options)
 write(6,'(a24,1x,L8)')      ' B-Bar stabilisation:    ',BBarStabilisation
#endif
 

!--------------------------------------------------------------------------------------------------
! sanity checks
 if (defgradTolerance <= 0.0_pReal)        call IO_error(301,ext_msg='defgradTolerance')
 if (iJacoStiffness < 1)                   call IO_error(301,ext_msg='iJacoStiffness')
 if (nMPstate < 1)                         call IO_error(301,ext_msg='nMPstate')
 if (subStepMinHomog <= 0.0_pReal)         call IO_error(301,ext_msg='subStepMinHomog')
 if (subStepSizeHomog <= 0.0_pReal)        call IO_error(301,ext_msg='subStepSizeHomog')
 if (stepIncreaseHomog <= 0.0_pReal)       call IO_error(301,ext_msg='stepIncreaseHomog')
 if (numerics_integrator <= 0 .or. numerics_integrator >= 6) &
                                           call IO_error(301,ext_msg='integrator')
 if (numerics_unitlength <= 0.0_pReal)     call IO_error(301,ext_msg='unitlength')
 if (absTol_RGC <= 0.0_pReal)              call IO_error(301,ext_msg='absTol_RGC')
 if (relTol_RGC <= 0.0_pReal)              call IO_error(301,ext_msg='relTol_RGC')
 if (absMax_RGC <= 0.0_pReal)              call IO_error(301,ext_msg='absMax_RGC')
 if (relMax_RGC <= 0.0_pReal)              call IO_error(301,ext_msg='relMax_RGC')
 if (pPert_RGC <= 0.0_pReal)               call IO_error(301,ext_msg='pPert_RGC')
 if (xSmoo_RGC <= 0.0_pReal)               call IO_error(301,ext_msg='xSmoo_RGC')
 if (viscPower_RGC < 0.0_pReal)            call IO_error(301,ext_msg='viscPower_RGC')
 if (viscModus_RGC < 0.0_pReal)            call IO_error(301,ext_msg='viscModus_RGC')
 if (refRelaxRate_RGC <= 0.0_pReal)        call IO_error(301,ext_msg='refRelaxRate_RGC')
 if (maxdRelax_RGC <= 0.0_pReal)           call IO_error(301,ext_msg='maxdRelax_RGC')
 if (maxVolDiscr_RGC <= 0.0_pReal)         call IO_error(301,ext_msg='maxVolDiscr_RGC')
 if (volDiscrMod_RGC < 0.0_pReal)          call IO_error(301,ext_msg='volDiscrMod_RGC')
 if (volDiscrPow_RGC <= 0.0_pReal)         call IO_error(301,ext_msg='volDiscrPw_RGC')
 if (residualStiffness < 0.0_pReal)        call IO_error(301,ext_msg='residualStiffness')
 if (itmax <= 1)                           call IO_error(301,ext_msg='itmax')
 if (itmin > itmax .or. itmin < 1)         call IO_error(301,ext_msg='itmin')
 if (maxCutBack < 0)                       call IO_error(301,ext_msg='maxCutBack')
 if (stagItMax < 0)                        call IO_error(301,ext_msg='maxStaggeredIter')
 if (err_struct_tolRel <= 0.0_pReal)       call IO_error(301,ext_msg='err_struct_tolRel')
 if (err_struct_tolAbs <= 0.0_pReal)       call IO_error(301,ext_msg='err_struct_tolAbs')
 if (err_thermal_tolabs <= 0.0_pReal)      call IO_error(301,ext_msg='err_thermal_tolabs')
 if (err_thermal_tolrel <= 0.0_pReal)      call IO_error(301,ext_msg='err_thermal_tolrel')
 if (err_damage_tolabs <= 0.0_pReal)       call IO_error(301,ext_msg='err_damage_tolabs')
 if (err_damage_tolrel <= 0.0_pReal)       call IO_error(301,ext_msg='err_damage_tolrel')
#ifdef Grid
 if (err_stress_tolrel <= 0.0_pReal)       call IO_error(301,ext_msg='err_stress_tolRel')
 if (err_stress_tolabs <= 0.0_pReal)       call IO_error(301,ext_msg='err_stress_tolAbs')
 if (err_div_tolRel < 0.0_pReal)           call IO_error(301,ext_msg='err_div_tolRel')
 if (err_div_tolAbs <= 0.0_pReal)          call IO_error(301,ext_msg='err_div_tolAbs')
 if (err_curl_tolRel < 0.0_pReal)          call IO_error(301,ext_msg='err_curl_tolRel')
 if (err_curl_tolAbs <= 0.0_pReal)         call IO_error(301,ext_msg='err_curl_tolAbs')
 if (polarAlpha <= 0.0_pReal .or. &
     polarAlpha >  2.0_pReal)              call IO_error(301,ext_msg='polarAlpha')
 if (polarBeta < 0.0_pReal .or. &
     polarBeta > 2.0_pReal)                call IO_error(301,ext_msg='polarBeta')
#endif

end subroutine numerics_init

end module numerics
