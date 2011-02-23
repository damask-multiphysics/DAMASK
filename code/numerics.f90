!* $Id$
!##############################################################
MODULE numerics
!##############################################################

use prec, only: pInt, pReal
implicit none

character(len=64), parameter :: numerics_configFile = 'numerics.config' ! name of configuration file
integer(pInt)                   iJacoStiffness, &                       ! frequency of stiffness update
                                iJacoLpresiduum, &                      ! frequency of Jacobian update of residuum in Lp
                                nHomog, &                               ! homogenization loop limit (only for debugging info, loop limit is determined by "subStepMinHomog")
                                nMPstate, &                             ! materialpoint state loop limit
                                nCryst, &                               ! crystallite loop limit (only for debugging info, loop limit is determined by "subStepMinCryst")
                                nState, &                               ! state loop limit
                                nStress, &                              ! stress loop limit
                                pert_method, &                          ! method used in perturbation technique for tangent
                                numerics_integrationMode                ! integration mode 1 = central solution ; integration mode 2 = perturbation
integer(pInt), dimension(2) ::  numerics_integrator                     ! method used for state integration (central & perturbed state)
real(pReal)                     relevantStrain, &                       ! strain increment considered significant (used by crystallite to determine whether strain inc is considered significant)
                                defgradTolerance, &                     ! deviation of deformation gradient that is still allowed (used by CPFEM to determine outdated ffn1)
                                pert_Fg, &                              ! strain perturbation for FEM Jacobi
                                subStepMinCryst, &                      ! minimum (relative) size of sub-step allowed during cutback in crystallite
                                subStepMinHomog, &                      ! minimum (relative) size of sub-step allowed during cutback in homogenization
                                subStepSizeCryst, &                     ! size of first substep when cutback in crystallite
                                subStepSizeHomog, &                     ! size of first substep when cutback in homogenization
                                stepIncreaseCryst, &                    ! increase of next substep size when previous substep converged in crystallite
                                stepIncreaseHomog, &                    ! increase of next substep size when previous substep converged in homogenization
                                rTol_crystalliteState, &                ! relative tolerance in crystallite state loop 
                                rTol_crystalliteTemperature, &          ! relative tolerance in crystallite temperature loop 
                                rTol_crystalliteStress, &               ! relative tolerance in crystallite stress loop
                                aTol_crystalliteStress, &               ! absolute tolerance in crystallite stress loop

!* RGC parameters: added <<<updated 17.12.2009>>>
                                absTol_RGC, &                           ! absolute tolerance of RGC residuum
                                relTol_RGC, &                           ! relative tolerance of RGC residuum
                                absMax_RGC, &                           ! absolute maximum of RGC residuum
                                relMax_RGC, &                           ! relative maximum of RGC residuum
                                pPert_RGC, &                            ! perturbation for computing RGC penalty tangent
                                xSmoo_RGC, &                            ! RGC penalty smoothing parameter (hyperbolic tangent)
                                viscPower_RGC, &                        ! power (sensitivity rate) of numerical viscosity in RGC scheme
                                viscModus_RGC, &                        ! stress modulus of RGC numerical viscosity
                                refRelaxRate_RGC, &                     ! reference relaxation rate in RGC viscosity
                                maxdRelax_RGC, &                        ! threshold of maximum relaxation vector increment (if exceed this then cutback)
                                maxVolDiscr_RGC, &                      ! threshold of maximum volume discrepancy allowed
                                volDiscrMod_RGC, &                      ! stiffness of RGC volume discrepancy (zero = without volume discrepancy constraint)
                                volDiscrPow_RGC, &                      ! powerlaw penalty for volume discrepancy
!* spectral parameters:
                                err_div_tol, &                          ! error of divergence in fourier space
                                err_stress_tol, &                       ! absolut stress error, will be computed from err_stress_tolrel (dont prescribe a value)
                                err_stress_tolrel, &                    ! factor to multiply with highest stress to get err_stress_tol
                                err_defgrad_tol                         ! tolerance for error of defgrad compared to prescribed defgrad
logical                         memory_efficient                          ! for fast execution (pre calculation of gamma_hat)
integer(pInt)                   itmax , &                               ! maximum number of iterations


!* Random seeding parameters
                                fixedSeed                               ! fixed seeding for pseudo-random number generator
!* OpenMP variable
!$ integer(pInt)                mpieNumThreadsInt                       ! value stored in environment variable MPIE_NUM_THREADS


CONTAINS
 
!*******************************************
!    initialization subroutine
!*******************************************
subroutine numerics_init()
  
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

  !*** input variables ***!
  
  !*** output variables ***!
  
  !*** local variables ***!
  integer(pInt), parameter ::                 fileunit = 300  
  integer(pInt), parameter ::                 maxNchunks = 2
  integer(pInt), dimension(1+2*maxNchunks) :: positions
  character(len=64)                           tag
  character(len=1024)                         line
  
! OpenMP variable
!$ character(len=4) mpieNumThreadsString                               !environment variable MPIE_NUMTHREADS

  write(6,*)
  write(6,*) '<<<+-  numerics init  -+>>>'
  write(6,*) '$Id$'
  write(6,*)
  
  ! initialize all parameters with standard values
  relevantStrain          = 1.0e-7_pReal
  defgradTolerance        = 1.0e-7_pReal
  iJacoStiffness          = 1_pInt
  iJacoLpresiduum         = 1_pInt
  pert_Fg                 = 1.0e-7_pReal
  pert_method             = 1
  nHomog                  = 20_pInt
  subStepMinHomog         = 1.0e-3_pReal
  subStepSizeHomog        = 0.25
  stepIncreaseHomog       = 1.5
  nMPstate                = 10_pInt
  nCryst                  = 20_pInt
  subStepMinCryst         = 1.0e-3_pReal
  subStepsizeCryst        = 0.25
  stepIncreaseCryst       = 1.5
  nState                  = 10_pInt
  nStress                 = 40_pInt
  rTol_crystalliteState   = 1.0e-6_pReal
  rTol_crystalliteTemperature = 1.0e-6_pReal
  rTol_crystalliteStress  = 1.0e-6_pReal
  aTol_crystalliteStress  = 1.0e-8_pReal            ! residuum is in Lp (hence strain on the order of 1e-8 here)
  numerics_integrator(1)  = 1                       ! fix-point iteration
  numerics_integrator(2)  = 1                       ! fix-point iteration
  
!* RGC parameters: added <<<updated 17.12.2009>>> with moderate setting
  absTol_RGC              = 1.0e+4
  relTol_RGC              = 1.0e-3
  absMax_RGC              = 1.0e+10
  relMax_RGC              = 1.0e+2
  pPert_RGC               = 1.0e-7
  xSmoo_RGC               = 1.0e-5
  viscPower_RGC           = 1.0e+0  ! Newton viscosity (linear model)
  viscModus_RGC           = 0.0e+0  ! No viscosity is applied
  refRelaxRate_RGC        = 1.0e-3
  maxdRelax_RGC           = 1.0e+0
  maxVolDiscr_RGC         = 1.0e-5  ! tolerance for volume discrepancy allowed
  volDiscrMod_RGC         = 1.0e+12
  volDiscrPow_RGC         = 5.0

!* spectral parameters:  
  err_div_tol             = 1.0e-4   ! proposed by Suquet, less strict criteria are usefull, e.g. 5e-3
  err_defgrad_tol         = 1.0e-3   ! relative tolerance for fullfillment of average deformation gradient (is usually passively fullfilled)
  err_stress_tolrel       = 0.01     ! relative tolerance for fullfillment of stress BC
  itmax                   = 20_pInt  ! Maximum iteration number
  memory_efficient        = .true.   ! Precalculate Gamma-operator (81 double per point)

!* Random seeding parameters: added <<<updated 27.08.2009>>>
  fixedSeed               = 0_pInt


!* determin number of threads from environment variable MPIE_NUM_THREADS
!$ call GetEnv('MPIE_NUM_THREADS',mpieNumThreadsString)                     ! get environment variable MPIE_NUM_THREADS...
!$ read(mpieNumThreadsString,'(i4)') mpieNumThreadsInt                      ! ...convert it to integer...
!$ if (mpieNumThreadsInt < 1) mpieNumThreadsInt = 1                         ! ...ensure that its at least one...
!$ call omp_set_num_threads(mpieNumThreadsInt)                              ! ...and use it as number of threads for parallel execution

  ! try to open the config file
  if(IO_open_file(fileunit,numerics_configFile)) then 
  
    write(6,*) '   ... using values from config file'
    write(6,*)
    
    line = ''
    ! read variables from config file and overwrite parameters
    do
      read(fileunit,'(a1024)',END=100) line
      if (IO_isBlank(line)) cycle                           ! skip empty lines
      positions = IO_stringPos(line,maxNchunks)
      tag = IO_lc(IO_stringValue(line,positions,1))         ! extract key
      select case(tag)
        case ('relevantstrain')
              relevantStrain = IO_floatValue(line,positions,2)
        case ('defgradtolerance')
              defgradTolerance = IO_floatValue(line,positions,2)
        case ('ijacostiffness')
              iJacoStiffness = IO_intValue(line,positions,2)
        case ('ijacolpresiduum')
              iJacoLpresiduum = IO_intValue(line,positions,2)
        case ('pert_fg')
              pert_Fg = IO_floatValue(line,positions,2)
        case ('pert_method')
              pert_method = IO_intValue(line,positions,2)
        case ('nhomog')
              nHomog = IO_intValue(line,positions,2)
        case ('nmpstate')
              nMPstate = IO_intValue(line,positions,2)
        case ('ncryst')
              nCryst = IO_intValue(line,positions,2)
        case ('nstate')
              nState = IO_intValue(line,positions,2)
        case ('nstress')
              nStress = IO_intValue(line,positions,2)
        case ('substepmincryst')
              subStepMinCryst = IO_floatValue(line,positions,2)
        case ('substepsizecryst')
              subStepSizeCryst = IO_floatValue(line,positions,2)
        case ('stepincreasecryst')
              stepIncreaseCryst = IO_floatValue(line,positions,2)
        case ('substepminhomog')
              subStepMinHomog = IO_floatValue(line,positions,2)
        case ('substepsizehomog')
              subStepSizeHomog = IO_floatValue(line,positions,2)
        case ('stepincreasehomog')
              stepIncreaseHomog = IO_floatValue(line,positions,2)
        case ('rtol_crystallitestate')
              rTol_crystalliteState = IO_floatValue(line,positions,2)
        case ('rtol_crystallitetemperature')
              rTol_crystalliteTemperature = IO_floatValue(line,positions,2)
        case ('rtol_crystallitestress')
              rTol_crystalliteStress = IO_floatValue(line,positions,2)
        case ('atol_crystallitestress')
              aTol_crystalliteStress = IO_floatValue(line,positions,2)
        case ('integrator')
              numerics_integrator(1) = IO_intValue(line,positions,2)
        case ('integratorstiffness')
              numerics_integrator(2) = IO_intValue(line,positions,2)

!* RGC parameters: 
        case ('atol_rgc')
              absTol_RGC = IO_floatValue(line,positions,2)
        case ('rtol_rgc')
              relTol_RGC = IO_floatValue(line,positions,2)
        case ('amax_rgc')
              absMax_RGC = IO_floatValue(line,positions,2)
        case ('rmax_rgc')
              relMax_RGC = IO_floatValue(line,positions,2)
        case ('perturbpenalty_rgc')
              pPert_RGC = IO_floatValue(line,positions,2)
        case ('relevantmismatch_rgc')
              xSmoo_RGC = IO_floatValue(line,positions,2)
        case ('viscositypower_rgc')
              viscPower_RGC = IO_floatValue(line,positions,2)
        case ('viscositymodulus_rgc')
              viscModus_RGC = IO_floatValue(line,positions,2)
        case ('refrelaxationrate_rgc')
              refRelaxRate_RGC = IO_floatValue(line,positions,2)
        case ('maxrelaxation_rgc')
              maxdRelax_RGC = IO_floatValue(line,positions,2)
        case ('maxvoldiscrepancy_rgc')
              maxVolDiscr_RGC = IO_floatValue(line,positions,2)
        case ('voldiscrepancymod_rgc')
              volDiscrMod_RGC = IO_floatValue(line,positions,2)
        case ('discrepancypower_rgc')
              volDiscrPow_RGC = IO_floatValue(line,positions,2)

!* spectral parameters
        case ('err_div_tol')
              err_div_tol = IO_floatValue(line,positions,2)
        case ('err_defgrad_tol')
              err_defgrad_tol = IO_floatValue(line,positions,2)
        case ('err_stress_tolrel')
              err_stress_tolrel = IO_floatValue(line,positions,2)
        case ('itmax')
              itmax = IO_intValue(line,positions,2)
        case ('memory_efficient')
              memory_efficient = IO_intValue(line,positions,2)  > 0_pInt

!* Random seeding parameters
        case ('fixed_seed')
              fixedSeed = IO_floatValue(line,positions,2)
      endselect
    enddo
    100 close(fileunit)
  
  ! no config file, so we use standard values
  else 
  
    write(6,*) '   ... using standard values'
    write(6,*)
    
  endif

  ! writing parameters to output file
  write(6,'(a24,x,e8.1)') 'relevantStrain:         ',relevantStrain
  write(6,'(a24,x,e8.1)') 'defgradTolerance:       ',defgradTolerance
  write(6,'(a24,x,i8)')   'iJacoStiffness:         ',iJacoStiffness
  write(6,'(a24,x,i8)')   'iJacoLpresiduum:        ',iJacoLpresiduum
  write(6,'(a24,x,e8.1)') 'pert_Fg:                ',pert_Fg
  write(6,'(a24,x,i8)')   'pert_method:            ',pert_method
  write(6,'(a24,x,i8)')   'nCryst:                 ',nCryst
  write(6,'(a24,x,e8.1)') 'subStepMinCryst:        ',subStepMinCryst
  write(6,'(a24,x,e8.1)') 'subStepSizeCryst:       ',subStepSizeCryst
  write(6,'(a24,x,e8.1)') 'stepIncreaseCryst:      ',stepIncreaseCryst
  write(6,'(a24,x,i8)')   'nState:                 ',nState
  write(6,'(a24,x,i8)')   'nStress:                ',nStress
  write(6,'(a24,x,e8.1)') 'rTol_crystalliteState:  ',rTol_crystalliteState
  write(6,'(a24,x,e8.1)') 'rTol_crystalliteTemp:   ',rTol_crystalliteTemperature
  write(6,'(a24,x,e8.1)') 'rTol_crystalliteStress: ',rTol_crystalliteStress
  write(6,'(a24,x,e8.1)') 'aTol_crystalliteStress: ',aTol_crystalliteStress
  write(6,'(a24,2(x,i8))')'integrator:             ',numerics_integrator
  write(6,*)

  write(6,'(a24,x,i8)')   'nHomog:                 ',nHomog
  write(6,'(a24,x,e8.1)') 'subStepMinHomog:        ',subStepMinHomog
  write(6,'(a24,x,e8.1)') 'subStepSizeHomog:       ',subStepSizeHomog
  write(6,'(a24,x,e8.1)') 'stepIncreaseHomog:      ',stepIncreaseHomog
  write(6,'(a24,x,i8)')   'nMPstate:               ',nMPstate
  write(6,*)

!* RGC parameters
  write(6,'(a24,x,e8.1)') 'aTol_RGC:               ',absTol_RGC
  write(6,'(a24,x,e8.1)') 'rTol_RGC:               ',relTol_RGC
  write(6,'(a24,x,e8.1)') 'aMax_RGC:               ',absMax_RGC
  write(6,'(a24,x,e8.1)') 'rMax_RGC:               ',relMax_RGC
  write(6,'(a24,x,e8.1)') 'perturbPenalty_RGC:     ',pPert_RGC
  write(6,'(a24,x,e8.1)') 'relevantMismatch_RGC:   ',xSmoo_RGC
  write(6,'(a24,x,e8.1)') 'viscosityrate_RGC:      ',viscPower_RGC
  write(6,'(a24,x,e8.1)') 'viscositymodulus_RGC:   ',viscModus_RGC
  write(6,'(a24,x,e8.1)') 'maxrelaxation_RGC:      ',maxdRelax_RGC
  write(6,'(a24,x,e8.1)') 'maxVolDiscrepancy_RGC:  ',maxVolDiscr_RGC
  write(6,'(a24,x,e8.1)') 'volDiscrepancyMod_RGC:  ',volDiscrMod_RGC
  write(6,'(a24,x,e8.1)') 'discrepancyPower_RGC:   ',volDiscrPow_RGC
  write(6,*)

!* spectral parameters
  write(6,'(a24,x,e8.1)') 'err_div_tol:             ',err_div_tol
  write(6,'(a24,x,e8.1)') 'err_defgrad_tol:         ',err_defgrad_tol
  write(6,'(a24,x,e8.1)') 'err_stress_tolrel:       ',err_stress_tolrel
  write(6,'(a24,x,i8)')   'itmax:                   ',itmax
  write(6,'(a24,x,L8)')   'memory_efficient:          ',memory_efficient
  write(6,*)

!* Random seeding parameters
  write(6,'(a24,x,i8)')   'fixed_seed:             ',fixedSeed
  write(6,*)

!* openMP parameter
!$  write(6,'(a24,x,i8)')   'number of threads:      ',OMP_get_max_threads()
!$  write(6,*)
  
  ! sanity check  
  if (relevantStrain <= 0.0_pReal)          call IO_error(260)
  if (defgradTolerance <= 0.0_pReal)        call IO_error(294)
  if (iJacoStiffness < 1_pInt)              call IO_error(261)
  if (iJacoLpresiduum < 1_pInt)             call IO_error(262)
  if (pert_Fg <= 0.0_pReal)                 call IO_error(263)
  if (pert_method <= 0_pInt .or. pert_method >= 4_pInt) &
                                            call IO_error(299)
  if (nHomog < 1_pInt)                      call IO_error(264)
  if (nMPstate < 1_pInt)                    call IO_error(279)  !! missing in IO !!
  if (nCryst < 1_pInt)                      call IO_error(265)
  if (nState < 1_pInt)                      call IO_error(266)
  if (nStress < 1_pInt)                     call IO_error(267)
  if (subStepMinCryst <= 0.0_pReal)         call IO_error(268)
  if (subStepSizeCryst <= 0.0_pReal)        call IO_error(268)
  if (stepIncreaseCryst <= 0.0_pReal)       call IO_error(268)
  if (subStepMinHomog <= 0.0_pReal)         call IO_error(268)
  if (subStepSizeHomog <= 0.0_pReal)        call IO_error(268)
  if (stepIncreaseHomog <= 0.0_pReal)       call IO_error(268)
  if (rTol_crystalliteState <= 0.0_pReal)   call IO_error(269)
  if (rTol_crystalliteTemperature <= 0.0_pReal) call IO_error(276) !! oops !!
  if (rTol_crystalliteStress <= 0.0_pReal)  call IO_error(270)
  if (aTol_crystalliteStress <= 0.0_pReal)  call IO_error(271)
  if (any(numerics_integrator <= 0_pInt) .or. any(numerics_integrator >= 6_pInt)) &
                                            call IO_error(298)

!* RGC parameters: added <<<updated 17.11.2009>>>
  if (absTol_RGC <= 0.0_pReal)              call IO_error(272)
  if (relTol_RGC <= 0.0_pReal)              call IO_error(273)
  if (absMax_RGC <= 0.0_pReal)              call IO_error(274)
  if (relMax_RGC <= 0.0_pReal)              call IO_error(275)
  if (pPert_RGC <= 0.0_pReal)               call IO_error(276)   !! oops !!
  if (xSmoo_RGC <= 0.0_pReal)               call IO_error(277)
  if (viscPower_RGC < 0.0_pReal)            call IO_error(278)
  if (viscModus_RGC < 0.0_pReal)            call IO_error(278)
  if (refRelaxRate_RGC <= 0.0_pReal)        call IO_error(278)
  if (maxdRelax_RGC <= 0.0_pReal)           call IO_error(288)
  if (maxVolDiscr_RGC <= 0.0_pReal)         call IO_error(289)
  if (volDiscrMod_RGC < 0.0_pReal)          call IO_error(289)
  if (volDiscrPow_RGC <= 0.0_pReal)         call IO_error(289)

!* spectral parameters
  if (err_div_tol <= 0.0_pReal)             call IO_error(49)
  if (err_defgrad_tol <= 0.0_pReal)         call IO_error(49)
  if (err_stress_tolrel <= 0.0_pReal)       call IO_error(49)
  if (itmax <= 1.0_pInt)                    call IO_error(49)
  
  if (fixedSeed <= 0_pInt)                  write(6,'(a)') 'Random is random!'
endsubroutine

END MODULE numerics
