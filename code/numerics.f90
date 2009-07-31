!##############################################################
MODULE numerics
!##############################################################

use prec, only: pInt, pReal
implicit none

character(len=64), parameter :: numerics_configFile = 'numerics.config' ! name of configuration file
integer(pInt)                   iJacoStiffness, &                       ! frequency of stiffness update
                                iJacoLpresiduum, &                      ! frequency of Jacobian update of residuum in Lp
                                nHomog, &                               ! homogenization loop limit
                                nCryst, &                               ! crystallite loop limit (only for debugging info, real loop limit is "subStepMin")
                                nState, &                               ! state loop limit
                                nStress                                 ! stress loop limit
real(pReal)                     relevantStrain, &                       ! strain increment considered significant
                                pert_Fg, &                              ! strain perturbation for FEM Jacobi
                                subStepMin, &                           ! minimum (relative) size of sub-step allowed during cutback in crystallite
                                rTol_crystalliteState, &                ! relative tolerance in crystallite state loop 
                                rTol_crystalliteTemperature, &          ! relative tolerance in crystallite temperature loop 
                                rTol_crystalliteStress, &               ! relative tolerance in crystallite stress loop
                                aTol_crystalliteStress, &               ! absolute tolerance in crystallite stress loop

!* RGC parameters: added <<<updated 31.07.2009>>>
                                absTol_RGC, &                           ! absolute tolerance of RGC residuum
                                relTol_RGC, &                           ! relative tolerance of RGC residuum
                                absMax_RGC, &                           ! absolute maximum of RGC residuum
                                relMax_RGC, &                           ! relative maximum of RGC residuum
                                pPert_RGC, &                            ! perturbation for computing RGC penalty tangent
                                xSmoo_RGC                               ! RGC penalty smoothing parameter (hyperbolic tangent)

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
  
  implicit none

  !*** input variables ***!
  
  !*** output variables ***!
  
  !*** local variables ***!
  integer(pInt), parameter ::                 fileunit = 300  
  integer(pInt), parameter ::                 maxNchunks = 2
  integer(pInt), dimension(1+2*maxNchunks) :: positions
  character(len=64)                           tag
  character(len=1024)                         line
  
  write(6,*)
  write(6,*) '<<<+-  numerics init  -+>>>'
  write(6,*)
  
  ! initialize all parameters with standard values
  relevantStrain          = 1.0e-7_pReal
  iJacoStiffness          = 1_pInt
  iJacoLpresiduum         = 1_pInt
  pert_Fg                 = 1.0e-6_pReal
  nHomog                  = 10_pInt
  nCryst                  = 20_pInt
  nState                  = 10_pInt
  nStress                 = 40_pInt
  subStepMin              = 1.0e-3_pReal
  rTol_crystalliteState   = 1.0e-6_pReal
  rTol_crystalliteTemperature = 1.0e-6_pReal
  rTol_crystalliteStress  = 1.0e-6_pReal
  aTol_crystalliteStress  = 1.0e-8_pReal

!* RGC parameters: added <<<updated 31.07.2009>>>
  absTol_RGC              = 1.0e+3
  relTol_RGC              = 1.0e-3
  absMax_RGC              = 1.0e+9
  relMax_RGC              = 1.0e+2
  pPert_RGC               = 1.0e-8
  xSmoo_RGC               = 1.0e-5

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
        case ('ijacostiffness')
              iJacoStiffness = IO_intValue(line,positions,2)
        case ('ijacolpresiduum')
              iJacoLpresiduum = IO_intValue(line,positions,2)
        case ('pert_fg')
              pert_Fg = IO_floatValue(line,positions,2)
        case ('nhomog')
              nHomog = IO_intValue(line,positions,2)
        case ('ncryst')
              nCryst = IO_intValue(line,positions,2)
        case ('nstate')
              nState = IO_intValue(line,positions,2)
        case ('nstress')
              nStress = IO_intValue(line,positions,2)
        case ('substepmin')
              subStepMin = IO_floatValue(line,positions,2)
        case ('rtol_crystallitestate')
              rTol_crystalliteState = IO_floatValue(line,positions,2)
        case ('rtol_crystallitetemperature')
              rTol_crystalliteTemperature = IO_floatValue(line,positions,2)
        case ('rtol_crystallitestress')
              rTol_crystalliteStress = IO_floatValue(line,positions,2)
        case ('atol_crystallitestress')
              aTol_crystalliteStress = IO_floatValue(line,positions,2)

!* RGC parameters: added <<<updated 31.07.2009>>>
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
  write(6,'(a24,x,i8)')   'iJacoStiffness:         ',iJacoStiffness
  write(6,'(a24,x,i8)')   'iJacoLpresiduum:        ',iJacoLpresiduum
  write(6,'(a24,x,e8.1)') 'pert_Fg:                ',pert_Fg
  write(6,'(a24,x,i8)')   'nHomog:                 ',nHomog
  write(6,'(a24,x,i8)')   'nCryst:                 ',nCryst
  write(6,'(a24,x,i8)')   'nState:                 ',nState
  write(6,'(a24,x,i8)')   'nStress:                ',nStress
  write(6,'(a24,x,e8.1)') 'subStepMin:             ',subStepMin
  write(6,'(a24,x,e8.1)') 'rTol_crystalliteState:  ',rTol_crystalliteState
  write(6,'(a24,x,e8.1)') 'rTol_crystalliteTemp:   ',rTol_crystalliteTemperature
  write(6,'(a24,x,e8.1)') 'rTol_crystalliteStress: ',rTol_crystalliteStress
  write(6,'(a24,x,e8.1)') 'aTol_crystalliteStress: ',aTol_crystalliteStress

!* RGC parameters: added <<<updated 31.07.2009>>>
  write(6,'(a24,x,e8.1)') 'aTol_RGC:             ',absTol_RGC
  write(6,'(a24,x,e8.1)') 'rTol_RGC:             ',relTol_RGC
  write(6,'(a24,x,e8.1)') 'aMax_RGC:             ',absMax_RGC
  write(6,'(a24,x,e8.1)') 'rMax_RGC:             ',relMax_RGC
  write(6,'(a24,x,e8.1)') 'perturbPenalty_RGC:   ',pPert_RGC
  write(6,'(a24,x,e8.1)') 'relevantMismatch_RGC: ',xSmoo_RGC
  write(6,*)
  
  ! sanity check  
  if (relevantStrain <= 0.0_pReal)          call IO_error(260)
  if (iJacoStiffness < 1_pInt)              call IO_error(261)
  if (iJacoLpresiduum < 1_pInt)             call IO_error(262)
  if (pert_Fg <= 0.0_pReal)                 call IO_error(263)
  if (nHomog < 1_pInt)                      call IO_error(264)
  if (nCryst < 1_pInt)                      call IO_error(265)
  if (nState < 1_pInt)                      call IO_error(266)
  if (nStress < 1_pInt)                     call IO_error(267)
  if (subStepMin <= 0.0_pReal)              call IO_error(268)
  if (rTol_crystalliteState <= 0.0_pReal)   call IO_error(269)
  if (rTol_crystalliteTemperature <= 0.0_pReal) call IO_error(276)
  if (rTol_crystalliteStress <= 0.0_pReal)  call IO_error(270)
  if (aTol_crystalliteStress <= 0.0_pReal)  call IO_error(271)

!* RGC parameters: added <<<updated 31.07.2009>>>
  if (absTol_RGC <= 0.0_pReal)              call IO_error(272)
  if (relTol_RGC <= 0.0_pReal)              call IO_error(273)
  if (absMax_RGC <= 0.0_pReal)              call IO_error(274)
  if (relMax_RGC <= 0.0_pReal)              call IO_error(275)
  if (pPert_RGC <= 0.0_pReal)               call IO_error(276)
  if (xSmoo_RGC <= 0.0_pReal)               call IO_error(277)
 
endsubroutine

END MODULE numerics
