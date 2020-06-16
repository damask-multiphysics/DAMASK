!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Managing of parameters related to numerics
!--------------------------------------------------------------------------------------------------
module numerics
 use prec
 use IO
 use YAML_types
 use YAML_parse

#ifdef PETSc
#include <petsc/finclude/petscsys.h>
   use petscsys
#endif
!$ use OMP_LIB

 implicit none
 private

 class(tNode), pointer, public :: &
   numerics_root
 integer, protected, public :: &
   iJacoStiffness             =  1, &                                                               !< frequency of stiffness update
   randomSeed                 =  0, &                                                               !< fixed seeding for pseudo-random number generator, Default 0: use random seed
   worldrank                  =  0, &                                                               !< MPI worldrank (/=0 for MPI simulations only)
   worldsize                  =  1                                                                  !< MPI worldsize (/=1 for MPI simulations only)
 integer(4), protected, public :: &
   DAMASK_NumThreadsInt       =  0                                                                  !< value stored in environment variable DAMASK_NUM_THREADS, set to zero if no OpenMP directive
 real(pReal), protected, public :: &
   defgradTolerance           =  1.0e-7_pReal, &                                                    !< deviation of deformation gradient that is still allowed (used by CPFEM to determine outdated ffn1)
   numerics_unitlength        =  1.0_pReal, &                                                       !< determines the physical length of one computational length unit
   charLength                 =  1.0_pReal, &                                                       !< characteristic length scale for gradient problems
   residualStiffness          =  1.0e-6_pReal                                                       !< non-zero residual damage

!--------------------------------------------------------------------------------------------------
! field parameters:
 real(pReal), protected, public :: &
   err_struct_tolAbs          =  1.0e-10_pReal, &                                                   !< absolute tolerance for mechanical equilibrium
   err_struct_tolRel          =  1.0e-4_pReal                                                       !< relative tolerance for mechanical equilibrium
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
   polarAlpha                 =  1.0_pReal, &                                                       !< polarization scheme parameter 0.0 < alpha < 2.0. alpha = 1.0 ==> AL scheme, alpha = 2.0 ==> accelerated scheme
   polarBeta                  =  1.0_pReal                                                          !< polarization scheme parameter 0.0 < beta < 2.0. beta = 1.0 ==> AL scheme, beta = 2.0 ==> accelerated scheme
 character(len=pStringLen), protected, public :: &
   petsc_options              = ''
#endif

!--------------------------------------------------------------------------------------------------
! Mesh parameters:
#ifdef Mesh
 integer, protected, public :: &
   integrationOrder           =  2, &                                                              !< order of quadrature rule required
   structOrder                =  2                                                                 !< order of displacement shape functions
 logical, protected, public :: &
   BBarStabilisation          = .false.
 character(len=pStringLen), protected, public :: &
   petsc_options           = ''
#endif

 public :: numerics_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief reads in parameters from numerics.config and sets openMP related parameters. Also does
! a sanity check
!--------------------------------------------------------------------------------------------------
subroutine numerics_init
!$ integer ::                                gotDAMASK_NUM_THREADS = 1
 integer :: i,j, ierr
 character(len=:), allocatable :: &
   numerics_input, &
   numerics_inFlow, &
   key
 class (tNode), pointer :: &
   num_grid, &
   num_mesh, &
   num_generic
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

 inquire(file='numerics.yaml', exist=fexist)
 
 fileExists: if (fexist) then
   write(6,'(a,/)') ' using values from config file'
   flush(6)
   numerics_input =  IO_read('numerics.yaml')
   numerics_inFlow = to_flow(numerics_input)
   numerics_root =>  parse_flow(numerics_inFlow,defaultVal=emptyDict)
!--------------------------------------------------------------------------------------------------
! spectral parameters
   num_grid => numerics_root%get('grid',defaultVal=emptyDict)
   do i=1,num_grid%length
     key = num_grid%getKey(i)
     select case(key) 
#ifdef Grid
       case ('err_div_tolabs')
         err_div_tolAbs = num_grid%get_asFloat(key)
       case ('err_div_tolrel')
         err_div_tolRel = num_grid%get_asFloat(key)
       case ('err_stress_tolrel')
         err_stress_tolrel = num_grid%get_asFloat(key)
       case ('err_stress_tolabs')
         err_stress_tolabs = num_grid%get_asFloat(key)
       case ('err_curl_tolabs')
         err_curl_tolAbs = num_grid%get_asFloat(key)
       case ('err_curl_tolrel')
         err_curl_tolRel = num_grid%get_asFloat(key)
       case ('polaralpha')
         polarAlpha = num_grid%get_asFloat(key)
       case ('polarbeta')
         polarBeta = num_grid%get_asFloat(key)
#endif
       case ('itmax')
         itmax = num_grid%get_asInt(key)
       case ('itmin')
         itmin = num_grid%get_asInt(key)
       case ('maxCutBack')
         maxCutBack = num_grid%get_asInt(key)
       case ('maxStaggeredIter')
         stagItMax = num_grid%get_asInt(key)
#ifdef PETSC
       case ('petsc_options')
         petsc_options = num_grid%get_asString(key)
#endif 
     endselect
   enddo
   
   num_generic => numerics_root%get('generic',defaultVal=emptyDict)
   do i=1,num_generic%length
     key = num_generic%getKey(i)
     select case(key)
       case ('defgradtolerance')
         defgradTolerance = num_generic%get_asFloat(key)
       case ('ijacostiffness')
         iJacoStiffness = num_generic%get_asInt(key)
       case ('unitlength')
         numerics_unitlength = num_generic%get_asFloat(key)

!--------------------------------------------------------------------------------------------------
! random seeding parameter
       case ('fixed_seed', 'random_seed')
         randomSeed = num_generic%get_asInt(key)

!--------------------------------------------------------------------------------------------------
! gradient parameter
       case ('charLength')
         charLength = num_generic%get_asFloat(key)
       case ('residualStiffness')
         residualStiffness = num_generic%get_asFloat(key)
!--------------------------------------------------------------------------------------------------
! field parameters
       case ('err_struct_tolabs')
         err_struct_tolAbs = num_generic%get_asFloat(key)
       case ('err_struct_tolrel')
         err_struct_tolRel = num_generic%get_asFloat(key)
     endselect
   enddo

#ifdef Mesh
   num_grid => numerics_root%get('mesh',defaultVal=emptyDict)
   do i=1,num_grid%length
     key = num_grid%getKey(i)
     select case(key) 
       case ('integrationorder')
         integrationorder = num_generic%get_asInt(key)
       case ('structorder')
         structorder = num_generic%get_asInt(key)
       case ('bbarstabilisation')
         BBarStabilisation = num_generic%get_asInt(key) > 0
     end select
   enddo
#endif

 else fileExists
   write(6,'(a,/)') ' using standard values'
   flush(6)
 endif fileExists

!--------------------------------------------------------------------------------------------------
! writing parameters to output
 write(6,'(a24,1x,es8.1)')  ' defgradTolerance:       ',defgradTolerance
 write(6,'(a24,1x,i8)')     ' iJacoStiffness:         ',iJacoStiffness
 write(6,'(a24,1x,es8.1,/)')' unitlength:             ',numerics_unitlength

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

!--------------------------------------------------------------------------------------------------
! spectral parameters
#ifdef Grid
 write(6,'(a24,1x,es8.1)')   ' err_stress_tolAbs:      ',err_stress_tolAbs
 write(6,'(a24,1x,es8.1)')   ' err_stress_tolRel:      ',err_stress_tolRel
 write(6,'(a24,1x,es8.1)')   ' err_div_tolAbs:         ',err_div_tolAbs
 write(6,'(a24,1x,es8.1)')   ' err_div_tolRel:         ',err_div_tolRel
 write(6,'(a24,1x,es8.1)')   ' err_curl_tolAbs:        ',err_curl_tolAbs
 write(6,'(a24,1x,es8.1)')   ' err_curl_tolRel:        ',err_curl_tolRel
 write(6,'(a24,1x,es8.1)')   ' polarAlpha:             ',polarAlpha
 write(6,'(a24,1x,es8.1)')   ' polarBeta:              ',polarBeta
#endif

!--------------------------------------------------------------------------------------------------
! spectral parameters
#ifdef Mesh
 write(6,'(a24,1x,i8)')      ' integrationOrder:       ',integrationOrder
 write(6,'(a24,1x,i8)')      ' structOrder:            ',structOrder
 write(6,'(a24,1x,L8)')      ' B-Bar stabilisation:    ',BBarStabilisation
#endif

#ifdef PETSC
 write(6,'(a24,1x,a)')       ' PETSc_options:          ',trim(petsc_options)
#endif

!--------------------------------------------------------------------------------------------------
! sanity checks
 if (defgradTolerance <= 0.0_pReal)        call IO_error(301,ext_msg='defgradTolerance')
 if (iJacoStiffness < 1)                   call IO_error(301,ext_msg='iJacoStiffness')
 if (numerics_unitlength <= 0.0_pReal)     call IO_error(301,ext_msg='unitlength')
 if (residualStiffness < 0.0_pReal)        call IO_error(301,ext_msg='residualStiffness')
 if (itmax <= 1)                           call IO_error(301,ext_msg='itmax')
 if (itmin > itmax .or. itmin < 1)         call IO_error(301,ext_msg='itmin')
 if (maxCutBack < 0)                       call IO_error(301,ext_msg='maxCutBack')
 if (stagItMax < 0)                        call IO_error(301,ext_msg='maxStaggeredIter')
 if (err_struct_tolRel <= 0.0_pReal)       call IO_error(301,ext_msg='err_struct_tolRel')
 if (err_struct_tolAbs <= 0.0_pReal)       call IO_error(301,ext_msg='err_struct_tolAbs')
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
