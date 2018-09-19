!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Chen Zhang, Michigan State University
!> @brief crystallite state integration functions and reporting of results
!--------------------------------------------------------------------------------------------------

module crystallite
 use prec, only: &
   pReal, &
   pInt

 implicit none

 private
 character(len=64),         dimension(:,:),          allocatable, private :: &
   crystallite_output                                                                               !< name of each post result output
 integer(pInt),                                                   public, protected :: &
   crystallite_maxSizePostResults                                                                   !< description not available
 integer(pInt),             dimension(:),            allocatable, public, protected :: &
   crystallite_sizePostResults                                                                      !< description not available
 integer(pInt),             dimension(:,:),          allocatable, private :: &
   crystallite_sizePostResult                                                                       !< description not available

 real(pReal),               dimension(:,:,:),        allocatable, public :: &
   crystallite_dt                                                                                   !< requested time increment of each grain
 real(pReal),               dimension(:,:,:),        allocatable, private :: &
   crystallite_subdt, &                                                                             !< substepped time increment of each grain
   crystallite_subFrac, &                                                                           !< already calculated fraction of increment
   crystallite_subStep                                                                              !< size of next integration step
 real(pReal),               dimension(:,:,:,:),      allocatable, public :: &
   crystallite_Tstar_v, &                                                                           !< current 2nd Piola-Kirchhoff stress vector (end of converged time step)
   crystallite_Tstar0_v, &                                                                          !< 2nd Piola-Kirchhoff stress vector at start of FE inc
   crystallite_partionedTstar0_v                                                                    !< 2nd Piola-Kirchhoff stress vector at start of homog inc
 real(pReal),               dimension(:,:,:,:),      allocatable, private :: &
   crystallite_subTstar0_v, &                                                                       !< 2nd Piola-Kirchhoff stress vector at start of crystallite inc
   crystallite_orientation, &                                                                       !< orientation as quaternion
   crystallite_orientation0, &                                                                      !< initial orientation as quaternion
   crystallite_rotation                                                                             !< grain rotation away from initial orientation as axis-angle (in degrees) in crystal reference frame
 real(pReal),               dimension(:,:,:,:,:),    allocatable, public :: &
   crystallite_Fp, &                                                                                !< current plastic def grad (end of converged time step)
   crystallite_Fp0, &                                                                               !< plastic def grad at start of FE inc
   crystallite_partionedFp0,&                                                                       !< plastic def grad at start of homog inc
   crystallite_Fi, &                                                                                !< current intermediate def grad (end of converged time step)
   crystallite_Fi0, &                                                                               !< intermediate def grad at start of FE inc
   crystallite_partionedFi0,&                                                                       !< intermediate def grad at start of homog inc
   crystallite_F0, &                                                                                !< def grad at start of FE inc
   crystallite_partionedF,  &                                                                       !< def grad to be reached at end of homog inc
   crystallite_partionedF0, &                                                                       !< def grad at start of homog inc
   crystallite_Lp, &                                                                                !< current plastic velocitiy grad (end of converged time step)
   crystallite_Lp0, &                                                                               !< plastic velocitiy grad at start of FE inc
   crystallite_partionedLp0,&                                                                       !< plastic velocity grad at start of homog inc
   crystallite_Li, &                                                                                !< current intermediate velocitiy grad (end of converged time step)
   crystallite_Li0, &                                                                               !< intermediate velocitiy grad at start of FE inc
   crystallite_partionedLi0,&                                                                       !< intermediate velocity grad at start of homog inc
   crystallite_Fe, &                                                                                !< current "elastic" def grad (end of converged time step)
   crystallite_P                                                                                    !< 1st Piola-Kirchhoff stress per grain
 real(pReal),                dimension(:,:,:,:,:),    allocatable, private :: &
   crystallite_invFp, &                                                                             !< inverse of current plastic def grad (end of converged time step)
   crystallite_subFp0,&                                                                             !< plastic def grad at start of crystallite inc
   crystallite_invFi, &                                                                             !< inverse of current intermediate def grad (end of converged time step)
   crystallite_subFi0,&                                                                             !< intermediate def grad at start of crystallite inc
   crystallite_subF,  &                                                                             !< def grad to be reached at end of crystallite inc
   crystallite_subF0, &                                                                             !< def grad at start of crystallite inc
   crystallite_subLp0,&                                                                             !< plastic velocity grad at start of crystallite inc
   crystallite_subLi0,&                                                                             !< intermediate velocity grad at start of crystallite inc
   crystallite_disorientation                                                                       !< disorientation between two neighboring ips (only calculated for single grain IPs)
 real(pReal),                dimension(:,:,:,:,:,:,:), allocatable, public :: &
   crystallite_dPdF, &                                                                              !< current individual dPdF per grain (end of converged time step)
   crystallite_dPdF0, &                                                                             !< individual dPdF per grain at start of FE inc
   crystallite_partioneddPdF0                                                                       !< individual dPdF per grain at start of homog inc
 logical,                    dimension(:,:,:),         allocatable, public :: &
   crystallite_requested                                                                            !< flag to request crystallite calculation
 logical,                    dimension(:,:,:),         allocatable, public, protected :: &
   crystallite_converged, &                                                                         !< convergence flag
   crystallite_localPlasticity                                                                      !< indicates this grain to have purely local constitutive law
 logical,                    dimension(:,:,:),         allocatable, private :: &
   crystallite_todo                                                                                 !< flag to indicate need for further computation
 logical,                    dimension(:,:),           allocatable, private :: &
   crystallite_clearToWindForward, &                                                                !< description not available
   crystallite_clearToCutback, &                                                                    !< description not available
   crystallite_syncSubFrac, &                                                                       !< description not available
   crystallite_syncSubFracCompleted, &                                                              !< description not available
   crystallite_neighborEnforcedCutback                                                              !< description not available

 enum, bind(c)
   enumerator :: undefined_ID, &
                 phase_ID, &
                 texture_ID, &
                 volume_ID, &
                 grainrotationx_ID, &
                 grainrotationy_ID, &
                 grainrotationz_ID, &
                 orientation_ID, &
                 grainrotation_ID, &
                 eulerangles_ID, &
                 defgrad_ID, &
                 fe_ID, &
                 fp_ID, &
                 fi_ID, &
                 lp_ID, &
                 li_ID, &
                 e_ID, &
                 ee_ID, &
                 p_ID, &
                 s_ID, &
                 elasmatrix_ID, &
                 neighboringip_ID, &
                 neighboringelement_ID
 end enum
 integer(kind(undefined_ID)),dimension(:,:),   allocatable,          private :: &
   crystallite_outputID                                                                             !< ID of each post result output
 procedure(), pointer :: integrateState

 public :: &
   crystallite_init, &
   crystallite_stressAndItsTangent, &
   crystallite_orientations, &
   crystallite_push33ToRef, &
   crystallite_postResults
 private :: &
   integrateState, &
   crystallite_integrateStateFPI, &
   crystallite_integrateStateEuler, &
   crystallite_integrateStateAdaptiveEuler, &
   crystallite_integrateStateRK4, &
   crystallite_integrateStateRKCK45, &
   crystallite_integrateStress, &
   crystallite_stateJump

contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates and initialize per grain variables
!--------------------------------------------------------------------------------------------------
subroutine crystallite_init
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use debug, only: &
   debug_info, &
   debug_reset, &
   debug_level, &
   debug_crystallite, &
   debug_levelBasic
 use numerics, only: &
   numerics_integrator, &
   worldrank, &
   usePingPong
 use math, only: &
   math_I3, &
   math_EulerToR, &
   math_inv33, &
   math_mul33xx33, &
   math_mul33x33
 use FEsolving, only:  &
   FEsolving_execElem, &
   FEsolving_execIP
 use mesh, only: &
   mesh_element, &
   mesh_NcpElems, &
   mesh_maxNips, &
   mesh_maxNipNeighbors
 use IO, only: &
   IO_timeStamp, &
   IO_stringValue, &
   IO_write_jobFile, &
   IO_error
 use material
 use config, only: &
  config_deallocate, &
  config_crystallite, &
  crystallite_name, &
  material_Nphase
 use constitutive, only: &
   constitutive_initialFi, &
   constitutive_microstructure                                                                      ! derived (shortcut) quantities of given state

 implicit none

 integer(pInt), parameter :: FILEUNIT=434_pInt
 integer(pInt) :: &
   c, &                                                                                             !< counter in integration point component loop
   i, &                                                                                             !< counter in integration point loop
   e, &                                                                                             !< counter in element loop
   o = 0_pInt, &                                                                                    !< counter in output loop
   r, &  
   ph, &                                                                                            !< counter in crystallite loop
   cMax, &                                                                                          !< maximum number of  integration point components
   iMax, &                                                                                          !< maximum number of integration points
   eMax, &                                                                                          !< maximum number of elements
   nMax, &                                                                                          !< maximum number of ip neighbors
   myNcomponents, &                                                                                 !< number of components at current IP
   mySize

 character(len=65536), dimension(:), allocatable :: str
 character(len=65536) :: &
   tag = ''

 write(6,'(/,a)')   ' <<<+-  crystallite init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 cMax = homogenization_maxNgrains
 iMax = mesh_maxNips
 eMax = mesh_NcpElems
 nMax = mesh_maxNipNeighbors


 allocate(crystallite_Tstar0_v(6,cMax,iMax,eMax),            source=0.0_pReal)
 allocate(crystallite_partionedTstar0_v(6,cMax,iMax,eMax),   source=0.0_pReal)
 allocate(crystallite_subTstar0_v(6,cMax,iMax,eMax),         source=0.0_pReal)
 allocate(crystallite_Tstar_v(6,cMax,iMax,eMax),             source=0.0_pReal)
 allocate(crystallite_P(3,3,cMax,iMax,eMax),                 source=0.0_pReal)
 allocate(crystallite_F0(3,3,cMax,iMax,eMax),                source=0.0_pReal)
 allocate(crystallite_partionedF0(3,3,cMax,iMax,eMax),       source=0.0_pReal)
 allocate(crystallite_partionedF(3,3,cMax,iMax,eMax),        source=0.0_pReal)
 allocate(crystallite_subF0(3,3,cMax,iMax,eMax),             source=0.0_pReal)
 allocate(crystallite_subF(3,3,cMax,iMax,eMax),              source=0.0_pReal)
 allocate(crystallite_Fp0(3,3,cMax,iMax,eMax),               source=0.0_pReal)
 allocate(crystallite_partionedFp0(3,3,cMax,iMax,eMax),      source=0.0_pReal)
 allocate(crystallite_subFp0(3,3,cMax,iMax,eMax),            source=0.0_pReal)
 allocate(crystallite_Fp(3,3,cMax,iMax,eMax),                source=0.0_pReal)
 allocate(crystallite_invFp(3,3,cMax,iMax,eMax),             source=0.0_pReal)
 allocate(crystallite_Fi0(3,3,cMax,iMax,eMax),               source=0.0_pReal)
 allocate(crystallite_partionedFi0(3,3,cMax,iMax,eMax),      source=0.0_pReal)
 allocate(crystallite_subFi0(3,3,cMax,iMax,eMax),            source=0.0_pReal)
 allocate(crystallite_Fi(3,3,cMax,iMax,eMax),                source=0.0_pReal)
 allocate(crystallite_invFi(3,3,cMax,iMax,eMax),             source=0.0_pReal)
 allocate(crystallite_Fe(3,3,cMax,iMax,eMax),                source=0.0_pReal)
 allocate(crystallite_Lp0(3,3,cMax,iMax,eMax),               source=0.0_pReal)
 allocate(crystallite_partionedLp0(3,3,cMax,iMax,eMax),      source=0.0_pReal)
 allocate(crystallite_subLp0(3,3,cMax,iMax,eMax),            source=0.0_pReal)
 allocate(crystallite_Lp(3,3,cMax,iMax,eMax),                source=0.0_pReal)
 allocate(crystallite_Li0(3,3,cMax,iMax,eMax),               source=0.0_pReal)
 allocate(crystallite_partionedLi0(3,3,cMax,iMax,eMax),      source=0.0_pReal)
 allocate(crystallite_subLi0(3,3,cMax,iMax,eMax),            source=0.0_pReal)
 allocate(crystallite_Li(3,3,cMax,iMax,eMax),                source=0.0_pReal)
 allocate(crystallite_dPdF(3,3,3,3,cMax,iMax,eMax),          source=0.0_pReal)
 allocate(crystallite_dPdF0(3,3,3,3,cMax,iMax,eMax),         source=0.0_pReal)
 allocate(crystallite_partioneddPdF0(3,3,3,3,cMax,iMax,eMax),source=0.0_pReal)
 allocate(crystallite_dt(cMax,iMax,eMax),                    source=0.0_pReal)
 allocate(crystallite_subdt(cMax,iMax,eMax),                 source=0.0_pReal)
 allocate(crystallite_subFrac(cMax,iMax,eMax),               source=0.0_pReal)
 allocate(crystallite_subStep(cMax,iMax,eMax),               source=0.0_pReal)
 allocate(crystallite_orientation(4,cMax,iMax,eMax),         source=0.0_pReal)
 allocate(crystallite_orientation0(4,cMax,iMax,eMax),        source=0.0_pReal)
 allocate(crystallite_rotation(4,cMax,iMax,eMax),            source=0.0_pReal)
 if (any(plasticState%nonLocal)) &
   allocate(crystallite_disorientation(4,nMax,cMax,iMax,eMax),source=0.0_pReal)
 allocate(crystallite_localPlasticity(cMax,iMax,eMax),       source=.true.)
 allocate(crystallite_requested(cMax,iMax,eMax),             source=.false.)
 allocate(crystallite_todo(cMax,iMax,eMax),                  source=.false.)
 allocate(crystallite_converged(cMax,iMax,eMax),             source=.true.)
 allocate(crystallite_clearToWindForward(iMax,eMax),         source=.true.)
 allocate(crystallite_syncSubFrac(iMax,eMax),                source=.false.)
 allocate(crystallite_syncSubFracCompleted(iMax,eMax),       source=.false.)
 allocate(crystallite_clearToCutback(iMax,eMax),             source=.true.)
 allocate(crystallite_neighborEnforcedCutback(iMax,eMax),    source=.false.)
 allocate(crystallite_output(maxval(crystallite_Noutput), &
                             size(config_crystallite))) ;       crystallite_output = ''
 allocate(crystallite_outputID(maxval(crystallite_Noutput), &
                             size(config_crystallite)),         source=undefined_ID)
 allocate(crystallite_sizePostResults(size(config_crystallite)),source=0_pInt)
 allocate(crystallite_sizePostResult(maxval(crystallite_Noutput), &
                                     size(config_crystallite)), source=0_pInt)

 select case(numerics_integrator(1))
   case(1_pInt)
     integrateState => crystallite_integrateStateFPI
   case(2_pInt)
     integrateState => crystallite_integrateStateEuler
   case(3_pInt)
     integrateState => crystallite_integrateStateAdaptiveEuler
   case(4_pInt)
     integrateState => crystallite_integrateStateRK4
   case(5_pInt)
     integrateState => crystallite_integrateStateRKCK45
 end select



 do c = 1_pInt, size(config_crystallite)
#if defined(__GFORTRAN__)
   str = ['GfortranBug86277']
   str = config_crystallite(c)%getStrings('(output)',defaultVal=str)
   if (str(1) == 'GfortranBug86277') str = [character(len=65536)::]
#else
   str = config_crystallite(c)%getStrings('(output)',defaultVal=[character(len=65536)::])
#endif
   do o = 1_pInt, size(str)
     crystallite_output(o,c) = str(o)
     outputName: select case(str(o))
           case ('phase') outputName
             crystallite_outputID(o,c) = phase_ID
           case ('texture') outputName
             crystallite_outputID(o,c) = texture_ID
           case ('volume') outputName
             crystallite_outputID(o,c) = volume_ID
           case ('grainrotationx') outputName
             crystallite_outputID(o,c) = grainrotationx_ID
           case ('grainrotationy') outputName
             crystallite_outputID(o,c) = grainrotationy_ID
           case ('grainrotationz') outputName
             crystallite_outputID(o,c) = grainrotationx_ID
           case ('orientation') outputName
             crystallite_outputID(o,c) = orientation_ID
           case ('grainrotation') outputName
             crystallite_outputID(o,c) = grainrotation_ID
           case ('eulerangles') outputName
             crystallite_outputID(o,c) = eulerangles_ID
           case ('defgrad','f') outputName
             crystallite_outputID(o,c) = defgrad_ID
           case ('fe') outputName
             crystallite_outputID(o,c) = fe_ID
           case ('fp') outputName
             crystallite_outputID(o,c) = fp_ID
           case ('fi') outputName
             crystallite_outputID(o,c) = fi_ID
           case ('lp') outputName
             crystallite_outputID(o,c) = lp_ID
           case ('li') outputName
             crystallite_outputID(o,c) = li_ID
           case ('e') outputName
             crystallite_outputID(o,c) = e_ID
           case ('ee') outputName
             crystallite_outputID(o,c) = ee_ID
           case ('p','firstpiola','1stpiola') outputName
             crystallite_outputID(o,c) = p_ID
           case ('s','tstar','secondpiola','2ndpiola') outputName
             crystallite_outputID(o,c) = s_ID
           case ('elasmatrix') outputName
             crystallite_outputID(o,c) = elasmatrix_ID
           case ('neighboringip') outputName
             crystallite_outputID(o,c) = neighboringip_ID
           case ('neighboringelement') outputName
             crystallite_outputID(o,c) = neighboringelement_ID
           case default outputName
             call IO_error(105_pInt,ext_msg=tag//' (Crystallite)')
         end select outputName
   enddo
 enddo


 do r = 1_pInt,size(config_crystallite)
   do o = 1_pInt,crystallite_Noutput(r)
     select case(crystallite_outputID(o,r))
       case(phase_ID,texture_ID,volume_ID,grainrotationx_ID,grainrotationy_ID,grainrotationz_ID)
         mySize = 1_pInt
       case(orientation_ID,grainrotation_ID)
         mySize = 4_pInt
       case(eulerangles_ID)
         mySize = 3_pInt
       case(defgrad_ID,fe_ID,fp_ID,fi_ID,lp_ID,li_ID,e_ID,ee_ID,p_ID,s_ID)
         mySize = 9_pInt
       case(elasmatrix_ID)
         mySize = 36_pInt
       case(neighboringip_ID,neighboringelement_ID)
         mySize = mesh_maxNipNeighbors
       case default
         mySize = 0_pInt
     end select
     crystallite_sizePostResult(o,r) = mySize
     crystallite_sizePostResults(r) = crystallite_sizePostResults(r) + mySize
   enddo
 enddo

 crystallite_maxSizePostResults = &
   maxval(crystallite_sizePostResults(microstructure_crystallite),microstructure_active)
 

!--------------------------------------------------------------------------------------------------
! write description file for crystallite output
 if (worldrank == 0_pInt) then
   call IO_write_jobFile(FILEUNIT,'outputCrystallite')

   do r = 1_pInt,size(config_crystallite)
     if (any(microstructure_crystallite(mesh_element(4,:)) == r)) then
       write(FILEUNIT,'(/,a,/)') '['//trim(crystallite_name(r))//']'
       do o = 1_pInt,crystallite_Noutput(r)
         write(FILEUNIT,'(a,i4)') trim(crystallite_output(o,r))//char(9),crystallite_sizePostResult(o,r)
       enddo
     endif
   enddo

   close(FILEUNIT)
 endif

 call config_deallocate('material.config/crystallite')

!--------------------------------------------------------------------------------------------------
! initialize
!$OMP PARALLEL DO PRIVATE(myNcomponents)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     myNcomponents = homogenization_Ngrains(mesh_element(3,e))
     forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), c = 1_pInt:myNcomponents)
       crystallite_Fp0(1:3,1:3,c,i,e) = math_EulerToR(material_EulerAngles(1:3,c,i,e))              ! plastic def gradient reflects init orientation
       crystallite_Fi0(1:3,1:3,c,i,e) = constitutive_initialFi(c,i,e)
       crystallite_F0(1:3,1:3,c,i,e)  = math_I3
       crystallite_localPlasticity(c,i,e) = phase_localPlasticity(material_phase(c,i,e))
       crystallite_Fe(1:3,1:3,c,i,e)  = math_inv33(math_mul33x33(crystallite_Fi0(1:3,1:3,c,i,e), &
                                                                 crystallite_Fp0(1:3,1:3,c,i,e)))   ! assuming that euler angles are given in internal strain free configuration
       crystallite_Fp(1:3,1:3,c,i,e)  = crystallite_Fp0(1:3,1:3,c,i,e)
       crystallite_Fi(1:3,1:3,c,i,e)  = crystallite_Fi0(1:3,1:3,c,i,e)
       crystallite_requested(c,i,e) = .true.
     endforall
   enddo
 !$OMP END PARALLEL DO

 if(any(.not. crystallite_localPlasticity) .and. .not. usePingPong) call IO_error(601_pInt)         ! exit if nonlocal but no ping-pong

 crystallite_partionedFp0 = crystallite_Fp0
 crystallite_partionedFi0 = crystallite_Fi0
 crystallite_partionedF0  = crystallite_F0
 crystallite_partionedF   = crystallite_F0

 call crystallite_orientations()
 crystallite_orientation0 = crystallite_orientation                                                 ! store initial orientations for calculation of grain rotations

 !$OMP PARALLEL DO PRIVATE(myNcomponents)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     myNcomponents = homogenization_Ngrains(mesh_element(3,e))
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do c = 1_pInt,myNcomponents
         call constitutive_microstructure(crystallite_orientation, &                                ! pass orientation to constitutive module
                                          crystallite_Fe(1:3,1:3,c,i,e), &
                                          crystallite_Fp(1:3,1:3,c,i,e), &
                                          c,i,e)                                                    ! update dependent state variables to be consistent with basic states
      enddo
     enddo
   enddo
 !$OMP END PARALLEL DO

 call crystallite_stressAndItsTangent(.true.)                                                       ! request elastic answers

!--------------------------------------------------------------------------------------------------
! debug output
 if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Fe:                    ', shape(crystallite_Fe)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Fp:                    ', shape(crystallite_Fp)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Fi:                    ', shape(crystallite_Fi)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Lp:                    ', shape(crystallite_Lp)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Li:                    ', shape(crystallite_Li)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_F0:                    ', shape(crystallite_F0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Fp0:                   ', shape(crystallite_Fp0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Fi0:                   ', shape(crystallite_Fi0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Lp0:                   ', shape(crystallite_Lp0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Li0:                   ', shape(crystallite_Li0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedF:            ', shape(crystallite_partionedF)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedF0:           ', shape(crystallite_partionedF0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedFp0:          ', shape(crystallite_partionedFp0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedFi0:          ', shape(crystallite_partionedFi0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedLp0:          ', shape(crystallite_partionedLp0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedLi0:          ', shape(crystallite_partionedLi0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subF:                  ', shape(crystallite_subF)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subF0:                 ', shape(crystallite_subF0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subFp0:                ', shape(crystallite_subFp0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subFi0:                ', shape(crystallite_subFi0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subLp0:                ', shape(crystallite_subLp0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subLi0:                ', shape(crystallite_subLi0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_P:                     ', shape(crystallite_P)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Tstar_v:               ', shape(crystallite_Tstar_v)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Tstar0_v:              ', shape(crystallite_Tstar0_v)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedTstar0_v:     ', shape(crystallite_partionedTstar0_v)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subTstar0_v:           ', shape(crystallite_subTstar0_v)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_dPdF:                  ', shape(crystallite_dPdF)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_dPdF0:                 ', shape(crystallite_dPdF0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_partioneddPdF0:        ', shape(crystallite_partioneddPdF0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_orientation:           ', shape(crystallite_orientation)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_orientation0:          ', shape(crystallite_orientation0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_rotation:              ', shape(crystallite_rotation)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_disorientation:        ', shape(crystallite_disorientation)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_dt:                    ', shape(crystallite_dt)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subdt:                 ', shape(crystallite_subdt)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subFrac:               ', shape(crystallite_subFrac)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subStep:               ', shape(crystallite_subStep)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_localPlasticity:       ', shape(crystallite_localPlasticity)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_requested:             ', shape(crystallite_requested)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_todo:                  ', shape(crystallite_todo)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_converged:             ', shape(crystallite_converged)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_sizePostResults:       ', shape(crystallite_sizePostResults)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_sizePostResult:        ', shape(crystallite_sizePostResult)
   write(6,'(/,a35,1x,i10)')    'Number of nonlocal grains:         ',count(.not. crystallite_localPlasticity)
   flush(6)
 endif

 call debug_info
 call debug_reset

end subroutine crystallite_init


!--------------------------------------------------------------------------------------------------
!> @brief calculate stress (P) and tangent (dPdF) for crystallites
!--------------------------------------------------------------------------------------------------
subroutine crystallite_stressAndItsTangent(updateJaco)
 use prec, only: &
   tol_math_check, &
   dNeq0
 use numerics, only: &
   subStepMinCryst, &
   subStepSizeCryst, &
   stepIncreaseCryst, &
   numerics_timeSyncing
 use debug, only: &
   debug_level, &
   debug_crystallite, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_g
 use IO, only: &
   IO_warning, &
   IO_error
 use math, only: &
   math_inv33, &
   math_identity2nd, &
   math_mul33x33, &
   math_mul66x6, &
   math_Mandel6to33, &
   math_Mandel33to6, &
   math_Plain3333to99, &
   math_Plain99to3333, &
   math_I3, &
   math_mul3333xx3333, &
   math_mul33xx33, &
   math_invert, &
   math_det33
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP
 use mesh, only: &
   mesh_element, &
   mesh_maxNips, &
   mesh_ipNeighborhood, &
   FE_NipNeighbors, &
   FE_geomtype, &
   FE_cellType
 use material, only: &
   homogenization_Ngrains, &
   plasticState, &
   sourceState, &
   phase_Nsources, &
   phaseAt, phasememberAt
 use constitutive, only:  &
   constitutive_SandItsTangents, &
   constitutive_LpAndItsTangents, &
   constitutive_LiAndItsTangents

 implicit none
 logical, intent(in) :: &
   updateJaco                                                                                       !< whether to update the Jacobian (stiffness) or not
 real(pReal) :: &
   formerSubStep, &
   subFracIntermediate
 real(pReal), dimension(3,3) :: &
   invFp, &                                                                                         ! inverse of the plastic deformation gradient
   Fe_guess, &                                                                                      ! guess for elastic deformation gradient
   Tstar                                                                                            ! 2nd Piola-Kirchhoff stress tensor
 integer(pInt) :: &
   NiterationCrystallite, &                                                                         ! number of iterations in crystallite loop
   c, &                                                                                             !< counter in integration point component loop
   i, &                                                                                             !< counter in integration point loop
   e, &                                                                                             !< counter in element loop
   n, startIP, endIP, &
   neighboring_e, &
   neighboring_i, &
   o, &
   p, &
   myNcomponents, &
   mySource
 ! local variables used for calculating analytic Jacobian
 real(pReal), dimension(3,3)     ::   temp_33
 real(pReal), dimension(3,3,3,3) ::   dSdFe, &
                                      dSdF, &
                                      dSdFi, &
                                      dLidS, &
                                      dLidFi, &
                                      dLpdS, &
                                      dLpdFi, &
                                      dFidS, &
                                      dFpinvdF, &
                                      rhs_3333, &
                                      lhs_3333, &
                                      temp_3333
 real(pReal), dimension(9,9)::        temp_99
 logical :: error


 if (iand(debug_level(debug_crystallite),debug_levelSelective) /= 0_pInt &
     .and. FEsolving_execElem(1) <= debug_e &
     .and.                          debug_e <= FEsolving_execElem(2)) then
     write(6,'(/,a,i8,1x,a,i8,a,1x,i2,1x,i3)')      '<< CRYST >> boundary values at el ip ipc ', &
       debug_e,'(',mesh_element(1,debug_e), ')',debug_i, debug_g
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> F  ', &
                                         transpose(crystallite_partionedF(1:3,1:3,debug_g,debug_i,debug_e))
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> F0 ', &
                                         transpose(crystallite_partionedF0(1:3,1:3,debug_g,debug_i,debug_e))
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Fp0', &
                                         transpose(crystallite_partionedFp0(1:3,1:3,debug_g,debug_i,debug_e))
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Fi0', &
                                         transpose(crystallite_partionedFi0(1:3,1:3,debug_g,debug_i,debug_e))
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Lp0', &
                                         transpose(crystallite_partionedLp0(1:3,1:3,debug_g,debug_i,debug_e))
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Li0', &
                                         transpose(crystallite_partionedLi0(1:3,1:3,debug_g,debug_i,debug_e))
 endif

!--------------------------------------------------------------------------------------------------
! initialize to starting condition
 crystallite_subStep = 0.0_pReal

 !$OMP PARALLEL DO PRIVATE(myNcomponents)
   elementLooping1: do e = FEsolving_execElem(1),FEsolving_execElem(2)
     myNcomponents = homogenization_Ngrains(mesh_element(3,e))
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e); do c = 1_pInt,myNcomponents
       if (crystallite_requested(c,i,e)) then
         plasticState    (phaseAt(c,i,e))%subState0(      :,phasememberAt(c,i,e)) = &
         plasticState    (phaseAt(c,i,e))%partionedState0(:,phasememberAt(c,i,e))

         do mySource = 1_pInt, phase_Nsources(phaseAt(c,i,e))
           sourceState(phaseAt(c,i,e))%p(mySource)%subState0(      :,phasememberAt(c,i,e)) = &
           sourceState(phaseAt(c,i,e))%p(mySource)%partionedState0(:,phasememberAt(c,i,e))
         enddo
         crystallite_subFp0(1:3,1:3,c,i,e) = crystallite_partionedFp0(1:3,1:3,c,i,e)                  ! ...plastic def grad
         crystallite_subLp0(1:3,1:3,c,i,e) = crystallite_partionedLp0(1:3,1:3,c,i,e)                  ! ...plastic velocity grad
         crystallite_subFi0(1:3,1:3,c,i,e) = crystallite_partionedFi0(1:3,1:3,c,i,e)                  ! ...intermediate def grad
         crystallite_subLi0(1:3,1:3,c,i,e) = crystallite_partionedLi0(1:3,1:3,c,i,e)                  ! ...intermediate velocity grad
         crystallite_dPdF0(1:3,1:3,1:3,1:3,c,i,e) = crystallite_partioneddPdF0(1:3,1:3,1:3,1:3,c,i,e) ! ...stiffness
         crystallite_subF0(1:3,1:3,c,i,e) = crystallite_partionedF0(1:3,1:3,c,i,e)                    ! ...def grad
         crystallite_subTstar0_v(1:6,c,i,e) = crystallite_partionedTstar0_v(1:6,c,i,e)                !...2nd PK stress
         crystallite_subFrac(c,i,e) = 0.0_pReal
         crystallite_subStep(c,i,e) = 1.0_pReal/subStepSizeCryst
         crystallite_todo(c,i,e) = .true.
         crystallite_converged(c,i,e) = .false.                                                       ! pretend failed step of twice the required size
       endif
     enddo; enddo
   enddo elementLooping1
 !$OMP END PARALLEL DO

 singleRun: if (FEsolving_execELem(1) == FEsolving_execElem(2) .and. &
     FEsolving_execIP(1,FEsolving_execELem(1))==FEsolving_execIP(2,FEsolving_execELem(1))) then
   startIP = FEsolving_execIP(1,FEsolving_execELem(1))
   endIP   = startIP
 else singleRun
   startIP = 1_pInt
   endIP = mesh_maxNips
 endif singleRun

 NiterationCrystallite = 0_pInt
 cutbackLooping: do while (any(crystallite_todo(:,startIP:endIP,FEsolving_execELem(1):FEsolving_execElem(2))))

   if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) &
     write(6,'(a,i6)') '<< CRYST >> crystallite iteration ',NiterationCrystallite

   timeSyncing1: if (any(.not. crystallite_localPlasticity) .and. numerics_timeSyncing) then

    ! Time synchronization can only be used for nonlocal calculations, and only there it makes sense.
    ! The idea is that in nonlocal calculations often the vast majority of the ips
    ! converges in one iteration whereas a small fraction of ips has to do a lot of cutbacks.
    ! Hence, we try to minimize the computational effort by just doing a lot of cutbacks
    ! in the vicinity of the "bad" ips and leave the easily converged volume more or less as it is.
    ! However, some synchronization of the time step has to be done at the border between "bad" ips
    ! and the ones that immediately converged.

     if (any(crystallite_syncSubFrac)) then

       ! Just did a time synchronization.
       ! If all synchronizers converged, then do nothing else than winding them forward.
       ! If any of the synchronizers did not converge, something went completely wrong
       ! and its not clear how to fix this, so all nonlocals become terminally ill.

       if (any(crystallite_syncSubFrac .and. .not. crystallite_converged(1,:,:))) then
         if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
               if (crystallite_syncSubFrac(i,e) .and. .not. crystallite_converged(1,i,e)) &
                 write(6,'(a,i8,1x,i2)') '<< CRYST >> time synchronization: failed at el,ip ',e,i
             enddo
           enddo
         endif
         crystallite_syncSubFrac = .false.
         where(.not. crystallite_localPlasticity)
           crystallite_substep = 0.0_pReal
           crystallite_todo = .false.
         endwhere
       else
         !$OMP PARALLEL DO
         do e = FEsolving_execElem(1),FEsolving_execElem(2)
           do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
             crystallite_clearToWindForward(i,e) = crystallite_localPlasticity(1,i,e) .or. crystallite_syncSubFrac(i,e)
             crystallite_clearToCutback(i,e) = crystallite_localPlasticity(1,i,e)
           enddo
         enddo
         !$OMP END PARALLEL DO
         if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) &
           write(6,'(a,i6)') '<< CRYST >> time synchronization: wind forward'
       endif

     elseif (any(crystallite_syncSubFracCompleted)) then

       ! Just completed a time synchronization.
       ! Make sure that the ips that synchronized their time step start non-converged

       do e = FEsolving_execElem(1),FEsolving_execElem(2)
         do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
           if (crystallite_syncSubFracCompleted(i,e)) crystallite_converged(1,i,e) = .false.
           crystallite_syncSubFracCompleted(i,e) = .false.
           crystallite_clearToWindForward(i,e) = crystallite_localPlasticity(1,i,e)
           crystallite_clearToCutback(i,e) = crystallite_localPlasticity(1,i,e) .or. .not. crystallite_converged(1,i,e)
         enddo
       enddo
       if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) &
         write(6,'(a,i6)') '<< CRYST >> time synchronization: done, proceed with cutback'
     else

       ! Normal calculation.
       ! If all converged and are at the end of the time increment, then just do a final wind forward.
       ! If all converged, but not all reached the end of the time increment, then we only wind
       ! those forward that are still on their way, all others have to wait.
       ! If some did not converge and all are still at the start of the time increment,
       ! then all non-convergers force their converged neighbors to also do a cutback.
       ! In case that some ips have already wound forward to an intermediate time (subfrac),
       ! then all those ips that converged in the first iteration, but now have a non-converged neighbor
       ! have to synchronize their time step to the same intermediate time. If such a synchronization
       ! takes place, all other ips have to wait and only the synchronizers do a cutback. In the next
       ! iteration those will do a wind forward while all others still wait.

       !$OMP PARALLEL DO
       do e = FEsolving_execElem(1),FEsolving_execElem(2)
         do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
           crystallite_clearToWindForward(i,e) = crystallite_localPlasticity(1,i,e)
           crystallite_clearToCutback(i,e) = crystallite_localPlasticity(1,i,e)
         enddo
       enddo
       !$OMP END PARALLEL DO
       if (all(crystallite_localPlasticity .or. crystallite_converged)) then
         if (all(crystallite_localPlasticity .or. crystallite_subStep + crystallite_subFrac >= 1.0_pReal)) then
           crystallite_clearToWindForward = .true.   ! final wind forward
           if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) &
             write(6,'(a,i6)') '<< CRYST >> final wind forward'
         else
           !$OMP PARALLEL DO
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
               crystallite_clearToWindForward(i,e) = crystallite_localPlasticity(1,i,e) .or. crystallite_subStep(1,i,e) < 1.0_pReal
             enddo
           enddo
           !$OMP END PARALLEL DO
           if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) &
             write(6,'(a,i6)') '<< CRYST >> wind forward'
         endif
       else
         subFracIntermediate = maxval(crystallite_subFrac, mask=.not.crystallite_localPlasticity)
         if (dNeq0(subFracIntermediate)) then
           crystallite_neighborEnforcedCutback = .false.  ! look for ips that require a cutback because of a nonconverged neighbor
           !$OMP PARALLEL
           !$OMP DO PRIVATE(neighboring_e,neighboring_i)
             do e = FEsolving_execElem(1),FEsolving_execElem(2)
               do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
                 if (.not. crystallite_localPlasticity(1,i,e) .and. crystallite_converged(1,i,e)) then
                   do n = 1_pInt,FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,e))))
                     neighboring_e = mesh_ipNeighborhood(1,n,i,e)
                     neighboring_i = mesh_ipNeighborhood(2,n,i,e)
                     if (neighboring_e > 0_pInt .and. neighboring_i > 0_pInt) then
                       if (.not. crystallite_localPlasticity(1,neighboring_i,neighboring_e) &
                           .and. .not. crystallite_converged(1,neighboring_i,neighboring_e)) then
                         crystallite_neighborEnforcedCutback(i,e) = .true.
#ifdef DEBUG
                         if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) &
                            write(6,'(a12,i5,1x,i2,a,i5,1x,i2)') '<< CRYST >> ', neighboring_e,neighboring_i, &
                                                                 ' enforced cutback at ',e,i
#endif
                         exit
                       endif
                     endif
                   enddo
                 endif
               enddo
             enddo
           !$OMP END DO
           !$OMP DO
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
               if(crystallite_neighborEnforcedCutback(i,e)) crystallite_converged(1,i,e) = .false.
             enddo
           enddo
           !$OMP END DO
           !$OMP END PARALLEL
         else
           crystallite_syncSubFrac = .false.  ! look for ips that have to do a time synchronization because of a nonconverged neighbor
           !$OMP PARALLEL
           !$OMP DO PRIVATE(neighboring_e,neighboring_i)
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
               if (.not. crystallite_localPlasticity(1,i,e) .and. dNeq0(crystallite_subFrac(1,i,e))) then
                 do n = 1_pInt,FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,e))))
                   neighboring_e = mesh_ipNeighborhood(1,n,i,e)
                   neighboring_i = mesh_ipNeighborhood(2,n,i,e)
                   if (neighboring_e > 0_pInt .and. neighboring_i > 0_pInt) then
                     if (.not. crystallite_localPlasticity(1,neighboring_i,neighboring_e) &
                         .and. .not. crystallite_converged(1,neighboring_i,neighboring_e)) then
                       crystallite_syncSubFrac(i,e) = .true.
#ifdef DEBUG
                       if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) &
                         write(6,'(a12,i5,1x,i2,a,i5,1x,i2)') '<< CRYST >> ',neighboring_e,neighboring_i, &
                                                             ' enforced time synchronization at ',e,i
#endif
                       exit
                     endif
                   endif
                 enddo
               endif
             enddo
           enddo
           !$OMP END DO
           !$OMP DO
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
               if(crystallite_syncSubFrac(i,e)) crystallite_converged(1,i,e) = .false.
             enddo
           enddo
           !$OMP END DO
           !$OMP END PARALLEL
         endif
         where(.not. crystallite_localPlasticity .and. crystallite_subStep < 1.0_pReal) &
           crystallite_converged = .false.
         if (any(crystallite_syncSubFrac)) then   ! have to do syncing now, so all wait except for the synchronizers which do a cutback
           !$OMP PARALLEL DO
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
               crystallite_clearToWindForward(i,e) = crystallite_localPlasticity(1,i,e)
               crystallite_clearToCutback(i,e) = crystallite_localPlasticity(1,i,e) .or. crystallite_syncSubFrac(i,e)
             enddo
           enddo
           !$OMP END PARALLEL DO
           if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) &
             write(6,'(a,i6)') '<< CRYST >> time synchronization: cutback'
         else
           !$OMP PARALLEL DO
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
               if(.not. crystallite_converged(1,i,e)) crystallite_clearToCutback(i,e) = .true.
             enddo
           enddo
           !$OMP END PARALLEL DO
           if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) &
             write(6,'(a,i6)') '<< CRYST >> cutback'
         endif
       endif
     endif

     ! Make sure that all cutbackers start with the same substep

     where(.not. crystallite_localPlasticity .and. .not. crystallite_converged) &
       crystallite_subStep = minval(crystallite_subStep, mask=.not. crystallite_localPlasticity &
                                                              .and. .not. crystallite_converged)

     ! Those that do neither wind forward nor cutback are not to do

     !$OMP PARALLEL DO
     elementLooping2: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
         if(.not. crystallite_clearToWindForward(i,e) .and. .not. crystallite_clearToCutback(i,e)) &
           crystallite_todo(1,i,e) = .false.
       enddo
     enddo elementLooping2
     !$OMP END PARALLEL DO

   endif timeSyncing1

   !$OMP PARALLEL DO PRIVATE(myNcomponents,formerSubStep)
     elementLooping3: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       myNcomponents = homogenization_Ngrains(mesh_element(3,e))
       do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                            ! iterate over IPs of this element to be processed
         do c = 1,myNcomponents
           ! --- wind forward ---

           if (crystallite_converged(c,i,e) .and. crystallite_clearToWindForward(i,e)) then
             formerSubStep = crystallite_subStep(c,i,e)
             crystallite_subFrac(c,i,e) = crystallite_subFrac(c,i,e) + crystallite_subStep(c,i,e)
             !$OMP FLUSH(crystallite_subFrac)
             crystallite_subStep(c,i,e) = min(1.0_pReal - crystallite_subFrac(c,i,e), &
                                              stepIncreaseCryst * crystallite_subStep(c,i,e))
             !$OMP FLUSH(crystallite_subStep)
             if (crystallite_subStep(c,i,e) > 0.0_pReal) then
               crystallite_subF0(1:3,1:3,c,i,e) = crystallite_subF(1:3,1:3,c,i,e)                    ! ...def grad
               !$OMP FLUSH(crystallite_subF0)
               crystallite_subLp0(1:3,1:3,c,i,e) = crystallite_Lp(1:3,1:3,c,i,e)                     ! ...plastic velocity gradient
               crystallite_subLi0(1:3,1:3,c,i,e) = crystallite_Li(1:3,1:3,c,i,e)                     ! ...intermediate velocity gradient
               crystallite_subFp0(1:3,1:3,c,i,e) = crystallite_Fp(1:3,1:3,c,i,e)                     ! ...plastic def grad
               crystallite_subFi0(1:3,1:3,c,i,e) = crystallite_Fi(1:3,1:3,c,i,e)                     ! ...intermediate def grad
               !if abbrevation, make c and p private in omp
               plasticState    (phaseAt(c,i,e))%subState0(:,phasememberAt(c,i,e)) = &
               plasticState    (phaseAt(c,i,e))%state(    :,phasememberAt(c,i,e))
               do mySource = 1_pInt, phase_Nsources(phaseAt(c,i,e))
                 sourceState(phaseAt(c,i,e))%p(mySource)%subState0(:,phasememberAt(c,i,e)) = &
                 sourceState(phaseAt(c,i,e))%p(mySource)%state(    :,phasememberAt(c,i,e))
               enddo
               crystallite_subTstar0_v(1:6,c,i,e) = crystallite_Tstar_v(1:6,c,i,e)                   ! ...2nd PK stress
               if (crystallite_syncSubFrac(i,e)) then                                                ! if we just did a synchronization of states, then we wind forward without any further time integration
                 crystallite_syncSubFracCompleted(i,e) = .true.
                 crystallite_syncSubFrac(i,e) = .false.
                 crystallite_todo(c,i,e) = .false.
               else
                 crystallite_todo(c,i,e) = .true.
               endif
               !$OMP FLUSH(crystallite_todo)
#ifdef DEBUG
               if (iand(debug_level(debug_crystallite),debug_levelBasic) /= 0_pInt &
                   .and. ((e == debug_e .and. i == debug_i .and. c == debug_g) &
                          .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) &
                 write(6,'(a,f12.8,a,f12.8,a,i8,1x,i2,1x,i3,/)') '<< CRYST >> winding forward from ', &
                   crystallite_subFrac(c,i,e)-formerSubStep,' to current crystallite_subfrac ', &
                   crystallite_subFrac(c,i,e),' in crystallite_stressAndItsTangent at el ip ipc ',e,i,c
#endif
             else                                                                                    ! this crystallite just converged for the entire timestep
               crystallite_todo(c,i,e) = .false.                                                     ! so done here
               !$OMP FLUSH(crystallite_todo)
             endif

           ! --- cutback ---

           elseif (.not. crystallite_converged(c,i,e) .and. crystallite_clearToCutback(i,e)) then
             if (crystallite_syncSubFrac(i,e)) then                                                  ! synchronize time
               crystallite_subStep(c,i,e) = subFracIntermediate
             else
               crystallite_subStep(c,i,e) = subStepSizeCryst * crystallite_subStep(c,i,e)            ! cut step in half and restore...
             endif
             !$OMP FLUSH(crystallite_subStep)
             crystallite_Fp(1:3,1:3,c,i,e) = crystallite_subFp0(1:3,1:3,c,i,e)                       ! ...plastic def grad
             !$OMP FLUSH(crystallite_Fp)
             crystallite_invFp(1:3,1:3,c,i,e) = math_inv33(crystallite_Fp(1:3,1:3,c,i,e))
             !$OMP FLUSH(crystallite_invFp)
             crystallite_Fi(1:3,1:3,c,i,e) = crystallite_subFi0(1:3,1:3,c,i,e)                       ! ...intermediate def grad
             !$OMP FLUSH(crystallite_Fi)
             crystallite_invFi(1:3,1:3,c,i,e) = math_inv33(crystallite_Fi(1:3,1:3,c,i,e))
             !$OMP FLUSH(crystallite_invFi)
             crystallite_Lp(1:3,1:3,c,i,e)    = crystallite_subLp0(1:3,1:3,c,i,e)                    ! ...plastic velocity grad
             crystallite_Li(1:3,1:3,c,i,e)    = crystallite_subLi0(1:3,1:3,c,i,e)                    ! ...intermediate velocity grad
             plasticState    (phaseAt(c,i,e))%state(    :,phasememberAt(c,i,e)) = &
             plasticState    (phaseAt(c,i,e))%subState0(:,phasememberAt(c,i,e))
             do mySource = 1_pInt, phase_Nsources(phaseAt(c,i,e))
               sourceState(phaseAt(c,i,e))%p(mySource)%state(    :,phasememberAt(c,i,e)) = &
               sourceState(phaseAt(c,i,e))%p(mySource)%subState0(:,phasememberAt(c,i,e))
             enddo
             crystallite_Tstar_v(1:6,c,i,e)   = crystallite_subTstar0_v(1:6,c,i,e)                   ! ...2nd PK stress

                                                                                                     ! cant restore dotState here, since not yet calculated in first cutback after initialization
             crystallite_todo(c,i,e) = crystallite_subStep(c,i,e) > subStepMinCryst                  ! still on track or already done (beyond repair)
             !$OMP FLUSH(crystallite_todo)
#ifdef DEBUG
             if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
                .and. ((e == debug_e .and. i == debug_i .and. c == debug_g) &
                       .or. .not. iand(debug_level(debug_crystallite),debug_levelSelective) /= 0_pInt)) then
               if (crystallite_todo(c,i,e)) then
                 write(6,'(a,f12.8,a,i8,1x,i2,1x,i3,/)') '<< CRYST >> cutback step in crystallite_stressAndItsTangent &
                                                        &with new crystallite_subStep: ',&
                                                       crystallite_subStep(c,i,e),' at el ip ipc ',e,i,c
               else
                 write(6,'(a,i8,1x,i2,1x,i3,/)') '<< CRYST >> reached minimum step size &
                                               &in crystallite_stressAndItsTangent at el ip ipc ',e,i,c
               endif
             endif
#endif
           endif

           ! --- prepare for integration ---

           if (crystallite_todo(c,i,e) .and. (crystallite_clearToWindForward(i,e) .or. crystallite_clearToCutback(i,e))) then
             crystallite_subF(1:3,1:3,c,i,e) = crystallite_subF0(1:3,1:3,c,i,e) &
                                             + crystallite_subStep(c,i,e) &
                                               * (crystallite_partionedF(1:3,1:3,c,i,e) &
                                             - crystallite_partionedF0(1:3,1:3,c,i,e))
             !$OMP FLUSH(crystallite_subF)
             crystallite_Fe(1:3,1:3,c,i,e) = math_mul33x33(math_mul33x33(crystallite_subF (1:3,1:3,c,i,e), &
                                                                         crystallite_invFp(1:3,1:3,c,i,e)), &
                                                           crystallite_invFi(1:3,1:3,c,i,e))
             crystallite_subdt(c,i,e) = crystallite_subStep(c,i,e) * crystallite_dt(c,i,e)
             crystallite_converged(c,i,e) = .false.                                                  ! start out non-converged
           endif

         enddo    ! grains
       enddo      ! IPs
     enddo elementLooping3
   !$OMP END PARALLEL DO

   timeSyncing2: if(numerics_timeSyncing) then
     if (any(.not. crystallite_localPlasticity .and. .not. crystallite_todo .and. .not. crystallite_converged &
             .and. crystallite_subStep <= subStepMinCryst)) then                                     ! no way of rescuing a nonlocal ip that violated the lower time step limit, ...
       if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
         elementLooping4: do e = FEsolving_execElem(1),FEsolving_execElem(2)
           myNcomponents = homogenization_Ngrains(mesh_element(3,e))
           do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
             do c = 1,myNcomponents
               if (.not. crystallite_localPlasticity(c,i,e) .and. .not. crystallite_todo(c,i,e) &
                   .and. .not. crystallite_converged(c,i,e) .and. crystallite_subStep(c,i,e) <= subStepMinCryst) &
                 write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> nonlocal violated minimum subStep at el ip ipc ',e,i,c
             enddo
           enddo
         enddo elementLooping4
       endif
       where(.not. crystallite_localPlasticity)
         crystallite_todo = .false.                                                                  ! ... so let all nonlocal ips die peacefully
         crystallite_subStep = 0.0_pReal
       endwhere
     endif
   endif timeSyncing2

   if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
     write(6,'(/,a,f8.5)') '<< CRYST >> min(subStep) ',minval(crystallite_subStep)
     write(6,'(a,f8.5)')   '<< CRYST >> max(subStep) ',maxval(crystallite_subStep)
     write(6,'(a,f8.5)')   '<< CRYST >> min(subFrac) ',minval(crystallite_subFrac)
     write(6,'(a,f8.5,/)') '<< CRYST >> max(subFrac) ',maxval(crystallite_subFrac)
     flush(6)
     if (iand(debug_level(debug_crystallite),debug_levelSelective) /= 0_pInt) then
       write(6,'(/,a,f8.5,1x,a,1x,f8.5,1x,a)') '<< CRYST >> subFrac + subStep = ',&
          crystallite_subFrac(debug_g,debug_i,debug_e),'+',crystallite_subStep(debug_g,debug_i,debug_e),'@selective'
       flush(6)
     endif
   endif

   ! --- integrate --- requires fully defined state array (basic + dependent state)

   if (any(crystallite_todo)) call integrateState()
   where(.not. crystallite_converged .and. crystallite_subStep > subStepMinCryst) &                  ! do not try non-converged & fully cutbacked any further
     crystallite_todo = .true.

   NiterationCrystallite = NiterationCrystallite + 1_pInt

 enddo cutbackLooping


! --+>> CHECK FOR NON-CONVERGED CRYSTALLITES <<+--

 elementLooping5: do e = FEsolving_execElem(1),FEsolving_execElem(2)
   myNcomponents = homogenization_Ngrains(mesh_element(3,e))
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                              ! iterate over IPs of this element to be processed
     do c = 1,myNcomponents
       if (.not. crystallite_converged(c,i,e)) then                                                ! respond fully elastically (might be not required due to becoming terminally ill anyway)
         if(iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) &
           write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3,/)') '<< CRYST >> no convergence: respond fully elastic at el (elFE) ip ipc ', &
             e,'(',mesh_element(1,e),')',i,c
         invFp = math_inv33(crystallite_partionedFp0(1:3,1:3,c,i,e))
         Fe_guess = math_mul33x33(math_mul33x33(crystallite_partionedF(1:3,1:3,c,i,e), invFp), &
                                  math_inv33(crystallite_partionedFi0(1:3,1:3,c,i,e)))
         call constitutive_SandItsTangents(Tstar,dSdFe,dSdFi,Fe_guess,crystallite_partionedFi0(1:3,1:3,c,i,e),c,i,e)
         crystallite_P(1:3,1:3,c,i,e) = math_mul33x33(math_mul33x33(crystallite_partionedF(1:3,1:3,c,i,e), invFp), &
                                                      math_mul33x33(Tstar,transpose(invFp)))
       endif
       if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
           .and. ((e == debug_e .and. i == debug_i .and. c == debug_g) &
                  .or. .not. iand(debug_level(debug_crystallite),debug_levelSelective) /= 0_pInt)) then
         write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> central solution of cryst_StressAndTangent at el ip ipc ',e,i,c
         write(6,'(/,a,/,3(12x,3(f12.4,1x)/))') '<< CRYST >> P / MPa', &
                                          transpose(crystallite_P(1:3,1:3,c,i,e))*1.0e-6_pReal
         write(6,'(a,/,3(12x,3(f14.9,1x)/))')   '<< CRYST >> Fp', &
                                          transpose(crystallite_Fp(1:3,1:3,c,i,e))
         write(6,'(a,/,3(12x,3(f14.9,1x)/))')   '<< CRYST >> Fi', &
                                          transpose(crystallite_Fi(1:3,1:3,c,i,e))
         write(6,'(a,/,3(12x,3(f14.9,1x)/),/)') '<< CRYST >> Lp', &
                                          transpose(crystallite_Lp(1:3,1:3,c,i,e))
         write(6,'(a,/,3(12x,3(f14.9,1x)/),/)') '<< CRYST >> Li', &
                                          transpose(crystallite_Li(1:3,1:3,c,i,e))
         flush(6)
       endif
     enddo
   enddo
 enddo elementLooping5


! --+>> STIFFNESS CALCULATION <<+--

 computeJacobian: if(updateJaco) then
   !$OMP PARALLEL DO PRIVATE(dSdF,dSdFe,dSdFi,dLpdS,dLpdFi,dFpinvdF,dLidS,dLidFi,dFidS,&
   !$OMP                     rhs_3333,lhs_3333,temp_99,temp_33,temp_3333,myNcomponents,error)
   elementLooping6: do e = FEsolving_execElem(1),FEsolving_execElem(2)
     myNcomponents = homogenization_Ngrains(mesh_element(3,e))
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                            ! iterate over IPs of this element to be processed
       do c = 1_pInt,myNcomponents
         call constitutive_SandItsTangents(temp_33,dSdFe,dSdFi,crystallite_Fe(1:3,1:3,c,i,e), &
                                          crystallite_Fi(1:3,1:3,c,i,e),c,i,e)                     ! call constitutive law to calculate elastic stress tangent

         call constitutive_LiAndItsTangents(temp_33,dLidS,dLidFi,crystallite_Tstar_v(1:6,c,i,e), &
                                           crystallite_Fi(1:3,1:3,c,i,e), &
                                           c,i,e)                                                  ! call constitutive law to calculate Li tangent in lattice configuration
         if (sum(abs(dLidS)) < tol_math_check) then
           dFidS = 0.0_pReal
         else
           temp_33 = math_inv33(crystallite_subFi0(1:3,1:3,c,i,e))
           lhs_3333 = 0.0_pReal; rhs_3333 = 0.0_pReal
           do o=1_pInt,3_pInt; do p=1_pInt,3_pInt
             lhs_3333(1:3,1:3,o,p) = lhs_3333(1:3,1:3,o,p) + &
                                     crystallite_subdt(c,i,e)*math_mul33x33(temp_33,dLidFi(1:3,1:3,o,p))
             lhs_3333(1:3,o,1:3,p) = lhs_3333(1:3,o,1:3,p) + &
                                     crystallite_invFi(1:3,1:3,c,i,e)*crystallite_invFi(p,o,c,i,e)
             rhs_3333(1:3,1:3,o,p) = rhs_3333(1:3,1:3,o,p) - &
                                     crystallite_subdt(c,i,e)*math_mul33x33(temp_33,dLidS(1:3,1:3,o,p))
           enddo; enddo
           call math_invert(9_pInt,math_Plain3333to99(lhs_3333),temp_99,error)
           if (error) then
             call IO_warning(warning_ID=600_pInt,el=e,ip=i,g=c, &
                             ext_msg='inversion error in analytic tangent calculation')
             dFidS = 0.0_pReal
           else
             dFidS = math_mul3333xx3333(math_Plain99to3333(temp_99),rhs_3333)
           endif
           dLidS = math_mul3333xx3333(dLidFi,dFidS) + dLidS
         endif

         call constitutive_LpAndItsTangents(temp_33,dLpdS,dLpdFi,crystallite_Tstar_v(1:6,c,i,e), &
                                           crystallite_Fi(1:3,1:3,c,i,e),c,i,e)                    ! call constitutive law to calculate Lp tangent in lattice configuration
         dLpdS = math_mul3333xx3333(dLpdFi,dFidS) + dLpdS

         temp_33   = transpose(math_mul33x33(crystallite_invFp(1:3,1:3,c,i,e), &
                                                    crystallite_invFi(1:3,1:3,c,i,e)))
         rhs_3333 = 0.0_pReal
         forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
           rhs_3333(p,o,1:3,1:3) = math_mul33x33(dSdFe(p,o,1:3,1:3),temp_33)

         temp_3333 = 0.0_pReal
         temp_33 = math_mul33x33(crystallite_subF(1:3,1:3,c,i,e), &
                                 math_inv33(crystallite_subFp0(1:3,1:3,c,i,e)))
         forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
           temp_3333(1:3,1:3,p,o) = math_mul33x33(math_mul33x33(temp_33,dLpdS(1:3,1:3,p,o)), &
                                                  crystallite_invFi(1:3,1:3,c,i,e))

         temp_33 = math_mul33x33(math_mul33x33(crystallite_subF(1:3,1:3,c,i,e), &
                                               crystallite_invFp(1:3,1:3,c,i,e)), &
                                 math_inv33(crystallite_subFi0(1:3,1:3,c,i,e)))
         forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
           temp_3333(1:3,1:3,p,o) = temp_3333(1:3,1:3,p,o) + math_mul33x33(temp_33,dLidS(1:3,1:3,p,o))

         lhs_3333 = crystallite_subdt(c,i,e)*math_mul3333xx3333(dSdFe,temp_3333) + &
                    math_mul3333xx3333(dSdFi,dFidS)

         call math_invert(9_pInt,math_identity2nd(9_pInt)+math_Plain3333to99(lhs_3333),temp_99,error)
         if (error) then
           call IO_warning(warning_ID=600_pInt,el=e,ip=i,g=c, &
                           ext_msg='inversion error in analytic tangent calculation')
           dSdF = rhs_3333
         else
           dSdF = math_mul3333xx3333(math_Plain99to3333(temp_99),rhs_3333)
         endif

         dFpinvdF = 0.0_pReal
         temp_3333 = math_mul3333xx3333(dLpdS,dSdF)
         forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
           dFpinvdF(1:3,1:3,p,o) = -crystallite_subdt(c,i,e)* &
                                    math_mul33x33(math_inv33(crystallite_subFp0(1:3,1:3,c,i,e)), &
                                                  math_mul33x33(temp_3333(1:3,1:3,p,o), &
                                                                crystallite_invFi(1:3,1:3,c,i,e)))

         crystallite_dPdF(1:3,1:3,1:3,1:3,c,i,e) = 0.0_pReal
         temp_33 = math_mul33x33(crystallite_invFp(1:3,1:3,c,i,e), &
                              math_mul33x33(math_Mandel6to33(crystallite_Tstar_v(1:6,c,i,e)), &
                                            transpose(crystallite_invFp(1:3,1:3,c,i,e))))
         forall(p=1_pInt:3_pInt) &
           crystallite_dPdF(p,1:3,p,1:3,c,i,e) = transpose(temp_33)

         temp_33 = math_mul33x33(math_Mandel6to33(crystallite_Tstar_v(1:6,c,i,e)), &
                              transpose(crystallite_invFp(1:3,1:3,c,i,e)))
         forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
           crystallite_dPdF(1:3,1:3,p,o,c,i,e) = crystallite_dPdF(1:3,1:3,p,o,c,i,e) + &
             math_mul33x33(math_mul33x33(crystallite_subF(1:3,1:3,c,i,e),dFpinvdF(1:3,1:3,p,o)),temp_33)

         temp_33 = math_mul33x33(crystallite_subF(1:3,1:3,c,i,e), &
                              crystallite_invFp(1:3,1:3,c,i,e))
         forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
           crystallite_dPdF(1:3,1:3,p,o,c,i,e) = crystallite_dPdF(1:3,1:3,p,o,c,i,e) + &
             math_mul33x33(math_mul33x33(temp_33,dSdF(1:3,1:3,p,o)), &
                           transpose(crystallite_invFp(1:3,1:3,c,i,e)))

         temp_33 = math_mul33x33(math_mul33x33(crystallite_subF(1:3,1:3,c,i,e), &
                                            crystallite_invFp(1:3,1:3,c,i,e)), &
                              math_Mandel6to33(crystallite_Tstar_v(1:6,c,i,e)))
         forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
           crystallite_dPdF(1:3,1:3,p,o,c,i,e) = crystallite_dPdF(1:3,1:3,p,o,c,i,e) + &
             math_mul33x33(temp_33,transpose(dFpinvdF(1:3,1:3,p,o)))

     enddo; enddo
   enddo elementLooping6
   !$OMP END PARALLEL DO
 endif computeJacobian
!why not OMP?

end subroutine crystallite_stressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with 4th order explicit Runge Kutta method
!--------------------------------------------------------------------------------------------------
subroutine crystallite_integrateStateRK4()
 use, intrinsic :: &
   IEEE_arithmetic
 use debug, only: &
#ifdef DEBUG
   debug_e, &
   debug_i, &
   debug_g, &
#endif
   debug_level, &
   debug_crystallite, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP
 use mesh, only: &
   mesh_element, &
   mesh_NcpElems
 use material, only: &
   homogenization_Ngrains, &
   plasticState, &
   sourceState, &
   phase_Nsources, &
   phaseAt, phasememberAt
 use config, only: &
   material_Nphase
 use constitutive, only: &
   constitutive_collectDotState, &
   constitutive_microstructure

 implicit none
 real(pReal), dimension(4), parameter :: &
   TIMESTEPFRACTION = [0.5_pReal, 0.5_pReal, 1.0_pReal, 1.0_pReal]                                   ! factor giving the fraction of the original timestep used for Runge Kutta Integration
 real(pReal), dimension(4), parameter :: &
   WEIGHT = [1.0_pReal, 2.0_pReal, 2.0_pReal, 1.0_pReal/6.0_pReal]                                   ! weight of slope used for Runge Kutta integration (final weight divided by 6)

 integer(pInt) ::                              e, &                                                  ! element index in element loop
                                               i, &                                                  ! integration point index in ip loop
                                               g, &                                                  ! grain index in grain loop
                                               p, &                                                  ! phase loop
                                               c, &
                                               n, &
                                               mySource, &
                                               mySizePlasticDotState, &
                                               mySizeSourceDotState
 integer(pInt), dimension(2) ::                eIter                                                 ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) ::  iIter, &                                              ! bounds for ip iteration
                                               gIter                                                 ! bounds for grain iteration
 logical   ::                                  NaN, &
                                               singleRun                                             ! flag indicating computation for single (g,i,e) triple

 eIter = FEsolving_execElem(1:2)
 do e = eIter(1),eIter(2)
   iIter(1:2,e) = FEsolving_execIP(1:2,e)
   gIter(1:2,e) = [ 1_pInt,homogenization_Ngrains(mesh_element(3,e))]
 enddo

 singleRun = (eIter(1) == eIter(2) .and. iIter(1,eIter(1)) == iIter(2,eIter(2)))

!--------------------------------------------------------------------------------------------------
! initialize dotState
 if (.not. singleRun) then
   do p = 1_pInt, material_Nphase
     plasticState(p)%RK4dotState = 0.0_pReal
     do mySource = 1_pInt, phase_Nsources(p)
       sourceState(p)%p(mySource)%RK4dotState = 0.0_pReal
     enddo
   enddo
 else
   e = eIter(1)
   i = iIter(1,e)
   do g = gIter(1,e), gIter(2,e)
     plasticState(phaseAt(g,i,e))%RK4dotState(:,phasememberAt(g,i,e)) = 0.0_pReal
     do mySource = 1_pInt, phase_Nsources(phaseAt(g,i,e))
       sourceState(phaseAt(g,i,e))%p(mySource)%RK4dotState(:,phasememberAt(g,i,e)) = 0.0_pReal
     enddo
   enddo
 endif

!--------------------------------------------------------------------------------------------------
! first Runge-Kutta step
 !$OMP PARALLEL
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                 ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) &
       call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), &
                                         crystallite_Fe, &
                                         crystallite_Fi(1:3,1:3,g,i,e), &
                                         crystallite_Fp, &
                                         crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
   enddo; enddo; enddo
 !$OMP ENDDO

 !$OMP DO PRIVATE(p,c,NaN)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                 ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       c = phasememberAt(g,i,e)
       p = phaseAt(g,i,e)
       NaN = any(IEEE_is_NaN(plasticState(p)%dotState(:,c)))
       do mySource = 1_pInt, phase_Nsources(p)
         NaN = NaN .or. any(IEEE_is_NaN(sourceState(p)%p(mySource)%dotState(:,c)))
       enddo
       if (NaN) then                                                                                    ! NaN occured in any dotState
         if (.not. crystallite_localPlasticity(g,i,e)) then                                             ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                      ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         else                                                                                           ! if broken local...
           crystallite_todo(g,i,e) = .false.                                                            ! ... skip this one next time
         endif
       endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 !$OMP END PARALLEL
!--------------------------------------------------------------------------------------------------
! --- SECOND TO FOURTH RUNGE KUTTA STEP PLUS FINAL INTEGRATION ---

 do n = 1_pInt,4_pInt
   ! --- state update ---

   !$OMP PARALLEL
   !$OMP DO PRIVATE(p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         p = phaseAt(g,i,e)
         c = phasememberAt(g,i,e)
         plasticState(p)%RK4dotState(:,c) = plasticState(p)%RK4dotState(:,c) &
                                          + weight(n)*plasticState(p)%dotState(:,c)
         do mySource = 1_pInt, phase_Nsources(p)
           sourceState(p)%p(mySource)%RK4dotState(:,c) = sourceState(p)%p(mySource)%RK4dotState(:,c) &
                                                       + weight(n)*sourceState(p)%p(mySource)%dotState(:,c)
         enddo
       endif
     enddo; enddo; enddo
   !$OMP ENDDO

   !$OMP DO PRIVATE(mySizePlasticDotState,mySizeSourceDotState,p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then

         p = phaseAt(g,i,e)
         c = phasememberAt(g,i,e)
         mySizePlasticDotState = plasticState(p)%sizeDotState
         plasticState(p)%state    (1:mySizePlasticDotState,c) = &
         plasticState(p)%subState0(1:mySizePlasticDotState,c) &
       + plasticState(p)%dotState (1:mySizePlasticDotState,c) &
       * crystallite_subdt(g,i,e) * timeStepFraction(n)
         do mySource = 1_pInt, phase_Nsources(p)
           mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
           sourceState(p)%p(mySource)%state    (1:mySizeSourceDotState,c)  = &
           sourceState(p)%p(mySource)%subState0(1:mySizeSourceDotState,c) &
         + sourceState(p)%p(mySource)%dotState (1:mySizeSourceDotState,c) &
         * crystallite_subdt(g,i,e) * timeStepFraction(n)
         enddo

#ifdef DEBUG
         if (n == 4 &
             .and. iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
             .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                    .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then ! final integration step

           write(6,'(a,i8,1x,i2,1x,i3,/)')       '<< CRYST >> updateState at el ip g ',e,i,g
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> dotState',  plasticState(p)%dotState(1:mySizePlasticDotState,c)
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> new state', plasticState(p)%state(1:mySizePlasticDotState,c)
         endif
#endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- state jump ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         crystallite_todo(g,i,e) = crystallite_stateJump(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then              ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                           ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- update dependent states ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) &
         call constitutive_microstructure(crystallite_orientation,       &
                                          crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), &
                                          g, i, e)                                                         ! update dependent state variables to be consistent with basic states
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- stress integration ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         crystallite_todo(g,i,e) = crystallite_integrateStress(g,i,e,timeStepFraction(n))                  ! fraction of original times step
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then            ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- dot state and RK dot state---

   first3steps: if (n < 4) then
     !$OMP DO
       do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                ! iterate over elements, ips and grains
         if (crystallite_todo(g,i,e)) &
           call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), &
                                             crystallite_Fe, &
                                             crystallite_Fi(1:3,1:3,g,i,e), &
                                             crystallite_Fp, &
                                             timeStepFraction(n)*crystallite_subdt(g,i,e), &               ! fraction of original timestep
                                             crystallite_subFrac, g,i,e)
       enddo; enddo; enddo
     !$OMP ENDDO

     !$OMP DO PRIVATE(p,c,NaN)
       do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                ! iterate over elements, ips and grains
         !$OMP FLUSH(crystallite_todo)
         if (crystallite_todo(g,i,e)) then

           p = phaseAt(g,i,e)
           c = phasememberAt(g,i,e)
           NaN = any(IEEE_is_NaN(plasticState(p)%dotState(:,c)))
           do mySource = 1_pInt, phase_Nsources(p)
             NaN = NaN .or. any(IEEE_is_NaN(sourceState(p)%p(mySource)%dotState(:,c)))
           enddo
           if (NaN) then                                                                                   ! NaN occured in any dotState
             if (.not. crystallite_localPlasticity(g,i,e)) then                                            ! if broken non-local...
               !$OMP CRITICAL (checkTodo)
                 crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                     ! ...all non-locals skipped
               !$OMP END CRITICAL (checkTodo)
             else                                                                                          ! if broken local...
               crystallite_todo(g,i,e) = .false.                                                           ! ... skip this one next time
             endif
           endif
         endif
       enddo; enddo; enddo
     !$OMP ENDDO
   endif first3steps
 !$OMP END PARALLEL

 enddo


 ! --- SET CONVERGENCE FLAG ---

 do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                     ! iterate over elements, ips and grains
   crystallite_converged(g,i,e) = crystallite_todo(g,i,e) .or. crystallite_converged(g,i,e)               ! if still "to do" then converged per definitionem
 enddo; enddo; enddo


 ! --- CHECK NONLOCAL CONVERGENCE ---

 if (.not. singleRun) then                                                                                 ! if not requesting Integration of just a single IP
   if (any(.not. crystallite_converged .and. .not. crystallite_localPlasticity)) then                      ! any non-local not yet converged (or broken)...
     crystallite_converged = crystallite_converged .and. crystallite_localPlasticity                       ! ...restart all non-local as not converged
   endif
 endif

end subroutine crystallite_integrateStateRK4


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with 5th order Runge-Kutta Cash-Karp method with
!> adaptive step size  (use 5th order solution to advance = "local extrapolation")
!--------------------------------------------------------------------------------------------------
subroutine crystallite_integrateStateRKCK45()
 use, intrinsic :: &
   IEEE_arithmetic
 use debug, only: &
#ifdef DEBUG
   debug_e, &
   debug_i, &
   debug_g, &
#endif
   debug_level, &
   debug_crystallite, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective
 use numerics, only: &
   rTol_crystalliteState
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP
 use mesh, only: &
   mesh_element, &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_Ngrains, &
   plasticState, &
   sourceState, &
   phase_Nsources, &
   phaseAt, phasememberAt, &
   homogenization_maxNgrains
 use constitutive, only: &
   constitutive_collectDotState, &
   constitutive_plasticity_maxSizeDotState, &
   constitutive_source_maxSizeDotState, &
   constitutive_microstructure

 implicit none
 real(pReal), dimension(5,5), parameter :: &
   A = reshape([&
     .2_pReal,  .075_pReal,  .3_pReal, -11.0_pReal/54.0_pReal,  1631.0_pReal/55296.0_pReal, &
     .0_pReal,  .225_pReal, -.9_pReal,   2.5_pReal,             175.0_pReal/512.0_pReal, &
     .0_pReal,   .0_pReal,  1.2_pReal, -70.0_pReal/27.0_pReal,  575.0_pReal/13824.0_pReal, &
     .0_pReal,   .0_pReal,   .0_pReal,  35.0_pReal/27.0_pReal,  44275.0_pReal/110592.0_pReal, &
     .0_pReal,   .0_pReal,   .0_pReal,    .0_pReal,             253.0_pReal/4096.0_pReal], &
     [5,5], order=[2,1])                                                                            !< coefficients in Butcher tableau (used for preliminary integration in stages 2 to 6)

 real(pReal), dimension(6), parameter :: &
   B = &
     [37.0_pReal/378.0_pReal, .0_pReal, 250.0_pReal/621.0_pReal, &
     125.0_pReal/594.0_pReal, .0_pReal, 512.0_pReal/1771.0_pReal], &                                !< coefficients in Butcher tableau (used for final integration and error estimate)
   DB = B - &
     [2825.0_pReal/27648.0_pReal, .0_pReal,                   18575.0_pReal/48384.0_pReal,&
     13525.0_pReal/55296.0_pReal, 277.0_pReal/14336.0_pReal,  0.25_pReal]                           !< coefficients in Butcher tableau (used for final integration and error estimate)

 real(pReal), dimension(5), parameter :: &
   C = [0.2_pReal, 0.3_pReal, 0.6_pReal, 1.0_pReal, 0.875_pReal]                                    !< coefficients in Butcher tableau (fractions of original time step in stages 2 to 6)

 integer(pInt) :: &
   e, &                                                                                             ! element index in element loop
   i, &                                                                                             ! integration point index in ip loop
   g, &                                                                                             ! grain index in grain loop
   stage, &                                                                                         ! stage index in integration stage loop
   s, &                                                                                             ! state index
   n, &
   p, &
   cc, &
   mySource, &
   mySizePlasticDotState, &                                                                         ! size of dot States
   mySizeSourceDotState
 integer(pInt), dimension(2) :: &
   eIter                                                                                            ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) :: &
   iIter, &                                                                                         ! bounds for ip iteration
   gIter                                                                                            ! bounds for grain iteration

 real(pReal), dimension(constitutive_plasticity_maxSizeDotState,            &
                        homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   plasticStateResiduum, &                                                                          ! residuum from evolution in microstructure
   relPlasticStateResiduum                                                                          ! relative residuum from evolution in microstructure
 real(pReal), dimension(constitutive_source_maxSizeDotState, &
                        maxval(phase_Nsources), &
                        homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   sourceStateResiduum, &                                                                           ! residuum from evolution in microstructure
   relSourceStateResiduum                                                                           ! relative residuum from evolution in microstructure
 logical :: &
   NaN, &
   singleRun                                                                                        ! flag indicating computation for single (g,i,e) triple

 eIter = FEsolving_execElem(1:2)
 if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) &
   write(6,'(a,1x,i1)') '<< CRYST >> Runge--Kutta step',1

 ! --- LOOP ITERATOR FOR ELEMENT, GRAIN, IP ---
 do e = eIter(1),eIter(2)
   iIter(1:2,e) = FEsolving_execIP(1:2,e)
   gIter(1:2,e) = [ 1_pInt,homogenization_Ngrains(mesh_element(3,e))]
 enddo

 singleRun = (eIter(1) == eIter(2) .and. iIter(1,eIter(1)) == iIter(2,eIter(2)))



 ! --- FIRST RUNGE KUTTA STEP ---

 !$OMP PARALLEL
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) &
       call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), &
                                         crystallite_Fe, &
                                         crystallite_Fi(1:3,1:3,g,i,e), &
                                         crystallite_Fp, &
                                         crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
   enddo; enddo; enddo
 !$OMP ENDDO
 !$OMP DO PRIVATE(p,cc,NaN)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       cc = phasememberAt(g,i,e)
       p = phaseAt(g,i,e)
       NaN = any(IEEE_is_NaN(plasticState(p)%dotState(:,cc)))
       do mySource = 1_pInt, phase_Nsources(p)
         NaN = NaN .or. any(IEEE_is_NaN(sourceState(p)%p(mySource)%dotState(:,cc)))
       enddo
       if (NaN) then                                                                                       ! NaN occured in any dotState
         if (.not. crystallite_localPlasticity(g,i,e)) then                                                ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         else                                                                                              ! if broken local...
           crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
         endif
       endif
     endif
  enddo; enddo; enddo
 !$OMP ENDDO
 !$OMP END PARALLEL


 ! --- SECOND TO SIXTH RUNGE KUTTA STEP ---

 do stage = 1_pInt,5_pInt

   ! --- state update ---

   !$OMP PARALLEL
   !$OMP DO PRIVATE(p,cc)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         p = phaseAt(g,i,e)
         cc = phasememberAt(g,i,e)
         plasticState(p)%RKCK45dotState(stage,:,cc) = plasticState(p)%dotState(:,cc)                       ! store Runge-Kutta dotState
         do mySource = 1_pInt, phase_Nsources(p)
           sourceState(p)%p(mySource)%RKCK45dotState(stage,:,cc) = sourceState(p)%p(mySource)%dotState(:,cc)
         enddo
       endif
     enddo; enddo; enddo
   !$OMP ENDDO

   !$OMP DO PRIVATE(p,cc,n)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         p = phaseAt(g,i,e)
         cc = phasememberAt(g,i,e)

         plasticState(p)%dotState(:,cc) = A(1,stage) * plasticState(p)%RKCK45dotState(1,:,cc)
         do mySource = 1_pInt, phase_Nsources(p)
           sourceState(p)%p(mySource)%dotState(:,cc) = A(1,stage) * sourceState(p)%p(mySource)%RKCK45dotState(1,:,cc)
         enddo
         do n = 2_pInt, stage
           plasticState(p)%dotState(:,cc) = &
           plasticState(p)%dotState(:,cc) + A(n,stage) * plasticState(p)%RKCK45dotState(n,:,cc)
           do mySource = 1_pInt, phase_Nsources(p)
             sourceState(p)%p(mySource)%dotState(:,cc) = &
             sourceState(p)%p(mySource)%dotState(:,cc) + A(n,stage) * sourceState(p)%p(mySource)%RKCK45dotState(n,:,cc)
           enddo
         enddo
       endif
     enddo; enddo; enddo
   !$OMP ENDDO

   !$OMP DO PRIVATE(mySizePlasticDotState,mySizeSourceDotState,p,cc)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         p = phaseAt(g,i,e)
         cc = phasememberAt(g,i,e)
         mySizePlasticDotState = plasticState(p)%sizeDotState
         plasticState    (p)%state    (1:mySizePlasticDotState,    cc) = &
         plasticState    (p)%subState0(1:mySizePlasticDotState,    cc) &
       + plasticState    (p)%dotState (1:mySizePlasticDotState,    cc) &
       * crystallite_subdt(g,i,e)
         do mySource = 1_pInt, phase_Nsources(p)
           mySizeSourceDotState  = sourceState(p)%p(mySource)%sizeDotState
           sourceState(p)%p(mySource)%state    (1:mySizeSourceDotState,cc) = &
           sourceState(p)%p(mySource)%subState0(1:mySizeSourceDotState,cc) &
         + sourceState(p)%p(mySource)%dotState (1:mySizeSourceDotState,cc) &
         * crystallite_subdt(g,i,e)
         enddo
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- state jump ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                   ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         crystallite_todo(g,i,e) = crystallite_stateJump(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then             ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                          ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- update dependent states ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) &
         call constitutive_microstructure(crystallite_orientation,       &
                                          crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), &
                                          g, i, e)                                                         ! update dependent state variables to be consistent with basic states
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- stress integration ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         crystallite_todo(g,i,e) =  crystallite_integrateStress(g,i,e,C(stage))                            ! fraction of original time step
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then            ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- dot state and RK dot state---
#ifdef DEBUG
   if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) &
     write(6,'(a,1x,i1)') '<< CRYST >> Runge--Kutta step',stage+1_pInt
#endif
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) &
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), &
                                           crystallite_Fe, &
                                           crystallite_Fi(1:3,1:3,g,i,e), &
                                           crystallite_Fp, &
                                           C(stage)*crystallite_subdt(g,i,e), & ! fraction of original timestep
                                           crystallite_subFrac, g,i,e)
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO PRIVATE(p,cc,NaN)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then

         p = phaseAt(g,i,e)
         cc = phasememberAt(g,i,e)
         NaN = any(IEEE_is_NaN(plasticState(p)%dotState(:,cc)))
         do mySource = 1_pInt, phase_Nsources(p)
           NaN = NaN .or. any(IEEE_is_NaN(sourceState(p)%p(mySource)%dotState(:,cc)))
         enddo
         if (NaN) then                                                                                   ! NaN occured in any dotState
           if (.not. crystallite_localPlasticity(g,i,e)) then                                              ! if broken non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                       ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           else                                                                                            ! if broken local...
             crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
           endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP END PARALLEL

 enddo


!--------------------------------------------------------------------------------------------------
! --- STATE UPDATE WITH ERROR ESTIMATE FOR STATE ---

 relPlasticStateResiduum = 0.0_pReal
 relSourceStateResiduum = 0.0_pReal
 !$OMP PARALLEL
 !$OMP DO PRIVATE(p,cc)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       p = phaseAt(g,i,e)
       cc = phasememberAt(g,i,e)
       plasticState(p)%RKCK45dotState(6,:,cc) = plasticState (p)%dotState(:,cc)                            ! store Runge-Kutta dotState
       do mySource = 1_pInt, phase_Nsources(p)
         sourceState(p)%p(mySource)%RKCK45dotState(6,:,cc) = sourceState(p)%p(mySource)%dotState(:,cc)     ! store Runge-Kutta dotState
       enddo
     endif
   enddo; enddo; enddo
 !$OMP ENDDO

 !$OMP DO PRIVATE(mySizePlasticDotState,mySizeSourceDotState,p,cc)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       p = phaseAt(g,i,e)
       cc = phasememberAt(g,i,e)

       ! --- absolute residuum in state  ---
       mySizePlasticDotState = plasticState(p)%sizeDotState
       plasticStateResiduum(1:mySizePlasticDotState,g,i,e) = &
         matmul(transpose(plasticState(p)%RKCK45dotState(1:6,1:mySizePlasticDotState,cc)),DB) &
       * crystallite_subdt(g,i,e)
       do mySource = 1_pInt, phase_Nsources(p)
         mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
         sourceStateResiduum(1:mySizeSourceDotState,mySource,g,i,e) = &
           matmul(transpose(sourceState(p)%p(mySource)%RKCK45dotState(1:6,1:mySizeSourceDotState,cc)),DB) &
         * crystallite_subdt(g,i,e)
       enddo

       ! --- dot state ---
       plasticState(p)%dotState(:,cc) =  &
         matmul(transpose(plasticState(p)%RKCK45dotState(1:6,1:mySizePlasticDotState,cc)), B)
       do mySource = 1_pInt, phase_Nsources(p)
         mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
         sourceState(p)%p(mySource)%dotState(:,cc)  = &
           matmul(transpose(sourceState(p)%p(mySource)%RKCK45dotState(1:6,1:mySizeSourceDotState,cc)),B)
       enddo
     endif
   enddo; enddo; enddo
 !$OMP ENDDO

 ! --- state and update ---

 !$OMP DO PRIVATE(mySizePlasticDotState,mySizeSourceDotState,p,cc)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then

         p = phaseAt(g,i,e)
         cc = phasememberAt(g,i,e)
         mySizePlasticDotState = plasticState(p)%sizeDotState
         plasticState(p)%state    (1:mySizePlasticDotState,cc) = &
         plasticState(p)%subState0(1:mySizePlasticDotState,cc) &
       + plasticState(p)%dotState (1:mySizePlasticDotState,cc) &
       * crystallite_subdt(g,i,e)
         do mySource = 1_pInt, phase_Nsources(p)
           mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
           sourceState(p)%p(mySource)%state    (1:mySizeSourceDotState,cc) = &
           sourceState(p)%p(mySource)%subState0(1:mySizeSourceDotState,cc) &
         + sourceState(p)%p(mySource)%dotState (1:mySizeSourceDotState,cc)&
         * crystallite_subdt(g,i,e)
         enddo
     endif
   enddo; enddo; enddo
 !$OMP ENDDO

 ! --- relative residui and state convergence ---

 !$OMP DO PRIVATE(mySizePlasticDotState,mySizeSourceDotState,p,cc,s)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       p  = phaseAt(g,i,e)
       cc = phasememberAt(g,i,e)
       mySizePlasticDotState = plasticState(p)%sizeDotState
       forall (s = 1_pInt:mySizePlasticDotState,    abs(plasticState(p)%state(s,cc)) > 0.0_pReal) &
         relPlasticStateResiduum(s,g,i,e) = &
            plasticStateResiduum(s,g,i,e) / plasticState(p)%state(s,cc)

       do mySource = 1_pInt, phase_Nsources(p)
         mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
         forall (s = 1_pInt:mySizeSourceDotState,abs(sourceState(p)%p(mySource)%state(s,cc)) > 0.0_pReal) &
           relSourceStateResiduum(s,mySource,g,i,e) = &
              sourceStateResiduum(s,mySource,g,i,e) / sourceState(p)%p(mySource)%state(s,cc)
       enddo
       !$OMP FLUSH(relPlasticStateResiduum)
       !$OMP FLUSH(relSourceStateResiduum)
! @Martin: do we need flushing? why..?
       crystallite_todo(g,i,e) = all(abs(relPlasticStateResiduum(1:mySizePlasticDotState,g,i,e)) < &
                                     rTol_crystalliteState .or. &
                                     abs(plasticStateResiduum(1:mySizePlasticDotState,g,i,e)) < &
                                     plasticState(p)%aTolState(1:mySizePlasticDotState))
       do mySource = 1_pInt, phase_Nsources(p)
         mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
         crystallite_todo(g,i,e) = crystallite_todo(g,i,e) .and. &
                                   all(abs(relSourceStateResiduum(1:mySizeSourceDotState,mySource,g,i,e)) < &
                                       rTol_crystalliteState .or. &
                                       abs(sourceStateResiduum(1:mySizeSourceDotState,mySource,g,i,e)) < &
                                       sourceState(p)%p(mySource)%aTolState(1:mySizeSourceDotState))
       enddo

#ifdef DEBUG
       if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt&
           .and. ((e == debug_e .and. i == debug_i .and. g == debug_g)&
                  .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
         write(6,'(a,i8,1x,i3,1x,i3,/)') '<< CRYST >> updateState at el ip ipc ',e,i,g
         write(6,'(a,/,(12x,12(f12.1,1x)),/)') '<< CRYST >> absolute residuum tolerance', &
               plasticStateResiduum(1:mySizePlasticDotState,g,i,e) / plasticState(p)%aTolState(1:mySizePlasticDotState)
         write(6,'(a,/,(12x,12(f12.1,1x)),/)') '<< CRYST >> relative residuum tolerance', &
               relPlasticStateResiduum(1:mySizePlasticDotState,g,i,e) / rTol_crystalliteState
         write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> dotState', &
              plasticState(p)%dotState(1:mySizePlasticDotState,cc)
         write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> new state', &
              plasticState(p)%state(1:mySizePlasticDotState,cc)
       endif
#endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO


 ! --- STATE JUMP ---

 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       crystallite_todo(g,i,e) = crystallite_stateJump(g,i,e)
       !$OMP FLUSH(crystallite_todo)
       if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then              ! if broken non-local...
         !$OMP CRITICAL (checkTodo)
           crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                           ! ...all non-locals skipped
         !$OMP END CRITICAL (checkTodo)
       endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO


!--------------------------------------------------------------------------------------------------
! --- UPDATE DEPENDENT STATES IF RESIDUUM BELOW TOLERANCE ---
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) &
       call constitutive_microstructure(crystallite_orientation,       &
                                        crystallite_Fe(1:3,1:3,g,i,e), &
                                        crystallite_Fp(1:3,1:3,g,i,e), &
                                        g, i, e)                                                           ! update dependent state variables to be consistent with basic states
  enddo; enddo; enddo
 !$OMP ENDDO


!--------------------------------------------------------------------------------------------------
! --- FINAL STRESS INTEGRATION STEP IF RESIDUUM BELOW TOLERANCE ---
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       crystallite_todo(g,i,e) = crystallite_integrateStress(g,i,e)
       !$OMP FLUSH(crystallite_todo)
       if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then              ! if broken non-local...
         !$OMP CRITICAL (checkTodo)
           crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                           ! ...all non-locals skipped
         !$OMP END CRITICAL (checkTodo)
       endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO


!--------------------------------------------------------------------------------------------------
! --- SET CONVERGENCE FLAG ---
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     crystallite_converged(g,i,e) = crystallite_todo(g,i,e) .or. crystallite_converged(g,i,e)              ! if still "to do" then converged per definition
   enddo; enddo; enddo
 !$OMP ENDDO

 !$OMP END PARALLEL


 ! --- nonlocal convergence check ---

 if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) &
   write(6,'(a,i8,a,i2,/)') '<< CRYST >> ', count(crystallite_converged(:,:,:)), ' grains converged'    ! if not requesting Integration of just a single IP
 if ((.not. singleRun) .and. any(.not. crystallite_converged .and. .not. crystallite_localPlasticity)) &  ! any non-local not yet converged (or broken)...
   crystallite_converged = crystallite_converged .and. crystallite_localPlasticity                       ! ...restart all non-local as not converged

end subroutine crystallite_integrateStateRKCK45


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with 1st order Euler method with adaptive step size
!--------------------------------------------------------------------------------------------------
subroutine crystallite_integrateStateAdaptiveEuler()
 use, intrinsic :: &
   IEEE_arithmetic
 use debug, only: &
#ifdef DEBUG
   debug_e, &
   debug_i, &
   debug_g, &
#endif
   debug_level, &
   debug_crystallite, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective
 use numerics, only: &
   rTol_crystalliteState
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP
 use mesh, only: &
   mesh_element, &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_Ngrains, &
   plasticState, &
   sourceState, &
   phaseAt, phasememberAt, &
   phase_Nsources, &
   homogenization_maxNgrains
 use constitutive, only: &
   constitutive_collectDotState, &
   constitutive_microstructure, &
   constitutive_plasticity_maxSizeDotState, &
   constitutive_source_maxSizeDotState

 implicit none
 integer(pInt) :: &
   e, &                                                                                             ! element index in element loop
   i, &                                                                                             ! integration point index in ip loop
   g, &                                                                                             ! grain index in grain loop
   s, &                                                                                             ! state index
   p, &
   c, &
   mySource, &
   mySizePlasticDotState, &                                                                         ! size of dot states
   mySizeSourceDotState
 integer(pInt), dimension(2) :: &
   eIter                                                                                            ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) :: &
   iIter, &                                                                                         ! bounds for ip iteration
   gIter                                                                                            ! bounds for grain iteration
 real(pReal), dimension(constitutive_plasticity_maxSizeDotState,            &
                        homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   plasticStateResiduum, &                                                                          ! residuum from evolution in micrstructure
   relPlasticStateResiduum                                                                          ! relative residuum from evolution in microstructure
 real(pReal), dimension(constitutive_source_maxSizeDotState,&
                        maxval(phase_Nsources), &
                        homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   sourceStateResiduum, &                                                                           ! residuum from evolution in micrstructure
   relSourceStateResiduum                                                                           ! relative residuum from evolution in microstructure

 logical :: &
   converged, &
   NaN, &
   singleRun                                                                                        ! flag indicating computation for single (g,i,e) triple


 ! --- LOOP ITERATOR FOR ELEMENT, GRAIN, IP ---
 eIter = FEsolving_execElem(1:2)
 do e = eIter(1),eIter(2)
   iIter(1:2,e) = FEsolving_execIP(1:2,e)
   gIter(1:2,e) = [ 1_pInt,homogenization_Ngrains(mesh_element(3,e))]
 enddo

 singleRun = (eIter(1) == eIter(2) .and. iIter(1,eIter(1)) == iIter(2,eIter(2)))


 plasticStateResiduum = 0.0_pReal
 relPlasticStateResiduum = 0.0_pReal
 sourceStateResiduum = 0.0_pReal
 relSourceStateResiduum = 0.0_pReal


 !$OMP PARALLEL
   ! --- DOT STATE (EULER INTEGRATION) ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) &
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), &
                                           crystallite_Fe, &
                                           crystallite_Fi(1:3,1:3,g,i,e), &
                                           crystallite_Fp, &
                                           crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
    enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO PRIVATE(p,c,NaN)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         p = phaseAt(g,i,e)
         c = phasememberAt(g,i,e)
         NaN = any(IEEE_is_NaN(plasticState(p)%dotState(:,c)))
         do mySource = 1_pInt, phase_Nsources(p)
           NaN = NaN .or. any(IEEE_is_NaN(sourceState(p)%p(mySource)%dotState(:,c)))
         enddo
         if (NaN) then                                                                                       ! NaN occured in any dotState
           if (.not. crystallite_localPlasticity(g,i,e)) then                                                ! if broken non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           else                                                                                              ! if broken local...
             crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
           endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- STATE UPDATE (EULER INTEGRATION) ---

   !$OMP DO PRIVATE(mySizePlasticDotState,mySizeSourceDotState,p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         p = phaseAt(g,i,e)
         c = phasememberAt(g,i,e)
         mySizePlasticDotState = plasticState(p)%sizeDotState
         plasticStateResiduum(1:mySizePlasticDotState,g,i,e) = &
       - 0.5_pReal &
       * plasticState(p)%dotstate(1:mySizePlasticDotState,c) &
       * crystallite_subdt(g,i,e)                                                                            ! contribution to absolute residuum in state
         plasticState(p)%state   (1:mySizePlasticDotState,c) = &
         plasticState(p)%state   (1:mySizePlasticDotState,c) &
       + plasticState(p)%dotstate(1:mySizePlasticDotState,c) &
       * crystallite_subdt(g,i,e)
         do mySource = 1_pInt, phase_Nsources(p)
           mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
           sourceStateResiduum(1:mySizeSourceDotState,mySource,g,i,e) = &
         - 0.5_pReal &
         * sourceState(p)%p(mySource)%dotstate(1:mySizeSourceDotState,c) &
         * crystallite_subdt(g,i,e)                                                                         ! contribution to absolute residuum in state
           sourceState(p)%p(mySource)%state   (1:mySizeSourceDotState,c) = &
           sourceState(p)%p(mySource)%state   (1:mySizeSourceDotState,c) &
         + sourceState(p)%p(mySource)%dotstate(1:mySizeSourceDotState,c) &
         * crystallite_subdt(g,i,e)
         enddo
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- STATE JUMP ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         crystallite_todo(g,i,e) = crystallite_stateJump(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then            ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- UPDATE DEPENDENT STATES (EULER INTEGRATION) ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) &
         call constitutive_microstructure(crystallite_orientation,       &
                                          crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), &
                                          g, i, e)                                                         ! update dependent state variables to be consistent with basic states
     enddo; enddo; enddo
   !$OMP ENDDO
 !$OMP END PARALLEL


 ! --- STRESS INTEGRATION (EULER INTEGRATION) ---

 !$OMP PARALLEL DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       crystallite_todo(g,i,e) = crystallite_integrateStress(g,i,e)
       !$OMP FLUSH(crystallite_todo)
       if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then              ! if broken non-local...
         !$OMP CRITICAL (checkTodo)
           crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                           ! ...all non-locals skipped
         !$OMP END CRITICAL (checkTodo)
       endif
     endif
   enddo; enddo; enddo
 !$OMP END PARALLEL DO

   !$OMP PARALLEL
   ! --- DOT STATE (HEUN METHOD) ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) &
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), &
                                           crystallite_Fe, &
                                           crystallite_Fi(1:3,1:3,g,i,e), &
                                           crystallite_Fp, &
                                           crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO PRIVATE(p,c,NaN)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         p = phaseAt(g,i,e)
         c = phasememberAt(g,i,e)
         NaN = any(IEEE_is_NaN(plasticState(p)%dotState(:,c)))
         do mySource = 1_pInt, phase_Nsources(p)
           NaN = NaN .or. any(IEEE_is_NaN(sourceState(p)%p(mySource)%dotState(:,c)))
         enddo
         if (NaN) then                                                                                     ! NaN occured in any dotState
           if (.not. crystallite_localPlasticity(g,i,e)) then                                              ! if broken non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                       ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           else                                                                                            ! if broken local...
             crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
           endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- ERROR ESTIMATE FOR STATE  (HEUN METHOD) ---

   !$OMP SINGLE
   relPlasticStateResiduum = 0.0_pReal
   relSourceStateResiduum = 0.0_pReal
   !$OMP END SINGLE

   !$OMP DO PRIVATE(mySizePlasticDotState,mySizeSourceDotState,converged,p,c,s)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                   ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         p = phaseAt(g,i,e)
         c = phasememberAt(g,i,e)
         ! --- contribution of heun step to absolute residui ---
         mySizePlasticDotState = plasticState(p)%sizeDotState
         plasticStateResiduum(1:mySizePlasticDotState,g,i,e) = &
         plasticStateResiduum(1:mySizePlasticDotState,g,i,e) &
       + 0.5_pReal * plasticState(p)%dotState(:,c) &
       * crystallite_subdt(g,i,e)                                                                           ! contribution to absolute residuum in state
         do mySource = 1_pInt, phase_Nsources(p)
           mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
           sourceStateResiduum(1:mySizeSourceDotState,mySource,g,i,e) = &
           sourceStateResiduum(1:mySizeSourceDotState,mySource,g,i,e) &
         + 0.5_pReal * sourceState(p)%p(mySource)%dotState(:,c) &
         * crystallite_subdt(g,i,e)                                                                         ! contribution to absolute residuum in state
         enddo
         !$OMP FLUSH(plasticStateResiduum)
         !$OMP FLUSH(sourceStateResiduum)

         ! --- relative residui ---
         forall (s = 1_pInt:mySizePlasticDotState, abs(plasticState(p)%dotState(s,c)) > 0.0_pReal) &
           relPlasticStateResiduum(s,g,i,e) = &
              plasticStateResiduum(s,g,i,e) / plasticState(p)%dotState(s,c)
         do mySource = 1_pInt, phase_Nsources(p)
           mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
           forall (s = 1_pInt:mySizeSourceDotState,abs(sourceState(p)%p(mySource)%dotState(s,c)) > 0.0_pReal) &
             relSourceStateResiduum(s,mySource,g,i,e) = &
                sourceStateResiduum(s,mySource,g,i,e) / sourceState(p)%p(mySource)%dotState(s,c)
         enddo
         !$OMP FLUSH(relPlasticStateResiduum)
         !$OMP FLUSH(relSourceStateResiduum)

#ifdef DEBUG

         if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
             .and. ((e == debug_e .and. i == debug_i .and. g == debug_g)&
                     .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
           write(6,'(a,i8,1x,i2,1x,i3,/)')       '<< CRYST >> updateState at el ip g ',e,i,g
           write(6,'(a,/,(12x,12(f12.1,1x)),/)') '<< CRYST >> absolute residuum tolerance', &
                 plasticStateResiduum(1:mySizePlasticDotState,g,i,e) / plasticState(p)%aTolState(1:mySizePlasticDotState)
           write(6,'(a,/,(12x,12(f12.1,1x)),/)') '<< CRYST >> relative residuum tolerance', &
                 relPlasticStateResiduum(1:mySizePlasticDotState,g,i,e) / rTol_crystalliteState
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> dotState', plasticState(p)%dotState(1:mySizePlasticDotState,c) &
                 - 2.0_pReal * plasticStateResiduum(1:mySizePlasticDotState,g,i,e) / crystallite_subdt(g,i,e)     ! calculate former dotstate from higher order solution and state residuum
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> new state', plasticState(p)%state(1:mySizePlasticDotState,c)
         endif
#endif

         ! --- converged ? ---
         converged = all(abs(relPlasticStateResiduum(1:mySizePlasticDotState,g,i,e)) < &
                         rTol_crystalliteState .or. &
                         abs(plasticStateResiduum(1:mySizePlasticDotState,g,i,e)) < &
                         plasticState(p)%aTolState(1:mySizePlasticDotState))
         do mySource = 1_pInt, phase_Nsources(p)
           mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
           converged = converged .and. &
                       all(abs(relSourceStateResiduum(1:mySizeSourceDotState,mySource,g,i,e)) < &
                       rTol_crystalliteState .or. &
                       abs(sourceStateResiduum(1:mySizeSourceDotState,mySource,g,i,e)) < &
                       sourceState(p)%p(mySource)%aTolState(1:mySizeSourceDotState))
         enddo
         if (converged) then
           crystallite_converged(g,i,e) = .true.                                                             ! ... converged per definitionem
           if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
             !$OMP CRITICAL (distributionState)
             !$OMP END CRITICAL (distributionState)
           endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 !$OMP END PARALLEL


 ! --- NONLOCAL CONVERGENCE CHECK ---

 if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) &
   write(6,'(a,i8,a,i2,/)') '<< CRYST >> ', count(crystallite_converged(:,:,:)), ' grains converged'
 if ((.not. singleRun) .and. any(.not. crystallite_converged .and. .not. crystallite_localPlasticity)) &     ! any non-local not yet converged (or broken)...
    crystallite_converged = crystallite_converged .and. crystallite_localPlasticity                          ! ...restart all non-local as not converged


end subroutine crystallite_integrateStateAdaptiveEuler


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, and state with 1st order explicit Euler method
!--------------------------------------------------------------------------------------------------
subroutine crystallite_integrateStateEuler()
 use, intrinsic :: &
   IEEE_arithmetic
 use debug, only: &
#ifdef DEBUG
   debug_e, &
   debug_i, &
   debug_g, &
#endif
   debug_level, &
   debug_crystallite, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective
 use numerics, only: &
   numerics_timeSyncing
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP
 use mesh, only: &
   mesh_element, &
   mesh_NcpElems
 use material, only: &
   plasticState, &
   sourceState, &
   phaseAt, phasememberAt, &
   phase_Nsources, &
   homogenization_Ngrains
 use constitutive, only: &
   constitutive_collectDotState, &
   constitutive_microstructure

 implicit none

 integer(pInt) :: &
   e, &                                                                                             ! element index in element loop
   i, &                                                                                             ! integration point index in ip loop
   g, &                                                                                             ! grain index in grain loop
   p, &                                                                                             ! phase loop
   c, &
   mySource, &
   mySizePlasticDotState, &
   mySizeSourceDotState
 integer(pInt), dimension(2) :: &
   eIter                                                                                            ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) :: &
   iIter, &                                                                                         ! bounds for ip iteration
   gIter                                                                                            ! bounds for grain iteration
 logical :: &
   NaN, &
   singleRun                                                                                        ! flag indicating computation for single (g,i,e) triple


eIter = FEsolving_execElem(1:2)
 do e = eIter(1),eIter(2)
   iIter(1:2,e) = FEsolving_execIP(1:2,e)
   gIter(1:2,e) = [ 1_pInt,homogenization_Ngrains(mesh_element(3,e))]
 enddo

 singleRun = (eIter(1) == eIter(2) .and. iIter(1,eIter(1)) == iIter(2,eIter(2)))

 !$OMP PARALLEL

   ! --- DOT STATE  ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) &
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), &
                                           crystallite_Fe, &
                                           crystallite_Fi(1:3,1:3,g,i,e), &
                                           crystallite_Fp, &
                                           crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO PRIVATE(p,c,NaN)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         c = phasememberAt(g,i,e)
         p = phaseAt(g,i,e)
         NaN = any(IEEE_is_NaN(plasticState(p)%dotState(:,c)))
         do mySource = 1_pInt, phase_Nsources(p)
           NaN = NaN .or. any(IEEE_is_NaN(sourceState(p)%p(mySource)%dotState(:,c)))
         enddo
         if (NaN) then                                                                                       ! NaN occured in any dotState
           if (.not. crystallite_localPlasticity(g,i,e) .and. .not. numerics_timeSyncing) then               ! if broken non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           else                                                                                              ! if broken local...
             crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
           endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- UPDATE STATE  ---

   !$OMP DO PRIVATE(mySizePlasticDotState,mySizeSourceDotState,p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         p = phaseAt(g,i,e)
         c = phasememberAt(g,i,e)
         mySizePlasticDotState = plasticState(p)%sizeDotState
         plasticState(p)%state(   1:mySizePlasticDotState,c) = &
         plasticState(p)%state(   1:mySizePlasticDotState,c) &
       + plasticState(p)%dotState(1:mySizePlasticDotState,c) &
       * crystallite_subdt(g,i,e)
         do mySource = 1_pInt, phase_Nsources(p)
           mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
           sourceState(p)%p(mySource)%state(   1:mySizeSourceDotState,c) = &
           sourceState(p)%p(mySource)%state(   1:mySizeSourceDotState,c) &
         + sourceState(p)%p(mySource)%dotState(1:mySizeSourceDotState,c) &
         * crystallite_subdt(g,i,e)
         enddo

#ifdef DEBUG
         if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
             .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                     .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
           p = phaseAt(g,i,e)
           c = phasememberAt(g,i,e)
           write(6,'(a,i8,1x,i2,1x,i3,/)')       '<< CRYST >> update state at el ip g ',e,i,g
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> dotState',  plasticState(p)%dotState(1:mySizePlasticDotState,c)
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> new state', plasticState(p)%state   (1:mySizePlasticDotState,c)
         endif
#endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- STATE JUMP ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                   ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         crystallite_todo(g,i,e) = crystallite_stateJump(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e) &                 ! if broken non-local...
             .and. .not. numerics_timeSyncing) then
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                          ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- UPDATE DEPENDENT STATES ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) &
         call constitutive_microstructure(crystallite_orientation,       &
                                          crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), &
                                          g, i, e)                                                           ! update dependent state variables to be consistent with basic states
   enddo; enddo; enddo
   !$OMP ENDDO
  !$OMP END PARALLEL


 !$OMP PARALLEL
 ! --- STRESS INTEGRATION ---

 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
       crystallite_todo(g,i,e) = crystallite_integrateStress(g,i,e)
       !$OMP FLUSH(crystallite_todo)
       if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e) &                  ! if broken non-local...
           .and. .not. numerics_timeSyncing) then
         !$OMP CRITICAL (checkTodo)
           crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                           ! ...all non-locals skipped
         !$OMP END CRITICAL (checkTodo)
       endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO


 ! --- SET CONVERGENCE FLAG ---

 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     crystallite_converged(g,i,e) = crystallite_todo(g,i,e) .or. crystallite_converged(g,i,e)              ! if still "to do" then converged per definitionem
   enddo; enddo; enddo
 !$OMP ENDDO

 !$OMP END PARALLEL


 ! --- CHECK NON-LOCAL CONVERGENCE ---

 if (.not. singleRun) then                                                                                 ! if not requesting Integration of just a single IP
   if (any(.not. crystallite_converged .and. .not. crystallite_localPlasticity) &                          ! any non-local not yet converged (or broken)...
       .and. .not. numerics_timeSyncing) &
     crystallite_converged = crystallite_converged .and. crystallite_localPlasticity                       ! ...restart all non-local as not converged
 endif

end subroutine crystallite_integrateStateEuler


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with adaptive 1st order explicit Euler method
!> using Fixed Point Iteration to adapt the stepsize
!--------------------------------------------------------------------------------------------------
subroutine crystallite_integrateStateFPI()
 use, intrinsic :: &
   IEEE_arithmetic
 use debug, only: &
#ifdef DEBUG
   debug_e, &
   debug_i, &
   debug_g, &
#endif
   debug_level,&
   debug_crystallite, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective
 use numerics, only: &
   nState, &
   rTol_crystalliteState
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP
 use mesh, only: &
   mesh_element, &
   mesh_NcpElems
 use material, only: &
   plasticState, &
   sourceState, &
   phaseAt, phasememberAt, &
   phase_Nsources, &
   homogenization_Ngrains
 use constitutive, only: &
   constitutive_collectDotState, &
   constitutive_microstructure, &
   constitutive_plasticity_maxSizeDotState, &
   constitutive_source_maxSizeDotState

 implicit none

 integer(pInt) :: &
   NiterationState, &                                                                               !< number of iterations in state loop
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g, &                                                                                             !< grain index in grain loop
   p, &
   c, &
   mySource, &
   mySizePlasticDotState, &                                                                         ! size of dot states
   mySizeSourceDotState
 integer(pInt), dimension(2) :: &
   eIter                                                                                            ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) :: &
   iIter, &                                                                                         ! bounds for ip iteration
   gIter                                                                                            ! bounds for grain iteration
 real(pReal) :: &
   dot_prod12, &
   dot_prod22, &
   plasticStateDamper, &                                                                            ! damper for integration of state
   sourceStateDamper
 real(pReal), dimension(constitutive_plasticity_maxSizeDotState) :: &
   plasticStateResiduum, &
   tempPlasticState
 real(pReal), dimension(constitutive_source_maxSizeDotState, maxval(phase_Nsources)) :: &
   sourceStateResiduum, &                                                                           ! residuum from evolution in micrstructure
   tempSourceState
 logical :: &
   converged, &
   NaN, &
   singleRun, &                                                                                     ! flag indicating computation for single (g,i,e) triple
   doneWithIntegration

 eIter = FEsolving_execElem(1:2)
 do e = eIter(1),eIter(2)
   iIter(1:2,e) = FEsolving_execIP(1:2,e)
   gIter(1:2,e) = [ 1_pInt,homogenization_Ngrains(mesh_element(3,e))]
 enddo

 singleRun = (eIter(1) == eIter(2) .and. iIter(1,eIter(1)) == iIter(2,eIter(2)))

 if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) &
   write(6,'(a,i8,a)') '<< CRYST >> ', count(crystallite_todo(:,:,:)),' grains todo at start of state integration'

!--------------------------------------------------------------------------------------------------
! initialize dotState
 if (.not. singleRun) then
   forall(p = 1_pInt:size(plasticState))
     plasticState(p)%previousDotState  = 0.0_pReal
     plasticState(p)%previousDotState2 = 0.0_pReal
   end forall
   do p = 1_pInt, size(sourceState); do mySource = 1_pInt, phase_Nsources(p)
     sourceState(p)%p(mySource)%previousDotState  = 0.0_pReal
     sourceState(p)%p(mySource)%previousDotState2 = 0.0_pReal
   enddo; enddo
 else
   e = eIter(1)
   i = iIter(1,e)
   do g = gIter(1,e), gIter(2,e)
     p = phaseAt(g,i,e)
     c = phasememberAt(g,i,e)
     plasticState(p)%previousDotState (:,c) = 0.0_pReal
     plasticState(p)%previousDotState2(:,c) = 0.0_pReal
     do mySource = 1_pInt, phase_Nsources(p)
       sourceState(p)%p(mySource)%previousDotState (:,c) = 0.0_pReal
       sourceState(p)%p(mySource)%previousDotState2(:,c) = 0.0_pReal
     enddo
   enddo
 endif

 ! --+>> PREGUESS FOR STATE <<+--

 ! --- DOT STATES ---

 !$OMP PARALLEL
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)             ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) &
       call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), &
                                         crystallite_Fe, &
                                         crystallite_Fi(1:3,1:3,g,i,e), &
                                         crystallite_Fp, &
                                         crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
   enddo; enddo; enddo

 !$OMP ENDDO
 !$OMP DO PRIVATE(p,c,NaN)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       p = phaseAt(g,i,e)
       c = phasememberAt(g,i,e)
       NaN = any(IEEE_is_NaN(plasticState(p)%dotState(:,c)))
       do mySource = 1_pInt, phase_Nsources(p)
         NaN = NaN .or. any(IEEE_is_NaN(sourceState(p)%p(mySource)%dotState(:,c)))
       enddo
       if (NaN) then                                                                                       ! NaN occured in any dotState
         if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) &
           write(6,*) '<< CRYST >> dotstate ',plasticState(p)%dotState(:,c)
         if (.not. crystallite_localPlasticity(g,i,e)) then                                                ! if broken is a non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals done (and broken)
           !$OMP END CRITICAL (checkTodo)
         else                                                                                              ! broken one was local...
           crystallite_todo(g,i,e) = .false.                                                               ! ... done (and broken)
         endif
       endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO

 ! --- UPDATE STATE  ---
 if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) &
   write(6,'(a,i8,a)') '<< CRYST >> ', count(crystallite_todo(:,:,:)),' grains todo after preguess of state'


 !$OMP DO PRIVATE(mySizePlasticDotState,mySizeSourceDotState,p,c)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       p = phaseAt(g,i,e)
       c = phasememberAt(g,i,e)
       mySizePlasticDotState = plasticState(p)%sizeDotState
       plasticState(p)%state(1:mySizePlasticDotState,c) = &
         plasticState(p)%subState0(1:mySizePlasticDotState,c) &
       + plasticState(p)%dotState (1:mySizePlasticDotState,c) &
       * crystallite_subdt(g,i,e)
       do mySource = 1_pInt, phase_Nsources(p)
         mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
         sourceState(p)%p(mySource)%state(1:mySizeSourceDotState,c) = &
           sourceState(p)%p(mySource)%subState0(1:mySizeSourceDotState,c) &
         + sourceState(p)%p(mySource)%dotState (1:mySizeSourceDotState,c) &
         * crystallite_subdt(g,i,e)
       enddo
     endif
   enddo; enddo; enddo
 !$OMP ENDDO

 !$OMP END PARALLEL

 ! --+>> STATE LOOP <<+--

 NiterationState = 0_pInt
 doneWithIntegration = .false.
 crystalliteLooping: do while (.not. doneWithIntegration .and. NiterationState < nState)
   NiterationState = NiterationState + 1_pInt

   !$OMP PARALLEL

   ! --- UPDATE DEPENDENT STATES ---

   !$OMP DO PRIVATE(p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) &
         call constitutive_microstructure(crystallite_orientation,       &
                                          crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), &
                                          g, i, e)                                                           ! update dependent state variables to be consistent with basic states
       p = phaseAt(g,i,e)
       c = phasememberAt(g,i,e)
       plasticState(p)%previousDotState2(:,c) = plasticState(p)%previousDotState(:,c)
       plasticState(p)%previousDotState (:,c) = plasticState(p)%dotState(:,c)
       do mySource = 1_pInt, phase_Nsources(p)
         sourceState(p)%p(mySource)%previousDotState2(:,c) = sourceState(p)%p(mySource)%previousDotState(:,c)
         sourceState(p)%p(mySource)%previousDotState (:,c) = sourceState(p)%p(mySource)%dotState(:,c)
       enddo
     enddo; enddo; enddo
   !$OMP ENDDO

   ! --- STRESS INTEGRATION ---

   if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) &
     write(6,'(a,i8,a)') '<< CRYST >> ', count(crystallite_todo(:,:,:)),' grains todo before stress integration'

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         crystallite_todo(g,i,e) = crystallite_integrateStress(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then            ! broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ... then all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO

   !$OMP SINGLE
   !$OMP CRITICAL (write2out)
   if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) &
     write(6,'(a,i8,a)') '<< CRYST >> ', count(crystallite_todo(:,:,:)),' grains todo after stress integration'
   !$OMP END CRITICAL (write2out)
   !$OMP END SINGLE


   ! --- DOT STATE  ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) &
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), &
                                           crystallite_Fe, &
                                           crystallite_Fi(1:3,1:3,g,i,e), &
                                           crystallite_Fp, &
                                           crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
     enddo; enddo; enddo
   !$OMP ENDDO

   !$OMP DO PRIVATE(p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         p = phaseAt(g,i,e)
         c = phasememberAt(g,i,e)
         NaN = any(IEEE_is_NaN(plasticState(p)%dotState(:,c)))
         do mySource = 1_pInt, phase_Nsources(p)
           NaN = NaN .or. any(IEEE_is_NaN(sourceState(p)%p(mySource)%dotState(:,c)))
         enddo
         if (NaN) then                                                                                       ! NaN occured in any dotState
           crystallite_todo(g,i,e) = .false.                                                                 ! ... skip me next time
           if (.not. crystallite_localPlasticity(g,i,e)) then                                                ! if me is non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           endif
         endif

       endif

     enddo; enddo; enddo
   !$OMP ENDDO

   ! --- UPDATE STATE  ---

   !$OMP DO PRIVATE(dot_prod12,dot_prod22, &
   !$OMP&           mySizePlasticDotState,mySizeSourceDotState, &
   !$OMP&           plasticStateResiduum,sourceStateResiduum, &
   !$OMP&           plasticStatedamper,sourceStateDamper, &
   !$OMP&           tempPlasticState,tempSourceState,converged,p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)           ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then

         p = phaseAt(g,i,e)
         c = phasememberAt(g,i,e)
         dot_prod12 = dot_product(  plasticState(p)%dotState         (:,c) &
                                  - plasticState(p)%previousDotState (:,c), &
                                    plasticState(p)%previousDotState (:,c) &
                                  - plasticState(p)%previousDotState2(:,c))
         dot_prod22 = dot_product(  plasticState(p)%previousDotState (:,c) &
                                  - plasticState(p)%previousDotState2(:,c), &
                                    plasticState(p)%previousDotState (:,c) &
                                  - plasticState(p)%previousDotState2(:,c))
         if (      dot_prod22 > 0.0_pReal &
             .and. (     dot_prod12 < 0.0_pReal &
                    .or. dot_product(plasticState(p)%dotState(:,c), &
                                     plasticState(p)%previousDotState(:,c)) < 0.0_pReal) ) then
           plasticStateDamper = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)
         else
           plasticStateDamper = 1.0_pReal
         endif
         ! --- get residui ---

         mySizePlasticDotState = plasticState(p)%sizeDotState
         plasticStateResiduum(1:mySizePlasticDotState) = &
           plasticState(p)%state(1:mySizePlasticDotState,c)      &
         - plasticState(p)%subState0(1:mySizePlasticDotState,c)  &
         - (  plasticState(p)%dotState(1:mySizePlasticDotState,c) * plasticStateDamper &
            + plasticState(p)%previousDotState(1:mySizePlasticDotState,c) &
            * (1.0_pReal - plasticStateDamper)) * crystallite_subdt(g,i,e)

         ! --- correct state with residuum ---
         tempPlasticState(1:mySizePlasticDotState) = &
           plasticState(p)%state(1:mySizePlasticDotState,c) &
         - plasticStateResiduum(1:mySizePlasticDotState)                              ! need to copy to local variable, since we cant flush a pointer in openmp

         ! --- store corrected dotState --- (cannot do this before state update, because not sure how to flush pointers in openmp)

         plasticState(p)%dotState(:,c) = plasticState(p)%dotState(:,c) * plasticStateDamper &
                                       + plasticState(p)%previousDotState(:,c) &
                                       * (1.0_pReal - plasticStateDamper)

         do mySource = 1_pInt, phase_Nsources(p)
           mySizeSourceDotState  = sourceState(p)%p(mySource)%sizeDotState
           dot_prod12 = dot_product(  sourceState(p)%p(mySource)%dotState         (:,c) &
                                    - sourceState(p)%p(mySource)%previousDotState (:,c), &
                                      sourceState(p)%p(mySource)%previousDotState (:,c) &
                                    - sourceState(p)%p(mySource)%previousDotState2(:,c))
           dot_prod22 = dot_product(  sourceState(p)%p(mySource)%previousDotState (:,c) &
                                    - sourceState(p)%p(mySource)%previousDotState2(:,c), &
                                      sourceState(p)%p(mySource)%previousDotState (:,c) &
                                    - sourceState(p)%p(mySource)%previousDotState2(:,c))

           if (      dot_prod22 > 0.0_pReal &
               .and. (     dot_prod12 < 0.0_pReal &
                      .or. dot_product(sourceState(p)%p(mySource)%dotState(:,c), &
                                       sourceState(p)%p(mySource)%previousDotState(:,c)) < 0.0_pReal) ) then
             sourceStateDamper = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)
           else
             sourceStateDamper = 1.0_pReal
           endif
         ! --- get residui ---
           mySizeSourceDotState  = sourceState(p)%p(mySource)%sizeDotState
           sourceStateResiduum(1:mySizeSourceDotState,mySource) = &
             sourceState(p)%p(mySource)%state(1:mySizeSourceDotState,c)      &
           - sourceState(p)%p(mySource)%subState0(1:mySizeSourceDotState,c)  &
           - (  sourceState(p)%p(mySource)%dotState(1:mySizeSourceDotState,c) * sourceStateDamper &
              + sourceState(p)%p(mySource)%previousDotState(1:mySizeSourceDotState,c) &
              * (1.0_pReal - sourceStateDamper)) * crystallite_subdt(g,i,e)

         ! --- correct state with residuum ---
           tempSourceState(1:mySizeSourceDotState,mySource) = &
             sourceState(p)%p(mySource)%state(1:mySizeSourceDotState,c) &
           - sourceStateResiduum(1:mySizeSourceDotState,mySource)                     ! need to copy to local variable, since we cant flush a pointer in openmp

         ! --- store corrected dotState --- (cannot do this before state update, because not sure how to flush pointers in openmp)
           sourceState(p)%p(mySource)%dotState(:,c) = &
             sourceState(p)%p(mySource)%dotState(:,c) * sourceStateDamper &
           + sourceState(p)%p(mySource)%previousDotState(:,c) &
           * (1.0_pReal - sourceStateDamper)
         enddo

#ifdef DEBUG
         if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
             .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                    .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
           write(6,'(a,i8,1x,i2,1x,i3,/)')       '<< CRYST >> update state at el ip g ',e,i,g
           write(6,'(a,f6.1,/)')                 '<< CRYST >> plasticstatedamper ',plasticStatedamper
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> plastic state residuum',&
                                                  abs(plasticStateResiduum(1:mySizePlasticDotState))
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> abstol dotstate',plasticState(p)%aTolState(1:mySizePlasticDotState)
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> reltol dotstate',rTol_crystalliteState* &
                                                  abs(tempPlasticState(1:mySizePlasticDotState))
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> new state',tempPlasticState(1:mySizePlasticDotState)
         endif
#endif

         ! --- converged ? ---
         converged = all(    abs(plasticStateResiduum(1:mySizePlasticDotState)) < &
                             plasticState(p)%aTolState(1:mySizePlasticDotState) &
                        .or. abs(plasticStateResiduum(1:mySizePlasticDotState)) < &
                             rTol_crystalliteState * abs(tempPlasticState(1:mySizePlasticDotState)))
         do mySource = 1_pInt, phase_Nsources(p)
           mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
           converged = converged .and. &
                       all(    abs(sourceStateResiduum(1:mySizeSourceDotState,mySource)) < &
                               sourceState(p)%p(mySource)%aTolState(1:mySizeSourceDotState) &
                          .or. abs(sourceStateResiduum(1:mySizeSourceDotState,mySource)) < &
                               rTol_crystalliteState * abs(tempSourceState(1:mySizeSourceDotState,mySource)))
         enddo
         if (converged) crystallite_converged(g,i,e) = .true.                                                                   ! ... converged per definition

         plasticState(p)%state(1:mySizePlasticDotState,c) = &
           tempPlasticState(1:mySizePlasticDotState)
         do mySource = 1_pInt, phase_Nsources(p)
           mySizeSourceDotState = sourceState(p)%p(mySource)%sizeDotState
           sourceState(p)%p(mySource)%state(1:mySizeSourceDotState,c) = &
             tempSourceState(1:mySizeSourceDotState,mySource)
         enddo
        endif
     enddo; enddo; enddo
   !$OMP ENDDO
   ! --- STATE JUMP ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. crystallite_converged(g,i,e)) then                                  ! converged and still alive...
         crystallite_todo(g,i,e) = crystallite_stateJump(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e)) then                                                             ! if state jump fails, then convergence is broken
           crystallite_converged(g,i,e) = .false.
           if (.not. crystallite_localPlasticity(g,i,e)) then                                                ! if broken non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP END PARALLEL

   if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) &
     write(6,'(a,i8,a,i2)') '<< CRYST >> ', count(crystallite_converged(:,:,:)), &
                            ' grains converged after state integration #', NiterationState


   ! --- NON-LOCAL CONVERGENCE CHECK ---

   if (.not. singleRun) then                                                                               ! if not requesting Integration of just a single IP
     if (any(.not. crystallite_converged .and. .not. crystallite_localPlasticity)) &                       ! any non-local not yet converged (or broken)...
       crystallite_converged = crystallite_converged .and. crystallite_localPlasticity                     ! ...restart all non-local as not converged
   endif

   if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) then
     write(6,'(a,i8,a)')      '<< CRYST >> ', count(crystallite_converged(:,:,:)), &
                              ' grains converged after non-local check'
     write(6,'(a,i8,a,i2,/)') '<< CRYST >> ', count(crystallite_todo(:,:,:)), &
                              ' grains todo after state integration #', NiterationState
   endif
   ! --- CHECK IF DONE WITH INTEGRATION ---

   doneWithIntegration = .true.
   elemLoop: do e = eIter(1),eIter(2)
     do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         doneWithIntegration = .false.
         exit elemLoop
       endif
     enddo; enddo
   enddo elemLoop

 enddo crystalliteLooping
end subroutine crystallite_integrateStateFPI


!--------------------------------------------------------------------------------------------------
!> @brief calculates a jump in the state according to the current state and the current stress
!> returns true, if state jump was successfull or not needed. false indicates NaN in delta state
!--------------------------------------------------------------------------------------------------
logical function crystallite_stateJump(ipc,ip,el)
 use, intrinsic :: &
   IEEE_arithmetic
 use prec, only: &
   dNeq0
#ifdef DEBUG
 use debug, only: &
   debug_e, &
   debug_i, &
   debug_g, &
   debug_level, &
   debug_crystallite, &
   debug_levelExtensive, &
   debug_levelSelective
#endif
 use material, only: &
   plasticState, &
   sourceState, &
   phase_Nsources, &
   phaseAt, phasememberAt
 use constitutive, only: &
   constitutive_collectDeltaState

 implicit none
 integer(pInt), intent(in):: &
   el, &                       ! element index
   ip, &                       ! integration point index
   ipc                         ! grain index

 integer(pInt) :: &
   c, &
   p, &
   mySource, &
   myOffsetPlasticDeltaState, &
   myOffsetSourceDeltaState, &
   mySizePlasticDeltaState, &
   mySizeSourceDeltaState

 c = phasememberAt(ipc,ip,el)
 p = phaseAt(ipc,ip,el)

 call constitutive_collectDeltaState(crystallite_Tstar_v(1:6,ipc,ip,el), &
                                     crystallite_Fe(1:3,1:3,ipc,ip,el), &
                                     crystallite_Fi(1:3,1:3,ipc,ip,el), &
                                     ipc,ip,el)

 myOffsetPlasticDeltaState = plasticState(p)%offsetDeltaState
 mySizePlasticDeltaState   = plasticState(p)%sizeDeltaState

 if( any(IEEE_is_NaN(plasticState(p)%deltaState(1:mySizePlasticDeltaState,c)))) then                                       ! NaN occured in deltaState
   crystallite_stateJump = .false.
   return
 endif

 plasticState(p)%state(myOffsetPlasticDeltaState + 1_pInt                 : &
                       myOffsetPlasticDeltaState + mySizePlasticDeltaState,c) = &
 plasticState(p)%state(myOffsetPlasticDeltaState + 1_pInt                 : &
                       myOffsetPlasticDeltaState + mySizePlasticDeltaState,c) + &
    plasticState(p)%deltaState(1:mySizePlasticDeltaState,c)

 do mySource = 1_pInt, phase_Nsources(p)
   myOffsetSourceDeltaState = sourceState(p)%p(mySource)%offsetDeltaState
   mySizeSourceDeltaState   = sourceState(p)%p(mySource)%sizeDeltaState
   if (any(IEEE_is_NaN(sourceState(p)%p(mySource)%deltaState(1:mySizeSourceDeltaState,c)))) then   ! NaN occured in deltaState
     crystallite_stateJump = .false.
     return
   endif
   sourceState(p)%p(mySource)%state(myOffsetSourceDeltaState + 1_pInt                : &
                                    myOffsetSourceDeltaState + mySizeSourceDeltaState,c) = &
   sourceState(p)%p(mySource)%state(myOffsetSourceDeltaState + 1_pInt                : &
                                    myOffsetSourceDeltaState + mySizeSourceDeltaState,c) + &
     sourceState(p)%p(mySource)%deltaState(1:mySizeSourceDeltaState,c)
 enddo

#ifdef DEBUG
 if (any(dNeq0(plasticState(p)%deltaState(1:mySizePlasticDeltaState,c))) &
     .and. iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
     .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
             .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
   write(6,'(a,i8,1x,i2,1x,i3, /)') '<< CRYST >> update state at el ip ipc ',el,ip,ipc
   write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> deltaState', plasticState(p)%deltaState(1:mySizePlasticDeltaState,c)
   write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> new state', &
     plasticState(p)%state(myOffsetPlasticDeltaState + 1_pInt                : &
                           myOffsetPlasticDeltaState + mySizePlasticDeltaState,c)
 endif
#endif

 crystallite_stateJump = .true.

end function crystallite_stateJump


!--------------------------------------------------------------------------------------------------
!> @brief Map 2nd order tensor to reference config
!--------------------------------------------------------------------------------------------------
function crystallite_push33ToRef(ipc,ip,el, tensor33)
 use math, only: &
  math_mul33x33, &
  math_inv33, &
  math_EulerToR
 use material, only: &
  material_EulerAngles

 implicit none
 real(pReal), dimension(3,3) :: crystallite_push33ToRef
 real(pReal), dimension(3,3), intent(in) :: tensor33
 real(pReal), dimension(3,3)             :: T
 integer(pInt), intent(in):: &
   el, &                      ! element index
   ip, &                      ! integration point index
   ipc                         ! grain index

 T = math_mul33x33(math_EulerToR(material_EulerAngles(1:3,ipc,ip,el)), &
                   transpose(math_inv33(crystallite_subF(1:3,1:3,ipc,ip,el))))
 crystallite_push33ToRef = math_mul33x33(transpose(T),math_mul33x33(tensor33,T))

end function crystallite_push33ToRef


!--------------------------------------------------------------------------------------------------
!> @brief calculation of stress (P) with time integration based on a residuum in Lp and
!> intermediate acceleration of the Newton-Raphson correction
!--------------------------------------------------------------------------------------------------
logical function crystallite_integrateStress(&
      ipc,&                                                                                       ! grain number
      ip,&                                                                                        ! integration point number
      el,&                                                                                        ! element number
      timeFraction &
      )
 use, intrinsic :: &
   IEEE_arithmetic
 use prec, only:         pLongInt, &
                         tol_math_check, &
                         dEq0
 use numerics, only:     nStress, &
                         aTol_crystalliteStress, &
                         rTol_crystalliteStress, &
                         iJacoLpresiduum, &
                         subStepSizeLp, &
                         subStepSizeLi
 use debug, only:        debug_level, &
#ifdef DEBUG
                         debug_e, &
                         debug_i, &
                         debug_g, &
#endif
                         debug_crystallite, &
                         debug_levelBasic, &
                         debug_levelExtensive, &
                         debug_levelSelective

 use constitutive, only: constitutive_LpAndItsTangents, &
                         constitutive_LiAndItsTangents, &
                         constitutive_SandItsTangents
 use math, only:         math_mul33x33, &
                         math_mul33xx33, &
                         math_mul3333xx3333, &
                         math_mul66x6, &
                         math_mul99x99, &
                         math_inv33, &
                         math_invert, &
                         math_det33, &
                         math_I3, &
                         math_identity2nd, &
                         math_Mandel66to3333, &
                         math_Mandel6to33, &
                         math_Mandel33to6, &
                         math_Plain3333to99, &
                         math_Plain33to9, &
                         math_Plain9to33, &
                         math_Plain99to3333
#ifdef DEBUG
 use mesh, only:         mesh_element
#endif

 implicit none
 integer(pInt), intent(in)::         el, &                                                           ! element index
                                     ip, &                                                           ! integration point index
                                     ipc                                                             ! grain index
 real(pReal), optional, intent(in) :: timeFraction                                                   ! fraction of timestep

 !*** local variables ***!
 real(pReal), dimension(3,3)::       Fg_new, &                                                       ! deformation gradient at end of timestep
                                     Fp_current, &                                                   ! plastic deformation gradient at start of timestep
                                     Fi_current, &                                                   ! intermediate deformation gradient at start of timestep
                                     Fp_new, &                                                       ! plastic deformation gradient at end of timestep
                                     Fe_new, &                                                       ! elastic deformation gradient at end of timestep
                                     invFp_new, &                                                    ! inverse of Fp_new
                                     Fi_new, &                                                       ! gradient of intermediate deformation stages
                                     invFi_new, &
                                     invFp_current, &                                                ! inverse of Fp_current
                                     invFi_current, &                                                ! inverse of Fp_current
                                     Lpguess, &                                                      ! current guess for plastic velocity gradient
                                     Lpguess_old, &                                                  ! known last good guess for plastic velocity gradient
                                     Lp_constitutive, &                                              ! plastic velocity gradient resulting from constitutive law
                                     residuumLp, &                                                   ! current residuum of plastic velocity gradient
                                     residuumLp_old, &                                               ! last residuum of plastic velocity gradient
                                     deltaLp, &                                                      ! direction of next guess
                                     Liguess, &                                                      ! current guess for intermediate velocity gradient
                                     Liguess_old, &                                                  ! known last good guess for intermediate velocity gradient
                                     Li_constitutive, &                                              ! intermediate velocity gradient resulting from constitutive law
                                     residuumLi, &                                                   ! current residuum of intermediate velocity gradient
                                     residuumLi_old, &                                               ! last residuum of intermediate velocity gradient
                                     deltaLi, &                                                      ! direction of next guess
                                     Tstar, &                                                        ! 2nd Piola-Kirchhoff Stress in plastic (lattice) configuration
                                     A, &
                                     B, &
                                     Fe, &                                                           ! elastic deformation gradient
                                     temp_33
 real(pReal), dimension(6)::         Tstar_v                                                         ! 2nd Piola-Kirchhoff Stress in Mandel-Notation
 real(pReal), dimension(9)::         work                                                            ! needed for matrix inversion by LAPACK
 integer(pInt), dimension(9) ::      ipiv                                                            ! needed for matrix inversion by LAPACK
 real(pReal), dimension(9,9) ::      dRLp_dLp, &                                                     ! partial derivative of residuum (Jacobian for NEwton-Raphson scheme)
                                     dRLp_dLp2, &                                                    ! working copy of dRdLp
                                     dRLi_dLi                                                        ! partial derivative of residuumI (Jacobian for NEwton-Raphson scheme)
 real(pReal), dimension(3,3,3,3)::   dS_dFe, &                                                   ! partial derivative of 2nd Piola-Kirchhoff stress
                                     dS_dFi, &
                                     dFe_dLp, &                                                  ! partial derivative of elastic deformation gradient
                                     dFe_dLi, &
                                     dFi_dLi, &
                                     dLp_dFi, &
                                     dLi_dFi, &
                                     dLp_dS, &
                                     dLi_dS
 real(pReal)                         detInvFi, &                                                     ! determinant of InvFi
                                     steplengthLp, &
                                     steplengthLi, &
                                     dt, &                                                           ! time increment
                                     aTolLp, &
                                     aTolLi
 integer(pInt)                       NiterationStressLp, &                                           ! number of stress integrations
                                     NiterationStressLi, &                                           ! number of inner stress integrations
                                     ierr, &                                                         ! error indicator for LAPACK
                                     o, &
                                     p, &
                                     jacoCounterLp, &
                                     jacoCounterLi                                                    ! counters to check for Jacobian update
 external :: &
   dgesv

 !* be pessimistic
 crystallite_integrateStress = .false.
#ifdef DEBUG
 if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
     .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
            .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) &
 write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> integrateStress at el ip ipc ',el,ip,ipc
#endif

 !* only integrate over fraction of timestep?

 if (present(timeFraction)) then
   dt = crystallite_subdt(ipc,ip,el) * timeFraction
   Fg_new = crystallite_subF0(1:3,1:3,ipc,ip,el) &
          + (crystallite_subF(1:3,1:3,ipc,ip,el) - crystallite_subF0(1:3,1:3,ipc,ip,el)) * timeFraction
 else
   dt = crystallite_subdt(ipc,ip,el)
   Fg_new = crystallite_subF(1:3,1:3,ipc,ip,el)
 endif


 !* feed local variables

 Fp_current  =   crystallite_subFp0(1:3,1:3,ipc,ip,el)                                                  ! "Fp_current" is only used as temp var here...
 Lpguess     =   crystallite_Lp    (1:3,1:3,ipc,ip,el)                                                  ! ... and take it as first guess
 Fi_current  =   crystallite_subFi0(1:3,1:3,ipc,ip,el)                                                  ! intermediate configuration, assume decomposition as F = Fe Fi Fp
 Liguess     =   crystallite_Li    (1:3,1:3,ipc,ip,el)                                                  ! ... and take it as first guess
 Liguess_old =   Liguess


 !* inversion of Fp_current...

 invFp_current = math_inv33(Fp_current)
 failedInversionFp: if (all(dEq0(invFp_current))) then
#ifdef DEBUG
   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3)') '<< CRYST >> integrateStress failed on inversion of Fp_current at el (elFE) ip ipc ',&
       el,'(',mesh_element(1,el),')',ip,ipc
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) > 0_pInt) &
       write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fp_current',transpose(Fp_current(1:3,1:3))
   endif
#endif
   return
 endif failedInversionFp
 A = math_mul33x33(Fg_new,invFp_current)                                                            ! intermediate tensor needed later to calculate dFe_dLp

 !* inversion of Fi_current...

 invFi_current = math_inv33(Fi_current)
 failedInversionFi: if (all(dEq0(invFi_current))) then
#ifdef DEBUG
   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3)') '<< CRYST >> integrateStress failed on inversion of Fi_current at el (elFE) ip ipc ',&
       el,'(',mesh_element(1,el),')',ip,ipc
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) > 0_pInt) &
       write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fp_current',transpose(Fi_current(1:3,1:3))
   endif
#endif
   return
 endif failedInversionFi

 !* start LpLoop with normal step length

 NiterationStressLi = 0_pInt
 jacoCounterLi      = 0_pInt
 steplengthLi       = 1.0_pReal
 residuumLi_old     = 0.0_pReal

 LiLoop: do
   NiterationStressLi = NiterationStressLi + 1_pInt
   IloopsExeced: if (NiterationStressLi > nStress) then
#ifdef DEBUG
     if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) &
       write(6,'(a,i3,a,i8,1x,a,i8,a,1x,i2,1x,i3,/)') '<< CRYST >> integrateStress reached inelastic loop limit',nStress, &
       ' at el (elFE) ip ipc ', el,'(',mesh_element(1,el),')',ip,ipc
#endif
     return
   endif IloopsExeced

   invFi_new = math_mul33x33(invFi_current,math_I3 - dt*Liguess)
   Fi_new    = math_inv33(invFi_new)
   detInvFi  = math_det33(invFi_new)

   NiterationStressLp = 0_pInt
   jacoCounterLp      = 0_pInt
   steplengthLp       = 1.0_pReal
   residuumLp_old     = 0.0_pReal
   Lpguess_old        = Lpguess

   LpLoop: do                                    ! inner stress integration loop for consistency with Fi
     NiterationStressLp = NiterationStressLp + 1_pInt
     loopsExeced: if (NiterationStressLp > nStress) then
#ifdef DEBUG
     if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) &
         write(6,'(a,i3,a,i8,1x,a,i8,a,1x,i2,1x,i3,/)') '<< CRYST >> integrateStress reached loop limit',nStress, &
         ' at el (elFE) ip ipc ', el,'(',mesh_element(1,el),')',ip,ipc
#endif
       return
     endif loopsExeced

     !* calculate (elastic) 2nd Piola--Kirchhoff stress tensor and its tangent from constitutive law

     B  = math_I3 - dt*Lpguess
     Fe = math_mul33x33(math_mul33x33(A,B), invFi_new)                                                  ! current elastic deformation tensor
     call constitutive_SandItsTangents(Tstar, dS_dFe, dS_dFi, &
                                      Fe, Fi_new, ipc, ip, el)                                          ! call constitutive law to calculate 2nd Piola-Kirchhoff stress and its derivative in unloaded configuration
     Tstar_v = math_Mandel33to6(Tstar)

     !* calculate plastic velocity gradient and its tangent from constitutive law

#ifdef DEBUG
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
         .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
       write(6,'(a,i3,/)')                  '<< CRYST >> stress iteration ', NiterationStressLp
       write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST >> Lpguess', transpose(Lpguess)
       write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST >> Fi', transpose(Fi_new)
       write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST >> Fe', transpose(Fe)
       write(6,'(a,/,6(e20.10,1x))')         '<< CRYST >> Tstar', Tstar_v
     endif
#endif
     call constitutive_LpAndItsTangents(Lp_constitutive, dLp_dS, dLp_dFi, &
                                       Tstar_v, Fi_new, ipc, ip, el)

#ifdef DEBUG
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
         .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
       write(6,'(a,/,3(12x,3(e20.10,1x)/))') '<< CRYST >> Lp_constitutive', transpose(Lp_constitutive)
     endif
#endif


     !* update current residuum and check for convergence of loop

     aTolLp = max(rTol_crystalliteStress * max(norm2(Lpguess),norm2(Lp_constitutive)), &             ! absolute tolerance from largest acceptable relative error
                  aTol_crystalliteStress)                                                            ! minimum lower cutoff
     residuumLp = Lpguess - Lp_constitutive

     if (any(IEEE_is_NaN(residuumLp))) then                                                          ! NaN in residuum...
#ifdef DEBUG
       if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) &
         write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3,a,i3,a)') '<< CRYST >> integrateStress encountered NaN at el (elFE) ip ipc ', &
           el,'(',mesh_element(1,el),')',ip,ipc, &
                                                        ' ; iteration ', NiterationStressLp,&
                                                        ' >> returning..!'
#endif
       return                                                                                         ! ...me = .false. to inform integrator about problem
     elseif (norm2(residuumLp) < aTolLp) then                                                         ! converged if below absolute tolerance
       exit LpLoop                                                                                    ! ...leave iteration loop
     elseif (     NiterationStressLp == 1_pInt &
             .or. norm2(residuumLp) < norm2(residuumLp_old)) then                                     ! not converged, but improved norm of residuum (always proceed in first iteration)...
       residuumLp_old = residuumLp                                                                    ! ...remember old values and...
       Lpguess_old    = Lpguess
       steplengthLp   = 1.0_pReal                                                                     ! ...proceed with normal step length (calculate new search direction)
     else                                                                                             ! not converged and residuum not improved...
       steplengthLp = subStepSizeLp * steplengthLp                                                    ! ...try with smaller step length in same direction
       Lpguess    = Lpguess_old + steplengthLp * deltaLp
#ifdef DEBUG
       if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
           .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                  .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
         write(6,'(a,1x,f7.4)') '<< CRYST >> linear search for Lpguess with step', steplengthLp
       endif
#endif
       cycle LpLoop
     endif


     !* calculate Jacobian for correction term

     if (mod(jacoCounterLp, iJacoLpresiduum) == 0_pInt) then
       dFe_dLp = 0.0_pReal
       forall(o=1_pInt:3_pInt,p=1_pInt:3_pInt) &
         dFe_dLp(o,1:3,p,1:3) = A(o,p)*transpose(invFi_new)                                ! dFe_dLp(i,j,k,l) = -dt * A(i,k) invFi(l,j)
       dFe_dLp = - dt * dFe_dLp
       dRLp_dLp    =   math_identity2nd(9_pInt) &
                     - math_Plain3333to99(math_mul3333xx3333(math_mul3333xx3333(dLp_dS,dS_dFe),dFe_dLp))
#ifdef DEBUG
       if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
           .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                  .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
         write(6,'(a,/,9(12x,9(e12.4,1x)/))') '<< CRYST >> dLp_dS', math_Plain3333to99(dLp_dS)
         write(6,'(a,1x,e20.10)') '<< CRYST >> dLp_dS norm', norm2(math_Plain3333to99(dLp_dS))
         write(6,'(a,/,9(12x,9(e12.4,1x)/))') '<< CRYST >> dRLp_dLp', dRLp_dLp - math_identity2nd(9_pInt)
         write(6,'(a,1x,e20.10)') '<< CRYST >> dRLp_dLp norm', norm2(dRLp_dLp - math_identity2nd(9_pInt))
       endif
#endif
       dRLp_dLp2   = dRLp_dLp                                                                         ! will be overwritten in first call to LAPACK routine
       work = math_plain33to9(residuumLp)
       call dgesv(9,1,dRLp_dLp2,9,ipiv,work,9,ierr)                                                   ! solve dRLp/dLp * delta Lp = -res for delta Lp
       if (ierr /= 0_pInt) then
#ifdef DEBUG
         if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
           write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3)') '<< CRYST >> integrateStress failed on dR/dLp inversion at el (elFE) ip ipc ', &
             el,'(',mesh_element(1,el),')',ip,ipc
           if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
               .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g)&
                      .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
             write(6,*)
             write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dR_dLp',transpose(dRLp_dLp)
             write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dFe_dLp',transpose(math_Plain3333to99(dFe_dLp))
             write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dS_dFe_constitutive',transpose(math_Plain3333to99(dS_dFe))
             write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dLp_dS_constitutive',transpose(math_Plain3333to99(dLp_dS))
             write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> A',transpose(A)
             write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> B',transpose(B)
             write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Lp_constitutive',transpose(Lp_constitutive)
             write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Lpguess',transpose(Lpguess)
           endif
         endif
#endif
         return
       endif
       deltaLp = - math_plain9to33(work)
     endif
     jacoCounterLp = jacoCounterLp + 1_pInt                                                         ! increase counter for jaco update

     Lpguess = Lpguess + steplengthLp * deltaLp

   enddo LpLoop

   !* calculate intermediate velocity gradient and its tangent from constitutive law

   call constitutive_LiAndItsTangents(Li_constitutive, dLi_dS, dLi_dFi, &
                                     Tstar_v, Fi_new, ipc, ip, el)

#ifdef DEBUG
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
         .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
       write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Li_constitutive', transpose(Li_constitutive)
       write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Liguess', transpose(Liguess)
     endif
#endif
   !* update current residuum and check for convergence of loop

   aTolLi = max(rTol_crystalliteStress * max(norm2(Liguess),norm2(Li_constitutive)), &              ! absolute tolerance from largest acceptable relative error
                aTol_crystalliteStress)                                                             ! minimum lower cutoff
   residuumLi = Liguess - Li_constitutive
   if (any(IEEE_is_NaN(residuumLi))) then                                                           ! NaN in residuum...
     return                                                                                         ! ...me = .false. to inform integrator about problem
   elseif (norm2(residuumLi) < aTolLi) then                                                         ! converged if below absolute tolerance
     exit LiLoop                                                                                    ! ...leave iteration loop
   elseif (     NiterationStressLi == 1_pInt &
           .or. norm2(residuumLi) < norm2(residuumLi_old)) then                                     ! not converged, but improved norm of residuum (always proceed in first iteration)...
     residuumLi_old = residuumLi                                                                    ! ...remember old values and...
     Liguess_old    = Liguess
     steplengthLi   = 1.0_pReal                                                                     ! ...proceed with normal step length (calculate new search direction)
   else                                                                                             ! not converged and residuum not improved...
     steplengthLi   = subStepSizeLi * steplengthLi                                                  ! ...try with smaller step length in same direction
     Liguess        = Liguess_old + steplengthLi * deltaLi
     cycle LiLoop
   endif

   !* calculate Jacobian for correction term

   if (mod(jacoCounterLi, iJacoLpresiduum) == 0_pInt) then
     temp_33     = math_mul33x33(math_mul33x33(A,B),invFi_current)
     dFe_dLi = 0.0_pReal
     dFi_dLi = 0.0_pReal
     forall(o=1_pInt:3_pInt,p=1_pInt:3_pInt)
       dFe_dLi(1:3,o,1:3,p) = -dt*math_I3(o,p)*temp_33                                          ! dFe_dLp(i,j,k,l) = -dt * A(i,k) invFi(l,j)
       dFi_dLi(1:3,o,1:3,p) = -dt*math_I3(o,p)*invFi_current
     end forall
     forall(o=1_pInt:3_pInt,p=1_pInt:3_pInt) &
       dFi_dLi(1:3,1:3,o,p) = math_mul33x33(math_mul33x33(Fi_new,dFi_dLi(1:3,1:3,o,p)),Fi_new)

     dRLi_dLi  = math_identity2nd(9_pInt) &
               - math_Plain3333to99(math_mul3333xx3333(dLi_dS, math_mul3333xx3333(dS_dFe, dFe_dLi) + &
                                                                   math_mul3333xx3333(dS_dFi, dFi_dLi)))  &
               - math_Plain3333to99(math_mul3333xx3333(dLi_dFi, dFi_dLi))
     work = math_plain33to9(residuumLi)
     call dgesv(9,1,dRLi_dLi,9,ipiv,work,9,ierr)                                                    ! solve dRLi/dLp * delta Li = -res for delta Li
     if (ierr /= 0_pInt) then
#ifdef DEBUG
       if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
         write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3)') '<< CRYST >> integrateStress failed on dR/dLi inversion at el (elFE) ip ipc ', &
               el,'(',mesh_element(1,el),')',ip,ipc
         if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
             .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g)&
                    .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
           write(6,*)
           write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dR_dLi',transpose(dRLi_dLi)
           write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dFe_dLi',transpose(math_Plain3333to99(dFe_dLi))
           write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dS_dFi_constitutive',transpose(math_Plain3333to99(dS_dFi))
           write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dLi_dS_constitutive',transpose(math_Plain3333to99(dLi_dS))
           write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Li_constitutive',transpose(Li_constitutive)
           write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Liguess',transpose(Liguess)
         endif
       endif
#endif
       return
     endif

     deltaLi = - math_plain9to33(work)
   endif
   jacoCounterLi = jacoCounterLi + 1_pInt                                                           ! increase counter for jaco update

   Liguess = Liguess + steplengthLi * deltaLi
 enddo LiLoop

 !* calculate new plastic and elastic deformation gradient

 invFp_new = math_mul33x33(invFp_current,B)
 invFp_new = invFp_new / math_det33(invFp_new)**(1.0_pReal/3.0_pReal)                               ! regularize by det
 Fp_new = math_inv33(invFp_new)
 failedInversionInvFp: if (all(dEq0(Fp_new))) then
#ifdef DEBUG
   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3,a,i3)') '<< CRYST >> integrateStress failed on invFp_new inversion at el (elFE) ip ipc ',&
       el,'(',mesh_element(1,el),')',ip,ipc, ' ; iteration ', NiterationStressLp
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
         .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) &
       write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> invFp_new',transpose(invFp_new)
   endif
#endif
   return
 endif failedInversionInvFp
 Fe_new = math_mul33x33(math_mul33x33(Fg_new,invFp_new),invFi_new)    ! calc resulting Fe

 !* calculate 1st Piola-Kirchhoff stress

 crystallite_P(1:3,1:3,ipc,ip,el) = math_mul33x33(math_mul33x33(Fg_new,invFp_new), &
                                              math_mul33x33(math_Mandel6to33(Tstar_v), &
                                                            transpose(invFp_new)))

 !* store local values in global variables

 crystallite_Lp(1:3,1:3,ipc,ip,el)    = Lpguess
 crystallite_Li(1:3,1:3,ipc,ip,el)    = Liguess
 crystallite_Tstar_v(1:6,ipc,ip,el)   = Tstar_v
 crystallite_Fp(1:3,1:3,ipc,ip,el)    = Fp_new
 crystallite_Fi(1:3,1:3,ipc,ip,el)    = Fi_new
 crystallite_Fe(1:3,1:3,ipc,ip,el)    = Fe_new
 crystallite_invFp(1:3,1:3,ipc,ip,el) = invFp_new
 crystallite_invFi(1:3,1:3,ipc,ip,el) = invFi_new

 !* set return flag to true

 crystallite_integrateStress = .true.
#ifdef DEBUG
 if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt &
     .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
             .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> P / MPa',transpose(crystallite_P(1:3,1:3,ipc,ip,el))*1.0e-6_pReal
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Cauchy / MPa', &
              math_mul33x33(crystallite_P(1:3,1:3,ipc,ip,el), transpose(Fg_new)) * 1.0e-6_pReal / math_det33(Fg_new)
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fe Lp Fe^-1', &
              transpose(math_mul33x33(Fe_new, math_mul33x33(crystallite_Lp(1:3,1:3,ipc,ip,el), math_inv33(Fe_new))))    ! transpose to get correct print out order
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fp',transpose(crystallite_Fp(1:3,1:3,ipc,ip,el))
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fi',transpose(crystallite_Fi(1:3,1:3,ipc,ip,el))
 endif
#endif

end function crystallite_integrateStress


!--------------------------------------------------------------------------------------------------
!> @brief calculates orientations and disorientations (in case of single grain ips)
!--------------------------------------------------------------------------------------------------
subroutine crystallite_orientations
 use math, only: &
   math_rotationalPart33, &
   math_RtoQ, &
   math_qConj
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP
 use material, only: &
   material_phase, &
   homogenization_Ngrains, &
   plasticState
 use mesh, only: &
   mesh_element, &
   mesh_ipNeighborhood, &
   FE_NipNeighbors, &
   FE_geomtype, &
   FE_celltype
 use lattice, only: &
   lattice_qDisorientation, &
   lattice_structure
 use plastic_nonlocal, only: &
   plastic_nonlocal_updateCompatibility

 implicit none
 integer(pInt) &
   c, &                                                                                             !< counter in integration point component loop
   i, &                                                                                             !< counter in integration point loop
   e, &                                                                                             !< counter in element loop
   n, &                                                                                             !< counter in neighbor loop
   neighboring_e, &                                                                                 !< neighbor element
   neighboring_i, &                                                                                 !< neighbor integration point
   myPhase, &                    ! phase
   neighboringPhase
 real(pReal), dimension(4) :: &
   orientation

 ! --- CALCULATE ORIENTATION AND LATTICE ROTATION ---

!$OMP PARALLEL DO PRIVATE(orientation)
 do e = FEsolving_execElem(1),FEsolving_execElem(2)
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do c = 1_pInt,homogenization_Ngrains(mesh_element(3,e))
! somehow this subroutine is not threadsafe, so need critical statement here; not clear, what exactly the problem is
!$OMP CRITICAL (polarDecomp)
       orientation = math_RtoQ(transpose(math_rotationalPart33(crystallite_Fe(1:3,1:3,c,i,e))))
!$OMP END CRITICAL (polarDecomp)
       crystallite_rotation(1:4,c,i,e) = lattice_qDisorientation(crystallite_orientation0(1:4,c,i,e), &! active rotation from initial
                                                                  orientation)                         ! to current orientation (with no symmetry)
       crystallite_orientation(1:4,c,i,e) = orientation
 enddo; enddo; enddo
!$OMP END PARALLEL DO
 
  
 ! --- UPDATE SOME ADDITIONAL VARIABLES THAT ARE NEEDED FOR NONLOCAL MATERIAL ---
 ! --- we use crystallite_orientation from above, so need a separate loop
  
 nonlocalPresent: if (any(plasticState%nonLocal)) then
!$OMP PARALLEL DO PRIVATE(myPhase,neighboring_e,neighboring_i,neighboringPhase)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       myPhase = material_phase(1,i,e)                                                              ! get my phase (non-local models make no sense with more than one grain per material point)
       if (plasticState(myPhase)%nonLocal) then                                                     ! if nonlocal model
         ! --- calculate disorientation between me and my neighbor ---

         do n = 1_pInt,FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,e))))                 ! loop through my neighbors
           neighboring_e = mesh_ipNeighborhood(1,n,i,e)
           neighboring_i = mesh_ipNeighborhood(2,n,i,e)
           if (neighboring_e > 0 .and. neighboring_i > 0) then                                      ! if neighbor exists
             neighboringPhase = material_phase(1,neighboring_i,neighboring_e)                       ! get my neighbor's phase
             if (plasticState(neighboringPhase)%nonLocal) then                                      ! neighbor got also nonlocal plasticity
               if (lattice_structure(myPhase) == lattice_structure(neighboringPhase)) then          ! if my neighbor has same crystal structure like me
                 crystallite_disorientation(:,n,1,i,e) = &
                   lattice_qDisorientation( crystallite_orientation(1:4,1,i,e), &
                                            crystallite_orientation(1:4,1,neighboring_i,neighboring_e), &
                                            lattice_structure(myPhase))                             ! calculate disorientation for given symmetry
               else                                                                                 ! for neighbor with different phase
                 crystallite_disorientation(:,n,1,i,e) = [0.0_pReal, 1.0_pReal, 0.0_pReal, 0.0_pReal]! 180 degree rotation about 100 axis
               endif
             else                                                                                   ! for neighbor with local plasticity
               crystallite_disorientation(:,n,1,i,e) = [-1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal]! homomorphic identity
             endif
           else                                                                                     ! no existing neighbor
             crystallite_disorientation(:,n,1,i,e) = [-1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal]  ! homomorphic identity
           endif
         enddo


         ! --- calculate compatibility and transmissivity between me and my neighbor ---

         call plastic_nonlocal_updateCompatibility(crystallite_orientation,i,e)

       endif
   enddo; enddo
!$OMP END PARALLEL DO
 endif nonlocalPresent

end subroutine crystallite_orientations

!--------------------------------------------------------------------------------------------------
!> @brief return results of particular grain
!--------------------------------------------------------------------------------------------------
function crystallite_postResults(ipc, ip, el)
 use math, only: &
   math_qToEuler, &
   math_qToEulerAxisAngle, &
   math_mul33x33, &
   math_det33, &
   math_I3, &
   inDeg, &
   math_Mandel6to33, &
   math_qMul, &
   math_qConj
 use mesh, only: &
   mesh_element, &
   mesh_ipVolume, &
   mesh_maxNipNeighbors, &
   mesh_ipNeighborhood, &
   FE_NipNeighbors, &
   FE_geomtype, &
   FE_celltype
 use material, only: &
   plasticState, &
   sourceState, &
   microstructure_crystallite, &
   crystallite_Noutput, &
   material_phase, &
   material_texture, &
   homogenization_Ngrains
 use constitutive, only: &
   constitutive_homogenizedC, &
   constitutive_postResults

 implicit none
 integer(pInt), intent(in):: &
   el, &                         !< element index
   ip, &                         !< integration point index
   ipc                           !< grain index

 real(pReal), dimension(1+crystallite_sizePostResults(microstructure_crystallite(mesh_element(4,el))) + &
                        1+plasticState(material_phase(ipc,ip,el))%sizePostResults + &
                          sum(sourceState(material_phase(ipc,ip,el))%p(:)%sizePostResults)) :: &
   crystallite_postResults
 real(pReal), dimension(3,3) :: &
   Ee
 real(pReal), dimension(4) :: &
   rotation
 real(pReal) :: &
   detF
 integer(pInt) :: &
   o, &
   c, &
   crystID, &
   mySize, &
   n


 crystID = microstructure_crystallite(mesh_element(4,el))

 crystallite_postResults = 0.0_pReal
 c = 0_pInt
 crystallite_postResults(c+1) = real(crystallite_sizePostResults(crystID),pReal)                  ! size of results from cryst
 c = c + 1_pInt

 do o = 1_pInt,crystallite_Noutput(crystID)
   mySize = 0_pInt
   select case(crystallite_outputID(o,crystID))
     case (phase_ID)
       mySize = 1_pInt
       crystallite_postResults(c+1) = real(material_phase(ipc,ip,el),pReal)                         ! phaseID of grain
     case (texture_ID)
       mySize = 1_pInt
       crystallite_postResults(c+1) = real(material_texture(ipc,ip,el),pReal)                       ! textureID of grain
     case (volume_ID)
       mySize = 1_pInt
       detF = math_det33(crystallite_partionedF(1:3,1:3,ipc,ip,el))                                 ! V_current = det(F) * V_reference
       crystallite_postResults(c+1) = detF * mesh_ipVolume(ip,el) &
                                           / real(homogenization_Ngrains(mesh_element(3,el)),pReal) ! grain volume (not fraction but absolute)
     case (orientation_ID)
       mySize = 4_pInt
       crystallite_postResults(c+1:c+mySize) = crystallite_orientation(1:4,ipc,ip,el)               ! grain orientation as quaternion
     case (eulerangles_ID)
       mySize = 3_pInt
       crystallite_postResults(c+1:c+mySize) = inDeg &
                                             * math_qToEuler(crystallite_orientation(1:4,ipc,ip,el)) ! grain orientation as Euler angles in degree
     case (grainrotation_ID)
       mySize = 4_pInt
       crystallite_postResults(c+1:c+mySize) = &
           math_qToEulerAxisAngle(crystallite_rotation(1:4,ipc,ip,el))                              ! grain rotation away from initial orientation as axis-angle in sample reference coordinates
       crystallite_postResults(c+4) = inDeg * crystallite_postResults(c+4)                          ! angle in degree
     case (grainrotationx_ID)
       mySize = 1_pInt
       rotation = math_qToEulerAxisAngle(crystallite_rotation(1:4,ipc,ip,el))                       ! grain rotation away from initial orientation as axis-angle in sample reference coordinates
       crystallite_postResults(c+1) = inDeg * rotation(1) * rotation(4)                             ! angle in degree
     case (grainrotationy_ID)
       mySize = 1_pInt
       rotation = math_qToEulerAxisAngle(crystallite_rotation(1:4,ipc,ip,el))                       ! grain rotation away from initial orientation as axis-angle in sample reference coordinates
       crystallite_postResults(c+1) = inDeg * rotation(2) * rotation(4)                             ! angle in degree
     case (grainrotationz_ID)
       mySize = 1_pInt
       rotation = math_qToEulerAxisAngle(crystallite_rotation(1:4,ipc,ip,el))                       ! grain rotation away from initial orientation as axis-angle in sample reference coordinates
       crystallite_postResults(c+1) = inDeg * rotation(3) * rotation(4)                             ! angle in degree

! remark: tensor output is of the form 11,12,13, 21,22,23, 31,32,33
! thus row index i is slow, while column index j is fast. reminder: "row is slow"

     case (defgrad_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_partionedF(1:3,1:3,ipc,ip,el)),[mySize])
     case (e_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = 0.5_pReal * reshape((math_mul33x33( &
                                               transpose(crystallite_partionedF(1:3,1:3,ipc,ip,el)), &
                                               crystallite_partionedF(1:3,1:3,ipc,ip,el)) - math_I3),[mySize])
     case (fe_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_Fe(1:3,1:3,ipc,ip,el)),[mySize])
     case (ee_ID)
       Ee = 0.5_pReal *(math_mul33x33(transpose(crystallite_Fe(1:3,1:3,ipc,ip,el)), &
                                               crystallite_Fe(1:3,1:3,ipc,ip,el)) - math_I3)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = reshape(Ee,[mySize])
     case (fp_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_Fp(1:3,1:3,ipc,ip,el)),[mySize])
     case (fi_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_Fi(1:3,1:3,ipc,ip,el)),[mySize])
     case (lp_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_Lp(1:3,1:3,ipc,ip,el)),[mySize])
     case (li_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_Li(1:3,1:3,ipc,ip,el)),[mySize])
     case (p_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(transpose(crystallite_P(1:3,1:3,ipc,ip,el)),[mySize])
     case (s_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(math_Mandel6to33(crystallite_Tstar_v(1:6,ipc,ip,el)),[mySize])
     case (elasmatrix_ID)
       mySize = 36_pInt
       crystallite_postResults(c+1:c+mySize) = reshape(constitutive_homogenizedC(ipc,ip,el),[mySize])
     case(neighboringelement_ID)
       mySize = mesh_maxNipNeighbors
       crystallite_postResults(c+1:c+mySize) = 0.0_pReal
       forall (n = 1_pInt:FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,el))))) &
         crystallite_postResults(c+n) = real(mesh_ipNeighborhood(1,n,ip,el),pReal)
     case(neighboringip_ID)
       mySize = mesh_maxNipNeighbors
       crystallite_postResults(c+1:c+mySize) = 0.0_pReal
       forall (n = 1_pInt:FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,el))))) &
         crystallite_postResults(c+n) = real(mesh_ipNeighborhood(2,n,ip,el),pReal)
   end select
   c = c + mySize
 enddo

 crystallite_postResults(c+1) = real(plasticState(material_phase(ipc,ip,el))%sizePostResults,pReal)             ! size of constitutive results
 c = c + 1_pInt
 if (size(crystallite_postResults)-c > 0_pInt) &
   crystallite_postResults(c+1:size(crystallite_postResults)) = &
      constitutive_postResults(crystallite_Tstar_v(1:6,ipc,ip,el), crystallite_Fe, &
                               ipc, ip, el)

end function crystallite_postResults

end module crystallite
