!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters,     Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords,  Max-Planck-Institut für Eisenforschung GmbH
!> @author Chen Zhang,       Michigan State University
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
   crystallite_subFe0,&                                                                             !< "elastic" def grad at start of crystallite inc
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
 real(pReal),                dimension(:,:,:,:,:,:,:), allocatable, private :: &
   crystallite_fallbackdPdF                                                                         !< dPdF fallback for non-converged grains (elastic prediction)
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

 public :: &
   crystallite_init, &
   crystallite_stressAndItsTangent, &
   crystallite_orientations, &
   crystallite_push33ToRef, &
   crystallite_postResults
 private :: &
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
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_info, &
   debug_reset, &
   debug_level, &
   debug_crystallite, &
   debug_levelBasic
 use numerics, only: &
   worldrank, &
   usePingPong
 use math, only: &
   math_I3, &
   math_EulerToR, &
   math_inv33, &
   math_transpose33, &
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
   IO_read, &
   IO_timeStamp, &
   IO_open_jobFile_stat, &
   IO_open_file, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_write_jobFile, &
   IO_error, &
   IO_EOF
 use material
 use constitutive, only: &
   constitutive_initialFi, &
   constitutive_microstructure                                                                     ! derived (shortcut) quantities of given state

 implicit none
 integer(pInt), parameter :: &
   FILEUNIT = 200_pInt

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: &
   c, &                                                                                             !< counter in integration point component loop
   i, &                                                                                             !< counter in integration point loop
   e, &                                                                                             !< counter in element loop
   o, &                                                                                             !< counter in output loop
   r, &                                                                                             !< counter in crystallite loop
   cMax, &                                                                                          !< maximum number of  integration point components
   iMax, &                                                                                          !< maximum number of integration points
   eMax, &                                                                                          !< maximum number of elements
   nMax, &                                                                                          !< maximum number of ip neighbors
   myNcomponents, &                                                                                 !< number of components at current IP
   section = 0_pInt, &
   j, &
   p, &
   mySize

 character(len=65536) :: &
   tag = '', &
   line= ''

 mainProcess: if (worldrank == 0) then
   write(6,'(/,a)')   ' <<<+-  crystallite init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

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
 allocate(crystallite_subFe0(3,3,cMax,iMax,eMax),            source=0.0_pReal)
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
 allocate(crystallite_fallbackdPdF(3,3,3,3,cMax,iMax,eMax),  source=0.0_pReal)
 allocate(crystallite_dt(cMax,iMax,eMax),                    source=0.0_pReal)
 allocate(crystallite_subdt(cMax,iMax,eMax),                 source=0.0_pReal)
 allocate(crystallite_subFrac(cMax,iMax,eMax),               source=0.0_pReal)
 allocate(crystallite_subStep(cMax,iMax,eMax),               source=0.0_pReal)
 allocate(crystallite_orientation(4,cMax,iMax,eMax),         source=0.0_pReal)
 allocate(crystallite_orientation0(4,cMax,iMax,eMax),        source=0.0_pReal)
 allocate(crystallite_rotation(4,cMax,iMax,eMax),            source=0.0_pReal)
 allocate(crystallite_disorientation(4,nMax,cMax,iMax,eMax), source=0.0_pReal)
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
                             material_Ncrystallite)) ;       crystallite_output = ''
 allocate(crystallite_outputID(maxval(crystallite_Noutput), &
                             material_Ncrystallite),         source=undefined_ID)
 allocate(crystallite_sizePostResults(material_Ncrystallite),source=0_pInt)
 allocate(crystallite_sizePostResult(maxval(crystallite_Noutput), &
                                     material_Ncrystallite), source=0_pInt)

 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  !  no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ...open material.config file

 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partCrystallite)   ! wind forward to <crystallite>
   line = IO_read(FILEUNIT)
 enddo

 do while (trim(line) /= IO_EOF)                                                                    ! read through sections of crystallite part
   line = IO_read(FILEUNIT)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(FILEUNIT, .true.)                                                               ! reset IO_read
     exit
   endif
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     o = 0_pInt                                                                                     ! reset output counter
     cycle                                                                                          ! skip to next line
   endif
   if (section > 0_pInt) then
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                              ! extract key
     select case(tag)
       case ('(output)')
         o = o + 1_pInt
         crystallite_output(o,section) = IO_lc(IO_stringValue(line,chunkPos,2_pInt))
         outputName: select case(crystallite_output(o,section))
           case ('phase') outputName
             crystallite_outputID(o,section) = phase_ID
           case ('texture') outputName
             crystallite_outputID(o,section) = texture_ID
           case ('volume') outputName
             crystallite_outputID(o,section) = volume_ID
           case ('grainrotationx') outputName
             crystallite_outputID(o,section) = grainrotationx_ID
           case ('grainrotationy') outputName
             crystallite_outputID(o,section) = grainrotationy_ID
           case ('grainrotationz') outputName
             crystallite_outputID(o,section) = grainrotationx_ID
           case ('orientation') outputName
             crystallite_outputID(o,section) = orientation_ID
           case ('grainrotation') outputName
             crystallite_outputID(o,section) = grainrotation_ID
           case ('eulerangles') outputName
             crystallite_outputID(o,section) = eulerangles_ID
           case ('defgrad','f') outputName
             crystallite_outputID(o,section) = defgrad_ID
           case ('fe') outputName
             crystallite_outputID(o,section) = fe_ID
           case ('fp') outputName
             crystallite_outputID(o,section) = fp_ID
           case ('fi') outputName
             crystallite_outputID(o,section) = fi_ID
           case ('lp') outputName
             crystallite_outputID(o,section) = lp_ID
           case ('li') outputName
             crystallite_outputID(o,section) = li_ID
           case ('e') outputName
             crystallite_outputID(o,section) = e_ID
           case ('ee') outputName
             crystallite_outputID(o,section) = ee_ID
           case ('p','firstpiola','1stpiola') outputName
             crystallite_outputID(o,section) = p_ID
           case ('s','tstar','secondpiola','2ndpiola') outputName
             crystallite_outputID(o,section) = s_ID
           case ('elasmatrix') outputName
             crystallite_outputID(o,section) = elasmatrix_ID
           case ('neighboringip') outputName
             crystallite_outputID(o,section) = neighboringip_ID
           case ('neighboringelement') outputName
             crystallite_outputID(o,section) = neighboringelement_ID
           case default outputName
             call IO_error(105_pInt,ext_msg=IO_stringValue(line,chunkPos,2_pInt)//' (Crystallite)')
         end select outputName
     end select
   endif
 enddo

 close(FILEUNIT)

 do r = 1_pInt,material_Ncrystallite
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

   do r = 1_pInt,material_Ncrystallite
     if (any(microstructure_crystallite(mesh_element(4,:)) == r)) then
       write(FILEUNIT,'(/,a,/)') '['//trim(crystallite_name(r))//']'
       do o = 1_pInt,crystallite_Noutput(r)
         write(FILEUNIT,'(a,i4)') trim(crystallite_output(o,r))//char(9),crystallite_sizePostResult(o,r)
       enddo
     endif
   enddo

   close(FILEUNIT)
 endif

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
 crystallite_fallbackdPdF = crystallite_dPdF                                                        ! use initial elastic stiffness as fallback

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
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subFe0:                ', shape(crystallite_subFe0)
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
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_fallbackdPdF:          ', shape(crystallite_fallbackdPdF)
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
   tol_math_check
 use numerics, only: &
   subStepMinCryst, &
   subStepSizeCryst, &
   stepIncreaseCryst, &
   pert_Fg, &
   pert_method, &
   nCryst, &
   numerics_integrator, &
   numerics_integrationMode, &
   numerics_timeSyncing, &
   analyticJaco
 use debug, only: &
   debug_level, &
   debug_crystallite, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_g, &
   debug_CrystalliteLoopDistribution
 use IO, only: &
   IO_warning, &
   IO_error
 use math, only: &
   math_inv33, &
   math_identity2nd, &
   math_transpose33, &
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
   mesh_NcpElems, &
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
   phaseAt, phasememberAt, &
   homogenization_maxNgrains
 use constitutive, only:  &
   constitutive_TandItsTangent, &
   constitutive_LpAndItsTangent, &
   constitutive_LiAndItsTangent

 implicit none
 logical, intent(in) :: &
   updateJaco                                                                                       !< whether to update the Jacobian (stiffness) or not
 real(pReal) :: &
   myPert, &                                                                                        ! perturbation with correct sign
   formerSubStep, &
   subFracIntermediate
 real(pReal), dimension(3,3) :: &
   invFp, &                                                                                         ! inverse of the plastic deformation gradient
   Fe_guess, &                                                                                      ! guess for elastic deformation gradient
   Tstar                                                                                            ! 2nd Piola-Kirchhoff stress tensor
 real(pReal), allocatable, dimension(:,:,:,:,:,:,:) :: &
   dPdF_perturbation1, &
   dPdF_perturbation2
 real(pReal), allocatable, dimension(:,:,:,:,:) :: &
   F_backup, &
   Fp_backup, &
   InvFp_backup, &
   Fi_backup, &
   InvFi_backup, &
   Fe_backup, &
   Lp_backup, &
   Li_backup, &
   P_backup
 real(pReal), allocatable, dimension(:,:,:,:) :: &
   Tstar_v_backup
 logical,     allocatable, dimension(:,:,:) :: &
   convergenceFlag_backup
 integer(pInt) :: &
   NiterationCrystallite, &                                                                         ! number of iterations in crystallite loop
   c, &                                                                                             !< counter in integration point component loop
   i, &                                                                                             !< counter in integration point loop
   e, &                                                                                             !< counter in element loop
   k, &
   l, &
   n, startIP, endIP, &
   neighboring_e, &
   neighboring_i, &
   o, &
   p, &
   perturbation , &                                                                                 ! loop counter for forward,backward perturbation mode
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
                                         math_transpose33(crystallite_partionedF(1:3,1:3,debug_g,debug_i,debug_e))
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> F0 ', &
                                         math_transpose33(crystallite_partionedF0(1:3,1:3,debug_g,debug_i,debug_e))
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Fp0', &
                                         math_transpose33(crystallite_partionedFp0(1:3,1:3,debug_g,debug_i,debug_e))
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Fi0', &
                                         math_transpose33(crystallite_partionedFi0(1:3,1:3,debug_g,debug_i,debug_e))
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Lp0', &
                                         math_transpose33(crystallite_partionedLp0(1:3,1:3,debug_g,debug_i,debug_e))
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Li0', &
                                         math_transpose33(crystallite_partionedLi0(1:3,1:3,debug_g,debug_i,debug_e))
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
         crystallite_subFe0(1:3,1:3,c,i,e) = math_mul33x33(math_mul33x33(crystallite_subF0(1:3,1:3,c,i,e), &
                                                                         math_inv33(crystallite_subFp0(1:3,1:3,c,i,e))), &
                                                           math_inv33(crystallite_subFi0(1:3,1:3,c,i,e)))! only needed later on for stiffness calculation
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
 numerics_integrationMode = 1_pInt
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
         if (abs(subFracIntermediate) > tiny(0.0_pReal)) then
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
#ifndef _OPENMP
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
               if (.not. crystallite_localPlasticity(1,i,e) .and. abs(crystallite_subFrac(1,i,e)) > tiny(0.0_pReal)) then
                 do n = 1_pInt,FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,e))))
                   neighboring_e = mesh_ipNeighborhood(1,n,i,e)
                   neighboring_i = mesh_ipNeighborhood(2,n,i,e)
                   if (neighboring_e > 0_pInt .and. neighboring_i > 0_pInt) then
                     if (.not. crystallite_localPlasticity(1,neighboring_i,neighboring_e) &
                         .and. .not. crystallite_converged(1,neighboring_i,neighboring_e)) then
                       crystallite_syncSubFrac(i,e) = .true.
#ifndef _OPENMP
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
               crystallite_subFe0(1:3,1:3,c,i,e) = math_mul33x33(math_mul33x33(crystallite_subF (1:3,1:3,c,i,e), &
                                                                               crystallite_invFp(1:3,1:3,c,i,e)), &
                                                                 crystallite_invFi(1:3,1:3,c,i,e))  ! only needed later on for stiffness calculation
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
#ifndef _OPENMP
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
               if (iand(debug_level(debug_crystallite),debug_levelBasic) /= 0_pInt &
                   .and. formerSubStep > 0.0_pReal) then
                 !$OMP CRITICAL (distributionCrystallite)
                   debug_CrystalliteLoopDistribution(min(nCryst+1_pInt,NiterationCrystallite)) = &
                     debug_CrystalliteLoopDistribution(min(nCryst+1_pInt,NiterationCrystallite)) + 1_pInt
                 !$OMP END CRITICAL (distributionCrystallite)
               endif
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
#ifndef _OPENMP
             if (iand(debug_level(debug_crystallite),debug_levelBasic) /= 0_pInt) then
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
             .and. crystallite_subStep <= subStepMinCryst)) then                                      ! no way of rescuing a nonlocal ip that violated the lower time step limit, ...
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
         crystallite_todo = .false.                                                                       ! ... so let all nonlocal ips die peacefully
         crystallite_subStep = 0.0_pReal
       endwhere
     endif
   endif timeSyncing2

   if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
     write(6,'(/,a,e12.5)') '<< CRYST >> min(subStep) ',minval(crystallite_subStep)
     write(6,'(a,e12.5)')   '<< CRYST >> max(subStep) ',maxval(crystallite_subStep)
     write(6,'(a,e12.5)')   '<< CRYST >> min(subFrac) ',minval(crystallite_subFrac)
     write(6,'(a,e12.5,/)') '<< CRYST >> max(subFrac) ',maxval(crystallite_subFrac)
     flush(6)
   endif

   ! --- integrate --- requires fully defined state array (basic + dependent state)

   if (any(crystallite_todo)) then
     select case(numerics_integrator(numerics_integrationMode))
       case(1_pInt)
         call crystallite_integrateStateFPI()
       case(2_pInt)
         call crystallite_integrateStateEuler()
       case(3_pInt)
         call crystallite_integrateStateAdaptiveEuler()
       case(4_pInt)
         call crystallite_integrateStateRK4()
       case(5_pInt)
         call crystallite_integrateStateRKCK45()
     end select
   endif

   where(.not. crystallite_converged .and. crystallite_subStep > subStepMinCryst) &                  ! do not try non-converged & fully cutbacked any further
     crystallite_todo = .true.

   NiterationCrystallite = NiterationCrystallite + 1_pInt

 enddo cutbackLooping


! --+>> CHECK FOR NON-CONVERGED CRYSTALLITES <<+--

 elementLooping5: do e = FEsolving_execElem(1),FEsolving_execElem(2)
   myNcomponents = homogenization_Ngrains(mesh_element(3,e))
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                  ! iterate over IPs of this element to be processed
     do c = 1,myNcomponents
       if (.not. crystallite_converged(c,i,e)) then                                                    ! respond fully elastically (might be not required due to becoming terminally ill anyway)
         if(iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) &
           write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3,/)') '<< CRYST >> no convergence: respond fully elastic at el (elFE) ip ipc ', &
             e,'(',mesh_element(1,e),')',i,c
         invFp = math_inv33(crystallite_partionedFp0(1:3,1:3,c,i,e))
         Fe_guess = math_mul33x33(math_mul33x33(crystallite_partionedF(1:3,1:3,c,i,e), invFp), &
                                  math_inv33(crystallite_partionedFi0(1:3,1:3,c,i,e)))
         call constitutive_TandItsTangent(Tstar,dSdFe,dSdFi,Fe_guess,crystallite_partionedFi0(1:3,1:3,c,i,e),c,i,e)
         crystallite_P(1:3,1:3,c,i,e) = math_mul33x33(math_mul33x33(crystallite_partionedF(1:3,1:3,c,i,e), invFp), &
                                                      math_mul33x33(Tstar,transpose(invFp)))
       endif
       if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
           .and. ((e == debug_e .and. i == debug_i .and. c == debug_g) &
                  .or. .not. iand(debug_level(debug_crystallite),debug_levelSelective) /= 0_pInt)) then
         write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> central solution of cryst_StressAndTangent at el ip ipc ',e,i,c
         write(6,'(/,a,/,3(12x,3(f12.4,1x)/))') '<< CRYST >> P / MPa', &
                                          math_transpose33(crystallite_P(1:3,1:3,c,i,e))*1.0e-6_pReal
         write(6,'(a,/,3(12x,3(f14.9,1x)/))')   '<< CRYST >> Fp', &
                                          math_transpose33(crystallite_Fp(1:3,1:3,c,i,e))
         write(6,'(a,/,3(12x,3(f14.9,1x)/))')   '<< CRYST >> Fi', &
                                          math_transpose33(crystallite_Fi(1:3,1:3,c,i,e))
         write(6,'(a,/,3(12x,3(f14.9,1x)/),/)') '<< CRYST >> Lp', &
                                          math_transpose33(crystallite_Lp(1:3,1:3,c,i,e))
         write(6,'(a,/,3(12x,3(f14.9,1x)/),/)') '<< CRYST >> Li', &
                                          math_transpose33(crystallite_Li(1:3,1:3,c,i,e))
         flush(6)
       endif
     enddo
   enddo
 enddo elementLooping5


! --+>> STIFFNESS CALCULATION <<+--

 computeJacobian: if(updateJaco) then
   jacobianMethod: if (analyticJaco) then

     ! --- ANALYTIC JACOBIAN ---

     !$OMP PARALLEL DO PRIVATE(dSdF,dSdFe,dSdFi,dLpdS,dLpdFi,dFpinvdF,dLidS,dLidFi,dFidS,&
     !$OMP                     rhs_3333,lhs_3333,temp_99,temp_33,temp_3333,myNcomponents,error)
       elementLooping6: do e = FEsolving_execElem(1),FEsolving_execElem(2)
         myNcomponents = homogenization_Ngrains(mesh_element(3,e))
         do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                            ! iterate over IPs of this element to be processed
           do c = 1_pInt,myNcomponents
             call constitutive_TandItsTangent(temp_33,dSdFe,dSdFi,crystallite_Fe(1:3,1:3,c,i,e), &
                                              crystallite_Fi(1:3,1:3,c,i,e),c,i,e)                     ! call constitutive law to calculate elastic stress tangent

             call constitutive_LiAndItsTangent(temp_33,dLidS,dLidFi,crystallite_Tstar_v(1:6,c,i,e), &
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

             call constitutive_LpAndItsTangent(temp_33,dLpdS,dLpdFi,crystallite_Tstar_v(1:6,c,i,e), &
                                               crystallite_Fi(1:3,1:3,c,i,e),c,i,e)                    ! call constitutive law to calculate Lp tangent in lattice configuration
             dLpdS = math_mul3333xx3333(dLpdFi,dFidS) + dLpdS

             temp_33   = math_transpose33(math_mul33x33(crystallite_invFp(1:3,1:3,c,i,e), &
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
                                                math_transpose33(crystallite_invFp(1:3,1:3,c,i,e))))
             forall(p=1_pInt:3_pInt) &
               crystallite_dPdF(p,1:3,p,1:3,c,i,e) = math_transpose33(temp_33)

             temp_33 = math_mul33x33(math_Mandel6to33(crystallite_Tstar_v(1:6,c,i,e)), &
                                  math_transpose33(crystallite_invFp(1:3,1:3,c,i,e)))
             forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
               crystallite_dPdF(1:3,1:3,p,o,c,i,e) = crystallite_dPdF(1:3,1:3,p,o,c,i,e) + &
                 math_mul33x33(math_mul33x33(crystallite_subF(1:3,1:3,c,i,e),dFpinvdF(1:3,1:3,p,o)),temp_33)

             temp_33 = math_mul33x33(crystallite_subF(1:3,1:3,c,i,e), &
                                  crystallite_invFp(1:3,1:3,c,i,e))
             forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
               crystallite_dPdF(1:3,1:3,p,o,c,i,e) = crystallite_dPdF(1:3,1:3,p,o,c,i,e) + &
                 math_mul33x33(math_mul33x33(temp_33,dSdF(1:3,1:3,p,o)), &
                               math_transpose33(crystallite_invFp(1:3,1:3,c,i,e)))

             temp_33 = math_mul33x33(math_mul33x33(crystallite_subF(1:3,1:3,c,i,e), &
                                                crystallite_invFp(1:3,1:3,c,i,e)), &
                                  math_Mandel6to33(crystallite_Tstar_v(1:6,c,i,e)))
             forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
               crystallite_dPdF(1:3,1:3,p,o,c,i,e) = crystallite_dPdF(1:3,1:3,p,o,c,i,e) + &
                 math_mul33x33(temp_33,math_transpose33(dFpinvdF(1:3,1:3,p,o)))

         enddo; enddo
       enddo elementLooping6
     !$OMP END PARALLEL DO

   else jacobianMethod

     ! --- STANDARD (PERTURBATION METHOD) FOR JACOBIAN ---

     numerics_integrationMode = 2_pInt

     ! --- BACKUP ---
     allocate(dPdF_perturbation1(3,3,3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source = 0.0_pReal)
     allocate(dPdF_perturbation2(3,3,3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source = 0.0_pReal)
     allocate(F_backup          (3,3,    homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source = 0.0_pReal)
     allocate(Fp_backup         (3,3,    homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source = 0.0_pReal)
     allocate(InvFp_backup      (3,3,    homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source = 0.0_pReal)
     allocate(Fi_backup         (3,3,    homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source = 0.0_pReal)
     allocate(InvFi_backup      (3,3,    homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source = 0.0_pReal)
     allocate(Fe_backup         (3,3,    homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source = 0.0_pReal)
     allocate(Lp_backup         (3,3,    homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source = 0.0_pReal)
     allocate(Li_backup         (3,3,    homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source = 0.0_pReal)
     allocate(P_backup          (3,3,    homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source = 0.0_pReal)
     allocate(Tstar_v_backup    (6,      homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source = 0.0_pReal)
     allocate(convergenceFlag_backup    (homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source = .false.)

     !$OMP PARALLEL DO PRIVATE(myNcomponents)
       elementLooping7: do e = FEsolving_execElem(1),FEsolving_execElem(2)
         myNcomponents = homogenization_Ngrains(mesh_element(3,e))
         do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e); do c = 1,myNcomponents

           plasticState    (phaseAt(c,i,e))%state_backup(:,phasememberAt(c,i,e)) = &
           plasticState    (phaseAt(c,i,e))%state(       :,phasememberAt(c,i,e))
           do mySource = 1_pInt, phase_Nsources(phaseAt(c,i,e))
             sourceState(phaseAt(c,i,e))%p(mySource)%state_backup(:,phasememberAt(c,i,e)) = &
             sourceState(phaseAt(c,i,e))%p(mySource)%state(       :,phasememberAt(c,i,e))
           enddo

           plasticState    (phaseAt(c,i,e))%dotState_backup(:,phasememberAt(c,i,e)) = &
           plasticState    (phaseAt(c,i,e))%dotState(       :,phasememberAt(c,i,e))
           do mySource = 1_pInt, phase_Nsources(phaseAt(c,i,e))
             sourceState(phaseAt(c,i,e))%p(mySource)%dotState_backup(:,phasememberAt(c,i,e)) = &
             sourceState(phaseAt(c,i,e))%p(mySource)%dotState(       :,phasememberAt(c,i,e))
           enddo

           F_backup(1:3,1:3,c,i,e)       = crystallite_subF(1:3,1:3,c,i,e)                            ! ... and kinematics
           Fp_backup(1:3,1:3,c,i,e)      = crystallite_Fp(1:3,1:3,c,i,e)
           InvFp_backup(1:3,1:3,c,i,e)   = crystallite_invFp(1:3,1:3,c,i,e)
           Fi_backup(1:3,1:3,c,i,e)      = crystallite_Fi(1:3,1:3,c,i,e)
           InvFi_backup(1:3,1:3,c,i,e)   = crystallite_invFi(1:3,1:3,c,i,e)
           Fe_backup(1:3,1:3,c,i,e)      = crystallite_Fe(1:3,1:3,c,i,e)
           Lp_backup(1:3,1:3,c,i,e)      = crystallite_Lp(1:3,1:3,c,i,e)
           Li_backup(1:3,1:3,c,i,e)      = crystallite_Li(1:3,1:3,c,i,e)
           Tstar_v_backup(1:6,c,i,e)     = crystallite_Tstar_v(1:6,c,i,e)
           P_backup(1:3,1:3,c,i,e)       = crystallite_P(1:3,1:3,c,i,e)
           convergenceFlag_backup(c,i,e) = crystallite_converged(c,i,e)
         enddo; enddo
       enddo elementLooping7
     !$END PARALLEL DO
     ! --- CALCULATE STATE AND STRESS FOR PERTURBATION ---

     dPdF_perturbation1 = crystallite_dPdF0                                                            ! initialize stiffness with known good values from last increment
     dPdF_perturbation2 = crystallite_dPdF0                                                            ! initialize stiffness with known good values from last increment
     pertubationLoop: do perturbation = 1,2                                                                             ! forward and backward perturbation
       if (iand(pert_method,perturbation) > 0_pInt) then                                               ! mask for desired direction
         myPert = -pert_Fg * (-1.0_pReal)**perturbation                                                ! set perturbation step
         do k = 1,3; do l = 1,3                                                                        ! ...alter individual components
           if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
               .and. ((e == debug_e .and. i == debug_i .and. c == debug_g) &
                      .or. .not. iand(debug_level(debug_crystallite),debug_levelSelective) /= 0_pInt)) &
               write(6,'(a,2(1x,i1),1x,a,/)') '<< CRYST >> [[[[[[ Stiffness perturbation',k,l,']]]]]]'
           ! --- INITIALIZE UNPERTURBED STATE ---

           select case(numerics_integrator(numerics_integrationMode))
             case(1_pInt)
!why not OMP?                                                                           ! Fix-point method: restore to last converged state at end of subinc, since this is probably closest to perturbed state
               do e = FEsolving_execElem(1),FEsolving_execElem(2)
                 myNcomponents = homogenization_Ngrains(mesh_element(3,e))
                 do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e); do c = 1,myNcomponents

                   plasticState    (phaseAt(c,i,e))%state(       :,phasememberAt(c,i,e)) = &
                   plasticState    (phaseAt(c,i,e))%state_backup(:,phasememberAt(c,i,e))
                   do mySource = 1_pInt, phase_Nsources(phaseAt(c,i,e))
                     sourceState(phaseAt(c,i,e))%p(mySource)%state(       :,phasememberAt(c,i,e)) = &
                     sourceState(phaseAt(c,i,e))%p(mySource)%state_backup(:,phasememberAt(c,i,e))
                   enddo

                   plasticState    (phaseAt(c,i,e))%dotState(       :,phasememberAt(c,i,e)) = &
                   plasticState    (phaseAt(c,i,e))%dotState_backup(:,phasememberAt(c,i,e))
                   do mySource = 1_pInt, phase_Nsources(phaseAt(c,i,e))
                     sourceState(phaseAt(c,i,e))%p(mySource)%dotState(       :,phasememberAt(c,i,e)) = &
                     sourceState(phaseAt(c,i,e))%p(mySource)%dotState_backup(:,phasememberAt(c,i,e))
                   enddo

                   crystallite_Fp(1:3,1:3,c,i,e)    = Fp_backup(1:3,1:3,c,i,e)
                   crystallite_invFp(1:3,1:3,c,i,e) = InvFp_backup(1:3,1:3,c,i,e)
                   crystallite_Fi(1:3,1:3,c,i,e)    = Fi_backup(1:3,1:3,c,i,e)
                   crystallite_invFi(1:3,1:3,c,i,e) = InvFi_backup(1:3,1:3,c,i,e)
                   crystallite_Fe(1:3,1:3,c,i,e)    = Fe_backup(1:3,1:3,c,i,e)
                   crystallite_Lp(1:3,1:3,c,i,e)    = Lp_backup(1:3,1:3,c,i,e)
                   crystallite_Li(1:3,1:3,c,i,e)    = Li_backup(1:3,1:3,c,i,e)
                   crystallite_Tstar_v(1:6,c,i,e)   = Tstar_v_backup(1:6,c,i,e)
                 enddo; enddo
               enddo
             case(2_pInt,3_pInt)                                                                    ! explicit Euler methods: nothing to restore (except for F), since we are only doing a stress integration step
             case(4_pInt,5_pInt)
!why not OMP?                                                                   ! explicit Runge-Kutta methods: restore to start of subinc, since we are doing a full integration of state and stress
               do e = FEsolving_execElem(1),FEsolving_execElem(2)
                 myNcomponents = homogenization_Ngrains(mesh_element(3,e))
                 do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e); do c = 1,myNcomponents

                   plasticState    (phaseAt(c,i,e))%state(    :,phasememberAt(c,i,e)) = &
                   plasticState    (phaseAt(c,i,e))%subState0(:,phasememberAt(c,i,e))
                   do mySource = 1_pInt, phase_Nsources(phaseAt(c,i,e))
                     sourceState(phaseAt(c,i,e))%p(mySource)%state(    :,phasememberAt(c,i,e)) = &
                     sourceState(phaseAt(c,i,e))%p(mySource)%subState0(:,phasememberAt(c,i,e))
                   enddo

                   plasticState    (phaseAt(c,i,e))%dotState(       :,phasememberAt(c,i,e)) = &
                   plasticState    (phaseAt(c,i,e))%dotState_backup(:,phasememberAt(c,i,e))
                   do mySource = 1_pInt, phase_Nsources(phaseAt(c,i,e))
                     sourceState(phaseAt(c,i,e))%p(mySource)%dotState(       :,phasememberAt(c,i,e)) = &
                     sourceState(phaseAt(c,i,e))%p(mySource)%dotState_backup(:,phasememberAt(c,i,e))
                   enddo

                   crystallite_Fp(1:3,1:3,c,i,e)    = crystallite_subFp0(1:3,1:3,c,i,e)
                   crystallite_Fi(1:3,1:3,c,i,e)    = crystallite_subFi0(1:3,1:3,c,i,e)
                   crystallite_Fe(1:3,1:3,c,i,e)    = crystallite_subFe0(1:3,1:3,c,i,e)
                   crystallite_Lp(1:3,1:3,c,i,e)    = crystallite_subLp0(1:3,1:3,c,i,e)
                   crystallite_Li(1:3,1:3,c,i,e)    = crystallite_subLi0(1:3,1:3,c,i,e)
                   crystallite_Tstar_v(1:6,c,i,e)   = crystallite_subTstar0_v(1:6,c,i,e)
                 enddo; enddo
               enddo
           end select

           ! --- PERTURB EITHER FORWARD OR BACKWARD ---
!why not OMP?
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             myNcomponents = homogenization_Ngrains(mesh_element(3,e))
             do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
               do c = 1,myNcomponents
                 crystallite_subF(1:3,1:3,c,i,e) = F_backup(1:3,1:3,c,i,e)
                 crystallite_subF(k,l,c,i,e) = crystallite_subF(k,l,c,i,e) + myPert
                 crystallite_todo(c,i,e) = crystallite_requested(c,i,e) &
                                          .and. convergenceFlag_backup(c,i,e)
                 if (crystallite_todo(c,i,e)) crystallite_converged(c,i,e) = .false.                ! start out non-converged
           enddo; enddo; enddo


           select case(numerics_integrator(numerics_integrationMode))
             case(1_pInt)
               call crystallite_integrateStateFPI()
             case(2_pInt)
               call crystallite_integrateStateEuler()
             case(3_pInt)
               call crystallite_integrateStateAdaptiveEuler()
             case(4_pInt)
               call crystallite_integrateStateRK4()
             case(5_pInt)
               call crystallite_integrateStateRKCK45()
           end select
           !why not OMP?
           elementLooping8: do e = FEsolving_execElem(1),FEsolving_execElem(2)
             myNcomponents = homogenization_Ngrains(mesh_element(3,e))
             select case(perturbation)
               case(1_pInt)
                 forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), c = 1:myNcomponents, &
                         crystallite_requested(c,i,e) .and. crystallite_converged(c,i,e)) &         ! converged state warrants stiffness update
                   dPdF_perturbation1(1:3,1:3,k,l,c,i,e) = &
                                  (crystallite_P(1:3,1:3,c,i,e) - P_backup(1:3,1:3,c,i,e)) / myPert ! tangent dP_ij/dFg_kl
               case(2_pInt)
                 forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), c = 1:myNcomponents, &
                         crystallite_requested(c,i,e) .and. crystallite_converged(c,i,e)) &         ! converged state warrants stiffness update
                   dPdF_perturbation2(1:3,1:3,k,l,c,i,e) = &
                                  (crystallite_P(1:3,1:3,c,i,e) - P_backup(1:3,1:3,c,i,e)) / myPert ! tangent dP_ij/dFg_kl
               end select
           enddo elementLooping8

         enddo; enddo                                                                               ! k,l component perturbation loop

       endif
     enddo pertubationLoop

     ! --- STIFFNESS ACCORDING TO PERTURBATION METHOD AND CONVERGENCE ---

     elementLooping9: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       myNcomponents = homogenization_Ngrains(mesh_element(3,e))
       select case(pert_method)
         case(1_pInt)
           forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), c = 1:myNcomponents, &
                   crystallite_requested(c,i,e) .and. convergenceFlag_backup(c,i,e)) &                 ! perturbation mode 1: central solution converged
             crystallite_dPdF(1:3,1:3,1:3,1:3,c,i,e) = dPdF_perturbation1(1:3,1:3,1:3,1:3,c,i,e)
         case(2_pInt)
           forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), c = 1:myNcomponents, &
                   crystallite_requested(c,i,e) .and. convergenceFlag_backup(c,i,e)) &                 ! perturbation mode 2: central solution converged
             crystallite_dPdF(1:3,1:3,1:3,1:3,c,i,e) = dPdF_perturbation2(1:3,1:3,1:3,1:3,c,i,e)
         case(3_pInt)
           forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), c = 1:myNcomponents, &
                   crystallite_requested(c,i,e) .and. convergenceFlag_backup(c,i,e)) &                 ! perturbation mode 3: central solution converged
             crystallite_dPdF(1:3,1:3,1:3,1:3,c,i,e) = 0.5_pReal* (  dPdF_perturbation1(1:3,1:3,1:3,1:3,c,i,e) &
                                                                   + dPdF_perturbation2(1:3,1:3,1:3,1:3,c,i,e))
       end select
       forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), c = 1:myNcomponents, &
               crystallite_requested(c,i,e) .and. .not. convergenceFlag_backup(c,i,e)) &               ! for any pertubation mode: if central solution did not converge...
         crystallite_dPdF(1:3,1:3,1:3,1:3,c,i,e) = crystallite_fallbackdPdF(1:3,1:3,1:3,1:3,c,i,e)     ! ...use (elastic) fallback
     enddo elementLooping9

     ! --- RESTORE ---
!why not OMP?
     elementLooping10: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       myNcomponents = homogenization_Ngrains(mesh_element(3,e))
       do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e); do c = 1,myNcomponents

         plasticState    (phaseAt(c,i,e))%state(       :,phasememberAt(c,i,e)) = &
         plasticState    (phaseAt(c,i,e))%state_backup(:,phasememberAt(c,i,e))
         do mySource = 1_pInt, phase_Nsources(phaseAt(c,i,e))
           sourceState(phaseAt(c,i,e))%p(mySource)%state(       :,phasememberAt(c,i,e)) = &
           sourceState(phaseAt(c,i,e))%p(mySource)%state_backup(:,phasememberAt(c,i,e))
         enddo

         plasticState    (phaseAt(c,i,e))%dotState(       :,phasememberAt(c,i,e)) = &
         plasticState    (phaseAt(c,i,e))%dotState_backup(:,phasememberAt(c,i,e))
         do mySource = 1_pInt, phase_Nsources(phaseAt(c,i,e))
           sourceState(phaseAt(c,i,e))%p(mySource)%dotState(       :,phasememberAt(c,i,e)) = &
           sourceState(phaseAt(c,i,e))%p(mySource)%dotState_backup(:,phasememberAt(c,i,e))
         enddo

         crystallite_subF(1:3,1:3,c,i,e)  = F_backup(1:3,1:3,c,i,e)
         crystallite_Fp(1:3,1:3,c,i,e)    = Fp_backup(1:3,1:3,c,i,e)
         crystallite_invFp(1:3,1:3,c,i,e) = InvFp_backup(1:3,1:3,c,i,e)
         crystallite_Fi(1:3,1:3,c,i,e)    = Fi_backup(1:3,1:3,c,i,e)
         crystallite_invFi(1:3,1:3,c,i,e) = InvFi_backup(1:3,1:3,c,i,e)
         crystallite_Fe(1:3,1:3,c,i,e)    = Fe_backup(1:3,1:3,c,i,e)
         crystallite_Lp(1:3,1:3,c,i,e)    = Lp_backup(1:3,1:3,c,i,e)
         crystallite_Li(1:3,1:3,c,i,e)    = Li_backup(1:3,1:3,c,i,e)
         crystallite_Tstar_v(1:6,c,i,e)   = Tstar_v_backup(1:6,c,i,e)
         crystallite_P(1:3,1:3,c,i,e)     = P_backup(1:3,1:3,c,i,e)
         crystallite_converged(c,i,e)     = convergenceFlag_backup(c,i,e)
       enddo; enddo
     enddo elementLooping10

     deallocate(dPdF_perturbation1)
     deallocate(dPdF_perturbation2)
     deallocate(F_backup          )
     deallocate(Fp_backup         )
     deallocate(InvFp_backup      )
     deallocate(Fi_backup         )
     deallocate(InvFi_backup      )
     deallocate(Fe_backup         )
     deallocate(Lp_backup         )
     deallocate(Li_backup         )
     deallocate(P_backup          )
     deallocate(Tstar_v_backup    )
     deallocate(convergenceFlag_backup)

   endif jacobianMethod
 endif computeJacobian
!why not OMP?

end subroutine crystallite_stressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with 4th order explicit Runge Kutta method
!--------------------------------------------------------------------------------------------------
subroutine crystallite_integrateStateRK4()
 use prec, only: &
   prec_isNaN
 use numerics, only: &
   numerics_integrationMode
 use debug, only: &
   debug_level, &
   debug_crystallite, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_g, &
   debug_StateLoopDistribution
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
   material_Nphase, &
   phaseAt, phasememberAt
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
       NaN = any(prec_isNaN(plasticState(p)%dotState(:,c)))
       do mySource = 1_pInt, phase_Nsources(p)
         NaN = NaN .or. any(prec_isNaN(sourceState(p)%p(mySource)%dotState(:,c)))
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

#ifndef _OPENMP
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
         !***dirty way to pass orientation information
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
           NaN = any(prec_isNaN(plasticState(p)%dotState(:,c)))
           do mySource = 1_pInt, phase_Nsources(p)
             NaN = NaN .or. any(prec_isNaN(sourceState(p)%p(mySource)%dotState(:,c)))
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
   if (crystallite_todo(g,i,e)) then
     crystallite_converged(g,i,e) = .true.                                                                ! if still "to do" then converged per definitionem
     if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
       !$OMP CRITICAL (distributionState)
         debug_StateLoopDistribution(4,numerics_integrationMode) = &
           debug_StateLoopDistribution(4,numerics_integrationMode) + 1_pInt
       !$OMP END CRITICAL (distributionState)
     endif
   endif
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
 use prec, only: &
   prec_isNaN
 use debug, only: &
   debug_level, &
   debug_crystallite, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_g, &
   debug_StateLoopDistribution
 use numerics, only: &
   rTol_crystalliteState, &
   numerics_integrationMode
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
       NaN = any(prec_isNaN(plasticState(p)%dotState(:,cc)))
       do mySource = 1_pInt, phase_Nsources(p)
         NaN = NaN .or. any(prec_isNaN(sourceState(p)%p(mySource)%dotState(:,cc)))
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
         !***dirty way to pass orientations to constitutive_microstructure
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
#ifndef _OPENMP
   if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) &
     write(6,'(a,1x,i1)') '<< CRYST >> Runge--Kutta step',stage+1_pInt
#endif
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) &
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), &
                                           crystallite_Fe, &
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
         NaN = any(prec_isNaN(plasticState(p)%dotState(:,cc)))
         do mySource = 1_pInt, phase_Nsources(p)
           NaN = NaN .or. any(prec_isNaN(sourceState(p)%p(mySource)%dotState(:,cc)))
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

#ifndef _OPENMP
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
       !***dirty way to pass orientations to constitutive_microstructure
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
     if (crystallite_todo(g,i,e)) then
       crystallite_converged(g,i,e) = .true.                                                               ! if still "to do" then converged per definition
       if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
         !$OMP CRITICAL (distributionState)
           debug_StateLoopDistribution(6,numerics_integrationMode) = &
             debug_StateLoopDistribution(6,numerics_integrationMode) + 1_pInt
         !$OMP END CRITICAL (distributionState)
       endif
     endif
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
 use prec, only: &
   prec_isNaN
 use debug, only: &
   debug_level, &
   debug_crystallite, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_g, &
   debug_StateLoopDistribution
 use numerics, only: &
   rTol_crystalliteState, &
   numerics_integrationMode
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

 integrationMode: if (numerics_integrationMode == 1_pInt) then

 !$OMP PARALLEL
   ! --- DOT STATE (EULER INTEGRATION) ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) &
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), &
                                           crystallite_Fe, &
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
         NaN = any(prec_isNaN(plasticState(p)%dotState(:,c)))
         do mySource = 1_pInt, phase_Nsources(p)
           NaN = NaN .or. any(prec_isNaN(sourceState(p)%p(mySource)%dotState(:,c)))
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
         !***dirty way to pass orientations to constitutive_microstructure
         call constitutive_microstructure(crystallite_orientation,       &
                                          crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), &
                                          g, i, e)                                                         ! update dependent state variables to be consistent with basic states
     enddo; enddo; enddo
   !$OMP ENDDO
 !$OMP END PARALLEL
 endif integrationMode


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


 if (numerics_integrationMode == 1_pInt) then

   !$OMP PARALLEL
   ! --- DOT STATE (HEUN METHOD) ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) &
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), &
                                           crystallite_Fe, &
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
         NaN = any(prec_isNaN(plasticState(p)%dotState(:,c)))
         do mySource = 1_pInt, phase_Nsources(p)
           NaN = NaN .or. any(prec_isNaN(sourceState(p)%p(mySource)%dotState(:,c)))
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

#ifndef _OPENMP

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
               debug_StateLoopDistribution(2,numerics_integrationMode) = &
                 debug_StateLoopDistribution(2,numerics_integrationMode) + 1_pInt
             !$OMP END CRITICAL (distributionState)
           endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 !$OMP END PARALLEL

 elseif (numerics_integrationMode > 1) then ! stiffness calculation

   !$OMP PARALLEL DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         crystallite_converged(g,i,e) = .true.                                                               ! ... converged per definitionem
         if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
           !$OMP CRITICAL (distributionState)
             debug_StateLoopDistribution(2,numerics_integrationMode) = &
               debug_StateLoopDistribution(2,numerics_integrationMode) + 1_pInt
           !$OMP END CRITICAL (distributionState)
         endif
       endif
     enddo; enddo; enddo
   !$OMP END PARALLEL DO

 endif



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
 use prec, only: &
   prec_isNaN
 use debug, only: &
   debug_level, &
   debug_crystallite, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_g, &
   debug_StateLoopDistribution
 use numerics, only: &
   numerics_integrationMode, &
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

 if (numerics_integrationMode == 1_pInt) then
 !$OMP PARALLEL

   ! --- DOT STATE  ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) &
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), &
                                           crystallite_Fe, &
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
         NaN = any(prec_isNaN(plasticState(p)%dotState(:,c)))
         do mySource = 1_pInt, phase_Nsources(p)
           NaN = NaN .or. any(prec_isNaN(sourceState(p)%p(mySource)%dotState(:,c)))
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

#ifndef _OPENMP
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
         !***dirty way to pass orientations to constitutive_microstructure
         call constitutive_microstructure(crystallite_orientation,       &
                                          crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), &
                                          g, i, e)                                                           ! update dependent state variables to be consistent with basic states
   enddo; enddo; enddo
   !$OMP ENDDO
  !$OMP END PARALLEL
 endif


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
     if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
       crystallite_converged(g,i,e) = .true.                                                               ! if still "to do" then converged per definitionem
       if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
         !$OMP CRITICAL (distributionState)
           debug_StateLoopDistribution(1,numerics_integrationMode) = &
             debug_StateLoopDistribution(1,numerics_integrationMode) + 1_pInt
         !$OMP END CRITICAL (distributionState)
       endif
     endif
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
 use prec, only: &
   prec_isNaN
 use debug, only: &
   debug_e, &
   debug_i, &
   debug_g, &
   debug_level,&
   debug_crystallite, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_StateLoopDistribution
 use numerics, only: &
   nState, &
   numerics_integrationMode, &
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
       NaN = any(prec_isNaN(plasticState(p)%dotState(:,c)))
       do mySource = 1_pInt, phase_Nsources(p)
         NaN = NaN .or. any(prec_isNaN(sourceState(p)%p(mySource)%dotState(:,c)))
       enddo
       if (NaN) then                                                                                       ! NaN occured in any dotState
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
         !***dirty way to pass orientations to constitutive_micrsotructure
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
         NaN = any(prec_isNaN(plasticState(p)%dotState(:,c)))
         do mySource = 1_pInt, phase_Nsources(p)
           NaN = NaN .or. any(prec_isNaN(sourceState(p)%p(mySource)%dotState(:,c)))
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

#ifndef _OPENMP
         if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
             .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                    .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
           write(6,'(a,i8,1x,i2,1x,i3,/)')       '<< CRYST >> update state at el ip g ',e,i,g
           write(6,'(a,f6.1,/)')                 '<< CRYST >> plasticstatedamper ',plasticStatedamper
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> plastic state residuum',plasticStateResiduum(1:mySizePlasticDotState)
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
         if (converged) then
           crystallite_converged(g,i,e) = .true.                                                                   ! ... converged per definition

           if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
             !$OMP CRITICAL (distributionState)
               debug_StateLoopDistribution(NiterationState,numerics_integrationMode) = &
                 debug_StateLoopDistribution(NiterationState,numerics_integrationMode) + 1_pInt
             !$OMP END CRITICAL (distributionState)
           endif
         endif
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
     write(6,'(a,i8,a,i2,/)') '<< CRYST >> ', count(crystallite_converged(:,:,:)), &
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
 use prec, only: &
   prec_isNaN
 use debug, only: &
   debug_level, &
   debug_crystallite, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_g
 use material, only: &
   plasticState, &
   sourceState, &
   phase_Nsources, &
   phaseAt, phasememberAt
 use constitutive, only: &
   constitutive_collectDeltaState

 implicit none
 integer(pInt), intent(in):: &
   el, &                      ! element index
   ip, &                      ! integration point index
   ipc                         ! grain index

 integer(pInt) :: &
   c, &
   p, &
   mySource, &
   mySizePlasticDeltaState, &
   mySizeSourceDeltaState

 c= phasememberAt(ipc,ip,el)
 p = phaseAt(ipc,ip,el)
 call constitutive_collectDeltaState(crystallite_Tstar_v(1:6,ipc,ip,el), crystallite_Fe(1:3,1:3,ipc,ip,el), ipc,ip,el)
 mySizePlasticDeltaState = plasticState(p)%sizeDeltaState
 if( any(prec_isNaN(plasticState(p)%deltaState(:,c)))) then                                         ! NaN occured in deltaState
   crystallite_stateJump = .false.
   return
 endif
 plasticState(p)%state(1:mySizePlasticDeltaState,c) = plasticState(p)%state(1:mySizePlasticDeltaState,c) + &
                                                      plasticState(p)%deltaState(1:mySizePlasticDeltaState,c)
 do mySource = 1_pInt, phase_Nsources(p)
   mySizeSourceDeltaState = sourceState(p)%p(mySource)%sizeDeltaState
   if( any(prec_isNaN(sourceState(p)%p(mySource)%deltaState(:,c)))) then                            ! NaN occured in deltaState
     crystallite_stateJump = .false.
     return
   endif
   sourceState(p)%p(mySource)%state(1:mySizeSourceDeltaState,c) = &
     sourceState(p)%p(mySource)%state(1:mySizeSourceDeltaState,c) + &
     sourceState(p)%p(mySource)%deltaState(1:mySizeSourceDeltaState,c)
 enddo

#ifndef _OPENMP
 if (any(plasticState(p)%deltaState(1:mySizePlasticDeltaState,c) /= 0.0_pReal) &
     .and. iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
     .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
             .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
   write(6,'(a,i8,1x,i2,1x,i3, /)') '<< CRYST >> update state at el ip ipc ',el,ip,ipc
   write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> deltaState', plasticState(p)%deltaState(1:mySizePlasticDeltaState,c)
   write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> new state',  plasticState(p)%state     (1:mySizePlasticDeltaState,c)
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
  math_transpose33, &
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
                   math_transpose33(math_inv33(crystallite_subF(1:3,1:3,ipc,ip,el))))
 crystallite_push33ToRef = math_mul33x33(math_transpose33(T),math_mul33x33(tensor33,T))

end function crystallite_push33ToRef


!--------------------------------------------------------------------------------------------------
!> @brief calculation of stress (P) with time integration based on a residuum in Lp and
!> intermediate acceleration of the Newton-Raphson correction
!--------------------------------------------------------------------------------------------------
logical function crystallite_integrateStress(&
      ipc,&          ! grain number
      ip,&          ! integration point number
      el,&          ! element number
      timeFraction &
      )
 use prec, only:         pLongInt, &
                         tol_math_check, &
                         prec_isNaN
 use numerics, only:     nStress, &
                         aTol_crystalliteStress, &
                         rTol_crystalliteStress, &
                         iJacoLpresiduum, &
                         numerics_integrationMode
 use debug, only:        debug_level, &
                         debug_crystallite, &
                         debug_levelBasic, &
                         debug_levelExtensive, &
                         debug_levelSelective, &
                         debug_e, &
                         debug_i, &
                         debug_g, &
                         debug_cumLpCalls, &
                         debug_cumLpTicks, &
                         debug_StressLoopLpDistribution, &
                         debug_StressLoopLiDistribution
 use constitutive, only: constitutive_LpAndItsTangent, &
                         constitutive_LiAndItsTangent, &
                         constitutive_TandItsTangent
 use math, only:         math_mul33x33, &
                         math_mul33xx33, &
                         math_mul3333xx3333, &
                         math_mul66x6, &
                         math_mul99x99, &
                         math_transpose33, &
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
 use mesh, only:         mesh_element

 implicit none
 integer(pInt), intent(in)::         el, &                          ! element index
                                     ip, &                          ! integration point index
                                     ipc                             ! grain index
 real(pReal), optional, intent(in) :: timeFraction                 ! fraction of timestep

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
 real(pReal), dimension(3,3,3,3)::   dT_dFe3333, &                                                   ! partial derivative of 2nd Piola-Kirchhoff stress
                                     dT_dFi3333, &
                                     dFe_dLp3333, &                                                  ! partial derivative of elastic deformation gradient
                                     dFe_dLi3333, &
                                     dFi_dLi3333, &
                                     dLp_dFi3333, &
                                     dLi_dFi3333, &
                                     dLp_dT3333, &
                                     dLi_dT3333
 real(pReal)                         detInvFi, &                                                     ! determinant of InvFi
                                     steplengthLp0, &
                                     steplengthLp, &
                                     steplengthLi0, &
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
 integer(pLongInt)                   tick, &
                                     tock, &
                                     tickrate, &
                                     maxticks

 external :: &
#if(FLOAT==8)
   dgesv
#elif(FLOAT==4)
   sgesv
#endif

 !* be pessimistic
 crystallite_integrateStress = .false.
#ifndef _OPENMP
 if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
     .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
            .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
   write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> integrateStress at el ip ipc ',el,ip,ipc
 endif
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
 if (all(abs(invFp_current) <= tiny(0.0_pReal))) then                                               ! math_inv33 returns zero when failed, avoid floating point comparison
#ifndef _OPENMP
   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3)') '<< CRYST >> integrateStress failed on inversion of Fp_current at el (elFE) ip g ',&
       el,'(',mesh_element(1,el),')',ip,ipc
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) > 0_pInt) &
       write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fp_current',math_transpose33(Fp_current(1:3,1:3))
   endif
#endif
   return
 endif
 A = math_mul33x33(Fg_new,invFp_current)                                                            ! intermediate tensor needed later to calculate dFe_dLp

 !* inversion of Fi_current...

 invFi_current = math_inv33(Fi_current)
 if (all(abs(invFi_current) <= tiny(0.0_pReal))) then                                               ! math_inv33 returns zero when failed, avoid floating point comparison
#ifndef _OPENMP
   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3)') '<< CRYST >> integrateStress failed on inversion of Fi_current at el (elFE) ip ipc ',&
       el,'(',mesh_element(1,el),')',ip,ipc
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) > 0_pInt) &
       write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fp_current',math_transpose33(Fi_current(1:3,1:3))
   endif
#endif
   return
 endif

 !* start LpLoop with normal step length

 NiterationStressLi = 0_pInt
 jacoCounterLi      = 0_pInt
 steplengthLi0      = 1.0_pReal
 steplengthLi       = steplengthLi0
 residuumLi_old     = 0.0_pReal

 LiLoop: do
   NiterationStressLi = NiterationStressLi + 1_pInt
   IloopsExeced: if (NiterationStressLi > nStress) then
#ifndef _OPENMP
     if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) &
       write(6,'(a,i3,a,i8,1x,a,i8,a,1x,i2,1x,i3,/)') '<< CRYST >> integrateStress reached inelastic loop limit',nStress, &
       ' at el (elFE) ip ipc ', el,mesh_element(1,el),ip,ipc
#endif
     return
   endif IloopsExeced

   invFi_new = math_mul33x33(invFi_current,math_I3 - dt*Liguess)
   Fi_new    = math_inv33(invFi_new)
   detInvFi  = math_det33(invFi_new)

   NiterationStressLp = 0_pInt
   jacoCounterLp      = 0_pInt
   steplengthLp0      = 1.0_pReal
   steplengthLp       = steplengthLp0
   residuumLp_old     = 0.0_pReal
   Lpguess_old        = Lpguess

   LpLoop: do                                    ! inner stress integration loop for consistency with Fi
     NiterationStressLp = NiterationStressLp + 1_pInt
     loopsExeced: if (NiterationStressLp > nStress) then
#ifndef _OPENMP
     if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) &
         write(6,'(a,i3,a,i8,1x,a,i8,a,1x,i2,1x,i3,/)') '<< CRYST >> integrateStress reached loop limit',nStress, &
         ' at el (elFE) ip ipc ', el,mesh_element(1,el),ip,ipc
#endif
       return
     endif loopsExeced

     !* calculate (elastic) 2nd Piola--Kirchhoff stress tensor and its tangent from constitutive law

     B  = math_I3 - dt*Lpguess
     Fe = math_mul33x33(math_mul33x33(A,B), invFi_new)                                                  ! current elastic deformation tensor
     call constitutive_TandItsTangent(Tstar, dT_dFe3333, dT_dFi3333, Fe, Fi_new, ipc, ip, el)               ! call constitutive law to calculate 2nd Piola-Kirchhoff stress and its derivative in unloaded configuration
     Tstar_v = math_Mandel33to6(Tstar)

     !* calculate plastic velocity gradient and its tangent from constitutive law

     if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
       call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
     endif

     call constitutive_LpAndItsTangent(Lp_constitutive, dLp_dT3333, dLp_dFi3333, &
                                       Tstar_v, Fi_new, ipc, ip, el)

     if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
       call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
       !$OMP CRITICAL (debugTimingLpTangent)
         debug_cumLpCalls = debug_cumLpCalls + 1_pInt
         debug_cumLpTicks = debug_cumLpTicks + tock-tick
         !$OMP FLUSH (debug_cumLpTicks)
         if (tock < tick) debug_cumLpTicks = debug_cumLpTicks + maxticks
       !$OMP END CRITICAL (debugTimingLpTangent)
     endif

#ifndef _OPENMP
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
         .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
       write(6,'(a,i3,/)') '<< CRYST >> stress iteration ', NiterationStressLp
       write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Lp_constitutive', math_transpose33(Lp_constitutive)
       write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Lpguess', math_transpose33(Lpguess)
     endif
#endif


     !* update current residuum and check for convergence of loop

     aTolLp = max(rTol_crystalliteStress * max(norm2(Lpguess),norm2(Lp_constitutive)), &             ! absolute tolerance from largest acceptable relative error
                  aTol_crystalliteStress)                                                            ! minimum lower cutoff
     residuumLp = Lpguess - Lp_constitutive

     if (any(prec_isNaN(residuumLp))) then                                                           ! NaN in residuum...
#ifndef _OPENMP
       if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) &
         write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3,a,i3,a)') '<< CRYST >> integrateStress encountered NaN at el (elFE) ip ipc ', &
           el,mesh_element(1,el),ip,ipc, &
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
       steplengthLp   = steplengthLp0                                                                 ! ...proceed with normal step length (calculate new search direction)
     else                                                                                             ! not converged and residuum not improved...
       steplengthLp = 0.5_pReal * steplengthLp                                                        ! ...try with smaller step length in same direction
       Lpguess    = Lpguess_old + steplengthLp * deltaLp
       cycle LpLoop
     endif


     !* calculate Jacobian for correction term

     if (mod(jacoCounterLp, iJacoLpresiduum) == 0_pInt) then
       dFe_dLp3333 = 0.0_pReal
       forall(o=1_pInt:3_pInt,p=1_pInt:3_pInt) &
         dFe_dLp3333(o,1:3,p,1:3) = A(o,p)*math_transpose33(invFi_new)                                ! dFe_dLp(i,j,k,l) = -dt * A(i,k) invFi(l,j)
       dFe_dLp3333 = - dt * dFe_dLp3333
       dRLp_dLp    =   math_identity2nd(9_pInt) &
                     - math_Plain3333to99(math_mul3333xx3333(math_mul3333xx3333(dLp_dT3333,dT_dFe3333),dFe_dLp3333))
       dRLp_dLp2   = dRLp_dLp                                                                         ! will be overwritten in first call to LAPACK routine
       work = math_plain33to9(residuumLp)
#if(FLOAT==8)
       call dgesv(9,1,dRLp_dLp2,9,ipiv,work,9,ierr)                                                   ! solve dRLp/dLp * delta Lp = -res for delta Lp
#elif(FLOAT==4)
       call sgesv(9,1,dRLp_dLp2,9,ipiv,work,9,ierr)                                                   ! solve dRLp/dLp * delta Lp = -res for delta Lp
#endif
       if (ierr /= 0_pInt) then
#ifndef _OPENMP
         if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
           write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3,a,i3)') '<< CRYST >> integrateStress failed on dR/dLp inversion at el ip ipc ', &
             el,mesh_element(1,el),ip,ipc
           if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
               .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g)&
                      .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
             write(6,*)
             write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dR_dLp',transpose(dRLp_dLp)
             write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dFe_dLp',transpose(math_Plain3333to99(dFe_dLp3333))
             write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dT_dFe_constitutive',transpose(math_Plain3333to99(dT_dFe3333))
             write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dLp_dT_constitutive',transpose(math_Plain3333to99(dLp_dT3333))
             write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> A',math_transpose33(A)
             write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> B',math_transpose33(B)
             write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Lp_constitutive',math_transpose33(Lp_constitutive)
             write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Lpguess',math_transpose33(Lpguess)
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

   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     !$OMP CRITICAL (distributionStress)
      debug_StressLoopLpDistribution(NiterationStressLp,numerics_integrationMode) = &
        debug_StressLoopLpDistribution(NiterationStressLp,numerics_integrationMode) + 1_pInt
     !$OMP END CRITICAL (distributionStress)
   endif

   !* calculate intermediate velocity gradient and its tangent from constitutive law

   call constitutive_LiAndItsTangent(Li_constitutive, dLi_dT3333, dLi_dFi3333, &
                                     Tstar_v, Fi_new, ipc, ip, el)

#ifndef _OPENMP
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
         .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
       write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Li_constitutive', math_transpose33(Li_constitutive)
       write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Liguess', math_transpose33(Liguess)
     endif
#endif
   !* update current residuum and check for convergence of loop

   aTolLi = max(rTol_crystalliteStress * max(norm2(Liguess),norm2(Li_constitutive)), &              ! absolute tolerance from largest acceptable relative error
                aTol_crystalliteStress)                                                             ! minimum lower cutoff
   residuumLi = Liguess - Li_constitutive
   if (any(prec_isNaN(residuumLi))) then                                                            ! NaN in residuum...
     return                                                                                         ! ...me = .false. to inform integrator about problem
   elseif (norm2(residuumLi) < aTolLi) then                                                         ! converged if below absolute tolerance
     exit LiLoop                                                                                    ! ...leave iteration loop
   elseif (     NiterationStressLi == 1_pInt &
           .or. norm2(residuumLi) < norm2(residuumLi_old)) then                                     ! not converged, but improved norm of residuum (always proceed in first iteration)...
     residuumLi_old = residuumLi                                                                    ! ...remember old values and...
     Liguess_old    = Liguess
     steplengthLi   = steplengthLi0                                                                 ! ...proceed with normal step length (calculate new search direction)
   else                                                                                             ! not converged and residuum not improved...
     steplengthLi   = 0.5_pReal * steplengthLi                                                      ! ...try with smaller step length in same direction
     Liguess        = Liguess_old + steplengthLi * deltaLi
     cycle LiLoop
   endif

   !* calculate Jacobian for correction term

   if (mod(jacoCounterLi, iJacoLpresiduum) == 0_pInt) then
     temp_33     = math_mul33x33(math_mul33x33(A,B),invFi_current)
     dFe_dLi3333 = 0.0_pReal
     dFi_dLi3333 = 0.0_pReal
     forall(o=1_pInt:3_pInt,p=1_pInt:3_pInt)
       dFe_dLi3333(1:3,o,1:3,p) = -dt*math_I3(o,p)*temp_33                                          ! dFe_dLp(i,j,k,l) = -dt * A(i,k) invFi(l,j)
       dFi_dLi3333(1:3,o,1:3,p) = -dt*math_I3(o,p)*invFi_current
     end forall
     forall(o=1_pInt:3_pInt,p=1_pInt:3_pInt) &
       dFi_dLi3333(1:3,1:3,o,p) = math_mul33x33(math_mul33x33(Fi_new,dFi_dLi3333(1:3,1:3,o,p)),Fi_new)

     dRLi_dLi  = math_identity2nd(9_pInt) &
               - math_Plain3333to99(math_mul3333xx3333(dLi_dT3333, math_mul3333xx3333(dT_dFe3333, dFe_dLi3333) + &
                                                                   math_mul3333xx3333(dT_dFi3333, dFi_dLi3333)))  &
               - math_Plain3333to99(math_mul3333xx3333(dLi_dFi3333, dFi_dLi3333))
     work = math_plain33to9(residuumLi)
#if(FLOAT==8)
     call dgesv(9,1,dRLi_dLi,9,ipiv,work,9,ierr)                                                    ! solve dRLi/dLp * delta Li = -res for delta Li
#elif(FLOAT==4)
     call sgesv(9,1,dRLi_dLi,9,ipiv,work,9,ierr)                                                    ! solve dRLi/dLp * delta Li = -res for delta Li
#endif
       if (ierr /= 0_pInt) then
#ifndef _OPENMP
         if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
           write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3,a,i3)') '<< CRYST >> integrateStress failed on dR/dLi inversion at el ip ipc ', &
             el,mesh_element(1,el),ip,ipc
           if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
               .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g)&
                      .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
             write(6,*)
             write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dR_dLi',transpose(dRLi_dLi)
             write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dFe_dLi',transpose(math_Plain3333to99(dFe_dLi3333))
             write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dT_dFi_constitutive',transpose(math_Plain3333to99(dT_dFi3333))
             write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dLi_dT_constitutive',transpose(math_Plain3333to99(dLi_dT3333))
             write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Li_constitutive',math_transpose33(Li_constitutive)
             write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Liguess',math_transpose33(Liguess)
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

 if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
   !$OMP CRITICAL (distributionStress)
    debug_StressLoopLiDistribution(NiterationStressLi,numerics_integrationMode) = &
      debug_StressLoopLiDistribution(NiterationStressLi,numerics_integrationMode) + 1_pInt
   !$OMP END CRITICAL (distributionStress)
 endif

 !* calculate new plastic and elastic deformation gradient

 invFp_new = math_mul33x33(invFp_current,B)
 invFp_new = invFp_new / math_det33(invFp_new)**(1.0_pReal/3.0_pReal)                               ! regularize by det
 Fp_new = math_inv33(invFp_new)
 if (all(abs(Fp_new)<= tiny(0.0_pReal))) then                                                       ! math_inv33 returns zero when failed, avoid floating point comparison
#ifndef _OPENMP
   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3,a,i3)') '<< CRYST >> integrateStress failed on invFp_new inversion at el ip ipc ',&
       el,mesh_element(1,el),ip,ipc, ' ; iteration ', NiterationStressLp
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
         .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) &
       write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> invFp_new',math_transpose33(invFp_new)
   endif
#endif
   return
 endif
 Fe_new = math_mul33x33(math_mul33x33(Fg_new,invFp_new),invFi_new)    ! calc resulting Fe

 !* calculate 1st Piola-Kirchhoff stress

 crystallite_P(1:3,1:3,ipc,ip,el) = math_mul33x33(math_mul33x33(Fg_new,invFp_new), &
                                              math_mul33x33(math_Mandel6to33(Tstar_v), &
                                                            math_transpose33(invFp_new)))

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
#ifndef _OPENMP
 if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt &
     .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
             .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> P / MPa',math_transpose33(crystallite_P(1:3,1:3,ipc,ip,el))*1.0e-6_pReal
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Cauchy / MPa', &
              math_mul33x33(crystallite_P(1:3,1:3,ipc,ip,el), math_transpose33(Fg_new)) * 1.0e-6_pReal / math_det33(Fg_new)
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fe Lp Fe^-1', &
              math_transpose33(math_mul33x33(Fe_new, math_mul33x33(crystallite_Lp(1:3,1:3,ipc,ip,el), math_inv33(Fe_new))))    ! transpose to get correct print out order
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fp',math_transpose33(crystallite_Fp(1:3,1:3,ipc,ip,el))
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fi',math_transpose33(crystallite_Fi(1:3,1:3,ipc,ip,el))
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
         orientation = math_RtoQ(transpose(math_rotationalPart33(crystallite_Fe(1:3,1:3,c,i,e))))          ! rotational part from polar decomposition as quaternion
         !$OMP END CRITICAL (polarDecomp)
         crystallite_rotation(1:4,c,i,e) = lattice_qDisorientation(crystallite_orientation0(1:4,c,i,e), &  ! active rotation from ori0
                                                                orientation)                               ! to current orientation (with no symmetry)
         crystallite_orientation(1:4,c,i,e) = orientation
  enddo; enddo; enddo
 !$OMP END PARALLEL DO


 ! --- UPDATE SOME ADDITIONAL VARIABLES THAT ARE NEEDED FOR NONLOCAL MATERIAL ---
 ! --- we use crystallite_orientation from above, so need a separate loop

 !$OMP PARALLEL DO PRIVATE(myPhase,neighboring_e,neighboring_i,neighboringPhase)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       myPhase = material_phase(1,i,e)                                                                     ! get my phase (non-local models make no sense with more than one grain per material point)
       if (plasticState(myPhase)%nonLocal) then                                                            ! if nonlocal model
         ! --- calculate disorientation between me and my neighbor ---

         do n = 1_pInt,FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,e))))                        ! loop through my neighbors
           neighboring_e = mesh_ipNeighborhood(1,n,i,e)
           neighboring_i = mesh_ipNeighborhood(2,n,i,e)
           if (neighboring_e > 0 .and. neighboring_i > 0) then                                             ! if neighbor exists
             neighboringPhase = material_phase(1,neighboring_i,neighboring_e)                              ! get my neighbor's phase
             if (plasticState(neighboringPhase)%nonLocal) then                                             ! neighbor got also nonlocal plasticity
               if (lattice_structure(myPhase) == lattice_structure(neighboringPhase)) then                 ! if my neighbor has same crystal structure like me
                 crystallite_disorientation(:,n,1,i,e) = &
                   lattice_qDisorientation( crystallite_orientation(1:4,1,i,e), &
                                            crystallite_orientation(1:4,1,neighboring_i,neighboring_e), &
                                            lattice_structure(myPhase))                                    ! calculate disorientation for given symmetry
               else                                                                                        ! for neighbor with different phase
                 crystallite_disorientation(:,n,1,i,e) = [0.0_pReal, 1.0_pReal, 0.0_pReal, 0.0_pReal]      ! 180 degree rotation about 100 axis
               endif
             else                                                                                          ! for neighbor with local plasticity
               crystallite_disorientation(:,n,1,i,e) = [-1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal]       ! homomorphic identity
             endif
           else                                                                                            ! no existing neighbor
             crystallite_disorientation(:,n,1,i,e) = [-1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal]         ! homomorphic identity
           endif
         enddo


         ! --- calculate compatibility and transmissivity between me and my neighbor ---

         call plastic_nonlocal_updateCompatibility(crystallite_orientation,i,e)

       endif
   enddo; enddo
 !$OMP END PARALLEL DO

end subroutine crystallite_orientations

!--------------------------------------------------------------------------------------------------
!> @brief return results of particular grain
!--------------------------------------------------------------------------------------------------
function crystallite_postResults(ipc, ip, el)
 use math, only: &
   math_qToEuler, &
   math_qToEulerAxisAngle, &
   math_mul33x33, &
   math_transpose33, &
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
                                           / homogenization_Ngrains(mesh_element(3,el))             ! grain volume (not fraction but absolute)
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
         reshape(math_transpose33(crystallite_partionedF(1:3,1:3,ipc,ip,el)),[mySize])
     case (e_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = 0.5_pReal * reshape((math_mul33x33( &
                                               math_transpose33(crystallite_partionedF(1:3,1:3,ipc,ip,el)), &
                                               crystallite_partionedF(1:3,1:3,ipc,ip,el)) - math_I3),[mySize])
     case (fe_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(math_transpose33(crystallite_Fe(1:3,1:3,ipc,ip,el)),[mySize])
     case (ee_ID)
       Ee = 0.5_pReal *(math_mul33x33(math_transpose33(crystallite_Fe(1:3,1:3,ipc,ip,el)), &
                                               crystallite_Fe(1:3,1:3,ipc,ip,el)) - math_I3)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = reshape(Ee,[mySize])
     case (fp_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(math_transpose33(crystallite_Fp(1:3,1:3,ipc,ip,el)),[mySize])
     case (fi_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(math_transpose33(crystallite_Fi(1:3,1:3,ipc,ip,el)),[mySize])
     case (lp_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(math_transpose33(crystallite_Lp(1:3,1:3,ipc,ip,el)),[mySize])
     case (li_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(math_transpose33(crystallite_Li(1:3,1:3,ipc,ip,el)),[mySize])
     case (p_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(math_transpose33(crystallite_P(1:3,1:3,ipc,ip,el)),[mySize])
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
