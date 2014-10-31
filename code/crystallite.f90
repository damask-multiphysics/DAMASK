!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
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
   crystallite_F0, &                                                                                !< def grad at start of FE inc
   crystallite_partionedF,  &                                                                       !< def grad to be reached at end of homog inc
   crystallite_partionedF0, &                                                                       !< def grad at start of homog inc
   crystallite_Lp, &                                                                                !< current plastic velocitiy grad (end of converged time step)
   crystallite_Lp0, &                                                                               !< plastic velocitiy grad at start of FE inc
   crystallite_partionedLp0,&                                                                       !< plastic velocity grad at start of homog inc
   crystallite_P                                                                                    !< 1st Piola-Kirchhoff stress per grain
 real(pReal),                dimension(:,:,:,:,:),    allocatable, private :: &
   crystallite_Fe, &                                                                                !< current "elastic" def grad (end of converged time step)
   crystallite_subFe0,&                                                                             !< "elastic" def grad at start of crystallite inc
   crystallite_invFp, &                                                                             !< inverse of current plastic def grad (end of converged time step)
   crystallite_subFp0,&                                                                             !< plastic def grad at start of crystallite inc
   crystallite_subF,  &                                                                             !< def grad to be reached at end of crystallite inc
   crystallite_subF0, &                                                                             !< def grad at start of crystallite inc
   crystallite_subLp0,&                                                                             !< plastic velocity grad at start of crystallite inc
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
                 lp_ID, &
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
 use lattice, only: &
   lattice_structure
 use constitutive, only: &
   constitutive_microstructure                                                                     ! derived (shortcut) quantities of given state

 implicit none
 integer(pInt), parameter :: &
   FILEUNIT = 200_pInt, &
   MAXNCHUNKS = 2_pInt

 integer(pInt), dimension(1+2*MAXNCHUNKS) :: positions
 integer(pInt) :: &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e, &                                                                                             !< element number
   gMax, &                                                                                          !< maximum number of grains
   iMax, &                                                                                          !< maximum number of integration points
   eMax, &                                                                                          !< maximum number of elements
   nMax, &                                                                                          !< maximum number of ip neighbors
   myNgrains, &                                                                                     !< number of grains in current IP
   section = 0_pInt, &
   j, &
   p, &
   output = 0_pInt, &
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

 gMax = homogenization_maxNgrains
 iMax = mesh_maxNips
 eMax = mesh_NcpElems
 nMax = mesh_maxNipNeighbors


 allocate(crystallite_Tstar0_v(6,gMax,iMax,eMax),            source=0.0_pReal)
 allocate(crystallite_partionedTstar0_v(6,gMax,iMax,eMax),   source=0.0_pReal)
 allocate(crystallite_subTstar0_v(6,gMax,iMax,eMax),         source=0.0_pReal)
 allocate(crystallite_Tstar_v(6,gMax,iMax,eMax),             source=0.0_pReal)
 allocate(crystallite_P(3,3,gMax,iMax,eMax),                 source=0.0_pReal)
 allocate(crystallite_F0(3,3,gMax,iMax,eMax),                source=0.0_pReal)
 allocate(crystallite_partionedF0(3,3,gMax,iMax,eMax),       source=0.0_pReal)
 allocate(crystallite_partionedF(3,3,gMax,iMax,eMax),        source=0.0_pReal)
 allocate(crystallite_subF0(3,3,gMax,iMax,eMax),             source=0.0_pReal)
 allocate(crystallite_subF(3,3,gMax,iMax,eMax),              source=0.0_pReal)
 allocate(crystallite_Fp0(3,3,gMax,iMax,eMax),               source=0.0_pReal)
 allocate(crystallite_partionedFp0(3,3,gMax,iMax,eMax),      source=0.0_pReal)
 allocate(crystallite_subFp0(3,3,gMax,iMax,eMax),            source=0.0_pReal)
 allocate(crystallite_Fp(3,3,gMax,iMax,eMax),                source=0.0_pReal)
 allocate(crystallite_invFp(3,3,gMax,iMax,eMax),             source=0.0_pReal)
 allocate(crystallite_Fe(3,3,gMax,iMax,eMax),                source=0.0_pReal)
 allocate(crystallite_subFe0(3,3,gMax,iMax,eMax),            source=0.0_pReal)
 allocate(crystallite_Lp0(3,3,gMax,iMax,eMax),               source=0.0_pReal)
 allocate(crystallite_partionedLp0(3,3,gMax,iMax,eMax),      source=0.0_pReal)
 allocate(crystallite_subLp0(3,3,gMax,iMax,eMax),            source=0.0_pReal)
 allocate(crystallite_Lp(3,3,gMax,iMax,eMax),                source=0.0_pReal)
 allocate(crystallite_dPdF(3,3,3,3,gMax,iMax,eMax),          source=0.0_pReal)
 allocate(crystallite_dPdF0(3,3,3,3,gMax,iMax,eMax),         source=0.0_pReal)
 allocate(crystallite_partioneddPdF0(3,3,3,3,gMax,iMax,eMax),source=0.0_pReal)
 allocate(crystallite_fallbackdPdF(3,3,3,3,gMax,iMax,eMax),  source=0.0_pReal)
 allocate(crystallite_dt(gMax,iMax,eMax),                    source=0.0_pReal)
 allocate(crystallite_subdt(gMax,iMax,eMax),                 source=0.0_pReal)
 allocate(crystallite_subFrac(gMax,iMax,eMax),               source=0.0_pReal)
 allocate(crystallite_subStep(gMax,iMax,eMax),               source=0.0_pReal)
 allocate(crystallite_orientation(4,gMax,iMax,eMax),         source=0.0_pReal)
 allocate(crystallite_orientation0(4,gMax,iMax,eMax),        source=0.0_pReal)
 allocate(crystallite_rotation(4,gMax,iMax,eMax),            source=0.0_pReal)
 allocate(crystallite_disorientation(4,nMax,gMax,iMax,eMax), source=0.0_pReal)
 allocate(crystallite_localPlasticity(gMax,iMax,eMax),       source=.true.)
 allocate(crystallite_requested(gMax,iMax,eMax),             source=.false.)
 allocate(crystallite_todo(gMax,iMax,eMax),                  source=.false.)
 allocate(crystallite_converged(gMax,iMax,eMax),             source=.true.)
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
 rewind(FILEUNIT)
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
     output = 0_pInt                                                                                ! reset output counter
     cycle                                                                                          ! skip to next line
   endif
   if (section > 0_pInt) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         output = output + 1_pInt
         crystallite_output(output,section) = IO_lc(IO_stringValue(line,positions,2_pInt))
         select case(crystallite_output(output,section))
           case ('phase')
             crystallite_outputID(output,section) = phase_ID
           case ('texture')
             crystallite_outputID(output,section) = texture_ID
           case ('volume')
             crystallite_outputID(output,section) = volume_ID
           case ('grainrotationx')
             crystallite_outputID(output,section) = grainrotationx_ID
           case ('grainrotationy')
             crystallite_outputID(output,section) = grainrotationy_ID
           case ('grainrotationz')
             crystallite_outputID(output,section) = grainrotationx_ID
           case ('orientation')
             crystallite_outputID(output,section) = orientation_ID
           case ('grainrotation')
             crystallite_outputID(output,section) = grainrotation_ID
           case ('eulerangles')
             crystallite_outputID(output,section) = eulerangles_ID
           case ('defgrad','f')
             crystallite_outputID(output,section) = defgrad_ID
           case ('fe')
             crystallite_outputID(output,section) = fe_ID
           case ('fp')
             crystallite_outputID(output,section) = fp_ID
           case ('lp')
             crystallite_outputID(output,section) = lp_ID
           case ('e')
             crystallite_outputID(output,section) = e_ID
           case ('ee')
             crystallite_outputID(output,section) = ee_ID
           case ('p','firstpiola','1piola')
             crystallite_outputID(output,section) = p_ID
           case ('s','tstar','secondpiola','2ndpiola')
             crystallite_outputID(output,section) = s_ID
           case ('elasmatrix')
             crystallite_outputID(output,section) = elasmatrix_ID
           case ('neighboringip')
             crystallite_outputID(output,section) = neighboringip_ID
           case ('neighboringelement')
             crystallite_outputID(output,section) = neighboringelement_ID
           case default
             call IO_error(105_pInt,ext_msg=IO_stringValue(line,positions,2_pInt)//' (Crystallite)')
         end select
     end select
   endif
 enddo

 close(FILEUNIT)

 do i = 1_pInt,material_Ncrystallite
   do j = 1_pInt,crystallite_Noutput(i)
     select case(crystallite_outputID(j,i))
       case(phase_ID,texture_ID,volume_ID,grainrotationx_ID,grainrotationy_ID,grainrotationz_ID)
         mySize = 1_pInt
       case(orientation_ID,grainrotation_ID)                                                          ! orientation as quaternion, or deviation from initial grain orientation in axis-angle form (angle in degrees)
         mySize = 4_pInt
       case(eulerangles_ID)
         mySize = 3_pInt
       case(defgrad_ID,fe_ID,fp_ID,lp_ID,e_ID,ee_ID,p_ID,s_ID)
         mySize = 9_pInt
       case(elasmatrix_ID)
         mySize = 36_pInt
       case(neighboringip_ID,neighboringelement_ID)
         mySize = mesh_maxNipNeighbors
       case default
         mySize = 0_pInt
     end select

     outputFound: if (mySize > 0_pInt) then
       crystallite_sizePostResult(j,i) = mySize
       crystallite_sizePostResults(i) = crystallite_sizePostResults(i) + mySize
     endif outputFound
   enddo
 enddo

 crystallite_maxSizePostResults = 0_pInt
 do j = 1_pInt,material_Nmicrostructure
   if (microstructure_active(j)) &
     crystallite_maxSizePostResults = max(crystallite_maxSizePostResults,&
                                          crystallite_sizePostResults(microstructure_crystallite(j)))
 enddo

!--------------------------------------------------------------------------------------------------
! write description file for crystallite output
 call IO_write_jobFile(FILEUNIT,'outputCrystallite')

 do p = 1_pInt,material_Ncrystallite
   write(FILEUNIT,'(/,a,/)') '['//trim(crystallite_name(p))//']'
   do e = 1_pInt,crystallite_Noutput(p)
     write(FILEUNIT,'(a,i4)') trim(crystallite_output(e,p))//char(9),crystallite_sizePostResult(e,p)
   enddo
 enddo

 close(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! initialize
!$OMP PARALLEL DO PRIVATE(myNgrains)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                 ! iterate over all cp elements
     myNgrains = homogenization_Ngrains(mesh_element(3,e))                                            ! look up homogenization-->grainCount
     forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1_pInt:myNgrains)
       crystallite_Fp0(1:3,1:3,g,i,e) = math_EulerToR(material_EulerAngles(1:3,g,i,e))                ! plastic def gradient reflects init orientation
       crystallite_F0(1:3,1:3,g,i,e)  = math_I3
       crystallite_localPlasticity(g,i,e) = phase_localPlasticity(material_phase(g,i,e))
       crystallite_Fe(1:3,1:3,g,i,e)  = math_transpose33(crystallite_Fp0(1:3,1:3,g,i,e))
       crystallite_Fp(1:3,1:3,g,i,e)  = crystallite_Fp0(1:3,1:3,g,i,e)
       crystallite_requested(g,i,e) = .true.
     endforall
   enddo
 !$OMP END PARALLEL DO

 if(any(.not. crystallite_localPlasticity) .and. .not. usePingPong) call IO_error(601_pInt)         ! exit if nonlocal but no ping-pong

 crystallite_partionedFp0 = crystallite_Fp0
 crystallite_partionedF0 = crystallite_F0
 crystallite_partionedF = crystallite_F0

 call crystallite_orientations()
 crystallite_orientation0 = crystallite_orientation                                                 ! store initial orientations for calculation of grain rotations


 !$OMP PARALLEL DO PRIVATE(myNgrains)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     myNgrains = homogenization_Ngrains(mesh_element(3,e))
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1_pInt,myNgrains
         call constitutive_microstructure( &
                          crystallite_Tstar_v(1:6,g,i,e), &
                          crystallite_Fe(1:3,1:3,g,i,e), &
                          crystallite_Fp(1:3,1:3,g,i,e),g,i,e)                                      ! update dependent state variables to be consistent with basic states
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
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Lp:                    ', shape(crystallite_Lp)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_F0:                    ', shape(crystallite_F0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Fp0:                   ', shape(crystallite_Fp0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_Lp0:                   ', shape(crystallite_Lp0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedF:            ', shape(crystallite_partionedF)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedF0:           ', shape(crystallite_partionedF0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedFp0:          ', shape(crystallite_partionedFp0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedLp0:          ', shape(crystallite_partionedLp0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subF:                  ', shape(crystallite_subF)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subF0:                 ', shape(crystallite_subF0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subFe0:                ', shape(crystallite_subFe0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subFp0:                ', shape(crystallite_subFp0)
   write(6,'(a35,1x,7(i8,1x))') 'crystallite_subLp0:                ', shape(crystallite_subLp0)
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
   relevantStrain, &
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
   damageState, &
   thermalState, &
   mappingConstitutive, &
   homogenization_maxNgrains
 use constitutive, only:  &
   constitutive_TandItsTangent, &
   constitutive_LpAndItsTangent, &
   constitutive_LiAndItsTangent, &
   constitutive_getFi, &
   constitutive_getFi0, &
   constitutive_getPartionedFi0

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
 real(pReal), dimension(3,3,3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   dPdF_perturbation1, &
   dPdF_perturbation2
 real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   F_backup, &
   Fp_backup, &
   InvFp_backup, &
   Fe_backup, &
   Lp_backup, &
   P_backup
 real(pReal), dimension(6,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   Tstar_v_backup
 integer(pInt) :: &
   NiterationCrystallite, &                                                                         ! number of iterations in crystallite loop
   e, &                                                                                             ! element index
   i, &                                                                                             ! integration point index
   g, &                                                                                             ! grain index
   k, &
   l, &
   n, startIP, endIP, &
   neighboring_e, &
   neighboring_i, &
   o, &
   p, &
   perturbation , &                                                                                 ! loop counter for forward,backward perturbation mode
   myNgrains
 logical, dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                                       convergenceFlag_backup
 ! local variables used for calculating analytic Jacobian
 real(pReal)                     ::   detInvFi
 real(pReal), dimension(3,3)     ::   junk, &
                                      Fi, &
                                      invFi, &
                                      Fi0, &
                                      invFi0
 real(pReal), dimension(3,3,3,3) ::   dSdFe, &
                                      dSdF, &
                                      junk2, &
                                      dLidS, &
                                      dLpdS, &
                                      dFpinvdF, &
                                      rhs_3333, &
                                      lhs_3333, &
                                      temp_3333
 real(pReal), dimension(9,9)::       temp_99
 logical :: error


 if (iand(debug_level(debug_crystallite),debug_levelSelective) /= 0_pInt &
     .and. FEsolving_execElem(1) <= debug_e &
     .and.                          debug_e <= FEsolving_execElem(2)) then
     write(6,'(/,a,i8,1x,a,i8,a,1x,i2,1x,i3)')      '<< CRYST >> values at el (elFE) ip g ', &
       debug_e,'(',mesh_element(1,debug_e), ')',debug_i, debug_g
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> F0 ', &
                                         math_transpose33(crystallite_partionedF0(1:3,1:3,debug_g,debug_i,debug_e))
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Fp0', &
                                         math_transpose33(crystallite_partionedFp0(1:3,1:3,debug_g,debug_i,debug_e))
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Lp0', &
                                         math_transpose33(crystallite_partionedLp0(1:3,1:3,debug_g,debug_i,debug_e))
   write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> F  ', &
                                         math_transpose33(crystallite_partionedF(1:3,1:3,debug_g,debug_i,debug_e))
 endif

!--------------------------------------------------------------------------------------------------
! initialize to starting condition
 crystallite_subStep = 0.0_pReal

 !$OMP PARALLEL DO PRIVATE(myNgrains)
   elementLooping1: do e = FEsolving_execElem(1),FEsolving_execElem(2)
     myNgrains = homogenization_Ngrains(mesh_element(3,e))
     forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             g = 1_pInt:myNgrains, crystallite_requested(g,i,e))
       plasticState(mappingConstitutive(2,g,i,e))%subState0(      :,mappingConstitutive(1,g,i,e)) = &
       plasticState(mappingConstitutive(2,g,i,e))%partionedState0(:,mappingConstitutive(1,g,i,e))
       damageState( mappingConstitutive(2,g,i,e))%subState0(      :,mappingConstitutive(1,g,i,e)) = &
       damageState( mappingConstitutive(2,g,i,e))%partionedState0(:,mappingConstitutive(1,g,i,e))
       thermalState(mappingConstitutive(2,g,i,e))%subState0(      :,mappingConstitutive(1,g,i,e)) = &
       thermalState(mappingConstitutive(2,g,i,e))%partionedState0(:,mappingConstitutive(1,g,i,e))
       crystallite_subFp0(1:3,1:3,g,i,e) = crystallite_partionedFp0(1:3,1:3,g,i,e)                  ! ...plastic def grad
       crystallite_subLp0(1:3,1:3,g,i,e) = crystallite_partionedLp0(1:3,1:3,g,i,e)                  ! ...plastic velocity grad
       crystallite_dPdF0(1:3,1:3,1:3,1:3,g,i,e) = crystallite_partioneddPdF0(1:3,1:3,1:3,1:3,g,i,e) ! ...stiffness
       crystallite_subF0(1:3,1:3,g,i,e) = crystallite_partionedF0(1:3,1:3,g,i,e)                    ! ...def grad
       crystallite_subTstar0_v(1:6,g,i,e) = crystallite_partionedTstar0_v(1:6,g,i,e)                !...2nd PK stress
       crystallite_subFe0(1:3,1:3,g,i,e) = math_mul33x33(math_mul33x33(crystallite_subF0(1:3,1:3,g,i,e), &
                                             math_inv33(crystallite_subFp0(1:3,1:3,g,i,e))), &
                                             math_inv33(constitutive_getFi0(g,i,e)))                ! only needed later on for stiffness calculation
       crystallite_subFrac(g,i,e) = 0.0_pReal
       crystallite_subStep(g,i,e) = 1.0_pReal/subStepSizeCryst
       crystallite_todo(g,i,e) = .true.
       crystallite_converged(g,i,e) = .false.                                                       ! pretend failed step of twice the required size
     endforall
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
   timeSyncing1: if (any(.not. crystallite_localPlasticity) .and. numerics_timeSyncing) then

    ! Time synchronization can only be used for nonlocal calculations, and only there it makes sense.
    ! The idea is that in nonlocal calculations often the vast majority of the ips
    ! converges in one iteration whereas a small fraction of ips has to do a lot of cutbacks.
    ! Hence, we try to minimize the computational effort by just doing a lot of cutbacks
    ! in the vicinity of the "bad" ips and leave the easily converged volume more or less as it is.
    ! However, some synchronization of the time step has to be done at the border between "bad" ips
    ! and the ones that immediately converged.

     if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) &
       write(6,'(a,i6)') '<< CRYST >> crystallite iteration ',NiterationCrystallite

     if (any(crystallite_syncSubFrac)) then

       ! Just did a time synchronization.
       ! If all synchronizers converged, then do nothing else than winding them forward.
       ! If any of the synchronizers did not converge, something went completely wrong
       ! and its not clear how to fix this, so all nonlocals become terminally ill.

       if (any(crystallite_syncSubFrac .and. .not. crystallite_converged(1,:,:))) then
         if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             myNgrains = homogenization_Ngrains(mesh_element(3,e))
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
         !$OMP PARALLEL DO PRIVATE(myNgrains)
         do e = FEsolving_execElem(1),FEsolving_execElem(2)
           myNgrains = homogenization_Ngrains(mesh_element(3,e))
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
         myNgrains = homogenization_Ngrains(mesh_element(3,e))
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

       !$OMP PARALLEL DO PRIVATE(myNgrains)
       do e = FEsolving_execElem(1),FEsolving_execElem(2)
         myNgrains = homogenization_Ngrains(mesh_element(3,e))
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
           !$OMP PARALLEL DO PRIVATE(myNgrains)
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             myNgrains = homogenization_Ngrains(mesh_element(3,e))
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
         if (subFracIntermediate == 0.0_pReal) then
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
           !$OMP DO PRIVATE(myNgrains)
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             myNgrains = homogenization_Ngrains(mesh_element(3,e))
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
               if (.not. crystallite_localPlasticity(1,i,e) .and. crystallite_subFrac(1,i,e) == 0.0_pReal) then
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
           !$OMP DO PRIVATE(myNgrains)
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             myNgrains = homogenization_Ngrains(mesh_element(3,e))
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
           !$OMP PARALLEL DO PRIVATE(myNgrains)
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             myNgrains = homogenization_Ngrains(mesh_element(3,e))
             do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
               crystallite_clearToWindForward(i,e) = crystallite_localPlasticity(1,i,e)
               crystallite_clearToCutback(i,e) = crystallite_localPlasticity(1,i,e) .or. crystallite_syncSubFrac(i,e)
             enddo
           enddo
           !$OMP END PARALLEL DO
           if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) &
             write(6,'(a,i6)') '<< CRYST >> time synchronization: cutback'
         else
           !$OMP PARALLEL DO PRIVATE(myNgrains)
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             myNgrains = homogenization_Ngrains(mesh_element(3,e))
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

     !$OMP PARALLEL DO PRIVATE(myNgrains)
     elementLooping2: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       myNgrains = homogenization_Ngrains(mesh_element(3,e))
       do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
         if(.not. crystallite_clearToWindForward(i,e) .and. .not. crystallite_clearToCutback(i,e)) &
           crystallite_todo(1,i,e) = .false.
       enddo
     enddo elementLooping2
     !$OMP END PARALLEL DO

   endif timeSyncing1
   !$OMP PARALLEL DO PRIVATE(myNgrains,formerSubStep)
     elementLooping3: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       myNgrains = homogenization_Ngrains(mesh_element(3,e))
       do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                            ! iterate over IPs of this element to be processed
         do g = 1,myNgrains
           ! --- wind forward ---

           if (crystallite_converged(g,i,e) .and. crystallite_clearToWindForward(i,e)) then
             formerSubStep = crystallite_subStep(g,i,e)
             crystallite_subFrac(g,i,e) = crystallite_subFrac(g,i,e) + crystallite_subStep(g,i,e)
             !$OMP FLUSH(crystallite_subFrac)
             crystallite_subStep(g,i,e) = min(1.0_pReal - crystallite_subFrac(g,i,e), &
                                              stepIncreaseCryst * crystallite_subStep(g,i,e))
             !$OMP FLUSH(crystallite_subStep)
             if (crystallite_subStep(g,i,e) > 0.0_pReal) then
               crystallite_subF0(1:3,1:3,g,i,e) = crystallite_subF(1:3,1:3,g,i,e)                    ! ...def grad
               !$OMP FLUSH(crystallite_subF0)
               crystallite_subLp0(1:3,1:3,g,i,e) = crystallite_Lp(1:3,1:3,g,i,e)                     ! ...plastic velocity gradient
               crystallite_subFp0(1:3,1:3,g,i,e) = crystallite_Fp(1:3,1:3,g,i,e)                     ! ...plastic def grad
               crystallite_subFe0(1:3,1:3,g,i,e) = math_mul33x33(math_mul33x33(crystallite_subF (1:3,1:3,g,i,e), &
                                                                               crystallite_invFp(1:3,1:3,g,i,e)), &
                                                                 math_inv33(constitutive_getFi(g,i,e)))  ! only needed later on for stiffness calculation
               !if abbrevation, make c and p private in omp
               plasticState(mappingConstitutive(2,g,i,e))%subState0(:,mappingConstitutive(1,g,i,e)) = &
               plasticState(mappingConstitutive(2,g,i,e))%state(    :,mappingConstitutive(1,g,i,e))
               damageState( mappingConstitutive(2,g,i,e))%subState0(:,mappingConstitutive(1,g,i,e)) = &
               damageState( mappingConstitutive(2,g,i,e))%state(    :,mappingConstitutive(1,g,i,e))
               thermalState(mappingConstitutive(2,g,i,e))%subState0(:,mappingConstitutive(1,g,i,e)) = &
               thermalState(mappingConstitutive(2,g,i,e))%state(    :,mappingConstitutive(1,g,i,e))
               crystallite_subTstar0_v(1:6,g,i,e) = crystallite_Tstar_v(1:6,g,i,e)                   ! ...2nd PK stress
               if (crystallite_syncSubFrac(i,e)) then                                                ! if we just did a synchronization of states, then we wind forward without any further time integration
                 crystallite_syncSubFracCompleted(i,e) = .true.
                 crystallite_syncSubFrac(i,e) = .false.
                 crystallite_todo(g,i,e) = .false.
               else
                 crystallite_todo(g,i,e) = .true.
               endif
               !$OMP FLUSH(crystallite_todo)
#ifndef _OPENMP
               if (iand(debug_level(debug_crystallite),debug_levelBasic) /= 0_pInt &
                   .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                          .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) &
                 write(6,'(a,f12.8,a,f12.8,a,i8,1x,i2,1x,i3,/)') '<< CRYST >> winding forward from ', &
                   crystallite_subFrac(g,i,e)-formerSubStep,' to current crystallite_subfrac ', &
                   crystallite_subFrac(g,i,e),' in crystallite_stressAndItsTangent at el ip g ',e,i,g
#endif
             else                                                                                    ! this crystallite just converged for the entire timestep
               crystallite_todo(g,i,e) = .false.                                                     ! so done here
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

           elseif (.not. crystallite_converged(g,i,e) .and. crystallite_clearToCutback(i,e)) then
             if (crystallite_syncSubFrac(i,e)) then                                                  ! synchronize time
               crystallite_subStep(g,i,e) = subFracIntermediate
             else
               crystallite_subStep(g,i,e) = subStepSizeCryst * crystallite_subStep(g,i,e)            ! cut step in half and restore...
             endif
             !$OMP FLUSH(crystallite_subStep)
             crystallite_Fp(1:3,1:3,g,i,e) = crystallite_subFp0(1:3,1:3,g,i,e)                       ! ...plastic def grad
             !$OMP FLUSH(crystallite_Fp)
             crystallite_invFp(1:3,1:3,g,i,e) = math_inv33(crystallite_Fp(1:3,1:3,g,i,e))
             !$OMP FLUSH(crystallite_invFp)
             crystallite_Lp(1:3,1:3,g,i,e)    = crystallite_subLp0(1:3,1:3,g,i,e)                    ! ...plastic velocity grad
             plasticState(mappingConstitutive(2,g,i,e))%state(    :,mappingConstitutive(1,g,i,e)) = &
             plasticState(mappingConstitutive(2,g,i,e))%subState0(:,mappingConstitutive(1,g,i,e))
             damageState( mappingConstitutive(2,g,i,e))%state(    :,mappingConstitutive(1,g,i,e)) = &
             damageState( mappingConstitutive(2,g,i,e))%subState0(:,mappingConstitutive(1,g,i,e))
             thermalState(mappingConstitutive(2,g,i,e))%state(    :,mappingConstitutive(1,g,i,e)) = &
             thermalState(mappingConstitutive(2,g,i,e))%subState0(:,mappingConstitutive(1,g,i,e))
             crystallite_Tstar_v(1:6,g,i,e)   = crystallite_subTstar0_v(1:6,g,i,e)                   ! ...2nd PK stress

                                                                                                     ! cant restore dotState here, since not yet calculated in first cutback after initialization
             crystallite_todo(g,i,e) = crystallite_subStep(g,i,e) > subStepMinCryst                  ! still on track or already done (beyond repair)
             !$OMP FLUSH(crystallite_todo)
#ifndef _OPENMP
             if (iand(debug_level(debug_crystallite),debug_levelBasic) /= 0_pInt) then
               if (crystallite_todo(g,i,e)) then
                 write(6,'(a,f12.8,a,i8,1x,i2,1x,i3,/)') '<< CRYST >> cutback step in crystallite_stressAndItsTangent &
                                                        &with new crystallite_subStep: ',&
                                                       crystallite_subStep(g,i,e),' at el ip g ',e,i,g
               else
                 write(6,'(a,i8,1x,i2,1x,i3,/)') '<< CRYST >> reached minimum step size &
                                               &in crystallite_stressAndItsTangent at el ip g ',e,i,g
               endif
             endif
#endif
           endif

           ! --- prepare for integration ---

           if (crystallite_todo(g,i,e) .and. (crystallite_clearToWindForward(i,e) .or. crystallite_clearToCutback(i,e))) then
             crystallite_subF(1:3,1:3,g,i,e) = crystallite_subF0(1:3,1:3,g,i,e) &
                                             + crystallite_subStep(g,i,e) &
                                               * (crystallite_partionedF(1:3,1:3,g,i,e) &
                                             - crystallite_partionedF0(1:3,1:3,g,i,e))
             !$OMP FLUSH(crystallite_subF)
             crystallite_Fe(1:3,1:3,g,i,e) = math_mul33x33(crystallite_subF(1:3,1:3,g,i,e), crystallite_invFp(1:3,1:3,g,i,e))
             crystallite_Fe(1:3,1:3,g,i,e) = math_mul33x33(math_mul33x33(crystallite_subF (1:3,1:3,g,i,e), &
                                                                         crystallite_invFp(1:3,1:3,g,i,e)), &
                                                           math_inv33(constitutive_getFi(g,i,e)))
             crystallite_subdt(g,i,e) = crystallite_subStep(g,i,e) * crystallite_dt(g,i,e)
             crystallite_converged(g,i,e) = .false.                                                  ! start out non-converged
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
           myNgrains = homogenization_Ngrains(mesh_element(3,e))
           do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
             do g = 1,myNgrains
               if (.not. crystallite_localPlasticity(g,i,e) .and. .not. crystallite_todo(g,i,e) &
                   .and. .not. crystallite_converged(g,i,e) .and. crystallite_subStep(g,i,e) <= subStepMinCryst) &
                 write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> nonlocal violated minimum subStep at el,ip,g ',e,i,g
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
   myNgrains = homogenization_Ngrains(mesh_element(3,e))
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                  ! iterate over IPs of this element to be processed
     do g = 1,myNgrains
       if (.not. crystallite_converged(g,i,e)) then                                                    ! respond fully elastically (might be not required due to becoming terminally ill anyway)
         if(iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) &
           write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3,/)') '<< CRYST >> no convergence: respond fully elastic at el (elFE) ip g ', &
             e,'(',mesh_element(1,e),')',i,g
         invFp = math_inv33(crystallite_partionedFp0(1:3,1:3,g,i,e))
         Fe_guess = math_mul33x33(math_mul33x33(crystallite_partionedF(1:3,1:3,g,i,e), invFp), &
                                  math_inv33(constitutive_getPartionedFi0(g,i,e)))
         call constitutive_TandItsTangent(Tstar, junk2, Fe_guess,g,i,e)
         crystallite_P(1:3,1:3,g,i,e) = math_mul33x33(math_mul33x33(Fe_guess,constitutive_getPartionedFi0(g,i,e)), &
                                                      math_mul33x33(Tstar,transpose(invFp)))* &
                                        math_det33(constitutive_getPartionedFi0(g,i,e))              
       endif
       if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
           .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                  .or. .not. iand(debug_level(debug_crystallite),debug_levelSelective) /= 0_pInt)) then
         write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> central solution of cryst_StressAndTangent at el ip g ',e,i,g
         write(6,'(/,a,/,3(12x,3(f12.4,1x)/))') '<< CRYST >> P / MPa', &
                                          math_transpose33(crystallite_P(1:3,1:3,g,i,e))*1.0e-6_pReal
         write(6,'(a,/,3(12x,3(f14.9,1x)/))')   '<< CRYST >> Fp', &
                                          math_transpose33(crystallite_Fp(1:3,1:3,g,i,e))
         write(6,'(a,/,3(12x,3(f14.9,1x)/),/)') '<< CRYST >> Lp', &
                                          math_transpose33(crystallite_Lp(1:3,1:3,g,i,e))
         flush(6)
       endif
     enddo
   enddo
 enddo elementLooping5


! --+>> STIFFNESS CALCULATION <<+--

 computeJacobian: if(updateJaco) then
   jacobianMethod: if (analyticJaco) then

     ! --- ANALYTIC JACOBIAN ---

     !$OMP PARALLEL DO PRIVATE(dFedF,dSdF,dSdFe,dLpdS,dFpinvdF,rhs_3333,lhs_3333,temp_99,junk,dLidS,Fi,invFi,Fi0,invFi0,detInvFi,temp_3333,myNgrains)
       elementLooping6: do e = FEsolving_execElem(1),FEsolving_execElem(2)
         myNgrains = homogenization_Ngrains(mesh_element(3,e))
         do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                            ! iterate over IPs of this element to be processed
           do g = 1_pInt,myNgrains
             call constitutive_TandItsTangent(junk,dSdFe,crystallite_Fe(1:3,1:3,g,i,e),g,i,e)      ! call constitutive law to calculate 2nd Piola-Kirchhoff stress and its derivative
             call constitutive_LpAndItsTangent(junk,temp_99,crystallite_Tstar_v(1:6,g,i,e), &
                                               g,i,e)
             dLpdS = reshape(temp_99,shape=[3,3,3,3])
             call constitutive_LiAndItsTangent(junk,temp_99,crystallite_Tstar_v(1:6,g,i,e), &
                                               crystallite_Lp(1:3,1:3,g,i,e), g,i,e)
             dLidS    = reshape(temp_99,shape=[3,3,3,3])
             Fi       = constitutive_getFi(g,i,e)
             invFi    = math_inv33(Fi)
             Fi0      = constitutive_getFi0(g,i,e)
             invFi0   = math_inv33(Fi0)
             detInvFi = math_det33(invFi)
             
             junk   = math_transpose33(math_mul33x33(crystallite_invFp(1:3,1:3,g,i,e),invFi))  
             rhs_3333 = 0.0_pReal
             forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
               rhs_3333(p,o,1:3,1:3) = math_mul33x33(dSdFe(p,o,1:3,1:3),junk)
             
             temp_3333 = 0.0_pReal
             junk = math_mul33x33(crystallite_subF(1:3,1:3,g,i,e), &
                                  math_inv33(crystallite_subFp0(1:3,1:3,g,i,e)))
             forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
               temp_3333(1:3,1:3,p,o) = math_mul33x33(math_mul33x33(junk,dLpdS(1:3,1:3,p,o)),invFi)
             
             junk = math_mul33x33(math_mul33x33(crystallite_subF(1:3,1:3,g,i,e), &
                                                crystallite_invFp(1:3,1:3,g,i,e)), invFi0)
             forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
               temp_3333(1:3,1:3,p,o) = temp_3333(1:3,1:3,p,o) + math_mul33x33(junk,dLidS(1:3,1:3,p,o))
             
             lhs_3333 = crystallite_subdt(g,i,e)*math_mul3333xx3333(dSdFe,temp_3333)

             call math_invert(9_pInt,math_identity2nd(9_pInt)+reshape(lhs_3333,shape=[9,9]),temp_99,error)
             if (error) call IO_error(error_ID=400_pInt,ext_msg='analytic tangent inversion')
             dSdF = math_mul3333xx3333(reshape(temp_99,shape=[3,3,3,3]),rhs_3333)
             
             dFpinvdF = 0.0_pReal
             temp_3333 = math_mul3333xx3333(dLpdS,dSdF)
             forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
               dFpinvdF(1:3,1:3,p,o) = -crystallite_subdt(g,i,e)* &
                                        math_mul33x33(math_inv33(crystallite_subFp0(1:3,1:3,g,i,e)), &
                                                      math_mul33x33(temp_3333(1:3,1:3,p,o),invFi))
             
             crystallite_dPdF(1:3,1:3,1:3,1:3,g,i,e) = 0.0_pReal
             junk = math_mul33x33(crystallite_invFp(1:3,1:3,g,i,e), &
                                  math_mul33x33(math_Mandel6to33(crystallite_Tstar_v(1:6,g,i,e)), &
                                                math_transpose33(crystallite_invFp(1:3,1:3,g,i,e))))/detInvFi
             forall(p=1_pInt:3_pInt) &
               crystallite_dPdF(p,1:3,p,1:3,g,i,e) = math_transpose33(junk)

             junk = math_mul33x33(math_Mandel6to33(crystallite_Tstar_v(1:6,g,i,e)), &
                                  math_transpose33(crystallite_invFp(1:3,1:3,g,i,e)))/detInvFi
             forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
               crystallite_dPdF(1:3,1:3,p,o,g,i,e) = crystallite_dPdF(1:3,1:3,p,o,g,i,e) + &
                 math_mul33x33(math_mul33x33(crystallite_subF(1:3,1:3,g,i,e),dFpinvdF(1:3,1:3,p,o)),junk)

             junk = math_mul33x33(crystallite_subF(1:3,1:3,g,i,e), &
                                  crystallite_invFp(1:3,1:3,g,i,e))/detInvFi
             forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
               crystallite_dPdF(1:3,1:3,p,o,g,i,e) = crystallite_dPdF(1:3,1:3,p,o,g,i,e) + &
                 math_mul33x33(math_mul33x33(junk,dSdF(1:3,1:3,p,o)), &
                               math_transpose33(crystallite_invFp(1:3,1:3,g,i,e)))

             junk = math_mul33x33(math_mul33x33(crystallite_subF(1:3,1:3,g,i,e), &
                                                crystallite_invFp(1:3,1:3,g,i,e)), &
                                  math_Mandel6to33(crystallite_Tstar_v(1:6,g,i,e)))/detInvFi
             forall(p=1_pInt:3_pInt, o=1_pInt:3_pInt) &
               crystallite_dPdF(1:3,1:3,p,o,g,i,e) = crystallite_dPdF(1:3,1:3,p,o,g,i,e) + &
                 math_mul33x33(junk,math_transpose33(dFpinvdF(1:3,1:3,p,o)))

         enddo; enddo
       enddo elementLooping6
     !$OMP END PARALLEL DO

   else jacobianMethod

     ! --- STANDARD (PERTURBATION METHOD) FOR JACOBIAN ---

     numerics_integrationMode = 2_pInt

     ! --- BACKUP ---
     !$OMP PARALLEL DO PRIVATE(myNgrains)
       elementLooping7: do e = FEsolving_execElem(1),FEsolving_execElem(2)
         myNgrains = homogenization_Ngrains(mesh_element(3,e))
         forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains)

           plasticState(mappingConstitutive(2,g,i,e))%state_backup(:,mappingConstitutive(1,g,i,e)) = &
           plasticState(mappingConstitutive(2,g,i,e))%state(       :,mappingConstitutive(1,g,i,e))
           damageState( mappingConstitutive(2,g,i,e))%state_backup(:,mappingConstitutive(1,g,i,e)) = &
           damageState( mappingConstitutive(2,g,i,e))%state(       :,mappingConstitutive(1,g,i,e))
           thermalState(mappingConstitutive(2,g,i,e))%state_backup(:,mappingConstitutive(1,g,i,e)) = &
           thermalState(mappingConstitutive(2,g,i,e))%state(       :,mappingConstitutive(1,g,i,e))

           plasticState(mappingConstitutive(2,g,i,e))%dotState_backup(:,mappingConstitutive(1,g,i,e)) = &
           plasticState(mappingConstitutive(2,g,i,e))%dotState(       :,mappingConstitutive(1,g,i,e))
           damageState( mappingConstitutive(2,g,i,e))%dotState_backup(:,mappingConstitutive(1,g,i,e)) = &
           damageState( mappingConstitutive(2,g,i,e))%dotState(       :,mappingConstitutive(1,g,i,e))
           thermalState(mappingConstitutive(2,g,i,e))%dotState_backup(:,mappingConstitutive(1,g,i,e)) = &
           thermalState(mappingConstitutive(2,g,i,e))%dotState(       :,mappingConstitutive(1,g,i,e))

           F_backup(1:3,1:3,g,i,e)       = crystallite_subF(1:3,1:3,g,i,e)                            ! ... and kinematics
           Fp_backup(1:3,1:3,g,i,e)      = crystallite_Fp(1:3,1:3,g,i,e)
           InvFp_backup(1:3,1:3,g,i,e)   = crystallite_invFp(1:3,1:3,g,i,e)
           Fe_backup(1:3,1:3,g,i,e)      = crystallite_Fe(1:3,1:3,g,i,e)
           Lp_backup(1:3,1:3,g,i,e)      = crystallite_Lp(1:3,1:3,g,i,e)
           Tstar_v_backup(1:6,g,i,e)     = crystallite_Tstar_v(1:6,g,i,e)
           P_backup(1:3,1:3,g,i,e)       = crystallite_P(1:3,1:3,g,i,e)
           convergenceFlag_backup(g,i,e) = crystallite_converged(g,i,e)
         endforall
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
               .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                      .or. .not. iand(debug_level(debug_crystallite),debug_levelSelective) /= 0_pInt)) &
               write(6,'(a,2(1x,i1),1x,a,/)') '<< CRYST >> [[[[[[ Stiffness perturbation',k,l,']]]]]]'
           ! --- INITIALIZE UNPERTURBED STATE ---

           select case(numerics_integrator(numerics_integrationMode))
             case(1_pInt)
!why not OMP?                                                                           ! Fix-point method: restore to last converged state at end of subinc, since this is probably closest to perturbed state
               do e = FEsolving_execElem(1),FEsolving_execElem(2)
                 myNgrains = homogenization_Ngrains(mesh_element(3,e))
                 forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains)

                   plasticState(mappingConstitutive(2,g,i,e))%state(       :,mappingConstitutive(1,g,i,e)) = &
                   plasticState(mappingConstitutive(2,g,i,e))%state_backup(:,mappingConstitutive(1,g,i,e))  
                   damageState( mappingConstitutive(2,g,i,e))%state(       :,mappingConstitutive(1,g,i,e)) = &
                   damageState( mappingConstitutive(2,g,i,e))%state_backup(:,mappingConstitutive(1,g,i,e))  
                   thermalState(mappingConstitutive(2,g,i,e))%state(       :,mappingConstitutive(1,g,i,e)) = &
                   thermalState(mappingConstitutive(2,g,i,e))%state_backup(:,mappingConstitutive(1,g,i,e))  

                   plasticState(mappingConstitutive(2,g,i,e))%dotState(       :,mappingConstitutive(1,g,i,e)) = &
                   plasticState(mappingConstitutive(2,g,i,e))%dotState_backup(:,mappingConstitutive(1,g,i,e))
                   damageState( mappingConstitutive(2,g,i,e))%dotState(       :,mappingConstitutive(1,g,i,e)) = &
                   damageState( mappingConstitutive(2,g,i,e))%dotState_backup(:,mappingConstitutive(1,g,i,e))
                   thermalState(mappingConstitutive(2,g,i,e))%dotState(       :,mappingConstitutive(1,g,i,e)) = &
                   thermalState(mappingConstitutive(2,g,i,e))%dotState_backup(:,mappingConstitutive(1,g,i,e))

                   crystallite_Fp(1:3,1:3,g,i,e)    = Fp_backup(1:3,1:3,g,i,e) 
                   crystallite_invFp(1:3,1:3,g,i,e) = InvFp_backup(1:3,1:3,g,i,e)
                   crystallite_Fe(1:3,1:3,g,i,e)    = Fe_backup(1:3,1:3,g,i,e)
                   crystallite_Lp(1:3,1:3,g,i,e)    = Lp_backup(1:3,1:3,g,i,e)
                   crystallite_Tstar_v(1:6,g,i,e)   = Tstar_v_backup(1:6,g,i,e)
                 endforall
               enddo
             case(2_pInt,3_pInt)                                                                    ! explicit Euler methods: nothing to restore (except for F), since we are only doing a stress integration step
             case(4_pInt,5_pInt)
!why not OMP?                                                                   ! explicit Runge-Kutta methods: restore to start of subinc, since we are doing a full integration of state and stress
               do e = FEsolving_execElem(1),FEsolving_execElem(2)
                 myNgrains = homogenization_Ngrains(mesh_element(3,e))
                 forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains)

                   plasticState(mappingConstitutive(2,g,i,e))%state(    :,mappingConstitutive(1,g,i,e)) = &
                   plasticState(mappingConstitutive(2,g,i,e))%subState0(:,mappingConstitutive(1,g,i,e))
                   damageState( mappingConstitutive(2,g,i,e))%state(    :,mappingConstitutive(1,g,i,e)) = &
                   damageState( mappingConstitutive(2,g,i,e))%subState0(:,mappingConstitutive(1,g,i,e))
                   thermalState(mappingConstitutive(2,g,i,e))%state(    :,mappingConstitutive(1,g,i,e)) = &
                   thermalState(mappingConstitutive(2,g,i,e))%subState0(:,mappingConstitutive(1,g,i,e))

                   plasticState(mappingConstitutive(2,g,i,e))%dotState(       :,mappingConstitutive(1,g,i,e)) = &
                   plasticState(mappingConstitutive(2,g,i,e))%dotState_backup(:,mappingConstitutive(1,g,i,e))
                   damageState( mappingConstitutive(2,g,i,e))%dotState(       :,mappingConstitutive(1,g,i,e)) = &
                   damageState( mappingConstitutive(2,g,i,e))%dotState_backup(:,mappingConstitutive(1,g,i,e))
                   thermalState(mappingConstitutive(2,g,i,e))%dotState(       :,mappingConstitutive(1,g,i,e)) = &
                   thermalState(mappingConstitutive(2,g,i,e))%dotState_backup(:,mappingConstitutive(1,g,i,e))

                   crystallite_Fp(1:3,1:3,g,i,e)    = crystallite_subFp0(1:3,1:3,g,i,e) 
                   crystallite_Fe(1:3,1:3,g,i,e)    = crystallite_subFe0(1:3,1:3,g,i,e)
                   crystallite_Lp(1:3,1:3,g,i,e)    = crystallite_subLp0(1:3,1:3,g,i,e)
                   crystallite_Tstar_v(1:6,g,i,e)   = crystallite_subTstar0_v(1:6,g,i,e)
                 endforall
               enddo
           end select

           ! --- PERTURB EITHER FORWARD OR BACKWARD ---
!why not OMP?
           do e = FEsolving_execElem(1),FEsolving_execElem(2)
             myNgrains = homogenization_Ngrains(mesh_element(3,e))
             do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
               do g = 1,myNgrains
                 crystallite_subF(1:3,1:3,g,i,e) = F_backup(1:3,1:3,g,i,e)
                 crystallite_subF(k,l,g,i,e) = crystallite_subF(k,l,g,i,e) + myPert
                 crystallite_todo(g,i,e) = crystallite_requested(g,i,e) &
                                          .and. convergenceFlag_backup(g,i,e)
                 if (crystallite_todo(g,i,e)) crystallite_converged(g,i,e) = .false.                ! start out non-converged
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
             myNgrains = homogenization_Ngrains(mesh_element(3,e))
             select case(perturbation)
               case(1_pInt)
                 forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains, &
                         crystallite_requested(g,i,e) .and. crystallite_converged(g,i,e)) &         ! converged state warrants stiffness update
                   dPdF_perturbation1(1:3,1:3,k,l,g,i,e) = &
                                  (crystallite_P(1:3,1:3,g,i,e) - P_backup(1:3,1:3,g,i,e)) / myPert ! tangent dP_ij/dFg_kl
               case(2_pInt)
                 forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains, &
                         crystallite_requested(g,i,e) .and. crystallite_converged(g,i,e)) &         ! converged state warrants stiffness update
                   dPdF_perturbation2(1:3,1:3,k,l,g,i,e) = &
                                  (crystallite_P(1:3,1:3,g,i,e) - P_backup(1:3,1:3,g,i,e)) / myPert ! tangent dP_ij/dFg_kl
               end select
           enddo elementLooping8

         enddo; enddo                                                                               ! k,l component perturbation loop

       endif
     enddo pertubationLoop

     ! --- STIFFNESS ACCORDING TO PERTURBATION METHOD AND CONVERGENCE ---

     elementLooping9: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       myNgrains = homogenization_Ngrains(mesh_element(3,e))
       select case(pert_method)
         case(1_pInt)
           forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains, &
                   crystallite_requested(g,i,e) .and. convergenceFlag_backup(g,i,e)) &                 ! perturbation mode 1: central solution converged
             crystallite_dPdF(1:3,1:3,1:3,1:3,g,i,e) = dPdF_perturbation1(1:3,1:3,1:3,1:3,g,i,e)
         case(2_pInt)
           forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains, &
                   crystallite_requested(g,i,e) .and. convergenceFlag_backup(g,i,e)) &                 ! perturbation mode 2: central solution converged
             crystallite_dPdF(1:3,1:3,1:3,1:3,g,i,e) = dPdF_perturbation2(1:3,1:3,1:3,1:3,g,i,e)
         case(3_pInt)
           forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains, &
                   crystallite_requested(g,i,e) .and. convergenceFlag_backup(g,i,e)) &                 ! perturbation mode 3: central solution converged
             crystallite_dPdF(1:3,1:3,1:3,1:3,g,i,e) = 0.5_pReal* (  dPdF_perturbation1(1:3,1:3,1:3,1:3,g,i,e) &
                                                                   + dPdF_perturbation2(1:3,1:3,1:3,1:3,g,i,e))
       end select
       forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains, &
               crystallite_requested(g,i,e) .and. .not. convergenceFlag_backup(g,i,e)) &               ! for any pertubation mode: if central solution did not converge...
         crystallite_dPdF(1:3,1:3,1:3,1:3,g,i,e) = crystallite_fallbackdPdF(1:3,1:3,1:3,1:3,g,i,e)     ! ...use (elastic) fallback
     enddo elementLooping9

     ! --- RESTORE ---
!why not OMP?
     elementLooping10: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       myNgrains = homogenization_Ngrains(mesh_element(3,e))
       forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains)

         plasticState(mappingConstitutive(2,g,i,e))%state(       :,mappingConstitutive(1,g,i,e)) = &
         plasticState(mappingConstitutive(2,g,i,e))%state_backup(:,mappingConstitutive(1,g,i,e))
         damageState( mappingConstitutive(2,g,i,e))%state(       :,mappingConstitutive(1,g,i,e)) = &
         damageState( mappingConstitutive(2,g,i,e))%state_backup(:,mappingConstitutive(1,g,i,e))
         thermalState(mappingConstitutive(2,g,i,e))%state(       :,mappingConstitutive(1,g,i,e)) = &
         thermalState(mappingConstitutive(2,g,i,e))%state_backup(:,mappingConstitutive(1,g,i,e))

         plasticState(mappingConstitutive(2,g,i,e))%dotState(       :,mappingConstitutive(1,g,i,e)) = &
         plasticState(mappingConstitutive(2,g,i,e))%dotState_backup(:,mappingConstitutive(1,g,i,e))
         damageState( mappingConstitutive(2,g,i,e))%dotState(       :,mappingConstitutive(1,g,i,e)) = &
         damageState( mappingConstitutive(2,g,i,e))%dotState_backup(:,mappingConstitutive(1,g,i,e))
         thermalState(mappingConstitutive(2,g,i,e))%dotState(       :,mappingConstitutive(1,g,i,e)) = &
         thermalState(mappingConstitutive(2,g,i,e))%dotState_backup(:,mappingConstitutive(1,g,i,e))

         crystallite_subF(1:3,1:3,g,i,e)  = F_backup(1:3,1:3,g,i,e)
         crystallite_Fp(1:3,1:3,g,i,e)    = Fp_backup(1:3,1:3,g,i,e)
         crystallite_invFp(1:3,1:3,g,i,e) = InvFp_backup(1:3,1:3,g,i,e)
         crystallite_Fe(1:3,1:3,g,i,e)    = Fe_backup(1:3,1:3,g,i,e)
         crystallite_Lp(1:3,1:3,g,i,e)    = Lp_backup(1:3,1:3,g,i,e)
         crystallite_Tstar_v(1:6,g,i,e)   = Tstar_v_backup(1:6,g,i,e)
         crystallite_P(1:3,1:3,g,i,e)     = P_backup(1:3,1:3,g,i,e)
         crystallite_converged(g,i,e)     = convergenceFlag_backup(g,i,e)
       endforall
     enddo elementLooping10

   endif jacobianMethod
 endif computeJacobian
!why not OMP?

end subroutine crystallite_stressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with 4th order explicit Runge Kutta method
!--------------------------------------------------------------------------------------------------
subroutine crystallite_integrateStateRK4()
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
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_Ngrains, &
   plasticState, &
   damageState, &
   thermalState, &
   mappingConstitutive, &
   homogenization_maxNgrains
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
                                               mySizePlasticDotState, &
                                               mySizeDamageDotState, &
                                               mySizeThermalDotState
 integer(pInt), dimension(2) ::                eIter                                                 ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) ::  iIter, &                                              ! bounds for ip iteration
                                               gIter                                                 ! bounds for grain iteration
 logical   ::                                  singleRun                                             ! flag indicating computation for single (g,i,e) triple

 eIter = FEsolving_execElem(1:2)
 do e = eIter(1),eIter(2)
   iIter(1:2,e) = FEsolving_execIP(1:2,e)
   gIter(1:2,e) = [ 1_pInt,homogenization_Ngrains(mesh_element(3,e))]
 enddo

 singleRun = (eIter(1) == eIter(2) .and. iIter(1,eIter(1)) == iIter(2,eIter(2)))

!--------------------------------------------------------------------------------------------------
! initialize dotState
 if (.not. singleRun) then
   forall(p = 1_pInt:size(plasticState)) plasticState(p)%RK4dotState = 0.0_pReal
   forall(p = 1_pInt:size(damageState))  damageState(p)%RK4dotState  = 0.0_pReal
   forall(p = 1_pInt:size(thermalState)) thermalState(p)%RK4dotState = 0.0_pReal
 else
   e = eIter(1)
   i = iIter(1,e)
   do g = iIter(1,e), iIter(2,e)
     plasticState(mappingConstitutive(2,g,i,e))%RK4dotState(:,mappingConstitutive(1,g,i,e)) = 0.0_pReal  
     damageState( mappingConstitutive(2,g,i,e))%RK4dotState(:,mappingConstitutive(1,g,i,e)) = 0.0_pReal  
     thermalState(mappingConstitutive(2,g,i,e))%RK4dotState(:,mappingConstitutive(1,g,i,e)) = 0.0_pReal   
   enddo
 endif

!--------------------------------------------------------------------------------------------------
! first Runge-Kutta step
 !$OMP PARALLEL
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                 ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) &
       call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Lp(1:3,1:3,g,i,e), &
                                         crystallite_Fe, &
                                         crystallite_Fp, &
                                         crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
   enddo; enddo; enddo
 !$OMP ENDDO

 !$OMP DO PRIVATE(p,c)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                 ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       c = mappingConstitutive(1,g,i,e) 
       p = mappingConstitutive(2,g,i,e) 
       if ( any(plasticState(p)%dotState(:,c) /= plasticState(p)%dotState(:,c)) .or.&
            any(damageState(p)%dotState(:,c)  /= damageState(p)%dotState(:,c))  .or.&
            any(thermalState(p)%dotState(:,c) /= thermalState(p)%dotState(:,c))) then                   ! NaN occured in dotState
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
         p = mappingConstitutive(2,g,i,e)
         c = mappingConstitutive(1,g,i,e)

         plasticState(p)%RK4dotState(:,c) = plasticState(p)%RK4dotState(:,c) &
                                          + weight(n)*plasticState(p)%dotState(:,c) 
         damageState(p)%RK4dotState(:,c)  = damageState(p)%RK4dotState(:,c) &
                                          + weight(n)*damageState(p)%dotState(:,c) 
         thermalState(p)%RK4dotState(:,c) = thermalState(p)%RK4dotState(:,c) &
                                          + weight(n)*thermalState(p)%dotState(:,c) 
       endif
     enddo; enddo; enddo
   !$OMP ENDDO

   !$OMP DO PRIVATE(mySizePlasticDotState,mySizeDamageDotState,mySizeThermalDotState,p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then

         p = mappingConstitutive(2,g,i,e)
         c = mappingConstitutive(1,g,i,e)
         mySizePlasticDotState = plasticState(p)%sizeDotState
         mySizeDamageDotState  = damageState(p)%sizeDotState
         mySizeThermalDotState = thermalState(p)%sizeDotState
         plasticState(p)%state(1:mySizePlasticDotState,c) = plasticState(p)%subState0(1:mySizePlasticDotState,c) &
                                                          + plasticState(p)%dotState (1:mySizePlasticDotState,c) &
                                                          * crystallite_subdt(g,i,e) * timeStepFraction(n)
         damageState( p)%state(1:mySizeDamageDotState,c)  = damageState(p)%subState0(1:mySizeDamageDotState,c) &
                                                          + damageState(p)%dotState (1:mySizeDamageDotState,c) &
                                                          * crystallite_subdt(g,i,e) * timeStepFraction(n)
         thermalState(p)%state(1:mySizeThermalDotState,c) = thermalState(p)%subState0(1:mySizeThermalDotState,c) &
                                                          + thermalState(p)%dotState (1:mySizeThermalDotState,c) &
                                                          * crystallite_subdt(g,i,e) * timeStepFraction(n)

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
         call constitutive_microstructure(crystallite_Tstar_v(1:6,g,i,e), &
                                          crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), g, i, e)                          ! update dependent state variables to be consistent with basic states
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
           call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Lp(1:3,1:3,g,i,e), &
                                             crystallite_Fe, &
                                             crystallite_Fp, &
                                             timeStepFraction(n)*crystallite_subdt(g,i,e), &               ! fraction of original timestep
                                             crystallite_subFrac, g,i,e)
       enddo; enddo; enddo
     !$OMP ENDDO

     !$OMP DO PRIVATE(p,c)
       do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                ! iterate over elements, ips and grains
         !$OMP FLUSH(crystallite_todo)
         if (crystallite_todo(g,i,e)) then

           p = mappingConstitutive(2,g,i,e) 
           c = mappingConstitutive(1,g,i,e)  
           if ( any(plasticState(p)%dotState(:,c) /= plasticState(p)%dotState(:,c)) .or.&
                any(damageState(p)%dotState(:,c)  /=  damageState(p)%dotState(:,c)) .or.&
                any(thermalState(p)%dotState(:,c) /= thermalState(p)%dotState(:,c))) then             ! NaN occured in dotState
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

 do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
   if (crystallite_todo(g,i,e)) then
     crystallite_converged(g,i,e) = .true.                                                               ! if still "to do" then converged per definitionem
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
   damageState, &
   thermalState, &
   mappingConstitutive, &
   homogenization_maxNgrains
 use constitutive, only: &
   constitutive_collectDotState, &
   constitutive_maxSizeDotState, &
   constitutive_damage_maxSizeDotState, &
   constitutive_thermal_maxSizeDotState, &
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
   mySizePlasticDotState, &                                                                         ! size of dot States
   mySizeDamageDotState, &
   mySizeThermalDotState
 integer(pInt), dimension(2) :: &
   eIter                                                                                            ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) :: &
   iIter, &                                                                                         ! bounds for ip iteration
   gIter                                                                                            ! bounds for grain iteration

 real(pReal), dimension(constitutive_maxSizeDotState,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   stateResiduum, &                                                                                 ! residuum from evolution in microstructure
   relStateResiduum                                                                                 ! relative residuum from evolution in microstructure
 real(pReal), dimension(constitutive_damage_maxSizeDotState,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   damageStateResiduum, &                                                                           ! residuum from evolution in microstructure
   relDamageStateResiduum                                                                           ! relative residuum from evolution in microstructure
 real(pReal), dimension(constitutive_thermal_maxSizeDotState,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   thermalStateResiduum, &                                                                          ! residuum from evolution in microstructure
   relThermalStateResiduum                                                                          ! relative residuum from evolution in microstructure
 logical :: &
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
       call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Lp(1:3,1:3,g,i,e), &
                                         crystallite_Fe, &
                                         crystallite_Fp, &
                                         crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
   enddo; enddo; enddo
 !$OMP ENDDO
 !$OMP DO PRIVATE(p,cc)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       cc = mappingConstitutive(1,g,i,e) 
       p = mappingConstitutive(2,g,i,e) 
       if ( any(plasticState(p)%dotState(:,cc) /= plasticState(p)%dotState(:,cc)) .or.&
            any(damageState(p)%dotState(:,cc)  /= damageState(p)%dotState(:,cc))  .or.&
            any(thermalState(p)%dotState(:,cc) /= thermalState(p)%dotState(:,cc))) then                    ! NaN occured in dotState
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
         p = mappingConstitutive(2,g,i,e)
         cc = mappingConstitutive(1,g,i,e)
         plasticState(p)%RKCK45dotState(stage,:,cc) = plasticState(p)%dotState(:,cc)                       ! store Runge-Kutta dotState
         damageState(p)%RKCK45dotState(stage,:,cc)  = damageState(p)%dotState(:,cc)                        ! store Runge-Kutta dotState
         thermalState(p)%RKCK45dotState(stage,:,cc) = thermalState(p)%dotState(:,cc)                       ! store Runge-Kutta dotState
       endif
     enddo; enddo; enddo
   !$OMP ENDDO

   !$OMP DO PRIVATE(p,cc,n)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         p = mappingConstitutive(2,g,i,e)
         cc = mappingConstitutive(1,g,i,e)

         plasticState(p)%dotState(:,cc) = A(1,stage) * plasticState(p)%RKCK45dotState(1,:,cc)
         damageState( p)%dotState(:,cc) = A(1,stage) * damageState( p)%RKCK45dotState(1,:,cc)
         thermalState(p)%dotState(:,cc) = A(1,stage) * thermalState(p)%RKCK45dotState(1,:,cc)
         do n = 2_pInt, stage
           plasticState(p)%dotState(:,cc) = &
           plasticState(p)%dotState(:,cc) + A(n,stage) * plasticState(p)%RKCK45dotState(n,:,cc)
           damageState( p)%dotState(:,cc) = &
           damageState( p)%dotState(:,cc) + A(n,stage) * damageState( p)%RKCK45dotState(n,:,cc)
           thermalState(p)%dotState(:,cc) = &
           thermalState(p)%dotState(:,cc) + A(n,stage) * thermalState(p)%RKCK45dotState(n,:,cc)
         enddo
       endif
     enddo; enddo; enddo
   !$OMP ENDDO

   !$OMP DO PRIVATE(mySizePlasticDotState,mySizeDamageDotState,mySizeThermalDotState,p,cc)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         p = mappingConstitutive(2,g,i,e)
         cc = mappingConstitutive(1,g,i,e)
         mySizePlasticDotState = plasticState(p)%sizeDotState
         mySizeDamageDotState  = damageState( p)%sizeDotState
         mySizeThermalDotState = thermalState(p)%sizeDotState
         plasticState(p)%state(1:mySizePlasticDotState,cc) = plasticState(p)%subState0(1:mySizePlasticDotState,cc) &
                                                           + plasticState(p)%dotState (1:mySizePlasticDotState,cc) &
                                                           * crystallite_subdt(g,i,e)
         damageState(p)%state(1:mySizeDamageDotState,cc)   = damageState(p)%subState0( 1:mySizeDamageDotState,cc) &
                                                           + damageState(p)%dotState ( 1:mySizeDamageDotState,cc) &
                                                           * crystallite_subdt(g,i,e)
         thermalState(p)%state(1:mySizeThermalDotState,cc) = thermalState(p)%subState0(1:mySizeThermalDotState,cc) &
                                                           + thermalState(p)%dotState (1:mySizeThermalDotState,cc) &
                                                           * crystallite_subdt(g,i,e)
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
         call constitutive_microstructure(crystallite_Tstar_v(1:6,g,i,e), &
                                          crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), g, i, e)                           ! update dependent state variables to be consistent with basic states
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
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Lp(1:3,1:3,g,i,e), &
                                           crystallite_Fe, &
                                           crystallite_Fp, &
                                           C(stage)*crystallite_subdt(g,i,e), & ! fraction of original timestep
                                           crystallite_subFrac, g,i,e)
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO PRIVATE(p,cc)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then

         p = mappingConstitutive(2,g,i,e) 
         cc = mappingConstitutive(1,g,i,e)  
       if ( any(plasticState(p)%dotState(:,cc) /= plasticState(p)%dotState(:,cc)) .or.&
            any(damageState(p)%dotState(:,cc)  /= damageState(p)%dotState(:,cc))  .or.&
            any(thermalState(p)%dotState(:,cc) /= thermalState(p)%dotState(:,cc))) then                    ! NaN occured in dotState
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

 relStateResiduum = 0.0_pReal
 relDamageStateResiduum  = 0.0_pReal
 relThermalStateResiduum = 0.0_pReal
 !$OMP PARALLEL
 !$OMP DO PRIVATE(p,cc)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       p = mappingConstitutive(2,g,i,e) 
       cc = mappingConstitutive(1,g,i,e)  
       plasticState(p)%RKCK45dotState(6,:,cc) = plasticState(p)%dotState(:,cc)            ! store Runge-Kutta dotState
       damageState( p)%RKCK45dotState(6,:,cc) = damageState( p)%dotState(:,cc)            ! store Runge-Kutta dotState
       thermalState(p)%RKCK45dotState(6,:,cc) = thermalState(p)%dotState(:,cc)            ! store Runge-Kutta dotState
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
       
 !$OMP DO PRIVATE(mySizePlasticDotState,mySizeDamageDotState,mySizeThermalDotState,p,cc)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       p = mappingConstitutive(2,g,i,e)
       cc = mappingConstitutive(1,g,i,e)
       mySizePlasticDotState = plasticState(p)%sizeDotState
       mySizeDamageDotState  = damageState( p)%sizeDotState
       mySizeThermalDotState = thermalState(p)%sizeDotState

       ! --- absolute residuum in state  ---
       ! NEED TO DO THE ADDITION IN THIS LENGTHY WAY BECAUSE OF PARALLELIZATION
       ! CAN'T USE A REDUCTION CLAUSE ON A POINTER OR USER DEFINED TYPE
       
       stateResiduum(1:mySizePlasticDotState,g,i,e) = &
                 (   DB(1) * plasticState(p)%RKCK45dotState(1,1:mySizePlasticDotState,cc) &
                   + DB(2) * plasticState(p)%RKCK45dotState(2,1:mySizePlasticDotState,cc) &
                   + DB(3) * plasticState(p)%RKCK45dotState(3,1:mySizePlasticDotState,cc) &
                   + DB(4) * plasticState(p)%RKCK45dotState(4,1:mySizePlasticDotState,cc) &
                   + DB(5) * plasticState(p)%RKCK45dotState(5,1:mySizePlasticDotState,cc) &
                   + DB(6) * plasticState(p)%RKCK45dotState(6,1:mySizePlasticDotState,cc) &
                 ) * crystallite_subdt(g,i,e)
       damageStateResiduum(1:mySizeDamageDotState,g,i,e) = &
                 (   DB(1) * damageState(p)%RKCK45dotState(1,1:mySizeDamageDotState,cc) &
                   + DB(2) * damageState(p)%RKCK45dotState(2,1:mySizeDamageDotState,cc) &
                   + DB(3) * damageState(p)%RKCK45dotState(3,1:mySizeDamageDotState,cc) &
                   + DB(4) * damageState(p)%RKCK45dotState(4,1:mySizeDamageDotState,cc) &
                   + DB(5) * damageState(p)%RKCK45dotState(5,1:mySizeDamageDotState,cc) &
                   + DB(6) * damageState(p)%RKCK45dotState(6,1:mySizeDamageDotState,cc) &
                 ) * crystallite_subdt(g,i,e)
       thermalStateResiduum(1:mySizethermalDotState,g,i,e) = &
                 (   DB(1) * thermalState(p)%RKCK45dotState(1,1:mySizeThermalDotState,cc) &
                   + DB(2) * thermalState(p)%RKCK45dotState(2,1:mySizeThermalDotState,cc) &
                   + DB(3) * thermalState(p)%RKCK45dotState(3,1:mySizeThermalDotState,cc) &
                   + DB(4) * thermalState(p)%RKCK45dotState(4,1:mySizeThermalDotState,cc) &
                   + DB(5) * thermalState(p)%RKCK45dotState(5,1:mySizeThermalDotState,cc) &
                   + DB(6) * thermalState(p)%RKCK45dotState(6,1:mySizeThermalDotState,cc) &
                 ) * crystallite_subdt(g,i,e)
 
       ! --- dot state ---

       plasticState(p)%dotState (:,cc) = B(1) * plasticState(p)%RKCK45dotState(1,:,cc) &
                                       + B(2) * plasticState(p)%RKCK45dotState(2,:,cc) &
                                       + B(3) * plasticState(p)%RKCK45dotState(3,:,cc) &
                                       + B(4) * plasticState(p)%RKCK45dotState(4,:,cc) &
                                       + B(5) * plasticState(p)%RKCK45dotState(5,:,cc) &
                                       + B(6) * plasticState(p)%RKCK45dotState(6,:,cc) 
       damageState(p)%dotState (:,cc)  = B(1) * damageState( p)%RKCK45dotState(1,:,cc) &
                                       + B(2) * damageState( p)%RKCK45dotState(2,:,cc) &
                                       + B(3) * damageState( p)%RKCK45dotState(3,:,cc) &
                                       + B(4) * damageState( p)%RKCK45dotState(4,:,cc) &
                                       + B(5) * damageState( p)%RKCK45dotState(5,:,cc) &
                                       + B(6) * damageState( p)%RKCK45dotState(6,:,cc) 
       thermalState(p)%dotState (:,cc) = B(1) * thermalState(p)%RKCK45dotState(1,:,cc) &
                                       + B(2) * thermalState(p)%RKCK45dotState(2,:,cc) &
                                       + B(3) * thermalState(p)%RKCK45dotState(3,:,cc) &
                                       + B(4) * thermalState(p)%RKCK45dotState(4,:,cc) &
                                       + B(5) * thermalState(p)%RKCK45dotState(5,:,cc) &
                                       + B(6) * thermalState(p)%RKCK45dotState(6,:,cc) 
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 
 ! --- state and update ---      
 
 !$OMP DO PRIVATE(mySizePlasticDotState,mySizeDamageDotState,mySizeThermalDotState,p,cc)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then

         p = mappingConstitutive(2,g,i,e)
         cc = mappingConstitutive(1,g,i,e)
         mySizePlasticDotState = plasticState(p)%sizeDotState
         mySizeDamageDotState  = damageState( p)%sizeDotState
         mySizeThermalDotState = thermalState(p)%sizeDotState
         plasticState(p)%state(1:mySizePlasticDotState,cc) = plasticState(p)%subState0(1:mySizePlasticDotState,cc)&
                                                           + plasticState(p)%dotState (1:mySizePlasticDotState,cc)&
                                                           * crystallite_subdt(g,i,e)
         damageState(p)%state (1:mySizeDamageDotState,cc)  = damageState(p)%subState0(1:mySizeDamageDotState,cc)&
                                                           + damageState(p)%dotState (1:mySizeDamageDotState,cc)&
                                                           * crystallite_subdt(g,i,e)
         thermalState(p)%state(1:mySizeThermalDotState,cc) = thermalState(p)%subState0(1:mySizeThermalDotState,cc)&
                                                           + thermalState(p)%dotState (1:mySizeThermalDotState,cc)&
                                                           * crystallite_subdt(g,i,e)
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 
 ! --- relative residui and state convergence ---      
 
 !$OMP DO PRIVATE(mySizePlasticDotState,mySizeDamageDotState,mySizeThermalDotState,p,cc,s)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       p = mappingConstitutive(2,g,i,e)
       cc = mappingConstitutive(1,g,i,e)
       mySizePlasticDotState = plasticState(p)%sizeDotState
       mySizeDamageDotState  = damageState( p)%sizeDotState
       mySizeThermalDotState = thermalState(p)%sizeDotState
       forall (s = 1_pInt:mySizePlasticDotState, abs(plasticState(p)%state(s,cc)) > 0.0_pReal) &
         relStateResiduum(s,g,i,e) = stateResiduum(s,g,i,e) / plasticState(p)%state(s,cc)
       forall (s = 1_pInt:mySizeDamageDotState,  abs(damageState( p)%state(s,cc)) > 0.0_pReal) &
         relDamageStateResiduum(s,g,i,e) = damageStateResiduum(s,g,i,e) / damageState(p)%state(s,cc)
       forall (s = 1_pInt:mySizeThermalDotState, abs(thermalState(p)%state(s,cc)) > 0.0_pReal) &
         relThermalStateResiduum(s,g,i,e) = thermalStateResiduum(s,g,i,e) / thermalState(p)%state(s,cc)
       !$OMP FLUSH(relStateResiduum)
       !$OMP FLUSH(relDamageStateResiduum)
       !$OMP FLUSH(relThermalStateResiduum)
! @Martin: do we need flushing? why..?
       crystallite_todo(g,i,e) = &
           (      all(abs(relStateResiduum(1:mySizePlasticDotState,g,i,e)) < &
                      rTol_crystalliteState .or. &
                      abs(stateResiduum(1:mySizePlasticDotState,g,i,e)) < &
                      plasticState(p)%aTolState(1:mySizePlasticDotState)) &
            .and. all(abs(relDamageStateResiduum(1:mySizeDamageDotState,g,i,e)) < &
                      rTol_crystalliteState .or. &
                      abs(damageStateResiduum(1:mySizeDamageDotState,g,i,e)) < &
                      damageState(p)%aTolState(1:mySizeDamageDotState)) &
            .and. all(abs(relThermalStateResiduum(1:mySizeThermalDotState,g,i,e)) < &
                      rTol_crystalliteState .or. &
                      abs(thermalStateResiduum(1:mySizeThermalDotState,g,i,e)) < &
                      thermalState(p)%aTolState(1:mySizeThermalDotState)))

#ifndef _OPENMP
       if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt&
           .and. ((e == debug_e .and. i == debug_i .and. g == debug_g)&
                  .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
         write(6,'(a,i8,1x,i3,1x,i3,/)') '<< CRYST >> updateState at el ip g ',e,i,g
         write(6,'(a,/,(12x,12(f12.1,1x)),/)') '<< CRYST >> absolute residuum tolerance', &
               stateResiduum(1:mySizePlasticDotState,g,i,e) / plasticState(p)%aTolState(1:mySizePlasticDotState)
         write(6,'(a,/,(12x,12(f12.1,1x)),/)') '<< CRYST >> relative residuum tolerance', &
               relStateResiduum(1:mySizePlasticDotState,g,i,e) / rTol_crystalliteState
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
       call constitutive_microstructure(crystallite_Tstar_v(1:6,g,i,e), &
                                        crystallite_Fe(1:3,1:3,g,i,e), &
                                        crystallite_Fp(1:3,1:3,g,i,e), g, i, e)                            ! update dependent state variables to be consistent with basic states
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
   damageState, &
   thermalState, &
   mappingConstitutive, &
   homogenization_maxNgrains
 use constitutive, only: &
   constitutive_collectDotState, &
   constitutive_microstructure, &
   constitutive_maxSizeDotState, &
   constitutive_damage_maxSizeDotState, &
   constitutive_thermal_maxSizeDotState

 implicit none
 integer(pInt) :: &
   e, &                                                                                             ! element index in element loop
   i, &                                                                                             ! integration point index in ip loop
   g, &                                                                                             ! grain index in grain loop
   s, &                                                                                             ! state index
   p, &
   c, &
   mySizePlasticDotState, &                                                                         ! size of dot states
   mySizeDamageDotState, &
   mySizeThermalDotState
 integer(pInt), dimension(2) :: &
   eIter                                                                                            ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) :: &
   iIter, &                                                                                         ! bounds for ip iteration
   gIter                                                                                            ! bounds for grain iteration
 real(pReal), dimension(constitutive_maxSizeDotState,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   stateResiduum, &                                                                                 ! residuum from evolution in micrstructure
   relStateResiduum                                                                                 ! relative residuum from evolution in microstructure
 real(pReal), dimension(constitutive_damage_maxSizeDotState,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   damageStateResiduum, &                                                                           ! residuum from evolution in micrstructure
   relDamageStateResiduum                                                                           ! relative residuum from evolution in microstructure
 real(pReal), dimension(constitutive_thermal_maxSizeDotState,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   thermalStateResiduum, &                                                                          ! residuum from evolution in micrstructure
   relThermalStateResiduum                                                                          ! relative residuum from evolution in microstructure

 logical :: &
   singleRun                                                                                        ! flag indicating computation for single (g,i,e) triple


 ! --- LOOP ITERATOR FOR ELEMENT, GRAIN, IP ---
 eIter = FEsolving_execElem(1:2)
 do e = eIter(1),eIter(2)
   iIter(1:2,e) = FEsolving_execIP(1:2,e)
   gIter(1:2,e) = [ 1_pInt,homogenization_Ngrains(mesh_element(3,e))]
 enddo

 singleRun = (eIter(1) == eIter(2) .and. iIter(1,eIter(1)) == iIter(2,eIter(2)))


 stateResiduum    = 0.0_pReal
 relStateResiduum = 0.0_pReal
 damageStateResiduum     = 0.0_pReal
 relDamageStateResiduum  = 0.0_pReal
 thermalStateResiduum    = 0.0_pReal
 relThermalStateResiduum = 0.0_pReal


 integrationMode: if (numerics_integrationMode == 1_pInt) then

 !$OMP PARALLEL
   ! --- DOT STATE (EULER INTEGRATION) ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) &
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Lp(1:3,1:3,g,i,e), &
                                           crystallite_Fe, &
                                           crystallite_Fp, &
                                           crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
    enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO PRIVATE(p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         p = mappingConstitutive(2,g,i,e) 
         c = mappingConstitutive(1,g,i,e)  
         if ( any(plasticState(p)%dotState(:,c) /= plasticState(p)%dotState(:,c)) .or. &
              any(damageState( p)%dotState(:,c) /= damageState( p)%dotState(:,c)) .or. &
              any(thermalState(p)%dotState(:,c) /= thermalState(p)%dotState(:,c))) then                      ! NaN occured in dotState
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
 
   !$OMP DO PRIVATE(mySizePlasticDotState,mySizeDamageDotState,mySizeThermalDotState,p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         p = mappingConstitutive(2,g,i,e)
         c = mappingConstitutive(1,g,i,e)
         mySizePlasticDotState = plasticState(p)%sizeDotState
         mySizeDamageDotState  = damageState( p)%sizeDotState
         mySizeThermalDotState = thermalState(p)%sizeDotState
         stateResiduum(1:mySizePlasticDotState,g,i,e) = - 0.5_pReal &
                                                        * plasticState(p)%dotstate(1:mySizePlasticDotState,c) &
                                                        * crystallite_subdt(g,i,e)            ! contribution to absolute residuum in state
         damageStateResiduum(1:mySizeDamageDotState,g,i,e) = - 0.5_pReal &
                                                             * damageState(p)%dotstate(1:mySizeDamageDotState,c) &
                                                             * crystallite_subdt(g,i,e)            ! contribution to absolute residuum in state
         thermalStateResiduum(1:mySizeThermalDotState,g,i,e) = - 0.5_pReal &
                                                               * thermalState(p)%dotstate(1:mySizeThermalDotState,c) &
                                                               * crystallite_subdt(g,i,e)            ! contribution to absolute residuum in state
         plasticState(p)%state(1:mySizePlasticDotState,c) = plasticState(p)%state(1:mySizePlasticDotState,c) &
                                                          + plasticState(p)%dotstate(1:mySizePlasticDotState,c) &
                                                          * crystallite_subdt(g,i,e)
         damageState(p)%state(1:mySizeDamageDotState,c)   = damageState(p)%state(1:mySizeDamageDotState,c) &
                                                          + damageState(p)%dotstate(1:mySizeDamageDotState,c) &
                                                          * crystallite_subdt(g,i,e)
         thermalState(p)%state(1:mySizeThermalDotState,c) = thermalState(p)%state(1:mySizeThermalDotState,c) &
                                                          + thermalState(p)%dotstate(1:mySizeThermalDotState,c) &
                                                          * crystallite_subdt(g,i,e)
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
         call constitutive_microstructure(crystallite_Tstar_v(1:6,g,i,e), & 
                                          crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), g, i, e)                          ! update dependent state variables to be consistent with basic states
     enddo; enddo; enddo
   !$OMP ENDDO
 !$OMP END PARALLEL
 endif integrationMode


 ! --- STRESS INTEGRATION (EULER INTEGRATION) ---

 !$OMP PARALLEL DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       crystallite_todo(g,i,e) = crystallite_integrateStress(g,i,e)
       !$OMP FLUSH(crystallite_todo)
       if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then            ! if broken non-local...
         !$OMP CRITICAL (checkTodo)
           crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
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
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Lp(1:3,1:3,g,i,e), &
                                           crystallite_Fe, &
                                           crystallite_Fp, &
                                           crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO PRIVATE(p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         p = mappingConstitutive(2,g,i,e)
         c = mappingConstitutive(1,g,i,e)
         if ( any(plasticState(p)%dotState(:,c)  /= plasticState(p)%dotState(:,c)) .or.&
              any(damageState( p)%dotState(:,c)  /= damageState( p)%dotState(:,c)) .or.&
              any(thermalState(p)%dotState(:,c)  /= thermalState(p)%dotState(:,c))) then             ! NaN occured in dotState
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
   relStateResiduum = 0.0_pReal
   relDamageStateResiduum  = 0.0_pReal
   relThermalStateResiduum = 0.0_pReal
   !$OMP END SINGLE

   !$OMP DO PRIVATE(mySizePlasticDotState,mySizeDamageDotState,mySizeThermalDotState,p,c,s)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                   ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         p = mappingConstitutive(2,g,i,e)
         c = mappingConstitutive(1,g,i,e)
         mySizePlasticDotState = plasticState(p)%sizeDotState
         mySizeDamageDotState  = damageState(p)%sizeDotState
         mySizeThermalDotState = thermalState(p)%sizeDotState
         ! --- contribution of heun step to absolute residui ---

         stateResiduum(1:mySizePlasticDotState,g,i,e) = stateResiduum(1:mySizePlasticDotState,g,i,e) &
                                                      + 0.5_pReal * plasticState(p)%dotState(:,c) &
                                                      * crystallite_subdt(g,i,e)                         ! contribution to absolute residuum in state
         damageStateResiduum(1:mySizeDamageDotState,g,i,e) = damageStateResiduum(1:mySizeDamageDotState,g,i,e) &
                                                           + 0.5_pReal * damageState(p)%dotState(:,c) &
                                                           * crystallite_subdt(g,i,e)                         ! contribution to absolute residuum in state
         thermalStateResiduum(1:mySizeThermalDotState,g,i,e) = thermalStateResiduum(1:mySizeThermalDotState,g,i,e) &
                                                             + 0.5_pReal * thermalState(p)%dotState(:,c) &
                                                             * crystallite_subdt(g,i,e)                      ! contribution to absolute residuum in state

         !$OMP FLUSH(stateResiduum)
         !$OMP FLUSH(damageStateResiduum)
         !$OMP FLUSH(thermalStateResiduum)
 
         ! --- relative residui ---  
         forall (s = 1_pInt:mySizePlasticDotState, abs(plasticState(p)%dotState(s,c)) > 0.0_pReal) &
           relStateResiduum(s,g,i,e) = stateResiduum(s,g,i,e) / plasticState(p)%dotState(s,c)
         forall (s = 1_pInt:mySizeDamageDotState, abs(damageState(p)%dotState(s,c)) > 0.0_pReal) &
           relDamageStateResiduum(s,g,i,e) = damageStateResiduum(s,g,i,e) / damageState(p)%dotState(s,c)
         forall (s = 1_pInt:mySizeThermalDotState, abs(thermalState(p)%dotState(s,c)) > 0.0_pReal) &
           relThermalStateResiduum(s,g,i,e) = thermalStateResiduum(s,g,i,e) / thermalState(p)%dotState(s,c)
         !$OMP FLUSH(relStateResiduum)
         !$OMP FLUSH(relDamageStateResiduum)
         !$OMP FLUSH(relthermalStateResiduum)
 
#ifndef _OPENMP 

         if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
             .and. ((e == debug_e .and. i == debug_i .and. g == debug_g)&
                     .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
           write(6,'(a,i8,1x,i2,1x,i3,/)')       '<< CRYST >> updateState at el ip g ',e,i,g
           write(6,'(a,/,(12x,12(f12.1,1x)),/)') '<< CRYST >> absolute residuum tolerance', &
                 stateResiduum(1:mySizePlasticDotState,g,i,e) / plasticState(p)%aTolState(1:mySizePlasticDotState)
           write(6,'(a,/,(12x,12(f12.1,1x)),/)') '<< CRYST >> relative residuum tolerance', &
                 relStateResiduum(1:mySizePlasticDotState,g,i,e) / rTol_crystalliteState
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> dotState', plasticState(p)%dotState(1:mySizePlasticDotState,c) &
                 - 2.0_pReal * stateResiduum(1:mySizePlasticDotState,g,i,e) / crystallite_subdt(g,i,e)     ! calculate former dotstate from higher order solution and state residuum
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> new state', plasticState(p)%state(1:mySizePlasticDotState,c)
         endif
#endif

         ! --- converged ? ---
         if (      all(abs(relStateResiduum(1:mySizePlasticDotState,g,i,e)) < &
                       rTol_crystalliteState .or. &
                       abs(stateResiduum(1:mySizePlasticDotState,g,i,e)) < &
                       plasticState(p)%aTolState(1:mySizePlasticDotState)) &
             .and. all(abs(relDamageStateResiduum(1:mySizeDamageDotState,g,i,e)) < &
                       rTol_crystalliteState .or. &
                       abs(damageStateResiduum(1:mySizeDamageDotState,g,i,e)) < &
                       damageState(p)%aTolState(1:mySizeDamageDotState)) &
             .and. all(abs(relThermalStateResiduum(1:mySizeThermalDotState,g,i,e)) < &
                       rTol_crystalliteState .or. &
                       abs(thermalStateResiduum(1:mySizeThermalDotState,g,i,e)) < &
                       thermalState(p)%aTolState(1:mySizeThermalDotState))) then
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
   damageState, &
   thermalState, &
   mappingConstitutive, &
   homogenization_Ngrains
 use constitutive, only: &
   constitutive_collectDotState, &
   constitutive_microstructure

 implicit none

 integer(pInt) :: &
   e, &                                                                                             ! element index in element loop
   i, &                                                                                             ! integration point index in ip loop
   g, &                                                                                             ! grain index in grain loop
   p, &                                                  ! phase loop
   c, &
   mySizePlasticDotState, &
   mySizeDamageDotState, &
   mySizeThermalDotState
 integer(pInt), dimension(2) :: &
   eIter                                                                                            ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) :: &
   iIter, &                                                                                         ! bounds for ip iteration
   gIter                                                                                            ! bounds for grain iteration
 logical :: &
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
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Lp(1:3,1:3,g,i,e), &
                                           crystallite_Fe, &
                                           crystallite_Fp, &
                                           crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO PRIVATE(p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         c = mappingConstitutive(1,g,i,e) 
         p = mappingConstitutive(2,g,i,e) 
         if ( any(plasticState(p)%dotState(:,c)  /= plasticState(p)%dotState(:,c)) .or. &
              any(damageState( p)%dotState(:,c)  /= damageState( p)%dotState(:,c)) .or. &
              any(thermalState(p)%dotState(:,c)  /= thermalState(p)%dotState(:,c))) then                     ! NaN occured in dotState
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
 
   !$OMP DO PRIVATE(mySizePlasticDotState,mySizeDamageDotState,mySizeThermalDotState,p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         p = mappingConstitutive(2,g,i,e)
         c = mappingConstitutive(1,g,i,e)
         mySizePlasticDotState = plasticState(p)%sizeDotState
         mySizeDamageDotState  = damageState( p)%sizeDotState
         mySizeThermalDotState = thermalState(p)%sizeDotState
         plasticState(p)%state(1:mySizePlasticDotState,c) = plasticState(p)%state(   1:mySizePlasticDotState,c) &
                                                          + plasticState(p)%dotState(1:mySizePlasticDotState,c) &
                                                          * crystallite_subdt(g,i,e)
         damageState( p)%state(1:mySizeDamageDotState,c)  = damageState( p)%state(   1:mySizeDamageDotState,c) &
                                                          + damageState( p)%dotState(1:mySizeDamageDotState,c) &
                                                          * crystallite_subdt(g,i,e)
         thermalState(p)%state(1:mySizeThermalDotState,c) = thermalState(p)%state(   1:mySizeThermalDotState,c) &
                                                          + thermalState(p)%dotState(1:mySizeThermalDotState,c) &
                                                          * crystallite_subdt(g,i,e)

#ifndef _OPENMP
         if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
             .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                     .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
           p = mappingConstitutive(2,g,i,e)
           c = mappingConstitutive(1,g,i,e)
           write(6,'(a,i8,1x,i2,1x,i3,/)')       '<< CRYST >> update state at el ip g ',e,i,g
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> dotState',  plasticState(p)%dotState (1:mySizePlasticDotState,c)
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> new state', plasticState(p)%state    (1:mySizePlasticDotState,c)
         endif
#endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- STATE JUMP ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         crystallite_todo(g,i,e) = crystallite_stateJump(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e) &                ! if broken non-local...
             .and. .not. numerics_timeSyncing) then
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO


   ! --- UPDATE DEPENDENT STATES ---

   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) &
         call constitutive_microstructure(crystallite_Tstar_v(1:6,g,i,e), &
                                          crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), g, i, e)                            ! update dependent state variables to be consistent with basic states
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
   damageState, &
   thermalState, &
   mappingConstitutive, &
   homogenization_Ngrains
 use constitutive, only: &
   constitutive_collectDotState, &
   constitutive_microstructure, &
   constitutive_maxSizeDotState, &
   constitutive_damage_maxSizeDotState, &
   constitutive_thermal_maxSizeDotState

 implicit none

 integer(pInt) :: &
   NiterationState, &                                                                               !< number of iterations in state loop
   e, &                                                                                             !< element index in element loop
   i, &                                                                                             !< integration point index in ip loop
   g, &                                                                                             !< grain index in grain loop
   p, &
   c, &
   mySizePlasticDotState, &                                                                         ! size of dot states
   mySizeDamageDotState, &
   mySizeThermalDotState
 integer(pInt), dimension(2) :: &
   eIter                                                                                            ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) :: &
   iIter, &                                                                                         ! bounds for ip iteration
   gIter                                                                                            ! bounds for grain iteration
 real(pReal) :: &
   dot_prod12, &
   dot_prod22, &
   stateDamper, &                                                                                   ! damper for integration of state
   damageStateDamper, &
   thermalStateDamper
 real(pReal), dimension(constitutive_maxSizeDotState) :: &
   stateResiduum, &
   tempState
 real(pReal), dimension(constitutive_damage_maxSizeDotState) :: &
   damageStateResiduum, &                                                                           ! residuum from evolution in micrstructure
   tempDamageState
 real(pReal), dimension(constitutive_thermal_maxSizeDotState) :: &
   thermalStateResiduum, &                                                                          ! residuum from evolution in micrstructure
   tempThermalState
 logical :: &
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
   forall(p = 1_pInt:size(damageState)) 
     damageState(p)%previousDotState  = 0.0_pReal
     damageState(p)%previousDotState2 = 0.0_pReal
   end forall
   forall(p = 1_pInt:size(thermalState)) 
     thermalState(p)%previousDotState  = 0.0_pReal
     thermalState(p)%previousDotState2 = 0.0_pReal
   end forall
 else
   e = eIter(1)
   i = iIter(1,e)
   do g = gIter(1,e), gIter(2,e)
     p = mappingConstitutive(2,g,i,e) 
     c = mappingConstitutive(1,g,i,e)  
     plasticState(p)%previousDotState (:,c) = 0.0_pReal   
     plasticState(p)%previousDotState2(:,c) = 0.0_pReal   
     damageState( p)%previousDotState (:,c) = 0.0_pReal   
     damageState( p)%previousDotState2(:,c) = 0.0_pReal   
     thermalState(p)%previousDotState (:,c) = 0.0_pReal   
     thermalState(p)%previousDotState2(:,c) = 0.0_pReal   
   enddo
 endif

 ! --+>> PREGUESS FOR STATE <<+--

 ! --- DOT STATES ---

 !$OMP PARALLEL
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)             ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) &
       call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Lp(1:3,1:3,g,i,e), &
                                         crystallite_Fe, &
                                         crystallite_Fp, &
                                         crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
   enddo; enddo; enddo

 !$OMP ENDDO
 !$OMP DO PRIVATE(p,c)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       p = mappingConstitutive(2,g,i,e) 
       c = mappingConstitutive(1,g,i,e)  
       if ( any(plasticState(p)%dotState(:,c) /= plasticState(p)%dotState(:,c)) .or. &
            any(damageState( p)%dotState(:,c) /= damageState( p)%dotState(:,c)) .or. &
            any(thermalState(p)%dotState(:,c) /= thermalState(p)%dotState(:,c))) then                      !NaN occured in dotState

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
 
 !$OMP DO PRIVATE(mySizePlasticDotState,mySizeDamageDotState,mySizeThermalDotState,p,c)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       p = mappingConstitutive(2,g,i,e)
       c = mappingConstitutive(1,g,i,e)
       mySizePlasticDotState = plasticState(p)%sizeDotState
       mySizeDamageDotState  = damageState( p)%sizeDotState
       mySizeThermalDotState = thermalState(p)%sizeDotState
       plasticState(p)%state(1:mySizePlasticDotState,c) = plasticState(p)%subState0(1:mySizePlasticDotState,c) &
                                                        + plasticState(p)%dotState (1:mySizePlasticDotState,c) &
                                                        * crystallite_subdt(g,i,e)
       damageState( p)%state(1:mySizeDamageDotState,c)  = damageState( p)%subState0(1:mySizeDamageDotState,c) &
                                                        + damageState( p)%dotState (1:mySizeDamageDotState,c) &
                                                        * crystallite_subdt(g,i,e)
       thermalState(p)%state(1:mySizeThermalDotState,c) = thermalState(p)%subState0(1:mySizeThermalDotState,c) &
                                                        + thermalState(p)%dotState (1:mySizeThermalDotState,c) &
                                                        * crystallite_subdt(g,i,e)
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
         call constitutive_microstructure(crystallite_Tstar_v(1:6,g,i,e), &
                                          crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), g, i, e)                            ! update dependent state variables to be consistent with basic states
       p = mappingConstitutive(2,g,i,e) 
       c = mappingConstitutive(1,g,i,e)  
       plasticState(p)%previousDotState2(:,c) = plasticState(p)%previousDotState(:,c) 
       damageState( p)%previousDotState2(:,c) = damageState( p)%previousDotState(:,c) 
       thermalState(p)%previousDotState2(:,c) = thermalState(p)%previousDotState(:,c) 
       plasticState(p)%previousDotState (:,c) = plasticState(p)%dotState(:,c)          
       damageState( p)%previousDotState (:,c) = damageState( p)%dotState(:,c)          
       thermalState(p)%previousDotState (:,c) = thermalState(p)%dotState(:,c)          
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
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Lp(1:3,1:3,g,i,e), &
                                           crystallite_Fe, &
                                           crystallite_Fp, &
                                           crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
     enddo; enddo; enddo
   !$OMP ENDDO

   !$OMP DO PRIVATE(p,c)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         p = mappingConstitutive(2,g,i,e) 
         c = mappingConstitutive(1,g,i,e)  
         if ( any(plasticState(p)%dotState(:,c) /= plasticState(p)%dotState(:,c)) .or. &
              any(damageState( p)%dotState(:,c) /= damageState( p)%dotState(:,c)) .or. &
              any(thermalState(p)%dotState(:,c) /= thermalState(p)%dotState(:,c))) then                      ! NaN occured in dotState
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
   !$OMP&           mySizePlasticDotState,mySizeDamageDotState,mySizeThermalDotState, &
   !$OMP&           damageStateResiduum,thermalStateResiduum,damageStateDamper,thermalStateDamper, &
   !$OMP&           tempDamageState,tempThermalState,p,c, &
   !$OMP&           statedamper,stateResiduum,tempState)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)           ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then

         p = mappingConstitutive(2,g,i,e)
         c = mappingConstitutive(1,g,i,e)
         mySizePlasticDotState = plasticState(p)%sizeDotState
         mySizeDamageDotState  = damageState( p)%sizeDotState
         mySizeThermalDotState = thermalState(p)%sizeDotState

         dot_prod12 = dot_product(plasticState(p)%dotState(:,c)         - plasticState(p)%previousDotState(:,c), &
                                  plasticState(p)%previousDotState(:,c) - plasticState(p)%previousDotState2(:,c))
         dot_prod22 = dot_product(plasticState(p)%previousDotState(:,c) - plasticState(p)%previousDotState2(:,c), &
                                  plasticState(p)%previousDotState(:,c) - plasticState(p)%previousDotState2(:,c))
         if (      dot_prod22 > 0.0_pReal &
             .and. (     dot_prod12 < 0.0_pReal &
                    .or. dot_product(plasticState(p)%dotState(:,c), &
                                     plasticState(p)%previousDotState(:,c)) < 0.0_pReal) ) then
           stateDamper = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)
         else
           stateDamper = 1.0_pReal
         endif

         dot_prod12 = dot_product(damageState(p)%dotState(:,c)         - damageState(p)%previousDotState(:,c), &
                                  damageState(p)%previousDotState(:,c) - damageState(p)%previousDotState2(:,c))
         dot_prod22 = dot_product(damageState(p)%previousDotState(:,c) - damageState(p)%previousDotState2(:,c), &
                                  damageState(p)%previousDotState(:,c) - damageState(p)%previousDotState2(:,c))
         if (      dot_prod22 > 0.0_pReal &
             .and. (     dot_prod12 < 0.0_pReal &
                    .or. dot_product(damageState(p)%dotState(:,c), &
                                     damageState(p)%previousDotState(:,c)) < 0.0_pReal) ) then
           damageStateDamper = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)
         else
           damageStateDamper = 1.0_pReal
         endif        

         dot_prod12 = dot_product(thermalState(p)%dotState(:,c)         - thermalState(p)%previousDotState(:,c), &
                                  thermalState(p)%previousDotState(:,c) - thermalState(p)%previousDotState2(:,c))
         dot_prod22 = dot_product(thermalState(p)%previousDotState(:,c) - thermalState(p)%previousDotState2(:,c), &
                                  thermalState(p)%previousDotState(:,c) - thermalState(p)%previousDotState2(:,c))
         if (      dot_prod22 > 0.0_pReal &
             .and. (     dot_prod12 < 0.0_pReal &
                    .or. dot_product(thermalState(p)%dotState(:,c), &
                                     thermalState(p)%previousDotState(:,c)) < 0.0_pReal) ) then
           thermalStateDamper = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)
         else
           thermalStateDamper = 1.0_pReal
         endif
         ! --- get residui ---

         stateResiduum(1:mySizePlasticDotState) = plasticState(p)%state(1:mySizePlasticDotState,c)      &
                                                - plasticState(p)%subState0(1:mySizePlasticDotState,c)  &
                                                - (plasticState(p)%dotState(1:mySizePlasticDotState,c) * stateDamper &
                                                 + plasticState(p)%previousDotState(1:mySizePlasticDotState,c) &
                                                * (1.0_pReal - stateDamper)) * crystallite_subdt(g,i,e)

         damageStateResiduum(1:mySizeDamageDotState) = damageState(p)%state(1:mySizeDamageDotState,c)      &
                                                     - damageState(p)%subState0(1:mySizeDamageDotState,c)  &
                                                     - (damageState(p)%dotState(1:mySizeDamageDotState,c) * damageStateDamper &
                                                      + damageState(p)%previousDotState(1:mySizeDamageDotState,c) &
                                                     * (1.0_pReal - damageStatedamper)) * crystallite_subdt(g,i,e)

         thermalStateResiduum(1:mySizeThermalDotState) = thermalState(p)%state(1:mySizeThermalDotState,c)      &
                                                       - thermalState(p)%subState0(1:mySizeThermalDotState,c)  &
                                                       - (thermalState(p)%dotState(1:mySizeThermalDotState,c) * thermalStateDamper &
                                                        + thermalState(p)%previousDotState(1:mySizeThermalDotState,c) &
                                                       * (1.0_pReal - thermalStateDamper)) * crystallite_subdt(g,i,e)
         ! --- correct state with residuum ---
         tempState(1:mySizePlasticDotState) = plasticState(p)%state(1:mySizePlasticDotState,c) &
                                            - stateResiduum(1:mySizePlasticDotState)                              ! need to copy to local variable, since we cant flush a pointer in openmp
         tempDamageState(1:mySizeDamageDotState) = damageState(p)%state(1:mySizeDamageDotState,c) &
                                                 - damageStateResiduum(1:mySizeDamageDotState)                              ! need to copy to local variable, since we cant flush a pointer in openmp
         tempThermalState(1:mySizeThermalDotState) = thermalState(p)%state(1:mySizeThermalDotState,c) &
                                                   - thermalStateResiduum(1:mySizeThermalDotState)                              ! need to copy to local variable, since we cant flush a pointer in openmp
#ifndef _OPENMP
         if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
             .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                    .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
           write(6,'(a,i8,1x,i2,1x,i3,/)')       '<< CRYST >> update state at el ip g ',e,i,g
           write(6,'(a,f6.1,/)')                 '<< CRYST >> statedamper ',statedamper
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> state residuum',stateResiduum(1:mySizePlasticDotState)
           write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> new state',tempState(1:mySizePlasticDotState)
         endif
#endif
         ! --- store corrected dotState --- (cannot do this before state update, because not sure how to flush pointers in openmp)

         plasticState(p)%dotState(:,c) = plasticState(p)%dotState(:,c) * stateDamper &
                                       + plasticState(p)%previousDotState(:,c) &
                                       * (1.0_pReal - stateDamper)
         damageState( p)%dotState(:,c) = damageState(p)%dotState(:,c) * damageStateDamper &
                                       + damageState(p)%previousDotState(:,c) &
                                       * (1.0_pReal - damageStateDamper)
         thermalState(p)%dotState(:,c) = thermalState(p)%dotState(:,c) * thermalStateDamper &
                                       + thermalState(p)%previousDotState(:,c) &
                                       * (1.0_pReal - thermalStateDamper)
         ! --- converged ? ---
         if ( all(    abs(stateResiduum(1:mySizePlasticDotState))        < plasticState(p)%aTolState(1:mySizePlasticDotState) &
                 .or. abs(stateResiduum(1:mySizePlasticDotState))        < rTol_crystalliteState &
                                                                         * abs(tempState(1:mySizePlasticDotState)) ) &
        .and. all(    abs(damageStateResiduum(1:mySizeDamageDotState))   < damageState(p)%aTolState(1:mySizeDamageDotState) &
                 .or. abs(damageStateResiduum(1:mySizeDamageDotState))   < rTol_crystalliteState &
                                                                         * abs(tempDamageState(1:mySizeDamageDotState)) ) &
        .and. all(    abs(thermalStateResiduum(1:mySizeThermalDotState)) < thermalState(p)%aTolState(1:mySizeThermalDotState) &
                 .or. abs(thermalStateResiduum(1:mySizeThermalDotState)) < rTol_crystalliteState &
                                                                         * abs(tempThermalState(1:mySizeThermalDotState)) )) then
           crystallite_converged(g,i,e) = .true.                                                                   ! ... converged per definition

           if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
             !$OMP CRITICAL (distributionState)
               debug_StateLoopDistribution(NiterationState,numerics_integrationMode) = &
                 debug_StateLoopDistribution(NiterationState,numerics_integrationMode) + 1_pInt
             !$OMP END CRITICAL (distributionState)
           endif
         endif
         plasticState(p)%state(1:mySizePlasticDotState,c) = tempState(1:mySizePlasticDotState)
         damageState( p)%state(1:mySizeDamageDotState, c) = tempDamageState(1:mySizeDamageDotState)
         thermalState(p)%state(1:mySizeThermalDotState,c) = tempThermalState(1:mySizeThermalDotState)
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
logical function crystallite_stateJump(g,i,e)
 use debug, only: &
   debug_level, &
   debug_crystallite, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_g
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP
 use mesh, only: &
   mesh_element, &
   mesh_NcpElems
 use material, only: &
   plasticState, &
   damageState, &
   thermalState, &
   mappingConstitutive, &
   homogenization_Ngrains
 use constitutive, only: &
   constitutive_collectDeltaState

 implicit none
 integer(pInt), intent(in):: &
   e, &                      ! element index
   i, &                      ! integration point index
   g                         ! grain index

 integer(pInt) :: &
   c, &
   p, &
   mySizePlasticDotState

 c = mappingConstitutive(1,g,i,e)
 p = mappingConstitutive(2,g,i,e)
 if (constitutive_collectDeltaState(crystallite_Tstar_v(1:6,g,i,e), g,i,e)) then
   mySizePlasticDotState = plasticState(p)%sizeDotState
   if( any(plasticState(p)%deltaState(:,c) /= plasticState(p)%deltaState(:,c))) then    ! NaN occured in deltaState
     crystallite_stateJump = .false.
     return
   endif
   plasticState(p)%state(1:mySizePlasticDotState,c) = plasticState(p)%state(1:mySizePlasticDotState,c) + &
                                                      plasticState(p)%deltaState(1:mySizePlasticDotState,c)

#ifndef _OPENMP
   p = mappingConstitutive(2,g,i,e)
   c = mappingConstitutive(1,g,i,e)
   if (any(plasticState(p)%deltaState(1:mySizePlasticDotState,c) /= 0.0_pReal) &
       .and. iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
       .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
               .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
     write(6,'(a,i8,1x,i2,1x,i3, /)') '<< CRYST >> update state at el ip g ',e,i,g
     write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> deltaState', plasticState(p)%deltaState(1:mySizePlasticDotState,c)
     write(6,'(a,/,(12x,12(e12.5,1x)),/)') '<< CRYST >> new state',  plasticState(p)%state     (1:mySizePlasticDotState,c)
   endif
#endif
 endif

 crystallite_stateJump = .true.

end function crystallite_stateJump
!--------------------------------------------------------------------------------------------------
!> @brief Map 2nd order tensor to reference config
!--------------------------------------------------------------------------------------------------
function crystallite_push33ToRef(g,i,e, tensor33)
 use math, only: &
  math_inv33

 implicit none
 real(pReal), dimension(3,3) :: crystallite_push33ToRef
 real(pReal), dimension(3,3), intent(in) :: tensor33
 real(pReal), dimension(3,3)              :: invFp
 integer(pInt), intent(in):: &
   e, &                      ! element index
   i, &                      ! integration point index
   g                         ! grain index

 invFp = math_inv33(crystallite_Fp(1:3,1:3,g,i,e))
 crystallite_push33ToRef = matmul(invFp,matmul(tensor33,transpose(invFp)))

end function crystallite_push33ToRef
!--------------------------------------------------------------------------------------------------
!> @brief calculation of stress (P) with time integration based on a residuum in Lp and
!> intermediate acceleration of the Newton-Raphson correction
!--------------------------------------------------------------------------------------------------
logical function crystallite_integrateStress(&
      g,&          ! grain number
      i,&          ! integration point number
      e,&          ! element number
      timeFraction &
      )
 use prec, only:         pLongInt
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
                         debug_StressLoopDistribution
 use constitutive, only: constitutive_LpAndItsTangent, &
                         constitutive_LiAndItsTangent, &
                         constitutive_getFi0, &
                         constitutive_TandItsTangent
 use math, only:         math_mul33x33, &
                         math_mul33xx33, &
                         math_mul3333xx3333, &
                         math_mul66x6, &
                         math_mul99x99, &
                         math_transpose33, &
                         math_inv33, &
                         math_invert33, &
                         math_invert, &
                         math_det33, &
                         math_norm33, &
                         math_I3, &
                         math_identity2nd, &
                         math_Mandel66to3333, &
                         math_Mandel6to33, &
                         math_Mandel33to6, &
                         math_Plain3333to99, &
                         math_Plain33to9, &
                         math_Plain9to33
 use mesh, only:         mesh_element

 implicit none
 integer(pInt), intent(in)::         e, &                          ! element index
                                     i, &                          ! integration point index
                                     g                             ! grain index
 real(pReal), optional, intent(in) :: timeFraction                 ! fraction of timestep

 !*** local variables ***!
 real(pReal), dimension(3,3)::       Fg_new, &                                                       ! deformation gradient at end of timestep
                                     Fp_current, &                                                   ! plastic deformation gradient at start of timestep
                                     Fi_current, &                                                   ! intermediate deformation gradient at start of timestep
                                     Fp_new, &                                                       ! plastic deformation gradient at end of timestep
                                     Fe_new, &                                                       ! elastic deformation gradient at end of timestep
                                     invFp_new, &                                                    ! inverse of Fp_new
                                     invFp_current, &                                                ! inverse of Fp_current
                                     invFi_current, &                                                ! inverse of Fp_current
                                     Lpguess, &                                                      ! current guess for plastic velocity gradient
                                     Lpguess_old, &                                                  ! known last good guess for plastic velocity gradient
                                     Lp_constitutive, &                                              ! plastic velocity gradient resulting from constitutive law
                                     residuum, &                                                     ! current residuum of plastic velocity gradient
                                     residuum_old, &                                                 ! last residuum of plastic velocity gradient
                                     deltaLp, &                                                      ! direction of next guess
                                     Liguess, &                                                      ! current guess for plastic velocity gradient
                                     Liguess_old, &                                                  ! known last good guess for plastic velocity gradient
                                     Li_constitutive, &                                              ! plastic velocity gradient resulting from constitutive law
                                     residuumI, &                                                    ! current residuum of plastic velocity gradient
                                     residuumI_old, &                                                ! last residuum of plastic velocity gradient
                                     deltaLi, &                                                      ! direction of next guess
                                     Tstar, &                                                        ! 2nd Piola-Kirchhoff Stress in plastic (lattice) configuration
                                     Tstar_unloaded, &                                               ! 2nd Piola-Kirchhoff Stress in unloaded configuration
                                     A, &
                                     B, &
                                     Fe, &                                                           ! elastic deformation gradient
                                     Fi, &                                                           ! gradient of intermediate deformation stages
                                     invFi, &
                                     temp_33
 real(pReal), dimension(6)::         Tstar_v                                                         ! 2nd Piola-Kirchhoff Stress in Mandel-Notation
 real(pReal), dimension(9)::         work                                                            ! needed for matrix inversion by LAPACK
 integer(pInt), dimension(9) ::      ipiv                                                            ! needed for matrix inversion by LAPACK
 real(pReal), dimension(9,9) ::      dLp_dT_constitutive, &                                          ! partial derivative of plastic velocity gradient calculated by constitutive law
                                     dLi_dT_constitutive, &                                          ! partial derivative of intermediate velocity gradient calculated by constitutive law
                                     dT_dFe_constitutive, &                                          ! partial derivative of 2nd Piola-Kirchhoff stress calculated by constitutive law
                                     dFe_dLp, &                                                      ! partial derivative of elastic deformation gradient
                                     dR_dLp, &                                                       ! partial derivative of residuum (Jacobian for NEwton-Raphson scheme)
                                     dR_dLp2, &                                                      ! working copy of dRdLp
                                     dRI_dLp, &                                                      ! partial derivative of residuumI (Jacobian for NEwton-Raphson scheme)
                                     dRI_dLp2                                                        ! working copy of dRIdLp
 real(pReal), dimension(3,3,3,3)::   dT_dFe3333, &                                                   ! partial derivative of 2nd Piola-Kirchhoff stress
                                     dFe_dLp3333, &                                                  ! partial derivative of elastic deformation gradient
                                     dInvFi_dLp3333, &
                                     dFe_dInvFi3333, &
                                     dT_dInvFi3333, &
                                     temp_3333
 real(pReal)                         det, &                                                          ! determinant
                                     detInvFi, &
                                     steplength0, &
                                     steplength, &
                                     steplengthI0, &
                                     steplengthI, &
                                     dt, &                                                           ! time increment
                                     aTol
 logical                             error                                                           ! flag indicating an error
 integer(pInt)                       NiterationStress, &                                             ! number of stress integrations
                                     NiterationStressI, &                                            ! number of inner stress integrations
                                     ierr, &                                                         ! error indicator for LAPACK
                                     o, &
                                     p, &
                                     jacoCounter, &
                                     jacoICounter                                                    ! counters to check for Jacobian update
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
     .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
            .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
   write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> integrateStress at el ip g ',e,i,g
 endif
#endif


 !* only integrate over fraction of timestep?

 if (present(timeFraction)) then
   dt = crystallite_subdt(g,i,e) * timeFraction
   Fg_new = crystallite_subF0(1:3,1:3,g,i,e) &
          + (crystallite_subF(1:3,1:3,g,i,e) - crystallite_subF0(1:3,1:3,g,i,e)) * timeFraction
 else
   dt = crystallite_subdt(g,i,e)
   Fg_new = crystallite_subF(1:3,1:3,g,i,e)
 endif


 !* feed local variables

 Fp_current =   crystallite_subFp0(1:3,1:3,g,i,e)                   ! "Fp_current" is only used as temp var here...
 Lpguess_old =  crystallite_Lp(1:3,1:3,g,i,e)                       ! consider present Lp good (i.e. worth remembering) ...
 Lpguess =      crystallite_Lp(1:3,1:3,g,i,e)                       ! ... and take it as first guess


 !* inversion of Fp_current...

 invFp_current = math_inv33(Fp_current)
 if (all(invFp_current == 0.0_pReal)) then                          ! ... failed?
#ifndef _OPENMP
   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3)') '<< CRYST >> integrateStress failed on inversion of Fp_current at el (elFE) ip g ',&
       e,'(',mesh_element(1,e),')',i,g
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) > 0_pInt) &
       write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fp_current',math_transpose33(Fp_current(1:3,1:3))
   endif
#endif
   return
 endif
 A = math_mul33x33(Fg_new,invFp_current)                                    ! intermediate tensor needed later to calculate dFe_dLp

 Fi_current = constitutive_getFi0(g,i,e)                                    ! intermediate configuration, assume decomposition as F = Fe Fi Fp
 invFi_current = math_inv33(Fi_current)
 if (all(invFi_current == 0.0_pReal)) then                          ! ... failed?
#ifndef _OPENMP
   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3)') '<< CRYST >> integrateStress failed on inversion of Fi_current at el (elFE) ip g ',&
       e,'(',mesh_element(1,e),')',i,g
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) > 0_pInt) &
       write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fp_current',math_transpose33(Fi_current(1:3,1:3))
   endif
#endif
   return
 endif

 !* start LpLoop with normal step length

 NiterationStress = 0_pInt
 jacoCounter      = 0_pInt
 steplength0      = 1.0_pReal
 steplength       = steplength0
 residuum_old     = 0.0_pReal

 LpLoop: do
   NiterationStress = NiterationStress + 1_pInt
   loopsExeced: if (NiterationStress > nStress) then
#ifndef _OPENMP
     if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) &
       write(6,'(a,i3,a,i8,1x,a,i8,a,1x,i2,1x,i3,/)') '<< CRYST >> integrateStress reached loop limit',nStress, &
       ' at el (elFE) ip g ', e,mesh_element(1,e),i,g
#endif
     return
   endif loopsExeced
   B                 = math_I3 - dt*Lpguess
   NiterationStressI = 0_pInt
   Liguess           = 0.0_pReal
   Liguess_old       = Liguess
   jacoICounter      = 0_pInt
   steplengthI0      = 1.0_pReal
   steplengthI       = steplengthI0
   residuumI_old     = 0.0_pReal
   LiLoop: do                                    ! inner stress integration loop for consistency with Fi 
     NiterationStressI = NiterationStressI + 1_pInt
     IloopsExeced: if (NiterationStressI > nStress) then
#ifndef _OPENMP
     if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) &
         write(6,'(a,i3,a,i8,1x,a,i8,a,1x,i2,1x,i3,/)') '<< CRYST >> integrateStress reached inelastic loop limit',nStress, &
         ' at el (elFE) ip g ', e,mesh_element(1,e),i,g
#endif
       return
     endif IloopsExeced
     
     !* calculate (elastic) 2nd Piola--Kirchhoff stress tensor and its tangent from constitutive law

     invFi = math_mul33x33(invFi_current,math_I3 - dt*Liguess)
     detInvFi = math_det33(invFi)
     Fi = math_inv33(invFi)
     Fe = math_mul33x33(math_mul33x33(A,B),invFi)                                                       ! current elastic deformation tensor
     call constitutive_TandItsTangent(Tstar_unloaded, dT_dFe3333, Fe, g,i,e)                            ! call constitutive law to calculate 2nd Piola-Kirchhoff stress and its derivative in unloaded configuration
     Tstar = math_mul33x33(invFi, &
               math_mul33x33(Tstar_unloaded,math_transpose33(invFi)))/detInvFi                          ! push Tstar forward from unloaded to plastic (lattice) configuration 
     Tstar_v = math_Mandel33to6(Tstar)

     !* calculate intermediate velocity gradient and its tangent from constitutive law

     call constitutive_LiAndItsTangent(Li_constitutive, dLi_dT_constitutive, Tstar_v, Lpguess, &
                                       g, i, e)

     !* update current residuum and check for convergence of loop

     aTol = max(rTol_crystalliteStress * max(math_norm33(Liguess),math_norm33(Li_constitutive)), &      ! absolute tolerance from largest acceptable relative error
                aTol_crystalliteStress)                                                                 ! minimum lower cutoff
     residuumI = Liguess - Li_constitutive
     if (any(residuumI /= residuumI)) then                                                              ! NaN in residuum...
       return                                                                                           ! ...me = .false. to inform integrator about problem
     elseif (math_norm33(residuumI) < aTol) then                                                        ! converged if below absolute tolerance
       exit LiLoop                                                                                      ! ...leave iteration loop
     elseif (math_norm33(residuumI) < math_norm33(residuumI_old) .or. NiterationStress == 1_pInt ) then ! not converged, but improved norm of residuum (always proceed in first iteration)...
       residuumI_old = residuumI                                                                        ! ...remember old values and...
       Liguess_old = Liguess
       steplengthI = steplengthI0                                                                       ! ...proceed with normal step length (calculate new search direction)
     else                                                                                               ! not converged and residuum not improved...
       steplengthI = 0.5_pReal * steplengthI                                                            ! ...try with smaller step length in same direction
       Liguess = Liguess_old + steplengthI * deltaLi
       cycle LiLoop
     endif

     !* calculate Jacobian for correction term

     if (mod(jacoICounter, iJacoLpresiduum) == 0_pInt) then
       dInvFi_dLp3333 = 0.0_pReal
       do o=1_pInt,3_pInt
         dInvFi_dLp3333(1:3,o,1:3,o) = invFi_current
       enddo 
       dInvFi_dLp3333 = -dt * dInvFi_dLp3333 
       dFe_dInvFi3333 = 0.0_pReal
       temp_33 = math_mul33x33(A,B)
       do o=1_pInt,3_pInt
         dFe_dInvFi3333(1:3,o,1:3,o) = temp_33
       enddo
       dT_dInvFi3333 = 0.0_pReal
       do o=1_pInt,3_pInt; do p=1_pInt,3_pInt
         dT_dInvFi3333(1:3,1:3,p,o) = -Tstar * Fi(o,p)
       enddo; enddo
       temp_33 = math_mul33x33(invFi,Tstar_unloaded)/detInvFi         
       do o=1_pInt,3_pInt
         dT_dInvFi3333(1:3,o,o,1:3) = dT_dInvFi3333(1:3,o,o,1:3) + temp_33
         dT_dInvFi3333(o,1:3,o,1:3) = dT_dInvFi3333(o,1:3,o,1:3) + temp_33
       enddo
       temp_3333 = math_mul3333xx3333(dT_dFe3333,dFe_dInvFi3333)
       do o=1_pInt,3_pInt; do p=1_pInt,3_pInt
         dT_dInvFi3333(1:3,1:3,o,p) = dT_dInvFi3333(1:3,1:3,o,p) + &
                                      math_mul33x33(math_mul33x33(invFi,temp_3333(1:3,1:3,o,p)), &
                                                    math_transpose33(invFi))/detInvFi
       enddo; enddo
       dRI_dLp = math_identity2nd(9_pInt) - &
                math_mul99x99(dLi_dT_constitutive,math_Plain3333to99(math_mul3333xx3333(dT_dInvFi3333,dInvFi_dLp3333)))
       dRI_dLp2 = dRI_dLp                                                                               ! will be overwritten in first call to LAPACK routine
       work = math_plain33to9(residuumI)
#if(FLOAT==8)
       call dgesv(9,1,dRI_dLp2,9,ipiv,work,9,ierr)                                                      ! solve dR/dLp * delta Lp = -res for dR/dLp
#elif(FLOAT==4)
       call sgesv(9,1,dRI_dLp2,9,ipiv,work,9,ierr)                                                      ! solve dR/dLp * delta Lp = -res for dR/dLp
#endif
       if (ierr /= 0_pInt) then
         return
       endif
       deltaLi = - math_plain9to33(work)
     endif
     jacoICounter = jacoICounter + 1_pInt                                                               ! increase counter for jaco update

     Liguess = Liguess + steplengthI * deltaLi
   enddo LiLoop  

   !* calculate plastic velocity gradient and its tangent from constitutive law

   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
   endif

   call constitutive_LpAndItsTangent(Lp_constitutive, dLp_dT_constitutive, Tstar_v, &
                                     g, i, e)

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
       .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
              .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
     write(6,'(a,i3,/)') '<< CRYST >> iteration ', NiterationStress
     write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Lp_constitutive', math_transpose33(Lp_constitutive)
     write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Lpguess', math_transpose33(Lpguess)
   endif
#endif


   !* update current residuum and check for convergence of loop

   aTol = max(rTol_crystalliteStress * max(math_norm33(Lpguess),math_norm33(Lp_constitutive)), &    ! absolute tolerance from largest acceptable relative error
              aTol_crystalliteStress)                                                               ! minimum lower cutoff
   residuum = Lpguess - Lp_constitutive

   if (any(residuum /= residuum)) then                                                              ! NaN in residuum...
#ifndef _OPENMP
     if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) &
       write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3,a,i3,a)') '<< CRYST >> integrateStress encountered NaN at el (elFE) ip g ', &
         e,mesh_element(1,e),i,g, &
                                                      ' ; iteration ', NiterationStress,&
                                                      ' >> returning..!'
#endif
     return                                                                                         ! ...me = .false. to inform integrator about problem
   elseif (math_norm33(residuum) < aTol) then                                                       ! converged if below absolute tolerance
     exit LpLoop                                                                                    ! ...leave iteration loop
   elseif (math_norm33(residuum) < math_norm33(residuum_old) .or. NiterationStress == 1_pInt ) then ! not converged, but improved norm of residuum (always proceed in first iteration)...
     residuum_old = residuum                                                                        ! ...remember old values and...
     Lpguess_old = Lpguess
     steplength = steplength0                                                                       ! ...proceed with normal step length (calculate new search direction)
   else                                                                                             ! not converged and residuum not improved...
     steplength = 0.5_pReal * steplength                                                            ! ...try with smaller step length in same direction
     Lpguess = Lpguess_old + steplength * deltaLp
     cycle LpLoop
   endif


   !* calculate Jacobian for correction term

   if (mod(jacoCounter, iJacoLpresiduum) == 0_pInt) then
     dFe_dLp3333 = 0.0_pReal
     do o=1_pInt,3_pInt; do p=1_pInt,3_pInt
       dFe_dLp3333(o,1:3,p,1:3) = A(o,p)*math_transpose33(invFi)                                    ! dFe_dLp(i,j,k,l) = -dt * A(i,k) invFi(l,j)
     enddo; enddo
     dFe_dLp3333 = -dt * dFe_dLp3333
     dFe_dLp = math_Plain3333to99(dFe_dLp3333)
     dT_dFe_constitutive = math_Plain3333to99(dT_dFe3333)
     dR_dLp = math_identity2nd(9_pInt) - &
              math_mul99x99(dLp_dT_constitutive, math_mul99x99(dT_dFe_constitutive , dFe_dLp))
     dR_dLp2 = dR_dLp                                                                               ! will be overwritten in first call to LAPACK routine
     work = math_plain33to9(residuum)
#if(FLOAT==8)
     call dgesv(9,1,dR_dLp2,9,ipiv,work,9,ierr)                                                     ! solve dR/dLp * delta Lp = -res for dR/dLp
#elif(FLOAT==4)
     call sgesv(9,1,dR_dLp2,9,ipiv,work,9,ierr)                                                     ! solve dR/dLp * delta Lp = -res for dR/dLp
#endif
     if (ierr /= 0_pInt) then
#ifndef _OPENMP
       if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
         write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3,a,i3)') '<< CRYST >> integrateStress failed on dR/dLp inversion at el ip g ', &
           e,mesh_element(1,e),i,g
         if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
             .and. ((e == debug_e .and. i == debug_i .and. g == debug_g)&
                    .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
           write(6,*)
           write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dR_dLp',transpose(dR_dLp)
           write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dFe_dLp',transpose(dFe_dLp)
           write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dT_dFe_constitutive',transpose(dT_dFe_constitutive)
           write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dLp_dT_constitutive',transpose(dLp_dT_constitutive)
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
   jacoCounter = jacoCounter + 1_pInt                             ! increase counter for jaco update

   Lpguess = Lpguess + steplength * deltaLp

 enddo LpLoop


 !* calculate new plastic and elastic deformation gradient

 invFp_new = math_mul33x33(invFp_current,B)
 invFp_new = invFp_new / math_det33(invFp_new)**(1.0_pReal/3.0_pReal)  ! regularize by det
 call math_invert33(invFp_new,Fp_new,det,error)
 if (error .or. any(Fp_new /= Fp_new)) then
#ifndef _OPENMP
   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     write(6,'(a,i8,1x,a,i8,a,1x,i2,1x,i3,a,i3)') '<< CRYST >> integrateStress failed on invFp_new inversion at el ip g ',&
       e,mesh_element(1,e),i,g, ' ; iteration ', NiterationStress
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
         .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) &
       write(6,'(/,a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> invFp_new',math_transpose33(invFp_new)
   endif
#endif
   return
 endif
 Fe_new = math_mul33x33(math_mul33x33(Fg_new,invFp_new),invFi)    ! calc resulting Fe


 !* calculate 1st Piola-Kirchhoff stress

 crystallite_P(1:3,1:3,g,i,e) = math_mul33x33(math_mul33x33(Fe_new,Fi), &
                                              math_mul33x33(math_Mandel6to33(Tstar_v), &
                                                            math_transpose33(invFp_new)))/detInvFi


 !* store local values in global variables

 crystallite_Lp(1:3,1:3,g,i,e) = Lpguess
 crystallite_Tstar_v(1:6,g,i,e) = Tstar_v
 crystallite_Fp(1:3,1:3,g,i,e) = Fp_new
 crystallite_Fe(1:3,1:3,g,i,e) = Fe_new
 crystallite_invFp(1:3,1:3,g,i,e) = invFp_new


 !* set return flag to true

 crystallite_integrateStress = .true.
#ifndef _OPENMP
 if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt &
     .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
             .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> P / MPa',math_transpose33(crystallite_P(1:3,1:3,g,i,e))*1.0e-6_pReal
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Cauchy / MPa', &
              math_mul33x33(crystallite_P(1:3,1:3,g,i,e), math_transpose33(Fg_new)) * 1.0e-6_pReal / math_det33(Fg_new)
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fe Lp Fe^-1', &
              math_transpose33(math_mul33x33(Fe_new, math_mul33x33(crystallite_Lp(1:3,1:3,g,i,e), math_inv33(Fe_new))))    ! transpose to get correct print out order
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fp',math_transpose33(crystallite_Fp(1:3,1:3,g,i,e))
 endif
#endif

 if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
   !$OMP CRITICAL (distributionStress)
    debug_StressLoopDistribution(NiterationStress,numerics_integrationMode) = &
      debug_StressLoopDistribution(NiterationStress,numerics_integrationMode) + 1_pInt
   !$OMP END CRITICAL (distributionStress)
 endif

end function crystallite_integrateStress


!--------------------------------------------------------------------------------------------------
!> @brief calculates orientations and disorientations (in case of single grain ips)
!--------------------------------------------------------------------------------------------------
subroutine crystallite_orientations
 use math, only: &
   math_pDecomposition, &
   math_RtoQ, &
   math_qConj
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP
 use IO, only: &
   IO_warning
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
 use constitutive_nonlocal, only: &
   constitutive_nonlocal_updateCompatibility


 implicit none
 integer(pInt) &
   e, &                          ! element index
   i, &                          ! integration point index
   g, &                          ! grain index
   n, &                          ! neighbor index
   neighboring_e, &              ! element index of my neighbor
   neighboring_i, &              ! integration point index of my neighbor
   myPhase, &                    ! phase
   neighboringPhase
 real(pReal), dimension(3,3) :: &
   U, &
   R
 real(pReal), dimension(4) :: &
   orientation
 logical &
   error

 ! --- CALCULATE ORIENTATION AND LATTICE ROTATION ---

 !$OMP PARALLEL DO PRIVATE(error,U,R,orientation)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1_pInt,homogenization_Ngrains(mesh_element(3,e))
! somehow this subroutine is not threadsafe, so need critical statement here; not clear, what exactly the problem is
         !$OMP CRITICAL (polarDecomp)
         call math_pDecomposition(crystallite_Fe(1:3,1:3,g,i,e), U, R, error)                       ! polar decomposition of Fe
         !$OMP END CRITICAL (polarDecomp)
         if (error) then
           call IO_warning(650_pInt, e, i, g)
           orientation = [1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal]                               ! fake orientation
         else
           orientation = math_RtoQ(transpose(R))
         endif
         crystallite_rotation(1:4,g,i,e) = lattice_qDisorientation(crystallite_orientation0(1:4,g,i,e), &  ! active rotation from ori0
                                                                orientation)                               ! to current orientation (with no symmetry)
         crystallite_orientation(1:4,g,i,e) = orientation
       enddo
     enddo
   enddo
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
                                            lattice_structure(myPhase))                              ! calculate disorientation for given symmetry
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

         call constitutive_nonlocal_updateCompatibility(crystallite_orientation,i,e)

       endif
     enddo
   enddo
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
   mesh_ipCoordinates, &
   mesh_maxNipNeighbors, &
   mesh_ipNeighborhood, &
   FE_NipNeighbors, &
   FE_geomtype, &
   FE_celltype
 use material, only: &
   plasticState, &
   damageState, &
   thermalState, &
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

#ifdef multiphysicsOut
 real(pReal), dimension(1+crystallite_sizePostResults(microstructure_crystallite(mesh_element(4,el))) + &
                        1+plasticState(material_phase(ipc,ip,el))%sizePostResults + &
                          damageState( material_phase(ipc,ip,el))%sizePostResults + &
                          thermalState(material_phase(ipc,ip,el))%sizePostResults) :: & 
   crystallite_postResults
#else
 real(pReal), dimension(1+crystallite_sizePostResults(microstructure_crystallite(mesh_element(4,el)))+ &
                        1+plasticState(material_phase(ipc,ip,el))%sizePostResults) :: &
   crystallite_postResults
#endif
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
     case (lp_ID)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = &
         reshape(math_transpose33(crystallite_Lp(1:3,1:3,ipc,ip,el)),[mySize])
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
