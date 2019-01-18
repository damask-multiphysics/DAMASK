!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Denny Tjahjanto, Max-Planck-Institut für Eisenforschung GmbH
!> @brief homogenization manager, organizing deformation partitioning and stress homogenization
!--------------------------------------------------------------------------------------------------
module homogenization
 use prec, only: &
   pInt, &
   pReal

!--------------------------------------------------------------------------------------------------
! General variables for the homogenization at a  material point
 implicit none
 private
   real(pReal),   dimension(:,:,:,:),   allocatable, public :: &
   materialpoint_F0, &                                                                              !< def grad of IP at start of FE increment
   materialpoint_F, &                                                                               !< def grad of IP to be reached at end of FE increment
   materialpoint_P                                                                                  !< first P--K stress of IP
 real(pReal),   dimension(:,:,:,:,:,:), allocatable, public ::  &
   materialpoint_dPdF                                                                               !< tangent of first P--K stress at IP
 real(pReal),   dimension(:,:,:),       allocatable, public :: &
   materialpoint_results                                                                            !< results array of material point
 integer(pInt),                                      public, protected  :: &
   materialpoint_sizeResults, &
   homogenization_maxSizePostResults, &
   thermal_maxSizePostResults, &
   damage_maxSizePostResults

 real(pReal),   dimension(:,:,:,:),     allocatable, private :: &
   materialpoint_subF0, &                                                                           !< def grad of IP at beginning of homogenization increment
   materialpoint_subF                                                                               !< def grad of IP to be reached at end of homog inc
 real(pReal),   dimension(:,:),         allocatable, private :: &
   materialpoint_subFrac, &
   materialpoint_subStep, &
   materialpoint_subdt
 logical,       dimension(:,:),         allocatable, private :: &
   materialpoint_requested, &
   materialpoint_converged
 logical,       dimension(:,:,:),       allocatable, private :: &
   materialpoint_doneAndHappy

 public ::  &
   homogenization_init, &
   materialpoint_stressAndItsTangent, &
   materialpoint_postResults
 private :: &
   homogenization_partitionDeformation, &
   homogenization_updateState, &
   homogenization_averageStressAndItsTangent, &
   homogenization_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!--------------------------------------------------------------------------------------------------
subroutine homogenization_init
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use math, only: &
   math_I3
 use debug, only: &
   debug_level, &
   debug_homogenization, &
   debug_levelBasic, &
   debug_e, &
   debug_g
 use mesh, only: &
   mesh_maxNips, &
   mesh_NcpElems, &
   mesh_element, &
   FE_Nips, &
   FE_geomtype
 use constitutive, only: &
   constitutive_plasticity_maxSizePostResults, &
   constitutive_source_maxSizePostResults
 use crystallite, only: &
   crystallite_maxSizePostResults
 use config, only: &
  material_configFile, &
  material_localFileExt, &
  config_deallocate, &
  config_homogenization, &
  homogenization_name
 use material
 use homogenization_none
 use homogenization_isostrain
 use homogenization_RGC
 use thermal_isothermal
 use thermal_adiabatic
 use thermal_conduction
 use damage_none
 use damage_local
 use damage_nonlocal
 use IO
 use numerics, only: &
   worldrank

 implicit none
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt) :: e,i,p
 integer(pInt), dimension(:,:), pointer :: thisSize
 integer(pInt), dimension(:)  , pointer :: thisNoutput
 character(len=64), dimension(:,:), pointer :: thisOutput
 character(len=32) :: outputName                                                                    !< name of output, intermediate fix until HDF5 output is ready
 logical :: valid


!--------------------------------------------------------------------------------------------------
! open material.config
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file

!--------------------------------------------------------------------------------------------------
! parse homogenization from config file 
 if (any(homogenization_type == HOMOGENIZATION_NONE_ID)) &
   call homogenization_none_init()
 if (any(homogenization_type == HOMOGENIZATION_ISOSTRAIN_ID)) &
   call homogenization_isostrain_init(FILEUNIT)
 if (any(homogenization_type == HOMOGENIZATION_RGC_ID)) &
   call homogenization_RGC_init(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! parse thermal from config file
 call IO_checkAndRewind(FILEUNIT)
 if (any(thermal_type == THERMAL_isothermal_ID)) &
   call thermal_isothermal_init()
 if (any(thermal_type == THERMAL_adiabatic_ID)) &
   call thermal_adiabatic_init(FILEUNIT)
 if (any(thermal_type == THERMAL_conduction_ID)) &
   call thermal_conduction_init(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! parse damage from config file
 call IO_checkAndRewind(FILEUNIT)
 if (any(damage_type == DAMAGE_none_ID)) &
   call damage_none_init()
 if (any(damage_type == DAMAGE_local_ID)) &
   call damage_local_init(FILEUNIT)
 if (any(damage_type == DAMAGE_nonlocal_ID)) &
   call damage_nonlocal_init(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! write description file for homogenization output
 mainProcess2: if (worldrank == 0) then
   call IO_write_jobFile(FILEUNIT,'outputHomogenization')
   do p = 1,size(config_homogenization)
     if (any(material_homog == p)) then
       i = homogenization_typeInstance(p)                                                               ! which instance of this homogenization type
       valid = .true.                                                                                   ! assume valid
       select case(homogenization_type(p))                                                              ! split per homogenization type
         case (HOMOGENIZATION_NONE_ID)
           outputName = HOMOGENIZATION_NONE_label
           thisNoutput => null()
           thisOutput => null()
           thisSize   => null()
         case (HOMOGENIZATION_ISOSTRAIN_ID)
           outputName = HOMOGENIZATION_ISOSTRAIN_label
           thisNoutput => homogenization_isostrain_Noutput
           thisOutput => homogenization_isostrain_output
           thisSize   => homogenization_isostrain_sizePostResult
         case (HOMOGENIZATION_RGC_ID)
           outputName = HOMOGENIZATION_RGC_label
           thisNoutput => homogenization_RGC_Noutput
           thisOutput => homogenization_RGC_output
           thisSize   => homogenization_RGC_sizePostResult
         case default
           valid = .false.
       end select
       write(FILEUNIT,'(/,a,/)')  '['//trim(homogenization_name(p))//']'
       if (valid) then
         write(FILEUNIT,'(a)') '(type)'//char(9)//trim(outputName)
         write(FILEUNIT,'(a,i4)') '(ngrains)'//char(9),homogenization_Ngrains(p)
         if (homogenization_type(p) /= HOMOGENIZATION_NONE_ID) then
           do e = 1,thisNoutput(i)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
           enddo
         endif
       endif
       i = thermal_typeInstance(p)                                                                      ! which instance of this thermal type
       valid = .true.                                                                                   ! assume valid
       select case(thermal_type(p))                                                                     ! split per thermal type
         case (THERMAL_isothermal_ID)
           outputName = THERMAL_isothermal_label
           thisNoutput => null()
           thisOutput => null()
           thisSize   => null()
         case (THERMAL_adiabatic_ID)
           outputName = THERMAL_adiabatic_label
           thisNoutput => thermal_adiabatic_Noutput
           thisOutput => thermal_adiabatic_output
           thisSize   => thermal_adiabatic_sizePostResult
         case (THERMAL_conduction_ID)
           outputName = THERMAL_conduction_label
           thisNoutput => thermal_conduction_Noutput
           thisOutput => thermal_conduction_output
           thisSize   => thermal_conduction_sizePostResult
         case default
           valid = .false.
       end select
       if (valid) then
         write(FILEUNIT,'(a)') '(thermal)'//char(9)//trim(outputName)
         if (thermal_type(p) /= THERMAL_isothermal_ID) then
           do e = 1,thisNoutput(i)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
           enddo
         endif
       endif
       i = damage_typeInstance(p)                                                                       ! which instance of this damage type
       valid = .true.                                                                                   ! assume valid
       select case(damage_type(p))                                                                      ! split per damage type
         case (DAMAGE_none_ID)
           outputName = DAMAGE_none_label
           thisNoutput => null()
           thisOutput => null()
           thisSize   => null()
         case (DAMAGE_local_ID)
           outputName = DAMAGE_local_label
           thisNoutput => damage_local_Noutput
           thisOutput => damage_local_output
           thisSize   => damage_local_sizePostResult
         case (DAMAGE_nonlocal_ID)
           outputName = DAMAGE_nonlocal_label
           thisNoutput => damage_nonlocal_Noutput
           thisOutput => damage_nonlocal_output
           thisSize   => damage_nonlocal_sizePostResult
         case default
           valid = .false.
       end select
       if (valid) then
         write(FILEUNIT,'(a)') '(damage)'//char(9)//trim(outputName)
         if (damage_type(p) /= DAMAGE_none_ID) then
           do e = 1,thisNoutput(i)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
           enddo
         endif
       endif
     endif
   enddo
   close(FILEUNIT)
 endif mainProcess2

 call config_deallocate('material.config/homogenization')

!--------------------------------------------------------------------------------------------------
! allocate and initialize global variables
 allocate(materialpoint_dPdF(3,3,3,3,mesh_maxNips,mesh_NcpElems),       source=0.0_pReal)
 allocate(materialpoint_F0(3,3,mesh_maxNips,mesh_NcpElems),             source=0.0_pReal)
 materialpoint_F0 = spread(spread(math_I3,3,mesh_maxNips),4,mesh_NcpElems)                          ! initialize to identity
 allocate(materialpoint_F(3,3,mesh_maxNips,mesh_NcpElems),              source=0.0_pReal)
 materialpoint_F = materialpoint_F0                                                                 ! initialize to identity
 allocate(materialpoint_subF0(3,3,mesh_maxNips,mesh_NcpElems),          source=0.0_pReal)
 allocate(materialpoint_subF(3,3,mesh_maxNips,mesh_NcpElems),           source=0.0_pReal)
 allocate(materialpoint_P(3,3,mesh_maxNips,mesh_NcpElems),              source=0.0_pReal)
 allocate(materialpoint_subFrac(mesh_maxNips,mesh_NcpElems),            source=0.0_pReal)
 allocate(materialpoint_subStep(mesh_maxNips,mesh_NcpElems),            source=0.0_pReal)
 allocate(materialpoint_subdt(mesh_maxNips,mesh_NcpElems),              source=0.0_pReal)
 allocate(materialpoint_requested(mesh_maxNips,mesh_NcpElems),          source=.false.)
 allocate(materialpoint_converged(mesh_maxNips,mesh_NcpElems),          source=.true.)
 allocate(materialpoint_doneAndHappy(2,mesh_maxNips,mesh_NcpElems),     source=.true.)

!--------------------------------------------------------------------------------------------------
! allocate and initialize global state and postresutls variables
 homogenization_maxSizePostResults = 0_pInt
 thermal_maxSizePostResults        = 0_pInt
 damage_maxSizePostResults         = 0_pInt
 do p = 1,size(config_homogenization)
   homogenization_maxSizePostResults = max(homogenization_maxSizePostResults,homogState       (p)%sizePostResults)
   thermal_maxSizePostResults        = max(thermal_maxSizePostResults,       thermalState     (p)%sizePostResults)
   damage_maxSizePostResults         = max(damage_maxSizePostResults        ,damageState      (p)%sizePostResults)
 enddo

 materialpoint_sizeResults = 1 &                                                                    ! grain count
                           + 1 + homogenization_maxSizePostResults &                                ! homogSize & homogResult
                               + thermal_maxSizePostResults        &
                               + damage_maxSizePostResults         &
                           + homogenization_maxNgrains * (1 + crystallite_maxSizePostResults &      ! crystallite size & crystallite results
                                                        + 1 + constitutive_plasticity_maxSizePostResults &     ! constitutive size & constitutive results
                                                            + constitutive_source_maxSizePostResults)
 allocate(materialpoint_results(materialpoint_sizeResults,mesh_maxNips,mesh_NcpElems))

 write(6,'(/,a)')   ' <<<+-  homogenization init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt) then
#ifdef TODO
   write(6,'(a32,1x,7(i8,1x))')   'homogenization_state0:          ', shape(homogenization_state0)
   write(6,'(a32,1x,7(i8,1x))')   'homogenization_subState0:       ', shape(homogenization_subState0)
   write(6,'(a32,1x,7(i8,1x))')   'homogenization_state:           ', shape(homogenization_state)
#endif
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_dPdF:             ', shape(materialpoint_dPdF)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_F0:               ', shape(materialpoint_F0)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_F:                ', shape(materialpoint_F)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subF0:            ', shape(materialpoint_subF0)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subF:             ', shape(materialpoint_subF)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_P:                ', shape(materialpoint_P)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subFrac:          ', shape(materialpoint_subFrac)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subStep:          ', shape(materialpoint_subStep)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subdt:            ', shape(materialpoint_subdt)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_requested:        ', shape(materialpoint_requested)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_converged:        ', shape(materialpoint_converged)
   write(6,'(a32,1x,7(i8,1x),/)') 'materialpoint_doneAndHappy:     ', shape(materialpoint_doneAndHappy)
   write(6,'(a32,1x,7(i8,1x))')   'maxSizePostResults: ', homogenization_maxSizePostResults
 endif
 flush(6)

 if (debug_g < 1 .or. debug_g > homogenization_Ngrains(mesh_element(3,debug_e))) &
   call IO_error(602_pInt,ext_msg='constituent', el=debug_e, g=debug_g)

end subroutine homogenization_init


!--------------------------------------------------------------------------------------------------
!> @brief  parallelized calculation of stress and corresponding tangent at material points
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_stressAndItsTangent(updateJaco,dt)
 use numerics, only: &
   subStepMinHomog, &
   subStepSizeHomog, &
   stepIncreaseHomog, &
   nMPstate
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP, &
   terminallyIll
 use mesh, only: &
   mesh_element
 use material, only: &
   plasticState, &
   sourceState, &
   homogState, &
   thermalState, &
   damageState, &
   phase_Nsources, &
   mappingHomogenization, &
   phaseAt, phasememberAt, &
   homogenization_Ngrains
 use crystallite, only: &
   crystallite_F0, &
   crystallite_Fp0, &
   crystallite_Fp, &
   crystallite_Fi0, &
   crystallite_Fi, &
   crystallite_Lp0, &
   crystallite_Lp, &
   crystallite_Li0, &
   crystallite_Li, &
   crystallite_dPdF, &
   crystallite_Tstar0_v, &
   crystallite_Tstar_v, &
   crystallite_partionedF0, &
   crystallite_partionedF, &
   crystallite_partionedFp0, &
   crystallite_partionedLp0, &
   crystallite_partionedFi0, &
   crystallite_partionedLi0, &
   crystallite_partionedTstar0_v, &
   crystallite_dt, &
   crystallite_requested, &
   crystallite_stress, &
   crystallite_stressTangent, &
   crystallite_orientations
#ifdef DEBUG
 use debug, only: &
   debug_level, &
   debug_homogenization, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_e, &
   debug_i
#endif

 implicit none
 real(pReal), intent(in) :: dt                                                                      !< time increment
 logical,     intent(in) :: updateJaco                                                              !< initiating Jacobian update
 integer(pInt) :: &
   NiterationHomog, &
   NiterationMPstate, &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e, &                                                                                             !< element number
   mySource, &
   myNgrains

#ifdef DEBUG
 if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt) then
     write(6,'(/a,i5,1x,i2)') '<< HOMOG >> Material Point start at el ip ', debug_e, debug_i

     write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< HOMOG >> F0', &
                                     transpose(materialpoint_F0(1:3,1:3,debug_i,debug_e))
     write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< HOMOG >> F', &
                                     transpose(materialpoint_F(1:3,1:3,debug_i,debug_e))
 endif
#endif

!--------------------------------------------------------------------------------------------------
! initialize restoration points of ...
 do e = FEsolving_execElem(1),FEsolving_execElem(2)
   myNgrains = homogenization_Ngrains(mesh_element(3,e))
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e); do g = 1,myNgrains

     plasticState    (phaseAt(g,i,e))%partionedState0(:,phasememberAt(g,i,e)) = &
     plasticState    (phaseAt(g,i,e))%state0(         :,phasememberAt(g,i,e))
     do mySource = 1_pInt, phase_Nsources(phaseAt(g,i,e))
       sourceState(phaseAt(g,i,e))%p(mySource)%partionedState0(:,phasememberAt(g,i,e)) = &
       sourceState(phaseAt(g,i,e))%p(mySource)%state0(         :,phasememberAt(g,i,e))
     enddo

     crystallite_partionedFp0(1:3,1:3,g,i,e) = crystallite_Fp0(1:3,1:3,g,i,e)                       ! ...plastic def grads
     crystallite_partionedLp0(1:3,1:3,g,i,e) = crystallite_Lp0(1:3,1:3,g,i,e)                       ! ...plastic velocity grads
     crystallite_partionedFi0(1:3,1:3,g,i,e) = crystallite_Fi0(1:3,1:3,g,i,e)                       ! ...intermediate def grads
     crystallite_partionedLi0(1:3,1:3,g,i,e) = crystallite_Li0(1:3,1:3,g,i,e)                       ! ...intermediate velocity grads
     crystallite_partionedF0(1:3,1:3,g,i,e) = crystallite_F0(1:3,1:3,g,i,e)                         ! ...def grads
     crystallite_partionedTstar0_v(1:6,g,i,e) = crystallite_Tstar0_v(1:6,g,i,e)                     ! ...2nd PK stress

   enddo; enddo
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e))
     materialpoint_subF0(1:3,1:3,i,e) = materialpoint_F0(1:3,1:3,i,e)                               ! ...def grad
     materialpoint_subFrac(i,e) = 0.0_pReal
     materialpoint_subStep(i,e) = 1.0_pReal/subStepSizeHomog                                        ! <<added to adopt flexibility in cutback size>>
     materialpoint_converged(i,e) = .false.                                                         ! pretend failed step of twice the required size
     materialpoint_requested(i,e) = .true.                                                          ! everybody requires calculation
   endforall
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
     homogState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
       homogState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
       homogState(mappingHomogenization(2,i,e))%State0(   :,mappingHomogenization(1,i,e))           ! ...internal homogenization state
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
     thermalState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
       thermalState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
       thermalState(mappingHomogenization(2,i,e))%State0(   :,mappingHomogenization(1,i,e))         ! ...internal thermal state
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
     damageState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
       damageState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
       damageState(mappingHomogenization(2,i,e))%State0(   :,mappingHomogenization(1,i,e))          ! ...internal damage state
 enddo
 NiterationHomog = 0_pInt

 cutBackLooping: do while (.not. terminallyIll .and. &
      any(materialpoint_subStep(:,FEsolving_execELem(1):FEsolving_execElem(2)) > subStepMinHomog))

   !$OMP PARALLEL DO PRIVATE(myNgrains)
   elementLooping1: do e = FEsolving_execElem(1),FEsolving_execElem(2)
     myNgrains = homogenization_Ngrains(mesh_element(3,e))
     IpLooping1: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)

       converged: if ( materialpoint_converged(i,e) ) then
#ifdef DEBUG
         if (iand(debug_level(debug_homogenization), debug_levelExtensive) /= 0_pInt &
            .and. ((e == debug_e .and. i == debug_i) &
                   .or. .not. iand(debug_level(debug_homogenization),debug_levelSelective) /= 0_pInt)) then
           write(6,'(a,1x,f12.8,1x,a,1x,f12.8,1x,a,i8,1x,i2/)') '<< HOMOG >> winding forward from', &
             materialpoint_subFrac(i,e), 'to current materialpoint_subFrac', &
             materialpoint_subFrac(i,e)+materialpoint_subStep(i,e),'in materialpoint_stressAndItsTangent at el ip',e,i
         endif
#endif

!---------------------------------------------------------------------------------------------------
! calculate new subStep and new subFrac
         materialpoint_subFrac(i,e) = materialpoint_subFrac(i,e) + materialpoint_subStep(i,e)
         materialpoint_subStep(i,e) = min(1.0_pReal-materialpoint_subFrac(i,e), &
                                          stepIncreaseHomog*materialpoint_subStep(i,e))             ! introduce flexibility for step increase/acceleration

         steppingNeeded: if (materialpoint_subStep(i,e) > subStepMinHomog) then

           ! wind forward grain starting point of...
           crystallite_partionedF0(1:3,1:3,1:myNgrains,i,e) =  &
              crystallite_partionedF(1:3,1:3,1:myNgrains,i,e)                                       ! ...def grads

           crystallite_partionedFp0(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_Fp(1:3,1:3,1:myNgrains,i,e)                                                ! ...plastic def grads

           crystallite_partionedLp0(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_Lp(1:3,1:3,1:myNgrains,i,e)                                                ! ...plastic velocity grads

           crystallite_partionedFi0(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_Fi(1:3,1:3,1:myNgrains,i,e)                                                ! ...intermediate def grads

           crystallite_partionedLi0(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_Li(1:3,1:3,1:myNgrains,i,e)                                                ! ...intermediate velocity grads

           crystallite_partionedTstar0_v(1:6,1:myNgrains,i,e) = &
             crystallite_Tstar_v(1:6,1:myNgrains,i,e)                                               ! ...2nd PK stress

           do g = 1,myNgrains
             plasticState    (phaseAt(g,i,e))%partionedState0(:,phasememberAt(g,i,e)) = &
             plasticState    (phaseAt(g,i,e))%state(          :,phasememberAt(g,i,e))
             do mySource = 1_pInt, phase_Nsources(phaseAt(g,i,e))
               sourceState(phaseAt(g,i,e))%p(mySource)%partionedState0(:,phasememberAt(g,i,e)) = &
               sourceState(phaseAt(g,i,e))%p(mySource)%state(          :,phasememberAt(g,i,e))
             enddo
           enddo

           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             homogState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               homogState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
               homogState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e))   ! ...internal homogenization state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             thermalState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               thermalState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
               thermalState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e)) ! ...internal thermal state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             damageState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               damageState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
               damageState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e))  ! ...internal damage state
           materialpoint_subF0(1:3,1:3,i,e) = materialpoint_subF(1:3,1:3,i,e)                       ! ...def grad
         endif steppingNeeded

       else converged
         if ( (myNgrains == 1_pInt .and. materialpoint_subStep(i,e) <= 1.0 ) .or. &                 ! single grain already tried internal subStepping in crystallite
              subStepSizeHomog * materialpoint_subStep(i,e) <=  subStepMinHomog ) then              ! would require too small subStep
                                                                                                    ! cutback makes no sense
           !$OMP FLUSH(terminallyIll)
           if (.not. terminallyIll) then                                                            ! so first signals terminally ill...
             !$OMP CRITICAL (write2out)
               write(6,*) 'Integration point ', i,' at element ', e, ' terminally ill'
             !$OMP END CRITICAL (write2out)
           endif
           !$OMP CRITICAL (setTerminallyIll)
             terminallyIll = .true.                                                                 ! ...and kills all others
           !$OMP END CRITICAL (setTerminallyIll)
         else                                                                                       ! cutback makes sense
           materialpoint_subStep(i,e) = subStepSizeHomog * materialpoint_subStep(i,e)               ! crystallite had severe trouble, so do a significant cutback

#ifdef DEBUG
           if (iand(debug_level(debug_homogenization), debug_levelExtensive) /= 0_pInt &
              .and. ((e == debug_e .and. i == debug_i) &
                    .or. .not. iand(debug_level(debug_homogenization), debug_levelSelective) /= 0_pInt)) then
             write(6,'(a,1x,f12.8,a,i8,1x,i2/)') &
               '<< HOMOG >> cutback step in materialpoint_stressAndItsTangent with new materialpoint_subStep:',&
               materialpoint_subStep(i,e),' at el ip',e,i
           endif
#endif

!--------------------------------------------------------------------------------------------------
! restore...
           crystallite_Fp(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_partionedFp0(1:3,1:3,1:myNgrains,i,e)                                      ! ...plastic def grads
           crystallite_Lp(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_partionedLp0(1:3,1:3,1:myNgrains,i,e)                                      ! ...plastic velocity grads
           crystallite_Fi(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_partionedFi0(1:3,1:3,1:myNgrains,i,e)                                      ! ...intermediate def grads
           crystallite_Li(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_partionedLi0(1:3,1:3,1:myNgrains,i,e)                                      ! ...intermediate velocity grads
           crystallite_Tstar_v(1:6,1:myNgrains,i,e) = &
              crystallite_partionedTstar0_v(1:6,1:myNgrains,i,e)                                    ! ...2nd PK stress
           do g = 1, myNgrains
             plasticState    (phaseAt(g,i,e))%state(          :,phasememberAt(g,i,e)) = &
             plasticState    (phaseAt(g,i,e))%partionedState0(:,phasememberAt(g,i,e))
             do mySource = 1_pInt, phase_Nsources(phaseAt(g,i,e))
               sourceState(phaseAt(g,i,e))%p(mySource)%state(          :,phasememberAt(g,i,e)) = &
               sourceState(phaseAt(g,i,e))%p(mySource)%partionedState0(:,phasememberAt(g,i,e))
             enddo
           enddo
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             homogState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               homogState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e)) = &
               homogState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e))   ! ...internal homogenization state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             thermalState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               thermalState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e)) = &
               thermalState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) ! ...internal thermal state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             damageState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               damageState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e)) = &
               damageState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e))  ! ...internal damage state
         endif
       endif converged

       if (materialpoint_subStep(i,e) > subStepMinHomog) then
         materialpoint_requested(i,e) = .true.
         materialpoint_subF(1:3,1:3,i,e) = materialpoint_subF0(1:3,1:3,i,e) &
                                         + materialpoint_subStep(i,e) * (materialpoint_F(1:3,1:3,i,e) &
                                         - materialpoint_F0(1:3,1:3,i,e))
         materialpoint_subdt(i,e) = materialpoint_subStep(i,e) * dt
         materialpoint_doneAndHappy(1:2,i,e) = [.false.,.true.]
       endif
     enddo IpLooping1
   enddo elementLooping1
   !$OMP END PARALLEL DO

   NiterationMPstate = 0_pInt

   convergenceLooping: do while (.not. terminallyIll .and. &
             any(            materialpoint_requested(:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                 .and. .not. materialpoint_doneAndHappy(1,:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                ) .and. &
             NiterationMPstate < nMPstate)
     NiterationMPstate = NiterationMPstate + 1

!--------------------------------------------------------------------------------------------------
! deformation partitioning
! based on materialpoint_subF0,.._subF,crystallite_partionedF0, and homogenization_state,
! results in crystallite_partionedF
     !$OMP PARALLEL DO PRIVATE(myNgrains)
     elementLooping2: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       myNgrains = homogenization_Ngrains(mesh_element(3,e))
       IpLooping2: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
         if (      materialpoint_requested(i,e) .and. &                                             ! process requested but...
             .not. materialpoint_doneAndHappy(1,i,e)) then                                          ! ...not yet done material points
           call homogenization_partitionDeformation(i,e)                                            ! partition deformation onto constituents
           crystallite_dt(1:myNgrains,i,e) = materialpoint_subdt(i,e)                               ! propagate materialpoint dt to grains
           crystallite_requested(1:myNgrains,i,e) = .true.                                          ! request calculation for constituents
         else
           crystallite_requested(1:myNgrains,i,e) = .false.                                         ! calculation for constituents not required anymore
         endif
       enddo IpLooping2
     enddo elementLooping2
     !$OMP END PARALLEL DO

!--------------------------------------------------------------------------------------------------
! crystallite integration
! based on crystallite_partionedF0,.._partionedF
! incrementing by crystallite_dt
     materialpoint_converged = crystallite_stress() !ToDo: MD not sure if that is the best logic

!--------------------------------------------------------------------------------------------------
! state update
     !$OMP PARALLEL DO
     elementLooping3: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       IpLooping3: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
         if (      materialpoint_requested(i,e) .and. &
             .not. materialpoint_doneAndHappy(1,i,e)) then
           if (.not. materialpoint_converged(i,e)) then
             materialpoint_doneAndHappy(1:2,i,e) = [.true.,.false.]
           else
             materialpoint_doneAndHappy(1:2,i,e) = homogenization_updateState(i,e)
             materialpoint_converged(i,e) = all(materialpoint_doneAndHappy(1:2,i,e))                  ! converged if done and happy
           endif
         endif
       enddo IpLooping3
     enddo elementLooping3
     !$OMP END PARALLEL DO

   enddo convergenceLooping

   NiterationHomog = NiterationHomog + 1_pInt

 enddo cutBackLooping
 
 if(updateJaco) call crystallite_stressTangent

 if (.not. terminallyIll ) then
   call crystallite_orientations()                                                                  ! calculate crystal orientations
   !$OMP PARALLEL DO
   elementLooping4: do e = FEsolving_execElem(1),FEsolving_execElem(2)
     IpLooping4: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       call homogenization_averageStressAndItsTangent(i,e)
     enddo IpLooping4
   enddo elementLooping4
   !$OMP END PARALLEL DO
 else
   write(6,'(/,a,/)') '<< HOMOG >> Material Point terminally ill'
 endif

end subroutine materialpoint_stressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief parallelized calculation of result array at material points
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_postResults
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP
 use mesh, only: &
   mesh_element
 use material, only: &
   mappingHomogenization, &
   homogState, &
   thermalState, &
   damageState, &
   plasticState, &
   sourceState, &
   material_phase, &
   homogenization_Ngrains, &
   microstructure_crystallite
 use crystallite, only: &
   crystallite_sizePostResults, &
   crystallite_postResults

 implicit none
 integer(pInt) :: &
   thePos, &
   theSize, &
   myNgrains, &
   myCrystallite, &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number

 !$OMP PARALLEL DO PRIVATE(myNgrains,myCrystallite,thePos,theSize)
   elementLooping: do e = FEsolving_execElem(1),FEsolving_execElem(2)
     myNgrains = homogenization_Ngrains(mesh_element(3,e))
     myCrystallite = microstructure_crystallite(mesh_element(4,e))
     IpLooping: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       thePos = 0_pInt

       theSize = homogState       (mappingHomogenization(2,i,e))%sizePostResults &
               + thermalState     (mappingHomogenization(2,i,e))%sizePostResults &
               + damageState      (mappingHomogenization(2,i,e))%sizePostResults
       materialpoint_results(thePos+1,i,e) = real(theSize,pReal)                                    ! tell size of homogenization results
       thePos = thePos + 1_pInt

       if (theSize > 0_pInt) then                                                                   ! any homogenization results to mention?
         materialpoint_results(thePos+1:thePos+theSize,i,e) = homogenization_postResults(i,e)       ! tell homogenization results
         thePos = thePos + theSize
       endif

       materialpoint_results(thePos+1,i,e) = real(myNgrains,pReal)                                  ! tell number of grains at materialpoint
       thePos = thePos + 1_pInt

       grainLooping :do g = 1,myNgrains
         theSize = 1 + crystallite_sizePostResults(myCrystallite) + &
                   1 + plasticState    (material_phase(g,i,e))%sizePostResults + &                    !ToDo
                       sum(sourceState(material_phase(g,i,e))%p(:)%sizePostResults)
         materialpoint_results(thePos+1:thePos+theSize,i,e) = crystallite_postResults(g,i,e)        ! tell crystallite results
         thePos = thePos + theSize
       enddo grainLooping
     enddo IpLooping
   enddo elementLooping
 !$OMP END PARALLEL DO

end subroutine materialpoint_postResults


!--------------------------------------------------------------------------------------------------
!> @brief  partition material point def grad onto constituents
!--------------------------------------------------------------------------------------------------
subroutine homogenization_partitionDeformation(ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_type, &
   homogenization_maxNgrains, &
   HOMOGENIZATION_NONE_ID, &
   HOMOGENIZATION_ISOSTRAIN_ID, &
   HOMOGENIZATION_RGC_ID
 use crystallite, only: &
   crystallite_partionedF
 use homogenization_isostrain, only: &
   homogenization_isostrain_partitionDeformation
 use homogenization_RGC, only: &
   homogenization_RGC_partitionDeformation

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number

 chosenHomogenization: select case(homogenization_type(mesh_element(3,el)))

   case (HOMOGENIZATION_NONE_ID) chosenHomogenization
     crystallite_partionedF(1:3,1:3,1:homogenization_maxNgrains,ip,el) = 0.0_pReal
     crystallite_partionedF(1:3,1:3,1:1,ip,el) = &
       spread(materialpoint_subF(1:3,1:3,ip,el),3,1)

   case (HOMOGENIZATION_ISOSTRAIN_ID) chosenHomogenization
     call homogenization_isostrain_partitionDeformation(&
                          crystallite_partionedF(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                          materialpoint_subF(1:3,1:3,ip,el),&
                          el)
   case (HOMOGENIZATION_RGC_ID) chosenHomogenization
     call homogenization_RGC_partitionDeformation(&
                         crystallite_partionedF(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                         materialpoint_subF(1:3,1:3,ip,el),&
                         ip, &
                         el)
 end select chosenHomogenization

end subroutine homogenization_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief update the internal state of the homogenization scheme and tell whether "done" and
!> "happy" with result
!--------------------------------------------------------------------------------------------------
function homogenization_updateState(ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_type, &
   thermal_type, &
   damage_type, &
   homogenization_maxNgrains, &
   HOMOGENIZATION_RGC_ID, &
   THERMAL_adiabatic_ID, &
   DAMAGE_local_ID
 use crystallite, only: &
   crystallite_P, &
   crystallite_dPdF, &
   crystallite_partionedF,&
   crystallite_partionedF0
 use homogenization_RGC, only: &
   homogenization_RGC_updateState
 use thermal_adiabatic, only: &
   thermal_adiabatic_updateState
 use damage_local, only: &
   damage_local_updateState

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number
 logical, dimension(2) :: homogenization_updateState

 homogenization_updateState = .true.
 chosenHomogenization: select case(homogenization_type(mesh_element(3,el)))
   case (HOMOGENIZATION_RGC_ID) chosenHomogenization
     homogenization_updateState = &
       homogenization_updateState .and. &
        homogenization_RGC_updateState(crystallite_P(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                                       crystallite_partionedF(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                                       crystallite_partionedF0(1:3,1:3,1:homogenization_maxNgrains,ip,el),&
                                       materialpoint_subF(1:3,1:3,ip,el),&
                                       materialpoint_subdt(ip,el), &
                                       crystallite_dPdF(1:3,1:3,1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                                       ip, &
                                       el)
 end select chosenHomogenization

 chosenThermal: select case (thermal_type(mesh_element(3,el)))
   case (THERMAL_adiabatic_ID) chosenThermal
     homogenization_updateState = &
       homogenization_updateState .and. &
       thermal_adiabatic_updateState(materialpoint_subdt(ip,el), &
                                     ip, &
                                     el)
 end select chosenThermal

 chosenDamage: select case (damage_type(mesh_element(3,el)))
   case (DAMAGE_local_ID) chosenDamage
     homogenization_updateState = &
       homogenization_updateState .and. &
       damage_local_updateState(materialpoint_subdt(ip,el), &
                                ip, &
                                el)
 end select chosenDamage

end function homogenization_updateState


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities
!--------------------------------------------------------------------------------------------------
subroutine homogenization_averageStressAndItsTangent(ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_type, &
   homogenization_maxNgrains, &
   HOMOGENIZATION_NONE_ID, &
   HOMOGENIZATION_ISOSTRAIN_ID, &
   HOMOGENIZATION_RGC_ID
 use crystallite, only: &
   crystallite_P,crystallite_dPdF
 use homogenization_isostrain, only: &
   homogenization_isostrain_averageStressAndItsTangent
 use homogenization_RGC, only: &
   homogenization_RGC_averageStressAndItsTangent

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number

 chosenHomogenization: select case(homogenization_type(mesh_element(3,el)))
   case (HOMOGENIZATION_NONE_ID) chosenHomogenization
       materialpoint_P(1:3,1:3,ip,el) = sum(crystallite_P(1:3,1:3,1:1,ip,el),3)
       materialpoint_dPdF(1:3,1:3,1:3,1:3,ip,el) &
        = sum(crystallite_dPdF(1:3,1:3,1:3,1:3,1:1,ip,el),5)

   case (HOMOGENIZATION_ISOSTRAIN_ID) chosenHomogenization
     call homogenization_isostrain_averageStressAndItsTangent(&
       materialpoint_P(1:3,1:3,ip,el), &
       materialpoint_dPdF(1:3,1:3,1:3,1:3,ip,el),&
       crystallite_P(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
       crystallite_dPdF(1:3,1:3,1:3,1:3,1:homogenization_maxNgrains,ip,el), &
       el)
   case (HOMOGENIZATION_RGC_ID) chosenHomogenization
     call homogenization_RGC_averageStressAndItsTangent(&
       materialpoint_P(1:3,1:3,ip,el), &
       materialpoint_dPdF(1:3,1:3,1:3,1:3,ip,el),&
       crystallite_P(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
       crystallite_dPdF(1:3,1:3,1:3,1:3,1:homogenization_maxNgrains,ip,el), &
       el)
 end select chosenHomogenization

end subroutine homogenization_averageStressAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion. call only,
!> if homogenization_sizePostResults(i,e) > 0 !!
!--------------------------------------------------------------------------------------------------
function homogenization_postResults(ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   mappingHomogenization, &
   homogState, &
   thermalState, &
   damageState, &
   homogenization_type, &
   thermal_type, &
   damage_type, &
   HOMOGENIZATION_NONE_ID, &
   HOMOGENIZATION_ISOSTRAIN_ID, &
   HOMOGENIZATION_RGC_ID, &
   THERMAL_isothermal_ID, &
   THERMAL_adiabatic_ID, &
   THERMAL_conduction_ID, &
   DAMAGE_none_ID, &
   DAMAGE_local_ID, &
   DAMAGE_nonlocal_ID
 use homogenization_isostrain, only: &
   homogenization_isostrain_postResults
 use homogenization_RGC, only: &
   homogenization_RGC_postResults
 use thermal_adiabatic, only: &
   thermal_adiabatic_postResults
 use thermal_conduction, only: &
   thermal_conduction_postResults
 use damage_local, only: &
   damage_local_postResults
 use damage_nonlocal, only: &
   damage_nonlocal_postResults

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number
 real(pReal), dimension(  homogState       (mappingHomogenization(2,ip,el))%sizePostResults &
                        + thermalState     (mappingHomogenization(2,ip,el))%sizePostResults &
                        + damageState      (mappingHomogenization(2,ip,el))%sizePostResults) :: &
   homogenization_postResults
 integer(pInt) :: &
   startPos, endPos

 homogenization_postResults = 0.0_pReal

 startPos = 1_pInt
 endPos   = homogState(mappingHomogenization(2,ip,el))%sizePostResults
 chosenHomogenization: select case (homogenization_type(mesh_element(3,el)))
   case (HOMOGENIZATION_NONE_ID) chosenHomogenization

   case (HOMOGENIZATION_ISOSTRAIN_ID) chosenHomogenization
     homogenization_postResults(startPos:endPos) = &
       homogenization_isostrain_postResults(&
                                  ip, &
                                  el, &
                                  materialpoint_P(1:3,1:3,ip,el), &
                                  materialpoint_F(1:3,1:3,ip,el))
   case (HOMOGENIZATION_RGC_ID) chosenHomogenization
     homogenization_postResults(startPos:endPos) = &
       homogenization_RGC_postResults(&
                                  ip, &
                                  el, &
                                  materialpoint_P(1:3,1:3,ip,el), &
                                  materialpoint_F(1:3,1:3,ip,el))
 end select chosenHomogenization

 startPos = endPos + 1_pInt
 endPos   = endPos + thermalState(mappingHomogenization(2,ip,el))%sizePostResults
 chosenThermal: select case (thermal_type(mesh_element(3,el)))
   case (THERMAL_isothermal_ID) chosenThermal

   case (THERMAL_adiabatic_ID) chosenThermal
     homogenization_postResults(startPos:endPos) = &
       thermal_adiabatic_postResults(ip, el)
   case (THERMAL_conduction_ID) chosenThermal
     homogenization_postResults(startPos:endPos) = &
       thermal_conduction_postResults(ip, el)
 end select chosenThermal

 startPos = endPos + 1_pInt
 endPos   = endPos + damageState(mappingHomogenization(2,ip,el))%sizePostResults
 chosenDamage: select case (damage_type(mesh_element(3,el)))
   case (DAMAGE_none_ID) chosenDamage

   case (DAMAGE_local_ID) chosenDamage
     homogenization_postResults(startPos:endPos) = &
       damage_local_postResults(ip, el)

   case (DAMAGE_nonlocal_ID) chosenDamage
     homogenization_postResults(startPos:endPos) = &
       damage_nonlocal_postResults(ip, el)
 end select chosenDamage

end function homogenization_postResults

end module homogenization
