!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief elasticity, plasticity, internal microstructure state
!--------------------------------------------------------------------------------------------------
module constitutive
 use prec, only: &
   pInt

 implicit none
 private
 integer(pInt), public, protected :: &
   constitutive_plasticity_maxSizePostResults, &
   constitutive_plasticity_maxSizeDotState, &
   constitutive_source_maxSizePostResults, &
   constitutive_source_maxSizeDotState

 public :: &
   constitutive_init, &
   constitutive_homogenizedC, &
   constitutive_microstructure, &
   constitutive_LpAndItsTangent, &
   constitutive_LiAndItsTangent, &
   constitutive_initialFi, &
   constitutive_TandItsTangent, &
   constitutive_collectDotState, &
   constitutive_collectDeltaState, &
   constitutive_postResults

 private :: &
   constitutive_hooke_TandItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates arrays pointing to array of the various constitutive modules
!--------------------------------------------------------------------------------------------------
subroutine constitutive_init()
#ifdef HDF
 use hdf5, only: &
   HID_T
 use IO, only : &
   HDF5_mappingConstitutive
#endif

 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use prec, only: &
   pReal
 use debug, only: &
   debug_constitutive, &
   debug_levelBasic
 use numerics, only: &
   worldrank
 use IO, only: &
   IO_error, &
   IO_open_file, &
   IO_open_jobFile_stat, &
   IO_write_jobFile, &
   IO_write_jobIntFile, &
   IO_timeStamp
 use mesh, only: &
   FE_geomtype
 use material, only: &
   material_phase, &
   material_Nphase, &
   material_localFileExt, &
   material_configFile, &
   phase_name, &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Nsources, &
   phase_source, &
   phase_kinematics, &
   ELASTICITY_hooke_ID, &
   PLASTICITY_none_ID, &
   PLASTICITY_isotropic_ID, &
   PLASTICITY_j2_ID, &
   PLASTICITY_phenopowerlaw_ID, &
   PLASTICITY_phenoplus_ID, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_titanmod_ID, &
   PLASTICITY_nonlocal_ID ,&
   SOURCE_thermal_dissipation_ID, &
   SOURCE_thermal_externalheat_ID, &
   SOURCE_damage_isoBrittle_ID, &
   SOURCE_damage_isoDuctile_ID, &
   SOURCE_damage_anisoBrittle_ID, &
   SOURCE_damage_anisoDuctile_ID, &
   SOURCE_vacancy_phenoplasticity_ID, &
   SOURCE_vacancy_irradiation_ID, &
   SOURCE_vacancy_thermalfluc_ID, &
   KINEMATICS_cleavage_opening_ID, &
   KINEMATICS_slipplane_opening_ID, &
   KINEMATICS_thermal_expansion_ID, &
   KINEMATICS_vacancy_strain_ID, &
   KINEMATICS_hydrogen_strain_ID, &
   ELASTICITY_HOOKE_label, &
   PLASTICITY_NONE_label, &
   PLASTICITY_ISOTROPIC_label, &
   PLASTICITY_J2_label, &
   PLASTICITY_PHENOPOWERLAW_label, &
   PLASTICITY_PHENOPLUS_label, &
   PLASTICITY_DISLOTWIN_label, &
   PLASTICITY_DISLOUCLA_label, &
   PLASTICITY_TITANMOD_label, &
   PLASTICITY_NONLOCAL_label, &
   SOURCE_thermal_dissipation_label, &
   SOURCE_thermal_externalheat_label, &
   SOURCE_damage_isoBrittle_label, &
   SOURCE_damage_isoDuctile_label, &
   SOURCE_damage_anisoBrittle_label, &
   SOURCE_damage_anisoDuctile_label, &
   SOURCE_vacancy_phenoplasticity_label, &
   SOURCE_vacancy_irradiation_label, &
   SOURCE_vacancy_thermalfluc_label, &
   plasticState, &
   sourceState

 use plastic_none
 use plastic_isotropic
 use plastic_j2
 use plastic_phenopowerlaw
 use plastic_phenoplus
 use plastic_dislotwin
 use plastic_disloucla
 use plastic_titanmod
 use plastic_nonlocal
 use source_thermal_dissipation
 use source_thermal_externalheat
 use source_damage_isoBrittle
 use source_damage_isoDuctile
 use source_damage_anisoBrittle
 use source_damage_anisoDuctile
 use source_vacancy_phenoplasticity
 use source_vacancy_irradiation
 use source_vacancy_thermalfluc
 use kinematics_cleavage_opening
 use kinematics_slipplane_opening
 use kinematics_thermal_expansion
 use kinematics_vacancy_strain
 use kinematics_hydrogen_strain

 implicit none
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt) :: &
  e, &                                                                                              !< maximum number of elements
  phase, &
  mySource, &
  instance

 integer(pInt), dimension(:,:), pointer :: thisSize
 integer(pInt), dimension(:)  , pointer :: thisNoutput
 character(len=64), dimension(:,:), pointer :: thisOutput
 character(len=32) :: outputName                                                                    !< name of output, intermediate fix until HDF5 output is ready
 logical :: knownPlasticity, knownSource, nonlocalConstitutionPresent
 nonlocalConstitutionPresent = .false.

!--------------------------------------------------------------------------------------------------
! parse plasticities from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(phase_plasticity == PLASTICITY_NONE_ID))          call plastic_none_init
 if (any(phase_plasticity == PLASTICITY_ISOTROPIC_ID))     call plastic_isotropic_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_J2_ID))            call plastic_j2_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_PHENOPOWERLAW_ID)) call plastic_phenopowerlaw_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_PHENOPLUS_ID))     call plastic_phenoplus_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_DISLOTWIN_ID))     call plastic_dislotwin_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_DISLOUCLA_ID))     call plastic_disloucla_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_TITANMOD_ID))      call plastic_titanmod_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_NONLOCAL_ID)) then
  call plastic_nonlocal_init(FILEUNIT)
  call plastic_nonlocal_stateInit()
 endif
 close(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! parse source mechanisms from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(phase_source == SOURCE_thermal_dissipation_ID))     call source_thermal_dissipation_init(FILEUNIT)
 if (any(phase_source == SOURCE_thermal_externalheat_ID))    call source_thermal_externalheat_init(FILEUNIT)
 if (any(phase_source == SOURCE_damage_isoBrittle_ID))       call source_damage_isoBrittle_init(FILEUNIT)
 if (any(phase_source == SOURCE_damage_isoDuctile_ID))       call source_damage_isoDuctile_init(FILEUNIT)
 if (any(phase_source == SOURCE_damage_anisoBrittle_ID))     call source_damage_anisoBrittle_init(FILEUNIT)
 if (any(phase_source == SOURCE_damage_anisoDuctile_ID))     call source_damage_anisoDuctile_init(FILEUNIT)
 if (any(phase_source == SOURCE_vacancy_phenoplasticity_ID)) call source_vacancy_phenoplasticity_init(FILEUNIT)
 if (any(phase_source == SOURCE_vacancy_irradiation_ID))     call source_vacancy_irradiation_init(FILEUNIT)
 if (any(phase_source == SOURCE_vacancy_thermalfluc_ID))     call source_vacancy_thermalfluc_init(FILEUNIT)
 close(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! parse kinematic mechanisms from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(phase_kinematics == KINEMATICS_cleavage_opening_ID))  call kinematics_cleavage_opening_init(FILEUNIT)
 if (any(phase_kinematics == KINEMATICS_slipplane_opening_ID)) call kinematics_slipplane_opening_init(FILEUNIT)
 if (any(phase_kinematics == KINEMATICS_thermal_expansion_ID)) call kinematics_thermal_expansion_init(FILEUNIT)
 if (any(phase_kinematics == KINEMATICS_vacancy_strain_ID))    call kinematics_vacancy_strain_init(FILEUNIT)
 if (any(phase_kinematics == KINEMATICS_hydrogen_strain_ID))   call kinematics_hydrogen_strain_init(FILEUNIT)
 close(FILEUNIT)

 mainProcess: if (worldrank == 0) then
   write(6,'(/,a)')   ' <<<+-  constitutive init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

!--------------------------------------------------------------------------------------------------
! write description file for constitutive phase output
 if (worldrank == 0_pInt) then
   call IO_write_jobFile(FILEUNIT,'outputConstitutive')
   do phase = 1_pInt,material_Nphase
     if (any(material_phase == phase)) then                                                             ! is this phase active?
       instance = phase_plasticityInstance(phase)                                                       ! which instance of a plasticity is present phase
       knownPlasticity = .true.                                                                         ! assume valid
       select case(phase_plasticity(phase))                                                             ! split per constititution
         case (PLASTICITY_NONE_ID)
           outputName = PLASTICITY_NONE_label
           thisNoutput => null()
           thisOutput => null()                                                                         ! plastic_none_output
           thisSize   => null()                                                                         ! plastic_none_sizePostResult
         case (PLASTICITY_ISOTROPIC_ID)
           outputName = PLASTICITY_ISOTROPIC_label
           thisNoutput => plastic_isotropic_Noutput
           thisOutput => plastic_isotropic_output
           thisSize   => plastic_isotropic_sizePostResult
         case (PLASTICITY_J2_ID)
           outputName = PLASTICITY_J2_label
           thisNoutput => plastic_j2_Noutput
           thisOutput => plastic_j2_output
           thisSize   => plastic_j2_sizePostResult
         case (PLASTICITY_PHENOPOWERLAW_ID)
           outputName = PLASTICITY_PHENOPOWERLAW_label
           thisNoutput => plastic_phenopowerlaw_Noutput
           thisOutput => plastic_phenopowerlaw_output
           thisSize   => plastic_phenopowerlaw_sizePostResult
         case (PLASTICITY_PHENOPLUS_ID)
           outputName = PLASTICITY_PHENOPLUS_label
           thisNoutput => plastic_phenoplus_Noutput
           thisOutput => plastic_phenoplus_output
           thisSize   => plastic_phenoplus_sizePostResult
         case (PLASTICITY_DISLOTWIN_ID)
           outputName = PLASTICITY_DISLOTWIN_label
           thisNoutput => plastic_dislotwin_Noutput
           thisOutput => plastic_dislotwin_output
           thisSize   => plastic_dislotwin_sizePostResult
         case (PLASTICITY_DISLOUCLA_ID)
           outputName = PLASTICITY_DISLOUCLA_label
           thisNoutput => plastic_disloucla_Noutput
           thisOutput => plastic_disloucla_output
           thisSize   => plastic_disloucla_sizePostResult
         case (PLASTICITY_TITANMOD_ID)
           outputName = PLASTICITY_TITANMOD_label
           thisNoutput => plastic_titanmod_Noutput
           thisOutput => plastic_titanmod_output
           thisSize   => plastic_titanmod_sizePostResult
         case (PLASTICITY_NONLOCAL_ID)
           outputName = PLASTICITY_NONLOCAL_label
           thisNoutput => plastic_nonlocal_Noutput
           thisOutput => plastic_nonlocal_output
           thisSize   => plastic_nonlocal_sizePostResult
         case default
           knownPlasticity = .false.
       end select
       write(FILEUNIT,'(/,a,/)') '['//trim(phase_name(phase))//']'
       if (knownPlasticity) then
         write(FILEUNIT,'(a)') '(plasticity)'//char(9)//trim(outputName)
         if (phase_plasticity(phase) /= PLASTICITY_NONE_ID) then
           do e = 1_pInt,thisNoutput(instance)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(e,instance))//char(9),thisSize(e,instance)
           enddo
         endif
       endif
       do mySource = 1_pInt, phase_Nsources(phase)
         knownSource = .true.
         select case (phase_source(mySource,phase))
           case (SOURCE_thermal_dissipation_ID)
             instance = source_thermal_dissipation_instance(phase)
             outputName = SOURCE_thermal_dissipation_label
             thisNoutput => source_thermal_dissipation_Noutput
             thisOutput => source_thermal_dissipation_output
             thisSize   => source_thermal_dissipation_sizePostResult
           case (SOURCE_thermal_externalheat_ID)
             instance = source_thermal_externalheat_instance(phase)
             outputName = SOURCE_thermal_externalheat_label
             thisNoutput => source_thermal_externalheat_Noutput
             thisOutput => source_thermal_externalheat_output
             thisSize   => source_thermal_externalheat_sizePostResult
           case (SOURCE_damage_isoBrittle_ID)
             instance = source_damage_isoBrittle_instance(phase)
             outputName = SOURCE_damage_isoBrittle_label
             thisNoutput => source_damage_isoBrittle_Noutput
             thisOutput => source_damage_isoBrittle_output
             thisSize   => source_damage_isoBrittle_sizePostResult
           case (SOURCE_damage_isoDuctile_ID)
             instance = source_damage_isoDuctile_instance(phase)
             outputName = SOURCE_damage_isoDuctile_label
             thisNoutput => source_damage_isoDuctile_Noutput
             thisOutput => source_damage_isoDuctile_output
             thisSize   => source_damage_isoDuctile_sizePostResult
           case (SOURCE_damage_anisoBrittle_ID)
             instance = source_damage_anisoBrittle_instance(phase)
             outputName = SOURCE_damage_anisoBrittle_label
             thisNoutput => source_damage_anisoBrittle_Noutput
             thisOutput => source_damage_anisoBrittle_output
             thisSize   => source_damage_anisoBrittle_sizePostResult
           case (SOURCE_damage_anisoDuctile_ID)
             instance = source_damage_anisoDuctile_instance(phase)
             outputName = SOURCE_damage_anisoDuctile_label
             thisNoutput => source_damage_anisoDuctile_Noutput
             thisOutput => source_damage_anisoDuctile_output
             thisSize   => source_damage_anisoDuctile_sizePostResult
           case (SOURCE_vacancy_phenoplasticity_ID)
             instance = source_vacancy_phenoplasticity_instance(phase)
             outputName = SOURCE_vacancy_phenoplasticity_label
             thisNoutput => source_vacancy_phenoplasticity_Noutput
             thisOutput => source_vacancy_phenoplasticity_output
             thisSize   => source_vacancy_phenoplasticity_sizePostResult
           case (SOURCE_vacancy_irradiation_ID)
             instance = source_vacancy_irradiation_instance(phase)
             outputName = SOURCE_vacancy_irradiation_label
             thisNoutput => source_vacancy_irradiation_Noutput
             thisOutput => source_vacancy_irradiation_output
             thisSize   => source_vacancy_irradiation_sizePostResult
           case (SOURCE_vacancy_thermalfluc_ID)
             instance = source_vacancy_thermalfluc_instance(phase)
             outputName = SOURCE_vacancy_thermalfluc_label
             thisNoutput => source_vacancy_thermalfluc_Noutput
             thisOutput => source_vacancy_thermalfluc_output
             thisSize   => source_vacancy_thermalfluc_sizePostResult
           case default
             knownSource = .false.
         end select
         if (knownSource) then
           write(FILEUNIT,'(a)') '(source)'//char(9)//trim(outputName)
           do e = 1_pInt,thisNoutput(instance)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(e,instance))//char(9),thisSize(e,instance)
           enddo
         endif
       enddo
     endif
   enddo
   close(FILEUNIT)
 endif

 constitutive_plasticity_maxSizeDotState = 0_pInt
 constitutive_plasticity_maxSizePostResults = 0_pInt
 constitutive_source_maxSizeDotState = 0_pInt
 constitutive_source_maxSizePostResults = 0_pInt

 PhaseLoop2:do phase = 1_pInt,material_Nphase
   plasticState    (phase)%partionedState0 = plasticState    (phase)%State0
   plasticState    (phase)%State           = plasticState    (phase)%State0
   forall(mySource = 1_pInt:phase_Nsources(phase)) &
     sourceState(phase)%p(mySource)%partionedState0 = sourceState(phase)%p(mySource)%State0
   forall(mySource = 1_pInt:phase_Nsources(phase)) &
     sourceState(phase)%p(mySource)%State           = sourceState(phase)%p(mySource)%State0

   constitutive_plasticity_maxSizeDotState    = max(constitutive_plasticity_maxSizeDotState,    &
                                                    plasticState(phase)%sizeDotState)
   constitutive_plasticity_maxSizePostResults = max(constitutive_plasticity_maxSizePostResults, &
                                                    plasticState(phase)%sizePostResults)
   constitutive_source_maxSizeDotState        = max(constitutive_source_maxSizeDotState, &
                                                    maxval(sourceState(phase)%p(:)%sizeDotState))
   constitutive_source_maxSizePostResults     = max(constitutive_source_maxSizePostResults, &
                                                    maxval(sourceState(phase)%p(:)%sizePostResults))
 enddo PhaseLoop2

#ifdef HDF
 call  HDF5_mappingConstitutive(mappingConstitutive)
 do phase = 1_pInt,material_Nphase
   instance = phase_plasticityInstance(phase)                                                       ! which instance of a plasticity is present phase
   select case(phase_plasticity(phase))                                                             ! split per constititution
     case (PLASTICITY_NONE_ID)
     case (PLASTICITY_ISOTROPIC_ID)
     case (PLASTICITY_J2_ID)
   end select
 enddo
#endif

#ifdef TODO
!--------------------------------------------------------------------------------------------------
! report
 constitutive_maxSizeState       = maxval(constitutive_sizeState)
 constitutive_plasticity_maxSizeDotState    = maxval(constitutive_sizeDotState)

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) then
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_state0:          ', shape(constitutive_state0)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_partionedState0: ', shape(constitutive_partionedState0)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_subState0:       ', shape(constitutive_subState0)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_state:           ', shape(constitutive_state)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_aTolState:       ', shape(constitutive_aTolState)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_dotState:        ', shape(constitutive_dotState)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_deltaState:      ', shape(constitutive_deltaState)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_sizeState:       ', shape(constitutive_sizeState)
   write(6,'(a32,1x,7(i8,1x))')   'constitutive_sizeDotState:    ', shape(constitutive_sizeDotState)
   write(6,'(a32,1x,7(i8,1x),/)') 'constitutive_sizePostResults: ', shape(constitutive_sizePostResults)
   write(6,'(a32,1x,7(i8,1x))')   'maxSizeState:       ', constitutive_maxSizeState
   write(6,'(a32,1x,7(i8,1x))')   'maxSizeDotState:    ', constitutive_plasticity_maxSizeDotState
   write(6,'(a32,1x,7(i8,1x))')   'maxSizePostResults: ', constitutive_plasticity_maxSizePostResults
 endif
 flush(6)
#endif


end subroutine constitutive_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenize elasticity matrix
!--------------------------------------------------------------------------------------------------
function constitutive_homogenizedC(ipc,ip,el)
 use prec, only: &
   pReal
 use material, only: &
   phase_plasticity, &
   material_phase, &
   PLASTICITY_TITANMOD_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOUCLA_ID
 use plastic_titanmod, only: &
   plastic_titanmod_homogenizedC
 use plastic_dislotwin, only: &
   plastic_dislotwin_homogenizedC
 use plastic_disloucla, only: &
   plastic_disloucla_homogenizedC
 use lattice, only: &
   lattice_C66

 implicit none
 real(pReal), dimension(6,6) :: constitutive_homogenizedC
 integer(pInt), intent(in) :: &
   ipc, &                                                                                            !< grain number
   ip, &                                                                                             !< integration point number
   el                                                                                                !< element number

 select case (phase_plasticity(material_phase(ipc,ip,el)))

   case (PLASTICITY_DISLOTWIN_ID)
     constitutive_homogenizedC = plastic_dislotwin_homogenizedC(ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID)
     constitutive_homogenizedC = plastic_disloucla_homogenizedC(ipc,ip,el)
   case (PLASTICITY_TITANMOD_ID)
     constitutive_homogenizedC = plastic_titanmod_homogenizedC (ipc,ip,el)
   case default
     constitutive_homogenizedC = lattice_C66(1:6,1:6,material_phase (ipc,ip,el))

 end select

end function constitutive_homogenizedC

!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different constitutive models
!--------------------------------------------------------------------------------------------------
subroutine constitutive_microstructure(orientations, Fe, Fp, ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   phase_plasticity, &
   material_phase, &
   material_homog, &
   temperature, &
   thermalMapping, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_titanmod_ID, &
   PLASTICITY_nonlocal_ID, &
   PLASTICITY_phenoplus_ID
 use plastic_titanmod, only: &
   plastic_titanmod_microstructure
 use plastic_nonlocal, only: &
   plastic_nonlocal_microstructure
 use plastic_dislotwin, only: &
   plastic_dislotwin_microstructure
 use plastic_disloucla, only: &
   plastic_disloucla_microstructure
 use plastic_phenoplus, only: &
   plastic_phenoplus_microstructure

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in), dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fp                                                                                               !< plastic deformation gradient
 integer(pInt) :: &
   phase, homog, offset
 real(pReal),   intent(in), dimension(:,:,:,:) :: &
   orientations                                                                                     !< crystal orientation in quaternions

 phase = material_phase(ipc,ip,el)
 homog = material_homog(    ip,el)
 offset = thermalMapping(homog)%p(ip,el)
 select case (phase_plasticity(phase))

   case (PLASTICITY_DISLOTWIN_ID)
     call plastic_dislotwin_microstructure(temperature(homog)%p(offset),ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID)
     call plastic_disloucla_microstructure(temperature(homog)%p(offset),ipc,ip,el)
   case (PLASTICITY_TITANMOD_ID)
     call plastic_titanmod_microstructure (temperature(homog)%p(offset),ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID)
     call plastic_nonlocal_microstructure (Fe,Fp,ip,el)
   case (PLASTICITY_PHENOPLUS_ID)
     call plastic_phenoplus_microstructure(orientations,ipc,ip,el)

 end select

end subroutine constitutive_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LpAndItsTangent(Lp, dLp_dTstar3333, dLp_dFi3333, Tstar_v, Fi, ipc, ip, el)
 use prec, only: &
   pReal
 use math, only: &
   math_transpose33, &
   math_mul33x33, &
   math_Mandel6to33, &
   math_Mandel33to6, &
   math_Plain99to3333
 use material, only: &
   phase_plasticity, &
   material_phase, &
   material_homog, &
   temperature, &
   thermalMapping, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_ISOTROPIC_ID, &
   PLASTICITY_J2_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_PHENOPLUS_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOUCLA_ID, &
   PLASTICITY_TITANMOD_ID, &
   PLASTICITY_NONLOCAL_ID
 use plastic_isotropic, only: &
   plastic_isotropic_LpAndItsTangent
 use plastic_j2, only: &
   plastic_j2_LpAndItsTangent
 use plastic_phenopowerlaw, only: &
   plastic_phenopowerlaw_LpAndItsTangent
 use plastic_phenoplus, only: &
   plastic_phenoplus_LpAndItsTangent
 use plastic_dislotwin, only: &
   plastic_dislotwin_LpAndItsTangent
 use plastic_disloucla, only: &
   plastic_disloucla_LpAndItsTangent
 use plastic_titanmod, only: &
   plastic_titanmod_LpAndItsTangent
 use plastic_nonlocal, only: &
   plastic_nonlocal_LpAndItsTangent

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dLp_dTstar3333, &                                                                                !< derivative of Lp with respect to Tstar (4th-order tensor)
   dLp_dFi3333                                                                                      !< derivative of Lp with respect to Fi (4th-order tensor)
 real(pReal), dimension(6) :: &
   Mstar_v                                                                                          !< Mandel stress work conjugate with Lp
 real(pReal), dimension(9,9) :: &
   dLp_dMstar                                                                                       !< derivative of Lp with respect to Mstar (4th-order tensor)
 real(pReal), dimension(3,3) :: &
   temp_33
 integer(pInt) :: &
   i, j, phase, homog, offset

 phase = material_phase(ipc,ip,el)
 homog = material_homog(    ip,el)
 offset = thermalMapping(homog)%p(ip,el)
 Mstar_v = math_Mandel33to6(math_mul33x33(math_mul33x33(math_transpose33(Fi),Fi), &
                            math_Mandel6to33(Tstar_v)))
 select case (phase_plasticity(phase))

   case (PLASTICITY_NONE_ID)
     Lp = 0.0_pReal
     dLp_dMstar = 0.0_pReal
   case (PLASTICITY_ISOTROPIC_ID)
     call plastic_isotropic_LpAndItsTangent(Lp,dLp_dMstar,Mstar_v,ipc,ip,el)
   case (PLASTICITY_J2_ID)
     call plastic_j2_LpAndItsTangent(Lp,dLp_dMstar,Mstar_v,ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID)
     call plastic_phenopowerlaw_LpAndItsTangent(Lp,dLp_dMstar,Mstar_v,ipc,ip,el)
   case (PLASTICITY_PHENOPLUS_ID)
     call plastic_phenoplus_LpAndItsTangent(Lp,dLp_dMstar,Mstar_v,ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID)
     call plastic_nonlocal_LpAndItsTangent(Lp,dLp_dMstar,Mstar_v, &
                                           temperature(homog)%p(offset), &
                                           ip,el)
   case (PLASTICITY_DISLOTWIN_ID)
     call plastic_dislotwin_LpAndItsTangent(Lp,dLp_dMstar,Mstar_v, &
                                            temperature(homog)%p(offset), &
                                            ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID)
     call plastic_disloucla_LpAndItsTangent(Lp,dLp_dMstar,Mstar_v, &
                                            temperature(homog)%p(offset), &
                                            ipc,ip,el)
   case (PLASTICITY_TITANMOD_ID)
     call plastic_titanmod_LpAndItsTangent(Lp,dLp_dMstar,Mstar_v, &
                                           temperature(homog)%p(offset), &
                                           ipc,ip,el)

 end select

 dLp_dTstar3333 = math_Plain99to3333(dLp_dMstar)
 temp_33 = math_mul33x33(Fi,math_Mandel6to33(Tstar_v))
 do i = 1_pInt, 3_pInt; do j = 1_pInt, 3_pInt
   dLp_dFi3333(i,j,1:3,1:3) = math_mul33x33(temp_33,math_transpose33(dLp_dTstar3333(i,j,1:3,1:3))) + &
                              math_mul33x33(math_mul33x33(Fi,dLp_dTstar3333(i,j,1:3,1:3)),math_Mandel6to33(Tstar_v))
 enddo; enddo
 temp_33 = math_mul33x33(math_transpose33(Fi),Fi)
 do i = 1_pInt, 3_pInt; do j = 1_pInt, 3_pInt
   dLp_dTstar3333(i,j,1:3,1:3) = math_mul33x33(temp_33,dLp_dTstar3333(i,j,1:3,1:3))
 enddo; enddo

end subroutine constitutive_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LiAndItsTangent(Li, dLi_dTstar3333, dLi_dFi3333, Tstar_v, Fi, ipc, ip, el)
 use prec, only: &
   pReal
 use math, only: &
   math_I3, &
   math_inv33, &
   math_det33, &
   math_transpose33, &
   math_mul33x33
 use material, only: &
   phase_plasticity, &
   material_phase, &
   material_homog, &
   mappingConstitutive, &
   phase_kinematics, &
   phase_Nkinematics, &
   PLASTICITY_isotropic_ID, &
   KINEMATICS_cleavage_opening_ID, &
   KINEMATICS_slipplane_opening_ID, &
   KINEMATICS_thermal_expansion_ID, &
   KINEMATICS_vacancy_strain_ID, &
   KINEMATICS_hydrogen_strain_ID
 use plastic_isotropic, only: &
   plastic_isotropic_LiAndItsTangent
 use kinematics_cleavage_opening, only: &
   kinematics_cleavage_opening_LiAndItsTangent
 use kinematics_slipplane_opening, only: &
   kinematics_slipplane_opening_LiAndItsTangent
 use kinematics_thermal_expansion, only: &
   kinematics_thermal_expansion_LiAndItsTangent
 use kinematics_vacancy_strain, only: &
   kinematics_vacancy_strain_LiAndItsTangent
 use kinematics_hydrogen_strain, only: &
   kinematics_hydrogen_strain_LiAndItsTangent

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   Li                                                                                               !< intermediate velocity gradient
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dLi_dTstar3333, &                                                                                !< derivative of Li with respect to Tstar (4th-order tensor)
   dLi_dFi3333
 real(pReal), dimension(3,3) :: &
   my_Li                                                                                            !< intermediate velocity gradient
 real(pReal), dimension(3,3,3,3) :: &
   my_dLi_dTstar
 real(pReal), dimension(3,3) :: &
   FiInv, &
   temp_33
 real(pReal) :: &
   detFi
 integer(pInt) :: &
   i, j, kinematics, phase, homog

 phase = material_phase(ipc,ip,el)
 homog = material_homog(    ip,el)

 Li = 0.0_pReal
 dLi_dTstar3333  = 0.0_pReal
 dLi_dFi3333     = 0.0_pReal

 select case (phase_plasticity(phase))

   case (PLASTICITY_isotropic_ID)
     call plastic_isotropic_LiAndItsTangent(my_Li, my_dLi_dTstar, Tstar_v, ipc, ip, el)
   
   case default
     my_Li = 0.0_pReal
     my_dLi_dTstar = 0.0_pReal
 end select
 Li = Li + my_Li
 dLi_dTstar3333 = dLi_dTstar3333 + my_dLi_dTstar

 do kinematics = 1_pInt, phase_Nkinematics(material_phase(ipc,ip,el))
   select case (phase_kinematics(kinematics,material_phase(ipc,ip,el)))
     case (KINEMATICS_cleavage_opening_ID)
       call kinematics_cleavage_opening_LiAndItsTangent(my_Li, my_dLi_dTstar, Tstar_v, ipc, ip, el)

     case (KINEMATICS_slipplane_opening_ID)
       call kinematics_slipplane_opening_LiAndItsTangent(my_Li, my_dLi_dTstar, Tstar_v, ipc, ip, el)

     case (KINEMATICS_thermal_expansion_ID)
       call kinematics_thermal_expansion_LiAndItsTangent(my_Li, my_dLi_dTstar, ipc, ip, el)

     case (KINEMATICS_vacancy_strain_ID)
       call kinematics_vacancy_strain_LiAndItsTangent(my_Li, my_dLi_dTstar, ipc, ip, el)

     case (KINEMATICS_hydrogen_strain_ID)
       call kinematics_hydrogen_strain_LiAndItsTangent(my_Li, my_dLi_dTstar, ipc, ip, el)

     case default
       my_Li = 0.0_pReal
       my_dLi_dTstar = 0.0_pReal
   end select
   Li = Li + my_Li
   dLi_dTstar3333 = dLi_dTstar3333 + my_dLi_dTstar
 enddo

 FiInv = math_inv33(Fi)
 detFi = math_det33(Fi)
 Li = math_mul33x33(math_mul33x33(Fi,Li),FiInv)*detFi                                               !< push forward to intermediate configuration
 temp_33 = math_mul33x33(FiInv,Li)
 do i = 1_pInt, 3_pInt; do j = 1_pInt, 3_pInt
   dLi_dTstar3333(1:3,1:3,i,j) = math_mul33x33(math_mul33x33(Fi,dLi_dTstar3333(1:3,1:3,i,j)),FiInv)*detFi
   dLi_dFi3333   (1:3,1:3,i,j) = dLi_dFi3333(1:3,1:3,i,j) + Li*FiInv(j,i)
   dLi_dFi3333   (1:3,i,1:3,j) = dLi_dFi3333(1:3,i,1:3,j) + math_I3*temp_33(j,i) + Li*FiInv(j,i)
 enddo; enddo

end subroutine constitutive_LiAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief  collects initial intermediate deformation gradient
!--------------------------------------------------------------------------------------------------
pure function constitutive_initialFi(ipc, ip, el)
 use prec, only: &
   pReal
 use math, only: &
   math_I3, &
   math_inv33, &
   math_mul33x33
 use material, only: &
   phase_kinematics, &
   phase_Nkinematics, &
   material_phase, &
   KINEMATICS_thermal_expansion_ID, &
   KINEMATICS_vacancy_strain_ID, &
   KINEMATICS_hydrogen_strain_ID
 use kinematics_thermal_expansion, only: &
   kinematics_thermal_expansion_initialStrain
 use kinematics_vacancy_strain, only: &
   kinematics_vacancy_strain_initialStrain
 use kinematics_hydrogen_strain, only: &
   kinematics_hydrogen_strain_initialStrain

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   constitutive_initialFi                                                                           !< composite initial intermediate deformation gradient
 integer(pInt) :: &
   kinematics

 constitutive_initialFi = math_I3

 do kinematics = 1_pInt, phase_Nkinematics(material_phase(ipc,ip,el))                               !< Warning: small initial strain assumption
   select case (phase_kinematics(kinematics,material_phase(ipc,ip,el)))
     case (KINEMATICS_thermal_expansion_ID)
       constitutive_initialFi = &
         constitutive_initialFi + kinematics_thermal_expansion_initialStrain(ipc, ip, el)

     case (KINEMATICS_vacancy_strain_ID)
       constitutive_initialFi = &
         constitutive_initialFi + kinematics_vacancy_strain_initialStrain(ipc, ip, el)

     case (KINEMATICS_hydrogen_strain_ID)
       constitutive_initialFi = &
         constitutive_initialFi + kinematics_hydrogen_strain_initialStrain(ipc, ip, el)

   end select
enddo

end function constitutive_initialFi


!--------------------------------------------------------------------------------------------------
!> @brief returns the 2nd Piola-Kirchhoff stress tensor and its tangent with respect to
!> the elastic deformation gradient depending on the selected elastic law (so far no case switch
!! because only hooke is implemented
!--------------------------------------------------------------------------------------------------
subroutine constitutive_TandItsTangent(T, dT_dFe, dT_dFi, Fe, Fi, ipc, ip, el)
 use prec, only: &
   pReal

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   T                                                                                                !< 2nd Piola-Kirchhoff stress tensor
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dT_dFe, &                                                                                        !< derivative of 2nd P-K stress with respect to elastic deformation gradient
   dT_dFi                                                                                           !< derivative of 2nd P-K stress with respect to intermediate deformation gradient

 call constitutive_hooke_TandItsTangent(T, dT_dFe, dT_dFi, Fe, Fi, ipc, ip, el)


end subroutine constitutive_TandItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief returns the 2nd Piola-Kirchhoff stress tensor and its tangent with respect to
!> the elastic deformation gradient using hookes law
!--------------------------------------------------------------------------------------------------
subroutine constitutive_hooke_TandItsTangent(T, dT_dFe, dT_dFi, Fe, Fi, ipc, ip, el)
 use prec, only: &
   pReal
 use math, only : &
   math_mul3x3, &
   math_mul33x33, &
   math_mul3333xx33, &
   math_Mandel66to3333, &
   math_transpose33, &
   math_trace33, &
   math_I3
 use material, only: &
   material_phase, &
   material_homog, &
   phase_NstiffnessDegradations, &
   phase_stiffnessDegradation, &
   damage, &
   damageMapping, &
   porosity, &
   porosityMapping, &
   STIFFNESS_DEGRADATION_damage_ID, &
   STIFFNESS_DEGRADATION_porosity_ID

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   T                                                                                                !< 2nd Piola-Kirchhoff stress tensor in lattice configuration
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dT_dFe, &                                                                                        !< derivative of 2nd P-K stress with respect to elastic deformation gradient
   dT_dFi                                                                                           !< derivative of 2nd P-K stress with respect to intermediate deformation gradient

 integer(pInt) :: i, j, phase, homog
 real(pReal), dimension(3,3) :: E
 real(pReal), dimension(3,3,3,3) :: C

 phase = material_phase(ipc,ip,el)
 homog = material_homog(ip,el)
 C = math_Mandel66to3333(constitutive_homogenizedC(ipc,ip,el))
 do i = 1_pInt, phase_NstiffnessDegradations(phase)
   select case(phase_stiffnessDegradation(i,phase))
     case (STIFFNESS_DEGRADATION_damage_ID)
       C = damage(homog)%p(damageMapping(homog)%p(ip,el))* &
           damage(homog)%p(damageMapping(homog)%p(ip,el))* &
           C

     case (STIFFNESS_DEGRADATION_porosity_ID)
       C = porosity(homog)%p(porosityMapping(homog)%p(ip,el))* &
           porosity(homog)%p(porosityMapping(homog)%p(ip,el))* &
           C
   end select
 enddo

 E = 0.5_pReal*(math_mul33x33(math_transpose33(Fe),Fe)-math_I3)                                     !< Green-Lagrange strain in unloaded configuration
 T = math_mul3333xx33(C,math_mul33x33(math_mul33x33(math_transpose33(Fi),E),Fi))                    !< 2PK stress in lattice configuration in work conjugate with GL strain pulled back to lattice configuration

 dT_dFe = 0.0_pReal
 forall (i=1_pInt:3_pInt, j=1_pInt:3_pInt)
   dT_dFe(i,j,1:3,1:3) = &
     math_mul33x33(Fe,math_mul33x33(math_mul33x33(Fi,C(i,j,1:3,1:3)),math_transpose33(Fi)))         !< dT_ij/dFe_kl = C_ijmn * Fi_lm * Fi_on * Fe_ko
   dT_dFi(i,j,1:3,1:3) = 2.0_pReal*math_mul33x33(math_mul33x33(E,Fi),C(i,j,1:3,1:3))                !< dT_ij/dFi_kl = C_ijln * E_km * Fe_mn
 end forall

end subroutine constitutive_hooke_TandItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine constitutive_collectDotState(Tstar_v, FeArray, FpArray, subdt, subfracArray,ipc, ip, el)
 use prec, only: &
   pReal, &
   pLongInt
 use debug, only: &
   debug_cumDotStateCalls, &
   debug_cumDotStateTicks, &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   phase_plasticity, &
   phase_source, &
   phase_Nsources, &
   material_phase, &
   material_homog, &
   temperature, &
   thermalMapping, &
   homogenization_maxNgrains, &
   PLASTICITY_none_ID, &
   PLASTICITY_isotropic_ID, &
   PLASTICITY_j2_ID, &
   PLASTICITY_phenopowerlaw_ID, &
   PLASTICITY_phenoplus_ID, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_titanmod_ID, &
   PLASTICITY_nonlocal_ID, &
   SOURCE_damage_isoDuctile_ID, &
   SOURCE_damage_anisoBrittle_ID, &
   SOURCE_damage_anisoDuctile_ID, &
   SOURCE_thermal_externalheat_ID
 use plastic_isotropic, only:  &
   plastic_isotropic_dotState
 use plastic_j2, only:  &
   plastic_j2_dotState
 use plastic_phenopowerlaw, only: &
   plastic_phenopowerlaw_dotState
 use plastic_phenoplus, only: &
   plastic_phenoplus_dotState
 use plastic_dislotwin, only: &
   plastic_dislotwin_dotState
 use plastic_disloucla, only: &
   plastic_disloucla_dotState
 use plastic_titanmod, only: &
   plastic_titanmod_dotState
 use plastic_nonlocal, only: &
   plastic_nonlocal_dotState
 use source_damage_isoDuctile, only: &
   source_damage_isoDuctile_dotState
 use source_damage_anisoBrittle, only: &
   source_damage_anisoBrittle_dotState
 use source_damage_anisoDuctile, only: &
   source_damage_anisoDuctile_dotState
 use source_thermal_externalheat, only: &
   source_thermal_externalheat_dotState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  intent(in) :: &
   subdt                                                                                            !< timestep
 real(pReal),  intent(in), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   subfracArray                                                                                     !< subfraction of timestep
 real(pReal),  intent(in), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   FeArray, &                                                                                       !< elastic deformation gradient
   FpArray                                                                                          !< plastic deformation gradient
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 integer(pLongInt) :: &
   tick, tock, &
   tickrate, &
   maxticks
 integer(pInt) :: &
   phase, homog, offset, mySource

 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) &
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)

 phase = material_phase(ipc,ip,el)
 homog = material_homog(    ip,el)
 offset = thermalMapping(homog)%p(ip,el)
 select case (phase_plasticity(phase))
   case (PLASTICITY_ISOTROPIC_ID)
     call plastic_isotropic_dotState           (Tstar_v,ipc,ip,el)
   case (PLASTICITY_J2_ID)
     call plastic_j2_dotState           (Tstar_v,ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID)
     call plastic_phenopowerlaw_dotState(Tstar_v,ipc,ip,el)
   case (PLASTICITY_PHENOPLUS_ID)
     call plastic_phenoplus_dotState(Tstar_v,ipc,ip,el)
   case (PLASTICITY_DISLOTWIN_ID)
     call plastic_dislotwin_dotState    (Tstar_v,temperature(homog)%p(offset), &
                                         ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID)
     call plastic_disloucla_dotState    (Tstar_v,temperature(homog)%p(offset), &
                                         ipc,ip,el)
   case (PLASTICITY_TITANMOD_ID)
     call plastic_titanmod_dotState     (Tstar_v,temperature(homog)%p(offset), &
                                         ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID)
     call plastic_nonlocal_dotState     (Tstar_v,FeArray,FpArray,temperature(homog)%p(offset), &
                                         subdt,subfracArray,ip,el)
 end select

 do mySource = 1_pInt, phase_Nsources(phase)
   select case (phase_source(mySource,phase))
     case (SOURCE_damage_anisoBrittle_ID)
       call source_damage_anisoBrittle_dotState (Tstar_v, ipc, ip, el)
     case (SOURCE_damage_isoDuctile_ID)
       call source_damage_isoDuctile_dotState   (         ipc, ip, el)
     case (SOURCE_damage_anisoDuctile_ID)
       call source_damage_anisoDuctile_dotState (         ipc, ip, el)
     case (SOURCE_thermal_externalheat_ID)
       call source_thermal_externalheat_dotState(         ipc, ip, el)

   end select
 enddo

 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) then
   call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
   !$OMP CRITICAL (debugTimingDotState)
     debug_cumDotStateCalls = debug_cumDotStateCalls + 1_pInt
     debug_cumDotStateTicks = debug_cumDotStateTicks + tock-tick
     !$OMP FLUSH (debug_cumDotStateTicks)
     if (tock < tick) debug_cumDotStateTicks  = debug_cumDotStateTicks + maxticks
   !$OMP END CRITICAL (debugTimingDotState)
 endif
end subroutine constitutive_collectDotState

!--------------------------------------------------------------------------------------------------
!> @brief for constitutive models having an instantaneous change of state (so far, only nonlocal)
!> will return false if delta state is not needed/supported by the constitutive model
!--------------------------------------------------------------------------------------------------
subroutine constitutive_collectDeltaState(Tstar_v, Fe, ipc, ip, el)
 use prec, only: &
   pReal, &
   pLongInt
 use debug, only: &
   debug_cumDeltaStateCalls, &
   debug_cumDeltaStateTicks, &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use material, only: &
   phase_plasticity, &
   phase_source, &
   phase_Nsources, &
   material_phase, &
   PLASTICITY_NONLOCAL_ID, &
   SOURCE_damage_isoBrittle_ID, &
   SOURCE_vacancy_irradiation_ID, &
   SOURCE_vacancy_thermalfluc_ID
 use plastic_nonlocal, only: &
   plastic_nonlocal_deltaState
 use source_damage_isoBrittle, only: &
   source_damage_isoBrittle_deltaState
 use source_vacancy_irradiation, only: &
   source_vacancy_irradiation_deltaState
 use source_vacancy_thermalfluc, only: &
   source_vacancy_thermalfluc_deltaState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(in), dimension(3,3) :: &
   Fe                                                                                               !< elastic deformation gradient
 integer(pInt) :: &
   mySource
 integer(pLongInt) :: &
   tick, tock, &
   tickrate, &
   maxticks

 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) &
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)

 select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_NONLOCAL_ID)
     call plastic_nonlocal_deltaState(Tstar_v,ip,el)

 end select

 do mySource = 1_pInt, phase_Nsources(material_phase(ipc,ip,el))
   select case (phase_source(mySource,material_phase(ipc,ip,el)))
     case (SOURCE_damage_isoBrittle_ID)
       call source_damage_isoBrittle_deltaState  (constitutive_homogenizedC(ipc,ip,el), Fe, &
                                                  ipc, ip, el)
     case (SOURCE_vacancy_irradiation_ID)
       call source_vacancy_irradiation_deltaState(ipc, ip, el)
     case (SOURCE_vacancy_thermalfluc_ID)
       call source_vacancy_thermalfluc_deltaState(ipc, ip, el)

   end select
 enddo

 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) then
   call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
   !$OMP CRITICAL (debugTimingDeltaState)
     debug_cumDeltaStateCalls = debug_cumDeltaStateCalls + 1_pInt
     debug_cumDeltaStateTicks = debug_cumDeltaStateTicks + tock-tick
     !$OMP FLUSH (debug_cumDeltaStateTicks)
     if (tock < tick) debug_cumDeltaStateTicks  = debug_cumDeltaStateTicks + maxticks
   !$OMP END CRITICAL (debugTimingDeltaState)
 endif

end subroutine constitutive_collectDeltaState


!--------------------------------------------------------------------------------------------------
!> @brief returns array of constitutive results
!--------------------------------------------------------------------------------------------------
function constitutive_postResults(Tstar_v, FeArray, ipc, ip, el)
 use prec, only: &
   pReal
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   plasticState, &
   sourceState, &
   phase_plasticity, &
   phase_source, &
   phase_Nsources, &
   material_phase, &
   material_homog, &
   temperature, &
   thermalMapping, &
   homogenization_maxNgrains, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_ISOTROPIC_ID, &
   PLASTICITY_J2_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_PHENOPLUS_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOUCLA_ID, &
   PLASTICITY_TITANMOD_ID, &
   PLASTICITY_NONLOCAL_ID, &
   SOURCE_damage_isoBrittle_ID, &
   SOURCE_damage_isoDuctile_ID, &
   SOURCE_damage_anisoBrittle_ID, &
   SOURCE_damage_anisoDuctile_ID
 use plastic_isotropic, only: &
   plastic_isotropic_postResults
 use plastic_j2, only: &
   plastic_j2_postResults
 use plastic_phenopowerlaw, only: &
   plastic_phenopowerlaw_postResults
 use plastic_phenoplus, only: &
   plastic_phenoplus_postResults
 use plastic_dislotwin, only: &
   plastic_dislotwin_postResults
 use plastic_disloucla, only: &
   plastic_disloucla_postResults
 use plastic_titanmod, only: &
   plastic_titanmod_postResults
 use plastic_nonlocal, only: &
   plastic_nonlocal_postResults
 use source_damage_isoBrittle, only: &
   source_damage_isoBrittle_postResults
 use source_damage_isoDuctile, only: &
   source_damage_isoDuctile_postResults
 use source_damage_anisoBrittle, only: &
   source_damage_anisoBrittle_postResults
 use source_damage_anisoDuctile, only: &
   source_damage_anisoDuctile_postResults

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%sizePostResults + &
                        sum(sourceState(material_phase(ipc,ip,el))%p(:)%sizePostResults)) :: &
   constitutive_postResults
 real(pReal),  intent(in), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   FeArray                                                                                          !< elastic deformation gradient
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 integer(pInt) :: &
   startPos, endPos, phase, homog, offset, mySource

 constitutive_postResults = 0.0_pReal

 phase = material_phase(ipc,ip,el)
 homog = material_homog(    ip,el)
 offset = thermalMapping(homog)%p(ip,el)

 startPos = 1_pInt
 endPos = plasticState(material_phase(ipc,ip,el))%sizePostResults
 select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_TITANMOD_ID)
     constitutive_postResults(startPos:endPos) = plastic_titanmod_postResults(ipc,ip,el)
   case (PLASTICITY_ISOTROPIC_ID)
     constitutive_postResults(startPos:endPos) = plastic_isotropic_postResults(Tstar_v,ipc,ip,el)
   case (PLASTICITY_J2_ID)
     constitutive_postResults(startPos:endPos) = plastic_j2_postResults(Tstar_v,ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID)
     constitutive_postResults(startPos:endPos) = &
       plastic_phenopowerlaw_postResults(Tstar_v,ipc,ip,el)
   case (PLASTICITY_PHENOPLUS_ID)
     constitutive_postResults(startPos:endPos) = &
       plastic_phenoplus_postResults(Tstar_v,ipc,ip,el)
   case (PLASTICITY_DISLOTWIN_ID)
     constitutive_postResults(startPos:endPos) = &
       plastic_dislotwin_postResults(Tstar_v,temperature(homog)%p(offset),ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID)
     constitutive_postResults(startPos:endPos) = &
       plastic_disloucla_postResults(Tstar_v,temperature(homog)%p(offset),ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID)
     constitutive_postResults(startPos:endPos) = &
       plastic_nonlocal_postResults (Tstar_v,FeArray,ip,el)
 end select

 do mySource = 1_pInt, phase_Nsources(phase)
   startPos = endPos + 1_pInt
   endPos = endPos + sourceState(material_phase(ipc,ip,el))%p(mySource)%sizePostResults
   select case (phase_source(mySource,material_phase(ipc,ip,el)))
     case (SOURCE_damage_isoBrittle_ID)
       constitutive_postResults(startPos:endPos) = source_damage_isoBrittle_postResults(ipc, ip, el)
     case (SOURCE_damage_isoDuctile_ID)
       constitutive_postResults(startPos:endPos) = source_damage_isoDuctile_postResults(ipc, ip, el)
     case (SOURCE_damage_anisoBrittle_ID)
       constitutive_postResults(startPos:endPos) = source_damage_anisoBrittle_postResults(ipc, ip, el)
     case (SOURCE_damage_anisoDuctile_ID)
       constitutive_postResults(startPos:endPos) = source_damage_anisoDuctile_postResults(ipc, ip, el)
   end select
 enddo

end function constitutive_postResults

end module constitutive
