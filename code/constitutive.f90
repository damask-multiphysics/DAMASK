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
   constitutive_maxSizePostResults, &
   constitutive_maxSizeDotState, &
   constitutive_damage_maxSizePostResults, &
   constitutive_damage_maxSizeDotState, &
   constitutive_thermal_maxSizePostResults, &
   constitutive_thermal_maxSizeDotState, &
   constitutive_vacancy_maxSizePostResults, &
   constitutive_vacancy_maxSizeDotState

 public :: & 
   constitutive_init, &
   constitutive_homogenizedC, &
   constitutive_damagedC, &
   constitutive_microstructure, &
   constitutive_LpAndItsTangent, &
   constitutive_LiAndItsTangent, &
   constitutive_getFi, &
   constitutive_getUndamagedFi, &
   constitutive_putFi, &
   constitutive_getFi0, &
   constitutive_getPartionedFi0, &
   constitutive_TandItsTangent, &
   constitutive_collectDotState, &
   constitutive_collectDeltaState, &
   constitutive_getLocalDamage, &
   constitutive_putLocalDamage, & 
   constitutive_getDamage, &
   constitutive_getSlipDamage, &
   constitutive_getDamageDiffusion33, &
   constitutive_getAdiabaticTemperature, &
   constitutive_putAdiabaticTemperature, &
   constitutive_getTemperature, &
   constitutive_getHeatGeneration, &
   constitutive_getLocalVacancyConcentration, &
   constitutive_putLocalVacancyConcentration, &
   constitutive_getVacancyConcentration, &
   constitutive_getVacancyDiffusion33, &
   constitutive_getVacancyMobility33, &
   constitutive_getVacancyPotentialDrivingForce, &
   constitutive_postResults
 
 private :: &
   constitutive_hooke_TandItsTangent
 
contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates arrays pointing to array of the various constitutive modules
!--------------------------------------------------------------------------------------------------
subroutine constitutive_init(temperature_init)
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
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use numerics, only: &
   worldrank, &
   numerics_integrator
 use IO, only: &
   IO_error, &
   IO_open_file, &
   IO_open_jobFile_stat, &
   IO_write_jobFile, &
   IO_write_jobIntFile, &
   IO_timeStamp
 use mesh, only: &
   mesh_maxNips, &
   mesh_NcpElems, &
   mesh_element, &
   FE_Nips, &
   FE_geomtype
 use material, only: &
   material_phase, &
   material_Nphase, &
   material_localFileExt, &    
   material_configFile, &    
   phase_name, &
   phase_elasticity, &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_damage, &
   phase_damageInstance, &
   phase_thermal, &
   phase_thermalInstance, &
   phase_vacancy, &
   phase_vacancyInstance, &
   phase_Noutput, &
   homogenization_Ngrains, &
   homogenization_maxNgrains, &
   ELASTICITY_hooke_ID, &
   PLASTICITY_none_ID, &
   PLASTICITY_j2_ID, &
   PLASTICITY_phenopowerlaw_ID, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_dislokmc_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_titanmod_ID, &
   PLASTICITY_nonlocal_ID ,&
   ELASTICITY_HOOKE_label, &
   PLASTICITY_NONE_label, &
   PLASTICITY_J2_label, &
   PLASTICITY_PHENOPOWERLAW_label, &
   PLASTICITY_DISLOTWIN_label, &
   PLASTICITY_DISLOKMC_label, &
   PLASTICITY_DISLOUCLA_label, &
   PLASTICITY_TITANMOD_label, &
   PLASTICITY_NONLOCAL_label, &
   LOCAL_DAMAGE_none_ID, &
   LOCAL_DAMAGE_isoBrittle_ID, &
   LOCAL_DAMAGE_isoDuctile_ID, &
   LOCAL_DAMAGE_anisoBrittle_ID, &
   LOCAL_DAMAGE_anisoDuctile_ID, &
   LOCAL_DAMAGE_gurson_ID, &
   LOCAL_DAMAGE_phaseField_ID, &
   LOCAL_THERMAL_isothermal_ID, &
   LOCAL_THERMAL_adiabatic_ID, &
   LOCAL_VACANCY_constant_ID, &
   LOCAL_VACANCY_generation_ID, &
   LOCAL_DAMAGE_none_LABEL, &
   LOCAL_DAMAGE_isoBrittle_LABEL, &
   LOCAL_DAMAGE_isoDuctile_LABEL, &
   LOCAL_DAMAGE_anisoBrittle_LABEL, &
   LOCAL_DAMAGE_anisoDuctile_LABEL, &
   LOCAL_DAMAGE_gurson_LABEL, &
   LOCAL_DAMAGE_phaseField_label, &
   LOCAL_THERMAL_isothermal_label, &
   LOCAL_THERMAL_adiabatic_label, &
   LOCAL_VACANCY_constant_label, &
   LOCAL_VACANCY_generation_label, &
   plasticState, &
   damageState, &
   thermalState, &
   vacancyState, &
   mappingConstitutive
 

 use plastic_none
 use plastic_j2
 use plastic_phenopowerlaw
 use plastic_dislotwin
 use plastic_dislokmc
 use plastic_disloucla
 use plastic_titanmod
 use plastic_nonlocal
 use damage_none
 use damage_isoBrittle
 use damage_isoDuctile
 use damage_anisoDuctile
 use damage_anisoBrittle
 use damage_gurson
 use damage_phaseField
 use thermal_isothermal
 use thermal_adiabatic
 use vacancy_constant
 use vacancy_generation

 implicit none
 real(pReal), intent(in)  :: temperature_init                                                       !< initial temperature
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt) :: &
  e, &                                                                                              !< maximum number of elements
  phase, &
  instance

 integer(pInt), dimension(:,:), pointer :: thisSize
 integer(pInt), dimension(:)  , pointer :: thisNoutput
 character(len=64), dimension(:,:), pointer :: thisOutput
 character(len=32) :: outputName                                                                    !< name of output, intermediate fix until HDF5 output is ready
 logical :: knownPlasticity, knownDamage, knownThermal, knownVacancy, nonlocalConstitutionPresent
 nonlocalConstitutionPresent = .false.
 
!--------------------------------------------------------------------------------------------------
! parse plasticities from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(phase_plasticity == PLASTICITY_NONE_ID))          call plastic_none_init
 if (any(phase_plasticity == PLASTICITY_J2_ID))            call plastic_j2_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_PHENOPOWERLAW_ID)) call plastic_phenopowerlaw_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_DISLOTWIN_ID))     call plastic_dislotwin_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_DISLOKMC_ID))      call plastic_dislokmc_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_DISLOUCLA_ID))     call plastic_disloucla_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_TITANMOD_ID))      call plastic_titanmod_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_NONLOCAL_ID)) then
  call plastic_nonlocal_init(FILEUNIT)
  call plastic_nonlocal_stateInit()
 endif
 close(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! parse damage from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(phase_damage == LOCAL_DAMAGE_none_ID))             call damage_none_init
 if (any(phase_damage == LOCAL_DAMAGE_isoBrittle_ID))       call damage_isoBrittle_init(FILEUNIT)
 if (any(phase_damage == LOCAL_DAMAGE_isoductile_ID))       call damage_isoDuctile_init(FILEUNIT)
 if (any(phase_damage == LOCAL_DAMAGE_anisoBrittle_ID))     call damage_anisoBrittle_init(FILEUNIT)
 if (any(phase_damage == LOCAL_DAMAGE_anisoductile_ID))     call damage_anisoDuctile_init(FILEUNIT)
 if (any(phase_damage == LOCAL_DAMAGE_gurson_ID))           call damage_gurson_init(FILEUNIT)
 if (any(phase_damage == LOCAL_DAMAGE_phaseField_ID))       call damage_phaseField_init(FILEUNIT)
 close(FILEUNIT)
 
!--------------------------------------------------------------------------------------------------
! parse thermal from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(phase_thermal == LOCAL_THERMAL_isothermal_ID))    call thermal_isothermal_init(temperature_init)
 if (any(phase_thermal == LOCAL_THERMAL_adiabatic_ID))     call thermal_adiabatic_init(FILEUNIT,temperature_init)
 close(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! parse vacancy model from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(phase_vacancy == LOCAL_VACANCY_constant_ID))      call vacancy_constant_init
 if (any(phase_vacancy == LOCAL_VACANCY_generation_ID))    call vacancy_generation_init(FILEUNIT)
 close(FILEUNIT)

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  constitutive init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
!--------------------------------------------------------------------------------------------------
! write description file for constitutive phase output
 call IO_write_jobFile(FILEUNIT,'outputConstitutive') 
 do phase = 1_pInt,material_Nphase
   instance = phase_plasticityInstance(phase)                                                       ! which instance of a plasticity is present phase
   knownPlasticity = .true.                                                                         ! assume valid
   select case(phase_plasticity(phase))                                                             ! split per constititution
     case (PLASTICITY_NONE_ID)
       outputName = PLASTICITY_NONE_label
       thisNoutput => null()
       thisOutput => null()                                                                         ! plastic_none_output
       thisSize   => null()                                                                         ! plastic_none_sizePostResult
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
     case (PLASTICITY_DISLOTWIN_ID)
       outputName = PLASTICITY_DISLOTWIN_label
       thisNoutput => plastic_dislotwin_Noutput
       thisOutput => plastic_dislotwin_output
       thisSize   => plastic_dislotwin_sizePostResult
     case (PLASTICITY_DISLOKMC_ID)
       outputName = PLASTICITY_DISLOKMC_label
       thisNoutput => plastic_dislokmc_Noutput
       thisOutput => plastic_dislokmc_output
       thisSize   => plastic_dislokmc_sizePostResult
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
#ifdef multiphysicsOut
   instance = phase_damageInstance(phase)                                                           ! which instance of a plasticity is present phase
   knownDamage = .true.
   select case(phase_damage(phase))                                                                 ! split per constititution
     case (LOCAL_DAMAGE_none_ID)
       outputName = LOCAL_DAMAGE_NONE_label
       thisNoutput => null()
       thisOutput => null()
       thisSize   => null()
     case (LOCAL_DAMAGE_isoBrittle_ID)
       outputName = LOCAL_DAMAGE_isoBrittle_LABEL
       thisNoutput => damage_isoBrittle_Noutput
       thisOutput => damage_isoBrittle_output
       thisSize   => damage_isoBrittle_sizePostResult
     case (LOCAL_DAMAGE_isoDuctile_ID)
       outputName = LOCAL_DAMAGE_isoDuctile_LABEL
       thisNoutput => damage_isoDuctile_Noutput
       thisOutput => damage_isoDuctile_output
       thisSize   => damage_isoDuctile_sizePostResult
     case (LOCAL_DAMAGE_anisoBrittle_ID)
       outputName = LOCAL_DAMAGE_anisoBrittle_label
       thisNoutput => damage_anisoBrittle_Noutput
       thisOutput => damage_anisoBrittle_output
       thisSize   => damage_anisoBrittle_sizePostResult
     case (LOCAL_DAMAGE_anisoDuctile_ID)
       outputName = LOCAL_DAMAGE_anisoDuctile_LABEL
       thisNoutput => damage_anisoDuctile_Noutput
       thisOutput => damage_anisoDuctile_output
       thisSize   => damage_anisoDuctile_sizePostResult
     case (LOCAL_DAMAGE_gurson_ID)
       outputName = LOCAL_DAMAGE_gurson_label
       thisNoutput => damage_gurson_Noutput
       thisOutput => damage_gurson_output
       thisSize   => damage_gurson_sizePostResult
     case (LOCAL_DAMAGE_phaseField_ID)
       outputName = LOCAL_DAMAGE_phaseField_label
       thisNoutput => damage_phaseField_Noutput
       thisOutput => damage_phaseField_output
       thisSize   => damage_phaseField_sizePostResult
     case default
       knownDamage = .false.
   end select   
   if (knownDamage) then
     write(FILEUNIT,'(a)') '(damage)'//char(9)//trim(outputName)
     if (phase_damage(phase) /= LOCAL_DAMAGE_none_ID) then
       do e = 1_pInt,thisNoutput(instance)
         write(FILEUNIT,'(a,i4)') trim(thisOutput(e,instance))//char(9),thisSize(e,instance)
       enddo
     endif
   endif
   instance = phase_thermalInstance(phase)                                                              ! which instance is present phase
   knownThermal = .true.
   select case(phase_thermal(phase))                                                                 ! split per constititution
     case (LOCAL_THERMAL_isothermal_ID)
       outputName = LOCAL_THERMAL_ISOTHERMAL_label
       thisNoutput => null()
       thisOutput => null()
       thisSize   => null()
     case (LOCAL_THERMAL_adiabatic_ID)
       outputName = LOCAL_THERMAL_ADIABATIC_label
       thisNoutput => thermal_adiabatic_Noutput
       thisOutput => thermal_adiabatic_output
       thisSize   => thermal_adiabatic_sizePostResult
     case default
       knownThermal = .false.
   end select   
   if (knownThermal) then
     write(FILEUNIT,'(a)') '(thermal)'//char(9)//trim(outputName)
     if (phase_thermal(phase) /= LOCAL_THERMAL_isothermal_ID) then
       do e = 1_pInt,thisNoutput(instance)
         write(FILEUNIT,'(a,i4)') trim(thisOutput(e,instance))//char(9),thisSize(e,instance)
       enddo
     endif
   endif
   instance = phase_vacancyInstance(phase)                                                              ! which instance is present phase
   knownVacancy = .true.
   select case(phase_vacancy(phase))                                                                 ! split per constititution
     case (LOCAL_VACANCY_constant_ID)
       outputName = LOCAL_VACANCY_constant_label
       thisNoutput => null()
       thisOutput => null()
       thisSize   => null()
     case (LOCAL_VACANCY_generation_ID)
       outputName = LOCAL_VACANCY_generation_label
       thisNoutput => vacancy_generation_Noutput
       thisOutput => vacancy_generation_output
       thisSize   => vacancy_generation_sizePostResult
     case default
       knownVacancy = .false.
   end select   
   if (knownVacancy) then
     write(FILEUNIT,'(a)') '(vacancy)'//char(9)//trim(outputName)
     if (phase_vacancy(phase) /= LOCAL_VACANCY_constant_ID) then
       do e = 1_pInt,thisNoutput(instance)
         write(FILEUNIT,'(a,i4)') trim(thisOutput(e,instance))//char(9),thisSize(e,instance)
       enddo
     endif
   endif
#endif
 enddo
 close(FILEUNIT)
 
 constitutive_maxSizeDotState = 0_pInt
 constitutive_maxSizePostResults = 0_pInt
 constitutive_damage_maxSizePostResults = 0_pInt
 constitutive_damage_maxSizeDotState = 0_pInt
 constitutive_thermal_maxSizePostResults = 0_pInt
 constitutive_thermal_maxSizeDotState = 0_pInt
 constitutive_vacancy_maxSizePostResults = 0_pInt
 constitutive_vacancy_maxSizeDotState = 0_pInt

 PhaseLoop2:do phase = 1_pInt,material_Nphase
   plasticState(phase)%partionedState0 = plasticState(phase)%State0
   plasticState(phase)%State = plasticState(phase)%State0
   constitutive_maxSizeDotState = max(constitutive_maxSizeDotState, plasticState(phase)%sizeDotState)
   constitutive_maxSizePostResults = max(constitutive_maxSizePostResults, plasticState(phase)%sizePostResults)
   damageState(phase)%partionedState0 = damageState(phase)%State0
   damageState(phase)%State = damageState(phase)%State0
   constitutive_damage_maxSizeDotState = max(constitutive_damage_maxSizeDotState, damageState(phase)%sizeDotState)
   constitutive_damage_maxSizePostResults = max(constitutive_damage_maxSizePostResults, damageState(phase)%sizePostResults)
   thermalState(phase)%partionedState0 = thermalState(phase)%State0
   thermalState(phase)%State = thermalState(phase)%State0
   constitutive_thermal_maxSizeDotState = max(constitutive_thermal_maxSizeDotState, thermalState(phase)%sizeDotState)
   constitutive_thermal_maxSizePostResults = max(constitutive_thermal_maxSizePostResults, thermalState(phase)%sizePostResults)
   vacancyState(phase)%partionedState0 = vacancyState(phase)%State0
   vacancyState(phase)%State = vacancyState(phase)%State0
   constitutive_vacancy_maxSizeDotState = max(constitutive_vacancy_maxSizeDotState, vacancyState(phase)%sizeDotState)
   constitutive_vacancy_maxSizePostResults = max(constitutive_vacancy_maxSizePostResults, vacancyState(phase)%sizePostResults)
 enddo PhaseLoop2

#ifdef HDF
 call  HDF5_mappingConstitutive(mappingConstitutive)
 do phase = 1_pInt,material_Nphase
   instance = phase_plasticityInstance(phase)                                                       ! which instance of a plasticity is present phase
   select case(phase_plasticity(phase))                                                             ! split per constititution
     case (PLASTICITY_NONE_ID)
     case (PLASTICITY_J2_ID)
   end select
 enddo
#endif

#ifdef TODO
!--------------------------------------------------------------------------------------------------
! report
 constitutive_maxSizeState       = maxval(constitutive_sizeState)
 constitutive_maxSizeDotState    = maxval(constitutive_sizeDotState)
 
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
   write(6,'(a32,1x,7(i8,1x))')   'maxSizeDotState:    ', constitutive_maxSizeDotState
   write(6,'(a32,1x,7(i8,1x))')   'maxSizePostResults: ', constitutive_maxSizePostResults
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
   PLASTICITY_DISLOKMC_ID, &
   PLASTICITY_DISLOUCLA_ID, &
   plasticState,&
   mappingConstitutive

 use plastic_titanmod, only: &
   plastic_titanmod_homogenizedC
 use plastic_dislotwin, only: &
   plastic_dislotwin_homogenizedC
 use plastic_dislokmc, only: &
   plastic_dislokmc_homogenizedC
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
   case (PLASTICITY_DISLOKMC_ID)
     constitutive_homogenizedC = plastic_dislokmc_homogenizedC(ipc,ip,el) 
   case (PLASTICITY_DISLOUCLA_ID)
     constitutive_homogenizedC = plastic_disloucla_homogenizedC(ipc,ip,el)  
   case (PLASTICITY_TITANMOD_ID)
     constitutive_homogenizedC = plastic_titanmod_homogenizedC (ipc,ip,el)
   case default
     constitutive_homogenizedC = lattice_C66(1:6,1:6,material_phase (ipc,ip,el))
     
 end select

end function constitutive_homogenizedC

!--------------------------------------------------------------------------------------------------
!> @brief returns the damaged elasticity matrix if relevant 
!--------------------------------------------------------------------------------------------------
function constitutive_damagedC(ipc,ip,el)
 use prec, only: &
   pReal 
 use material, only: &
   material_phase, &
   LOCAL_DAMAGE_isoBrittle_ID, &
   LOCAL_DAMAGE_isoDuctile_ID, &
   LOCAL_DAMAGE_phaseField_ID, &   
   phase_damage
 use damage_isoBrittle, only: &
   damage_isoBrittle_getDamagedC66
 use damage_isoDuctile, only: &
   damage_isoDuctile_getDamagedC66  
 use damage_phaseField, only: &
   damage_phaseField_getDamagedC66  

 implicit none
 real(pReal), dimension(6,6) :: constitutive_damagedC
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el 
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_isoBrittle_ID)
     constitutive_damagedC = damage_isoBrittle_getDamagedC66(constitutive_homogenizedC(ipc,ip,el), &
                                                             ipc,ip,el)
                                                             
   case (LOCAL_DAMAGE_isoDuctile_ID)
     constitutive_damagedC = damage_isoDuctile_getDamagedC66(constitutive_homogenizedC(ipc,ip,el), &
                                                             ipc,ip,el)
                                                             
   case (LOCAL_DAMAGE_phaseField_ID)
     constitutive_damagedC = damage_phaseField_getDamagedC66(constitutive_homogenizedC(ipc,ip,el), &
                                                             ipc,ip,el)
   case default
     constitutive_damagedC = constitutive_homogenizedC(ipc,ip,el)

 end select 
   
end function constitutive_damagedC

!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different constitutive models
!--------------------------------------------------------------------------------------------------
subroutine constitutive_microstructure(Tstar_v, Fe, Fp, subdt, ipc, ip, el)
 use prec, only: &
   pReal 
 use material, only: &
   phase_plasticity, &
   phase_damage, &
   phase_vacancy, &
   material_phase, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_dislokmc_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_titanmod_ID, &
   PLASTICITY_nonlocal_ID, &
   LOCAL_DAMAGE_isoBrittle_ID, &
   LOCAL_DAMAGE_isoDuctile_ID, &
   LOCAL_DAMAGE_anisoBrittle_ID, &
   LOCAL_DAMAGE_anisoDuctile_ID, &
   LOCAL_DAMAGE_gurson_ID, &
   LOCAL_DAMAGE_phaseField_ID, &
   LOCAL_VACANCY_generation_ID
 use plastic_titanmod, only: &
   plastic_titanmod_microstructure
 use plastic_nonlocal, only: &
   plastic_nonlocal_microstructure
 use plastic_dislotwin, only: &
   plastic_dislotwin_microstructure
 use plastic_dislokmc, only: &
   plastic_dislokmc_microstructure
 use plastic_disloucla, only: &
   plastic_disloucla_microstructure
 use damage_isoBrittle, only: &
   damage_isoBrittle_microstructure
 use damage_isoDuctile, only: &
   damage_isoDuctile_microstructure
 use damage_anisoBrittle, only: &
   damage_anisoBrittle_microstructure
 use damage_anisoDuctile, only: &
   damage_anisoDuctile_microstructure
 use damage_gurson, only: &
   damage_gurson_microstructure
 use damage_phaseField, only: &
   damage_phaseField_microstructure
 use vacancy_generation, only: &
   vacancy_generation_microstructure    

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal),   intent(in), dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fp                                                                                               !< plastic deformation gradient
 real(pReal),   intent(in) :: &
   subdt                                                                                            !< timestep

 select case (phase_plasticity(material_phase(ipc,ip,el)))
       
   case (PLASTICITY_DISLOTWIN_ID)
     call plastic_dislotwin_microstructure(constitutive_getTemperature(ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_DISLOKMC_ID)
     call plastic_dislokmc_microstructure(constitutive_getTemperature(ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID)
     call plastic_disloucla_microstructure(constitutive_getTemperature(ipc,ip,el),ipc,ip,el)  
   case (PLASTICITY_TITANMOD_ID)
     call plastic_titanmod_microstructure (constitutive_getTemperature(ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID)
     call plastic_nonlocal_microstructure (Fe,Fp,          ip,el)

 end select
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_isoBrittle_ID)
     call damage_isoBrittle_microstructure(constitutive_homogenizedC(ipc,ip,el), Fe, subdt, &
                                           ipc, ip, el)
   case (LOCAL_DAMAGE_isoDuctile_ID)
     call damage_isoDuctile_microstructure(subdt, ipc, ip, el)
   case (LOCAL_DAMAGE_anisoBrittle_ID)
     call damage_anisoBrittle_microstructure(Tstar_v, subdt, ipc, ip, el)
   case (LOCAL_DAMAGE_anisoDuctile_ID)
     call damage_anisoDuctile_microstructure(subdt, ipc, ip, el)
   case (LOCAL_DAMAGE_gurson_ID)
     call damage_gurson_microstructure(ipc, ip, el)
   case (LOCAL_DAMAGE_phaseField_ID)
     call damage_phaseField_microstructure(constitutive_homogenizedC(ipc,ip,el), Fe, &
                                           constitutive_getVacancyConcentration(ipc, ip, el), &
                                           subdt, ipc, ip, el)

 end select

 select case (phase_vacancy(material_phase(ipc,ip,el)))
   case (LOCAL_VACANCY_generation_ID)
     call vacancy_generation_microstructure(constitutive_homogenizedC(ipc,ip,el), Fe, &
                                            constitutive_getTemperature(ipc,ip,el), &
                                            subdt,ipc,ip,el)

 end select

end subroutine constitutive_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LpAndItsTangent(Lp, dLp_dTstar, Tstar_v, ipc, ip, el)
 use prec, only: &
   pReal 
 use math, only: &
   math_identity2nd
 use material, only: &
   phase_plasticity, &
   phase_plasticityInstance, &
   material_phase, &
   plasticState,&
   mappingConstitutive, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_J2_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOKMC_ID, &
   PLASTICITY_DISLOUCLA_ID, &
   PLASTICITY_TITANMOD_ID, &
   PLASTICITY_NONLOCAL_ID
 use plastic_j2, only: &
   plastic_j2_LpAndItsTangent
 use plastic_phenopowerlaw, only: &
   plastic_phenopowerlaw_LpAndItsTangent, &
   plastic_phenopowerlaw_totalNslip
 use plastic_dislotwin, only: &
   plastic_dislotwin_LpAndItsTangent, &
   plastic_dislotwin_totalNslip
 use plastic_dislokmc, only: &
   plastic_dislokmc_LpAndItsTangent, &
   plastic_dislokmc_totalNslip
 use plastic_disloucla, only: &
   plastic_disloucla_LpAndItsTangent, &
   plastic_disloucla_totalNslip
 use plastic_titanmod, only: &
   plastic_titanmod_LpAndItsTangent, &
   plastic_titanmod_totalNslip
 use plastic_nonlocal, only: &
   plastic_nonlocal_LpAndItsTangent, &
   totalNslip
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(out), dimension(3,3) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal),   intent(out), dimension(9,9) :: &
   dLp_dTstar                                                                                       !< derivative of Lp with respect to Tstar (4th-order tensor)
 integer(pInt) :: &
   nSlip

 select case (phase_plasticity(material_phase(ipc,ip,el)))
 
   case (PLASTICITY_NONE_ID)
     Lp = 0.0_pReal
     dLp_dTstar = 0.0_pReal
   case (PLASTICITY_J2_ID)
     nSlip = 1_pInt
     call plastic_j2_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v, &
                                     nSlip,constitutive_getSlipDamage(nSlip,Tstar_v,ipc,ip,el), &
                                     ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID)
     nSlip = plastic_phenopowerlaw_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))
     call plastic_phenopowerlaw_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v, &
                                                nSlip,constitutive_getSlipDamage(nSlip,Tstar_v,ipc,ip,el), &
                                                ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID)
     nSlip = totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))
     call plastic_nonlocal_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v, &
                                           constitutive_getTemperature(ipc,ip,el), &
                                           nSlip,constitutive_getSlipDamage(nSlip,Tstar_v,ipc,ip,el), &
                                           ipc,ip,el)
   case (PLASTICITY_DISLOTWIN_ID)
     nSlip = plastic_dislotwin_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))
     call plastic_dislotwin_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v, &
                                            constitutive_getTemperature(ipc,ip,el), &
                                            nSlip,constitutive_getSlipDamage(nSlip,Tstar_v,ipc,ip,el), &
                                            ipc,ip,el)
   case (PLASTICITY_DISLOKMC_ID)
     nSlip = plastic_dislokmc_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))
     call plastic_dislokmc_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v, &
                                           constitutive_getTemperature(ipc,ip,el), &
                                           nSlip,constitutive_getSlipDamage(nSlip,Tstar_v,ipc,ip,el), &
                                           ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID)
     nSlip = plastic_disloucla_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))
     call plastic_disloucla_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v, &
                                           constitutive_getTemperature(ipc,ip,el), &
                                           nSlip,constitutive_getSlipDamage(nSlip,Tstar_v,ipc,ip,el), &
                                           ipc,ip,el)  
   case (PLASTICITY_TITANMOD_ID)
     nSlip = plastic_titanmod_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))
     call plastic_titanmod_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v, &
                                           constitutive_getTemperature(ipc,ip,el), &
                                           nSlip,constitutive_getSlipDamage(nSlip,Tstar_v,ipc,ip,el), &
                                           ipc,ip,el)

 end select
 
end subroutine constitutive_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LiAndItsTangent(Li, dLi_dTstar, Tstar_v, Lp, ipc, ip, el)
 use prec, only: &
   pReal 
 use material, only: &
   phase_damage, &
   phase_thermal, &
   material_phase, &
   LOCAL_DAMAGE_anisoBrittle_ID, &
   LOCAL_THERMAL_adiabatic_ID
 use damage_anisoBrittle, only: &
   damage_anisoBrittle_LdAndItsTangent
 use thermal_adiabatic, only: &
   thermal_adiabatic_LTAndItsTangent
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(in),  dimension(3,3) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   Li                                                                                               !< intermediate velocity gradient
 real(pReal),   intent(out), dimension(9,9) :: &
   dLi_dTstar                                                                                       !< derivative of Li with respect to Tstar (2nd-order tensor)
 real(pReal), dimension(3,3) :: &
   Li_temp                                                                                          !< intermediate velocity gradient
 real(pReal), dimension(9,9) :: &
   dLi_dTstar_temp                                                                                  !< derivative of Li with respect to Tstar (4th-order tensor)

 Li = 0.0_pReal
 dLi_dTstar = 0.0_pReal
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_anisoBrittle_ID)
     call damage_anisoBrittle_LdAndItsTangent(Li_temp, dLi_dTstar_temp, Tstar_v, ipc, ip, el)
     Li = Li + Li_temp
     dLi_dTstar = dLi_dTstar + dLi_dTstar_temp
 
 end select

 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_adiabatic_ID)
     call thermal_adiabatic_LTAndItsTangent(Li_temp, dLi_dTstar_temp, Tstar_v, Lp, ipc, ip, el)
     Li = Li + Li_temp
     dLi_dTstar = dLi_dTstar + dLi_dTstar_temp

 end select
 
end subroutine constitutive_LiAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the intermediate deformation gradient  
!--------------------------------------------------------------------------------------------------
pure function constitutive_getFi(ipc, ip, el)
 use prec, only: &
   pReal 
 use math, only: &
   math_I3, &
   math_mul33x33
 use material, only: &
   phase_damage, &
   phase_thermal, &
   material_phase, &
   LOCAL_DAMAGE_anisoBrittle_ID, &
   LOCAL_THERMAL_adiabatic_ID
 use damage_anisoBrittle, only: &
   damage_anisoBrittle_getFd
 use thermal_adiabatic, only: &
   thermal_adiabatic_getFT
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   dimension(3,3) :: &
   constitutive_getFi                                                                              !< intermediate deformation gradient

 constitutive_getFi = math_I3
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_anisoBrittle_ID)
     constitutive_getFi = math_mul33x33(constitutive_getFi,damage_anisoBrittle_getFd (ipc, ip, el))
 
 end select

 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_adiabatic_ID)
     constitutive_getFi = math_mul33x33(constitutive_getFi,thermal_adiabatic_getFT (ipc, ip, el))

 end select
 
end function constitutive_getFi

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the undamaged 
!>  intermediate deformation gradient  
!--------------------------------------------------------------------------------------------------
pure function constitutive_getUndamagedFi(ipc, ip, el)
 use prec, only: &
   pReal 
 use math, only: &
   math_I3, &
   math_mul33x33
 use material, only: &
   phase_thermal, &
   material_phase, &
   LOCAL_THERMAL_adiabatic_ID
 use thermal_adiabatic, only: &
   thermal_adiabatic_getFT
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   dimension(3,3) :: &
   constitutive_getUndamagedFi                                                                              !< intermediate deformation gradient

 constitutive_getUndamagedFi = math_I3

 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_adiabatic_ID)
     constitutive_getUndamagedFi = math_mul33x33(constitutive_getUndamagedFi,thermal_adiabatic_getFT (ipc, ip, el))

 end select
 
end function constitutive_getUndamagedFi


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the intermediate deformation gradient  
!--------------------------------------------------------------------------------------------------
subroutine constitutive_putFi(Tstar_v, Lp, dt, ipc, ip, el)
 use prec, only: &
   pReal 
 use material, only: &
   phase_damage, &
   phase_thermal, &
   material_phase, &
   LOCAL_DAMAGE_anisoBrittle_ID, &
   LOCAL_THERMAL_adiabatic_ID
 use damage_anisoBrittle, only: &
   damage_anisoBrittle_putFd
 use thermal_adiabatic, only: &
   thermal_adiabatic_putFT
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(in),  dimension(3,3) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal),   intent(in) :: &
   dt
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_anisoBrittle_ID)
     call damage_anisoBrittle_putFd (Tstar_v, dt, ipc, ip, el)
 
 end select

 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_adiabatic_ID)
     call thermal_adiabatic_putFT (Tstar_v, Lp, dt, ipc, ip, el)

 end select
 
end subroutine constitutive_putFi


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the intermediate deformation gradient  
!--------------------------------------------------------------------------------------------------
pure function constitutive_getFi0(ipc, ip, el)
 use prec, only: &
   pReal 
 use math, only: &
   math_I3, &
   math_mul33x33
 use material, only: &
   phase_damage, &
   phase_thermal, &
   material_phase, &
   LOCAL_DAMAGE_anisoBrittle_ID, &
   LOCAL_THERMAL_adiabatic_ID
 use damage_anisoBrittle, only: &
   damage_anisoBrittle_getFd0
 use thermal_adiabatic, only: &
   thermal_adiabatic_getFT0
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   dimension(3,3) :: &
   constitutive_getFi0                                                                              !< intermediate deformation gradient

 constitutive_getFi0 = math_I3
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_anisoBrittle_ID)
     constitutive_getFi0 = math_mul33x33(constitutive_getFi0,damage_anisoBrittle_getFd0 (ipc, ip, el))
 
 end select

 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_adiabatic_ID)
     constitutive_getFi0 = math_mul33x33(constitutive_getFi0,thermal_adiabatic_getFT0 (ipc, ip, el))

 end select
 
end function constitutive_getFi0


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the intermediate deformation gradient  
!--------------------------------------------------------------------------------------------------
pure function constitutive_getPartionedFi0(ipc, ip, el)
 use prec, only: &
   pReal 
 use math, only: &
   math_I3, &
   math_mul33x33
 use material, only: &
   phase_damage, &
   phase_thermal, &
   material_phase, &
   LOCAL_DAMAGE_anisoBrittle_ID, &
   LOCAL_THERMAL_adiabatic_ID
 use damage_anisoBrittle, only: &
   damage_anisoBrittle_getPartionedFd0
 use thermal_adiabatic, only: &
   thermal_adiabatic_getPartionedFT0
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   dimension(3,3) :: &
   constitutive_getPartionedFi0                                                                     !< intermediate deformation gradient

 constitutive_getPartionedFi0 = math_I3
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_anisoBrittle_ID)
     constitutive_getPartionedFi0 = math_mul33x33(constitutive_getPartionedFi0, &
                                                  damage_anisoBrittle_getPartionedFd0(ipc, ip, el))
 
 end select

 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_adiabatic_ID)
     constitutive_getPartionedFi0 = math_mul33x33(constitutive_getPartionedFi0, &
                                                  thermal_adiabatic_getPartionedFT0(ipc, ip, el))

 end select
 
end function constitutive_getPartionedFi0


!--------------------------------------------------------------------------------------------------
!> @brief returns the 2nd Piola-Kirchhoff stress tensor and its tangent with respect to 
!> the elastic deformation gradient depending on the selected elastic law (so far no case switch
!! because only hooke is implemented
!--------------------------------------------------------------------------------------------------
subroutine constitutive_TandItsTangent(T, dT_dFe, Fe, ipc, ip, el)
 use prec, only: &
   pReal

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fe                                                                                               !< elastic deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   T                                                                                                !< 2nd Piola-Kirchhoff stress tensor
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dT_dFe                                                                                           !< derivative of 2nd P-K stress with respect to elastic deformation gradient
 
 call constitutive_hooke_TandItsTangent(T, dT_dFe, Fe, ipc, ip, el)

 
end subroutine constitutive_TandItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief returns the 2nd Piola-Kirchhoff stress tensor and its tangent with respect to 
!> the elastic deformation gradient using hookes law
!--------------------------------------------------------------------------------------------------
subroutine constitutive_hooke_TandItsTangent(T, dT_dFe, Fe, ipc, ip, el)
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
   mappingConstitutive
 use lattice, only: &
   lattice_referenceTemperature, &
   lattice_thermalExpansion33

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fe                                                                                               !< elastic deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   T                                                                                                !< 2nd Piola-Kirchhoff stress tensor
 real(pReal),   intent(out), dimension(3,3,3,3) :: & 
   dT_dFe                                                                                           !< dT/dFe
 
 integer(pInt) :: i, j, k, l
 real(pReal), dimension(3,3,3,3) :: C

 C = math_Mandel66to3333(constitutive_damagedC(ipc,ip,el))
 T = math_mul3333xx33(C,0.5_pReal*(math_mul33x33(math_transpose33(Fe),Fe)-math_I3))
 
 dT_dFe = 0.0_pReal
 forall (i=1_pInt:3_pInt, j=1_pInt:3_pInt, k=1_pInt:3_pInt, l=1_pInt:3_pInt) &
   dT_dFe(i,j,k,l) = sum(C(i,j,l,1:3)*Fe(k,1:3))                                                     ! dT*_ij/dFe_kl
 
end subroutine constitutive_hooke_TandItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure 
!--------------------------------------------------------------------------------------------------
subroutine constitutive_collectDotState(Tstar_v, Lp, FeArray, FpArray, subdt, subfracArray,&
                                                                                        ipc, ip, el)
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
   phase_plasticityInstance, &
   phase_damage, &
   phase_vacancy, &
   material_phase, &
   homogenization_maxNgrains, &
   PLASTICITY_none_ID, &
   PLASTICITY_j2_ID, &
   PLASTICITY_phenopowerlaw_ID, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_dislokmc_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_titanmod_ID, &
   PLASTICITY_nonlocal_ID, &
   LOCAL_DAMAGE_gurson_ID
 use plastic_j2, only:  &
   plastic_j2_dotState
 use plastic_phenopowerlaw, only: &
   plastic_phenopowerlaw_dotState, &
   plastic_phenopowerlaw_totalNslip
 use plastic_dislotwin, only: &
   plastic_dislotwin_dotState, &
   plastic_dislotwin_totalNslip
 use plastic_dislokmc, only: &
   plastic_dislokmc_dotState, &
   plastic_dislokmc_totalNslip
 use plastic_disloucla, only: &
   plastic_disloucla_dotState, &
   plastic_disloucla_totalNslip
 use plastic_titanmod, only: &
   plastic_titanmod_dotState, &
   plastic_titanmod_totalNslip
 use plastic_nonlocal, only: &
   plastic_nonlocal_dotState, &
   totalNslip
 use damage_gurson, only: &
   damage_gurson_dotState

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
 real(pReal),  intent(in), dimension(3,3) :: &
   Lp                                                                                               !< plastic velocity gradient
 integer(pLongInt) :: &
   tick, tock, & 
   tickrate, &
   maxticks
 integer(pInt) :: &
   nSlip
 
 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) &
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
 
 select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_J2_ID)
     nSlip = 1_pInt
     call plastic_j2_dotState           (Tstar_v,nSlip,constitutive_getSlipDamage(nSlip,Tstar_v,ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID)
     nSlip = plastic_phenopowerlaw_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))
     call plastic_phenopowerlaw_dotState(Tstar_v,nSlip,constitutive_getSlipDamage(nSlip,Tstar_v,ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_DISLOTWIN_ID)
     nSlip = plastic_dislotwin_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))
     call plastic_dislotwin_dotState    (Tstar_v,constitutive_getTemperature(ipc,ip,el), &
                                         nSlip,constitutive_getSlipDamage(nSlip,Tstar_v,ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_DISLOKMC_ID)
     nSlip = plastic_dislokmc_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))
     call plastic_dislokmc_dotState     (Tstar_v,constitutive_getTemperature(ipc,ip,el), &
                                         nSlip,constitutive_getSlipDamage(nSlip,Tstar_v,ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID)
     nSlip = plastic_disloucla_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))
     call plastic_disloucla_dotState     (Tstar_v,constitutive_getTemperature(ipc,ip,el), &
                                         nSlip,constitutive_getSlipDamage(nSlip,Tstar_v,ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_TITANMOD_ID)
     call plastic_titanmod_dotState     (Tstar_v,constitutive_getTemperature(ipc,ip,el), &
                                         ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID)
     nSlip = totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))
     call plastic_nonlocal_dotState     (Tstar_v,FeArray,FpArray,constitutive_getTemperature(ipc,ip,el), &
                                         nSlip,constitutive_getSlipDamage(nSlip,Tstar_v,ipc,ip,el),subdt,subfracArray,ip,el)
 end select
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_gurson_ID)
     call damage_gurson_dotState(Tstar_v, Lp, ipc, ip, el)
 end select

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
logical function constitutive_collectDeltaState(Tstar_v, ipc, ip, el)
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
   material_phase, &
   plasticState, &
   mappingConstitutive, &
   PLASTICITY_NONLOCAL_ID 
 use plastic_nonlocal, only: &
   plastic_nonlocal_deltaState
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress     
 integer(pLongInt) :: &
   tick, tock, & 
   tickrate, &
   maxticks

 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) &
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)

 select case (phase_plasticity(material_phase(ipc,ip,el)))

   case (PLASTICITY_NONLOCAL_ID)
     constitutive_collectDeltaState = .true.
     call plastic_nonlocal_deltaState(Tstar_v,ip,el)
   case default
     constitutive_collectDeltaState = .false.

 end select

 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) then
   call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
   !$OMP CRITICAL (debugTimingDeltaState)
     debug_cumDeltaStateCalls = debug_cumDeltaStateCalls + 1_pInt
     debug_cumDeltaStateTicks = debug_cumDeltaStateTicks + tock-tick
     !$OMP FLUSH (debug_cumDeltaStateTicks)
     if (tock < tick) debug_cumDeltaStateTicks  = debug_cumDeltaStateTicks + maxticks
   !$OMP END CRITICAL (debugTimingDeltaState)
 endif

end function constitutive_collectDeltaState


!--------------------------------------------------------------------------------------------------
!> @brief Returns the local(regularised)  damage 
!--------------------------------------------------------------------------------------------------
function constitutive_getLocalDamage(ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_DAMAGE_none_ID, &
   LOCAL_DAMAGE_isoBrittle_ID, &
   LOCAL_DAMAGE_isoDuctile_ID, &
   LOCAL_DAMAGE_anisoBrittle_ID, &
   LOCAL_DAMAGE_anisoDuctile_ID, &
   LOCAL_DAMAGE_gurson_ID, &
   LOCAL_DAMAGE_phaseField_ID, &
   phase_damage
 use damage_isoBrittle, only: &
   damage_isoBrittle_getLocalDamage
 use damage_isoDuctile, only: &
   damage_isoDuctile_getLocalDamage
 use damage_anisoBrittle, only: &
   damage_anisoBrittle_getLocalDamage
 use damage_anisoDuctile, only: &
   damage_anisoDuctile_getLocalDamage
 use damage_gurson, only: &
   damage_gurson_getLocalDamage
 use damage_phaseField, only: &
   damage_phaseField_getLocalDamage

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: constitutive_getLocalDamage
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_none_ID)
     constitutive_getLocalDamage = 1.0_pReal
     
   case (LOCAL_DAMAGE_isoBrittle_ID)
     constitutive_getLocalDamage = damage_isoBrittle_getLocalDamage(ipc, ip, el)
   
   case (LOCAL_DAMAGE_isoDuctile_ID)
     constitutive_getLocalDamage = damage_isoDuctile_getLocalDamage(ipc, ip, el)
     
   case (LOCAL_DAMAGE_anisoBrittle_ID)
     constitutive_getLocalDamage = damage_anisoBrittle_getLocalDamage(ipc, ip, el)

   case (LOCAL_DAMAGE_anisoDuctile_ID)     
     constitutive_getLocalDamage = damage_anisoDuctile_getLocalDamage(ipc, ip, el)
     
   case (LOCAL_DAMAGE_gurson_ID)
     constitutive_getLocalDamage = damage_gurson_getLocalDamage(ipc, ip, el)

   case (LOCAL_DAMAGE_phaseField_ID)
     constitutive_getLocalDamage = damage_phaseField_getLocalDamage(ipc, ip, el)

 end select

end function constitutive_getLocalDamage

!--------------------------------------------------------------------------------------------------
!> @brief Returns the local(unregularised)  damage 
!--------------------------------------------------------------------------------------------------
subroutine constitutive_putLocalDamage(ipc, ip, el, localDamage)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_DAMAGE_isoBrittle_ID, &
   LOCAL_DAMAGE_isoDuctile_ID, &
   LOCAL_DAMAGE_anisoBrittle_ID, &
   LOCAL_DAMAGE_anisoDuctile_ID, &
   LOCAL_DAMAGE_gurson_ID, &
   LOCAL_DAMAGE_phaseField_ID, &
   phase_damage
 use damage_isoBrittle, only: &
   damage_isoBrittle_putLocalDamage
 use damage_isoDuctile, only: &
   damage_isoDuctile_putLocalDamage
 use damage_anisoBrittle, only: &
   damage_anisoBrittle_putLocalDamage
 use damage_anisoDuctile, only: &
   damage_anisoDuctile_putLocalDamage
 use damage_gurson, only: &
   damage_gurson_putLocalDamage
 use damage_phaseField, only: &
   damage_phaseField_putLocalDamage

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   localDamage
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_isoBrittle_ID)
     call damage_isoBrittle_putLocalDamage(ipc, ip, el, localDamage)
   
   case (LOCAL_DAMAGE_isoDuctile_ID)
     call damage_isoDuctile_putLocalDamage(ipc, ip, el, localDamage)
     
   case (LOCAL_DAMAGE_anisoBrittle_ID)
     call damage_anisoBrittle_putLocalDamage(ipc, ip, el, localDamage)

   case (LOCAL_DAMAGE_anisoDuctile_ID)
     call damage_anisoDuctile_putLocalDamage(ipc, ip, el, localDamage)
     
   case (LOCAL_DAMAGE_gurson_ID)
     call damage_gurson_putLocalDamage(ipc, ip, el, localDamage)
 
   case (LOCAL_DAMAGE_phaseField_ID)
     call damage_phaseField_putLocalDamage(ipc, ip, el, localDamage)

 end select

end subroutine constitutive_putLocalDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns nonlocal (regularised) damage
!--------------------------------------------------------------------------------------------------
function constitutive_getDamage(ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_DAMAGE_none_ID, &
   LOCAL_DAMAGE_isoBrittle_ID, &
   LOCAL_DAMAGE_isoDuctile_ID, &
   LOCAL_DAMAGE_anisoBrittle_ID, &
   LOCAL_DAMAGE_anisoDuctile_ID, &
   LOCAL_DAMAGE_gurson_ID, &
   LOCAL_DAMAGE_phaseField_ID, &
   phase_damage
 use damage_isoBrittle, only: &
   damage_isoBrittle_getDamage
 use damage_isoDuctile, only: &
   damage_isoDuctile_getDamage
 use damage_anisoBrittle, only: &
   damage_anisoBrittle_getDamage
 use damage_anisoDuctile, only: &
   damage_anisoDuctile_getDamage
 use damage_gurson, only: &
   damage_gurson_getDamage
 use damage_phaseField, only: &
   damage_phaseField_getDamage

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: constitutive_getDamage
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_none_ID)
     constitutive_getDamage = 1.0_pReal
     
   case (LOCAL_DAMAGE_isoBrittle_ID)
     constitutive_getDamage = damage_isoBrittle_getDamage(ipc, ip, el)
   
   case (LOCAL_DAMAGE_isoDuctile_ID)
     constitutive_getDamage = damage_isoDuctile_getDamage(ipc, ip, el)
     
   case (LOCAL_DAMAGE_anisoBrittle_ID)
     constitutive_getDamage = damage_anisoBrittle_getDamage(ipc, ip, el)

   case (LOCAL_DAMAGE_anisoDuctile_ID)
     constitutive_getDamage = damage_anisoDuctile_getDamage(ipc, ip, el)
     
   case (LOCAL_DAMAGE_gurson_ID)
     constitutive_getDamage = damage_gurson_getDamage(ipc, ip, el)

   case (LOCAL_DAMAGE_phaseField_ID)
     constitutive_getDamage = damage_phaseField_getDamage(ipc, ip, el)

 end select

end function constitutive_getDamage

!--------------------------------------------------------------------------------------------------
!> @brief Returns the damage on each slip system
!--------------------------------------------------------------------------------------------------
function constitutive_getSlipDamage(nSlip, Tstar_v, ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_DAMAGE_isoDuctile_ID, &
   LOCAL_DAMAGE_anisoDuctile_ID, &
   LOCAL_DAMAGE_gurson_ID, &
   phase_damage
 use damage_isoDuctile, only: &
   damage_isoDuctile_getSlipDamage
 use damage_anisoDuctile, only: &
   damage_anisoDuctile_getSlipDamage
 use damage_gurson, only: &
   damage_gurson_getSlipDamage

 implicit none
 integer(pInt), intent(in) :: &
   nSlip, &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal) :: &
   constitutive_getSlipDamage(nSlip)
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_isoDuctile_ID)
     constitutive_getSlipDamage = damage_isoDuctile_getSlipDamage(ipc, ip, el)
     
   case (LOCAL_DAMAGE_anisoDuctile_ID)
     constitutive_getSlipDamage = damage_anisoDuctile_getSlipDamage(ipc, ip, el)
     
   case (LOCAL_DAMAGE_gurson_ID)
     constitutive_getSlipDamage = damage_gurson_getSlipDamage(Tstar_v, ipc, ip, el)

   case default
     constitutive_getSlipDamage = 1.0_pReal
 end select

end function constitutive_getSlipDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns damage diffusion tensor
!--------------------------------------------------------------------------------------------------
function constitutive_getDamageDiffusion33(ipc, ip, el)
 use prec, only: &
   pReal
 use lattice, only: &
   lattice_DamageDiffusion33
 use material, only: &
   material_phase, &
   phase_damage, &
   LOCAL_DAMAGE_isoBrittle_ID, &
   LOCAL_DAMAGE_phaseField_ID
 use damage_isoBrittle, only: &
   damage_isoBrittle_getDamageDiffusion33
 use damage_phaseField, only: &
   damage_phaseField_getDamageDiffusion33

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   constitutive_getDamageDiffusion33
 
 constitutive_getDamageDiffusion33 = lattice_DamageDiffusion33(1:3,1:3,material_phase(ipc,ip,el))
 select case(phase_damage(material_phase(ipc,ip,el)))                                                   
   case (LOCAL_DAMAGE_isoBrittle_ID)
    constitutive_getDamageDiffusion33 = damage_isoBrittle_getDamageDiffusion33(ipc, ip, el)
   case (LOCAL_DAMAGE_phaseField_ID)
    constitutive_getDamageDiffusion33 = damage_phaseField_getDamageDiffusion33(ipc, ip, el)
    
 end select

end function constitutive_getDamageDiffusion33

!--------------------------------------------------------------------------------------------------
!> @brief returns local (unregularised) temperature
!--------------------------------------------------------------------------------------------------
function constitutive_getAdiabaticTemperature(ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_THERMAL_isothermal_ID, &
   LOCAL_THERMAL_adiabatic_ID, &
   phase_thermal, &
   phase_thermalInstance
 use thermal_isothermal, only: &
   thermal_isothermal_temperature
 use thermal_adiabatic, only: &
   thermal_adiabatic_getTemperature

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: constitutive_getAdiabaticTemperature
 
 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_isothermal_ID)
     constitutive_getAdiabaticTemperature = &
       thermal_isothermal_temperature(phase_thermalInstance(material_phase(ipc,ip,el)))
     
   case (LOCAL_THERMAL_adiabatic_ID)
     constitutive_getAdiabaticTemperature = thermal_adiabatic_getTemperature(ipc, ip, el)
 end select

end function constitutive_getAdiabaticTemperature

!--------------------------------------------------------------------------------------------------
!> @brief assigns the local/nonlocal value of temperature to local thermal state
!--------------------------------------------------------------------------------------------------
subroutine constitutive_putAdiabaticTemperature(ipc, ip, el, localTemperature)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_THERMAL_adiabatic_ID, &
   phase_thermal
 use thermal_adiabatic, only: &
   thermal_adiabatic_putTemperature

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   localTemperature
 
 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_adiabatic_ID)
     call thermal_adiabatic_putTemperature(ipc, ip, el, localTemperature)

 end select

end subroutine constitutive_putAdiabaticTemperature

!--------------------------------------------------------------------------------------------------
!> @brief returns nonlocal (regularised) temperature
!--------------------------------------------------------------------------------------------------
function constitutive_getTemperature(ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   mappingHomogenization, &
   material_phase, &
   fieldThermal, &
   field_thermal_type, &
   FIELD_THERMAL_local_ID, &
   FIELD_THERMAL_nonlocal_ID, &
   material_homog
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: constitutive_getTemperature
 
 select case(field_thermal_type(material_homog(ip,el)))                                                   
   case (FIELD_THERMAL_local_ID)
    constitutive_getTemperature = constitutive_getAdiabaticTemperature(ipc, ip, el)      
    
   case (FIELD_THERMAL_nonlocal_ID)
    constitutive_getTemperature = fieldThermal(material_homog(ip,el))% &
      field(1,mappingHomogenization(1,ip,el))                                                     ! Taylor type 

 end select

end function constitutive_getTemperature

!--------------------------------------------------------------------------------------------------
!> @brief returns heat generation rate
!--------------------------------------------------------------------------------------------------
function constitutive_getHeatGeneration(Tstar_v, Lp, ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_THERMAL_isothermal_ID, &
   LOCAL_THERMAL_adiabatic_ID, &
   phase_thermal
 use thermal_adiabatic, only: &
   thermal_adiabatic_getHeatGeneration

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(in),  dimension(3,3) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal) :: constitutive_getHeatGeneration
 
 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_isothermal_ID)
     constitutive_getHeatGeneration = 0.0_pReal
     
   case (LOCAL_THERMAL_adiabatic_ID)
     constitutive_getHeatGeneration = thermal_adiabatic_getHeatGeneration(Tstar_v, Lp)
 end select

end function constitutive_getHeatGeneration

!--------------------------------------------------------------------------------------------------
!> @brief returns local vacancy concentration
!--------------------------------------------------------------------------------------------------
function constitutive_getLocalVacancyConcentration(ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_VACANCY_constant_ID, &
   LOCAL_VACANCY_generation_ID, &
   phase_vacancy
 use vacancy_generation, only: &
   vacancy_generation_getLocalConcentration
 use lattice, only: &
   lattice_equilibriumVacancyConcentration

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: constitutive_getLocalVacancyConcentration
 
 select case (phase_vacancy(material_phase(ipc,ip,el)))
   case (LOCAL_VACANCY_constant_ID)
     constitutive_getLocalVacancyConcentration = &
       lattice_equilibriumVacancyConcentration(material_phase(ipc,ip,el))
     
   case (LOCAL_VACANCY_generation_ID)
     constitutive_getLocalVacancyConcentration = vacancy_generation_getLocalConcentration(ipc, ip, el)
 end select

end function constitutive_getLocalVacancyConcentration

!--------------------------------------------------------------------------------------------------
!> @brief Puts local vacancy concentration 
!--------------------------------------------------------------------------------------------------
subroutine constitutive_putLocalVacancyConcentration(ipc, ip, el, localVacancyConcentration)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_VACANCY_generation_ID, &
   phase_vacancy
 use vacancy_generation, only: &
   vacancy_generation_putLocalConcentration

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   localVacancyConcentration
 
 select case (phase_vacancy(material_phase(ipc,ip,el)))
   case (LOCAL_VACANCY_generation_ID)
     call vacancy_generation_putLocalConcentration(ipc, ip, el, localVacancyConcentration)

 end select

end subroutine constitutive_putLocalVacancyConcentration

!--------------------------------------------------------------------------------------------------
!> @brief returns nonlocal vacancy concentration
!--------------------------------------------------------------------------------------------------
function constitutive_getVacancyConcentration(ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_VACANCY_constant_ID, &
   LOCAL_VACANCY_generation_ID, &
   phase_vacancy
 use vacancy_generation, only: &
   vacancy_generation_getConcentration
 use lattice, only: &
   lattice_equilibriumVacancyConcentration
 implicit none

 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: constitutive_getVacancyConcentration
 
 select case (phase_vacancy(material_phase(ipc,ip,el)))
   case (LOCAL_VACANCY_constant_ID)
     constitutive_getVacancyConcentration = &
       lattice_equilibriumVacancyConcentration(material_phase(ipc,ip,el))
     
   case (LOCAL_VACANCY_generation_ID)
     constitutive_getVacancyConcentration = vacancy_generation_getConcentration(ipc, ip, el)
 end select

end function constitutive_getVacancyConcentration

!--------------------------------------------------------------------------------------------------
!> @brief returns vacancy diffusion tensor
!--------------------------------------------------------------------------------------------------
function constitutive_getVacancyDiffusion33(ipc, ip, el)
 use prec, only: &
   pReal
 use lattice, only: &
   lattice_VacancyDiffusion33
 use material, only: &
   material_phase, &
   LOCAL_VACANCY_generation_ID, &
   phase_vacancy
 use vacancy_generation, only: &
   vacancy_generation_getVacancyDiffusion33

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   constitutive_getVacancyDiffusion33
 
 select case(phase_vacancy(material_phase(ipc,ip,el)))                                                   
   case (LOCAL_VACANCY_generation_ID)
    constitutive_getVacancyDiffusion33 = &
      vacancy_generation_getVacancyDiffusion33(ipc,ip,el)
    
 end select

end function constitutive_getVacancyDiffusion33

!--------------------------------------------------------------------------------------------------
!> @brief returns vacancy diffusion tensor
!--------------------------------------------------------------------------------------------------
function constitutive_getVacancyMobility33(ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_VACANCY_generation_ID, &
   phase_vacancy
 use vacancy_generation, only: &
   vacancy_generation_getVacancyMobility33

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   constitutive_getVacancyMobility33

 integer(pInt) :: &
   nSlip
 
 select case(phase_vacancy(material_phase(ipc,ip,el)))                                                   
   case (LOCAL_VACANCY_generation_ID)
    constitutive_getVacancyMobility33 = &
      vacancy_generation_getVacancyMobility33(nSlip,constitutive_getTemperature(ipc,ip,el), &
                                               ipc,ip,el)
    
 end select

end function constitutive_getVacancyMobility33

!--------------------------------------------------------------------------------------------------
!> @brief returns vacancy chemical potential driving force
!--------------------------------------------------------------------------------------------------
real(pReal) function constitutive_getVacancyPotentialDrivingForce(ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_VACANCY_generation_ID, &
   phase_vacancy
 use vacancy_generation, only: &
   vacancy_generation_getVacancyPotentialDrivingForce

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 
 select case(phase_vacancy(material_phase(ipc,ip,el)))                                                   
   case (LOCAL_VACANCY_generation_ID)
    constitutive_getVacancyPotentialDrivingForce = &
      vacancy_generation_getVacancyPotentialDrivingForce(constitutive_getDamage(ipc, ip, el), &
                                                         ipc,ip,el)
    
 end select

end function constitutive_getVacancyPotentialDrivingForce


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
   damageState, &
   thermalState, &
   vacancyState, &
   phase_plasticity, &
   phase_damage, &
   phase_thermal, &
   phase_vacancy, &
   material_phase, &
   homogenization_maxNgrains, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_J2_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOKMC_ID, &
   PLASTICITY_DISLOUCLA_ID, &
   PLASTICITY_TITANMOD_ID, &
   PLASTICITY_NONLOCAL_ID, &
   LOCAL_DAMAGE_isoBrittle_ID, &
   LOCAL_DAMAGE_isoDuctile_ID, &
   LOCAL_DAMAGE_anisoBrittle_ID, &
   LOCAL_DAMAGE_anisoDuctile_ID, &
   LOCAL_DAMAGE_gurson_ID, &
   LOCAL_DAMAGE_phaseField_ID, &
   LOCAL_THERMAL_ADIABATIC_ID, &
   LOCAL_VACANCY_generation_ID
 use plastic_j2, only: &
#ifdef HDF
   plastic_j2_postResults2,&
#endif
   plastic_j2_postResults
 use plastic_phenopowerlaw, only: &
   plastic_phenopowerlaw_postResults
 use plastic_dislotwin, only: &
   plastic_dislotwin_postResults
 use plastic_dislokmc, only: &
   plastic_dislokmc_postResults
 use plastic_disloucla, only: &
   plastic_disloucla_postResults
 use plastic_titanmod, only: &
   plastic_titanmod_postResults
 use plastic_nonlocal, only: &
   plastic_nonlocal_postResults
#ifdef multiphysicsOut
 use damage_isoBrittle, only: &
   damage_isoBrittle_postResults
 use damage_isoDuctile, only: &
   damage_isoDuctile_postResults
 use damage_anisoBrittle, only: &
   damage_anisoBrittle_postResults
 use damage_anisoDuctile, only: &
   damage_anisoDuctile_postResults
 use damage_gurson, only: &
   damage_gurson_postResults
 use damage_phaseField, only: &
   damage_phaseField_postResults
 use thermal_adiabatic, only: &
   thermal_adiabatic_postResults
 use vacancy_generation, only: &
   vacancy_generation_postResults
#endif

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
#ifdef multiphysicsOut
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%sizePostResults + &
                        damageState( material_phase(ipc,ip,el))%sizePostResults + &
                        thermalState(material_phase(ipc,ip,el))%sizePostResults + &
                        vacancyState(material_phase(ipc,ip,el))%sizePostResults) :: & 
   constitutive_postResults
#else
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%sizePostResults) :: &
   constitutive_postResults
#endif
 real(pReal),  intent(in), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   FeArray                                                                                          !< elastic deformation gradient
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 integer(pInt) :: &
   startPos, endPos
 
 constitutive_postResults = 0.0_pReal
 
 startPos = 1_pInt
 endPos = plasticState(material_phase(ipc,ip,el))%sizePostResults
 select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_TITANMOD_ID)
     constitutive_postResults(startPos:endPos) = plastic_titanmod_postResults(ipc,ip,el)
   case (PLASTICITY_J2_ID)
     constitutive_postResults(startPos:endPos) = plastic_j2_postResults(Tstar_v,ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID)
     constitutive_postResults(startPos:endPos) = &
       plastic_phenopowerlaw_postResults(Tstar_v,ipc,ip,el)
   case (PLASTICITY_DISLOTWIN_ID)
     constitutive_postResults(startPos:endPos) = &
       plastic_dislotwin_postResults(Tstar_v,constitutive_getTemperature(ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_DISLOKMC_ID)
     constitutive_postResults(startPos:endPos) = &
       plastic_dislokmc_postResults(Tstar_v,constitutive_getTemperature(ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID)
     constitutive_postResults(startPos:endPos) = &
       plastic_disloucla_postResults(Tstar_v,constitutive_getTemperature(ipc,ip,el),ipc,ip,el)  
   case (PLASTICITY_NONLOCAL_ID)
     constitutive_postResults(startPos:endPos) = &
       plastic_nonlocal_postResults (Tstar_v,FeArray,ip,el)
 end select

#ifdef multiphysicsOut
 startPos = endPos + 1_pInt
 endPos = endPos + damageState(material_phase(ipc,ip,el))%sizePostResults
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_isoBrittle_ID)
     constitutive_postResults(startPos:endPos) = damage_isoBrittle_postResults(ipc, ip, el)
   case (LOCAL_DAMAGE_isoDuctile_ID)
     constitutive_postResults(startPos:endPos) = damage_isoDuctile_postResults(ipc, ip, el)
   case (LOCAL_DAMAGE_anisoBrittle_ID)
     constitutive_postResults(startPos:endPos) = damage_anisoBrittle_postResults(ipc, ip, el)
   case (LOCAL_DAMAGE_anisoDuctile_ID)
     constitutive_postResults(startPos:endPos) = damage_anisoDuctile_postResults(ipc, ip, el)
   case (LOCAL_DAMAGE_gurson_ID)
     constitutive_postResults(startPos:endPos) = damage_gurson_postResults(ipc, ip, el)
   case (LOCAL_DAMAGE_phaseField_ID)
     constitutive_postResults(startPos:endPos) = damage_phaseField_postResults(ipc, ip, el)
 end select

 startPos = endPos + 1_pInt
 endPos = endPos + thermalState(material_phase(ipc,ip,el))%sizePostResults
 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_ADIABATIC_ID)
     constitutive_postResults(startPos:endPos) = thermal_adiabatic_postResults(ipc, ip, el)
 end select

 startPos = endPos + 1_pInt
 endPos = endPos + vacancyState(material_phase(ipc,ip,el))%sizePostResults
 select case (phase_vacancy(material_phase(ipc,ip,el)))
   case (LOCAL_VACANCY_generation_ID)
     constitutive_postResults(startPos:endPos) = vacancy_generation_postResults(ipc, ip, el)
 end select
#endif
  
end function constitutive_postResults


end module constitutive
