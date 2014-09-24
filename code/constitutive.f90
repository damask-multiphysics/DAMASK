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
   constitutive_thermal_maxSizeDotState

 public :: & 
   constitutive_init, &
   constitutive_homogenizedC, &
   constitutive_microstructure, &
   constitutive_LpAndItsTangent, &
   constitutive_TandItsTangent, &
   constitutive_collectDotState, &
   constitutive_collectDeltaState, &
   constitutive_getLocalDamage, &
   constitutive_putLocalDamage, & 
   constitutive_getNonLocalDamage, &
   constitutive_getAdiabaticThermal, &
   constitutive_putAdiabaticThermal, &
   constitutive_getConductionThermal, &
   constitutive_postResults
 
 private :: &
   constitutive_hooke_TandItsTangent
 
contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates arrays pointing to array of the various constitutive modules
!--------------------------------------------------------------------------------------------------
subroutine constitutive_init
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
   phase_Noutput, &
   homogenization_Ngrains, &
   homogenization_maxNgrains, &
   ELASTICITY_HOOKE_ID, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_J2_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOKMC_ID, &
   PLASTICITY_TITANMOD_ID, &
   PLASTICITY_NONLOCAL_ID ,&
   ELASTICITY_HOOKE_label, &
   PLASTICITY_NONE_label, &
   PLASTICITY_J2_label, &
   PLASTICITY_PHENOPOWERLAW_label, &
   PLASTICITY_DISLOTWIN_label, &
   PLASTICITY_DISLOKMC_label, &
   PLASTICITY_TITANMOD_label, &
   PLASTICITY_NONLOCAL_label, &
   LOCAL_DAMAGE_NONE_ID, &
   LOCAL_DAMAGE_BRITTLE_ID, &
   LOCAL_THERMAL_none_ID, &
   LOCAL_THERMAL_HEATGEN_ID, &
   LOCAL_DAMAGE_NONE_label, &
   LOCAL_DAMAGE_BRITTLE_label, &
   LOCAL_THERMAL_none_label, &
   LOCAL_THERMAL_HEATGEN_label, &
   plasticState, &
   damageState, &
   thermalState, &
   mappingConstitutive
 

 use constitutive_none
 use constitutive_j2
 use constitutive_phenopowerlaw
 use constitutive_dislotwin
 use constitutive_dislokmc
 use constitutive_titanmod
 use constitutive_nonlocal
 use damage_none
 use damage_brittle
 use thermal_none
 use thermal_adiabatic
 implicit none
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt) :: &
  e, &                                                                                        !< maximum number of elements
  phase, &
  instance

 integer(pInt), dimension(:,:), pointer :: thisSize
 character(len=64), dimension(:,:), pointer :: thisOutput
 character(len=32) :: outputName                                                                    !< name of output, intermediate fix until HDF5 output is ready
 logical :: knownPlasticity, knownDamage, knownThermal, nonlocalConstitutionPresent
 nonlocalConstitutionPresent = .false.
 
!--------------------------------------------------------------------------------------------------
! parse plasticities from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(phase_plasticity == PLASTICITY_NONE_ID))          call constitutive_none_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_J2_ID))            call constitutive_j2_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_PHENOPOWERLAW_ID)) call constitutive_phenopowerlaw_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_DISLOTWIN_ID))     call constitutive_dislotwin_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_DISLOKMC_ID))      call constitutive_dislokmc_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_TITANMOD_ID))      call constitutive_titanmod_init(FILEUNIT)
 if (any(phase_plasticity == PLASTICITY_NONLOCAL_ID)) then
  call constitutive_nonlocal_init(FILEUNIT)
  call constitutive_nonlocal_stateInit()
 endif
 close(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! parse damage from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(phase_damage == LOCAL_DAMAGE_NONE_ID))       call damage_none_init(FILEUNIT)
 if (any(phase_damage == LOCAL_DAMAGE_BRITTLE_ID))    call damage_brittle_init(FILEUNIT)
 close(FILEUNIT)
 
!--------------------------------------------------------------------------------------------------
! parse thermal from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(phase_thermal == LOCAL_THERMAL_none_ID))       call thermal_none_init(FILEUNIT)
! if (any(phase_thermal == LOCAL_THERMAL_HEATGEN_ID)) call thermal_heatgen_init(FILEUNIT)
 close(FILEUNIT)

 write(6,'(/,a)')   ' <<<+-  constitutive init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
!--------------------------------------------------------------------------------------------------
! write description file for constitutive phase output
 call IO_write_jobFile(FILEUNIT,'outputConstitutive') 
 do phase = 1_pInt,material_Nphase
   instance = phase_plasticityInstance(phase)                                                       ! which instance of a plasticity is present phase
   knownPlasticity = .true.                                                                         ! assume valid
   select case(phase_plasticity(phase))                                                             ! split per constititution
     case (PLASTICITY_NONE_ID)
       outputName = PLASTICITY_NONE_label
       thisOutput => null()                                                                         ! constitutive_none_output
       thisSize   => null()                                                                         ! constitutive_none_sizePostResult
     case (PLASTICITY_J2_ID)
       outputName = PLASTICITY_J2_label
       thisOutput => constitutive_j2_output
       thisSize   => constitutive_j2_sizePostResult
     case (PLASTICITY_PHENOPOWERLAW_ID)
       outputName = PLASTICITY_PHENOPOWERLAW_label
       thisOutput => constitutive_phenopowerlaw_output
       thisSize   => constitutive_phenopowerlaw_sizePostResult
     case (PLASTICITY_DISLOTWIN_ID)
       outputName = PLASTICITY_DISLOTWIN_label
       thisOutput => constitutive_dislotwin_output
       thisSize   => constitutive_dislotwin_sizePostResult
     case (PLASTICITY_DISLOKMC_ID)
       outputName = PLASTICITY_DISLOKMC_label
       thisOutput => constitutive_dislokmc_output
       thisSize   => constitutive_dislokmc_sizePostResult
     case (PLASTICITY_TITANMOD_ID)
       outputName = PLASTICITY_TITANMOD_label
       thisOutput => constitutive_titanmod_output
       thisSize   => constitutive_titanmod_sizePostResult
     case (PLASTICITY_NONLOCAL_ID)
       outputName = PLASTICITY_NONLOCAL_label
       thisOutput => constitutive_nonlocal_output
       thisSize   => constitutive_nonlocal_sizePostResult
     case default
       knownPlasticity = .false.
   end select   
   write(FILEUNIT,'(/,a,/)') '['//trim(phase_name(phase))//']'
   if (knownPlasticity) then
     write(FILEUNIT,'(a)') '(plasticity)'//char(9)//trim(outputName)
     if (phase_plasticity(phase) /= PLASTICITY_NONE_ID) then
       do e = 1_pInt,phase_Noutput(phase)
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
       thisOutput => null()
       thisSize   => null()
     case (LOCAL_DAMAGE_BRITTLE_ID)
       outputName = LOCAL_DAMAGE_BRITTLE_label
       thisOutput => damage_brittle_output
       thisSize   => damage_brittle_sizePostResult
     case default
       knownDamage = .false.
   end select   
   if (knownDamage) then
     write(FILEUNIT,'(a)') '(damage)'//char(9)//trim(outputName)
     if (phase_damage(phase) /= LOCAL_DAMAGE_none_ID) then
       do e = 1_pInt,phase_Noutput(phase)
         write(FILEUNIT,'(a,i4)') trim(thisOutput(e,instance))//char(9),thisSize(e,instance)
       enddo
     endif
   endif
   instance = phase_thermalInstance(phase)                                                              ! which instance is present phase
   knownThermal = .true.
   select case(phase_thermal(phase))                                                                 ! split per constititution
     case (LOCAL_THERMAL_none_ID)
       outputName = LOCAL_THERMAL_NONE_label
       thisOutput => null()
       thisSize   => null()
     case (LOCAL_THERMAL_heatgen_ID)
       outputName = LOCAL_THERMAL_HEATGEN_label
       thisOutput => null()
       thisSize   => null()
     case default
       knownThermal = .false.
   end select   
   if (knownThermal) then
     write(FILEUNIT,'(a)') '(thermal)'//char(9)//trim(outputName)
     if (phase_thermal(phase) /= LOCAL_THERMAL_none_ID) then
       do e = 1_pInt,phase_Noutput(phase)
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
   plasticState,&
   mappingConstitutive

 use constitutive_titanmod, only: &
   constitutive_titanmod_homogenizedC
 use constitutive_dislotwin, only: &
   constitutive_dislotwin_homogenizedC
 use constitutive_dislokmc, only: &
   constitutive_dislokmc_homogenizedC
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
     constitutive_homogenizedC = constitutive_dislotwin_homogenizedC(ipc,ip,el)
   case (PLASTICITY_DISLOKMC_ID)
     constitutive_homogenizedC = constitutive_dislokmc_homogenizedC(ipc,ip,el) 
   case (PLASTICITY_TITANMOD_ID)
     constitutive_homogenizedC = constitutive_titanmod_homogenizedC (ipc,ip,el)
   case default
     constitutive_homogenizedC = lattice_C66(1:6,1:6,material_phase (ipc,ip,el))
     
 end select

end function constitutive_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different constitutive models
!--------------------------------------------------------------------------------------------------
subroutine constitutive_microstructure(temperature, Tstar_v, Fe, Fp, ipc, ip, el)
 use prec, only: &
   pReal 
 use material, only: &
   phase_plasticity, &
   phase_damage, &
   material_phase, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOKMC_ID, &
   PLASTICITY_TITANMOD_ID, &
   PLASTICITY_NONLOCAL_ID, &
   LOCAL_DAMAGE_BRITTLE_ID

 use constitutive_titanmod, only: &
   constitutive_titanmod_microstructure
 use constitutive_nonlocal, only: &
   constitutive_nonlocal_microstructure
 use constitutive_dislotwin, only: &
   constitutive_dislotwin_microstructure
 use constitutive_dislokmc, only: &
   constitutive_dislokmc_microstructure
 use damage_brittle, only: &
   damage_brittle_microstructure

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   temperature
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal),   intent(in), dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fp                                                                                               !< plastic deformation gradient
 real(pReal) :: damage, Tstar_v_effective(6)
 
 damage = constitutive_getNonlocalDamage(ipc,ip,el)
 Tstar_v_effective = Tstar_v/(damage*damage)

 select case (phase_plasticity(material_phase(ipc,ip,el)))
       
   case (PLASTICITY_DISLOTWIN_ID)
     call constitutive_dislotwin_microstructure(temperature,ipc,ip,el)
   case (PLASTICITY_DISLOKMC_ID)
     call constitutive_dislokmc_microstructure(temperature,ipc,ip,el)
   case (PLASTICITY_TITANMOD_ID)
     call constitutive_titanmod_microstructure (temperature,ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID)
     call constitutive_nonlocal_microstructure (Fe,Fp,          ip,el)

 end select
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_BRITTLE_ID)
     call damage_brittle_microstructure(Tstar_v_effective, Fe, ipc, ip, el)

 end select

end subroutine constitutive_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LpAndItsTangent(Lp, dLp_dTstar, Tstar_v, temperature, ipc, ip, el)
 use prec, only: &
   pReal 
 use math, only: &
   math_identity2nd
 use material, only: &
   phase_plasticity, &
   material_phase, &
   plasticState,&
   mappingConstitutive, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_J2_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOKMC_ID, &
   PLASTICITY_TITANMOD_ID, &
   PLASTICITY_NONLOCAL_ID
 use constitutive_j2, only: &
   constitutive_j2_LpAndItsTangent
 use constitutive_phenopowerlaw, only: &
   constitutive_phenopowerlaw_LpAndItsTangent
 use constitutive_dislotwin, only: &
   constitutive_dislotwin_LpAndItsTangent
 use constitutive_dislokmc, only: &
   constitutive_dislokmc_LpAndItsTangent
 use constitutive_titanmod, only: &
   constitutive_titanmod_LpAndItsTangent
 use constitutive_nonlocal, only: &
   constitutive_nonlocal_LpAndItsTangent
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   Temperature
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(out), dimension(3,3) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal),   intent(out), dimension(9,9) :: &
   dLp_dTstar                                                                                       !< derivative of Lp with respect to Tstar (4th-order tensor)
 real(pReal) :: damage, Tstar_v_effective(6)
 
 damage = constitutive_getNonlocalDamage(ipc,ip,el)
 Tstar_v_effective = Tstar_v/(damage*damage)
 select case (phase_plasticity(material_phase(ipc,ip,el)))
 
   case (PLASTICITY_NONE_ID)
     Lp = 0.0_pReal
     dLp_dTstar = 0.0_pReal
   case (PLASTICITY_J2_ID)
     call constitutive_j2_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v_effective,ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID)
     call constitutive_phenopowerlaw_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v_effective,ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID)
     call constitutive_nonlocal_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v_effective,temperature,ip,el)
   case (PLASTICITY_DISLOTWIN_ID)
     call constitutive_dislotwin_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v_effective,temperature,ipc,ip,el)
   case (PLASTICITY_DISLOKMC_ID)
     call constitutive_dislokmc_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v_effective,temperature,ipc,ip,el)
   case (PLASTICITY_TITANMOD_ID)
     call constitutive_titanmod_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v_effective,temperature,ipc,ip,el)

 end select
 
end subroutine constitutive_LpAndItsTangent



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
 real(pReal)   :: damage
 real(pReal), dimension(3,3,3,3) :: C

 damage = constitutive_getNonlocalDamage(ipc,ip,el)
 C = damage*damage*math_Mandel66to3333(constitutive_homogenizedC(ipc,ip,el))
 T = math_mul3333xx33(C,0.5_pReal*(math_mul33x33(math_transpose33(Fe),Fe)-math_I3) - &
                        lattice_thermalExpansion33(1:3,1:3,mappingConstitutive(2,ipc,ip,el))* &
                        (constitutive_getConductionThermal(ipc,ip,el) - &
                         lattice_referenceTemperature(mappingConstitutive(2,ipc,ip,el))))
 
 dT_dFe = 0.0_pReal
 forall (i=1_pInt:3_pInt, j=1_pInt:3_pInt, k=1_pInt:3_pInt, l=1_pInt:3_pInt) &
   dT_dFe(i,j,k,l) = sum(C(i,j,l,1:3)*Fe(k,1:3))                                                     ! dT*_ij/dFe_kl
 
end subroutine constitutive_hooke_TandItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure 
!--------------------------------------------------------------------------------------------------
subroutine constitutive_collectDotState(Tstar_v, FeArray, FpArray, Temperature, subdt, subfracArray,&
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
   phase_damage, &
   phase_thermal, &
   material_phase, &
   homogenization_maxNgrains, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_J2_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOKMC_ID, &
   PLASTICITY_TITANMOD_ID, &
   PLASTICITY_NONLOCAL_ID, &
   LOCAL_DAMAGE_BRITTLE_ID, &
   LOCAL_THERMAL_HEATGEN_ID
 use constitutive_j2, only:  &
   constitutive_j2_dotState
 use constitutive_phenopowerlaw, only: &
   constitutive_phenopowerlaw_dotState
 use constitutive_dislotwin, only: &
   constitutive_dislotwin_dotState
 use constitutive_dislokmc, only: &
   constitutive_dislokmc_dotState
 use constitutive_titanmod, only: &
   constitutive_titanmod_dotState
 use constitutive_nonlocal, only: &
   constitutive_nonlocal_dotState
 use damage_brittle, only: &
   damage_brittle_dotState
 use thermal_adiabatic, only: &
   thermal_adiabatic_dotState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  intent(in) :: &
   Temperature, &
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
 
 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) &
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
 
 select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_J2_ID)
     call constitutive_j2_dotState           (Tstar_v,ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID)
     call constitutive_phenopowerlaw_dotState(Tstar_v,ipc,ip,el)
   case (PLASTICITY_DISLOTWIN_ID)
     call constitutive_dislotwin_dotState    (Tstar_v,Temperature,ipc,ip,el)
   case (PLASTICITY_DISLOKMC_ID)
     call constitutive_dislokmc_dotState    (Tstar_v,Temperature,ipc,ip,el)
   case (PLASTICITY_TITANMOD_ID)
     call constitutive_titanmod_dotState     (Tstar_v,Temperature,ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID)
     call constitutive_nonlocal_dotState     (Tstar_v,FeArray,FpArray,Temperature, subdt, &
                                              subfracArray,ip,el)
 end select
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_BRITTLE_ID)
     call damage_brittle_dotState(ipc, ip, el)

 end select

 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_HEATGEN_ID)
!     call thermal_adiabatic_dotState(Tstar_v, Lp, ipc, ip, el)
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
 use constitutive_nonlocal, only: &
   constitutive_nonlocal_deltaState
 
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
     call constitutive_nonlocal_deltaState(Tstar_v,ip,el)
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
!> @brief Returns the local(unregularised)  damage 
!--------------------------------------------------------------------------------------------------
function constitutive_getLocalDamage(ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_DAMAGE_none_ID, &
   LOCAL_DAMAGE_BRITTLE_ID, &
   phase_damage
 use damage_brittle, only: &
   constitutive_brittle_getDamage

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: constitutive_getLocalDamage
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_none_ID)
     constitutive_getLocalDamage = 1.0_pReal
     
   case (LOCAL_DAMAGE_BRITTLE_ID)
     constitutive_getLocalDamage = constitutive_brittle_getDamage(ipc, ip, el)
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
   LOCAL_DAMAGE_BRITTLE_ID, &
   phase_damage
 use damage_brittle, only: &
   constitutive_brittle_putDamage

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   localDamage
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_BRITTLE_ID)
     call constitutive_brittle_putDamage(ipc, ip, el, localDamage)

 end select

end subroutine constitutive_putLocalDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns nonlocal (regularised) damage
!--------------------------------------------------------------------------------------------------
function constitutive_getNonlocalDamage(ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
    material_homog, &
    mappingHomogenization, &
    fieldDamage, &
    field_damage_type, &
    FIELD_DAMAGE_LOCAL_ID, &
    FIELD_DAMAGE_NONLOCAL_ID

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: constitutive_getNonlocalDamage
 
 select case(field_damage_type(material_homog(ip,el)))                                                   
 
   case (FIELD_DAMAGE_LOCAL_ID)
    constitutive_getNonlocalDamage = constitutive_getLocalDamage(ipc, ip, el)
    
   case (FIELD_DAMAGE_NONLOCAL_ID)
    constitutive_getNonlocalDamage =    fieldDamage(material_homog(ip,el))% &
      field(1,mappingHomogenization(1,ip,el))                                                     ! Taylor type 

 end select

end function constitutive_getNonlocalDamage
!--------------------------------------------------------------------------------------------------
!> @brief returns local (unregularised) temperature
!--------------------------------------------------------------------------------------------------
function constitutive_getAdiabaticThermal(ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_THERMAL_none_ID, &
   LOCAL_THERMAL_HEATGEN_ID, &
   phase_thermal
 use thermal_adiabatic, only: &
   constitutive_heatgen_getThermal
 use lattice, only: &
   lattice_referenceTemperature

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: constitutive_getAdiabaticThermal
 
 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_none_ID)
     constitutive_getAdiabaticThermal = lattice_referenceTemperature(material_phase(ipc,ip,el))
     
   case (LOCAL_THERMAL_HEATGEN_ID)
     constitutive_getAdiabaticThermal = constitutive_heatgen_getThermal(ipc, ip, el)
 end select

end function constitutive_getAdiabaticThermal

!--------------------------------------------------------------------------------------------------
!> @brief Returns the local(unregularised)  damage 
!--------------------------------------------------------------------------------------------------
subroutine constitutive_putAdiabaticThermal(ipc, ip, el, localTemperature)
 use prec, only: &
   pReal
 use material, only: &
   material_phase, &
   LOCAL_THERMAL_HEATGEN_ID, &
   phase_thermal
 use thermal_adiabatic, only: &
   constitutive_heatgen_putThermal

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   localTemperature
 
 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_HEATGEN_ID)
     call constitutive_heatgen_putThermal(ipc, ip, el, localTemperature)

 end select

end subroutine constitutive_putAdiabaticThermal

!--------------------------------------------------------------------------------------------------
!> @brief returns nonlocal (regularised) temperature
!--------------------------------------------------------------------------------------------------
function constitutive_getConductionThermal(ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   mappingHomogenization, &
   material_phase, &
   fieldThermal, &
   field_thermal_type, &
   FIELD_THERMAL_ADIABATIC_ID, &
   FIELD_THERMAL_CONDUCTION_ID, &
   material_homog
 use lattice, only: &
   lattice_referenceTemperature
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: constitutive_getConductionThermal
 
   select case(field_thermal_type(material_homog(ip,el)))                                                   
   
     case (FIELD_THERMAL_ADIABATIC_ID)
      constitutive_getConductionThermal = constitutive_getAdiabaticThermal(ipc, ip, el)      
      
     case (FIELD_THERMAL_CONDUCTION_ID)
      constitutive_getConductionThermal =    fieldThermal(material_homog(ip,el))% &
        field(1,mappingHomogenization(1,ip,el))                                                     ! Taylor type 

   end select

end function constitutive_getConductionThermal
!--------------------------------------------------------------------------------------------------
!> @brief returns array of constitutive results
!--------------------------------------------------------------------------------------------------
function constitutive_postResults(Tstar_v, FeArray, temperature, ipc, ip, el)
 use prec, only: &
   pReal 
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   plasticState, &
   damageState, &
   thermalState, &
   phase_plasticity, &
   phase_damage, &
   phase_thermal, &
   material_phase, &
   homogenization_maxNgrains, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_J2_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOKMC_ID, &
   PLASTICITY_TITANMOD_ID, &
   PLASTICITY_NONLOCAL_ID, &
   LOCAL_DAMAGE_BRITTLE_ID, &
   LOCAL_THERMAL_HEATGEN_ID
 use constitutive_j2, only: &
#ifdef HDF
   constitutive_j2_postResults2,&
#endif
   constitutive_j2_postResults
 use constitutive_phenopowerlaw, only: &
   constitutive_phenopowerlaw_postResults
 use constitutive_dislotwin, only: &
   constitutive_dislotwin_postResults
 use constitutive_dislokmc, only: &
   constitutive_dislokmc_postResults
 use constitutive_titanmod, only: &
   constitutive_titanmod_postResults
 use constitutive_nonlocal, only: &
   constitutive_nonlocal_postResults
#ifdef multiphysicsOut
 use damage_brittle, only: &
   damage_brittle_postResults
! use thermal_adiabatic, only: &
!   thermal_adiabatic_postResults
#endif

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
#ifdef multiphysicsOut
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%sizePostResults + &
                        damageState( material_phase(ipc,ip,el))%sizePostResults + &
                        thermalState(material_phase(ipc,ip,el))%sizePostResults) :: & 
   constitutive_postResults
#else
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%sizePostResults) :: &
   constitutive_postResults
#endif
 real(pReal),  intent(in) :: &
   temperature
 real(pReal),  intent(in), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   FeArray                                                                                          !< elastic deformation gradient
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal) :: damage, Tstar_v_effective(6)
 integer(pInt) :: startPos, endPos
 
 damage = constitutive_getNonlocalDamage(ipc,ip,el)
 Tstar_v_effective = damage*damage*Tstar_v

 constitutive_postResults = 0.0_pReal
 
 startPos = 1_pInt
 endPos = plasticState(material_phase(ipc,ip,el))%sizePostResults
 select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_TITANMOD_ID)
     constitutive_postResults(startPos:endPos) = constitutive_titanmod_postResults(ipc,ip,el)
   case (PLASTICITY_J2_ID)
     constitutive_postResults(startPos:endPos) = constitutive_j2_postResults(Tstar_v_effective,ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID)
     constitutive_postResults(startPos:endPos) = &
       constitutive_phenopowerlaw_postResults(Tstar_v_effective,ipc,ip,el)
   case (PLASTICITY_DISLOTWIN_ID)
     constitutive_postResults(startPos:endPos) = &
       constitutive_dislotwin_postResults(Tstar_v_effective,Temperature,ipc,ip,el)
   case (PLASTICITY_DISLOKMC_ID)
     constitutive_postResults(startPos:endPos) = &
       constitutive_dislokmc_postResults(Tstar_v_effective,Temperature,ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID)
     constitutive_postResults(startPos:endPos) = &
       constitutive_nonlocal_postResults (Tstar_v_effective,FeArray,ip,el)
 end select

#ifdef multiphysicsOut
 startPos = endPos + 1_pInt
 endPos = endPos + damageState(material_phase(ipc,ip,el))%sizePostResults
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_BRITTLE_ID)
     constitutive_postResults(startPos:endPos) = damage_brittle_postResults(ipc, ip, el)
 end select

 startPos = endPos + 1_pInt
 endPos = endPos + thermalState(material_phase(ipc,ip,el))%sizePostResults
 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (LOCAL_THERMAL_HEATGEN_ID)
!     constitutive_postResults(startPos:endPos) = thermal_adiabatic_postResults(ipc, ip, el)
 end select
#endif
  
end function constitutive_postResults


end module constitutive
