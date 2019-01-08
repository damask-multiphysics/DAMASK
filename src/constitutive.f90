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
   constitutive_LpAndItsTangents, &
   constitutive_LiAndItsTangents, &
   constitutive_initialFi, &
   constitutive_SandItsTangents, &
   constitutive_collectDotState, &
   constitutive_collectDeltaState, &
   constitutive_postResults, &
   constitutive_results

 private :: &
   constitutive_hooke_SandItsTangents

contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates arrays pointing to array of the various constitutive modules
!--------------------------------------------------------------------------------------------------
subroutine constitutive_init()
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
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
   IO_checkAndRewind, &
   IO_open_jobFile_stat, &
   IO_write_jobFile, &
   IO_write_jobIntFile, &
   IO_timeStamp
 use config, only: &
   config_phase
 use mesh, only: &
   FE_geomtype
 use config, only: &
   material_Nphase, &
   material_localFileExt, &
   phase_name, &
   material_configFile, &
   config_deallocate
 use material, only: &
   material_phase, &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Nsources, &
   phase_source, &
   phase_kinematics, &
   ELASTICITY_hooke_ID, &
   PLASTICITY_none_ID, &
   PLASTICITY_isotropic_ID, &
   PLASTICITY_phenopowerlaw_ID, &
   PLASTICITY_kinehardening_ID, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_nonlocal_ID ,&
   SOURCE_thermal_dissipation_ID, &
   SOURCE_thermal_externalheat_ID, &
   SOURCE_damage_isoBrittle_ID, &
   SOURCE_damage_isoDuctile_ID, &
   SOURCE_damage_anisoBrittle_ID, &
   SOURCE_damage_anisoDuctile_ID, &
   KINEMATICS_cleavage_opening_ID, &
   KINEMATICS_slipplane_opening_ID, &
   KINEMATICS_thermal_expansion_ID, &
   ELASTICITY_HOOKE_label, &
   PLASTICITY_NONE_label, &
   PLASTICITY_ISOTROPIC_label, &
   PLASTICITY_PHENOPOWERLAW_label, &
   PLASTICITY_KINEHARDENING_label, &
   PLASTICITY_DISLOTWIN_label, &
   PLASTICITY_DISLOUCLA_label, &
   PLASTICITY_NONLOCAL_label, &
   SOURCE_thermal_dissipation_label, &
   SOURCE_thermal_externalheat_label, &
   SOURCE_damage_isoBrittle_label, &
   SOURCE_damage_isoDuctile_label, &
   SOURCE_damage_anisoBrittle_label, &
   SOURCE_damage_anisoDuctile_label, &
   plasticState, &
   sourceState

 use plastic_none
 use plastic_isotropic
 use plastic_phenopowerlaw
 use plastic_kinehardening
 use plastic_dislotwin
 use plastic_disloucla
 use plastic_nonlocal
 use source_thermal_dissipation
 use source_thermal_externalheat
 use source_damage_isoBrittle
 use source_damage_isoDuctile
 use source_damage_anisoBrittle
 use source_damage_anisoDuctile
 use kinematics_cleavage_opening
 use kinematics_slipplane_opening
 use kinematics_thermal_expansion

 implicit none
 integer(pInt), parameter :: FILEUNIT = 204_pInt
 integer(pInt) :: &
   o, &                                                                                             !< counter in output loop
   ph, &                                                                                            !< counter in phase loop
   s, &                                                                                             !< counter in source loop
   ins                                                                                              !< instance of plasticity/source

 integer(pInt), dimension(:,:), pointer :: thisSize
 character(len=64), dimension(:,:), pointer :: thisOutput
 character(len=32) :: outputName                                                                    !< name of output, intermediate fix until HDF5 output is ready
 logical :: knownPlasticity, knownSource, nonlocalConstitutionPresent
 nonlocalConstitutionPresent = .false.

!--------------------------------------------------------------------------------------------------
! open material.config
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file

!--------------------------------------------------------------------------------------------------
! parse plasticities from config file
 if (any(phase_plasticity == PLASTICITY_NONE_ID))          call plastic_none_init
 if (any(phase_plasticity == PLASTICITY_ISOTROPIC_ID))     call plastic_isotropic_init
 if (any(phase_plasticity == PLASTICITY_PHENOPOWERLAW_ID)) call plastic_phenopowerlaw_init
 if (any(phase_plasticity == PLASTICITY_KINEHARDENING_ID)) call plastic_kinehardening_init
 if (any(phase_plasticity == PLASTICITY_DISLOTWIN_ID))     call plastic_dislotwin_init
 if (any(phase_plasticity == PLASTICITY_DISLOUCLA_ID))     call plastic_disloucla_init
 if (any(phase_plasticity == PLASTICITY_NONLOCAL_ID)) then
  call plastic_nonlocal_init(FILEUNIT)
  call plastic_nonlocal_stateInit()
 endif

!--------------------------------------------------------------------------------------------------
! parse source mechanisms from config file
 call IO_checkAndRewind(FILEUNIT)
 if (any(phase_source == SOURCE_thermal_dissipation_ID))     call source_thermal_dissipation_init(FILEUNIT)
 if (any(phase_source == SOURCE_thermal_externalheat_ID))    call source_thermal_externalheat_init(FILEUNIT)
 if (any(phase_source == SOURCE_damage_isoBrittle_ID))       call source_damage_isoBrittle_init(FILEUNIT)
 if (any(phase_source == SOURCE_damage_isoDuctile_ID))       call source_damage_isoDuctile_init(FILEUNIT)
 if (any(phase_source == SOURCE_damage_anisoBrittle_ID))     call source_damage_anisoBrittle_init(FILEUNIT)
 if (any(phase_source == SOURCE_damage_anisoDuctile_ID))     call source_damage_anisoDuctile_init(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! parse kinematic mechanisms from config file
 call IO_checkAndRewind(FILEUNIT)
 if (any(phase_kinematics == KINEMATICS_cleavage_opening_ID))  call kinematics_cleavage_opening_init(FILEUNIT)
 if (any(phase_kinematics == KINEMATICS_slipplane_opening_ID)) call kinematics_slipplane_opening_init(FILEUNIT)
 if (any(phase_kinematics == KINEMATICS_thermal_expansion_ID)) call kinematics_thermal_expansion_init(FILEUNIT)
 close(FILEUNIT)

 call config_deallocate('material.config/phase')

 write(6,'(/,a)')   ' <<<+-  constitutive init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 mainProcess: if (worldrank == 0) then
!--------------------------------------------------------------------------------------------------
! write description file for constitutive output
   call IO_write_jobFile(FILEUNIT,'outputConstitutive')
   PhaseLoop: do ph = 1_pInt,material_Nphase
     activePhase: if (any(material_phase == ph)) then
       ins = phase_plasticityInstance(ph)
       knownPlasticity = .true.                                                                     ! assume valid
       plasticityType: select case(phase_plasticity(ph))
         case (PLASTICITY_NONE_ID) plasticityType
           outputName = PLASTICITY_NONE_label
           thisOutput => null()
           thisSize   => null()
         case (PLASTICITY_ISOTROPIC_ID) plasticityType
           outputName = PLASTICITY_ISOTROPIC_label
           thisOutput => plastic_isotropic_output
           thisSize   => plastic_isotropic_sizePostResult
         case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
           outputName = PLASTICITY_PHENOPOWERLAW_label
           thisOutput => plastic_phenopowerlaw_output
           thisSize   => plastic_phenopowerlaw_sizePostResult
         case (PLASTICITY_KINEHARDENING_ID) plasticityType
           outputName = PLASTICITY_KINEHARDENING_label
           thisOutput => plastic_kinehardening_output
           thisSize   => plastic_kinehardening_sizePostResult  
         case (PLASTICITY_DISLOTWIN_ID) plasticityType
           outputName = PLASTICITY_DISLOTWIN_label
           thisOutput => plastic_dislotwin_output
           thisSize   => plastic_dislotwin_sizePostResult
         case (PLASTICITY_DISLOUCLA_ID) plasticityType
           outputName = PLASTICITY_DISLOUCLA_label
           thisOutput => plastic_disloucla_output
           thisSize   => plastic_disloucla_sizePostResult
         case (PLASTICITY_NONLOCAL_ID) plasticityType
           outputName = PLASTICITY_NONLOCAL_label
           thisOutput => plastic_nonlocal_output
           thisSize   => plastic_nonlocal_sizePostResult
         case default plasticityType
           knownPlasticity = .false.
       end select plasticityType
       write(FILEUNIT,'(/,a,/)') '['//trim(phase_name(ph))//']'
       if (knownPlasticity) then
         write(FILEUNIT,'(a)') '(plasticity)'//char(9)//trim(outputName)
         if (phase_plasticity(ph) /= PLASTICITY_NONE_ID) then
           OutputPlasticityLoop: do o = 1_pInt,size(thisOutput(:,ins))
             if(len(trim(thisOutput(o,ins))) > 0_pInt) &
               write(FILEUNIT,'(a,i4)') trim(thisOutput(o,ins))//char(9),thisSize(o,ins)
           enddo OutputPlasticityLoop
         endif
       endif
       
       SourceLoop: do s = 1_pInt, phase_Nsources(ph)
         knownSource = .true.                                                                       ! assume valid
         sourceType: select case (phase_source(s,ph))
           case (SOURCE_thermal_dissipation_ID) sourceType
             ins = source_thermal_dissipation_instance(ph)
             outputName = SOURCE_thermal_dissipation_label
             thisOutput => source_thermal_dissipation_output
             thisSize   => source_thermal_dissipation_sizePostResult
           case (SOURCE_thermal_externalheat_ID) sourceType
             ins = source_thermal_externalheat_instance(ph)
             outputName = SOURCE_thermal_externalheat_label
             thisOutput => source_thermal_externalheat_output
             thisSize   => source_thermal_externalheat_sizePostResult
           case (SOURCE_damage_isoBrittle_ID) sourceType
             ins = source_damage_isoBrittle_instance(ph)
             outputName = SOURCE_damage_isoBrittle_label
             thisOutput => source_damage_isoBrittle_output
             thisSize   => source_damage_isoBrittle_sizePostResult
           case (SOURCE_damage_isoDuctile_ID) sourceType
             ins = source_damage_isoDuctile_instance(ph)
             outputName = SOURCE_damage_isoDuctile_label
             thisOutput => source_damage_isoDuctile_output
             thisSize   => source_damage_isoDuctile_sizePostResult
           case (SOURCE_damage_anisoBrittle_ID) sourceType
             ins = source_damage_anisoBrittle_instance(ph)
             outputName = SOURCE_damage_anisoBrittle_label
             thisOutput => source_damage_anisoBrittle_output
             thisSize   => source_damage_anisoBrittle_sizePostResult
           case (SOURCE_damage_anisoDuctile_ID) sourceType
             ins = source_damage_anisoDuctile_instance(ph)
             outputName = SOURCE_damage_anisoDuctile_label
             thisOutput => source_damage_anisoDuctile_output
             thisSize   => source_damage_anisoDuctile_sizePostResult
           case default sourceType
             knownSource = .false.
         end select sourceType
         if (knownSource) then
           write(FILEUNIT,'(a)') '(source)'//char(9)//trim(outputName)
           OutputSourceLoop: do o = 1_pInt,size(thisOutput(:,ins))
             if(len(trim(thisOutput(o,ins))) > 0_pInt) &
               write(FILEUNIT,'(a,i4)') trim(thisOutput(o,ins))//char(9),thisSize(o,ins)
           enddo OutputSourceLoop
         endif
       enddo SourceLoop
     endif activePhase
   enddo PhaseLoop
   close(FILEUNIT)
 endif mainProcess

 constitutive_plasticity_maxSizeDotState = 0_pInt
 constitutive_plasticity_maxSizePostResults = 0_pInt
 constitutive_source_maxSizeDotState = 0_pInt
 constitutive_source_maxSizePostResults = 0_pInt

 PhaseLoop2:do ph = 1_pInt,material_Nphase
!--------------------------------------------------------------------------------------------------
! partition and inititalize state
   plasticState(ph)%partionedState0 = plasticState(ph)%state0
   plasticState(ph)%state           = plasticState(ph)%partionedState0
   forall(s = 1_pInt:phase_Nsources(ph))
     sourceState(ph)%p(s)%partionedState0 = sourceState(ph)%p(s)%state0
     sourceState(ph)%p(s)%state           = sourceState(ph)%p(s)%partionedState0
   end forall
!--------------------------------------------------------------------------------------------------
! determine max size of state and output
   constitutive_plasticity_maxSizeDotState    = max(constitutive_plasticity_maxSizeDotState,    &
                                                    plasticState(ph)%sizeDotState)
   constitutive_plasticity_maxSizePostResults = max(constitutive_plasticity_maxSizePostResults, &
                                                    plasticState(ph)%sizePostResults)
   constitutive_source_maxSizeDotState        = max(constitutive_source_maxSizeDotState, &
                                                    maxval(sourceState(ph)%p(:)%sizeDotState))
   constitutive_source_maxSizePostResults     = max(constitutive_source_maxSizePostResults, &
                                                    maxval(sourceState(ph)%p(:)%sizePostResults))
 enddo PhaseLoop2


end subroutine constitutive_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenize elasticity matrix
!> ToDo: homogenizedC66 would be more consistent
!--------------------------------------------------------------------------------------------------
function constitutive_homogenizedC(ipc,ip,el)
 use prec, only: &
   pReal
 use material, only: &
   phase_plasticity, &
   material_phase, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOUCLA_ID
 use plastic_dislotwin, only: &
   plastic_dislotwin_homogenizedC
 use lattice, only: &
   lattice_C66

 implicit none
 real(pReal), dimension(6,6) :: constitutive_homogenizedC
 integer(pInt), intent(in) :: &
   ipc, &                                                                                            !< component-ID of integration point
   ip, &                                                                                             !< integration point
   el                                                                                                !< element

 plasticityType: select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     constitutive_homogenizedC = plastic_dislotwin_homogenizedC(ipc,ip,el)
   case default plasticityType
     constitutive_homogenizedC = lattice_C66(1:6,1:6,material_phase (ipc,ip,el))
 end select plasticityType

end function constitutive_homogenizedC

!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different constitutive models
!--------------------------------------------------------------------------------------------------
subroutine constitutive_microstructure(orientations, Fe, Fp, ipc, ip, el)
 use prec, only: &
   pReal
 use material, only: &
   phasememberAt, &
   phase_plasticity, &
   phase_plasticityInstance, &
   material_phase, &
   material_homogenizationAt, &
   temperature, &
   thermalMapping, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_nonlocal_ID
 use plastic_nonlocal, only: &
   plastic_nonlocal_microstructure
 use plastic_dislotwin, only: &
   plastic_dislotwin_microstructure
 use plastic_disloUCLA, only: &
   plastic_disloUCLA_dependentState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in), dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fp                                                                                               !< plastic deformation gradient
 integer(pInt) :: &
   ho, &                                                                                            !< homogenization
   tme, &                                                                                           !< thermal member position
   instance, of
 real(pReal),   intent(in), dimension(:,:,:,:) :: &
   orientations                                                                                     !< crystal orientations as quaternions

 ho = material_homogenizationAt(el)
 tme = thermalMapping(ho)%p(ip,el)

 plasticityType: select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     call plastic_dislotwin_microstructure(temperature(ho)%p(tme),ipc,ip,el)
   case (PLASTICITY_DISLOUCLA_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     call plastic_disloUCLA_dependentState(instance,of)
   case (PLASTICITY_NONLOCAL_ID) plasticityType
     call plastic_nonlocal_microstructure (Fe,Fp,ip,el)
 end select plasticityType

end subroutine constitutive_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LpAndItsTangents(Lp, dLp_dS, dLp_dFi, S6, Fi, ipc, ip, el)
 use prec, only: &
   pReal
 use math, only: &
   math_mul33x33, &
   math_Mandel6to33, &
   math_Mandel33to6, &
   math_Plain99to3333
 use material, only: &
   phasememberAt, &
   phase_plasticity, &
   phase_plasticityInstance, &
   material_phase, &
   material_homogenizationAt, &
   temperature, &
   thermalMapping, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_ISOTROPIC_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_KINEHARDENING_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOUCLA_ID, &
   PLASTICITY_NONLOCAL_ID
 use plastic_isotropic, only: &
   plastic_isotropic_LpAndItsTangent
 use plastic_phenopowerlaw, only: &
   plastic_phenopowerlaw_LpAndItsTangent
 use plastic_kinehardening, only: &
   plastic_kinehardening_LpAndItsTangent  
 use plastic_dislotwin, only: &
   plastic_dislotwin_LpAndItsTangent
 use plastic_disloucla, only: &
   plastic_disloucla_LpAndItsTangent
 use plastic_nonlocal, only: &
   plastic_nonlocal_LpAndItsTangent

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in),  dimension(6) :: &
   S6                                                                                               !< 2nd Piola-Kirchhoff stress (vector notation)
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dLp_dS, &
   dLp_dFi                                                                                          !< derivative of Lp with respect to Fi
 real(pReal), dimension(3,3,3,3) :: &
   dLp_dMp                                                                                          !< derivative of Lp with respect to Mandel stress
 real(pReal), dimension(9,9) :: &
   dLp_dMp99                                                                                        !< derivative of Lp with respect to Mstar (matrix notation)
 real(pReal), dimension(3,3) :: &
   Mp, &                                                                                            !< Mandel stress work conjugate with Lp
   S                                                                                                !< 2nd Piola-Kirchhoff stress
 integer(pInt) :: &
   ho, &                                                                                            !< homogenization
   tme                                                                                              !< thermal member position
 integer(pInt) :: &
   i, j, instance, of

 ho = material_homogenizationAt(el)
 tme = thermalMapping(ho)%p(ip,el)

 S  = math_Mandel6to33(S6)
 Mp  = math_mul33x33(math_mul33x33(transpose(Fi),Fi),S)

 plasticityType: select case (phase_plasticity(material_phase(ipc,ip,el)))

   case (PLASTICITY_NONE_ID) plasticityType
     Lp = 0.0_pReal
     dLp_dMp = 0.0_pReal

   case (PLASTICITY_ISOTROPIC_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     call plastic_isotropic_LpAndItsTangent       (Lp,dLp_dMp,Mp,instance,of)

   case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     call plastic_phenopowerlaw_LpAndItsTangent   (Lp,dLp_dMp,Mp,instance,of)

   case (PLASTICITY_KINEHARDENING_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     call plastic_kinehardening_LpAndItsTangent   (Lp,dLp_dMp, Mp,instance,of)

   case (PLASTICITY_NONLOCAL_ID) plasticityType
     call plastic_nonlocal_LpAndItsTangent        (Lp,dLp_dMp99, math_Mandel33to6(Mp), &
                                                   temperature(ho)%p(tme),ip,el)
     dLp_dMp = math_Plain99to3333(dLp_dMp99)                                                        ! ToDo: We revert here the last statement in plastic_xx_LpAndItsTanget

   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     call plastic_dislotwin_LpAndItsTangent       (Lp,dLp_dMp,Mp,temperature(ho)%p(tme),instance,of)

   case (PLASTICITY_DISLOUCLA_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     call plastic_disloucla_LpAndItsTangent       (Lp,dLp_dMp,Mp,temperature(ho)%p(tme),instance,of)

 end select plasticityType

#ifdef __INTEL_COMPILER
 forall(i = 1_pInt:3_pInt, j = 1_pInt:3_pInt)
#else
 do concurrent(i = 1_pInt:3_pInt, j = 1_pInt:3_pInt)
#endif
   dLp_dFi(i,j,1:3,1:3) = math_mul33x33(math_mul33x33(Fi,S),transpose(dLp_dMp(i,j,1:3,1:3))) + &
                          math_mul33x33(math_mul33x33(Fi,dLp_dMp(i,j,1:3,1:3)),S)
   dLp_dS(i,j,1:3,1:3)  = math_mul33x33(math_mul33x33(transpose(Fi),Fi),dLp_dMp(i,j,1:3,1:3))       ! ToDo: @PS: why not:   dLp_dMp:(FiT Fi)
#ifdef __INTEL_COMPILER
 end forall
#else
 enddo
#endif

end subroutine constitutive_LpAndItsTangents


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
! ToDo: MD: S is Mi?
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LiAndItsTangents(Li, dLi_dS, dLi_dFi, S6, Fi, ipc, ip, el)
 use prec, only: &
   pReal
 use math, only: &
   math_I3, &
   math_inv33, &
   math_det33, &
   math_mul33x33, &
   math_Mandel6to33
 use material, only: &
   phasememberAt, &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_plasticity, &
   material_phase, &
   phase_kinematics, &
   phase_Nkinematics, &
   PLASTICITY_isotropic_ID, &
   KINEMATICS_cleavage_opening_ID, &
   KINEMATICS_slipplane_opening_ID, &
   KINEMATICS_thermal_expansion_ID
 use plastic_isotropic, only: &
   plastic_isotropic_LiAndItsTangent
 use kinematics_cleavage_opening, only: &
   kinematics_cleavage_opening_LiAndItsTangent
 use kinematics_slipplane_opening, only: &
   kinematics_slipplane_opening_LiAndItsTangent
 use kinematics_thermal_expansion, only: &
   kinematics_thermal_expansion_LiAndItsTangent

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in),  dimension(6) :: &
   S6                                                                                               !< 2nd Piola-Kirchhoff stress (vector notation)
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   Li                                                                                               !< intermediate velocity gradient
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dLi_dS, &                                                                                        !< derivative of Li with respect to S
   dLi_dFi
 
 real(pReal), dimension(3,3) :: &
   my_Li, &                                                                                            !< intermediate velocity gradient
   FiInv, &
   temp_33
 real(pReal), dimension(3,3,3,3) :: &
   my_dLi_dS
 real(pReal) :: &
   detFi
 integer(pInt) :: &
   k, i, j, &
   instance, of

 Li = 0.0_pReal
 dLi_dS  = 0.0_pReal
 dLi_dFi = 0.0_pReal

 plasticityType: select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_isotropic_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     call plastic_isotropic_LiAndItsTangent(my_Li, my_dLi_dS, math_Mandel6to33(S6),instance,of)
   case default plasticityType
     my_Li = 0.0_pReal
     my_dLi_dS = 0.0_pReal
 end select plasticityType

 Li = Li + my_Li
 dLi_dS = dLi_dS + my_dLi_dS

 KinematicsLoop: do k = 1_pInt, phase_Nkinematics(material_phase(ipc,ip,el))
   kinematicsType: select case (phase_kinematics(k,material_phase(ipc,ip,el)))
     case (KINEMATICS_cleavage_opening_ID) kinematicsType
       call kinematics_cleavage_opening_LiAndItsTangent(my_Li, my_dLi_dS, S6, ipc, ip, el)
     case (KINEMATICS_slipplane_opening_ID) kinematicsType
       call kinematics_slipplane_opening_LiAndItsTangent(my_Li, my_dLi_dS, S6, ipc, ip, el)
     case (KINEMATICS_thermal_expansion_ID) kinematicsType
       call kinematics_thermal_expansion_LiAndItsTangent(my_Li, my_dLi_dS, ipc, ip, el)
     case default kinematicsType
       my_Li = 0.0_pReal
       my_dLi_dS = 0.0_pReal
   end select kinematicsType
   Li = Li + my_Li
   dLi_dS = dLi_dS + my_dLi_dS
 enddo KinematicsLoop

 FiInv = math_inv33(Fi)
 detFi = math_det33(Fi)
 Li = math_mul33x33(math_mul33x33(Fi,Li),FiInv)*detFi                                               !< push forward to intermediate configuration
 temp_33 = math_mul33x33(FiInv,Li)

 do i = 1_pInt,3_pInt; do j = 1_pInt,3_pInt
   dLi_dS(1:3,1:3,i,j)  = math_mul33x33(math_mul33x33(Fi,dLi_dS(1:3,1:3,i,j)),FiInv)*detFi
   dLi_dFi(1:3,1:3,i,j) = dLi_dFi(1:3,1:3,i,j) + Li*FiInv(j,i)
   dLi_dFi(1:3,i,1:3,j) = dLi_dFi(1:3,i,1:3,j) + math_I3*temp_33(j,i) + Li*FiInv(j,i)
 end do; end do

end subroutine constitutive_LiAndItsTangents


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
   KINEMATICS_thermal_expansion_ID
 use kinematics_thermal_expansion, only: &
   kinematics_thermal_expansion_initialStrain

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(3,3) :: &
   constitutive_initialFi                                                                           !< composite initial intermediate deformation gradient
 integer(pInt) :: &
   k                                                                                                !< counter in kinematics loop

 constitutive_initialFi = math_I3

 KinematicsLoop: do k = 1_pInt, phase_Nkinematics(material_phase(ipc,ip,el))                        !< Warning: small initial strain assumption
   kinematicsType: select case (phase_kinematics(k,material_phase(ipc,ip,el)))
     case (KINEMATICS_thermal_expansion_ID) kinematicsType
       constitutive_initialFi = &
         constitutive_initialFi + kinematics_thermal_expansion_initialStrain(ipc, ip, el)
   end select kinematicsType
 enddo KinematicsLoop

end function constitutive_initialFi


!--------------------------------------------------------------------------------------------------
!> @brief returns the 2nd Piola-Kirchhoff stress tensor and its tangent with respect to
!> the elastic/intermediate deformation gradients depending on the selected elastic law 
!! (so far no case switch because only Hooke is implemented)
!--------------------------------------------------------------------------------------------------
subroutine constitutive_SandItsTangents(S, dS_dFe, dS_dFi, Fe, Fi, ipc, ip, el)
 use prec, only: &
   pReal

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   S                                                                                                !< 2nd Piola-Kirchhoff stress tensor
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dS_dFe, &                                                                                        !< derivative of 2nd P-K stress with respect to elastic deformation gradient
   dS_dFi                                                                                           !< derivative of 2nd P-K stress with respect to intermediate deformation gradient

 call constitutive_hooke_SandItsTangents(S, dS_dFe, dS_dFi, Fe, Fi, ipc, ip, el)


end subroutine constitutive_SandItsTangents


!--------------------------------------------------------------------------------------------------
!> @brief returns the 2nd Piola-Kirchhoff stress tensor and its tangent with respect to
!> the elastic and intermeidate deformation gradients using Hookes law
!--------------------------------------------------------------------------------------------------
subroutine constitutive_hooke_SandItsTangents(S, dS_dFe, dS_dFi, Fe, Fi, ipc, ip, el)
 use prec, only: &
   pReal
 use math, only : &
   math_mul33x33, &
   math_mul3333xx33, &
   math_Mandel66to3333, &
   math_I3
 use material, only: &
   material_phase, &
   material_homogenizationAt, &
   phase_NstiffnessDegradations, &
   phase_stiffnessDegradation, &
   damage, &
   damageMapping, &
   STIFFNESS_DEGRADATION_damage_ID

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   S                                                                                                !< 2nd Piola-Kirchhoff stress tensor in lattice configuration
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dS_dFe, &                                                                                        !< derivative of 2nd P-K stress with respect to elastic deformation gradient
   dS_dFi                                                                                           !< derivative of 2nd P-K stress with respect to intermediate deformation gradient
 real(pReal), dimension(3,3) :: E
 real(pReal), dimension(3,3,3,3) :: C
 integer(pInt) :: &
   ho, &                                                                                            !< homogenization
   d                                                                                                !< counter in degradation loop
 integer(pInt) :: &
   i, j

 ho = material_homogenizationAt(el)
 C = math_Mandel66to3333(constitutive_homogenizedC(ipc,ip,el))

 DegradationLoop: do d = 1_pInt, phase_NstiffnessDegradations(material_phase(ipc,ip,el))
   degradationType: select case(phase_stiffnessDegradation(d,material_phase(ipc,ip,el)))
     case (STIFFNESS_DEGRADATION_damage_ID) degradationType
       C = C * damage(ho)%p(damageMapping(ho)%p(ip,el))**2_pInt
   end select degradationType
 enddo DegradationLoop

 E = 0.5_pReal*(math_mul33x33(transpose(Fe),Fe)-math_I3)                                            !< Green-Lagrange strain in unloaded configuration
 S = math_mul3333xx33(C,math_mul33x33(math_mul33x33(transpose(Fi),E),Fi))                           !< 2PK stress in lattice configuration in work conjugate with GL strain pulled back to lattice configuration

 dS_dFe = 0.0_pReal
 forall (i=1_pInt:3_pInt, j=1_pInt:3_pInt)
   dS_dFe(i,j,1:3,1:3) = &
     math_mul33x33(Fe,math_mul33x33(math_mul33x33(Fi,C(i,j,1:3,1:3)),transpose(Fi)))                !< dS_ij/dFe_kl = C_ijmn * Fi_lm * Fi_on * Fe_ko
   dS_dFi(i,j,1:3,1:3) = 2.0_pReal*math_mul33x33(math_mul33x33(E,Fi),C(i,j,1:3,1:3))                !< dS_ij/dFi_kl = C_ijln * E_km * Fe_mn
 end forall

end subroutine constitutive_hooke_SandItsTangents


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine constitutive_collectDotState(S6, FeArray, Fi, FpArray, subdt, subfracArray,ipc, ip, el)
 use prec, only: &
   pReal, &
   pLongInt
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use math, only: &
   math_mul33x33, &
   math_Mandel6to33, &
   math_Mandel33to6, &
   math_mul33x33
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   phasememberAt, &
   phase_plasticityInstance, &
   phase_plasticity, &
   phase_source, &
   phase_Nsources, &
   material_phase, &
   material_homogenizationAt, &
   temperature, &
   thermalMapping, &
   homogenization_maxNgrains, &
   PLASTICITY_none_ID, &
   PLASTICITY_isotropic_ID, &
   PLASTICITY_phenopowerlaw_ID, &
   PLASTICITY_kinehardening_ID, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_nonlocal_ID, &
   SOURCE_damage_isoDuctile_ID, &
   SOURCE_damage_anisoBrittle_ID, &
   SOURCE_damage_anisoDuctile_ID, &
   SOURCE_thermal_externalheat_ID
 use plastic_isotropic, only:  &
   plastic_isotropic_dotState
 use plastic_phenopowerlaw, only: &
   plastic_phenopowerlaw_dotState
 use plastic_kinehardening, only: &
   plastic_kinehardening_dotState  
 use plastic_dislotwin, only: &
   plastic_dislotwin_dotState
 use plastic_disloucla, only: &
   plastic_disloucla_dotState
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
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),  intent(in) :: &
   subdt                                                                                            !< timestep
 real(pReal),  intent(in), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   subfracArray                                                                                     !< subfraction of timestep
 real(pReal),  intent(in), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   FeArray, &                                                                                       !< elastic deformation gradient
   FpArray                                                                                          !< plastic deformation gradient
 real(pReal),  intent(in), dimension(3,3) :: &
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),  intent(in), dimension(6) :: &
   S6                                                                                               !< 2nd Piola Kirchhoff stress (vector notation)
 real(pReal),              dimension(3,3) :: &
   Mp
 integer(pInt) :: &
   ho, &                                                                                            !< homogenization
   tme, &                                                                                           !< thermal member position
   s, &                                                                                                !< counter in source loop
   instance, of

 ho = material_homogenizationAt(el)
 tme = thermalMapping(ho)%p(ip,el)

 Mp  = math_mul33x33(math_mul33x33(transpose(Fi),Fi),math_Mandel6to33(S6))

 plasticityType: select case (phase_plasticity(material_phase(ipc,ip,el)))

   case (PLASTICITY_ISOTROPIC_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     call plastic_isotropic_dotState    (Mp,instance,of)

   case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     call plastic_phenopowerlaw_dotState(Mp,instance,of)

   case (PLASTICITY_KINEHARDENING_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     call plastic_kinehardening_dotState(Mp,instance,of)

   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     call plastic_dislotwin_dotState    (Mp,temperature(ho)%p(tme),instance,of)

   case (PLASTICITY_DISLOUCLA_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     call plastic_disloucla_dotState    (Mp,temperature(ho)%p(tme),instance,of)

   case (PLASTICITY_NONLOCAL_ID) plasticityType
     call plastic_nonlocal_dotState     (math_Mandel33to6(Mp),FeArray,FpArray,temperature(ho)%p(tme), &
                                         subdt,subfracArray,ip,el)
 end select plasticityType

 SourceLoop: do s = 1_pInt, phase_Nsources(material_phase(ipc,ip,el))

   sourceType: select case (phase_source(s,material_phase(ipc,ip,el)))

     case (SOURCE_damage_anisoBrittle_ID) sourceType
       call source_damage_anisoBrittle_dotState (S6, ipc, ip, el) !< correct stress?

     case (SOURCE_damage_isoDuctile_ID) sourceType
       call source_damage_isoDuctile_dotState   (         ipc, ip, el)

     case (SOURCE_damage_anisoDuctile_ID) sourceType
       call source_damage_anisoDuctile_dotState (         ipc, ip, el)

     case (SOURCE_thermal_externalheat_ID) sourceType
       call source_thermal_externalheat_dotState(         ipc, ip, el)

   end select sourceType

 enddo SourceLoop

end subroutine constitutive_collectDotState

!--------------------------------------------------------------------------------------------------
!> @brief for constitutive models having an instantaneous change of state
!> will return false if delta state is not needed/supported by the constitutive model
!--------------------------------------------------------------------------------------------------
subroutine constitutive_collectDeltaState(S6, Fe, Fi, ipc, ip, el)
 use prec, only: &
   pReal, &
   pLongInt
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use math, only: &
   math_Mandel6to33, &
   math_Mandel33to6, &
   math_mul33x33
 use material, only: &
   phasememberAt, &
   phase_plasticityInstance, &
   phase_plasticity, &
   phase_source, &
   phase_Nsources, &
   material_phase, &
   PLASTICITY_KINEHARDENING_ID, &
   PLASTICITY_NONLOCAL_ID, &
   SOURCE_damage_isoBrittle_ID
 use plastic_kinehardening, only: &
   plastic_kinehardening_deltaState   
 use plastic_nonlocal, only: &
   plastic_nonlocal_deltaState
 use source_damage_isoBrittle, only: &
   source_damage_isoBrittle_deltaState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in),  dimension(6) :: &
   S6                                                                                               !< 2nd Piola Kirchhoff stress (vector notation)
 real(pReal),   intent(in), dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),               dimension(3,3) :: &
   Mp
 integer(pInt) :: &
   s, &                                                                                                !< counter in source loop
   instance, of

 Mp  = math_mul33x33(math_mul33x33(transpose(Fi),Fi),math_Mandel6to33(S6))

 plasticityType: select case (phase_plasticity(material_phase(ipc,ip,el)))

   case (PLASTICITY_KINEHARDENING_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     call plastic_kinehardening_deltaState(Mp,instance,of)

   case (PLASTICITY_NONLOCAL_ID) plasticityType
     call plastic_nonlocal_deltaState(math_Mandel33to6(Mp),ip,el)

 end select plasticityType

 sourceLoop: do s = 1_pInt, phase_Nsources(material_phase(ipc,ip,el))

    sourceType: select case (phase_source(s,material_phase(ipc,ip,el)))

     case (SOURCE_damage_isoBrittle_ID) sourceType
       call source_damage_isoBrittle_deltaState  (constitutive_homogenizedC(ipc,ip,el), Fe, &
                                                  ipc, ip, el)

   end select sourceType

 enddo SourceLoop

end subroutine constitutive_collectDeltaState


!--------------------------------------------------------------------------------------------------
!> @brief returns array of constitutive results
!--------------------------------------------------------------------------------------------------
function constitutive_postResults(S6, Fi, FeArray, ipc, ip, el)
 use prec, only: &
   pReal
 use math, only: &
  math_Mandel6to33, &
  math_mul33x33
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   phasememberAt, &
   phase_plasticityInstance, &
   plasticState, &
   sourceState, &
   phase_plasticity, &
   phase_source, &
   phase_Nsources, &
   material_phase, &
   material_homogenizationAt, &
   temperature, &
   thermalMapping, &
   homogenization_maxNgrains, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_ISOTROPIC_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_KINEHARDENING_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOUCLA_ID, &
   PLASTICITY_NONLOCAL_ID, &
   SOURCE_damage_isoBrittle_ID, &
   SOURCE_damage_isoDuctile_ID, &
   SOURCE_damage_anisoBrittle_ID, &
   SOURCE_damage_anisoDuctile_ID
 use plastic_isotropic, only: &
   plastic_isotropic_postResults
 use plastic_phenopowerlaw, only: &
   plastic_phenopowerlaw_postResults
 use plastic_kinehardening, only: &
   plastic_kinehardening_postResults
 use plastic_dislotwin, only: &
   plastic_dislotwin_postResults
 use plastic_disloucla, only: &
   plastic_disloucla_postResults
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
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%sizePostResults + &
                        sum(sourceState(material_phase(ipc,ip,el))%p(:)%sizePostResults)) :: &
   constitutive_postResults
 real(pReal),  intent(in), dimension(3,3) :: &
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),  intent(in), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   FeArray                                                                                          !< elastic deformation gradient
 real(pReal),  intent(in), dimension(6) :: &
   S6                                                                                               !< 2nd Piola Kirchhoff stress (vector notation)
 real(pReal), dimension(3,3) :: &
   Mp                                                                                               !< Mandel stress
 integer(pInt) :: &
   startPos, endPos
 integer(pInt) :: &
   ho, &                                                                                            !< homogenization
   tme, &                                                                                           !< thermal member position
   s, of, instance                                                                                  !< counter in source loop

 constitutive_postResults = 0.0_pReal

 Mp  = math_mul33x33(math_mul33x33(transpose(Fi),Fi),math_Mandel6to33(S6))

 ho = material_homogenizationAt(el)
 tme = thermalMapping(ho)%p(ip,el)

 startPos = 1_pInt
 endPos = plasticState(material_phase(ipc,ip,el))%sizePostResults

 plasticityType: select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_ISOTROPIC_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     constitutive_postResults(startPos:endPos) = &
       plastic_isotropic_postResults(Mp,instance,of)

   case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     constitutive_postResults(startPos:endPos) = &
       plastic_phenopowerlaw_postResults(Mp,instance,of)

   case (PLASTICITY_KINEHARDENING_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     constitutive_postResults(startPos:endPos) = &
       plastic_kinehardening_postResults(Mp,instance,of)

   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     constitutive_postResults(startPos:endPos) = &
       plastic_dislotwin_postResults(Mp,temperature(ho)%p(tme),instance,of)

   case (PLASTICITY_DISLOUCLA_ID) plasticityType
     of = phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phase(ipc,ip,el))
     constitutive_postResults(startPos:endPos) = &
       plastic_disloucla_postResults(Mp,temperature(ho)%p(tme),instance,of)

   case (PLASTICITY_NONLOCAL_ID) plasticityType
     constitutive_postResults(startPos:endPos) = &
       plastic_nonlocal_postResults (S6,FeArray,ip,el)
 end select plasticityType

 SourceLoop: do s = 1_pInt, phase_Nsources(material_phase(ipc,ip,el))
   startPos = endPos + 1_pInt
   endPos = endPos + sourceState(material_phase(ipc,ip,el))%p(s)%sizePostResults
   sourceType: select case (phase_source(s,material_phase(ipc,ip,el)))
     case (SOURCE_damage_isoBrittle_ID) sourceType
       constitutive_postResults(startPos:endPos) = source_damage_isoBrittle_postResults(ipc, ip, el)
     case (SOURCE_damage_isoDuctile_ID) sourceType
       constitutive_postResults(startPos:endPos) = source_damage_isoDuctile_postResults(ipc, ip, el)
     case (SOURCE_damage_anisoBrittle_ID) sourceType
       constitutive_postResults(startPos:endPos) = source_damage_anisoBrittle_postResults(ipc, ip, el)
     case (SOURCE_damage_anisoDuctile_ID) sourceType
       constitutive_postResults(startPos:endPos) = source_damage_anisoDuctile_postResults(ipc, ip, el)
   end select sourceType
 enddo SourceLoop

end function constitutive_postResults


!--------------------------------------------------------------------------------------------------
!> @brief writes constitutive results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine constitutive_results()
 use material, only: &
   PLASTICITY_ISOTROPIC_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_KINEHARDENING_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOUCLA_ID, &
   PLASTICITY_NONLOCAL_ID
#if defined(PETSc) || defined(DAMASKHDF5)
 use results
 use HDF5_utilities
 use config, only: &
   config_name_phase => phase_name                                                                  ! anticipate logical name
   
 use material, only: &
   phase_plasticityInstance, &
   material_phase_plasticity_type => phase_plasticity
 use plastic_phenopowerlaw, only: &
   plastic_phenopowerlaw_results
 
 implicit none
 integer(pInt) :: p  
 call HDF5_closeGroup(results_addGroup('current/phase'))                                              
 do p=1,size(config_name_phase)                                                                           
   call HDF5_closeGroup(results_addGroup('current/phase/'//trim(config_name_phase(p))))
   if (material_phase_plasticity_type(p) == PLASTICITY_PHENOPOWERLAW_ID) then
     call plastic_phenopowerlaw_results(phase_plasticityInstance(p),'current/phase/'//trim(config_name_phase(p)))
   endif
 enddo      

#endif


end subroutine constitutive_results

end module constitutive
