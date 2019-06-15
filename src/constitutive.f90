!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief elasticity, plasticity, internal microstructure state
!--------------------------------------------------------------------------------------------------
module constitutive
 use math
 use debug
 use numerics
 use IO
 use config
 use material
 use results
 use HDF5_utilities
 use lattice
 use mesh
 use discretization
 use plastic_none
 use plastic_isotropic
 use plastic_phenopowerlaw
 use plastic_kinehardening
 use plastic_dislotwin
 use plastic_disloucla
 use plastic_nonlocal
 use geometry_plastic_nonlocal
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
 private
 
 integer, public, protected :: &
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
subroutine constitutive_init

 integer, parameter :: FILEUNIT = 204
 integer :: &
   o, &                                                                                             !< counter in output loop
   ph, &                                                                                            !< counter in phase loop
   s, &                                                                                             !< counter in source loop
   ins                                                                                              !< instance of plasticity/source

 integer, dimension(:,:), pointer :: thisSize
 character(len=64), dimension(:,:), pointer :: thisOutput
 character(len=32) :: outputName                                                                    !< name of output, intermediate fix until HDF5 output is ready
 logical :: knownPlasticity, knownSource, nonlocalConstitutionPresent
 nonlocalConstitutionPresent = .false.

!--------------------------------------------------------------------------------------------------
! initialized plasticity
 if (any(phase_plasticity == PLASTICITY_NONE_ID))          call plastic_none_init
 if (any(phase_plasticity == PLASTICITY_ISOTROPIC_ID))     call plastic_isotropic_init
 if (any(phase_plasticity == PLASTICITY_PHENOPOWERLAW_ID)) call plastic_phenopowerlaw_init
 if (any(phase_plasticity == PLASTICITY_KINEHARDENING_ID)) call plastic_kinehardening_init
 if (any(phase_plasticity == PLASTICITY_DISLOTWIN_ID))     call plastic_dislotwin_init
 if (any(phase_plasticity == PLASTICITY_DISLOUCLA_ID))     call plastic_disloucla_init
 if (any(phase_plasticity == PLASTICITY_NONLOCAL_ID)) then
   call plastic_nonlocal_init
 else
   call geometry_plastic_nonlocal_disable
 endif
!--------------------------------------------------------------------------------------------------
! initialize source mechanisms
 if (any(phase_source == SOURCE_thermal_dissipation_ID))     call source_thermal_dissipation_init
 if (any(phase_source == SOURCE_thermal_externalheat_ID))    call source_thermal_externalheat_init
 if (any(phase_source == SOURCE_damage_isoBrittle_ID))       call source_damage_isoBrittle_init
 if (any(phase_source == SOURCE_damage_isoDuctile_ID))       call source_damage_isoDuctile_init
 if (any(phase_source == SOURCE_damage_anisoBrittle_ID))     call source_damage_anisoBrittle_init
 if (any(phase_source == SOURCE_damage_anisoDuctile_ID))     call source_damage_anisoDuctile_init
 
!--------------------------------------------------------------------------------------------------
! initialize kinematic mechanisms
 if (any(phase_kinematics == KINEMATICS_cleavage_opening_ID))  call kinematics_cleavage_opening_init
 if (any(phase_kinematics == KINEMATICS_slipplane_opening_ID)) call kinematics_slipplane_opening_init
 if (any(phase_kinematics == KINEMATICS_thermal_expansion_ID)) call kinematics_thermal_expansion_init

 write(6,'(/,a)')   ' <<<+-  constitutive init  -+>>>'

 mainProcess: if (worldrank == 0) then
!--------------------------------------------------------------------------------------------------
! write description file for constitutive output
   call IO_write_jobFile(FILEUNIT,'outputConstitutive')
   PhaseLoop: do ph = 1,material_Nphase
     activePhase: if (any(material_phaseAt == ph)) then
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
           OutputPlasticityLoop: do o = 1,size(thisOutput(:,ins))
             if(len(trim(thisOutput(o,ins))) > 0) &
               write(FILEUNIT,'(a,i4)') trim(thisOutput(o,ins))//char(9),thisSize(o,ins)
           enddo OutputPlasticityLoop
         endif
       endif
       
       SourceLoop: do s = 1, phase_Nsources(ph)
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
           OutputSourceLoop: do o = 1,size(thisOutput(:,ins))
             if(len(trim(thisOutput(o,ins))) > 0) &
               write(FILEUNIT,'(a,i4)') trim(thisOutput(o,ins))//char(9),thisSize(o,ins)
           enddo OutputSourceLoop
         endif
       enddo SourceLoop
     endif activePhase
   enddo PhaseLoop
   close(FILEUNIT)
 endif mainProcess

 constitutive_plasticity_maxSizeDotState = 0
 constitutive_plasticity_maxSizePostResults = 0
 constitutive_source_maxSizeDotState = 0
 constitutive_source_maxSizePostResults = 0

 PhaseLoop2:do ph = 1,material_Nphase
!--------------------------------------------------------------------------------------------------
! partition and inititalize state
   plasticState(ph)%partionedState0 = plasticState(ph)%state0
   plasticState(ph)%state           = plasticState(ph)%partionedState0
   forall(s = 1:phase_Nsources(ph))
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

 real(pReal), dimension(6,6) :: constitutive_homogenizedC
 integer, intent(in) :: &
   ipc, &                                                                                            !< component-ID of integration point
   ip, &                                                                                             !< integration point
   el                                                                                                !< element

 plasticityType: select case (phase_plasticity(material_phaseAt(ipc,el)))
   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     constitutive_homogenizedC = plastic_dislotwin_homogenizedC(ipc,ip,el)
   case default plasticityType
     constitutive_homogenizedC = lattice_C66(1:6,1:6,material_phaseAt(ipc,el))
 end select plasticityType

end function constitutive_homogenizedC

!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different constitutive models
!--------------------------------------------------------------------------------------------------
subroutine constitutive_microstructure(Fe, Fp, ipc, ip, el)

 integer, intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in), dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fp                                                                                               !< plastic deformation gradient
 integer :: &
   ho, &                                                                                            !< homogenization
   tme, &                                                                                           !< thermal member position
   instance, of

 ho = material_homogenizationAt(el)
 tme = thermalMapping(ho)%p(ip,el)

 plasticityType: select case (phase_plasticity(material_phaseAt(ipc,el)))
   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_dislotwin_dependentState(temperature(ho)%p(tme),instance,of)
   case (PLASTICITY_DISLOUCLA_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_disloUCLA_dependentState(instance,of)
   case (PLASTICITY_NONLOCAL_ID) plasticityType
     call plastic_nonlocal_dependentState (Fe,Fp,ip,el)
 end select plasticityType

end subroutine constitutive_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
! ToDo: Discuss wheter it makes sense if crystallite handles the configuration conversion, i.e.
! Mp in, dLp_dMp out
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LpAndItsTangents(Lp, dLp_dS, dLp_dFi, &
                                         S, Fi, ipc, ip, el)

 integer, intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in),  dimension(3,3) :: &
   S, &                                                                                             !< 2nd Piola-Kirchhoff stress
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dLp_dS, &
   dLp_dFi                                                                                          !< derivative of Lp with respect to Fi
 real(pReal), dimension(3,3,3,3) :: &
   dLp_dMp                                                                                          !< derivative of Lp with respect to Mandel stress
 real(pReal), dimension(3,3) :: &
   Mp                                                                                               !< Mandel stress work conjugate with Lp
 integer :: &
   ho, &                                                                                            !< homogenization
   tme                                                                                              !< thermal member position
 integer :: &
   i, j, instance, of

 ho = material_homogenizationAt(el)
 tme = thermalMapping(ho)%p(ip,el)

 Mp  = matmul(matmul(transpose(Fi),Fi),S)

 plasticityType: select case (phase_plasticity(material_phaseAt(ipc,el)))

   case (PLASTICITY_NONE_ID) plasticityType
     Lp = 0.0_pReal
     dLp_dMp = 0.0_pReal

   case (PLASTICITY_ISOTROPIC_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_isotropic_LpAndItsTangent       (Lp,dLp_dMp,Mp,instance,of)

   case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_phenopowerlaw_LpAndItsTangent   (Lp,dLp_dMp,Mp,instance,of)

   case (PLASTICITY_KINEHARDENING_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_kinehardening_LpAndItsTangent   (Lp,dLp_dMp, Mp,instance,of)

   case (PLASTICITY_NONLOCAL_ID) plasticityType
     call plastic_nonlocal_LpAndItsTangent        (Lp,dLp_dMp,Mp, &
                                                   temperature(ho)%p(tme),geometry_plastic_nonlocal_IPvolume0(ip,el),ip,el)

   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_dislotwin_LpAndItsTangent       (Lp,dLp_dMp,Mp,temperature(ho)%p(tme),instance,of)

   case (PLASTICITY_DISLOUCLA_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_disloucla_LpAndItsTangent       (Lp,dLp_dMp,Mp,temperature(ho)%p(tme),instance,of)

 end select plasticityType

 do i=1,3; do j=1,3
   dLp_dFi(i,j,1:3,1:3) = matmul(matmul(Fi,S),transpose(dLp_dMp(i,j,1:3,1:3))) + &
                          matmul(matmul(Fi,dLp_dMp(i,j,1:3,1:3)),S)
   dLp_dS(i,j,1:3,1:3)  = matmul(matmul(transpose(Fi),Fi),dLp_dMp(i,j,1:3,1:3))                     ! ToDo: @PS: why not:   dLp_dMp:(FiT Fi)
 enddo; enddo

end subroutine constitutive_LpAndItsTangents


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
! ToDo: MD: S is Mi?
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LiAndItsTangents(Li, dLi_dS, dLi_dFi, &
                                         S, Fi, ipc, ip, el)

 integer, intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in),  dimension(3,3) :: &
   S                                                                                                !< 2nd Piola-Kirchhoff stress
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
 integer :: &
   k, i, j, &
   instance, of

 Li = 0.0_pReal
 dLi_dS  = 0.0_pReal
 dLi_dFi = 0.0_pReal

 plasticityType: select case (phase_plasticity(material_phaseAt(ipc,el)))
   case (PLASTICITY_isotropic_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_isotropic_LiAndItsTangent(my_Li, my_dLi_dS, S ,instance,of)
   case default plasticityType
     my_Li = 0.0_pReal
     my_dLi_dS = 0.0_pReal
 end select plasticityType

 Li = Li + my_Li
 dLi_dS = dLi_dS + my_dLi_dS

 KinematicsLoop: do k = 1, phase_Nkinematics(material_phaseAt(ipc,el))
   kinematicsType: select case (phase_kinematics(k,material_phaseAt(ipc,el)))
     case (KINEMATICS_cleavage_opening_ID) kinematicsType
       call kinematics_cleavage_opening_LiAndItsTangent(my_Li, my_dLi_dS, S, ipc, ip, el)
     case (KINEMATICS_slipplane_opening_ID) kinematicsType
       call kinematics_slipplane_opening_LiAndItsTangent(my_Li, my_dLi_dS, S, ipc, ip, el)
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
 Li = matmul(matmul(Fi,Li),FiInv)*detFi                                               !< push forward to intermediate configuration
 temp_33 = matmul(FiInv,Li)

 do i = 1,3; do j = 1,3
   dLi_dS(1:3,1:3,i,j)  = matmul(matmul(Fi,dLi_dS(1:3,1:3,i,j)),FiInv)*detFi
   dLi_dFi(1:3,1:3,i,j) = dLi_dFi(1:3,1:3,i,j) + Li*FiInv(j,i)
   dLi_dFi(1:3,i,1:3,j) = dLi_dFi(1:3,i,1:3,j) + math_I3*temp_33(j,i) + Li*FiInv(j,i)
 end do; end do

end subroutine constitutive_LiAndItsTangents


!--------------------------------------------------------------------------------------------------
!> @brief  collects initial intermediate deformation gradient
!--------------------------------------------------------------------------------------------------
pure function constitutive_initialFi(ipc, ip, el)

 integer, intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(3,3) :: &
   constitutive_initialFi                                                                           !< composite initial intermediate deformation gradient
 integer :: &
   k                                                                                                !< counter in kinematics loop
 integer :: &
   phase, &
   homog, offset

 constitutive_initialFi = math_I3
 phase = material_phaseAt(ipc,el)

 KinematicsLoop: do k = 1, phase_Nkinematics(phase)                                            !< Warning: small initial strain assumption
   kinematicsType: select case (phase_kinematics(k,phase))
     case (KINEMATICS_thermal_expansion_ID) kinematicsType
       homog = material_homogenizationAt(el)
       offset = thermalMapping(homog)%p(ip,el)
       constitutive_initialFi = &
         constitutive_initialFi + kinematics_thermal_expansion_initialStrain(homog,phase,offset)
   end select kinematicsType
 enddo KinematicsLoop

end function constitutive_initialFi


!--------------------------------------------------------------------------------------------------
!> @brief returns the 2nd Piola-Kirchhoff stress tensor and its tangent with respect to
!> the elastic/intermediate deformation gradients depending on the selected elastic law 
!! (so far no case switch because only Hooke is implemented)
!--------------------------------------------------------------------------------------------------
subroutine constitutive_SandItsTangents(S, dS_dFe, dS_dFi, Fe, Fi, ipc, ip, el)

 integer, intent(in) :: &
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
subroutine constitutive_hooke_SandItsTangents(S, dS_dFe, dS_dFi, &
                                              Fe, Fi, ipc, ip, el)

 integer, intent(in) :: &
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
 integer :: &
   ho, &                                                                                            !< homogenization
   d                                                                                                !< counter in degradation loop
 integer :: &
   i, j

 ho = material_homogenizationAt(el)
 C = math_66toSym3333(constitutive_homogenizedC(ipc,ip,el))

 DegradationLoop: do d = 1, phase_NstiffnessDegradations(material_phaseAt(ipc,el))
   degradationType: select case(phase_stiffnessDegradation(d,material_phaseAt(ipc,el)))
     case (STIFFNESS_DEGRADATION_damage_ID) degradationType
       C = C * damage(ho)%p(damageMapping(ho)%p(ip,el))**2
   end select degradationType
 enddo DegradationLoop

 E = 0.5_pReal*(matmul(transpose(Fe),Fe)-math_I3)                                            !< Green-Lagrange strain in unloaded configuration
 S = math_mul3333xx33(C,matmul(matmul(transpose(Fi),E),Fi))                           !< 2PK stress in lattice configuration in work conjugate with GL strain pulled back to lattice configuration

 dS_dFe = 0.0_pReal
 forall (i=1:3, j=1:3)
   dS_dFe(i,j,1:3,1:3) = &
     matmul(Fe,matmul(matmul(Fi,C(i,j,1:3,1:3)),transpose(Fi)))                !< dS_ij/dFe_kl = C_ijmn * Fi_lm * Fi_on * Fe_ko
   dS_dFi(i,j,1:3,1:3) = 2.0_pReal*matmul(matmul(E,Fi),C(i,j,1:3,1:3))                !< dS_ij/dFi_kl = C_ijln * E_km * Fe_mn
 end forall

end subroutine constitutive_hooke_SandItsTangents


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine constitutive_collectDotState(S, FeArray, Fi, FpArray, subdt, ipc, ip, el)

 integer, intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),  intent(in) :: &
   subdt                                                                                            !< timestep
 real(pReal),  intent(in), dimension(3,3,homogenization_maxNgrains,discretization_nIP,discretization_nElem) :: &
   FeArray, &                                                                                       !< elastic deformation gradient
   FpArray                                                                                          !< plastic deformation gradient
 real(pReal),  intent(in), dimension(3,3) :: &
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),  intent(in), dimension(3,3) :: &
   S                                                                                                !< 2nd Piola Kirchhoff stress (vector notation)
 real(pReal),              dimension(3,3) :: &
   Mp
 integer :: &
   ho, &                                                                                            !< homogenization
   tme, &                                                                                           !< thermal member position
   i, &                                                                                             !< counter in source loop
   instance, of

 ho = material_homogenizationAt(el)
 tme = thermalMapping(ho)%p(ip,el)

 Mp  = matmul(matmul(transpose(Fi),Fi),S)

 plasticityType: select case (phase_plasticity(material_phaseAt(ipc,el)))

   case (PLASTICITY_ISOTROPIC_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_isotropic_dotState    (Mp,instance,of)

   case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_phenopowerlaw_dotState(Mp,instance,of)

   case (PLASTICITY_KINEHARDENING_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_kinehardening_dotState(Mp,instance,of)

   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_dislotwin_dotState    (Mp,temperature(ho)%p(tme),instance,of)

   case (PLASTICITY_DISLOUCLA_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_disloucla_dotState    (Mp,temperature(ho)%p(tme),instance,of)

   case (PLASTICITY_NONLOCAL_ID) plasticityType
     call plastic_nonlocal_dotState     (Mp,FeArray,FpArray,temperature(ho)%p(tme), &
                                         subdt,ip,el)
 end select plasticityType

 SourceLoop: do i = 1, phase_Nsources(material_phaseAt(ipc,el))

   sourceType: select case (phase_source(i,material_phaseAt(ipc,el)))

     case (SOURCE_damage_anisoBrittle_ID) sourceType
       call source_damage_anisoBrittle_dotState (S, ipc, ip, el) !< correct stress?

     case (SOURCE_damage_isoDuctile_ID) sourceType
       call source_damage_isoDuctile_dotState   (         ipc, ip, el)

     case (SOURCE_damage_anisoDuctile_ID) sourceType
       call source_damage_anisoDuctile_dotState (         ipc, ip, el)

     case (SOURCE_thermal_externalheat_ID) sourceType
       of = material_phasememberAt(ipc,ip,el)
       call source_thermal_externalheat_dotState(material_phaseAt(ipc,el),of)

   end select sourceType

 enddo SourceLoop

end subroutine constitutive_collectDotState

!--------------------------------------------------------------------------------------------------
!> @brief for constitutive models having an instantaneous change of state
!> will return false if delta state is not needed/supported by the constitutive model
!--------------------------------------------------------------------------------------------------
subroutine constitutive_collectDeltaState(S, Fe, Fi, ipc, ip, el)

 integer, intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in), dimension(3,3) :: &
   S, &                                                                                             !< 2nd Piola Kirchhoff stress
   Fe, &                                                                                            !< elastic deformation gradient
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),               dimension(3,3) :: &
   Mp
 integer :: &
   i, &
   instance, of

 Mp  = matmul(matmul(transpose(Fi),Fi),S)

 plasticityType: select case (phase_plasticity(material_phaseAt(ipc,el)))

   case (PLASTICITY_KINEHARDENING_ID) plasticityType
     of = material_phasememberAt(ipc,ip,el)
     instance = phase_plasticityInstance(material_phaseAt(ipc,el))
     call plastic_kinehardening_deltaState(Mp,instance,of)

   case (PLASTICITY_NONLOCAL_ID) plasticityType
     call plastic_nonlocal_deltaState(Mp,ip,el)

 end select plasticityType

 sourceLoop: do i = 1, phase_Nsources(material_phaseAt(ipc,el))

    sourceType: select case (phase_source(i,material_phaseAt(ipc,el)))

     case (SOURCE_damage_isoBrittle_ID) sourceType
       call source_damage_isoBrittle_deltaState  (constitutive_homogenizedC(ipc,ip,el), Fe, &
                                                  ipc, ip, el)

   end select sourceType

 enddo SourceLoop

end subroutine constitutive_collectDeltaState


!--------------------------------------------------------------------------------------------------
!> @brief returns array of constitutive results
!--------------------------------------------------------------------------------------------------
function constitutive_postResults(S, Fi, ipc, ip, el)

 integer, intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(plasticState(material_phaseAt(ipc,el))%sizePostResults + &
                        sum(sourceState(material_phaseAt(ipc,el))%p(:)%sizePostResults)) :: &
   constitutive_postResults
 real(pReal),  intent(in), dimension(3,3) :: &
   Fi                                                                                               !< intermediate deformation gradient
 real(pReal),  intent(in), dimension(3,3) :: &
   S                                                                                                !< 2nd Piola Kirchhoff stress
 real(pReal), dimension(3,3) :: &
   Mp                                                                                               !< Mandel stress
 integer :: &
   startPos, endPos
 integer :: &
   ho, &                                                                                            !< homogenization
   tme, &                                                                                           !< thermal member position
   i, of, instance                                                                                  !< counter in source loop

 constitutive_postResults = 0.0_pReal

 Mp  = matmul(matmul(transpose(Fi),Fi),S)

 ho = material_homogenizationAt(el)
 tme = thermalMapping(ho)%p(ip,el)

 startPos = 1
 endPos = plasticState(material_phaseAt(ipc,el))%sizePostResults

 of = material_phasememberAt(ipc,ip,el)
 instance = phase_plasticityInstance(material_phaseAt(ipc,el))

 plasticityType: select case (phase_plasticity(material_phaseAt(ipc,el)))
   case (PLASTICITY_ISOTROPIC_ID) plasticityType
     constitutive_postResults(startPos:endPos) = &
       plastic_isotropic_postResults(Mp,instance,of)

   case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
     constitutive_postResults(startPos:endPos) = &
       plastic_phenopowerlaw_postResults(Mp,instance,of)

   case (PLASTICITY_KINEHARDENING_ID) plasticityType
     constitutive_postResults(startPos:endPos) = &
       plastic_kinehardening_postResults(Mp,instance,of)

   case (PLASTICITY_DISLOTWIN_ID) plasticityType
     constitutive_postResults(startPos:endPos) = &
       plastic_dislotwin_postResults(Mp,temperature(ho)%p(tme),instance,of)

   case (PLASTICITY_DISLOUCLA_ID) plasticityType
     constitutive_postResults(startPos:endPos) = &
       plastic_disloucla_postResults(Mp,temperature(ho)%p(tme),instance,of)

   case (PLASTICITY_NONLOCAL_ID) plasticityType
     constitutive_postResults(startPos:endPos) = &
       plastic_nonlocal_postResults (material_phaseAt(ipc,el),instance,of)

 end select plasticityType

 SourceLoop: do i = 1, phase_Nsources(material_phaseAt(ipc,el))
   startPos = endPos + 1
   endPos = endPos + sourceState(material_phaseAt(ipc,el))%p(i)%sizePostResults
   of = material_phasememberAt(ipc,ip,el)
   sourceType: select case (phase_source(i,material_phaseAt(ipc,el)))
     case (SOURCE_damage_isoBrittle_ID) sourceType
       constitutive_postResults(startPos:endPos) = source_damage_isoBrittle_postResults(material_phaseAt(ipc,el),of)
     case (SOURCE_damage_isoDuctile_ID) sourceType
       constitutive_postResults(startPos:endPos) = source_damage_isoDuctile_postResults(material_phaseAt(ipc,el),of)
     case (SOURCE_damage_anisoBrittle_ID) sourceType
       constitutive_postResults(startPos:endPos) = source_damage_anisoBrittle_postResults(material_phaseAt(ipc,el),of)
     case (SOURCE_damage_anisoDuctile_ID) sourceType
       constitutive_postResults(startPos:endPos) = source_damage_anisoDuctile_postResults(material_phaseAt(ipc,el),of)
   end select sourceType

 enddo SourceLoop

end function constitutive_postResults


!--------------------------------------------------------------------------------------------------
!> @brief writes constitutive results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine constitutive_results
        
  integer :: p
  character(len=256) :: group
#if defined(PETSc) || defined(DAMASK_HDF5)                                             
  do p=1,size(phase_name)
    group = trim('current/constituent')//'/'//trim(phase_name(p))
    call HDF5_closeGroup(results_addGroup(group))
    
    group = trim(group)//'/plastic'
    
    call HDF5_closeGroup(results_addGroup(group))  
    select case(phase_plasticity(p))
    
      case(PLASTICITY_ISOTROPIC_ID)
        call plastic_isotropic_results(phase_plasticityInstance(p),group) 
        
      case(PLASTICITY_PHENOPOWERLAW_ID)
        call plastic_phenopowerlaw_results(phase_plasticityInstance(p),group) 
         
      case(PLASTICITY_KINEHARDENING_ID)
        call plastic_kinehardening_results(phase_plasticityInstance(p),group) 

      case(PLASTICITY_DISLOTWIN_ID)
        call plastic_dislotwin_results(phase_plasticityInstance(p),group) 
      
      case(PLASTICITY_DISLOUCLA_ID)
        call plastic_disloUCLA_results(phase_plasticityInstance(p),group) 
       
      case(PLASTICITY_NONLOCAL_ID)
        call plastic_nonlocal_results(phase_plasticityInstance(p),group) 
    end select
  
 enddo   
#endif


end subroutine constitutive_results


end module constitutive
