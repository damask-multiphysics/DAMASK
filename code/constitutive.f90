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
   constitutive_maxSizeDotState

 public :: & 
   constitutive_init, &
   constitutive_homogenizedC, &
   constitutive_microstructure, &
   constitutive_LpAndItsTangent, &
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
   plasticState, &
   mappingConstitutive
 

 use constitutive_none
 use constitutive_j2
 use constitutive_phenopowerlaw
 use constitutive_dislotwin
 use constitutive_dislokmc
 use constitutive_titanmod
 use constitutive_nonlocal
 implicit none
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt) :: &
  e, &                                                                                        !< maximum number of elements
  phase, &
  instance

 integer(pInt), dimension(:,:), pointer :: thisSize
 character(len=64), dimension(:,:), pointer :: thisOutput
 character(len=32) :: outputName                                                                    !< name of output, intermediate fix until HDF5 output is ready
 logical :: knownPlasticity, nonlocalConstitutionPresent
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
 enddo
 close(FILEUNIT)
 
 constitutive_maxSizeDotState = 0_pInt
 constitutive_maxSizePostResults = 0_pInt

 PhaseLoop2:do phase = 1_pInt,material_Nphase
   plasticState(phase)%partionedState0 = plasticState(phase)%State0
   plasticState(phase)%State = plasticState(phase)%State0
   constitutive_maxSizeDotState = max(constitutive_maxSizeDotState, plasticState(phase)%sizeDotState)
   constitutive_maxSizePostResults = max(constitutive_maxSizePostResults, plasticState(phase)%sizePostResults)
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
! write out state size file
 call IO_write_jobIntFile(777,'sizeStateConst', size(constitutive_sizeState))
 write (777,rec=1) constitutive_sizeState
 close(777)

!--------------------------------------------------------------------------------------------------
! report
 constitutive_maxSizeState       = maxval(constitutive_sizeState)
 constitutive_maxSizeDotState    = maxval(constitutive_sizeDotState)
 constitutive_maxSizePostResults = maxval(constitutive_sizePostResults)
 
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
subroutine constitutive_microstructure(temperature, Fe, Fp, ipc, ip, el)
 use prec, only: &
   pReal 
 use material, only: &
   phase_plasticity, &
   material_phase, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOKMC_ID, &
   PLASTICITY_TITANMOD_ID, &
   PLASTICITY_NONLOCAL_ID, &
   plasticState, &
   mappingConstitutive

 use constitutive_titanmod, only: &
   constitutive_titanmod_microstructure
 use constitutive_nonlocal, only: &
   constitutive_nonlocal_microstructure
 use constitutive_dislotwin, only: &
   constitutive_dislotwin_microstructure
 use constitutive_dislokmc, only: &
   constitutive_dislokmc_microstructure

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   temperature
 real(pReal),   intent(in), dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fp                                                                                               !< plastic deformation gradient

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
 use constitutive_damage, only: &
   constitutive_damageValue  
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
 
 select case (phase_plasticity(material_phase(ipc,ip,el)))
 
   case (PLASTICITY_NONE_ID)
     Lp = 0.0_pReal
     dLp_dTstar = 0.0_pReal
   case (PLASTICITY_J2_ID)
     call constitutive_j2_LpAndItsTangent           (Lp,dLp_dTstar,Tstar_v,constitutive_damageValue(ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID)
     call constitutive_phenopowerlaw_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,constitutive_damageValue(ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID)
     call constitutive_nonlocal_LpAndItsTangent     (Lp,dLp_dTstar,Tstar_v,temperature,    ip,el)
   case (PLASTICITY_DISLOTWIN_ID)
     call constitutive_dislotwin_LpAndItsTangent    (Lp,dLp_dTstar,Tstar_v,temperature,ipc,ip,el)
   case (PLASTICITY_DISLOKMC_ID)
     call constitutive_dislokmc_LpAndItsTangent    (Lp,dLp_dTstar,Tstar_v,temperature,ipc,ip,el)
   case (PLASTICITY_TITANMOD_ID)
     call constitutive_titanmod_LpAndItsTangent     (Lp,dLp_dTstar,Tstar_v,temperature,ipc,ip,el)

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
 use constitutive_damage, only: &
   constitutive_damageValue  
 use constitutive_thermal, only: &
   constitutive_temperature  

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
 real(pReal)   :: damage, pressure
 real(pReal), dimension(3,3,3,3) :: C

 C = math_Mandel66to3333(constitutive_homogenizedC(ipc,ip,el))
 T = math_mul3333xx33(C,0.5_pReal*(math_mul33x33(math_transpose33(Fe),Fe)-math_I3) - &
                        lattice_thermalExpansion33(1:3,1:3,mappingConstitutive(2,ipc,ip,el))* &
                        (constitutive_temperature(ipc,ip,el) - &
                         lattice_referenceTemperature(mappingConstitutive(2,ipc,ip,el))))
 
 pressure = math_trace33(T)/3.0_pReal
 damage = constitutive_damageValue(ipc,ip,el)
 T = damage*T
 if (pressure < 0.0_pReal) T = T + (1.0_pReal - damage)*pressure*math_I3

 dT_dFe = 0.0_pReal
 forall (i=1_pInt:3_pInt, j=1_pInt:3_pInt, k=1_pInt:3_pInt, l=1_pInt:3_pInt) &
   dT_dFe(i,j,k,l) = damage*sum(C(i,j,l,1:3)*Fe(k,1:3))                                            ! dT*_ij/dFe_kl
 
 if (pressure < 0.0_pReal) then
   do i=1_pInt, 3_pInt; do k=1_pInt,3_pInt; do l=1_pInt,3_pInt; do j=1_pInt,3_pInt
     dT_dFe(i,i,k,l) = dT_dFe(i,i,k,l) + &
                       (1.0_pReal - damage)*math_trace33(C(1:3,1:3,l,j))*Fe(k,j)/3.0_pReal
   enddo; enddo; enddo; enddo
 endif  
 
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
   plasticState, &
   mappingConstitutive, &  
   material_phase, &
   homogenization_maxNgrains, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_J2_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOKMC_ID, &
   PLASTICITY_TITANMOD_ID, &
   PLASTICITY_NONLOCAL_ID
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
#ifdef NEWSTATE   
   plasticState, &
   mappingConstitutive, &
#endif
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
   mappingConstitutive, &
   phase_plasticity, &
   material_phase, &
   homogenization_maxNgrains, &
   PLASTICITY_NONE_ID, &
   PLASTICITY_J2_ID, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   PLASTICITY_DISLOTWIN_ID, &
   PLASTICITY_DISLOKMC_ID, &
   PLASTICITY_TITANMOD_ID, &
   PLASTICITY_NONLOCAL_ID
 use constitutive_damage, only: &
   constitutive_damageValue  
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
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%sizePostResults) :: &
   constitutive_postResults
 real(pReal),  intent(in) :: &
   temperature
 real(pReal),  intent(in), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   FeArray                                                                                          !< elastic deformation gradient
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)

 constitutive_postResults = 0.0_pReal
 
 select case (phase_plasticity(material_phase(ipc,ip,el)))
   case (PLASTICITY_TITANMOD_ID)
     constitutive_postResults = constitutive_titanmod_postResults             (ipc,ip,el)
   case (PLASTICITY_J2_ID)
     constitutive_postResults= constitutive_j2_postResults            (Tstar_v,constitutive_damageValue(ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_PHENOPOWERLAW_ID)
     constitutive_postResults = constitutive_phenopowerlaw_postResults(Tstar_v,constitutive_damageValue(ipc,ip,el),ipc,ip,el)
   case (PLASTICITY_DISLOTWIN_ID)
     constitutive_postResults = constitutive_dislotwin_postResults(Tstar_v,Temperature,ipc,ip,el)
   case (PLASTICITY_DISLOKMC_ID)
     constitutive_postResults = constitutive_dislokmc_postResults(Tstar_v,Temperature,ipc,ip,el)
   case (PLASTICITY_NONLOCAL_ID)
     constitutive_postResults = constitutive_nonlocal_postResults (Tstar_v,FeArray,        ip,el)
 end select
  
end function constitutive_postResults


end module constitutive
