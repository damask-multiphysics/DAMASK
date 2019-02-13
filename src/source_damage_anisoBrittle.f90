!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic brittle damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_damage_anisoBrittle
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   source_damage_anisoBrittle_offset, &                                                                         !< which source is my current source mechanism?
   source_damage_anisoBrittle_instance                                                                          !< instance of source mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public  :: &
   source_damage_anisoBrittle_sizePostResult                                                                    !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public  :: &
   source_damage_anisoBrittle_output                                                                            !< name of each post result output
   
 integer(pInt),                       dimension(:,:),         allocatable,         private :: &
   source_damage_anisoBrittle_Ncleavage                                                                         !< number of cleavage systems per family

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 damage_drivingforce_ID
 end enum                                                


 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   real(pReal) :: &
     aTol, &
     sdot_0, &
     N
   real(pReal), dimension(:), allocatable :: &
     critDisp, &
     critLoad
   integer(pInt) :: &
     totalNcleavage
   integer(pInt), dimension(:), allocatable :: &
     Ncleavage
   integer(kind(undefined_ID)), allocatable, dimension(:) :: &
     outputID                                                                                       !< ID of each post result output
 end type tParameters

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)


 public :: &
   source_damage_anisoBrittle_init, &
   source_damage_anisoBrittle_dotState, &
   source_damage_anisobrittle_getRateAndItsTangent, &
   source_damage_anisoBrittle_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisoBrittle_init
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   pStringLen
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use IO, only: &
   IO_error
    use math, only: &
   math_expand
 use material, only: &
   material_allocateSourceState, &
   phase_source, &
   phase_Nsources, &
   phase_Noutput, &
   SOURCE_damage_anisoBrittle_label, &
   SOURCE_damage_anisoBrittle_ID, &
   material_phase, &
   sourceState
 use config, only: &
   config_phase, &
   material_Nphase, &
   MATERIAL_partPhase
 use lattice, only: &
   lattice_maxNcleavageFamily

 implicit none

 integer(pInt) :: Ninstance,phase,instance,source,sourceOffset
 integer(pInt) :: NofMyPhase,p   ,i
 integer(pInt),          dimension(0), parameter :: emptyIntArray    = [integer(pInt)::]
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
 integer(kind(undefined_ID)) :: &
   outputID

 character(len=pStringLen) :: &
   extmsg = ''
 character(len=65536), dimension(:), allocatable :: &
   outputs

 write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_DAMAGE_ANISOBRITTLE_LABEL//' init  -+>>>'
#include "compilation_info.f90"

 Ninstance = int(count(phase_source == SOURCE_damage_anisoBrittle_ID),pInt)
 if (Ninstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',Ninstance
 
 allocate(source_damage_anisoBrittle_offset(material_Nphase), source=0_pInt)
 allocate(source_damage_anisoBrittle_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   source_damage_anisoBrittle_instance(phase) = count(phase_source(:,1:phase) == source_damage_anisoBrittle_ID)
   do source = 1, phase_Nsources(phase)
     if (phase_source(source,phase) == source_damage_anisoBrittle_ID) &
       source_damage_anisoBrittle_offset(phase) = source
   enddo    
 enddo
 
 allocate(source_damage_anisoBrittle_sizePostResult(maxval(phase_Noutput),Ninstance), source=0_pInt)
 allocate(source_damage_anisoBrittle_output(maxval(phase_Noutput),Ninstance))
          source_damage_anisoBrittle_output = ''

 allocate(source_damage_anisoBrittle_Ncleavage(lattice_maxNcleavageFamily,Ninstance), source=0_pInt)

 allocate(param(Ninstance))
 
 do p=1, size(config_phase)
   if (all(phase_source(:,p) /= SOURCE_DAMAGE_ANISOBRITTLE_ID)) cycle
   associate(prm => param(source_damage_anisoBrittle_instance(p)), &
             config => config_phase(p))
             
   prm%aTol      = config%getFloat('anisobrittle_atol',defaultVal = 1.0e-3_pReal)

   prm%N         = config%getFloat('anisobrittle_ratesensitivity')
   prm%sdot_0    = config%getFloat('anisobrittle_sdot0')
   
   ! sanity checks
   if (prm%aTol      < 0.0_pReal) extmsg = trim(extmsg)//' anisobrittle_atol'
   
   if (prm%N        <= 0.0_pReal) extmsg = trim(extmsg)//' anisobrittle_ratesensitivity'
   if (prm%sdot_0   <= 0.0_pReal) extmsg = trim(extmsg)//' anisobrittle_sdot0'
   
   prm%Ncleavage = config%getInts('ncleavage',defaultVal=emptyIntArray)

   prm%critDisp = config%getFloats('anisobrittle_criticaldisplacement',requiredSize=size(prm%Ncleavage))
   prm%critLoad = config%getFloats('anisobrittle_criticalload',        requiredSize=size(prm%Ncleavage))

     ! expand: family => system
     prm%critDisp  = math_expand(prm%critDisp, prm%Ncleavage)
     prm%critLoad  = math_expand(prm%critLoad, prm%Ncleavage)
     
     if (any(prm%critLoad < 0.0_pReal))     extmsg = trim(extmsg)//' anisobrittle_criticalload'
    if (any(prm%critDisp < 0.0_pReal))     extmsg = trim(extmsg)//' anisobrittle_criticaldisplacement'  
!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
   if (extmsg /= '') &
     call IO_error(211_pInt,ext_msg=trim(extmsg)//'('//SOURCE_DAMAGE_ANISOBRITTLE_LABEL//')')

!--------------------------------------------------------------------------------------------------
!  output pararameters
   outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   do i=1_pInt, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
     
       case ('anisobrittle_drivingforce')
         source_damage_anisoBrittle_sizePostResult(i,source_damage_anisoBrittle_instance(p)) = 1_pInt
         source_damage_anisoBrittle_output(i,source_damage_anisoBrittle_instance(p)) = outputs(i)
         prm%outputID = [prm%outputID, damage_drivingforce_ID]

     end select

   enddo

   end associate
   
   phase = p
   NofMyPhase=count(material_phase==phase)
   instance = source_damage_anisoBrittle_instance(phase)
   sourceOffset = source_damage_anisoBrittle_offset(phase)


   call material_allocateSourceState(phase,sourceOffset,NofMyPhase,1_pInt)
   sourceState(phase)%p(sourceOffset)%sizePostResults = sum(source_damage_anisoBrittle_sizePostResult(:,instance))
   sourceState(phase)%p(sourceOffset)%aTolState=param(instance)%aTol


   source_damage_anisoBrittle_Ncleavage(1:size(param(instance)%Ncleavage),instance) = param(instance)%Ncleavage
 enddo

 
end subroutine source_damage_anisoBrittle_init

!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisoBrittle_dotState(S, ipc, ip, el)
 use math, only: &
   math_mul33xx33
 use material, only: &
   phaseAt, phasememberAt, &
   sourceState, &
   material_homog, &
   damage, &
   damageMapping
 use lattice, only: &
   lattice_Scleavage, &
   lattice_maxNcleavageFamily, &
   lattice_NcleavageSystem

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),  intent(in), dimension(3,3) :: &
   S
 integer(pInt) :: &
   phase, &
   constituent, &
   instance, &
   sourceOffset, &
   damageOffset, &
   homog, &
   f, i, index_myFamily, index
 real(pReal) :: &
   traction_d, traction_t, traction_n, traction_crit

 phase = phaseAt(ipc,ip,el)
 constituent = phasememberAt(ipc,ip,el)
 instance = source_damage_anisoBrittle_instance(phase)
 sourceOffset = source_damage_anisoBrittle_offset(phase)
 homog = material_homog(ip,el)
 damageOffset = damageMapping(homog)%p(ip,el)
 
 sourceState(phase)%p(sourceOffset)%dotState(1,constituent) = 0.0_pReal
 
 index = 1_pInt
 do f = 1_pInt,lattice_maxNcleavageFamily
   index_myFamily = sum(lattice_NcleavageSystem(1:f-1_pInt,phase))                                  ! at which index starts my family
   do i = 1_pInt,source_damage_anisoBrittle_Ncleavage(f,instance)                                   ! process each (active) cleavage system in family
     traction_d    = math_mul33xx33(S,lattice_Scleavage(1:3,1:3,1,index_myFamily+i,phase))
     traction_t    = math_mul33xx33(S,lattice_Scleavage(1:3,1:3,2,index_myFamily+i,phase))
     traction_n    = math_mul33xx33(S,lattice_Scleavage(1:3,1:3,3,index_myFamily+i,phase))
     
     traction_crit = param(instance)%critLoad(index)* &
                     damage(homog)%p(damageOffset)*damage(homog)%p(damageOffset)
     sourceState(phase)%p(sourceOffset)%dotState(1,constituent) = &
       sourceState(phase)%p(sourceOffset)%dotState(1,constituent) + &
       param(instance)%sdot_0* &
       ((max(0.0_pReal, abs(traction_d) - traction_crit)/traction_crit)**param(instance)%N + &
        (max(0.0_pReal, abs(traction_t) - traction_crit)/traction_crit)**param(instance)%N + &
        (max(0.0_pReal, abs(traction_n) - traction_crit)/traction_crit)**param(instance)%N)/ &
       param(instance)%critDisp(index)

   index = index + 1_pInt
   enddo
 enddo

end subroutine source_damage_anisoBrittle_dotState

!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisobrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)
 use material, only: &
   sourceState

 implicit none
 integer(pInt), intent(in) :: &
   phase, &
   constituent
 real(pReal),  intent(in) :: &
   phi
 real(pReal),  intent(out) :: &
   localphiDot, &
   dLocalphiDot_dPhi
 integer(pInt) :: &
   sourceOffset

 sourceOffset = source_damage_anisoBrittle_offset(phase)
 
 localphiDot = 1.0_pReal - &
               sourceState(phase)%p(sourceOffset)%state(1,constituent)*phi
 
 dLocalphiDot_dPhi = -sourceState(phase)%p(sourceOffset)%state(1,constituent)
 
end subroutine source_damage_anisobrittle_getRateAndItsTangent
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of local damage results
!--------------------------------------------------------------------------------------------------
function source_damage_anisoBrittle_postResults(phase, constituent)
 use material, only: &
   sourceState

 implicit none
 integer(pInt), intent(in) :: &
   phase, &
   constituent
 real(pReal), dimension(sum(source_damage_anisoBrittle_sizePostResult(:, &
                          source_damage_anisoBrittle_instance(phase)))) :: &
   source_damage_anisoBrittle_postResults

 integer(pInt) :: &
   instance, sourceOffset, o, c
   
 instance = source_damage_anisoBrittle_instance(phase)
 sourceOffset = source_damage_anisoBrittle_offset(phase)

 c = 0_pInt

 do o = 1_pInt,size(param(instance)%outputID)
    select case(param(instance)%outputID(o))
      case (damage_drivingforce_ID)
        source_damage_anisoBrittle_postResults(c+1_pInt) = &
          sourceState(phase)%p(sourceOffset)%state(1,constituent)
        c = c + 1_pInt

    end select
 enddo
end function source_damage_anisoBrittle_postResults

end module source_damage_anisoBrittle
