!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic ductile damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_damage_anisoDuctile
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   source_damage_anisoDuctile_offset, &                                                                         !< which source is my current damage mechanism?
   source_damage_anisoDuctile_instance                                                                          !< instance of damage source mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public  :: &
   source_damage_anisoDuctile_sizePostResult                                                                    !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public  :: &
   source_damage_anisoDuctile_output                                                                            !< name of each post result output
   

 integer(pInt),                       dimension(:,:),         allocatable,         private :: &
   source_damage_anisoDuctile_Nslip                                                                         !< number of slip systems per family
   
 enum, bind(c) 
   enumerator :: undefined_ID, &
                 damage_drivingforce_ID
 end enum 


 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   real(pReal) :: &
     aTol, &
     N
   real(pReal), dimension(:), allocatable :: &
     critPlasticStrain
   integer(pInt) :: &
     totalNslip
   integer(pInt), dimension(:), allocatable :: &
     Nslip
   integer(kind(undefined_ID)), allocatable, dimension(:) :: &
     outputID
 end type tParameters

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)


 public :: &
   source_damage_anisoDuctile_init, &
   source_damage_anisoDuctile_dotState, &
   source_damage_anisoDuctile_getRateAndItsTangent, &
   source_damage_anisoDuctile_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisoDuctile_init
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
   SOURCE_damage_anisoDuctile_label, &
   SOURCE_damage_anisoDuctile_ID, &
   material_phase, &  
   sourceState
 use config, only: &
   config_phase, &
   material_Nphase, &
   MATERIAL_partPhase
 use lattice, only: &
   lattice_maxNslipFamily
   
 implicit none

 integer(pInt) :: Ninstance,phase,instance,source,sourceOffset
 integer(pInt) :: NofMyPhase,p ,i

 integer(pInt),          dimension(0), parameter :: emptyIntArray    = [integer(pInt)::]
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
 integer(kind(undefined_ID)) :: &
   outputID

 character(len=pStringLen) :: &
   extmsg = ''
 character(len=65536), dimension(:), allocatable :: &
   outputs

 write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_DAMAGE_ANISODUCTILE_LABEL//' init  -+>>>'
#include "compilation_info.f90"

 Ninstance = int(count(phase_source == SOURCE_damage_anisoDuctile_ID),pInt)
 if (Ninstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',Ninstance
 
 allocate(source_damage_anisoDuctile_offset(material_Nphase), source=0_pInt)
 allocate(source_damage_anisoDuctile_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   source_damage_anisoDuctile_instance(phase) = count(phase_source(:,1:phase) == source_damage_anisoDuctile_ID)
   do source = 1, phase_Nsources(phase)
     if (phase_source(source,phase) == source_damage_anisoDuctile_ID) &
       source_damage_anisoDuctile_offset(phase) = source
   enddo    
 enddo
   
 allocate(source_damage_anisoDuctile_sizePostResult(maxval(phase_Noutput),Ninstance),source=0_pInt)
 allocate(source_damage_anisoDuctile_output(maxval(phase_Noutput),Ninstance))
          source_damage_anisoDuctile_output = ''

 allocate(source_damage_anisoDuctile_Nslip(lattice_maxNslipFamily,Ninstance),        source=0_pInt)

 allocate(param(Ninstance))
 
 do p=1, size(config_phase)
   if (all(phase_source(:,p) /= SOURCE_DAMAGE_ANISODUCTILE_ID)) cycle
   associate(prm => param(source_damage_anisoDuctile_instance(p)), &
             config => config_phase(p))
             
   prm%aTol   = config%getFloat('anisoductile_atol',defaultVal = 1.0e-3_pReal)

   prm%N      = config%getFloat('anisoductile_ratesensitivity')
   
   ! sanity checks
   if (prm%aTol                 < 0.0_pReal) extmsg = trim(extmsg)//' anisoductile_atol'
   
   if (prm%N                   <= 0.0_pReal) extmsg = trim(extmsg)//' anisoductile_ratesensitivity'
   
   prm%Nslip  = config%getInts('nslip',defaultVal=emptyIntArray)
   
   prm%critPlasticStrain = config%getFloats('anisoductile_criticalplasticstrain',requiredSize=size(prm%Nslip))

     ! expand: family => system
     prm%critPlasticStrain   = math_expand(prm%critPlasticStrain,  prm%Nslip)
     
          if (any(prm%critPlasticStrain < 0.0_pReal))     extmsg = trim(extmsg)//' anisoductile_criticalplasticstrain'
   
!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
   if (extmsg /= '') &
     call IO_error(211_pInt,ext_msg=trim(extmsg)//'('//SOURCE_DAMAGE_ANISODUCTILE_LABEL//')')

!--------------------------------------------------------------------------------------------------
!  output pararameters
   outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   do i=1_pInt, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
     
       case ('anisoductile_drivingforce')
         source_damage_anisoDuctile_sizePostResult(i,source_damage_anisoDuctile_instance(p)) = 1_pInt
         source_damage_anisoDuctile_output(i,source_damage_anisoDuctile_instance(p)) = outputs(i)
         prm%outputID = [prm%outputID, damage_drivingforce_ID]

     end select

   enddo

   end associate
   
   phase = p
   
   NofMyPhase=count(material_phase==phase)
   instance = source_damage_anisoDuctile_instance(phase)
   sourceOffset = source_damage_anisoDuctile_offset(phase)

   call material_allocateSourceState(phase,sourceOffset,NofMyPhase,1_pInt,1_pInt,0_pInt)
   sourceState(phase)%p(sourceOffset)%sizePostResults = sum(source_damage_anisoDuctile_sizePostResult(:,instance))
   sourceState(phase)%p(sourceOffset)%aTolState=param(instance)%aTol
   
   source_damage_anisoDuctile_Nslip(1:size(param(instance)%Nslip),instance) = param(instance)%Nslip
   
 enddo
  
end subroutine source_damage_anisoDuctile_init

!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisoDuctile_dotState(ipc, ip, el)
 use material, only: &
   phaseAt, phasememberAt, &
   plasticState, &
   sourceState, &
   material_homog, &
   damage, &
   damageMapping
 use lattice, only: &
   lattice_maxNslipFamily

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 integer(pInt) :: &
   phase, &
   constituent, &
   sourceOffset, &
   homog, damageOffset, &
   instance, &
   index, f, i

 phase = phaseAt(ipc,ip,el)
 constituent = phasememberAt(ipc,ip,el)
 instance = source_damage_anisoDuctile_instance(phase)
 sourceOffset = source_damage_anisoDuctile_offset(phase)
 homog = material_homog(ip,el)
 damageOffset = damageMapping(homog)%p(ip,el)

 index = 1_pInt
 sourceState(phase)%p(sourceOffset)%dotState(1,constituent) = 0.0_pReal
 do f = 1_pInt,lattice_maxNslipFamily
   do i = 1_pInt,source_damage_anisoDuctile_Nslip(f,instance)                                              ! process each (active) slip system in family
     sourceState(phase)%p(sourceOffset)%dotState(1,constituent) = &
       sourceState(phase)%p(sourceOffset)%dotState(1,constituent) + &
       plasticState(phase)%slipRate(index,constituent)/ &
       ((damage(homog)%p(damageOffset))**param(instance)%N)/param(instance)%critPlasticStrain(index) 
     
     index = index + 1_pInt
   enddo
 enddo
 
end subroutine source_damage_anisoDuctile_dotState

!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisoDuctile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)
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

 sourceOffset = source_damage_anisoDuctile_offset(phase)
 
 localphiDot = 1.0_pReal - &
               sourceState(phase)%p(sourceOffset)%state(1,constituent)* &
               phi
 
 dLocalphiDot_dPhi = -sourceState(phase)%p(sourceOffset)%state(1,constituent)
 
end subroutine source_damage_anisoDuctile_getRateAndItsTangent
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of local damage results
!--------------------------------------------------------------------------------------------------
function source_damage_anisoDuctile_postResults(phase, constituent)
 use material, only: &
   sourceState

 implicit none
 integer(pInt), intent(in) :: &
   phase, &
   constituent
 real(pReal), dimension(sum(source_damage_anisoDuctile_sizePostResult(:, &
                          source_damage_anisoDuctile_instance(phase)))) :: &
   source_damage_anisoDuctile_postResults

 integer(pInt) :: &
   instance, sourceOffset, o, c
   
 instance = source_damage_anisoDuctile_instance(phase)
 sourceOffset = source_damage_anisoDuctile_offset(phase)

 c = 0_pInt

 do o = 1_pInt,size(param(instance)%outputID)
    select case(param(instance)%outputID(o))
      case (damage_drivingforce_ID)
        source_damage_anisoDuctile_postResults(c+1_pInt) = &
          sourceState(phase)%p(sourceOffset)%state(1,constituent)
        c = c + 1_pInt

    end select
 enddo
end function source_damage_anisoDuctile_postResults

end module source_damage_anisoDuctile
