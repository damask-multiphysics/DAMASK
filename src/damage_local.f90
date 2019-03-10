!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for locally evolving damage field
!--------------------------------------------------------------------------------------------------
module damage_local
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   damage_local_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   damage_local_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   damage_local_Noutput                                                                   !< number of outputs per instance of this damage 

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 damage_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   damage_local_outputID                                                                  !< ID of each post result output
   
 type, private :: tParameters
   integer(kind(undefined_ID)),         dimension(:),   allocatable   :: &
     outputID
 end type tParameters
 
 type(tparameters),          dimension(:), allocatable, private :: &
   param
   
 public :: &
   damage_local_init, &
   damage_local_updateState, &
   damage_local_postResults
 private :: &
   damage_local_getSourceAndItsTangent  

contains

!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine damage_local_init
 use material, only: &
   damage_type, &
   damage_typeInstance, &
   homogenization_Noutput, &
   DAMAGE_local_label, &
   DAMAGE_local_ID, &
   material_homogenizationAt, & 
   mappingHomogenization, & 
   damageState, &
   damageMapping, &
   damage, &
   damage_initialPhi
 use config, only: &
   config_homogenization
 
 implicit none

 integer(pInt) :: maxNinstance,homog,instance,o,i
 integer(pInt) :: sizeState
 integer(pInt) :: NofMyHomog, h
  integer(kind(undefined_ID)) :: &
   outputID
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
  character(len=65536), dimension(:), allocatable :: &
   outputs
 
 write(6,'(/,a)')   ' <<<+-  damage_'//DAMAGE_local_label//' init  -+>>>'

 maxNinstance = int(count(damage_type == DAMAGE_local_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 allocate(damage_local_sizePostResult (maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(damage_local_output         (maxval(homogenization_Noutput),maxNinstance))
          damage_local_output = ''
 allocate(damage_local_outputID       (maxval(homogenization_Noutput),maxNinstance),source=undefined_ID)
 allocate(damage_local_Noutput        (maxNinstance),                               source=0_pInt) 
 
 allocate(param(maxNinstance))
  
 do h = 1, size(damage_type)
   if (damage_type(h) /= DAMAGE_LOCAL_ID) cycle
   associate(prm => param(damage_typeInstance(h)), &
             config => config_homogenization(h))
             

   outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   
   do i=1, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
     
           case ('damage')
           damage_local_output(i,damage_typeInstance(h)) = outputs(i)
             damage_local_Noutput(instance) = damage_local_Noutput(instance) + 1
            damage_local_sizePostResult(i,damage_typeInstance(h)) = 1
       prm%outputID = [prm%outputID , damage_ID]
          end select
     
   enddo


   homog = h

     NofMyHomog = count(material_homogenizationAt == homog)
     instance = damage_typeInstance(homog)


! allocate state arrays
     sizeState = 1_pInt
     damageState(homog)%sizeState = sizeState
     damageState(homog)%sizePostResults = sum(damage_local_sizePostResult(:,instance))
     allocate(damageState(homog)%state0   (sizeState,NofMyHomog), source=damage_initialPhi(homog))
     allocate(damageState(homog)%subState0(sizeState,NofMyHomog), source=damage_initialPhi(homog))
     allocate(damageState(homog)%state    (sizeState,NofMyHomog), source=damage_initialPhi(homog))

     nullify(damageMapping(homog)%p)
     damageMapping(homog)%p => mappingHomogenization(1,:,:)
     deallocate(damage(homog)%p)
     damage(homog)%p => damageState(homog)%state(1,:)
     
   end associate
 enddo


end subroutine damage_local_init

!--------------------------------------------------------------------------------------------------
!> @brief  calculates local change in damage field   
!--------------------------------------------------------------------------------------------------
function damage_local_updateState(subdt, ip, el)
 use numerics, only: &
   residualStiffness, &
   err_damage_tolAbs, &
   err_damage_tolRel
 use material, only: &
   material_homogenizationAt, &
   mappingHomogenization, &
   damageState
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   subdt
 logical,                    dimension(2)  :: &
   damage_local_updateState
 integer(pInt) :: &
   homog, &
   offset
 real(pReal) :: &
   phi, phiDot, dPhiDot_dPhi  
 
 homog  = material_homogenizationAt(el)
 offset = mappingHomogenization(1,ip,el)
 phi = damageState(homog)%subState0(1,offset)
 call damage_local_getSourceAndItsTangent(phiDot, dPhiDot_dPhi, phi, ip, el)
 phi = max(residualStiffness,min(1.0_pReal,phi + subdt*phiDot))
 
 damage_local_updateState = [     abs(phi - damageState(homog)%state(1,offset)) &
                               <= err_damage_tolAbs &
                             .or. abs(phi - damageState(homog)%state(1,offset)) &
                               <= err_damage_tolRel*abs(damageState(homog)%state(1,offset)), &
                             .true.]

 damageState(homog)%state(1,offset) = phi  

end function damage_local_updateState

!--------------------------------------------------------------------------------------------------
!> @brief  calculates homogenized local damage driving forces  
!--------------------------------------------------------------------------------------------------
subroutine damage_local_getSourceAndItsTangent(phiDot, dPhiDot_dPhi, phi, ip, el)
 use material, only: &
   homogenization_Ngrains, &
   material_homogenizationAt, &
   phaseAt, &
   phasememberAt, &
   phase_source, &
   phase_Nsources, &
   SOURCE_damage_isoBrittle_ID, &
   SOURCE_damage_isoDuctile_ID, &
   SOURCE_damage_anisoBrittle_ID, &
   SOURCE_damage_anisoDuctile_ID
 use source_damage_isoBrittle, only: &
   source_damage_isobrittle_getRateAndItsTangent
 use source_damage_isoDuctile, only: &
   source_damage_isoductile_getRateAndItsTangent
 use source_damage_anisoBrittle, only: &
   source_damage_anisobrittle_getRateAndItsTangent
 use source_damage_anisoDuctile, only: &
   source_damage_anisoductile_getRateAndItsTangent
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   phi
 integer(pInt) :: &
   phase, &
   grain, &
   source, &
   constituent
 real(pReal) :: &
   phiDot, dPhiDot_dPhi, localphiDot, dLocalphiDot_dPhi  

 phiDot = 0.0_pReal
 dPhiDot_dPhi = 0.0_pReal
 do grain = 1, homogenization_Ngrains(material_homogenizationAt(el))
   phase = phaseAt(grain,ip,el)
   constituent = phasememberAt(grain,ip,el)
   do source = 1, phase_Nsources(phase)
     select case(phase_source(source,phase))                                                   
       case (SOURCE_damage_isoBrittle_ID)
        call source_damage_isobrittle_getRateAndItsTangent  (localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

       case (SOURCE_damage_isoDuctile_ID)
        call source_damage_isoductile_getRateAndItsTangent  (localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

       case (SOURCE_damage_anisoBrittle_ID)
        call source_damage_anisobrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

       case (SOURCE_damage_anisoDuctile_ID)
        call source_damage_anisoductile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

       case default
        localphiDot = 0.0_pReal
        dLocalphiDot_dPhi = 0.0_pReal

     end select
     phiDot = phiDot + localphiDot
     dPhiDot_dPhi = dPhiDot_dPhi + dLocalphiDot_dPhi
   enddo  
 enddo
 
 phiDot = phiDot/real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)
 dPhiDot_dPhi = dPhiDot_dPhi/real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)
 
end subroutine damage_local_getSourceAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief return array of damage results
!--------------------------------------------------------------------------------------------------
function damage_local_postResults(ip,el)
 use material, only: &
   material_homogenizationAt, &
   damage_typeInstance, &
   damageMapping, &
   damage

 implicit none
 integer(pInt),              intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(sum(damage_local_sizePostResult(:,damage_typeInstance(material_homogenizationAt(el))))) :: &
   damage_local_postResults

 integer(pInt) :: &
   instance, homog, offset, o, c
   
 homog     = material_homogenizationAt(el)
 offset    = damageMapping(homog)%p(ip,el)
 instance  = damage_typeInstance(homog)
 associate(prm => param(instance))
 c = 0_pInt

 outputsLoop: do o = 1_pInt,size(prm%outputID)
   select case(prm%outputID(o))
 
      case (damage_ID)
        damage_local_postResults(c+1_pInt) = damage(homog)%p(offset)
        c = c + 1
    end select
 enddo outputsLoop
 
 end associate
end function damage_local_postResults

end module damage_local
