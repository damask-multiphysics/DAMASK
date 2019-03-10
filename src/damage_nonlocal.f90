!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for non-locally evolving damage field
!> @details to be done
!--------------------------------------------------------------------------------------------------
module damage_nonlocal
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   damage_nonlocal_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   damage_nonlocal_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   damage_nonlocal_Noutput                                                                   !< number of outputs per instance of this damage 

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 damage_ID
 end enum

 type, private :: tParameters
   integer(kind(undefined_ID)),         dimension(:),   allocatable   :: &
     outputID
 end type tParameters
 
 type(tparameters),          dimension(:), allocatable, private :: &
   param

 public :: &
   damage_nonlocal_init, &
   damage_nonlocal_getSourceAndItsTangent, &
   damage_nonlocal_getDiffusion33, &
   damage_nonlocal_getMobility, &
   damage_nonlocal_putNonLocalDamage, &
   damage_nonlocal_postResults

contains

!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine damage_nonlocal_init
 use material, only: &
   damage_type, &
   damage_typeInstance, &
   homogenization_Noutput, &
   DAMAGE_nonlocal_label, &
   DAMAGE_nonlocal_ID, &
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

 write(6,'(/,a)')   ' <<<+-  damage_'//DAMAGE_nonlocal_label//' init  -+>>>'
 
 maxNinstance = int(count(damage_type == DAMAGE_nonlocal_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 allocate(damage_nonlocal_sizePostResult (maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(damage_nonlocal_output         (maxval(homogenization_Noutput),maxNinstance))
          damage_nonlocal_output = ''
 allocate(damage_nonlocal_Noutput        (maxNinstance),                               source=0_pInt) 

 allocate(param(maxNinstance))
  
 do h = 1, size(damage_type)
   if (damage_type(h) /= DAMAGE_NONLOCAL_ID) cycle
   associate(prm => param(damage_typeInstance(h)), &
             config => config_homogenization(h))
             
   instance = damage_typeInstance(h)
   outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   
   do i=1, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
     
           case ('damage')
           damage_nonlocal_output(i,damage_typeInstance(h)) = outputs(i)
             damage_nonlocal_Noutput(instance) = damage_nonlocal_Noutput(instance) + 1
            damage_nonlocal_sizePostResult(i,damage_typeInstance(h)) = 1
       prm%outputID = [prm%outputID , damage_ID]
          end select
     
   enddo

   homog = h

     NofMyHomog = count(material_homogenizationAt == homog)
     instance = damage_typeInstance(homog)


! allocate state arrays
     sizeState = 1_pInt
     damageState(homog)%sizeState = sizeState
     damageState(homog)%sizePostResults = sum(damage_nonlocal_sizePostResult(:,instance))
     allocate(damageState(homog)%state0   (sizeState,NofMyHomog), source=damage_initialPhi(homog))
     allocate(damageState(homog)%subState0(sizeState,NofMyHomog), source=damage_initialPhi(homog))
     allocate(damageState(homog)%state    (sizeState,NofMyHomog), source=damage_initialPhi(homog))

     nullify(damageMapping(homog)%p)
     damageMapping(homog)%p => mappingHomogenization(1,:,:)
     deallocate(damage(homog)%p)
     damage(homog)%p => damageState(homog)%state(1,:)
     
   end associate
 enddo
end subroutine damage_nonlocal_init

!--------------------------------------------------------------------------------------------------
!> @brief  calculates homogenized damage driving forces  
!--------------------------------------------------------------------------------------------------
subroutine damage_nonlocal_getSourceAndItsTangent(phiDot, dPhiDot_dPhi, phi, ip, el)
 use material, only: &
   homogenization_Ngrains, &
   mappingHomogenization, &
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
 do grain = 1, homogenization_Ngrains(mappingHomogenization(2,ip,el))
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
 
 phiDot = phiDot/real(homogenization_Ngrains(mappingHomogenization(2,ip,el)),pReal)
 dPhiDot_dPhi = dPhiDot_dPhi/real(homogenization_Ngrains(mappingHomogenization(2,ip,el)),pReal)
 
end subroutine damage_nonlocal_getSourceAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized non local damage diffusion tensor in reference configuration
!--------------------------------------------------------------------------------------------------
function damage_nonlocal_getDiffusion33(ip,el)
 use numerics, only: &
   charLength
 use lattice, only: &
   lattice_DamageDiffusion33
 use material, only: &
   homogenization_Ngrains, &
   material_phase, &
   mappingHomogenization
 use crystallite, only: &
   crystallite_push33ToRef

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   damage_nonlocal_getDiffusion33
 integer(pInt) :: &
   homog, &
   grain
   
 homog  = mappingHomogenization(2,ip,el)
 damage_nonlocal_getDiffusion33 = 0.0_pReal  
 do grain = 1, homogenization_Ngrains(homog)
   damage_nonlocal_getDiffusion33 = damage_nonlocal_getDiffusion33 + &
     crystallite_push33ToRef(grain,ip,el,lattice_DamageDiffusion33(1:3,1:3,material_phase(grain,ip,el)))
 enddo

 damage_nonlocal_getDiffusion33 = &
   charLength**2_pInt*damage_nonlocal_getDiffusion33/real(homogenization_Ngrains(homog),pReal)
 
end function damage_nonlocal_getDiffusion33
 
!--------------------------------------------------------------------------------------------------
!> @brief Returns homogenized nonlocal damage mobility 
!--------------------------------------------------------------------------------------------------
real(pReal) function damage_nonlocal_getMobility(ip,el)
 use mesh, only: &
   mesh_element
 use lattice, only: &
   lattice_damageMobility
 use material, only: &
   material_phase, &
   homogenization_Ngrains

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 integer(pInt) :: &
   ipc
 
 damage_nonlocal_getMobility = 0.0_pReal
                                                
 do ipc = 1, homogenization_Ngrains(mesh_element(3,el))
   damage_nonlocal_getMobility = damage_nonlocal_getMobility + lattice_DamageMobility(material_phase(ipc,ip,el))
 enddo

 damage_nonlocal_getMobility = damage_nonlocal_getMobility/&
                               real(homogenization_Ngrains(mesh_element(3,el)),pReal)

end function damage_nonlocal_getMobility

!--------------------------------------------------------------------------------------------------
!> @brief updated nonlocal damage field with solution from damage phase field PDE
!--------------------------------------------------------------------------------------------------
subroutine damage_nonlocal_putNonLocalDamage(phi,ip,el)
 use material, only: &
   material_homogenizationAt, &
   damageMapping, &
   damage

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   phi
 integer(pInt) :: &
   homog, &
   offset
 
 homog  = material_homogenizationAt(el)
 offset = damageMapping(homog)%p(ip,el)
 damage(homog)%p(offset) = phi

end subroutine damage_nonlocal_putNonLocalDamage
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of damage results
!--------------------------------------------------------------------------------------------------
function damage_nonlocal_postResults(ip,el)
 use material, only: &
   mappingHomogenization, &
   damage_typeInstance, &
   damageMapping, &
   damage

 implicit none
 integer(pInt),              intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(sum(damage_nonlocal_sizePostResult(:,damage_typeInstance(mappingHomogenization(2,ip,el))))) :: &
   damage_nonlocal_postResults

 integer(pInt) :: &
   instance, homog, offset, o, c
   
 homog     = mappingHomogenization(2,ip,el)
 offset    = damageMapping(homog)%p(ip,el)
 instance  = damage_typeInstance(homog)
 associate(prm => param(instance))
 c = 0_pInt

 outputsLoop: do o = 1_pInt,size(prm%outputID)
   select case(prm%outputID(o))
 
      case (damage_ID)
        damage_nonlocal_postResults(c+1_pInt) = damage(homog)%p(offset)
        c = c + 1
    end select
 enddo outputsLoop

 end associate
end function damage_nonlocal_postResults

end module damage_nonlocal
