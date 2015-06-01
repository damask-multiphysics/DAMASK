!--------------------------------------------------------------------------------------------------
! $Id$
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
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   damage_nonlocal_sizePostResults                                                           !< cumulative size of post results

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
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   damage_nonlocal_outputID                                                                  !< ID of each post result output


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
subroutine damage_nonlocal_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use IO, only: &
   IO_read, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue, &
   IO_intValue, &
   IO_warning, &
   IO_error, &
   IO_timeStamp, &
   IO_EOF
 use material, only: &
   damage_type, &
   damage_typeInstance, &
   homogenization_Noutput, &
   DAMAGE_nonlocal_label, &
   DAMAGE_nonlocal_ID, &
   material_homog, & 
   mappingHomogenization, & 
   damageState, &
   damageMapping, &
   damage, &
   material_partHomogenization
 use numerics,only: &
   worldrank

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), parameter :: MAXNCHUNKS = 7_pInt
 integer(pInt), dimension(1+2*MAXNCHUNKS) :: positions
 integer(pInt) :: maxNinstance,mySize=0_pInt,section,instance,o
 integer(pInt) :: sizeState
 integer(pInt) :: NofMyHomog   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  damage_'//DAMAGE_nonlocal_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(damage_type == DAMAGE_nonlocal_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 allocate(damage_nonlocal_sizePostResults(maxNinstance),                               source=0_pInt)
 allocate(damage_nonlocal_sizePostResult (maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(damage_nonlocal_output         (maxval(homogenization_Noutput),maxNinstance))
          damage_nonlocal_output = ''
 allocate(damage_nonlocal_outputID       (maxval(homogenization_Noutput),maxNinstance),source=undefined_ID)
 allocate(damage_nonlocal_Noutput        (maxNinstance),                               source=0_pInt) 

 rewind(fileUnit)
 section = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partHomogenization)! wind forward to <homogenization>
   line = IO_read(fileUnit)
 enddo
 
 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of homog part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif   
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next homog section
     section = section + 1_pInt                                                                     ! advance homog section counter
     cycle                                                                                          ! skip to next line
   endif

   if (section > 0_pInt ) then; if (damage_type(section) == DAMAGE_nonlocal_ID) then             ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran

     instance = damage_typeInstance(section)                                                       ! which instance of my damage is present homog
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('damage')
             damage_nonlocal_Noutput(instance) = damage_nonlocal_Noutput(instance) + 1_pInt
             damage_nonlocal_outputID(damage_nonlocal_Noutput(instance),instance) = damage_ID
             damage_nonlocal_output(damage_nonlocal_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select

     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do section = 1_pInt, size(damage_type)
   if (damage_type(section) == DAMAGE_nonlocal_ID) then
     NofMyHomog=count(material_homog==section)
     instance = damage_typeInstance(section)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,damage_nonlocal_Noutput(instance)
       select case(damage_nonlocal_outputID(o,instance))
         case(damage_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          damage_nonlocal_sizePostResult(o,instance) = mySize
          damage_nonlocal_sizePostResults(instance)  = damage_nonlocal_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop

! allocate state arrays
     sizeState = 0_pInt
     damageState(section)%sizeState = sizeState
     damageState(section)%sizePostResults = damage_nonlocal_sizePostResults(instance)
     allocate(damageState(section)%state0   (sizeState,NofMyHomog))
     allocate(damageState(section)%subState0(sizeState,NofMyHomog))
     allocate(damageState(section)%state    (sizeState,NofMyHomog))

     nullify(damageMapping(section)%p)
     damageMapping(section)%p => mappingHomogenization(1,:,:)
     deallocate(damage(section)%p)
     allocate(damage(section)%p(NofMyHomog), source=1.0_pReal)
     
   endif
 
 enddo initializeInstances
end subroutine damage_nonlocal_init

!--------------------------------------------------------------------------------------------------
!> @brief  calculates homogenized damage driving forces  
!--------------------------------------------------------------------------------------------------
subroutine damage_nonlocal_getSourceAndItsTangent(phiDot, dPhiDot_dPhi, phi, ip, el)
 use material, only: &
   homogenization_Ngrains, &
   mappingHomogenization, &
   mappingConstitutive, &
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
   source
 real(pReal) :: &
   phiDot, dPhiDot_dPhi, localphiDot, dLocalphiDot_dPhi  

 phiDot = 0.0_pReal
 dPhiDot_dPhi = 0.0_pReal
 do grain = 1, homogenization_Ngrains(mappingHomogenization(2,ip,el))
   phase = mappingConstitutive(2,grain,ip,el)
   do source = 1_pInt, phase_Nsources(phase)
     select case(phase_source(source,phase))                                                   
       case (SOURCE_damage_isoBrittle_ID)
        call source_damage_isobrittle_getRateAndItsTangent  (localphiDot, dLocalphiDot_dPhi, phi, grain, ip,  el)

       case (SOURCE_damage_isoDuctile_ID)
        call source_damage_isoductile_getRateAndItsTangent  (localphiDot, dLocalphiDot_dPhi, phi, grain, ip,  el)

       case (SOURCE_damage_anisoBrittle_ID)
        call source_damage_anisobrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, grain, ip,  el)

       case (SOURCE_damage_anisoDuctile_ID)
        call source_damage_anisoductile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, grain, ip,  el)

       case default
        localphiDot = 0.0_pReal
        dLocalphiDot_dPhi = 0.0_pReal

     end select
     phiDot = phiDot + localphiDot
     dPhiDot_dPhi = dPhiDot_dPhi + dLocalphiDot_dPhi
   enddo  
 enddo
 
 phiDot = phiDot/homogenization_Ngrains(mappingHomogenization(2,ip,el))
 dPhiDot_dPhi = dPhiDot_dPhi/homogenization_Ngrains(mappingHomogenization(2,ip,el))
 
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
   charLength*charLength* &
   damage_nonlocal_getDiffusion33/ &
   homogenization_Ngrains(homog)
 
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

 damage_nonlocal_getMobility = damage_nonlocal_getMobility /homogenization_Ngrains(mesh_element(3,el))

end function damage_nonlocal_getMobility

!--------------------------------------------------------------------------------------------------
!> @brief updated nonlocal damage field with solution from damage phase field PDE
!--------------------------------------------------------------------------------------------------
subroutine damage_nonlocal_putNonLocalDamage(phi,ip,el)
 use material, only: &
   material_homog, &
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
 
 homog  = material_homog(ip,el)
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
   damage

 implicit none
 integer(pInt),              intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(damage_nonlocal_sizePostResults(damage_typeInstance(mappingHomogenization(2,ip,el)))) :: &
   damage_nonlocal_postResults

 integer(pInt) :: &
   instance, homog, offset, o, c
   
 homog     = mappingHomogenization(2,ip,el)
 offset    = mappingHomogenization(1,ip,el)
 instance  = damage_typeInstance(homog)

 c = 0_pInt
 damage_nonlocal_postResults = 0.0_pReal

 do o = 1_pInt,damage_nonlocal_Noutput(instance)
    select case(damage_nonlocal_outputID(o,instance))
 
      case (damage_ID)
        damage_nonlocal_postResults(c+1_pInt) = damage(homog)%p(offset)
        c = c + 1
    end select
 enddo
end function damage_nonlocal_postResults

end module damage_nonlocal
