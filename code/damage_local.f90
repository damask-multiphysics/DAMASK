!--------------------------------------------------------------------------------------------------
! $Id$
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
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   damage_local_sizePostResults                                                           !< cumulative size of post results

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
 
 public :: &
   damage_local_init, &
   damage_local_updateState, &
   damage_local_postResults
 private :: &
   damage_local_getSourceAndItsTangent  

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine damage_local_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                         ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
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
   DAMAGE_local_label, &
   DAMAGE_local_ID, &
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
 integer(pInt) :: maxNinstance,mySize=0_pInt,homog,instance,o
 integer(pInt) :: sizeState
 integer(pInt) :: NofMyHomog   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  damage_'//DAMAGE_local_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 maxNinstance = int(count(damage_type == DAMAGE_local_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 allocate(damage_local_sizePostResults(maxNinstance),                               source=0_pInt)
 allocate(damage_local_sizePostResult (maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(damage_local_output         (maxval(homogenization_Noutput),maxNinstance))
          damage_local_output = ''
 allocate(damage_local_outputID       (maxval(homogenization_Noutput),maxNinstance),source=undefined_ID)
 allocate(damage_local_Noutput        (maxNinstance),                               source=0_pInt) 

 rewind(fileUnit)
 homog = 0_pInt
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
     homog = homog + 1_pInt                                                                         ! advance homog section counter
     cycle                                                                                          ! skip to next line
   endif

   if (homog > 0_pInt ) then; if (damage_type(homog) == DAMAGE_local_ID) then                       ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran

     instance = damage_typeInstance(homog)                                                          ! which instance of my damage is present homog
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('damage')
             damage_local_Noutput(instance) = damage_local_Noutput(instance) + 1_pInt
             damage_local_outputID(damage_local_Noutput(instance),instance) = damage_ID
             damage_local_output(damage_local_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select

     end select
   endif; endif
 enddo parsingFile

 initializeInstances: do homog = 1_pInt, size(damage_type)
   
   myhomog: if (damage_type(homog) == DAMAGE_local_ID) then
     NofMyHomog = count(material_homog == homog)
     instance = damage_typeInstance(homog)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,damage_local_Noutput(instance)
       select case(damage_local_outputID(o,instance))
         case(damage_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          damage_local_sizePostResult(o,instance) = mySize
          damage_local_sizePostResults(instance)  = damage_local_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop

! allocate state arrays
     sizeState = 1_pInt
     damageState(homog)%sizeState = sizeState
     damageState(homog)%sizePostResults = damage_local_sizePostResults(instance)
     allocate(damageState(homog)%state0   (sizeState,NofMyHomog))
     allocate(damageState(homog)%subState0(sizeState,NofMyHomog))
     allocate(damageState(homog)%state    (sizeState,NofMyHomog))

     nullify(damageMapping(homog)%p)
     damageMapping(homog)%p => mappingHomogenization(1,:,:)
     deallocate(damage(homog)%p)
     damage(homog)%p => damageState(homog)%state(1,:)
     
   endif myhomog
 enddo initializeInstances


end subroutine damage_local_init

!--------------------------------------------------------------------------------------------------
!> @brief  calculates local change in damage field   
!--------------------------------------------------------------------------------------------------
function damage_local_updateState(subdt, ip, el)
 use numerics, only: &
   err_damage_tolAbs, &
   err_damage_tolRel
 use material, only: &
   mappingHomogenization, &
   damageState
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   subdt
 logical,                    dimension(2)                             :: &
   damage_local_updateState
 integer(pInt) :: &
   homog, &
   offset
 real(pReal) :: &
   phi, phiDot, dPhiDot_dPhi  
 
 homog  = mappingHomogenization(2,ip,el)
 offset = mappingHomogenization(1,ip,el)
 phi = damageState(homog)%subState0(1,offset)
 call damage_local_getSourceAndItsTangent(phiDot, dPhiDot_dPhi, phi, ip, el)
 phi = phi + subdt*phiDot
 
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
   do source = 1, phase_Nsources(phase)
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
 
end subroutine damage_local_getSourceAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief return array of damage results
!--------------------------------------------------------------------------------------------------
function damage_local_postResults(ip,el)
 use material, only: &
   mappingHomogenization, &
   damage_typeInstance, &
   damageMapping, &
   damage

 implicit none
 integer(pInt),              intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(damage_local_sizePostResults(damage_typeInstance(mappingHomogenization(2,ip,el)))) :: &
   damage_local_postResults

 integer(pInt) :: &
   instance, homog, offset, o, c
   
 homog     = mappingHomogenization(2,ip,el)
 offset    = damageMapping(homog)%p(ip,el)
 instance  = damage_typeInstance(homog)

 c = 0_pInt
 damage_local_postResults = 0.0_pReal

 do o = 1_pInt,damage_local_Noutput(instance)
    select case(damage_local_outputID(o,instance))
 
      case (damage_ID)
        damage_local_postResults(c+1_pInt) = damage(homog)%p(offset)
        c = c + 1
    end select
 enddo
end function damage_local_postResults

end module damage_local
