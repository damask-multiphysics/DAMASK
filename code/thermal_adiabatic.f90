!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for adiabatic temperature evolution
!> @details to be done
!--------------------------------------------------------------------------------------------------
module thermal_adiabatic
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   thermal_adiabatic_sizePostResults                                                           !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   thermal_adiabatic_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   thermal_adiabatic_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   thermal_adiabatic_Noutput                                                                   !< number of outputs per instance of this thermal model 

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 temperature_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   thermal_adiabatic_outputID                                                                  !< ID of each post result output


 public :: &
   thermal_adiabatic_init, &
   thermal_adiabatic_updateState, &
   thermal_adiabatic_getSourceAndItsTangent, &
   thermal_adiabatic_getSpecificHeat, &
   thermal_adiabatic_getMassDensity, &
   thermal_adiabatic_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine thermal_adiabatic_init(fileUnit)
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
   thermal_type, &
   thermal_typeInstance, &
   homogenization_Noutput, &
   THERMAL_ADIABATIC_label, &
   THERMAL_adiabatic_ID, &
   material_homog, & 
   mappingHomogenization, & 
   thermalState, &
   thermalMapping, &
   thermal_initialT, &
   temperature, &
   temperatureRate, &
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
   write(6,'(/,a)')   ' <<<+-  thermal_'//THERMAL_ADIABATIC_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(thermal_type == THERMAL_adiabatic_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 allocate(thermal_adiabatic_sizePostResults(maxNinstance),                               source=0_pInt)
 allocate(thermal_adiabatic_sizePostResult (maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(thermal_adiabatic_output         (maxval(homogenization_Noutput),maxNinstance))
          thermal_adiabatic_output = ''
 allocate(thermal_adiabatic_outputID       (maxval(homogenization_Noutput),maxNinstance),source=undefined_ID)
 allocate(thermal_adiabatic_Noutput        (maxNinstance),                               source=0_pInt) 

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

   if (section > 0_pInt ) then; if (thermal_type(section) == THERMAL_adiabatic_ID) then             ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran

     instance = thermal_typeInstance(section)                                                       ! which instance of my thermal is present homog
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('temperature')
             thermal_adiabatic_Noutput(instance) = thermal_adiabatic_Noutput(instance) + 1_pInt
             thermal_adiabatic_outputID(thermal_adiabatic_Noutput(instance),instance) = temperature_ID
             thermal_adiabatic_output(thermal_adiabatic_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select

     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do section = 1_pInt, size(thermal_type)
   if (thermal_type(section) == THERMAL_adiabatic_ID) then
     NofMyHomog=count(material_homog==section)
     instance = thermal_typeInstance(section)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,thermal_adiabatic_Noutput(instance)
       select case(thermal_adiabatic_outputID(o,instance))
         case(temperature_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          thermal_adiabatic_sizePostResult(o,instance) = mySize
          thermal_adiabatic_sizePostResults(instance)  = thermal_adiabatic_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop

! allocate state arrays
     sizeState = 1_pInt
     thermalState(section)%sizeState = sizeState
     thermalState(section)%sizePostResults = thermal_adiabatic_sizePostResults(instance)
     allocate(thermalState(section)%state0   (sizeState,NofMyHomog), source=thermal_initialT(section))
     allocate(thermalState(section)%subState0(sizeState,NofMyHomog), source=thermal_initialT(section))
     allocate(thermalState(section)%state    (sizeState,NofMyHomog), source=thermal_initialT(section))

     nullify(thermalMapping(section)%p)
     thermalMapping(section)%p => mappingHomogenization(1,:,:)
     deallocate(temperature(section)%p)
     temperature(section)%p => thermalState(section)%state(1,:)
     deallocate(temperatureRate(section)%p)
     allocate  (temperatureRate(section)%p(NofMyHomog), source=0.0_pReal)
     
   endif
 
 enddo initializeInstances
end subroutine thermal_adiabatic_init

!--------------------------------------------------------------------------------------------------
!> @brief  calculates adiabatic change in temperature based on local heat generation model  
!--------------------------------------------------------------------------------------------------
function thermal_adiabatic_updateState(subdt, ip, el)
 use numerics, only: &
   err_thermal_tolAbs, &
   err_thermal_tolRel
 use material, only: &
   mappingHomogenization, &
   thermalState, &
   temperature, &
   temperatureRate, &
   thermalMapping
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   subdt
 logical,                    dimension(2)                             :: &
   thermal_adiabatic_updateState
 integer(pInt) :: &
   homog, &
   offset
 real(pReal) :: &
   T, Tdot, dTdot_dT  

 homog  = mappingHomogenization(2,ip,el)
 offset = mappingHomogenization(1,ip,el)
 
 T = thermalState(homog)%subState0(1,offset)
 call thermal_adiabatic_getSourceAndItsTangent(Tdot, dTdot_dT, T, ip, el)
 T = T + subdt*Tdot/(thermal_adiabatic_getSpecificHeat(ip,el)*thermal_adiabatic_getMassDensity(ip,el))
 
 thermal_adiabatic_updateState = [     abs(T - thermalState(homog)%state(1,offset)) &
                                    <= err_thermal_tolAbs &
                                  .or. abs(T - thermalState(homog)%state(1,offset)) &
                                    <= err_thermal_tolRel*abs(thermalState(homog)%state(1,offset)), &
                                  .true.]

 temperature    (homog)%p(thermalMapping(homog)%p(ip,el)) = T  
 temperatureRate(homog)%p(thermalMapping(homog)%p(ip,el)) = &
   (thermalState(homog)%state(1,offset) - thermalState(homog)%subState0(1,offset))/(subdt+tiny(0.0_pReal))
 
end function thermal_adiabatic_updateState

!--------------------------------------------------------------------------------------------------
!> @brief returns heat generation rate
!--------------------------------------------------------------------------------------------------
subroutine thermal_adiabatic_getSourceAndItsTangent(Tdot, dTdot_dT, T, ip, el)
 use math, only: &
   math_Mandel6to33
 use material, only: &
   homogenization_Ngrains, &
   mappingHomogenization, &
   mappingConstitutive, &
   thermal_typeInstance, &
   phase_Nsources, &
   phase_source, &
   SOURCE_thermal_dissipation_ID, &
   SOURCE_thermal_externalheat_ID
 use source_thermal_dissipation, only: &
   source_thermal_dissipation_getRateAndItsTangent
 use source_thermal_externalheat, only: &
   source_thermal_externalheat_getRateAndItsTangent
 use crystallite, only: &
   crystallite_Tstar_v, &
   crystallite_Lp  

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), intent(in) :: &
   T
 real(pReal), intent(out) :: &
   Tdot, dTdot_dT
 real(pReal) :: &
   my_Tdot, my_dTdot_dT
 integer(pInt) :: &
   phase, &
   homog, &
   offset, &
   instance, &
   grain, &
   source
   
 homog  = mappingHomogenization(2,ip,el)
 offset = mappingHomogenization(1,ip,el)
 instance = thermal_typeInstance(homog)
  
 Tdot = 0.0_pReal
 dTdot_dT = 0.0_pReal
 do grain = 1, homogenization_Ngrains(homog)
   phase = mappingConstitutive(2,grain,ip,el)
   do source = 1, phase_Nsources(phase)
     select case(phase_source(source,phase))                                                   
       case (SOURCE_thermal_dissipation_ID)
        call source_thermal_dissipation_getRateAndItsTangent(my_Tdot, my_dTdot_dT, &
                                                             crystallite_Tstar_v(1:6,grain,ip,el), &
                                                             crystallite_Lp(1:3,1:3,grain,ip,el), &
                                                             grain, ip, el)

       case (SOURCE_thermal_externalheat_ID)
        call source_thermal_externalheat_getRateAndItsTangent(my_Tdot, my_dTdot_dT, &
                                                              grain, ip, el)

       case default
        my_Tdot = 0.0_pReal
        my_dTdot_dT = 0.0_pReal
     end select
     Tdot = Tdot + my_Tdot
     dTdot_dT = dTdot_dT + my_dTdot_dT
   enddo  
 enddo
 
 Tdot = Tdot/homogenization_Ngrains(homog)
 dTdot_dT = dTdot_dT/homogenization_Ngrains(homog)
 
end subroutine thermal_adiabatic_getSourceAndItsTangent
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized specific heat capacity
!--------------------------------------------------------------------------------------------------
function thermal_adiabatic_getSpecificHeat(ip,el)
 use lattice, only: &
   lattice_specificHeat
 use material, only: &
   homogenization_Ngrains, &
   mappingHomogenization, &
   material_phase
 use mesh, only: &
   mesh_element
 use crystallite, only: &
   crystallite_push33ToRef

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: &
   thermal_adiabatic_getSpecificHeat
 integer(pInt) :: &
   homog, grain
  
 thermal_adiabatic_getSpecificHeat = 0.0_pReal
 
 homog  = mappingHomogenization(2,ip,el)
  
 do grain = 1, homogenization_Ngrains(mesh_element(3,el))
   thermal_adiabatic_getSpecificHeat = thermal_adiabatic_getSpecificHeat + &
    lattice_specificHeat(material_phase(grain,ip,el))
 enddo

 thermal_adiabatic_getSpecificHeat = &
   thermal_adiabatic_getSpecificHeat/ &
   homogenization_Ngrains(mesh_element(3,el))
 
end function thermal_adiabatic_getSpecificHeat
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized mass density
!--------------------------------------------------------------------------------------------------
function thermal_adiabatic_getMassDensity(ip,el)
 use lattice, only: &
   lattice_massDensity
 use material, only: &
   homogenization_Ngrains, &
   mappingHomogenization, &
   material_phase
 use mesh, only: &
   mesh_element
 use crystallite, only: &
   crystallite_push33ToRef

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: &
   thermal_adiabatic_getMassDensity
 integer(pInt) :: &
   homog, grain
  
 thermal_adiabatic_getMassDensity = 0.0_pReal
 
 homog  = mappingHomogenization(2,ip,el)
  
 do grain = 1, homogenization_Ngrains(mesh_element(3,el))
   thermal_adiabatic_getMassDensity = thermal_adiabatic_getMassDensity + &
    lattice_massDensity(material_phase(grain,ip,el))
 enddo

 thermal_adiabatic_getMassDensity = &
   thermal_adiabatic_getMassDensity/ &
   homogenization_Ngrains(mesh_element(3,el))
 
end function thermal_adiabatic_getMassDensity
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of thermal results
!--------------------------------------------------------------------------------------------------
function thermal_adiabatic_postResults(ip,el)
 use material, only: &
   mappingHomogenization, &
   thermal_typeInstance, &
   thermalMapping, &
   temperature

 implicit none
 integer(pInt),              intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(thermal_adiabatic_sizePostResults(thermal_typeInstance(mappingHomogenization(2,ip,el)))) :: &
   thermal_adiabatic_postResults

 integer(pInt) :: &
   instance, homog, offset, o, c
   
 homog     = mappingHomogenization(2,ip,el)
 offset    = thermalMapping(homog)%p(ip,el)
 instance  = thermal_typeInstance(homog)

 c = 0_pInt
 thermal_adiabatic_postResults = 0.0_pReal

 do o = 1_pInt,thermal_adiabatic_Noutput(instance)
    select case(thermal_adiabatic_outputID(o,instance))
 
      case (temperature_ID)
        thermal_adiabatic_postResults(c+1_pInt) = temperature(homog)%p(offset)
        c = c + 1
    end select
 enddo
end function thermal_adiabatic_postResults

end module thermal_adiabatic
