!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for temperature evolution from heat conduction
!> @details to be done
!--------------------------------------------------------------------------------------------------
module thermal_conduction
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   thermal_conduction_sizePostResults                                                           !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   thermal_conduction_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   thermal_conduction_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   thermal_conduction_Noutput                                                                   !< number of outputs per instance of this damage 

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 temperature_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   thermal_conduction_outputID                                                                  !< ID of each post result output


 public :: &
   thermal_conduction_init, &
   thermal_conduction_getSourceAndItsTangent, &
   thermal_conduction_getConductivity33, &
   thermal_conduction_getSpecificHeat, &
   thermal_conduction_getMassDensity, &
   thermal_conduction_putTemperatureAndItsRate, &
   thermal_conduction_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine thermal_conduction_init(fileUnit)
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
   THERMAL_conduction_label, &
   THERMAL_conduction_ID, &
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

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: maxNinstance,mySize=0_pInt,section,instance,o
 integer(pInt) :: sizeState
 integer(pInt) :: NofMyHomog   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  thermal_'//THERMAL_CONDUCTION_label//' init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(thermal_type == THERMAL_conduction_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 allocate(thermal_conduction_sizePostResults(maxNinstance),                               source=0_pInt)
 allocate(thermal_conduction_sizePostResult (maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(thermal_conduction_output         (maxval(homogenization_Noutput),maxNinstance))
          thermal_conduction_output = ''
 allocate(thermal_conduction_outputID       (maxval(homogenization_Noutput),maxNinstance),source=undefined_ID)
 allocate(thermal_conduction_Noutput        (maxNinstance),                               source=0_pInt) 

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

   if (section > 0_pInt ) then; if (thermal_type(section) == THERMAL_conduction_ID) then             ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran

     instance = thermal_typeInstance(section)                                                       ! which instance of my thermal is present homog
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case ('temperature')
             thermal_conduction_Noutput(instance) = thermal_conduction_Noutput(instance) + 1_pInt
             thermal_conduction_outputID(thermal_conduction_Noutput(instance),instance) = temperature_ID
             thermal_conduction_output(thermal_conduction_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
          end select

     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do section = 1_pInt, size(thermal_type)
   if (thermal_type(section) == THERMAL_conduction_ID) then
     NofMyHomog=count(material_homog==section)
     instance = thermal_typeInstance(section)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,thermal_conduction_Noutput(instance)
       select case(thermal_conduction_outputID(o,instance))
         case(temperature_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          thermal_conduction_sizePostResult(o,instance) = mySize
          thermal_conduction_sizePostResults(instance)  = thermal_conduction_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop

! allocate state arrays
     sizeState = 0_pInt
     thermalState(section)%sizeState = sizeState
     thermalState(section)%sizePostResults = thermal_conduction_sizePostResults(instance)
     allocate(thermalState(section)%state0   (sizeState,NofMyHomog))
     allocate(thermalState(section)%subState0(sizeState,NofMyHomog))
     allocate(thermalState(section)%state    (sizeState,NofMyHomog))

     nullify(thermalMapping(section)%p)
     thermalMapping(section)%p => mappingHomogenization(1,:,:)
     deallocate(temperature    (section)%p)
     allocate  (temperature    (section)%p(NofMyHomog), source=thermal_initialT(section))
     deallocate(temperatureRate(section)%p)
     allocate  (temperatureRate(section)%p(NofMyHomog), source=0.0_pReal)
     
   endif
 
 enddo initializeInstances
end subroutine thermal_conduction_init

!--------------------------------------------------------------------------------------------------
!> @brief returns heat generation rate
!--------------------------------------------------------------------------------------------------
subroutine thermal_conduction_getSourceAndItsTangent(Tdot, dTdot_dT, T, ip, el)
 use math, only: &
   math_Mandel6to33
 use material, only: &
   homogenization_Ngrains, &
   mappingHomogenization, &
   phaseAt, phasememberAt, &
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
   phase = phaseAt(grain,ip,el)
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
 
end subroutine thermal_conduction_getSourceAndItsTangent
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized thermal conductivity in reference configuration
!--------------------------------------------------------------------------------------------------
function thermal_conduction_getConductivity33(ip,el)
 use lattice, only: &
   lattice_thermalConductivity33
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
 real(pReal), dimension(3,3) :: &
   thermal_conduction_getConductivity33
 integer(pInt) :: &
   homog, &
   grain
   
 homog  = mappingHomogenization(2,ip,el)
  
 thermal_conduction_getConductivity33 = 0.0_pReal
 do grain = 1, homogenization_Ngrains(mesh_element(3,el))
   thermal_conduction_getConductivity33 = thermal_conduction_getConductivity33 + &
    crystallite_push33ToRef(grain,ip,el,lattice_thermalConductivity33(:,:,material_phase(grain,ip,el)))
 enddo

 thermal_conduction_getConductivity33 = &
   thermal_conduction_getConductivity33/ &
   homogenization_Ngrains(mesh_element(3,el))
 
end function thermal_conduction_getConductivity33
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized specific heat capacity
!--------------------------------------------------------------------------------------------------
function thermal_conduction_getSpecificHeat(ip,el)
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
   thermal_conduction_getSpecificHeat
 integer(pInt) :: &
   homog, grain
  
 thermal_conduction_getSpecificHeat = 0.0_pReal
 
 homog  = mappingHomogenization(2,ip,el)
  
 do grain = 1, homogenization_Ngrains(mesh_element(3,el))
   thermal_conduction_getSpecificHeat = thermal_conduction_getSpecificHeat + &
    lattice_specificHeat(material_phase(grain,ip,el))
 enddo

 thermal_conduction_getSpecificHeat = &
   thermal_conduction_getSpecificHeat/ &
   homogenization_Ngrains(mesh_element(3,el))
 
end function thermal_conduction_getSpecificHeat
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized mass density
!--------------------------------------------------------------------------------------------------
function thermal_conduction_getMassDensity(ip,el)
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
   thermal_conduction_getMassDensity
 integer(pInt) :: &
   homog, grain
  
 thermal_conduction_getMassDensity = 0.0_pReal
 
 homog  = mappingHomogenization(2,ip,el)
  
 do grain = 1, homogenization_Ngrains(mesh_element(3,el))
   thermal_conduction_getMassDensity = thermal_conduction_getMassDensity + &
    lattice_massDensity(material_phase(grain,ip,el))
 enddo

 thermal_conduction_getMassDensity = &
   thermal_conduction_getMassDensity/ &
   homogenization_Ngrains(mesh_element(3,el))
 
end function thermal_conduction_getMassDensity
 
!--------------------------------------------------------------------------------------------------
!> @brief updates thermal state with solution from heat conduction PDE
!--------------------------------------------------------------------------------------------------
subroutine thermal_conduction_putTemperatureAndItsRate(T,Tdot,ip,el)
 use material, only: &
   mappingHomogenization, &
   temperature, &
   temperatureRate, &
   thermalMapping

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   T, &
   Tdot
 integer(pInt) :: &
   homog, &
   offset  
 
 homog  = mappingHomogenization(2,ip,el)
 offset = thermalMapping(homog)%p(ip,el)
 temperature    (homog)%p(offset) = T
 temperatureRate(homog)%p(offset) = Tdot

end subroutine thermal_conduction_putTemperatureAndItsRate
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of thermal results
!--------------------------------------------------------------------------------------------------
function thermal_conduction_postResults(ip,el)
 use material, only: &
   mappingHomogenization, &
   thermal_typeInstance, &
   temperature, &
   thermalMapping

 implicit none
 integer(pInt),              intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(thermal_conduction_sizePostResults(thermal_typeInstance(mappingHomogenization(2,ip,el)))) :: &
   thermal_conduction_postResults

 integer(pInt) :: &
   instance, homog, offset, o, c
   
 homog     = mappingHomogenization(2,ip,el)
 offset    = thermalMapping(homog)%p(ip,el)
 instance  = thermal_typeInstance(homog)

 c = 0_pInt
 thermal_conduction_postResults = 0.0_pReal

 do o = 1_pInt,thermal_conduction_Noutput(instance)
    select case(thermal_conduction_outputID(o,instance))
 
      case (temperature_ID)
        thermal_conduction_postResults(c+1_pInt) = temperature(homog)%p(offset)
        c = c + 1
    end select
 enddo
end function thermal_conduction_postResults

end module thermal_conduction
