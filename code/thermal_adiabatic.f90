!--------------------------------------------------------------------------------------------------
! $Id: thermal_adiabatic.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating dislocation and twinning physics
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
   thermal_adiabatic_Noutput                                                                   !< number of outputs per instance of this damage 

 real(pReal),                         dimension(:),     allocatable,         public :: &
   thermal_adiabatic_aTol

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 temperature_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   thermal_adiabatic_outputID                                                                  !< ID of each post result output


 public :: &
   thermal_adiabatic_init, &
   thermal_adiabatic_stateInit, &
   thermal_adiabatic_aTolState, &
   thermal_adiabatic_dotState, &
   thermal_adiabatic_getTemperature, &
   thermal_adiabatic_putTemperature, &
   thermal_adiabatic_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine thermal_adiabatic_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use mesh, only: &
   mesh_maxNips, &
   mesh_NcpElems
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
   homogenization_maxNgrains, &
   phase_thermal, &
   phase_thermalInstance, &
   phase_Noutput, &
   LOCAL_THERMAL_ADIABATIC_label, &
   LOCAL_THERMAL_adiabatic_ID, &
   material_phase, &  
   thermalState, &
   MATERIAL_partPhase
 use numerics,only: &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), parameter :: MAXNCHUNKS = 7_pInt
 integer(pInt), dimension(1+2*MAXNCHUNKS) :: positions
 integer(pInt) :: maxNinstance,mySize=0_pInt,phase,instance,o
 integer(pInt) :: sizeState, sizeDotState
 integer(pInt) :: NofMyPhase   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 write(6,'(/,a)')   ' <<<+-  thermal_'//LOCAL_THERMAL_ADIABATIC_label//' init  -+>>>'
 write(6,'(a)')     ' $Id: thermal_adiabatic.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_thermal == LOCAL_THERMAL_adiabatic_ID),pInt)
 if (maxNinstance == 0_pInt) return
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(thermal_adiabatic_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(thermal_adiabatic_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(thermal_adiabatic_output(maxval(phase_Noutput),maxNinstance))
          thermal_adiabatic_output = ''
 allocate(thermal_adiabatic_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(thermal_adiabatic_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(thermal_adiabatic_aTol(maxNinstance),                       source=0.0_pReal) 

 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= MATERIAL_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo
 
 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif   
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next phase section
     phase = phase + 1_pInt                                                                         ! advance phase section counter
     cycle                                                                                          ! skip to next line
   endif

   if (phase > 0_pInt ) then; if (phase_thermal(phase) == LOCAL_THERMAL_adiabatic_ID) then               ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran

     instance = phase_thermalInstance(phase)                                                     ! which instance of my thermal is present phase
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

       case ('atol_adiabatic')
         thermal_adiabatic_aTol(instance) = IO_floatValue(line,positions,2_pInt)

     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do phase = 1_pInt, size(phase_thermal)
   if (phase_thermal(phase) == LOCAL_THERMAL_adiabatic_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_thermalInstance(phase)

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
! Determine size of state array
     sizeDotState              =   1_pInt
     sizeState                 =   1_pInt
     thermalState(phase)%sizeState = sizeState
     thermalState(phase)%sizeDotState = sizeDotState
     thermalState(phase)%sizePostResults = thermal_adiabatic_sizePostResults(instance)
     allocate(thermalState(phase)%aTolState           (sizeState),                source=0.0_pReal)
     allocate(thermalState(phase)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(thermalState(phase)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(thermalState(phase)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(thermalState(phase)%state               (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(thermalState(phase)%state_backup        (sizeState,NofMyPhase),     source=0.0_pReal)

     allocate(thermalState(phase)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(thermalState(phase)%deltaState          (sizeDotState,NofMyPhase),     source=0.0_pReal)
     allocate(thermalState(phase)%dotState_backup     (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(thermalState(phase)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
       allocate(thermalState(phase)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(thermalState(phase)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(thermalState(phase)%RKCK45dotState    (6,sizeDotState,NofMyPhase),source=0.0_pReal)

     call thermal_adiabatic_stateInit(phase,instance)
     call thermal_adiabatic_aTolState(phase,instance)
   endif
 
 enddo initializeInstances
end subroutine thermal_adiabatic_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant  NEW state values for a given instance of this thermal
!--------------------------------------------------------------------------------------------------
subroutine thermal_adiabatic_stateInit(phase,instance)
 use material, only: &
   thermalState
 use lattice, only: &
  lattice_referenceTemperature
 
 implicit none
 integer(pInt),              intent(in) :: instance                                                 !< number specifying the instance of the thermal
 integer(pInt),              intent(in) :: phase                                                    !< number specifying the phase of the thermal

 real(pReal), dimension(thermalState(phase)%sizeState) :: tempState

 tempState(1) = lattice_referenceTemperature(phase)
 thermalState(phase)%state = spread(tempState,2,size(thermalState(phase)%state(1,:)))
 thermalState(phase)%state0 = thermalState(phase)%state
 thermalState(phase)%partionedState0 = thermalState(phase)%state
end subroutine thermal_adiabatic_stateInit

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this thermal
!--------------------------------------------------------------------------------------------------
subroutine thermal_adiabatic_aTolState(phase,instance)
 use material, only: &
  thermalState

 implicit none
 integer(pInt), intent(in) ::  &
   phase, &
   instance                                                                                         ! number specifying the current instance of the thermal
 real(pReal), dimension(thermalState(phase)%sizeState) :: tempTol

 tempTol = thermal_adiabatic_aTol
 thermalState(phase)%aTolState = tempTol
end subroutine thermal_adiabatic_aTolState
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine thermal_adiabatic_dotState(Tstar_v, Lp, ipc, ip, el)
 use lattice, only: &
   lattice_massDensity, &
   lattice_specificHeat
 use material, only: &
   mappingConstitutive, &
   phase_thermalInstance, &
   thermalState
 use math, only: &
   math_Mandel6to33

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal),  intent(in), dimension(3,3) :: &
   Lp
 integer(pInt) :: &
   instance, phase, constituent 

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_thermalInstance(phase)
 
 thermalState(phase)%dotState(1,constituent) = &
    0.95_pReal &
  * sum(abs(math_Mandel6to33(Tstar_v)*Lp)) &
  / (lattice_massDensity(phase)*lattice_specificHeat(phase))
  
end subroutine thermal_adiabatic_dotState

!--------------------------------------------------------------------------------------------------
!> @brief returns temperature based on local damage model state layout 
!--------------------------------------------------------------------------------------------------
function thermal_adiabatic_getTemperature(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   ThermalState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: thermal_adiabatic_getTemperature
 
 thermal_adiabatic_getTemperature = &
   thermalState(mappingConstitutive(2,ipc,ip,el))%state(1,mappingConstitutive(1,ipc,ip,el))
 
end function thermal_adiabatic_getTemperature
 
!--------------------------------------------------------------------------------------------------
!> @brief returns temperature based on local damage model state layout 
!--------------------------------------------------------------------------------------------------
subroutine thermal_adiabatic_putTemperature(ipc, ip, el, localTemperature)
 use material, only: &
   mappingConstitutive, &
   ThermalState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   localTemperature
 
 thermalState(mappingConstitutive(2,ipc,ip,el))%state(1,mappingConstitutive(1,ipc,ip,el))= &
   localTemperature
 
end subroutine thermal_adiabatic_putTemperature
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function thermal_adiabatic_postResults(ipc,ip,el)
 use material, only: &
   mappingConstitutive, &
   phase_thermalInstance, &
   thermalState

 implicit none
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(thermal_adiabatic_sizePostResults(phase_thermalInstance(mappingConstitutive(2,ipc,ip,el)))) :: &
   thermal_adiabatic_postResults

 integer(pInt) :: &
   instance, phase, constituent, o, c
   
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance  = phase_thermalInstance(phase)

 c = 0_pInt
 thermal_adiabatic_postResults = 0.0_pReal

 do o = 1_pInt,thermal_adiabatic_Noutput(instance)
    select case(thermal_adiabatic_outputID(o,instance))
 
      case (temperature_ID)
        thermal_adiabatic_postResults(c+1_pInt) = thermalState(phase)%state(1,constituent)
        c = c + 1
    end select
 enddo
end function thermal_adiabatic_postResults

end module thermal_adiabatic
