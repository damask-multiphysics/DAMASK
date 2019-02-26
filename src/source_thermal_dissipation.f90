!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for thermal source due to plastic dissipation
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_thermal_dissipation
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   source_thermal_dissipation_offset, &                                                                 !< which source is my current thermal dissipation mechanism?
   source_thermal_dissipation_instance                                                                  !< instance of thermal dissipation source mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   source_thermal_dissipation_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   source_thermal_dissipation_output                                                                    !< name of each post result output

 real(pReal),                         dimension(:),           allocatable,        private :: &
   source_thermal_dissipation_coldworkCoeff


 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   real(pReal) :: &
     coldworkCoeff
 end type tParameters

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)


 public :: &
   source_thermal_dissipation_init, &
   source_thermal_dissipation_getRateAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_thermal_dissipation_init
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use material, only: &
   material_allocateSourceState, &
   phase_source, &
   phase_Nsources, &
   phase_Noutput, &
   SOURCE_thermal_dissipation_label, &
   SOURCE_thermal_dissipation_ID, &
   material_phase, &  
   sourceState
 use config, only: &
   config_phase, &
   material_Nphase, &
   MATERIAL_partPhase

 implicit none
 integer(pInt) :: Ninstance,instance,source,sourceOffset
 integer(pInt) :: NofMyPhase,p   

 write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_thermal_dissipation_label//' init  -+>>>'

 
 Ninstance = int(count(phase_source == SOURCE_thermal_dissipation_ID),pInt)
 if (Ninstance == 0_pInt) return
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',Ninstance
 
 allocate(source_thermal_dissipation_offset(material_Nphase), source=0_pInt)
 allocate(source_thermal_dissipation_instance(material_Nphase), source=0_pInt)
 do p = 1, material_Nphase
   source_thermal_dissipation_instance(p) = count(phase_source(:,1:p) == SOURCE_thermal_dissipation_ID)
   do source = 1, phase_Nsources(p)
     if (phase_source(source,p) == SOURCE_thermal_dissipation_ID) &
       source_thermal_dissipation_offset(p) = source
   enddo    
 enddo
   
 allocate(source_thermal_dissipation_sizePostResult(maxval(phase_Noutput),Ninstance),source=0_pInt)
 allocate(source_thermal_dissipation_output  (maxval(phase_Noutput),Ninstance))
          source_thermal_dissipation_output = ''

 allocate(source_thermal_dissipation_coldworkCoeff(Ninstance),                       source=0.0_pReal) 

 do p=1, size(config_phase)
   if (all(phase_source(:,p) /= SOURCE_THERMAL_DISSIPATION_ID)) cycle
   instance = source_thermal_dissipation_instance(p)
   source_thermal_dissipation_coldworkCoeff(instance) = config_phase(p)%getFloat('dissipation_coldworkcoeff')
   NofMyPhase=count(material_phase==p)
   sourceOffset = source_thermal_dissipation_offset(p)

   call material_allocateSourceState(p,sourceOffset,NofMyPhase,0_pInt,0_pInt,0_pInt)
 
 enddo
  
end subroutine source_thermal_dissipation_init


!--------------------------------------------------------------------------------------------------
!> @brief returns local vacancy generation rate 
!--------------------------------------------------------------------------------------------------
subroutine source_thermal_dissipation_getRateAndItsTangent(TDot, dTDOT_dT, Tstar, Lp, phase)

 implicit none
 integer(pInt), intent(in) :: &
   phase
 real(pReal),  intent(in), dimension(3,3) :: &
   Tstar
 real(pReal),  intent(in), dimension(3,3) :: &
   Lp
 real(pReal),  intent(out) :: &
   TDot, &
   dTDOT_dT
 integer(pInt) :: &
   instance

 instance = source_thermal_dissipation_instance(phase)
 
 TDot = source_thermal_dissipation_coldworkCoeff(instance)*sum(abs(Tstar*Lp))
 dTDOT_dT = 0.0_pReal       
 
end subroutine source_thermal_dissipation_getRateAndItsTangent
 
end module source_thermal_dissipation
