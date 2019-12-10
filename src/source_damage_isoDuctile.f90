!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating isotropic ductile damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_damage_isoDuctile
 use prec
 use debug
 use IO
 use discretization
 use material
 use config
 use results

 implicit none
 private
 integer,                       dimension(:),           allocatable,         public, protected :: &
   source_damage_isoDuctile_offset, &                                                                 !< which source is my current damage mechanism?
   source_damage_isoDuctile_instance                                                                  !< instance of damage source mechanism

 character(len=64),             dimension(:,:),         allocatable, target, public :: &
   source_damage_isoDuctile_output                                                                    !< name of each post result output


 enum, bind(c) 
   enumerator :: undefined_ID, &
                 damage_drivingforce_ID
 end enum                                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 ToDo

 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   real(pReal) :: &
     critPlasticStrain, &
     N, &
     aTol
   integer(kind(undefined_ID)), allocatable, dimension(:) :: &
     outputID
 end type tParameters

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)


 public :: &
   source_damage_isoDuctile_init, &
   source_damage_isoDuctile_dotState, &
   source_damage_isoDuctile_getRateAndItsTangent, &
   source_damage_isoDuctile_Results

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoDuctile_init

 integer :: Ninstance,phase,instance,source,sourceOffset
 integer :: NofMyPhase,p,i
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
 integer(kind(undefined_ID)) :: &
   outputID

 character(len=pStringLen) :: &
   extmsg = ''
 character(len=65536), dimension(:), allocatable :: &
   outputs

 write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_DAMAGE_ISODUCTILE_LABEL//' init  -+>>>'

 Ninstance = count(phase_source == SOURCE_damage_isoDuctile_ID)
 if (Ninstance == 0) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
   write(6,'(a16,1x,i5,/)') '# instances:',Ninstance
 
 allocate(source_damage_isoDuctile_offset(material_Nphase), source=0)
 allocate(source_damage_isoDuctile_instance(material_Nphase), source=0)
 do phase = 1, material_Nphase
   source_damage_isoDuctile_instance(phase) = count(phase_source(:,1:phase) == source_damage_isoDuctile_ID)
   do source = 1, phase_Nsources(phase)
     if (phase_source(source,phase) == source_damage_isoDuctile_ID) &
       source_damage_isoDuctile_offset(phase) = source
   enddo    
 enddo
   
 allocate(source_damage_isoDuctile_output(maxval(phase_Noutput),Ninstance))
          source_damage_isoDuctile_output = ''

 allocate(param(Ninstance))
 
 do p=1, size(config_phase)
   if (all(phase_source(:,p) /= SOURCE_DAMAGE_ISODUCTILE_ID)) cycle
   associate(prm => param(source_damage_isoDuctile_instance(p)), &
             config => config_phase(p))
             
   prm%aTol              = config%getFloat('isoductile_atol',defaultVal = 1.0e-3_pReal)

   prm%N                 = config%getFloat('isoductile_ratesensitivity')
   prm%critPlasticStrain = config%getFloat('isoductile_criticalplasticstrain')
   
   ! sanity checks
   if (prm%aTol                 < 0.0_pReal) extmsg = trim(extmsg)//' isoductile_atol'
   
   if (prm%N                   <= 0.0_pReal) extmsg = trim(extmsg)//' isoductile_ratesensitivity'
   if (prm%critPlasticStrain   <= 0.0_pReal) extmsg = trim(extmsg)//' isoductile_criticalplasticstrain'
   
!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
   if (extmsg /= '') &
     call IO_error(211,ext_msg=trim(extmsg)//'('//SOURCE_DAMAGE_ISODUCTILE_LABEL//')')

!--------------------------------------------------------------------------------------------------
!  output pararameters
   outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   do i=1, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
     
       case ('isoductile_drivingforce')
         source_damage_isoDuctile_output(i,source_damage_isoDuctile_instance(p)) = outputs(i)
         prm%outputID = [prm%outputID, damage_drivingforce_ID]

     end select

   enddo

   end associate
   
   phase = p
   NofMyPhase=count(material_phaseAt==phase) * discretization_nIP
   instance = source_damage_isoDuctile_instance(phase)
   sourceOffset = source_damage_isoDuctile_offset(phase)

   call material_allocateSourceState(phase,sourceOffset,NofMyPhase,1,1,0)
   sourceState(phase)%p(sourceOffset)%aTolState=param(instance)%aTol
 
 enddo
 
end subroutine source_damage_isoDuctile_init

!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoDuctile_dotState(ipc, ip, el)

 integer, intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 integer :: &
   phase, constituent, instance, homog, sourceOffset, damageOffset

 phase = material_phaseAt(ipc,el)
 constituent = material_phasememberAt(ipc,ip,el)
 instance = source_damage_isoDuctile_instance(phase)
 sourceOffset = source_damage_isoDuctile_offset(phase)
 homog = material_homogenizationAt(el)
 damageOffset = damageMapping(homog)%p(ip,el)

 sourceState(phase)%p(sourceOffset)%dotState(1,constituent) = &
   sum(plasticState(phase)%slipRate(:,constituent))/ &
   ((damage(homog)%p(damageOffset))**param(instance)%N)/ & 
   param(instance)%critPlasticStrain 

end subroutine source_damage_isoDuctile_dotState
 
!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoDuctile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

 integer, intent(in) :: &
   phase, &
   constituent
 real(pReal),  intent(in) :: &
   phi
 real(pReal),  intent(out) :: &
   localphiDot, &
   dLocalphiDot_dPhi
 integer :: &
   sourceOffset

 sourceOffset = source_damage_isoDuctile_offset(phase)
 
 localphiDot = 1.0_pReal &
             - sourceState(phase)%p(sourceOffset)%state(1,constituent) * phi
 
 dLocalphiDot_dPhi = -sourceState(phase)%p(sourceOffset)%state(1,constituent)
 
end subroutine source_damage_isoDuctile_getRateAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoDuctile_results(phase,group)

  integer, intent(in) :: phase
  character(len=*), intent(in) :: group
#if defined(PETSc) || defined(DAMASK_HDF5)  
  integer :: sourceOffset, o, instance
   
  instance     = source_damage_isoDuctile_instance(phase)
  sourceOffset = source_damage_isoDuctile_offset(phase)

   associate(prm => param(instance), stt => sourceState(phase)%p(sourceOffset)%state)
   outputsLoop: do o = 1,size(prm%outputID)
     select case(prm%outputID(o))
       case (damage_drivingforce_ID)
         call results_writeDataset(group,stt,'tbd','driving force','tbd')
     end select
   enddo outputsLoop
   end associate
#endif

end subroutine source_damage_isoDuctile_results


end module source_damage_isoDuctile
