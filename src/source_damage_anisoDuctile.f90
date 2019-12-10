!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic ductile damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_damage_anisoDuctile
  use prec
  use debug
  use IO
  use math
  use discretization
  use material
  use config
  use results
 
  implicit none
  private
   
  integer,                       dimension(:),           allocatable,         public, protected :: &
    source_damage_anisoDuctile_offset, &                                                            !< which source is my current damage mechanism?
    source_damage_anisoDuctile_instance                                                             !< instance of damage source mechanism
 
  character(len=64),             dimension(:,:),         allocatable, target, public  :: &
    source_damage_anisoDuctile_output                                                               !< name of each post result output
    
 
  enum, bind(c) 
    enumerator :: undefined_ID, &
                  damage_drivingforce_ID
  end enum 
 
 
  type, private :: tParameters                                                                      !< container type for internal constitutive parameters
    real(pReal) :: &
      aTol, &
      N
    real(pReal), dimension(:), allocatable :: &
      critPlasticStrain
    integer :: &
      totalNslip
    integer, dimension(:), allocatable :: &
      Nslip
    integer(kind(undefined_ID)), allocatable, dimension(:) :: &
      outputID
  end type tParameters
 
  type(tParameters), dimension(:), allocatable, private :: param                                    !< containers of constitutive parameters (len Ninstance)
 
 
  public :: &
    source_damage_anisoDuctile_init, &
    source_damage_anisoDuctile_dotState, &
    source_damage_anisoDuctile_getRateAndItsTangent, &
    source_damage_anisoDuctile_results

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisoDuctile_init
   
  integer :: Ninstance,phase,instance,source,sourceOffset
  integer :: NofMyPhase,p ,i
 
  integer,              dimension(0), parameter :: emptyIntArray    = [integer::]
  character(len=65536), dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
  integer(kind(undefined_ID)) :: &
    outputID
 
  character(len=pStringLen) :: &
    extmsg = ''
  character(len=65536), dimension(:), allocatable :: &
    outputs
 
  write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_DAMAGE_ANISODUCTILE_LABEL//' init  -+>>>'
 
  Ninstance = count(phase_source == SOURCE_damage_anisoDuctile_ID)
  if (Ninstance == 0) return
  
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance
  
  allocate(source_damage_anisoDuctile_offset(size(config_phase)), source=0)
  allocate(source_damage_anisoDuctile_instance(size(config_phase)), source=0)
  do phase = 1, size(config_phase)
    source_damage_anisoDuctile_instance(phase) = count(phase_source(:,1:phase) == source_damage_anisoDuctile_ID)
    do source = 1, phase_Nsources(phase)
      if (phase_source(source,phase) == source_damage_anisoDuctile_ID) &
        source_damage_anisoDuctile_offset(phase) = source
    enddo    
  enddo
    
  allocate(source_damage_anisoDuctile_output(maxval(phase_Noutput),Ninstance))
           source_damage_anisoDuctile_output = ''
 
 
  allocate(param(Ninstance))
  
  do p=1, size(config_phase)
    if (all(phase_source(:,p) /= SOURCE_DAMAGE_ANISODUCTILE_ID)) cycle
    associate(prm => param(source_damage_anisoDuctile_instance(p)), &
              config => config_phase(p))
              
    prm%aTol   = config%getFloat('anisoductile_atol',defaultVal = 1.0e-3_pReal)
 
    prm%N      = config%getFloat('anisoductile_ratesensitivity')
    prm%totalNslip = sum(prm%Nslip)
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
    if (extmsg /= '')  call IO_error(211,ext_msg=trim(extmsg)//'('//SOURCE_DAMAGE_ANISODUCTILE_LABEL//')')

!--------------------------------------------------------------------------------------------------
!  output pararameters
    outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
    allocate(prm%outputID(0))
    do i=1, size(outputs)
      outputID = undefined_ID
      select case(outputs(i))
      
        case ('anisoductile_drivingforce')
          source_damage_anisoDuctile_output(i,source_damage_anisoDuctile_instance(p)) = outputs(i)
          prm%outputID = [prm%outputID, damage_drivingforce_ID]
 
      end select
 
    enddo
 
    end associate
    
    phase = p
    
    NofMyPhase=count(material_phaseAt==phase) * discretization_nIP
    instance = source_damage_anisoDuctile_instance(phase)
    sourceOffset = source_damage_anisoDuctile_offset(phase)
 
    call material_allocateSourceState(phase,sourceOffset,NofMyPhase,1,1,0)
      sourceState(phase)%p(sourceOffset)%aTolState=param(instance)%aTol
    
  enddo
  
end subroutine source_damage_anisoDuctile_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisoDuctile_dotState(ipc, ip, el)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element
  integer :: &
    phase, &
    constituent, &
    sourceOffset, &
    homog, damageOffset, &
    instance, &
    i
 
  phase = material_phaseAt(ipc,el)
  constituent = material_phasememberAt(ipc,ip,el)
  instance = source_damage_anisoDuctile_instance(phase)
  sourceOffset = source_damage_anisoDuctile_offset(phase)
  homog = material_homogenizationAt(el)
  damageOffset = damageMapping(homog)%p(ip,el)
 
 
  do i = 1, param(instance)%totalNslip
      sourceState(phase)%p(sourceOffset)%dotState(1,constituent) = &
        sourceState(phase)%p(sourceOffset)%dotState(1,constituent) + &
        plasticState(phase)%slipRate(i,constituent)/ &
        ((damage(homog)%p(damageOffset))**param(instance)%N)/param(instance)%critPlasticStrain(i) 
  enddo
 
end subroutine source_damage_anisoDuctile_dotState


!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisoDuctile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

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
 
  sourceOffset = source_damage_anisoDuctile_offset(phase)
  
  localphiDot = 1.0_pReal &
              - sourceState(phase)%p(sourceOffset)%state(1,constituent) * phi
  
  dLocalphiDot_dPhi = -sourceState(phase)%p(sourceOffset)%state(1,constituent)
 
end subroutine source_damage_anisoDuctile_getRateAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisoDuctile_results(phase,group)

  integer, intent(in) :: phase
  character(len=*), intent(in) :: group
#if defined(PETSc) || defined(DAMASK_HDF5)  
  integer :: sourceOffset, o, instance
   
  instance     = source_damage_anisoDuctile_instance(phase)
  sourceOffset = source_damage_anisoDuctile_offset(phase)

   associate(prm => param(instance), stt => sourceState(phase)%p(sourceOffset)%state)
   outputsLoop: do o = 1,size(prm%outputID)
     select case(prm%outputID(o))
       case (damage_drivingforce_ID)
         call results_writeDataset(group,stt,'tbd','driving force','tbd')
     end select
   enddo outputsLoop
   end associate
#endif

end subroutine source_damage_anisoDuctile_results

end module source_damage_anisoDuctile
