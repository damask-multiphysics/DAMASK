!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating isotropic brittle damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_damage_isoBrittle
  use prec
  use debug
  use IO
  use math
  use material
  use config

  implicit none
  private
  integer,                       dimension(:),           allocatable,         public, protected :: &
    source_damage_isoBrittle_offset, &
    source_damage_isoBrittle_instance
  integer,                       dimension(:,:),         allocatable, target, public :: &
    source_damage_isoBrittle_sizePostResult
  character(len=64),             dimension(:,:),         allocatable, target, public :: &
    source_damage_isoBrittle_output

  enum, bind(c) 
    enumerator :: &
      undefined_ID, &
      damage_drivingforce_ID
  end enum


  type, private :: tParameters                                                                      !< container type for internal constitutive parameters
    real(pReal) :: &
      critStrainEnergy, &
      N, &
      aTol
    integer(kind(undefined_ID)), allocatable, dimension(:) :: &
      outputID
  end type tParameters

  type(tParameters), dimension(:), allocatable  :: param                                            !< containers of constitutive parameters (len Ninstance)


  public :: &
    source_damage_isoBrittle_init, &
    source_damage_isoBrittle_deltaState, &
    source_damage_isoBrittle_getRateAndItsTangent, &
    source_damage_isoBrittle_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoBrittle_init

  integer :: Ninstance,phase,instance,source,sourceOffset
  integer :: NofMyPhase,p,i   
  character(len=65536), dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
  integer(kind(undefined_ID)) :: &
    outputID
 
  character(len=pStringLen) :: &
    extmsg = ''
  character(len=65536), dimension(:), allocatable :: &
    outputs
  
  write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_DAMAGE_ISOBRITTLE_LABEL//' init  -+>>>'
 
  Ninstance = count(phase_source == SOURCE_damage_isoBrittle_ID)
  if (Ninstance == 0) return
  
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance
  
  allocate(source_damage_isoBrittle_offset(material_Nphase), source=0)
  allocate(source_damage_isoBrittle_instance(material_Nphase), source=0)
  do phase = 1, material_Nphase
    source_damage_isoBrittle_instance(phase) = count(phase_source(:,1:phase) == source_damage_isoBrittle_ID)
    do source = 1, phase_Nsources(phase)
      if (phase_source(source,phase) == source_damage_isoBrittle_ID) &
        source_damage_isoBrittle_offset(phase) = source
    enddo    
  enddo
    
  allocate(source_damage_isoBrittle_sizePostResult(maxval(phase_Noutput),Ninstance),source=0)
  allocate(source_damage_isoBrittle_output(maxval(phase_Noutput),Ninstance))
           source_damage_isoBrittle_output = ''
 
  allocate(param(Ninstance))
  
  do p=1, size(config_phase)
    if (all(phase_source(:,p) /= SOURCE_DAMAGE_ISOBRITTLE_ID)) cycle
    associate(prm => param(source_damage_isoBrittle_instance(p)), &
              config => config_phase(p))
              
    prm%aTol             = config%getFloat('isobrittle_atol',defaultVal = 1.0e-3_pReal)
    
    prm%N                = config%getFloat('isobrittle_n')
    prm%critStrainEnergy = config%getFloat('isobrittle_criticalstrainenergy')
    
    ! sanity checks
    if (prm%aTol                < 0.0_pReal) extmsg = trim(extmsg)//' isobrittle_atol'
    
    if (prm%N                  <= 0.0_pReal) extmsg = trim(extmsg)//' isobrittle_n'
    if (prm%critStrainEnergy   <= 0.0_pReal) extmsg = trim(extmsg)//' isobrittle_criticalstrainenergy'
   
!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') &
      call IO_error(211,ext_msg=trim(extmsg)//'('//SOURCE_DAMAGE_ISOBRITTLE_LABEL//')')

!--------------------------------------------------------------------------------------------------
!  output pararameters
    outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
    allocate(prm%outputID(0))
    do i=1, size(outputs)
      outputID = undefined_ID
      select case(outputs(i))
      
        case ('isobrittle_drivingforce')
          source_damage_isoBrittle_sizePostResult(i,source_damage_isoBrittle_instance(p)) = 1
          source_damage_isoBrittle_output(i,source_damage_isoBrittle_instance(p)) = outputs(i)
          prm%outputID = [prm%outputID, damage_drivingforce_ID]
            
      end select

    enddo

    end associate
    
   phase = p
   
   NofMyPhase=count(material_phase==phase)
   instance = source_damage_isoBrittle_instance(phase)
   sourceOffset = source_damage_isoBrittle_offset(phase)
      
   call material_allocateSourceState(phase,sourceOffset,NofMyPhase,1,1,1)
   sourceState(phase)%p(sourceOffset)%sizePostResults = sum(source_damage_isoBrittle_sizePostResult(:,instance))
   sourceState(phase)%p(sourceOffset)%aTolState=param(instance)%aTol
  
  enddo
 
end subroutine source_damage_isoBrittle_init

!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoBrittle_deltaState(C, Fe, ipc, ip, el)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element
  real(pReal),  intent(in), dimension(3,3) :: &
    Fe
  real(pReal),  intent(in), dimension(6,6) :: &
    C
  integer :: &
    phase, constituent, instance, sourceOffset
  real(pReal) :: &
    strain(6), &
    strainenergy

  phase = phaseAt(ipc,ip,el)                                                                        !< phase ID at ipc,ip,el
  constituent = phasememberAt(ipc,ip,el)                                                            !< state array offset for phase ID at ipc,ip,el
  ! ToDo: capability for multiple instances of SAME source within given phase. Needs Ninstance loop from here on!
  instance = source_damage_isoBrittle_instance(phase)                                               !< instance of damage_isoBrittle source
  sourceOffset = source_damage_isoBrittle_offset(phase)

                                    
  strain = 0.5_pReal*math_sym33to6(matmul(transpose(Fe),Fe)-math_I3)

  strainenergy = 2.0_pReal*sum(strain*matmul(C,strain))/param(instance)%critStrainEnergy
  ! ToDo: check strainenergy = 2.0_pReal*dot_product(strain,matmul(C,strain))/param(instance)%critStrainEnergy
  
  if (strainenergy > sourceState(phase)%p(sourceOffset)%subState0(1,constituent)) then
    sourceState(phase)%p(sourceOffset)%deltaState(1,constituent) = &
      strainenergy - sourceState(phase)%p(sourceOffset)%state(1,constituent)
  else
    sourceState(phase)%p(sourceOffset)%deltaState(1,constituent) = &
      sourceState(phase)%p(sourceOffset)%subState0(1,constituent) - &
      sourceState(phase)%p(sourceOffset)%state(1,constituent)
  endif
 
end subroutine source_damage_isoBrittle_deltaState
 
!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoBrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

  integer, intent(in) :: &
    phase, &
    constituent
  real(pReal),  intent(in) :: &
    phi
  real(pReal),  intent(out) :: &
    localphiDot, &
    dLocalphiDot_dPhi
  integer :: &
    instance, sourceOffset

  instance = source_damage_isoBrittle_instance(phase)
  sourceOffset = source_damage_isoBrittle_offset(phase)
  
  localphiDot = (1.0_pReal - phi)**(param(instance)%N - 1.0_pReal) - &
                phi*sourceState(phase)%p(sourceOffset)%state(1,constituent)
  dLocalphiDot_dPhi = - (param(instance)%N - 1.0_pReal)* &
                        (1.0_pReal - phi)**max(0.0_pReal,param(instance)%N - 2.0_pReal) &
                      - sourceState(phase)%p(sourceOffset)%state(1,constituent)
  
end subroutine source_damage_isoBrittle_getRateAndItsTangent
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of local damage results
!--------------------------------------------------------------------------------------------------
function source_damage_isoBrittle_postResults(phase, constituent)

  integer, intent(in) :: &
    phase, &
    constituent
  real(pReal), dimension(sum(source_damage_isoBrittle_sizePostResult(:, &
                          source_damage_isoBrittle_instance(phase)))) :: &
    source_damage_isoBrittle_postResults

  integer :: &
    instance, sourceOffset, o, c
    
  instance = source_damage_isoBrittle_instance(phase)
  sourceOffset = source_damage_isoBrittle_offset(phase)

  c = 0

  do o = 1,size(param(instance)%outputID)
     select case(param(instance)%outputID(o))
       case (damage_drivingforce_ID)
         source_damage_isoBrittle_postResults(c+1) = sourceState(phase)%p(sourceOffset)%state(1,constituent)
         c = c + 1

     end select
  enddo
end function source_damage_isoBrittle_postResults

end module source_damage_isoBrittle
