!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Denny Tjahjanto, Max-Planck-Institut für Eisenforschung GmbH
!> @brief homogenization manager, organizing deformation partitioning and stress homogenization
!--------------------------------------------------------------------------------------------------
module homogenization
  use prec
  use IO
  use config
  use debug
  use math
  use material
  use numerics
  use constitutive
  use crystallite
  use FEsolving
  use discretization
  use thermal_isothermal
  use thermal_adiabatic
  use thermal_conduction
  use damage_none
  use damage_local
  use damage_nonlocal
  use results
  use HDF5_utilities
 
  implicit none
  private

!--------------------------------------------------------------------------------------------------
! General variables for the homogenization at a  material point
  real(pReal),   dimension(:,:,:,:),   allocatable, public :: &
    materialpoint_F0, &                                                                             !< def grad of IP at start of FE increment
    materialpoint_F, &                                                                              !< def grad of IP to be reached at end of FE increment
    materialpoint_P                                                                                 !< first P--K stress of IP
  real(pReal),   dimension(:,:,:,:,:,:), allocatable, public ::  &
    materialpoint_dPdF                                                                              !< tangent of first P--K stress at IP
  real(pReal),   dimension(:,:,:),       allocatable, public :: &
    materialpoint_results                                                                           !< results array of material point
  integer,                                            public, protected  :: &
    materialpoint_sizeResults, &
    thermal_maxSizePostResults, &
    damage_maxSizePostResults

  real(pReal),   dimension(:,:,:,:),     allocatable :: &
    materialpoint_subF0, &                                                                          !< def grad of IP at beginning of homogenization increment
    materialpoint_subF                                                                              !< def grad of IP to be reached at end of homog inc
  real(pReal),   dimension(:,:),         allocatable :: &
    materialpoint_subFrac, &
    materialpoint_subStep, &
    materialpoint_subdt
  logical,       dimension(:,:),         allocatable :: &
    materialpoint_requested, &
    materialpoint_converged
  logical,       dimension(:,:,:),       allocatable :: &
    materialpoint_doneAndHappy
    
  interface

    module subroutine mech_none_init
    end subroutine mech_none_init
    
    module subroutine mech_isostrain_init
    end subroutine mech_isostrain_init
    
    module subroutine mech_RGC_init
    end subroutine mech_RGC_init
    
    
    module subroutine mech_isostrain_partitionDeformation(F,avgF)
      real(pReal),   dimension (:,:,:), intent(out) :: F                                            !< partitioned deformation gradient
      real(pReal),   dimension (3,3),   intent(in)  :: avgF                                         !< average deformation gradient at material point
    end subroutine mech_isostrain_partitionDeformation
    
    module subroutine mech_RGC_partitionDeformation(F,avgF,instance,of)
      real(pReal),   dimension (:,:,:), intent(out) :: F                                            !< partitioned deformation gradient
      real(pReal),   dimension (3,3),   intent(in)  :: avgF                                         !< average deformation gradient at material point
      integer,                          intent(in)  :: &
        instance, &
        of
    end subroutine mech_RGC_partitionDeformation
    
    
    module subroutine mech_isostrain_averageStressAndItsTangent(avgP,dAvgPdAvgF,P,dPdF,instance)
      real(pReal),   dimension (3,3),       intent(out) :: avgP                                     !< average stress at material point
      real(pReal),   dimension (3,3,3,3),   intent(out) :: dAvgPdAvgF                               !< average stiffness at material point
  
      real(pReal),   dimension (:,:,:),     intent(in)  :: P                                        !< partitioned stresses
      real(pReal),   dimension (:,:,:,:,:), intent(in)  :: dPdF                                     !< partitioned stiffnesses
      integer,                              intent(in)  :: instance 
    end subroutine mech_isostrain_averageStressAndItsTangent
    
    module subroutine mech_RGC_averageStressAndItsTangent(avgP,dAvgPdAvgF,P,dPdF,instance)
      real(pReal),   dimension (3,3),       intent(out) :: avgP                                     !< average stress at material point
      real(pReal),   dimension (3,3,3,3),   intent(out) :: dAvgPdAvgF                               !< average stiffness at material point
  
      real(pReal),   dimension (:,:,:),     intent(in)  :: P                                        !< partitioned stresses
      real(pReal),   dimension (:,:,:,:,:), intent(in)  :: dPdF                                     !< partitioned stiffnesses
      integer,                              intent(in)  :: instance 
    end subroutine mech_RGC_averageStressAndItsTangent
    
    
    module function mech_RGC_updateState(P,F,F0,avgF,dt,dPdF,ip,el)
      logical, dimension(2) :: mech_RGC_updateState
      real(pReal), dimension(:,:,:),     intent(in)    :: & 
        P,&                                                                                         !< partitioned stresses
        F,&                                                                                         !< partitioned deformation gradients
        F0                                                                                          !< partitioned initial deformation gradients
     real(pReal), dimension(:,:,:,:,:), intent(in) :: dPdF                                          !< partitioned stiffnesses
     real(pReal), dimension(3,3),       intent(in) :: avgF                                          !< average F
     real(pReal),                       intent(in) :: dt                                            !< time increment
     integer,                           intent(in) :: &
      ip, &                                                                                         !< integration point number
      el                                                                                            !< element number
    end function mech_RGC_updateState


    module subroutine mech_RGC_results(instance,group)
      integer,          intent(in) :: instance                                                      !< homogenization instance
      character(len=*), intent(in) :: group                                                         !< group name in HDF5 file
    end subroutine mech_RGC_results
   
  end interface

  public ::  &
    homogenization_init, &
    materialpoint_stressAndItsTangent, &
    materialpoint_postResults, &
    homogenization_results

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!--------------------------------------------------------------------------------------------------
subroutine homogenization_init

  integer, parameter :: FILEUNIT = 200
  integer :: e,i,p
  integer, dimension(:,:), pointer :: thisSize
  integer, dimension(:)  , pointer :: thisNoutput
  character(len=64), dimension(:,:), pointer :: thisOutput
  character(len=32) :: outputName                                                                   !< name of output, intermediate fix until HDF5 output is ready
  logical :: valid

  if (any(homogenization_type == HOMOGENIZATION_NONE_ID))      call mech_none_init
  if (any(homogenization_type == HOMOGENIZATION_ISOSTRAIN_ID)) call mech_isostrain_init
  if (any(homogenization_type == HOMOGENIZATION_RGC_ID))       call mech_RGC_init

  if (any(thermal_type == THERMAL_isothermal_ID)) call thermal_isothermal_init
  if (any(thermal_type == THERMAL_adiabatic_ID))  call thermal_adiabatic_init
  if (any(thermal_type == THERMAL_conduction_ID)) call thermal_conduction_init

  if (any(damage_type == DAMAGE_none_ID))      call damage_none_init
  if (any(damage_type == DAMAGE_local_ID))     call damage_local_init
  if (any(damage_type == DAMAGE_nonlocal_ID))  call damage_nonlocal_init

!--------------------------------------------------------------------------------------------------
! write description file for homogenization output
  mainProcess: if (worldrank == 0) then
    call IO_write_jobFile(FILEUNIT,'outputHomogenization')
    do p = 1,size(config_homogenization)
      if (any(material_homogenizationAt == p)) then
        write(FILEUNIT,'(/,a,/)')  '['//trim(config_name_homogenization(p))//']'
        write(FILEUNIT,'(a)') '(type) n/a'
        write(FILEUNIT,'(a,i4)') '(ngrains)'//char(9),homogenization_Ngrains(p)
        
        i = thermal_typeInstance(p)                                                                 ! which instance of this thermal type
        valid = .true.                                                                              ! assume valid
        select case(thermal_type(p))                                                                ! split per thermal type
          case (THERMAL_isothermal_ID)
            outputName = THERMAL_isothermal_label
            thisNoutput => null()
            thisOutput => null()
            thisSize   => null()
          case (THERMAL_adiabatic_ID)
            outputName = THERMAL_adiabatic_label
            thisNoutput => thermal_adiabatic_Noutput
            thisOutput => thermal_adiabatic_output
            thisSize   => thermal_adiabatic_sizePostResult
          case (THERMAL_conduction_ID)
            outputName = THERMAL_conduction_label
            thisNoutput => thermal_conduction_Noutput
            thisOutput => thermal_conduction_output
            thisSize   => thermal_conduction_sizePostResult
          case default
            valid = .false.
        end select
        if (valid) then
          write(FILEUNIT,'(a)') '(thermal)'//char(9)//trim(outputName)
          if (thermal_type(p) /= THERMAL_isothermal_ID) then
            do e = 1,thisNoutput(i)
              write(FILEUNIT,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
            enddo
          endif
        endif
        
        i = damage_typeInstance(p)                                                                  ! which instance of this damage type
        valid = .true.                                                                              ! assume valid
        select case(damage_type(p))                                                                 ! split per damage type
          case (DAMAGE_none_ID)
            outputName = DAMAGE_none_label
            thisNoutput => null()
            thisOutput => null()
            thisSize   => null()
          case (DAMAGE_local_ID)
            outputName = DAMAGE_local_label
            thisNoutput => damage_local_Noutput
            thisOutput => damage_local_output
            thisSize   => damage_local_sizePostResult
          case (DAMAGE_nonlocal_ID)
            outputName = DAMAGE_nonlocal_label
            thisNoutput => damage_nonlocal_Noutput
            thisOutput => damage_nonlocal_output
            thisSize   => damage_nonlocal_sizePostResult
          case default
            valid = .false.
        end select
        if (valid) then
          write(FILEUNIT,'(a)') '(damage)'//char(9)//trim(outputName)
          if (damage_type(p) /= DAMAGE_none_ID) then
            do e = 1,thisNoutput(i)
              write(FILEUNIT,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
            enddo
          endif
        endif
      endif
    enddo
    close(FILEUNIT)
  endif mainProcess

  call config_deallocate('material.config/homogenization')

!--------------------------------------------------------------------------------------------------
! allocate and initialize global variables
  allocate(materialpoint_dPdF(3,3,3,3,discretization_nIP,discretization_nElem),       source=0.0_pReal)
  allocate(materialpoint_F0(3,3,discretization_nIP,discretization_nElem),             source=0.0_pReal)
  materialpoint_F0 = spread(spread(math_I3,3,discretization_nIP),4,discretization_nElem)            ! initialize to identity
  allocate(materialpoint_F(3,3,discretization_nIP,discretization_nElem),              source=0.0_pReal)
  materialpoint_F = materialpoint_F0                                                                ! initialize to identity
  allocate(materialpoint_subF0(3,3,discretization_nIP,discretization_nElem),          source=0.0_pReal)
  allocate(materialpoint_subF(3,3,discretization_nIP,discretization_nElem),           source=0.0_pReal)
  allocate(materialpoint_P(3,3,discretization_nIP,discretization_nElem),              source=0.0_pReal)
  allocate(materialpoint_subFrac(discretization_nIP,discretization_nElem),            source=0.0_pReal)
  allocate(materialpoint_subStep(discretization_nIP,discretization_nElem),            source=0.0_pReal)
  allocate(materialpoint_subdt(discretization_nIP,discretization_nElem),              source=0.0_pReal)
  allocate(materialpoint_requested(discretization_nIP,discretization_nElem),          source=.false.)
  allocate(materialpoint_converged(discretization_nIP,discretization_nElem),          source=.true.)
  allocate(materialpoint_doneAndHappy(2,discretization_nIP,discretization_nElem),     source=.true.)

!--------------------------------------------------------------------------------------------------
! allocate and initialize global state and postresutls variables
  thermal_maxSizePostResults        = 0
  damage_maxSizePostResults         = 0
  do p = 1,size(config_homogenization)
    thermal_maxSizePostResults        = max(thermal_maxSizePostResults,       thermalState     (p)%sizePostResults)
    damage_maxSizePostResults         = max(damage_maxSizePostResults        ,damageState      (p)%sizePostResults)
  enddo

  materialpoint_sizeResults = 1 &                                                                   ! grain count
                            + 1 + thermal_maxSizePostResults        &
                                + damage_maxSizePostResults         &
                            + homogenization_maxNgrains * (1 &     ! crystallite size
                                                         + 1 + constitutive_plasticity_maxSizePostResults & ! constitutive size & constitutive results
                                                             + constitutive_source_maxSizePostResults)
  allocate(materialpoint_results(materialpoint_sizeResults,discretization_nIP,discretization_nElem))

  write(6,'(/,a)')   ' <<<+-  homogenization init  -+>>>'

  if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0) then
    write(6,'(a32,1x,7(i8,1x))')   'materialpoint_dPdF:             ', shape(materialpoint_dPdF)
    write(6,'(a32,1x,7(i8,1x))')   'materialpoint_F0:               ', shape(materialpoint_F0)
    write(6,'(a32,1x,7(i8,1x))')   'materialpoint_F:                ', shape(materialpoint_F)
    write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subF0:            ', shape(materialpoint_subF0)
    write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subF:             ', shape(materialpoint_subF)
    write(6,'(a32,1x,7(i8,1x))')   'materialpoint_P:                ', shape(materialpoint_P)
    write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subFrac:          ', shape(materialpoint_subFrac)
    write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subStep:          ', shape(materialpoint_subStep)
    write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subdt:            ', shape(materialpoint_subdt)
    write(6,'(a32,1x,7(i8,1x))')   'materialpoint_requested:        ', shape(materialpoint_requested)
    write(6,'(a32,1x,7(i8,1x))')   'materialpoint_converged:        ', shape(materialpoint_converged)
    write(6,'(a32,1x,7(i8,1x),/)') 'materialpoint_doneAndHappy:     ', shape(materialpoint_doneAndHappy)
  endif
  flush(6)

  if (debug_g < 1 .or. debug_g > homogenization_Ngrains(material_homogenizationAt(debug_e))) &
    call IO_error(602,ext_msg='constituent', el=debug_e, g=debug_g)

end subroutine homogenization_init


!--------------------------------------------------------------------------------------------------
!> @brief  parallelized calculation of stress and corresponding tangent at material points
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_stressAndItsTangent(updateJaco,dt)

  real(pReal), intent(in) :: dt                                                                     !< time increment
  logical,     intent(in) :: updateJaco                                                             !< initiating Jacobian update
  integer :: &
    NiterationHomog, &
    NiterationMPstate, &
    g, &                                                                                            !< grain number
    i, &                                                                                            !< integration point number
    e, &                                                                                            !< element number
    mySource, &
    myNgrains

#ifdef DEBUG
  if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0) then
    write(6,'(/a,i5,1x,i2)') '<< HOMOG >> Material Point start at el ip ', debug_e, debug_i

    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< HOMOG >> F0', &
                                    transpose(materialpoint_F0(1:3,1:3,debug_i,debug_e))
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< HOMOG >> F', &
                                    transpose(materialpoint_F(1:3,1:3,debug_i,debug_e))
  endif
#endif

!--------------------------------------------------------------------------------------------------
! initialize restoration points of ...
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    myNgrains = homogenization_Ngrains(material_homogenizationAt(e))
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e);
      do g = 1,myNgrains

        plasticState    (material_phaseAt(g,e))%partionedState0(:,material_phasememberAt(g,i,e)) = &
        plasticState    (material_phaseAt(g,e))%state0(         :,material_phasememberAt(g,i,e))
        do mySource = 1, phase_Nsources(material_phaseAt(g,e))
          sourceState(material_phaseAt(g,e))%p(mySource)%partionedState0(:,material_phasememberAt(g,i,e)) = &
          sourceState(material_phaseAt(g,e))%p(mySource)%state0(         :,material_phasememberAt(g,i,e))
        enddo

        crystallite_partionedFp0(1:3,1:3,g,i,e) = crystallite_Fp0(1:3,1:3,g,i,e)
        crystallite_partionedLp0(1:3,1:3,g,i,e) = crystallite_Lp0(1:3,1:3,g,i,e)
        crystallite_partionedFi0(1:3,1:3,g,i,e) = crystallite_Fi0(1:3,1:3,g,i,e)
        crystallite_partionedLi0(1:3,1:3,g,i,e) = crystallite_Li0(1:3,1:3,g,i,e)
        crystallite_partionedF0(1:3,1:3,g,i,e)  = crystallite_F0(1:3,1:3,g,i,e)
        crystallite_partionedS0(1:3,1:3,g,i,e)  = crystallite_S0(1:3,1:3,g,i,e)

      enddo


      materialpoint_subF0(1:3,1:3,i,e) = materialpoint_F0(1:3,1:3,i,e)
      materialpoint_subFrac(i,e) = 0.0_pReal
      materialpoint_subStep(i,e) = 1.0_pReal/subStepSizeHomog                                       ! <<added to adopt flexibility in cutback size>>
      materialpoint_converged(i,e) = .false.                                                        ! pretend failed step of twice the required size
      materialpoint_requested(i,e) = .true.                                                         ! everybody requires calculation
    
      if (homogState(material_homogenizationAt(e))%sizeState > 0) &
          homogState(material_homogenizationAt(e))%subState0(:,mappingHomogenization(1,i,e)) = &
          homogState(material_homogenizationAt(e))%State0(   :,mappingHomogenization(1,i,e))        ! ...internal homogenization state

      if (thermalState(material_homogenizationAt(e))%sizeState > 0) &
          thermalState(material_homogenizationAt(e))%subState0(:,mappingHomogenization(1,i,e)) = &
          thermalState(material_homogenizationAt(e))%State0(   :,mappingHomogenization(1,i,e))      ! ...internal thermal state
        
      if (damageState(material_homogenizationAt(e))%sizeState > 0) &
          damageState(material_homogenizationAt(e))%subState0(:,mappingHomogenization(1,i,e)) = &
          damageState(material_homogenizationAt(e))%State0(   :,mappingHomogenization(1,i,e))       ! ...internal damage state
    enddo
  enddo
  
  NiterationHomog = 0

  cutBackLooping: do while (.not. terminallyIll .and. &
       any(materialpoint_subStep(:,FEsolving_execELem(1):FEsolving_execElem(2)) > subStepMinHomog))

    !$OMP PARALLEL DO PRIVATE(myNgrains)
    elementLooping1: do e = FEsolving_execElem(1),FEsolving_execElem(2)
      myNgrains = homogenization_Ngrains(material_homogenizationAt(e))
      IpLooping1: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)

        converged: if (materialpoint_converged(i,e)) then
#ifdef DEBUG
          if (iand(debug_level(debug_homogenization), debug_levelExtensive) /= 0 &
             .and. ((e == debug_e .and. i == debug_i) &
                    .or. .not. iand(debug_level(debug_homogenization),debug_levelSelective) /= 0)) then
            write(6,'(a,1x,f12.8,1x,a,1x,f12.8,1x,a,i8,1x,i2/)') '<< HOMOG >> winding forward from', &
              materialpoint_subFrac(i,e), 'to current materialpoint_subFrac', &
              materialpoint_subFrac(i,e)+materialpoint_subStep(i,e),'in materialpoint_stressAndItsTangent at el ip',e,i
          endif
#endif

!---------------------------------------------------------------------------------------------------
! calculate new subStep and new subFrac
          materialpoint_subFrac(i,e) = materialpoint_subFrac(i,e) + materialpoint_subStep(i,e)
          materialpoint_subStep(i,e) = min(1.0_pReal-materialpoint_subFrac(i,e), &
                                           stepIncreaseHomog*materialpoint_subStep(i,e))            ! introduce flexibility for step increase/acceleration

          steppingNeeded: if (materialpoint_subStep(i,e) > subStepMinHomog) then

            ! wind forward grain starting point of...
            crystallite_partionedF0  (1:3,1:3,1:myNgrains,i,e) =  &
               crystallite_partionedF(1:3,1:3,1:myNgrains,i,e)

            crystallite_partionedFp0 (1:3,1:3,1:myNgrains,i,e) = &
              crystallite_Fp         (1:3,1:3,1:myNgrains,i,e)

            crystallite_partionedLp0 (1:3,1:3,1:myNgrains,i,e) = &
              crystallite_Lp         (1:3,1:3,1:myNgrains,i,e)

            crystallite_partionedFi0 (1:3,1:3,1:myNgrains,i,e) = &
              crystallite_Fi         (1:3,1:3,1:myNgrains,i,e)

            crystallite_partionedLi0 (1:3,1:3,1:myNgrains,i,e) = &
              crystallite_Li         (1:3,1:3,1:myNgrains,i,e)

            crystallite_partionedS0  (1:3,1:3,1:myNgrains,i,e) = &
              crystallite_S          (1:3,1:3,1:myNgrains,i,e)

            do g = 1,myNgrains
              plasticState    (material_phaseAt(g,e))%partionedState0(:,material_phasememberAt(g,i,e)) = &
              plasticState    (material_phaseAt(g,e))%state          (:,material_phasememberAt(g,i,e))
              do mySource = 1, phase_Nsources(material_phaseAt(g,e))
                sourceState(material_phaseAt(g,e))%p(mySource)%partionedState0(:,material_phasememberAt(g,i,e)) = &
                sourceState(material_phaseAt(g,e))%p(mySource)%state          (:,material_phasememberAt(g,i,e))
              enddo
            enddo

            if(homogState(material_homogenizationAt(e))%sizeState > 0) &
                homogState(material_homogenizationAt(e))%subState0(:,mappingHomogenization(1,i,e)) = &
                homogState(material_homogenizationAt(e))%State    (:,mappingHomogenization(1,i,e))
            if(thermalState(material_homogenizationAt(e))%sizeState > 0) &
                thermalState(material_homogenizationAt(e))%subState0(:,mappingHomogenization(1,i,e)) = &
                thermalState(material_homogenizationAt(e))%State    (:,mappingHomogenization(1,i,e))
            if(damageState(material_homogenizationAt(e))%sizeState > 0) &
                damageState(material_homogenizationAt(e))%subState0(:,mappingHomogenization(1,i,e)) = &
                damageState(material_homogenizationAt(e))%State    (:,mappingHomogenization(1,i,e))
                
            materialpoint_subF0(1:3,1:3,i,e) = materialpoint_subF(1:3,1:3,i,e)

          endif steppingNeeded

        else converged
          if ( (myNgrains == 1 .and. materialpoint_subStep(i,e) <= 1.0 ) .or. &                     ! single grain already tried internal subStepping in crystallite
               subStepSizeHomog * materialpoint_subStep(i,e) <=  subStepMinHomog ) then             ! would require too small subStep
                                                                                                    ! cutback makes no sense
            !$OMP FLUSH(terminallyIll)
            if (.not. terminallyIll) then                                                           ! so first signals terminally ill...
              !$OMP CRITICAL (write2out)
                write(6,*) 'Integration point ', i,' at element ', e, ' terminally ill'
              !$OMP END CRITICAL (write2out)
            endif
            !$OMP CRITICAL (setTerminallyIll)
              terminallyIll = .true.                                                                ! ...and kills all others
            !$OMP END CRITICAL (setTerminallyIll)
          else                                                                                      ! cutback makes sense
            materialpoint_subStep(i,e) = subStepSizeHomog * materialpoint_subStep(i,e)              ! crystallite had severe trouble, so do a significant cutback

#ifdef DEBUG
            if (iand(debug_level(debug_homogenization), debug_levelExtensive) /= 0 &
               .and. ((e == debug_e .and. i == debug_i) &
                     .or. .not. iand(debug_level(debug_homogenization), debug_levelSelective) /= 0)) then
              write(6,'(a,1x,f12.8,a,i8,1x,i2/)') &
                '<< HOMOG >> cutback step in materialpoint_stressAndItsTangent with new materialpoint_subStep:',&
                materialpoint_subStep(i,e),' at el ip',e,i
            endif
#endif

!--------------------------------------------------------------------------------------------------
! restore...
            if (materialpoint_subStep(i,e) < 1.0_pReal) then                                        ! protect against fake cutback from \Delta t = 2 to 1. Maybe that "trick" is not necessary anymore at all? I.e. start with \Delta t = 1
              crystallite_Lp(1:3,1:3,1:myNgrains,i,e) = &
                crystallite_partionedLp0(1:3,1:3,1:myNgrains,i,e)
              crystallite_Li(1:3,1:3,1:myNgrains,i,e) = &
                crystallite_partionedLi0(1:3,1:3,1:myNgrains,i,e)
            endif                                                                                   ! maybe protecting everything from overwriting (not only L) makes even more sense
            crystallite_Fp(1:3,1:3,1:myNgrains,i,e) = &
              crystallite_partionedFp0(1:3,1:3,1:myNgrains,i,e)
            crystallite_Fi(1:3,1:3,1:myNgrains,i,e) = &
              crystallite_partionedFi0(1:3,1:3,1:myNgrains,i,e)
            crystallite_S(1:3,1:3,1:myNgrains,i,e) = &
               crystallite_partionedS0(1:3,1:3,1:myNgrains,i,e)
            do g = 1, myNgrains
              plasticState    (material_phaseAt(g,e))%state(          :,material_phasememberAt(g,i,e)) = &
              plasticState    (material_phaseAt(g,e))%partionedState0(:,material_phasememberAt(g,i,e))
              do mySource = 1, phase_Nsources(material_phaseAt(g,e))
                sourceState(material_phaseAt(g,e))%p(mySource)%state(          :,material_phasememberAt(g,i,e)) = &
                sourceState(material_phaseAt(g,e))%p(mySource)%partionedState0(:,material_phasememberAt(g,i,e))
              enddo
            enddo
            if(homogState(material_homogenizationAt(e))%sizeState > 0) &
                homogState(material_homogenizationAt(e))%State(    :,mappingHomogenization(1,i,e)) = &
                homogState(material_homogenizationAt(e))%subState0(:,mappingHomogenization(1,i,e))
            if(thermalState(material_homogenizationAt(e))%sizeState > 0) &
                thermalState(material_homogenizationAt(e))%State(    :,mappingHomogenization(1,i,e)) = &
                thermalState(material_homogenizationAt(e))%subState0(:,mappingHomogenization(1,i,e))
            if(damageState(material_homogenizationAt(e))%sizeState > 0) &
                damageState(material_homogenizationAt(e))%State(    :,mappingHomogenization(1,i,e)) = &
                damageState(material_homogenizationAt(e))%subState0(:,mappingHomogenization(1,i,e))
          endif
        endif converged

        if (materialpoint_subStep(i,e) > subStepMinHomog) then
          materialpoint_requested(i,e) = .true.
          materialpoint_subF(1:3,1:3,i,e) = materialpoint_subF0(1:3,1:3,i,e) &
                                          + materialpoint_subStep(i,e) * (materialpoint_F(1:3,1:3,i,e) &
                                          - materialpoint_F0(1:3,1:3,i,e))
          materialpoint_subdt(i,e) = materialpoint_subStep(i,e) * dt
          materialpoint_doneAndHappy(1:2,i,e) = [.false.,.true.]
        endif
      enddo IpLooping1
    enddo elementLooping1
    !$OMP END PARALLEL DO

    NiterationMPstate = 0

    convergenceLooping: do while (.not. terminallyIll .and. &
              any(            materialpoint_requested(:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                  .and. .not. materialpoint_doneAndHappy(1,:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                 ) .and. &
              NiterationMPstate < nMPstate)
      NiterationMPstate = NiterationMPstate + 1

!--------------------------------------------------------------------------------------------------
! deformation partitioning
! based on materialpoint_subF0,.._subF,crystallite_partionedF0, and homogenization_state,
! results in crystallite_partionedF
      !$OMP PARALLEL DO PRIVATE(myNgrains)
      elementLooping2: do e = FEsolving_execElem(1),FEsolving_execElem(2)
        myNgrains = homogenization_Ngrains(material_homogenizationAt(e))
        IpLooping2: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
          if (      materialpoint_requested(i,e) .and. &                                            ! process requested but...
              .not. materialpoint_doneAndHappy(1,i,e)) then                                         ! ...not yet done material points
            call partitionDeformation(i,e)                                                          ! partition deformation onto constituents
            crystallite_dt(1:myNgrains,i,e) = materialpoint_subdt(i,e)                              ! propagate materialpoint dt to grains
            crystallite_requested(1:myNgrains,i,e) = .true.                                         ! request calculation for constituents
          else
            crystallite_requested(1:myNgrains,i,e) = .false.                                        ! calculation for constituents not required anymore
          endif
        enddo IpLooping2
      enddo elementLooping2
      !$OMP END PARALLEL DO

!--------------------------------------------------------------------------------------------------
! crystallite integration
! based on crystallite_partionedF0,.._partionedF
! incrementing by crystallite_dt
    
      materialpoint_converged = crystallite_stress() !ToDo: MD not sure if that is the best logic

!--------------------------------------------------------------------------------------------------
! state update
     !$OMP PARALLEL DO
      elementLooping3: do e = FEsolving_execElem(1),FEsolving_execElem(2)
        IpLooping3: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
          if (      materialpoint_requested(i,e) .and. &
              .not. materialpoint_doneAndHappy(1,i,e)) then
            if (.not. materialpoint_converged(i,e)) then
              materialpoint_doneAndHappy(1:2,i,e) = [.true.,.false.]
            else
              materialpoint_doneAndHappy(1:2,i,e) = updateState(i,e)
              materialpoint_converged(i,e) = all(materialpoint_doneAndHappy(1:2,i,e))               ! converged if done and happy
            endif
          endif
        enddo IpLooping3
      enddo elementLooping3
      !$OMP END PARALLEL DO
 
    enddo convergenceLooping
 
    NiterationHomog = NiterationHomog + 1
 
  enddo cutBackLooping
  
  if(updateJaco) call crystallite_stressTangent
 
  if (.not. terminallyIll ) then
    call crystallite_orientations()                                                                 ! calculate crystal orientations
    !$OMP PARALLEL DO
    elementLooping4: do e = FEsolving_execElem(1),FEsolving_execElem(2)
      IpLooping4: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
        call averageStressAndItsTangent(i,e)
      enddo IpLooping4
    enddo elementLooping4
    !$OMP END PARALLEL DO
  else
    write(6,'(/,a,/)') '<< HOMOG >> Material Point terminally ill'
  endif

end subroutine materialpoint_stressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief parallelized calculation of result array at material points
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_postResults

  integer :: &
    thePos, &
    theSize, &
    myNgrains, &
    myCrystallite, &
    g, &                                                                                            !< grain number
    i, &                                                                                            !< integration point number
    e                                                                                               !< element number

  !$OMP PARALLEL DO PRIVATE(myNgrains,myCrystallite,thePos,theSize)
  elementLooping: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    myNgrains = homogenization_Ngrains(material_homogenizationAt(e))
    myCrystallite = microstructure_crystallite(discretization_microstructureAt(e))
    IpLooping: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      thePos = 0

      theSize = thermalState     (material_homogenizationAt(e))%sizePostResults &
              + damageState      (material_homogenizationAt(e))%sizePostResults
      materialpoint_results(thePos+1,i,e) = real(theSize,pReal)                                    ! tell size of homogenization results
      thePos = thePos + 1

      if (theSize > 0) then                                                                         ! any homogenization results to mention?
        materialpoint_results(thePos+1:thePos+theSize,i,e) = postResults(i,e)
        thePos = thePos + theSize
      endif

      materialpoint_results(thePos+1,i,e) = real(myNgrains,pReal)                                   ! tell number of grains at materialpoint
      thePos = thePos + 1

      grainLooping :do g = 1,myNgrains
        theSize = 1 + &
                  1 + plasticState    (material_phaseAt(g,e))%sizePostResults + &
                      sum(sourceState(material_phaseAt(g,e))%p(:)%sizePostResults)
        materialpoint_results(thePos+1:thePos+theSize,i,e) = crystallite_postResults(g,i,e)        ! tell crystallite results
        thePos = thePos + theSize
      enddo grainLooping
    enddo IpLooping
  enddo elementLooping
 !$OMP END PARALLEL DO

end subroutine materialpoint_postResults


!--------------------------------------------------------------------------------------------------
!> @brief  partition material point def grad onto constituents
!--------------------------------------------------------------------------------------------------
subroutine partitionDeformation(ip,el)

 integer, intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number

 chosenHomogenization: select case(homogenization_type(material_homogenizationAt(el)))

   case (HOMOGENIZATION_NONE_ID) chosenHomogenization
     crystallite_partionedF(1:3,1:3,1,ip,el) = materialpoint_subF(1:3,1:3,ip,el)

   case (HOMOGENIZATION_ISOSTRAIN_ID) chosenHomogenization
     call mech_isostrain_partitionDeformation(&
                          crystallite_partionedF(1:3,1:3,1:homogenization_Ngrains(material_homogenizationAt(el)),ip,el), &
                          materialpoint_subF(1:3,1:3,ip,el))

   case (HOMOGENIZATION_RGC_ID) chosenHomogenization
     call mech_RGC_partitionDeformation(&
                         crystallite_partionedF(1:3,1:3,1:homogenization_Ngrains(material_homogenizationAt(el)),ip,el), &
                         materialpoint_subF(1:3,1:3,ip,el),&
                         ip, &
                         el)
 end select chosenHomogenization

end subroutine partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief update the internal state of the homogenization scheme and tell whether "done" and
!> "happy" with result
!--------------------------------------------------------------------------------------------------
function updateState(ip,el)

 integer, intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number
 logical, dimension(2) :: updateState

 updateState = .true.
 chosenHomogenization: select case(homogenization_type(material_homogenizationAt(el)))
   case (HOMOGENIZATION_RGC_ID) chosenHomogenization
     updateState = &
       updateState .and. &
        mech_RGC_updateState(crystallite_P(1:3,1:3,1:homogenization_Ngrains(material_homogenizationAt(el)),ip,el), &
                             crystallite_partionedF(1:3,1:3,1:homogenization_Ngrains(material_homogenizationAt(el)),ip,el), &
                             crystallite_partionedF0(1:3,1:3,1:homogenization_Ngrains(material_homogenizationAt(el)),ip,el),&
                             materialpoint_subF(1:3,1:3,ip,el),&
                             materialpoint_subdt(ip,el), &
                             crystallite_dPdF(1:3,1:3,1:3,1:3,1:homogenization_Ngrains(material_homogenizationAt(el)),ip,el), &
                             ip, &
                             el)
 end select chosenHomogenization

 chosenThermal: select case (thermal_type(material_homogenizationAt(el)))
   case (THERMAL_adiabatic_ID) chosenThermal
     updateState = &
       updateState .and. &
       thermal_adiabatic_updateState(materialpoint_subdt(ip,el), &
                                     ip, &
                                     el)
 end select chosenThermal

 chosenDamage: select case (damage_type(material_homogenizationAt(el)))
   case (DAMAGE_local_ID) chosenDamage
     updateState = &
       updateState .and. &
       damage_local_updateState(materialpoint_subdt(ip,el), &
                                ip, &
                                el)
 end select chosenDamage

end function updateState


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities
!--------------------------------------------------------------------------------------------------
subroutine averageStressAndItsTangent(ip,el)

 integer, intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number

 chosenHomogenization: select case(homogenization_type(material_homogenizationAt(el)))
   case (HOMOGENIZATION_NONE_ID) chosenHomogenization
       materialpoint_P(1:3,1:3,ip,el)            = crystallite_P(1:3,1:3,1,ip,el)
       materialpoint_dPdF(1:3,1:3,1:3,1:3,ip,el) = crystallite_dPdF(1:3,1:3,1:3,1:3,1,ip,el)

   case (HOMOGENIZATION_ISOSTRAIN_ID) chosenHomogenization
     call mech_isostrain_averageStressAndItsTangent(&
       materialpoint_P(1:3,1:3,ip,el), &
       materialpoint_dPdF(1:3,1:3,1:3,1:3,ip,el),&
       crystallite_P(1:3,1:3,1:homogenization_Ngrains(material_homogenizationAt(el)),ip,el), &
       crystallite_dPdF(1:3,1:3,1:3,1:3,1:homogenization_Ngrains(material_homogenizationAt(el)),ip,el), &
       homogenization_typeInstance(material_homogenizationAt(el)))

   case (HOMOGENIZATION_RGC_ID) chosenHomogenization
     call mech_RGC_averageStressAndItsTangent(&
       materialpoint_P(1:3,1:3,ip,el), &
       materialpoint_dPdF(1:3,1:3,1:3,1:3,ip,el),&
       crystallite_P(1:3,1:3,1:homogenization_Ngrains(material_homogenizationAt(el)),ip,el), &
       crystallite_dPdF(1:3,1:3,1:3,1:3,1:homogenization_Ngrains(material_homogenizationAt(el)),ip,el), &
       homogenization_typeInstance(material_homogenizationAt(el)))
 end select chosenHomogenization

end subroutine averageStressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion. call only,
!> if homogenization_sizePostResults(i,e) > 0 !!
!--------------------------------------------------------------------------------------------------
function postResults(ip,el)

 integer, intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number
 real(pReal), dimension(  thermalState     (material_homogenizationAt(el))%sizePostResults &
                        + damageState      (material_homogenizationAt(el))%sizePostResults) :: &
   postResults
 integer :: &
   startPos, endPos ,&
   homog


 postResults = 0.0_pReal
 startPos = 1
 endPos   = thermalState(material_homogenizationAt(el))%sizePostResults
 chosenThermal: select case (thermal_type(material_homogenizationAt(el)))

   case (THERMAL_adiabatic_ID) chosenThermal
     homog = material_homogenizationAt(el)
     postResults(startPos:endPos) = &
       thermal_adiabatic_postResults(homog,thermal_typeInstance(homog),thermalMapping(homog)%p(ip,el))
   case (THERMAL_conduction_ID) chosenThermal
     homog = material_homogenizationAt(el)
     postResults(startPos:endPos) = &
       thermal_conduction_postResults(homog,thermal_typeInstance(homog),thermalMapping(homog)%p(ip,el))

 end select chosenThermal

 startPos = endPos + 1
 endPos   = endPos + damageState(material_homogenizationAt(el))%sizePostResults
 chosenDamage: select case (damage_type(material_homogenizationAt(el)))

   case (DAMAGE_local_ID) chosenDamage
     postResults(startPos:endPos) = damage_local_postResults(ip, el)
   case (DAMAGE_nonlocal_ID) chosenDamage
     postResults(startPos:endPos) = damage_nonlocal_postResults(ip, el)
     
 end select chosenDamage

end function postResults


!--------------------------------------------------------------------------------------------------
!> @brief writes homogenization results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_results
#if defined(PETSc) || defined(DAMASK_HDF5)
  use material, only: &
    material_homogenization_type => homogenization_type
    
  integer :: p
  character(len=256) :: group
  
  !real(pReal), dimension(:,:,:), allocatable :: temp
                                             
  do p=1,size(config_name_homogenization)
    group = trim('current/materialpoint')//'/'//trim(config_name_homogenization(p))
    call HDF5_closeGroup(results_addGroup(group))
    
    group = trim(group)//'/mech'
    
    call HDF5_closeGroup(results_addGroup(group))  
    select case(material_homogenization_type(p))
      case(HOMOGENIZATION_rgc_ID)
        call mech_RGC_results(homogenization_typeInstance(p),group)
    end select
    
    group = trim('current/materialpoint')//'/'//trim(config_name_homogenization(p))//'/generic'
    call HDF5_closeGroup(results_addGroup(group))
    
    !temp = reshape(materialpoint_F,[3,3,discretization_nIP*discretization_nElem])
    !call results_writeDataset(group,temp,'F',&
    !                          'deformation gradient','1')  
    !temp = reshape(materialpoint_P,[3,3,discretization_nIP*discretization_nElem])
    !call results_writeDataset(group,temp,'P',&
    !                          '1st Piola-Kirchoff stress','Pa')  

 enddo   
#endif
end subroutine homogenization_results

end module homogenization
