!----------------------------------------------------------------------------------------------------
!> @brief internal microstructure state for all damage sources and kinematics constitutive models
!----------------------------------------------------------------------------------------------------
submodule(constitutive) constitutive_damage
  enum, bind(c); enumerator :: &
    DAMAGE_UNDEFINED_ID, &
    DAMAGE_ISOBRITTLE_ID, &
    DAMAGE_ISODUCTILE_ID, &
    DAMAGE_ANISOBRITTLE_ID, &
    DAMAGE_ANISODUCTILE_ID
  end enum


  type :: tDataContainer
    real(pReal), dimension(:), allocatable :: phi, d_phi_d_dot_phi
  end type tDataContainer

  integer(kind(DAMAGE_UNDEFINED_ID)),     dimension(:,:), allocatable :: &
    phase_source                                                                                !< active sources mechanisms of each phase

  type(tDataContainer), dimension(:), allocatable :: current

  interface

  module function source_damage_anisoBrittle_init(source_length) result(mySources)
    integer, intent(in) :: source_length
    logical, dimension(:,:), allocatable :: mySources
  end function source_damage_anisoBrittle_init

  module function source_damage_anisoDuctile_init(source_length) result(mySources)
    integer, intent(in) :: source_length
    logical, dimension(:,:), allocatable :: mySources
  end function source_damage_anisoDuctile_init

  module function source_damage_isoBrittle_init(source_length) result(mySources)
    integer, intent(in) :: source_length
    logical, dimension(:,:), allocatable :: mySources
  end function source_damage_isoBrittle_init

  module function source_damage_isoDuctile_init(source_length) result(mySources)
    integer, intent(in) :: source_length
    logical, dimension(:,:), allocatable :: mySources
  end function source_damage_isoDuctile_init

  module function kinematics_cleavage_opening_init(kinematics_length) result(myKinematics)
    integer, intent(in) :: kinematics_length
    logical, dimension(:,:), allocatable :: myKinematics
  end function kinematics_cleavage_opening_init

  module function kinematics_slipplane_opening_init(kinematics_length) result(myKinematics)
    integer, intent(in) :: kinematics_length
    logical, dimension(:,:), allocatable :: myKinematics
  end function kinematics_slipplane_opening_init

    module subroutine source_damage_isoBrittle_deltaState(C, Fe, ph, me)
      integer, intent(in) :: ph,me
      real(pReal),  intent(in), dimension(3,3) :: &
        Fe
      real(pReal),  intent(in), dimension(6,6) :: &
        C
    end subroutine source_damage_isoBrittle_deltaState


    module subroutine source_damage_anisoBrittle_dotState(S, co, ip, el)
      integer, intent(in) :: &
        co, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
      real(pReal),  intent(in), dimension(3,3) :: &
        S
    end subroutine source_damage_anisoBrittle_dotState

    module subroutine source_damage_anisoDuctile_dotState(co, ip, el)
      integer, intent(in) :: &
        co, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
    end subroutine source_damage_anisoDuctile_dotState

    module subroutine source_damage_isoDuctile_dotState(co, ip, el)
      integer, intent(in) :: &
        co, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
    end subroutine source_damage_isoDuctile_dotState


  module subroutine source_damage_anisobrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)
    integer, intent(in) :: &
      phase, &                                                                                      !< phase ID of element
      constituent                                                                                   !< position of element within its phase instance
    real(pReal),  intent(in) :: &
      phi                                                                                           !< damage parameter
    real(pReal),  intent(out) :: &
      localphiDot, &
      dLocalphiDot_dPhi
  end subroutine source_damage_anisoBrittle_getRateAndItsTangent

  module subroutine source_damage_anisoDuctile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)
    integer, intent(in) :: &
      phase, &                                                                                      !< phase ID of element
      constituent                                                                                   !< position of element within its phase instance
    real(pReal),  intent(in) :: &
      phi                                                                                           !< damage parameter
    real(pReal),  intent(out) :: &
      localphiDot, &
      dLocalphiDot_dPhi
  end subroutine source_damage_anisoDuctile_getRateAndItsTangent

  module subroutine source_damage_isoBrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)
    integer, intent(in) :: &
      phase, &                                                                                      !< phase ID of element
      constituent                                                                                   !< position of element within its phase instance
    real(pReal),  intent(in) :: &
      phi                                                                                           !< damage parameter
    real(pReal),  intent(out) :: &
      localphiDot, &
      dLocalphiDot_dPhi
  end subroutine source_damage_isoBrittle_getRateAndItsTangent

  module subroutine source_damage_isoDuctile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)
    integer, intent(in) :: &
      phase, &                                                                                      !< phase ID of element
      constituent                                                                                   !< position of element within its phase instance
    real(pReal),  intent(in) :: &
      phi                                                                                           !< damage parameter
    real(pReal),  intent(out) :: &
      localphiDot, &
      dLocalphiDot_dPhi
  end subroutine source_damage_isoDuctile_getRateAndItsTangent

  module subroutine source_damage_anisoBrittle_results(phase,group)
    integer,          intent(in) :: phase
    character(len=*), intent(in) :: group
  end subroutine source_damage_anisoBrittle_results

  module subroutine source_damage_anisoDuctile_results(phase,group)
    integer,          intent(in) :: phase
    character(len=*), intent(in) :: group
  end subroutine source_damage_anisoDuctile_results

  module subroutine source_damage_isoBrittle_results(phase,group)
    integer,          intent(in) :: phase
    character(len=*), intent(in) :: group
  end subroutine source_damage_isoBrittle_results

  module subroutine source_damage_isoDuctile_results(phase,group)
    integer,          intent(in) :: phase
    character(len=*), intent(in) :: group
  end subroutine source_damage_isoDuctile_results

 end interface

contains

!----------------------------------------------------------------------------------------------
!< @brief initialize damage sources and kinematics mechanism
!----------------------------------------------------------------------------------------------
module subroutine damage_init

  integer :: &
    ph, &                                                                                           !< counter in phase loop
    Nconstituents
  class(tNode), pointer :: &
   phases, &
   phase, &
   sources, &
   kinematics

  phases => config_material%get('phase')

  allocate(current(phases%length))

  allocate(damageState (phases%length))
  allocate(phase_Nsources(phases%length),source = 0)           ! same for kinematics

  do ph = 1,phases%length

    Nconstituents = count(material_phaseAt == ph) * discretization_nIPs

    allocate(current(ph)%phi(Nconstituents),source=1.0_pReal)
    allocate(current(ph)%d_phi_d_dot_phi(Nconstituents),source=0.0_pReal)

    phase => phases%get(ph)
    sources => phase%get('source',defaultVal=emptyList)
    phase_Nsources(ph) = sources%length
    allocate(damageState(ph)%p(phase_Nsources(ph)))
  enddo

  allocate(phase_source(maxval(phase_Nsources),phases%length), source = DAMAGE_UNDEFINED_ID)

! initialize source mechanisms
  if(maxval(phase_Nsources) /= 0) then
    where(source_damage_isoBrittle_init   (maxval(phase_Nsources))) phase_source = DAMAGE_ISOBRITTLE_ID
    where(source_damage_isoDuctile_init   (maxval(phase_Nsources))) phase_source = DAMAGE_ISODUCTILE_ID
    where(source_damage_anisoBrittle_init (maxval(phase_Nsources))) phase_source = DAMAGE_ANISOBRITTLE_ID
    where(source_damage_anisoDuctile_init (maxval(phase_Nsources))) phase_source = DAMAGE_ANISODUCTILE_ID
  endif

!--------------------------------------------------------------------------------------------------
! initialize kinematic mechanisms
  allocate(phase_Nkinematics(phases%length),source = 0)
  do ph = 1,phases%length
    phase => phases%get(ph)
    kinematics => phase%get('kinematics',defaultVal=emptyList)
    phase_Nkinematics(ph) = kinematics%length
  enddo

  allocate(phase_kinematics(maxval(phase_Nkinematics),phases%length), source = KINEMATICS_undefined_ID)

  if(maxval(phase_Nkinematics) /= 0) then
    where(kinematics_cleavage_opening_init(maxval(phase_Nkinematics)))  phase_kinematics = KINEMATICS_cleavage_opening_ID
    where(kinematics_slipplane_opening_init(maxval(phase_Nkinematics))) phase_kinematics = KINEMATICS_slipplane_opening_ID
  endif

end subroutine damage_init


!----------------------------------------------------------------------------------------------
!< @brief returns local part of nonlocal damage driving force
!----------------------------------------------------------------------------------------------
module subroutine constitutive_damage_getRateAndItsTangents(phiDot, dPhiDot_dPhi, phi, ip, el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal), intent(in) :: &
    phi                                                                                             !< damage parameter
  real(pReal), intent(inout) :: &
    phiDot, &
    dPhiDot_dPhi

  real(pReal) :: &
    localphiDot, &
    dLocalphiDot_dPhi
  integer :: &
    phase, &
    grain, &
    source, &
    constituent

   phiDot = 0.0_pReal
   dPhiDot_dPhi = 0.0_pReal

   do grain = 1, homogenization_Nconstituents(material_homogenizationAt(el))
     phase = material_phaseAt(grain,el)
     constituent = material_phasememberAt(grain,ip,el)
     do source = 1, phase_Nsources(phase)
       select case(phase_source(source,phase))
         case (DAMAGE_ISOBRITTLE_ID)
           call source_damage_isobrittle_getRateAndItsTangent  (localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

         case (DAMAGE_ISODUCTILE_ID)
           call source_damage_isoductile_getRateAndItsTangent  (localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

         case (DAMAGE_ANISOBRITTLE_ID)
           call source_damage_anisobrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

         case (DAMAGE_ANISODUCTILE_ID)
           call source_damage_anisoductile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

         case default
         localphiDot = 0.0_pReal
         dLocalphiDot_dPhi = 0.0_pReal

      end select
      phiDot = phiDot + localphiDot
      dPhiDot_dPhi = dPhiDot_dPhi + dLocalphiDot_dPhi
    enddo
  enddo

end subroutine constitutive_damage_getRateAndItsTangents



!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with adaptive 1st order explicit Euler method
!> using Fixed Point Iteration to adapt the stepsize
!--------------------------------------------------------------------------------------------------
module function integrateDamageState(dt,co,ip,el) result(broken)

  real(pReal), intent(in) :: dt
  integer, intent(in) :: &
    el, &                                                                                            !< element index in element loop
    ip, &                                                                                            !< integration point index in ip loop
    co                                                                                               !< grain index in grain loop
  logical :: broken

  integer :: &
    NiterationState, &                                                                              !< number of iterations in state loop
    ph, &
    me, &
    so
  integer, dimension(maxval(phase_Nsources)) :: &
    size_so
  real(pReal) :: &
    zeta
  real(pReal), dimension(constitutive_source_maxSizeDotState) :: &
    r                                                                                               ! state residuum
  real(pReal), dimension(constitutive_source_maxSizeDotState,2,maxval(phase_Nsources)) :: source_dotState
  logical :: &
    converged_


  ph = material_phaseAt(co,el)
  me = material_phaseMemberAt(co,ip,el)

  converged_ = .true.
  broken = constitutive_damage_collectDotState(co,ip,el,ph,me)
  if(broken) return

  do so = 1, phase_Nsources(ph)
    size_so(so) = damageState(ph)%p(so)%sizeDotState
    damageState(ph)%p(so)%state(1:size_so(so),me) = damageState(ph)%p(so)%subState0(1:size_so(so),me) &
                                                  + damageState(ph)%p(so)%dotState (1:size_so(so),me) * dt
    source_dotState(1:size_so(so),2,so) = 0.0_pReal
  enddo

  iteration: do NiterationState = 1, num%nState

    do so = 1, phase_Nsources(ph)
      if(nIterationState > 1) source_dotState(1:size_so(so),2,so) = source_dotState(1:size_so(so),1,so)
      source_dotState(1:size_so(so),1,so) = damageState(ph)%p(so)%dotState(:,me)
    enddo

    broken = constitutive_damage_collectDotState(co,ip,el,ph,me)
    if(broken) exit iteration

    do so = 1, phase_Nsources(ph)
      zeta = damper(damageState(ph)%p(so)%dotState(:,me), &
                    source_dotState(1:size_so(so),1,so),&
                    source_dotState(1:size_so(so),2,so))
      damageState(ph)%p(so)%dotState(:,me) = damageState(ph)%p(so)%dotState(:,me) * zeta &
                                        + source_dotState(1:size_so(so),1,so)* (1.0_pReal - zeta)
      r(1:size_so(so)) = damageState(ph)%p(so)%state    (1:size_so(so),me)  &
                       - damageState(ph)%p(so)%subState0(1:size_so(so),me)  &
                       - damageState(ph)%p(so)%dotState (1:size_so(so),me) * dt
      damageState(ph)%p(so)%state(1:size_so(so),me) = damageState(ph)%p(so)%state(1:size_so(so),me) &
                                                - r(1:size_so(so))
      converged_ = converged_  .and. converged(r(1:size_so(so)), &
                                               damageState(ph)%p(so)%state(1:size_so(so),me), &
                                               damageState(ph)%p(so)%atol(1:size_so(so)))
    enddo

    if(converged_) then
      broken = constitutive_damage_deltaState(mech_F_e(ph,me),ph,me)
      exit iteration
    endif

  enddo iteration

  broken = broken .or. .not. converged_


  contains

  !--------------------------------------------------------------------------------------------------
  !> @brief calculate the damping for correction of state and dot state
  !--------------------------------------------------------------------------------------------------
  real(pReal) pure function damper(current,previous,previous2)

  real(pReal), dimension(:), intent(in) ::&
    current, previous, previous2

  real(pReal) :: dot_prod12, dot_prod22

  dot_prod12 = dot_product(current  - previous,  previous - previous2)
  dot_prod22 = dot_product(previous - previous2, previous - previous2)
  if ((dot_product(current,previous) < 0.0_pReal .or. dot_prod12 < 0.0_pReal) .and. dot_prod22 > 0.0_pReal) then
    damper = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)
  else
    damper = 1.0_pReal
  endif

  end function damper

end function integrateDamageState


!----------------------------------------------------------------------------------------------
!< @brief writes damage sources results to HDF5 output file
!----------------------------------------------------------------------------------------------
module subroutine damage_results(group,ph)

  character(len=*), intent(in) :: group
  integer,          intent(in) :: ph

  integer :: so

  sourceLoop: do so = 1, phase_Nsources(ph)

  if (phase_source(so,ph) /= DAMAGE_UNDEFINED_ID) &
    call results_closeGroup(results_addGroup(group//'sources/')) ! should be 'damage'

    sourceType: select case (phase_source(so,ph))

      case (DAMAGE_ISOBRITTLE_ID) sourceType
        call source_damage_isoBrittle_results(ph,group//'sources/')

      case (DAMAGE_ISODUCTILE_ID) sourceType
        call source_damage_isoDuctile_results(ph,group//'sources/')

      case (DAMAGE_ANISOBRITTLE_ID) sourceType
        call source_damage_anisoBrittle_results(ph,group//'sources/')

      case (DAMAGE_ANISODUCTILE_ID) sourceType
        call source_damage_anisoDuctile_results(ph,group//'sources/')

    end select sourceType

  enddo SourceLoop

end subroutine damage_results


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
function constitutive_damage_collectDotState(co,ip,el,ph,me) result(broken)

  integer, intent(in) :: &
    co, &                                                                                          !< component-ID me integration point
    ip, &                                                                                           !< integration point
    el, &                                                                                           !< element
    ph, &
    me
  integer :: &
    so                                                                                               !< counter in source loop
  logical :: broken


  broken = .false.

  SourceLoop: do so = 1, phase_Nsources(ph)

    sourceType: select case (phase_source(so,ph))

      case (DAMAGE_ISODUCTILE_ID) sourceType
        call source_damage_isoDuctile_dotState(co, ip, el)

      case (DAMAGE_ANISODUCTILE_ID) sourceType
        call source_damage_anisoDuctile_dotState(co, ip, el)

      case (DAMAGE_ANISOBRITTLE_ID) sourceType
        call source_damage_anisoBrittle_dotState(mech_S(ph,me),co, ip, el) ! correct stress?

    end select sourceType

    broken = broken .or. any(IEEE_is_NaN(damageState(ph)%p(so)%dotState(:,me)))

  enddo SourceLoop

end function constitutive_damage_collectDotState



!--------------------------------------------------------------------------------------------------
!> @brief for constitutive models having an instantaneous change of state
!> will return false if delta state is not needed/supported by the constitutive model
!--------------------------------------------------------------------------------------------------
function constitutive_damage_deltaState(Fe, ph, me) result(broken)

  integer, intent(in) :: &
    ph, &
    me
  real(pReal),   intent(in), dimension(3,3) :: &
    Fe                                                                                              !< elastic deformation gradient
  integer :: &
    so, &
    myOffset, &
    mySize
  logical :: &
    broken


  broken = .false.

  sourceLoop: do so = 1, phase_Nsources(ph)

     sourceType: select case (phase_source(so,ph))

      case (DAMAGE_ISOBRITTLE_ID) sourceType
        call source_damage_isoBrittle_deltaState(constitutive_homogenizedC(ph,me), Fe, ph,me)
        broken = any(IEEE_is_NaN(damageState(ph)%p(so)%deltaState(:,me)))
        if(.not. broken) then
          myOffset = damageState(ph)%p(so)%offsetDeltaState
          mySize   = damageState(ph)%p(so)%sizeDeltaState
          damageState(ph)%p(so)%state(myOffset + 1: myOffset + mySize,me) = &
          damageState(ph)%p(so)%state(myOffset + 1: myOffset + mySize,me) + damageState(ph)%p(so)%deltaState(1:mySize,me)
        endif

    end select sourceType

  enddo SourceLoop

end function constitutive_damage_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief checks if a source mechanism is active or not
!--------------------------------------------------------------------------------------------------
function source_active(source_label,src_length)  result(active_source)

  character(len=*), intent(in)         :: source_label                                              !< name of source mechanism
  integer,          intent(in)         :: src_length                                                !< max. number of sources in system
  logical, dimension(:,:), allocatable :: active_source

  class(tNode), pointer :: &
    phases, &
    phase, &
    sources, &
    src
  integer :: p,s

  phases => config_material%get('phase')
  allocate(active_source(src_length,phases%length), source = .false. )
  do p = 1, phases%length
    phase => phases%get(p)
    sources => phase%get('source',defaultVal=emptyList)
    do s = 1, sources%length
      src => sources%get(s)
      if(src%get_asString('type') == source_label) active_source(s,p) = .true.
    enddo
  enddo


end function source_active


!----------------------------------------------------------------------------------------------
!< @brief Set damage parameter
!----------------------------------------------------------------------------------------------
module subroutine constitutive_damage_set_phi(phi,co,ce)

  real(pReal), intent(in) :: phi
  integer, intent(in) :: ce, co


  current(material_phaseAt2(co,ce))%phi(material_phaseMemberAt2(co,ce)) = phi

end subroutine constitutive_damage_set_phi


module function constitutive_damage_get_phi(co,ip,el) result(phi)
  
  integer, intent(in) :: co, ip, el
  real(pReal) :: phi

  phi = current(material_phaseAt(co,el))%phi(material_phaseMemberAt(co,ip,el))

end function constitutive_damage_get_phi


end submodule constitutive_damage
