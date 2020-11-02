!----------------------------------------------------------------------------------------------------
!> @brief internal microstructure state for all damage sources and kinematics constitutive models  
!----------------------------------------------------------------------------------------------------
submodule(constitutive) constitutive_damage

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
    ph                                                                                              !< counter in phase loop
  class(tNode), pointer :: &
   phases, &
   phase, &
   sources, &
   kinematics

  phases => config_material%get('phase')

  allocate(sourceState (phases%length))
  allocate(phase_Nsources(phases%length),source = 0)           ! same for kinematics

  do ph = 1,phases%length
    phase => phases%get(ph)
    sources => phase%get('source',defaultVal=emptyList)
    phase_Nsources(ph) = sources%length
    allocate(sourceState(ph)%p(phase_Nsources(ph)))
  enddo

  allocate(phase_source(maxval(phase_Nsources),phases%length), source = SOURCE_undefined_ID) 

! initialize source mechanisms
  if(maxval(phase_Nsources) /= 0) then
    where(source_damage_isoBrittle_init   (maxval(phase_Nsources))) phase_source = SOURCE_damage_isoBrittle_ID
    where(source_damage_isoDuctile_init   (maxval(phase_Nsources))) phase_source = SOURCE_damage_isoDuctile_ID
    where(source_damage_anisoBrittle_init (maxval(phase_Nsources))) phase_source = SOURCE_damage_anisoBrittle_ID
    where(source_damage_anisoDuctile_init (maxval(phase_Nsources))) phase_source = SOURCE_damage_anisoDuctile_ID
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
         case (SOURCE_damage_isoBrittle_ID)
           call source_damage_isobrittle_getRateAndItsTangent  (localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

         case (SOURCE_damage_isoDuctile_ID)
           call source_damage_isoductile_getRateAndItsTangent  (localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

         case (SOURCE_damage_anisoBrittle_ID)
           call source_damage_anisobrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

         case (SOURCE_damage_anisoDuctile_ID)
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


!----------------------------------------------------------------------------------------------
!< @brief writes damage sources results to HDF5 output file
!----------------------------------------------------------------------------------------------
module subroutine damage_results

  integer :: p,i
  character(len=pStringLen) :: group

  do p = 1, size(material_name_phase)

    sourceLoop: do i = 1, phase_Nsources(p)
    group = trim('current/constituent')//'/'//trim(material_name_phase(p))
    group = trim(group)//'/sources'
    call results_closeGroup(results_addGroup(group))

      sourceType: select case (phase_source(i,p))

        case (SOURCE_damage_anisoBrittle_ID) sourceType
          call source_damage_anisoBrittle_results(p,group)
        case (SOURCE_damage_anisoDuctile_ID) sourceType
          call source_damage_anisoDuctile_results(p,group)
        case (SOURCE_damage_isoBrittle_ID) sourceType
          call source_damage_isoBrittle_results(p,group)
        case (SOURCE_damage_isoDuctile_ID) sourceType
          call source_damage_isoDuctile_results(p,group)
      end select sourceType

    enddo SourceLoop
  enddo

end subroutine damage_results


end submodule constitutive_damage
