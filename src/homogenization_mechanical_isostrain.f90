!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Isostrain (full constraint Taylor assuption) homogenization scheme
!--------------------------------------------------------------------------------------------------
submodule(homogenization:mechanical) isostrain

  enum, bind(c); enumerator :: &
    parallel_ID, &
    average_ID
  end enum

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    integer :: &
      N_constituents
    integer(kind(average_ID)) :: &
      mapping
  end type

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstances)


contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_isostrain_init

  integer :: &
    h, &
    Nmaterialpoints
  class(tNode), pointer :: &
    material_homogenization, &
    homog, &
    homogMech

  print'(/,a)', ' <<<+-  homogenization:mechanical:isostrain init  -+>>>'

  print'(a,i2)', ' # instances: ',count(homogenization_type == HOMOGENIZATION_ISOSTRAIN_ID); flush(IO_STDOUT)

  material_homogenization => config_material%get('homogenization')
  allocate(param(material_homogenization%length))                                                   ! one container of parameters per homog

  do h = 1, size(homogenization_type)
    if (homogenization_type(h) /= HOMOGENIZATION_ISOSTRAIN_ID) cycle
    homog => material_homogenization%get(h)
    homogMech => homog%get('mechanical')
    associate(prm => param(h))

    prm%N_constituents = homogenization_Nconstituents(h)
    select case(homogMech%get_asString('mapping',defaultVal = 'sum'))
      case ('sum')
        prm%mapping = parallel_ID
      case ('avg')
        prm%mapping = average_ID
      case default
        call IO_error(211,ext_msg='sum'//' (mechanical_isostrain)')
    end select

    Nmaterialpoints = count(material_homogenizationAt == h)
    homogState(h)%sizeState       = 0
    allocate(homogState(h)%state0   (0,Nmaterialpoints))
    allocate(homogState(h)%state    (0,Nmaterialpoints))

    end associate

  enddo

end subroutine mechanical_isostrain_init


!--------------------------------------------------------------------------------------------------
!> @brief partitions the deformation gradient onto the constituents
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_isostrain_partitionDeformation(F,avgF)

  real(pReal),   dimension (:,:,:), intent(out) :: F                                                !< partitioned deformation gradient

  real(pReal),   dimension (3,3),   intent(in)  :: avgF                                             !< average deformation gradient at material point

  F = spread(avgF,3,size(F,3))

end subroutine mechanical_isostrain_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_isostrain_averageStressAndItsTangent(avgP,dAvgPdAvgF,P,dPdF,ho)

  real(pReal),   dimension (3,3),       intent(out) :: avgP                                         !< average stress at material point
  real(pReal),   dimension (3,3,3,3),   intent(out) :: dAvgPdAvgF                                   !< average stiffness at material point

  real(pReal),   dimension (:,:,:),     intent(in)  :: P                                            !< partitioned stresses
  real(pReal),   dimension (:,:,:,:,:), intent(in)  :: dPdF                                         !< partitioned stiffnesses
  integer,                              intent(in)  :: ho

  associate(prm => param(ho))

  select case (prm%mapping)
    case (parallel_ID)
      avgP       = sum(P,3)
      dAvgPdAvgF = sum(dPdF,5)
    case (average_ID)
      avgP       = sum(P,3)   /real(prm%N_constituents,pReal)
      dAvgPdAvgF = sum(dPdF,5)/real(prm%N_constituents,pReal)
  end select

  end associate

end subroutine mechanical_isostrain_averageStressAndItsTangent

end submodule isostrain
