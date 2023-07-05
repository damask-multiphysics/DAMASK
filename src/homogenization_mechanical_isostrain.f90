!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Isostrain (full constraint Taylor assuption) homogenization scheme
!--------------------------------------------------------------------------------------------------
submodule(homogenization:mechanical) isostrain

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
module subroutine isostrain_init

  integer :: &
    ho, &
    Nmembers

  print'(/,1x,a)', '<<<+-  homogenization:mechanical:isostrain init  -+>>>'

  print'(/,a,i0)', ' # homogenizations: ',count(mechanical_type == MECHANICAL_ISOSTRAIN_ID)
  flush(IO_STDOUT)

  do ho = 1, size(mechanical_type)
    if (mechanical_type(ho) /= MECHANICAL_ISOSTRAIN_ID) cycle

    Nmembers = count(material_ID_homogenization == ho)
    homogState(ho)%sizeState = 0
    allocate(homogState(ho)%state0(0,Nmembers))
    allocate(homogState(ho)%state (0,Nmembers))

  end do

end subroutine isostrain_init


!--------------------------------------------------------------------------------------------------
!> @brief partitions the deformation gradient onto the constituents
!--------------------------------------------------------------------------------------------------
module subroutine isostrain_partitionDeformation(F,avgF)

  real(pREAL),   dimension (:,:,:), intent(out) :: F                                                !< partitioned deformation gradient

  real(pREAL),   dimension (3,3),   intent(in)  :: avgF                                             !< average deformation gradient at material point


  F = spread(avgF,3,size(F,3))

end subroutine isostrain_partitionDeformation

end submodule isostrain
