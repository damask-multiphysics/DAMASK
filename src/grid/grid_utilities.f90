!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Utilities used by the grid solver variants
!--------------------------------------------------------------------------------------------------
module grid_utilities
#include <petsc/finclude/petscsys.h>
  use PETScSys
#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
  use MPI_f08
#endif

  use prec
  use parallelization
  use math
  use rotations
  use IO
  use discretization_grid
  use discretization
  use spectral_utilities
  use homogenization
  use constants

#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
  implicit none(type,external)
#else
  implicit none
#endif
  private

!--------------------------------------------------------------------------------------------------
! derived types
  type, public :: tBoundaryCondition                                                                !< set of parameters defining a boundary condition
    real(pREAL), dimension(3,3)   :: values = 0.0_pREAL
    logical,     dimension(3,3)   :: mask   = .true.
    character(len=:), allocatable :: myType
  end type tBoundaryCondition

  type, public :: tSolutionParams
    real(pREAL), dimension(3,3) :: stress_BC
    logical, dimension(3,3)     :: stress_mask
    type(tRotation)             :: rotation_BC
    real(pREAL) :: Delta_t
  end type tSolutionParams

  public :: &
    utilities_maskedCompliance, &
    utilities_constitutiveResponse, &
    utilities_calculateRate, &
    utilities_forwardTensorField

contains


!--------------------------------------------------------------------------------------------------
!> @brief Calculate masked compliance tensor used to adjust F to fullfill stress BC.
!--------------------------------------------------------------------------------------------------
function utilities_maskedCompliance(rot_BC,mask_stress,C)

  real(pREAL),                dimension(3,3,3,3) :: utilities_maskedCompliance                      !< masked compliance
  real(pREAL),    intent(in), dimension(3,3,3,3) :: C                                               !< current average stiffness
  type(tRotation), intent(in)                    :: rot_BC                                          !< rotation of load frame
  logical,        intent(in), dimension(3,3)     :: mask_stress                                     !< mask of stress BC

  integer :: i, j
  logical, dimension(9)   :: mask_stressVector
  logical, dimension(9,9) :: mask
  real(pREAL), dimension(9,9) :: temp99_real
  integer :: size_reduced = 0
  real(pREAL),              dimension(:,:), allocatable ::  &
    s_reduced, &                                                                                    !< reduced compliance matrix (depending on number of stress BC)
    c_reduced, &                                                                                    !< reduced stiffness (depending on number of stress BC)
    sTimesC                                                                                         !< temp variable to check inversion
  logical :: errmatinv
  character(len=pSTRLEN):: formatString


  mask_stressVector = .not. reshape(transpose(mask_stress), [9])
  size_reduced = count(mask_stressVector)
  if (size_reduced > 0) then
    temp99_real = math_3333to99(rot_BC%rotate(C))

    do i = 1,9; do j = 1,9
      mask(i,j) = mask_stressVector(i) .and. mask_stressVector(j)
    end do; end do
    c_reduced = reshape(pack(temp99_Real,mask),[size_reduced,size_reduced])

    allocate(s_reduced,mold = c_reduced)
    call math_invert(s_reduced, errmatinv, c_reduced)                                               ! invert reduced stiffness
    if (any(IEEE_is_NaN(s_reduced))) errmatinv = .true.

!--------------------------------------------------------------------------------------------------
! check if inversion was successful
    sTimesC = matmul(c_reduced,s_reduced)
    errmatinv = errmatinv .or. any(dNeq(sTimesC,math_eye(size_reduced),1.0e-12_pREAL))
    if (errmatinv) then
      write(formatString, '(i2)') size_reduced
      formatString = '(/,1x,a,/,'//trim(formatString)//'('//trim(formatString)//'(2x,es9.2,1x)/))'
      print trim(formatString), 'C * S (load) ', transpose(matmul(c_reduced,s_reduced))
      print trim(formatString), 'C (load) ', transpose(c_reduced)
      print trim(formatString), 'S (load) ', transpose(s_reduced)
      if (errmatinv) error stop 'matrix inversion error'
    end if
    temp99_real = reshape(unpack(reshape(s_reduced,[size_reduced**2]),reshape(mask,[81]),0.0_pREAL),[9,9])
  else
    temp99_real = 0.0_pREAL
  end if

  utilities_maskedCompliance = math_99to3333(temp99_Real)

end function utilities_maskedCompliance


!--------------------------------------------------------------------------------------------------
!> @brief Calculate constitutive response.
!--------------------------------------------------------------------------------------------------
subroutine utilities_constitutiveResponse(status, P,P_av,C_volAvg,C_minmaxAvg,&
                                          F,Delta_t,rotation_BC)

  integer(kind(STATUS_OK)),  intent(out)                            :: status
  real(pREAL),    intent(out), dimension(3,3,3,3)                   :: C_volAvg, C_minmaxAvg        !< average stiffness
  real(pREAL),    intent(out), dimension(3,3)                       :: P_av                         !< average PK stress
  real(pREAL),    intent(out), dimension(3,3,cells(1),cells(2),cells3) :: P                         !< PK stress
  real(pREAL),    intent(in),  dimension(3,3,cells(1),cells(2),cells3) :: F                         !< deformation gradient target
  real(pREAL),    intent(in)                                        :: Delta_t                      !< loading time
  type(tRotation), intent(in),  optional                            :: rotation_BC                  !< rotation of load frame

  integer :: ce
  integer(MPI_INTEGER_KIND) :: err_MPI
  real(pREAL), dimension(3,3,3,3) :: dPdF_max,      dPdF_min
  real(pREAL)                     :: dPdF_norm_max, dPdF_norm_min
  real(pREAL), dimension(2) :: valueAndRank                                                         !< pair of min/max norm of dPdF to synchronize min/max of dPdF


  print'(/,1x,a)', '... evaluating constitutive response ......................................'
  flush(IO_STDOUT)

  homogenization_F  = reshape(F,[3,3,product(cells(1:2))*cells3])                                   ! set materialpoint target F to estimated field

  call homogenization_mechanical_response(status,Delta_t,1,product(cells(1:2))*cells3)              ! calculate P field

  P = reshape(homogenization_P, [3,3,cells(1),cells(2),cells3])
  P_av = sum(sum(sum(P,dim=5),dim=4),dim=3) * wgt
  call MPI_Allreduce(MPI_IN_PLACE,P_av,9_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  if (present(rotation_BC)) then
    if (any(dNeq(rotation_BC%asQuaternion(), real([1.0, 0.0, 0.0, 0.0],pREAL)))) &
      print'(/,1x,a,/,2(3(2x,f12.4,1x)/),3(2x,f12.4,1x))', &
      'Piola--Kirchhoff stress (lab) / MPa =', transpose(P_av)*1.e-6_pREAL
    P_av = rotation_BC%rotate(P_av)
  end if
  print'(/,1x,a,/,2(3(2x,f12.4,1x)/),3(2x,f12.4,1x))', &
    'Piola--Kirchhoff stress       / MPa =', transpose(P_av)*1.e-6_pREAL
  flush(IO_STDOUT)

  dPdF_max = 0.0_pREAL
  dPdF_norm_max = 0.0_pREAL
  dPdF_min = huge(1.0_pREAL)
  dPdF_norm_min = huge(1.0_pREAL)
  do ce = 1, product(cells(1:2))*cells3
    if (dPdF_norm_max < sum(homogenization_dPdF(1:3,1:3,1:3,1:3,ce)**2)) then
      dPdF_max = homogenization_dPdF(1:3,1:3,1:3,1:3,ce)
      dPdF_norm_max = sum(homogenization_dPdF(1:3,1:3,1:3,1:3,ce)**2)
    end if
    if (dPdF_norm_min > sum(homogenization_dPdF(1:3,1:3,1:3,1:3,ce)**2)) then
      dPdF_min = homogenization_dPdF(1:3,1:3,1:3,1:3,ce)
      dPdF_norm_min = sum(homogenization_dPdF(1:3,1:3,1:3,1:3,ce)**2)
    end if
  end do

  valueAndRank = [dPdF_norm_max,real(worldrank,pREAL)]
  call MPI_Allreduce(MPI_IN_PLACE,valueAndRank,1_MPI_INTEGER_KIND,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Bcast(dPdF_max,81_MPI_INTEGER_KIND,MPI_DOUBLE,int(valueAndRank(2),MPI_INTEGER_KIND),MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

  valueAndRank = [dPdF_norm_min,real(worldrank,pREAL)]
  call MPI_Allreduce(MPI_IN_PLACE,valueAndRank,1_MPI_INTEGER_KIND,MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Bcast(dPdF_min,81_MPI_INTEGER_KIND,MPI_DOUBLE,int(valueAndRank(2),MPI_INTEGER_KIND),MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

  C_minmaxAvg = 0.5_pREAL*(dPdF_max + dPdF_min)

  C_volAvg = sum(homogenization_dPdF,dim=5)
  call MPI_Allreduce(MPI_IN_PLACE,C_volAvg,81_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  C_volAvg = C_volAvg * wgt


end subroutine utilities_constitutiveResponse


!--------------------------------------------------------------------------------------------------
!> @brief Calculate forward rate, either as local guess or as homogeneous add on.
!--------------------------------------------------------------------------------------------------
pure function utilities_calculateRate(heterogeneous,field0,field,dt,avRate)

  real(pREAL), intent(in), dimension(3,3) :: &
    avRate                                                                                          !< homogeneous addon
  real(pREAL), intent(in) :: &
    dt                                                                                              !< Delta_t between field0 and field
  logical, intent(in) :: &
    heterogeneous                                                                                   !< calculate field of rates
  real(pREAL), intent(in), dimension(3,3,cells(1),cells(2),cells3) :: &
    field0, &                                                                                       !< data of previous step
    field                                                                                           !< data of current step
  real(pREAL),             dimension(3,3,cells(1),cells(2),cells3) :: &
    utilities_calculateRate


  utilities_calculateRate = merge((field-field0) / dt, &
                                  spread(spread(spread(avRate,3,cells(1)),4,cells(2)),5,cells3), &
                                  heterogeneous)

end function utilities_calculateRate


!--------------------------------------------------------------------------------------------------
!> @brief forwards a field with a pointwise given rate, if aim is given,
!> ensures that the average matches the aim
!--------------------------------------------------------------------------------------------------
function utilities_forwardTensorField(Delta_t,field_lastInc,rate,aim)

  real(pREAL), intent(in) :: &
    Delta_t                                                                                         !< Delta_t of current step
  real(pREAL), intent(in),           dimension(3,3,cells(1),cells(2),cells3) :: &
    field_lastInc, &                                                                                !< initial field
    rate                                                                                            !< rate by which to forward
  real(pREAL), intent(in), optional, dimension(3,3) :: &
    aim                                                                                             !< average field value aim

  real(pREAL),                       dimension(3,3,cells(1),cells(2),cells3) :: &
    utilities_forwardTensorField
  real(pREAL),                       dimension(3,3) :: fieldDiff                                    !< <a + adot*t> - aim
  integer(MPI_INTEGER_KIND) :: err_MPI


  utilities_forwardTensorField = field_lastInc + rate*Delta_t
  if (present(aim)) then                                                                            !< correct to match average
    fieldDiff = sum(sum(sum(utilities_forwardTensorField,dim=5),dim=4),dim=3)*wgt
    call MPI_Allreduce(MPI_IN_PLACE,fieldDiff,9_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)
    fieldDiff = fieldDiff - aim
    utilities_forwardTensorField = utilities_forwardTensorField &
                                 - spread(spread(spread(fieldDiff,3,cells(1)),4,cells(2)),5,cells3)
  end if

end function utilities_forwardTensorField

end module grid_utilities
