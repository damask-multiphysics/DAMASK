!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Denny Tjahjanto, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Relaxed grain cluster (RGC) homogenization scheme
!> N_constituents is defined as p x q x r (cluster)
!--------------------------------------------------------------------------------------------------
submodule(homogenization:mechanical) RGC
  use rotations
  use crystal

  type :: tParameters
    integer, dimension(:), allocatable :: &
      N_constituents
    real(pREAL) :: &
      xi_alpha, &
      c_Alpha
    real(pREAL), dimension(:), allocatable :: &
      D_alpha, &
      a_g
    character(len=pSTRLEN), allocatable, dimension(:) :: &
      output
  end type tParameters

  type :: tRGCstate
    real(pREAL), pointer,     dimension(:,:) :: &
      relaxationVector
  end type tRGCstate

  type :: tRGCdependentState
    real(pREAL), allocatable,     dimension(:) :: &
      volumeDiscrepancy, &
      relaxationRate_avg, &
      relaxationRate_max
    real(pREAL), allocatable,     dimension(:,:) :: &
      mismatch
    real(pREAL), allocatable,     dimension(:,:,:) :: &
      orientation
  end type tRGCdependentState

  type :: tNumerics_RGC
    real(pREAL) :: &
      atol, &                                                                                       !< absolute tolerance of RGC residuum
      rtol, &                                                                                       !< relative tolerance of RGC residuum
      absMax, &                                                                                     !< absolute maximum of RGC residuum
      relMax, &                                                                                     !< relative maximum of RGC residuum
      pPert, &                                                                                      !< perturbation for computing RGC penalty tangent
      xSmoo, &                                                                                      !< RGC penalty smoothing parameter (hyperbolic tangent)
      viscPower, &                                                                                  !< power (sensitivity rate) of numerical viscosity in RGC scheme, Default 1.0e0: Newton viscosity (linear model)
      viscModus, &                                                                                  !< stress modulus of RGC numerical viscosity, Default 0.0e0: No viscosity is applied
      refRelaxRate, &                                                                               !< reference relaxation rate in RGC viscosity
      maxdRelax, &                                                                                  !< threshold of maximum relaxation vector increment (if exceed this then cutback)
      maxVolDiscr, &                                                                                !< threshold of maximum volume discrepancy allowed
      volDiscrMod, &                                                                                !< stiffness of RGC volume discrepancy (zero = without volume discrepancy constraint)
      volDiscrPow                                                                                   !< powerlaw penalty for volume discrepancy
  end type tNumerics_RGC

  type(tparameters),          dimension(:), allocatable :: &
    param
  type(tRGCstate),            dimension(:), allocatable :: &
    state, &
    state0
  type(tRGCdependentState),   dimension(:), allocatable :: &
    dependentState
  type(tNumerics_RGC) :: &
    num                                                                                             ! numerics parameters. Better name?

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all necessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
module subroutine RGC_init()

  integer :: &
    ho, &
    Nmembers, &
    sizeState, nIntFaceTot

  class(tDict), pointer :: &
    num_homogenization, &
    num_mechanical, &
    num_RGC, &                                                                                      ! pointer to RGC numerics data
    material_homogenization, &
    homog, &
    homogMech

  print'(/,1x,a)', '<<<+-  homogenization:mechanical:RGC init  -+>>>'

  print'(/,a,i0)', ' # homogenizations: ',count(mechanical_type == MECHANICAL_RGC_ID)
  flush(IO_STDOUT)

  print'(/,1x,a)', 'D.D. Tjahjanto et al., International Journal of Material Forming 2(1):939–942, 2009'
  print'(  1x,a)', 'https://doi.org/10.1007/s12289-009-0619-1'//IO_EOL

  print'(/,1x,a)', 'D.D. Tjahjanto et al., Modelling and Simulation in Materials Science and Engineering 18:015006, 2010'
  print'(  1x,a)', 'https://doi.org/10.1088/0965-0393/18/1/015006'//IO_EOL


  material_homogenization => config_material%get_dict('homogenization')
  allocate(param(material_homogenization%length))
  allocate(state(material_homogenization%length))
  allocate(state0(material_homogenization%length))
  allocate(dependentState(material_homogenization%length))

  num_homogenization => config_numerics%get_dict('homogenization',defaultVal=emptyDict)
  num_mechanical => num_homogenization%get_dict('mechanical',defaultVal=emptyDict)
  num_RGC => num_mechanical%get_dict('RGC',defaultVal=emptyDict)

  num%atol         =  num_RGC%get_asReal('eps_abs_P',          defaultVal=1.0e+4_pREAL)
  num%rtol         =  num_RGC%get_asReal('eps_rel_P',          defaultVal=1.0e-3_pREAL)
  num%absMax       =  num_RGC%get_asReal('eps_abs_max',        defaultVal=1.0e+10_pREAL)
  num%relMax       =  num_RGC%get_asReal('eps_rel_max',        defaultVal=1.0e+2_pREAL)
  num%pPert        =  num_RGC%get_asReal('Delta_a',            defaultVal=1.0e-7_pREAL)
  num%xSmoo        =  num_RGC%get_asReal('relevant_mismatch',  defaultVal=1.0e-5_pREAL)
  num%viscPower    =  num_RGC%get_asReal('viscosity_exponent', defaultVal=1.0e+0_pREAL)
  num%viscModus    =  num_RGC%get_asReal('viscosity_modulus',  defaultVal=0.0e+0_pREAL)
  num%refRelaxRate =  num_RGC%get_asReal('dot_a_ref',          defaultVal=1.0e-3_pREAL)
  num%maxdRelax    =  num_RGC%get_asReal('dot_a_max',          defaultVal=1.0e+0_pREAL)
  num%maxVolDiscr  =  num_RGC%get_asReal('Delta_V_max',        defaultVal=1.0e-5_pREAL)
  num%volDiscrMod  =  num_RGC%get_asReal('Delta_V_modulus',    defaultVal=1.0e+12_pREAL)
  num%volDiscrPow  =  num_RGC%get_asReal('Delta_V_exponent',   defaultVal=5.0_pREAL)

  if (num%atol <= 0.0_pREAL)         call IO_error(301,ext_msg='eps_abs_P')
  if (num%rtol <= 0.0_pREAL)         call IO_error(301,ext_msg='eps_rel_P')
  if (num%absMax <= 0.0_pREAL)       call IO_error(301,ext_msg='eps_abs_max')
  if (num%relMax <= 0.0_pREAL)       call IO_error(301,ext_msg='eps_rel_max')
  if (num%pPert <= 0.0_pREAL)        call IO_error(301,ext_msg='Delta_a')
  if (num%xSmoo <= 0.0_pREAL)        call IO_error(301,ext_msg='relevant_mismatch')
  if (num%viscPower < 0.0_pREAL)     call IO_error(301,ext_msg='viscosity_exponent')
  if (num%viscModus < 0.0_pREAL)     call IO_error(301,ext_msg='viscosity_modulus')
  if (num%refRelaxRate <= 0.0_pREAL) call IO_error(301,ext_msg='dot_a_ref')
  if (num%maxdRelax <= 0.0_pREAL)    call IO_error(301,ext_msg='dot_a_max')
  if (num%maxVolDiscr <= 0.0_pREAL)  call IO_error(301,ext_msg='Delta_V_max')
  if (num%volDiscrMod < 0.0_pREAL)   call IO_error(301,ext_msg='Delta_V_modulus')
  if (num%volDiscrPow <= 0.0_pREAL)  call IO_error(301,ext_msg='Delta_V_exponent')


  do ho = 1, size(mechanical_type)
    if (mechanical_type(ho) /= MECHANICAL_RGC_ID) cycle
    homog => material_homogenization%get_dict(ho)
    homogMech => homog%get_dict('mechanical')
    associate(prm => param(ho), &
              stt => state(ho), &
              st0 => state0(ho), &
              dst => dependentState(ho))

#if defined (__GFORTRAN__)
    prm%output = output_as1dStr(homogMech)
#else
    prm%output = homogMech%get_as1dStr('output',defaultVal=emptyStrArray)
#endif

    prm%N_constituents = homogMech%get_as1dInt('cluster_size',requiredSize=3)
    if (homogenization_Nconstituents(ho) /= product(prm%N_constituents)) &
      call IO_error(211,ext_msg='N_constituents (RGC)')

    prm%xi_alpha = homogMech%get_asReal('xi_alpha')
    prm%c_alpha  = homogMech%get_asReal('c_alpha')

    prm%D_alpha  = homogMech%get_as1dReal('D_alpha', requiredSize=3)
    prm%a_g      = homogMech%get_as1dReal('a_g',     requiredSize=3)

    Nmembers = count(material_ID_homogenization == ho)
    nIntFaceTot = 3*(  (prm%N_constituents(1)-1)*prm%N_constituents(2)*prm%N_constituents(3) &
                      + prm%N_constituents(1)*(prm%N_constituents(2)-1)*prm%N_constituents(3) &
                      + prm%N_constituents(1)*prm%N_constituents(2)*(prm%N_constituents(3)-1))
    sizeState = nIntFaceTot

    homogState(ho)%sizeState = sizeState
    allocate(homogState(ho)%state0   (sizeState,Nmembers), source=0.0_pREAL)
    allocate(homogState(ho)%state    (sizeState,Nmembers), source=0.0_pREAL)

    stt%relaxationVector => homogState(ho)%state(1:nIntFaceTot,:)
    st0%relaxationVector => homogState(ho)%state0(1:nIntFaceTot,:)

    allocate(dst%volumeDiscrepancy(   Nmembers), source=0.0_pREAL)
    allocate(dst%relaxationRate_avg(  Nmembers), source=0.0_pREAL)
    allocate(dst%relaxationRate_max(  Nmembers), source=0.0_pREAL)
    allocate(dst%mismatch(          3,Nmembers), source=0.0_pREAL)

!--------------------------------------------------------------------------------------------------
! assigning cluster orientations
    dependentState(ho)%orientation = spread(eu2om(prm%a_g*inRad),3,Nmembers)
    !dst%orientation = spread(eu2om(prm%a_g*inRad),3,Nmembers) ifort version 18.0.1 crashes (for whatever reason)

    end associate

  end do

end subroutine RGC_init


!--------------------------------------------------------------------------------------------------
!> @brief partitions the deformation gradient onto the constituents
!--------------------------------------------------------------------------------------------------
module subroutine RGC_partitionDeformation(F,avgF,ce)

  real(pREAL),   dimension (:,:,:), intent(out) :: F                                                !< partitioned F  per grain

  real(pREAL),   dimension (3,3),   intent(in)  :: avgF                                             !< averaged F
  integer,                          intent(in)  :: &
    ce

  real(pREAL), dimension(3) :: aVect,nVect
  integer,     dimension(4) :: intFace
  integer,     dimension(3) :: iGrain3
  integer ::  iGrain,iFace,i,j,ho,en

  associate(prm => param(material_ID_homogenization(ce)))

  ho = material_ID_homogenization(ce)
  en = material_entry_homogenization(ce)
!--------------------------------------------------------------------------------------------------
! compute the deformation gradient of individual grains due to relaxations
  F = 0.0_pREAL
  do iGrain = 1,product(prm%N_constituents)
    iGrain3 = grain1to3(iGrain,prm%N_constituents)
    do iFace = 1,6
      intFace = getInterface(iFace,iGrain3)                                                         ! identifying 6 interfaces of each grain
      aVect = relaxationVector(intFace,ho,en)                                                       ! get the relaxation vectors for each interface from global relaxation vector array
      nVect = interfaceNormal(intFace,ho,en)
      forall (i=1:3,j=1:3) &
        F(i,j,iGrain) = F(i,j,iGrain) + aVect(i)*nVect(j)                                           ! calculating deformation relaxations due to interface relaxation
    end do
    F(1:3,1:3,iGrain) = F(1:3,1:3,iGrain) + avgF                                                    ! resulting relaxed deformation gradient
  end do

  end associate

end subroutine RGC_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief update the internal state of the homogenization scheme and tell whether "done" and
! "happy" with result
!--------------------------------------------------------------------------------------------------
module function RGC_updateState(P,F,avgF,dt,dPdF,ce) result(doneAndHappy)
      logical, dimension(2) :: doneAndHappy
      real(pREAL), dimension(:,:,:),     intent(in)    :: &
        P,&                                                                                         !< partitioned stresses
        F                                                                                           !< partitioned deformation gradients
      real(pREAL), dimension(:,:,:,:,:), intent(in) :: dPdF                                         !< partitioned stiffnesses
      real(pREAL), dimension(3,3),       intent(in) :: avgF                                         !< average F
      real(pREAL),                       intent(in) :: dt                                           !< time increment
      integer,                           intent(in) :: &
        ce                                                                                          !< cell

  integer, dimension(4) :: intFaceN,intFaceP,faceID
  integer, dimension(3) :: nGDim,iGr3N,iGr3P
  integer :: ho,iNum,i,j,nIntFaceTot,iGrN,iGrP,iMun,iFace,k,l,ipert,nGrain, en
  real(pREAL), dimension(3,3,size(P,3)) :: R,pF,pR,D,pD
  real(pREAL), dimension(3,size(P,3))   :: NN,devNull
  real(pREAL), dimension(3)             :: normP,normN,mornP,mornN
  real(pREAL) :: residMax,stresMax
  logical :: error
  real(pREAL), dimension(:,:), allocatable :: tract,jmatrix,jnverse,smatrix,pmatrix,rmatrix
  real(pREAL), dimension(:), allocatable   :: resid,relax,p_relax,p_resid,drelax

  zeroTimeStep: if (dEq0(dt)) then
    doneAndHappy = .true.                                                                   ! pretend everything is fine and return
    return
  end if zeroTimeStep

  ho  = material_ID_homogenization(ce)
  en = material_entry_homogenization(ce)

  associate(stt => state(ho), st0 => state0(ho), dst => dependentState(ho), prm => param(ho))

!--------------------------------------------------------------------------------------------------
! get the dimension of the cluster (grains and interfaces)
  nGDim  = prm%N_constituents
  nGrain = product(nGDim)
  nIntFaceTot = (nGDim(1)-1)*nGDim(2)*nGDim(3) &
              + nGDim(1)*(nGDim(2)-1)*nGDim(3) &
              + nGDim(1)*nGDim(2)*(nGDim(3)-1)

!--------------------------------------------------------------------------------------------------
! allocate the size of the global relaxation arrays/jacobian matrices depending on the size of the cluster
  allocate(resid(3*nIntFaceTot), source=0.0_pREAL)
  allocate(tract(nIntFaceTot,3), source=0.0_pREAL)
  relax  = stt%relaxationVector(:,en)
  drelax = stt%relaxationVector(:,en) - st0%relaxationVector(:,en)

!--------------------------------------------------------------------------------------------------
! computing interface mismatch and stress penalty tensor for all interfaces of all grains
  call stressPenalty(R,NN,avgF,F,ho,en)

!--------------------------------------------------------------------------------------------------
! calculating volume discrepancy and stress penalty related to overall volume discrepancy
  call volumePenalty(D,dst%volumeDiscrepancy(en),avgF,F,nGrain)

!------------------------------------------------------------------------------------------------
! computing the residual stress from the balance of traction at all (interior) interfaces
  do iNum = 1,nIntFaceTot
    faceID = interface1to4(iNum,param(ho)%N_constituents)                                           ! identifying the interface ID in local coordinate system (4-dimensional index)

!--------------------------------------------------------------------------------------------------
! identify the left/bottom/back grain (-|N)
    iGr3N = faceID(2:4)                                                                             ! identifying the grain ID in local coordinate system (3-dimensional index)
    iGrN = grain3to1(iGr3N,param(ho)%N_constituents)                                                ! translate the local grain ID into global coordinate system (1-dimensional index)
    intFaceN = getInterface(2*faceID(1),iGr3N)
    normN = interfaceNormal(intFaceN,ho,en)

!--------------------------------------------------------------------------------------------------
! identify the right/up/front grain (+|P)
    iGr3P = iGr3N
    iGr3P(faceID(1)) = iGr3N(faceID(1))+1                                                           ! identifying the grain ID in local coordinate system (3-dimensional index)
    iGrP = grain3to1(iGr3P,param(ho)%N_constituents)                                                ! translate the local grain ID into global coordinate system (1-dimensional index)
    intFaceP = getInterface(2*faceID(1)-1,iGr3P)
    normP = interfaceNormal(intFaceP,ho,en)

!--------------------------------------------------------------------------------------------------
! compute the residual of traction at the interface (in local system, 4-dimensional index)
    do i = 1,3
      tract(iNum,i) = sign(num%viscModus*(abs(drelax(i+3*(iNum-1)))/(num%refRelaxRate*dt))**num%viscPower, &
                           drelax(i+3*(iNum-1)))                                                    ! contribution from the relaxation viscosity
      do j = 1,3
        tract(iNum,i) = tract(iNum,i) + (P(i,j,iGrP) + R(i,j,iGrP) + D(i,j,iGrP))*normP(j) &        ! contribution from material stress P, mismatch penalty R, and volume penalty D projected into the interface
                                      + (P(i,j,iGrN) + R(i,j,iGrN) + D(i,j,iGrN))*normN(j)
        resid(i+3*(iNum-1)) = tract(iNum,i)                                                         ! translate the local residual into global 1-dimensional residual array
      end do
    end do

  end do

!--------------------------------------------------------------------------------------------------
! convergence check for stress residual
  stresMax = maxval(abs(P))                                                                         ! get the maximum of first Piola-Kirchhoff (material) stress
  residMax = maxval(abs(tract))                                                                     ! get the maximum of the residual

  doneAndHappy = .false.

!--------------------------------------------------------------------------------------------------
!  If convergence reached => done and happy
  if (residMax < num%rtol*stresMax .or. residMax < num%atol) then
    doneAndHappy = .true.

    dst%mismatch(1:3,en)       = sum(NN,2)/real(nGrain,pREAL)
    dst%relaxationRate_avg(en) = sum(abs(drelax))/dt/real(3*nIntFaceTot,pREAL)
    dst%relaxationRate_max(en) = maxval(abs(drelax))/dt

    return

!--------------------------------------------------------------------------------------------------
! if residual blows-up => done but unhappy
  elseif (residMax > num%relMax*stresMax .or. residMax > num%absMax) then                           ! try to restart when residual blows up exceeding maximum bound
    doneAndHappy = [.true.,.false.]                                                         ! with direct cut-back
   return
  end if

!---------------------------------------------------------------------------------------------------
! construct the global Jacobian matrix for updating the global relaxation vector array when
! convergence is not yet reached ...

!--------------------------------------------------------------------------------------------------
! ... of the constitutive stress tangent, assembled from dPdF or material constitutive model "smatrix"
  allocate(smatrix(3*nIntFaceTot,3*nIntFaceTot), source=0.0_pREAL)
  do iNum = 1,nIntFaceTot
    faceID = interface1to4(iNum,param(ho)%N_constituents)                                           ! assembling of local dPdF into global Jacobian matrix

!--------------------------------------------------------------------------------------------------
! identify the left/bottom/back grain (-|N)
    iGr3N = faceID(2:4)                                                                             ! identifying the grain ID in local coordinate sytem
    iGrN = grain3to1(iGr3N,param(ho)%N_constituents)                                                ! translate into global grain ID
    intFaceN = getInterface(2*faceID(1),iGr3N)                                                      ! identifying the connecting interface in local coordinate system
    normN = interfaceNormal(intFaceN,ho,en)
    do iFace = 1,6
      intFaceN = getInterface(iFace,iGr3N)                                                          ! identifying all interfaces that influence relaxation of the above interface
      mornN = interfaceNormal(intFaceN,ho,en)
      iMun = interface4to1(intFaceN,param(ho)%N_constituents)                                       ! translate the interfaces ID into local 4-dimensional index
      if (iMun > 0) then                                                                            ! get the corresponding tangent
        do i=1,3; do j=1,3; do k=1,3; do l=1,3
          smatrix(3*(iNum-1)+i,3*(iMun-1)+j) = smatrix(3*(iNum-1)+i,3*(iMun-1)+j) &
                                             + dPdF(i,k,j,l,iGrN)*normN(k)*mornN(l)
        end do;end do;end do;end do
! projecting the material tangent dPdF into the interface
! to obtain the Jacobian matrix contribution of dPdF
      end if
    end do

!--------------------------------------------------------------------------------------------------
! identify the right/up/front grain (+|P)
    iGr3P = iGr3N
    iGr3P(faceID(1)) = iGr3N(faceID(1))+1                                                           ! identifying the grain ID in local coordinate sytem
    iGrP = grain3to1(iGr3P,param(ho)%N_constituents)                                                ! translate into global grain ID
    intFaceP = getInterface(2*faceID(1)-1,iGr3P)                                                    ! identifying the connecting interface in local coordinate system
    normP = interfaceNormal(intFaceP,ho,en)
    do iFace = 1,6
      intFaceP = getInterface(iFace,iGr3P)                                                          ! identifying all interfaces that influence relaxation of the above interface
      mornP = interfaceNormal(intFaceP,ho,en)
      iMun = interface4to1(intFaceP,param(ho)%N_constituents)                                       ! translate the interfaces ID into local 4-dimensional index
      if (iMun > 0) then                                                                            ! get the corresponding tangent
        do i=1,3; do j=1,3; do k=1,3; do l=1,3
          smatrix(3*(iNum-1)+i,3*(iMun-1)+j) = smatrix(3*(iNum-1)+i,3*(iMun-1)+j) &
                                             + dPdF(i,k,j,l,iGrP)*normP(k)*mornP(l)
        end do;end do;end do;end do
      end if
    end do
  end do

!--------------------------------------------------------------------------------------------------
! ... of the stress penalty tangent (mismatch penalty and volume penalty, computed using numerical
! perturbation method) "pmatrix"
  allocate(pmatrix(3*nIntFaceTot,3*nIntFaceTot), source=0.0_pREAL)
  allocate(p_relax(3*nIntFaceTot),               source=0.0_pREAL)
  allocate(p_resid(3*nIntFaceTot),               source=0.0_pREAL)

  do ipert = 1,3*nIntFaceTot
    p_relax = relax
    p_relax(ipert) = relax(ipert) + num%pPert                                                       ! perturb the relaxation vector
    stt%relaxationVector(:,en) = p_relax
    call grainDeformation(pF,avgF,ho,en)                                                            ! rain deformation from perturbed state
    call stressPenalty(pR,DevNull,      avgF,pF,ho,en)                                              ! stress penalty due to interface mismatch from perturbed state
    call volumePenalty(pD,devNull(1,1), avgF,pF,nGrain)                                             ! stress penalty due to volume discrepancy from perturbed state

!--------------------------------------------------------------------------------------------------
! computing the global stress residual array from the perturbed state
    p_resid = 0.0_pREAL
    do iNum = 1,nIntFaceTot
      faceID = interface1to4(iNum,param(ho)%N_constituents)                                         ! identifying the interface ID in local coordinate system (4-dimensional index)

!--------------------------------------------------------------------------------------------------
! identify the left/bottom/back grain (-|N)
      iGr3N = faceID(2:4)                                                                           ! identify the grain ID in local coordinate system (3-dimensional index)
      iGrN = grain3to1(iGr3N,param(ho)%N_constituents)                                              ! translate the local grain ID into global coordinate system (1-dimensional index)
      intFaceN = getInterface(2*faceID(1),iGr3N)                                                    ! identify the interface ID of the grain
      normN = interfaceNormal(intFaceN,ho,en)

!--------------------------------------------------------------------------------------------------
! identify the right/up/front grain (+|P)
      iGr3P = iGr3N
      iGr3P(faceID(1)) = iGr3N(faceID(1))+1                                                         ! identify the grain ID in local coordinate system (3-dimensional index)
      iGrP = grain3to1(iGr3P,param(ho)%N_constituents)                                              ! translate the local grain ID into global coordinate system (1-dimensional index)
      intFaceP = getInterface(2*faceID(1)-1,iGr3P)                                                  ! identify the interface ID of the grain
      normP = interfaceNormal(intFaceP,ho,en)

!--------------------------------------------------------------------------------------------------
! compute the residual stress (contribution of mismatch and volume penalties) from perturbed state
! at all interfaces
      do i = 1,3; do j = 1,3
        p_resid(i+3*(iNum-1)) = p_resid(i+3*(iNum-1)) + (pR(i,j,iGrP) - R(i,j,iGrP))*normP(j) &
                                                      + (pR(i,j,iGrN) - R(i,j,iGrN))*normN(j) &
                                                      + (pD(i,j,iGrP) - D(i,j,iGrP))*normP(j) &
                                                      + (pD(i,j,iGrN) - D(i,j,iGrN))*normN(j)
      end do; end do
    end do
    pmatrix(:,ipert) = p_resid/num%pPert
  end do


!--------------------------------------------------------------------------------------------------
! ... of the numerical viscosity traction "rmatrix"
  allocate(rmatrix(3*nIntFaceTot,3*nIntFaceTot),source=0.0_pREAL)
  do i=1,3*nIntFaceTot
    rmatrix(i,i) = num%viscModus*num%viscPower/(num%refRelaxRate*dt)* &                             ! tangent due to numerical viscosity traction appears
                   (abs(drelax(i))/(num%refRelaxRate*dt))**(num%viscPower - 1.0_pREAL)              ! only in the main diagonal term
  end do


!--------------------------------------------------------------------------------------------------
! The overall Jacobian matrix summarizing contributions of smatrix, pmatrix, rmatrix
  allocate(jmatrix(3*nIntFaceTot,3*nIntFaceTot)); jmatrix = smatrix + pmatrix + rmatrix

!--------------------------------------------------------------------------------------------------
! computing the update of the state variable (relaxation vectors) using the Jacobian matrix
  allocate(jnverse(3*nIntFaceTot,3*nIntFaceTot),source=0.0_pREAL)
  call math_invert(jnverse,error,jmatrix)

!--------------------------------------------------------------------------------------------------
! calculate the state update (global relaxation vectors) for the next Newton-Raphson iteration
  drelax = 0.0_pREAL
  do i = 1,3*nIntFaceTot;do j = 1,3*nIntFaceTot
    drelax(i) = drelax(i) - jnverse(i,j)*resid(j)                                                   ! Calculate the correction for the state variable
  end do; end do
  stt%relaxationVector(:,en) = relax + drelax                                                       ! Updateing the state variable for the next iteration
  if (any(abs(drelax) > num%maxdRelax)) then                                                        ! Forcing cutback when the incremental change of relaxation vector becomes too large
    doneAndHappy = [.true.,.false.]
    !$OMP CRITICAL (write2out)
    print'(a,i3,a,i3,a)',' RGC_updateState: enforces cutback'
    print'(a,e15.8)',' due to large relaxation change = ',maxval(abs(drelax))
    flush(IO_STDOUT)
    !$OMP END CRITICAL (write2out)
  end if

  end associate

  contains
  !------------------------------------------------------------------------------------------------
  !> @brief calculate stress-like penalty due to deformation mismatch
  !------------------------------------------------------------------------------------------------
  subroutine stressPenalty(rPen,nMis,avgF,fDef,ho,en)

    real(pREAL),   dimension (:,:,:), intent(out) :: rPen                                           !< stress-like penalty
    real(pREAL),   dimension (:,:),   intent(out) :: nMis                                           !< total amount of mismatch

    real(pREAL),   dimension (:,:,:), intent(in)  :: fDef                                           !< deformation gradients
    real(pREAL),   dimension (3,3),   intent(in)  :: avgF                                           !< initial effective stretch tensor
    integer,                          intent(in)  :: ho, en

    integer, dimension (4)   :: intFace
    integer, dimension (3)   :: iGrain3,iGNghb3,nGDim
    real(pREAL),   dimension (3,3) :: gDef,nDef
    real(pREAL),   dimension (3)   :: nVect,surfCorr
    integer :: iGrain,iGNghb,iFace,i,j,k,l
    real(pREAL) :: muGrain,muGNghb,nDefNorm
    real(pREAL), parameter  :: &
      nDefToler = 1.0e-10_pREAL, &
      b = 2.5e-10_pREAL                                                                             ! Length of Burgers vector

    nGDim = param(ho)%N_constituents
    rPen = 0.0_pREAL
    nMis = 0.0_pREAL

    !----------------------------------------------------------------------------------------------
    ! get the correction factor the modulus of penalty stress representing the evolution of area of
    ! the interfaces due to deformations

    surfCorr = surfaceCorrection(avgF,ho,en)

    associate(prm => param(ho))

   !-----------------------------------------------------------------------------------------------
   ! computing the mismatch and penalty stress tensor of all grains
   grainLoop: do iGrain = 1,product(prm%N_constituents)
     muGrain = equivalentMu(iGrain,ce)
     iGrain3 = grain1to3(iGrain,prm%N_constituents)                                                 ! get the grain ID in local 3-dimensional index (x,y,z)-position

     interfaceLoop: do iFace = 1,6
       intFace = getInterface(iFace,iGrain3)                                                        ! get the 4-dimensional index of the interface in local numbering system of the grain
       nVect = interfaceNormal(intFace,ho,en)
       iGNghb3 = iGrain3                                                                            ! identify the neighboring grain across the interface
       iGNghb3(abs(intFace(1))) = iGNghb3(abs(intFace(1))) &
                                + int(real(intFace(1),pREAL)/real(abs(intFace(1)),pREAL))
       where(iGNghb3 < 1)    iGNghb3 = nGDim
       where(iGNghb3 >nGDim) iGNghb3 = 1
       iGNghb  = grain3to1(iGNghb3,prm%N_constituents)                                              ! get the ID of the neighboring grain
       muGNghb = equivalentMu(iGNghb,ce)
       gDef = 0.5_pREAL*(fDef(1:3,1:3,iGNghb) - fDef(1:3,1:3,iGrain))                               ! difference/jump in deformation gradeint across the neighbor

       !-------------------------------------------------------------------------------------------
       ! compute the mismatch tensor of all interfaces
       nDefNorm = 0.0_pREAL
       nDef = 0.0_pREAL
       do i = 1,3; do j = 1,3
         do k = 1,3; do l = 1,3
           nDef(i,j) = nDef(i,j) - nVect(k)*gDef(i,l)*math_LeviCivita(j,k,l)                        ! compute the interface mismatch tensor from the jump of deformation gradient
         end do; end do
         nDefNorm = nDefNorm + nDef(i,j)**2                                                         ! compute the norm of the mismatch tensor
       end do; end do
       nDefNorm = max(nDefToler,sqrt(nDefNorm))                                                     ! approximation to zero mismatch if mismatch is zero (singularity)
       nMis(abs(intFace(1)),iGrain) = nMis(abs(intFace(1)),iGrain) + nDefNorm                       ! total amount of mismatch experienced by the grain (at all six interfaces)


       !-------------------------------------------------------------------------------------------
       ! compute the stress penalty of all interfaces
       do i = 1,3; do j = 1,3; do k = 1,3; do l = 1,3
         rPen(i,j,iGrain) = rPen(i,j,iGrain) + 0.5_pREAL*(muGrain*b + muGNghb*b)*prm%xi_alpha &
                                                *surfCorr(abs(intFace(1)))/prm%D_alpha(abs(intFace(1))) &
                                                *cosh(prm%c_alpha*nDefNorm) &
                                                *0.5_pREAL*nVect(l)*nDef(i,k)/nDefNorm*math_LeviCivita(k,l,j) &
                                                *tanh(nDefNorm/num%xSmoo)
       end do; end do;end do; end do
     end do interfaceLoop


   end do grainLoop

   end associate

  end subroutine stressPenalty


  !------------------------------------------------------------------------------------------------
  !> @brief calculate stress-like penalty due to volume discrepancy
  !------------------------------------------------------------------------------------------------
  subroutine volumePenalty(vPen,vDiscrep,fAvg,fDef,nGrain)

    real(pREAL), dimension (:,:,:), intent(out) :: vPen                                             ! stress-like penalty due to volume
    real(pREAL),                    intent(out) :: vDiscrep                                         ! total volume discrepancy

    real(pREAL), dimension (:,:,:), intent(in)  :: fDef                                             ! deformation gradients
    real(pREAL), dimension (3,3),   intent(in)  :: fAvg                                             ! overall deformation gradient
    integer,                        intent(in) :: &
      Ngrain

    real(pREAL), dimension(size(vPen,3)) :: gVol
    integer :: i

    !----------------------------------------------------------------------------------------------
    ! compute the volumes of grains and of cluster
    vDiscrep = math_det33(fAvg)                                                                     ! compute the volume of the cluster
    do i = 1,nGrain
      gVol(i) = math_det33(fDef(1:3,1:3,i))                                                         ! compute the volume of individual grains
      vDiscrep     = vDiscrep - gVol(i)/real(nGrain,pREAL)                                          ! calculate the difference/dicrepancy between
                                                                                                    ! the volume of the cluster and the the total volume of grains
    end do

    !----------------------------------------------------------------------------------------------
    ! calculate the stress and penalty due to volume discrepancy
    vPen      = 0.0_pREAL
    do i = 1,nGrain
      vPen(:,:,i) = -real(nGrain,pREAL)**(-1)*num%volDiscrMod*num%volDiscrPow/num%maxVolDiscr &
                  * sign((abs(vDiscrep)/num%maxVolDiscr)**(num%volDiscrPow - 1.0_pREAL),vDiscrep) &
                  * gVol(i)*transpose(math_inv33(fDef(:,:,i)))
    end do

  end subroutine volumePenalty


  !--------------------------------------------------------------------------------------------------
  !> @brief compute the correction factor accouted for surface evolution (area change) due to
  ! deformation
  !--------------------------------------------------------------------------------------------------
  function surfaceCorrection(avgF,ho,en)

    real(pREAL), dimension(3)               :: surfaceCorrection

    real(pREAL), dimension(3,3), intent(in) :: avgF                                                 !< average F
    integer,                     intent(in) :: &
      ho, &
      en
    real(pREAL), dimension(3,3)             :: invC
    real(pREAL), dimension(3)               :: nVect
    real(pREAL)  :: detF
    integer :: i,j,iBase
    logical     ::  error

    call math_invert33(invC,detF,error,matmul(transpose(avgF),avgF))

    surfaceCorrection = 0.0_pREAL
    do iBase = 1,3
      nVect = interfaceNormal([iBase,1,1,1],ho,en)
      do i = 1,3; do j = 1,3
        surfaceCorrection(iBase) = surfaceCorrection(iBase) + invC(i,j)*nVect(i)*nVect(j)           ! compute the component of (the inverse of) the stretch in the direction of the normal
      end do; end do
      surfaceCorrection(iBase) = sqrt(surfaceCorrection(iBase))*detF                                ! get the surface correction factor (area contraction/enlargement)
    end do

  end function surfaceCorrection


  !-------------------------------------------------------------------------------------------------
  !> @brief compute the equivalent shear and bulk moduli from the elasticity tensor
  !-------------------------------------------------------------------------------------------------
  real(pREAL) function equivalentMu(co,ce)

    integer, intent(in)    :: &
      co,&
      ce

    real(pREAL), dimension(6,6) :: C

    C = phase_homogenizedC66(material_ID_phase(co,ce),material_entry_phase(co,ce))                  ! damage not included!

    equivalentMu = crystal_isotropic_mu(C,'isostrain')

  end function equivalentMu


  !-------------------------------------------------------------------------------------------------
  !> @brief calculating the grain deformation gradient (the same with
  ! homogenization_RGC_partitionDeformation, but used only for perturbation scheme)
  !-------------------------------------------------------------------------------------------------
  subroutine grainDeformation(F, avgF, ho, en)

    real(pREAL),   dimension(:,:,:), intent(out) :: F                                               !< partitioned F  per grain

    real(pREAL),   dimension(:,:),   intent(in)  :: avgF                                            !< averaged F
    integer,                         intent(in)  :: &
      ho, &
      en

    real(pREAL),   dimension(3) :: aVect,nVect
    integer,       dimension(4) :: intFace
    integer,       dimension(3) :: iGrain3
    integer :: iGrain,iFace,i,j

    !-----------------------------------------------------------------------------------------------
    ! compute the deformation gradient of individual grains due to relaxations

    associate (prm => param(ho))

    F = 0.0_pREAL
    do iGrain = 1,product(prm%N_constituents)
      iGrain3 = grain1to3(iGrain,prm%N_constituents)
      do iFace = 1,6
        intFace = getInterface(iFace,iGrain3)
        aVect   = relaxationVector(intFace,ho,en)
        nVect   = interfaceNormal(intFace,ho,en)
        forall (i=1:3,j=1:3) &
          F(i,j,iGrain) = F(i,j,iGrain) + aVect(i)*nVect(j)                                         ! effective relaxations
      end do
      F(1:3,1:3,iGrain) = F(1:3,1:3,iGrain) + avgF                                                  ! relaxed deformation gradient
    end do

    end associate

  end subroutine grainDeformation

end function RGC_updateState


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
module subroutine RGC_result(ho,group)

  integer,          intent(in) :: ho
  character(len=*), intent(in) :: group

  integer :: o

  associate(stt => state(ho), dst => dependentState(ho), prm => param(ho))
    outputsLoop: do o = 1,size(prm%output)
      select case(trim(prm%output(o)))
        case('M')
          call result_writeDataset(dst%mismatch,group,trim(prm%output(o)), &
                                   'average mismatch tensor','1')
        case('Delta_V')
          call result_writeDataset(dst%volumeDiscrepancy,group,trim(prm%output(o)), &
                                 'volume discrepancy','m³')
        case('max_dot_a')
          call result_writeDataset(dst%relaxationrate_max,group,trim(prm%output(o)), &
                                 'maximum relaxation rate','m/s')
        case('avg_dot_a')
          call result_writeDataset(dst%relaxationrate_avg,group,trim(prm%output(o)), &
                                   'average relaxation rate','m/s')
      end select
    end do outputsLoop
  end associate

end subroutine RGC_result


!--------------------------------------------------------------------------------------------------
!> @brief collect relaxation vectors of an interface
!--------------------------------------------------------------------------------------------------
pure function relaxationVector(intFace,ho,en)

  real(pREAL), dimension (3)            :: relaxationVector

  integer,                   intent(in) :: ho,en
  integer,     dimension(4), intent(in) :: intFace                                                  !< set of interface ID in 4D array (normal and position)

  integer :: iNum

!--------------------------------------------------------------------------------------------------
! collect the interface relaxation vector from the global state array

  associate (prm => param(ho), &
             stt => state(ho))

  iNum = interface4to1(intFace,prm%N_constituents)                                                  ! identify the position of the interface in global state array
  if (iNum > 0) then
    relaxationVector = stt%relaxationVector((3*iNum-2):(3*iNum),en)
  else
    relaxationVector = 0.0_pREAL
  end if

  end associate

end function relaxationVector


!--------------------------------------------------------------------------------------------------
!> @brief identify the normal of an interface
!--------------------------------------------------------------------------------------------------
pure function interfaceNormal(intFace,ho,en) result(n)

  real(pREAL), dimension(3)             :: n
  integer,     dimension(4), intent(in) :: intFace                                                  !< interface ID in 4D array (normal and position)
  integer,                   intent(in) :: &
    ho, &
    en


  associate (dst => dependentState(ho))

    n = 0.0_pREAL
    n(abs(intFace(1))) = real(intFace(1)/abs(intFace(1)),pREAL)                                     ! get the normal vector w.r.t. cluster axis

    n = matmul(dst%orientation(1:3,1:3,en),n)                                                       ! map the normal vector into sample coordinate system (basis)

  end associate

end function interfaceNormal


!--------------------------------------------------------------------------------------------------
!> @brief collect six faces of a grain in 4D (normal and position)
!--------------------------------------------------------------------------------------------------
pure function getInterface(iFace,iGrain3) result(i)

  integer, dimension(4)             :: i
  integer, dimension(3), intent(in) :: iGrain3                                                      !< grain ID in 3D array
  integer,               intent(in) :: iFace                                                        !< face index (1..6) mapped like (-e1,-e2,-e3,+e1,+e2,+e3) or iDir = (-1,-2,-3,1,2,3)

  integer :: iDir                                                                                   !< direction of interface normal


  iDir = (int(real(iFace-1,pREAL)/2.0_pREAL)+1)*(-1)**iFace
  i = [iDir,iGrain3]
  if (iDir < 0) i(1-iDir) = i(1-iDir)-1                                       ! to have a correlation with coordinate/position in real space

end function getInterface


!--------------------------------------------------------------------------------------------------
!> @brief map grain ID from in 1D (global array) to in 3D (local position)
!--------------------------------------------------------------------------------------------------
pure function grain1to3(grain1,nGDim)

  integer, dimension(3)             :: grain1to3

  integer,               intent(in) :: grain1                                                       !< grain ID in 1D array
  integer, dimension(3), intent(in) :: nGDim

  grain1to3 = 1 + [mod((grain1-1), nGDim(1)), &
                   mod((grain1-1)/ nGDim(1),nGDim(2)), &
                       (grain1-1)/(nGDim(1)*nGDim(2))]

end function grain1to3


!--------------------------------------------------------------------------------------------------
!> @brief map grain ID from in 3D (local position) to in 1D (global array)
!--------------------------------------------------------------------------------------------------
integer pure function grain3to1(grain3,nGDim)

  integer, dimension(3), intent(in) :: grain3                                                       !< grain ID in 3D array (pos.x,pos.y,pos.z)
  integer, dimension(3), intent(in) :: nGDim

  grain3to1 = grain3(1) &
            + nGDim(1)*(grain3(2)-1) &
            + nGDim(1)*nGDim(2)*(grain3(3)-1)

end function grain3to1


!--------------------------------------------------------------------------------------------------
!> @brief maps interface ID from 4D (normal and local position) into 1D (global array)
!--------------------------------------------------------------------------------------------------
integer pure function interface4to1(iFace4D, nGDim)

  integer, dimension(4), intent(in) :: iFace4D                                                      !< interface ID in 4D array (n.dir,pos.x,pos.y,pos.z)
  integer, dimension(3), intent(in) :: nGDim


  select case(abs(iFace4D(1)))

    case(1)
      if ((iFace4D(2) == 0) .or. (iFace4D(2) == nGDim(1))) then
        interface4to1 = 0
      else
        interface4to1 = iFace4D(3) + nGDim(2)*(iFace4D(4)-1) &
                      + nGDim(2)*nGDim(3)*(iFace4D(2)-1)
      end if

    case(2)
      if ((iFace4D(3) == 0) .or. (iFace4D(3) == nGDim(2))) then
        interface4to1 = 0
      else
        interface4to1 = iFace4D(4) + nGDim(3)*(iFace4D(2)-1) &
                      + nGDim(3)*nGDim(1)*(iFace4D(3)-1) &
                      + (nGDim(1)-1)*nGDim(2)*nGDim(3)                                              ! total # of interfaces normal || e1
      end if

    case(3)
      if ((iFace4D(4) == 0) .or. (iFace4D(4) == nGDim(3))) then
        interface4to1 = 0
      else
        interface4to1 = iFace4D(2) + nGDim(1)*(iFace4D(3)-1) &
                      + nGDim(1)*nGDim(2)*(iFace4D(4)-1) &
                      + (nGDim(1)-1)*nGDim(2)*nGDim(3) &                                            ! total # of interfaces normal || e1
                      + nGDim(1)*(nGDim(2)-1)*nGDim(3)                                              ! total # of interfaces normal || e2
      end if

    case default
      interface4to1 = -1

  end select

end function interface4to1


!--------------------------------------------------------------------------------------------------
!> @brief maps interface ID from 1D (global array) into 4D (normal and local position)
!--------------------------------------------------------------------------------------------------
pure function interface1to4(iFace1D, nGDim)

  integer, dimension(4)             :: interface1to4

  integer,               intent(in) :: iFace1D                                                      !< interface ID in 1D array
  integer, dimension(3), intent(in) :: nGDim
  integer, dimension(3)             :: nIntFace

!--------------------------------------------------------------------------------------------------
! compute the total number of interfaces, which ...
  nIntFace = [(nGDim(1)-1)*nGDim(2)*nGDim(3), &                                                     ! ... normal || e1
              nGDim(1)*(nGDim(2)-1)*nGDim(3), &                                                     ! ... normal || e2
              nGDim(1)*nGDim(2)*(nGDim(3)-1)]                                                       ! ... normal || e3

!--------------------------------------------------------------------------------------------------
! get the corresponding interface ID in 4D (normal and local position)
  if (iFace1D > 0 .and. iFace1D <= nIntFace(1)) then                                                ! interface with normal || e1
    interface1to4(1) = 1
    interface1to4(3) = mod((iFace1D-1),nGDim(2))+1
    interface1to4(4) = mod(int(real(iFace1D-1,pREAL)/real(nGDim(2),pREAL)),nGDim(3))+1
    interface1to4(2) = int(real(iFace1D-1,pREAL)/real(nGDim(2),pREAL)/real(nGDim(3),pREAL))+1
  elseif (iFace1D > nIntFace(1) .and. iFace1D <= (nIntFace(2) + nIntFace(1))) then                  ! interface with normal || e2
    interface1to4(1) = 2
    interface1to4(4) = mod((iFace1D-nIntFace(1)-1),nGDim(3))+1
    interface1to4(2) = mod(int(real(iFace1D-nIntFace(1)-1,pREAL)/real(nGDim(3),pREAL)),nGDim(1))+1
    interface1to4(3) = int(real(iFace1D-nIntFace(1)-1,pREAL)/real(nGDim(3),pREAL)/real(nGDim(1),pREAL))+1
  elseif (iFace1D > nIntFace(2) + nIntFace(1) .and. iFace1D <= (nIntFace(3) + nIntFace(2) + nIntFace(1))) then ! interface with normal || e3
    interface1to4(1) = 3
    interface1to4(2) = mod((iFace1D-nIntFace(2)-nIntFace(1)-1),nGDim(1))+1
    interface1to4(3) = mod(int(real(iFace1D-nIntFace(2)-nIntFace(1)-1,pREAL)/real(nGDim(1),pREAL)),nGDim(2))+1
    interface1to4(4) = int(real(iFace1D-nIntFace(2)-nIntFace(1)-1,pREAL)/real(nGDim(1),pREAL)/real(nGDim(2),pREAL))+1
  end if

end function interface1to4


end submodule RGC
