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
  use lattice

  type :: tParameters
    integer, dimension(:), allocatable :: &
      N_constituents
    real(pReal) :: &
      xi_alpha, &
      c_Alpha
    real(pReal), dimension(:), allocatable :: &
      D_alpha, &
      a_g
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type :: tRGCstate
    real(pReal), pointer,     dimension(:,:) :: &
      relaxationVector
  end type tRGCstate

  type :: tRGCdependentState
    real(pReal), allocatable,     dimension(:) :: &
      volumeDiscrepancy, &
      relaxationRate_avg, &
      relaxationRate_max
    real(pReal), allocatable,     dimension(:,:) :: &
      mismatch
    real(pReal), allocatable,     dimension(:,:,:) :: &
      orientation
  end type tRGCdependentState

  type :: tNumerics_RGC
    real(pReal) :: &
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
module subroutine RGC_init(num_homogMech)

  class(tNode), pointer, intent(in) :: &
    num_homogMech                                                                                   !< pointer to mechanical homogenization numerics data

  integer :: &
    ho, &
    Nmembers, &
    sizeState, nIntFaceTot

  class (tNode), pointer :: &
    num_RGC, &                                                                                      ! pointer to RGC numerics data
    material_homogenization, &
    homog, &
    homogMech

  print'(/,a)', ' <<<+-  homogenization:mechanical:RGC init  -+>>>'

  print'(a,i0)', ' # homogenizations: ',count(homogenization_type == HOMOGENIZATION_RGC_ID)
  flush(IO_STDOUT)

  print*, 'D.D. Tjahjanto et al., International Journal of Material Forming 2(1):939–942, 2009'
  print*, 'https://doi.org/10.1007/s12289-009-0619-1'//IO_EOL

  print*, 'D.D. Tjahjanto et al., Modelling and Simulation in Materials Science and Engineering 18:015006, 2010'
  print*, 'https://doi.org/10.1088/0965-0393/18/1/015006'//IO_EOL


  material_homogenization => config_material%get('homogenization')
  allocate(param(material_homogenization%length))
  allocate(state(material_homogenization%length))
  allocate(state0(material_homogenization%length))
  allocate(dependentState(material_homogenization%length))

  num_RGC => num_homogMech%get('RGC',defaultVal=emptyDict)

  num%atol         =  num_RGC%get_asFloat('atol',              defaultVal=1.0e+4_pReal)
  num%rtol         =  num_RGC%get_asFloat('rtol',              defaultVal=1.0e-3_pReal)
  num%absMax       =  num_RGC%get_asFloat('amax',              defaultVal=1.0e+10_pReal)
  num%relMax       =  num_RGC%get_asFloat('rmax',              defaultVal=1.0e+2_pReal)
  num%pPert        =  num_RGC%get_asFloat('perturbpenalty',    defaultVal=1.0e-7_pReal)
  num%xSmoo        =  num_RGC%get_asFloat('relvantmismatch',   defaultVal=1.0e-5_pReal)
  num%viscPower    =  num_RGC%get_asFloat('viscositypower',    defaultVal=1.0e+0_pReal)
  num%viscModus    =  num_RGC%get_asFloat('viscositymodulus',  defaultVal=0.0e+0_pReal)
  num%refRelaxRate =  num_RGC%get_asFloat('refrelaxationrate', defaultVal=1.0e-3_pReal)
  num%maxdRelax    =  num_RGC%get_asFloat('maxrelaxationrate', defaultVal=1.0e+0_pReal)
  num%maxVolDiscr  =  num_RGC%get_asFloat('maxvoldiscrepancy', defaultVal=1.0e-5_pReal)
  num%volDiscrMod  =  num_RGC%get_asFloat('voldiscrepancymod', defaultVal=1.0e+12_pReal)
  num%volDiscrPow  =  num_RGC%get_asFloat('dicrepancypower',   defaultVal=5.0_pReal)

  if (num%atol <= 0.0_pReal)         call IO_error(301,ext_msg='absTol_RGC')
  if (num%rtol <= 0.0_pReal)         call IO_error(301,ext_msg='relTol_RGC')
  if (num%absMax <= 0.0_pReal)       call IO_error(301,ext_msg='absMax_RGC')
  if (num%relMax <= 0.0_pReal)       call IO_error(301,ext_msg='relMax_RGC')
  if (num%pPert <= 0.0_pReal)        call IO_error(301,ext_msg='pPert_RGC')
  if (num%xSmoo <= 0.0_pReal)        call IO_error(301,ext_msg='xSmoo_RGC')
  if (num%viscPower < 0.0_pReal)     call IO_error(301,ext_msg='viscPower_RGC')
  if (num%viscModus < 0.0_pReal)     call IO_error(301,ext_msg='viscModus_RGC')
  if (num%refRelaxRate <= 0.0_pReal) call IO_error(301,ext_msg='refRelaxRate_RGC')
  if (num%maxdRelax <= 0.0_pReal)    call IO_error(301,ext_msg='maxdRelax_RGC')
  if (num%maxVolDiscr <= 0.0_pReal)  call IO_error(301,ext_msg='maxVolDiscr_RGC')
  if (num%volDiscrMod < 0.0_pReal)   call IO_error(301,ext_msg='volDiscrMod_RGC')
  if (num%volDiscrPow <= 0.0_pReal)  call IO_error(301,ext_msg='volDiscrPw_RGC')


  do ho = 1, size(homogenization_type)
    if (homogenization_type(ho) /= HOMOGENIZATION_RGC_ID) cycle
    homog => material_homogenization%get(ho)
    homogMech => homog%get('mechanical')
    associate(prm => param(ho), &
              stt => state(ho), &
              st0 => state0(ho), &
              dst => dependentState(ho))

#if defined (__GFORTRAN__)
    prm%output = output_as1dString(homogMech)
#else
    prm%output = homogMech%get_as1dString('output',defaultVal=emptyStringArray)
#endif

    prm%N_constituents = homogMech%get_as1dInt('cluster_size',requiredSize=3)
    if (homogenization_Nconstituents(ho) /= product(prm%N_constituents)) &
      call IO_error(211,ext_msg='N_constituents (RGC)')

    prm%xi_alpha = homogMech%get_asFloat('xi_alpha')
    prm%c_alpha  = homogMech%get_asFloat('c_alpha')

    prm%D_alpha  = homogMech%get_as1dFloat('D_alpha', requiredSize=3)
    prm%a_g      = homogMech%get_as1dFloat('a_g',     requiredSize=3)

    Nmembers = count(material_homogenizationID == ho)
    nIntFaceTot = 3*(  (prm%N_constituents(1)-1)*prm%N_constituents(2)*prm%N_constituents(3) &
                      + prm%N_constituents(1)*(prm%N_constituents(2)-1)*prm%N_constituents(3) &
                      + prm%N_constituents(1)*prm%N_constituents(2)*(prm%N_constituents(3)-1))
    sizeState = nIntFaceTot

    homogState(ho)%sizeState = sizeState
    allocate(homogState(ho)%state0   (sizeState,Nmembers), source=0.0_pReal)
    allocate(homogState(ho)%state    (sizeState,Nmembers), source=0.0_pReal)

    stt%relaxationVector   => homogState(ho)%state(1:nIntFaceTot,:)
    st0%relaxationVector   => homogState(ho)%state0(1:nIntFaceTot,:)

    allocate(dst%volumeDiscrepancy(   Nmembers), source=0.0_pReal)
    allocate(dst%relaxationRate_avg(  Nmembers), source=0.0_pReal)
    allocate(dst%relaxationRate_max(  Nmembers), source=0.0_pReal)
    allocate(dst%mismatch(          3,Nmembers), source=0.0_pReal)

!--------------------------------------------------------------------------------------------------
! assigning cluster orientations
    dependentState(ho)%orientation = spread(eu2om(prm%a_g*inRad),3,Nmembers)
    !dst%orientation = spread(eu2om(prm%a_g*inRad),3,Nmembers) ifort version 18.0.1 crashes (for whatever reason)

    end associate

  enddo

end subroutine RGC_init


!--------------------------------------------------------------------------------------------------
!> @brief partitions the deformation gradient onto the constituents
!--------------------------------------------------------------------------------------------------
module subroutine RGC_partitionDeformation(F,avgF,ce)

  real(pReal),   dimension (:,:,:), intent(out) :: F                                                !< partitioned F  per grain

  real(pReal),   dimension (3,3),   intent(in)  :: avgF                                             !< averaged F
  integer,                          intent(in)  :: &
    ce

  real(pReal), dimension(3) :: aVect,nVect
  integer,     dimension(4) :: intFace
  integer,     dimension(3) :: iGrain3
  integer ::  iGrain,iFace,i,j,ho,en

  associate(prm => param(material_homogenizationID(ce)))

  ho = material_homogenizationID(ce)
  en = material_homogenizationEntry(ce)
!--------------------------------------------------------------------------------------------------
! compute the deformation gradient of individual grains due to relaxations
  F = 0.0_pReal
  do iGrain = 1,product(prm%N_constituents)
    iGrain3 = grain1to3(iGrain,prm%N_constituents)
    do iFace = 1,6
      intFace = getInterface(iFace,iGrain3)                                                         ! identifying 6 interfaces of each grain
      aVect = relaxationVector(intFace,ho,en)                                                       ! get the relaxation vectors for each interface from global relaxation vector array
      nVect = interfaceNormal(intFace,ho,en)
      forall (i=1:3,j=1:3) &
        F(i,j,iGrain) = F(i,j,iGrain) + aVect(i)*nVect(j)                                           ! calculating deformation relaxations due to interface relaxation
    enddo
    F(1:3,1:3,iGrain) = F(1:3,1:3,iGrain) + avgF                                                    ! resulting relaxed deformation gradient
  enddo

  end associate

end subroutine RGC_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief update the internal state of the homogenization scheme and tell whether "done" and
! "happy" with result
!--------------------------------------------------------------------------------------------------
module function RGC_updateState(P,F,avgF,dt,dPdF,ce) result(doneAndHappy)
      logical, dimension(2) :: doneAndHappy
      real(pReal), dimension(:,:,:),     intent(in)    :: &
        P,&                                                                                         !< partitioned stresses
        F                                                                                           !< partitioned deformation gradients
      real(pReal), dimension(:,:,:,:,:), intent(in) :: dPdF                                         !< partitioned stiffnesses
      real(pReal), dimension(3,3),       intent(in) :: avgF                                         !< average F
      real(pReal),                       intent(in) :: dt                                           !< time increment
      integer,                           intent(in) :: &
        ce                                                                                          !< cell

  integer, dimension(4) :: intFaceN,intFaceP,faceID
  integer, dimension(3) :: nGDim,iGr3N,iGr3P
  integer :: ho,iNum,i,j,nIntFaceTot,iGrN,iGrP,iMun,iFace,k,l,ipert,nGrain, en
  real(pReal), dimension(3,3,size(P,3)) :: R,pF,pR,D,pD
  real(pReal), dimension(3,size(P,3))   :: NN,devNull
  real(pReal), dimension(3)             :: normP,normN,mornP,mornN
  real(pReal) :: residMax,stresMax
  logical :: error
  real(pReal), dimension(:,:), allocatable :: tract,jmatrix,jnverse,smatrix,pmatrix,rmatrix
  real(pReal), dimension(:), allocatable   :: resid,relax,p_relax,p_resid,drelax

  zeroTimeStep: if(dEq0(dt)) then
    doneAndHappy = .true.                                                                   ! pretend everything is fine and return
    return
  endif zeroTimeStep

  ho  = material_homogenizationID(ce)
  en = material_homogenizationEntry(ce)

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
  allocate(resid(3*nIntFaceTot), source=0.0_pReal)
  allocate(tract(nIntFaceTot,3), source=0.0_pReal)
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
      enddo
    enddo

  enddo

!--------------------------------------------------------------------------------------------------
! convergence check for stress residual
  stresMax = maxval(abs(P))                                                                         ! get the maximum of first Piola-Kirchhoff (material) stress
  residMax = maxval(abs(tract))                                                                     ! get the maximum of the residual

  doneAndHappy = .false.

!--------------------------------------------------------------------------------------------------
!  If convergence reached => done and happy
  if (residMax < num%rtol*stresMax .or. residMax < num%atol) then
    doneAndHappy = .true.

    dst%mismatch(1:3,en)       = sum(NN,2)/real(nGrain,pReal)
    dst%relaxationRate_avg(en) = sum(abs(drelax))/dt/real(3*nIntFaceTot,pReal)
    dst%relaxationRate_max(en) = maxval(abs(drelax))/dt

    return

!--------------------------------------------------------------------------------------------------
! if residual blows-up => done but unhappy
  elseif (residMax > num%relMax*stresMax .or. residMax > num%absMax) then                           ! try to restart when residual blows up exceeding maximum bound
    doneAndHappy = [.true.,.false.]                                                         ! with direct cut-back
   return
  endif

!---------------------------------------------------------------------------------------------------
! construct the global Jacobian matrix for updating the global relaxation vector array when
! convergence is not yet reached ...

!--------------------------------------------------------------------------------------------------
! ... of the constitutive stress tangent, assembled from dPdF or material constitutive model "smatrix"
  allocate(smatrix(3*nIntFaceTot,3*nIntFaceTot), source=0.0_pReal)
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
        enddo;enddo;enddo;enddo
! projecting the material tangent dPdF into the interface
! to obtain the Jacobian matrix contribution of dPdF
      endif
    enddo

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
        enddo;enddo;enddo;enddo
      endif
    enddo
  enddo

!--------------------------------------------------------------------------------------------------
! ... of the stress penalty tangent (mismatch penalty and volume penalty, computed using numerical
! perturbation method) "pmatrix"
  allocate(pmatrix(3*nIntFaceTot,3*nIntFaceTot), source=0.0_pReal)
  allocate(p_relax(3*nIntFaceTot),               source=0.0_pReal)
  allocate(p_resid(3*nIntFaceTot),               source=0.0_pReal)

  do ipert = 1,3*nIntFaceTot
    p_relax = relax
    p_relax(ipert) = relax(ipert) + num%pPert                                                       ! perturb the relaxation vector
    stt%relaxationVector(:,en) = p_relax
    call grainDeformation(pF,avgF,ho,en)                                                            ! rain deformation from perturbed state
    call stressPenalty(pR,DevNull,      avgF,pF,ho,en)                                              ! stress penalty due to interface mismatch from perturbed state
    call volumePenalty(pD,devNull(1,1), avgF,pF,nGrain)                                             ! stress penalty due to volume discrepancy from perturbed state

!--------------------------------------------------------------------------------------------------
! computing the global stress residual array from the perturbed state
    p_resid = 0.0_pReal
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
      enddo; enddo
    enddo
    pmatrix(:,ipert) = p_resid/num%pPert
  enddo


!--------------------------------------------------------------------------------------------------
! ... of the numerical viscosity traction "rmatrix"
  allocate(rmatrix(3*nIntFaceTot,3*nIntFaceTot),source=0.0_pReal)
  do i=1,3*nIntFaceTot
    rmatrix(i,i) = num%viscModus*num%viscPower/(num%refRelaxRate*dt)* &                             ! tangent due to numerical viscosity traction appears
                   (abs(drelax(i))/(num%refRelaxRate*dt))**(num%viscPower - 1.0_pReal)              ! only in the main diagonal term
  enddo


!--------------------------------------------------------------------------------------------------
! The overall Jacobian matrix summarizing contributions of smatrix, pmatrix, rmatrix
  allocate(jmatrix(3*nIntFaceTot,3*nIntFaceTot)); jmatrix = smatrix + pmatrix + rmatrix

!--------------------------------------------------------------------------------------------------
! computing the update of the state variable (relaxation vectors) using the Jacobian matrix
  allocate(jnverse(3*nIntFaceTot,3*nIntFaceTot),source=0.0_pReal)
  call math_invert(jnverse,error,jmatrix)

!--------------------------------------------------------------------------------------------------
! calculate the state update (global relaxation vectors) for the next Newton-Raphson iteration
  drelax = 0.0_pReal
  do i = 1,3*nIntFaceTot;do j = 1,3*nIntFaceTot
    drelax(i) = drelax(i) - jnverse(i,j)*resid(j)                                                   ! Calculate the correction for the state variable
  enddo; enddo
  stt%relaxationVector(:,en) = relax + drelax                                                       ! Updateing the state variable for the next iteration
  if (any(abs(drelax) > num%maxdRelax)) then                                                        ! Forcing cutback when the incremental change of relaxation vector becomes too large
    doneAndHappy = [.true.,.false.]
    !$OMP CRITICAL (write2out)
    print'(a,i3,a,i3,a)',' RGC_updateState: enforces cutback'
    print'(a,e15.8)',' due to large relaxation change = ',maxval(abs(drelax))
    flush(IO_STDOUT)
    !$OMP END CRITICAL (write2out)
  endif

  end associate

  contains
  !------------------------------------------------------------------------------------------------
  !> @brief calculate stress-like penalty due to deformation mismatch
  !------------------------------------------------------------------------------------------------
  subroutine stressPenalty(rPen,nMis,avgF,fDef,ho,en)

    real(pReal),   dimension (:,:,:), intent(out) :: rPen                                           !< stress-like penalty
    real(pReal),   dimension (:,:),   intent(out) :: nMis                                           !< total amount of mismatch

    real(pReal),   dimension (:,:,:), intent(in)  :: fDef                                           !< deformation gradients
    real(pReal),   dimension (3,3),   intent(in)  :: avgF                                           !< initial effective stretch tensor
    integer,                          intent(in)  :: ho, en

    integer, dimension (4)   :: intFace
    integer, dimension (3)   :: iGrain3,iGNghb3,nGDim
    real(pReal),   dimension (3,3) :: gDef,nDef
    real(pReal),   dimension (3)   :: nVect,surfCorr
    integer :: iGrain,iGNghb,iFace,i,j,k,l
    real(pReal) :: muGrain,muGNghb,nDefNorm
    real(pReal), parameter  :: &
      nDefToler = 1.0e-10_pReal, &
      b = 2.5e-10_pReal                                                                             ! Length of Burgers vector

    nGDim = param(ho)%N_constituents
    rPen = 0.0_pReal
    nMis = 0.0_pReal

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
                                + int(real(intFace(1),pReal)/real(abs(intFace(1)),pReal))
       where(iGNghb3 < 1)    iGNghb3 = nGDim
       where(iGNghb3 >nGDim) iGNghb3 = 1
       iGNghb  = grain3to1(iGNghb3,prm%N_constituents)                                              ! get the ID of the neighboring grain
       muGNghb = equivalentMu(iGNghb,ce)
       gDef = 0.5_pReal*(fDef(1:3,1:3,iGNghb) - fDef(1:3,1:3,iGrain))                               ! difference/jump in deformation gradeint across the neighbor

       !-------------------------------------------------------------------------------------------
       ! compute the mismatch tensor of all interfaces
       nDefNorm = 0.0_pReal
       nDef = 0.0_pReal
       do i = 1,3; do j = 1,3
         do k = 1,3; do l = 1,3
           nDef(i,j) = nDef(i,j) - nVect(k)*gDef(i,l)*math_LeviCivita(j,k,l)                        ! compute the interface mismatch tensor from the jump of deformation gradient
         enddo; enddo
         nDefNorm = nDefNorm + nDef(i,j)**2.0_pReal                                                 ! compute the norm of the mismatch tensor
       enddo; enddo
       nDefNorm = max(nDefToler,sqrt(nDefNorm))                                                     ! approximation to zero mismatch if mismatch is zero (singularity)
       nMis(abs(intFace(1)),iGrain) = nMis(abs(intFace(1)),iGrain) + nDefNorm                       ! total amount of mismatch experienced by the grain (at all six interfaces)


       !-------------------------------------------------------------------------------------------
       ! compute the stress penalty of all interfaces
       do i = 1,3; do j = 1,3; do k = 1,3; do l = 1,3
         rPen(i,j,iGrain) = rPen(i,j,iGrain) + 0.5_pReal*(muGrain*b + muGNghb*b)*prm%xi_alpha &
                                                *surfCorr(abs(intFace(1)))/prm%D_alpha(abs(intFace(1))) &
                                                *cosh(prm%c_alpha*nDefNorm) &
                                                *0.5_pReal*nVect(l)*nDef(i,k)/nDefNorm*math_LeviCivita(k,l,j) &
                                                *tanh(nDefNorm/num%xSmoo)
       enddo; enddo;enddo; enddo
     enddo interfaceLoop


   enddo grainLoop

   end associate

  end subroutine stressPenalty


  !------------------------------------------------------------------------------------------------
  !> @brief calculate stress-like penalty due to volume discrepancy
  !------------------------------------------------------------------------------------------------
  subroutine volumePenalty(vPen,vDiscrep,fAvg,fDef,nGrain)

    real(pReal), dimension (:,:,:), intent(out) :: vPen                                             ! stress-like penalty due to volume
    real(pReal),                    intent(out) :: vDiscrep                                         ! total volume discrepancy

    real(pReal), dimension (:,:,:), intent(in)  :: fDef                                             ! deformation gradients
    real(pReal), dimension (3,3),   intent(in)  :: fAvg                                             ! overall deformation gradient
    integer,                        intent(in) :: &
      Ngrain

    real(pReal), dimension(size(vPen,3)) :: gVol
    integer :: i

    !----------------------------------------------------------------------------------------------
    ! compute the volumes of grains and of cluster
    vDiscrep = math_det33(fAvg)                                                                     ! compute the volume of the cluster
    do i = 1,nGrain
      gVol(i) = math_det33(fDef(1:3,1:3,i))                                                         ! compute the volume of individual grains
      vDiscrep     = vDiscrep - gVol(i)/real(nGrain,pReal)                                          ! calculate the difference/dicrepancy between
                                                                                                    ! the volume of the cluster and the the total volume of grains
    enddo

    !----------------------------------------------------------------------------------------------
    ! calculate the stress and penalty due to volume discrepancy
    vPen      = 0.0_pReal
    do i = 1,nGrain
      vPen(:,:,i) = -1.0_pReal/real(nGrain,pReal)*num%volDiscrMod*num%volDiscrPow/num%maxVolDiscr* &
                         sign((abs(vDiscrep)/num%maxVolDiscr)**(num%volDiscrPow - 1.0),vDiscrep)* &
                         gVol(i)*transpose(math_inv33(fDef(:,:,i)))
    enddo

  end subroutine volumePenalty


  !--------------------------------------------------------------------------------------------------
  !> @brief compute the correction factor accouted for surface evolution (area change) due to
  ! deformation
  !--------------------------------------------------------------------------------------------------
  function surfaceCorrection(avgF,ho,en)

    real(pReal), dimension(3)               :: surfaceCorrection

    real(pReal), dimension(3,3), intent(in) :: avgF                                                 !< average F
    integer,                     intent(in) :: &
      ho, &
      en
    real(pReal), dimension(3,3)             :: invC
    real(pReal), dimension(3)               :: nVect
    real(pReal)  :: detF
    integer :: i,j,iBase
    logical     ::  error

    call math_invert33(invC,detF,error,matmul(transpose(avgF),avgF))

    surfaceCorrection = 0.0_pReal
    do iBase = 1,3
      nVect = interfaceNormal([iBase,1,1,1],ho,en)
      do i = 1,3; do j = 1,3
        surfaceCorrection(iBase) = surfaceCorrection(iBase) + invC(i,j)*nVect(i)*nVect(j)           ! compute the component of (the inverse of) the stretch in the direction of the normal
      enddo; enddo
      surfaceCorrection(iBase) = sqrt(surfaceCorrection(iBase))*detF                                ! get the surface correction factor (area contraction/enlargement)
    enddo

  end function surfaceCorrection


  !-------------------------------------------------------------------------------------------------
  !> @brief compute the equivalent shear and bulk moduli from the elasticity tensor
  !-------------------------------------------------------------------------------------------------
  real(pReal) function equivalentMu(grainID,ce)

    integer, intent(in)    :: &
      grainID,&
      ce

    real(pReal), dimension(6,6) :: C


    C = phase_homogenizedC(material_phaseID(grainID,ce),material_phaseEntry(grainID,ce))
    equivalentMu = lattice_equivalent_mu(C,'voigt')

  end function equivalentMu


  !-------------------------------------------------------------------------------------------------
  !> @brief calculating the grain deformation gradient (the same with
  ! homogenization_RGC_partitionDeformation, but used only for perturbation scheme)
  !-------------------------------------------------------------------------------------------------
  subroutine grainDeformation(F, avgF, ho, en)

    real(pReal),   dimension(:,:,:), intent(out) :: F                                               !< partitioned F  per grain

    real(pReal),   dimension(:,:),   intent(in)  :: avgF                                            !< averaged F
    integer,                          intent(in)  :: &
      ho, &
      en

    real(pReal),   dimension(3) :: aVect,nVect
    integer,       dimension(4) :: intFace
    integer,       dimension(3) :: iGrain3
    integer :: iGrain,iFace,i,j

    !-----------------------------------------------------------------------------------------------
    ! compute the deformation gradient of individual grains due to relaxations

    associate (prm => param(ho))

    F = 0.0_pReal
    do iGrain = 1,product(prm%N_constituents)
      iGrain3 = grain1to3(iGrain,prm%N_constituents)
      do iFace = 1,6
        intFace = getInterface(iFace,iGrain3)
        aVect   = relaxationVector(intFace,ho,en)
        nVect   = interfaceNormal(intFace,ho,en)
        forall (i=1:3,j=1:3) &
          F(i,j,iGrain) = F(i,j,iGrain) + aVect(i)*nVect(j)                                         ! effective relaxations
      enddo
      F(1:3,1:3,iGrain) = F(1:3,1:3,iGrain) + avgF                                                  ! relaxed deformation gradient
    enddo

    end associate

  end subroutine grainDeformation

end function RGC_updateState


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
module subroutine RGC_results(ho,group)

  integer,          intent(in) :: ho
  character(len=*), intent(in) :: group

  integer :: o

  associate(stt => state(ho), dst => dependentState(ho), prm => param(ho))
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case('M')
        call results_writeDataset(dst%mismatch,group,trim(prm%output(o)), &
                                  'average mismatch tensor','1')
      case('Delta_V')
        call results_writeDataset(dst%volumeDiscrepancy,group,trim(prm%output(o)), &
                                  'volume discrepancy','m³')
      case('max_a_dot')
        call results_writeDataset(dst%relaxationrate_max,group,trim(prm%output(o)), &
                                  'maximum relaxation rate','m/s')
      case('avg_a_dot')
        call results_writeDataset(dst%relaxationrate_avg,group,trim(prm%output(o)), &
                                  'average relaxation rate','m/s')
    end select
  enddo outputsLoop
  end associate

end subroutine RGC_results


!--------------------------------------------------------------------------------------------------
!> @brief collect relaxation vectors of an interface
!--------------------------------------------------------------------------------------------------
pure function relaxationVector(intFace,ho,en)

  real(pReal), dimension (3)            :: relaxationVector

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
    relaxationVector = 0.0_pReal
  endif

  end associate

end function relaxationVector


!--------------------------------------------------------------------------------------------------
!> @brief identify the normal of an interface
!--------------------------------------------------------------------------------------------------
pure function interfaceNormal(intFace,ho,en)

  real(pReal), dimension(3)             :: interfaceNormal

  integer,     dimension(4), intent(in) :: intFace                                                  !< interface ID in 4D array (normal and position)
  integer,                   intent(in) :: &
    ho, &
    en

  integer :: nPos
  associate (dst => dependentState(ho))

!--------------------------------------------------------------------------------------------------
! get the normal of the interface, identified from the value of intFace(1)
  interfaceNormal = 0.0_pReal
  nPos = abs(intFace(1))                                                                            ! identify the position of the interface in global state array
  interfaceNormal(nPos) = real(intFace(1)/abs(intFace(1)),pReal)                                    ! get the normal vector w.r.t. cluster axis

  interfaceNormal = matmul(dst%orientation(1:3,1:3,en),interfaceNormal)                             ! map the normal vector into sample coordinate system (basis)

  end associate

end function interfaceNormal


!--------------------------------------------------------------------------------------------------
!> @brief collect six faces of a grain in 4D (normal and position)
!--------------------------------------------------------------------------------------------------
pure function getInterface(iFace,iGrain3)

  integer, dimension(4)             :: getInterface

  integer, dimension(3), intent(in) :: iGrain3                                                      !< grain ID in 3D array
  integer,               intent(in) :: iFace                                                        !< face index (1..6) mapped like (-e1,-e2,-e3,+e1,+e2,+e3) or iDir = (-1,-2,-3,1,2,3)

  integer :: iDir                                                                                   !< direction of interface normal

 iDir = (int(real(iFace-1,pReal)/2.0_pReal)+1)*(-1)**iFace
 getInterface(1) = iDir

!--------------------------------------------------------------------------------------------------
! identify the interface position by the direction of its normal
  getInterface(2:4) = iGrain3
  if (iDir < 0) getInterface(1-iDir) = getInterface(1-iDir)-1                                       ! to have a correlation with coordinate/position in real space

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
      endif

    case(2)
      if ((iFace4D(3) == 0) .or. (iFace4D(3) == nGDim(2))) then
        interface4to1 = 0
      else
        interface4to1 = iFace4D(4) + nGDim(3)*(iFace4D(2)-1) &
                      + nGDim(3)*nGDim(1)*(iFace4D(3)-1) &
                      + (nGDim(1)-1)*nGDim(2)*nGDim(3)                                              ! total # of interfaces normal || e1
      endif

    case(3)
      if ((iFace4D(4) == 0) .or. (iFace4D(4) == nGDim(3))) then
        interface4to1 = 0
      else
        interface4to1 = iFace4D(2) + nGDim(1)*(iFace4D(3)-1) &
                      + nGDim(1)*nGDim(2)*(iFace4D(4)-1) &
                      + (nGDim(1)-1)*nGDim(2)*nGDim(3) &                                            ! total # of interfaces normal || e1
                      + nGDim(1)*(nGDim(2)-1)*nGDim(3)                                              ! total # of interfaces normal || e2
      endif

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
    interface1to4(4) = mod(int(real(iFace1D-1,pReal)/real(nGDim(2),pReal)),nGDim(3))+1
    interface1to4(2) = int(real(iFace1D-1,pReal)/real(nGDim(2),pReal)/real(nGDim(3),pReal))+1
  elseif (iFace1D > nIntFace(1) .and. iFace1D <= (nIntFace(2) + nIntFace(1))) then                  ! interface with normal || e2
    interface1to4(1) = 2
    interface1to4(4) = mod((iFace1D-nIntFace(1)-1),nGDim(3))+1
    interface1to4(2) = mod(int(real(iFace1D-nIntFace(1)-1,pReal)/real(nGDim(3),pReal)),nGDim(1))+1
    interface1to4(3) = int(real(iFace1D-nIntFace(1)-1,pReal)/real(nGDim(3),pReal)/real(nGDim(1),pReal))+1
  elseif (iFace1D > nIntFace(2) + nIntFace(1) .and. iFace1D <= (nIntFace(3) + nIntFace(2) + nIntFace(1))) then ! interface with normal || e3
    interface1to4(1) = 3
    interface1to4(2) = mod((iFace1D-nIntFace(2)-nIntFace(1)-1),nGDim(1))+1
    interface1to4(3) = mod(int(real(iFace1D-nIntFace(2)-nIntFace(1)-1,pReal)/real(nGDim(1),pReal)),nGDim(2))+1
    interface1to4(4) = int(real(iFace1D-nIntFace(2)-nIntFace(1)-1,pReal)/real(nGDim(1),pReal)/real(nGDim(2),pReal))+1
  endif

end function interface1to4


end submodule RGC
