!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Denny Tjahjanto, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Relaxed grain cluster (RGC) homogenization scheme
!> Nconstituents is defined as p x q x r (cluster)
!--------------------------------------------------------------------------------------------------
module homogenization_RGC
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),              dimension(:,:),     allocatable,target, public :: &
   homogenization_RGC_sizePostResult
 character(len=64),          dimension(:,:),     allocatable,target, public :: &
   homogenization_RGC_output                                                                        ! name of each post result output

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 constitutivework_ID, &
                 penaltyenergy_ID, &
                 volumediscrepancy_ID, &
                 averagerelaxrate_ID,&
                 maximumrelaxrate_ID,&
                 ipcoords_ID,&
                 magnitudemismatch_ID,&
                 avgdefgrad_ID,&
                 avgfirstpiola_ID
 end enum
 
 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   integer(pInt), dimension(:), allocatable :: &
     Nconstituents
   real(pReal) :: &
     xiAlpha, &
     ciAlpha
   real(pReal), dimension(:), allocatable :: &
     dAlpha, &
     angles
   integer(kind(undefined_ID)),         dimension(:),   allocatable   :: &
     outputID                                                                                       !< ID of each post result output
 end type

 type, private :: tRGCstate
   real(pReal), pointer,     dimension(:) :: &
     work, &
     penaltyEnergy, &
     volumeDiscrepancy, &
     relaxationRate_avg, &
     relaxationRate_max
   real(pReal), pointer,     dimension(:,:) :: &
     relaxationVector, &
     mismatch
 end type tRGCstate

 type, private :: tRGCdependentState
   real(pReal), allocatable,     dimension(:,:,:) :: &
     orientation
 end type tRGCdependentState

 type(tparameters), dimension(:), allocatable, private :: param                                     !< containers of parameters (len Ninstance)
 type(tRGCstate),   dimension(:), allocatable, private :: state
 type(tRGCdependentState),   dimension(:), allocatable, private :: dependentState

 public :: &
   homogenization_RGC_init, &
   homogenization_RGC_partitionDeformation, &
   homogenization_RGC_averageStressAndItsTangent, &
   homogenization_RGC_updateState, &
   homogenization_RGC_postResults
 private :: &
   relaxationVector, &
   interfaceNormal, &
   getInterface, &
   grain1to3, &
   grain3to1, &
   interface4to1, &
   interface1to4

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all necessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_RGC_init()
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   pReal, &
   pInt 
 use debug, only: &
  debug_level, &
  debug_homogenization, &
  debug_levelBasic, &
  debug_levelExtensive
 use math, only: &
   math_EulerToR,&
   INRAD
 use mesh, only: &
   mesh_NcpElems,&
   mesh_NipsPerElem
 use IO
 use material
 use config

 implicit none
 integer(pInt) :: &
   NofMyHomog, &
   h, &
   outputSize, &
   instance, &
   sizeHState, nIntFaceTot
 integer(pInt) :: maxNinstance, i,j,e, of
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
 integer(kind(undefined_ID))                 :: &
   outputID                                                                                     !< ID of each post result output
 character(len=65536), dimension(:), allocatable :: outputs
 
 write(6,'(/,a)')   ' <<<+-  homogenization_'//HOMOGENIZATION_RGC_label//' init  -+>>>'
 write(6,'(/,a)')   ' Tjahjanto et al., International Journal of Material Forming, 2(1):939–942, 2009'
 write(6,'(a)')     ' https://doi.org/10.1007/s12289-009-0619-1'
 write(6,'(/,a)')   ' Tjahjanto et al., Modelling and Simulation in Materials Science and Engineering, 18:015006, 2010'
 write(6,'(a)')     ' https://doi.org/10.1088/0965-0393/18/1/015006'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(homogenization_type == HOMOGENIZATION_RGC_ID),pInt)
 if (maxNinstance == 0_pInt) return
 if (iand(debug_level(debug_HOMOGENIZATION),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 allocate(param(maxNinstance))                                                                      ! one container of parameters per instance
 allocate(state(maxNinstance))                                                                      ! one container per instance
 allocate(dependentState(maxNinstance))                                                             ! one container per instance

 allocate(homogenization_RGC_output(maxval(homogenization_Noutput),maxNinstance))
                                                              homogenization_RGC_output=''
 allocate(homogenization_RGC_sizePostResult(maxval(homogenization_Noutput),maxNinstance),&
                                                                                 source=0_pInt)

 do h = 1_pInt, size(homogenization_type)
   if (homogenization_type(h) /= HOMOGENIZATION_RGC_ID) cycle
   instance = homogenization_typeInstance(h)
   associate(prm => param(instance))

   prm%Nconstituents = config_homogenization(h)%getInts('clustersize',requiredShape=[3])
   if (homogenization_Ngrains(h) /= product(prm%Nconstituents)) &
     call IO_error(211_pInt,ext_msg='clustersize ('//HOMOGENIZATION_RGC_label//')')
   prm%xiAlpha = config_homogenization(h)%getFloat('scalingparameter')
   prm%ciAlpha = config_homogenization(h)%getFloat('overproportionality')
   prm%dAlpha  = config_homogenization(h)%getFloats('grainsize',requiredShape=[3])
   prm%angles  = config_homogenization(h)%getFloats('clusterorientation',requiredShape=[3])

   outputs = config_homogenization(h)%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(prm%outputID(0))

   do i=1_pInt, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
       case('constitutivework')
         outputID = constitutivework_ID
         outputSize = 1_pInt
       case('penaltyenergy')
         outputID = penaltyenergy_ID
         outputSize = 1_pInt
       case('volumediscrepancy')
         outputID = volumediscrepancy_ID
         outputSize = 1_pInt
       case('averagerelaxrate')
         outputID = averagerelaxrate_ID
         outputSize = 1_pInt
       case('maximumrelaxrate')
         outputID = maximumrelaxrate_ID
         outputSize = 1_pInt
       case('magnitudemismatch')
         outputID = magnitudemismatch_ID
         outputSize = 3_pInt
       case default
     end select
     if (outputID /= undefined_ID) then
       homogenization_RGC_output(i,instance) = outputs(i)
       homogenization_RGC_sizePostResult(i,instance) = outputSize
       prm%outputID = [prm%outputID , outputID]
     endif
   enddo

   if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt) then
     write(6,'(a15,1x,i4,/)')     'instance:  ', instance
     write(6,'(a25,3(1x,i8))')    'cluster size:         ',(prm%Nconstituents(j),j=1_pInt,3_pInt)
     write(6,'(a25,1x,e10.3)')    'scaling parameter:    ', prm%xiAlpha
     write(6,'(a25,1x,e10.3)')    'over-proportionality: ', prm%ciAlpha
     write(6,'(a25,3(1x,e10.3))') 'grain size:           ',(prm%dAlpha(j),j=1_pInt,3_pInt)
     write(6,'(a25,3(1x,e10.3))') 'cluster orientation:  ',(prm%angles(j),j=1_pInt,3_pInt)
   endif

   NofMyHomog = count(material_homog == h)
   nIntFaceTot = 3_pInt*(  (prm%Nconstituents(1)-1_pInt)*prm%Nconstituents(2)*prm%Nconstituents(3) &
                         + prm%Nconstituents(1)*(prm%Nconstituents(2)-1_pInt)*prm%Nconstituents(3) &
                         + prm%Nconstituents(1)*prm%Nconstituents(2)*(prm%Nconstituents(3)-1_pInt))
   sizeHState = nIntFaceTot &
       + 8_pInt   ! (1) Average constitutive work, (2-4) Overall mismatch, (5) Average penalty energy, 
                  ! (6) Volume discrepancy, (7) Avg relaxation rate component, (8) Max relaxation rate component
                  
   homogState(h)%sizeState = sizeHState
   homogState(h)%sizePostResults = sum(homogenization_RGC_sizePostResult(:,instance))
   allocate(homogState(h)%state0   (sizeHState,NofMyHomog), source=0.0_pReal)
   allocate(homogState(h)%subState0(sizeHState,NofMyHomog), source=0.0_pReal)
   allocate(homogState(h)%state    (sizeHState,NofMyHomog), source=0.0_pReal)

   state(instance)%relaxationVector   => homogState(h)%state(1:nIntFaceTot,:)
   state(instance)%work               => homogState(h)%state(nIntFaceTot+1,:)
   state(instance)%mismatch           => homogState(h)%state(nIntFaceTot+2:nIntFaceTot+4,:)
   state(instance)%penaltyEnergy      => homogState(h)%state(nIntFaceTot+5,:)
   state(instance)%volumeDiscrepancy  => homogState(h)%state(nIntFaceTot+6,:)
   state(instance)%relaxationRate_avg => homogState(h)%state(nIntFaceTot+7,:)
   state(instance)%relaxationRate_max => homogState(h)%state(nIntFaceTot+8,:)

   allocate(dependentState(instance)%orientation(3,3,NofMyHomog))


!--------------------------------------------------------------------------------------------------
! * assigning cluster orientations
   elementLooping: do e = 1_pInt,mesh_NcpElems
     if (homogenization_typeInstance(material_homogenizationAt(e)) == instance .and. NofMyHomog > 0_pInt) then
       do i = 1_pInt,mesh_NipsPerElem
         of = mappingHomogenization(1,i,e)
         dependentState(instance)%orientation(1:3,1:3,of) = math_EulerToR(prm%angles*inRad)
       enddo
     endif
   enddo elementLooping
   end associate
 enddo  
 
end subroutine homogenization_RGC_init


!--------------------------------------------------------------------------------------------------
!> @brief partitions the deformation gradient onto the constituents
!--------------------------------------------------------------------------------------------------
subroutine homogenization_RGC_partitionDeformation(F,avgF,instance,of)
 use debug, only: &
   debug_level, &
   debug_homogenization, &
   debug_levelExtensive
 use material, only: &
   homogenization_maxNgrains

 implicit none
 real(pReal),   dimension (3,3,homogenization_maxNgrains), intent(out) :: F                         !< partioned F  per grain
 real(pReal),   dimension (3,3),                           intent(in)  :: avgF                      !< averaged F
 integer(pInt),                                            intent(in)  :: &
   instance, &
   of                                                                                               !< element number
 real(pReal),   dimension (3) :: aVect,nVect
 integer(pInt), dimension (4) :: intFace
 integer(pInt), dimension (3) :: iGrain3
 integer(pInt) ::  iGrain,iFace,i,j

!--------------------------------------------------------------------------------------------------
! compute the deformation gradient of individual grains due to relaxations
 associate(prm => param(instance))
 F = 0.0_pReal
 do iGrain = 1_pInt,product(prm%Nconstituents)
   iGrain3 = grain1to3(iGrain,prm%Nconstituents)
   do iFace = 1_pInt,6_pInt
     intFace = getInterface(iFace,iGrain3)                                                          ! identifying 6 interfaces of each grain
     aVect = relaxationVector(intFace,instance,of)                                                  ! get the relaxation vectors for each interface from global relaxation vector array
     nVect = interfaceNormal(intFace,instance,of)
     forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) &
       F(i,j,iGrain) = F(i,j,iGrain) + aVect(i)*nVect(j)                                            ! calculating deformation relaxations due to interface relaxation
   enddo
   F(1:3,1:3,iGrain) = F(1:3,1:3,iGrain) + avgF                                                     ! resulting relaxed deformation gradient

!--------------------------------------------------------------------------------------------------
! debugging the grain deformation gradients
   if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt) then
     !$OMP CRITICAL (write2out)
     write(6,'(1x,a32,1x,i3)')'Deformation gradient of grain: ',iGrain
     do i = 1_pInt,3_pInt
       write(6,'(1x,3(e15.8,1x))')(F(i,j,iGrain), j = 1_pInt,3_pInt)
     enddo
     write(6,*)' '
     flush(6)
     !$OMP END CRITICAL (write2out)
   endif

 enddo
 end associate

end subroutine homogenization_RGC_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief update the internal state of the homogenization scheme and tell whether "done" and 
! "happy" with result
!--------------------------------------------------------------------------------------------------
function homogenization_RGC_updateState(P,F,F0,avgF,dt,dPdF,ip,el)
 use prec, only: &
   dEq0
 use debug, only: &
   debug_level, &
   debug_homogenization,&
   debug_levelExtensive, &
   debug_e, &
   debug_i
 use math, only: &
   math_invert
 use material, only: &
   material_homogenizationAt, &
   homogenization_maxNgrains, &
   homogenization_typeInstance, &
   homogState, &
   mappingHomogenization, &   
   homogenization_Ngrains
 use numerics, only: &
   absTol_RGC, &
   relTol_RGC, &
   absMax_RGC, &
   relMax_RGC, &
   pPert_RGC, &
   maxdRelax_RGC, &
   viscPower_RGC, &
   viscModus_RGC, &
   refRelaxRate_RGC

 implicit none

 real(pReal), dimension (3,3,homogenization_maxNgrains),     intent(in)    :: & 
   P,&                                                                                              !< array of P
   F,&                                                                                              !< array of F
   F0                                                                                               !< array of initial F
 real(pReal), dimension (3,3,3,3,homogenization_maxNgrains), intent(in) :: dPdF                     !< array of current grain stiffness
 real(pReal), dimension (3,3),                               intent(in) :: avgF                     !< average F
 real(pReal),                                                intent(in) :: dt                       !< time increment
 integer(pInt),                                              intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number

 logical, dimension(2)        :: homogenization_RGC_updateState
 integer(pInt), dimension (4) :: intFaceN,intFaceP,faceID
 integer(pInt), dimension (3) :: nGDim,iGr3N,iGr3P,stresLoc
 integer(pInt), dimension (2) :: residLoc
 integer(pInt) instance,iNum,i,j,nIntFaceTot,iGrN,iGrP,iMun,iFace,k,l,ipert,iGrain,nGrain, of
 real(pReal), dimension (3,3,homogenization_maxNgrains) :: R,pF,pR,D,pD
 real(pReal), dimension (3,homogenization_maxNgrains)   :: NN,pNN
 real(pReal), dimension (3)                             :: normP,normN,mornP,mornN
 real(pReal) :: residMax,stresMax,volDiscrep
 logical error

 real(pReal), dimension(:,:), allocatable :: tract,jmatrix,jnverse,smatrix,pmatrix,rmatrix
 real(pReal), dimension(:), allocatable   :: resid,relax,p_relax,p_resid,drelax
 
 zeroTimeStep: if(dEq0(dt)) then
   homogenization_RGC_updateState = .true.                                                          ! pretend everything is fine and return
   return                                                                    
 endif zeroTimeStep

!--------------------------------------------------------------------------------------------------
! get the dimension of the cluster (grains and interfaces)
 instance  = homogenization_typeInstance(material_homogenizationAt(el))
 of = mappingHomogenization(1,ip,el)
 nGDim  = param(instance)%Nconstituents
 nGrain = homogenization_Ngrains(material_homogenizationAt(el))
 nIntFaceTot = (nGDim(1)-1_pInt)*nGDim(2)*nGDim(3) &
             + nGDim(1)*(nGDim(2)-1_pInt)*nGDim(3) &
             + nGDim(1)*nGDim(2)*(nGDim(3)-1_pInt)

!--------------------------------------------------------------------------------------------------
! allocate the size of the global relaxation arrays/jacobian matrices depending on the size of the cluster
 allocate(resid(3_pInt*nIntFaceTot),   source=0.0_pReal)
 allocate(tract(nIntFaceTot,3),        source=0.0_pReal)
 relax  = homogState(mappingHomogenization(2,ip,el))%state (1:3_pInt*nIntFaceTot,of)
 drelax = relax &
        - homogState(mappingHomogenization(2,ip,el))%state0(1:3_pInt*nIntFaceTot,of)
!--------------------------------------------------------------------------------------------------
! debugging the obtained state
 if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Obtained state: '
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,2(e15.8,1x))')homogState(mappingHomogenization(2,ip,el))%state(i,of)
   enddo
   write(6,*)' '
   !$OMP END CRITICAL (write2out)
 endif

!--------------------------------------------------------------------------------------------------
! computing interface mismatch and stress penalty tensor for all interfaces of all grains
 call stressPenalty(R,NN,avgF,F,ip,el,instance)

!--------------------------------------------------------------------------------------------------
! calculating volume discrepancy and stress penalty related to overall volume discrepancy 
 call volumePenalty(D,volDiscrep,F,avgF,ip,el)

!--------------------------------------------------------------------------------------------------
! debugging the mismatch, stress and penalties of grains
 if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   do iGrain = 1_pInt,nGrain
     write(6,'(1x,a30,1x,i3,1x,a4,3(1x,e15.8))')'Mismatch magnitude of grain(',iGrain,') :',&
       NN(1,iGrain),NN(2,iGrain),NN(3,iGrain)
     write(6,'(/,1x,a30,1x,i3)')'Stress and penalties of grain: ',iGrain
     do i = 1_pInt,3_pInt
       write(6,'(1x,3(e15.8,1x),1x,3(e15.8,1x),1x,3(e15.8,1x))')(P(i,j,iGrain), j = 1_pInt,3_pInt), &
                                                                (R(i,j,iGrain), j = 1_pInt,3_pInt), &
                                                                (D(i,j,iGrain), j = 1_pInt,3_pInt)
     enddo
     write(6,*)' '
   enddo
   !$OMP END CRITICAL (write2out)
 endif

!------------------------------------------------------------------------------------------------
! computing the residual stress from the balance of traction at all (interior) interfaces
 do iNum = 1_pInt,nIntFaceTot
   faceID = interface1to4(iNum,param(instance)%Nconstituents)                                         ! identifying the interface ID in local coordinate system (4-dimensional index)

!--------------------------------------------------------------------------------------------------
! identify the left/bottom/back grain (-|N)
   iGr3N = faceID(2:4)                                                                              ! identifying the grain ID in local coordinate system (3-dimensional index)
   iGrN = grain3to1(iGr3N,param(instance)%Nconstituents)                                              ! translate the local grain ID into global coordinate system (1-dimensional index)
   intFaceN = getInterface(2_pInt*faceID(1),iGr3N)
   normN = interfaceNormal(intFaceN,instance,of)
   
!--------------------------------------------------------------------------------------------------
! identify the right/up/front grain (+|P)
   iGr3P = iGr3N
   iGr3P(faceID(1)) = iGr3N(faceID(1))+1_pInt                                                       ! identifying the grain ID in local coordinate system (3-dimensional index)
   iGrP = grain3to1(iGr3P,param(instance)%Nconstituents)                                              ! translate the local grain ID into global coordinate system (1-dimensional index)
   intFaceP = getInterface(2_pInt*faceID(1)-1_pInt,iGr3P)
   normP = interfaceNormal(intFaceP,instance,of)

!--------------------------------------------------------------------------------------------------
! compute the residual of traction at the interface (in local system, 4-dimensional index)
   do i = 1_pInt,3_pInt
     tract(iNum,i) = sign(viscModus_RGC*(abs(drelax(i+3*(iNum-1_pInt)))/(refRelaxRate_RGC*dt))**viscPower_RGC, &
                          drelax(i+3*(iNum-1_pInt)))                                                ! contribution from the relaxation viscosity
     do j = 1_pInt,3_pInt
       tract(iNum,i) = tract(iNum,i) + (P(i,j,iGrP) + R(i,j,iGrP) + D(i,j,iGrP))*normP(j) &         ! contribution from material stress P, mismatch penalty R, and volume penalty D projected into the interface
                                     + (P(i,j,iGrN) + R(i,j,iGrN) + D(i,j,iGrN))*normN(j)
       resid(i+3_pInt*(iNum-1_pInt)) = tract(iNum,i)                                                ! translate the local residual into global 1-dimensional residual array
     enddo
   enddo
   
!--------------------------------------------------------------------------------------------------
! debugging the residual stress
   if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt) then
     !$OMP CRITICAL (write2out)
     write(6,'(1x,a30,1x,i3)')'Traction at interface: ',iNum
     write(6,'(1x,3(e15.8,1x))')(tract(iNum,j), j = 1_pInt,3_pInt)
     write(6,*)' '
     !$OMP END CRITICAL (write2out)
   endif
 enddo

!--------------------------------------------------------------------------------------------------
! convergence check for stress residual
 stresMax = maxval(abs(P))                                                                          ! get the maximum of first Piola-Kirchhoff (material) stress
 stresLoc = int(maxloc(abs(P)),pInt)                                                                ! get the location of the maximum stress
 residMax = maxval(abs(tract))                                                                      ! get the maximum of the residual
 residLoc = int(maxloc(abs(tract)),pInt)                                                            ! get the position of the maximum residual

!--------------------------------------------------------------------------------------------------
!  Debugging the convergent criteria
 if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt &
     .and. debug_e == el .and. debug_i == ip) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a)')' '
   write(6,'(1x,a,1x,i2,1x,i4)')'RGC residual check ...',ip,el
   write(6,'(1x,a15,1x,e15.8,1x,a7,i3,1x,a12,i2,i2)')'Max stress: ',stresMax, &
              '@ grain',stresLoc(3),'in component',stresLoc(1),stresLoc(2)
   write(6,'(1x,a15,1x,e15.8,1x,a7,i3,1x,a12,i2)')'Max residual: ',residMax, &
              '@ iface',residLoc(1),'in direction',residLoc(2)
   flush(6)
   !$OMP END CRITICAL (write2out)
 endif
 
 homogenization_RGC_updateState = .false.
 
!--------------------------------------------------------------------------------------------------
!  If convergence reached => done and happy
 if (residMax < relTol_RGC*stresMax .or. residMax < absTol_RGC) then 
   homogenization_RGC_updateState = .true.
   
    if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt &
        .and. debug_e == el .and. debug_i == ip) then 
     !$OMP CRITICAL (write2out)
     write(6,'(1x,a55,/)')'... done and happy'
     flush(6)
     !$OMP END CRITICAL (write2out)
   endif

!--------------------------------------------------------------------------------------------------
! compute/update the state for postResult, i.e., all energy densities computed by time-integration
   do iGrain = 1_pInt,homogenization_Ngrains(material_homogenizationAt(el))                             ! time-integration loop for work and energy
     do i = 1_pInt,3_pInt;do j = 1_pInt,3_pInt
       state(instance)%work(of)          = state(instance)%work(of) &
                                         + P(i,j,iGrain)*(F(i,j,iGrain) - F0(i,j,iGrain))/real(nGrain,pReal)
       state(instance)%penaltyEnergy(of) = state(instance)%penaltyEnergy(of) &
                                         + R(i,j,iGrain)*(F(i,j,iGrain) - F0(i,j,iGrain))/real(nGrain,pReal)
     enddo; enddo
   enddo

   state(instance)%mismatch(1:3,of)       = sum(NN,2)/real(nGrain,pReal)                            ! the overall mismatch of all interface normals
   state(instance)%volumeDiscrepancy(of)  = volDiscrep
   state(instance)%relaxationRate_avg(of) = sum(abs(drelax))/dt/real(3_pInt*nIntFaceTot,pReal)
   state(instance)%relaxationRate_max(of) = maxval(abs(drelax))/dt

   if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt &
        .and. debug_e == el .and. debug_i == ip) then
     !$OMP CRITICAL (write2out)
     write(6,'(1x,a30,1x,e15.8)')   'Constitutive work: ',state(instance)%work(of)
     write(6,'(1x,a30,3(1x,e15.8))')'Magnitude mismatch: ',state(instance)%mismatch(1,of), &
                                                           state(instance)%mismatch(2,of), &
                                                           state(instance)%mismatch(3,of)
     write(6,'(1x,a30,1x,e15.8)')   'Penalty energy: ',          state(instance)%penaltyEnergy(of)
     write(6,'(1x,a30,1x,e15.8,/)') 'Volume discrepancy: ',      state(instance)%volumeDiscrepancy(of)
     write(6,'(1x,a30,1x,e15.8)')   'Maximum relaxation rate: ', state(instance)%relaxationRate_max(of)
     write(6,'(1x,a30,1x,e15.8,/)') 'Average relaxation rate: ', state(instance)%relaxationRate_avg(of)
     flush(6)
     !$OMP END CRITICAL (write2out)
   endif
   
   return

!--------------------------------------------------------------------------------------------------
! if residual blows-up => done but unhappy
 elseif (residMax > relMax_RGC*stresMax .or. residMax > absMax_RGC) then                            ! try to restart when residual blows up exceeding maximum bound
   homogenization_RGC_updateState = [.true.,.false.]                                                ! with direct cut-back

   if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt &
       .and. debug_e == el .and. debug_i == ip) then
     !$OMP CRITICAL (write2out)
     write(6,'(1x,a55,/)')'... broken'
     flush(6)
     !$OMP END CRITICAL (write2out)
   endif
   
   return 
 else                                                                                               ! proceed with computing the Jacobian and state update
   if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt &
     .and. debug_e == el .and. debug_i == ip) then
     !$OMP CRITICAL (write2out)
     write(6,'(1x,a55,/)')'... not yet done'
     flush(6)
     !$OMP END CRITICAL (write2out)
   endif
   
 endif

!---------------------------------------------------------------------------------------------------
! construct the global Jacobian matrix for updating the global relaxation vector array when 
! convergence is not yet reached ...

!--------------------------------------------------------------------------------------------------
! ... of the constitutive stress tangent, assembled from dPdF or material constitutive model "smatrix"
 allocate(smatrix(3*nIntFaceTot,3*nIntFaceTot), source=0.0_pReal)
 do iNum = 1_pInt,nIntFaceTot
   faceID = interface1to4(iNum,param(instance)%Nconstituents)                                       ! assembling of local dPdF into global Jacobian matrix
   
!--------------------------------------------------------------------------------------------------
! identify the left/bottom/back grain (-|N)
   iGr3N = faceID(2:4)                                                                              ! identifying the grain ID in local coordinate sytem
   iGrN = grain3to1(iGr3N,param(instance)%Nconstituents)                                            ! translate into global grain ID
   intFaceN = getInterface(2_pInt*faceID(1),iGr3N)                                                  ! identifying the connecting interface in local coordinate system
   normN = interfaceNormal(intFaceN,instance,of)
   do iFace = 1_pInt,6_pInt
     intFaceN = getInterface(iFace,iGr3N)                                                           ! identifying all interfaces that influence relaxation of the above interface
     mornN = interfaceNormal(intFaceN,instance,of)
     iMun = interface4to1(intFaceN,param(instance)%Nconstituents)                                     ! translate the interfaces ID into local 4-dimensional index
     if (iMun > 0) then                                                                             ! get the corresponding tangent
       do i=1_pInt,3_pInt; do j=1_pInt,3_pInt; do k=1_pInt,3_pInt; do l=1_pInt,3_pInt
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
   iGr3P(faceID(1)) = iGr3N(faceID(1))+1_pInt                                                       ! identifying the grain ID in local coordinate sytem
   iGrP = grain3to1(iGr3P,param(instance)%Nconstituents)                                            ! translate into global grain ID
   intFaceP = getInterface(2_pInt*faceID(1)-1_pInt,iGr3P)                                           ! identifying the connecting interface in local coordinate system
   normP = interfaceNormal(intFaceP,instance,of)
   do iFace = 1_pInt,6_pInt
     intFaceP = getInterface(iFace,iGr3P)                                                           ! identifying all interfaces that influence relaxation of the above interface
     mornP = interfaceNormal(intFaceP,instance,of)
     iMun = interface4to1(intFaceP,param(instance)%Nconstituents)                                   ! translate the interfaces ID into local 4-dimensional index
     if (iMun > 0_pInt) then                                                                        ! get the corresponding tangent
       do i=1_pInt,3_pInt; do j=1_pInt,3_pInt; do k=1_pInt,3_pInt; do l=1_pInt,3_pInt
         smatrix(3*(iNum-1)+i,3*(iMun-1)+j) = smatrix(3*(iNum-1)+i,3*(iMun-1)+j) &
                                            + dPdF(i,k,j,l,iGrP)*normP(k)*mornP(l)
       enddo;enddo;enddo;enddo
     endif
   enddo
 enddo
 
!--------------------------------------------------------------------------------------------------
! debugging the global Jacobian matrix of stress tangent
 if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Jacobian matrix of stress'
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,100(e11.4,1x))')(smatrix(i,j), j = 1_pInt,3_pInt*nIntFaceTot)
   enddo
   write(6,*)' '
   flush(6)
   !$OMP END CRITICAL (write2out)
 endif
 
!--------------------------------------------------------------------------------------------------
! ... of the stress penalty tangent (mismatch penalty and volume penalty, computed using numerical 
! perturbation method) "pmatrix"
 allocate(pmatrix(3*nIntFaceTot,3*nIntFaceTot), source=0.0_pReal)
 allocate(p_relax(3*nIntFaceTot),               source=0.0_pReal)
 allocate(p_resid(3*nIntFaceTot),               source=0.0_pReal)

 do ipert = 1_pInt,3_pInt*nIntFaceTot
   p_relax = relax
   p_relax(ipert) = relax(ipert) + pPert_RGC                                                        ! perturb the relaxation vector
   homogState(mappingHomogenization(2,ip,el))%state(1:3*nIntFaceTot,of) = p_relax
   call grainDeformation(pF,avgF,instance,of)                                                       ! rain deformation from perturbed state
   call stressPenalty(pR,pNN,avgF,pF,ip,el,instance)                                                ! stress penalty due to interface mismatch from perturbed state
   call volumePenalty(pD,volDiscrep,pF,avgF,ip,el)                                                  ! stress penalty due to volume discrepancy from perturbed state

!--------------------------------------------------------------------------------------------------
! computing the global stress residual array from the perturbed state
   p_resid = 0.0_pReal
   do iNum = 1_pInt,nIntFaceTot
     faceID = interface1to4(iNum,param(instance)%Nconstituents)                                     ! identifying the interface ID in local coordinate system (4-dimensional index)

!--------------------------------------------------------------------------------------------------
! identify the left/bottom/back grain (-|N)
     iGr3N = faceID(2:4)                                                                            ! identify the grain ID in local coordinate system (3-dimensional index)
     iGrN = grain3to1(iGr3N,param(instance)%Nconstituents)                                          ! translate the local grain ID into global coordinate system (1-dimensional index)
     intFaceN = getInterface(2_pInt*faceID(1),iGr3N)                                                ! identify the interface ID of the grain
     normN = interfaceNormal(intFaceN,instance,of)
 
!--------------------------------------------------------------------------------------------------
! identify the right/up/front grain (+|P)
     iGr3P = iGr3N
     iGr3P(faceID(1)) = iGr3N(faceID(1))+1_pInt                                                     ! identify the grain ID in local coordinate system (3-dimensional index)
     iGrP = grain3to1(iGr3P,param(instance)%Nconstituents)                                          ! translate the local grain ID into global coordinate system (1-dimensional index)
     intFaceP = getInterface(2_pInt*faceID(1)-1_pInt,iGr3P)                                         ! identify the interface ID of the grain
     normP = interfaceNormal(intFaceP,instance,of)
 
!--------------------------------------------------------------------------------------------------
! compute the residual stress (contribution of mismatch and volume penalties) from perturbed state 
! at all interfaces
     do i = 1_pInt,3_pInt; do j = 1_pInt,3_pInt
       p_resid(i+3*(iNum-1)) = p_resid(i+3*(iNum-1)) + (pR(i,j,iGrP) - R(i,j,iGrP))*normP(j) &
                                                     + (pR(i,j,iGrN) - R(i,j,iGrN))*normN(j) &
                                                     + (pD(i,j,iGrP) - D(i,j,iGrP))*normP(j) &
                                                     + (pD(i,j,iGrN) - D(i,j,iGrN))*normN(j)
     enddo; enddo
   enddo
   pmatrix(:,ipert) = p_resid/pPert_RGC
 enddo
 
!--------------------------------------------------------------------------------------------------
! debugging the global Jacobian matrix of penalty tangent
 if (iand(debug_level(debug_homogenization), debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Jacobian matrix of penalty'
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,100(e11.4,1x))')(pmatrix(i,j), j = 1_pInt,3_pInt*nIntFaceTot)
   enddo
   write(6,*)' '
   flush(6)
   !$OMP END CRITICAL (write2out)
 endif
 
!--------------------------------------------------------------------------------------------------
! ... of the numerical viscosity traction "rmatrix"
 allocate(rmatrix(3*nIntFaceTot,3*nIntFaceTot),source=0.0_pReal)
 forall (i=1_pInt:3_pInt*nIntFaceTot) &
   rmatrix(i,i) = viscModus_RGC*viscPower_RGC/(refRelaxRate_RGC*dt)* &                              ! tangent due to numerical viscosity traction appears
                  (abs(drelax(i))/(refRelaxRate_RGC*dt))**(viscPower_RGC - 1.0_pReal)               ! only in the main diagonal term 
                                                                           
                                                                           

!--------------------------------------------------------------------------------------------------
! debugging the global Jacobian matrix of numerical viscosity tangent
 if (iand(debug_level(debug_homogenization), debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Jacobian matrix of penalty'
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,100(e11.4,1x))')(rmatrix(i,j), j = 1_pInt,3_pInt*nIntFaceTot)
   enddo
   write(6,*)' '
   flush(6)
   !$OMP END CRITICAL (write2out)
 endif

!--------------------------------------------------------------------------------------------------
! The overall Jacobian matrix summarizing contributions of smatrix, pmatrix, rmatrix
 allocate(jmatrix(3*nIntFaceTot,3*nIntFaceTot)); jmatrix = smatrix + pmatrix + rmatrix
 
 if (iand(debug_level(debug_homogenization), debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Jacobian matrix (total)'
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,100(e11.4,1x))')(jmatrix(i,j), j = 1_pInt,3_pInt*nIntFaceTot)
   enddo
   write(6,*)' '
   flush(6)
   !$OMP END CRITICAL (write2out)
 endif

!--------------------------------------------------------------------------------------------------
! computing the update of the state variable (relaxation vectors) using the Jacobian matrix
 allocate(jnverse(3_pInt*nIntFaceTot,3_pInt*nIntFaceTot),source=0.0_pReal)
 call math_invert(size(jmatrix,1),jmatrix,jnverse,error)                                            ! Compute the inverse of the overall Jacobian matrix
 
!--------------------------------------------------------------------------------------------------
! debugging the inverse Jacobian matrix
 if (iand(debug_level(debug_homogenization), debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Jacobian inverse'
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,100(e11.4,1x))')(jnverse(i,j), j = 1_pInt,3_pInt*nIntFaceTot)
   enddo
   write(6,*)' '
   flush(6)
   !$OMP END CRITICAL (write2out)
 endif

!--------------------------------------------------------------------------------------------------
! calculate the state update (global relaxation vectors) for the next Newton-Raphson iteration
 drelax = 0.0_pReal
 do i = 1_pInt,3_pInt*nIntFaceTot;do j = 1_pInt,3_pInt*nIntFaceTot
   drelax(i) = drelax(i) - jnverse(i,j)*resid(j)                                                    ! Calculate the correction for the state variable
 enddo; enddo
 relax = relax + drelax                                                                             ! Updateing the state variable for the next iteration
 homogState(mappingHomogenization(2,ip,el))%state(1:3*nIntFaceTot,of) = relax
 if (any(abs(drelax) > maxdRelax_RGC)) then                                                         ! Forcing cutback when the incremental change of relaxation vector becomes too large
   homogenization_RGC_updateState = [.true.,.false.]
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a,1x,i3,1x,a,1x,i3,1x,a)')'RGC_updateState: ip',ip,'| el',el,'enforces cutback'
   write(6,'(1x,a,1x,e15.8)')'due to large relaxation change =',maxval(abs(drelax))
   flush(6)
   !$OMP END CRITICAL (write2out)
 endif

!--------------------------------------------------------------------------------------------------
! debugging the return state
 if (iand(debug_homogenization, debug_levelExtensive) > 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Returned state: '
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,2(e15.8,1x))')homogState(mappingHomogenization(2,ip,el))%state(i,of)
   enddo
   write(6,*)' '
   flush(6)
   !$OMP END CRITICAL (write2out)
 endif

 contains
 !--------------------------------------------------------------------------------------------------
 !> @brief calculate stress-like penalty due to deformation mismatch
 !--------------------------------------------------------------------------------------------------
 subroutine stressPenalty(rPen,nMis,avgF,fDef,ip,el,instance)
  use math, only: &
    math_civita
  use numerics, only: &
    xSmoo_RGC
  
  implicit none
  real(pReal),   dimension (3,3,homogenization_maxNgrains), intent(out) :: rPen                      !< stress-like penalty
  real(pReal),   dimension (3,homogenization_maxNgrains),   intent(out) :: nMis                      !< total amount of mismatch
  real(pReal),   dimension (3,3,homogenization_maxNgrains), intent(in)  :: fDef                      !< deformation gradients
  real(pReal),   dimension (3,3),                           intent(in)  :: avgF                      !< initial effective stretch tensor
  integer(pInt),                                            intent(in)  :: ip,el,instance
  integer(pInt), dimension (4)   :: intFace
  integer(pInt), dimension (3)   :: iGrain3,iGNghb3,nGDim
  real(pReal),   dimension (3,3) :: gDef,nDef
  real(pReal),   dimension (3)   :: nVect,surfCorr
  real(pReal),   dimension (2)   :: Gmoduli
  integer(pInt) :: iGrain,iGNghb,iFace,i,j,k,l,of
  real(pReal)   :: muGrain,muGNghb,nDefNorm,bgGrain,bgGNghb
  real(pReal),                                               parameter  :: nDefToler = 1.0e-10_pReal
  logical :: debugActive

  debugActive = iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt &
                .and. debug_e == el .and. debug_i == ip
 
  nGDim = param(instance)%Nconstituents
  rPen = 0.0_pReal
  nMis = 0.0_pReal
 
 !--------------------------------------------------------------------------------------------------
 ! get the correction factor the modulus of penalty stress representing the evolution of area of 
 ! the interfaces due to deformations

  of = mappingHomogenization(1,ip,el)
  surfCorr = surfaceCorrection(avgF,instance,of)
  associate(prm => param(instance))

  if (debugActive) then
    write(6,'(1x,a20,2(1x,i3))')'Correction factor: ',ip,el
    write(6,*) surfCorr
  endif
 
 !--------------------------------------------------------------------------------------------------
 ! computing the mismatch and penalty stress tensor of all grains 
  do iGrain = 1_pInt,product(prm%Nconstituents)
    Gmoduli = equivalentModuli(iGrain,ip,el)
    muGrain = Gmoduli(1)                                                                            ! collecting the equivalent shear modulus of grain
    bgGrain = Gmoduli(2)                                                                            ! and the lengthh of Burgers vector
    iGrain3 = grain1to3(iGrain,prm%Nconstituents)                                                   ! get the grain ID in local 3-dimensional index (x,y,z)-position
 
 !* Looping over all six interfaces of each grain
    do iFace = 1_pInt,6_pInt
      intFace = getInterface(iFace,iGrain3)                                                         ! get the 4-dimensional index of the interface in local numbering system of the grain
      nVect = interfaceNormal(intFace,instance,of)
      iGNghb3 = iGrain3                                                                             ! identify the neighboring grain across the interface
      iGNghb3(abs(intFace(1))) = iGNghb3(abs(intFace(1))) &
                               + int(real(intFace(1),pReal)/real(abs(intFace(1)),pReal),pInt)
      where(iGNghb3 < 1)    iGNghb3 = nGDim
      where(iGNghb3 >nGDim) iGNghb3 = 1_pInt
      iGNghb  = grain3to1(iGNghb3,prm%Nconstituents)                                                ! get the ID of the neighboring grain
      Gmoduli = equivalentModuli(iGNghb,ip,el)                                                      ! collect the shear modulus and Burgers vector of the neighbor
      muGNghb = Gmoduli(1)
      bgGNghb = Gmoduli(2)
      gDef = 0.5_pReal*(fDef(1:3,1:3,iGNghb) - fDef(1:3,1:3,iGrain))                                ! difference/jump in deformation gradeint across the neighbor
 
 !--------------------------------------------------------------------------------------------------
 ! compute the mismatch tensor of all interfaces
      nDefNorm = 0.0_pReal
      nDef = 0.0_pReal
      do i = 1_pInt,3_pInt; do j = 1_pInt,3_pInt
        do k = 1_pInt,3_pInt; do l = 1_pInt,3_pInt
          nDef(i,j) = nDef(i,j) - nVect(k)*gDef(i,l)*math_civita(j,k,l)                             ! compute the interface mismatch tensor from the jump of deformation gradient
        enddo; enddo
        nDefNorm = nDefNorm + nDef(i,j)**2.0_pReal                                                  ! compute the norm of the mismatch tensor
      enddo; enddo
      nDefNorm = max(nDefToler,sqrt(nDefNorm))                                                      ! approximation to zero mismatch if mismatch is zero (singularity)
      nMis(abs(intFace(1)),iGrain) = nMis(abs(intFace(1)),iGrain) + nDefNorm                        ! total amount of mismatch experienced by the grain (at all six interfaces)
 
      if (debugActive) then
        write(6,'(1x,a20,i2,1x,a20,1x,i3)')'Mismatch to face: ',intFace(1),'neighbor grain: ',iGNghb
        write(6,*) transpose(nDef)
        write(6,'(1x,a20,e11.4)')'with magnitude: ',nDefNorm
      endif
 
 !--------------------------------------------------------------------------------------------------
 ! compute the stress penalty of all interfaces
      do i = 1_pInt,3_pInt; do j = 1_pInt,3_pInt; do k = 1_pInt,3_pInt; do l = 1_pInt,3_pInt
        rPen(i,j,iGrain) = rPen(i,j,iGrain) + 0.5_pReal*(muGrain*bgGrain + muGNghb*bgGNghb)*prm%xiAlpha &
                                               *surfCorr(abs(intFace(1)))/prm%dAlpha(abs(intFace(1))) &
                                               *cosh(prm%ciAlpha*nDefNorm) &
                                               *0.5_pReal*nVect(l)*nDef(i,k)/nDefNorm*math_civita(k,l,j) &
                                               *tanh(nDefNorm/xSmoo_RGC)
      enddo; enddo;enddo; enddo
    enddo
    
    if (debugActive) then
     write(6,'(1x,a20,i2)')'Penalty of grain: ',iGrain
     write(6,*) transpose(rPen(1:3,1:3,iGrain))
    endif
 
  enddo
  end associate
 
 end subroutine stressPenalty
 
 
 !--------------------------------------------------------------------------------------------------
 !> @brief calculate stress-like penalty due to volume discrepancy 
 !--------------------------------------------------------------------------------------------------
 subroutine volumePenalty(vPen,vDiscrep,fDef,fAvg,ip,el) 
  use math, only: &
    math_det33, &
    math_inv33
  use numerics, only: &
    maxVolDiscr_RGC,&
    volDiscrMod_RGC,&
    volDiscrPow_RGC
 
  implicit none
  real(pReal), dimension (3,3,homogenization_maxNgrains), intent(out) :: vPen                        ! stress-like penalty due to volume
  real(pReal), intent(out)                  :: vDiscrep                                              ! total volume discrepancy
  real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in)  :: fDef                        ! deformation gradients
  real(pReal), dimension (3,3), intent(in)  :: fAvg                                                  ! overall deformation gradient
  integer(pInt), intent(in)                 :: ip,&                                                  ! integration point
    el 
  real(pReal), dimension (homogenization_maxNgrains) :: gVol
  integer(pInt) :: iGrain,nGrain
  logical :: debugActive

  debugActive = iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt &
                .and. debug_e == el .and. debug_i == ip

  nGrain = homogenization_Ngrains(material_homogenizationAt(el))
   
 !--------------------------------------------------------------------------------------------------
 ! compute the volumes of grains and of cluster
  vDiscrep = math_det33(fAvg)                                                                        ! compute the volume of the cluster
  do iGrain = 1_pInt,nGrain
    gVol(iGrain) = math_det33(fDef(1:3,1:3,iGrain))                                                  ! compute the volume of individual grains
    vDiscrep     = vDiscrep - gVol(iGrain)/real(nGrain,pReal)                                        ! calculate the difference/dicrepancy between
                                                                                                     ! the volume of the cluster and the the total volume of grains
  enddo
 
 !--------------------------------------------------------------------------------------------------
 ! calculate the stress and penalty due to volume discrepancy
  vPen      = 0.0_pReal
  do iGrain = 1_pInt,nGrain
    vPen(:,:,iGrain) = -1.0_pReal/real(nGrain,pReal)*volDiscrMod_RGC*volDiscrPow_RGC/maxVolDiscr_RGC* &
                       sign((abs(vDiscrep)/maxVolDiscr_RGC)**(volDiscrPow_RGC - 1.0),vDiscrep)* &
                       gVol(iGrain)*transpose(math_inv33(fDef(:,:,iGrain)))
 
    if (debugActive) then
      write(6,'(1x,a30,i2)')'Volume penalty of grain: ',iGrain
      write(6,*) transpose(vPen(:,:,iGrain))
    endif
  enddo

 end subroutine volumePenalty


 !--------------------------------------------------------------------------------------------------
 !> @brief compute the correction factor accouted for surface evolution (area change) due to 
 ! deformation
 !--------------------------------------------------------------------------------------------------
 function surfaceCorrection(avgF,instance,of)
  use math, only: &
    math_invert33, &
    math_mul33x33
  
  implicit none
  real(pReal), dimension(3)               :: surfaceCorrection
  real(pReal), dimension(3,3), intent(in) :: avgF                                                    !< average F
  integer(pInt),               intent(in) :: &
    instance, &
    of
  real(pReal), dimension(3,3)             :: invC
  real(pReal), dimension(3)               :: nVect
  real(pReal)  :: detF
  integer(pInt) :: i,j,iBase
  logical     ::  error
 
  call math_invert33(math_mul33x33(transpose(avgF),avgF),invC,detF,error)
 
  surfaceCorrection = 0.0_pReal
  do iBase = 1_pInt,3_pInt
    nVect = interfaceNormal([iBase,1_pInt,1_pInt,1_pInt],instance,of)
    do i = 1_pInt,3_pInt; do j = 1_pInt,3_pInt
      surfaceCorrection(iBase) = surfaceCorrection(iBase) + invC(i,j)*nVect(i)*nVect(j)              ! compute the component of (the inverse of) the stretch in the direction of the normal
    enddo; enddo
    surfaceCorrection(iBase) = sqrt(surfaceCorrection(iBase))*detF                                   ! get the surface correction factor (area contraction/enlargement)
  enddo
 
 end function surfaceCorrection


 !--------------------------------------------------------------------------------------------------
 !> @brief compute the equivalent shear and bulk moduli from the elasticity tensor
 !--------------------------------------------------------------------------------------------------
 function equivalentModuli(grainID,ip,el)
  use constitutive, only: &
    constitutive_homogenizedC
 
  implicit none
  real(pReal), dimension(2)    :: equivalentModuli
  integer(pInt), intent(in)    :: &
    grainID,&
    ip, &                                                                                            !< integration point number
    el                                                                                               !< element number
  real(pReal), dimension(6,6) :: elasTens
  real(pReal) :: &
    cEquiv_11, &
    cEquiv_12, &
    cEquiv_44
 
  elasTens = constitutive_homogenizedC(grainID,ip,el)
 
 !--------------------------------------------------------------------------------------------------
 ! compute the equivalent shear modulus after Turterltaub and Suiker, JMPS (2005)
  cEquiv_11 = (elasTens(1,1) + elasTens(2,2) + elasTens(3,3))/3.0_pReal
  cEquiv_12 = (elasTens(1,2) + elasTens(2,3) + elasTens(3,1) + &
               elasTens(1,3) + elasTens(2,1) + elasTens(3,2))/6.0_pReal
  cEquiv_44 = (elasTens(4,4) + elasTens(5,5) + elasTens(6,6))/3.0_pReal
  equivalentModuli(1) = 0.2_pReal*(cEquiv_11 - cEquiv_12) + 0.6_pReal*cEquiv_44
 
 !--------------------------------------------------------------------------------------------------
 ! obtain the length of Burgers vector (could be model dependend)
  equivalentModuli(2) = 2.5e-10_pReal
 
 end function equivalentModuli
 
 
 !--------------------------------------------------------------------------------------------------
 !> @brief calculating the grain deformation gradient (the same with 
 ! homogenization_RGC_partionDeformation, but used only for perturbation scheme)
 !--------------------------------------------------------------------------------------------------
 subroutine grainDeformation(F, avgF, instance, of)
  
  implicit none
  real(pReal),   dimension (3,3,homogenization_maxNgrains), intent(out) :: F                         !< partioned F per grain
  real(pReal),   dimension (3,3),                           intent(in)  :: avgF                      !< 
  integer(pInt),                                            intent(in)  :: &
    instance, &
    of
  real(pReal),   dimension (3) :: aVect,nVect
  integer(pInt), dimension (4) :: intFace
  integer(pInt), dimension (3) :: iGrain3
  integer(pInt) :: iGrain,iFace,i,j
 
  !--------------------------------------------------------------------------------------------------
  ! compute the deformation gradient of individual grains due to relaxations
  F = 0.0_pReal
  do iGrain = 1_pInt,product(param(instance)%Nconstituents)
    iGrain3 = grain1to3(iGrain,param(instance)%Nconstituents)
    do iFace = 1_pInt,6_pInt
      intFace = getInterface(iFace,iGrain3)
      aVect   = relaxationVector(intFace,instance,of)
      nVect   = interfaceNormal(intFace,instance,of)
      forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) &
        F(i,j,iGrain) = F(i,j,iGrain) + aVect(i)*nVect(j)                                           ! effective relaxations
    enddo
    F(1:3,1:3,iGrain) = F(1:3,1:3,iGrain) + avgF                                                    ! relaxed deformation gradient
  enddo
  
 end subroutine grainDeformation
 
end function homogenization_RGC_updateState


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities 
!--------------------------------------------------------------------------------------------------
subroutine homogenization_RGC_averageStressAndItsTangent(avgP,dAvgPdAvgF,P,dPdF,instance)
 use material, only: &
   homogenization_maxNgrains

 implicit none
 real(pReal), dimension (3,3),                               intent(out) :: avgP                    !< average stress at material point
 real(pReal), dimension (3,3,3,3),                           intent(out) :: dAvgPdAvgF              !< average stiffness at material point
 real(pReal), dimension (3,3,homogenization_maxNgrains),     intent(in)  :: P                       !< array of current grain stresses
 real(pReal), dimension (3,3,3,3,homogenization_maxNgrains), intent(in)  :: dPdF                    !< array of current grain stiffnesses
 integer(pInt),                                              intent(in)  :: instance

 avgP       = sum(P,3)   /real(product(param(instance)%Nconstituents),pReal)
 dAvgPdAvgF = sum(dPdF,5)/real(product(param(instance)%Nconstituents),pReal)

end subroutine homogenization_RGC_averageStressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion 
!--------------------------------------------------------------------------------------------------
pure function homogenization_RGC_postResults(ip,el) result(postResults)
 use material, only: &
   material_homogenizationAt, &
   homogenization_typeInstance,&
   mappingHomogenization
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number

 integer(pInt) instance,o,c,of
 real(pReal), dimension(sum(homogenization_RGC_sizePostResult(:,homogenization_typeInstance(material_homogenizationAt(el))))) :: &
   postResults

 instance = homogenization_typeInstance(material_homogenizationAt(el))
 associate(prm => param(instance))
 of = mappingHomogenization(1,ip,el)

 c = 0_pInt
 postResults = 0.0_pReal
 outputsLoop: do o = 1_pInt,size(prm%outputID)
   select case(prm%outputID(o))
     case (constitutivework_ID)
       postResults(c+1) = state(instance)%work(of)
       c = c + 1_pInt
     case (magnitudemismatch_ID)
       postResults(c+1:c+3) = state(instance)%mismatch(1:3,of)
       c = c + 3_pInt
     case (penaltyenergy_ID)
       postResults(c+1) = state(instance)%penaltyEnergy(of)
       c = c + 1_pInt
     case (volumediscrepancy_ID)
       postResults(c+1) = state(instance)%volumeDiscrepancy(of)
       c = c + 1_pInt
     case (averagerelaxrate_ID)
       postResults(c+1) = state(instance)%relaxationrate_avg(of)
       c = c + 1_pInt
     case (maximumrelaxrate_ID)
       postResults(c+1) = state(instance)%relaxationrate_max(of)
       c = c + 1_pInt
   end select
 enddo outputsLoop
 end associate
end function homogenization_RGC_postResults


!--------------------------------------------------------------------------------------------------
!> @brief collect relaxation vectors of an interface
!--------------------------------------------------------------------------------------------------
pure function relaxationVector(intFace,instance,of)

 implicit none
 integer(pInt),                intent(in) :: instance,of
 real(pReal),   dimension (3)             :: relaxationVector
 integer(pInt), dimension (4), intent(in) :: intFace                                                !< set of interface ID in 4D array (normal and position)
 integer(pInt) :: &
   iNum

!--------------------------------------------------------------------------------------------------
! collect the interface relaxation vector from the global state array
 relaxationVector = 0.0_pReal
 iNum = interface4to1(intFace,param(instance)%Nconstituents)                                        ! identify the position of the interface in global state array
 if (iNum > 0_pInt) relaxationVector = state(instance)%relaxationVector((3*iNum-2):(3*iNum),of)

end function relaxationVector


!--------------------------------------------------------------------------------------------------
!> @brief identify the normal of an interface 
!--------------------------------------------------------------------------------------------------
pure function interfaceNormal(intFace,instance,of)
 use math, only: &
   math_mul33x3
 
 implicit none
 real(pReal),   dimension (3)  ::            interfaceNormal
 integer(pInt), dimension (4), intent(in) :: intFace                                                !< interface ID in 4D array (normal and position)
 integer(pInt),                intent(in) :: &
   instance, &
   of
 integer(pInt) ::                            nPos

!--------------------------------------------------------------------------------------------------
! get the normal of the interface, identified from the value of intFace(1)
 interfaceNormal = 0.0_pReal
 nPos = abs(intFace(1))                                                                             ! identify the position of the interface in global state array
 interfaceNormal(nPos) = real(intFace(1)/abs(intFace(1)),pReal)                                     ! get the normal vector w.r.t. cluster axis

 interfaceNormal = math_mul33x3(dependentState(instance)%orientation(1:3,1:3,of),interfaceNormal)   ! map the normal vector into sample coordinate system (basis)

end function interfaceNormal


!--------------------------------------------------------------------------------------------------
!> @brief collect six faces of a grain in 4D (normal and position)
!--------------------------------------------------------------------------------------------------
pure function getInterface(iFace,iGrain3)

 implicit none
 integer(pInt), dimension (4) ::             getInterface
 integer(pInt), dimension (3), intent(in) :: iGrain3                                                !< grain ID in 3D array
 integer(pInt),                intent(in) :: iFace                                                  !< face index (1..6) mapped like (-e1,-e2,-e3,+e1,+e2,+e3) or iDir = (-1,-2,-3,1,2,3)
 integer(pInt) ::                            iDir
 
!* Direction of interface normal
 iDir = (int(real(iFace-1_pInt,pReal)/2.0_pReal,pInt)+1_pInt)*(-1_pInt)**iFace
 getInterface(1) = iDir
 
!--------------------------------------------------------------------------------------------------
! identify the interface position by the direction of its normal
 getInterface(2:4) = iGrain3
 if (iDir < 0_pInt) getInterface(1_pInt-iDir) = getInterface(1_pInt-iDir)-1_pInt                    ! to have a correlation with coordinate/position in real space

end function getInterface


!--------------------------------------------------------------------------------------------------
!> @brief map grain ID from in 1D (global array) to in 3D (local position)
!--------------------------------------------------------------------------------------------------
pure function grain1to3(grain1,nGDim)
 
 implicit none
 integer(pInt), dimension(3)             :: grain1to3
 integer(pInt),               intent(in) :: grain1                                                  !< grain ID in 1D array
 integer(pInt), dimension(3), intent(in) :: nGDim

 grain1to3 = 1_pInt + [mod((grain1-1_pInt),nGDim(1)), &
                       mod((grain1-1_pInt)/nGDim(1),nGDim(2)), &
                       (grain1-1_pInt)/(nGDim(1)*nGDim(2))]

end function grain1to3


!--------------------------------------------------------------------------------------------------
!> @brief map grain ID from in 3D (local position) to in 1D (global array)
!--------------------------------------------------------------------------------------------------
integer(pInt) pure function grain3to1(grain3,nGDim)

 implicit none
 integer(pInt), dimension(3), intent(in) :: grain3                                                  !< grain ID in 3D array (pos.x,pos.y,pos.z)
 integer(pInt), dimension(3), intent(in) :: nGDim

 grain3to1 = grain3(1) &
           + nGDim(1)*(grain3(2)-1_pInt) &
           + nGDim(1)*nGDim(2)*(grain3(3)-1_pInt)

end function grain3to1


!--------------------------------------------------------------------------------------------------
!> @brief maps interface ID from 4D (normal and local position) into 1D (global array)
!--------------------------------------------------------------------------------------------------
integer(pInt) pure function interface4to1(iFace4D, nGDim)
 
 implicit none
 integer(pInt), dimension(4), intent(in) :: iFace4D                                                 !< interface ID in 4D array (n.dir,pos.x,pos.y,pos.z)
 integer(pInt), dimension(3), intent(in) :: nGDim

 
 select case(abs(iFace4D(1)))

   case(1_pInt)
     if ((iFace4D(2) == 0_pInt) .or. (iFace4D(2) == nGDim(1))) then
       interface4to1 = 0_pInt
     else
       interface4to1 = iFace4D(3) + nGDim(2)*(iFace4D(4)-1_pInt) &
                     + nGDim(2)*nGDim(3)*(iFace4D(2)-1_pInt)
     endif

   case(2_pInt)
     if ((iFace4D(3) == 0_pInt) .or. (iFace4D(3) == nGDim(2))) then
       interface4to1 = 0_pInt
     else
       interface4to1 = iFace4D(4) + nGDim(3)*(iFace4D(2)-1_pInt) &
                     + nGDim(3)*nGDim(1)*(iFace4D(3)-1_pInt) &
                     + (nGDim(1)-1_pInt)*nGDim(2)*nGDim(3)                                          ! total number of interfaces normal //e1
     endif

   case(3_pInt)
     if ((iFace4D(4) == 0_pInt) .or. (iFace4D(4) == nGDim(3))) then
       interface4to1 = 0_pInt
     else
       interface4to1 = iFace4D(2) + nGDim(1)*(iFace4D(3)-1_pInt) &
                     + nGDim(1)*nGDim(2)*(iFace4D(4)-1_pInt) &
                     + (nGDim(1)-1_pInt)*nGDim(2)*nGDim(3) &                                        ! total number of interfaces normal //e1
                     + nGDim(1)*(nGDim(2)-1_pInt)*nGDim(3)                                          ! total number of interfaces normal //e2
     endif

   case default
     interface4to1 = -1_pInt
     
 end select

end function interface4to1


!--------------------------------------------------------------------------------------------------
!> @brief maps interface ID from 1D (global array) into 4D (normal and local position)
!--------------------------------------------------------------------------------------------------
pure function interface1to4(iFace1D, nGDim)
 
 implicit none
 integer(pInt), dimension(4)             :: interface1to4
 integer(pInt),               intent(in) :: iFace1D                                                            !< interface ID in 1D array
 integer(pInt), dimension(3), intent(in) :: nGDim
 integer(pInt), dimension (3)            :: nIntFace

!--------------------------------------------------------------------------------------------------
! compute the total number of interfaces, which ...
 nIntFace = [(nGDim(1)-1_pInt)*nGDim(2)*nGDim(3), &                                                 ! ... normal //e1
             nGDim(1)*(nGDim(2)-1_pInt)*nGDim(3), &                                                 ! ... normal //e2
             nGDim(1)*nGDim(2)*(nGDim(3)-1_pInt)]                                                   ! ... normal //e3

!--------------------------------------------------------------------------------------------------
! get the corresponding interface ID in 4D (normal and local position)
 if (iFace1D > 0 .and. iFace1D <= nIntFace(1)) then                                                 ! interface with normal //e1
   interface1to4(1) = 1_pInt
   interface1to4(3) = mod((iFace1D-1_pInt),nGDim(2))+1_pInt
   interface1to4(4) = mod(&
                                             int(&
                                                 real(iFace1D-1_pInt,pReal)/&
                                                 real(nGDim(2),pReal)&
                                                 ,pInt)&
                                             ,nGDim(3))+1_pInt
   interface1to4(2) = int(&
                                             real(iFace1D-1_pInt,pReal)/&
                                             real(nGDim(2),pReal)/&
                                             real(nGDim(3),pReal)&
                                             ,pInt)+1_pInt
 elseif (iFace1D > nIntFace(1) .and. iFace1D <= (nIntFace(2) + nIntFace(1))) then                   ! interface with normal //e2
   interface1to4(1) = 2_pInt
   interface1to4(4) = mod((iFace1D-nIntFace(1)-1_pInt),nGDim(3))+1_pInt
   interface1to4(2) = mod(&
                                             int(&
                                                 real(iFace1D-nIntFace(1)-1_pInt,pReal)/&
                                                 real(nGDim(3),pReal)&
                                                 ,pInt)&
                                              ,nGDim(1))+1_pInt
   interface1to4(3) = int(&
                                             real(iFace1D-nIntFace(1)-1_pInt,pReal)/&
                                             real(nGDim(3),pReal)/&
                                             real(nGDim(1),pReal)&
                                             ,pInt)+1_pInt
 elseif (iFace1D > nIntFace(2) + nIntFace(1) .and. iFace1D <= (nIntFace(3) + nIntFace(2) + nIntFace(1))) then ! interface with normal //e3
   interface1to4(1) = 3_pInt
   interface1to4(2) = mod((iFace1D-nIntFace(2)-nIntFace(1)-1_pInt),nGDim(1))+1_pInt
   interface1to4(3) = mod(&
                                             int(&
                                                 real(iFace1D-nIntFace(2)-nIntFace(1)-1_pInt,pReal)/&
                                                 real(nGDim(1),pReal)&
                                                 ,pInt)&
                                             ,nGDim(2))+1_pInt
   interface1to4(4) = int(&
                                             real(iFace1D-nIntFace(2)-nIntFace(1)-1_pInt,pReal)/&
                                             real(nGDim(1),pReal)/&
                                             real(nGDim(2),pReal)&
                                             ,pInt)+1_pInt
 endif

end function interface1to4


end module homogenization_RGC
