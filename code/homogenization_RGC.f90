!* $Id$
!*****************************************************
!*      Module: HOMOGENIZATION_RGC                   *
!*****************************************************
!* contains:                                         *
!*****************************************************

!	[rgc]
!	type            rgc
!	Ngrains         p x q x r (cluster)
!   (output)        Ngrains

MODULE homogenization_RGC

!*** Include other modules ***
 use prec, only: pReal,pInt
 implicit none

 character (len=*), parameter :: homogenization_RGC_label = 'rgc'
 
 integer(pInt),     dimension(:),     allocatable :: homogenization_RGC_sizeState, &
                                                     homogenization_RGC_sizePostResults
 integer(pInt),     dimension(:,:),   allocatable :: homogenization_RGC_Ngrains
 real(pReal),       dimension(:,:),   allocatable :: homogenization_RGC_xiAlpha, &
                                                     homogenization_RGC_ciAlpha
 character(len=64), dimension(:,:),   allocatable :: homogenization_RGC_output

CONTAINS
!****************************************
!* - homogenization_RGC_init
!* - homogenization_RGC_stateInit
!* - homogenization_RGC_deformationPartition
!* - homogenization_RGC_stateUpdate
!* - homogenization_RGC_averageStressAndItsTangent
!* - homogenization_RGC_postResults
!****************************************


!**************************************
!*      Module initialization         *
!**************************************
subroutine homogenization_RGC_init(&
   file  &    ! file pointer to material configuration
  )

 use prec, only: pInt, pReal
 use math, only: math_Mandel3333to66, math_Voigt66to3333
 use mesh, only: mesh_maxNips,mesh_NcpElems
 use IO
 use material
 integer(pInt), intent(in) :: file
 integer(pInt), parameter  :: maxNchunks = 4
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) section, maxNinstance, i,j,k,l, output
 character(len=64) tag
 character(len=1024) line
 
    write(6,*)
    write(6,'(a20,a20,a12)') '<<<+-  homogenization',homogenization_RGC_label,' init  -+>>>'
    write(6,*) '$Id$'
    write(6,*)

 maxNinstance = count(homogenization_type == homogenization_RGC_label)
 if (maxNinstance == 0) return

 allocate(homogenization_RGC_sizeState(maxNinstance));       homogenization_RGC_sizeState = 0_pInt
 allocate(homogenization_RGC_sizePostResults(maxNinstance)); homogenization_RGC_sizePostResults = 0_pInt
 allocate(homogenization_RGC_Ngrains(3,maxNinstance));       homogenization_RGC_Ngrains = 0_pInt
 allocate(homogenization_RGC_ciAlpha(3,maxNinstance));       homogenization_RGC_ciAlpha = 0.0_pReal
 allocate(homogenization_RGC_xiAlpha(3,maxNinstance));       homogenization_RGC_xiAlpha = 0.0_pReal
 allocate(homogenization_RGC_output(maxval(homogenization_Noutput),maxNinstance)); homogenization_RGC_output = ''
 
 rewind(file)
 line = ''
 section = 0
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= material_partHomogenization)     ! wind forward to <homogenization>
   read(file,'(a1024)',END=100) line
 enddo

 do                                                       ! read thru sections of phase part
   read(file,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                ! next section
     section = section + 1
     output = 0                                           ! reset output counter
   endif
   if (section > 0 .and. homogenization_type(section) == homogenization_RGC_label) then  ! one of my sections
     i = homogenization_typeInstance(section)             ! which instance of my type is present homogenization
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     select case(tag)
       case ('(output)')
         output = output + 1
         homogenization_RGC_output(output,i) = IO_lc(IO_stringValue(line,positions,2))
       case ('clustersize')
              homogenization_RGC_Ngrains(1,i) = IO_intValue(line,positions,2)
              homogenization_RGC_Ngrains(2,i) = IO_intValue(line,positions,3)
              homogenization_RGC_Ngrains(3,i) = IO_intValue(line,positions,4)
       case ('grainsizeparameter')
              homogenization_RGC_xiAlpha(1,i) = IO_floatValue(line,positions,2)
              homogenization_RGC_xiAlpha(2,i) = IO_floatValue(line,positions,3)
              homogenization_RGC_xiAlpha(3,i) = IO_floatValue(line,positions,4)
       case ('overproportionality')
              homogenization_RGC_ciAlpha(1,i) = IO_floatValue(line,positions,2)
              homogenization_RGC_ciAlpha(2,i) = IO_floatValue(line,positions,3)
              homogenization_RGC_ciAlpha(3,i) = IO_floatValue(line,positions,4)
     end select
   endif
 enddo

100 do i = 1,maxNinstance                                        ! sanity checks
 enddo

 do i = 1,maxNinstance
   do j = 1,maxval(homogenization_Noutput)
     select case(homogenization_RGC_output(j,i))
       case('constitutivework')
         homogenization_RGC_sizePostResults(i) = &
         homogenization_RGC_sizePostResults(i) + 1
       case('penaltyenergy')
         homogenization_RGC_sizePostResults(i) = &
         homogenization_RGC_sizePostResults(i) + 1
       case('magnitudemismatch')
         homogenization_RGC_sizePostResults(i) = &
         homogenization_RGC_sizePostResults(i) + 1
     end select
   enddo

   homogenization_RGC_sizeState(i) &
       = 3*(homogenization_RGC_Ngrains(1,i)-1)*homogenization_RGC_Ngrains(2,i)*homogenization_RGC_Ngrains(3,i) &
         + 3*homogenization_RGC_Ngrains(1,i)*(homogenization_RGC_Ngrains(2,i)-1)*homogenization_RGC_Ngrains(3,i) &
         + 3*homogenization_RGC_Ngrains(1,i)*homogenization_RGC_Ngrains(2,i)*(homogenization_RGC_Ngrains(3,i)-1) &
         + homogenization_RGC_sizePostResults(i)
 enddo

 return

endsubroutine


!*********************************************************************
!* initial homogenization state                                      *
!*********************************************************************
function homogenization_RGC_stateInit(myInstance)
 use prec, only: pReal,pInt
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: myInstance
 real(pReal), dimension(homogenization_RGC_sizeState(myInstance)) :: homogenization_RGC_stateInit

!* Open a debugging file 
 open(1978,file='homogenization_RGC_debugging.out')
 homogenization_RGC_stateInit = 0.0_pReal

 return
 
endfunction


!********************************************************************
! partition material point def grad onto constituents
!********************************************************************
subroutine homogenization_RGC_partitionDeformation(&
   F, &             ! partioned def grad per grain
!
   F0, &            ! initial partioned def grad per grain
   avgF, &          ! my average def grad
   state, &         ! my state
   ip, &            ! my integration point
   el  &            ! my element
  )
 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element,mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,homogenization_Ngrains,homogenization_typeInstance
 implicit none

!* Definition of variables
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(out) :: F
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in)  :: F0
 real(pReal), dimension (3,3), intent(in) :: avgF
 type(p_vec), intent(in) :: state
 integer(pInt), intent(in) :: ip,el
!
 real(pReal), dimension (3)   :: aVect,nVect
 integer(pInt), dimension (4) :: intFace
 integer(pInt), dimension (3) :: iGrain3
 integer(pInt) homID, iGrain,iFace,i,j
 logical RGCdebug
!
 integer(pInt), parameter :: nFace = 6
 
 RGCdebug = (el == 1 .and. ip == 1)

!* Debugging the overall deformation gradient
 if (RGCdebug) then
   write(1978,'(x,a32,i3,i3)')'Overall deformation gradient: '
   do i = 1,3
     write(1978,'(x,3(e10.4,x))')(avgF(i,j), j = 1,3)
   enddo
   write(1978,*)' '
 endif

!* Compute the deformation gradient of individual grains due to relaxations
 homID = homogenization_typeInstance(mesh_element(3,el))
 F = 0.0_pReal
 do iGrain = 1,homogenization_Ngrains(mesh_element(3,el))
   call homogenization_RGC_grain1to3(iGrain3,iGrain,homID)
   do iFace = 1,nFace
     call homogenization_RGC_getInterface(intFace,iFace,iGrain3)
     call homogenization_RGC_relaxationVector(aVect,intFace,state,homID)
     call homogenization_RGC_interfaceNormal(nVect,intFace)
     forall (i=1:3,j=1:3) &
     F(i,j,iGrain) = F(i,j,iGrain) + aVect(i)*nVect(j)               ! effective relaxations
   enddo
   F(:,:,iGrain) = F(:,:,iGrain) + avgF(:,:)                         ! relaxed deformation gradient
!* Debugging the grain deformation gradients
   if (RGCdebug) then
     write(1978,'(x,a32,x,i3)')'Deformation gradient of grain: ',iGrain
     do i = 1,3
       write(1978,'(x,3(e10.4,x))')(F(i,j,iGrain), j = 1,3)
     enddo
     write(1978,*)' '
   endif
 enddo
 
 return

endsubroutine

!********************************************************************
! update the internal state of the homogenization scheme
! and tell whether "done" and "happy" with result
!********************************************************************
function homogenization_RGC_updateState(&
   state, &         ! my state
!
   P, &             ! array of current grain stresses
   F, &             ! array of current grain deformation gradients
   F0, &            ! array of initial grain deformation gradients
   avgF, &          ! average deformation gradient
   dPdF, &          ! array of current grain stiffnesses
   ip, &            ! my integration point
   el  &            ! my element
  )

 use prec, only: pReal,pInt,p_vec
 use math, only: math_invert
 use mesh, only: mesh_element,mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,homogenization_typeInstance,homogenization_Ngrains
 use numerics, only: absTol_RGC,relTol_RGC,absMax_RGC,relMax_RGC,pPert_RGC

 implicit none

!* Definition of variables
 type(p_vec), intent(inout) :: state
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in) :: P,F,F0
 real(pReal), dimension (3,3,3,3,homogenization_maxNgrains), intent(in) :: dPdF
 real(pReal), dimension (3,3), intent(in) :: avgF
 integer(pInt), intent(in)    :: ip,el
 logical, dimension(2)        :: homogenization_RGC_updateState
 integer(pInt), dimension (4) :: intFaceN,intFaceP,faceID
 integer(pInt), dimension (3) :: nGDim,iGr3N,iGr3P,stresLoc
 integer(pInt), dimension (2) :: residLoc
 integer(pInt) homID,i1,i2,i3,iNum,i,j,nIntFaceTot,iGrN,iGrP,iMun,iFace,k,l,ival,ipert,iGrain
 real(pReal), dimension (3,3,homogenization_maxNgrains) :: R,pF,pR
 real(pReal), dimension (homogenization_maxNgrains)     :: NN,pNN
 real(pReal), dimension (3)   :: normP,normN,mornP,mornN
 real(pReal) residMax,stresMax,constitutiveWork,penaltyEnergy
 logical error,RGCdebug,RGCdebugJacobi
!
 integer(pInt), parameter :: nFace = 6
!
 real(pReal), dimension(:,:), allocatable :: tract,jmatrix,jnverse,smatrix,pmatrix
 real(pReal), dimension(:), allocatable   :: resid,relax,p_relax,p_resid

 RGCdebug = (el == 1 .and. ip == 1)
 RGCdebugJacobi = .false.

!* Get the dimension of the cluster (grains and interfaces)
 homID = homogenization_typeInstance(mesh_element(3,el))
 nGDim = homogenization_RGC_Ngrains(:,homID)
 nIntFaceTot = (nGDim(1)-1)*nGDim(2)*nGDim(3) + nGDim(1)*(nGDim(2)-1)*nGDim(3) &
               + nGDim(1)*nGDim(2)*(nGDim(3)-1)

!* Allocate the size of the arrays/matrices depending on the size of the cluster
 allocate(resid(3*nIntFaceTot)); resid = 0.0_pReal
 allocate(tract(nIntFaceTot,3)); tract = 0.0_pReal
 allocate(relax(3*nIntFaceTot)); relax = state%p(1:3*nIntFaceTot)
!* Debugging the obtained state
 if (RGCdebug) then
   write(1978,'(x,a30)')'Obtained state: '
   do i = 1,3*nIntFaceTot
     write(1978,'(x,2(e10.4,x))')state%p(i)
   enddo
   write(1978,*)' '
 endif

!* Stress-like penalty related to mismatch or incompatibility at interfaces
 call homogenization_RGC_stressPenalty(R,NN,F,ip,el,homID)
!* Debugging the mismatch, stress and penalty of grains
 if (RGCdebug) then
   do iGrain = 1,homogenization_Ngrains(mesh_element(3,el))
     write(1978,'(x,a30,x,i3,x,a4,x,e10.4)')'Mismatch magnitude of grain(',iGrain,') :',NN(iGrain)
     write(1978,*)' '
     write(1978,'(x,a30,x,i3)')'Stress and penalty of grain: ',iGrain
     do i = 1,3
       write(1978,'(x,3(e10.4,x),x,3(e10.4,x))')(P(i,j,iGrain), j = 1,3),(R(i,j,iGrain), j = 1,3)
     enddo
     write(1978,*)' '
   enddo
 endif

!* Compute the residual stress from the balance of traction at all (interior) interfaces
 do iNum = 1,nIntFaceTot
   call homogenization_RGC_interface1to4(faceID,iNum,homID)
!* Identify the left/bottom/back grain (-|N)
   iGr3N = faceID(2:4)
   call homogenization_RGC_grain3to1(iGrN,iGr3N,homID)
   call homogenization_RGC_getInterface(intFaceN,2*faceID(1),iGr3N)
   call homogenization_RGC_interfaceNormal(normN,intFaceN)            ! get the interface normal
!* Identify the right/up/front grain (+|P)
   iGr3P = iGr3N
   iGr3P(faceID(1)) = iGr3N(faceID(1))+1
   call homogenization_RGC_grain3to1(iGrP,iGr3P,homID)
   call homogenization_RGC_getInterface(intFaceP,2*faceID(1)-1,iGr3P)
   call homogenization_RGC_interfaceNormal(normP,intFaceP)            ! get the interface normal
   do i = 1,3                                                         ! compute the traction balance at the interface
   do j = 1,3
     tract(iNum,i) = tract(iNum,i) + (P(i,j,iGrP) + R(i,j,iGrP))*normP(j) &
                                   + (P(i,j,iGrN) + R(i,j,iGrN))*normN(j)
     resid(i+3*(iNum-1)) = tract(iNum,i)                              ! map into 1D residual array
   enddo
   enddo
!* Debugging the residual stress
   if (RGCdebug) then
     write(1978,'(x,a30,x,i3)')'Traction difference: ',iNum
     write(1978,'(x,3(e10.4,x))')(tract(iNum,j), j = 1,3)
     write(1978,*)' '
   endif
 enddo

!* Convergence check for stress residual
 stresMax = maxval(P)
 stresLoc = maxloc(P)
 residMax = maxval(tract)
 residLoc = maxloc(tract)
!* Debugging the convergent criteria
 if (RGCdebug) then
   write(1978,'(x,a)')' '
   write(1978,'(x,a)')'Residual check ...'
   write(1978,'(x,a15,x,e10.4,x,a7,i3,x,a12,i2,i2)')'Max stress: ',stresMax, &
              '@ grain',stresLoc(3),'in component',stresLoc(1),stresLoc(2)
   write(1978,'(x,a15,x,e10.4,x,a7,i3,x,a12,i2)')'Max residual: ',residMax, &
              '@ iface',residLoc(1),'in direction',residLoc(2)
 endif
 homogenization_RGC_updateState = .false.
!*** If convergence reached => done and happy
 if (residMax < relTol_RGC*stresMax .or. residMax < absTol_RGC) then 
   homogenization_RGC_updateState = .true.
   if (RGCdebug) then 
     write(1978,'(x,a55)')'... done and happy'
     write(1978,*)' '
   endif
   write(6,'(x,a,x,i3,x,a6,x,i3,x,a12)')'RGC_updateState: ip',ip,'| el',el,'converged :)'
!* Then compute/update the state for postResult, i.e., ...
!* ... the (bulk) constitutive work and penalty energy
   constitutiveWork = state%p(3*nIntFaceTot+1)
   penaltyEnergy    = state%p(3*nIntFaceTot+2)
   do iGrain = 1,homogenization_Ngrains(mesh_element(3,el))
     do i = 1,3
     do j = 1,3
       constitutiveWork = constitutiveWork + P(i,j,iGrain)*(F(i,j,iGrain) - F0(i,j,iGrain))
       penaltyEnergy    = penaltyEnergy    + R(i,j,iGrain)*(F(i,j,iGrain) - F0(i,j,iGrain))
     enddo
     enddo
   enddo
   state%p(3*nIntFaceTot+1) = constitutiveWork
   state%p(3*nIntFaceTot+2) = penaltyEnergy
!* ... the overall mismatch
   state%p(3*nIntFaceTot+3) = sum(NN)
   if (RGCdebug) then
     write(1978,'(x,a30,x,e10.4)')'Constitutive work: ',constitutiveWork
     write(1978,'(x,a30,x,e10.4)')'Penalty energy: ',penaltyEnergy
     write(1978,'(x,a30,x,e10.4)')'Magnitude mismatch: ',sum(NN)
     write(1978,*)' '
   endif
   deallocate(tract,resid,relax)
   return
!*** If residual blows-up => done but unhappy
 elseif (residMax > relMax_RGC*stresMax .or. residMax > absMax_RGC) then
   homogenization_RGC_updateState = (/.true.,.false./)
   if (RGCdebug) then
     write(1978,'(x,a55)')'... broken'
     write(1978,*)' '
   endif
   write(6,'(x,a,x,i3,x,a6,x,i3,x,a9)')'RGC_updateState: ip',ip,'| el',el,'broken :('
   deallocate(tract,resid,relax)
   return 
!*** Otherwise, proceed with computing the Jacobian and state update
 else
   if (RGCdebug) then
     write(1978,'(x,a55)')'... not yet done'
     write(1978,*)' '
   endif
 endif

!* Construct the Jacobian matrix of the constitutive stress tangent from dPdF
 allocate(smatrix(3*nIntFaceTot,3*nIntFaceTot)); smatrix = 0.0_pReal
 do iNum = 1,nIntFaceTot
   call homogenization_RGC_interface1to4(faceID,iNum,homID)
!* Identify the left/bottom/back grain (-|N)
   iGr3N = faceID(2:4)
   call homogenization_RGC_grain3to1(iGrN,iGr3N,homID)
   call homogenization_RGC_getInterface(intFaceN,2*faceID(1),iGr3N)
   call homogenization_RGC_interfaceNormal(normN,intFaceN)            ! get the interface normal
   do iFace = 1,nFace
     call homogenization_RGC_getInterface(intFaceN,iFace,iGr3N)
     call homogenization_RGC_interfaceNormal(mornN,intFaceN)          ! get influencing interfaces normal
     call homogenization_RGC_interface4to1(iMun,intFaceN,homID)
     if (iMun .gt. 0) then                                            ! get the tangent
       forall(i=1:3,j=1:3,k=1:3,l=1:3) &
       smatrix(3*(iNum-1)+i,3*(iMun-1)+j) = smatrix(3*(iNum-1)+i,3*(iMun-1)+j) + dPdF(i,k,j,l,iGrN)*normN(k)*mornN(l)
     endif
   enddo
!* Identify the right/up/front grain (+|P)
   iGr3P = iGr3N
   iGr3P(faceID(1)) = iGr3N(faceID(1))+1
   call homogenization_RGC_grain3to1(iGrP,iGr3P,homID)
   call homogenization_RGC_getInterface(intFaceP,2*faceID(1)-1,iGr3P)
   call homogenization_RGC_interfaceNormal(normP,intFaceP)            ! get the interface normal
   do iFace = 1,nFace
     call homogenization_RGC_getInterface(intFaceP,iFace,iGr3P)
     call homogenization_RGC_interfaceNormal(mornP,intFaceP)          ! get influencing interfaces normal
     call homogenization_RGC_interface4to1(iMun,intFaceP,homID)
     if (iMun .gt. 0) then                                            ! get the tangent
       forall(i=1:3,j=1:3,k=1:3,l=1:3) &
       smatrix(3*(iNum-1)+i,3*(iMun-1)+j) = smatrix(3*(iNum-1)+i,3*(iMun-1)+j) + dPdF(i,k,j,l,iGrP)*normP(k)*mornP(l)
     endif
   enddo
 enddo
!* Debugging the global Jacobian matrix of stress tangent
 if (RGCdebugJacobi) then
   write(1978,'(x,a30)')'Jacobian matrix of stress'
   do i = 1,3*nIntFaceTot
     write(1978,'(x,100(e10.4,x))')(smatrix(i,j), j = 1,3*nIntFaceTot)
   enddo
   write(1978,*)' '
 endif

!* Compute the Jacobian of the stress-like penalty (penalty tangent) using perturbation technique
 allocate(pmatrix(3*nIntFaceTot,3*nIntFaceTot)); pmatrix = 0.0_pReal
 allocate(p_relax(3*nIntFaceTot));               p_relax = 0.0_pReal
 allocate(p_resid(3*nIntFaceTot));               p_resid = 0.0_pReal
 do ipert = 1,3*nIntFaceTot
   p_relax = relax
   p_relax(ipert) = relax(ipert) + pPert_RGC                            ! perturb the relaxation vector
   state%p(1:3*nIntFaceTot) = p_relax
   call homogenization_RGC_grainDeformation(pF,F0,avgF,state,el)
   call homogenization_RGC_stressPenalty(pR,pNN,pF,ip,el,homID)
   p_resid = 0.0_pReal
   do iNum = 1,nIntFaceTot
     call homogenization_RGC_interface1to4(faceID,iNum,homID)
!* Identify the left/bottom/back grain (-|N)
     iGr3N = faceID(2:4)
     call homogenization_RGC_grain3to1(iGrN,iGr3N,homID)
     call homogenization_RGC_getInterface(intFaceN,2*faceID(1),iGr3N)
     call homogenization_RGC_interfaceNormal(normN,intFaceN)            ! get the corresponding normal
!* Identify the right/up/front grain (+|P)
     iGr3P = iGr3N
     iGr3P(faceID(1)) = iGr3N(faceID(1))+1
     call homogenization_RGC_grain3to1(iGrP,iGr3P,homID)
     call homogenization_RGC_getInterface(intFaceP,2*faceID(1)-1,iGr3P)
     call homogenization_RGC_interfaceNormal(normP,intFaceP)            ! get the corresponding normal
!* Compute the perturbed traction at interface
     do i = 1,3
     do j = 1,3
       p_resid(i+3*(iNum-1)) = p_resid(i+3*(iNum-1)) + (pR(i,j,iGrP) - R(i,j,iGrP))*normP(j) &
                                                     + (pR(i,j,iGrN) - R(i,j,iGrN))*normN(j)
     enddo
     enddo
   enddo
   pmatrix(:,ipert) = p_resid/pPert_RGC
 enddo
!* Debugging the global Jacobian matrix of penalty tangent
 if (RGCdebugJacobi) then
   write(1978,'(x,a30)')'Jacobian matrix of penalty'
   do i = 1,3*nIntFaceTot
     write(1978,'(x,100(e10.4,x))')(pmatrix(i,j), j = 1,3*nIntFaceTot)
   enddo
   write(1978,*)' '
 endif

!* The overall Jacobian matrix (due to constitutive and penalty tangents)
 allocate(jmatrix(3*nIntFaceTot,3*nIntFaceTot)); jmatrix = smatrix + pmatrix
 if (RGCdebugJacobi) then
   write(1978,'(x,a30)')'Jacobian matrix (total)'
   do i = 1,3*nIntFaceTot
     write(1978,'(x,100(e10.4,x))')(jmatrix(i,j), j = 1,3*nIntFaceTot)
   enddo
   write(1978,*)' '
 endif
 allocate(jnverse(3*nIntFaceTot,3*nIntFaceTot)); jnverse = 0.0_pReal
 call math_invert(3*nIntFaceTot,jmatrix,jnverse,ival,error)
!* Debugging the inverse Jacobian matrix
 if (RGCdebugJacobi) then
   write(1978,'(x,a30)')'Jacobian inverse'
   do i = 1,3*nIntFaceTot
     write(1978,'(x,100(e10.4,x))')(jnverse(i,j), j = 1,3*nIntFaceTot)
   enddo
   write(1978,*)' '
 endif

!* Calculate the state update (i.e., new relaxation vectors)
 forall(i=1:3*nIntFaceTot,j=1:3*nIntFaceTot) relax(i) = relax(i) - jnverse(i,j)*resid(j)
 state%p(1:3*nIntFaceTot) = relax
!* Debugging the return state
 if (RGCdebugJacobi) then
   write(1978,'(x,a30)')'Returned state: '
   do i = 1,3*nIntFaceTot
     write(1978,'(x,2(e10.4,x))')state%p(i)
   enddo
   write(1978,*)' '
 endif

 deallocate(tract,resid,jmatrix,jnverse,relax,pmatrix,smatrix,p_relax,p_resid)
 return 
 
endfunction

!********************************************************************
! derive average stress and stiffness from constituent quantities
!********************************************************************
subroutine homogenization_RGC_averageStressAndItsTangent(&
   avgP, &          ! average stress at material point
   dAvgPdAvgF, &    ! average stiffness at material point
!
   P, &             ! array of current grain stresses
   dPdF, &          ! array of current grain stiffnesses
   ip, &            ! my integration point
   el  &            ! my element
  )

 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element,mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains, homogenization_Ngrains,homogenization_typeInstance
 use math, only: math_Plain3333to99
 implicit none

!* Definition of variables
 real(pReal), dimension (3,3), intent(out) :: avgP
 real(pReal), dimension (3,3,3,3), intent(out) :: dAvgPdAvgF
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in) :: P
 real(pReal), dimension (3,3,3,3,homogenization_maxNgrains), intent(in) :: dPdF
 real(pReal), dimension (9,9) :: dPdF99
 integer(pInt), intent(in) :: ip,el
!
 logical homogenization_RGC_stateUpdate,RGCdebug
 integer(pInt) homID, i, j, Ngrains, iGrain

 RGCdebug = .false. !(ip == 1 .and. el == 1)

 homID = homogenization_typeInstance(mesh_element(3,el))
!* Debugging the grain tangent
 if (RGCdebug) then
   do iGrain = 1,homogenization_Ngrains(mesh_element(3,el))
     dPdF99 = math_Plain3333to99(dPdF(:,:,:,:,iGrain))
     write(1978,'(x,a30,x,i3)')'Stress tangent of grain: ',iGrain
     do i = 1,9
       write(1978,'(x,(e10.4,x))') (dPdF99(i,j), j = 1,9)
     enddo
     write(1978,*)' '
   enddo
 endif
 Ngrains = homogenization_Ngrains(mesh_element(3,el))
 avgP = sum(P,3)/dble(Ngrains)
 dAvgPdAvgF = sum(dPdF,5)/dble(Ngrains)

 return

endsubroutine

!********************************************************************
! derive average stress and stiffness from constituent quantities
!********************************************************************
function homogenization_RGC_averageTemperature(&
   Temperature, &   ! temperature
   ip, &            ! my integration point
   el  &            ! my element
  )

 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element,mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains, homogenization_Ngrains
 implicit none

!* Definition of variables
 real(pReal), dimension (homogenization_maxNgrains), intent(in) :: Temperature
 integer(pInt), intent(in) :: ip,el
 real(pReal) homogenization_RGC_averageTemperature
 integer(pInt) homID, i, Ngrains

! homID = homogenization_typeInstance(mesh_element(3,el))        ! <<not required at the moment>>
 Ngrains = homogenization_Ngrains(mesh_element(3,el))
 homogenization_RGC_averageTemperature = sum(Temperature(1:Ngrains))/dble(Ngrains)

 return

endfunction

!********************************************************************
! return array of homogenization results for post file inclusion
!********************************************************************
pure function homogenization_RGC_postResults(&
   state, &         ! my state
   ip, &            ! my integration point
   el  &            ! my element
  )

 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element
 use material, only: homogenization_typeInstance,homogenization_Noutput,homogenization_Ngrains
 implicit none

!* Definition of variables
 type(p_vec), intent(in) :: state
 integer(pInt), intent(in) :: ip,el
!
 integer(pInt) homID,o,c,nIntFaceTot
 real(pReal), dimension(homogenization_RGC_sizePostResults(homogenization_typeInstance(mesh_element(3,el)))) :: &
   homogenization_RGC_postResults

 homID = homogenization_typeInstance(mesh_element(3,el))
 nIntFaceTot = (homogenization_RGC_Ngrains(1,homID)-1)*homogenization_RGC_Ngrains(2,homID)*homogenization_RGC_Ngrains(3,homID) + & 
               homogenization_RGC_Ngrains(1,homID)*(homogenization_RGC_Ngrains(2,homID)-1)*homogenization_RGC_Ngrains(3,homID) + &
               homogenization_RGC_Ngrains(1,homID)*homogenization_RGC_Ngrains(2,homID)*(homogenization_RGC_Ngrains(3,homID)-1)

 c = 0_pInt
 homogenization_RGC_postResults = 0.0_pReal
 do o = 1,homogenization_Noutput(mesh_element(3,el))
   select case(homogenization_RGC_output(o,homID))
     case('constitutivework')
       homogenization_RGC_postResults(c+1) = state%p(3*nIntFaceTot+1)
       c = c + 1
     case('penaltyenergy')
       homogenization_RGC_postResults(c+1) = state%p(3*nIntFaceTot+2)
       c = c + 1
     case('magnitudemismatch')
       homogenization_RGC_postResults(c+1) = state%p(3*nIntFaceTot+3)
       c = c + 1
   end select
 enddo
 
 return

endfunction

!********************************************************************
! subroutine to calculate stress-like penalty due to mismatch
!********************************************************************
subroutine homogenization_RGC_stressPenalty(&
   rPen, &          ! stress-like penalty
   nMis, &          ! total amount of mismatch
!
   fDef, &          ! relaxation vectors
   ip, &            ! integration point
   el, &            ! element
   homID &          ! homogenization ID
  )
 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element
 use constitutive, only: constitutive_homogenizedC
 use math, only: math_civita
 use material, only: homogenization_maxNgrains,homogenization_Ngrains
 use numerics, only: xSmoo_RGC
 implicit none

!* Definition of variables
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(out) :: rPen
 real(pReal), dimension (homogenization_maxNgrains), intent(out)     :: nMis
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in)  :: fDef
 integer(pInt), intent(in)    :: ip,el
 integer(pInt), dimension (4) :: intFace
 integer(pInt), dimension (3) :: iGrain3,iGNghb3,nGDim
 real(pReal), dimension (3,3) :: gDef,nDef
 real(pReal), dimension (3)   :: nVect
 integer(pInt) homID,iGrain,iGNghb,iFace,i,j,k,l
 real(pReal) muGrain,muGNghb,nDefNorm
! 
 integer(pInt), parameter :: nFace = 6
 real(pReal), parameter   :: nDefToler = 1.0e-10

 nGDim = homogenization_RGC_Ngrains(:,homID)

 rPen = 0.0_pReal
 nMis = 0.0_pReal
 do iGrain = 1,homogenization_Ngrains(mesh_element(3,el))
   call homogenization_RGC_equivalentShearMod(muGrain,constitutive_homogenizedC(iGrain,ip,el))
   call homogenization_RGC_grain1to3(iGrain3,iGrain,homID)
!* Compute the mismatch tensor at all six interfaces
   do iFace = 1,nFace
     call homogenization_RGC_getInterface(intFace,iFace,iGrain3)
     call homogenization_RGC_interfaceNormal(nVect,intFace)             ! get the interface normal
     iGNghb3 = iGrain3                                                  ! identify the grain neighbor
     iGNghb3(abs(intFace(1))) = iGNghb3(abs(intFace(1))) + int(dble(intFace(1))/dble(abs(intFace(1))))
!* The grain periodicity along e1
     if (iGNghb3(1) < 1)        iGNghb3(1) = nGDim(1)
     if (iGNghb3(1) > nGDim(1)) iGNghb3(1) = 1
!* The grain periodicity along e2
     if (iGNghb3(2) < 1)        iGNghb3(2) = nGDim(2)
     if (iGNghb3(2) > nGDim(2)) iGNghb3(2) = 1
!* The grain periodicity along e3
     if (iGNghb3(3) < 1)        iGNghb3(3) = nGDim(3)
     if (iGNghb3(3) > nGDim(3)) iGNghb3(3) = 1
     call homogenization_RGC_grain3to1(iGNghb,iGNghb3,homID)            ! get the grain neighbor
     call homogenization_RGC_equivalentShearMod(muGNghb,constitutive_homogenizedC(iGNghb,ip,el))
     gDef = 0.5_pReal*(fDef(:,:,iGNghb) - fDef(:,:,iGrain))             ! difference in F with the neighbor
     nDefNorm = 0.0_pReal
     nDef = 0.0_pReal
     do i = 1,3
     do j = 1,3
       do k = 1,3
       do l = 1,3
         nDef(i,j) = nDef(i,j) - nVect(k)*gDef(i,l)*math_civita(j,k,l)  ! compute the interface mismatch tensor
       enddo
       enddo
       nDefNorm = nDefNorm + nDef(i,j)*nDef(i,j)
     enddo
     enddo
     nDefNorm = max(nDefToler,sqrt(nDefNorm))                           ! zero mismatch approximation if too small
!* Debugging the mismatch tensor
!      if (ip == 1 .and. el == 1) then
!        write(1978,'(x,a20,i2,x,a20,x,i3)')'Mismatch to face: ',intFace(1),'neighbor grain: ',iGNghb
!        do i = 1,3
!          write(1978,'(x,3(e10.4,x))')(nDef(i,j), j = 1,3)
!        enddo
!        write(1978,'(x,a20,e10.4))')'with magnitude: ',nDefNorm
!      endif
!* Compute the stress-like penalty from all six interfaces
     do i = 1,3
     do j = 1,3
       do k = 1,3
       do l = 1,3
         rPen(i,j,iGrain) = rPen(i,j,iGrain) + 0.5_pReal*(muGrain + muGNghb)/homogenization_RGC_xiAlpha(abs(intFace(1)),homID) &
                                               *cosh(homogenization_RGC_ciAlpha(abs(intFace(1)),homID)*nDefNorm) &
                                               *0.5_pReal*nVect(l)*nDef(i,k)/nDefNorm*math_civita(k,l,j) &
                                               *tanh(nDefNorm/xSmoo_RGC)
       enddo
       enddo
     enddo
     enddo
!* Total amount of mismatch experienced by the grain (at all six interfaces)
     nMis(iGrain) = nMis(iGrain) + nDefNorm
   enddo
!* Debugging the stress-like penalty
!    if (ip == 1 .and. el == 1) then
!      write(1978,'(x,a20,i2)')'Penalty of grain: ',iGrain
!      do i = 1,3
!        write(1978,'(x,3(e10.4,x))')(rPen(i,j,iGrain), j = 1,3)
!      enddo
!    endif
 enddo

 return

endsubroutine

!********************************************************************
! subroutine to compute the equivalent shear modulus from the elasticity tensor
!********************************************************************
subroutine homogenization_RGC_equivalentShearMod(&
   shearMod, &       ! equivalent (isotropic) shear modulus
!
   elasTens &        ! elasticity tensor in Mandel notation
  )

 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element,mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_typeInstance
 implicit none

!* Definition of variables
 real(pReal), dimension (6,6), intent(in) :: elasTens
 real(pReal), intent(out)                 :: shearMod
 real(pReal) cEquiv_11,cEquiv_12,cEquiv_44

!* Compute the equivalent shear modulus using Turterltaub and Suiker, JMPS (2005)
 cEquiv_11 = (elasTens(1,1) + elasTens(2,2) + elasTens(3,3))/3.0_pReal
 cEquiv_12 = (elasTens(1,2) + elasTens(2,3) + elasTens(3,1) + &
              elasTens(1,3) + elasTens(2,1) + elasTens(3,2))/6.0_pReal
 cEquiv_44 = (elasTens(4,4) + elasTens(5,5) + elasTens(6,6))/3.0_pReal
 shearMod = 0.2_pReal*(cEquiv_11 - cEquiv_12) + 0.6_pReal*cEquiv_44

 return

endsubroutine

!********************************************************************
! subroutine to collect relaxation vectors of an interface
!********************************************************************
subroutine homogenization_RGC_relaxationVector(&
   aVect, &          ! relaxation vector of the interface
!
   intFace, &        ! set of interface ID in 4D array (normal and position)
   state, &          ! set of global relaxation vectors
   homID &           ! homogenization ID
  )

 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element,mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_typeInstance
 implicit none

!* Definition of variables
 real(pReal), dimension (3), intent(out)  :: aVect
 integer(pInt), dimension (4), intent(in) :: intFace
 type(p_vec), intent(in)      :: state
 integer(pInt), dimension (3) :: nGDim
 integer(pInt) iNum,homID

!* Collect the interface relaxation vector from the global state array
 aVect = 0.0_pReal
 nGDim = homogenization_RGC_Ngrains(:,homID)
 call homogenization_RGC_interface4to1(iNum,intFace,homID)     ! Get the position in global state array
 if (iNum .gt. 0_pInt) aVect = state%p((3*iNum-2):(3*iNum))    ! Collect the corresponding entries

 return

endsubroutine

!********************************************************************
! subroutine to identify the normal of an interface
!********************************************************************
subroutine homogenization_RGC_interfaceNormal(&
   nVect, &          ! interface normal
!
   intFace &         ! interface ID in 4D array (normal and position)
  )

 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element,mesh_NcpElems,mesh_maxNips
 implicit none

!* Definition of variables
 real(pReal), dimension (3), intent(out)  :: nVect
 integer(pInt), dimension (4), intent(in) :: intFace
 integer(pInt) nPos

!* Get the normal of the interface, identified from the value of intFace(1)
 nVect = 0.0_pReal
 nPos = abs(intFace(1))
 nVect(nPos) = intFace(1)/abs(intFace(1))

 return

endsubroutine

!********************************************************************
! subroutine to collect six faces of a grain in 4D (normal and position)
!********************************************************************
subroutine homogenization_RGC_getInterface(&
   intFace, &       ! interface ID in 4D (normal and position)
!
   iFace, &         ! number of faces of grain
   iGrain3 &        ! grain ID in 3D array
  )
 use prec, only: pReal,pInt,p_vec
 implicit none

!* Definition of variables
 integer(pInt), dimension (4), intent(out)   :: intFace
 integer(pInt), dimension (3), intent(in)    :: iGrain3
 integer(pInt), intent(in) :: iFace
 integer(pInt) iDir
 
!* Direction of interface normal
 iDir = (int(dble(iFace-1)/2.0_pReal)+1)*(-1_pInt)**iFace
 intFace(1) = iDir
 
!* Identify the interface position by the direction of its normal
 intFace(2:4) = iGrain3(:)
 if (iDir .eq. -1_pInt) intFace(2) = intFace(2)-1
 if (iDir .eq. -2_pInt) intFace(3) = intFace(3)-1
 if (iDir .eq. -3_pInt) intFace(4) = intFace(4)-1

 return

endsubroutine

!********************************************************************
! subroutine to map grain ID from in 1D (array) to in 3D (position)
!********************************************************************
subroutine homogenization_RGC_grain1to3(&
   grain3, &        ! grain ID in 3D array (pos.x,pos.y,pos.z)
!
   grain1, &        ! grain ID in 1D array
   homID &          ! homogenization ID
  )

 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element,mesh_NcpElems,mesh_maxNips
 implicit none

!* Definition of variables
 integer(pInt), dimension (3), intent(out) :: grain3
 integer(pInt), intent(in)    :: grain1,homID
 integer(pInt), dimension (3) :: nGDim

!* Get the grain position
 nGDim = homogenization_RGC_Ngrains(:,homID)
 grain3(3) = int(dble(grain1-1)/dble(nGDim(1))/dble(nGDim(2)))+1
 grain3(2) = mod(int(dble(grain1-1)/dble(nGDim(1))),nGDim(2))+1
 grain3(1) = mod((grain1-1),nGDim(1))+1

 return

endsubroutine

!********************************************************************
! subroutine to map grain ID from in 3D (position) to in 1D (array)
!********************************************************************
subroutine homogenization_RGC_grain3to1(&
   grain1, &       ! grain ID in 1D array
!
   grain3, &       ! grain ID in 3D array (pos.x,pos.y,pos.z)
   homID &         ! homogenization ID
  )

 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element,mesh_NcpElems,mesh_maxNips
 implicit none

!* Definition of variables
 integer(pInt), dimension (3), intent(in) :: grain3
 integer(pInt), intent(out)   :: grain1
 integer(pInt), dimension (3) :: nGDim
 integer(pInt) homID

!* Get the grain ID
 nGDim = homogenization_RGC_Ngrains(:,homID)
 grain1 = grain3(1) + nGDim(1)*(grain3(2)-1) + nGDim(1)*nGDim(2)*(grain3(3)-1)

 return

endsubroutine

!********************************************************************
! subroutine to map interface ID from 4D (normal and position) into 1D (array)
!********************************************************************
subroutine homogenization_RGC_interface4to1(&
   iFace1D, &        ! interface ID in 1D array
!
   iFace4D, &        ! interface ID in 4D array (n.dir,pos.x,pos.y,pos.z)
   homID &           ! homogenization ID
  )

 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element,mesh_NcpElems,mesh_maxNips
 implicit none

!* Definition of variables
 integer(pInt), dimension (4), intent(in) :: iFace4D
 integer(pInt), intent(out)   :: iFace1D
 integer(pInt), dimension (3) :: nGDim,nIntFace
 integer(pInt) homID

 nGDim = homogenization_RGC_Ngrains(:,homID)
!* Get the number of interfaces, which ...
 nIntFace(1) = (nGDim(1)-1)*nGDim(2)*nGDim(3)      ! ... normal //e1
 nIntFace(2) = nGDim(1)*(nGDim(2)-1)*nGDim(3)      ! ... normal //e2
 nIntFace(3) = nGDim(1)*nGDim(2)*(nGDim(3)-1)      ! ... normal //e3

!* For interface with normal //e1
 if (abs(iFace4D(1)) == 1_pInt) then
   iFace1D = iFace4D(3) + nGDim(2)*(iFace4D(4)-1) + nGDim(2)*nGDim(3)*(iFace4D(2)-1)
   if ((iFace4D(2) == 0_pInt) .or. (iFace4D(2) == nGDim(1))) iFace1D = 0_pInt
!* For interface with normal //e2 
 elseif (abs(iFace4D(1)) == 2_pInt) then
   iFace1D = iFace4D(4) + nGDim(3)*(iFace4D(2)-1) + nGDim(3)*nGDim(1)*(iFace4D(3)-1) + nIntFace(1)
   if ((iFace4D(3) == 0_pInt) .or. (iFace4D(3) == nGDim(2))) iFace1D = 0_pInt
!* For interface with normal //e3 
 elseif (abs(iFace4D(1)) == 3_pInt) then
   iFace1D = iFace4D(2) + nGDim(1)*(iFace4D(3)-1) + nGDim(1)*nGDim(2)*(iFace4D(4)-1) + nIntFace(1) + nIntFace(2)
   if ((iFace4D(4) == 0_pInt) .or. (iFace4D(4) == nGDim(3))) iFace1D = 0_pInt
 endif

 return

endsubroutine

!********************************************************************
! subroutine to map interface ID from 1D (array) into 4D (normal and position)
!********************************************************************
subroutine homogenization_RGC_interface1to4(&
   iFace4D, &        ! interface ID in 4D array (n.dir,pos.x,pos.y,pos.z)
!
   iFace1D, &        ! interface ID in 1D array
   homID &           ! homogenization ID
  )

 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element,mesh_NcpElems,mesh_maxNips
 implicit none

!* Definition of variables
 integer(pInt), dimension (4), intent(out) :: iFace4D
 integer(pInt), intent(in)    :: iFace1D
 integer(pInt), dimension (3) :: nGDim,nIntFace
 integer(pInt) homID

 nGDim = homogenization_RGC_Ngrains(:,homID)
!* Get the number of interfaces, which ...
 nIntFace(1) = (nGDim(1)-1)*nGDim(2)*nGDim(3)      ! ... normal //e1
 nIntFace(2) = nGDim(1)*(nGDim(2)-1)*nGDim(3)      ! ... normal //e2
 nIntFace(3) = nGDim(1)*nGDim(2)*(nGDim(3)-1)      ! ... normal //e3

!* For interface ID between 1 and nIntFace(1)
 if (iFace1D > 0 .and. iFace1D <= nIntFace(1)) then
   iFace4D(1) = 1
   iFace4D(3) = mod((iFace1D-1),nGDim(2))+1
   iFace4D(4) = mod(int(dble(iFace1D-1)/dble(nGDim(2))),nGDim(3))+1
   iFace4D(2) = int(dble(iFace1D-1)/dble(nGDim(2))/dble(nGDim(3)))+1
!* For interface ID between nIntFace(1) and nIntFace(1) + nIntFace(2)
 elseif (iFace1D > nIntFace(1) .and. iFace1D <= (nIntFace(2) + nIntFace(1))) then
   iFace4D(1) = 2
   iFace4D(4) = mod((iFace1D-nIntFace(1)-1),nGDim(3))+1
   iFace4D(2) = mod(int(dble(iFace1D-nIntFace(1)-1)/dble(nGDim(3))),nGDim(1))+1
   iFace4D(3) = int(dble(iFace1D-nIntFace(1)-1)/dble(nGDim(3))/dble(nGDim(1)))+1
!* For interface ID between nIntFace(1) + nIntFace(2) and nIntFace(1) + nIntFace(2) + nIntFace(3)
 elseif (iFace1D > nIntFace(2) + nIntFace(1) .and. iFace1D <= (nIntFace(3) + nIntFace(2) + nIntFace(1))) then
   iFace4D(1) = 3
   iFace4D(2) = mod((iFace1D-nIntFace(2)-nIntFace(1)-1),nGDim(1))+1
   iFace4D(3) = mod(int(dble(iFace1D-nIntFace(2)-nIntFace(1)-1)/dble(nGDim(1))),nGDim(2))+1
   iFace4D(4) = int(dble(iFace1D-nIntFace(2)-nIntFace(1)-1)/dble(nGDim(1))/dble(nGDim(2)))+1
 endif

 return

endsubroutine

!********************************************************************
! calculating the grain deformation gradient
!********************************************************************
subroutine homogenization_RGC_grainDeformation(&
   F, &             ! partioned def grad per grain
!
   F0, &            ! initial partioned def grad per grain
   avgF, &          ! my average def grad
   state, &         ! my state
   el  &            ! my element
  )
 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element
 use material, only: homogenization_maxNgrains,homogenization_Ngrains,homogenization_typeInstance
 implicit none

!* Definition of variables
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(out) :: F
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in)  :: F0
 real(pReal), dimension (3,3), intent(in) :: avgF
 type(p_vec), intent(in) :: state
 integer(pInt), intent(in) :: el
!
 real(pReal), dimension (3)   :: aVect,nVect
 integer(pInt), dimension (4) :: intFace
 integer(pInt), dimension (3) :: iGrain3
 integer(pInt) homID, iGrain,iFace,i,j
!
 integer(pInt), parameter :: nFace = 6

!* Compute the deformation gradient of individual grains due to relaxations
 homID = homogenization_typeInstance(mesh_element(3,el))
 F = 0.0_pReal
 do iGrain = 1,homogenization_Ngrains(mesh_element(3,el))
   call homogenization_RGC_grain1to3(iGrain3,iGrain,homID)
   do iFace = 1,nFace
     call homogenization_RGC_getInterface(intFace,iFace,iGrain3)
     call homogenization_RGC_relaxationVector(aVect,intFace,state,homID)
     call homogenization_RGC_interfaceNormal(nVect,intFace)
     forall (i=1:3,j=1:3) &
     F(i,j,iGrain) = F(i,j,iGrain) + aVect(i)*nVect(j)               ! effective relaxations
   enddo
   F(:,:,iGrain) = F(:,:,iGrain) + avgF(:,:)                         ! relaxed deformation gradient
 enddo
 
 return

endsubroutine

END MODULE
