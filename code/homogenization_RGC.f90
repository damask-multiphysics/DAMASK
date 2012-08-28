! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!##############################################################
! $Id$
!*****************************************************
!*      Module: HOMOGENIZATION_RGC                   *
!*****************************************************
!* contains:                                         *
!*****************************************************

! [rgc]
!  type            rgc
!  Ngrains         p x q x r (cluster)
!   (output)        Ngrains

MODULE homogenization_RGC

!*** Include other modules ***
 use prec, only: pReal,pInt

 implicit none
 character (len=*), parameter :: homogenization_RGC_label = 'rgc'
 
 integer(pInt),     dimension(:),       allocatable :: homogenization_RGC_sizeState, &
                                                       homogenization_RGC_sizePostResults
 integer(pInt),     dimension(:,:),     allocatable,target :: homogenization_RGC_sizePostResult
 integer(pInt),     dimension(:,:),     allocatable :: homogenization_RGC_Ngrains
 real(pReal),       dimension(:,:),     allocatable :: homogenization_RGC_dAlpha, &
                                                       homogenization_RGC_angles
 real(pReal),       dimension(:,:,:,:), allocatable :: homogenization_RGC_orientation
 real(pReal),       dimension(:),       allocatable :: homogenization_RGC_xiAlpha, &
                                                       homogenization_RGC_ciAlpha
 character(len=64), dimension(:,:),     allocatable,target :: homogenization_RGC_output             ! name of each post result output

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
   myFile  &    ! file pointer to material configuration
  )

 use, intrinsic :: iso_fortran_env                                ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: debug_level, &
                  debug_homogenization, &
                  debug_levelBasic, &
                  debug_levelExtensive
 use math,  only: math_Mandel3333to66,&
                  math_Voigt66to3333, &
                  math_I3, &
                  math_sampleRandomOri,&
                  math_EulerToR,&
                  INRAD
 use mesh,  only: mesh_maxNips,mesh_NcpElems,mesh_element,FE_Nips
 use IO
 use material

 implicit none
 integer(pInt), intent(in) :: myFile
 integer(pInt), parameter  :: maxNchunks = 4_pInt
 integer(pInt), dimension(1_pInt+2_pInt*maxNchunks) :: positions
 integer(pInt) section, maxNinstance, i,j,e, output, mySize, myInstance
 character(len=64) tag
 character(len=1024) line

!$OMP CRITICAL (write2out)
   write(6,*)
   write(6,*) '<<<+-  homogenization_',trim(homogenization_RGC_label),' init  -+>>>'
   write(6,*) '$Id$'
#include "compilation_info.f90"
!$OMP END CRITICAL (write2out)

 maxNinstance = int(count(homogenization_type == homogenization_RGC_label),pInt)
 if (maxNinstance == 0) return

 allocate(homogenization_RGC_sizeState(maxNinstance));       homogenization_RGC_sizeState = 0_pInt
 allocate(homogenization_RGC_sizePostResults(maxNinstance)); homogenization_RGC_sizePostResults = 0_pInt
 allocate(homogenization_RGC_Ngrains(3,maxNinstance));       homogenization_RGC_Ngrains = 0_pInt
 allocate(homogenization_RGC_ciAlpha(maxNinstance));         homogenization_RGC_ciAlpha = 0.0_pReal
 allocate(homogenization_RGC_xiAlpha(maxNinstance));         homogenization_RGC_xiAlpha = 0.0_pReal
 allocate(homogenization_RGC_dAlpha(3,maxNinstance));        homogenization_RGC_dAlpha = 0.0_pReal
 allocate(homogenization_RGC_angles(3,maxNinstance));        homogenization_RGC_angles = 400.0_pReal
 allocate(homogenization_RGC_output(maxval(homogenization_Noutput),maxNinstance)); homogenization_RGC_output = ''
 allocate(homogenization_RGC_sizePostResult(maxval(homogenization_Noutput),maxNinstance))
  homogenization_RGC_sizePostResult = 0_pInt
 allocate(homogenization_RGC_orientation(3,3,mesh_maxNips,mesh_NcpElems))
 forall (i = 1_pInt:mesh_maxNips,e = 1_pInt:mesh_NcpElems)
   homogenization_RGC_orientation(:,:,i,e) = math_I3
 end forall
 
 rewind(myFile)
 line = ''
 section = 0_pInt
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= material_partHomogenization)     ! wind forward to <homogenization>
   read(myFile,'(a1024)',END=100) line
 enddo

 do                                                       ! read thru sections of phase part
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                ! next section
     section = section + 1_pInt
     output = 0_pInt                                           ! reset output counter
   endif
   if (section > 0_pInt .and. homogenization_type(section) == homogenization_RGC_label) then  ! one of my sections
     i = homogenization_typeInstance(section)             ! which instance of my type is present homogenization
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))        ! extract key
     select case(tag)
       case ('(output)')
         output = output + 1_pInt
         homogenization_RGC_output(output,i) = IO_lc(IO_stringValue(line,positions,2_pInt))
       case ('clustersize')
              homogenization_RGC_Ngrains(1,i) = IO_intValue(line,positions,2_pInt)
              homogenization_RGC_Ngrains(2,i) = IO_intValue(line,positions,3_pInt)
              homogenization_RGC_Ngrains(3,i) = IO_intValue(line,positions,4_pInt)
       case ('scalingparameter')
              homogenization_RGC_xiAlpha(i) = IO_floatValue(line,positions,2_pInt)
       case ('overproportionality')
              homogenization_RGC_ciAlpha(i) = IO_floatValue(line,positions,2_pInt)
       case ('grainsize')
              homogenization_RGC_dAlpha(1,i) = IO_floatValue(line,positions,2_pInt)
              homogenization_RGC_dAlpha(2,i) = IO_floatValue(line,positions,3_pInt)
              homogenization_RGC_dAlpha(3,i) = IO_floatValue(line,positions,4_pInt)
       case ('clusterorientation')
              homogenization_RGC_angles(1,i) = IO_floatValue(line,positions,2_pInt)
              homogenization_RGC_angles(2,i) = IO_floatValue(line,positions,3_pInt)
              homogenization_RGC_angles(3,i) = IO_floatValue(line,positions,4_pInt)
     end select
   endif
 enddo

!*** assigning cluster orientations
 do e = 1_pInt,mesh_NcpElems
   if (homogenization_type(mesh_element(3,e)) == homogenization_RGC_label) then
     myInstance = homogenization_typeInstance(mesh_element(3,e))
     if (all (homogenization_RGC_angles(:,myInstance) >= 399.9_pReal)) then
       homogenization_RGC_orientation(:,:,1,e) = math_EulerToR(math_sampleRandomOri())
       do i = 1_pInt,FE_Nips(mesh_element(2,e))
         if (microstructure_elemhomo(mesh_element(4,e))) then
           homogenization_RGC_orientation(:,:,i,e) = homogenization_RGC_orientation(:,:,1,e)
         else
           homogenization_RGC_orientation(:,:,i,e) = math_EulerToR(math_sampleRandomOri())
         endif
       enddo
     else
       do i = 1_pInt,FE_Nips(mesh_element(2,e))
         homogenization_RGC_orientation(:,:,i,e) = math_EulerToR(homogenization_RGC_angles(:,myInstance)*inRad)
       enddo
     endif
   endif
 enddo

100   if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt) then
  !$OMP CRITICAL (write2out)
   do i = 1_pInt,maxNinstance
     write(6,'(a15,1x,i4)')  'instance:  ', i
     write(6,*)
     write(6,'(a25,3(1x,i8))')    'cluster size:         ',(homogenization_RGC_Ngrains(j,i),j=1_pInt,3_pInt)
     write(6,'(a25,1x,e10.3)')    'scaling parameter:    ', homogenization_RGC_xiAlpha(i)
     write(6,'(a25,1x,e10.3)')    'over-proportionality: ', homogenization_RGC_ciAlpha(i)
     write(6,'(a25,3(1x,e10.3))') 'grain size:           ',(homogenization_RGC_dAlpha(j,i),j=1_pInt,3_pInt)
     write(6,'(a25,3(1x,e10.3))') 'cluster orientation:  ',(homogenization_RGC_angles(j,i),j=1_pInt,3_pInt)
   enddo
  !$OMP END CRITICAL (write2out)
 endif
 
 do i = 1_pInt,maxNinstance
   do j = 1_pInt,maxval(homogenization_Noutput)
     select case(homogenization_RGC_output(j,i))
       case('constitutivework')
         mySize = 1_pInt
       case('magnitudemismatch')
         mySize = 3_pInt
       case('penaltyenergy')
         mySize = 1_pInt
       case('volumediscrepancy')
         mySize = 1_pInt
       case('averagerelaxrate')
         mySize = 1_pInt
       case('maximumrelaxrate')
         mySize = 1_pInt
       case default
         mySize = 0_pInt
     end select

     if (mySize > 0_pInt) then                               ! any meaningful output found
         homogenization_RGC_sizePostResult(j,i) = mySize
         homogenization_RGC_sizePostResults(i) = &
         homogenization_RGC_sizePostResults(i) + mySize
     endif
   enddo



   homogenization_RGC_sizeState(i) &
       = 3_pInt*(homogenization_RGC_Ngrains(1,i)-1_pInt)*homogenization_RGC_Ngrains(2,i)*homogenization_RGC_Ngrains(3,i) &
         + 3_pInt*homogenization_RGC_Ngrains(1,i)*(homogenization_RGC_Ngrains(2,i)-1_pInt)*homogenization_RGC_Ngrains(3,i) &
         + 3_pInt*homogenization_RGC_Ngrains(1,i)*homogenization_RGC_Ngrains(2,i)*(homogenization_RGC_Ngrains(3,i)-1_pInt) &
         + 8_pInt   ! (1) Average constitutive work, (2-4) Overall mismatch, (5) Average penalty energy, 
                    ! (6) Volume discrepancy, (7) Avg relaxation rate component, (8) Max relaxation rate component
 enddo

endsubroutine


!*********************************************************************
!* initial homogenization state                                      *
!*********************************************************************
function homogenization_RGC_stateInit(myInstance)

 implicit none
!* Definition of variables
 integer(pInt), intent(in) :: myInstance
 real(pReal), dimension(homogenization_RGC_sizeState(myInstance)) :: homogenization_RGC_stateInit

!* Open a debugging file 
!  open(1978,file='homogenization_RGC_debugging.out',status='unknown')
 homogenization_RGC_stateInit = 0.0_pReal
 
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
 use prec,  only: p_vec
 use debug, only: debug_level, &
                  debug_homogenization, &
                  debug_levelExtensive
 use mesh,  only: mesh_element
 use material,  only: homogenization_maxNgrains,homogenization_Ngrains,homogenization_typeInstance
 use FEsolving, only: theInc,cycleCounter

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
!
 integer(pInt), parameter :: nFace = 6_pInt
 

!* Debugging the overall deformation gradient
 if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a,i3,a,i3,a)')'========== Increment: ',theInc,' Cycle: ',cycleCounter,' =========='
   write(6,'(1x,a32)')'Overall deformation gradient: '
   do i = 1_pInt,3_pInt
     write(6,'(1x,3(e15.8,1x))')(avgF(i,j), j = 1_pInt,3_pInt)
   enddo
   write(6,*)' '
   call flush(6)
   !$OMP END CRITICAL (write2out)
 endif

!* Compute the deformation gradient of individual grains due to relaxations
 homID = homogenization_typeInstance(mesh_element(3,el))
 F = 0.0_pReal
 do iGrain = 1_pInt,homogenization_Ngrains(mesh_element(3,el))
   iGrain3 = homogenization_RGC_grain1to3(iGrain,homID)
   do iFace = 1_pInt,nFace
     intFace = homogenization_RGC_getInterface(iFace,iGrain3)        ! identifying 6 interfaces of each grain
     aVect = homogenization_RGC_relaxationVector(intFace,state,homID)! get the relaxation vectors for each interface from global relaxation vector array
     nVect = homogenization_RGC_interfaceNormal(intFace,ip,el)       ! get the normal of each interface
     forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) &
     F(i,j,iGrain) = F(i,j,iGrain) + aVect(i)*nVect(j)               ! calculating deformation relaxations due to interface relaxation
   enddo
   F(:,:,iGrain) = F(:,:,iGrain) + avgF(:,:)                         ! resulting relaxed deformation gradient

!* Debugging the grain deformation gradients
   if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt) then
     !$OMP CRITICAL (write2out)
     write(6,'(1x,a32,1x,i3)')'Deformation gradient of grain: ',iGrain
     do i = 1_pInt,3_pInt
       write(6,'(1x,3(e15.8,1x))')(F(i,j,iGrain), j = 1_pInt,3_pInt)
     enddo
     write(6,*)' '
     call flush(6)
     !$OMP END CRITICAL (write2out)
   endif

 enddo

endsubroutine

!********************************************************************
! update the internal state of the homogenization scheme
! and tell whether "done" and "happy" with result
!********************************************************************
function homogenization_RGC_updateState(&
   state, &         ! my state
   state0, &        ! my state at the beginning of increment
!
   P, &             ! array of current grain stresses
   F, &             ! array of current grain deformation gradients
   F0, &            ! array of initial grain deformation gradients
   avgF, &          ! average deformation gradient
   dt, &            ! time increment
   dPdF, &          ! array of current grain stiffnesses
   ip, &            ! my integration point
   el  &            ! my element
  )

 use prec,  only: pReal,pInt,p_vec
 use debug, only: debug_level, &
                  debug_homogenization,&
                  debug_levelExtensive, &
                  debug_e, &
                  debug_i
 use math,  only: math_invert
 use mesh,  only: mesh_element
 use material, only: homogenization_maxNgrains,homogenization_typeInstance, &
                     homogenization_Ngrains
 use numerics, only: absTol_RGC,relTol_RGC,absMax_RGC,relMax_RGC,pPert_RGC, &
                     maxdRelax_RGC,viscPower_RGC,viscModus_RGC,refRelaxRate_RGC

 implicit none

!* Definition of variables
 type(p_vec), intent(inout) :: state
 type(p_vec), intent(in)    :: state0
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in)     :: P,F,F0
 real(pReal), dimension (3,3,3,3,homogenization_maxNgrains), intent(in) :: dPdF
 real(pReal), dimension (3,3), intent(in) :: avgF
 real(pReal), intent(in)                  :: dt
 integer(pInt), intent(in)                :: ip,el
! 
 logical, dimension(2)        :: homogenization_RGC_updateState
 integer(pInt), dimension (4) :: intFaceN,intFaceP,faceID
 integer(pInt), dimension (3) :: nGDim,iGr3N,iGr3P,stresLoc
 integer(pInt), dimension (2) :: residLoc
 integer(pInt) homID,iNum,i,j,nIntFaceTot,iGrN,iGrP,iMun,iFace,k,l,ipert,iGrain,nGrain
 real(pReal), dimension (3,3,homogenization_maxNgrains) :: R,pF,pR,D,pD
 real(pReal), dimension (3,homogenization_maxNgrains)   :: NN,pNN
 real(pReal), dimension (3)   :: normP,normN,mornP,mornN
 real(pReal) residMax,stresMax,constitutiveWork,penaltyEnergy,volDiscrep
 logical error
!
 integer(pInt), parameter :: nFace = 6_pInt
!
 real(pReal), dimension(:,:), allocatable :: tract,jmatrix,jnverse,smatrix,pmatrix,rmatrix
 real(pReal), dimension(:), allocatable   :: resid,relax,p_relax,p_resid,drelax

!* -------------------------------------------------------------------------------------------------------------
!*** Initialization of RGC update state calculation
!* Get the dimension of the cluster (grains and interfaces)
 homID  = homogenization_typeInstance(mesh_element(3,el))
 nGDim  = homogenization_RGC_Ngrains(:,homID)
 nGrain = homogenization_Ngrains(mesh_element(3,el))
 nIntFaceTot = (nGDim(1)-1_pInt)*nGDim(2)*nGDim(3) + nGDim(1)*(nGDim(2)-1_pInt)*nGDim(3) &
               + nGDim(1)*nGDim(2)*(nGDim(3)-1_pInt)

!* Allocate the size of the global relaxation arrays/jacobian matrices depending on the size of the cluster
 allocate(resid(3_pInt*nIntFaceTot));  resid = 0.0_pReal
 allocate(tract(nIntFaceTot,3));  tract = 0.0_pReal
 allocate(relax(3_pInt*nIntFaceTot));  relax = state%p(1:3_pInt*nIntFaceTot)
 allocate(drelax(3_pInt*nIntFaceTot))
   drelax = state%p(1:3_pInt*nIntFaceTot) - state0%p(1:3_pInt*nIntFaceTot)

!* Debugging the obtained state
 if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Obtained state: '
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,2(e15.8,1x))')state%p(i)
   enddo
   write(6,*)' '
   !$OMP END CRITICAL (write2out)
 endif

!* Computing interface mismatch and stress penalty tensor for all interfaces of all grains
 call homogenization_RGC_stressPenalty(R,NN,avgF,F,ip,el,homID)

!* Calculating volume discrepancy and stress penalty related to overall volume discrepancy 
 call homogenization_RGC_volumePenalty(D,volDiscrep,F,avgF,ip,el,homID)

!* Debugging the mismatch, stress and penalties of grains
 if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   do iGrain = 1_pInt,nGrain
     write(6,'(1x,a30,1x,i3,1x,a4,3(1x,e15.8))')'Mismatch magnitude of grain(',iGrain,') :',NN(1,iGrain),NN(2,iGrain),NN(3,iGrain)
     write(6,*)' '
     write(6,'(1x,a30,1x,i3)')'Stress and penalties of grain: ',iGrain
     do i = 1_pInt,3_pInt
       write(6,'(1x,3(e15.8,1x),1x,3(e15.8,1x),1x,3(e15.8,1x))')(P(i,j,iGrain), j = 1_pInt,3_pInt), &
                                                          (R(i,j,iGrain), j = 1_pInt,3_pInt), &
                                                          (D(i,j,iGrain), j = 1_pInt,3_pInt)
     enddo
     write(6,*)' '
   enddo
   !$OMP END CRITICAL (write2out)
 endif
!* End of initialization

!* -------------------------------------------------------------------------------------------------------------
!*** Computing the residual stress from the balance of traction at all (interior) interfaces
 do iNum = 1_pInt,nIntFaceTot
   faceID = homogenization_RGC_interface1to4(iNum,homID)                 ! identifying the interface ID in local coordinate system (4-dimensional index)
   
!* Identify the left/bottom/back grain (-|N)
   iGr3N = faceID(2:4)                                                   ! identifying the grain ID in local coordinate system (3-dimensional index)
   iGrN = homogenization_RGC_grain3to1(iGr3N,homID)                      ! translate the local grain ID into global coordinate system (1-dimensional index)
   intFaceN = homogenization_RGC_getInterface(2_pInt*faceID(1),iGr3N)
   normN = homogenization_RGC_interfaceNormal(intFaceN,ip,el)            ! get the interface normal
   
!* Identify the right/up/front grain (+|P)
   iGr3P = iGr3N
   iGr3P(faceID(1)) = iGr3N(faceID(1))+1_pInt                                 ! identifying the grain ID in local coordinate system (3-dimensional index)
   iGrP = homogenization_RGC_grain3to1(iGr3P,homID)                      ! translate the local grain ID into global coordinate system (1-dimensional index)
   intFaceP = homogenization_RGC_getInterface(2_pInt*faceID(1)-1_pInt,iGr3P)
   normP = homogenization_RGC_interfaceNormal(intFaceP,ip,el)            ! get the interface normal

!* Compute the residual of traction at the interface (in local system, 4-dimensional index)
   do i = 1_pInt,3_pInt
     tract(iNum,i) = sign(viscModus_RGC*(abs(drelax(i+3*(iNum-1_pInt)))/(refRelaxRate_RGC*dt))**viscPower_RGC, &
                          drelax(i+3*(iNum-1_pInt)))                          ! contribution from the relaxation viscosity
     do j = 1_pInt,3_pInt
       tract(iNum,i) = tract(iNum,i) + (P(i,j,iGrP) + R(i,j,iGrP) + D(i,j,iGrP))*normP(j) &
                                     + (P(i,j,iGrN) + R(i,j,iGrN) + D(i,j,iGrN))*normN(j)
                                                                         ! contribution from material stress P, mismatch penalty R, and volume penalty D
                                                                         ! projected into the interface
       resid(i+3_pInt*(iNum-1_pInt)) = tract(iNum,i)                               ! translate the local residual into global 1-dimensional residual array
     enddo
   enddo
   
!* Debugging the residual stress
   if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt) then
     !$OMP CRITICAL (write2out)
     write(6,'(1x,a30,1x,i3)')'Traction at interface: ',iNum
     write(6,'(1x,3(e15.8,1x))')(tract(iNum,j), j = 1_pInt,3_pInt)
     write(6,*)' '
     !$OMP END CRITICAL (write2out)
   endif
 enddo
!* End of residual stress calculation

!* -------------------------------------------------------------------------------------------------------------
!*** Convergence check for stress residual
 stresMax = maxval(abs(P))                                                        ! get the maximum of first Piola-Kirchhoff (material) stress
 stresLoc = int(maxloc(abs(P)),pInt)                                              ! get the location of the maximum stress
 residMax = maxval(abs(tract))                                                    ! get the maximum of the residual
 residLoc = int(maxloc(abs(tract)),pInt)                                          ! get the position of the maximum residual

 !* Debugging the convergent criteria
 if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt &
     .and. debug_e == el .and. debug_i == ip) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a)')' '
   write(6,'(1x,a,1x,i2,1x,i4)')'RGC residual check ...',ip,el
   write(6,'(1x,a15,1x,e15.8,1x,a7,i3,1x,a12,i2,i2)')'Max stress: ',stresMax, &
              '@ grain',stresLoc(3),'in component',stresLoc(1),stresLoc(2)
   write(6,'(1x,a15,1x,e15.8,1x,a7,i3,1x,a12,i2)')'Max residual: ',residMax, &
              '@ iface',residLoc(1),'in direction',residLoc(2)
   call flush(6)
   !$OMP END CRITICAL (write2out)
 endif
 
 homogenization_RGC_updateState = .false.
!* If convergence reached => done and happy
 if (residMax < relTol_RGC*stresMax .or. residMax < absTol_RGC) then 
   homogenization_RGC_updateState = .true.
   
    if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt &
        .and. debug_e == el .and. debug_i == ip) then 
     !$OMP CRITICAL (write2out)
     write(6,'(1x,a55)')'... done and happy'
     write(6,*)' '
     call flush(6)
     !$OMP END CRITICAL (write2out)
   endif
!    write(6,'(1x,a,1x,i3,1x,a6,1x,i3,1x,a12)')'RGC_updateState: ip',ip,'| el',el,'converged :)'

!* Then compute/update the state for postResult, i.e., ...
!* ... all energy densities computed by time-integration
   constitutiveWork = state%p(3*nIntFaceTot+1)
   penaltyEnergy    = state%p(3*nIntFaceTot+5)
   do iGrain = 1_pInt,homogenization_Ngrains(mesh_element(3,el))             ! time-integration loop for the calculating the work and energy
     do i = 1_pInt,3_pInt
     do j = 1_pInt,3_pInt
       constitutiveWork = constitutiveWork + P(i,j,iGrain)*(F(i,j,iGrain) - F0(i,j,iGrain))/real(nGrain,pReal)
       penaltyEnergy    = penaltyEnergy    + R(i,j,iGrain)*(F(i,j,iGrain) - F0(i,j,iGrain))/real(nGrain,pReal)
     enddo
     enddo
   enddo
   state%p(3*nIntFaceTot+1) = constitutiveWork                          ! the bulk mechanical/constitutive work
   state%p(3*nIntFaceTot+2) = sum(NN(1,:))/real(nGrain,pReal)                 ! the overall mismatch of all interface normal to e1-direction
   state%p(3*nIntFaceTot+3) = sum(NN(2,:))/real(nGrain,pReal)                 ! the overall mismatch of all interface normal to e2-direction
   state%p(3*nIntFaceTot+4) = sum(NN(3,:))/real(nGrain,pReal)                 ! the overall mismatch of all interface normal to e3-direction
   state%p(3*nIntFaceTot+5) = penaltyEnergy                             ! the overall penalty energy
   state%p(3*nIntFaceTot+6) = volDiscrep                                ! the overall volume discrepancy
   state%p(3*nIntFaceTot+7) = sum(abs(drelax))/dt/real(3_pInt*nIntFaceTot,pReal)   ! the average rate of relaxation vectors
   state%p(3*nIntFaceTot+8) = maxval(abs(drelax))/dt                    ! the maximum rate of relaxation vectors

   if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt &
        .and. debug_e == el .and. debug_i == ip) then
     !$OMP CRITICAL (write2out)
     write(6,'(1x,a30,1x,e15.8)')'Constitutive work: ',constitutiveWork
     write(6,'(1x,a30,3(1x,e15.8))')'Magnitude mismatch: ',sum(NN(1,:))/real(nGrain,pReal), &
                                                           sum(NN(2,:))/real(nGrain,pReal), &
                                                           sum(NN(3,:))/real(nGrain,pReal)
     write(6,'(1x,a30,1x,e15.8)')'Penalty energy: ',penaltyEnergy
     write(6,'(1x,a30,1x,e15.8)')'Volume discrepancy: ',volDiscrep
     write(6,*)''
     write(6,'(1x,a30,1x,e15.8)')'Maximum relaxation rate: ',maxval(abs(drelax))/dt
     write(6,'(1x,a30,1x,e15.8)')'Average relaxation rate: ',sum(abs(drelax))/dt/real(3_pInt*nIntFaceTot,pReal)
     write(6,*)''
     call flush(6)
     !$OMP END CRITICAL (write2out)
   endif
   
   deallocate(tract,resid,relax,drelax)
   return
   
!* If residual blows-up => done but unhappy
 elseif (residMax > relMax_RGC*stresMax .or. residMax > absMax_RGC) then
!* Try to restart when residual blows up exceeding maximum bound
   homogenization_RGC_updateState = (/.true.,.false./)                  ! with direct cut-back

   if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt &
       .and. debug_e == el .and. debug_i == ip) then
     !$OMP CRITICAL (write2out)
     write(6,'(1x,a55)')'... broken'
     write(6,*)' '
     call flush(6)
     !$OMP END CRITICAL (write2out)
   endif
   
   deallocate(tract,resid,relax,drelax)
   return 

!* Otherwise, proceed with computing the Jacobian and state update
 else

   if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt &
     .and. debug_e == el .and. debug_i == ip) then
     !$OMP CRITICAL (write2out)
     write(6,'(1x,a55)')'... not yet done'
     write(6,*)' '
     call flush(6)
     !$OMP END CRITICAL (write2out)
   endif
   
 endif
!*** End of convergence check for residual stress

!* -------------------------------------------------------------------------------------------------------------
!*** Construct the global Jacobian matrix for updating the global relaxation vector array when convergence is not yet reached ...
!* ... of the constitutive stress tangent, 
!* assembled from dPdF or material constitutive model "smatrix"
 allocate(smatrix(3*nIntFaceTot,3*nIntFaceTot)); smatrix = 0.0_pReal
 do iNum = 1_pInt,nIntFaceTot
   faceID = homogenization_RGC_interface1to4(iNum,homID)                 ! assembling of local dPdF into global Jacobian matrix
   
!* Identify the left/bottom/back grain (-|N)
   iGr3N = faceID(2:4)                                                   ! identifying the grain ID in local coordinate sytem
   iGrN = homogenization_RGC_grain3to1(iGr3N,homID)                      ! translate into global grain ID
   intFaceN = homogenization_RGC_getInterface(2_pInt*faceID(1),iGr3N)         ! identifying the connecting interface in local coordinate system
   normN = homogenization_RGC_interfaceNormal(intFaceN,ip,el)            ! get the interface normal
   do iFace = 1_pInt,nFace
     intFaceN = homogenization_RGC_getInterface(iFace,iGr3N)             ! identifying all interfaces that influence relaxation of the above interface
     mornN = homogenization_RGC_interfaceNormal(intFaceN,ip,el)          ! get normal of the interfaces
     iMun = homogenization_RGC_interface4to1(intFaceN,homID)             ! translate the interfaces ID into local 4-dimensional index
     if (iMun .gt. 0) then                                               ! get the corresponding tangent
       do i=1_pInt,3_pInt; do j=1_pInt,3_pInt; do k=1_pInt,3_pInt; do l=1_pInt,3_pInt
         smatrix(3*(iNum-1)+i,3*(iMun-1)+j) = smatrix(3*(iNum-1)+i,3*(iMun-1)+j) + dPdF(i,k,j,l,iGrN)*normN(k)*mornN(l)
       enddo;enddo;enddo;enddo
                                                                         ! projecting the material tangent dPdF into the interface
                                                                         ! to obtain the Jacobian matrix contribution of dPdF
     endif
   enddo
   
!* Identify the right/up/front grain (+|P)
   iGr3P = iGr3N
   iGr3P(faceID(1)) = iGr3N(faceID(1))+1_pInt                                 ! identifying the grain ID in local coordinate sytem
   iGrP = homogenization_RGC_grain3to1(iGr3P,homID)                      ! translate into global grain ID
   intFaceP = homogenization_RGC_getInterface(2_pInt*faceID(1)-1_pInt,iGr3P)       ! identifying the connecting interface in local coordinate system
   normP = homogenization_RGC_interfaceNormal(intFaceP,ip,el)            ! get the interface normal
   do iFace = 1_pInt,nFace
     intFaceP = homogenization_RGC_getInterface(iFace,iGr3P)             ! identifying all interfaces that influence relaxation of the above interface
     mornP = homogenization_RGC_interfaceNormal(intFaceP,ip,el)          ! get normal of the interfaces
     iMun = homogenization_RGC_interface4to1(intFaceP,homID)             ! translate the interfaces ID into local 4-dimensional index
     if (iMun .gt. 0) then                                               ! get the corresponding tangent
       do i=1_pInt,3_pInt; do j=1_pInt,3_pInt; do k=1_pInt,3_pInt; do l=1_pInt,3_pInt
         smatrix(3*(iNum-1)+i,3*(iMun-1)+j) = smatrix(3*(iNum-1)+i,3*(iMun-1)+j) + dPdF(i,k,j,l,iGrP)*normP(k)*mornP(l)
       enddo;enddo;enddo;enddo
     endif
   enddo
 enddo
 
!* Debugging the global Jacobian matrix of stress tangent
 if (iand(debug_level(debug_homogenization),debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Jacobian matrix of stress'
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,100(e11.4,1x))')(smatrix(i,j), j = 1_pInt,3_pInt*nIntFaceTot)
   enddo
   write(6,*)' '
   call flush(6)
   !$OMP END CRITICAL (write2out)
 endif

!* ... of the stress penalty tangent (mismatch penalty and volume penalty, 
!*     computed using numerical perturbation method) "pmatrix"
 allocate(pmatrix(3*nIntFaceTot,3*nIntFaceTot)); pmatrix = 0.0_pReal
 allocate(p_relax(3*nIntFaceTot));               p_relax = 0.0_pReal
 allocate(p_resid(3*nIntFaceTot));               p_resid = 0.0_pReal
 do ipert = 1_pInt,3_pInt*nIntFaceTot
   p_relax = relax
   p_relax(ipert) = relax(ipert) + pPert_RGC                               ! perturb the relaxation vector
   state%p(1:3*nIntFaceTot) = p_relax
   call homogenization_RGC_grainDeformation(pF,F0,avgF,state,ip,el)        ! compute the grains deformation from perturbed state
   call homogenization_RGC_stressPenalty(pR,pNN,avgF,pF,ip,el,homID)       ! compute stress penalty due to interface mismatch from perturbed state
   call homogenization_RGC_volumePenalty(pD,volDiscrep,pF,avgF,ip,el,homID)! compute stress penalty due to volume discrepancy from perturbed state

!* Computing the global stress residual array from the perturbed state
   p_resid = 0.0_pReal
   do iNum = 1_pInt,nIntFaceTot
     faceID = homogenization_RGC_interface1to4(iNum,homID)                 ! identifying the interface ID in local coordinate system (4-dimensional index)

!* Identify the left/bottom/back grain (-|N)
     iGr3N = faceID(2:4)                                                   ! identifying the grain ID in local coordinate system (3-dimensional index)
     iGrN = homogenization_RGC_grain3to1(iGr3N,homID)                      ! translate the local grain ID into global coordinate system (1-dimensional index)
     intFaceN = homogenization_RGC_getInterface(2_pInt*faceID(1),iGr3N)         ! identifying the interface ID of the grain
     normN = homogenization_RGC_interfaceNormal(intFaceN,ip,el)            ! get the corresponding interface normal
 
!* Identify the right/up/front grain (+|P)
     iGr3P = iGr3N
     iGr3P(faceID(1)) = iGr3N(faceID(1))+1_pInt                                 ! identifying the grain ID in local coordinate system (3-dimensional index)
     iGrP = homogenization_RGC_grain3to1(iGr3P,homID)                      ! translate the local grain ID into global coordinate system (1-dimensional index)
     intFaceP = homogenization_RGC_getInterface(2_pInt*faceID(1)-1_pInt,iGr3P)       ! identifying the interface ID of the grain
     normP = homogenization_RGC_interfaceNormal(intFaceP,ip,el)            ! get the corresponding normal
 
!* Compute the residual stress (contribution of mismatch and volume penalties) from perturbed state at all interfaces
     do i = 1_pInt,3_pInt
     do j = 1_pInt,3_pInt
       p_resid(i+3*(iNum-1)) = p_resid(i+3*(iNum-1)) + (pR(i,j,iGrP) - R(i,j,iGrP))*normP(j) &
                                                     + (pR(i,j,iGrN) - R(i,j,iGrN))*normN(j) &
                                                     + (pD(i,j,iGrP) - D(i,j,iGrP))*normP(j) &
                                                     + (pD(i,j,iGrN) - D(i,j,iGrN))*normN(j)
     enddo
     enddo
   enddo
   pmatrix(:,ipert) = p_resid/pPert_RGC
 enddo
 
!* Debugging the global Jacobian matrix of penalty tangent
 if (iand(debug_level(debug_homogenization), debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Jacobian matrix of penalty'
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,100(e11.4,1x))')(pmatrix(i,j), j = 1_pInt,3_pInt*nIntFaceTot)
   enddo
   write(6,*)' '
   call flush(6)
   !$OMP END CRITICAL (write2out)
 endif
 
!* ... of the numerical viscosity traction "rmatrix"
 allocate(rmatrix(3*nIntFaceTot,3*nIntFaceTot)); rmatrix = 0.0_pReal
 forall (i=1_pInt:3_pInt*nIntFaceTot) &
   rmatrix(i,i) = viscModus_RGC*viscPower_RGC/(refRelaxRate_RGC*dt)* &
                  (abs(drelax(i))/(refRelaxRate_RGC*dt))**(viscPower_RGC - 1.0_pReal)
                                                                           ! tangent due to numerical viscosity traction appears
                                                                           ! only in the main diagonal term 

!* Debugging the global Jacobian matrix of numerical viscosity tangent
 if (iand(debug_level(debug_homogenization), debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Jacobian matrix of penalty'
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,100(e11.4,1x))')(rmatrix(i,j), j = 1_pInt,3_pInt*nIntFaceTot)
   enddo
   write(6,*)' '
   call flush(6)
   !$OMP END CRITICAL (write2out)
 endif

!* The overall Jacobian matrix summarizing contributions of smatrix, pmatrix, rmatrix
 allocate(jmatrix(3*nIntFaceTot,3*nIntFaceTot)); jmatrix = smatrix + pmatrix + rmatrix
 
 if (iand(debug_level(debug_homogenization), debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Jacobian matrix (total)'
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,100(e11.4,1x))')(jmatrix(i,j), j = 1_pInt,3_pInt*nIntFaceTot)
   enddo
   write(6,*)' '
   call flush(6)
   !$OMP END CRITICAL (write2out)
 endif
 
!*** End of construction and assembly of Jacobian matrix

!* -------------------------------------------------------------------------------------------------------------
!*** Computing the update of the state variable (relaxation vectors) using the Jacobian matrix
 allocate(jnverse(3_pInt*nIntFaceTot,3_pInt*nIntFaceTot)); jnverse = 0.0_pReal
 call math_invert(size(jmatrix,1),jmatrix,jnverse,error)                              ! Compute the inverse of the overall Jacobian matrix
 
!* Debugging the inverse Jacobian matrix
 if (iand(debug_level(debug_homogenization), debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Jacobian inverse'
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,100(e11.4,1x))')(jnverse(i,j), j = 1_pInt,3_pInt*nIntFaceTot)
   enddo
   write(6,*)' '
   call flush(6)
   !$OMP END CRITICAL (write2out)
 endif

!* Calculate the state update (global relaxation vectors) for the next Newton-Raphson iteration
 drelax = 0.0_pReal
 do i = 1_pInt,3_pInt*nIntFaceTot
   do j = 1_pInt,3_pInt*nIntFaceTot
     drelax(i) = drelax(i) - jnverse(i,j)*resid(j)                     ! Calculate the correction for the state variable
   enddo
 enddo
 relax = relax + drelax                                                ! Updateing the state variable for the next iteration
 state%p(1:3*nIntFaceTot) = relax
 if (any(abs(drelax(:)) > maxdRelax_RGC)) then                         ! Forcing cutback when the incremental change of relaxation vector becomes too large
   homogenization_RGC_updateState = (/.true.,.false./)
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a,1x,i3,1x,a,1x,i3,1x,a)')'RGC_updateState: ip',ip,'| el',el,'enforces cutback'
   write(6,'(1x,a,1x,e15.8)')'due to large relaxation change =',maxval(abs(drelax))
   call flush(6)
   !$OMP END CRITICAL (write2out)
 endif

 !* Debugging the return state
 if (iand(debug_homogenization, debug_levelExtensive) > 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(1x,a30)')'Returned state: '
   do i = 1_pInt,3_pInt*nIntFaceTot
     write(6,'(1x,2(e15.8,1x))')state%p(i)
   enddo
   write(6,*)' '
   call flush(6)
   !$OMP END CRITICAL (write2out)
 endif

 deallocate(tract,resid,jmatrix,jnverse,relax,drelax,pmatrix,smatrix,p_relax,p_resid)
!*** End of calculation of state update
 
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

 use prec,  only: pReal,pInt,p_vec
 use debug, only: debug_level, &
                  debug_homogenization,&
                  debug_levelExtensive
 use mesh,  only: mesh_element
 use material, only: homogenization_maxNgrains,homogenization_Ngrains,homogenization_typeInstance
 use math, only: math_Plain3333to99
 
 implicit none
 real(pReal), dimension (3,3), intent(out) :: avgP
 real(pReal), dimension (3,3,3,3), intent(out) :: dAvgPdAvgF
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in) :: P
 real(pReal), dimension (3,3,3,3,homogenization_maxNgrains), intent(in) :: dPdF
 real(pReal), dimension (9,9) :: dPdF99
 integer(pInt), intent(in) :: ip,el
!
 integer(pInt) homID, i, j, Ngrains, iGrain

 homID = homogenization_typeInstance(mesh_element(3,el))
 Ngrains = homogenization_Ngrains(mesh_element(3,el))
 
!* Debugging the grain tangent
 if (iand(debug_level(debug_homogenization), debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   do iGrain = 1_pInt,Ngrains
     dPdF99 = math_Plain3333to99(dPdF(1:3,1:3,1:3,1:3,iGrain))
     write(6,'(1x,a30,1x,i3)')'Stress tangent of grain: ',iGrain
     do i = 1_pInt,9_pInt
       write(6,'(1x,(e15.8,1x))') (dPdF99(i,j), j = 1_pInt,9_pInt)
     enddo
     write(6,*)' '
   enddo
   call flush(6)
   !$OMP END CRITICAL (write2out)
 endif
 
!* Computing the average first Piola-Kirchhoff stress P and the average tangent dPdF
 avgP = sum(P,3)/real(Ngrains,pReal)
 dAvgPdAvgF = sum(dPdF,5)/real(Ngrains,pReal)

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
 use mesh, only: mesh_element
 use material, only: homogenization_maxNgrains, homogenization_Ngrains
 
 implicit none
 real(pReal), dimension (homogenization_maxNgrains), intent(in) :: Temperature
 integer(pInt), intent(in) :: ip,el
 real(pReal) homogenization_RGC_averageTemperature
 integer(pInt) :: Ngrains

!* Computing the average temperature
 Ngrains = homogenization_Ngrains(mesh_element(3,el))
 homogenization_RGC_averageTemperature = sum(Temperature(1:Ngrains))/real(Ngrains,pReal)

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
 use material, only: homogenization_typeInstance,homogenization_Noutput
 
 implicit none
 type(p_vec), intent(in) :: state
 integer(pInt), intent(in) :: ip,el
!
 integer(pInt) homID,o,c,nIntFaceTot
 real(pReal), dimension(homogenization_RGC_sizePostResults(homogenization_typeInstance(mesh_element(3,el)))) :: &
   homogenization_RGC_postResults

 homID = homogenization_typeInstance(mesh_element(3,el))
 nIntFaceTot=(homogenization_RGC_Ngrains(1,homID)-1_pInt)*homogenization_RGC_Ngrains(2,homID)*homogenization_RGC_Ngrains(3,homID)& 
            + homogenization_RGC_Ngrains(1,homID)*(homogenization_RGC_Ngrains(2,homID)-1_pInt)*homogenization_RGC_Ngrains(3,homID)&
            + homogenization_RGC_Ngrains(1,homID)*homogenization_RGC_Ngrains(2,homID)*(homogenization_RGC_Ngrains(3,homID)-1_pInt)

 c = 0_pInt
 homogenization_RGC_postResults = 0.0_pReal
 do o = 1_pInt,homogenization_Noutput(mesh_element(3,el))
   select case(homogenization_RGC_output(o,homID))
     case('constitutivework')
       homogenization_RGC_postResults(c+1) = state%p(3*nIntFaceTot+1)
       c = c + 1_pInt
     case('magnitudemismatch')
       homogenization_RGC_postResults(c+1) = state%p(3*nIntFaceTot+2)
       homogenization_RGC_postResults(c+2) = state%p(3*nIntFaceTot+3)
       homogenization_RGC_postResults(c+3) = state%p(3*nIntFaceTot+4)
       c = c + 3_pInt
     case('penaltyenergy')
       homogenization_RGC_postResults(c+1) = state%p(3*nIntFaceTot+5)
       c = c + 1_pInt
     case('volumediscrepancy')
       homogenization_RGC_postResults(c+1) = state%p(3*nIntFaceTot+6)
       c = c + 1_pInt
     case('averagerelaxrate')
       homogenization_RGC_postResults(c+1) = state%p(3*nIntFaceTot+7)
       c = c + 1_pInt
     case('maximumrelaxrate')
       homogenization_RGC_postResults(c+1) = state%p(3*nIntFaceTot+8)
       c = c + 1_pInt
   end select
 enddo

endfunction

!********************************************************************
! subroutine to calculate stress-like penalty due to deformation mismatch
!********************************************************************
subroutine homogenization_RGC_stressPenalty(&
   rPen, &          ! stress-like penalty
   nMis, &          ! total amount of mismatch
!
   avgF, &         ! initial effective stretch tensor
   fDef, &          ! deformation gradients
   ip, &            ! integration point
   el, &            ! element
   homID &          ! homogenization ID
  )
 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element
 use constitutive, only: constitutive_homogenizedC
 use math, only: math_civita,math_invert33
 use material, only: homogenization_maxNgrains,homogenization_Ngrains
 use numerics, only: xSmoo_RGC
 
 implicit none
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(out) :: rPen
 real(pReal), dimension (3,homogenization_maxNgrains), intent(out)   :: nMis
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in)  :: fDef
 real(pReal), dimension (3,3), intent(in)                            :: avgF
 integer(pInt), intent(in)    :: ip,el
 integer(pInt), dimension (4) :: intFace
 integer(pInt), dimension (3) :: iGrain3,iGNghb3,nGDim
 real(pReal), dimension (3,3) :: gDef,nDef
 real(pReal), dimension (3)   :: nVect,surfCorr
 real(pReal), dimension (2)   :: Gmoduli
 integer(pInt) homID,iGrain,iGNghb,iFace,i,j,k,l
 real(pReal) muGrain,muGNghb,nDefNorm,bgGrain,bgGNghb
! 
 integer(pInt), parameter :: nFace = 6_pInt
 real(pReal), parameter   :: nDefToler = 1.0e-10_pReal

 nGDim = homogenization_RGC_Ngrains(:,homID)

 rPen = 0.0_pReal
 nMis = 0.0_pReal

!* Get the correction factor the modulus of penalty stress representing the evolution of area of the interfaces due to deformations
 surfCorr = homogenization_RGC_surfaceCorrection(avgF,ip,el)

!* Debugging the surface correction factor
!  if (ip == 1 .and. el == 1) then
!    write(6,'(1x,a20,2(1x,i3))')'Correction factor: ',ip,el
!    write(6,'(1x,3(e11.4,1x))')(surfCorr(i), i = 1,3)
!  endif

!* -------------------------------------------------------------------------------------------------------------
!*** Computing the mismatch and penalty stress tensor of all grains 
 do iGrain = 1_pInt,homogenization_Ngrains(mesh_element(3,el))
   Gmoduli = homogenization_RGC_equivalentModuli(iGrain,ip,el)
   muGrain = Gmoduli(1)                                               ! collecting the equivalent shear modulus of grain
   bgGrain = Gmoduli(2)                                               ! and the lengthh of Burgers vector
   iGrain3 = homogenization_RGC_grain1to3(iGrain,homID)               ! get the grain ID in local 3-dimensional index (x,y,z)-position

!* Looping over all six interfaces of each grain
   do iFace = 1_pInt,nFace
     intFace = homogenization_RGC_getInterface(iFace,iGrain3)         ! get the 4-dimensional index of the interface in local numbering system of the grain
     nVect = homogenization_RGC_interfaceNormal(intFace,ip,el)        ! get the interface normal
     iGNghb3 = iGrain3                                                ! identify the neighboring grain across the interface
     iGNghb3(abs(intFace(1))) = iGNghb3(abs(intFace(1))) + int(real(intFace(1),pReal)/real(abs(intFace(1)),pReal),pInt)
     if (iGNghb3(1) < 1)        iGNghb3(1) = nGDim(1)                 ! with periodicity along e1 direction
     if (iGNghb3(1) > nGDim(1)) iGNghb3(1) = 1_pInt
     if (iGNghb3(2) < 1)        iGNghb3(2) = nGDim(2)                 ! with periodicity along e2 direction
     if (iGNghb3(2) > nGDim(2)) iGNghb3(2) = 1_pInt
     if (iGNghb3(3) < 1)        iGNghb3(3) = nGDim(3)                 ! with periodicity along e3 direction
     if (iGNghb3(3) > nGDim(3)) iGNghb3(3) = 1_pInt
     iGNghb  = homogenization_RGC_grain3to1(iGNghb3,homID)            ! get the ID of the neighboring grain
     Gmoduli = homogenization_RGC_equivalentModuli(iGNghb,ip,el)      ! collecting the shear modulus and Burgers vector of the neighbor
     muGNghb = Gmoduli(1)
     bgGNghb = Gmoduli(2)
     gDef = 0.5_pReal*(fDef(:,:,iGNghb) - fDef(:,:,iGrain))           ! compute the difference/jump in deformation gradeint across the neighbor

!* Compute the mismatch tensor of all interfaces
     nDefNorm = 0.0_pReal
     nDef = 0.0_pReal
     do i = 1_pInt,3_pInt
     do j = 1_pInt,3_pInt
       do k = 1_pInt,3_pInt
       do l = 1_pInt,3_pInt
         nDef(i,j) = nDef(i,j) - nVect(k)*gDef(i,l)*math_civita(j,k,l)! compute the interface mismatch tensor from the jump of deformation gradient
       enddo
       enddo
       nDefNorm = nDefNorm + nDef(i,j)*nDef(i,j)                      ! compute the norm of the mismatch tensor
     enddo
     enddo
     nDefNorm = max(nDefToler,sqrt(nDefNorm))                         ! approximation to zero mismatch if mismatch is zero (singularity)
     nMis(abs(intFace(1)),iGrain) = nMis(abs(intFace(1)),iGrain) + nDefNorm
                                                                      ! total amount of mismatch experienced by the grain (at all six interfaces)

!* Debugging the mismatch tensor
!      if (ip == 1 .and. el == 1) then
!        write(6,'(1x,a20,i2,1x,a20,1x,i3)')'Mismatch to face: ',intFace(1),'neighbor grain: ',iGNghb
!        do i = 1,3
!          write(6,'(1x,3(e11.4,1x))')(nDef(i,j), j = 1,3)
!        enddo
!        write(6,'(1x,a20,e11.4))')'with magnitude: ',nDefNorm
!      endif

!* Compute the stress penalty of all interfaces
     do i = 1_pInt,3_pInt
     do j = 1_pInt,3_pInt
       do k = 1_pInt,3_pInt
       do l = 1_pInt,3_pInt
         rPen(i,j,iGrain) = rPen(i,j,iGrain) + 0.5_pReal*(muGrain*bgGrain + muGNghb*bgGNghb)*homogenization_RGC_xiAlpha(homID) &
                                               *surfCorr(abs(intFace(1)))/homogenization_RGC_dAlpha(abs(intFace(1)),homID) &
                                               *cosh(homogenization_RGC_ciAlpha(homID)*nDefNorm) &
                                               *0.5_pReal*nVect(l)*nDef(i,k)/nDefNorm*math_civita(k,l,j) &
                                               *tanh(nDefNorm/xSmoo_RGC)
       enddo
       enddo
     enddo
     enddo
   enddo
   
!* Debugging the stress-like penalty
!    if (ip == 1 .and. el == 1) then
!      write(6,'(1x,a20,i2)')'Penalty of grain: ',iGrain
!      do i = 1,3
!        write(6,'(1x,3(e11.4,1x))')(rPen(i,j,iGrain), j = 1,3)
!      enddo
!    endif

 enddo
!*** End of mismatch and penalty stress tensor calculation

endsubroutine

!********************************************************************
! subroutine to calculate stress-like penalty due to volume discrepancy
!********************************************************************
subroutine homogenization_RGC_volumePenalty(&
   vPen, &          ! stress-like penalty due to volume
   vDiscrep, &      ! total volume discrepancy
!
   fDef, &          ! deformation gradients
   fAvg, &          ! overall deformation gradient
   ip, &            ! integration point
   el, &            ! element
   homID &          ! homogenization ID
  )
 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element
 use math, only: math_det33,math_inv33
 use material, only: homogenization_maxNgrains,homogenization_Ngrains
 use numerics, only: maxVolDiscr_RGC,volDiscrMod_RGC,volDiscrPow_RGC

 implicit none
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(out) :: vPen
 real(pReal), intent(out)                  :: vDiscrep
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in)  :: fDef
 real(pReal), dimension (3,3), intent(in)  :: fAvg
 integer(pInt), intent(in)                 :: ip,el
 real(pReal), dimension (homogenization_maxNgrains) :: gVol
 integer(pInt) homID,iGrain,nGrain
! 
 nGrain = homogenization_Ngrains(mesh_element(3,el))

!* Compute the volumes of grains and of cluster
 vDiscrep = math_det33(fAvg)                                           ! compute the volume of the cluster
 do iGrain = 1_pInt,nGrain
   gVol(iGrain) = math_det33(fDef(:,:,iGrain))                         ! compute the volume of individual grains
   vDiscrep     = vDiscrep - gVol(iGrain)/real(nGrain,pReal)                  ! calculate the difference/dicrepancy between
                                                                        ! the volume of the cluster and the the total volume of grains
 enddo

!* Calculate the stress and penalty due to volume discrepancy
 vPen      = 0.0_pReal
 do iGrain = 1_pInt,nGrain
   vPen(:,:,iGrain) = -1.0_pReal/real(nGrain,pReal)*volDiscrMod_RGC*volDiscrPow_RGC/maxVolDiscr_RGC* &
                      sign((abs(vDiscrep)/maxVolDiscr_RGC)**(volDiscrPow_RGC - 1.0),vDiscrep)* &
                      gVol(iGrain)*transpose(math_inv33(fDef(:,:,iGrain)))

!* Debugging the stress-like penalty of volume discrepancy
!    if (ip == 1 .and. el == 1) then
!      write(6,'(1x,a30,i2)')'Volume penalty of grain: ',iGrain
!      do i = 1,3
!        write(6,'(1x,3(e11.4,1x))')(vPen(i,j,iGrain), j = 1,3)
!      enddo
!    endif

 enddo

endsubroutine

!********************************************************************
! subroutine to compute the correction factor due to surface area evolution
!********************************************************************
function homogenization_RGC_surfaceCorrection(&
   avgF, &        ! average deformation gradient
   ip, &          ! my IP
   el &           ! my element
  )

 use prec, only: pReal,pInt,p_vec
 use math, only: math_invert33,math_mul33x33
 
 implicit none
 real(pReal), dimension(3,3), intent(in) :: avgF
 real(pReal), dimension(3)               :: homogenization_RGC_surfaceCorrection
 integer(pInt), intent(in)               :: ip,el
 real(pReal), dimension(3,3)             :: invC,avgC
 real(pReal), dimension(3)               :: nVect
 real(pReal)   detF
 integer(pInt), dimension(4)             :: intFace
 integer(pInt) i,j,iBase
 logical       error

!* Compute the correction factor accouted for surface evolution (area change) due to deformation
 avgC = 0.0_pReal
 avgC = math_mul33x33(transpose(avgF),avgF)
 invC = 0.0_pReal
 call math_invert33(avgC,invC,detF,error)
 homogenization_RGC_surfaceCorrection = 0.0_pReal
 do iBase = 1_pInt,3_pInt
   intFace = (/iBase,1_pInt,1_pInt,1_pInt/)
   nVect = homogenization_RGC_interfaceNormal(intFace,ip,el)             ! get the normal of the interface
   do i = 1_pInt,3_pInt
   do j = 1_pInt,3_pInt
     homogenization_RGC_surfaceCorrection(iBase) = &                     ! compute the component of (the inverse of) the stretch in the direction of the normal
       homogenization_RGC_surfaceCorrection(iBase) + invC(i,j)*nVect(i)*nVect(j)
   enddo
   enddo
   homogenization_RGC_surfaceCorrection(iBase) = &                       ! get the surface correction factor (area contraction/enlargement)
     sqrt(homogenization_RGC_surfaceCorrection(iBase))*detF
 enddo

endfunction

!********************************************************************
! subroutine to compute the equivalent shear and bulk moduli from the elasticity tensor
!********************************************************************
function homogenization_RGC_equivalentModuli(&
   grainID, &        ! grain ID
   ip, &             ! IP number
   el &              ! element number
  )

 use prec, only: pReal,pInt,p_vec
 use constitutive, only: constitutive_homogenizedC,constitutive_averageBurgers

 implicit none
 integer(pInt), intent(in)    :: grainID,ip,el
 real(pReal), dimension (6,6) :: elasTens
 real(pReal), dimension(2)    :: homogenization_RGC_equivalentModuli
 real(pReal) cEquiv_11,cEquiv_12,cEquiv_44

 elasTens = constitutive_homogenizedC(grainID,ip,el)

!* Compute the equivalent shear modulus after Turterltaub and Suiker, JMPS (2005)
 cEquiv_11 = (elasTens(1,1) + elasTens(2,2) + elasTens(3,3))/3.0_pReal
 cEquiv_12 = (elasTens(1,2) + elasTens(2,3) + elasTens(3,1) + &
              elasTens(1,3) + elasTens(2,1) + elasTens(3,2))/6.0_pReal
 cEquiv_44 = (elasTens(4,4) + elasTens(5,5) + elasTens(6,6))/3.0_pReal
 homogenization_RGC_equivalentModuli(1) = 0.2_pReal*(cEquiv_11 - cEquiv_12) + 0.6_pReal*cEquiv_44

!* Obtain the length of Burgers vector
 homogenization_RGC_equivalentModuli(2) = constitutive_averageBurgers(grainID,ip,el)

endfunction

!********************************************************************
! subroutine to collect relaxation vectors of an interface
!********************************************************************
function homogenization_RGC_relaxationVector(&
   intFace, &        ! set of interface ID in 4D array (normal and position)
   state, &          ! set of global relaxation vectors
   homID &           ! homogenization ID
  )

 use prec, only: pReal,pInt,p_vec
 
 implicit none
 real(pReal), dimension (3)               :: homogenization_RGC_relaxationVector
 integer(pInt), dimension (4), intent(in) :: intFace
 type(p_vec), intent(in)      :: state
 integer(pInt), dimension (3) :: nGDim
 integer(pInt) iNum,homID

!* Collect the interface relaxation vector from the global state array
 homogenization_RGC_relaxationVector = 0.0_pReal
 nGDim = homogenization_RGC_Ngrains(:,homID)
 iNum = homogenization_RGC_interface4to1(intFace,homID)                  ! identify the position of the interface in global state array
 if (iNum .gt. 0_pInt) homogenization_RGC_relaxationVector = state%p((3*iNum-2):(3*iNum))    
                                                                         ! get the corresponding entries

endfunction

!********************************************************************
! subroutine to identify the normal of an interface
!********************************************************************
function homogenization_RGC_interfaceNormal(&
   intFace, &        ! interface ID in 4D array (normal and position)
   ip, &             ! my IP
   el &              ! my element
  )

 use prec, only: pReal,pInt,p_vec
 use math, only: math_mul33x3
 
 implicit none
 real(pReal), dimension (3)               :: homogenization_RGC_interfaceNormal
 integer(pInt), dimension (4), intent(in) :: intFace
 integer(pInt), intent(in)                :: ip,el
 integer(pInt) nPos

!* Get the normal of the interface, identified from the value of intFace(1)
 homogenization_RGC_interfaceNormal = 0.0_pReal
 nPos = abs(intFace(1))                                                                             ! identify the position of the interface in global state array
 homogenization_RGC_interfaceNormal(nPos) = real(intFace(1)/abs(intFace(1)),pReal)                  ! get the normal vector w.r.t. cluster axis

 homogenization_RGC_interfaceNormal = &
   math_mul33x3(homogenization_RGC_orientation(:,:,ip,el),homogenization_RGC_interfaceNormal)
                                                                                                    ! map the normal vector into sample coordinate system (basis)

!  if (ip == 1 .and. el == 1) then
!    write(6,'(1x,a32,3(1x,i3))')'Interface normal: ',intFace(1)
!    write(6,'(1x,3(e15.8,1x))')(nVect(i), i = 1,3)
!    write(6,*)' '
!    call flush(6)
!  endif

endfunction

!********************************************************************
! subroutine to collect six faces of a grain in 4D (normal and position)
!********************************************************************
function homogenization_RGC_getInterface(&
   iFace, &         ! face index (1..6) mapped like (-e1,-e2,-e3,+e1,+e2,+e3) or iDir = (-1,-2,-3,1,2,3)
   iGrain3 &        ! grain ID in 3D array
  )
 use prec, only: pReal,pInt,p_vec
 
 implicit none
 integer(pInt), dimension (4) :: homogenization_RGC_getInterface
 integer(pInt), dimension (3), intent(in)    :: iGrain3
 integer(pInt), intent(in) :: iFace
 integer(pInt) iDir
 
!* Direction of interface normal
 iDir = (int(real(iFace-1_pInt,pReal)/2.0_pReal,pInt)+1_pInt)*(-1_pInt)**iFace
 homogenization_RGC_getInterface(1) = iDir
 
!* Identify the interface position by the direction of its normal
 homogenization_RGC_getInterface(2:4) = iGrain3(:)
 if (iDir < 0_pInt) &                                                   ! to have a correlation with coordinate/position in real space
   homogenization_RGC_getInterface(1_pInt-iDir) = homogenization_RGC_getInterface(1_pInt-iDir)-1_pInt

endfunction

!********************************************************************
! subroutine to map grain ID from in 1D (global array) to in 3D (local position)
!********************************************************************
function homogenization_RGC_grain1to3(&
   grain1, &        ! grain ID in 1D array
   homID &          ! homogenization ID
  )

 use prec, only: pInt,p_vec
 
 implicit none
 integer(pInt), dimension (3) :: homogenization_RGC_grain1to3
 integer(pInt), intent(in)    :: grain1,homID
 integer(pInt), dimension (3) :: nGDim

!* Get the grain position
 nGDim = homogenization_RGC_Ngrains(:,homID)
 homogenization_RGC_grain1to3(3) = 1_pInt+(grain1-1_pInt)/(nGDim(1)*nGDim(2))
 homogenization_RGC_grain1to3(2) = 1_pInt+mod((grain1-1_pInt)/nGDim(1),nGDim(2))
 homogenization_RGC_grain1to3(1) = 1_pInt+mod((grain1-1_pInt),nGDim(1))

endfunction

!********************************************************************
! subroutine to map grain ID from in 3D (local position) to in 1D (global array)
!********************************************************************
function homogenization_RGC_grain3to1(&
   grain3, &       ! grain ID in 3D array (pos.x,pos.y,pos.z)
   homID &         ! homogenization ID
  )

 use prec, only: pInt,p_vec

 implicit none
 integer(pInt), dimension (3), intent(in) :: grain3
 integer(pInt)                :: homogenization_RGC_grain3to1
 integer(pInt), dimension (3) :: nGDim
 integer(pInt) homID

!* Get the grain ID
 nGDim = homogenization_RGC_Ngrains(:,homID)
 homogenization_RGC_grain3to1 = grain3(1) + nGDim(1)*(grain3(2)-1_pInt) + nGDim(1)*nGDim(2)*(grain3(3)-1_pInt)

endfunction

!********************************************************************
! subroutine to map interface ID from 4D (normal and local position) into 1D (global array)
!********************************************************************
function homogenization_RGC_interface4to1(&
   iFace4D, &        ! interface ID in 4D array (n.dir,pos.x,pos.y,pos.z)
   homID &           ! homogenization ID
  )

 use prec, only: pInt,p_vec
 
 implicit none
 integer(pInt), dimension (4), intent(in) :: iFace4D
 integer(pInt)                :: homogenization_RGC_interface4to1
 integer(pInt), dimension (3) :: nGDim,nIntFace
 integer(pInt) homID

 nGDim = homogenization_RGC_Ngrains(:,homID)
!* Compute the total number of interfaces, which ...
 nIntFace(1) = (nGDim(1)-1_pInt)*nGDim(2)*nGDim(3)                           ! ... normal //e1
 nIntFace(2) = nGDim(1)*(nGDim(2)-1_pInt)*nGDim(3)                           ! ... normal //e2
 nIntFace(3) = nGDim(1)*nGDim(2)*(nGDim(3)-1_pInt)                           ! ... normal //e3

!* Get the corresponding interface ID in 1D global array
 if (abs(iFace4D(1)) == 1_pInt) then                                    ! ... interface with normal //e1
   homogenization_RGC_interface4to1 = iFace4D(3) + nGDim(2)*(iFace4D(4)-1_pInt) &
                                      + nGDim(2)*nGDim(3)*(iFace4D(2)-1_pInt)
   if ((iFace4D(2) == 0_pInt) .or. (iFace4D(2) == nGDim(1))) homogenization_RGC_interface4to1 = 0_pInt
 elseif (abs(iFace4D(1)) == 2_pInt) then                                ! ... interface with normal //e2
   homogenization_RGC_interface4to1 = iFace4D(4) + nGDim(3)*(iFace4D(2)-1_pInt) &
                                      + nGDim(3)*nGDim(1)*(iFace4D(3)-1_pInt) + nIntFace(1)
   if ((iFace4D(3) == 0_pInt) .or. (iFace4D(3) == nGDim(2))) homogenization_RGC_interface4to1 = 0_pInt
 elseif (abs(iFace4D(1)) == 3_pInt) then                                ! ... interface with normal //e3
   homogenization_RGC_interface4to1 = iFace4D(2) + nGDim(1)*(iFace4D(3)-1_pInt) &
                                      + nGDim(1)*nGDim(2)*(iFace4D(4)-1_pInt) + nIntFace(1) + nIntFace(2)
   if ((iFace4D(4) == 0_pInt) .or. (iFace4D(4) == nGDim(3))) homogenization_RGC_interface4to1 = 0_pInt
 endif

endfunction

!********************************************************************
! subroutine to map interface ID from 1D (global array) into 4D (normal and local position)
!********************************************************************
function homogenization_RGC_interface1to4(&
   iFace1D, &        ! interface ID in 1D array
   homID &           ! homogenization ID
  )

 use prec, only: pReal,pInt,p_vec
 
 implicit none
 integer(pInt), dimension (4) :: homogenization_RGC_interface1to4
 integer(pInt), intent(in)    :: iFace1D
 integer(pInt), dimension (3) :: nGDim,nIntFace
 integer(pInt) homID

 nGDim = homogenization_RGC_Ngrains(:,homID)
!* Compute the total number of interfaces, which ...
 nIntFace(1) = (nGDim(1)-1_pInt)*nGDim(2)*nGDim(3)      ! ... normal //e1
 nIntFace(2) = nGDim(1)*(nGDim(2)-1_pInt)*nGDim(3)      ! ... normal //e2
 nIntFace(3) = nGDim(1)*nGDim(2)*(nGDim(3)-1_pInt)      ! ... normal //e3

!* Get the corresponding interface ID in 4D (normal and local position)
 if (iFace1D > 0 .and. iFace1D <= nIntFace(1)) then                              ! ... interface with normal //e1
   homogenization_RGC_interface1to4(1) = 1_pInt
   homogenization_RGC_interface1to4(3) = mod((iFace1D-1_pInt),nGDim(2))+1_pInt
   homogenization_RGC_interface1to4(4) = mod(&
                                             int(&
                                                 real(iFace1D-1_pInt,pReal)/&
                                                 real(nGDim(2),pReal)&
                                                 ,pInt)&
                                             ,nGDim(3))+1_pInt
   homogenization_RGC_interface1to4(2) = int(&
                                             real(iFace1D-1_pInt,pReal)/&
                                             real(nGDim(2),pReal)/&
                                             real(nGDim(3),pReal)&
                                             ,pInt)+1_pInt
 elseif (iFace1D > nIntFace(1) .and. iFace1D <= (nIntFace(2) + nIntFace(1))) then  ! ... interface with normal //e2
   homogenization_RGC_interface1to4(1) = 2_pInt
   homogenization_RGC_interface1to4(4) = mod((iFace1D-nIntFace(1)-1_pInt),nGDim(3))+1_pInt
   homogenization_RGC_interface1to4(2) = mod(&
                                             int(&
                                                 real(iFace1D-nIntFace(1)-1_pInt,pReal)/&
                                                 real(nGDim(3),pReal)&
                                                 ,pInt)&
                                              ,nGDim(1))+1_pInt
   homogenization_RGC_interface1to4(3) = int(&
                                             real(iFace1D-nIntFace(1)-1_pInt,pReal)/&
                                             real(nGDim(3),pReal)/&
                                             real(nGDim(1),pReal)&
                                             ,pInt)+1_pInt
 elseif (iFace1D > nIntFace(2) + nIntFace(1) .and. iFace1D <= (nIntFace(3) + nIntFace(2) + nIntFace(1))) then ! ... interface with normal //e3
   homogenization_RGC_interface1to4(1) = 3_pInt
   homogenization_RGC_interface1to4(2) = mod((iFace1D-nIntFace(2)-nIntFace(1)-1_pInt),nGDim(1))+1_pInt
   homogenization_RGC_interface1to4(3) = mod(&
                                             int(&
                                                 real(iFace1D-nIntFace(2)-nIntFace(1)-1_pInt,pReal)/&
                                                 real(nGDim(1),pReal)&
                                                 ,pInt)&
                                             ,nGDim(2))+1_pInt
   homogenization_RGC_interface1to4(4) = int(&
                                             real(iFace1D-nIntFace(2)-nIntFace(1)-1_pInt,pReal)/&
                                             real(nGDim(1),pReal)/&
                                             real(nGDim(2),pReal)&
                                             ,pInt)+1_pInt
 endif

endfunction

!********************************************************************
! calculating the grain deformation gradient
! (the same with homogenization_RGC_partionDeformation,
!  but used only for perturbation scheme)
!********************************************************************
subroutine homogenization_RGC_grainDeformation(&
   F, &             ! partioned def grad per grain
!
   F0, &            ! initial partioned def grad per grain
   avgF, &          ! my average def grad
   state, &         ! my state
   ip, &            ! my IP
   el  &            ! my element
  )
 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element
 use material, only: homogenization_maxNgrains,homogenization_Ngrains,homogenization_typeInstance
 
 implicit none
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(out) :: F
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in)  :: F0
 real(pReal), dimension (3,3), intent(in) :: avgF
 type(p_vec), intent(in) :: state
 integer(pInt), intent(in) :: el,ip
!
 real(pReal), dimension (3)   :: aVect,nVect
 integer(pInt), dimension (4) :: intFace
 integer(pInt), dimension (3) :: iGrain3
 integer(pInt) homID, iGrain,iFace,i,j
!
 integer(pInt), parameter :: nFace = 6_pInt

!* Compute the deformation gradient of individual grains due to relaxations
 homID = homogenization_typeInstance(mesh_element(3,el))
 F = 0.0_pReal
 do iGrain = 1_pInt,homogenization_Ngrains(mesh_element(3,el))
   iGrain3 = homogenization_RGC_grain1to3(iGrain,homID)
   do iFace = 1_pInt,nFace
     intFace = homogenization_RGC_getInterface(iFace,iGrain3)
     aVect = homogenization_RGC_relaxationVector(intFace,state,homID)
     nVect = homogenization_RGC_interfaceNormal(intFace,ip,el)
     forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) &
     F(i,j,iGrain) = F(i,j,iGrain) + aVect(i)*nVect(j)               ! effective relaxations
   enddo
   F(:,:,iGrain) = F(:,:,iGrain) + avgF(:,:)                         ! relaxed deformation gradient
 enddo

endsubroutine

END MODULE
