!************************************
!*      Module: CONSTITUTIVE        *
!************************************
!* contains:                        *
!* - constitutive equations         *
!* - Schmid matrices calculation    *
!* - Hardening matrices definition  *
!* - Parameters definition          *
!* - orientations?                  *
!************************************

MODULE constitutive

!*** Include other modules ***
use prec, only: pReal,pInt
! NB: 'only'-commend may not be needed
implicit none

!*****************************
!*     Material parameters   *
!*****************************
!* Character *
character*80, allocatble :: TCfile(:),ODFfile(:)
! NB: orientation files  TCfile(number of material)

!* Integer *
integer(pInt) constitutive_Nmats
! NB: Number of materials (read in material file)
integer(pInt), allocatable :: constitutive_crystal_structure(:)
! NB: crystal_structure(number of material)=1-3
integer(pInt) constitutive_Nslip(3)
! NB: Number of systems for each crystal structure (3)
! NB: not forget a MaxSlip variable that give the effective number
!     of slip system I have  

!* Real *
real(pReal), allocatable :: constitutive_Cslip_66(:,:,:)
! NB: Cslip_66(1:6,1:6,number of materials)
real(pReal), allocatable :: constitutive_s0_slip(:)
real(pReal), allocatable :: constitutive_gdot0_slip(:)
real(pReal), allocatable :: constitutive_n_slip(:)
real(pReal), allocatable :: constitutive_h0(:)
real(pReal), allocatable :: constitutive_w0(:)
real(pReal), allocatable :: constitutive_s_sat(:)
! NB: Parameters(number of materials)
real(pReal), allocatable :: constitutive_hardening_matrix(:,:,:)
! NB: hardening_matrix(48,48,3)
real(pReal), parameter :: constitutive_latent_hardening=1.4_pReal
real(pReal) constitutive_sn(3,48,3),constitutive_sd(3,48,3)
! NB: slip normale and slip direction for 3 crystal structures
!     Is 48 always the maximum number of systems?
real(pReal) constitutive_Sslip(3,3,48,3),constitutive_Sslip_v(6,48,3)
! NB: Schmid matrices and corresponding Schmid vectors

!*** Slip systems for FCC structures (1) ***
constitutive_Nslip(1)=12_pInt
!* System {111}<110>  Sort according Eisenlohr&Hantcherli
data constitutive_sd(:, 1,1)/ 0, 1,-1/ ; data constitutive_sn(:, 1,1)/ 1, 1, 1/
data constitutive_sd(:, 2,1)/-1, 0, 1/ ; data constitutive_sn(:, 2,1)/ 1, 1, 1/
data constitutive_sd(:, 3,1)/ 1,-1, 0/ ; data constitutive_sn(:, 3,1)/ 1, 1, 1/
data constitutive_sd(:, 4,1)/ 0,-1,-1/ ; data constitutive_sn(:, 4,1)/-1,-1, 1/
data constitutive_sd(:, 5,1)/ 1, 0, 1/ ; data constitutive_sn(:, 5,1)/-1,-1, 1/
data constitutive_sd(:, 6,1)/-1, 1, 0/ ; data constitutive_sn(:, 6,1)/-1,-1, 1/
data constitutive_sd(:, 7,1)/ 0,-1, 1/ ; data constitutive_sn(:, 7,1)/ 1,-1,-1/
data constitutive_sd(:, 8,1)/-1, 0,-1/ ; data constitutive_sn(:, 8,1)/ 1,-1,-1/
data constitutive_sd(:, 9,1)/ 1, 1, 0/ ; data constitutive_sn(:, 9,1)/ 1,-1,-1/
data constitutive_sd(:,10,1)/ 0, 1, 1/ ; data constitutive_sn(:,10,1)/-1, 1,-1/
data constitutive_sd(:,11,1)/ 1, 0,-1/ ; data constitutive_sn(:,11,1)/-1, 1,-1/
data constitutive_sd(:,12,1)/-1,-1, 0/ ; data constitutive_sn(:,12,1)/-1, 1,-1/

!*** Slip systems for BCC structures (2) ***
constitutive_Nslip(2)=48_pInt
!* System {110}<111>
!* Sort?
data constitutive_sd(:, 1,2)/ 1,-1, 1/ ; data constitutive_sn(:, 1,2)/ 0, 1, 1/
data constitutive_sd(:, 2,2)/-1,-1, 1/ ; data constitutive_sn(:, 2,2)/ 0, 1, 1/
data constitutive_sd(:, 3,2)/ 1, 1, 1/ ; data constitutive_sn(:, 3,2)/ 0,-1, 1/
data constitutive_sd(:, 4,2)/-1, 1, 1/ ; data constitutive_sn(:, 4,2)/ 0,-1, 1/
data constitutive_sd(:, 5,2)/-1, 1, 1/ ; data constitutive_sn(:, 5,2)/ 1, 0, 1/
data constitutive_sd(:, 6,2)/-1,-1, 1/ ; data constitutive_sn(:, 6,2)/ 1, 0, 1/
data constitutive_sd(:, 7,2)/ 1, 1, 1/ ; data constitutive_sn(:, 7,2)/-1, 0, 1/
data constitutive_sd(:, 8,2)/ 1,-1, 1/ ; data constitutive_sn(:, 8,2)/-1, 0, 1/
data constitutive_sd(:, 9,2)/-1, 1, 1/ ; data constitutive_sn(:, 9,2)/ 1, 1, 0/
data constitutive_sd(:,10,2)/-1, 1,-1/ ; data constitutive_sn(:,10,2)/ 1, 1, 0/
data constitutive_sd(:,11,2)/ 1, 1, 1/ ; data constitutive_sn(:,11,2)/-1, 1, 0/
data constitutive_sd(:,12,2)/ 1, 1,-1/ ; data constitutive_sn(:,12,2)/-1, 1, 0/
!* System {112}<111>
!* Sort?
data constitutive_sd(:,13,2)/-1, 1, 1/ ; data constitutive_sn(:,13,2)/ 2, 1, 1/
data constitutive_sd(:,14,2)/ 1, 1, 1/ ; data constitutive_sn(:,14,2)/-2, 1, 1/
data constitutive_sd(:,15,2)/ 1, 1,-1/ ; data constitutive_sn(:,15,2)/ 2,-1, 1/
data constitutive_sd(:,16,2)/ 1,-1, 1/ ; data constitutive_sn(:,16,2)/ 2, 1,-1/
data constitutive_sd(:,17,2)/ 1,-1, 1/ ; data constitutive_sn(:,17,2)/ 1, 2, 1/
data constitutive_sd(:,18,2)/ 1, 1,-1/ ; data constitutive_sn(:,18,2)/-1, 2, 1/
data constitutive_sd(:,19,2)/ 1, 1, 1/ ; data constitutive_sn(:,19,2)/ 1,-2, 1/
data constitutive_sd(:,20,2)/-1, 1, 1/ ; data constitutive_sn(:,20,2)/ 1, 2,-1/
data constitutive_sd(:,21,2)/ 1, 1,-1/ ; data constitutive_sn(:,21,2)/ 1, 1, 2/
data constitutive_sd(:,22,2)/ 1,-1, 1/ ; data constitutive_sn(:,22,2)/-1, 1, 2/
data constitutive_sd(:,23,2)/-1, 1, 1/ ; data constitutive_sn(:,23,2)/ 1,-1, 2/
data constitutive_sd(:,24,2)/ 1, 1, 1/ ; data constitutive_sn(:,24,2)/ 1, 1,-2/
!* System {123}<111>
!* Sort?
data constitutive_sd(:,25,2)/ 1, 1,-1/ ; data constitutive_sn(:,25,2)/ 1, 2, 3/
data constitutive_sd(:,26,2)/ 1,-1, 1/ ; data constitutive_sn(:,26,2)/-1, 2, 3/
data constitutive_sd(:,27,2)/-1, 1, 1/ ; data constitutive_sn(:,27,2)/ 1,-2, 3/
data constitutive_sd(:,28,2)/ 1, 1, 1/ ; data constitutive_sn(:,28,2)/ 1, 2,-3/
data constitutive_sd(:,29,2)/ 1,-1, 1/ ; data constitutive_sn(:,29,2)/ 1, 3, 2/
data constitutive_sd(:,30,2)/ 1, 1,-1/ ; data constitutive_sn(:,30,2)/-1, 3, 2/
data constitutive_sd(:,31,2)/ 1, 1, 1/ ; data constitutive_sn(:,31,2)/ 1,-3, 2/
data constitutive_sd(:,32,2)/-1, 1, 1/ ; data constitutive_sn(:,32,2)/ 1, 3,-2/
data constitutive_sd(:,33,2)/ 1, 1,-1/ ; data constitutive_sn(:,33,2)/ 2, 1, 3/
data constitutive_sd(:,34,2)/ 1,-1, 1/ ; data constitutive_sn(:,34,2)/-2, 1, 3/
data constitutive_sd(:,35,2)/-1, 1, 1/ ; data constitutive_sn(:,35,2)/ 2,-1, 3/
data constitutive_sd(:,36,2)/ 1, 1, 1/ ; data constitutive_sn(:,36,2)/ 2, 1,-3/
data constitutive_sd(:,37,2)/ 1,-1, 1/ ; data constitutive_sn(:,37,2)/ 2, 3, 1/
data constitutive_sd(:,38,2)/ 1, 1,-1/ ; data constitutive_sn(:,38,2)/-2, 3, 1/
data constitutive_sd(:,39,2)/ 1, 1, 1/ ; data constitutive_sn(:,39,2)/ 2,-3, 1/
data constitutive_sd(:,40,2)/-1, 1, 1/ ; data constitutive_sn(:,40,2)/ 2, 3,-1/
data constitutive_sd(:,41,2)/-1, 1, 1/ ; data constitutive_sn(:,41,2)/ 3, 1, 2/
data constitutive_sd(:,42,2)/ 1, 1, 1/ ; data constitutive_sn(:,42,2)/-3, 1, 2/
data constitutive_sd(:,43,2)/ 1, 1,-1/ ; data constitutive_sn(:,43,2)/ 3,-1, 2/
data constitutive_sd(:,44,2)/ 1,-1, 1/ ; data constitutive_sn(:,44,2)/ 3, 1,-2/
data constitutive_sd(:,45,2)/-1, 1, 1/ ; data constitutive_sn(:,45,2)/ 3, 2, 1/
data constitutive_sd(:,46,2)/ 1, 1, 1/ ; data constitutive_sn(:,46,2)/-3, 2, 1/
data constitutive_sd(:,47,2)/ 1, 1,-1/ ; data constitutive_sn(:,47,2)/ 3,-2, 1/
data constitutive_sd(:,48,2)/ 1,-1, 1/ ; data constitutive_sn(:,48,2)/ 3, 2,-1/

!*** Slip systems for HCP structures (3) ***
constitutive_Nslip(3)=12_pInt
!* Basal systems {0001}<1120>
!* 1- (0 0 0 1)[-2  1  1  0]
!* 2- (0 0 0 1)[ 1 -2  1  0]
!* 3- (0 0 0 1)[ 1  1 -2  0]
!* Plane (hkil)->(hkl)
!* Direction [uvtw]->[(u-t) (v-t) w]
!* Automatical transformation from Bravais to Miller
!* not done for the moment
!* Sort?
data constitutive_sd(:, 1,3)/-1, 0, 0/ ; data constitutive_sn(:, 1,3)/ 0, 0, 1/
data constitutive_sd(:, 2,3)/ 0,-1, 0/ ; data constitutive_sn(:, 2,3)/ 0, 0, 1/
data constitutive_sd(:, 3,3)/ 1, 1, 0/ ; data constitutive_sn(:, 3,3)/ 0, 0, 1/
!* 1st type prismatic systems {1010}<1120> 
!* 1- ( 0  1 -1  0)[-2  1  1  0]
!* 2- ( 1  0 -1  0)[ 1 -2  1  0]
!* 3- (-1  1  0  0)[ 1  1 -2  0]
!* Sort?
data constitutive_sd(:, 4,3)/-1, 0, 0/ ; data constitutive_sn(:, 4,3)/ 0, 1, 0/
data constitutive_sd(:, 5,3)/ 0,-1, 0/ ; data constitutive_sn(:, 5,3)/ 1, 0, 0/
data constitutive_sd(:, 6,3)/ 1, 1, 0/ ; data constitutive_sn(:, 6,3)/-1, 1, 0/
!* 1st type 1st order pyramidal systems {1011}<1120>
!* 1- ( 0 -1  1  1)[-2  1  1  0]
!* 2- ( 0  1 -1  1)[-2  1  1  0]
!* 3- (-1  0  1  1)[ 1 -2  1  0]
!* 4- ( 1  0 -1  1)[ 1 -2  1  0]
!* 5- (-1  1  0  1)[ 1  1 -2  0]
!* 6- ( 1 -1  0  1)[ 1  1 -2  0]
!* Sort?
data constitutive_sd(:, 7,3)/-1, 0, 0/ ; data constitutive_sn(:, 7,3)/ 0,-1, 1/
data constitutive_sd(:, 8,3)/ 0,-1, 0/ ; data constitutive_sn(:, 8,3)/ 0, 1, 1/
data constitutive_sd(:, 9,3)/ 1, 1, 0/ ; data constitutive_sn(:, 9,3)/-1, 0, 1/
data constitutive_sd(:,10,3)/-1, 0, 0/ ; data constitutive_sn(:,10,3)/ 1, 0, 1/
data constitutive_sd(:,11,3)/ 0,-1, 0/ ; data constitutive_sn(:,11,3)/-1, 1, 1/
data constitutive_sd(:,12,3)/ 1, 1, 0/ ; data constitutive_sn(:,12,3)/ 1,-1, 1/


CONTAINS
!****************************************
!* - constitutive_init                  *
!* - constitutive_calc_SchmidM          *
!* - constitutive_calc_HardeningM       *
!* - constitutive_parse_materialDat     *
!* - orientation reading????            *
!* - constitutive_calc_SlipRates        *
!* - constitutive_calc_Hardening        *
!* - consistutive_calc_PlasVeloGradient *
!* - CPFEM_CauchyStress???????          *
!****************************************


subroutine constitutive_init()
!**************************************
!***      Module initialization     ***
!**************************************
call constitutive_calc_SchmidM()
call constitutive_calc_hardeningM()
call constitutive_parse_materialDat()
end subroutine
 

subroutine constitutive_calc_SchmidM()
!**************************************
!*** Calculation of Schmid matrices ***
!**************************************
use prec, only: pReal,pInt
implicit none

!* Definition of variables
integer(pInt) i,j,k,l
real(pReal) invNorm

!* Iteration over the crystal structures 
do l=1,3  
!* Iteration over the systems
   do k=1,constitutive_Nslip(l)
!* Defintion of Schmid matrix  
      forall (i=1:3,j=1:3) 
	         constitutive_Sslip(i,j,k,l)=constitutive_sd(i,k,l)*constitutive_sn(j,k,l)
      endforall
!* Normalization of Schmid matrix
      invNorm=dsqrt(1.0_pReal/((constitutive_sn(1,k,l)**2+constitutive_sn(2,k,l)**2+constitutive_sn(3,k,l)**2)*(constitutive_sd(1,k,l)**2+constitutive_sd(2,k,l)**2+constitutive_sd(3,k,l)**2)))
      constitutive_Sslip(:,:,k,l)=constitutive_Sslip(:,:,k,l)*invNorm
!* Vectorization of normalized Schmid matrix
!* according MARC component order 11,22,33,12,23,13
      constitutive_Sslip_v(1,k,l)=constitutive_Sslip(1,1,k,l)
      constitutive_Sslip_v(2,k,l)=constitutive_Sslip(2,2,k,l)
      constitutive_Sslip_v(3,k,l)=constitutive_Sslip(3,3,k,l)
      constitutive_Sslip_v(4,k,l)=constitutive_Sslip(1,2,k,l)+constitutive_Sslip(2,1,k,l)
      constitutive_Sslip_v(5,k,l)=constitutive_Sslip(2,3,k,l)+constitutive_Sslip(3,3,k,l)
      constitutive_Sslip_v(6,k,l)=constitutive_Sslip(1,3,k,l)+constitutive_Sslip(3,1,k,l)
   enddo
enddo

end subroutine


subroutine constitutive_calc_HardeningM()
!****************************************
!*** Hardening matrix (see Kalidindi) ***
!****************************************
use prec, only: pReal,pInt
implicit none

!* Definition of variables
integer(pInt) i,j,k,l

!* Initialization of the hardening matrix
constitutive_hardening_matrix=constitutive_latent_hardening
!* Iteration over the crystal structures 
do l=1,3
   select case(l) 
!* Hardening matrix for FCC structures    
   case (1)
   do k=1,10,3
      forall (i=1:3,j=1:3)
             constitutive_hardening_matrix(k-1+i,k-1+j,l)=1.0_pReal
      endforall
   enddo
!* Hardening matrix for BCC structures
   case (2)
   do k=1,11,2
      forall (i=1:2,j=1:2)
             constitutive_hardening_matrix(k-1+i,k-1+j,l)=1.0_pReal
      endforall
   enddo
   do k=13,48
      constitutive_hardening_matrix(k,k,l)=1.0_pReal
   enddo
!* Hardening matrix for HCP structures
   case (3)
   forall (i=1:3,j=1:3)
          constitutive_hardening_matrix(i,j,l)=1.0_pReal
   endforall
   do k=4,12
      constitutive_hardening_matrix(k,k,l)=1.0_ZdRe
   enddo
   end select
enddo

end subroutine






!* NOT YET IMPLEMENTED *!
subroutine constitutive_parse_materialDat()
!****************************************
!***      Reading parameter files     ***
!****************************************
use prec, only: pReal,pInt
implicit none

!* Definition of variables
character*80 line
integer(pIn) i,j,k,l,positions(4)

! MISSING: needs to be 2 pass
! first pass to count Nmats and allocate
! 2nd pass to read actual parameters

 write(6,*) '## constitutive_parse_materialDat ##'
 write(6,*)

 constitutive_Nmats = 1
 open(200,FILE='material.mpie',ACTION='READ',STATUS='OLD',ERR=100)
 read(200,610,ERR=200,END=200) line

 IF( line(1:1).ne.'[' )THEN
   WRITE(6,*) 'Problem with mat file: no mat. in 1st line'
 ELSE
   WRITE(6,*) 'Reading mat. data'
   DO WHILE( .true. )
     READ(200,610,END=220) line
     IF( line(1:1).eq.'[' )THEN
  constitutive_Nmats = constitutive_Nmats+1
     ELSE
  positions = IO_stringPos(line,2) ! parse 2 parts
  SELECT CASE (IO_stringValue(line,positions,1))
  CASE ('s0_slip')
    s0_slip(mat) = IO_floatValue(line,positions,2)
  CASE ('g0_slip')
    g0_slip(mat) = IO_floatValue(line,positions,2)
  CASE ('n_slip')
    n_slip(mat) = IO_intValue(line,positions,2)
  CASE ('h0')
    h0(mat) = IO_floatValue(line,positions,2)
  CASE ('w0')
    w0(mat) = IO_floatValue(line,positions,2)
  CASE ('tauc_sat')
    tauc_sat(mat) = IO_floatValue(line,positions,2)
  CASE ('C11')
    C11(mat) = IO_floatValue(line,positions,2)
  CASE ('C12')
    C12(mat) = IO_floatValue(line,positions,2)
  CASE ('C44')
    C44(mat) = IO_floatValue(line,positions,2)
  CASE ('TCfile')
    TCfile(mat) = IO_stringValue(line,positions,2)
  CASE ('ODFfile')
    ODFfile(mat) = IO_stringValue(line,positions,2)
  CASE ('Ngrains')
    Ngrains(mat) = IO_intValue(line,positions,2)

  CASE DEFAULT
    WRITE(6,*) 'Unknown mat. parameter ',line
     END IF
   END DO
 END IF

  220 continue
 close(200)

! ** Defintion of stiffness matrices **
! MISSING: this needs to be iterated over the materials 
 Cslip_66 = 0.0_pRe
 do i=1,3
   do j=1,3
     Cslip_66(i,j)   = C12
   enddo
   Cslip_66(i,i)   = C11
   Cslip_66(i+3,i+3) = C44
 enddo

 Cslip_3333(:,:,:,:) = math_66to3333(Cslip_66(:,:))   

! *** Transformation to get the MARC order ***
! ***    11,22,33,12,23,13      ***
! MISSING this should be outsourced to FEM-spec

 temp=Cslip_66(4,:)
 Cslip_66(4,:)=Cslip_66(6,:)    
 Cslip_66(6,:)=Cslip_66(5,:)    
 Cslip_66(5,:)=temp
 temp=Cslip_66(:,4)
 Cslip_66(:,4)=2.0d0*Cslip_66(:,6)    
 Cslip_66(:,6)=2.0d0*Cslip_66(:,5)    
 Cslip_66(:,5)=2.0d0*temp 


!    *** Output to MARC output file ***  
 write(6,*) 'Material data:'
 write(6,*) 'Slip parameter:(s0_slip,g0_slip,n_slip)'
 write(6,*) s0_slip,g0_slip,n_slip
 write(6,*) 'Slip hardening parameter:(h0,tauc_sat,w0)'
 write(6,*) h0,tauc_sat,w0
 write(6,*) 'Elasticity matrix:'
 write(6,*) Cslip_66(1,:)
 write(6,*) Cslip_66(2,:)
 write(6,*) Cslip_66(3,:)
 write(6,*) Cslip_66(4,:)/2.0d0
 write(6,*) Cslip_66(5,:)/2.0d0
 write(6,*) Cslip_66(6,:)/2.0d0 
 write(6,*)
 call flush(6)

! END OF MISSING mat iterations

 return
 100 call _error(110)
 200 call _error(210)
 end


!* NOT YET IMPLEMENTED *!
 subroutine READ_ORIENTATIONS
!***********************************************************************
!*** This routine reads orientations from 'orientations.mpie'    ***
!***********************************************************************
 use mpie
 use Zahlendarstellung, only: ZdRe,ZdIn
 implicit none

!    *** Definition of variables ***
 integer(ZdIn) i,j

!    *** Read 'orientations.mpie' file ***
 open(100,FILE='orientations.mpie',ACTION='READ',STATUS='OLD',
     &     ERR=100)
 read(100,*,ERR=200,END=200)

!    *** Read number of states, maximum of components over the states ***
 read(100,*,ERR=200,END=200) mpie_nmat,mpie_norimx
!    *** Allocate memory for the arrays ***
 allocate(mpie_mat(mpie_nmat,2+7*mpie_norimx))
 allocate(mpie_cko(mpie_nmat,4:35,3,0:35,2))
 allocate(mpie_ckofile(mpie_nmat,80))
 allocate(mpie_odfmax(mpie_nmat))
 mpie_mat=0.0_ZdRe
 mpie_cko=0.0_ZdRe
 mpie_ckofile=''
 mpie_odfmax=0.0_ZdRe

!    *** Read the different states ***
 do i=1,mpie_nmat
    read(100,*,ERR=200,END=200)
!    *** Number of component and symmetry ***
    read(100,*,ERR=200,END=200) mpie_mat(i,1),mpie_mat(i,2)
!    *** If symmetry = 2, use direct ODF sampling,i.e. read coefficience ***
    if (mpie_mat(i,2)==2_ZdIn) then
  read(100,'(80A)',ERR=200,END=201) mpie_ckofile(i,:)
 201  call mpie_read_ckofile(mpie_cko(i,:,:,:,:),
     &         mpie_ckofile(i,:))
  call mpie_odf_max(mpie_cko(i,:,:,:,:),mpie_odfmax(i))
!    *** Set volume fraction to inverse of orientation number for each orientation ***
  do j=1,int(mpie_mat(i,1),ZdIn)
     mpie_mat(i,2+7*j)=1/mpie_mat(i,1)
  enddo
    else
!    *** Read for every component:          ***
!    *** gauss: euler angles (phi1, PHI, phi2), dummy, scatter, volume fraction ***
!    *** fiber: alpha1, alpha2, beta1, beta2, scatter, volume fraction    ***
  do j=1,int(mpie_mat(i,1),ZdIn)
     read(100,*,ERR=200,END=200) mpie_mat(i,7*j-4),
     &    mpie_mat(i,7*j-3),mpie_mat(i,7*j-2),
     &    mpie_mat(i,7*j-1),mpie_mat(i,7*j),
     &    mpie_mat(i,7*j+1),mpie_mat(i,7*j+2)
  enddo
     endif
 enddo
 close(100)

!    *** Output to MARC output file ***
 write(6,*) 'MPIE Material Routine Ver. 0.1 by L. Hantcherli'
 write(6,*)
 write(6,*) 'Orientations data:'
 write(6,*) 'Number of materials:  ', mpie_nmat
 write(6,*) 'Maximum number of components: ', mpie_norimx
 write(6,*)
 do i=1,mpie_nmat
    write(6,*) 'State', i
    if (mpie_mat(i,2)==2_ZdIn) then
  write(6,*) mpie_ckofile(i,:),mpie_mat(i,9),mpie_odfmax(i)
    else
  write(6,*) mpie_mat(i,:)
    endif
    write(6,*)
 enddo
 call flush(6)
 return
 100 call _error(100)
 200 call _error(200)
 end





subroutine constitutive_calc_SlipRates(matID,tau_slip,tauc_slip,gdot_slip,dgdot_dtaucslip)
!*********************************************************************
!* This subroutine contains the constitutive equation for the slip   *
!* rate on each slip system                                          *
!* INPUT:                                                            *
!*  - matID           : material identifier                          *
!*  - tau_slip        : applied shear stress on each slip system     *
!*  - tauc_slip       : critical shear stress on each slip system    *
!* OUTPUT:                                                           *
!*  - gdot_slip       : slip rate on each slip system                *
!*  - dgdot_dtaucslip : derivative of slip rate on each slip system  *
!*********************************************************************
use prec, only: pReal,pInt
implicit none
 
!* Definition of variables
integer(pInt) matID,i
real(pReal) tau_slip(constitutive_Nslip(constitutive_crystal_structure(matID)))
real(pReal) tauc_slip_new(constitutive_Nslip(constitutive_crystal_structure(matID)))
real(pReal) gdot_slip(constitutive_Nslip(constitutive_crystal_structure(matID)))
real(pReal) dgdot_dtaucslip(constitutive_Nslip(constitutive_crystal_structure(matID)))

!* Iteration over the systems 
do i=1,constitutive_Nslip(constitutive_crystal_structure(matID))
   gdot_slip(i)=constitutive_gdot0_slip(matID)*(abs(tau_slip(i))/tauc_slip(i))**constitutive_n_slip(matID)*sign(1.0_pReal,tau_slip(i))
   dgdot_dtaucslip(i)=constitutive_gdot0_slip(matID)*(abs(tau_slip(i))/tauc_slip(i))**(constitutive_n_slip(matID)-1.0_pReal)*constitutive_n_slip(matID)/tauc_slip(i)
enddo

return
end subroutine


subroutine constitutive_calc_Hardening(matID,tauc_slip,gdot_slip,dtauc_slip)
!*********************************************************************
!* This subroutine calculates the increment in critical shear stress *
!* due to plastic deformation on each slip system                    *
!* INPUT:                                                            *
!*  - matID      : material identifier                               *
!*  - tauc_slip  : critical shear stress on each slip system         *
!*  - gdot_slip  : slip rate on each slip system                     *
!* OUTPUT:                                                           *
!*  - dtauc_slip : increment of hardening due to slip on each system *
!*********************************************************************
use prec, only: pReal,pInt
implicit none

!* Definition of variables
integer(pInt) matID,i,j
real(pReal) tauc_slip(constitutive_Nslip(constitutive_crystal_structure(matID)))
real(pReal) gdot_slip(constitutive_Nslip(constitutive_crystal_structure(matID)))
real(pReal) dtauc_slip(constitutive_Nslip(constitutive_crystal_structure(matID)))
real(pReal) self_hardening(constitutive_Nslip(constitutive_crystal_structure(matID)))

!* Self-Hardening of each system
do i=1,constitutive_Nslip(constitutive_crystal_structure(matID))
   self_hardening(i)=constitutive_h0(matID)*(1.0_pReal-tauc_slip(i)/constitutive_s_sat(matID))**constitutive_w0(matID)*abs(gdot_slip(i))
enddo

!* Hardening for all systems
i=Nslip(crystal_structure(matID))
j=crystal_structure(matID)
dtauc_slip=matmul(constitutive_hardening_matrix(1:i,1:i,j),self_hardening)

return
end subroutine


subroutine constitutive_calc_PlasVeloGradient(dt,tau_slip,tauc_slip_new,Lp)
!*********************************************************************
!* This subroutine calculates the plastic velocity gradient given    *
!* the slip rates                                                    *
!* INPUT:                                                            *
!*  - matID      : material identifier                               *
!*  - dt         : time step                                         *
!*  - tau_slip   : applied shear stress on each slip system          *
!*  - tauc_slip  : critical shear stress on each slip system         *
!* OUTPUT:                                                           *
!*  - Lp         : plastic velocity gradient                         *
!*********************************************************************
use prec, only: pReal,pInt
implicit none

!* Definition of variables
integer(pInt) matID,i
real(pReal) dt,Lp(3,3)
real(pReal) tau_slip(constitutive_Nslip(constitutive_crystal_structure(matID)))
real(pReal) tauc_slip_new(constitutive_Nslip(constitutive_crystal_structure(matID)))
real(pReal) gdot_slip(constitutive_Nslip(constitutive_crystal_structure(matID)))
 
!* Calculation of Lp
Lp=0.0_pReal
do i=1,constitutive_Nslip(constitutive_crystal_structure(matID))
   gdot_slip(i)=constitutive_gdot0_slip(matID)*(abs(tau_slip(i))/tauc_slip(i))**constitutive_n_slip(matID)*sign(1.0_pReal,tau_slip(i))
   Lp=Lp+gdot_slip(i)*constitutive_Sslip(:,:,i,constitutive_crystal_structure(matID))
enddo

return
end subroutine


!function CPFEM_Cauchy(Estar_v,Fe,C66)
! ***************************************************************
! Subroutine calculates the cauchy from the elastic strain tensor
! Input: Estar_v : elastic strain tensor (in vector form)
!   Fe    : elastic deformation gradient
!   C66    : Stiffness Tensor
! Output: cs    : cauchy stress
! Local: Tstar_v,Tstar,mm,det
! ***************************************************************
!use math
!use prec
!implicit none

!real(pRe) Estar_v(6),Fe(3,3),C66(6,6),CPFEM_Cauchy(6)
!real(pRe) det,mm(3,3),Tstar(3,3)
!integer(pIn) i

!det = math_det(Fe)
!Tstar = math_6to33(matmul(C66,Estar_v))
!mm=matmul(matmul(Fe,Tstar),transpose(Fe))/det
!CPFEM_Cauchy = math_33to6(mm)

!return
!end function
 
END MODULE
