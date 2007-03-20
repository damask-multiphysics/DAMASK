
!    ---------------------------
 MODULE constitutive 
!    ---------------------------
!    *** constitutive equations ***

 use prec, only: pRe,pIn
 implicit none

! ***************************
! *** Material parameters ***
! ***************************  
 real(pRe), allocatable :: Cslip_66(:,:,:),Cslip_3333(:,:,:,:,:)
 real(pRe), allocatable :: s0_slip(:),gdot0_slip(:)
 real(pRe), allocatable :: h0(:),w0(:),s_sat(:),q0(:),n_slip(:)
 real(pRe), allocatable :: hardening_matrix(:,:,:)
 character*80, allocatable :: TCfile(:), ODFfile(:)
 real(pRe), parameter :: latent=1.4_pRe
 integer(pIn), parameter :: Nslip(3)
 integer(pIn) Nmats

 real(pRe) sn(3,48,3),sd(3,48,3)
 real(pRe) Sslip(3,48,3,3),Sslip_v(3,48,6)
!    *** Vectors n and d for each fcc slip systems ***
! MISSING needs to be generalized to fcc and bcc (and hcp?)
! 1: fcc, 2: bcc, 3: hcp
! the respective crystal structure has to be defined
! via material parameter 'crystal_structure' in [material]
 data sd( 1,:)/ 0, 1,-1/ ; data sn( 1,:)/ 1, 1, 1/
 data sd( 2,:)/-1, 0, 1/ ; data sn( 2,:)/ 1, 1, 1/
 data sd( 3,:)/ 1,-1, 0/ ; data sn( 3,:)/ 1, 1, 1/
 data sd( 4,:)/ 0,-1,-1/ ; data sn( 4,:)/-1,-1, 1/
 data sd( 5,:)/ 1, 0, 1/ ; data sn( 5,:)/-1,-1, 1/
 data sd( 6,:)/-1, 1, 0/ ; data sn( 6,:)/-1,-1, 1/
 data sd( 7,:)/ 0,-1, 1/ ; data sn( 7,:)/ 1,-1,-1/
 data sd( 8,:)/-1, 0,-1/ ; data sn( 8,:)/ 1,-1,-1/
 data sd( 9,:)/ 1, 1, 0/ ; data sn( 9,:)/ 1,-1,-1/
 data sd(10,:)/ 0, 1, 1/ ; data sn(10,:)/-1, 1,-1/
 data sd(11,:)/ 1, 0,-1/ ; data sn(11,:)/-1, 1,-1/
 data sd(12,:)/-1,-1, 0/ ; data sn(12,:)/-1, 1,-1/

 contains
 
! **************************************
! ***     module Init      ***
! **************************************
 subroutine constitutive_init()
 
 call constitutive_calc_SchmidM()
 call constitutive_calc_hardeningM()
 call constitutive_parse_materialDat()
 
 end subroutine
 
! **************************************
! *** Calculation of Schmid matrices ***
! **************************************
 subroutine constitutive_calc_SchmidM()

 use prec, only: pRe,pIn
 implicit none

 integer(pIn) i,j,k,l
 real(pRe) invNorm
 
 do j=1,3  ! iterate over crystal system
   do i=1,Nslip(j)  ! iterate over slip systems
     do k=1,3
  do l=1,3
    Sslip(j,i,k,l)=sd(j,i,k)*sn(j,i,l)
  enddo
     enddo
     invNorm = dsqrt(1.0_pRe/
     &    (sn(j,i,1)**2+sn(j,1,2)**2+sn(j,i,3)**2)/
     &    (sd(j,i,1)**2+sd(j,1,2)**2+sd(j,i,3)**2))
     Sslip(j,i,:,:) = Sslip(j,i,:,:)*invNorm
     Sslip_v(j,i,1)=Sslip(j,i,1,1)
     Sslip_v(j,i,2)=Sslip(j,i,2,2)
     Sslip_v(j,i,3)=Sslip(j,i,3,3)
     Sslip_v(j,i,4)=Sslip(j,i,1,2)+Sslip(j,i,2,1)
     Sslip_v(j,i,5)=Sslip(j,i,2,3)+Sslip(j,i,3,2)
     Sslip_v(j,i,6)=Sslip(j,i,1,3)+Sslip(j,i,3,1)
   enddo
 enddo
 end subroutine
 
! ****************************************
! *** Hardening matrix (see Kalidindi) ***
! ****************************************
 subroutine constitutive_calc_hardeningM()

 use prec, only: pRe,pIn
 implicit none

 integer(pIn) i,j,k,l

! MISSING iteration over crystal systems
! PE does not understand the j,k looping

 hardening_matrix=latent
 do i=1,10,3
   do j=1,3
     do k=1,3
  hardening_matrix(i-1+j,i-1+k)=1.0_ZdRe
     enddo
   enddo
 enddo

! ****************************************
! ***    Reading 'material.mpie'  ***     
! ****************************************
 subroutine constitutive_parse_materialDat()

 use prec, only: pRe,pIn
 implicit none

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



 subroutine slip_rate (tau_slip,tauc_slip_new,gdot_slip,
     &       dgdot_dtaucslip)
C ********************************************************************
C Subroutine contains the constitutive equation for the slip
C    rate on each slip system
C Input: tau_slip      : shear stress on each slip system 
C   tauc_slip_new   : critical shear stress on each slip system
C Output: gdot_slip      : slip rate on each slip system
C    dgdot_dtaucslip: derivative of slip rate on each slip system
C ********************************************************************
 use mpie
 use Zahlendarstellung
 implicit none

 real(ZdRe) tau_slip(nslip),tauc_slip_new(nslip),
     &      gdot_slip(nslip),dgdot_dtaucslip(nslip)
 integer(ZdIn) i
 
 do i=1,nslip
   gdot_slip(i)=g0_slip*(abs(tau_slip(i))/tauc_slip_new(i))
     &    **n_slip*sign(1.0_ZdRe,tau_slip(i))
   dgdot_dtaucslip(i)=g0_slip*(abs(tau_slip(i))/tauc_slip_new(i))
     &     **(n_slip-1) *n_slip/tauc_slip_new(i)
 enddo

 return
 end

 subroutine hardening (tauc_slip_new,gdot_slip,dtauc_slip)
C *********************************************************************
C Subroutine calculates the increment in critical shear stress due 
C   to plastic deformation on each slip system
C Input: tauc_slip_new :critical shear stress needed for slip on each 
C        slip system
C   gdot_slip    :slip rate on each slip system
C Output: dtauc_slip   :increment of hardening due to slip on each 
C        slip system
C Local : selfhr
C *********************************************************************
 use mpie
 use Zahlendarstellung
 implicit none

 real(ZdRe) tauc_slip_new(nslip),gdot_slip(nslip),
     &      dtauc_slip(nslip)
 real(ZdRe) selfhr(nslip)
 integer(ZdIn) i

 do i=1,nslip
   selfhr(i)=h0*(1.0_ZdRe-tauc_slip_new(i)/
     &  tauc_sat)**w0
     &  *abs(gdot_slip(i))
 enddo
 dtauc_slip=matmul(hardening_matrix,selfhr)

 return
 end


 subroutine plastic_vel_grad(dt,tau_slip,tauc_slip_new,Lp)
C *************************************************************
C Subroutine calculates the plastic velocity gradient given the
C   slip rates
C Input: dt     : time step
C   tau_slip    : shear stress on each slip system on each 
C         slip system
C   tauc_slip_new : critical shear stress needed for slip on each 
C         slip system
C Output: Lp     : plastic velocity gradient
C    gdot_slip    : slip rate on each slip system
C *************************************************************
 use mpie
 use Zahlendarstellung
 implicit none

 real(ZdRe) dt,tau_slip(nslip),tauc_slip_new(nslip),
     &      Lp(3,3),gdot_slip(nslip)
 integer(ZdIn) i

 Lp=0
 do i=1,nslip
   gdot_slip(i)=g0_slip*(abs(tau_slip(i))/tauc_slip_new(i))
     &     **n_slip*sign(1.0_ZdRe,tau_slip(i))
   Lp=Lp+gdot_slip(i)*Sslip(i,:,:)
 enddo

 return
 end


 function CPFEM_Cauchy(Estar_v,Fe,C66)
C ***************************************************************
C Subroutine calculates the cauchy from the elastic strain tensor
C Input: Estar_v : elastic strain tensor (in vector form)
C   Fe    : elastic deformation gradient
C   C66    : Stiffness Tensor
C Output: cs    : cauchy stress
C Local: Tstar_v,Tstar,mm,det
C ***************************************************************
 use math
 use prec
 implicit none

 real(pRe) Estar_v(6),Fe(3,3),C66(6,6),CPFEM_Cauchy(6)
 real(pRe) det,mm(3,3),Tstar(3,3)
 integer(pIn) i

 det = math_det(Fe)
 Tstar = math_6to33(matmul(C66,Estar_v))
 mm=matmul(matmul(Fe,Tstar),transpose(Fe))/det
 CPFEM_Cauchy = math_33to6(mm)

 return
 end function
 
 end module
