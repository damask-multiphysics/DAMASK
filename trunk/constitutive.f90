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
implicit none

! MISSING this should be outsourced to FEM-spec
! *** Transformation to get the MARC order ***
! ***    11,22,33,12,23,13      ***
!temp=Cslip_66(4,:)
!Cslip_66(4,:)=Cslip_66(6,:)    
!Cslip_66(6,:)=Cslip_66(5,:)    
!Cslip_66(5,:)=temp
!temp=Cslip_66(:,4)
!Cslip_66(:,4)=2.0d0*Cslip_66(:,6)    
!Cslip_66(:,6)=2.0d0*Cslip_66(:,5)    
!Cslip_66(:,5)=2.0d0*temp 


!***********************************************
!* Definition of crystal structures properties *
!***********************************************
!* Number of crystal structures (1-FCC,2-BCC,3-HCP)
integer(pInt), parameter :: constitutive_MaxCrystalStructure = 3
!* Total number of slip systems per crystal structure
!* (as to be changed according the definition of slip systems)
integer(pInt), dimension(constitutive_MaxCrystalStructure), parameter :: constitutive_MaxNslipOfStructure = &
reshape((/12,48,12/),(/constitutive_MaxCrystalStructure/))
!* Maximum number of slip systems over crystal structures
integer(pInt), parameter :: constitutive_MaxMaxNslipOfStructure = 48 
!* Slip direction, slip normales and Schmid matrices
real(pReal), dimension(3,3,constitutive_MaxMaxNslipOfStructure,constitutive_MaxCrystalStructure) :: constitutive_Sslip
real(pReal), dimension(6,constitutive_MaxMaxNslipOfStructure,constitutive_MaxCrystalStructure) :: constitutive_Sslip_v
real(pReal), dimension(3,constitutive_MaxMaxNslipOfStructure,constitutive_MaxCrystalStructure) :: constitutive_sn
real(pReal), dimension(3,constitutive_MaxMaxNslipOfStructure,constitutive_MaxCrystalStructure) :: constitutive_sd

!*** Slip systems for FCC structures (1) ***
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
!* Basal systems {0001}<1120> (independent of c/a-ratio)
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
!* 1st type prismatic systems {1010}<1120>  (independent of c/a-ratio)
!* 1- ( 0  1 -1  0)[-2  1  1  0]
!* 2- ( 1  0 -1  0)[ 1 -2  1  0]
!* 3- (-1  1  0  0)[ 1  1 -2  0]
!* Sort?
data constitutive_sd(:, 4,3)/-1, 0, 0/ ; data constitutive_sn(:, 4,3)/ 0, 1, 0/
data constitutive_sd(:, 5,3)/ 0,-1, 0/ ; data constitutive_sn(:, 5,3)/ 1, 0, 0/
data constitutive_sd(:, 6,3)/ 1, 1, 0/ ; data constitutive_sn(:, 6,3)/-1, 1, 0/
!* 1st type 1st order pyramidal systems {1011}<1120> 
!* plane normales depend on the c/a-ratio
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

!* Slip-slip interactions matrices
!* (defined for the moment as crystal structure property and not as material property)
!* (may be changed in the future)
real(pReal), dimension(constitutive_MaxMaxNslipOfStructure,constitutive_MaxMaxNslipOfStructure,constitutive_MaxCrystalStructure) :: constitutive_hardening_matrix
real(pReal), parameter :: constitutive_latent_hardening=1.4

!*************************************
!* Definition of material properties *
!*************************************
!* Number of materials
integer(pInt) material_maxN
!* Crystal structure and number of selected slip systems per material
integer(pInt), dimension(:)  , allocatable :: material_CrystalStructure
integer(pInt), dimension(:)  , allocatable :: material_Nslip
!* Maximum number of selected slip systems over materials
integer(pInt) material_MaxNslip
!* Elastic constants and matrices
real(pReal), dimension(:)    , allocatable :: material_C11
real(pReal), dimension(:)    , allocatable :: material_C12
real(pReal), dimension(:)    , allocatable :: material_C13
real(pReal), dimension(:)    , allocatable :: material_C33
real(pReal), dimension(:)    , allocatable :: material_C44
real(pReal), dimension(:,:,:), allocatable :: material_Cslip_66
!* Visco-plastic material parameters
real(pReal), dimension(:)    , allocatable :: material_s0_slip
real(pReal), dimension(:)    , allocatable :: material_gdot0_slip
real(pReal), dimension(:)    , allocatable :: material_n_slip
real(pReal), dimension(:)    , allocatable :: material_h0
real(pReal), dimension(:)    , allocatable :: material_s_sat
real(pReal), dimension(:)    , allocatable :: material_w0

!************************************
!* Definition of texture properties *
!************************************
!* Number of textures, maximum number of Gauss and Fiber components
integer(pInt) texture_maxN
integer(pInt) texture_maxNGauss
integer(pInt) texture_maxNFiber
!* Textures definition
character(len=80), dimension(:), allocatable :: texture_ODFfile
character(len=80), dimension(:), allocatable :: texture_symmetry
integer(pInt), dimension(:)    , allocatable :: texture_Ngrains
integer(pInt), dimension(:)    , allocatable :: texture_NGauss
integer(pInt),dimension(:)     , allocatable :: texture_NFiber
real(pReal), dimension(:,:,:)  , allocatable :: texture_Gauss
real(pReal), dimension(:,:,:)  , allocatable :: texture_Fiber

!************************************
!*         State variables          *
!************************************
!* integer(pInt) constitutive_maxNstate???
real(pReal), dimension(:,:,:,:), allocatable :: constitutive_state_old
real(pReal), dimension(:,:,:,:), allocatable :: constitutive_state_new

!************************************
!*             Results              *
!************************************
integer(pInt) constitutive_MaxNresults
integer(pInt), dimension(:,:,:), allocatable :: constitutive_Nresults
real(pReal), dimension(:,:,:,:), allocatable :: constitutive_results
!* IS MISSING : allocation

!************************************
!*             Other                *
!************************************
integer(pInt) constitutive_maxNgrains
integer(pInt), dimension(:,:)    , allocatable :: constitutive_Ngrains
integer(pInt), dimension(:,:,:)  , allocatable :: constitutive_matID
real(pReal), dimension(:,:,:)    , allocatable :: constitutive_matVolFrac
integer(pInt), dimension(:,:,:)  , allocatable :: constitutive_texID
real(pReal), dimension(:,:,:)    , allocatable :: constitutive_texVolFrac


CONTAINS
!****************************************
!* - constitutive_Init                  
!* - constitutive_SchmidMatrices          
!* - constitutive_HardeningMatrices
!* - constitutive_CountSections 
!* - constitutive_Parse_UnknownPart 
!* - constitutive_Parse_MaterialPart
!* - constitutive_Parse_TexturePart     
!* - constitutive_Parse_MatTexDat
!* - constitutive_Assignement                          
!* - constitutive_LpAndItsTangent     
!* - consistutive_DotState          
!****************************************


subroutine constitutive_Init()
!**************************************
!*      Module initialization         *
!**************************************
call constitutive_SchmidMatrices()
call constitutive_HardeningMatrices()
call constitutive_Parse_MatTexDat('mattex.mpie')
call constitutive_Assignement()
end subroutine
 

subroutine constitutive_SchmidMatrices()
!**************************************
!*   Calculation of Schmid matrices   *
!**************************************
use prec, only: pReal,pInt
implicit none

!* Definition of variables
integer(pInt) i,j,k,l
real(pReal) invNorm

!* Iteration over the crystal structures 
do l=1,3  
!* Iteration over the systems
   do k=1,constitutive_MaxNslipOfStructure(l)
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


subroutine constitutive_HardeningMatrices()
!****************************************
!* Hardening matrix (see Kalidindi)     *
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
      constitutive_hardening_matrix(k,k,l)=1.0_pReal
   enddo
   end select
enddo

end subroutine


subroutine constitutive_CountSections(file,count,part)
!*********************************************************************
!* This subroutine reads a "part" from the input file until the next *
!* part is reached and counts the number of "sections" in the part   *
!* INPUT:                                                            *
!*  - file  : file ID                                                *
!* OUTPUT:                                                           *
!*  - part  : name of the next "part"                                *
!*  - count : number of sections inside the current "part"           *
!*********************************************************************
use prec, only: pInt
use IO,   only: IO_stringPos,IO_stringValue,IO_lc
implicit none

!* Definition of variables
character(len=80) part,line,tag
integer(pInt) file,count,pos
integer(pInt), dimension(3) :: positions

count=0
part=''

do while(.true.)
   read(file,'(a80)',END=100) line
   positions=IO_stringPos(line,1)
   tag=IO_lc(IO_stringValue(line,positions,1))
   if (tag(1:1)=='<'.AND.tag(len_trim(tag):len_trim(tag))=='>') then  
      part=tag(2:len_trim(tag)-1)
	  exit
   elseif (tag(1:1)=='[') then
      count=count+1
   endif
enddo
100 return

end subroutine


subroutine constitutive_CountGaussAndFiber(file,count,part)
!*********************************************************************
!*********************************************************************
use prec, only: pInt
use IO,   only: IO_stringPos,IO_stringValue,IO_lc
implicit none

!* Definition of variables
character(len=80) line,tag,part
integer(pInt) file,count,pos
integer(pInt), dimension(3) :: positions

part=''

do while(.true.)
   read(file,'(a80)',END=100) line
   positions=IO_stringPos(line,1)
   tag=IO_lc(IO_stringValue(line,positions,1))
   if (tag(1:1)=='<'.AND.tag(len_trim(tag):len_trim(tag))=='>') then 
      part=tag(2:len_trim(tag)-1)
	  exit
   elseif (tag(1:1)=='[') then  
	  count=count+1
   elseif (tag(2:len_trim(tag)-1)=='gauss') then
      texture_NGauss(count)=texture_NGauss(count)+1
   elseif (tag(2:len_trim(tag)-1)=='fiber') then
      texture_NFiber(count)=texture_NFiber(count)+1
   endif
enddo
100 return

end subroutine


character(len=80) function constitutive_Parse_UnknownPart(file)
!*********************************************************************
!* read an unknown "part" from the input file until                  *
!* the next part is reached                                          *
!* INPUT:                                                            *
!*  - file  : file ID                                                *
!*********************************************************************
use prec, only: pInt
use IO, only: IO_stringPos,IO_stringValue,IO_lc
implicit none

!* Definition of variables
character(len=80) line,tag
integer(pInt), parameter :: maxNchunks = 1
integer(pInt) file
integer(pInt), dimension(1+2*maxNchunks) :: positions 

constitutive_parse_unknownPart=''

do while(.true.)
   read(file,'(a80)',END=100) line
   positions=IO_stringPos(line,maxNchunks)
   tag=IO_lc(IO_stringValue(line,positions,1))
   if (tag(1:1)=='<'.AND.tag(len_trim(tag):len_trim(tag))=='>') then
      constitutive_Parse_UnknownPart=tag(2:len_trim(tag)-1)
	  exit
   endif
enddo
100 return

end function


character(len=80) function constitutive_Parse_MaterialPart(file)
!*********************************************************************
!* This function reads a material "part" from the input file until   *
!* the next part is reached                                          *
!* INPUT:                                                            *
!*  - file  : file ID                                                *
!*********************************************************************
use prec, only: pInt
use IO
implicit none

!* Definition of variables
character(len=80) line,tag
integer(pInt), parameter :: maxNchunks = 2! may be more than 2 chunks ..?
integer(pInt) file,section
integer(pInt), dimension(1+2*maxNchunks) :: positions 

section = 0
constitutive_parse_materialPart = ''

do while(.true.)
   read(file,'(a80)',END=100) line
   positions=IO_stringPos(line,maxNchunks) ! parse leading chunks
   tag=IO_lc(IO_stringValue(line,positions,1))
   if (tag(1:1)=='#') then  ! skip comment line
      cycle
   elseif (tag(1:1)=='<'.AND.tag(len_trim(tag):len_trim(tag))=='>') then
      constitutive_parse_materialPart=tag(2:len_trim(tag)-1)
	  exit
   elseif (tag(1:1)=='[') then
      section=section+1
   else
      if (section>0) then
         select case(tag)
	     case ('crystal_structure') 
              material_CrystalStructure(section)=IO_intValue(line,positions,2)
	     case ('nslip')
		      material_Nslip(section)=IO_intValue(line,positions,2)
		 case ('c11')
              material_C11(section)=IO_floatValue(line,positions,2)
		 case ('c12')
              material_C12(section)=IO_floatValue(line,positions,2)
		 case ('c13')
              material_C13(section)=IO_floatValue(line,positions,2)
		 case ('c33')
              material_C33(section)=IO_floatValue(line,positions,2)
		 case ('c44')
              material_C44(section)=IO_floatValue(line,positions,2)
         case ('s0_slip')
              material_s0_slip(section)=IO_floatValue(line,positions,2)
		 case ('gdot0_slip')
              material_gdot0_slip(section)=IO_floatValue(line,positions,2)
		 case ('n_slip')
              material_n_slip(section)=IO_floatValue(line,positions,2)
		 case ('h0')
              material_h0(section)=IO_floatValue(line,positions,2)
		 case ('s_sat')
              material_s_sat(section)=IO_floatValue(line,positions,2)
		 case ('w0')
              material_w0(section)=IO_floatValue(line,positions,2)
         end select
      endif
   endif
enddo
100 return

end function


character(len=80) function constitutive_Parse_TexturePart(file)
!*********************************************************************
!* This function reads a texture "part" from the input file until    *
!* the next part is reached                                          *
!* INPUT:                                                            *
!*  - file  : file ID                                                *
!*********************************************************************
use prec, only: pInt
use IO
implicit none

!* Definition of variables
character(len=80) line,tag
integer(pInt), parameter :: maxNchunks = 13 ! may be more than 10 chunks ..?
integer(pInt) file,pos,section,gaussCount,fiberCount,i
integer(pInt), dimension(1+2*maxNchunks) :: positions 

section = 0
gaussCount = 0
fiberCount = 0
constitutive_parse_texturePart = ''

do while(.true.)
   read(file,'(a80)',END=100) line
   positions=IO_stringPos(line,maxNchunks)  ! parse leading chunks
   tag=IO_lc(IO_stringValue(line,positions,1))
   if (tag(1:1)=='#') then  ! skip comment line
      cycle
   elseif (tag(1:1)=='<'.AND.tag(len_trim(tag):len_trim(tag))=='>') then
      constitutive_parse_texturePart=tag(2:len_trim(tag)-1)
	  exit
   elseif (tag(1:1)=='[') then
      section=section+1
	  gaussCount=0
	  fiberCount=0
   else
      if (section>0) then
         select case(tag)
  		 case ('hybridIA')
              texture_ODFfile(section)=IO_stringValue(line,positions,2)
		 case ('(gauss)')
		      gaussCount=gaussCount+1
		      do i=2,10,2
			     tag=IO_lc(IO_stringValue(line,positions,i))
				 select case (tag)
				 case('phi1')
				     texture_Gauss(1,gaussCount,section)=IO_floatValue(line,positions,i+1)
				 case('phi')
				     texture_Gauss(2,gaussCount,section)=IO_floatValue(line,positions,i+1)
				 case('phi2')
				     texture_Gauss(3,gaussCount,section)=IO_floatValue(line,positions,i+1)
				 case('scatter')
				     texture_Gauss(5,gaussCount,section)=IO_floatValue(line,positions,i+1)
				 case('fraction')
				     texture_Gauss(6,gaussCount,section)=IO_floatValue(line,positions,i+1)
				 end select
              enddo
		 case ('(fiber)')
              fiberCount=fiberCount+1
		      do i=2,12,2
			     tag=IO_lc(IO_stringValue(line,positions,i))
				 select case (tag)
				 case('alpha1')
				     texture_fiber(1,fiberCount,section)=IO_floatValue(line,positions,i+1)
				 case('alpha2')
				     texture_fiber(2,fiberCount,section)=IO_floatValue(line,positions,i+1)
				 case('beta1')
				     texture_fiber(3,fiberCount,section)=IO_floatValue(line,positions,i+1)
                 case('beta2')
				     texture_fiber(4,fiberCount,section)=IO_floatValue(line,positions,i+1)
				 case('scatter')
				     texture_fiber(5,fiberCount,section)=IO_floatValue(line,positions,i+1)
				 case('fraction')
				     texture_fiber(6,fiberCount,section)=IO_floatValue(line,positions,i+1)
				 end select
              enddo
		 case ('ngrains')
		      texture_Ngrains(section)=IO_intValue(line,positions,2)
         case ('symmetry')
		      texture_symmetry(section)=IO_stringValue(line,positions,2)
         end select
      endif
   endif
enddo
100 return

end function


subroutine constitutive_Parse_MatTexDat(filename)
!*********************************************************************
!* This function reads the material and texture input file           *
!* INPUT:                                                            *
!*  - filename : name of input file                                  *
!*********************************************************************
use prec, only: pReal,pInt
use IO
implicit none

!* Definition of variables
character(len=*) filename
character(len=80) part,formerPart
integer(pInt) sectionCount,dummy,i,j,m

!* First reading: number of materials and textures
!* determine material_maxN and texture_maxN
open(1,FILE=filename,ACTION='READ',STATUS='OLD',ERR=100)
part = '_dummy_'
do while (part/='')
   formerPart = part
   call constitutive_CountSections(1,sectionCount,part)
   select case (formerPart)
   case ('materials')
        material_maxN = sectionCount
   case ('textures')
        texture_maxN = sectionCount
   end select
enddo
close(1)
!* Arrays allocation
allocate(texture_NGauss(texture_maxN)) ; texture_NGauss=0_pInt
allocate(texture_NFiber(texture_maxN)) ; texture_NFiber=0_pInt

!* Second reading: number of Gauss and Fiber
!* determine material_maxN and texture_maxN
open(1,FILE=filename,ACTION='READ',STATUS='OLD',ERR=100)
part = '_dummy_'
sectionCount = 0
do while (part/='')
   select case (part)
   case ('textures') 
        call constitutive_CountGaussAndFiber(1,sectionCount,part)
   case default
        call constitutive_CountSections(1,dummy,part)
   end select
enddo
close(1)
!* Arrays allocation
texture_maxNGauss=maxval(texture_NGauss)
texture_maxNFiber=maxval(texture_NFiber)
allocate(material_CrystalStructure(material_maxN))		  ; material_CrystalStructure=0_pInt
allocate(material_Nslip(material_maxN))					  ; material_Nslip=0_pInt
allocate(material_C11(material_maxN))					  ; material_C11=0.0_pReal
allocate(material_C12(material_maxN))					  ; material_C12=0.0_pReal
allocate(material_C13(material_maxN))					  ; material_C13=0.0_pReal
allocate(material_C33(material_maxN))					  ; material_C33=0.0_pReal
allocate(material_C44(material_maxN))					  ; material_C44=0.0_pReal
allocate(material_s0_slip(material_maxN))				  ; material_s0_slip=0.0_pReal
allocate(material_gdot0_slip(material_maxN))			  ; material_gdot0_slip=0.0_pReal
allocate(material_n_slip(material_maxN))				  ; material_n_slip=0.0_pReal
allocate(material_h0(material_maxN))				      ; material_h0=0.0_pReal
allocate(material_s_sat(material_maxN))                   ; material_s_sat=0.0_pReal
allocate(material_w0(material_maxN))                      ; material_w0=0.0_pReal
allocate(texture_ODFfile(texture_maxN))                   ; texture_ODFfile=''
allocate(texture_Ngrains(texture_maxN))                   ; texture_Ngrains=0_pInt
allocate(texture_symmetry(texture_maxN))                  ; texture_symmetry=''
allocate(texture_Gauss(6,texture_maxNGauss,texture_maxN)) ; texture_Gauss=0.0_pReal
allocate(texture_Fiber(6,texture_maxNGauss,texture_maxN)) ; texture_Fiber=0.0_pReal

!* Third reading: materials and textures are stored
open(1,FILE=filename,ACTION='READ',STATUS='OLD',ERR=100)
part='_dummy_'
do while (part/='')
   select case (part)
   case ('materials')
	    part=constitutive_Parse_MaterialPart(1)
   case ('textures')
	    part=constitutive_Parse_TexturePart(1)
   case default
        part=constitutive_Parse_UnknownPart(1)
   end select
enddo
close(1)


!do m=1,material_maxN
!  material_Cslip_66(:,:,m) = 0.0_pReal
!  select case (material_crystal_structure)
!    case (1:2) ! cubic structure
!      do i=1,3
!        do j=1,3
!          material_Cslip_66(i,j,m)   = C12
!        enddo
!        material_Cslip_66(i,i,m)     = C11
!        material_Cslip_66(i+3,i+3,m) = C44
!      enddo
!    case (3)   ! hcp structure MISSING correct
!      do i=1,3
!        do j=1,3
!          material_Cslip_66(i,j,m)   = C12
!        enddo
!        material_Cslip_66(i,i,m)     = C11
!        material_Cslip_66(i+3,i+3,m) = C44
!      enddo
!  end select
!  material_Cslip_3333(:,:,:,:,m) = math_66to3333(Cslip_66(:,:,m))   
!end do


! MISSING some consistency checks may be..?

return
100 call IO_error(110) ! corrupt matarials_textures file
end subroutine


subroutine constitutive_Assignement()
!*********************************************************************
!* This subroutine assign material parameters according to ipc,ip,el *
!*********************************************************************
use prec, only: pReal,pInt
use mesh, only: mesh_NcpElems,FE_Nips,mesh_maxNips,mesh_element
!use CPFEM, only: CPFEM_Fp_old
implicit none

!* Definition of variables
integer(pInt) i,j,k,l

!* Allocate arrays
constitutive_maxNgrains=maxval(texture_Ngrains)
allocate(constitutive_Ngrains(mesh_maxNips,mesh_NcpElems))
constitutive_Ngrains=0_pInt
allocate(constitutive_matID(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
constitutive_matID=0_pInt
allocate(constitutive_texID(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
constitutive_texID=0_pInt
allocate(constitutive_MatVolFrac(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
constitutive_MatVolFrac=0.0_pReal
allocate(constitutive_TexVolFrac(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
constitutive_TexVolFrac=0.0_pReal
allocate(constitutive_state_old(material_maxNslip,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
constitutive_state_old=0.0_pReal
allocate(constitutive_state_new(material_maxNslip,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
constitutive_state_new=0.0_pReal

!* Assignement
do i=1,mesh_NcpElems
   do j=1,FE_Nips(mesh_element(2,i))
      !* Multiplicity of orientations per texture 
      constitutive_Ngrains(j,i)=texture_Ngrains(mesh_element(4,i))
      do k=1,constitutive_Ngrains(j,i)
	     !* MaterialID and TextureID 
         constitutive_matID(k,j,i)=mesh_element(3,i)
		 constitutive_texID(k,j,i)=mesh_element(4,i)
         constitutive_MatVolFrac(k,j,i)=1.0_pReal
!		 constitutive_TexVolFrac(k,j,i)=texture_VolFrac([gauss],mesh_element(4,i))
		 !* Initialization of state variables 
		 do l=1,material_Nslip(constitutive_matID(k,j,i))
		    constitutive_state_old(l,k,j,i)=material_s0_slip(constitutive_matID(k,j,i))
		    constitutive_state_new(l,k,j,i)=material_s0_slip(constitutive_matID(k,j,i))
		 enddo 
	  enddo
   enddo
enddo

end subroutine


!subroutine constitutive_InitFp(CPFEM_Fp_old)
!*********************************************************************
!* This function reads the material and texture input file           *
!* INPUT:                                                            *
!*  - CPFEM_Fp_old : old plastic deformation gradient                *
!*********************************************************************
!use prec, only: pReal,pInt
!use CPFEM, only: CPFEM_Fp_old
!implicit none

!* Definition of variables

!* Initialization of Fp_old with starting orientation


!end subroutine


subroutine constitutive_LpAndItsTangent(Tstar_v,ipc,ip,el,Lp,dLp_dTstar)
!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating the velocity gradient                                 *       
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor            *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:
!*  - Lp              : plastic velocity gradient                    *
!*  - dLp_dTstar      : derivative of Lp                             *
!*********************************************************************
use prec, only: pReal,pInt
implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
integer(pInt) matID,i,j,k
real(pReal) Tstar_v(6)
real(pReal) Lp(3,3)
real(pReal) dLp_dTstar(6,6)
real(pReal) dLpT_dTstar(6,6)
real(pReal), dimension(constitutive_matID(ipc,ip,el)) :: gdot_slip
real(pReal), dimension(constitutive_matID(ipc,ip,el)) :: dgdot_dtauslip
real(preal), dimension(constitutive_matID(ipc,ip,el)) :: tau_slip

!* Get the material-ID from the triplet(ipc,ip,el)
matID=constitutive_matID(ipc,ip,el)

!* Calculation of Lp
Lp=0.0_pReal
do i=1,material_Nslip(matID)
   tau_slip(i)=dot_product(Tstar_v,constitutive_Sslip_v(:,i,material_CrystalStructure(matID)))
   gdot_slip(i)=material_gdot0_slip(matID)*(abs(tau_slip(i))/constitutive_state_new(i,ipc,ip,el))**material_n_slip(matID)*sign(1.0_pReal,tau_slip(i))
   Lp=Lp+gdot_slip(i)*constitutive_Sslip(:,:,i,material_CrystalStructure(matID))
enddo

!* Calculation of the tangent of Lp
dLp_dTstar=0.0_pReal
do i=1,material_Nslip(matID)
   dgdot_dtauslip(i)=material_gdot0_slip(matID)*(abs(tau_slip(i))/constitutive_state_new(i,ipc,ip,el))**(material_n_slip(matID)-1.0_pReal)*material_n_slip(matID)/constitutive_state_new(i,ipc,ip,el)
   forall (j=1:6,k=1:6)
          dLp_dTstar(j,k)=dLp_dTstar(j,k)+constitutive_Sslip_v(j,i,material_CrystalStructure(matID))*constitutive_Sslip_v(k,i,material_CrystalStructure(matID))*dgdot_dtauslip(i) 
   endforall
enddo

return
end subroutine


function constitutive_DotState(Tstar_v,ipc,ip,el)
!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating the velocity gradient                                 *       
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor            *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - constitutive_DotState : evolution of state variable            *
!*********************************************************************
use prec, only: pReal,pInt
implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
integer(pInt) matID,i
real(pReal) Tstar_v(6)
real(pReal), dimension(constitutive_matID(ipc,ip,el)) :: constitutive_DotState
real(pReal), dimension(constitutive_matID(ipc,ip,el)) :: gdot_slip
real(pReal), dimension(constitutive_matID(ipc,ip,el)) :: tau_slip
real(pReal), dimension(constitutive_matID(ipc,ip,el)) :: self_hardening
 
!* Get the material-ID from the triplet(ipc,ip,el)
matID=constitutive_matID(ipc,ip,el)

!* Self-Hardening of each system
do i=1,material_Nslip(matID)
   tau_slip(i)=dot_product(Tstar_v,constitutive_Sslip_v(:,i,material_CrystalStructure(matID)))
   gdot_slip(i)=material_gdot0_slip(matID)*(abs(tau_slip(i))/constitutive_state_new(i,ipc,ip,el))**material_n_slip(matID)*sign(1.0_pReal,tau_slip(i))
   self_hardening(i)=material_h0(matID)*(1.0_pReal-constitutive_state_new(i,ipc,ip,el)/material_s_sat(matID))**material_w0(matID)*abs(gdot_slip(i))
enddo

!* Hardening for all systems
constitutive_DotState=matmul(constitutive_hardening_matrix(1:material_Nslip(matID),1:material_Nslip(matID),material_CrystalStructure(matID)),self_hardening)

return
end function


END MODULE
