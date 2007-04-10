! QUESTION fileunit 1  may run into trouble?


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

! MISSING consistency check after reading 'mattex.mpie'
character(len=300), parameter :: mattexFile = 'mattex.mpie'

!*************************************
!* Definition of material properties *
!*************************************
!* Number of materials
integer(pInt) material_maxN
!* Crystal structure and number of selected slip systems per material
integer(pInt), dimension(:)  , allocatable :: material_CrystalStructure
integer(pInt), dimension(:)  , allocatable :: material_Nslip
!* Maximum number of selected slip systems over materials
integer(pInt) material_maxNslip
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
integer(pInt),dimension(:)     , allocatable :: texture_NRandom
integer(pInt),dimension(:)     , allocatable :: texture_totalNgrains
real(pReal), dimension(:,:,:)  , allocatable :: texture_Gauss
real(pReal), dimension(:,:,:)  , allocatable :: texture_Fiber
real(pReal), dimension(:,:,:,:), allocatable :: constitutive_EulerAngles

!************************************
!*             Grains               *
!************************************
integer(pInt) constitutive_maxNgrains
integer(pInt), dimension(:,:)    , allocatable :: constitutive_Ngrains
integer(pInt), dimension(:,:,:)  , allocatable :: constitutive_matID
real(pReal), dimension(:,:,:)    , allocatable :: constitutive_matVolFrac
integer(pInt), dimension(:,:,:)  , allocatable :: constitutive_texID
real(pReal), dimension(:,:,:)    , allocatable :: constitutive_texVolFrac

!************************************
!*         State variables          *
!************************************
integer(pInt) constitutive_maxNstatevars
integer(pInt), dimension(:,:,:), allocatable :: constitutive_Nstatevars
real(pReal), dimension(:,:,:,:), allocatable :: constitutive_state_old
real(pReal), dimension(:,:,:,:), allocatable :: constitutive_state_new

!************************************
!*             Results              *
!************************************
integer(pInt) constitutive_maxNresults
integer(pInt), dimension(:,:,:), allocatable :: constitutive_Nresults
real(pReal), dimension(:,:,:,:), allocatable :: constitutive_results
!* MISSING : allocation


!***********************************************
!*      Crystal structures                     *
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

!***********************************************
!*      slip-slip interaction                  *
!***********************************************
!* (defined for the moment as crystal structure property and not as material property)
!* (may be changed in the future)
real(pReal), dimension(constitutive_MaxMaxNslipOfStructure,constitutive_MaxMaxNslipOfStructure,&
                       constitutive_MaxCrystalStructure) :: constitutive_HardeningMatrix
real(pReal), parameter :: constitutive_LatentHardening=1.4_pReal


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
!* - constitutive_Assignment
!* - constitutive_HomogenizedC                          
!* - constitutive_LpAndItsTangent     
!* - consistutive_DotState          
!****************************************


subroutine constitutive_Init()
!**************************************
!*      Module initialization         *
!**************************************
call constitutive_SchmidMatrices()
call constitutive_HardeningMatrices()
call constitutive_Parse_MatTexDat(mattexFile)
call constitutive_Assignment()
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
      invNorm=dsqrt(1.0_pReal/((constitutive_sn(1,k,l)**2+constitutive_sn(2,k,l)**2+constitutive_sn(3,k,l)**2)*&
	          (constitutive_sd(1,k,l)**2+constitutive_sd(2,k,l)**2+constitutive_sd(3,k,l)**2)))
      constitutive_Sslip(:,:,k,l)=constitutive_Sslip(:,:,k,l)*invNorm
!* Vectorization of normalized Schmid matrix
      constitutive_Sslip_v(1,k,l)=constitutive_Sslip(1,1,k,l)
      constitutive_Sslip_v(2,k,l)=constitutive_Sslip(2,2,k,l)
      constitutive_Sslip_v(3,k,l)=constitutive_Sslip(3,3,k,l)
	  !* be compatible with Mandel notation of Tstar
      constitutive_Sslip_v(4,k,l)=(constitutive_Sslip(1,2,k,l)+constitutive_Sslip(2,1,k,l))/dsqrt(2.0_pReal)
      constitutive_Sslip_v(5,k,l)=(constitutive_Sslip(2,3,k,l)+constitutive_Sslip(3,3,k,l))/dsqrt(2.0_pReal)
      constitutive_Sslip_v(6,k,l)=(constitutive_Sslip(1,3,k,l)+constitutive_Sslip(3,1,k,l))/dsqrt(2.0_pReal)
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
constitutive_HardeningMatrix=constitutive_LatentHardening
!* Iteration over the crystal structures 
do l=1,3
   select case(l) 
!* Hardening matrix for FCC structures    
   case (1)
      forall (k=1:10:3,i=0:2,j=0:2)
         constitutive_HardeningMatrix(k+i,k+j,l)=1.0_pReal
      endforall
!* Hardening matrix for BCC structures
   case (2)
      forall (k=1:11:2,i=0:1,j=0:1)
             constitutive_HardeningMatrix(k+i,k+j,l)=1.0_pReal
      endforall
      forall (k=13:48)
         constitutive_HardeningMatrix(k,k,l)=1.0_pReal
      endforall
!* Hardening matrix for HCP structures
   case (3)
      forall (i=1:3,j=1:3)
          constitutive_HardeningMatrix(i,j,l)=1.0_pReal
      endforall
      forall (k=4:12)
         constitutive_HardeningMatrix(k,k,l)=1.0_pReal
      endforall
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

do
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


character(len=80) function constitutive_assignNGaussAndFiber(file)
!*********************************************************************
!*********************************************************************
use prec, only: pInt
use IO,   only: IO_stringPos,IO_stringValue,IO_lc
implicit none

!* Definition of variables
character(len=80) line,tag
integer(pInt) file,section,pos
integer(pInt), dimension(3) :: positions

constitutive_assignNGaussAndFiber=''
section = 0_pInt

do
   read(file,'(a80)',END=100) line
   positions=IO_stringPos(line,1)
   tag=IO_lc(IO_stringValue(line,positions,1))
   if (tag(1:1)=='<'.AND.tag(len_trim(tag):len_trim(tag))=='>') then 
      constitutive_assignNGaussAndFiber=tag(2:len_trim(tag)-1)
	  exit
   elseif (tag(1:1)=='[') then  
	  section=section+1
      texture_NGauss(section) = 0_pInt
      texture_NFiber(section) = 0_pInt
   elseif (tag=='(gauss)') then
      texture_NGauss(section)=texture_NGauss(section)+1
   elseif (tag=='(fiber)') then
      texture_NFiber(section)=texture_NFiber(section)+1
   endif
enddo
100 return

end function


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

do
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
integer(pInt), parameter :: maxNchunks = 2 
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
use math, only: inRad
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
				     texture_Gauss(1,gaussCount,section)=IO_floatValue(line,positions,i+1)*inRad
				 case('phi')
				     texture_Gauss(2,gaussCount,section)=IO_floatValue(line,positions,i+1)*inRad
				 case('phi2')
				     texture_Gauss(3,gaussCount,section)=IO_floatValue(line,positions,i+1)*inRad
				 case('scatter')
				     texture_Gauss(4,gaussCount,section)=IO_floatValue(line,positions,i+1)
				 case('fraction')
				     texture_Gauss(5,gaussCount,section)=IO_floatValue(line,positions,i+1)
				 end select
              enddo
		 case ('(fiber)')
              fiberCount=fiberCount+1
		      do i=2,12,2
			     tag=IO_lc(IO_stringValue(line,positions,i))
				 select case (tag)
				 case('alpha1')
				     texture_fiber(1,fiberCount,section)=IO_floatValue(line,positions,i+1)*inRad
				 case('alpha2')
				     texture_fiber(2,fiberCount,section)=IO_floatValue(line,positions,i+1)*inRad
				 case('beta1')
				     texture_fiber(3,fiberCount,section)=IO_floatValue(line,positions,i+1)*inRad
                 case('beta2')
				     texture_fiber(4,fiberCount,section)=IO_floatValue(line,positions,i+1)*inRad
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
use IO, only: IO_error
use math, only: math_Mandel3333to66, math_Voigt66to3333
implicit none

!* Definition of variables
character(len=*) filename
character(len=80) part,formerPart
integer(pInt) sectionCount,i,j,k

!-----------------------------
!* First reading: number of materials and textures
!-----------------------------
!* determine material_maxN and texture_maxN from last respective parts
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
!* Array allocation
allocate(material_CrystalStructure(material_maxN))		  ; material_CrystalStructure=0_pInt
allocate(material_Nslip(material_maxN))					  ; material_Nslip=0_pInt
allocate(material_C11(material_maxN))					  ; material_C11=0.0_pReal
allocate(material_C12(material_maxN))					  ; material_C12=0.0_pReal
allocate(material_C13(material_maxN))					  ; material_C13=0.0_pReal
allocate(material_C33(material_maxN))					  ; material_C33=0.0_pReal
allocate(material_C44(material_maxN))					  ; material_C44=0.0_pReal
allocate(material_Cslip_66(6,6,material_maxN))            ; material_Cslip_66=0.0_pReal
allocate(material_s0_slip(material_maxN))				  ; material_s0_slip=0.0_pReal
allocate(material_gdot0_slip(material_maxN))			  ; material_gdot0_slip=0.0_pReal
allocate(material_n_slip(material_maxN))				  ; material_n_slip=0.0_pReal
allocate(material_h0(material_maxN))				      ; material_h0=0.0_pReal
allocate(material_s_sat(material_maxN))                   ; material_s_sat=0.0_pReal
allocate(material_w0(material_maxN))                      ; material_w0=0.0_pReal
allocate(texture_ODFfile(texture_maxN))                   ; texture_ODFfile=''
allocate(texture_Ngrains(texture_maxN))                   ; texture_Ngrains=0_pInt
allocate(texture_symmetry(texture_maxN))                  ; texture_symmetry=''
allocate(texture_NGauss(texture_maxN)) ; texture_NGauss=0_pInt
allocate(texture_NFiber(texture_maxN)) ; texture_NFiber=0_pInt
allocate(texture_NRandom(texture_maxN)) ; texture_NRandom=0_pInt

!-----------------------------
!* Second reading: number of Gauss and Fiber
!-----------------------------
rewind(1)
part = '_dummy_'
do while (part/='')
   select case (part)
   case ('textures') 
        part = constitutive_assignNGaussAndFiber(1)
   case default
        part = constitutive_Parse_UnknownPart(1)
   end select
enddo
!* Array allocation
texture_maxNGauss=maxval(texture_NGauss)
texture_maxNFiber=maxval(texture_NFiber)
allocate(texture_Gauss(5,texture_maxNGauss,texture_maxN)) ; texture_Gauss=0.0_pReal
allocate(texture_Fiber(6,texture_maxNFiber,texture_maxN)) ; texture_Fiber=0.0_pReal

!-----------------------------
!* Third reading: materials and textures are stored
!-----------------------------
rewind(1)
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


!* Construction of the elasticity matrices
do i=1,material_maxN
   select case (material_CrystalStructure(i))
   case(1:2) ! cubic(s)
       do k=1,3
          do j=1,3
             material_Cslip_66(k,j,i)=material_C12(i)
          enddo
          material_Cslip_66(k,k,i)=material_C11(i)
          material_Cslip_66(k+3,k+3,i)=material_C44(i)
       enddo
   case(3)   ! hcp
        material_Cslip_66(1,1,i)=material_C11(i)
        material_Cslip_66(2,2,i)=material_C11(i)
		material_Cslip_66(3,3,i)=material_C33(i)
		material_Cslip_66(1,2,i)=material_C12(i)
		material_Cslip_66(2,1,i)=material_C12(i)
		material_Cslip_66(1,3,i)=material_C13(i)
		material_Cslip_66(3,1,i)=material_C13(i)
		material_Cslip_66(2,3,i)=material_C13(i)
		material_Cslip_66(3,2,i)=material_C13(i)
        material_Cslip_66(4,4,i)=material_C44(i)
        material_Cslip_66(5,5,i)=material_C44(i)
        material_Cslip_66(6,6,i)=0.5_pReal*(material_C11(i)-material_C12(i))  
   end select
   material_Cslip_66(:,:,i) = math_Mandel3333to66(math_Voigt66to3333(material_Cslip_66(:,:,i)))
enddo


! MISSING some consistency checks may be..?
! if ODFfile present then set NGauss NFiber =0
return
100 call IO_error(110) ! corrupt materials_textures file
end subroutine


subroutine constitutive_Assignment()
!*********************************************************************
!* This subroutine assign material parameters according to ipc,ip,el *
!*********************************************************************
use prec, only: pReal,pInt
use math, only: math_sampleGaussOri,math_sampleFiberOri,math_sampleRandomOri,math_symmetricEulers,math_EulerToR
use mesh, only: mesh_NcpElems,FE_Nips,FE_mapElemtype,mesh_maxNips,mesh_element
use IO,   only: IO_hybridIA
use CPFEM,only: CPFEM_Fp_old
implicit none

!* Definition of variables
integer(pInt) e,i,j,k,l,m,o,g,s,Ngrains
integer(pInt) matID,texID
integer(pInt), dimension(:,:,:), allocatable :: hybridIA_population
integer(pInt), dimension(texture_maxN) :: Ncomponents,Nsym,multiplicity,sumVolfrac,ODFmap,sampleCount
real(pReal), dimension(3,4*(1+texture_maxNGauss+texture_maxNfiber)) :: Euler
real(pReal), dimension(4*(1+texture_maxNGauss+texture_maxNfiber)) :: texVolfrac

! process textures
o = 0_pInt       ! ODF counter
ODFmap = 0_pInt  ! blank mapping
sampleCount = 0_pInt  ! count orientations assigned per texture

do texID=1,texture_maxN
   if (texture_ODFfile(texID)=='') then
      sumVolfrac(texID) = sum(texture_gauss(5,:,texID))+sum(texture_fiber(6,:,texID))
      if (sumVolfrac(texID)<1.0_pReal) texture_NRandom(texID) = 1_pInt  ! check whether random component missing
      select case (texture_symmetry(texID))                             ! set symmetry factor
	  case ('orthotropic')
		   Nsym(texID) = 4_pInt
	  case ('monoclinic')
		   Nsym(texID) = 2_pInt
	  case default
	 	   Nsym(texID) = 1_pInt
	  end select
      Ncomponents(texID) = texture_NGauss(texID)+texture_NFiber(texID)+texture_NRandom(texID)
   else      ! hybrid IA
      o = o+1
	  ODFmap(texID) = o             ! remember mapping
      Ncomponents(texID) = 1_pInt   ! single "component"
      Nsym(texID) = 1_pInt          ! no symmetry (use full ODF instead)
   endif
! adjust multiplicity and number of grains per IP of components
   multiplicity(texID) = max(1_pInt,texture_Ngrains(texID)/Ncomponents(texID)/Nsym(texID))
   if (mod(texture_Ngrains(texID),Ncomponents(texID)*Nsym(texID)) /= 0_pInt) then
      texture_Ngrains(texID) = multiplicity(texID)*Ncomponents(texID)*Nsym(texID)
      write (6,*) 'changed Ngrains to',texture_Ngrains(texID),' for texture',texID
   endif
enddo

!* publish globals
constitutive_maxNgrains = maxval(texture_Ngrains)
constitutive_maxNstatevars = material_maxNslip

!* calc texture_totalNgrains
allocate(texture_totalNgrains(texture_maxN)) ; texture_totalNgrains=0_pInt
do i=1,mesh_NcpElems
   texID = mesh_element(4,i)
   texture_totalNgrains(texID) = texture_totalNgrains(texID) + FE_Nips(FE_mapElemtype(mesh_element(2,i)))*texture_Ngrains(texID)
enddo

! generate hybridIA samplings for ODFfile textures to later draw from these populations
allocate(hybridIA_population(3,maxval(texture_totalNgrains,ODFmap /= 0),o))
do texID = 1,texture_maxN
   if (ODFmap(texID) > 0) &
      hybridIA_population(:,:,ODFmap(texID)) = IO_hybridIA(texture_totalNgrains(texID),texture_ODFfile(texID))
enddo

!* Array allocation
allocate(constitutive_Ngrains(mesh_maxNips,mesh_NcpElems)) ; constitutive_Ngrains=0_pInt
allocate(constitutive_matID(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_matID=0_pInt
allocate(constitutive_texID(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_texID=0_pInt
allocate(constitutive_MatVolFrac(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_MatVolFrac=0.0_pReal
allocate(constitutive_TexVolFrac(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_TexVolFrac=0.0_pReal
allocate(constitutive_Nstatevars(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_Nstatevars=0_pInt
allocate(constitutive_state_old(constitutive_maxNstatevars,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
constitutive_state_old=0.0_pReal
allocate(constitutive_state_new(constitutive_maxNstatevars,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
constitutive_state_new=0.0_pReal
allocate(constitutive_Nresults(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_Nresults=0_pInt

!* Assignment of all grains in all IPs of all cp-elements
do e=1,mesh_NcpElems
   matID=mesh_element(3,e)
   texID=mesh_element(4,e)
   do i=1,FE_Nips(FE_mapElemtype(mesh_element(2,e)))
      g = 0_pInt     ! grain counter
	  do m = 1,multiplicity(texID)
		 o = 0_pInt  ! component counter
         if (texture_ODFfile(texID)=='') then
		    do k = 1,texture_nGauss(texID)            ! *** gauss ***
			   o = o+1
               Euler(:,o)    = math_sampleGaussOri(texture_Gauss(1:3,k,texID),texture_Gauss(4,k,texID))
		       texVolFrac(o) = texture_Gauss(5,k,texID)
			enddo
		    do k = 1,texture_nFiber(texID)            ! *** fiber ***
			   o = o+1
               Euler(:,o)    = math_sampleFiberOri(texture_Fiber(1:2,k,texID),texture_Fiber(3:4,k,texID),texture_Fiber(5,k,texID))
		       texVolFrac(o) = texture_Fiber(6,k,texID)
			enddo
		    do k = 1,texture_nRandom(texID)           ! *** random ***
			   o = o+1
               Euler(:,o)    = math_sampleRandomOri()
		       texVolfrac(o) = 1.0_pReal-sumVolfrac(texID)
			enddo
		 else                                         ! *** hybrid IA ***
		    o = 1 ! only singular orientation, i.e. single "component"
			Euler(:,o) = hybridIA_population(:,1+sampleCount(texID),ODFmap(texID))
            texVolfrac(o) = 1.0_pReal
		 endif
		 if (Nsym(texID) > 1) then   ! symmetry generates additional orientations
		    forall (k=1:o)
		       Euler(:,1+o+(Nsym(texID)-1)*(k-1):3+o+(Nsym(texID)-1)*(k-1)) = &
		       math_symmetricEulers(texture_symmetry(texID),Euler(:,k))
		       texVolfrac(1+o+(Nsym(texID)-1)*(k-1):3+o+(Nsym(texID)-1)*(k-1)) = texVolfrac(k)
		    end forall
         endif
         do s = 1,Nsym(texID)*o   ! loop over orientations to be assigned to ip (ex multiplicity)
		    g = g+1               ! next "grain"
			sampleCount(texID) = sampleCount(texID)+1  ! next member of population
            constitutive_matID(g,i,e) = matID          ! copy matID of element
		    constitutive_texID(g,i,e) = texID          ! copy texID of element
		    constitutive_MatVolFrac(g,i,e) = 1.0_pReal ! singular material (so far)
		    constitutive_TexVolFrac(g,i,e) = texVolfrac(s)/multiplicity(texID)/Nsym(texID)
		    constitutive_Nstatevars(g,i,e) = material_Nslip(matID) ! number of state variables (i.e. tau_c of each slip system)
		    CPFEM_Fp_old(:,:,g,i,e) = math_EulerToR(Euler(:,s))    ! set plastic deformation gradient at t_0
            forall (l=1:constitutive_Nstatevars(g,i,e))  ! initialize state variables
               constitutive_state_old(l,g,i,e) = material_s0_slip(matID)
               constitutive_state_new(l,g,i,e) = material_s0_slip(matID)
            end forall
		 enddo  ! components
	  enddo  ! multiplicity
   enddo ! ip
enddo ! cp_element


end subroutine


function constitutive_HomogenizedC(ipc,ip,el)
!*********************************************************************
!* This function returns the homogenized elacticity matrix           *
!* INPUT:                                                            *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
use prec, only: pReal,pInt
implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
real(pReal), dimension(6,6) :: constitutive_homogenizedC

!* Homogenization scheme
constitutive_homogenizedC=constitutive_MatVolFrac(ipc,ip,el)*material_Cslip_66(:,:,constitutive_matID(ipc,ip,el)) 

return
end function


subroutine constitutive_LpAndItsTangent(Lp,dLp_dTstar, Tstar_v,state,ipc,ip,el)
!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating the velocity gradient                                 *       
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - state           : current microstructure                       *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - Lp              : plastic velocity gradient                    *
!*  - dLp_dTstar      : derivative of Lp (4th-order tensor)          *
!*********************************************************************
use prec, only: pReal,pInt
implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
integer(pInt) matID,i,k,l,m,n
real(pReal) Tstar_v(6)
real(pReal) Lp(3,3)
real(pReal) dLp_dTstar(3,3,3,3)
real(pReal), dimension(constitutive_Nstatevars(ipc,ip,el)) :: state,gdot_slip,dgdot_dtauslip,tau_slip

!* Get the material-ID from the triplet(ipc,ip,el)
matID = constitutive_matID(ipc,ip,el)

!* Calculation of Lp
Lp = 0.0_pReal
do i=1,material_Nslip(matID)
   tau_slip(i)=dot_product(Tstar_v,constitutive_Sslip_v(:,i,material_CrystalStructure(matID)))
   gdot_slip(i)=material_gdot0_slip(matID)*(abs(tau_slip(i))/state(i))**&
                material_n_slip(matID)*sign(1.0_pReal,tau_slip(i))
   Lp=Lp+gdot_slip(i)*constitutive_Sslip(:,:,i,material_CrystalStructure(matID))
enddo

!* Calculation of the tangent of Lp
dLp_dTstar=0.0_pReal
do i=1,material_Nslip(matID)
   dgdot_dtauslip(i)=material_gdot0_slip(matID)*(abs(tau_slip(i))/state(i))**&
                     (material_n_slip(matID)-1.0_pReal)*material_n_slip(matID)/constitutive_state_new(i,ipc,ip,el)
   forall (k=1:3,l=1:3,m=1:3,n=1:3)
          dLp_dTstar(k,l,m,n)=dLp_dTstar(k,l,m,n)+constitutive_Sslip(k,l,i,material_CrystalStructure(matID))*&
		                      constitutive_Sslip(m,n,i,material_CrystalStructure(matID))*dgdot_dtauslip(i) 
   endforall
enddo

return
end subroutine


function constitutive_DotState(Tstar_v,state,ipc,ip,el)
!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating the velocity gradient                                 *       
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - state           : current microstructure                       *
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
real(pReal), dimension(6) :: Tstar_v
real(pReal), dimension(constitutive_Nstatevars(ipc,ip,el)) :: constitutive_DotState,&
                                                              state,gdot_slip,tau_slip,self_hardening
 
!* Get the material-ID from the triplet(ipc,ip,el)
matID = constitutive_matID(ipc,ip,el)

!* Self-Hardening of each system
do i=1,constitutive_Nstatevars(ipc,ip,el)
   tau_slip(i)=dot_product(Tstar_v,constitutive_Sslip_v(:,i,material_CrystalStructure(matID)))
   gdot_slip(i)=material_gdot0_slip(matID)*(abs(tau_slip(i))/state(i))**&
                material_n_slip(matID)*sign(1.0_pReal,tau_slip(i))
   self_hardening(i)=material_h0(matID)*(1.0_pReal-state(i)/&
                material_s_sat(matID))**material_w0(matID)*abs(gdot_slip(i))
enddo

!* Hardening for all systems
constitutive_DotState=matmul(constitutive_HardeningMatrix(1:material_Nslip(matID),1:material_Nslip(matID),&
                      material_CrystalStructure(matID)),self_hardening)

return
end function


END MODULE
