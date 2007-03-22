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

!*****************************
!*     Module parameters     *
!*****************************
!* Character *
character(len=80), allocatable :: constitutive_ODFfile(:)
! NB: ODFfile(number of texture)
character(len=80), allocatable :: constitutive_symmetry(:)
! NB: symmetry(number of texture)

!* Integer *
integer(pInt) constitutive_Nmats
! NB: Number of materials (read in material file)
integer(pInt) constitutive_Ntexts
! NB: Number of textures (read in material file)
integer(pInt), allocatable :: constitutive_crystal_structure(:)
! NB: crystal_structure(number of material)=1-3
integer(pInt), allocatable :: constitutive_Nslip(:)
! NB: Number of systems for each material
integer(pInt) constitutive_Nslip_max(3)
! NB: Number of defines slip systems
integer(pInt), allocatable :: constitutive_Ngrains(:)
! NB: Ngrains(number of texture)
 
!* Real *
real(pReal), allocatable :: constitutive_C11(:)
real(pReal), allocatable :: constitutive_C12(:)
real(pReal), allocatable :: constitutive_C13(:)
real(pReal), allocatable :: constitutive_C33(:)
real(pReal), allocatable :: constitutive_C44(:)
real(pReal), allocatable :: constitutive_Cslip_66(:,:,:)
! NB: Cslip_66(1:6,1:6,number of materials)
real(pReal), allocatable :: constitutive_s0_slip(:)
real(pReal), allocatable :: constitutive_gdot0_slip(:)
real(pReal), allocatable :: constitutive_n_slip(:)
real(pReal), allocatable :: constitutive_h0(:)
real(pReal), allocatable :: constitutive_s_sat(:)
real(pReal), allocatable :: constitutive_w0(:)
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
data constitutive_Nslip_max(1)/12/
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
data constitutive_Nslip_max(2)/48/
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
data constitutive_Nslip_max(3)/12/
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


CONTAINS
!****************************************
!* - constitutive_init                  *
!* - constitutive_calc_SchmidM          *
!* - constitutive_calc_HardeningM       *
!* - constitutive_parse_MatTexDat       *            
!* - constitutive_calc_SlipRates        *
!* - constitutive_calc_Hardening        *
!* - consistutive_calc_PlasVeloGradient *
!* - CPFEM_CauchyStress???????          *
!****************************************


subroutine constitutive_init()
!**************************************
!*      Module initialization         *
!**************************************
call constitutive_calc_SchmidM()
call constitutive_calc_hardeningM()
call constitutive_parse_MatTexDat('materials_textures.mpie')
end subroutine
 

subroutine constitutive_calc_SchmidM()
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
   do k=1,constitutive_Nslip_max(l)
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


!***********************************************************
!***********************************************************
SUBROUTINE constitutive_countSections(file,count,part)

use prec, only: pInt
use IO, only: IO_stringPos,IO_stringValue,IO_lc
implicit none

integer(pInt) file,count,pos
integer(pInt), dimension(3) :: positions
character(len=80) part,line,tag

count = 0
part = ''

do
  read(unit=file,fmt='(a80)',end=100) line
  positions = IO_stringPos(line,1)
  tag = IO_lc(IO_stringValue(line,positions,1))
  if (tag(1:1)=='<' .and. tag(len(tag):len(tag)=='>') then
    part = tag(2:len(tag)-1)
	exit
  elseif (tag(1:1)=='[' .and. tag(len(tag):len(tag)==']') then
    count = count+1
  endif
end do
100 return

END SUBROUTINE


!***********************************************************
!***********************************************************
character(len=80) function constitutive_parse_unknownPart(file)

use prec, only: pInt
use IO, only: IO_stringPos,IO_stringValue,IO_lc
implicit none

integer(pInt), parameter :: maxNchunks = 1
integer(pInt) file
integer(pInt), dimension(1+2*maxNchunks) :: positions 
character(len=80) line,tag

constitutive_parse_unknownPart = ''

do
  read(unit=file,fmt='(a80)',end=100) line
  positions = IO_stringPos(line,maxNchunks)
  tag = IO_lc(IO_stringValue(line,positions,1))
  if (tag(1:1)=='<' .and. tag(len(tag):len(tag)=='>') then
    constitutive_parse_unknownPart = tag(2:len(tag)-1)
	exit
  endif
end do

100 return

END FUNCTION


!***********************************************************
!***********************************************************
character(len=80) function constitutive_parse_materialPart(file)

use prec, only: pInt
use IO, only: IO_stringPos,IO_stringValue,IO_lc
implicit none

integer(pInt), parameter :: maxNchunks = 2! may be more than 2 chunks ..?
integer(pInt) file,section
integer(pInt), dimension(1+2*maxNchunks) :: positions 
character(len=80) line,tag

section = 0
constitutive_parse_materialPart = ''

do
  read(unit=file,fmt='(a80)',end=100) line
  positions = IO_stringPos(line,maxNchunks) ! parse leading chunks
  tag = IO_lc(IO_stringValue(line,positions,1))
  if (tag(1:1)=='#') then  ! skip comment line
    cycle
  elseif (tag(1:1)=='<' .and. tag(len(tag):len(tag)=='>') then
    constitutive_parse_materialPart = tag(2:len(tag)-1)
	exit
  elseif (tag(1:1)=='[' .and. tag(len(tag):len(tag)==']') then
    section = section+1
  else
    if (section>0) then
      select case(tag)
	    case ('crystal_structure') ! crystal structure
          constitutive_crystal_structure(section)=IO_intValue(line,positions,2)
	    case ('nslip')
		  constitutive_Nslip(section)=IO_intValue(line,positions,2)
		case ('C11')
          constitutive_C11(section)=IO_floatValue(line,positions,2)
		case ('C12')
          constitutive_C12(section)=IO_floatValue(line,positions,2)
		case ('C13')
          constitutive_C13(section)=IO_floatValue(line,positions,2)
		case ('C33')
          constitutive_C33(section)=IO_floatValue(line,positions,2)
		case ('C44')
          constitutive_C44(section)=IO_floatValue(line,positions,2)
        case ('s0_slip')
          constitutive_s0_slip(section)=IO_floatValue(line,positions,2)
		case ('gdot0_slip')
          constitutive_gdot0_slip(section)=IO_floatValue(line,positions,2)
		case ('n_slip')
          constitutive_n_slip(section)=IO_floatValue(line,positions,2)
		case ('h0')
          constitutive_h0(section)=IO_floatValue(line,positions,2)
		case ('s_sat')
          constitutive_s_sat(section)=IO_floatValue(line,positions,2)
		case ('w0')
          constitutive_w0(section)=IO_floatValue(line,positions,2)
        case default
          write(6,*) 'Unknown material parameter ',line
      end select
    endif
  endif
end do

100 return

END FUNCTION


!***********************************************************
!***********************************************************
character(len=80) function constitutive_parse_texturePart(file)

use prec, only: pInt
use IO, only: IO_stringPos,IO_stringValue,IO_lc
implicit none

integer(pInt), parameter :: maxNchunks = 10 ! may be more than 10 chunks ..?
integer(pInt) file,pos,section
integer(pInt), dimension(1+2*maxNchunks) :: positions 
character(len=80) line,tag

section = 0
constitutive_parse_texturePart = ''

do
  read(unit=file,fmt='(a80)',end=100) line
  positions = IO_stringPos(line,maxNchunks)  ! parse leading chunks
  tag = IO_lc(IO_stringValue(line,positions,1))
  if (tag(1:1)=='#') then  ! skip comment line
    cycle
  elseif (tag(1:1)=='<' .and. tag(len(tag):len(tag)=='>') then
    constitutive_parse_texturePart = tag(2:len(tag)-1)
	exit
  elseif (tag(1:1)=='[' .and. tag(len(tag):len(tag)==']') then
    section = section+1
  else
    if (section>0) then
      select case(tag)
  		case ('hybridIA')
          constitutive_ODFfile(section)=IO_stringValue(line,positions,2)
		case ('gauss')
		case ('fiber')
		case ('ngrains')
		  constitutive_Ngrains(section)=IO_intValue(line,positions,2)
        case ('symmetry')
		  constitutive_symmetry(section)=IO_stringValue(line,positions,2)
        case default
          write(6,*) 'Unknown texture parameter ',line
      end select
    endif
  endif
end do

100 return

END FUNCTION



subroutine constitutive_parse_MatTexDat(filename)
!***********************************************************
!* Reading material parameters and texture components file *
!***********************************************************
use prec, only: pReal,pInt
use IO
implicit none

character(len=*) filename
character(len=80) part,formerPart
integer(pInt) sectionCount,i,j,m


open(200,FILE=filename,ACTION='READ',STATUS='OLD',ERR=100)

part = '_dummy_'
do while (part/='')
  formerPart = part
  call constitutive_countSections(200,sectionCount,part)
  select case (formerPart)
    case ('materials')
      materials_maxN = sectionCount
    case ('textures')
      textures_maxN = sectionCount
  end select
end do

allocate(constitutive_ODFfile(constitutive_Ntexts))          ; constitutive_ODFfile=''
allocate(constitutive_Ngrains(constitutive_Ntexts))          ; constitutive_Ngrains=0_pInt
allocate(constitutive_symmetry(constitutive_Ntexts))         ; constitutive_symmetry=''
allocate(constitutive_crystal_structure(constitutive_Nmats)) ; constitutive_crystal_structure=0_pInt
allocate(constitutive_Nslip(constitutive_Nmats))             ; constitutive_Nslip=0_pInt
allocate(constitutive_C11(constitutive_Nmats))               ; constitutive_C11=0.0_pReal
allocate(constitutive_C12(constitutive_Nmats))               ; constitutive_C12=0.0_pReal
allocate(constitutive_C13(constitutive_Nmats))               ; constitutive_C13=0.0_pReal
allocate(constitutive_C33(constitutive_Nmats))               ; constitutive_C33=0.0_pReal
allocate(constitutive_C44(constitutive_Nmats))               ; constitutive_C44=0.0_pReal
allocate(constitutive_s0_slip(constitutive_Nmats))           ; constitutive_s0_slip=0.0_pReal
allocate(constitutive_gdot0_slip(constitutive_Nmats))        ; constitutive_gdot0_slip=0.0_pReal
allocate(constitutive_n_slip(constitutive_Nmats))            ; constitutive_n_slip=0.0_pReal
allocate(constitutive_h0(constitutive_Nmats))                ; constitutive_h0=0.0_pReal
allocate(constitutive_s_sat(constitutive_Nmats))             ; constitutive_s_sat=0.0_pReal
allocate(constitutive_w0(constitutive_Nmats))                ; constitutive_w0=0.0_pReal

part = '_dummy_'
do while (part/='')
  select case (part)
    case ('materials')
	  part = constitutive_parse_materialPart(200)
    case ('textures')
	  part = constitutive_parse_texturePart(200)
	case default
      part = constitutive_parse_unknownPart(200)
  end select
end do
close(200)


do m=1,material_maxN
  material_Cslip_66(:,:,m) = 0.0_pReal
  select case (material_crystal_structure)
    case (1:2) ! cubic structure
      do i=1,3
        do j=1,3
          material_Cslip_66(i,j,m)   = C12
        enddo
        material_Cslip_66(i,i,m)     = C11
        material_Cslip_66(i+3,i+3,m) = C44
      enddo
    case (3)   ! hcp structure MISSING correct
      do i=1,3
        do j=1,3
          material_Cslip_66(i,j,m)   = C12
        enddo
        material_Cslip_66(i,i,m)     = C11
        material_Cslip_66(i+3,i+3,m) = C44
      enddo
  end select
  material_Cslip_3333(:,:,:,:,m) = math_66to3333(Cslip_66(:,:,m))   
end do

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

! MISSING some consistency checks may be..?

return
100 call IO_error(110) ! corrupt matarials_textures file
end subroutine


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
real(pReal) tau_slip(constitutive_Nslip(matID))
real(pReal) tauc_slip(constitutive_Nslip(matID))
real(pReal) gdot_slip(constitutive_Nslip(matID))
real(pReal) dgdot_dtaucslip(constitutive_Nslip(matID))

!* Iteration over the systems 
do i=1,constitutive_Nslip(matID)
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
real(pReal) tauc_slip(constitutive_Nslip(matID))
real(pReal) gdot_slip(constitutive_Nslip(matID))
real(pReal) dtauc_slip(constitutive_Nslip(matID))
real(pReal) self_hardening(constitutive_Nslip(matID))

!* Self-Hardening of each system
do i=1,constitutive_Nslip(matID)
   self_hardening(i)=constitutive_h0(matID)*(1.0_pReal-tauc_slip(i)/constitutive_s_sat(matID))**constitutive_w0(matID)*abs(gdot_slip(i))
enddo

!* Hardening for all systems
i=constitutive_Nslip(matID)
j=constitutive_crystal_structure(matID)
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
real(pReal) tau_slip(constitutive_Nslip(matID))
real(pReal) tauc_slip_new(constitutive_Nslip(matID))
real(pReal) gdot_slip(constitutive_Nslip(matID))
 
!* Calculation of Lp
Lp=0.0_pReal
do i=1,constitutive_Nslip(matID)
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
