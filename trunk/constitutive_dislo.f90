
!************************************
!*      Module: CONSTITUTIVE        *
!************************************
!* contains:                        *
!* - constitutive equations         *
!* - parameters definition          *
!* - orientations                   *
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
!* Physical parameter, attack_frequency != Debye frequency
real(pReal), parameter :: attack_frequency = 1.0e10_pReal  
!* Physical parameter, Boltzman constant in J/Kelvin
real(pReal), parameter :: Kb = 1.38e-23_pReal

!*************************************
!* Definition of material properties *
!*************************************
!* Number of materials
integer(pInt) material_maxN
!* Crystal structure and number of selected slip or twin systems per material
integer(pInt), dimension(:)        , allocatable :: material_CrystalStructure
integer(pInt), dimension(:)        , allocatable :: material_Nslip
integer(pInt), dimension(:)        , allocatable :: material_Ntwin
!* Maximum number of selected slip or twin systems over materials
integer(pInt) material_maxNslip
integer(pInt) material_maxNtwin
!* Elastic constants and matrices
real(pReal), dimension(:)          , allocatable :: material_C11
real(pReal), dimension(:)          , allocatable :: material_C12
real(pReal), dimension(:)          , allocatable :: material_C13
real(pReal), dimension(:)          , allocatable :: material_C33
real(pReal), dimension(:)          , allocatable :: material_C44
real(pReal), dimension(:)          , allocatable :: material_Gmod
real(pReal), dimension(:,:,:)      , allocatable :: material_Cslip_66
real(pReal), dimension(:,:,:,:,:)  , allocatable :: material_Cslip_3333
real(preal), dimension(:,:,:,:)    , allocatable :: material_Ctwin_66
real(pReal), dimension(:,:,:,:,:,:), allocatable :: material_Ctwin_3333
!* Visco-plastic material parameters
real(pReal), dimension(:)          , allocatable :: material_rho0
real(pReal), dimension(:)          , allocatable :: material_bg
real(pReal), dimension(:)          , allocatable :: material_Qedge
real(pReal), dimension(:)          , allocatable :: material_tau0
real(pReal), dimension(:)          , allocatable :: material_c1
real(pReal), dimension(:)          , allocatable :: material_c2
real(pReal), dimension(:)          , allocatable :: material_c3
real(pReal), dimension(:)          , allocatable :: material_c4
real(pReal), dimension(:)          , allocatable :: material_c5
real(pReal), dimension(:,:)        , allocatable :: material_SlipIntCoeff

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
real(pReal), dimension(:)      , allocatable :: constitutive_passing_stress
real(pReal), dimension(:)      , allocatable :: constitutive_jump_width
real(pReal), dimension(:)      , allocatable :: constitutive_activation_volume
real(pReal), dimension(:)      , allocatable :: constitutive_rho_m
real(pReal), dimension(:)      , allocatable :: constitutive_rho_f
real(pReal), dimension(:)      , allocatable :: constitutive_rho_p
real(pReal), dimension(:)      , allocatable :: constitutive_g0_slip

!************************************
!*      Interaction matrices        *
!************************************
real(pReal), dimension(:,:,:), allocatable :: constitutive_Pforest
real(pReal), dimension(:,:,:), allocatable :: constitutive_Pparallel

!************************************
!*             Results              *
!************************************
integer(pInt) constitutive_maxNresults
integer(pInt), dimension(:,:,:), allocatable :: constitutive_Nresults
real(pReal), dimension(:,:,:,:), allocatable :: constitutive_results



CONTAINS
!****************************************
!* - constitutive_Init
!* - constitutive_CountSections
!* - constitutive_Parse_UnknownPart
!* - constitutive_Parse_MaterialPart
!* - constitutive_Parse_TexturePart
!* - constitutive_Parse_MatTexDat
!* - constitutive_Assignment
!* - constitutive_HomogenizedC
!* - constitutive_Microstructure
!* - constitutive_LpAndItsTangent
!* - consistutive_DotState
!****************************************


subroutine constitutive_Init()
!**************************************
!*      Module initialization         *
!**************************************
call constitutive_Parse_MatTexDat(mattexFile)
call constitutive_Assignment()
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
integer(pInt) file,count
integer(pInt), dimension(3) :: positions

count=0
part=''

do
   read(file,'(a80)',END=100) line
   positions=IO_stringPos(line,1)
   tag=IO_lc(IO_stringValue(line,positions,1))
   if (tag(1:1)=='#' .OR. positions(1)==0) then  ! skip comment and empty lines
      cycle
   elseif (tag(1:1)=='<'.AND.tag(len_trim(tag):len_trim(tag))=='>') then
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
integer(pInt) file,section
integer(pInt), dimension(3) :: positions

constitutive_assignNGaussAndFiber=''
section = 0_pInt

do
   read(file,'(a80)',END=100) line
   positions=IO_stringPos(line,1)
   tag=IO_lc(IO_stringValue(line,positions,1))
   if (tag(1:1)=='#' .OR. positions(1)==0) then  ! skip comment and empty lines
      cycle
   elseif (tag(1:1)=='<'.AND.tag(len_trim(tag):len_trim(tag))=='>') then
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
   if (tag(1:1)=='#' .OR. positions(1)==0) then  ! skip comment and empty lines
      cycle
   elseif (tag(1:1)=='<'.AND.tag(len_trim(tag):len_trim(tag))=='>') then
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
use prec, only: pInt,pReal
use IO
implicit none

!* Definition of variables
character(len=80) line,tag
integer(pInt), parameter :: maxNchunks = 7
integer(pInt) i,file,section
integer(pInt), dimension(1+2*maxNchunks) :: positions

section = 0
constitutive_parse_materialPart = ''

do while(.true.)
   read(file,'(a80)',END=100) line
   positions=IO_stringPos(line,maxNchunks) ! parse leading chunks
   tag=IO_lc(IO_stringValue(line,positions,1))
   if (tag(1:1)=='#' .OR. positions(1)==0) then  ! skip comment and empty lines
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
	     case ('ntwin')
		      material_Ntwin(section)=IO_intValue(line,positions,2)
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
         case ('rho0') 
              material_rho0(section)=IO_floatValue(line,positions,2)
	     case ('interaction_coefficients') 
		      do i=1,6
              material_SlipIntCoeff(i,section)=IO_floatValue(line,positions,i+1)
			  enddo
		 case ('burgers') 
              material_bg(section)=IO_floatValue(line,positions,2)	  
		 case ('Qedge') 
              material_Qedge(section)=IO_floatValue(line,positions,2)
		 case ('tau0')
              material_tau0(section)=IO_floatValue(line,positions,2)
		 case ('c1')
              material_c1(section)=IO_floatValue(line,positions,2)
		 case ('c2')
              material_c2(section)=IO_floatValue(line,positions,2)
		 case ('c3')
              material_c3(section)=IO_floatValue(line,positions,2)
		 case ('c4')
              material_c4(section)=IO_floatValue(line,positions,2)
		 case ('c5')
              material_c5(section)=IO_floatValue(line,positions,2)
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
integer(pInt) file,section,gaussCount,fiberCount,i
integer(pInt), dimension(1+2*maxNchunks) :: positions

section = 0
gaussCount = 0
fiberCount = 0
constitutive_parse_texturePart = ''

do while(.true.)
   read(file,'(a80)',END=100) line
   positions=IO_stringPos(line,maxNchunks)  ! parse leading chunks
   tag=IO_lc(IO_stringValue(line,positions,1))
   if (tag(1:1)=='#' .OR. positions(1)==0) then  ! skip comment and empty lines
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
				     texture_Gauss(4,gaussCount,section)=IO_floatValue(line,positions,i+1)*inRad
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
				     texture_fiber(5,fiberCount,section)=IO_floatValue(line,positions,i+1)*inRad
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
use IO, only: IO_error, IO_open_file
use math, only: math_Mandel3333to66, math_Voigt66to3333
use crystal, only: crystal_MaxMaxNslipOfStructure,crystal_MaxMaxNtwinOfStructure
implicit none

!* Definition of variables
character(len=*) filename
character(len=80) part,formerPart
integer(pInt) sectionCount,i,j,k, fileunit

! set fileunit
fileunit=200
!-----------------------------
!* First reading: number of materials and textures
!-----------------------------
!* determine material_maxN and texture_maxN from last respective parts
if(IO_open_file(fileunit,filename)==.false.) goto 100
part = '_dummy_'
do while (part/='')
   formerPart = part
   call constitutive_CountSections(fileunit,sectionCount,part)
   select case (formerPart)
   case ('materials')
        material_maxN = sectionCount
   case ('textures')
        texture_maxN = sectionCount
   end select
enddo
!* Array allocation
allocate(material_CrystalStructure(material_maxN))		                            ; material_CrystalStructure=0_pInt
allocate(material_Nslip(material_maxN))					                            ; material_Nslip=0_pInt
allocate(material_C11(material_maxN))					                            ; material_C11=0.0_pReal
allocate(material_C12(material_maxN))					                            ; material_C12=0.0_pReal
allocate(material_C13(material_maxN))					                            ; material_C13=0.0_pReal
allocate(material_C33(material_maxN))					                            ; material_C33=0.0_pReal
allocate(material_C44(material_maxN))					                            ; material_C44=0.0_pReal
allocate(material_Gmod(material_maxN))					                            ; material_Gmod=0.0_pReal
allocate(material_Cslip_66(6,6,material_maxN))                                      ; material_Cslip_66=0.0_pReal
allocate(material_Cslip_3333(3,3,3,3,material_maxN))                                ; material_Cslip_3333=0.0_pReal
allocate(material_Ctwin_66(6,6,crystal_MaxMaxNtwinOfStructure,material_maxN))       ; material_Ctwin_66=0.0_pReal
allocate(material_Ctwin_3333(3,3,3,3,crystal_MaxMaxNtwinOfStructure,material_maxN)) ; material_Ctwin_3333=0.0_pReal
allocate(material_rho0(material_maxN))				                                ; material_rho0=0.0_pReal
allocate(material_SlipIntCoeff(crystal_MaxMaxNslipOfStructure,material_maxN))       ; material_SlipIntCoeff=0.0_pReal 
allocate(material_bg(material_maxN))			                                    ; material_bg=0.0_pReal
allocate(material_Qedge(material_maxN))				                                ; material_Qedge=0.0_pReal
allocate(material_tau0(material_maxN))				                                ; material_tau0=0.0_pReal
allocate(material_c1(material_maxN))				                                ; material_c1=0.0_pReal
allocate(material_c2(material_maxN))                                                ; material_c2=0.0_pReal
allocate(material_c3(material_maxN))                                                ; material_c3=0.0_pReal
allocate(material_c4(material_maxN))                                                ; material_c4=0.0_pReal
allocate(material_c5(material_maxN))                                                ; material_c5=0.0_pReal
allocate(texture_ODFfile(texture_maxN))                                             ; texture_ODFfile=''
allocate(texture_Ngrains(texture_maxN))                                             ; texture_Ngrains=0_pInt
allocate(texture_symmetry(texture_maxN))                                            ; texture_symmetry=''
allocate(texture_NGauss(texture_maxN))                                              ; texture_NGauss=0_pInt
allocate(texture_NFiber(texture_maxN))                                              ; texture_NFiber=0_pInt
allocate(texture_NRandom(texture_maxN))                                             ; texture_NRandom=0_pInt

!-----------------------------
!* Second reading: number of Gauss and Fiber
!-----------------------------
rewind(fileunit)
part = '_dummy_'
do while (part/='')
   select case (part)
   case ('textures')
        part = constitutive_assignNGaussAndFiber(fileunit)
   case default
        part = constitutive_Parse_UnknownPart(fileunit)
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
rewind(fileunit)
part='_dummy_'
do while (part/='')
   select case (part)
   case ('materials')
	    part=constitutive_Parse_MaterialPart(fileunit)
   case ('textures')
	    part=constitutive_Parse_TexturePart(fileunit)
   case default
        part=constitutive_Parse_UnknownPart(fileunit)
   end select
enddo
close(fileunit)

!* Construction of the elasticity matrices
do i=1,material_maxN
   select case (material_CrystalStructure(i))
   case(1:2) ! cubic(s) 
       material_Gmod(i)=material_C44(i)
       forall(k=1:3)
          forall(j=1:3)
             material_Cslip_66(k,j,i)=material_C12(i)
          endforall
          material_Cslip_66(k,k,i)=material_C11(i)
          material_Cslip_66(k+3,k+3,i)=material_C44(i)
       endforall
   case(3)   ! hcp
        material_Gmod(i)=material_C44(i)
        !* MISSING: Warning message: C44 in hexagonal structures?
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
   material_Cslip_3333(:,:,:,:,i) = math_Voigt66to3333(material_Cslip_66(:,:,i))
   material_Cslip_66(:,:,i) = math_Mandel3333to66(material_Cslip_3333(:,:,:,:,i))
enddo

! MISSING some consistency checks may be..?
! if ODFfile present then set NGauss NFiber =0
return
100 call IO_error(200) ! corrupt materials_textures file
end subroutine


subroutine constitutive_Assignment()
!*********************************************************************
!* This subroutine assign material parameters according to ipc,ip,el *
!*********************************************************************
use prec, only: pReal,pInt
use math, only: math_sampleGaussOri,math_sampleFiberOri,math_sampleRandomOri,math_symmetricEulers,math_EulerToR,&
                math_Mandel3333to66
use mesh, only: mesh_NcpElems,FE_Nips,FE_mapElemtype,mesh_maxNips,mesh_element
use IO,   only: IO_hybridIA
use crystal, only: crystal_SlipIntType,crystal_sn,crystal_st,crystal_Qtwin
implicit none

!* Definition of variables
integer(pInt) e,g,i,j,k,l,m,n,o,p,q,r,s
integer(pInt) matID,texID
real(pReal) x,y
integer(pInt), dimension(:,:,:), allocatable :: hybridIA_population
integer(pInt), dimension(texture_maxN) :: Ncomponents,Nsym,multiplicity,ODFmap,sampleCount
real(pReal), dimension(3,4*(1+texture_maxNGauss+texture_maxNfiber)) :: Euler
real(pReal), dimension(4*(1+texture_maxNGauss+texture_maxNfiber)) :: texVolfrac
real(pReal), dimension(texture_maxN) :: sumVolfrac

! process textures
o = 0_pInt       ! ODF counter
ODFmap = 0_pInt  ! blank mapping
sampleCount = 0_pInt  ! count orientations assigned per texture

do texID=1,texture_maxN
   select case (texture_symmetry(texID))          ! set symmetry factor
   case ('orthotropic')
     Nsym(texID) = 4_pInt
   case ('monoclinic')
     Nsym(texID) = 2_pInt
   case default
     Nsym(texID) = 1_pInt
   end select
   if (texture_ODFfile(texID)=='') then           ! texture components
      sumVolfrac(texID) = sum(texture_gauss(5,:,texID))+sum(texture_fiber(6,:,texID))
      if (sumVolfrac(texID)<1.0_pReal) texture_NRandom(texID) = 1_pInt  ! check whether random component missing
      Ncomponents(texID) = texture_NGauss(texID)+texture_NFiber(texID)+texture_NRandom(texID)
   else      ! hybrid IA
      o = o+1
	  ODFmap(texID) = o             ! remember mapping
      Ncomponents(texID) = 1_pInt   ! single "component"
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
material_maxNslip       = maxval(material_Nslip)	! max of slip systems among materials present
material_maxNtwin       = maxval(material_Ntwin)    ! max of twin systems among materials present 
constitutive_maxNstatevars = maxval(material_Nslip) + maxval(material_Ntwin)
constitutive_maxNresults = 1_pInt

!* calc texture_totalNgrains
allocate(texture_totalNgrains(texture_maxN)) ; texture_totalNgrains=0_pInt
do i=1,mesh_NcpElems
   texID = mesh_element(4,i)
   texture_totalNgrains(texID) = texture_totalNgrains(texID) + FE_Nips(FE_mapElemtype(mesh_element(2,i)))*texture_Ngrains(texID)
enddo

! generate hybridIA samplings for ODFfile textures to later draw from these populations
allocate(hybridIA_population(3,maxval(texture_totalNgrains/Nsym,ODFmap /= 0),o))
do texID = 1,texture_maxN
   if (ODFmap(texID) > 0) &
      hybridIA_population(:,:,ODFmap(texID)) = IO_hybridIA(texture_totalNgrains(texID)/Nsym(texID),texture_ODFfile(texID))
enddo

!* Array allocation
allocate(constitutive_Ngrains(mesh_maxNips,mesh_NcpElems)) ; constitutive_Ngrains=0_pInt
allocate(constitutive_matID(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_matID=0_pInt
allocate(constitutive_texID(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_texID=0_pInt
allocate(constitutive_MatVolFrac(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_MatVolFrac=0.0_pReal
allocate(constitutive_TexVolFrac(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_TexVolFrac=0.0_pReal
allocate(constitutive_EulerAngles(3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_EulerAngles=0.0_pReal
allocate(constitutive_Nresults(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_Nresults=0_pInt
allocate(constitutive_results(constitutive_maxNresults,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
constitutive_results=0.0_pReal
allocate(constitutive_Nstatevars(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_Nstatevars=0_pInt
allocate(constitutive_state_old(constitutive_maxNstatevars,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
constitutive_state_old=0.0_pReal
allocate(constitutive_state_new(constitutive_maxNstatevars,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
constitutive_state_new=0.0_pReal
allocate(constitutive_Pforest(material_maxNslip,material_maxNslip,material_maxN)) 
constitutive_Pforest=0.0_pReal
allocate(constitutive_Pparallel(material_maxNslip,material_maxNslip,material_maxN))
constitutive_Pparallel=0.0_pReal
allocate(constitutive_rho_p(material_maxNslip))             ; constitutive_rho_p=0.0_pReal
allocate(constitutive_rho_f(material_maxNslip))             ; constitutive_rho_f=0.0_pReal
allocate(constitutive_rho_m(material_maxNslip))             ; constitutive_rho_m=0.0_pReal
allocate(constitutive_passing_stress(material_maxNslip))    ; constitutive_passing_stress=0.0_pReal
allocate(constitutive_jump_width(material_maxNslip))        ; constitutive_jump_width=0.0_pReal
allocate(constitutive_activation_volume(material_maxNslip)) ; constitutive_activation_volume=0.0_pReal
allocate(constitutive_g0_slip(material_maxNslip))           ; constitutive_g0_slip=0.0_pReal

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
            constitutive_matID(g,i,e) = matID          ! copy matID of element
		    constitutive_texID(g,i,e) = texID          ! copy texID of element
		    constitutive_MatVolFrac(g,i,e) = 1.0_pReal ! singular material (so far)
		    constitutive_TexVolFrac(g,i,e) = texVolfrac(s)/multiplicity(texID)/Nsym(texID)
		    constitutive_Nstatevars(g,i,e) = material_Nslip(matID) ! number of state variables (i.e. tau_c of each slip system)
		    constitutive_Nresults(g,i,e)   = 0         ! number of constitutive results
		    constitutive_EulerAngles(:,g,i,e) = Euler(:,s)    ! store initial orientation
            forall (l=1:constitutive_Nstatevars(g,i,e))  ! initialize state variables
               constitutive_state_old(l,g,i,e) = material_rho0(matID)
               constitutive_state_new(l,g,i,e) = material_rho0(matID)
            end forall
		 enddo  ! components
         sampleCount(texID) = sampleCount(texID)+1     ! next member of hybrid IA population
	  enddo  ! multiplicity
   enddo ! ip
enddo ! cp_element

!* Construction of the rotated elasticity matrices for twinning
do i=1,material_maxN
   do j=1,material_Ntwin(i)
      do k=1,3
      do l=1,3
      do m=1,3
      do n=1,3
         material_Ctwin_3333(k,l,m,n,j,i)=0.0_pReal
         do p=1,3
	     do q=1,3
	     do r=1,3
	     do s=1,3
	        material_Ctwin_3333(k,l,m,n,j,i)=material_Ctwin_3333(k,l,m,n,j,i)+material_Cslip_3333(p,q,r,s,i)*&
		                                     crystal_Qtwin(k,p,j,material_CrystalStructure(i))*&
											 crystal_Qtwin(l,q,j,material_CrystalStructure(i))*&
											 crystal_Qtwin(m,r,j,material_CrystalStructure(i))*&
											 crystal_Qtwin(n,s,j,material_CrystalStructure(i))
	     enddo
	     enddo
	     enddo
	     enddo
      enddo
      enddo
      enddo
      enddo
      material_Ctwin_66(:,:,j,i) = math_Mandel3333to66(material_Ctwin_3333(:,:,:,:,j,i))
   enddo
enddo

!* Construction of the hardening matrices
do i=1,material_maxN
!* Iteration over the systems
   do j=1,material_Nslip(i)
   do k=1,material_Nslip(i)
!* Projection of the dislocation *
	  x=dot_product(crystal_sn(:,j,i),crystal_st(:,k,i))
	  y=1.0_pReal-x**(2.0_pReal)
!* Interaction matrix *
      constitutive_Pforest(j,k,i)=abs(x)*material_SlipIntCoeff(crystal_SlipIntType(j,k,i),i)
	  if (y>0.0_pReal) then
	     constitutive_Pparallel(j,k,i)=sqrt(y)*material_SlipIntCoeff(crystal_SlipIntType(j,k,i),i)
	  else
	     constitutive_Pparallel(j,k,i)=0.0_pReal
	  endif
   enddo
   enddo
enddo

end subroutine


function constitutive_HomogenizedC(state,ipc,ip,el)
!*********************************************************************
!* This function returns the homogenized elacticity matrix           *
!* INPUT:                                                            *
!*  - state           : state variables                              *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
use prec, only: pReal,pInt
implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
real(pReal), dimension(6,6) :: constitutive_homogenizedC
real(pReal), dimension(constitutive_Nstatevars(ipc,ip,el)) :: state

!* Homogenization scheme
constitutive_homogenizedC=material_Cslip_66(:,:,constitutive_matID(ipc,ip,el))

return
end function


subroutine constitutive_Microstructure(state,Tp,ipc,ip,el)
!*********************************************************************
!* This function calculates from state needed variables              *
!* INPUT:                                                            *
!*  - state           : state variables                              *
!*  - Tp              : temperature                                  *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
use prec, only: pReal,pInt
implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
integer(pInt) matID,i
real(pReal) Tp
real(pReal), dimension(constitutive_Nstatevars(ipc,ip,el)) :: state

!* Get the material-ID from the triplet(ipc,ip,el)
matID = constitutive_matID(ipc,ip,el)

!* Quantities derivated from state
constitutive_rho_f=matmul(constitutive_Pforest  (1:material_Nslip(matID),1:material_Nslip(matID),matID),state)
constitutive_rho_p=matmul(constitutive_Pparallel(1:material_Nslip(matID),1:material_Nslip(matID),matID),state)	
do i=1,material_Nslip(matID)
   constitutive_passing_stress(i)=material_tau0(matID)+material_c1(matID)*material_Gmod(matID)*material_bg(matID)*&
                                  sqrt(constitutive_rho_p(i))
   constitutive_jump_width(i)=material_c2(matID)/sqrt(constitutive_rho_f(i))
   constitutive_activation_volume(i)=material_c3(matID)*constitutive_jump_width(i)*material_bg(matID)**2.0_pReal 
   constitutive_rho_m(i)=(2.0_pReal*Kb*Tp*sqrt(constitutive_rho_p(i)))/&
                         (material_c1(matID)*material_c3(matID)*material_Gmod(matID)*constitutive_jump_width(i)*&
						 material_bg(matID)**3.0_pReal) 
   constitutive_g0_slip(i)=constitutive_rho_m(i)*material_bg(matID)*attack_frequency*constitutive_jump_width(i)*&
                           exp(-(material_Qedge(matID)+constitutive_passing_stress(i)*constitutive_activation_volume(i))/&
						   (Kb*Tp))
enddo

end subroutine


subroutine constitutive_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,state,Tp,ipc,ip,el)
!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating the velocity gradient                                 *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - state           : current microstructure                       *
!*  - Tp              : temperature                                  *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - Lp              : plastic velocity gradient                    *
!*  - dLp_dTstar      : derivative of Lp (4th-order tensor)          *
!*********************************************************************
use prec, only: pReal,pInt
use crystal, only: crystal_Sslip,crystal_Sslip_v
implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
integer(pInt) matID,i,k,l,m,n
real(pReal) Tp
real(pReal), dimension(6) :: Tstar_v
real(pReal), dimension(3,3) :: Lp
real(pReal), dimension(3,3,3,3) :: dLp_dTstar
real(pReal), dimension(constitutive_Nstatevars(ipc,ip,el)) :: state,gdot_slip,dgdot_dtauslip,tau_slip

!* Get the material-ID from the triplet(ipc,ip,el)
matID = constitutive_matID(ipc,ip,el)

!* Calculation of Lp
Lp = 0.0_pReal
do i=1,material_Nslip(matID)
   tau_slip(i)=dot_product(Tstar_v,crystal_Sslip_v(:,i,material_CrystalStructure(matID)))
   gdot_slip(i)=constitutive_g0_slip(i)*sinh((abs(tau_slip(i))*constitutive_activation_volume(i))/(Kb*Tp))*&
                sign(1.0_pReal,tau_slip(i))   
   Lp=Lp+gdot_slip(i)*crystal_Sslip(:,:,i,material_CrystalStructure(matID))
enddo

!* Lp twin

!* Calculation of the tangent of Lp
dLp_dTstar=0.0_pReal
do i=1,material_Nslip(matID)
   dgdot_dtauslip(i)=((constitutive_g0_slip(i)*constitutive_activation_volume(i))/(Kb*Tp))*&
                     cosh((abs(tau_slip(i))*constitutive_activation_volume(i))/(Kb*Tp))  
   forall (k=1:3,l=1:3,m=1:3,n=1:3)
          dLp_dTstar(k,l,m,n) = dLp_dTstar(k,l,m,n)+ &
                                dgdot_dtauslip(i)*crystal_Sslip(k,l,i,material_CrystalStructure(matID))* &
		                        (crystal_Sslip(m,n,i,material_CrystalStructure(matID))+ &
		                         crystal_Sslip(n,m,i,material_CrystalStructure(matID)))/2.0_pReal ! force m,n symmetry
   endforall
enddo

return
end subroutine


function constitutive_dotState(Tstar_v,state,Tp,ipc,ip,el)
!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating rate of change of microstructure                      *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - state           : current microstructure                       *
!*  - Tp              : temperature                                  *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - constitutive_DotState : evolution of state variable            *
!*********************************************************************
use prec, only: pReal,pInt
use crystal, only: crystal_Sslip_v
implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
integer(pInt) matID,i
real(pReal) Tp,tau_slip,gdot_slip,lock,recovery
real(pReal), dimension(6) :: Tstar_v
real(pReal), dimension(constitutive_Nstatevars(ipc,ip,el)) :: constitutive_dotState,state

!* Get the material-ID from the triplet(ipc,ip,el)
matID = constitutive_matID(ipc,ip,el)

!* Hardening of each system
do i=1,material_Nslip(matID)
   tau_slip = dot_product(Tstar_v,crystal_Sslip_v(:,i,material_CrystalStructure(matID)))
   gdot_slip = constitutive_g0_slip(i)*sinh((abs(tau_slip)*constitutive_activation_volume(i))/(Kb*Tp))*&
               sign(1.0_pReal,tau_slip)
   lock=(material_c4(matID)*sqrt(constitutive_rho_f(i))*abs(gdot_slip))/material_bg(matID)
   recovery=material_c5(matID)*state(i)*abs(gdot_slip)
   constitutive_dotState(i)=lock-recovery
enddo

return
end function


function constitutive_post_results(Tstar_v,state,dt,Tp,ipc,ip,el)
!*********************************************************************
!* return array of constitutive results                              *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - state           : current microstructure                       *
!*  - dt              : current time increment                       *
!*  - Tp              : temperature                                  *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
use prec, only: pReal,pInt
use crystal, only: crystal_Sslip_v
implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
integer(pInt) matID,i
real(pReal) dt,Tp,tau_slip
real(pReal), dimension(6) :: Tstar_v
real(pReal), dimension(constitutive_Nstatevars(ipc,ip,el)) :: state
real(pReal), dimension(constitutive_Nresults(ipc,ip,el))   :: constitutive_post_results

!* Get the material-ID from the triplet(ipc,ip,el)
matID = constitutive_matID(ipc,ip,el)

if(constitutive_Nresults(ipc,ip,el)==0) return

do i=1,material_Nslip(matID)
   constitutive_post_results(i) = state(i)
   tau_slip=dot_product(Tstar_v,crystal_Sslip_v(:,i,material_CrystalStructure(matID)))
   constitutive_post_results(i+material_Nslip(matID)) = &
     dt*constitutive_g0_slip(i)*sinh((abs(tau_slip)*constitutive_activation_volume(i))/(Kb*Tp))*&
	 sign(1.0_pReal,tau_slip)
enddo
return

end function

END MODULE