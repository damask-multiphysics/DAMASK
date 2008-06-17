
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
!* Physical parameter, Boltzmann constant in J/Kelvin
real(pReal), parameter :: kB = 1.38e-23_pReal
!* Physical parameter, Avogadro number in 1/mol
real(pReal), parameter :: avogadro = 6.022e23_pReal 
!* Physical parameter, Gaz constant in J.mol/Kelvin
real(pReal), parameter :: Rgaz = 8.314_pReal    

!*************************************
!* Definition of material properties *
!*************************************
!* Number of materials
integer(pInt) material_maxN
!* Lattice structure and number of selected slip or twin systems per material
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
real(preal), dimension(:,:,:,:)    , allocatable :: material_Ctwin_66
!* Visco-plastic material parameters
real(pReal) material_sfe
real(pReal), dimension(:)          , allocatable :: material_rho0
real(pReal), dimension(:)          , allocatable :: material_bg
real(pReal), dimension(:)          , allocatable :: material_Qedge
real(pReal), dimension(:)          , allocatable :: material_GrainSize
real(pReal), dimension(:)          , allocatable :: material_StackSize
real(pReal), dimension(:)          , allocatable :: material_ActivationLength
real(pReal), dimension(:)          , allocatable :: material_TwinSaturation
real(pReal), dimension(:)          , allocatable :: material_SiteScaling
real(pReal), dimension(:)          , allocatable :: material_c1
real(pReal), dimension(:)          , allocatable :: material_c2
real(pReal), dimension(:)          , allocatable :: material_c3
real(pReal), dimension(:)          , allocatable :: material_c4
real(pReal), dimension(:)          , allocatable :: material_c5
real(pReal), dimension(:)          , allocatable :: material_c6
real(pReal), dimension(:)          , allocatable :: material_c7
real(pReal), dimension(:)          , allocatable :: material_q1
real(pReal), dimension(:)          , allocatable :: material_q2
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
!*        Kinetics variables        *
!************************************
real(pReal), dimension(:)      , allocatable :: constitutive_tau_slip
real(pReal), dimension(:)      , allocatable :: constitutive_tau_twin
real(pReal), dimension(:)      , allocatable :: constitutive_gdot_slip
real(pReal), dimension(:)      , allocatable :: constitutive_fdot_twin
real(pReal), dimension(:)      , allocatable :: constitutive_dgdot_dtauslip
real(pReal), dimension(:)      , allocatable :: constitutive_dfdot_dtautwin
real(pReal), dimension(:,:)    , allocatable :: constitutive_dfdot_dtauslip
real(pReal), dimension(:)      , allocatable :: constitutive_locks
real(pReal), dimension(:)      , allocatable :: constitutive_grainboundaries
real(pReal), dimension(:)      , allocatable :: constitutive_twinboundaries
real(pReal), dimension(:)      , allocatable :: constitutive_recovery

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
real(pReal), dimension(:)      , allocatable :: constitutive_twin_volume
real(pReal), dimension(:)      , allocatable :: constitutive_inv_intertwin_len_s
real(pReal), dimension(:)      , allocatable :: constitutive_inv_intertwin_len_t
real(pReal), dimension(:)      , allocatable :: constitutive_twin_mfp

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
!$OMP CRITICAL (write2out)
         select case(tag)
	     case ('lattice_structure')
              material_CrystalStructure(section)=IO_intValue(line,positions,2)
			  write(6,*) 'lattice_structure', material_CrystalStructure(section)
	     case ('nslip')
		      material_Nslip(section)=IO_intValue(line,positions,2)
			  write(6,*) 'nslip', material_Nslip(section)
	     case ('ntwin')
		      material_Ntwin(section)=IO_intValue(line,positions,2)
			  write(6,*) 'ntwin', material_Ntwin(section)
		 case ('c11')
              material_C11(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'c11', material_C11(section)
		 case ('c12')
              material_C12(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'c12', material_C12(section)
		 case ('c13')
              material_C13(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'c13', material_C13(section)
		 case ('c33')
              material_C33(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'c33', material_C33(section)
		 case ('c44')
              material_C44(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'c44', material_C44(section)
         case ('rho0') 
              material_rho0(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'rho0', material_rho0(section)
	     case ('interaction_coefficients') 
		      do i=1,6
              material_SlipIntCoeff(i,section)=IO_floatValue(line,positions,i+1)
              write(6,*) 'interaction_coefficients', material_SlipIntCoeff(i,section)
			  enddo
		 case ('burgers') 
              material_bg(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'burgers', material_bg(section) 
		 case ('qedge') 
              material_Qedge(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'Qedge', material_Qedge(section)
	     case ('grain_size')
              material_GrainSize(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'grain_size', material_GrainSize(section)
	     case ('stack_size')
              material_StackSize(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'stack_size', material_StackSize(section)
		 case ('f_sat')
              material_TwinSaturation(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'twin saturation', material_TwinSaturation(section)
         case ('site_scaling')
              material_SiteScaling(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'site_scaling', material_SiteScaling(section)
		 case ('c1')
              material_c1(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'c1', material_c1(section)
		 case ('c2')
              material_c2(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'c2', material_c2(section)
		 case ('c3')
              material_c3(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'c3', material_c3(section)
		 case ('c4')
              material_c4(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'c4', material_c4(section)
		 case ('c5')
              material_c5(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'c5', material_c5(section)
		 case ('c6')
              material_c6(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'c6', material_c6(section)
	     case ('c7')
              material_c7(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'c7', material_c7(section)
		 case ('q1')
              material_q1(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'q1', material_q1(section)
	     case ('q2')
              material_q2(section)=IO_floatValue(line,positions,2)
			  write(6,*) 'q2', material_q2(section)
         end select
!$OMP END CRITICAL (write2out)
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
use lattice, only: lattice_MaxMaxNslipOfStructure,lattice_MaxMaxNtwinOfStructure
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
allocate(material_Ntwin(material_maxN))                                             ; material_Ntwin=0_pInt
allocate(material_C11(material_maxN))					                            ; material_C11=0.0_pReal
allocate(material_C12(material_maxN))					                            ; material_C12=0.0_pReal
allocate(material_C13(material_maxN))					                            ; material_C13=0.0_pReal
allocate(material_C33(material_maxN))					                            ; material_C33=0.0_pReal
allocate(material_C44(material_maxN))					                            ; material_C44=0.0_pReal
allocate(material_Gmod(material_maxN))					                            ; material_Gmod=0.0_pReal
allocate(material_Cslip_66(6,6,material_maxN))                                      ; material_Cslip_66=0.0_pReal
allocate(material_Ctwin_66(6,6,lattice_MaxMaxNtwinOfStructure,material_maxN))       ; material_Ctwin_66=0.0_pReal
allocate(material_rho0(material_maxN))				                                ; material_rho0=0.0_pReal
allocate(material_SlipIntCoeff(lattice_MaxMaxNslipOfStructure,material_maxN))       ; material_SlipIntCoeff=0.0_pReal 
allocate(material_bg(material_maxN))			                                    ; material_bg=0.0_pReal
allocate(material_Qedge(material_maxN))				                                ; material_Qedge=0.0_pReal
allocate(material_GrainSize(material_maxN))				                            ; material_GrainSize=0.0_pReal
allocate(material_StackSize(material_maxN))				                            ; material_StackSize=0.0_pReal
allocate(material_ActivationLength(material_maxN))                                  ; material_ActivationLength=0.0_pReal
allocate(material_TwinSaturation(material_maxN))                                    ; material_TwinSaturation=0.0_pReal
allocate(material_SiteScaling(material_maxN))                                       ; material_SiteScaling=0.0_pReal
allocate(material_c1(material_maxN))				                                ; material_c1=0.0_pReal
allocate(material_c2(material_maxN))                                                ; material_c2=0.0_pReal
allocate(material_c3(material_maxN))                                                ; material_c3=0.0_pReal
allocate(material_c4(material_maxN))                                                ; material_c4=0.0_pReal
allocate(material_c5(material_maxN))                                                ; material_c5=0.0_pReal
allocate(material_c6(material_maxN))                                                ; material_c6=0.0_pReal
allocate(material_c7(material_maxN))                                                ; material_c7=0.0_pReal
allocate(material_q1(material_maxN))                                                ; material_q1=0.0_pReal
allocate(material_q2(material_maxN))                                                ; material_q2=0.0_pReal
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
   material_Cslip_66(:,:,i) = math_Mandel3333to66(math_Voigt66to3333(material_Cslip_66(:,:,i)))
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
                math_Mandel3333to66,math_Mandel66to3333
use mesh, only: mesh_NcpElems,FE_Nips,mesh_maxNips,mesh_element
use IO,   only: IO_hybridIA
use lattice, only: lattice_SlipIntType,lattice_sn,lattice_st,lattice_Qtwin,lattice_Sslip_v,lattice_Sslip
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
real(pReal), dimension(3,3,3,3) :: C_3333,Ctwin_3333
real(pReal), dimension(3,3) :: Qtwin

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
!$OMP CRITICAL (write2out)
      write (6,*) 'changed Ngrains to',texture_Ngrains(texID),' for texture',texID
!$OMP END CRITICAL (write2out)
   endif
enddo

!* publish globals
constitutive_maxNgrains = maxval(texture_Ngrains)
material_maxNslip       = maxval(material_Nslip)	! max of slip systems among materials present
material_maxNtwin       = maxval(material_Ntwin)    ! max of twin systems among materials present 
constitutive_maxNstatevars = maxval(material_Nslip) + maxval(material_Ntwin)
! -----------------------------------------------------------------------------------------------------------------------
constitutive_maxNresults = 2_pInt
! -----------------------------------------------------------------------------------------------------------------------


!* calc texture_totalNgrains
allocate(texture_totalNgrains(texture_maxN)) ; texture_totalNgrains=0_pInt
do i=1,mesh_NcpElems
   texID = mesh_element(4,i)
   texture_totalNgrains(texID) = texture_totalNgrains(texID) + FE_Nips(mesh_element(2,i))*texture_Ngrains(texID)
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
allocate(constitutive_twin_volume(material_maxNtwin))       ; constitutive_twin_volume=0.0_pReal
allocate(constitutive_inv_intertwin_len_s(material_maxNslip)) ; constitutive_inv_intertwin_len_s=0.0_pReal
allocate(constitutive_inv_intertwin_len_t(material_maxNtwin)) ; constitutive_inv_intertwin_len_t=0.0_pReal
allocate(constitutive_twin_mfp(material_maxNtwin))          ; constitutive_twin_mfp=0.0_pReal
allocate(constitutive_tau_slip(material_maxNslip))          ; constitutive_tau_slip=0.0_pReal
allocate(constitutive_tau_twin(material_maxNtwin))          ; constitutive_tau_twin=0.0_pReal
allocate(constitutive_gdot_slip(material_maxNslip))         ; constitutive_gdot_slip=0.0_pReal
allocate(constitutive_fdot_twin(material_maxNtwin))         ; constitutive_fdot_twin=0.0_pReal
allocate(constitutive_dgdot_dtauslip(material_maxNslip))    ; constitutive_dgdot_dtauslip=0.0_pReal
allocate(constitutive_dfdot_dtautwin(material_maxNtwin))    ; constitutive_dfdot_dtautwin=0.0_pReal
allocate(constitutive_dfdot_dtauslip(material_maxNtwin,material_maxNslip)) ; constitutive_dfdot_dtauslip=0.0_pReal
allocate(constitutive_locks(material_maxNslip))             ; constitutive_locks=0.0_pReal
allocate(constitutive_grainboundaries(material_maxNslip))   ; constitutive_grainboundaries=0.0_pReal
allocate(constitutive_twinboundaries(material_maxNslip))    ; constitutive_twinboundaries=0.0_pReal
allocate(constitutive_recovery(material_maxNslip))          ; constitutive_locks=0.0_pReal

!* Assignment of all grains in all IPs of all cp-elements
do e=1,mesh_NcpElems
   matID=mesh_element(3,e)
   texID=mesh_element(4,e)
   do i=1,FE_Nips(mesh_element(2,e))
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
		    constitutive_Nstatevars(g,i,e) = material_Nslip(matID) + material_Ntwin(matID)! number of state variables (i.e. tau_c of each slip system)
! -----------------------------------------------------------------------------------------------------------------------
		    constitutive_Nresults(g,i,e)   = 2         ! number of constitutive results output by constitutive_post_results
! -----------------------------------------------------------------------------------------------------------------------
		    constitutive_EulerAngles(:,g,i,e) = Euler(:,s)    ! store initial orientation
            forall (l=1:material_Nslip(matID))  ! initialize state variables
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
   C_3333=math_Mandel66to3333(material_Cslip_66(:,:,i))
   do j=1,material_Ntwin(i)
      Qtwin=lattice_Qtwin(:,:,j,material_CrystalStructure(i))
      do k=1,3
      do l=1,3
      do m=1,3
      do n=1,3
         Ctwin_3333(k,l,m,n)=0.0_pReal
         do p=1,3
	     do q=1,3
	     do r=1,3
	     do s=1,3
	        Ctwin_3333(k,l,m,n)=Ctwin_3333(k,l,m,n)+C_3333(p,q,r,s)*Qtwin(k,p)*Qtwin(l,q)*Qtwin(m,r)*Qtwin(n,s)
	     enddo
	     enddo
	     enddo
	     enddo
      enddo
      enddo
      enddo
      enddo
	  !* Mapping back to 66-format of the matrices
      material_Ctwin_66(:,:,j,i) = math_Mandel3333to66(Ctwin_3333)
   enddo
enddo

!* Construction of the hardening matrices
do i=1,material_maxN
!* Iteration over the systems
   do j=1,material_Nslip(i)
   do k=1,material_Nslip(i)
!* Projection of the dislocation *
	  x=dot_product(lattice_sn(:,j,i),lattice_st(:,k,i))
	  y=1.0_pReal-x**(2.0_pReal)
!* Interaction matrix *
      constitutive_Pforest(j,k,i)=abs(x)*material_SlipIntCoeff(lattice_SlipIntType(j,k,i),i)
	  if (y>0.0_pReal) then
	     constitutive_Pparallel(j,k,i)=sqrt(y)*material_SlipIntCoeff(lattice_SlipIntType(j,k,i),i)
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
integer(pInt) matID,i,startIdxTwin
real(pReal), dimension(6,6) :: constitutive_homogenizedC
real(pReal), dimension(constitutive_Nstatevars(ipc,ip,el)) :: state

!* Get the material-ID from the triplet(ipc,ip,el)
matID = constitutive_matID(ipc,ip,el)
startIdxTwin = material_Nslip(matID)

!* Homogenization scheme
constitutive_homogenizedC=(1-sum(state((startIdxTwin+1):(startIdxTwin+material_Ntwin(matID)))))*&
                          material_Cslip_66(:,:,matID)
do i=1,material_Ntwin(matID)
   constitutive_homogenizedC=constitutive_homogenizedC+state(startIdxTwin+i)*material_Ctwin_66(:,:,i,matID) 
enddo

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
use math, only: pi
use lattice, only: lattice_TwinIntType,lattice_SlipTwinIntType
implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
integer(pInt) matID,i,j,startIdxTwin
real(pReal) Tp,Ftwin
real(pReal), dimension(constitutive_Nstatevars(ipc,ip,el)) :: state
real(pReal) x_fe,x_Mn,x_C,beta_mart,Tp_mart,f_mart,beta_aust,Tp_aust,f_aust
real(pReal) deltaG1,deltaG2,deltaG3,deltaG4,deltaG5

!* Get the material-ID from the triplet(ipc,ip,el)
matID = constitutive_matID(ipc,ip,el)
startIdxTwin = material_Nslip(matID)

!* Quantities derived from state - slip
!$OMP CRITICAL (evilmatmul)
constitutive_rho_f=matmul(constitutive_Pforest  (1:material_Nslip(matID),1:material_Nslip(matID),matID),state)
constitutive_rho_p=matmul(constitutive_Pparallel(1:material_Nslip(matID),1:material_Nslip(matID),matID),state)
!$OMP END CRITICAL (evilmatmul)	
do i=1,material_Nslip(matID)
   constitutive_passing_stress(i) = material_c1(matID)*material_Gmod(matID)*material_bg(matID)*&
                                    sqrt(constitutive_rho_p(i))

   constitutive_jump_width(i) = material_c2(matID)/sqrt(constitutive_rho_f(i))

   constitutive_activation_volume(i) = material_c3(matID)*constitutive_jump_width(i)*material_bg(matID)**2.0_pReal 

   constitutive_rho_m(i) = (2.0_pReal*kB*Tp*sqrt(constitutive_rho_p(i)))/&
                           (material_c1(matID)*material_c3(matID)*material_Gmod(matID)*constitutive_jump_width(i)*&
						    material_bg(matID)**3.0_pReal) 

   constitutive_g0_slip(i) = constitutive_rho_m(i)*material_bg(matID)*attack_frequency*constitutive_jump_width(i)*&
                             exp(-(material_Qedge(matID)+constitutive_passing_stress(i)*constitutive_activation_volume(i))/&
						     (kB*Tp))
enddo

!* Quantities derived from state - twin
Ftwin = sum(state((startIdxTwin+1):(startIdxTwin+material_Ntwin(matID))))
do i=1,material_Nslip(matID)
   !* Inverse of the average distance between 2 twins of the same familly
   constitutive_inv_intertwin_len_s(i)=0.0_pReal
   do j=1,material_Ntwin(matID)
      constitutive_inv_intertwin_len_s(i)=constitutive_inv_intertwin_len_s(i)+&
	                                     (lattice_SlipTwinIntType(i,j,material_CrystalStructure(matID))*state(startIdxTwin+j))/&
										 (2.0_pReal*material_StackSize(matID)*(1.0_pReal-Ftwin))
   enddo
enddo
do i=1,material_Ntwin(matID)
   !* Inverse of the average distance between 2 twins of the same familly
   constitutive_inv_intertwin_len_t(i)=0.0_pReal
   do j=1,material_Ntwin(matID)
      constitutive_inv_intertwin_len_t(i)=constitutive_inv_intertwin_len_t(i)+&
	                                     (lattice_TwinIntType(i,j,material_CrystalStructure(matID))*state(startIdxTwin+j))/&
										 (2.0_pReal*material_StackSize(matID)*(1.0_pReal-Ftwin))
   enddo
   constitutive_twin_mfp(i)=1.0_pReal/((1.0_pReal/material_GrainSize(matID))+constitutive_inv_intertwin_len_t(i))
   constitutive_twin_volume(i)=((4.0_pReal*pi)/3.0_pReal)*material_StackSize(matID)*constitutive_twin_mfp(i)**2.0_pReal
enddo

!* Stacking fault energy as function of temperature (see Allain PhD Thesis p40-42) *
x_Fe=0.774_pReal ! atomic %
x_Mn=0.218_pReal ! atomic %
x_C =0.027_pReal ! atomic %

!* Chemical contribution *
deltaG1=x_Fe*(4.309_pReal*Tp-2243.38_pReal)
deltaG2=x_Mn*(1.123_preal*Tp-1000.0_pReal)
deltaG3=x_Fe*x_Mn*(2873.0_pReal+717.0_pReal*(x_Fe-x_Mn)) 
deltaG4=1246.0_pReal*(1.0_pReal-exp(-24.29_pReal*x_C))-17175.0_pReal*x_Mn*x_C

!* Magnetical contribution *
beta_mart=0.62_pReal*x_Mn-4.0_pReal*x_C   ! Magnetic spin in µB
Tp_mart=580.0_pReal*x_Mn                  ! Néel Temperature
beta_aust=0.7_pReal*x_Fe+0.62_pReal*x_Mn-0.64_pReal*x_Fe*x_Mn-4.0_pReal*x_C ! Magnetic spin in µB
Tp_aust=669.27_pReal*(1.0_pReal-exp(-5.46_pReal*x_Mn))-2408.0_pReal*x_C*x_Fe-109.0_pReal ! Néel Temperature
if (Tp<=Tp_mart) then
   f_mart=1.0_pReal-(1.0_pReal/2.34_pReal)*((79.0_pReal*Tp_mart)/(140.0_pReal*0.28_pReal*Tp)+&
          (474.0_pReal/497.0_pReal)*((1.0_pReal/0.28_pReal)-1.0_pReal)*(((Tp/Tp_mart)**3.0_pReal)/6.0_pReal+&
          ((Tp/Tp_mart)**9.0_pReal)/135.0_pReal+((Tp/Tp_mart)**15.0_pReal)/600.0_pReal))
else
   f_mart=-(1.0_pReal/2.34_pReal)*(((Tp/Tp_mart)**-5.0_pReal)/10.0_pReal+((Tp/Tp_mart)**-15.0_pReal)/315.0_pReal+&
          ((Tp/Tp_mart)**-25.0_pReal)/1500.0_pReal)
endif
if (Tp<=Tp_aust) then
   f_aust=1.0_pReal-(1.0_pReal/2.34_pReal)*((79.0_pReal*Tp_aust)/(140.0_pReal*0.28_pReal*Tp)+&
          (474.0_pReal/497.0_pReal)*((1.0_pReal/0.28_pReal)-1.0_pReal)*(((Tp/Tp_aust)**3.0_pReal)/6.0_pReal+&
          ((Tp/Tp_aust)**9.0_pReal)/135.0_pReal+((Tp/Tp_aust)**15.0_pReal)/600.0_pReal))
else
   f_aust=-(1.0_pReal/2.34_pReal)*(((Tp/Tp_aust)**-5.0_pReal)/10.0_pReal+((Tp/Tp_aust)**-15.0_pReal)/315.0_pReal+&
          ((Tp/Tp_aust)**-25.0_pReal)/1500.0_pReal)
endif
deltaG5=Rgaz*Tp*(log(beta_mart+1.0_pReal)*f_mart-log(beta_aust+1.0_pReal)*f_aust)

!* Final expression *
material_sfe=(4.0_pReal/(sqrt(3.0_pReal)*avogadro*material_bg(matID)**2.0_pReal))*&
             (deltaG1+deltaG2+deltaG3+deltaG4+deltaG5)+2.0_pReal*0.0035_pReal

return	 
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
use lattice, only: lattice_Sslip,lattice_Sslip_v,lattice_Stwin,lattice_Stwin_v,lattice_TwinShear
use math, only: pi,math_Plain3333to99
implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
integer(pInt) matID,startIdxTwin,i,j,k,l,m,n
real(pReal) Tp,Ftwin
real(pReal), dimension(6) :: Tstar_v
real(pReal), dimension(3,3) :: Lp,Sslip,Stwin
real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
real(pReal), dimension(9,9) :: dLp_dTstar
real(pReal), dimension(constitutive_Nstatevars(ipc,ip,el)) :: state
real(pReal), dimension(material_Nslip(constitutive_matID(ipc,ip,el))) :: gdot_slip,dgdot_dtauslip,tau_slip
real(pReal), dimension(material_Ntwin(constitutive_matID(ipc,ip,el))) :: sfe_eff,fdot_twin,dfdot_dtautwin,tau_twin,tauc_twin
real(pReal), dimension(material_Ntwin(constitutive_matID(ipc,ip,el)),material_Nslip(constitutive_matID(ipc,ip,el))) :: dfdot_dtauslip

!* Get the material-ID from the triplet(ipc,ip,el)
matID = constitutive_matID(ipc,ip,el)
startIdxTwin = material_Nslip(matID)
Ftwin = sum(state((startIdxTwin+1):(startIdxTwin+material_Ntwin(matID))))

!* Calculation of Lp - slip
gdot_slip = 0.0_pReal
dgdot_dtauslip = 0.0_pReal	
Lp = 0.0_pReal	
do i=1,material_Nslip(matID)
   tau_slip(i)=dot_product(Tstar_v,lattice_Sslip_v(:,i,material_CrystalStructure(matID)))
   if (abs(tau_slip(i))>constitutive_passing_stress(i)) then
      gdot_slip(i) = constitutive_g0_slip(i)*sinh((tau_slip(i)*constitutive_activation_volume(i))/(kB*Tp))
      dgdot_dtauslip(i) = (constitutive_g0_slip(i)*constitutive_activation_volume(i))/(kB*Tp)*&
                          cosh((tau_slip(i)*constitutive_activation_volume(i))/(kB*Tp)) 
   endif       							  						       
   Lp=Lp+(1.0_pReal-Ftwin)*gdot_slip(i)*lattice_Sslip(:,:,i,material_CrystalStructure(matID))
enddo

!write(6,*) '##############'
!write(6,*) '##############'

!write(6,*) 'Schmid_1', lattice_Sslip_v(:,1,material_CrystalStructure(matID))
!write(6,*) 'Schmid_2', lattice_Sslip_v(:,2,material_CrystalStructure(matID))
!write(6,*) 'Schmid_3', lattice_Sslip_v(:,3,material_CrystalStructure(matID))
!write(6,*) 'Schmid_4', lattice_Sslip_v(:,4,material_CrystalStructure(matID))
!write(6,*) 'Schmid_5', lattice_Sslip_v(:,5,material_CrystalStructure(matID))
!write(6,*) 'Schmid_6', lattice_Sslip_v(:,6,material_CrystalStructure(matID))
!write(6,*) 'Schmid_7', lattice_Sslip_v(:,7,material_CrystalStructure(matID))
!write(6,*) 'Schmid_8', lattice_Sslip_v(:,8,material_CrystalStructure(matID))
!write(6,*) 'Schmid_9', lattice_Sslip_v(:,9,material_CrystalStructure(matID))
!write(6,*) 'Schmid_10', lattice_Sslip_v(:,10,material_CrystalStructure(matID))
!write(6,*) 'Schmid_11', lattice_Sslip_v(:,11,material_CrystalStructure(matID))
!write(6,*) 'Schmid_12', lattice_Sslip_v(:,12,material_CrystalStructure(matID))
!write(6,*) 'Tstar_v',Tstar_v
!write(6,*) 'state',state
!write(6,*) 'Tp', Tp
!write(6,*) 'ssd_f', constitutive_rho_f
!write(6,*) 'ssd_p', constitutive_rho_p
!write(6,*) 'jump_width', constitutive_jump_width
!write(6,*) 'activation_volume', constitutive_activation_volume 
!write(6,*) 'passing_stress', constitutive_passing_stress
!write(6,*) 'ssd_m', constitutive_rho_m
!write(6,*) 'g0_slip', constitutive_g0_slip
!write(6,*) 'tau_slip',tau_slip
!write(6,*) 'gdot_slip', gdot_slip
!write(6,*) 'dgdot_dtauslip', dgdot_dtauslip


!* Calculation of Lp - twin
sfe_eff = 0.0_pReal
fdot_twin = 0.0_pReal
dfdot_dtautwin = 0.0_pReal
dfdot_dtauslip = 0.0_pReal
do i=1,material_Ntwin(matID)
   tau_twin(i)=dot_product(Tstar_v,lattice_Stwin_v(:,i,material_CrystalStructure(matID)))
   if ((tau_twin(i) > 0.0_pReal).AND.(material_TwinSaturation(matID)-Ftwin>=0)) then
      sfe_eff(i)=material_sfe-(sqrt(3.0_pReal)/3.0_pReal)*material_q1(matID)*material_q2(matID)*material_bg(matID)*tau_twin(i)
      if (sfe_eff(i)<0.0_pReal) sfe_eff(i) = 0.0_pReal
	  fdot_twin(i) = (material_TwinSaturation(matID)-Ftwin)*&
	                 constitutive_twin_volume(i)*&
	                 ((2.0_pReal*sqrt(6.0_pReal)*material_SiteScaling(matID)*sum(abs(gdot_slip))*&
					 sum(state(1:material_Nslip(matID)))**1.5_pReal)/3.0_pReal)*&
	                 exp((-25.0_pReal*pi**3.0_pReal*material_Gmod(matID)**2.0_pReal*sfe_eff(i)**3.0_pReal)/&
	                     (3.0_pReal*Kb*Tp*(material_q2(matID)*tau_twin(i))**4.0_pReal))
					 				 								  
	  dfdot_dtautwin(i) = (material_TwinSaturation(matID)-Ftwin)*&
	                      constitutive_twin_volume(i)*&
	                      ((2.0_pReal*sqrt(6.0_pReal)*material_SiteScaling(matID)*sum(abs(gdot_slip))*&
					      sum(state(1:material_Nslip(matID)))**1.5_pReal)/3.0_pReal)*&	  
                          ((-25.0_pReal*pi**3.0_pReal*material_Gmod(matID)**2.0_pReal*sfe_eff(i)**2.0_pReal)/&
				          (3.0_pReal*Kb*Tp*material_q2(matID)**4.0_pReal*tau_twin(i)**5.0_pReal))*&
                          (-sqrt(3.0_pReal)*material_q1(matID)*material_q2(matID)*material_bg(matID)*tau_twin(i)-&
					      4.0_pReal*sfe_eff(i))*&
	                      exp((-25.0_pReal*pi**3.0_pReal*material_Gmod(matID)**2.0_pReal*sfe_eff(i)**3.0_pReal)/&
	                          (3.0_pReal*Kb*Tp*(material_q2(matID)*tau_twin(i))**4.0_pReal))
	  
	  do j=1,material_Nslip(matID)
	     if (gdot_slip(j)>0.0_pReal) then
	        dfdot_dtauslip(i,j) = (material_TwinSaturation(matID)-Ftwin)*&
	                              constitutive_twin_volume(i)*&
	                              ((2.0_pReal*sqrt(6.0_pReal)*material_SiteScaling(matID)*dgdot_dtauslip(j)*&
					              sum(state(1:material_Nslip(matID)))**1.5_pReal)/3.0_pReal)*&
	                              exp((-25.0_pReal*pi**3.0_pReal*material_Gmod(matID)**2.0_pReal*sfe_eff(i)**3.0_pReal)/&
	                              (3.0_pReal*Kb*Tp*(material_q2(matID)*tau_twin(i))**4.0_pReal))	                                    
	  	 else
	  	    dfdot_dtauslip(i,j) = (material_TwinSaturation(matID)-Ftwin)*&
	                              constitutive_twin_volume(i)*&
	                              ((2.0_pReal*sqrt(6.0_pReal)*material_SiteScaling(matID)*(-dgdot_dtauslip(j))*&
					              sum(state(1:material_Nslip(matID)))**1.5_pReal)/3.0_pReal)*&
	                              exp((-25.0_pReal*pi**3.0_pReal*material_Gmod(matID)**2.0_pReal*sfe_eff(i)**3.0_pReal)/&
	                              (3.0_pReal*Kb*Tp*(material_q2(matID)*tau_twin(i))**4.0_pReal))
	  	 endif
	  enddo
   endif
   Lp=Lp+state(material_Nslip(matID)+i)*lattice_TwinShear(material_CrystalStructure(matID))*constitutive_fdot_twin(i)*&
      lattice_Stwin(:,:,i,material_CrystalStructure(matID))
enddo


!write(6,*) 'twin_mfp', constitutive_twin_mfp
!write(6,*) 'twin_volume', constitutive_twin_volume
!write(6,*) 'tau_twin',tau_twin
!write(6,*) 'part1:',material_TwinSaturation(matID)-Ftwin
!write(6,*) 'part2:',constitutive_twin_volume
!write(6,*) 'part3:',((2.0_pReal*sqrt(6.0_pReal)*sum(abs(gdot_slip))*&
!					 sum(state(1:material_Nslip(matID)))**1.5_pReal)/3.0_pReal)
!do i=1,12
!write(6,*) 'part4:',exp((-25.0_pReal*pi**3.0_pReal*material_Gmod(matID)**2.0_pReal*sfe_eff(i)**3.0_pReal)/&
!	                     (3.0_pReal*Kb*Tp*(material_q2(matID)*tau_twin(i))**4.0_pReal))
!write(6,*) 'part5:',(-25.0_pReal*pi**3.0_pReal*material_Gmod(matID)**2.0_pReal*sfe_eff(i)**3.0_pReal)/&
!	                     (3.0_pReal*Kb*Tp*(material_q2(matID)*tau_twin(i))**4.0_pReal)
!enddo
!write(6,*) 'sfe', material_sfe
!write(6,*) 'sfe_eff',sfe_eff
!write(6,*) 'fdot_twin', fdot_twin
!write(6,*) 'dfdot_dtautwin', dfdot_dtautwin
!write(6,*) 'dfdot_dtauslip', dfdot_dtauslip

	            
!* Calculation of the tangent of Lp
dLp_dTstar3333=0.0_pReal
do i=1,material_Nslip(matID)
   Sslip = lattice_Sslip(:,:,i,material_CrystalStructure(matID))
   forall (k=1:3,l=1:3,m=1:3,n=1:3)
     dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n)+ &
                               (1.0_pReal-Ftwin)*dgdot_dtauslip(i)*Sslip(k,l)*Sslip(m,n) !force m,n symmetry
   endforall
enddo
do i=1,material_Ntwin(matID)
   Stwin = lattice_Stwin(:,:,i,material_CrystalStructure(matID)) 
   forall (k=1:3,l=1:3,m=1:3,n=1:3)
     dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n)+ &
                               state(material_Nslip(matID)+i)*lattice_TwinShear(material_CrystalStructure(matID))*&
                               dfdot_dtautwin(i)*Stwin(k,l)*Stwin(m,n) !force m,n symmetry
   endforall
   do j=1,material_Nslip(matID)
      Sslip = lattice_Sslip(:,:,j,material_CrystalStructure(matID))
      forall (k=1:3,l=1:3,m=1:3,n=1:3)
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n)+ &
                                   state(material_Nslip(matID)+i)*lattice_TwinShear(material_CrystalStructure(matID))*&
   								   dfdot_dtauslip(i,j)*Stwin(k,l)*Sslip(m,n) !force m,n symmetry
      endforall  
   enddo
enddo
dLp_dTstar = math_Plain3333to99(dLp_dTstar3333)

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
use math, only: pi
use lattice, only: lattice_Sslip_v,lattice_Stwin_v
implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
integer(pInt) matID,i,j,startIdxTwin
real(pReal) Tp,Ftwin
real(pReal), dimension(6) :: Tstar_v
real(pReal), dimension(constitutive_Nstatevars(ipc,ip,el)) :: constitutive_dotState,state
real(pReal), dimension(material_Nslip(constitutive_matID(ipc,ip,el))) :: gdot_slip,tau_slip
real(pReal), dimension(material_Ntwin(constitutive_matID(ipc,ip,el))) :: sfe_eff,fdot_twin,tau_twin,tauc_twin
real(pReal), dimension(material_Nslip(constitutive_matID(ipc,ip,el))) :: locks,grainboundaries,recovery
real(pReal), dimension(material_Nslip(constitutive_matID(ipc,ip,el))) :: twinboundaries

!* Get the material-ID from the triplet(ipc,ip,el)
matID = constitutive_matID(ipc,ip,el)
startIdxTwin = material_Nslip(matID)
Ftwin = sum(state((startIdxTwin+1):(startIdxTwin+material_Ntwin(matID))))
constitutive_dotState = 0.0_pReal

!* Dislocation density evolution
gdot_slip = 0.0_pReal
do i=1,material_Nslip(matID)
   tau_slip(i)=dot_product(Tstar_v,lattice_Sslip_v(:,i,material_CrystalStructure(matID)))
   if (abs(tau_slip(i))>constitutive_passing_stress(i)) &
      gdot_slip(i) = constitutive_g0_slip(i)*sinh((tau_slip(i)*constitutive_activation_volume(i))/(kB*Tp)) 
	    
   locks(i)           = (sqrt(constitutive_rho_f(i))*abs(gdot_slip(i)))/(material_c4(matID)*material_bg(matID)) 
   grainboundaries(i) = abs(gdot_slip(i))/(material_c5(matID)*material_bg(matID)*material_GrainSize(matID))
   twinboundaries(i)  = (abs(gdot_slip(i))*constitutive_inv_intertwin_len_s(i))/(material_c6(matID)*material_bg(matID))
   recovery(i)        = material_c7(matID)*state(i)*abs(gdot_slip(i))

   constitutive_dotState(i) = locks(i)+grainboundaries(i)+twinboundaries(i)-recovery(i)
enddo

!* Twin volume fraction evolution
fdot_twin = 0.0_pReal
do i=1,material_Ntwin(matID)
   tau_twin(i)=dot_product(Tstar_v,lattice_Stwin_v(:,i,material_CrystalStructure(matID)))
   if (tau_twin(i)>0.0_pReal) then
      sfe_eff(i)=material_sfe-(sqrt(3.0_pReal)/3.0_pReal)*material_q1(matID)*material_q2(matID)*material_bg(matID)*tau_twin(i)
      if (sfe_eff(i)<0.0_pReal) sfe_eff(i) = 0.0_pReal
	  fdot_twin(i) = (material_TwinSaturation(matID)-Ftwin)*&
	                 constitutive_twin_volume(i)*&
	                 ((2.0_pReal*sqrt(6.0_pReal)*material_SiteScaling(matID)*sum(abs(gdot_slip))*&
					 sum(state(1:material_Nslip(matID)))**1.5_pReal)/3.0_pReal)*&
	                 exp((-25.0_pReal*pi**3.0_pReal*material_Gmod(matID)**2.0_pReal*sfe_eff(i)**3.0_pReal)/&
	                     (3.0_pReal*Kb*Tp*(material_q2(matID)*tau_twin(i))**4.0_pReal))
   endif			 								  
   constitutive_dotState(startIdxTwin+i) = fdot_twin(i)
enddo

return
end function


function constitutive_post_results(Tstar_v,state,Tp,dt,ipc,ip,el)
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
!* constitutive_Nresults has to be set accordingly in _Assignment    *
!*********************************************************************
use prec, only: pReal,pInt
use lattice, only: lattice_Sslip_v
implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
integer(pInt) matID,i,startIdxTwin
real(pReal) dt,Tp,tau_slip
real(pReal), dimension(6) :: Tstar_v
real(pReal), dimension(constitutive_Nstatevars(ipc,ip,el)) :: state
real(pReal), dimension(constitutive_Nresults(ipc,ip,el))   :: constitutive_post_results

!* Get the material-ID from the triplet(ipc,ip,el)
matID = constitutive_matID(ipc,ip,el)
startIdxTwin = material_Nslip(matID)

if(constitutive_Nresults(ipc,ip,el)==0) return

constitutive_post_results=0.0_pReal
!do i=1,material_Nslip(matID)
!   constitutive_post_results(i) = state(i)
!enddo
!do i=1,material_Ntwin(matID)
!   constitutive_post_results(startIdxTwin+i) = state(startIdxTwin+i)
!enddo
constitutive_post_results(1) = sum(state(1:material_Nslip(matID)))
constitutive_post_results(2) = sum(state((startIdxTwin+1):(startIdxTwin+material_Ntwin(matID))))

return
end function

END MODULE