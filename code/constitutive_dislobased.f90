
!************************************
!*      Module: CONSTITUTIVE        *
!************************************
!* contains:                        *
!* - constitutive equations         *
!* - parameters definition          *
!* - orientations                   *
!************************************

!	[Alu]
!	constitution             dislobased
!	(output)                 dislodensity
!	(output)                 rateofshear
!	lattice_structure        1
!	Nslip                    12
!	
!	c11                      106.75e9
!	c12                      60.41e9
!	c44                      28.34e9
!	
!	burgers                  2.86e-10	# Burgers vector [m]
!	Qedge                    3e-19		# Activation energy for dislocation glide [J/K] (0.5*G*b^3)
!	Qsd                      2.4e-19		# Activation energy for self diffusion [J/K] (gamma-iron)
!	diff0                    1e-3		# prefactor vacancy diffusion coeffficent (gamma-iron)
!	interaction_coefficients 1.0 2.2 3.0 1.6 3.8 4.5		# Dislocation interaction coefficients
!	
!	rho0                     6.0e12		# Initial dislocation density [m/m^3]
!	
!	c1                       0.1		# Passing stress adjustment
!	c2                       2.0		# Jump width adjustment
!	c3                       1.0		# Activation volume adjustment
!	c4                       50.0		# Average slip distance adjustment for lock formation
!	c7                       8.0		# Athermal recovery adjustment
!	c8                       1.0e10		# Thermal recovery adjustment (plays no role for me)

MODULE constitutive_dislobased
!*** Include other modules ***
 use prec, only: pReal,pInt
 implicit none

 character (len=*), parameter :: constitutive_dislobased_label = 'dislobased'
 
 integer(pInt),   dimension(:),     allocatable :: constitutive_dislobased_sizeDotState, &
                                                   constitutive_dislobased_sizeState, &
                                                   constitutive_dislobased_sizePostResults
 integer(pInt),   dimension(:,:),   allocatable,target :: constitutive_dislobased_sizePostResult     ! size of each post result output
 character(len=64), dimension(:,:), allocatable,target :: constitutive_dislobased_output             ! name of each post result output
 character(len=32), dimension(:),   allocatable :: constitutive_dislobased_structureName
 integer(pInt),   dimension(:),     allocatable :: constitutive_dislobased_structure
 integer(pInt),   dimension(:),     allocatable :: constitutive_dislobased_Nslip
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_C11
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_C12
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_C13
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_C33
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_C44
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_Gmod
 real(pReal), dimension(:,:,:), allocatable :: constitutive_dislobased_Cslip_66
!* Visco-plastic constitutive_phenomenological parameters
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_rho0
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_bg
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_Qedge
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_Qsd
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_D0
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_c1
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_c2
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_c3
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_c4
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_c5
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_c6
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_c7
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_c8
 real(pReal), dimension(:),     allocatable :: constitutive_dislobased_CoverA
 real(pReal), dimension(:,:),   allocatable :: constitutive_dislobased_SlipIntCoeff
 real(pReal), dimension(:,:,:), allocatable :: constitutive_dislobased_Iparallel
 real(pReal), dimension(:,:,:), allocatable :: constitutive_dislobased_Iforest

!*************************************
!* Definition of material properties *
!*************************************
!* Physical parameter, attack_frequency != Debye frequency
real(pReal), parameter :: attack_frequency = 1.0e10_pReal  
!* Physical parameter, Boltzmann constant in J/Kelvin
real(pReal), parameter :: kB = 1.38e-23_pReal
!* Physical parameter, Avogadro number in 1/mol
real(pReal), parameter :: avogadro = 6.022e23_pReal 
!* Physical parameter, Gas constant in J.mol/Kelvin
real(pReal), parameter :: Rgaz = 8.314_pReal    

CONTAINS
!****************************************
!* - constitutive_init
!* - constitutive_homogenizedC
!* - constitutive_microstructure
!* - constitutive_LpAndItsTangent
!* - consistutive_dotState
!* - consistutive_postResults
!****************************************

subroutine constitutive_dislobased_init(file)
!**************************************
!*      Module initialization         *
!**************************************
 use prec, only: pInt, pReal
 use math, only: math_Mandel3333to66, math_Voigt66to3333, math_mul3x3
 use IO
 use material
 use lattice, only: lattice_sn, lattice_st, lattice_interactionSlipSlip, lattice_initializeStructure
 integer(pInt), intent(in) :: file
 integer(pInt), parameter :: maxNchunks = 7
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) section, maxNinstance, i,j,k,l, output, mySize
 character(len=64) tag
 character(len=1024) line
 real(pReal) x,y
 
 write(6,*)
 write(6,'(a20,a20,a12)') '<<<+-  constitutive_',constitutive_dislobased_label,' init  -+>>>'
 write(6,*)
 
 maxNinstance = count(phase_constitution == constitutive_dislobased_label)
 if (maxNinstance == 0) return

 allocate(constitutive_dislobased_sizeDotState(maxNinstance)) ;   constitutive_dislobased_sizeDotState = 0_pInt
 allocate(constitutive_dislobased_sizeState(maxNinstance)) ;      constitutive_dislobased_sizeState = 0_pInt
 allocate(constitutive_dislobased_sizePostResults(maxNinstance)); constitutive_dislobased_sizePostResults = 0_pInt
 allocate(constitutive_dislobased_sizePostResult(maxval(phase_Noutput), &
                                                 maxNinstance)) ; constitutive_dislobased_sizePostResult = 0_pInt
 allocate(constitutive_dislobased_output(maxval(phase_Noutput), &
                                         maxNinstance)) ;         constitutive_dislobased_output = ''
 allocate(constitutive_dislobased_structureName(maxNinstance)) ;  constitutive_dislobased_structureName = ''
 allocate(constitutive_dislobased_structure(maxNinstance)) ;      constitutive_dislobased_structure = 0_pInt
 allocate(constitutive_dislobased_Nslip(maxNinstance)) ;          constitutive_dislobased_Nslip = 0_pInt
 allocate(constitutive_dislobased_C11(maxNinstance)) ;            constitutive_dislobased_C11 = 0.0_pReal
 allocate(constitutive_dislobased_C12(maxNinstance)) ;            constitutive_dislobased_C12 = 0.0_pReal
 allocate(constitutive_dislobased_C13(maxNinstance)) ;            constitutive_dislobased_C13 = 0.0_pReal
 allocate(constitutive_dislobased_C33(maxNinstance)) ;            constitutive_dislobased_C33 = 0.0_pReal
 allocate(constitutive_dislobased_C44(maxNinstance)) ;            constitutive_dislobased_C44 = 0.0_pReal
 allocate(constitutive_dislobased_Gmod(maxNinstance)) ;           constitutive_dislobased_Gmod = 0.0_pReal
 allocate(constitutive_dislobased_Cslip_66(6,6,maxNinstance)) ;   constitutive_dislobased_Cslip_66 = 0.0_pReal
 allocate(constitutive_dislobased_rho0(maxNinstance))         ;   constitutive_dislobased_rho0 = 0.0_pReal
 allocate(constitutive_dislobased_bg(maxNinstance))           ;   constitutive_dislobased_bg = 0.0_pReal
 allocate(constitutive_dislobased_Qedge(maxNinstance))        ;   constitutive_dislobased_Qedge = 0.0_pReal
 allocate(constitutive_dislobased_Qsd(maxNinstance))          ;   constitutive_dislobased_Qsd = 0.0_pReal
 allocate(constitutive_dislobased_D0(maxNinstance))           ;   constitutive_dislobased_D0 = 0.0_pReal
 allocate(constitutive_dislobased_c1(maxNinstance))           ;   constitutive_dislobased_c1 = 0.0_pReal
 allocate(constitutive_dislobased_c2(maxNinstance))           ;   constitutive_dislobased_c2 = 0.0_pReal
 allocate(constitutive_dislobased_c3(maxNinstance))           ;   constitutive_dislobased_c3 = 0.0_pReal
 allocate(constitutive_dislobased_c4(maxNinstance))           ;   constitutive_dislobased_c4 = 0.0_pReal
 allocate(constitutive_dislobased_c5(maxNinstance))           ;   constitutive_dislobased_c5 = 0.0_pReal
 allocate(constitutive_dislobased_c6(maxNinstance))           ;   constitutive_dislobased_c6 = 0.0_pReal
 allocate(constitutive_dislobased_c7(maxNinstance))           ;   constitutive_dislobased_c7 = 0.0_pReal
 allocate(constitutive_dislobased_c8(maxNinstance))           ;   constitutive_dislobased_c8 = 0.0_pReal
 allocate(constitutive_dislobased_CoverA(maxNinstance))       ;   constitutive_dislobased_CoverA = 0.0_pReal
 allocate(constitutive_dislobased_SlipIntCoeff(6,maxNinstance)) ; constitutive_dislobased_SlipIntCoeff = 0.0_pReal

 rewind(file)
 line = ''
 section = 0
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= 'phase')     ! wind forward to <phase>
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
   if (section > 0 .and. phase_constitution(section) == constitutive_dislobased_label) then  ! one of my sections
     i = phase_constitutionInstance(section)     ! which instance of my constitution is present phase
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     select case(tag)
       case ('(output)')
         output = output + 1
         constitutive_dislobased_output(output,i) = IO_lc(IO_stringValue(line,positions,2))
       case ('lattice_structure')
              constitutive_dislobased_structureName(i) = IO_lc(IO_stringValue(line,positions,2))
       case ('covera_ratio')
              constitutive_dislobased_CoverA(i) = IO_floatValue(line,positions,2)
       case ('nslip')
              constitutive_dislobased_Nslip(i) = IO_intValue(line,positions,2)
       case ('c11')
              constitutive_dislobased_C11(i) = IO_floatValue(line,positions,2)
       case ('c12')
              constitutive_dislobased_C12(i) = IO_floatValue(line,positions,2)
       case ('c13')
              constitutive_dislobased_C13(i) = IO_floatValue(line,positions,2)
       case ('c33')
              constitutive_dislobased_C33(i) = IO_floatValue(line,positions,2)
       case ('c44')
              constitutive_dislobased_C44(i) = IO_floatValue(line,positions,2)
       case ('rho0')
              constitutive_dislobased_rho0(i) = IO_floatValue(line,positions,2)
       case ('burgers')
              constitutive_dislobased_bg(i) = IO_floatValue(line,positions,2)
       case ('qedge')
              constitutive_dislobased_Qedge(i) = IO_floatValue(line,positions,2)
       case ('qsd')
              constitutive_dislobased_Qsd(i) = IO_floatValue(line,positions,2)
       case ('diff0')
              constitutive_dislobased_D0(i) = IO_floatValue(line,positions,2)
       case ('c1')
              constitutive_dislobased_c1(i) = IO_floatValue(line,positions,2)
       case ('c2') 
              constitutive_dislobased_c2(i) = IO_floatValue(line,positions,2)
       case ('c3') 
              constitutive_dislobased_c3(i) = IO_floatValue(line,positions,2)
       case ('c4') 
              constitutive_dislobased_c4(i) = IO_floatValue(line,positions,2)
       case ('c5') 
              constitutive_dislobased_c5(i) = IO_floatValue(line,positions,2)
       case ('c6') 
              constitutive_dislobased_c6(i) = IO_floatValue(line,positions,2)
       case ('c7') 
              constitutive_dislobased_c7(i) = IO_floatValue(line,positions,2)
       case ('c8') 
              constitutive_dislobased_c8(i) = IO_floatValue(line,positions,2)
       case ('interaction_coefficients') 
         forall (j=1:6) &
           constitutive_dislobased_SlipIntCoeff(j,i) = IO_floatValue(line,positions,1+j)
     end select
   endif
 enddo

  
100 do i = 1,maxNinstance
   constitutive_dislobased_structure(i) = lattice_initializeStructure(constitutive_dislobased_structureName(i), &
                                                                      constitutive_dislobased_CoverA(i))
! sanity checks
   if (constitutive_dislobased_structure(i) < 1)           call IO_error(205)
   if (constitutive_dislobased_rho0(i) < 0.0_pReal)        call IO_error(220)
   if (constitutive_dislobased_bg(i) <= 0.0_pReal)         call IO_error(221)
   if (constitutive_dislobased_Qedge(i) <= 0.0_pReal)      call IO_error(222)
   if (constitutive_dislobased_Qsd(i) <= 0.0_pReal)        call IO_error(223)
   if (constitutive_dislobased_D0(i) <= 0.0_pReal)         call IO_error(224)
   if (constitutive_dislobased_Nslip(i) < 1)               call IO_error(225)
 enddo

 allocate(constitutive_dislobased_Iparallel(maxval(constitutive_dislobased_Nslip),&
                                            maxval(constitutive_dislobased_Nslip),&
                                            maxNinstance))

 allocate(constitutive_dislobased_Iforest(maxval(constitutive_dislobased_Nslip),&
                                          maxval(constitutive_dislobased_Nslip),&
                                          maxNinstance))

 do i = 1,maxNinstance
   do j = 1,maxval(phase_Noutput)
	 select case(constitutive_dislobased_output(j,i))
	   case('dislodensity')
		 mySize = constitutive_dislobased_Nslip(i)
	   case('rateofshear')
		 mySize = constitutive_dislobased_Nslip(i)
	   case default
		 mySize = 0_pInt
	 end select

	 if (mySize > 0_pInt) then                               ! any meaningful output found
	   constitutive_dislobased_sizePostResult(j,i) = mySize
	   constitutive_dislobased_sizePostResults(i) = &
	   constitutive_dislobased_sizePostResults(i) + mySize
	 endif
   enddo
   
   constitutive_dislobased_sizeDotState(i) = constitutive_dislobased_Nslip(i)
   constitutive_dislobased_sizeState(i)    = 8*constitutive_dislobased_Nslip(i)

   constitutive_dislobased_Gmod(i) = constitutive_dislobased_C44(i)
   select case (constitutive_dislobased_structure(i))
   case(1:2) ! cubic(s)
     forall(k=1:3)
       forall(j=1:3) &
         constitutive_dislobased_Cslip_66(k,j,i)     = constitutive_dislobased_C12(i)
         constitutive_dislobased_Cslip_66(k,k,i)     = constitutive_dislobased_C11(i)
         constitutive_dislobased_Cslip_66(k+3,k+3,i) = constitutive_dislobased_C44(i)
     end forall
   case(3:)   ! all hex
     constitutive_dislobased_Cslip_66(1,1,i) = constitutive_dislobased_C11(i)
     constitutive_dislobased_Cslip_66(2,2,i) = constitutive_dislobased_C11(i)
     constitutive_dislobased_Cslip_66(3,3,i) = constitutive_dislobased_C33(i)
     constitutive_dislobased_Cslip_66(1,2,i) = constitutive_dislobased_C12(i)
     constitutive_dislobased_Cslip_66(2,1,i) = constitutive_dislobased_C12(i)
     constitutive_dislobased_Cslip_66(1,3,i) = constitutive_dislobased_C13(i)
     constitutive_dislobased_Cslip_66(3,1,i) = constitutive_dislobased_C13(i)
     constitutive_dislobased_Cslip_66(2,3,i) = constitutive_dislobased_C13(i)
     constitutive_dislobased_Cslip_66(3,2,i) = constitutive_dislobased_C13(i)
     constitutive_dislobased_Cslip_66(4,4,i) = constitutive_dislobased_C44(i)
     constitutive_dislobased_Cslip_66(5,5,i) = constitutive_dislobased_C44(i)
     constitutive_dislobased_Cslip_66(6,6,i) = 0.5_pReal*(constitutive_dislobased_C11(i)- &
                                                          constitutive_dislobased_C12(i))
   end select
   constitutive_dislobased_Cslip_66(:,:,i) = &
     math_Mandel3333to66(math_Voigt66to3333(constitutive_dislobased_Cslip_66(:,:,i)))


   !* Construction of the hardening matrices
   !* Iteration over the systems
   do j = 1,constitutive_dislobased_Nslip(i)
      do k = 1,constitutive_dislobased_Nslip(i)
       !* Projection of the dislocation *
        x = math_mul3x3(lattice_sn(:,j,i),lattice_st(:,k,i))
        y = 1.0_pReal-x**(2.0_pReal)
       !* Interaction matrix *
        constitutive_dislobased_Iforest(j,k,i) = abs(x)*&
        constitutive_dislobased_SlipIntCoeff(lattice_interactionSlipSlip(j,k,constitutive_dislobased_structure(i)),i)
        if (y>0.0_pReal) &
          constitutive_dislobased_Iparallel(j,k,i) = sqrt(y)*&
          constitutive_dislobased_SlipIntCoeff(lattice_interactionSlipSlip(j,k,constitutive_dislobased_structure(i)),i)
      enddo
   enddo

 enddo

 return

end subroutine


function constitutive_dislobased_stateInit(myInstance)
!*********************************************************************
!* initial microstructural state                                     *
!*********************************************************************
 use prec, only: pReal,pInt
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: myInstance
 real(pReal), dimension(constitutive_dislobased_Nslip(myInstance)) :: constitutive_dislobased_stateInit

 constitutive_dislobased_stateInit = constitutive_dislobased_rho0(myInstance)

 return
end function

function constitutive_dislobased_homogenizedC(state,ipc,ip,el)
!*********************************************************************
!* homogenized elacticity matrix                                     *
!* INPUT:                                                            *
!*  - state           : state variables                              *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_constitutionInstance
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: ipc,ip,el
 integer(pInt) matID
 real(pReal), dimension(6,6) :: constitutive_dislobased_homogenizedC
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 
 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 constitutive_dislobased_homogenizedC = constitutive_dislobased_Cslip_66(:,:,matID)

 return

end function


subroutine constitutive_dislobased_microstructure(Temperature,state,ipc,ip,el)
!*********************************************************************
!* calculate derived quantities from state (not used here)           *
!* INPUT:                                                            *
!*  - Tp              : temperature                                  *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_constitutionInstance
 implicit none

!* Definition of variables
 integer(pInt) ipc,ip,el,matID,n,i
 real(pReal) Temperature
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state

 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 n = constitutive_dislobased_Nslip(matID)
 !* Quantities derived from state - slip
 !* State: 1    : n   rho
 !*        n+1  : 2n  rho_f
 !*        2n+1 : 3n  rho_p
 !*        3n+1 : 4n  passing stress
 !*        4n+1 : 5n  jump width
 !*        5n+1 : 6n  activation volume
 !*        6n+1 : 7n  rho_m
 !*        7n+1 : 8n  g0_slip
 !$OMP CRITICAL (evilmatmul)
 state(ipc,ip,el)%p((n+1):(2*n))   = matmul(constitutive_dislobased_Iforest  (1:n,1:n,matID),state(ipc,ip,el)%p(1:n))
 state(ipc,ip,el)%p((2*n+1):(3*n)) = matmul(constitutive_dislobased_Iparallel(1:n,1:n,matID),state(ipc,ip,el)%p(1:n))
 !$OMP END CRITICAL (evilmatmul)

 do i=1,n

   state(ipc,ip,el)%p(3*n+i) = &
   constitutive_dislobased_c1(matID)*constitutive_dislobased_Gmod(matID)*&
   constitutive_dislobased_bg(matID)*sqrt(state(ipc,ip,el)%p(2*n+i)) 
   
   state(ipc,ip,el)%p(4*n+i) = &
   constitutive_dislobased_c2(matID)/sqrt(state(ipc,ip,el)%p(n+i))

   state(ipc,ip,el)%p(5*n+i) = &
   constitutive_dislobased_c3(matID)*state(ipc,ip,el)%p(4*n+i)*constitutive_dislobased_bg(matID)**2.0_pReal

   state(ipc,ip,el)%p(6*n+i) = &
   (2.0_pReal*kB*Temperature*sqrt(state(ipc,ip,el)%p(2*n+i)))/&
   (constitutive_dislobased_c1(matID)*constitutive_dislobased_c3(matID)*constitutive_dislobased_Gmod(matID)*&
   state(ipc,ip,el)%p(4*n+i)*constitutive_dislobased_bg(matID)**3.0_pReal)
    
   state(ipc,ip,el)%p(7*n+i) = &
   state(ipc,ip,el)%p(6*n+i)*constitutive_dislobased_bg(matID)*attack_frequency*state(ipc,ip,el)%p(4*n+i)*&
   exp(-constitutive_dislobased_Qedge(matID)/(kB*Temperature))

 enddo

end subroutine


subroutine constitutive_dislobased_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,state,ipc,ip,el)
!*********************************************************************
!* plastic velocity gradient and its tangent                         *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - Lp              : plastic velocity gradient                    *
!*  - dLp_dTstar      : derivative of Lp (4th-rank tensor)           *
!*********************************************************************
 use prec, only: pReal,pInt,p_vec
 use math, only: math_Plain3333to99, math_mul6x6
 use lattice, only: lattice_Sslip,lattice_Sslip_v
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_constitutionInstance

 implicit none

!* Definition of variables
 integer(pInt) ipc,ip,el
 integer(pInt) matID,i,k,l,m,n
 real(pReal) Temperature
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(3,3) :: Lp
 real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
 real(pReal), dimension(9,9) :: dLp_dTstar
 real(pReal), dimension(constitutive_dislobased_Nslip(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip,dgdot_dtauslip,tau_slip

 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 n = constitutive_dislobased_Nslip(matID)

!* Calculation of Lp
 Lp = 0.0_pReal
 gdot_slip = 0.0_pReal
 do i = 1,constitutive_dislobased_Nslip(matID)
   tau_slip(i)  = math_mul6x6(Tstar_v,lattice_Sslip_v(:,i,constitutive_dislobased_structure(matID)))
   if ((abs(tau_slip(i))-state(ipc,ip,el)%p(3*n+i))>0) &
      gdot_slip(i) = state(ipc,ip,el)%p(7*n+i)*sign(1.0_pReal,tau_slip(i))*&
                     sinh(((abs(tau_slip(i))-state(ipc,ip,el)%p(3*n+i))*state(ipc,ip,el)%p(5*n+i))/(kB*Temperature)) 
				   
   Lp = Lp + gdot_slip(i)*lattice_Sslip(:,:,i,constitutive_dislobased_structure(matID))
 enddo

!* Calculation of the tangent of Lp
 dLp_dTstar3333 = 0.0_pReal
 dLp_dTstar = 0.0_pReal
 dgdot_dtauslip = 0.0_pReal
 do i = 1,constitutive_dislobased_Nslip(matID)
   if ((abs(tau_slip(i))-state(ipc,ip,el)%p(3*n+i))>0) &
      dgdot_dtauslip(i) = (state(ipc,ip,el)%p(7*n+i)*state(ipc,ip,el)%p(5*n+i))/(kB*Temperature)*&
                          cosh(((abs(tau_slip(i))-state(ipc,ip,el)%p(3*n+i))*state(ipc,ip,el)%p(5*n+i))/(kB*Temperature)) 
   forall (k=1:3,l=1:3,m=1:3,n=1:3) &
          dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
             dgdot_dtauslip(i)*lattice_Sslip(k,l,i,constitutive_dislobased_structure(matID))* &
                               lattice_Sslip(m,n,i,constitutive_dislobased_structure(matID))
 enddo
 dLp_dTstar = math_Plain3333to99(dLp_dTstar3333)

 return
end subroutine


function constitutive_dislobased_dotState(Tstar_v,Temperature,state,ipc,ip,el)
!*********************************************************************
!* rate of change of microstructure                                  *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - constitutive_dotState : evolution of state variable            *
!*********************************************************************
 use prec, only: pReal,pInt,p_vec
 use lattice, only: lattice_Sslip_v
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_constitutionInstance
 implicit none

!* Definition of variables
 integer(pInt) ipc,ip,el
 integer(pInt) matID,i,n
 real(pReal) Temperature,tau_slip,gdot_slip,locks,athermal_recovery,thermal_recovery
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(constitutive_dislobased_Nslip(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_dislobased_dotState

 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 n = constitutive_dislobased_Nslip(matID)

!* Dislocation density evolution
 constitutive_dislobased_dotState = 0.0_pReal
 do i = 1,n
   tau_slip = dot_product(Tstar_v,lattice_Sslip_v(:,i,constitutive_dislobased_structure(matID)))
   if (abs(tau_slip) > state(ipc,ip,el)%p(3*n+i)) then
      gdot_slip = state(ipc,ip,el)%p(7*n+i)*sign(1.0_pReal,tau_slip)*&
                  sinh(((abs(tau_slip)-state(ipc,ip,el)%p(3*n+i))*state(ipc,ip,el)%p(5*n+i))/(kB*Temperature))

      locks     = (sqrt(state(ipc,ip,el)%p(n+i))*abs(gdot_slip))/&
                  (constitutive_dislobased_c4(matID)*constitutive_dislobased_bg(matID))
				   
      athermal_recovery = constitutive_dislobased_c7(matID)*state(ipc,ip,el)%p(i)*abs(gdot_slip)
      
	  !thermal_recovery  = constitutive_dislobased_c8(matID)*abs(tau_slip)*state(ipc,ip,el)%p(i)**(2.0_pReal)*&
      !                    ((constitutive_dislobased_D0(matID)*constitutive_dislobased_bg(matID)**(3.0_pReal))/&
      !                    (kB*Temperature))*exp(-constitutive_dislobased_Qsd(matID)/(kB*Temperature))

      constitutive_dislobased_dotState(i) = locks - athermal_recovery
   endif
 enddo

 return
end function


!****************************************************************
!* calculates the rate of change of temperature                 *
!****************************************************************
pure function constitutive_dislobased_dotTemperature(Tstar_v,Temperature,state,ipc,ip,el)

  !*** variables and functions from other modules ***!
  use prec,     only: pReal,pInt,p_vec
  use lattice,  only: lattice_Sslip_v
  use mesh,     only: mesh_NcpElems,mesh_maxNips
  use material, only: homogenization_maxNgrains,material_phase,phase_constitutionInstance  
  implicit none

  !*** input variables ***!
  real(pReal), dimension(6), intent(in) ::  Tstar_v                   ! 2nd Piola Kirchhoff stress tensor in Mandel notation
  real(pReal), intent(in) ::                Temperature
  integer(pInt), intent(in)::               ipc, &                    ! grain number
                                            ip, &                     ! integration point number
                                            el                        ! element number
  type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state ! state of the current microstructure
  
  !*** output variables ***!
  real(pReal) constitutive_dislobased_dotTemperature                  ! rate of change of temparature
  
  ! calculate dotTemperature
  constitutive_dislobased_dotTemperature = 0.0_pReal

  return
endfunction



pure function constitutive_dislobased_postResults(Tstar_v,Temperature,dt,state,ipc,ip,el)
!*********************************************************************
!* return array of constitutive results                              *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - dt              : current time increment                       *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec, only: pReal,pInt,p_vec
 use math, only: math_mul6x6
 use lattice, only: lattice_Sslip_v
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase,phase_constitutionInstance,phase_Noutput
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: ipc,ip,el
 real(pReal), intent(in) :: dt,Temperature
 real(pReal), dimension(6), intent(in) :: Tstar_v
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
 integer(pInt) matID,o,i,c,n
 real(pReal) tau_slip, active_rate
 real(pReal), dimension(constitutive_dislobased_sizePostResults(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_dislobased_postResults

 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 n = constitutive_dislobased_Nslip(matID)
 c = 0_pInt
 constitutive_dislobased_postResults = 0.0_pReal

 do o = 1,phase_Noutput(material_phase(ipc,ip,el))
   select case(constitutive_dislobased_output(o,matID))
     case ('dislodensity')
       constitutive_dislobased_postResults(c+1:c+n) = state(ipc,ip,el)%p(1:n)
       c = c + n
     case ('rateofshear')
       do i = 1,n
         tau_slip = math_mul6x6(Tstar_v,lattice_Sslip_v(:,i,constitutive_dislobased_structure(matID)))
		 if ((abs(tau_slip)-state(ipc,ip,el)%p(3*n+i))>0) then
            constitutive_dislobased_postResults(c+i) = state(ipc,ip,el)%p(7*n+i)*sign(1.0_pReal,tau_slip)*&
            sinh(((abs(tau_slip)-state(ipc,ip,el)%p(3*n+i))*state(ipc,ip,el)%p(5*n+i))/(kB*Temperature))
	     else
		    constitutive_dislobased_postResults(c+i) = 0.0_pReal
		 endif
       enddo
       c = c + n
   end select
 enddo
 
 return

end function

END MODULE