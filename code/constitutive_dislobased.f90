!* $Id$
!************************************
!*      Module: CONSTITUTIVE        *
!************************************
!* contains:                        *
!* - constitutive equations         *
!* - parameters definition          *
!* - orientations                   *
!************************************

!	[TWIP steel FeMnC]

!	C11		175.0e9		# elastic constants in Pa
!	C12		115.0e9
!	C44		135.0e9
!	lattice_structure	fcc
!	Nslip			12
!	Ntwin			12
!	constitution		dislobased
!	(output)		dislocationdensity
!	(output)		shearrate_slip
!	(output)		mfp_slip				# mean free path
!	(output)		resolvedstress_slip
!	(output)		resistance_slip			# passing stress
!	(output)		volumefraction
!	(output)		shearrate_twin
!	(output)		mfp_twin				# mean free path
!	(output)		resolvedstress_twin
!	(output)		resistance_twin			# "nucleation barrier"

!	### dislocation density-based constitutive parameters ###
!	burgers			2.56e-10		# Burgers vector [m]
!	Qedge			5.5e-19			# Activation energy for dislocation glide [J/K] (0.5*G*b^3)
!	grainsize		2.0e-5			# Average grain size [m]
!	stacksize		5.0e-8			# Twin stack mean thickness [m]	
!	interaction_slipslip	1.0 2.2 3.0 1.6 3.8 4.5	# Dislocation interaction coefficients
!	interaction_sliptwin	0.0 1.0                 # Dislocation interaction coefficients
!	interaction_twintwin	0.0 1.0                 # Dislocation interaction coefficients
!	# dislocation glide
!	rho0			2.5e12			# Initial dislocation density [m/m³]
!	Cmfpslip		1.0			# Adjustable parameter controlling dislocation mean free path
!	Cactivolume		1.0			# Adjustable parameter controlling activation volume
!	Cthresholdslip		0.1			# Adjustable parameter controlling threshold stress for dislocation motion
!	Cstorage		0.02			# Adjustable parameter controlling dislocation storage
!	Carecovery		15.0			# Adjustable parameter controlling athermal recovery
!	# mechanical twinning
!	Ndot0			0.0			# Number of potential twin source per volume per time [1/m³.s] 
!	fmax			1.0			# Maximum admissible twin volume fraction
!	Cmfptwin		1.0			# Adjustable parameter controlling twin mean free path
!	Cthresholdtwin		1.0			# Adjustable parameter controlling threshold stress for deformation twinning


MODULE constitutive_dislobased

!*** Include other modules ***
 use prec, only: pReal,pInt
 implicit none

 character (len=*), parameter :: constitutive_dislobased_label = 'dislobased'

 integer(pInt),     dimension(:),           allocatable        :: constitutive_dislobased_sizeDotState, &
                                                                  constitutive_dislobased_sizeState, &
                                                                  constitutive_dislobased_sizePostResults 
 integer(pInt),     dimension(:,:),         allocatable,target :: constitutive_dislobased_sizePostResult
 character(len=64), dimension(:,:),         allocatable,target :: constitutive_dislobased_output
 
 character(len=32), dimension(:),           allocatable        :: constitutive_dislobased_structureName
 integer(pInt),     dimension(:),           allocatable        :: constitutive_dislobased_structure, &
                                                                  constitutive_dislobased_totalNslip, &
                                                                  constitutive_dislobased_totalNtwin
 integer(pInt),     dimension(:,:),         allocatable        :: constitutive_dislobased_Nslip, &
                                                                  constitutive_dislobased_Ntwin, &
                                                                  constitutive_dislobased_slipFamily, &
                                                                  constitutive_dislobased_twinFamily

 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_CoverA, &
                                                           constitutive_dislobased_C11, &
                                                           constitutive_dislobased_C12, &
                                                           constitutive_dislobased_C13, &
                                                           constitutive_dislobased_C33, &
                                                           constitutive_dislobased_C44, &
                                                           constitutive_dislobased_Gmod
 real(pReal),       dimension(:,:,:),       allocatable :: constitutive_dislobased_Cslip_66
 real(pReal),       dimension(:,:,:,:),     allocatable :: constitutive_dislobased_Ctwin_66
 real(pReal),       dimension(:,:,:,:,:),   allocatable :: constitutive_dislobased_Cslip_3333
 real(pReal),       dimension(:,:,:,:,:,:), allocatable :: constitutive_dislobased_Ctwin_3333
 real(pReal),       dimension(:,:),         allocatable :: constitutive_dislobased_rho0, &
                                                           constitutive_dislobased_Burgers, &
                                                           constitutive_dislobased_Qedge, &
                                                           constitutive_dislobased_stacksize, &
                                                           constitutive_dislobased_Ndot0, &

                                                           constitutive_dislobased_interaction_slipslip, &
                                                           constitutive_dislobased_interaction_sliptwin, &
                                                           constitutive_dislobased_interaction_twinslip, &
                                                           constitutive_dislobased_interaction_twintwin
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_grainsize, &
                                                           constitutive_dislobased_fmax, &
                                                           constitutive_dislobased_Cmfpslip, &
                                                           constitutive_dislobased_Cmfptwin, &
                                                           constitutive_dislobased_Cthresholdslip, &
                                                           constitutive_dislobased_Cthresholdtwin, &
                                                           constitutive_dislobased_Cactivolume, &
                                                           constitutive_dislobased_Carecovery, &
                                                           constitutive_dislobased_Cstorage


 real(pReal),       dimension(:,:,:),       allocatable :: constitutive_dislobased_parall_interaction, &
                                                           constitutive_dislobased_forest_interaction, &
                                                           constitutive_dislobased_hardeningMatrix_sliptwin, &
                                                           constitutive_dislobased_hardeningMatrix_twinslip, &
                                                           constitutive_dislobased_hardeningMatrix_twintwin 

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
real(pReal), parameter :: Rgas = 8.314_pReal    

CONTAINS
!****************************************
!* - constitutive_init
!* - constitutive_stateInit
!* - constitutive_homogenizedC
!* - constitutive_microstructure
!* - constitutive_LpAndItsTangent
!* - consistutive_dotState
!* - constitutive_dotTemperature
!* - consistutive_postResults
!****************************************


subroutine constitutive_dislobased_init(file)
!**************************************
!*      Module initialization         *
!**************************************
 use prec,    only: pInt,pReal
 use math,    only: math_Mandel3333to66,math_Voigt66to3333,math_mul3x3
 use IO
 use material
 use lattice

 integer(pInt), intent(in) :: file
 integer(pInt), parameter :: maxNchunks = 21
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) section,maxNinstance,i,j,k,l,m,n,o,p,q,r,s,output,mySize
 character(len=64) tag
 character(len=1024) line
 real(pReal) x,y

 write(6,*)
 write(6,'(a20,a20,a12)') '<<<+-  constitutive_',constitutive_dislobased_label,' init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)

 maxNinstance = count(phase_constitution == constitutive_dislobased_label)
 if (maxNinstance == 0) return

 allocate(constitutive_dislobased_sizeDotState(maxNinstance))       ; constitutive_dislobased_sizeDotState = 0_pInt
 allocate(constitutive_dislobased_sizeState(maxNinstance))          ; constitutive_dislobased_sizeState = 0_pInt
 allocate(constitutive_dislobased_sizePostResults(maxNinstance))    ; constitutive_dislobased_sizePostResults = 0_pInt
 allocate(constitutive_dislobased_sizePostResult(maxval(phase_Noutput), &
                                                    maxNinstance))  ; constitutive_dislobased_sizePostResult = 0_pInt
 allocate(constitutive_dislobased_output(maxval(phase_Noutput), &
                                         maxNinstance))             ; constitutive_dislobased_output = ''

 allocate(constitutive_dislobased_structureName(maxNinstance))      ; constitutive_dislobased_structureName = ''
 allocate(constitutive_dislobased_structure(maxNinstance))          ; constitutive_dislobased_structure = 0_pInt
 allocate(constitutive_dislobased_Nslip(lattice_maxNslipFamily,&
                                        maxNinstance))              ; constitutive_dislobased_Nslip = 0_pInt           
 allocate(constitutive_dislobased_Ntwin(lattice_maxNtwinFamily,&
                                        maxNinstance))              ; constitutive_dislobased_Ntwin = 0_pInt
 
 allocate(constitutive_dislobased_slipFamily(lattice_maxNslip,&
                                             maxNinstance))         ; constitutive_dislobased_slipFamily = 0_pInt           
 allocate(constitutive_dislobased_twinFamily(lattice_maxNtwin,&
                                             maxNinstance))         ; constitutive_dislobased_twinFamily = 0_pInt
 
 allocate(constitutive_dislobased_totalNslip(maxNinstance))         ; constitutive_dislobased_totalNslip = 0_pInt  
 allocate(constitutive_dislobased_totalNtwin(maxNinstance))         ; constitutive_dislobased_totalNtwin = 0_pInt 

 allocate(constitutive_dislobased_CoverA(maxNinstance))             ; constitutive_dislobased_CoverA = 0.0_pReal 
 allocate(constitutive_dislobased_C11(maxNinstance))                ; constitutive_dislobased_C11 = 0.0_pReal
 allocate(constitutive_dislobased_C12(maxNinstance))                ; constitutive_dislobased_C12 = 0.0_pReal
 allocate(constitutive_dislobased_C13(maxNinstance))                ; constitutive_dislobased_C13 = 0.0_pReal
 allocate(constitutive_dislobased_C33(maxNinstance))                ; constitutive_dislobased_C33 = 0.0_pReal
 allocate(constitutive_dislobased_C44(maxNinstance))                ; constitutive_dislobased_C44 = 0.0_pReal
 allocate(constitutive_dislobased_Gmod(maxNinstance))               ; constitutive_dislobased_Gmod = 0.0_pReal
 allocate(constitutive_dislobased_Cslip_66(6,6,maxNinstance))       ; constitutive_dislobased_Cslip_66 = 0.0_pReal
 allocate(constitutive_dislobased_Cslip_3333(3,3,3,3,maxNinstance)) ; constitutive_dislobased_Cslip_3333 = 0.0_pReal

 allocate(constitutive_dislobased_rho0(lattice_maxNslipFamily, &
                                       maxNinstance))               ; constitutive_dislobased_rho0 = 0.0_pReal
 allocate(constitutive_dislobased_Burgers(lattice_maxNslipFamily, &
                                          maxNinstance))            ; constitutive_dislobased_Burgers = 0.0_pReal
 allocate(constitutive_dislobased_Qedge(lattice_maxNslipFamily, &
                                        maxNinstance))              ; constitutive_dislobased_Qedge = 0.0_pReal
 allocate(constitutive_dislobased_grainsize(maxNinstance))          ; constitutive_dislobased_grainsize = 0.0_pReal
 allocate(constitutive_dislobased_stacksize(lattice_maxNtwinFamily, &
                                            maxNinstance))          ; constitutive_dislobased_stacksize = 0.0_pReal
 allocate(constitutive_dislobased_fmax(maxNinstance))               ; constitutive_dislobased_fmax = 0.0_pReal
 allocate(constitutive_dislobased_Ndot0(lattice_maxNtwinFamily, &
                                        maxNinstance))              ; constitutive_dislobased_Ndot0 = 0.0_pReal
 allocate(constitutive_dislobased_Cmfpslip(maxNinstance))           ; constitutive_dislobased_Cmfpslip = 0.0_pReal
 allocate(constitutive_dislobased_Cmfptwin(maxNinstance))           ; constitutive_dislobased_Cmfptwin = 0.0_pReal
 allocate(constitutive_dislobased_Cthresholdslip(maxNinstance))     ; constitutive_dislobased_Cthresholdslip = 0.0_pReal
 allocate(constitutive_dislobased_Cthresholdtwin(maxNinstance))     ; constitutive_dislobased_Cthresholdtwin = 0.0_pReal
 allocate(constitutive_dislobased_Cactivolume(maxNinstance))        ; constitutive_dislobased_Cactivolume = 0.0_pReal
 allocate(constitutive_dislobased_Carecovery(maxNinstance))         ; constitutive_dislobased_Carecovery = 0.0_pReal
 allocate(constitutive_dislobased_Cstorage(maxNinstance))           ; constitutive_dislobased_Cstorage = 0.0_pReal

 allocate(constitutive_dislobased_interaction_slipslip(lattice_maxNinteraction,&
                                                       maxNinstance)) ; constitutive_dislobased_interaction_slipslip = 0.0_pReal
 allocate(constitutive_dislobased_interaction_sliptwin(lattice_maxNinteraction,&
                                                       maxNinstance)) ; constitutive_dislobased_interaction_sliptwin = 0.0_pReal
 allocate(constitutive_dislobased_interaction_twinslip(lattice_maxNinteraction,&
                                                       maxNinstance)) ; constitutive_dislobased_interaction_twinslip = 0.0_pReal
 allocate(constitutive_dislobased_interaction_twintwin(lattice_maxNinteraction,&
                                                       maxNinstance)) ; constitutive_dislobased_interaction_twintwin = 0.0_pReal

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
       case ('nslip')
              forall (j = 1:lattice_maxNslipFamily) constitutive_dislobased_Nslip(j,i) = IO_intValue(line,positions,1+j)
       case ('ntwin')
              forall (j = 1:lattice_maxNtwinFamily) constitutive_dislobased_Ntwin(j,i) = IO_intValue(line,positions,1+j)
       case ('rho0')
              forall (j = 1:lattice_maxNslipFamily) constitutive_dislobased_rho0(j,i) = IO_floatValue(line,positions,1+j)
       case ('burgers')
              forall (j = 1:lattice_maxNslipFamily) constitutive_dislobased_Burgers(j,i) = IO_floatValue(line,positions,1+j)
       case ('qedge')
              forall (j = 1:lattice_maxNslipFamily) constitutive_dislobased_Qedge(j,i) = IO_floatValue(line,positions,1+j)
       case ('grainsize')
              constitutive_dislobased_grainsize(i) = IO_floatValue(line,positions,2)
       case ('stacksize')
              forall (j = 1:lattice_maxNtwinFamily) constitutive_dislobased_stacksize(j,i) = IO_floatValue(line,positions,1+j)
       case ('fmax')
              constitutive_dislobased_fmax(i) = IO_floatValue(line,positions,2)
       case ('ndot0')
              forall (j = 1:lattice_maxNtwinFamily) constitutive_dislobased_Ndot0(j,i) = IO_floatValue(line,positions,1+j)
       case ('cmfpslip')
              constitutive_dislobased_Cmfpslip(i) = IO_floatValue(line,positions,2)
       case ('cmfptwin') 
              constitutive_dislobased_Cmfptwin(i) = IO_floatValue(line,positions,2)
       case ('cthresholdslip')
              constitutive_dislobased_Cthresholdslip(i) = IO_floatValue(line,positions,2)
       case ('cthresholdtwin') 
              constitutive_dislobased_Cthresholdtwin(i) = IO_floatValue(line,positions,2)
       case ('cactivolume') 
              constitutive_dislobased_Cactivolume(i) = IO_floatValue(line,positions,2)
       case ('carecovery') 
              constitutive_dislobased_Carecovery(i) = IO_floatValue(line,positions,2)
       case ('cstorage') 
              constitutive_dislobased_Cstorage(i) = IO_floatValue(line,positions,2)
       case ('interaction_slipslip')
              forall (j = 1:lattice_maxNinteraction) &
			    constitutive_dislobased_interaction_slipslip(j,i) = IO_floatValue(line,positions,1+j)
       case ('interaction_sliptwin')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_dislobased_interaction_sliptwin(j,i) = IO_floatValue(line,positions,1+j)
       case ('interaction_twinslip')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_dislobased_interaction_twinslip(j,i) = IO_floatValue(line,positions,1+j)
       case ('interaction_twintwin')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_dislobased_interaction_twintwin(j,i) = IO_floatValue(line,positions,1+j)
     end select
   endif
 enddo

  
100 do i = 1,maxNinstance
   constitutive_dislobased_structure(i) = lattice_initializeStructure(constitutive_dislobased_structureName(i), &
                                                                      constitutive_dislobased_CoverA(i))
   constitutive_dislobased_Nslip(:,i) = min(lattice_NslipSystem(:,constitutive_dislobased_structure(i)),&   
                                            constitutive_dislobased_Nslip(:,i))
   constitutive_dislobased_Ntwin(:,i) = min(lattice_NtwinSystem(:,constitutive_dislobased_structure(i)),&
                                            constitutive_dislobased_Ntwin(:,i))
   constitutive_dislobased_totalNslip(i) = sum(constitutive_dislobased_Nslip(:,i))
   constitutive_dislobased_totalNtwin(i) = sum(constitutive_dislobased_Ntwin(:,i))

! sanity checks (still under construction)
   if (constitutive_dislobased_structure(i) < 1 .or. &    ! sanity checks                                                                         
       constitutive_dislobased_structure(i) > 3)           call IO_error(205)
   if (any(constitutive_dislobased_rho0(:,i) < 0.0_pReal))   call IO_error(220)
   if (any(constitutive_dislobased_Burgers(:,i) <= 0.0_pReal .and. &
           constitutive_dislobased_Nslip(:,i) > 0))        call IO_error(221)
   if (any(constitutive_dislobased_Qedge(:,i) <= 0.0_pReal .and. &
           constitutive_dislobased_Nslip(:,i) > 0))        call IO_error(222)
 enddo

 allocate(constitutive_dislobased_parall_interaction(maxval(constitutive_dislobased_totalNslip),&
                                                       maxval(constitutive_dislobased_totalNslip),&
													   maxNinstance))
 allocate(constitutive_dislobased_forest_interaction(maxval(constitutive_dislobased_totalNslip),&
                                                     maxval(constitutive_dislobased_totalNslip),&
													 maxNinstance))
 allocate(constitutive_dislobased_hardeningMatrix_sliptwin(maxval(constitutive_dislobased_totalNslip),&
                                                           maxval(constitutive_dislobased_totalNtwin),&
                                                           maxNinstance))
 allocate(constitutive_dislobased_hardeningMatrix_twinslip(maxval(constitutive_dislobased_totalNtwin),&
                                                           maxval(constitutive_dislobased_totalNslip),&
                                                           maxNinstance))
 allocate(constitutive_dislobased_hardeningMatrix_twintwin(maxval(constitutive_dislobased_totalNtwin),&
                                                           maxval(constitutive_dislobased_totalNtwin),&
                                                           maxNinstance))
 constitutive_dislobased_parall_interaction       = 0.0_pReal
 constitutive_dislobased_forest_interaction       = 0.0_pReal
 constitutive_dislobased_hardeningMatrix_sliptwin = 0.0_pReal
 constitutive_dislobased_hardeningMatrix_twinslip = 0.0_pReal
 constitutive_dislobased_hardeningMatrix_twintwin = 0.0_pReal

 allocate(constitutive_dislobased_Ctwin_66(6,6,maxval(constitutive_dislobased_totalNtwin),maxNinstance))
 constitutive_dislobased_Ctwin_66 = 0.0_pReal

 allocate(constitutive_dislobased_Ctwin_3333(3,3,3,3,maxval(constitutive_dislobased_totalNtwin),maxNinstance))
 constitutive_dislobased_Ctwin_3333 = 0.0_pReal

 do i = 1,maxNinstance   
   do j = 1,maxval(phase_Noutput)
     select case(constitutive_dislobased_output(j,i))
	   case('dislocationdensity', &
	        'shearrate_slip', &
	        'mfp_slip', &
	        'resolvedstress_slip', &
	        'resistance_slip' &
	        )
		 mySize = constitutive_dislobased_totalNslip(i)
	   case('volumefraction', &
	        'shearrate_twin', &
	        'mfp_twin', &
	        'resolvedstress_twin', &
	        'resistance_twin' &
	        )
		 mySize = constitutive_dislobased_totalNtwin(i)
	   case default
		 mySize = 0_pInt
     end select

	 if (mySize > 0_pInt) then                              ! any meaningful output found                               
	   constitutive_dislobased_sizePostResult(j,i) = mySize
	   constitutive_dislobased_sizePostResults(i)  = constitutive_dislobased_sizePostResults(i) + mySize
	 endif
   enddo

   constitutive_dislobased_sizeDotState(i) =    constitutive_dislobased_totalNslip(i) +   constitutive_dislobased_totalNtwin(i)
   constitutive_dislobased_sizeState(i)    = 10*constitutive_dislobased_totalNslip(i) + 5*constitutive_dislobased_totalNtwin(i)

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
   constitutive_dislobased_Cslip_66(:,:,i) = math_Mandel3333to66(math_Voigt66to3333(constitutive_dislobased_Cslip_66(:,:,i)))
   constitutive_dislobased_Cslip_3333(:,:,:,:,i) = math_Voigt66to3333(constitutive_dislobased_Cslip_66(:,:,i))

   !* Inverse lookup of my slip system family
   l = 0_pInt
   do j = 1,lattice_maxNslipFamily
      do k = 1,constitutive_dislobased_Nslip(j,i)
         l = l + 1
         constitutive_dislobased_slipFamily(l,i) = j
   enddo; enddo
   
   !* Inverse lookup of my twin system family
   l = 0_pInt
   do j = 1,lattice_maxNtwinFamily
      do k = 1,constitutive_dislobased_Ntwin(j,i)
         l = l + 1
         constitutive_dislobased_twinFamily(l,i) = j
   enddo; enddo
   
   !* Construction of the twin elasticity matrices
   do j=1,lattice_maxNtwinFamily
      do k=1,constitutive_dislobased_Ntwin(j,i)	   
         do l=1,3 ; do m=1,3 ; do n=1,3 ; do o=1,3 ; do p=1,3 ; do q=1,3 ; do r=1,3 ; do s=1,3
	     constitutive_dislobased_Ctwin_3333(l,m,n,o,sum(constitutive_dislobased_Nslip(1:j-1,i))+k,i) = &
		 constitutive_dislobased_Ctwin_3333(l,m,n,o,sum(constitutive_dislobased_Nslip(1:j-1,i))+k,i) + &
		 constitutive_dislobased_Cslip_3333(p,q,r,s,i)*&
		 lattice_Qtwin(l,p,sum(lattice_NslipSystem(1:j-1,constitutive_dislobased_structure(i)))+k,constitutive_dislobased_structure(i))* &
		 lattice_Qtwin(m,q,sum(lattice_NslipSystem(1:j-1,constitutive_dislobased_structure(i)))+k,constitutive_dislobased_structure(i))* &
		 lattice_Qtwin(n,r,sum(lattice_NslipSystem(1:j-1,constitutive_dislobased_structure(i)))+k,constitutive_dislobased_structure(i))* &
		 lattice_Qtwin(o,s,sum(lattice_NslipSystem(1:j-1,constitutive_dislobased_structure(i)))+k,constitutive_dislobased_structure(i))
	     enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
         constitutive_dislobased_Ctwin_66(:,:,k,i) = math_Mandel3333to66(constitutive_dislobased_Ctwin_3333(:,:,:,:,k,i))
	  enddo
   enddo

   !* Construction of the hardening matrices
   !* Iteration over the systems
   do j=1,lattice_maxNslipFamily
      do k=1,constitutive_dislobased_Nslip(j,i)
         do l=1,lattice_maxNslipFamily
            do m=1,constitutive_dislobased_Nslip(l,i)
               !* Projection of the dislocation *
               x = math_mul3x3(lattice_sn(:,sum(lattice_NslipSystem(1:j-1,constitutive_dislobased_structure(i)))+k, &
			                              constitutive_dislobased_structure(i)), &
			                   lattice_st(:,sum(lattice_NslipSystem(1:l-1,constitutive_dislobased_structure(i)))+m, &
							              constitutive_dislobased_structure(i)))
               y = 1.0_pReal-x**(2.0_pReal)
               !* Interaction matrix *
               constitutive_dislobased_forest_interaction(sum(constitutive_dislobased_Nslip(1:j-1,i))+k, &
				                                          sum(constitutive_dislobased_Nslip(1:l-1,i))+m,i) = &
			   abs(x)*constitutive_dislobased_interaction_slipslip(lattice_interactionSlipSlip( &
				                        sum(lattice_NslipSystem(1:j-1,constitutive_dislobased_structure(i)))+k, &
									    sum(lattice_NslipSystem(1:l-1,constitutive_dislobased_structure(i)))+m, &
									    constitutive_dislobased_structure(i)),i)
               if (y>0.0_pReal) &
                  constitutive_dislobased_parall_interaction(sum(constitutive_dislobased_Nslip(1:j-1,i))+k, &
				                                             sum(constitutive_dislobased_Nslip(1:l-1,i))+m,i) = &
				  sqrt(y)*constitutive_dislobased_interaction_slipslip(lattice_interactionSlipSlip( &
				                        sum(lattice_NslipSystem(1:j-1,constitutive_dislobased_structure(i)))+k, &
									    sum(lattice_NslipSystem(1:l-1,constitutive_dislobased_structure(i)))+m, &
									    constitutive_dislobased_structure(i)),i)
   enddo; enddo; enddo; enddo

   do j=1,lattice_maxNslipFamily
      do k=1,constitutive_dislobased_Nslip(j,i)
         do l=1,lattice_maxNtwinFamily
            do m=1,constitutive_dislobased_Ntwin(l,i)
               constitutive_dislobased_hardeningMatrix_sliptwin(sum(constitutive_dislobased_Nslip(1:j-1,i))+k,&
                                                                sum(constitutive_dislobased_Ntwin(1:l-1,i))+m,i) = &
               constitutive_dislobased_interaction_sliptwin(lattice_interactionSlipTwin( &
                                           sum(lattice_NslipSystem(1:j-1,constitutive_dislobased_structure(i)))+k, &
                                           sum(lattice_NtwinSystem(1:l-1,constitutive_dislobased_structure(i)))+m, &
                                           constitutive_dislobased_structure(i)),i)
   enddo; enddo; enddo; enddo

   do j=1,lattice_maxNtwinFamily
      do k=1,constitutive_dislobased_Ntwin(j,i)
         do l=1,lattice_maxNslipFamily
            do m=1,constitutive_dislobased_Nslip(l,i)
               constitutive_dislobased_hardeningMatrix_twinslip(sum(constitutive_dislobased_Ntwin(1:j-1,i))+k,&
                                                                sum(constitutive_dislobased_Nslip(1:l-1,i))+m,i) = &
               constitutive_dislobased_interaction_twinslip(lattice_interactionTwinSlip( &
                                           sum(lattice_NtwinSystem(1:j-1,constitutive_dislobased_structure(i)))+k, &
                                           sum(lattice_NslipSystem(1:l-1,constitutive_dislobased_structure(i)))+m, &
                                           constitutive_dislobased_structure(i)),i)
   enddo; enddo; enddo; enddo

   do j=1,lattice_maxNtwinFamily
      do k=1,constitutive_dislobased_Ntwin(j,i)
         do l=1,lattice_maxNtwinFamily
            do m=1,constitutive_dislobased_Ntwin(l,i)
               constitutive_dislobased_hardeningMatrix_twintwin(sum(constitutive_dislobased_Ntwin(1:j-1,i))+k,&
                                                                sum(constitutive_dislobased_Ntwin(1:l-1,i))+m,i) = &
               constitutive_dislobased_interaction_twintwin(lattice_interactionTwinTwin( &
                                           sum(lattice_NtwinSystem(1:j-1,constitutive_dislobased_structure(i)))+k, &
                                           sum(lattice_NtwinSystem(1:l-1,constitutive_dislobased_structure(i)))+m, &
                                           constitutive_dislobased_structure(i)), i )
   enddo; enddo; enddo; enddo

 enddo

 return
end subroutine


function constitutive_dislobased_stateInit(myInstance)
!*********************************************************************
!* initial microstructural state                                     *
!*********************************************************************
 use prec,    only: pReal,pInt
 use lattice, only: lattice_maxNslipFamily,lattice_maxNtwinFamily
 implicit none

 !* Definition of variables
 integer(pInt), intent(in) :: myInstance
 integer(pInt) i
 real(pReal), dimension(constitutive_dislobased_sizeState(myInstance)) :: constitutive_dislobased_stateInit
 
 constitutive_dislobased_stateInit = 0.0_pReal
 
 do i = 1,lattice_maxNslipFamily
    constitutive_dislobased_stateInit(1+sum(constitutive_dislobased_Nslip(1:i-1,myInstance)) : &
                                        sum(constitutive_dislobased_Nslip(1:i  ,myInstance))) = &
    constitutive_dislobased_rho0(i,myInstance)
 enddo
  
 return
end function


!*********************************************************************
!* relevant microstructural state                                    *
!*********************************************************************
pure function constitutive_dislobased_relevantState(myInstance)

use prec,     only: pReal, &
                    pInt
implicit none

!*** input variables
integer(pInt), intent(in) ::  myInstance                      ! number specifying the current instance of the constitution

!*** output variables
real(pReal), dimension(constitutive_dislobased_sizeState(myInstance)) :: &
                              constitutive_dislobased_relevantState ! relevant state values for the current instance of this constitution

!*** local variables

constitutive_dislobased_relevantState = 1.0e-200_pReal

endfunction


function constitutive_dislobased_homogenizedC(state,ipc,ip,el)
!*********************************************************************
!* calculates homogenized elacticity matrix                          *
!*  - state           : microstructure quantities                    *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec,     only: pReal,pInt,p_vec
 use mesh,     only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase,phase_constitutionInstance
 implicit none

 !* Definition of variables
 integer(pInt), intent(in) :: ipc,ip,el
 integer(pInt) matID,ns,nt,i
 real(pReal) sumf
 real(pReal), dimension(6,6) :: constitutive_dislobased_homogenizedC
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 
 !* Shortened notation
 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 ns = constitutive_dislobased_totalNslip(matID)
 nt = constitutive_dislobased_totalNtwin(matID)

 !* Total twin volume fraction
 sumf = sum(state(ipc,ip,el)%p((ns+1):(ns+nt)))                           ! safe for nt == 0

 !* Homogenized elasticity matrix
 constitutive_dislobased_homogenizedC = (1.0_pReal-sumf)*constitutive_dislobased_Cslip_66(:,:,matID)
 do i=1,nt
   constitutive_dislobased_homogenizedC = constitutive_dislobased_homogenizedC + &
                                          state(ipc,ip,el)%p(ns+i)*constitutive_dislobased_Ctwin_66(:,:,i,matID)
 enddo 

 return
end function


subroutine constitutive_dislobased_microstructure(Temperature,state,ipc,ip,el)
!*********************************************************************
!* calculates quantities characterizing the microstructure           *
!*  - Temperature     : temperature                                  *
!*  - state           : microstructure quantities                    *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec,     only: pReal,pInt,p_vec
 use math,     only: pi
 use mesh,     only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase,phase_constitutionInstance
 use lattice,  only: lattice_interactionSlipTwin,lattice_interactionTwinTwin
 implicit none

 !* Definition of variables
 integer(pInt), intent(in) :: ipc,ip,el
 integer(pInt) matID,ns,nt,i
 real(pReal) Temperature,sumf
 real(pReal), dimension(constitutive_dislobased_totalNtwin(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: fOverStacksize
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state

 !* Shortened notation
 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 ns = constitutive_dislobased_totalNslip(matID)
 nt = constitutive_dislobased_totalNtwin(matID)
 !* State: 1           :  ns         rho_ssd
 !* State: ns+1        :  ns+nt      f
 !* State: ns+nt+1     :  2*ns+nt    rho_forest
 !* State: 2*ns+nt+1   :  3*ns+nt    rho_parallel
 !* State: 3*ns+nt+1   :  4*ns+nt    1/lambda_slip
 !* State: 4*ns+nt+1   :  5*ns+nt    1/lambda_sliptwin
 !* State: 5*ns+nt+1   :  5*ns+2*nt  1/lambda_twin
 !* State: 5*ns+2*nt+1 :  6*ns+2*nt  mfp_slip
 !* State: 6*ns+2*nt+1 :  6*ns+3*nt  mfp_twin
 !* State: 6*ns+3*nt+1 :  7*ns+3*nt  threshold_stress_slip
 !* State: 7*ns+3*nt+1 :  7*ns+4*nt  threshold_stress_twin
 !* State: 7*ns+4*nt+1 :  8*ns+4*nt  activation volume
 !* State: 8*ns+4*nt+1 :  8*ns+5*nt  twin volume
 !* State: 8*ns+5*nt+1 :  9*ns+5*nt  rho_mobile
 !* State: 9*ns+5*nt+1 : 10*ns+5*nt  initial shear rate

 !* Total twin volume fraction
 sumf = sum(state(ipc,ip,el)%p((ns+1):(ns+nt)))                          ! safe for nt == 0

 !* rescaled twin volume fraction for topology
 forall (i = 1:nt) &
   fOverStacksize(i) = state(ipc,ip,el)%p(ns+i)/constitutive_dislobased_stacksize(constitutive_dislobased_twinFamily(i,matID),matID)

 !* Forest and parallel dislocation densities
 !$OMP CRITICAL (evilmatmul)
 state(ipc,ip,el)%p((ns+nt+1):(2*ns+nt))   = &
 matmul(constitutive_dislobased_forest_interaction(1:ns,1:ns,matID),state(ipc,ip,el)%p(1:ns)) 
 state(ipc,ip,el)%p((2*ns+nt+1):(3*ns+nt)) = &
 matmul(constitutive_dislobased_parall_interaction(1:ns,1:ns,matID),state(ipc,ip,el)%p(1:ns)) 
 !$OMP END CRITICAL (evilmatmul)
 
 !* 1/mean free distance between 2 forest dislocations seen by a moving dislocation
 forall (i=1:ns) state(ipc,ip,el)%p(3*ns+nt+i) = sqrt(state(ipc,ip,el)%p(ns+nt+i))

 !* 1/mean free distance between 2 twin stacks from different systems seen by a moving dislocation
 !$OMP CRITICAL (evilmatmul)
 state(ipc,ip,el)%p((4*ns+nt+1):(5*ns+nt)) = 0.0_pReal
 if (nt > 0_pInt) state(ipc,ip,el)%p((4*ns+nt+1):(5*ns+nt)) = &
 matmul(constitutive_dislobased_hardeningMatrix_sliptwin(1:ns,1:nt,matID),fOverStacksize(1:nt))/&
 (2.0_pReal*(1.0_pReal-sumf))
 !$OMP END CRITICAL (evilmatmul)

 !* 1/mean free distance between 2 twin stacks from different systems seen by a growing twin
 !$OMP CRITICAL (evilmatmul)
 if (nt > 0_pInt) state(ipc,ip,el)%p((5*ns+nt+1):(5*ns+2*nt)) = &
 matmul(constitutive_dislobased_hardeningMatrix_twintwin(1:nt,1:nt,matID),fOverStacksize(1:nt))/&
 (2.0_pReal*(1.0_pReal-sumf))
 !$OMP END CRITICAL (evilmatmul)

 !* mean free path between 2 obstacles seen by a moving dislocation
 do i=1,ns
    if (nt > 0_pInt) then
       state(ipc,ip,el)%p(5*ns+2*nt+i) = (constitutive_dislobased_Cmfpslip(matID)*constitutive_dislobased_grainsize(matID))/&
	   (1.0_pReal+constitutive_dislobased_grainsize(matID)*&
	   (state(ipc,ip,el)%p(3*ns+nt+i)+state(ipc,ip,el)%p(4*ns+nt+i)))  
	else
       state(ipc,ip,el)%p(5*ns+i) = (constitutive_dislobased_Cmfpslip(matID)*constitutive_dislobased_grainsize(matID))/&
	   (1.0_pReal+constitutive_dislobased_grainsize(matID)*(state(ipc,ip,el)%p(3*ns+i)))
	endif
 enddo

 !* mean free path between 2 obstacles seen by a growing twin
 forall (i=1:nt) &
    state(ipc,ip,el)%p(6*ns+2*nt+i) = (constitutive_dislobased_Cmfptwin(matID)*constitutive_dislobased_grainsize(matID))/&
	(1.0_pReal+constitutive_dislobased_grainsize(matID)*state(ipc,ip,el)%p(5*ns+nt+i))     

 !* threshold stress for dislocation motion
 forall (i=1:ns) &
   state(ipc,ip,el)%p(6*ns+3*nt+i) = constitutive_dislobased_Cthresholdslip(matID)*&
   constitutive_dislobased_Burgers(constitutive_dislobased_slipFamily(i,matID),matID)*&
   constitutive_dislobased_Gmod(matID)*sqrt(state(ipc,ip,el)%p(2*ns+nt+i))

 !* threshold stress for growing twin
 forall (i=1:nt) &
   state(ipc,ip,el)%p(7*ns+3*nt+i) = constitutive_dislobased_Cthresholdtwin(matID)*(sqrt(3.0_pReal)/3.0_pReal)*(&
   (0.0002_pReal*Temperature-0.0396_pReal)/constitutive_dislobased_Burgers(constitutive_dislobased_slipFamily(i,matID),matID)+&
   (constitutive_dislobased_Burgers(constitutive_dislobased_slipFamily(i,matID),matID)*&
    constitutive_dislobased_Gmod(matID))/state(ipc,ip,el)%p(5*ns+2*nt+i))

 !* activation volume for dislocation glide
 forall (i=1:ns) &
   state(ipc,ip,el)%p(7*ns+4*nt+i) = constitutive_dislobased_Cactivolume(matID)*&
   constitutive_dislobased_Burgers(constitutive_dislobased_slipFamily(i,matID),matID)**2*state(ipc,ip,el)%p(5*ns+2*nt+i)

 !* final twin volume after growth
 forall (i=1:nt) &
   state(ipc,ip,el)%p(8*ns+4*nt+i) = (pi/6.0_pReal)*&
                                     constitutive_dislobased_stacksize(constitutive_dislobased_twinFamily(i,matID),matID)*&
                                     state(ipc,ip,el)%p(6*ns+2*nt+i)*state(ipc,ip,el)%p(6*ns+2*nt+i)
 
 !* mobile dislocation densities
 forall (i=1:ns) &
   state(ipc,ip,el)%p(8*ns+5*nt+i) = (2.0_pReal*kB*Temperature*state(ipc,ip,el)%p(2*ns+nt+i))/&
                                     (state(ipc,ip,el)%p(6*ns+3*nt+i)*state(ipc,ip,el)%p(7*ns+4*nt+i))

 !* initial shear rate for slip
 forall (i=1:ns) &
   state(ipc,ip,el)%p(9*ns+5*nt+i) = state(ipc,ip,el)%p(8*ns+5*nt+i)*&
                                     constitutive_dislobased_Burgers(constitutive_dislobased_slipFamily(i,matID),matID)*&
                                     attack_frequency*state(ipc,ip,el)%p(5*ns+2*nt+i)*&
                                     exp(-constitutive_dislobased_Qedge(constitutive_dislobased_slipFamily(i,matID),matID)/&
                                     !   --------------------
                                          (kB*Temperature))
 
end subroutine


subroutine constitutive_dislobased_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,state,ipc,ip,el)
!*********************************************************************
!* calculates plastic velocity gradient and its tangent              *
!* INPUT:                                                            *
!*  - Temperature     : temperature                                  *
!*  - state           : microstructure quantities                    *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - Lp              : plastic velocity gradient                    *
!*  - dLp_dTstar      : derivative of Lp (4th-rank tensor)           *
!*********************************************************************
 use prec,     only: pReal,pInt,p_vec
 use math,     only: math_Plain3333to99
 use mesh,     only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase,phase_constitutionInstance
 use lattice,  only: lattice_Sslip,lattice_Sslip_v,lattice_Stwin,lattice_Stwin_v,lattice_maxNslipFamily,lattice_maxNtwinFamily, &
					 lattice_NslipSystem,lattice_NtwinSystem,lattice_shearTwin

 implicit none

 !* Definition of variables
 integer(pInt) ipc,ip,el
 integer(pInt) matID,structID,ns,nt,f,i,j,k,l,m,n,index_myFamily
 real(pReal) Temperature,sumf
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(3,3) :: Lp
 real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
 real(pReal), dimension(9,9) :: dLp_dTstar
 real(pReal), dimension(constitutive_dislobased_totalNslip(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip,dgdot_dtauslip,tau_slip
 real(pReal), dimension(constitutive_dislobased_totalNtwin(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   gdot_twin,dgdot_dtautwin,tau_twin

 !* Shortened notation
 matID    = phase_constitutionInstance(material_phase(ipc,ip,el))
 structID = constitutive_dislobased_structure(matID) 
 ns = constitutive_dislobased_totalNslip(matID)
 nt = constitutive_dislobased_totalNtwin(matID)

 !* Total twin volume fraction
 sumf = sum(state(ipc,ip,el)%p((ns+1):(ns+nt)))               ! safe for nt == 0

 Lp = 0.0_pReal
 dLp_dTstar3333 = 0.0_pReal
 dLp_dTstar = 0.0_pReal

 !* Dislocation glide part
 gdot_slip = 0.0_pReal
 dgdot_dtauslip = 0.0_pReal
 j = 0_pInt
 do f = 1,lattice_maxNslipFamily                              ! loop over all slip families
    index_myFamily = sum(lattice_NslipSystem(1:f-1,structID)) ! at which index starts my family
    do i = 1,constitutive_dislobased_Nslip(f,matID)           ! process each (active) slip system in family
       j = j+1_pInt

       !* Calculation of Lp
       tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,structID))       
       if ( abs(tau_slip(j)) > state(ipc,ip,el)%p(6*ns+3*nt+j) ) then

          gdot_slip(j) = state(ipc,ip,el)%p(9*ns+5*nt+j)*sign(1.0_pReal,tau_slip(j))*&
          sinh(((abs(tau_slip(j))-state(ipc,ip,el)%p(6*ns+3*nt+j))*state(ipc,ip,el)%p(7*ns+4*nt+j))/(kB*Temperature))

          dgdot_dtauslip(j) = (state(ipc,ip,el)%p(9*ns+5*nt+j)*state(ipc,ip,el)%p(7*ns+4*nt+j))/(kB*Temperature)*&
          cosh(((abs(tau_slip(j))-state(ipc,ip,el)%p(6*ns+3*nt+j))*state(ipc,ip,el)%p(7*ns+4*nt+j))/(kB*Temperature))

       endif
       Lp = Lp + (1.0_pReal - sumf)*gdot_slip(j)*lattice_Sslip(:,:,index_myFamily+i,structID)
				   
       !* Calculation of the tangent of Lp 						   
       forall (k=1:3,l=1:3,m=1:3,n=1:3) &
          dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + dgdot_dtauslip(j)*lattice_Sslip(k,l,index_myFamily+i,structID) &
                                                                               *lattice_Sslip(m,n,index_myFamily+i,structID) 
    enddo
 enddo

 !* Mechanical twinning part
 gdot_twin = 0.0_pReal
 dgdot_dtautwin = 0.0_pReal
 j = 0_pInt
 do f = 1,lattice_maxNtwinFamily                              ! loop over all slip families
    index_myFamily = sum(lattice_NtwinSystem(1:f-1,structID)) ! at which index starts my family
    do i = 1,constitutive_dislobased_Ntwin(f,matID)           ! process each (active) slip system in family
       j = j+1_pInt

       !* Calculation of Lp
       tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,structID))        
       if ( tau_twin(j) > 0.0_pReal ) then
	       
          gdot_twin(j) = (constitutive_dislobased_fmax(matID) - sumf)*lattice_shearTwin(index_myFamily+i,structID)*&
 	                     state(ipc,ip,el)%p(8*ns+4*nt+j)*constitutive_dislobased_Ndot0(f,matID)*&
                         exp(-(state(ipc,ip,el)%p(7*ns+3*nt+j)/tau_twin(j))**10.0_pReal) 
       
	      dgdot_dtautwin(j) = (gdot_twin(j)*10.0_pReal*state(ipc,ip,el)%p(7*ns+3*nt+j)**10.0_pReal)/(tau_twin(j)**11.0_pReal)

       endif 						   
       Lp = Lp + gdot_twin(j)*lattice_Stwin(:,:,index_myFamily+i,structID)

       !* Calculation of the tangent of Lp
       forall (k=1:3,l=1:3,m=1:3,n=1:3) &
          dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + dgdot_dtautwin(j)*lattice_Stwin(k,l,index_myFamily+i,structID) &
                                                                               *lattice_Stwin(m,n,index_myFamily+i,structID)
    enddo
 enddo

 dLp_dTstar = math_Plain3333to99(dLp_dTstar3333)

 return
end subroutine


function constitutive_dislobased_dotState(Tstar_v,Temperature,state,ipc,ip,el)
!*********************************************************************
!* rate of change of microstructure                                  *
!* INPUT:                                                            *
!*  - Temperature     : temperature                                  *
!*  - state           : microstructure quantities                    *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - constitutive_dotState : evolution of state variable            *
!*********************************************************************
 use prec,     only: pReal,pInt,p_vec
 use mesh,     only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_constitutionInstance
 use lattice,  only: lattice_Sslip,lattice_Sslip_v,lattice_Stwin,lattice_Stwin_v,lattice_maxNslipFamily,lattice_maxNtwinFamily, &
                     lattice_NslipSystem,lattice_NtwinSystem,lattice_shearTwin   
 implicit none

!* Definition of variables
 integer(pInt) ipc,ip,el
 integer(pInt) matID,structID,ns,nt,f,i,j,k,index_myFamily
 real(pReal) Temperature,sumf
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(constitutive_dislobased_totalNslip(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip,tau_slip,storage,arecovery
 real(pReal), dimension(constitutive_dislobased_totalNtwin(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   gdot_twin,tau_twin
 real(pReal), dimension(constitutive_dislobased_sizeDotState(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_dislobased_dotState
   
 !* Shortened notation
 matID    = phase_constitutionInstance(material_phase(ipc,ip,el))
 structID = constitutive_dislobased_structure(matID) 
 ns = constitutive_dislobased_totalNslip(matID)
 nt = constitutive_dislobased_totalNtwin(matID)

 !* Total twin volume fraction
 sumf = sum(state(ipc,ip,el)%p((ns+1):(ns+nt)))               ! safe for nt == 0

 constitutive_dislobased_dotState = 0.0_pReal

 !* Dislocation density evolution
 gdot_slip = 0.0_pReal
 j = 0_pInt
 do f = 1,lattice_maxNslipFamily                              ! loop over all slip families
    index_myFamily = sum(lattice_NslipSystem(1:f-1,structID)) ! at which index starts my family
    do i = 1,constitutive_dislobased_Nslip(f,matID)           ! process each (active) slip system in family
       j = j+1_pInt

       !* Calculation of Lp
       tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,structID))       
       if ( abs(tau_slip(j)) > state(ipc,ip,el)%p(6*ns+3*nt+j) ) then

          gdot_slip(j) = state(ipc,ip,el)%p(9*ns+5*nt+j)*sign(1.0_pReal,tau_slip(j))* &
                         sinh(((abs(tau_slip(j))-state(ipc,ip,el)%p(6*ns+3*nt+j))*state(ipc,ip,el)%p(7*ns+4*nt+j))/(kB*Temperature))

          storage(j) = (constitutive_dislobased_Cstorage(matID)*abs(gdot_slip(j)))/&
	                   (constitutive_dislobased_Burgers(f,matID)*state(ipc,ip,el)%p(5*ns+2*nt+j))
          
          arecovery(j) = constitutive_dislobased_Carecovery(matID)*state(ipc,ip,el)%p(j)*abs(gdot_slip(j))
		  
		  constitutive_dislobased_dotState(j) = storage(j) - arecovery(j)

       endif
    enddo
 enddo

 !* Twin volume fraction evolution
 gdot_twin = 0.0_pReal
 j = 0_pInt
 do f = 1,lattice_maxNtwinFamily                              ! loop over all twin families
    index_myFamily = sum(lattice_NtwinSystem(1:f-1,structID)) ! at which index starts my family
    do i = 1,constitutive_dislobased_Ntwin(f,matID)           ! process each (active) twin system in family
       j = j+1_pInt

       !* Calculation of Lp
       tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,structID))        
       if ( tau_twin(j) > 0.0_pReal ) &	       
          constitutive_dislobased_dotState(ns+j) = (constitutive_dislobased_fmax(matID) - sumf)* &
		  lattice_shearTwin(index_myFamily+i,structID)*state(ipc,ip,el)%p(8*ns+4*nt+j)*constitutive_dislobased_Ndot0(f,matID)*&
          exp(-(state(ipc,ip,el)%p(7*ns+3*nt+j)/tau_twin(j))**10.0_pReal)  						   
    enddo
 enddo

 return
end function


function constitutive_dislobased_dotTemperature(Tstar_v,Temperature,state,ipc,ip,el)
!*********************************************************************
!* rate of change of microstructure                                  *
!* INPUT:                                                            *
!*  - Temperature     : temperature                                  *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - constitutive_dotTemperature : evolution of Temperature         *
!*********************************************************************
 use prec,     only: pReal,pInt,p_vec
 use mesh,     only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains
 implicit none

!* Definition of variables
 integer(pInt) ipc,ip,el
 real(pReal) Temperature
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6) :: Tstar_v
 real(pReal) constitutive_dislobased_dotTemperature
     
 constitutive_dislobased_dotTemperature = 0.0_pReal
   
 return
end function


pure function constitutive_dislobased_postResults(Tstar_v,Temperature,dt,state,ipc,ip,el)
!*********************************************************************
!* return array of constitutive results                              *
!* INPUT:                                                            *
!*  - Temperature     : temperature                                  *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - dt              : current time increment                       *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec,     only: pReal,pInt,p_vec
 use mesh,     only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase,phase_constitutionInstance,phase_Noutput
 use lattice,  only: lattice_Sslip_v,lattice_Stwin_v,lattice_maxNslipFamily,lattice_maxNtwinFamily, &
					 lattice_NslipSystem,lattice_NtwinSystem,lattice_shearTwin  
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: ipc,ip,el
 real(pReal), intent(in) :: dt,Temperature
 real(pReal), dimension(6), intent(in) :: Tstar_v
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
 integer(pInt) matID,structID,ns,nt,f,o,i,c,j,index_myFamily
 real(pReal) sumf,tau
 real(pReal), dimension(constitutive_dislobased_sizePostResults(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_dislobased_postResults

 !* Shortened notation
 matID    = phase_constitutionInstance(material_phase(ipc,ip,el))
 structID = constitutive_dislobased_structure(matID) 
 ns = constitutive_dislobased_totalNslip(matID)
 nt = constitutive_dislobased_totalNtwin(matID)

 !* Total twin volume fraction
 sumf = sum(state(ipc,ip,el)%p((ns+1):(ns+nt)))                                  ! safe for nt == 0

 !* Required output 
 c = 0_pInt
 constitutive_dislobased_postResults = 0.0_pReal

 do o = 1,phase_Noutput(material_phase(ipc,ip,el))
   select case(constitutive_dislobased_output(o,matID))

     case ('dislocationdensity')
       constitutive_dislobased_postResults(c+1:c+ns) = state(ipc,ip,el)%p(1:ns)
       c = c + ns

     case ('shearrate_slip')
       j = 0_pInt
       do f = 1,lattice_maxNslipFamily                                           ! loop over all slip families
          index_myFamily = sum(lattice_NslipSystem(1:f-1,structID))              ! at which index starts my family
          do i = 1,constitutive_dislobased_Nslip(f,matID)                        ! process each (active) slip system in family
             j = j + 1_pInt
             tau = dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,structID))
             if ( abs(tau) > state(ipc,ip,el)%p(6*ns+3*nt+j) ) then
                constitutive_dislobased_postResults(c+j) = state(ipc,ip,el)%p(9*ns+5*nt+j)*sign(1.0_pReal,tau)* &
                sinh(((abs(tau)-state(ipc,ip,el)%p(6*ns+3*nt+j))*state(ipc,ip,el)%p(7*ns+4*nt+j))/(kB*Temperature))
			 else 
                constitutive_dislobased_postResults(c+j) = 0.0_pReal
			 endif
       enddo ; enddo
       c = c + ns

     case ('mfp_slip')
       constitutive_dislobased_postResults(c+1:c+ns) = state(ipc,ip,el)%p((5*ns+2*nt+1):(6*ns+2*nt))
       c = c + ns

     case ('resolvedstress_slip')
       j = 0_pInt
       do f = 1,lattice_maxNslipFamily                                           ! loop over all slip families
          index_myFamily = sum(lattice_NslipSystem(1:f-1,structID))              ! at which index starts my family
          do i = 1,constitutive_dislobased_Nslip(f,matID)                        ! process each (active) slip system in family
             j = j + 1_pInt
             constitutive_dislobased_postResults(c+j) = dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,structID))
	   enddo; enddo
       c = c + ns
   
     case ('resistance_slip')
       constitutive_dislobased_postResults(c+1:c+ns) = state(ipc,ip,el)%p((6*ns+3*nt+1):(7*ns+3*nt))
       c = c + ns

     case ('volumefraction')
       constitutive_dislobased_postResults(c+1:c+nt) = state(ipc,ip,el)%p((ns+1):(ns+nt))
       c = c + nt

     case ('shearrate_twin')
	   if (nt > 0_pInt) then 
       j = 0_pInt
       do f = 1,lattice_maxNtwinFamily                                           ! loop over all slip families
          index_myFamily = sum(lattice_NtwinSystem(1:f-1,structID))              ! at which index starts my family
          do i = 1,constitutive_dislobased_Ntwin(f,matID)                        ! process each (active) slip system in family
             j = j + 1_pInt
             tau = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,structID))
             if ( tau > 0.0_pReal ) then
                constitutive_dislobased_postResults(c+j) = (constitutive_dislobased_fmax(matID) - sumf)* &
				lattice_shearTwin(index_myFamily+i,structID)*state(ipc,ip,el)%p(8*ns+4*nt+j)* &
				constitutive_dislobased_Ndot0(f,matID)*exp(-(state(ipc,ip,el)%p(7*ns+3*nt+j)/tau)**10.0_pReal)
			 else 
                constitutive_dislobased_postResults(c+j) = 0.0_pReal
			 endif                        
       enddo ; enddo
	   endif
       c = c + nt

     case ('mfp_twin')
       constitutive_dislobased_postResults(c+1:c+nt) = state(ipc,ip,el)%p((6*ns+2*nt+1):(6*ns+3*nt))
       c = c + nt

     case ('resolvedstress_twin')
	   if (nt > 0_pInt) then
       j = 0_pInt
       do f = 1,lattice_maxNtwinFamily                                           ! loop over all slip families
          index_myFamily = sum(lattice_NtwinSystem(1:f-1,structID))              ! at which index starts my family
          do i = 1,constitutive_dislobased_Ntwin(f,matID)                        ! process each (active) slip system in family
             j = j + 1_pInt
             constitutive_dislobased_postResults(c+j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,structID))
	   enddo; enddo
	   endif
       c = c + nt

     case ('resistance_twin')
       constitutive_dislobased_postResults(c+1:c+nt) = state(ipc,ip,el)%p((7*ns+3*nt+1):(7*ns+4*nt))
       c = c + nt

   end select
 enddo
 
 return
end function

END MODULE