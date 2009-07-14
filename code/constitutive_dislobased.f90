
!************************************
!*      Module: CONSTITUTIVE        *
!************************************
!* contains:                        *
!* - constitutive equations         *
!* - parameters definition          *
!* - orientations                   *
!************************************

MODULE constitutive_dislobased
!*** Include other modules ***
 use prec, only: pReal,pInt
 implicit none

 character (len=*), parameter :: constitutive_dislobased_label = 'dislobased'
 integer(pInt),     dimension(:),     allocatable :: constitutive_dislobased_sizeDotState, &
                                                     constitutive_dislobased_sizeState, &
                                                     constitutive_dislobased_sizePostResults
 character(len=64), dimension(:,:),         allocatable :: constitutive_dislobased_output
 character(len=32), dimension(:),           allocatable :: constitutive_dislobased_structureName
 integer(pInt),     dimension(:),           allocatable :: constitutive_dislobased_structure
 integer(pInt),     dimension(:),           allocatable :: constitutive_dislobased_Nslip
 integer(pInt),     dimension(:),           allocatable :: constitutive_dislobased_Ntwin
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_C11
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_C12
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_C13
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_C33
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_C44
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_Gmod
 real(pReal),       dimension(:,:,:),       allocatable :: constitutive_dislobased_Cslip_66
 real(pReal),       dimension(:,:,:,:),     allocatable :: constitutive_dislobased_Ctwin_66
 real(pReal),       dimension(:,:,:,:,:),   allocatable :: constitutive_dislobased_Cslip_3333
 real(pReal),       dimension(:,:,:,:,:,:), allocatable :: constitutive_dislobased_Ctwin_3333
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_rho0
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_bg
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_Qedge
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_grainsize
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_stacksize
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_fmax
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_Ndot0
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_Cmfpslip
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_Cmfptwin
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_Cthresholdslip
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_Cthresholdtwin
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_Cactivolume
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_Cstorage
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_Carecovery
 real(pReal),       dimension(:),           allocatable :: constitutive_dislobased_CoverA
 real(pReal),       dimension(:,:),         allocatable :: constitutive_dislobased_SlipIntCoeff
 real(pReal),       dimension(:,:,:),       allocatable :: constitutive_dislobased_Iparallel
 real(pReal),       dimension(:,:,:),       allocatable :: constitutive_dislobased_Iforest

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
!* - constitutive_dotTemperature
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
 use lattice, only: lattice_sn,lattice_st,lattice_interactionSlipSlip,lattice_initializeStructure,lattice_Qtwin,lattice_tn
 integer(pInt), intent(in) :: file
 integer(pInt), parameter :: maxNchunks = 7
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) section,maxNinstance,i,j,k,l,m,n,o,p,q,r,s,output
 character(len=64) tag
 character(len=1024) line
 real(pReal) x,y
 
 maxNinstance = count(phase_constitution == constitutive_dislobased_label)
 if (maxNinstance == 0) return

 allocate(constitutive_dislobased_sizeDotState(maxNinstance))       ; constitutive_dislobased_sizeDotState = 0_pInt
 allocate(constitutive_dislobased_sizeState(maxNinstance))          ; constitutive_dislobased_sizeState = 0_pInt
 allocate(constitutive_dislobased_sizePostResults(maxNinstance))    ; constitutive_dislobased_sizePostResults = 0_pInt
 allocate(constitutive_dislobased_structureName(maxNinstance))      ; constitutive_dislobased_structureName = ''
 allocate(constitutive_dislobased_structure(maxNinstance))          ; constitutive_dislobased_structure = 0_pInt
 allocate(constitutive_dislobased_Nslip(maxNinstance))              ; constitutive_dislobased_Nslip = 0_pInt
 allocate(constitutive_dislobased_Ntwin(maxNinstance))              ; constitutive_dislobased_Ntwin = 0_pInt
 allocate(constitutive_dislobased_C11(maxNinstance))                ; constitutive_dislobased_C11 = 0.0_pReal
 allocate(constitutive_dislobased_C12(maxNinstance))                ; constitutive_dislobased_C12 = 0.0_pReal
 allocate(constitutive_dislobased_C13(maxNinstance))                ; constitutive_dislobased_C13 = 0.0_pReal
 allocate(constitutive_dislobased_C33(maxNinstance))                ; constitutive_dislobased_C33 = 0.0_pReal
 allocate(constitutive_dislobased_C44(maxNinstance))                ; constitutive_dislobased_C44 = 0.0_pReal
 allocate(constitutive_dislobased_Gmod(maxNinstance))               ; constitutive_dislobased_Gmod = 0.0_pReal
 allocate(constitutive_dislobased_Cslip_66(6,6,maxNinstance))       ; constitutive_dislobased_Cslip_66 = 0.0_pReal
 allocate(constitutive_dislobased_Cslip_3333(3,3,3,3,maxNinstance)) ; constitutive_dislobased_Ctwin_3333 = 0.0_pReal
 allocate(constitutive_dislobased_rho0(maxNinstance))               ; constitutive_dislobased_rho0 = 0.0_pReal
 allocate(constitutive_dislobased_bg(maxNinstance))                 ; constitutive_dislobased_bg = 0.0_pReal
 allocate(constitutive_dislobased_Qedge(maxNinstance))              ; constitutive_dislobased_Qedge = 0.0_pReal
 allocate(constitutive_dislobased_grainsize(maxNinstance))          ; constitutive_dislobased_grainsize = 0.0_pReal
 allocate(constitutive_dislobased_stacksize(maxNinstance))          ; constitutive_dislobased_stacksize = 0.0_pReal
 allocate(constitutive_dislobased_fmax(maxNinstance))               ; constitutive_dislobased_fmax = 0.0_pReal
 allocate(constitutive_dislobased_Ndot0(maxNinstance))              ; constitutive_dislobased_Ndot0 = 0.0_pReal
 allocate(constitutive_dislobased_Cmfpslip(maxNinstance))           ; constitutive_dislobased_Cmfpslip = 0.0_pReal
 allocate(constitutive_dislobased_Cmfptwin(maxNinstance))           ; constitutive_dislobased_Cmfptwin = 0.0_pReal
 allocate(constitutive_dislobased_Cthresholdslip(maxNinstance))     ; constitutive_dislobased_Cthresholdslip = 0.0_pReal
 allocate(constitutive_dislobased_Cthresholdtwin(maxNinstance))     ; constitutive_dislobased_Cthresholdtwin = 0.0_pReal
 allocate(constitutive_dislobased_Cactivolume(maxNinstance))        ; constitutive_dislobased_Cactivolume = 0.0_pReal
 allocate(constitutive_dislobased_Cstorage(maxNinstance))           ; constitutive_dislobased_Cstorage = 0.0_pReal
 allocate(constitutive_dislobased_Carecovery(maxNinstance))         ; constitutive_dislobased_Carecovery = 0.0_pReal
 allocate(constitutive_dislobased_CoverA(maxNinstance))             ; constitutive_dislobased_CoverA = 0.0_pReal
 allocate(constitutive_dislobased_SlipIntCoeff(6,maxNinstance))     ; constitutive_dislobased_SlipIntCoeff = 0.0_pReal
 allocate(constitutive_dislobased_output(maxval(phase_Noutput),maxNinstance)) ; constitutive_dislobased_output = ''

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
       case ('ntwin')
              constitutive_dislobased_Ntwin(i) = IO_intValue(line,positions,2)
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
       case ('grainsize')
              constitutive_dislobased_grainsize(i) = IO_floatValue(line,positions,2)
       case ('stacksize')
              constitutive_dislobased_stacksize(i) = IO_floatValue(line,positions,2)
       case ('fmax')
              constitutive_dislobased_fmax(i) = IO_floatValue(line,positions,2)
       case ('ndot0')
              constitutive_dislobased_Ndot0(i) = IO_floatValue(line,positions,2)
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
       case ('cstorage') 
              constitutive_dislobased_Cstorage(i) = IO_floatValue(line,positions,2)
       case ('carecovery') 
              constitutive_dislobased_Carecovery(i) = IO_floatValue(line,positions,2)
       case ('interaction_coefficients') 
         forall (j=2:min(7,positions(1))) &
           constitutive_dislobased_SlipIntCoeff(j-1,i) = IO_floatValue(line,positions,j)
     end select
   endif
 enddo

  
100 do i = 1,maxNinstance
   constitutive_dislobased_structure(i) = lattice_initializeStructure(constitutive_dislobased_structureName(i), &
                                                                      constitutive_dislobased_CoverA(i))
! sanity checks
   if (constitutive_dislobased_structure(i) < 1)           call IO_error(201)
   if (constitutive_dislobased_Nslip(i) < 1)               call IO_error(202)
   if (constitutive_dislobased_rho0(i) < 0.0_pReal)        call IO_error(220)
   if (constitutive_dislobased_bg(i) <= 0.0_pReal)         call IO_error(221)
   if (constitutive_dislobased_Qedge(i) <= 0.0_pReal)      call IO_error(222)
 enddo

 allocate(constitutive_dislobased_Iparallel(maxval(constitutive_dislobased_Nslip),maxval(constitutive_dislobased_Nslip),maxNinstance))
 constitutive_dislobased_Iparallel = 0.0_pReal
 allocate(constitutive_dislobased_Iforest(maxval(constitutive_dislobased_Nslip),maxval(constitutive_dislobased_Nslip),maxNinstance))
 constitutive_dislobased_Iforest = 0.0_pReal
 allocate(constitutive_dislobased_Ctwin_66(6,6,maxval(constitutive_dislobased_Ntwin),maxNinstance))
 constitutive_dislobased_Ctwin_66 = 0.0_pReal
 allocate(constitutive_dislobased_Ctwin_3333(3,3,3,3,maxval(constitutive_dislobased_Ntwin),maxNinstance))
 constitutive_dislobased_Ctwin_3333 = 0.0_pReal

 do i = 1,maxNinstance   
   constitutive_dislobased_sizeDotState(i) = constitutive_dislobased_Nslip(i)+constitutive_dislobased_Ntwin(i)
   constitutive_dislobased_sizeState(i)    = 10*constitutive_dislobased_Nslip(i)+5*constitutive_dislobased_Ntwin(i)

   do j = 1,maxval(phase_Noutput)
     select case(constitutive_dislobased_output(j,i))
       case('state_slip')
         constitutive_dislobased_sizePostResults(i) = constitutive_dislobased_sizePostResults(i) + constitutive_dislobased_Nslip(i)
       case('rateofshear_slip')
         constitutive_dislobased_sizePostResults(i) = constitutive_dislobased_sizePostResults(i) + constitutive_dislobased_Nslip(i)
       case('mfp_slip')
         constitutive_dislobased_sizePostResults(i) = constitutive_dislobased_sizePostResults(i) + constitutive_dislobased_Nslip(i)
       case('thresholdstress_slip')
         constitutive_dislobased_sizePostResults(i) = constitutive_dislobased_sizePostResults(i) + constitutive_dislobased_Nslip(i)
       case('state_twin')
         constitutive_dislobased_sizePostResults(i) = constitutive_dislobased_sizePostResults(i) + constitutive_dislobased_Ntwin(i)
       case('rateofshear_twin')
         constitutive_dislobased_sizePostResults(i) = constitutive_dislobased_sizePostResults(i) + constitutive_dislobased_Ntwin(i)
       case('mfp_twin')
         constitutive_dislobased_sizePostResults(i) = constitutive_dislobased_sizePostResults(i) + constitutive_dislobased_Ntwin(i)
       case('thresholdstress_twin')
         constitutive_dislobased_sizePostResults(i) = constitutive_dislobased_sizePostResults(i) + constitutive_dislobased_Ntwin(i)
     end select
   enddo

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
     
   !* Construction of the twin elasticity matrices
   !* Iteration over the systems
   constitutive_dislobased_Cslip_3333(:,:,:,:,i) = math_Voigt66to3333(constitutive_dislobased_Cslip_66(:,:,i))
   do j=1,constitutive_dislobased_Ntwin(i)
      do k=1,3
      do l=1,3
      do m=1,3
      do n=1,3
         do p=1,3
	     do q=1,3
	     do r=1,3
	     do s=1,3
	        constitutive_dislobased_Ctwin_3333(k,l,m,n,j,i) = constitutive_dislobased_Ctwin_3333(k,l,m,n,j,i) + &
			constitutive_dislobased_Cslip_3333(p,q,r,s,i)*&
			lattice_Qtwin(k,p,j,i)*lattice_Qtwin(l,q,j,i)*lattice_Qtwin(m,r,j,i)*lattice_Qtwin(n,s,j,i)
	     enddo
	     enddo
	     enddo
	     enddo
      enddo
      enddo
      enddo
      enddo
      constitutive_dislobased_Ctwin_66(:,:,j,i) = math_Mandel3333to66(constitutive_dislobased_Ctwin_3333(:,:,:,:,j,i))
   enddo

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
 real(pReal), dimension(constitutive_dislobased_sizeState(myInstance)) :: constitutive_dislobased_stateInit
 
 constitutive_dislobased_stateInit = 0.0_pReal 
 constitutive_dislobased_stateInit(1:constitutive_dislobased_Nslip(myInstance)) = constitutive_dislobased_rho0(myInstance)

 return
end function


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
 use material, only: homogenization_maxNgrains,material_phase, phase_constitutionInstance
 implicit none

 !* Definition of variables
 integer(pInt), intent(in) :: ipc,ip,el
 integer(pInt) matID,ns,nt,i
 real(pReal) sumf
 real(pReal), dimension(6,6) :: constitutive_dislobased_homogenizedC
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 
 !* Shortened notation
 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 ns = constitutive_dislobased_Nslip(matID)
 nt = constitutive_dislobased_Ntwin(matID)

 !* Total twin volume fraction
 sumf = 0.0_pReal
 if (nt > 0_pInt) sumf = sum(state(ipc,ip,el)%p((ns+1):(ns+nt)))

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
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state

 !* Shortened notation
 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 ns = constitutive_dislobased_Nslip(matID)
 nt = constitutive_dislobased_Ntwin(matID)
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
 sumf = 0.0_pReal
 if (nt > 0_pInt) sumf = sum(state(ipc,ip,el)%p((ns+1):(ns+nt)))

 !* Forest and parallel dislocation densities
 !$OMP CRITICAL (evilmatmul)
 state(ipc,ip,el)%p((ns+nt+1):(2*ns+nt))   = matmul(constitutive_dislobased_Iforest  (1:ns,1:ns,matID),state(ipc,ip,el)%p(1:ns)) 
 state(ipc,ip,el)%p((2*ns+nt+1):(3*ns+nt)) = matmul(constitutive_dislobased_Iparallel(1:ns,1:ns,matID),state(ipc,ip,el)%p(1:ns)) 
 !$OMP END CRITICAL (evilmatmul)
 
 !* 1/mean free distance between 2 forest dislocations seen by a moving dislocation
 do i=1,ns
   state(ipc,ip,el)%p(3*ns+nt+i) = sqrt(state(ipc,ip,el)%p(ns+nt+i))
 enddo 

 !* 1/mean free distance between 2 twin stacks from different systems seen by a moving dislocation
 !$OMP CRITICAL (evilmatmul)
 state(ipc,ip,el)%p((4*ns+nt+1):(5*ns+nt)) = 0.0_pReal
 if (nt > 0_pInt) state(ipc,ip,el)%p((4*ns+nt+1):(5*ns+nt)) = &
 matmul(lattice_interactionSlipTwin(1:ns,1:nt,constitutive_dislobased_structure(matID)),state(ipc,ip,el)%p((ns+1):(ns+nt)))/&
 (2.0_pReal*constitutive_dislobased_stacksize(matID)*(1.0_pReal-sumf)) 
 !$OMP END CRITICAL (evilmatmul)

 !* 1/mean free distance between 2 twin stacks from different systems seen by a growing twin
 !$OMP CRITICAL (evilmatmul)
 if (nt > 0_pInt) state(ipc,ip,el)%p((5*ns+nt+1):(5*ns+2*nt)) = &
 matmul(lattice_interactionTwinTwin(1:nt,1:nt,constitutive_dislobased_structure(matID)),state(ipc,ip,el)%p((ns+1):(ns+nt)))/&
 (2.0_pReal*constitutive_dislobased_stacksize(matID)*(1.0_pReal-sumf)) 
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
 do i=1,nt
    state(ipc,ip,el)%p(6*ns+2*nt+i) = (constitutive_dislobased_Cmfptwin(matID)*constitutive_dislobased_grainsize(matID))/&
	(1.0_pReal+constitutive_dislobased_grainsize(matID)*state(ipc,ip,el)%p(5*ns+nt+i))     
 enddo

 !* threshold stress for dislocation motion
 do i=1,ns
   state(ipc,ip,el)%p(6*ns+3*nt+i) = constitutive_dislobased_Cthresholdslip(matID)*&
   constitutive_dislobased_bg(matID)*constitutive_dislobased_Gmod(matID)*sqrt(state(ipc,ip,el)%p(2*ns+nt+i))
 enddo

 !* threshold stress for growing twin
 do i=1,nt
   state(ipc,ip,el)%p(7*ns+3*nt+i) = constitutive_dislobased_Cthresholdtwin(matID)*(sqrt(3.0_pReal)/3.0_pReal)*(&
   (0.0002_pReal*Temperature-0.0396_pReal)/constitutive_dislobased_bg(matID)+&
   (constitutive_dislobased_bg(matID)*constitutive_dislobased_Gmod(matID))/state(ipc,ip,el)%p(5*ns+2*nt+i))
 enddo

 !* activation volume for dislocation glide
 do i=1,ns
   state(ipc,ip,el)%p(7*ns+4*nt+i) = constitutive_dislobased_Cactivolume(matID)*&
   constitutive_dislobased_bg(matID)*constitutive_dislobased_bg(matID)*state(ipc,ip,el)%p(5*ns+2*nt+i)
 enddo

 !* final twin volume after growth
 do i=1,nt
   state(ipc,ip,el)%p(8*ns+4*nt+i) = (pi/6.0_pReal)*constitutive_dislobased_stacksize(matID)*&
   state(ipc,ip,el)%p(6*ns+2*nt+i)*state(ipc,ip,el)%p(6*ns+2*nt+i)
 enddo
 
 !* mobile dislocation densities
 do i=1,ns
   state(ipc,ip,el)%p(8*ns+5*nt+i) = (2.0_pReal*kB*Temperature*state(ipc,ip,el)%p(2*ns+nt+i))/&
   (state(ipc,ip,el)%p(6*ns+3*nt+i)*state(ipc,ip,el)%p(7*ns+4*nt+i))
 enddo

 !* initial shear rate for slip
 do i=1,ns
   state(ipc,ip,el)%p(9*ns+5*nt+i) = state(ipc,ip,el)%p(8*ns+5*nt+i)*constitutive_dislobased_bg(matID)*attack_frequency*&
   state(ipc,ip,el)%p(5*ns+2*nt+i)*exp(-constitutive_dislobased_Qedge(matID)/(kB*Temperature))
 enddo
 
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
 use material, only: homogenization_maxNgrains,material_phase, phase_constitutionInstance
 use lattice,  only: lattice_Sslip,lattice_Stwin,lattice_Sslip_v,lattice_Stwin_v,lattice_shearTwin
 implicit none

 !* Definition of variables
 integer(pInt) ipc,ip,el
 integer(pInt) matID,i,k,l,m,n,ns,nt
 real(pReal) Temperature,sumf
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(3,3) :: Lp
 real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
 real(pReal), dimension(9,9) :: dLp_dTstar
 real(pReal), dimension(constitutive_dislobased_Nslip(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip,dgdot_dtauslip,tau_slip
 real(pReal), dimension(constitutive_dislobased_Ntwin(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   gdot_twin,dgdot_dtautwin,tau_twin

 !* Shortened notation
 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 ns = constitutive_dislobased_Nslip(matID)
 nt = constitutive_dislobased_Ntwin(matID)

 !* Total twin volume fraction
 sumf = 0.0_pReal
 if (nt > 0_pInt) sumf = sum(state(ipc,ip,el)%p((ns+1):(ns+nt)))

 !* Calculation of Lp from dislocation glide
 Lp = 0.0_pReal
 gdot_slip = 0.0_pReal
 do i = 1,ns
    tau_slip(i) = dot_product(Tstar_v,lattice_Sslip_v(:,i,constitutive_dislobased_structure(matID)))

    if ( abs(tau_slip(i)) > state(ipc,ip,el)%p(6*ns+3*nt+i) ) &
       gdot_slip(i) = state(ipc,ip,el)%p(9*ns+5*nt+i)*sign(1.0_pReal,tau_slip(i))*&
       sinh(((abs(tau_slip(i))-state(ipc,ip,el)%p(6*ns+3*nt+i))*state(ipc,ip,el)%p(7*ns+4*nt+i))/(kB*Temperature)) 
				   
    Lp = Lp + (1.0_pReal - sumf)*gdot_slip(i)*lattice_Sslip(:,:,i,constitutive_dislobased_structure(matID))
 enddo

 !* Calculation of Lp from deformation twinning
 gdot_twin = 0.0_pReal
 do i = 1,nt
    tau_twin(i) = dot_product(Tstar_v,lattice_Stwin_v(:,i,constitutive_dislobased_structure(matID)))
 
    if ( tau_twin(i) > 0.0_pReal ) &    
       gdot_twin(i) = (constitutive_dislobased_fmax(matID) - sumf)*lattice_shearTwin(i,constitutive_dislobased_structure(matID))*&
 	   state(ipc,ip,el)%p(8*ns+4*nt+i)*constitutive_dislobased_Ndot0(matID)*&
       exp(-(state(ipc,ip,el)%p(7*ns+3*nt+i)/tau_twin(i))**10.0_pReal) 
 						   
     Lp = Lp + gdot_twin(i)*lattice_Stwin(:,:,i,constitutive_dislobased_structure(matID))
 enddo

 !* Calculation of the tangent of Lp from dislocation glide
 dLp_dTstar3333 = 0.0_pReal
 dLp_dTstar = 0.0_pReal
 dgdot_dtauslip = 0.0_pReal
 do i = 1,ns

    if ( abs(tau_slip(i)) > state(ipc,ip,el)%p(6*ns+3*nt+i) ) &
       dgdot_dtauslip(i) = (state(ipc,ip,el)%p(9*ns+5*nt+i)*state(ipc,ip,el)%p(7*ns+4*nt+i))/(kB*Temperature)*&
       cosh(((abs(tau_slip(i))-state(ipc,ip,el)%p(6*ns+3*nt+i))*state(ipc,ip,el)%p(7*ns+4*nt+i))/(kB*Temperature))
						   
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
           dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
           dgdot_dtauslip(i)*lattice_Sslip(k,l,i,constitutive_dislobased_structure(matID))&
                            *lattice_Sslip(m,n,i,constitutive_dislobased_structure(matID))
 enddo

 !* Calculation of the tangent of Lp from deformation twinning
 dgdot_dtautwin = 0.0_pReal
 do i = 1,nt

    if ( tau_twin(i) > 0.0_pReal ) &    
       dgdot_dtautwin(i) = (gdot_twin(i)*10.0_pReal*state(ipc,ip,el)%p(7*ns+3*nt+i)**10.0_pReal)/(tau_twin(i)**11.0_pReal)
						   
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
           dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
           dgdot_dtautwin(i)*lattice_Stwin(k,l,i,constitutive_dislobased_structure(matID)) &
                            *lattice_Stwin(m,n,i,constitutive_dislobased_structure(matID))
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
 use lattice,  only: lattice_Sslip_v,lattice_Stwin_v
 implicit none

!* Definition of variables
 integer(pInt) ipc,ip,el
 integer(pInt) matID,i,ns,nt
 real(pReal) Temperature,sumf,tau_slip,tau_twin,gdot_slip,gdot_twin,storage,arecovery
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(constitutive_dislobased_sizeDotState(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_dislobased_dotState
   
 !* Shortened notation
 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 ns = constitutive_dislobased_Nslip(matID)
 nt = constitutive_dislobased_Ntwin(matID)

 !* Total twin volume fraction
 sumf = 0.0_pReal
 if (nt > 0_pInt) sumf = sum(state(ipc,ip,el)%p((ns+1):(ns+nt)))

 !* Dislocation density evolution
 constitutive_dislobased_dotState = 0.0_pReal
 do i = 1,ns

   tau_slip = dot_product(Tstar_v,lattice_Sslip_v(:,i,constitutive_dislobased_structure(matID)))

   if ( abs(tau_slip) > state(ipc,ip,el)%p(6*ns+3*nt+i) ) then
      gdot_slip = state(ipc,ip,el)%p(9*ns+5*nt+i)*sign(1.0_pReal,tau_slip)*&
      sinh(((abs(tau_slip)-state(ipc,ip,el)%p(6*ns+3*nt+i))*state(ipc,ip,el)%p(7*ns+4*nt+i))/(kB*Temperature) ) 

      storage = (constitutive_dislobased_Cstorage(matID)*abs(gdot_slip))/&
	  (constitutive_dislobased_bg(matID)*state(ipc,ip,el)%p(5*ns+2*nt+i))
				   
      arecovery = constitutive_dislobased_Carecovery(matID)*state(ipc,ip,el)%p(i)*abs(gdot_slip)
      
      constitutive_dislobased_dotState(i) = storage - arecovery
   endif
 enddo

 !* Twin volume fraction evolution
 do i = 1,nt
 
   tau_twin = dot_product(Tstar_v,lattice_Stwin_v(:,i,constitutive_dislobased_structure(matID)))
 
    if ( tau_twin > 0.0_pReal ) &    
       constitutive_dislobased_dotState(ns+i) = (constitutive_dislobased_fmax(matID) - sumf)*&
 	   state(ipc,ip,el)%p(8*ns+4*nt+i)*constitutive_dislobased_Ndot0(matID)*&
       exp(-(state(ipc,ip,el)%p(7*ns+3*nt+i)/tau_twin)**10.0_pReal) 
 
 enddo

 return
end function


function constitutive_dislobased_dotTemperature(Tstar_v,Temperature,state,ipc,ip,el)
!*********************************************************************
!* rate of change of microstructure                                  *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - constitutive_dotTemperature : evolution of Temperature         *
!*********************************************************************
 use prec,     only: pReal,pInt,p_vec
 use lattice,  only: lattice_Sslip_v
 use mesh,     only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_constitutionInstance
 implicit none

!* Definition of variables
 integer(pInt) ipc,ip,el
 integer(pInt) matID,i,n
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
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - dt              : current time increment                       *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec,     only: pReal,pInt,p_vec
 use lattice,  only: lattice_Sslip_v,lattice_Stwin_v,lattice_shearTwin
 use mesh,     only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase,phase_constitutionInstance,phase_Noutput
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: ipc,ip,el
 real(pReal), intent(in) :: dt,Temperature
 real(pReal), dimension(6), intent(in) :: Tstar_v
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
 integer(pInt) matID,o,i,c,ns,nt
 real(pReal) sumf,tau_slip,tau_twin
 real(pReal), dimension(constitutive_dislobased_sizePostResults(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_dislobased_postResults

 !* Shortened notation
 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 ns = constitutive_dislobased_Nslip(matID)
 nt = constitutive_dislobased_Ntwin(matID)

 !* Total twin volume fraction
 sumf = 0.0_pReal
 if (nt > 0_pInt) sumf = sum(state(ipc,ip,el)%p((ns+1):(ns+nt)))

 !* Required output 
 c = 0_pInt
 constitutive_dislobased_postResults = 0.0_pReal

 do o = 1,phase_Noutput(material_phase(ipc,ip,el))
   select case(constitutive_dislobased_output(o,matID))

     case ('state_slip')
       constitutive_dislobased_postResults(c+1:c+ns) = state(ipc,ip,el)%p(1:ns)
       c = c + ns

     case ('rateofshear_slip')
       do i = 1,ns
         tau_slip = dot_product(Tstar_v,lattice_Sslip_v(:,i,constitutive_dislobased_structure(matID)))
       
         if ( abs(tau_slip) > state(ipc,ip,el)%p(6*ns+3*nt+i) ) then
            constitutive_dislobased_postResults(c+i) = state(ipc,ip,el)%p(9*ns+5*nt+i)*sign(1.0_pReal,tau_slip)*&
            sinh(((abs(tau_slip)-state(ipc,ip,el)%p(6*ns+3*nt+i))*state(ipc,ip,el)%p(7*ns+4*nt+i))/(kB*Temperature) ) 
	     else
	   	    constitutive_dislobased_postResults(c+i) = 0.0_pReal
	     endif
       enddo
       c = c + ns

     case ('mfp_slip')
       constitutive_dislobased_postResults(c+1:c+ns) = state(ipc,ip,el)%p((5*ns+2*nt+1):(6*ns+2*nt))
       c = c + ns
   
     case ('thresholdstress_slip')
       constitutive_dislobased_postResults(c+1:c+ns) = state(ipc,ip,el)%p((6*ns+3*nt+1):(7*ns+3*nt))
       c = c + ns

     case ('state_twin')
       if (nt > 0_pInt) constitutive_dislobased_postResults(c+1:c+nt) = state(ipc,ip,el)%p((ns+1):(ns+nt))
       c = c + nt

     case ('rateofshear_twin')
	   if (nt > 0_pInt) then 
          do i = 1,nt
             tau_twin = dot_product(Tstar_v,lattice_Stwin_v(:,i,constitutive_dislobased_structure(matID)))
       
             if ( tau_twin > 0.0_pReal ) then    
                constitutive_dislobased_postResults(c+i) = (constitutive_dislobased_fmax(matID) - sumf)*&
			    lattice_shearTwin(i,constitutive_dislobased_structure(matID))*&
 	            state(ipc,ip,el)%p(8*ns+4*nt+i)*constitutive_dislobased_Ndot0(matID)*&
                exp(-(state(ipc,ip,el)%p(7*ns+3*nt+i)/tau_twin)**10.0_pReal) 
	         else
	   	        constitutive_dislobased_postResults(c+i) = 0.0_pReal
	         endif
          enddo
	   endif
       c = c + nt

     case ('mfp_twin')
       if (nt > 0_pInt) constitutive_dislobased_postResults(c+1:c+nt) = state(ipc,ip,el)%p((6*ns+2*nt+1):(6*ns+3*nt))
       c = c + nt

     case ('thresholdstress_twin')
       if (nt > 0_pInt) constitutive_dislobased_postResults(c+1:c+nt) = state(ipc,ip,el)%p((7*ns+3*nt+1):(7*ns+4*nt))
       c = c + nt

   end select
 enddo
 
 return
end function

END MODULE