
!*****************************************************
!*      Module: CONSTITUTIVE_PHENOMENOLOGICAL        *
!*****************************************************
!* contains:                                         *
!* - constitutive equations                          *
!* - parameters definition                           *
!*****************************************************

!	[Alu]
!	constitution            phenomenological
!	(output)                slipresistance
!	(output)                rateofshear
!	lattice_structure       1
!	Nslip                   12
!	
!	c11                     106.75e9
!	c12                     60.41e9
!	c44                     28.34e9
!	
!	s0_slip                 31e6
!	gdot0_slip              0.001
!	n_slip                  20
!	h0                      75e6
!	s_sat                   63e6
!	w0                      2.25
!	latent_ratio            1.4

MODULE constitutive_phenomenological

!*** Include other modules ***
 use prec, only: pReal,pInt
 implicit none

 character (len=*), parameter :: constitutive_phenomenological_label = 'phenomenological'
 
 integer(pInt),   dimension(:),     allocatable :: constitutive_phenomenological_sizeDotState, &
                                                   constitutive_phenomenological_sizeState, &
                                                   constitutive_phenomenological_sizePostResults
 character(len=64), dimension(:,:), allocatable :: constitutive_phenomenological_output

 character(len=32), dimension(:),   allocatable :: constitutive_phenomenological_structureName
 integer(pInt),   dimension(:),     allocatable :: constitutive_phenomenological_structure
 integer(pInt),   dimension(:),     allocatable :: constitutive_phenomenological_Nslip

 real(pReal), dimension(:),     allocatable :: constitutive_phenomenological_CoverA
 real(pReal), dimension(:),     allocatable :: constitutive_phenomenological_C11
 real(pReal), dimension(:),     allocatable :: constitutive_phenomenological_C12
 real(pReal), dimension(:),     allocatable :: constitutive_phenomenological_C13
 real(pReal), dimension(:),     allocatable :: constitutive_phenomenological_C33
 real(pReal), dimension(:),     allocatable :: constitutive_phenomenological_C44
 real(pReal), dimension(:,:,:), allocatable :: constitutive_phenomenological_Cslip_66
!* Visco-plastic constitutive_phenomenological parameters
 real(pReal), dimension(:),     allocatable :: constitutive_phenomenological_s0_slip
 real(pReal), dimension(:),     allocatable :: constitutive_phenomenological_gdot0_slip
 real(pReal), dimension(:),     allocatable :: constitutive_phenomenological_n_slip
 real(pReal), dimension(:),     allocatable :: constitutive_phenomenological_h0
 real(pReal), dimension(:),     allocatable :: constitutive_phenomenological_s_sat
 real(pReal), dimension(:),     allocatable :: constitutive_phenomenological_w0
 real(pReal), dimension(:),     allocatable :: constitutive_phenomenological_latent
 real(pReal), dimension(:,:,:), allocatable :: constitutive_phenomenological_HardeningMatrix



CONTAINS
!****************************************
!* - constitutive_init
!* - constitutive_stateInit
!* - constitutive_homogenizedC
!* - constitutive_microstructure
!* - constitutive_LpAndItsTangent
!* - consistutive_dotState
!* - consistutive_postResults
!****************************************


subroutine constitutive_phenomenological_init(file)
!**************************************
!*      Module initialization         *
!**************************************
 use prec, only: pInt, pReal
 use math, only: math_Mandel3333to66, math_Voigt66to3333
 use IO
 use material

 use lattice, only: lattice_initializeStructure
 integer(pInt), intent(in) :: file
 integer(pInt), parameter :: maxNchunks = 7
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) section, maxNinstance, i,j,k, output
 character(len=64) tag
 character(len=1024) line
 
 maxNinstance = count(phase_constitution == constitutive_phenomenological_label)
 if (maxNinstance == 0) return

 allocate(constitutive_phenomenological_sizeDotState(maxNinstance)) ;   constitutive_phenomenological_sizeDotState = 0_pInt
 allocate(constitutive_phenomenological_sizeState(maxNinstance)) ;      constitutive_phenomenological_sizeState = 0_pInt
 allocate(constitutive_phenomenological_sizePostResults(maxNinstance)); constitutive_phenomenological_sizePostResults = 0_pInt
 allocate(constitutive_phenomenological_output(maxval(phase_Noutput), &
                                               maxNinstance)) ;         constitutive_phenomenological_output = ''

 allocate(constitutive_phenomenological_structureName(maxNinstance)) ;  constitutive_phenomenological_structureName = ''
 allocate(constitutive_phenomenological_structure(maxNinstance)) ;      constitutive_phenomenological_structure = 0_pInt
 allocate(constitutive_phenomenological_Nslip(maxNinstance)) ;          constitutive_phenomenological_Nslip = 0_pInt

 allocate(constitutive_phenomenological_CoverA(maxNinstance))       ;   constitutive_phenomenological_CoverA = 0.0_pReal
 allocate(constitutive_phenomenological_C11(maxNinstance)) ;            constitutive_phenomenological_C11 = 0.0_pReal
 allocate(constitutive_phenomenological_C12(maxNinstance)) ;            constitutive_phenomenological_C12 = 0.0_pReal
 allocate(constitutive_phenomenological_C13(maxNinstance)) ;            constitutive_phenomenological_C13 = 0.0_pReal
 allocate(constitutive_phenomenological_C33(maxNinstance)) ;            constitutive_phenomenological_C33 = 0.0_pReal
 allocate(constitutive_phenomenological_C44(maxNinstance)) ;            constitutive_phenomenological_C44 = 0.0_pReal
 allocate(constitutive_phenomenological_Cslip_66(6,6,maxNinstance)) ;   constitutive_phenomenological_Cslip_66 = 0.0_pReal
 allocate(constitutive_phenomenological_s0_slip(maxNinstance)) ;        constitutive_phenomenological_s0_slip = 0.0_pReal
 allocate(constitutive_phenomenological_gdot0_slip(maxNinstance)) ;     constitutive_phenomenological_gdot0_slip = 0.0_pReal
 allocate(constitutive_phenomenological_n_slip(maxNinstance)) ;         constitutive_phenomenological_n_slip = 0.0_pReal
 allocate(constitutive_phenomenological_h0(maxNinstance)) ;             constitutive_phenomenological_h0 = 0.0_pReal
 allocate(constitutive_phenomenological_s_sat(maxNinstance)) ;          constitutive_phenomenological_s_sat = 0.0_pReal
 allocate(constitutive_phenomenological_w0(maxNinstance)) ;             constitutive_phenomenological_w0 = 0.0_pReal
 allocate(constitutive_phenomenological_latent(maxNinstance)) ;         constitutive_phenomenological_latent = 0.0_pReal
 
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
   if (section > 0 .and. phase_constitution(section) == constitutive_phenomenological_label) then  ! one of my sections
     i = phase_constitutionInstance(section)     ! which instance of my constitution is present phase
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     select case(tag)
       case ('(output)')
         output = output + 1
         constitutive_phenomenological_output(output,i) = IO_lc(IO_stringValue(line,positions,2))
       case ('lattice_structure')
              constitutive_phenomenological_structureName(i) = IO_lc(IO_stringValue(line,positions,2))
       case ('nslip')
              constitutive_phenomenological_Nslip(i) = IO_intValue(line,positions,2)
	     case ('covera_ratio')
              constitutive_phenomenological_CoverA(i) = IO_floatValue(line,positions,2)
       case ('c11')
              constitutive_phenomenological_C11(i) = IO_floatValue(line,positions,2)
       case ('c12')
              constitutive_phenomenological_C12(i) = IO_floatValue(line,positions,2)
       case ('c13')
              constitutive_phenomenological_C13(i) = IO_floatValue(line,positions,2)
       case ('c33')
              constitutive_phenomenological_C33(i) = IO_floatValue(line,positions,2)
       case ('c44')
              constitutive_phenomenological_C44(i) = IO_floatValue(line,positions,2)
       case ('s0_slip')
              constitutive_phenomenological_s0_slip(i) = IO_floatValue(line,positions,2)
       case ('gdot0_slip')
              constitutive_phenomenological_gdot0_slip(i) = IO_floatValue(line,positions,2)
       case ('n_slip')
              constitutive_phenomenological_n_slip(i) = IO_floatValue(line,positions,2)
       case ('h0')
              constitutive_phenomenological_h0(i) = IO_floatValue(line,positions,2)
       case ('s_sat')
              constitutive_phenomenological_s_sat(i) = IO_floatValue(line,positions,2)
       case ('w0')
              constitutive_phenomenological_w0(i) = IO_floatValue(line,positions,2)
       case ('latent_ratio') 
           constitutive_phenomenological_latent(i) = IO_floatValue(line,positions,2)
     end select
   endif
 enddo

100 do i = 1,maxNinstance

   constitutive_phenomenological_structure(i) = lattice_initializeStructure(constitutive_phenomenological_structureName(i), &
                                                                      constitutive_phenomenological_CoverA(i))                                        ! sanity checks
   if (constitutive_phenomenological_structure(i) < 1 .or. &
       constitutive_phenomenological_structure(i) > 3)           call IO_error(201)
   if (constitutive_phenomenological_Nslip(i) < 1)               call IO_error(202)
   if (constitutive_phenomenological_s0_slip(i) < 0.0_pReal)     call IO_error(203)
   if (constitutive_phenomenological_gdot0_slip(i) <= 0.0_pReal) call IO_error(204)
   if (constitutive_phenomenological_n_slip(i) <= 0.0_pReal)     call IO_error(205)
   if (constitutive_phenomenological_h0(i) <= 0.0_pReal)         call IO_error(206)
   if (constitutive_phenomenological_s_sat(i) <= 0.0_pReal)      call IO_error(207)
   if (constitutive_phenomenological_w0(i) <= 0.0_pReal)         call IO_error(208)
   if (constitutive_phenomenological_latent(i) <= 0.0_pReal)      call IO_error(209)
 enddo

 allocate(constitutive_phenomenological_hardeningMatrix(maxval(constitutive_phenomenological_Nslip),&
                                                        maxval(constitutive_phenomenological_Nslip),&
                                                        maxNinstance))

 do i = 1,maxNinstance
   constitutive_phenomenological_sizeDotState(i) = constitutive_phenomenological_Nslip(i)
   constitutive_phenomenological_sizeState(i)    = constitutive_phenomenological_Nslip(i)

   do j = 1,maxval(phase_Noutput)
     select case(constitutive_phenomenological_output(j,i))
       case('slipresistance')
         constitutive_phenomenological_sizePostResults(i) = &
         constitutive_phenomenological_sizePostResults(i) + constitutive_phenomenological_Nslip(i)
       case('rateofshear')
         constitutive_phenomenological_sizePostResults(i) = &
         constitutive_phenomenological_sizePostResults(i) + constitutive_phenomenological_Nslip(i)
     end select
   enddo

   select case (constitutive_phenomenological_structure(i))
   case(1:2) ! cubic(s)
     forall(k=1:3)
       forall(j=1:3) &
         constitutive_phenomenological_Cslip_66(k,j,i) =   constitutive_phenomenological_C12(i)
       constitutive_phenomenological_Cslip_66(k,k,i) =     constitutive_phenomenological_C11(i)
       constitutive_phenomenological_Cslip_66(k+3,k+3,i) = constitutive_phenomenological_C44(i)
     end forall
   case(3)   ! hcp
     constitutive_phenomenological_Cslip_66(1,1,i) = constitutive_phenomenological_C11(i)
     constitutive_phenomenological_Cslip_66(2,2,i) = constitutive_phenomenological_C11(i)
     constitutive_phenomenological_Cslip_66(3,3,i) = constitutive_phenomenological_C33(i)
     constitutive_phenomenological_Cslip_66(1,2,i) = constitutive_phenomenological_C12(i)
     constitutive_phenomenological_Cslip_66(2,1,i) = constitutive_phenomenological_C12(i)
     constitutive_phenomenological_Cslip_66(1,3,i) = constitutive_phenomenological_C13(i)
     constitutive_phenomenological_Cslip_66(3,1,i) = constitutive_phenomenological_C13(i)
     constitutive_phenomenological_Cslip_66(2,3,i) = constitutive_phenomenological_C13(i)
     constitutive_phenomenological_Cslip_66(3,2,i) = constitutive_phenomenological_C13(i)
     constitutive_phenomenological_Cslip_66(4,4,i) = constitutive_phenomenological_C44(i)
     constitutive_phenomenological_Cslip_66(5,5,i) = constitutive_phenomenological_C44(i)
     constitutive_phenomenological_Cslip_66(6,6,i) = 0.5_pReal*(constitutive_phenomenological_C11(i)- &
                                                                constitutive_phenomenological_C12(i))
   end select
   constitutive_phenomenological_Cslip_66(:,:,i) = &
     math_Mandel3333to66(math_Voigt66to3333(constitutive_phenomenological_Cslip_66(:,:,i)))

   constitutive_phenomenological_hardeningMatrix(:,:,i) = constitutive_phenomenological_latent(i)
   forall (j = 1:constitutive_phenomenological_Nslip(i)) &
     constitutive_phenomenological_hardeningMatrix(j,j,i) = 1.0_pReal
 enddo

 return

end subroutine


function constitutive_phenomenological_stateInit(myInstance)
!*********************************************************************
!* initial microstructural state                                     *
!*********************************************************************
 use prec, only: pReal,pInt
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: myInstance
 real(pReal), dimension(constitutive_phenomenological_Nslip(myInstance)) :: constitutive_phenomenological_stateInit

 constitutive_phenomenological_stateInit = constitutive_phenomenological_s0_slip(myInstance)

 return
end function


function constitutive_phenomenological_homogenizedC(state,ipc,ip,el)
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
 real(pReal), dimension(6,6) :: constitutive_phenomenological_homogenizedC
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 
 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 constitutive_phenomenological_homogenizedC = constitutive_phenomenological_Cslip_66(:,:,matID)

 return

end function


subroutine constitutive_phenomenological_microstructure(Temperature,state,ipc,ip,el)
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
 integer(pInt) ipc,ip,el, matID
 real(pReal) Temperature
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state

 matID = phase_constitutionInstance(material_phase(ipc,ip,el))

end subroutine


subroutine constitutive_phenomenological_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,state,ipc,ip,el)
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
 use math, only: math_Plain3333to99
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
 real(pReal), dimension(constitutive_phenomenological_Nslip(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip,dgdot_dtauslip,tau_slip

 matID = phase_constitutionInstance(material_phase(ipc,ip,el))

!* Calculation of Lp
 Lp = 0.0_pReal
 do i = 1,constitutive_phenomenological_Nslip(matID)
   tau_slip(i)  = dot_product(Tstar_v,lattice_Sslip_v(:,i,constitutive_phenomenological_structure(matID)))

   gdot_slip(i) = constitutive_phenomenological_gdot0_slip(matID)*(abs(tau_slip(i))/state(ipc,ip,el)%p(i))**&
                  constitutive_phenomenological_n_slip(matID)*sign(1.0_pReal,tau_slip(i))
   
   Lp = Lp + gdot_slip(i)*lattice_Sslip(:,:,i,constitutive_phenomenological_structure(matID))
 enddo

!* Calculation of the tangent of Lp
 dLp_dTstar3333 = 0.0_pReal
 dLp_dTstar = 0.0_pReal
 do i = 1,constitutive_phenomenological_Nslip(matID)

   dgdot_dtauslip(i) = constitutive_phenomenological_gdot0_slip(matID)*(abs(tau_slip(i))/state(ipc,ip,el)%p(i))**&
                      (constitutive_phenomenological_n_slip(matID)-1.0_pReal)*&
                      constitutive_phenomenological_n_slip(matID)/state(ipc,ip,el)%p(i)

   forall (k=1:3,l=1:3,m=1:3,n=1:3) &
          dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
             dgdot_dtauslip(i)*lattice_Sslip(k,l,i,constitutive_phenomenological_structure(matID))* &
                               lattice_Sslip(m,n,i,constitutive_phenomenological_structure(matID))
 enddo
 dLp_dTstar = math_Plain3333to99(dLp_dTstar3333)

 return
end subroutine


function constitutive_phenomenological_dotState(Tstar_v,Temperature,state,ipc,ip,el)
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
 real(pReal) Temperature,tau_slip,gdot_slip
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(constitutive_phenomenological_Nslip(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_phenomenological_dotState,self_hardening

 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 n = constitutive_phenomenological_Nslip(matID)

!* Self-Hardening of each system
 do i = 1,n

   tau_slip = dot_product(Tstar_v,lattice_Sslip_v(:,i,constitutive_phenomenological_structure(matID)))

   gdot_slip = constitutive_phenomenological_gdot0_slip(matID)*(abs(tau_slip)/state(ipc,ip,el)%p(i))**&
               constitutive_phenomenological_n_slip(matID)*sign(1.0_pReal,tau_slip)

   self_hardening(i) = constitutive_phenomenological_h0(matID)*(1.0_pReal-state(ipc,ip,el)%p(i)/&
                       constitutive_phenomenological_s_sat(matID))**constitutive_phenomenological_w0(matID)*abs(gdot_slip)
 enddo

!$OMP CRITICAL (evilmatmul)
 constitutive_phenomenological_dotState = matmul(constitutive_phenomenological_hardeningMatrix(1:n,1:n,matID),self_hardening)
!$OMP END CRITICAL (evilmatmul)
 return

end function


function constitutive_phenomenological_dotTemperature(Tstar_v,Temperature,state,ipc,ip,el)
!*********************************************************************
!* rate of change of microstructure                                  *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - constitutive_dotTemperature : evolution of temperature         *
!*********************************************************************
 use prec, only: pReal,pInt,p_vec
 use lattice, only: lattice_Sslip_v
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_constitutionInstance
 implicit none

!* Definition of variables
 integer(pInt) ipc,ip,el
 integer(pInt) matID,i,n
 real(pReal) Temperature
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6) :: Tstar_v
 real(pReal) constitutive_phenomenological_dotTemperature

 constitutive_phenomenological_dotTemperature = 0.0_pReal

 return
end function


pure function constitutive_phenomenological_postResults(Tstar_v,Temperature,dt,state,ipc,ip,el)
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
 real(pReal) tau_slip
 real(pReal), dimension(constitutive_phenomenological_sizePostResults(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_phenomenological_postResults

 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 n = constitutive_phenomenological_Nslip(matID)
 c = 0_pInt
 constitutive_phenomenological_postResults = 0.0_pReal

 do o = 1,phase_Noutput(material_phase(ipc,ip,el))
   select case(constitutive_phenomenological_output(o,matID))

     case ('slipresistance')
       constitutive_phenomenological_postResults(c+1:c+n) = state(ipc,ip,el)%p(1:n)
       c = c + n

     case ('rateofshear')
       do i = 1,n
         tau_slip = dot_product(Tstar_v,lattice_Sslip_v(:,i,constitutive_phenomenological_structure(matID)))
         constitutive_phenomenological_postResults(c+i) = sign(1.0_pReal,tau_slip)*constitutive_phenomenological_gdot0_slip(matID)*&
                                                          (abs(tau_slip)/state(ipc,ip,el)%p(i))**&
                                                           constitutive_phenomenological_n_slip(matID)
       enddo
       c = c + n

   end select
 enddo
 
 return

end function

END MODULE
