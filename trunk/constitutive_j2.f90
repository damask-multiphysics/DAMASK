
!*****************************************************
!*      Module: CONSTITUTIVE_J2        				 *
!*****************************************************
!* contains:                                         *
!* - constitutive equations                          *
!* - parameters definition                           *
!*****************************************************

!	[Alu]
!	constitution            j2
!	(output)                flowstress
!	(output)                strainrate
!	c11                     110.9e9    # (3 C11 + 2 C12 + 2 C44) / 5  ... with C44 = C11-C12 !!
!	c12                     58.34e9    # (1 C11 + 4 C12 - 1 C44) / 5
!	taylorfactor            3
!	s0                      31e6
!	gdot0                   0.001
!	n                       20
!	h0                      75e6
!	s_sat                   63e6
!	w0                      2.25

MODULE constitutive_j2

!*** Include other modules ***
 use prec, only: pReal,pInt
 implicit none

 character (len=*), parameter :: constitutive_j2_label = 'j2'
 
 integer(pInt),   dimension(:),     allocatable :: constitutive_j2_sizeDotState, &
                                                   constitutive_j2_sizeState, &
                                                   constitutive_j2_sizePostResults
 character(len=64), dimension(:,:), allocatable :: constitutive_j2_output
 real(pReal), dimension(:),     allocatable :: constitutive_j2_C11
 real(pReal), dimension(:),     allocatable :: constitutive_j2_C12
 real(pReal), dimension(:,:,:), allocatable :: constitutive_j2_Cslip_66
!* Visco-plastic constitutive_j2 parameters
 real(pReal), dimension(:),     allocatable :: constitutive_j2_fTaylor
 real(pReal), dimension(:),     allocatable :: constitutive_j2_s0
 real(pReal), dimension(:),     allocatable :: constitutive_j2_gdot0
 real(pReal), dimension(:),     allocatable :: constitutive_j2_n
 real(pReal), dimension(:),     allocatable :: constitutive_j2_h0
 real(pReal), dimension(:),     allocatable :: constitutive_j2_s_sat
 real(pReal), dimension(:),     allocatable :: constitutive_j2_w0


CONTAINS
!****************************************
!* - constitutive_j2_init
!* - constitutive_j2_stateInit
!* - constitutive_j2_homogenizedC
!* - constitutive_j2_microstructure
!* - constitutive_j2_LpAndItsTangent
!* - consistutive_j2_dotState
!* - consistutive_j2_postResults
!****************************************


subroutine constitutive_j2_init(file)
!**************************************
!*      Module initialization         *
!**************************************
 use prec, only: pInt, pReal
 use math, only: math_Mandel3333to66, math_Voigt66to3333
 use IO
 use material
 integer(pInt), intent(in) :: file
 integer(pInt), parameter :: maxNchunks = 7
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) section, maxNinstance, i,j,k,l, output
 character(len=64) tag
 character(len=1024) line
 
 maxNinstance = count(phase_constitution == constitutive_j2_label)
 if (maxNinstance == 0) return

 allocate(constitutive_j2_sizeDotState(maxNinstance)) ;   constitutive_j2_sizeDotState = 0_pInt
 allocate(constitutive_j2_sizeState(maxNinstance)) ;      constitutive_j2_sizeState = 0_pInt
 allocate(constitutive_j2_sizePostResults(maxNinstance)); constitutive_j2_sizePostResults = 0_pInt
 allocate(constitutive_j2_output(maxval(phase_Noutput), &
                                 maxNinstance)) ;         constitutive_j2_output = ''
 allocate(constitutive_j2_C11(maxNinstance)) ;            constitutive_j2_C11 = 0.0_pReal
 allocate(constitutive_j2_C12(maxNinstance)) ;            constitutive_j2_C12 = 0.0_pReal
 allocate(constitutive_j2_Cslip_66(6,6,maxNinstance)) ;   constitutive_j2_Cslip_66 = 0.0_pReal
 allocate(constitutive_j2_fTaylor(maxNinstance)) ;        constitutive_j2_fTaylor = 0.0_pReal
 allocate(constitutive_j2_s0(maxNinstance)) ;             constitutive_j2_s0 = 0.0_pReal
 allocate(constitutive_j2_gdot0(maxNinstance)) ;          constitutive_j2_gdot0 = 0.0_pReal
 allocate(constitutive_j2_n(maxNinstance)) ;              constitutive_j2_n = 0.0_pReal
 allocate(constitutive_j2_h0(maxNinstance)) ;             constitutive_j2_h0 = 0.0_pReal
 allocate(constitutive_j2_s_sat(maxNinstance)) ;          constitutive_j2_s_sat = 0.0_pReal
 allocate(constitutive_j2_w0(maxNinstance)) ;             constitutive_j2_w0 = 0.0_pReal
 
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
   if (section > 0 .and. phase_constitution(section) == constitutive_j2_label) then  ! one of my sections
     i = phase_constitutionInstance(section)              ! which instance of my constitution is present phase
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     select case(tag)
       case ('(output)')
         output = output + 1
         constitutive_j2_output(output,i) = IO_lc(IO_stringValue(line,positions,2))
       case ('c11')
              constitutive_j2_C11(i) = IO_floatValue(line,positions,2)
       case ('c12')
              constitutive_j2_C12(i) = IO_floatValue(line,positions,2)
       case ('s0')
              constitutive_j2_s0(i) = IO_floatValue(line,positions,2)
       case ('gdot0')
              constitutive_j2_gdot0(i) = IO_floatValue(line,positions,2)
       case ('n')
              constitutive_j2_n(i) = IO_floatValue(line,positions,2)
       case ('h0')
              constitutive_j2_h0(i) = IO_floatValue(line,positions,2)
       case ('s_sat')
              constitutive_j2_s_sat(i) = IO_floatValue(line,positions,2)
       case ('w0')
              constitutive_j2_w0(i) = IO_floatValue(line,positions,2)
       case ('taylorfactor')
              constitutive_j2_fTaylor(i) = IO_floatValue(line,positions,2)
     end select
   endif
 enddo

100 do i = 1,maxNinstance                                        ! sanity checks
   if (constitutive_j2_s0(i) < 0.0_pReal)          call IO_error(203)
   if (constitutive_j2_gdot0(i) <= 0.0_pReal)      call IO_error(204)
   if (constitutive_j2_n(i) <= 0.0_pReal)          call IO_error(205)
   if (constitutive_j2_h0(i) <= 0.0_pReal)         call IO_error(206)
   if (constitutive_j2_s_sat(i) <= 0.0_pReal)      call IO_error(207)
   if (constitutive_j2_w0(i) <= 0.0_pReal)         call IO_error(208)
   if (constitutive_j2_fTaylor(i) <= 0.0_pReal)    call IO_error(240)
 enddo

 do i = 1,maxNinstance
   constitutive_j2_sizeDotState(i) = 1
   constitutive_j2_sizeState(i)    = 1

   do j = 1,maxval(phase_Noutput)
     select case(constitutive_j2_output(j,i))
       case('flowstress')
         constitutive_j2_sizePostResults(i) = &
         constitutive_j2_sizePostResults(i) + 1
       case('strainrate')
         constitutive_j2_sizePostResults(i) = &
         constitutive_j2_sizePostResults(i) + 1
     end select
   enddo

   forall(k=1:3)
     forall(j=1:3) &
       constitutive_j2_Cslip_66(k,j,i) =     constitutive_j2_C12(i)
       constitutive_j2_Cslip_66(k,k,i) =     constitutive_j2_C11(i)
       constitutive_j2_Cslip_66(k+3,k+3,i) = 0.5_pReal*(constitutive_j2_C11(i)-constitutive_j2_C12(i))
   end forall
   constitutive_j2_Cslip_66(:,:,i) = &
     math_Mandel3333to66(math_Voigt66to3333(constitutive_j2_Cslip_66(:,:,i)))

 enddo

 return

endsubroutine


function constitutive_j2_stateInit(myInstance)
!*********************************************************************
!* initial microstructural state                                     *
!*********************************************************************
 use prec, only: pReal,pInt
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: myInstance
 real(pReal), dimension(1) :: constitutive_j2_stateInit

 constitutive_j2_stateInit = constitutive_j2_s0(myInstance)

 return
endfunction


function constitutive_j2_homogenizedC(state,ipc,ip,el)
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
 real(pReal), dimension(6,6) :: constitutive_j2_homogenizedC
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 
 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 constitutive_j2_homogenizedC = constitutive_j2_Cslip_66(:,:,matID)

 return

endfunction


subroutine constitutive_j2_microstructure(Temperature,state,ipc,ip,el)
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

endsubroutine


subroutine constitutive_j2_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,state,ipc,ip,el)
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
 use math, only: math_mul6x6,math_Mandel6to33,math_Plain3333to99
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
 real(pReal), dimension(3,3) :: Tstar33
 real(pReal), dimension(3,3) :: Lp
 real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
 real(pReal), dimension(9,9) :: dLp_dTstar
 real(pReal) norm_Tstar, squarenorm_Tstar, factor
 
 matID = phase_constitutionInstance(material_phase(ipc,ip,el))

 Tstar33 = math_Mandel6to33(Tstar_v)
 squarenorm_Tstar = math_mul6x6(Tstar_v,Tstar_v)
 norm_Tstar = dsqrt(squarenorm_Tstar) 

!* Initialization of Lp and dLp_dTstar
 Lp = 0.0_pReal
 dLp_dTstar = 0.0_pReal

!* for Tstar==0 both Lp and dLp_dTstar are zero (if not n==1)
 if (norm_Tstar > 0) then
   !* Calculation of Lp
   Lp = dsqrt(1.5_pReal)*Tstar33/norm_Tstar*constitutive_j2_gdot0(matID)/constitutive_j2_fTaylor(matID)* &
        (dsqrt(1.5_pReal)/constitutive_j2_fTaylor(matID)*norm_Tstar/state(ipc,ip,el)%p(1))**constitutive_j2_n(matID)

   !* Calculation of the tangent of Lp
   factor = dsqrt(1.5_pReal)*constitutive_j2_gdot0(matID)/constitutive_j2_fTaylor(matID)* &
            (dsqrt(1.5_pReal)/constitutive_j2_fTaylor(matID)/state(ipc,ip,el)%p(1))**constitutive_j2_n(matID) * &
            norm_Tstar**(constitutive_j2_n(matID)-1.0_pReal)
   dLp_dTstar3333 = 0.0_pReal
   forall (k=1:3,l=1:3,m=1:3,n=1:3) &
     dLp_dTstar3333(k,l,m,n) = Tstar33(k,l)*Tstar33(m,n) * (constitutive_j2_n(matID)-1.0_pReal)/squarenorm_Tstar
   forall (k=1:3,l=1:3) &
     dLp_dTstar3333(k,l,k,l) = dLp_dTstar3333(k,l,k,l) + 1.0_pReal
   dLp_dTstar = math_Plain3333to99(factor * dLp_dTstar3333)
 end if
 
 return

endsubroutine


function constitutive_j2_dotState(Tstar_v,Temperature,state,ipc,ip,el)
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
 use math, only: math_mul6x6
 use lattice, only: lattice_Sslip_v
 use mesh, only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase, phase_constitutionInstance
 implicit none

!* Definition of variables
 integer(pInt) ipc,ip,el
 integer(pInt) matID
 real(pReal) Temperature
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(1) :: constitutive_j2_dotState
 real(pReal) norm_Tstar

 matID = phase_constitutionInstance(material_phase(ipc,ip,el))

 norm_Tstar = dsqrt(math_mul6x6(Tstar_v,Tstar_v))
 constitutive_j2_dotState = constitutive_j2_h0(matID)*(1.0_pReal-state(ipc,ip,el)%p(1)/constitutive_j2_s_sat(matID))** &
                            constitutive_j2_w0(matID) * &
                            constitutive_j2_gdot0(matID)/constitutive_j2_fTaylor(matID)* &
                            (dsqrt(1.5_pReal)/constitutive_j2_fTaylor(matID)*norm_Tstar/state(ipc,ip,el)%p(1))** &
                            constitutive_j2_n(matID) 

 return

endfunction


pure function constitutive_j2_postResults(Tstar_v,Temperature,dt,state,ipc,ip,el)
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
 real(pReal) norm_Tstar
 real(pReal), dimension(constitutive_j2_sizePostResults(phase_constitutionInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_j2_postResults

 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 norm_Tstar = dsqrt(math_mul6x6(Tstar_v,Tstar_v))
 c = 0_pInt
 constitutive_j2_postResults = 0.0_pReal

 do o = 1,phase_Noutput(material_phase(ipc,ip,el))
   select case(constitutive_j2_output(o,matID))
     case ('flowstress')
       constitutive_j2_postResults(c+1) = state(ipc,ip,el)%p(1)
       c = c + 1
     case ('strainrate')
       constitutive_j2_postResults(c+1) = constitutive_j2_gdot0(matID)/constitutive_j2_fTaylor(matID)* &
                                          (dsqrt(1.5_pReal)/constitutive_j2_fTaylor(matID)*norm_Tstar/state(ipc,ip,el)%p(1))** &
                                          constitutive_j2_n(matID)
       c = c + 1
   end select
 enddo
 
 return

endfunction

END MODULE
