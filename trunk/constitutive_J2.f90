
!*****************************************************
!*      Module: CONSTITUTIVE_J2        				 *
!*****************************************************
!* contains:                                         *
!* - constitutive equations                          *
!* - parameters definition                           *
!*****************************************************

MODULE constitutive_j2

!*** Include other modules ***
 use prec, only: pReal,pInt
 implicit none

 character (len=*), parameter :: constitutive_j2_label = 'j2'
 
 integer(pInt),   dimension(:),     allocatable :: constitutive_j2_sizeDotState, &
                                                   constitutive_j2_sizeState, &
                                                   constitutive_j2_sizePostResults
 character(len=64), dimension(:,:), allocatable :: constitutive_j2_output
 integer(pInt),   dimension(:),     allocatable :: constitutive_j2_structure
 real(pReal), dimension(:),     allocatable :: constitutive_j2_mu
 real(pReal), dimension(:),     allocatable :: constitutive_j2_lambda
 real(pReal), dimension(:,:,:), allocatable :: constitutive_j2_Cslip_66
!* Visco-plastic constitutive_j2 parameters
 real(pReal), dimension(:),     allocatable :: constitutive_j2_fTaylor
 real(pReal), dimension(:),     allocatable :: constitutive_j2_s0_slip
 real(pReal), dimension(:),     allocatable :: constitutive_j2_gdot0_slip
 real(pReal), dimension(:),     allocatable :: constitutive_j2_n_slip
 real(pReal), dimension(:),     allocatable :: constitutive_j2_h0
 real(pReal), dimension(:),     allocatable :: constitutive_j2_s_sat
 real(pReal), dimension(:),     allocatable :: constitutive_j2_w0


CONTAINS
!****************************************
!* - constitutive_init
!* - constitutive_homogenizedC
!* - constitutive_microstructure
!* - constitutive_LpAndItsTangent
!* - consistutive_dotState
!* - consistutive_postResults
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
 allocate(constitutive_j2_mu(maxNinstance)) ;             constitutive_j2_mu = 0.0_pReal
 allocate(constitutive_j2_lambda(maxNinstance)) ;         constitutive_j2_lambda = 0.0_pReal
 allocate(constitutive_j2_Cslip_66(6,6,maxNinstance)) ;   constitutive_j2_Cslip_66 = 0.0_pReal
 allocate(constitutive_j2_fTaylor(maxNinstance)) ;        constitutive_j2_fTaylor = 0.0_pReal
 allocate(constitutive_j2_s0_slip(maxNinstance)) ;        constitutive_j2_s0_slip = 0.0_pReal
 allocate(constitutive_j2_gdot0_slip(maxNinstance)) ;     constitutive_j2_gdot0_slip = 0.0_pReal
 allocate(constitutive_j2_n_slip(maxNinstance)) ;         constitutive_j2_n_slip = 0.0_pReal
 allocate(constitutive_j2_h0(maxNinstance)) ;             constitutive_j2_h0 = 0.0_pReal
 allocate(constitutive_j2_s_sat(maxNinstance)) ;          constitutive_j2_s_sat = 0.0_pReal
 allocate(constitutive_j2_w0(maxNinstance)) ;             constitutive_j2_w0 = 0.0_pReal
 allocate(constitutive_j2_latent(maxNinstance)) ;         constitutive_j2_latent = 1.0_pReal
 
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
     i = phase_constitutionInstance(section)     ! which instance of my constitution is present phase
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     select case(tag)
       case ('(output)')
         output = output + 1
         constitutive_j2_output(output,i) = IO_lc(IO_stringValue(line,positions,2))
       case ('nslip')
              constitutive_j2_Nslip(i) = IO_intValue(line,positions,2)
       case ('mu')
              constitutive_j2_mu(i) = IO_floatValue(line,positions,2)
       case ('lambda')
              constitutive_j2_lambda(i) = IO_floatValue(line,positions,2)
       case ('s0_slip')
              constitutive_j2_s0_slip(i) = IO_floatValue(line,positions,2)
       case ('gdot0_slip')
              constitutive_j2_gdot0_slip(i) = IO_floatValue(line,positions,2)
       case ('n_slip')
              constitutive_j2_n_slip(i) = IO_floatValue(line,positions,2)
       case ('h0')
              constitutive_j2_h0(i) = IO_floatValue(line,positions,2)
       case ('s_sat')
              constitutive_j2_s_sat(i) = IO_floatValue(line,positions,2)
       case ('w0')
              constitutive_j2_w0(i) = IO_floatValue(line,positions,2)
       case ('taylor')
              constitutive_j2_fTaylor(i) = IO_floatValue(line,positions,2)
     end select
   endif
 enddo

100 do i = 1,maxNinstance                                        ! sanity checks
   if (constitutive_j2_s0_slip(i) < 0.0_pReal)     call IO_error(203)
   if (constitutive_j2_gdot0_slip(i) <= 0.0_pReal) call IO_error(204)
   if (constitutive_j2_n_slip(i) <= 0.0_pReal)     call IO_error(205)
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
       constitutive_j2_Cslip_66(k,j,i) =   constitutive_j2_C12(i)
       constitutive_j2_Cslip_66(k,k,i) =     constitutive_j2_C11(i)
       constitutive_j2_Cslip_66(k+3,k+3,i) = constitutive_j2_C44(i)
   end forall
   constitutive_j2_Cslip_66(:,:,i) = &
     math_Mandel3333to66(math_Voigt66to3333(constitutive_j2_Cslip_66(:,:,i)))

 enddo

 return

end subroutine


function constitutive_j2_stateInit(ipc,ip,el)
!*********************************************************************
!* initial microstructural state                                     *
!*********************************************************************
 use prec, only: pReal,pInt
 use material, only: material_phase, phase_constitutionInstance
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: ipc,ip,el
 integer(pInt) matID
 real(pReal), dimension(1) :: &
   constitutive_j2_stateInit

 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 constitutive_j2_stateInit = constitutive_j2_s0_slip(matID)

 return
end function


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

end function


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

end subroutine


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
 real(pReal), dimension(3,3) :: Lp
 real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
 real(pReal), dimension(9,9) :: dLp_dTstar
 real(pReal), norm_Tstar
 matID = phase_constitutionInstance(material_phase(ipc,ip,el))

 norm_Tstar = dsqrt(math_mul6x6(Tstar_v,Tstar_v))
 
!* Calculation of Lp
 Lp = math_Mandel6to33(Tstar_v)/norm_Tstar*constitutive_j2_gdot0_slip(matID)/constitutive_j2_fTaylor(matID)* &
      (dsqrt(1.5_pReal)/constitutive_j2_fTaylor(matID)*norm_Tstar/state(ipc,ip,el)%p(1))**constitutive_j2_n_slip(matID)

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
 integer(pInt) matID,i,n
 real(pReal) Temperature,tau_slip,gdot_slip
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(1) :: &
   constitutive_phenomenological_dotState,self_hardening
 real(pReal), norm_Tstar

 matID = phase_constitutionInstance(material_phase(ipc,ip,el))

 norm_Tstar = dsqrt(math_mul6x6(Tstar_v,Tstar_v))
 constitutive_j2_dotState = constitutive_j2_gdot0_slip(matID)/constitutive_j2_fTaylor(matID)* &
                            (dsqrt(1.5_pReal)/constitutive_j2_fTaylor(matID)*norm_Tstar/state(ipc,ip,el)%p(1))** &
							constitutive_j2_n_slip(matID) * &
							constitutive_j2_h0(matID)*(1.0_pReal-state(ipc,ip,el)%p(1)/constitutive_j2_s0_sat(matID))** &
							constitutive_j2_w0(matID)

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
   constitutive_phenomenological_postResults

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
	   constitutive_j2_postResults(c+1) = constitutive_j2_gdot0_slip(matID)/constitutive_j2_fTaylor(matID)* &
                                          (dsqrt(1.5_pReal)/constitutive_j2_fTaylor(matID)*norm_Tstar/state(ipc,ip,el)%p(1))** &
										  constitutive_j2_n_slip(matID)
       c = c + 1
   end select
 enddo
 
 return

end function

END MODULE
