! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!##############################################################
!* $Id$
!*****************************************************
!*      Module: CONSTITUTIVE_J2                      *
!*****************************************************
!* contains:                                         *
!* - constitutive equations                          *
!* - parameters definition                           *
!*****************************************************

! [Alu]
! constitution            j2
! (output)                flowstress
! (output)                strainrate
! c11                     110.9e9    # (3 C11 + 2 C12 + 2 C44) / 5  ... with C44 = C11-C12 !!
! c12                     58.34e9    # (1 C11 + 4 C12 - 1 C44) / 5
! taylorfactor            3
! tau0                    31e6
! gdot0                   0.001
! n                       20
! h0                      75e6
! tausat                  63e6
! w0                      2.25

MODULE constitutive_j2

!*** Include other modules ***
 use prec, only: pReal,pInt
 implicit none

 character (len=*), parameter :: constitutive_j2_label = 'j2'
 
 integer(pInt),   dimension(:),     allocatable :: constitutive_j2_sizeDotState, &
                                                   constitutive_j2_sizeState, &
                                                   constitutive_j2_sizePostResults
 integer(pInt),   dimension(:,:),   allocatable,target :: constitutive_j2_sizePostResult     ! size of each post result output
 character(len=64), dimension(:,:), allocatable,target :: constitutive_j2_output             ! name of each post result output
 real(pReal), dimension(:),     allocatable :: constitutive_j2_C11
 real(pReal), dimension(:),     allocatable :: constitutive_j2_C12
 real(pReal), dimension(:,:,:), allocatable :: constitutive_j2_Cslip_66
!* Visco-plastic constitutive_j2 parameters
 real(pReal), dimension(:),     allocatable :: constitutive_j2_fTaylor
 real(pReal), dimension(:),     allocatable :: constitutive_j2_tau0
 real(pReal), dimension(:),     allocatable :: constitutive_j2_gdot0
 real(pReal), dimension(:),     allocatable :: constitutive_j2_n
 real(pReal), dimension(:),     allocatable :: constitutive_j2_h0
 real(pReal), dimension(:),     allocatable :: constitutive_j2_tausat
 real(pReal), dimension(:),     allocatable :: constitutive_j2_w0
 real(pReal), dimension(:),     allocatable :: constitutive_j2_aTolResistance


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
 use debug, only: debug_verbosity
 integer(pInt), intent(in) :: file
 integer(pInt), parameter :: maxNchunks = 7
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) section, maxNinstance, i,j,k,l, output, mySize
 character(len=64) tag
 character(len=1024) line

 !$OMP CRITICAL (write2out)
   write(6,*)
   write(6,'(a20,a20,a12)') '<<<+-  constitutive_',constitutive_j2_label,' init  -+>>>'
   write(6,*) '$Id$'
   write(6,*)
 !$OMP END CRITICAL (write2out)
 
 maxNinstance = count(phase_constitution == constitutive_j2_label)
 if (maxNinstance == 0) return

 if (debug_verbosity > 0) then
   !$OMP CRITICAL (write2out)
     write(6,'(a16,x,i5)') '# instances:',maxNinstance
     write(6,*)
   !$OMP END CRITICAL (write2out)
 endif
 
 allocate(constitutive_j2_sizeDotState(maxNinstance)) ;                         constitutive_j2_sizeDotState = 0_pInt
 allocate(constitutive_j2_sizeState(maxNinstance)) ;                            constitutive_j2_sizeState = 0_pInt
 allocate(constitutive_j2_sizePostResults(maxNinstance));                       constitutive_j2_sizePostResults = 0_pInt
 allocate(constitutive_j2_sizePostResult(maxval(phase_Noutput), maxNinstance)); constitutive_j2_sizePostResult = 0_pInt
 allocate(constitutive_j2_output(maxval(phase_Noutput), maxNinstance)) ;        constitutive_j2_output = ''
 allocate(constitutive_j2_C11(maxNinstance)) ;                                  constitutive_j2_C11 = 0.0_pReal
 allocate(constitutive_j2_C12(maxNinstance)) ;                                  constitutive_j2_C12 = 0.0_pReal
 allocate(constitutive_j2_Cslip_66(6,6,maxNinstance)) ;                         constitutive_j2_Cslip_66 = 0.0_pReal
 allocate(constitutive_j2_fTaylor(maxNinstance)) ;                              constitutive_j2_fTaylor = 0.0_pReal
 allocate(constitutive_j2_tau0(maxNinstance)) ;                                 constitutive_j2_tau0 = 0.0_pReal
 allocate(constitutive_j2_gdot0(maxNinstance)) ;                                constitutive_j2_gdot0 = 0.0_pReal
 allocate(constitutive_j2_n(maxNinstance)) ;                                    constitutive_j2_n = 0.0_pReal
 allocate(constitutive_j2_h0(maxNinstance)) ;                                   constitutive_j2_h0 = 0.0_pReal
 allocate(constitutive_j2_tausat(maxNinstance)) ;                               constitutive_j2_tausat = 0.0_pReal
 allocate(constitutive_j2_w0(maxNinstance)) ;                                   constitutive_j2_w0 = 0.0_pReal
 allocate(constitutive_j2_aTolResistance(maxNinstance)) ;                       constitutive_j2_aTolResistance = 0.0_pReal
 
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
       case ('tau0')
              constitutive_j2_tau0(i) = IO_floatValue(line,positions,2)
       case ('gdot0')
              constitutive_j2_gdot0(i) = IO_floatValue(line,positions,2)
       case ('n')
              constitutive_j2_n(i) = IO_floatValue(line,positions,2)
       case ('h0')
              constitutive_j2_h0(i) = IO_floatValue(line,positions,2)
       case ('tausat')
              constitutive_j2_tausat(i) = IO_floatValue(line,positions,2)
       case ('w0')
              constitutive_j2_w0(i) = IO_floatValue(line,positions,2)
       case ('taylorfactor')
              constitutive_j2_fTaylor(i) = IO_floatValue(line,positions,2)
       case ('atol_resistance')
              constitutive_j2_aTolResistance(i) = IO_floatValue(line,positions,2)
     end select
   endif
 enddo

100 do i = 1,maxNinstance                                        ! sanity checks
   if (constitutive_j2_tau0(i) < 0.0_pReal)               call IO_error(210)
   if (constitutive_j2_gdot0(i) <= 0.0_pReal)             call IO_error(211)
   if (constitutive_j2_n(i) <= 0.0_pReal)                 call IO_error(212)
   if (constitutive_j2_tausat(i) <= 0.0_pReal)            call IO_error(213)
   if (constitutive_j2_w0(i) <= 0.0_pReal)                call IO_error(241)
   if (constitutive_j2_fTaylor(i) <= 0.0_pReal)           call IO_error(240)
   if (constitutive_j2_aTolResistance(i) <= 0.0_pReal)    call IO_error(242)
 enddo

 do i = 1,maxNinstance
   do j = 1,maxval(phase_Noutput)
     select case(constitutive_j2_output(j,i))
       case('flowstress')
       mySize = 1_pInt
       case('strainrate')
       mySize = 1_pInt
       case default
       mySize = 0_pInt
     end select
  
     if (mySize > 0_pInt) then                               ! any meaningful output found
       constitutive_j2_sizePostResult(j,i) = mySize
       constitutive_j2_sizePostResults(i) = &
       constitutive_j2_sizePostResults(i) + mySize
     endif
   enddo

   constitutive_j2_sizeDotState(i) = 1
   constitutive_j2_sizeState(i)    = 1

   forall(k=1:3)
     forall(j=1:3) &
       constitutive_j2_Cslip_66(k,j,i) =   constitutive_j2_C12(i)
     constitutive_j2_Cslip_66(k,k,i) =     constitutive_j2_C11(i)
     constitutive_j2_Cslip_66(k+3,k+3,i) = 0.5_pReal*(constitutive_j2_C11(i)-constitutive_j2_C12(i))
   end forall
   constitutive_j2_Cslip_66(1:6,1:6,i) = &
     math_Mandel3333to66(math_Voigt66to3333(constitutive_j2_Cslip_66(1:6,1:6,i)))

 enddo

 return

endsubroutine


!*********************************************************************
!* initial microstructural state                                     *
!*********************************************************************
pure function constitutive_j2_stateInit(myInstance)
  
  use prec, only: pReal,pInt
  implicit none
  
  integer(pInt), intent(in) :: myInstance
  real(pReal), dimension(1) :: constitutive_j2_stateInit
  
  constitutive_j2_stateInit = constitutive_j2_tau0(myInstance)

  return
endfunction


!*********************************************************************
!* relevant microstructural state                                    *
!*********************************************************************
pure function constitutive_j2_aTolState(myInstance)

use prec,     only: pReal, &
                    pInt
implicit none

!*** input variables
integer(pInt), intent(in) ::  myInstance                      ! number specifying the current instance of the constitution

!*** output variables
real(pReal), dimension(constitutive_j2_sizeState(myInstance)) :: &
                              constitutive_j2_aTolState   ! relevant state values for the current instance of this constitution

!*** local variables

constitutive_j2_aTolState = constitutive_j2_aTolResistance(myInstance)

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

 integer(pInt), intent(in) :: ipc,ip,el
 integer(pInt) matID
 real(pReal), dimension(6,6) :: constitutive_j2_homogenizedC
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: state
 
 matID = phase_constitutionInstance(material_phase(ipc,ip,el))
 constitutive_j2_homogenizedC = constitutive_j2_Cslip_66(1:6,1:6,matID)

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


!****************************************************************
!* calculates plastic velocity gradient and its tangent         *
!****************************************************************
pure subroutine constitutive_j2_LpAndItsTangent(Lp, dLp_dTstar_99, Tstar_dev_v, Temperature, state, g, ip, el)

  !*** variables and functions from other modules ***!
  use prec,     only: pReal, &
                      pInt, &
                      p_vec
  use math,     only: math_mul6x6, &
                      math_Mandel6to33, &
                      math_Plain3333to99, &
                      math_spectral1
  use lattice,  only: lattice_Sslip, &
                      lattice_Sslip_v
  use mesh,     only: mesh_NcpElems, &
                      mesh_maxNips
  use material, only: homogenization_maxNgrains, &
                      material_phase, &
                      phase_constitutionInstance

  implicit none

  !*** input variables ***!
  real(pReal), dimension(6), intent(in)::       Tstar_dev_v               ! deviatoric part of the 2nd Piola Kirchhoff stress tensor in Mandel notation
  real(pReal), intent(in)::                     Temperature
  integer(pInt), intent(in)::                   g, &                      ! grain number
                                                ip, &                     ! integration point number
                                                el                        ! element number
  type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in):: state ! state of the current microstructure

  !*** output variables ***!
  real(pReal), dimension(3,3), intent(out) ::   Lp                        ! plastic velocity gradient
  real(pReal), dimension(9,9), intent(out) ::   dLp_dTstar_99             ! derivative of Lp with respect to Tstar (9x9 matrix)

  !*** local variables ***!
  real(pReal), dimension(3,3) ::                Tstar_dev_33              ! deviatoric part of the 2nd Piola Kirchhoff stress tensor as 2nd order tensor
  real(pReal), dimension(3,3,3,3) ::            dLp_dTstar_3333           ! derivative of Lp with respect to Tstar as 4th order tensor
  real(pReal)                                   gamma_dot, &              ! strainrate
                                                norm_Tstar_dev, &         ! euclidean norm of Tstar_dev
                                                squarenorm_Tstar_dev      ! square of the euclidean norm of Tstar_dev
  integer(pInt)                                 matID, &
                                                k, &
                                                l, &
                                                m, &
                                                n
 
  matID = phase_constitutionInstance(material_phase(g,ip,el))

  ! convert Tstar to matrix and calculate euclidean norm
  Tstar_dev_33 = math_Mandel6to33(Tstar_dev_v)
  squarenorm_Tstar_dev = math_mul6x6(Tstar_dev_v,Tstar_dev_v)
  norm_Tstar_dev = sqrt(squarenorm_Tstar_dev) 

  ! Initialization of Lp and dLp_dTstar
  Lp = 0.0_pReal
  dLp_dTstar_99 = 0.0_pReal

  ! for Tstar==0 both Lp and dLp_dTstar are zero (if not n==1)
  if (norm_Tstar_dev > 0) then
   
    ! Calculation of gamma_dot
    gamma_dot = constitutive_j2_gdot0(matID) * (            sqrt(1.5_pReal) * norm_Tstar_dev & 
                                               / &!---------------------------------------------------
                                                 (constitutive_j2_fTaylor(matID) * state(g,ip,el)%p(1)) ) **constitutive_j2_n(matID)

    ! Calculation of Lp
    Lp = Tstar_dev_33/norm_Tstar_dev * gamma_dot/constitutive_j2_fTaylor(matID)

    !* Calculation of the tangent of Lp
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dTstar_3333(k,l,m,n) = (constitutive_j2_n(matID)-1.0_pReal) * Tstar_dev_33(k,l)*Tstar_dev_33(m,n) / squarenorm_Tstar_dev
    forall (k=1:3,l=1:3) &
      dLp_dTstar_3333(k,l,k,l) = dLp_dTstar_3333(k,l,k,l) + 1.0_pReal
    dLp_dTstar_99 = math_Plain3333to99(gamma_dot / constitutive_j2_fTaylor(matID) * dLp_dTstar_3333 / norm_Tstar_dev)
  end if

  return

endsubroutine


!****************************************************************
!* calculates the rate of change of microstructure              *
!****************************************************************
pure function constitutive_j2_dotState(Tstar_v, Temperature, state, g, ip, el)

  !*** variables and functions from other modules ***!
  use prec,     only: pReal, &
                      pInt, &
                      p_vec
  use math,     only: math_mul6x6
  use lattice,  only: lattice_Sslip_v
  use mesh,     only: mesh_NcpElems, &
                      mesh_maxNips
  use material, only: homogenization_maxNgrains, &
                      material_phase, &
                      phase_constitutionInstance
  
  implicit none

  !*** input variables ***!
  real(pReal), dimension(6), intent(in) ::  Tstar_v                   ! 2nd Piola Kirchhoff stress tensor in Mandel notation
  real(pReal), intent(in) ::                Temperature
  integer(pInt), intent(in)::               g, &                      ! grain number
                                            ip, &                     ! integration point number
                                            el                        ! element number
  type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state ! state of the current microstructure
  
  !*** output variables ***!
  real(pReal), dimension(1) ::              constitutive_j2_dotState  ! evolution of state variable
  
  !*** local variables ***!
  real(pReal), dimension(6) ::              Tstar_dev_v               ! deviatoric part of the 2nd Piola Kirchhoff stress tensor in Mandel notation
  real(pReal)                               gamma_dot, &              ! strainrate
                                            hardening, &              ! hardening coefficient
                                            norm_Tstar_dev            ! euclidean norm of Tstar_dev
  integer(pInt)                             matID

  matID = phase_constitutionInstance(material_phase(g,ip,el))

  ! deviatoric part of 2nd Piola-Kirchhoff stress
  Tstar_dev_v(1:3) = Tstar_v(1:3) - sum(Tstar_v(1:3))/3.0_pReal
  Tstar_dev_v(4:6) = Tstar_v(4:6)
  
  norm_Tstar_dev = sqrt(math_mul6x6(Tstar_dev_v,Tstar_dev_v))
  
  ! gamma_dot
  gamma_dot = constitutive_j2_gdot0(matID) * (            sqrt(1.5_pReal) * norm_Tstar_dev & 
                                             / &!---------------------------------------------------
                                               (constitutive_j2_fTaylor(matID) * state(g,ip,el)%p(1)) ) ** constitutive_j2_n(matID)
  
  ! hardening coefficient
  hardening = constitutive_j2_h0(matID) * &
                  ( 1.0_pReal - state(g,ip,el)%p(1) / constitutive_j2_tausat(matID) ) ** constitutive_j2_w0(matID)
  
  ! dotState
  constitutive_j2_dotState =  hardening * gamma_dot

  return

endfunction


!****************************************************************
!* calculates the rate of change of temperature                 *
!****************************************************************
pure function constitutive_j2_dotTemperature(Tstar_v, Temperature, state, g, ip, el)

  !*** variables and functions from other modules ***!
  use prec,     only: pReal,pInt,p_vec
  use lattice,  only: lattice_Sslip_v
  use mesh,     only: mesh_NcpElems,mesh_maxNips
  use material, only: homogenization_maxNgrains,material_phase,phase_constitutionInstance  
  implicit none

  !*** input variables ***!
  real(pReal), dimension(6), intent(in) ::  Tstar_v                   ! 2nd Piola Kirchhoff stress tensor in Mandel notation
  real(pReal), intent(in) ::                Temperature
  integer(pInt), intent(in)::               g, &                      ! grain number
                                            ip, &                     ! integration point number
                                            el                        ! element number
  type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state ! state of the current microstructure
  
  !*** output variables ***!
  real(pReal) constitutive_j2_dotTemperature                          ! rate of change of temperature
  
  ! calculate dotTemperature
  constitutive_j2_dotTemperature = 0.0_pReal

  return
endfunction


!*********************************************************************
!* return array of constitutive results                              *
!*********************************************************************
pure function constitutive_j2_postResults(Tstar_v, Temperature, dt, state, g, ip, el)

!*** variables and functions from other modules ***!
  use prec,     only: pReal, &
                      pInt, &
                      p_vec
  use math,     only: math_mul6x6
  use lattice,  only: lattice_Sslip_v
  use mesh,     only: mesh_NcpElems, &
                      mesh_maxNips
  use material, only: homogenization_maxNgrains, &
                      material_phase, &
                      phase_constitutionInstance, &
                      phase_Noutput

  implicit none

  !*** input variables ***!
  real(pReal), dimension(6), intent(in)::   Tstar_v                    ! 2nd Piola Kirchhoff stress tensor in Mandel notation
  real(pReal), intent(in)::                 Temperature, &
                                            dt                         ! current time increment
  integer(pInt), intent(in)::               g, &                       ! grain number
                                            ip, &                      ! integration point number
                                            el                         ! element number
  type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state ! state of the current microstructure
  
  !*** output variables ***!
  real(pReal), dimension(constitutive_j2_sizePostResults(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            constitutive_j2_postResults
  
  !*** local variables ***!
  real(pReal), dimension(6) ::              Tstar_dev_v                ! deviatoric part of the 2nd Piola Kirchhoff stress tensor in Mandel notation
  real(pReal)                               norm_Tstar_dev             ! euclidean norm of Tstar_dev
  integer(pInt)                             matID, &
                                            o, &
                                            c

  !*** global variables ***!
  ! constitutive_j2_gdot0
  ! constitutive_j2_fTaylor
  ! constitutive_j2_n
  
  
  matID = phase_constitutionInstance(material_phase(g,ip,el))
  
  ! calculate deviatoric part of 2nd Piola-Kirchhoff stress and its norm
  Tstar_dev_v(1:3) = Tstar_v(1:3) - sum(Tstar_v(1:3))/3.0_pReal
  Tstar_dev_v(4:6) = Tstar_v(4:6)
  norm_Tstar_dev = sqrt(math_mul6x6(Tstar_dev_v,Tstar_dev_v))
  
  c = 0_pInt
  constitutive_j2_postResults = 0.0_pReal

  do o = 1,phase_Noutput(material_phase(g,ip,el))
    select case(constitutive_j2_output(o,matID))
      case ('flowstress')
        constitutive_j2_postResults(c+1) = state(g,ip,el)%p(1)
        c = c + 1
      case ('strainrate')
        constitutive_j2_postResults(c+1) = &
              constitutive_j2_gdot0(matID) * (            sqrt(1.5_pReal) * norm_Tstar_dev & 
                                             / &!---------------------------------------------------
                                               (constitutive_j2_fTaylor(matID) * state(g,ip,el)%p(1)) ) ** constitutive_j2_n(matID)
        c = c + 1
    end select
  enddo
 
 return

endfunction

END MODULE
