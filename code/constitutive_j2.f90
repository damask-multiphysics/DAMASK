! Copyright 2011-13 Max-Planck-Institut für Eisenforschung GmbH
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
!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for isotropic (J2) plasticity
!> @details Isotropic (J2) Plasticity which resembles the phenopowerlaw plasticity without
!! resolving the stress on the slip systems. Will give the response of phenopowerlaw for an
!! untextured polycrystal
!--------------------------------------------------------------------------------------------------
module constitutive_j2
 use prec, only: &
   pReal,&
   pInt
 use lattice, only: &
   LATTICE_iso_ID
 
 implicit none
 private
 integer(pInt),     dimension(:),     allocatable,         public, protected :: &
   constitutive_j2_sizeDotState, &                                                                  !< number of dotStates
   constitutive_j2_sizeState, &                                                                     !< total number of microstructural variables
   constitutive_j2_sizePostResults                                                                  !< cumulative size of post results
   
 integer(pInt),     dimension(:,:),   allocatable, target, public :: &
   constitutive_j2_sizePostResult                                                                   !< size of each post result output
   
 character(len=64), dimension(:,:),   allocatable, target, public :: &
   constitutive_j2_output                                                                           !< name of each post result output
 
 integer(kind(LATTICE_iso_ID)), dimension(:), allocatable, public :: &
   constitutive_j2_structureID                                                                !< ID of the lattice structure 
 
 integer(pInt),     dimension(:),     allocatable,         private :: &
   constitutive_j2_Noutput                                                                          !< number of outputs per instance
   
 real(pReal),       dimension(:),     allocatable,         private :: &
   constitutive_j2_fTaylor, &                                                                       !< Taylor factor
   constitutive_j2_tau0, &                                                                          !< initial plastic stress
   constitutive_j2_gdot0, &                                                                         !< reference velocity
   constitutive_j2_n, &                                                                             !< Visco-plastic parameter
!--------------------------------------------------------------------------------------------------
! h0 as function of h0 = A + B log (gammadot) 
   constitutive_j2_h0, &
   constitutive_j2_h0_slopeLnRate, &
   constitutive_j2_tausat, &                                                                        !< final plastic stress
   constitutive_j2_a, &
   constitutive_j2_aTolResistance, &
!--------------------------------------------------------------------------------------------------
! tausat += (asinh((gammadot / SinhFitA)**(1 / SinhFitD)))**(1 / SinhFitC) / (SinhFitB * (gammadot / gammadot0)**(1/n))
   constitutive_j2_tausat_SinhFitA, &                                                               !< fitting parameter for normalized strain rate vs. stress function
   constitutive_j2_tausat_SinhFitB, &                                                               !< fitting parameter for normalized strain rate vs. stress function
   constitutive_j2_tausat_SinhFitC, &                                                               !< fitting parameter for normalized strain rate vs. stress function
   constitutive_j2_tausat_SinhFitD                                                                  !< fitting parameter for normalized strain rate vs. stress function

 real(pReal),       dimension(:,:,:), allocatable,          private :: &
   constitutive_j2_Cslip_66

 public  :: &
   constitutive_j2_init, &
   constitutive_j2_stateInit, &
   constitutive_j2_aTolState, &
   constitutive_j2_homogenizedC, &
   constitutive_j2_LpAndItsTangent, &
   constitutive_j2_dotState, &
   constitutive_j2_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine constitutive_j2_init(myFile)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use math, only: &
   math_Mandel3333to66, &
   math_Voigt66to3333
 use IO, only: &
   IO_read, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue, &
   IO_error, &
   IO_timeStamp
 use material
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use lattice  

 implicit none
 integer(pInt), intent(in) :: myFile
 
 integer(pInt), parameter :: MAXNCHUNKS = 7_pInt
 
 integer(pInt), dimension(1_pInt+2_pInt*MAXNCHUNKS) :: positions
 integer(pInt) :: section = 0_pInt, maxNinstance, i,o, mySize
 character(len=65536) :: &
   tag  = '', &
   line = ''                                                                                        ! to start initialized

 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_J2_label//' init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_J2_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(constitutive_j2_sizeDotState(maxNinstance))
          constitutive_j2_sizeDotState = 0_pInt
 allocate(constitutive_j2_sizeState(maxNinstance))
          constitutive_j2_sizeState = 0_pInt
 allocate(constitutive_j2_sizePostResults(maxNinstance))
          constitutive_j2_sizePostResults = 0_pInt
 allocate(constitutive_j2_sizePostResult(maxval(phase_Noutput), maxNinstance))
          constitutive_j2_sizePostResult = 0_pInt
 allocate(constitutive_j2_output(maxval(phase_Noutput), maxNinstance))
          constitutive_j2_output = ''
 allocate(constitutive_j2_Noutput(maxNinstance))
          constitutive_j2_Noutput = 0_pInt
 allocate(constitutive_j2_structureID(maxNinstance))
          constitutive_j2_structureID        = -1
 allocate(constitutive_j2_Cslip_66(6,6,maxNinstance))
          constitutive_j2_Cslip_66 = 0.0_pReal
 allocate(constitutive_j2_fTaylor(maxNinstance))
          constitutive_j2_fTaylor = 0.0_pReal
 allocate(constitutive_j2_tau0(maxNinstance))
          constitutive_j2_tau0 = 0.0_pReal
 allocate(constitutive_j2_gdot0(maxNinstance))
          constitutive_j2_gdot0 = 0.0_pReal
 allocate(constitutive_j2_n(maxNinstance))
          constitutive_j2_n = 0.0_pReal
 allocate(constitutive_j2_h0(maxNinstance))
          constitutive_j2_h0 = 0.0_pReal
 allocate(constitutive_j2_h0_slopeLnRate(maxNinstance))
          constitutive_j2_h0_slopeLnRate = 0.0_pReal
 allocate(constitutive_j2_tausat(maxNinstance))
          constitutive_j2_tausat = 0.0_pReal
 allocate(constitutive_j2_a(maxNinstance))
          constitutive_j2_a = 0.0_pReal
 allocate(constitutive_j2_aTolResistance(maxNinstance))
          constitutive_j2_aTolResistance = 0.0_pReal
 allocate(constitutive_j2_tausat_SinhFitA(maxNinstance))
          constitutive_j2_tausat_SinhFitA = 0.0_pReal
 allocate(constitutive_j2_tausat_SinhFitB(maxNinstance))
          constitutive_j2_tausat_SinhFitB = 0.0_pReal
 allocate(constitutive_j2_tausat_SinhFitC(maxNinstance))
          constitutive_j2_tausat_SinhFitC = 0.0_pReal
 allocate(constitutive_j2_tausat_SinhFitD(maxNinstance))
          constitutive_j2_tausat_SinhFitD = 0.0_pReal
 
 rewind(myFile)
 do while (trim(line) /= '#EOF#' .and. IO_lc(IO_getTag(line,'<','>')) /= 'phase')                   ! wind forward to <phase>
   line = IO_read(myFile)
 enddo
 
 do while (trim(line) /= '#EOF#')                                                                   ! read through sections of phase part
   line = IO_read(myFile)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt                                                                     ! advance section counter
     cycle                                                                                          ! skip to next line
   endif
   if (section > 0_pInt ) then                                                                      ! do not short-circuit here (.and. with next if-statement). It's not safe in Fortran
     if (phase_plasticity(section) == PLASTICITY_J2_ID) then                                        ! one of my sections
       i = phase_plasticityInstance(section)                                                        ! which instance of my plasticity is present phase
       positions = IO_stringPos(line,MAXNCHUNKS)
       tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                           ! extract key
       select case(tag)
         case ('plasticity','elasticity')
           cycle
         case ('(output)')
           constitutive_j2_Noutput(i) = constitutive_j2_Noutput(i) + 1_pInt
           constitutive_j2_output(constitutive_j2_Noutput(i),i) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
         case ('lattice_structure')
           select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
             case(LATTICE_iso_label)
               constitutive_j2_structureID(i) = LATTICE_iso_ID
             case(LATTICE_fcc_label)
               constitutive_j2_structureID(i) = LATTICE_fcc_ID
             case(LATTICE_bcc_label)
               constitutive_j2_structureID(i) = LATTICE_bcc_ID
             case(LATTICE_hex_label)
               constitutive_j2_structureID(i) = LATTICE_hex_ID
             case(LATTICE_ort_label)
               constitutive_j2_structureID(i) = LATTICE_ort_ID
           end select
         case ('c11')
           constitutive_j2_Cslip_66(1,1,i) = IO_floatValue(line,positions,2_pInt)
         case ('c12')
           constitutive_j2_Cslip_66(1,2,i) = IO_floatValue(line,positions,2_pInt)
         case ('c13')
           constitutive_j2_Cslip_66(1,3,i) = IO_floatValue(line,positions,2_pInt)
         case ('c22')
           constitutive_j2_Cslip_66(2,2,i) = IO_floatValue(line,positions,2_pInt)
         case ('c23')
           constitutive_j2_Cslip_66(2,3,i) = IO_floatValue(line,positions,2_pInt)
         case ('c33')
           constitutive_j2_Cslip_66(3,3,i) = IO_floatValue(line,positions,2_pInt)
         case ('c44')
           constitutive_j2_Cslip_66(4,4,i) = IO_floatValue(line,positions,2_pInt)
         case ('c55')
           constitutive_j2_Cslip_66(5,5,i) = IO_floatValue(line,positions,2_pInt)
         case ('c66')
           constitutive_j2_Cslip_66(6,6,i) = IO_floatValue(line,positions,2_pInt)
         case ('tau0')
           constitutive_j2_tau0(i)         = IO_floatValue(line,positions,2_pInt)
         case ('gdot0')
           constitutive_j2_gdot0(i)        = IO_floatValue(line,positions,2_pInt)
         case ('n')
           constitutive_j2_n(i)            = IO_floatValue(line,positions,2_pInt)
         case ('h0')
           constitutive_j2_h0(i)           = IO_floatValue(line,positions,2_pInt)
         case ('h0_slope','slopelnrate')
           constitutive_j2_h0_slopeLnRate(i)  = IO_floatValue(line,positions,2_pInt)
         case ('tausat')
           constitutive_j2_tausat(i)          = IO_floatValue(line,positions,2_pInt)
         case ('tausat_sinhfita')
           constitutive_j2_tausat_SinhFitA(i) = IO_floatValue(line,positions,2_pInt)
         case ('tausat_sinhfitb')
           constitutive_j2_tausat_SinhFitB(i) = IO_floatValue(line,positions,2_pInt)
         case ('tausat_sinhfitc')
           constitutive_j2_tausat_SinhFitC(i) = IO_floatValue(line,positions,2_pInt)
         case ('tausat_sinhfitd')
           constitutive_j2_tausat_SinhFitD(i) = IO_floatValue(line,positions,2_pInt)
         case ('a', 'w0')
           constitutive_j2_a(i)               = IO_floatValue(line,positions,2_pInt)
         case ('taylorfactor')
           constitutive_j2_fTaylor(i)         = IO_floatValue(line,positions,2_pInt)
         case ('atol_resistance')
           constitutive_j2_aTolResistance(i)  = IO_floatValue(line,positions,2_pInt)
         case default
           call IO_error(210_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_J2_label//')')
       end select
     endif
   endif
 enddo

 sanityChecks: do i = 1_pInt,maxNinstance
   if (constitutive_j2_tau0(i) < 0.0_pReal)            call IO_error(211_pInt,ext_msg='tau0 (' &
                                                            //PLASTICITY_J2_label//')')
   if (constitutive_j2_gdot0(i) <= 0.0_pReal)          call IO_error(211_pInt,ext_msg='gdot0 (' &
                                                            //PLASTICITY_J2_label//')')
   if (constitutive_j2_n(i) <= 0.0_pReal)              call IO_error(211_pInt,ext_msg='n (' &
                                                            //PLASTICITY_J2_label//')')
   if (constitutive_j2_tausat(i) <= 0.0_pReal)         call IO_error(211_pInt,ext_msg='tausat (' &
                                                            //PLASTICITY_J2_label//')')
   if (constitutive_j2_a(i) <= 0.0_pReal)              call IO_error(211_pInt,ext_msg='a (' &
                                                            //PLASTICITY_J2_label//')')
   if (constitutive_j2_fTaylor(i) <= 0.0_pReal)        call IO_error(211_pInt,ext_msg='taylorfactor (' &
                                                            //PLASTICITY_J2_label//')')
   if (constitutive_j2_aTolResistance(i) <= 0.0_pReal) call IO_error(211_pInt,ext_msg='aTol_resistance (' &
                                                            //PLASTICITY_J2_label//')')
 enddo sanityChecks

 instancesLoop: do i = 1_pInt,maxNinstance
   outputsLoop: do o = 1_pInt,constitutive_j2_Noutput(i)
     select case(constitutive_j2_output(o,i))
       case('flowstress')
         mySize = 1_pInt
       case('strainrate')
         mySize = 1_pInt
       case default
         call IO_error(212_pInt,ext_msg=constitutive_j2_output(o,i)//' ('//PLASTICITY_J2_label//')')
     end select
  
     if (mySize > 0_pInt) then                                                                      ! any meaningful output found
       constitutive_j2_sizePostResult(o,i) = mySize
       constitutive_j2_sizePostResults(i) = &
       constitutive_j2_sizePostResults(i) + mySize
     endif
   enddo outputsLoop 

   constitutive_j2_sizeDotState(i) = 1_pInt
   constitutive_j2_sizeState(i)    = 1_pInt
   
   constitutive_j2_Cslip_66(1:6,1:6,i) = lattice_symmetrizeC66(constitutive_j2_structureID(i),&
                                                      constitutive_j2_Cslip_66(1:6,1:6,i)) 
   constitutive_j2_Cslip_66(1:6,1:6,i) = &
     math_Mandel3333to66(math_Voigt66to3333(constitutive_j2_Cslip_66(1:6,1:6,i)))                   ! Literature data is Voigt, DAMASK uses Mandel

 enddo instancesLoop

end subroutine constitutive_j2_init


!--------------------------------------------------------------------------------------------------
!> @brief sets the initial microstructural state for a given instance of this plasticity
!> @details initial microstructural state is set to the value specified by tau0
!--------------------------------------------------------------------------------------------------
pure function constitutive_j2_stateInit(matID)
  
 implicit none
 real(pReal),   dimension(1)            :: constitutive_j2_stateInit
 integer(pInt),              intent(in) :: matID                                               !< number specifying the instance of the plasticity
 
 constitutive_j2_stateInit = constitutive_j2_tau0(matID)

end function constitutive_j2_stateInit


!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
pure function constitutive_j2_aTolState(matID)

 implicit none
 integer(pInt), intent(in) :: matID                                                           !< number specifying the instance of the plasticity

 real(pReal), dimension(constitutive_j2_sizeState(matID)) :: &
                              constitutive_j2_aTolState

 constitutive_j2_aTolState = constitutive_j2_aTolResistance(matID)

end function constitutive_j2_aTolState


!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenized elasticity matrix
!--------------------------------------------------------------------------------------------------
pure function constitutive_j2_homogenizedC(ipc,ip,el)
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains,&
   material_phase, &
   phase_plasticityInstance
 
 implicit none
 real(pReal), dimension(6,6) :: &
   constitutive_j2_homogenizedC
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 constitutive_j2_homogenizedC = constitutive_j2_Cslip_66(1:6,1:6,&
                                               phase_plasticityInstance(material_phase(ipc,ip,el)))

end function constitutive_j2_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
pure subroutine constitutive_j2_LpAndItsTangent(Lp,dLp_dTstar99,Tstar_v,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use math, only: &
   math_mul6x6, &
   math_Mandel6to33, &
   math_Plain3333to99, &
   math_deviatoric33, &
   math_mul33xx33
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3),                                                  intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(9,9),                                                  intent(out) :: &
   dLp_dTstar99                                                                                    !< derivative of Lp with respect to 2nd Piola Kirchhoff stress

 real(pReal), dimension(6),                                                    intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),                                                                intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
   state                                                                                            !< microstructure state

 real(pReal), dimension(3,3) :: &
   Tstar_dev_33                                                                                     !< deviatoric part of the 2nd Piola Kirchhoff stress tensor as 2nd order tensor
 real(pReal), dimension(3,3,3,3) :: &
   dLp_dTstar_3333                                                                                  !< derivative of Lp with respect to Tstar as 4th order tensor
 real(pReal) :: &
   gamma_dot, &                                                                                     !< strainrate
   norm_Tstar_dev, &                                                                                !< euclidean norm of Tstar_dev
   squarenorm_Tstar_dev                                                                             !< square of the euclidean norm of Tstar_dev
 integer(pInt) :: &
   matID, &
   k, l, m, n

 matID = phase_plasticityInstance(material_phase(ipc,ip,el))
 Tstar_dev_33 = math_deviatoric33(math_Mandel6to33(Tstar_v))                                        ! deviatoric part of 2nd Piola-Kirchhoff stress
 squarenorm_Tstar_dev = math_mul33xx33(Tstar_dev_33,Tstar_dev_33)
 norm_Tstar_dev = sqrt(squarenorm_Tstar_dev) 

 if (norm_Tstar_dev <= 0.0_pReal) then                                                              ! Tstar == 0 --> both Lp and dLp_dTstar are zero
   Lp = 0.0_pReal
   dLp_dTstar99 = 0.0_pReal
 else
   gamma_dot = constitutive_j2_gdot0(matID) * (            sqrt(1.5_pReal) * norm_Tstar_dev & 
             / &!----------------------------------------------------------------------------------
               (constitutive_j2_fTaylor(matID) * state(ipc,ip,el)%p(1)) ) **constitutive_j2_n(matID)

   Lp = Tstar_dev_33/norm_Tstar_dev * gamma_dot/constitutive_j2_fTaylor(matID)

!--------------------------------------------------------------------------------------------------
! Calculation of the tangent of Lp
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,l,m,n) = (constitutive_j2_n(matID)-1.0_pReal) * &
                                      Tstar_dev_33(k,l)*Tstar_dev_33(m,n) / squarenorm_Tstar_dev
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,l,k,l) = dLp_dTstar_3333(k,l,k,l) + 1.0_pReal
   dLp_dTstar99 = math_Plain3333to99(gamma_dot / constitutive_j2_fTaylor(matID) * &
                                      dLp_dTstar_3333 / norm_Tstar_dev)
 end if

end subroutine constitutive_j2_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
pure function constitutive_j2_dotState(Tstar_v,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use math, only: &
   math_mul6x6
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance
 
 implicit none
 real(pReal), dimension(1) :: &
   constitutive_j2_dotState
 real(pReal), dimension(6),                                                    intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),                                                                intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
   state                                                                                            !< microstructure state

 real(pReal), dimension(6) :: &
   Tstar_dev_v                                                                                      !< deviatoric part of the 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal) :: &
   gamma_dot, &                                                                                     !< strainrate
   hardening, &                                                                                     !< hardening coefficient
   saturation, &                                                                                    !< saturation resistance
   norm_Tstar_dev                                                                                   !< euclidean norm of Tstar_dev
 integer(pInt) :: &
   matID

 matID = phase_plasticityInstance(material_phase(ipc,ip,el))
!--------------------------------------------------------------------------------------------------
! norm of deviatoric part of 2nd Piola-Kirchhoff stress
 Tstar_dev_v(1:3) = Tstar_v(1:3) - sum(Tstar_v(1:3))/3.0_pReal
 Tstar_dev_v(4:6) = Tstar_v(4:6)
 norm_Tstar_dev = sqrt(math_mul6x6(Tstar_dev_v,Tstar_dev_v))

!--------------------------------------------------------------------------------------------------
! strain rate 
 gamma_dot = constitutive_j2_gdot0(matID) * (            sqrt(1.5_pReal) * norm_Tstar_dev & 
            / &!-----------------------------------------------------------------------------------
              (constitutive_j2_fTaylor(matID) * state(ipc,ip,el)%p(1)) ) ** constitutive_j2_n(matID)
 
!--------------------------------------------------------------------------------------------------
! hardening coefficient
 if (abs(gamma_dot) > 1e-12_pReal) then
   if (constitutive_j2_tausat_SinhFitA(matID) == 0.0_pReal) then
     saturation = constitutive_j2_tausat(matID)
   else
     saturation = (  constitutive_j2_tausat(matID) &
                   + ( log(  ( gamma_dot / constitutive_j2_tausat_SinhFitA(matID)&
                               )**(1.0_pReal / constitutive_j2_tausat_SinhFitD(matID))&
                            + sqrt(  ( gamma_dot / constitutive_j2_tausat_SinhFitA(matID) &
                                      )**(2.0_pReal / constitutive_j2_tausat_SinhFitD(matID)) &
                                   + 1.0_pReal ) &
                            ) & ! asinh(K) = ln(K + sqrt(K^2 +1))
                       )**(1.0_pReal / constitutive_j2_tausat_SinhFitC(matID)) &
                   / (  constitutive_j2_tausat_SinhFitB(matID) &
                      * (gamma_dot / constitutive_j2_gdot0(matID))**(1.0_pReal / constitutive_j2_n(matID)) &
                      ) &
                   )
   endif
   hardening = ( constitutive_j2_h0(matID) + constitutive_j2_h0_slopeLnRate(matID) * log(gamma_dot) ) &
               * abs( 1.0_pReal - state(ipc,ip,el)%p(1)/saturation )**constitutive_j2_a(matID) &
               * sign(1.0_pReal, 1.0_pReal - state(ipc,ip,el)%p(1)/saturation)
 else
   hardening = 0.0_pReal
 endif
 
 constitutive_j2_dotState = hardening * gamma_dot

end function constitutive_j2_dotState


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
pure function constitutive_j2_postResults(Tstar_v,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use math, only: &
   math_mul6x6
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance, &
   phase_Noutput

 implicit none
 real(pReal), dimension(6),                                                    intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),                                                                intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
   state                                                                                            !< microstructure state

 real(pReal), dimension(constitutive_j2_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
                                           constitutive_j2_postResults
 
 real(pReal), dimension(6) :: &
   Tstar_dev_v                                                                                      ! deviatoric part of the 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal) :: &
   norm_Tstar_dev                                                                                   ! euclidean norm of Tstar_dev
 integer(pInt) :: &
   matID, &
   o, &
   c
 
 matID = phase_plasticityInstance(material_phase(ipc,ip,el))
 
!--------------------------------------------------------------------------------------------------
! calculate deviatoric part of 2nd Piola-Kirchhoff stress and its norm
 Tstar_dev_v(1:3) = Tstar_v(1:3) - sum(Tstar_v(1:3))/3.0_pReal
 Tstar_dev_v(4:6) = Tstar_v(4:6)
 norm_Tstar_dev = sqrt(math_mul6x6(Tstar_dev_v,Tstar_dev_v))
 
 c = 0_pInt
 constitutive_j2_postResults = 0.0_pReal

 outputsLoop: do o = 1_pInt,phase_Noutput(material_phase(ipc,ip,el))
   select case(constitutive_j2_output(o,matID))
     case ('flowstress')
       constitutive_j2_postResults(c+1_pInt) = state(ipc,ip,el)%p(1)
       c = c + 1_pInt
     case ('strainrate')
       constitutive_j2_postResults(c+1_pInt) = &
                constitutive_j2_gdot0(matID) * (            sqrt(1.5_pReal) * norm_Tstar_dev & 
             / &!----------------------------------------------------------------------------------
              (constitutive_j2_fTaylor(matID) * state(ipc,ip,el)%p(1)) ) ** constitutive_j2_n(matID)
       c = c + 1_pInt
   end select
 enddo outputsLoop

end function constitutive_j2_postResults

end module constitutive_j2
