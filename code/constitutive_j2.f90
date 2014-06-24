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
#ifdef HDF
 use hdf5, only: &
   HID_T
#endif
 use prec, only: &
   pReal,&
   pInt
 
 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,         public, protected :: &
#ifndef NEWSTATE
   constitutive_j2_sizeDotState, &                                                                  !< number of dotStates
   constitutive_j2_sizeState, &                                                                     !< total number of microstructural variables
#endif
   constitutive_j2_sizePostResults                                                                  !< cumulative size of post results
   
 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   constitutive_j2_sizePostResult                                                                   !< size of each post result output
   
 character(len=64),                   dimension(:,:),   allocatable, target, public :: &
   constitutive_j2_output                                                                           !< name of each post result output
 
 integer(pInt),                       dimension(:),     allocatable,         private :: &
   constitutive_j2_Noutput                                                                          !< number of outputs per instance
 real(pReal),                         dimension(:),     allocatable,         private :: &
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

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 flowstress_ID, &
                 strainrate_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),   allocatable,          private :: & 
   constitutive_j2_outputID                                                                         !< ID of each post result output
 
#ifdef HDF 
 type constitutive_j2_tOutput
   real(pReal),                         dimension(:),     allocatable,         private :: &
     flowstress, &
     strainrate
   logical :: flowstressActive = .false., strainrateActive = .false.                ! if we can write the output block wise, this is not needed anymore because we can do an if(allocated(xxx))                                 
 end type constitutive_j2_tOutput
 type(constitutive_j2_tOutput), allocatable, dimension(:) :: constitutive_j2_Output2
integer(HID_T), allocatable, dimension(:) :: outID
#endif 
 public  :: &
   constitutive_j2_init, &
#ifndef NEWSTATE 
   constitutive_j2_stateInit, &
   constitutive_j2_aTolState, &
#endif
   constitutive_j2_LpAndItsTangent, &
   constitutive_j2_dotState, &
   constitutive_j2_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine constitutive_j2_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
#ifdef HDF 
 use hdf5
#endif 
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use numerics, only: &
   numerics_integrator
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
   IO_timeStamp, &
#ifdef HDF 
   tempResults, &
   HDF5_addGroup, &
   HDF5_addScalarDataset,&
#endif
   IO_EOF
 use material, only: &
   homogenization_maxNgrains, &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_J2_label, &
   PLASTICITY_J2_ID, &
   material_phase, &
#ifdef NEWSTATE
   plasticState, &
#endif
   MATERIAL_partPhase
   
 use lattice  

 implicit none
 integer(pInt), intent(in) :: fileUnit
 
 integer(pInt), parameter :: MAXNCHUNKS = 7_pInt
 
 integer(pInt), dimension(1_pInt+2_pInt*MAXNCHUNKS) :: positions
 integer(pInt) :: phase, maxNinstance, instance,o, mySize, myConstituents
 character(len=65536) :: &
   tag  = '', &
   line = ''
  integer(pInt) :: NofMyPhase
#ifdef HDF 
 character(len=5) :: &
   str1
 integer(HID_T) :: ID,ID2,ID4
#endif
#ifdef NEWSTATE
 integer(pInt) :: sizeDotState,sizeState
#endif
 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_J2_label//' init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_J2_ID),pInt)
 if (maxNinstance == 0_pInt) return

#ifdef HDF 
  allocate(constitutive_j2_Output2(maxNinstance))
  allocate(outID(maxNinstance))
#endif

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
#ifndef NEWSTATE 
 allocate(constitutive_j2_sizeDotState(maxNinstance),                         source=1_pInt)
 allocate(constitutive_j2_sizeState(maxNinstance),                            source=1_pInt)
#endif
 allocate(constitutive_j2_sizePostResults(maxNinstance),                      source=0_pInt)
 allocate(constitutive_j2_sizePostResult(maxval(phase_Noutput), maxNinstance),source=0_pInt)
 allocate(constitutive_j2_output(maxval(phase_Noutput), maxNinstance))
          constitutive_j2_output = ''
 allocate(constitutive_j2_outputID(maxval(phase_Noutput),maxNinstance),       source=undefined_ID)
 allocate(constitutive_j2_Noutput(maxNinstance),                              source=0_pInt)
 allocate(constitutive_j2_fTaylor(maxNinstance),                              source=0.0_pReal)
 allocate(constitutive_j2_tau0(maxNinstance),                                 source=0.0_pReal)
 allocate(constitutive_j2_gdot0(maxNinstance),                                source=0.0_pReal)
 allocate(constitutive_j2_n(maxNinstance),                                    source=0.0_pReal)
 allocate(constitutive_j2_h0(maxNinstance),                                   source=0.0_pReal)
 allocate(constitutive_j2_h0_slopeLnRate(maxNinstance),                       source=0.0_pReal)
 allocate(constitutive_j2_tausat(maxNinstance),                               source=0.0_pReal)
 allocate(constitutive_j2_a(maxNinstance),                                    source=0.0_pReal)
 allocate(constitutive_j2_aTolResistance(maxNinstance),                       source=0.0_pReal)
 allocate(constitutive_j2_tausat_SinhFitA(maxNinstance),                      source=0.0_pReal)
 allocate(constitutive_j2_tausat_SinhFitB(maxNinstance),                      source=0.0_pReal)
 allocate(constitutive_j2_tausat_SinhFitC(maxNinstance),                      source=0.0_pReal)
 allocate(constitutive_j2_tausat_SinhFitD(maxNinstance),                      source=0.0_pReal)
 
 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo
 
 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     myConstituents = 0_pInt
     phase = phase + 1_pInt                                                                         ! advance section counter
     if (phase_plasticity(phase) == PLASTICITY_J2_ID) then
       instance = phase_plasticityInstance(phase)
       myConstituents = count(material_phase==phase)
#ifdef HDF
       outID(instance)=HDF5_addGroup(str1,tempResults)
#endif
     endif
     cycle                                                                                          ! skip to next line
   endif
   if (myConstituents > 0_pInt ) then
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my plasticity is present phase
     positions = IO_stringPos(line,MAXNCHUNKS) 
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key

     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('flowstress')
             constitutive_j2_outputID(constitutive_j2_Noutput(instance),instance) = flowstress_ID
             constitutive_j2_Noutput(instance) = constitutive_j2_Noutput(instance) + 1_pInt
             constitutive_j2_output(constitutive_j2_Noutput(instance),instance) = &
                                                IO_lc(IO_stringValue(line,positions,2_pInt))
#ifdef HDF 
             call HDF5_addScalarDataset(outID(instance),myConstituents,'flowstress','MPa')
             allocate(constitutive_j2_Output2(instance)%flowstress(myConstituents))
             constitutive_j2_Output2(instance)%flowstressActive = .true.
#endif
           case ('strainrate')
             constitutive_j2_outputID(constitutive_j2_Noutput(instance),instance) = strainrate_ID
             constitutive_j2_Noutput(instance) = constitutive_j2_Noutput(instance) + 1_pInt
             constitutive_j2_output(constitutive_j2_Noutput(instance),instance) = &
                                                IO_lc(IO_stringValue(line,positions,2_pInt))
#ifdef HDF 
             call HDF5_addScalarDataset(outID(instance),myConstituents,'strainrate','1/s')
             allocate(constitutive_j2_Output2(instance)%strainrate(myConstituents))
             constitutive_j2_Output2(instance)%strainrateActive = .true.
#endif
           case default

         end select
       case ('tau0')
         constitutive_j2_tau0(instance)         = IO_floatValue(line,positions,2_pInt)
         if (constitutive_j2_tau0(instance) < 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_J2_label//')')
       case ('gdot0')
         constitutive_j2_gdot0(instance)        = IO_floatValue(line,positions,2_pInt)
         if (constitutive_j2_gdot0(instance) <= 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_J2_label//')')
       case ('n')
         constitutive_j2_n(instance)            = IO_floatValue(line,positions,2_pInt)
         if (constitutive_j2_n(instance) <= 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_J2_label//')')
       case ('h0')
         constitutive_j2_h0(instance)           = IO_floatValue(line,positions,2_pInt)
       case ('h0_slope','slopelnrate')
         constitutive_j2_h0_slopeLnRate(instance)  = IO_floatValue(line,positions,2_pInt)
       case ('tausat')
         constitutive_j2_tausat(instance)          = IO_floatValue(line,positions,2_pInt)
         if (constitutive_j2_tausat(instance) <= 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_J2_label//')')
       case ('tausat_sinhfita')
         constitutive_j2_tausat_SinhFitA(instance) = IO_floatValue(line,positions,2_pInt)
       case ('tausat_sinhfitb')
         constitutive_j2_tausat_SinhFitB(instance) = IO_floatValue(line,positions,2_pInt)
       case ('tausat_sinhfitc')
         constitutive_j2_tausat_SinhFitC(instance) = IO_floatValue(line,positions,2_pInt)
       case ('tausat_sinhfitd')
         constitutive_j2_tausat_SinhFitD(instance) = IO_floatValue(line,positions,2_pInt)
       case ('a', 'w0')
         constitutive_j2_a(instance)               = IO_floatValue(line,positions,2_pInt)
         if (constitutive_j2_a(instance) <= 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_J2_label//')')
       case ('taylorfactor')
         constitutive_j2_fTaylor(instance)         = IO_floatValue(line,positions,2_pInt)
         if (constitutive_j2_fTaylor(instance) <= 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_J2_label//')')
       case ('atol_resistance')
         constitutive_j2_aTolResistance(instance)  = IO_floatValue(line,positions,2_pInt)
         if (constitutive_j2_aTolResistance(instance) <= 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_J2_label//')')
       case default

     end select
   endif
 enddo parsingFile

 initializeInstances: do phase = 1_pInt, size(phase_plasticity)
   NofMyPhase=count(material_phase==phase)
   if (phase_plasticity(phase) == PLASTICITY_j2_ID .and. NofMyPhase/=0) then
     instance = phase_plasticityInstance(phase)
     outputsLoop: do o = 1_pInt,constitutive_j2_Noutput(instance)
       select case(constitutive_j2_outputID(o,instance))
         case(flowstress_ID,strainrate_ID)
           mySize = 1_pInt
         case default
       end select
  
       if (mySize > 0_pInt) then                                                                      ! any meaningful output found
         constitutive_j2_sizePostResult(o,instance) = mySize
         constitutive_j2_sizePostResults(instance) = &
         constitutive_j2_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop
#ifdef NEWSTATE
     sizeState    = 1 
     plasticState(phase)%sizeState = sizeState
     sizeDotState = sizeState
     plasticState(phase)%sizeDotState = sizeDotState
     allocate(plasticState(phase)%state0         (sizeState,NofMyPhase),source=constitutive_j2_tau0(instance))
     allocate(plasticState(phase)%partionedState0(sizeState,NofMyPhase),source=constitutive_j2_tau0(instance))
     allocate(plasticState(phase)%subState0      (sizeState,NofMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%state          (sizeState,NofMyPhase),source=constitutive_j2_tau0(instance))
     allocate(plasticState(phase)%state_backup   (sizeState,NofMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%aTolState      (NofMyPhase),source=constitutive_j2_aTolResistance(instance))
     allocate(plasticState(phase)%dotState            (sizeDotState,NofMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%dotState_backup     (sizeDotState,NofMyPhase),source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(phase)%previousDotState  (sizeDotState,NofMyPhase),source=0.0_pReal)
       allocate(plasticState(phase)%previousDotState2 (sizeDotState,NofMyPhase),source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(phase)%RK4dotState       (sizeDotState,NofMyPhase),source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(phase)%RKCK45dotState    (6,sizeDotState,NofMyPhase),source=0.0_pReal)
#endif
   endif
 enddo initializeInstances

end subroutine constitutive_j2_init

#ifndef NEWSTATE
!--------------------------------------------------------------------------------------------------
!> @brief sets the initial microstructural state for a given instance of this plasticity
!> @details initial microstructural state is set to the value specified by tau0
! not needed for new state
!--------------------------------------------------------------------------------------------------
pure function constitutive_j2_stateInit(instance)
  
 implicit none
 real(pReal),   dimension(1)            :: constitutive_j2_stateInit
 integer(pInt),              intent(in) :: instance                                                 !< number specifying the instance of the plasticity
 
 constitutive_j2_stateInit = constitutive_j2_tau0(instance)


end function constitutive_j2_stateInit


!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
! not needed for new state
!--------------------------------------------------------------------------------------------------
pure function constitutive_j2_aTolState(instance)

 implicit none
 real(pReal),   dimension(1)            :: constitutive_j2_aTolState
 integer(pInt),              intent(in) :: instance                                                 !< number specifying the instance of the plasticity


 constitutive_j2_aTolState = constitutive_j2_aTolResistance(instance)

end function constitutive_j2_aTolState

#endif

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
 real(pReal), dimension(3,3), intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(9,9), intent(out) :: &
   dLp_dTstar99                                                                                     !< derivative of Lp with respect to 2nd Piola Kirchhoff stress

 real(pReal), dimension(6),   intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),               intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
#ifdef NEWSTATE
 real(pReal), dimension(1), intent(in) :: &
   state
#else
 type(p_vec),               intent(in) :: &
   state                                                                                            !< microstructure state
#endif 

 real(pReal), dimension(3,3) :: &
   Tstar_dev_33                                                                                     !< deviatoric part of the 2nd Piola Kirchhoff stress tensor as 2nd order tensor
 real(pReal), dimension(3,3,3,3) :: &
   dLp_dTstar_3333                                                                                  !< derivative of Lp with respect to Tstar as 4th order tensor
 real(pReal) :: &
   gamma_dot, &                                                                                     !< strainrate
   norm_Tstar_dev, &                                                                                !< euclidean norm of Tstar_dev
   squarenorm_Tstar_dev                                                                             !< square of the euclidean norm of Tstar_dev
 integer(pInt) :: &
   instance, &
   k, l, m, n

 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 Tstar_dev_33 = math_deviatoric33(math_Mandel6to33(Tstar_v))                                        ! deviatoric part of 2nd Piola-Kirchhoff stress
 squarenorm_Tstar_dev = math_mul33xx33(Tstar_dev_33,Tstar_dev_33)
 norm_Tstar_dev = sqrt(squarenorm_Tstar_dev) 

 if (norm_Tstar_dev <= 0.0_pReal) then                                                              ! Tstar == 0 --> both Lp and dLp_dTstar are zero
   Lp = 0.0_pReal
   dLp_dTstar99 = 0.0_pReal
 else
#ifdef NEWSTATE
   gamma_dot = constitutive_j2_gdot0(instance) &
             * (sqrt(1.5_pReal) * norm_Tstar_dev / constitutive_j2_fTaylor(instance) / state(1)) &
                                                  **constitutive_j2_n(instance)
#else
   gamma_dot = constitutive_j2_gdot0(instance) &
             * (sqrt(1.5_pReal) * norm_Tstar_dev / constitutive_j2_fTaylor(instance) / state%p(1)) &
                                                  **constitutive_j2_n(instance)
#endif
   Lp = Tstar_dev_33/norm_Tstar_dev * gamma_dot/constitutive_j2_fTaylor(instance)

!--------------------------------------------------------------------------------------------------
! Calculation of the tangent of Lp
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,l,m,n) = (constitutive_j2_n(instance)-1.0_pReal) * &
                                      Tstar_dev_33(k,l)*Tstar_dev_33(m,n) / squarenorm_Tstar_dev
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,l,k,l) = dLp_dTstar_3333(k,l,k,l) + 1.0_pReal
   forall (k=1_pInt:3_pInt,m=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,k,m,m) = dLp_dTstar_3333(k,k,m,m) - 1.0_pReal/3.0_pReal
   dLp_dTstar99 = math_Plain3333to99(gamma_dot / constitutive_j2_fTaylor(instance) * &
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
real(pReal) :: &
   tempState
 real(pReal), dimension(6), intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),             intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
#ifdef NEWSTATE
 real(pReal), dimension(1), intent(in) :: &
   state
#else
 type(p_vec),               intent(in) :: &
   state                                                                                            !< microstructure state
#endif
 real(pReal), dimension(6) :: &
   Tstar_dev_v                                                                                      !< deviatoric part of the 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal) :: &
   gamma_dot, &                                                                                     !< strainrate
   hardening, &                                                                                     !< hardening coefficient
   saturation, &                                                                                    !< saturation resistance
   norm_Tstar_dev                                                                                   !< euclidean norm of Tstar_dev
 integer(pInt) :: &
   instance
#ifdef NEWSTATE
  tempState = state(1)
#else
  tempState = state%p(1)
#endif

 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
!--------------------------------------------------------------------------------------------------
! norm of deviatoric part of 2nd Piola-Kirchhoff stress
 Tstar_dev_v(1:3) = Tstar_v(1:3) - sum(Tstar_v(1:3))/3.0_pReal
 Tstar_dev_v(4:6) = Tstar_v(4:6)
 norm_Tstar_dev = sqrt(math_mul6x6(Tstar_dev_v,Tstar_dev_v))

!--------------------------------------------------------------------------------------------------
! strain rate 
 gamma_dot = constitutive_j2_gdot0(instance) * (            sqrt(1.5_pReal) * norm_Tstar_dev & 
            / &!-----------------------------------------------------------------------------------
              (constitutive_j2_fTaylor(instance) * tempState) ) ** constitutive_j2_n(instance)
 
!--------------------------------------------------------------------------------------------------
! hardening coefficient
 if (abs(gamma_dot) > 1e-12_pReal) then
   if (constitutive_j2_tausat_SinhFitA(instance) == 0.0_pReal) then
     saturation = constitutive_j2_tausat(instance)
   else
     saturation = (  constitutive_j2_tausat(instance) &
                   + ( log(  ( gamma_dot / constitutive_j2_tausat_SinhFitA(instance)&
                               )**(1.0_pReal / constitutive_j2_tausat_SinhFitD(instance))&
                            + sqrt(  ( gamma_dot / constitutive_j2_tausat_SinhFitA(instance) &
                                      )**(2.0_pReal / constitutive_j2_tausat_SinhFitD(instance)) &
                                   + 1.0_pReal ) &
                            ) & ! asinh(K) = ln(K + sqrt(K^2 +1))
                       )**(1.0_pReal / constitutive_j2_tausat_SinhFitC(instance)) &
                   / (  constitutive_j2_tausat_SinhFitB(instance) &
                      * (gamma_dot / constitutive_j2_gdot0(instance))**(1.0_pReal / constitutive_j2_n(instance)) &
                      ) &
                   )
   endif
   hardening = ( constitutive_j2_h0(instance) + constitutive_j2_h0_slopeLnRate(instance) * log(gamma_dot) ) &
               * abs( 1.0_pReal - tempState/saturation )**constitutive_j2_a(instance) &
               * sign(1.0_pReal, 1.0_pReal - tempState/saturation)
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
 real(pReal), dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal) :: &
  tempState
#ifdef NEWSTATE
 real(pReal), dimension(1), intent(in) :: &
   state
#else
 type(p_vec),               intent(in) :: &
   state                                                                                            !< microstructure state
#endif
 real(pReal), dimension(constitutive_j2_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
                                           constitutive_j2_postResults

 real(pReal), dimension(6) :: &
   Tstar_dev_v                                                                                      ! deviatoric part of the 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal) :: &
   norm_Tstar_dev                                                                                   ! euclidean norm of Tstar_dev
 integer(pInt) :: &
   instance, &
   o, &
   c

#ifdef NEWSTATE
  tempState = state(1)
#else
  tempState = state%p(1)
#endif
 
 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 
!--------------------------------------------------------------------------------------------------
! calculate deviatoric part of 2nd Piola-Kirchhoff stress and its norm
 Tstar_dev_v(1:3) = Tstar_v(1:3) - sum(Tstar_v(1:3))/3.0_pReal
 Tstar_dev_v(4:6) = Tstar_v(4:6)
 norm_Tstar_dev = sqrt(math_mul6x6(Tstar_dev_v,Tstar_dev_v))
 
 c = 0_pInt
 constitutive_j2_postResults = 0.0_pReal

 outputsLoop: do o = 1_pInt,constitutive_j2_Noutput(instance)
   select case(constitutive_j2_outputID(o,instance))
     case (flowstress_ID)
       constitutive_j2_postResults(c+1_pInt) = tempState
       c = c + 1_pInt
     case (strainrate_ID)
       constitutive_j2_postResults(c+1_pInt) = &
                constitutive_j2_gdot0(instance) * (            sqrt(1.5_pReal) * norm_Tstar_dev & 
             / &!----------------------------------------------------------------------------------
              (constitutive_j2_fTaylor(instance) * tempState) ) ** constitutive_j2_n(instance)
       c = c + 1_pInt
   end select
 enddo outputsLoop

end function constitutive_j2_postResults


end module constitutive_j2
