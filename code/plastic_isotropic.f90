!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for isotropic (ISOTROPIC) plasticity
!> @details Isotropic (ISOTROPIC) Plasticity which resembles the phenopowerlaw plasticity without
!! resolving the stress on the slip systems. Will give the response of phenopowerlaw for an
!! untextured polycrystal
!--------------------------------------------------------------------------------------------------
module plastic_isotropic
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
   plastic_isotropic_sizePostResults                                                                  !< cumulative size of post results
   
 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   plastic_isotropic_sizePostResult                                                                   !< size of each post result output
   
 character(len=64),                   dimension(:,:),   allocatable, target, public :: &
   plastic_isotropic_output                                                                           !< name of each post result output
 
 integer(pInt),                       dimension(:),     allocatable, target, public :: &
   plastic_isotropic_Noutput                                                                          !< number of outputs per instance

 logical,                             dimension(:),     allocatable,         private :: &
   plastic_isotropic_dilatation                                                                       !< flag to indicate dilatation contribution of plasticity

 real(pReal),                         dimension(:),     allocatable,         private :: &
   plastic_isotropic_fTaylor, &                                                                       !< Taylor factor
   plastic_isotropic_tau0, &                                                                          !< initial plastic stress
   plastic_isotropic_gdot0, &                                                                         !< reference velocity
   plastic_isotropic_n, &                                                                             !< Visco-plastic parameter
!--------------------------------------------------------------------------------------------------
! h0 as function of h0 = A + B log (gammadot) 
   plastic_isotropic_h0, &
   plastic_isotropic_h0_slopeLnRate, &
   plastic_isotropic_tausat, &                                                                        !< final plastic stress
   plastic_isotropic_a, &
   plastic_isotropic_aTolResistance, &
   plastic_isotropic_aTolShear, &
!--------------------------------------------------------------------------------------------------
! tausat += (asinh((gammadot / SinhFitA)**(1 / SinhFitD)))**(1 / SinhFitC) / (SinhFitB * (gammadot / gammadot0)**(1/n))
   plastic_isotropic_tausat_SinhFitA, &                                                               !< fitting parameter for normalized strain rate vs. stress function
   plastic_isotropic_tausat_SinhFitB, &                                                               !< fitting parameter for normalized strain rate vs. stress function
   plastic_isotropic_tausat_SinhFitC, &                                                               !< fitting parameter for normalized strain rate vs. stress function
   plastic_isotropic_tausat_SinhFitD                                                                  !< fitting parameter for normalized strain rate vs. stress function

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 flowstress_ID, &
                 strainrate_ID
 end enum
 integer(kind(undefined_ID)),           dimension(:,:),   allocatable,        private :: & 
   plastic_isotropic_outputID                                                                         !< ID of each post result output
 

#ifdef HDF 
 type plastic_isotropic_tOutput
   real(pReal),                         dimension(:),     allocatable,         private :: &
     flowstress, &
     strainrate
 logical :: flowstressActive = .false., strainrateActive = .false.                ! if we can write the output block wise, this is not needed anymore because we can do an if(allocated(xxx))                                 
 end type plastic_isotropic_tOutput
 type(plastic_isotropic_tOutput), allocatable, dimension(:) :: plastic_isotropic_Output2
integer(HID_T), allocatable, dimension(:) :: outID
#endif


 public  :: &
   plastic_isotropic_init, &
   plastic_isotropic_LpAndItsTangent, &
   plastic_isotropic_LiAndItsTangent, &
   plastic_isotropic_dotState, &
   plastic_isotropic_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
#ifdef HDF 
 use hdf5
#endif 
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use numerics, only: &
   analyticJaco, &
   worldrank, &
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
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_ISOTROPIC_label, &
   PLASTICITY_ISOTROPIC_ID, &
   material_phase, &
   plasticState, &
   MATERIAL_partPhase
   
 use lattice  

 implicit none
 integer(pInt), intent(in) :: fileUnit
 
 
 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: &
   o, &
   phase, & 
   maxNinstance, &
   instance, &
   mySize, &
   sizeDotState, &
   sizeState, &
   sizeDeltaState
 character(len=65536) :: &
   tag  = '', &
   line = ''
  integer(pInt) :: NofMyPhase

#ifdef HDF 
 character(len=5) :: &
   str1
 integer(HID_T) :: ID,ID2,ID4
#endif

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_ISOTROPIC_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_ISOTROPIC_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

#ifdef HDF 
  allocate(plastic_isotropic_Output2(maxNinstance))
  allocate(outID(maxNinstance))
#endif

 allocate(plastic_isotropic_sizePostResults(maxNinstance),                      source=0_pInt)
 allocate(plastic_isotropic_sizePostResult(maxval(phase_Noutput), maxNinstance),source=0_pInt)
 allocate(plastic_isotropic_output(maxval(phase_Noutput), maxNinstance))
          plastic_isotropic_output = ''
 allocate(plastic_isotropic_outputID(maxval(phase_Noutput),maxNinstance),       source=undefined_ID)
 allocate(plastic_isotropic_Noutput(maxNinstance),                              source=0_pInt)
 allocate(plastic_isotropic_fTaylor(maxNinstance),                              source=0.0_pReal)
 allocate(plastic_isotropic_tau0(maxNinstance),                                 source=0.0_pReal)
 allocate(plastic_isotropic_gdot0(maxNinstance),                                source=0.0_pReal)
 allocate(plastic_isotropic_n(maxNinstance),                                    source=0.0_pReal)
 allocate(plastic_isotropic_h0(maxNinstance),                                   source=0.0_pReal)
 allocate(plastic_isotropic_h0_slopeLnRate(maxNinstance),                       source=0.0_pReal)
 allocate(plastic_isotropic_tausat(maxNinstance),                               source=0.0_pReal)
 allocate(plastic_isotropic_a(maxNinstance),                                    source=0.0_pReal)
 allocate(plastic_isotropic_aTolResistance(maxNinstance),                       source=0.0_pReal)
 allocate(plastic_isotropic_aTolShear     (maxNinstance),                       source=0.0_pReal)
 allocate(plastic_isotropic_tausat_SinhFitA(maxNinstance),                      source=0.0_pReal)
 allocate(plastic_isotropic_tausat_SinhFitB(maxNinstance),                      source=0.0_pReal)
 allocate(plastic_isotropic_tausat_SinhFitC(maxNinstance),                      source=0.0_pReal)
 allocate(plastic_isotropic_tausat_SinhFitD(maxNinstance),                      source=0.0_pReal)
 allocate(plastic_isotropic_dilatation(maxNinstance),                           source=.false.)
 
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
     phase = phase + 1_pInt                                                                         ! advance section counter
     if (phase_plasticity(phase) == PLASTICITY_ISOTROPIC_ID) then
       instance = phase_plasticityInstance(phase)
#ifdef HDF
       outID(instance)=HDF5_addGroup(str1,tempResults)
#endif
     endif
     cycle                                                                                          ! skip to next line
   endif
   if (phase > 0_pInt ) then; if (phase_plasticity(phase) == PLASTICITY_ISOTROPIC_ID) then                 ! one of my phases. Do not short-circuit here (.and. between if-statements), it's not safe in Fortran
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my plasticity is present phase
     chunkPos = IO_stringPos(line) 
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key

     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case ('flowstress')
             plastic_isotropic_Noutput(instance) = plastic_isotropic_Noutput(instance) + 1_pInt
             plastic_isotropic_outputID(plastic_isotropic_Noutput(instance),instance) = flowstress_ID
             plastic_isotropic_output(plastic_isotropic_Noutput(instance),instance) = &
                                                IO_lc(IO_stringValue(line,chunkPos,2_pInt))
#ifdef HDF 
             call HDF5_addScalarDataset(outID(instance),myConstituents,'flowstress','MPa')
             allocate(plastic_isotropic_Output2(instance)%flowstress(myConstituents))
             plastic_isotropic_Output2(instance)%flowstressActive = .true.
#endif
           case ('strainrate')
             plastic_isotropic_Noutput(instance) = plastic_isotropic_Noutput(instance) + 1_pInt
             plastic_isotropic_outputID(plastic_isotropic_Noutput(instance),instance) = strainrate_ID
             plastic_isotropic_output(plastic_isotropic_Noutput(instance),instance) = &
                                                IO_lc(IO_stringValue(line,chunkPos,2_pInt))
#ifdef HDF 
             call HDF5_addScalarDataset(outID(instance),myConstituents,'strainrate','1/s')
             allocate(plastic_isotropic_Output2(instance)%strainrate(myConstituents))
             plastic_isotropic_Output2(instance)%strainrateActive = .true.
#endif
           case default

         end select
       case ('/dilatation/')
         plastic_isotropic_dilatation(instance)   = .true.
       case ('tau0')
         plastic_isotropic_tau0(instance)         = IO_floatValue(line,chunkPos,2_pInt)
         if (plastic_isotropic_tau0(instance) < 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_ISOTROPIC_label//')')
       case ('gdot0')
         plastic_isotropic_gdot0(instance)        = IO_floatValue(line,chunkPos,2_pInt)
         if (plastic_isotropic_gdot0(instance) <= 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_ISOTROPIC_label//')')
       case ('n')
         plastic_isotropic_n(instance)            = IO_floatValue(line,chunkPos,2_pInt)
         if (plastic_isotropic_n(instance) <= 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_ISOTROPIC_label//')')
       case ('h0')
         plastic_isotropic_h0(instance)           = IO_floatValue(line,chunkPos,2_pInt)
       case ('h0_slope','slopelnrate')
         plastic_isotropic_h0_slopeLnRate(instance)  = IO_floatValue(line,chunkPos,2_pInt)
       case ('tausat')
         plastic_isotropic_tausat(instance)          = IO_floatValue(line,chunkPos,2_pInt)
         if (plastic_isotropic_tausat(instance) <= 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_ISOTROPIC_label//')')
       case ('tausat_sinhfita')
         plastic_isotropic_tausat_SinhFitA(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('tausat_sinhfitb')
         plastic_isotropic_tausat_SinhFitB(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('tausat_sinhfitc')
         plastic_isotropic_tausat_SinhFitC(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('tausat_sinhfitd')
         plastic_isotropic_tausat_SinhFitD(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('a', 'w0')
         plastic_isotropic_a(instance)               = IO_floatValue(line,chunkPos,2_pInt)
         if (plastic_isotropic_a(instance) <= 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_ISOTROPIC_label//')')
       case ('taylorfactor')
         plastic_isotropic_fTaylor(instance)         = IO_floatValue(line,chunkPos,2_pInt)
         if (plastic_isotropic_fTaylor(instance) <= 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_ISOTROPIC_label//')')
       case ('atol_resistance')
         plastic_isotropic_aTolResistance(instance)  = IO_floatValue(line,chunkPos,2_pInt)
         if (plastic_isotropic_aTolResistance(instance) <= 0.0_pReal) &
           call IO_error(211_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_ISOTROPIC_label//')')
       case ('atol_shear')
         plastic_isotropic_aTolShear(instance)  = IO_floatValue(line,chunkPos,2_pInt)

       case default

     end select
   endif; endif
 enddo parsingFile

 initializeInstances: do phase = 1_pInt, size(phase_plasticity)
   myPhase: if (phase_plasticity(phase) == PLASTICITY_isotropic_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_plasticityInstance(phase)
!--------------------------------------------------------------------------------------------------
!  sanity checks
     if (plastic_isotropic_aTolShear(instance) <= 0.0_pReal) &
       plastic_isotropic_aTolShear(instance) = 1.0e-6_pReal                                         ! default absolute tolerance 1e-6

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,plastic_isotropic_Noutput(instance)
       select case(plastic_isotropic_outputID(o,instance))
         case(flowstress_ID,strainrate_ID)
           mySize = 1_pInt
         case default
       end select
  
       outputFound: if (mySize > 0_pInt) then
         plastic_isotropic_sizePostResult(o,instance) = mySize
         plastic_isotropic_sizePostResults(instance) = &
         plastic_isotropic_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop

!--------------------------------------------------------------------------------------------------
! allocate state arrays
     sizeState    = 2_pInt
     sizeDotState = sizeState
     sizeDeltaState = 0_pInt
     plasticState(phase)%sizeState = sizeState
     plasticState(phase)%sizeDotState = sizeDotState
     plasticState(phase)%sizeDeltaState = sizeDeltaState
     plasticState(phase)%sizePostResults = plastic_isotropic_sizePostResults(instance)
     plasticState(phase)%nSlip = 1
     plasticState(phase)%nTwin = 0
     plasticState(phase)%nTrans= 0
     allocate(plasticState(phase)%aTolState          (   sizeState))
     plasticState(phase)%aTolState(1) = plastic_isotropic_aTolResistance(instance)
     plasticState(phase)%aTolState(2) = plastic_isotropic_aTolShear(instance)
     allocate(plasticState(phase)%state0             (   sizeState,NofMyPhase))
     plasticState(phase)%state0(1,1:NofMyPhase) = plastic_isotropic_tau0(instance)
     plasticState(phase)%state0(2,1:NofMyPhase) = 0.0_pReal
     allocate(plasticState(phase)%partionedState0    (   sizeState,NofMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%subState0          (   sizeState,NofMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%state              (   sizeState,NofMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%dotState           (sizeDotState,NofMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%deltaState       (sizeDeltaState,NofMyPhase),source=0.0_pReal)
     if (.not. analyticJaco) then
       allocate(plasticState(phase)%state_backup     (   sizeState,NofMyPhase),source=0.0_pReal)
       allocate(plasticState(phase)%dotState_backup  (sizeDotState,NofMyPhase),source=0.0_pReal)
     endif
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(phase)%previousDotState (sizeDotState,NofMyPhase),source=0.0_pReal)
       allocate(plasticState(phase)%previousDotState2(sizeDotState,NofMyPhase),source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(phase)%RK4dotState      (sizeDotState,NofMyPhase),source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(phase)%RKCK45dotState (6,sizeDotState,NofMyPhase),source=0.0_pReal)
     plasticState(phase)%slipRate        => plasticState(phase)%dotState(2:2,1:NofMyPhase)
     plasticState(phase)%accumulatedSlip => plasticState(phase)%state   (2:2,1:NofMyPhase)
   endif myPhase
 enddo initializeInstances

end subroutine plastic_isotropic_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_LpAndItsTangent(Lp,dLp_dTstar99,Tstar_v,ipc,ip,el)
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_g
 use math, only: &
   math_mul6x6, &
   math_Mandel6to33, &
   math_Plain3333to99, &
   math_deviatoric33, &
   math_mul33xx33, &
   math_transpose33
 use material, only: &
   phaseAt, phasememberAt, &
   plasticState, &
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
   gamma_dot = plastic_isotropic_gdot0(instance) &
             * (sqrt(1.5_pReal) * norm_Tstar_dev / (plastic_isotropic_fTaylor(instance) * &
               plasticState(phaseAt(ipc,ip,el))%state(1,phasememberAt(ipc,ip,el)))) &
                                                  **plastic_isotropic_n(instance)

   Lp = Tstar_dev_33/norm_Tstar_dev * gamma_dot/plastic_isotropic_fTaylor(instance)

   if (iand(debug_level(debug_constitutive), debug_levelExtensive) /= 0_pInt &
       .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
              .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
     write(6,'(a,i8,1x,i2,1x,i3)') '<< CONST isotropic >> at el ip g ',el,ip,ipc
     write(6,'(/,a,/,3(12x,3(f12.4,1x)/))') '<< CONST isotropic >> Tstar (dev) / MPa', &
                                      math_transpose33(Tstar_dev_33(1:3,1:3))*1.0e-6_pReal
     write(6,'(/,a,/,f12.5)') '<< CONST isotropic >> norm Tstar / MPa', norm_Tstar_dev*1.0e-6_pReal
     write(6,'(/,a,/,f12.5)') '<< CONST isotropic >> gdot', gamma_dot
   end if
!--------------------------------------------------------------------------------------------------
! Calculation of the tangent of Lp
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,l,m,n) = (plastic_isotropic_n(instance)-1.0_pReal) * &
                                      Tstar_dev_33(k,l)*Tstar_dev_33(m,n) / squarenorm_Tstar_dev
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,l,k,l) = dLp_dTstar_3333(k,l,k,l) + 1.0_pReal
   forall (k=1_pInt:3_pInt,m=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,k,m,m) = dLp_dTstar_3333(k,k,m,m) - 1.0_pReal/3.0_pReal
   dLp_dTstar99 = math_Plain3333to99(gamma_dot / plastic_isotropic_fTaylor(instance) * &
                                      dLp_dTstar_3333 / norm_Tstar_dev)
 end if
end subroutine plastic_isotropic_LpAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_LiAndItsTangent(Li,dLi_dTstar_3333,Tstar_v,ipc,ip,el)
 use math, only: &
   math_mul6x6, &
   math_Mandel6to33, &
   math_Plain3333to99, &
   math_spherical33, &
   math_mul33xx33
 use material, only: &
   phaseAt, phasememberAt, &
   plasticState, &
   material_phase, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3), intent(out) :: &
   Li                                                                                               !< plastic velocity gradient

 real(pReal), dimension(6),   intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),               intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 real(pReal), dimension(3,3) :: &
   Tstar_sph_33                                                                                     !< sphiatoric part of the 2nd Piola Kirchhoff stress tensor as 2nd order tensor
 real(pReal), dimension(3,3,3,3), intent(out)  :: &
   dLi_dTstar_3333                                                                                  !< derivative of Li with respect to Tstar as 4th order tensor
 real(pReal) :: &
   gamma_dot, &                                                                                     !< strainrate
   norm_Tstar_sph, &                                                                                !< euclidean norm of Tstar_sph
   squarenorm_Tstar_sph                                                                             !< square of the euclidean norm of Tstar_sph
 integer(pInt) :: &
   instance, &
   k, l, m, n

 instance = phase_plasticityInstance(material_phase(ipc,ip,el))

 Tstar_sph_33 = math_spherical33(math_Mandel6to33(Tstar_v))                                         ! spherical part of 2nd Piola-Kirchhoff stress
 squarenorm_Tstar_sph = math_mul33xx33(Tstar_sph_33,Tstar_sph_33)
 norm_Tstar_sph = sqrt(squarenorm_Tstar_sph) 

 if (plastic_isotropic_dilatation(instance)) then
     if (norm_Tstar_sph <= 0.0_pReal) then                                                          ! Tstar == 0 --> both Li and dLi_dTstar are zero
       Li = 0.0_pReal
       dLi_dTstar_3333 = 0.0_pReal
     else
       gamma_dot = plastic_isotropic_gdot0(instance) &
                   * (sqrt(1.5_pReal) * norm_Tstar_sph / (plastic_isotropic_fTaylor(instance) * &
                   plasticState(phaseAt(ipc,ip,el))%state(1,phasememberAt(ipc,ip,el)))) &
                                                           **plastic_isotropic_n(instance)

       Li = Tstar_sph_33/norm_Tstar_sph * gamma_dot/plastic_isotropic_fTaylor(instance)

       !--------------------------------------------------------------------------------------------------
       ! Calculation of the tangent of Li
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLi_dTstar_3333(k,l,m,n) = (plastic_isotropic_n(instance)-1.0_pReal) * &
                                          Tstar_sph_33(k,l)*Tstar_sph_33(m,n) / squarenorm_Tstar_sph
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
         dLi_dTstar_3333(k,l,k,l) = dLi_dTstar_3333(k,l,k,l) + 1.0_pReal

       dLi_dTstar_3333 = gamma_dot / plastic_isotropic_fTaylor(instance) * &
                                          dLi_dTstar_3333 / norm_Tstar_sph
     endif
 endif
 
end subroutine plastic_isotropic_LiAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_dotState(Tstar_v,ipc,ip,el)
 use math, only: &
   math_mul6x6
 use material, only: &
   phaseAt, phasememberAt, &
   plasticState, &
   material_phase, &
   phase_plasticityInstance
 
 implicit none
 real(pReal), dimension(6), intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),             intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(6) :: &
   Tstar_dev_v                                                                                      !< deviatoric 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal) :: &
   gamma_dot, &                                                                                     !< strainrate
   hardening, &                                                                                     !< hardening coefficient
   saturation, &                                                                                    !< saturation resistance
   norm_Tstar_v                                                                                     !< euclidean norm of Tstar_dev
 integer(pInt) :: &
   instance, &                                                                                      !< instance of my instance (unique number of my constitutive model)
   of, &                                                                                            !< shortcut notation for offset position in state array
   ph                                                                                               !< shortcut notation for phase ID (unique number of all phases, regardless of constitutive model)

 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(material_phase(ipc,ip,el))

!--------------------------------------------------------------------------------------------------
! norm of (deviatoric) 2nd Piola-Kirchhoff stress
 if (plastic_isotropic_dilatation(instance)) then
   norm_Tstar_v = sqrt(math_mul6x6(Tstar_v,Tstar_v))
 else
   Tstar_dev_v(1:3) = Tstar_v(1:3) - sum(Tstar_v(1:3))/3.0_pReal
   Tstar_dev_v(4:6) = Tstar_v(4:6)
   norm_Tstar_v = sqrt(math_mul6x6(Tstar_dev_v,Tstar_dev_v))
 end if
!--------------------------------------------------------------------------------------------------
! strain rate 
 gamma_dot = plastic_isotropic_gdot0(instance) * ( sqrt(1.5_pReal) * norm_Tstar_v & 
            / &!-----------------------------------------------------------------------------------
           (plastic_isotropic_fTaylor(instance)*plasticState(ph)%state(1,of)) )**plastic_isotropic_n(instance)
 
!--------------------------------------------------------------------------------------------------
! hardening coefficient
 if (abs(gamma_dot) > 1e-12_pReal) then
   if (abs(plastic_isotropic_tausat_SinhFitA(instance)) <= tiny(0.0_pReal)) then
     saturation = plastic_isotropic_tausat(instance)
   else
     saturation = (  plastic_isotropic_tausat(instance) &
                   + ( log(  ( gamma_dot / plastic_isotropic_tausat_SinhFitA(instance)&
                               )**(1.0_pReal / plastic_isotropic_tausat_SinhFitD(instance))&
                            + sqrt(  ( gamma_dot / plastic_isotropic_tausat_SinhFitA(instance) &
                                      )**(2.0_pReal / plastic_isotropic_tausat_SinhFitD(instance)) &
                                   + 1.0_pReal ) &
                            ) & ! asinh(K) = ln(K + sqrt(K^2 +1))
                       )**(1.0_pReal / plastic_isotropic_tausat_SinhFitC(instance)) &
                   / (  plastic_isotropic_tausat_SinhFitB(instance) &
                      * (gamma_dot / plastic_isotropic_gdot0(instance))**(1.0_pReal / plastic_isotropic_n(instance)) &
                      ) &
                   )
   endif
   hardening = ( plastic_isotropic_h0(instance) + plastic_isotropic_h0_slopeLnRate(instance) * log(gamma_dot) ) &
               * abs( 1.0_pReal - plasticState(ph)%state(1,of)/saturation )**plastic_isotropic_a(instance) &
               * sign(1.0_pReal, 1.0_pReal - plasticState(ph)%state(1,of)/saturation)
 else
   hardening = 0.0_pReal
 endif

  plasticState(ph)%dotState(1,of) = hardening * gamma_dot
  plasticState(ph)%dotState(2,of) =             gamma_dot

end subroutine plastic_isotropic_dotState

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_isotropic_postResults(Tstar_v,ipc,ip,el)
 use math, only: &
   math_mul6x6
 use material, only: &
   material_phase, &
   plasticState, &
   phaseAt, phasememberAt, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(plastic_isotropic_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
                                           plastic_isotropic_postResults

 real(pReal), dimension(6) :: &
   Tstar_dev_v                                                                                      !< deviatoric 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal) :: &
   norm_Tstar_v                                                                                     ! euclidean norm of Tstar_dev
 integer(pInt) :: &
   instance, &                                                                                      !< instance of my instance (unique number of my constitutive model)
   of, &                                                                                            !< shortcut notation for offset position in state array
   ph, &                                                                                            !< shortcut notation for phase ID (unique number of all phases, regardless of constitutive model)
   c, &
   o

 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 
!--------------------------------------------------------------------------------------------------
! norm of (deviatoric) 2nd Piola-Kirchhoff stress
 if (plastic_isotropic_dilatation(instance)) then
   norm_Tstar_v = sqrt(math_mul6x6(Tstar_v,Tstar_v))
 else
   Tstar_dev_v(1:3) = Tstar_v(1:3) - sum(Tstar_v(1:3))/3.0_pReal
   Tstar_dev_v(4:6) = Tstar_v(4:6)
   norm_Tstar_v = sqrt(math_mul6x6(Tstar_dev_v,Tstar_dev_v))
 end if
 
 c = 0_pInt
 plastic_isotropic_postResults = 0.0_pReal

 outputsLoop: do o = 1_pInt,plastic_isotropic_Noutput(instance)
   select case(plastic_isotropic_outputID(o,instance))
     case (flowstress_ID)
       plastic_isotropic_postResults(c+1_pInt) = plasticState(ph)%state(1,of)
       c = c + 1_pInt
     case (strainrate_ID)
       plastic_isotropic_postResults(c+1_pInt) = &
                plastic_isotropic_gdot0(instance) * (            sqrt(1.5_pReal) * norm_Tstar_v & 
             / &!----------------------------------------------------------------------------------
              (plastic_isotropic_fTaylor(instance) * plasticState(ph)%state(1,of)) ) ** plastic_isotropic_n(instance)
       c = c + 1_pInt
   end select
 enddo outputsLoop

end function plastic_isotropic_postResults


end module plastic_isotropic
