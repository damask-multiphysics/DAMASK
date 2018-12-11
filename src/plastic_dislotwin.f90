!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Su Leen Wong, Max-Planck-Institut für Eisenforschung GmbH
!> @author Nan Jia, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating dislocation and twinning physics
!> @details to be done
!--------------------------------------------------------------------------------------------------
module plastic_dislotwin
 use prec, only: &
   pReal, &
   pInt
    
 implicit none
 private
 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   plastic_dislotwin_sizePostResult                                                            !< size of each post result output
 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   plastic_dislotwin_output                                                                    !< name of each post result output
 
 real(pReal),                                                 parameter,           private :: &
   kB = 1.38e-23_pReal                                                                         !< Boltzmann constant in J/Kelvin

 enum, bind(c) 
   enumerator :: & 
     undefined_ID, &
     edge_density_ID, &
     dipole_density_ID, &
     shear_rate_slip_ID, &
     accumulated_shear_slip_ID, &
     mfp_slip_ID, &
     resolved_stress_slip_ID, &
     threshold_stress_slip_ID, &
     edge_dipole_distance_ID, &
     stress_exponent_ID, &
     twin_fraction_ID, &
     shear_rate_twin_ID, &
     accumulated_shear_twin_ID, &
     mfp_twin_ID, &
     resolved_stress_twin_ID, &
     threshold_stress_twin_ID, &
     resolved_stress_shearband_ID, &
     shear_rate_shearband_ID, &
     stress_trans_fraction_ID, &
     strain_trans_fraction_ID
 end enum
 
 type, private :: tParameters
   real(pReal) :: &
     mu, &
     nu, &
     CAtomicVolume, &                                                                               !< atomic volume in Bugers vector unit
     D0, &                                                                                          !< prefactor for self-diffusion coefficient
     Qsd, &                                                                                         !< activation energy for dislocation climb
     GrainSize, &                                                                                   !<grain size
     pShearBand, &                                                                                  !< p-exponent in shear band velocity
     qShearBand, &                                                                                  !< q-exponent in shear band velocity
     MaxTwinFraction, &                                                                             !<max allowed total twin volume fraction
     CEdgeDipMinDistance, &                                                                         !<
     Cmfptwin, &                                                                                    !<
     Cthresholdtwin, &                                                                              !<
     SolidSolutionStrength, &                                                                       !<strength due to elements in solid solution
     L0_twin, &                                                                                     !< Length of twin nuclei in Burgers vectors
     L0_trans, &                                                                                    !< Length of trans nuclei in Burgers vectors
     xc_twin, &                                                                                     !< critical distance for formation of twin nucleus
     xc_trans, &                                                                                    !< critical distance for formation of trans nucleus
     VcrossSlip, &                                                                                  !< cross slip volume
     sbResistance, &                                                                                !< value for shearband resistance (might become an internal state variable at some point)
     sbVelocity, &                                                                                  !< value for shearband velocity_0
     sbQedge, &                                                                                     !< value for shearband systems Qedge
     SFE_0K, &                                                                                      !< stacking fault energy at zero K
     dSFE_dT, &                                                                                     !< temperature dependance of stacking fault energy
     aTolRho, &                                                                                     !< absolute tolerance for integration of dislocation density
     aTolTwinFrac, &                                                                                !< absolute tolerance for integration of twin volume fraction
     aTolTransFrac, &                                                                               !< absolute tolerance for integration of trans volume fraction
     deltaG, &                                                                                      !< Free energy difference between austensite and martensite
     Cmfptrans, &                                                                                   !<
     Cthresholdtrans, &                                                                             !<
     transStackHeight                                                                               !< Stack height of hex nucleus 
   real(pReal),                  dimension(:),            allocatable :: & 
     rho0, &                                                                                        !< initial unipolar dislocation density per slip system
     rhoDip0, &                                                                                     !< initial dipole dislocation density per slip system
     burgers_slip, &                                                                                !< absolute length of burgers vector [m] for each slip system
     burgers_twin, &                                                                                !< absolute length of burgers vector [m] for each slip system
     burgers_trans, &                                                                               !< absolute length of burgers vector [m] for each twin system
     Qedge,&                                                                                        !< activation energy for glide [J] for each slip system
     v0, &                                                                                          !< dislocation velocity prefactor [m/s] for each slip system
     tau_peierls,&                                                                                  !< Peierls stress [Pa] for each slip system
     Ndot0_twin, &                                                                                  !< twin nucleation rate [1/m³s] for each twin system
     Ndot0_trans, &                                                                                 !< trans nucleation rate [1/m³s] for each trans system
     twinsize, &                                                                                    !< twin thickness [m] for each twin system
     CLambdaSlip, &                                                                                 !< Adj. parameter for distance between 2 forest dislocations for each slip system
     lamellarsizePerTransSystem, &                                                                  !< martensite lamellar thickness [m] for each trans system and instance
     p, &                                                                                           !< p-exponent in glide velocity
     q, &                                                                                           !< q-exponent in glide velocity
     r, &                                                                                           !< r-exponent in twin nucleation rate
     s, &                                                                                           !< s-exponent in trans nucleation rate
     shear_twin, &                                                                                  !< characteristic shear for twins
     B                                                                                              !< drag coefficient
   real(pReal),                  dimension(:,:),            allocatable :: & 
     interaction_SlipSlip, &                                                                        !< coefficients for slip-slip interaction for each interaction type
     interaction_SlipTwin, &                                                                        !< coefficients for slip-twin interaction for each interaction type
     interaction_TwinSlip, &                                                                        !< coefficients for twin-slip interaction for each interaction type
     interaction_TwinTwin, &                                                                        !< coefficients for twin-twin interaction for each interaction type
     interaction_SlipTrans, &                                                                       !< coefficients for slip-trans interaction for each interaction type
     interaction_TransTrans                                                                         !< coefficients for trans-trans interaction for each interaction type
   integer(pInt),                  dimension(:,:),            allocatable :: & 
     fcc_twinNucleationSlipPair                                                                     ! ToDo: Better name? Is also use for trans
   real(pReal),                  dimension(:,:),            allocatable :: & 
     forestProjectionEdge, &
     C66
   real(pReal),                  dimension(:,:,:),            allocatable :: &
     Schmid_trans, &
     Schmid_slip, &
     Schmid_twin, &
     C66_twin, &
     C66_trans
   logical :: &
     dipoleFormation                                                                                !< flag indicating consideration of dipole formation

   integer(kind(undefined_ID)),         dimension(:),         allocatable :: &
     outputID                                                                                       !< ID of each post result output
   
   logical :: &
     fccTwinTransNucleation                                                                             !< twinning and transformation models are for fcc
   integer(pInt) :: & 
     totalNslip, &                                                                                          !< number of active slip systems for each family and instance
     totalNtwin, &                                                                                          !< number of active twin systems for each family and instance
    totalNtrans                                                                                            !< number of active transformation systems for each family and instance 
   integer(pInt),                  dimension(:),            allocatable :: & 
     Nslip, &                                                                                          !< number of active slip systems for each family and instance
     Ntwin, &                                                                                          !< number of active twin systems for each family and instance
     Ntrans                                                                                            !< number of active transformation systems for each family and instance
  end type 
  
  type(tParameters), dimension(:), allocatable, private :: param                                    !< containers of constitutive parameters (len Ninstance)
 
 
 type, private :: tDislotwinState
   real(pReal), pointer, dimension(:,:) :: &
     rhoEdge, &
     rhoEdgeDip, &
     accshear_slip, &
     twinFraction, &
     accshear_twin, &
     stressTransFraction, &
     strainTransFraction, &
     whole
 end type tDislotwinState

 type, private :: tDislotwinMicrostructure
   real(pReal), allocatable,     dimension(:,:) :: &
     invLambdaSlip, &
     invLambdaSlipTwin, &
     invLambdaTwin, &
     invLambdaSlipTrans, &
     invLambdaTrans, &
     mfp_slip, &
     mfp_twin, &
     mfp_trans, &
     threshold_stress_slip, &
     threshold_stress_twin, &
     threshold_stress_trans, &
     twinVolume, &
     martensiteVolume, &
     tau_r_twin, &                                                             !< stress to bring partial close together for each twin system and instance
     tau_r_trans                                                               !< stress to bring partial close together for each trans system and instance
 end type tDislotwinMicrostructure

 type(tDislotwinState), allocatable, dimension(:), private :: &
   state, &
   dotState
 type(tDislotwinMicrostructure), allocatable, dimension(:), private :: &
   microstructure

 public :: &
   plastic_dislotwin_init, &
   plastic_dislotwin_homogenizedC, &
   plastic_dislotwin_microstructure, &
   plastic_dislotwin_LpAndItsTangent, &
   plastic_dislotwin_dotState, &
   plastic_dislotwin_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_init
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   pStringLen, &
   dEq0, &
   dNeq0, &
   dNeq
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use math, only: &
   math_rotate_forward3333, &
   math_Mandel3333to66, &
   math_mul3x3, &
   math_expand,&
   PI
 use IO, only: &
   IO_warning, &
   IO_error, &
   IO_timeStamp
 use material, only: &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   material_allocatePlasticState, &
   PLASTICITY_DISLOTWIN_label, &
   PLASTICITY_DISLOTWIN_ID, &
   material_phase, &  
   plasticState
 use config, only: &
   MATERIAL_partPhase, &
   config_phase
 use lattice

 implicit none
 integer(pInt) :: Ninstance,&
                  f,j,i,k,o,p, &
                  offset_slip, index_myFamily, index_otherFamily, &
                  startIndex, endIndex, outputSize
 integer(pInt) :: sizeState, sizeDotState
 integer(pInt) :: NipcMyPhase   
     
 real(pReal),  allocatable, dimension(:,:) :: temp1
 integer(pInt), dimension(1,200), parameter :: lattice_ntranssystem = 12 ! HACK!!
 integer(pInt),          dimension(0), parameter :: emptyIntArray    = [integer(pInt)::]
 real(pReal),            dimension(0), parameter :: emptyRealArray   = [real(pReal)::]
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]

 integer(kind(undefined_ID)) :: &
   outputID                                                                                         !< ID of each post result output

 character(len=pStringLen) :: &
   structure = '',&
   extmsg = ''
 character(len=65536), dimension(:), allocatable :: &
  outputs

 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_DISLOTWIN_label//' init  -+>>>'
 write(6,'(/,a)')   ' A. Ma and F. Roters, Acta Materialia, 52(12):3603–3612, 2004'
 write(6,'(a)')     ' https://doi.org/10.1016/j.actamat.2004.04.012'
 write(6,'(/,a)')   ' F.Roters et al., Computational Materials Science, 39:91–95, 2007'
 write(6,'(a)')     ' https://doi.org/10.1016/j.commatsci.2006.04.014'
 write(6,'(/,a)')   ' Wong et al., Acta Materialia, 118:140–151, 2016'
 write(6,'(a,/)')   ' https://doi.org/10.1016/j.actamat.2016.07.032'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 Ninstance = int(count(phase_plasticity == PLASTICITY_DISLOTWIN_ID),pInt)
 if (Ninstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',Ninstance
 

 allocate(plastic_dislotwin_sizePostResult(maxval(phase_Noutput),Ninstance),source=0_pInt)
 allocate(plastic_dislotwin_output(maxval(phase_Noutput),Ninstance))
          plastic_dislotwin_output = ''

 allocate(param(Ninstance))
 allocate(state(Ninstance))
 allocate(dotState(Ninstance))
 allocate(microstructure(Ninstance))

 do p = 1_pInt, size(phase_plasticityInstance)
   if (phase_plasticity(p) /= PLASTICITY_DISLOTWIN_ID) cycle
   associate(prm => param(phase_plasticityInstance(p)), &
             dot => dotState(phase_plasticityInstance(p)), &
             stt => state(phase_plasticityInstance(p)), &
             mse => microstructure(phase_plasticityInstance(p)))

   ! This data is read in already in lattice
   prm%mu = lattice_mu(p)
   prm%nu = lattice_nu(p)
   prm%C66 = lattice_C66(1:6,1:6,p)

   structure          = config_phase(p)%getString('lattice_structure')


!--------------------------------------------------------------------------------------------------
! slip related parameters
   prm%Nslip =  config_phase(p)%getInts('nslip',defaultVal=emptyIntArray)
   prm%totalNslip = sum(prm%Nslip)
   slipActive: if (prm%totalNslip > 0_pInt) then

     prm%fccTwinTransNucleation = merge(.true., .false., lattice_structure(p) == LATTICE_FCC_ID) &
                                .and. (prm%Nslip(1) == 12_pInt)
     if(prm%fccTwinTransNucleation) &
       prm%fcc_twinNucleationSlipPair = lattice_fcc_twinNucleationSlipPair


     prm%Schmid_slip          = lattice_SchmidMatrix_slip(prm%Nslip,structure(1:3),&
                                                          config_phase(p)%getFloat('c/a',defaultVal=0.0_pReal))
     prm%interaction_SlipSlip = lattice_interaction_SlipSlip(prm%Nslip, &
                                                             config_phase(p)%getFloats('interaction_slipslip'), &
                                                             structure(1:3))

     prm%rho0                 = config_phase(p)%getFloats('rhoedge0',   requiredShape=shape(prm%Nslip))  !ToDo: rename to rho_0
     prm%rhoDip0              = config_phase(p)%getFloats('rhoedgedip0',requiredShape=shape(prm%Nslip))  !ToDo: rename to rho_dip_0
     prm%v0                   = config_phase(p)%getFloats('v0',         requiredShape=shape(prm%Nslip))
     prm%burgers_slip         = config_phase(p)%getFloats('slipburgers',requiredShape=shape(prm%Nslip))
     prm%Qedge                = config_phase(p)%getFloats('qedge',      requiredShape=shape(prm%Nslip))  !ToDo: rename (ask Karo)
     prm%CLambdaSlip          = config_phase(p)%getFloats('clambdaslip',requiredShape=shape(prm%Nslip))
     prm%p                    = config_phase(p)%getFloats('p_slip',     requiredShape=shape(prm%Nslip))
     prm%q                    = config_phase(p)%getFloats('q_slip',     requiredShape=shape(prm%Nslip))
     prm%B                    = config_phase(p)%getFloats('b',          requiredShape=shape(prm%Nslip), &
                                                          defaultVal=[(0.0_pReal, i=1,size(prm%Nslip))])
     prm%tau_peierls          = config_phase(p)%getFloats('tau_peierls',requiredShape=shape(prm%Nslip), &
                                                          defaultVal=[(0.0_pReal, i=1,size(prm%Nslip))])

     prm%CEdgeDipMinDistance  = config_phase(p)%getFloat('cedgedipmindistance')

     ! expand: family => system
     prm%rho0         = math_expand(prm%rho0,        prm%Nslip)
     prm%rhoDip0      = math_expand(prm%rhoDip0,     prm%Nslip)
     prm%v0           = math_expand(prm%v0,          prm%Nslip)
     prm%burgers_slip = math_expand(prm%burgers_slip,prm%Nslip)
     prm%Qedge        = math_expand(prm%Qedge,       prm%Nslip)
     prm%CLambdaSlip  = math_expand(prm%CLambdaSlip, prm%Nslip)
     prm%p            = math_expand(prm%p,           prm%Nslip)
     prm%q            = math_expand(prm%q,           prm%Nslip)
     prm%B            = math_expand(prm%B,           prm%Nslip)
     prm%tau_peierls  = math_expand(prm%tau_peierls, prm%Nslip)

     ! sanity checks
     if (any(prm%rho0         <  0.0_pReal))         extmsg = trim(extmsg)//'rho0 '
     if (any(prm%rhoDip0      <  0.0_pReal))         extmsg = trim(extmsg)//'rhoDip0 '
     if (any(prm%v0           <  0.0_pReal))         extmsg = trim(extmsg)//'v0 '
     if (any(prm%burgers_slip <= 0.0_pReal))         extmsg = trim(extmsg)//'burgers_slip '
     if (any(prm%Qedge        <= 0.0_pReal))         extmsg = trim(extmsg)//'Qedge '
     if (any(prm%CLambdaSlip  <= 0.0_pReal))         extmsg = trim(extmsg)//'CLambdaSlip '
     if (any(prm%B            <  0.0_pReal))         extmsg = trim(extmsg)//'B '
     if (any(prm%tau_peierls  <  0.0_pReal))         extmsg = trim(extmsg)//'tau_peierls '
     if (any(prm%p<=0.0_pReal .or. prm%p>1.0_pReal)) extmsg = trim(extmsg)//'p '
     if (any(prm%q< 1.0_pReal .or. prm%q>2.0_pReal)) extmsg = trim(extmsg)//'q '

   else slipActive
     allocate(prm%burgers_slip(0))
   endif slipActive

!--------------------------------------------------------------------------------------------------
! twin related parameters
   prm%Ntwin      = config_phase(p)%getInts('ntwin', defaultVal=emptyIntArray)
   prm%totalNtwin = sum(prm%Ntwin)
   if (prm%totalNtwin > 0_pInt) then
     prm%Schmid_twin          = lattice_SchmidMatrix_twin(prm%Ntwin,structure(1:3),&
                                                          config_phase(p)%getFloat('c/a',defaultVal=0.0_pReal))
     prm%interaction_TwinTwin = lattice_interaction_TwinTwin(prm%Ntwin,&
                                                             config_phase(p)%getFloats('interaction_twintwin'), &
                                                             structure(1:3))

     prm%burgers_twin = config_phase(p)%getFloats('twinburgers')
     prm%twinsize     = config_phase(p)%getFloats('twinsize')
     prm%r            = config_phase(p)%getFloats('r_twin')

     prm%xc_twin         = config_phase(p)%getFloat('xc_twin')
     prm%L0_twin         = config_phase(p)%getFloat('l0_twin')
     prm%MaxTwinFraction = config_phase(p)%getFloat('maxtwinfraction') ! ToDo: only used in postResults
     prm%Cthresholdtwin  = config_phase(p)%getFloat('cthresholdtwin', defaultVal=0.0_pReal)
     prm%Cmfptwin        = config_phase(p)%getFloat('cmfptwin',       defaultVal=0.0_pReal) ! ToDo: How to handle that???

     prm%shear_twin      = lattice_characteristicShear_Twin(prm%Ntwin,structure(1:3),&
                                                          config_phase(p)%getFloat('c/a',defaultVal=0.0_pReal))

     prm%C66_twin        = lattice_C66_twin(prm%Ntwin,prm%C66,structure(1:3),&
                                                          config_phase(p)%getFloat('c/a',defaultVal=0.0_pReal))

     if (.not. prm%fccTwinTransNucleation) then
       prm%Ndot0_twin = config_phase(p)%getFloats('ndot0_twin') 
       prm%Ndot0_twin = math_expand(prm%Ndot0_twin,prm%Ntwin)
     endif

     ! expand: family => system
     prm%burgers_twin = math_expand(prm%burgers_twin,prm%Ntwin)
     prm%twinsize     = math_expand(prm%twinsize,prm%Ntwin)
     prm%r            = math_expand(prm%r,prm%Ntwin)
     
   else
     allocate(prm%twinsize(0))
     allocate(prm%burgers_twin(0))
     allocate(prm%r(0))
   endif
  
!--------------------------------------------------------------------------------------------------
! transformation related parameters
   prm%Ntrans      = config_phase(p)%getInts('ntrans', defaultVal=emptyIntArray)
   prm%totalNtrans = sum(prm%Ntrans)
   if (prm%totalNtrans > 0_pInt) then
     prm%burgers_trans = config_phase(p)%getFloats('transburgers')
     prm%burgers_trans = math_expand(prm%burgers_trans,prm%Ntrans)
     
     prm%Cthresholdtrans  = config_phase(p)%getFloat('cthresholdtrans', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%transStackHeight = config_phase(p)%getFloat('transstackheight', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%Cmfptrans        = config_phase(p)%getFloat('cmfptrans', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%deltaG           = config_phase(p)%getFloat('deltag')
     prm%xc_trans         = config_phase(p)%getFloat('xc_trans', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%L0_trans         = config_phase(p)%getFloat('l0_trans')

     prm%interaction_TransTrans = lattice_interaction_TransTrans(prm%Ntrans,&
                                                             config_phase(p)%getFloats('interaction_transtrans'), &
                                                             structure(1:3))
     if (lattice_structure(p) /= LATTICE_fcc_ID) then
        prm%Ndot0_trans = config_phase(p)%getFloats('ndot0_trans')
        prm%Ndot0_trans = math_expand(prm%Ndot0_trans,prm%Ntrans)
     endif
     prm%lamellarsizePerTransSystem = config_phase(p)%getFloats('lamellarsize')
     prm%lamellarsizePerTransSystem = math_expand(prm%lamellarsizePerTransSystem,prm%Ntrans)
     prm%s = config_phase(p)%getFloats('s_trans',defaultVal=[0.0_pReal])
     prm%s = math_expand(prm%s,prm%Ntrans)
   else
     allocate(prm%lamellarsizePerTransSystem(0))
     allocate(prm%burgers_trans(0))
   endif
   
   if (sum(prm%Ntwin) > 0_pInt  .or. prm%totalNtrans > 0_pInt) then
     prm%SFE_0K     = config_phase(p)%getFloat('sfe_0k')
     prm%dSFE_dT    = config_phase(p)%getFloat('dsfe_dt')
     prm%VcrossSlip = config_phase(p)%getFloat('vcrossslip')
   endif
   
   if (prm%totalNslip > 0_pInt .and. prm%totalNtwin > 0_pInt) then
     prm%interaction_SlipTwin = lattice_interaction_SlipTwin(prm%Nslip,prm%Ntwin,&
                                                             config_phase(p)%getFloats('interaction_sliptwin'), &
                                                             structure(1:3)) 
     prm%interaction_TwinSlip = lattice_interaction_TwinSlip(prm%Ntwin,prm%Nslip,&
                                                             config_phase(p)%getFloats('interaction_twinslip'), &
                                                             structure(1:3))
     if (prm%fccTwinTransNucleation .and. prm%totalNtwin > 12_pInt) write(6,*) 'mist' ! ToDo: implement better test. The model will fail also if ntwin is [6,6]
   endif    

   if (prm%totalNslip > 0_pInt .and. prm%totalNtrans > 0_pInt) then  
     prm%interaction_SlipTrans = lattice_interaction_SlipTrans(prm%Nslip,prm%Ntrans,&
                                                               config_phase(p)%getFloats('interaction_sliptrans'), &
                                                               structure(1:3)) 
     if (prm%fccTwinTransNucleation .and. prm%totalNtrans > 12_pInt) write(6,*) 'mist' ! ToDo: implement better test. The model will fail also if ntrans is [6,6]
   endif    


   prm%aTolRho       = config_phase(p)%getFloat('atol_rho', defaultVal=0.0_pReal)
   prm%aTolTwinFrac  = config_phase(p)%getFloat('atol_twinfrac', defaultVal=0.0_pReal)
   prm%aTolTransFrac = config_phase(p)%getFloat('atol_transfrac', defaultVal=0.0_pReal)

   prm%CAtomicVolume = config_phase(p)%getFloat('catomicvolume')
   prm%GrainSize =  config_phase(p)%getFloat('grainsize')


   prm%D0 = config_phase(p)%getFloat('d0')
   prm%Qsd = config_phase(p)%getFloat('qsd')
   prm%SolidSolutionStrength = config_phase(p)%getFloat('solidsolutionstrength')
   if (config_phase(p)%keyExists('dipoleformationfactor')) call IO_error(1,ext_msg='use /nodipoleformation/')
   prm%dipoleformation = .not. config_phase(p)%keyExists('/nodipoleformation/')
   prm%sbVelocity   = config_phase(p)%getFloat('shearbandvelocity',defaultVal=0.0_pReal)
   if (prm%sbVelocity > 0.0_pReal) then  
     prm%sbResistance = config_phase(p)%getFloat('shearbandresistance')
     prm%sbQedge = config_phase(p)%getFloat('qedgepersbsystem')
     prm%pShearBand = config_phase(p)%getFloat('p_shearband')
     prm%qShearBand = config_phase(p)%getFloat('q_shearband')
   endif

       !if (Ndot0PerTwinFamily(f,p) < 0.0_pReal) &
        ! call IO_error(211_pInt,el=p,ext_msg='ndot0_twin ('//PLASTICITY_DISLOTWIN_label//')')

   if (prm%CAtomicVolume <= 0.0_pReal) &
     call IO_error(211_pInt,el=p,ext_msg='cAtomicVolume ('//PLASTICITY_DISLOTWIN_label//')')
   if (prm%D0 <= 0.0_pReal) &
     call IO_error(211_pInt,el=p,ext_msg='D0 ('//PLASTICITY_DISLOTWIN_label//')')
   if (prm%Qsd <= 0.0_pReal) &
     call IO_error(211_pInt,el=p,ext_msg='Qsd ('//PLASTICITY_DISLOTWIN_label//')')
   if (prm%totalNtwin > 0_pInt) then
     if (dEq0(prm%SFE_0K) .and. &
         dEq0(prm%dSFE_dT) .and. &
                         lattice_structure(p) == LATTICE_fcc_ID) &
       call IO_error(211_pInt,el=p,ext_msg='SFE0K ('//PLASTICITY_DISLOTWIN_label//')')
     if (prm%aTolRho <= 0.0_pReal) &
       call IO_error(211_pInt,el=p,ext_msg='aTolRho ('//PLASTICITY_DISLOTWIN_label//')')   
     if (prm%aTolTwinFrac <= 0.0_pReal) &
       call IO_error(211_pInt,el=p,ext_msg='aTolTwinFrac ('//PLASTICITY_DISLOTWIN_label//')')
   endif
   if (prm%totalNtrans > 0_pInt) then
     if (dEq0(prm%SFE_0K) .and. &
         dEq0(prm%dSFE_dT) .and. &
                         lattice_structure(p) == LATTICE_fcc_ID) &
       call IO_error(211_pInt,el=p,ext_msg='SFE0K ('//PLASTICITY_DISLOTWIN_label//')')
     if (prm%aTolTransFrac <= 0.0_pReal) &
       call IO_error(211_pInt,el=p,ext_msg='aTolTransFrac ('//PLASTICITY_DISLOTWIN_label//')')
   endif
   !if (prm%sbResistance < 0.0_pReal) &
   !  call IO_error(211_pInt,el=p,ext_msg='sbResistance ('//PLASTICITY_DISLOTWIN_label//')')
   !if (prm%sbVelocity < 0.0_pReal) &
   !  call IO_error(211_pInt,el=p,ext_msg='sbVelocity ('//PLASTICITY_DISLOTWIN_label//')')
   !if (prm%sbVelocity > 0.0_pReal .and. &
   !    prm%pShearBand <= 0.0_pReal) &
   !  call IO_error(211_pInt,el=p,ext_msg='pShearBand ('//PLASTICITY_DISLOTWIN_label//')')
   if (prm%sbVelocity > 0.0_pReal .and. &
       prm%qShearBand <= 0.0_pReal) &
     call IO_error(211_pInt,el=p,ext_msg='qShearBand ('//PLASTICITY_DISLOTWIN_label//')')
 
   outputs = config_phase(p)%getStrings('(output)', defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   do i= 1_pInt, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
       case ('edge_density')
         outputID = merge(edge_density_ID,undefined_ID,prm%totalNslip > 0_pInt)
         outputSize = prm%totalNslip
       case ('dipole_density')
         outputID = merge(dipole_density_ID,undefined_ID,prm%totalNslip > 0_pInt)
         outputSize = prm%totalNslip
       case ('shear_rate_slip','shearrate_slip')
         outputID = merge(shear_rate_slip_ID,undefined_ID,prm%totalNslip > 0_pInt)
         outputSize = prm%totalNslip
       case ('accumulated_shear_slip')
         outputID = merge(accumulated_shear_slip_ID,undefined_ID,prm%totalNslip > 0_pInt)
         outputSize = prm%totalNslip
       case ('mfp_slip')
         outputID = merge(mfp_slip_ID,undefined_ID,prm%totalNslip > 0_pInt)
         outputSize = prm%totalNslip
       case ('resolved_stress_slip')
         outputID = merge(resolved_stress_slip_ID,undefined_ID,prm%totalNslip > 0_pInt)
         outputSize = prm%totalNslip
       case ('threshold_stress_slip')
         outputID= merge(threshold_stress_slip_ID,undefined_ID,prm%totalNslip > 0_pInt)
         outputSize = prm%totalNslip
       case ('edge_dipole_distance')
         outputID = merge(edge_dipole_distance_ID,undefined_ID,prm%totalNslip > 0_pInt)
         outputSize = prm%totalNslip
       case ('stress_exponent')
         outputID = merge(stress_exponent_ID,undefined_ID,prm%totalNslip > 0_pInt)
         outputSize = prm%totalNslip

       case ('twin_fraction')
         outputID = merge(twin_fraction_ID,undefined_ID,prm%totalNtwin >0_pInt)
         outputSize = prm%totalNtwin
       case ('shear_rate_twin','shearrate_twin')
         outputID = merge(shear_rate_twin_ID,undefined_ID,prm%totalNtwin >0_pInt)
         outputSize = prm%totalNtwin
       case ('accumulated_shear_twin')
         outputID = merge(accumulated_shear_twin_ID,undefined_ID,prm%totalNtwin >0_pInt)
         outputSize = prm%totalNtwin
       case ('mfp_twin')
         outputID = merge(mfp_twin_ID,undefined_ID,prm%totalNtwin >0_pInt)
         outputSize = prm%totalNtwin
       case ('resolved_stress_twin')
         outputID = merge(resolved_stress_twin_ID,undefined_ID,prm%totalNtwin >0_pInt)
         outputSize = prm%totalNtwin
       case ('threshold_stress_twin')
         outputID = merge(threshold_stress_twin_ID,undefined_ID,prm%totalNtwin >0_pInt)
         outputSize = prm%totalNtwin
         
       case ('resolved_stress_shearband')
         outputID = resolved_stress_shearband_ID
         outputSize = 6_pInt
       case ('shear_rate_shearband','shearrate_shearband')
         outputID = shear_rate_shearband_ID
         outputSize = 6_pInt
         
       case ('stress_trans_fraction')
         outputID = stress_trans_fraction_ID
         outputSize = prm%totalNtrans
       case ('strain_trans_fraction')
         outputID = strain_trans_fraction_ID
         outputSize = prm%totalNtrans
        
     end select
        
     if (outputID /= undefined_ID) then
       plastic_dislotwin_output(i,phase_plasticityInstance(p)) = outputs(i)
       plastic_dislotwin_sizePostResult(i,phase_plasticityInstance(p)) = outputSize
       prm%outputID = [prm%outputID, outputID]
     endif
   enddo

   
!--------------------------------------------------------------------------------------------------
! allocate state arrays
   NipcMyPhase=count(material_phase==p)
   sizeDotState     = int(size(['rho         ','rhoDip      ','accshearslip']),pInt) * prm%totalNslip &
                    + int(size(['twinFraction','accsheartwin']),pInt)                * prm%totalNtwin &
                    + int(size(['stressTransFraction','strainTransFraction']),pInt)  * prm%totalNtrans
   sizeState        = sizeDotState

   call material_allocatePlasticState(p,NipcMyPhase,sizeState,sizeDotState,0_pInt, &
                                      prm%totalNslip,prm%totalNtwin,prm%totalNtrans)
   plasticState(p)%sizePostResults = sum(plastic_dislotwin_sizePostResult(:,phase_plasticityInstance(p)))

   ! ToDo: do later on
   offset_slip = 2_pInt*plasticState(p)%nslip
   plasticState(p)%slipRate        => &
       plasticState(p)%dotState(offset_slip+1:offset_slip+plasticState(p)%nslip,1:NipcMyPhase)
   plasticState(p)%accumulatedSlip => &
       plasticState(p)%state   (offset_slip+1:offset_slip+plasticState(p)%nslip,1:NipcMyPhase)


! DEPRECATED BEGIN
   allocate(temp1(prm%totalNslip,prm%totalNtrans),source =0.0_pReal)        
   allocate(prm%forestProjectionEdge(prm%totalNslip,prm%totalNslip),source   = 0.0_pReal)
   i = 0_pInt
   mySlipFamilies: do f = 1_pInt,size(prm%Nslip,1)
     index_myFamily = sum(prm%Nslip(1:f-1_pInt))
     slipSystemsLoop: do j = 1_pInt,prm%Nslip(f)
       i = i + 1_pInt
       do o = 1_pInt, size(prm%Nslip,1)
         index_otherFamily = sum(prm%Nslip(1:o-1_pInt))
         do k = 1_pInt,prm%Nslip(o)                                   ! loop over (active) systems in other family (slip)
           prm%forestProjectionEdge(index_myFamily+j,index_otherFamily+k) = &
             abs(math_mul3x3(lattice_sn(:,sum(lattice_NslipSystem(1:f-1,p))+j,p), &
                             lattice_st(:,sum(lattice_NslipSystem(1:o-1,p))+k,p)))
       enddo; enddo
     enddo slipSystemsLoop
   enddo mySlipFamilies 

   allocate(prm%C66_trans(6,6,prm%totalNtrans)     ,source=0.0_pReal)
   allocate(prm%Schmid_trans(3,3,prm%totalNtrans),source  = 0.0_pReal)
   i = 0_pInt
   transFamiliesLoop: do f = 1_pInt,size(prm%Ntrans,1)
     index_myFamily = sum(prm%Ntrans(1:f-1_pInt))                                                   ! index in truncated trans system list
     transSystemsLoop: do j = 1_pInt,prm%Ntrans(f)
       i = i + 1_pInt
       prm%Schmid_trans(1:3,1:3,i) = lattice_Strans(1:3,1:3,sum(lattice_Ntranssystem(1:f-1,p))+j,p)
     !* Rotate trans elasticity matrices
       index_otherFamily = sum(lattice_NtransSystem(1:f-1_pInt,p))                                  ! index in full lattice trans list
       prm%C66_trans(1:6,1:6,index_myFamily+j) = &
         math_Mandel3333to66(math_rotate_forward3333(lattice_trans_C3333(1:3,1:3,1:3,1:3,p),&
                                                     lattice_Qtrans(1:3,1:3,index_otherFamily+j,p)))
     enddo transSystemsLoop
   enddo transFamiliesLoop
! DEPRECATED END  

   startIndex=1_pInt
   endIndex=prm%totalNslip
   stt%rhoEdge=>plasticState(p)%state(startIndex:endIndex,:)
   stt%rhoEdge= spread(prm%rho0,2,NipcMyPhase)
   dot%rhoEdge=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolRho

   startIndex=endIndex+1
   endIndex=endIndex+prm%totalNslip
   stt%rhoEdgeDip=>plasticState(p)%state(startIndex:endIndex,:)
   stt%rhoEdgeDip= spread(prm%rhoDip0,2,NipcMyPhase)
   dot%rhoEdgeDip=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolRho
   
   startIndex=endIndex+1
   endIndex=endIndex+prm%totalNslip
   stt%accshear_slip=>plasticState(p)%state(startIndex:endIndex,:)
   dot%accshear_slip=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = 1.0e6_pReal
   
   startIndex=endIndex+1
   endIndex=endIndex+prm%totalNtwin
   stt%twinFraction=>plasticState(p)%state(startIndex:endIndex,:)
   dot%twinFraction=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolTwinFrac
   
   startIndex=endIndex+1
   endIndex=endIndex+prm%totalNtwin
   stt%accshear_twin=>plasticState(p)%state(startIndex:endIndex,:)
   dot%accshear_twin=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = 1.0e6_pReal
   
   startIndex=endIndex+1
   endIndex=endIndex+prm%totalNtrans
   stt%stressTransFraction=>plasticState(p)%state(startIndex:endIndex,:)
   dot%stressTransFraction=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolTransFrac
   
   startIndex=endIndex+1
   endIndex=endIndex+prm%totalNtrans
   stt%strainTransFraction=>plasticState(p)%state(startIndex:endIndex,:)
   dot%strainTransFraction=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolTransFrac

   plasticState(p)%state0 = plasticState(p)%state
   dot%whole => plasticState(p)%dotState

   allocate(mse%invLambdaSlip         (prm%totalNslip, NipcMyPhase),source=0.0_pReal)
   allocate(mse%invLambdaSlipTwin     (prm%totalNslip, NipcMyPhase),source=0.0_pReal)
   allocate(mse%invLambdaSlipTrans    (prm%totalNslip, NipcMyPhase),source=0.0_pReal)
   allocate(mse%mfp_slip              (prm%totalNslip, NipcMyPhase),source=0.0_pReal)
   allocate(mse%threshold_stress_slip (prm%totalNslip, NipcMyPhase),source=0.0_pReal)

   allocate(mse%invLambdaTwin         (prm%totalNtwin, NipcMyPhase),source=0.0_pReal)
   allocate(mse%mfp_twin              (prm%totalNtwin, NipcMyPhase),source=0.0_pReal)
   allocate(mse%threshold_stress_twin (prm%totalNtwin, NipcMyPhase),source=0.0_pReal)
   allocate(mse%tau_r_twin            (prm%totalNtwin, NipcMyPhase),source=0.0_pReal)
   allocate(mse%twinVolume            (prm%totalNtwin, NipcMyPhase),source=0.0_pReal)

   allocate(mse%invLambdaTrans        (prm%totalNtrans,NipcMyPhase),source=0.0_pReal)
   allocate(mse%mfp_trans             (prm%totalNtrans,NipcMyPhase),source=0.0_pReal)
   allocate(mse%threshold_stress_trans(prm%totalNtrans,NipcMyPhase),source=0.0_pReal)
   allocate(mse%tau_r_trans           (prm%totalNtrans,NipcMyPhase),source=0.0_pReal)
   allocate(mse%martensiteVolume      (prm%totalNtrans,NipcMyPhase),source=0.0_pReal)

   end associate
 enddo
 
end subroutine plastic_dislotwin_init

!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenized elasticity matrix
!--------------------------------------------------------------------------------------------------
function plastic_dislotwin_homogenizedC(ipc,ip,el)
 use material, only: &
  material_phase, &
  phase_plasticityInstance, &
  phasememberAt
 
 implicit none
 real(pReal), dimension(6,6) :: &
   plastic_dislotwin_homogenizedC
 integer(pInt), intent(in) :: &
   ipc, &                                                                                          !< component-ID of integration point
   ip, &                                                                                           !< integration point
   el                                                                                              !< element

 integer(pInt) :: i, &
                  of
 real(pReal) :: f_unrotated

 of = phasememberAt(ipc,ip,el)
 associate(prm => param(phase_plasticityInstance(material_phase(ipc,ip,el))),&
           stt => state(phase_plasticityInstance(material_phase(ipc,ip,el))))

 f_unrotated = 1.0_pReal &
             - sum(stt%twinFraction(1_pInt:prm%totalNtwin,of)) &
             - sum(stt%stressTransFraction(1_pInt:prm%totalNtrans,of)) &
             - sum(stt%strainTransFraction(1_pInt:prm%totalNtrans,of))

 plastic_dislotwin_homogenizedC = f_unrotated * prm%C66
 do i=1_pInt,prm%totalNtwin
    plastic_dislotwin_homogenizedC = plastic_dislotwin_homogenizedC &
                                   + stt%twinFraction(i,of)*prm%C66_twin(1:6,1:6,i)
 enddo
 do i=1_pInt,prm%totalNtrans
    plastic_dislotwin_homogenizedC = plastic_dislotwin_homogenizedC &
                                   +(stt%stressTransFraction(i,of)+stt%strainTransFraction(i,of))*&
                                                            prm%C66_trans(1:6,1:6,i)
 enddo
 end associate
 end function plastic_dislotwin_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_microstructure(temperature,ipc,ip,el)
 use math, only: &
   PI
 use material, only: &
   material_phase, &
   phase_plasticityInstance, &
   phasememberAt

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in) :: &
   temperature                                                                                      !< temperature at IP 

 integer(pInt) :: &
   i, &
   of
 real(pReal) :: &
   sumf_twin,SFE,sumf_trans
 real(pReal), dimension(:), allocatable :: &
   x0, &
   fOverStacksize, &
   ftransOverLamellarSize

 of = phasememberAt(ipc,ip,el)

 associate(prm => param(phase_plasticityInstance(material_phase(ipc,ip,el))),&
           stt => state(phase_plasticityInstance(material_phase(ipc,ip,el))),&
           mse => microstructure(phase_plasticityInstance(material_phase(ipc,ip,el))))

 sumf_twin  = sum(stt%twinFraction(1:prm%totalNtwin,of))
 sumf_trans = sum(stt%stressTransFraction(1:prm%totalNtrans,of)) &
            + sum(stt%strainTransFraction(1:prm%totalNtrans,of))

 sfe = prm%SFE_0K + prm%dSFE_dT * Temperature
 
 !* rescaled volume fraction for topology
 fOverStacksize         =  stt%twinFraction(1_pInt:prm%totalNtwin,of)/prm%twinsize  !ToDo: this is per system
 ftransOverLamellarSize =  sumf_trans/prm%lamellarsizePerTransSystem                !ToDo: But this not ...
 !Todo: Physically ok, but naming could be adjusted


 !* 1/mean free distance between 2 forest dislocations seen by a moving dislocation
 forall (i = 1_pInt:prm%totalNslip) &
   mse%invLambdaSlip(i,of) = &
     sqrt(dot_product((stt%rhoEdge(1_pInt:prm%totalNslip,of)+stt%rhoEdgeDip(1_pInt:prm%totalNslip,of)),&
                      prm%forestProjectionEdge(1:prm%totalNslip,i)))/prm%CLambdaSlip(i)

 !* 1/mean free distance between 2 twin stacks from different systems seen by a moving dislocation
 !$OMP CRITICAL (evilmatmul)
 if (prm%totalNtwin > 0_pInt .and. prm%totalNslip > 0_pInt) &
   mse%invLambdaSlipTwin(1_pInt:prm%totalNslip,of) = &
     matmul(prm%interaction_SlipTwin,fOverStacksize)/(1.0_pReal-sumf_twin)

 !* 1/mean free distance between 2 twin stacks from different systems seen by a growing twin

  !ToDo: needed? if (prm%totalNtwin > 0_pInt) &
  mse%invLambdaTwin(1_pInt:prm%totalNtwin,of) = &
     matmul(prm%interaction_TwinTwin,fOverStacksize)/(1.0_pReal-sumf_twin)


 !* 1/mean free distance between 2 martensite lamellar from different systems seen by a moving dislocation
 if (prm%totalNtrans > 0_pInt .and. prm%totalNslip > 0_pInt) &
   mse%invLambdaSlipTrans(1_pInt:prm%totalNslip,of) = &                                  ! ToDo: does not work if Ntrans is not 12
      matmul(prm%interaction_SlipTrans,ftransOverLamellarSize)/(1.0_pReal-sumf_trans)

 !* 1/mean free distance between 2 martensite stacks from different systems seen by a growing martensite (1/lambda_trans)
 !ToDo: needed? if (prm%totalNtrans > 0_pInt) &

   mse%invLambdaTrans(1_pInt:prm%totalNtrans,of) = &
     matmul(prm%interaction_TransTrans,ftransOverLamellarSize)/(1.0_pReal-sumf_trans)
 !$OMP END CRITICAL (evilmatmul)

 !* mean free path between 2 obstacles seen by a moving dislocation
 do i = 1_pInt,prm%totalNslip
    if ((prm%totalNtwin > 0_pInt) .or. (prm%totalNtrans > 0_pInt)) then              ! ToDo: This is too simplified
       mse%mfp_slip(i,of) = &
         prm%GrainSize/(1.0_pReal+prm%GrainSize*&
         (mse%invLambdaSlip(i,of) + mse%invLambdaSlipTwin(i,of) + mse%invLambdaSlipTrans(i,of)))
    else
       mse%mfp_slip(i,of) = &  
         prm%GrainSize/&
         (1.0_pReal+prm%GrainSize*(mse%invLambdaSlip(i,of))) !!!!!! correct?
    endif
 enddo

 !* mean free path between 2 obstacles seen by a growing twin/martensite
 mse%mfp_twin(:,of)  = prm%Cmfptwin*prm%GrainSize/ (1.0_pReal+prm%GrainSize*mse%invLambdaTwin(:,of))
 mse%mfp_trans(:,of) = prm%Cmfptrans*prm%GrainSize/(1.0_pReal+prm%GrainSize*mse%invLambdaTrans(:,of))

 !* threshold stress for dislocation motion
 forall (i = 1_pInt:prm%totalNslip) mse%threshold_stress_slip(i,of) = &
     prm%mu*prm%burgers_slip(i)*&
     sqrt(dot_product(stt%rhoEdge(1_pInt:prm%totalNslip,of)+stt%rhoEdgeDip(1_pInt:prm%totalNslip,of),&
                      prm%interaction_SlipSlip(i,1:prm%totalNslip)))

 !* threshold stress for growing twin/martensite
 if(prm%totalNtwin == prm%totalNslip) &
 mse%threshold_stress_twin(:,of) = prm%Cthresholdtwin* &
     (sfe/(3.0_pReal*prm%burgers_twin)+ 3.0_pReal*prm%burgers_twin*prm%mu/ &
           (prm%L0_twin*prm%burgers_slip)) ! slip burgers here correct?
 if(prm%totalNtrans == prm%totalNslip) &
   mse%threshold_stress_trans(:,of) = prm%Cthresholdtrans* &
         (sfe/(3.0_pReal*prm%burgers_trans) + 3.0_pReal*prm%burgers_trans*prm%mu/&
              (prm%L0_trans*prm%burgers_slip) + prm%transStackHeight*prm%deltaG/ (3.0_pReal*prm%burgers_trans) )    
 
  ! final volume after growth
 mse%twinVolume(:,of) = (PI/4.0_pReal)*prm%twinsize*mse%mfp_twin(:,of)**2.0_pReal
 mse%martensiteVolume(:,of) = (PI/4.0_pReal)*prm%lamellarsizePerTransSystem*mse%mfp_trans(:,of)**2.0_pReal

 !* equilibrium separation of partial dislocations (twin)
 x0 = prm%mu*prm%burgers_twin**2.0_pReal/(sfe*8.0_pReal*PI)*(2.0_pReal+prm%nu)/(1.0_pReal-prm%nu)
 mse%tau_r_twin(:,of) = prm%mu*prm%burgers_twin/(2.0_pReal*PI)*(1.0_pReal/(x0+prm%xc_twin)+cos(pi/3.0_pReal)/x0)

 !* equilibrium separation of partial dislocations (trans)
 x0 = prm%mu*prm%burgers_trans**2.0_pReal/(sfe*8.0_pReal*PI)*(2.0_pReal+prm%nu)/(1.0_pReal-prm%nu)
 mse%tau_r_trans(:,of) = prm%mu*prm%burgers_trans/(2.0_pReal*PI)*(1.0_pReal/(x0+prm%xc_trans)+cos(pi/3.0_pReal)/x0)

end associate
end subroutine plastic_dislotwin_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_LpAndItsTangent(Lp,dLp_dMp,Mp,Temperature,instance,of)
 use prec, only: &
   tol_math_check, &
   dNeq0
 use math, only: &
   math_eigenValuesVectorsSym, &
   math_tensorproduct33, &
   math_symmetric33, &
   math_mul33xx33, &
   math_mul33x3
 
 implicit none
 real(pReal), dimension(3,3),     intent(out) :: Lp
 real(pReal), dimension(3,3,3,3), intent(out) :: dLp_dMp
 real(pReal), dimension(3,3),     intent(in)  :: Mp
 integer(pInt),                   intent(in)  :: instance,of
 real(pReal),                     intent(in)  :: Temperature

 integer(pInt) :: i,k,l,m,n,s1,s2
 real(pReal) :: f_unrotated,StressRatio_p,&
                StressRatio_r,BoltzmannRatio,Ndot0_twin,stressRatio, &
    Ndot0_trans,StressRatio_s, &
    dgdot_dtau, &
    tau
 real(pReal), dimension(param(instance)%totalNslip) :: &
    gdot_slip,dgdot_dtau_slip
 real(pReal), dimension(param(instance)%totalNtwin) :: &
    gdot_twin,dgdot_dtau_twin
 real(pReal):: gdot_sb,gdot_trans
 real(pReal), dimension(3,3) :: eigVectors, Schmid_shearBand
 real(pReal), dimension(3)   :: eigValues, sb_s, sb_m
 logical :: error
 real(pReal), dimension(3,6), parameter :: &
   sb_sComposition = &
     reshape(real([&
        1, 0, 1, &
        1, 0,-1, &
        1, 1, 0, &
        1,-1, 0, &
        0, 1, 1, &
        0, 1,-1  &
        ],pReal),[ 3,6]), &
   sb_mComposition = &
     reshape(real([&
        1, 0,-1, &
        1, 0,+1, &
        1,-1, 0, &
        1, 1, 0, &
        0, 1,-1, &
        0, 1, 1  &
        ],pReal),[ 3,6])

 associate(prm => param(instance), stt => state(instance), mse => microstructure(instance))

 f_unrotated = 1.0_pReal &
             - sum(stt%twinFraction(1_pInt:prm%totalNtwin,of)) &
             - sum(stt%stressTransFraction(1_pInt:prm%totalNtrans,of)) &
             - sum(stt%strainTransFraction(1_pInt:prm%totalNtrans,of))

 Lp = 0.0_pReal
 dLp_dMp = 0.0_pReal 

 call kinetics_slip(prm,stt,mse,of,Mp,temperature,gdot_slip,dgdot_dtau_slip)
 slipContribution: do i = 1_pInt, prm%totalNslip
   Lp = Lp + gdot_slip(i)*prm%Schmid_slip(1:3,1:3,i)
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + dgdot_dtau_slip(i) * prm%Schmid_slip(k,l,i) * prm%Schmid_slip(m,n,i)
 enddo slipContribution
 
 !ToDo: Why do this before shear banding?
 Lp      = Lp      * f_unrotated
 dLp_dMp = dLp_dMp * f_unrotated
 
 shearBandingContribution: if(dNeq0(prm%sbVelocity)) then

   BoltzmannRatio = prm%sbQedge/(kB*Temperature)
   call math_eigenValuesVectorsSym(Mp,eigValues,eigVectors,error)

   do i = 1_pInt,6_pInt
     sb_s = 0.5_pReal*sqrt(2.0_pReal)*math_mul33x3(eigVectors,sb_sComposition(1:3,i))
     sb_m = 0.5_pReal*sqrt(2.0_pReal)*math_mul33x3(eigVectors,sb_mComposition(1:3,i))
     Schmid_shearBand = math_tensorproduct33(sb_s,sb_m)
     tau = math_mul33xx33(Mp,Schmid_shearBand)
   
     significantShearBandStress: if (abs(tau) > tol_math_check) then
       StressRatio_p       = (abs(tau)/prm%sbResistance)**prm%pShearBand
       gdot_sb = sign(prm%sbVelocity*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**prm%qShearBand), tau)
       dgdot_dtau = ((abs(gdot_sb)*BoltzmannRatio* prm%pShearBand*prm%qShearBand)/ prm%sbResistance) &
                  * (abs(tau)/prm%sbResistance)**(prm%pShearBand-1.0_pReal) &
                  * (1.0_pReal-StressRatio_p)**(prm%qShearBand-1.0_pReal)
 
       Lp = Lp + gdot_sb * Schmid_shearBand
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                          + dgdot_dtau * Schmid_shearBand(k,l) * Schmid_shearBand(m,n)
     endif significantShearBandStress
   enddo

 endif shearBandingContribution
 
 call kinetics_twin(prm,stt,mse,of,Mp,temperature,gdot_slip,gdot_twin,dgdot_dtau_twin)
 gdot_twin       = f_unrotated * gdot_twin
 dgdot_dtau_twin = f_unrotated * dgdot_dtau_twin
 twinContibution: do i = 1_pInt, prm%totalNtwin
   Lp = Lp + gdot_twin(i)*prm%Schmid_twin(1:3,1:3,i)
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + dgdot_dtau_twin(i)* prm%Schmid_twin(k,l,i)*prm%Schmid_twin(m,n,i)
 enddo twinContibution

 transConstribution: do i = 1_pInt, prm%totalNtrans

   tau = math_mul33xx33(Mp,prm%Schmid_trans(1:3,1:3,i))

   significantTransStress: if (tau > tol_math_check) then
     StressRatio_s = (mse%threshold_stress_trans(i,of)/tau)**prm%s(i)

     isFCCtrans: if (prm%fccTwinTransNucleation) then
       s1=prm%fcc_twinNucleationSlipPair(1,i)
       s2=prm%fcc_twinNucleationSlipPair(2,i)
       if (tau < mse%tau_r_trans(i,of)) then
         Ndot0_trans=(abs(gdot_slip(s1))*(stt%rhoEdge(s2,of)+stt%rhoEdgeDip(s2,of))+&  !!!!! correct?
                           abs(gdot_slip(s2))*(stt%rhoEdge(s1,of)+stt%rhoEdgeDip(s1,of)))/&
                          (prm%L0_trans*prm%burgers_slip(i))*&
                          (1.0_pReal-exp(-prm%VcrossSlip/(kB*Temperature)*(mse%tau_r_trans(i,of)-tau)))
       else
         Ndot0_trans=0.0_pReal
       end if
     else isFCCtrans
       Ndot0_trans=prm%Ndot0_trans(i)
     endif isFCCtrans

     gdot_trans      = mse%martensiteVolume(i,of) * Ndot0_trans*exp(-StressRatio_s)
     gdot_trans      = f_unrotated * gdot_trans
     dgdot_dtau = ((gdot_trans*prm%s(i))/tau)*StressRatio_s
     Lp = Lp + gdot_trans*prm%Schmid_trans(1:3,1:3,i)
     
     forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
       dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                        + dgdot_dtau * prm%Schmid_trans(k,l,i)* prm%Schmid_trans(m,n,i)
   endif significantTransStress

 enddo transConstribution

 end associate
 
end subroutine plastic_dislotwin_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_dotState(Mp,Temperature,instance,of)
 use prec, only: &
   tol_math_check, &
   dEq0
 use math, only: &
   math_mul33xx33, &
   math_Mandel6to33, &
   pi
 use material, only: &
   plasticState

 implicit none
 real(pReal), dimension(3,3),  intent(in):: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                  intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),                intent(in) :: &
   instance, &
   of

 integer(pInt) :: i,s1,s2
 real(pReal) :: f_unrotated,StressRatio_p,BoltzmannRatio,&
             EdgeDipMinDistance,AtomicVolume,VacancyDiffusion,StressRatio_r,Ndot0_twin,stressRatio,&
             Ndot0_trans,StressRatio_s,EdgeDipDistance, ClimbVelocity,DotRhoEdgeDipClimb,DotRhoEdgeDipAnnihilation, &
             DotRhoDipFormation,DotRhoMultiplication,DotRhoEdgeEdgeAnnihilation, &
            tau
 real(pReal), dimension(plasticState(instance)%Nslip) :: &
 gdot_slip

 associate(prm => param(instance),    stt => state(instance), &
           dot => dotstate(instance), mse => microstructure(instance))

 dot%whole(:,of) = 0.0_pReal

 f_unrotated = 1.0_pReal &
             - sum(stt%twinFraction(1_pInt:prm%totalNtwin,of)) &
             - sum(stt%stressTransFraction(1_pInt:prm%totalNtrans,of)) &
             - sum(stt%strainTransFraction(1_pInt:prm%totalNtrans,of))

 call kinetics_slip(prm,stt,mse,of,Mp,temperature,gdot_slip)
 slipState: do i = 1_pInt, prm%totalNslip
   tau = math_mul33xx33(Mp,prm%Schmid_slip(1:3,1:3,i))

   DotRhoMultiplication = abs(gdot_slip(i))/(prm%burgers_slip(i)*mse%mfp_slip(i,of))
   EdgeDipMinDistance = prm%CEdgeDipMinDistance*prm%burgers_slip(i)

   significantSlipStress2: if (dEq0(tau)) then
     DotRhoDipFormation = 0.0_pReal
   else significantSlipStress2
     EdgeDipDistance = (3.0_pReal*prm%mu*prm%burgers_slip(i))/(16.0_pReal*PI*abs(tau))
     if (EdgeDipDistance>mse%mfp_slip(i,of)) EdgeDipDistance = mse%mfp_slip(i,of)
     if (EdgeDipDistance<EdgeDipMinDistance) EdgeDipDistance = EdgeDipMinDistance
     if (prm%dipoleFormation) then
       DotRhoDipFormation = ((2.0_pReal*(EdgeDipDistance-EdgeDipMinDistance))/prm%burgers_slip(i)) &
                          * stt%rhoEdge(i,of)*abs(gdot_slip(i))
     else
       DotRhoDipFormation = 0.0_pReal
     endif
   endif significantSlipStress2
 
   !* Spontaneous annihilation of 2 single edge dislocations
   DotRhoEdgeEdgeAnnihilation = ((2.0_pReal*EdgeDipMinDistance)/prm%burgers_slip(i))*&
                                stt%rhoEdge(i,of)*abs(gdot_slip(i))
   !* Spontaneous annihilation of a single edge dislocation with a dipole constituent
   DotRhoEdgeDipAnnihilation = ((2.0_pReal*EdgeDipMinDistance)/prm%burgers_slip(i)) &
                             * stt%rhoEdgeDip(i,of)*abs(gdot_slip(i))
 
   !* Dislocation dipole climb
   AtomicVolume     = prm%CAtomicVolume*prm%burgers_slip(i)**(3.0_pReal) ! no need to calculate this over and over again
   VacancyDiffusion = prm%D0*exp(-prm%Qsd/(kB*Temperature))

   if (dEq0(tau)) then
     DotRhoEdgeDipClimb = 0.0_pReal
   else
     if (dEq0(EdgeDipDistance-EdgeDipMinDistance)) then
       DotRhoEdgeDipClimb = 0.0_pReal
     else
       ClimbVelocity = 3.0_pReal*prm%mu*VacancyDiffusion*AtomicVolume/ &
                       (2.0_pReal*pi*kB*Temperature*(EdgeDipDistance+EdgeDipMinDistance))
       DotRhoEdgeDipClimb = 4.0_pReal*ClimbVelocity*stt%rhoEdgeDip(i,of)/ &
                       (EdgeDipDistance-EdgeDipMinDistance)
     endif
   endif
   dot%rhoEdge(i,of)    = DotRhoMultiplication-DotRhoDipFormation-DotRhoEdgeEdgeAnnihilation
   dot%rhoEdgeDip(i,of) = DotRhoDipFormation-DotRhoEdgeDipAnnihilation-DotRhoEdgeDipClimb
   dot%accshear_slip(i,of) = abs(gdot_slip(i))
 enddo slipState
 
 twinState: do i = 1_pInt, prm%totalNtwin
   
   tau = math_mul33xx33(Mp,prm%Schmid_slip(1:3,1:3,i))
   
   significantTwinStress: if (tau > tol_math_check) then
     StressRatio_r = (mse%threshold_stress_twin(i,of)/tau)**prm%r(i)
     isFCCtwin: if (prm%fccTwinTransNucleation) then
       s1=prm%fcc_twinNucleationSlipPair(1,i)
       s2=prm%fcc_twinNucleationSlipPair(2,i)
       if (tau < mse%tau_r_twin(i,of)) then
         Ndot0_twin=(abs(gdot_slip(s1))*(stt%rhoEdge(s2,of)+stt%rhoEdgeDip(s2,of))+&
                     abs(gdot_slip(s2))*(stt%rhoEdge(s1,of)+stt%rhoEdgeDip(s1,of)))/&
                 (prm%L0_twin*prm%burgers_slip(i))*(1.0_pReal-exp(-prm%VcrossSlip/(kB*Temperature)*&
                 (mse%tau_r_twin(i,of)-tau)))
       else
         Ndot0_twin=0.0_pReal
       end if
     else isFCCtwin
       Ndot0_twin=prm%Ndot0_twin(i)
     endif isFCCtwin
     dot%twinFraction(i,of) = f_unrotated * mse%twinVolume(i,of)*Ndot0_twin*exp(-StressRatio_r)
     dot%accshear_twin(i,of) = dot%twinFraction(i,of) * prm%shear_twin(i)
   endif significantTwinStress
 
 enddo twinState

 transState: do i = 1_pInt, prm%totalNtrans

   tau = math_mul33xx33(Mp,prm%Schmid_trans(1:3,1:3,i))
  
  significantTransStress: if (tau > tol_math_check) then
     StressRatio_s = (mse%threshold_stress_trans(i,of)/tau)**prm%s(i)
     isFCCtrans: if (prm%fccTwinTransNucleation) then
       s1=prm%fcc_twinNucleationSlipPair(1,i)
       s2=prm%fcc_twinNucleationSlipPair(2,i)
       if (tau < mse%tau_r_trans(i,of)) then
         Ndot0_trans=(abs(gdot_slip(s1))*(stt%rhoEdge(s2,of)+stt%rhoEdgeDip(s2,of))+&
                      abs(gdot_slip(s2))*(stt%rhoEdge(s1,of)+stt%rhoEdgeDip(s1,of)))/&
                          (prm%L0_trans*prm%burgers_slip(i))*(1.0_pReal-exp(-prm%VcrossSlip/(kB*Temperature)*&
                          (mse%tau_r_trans(i,of)-tau)))
       else
         Ndot0_trans=0.0_pReal
       end if
     else isFCCtrans
       Ndot0_trans=prm%Ndot0_trans(i)
     endif isFCCtrans
     dot%strainTransFraction(i,of) = f_unrotated * &
                             mse%martensiteVolume(i,of)*Ndot0_trans*exp(-StressRatio_s)
        !* Dotstate for accumulated shear due to transformation
        !dot%accshear_trans(i,of) = dot%strainTransFraction(i,of) * &
        !                                                  lattice_sheartrans(index_myfamily+i,ph)
   endif significantTransStress
 
 enddo transState

 end associate
end subroutine plastic_dislotwin_dotState


!--------------------------------------------------------------------------------------------------
!> @brief calculates shear rates on slip systems
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_slip(prm,stt,mse,of,Mp,temperature,gdot_slip,dgdot_dtau_slip)
 use prec, only: &
  tol_math_check, &
  dNeq0
 use math, only: &
   math_mul33xx33

 implicit none
 type(tParameters), intent(in) :: &
   prm
 type(tDislotwinState), intent(in) :: &
   stt
 integer(pInt),     intent(in) :: &
   of
 type(tDislotwinMicrostructure), intent(in) :: &
   mse
 real(pReal), dimension(prm%totalNslip), intent(out) :: &
   gdot_slip
 real(pReal), dimension(prm%totalNslip), optional, intent(out) :: &
   dgdot_dtau_slip
 real(pReal), dimension(prm%totalNslip) :: &
   dgdot_dtau
 real(pReal), dimension(3,3), intent(in) :: &
   Mp
 real(pReal), intent(in) :: &
   temperature

 real, dimension(prm%totalNslip) :: &
   tau, &
   stressRatio, &
   StressRatio_p, &
   BoltzmannRatio, &
   v_wait_inverse, &                                  !< inverse of the effective velocity of a dislocation waiting at obstacles (unsigned)
   v_run_inverse, &                                   !< inverse of the velocity of a free moving dislocation (unsigned)
   dV_wait_inverse_dTau, &
   dV_run_inverse_dTau, &
   dV_dTau, &
   tau_eff                                            !< effective resolved stress
 integer(pInt) :: i 

 do i = 1_pInt, prm%totalNslip
   tau(i) = math_mul33xx33(Mp,prm%Schmid_slip(1:3,1:3,i))
 enddo
 
 tau_eff = abs(tau)-mse%threshold_stress_slip(:,of)
   
 significantStress: where(tau_eff > tol_math_check)
   stressRatio    = tau_eff/(prm%SolidSolutionStrength+prm%tau_peierls)
   StressRatio_p  = stressRatio** prm%p
   BoltzmannRatio = prm%Qedge/(kB*Temperature)
   v_wait_inverse = prm%v0**(-1.0_pReal) * exp(BoltzmannRatio*(1.0_pReal-StressRatio_p)** prm%q)
   v_run_inverse  = prm%B/(tau_eff*prm%burgers_slip)

   gdot_slip = sign(stt%rhoEdge(:,of)*prm%burgers_slip/(v_wait_inverse+v_run_inverse),tau)

   dV_wait_inverse_dTau = v_wait_inverse * prm%p * prm%q * BoltzmannRatio &
                        * (stressRatio**(prm%p-1.0_pReal)) &
                        * (1.0_pReal-StressRatio_p)**(prm%q-1.0_pReal) &
                        / (prm%SolidSolutionStrength+prm%tau_peierls)
   dV_run_inverse_dTau  = v_run_inverse/tau_eff
   dV_dTau              = (dV_wait_inverse_dTau+dV_run_inverse_dTau) &
                        / (v_wait_inverse+v_run_inverse)**2.0_pReal
   dgdot_dtau = dV_dTau*stt%rhoEdge(:,of)*prm%burgers_slip
 else where significantStress
   gdot_slip  = 0.0_pReal
   dgdot_dtau = 0.0_pReal
 end where significantStress
 
 if(present(dgdot_dtau_slip)) dgdot_dtau_slip = dgdot_dtau

end subroutine kinetics_slip


!--------------------------------------------------------------------------------------------------
!> @brief calculates shear rates on twin systems
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_twin(prm,stt,mse,of,Mp,temperature,gdot_slip,gdot_twin,dgdot_dtau_twin)
 use prec, only: &
   tol_math_check, &
   dNeq0
 use math, only: &
   math_mul33xx33

 implicit none
 type(tParameters), intent(in) :: &
   prm
 type(tDislotwinState), intent(in) :: &
   stt
 integer(pInt),     intent(in) :: &
   of
 type(tDislotwinMicrostructure), intent(in) :: &
   mse
 real(pReal), dimension(prm%totalNslip), intent(out) :: &
   gdot_slip
 real(pReal), dimension(prm%totalNtwin), intent(out) :: &
   gdot_twin
 real(pReal), dimension(prm%totalNtwin), optional, intent(out) :: &
   dgdot_dtau_twin
 real(pReal), dimension(3,3), intent(in) :: &
   Mp
 real(pReal), intent(in) :: &
   temperature

 real, dimension(prm%totalNtwin) :: &
   tau, &
   Ndot0_twin, &
   stressRatio_r, &
   dgdot_dtau

 integer(pInt) :: i,s1,s2

 do i = 1_pInt, prm%totalNtwin
   tau(i) = math_mul33xx33(Mp,prm%Schmid_twin(1:3,1:3,i))
   isFCC: if (prm%fccTwinTransNucleation) then
     s1=prm%fcc_twinNucleationSlipPair(1,i)
     s2=prm%fcc_twinNucleationSlipPair(2,i)
     if (tau(i) < mse%tau_r_twin(i,of)) then
       Ndot0_twin=(abs(gdot_slip(s1))*(stt%rhoEdge(s2,of)+stt%rhoEdgeDip(s2,of))+&
                   abs(gdot_slip(s2))*(stt%rhoEdge(s1,of)+stt%rhoEdgeDip(s1,of)))/&
               (prm%L0_twin*prm%burgers_slip(i))*&
               (1.0_pReal-exp(-prm%VcrossSlip/(kB*Temperature)*&
               (mse%tau_r_twin(i,of)-tau)))
     else
       Ndot0_twin=0.0_pReal
     end if
   else isFCC
     Ndot0_twin=prm%Ndot0_twin(i)
   endif isFCC
 enddo


 significantStress: where(tau > tol_math_check)
   StressRatio_r = (mse%threshold_stress_twin(:,of)/tau)**prm%r
   gdot_twin = prm%shear_twin * mse%twinVolume(:,of) * Ndot0_twin*exp(-StressRatio_r)
   dgdot_dtau = ((gdot_twin*prm%r)/tau)*StressRatio_r
 else where significantStress
   gdot_twin  = 0.0_pReal
   dgdot_dtau = 0.0_pReal
 end where significantStress

 if(present(dgdot_dtau_twin)) dgdot_dtau_twin = dgdot_dtau

end subroutine kinetics_twin

 
!--------------------------------------------------------------------------------------------------
!> @brief calculates shear rates on transformation systems
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_trans(prm,stt,mse,of,Mp,temperature,gdot_slip,gdot_trans,dgdot_dtau_trans)
 use prec, only: &
   tol_math_check, &
   dNeq0
 use math, only: &
   math_mul33xx33

 implicit none
 type(tParameters), intent(in) :: &
   prm
 type(tDislotwinState), intent(in) :: &
   stt
 integer(pInt),     intent(in) :: &
   of
 type(tDislotwinMicrostructure), intent(in) :: &
   mse
 real(pReal), dimension(prm%totalNslip), intent(out) :: &
   gdot_slip
 real(pReal), dimension(prm%totalNtrans), intent(out) :: &
   gdot_trans
 real(pReal), dimension(prm%totalNtrans), optional, intent(out) :: &
   dgdot_dtau_trans
 real(pReal), dimension(3,3), intent(in) :: &
   Mp
 real(pReal), intent(in) :: &
   temperature

 real, dimension(prm%totalNtrans) :: &
   tau, &
   Ndot0_trans, &
   stressRatio_r, &
   dgdot_dtau

 integer(pInt) :: i,s1,s2

 do i = 1_pInt, prm%totalNtrans
   tau(i) = math_mul33xx33(Mp,prm%Schmid_trans(1:3,1:3,i))
   isFCC: if (prm%fccTwinTransNucleation) then
     s1=prm%fcc_twinNucleationSlipPair(1,i)
     s2=prm%fcc_twinNucleationSlipPair(2,i)
     if (tau(i) < mse%tau_r_trans(i,of)) then
       Ndot0_trans=(abs(gdot_slip(s1))*(stt%rhoEdge(s2,of)+stt%rhoEdgeDip(s2,of))+&
                    abs(gdot_slip(s2))*(stt%rhoEdge(s1,of)+stt%rhoEdgeDip(s1,of)))/&
               (prm%L0_trans*prm%burgers_slip(i))*&                               ! burgers_slip correct?
               (1.0_pReal-exp(-prm%VcrossSlip/(kB*Temperature)*&
               (mse%tau_r_trans(i,of)-tau)))
     else
       Ndot0_trans=0.0_pReal
     end if
   else isFCC
     Ndot0_trans=prm%Ndot0_trans(i)
   endif isFCC
 enddo
!
!  
!     endif isFCCtrans
!     dot%strainTransFraction(i,of) = f_unrotated * &
!                             mse%martensiteVolume(i,of)*Ndot0_trans*exp(-StressRatio_s)
!        !* Dotstate for accumulated shear due to transformation
!        !dot%accshear_trans(i,of) = dot%strainTransFraction(i,of) * &
!        !                                                  lattice_sheartrans(index_myfamily+i,ph)
!   endif significantTransStress
! 
! enddo transState
!
!
! significantStress: where(tau > tol_math_check)
!   StressRatio_r = (mse%threshold_stress_twin(:,of)/tau)**prm%r
!   gdot_twin = prm%shear_twin * mse%twinVolume(:,of) * Ndot0_twin*exp(-StressRatio_r)
!   dgdot_dtau = ((gdot_twin*prm%r)/tau)*StressRatio_r
! else where significantStress
!   gdot_twin  = 0.0_pReal
!   dgdot_dtau = 0.0_pReal
! end where significantStress
!
! if(present(dgdot_dtau_twin)) dgdot_dtau_twin = dgdot_dtau
!
end subroutine kinetics_trans

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_dislotwin_postResults(Mp,Temperature,instance,of) result(postResults)
 use prec, only: &
   tol_math_check, &
   dEq0
 use math, only: &
   PI, &
   math_mul33xx33, &
   math_Mandel6to33

 implicit none
 real(pReal), dimension(3,3),intent(in) :: &
   Mp                                                                                               !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
   instance, &
   of

 real(pReal), dimension(sum(plastic_dislotwin_sizePostResult(:,instance))) :: &
   postResults

 integer(pInt) :: &
   o,c,j,&
   s1,s2
 real(pReal) :: sumf_twin,tau,StressRatio_p,StressRatio_pminus1,BoltzmannRatio,DotGamma0,StressRatio_r,Ndot0_twin,dgdot_dtauslip, &
   stressRatio
 real(pReal), dimension(param(instance)%totalNslip) :: &
   gdot_slip
 
 type(tParameters) :: prm
 type(tDislotwinState) :: stt
 type(tDislotwinMicrostructure) :: mse
 

 associate(prm => param(instance), stt => state(instance), mse => microstructure(instance))

 sumf_twin = sum(stt%twinFraction(1_pInt:prm%totalNtwin,of))
 
 c = 0_pInt
 postResults = 0.0_pReal
 do o = 1_pInt,size(prm%outputID)
   select case(prm%outputID(o))
 
     case (edge_density_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%rhoEdge(1_pInt:prm%totalNslip,of)
       c = c + prm%totalNslip
     case (dipole_density_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%rhoEdgeDip(1_pInt:prm%totalNslip,of)
       c = c + prm%totalNslip
     case (shear_rate_slip_ID)
       call kinetics_slip(prm,stt,mse,of,Mp,temperature,postResults(c+1:c+prm%totalNslip))
       c = c + prm%totalNslip
     case (accumulated_shear_slip_ID)
      postResults(c+1_pInt:c+prm%totalNslip)  = stt%accshear_slip(1_pInt:prm%totalNslip,of)
       c = c + prm%totalNslip
     case (mfp_slip_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = mse%mfp_slip(1_pInt:prm%totalNslip,of)
       c = c + prm%totalNslip
     case (resolved_stress_slip_ID)
       do j = 1_pInt, prm%totalNslip
         postResults(c+j) = math_mul33xx33(Mp,prm%Schmid_slip(1:3,1:3,j))
       enddo
       c = c + prm%totalNslip
     case (threshold_stress_slip_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = mse%threshold_stress_slip(1_pInt:prm%totalNslip,of)
       c = c + prm%totalNslip
     case (edge_dipole_distance_ID)
       do j = 1_pInt, prm%totalNslip
         postResults(c+j) = (3.0_pReal*prm%mu*prm%burgers_slip(j)) &
                          / (16.0_pReal*PI*abs(math_mul33xx33(Mp,prm%Schmid_slip(1:3,1:3,j))))
         postResults(c+j)=min(postResults(c+j),mse%mfp_slip(j,of))
 !       postResults(c+j)=max(postResults(c+j),&
 !                                                       plasticState(ph)%state(4*ns+2*nt+2*nr+j, of))
       enddo
       c = c + prm%totalNslip
 !     case (resolved_stress_shearband_ID)
 !       do j = 1_pInt,6_pInt                                                                       ! loop over all shearband families
 !          postResults(c+j) = dot_product(Tstar_v,sbSv(1:6,j,ipc,ip,el))
 !       enddo
 !       c = c + 6_pInt
 !     case (shear_rate_shearband_ID)
 !       do j = 1_pInt,6_pInt                                                                       ! loop over all shearbands
 !         tau = dot_product(Tstar_v,sbSv(1:6,j,ipc,ip,el))
 !         if (abs(tau) < tol_math_check) then
 !           StressRatio_p = 0.0_pReal
 !           StressRatio_pminus1 = 0.0_pReal
 !         else
 !           StressRatio_p = (abs(tau)/prm%sbResistance)**prm%pShearBand
 !           StressRatio_pminus1 = (abs(tau)/prm%sbResistance)**(prm%pShearBand-1.0_pReal)
 !         endif
 !         BoltzmannRatio = prm%sbQedge/(kB*Temperature)
 !         DotGamma0 = prm%sbVelocity
 !         postResults(c+j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**prm%qShearBand)*&
 !                            sign(1.0_pReal,tau)
 !       enddo 
 !      c = c + 6_pInt
     case (twin_fraction_ID)
       postResults(c+1_pInt:c+prm%totalNtwin) = stt%twinFraction(1_pInt:prm%totalNtwin,of)
       c = c + prm%totalNtwin
     case (shear_rate_twin_ID)
       do j = 1_pInt, prm%totalNslip
         tau = math_mul33xx33(Mp,prm%Schmid_slip(1:3,1:3,j))
         if((abs(tau)-mse%threshold_stress_slip(j,of)) > tol_math_check) then
           StressRatio_p = ((abs(tau)-mse%threshold_stress_slip(j,of))/&
                           (prm%SolidSolutionStrength+&
                            prm%tau_peierls(j)))&
                                              **prm%p(j)
           StressRatio_pminus1 = ((abs(tau)-mse%threshold_stress_slip(j,of))/&
                                 (prm%SolidSolutionStrength+&
                                  prm%tau_peierls(j)))&
                                              **(prm%p(j)-1.0_pReal)
           BoltzmannRatio = prm%Qedge(j)/(kB*Temperature)
           DotGamma0 =  stt%rhoEdge(j,of)*prm%burgers_slip(j)* prm%v0(j)
 
           gdot_slip(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**&
                          prm%q(j))*sign(1.0_pReal,tau)
         else
           gdot_slip(j) = 0.0_pReal
         endif
       enddo

       do j = 1_pInt, prm%totalNtwin
         tau = math_mul33xx33(Mp,prm%Schmid_twin(1:3,1:3,j))

         if ( tau > 0.0_pReal ) then
           isFCCtwin: if (prm%fccTwinTransNucleation) then
             s1=prm%fcc_twinNucleationSlipPair(1,j)
             s2=prm%fcc_twinNucleationSlipPair(2,j)
             if (tau < mse%tau_r_twin(j,of)) then
                 Ndot0_twin=(abs(gdot_slip(s1))*(stt%rhoEdge(s2,of)+stt%rhoEdgeDip(s2,of))+&
                        abs(gdot_slip(s2))*(stt%rhoEdge(s1,of)+stt%rhoEdgeDip(s1,of)))/&
                       (prm%L0_twin* prm%burgers_slip(j))*&
                       (1.0_pReal-exp(-prm%VcrossSlip/(kB*Temperature)* (mse%tau_r_twin(j,of)-tau)))
              else
               Ndot0_twin=0.0_pReal
              end if
            else isFCCtwin
              Ndot0_twin=prm%Ndot0_twin(j)
            endif isFCCtwin
            StressRatio_r = (mse%threshold_stress_twin(j,of)/tau) **prm%r(j)
            postResults(c+j) = (prm%MaxTwinFraction-sumf_twin)*prm%shear_twin(j) &
                             * mse%twinVolume(j,of)*Ndot0_twin*exp(-StressRatio_r)
         endif
       enddo
       c = c + prm%totalNtwin
     case (accumulated_shear_twin_ID)
       postResults(c+1_pInt:c+prm%totalNtwin) = stt%accshear_twin(1_pInt:prm%totalNtwin,of)
       c = c + prm%totalNtwin     
     case (mfp_twin_ID)
       postResults(c+1_pInt:c+prm%totalNtwin) = mse%mfp_twin(1_pInt:prm%totalNtwin,of)
       c = c + prm%totalNtwin
     case (resolved_stress_twin_ID)
       do j = 1_pInt, prm%totalNtwin
         postResults(c+j) = math_mul33xx33(Mp,prm%Schmid_twin(1:3,1:3,j))
       enddo
       c = c + prm%totalNtwin
     case (threshold_stress_twin_ID)
       postResults(c+1_pInt:c+prm%totalNtwin) = mse%threshold_stress_twin(1_pInt:prm%totalNtwin,of)
       c = c + prm%totalNtwin
     case (stress_exponent_ID)
       do j = 1_pInt, prm%totalNslip
         tau = math_mul33xx33(Mp,prm%Schmid_slip(1:3,1:3,j))
         if((abs(tau)-mse%threshold_stress_slip(j,of)) > tol_math_check) then
           StressRatio_p = ((abs(tau)-mse%threshold_stress_slip(j,of))/&
                           (prm%SolidSolutionStrength+&
                            prm%tau_peierls(j)))&
                                              **prm%p(j)
           StressRatio_pminus1 = ((abs(tau)-mse%threshold_stress_slip(j,of))/&
                                 (prm%SolidSolutionStrength+&
                                  prm%tau_peierls(j)))&
                                              **(prm%p(j)-1.0_pReal)
           BoltzmannRatio = prm%Qedge(j)/(kB*Temperature)
           DotGamma0 = stt%rhoEdge(j,of)*prm%burgers_slip(j)* prm%v0(j)

           gdot_slip(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**&
                        prm%q(j))*sign(1.0_pReal,tau)

           dgdot_dtauslip = abs(gdot_slip(j))*BoltzmannRatio*prm%p(j) *prm%q(j)/&
             (prm%SolidSolutionStrength+ prm%tau_peierls(j))*&
             StressRatio_pminus1*(1-StressRatio_p)**(prm%q(j)-1.0_pReal)
         else
           gdot_slip(j) = 0.0_pReal
           dgdot_dtauslip = 0.0_pReal
         endif
         postResults(c+j) = merge(0.0_pReal,(tau/gdot_slip(j))*dgdot_dtauslip,dEq0(gdot_slip(j)))
       enddo
       c = c + prm%totalNslip
     case (stress_trans_fraction_ID)
       postResults(c+1_pInt:c+prm%totalNtrans) = stt%stressTransFraction(1_pInt:prm%totalNtrans,of)
       c = c + prm%totalNtrans
     case (strain_trans_fraction_ID)
       postResults(c+1_pInt:c+prm%totalNtrans) = stt%strainTransFraction(1_pInt:prm%totalNtrans,of)
       c = c + prm%totalNtrans
   end select
 enddo
 end associate
end function plastic_dislotwin_postResults

end module plastic_dislotwin
