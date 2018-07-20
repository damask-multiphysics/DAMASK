!--------------------------------------------------------------------------------------------------
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
 real(pReal),                         dimension(:,:,:,:),     allocatable,         private :: &
   Ctwin66,&                                                                   !< twin elasticity matrix in Mandel notation for each instance
   Ctrans66                                                                   !< trans elasticity matrix in Mandel notation for each instance
 real(pReal),                         dimension(:,:),         allocatable,         private :: &
   tau_r_twin, &                                                             !< stress to bring partial close together for each twin system and instance
   tau_r_trans                                                            !< stress to bring partial close together for each trans system and instance
 real(pReal),                         dimension(:,:,:),       allocatable,         private :: &
   forestProjectionEdge, &                                                   !< matrix of forest projections of edge dislocations for each instance
   projectionMatrix_Trans                                                    !< matrix for projection of slip system shear on fault band (twin) systems for each instance
 real(pReal),                         dimension(:,:,:,:,:),   allocatable,         private :: &
   sbSv

 enum, bind(c) 
   enumerator :: undefined_ID, &
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
                 sb_eigenvalues_ID, &
                 sb_eigenvectors_ID, &
                 stress_trans_fraction_ID, &
                 strain_trans_fraction_ID, &
                 trans_fraction_ID
 end enum
 
  type,private :: tParameters
   integer(kind(undefined_ID)),         dimension(:),         allocatable,          private :: &
     outputID                                                                                       !< ID of each post result output
   
   real(pReal) :: &
     CAtomicVolume, &                                                                               !< atomic volume in Bugers vector unit
     D0, &                                                                                          !< prefactor for self-diffusion coefficient
     Qsd, &                                                                                         !< activation energy for dislocation climb
     GrainSize, &                                                                                  !<grain size
     pShearBand, &                                                                                 !< p-exponent in shear band velocity
     qShearBand, &                                                                                 !< q-exponent in shear band velocity
     MaxTwinFraction, &                                                                            !<max allowed total twin volume fraction
     CEdgeDipMinDistance, &                                                                        !<
     Cmfptwin, &                                                                                   !<
     Cthresholdtwin, &                                                                             !<
     SolidSolutionStrength, &                                                                      !<strength due to elements in solid solution
     L0_twin, &                                                                                    !< Length of twin nuclei in Burgers vectors
     L0_trans, &                                                                                   !< Length of trans nuclei in Burgers vectors
     xc_twin, &                                                                                    !< critical distance for formation of twin nucleus
     xc_trans, &                                                                                   !< critical distance for formation of trans nucleus
     VcrossSlip, &                                                                                 !< cross slip volume
     sbResistance, &                                                                               !< value for shearband resistance (might become an internal state variable at some point)
     sbVelocity, &                                                                                 !< value for shearband velocity_0
     sbQedge, &                                                                                    !< value for shearband systems Qedge
     SFE_0K, &                                                                                     !< stacking fault energy at zero K
     dSFE_dT, &                                                                                    !< temperature dependance of stacking fault energy
     dipoleFormationFactor, &                                                                      !< scaling factor for dipole formation: 0: off, 1: on. other values not useful
     aTolRho, &                                                                                    !< absolute tolerance for integration of dislocation density
     aTolTwinFrac, &                                                                               !< absolute tolerance for integration of twin volume fraction
     aTolTransFrac, &                                                                              !< absolute tolerance for integration of trans volume fraction
     deltaG, &                                                                                     !< Free energy difference between austensite and martensite
     Cmfptrans, &                                                                                  !<
     Cthresholdtrans, &                                                                            !<
     transStackHeight                                                                              !< Stack height of hex nucleus 
   integer(pInt),                  private :: & 
     totalNslip, &                                                                                          !< number of active slip systems for each family and instance
     totalNtwin, &                                                                                          !< number of active twin systems for each family and instance
    totalNtrans                                                                                            !< number of active transformation systems for each family and instance 
   integer(pInt),                  dimension(:),            allocatable,           private :: & 
     Nslip, &                                                                                          !< number of active slip systems for each family and instance
     Ntwin, &                                                                                          !< number of active twin systems for each family and instance
     Ntrans                                                                                            !< number of active transformation systems for each family and instance
   real(pReal),                  dimension(:),            allocatable,           private :: & 
     rho0, & !< initial unipolar dislocation density per slip system
     rhoDip0, & !< initial dipole dislocation density per slip system
     burgers_slip, &                                                   !< absolute length of burgers vector [m] for each slip systems
     burgers_twin, &                                                   !< absolute length of burgers vector [m] for each slip systems
     burgers_trans, &                                                  !< absolute length of burgers vector [m] for each twin family and instance
     Qedge,&                                                !< activation energy for glide [J] for each slip system and instance
     v0, &                                                    !dislocation velocity prefactor [m/s] for each slip system and instance
     tau_peierlsPerSlipFamily,&                                            !< Peierls stress [Pa] for each family and instance
     Ndot0_twin, &                                                 !< twin nucleation rate [1/m³s] for each twin system and instance
     Ndot0_trans, &                                                !< trans nucleation rate [1/m³s] for each trans system and instance
     twinsize, &                                                  !< twin thickness [m] for each twin system and instance
     CLambdaSlipPerSlipSystem, &                                               !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
     lamellarsizePerTransSystem, &                                             !< martensite lamellar thickness [m] for each trans system and instance
     p, &                                                         !< p-exponent in glide velocity
     q, &                                                         !< q-exponent in glide velocity
     r, &                                                         !< r-exponent in twin nucleation rate
     s                                                           !< s-exponent in trans nucleation rate
   real(pReal),                  dimension(:,:),            allocatable,           private :: & 
   interaction_SlipSlip, &                                                   !< coefficients for slip-slip interaction for each interaction type and instance
   interaction_SlipTwin, &                                                   !< coefficients for slip-twin interaction for each interaction type and instance
   interaction_TwinSlip, &                                                   !< coefficients for twin-slip interaction for each interaction type and instance
   interaction_TwinTwin, &                                                   !< coefficients for twin-twin interaction for each interaction type and instance
   interaction_SlipTrans, &                                                  !< coefficients for slip-trans interaction for each interaction type and instance
   interaction_TransSlip, &                                                  !< coefficients for trans-slip interaction for each interaction type and instance
   interaction_TransTrans                                                   !< coefficients for trans-trans interaction for each interaction type and instance
  end type 
  
  type(tParameters), dimension(:), allocatable, private,target :: param                                !< containers of constitutive parameters (len Ninstance)
 
 
 type, private :: tDislotwinState
  
   real(pReal), pointer,     dimension(:,:) :: &
     rhoEdge, &
     rhoEdgeDip, &
     accshear_slip, &
     twinFraction, &
     accshear_twin, &
     stressTransFraction, &
     strainTransFraction , &
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
     martensiteVolume
 end type

 type(tDislotwinState), allocatable, dimension(:), private :: &
   state, &
   dotState

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
subroutine plastic_dislotwin_init(fileUnit)
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   dEq0, &
   dNeq0, &
   dNeq
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use math, only: &
   math_Mandel3333to66, &
   math_Voigt66to3333, &
   math_mul3x3, &
   math_expand,&
   pi
 use mesh, only: &
   mesh_maxNips, &
   mesh_NcpElems
 use IO, only: &
   IO_read, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue, &
   IO_intValue, &
   IO_warning, &
   IO_error, &
   IO_timeStamp, &
   IO_EOF
 use material, only: &
   homogenization_maxNgrains, &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_DISLOTWIN_label, &
   PLASTICITY_DISLOTWIN_ID, &
   material_phase, &  
   plasticState
 use config, only: &
   MATERIAL_partPhase, &
   config_phase
 use lattice
 use numerics,only: &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt) :: maxNinstance,&
                  f,instance,j,i,k,l,m,n,o,p,q,r,s,p1, &
                  offset_slip, index_myFamily, index_otherFamily, &
                  startIndex, endIndex, outputSize
 integer(pInt) :: sizeState, sizeDotState, sizeDeltaState
 integer(pInt) :: NofMyPhase   
 integer(kind(undefined_ID)) outputID
 
  real(pReal),                         dimension(:,:,:,:,:), allocatable :: &
   Ctwin3333, &                                                                 !< twin elasticity matrix for each instance
   Ctrans3333                                                                 !< trans elasticity matrix for each instance
 
 real(pReal),  allocatable, dimension(:) :: &
     invLambdaSlip0,&
     MeanFreePathSlip0,&
     MeanFreePathTrans0,&
     MeanFreePathTwin0,&
     tauSlipThreshold0,&
     TwinVolume0,&
     MartensiteVolume0
     
      real(pReal),  allocatable, dimension(:,:) :: temp1,temp2,temp3

 character(len=65536) :: &
   tag  = ''
 
 character(len=65536), dimension(:), allocatable :: outputs
 integer(pInt), dimension(0), parameter :: emptyInt = [integer(pInt)::]
 real(pReal),   dimension(0), parameter :: emptyReal = [real(pReal)::]
 character(len=65536),   dimension(0), parameter :: emptyString = [character(len=65536)::]


 type(tParameters),pointer :: prm

  
 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_DISLOTWIN_label//' init  -+>>>'
 write(6,'(/,a)')   ' A. Ma and F. Roters, Acta Materialia, 52(12):3603–3612, 2004'
 write(6,'(a)')     ' https://doi.org/10.1016/j.actamat.2004.04.012'
 write(6,'(/,a)')   ' F.Roters et al., Computational Materials Science, 39:91–95, 2007'
 write(6,'(a)')     ' https://doi.org/10.1016/j.commatsci.2006.04.014'
 write(6,'(/,a)')   ' Wong et al., Acta Materialia, 118:140–151, 2016'
 write(6,'(a,/)')     ' https://doi.org/10.1016/j.actamat.2016.07.032'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_DISLOTWIN_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 

 allocate(plastic_dislotwin_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(plastic_dislotwin_output(maxval(phase_Noutput),maxNinstance))
          plastic_dislotwin_output = ''
 allocate(param(maxNinstance))
 allocate(state(maxNinstance))
 allocate(dotState(maxNinstance))
 allocate(sbSv(6,6,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), source=0.0_pReal)

 

 do p = 1_pInt, size(phase_plasticityInstance)
   if (phase_plasticity(p) /= PLASTICITY_DISLOTWIN_ID) cycle
   instance = phase_plasticityInstance(p)
   prm => param(instance)
     
   prm%Nslip =  config_phase(p)%getInts('nslip',defaultVal=emptyInt)
   if (size(prm%Nslip) > count(lattice_NslipSystem(:,p) > 0_pInt))       call IO_error(150_pInt,ext_msg='Nslip')
   if (any(lattice_NslipSystem(1:size(prm%Nslip),p)-prm%Nslip < 0_pInt)) call IO_error(150_pInt,ext_msg='Nslip')
   if (any(prm%Nslip < 0_pInt))                                          call IO_error(150_pInt,ext_msg='Nslip')
   prm%totalNslip = sum(prm%Nslip)

     if (prm%totalNslip > 0_pInt) then
       prm%rho0 = config_phase(p)%getFloats('rhoedge0')
       prm%rhoDip0 = config_phase(p)%getFloats('rhoedgedip0')

       prm%burgers_slip = config_phase(p)%getFloats('slipburgers')
       if (size(prm%burgers_slip) /= size(prm%Nslip)) call IO_error(150_pInt,ext_msg='slipburgers')
       prm%burgers_slip = math_expand(prm%burgers_slip,prm%Nslip)

       prm%Qedge = config_phase(p)%getFloats('qedge')
       prm%Qedge = math_expand(prm%Qedge,prm%Nslip)

       prm%v0 = config_phase(p)%getFloats('v0')
       prm%v0 = math_expand(prm%v0,prm%Nslip)

       prm%interaction_SlipSlip = spread(config_phase(p)%getFloats('interaction_slipslip'),2,1)     

       prm%CEdgeDipMinDistance      = config_phase(p)%getFloat('cedgedipmindistance')

       prm%CLambdaSlipPerSlipSystem = config_phase(p)%getFloats('clambdaslip')
       prm%CLambdaSlipPerSlipSystem= math_expand(prm%CLambdaSlipPerSlipSystem,prm%Nslip)

       prm%tau_peierlsPerSlipFamily = config_phase(p)%getFloats('tau_peierls',defaultVal=[0.0_pReal])

       prm%p = config_phase(p)%getFloats('p_slip')
       prm%q = config_phase(p)%getFloats('q_slip')
   endif

   prm%Ntwin =  config_phase(p)%getInts('ntwin', defaultVal=emptyInt)
   if (size(prm%Ntwin) > count(lattice_NtwinSystem(:,p) > 0_pInt))       call IO_error(150_pInt,ext_msg='Ntwin')
   if (any(lattice_NtwinSystem(1:size(prm%Ntwin),p)-prm%Ntwin < 0_pInt)) call IO_error(150_pInt,ext_msg='Ntwin')
   if (any(prm%Ntwin < 0_pInt))                                          call IO_error(150_pInt,ext_msg='Ntwin')
   prm%totalNtwin = sum(prm%Ntwin)
   
   if (prm%totalNtwin > 0_pInt) then
     prm%burgers_twin = config_phase(p)%getFloats('twinburgers')
     prm%burgers_twin = math_expand(prm%burgers_twin,prm%Ntwin)
     
     prm%xc_twin = config_phase(p)%getFloat('xc_twin')
     prm%Cthresholdtwin = config_phase(p)%getFloat('cthresholdtwin', defaultVal=0.0_pReal)
     prm%Cmfptwin       = config_phase(p)%getFloat('cmfptwin', defaultVal=0.0_pReal) ! ToDo: How to handle that???

     prm%interaction_TwinTwin = spread(config_phase(p)%getFloats('interaction_twintwin'),2,1)     
     if (lattice_structure(p) /= LATTICE_fcc_ID) then
       prm%Ndot0_twin = config_phase(p)%getFloats('ndot0_twin') 
       prm%Ndot0_twin = math_expand(prm%Ndot0_twin,prm%Ntwin)
     endif

     prm%twinsize = config_phase(p)%getFloats('twinsize')
     prm%twinsize= math_expand(prm%twinsize,prm%Ntwin)
     
     prm%r = config_phase(p)%getFloats('r_twin')
     prm%r = math_expand(prm%r,prm%Ntwin)
     
            
     prm%L0_twin = config_phase(p)%getFloat('l0_twin')
     

   endif
   
   prm%Ntrans            =  config_phase(p)%getInts('ntrans', defaultVal=emptyInt)
   prm%totalNtrans = sum(prm%Ntrans)
   !if (size > Nchunks_SlipFamilies + 1_pInt) call IO_error(150_pInt,ext_msg=extmsg)
   if (prm%totalNtrans > 0_pInt) then
     prm%burgers_trans = config_phase(p)%getFloats('transburgers')
     prm%burgers_trans = math_expand(prm%burgers_trans,prm%Ntrans)
     
     prm%Cthresholdtrans = config_phase(p)%getFloat('cthresholdtrans', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%transStackHeight = config_phase(p)%getFloat('transstackheight', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%Cmfptrans = config_phase(p)%getFloat('cmfptrans', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%deltaG = config_phase(p)%getFloat('deltag')
     prm%xc_trans = config_phase(p)%getFloat('xc_trans', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%L0_trans = config_phase(p)%getFloat('l0_trans')

     prm%interaction_TransTrans = spread(config_phase(p)%getFloats('interaction_transtrans'),2,1)     
     if (lattice_structure(p) /= LATTICE_fcc_ID) then
        prm%Ndot0_trans = config_phase(p)%getFloats('ndot0_trans')
        prm%Ndot0_trans = math_expand(prm%Ndot0_trans,prm%Ntrans)
     endif
     prm%lamellarsizePerTransSystem = config_phase(p)%getFloats('lamellarsize')
     prm%lamellarsizePerTransSystem = math_expand(prm%lamellarsizePerTransSystem,prm%Ntrans)
     prm%s = config_phase(p)%getFloats('s_trans',defaultVal=[0.0_pReal])
   endif
   
   if (sum(prm%Ntwin) > 0_pInt  .or. sum(prm%Ntrans) > 0_pInt) then
     prm%SFE_0K = config_phase(p)%getFloat('sfe_0k')
     prm%dSFE_dT = config_phase(p)%getFloat('dsfe_dt')
     prm%VcrossSlip = config_phase(p)%getFloat('vcrossslip')
   endif
   
   if (prm%totalNslip > 0_pInt .and. prm%totalNtwin > 0_pInt) then
    prm%interaction_SlipTwin = spread(config_phase(p)%getFloats('interaction_sliptwin'),2,1)     
    prm%interaction_TwinSlip = spread(config_phase(p)%getFloats('interaction_twinslip'),2,1)
    endif    

   if (prm%totalNslip > 0_pInt .and. prm%totalNtrans > 0_pInt) then
      prm%interaction_TransSlip = spread(config_phase(p)%getFloats('interaction_transslip'),2,1)     
      prm%interaction_SlipTrans = spread(config_phase(p)%getFloats('interaction_sliptrans'),2,1)     
   endif    


     prm%aTolRho       = config_phase(p)%getFloat('atol_rho', defaultVal=0.0_pReal)
     prm%aTolTwinFrac  = config_phase(p)%getFloat('atol_twinfrac', defaultVal=0.0_pReal)
     prm%aTolTransFrac = config_phase(p)%getFloat('atol_transfrac', defaultVal=0.0_pReal)

     prm%CAtomicVolume = config_phase(p)%getFloat('catomicvolume')
     prm%GrainSize =  config_phase(p)%getFloat('grainsize')
     prm%MaxTwinFraction = config_phase(p)%getFloat('maxtwinfraction') ! ToDo: only used in postResults

     prm%D0 = config_phase(p)%getFloat('d0')
     prm%Qsd = config_phase(p)%getFloat('qsd')
     prm%SolidSolutionStrength = config_phase(p)%getFloat('solidsolutionstrength')
     prm%dipoleFormationFactor= config_phase(p)%getFloat('dipoleformationfactor', defaultVal=1.0_pReal) ! ToDo: How to handle that???
     prm%sbVelocity   = config_phase(p)%getFloat('shearbandvelocity',defaultVal=0.0_pReal)
     if (prm%sbVelocity > 0.0_pReal) then  
       prm%sbResistance = config_phase(p)%getFloat('shearbandresistance')
       prm%sbQedge = config_phase(p)%getFloat('qedgepersbsystem')
       prm%pShearBand = config_phase(p)%getFloat('p_shearband')
       prm%qShearBand = config_phase(p)%getFloat('q_shearband')
     endif
 
     outputs = config_phase(p)%getStrings('(output)', defaultVal=emptyString)
     allocate(prm%outputID(0))
     do i= 1_pInt, size(outputs)
        outputID = undefined_ID
        select case(outputs(i))

           case ('edge_density')
             outputID = edge_density_ID
             outputSize = sum(prm%Nslip)

           case ('dipole_density')
             outputID = dipole_density_ID
             outputSize = sum(prm%Nslip)

           case ('shear_rate_slip','shearrate_slip')
             outputID = shear_rate_slip_ID
             outputSize = sum(prm%Nslip)

           case ('accumulated_shear_slip')
             outputID = accumulated_shear_slip_ID
             outputSize = sum(prm%Nslip)

           case ('mfp_slip')
             outputID = mfp_slip_ID
             outputSize = sum(prm%Nslip)

           case ('resolved_stress_slip')
             outputID = resolved_stress_slip_ID
             outputSize = sum(prm%Nslip)

           case ('threshold_stress_slip')
             outputID= threshold_stress_slip_ID
             outputSize = sum(prm%Nslip)

           case ('edge_dipole_distance')
             outputID = edge_dipole_distance_ID
             outputSize = sum(prm%Nslip)

           case ('stress_exponent')
             outputID = stress_exponent_ID
             outputSize = sum(prm%Nslip)

           case ('twin_fraction')
             outputID = twin_fraction_ID
             outputSize = sum(prm%Ntwin)

           case ('shear_rate_twin','shearrate_twin')
             outputID= shear_rate_twin_ID
             outputSize = sum(prm%Ntwin)

           case ('accumulated_shear_twin')
             outputID = accumulated_shear_twin_ID
             outputSize = sum(prm%Ntwin)
             
           case ('mfp_twin')
             outputID = mfp_twin_ID
             outputSize = sum(prm%Ntwin)
             
           case ('resolved_stress_twin')
             outputID = resolved_stress_twin_ID
             outputSize = sum(prm%Ntwin)
             
           case ('threshold_stress_twin')
             outputID = threshold_stress_twin_ID
             outputSize = sum(prm%Ntwin)
             
           case ('resolved_stress_shearband')
             outputID = resolved_stress_shearband_ID
             outputSize = 6_pInt
             
           case ('shear_rate_shearband','shearrate_shearband')
             outputID = shear_rate_shearband_ID
             outputSize = 6_pInt
             
           case ('sb_eigenvalues')
             outputID = sb_eigenvalues_ID
             outputSize = 3_pInt
             
           case ('sb_eigenvectors')
             outputID = sb_eigenvectors_ID
             outputSize = 3_pInt
             
           case ('stress_trans_fraction')
             outputID = stress_trans_fraction_ID
             outputSize = sum(prm%Ntrans)
             
           case ('strain_trans_fraction')
             outputID = strain_trans_fraction_ID
             outputSize = sum(prm%Ntrans)
             
           case ('trans_fraction','total_trans_fraction')
             outputID = trans_fraction_ID
             outputSize = sum(prm%Ntrans)
            
          end select
          
          if (outputID /= undefined_ID) then
          plastic_dislotwin_output(i,instance) = outputs(i)
          plastic_dislotwin_sizePostResult(i,instance) = outputSize
          prm%outputID = [prm%outputID , outputID]
        endif

 enddo

 
      do f = 1_pInt,lattice_maxNslipFamily
        !  if (rhoEdge0(f,instance) < 0.0_pReal) &
        !    call IO_error(211_pInt,el=instance,ext_msg='rhoEdge0 ('//PLASTICITY_DISLOTWIN_label//')')
        !  if (rhoEdgeDip0(f,instance) < 0.0_pReal) & 
         !   call IO_error(211_pInt,el=instance,ext_msg='rhoEdgeDip0 ('//PLASTICITY_DISLOTWIN_label//')')
        !  if (burgersPerSlipFamily(f,instance) <= 0.0_pReal) &
        !    call IO_error(211_pInt,el=instance,ext_msg='slipBurgers ('//PLASTICITY_DISLOTWIN_label//')')
          !if (v0PerSlipFamily(f,instance) <= 0.0_pReal) &
           ! call IO_error(211_pInt,el=instance,ext_msg='v0 ('//PLASTICITY_DISLOTWIN_label//')')
          !if (prm%tau_peierlsPerSlipFamily(f) < 0.0_pReal) &
           ! call IO_error(211_pInt,el=instance,ext_msg='tau_peierls ('//PLASTICITY_DISLOTWIN_label//')')
      enddo
      do f = 1_pInt,lattice_maxNtwinFamily
         ! if (burgersPerTwinFamily(f,instance) <= 0.0_pReal) &
         !   call IO_error(211_pInt,el=instance,ext_msg='twinburgers ('//PLASTICITY_DISLOTWIN_label//')')
          !if (Ndot0PerTwinFamily(f,instance) < 0.0_pReal) &
           ! call IO_error(211_pInt,el=instance,ext_msg='ndot0_twin ('//PLASTICITY_DISLOTWIN_label//')')
      enddo
      if (prm%CAtomicVolume <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='cAtomicVolume ('//PLASTICITY_DISLOTWIN_label//')')
      if (prm%D0 <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='D0 ('//PLASTICITY_DISLOTWIN_label//')')
      if (prm%Qsd <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='Qsd ('//PLASTICITY_DISLOTWIN_label//')')
      if (prm%totalNtwin > 0_pInt) then
        if (dEq0(prm%SFE_0K) .and. &
            dEq0(prm%dSFE_dT) .and. &
                            lattice_structure(p) == LATTICE_fcc_ID) &
          call IO_error(211_pInt,el=instance,ext_msg='SFE0K ('//PLASTICITY_DISLOTWIN_label//')')
        if (prm%aTolRho <= 0.0_pReal) &
          call IO_error(211_pInt,el=instance,ext_msg='aTolRho ('//PLASTICITY_DISLOTWIN_label//')')   
        if (prm%aTolTwinFrac <= 0.0_pReal) &
          call IO_error(211_pInt,el=instance,ext_msg='aTolTwinFrac ('//PLASTICITY_DISLOTWIN_label//')')
      endif
      if (prm%totalNtrans > 0_pInt) then
        if (dEq0(prm%SFE_0K) .and. &
            dEq0(prm%dSFE_dT) .and. &
                            lattice_structure(p) == LATTICE_fcc_ID) &
          call IO_error(211_pInt,el=instance,ext_msg='SFE0K ('//PLASTICITY_DISLOTWIN_label//')')
        if (prm%aTolTransFrac <= 0.0_pReal) &
          call IO_error(211_pInt,el=instance,ext_msg='aTolTransFrac ('//PLASTICITY_DISLOTWIN_label//')')
      endif
      !if (prm%sbResistance < 0.0_pReal) &
      !  call IO_error(211_pInt,el=instance,ext_msg='sbResistance ('//PLASTICITY_DISLOTWIN_label//')')
      !if (prm%sbVelocity < 0.0_pReal) &
      !  call IO_error(211_pInt,el=instance,ext_msg='sbVelocity ('//PLASTICITY_DISLOTWIN_label//')')
      !if (prm%sbVelocity > 0.0_pReal .and. &
      !    prm%pShearBand <= 0.0_pReal) &
      !  call IO_error(211_pInt,el=instance,ext_msg='pShearBand ('//PLASTICITY_DISLOTWIN_label//')')
      if (dNeq0(prm%dipoleFormationFactor) .and. &
          dNeq(prm%dipoleFormationFactor, 1.0_pReal)) &
        call IO_error(211_pInt,el=instance,ext_msg='dipoleFormationFactor ('//PLASTICITY_DISLOTWIN_label//')')
      if (prm%sbVelocity > 0.0_pReal .and. &
          prm%qShearBand <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='qShearBand ('//PLASTICITY_DISLOTWIN_label//')')

 enddo 
 
 ! ToDo: this should be stored somewhere else. Works only for the whole instance!!
! prm%totalNtwin should be the maximum over all totalNtwins!
 allocate(tau_r_twin(prm%totalNtwin, maxNinstance),              source=0.0_pReal)
 allocate(tau_r_trans(prm%totalNtrans, maxNinstance),            source=0.0_pReal)

 allocate(forestProjectionEdge(prm%totalNslip,prm%totalNslip,maxNinstance), source=0.0_pReal)
 allocate(projectionMatrix_Trans(prm%totalNtrans,prm%totalNslip,maxNinstance), source=0.0_pReal)
 
 allocate(Ctwin66(6,6,prm%totalNtwin,maxNinstance),              source=0.0_pReal)
 allocate(Ctrans66(6,6,prm%totalNtrans,maxNinstance),            source=0.0_pReal)

 allocate(Ctwin3333(3,3,3,3,prm%totalNtwin),        source=0.0_pReal)
 allocate(Ctrans3333(3,3,3,3,prm%totalNtrans),      source=0.0_pReal)
 

 initializeInstances: do p = 1_pInt, size(phase_plasticity)
    if (phase_plasticity(p) /= PLASTICITY_dislotwin_ID) cycle
     NofMyPhase=count(material_phase==p)
     instance = phase_plasticityInstance(p)
     prm => param(instance) 
!--------------------------------------------------------------------------------------------------
! allocate state arrays

     sizeDotState     = int(size(['rho         ','rhoDip      ','accshearslip']),pInt) * prm%totalNslip &
                      + int(size(['twinFraction','accsheartwin']),pInt) * prm%totalNtwin &
                      + int(size(['stressTransFraction','strainTransFraction']),pInt) * prm%totalNtrans
     sizeDeltaState   =  0_pInt
     sizeState        = sizeDotState &
                      + int(size(['invLambdaSlip     ','invLambdaSlipTwin ','invLambdaSlipTrans',&
                                  'meanFreePathSlip  ','tauSlipThreshold  ']),pInt) * prm%totalNslip &
                      + int(size(['invLambdaTwin   ','meanFreePathTwin','tauTwinThreshold',&
                                  'twinVolume      ']),pInt) * prm%totalNtwin &
                      + int(size(['invLambdaTrans   ','meanFreePathTrans','tauTransThreshold', &
                                  'martensiteVolume ']),pInt) * prm%totalNtrans
                
     plasticState(p)%sizeState = sizeState
     plasticState(p)%sizeDotState = sizeDotState
     plasticState(p)%sizeDeltaState = sizeDeltaState
     plasticState(p)%sizePostResults = sum(plastic_dislotwin_sizePostResult(:,instance))
     plasticState(p)%nSlip = prm%totalNslip
     plasticState(p)%nTwin = prm%totalNtwin
     plasticState(p)%nTrans= prm%totalNtrans
     allocate(plasticState(p)%aTolState           (sizeState),                source=0.0_pReal)
     allocate(plasticState(p)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(p)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(p)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(p)%state               (sizeState,NofMyPhase),     source=0.0_pReal)

     allocate(plasticState(p)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(plasticState(p)%deltaState        (sizeDeltaState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(p)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
       allocate(plasticState(p)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(p)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(p)%RKCK45dotState    (6,sizeDotState,NofMyPhase),source=0.0_pReal)
     offset_slip = 2_pInt*plasticState(p)%nslip
     plasticState(p)%slipRate        => &
         plasticState(p)%dotState(offset_slip+1:offset_slip+plasticState(p)%nslip,1:NofMyPhase)
     plasticState(p)%accumulatedSlip => &
         plasticState(p)%state   (offset_slip+1:offset_slip+plasticState(p)%nslip,1:NofMyPhase)

     allocate(temp1(prm%totalNslip,prm%totalNslip), source =0.0_pReal)
     allocate(temp2(prm%totalNslip,prm%totalNtwin), source =0.0_pReal)        
     allocate(temp3(prm%totalNslip,prm%totalNtrans),source =0.0_pReal)        

     mySlipFamilies: do f = 1_pInt,size(prm%Nslip,1)
       index_myFamily = sum(prm%Nslip(1:f-1_pInt))

       slipSystemsLoop: do j = 1_pInt,prm%Nslip(f)
         do o = 1_pInt, size(prm%Nslip,1)
           index_otherFamily = sum(prm%Nslip(1:o-1_pInt))
           do k = 1_pInt,prm%Nslip(o)                                   ! loop over (active) systems in other family (slip)
             forestProjectionEdge(index_myFamily+j,index_otherFamily+k,instance) = &
               abs(math_mul3x3(lattice_sn(:,sum(lattice_NslipSystem(1:f-1,p))+j,p), &
                               lattice_st(:,sum(lattice_NslipSystem(1:o-1,p))+k,p)))
             temp1(index_myFamily+j,index_otherFamily+k) = &
                   prm%interaction_SlipSlip(lattice_interactionSlipSlip( &
                                                                 sum(lattice_NslipSystem(1:f-1,p))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1,p))+k, &
                                                                 p),1 )
         enddo; enddo
  
         do o = 1_pInt,size(prm%Ntwin,1)
           index_otherFamily = sum(prm%Ntwin(1:o-1_pInt))
           do k = 1_pInt,prm%Ntwin(o)                                   ! loop over (active) systems in other family (twin)
             temp2(index_myFamily+j,index_otherFamily+k) = &
                   prm%interaction_SlipTwin(lattice_interactionSlipTwin( &
                                                                 sum(lattice_NslipSystem(1:f-1_pInt,p))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,p))+k, &
                                                                 p),1 )
         enddo; enddo

         do o = 1_pInt,size(prm%Ntrans,1)
           index_otherFamily = sum(prm%Ntrans(1:o-1_pInt))
           do k = 1_pInt,prm%Ntrans(o)                                  ! loop over (active) systems in other family (trans)
             temp3(index_myFamily+j,index_otherFamily+k) = &
                   prm%interaction_SlipTrans(lattice_interactionSlipTrans( &
                                                                 sum(lattice_NslipSystem(1:f-1_pInt,p))+j, &
                                                                 sum(lattice_NtransSystem(1:o-1_pInt,p))+k, &
                                                                 p),1 )
         enddo; enddo
  
       enddo slipSystemsLoop
     enddo mySlipFamilies
     prm%interaction_SlipSlip  = temp1; deallocate(temp1)    
     prm%interaction_SlipTwin  = temp2; deallocate(temp2)    
     prm%interaction_SlipTrans = temp3; deallocate(temp3)    

 
     allocate(temp1(prm%totalNtwin,prm%totalNslip), source =0.0_pReal)
     allocate(temp2(prm%totalNtwin,prm%totalNtwin), source =0.0_pReal)        

     Ctwin3333 = 0.0_pReal

     twinFamiliesLoop: do f = 1_pInt, size(prm%Ntwin,1)
       index_myFamily = sum(prm%Ntwin(1:f-1_pInt))                      ! index in truncated twin system list
       twinSystemsLoop: do j = 1_pInt,prm%Ntwin(f)
 
       !  nucleation rate prefactor,
       !  and twin size
       !* Rotate twin elasticity matrices
         index_otherFamily = sum(lattice_NtwinSystem(1:f-1_pInt,p))                             ! index in full lattice twin list
         do l = 1_pInt,3_pInt; do m = 1_pInt,3_pInt; do n = 1_pInt,3_pInt; do o = 1_pInt,3_pInt
           do p1 = 1_pInt,3_pInt; do q = 1_pInt,3_pInt; do r = 1_pInt,3_pInt; do s = 1_pInt,3_pInt
             Ctwin3333(l,m,n,o,index_myFamily+j) = &
             Ctwin3333(l,m,n,o,index_myFamily+j) + &
               lattice_C3333(p1,q,r,s,p) * &
               lattice_Qtwin(l,p1,index_otherFamily+j,p) * &
               lattice_Qtwin(m,q,index_otherFamily+j,p) * &
               lattice_Qtwin(n,r,index_otherFamily+j,p) * &
               lattice_Qtwin(o,s,index_otherFamily+j,p)
           enddo; enddo; enddo; enddo
         enddo; enddo; enddo; enddo
         Ctwin66(1:6,1:6,index_myFamily+j,instance) = &
           math_Mandel3333to66(Ctwin3333(1:3,1:3,1:3,1:3,index_myFamily+j))
 
      !* Interaction matrices
         do o = 1_pInt,size(prm%Nslip,1)
           index_otherFamily = sum(prm%Nslip(1:o-1_pInt))
           do k = 1_pInt,prm%Nslip(o)                                   ! loop over (active) systems in other family (slip)
             temp1(index_myFamily+j,index_otherFamily+k) = &
                   prm%interaction_TwinSlip(lattice_interactionTwinSlip( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,p))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1_pInt,p))+k, &
                                                                 p),1 )
         enddo; enddo
 
         do o = 1_pInt,size(prm%Ntwin,1)
           index_otherFamily = sum(prm%Ntwin(1:o-1_pInt))
           do k = 1_pInt,prm%Ntwin(o)                                   ! loop over (active) systems in other family (twin)
             temp2(index_myFamily+j,index_otherFamily+k) = &
                  prm%interaction_TwinTwin(lattice_interactionTwinTwin( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,p))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,p))+k, &
                                                                 p),1 )
         enddo; enddo
 
       enddo twinSystemsLoop
     enddo twinFamiliesLoop
     prm%interaction_TwinSlip  = temp1; deallocate(temp1)    
     prm%interaction_TwinTwin  = temp2; deallocate(temp2)    

 
     allocate(temp1(prm%totalNtrans,prm%totalNslip), source =0.0_pReal)
     allocate(temp2(prm%totalNtrans,prm%totalNtrans), source =0.0_pReal)        


     transFamiliesLoop: do f = 1_pInt,size(prm%Ntrans,1)
       index_myFamily = sum(prm%Ntrans(1:f-1_pInt))                      ! index in truncated trans system list
       transSystemsLoop: do j = 1_pInt,prm%Ntrans(f)

       Ctrans3333 = 0.0_pReal
       index_otherFamily = sum(lattice_NtransSystem(1:f-1_pInt,p))                               ! index in full lattice trans list
         do l = 1_pInt,3_pInt; do m = 1_pInt,3_pInt; do n = 1_pInt,3_pInt; do o = 1_pInt,3_pInt
           do p1 = 1_pInt,3_pInt; do q = 1_pInt,3_pInt; do r = 1_pInt,3_pInt; do s = 1_pInt,3_pInt
             Ctrans3333(l,m,n,o,index_myFamily+j) = &
             Ctrans3333(l,m,n,o,index_myFamily+j) + &
               lattice_trans_C3333(p1,q,r,s,p) * &
               lattice_Qtrans(l,p1,index_otherFamily+j,p) * &
               lattice_Qtrans(m,q,index_otherFamily+j,p) * &
               lattice_Qtrans(n,r,index_otherFamily+j,p) * &
               lattice_Qtrans(o,s,index_otherFamily+j,p)
           enddo; enddo; enddo; enddo
         enddo; enddo; enddo; enddo
         Ctrans66(1:6,1:6,index_myFamily+j,instance) = &
           math_Mandel3333to66(Ctrans3333(1:3,1:3,1:3,1:3,index_myFamily+j))

       !* Interaction matrices
         do o = 1_pInt,size(prm%Nslip,1)
           index_otherFamily = sum(prm%Nslip(1:o-1_pInt))
           do k = 1_pInt,prm%Nslip(o)                                   ! loop over (active) systems in other family (slip)
             temp1(index_myFamily+j,index_otherFamily+k) = &
                   prm%interaction_TransSlip(lattice_interactionTransSlip( &
                                                                 sum(lattice_NtransSystem(1:f-1_pInt,p))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1_pInt,p))+k, &
                                                                 p) ,1 )
         enddo; enddo

         do o = 1_pInt,size(prm%Ntrans,1)
           index_otherFamily = sum(prm%Ntrans(1:o-1_pInt))
           do k = 1_pInt,prm%Ntrans(o)                                  ! loop over (active) systems in other family (trans)
             temp2(index_myFamily+j,index_otherFamily+k) = &
                   prm%interaction_TransTrans(lattice_interactionTransTrans( &
                                                                 sum(lattice_NtransSystem(1:f-1_pInt,p))+j, &
                                                                 sum(lattice_NtransSystem(1:o-1_pInt,p))+k, &
                                                                 p),1 )
         enddo; enddo

       !* Projection matrices for shear from slip systems to fault-band (twin) systems for strain-induced martensite nucleation
       !  select case(trans_lattice_structure(p))
       !    case (LATTICE_bcc_ID)
       !      do o = 1_pInt,sum(prm%Ntrans,1)
       !        index_otherFamily = sum(prm%Nslip(1:o-1_pInt))
       !        do k = 1_pInt,prm%Nslip(o)                                   ! loop over (active) systems in other family (trans)
       !          temp3(index_myFamily+j,index_otherFamily+k) = &
       !                lattice_projectionTrans( sum(lattice_NtransSystem(1:f-1,p))+j, &
       !                                         sum(lattice_NslipSystem(1:o-1,p))+k, p)
       !      enddo; enddo
       !  end select 

       enddo transSystemsLoop
     enddo transFamiliesLoop

     
     prm%interaction_TransSlip  = temp1; deallocate(temp1)    
     prm%interaction_TransTrans  = temp2; deallocate(temp2)    

     startIndex=1_pInt
     endIndex=prm%totalNslip
     state(instance)%rhoEdge=>plasticState(p)%state(startIndex:endIndex,:)
     dotState(instance)%rhoEdge=>plasticState(p)%dotState(startIndex:endIndex,:)
     plasticState(p)%state0(startIndex:endIndex,:) = &
        spread(math_expand(prm%rho0,prm%Nslip),2,NofMyPhase)
     plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolRho

     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNslip
     state(instance)%rhoEdgeDip=>plasticState(p)%state(startIndex:endIndex,:)
     dotState(instance)%rhoEdgeDip=>plasticState(p)%dotState(startIndex:endIndex,:)
     plasticState(p)%state0(startIndex:endIndex,:) = &
        spread(math_expand(prm%rhoDip0,prm%Nslip),2,NofMyPhase)
     plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolRho
     
     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNslip
     state(instance)%accshear_slip=>plasticState(p)%state(startIndex:endIndex,:)
     dotState(instance)%accshear_slip=>plasticState(p)%dotState(startIndex:endIndex,:)
     plasticState(p)%aTolState(startIndex:endIndex) = 1.0e6_pReal
     
     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNtwin
     state(instance)%twinFraction=>plasticState(p)%state(startIndex:endIndex,:)
     dotState(instance)%twinFraction=>plasticState(p)%dotState(startIndex:endIndex,:)
     plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolTwinFrac
     
     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNtwin
     state(instance)%accshear_twin=>plasticState(p)%state(startIndex:endIndex,:)
     dotState(instance)%accshear_twin=>plasticState(p)%dotState(startIndex:endIndex,:)
     plasticState(p)%aTolState(startIndex:endIndex) = 1.0e6_pReal
     
     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNtrans
     state(instance)%stressTransFraction=>plasticState(p)%state(startIndex:endIndex,:)
     dotState(instance)%stressTransFraction=>plasticState(p)%dotState(startIndex:endIndex,:)
     plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolTransFrac
     
     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNtrans
     state(instance)%strainTransFraction=>plasticState(p)%state(startIndex:endIndex,:)
     dotState(instance)%strainTransFraction=>plasticState(p)%dotState(startIndex:endIndex,:)
     plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolTransFrac
     
     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNslip
     state(instance)%invLambdaSlip=>plasticState(p)%state(startIndex:endIndex,:)
     invLambdaSlip0 = spread(0.0_pReal,1,prm%totalNslip)
     forall (i = 1_pInt:prm%totalNslip) &
       invLambdaSlip0(i) = sqrt(dot_product(math_expand(prm%rho0,prm%Nslip)+ &
       math_expand(prm%rhoDip0,prm%Nslip),forestProjectionEdge(1:prm%totalNslip,i,instance)))/ &
                           prm%CLambdaSlipPerSlipSystem(i)
     plasticState(p)%state0(startIndex:endIndex,:) = &
       spread(math_expand(invLambdaSlip0,prm%Nslip),2, NofMyPhase)
     
     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNslip
     state(instance)%invLambdaSlipTwin=>plasticState(p)%state(startIndex:endIndex,:)
     plasticState(p)%state0(startIndex:endIndex,:) = 0.0_pReal

     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNtwin
     state(instance)%invLambdaTwin=>plasticState(p)%state(startIndex:endIndex,:)
     plasticState(p)%state0(startIndex:endIndex,:) = 0.0_pReal

     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNslip
     state(instance)%invLambdaSlipTrans=>plasticState(p)%state(startIndex:endIndex,:)
     plasticState(p)%state0(startIndex:endIndex,:) = 0.0_pReal

     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNtrans
     state(instance)%invLambdaTrans=>plasticState(p)%state(startIndex:endIndex,:)
     plasticState(p)%state0(startIndex:endIndex,:) = 0.0_pReal

     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNslip
     state(instance)%mfp_slip=>plasticState(p)%state(startIndex:endIndex,:)
     MeanFreePathSlip0 = prm%GrainSize/(1.0_pReal+invLambdaSlip0*prm%GrainSize)
     plasticState(p)%state0(startIndex:endIndex,:) = &
       spread(math_expand(MeanFreePathSlip0,prm%Nslip),2, NofMyPhase)
     
     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNtwin
     state(instance)%mfp_twin=>plasticState(p)%state(startIndex:endIndex,:)
     MeanFreePathTwin0 = spread(prm%GrainSize,1,prm%totalNtwin)
     plasticState(p)%state0(startIndex:endIndex,:) = &
       spread(math_expand(MeanFreePathTwin0,prm%Ntwin),2, NofMyPhase)

     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNtrans
     state(instance)%mfp_trans=>plasticState(p)%state(startIndex:endIndex,:)
     MeanFreePathTrans0 = spread(prm%GrainSize,1,prm%totalNtrans)
     plasticState(p)%state0(startIndex:endIndex,:) = &
       spread(math_expand(MeanFreePathTrans0,prm%Ntrans),2, NofMyPhase)

     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNslip
     state(instance)%threshold_stress_slip=>plasticState(p)%state(startIndex:endIndex,:)
     tauSlipThreshold0 = spread(0.0_pReal,1,prm%totalNslip)
     forall (i = 1_pInt:prm%totalNslip) tauSlipThreshold0(i) = &
       lattice_mu(p)*prm%burgers_slip(i) * &
       sqrt(dot_product(math_expand(prm%rho0 + prm%rhoDip0,prm%Nslip),&
                        prm%interaction_SlipSlip(i,1:prm%totalNslip)))
     plasticState(p)%state0(startIndex:endIndex,:) = &
       spread(math_expand(tauSlipThreshold0,prm%Nslip),2, NofMyPhase)

     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNtwin
     state(instance)%threshold_stress_twin=>plasticState(p)%state(startIndex:endIndex,:)

     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNtrans
     state(instance)%threshold_stress_trans=>plasticState(p)%state(startIndex:endIndex,:)

     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNtwin
     state(instance)%twinVolume=>plasticState(p)%state(startIndex:endIndex,:)
     TwinVolume0= spread(0.0_pReal,1,prm%totalNtwin)
     forall (i = 1_pInt:prm%totalNtwin) TwinVolume0(i) = &
       (PI/4.0_pReal)*prm%twinsize(i)*MeanFreePathTwin0(i)**2.0_pReal
     plasticState(p)%state0(startIndex:endIndex,:) = &
        spread(math_expand(TwinVolume0,prm%Ntwin),2, NofMyPhase)

     startIndex=endIndex+1
     endIndex=endIndex+prm%totalNtrans
     state(instance)%martensiteVolume=>plasticState(p)%state(startIndex:endIndex,:)
     MartensiteVolume0= spread(0.0_pReal,1,prm%totalNtrans)
     forall (i = 1_pInt:prm%totalNtrans) MartensiteVolume0(i) = &
       (PI/4.0_pReal)*prm%lamellarsizePerTransSystem(i)*MeanFreePathTrans0(i)**2.0_pReal
     plasticState(p)%state0(startIndex:endIndex,:) = &
       spread(math_expand(MartensiteVolume0,prm%Ntrans),2, NofMyPhase)

 enddo initializeInstances
end subroutine plastic_dislotwin_init




!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenized elasticity matrix
!--------------------------------------------------------------------------------------------------
function plastic_dislotwin_homogenizedC(ipc,ip,el)
 use material, only: &
  phase_plasticityInstance, &
  phaseAt, phasememberAt
 use lattice, only: &
  lattice_C66
 
  implicit none
  real(pReal), dimension(6,6) :: &
    plastic_dislotwin_homogenizedC
  integer(pInt), intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element
 type(tParameters)     :: prm
 type(tDislotwinState) :: ste

 integer(pInt) :: instance,i, &
                  ph, &
                  of
 real(pReal) :: sumf, sumftr

 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(ph)
 associate( prm => param(instance), ste =>state(instance))


 !* Total twin volume fraction
 sumf = sum(ste%twinFraction(1_pInt:prm%totalNtwin,of))           ! safe for prm%totalNtwin == 0

 !* Total transformed volume fraction
 sumftr = sum(ste%stressTransFraction(1_pInt:prm%totalNtrans,of)) + &
          sum(ste%strainTransFraction(1_pInt:prm%totalNtrans,of))

 !* Homogenized elasticity matrix
 plastic_dislotwin_homogenizedC = (1.0_pReal-sumf-sumftr)*lattice_C66(1:6,1:6,ph)
 do i=1_pInt,prm%totalNtwin
    plastic_dislotwin_homogenizedC = plastic_dislotwin_homogenizedC &
                   + ste%twinFraction(i,of)*Ctwin66(1:6,1:6,i,instance)
 enddo
 do i=1_pInt,prm%totalNtrans
    plastic_dislotwin_homogenizedC = plastic_dislotwin_homogenizedC &
                   + (ste%stressTransFraction(i,of) + ste%strainTransFraction(i,of))*&
                     Ctrans66(1:6,1:6,i,instance)
 enddo
 end associate
 end function plastic_dislotwin_homogenizedC
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_microstructure(temperature,ipc,ip,el)
 use math, only: &
   pi
 use material, only: &
   material_phase, &
   phase_plasticityInstance, &
   plasticState, &  
   phaseAt, phasememberAt
 use lattice, only: &
   lattice_mu, &
   lattice_nu

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in) :: &
   temperature                                                                                      !< temperature at IP 

 integer(pInt) :: &
   instance, &
   s, &
   ph, &
   of
 real(pReal) :: &
   sumf,sfe,sumftr
    real(pReal), dimension(:), allocatable :: &
   x0
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%Ntwin)  :: fOverStacksize
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%Ntrans) :: ftransOverLamellarSize
   
 type(tParameters):: prm
 type(tDislotwinState) :: ste

 
 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(ph)

 associate(prm => param(instance), ste => state(instance))
 !* Total twin volume fraction
 sumf   = sum(ste%twinFraction(1:prm%totalNtwin,of)) ! safe for prm%totalNtwin == 0
 sumftr = sum(ste%stressTransFraction(1:prm%totalNtrans,of)) + &
          sum(ste%strainTransFraction(1:prm%totalNtrans,of))

 !* Stacking fault energy
 sfe = prm%SFE_0K + prm%dSFE_dT * Temperature
 
 !* rescaled volume fraction for topology
 fOverStacksize         =  ste%twinFraction(1_pInt:prm%totalNtwin,of)/prm%twinsize
 ftransOverLamellarSize =  sumftr                                    /prm%lamellarsizePerTransSystem
 
 !* 1/mean free distance between 2 forest dislocations seen by a moving dislocation
 forall (s = 1_pInt:prm%totalNslip) &
   ste%invLambdaSlip(s,of) = &
     sqrt(dot_product((ste%rhoEdge(1_pInt:prm%totalNslip,of)+ste%rhoEdgeDip(1_pInt:prm%totalNslip,of)),&
                      forestProjectionEdge(1:prm%totalNslip,s,instance)))/prm%CLambdaSlipPerSlipSystem(s)

 !* 1/mean free distance between 2 twin stacks from different systems seen by a moving dislocation
 !$OMP CRITICAL (evilmatmul)
 if (prm%totalNtwin > 0_pInt .and. prm%totalNslip > 0_pInt) &
   ste%invLambdaSlipTwin(1_pInt:prm%totalNslip,of) = &
     matmul(prm%interaction_SlipTwin,fOverStacksize)/(1.0_pReal-sumf)

 !* 1/mean free distance between 2 twin stacks from different systems seen by a growing twin

  !ToDo: needed? if (prm%totalNtwin > 0_pInt) &
  ste%invLambdaTwin(1_pInt:prm%totalNtwin,of) = &
     matmul(prm%interaction_TwinTwin,fOverStacksize)/(1.0_pReal-sumf)


 !* 1/mean free distance between 2 martensite lamellar from different systems seen by a moving dislocation
 if (prm%totalNtrans > 0_pInt .and. prm%totalNslip > 0_pInt) &
   ste%invLambdaSlipTrans(1_pInt:prm%totalNslip,of) = &
      matmul(prm%interaction_SlipTrans,ftransOverLamellarSize)/(1.0_pReal-sumftr)

 !* 1/mean free distance between 2 martensite stacks from different systems seen by a growing martensite (1/lambda_trans)
 !ToDo: needed? if (prm%totalNtrans > 0_pInt) &

   ste%invLambdaTrans(1_pInt:prm%totalNtrans,of) = &
     matmul(prm%interaction_TransTrans,ftransOverLamellarSize)/(1.0_pReal-sumftr)
 !$OMP END CRITICAL (evilmatmul)

 !* mean free path between 2 obstacles seen by a moving dislocation
 do s = 1_pInt,prm%totalNslip
    if ((prm%totalNtwin > 0_pInt) .or. (prm%totalNtrans > 0_pInt)) then              ! ToDo: This is two simplified
       ste%mfp_slip(s,of) = &
         prm%GrainSize/(1.0_pReal+prm%GrainSize*&
         (ste%invLambdaSlip(s,of) + ste%invLambdaSlipTwin(s,of) + ste%invLambdaSlipTrans(s,of)))
    else
       ste%mfp_slip(s,of) = &  
         prm%GrainSize/&
         (1.0_pReal+prm%GrainSize*(ste%invLambdaSlip(s,of))) !!!!!! correct?
    endif
 enddo

 !* mean free path between 2 obstacles seen by a growing twin/martensite
 ste%mfp_twin(:,of)  = prm%Cmfptwin*prm%GrainSize/ (1.0_pReal+prm%GrainSize*ste%invLambdaTwin(:,of))
 ste%mfp_trans(:,of) = prm%Cmfptrans*prm%GrainSize/(1.0_pReal+prm%GrainSize*ste%invLambdaTrans(:,of))

 !* threshold stress for dislocation motion
 forall (s = 1_pInt:prm%totalNslip) ste%threshold_stress_slip(s,of) = &
     lattice_mu(ph)*prm%burgers_slip(s)*&
     sqrt(dot_product(ste%rhoEdge(1_pInt:prm%totalNslip,of)+ste%rhoEdgeDip(1_pInt:prm%totalNslip,of),&
                      prm%interaction_SlipSlip(s,1:prm%totalNslip)))

 !* threshold stress for growing twin/martensite
   ste%threshold_stress_twin(:,of) = prm%Cthresholdtwin* &
     (sfe/(3.0_pReal*prm%burgers_twin)+ 3.0_pReal*prm%burgers_twin*lattice_mu(ph)/ &
           (prm%L0_twin*prm%burgers_slip)) ! slip burgers here correct?
   ste%threshold_stress_trans(:,of) = prm%Cthresholdtrans* &
       (sfe/(3.0_pReal*prm%burgers_trans) + 3.0_pReal*prm%burgers_trans*lattice_mu(ph)/&
            (prm%L0_trans*prm%burgers_slip) + prm%transStackHeight*prm%deltaG/ (3.0_pReal*prm%burgers_trans) )    
 
  ! final volume after growth
  ste%twinVolume(:,of) = (PI/4.0_pReal)*prm%twinsize*ste%mfp_twin(:,of)**2.0_pReal
  ste%martensiteVolume(:,of) = (PI/4.0_pReal)*prm%lamellarsizePerTransSystem*ste%mfp_trans(:,of)**2.0_pReal

 !* equilibrium separation of partial dislocations (twin)
 x0 = lattice_mu(ph)*prm%burgers_twin**2.0_pReal/(sfe*8.0_pReal*PI)*(2.0_pReal+lattice_nu(ph))/(1.0_pReal-lattice_nu(ph))
 tau_r_twin(:,instance)= lattice_mu(ph)*prm%burgers_twin/(2.0_pReal*PI)*&     
        (1/(x0+prm%xc_twin)+cos(pi/3.0_pReal)/x0)
 
!* equilibrium separation of partial dislocations (trans)
 x0 = lattice_mu(ph)*prm%burgers_trans**2.0_pReal/(sfe*8.0_pReal*PI)*(2.0_pReal+lattice_nu(ph))/(1.0_pReal-lattice_nu(ph))
 tau_r_trans(:,instance)= lattice_mu(ph)*prm%burgers_trans/(2.0_pReal*PI)*&     
        (1/(x0+prm%xc_trans)+cos(pi/3.0_pReal)/x0)

end associate
end subroutine plastic_dislotwin_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_LpAndItsTangent(Lp,dLp_dTstar99,Tstar_v,Temperature,ipc,ip,el)
 use prec, only: &
   tol_math_check, &
   dNeq0
 use math, only: &
   math_Plain3333to99, &
   math_Mandel6to33, &
   math_Mandel33to6, &
   math_eigenValuesVectorsSym, &
   math_tensorproduct33, &
   math_symmetric33, &
   math_mul33x3
 use material, only: &
   material_phase, &
   plasticState, &
   phase_plasticityInstance, &
   phaseAt, phasememberAt
 use lattice, only: &
   lattice_Sslip, &
   lattice_Sslip_v, &
   lattice_Stwin, &
   lattice_Stwin_v, &
   lattice_Strans, &
   lattice_Strans_v, &
   lattice_maxNslipFamily,&
   lattice_maxNtwinFamily, &
   lattice_maxNtransFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_NtransSystem, &
   lattice_shearTwin, &
   lattice_structure, &
   lattice_fcc_twinNucleationSlipPair, &
   LATTICE_fcc_ID
 
 implicit none
 integer(pInt), intent(in)                  :: ipc,ip,el
 real(pReal), intent(in)                    :: Temperature
 real(pReal), dimension(6),   intent(in)    :: Tstar_v
 real(pReal), dimension(3,3), intent(out)   :: Lp
 real(pReal), dimension(9,9), intent(out)   :: dLp_dTstar99

 integer(pInt) :: instance,ph,of,f,i,j,k,l,m,n,index_myFamily,s1,s2
 real(pReal) :: sumf,sumftr,StressRatio_p,StressRatio_pminus1,StressRatio_r,BoltzmannRatio,DotGamma0,Ndot0_twin,stressRatio, &
    Ndot0_trans,StressRatio_s
 real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%Nslip) :: &
    gdot_slip,dgdot_dtauslip,tau_slip
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%Ntwin) :: &
    gdot_twin,dgdot_dtautwin,tau_twin
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%Ntrans) :: &
    gdot_trans,dgdot_dtautrans,tau_trans
 real(pReal), dimension(6) :: gdot_sb,dgdot_dtausb,tau_sb
 real(pReal), dimension(3,3) :: eigVectors, sb_Smatrix
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
        
  type(tParameters) :: prm
  
 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance  = phase_plasticityInstance(ph)


 Lp = 0.0_pReal
 dLp_dTstar3333 = 0.0_pReal
 
 
associate(prm => param(instance))
!--------------------------------------------------------------------------------------------------
! Dislocation glide part
 gdot_slip = 0.0_pReal
 dgdot_dtauslip = 0.0_pReal
 j = 0_pInt
 slipFamiliesLoop: do f = 1_pInt,size(prm%Nslip,1)
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph)) ! at which index starts my family
   slipSystemsLoop: do i = 1_pInt,prm%Nslip(f)
      j = j+1_pInt
 
      !* Calculation of Lp
      !* Resolved shear stress on slip system
      tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))

      if((abs(tau_slip(j))-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
      !* Stress ratios
        stressRatio =((abs(tau_slip(j))- state(instance)%threshold_stress_slip(j,of))/&
          (prm%SolidSolutionStrength+prm%tau_peierlsPerSlipFamily(f)))
        StressRatio_p       = stressRatio** prm%p(f)
        StressRatio_pminus1 = stressRatio**(prm%p(f)-1.0_pReal)
      !* Boltzmann ratio
        BoltzmannRatio = prm%Qedge(j)/(kB*Temperature)
      !* Initial shear rates
        DotGamma0 = &
          state(instance)%rhoEdge(j,of)*prm%burgers_slip(j)*&
          prm%v0(j)
 
        !* Shear rates due to slip
        gdot_slip(j) = DotGamma0 &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(f)) &
                     * sign(1.0_pReal,tau_slip(j))
 
        !* Derivatives of shear rates
        dgdot_dtauslip(j) = &
          abs(gdot_slip(j))*BoltzmannRatio*prm%p(f)&
          *prm%q(f)/&
          (prm%SolidSolutionStrength+prm%tau_peierlsPerSlipFamily(f))*&
          StressRatio_pminus1*(1-StressRatio_p)**(prm%q(f)-1.0_pReal)
      endif
 
      !* Plastic velocity gradient for dislocation glide
      Lp = Lp + gdot_slip(j)*lattice_Sslip(:,:,1,index_myFamily+i,ph)
 
      !* Calculation of the tangent of Lp
      forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtauslip(j)*&
                                  lattice_Sslip(k,l,1,index_myFamily+i,ph)*&
                                  lattice_Sslip(m,n,1,index_myFamily+i,ph)
   enddo slipSystemsLoop
 enddo slipFamiliesLoop
 
!--------------------------------------------------------------------------------------------------
! correct Lp and dLp_dTstar3333 for twinned and transformed fraction 
 !* Total twin volume fraction
 sumf = sum(state(instance)%twinFraction(1_pInt:prm%totalNtwin,of)) ! safe for prm%totalNtwin == 0

 !* Total transformed volume fraction
 sumftr = sum(state(instance)%stressTransFraction(1_pInt:prm%totalNtrans,of)) + &
          sum(state(instance)%strainTransFraction(1_pInt:prm%totalNtrans,of))
 Lp = Lp * (1.0_pReal - sumf - sumftr)
 dLp_dTstar3333 = dLp_dTstar3333 * (1.0_pReal - sumf - sumftr)
 
!--------------------------------------------------------------------------------------------------
! Shear banding (shearband) part
 if(dNeq0(prm%sbVelocity)) then
   gdot_sb = 0.0_pReal
   dgdot_dtausb = 0.0_pReal
   call math_eigenValuesVectorsSym(math_Mandel6to33(Tstar_v),eigValues,eigVectors,error)
   do j = 1_pInt,6_pInt
     sb_s = 0.5_pReal*sqrt(2.0_pReal)*math_mul33x3(eigVectors,sb_sComposition(1:3,j))
     sb_m = 0.5_pReal*sqrt(2.0_pReal)*math_mul33x3(eigVectors,sb_mComposition(1:3,j))
     sb_Smatrix = math_tensorproduct33(sb_s,sb_m)
     sbSv(1:6,j,ipc,ip,el) = math_Mandel33to6(math_symmetric33(sb_Smatrix))
   
     !* Calculation of Lp
     !* Resolved shear stress on shear banding system
     tau_sb(j) = dot_product(Tstar_v,sbSv(1:6,j,ipc,ip,el))
   
     !* Stress ratios
      if (abs(tau_sb(j)) < tol_math_check) then
        StressRatio_p = 0.0_pReal
        StressRatio_pminus1 = 0.0_pReal
      else
       StressRatio_p = (abs(tau_sb(j))/prm%sbResistance)&
                 **prm%pShearBand
       StressRatio_pminus1 = (abs(tau_sb(j))/prm%sbResistance)&
                 **(prm%pShearBand-1.0_pReal)
      endif

     !* Boltzmann ratio
     BoltzmannRatio = prm%sbQedge/(kB*Temperature)
     !* Initial shear rates
     DotGamma0 = prm%sbVelocity
 
     !* Shear rates due to shearband
     gdot_sb(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**&
                  prm%qShearBand)*sign(1.0_pReal,tau_sb(j))
                  
     !* Derivatives of shear rates
     dgdot_dtausb(j) = &
       ((abs(gdot_sb(j))*BoltzmannRatio*&
       prm%pShearBand*prm%qShearBand)/&
       prm%sbResistance)*&
       StressRatio_pminus1*(1_pInt-StressRatio_p)**(prm%qShearBand-1.0_pReal)
 
     !* Plastic velocity gradient for shear banding
     Lp = Lp + gdot_sb(j)*sb_Smatrix
 
     !* Calculation of the tangent of Lp
     forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
       dLp_dTstar3333(k,l,m,n) = &
       dLp_dTstar3333(k,l,m,n) + dgdot_dtausb(j)*&
                                 sb_Smatrix(k,l)*&
                                 sb_Smatrix(m,n)
   enddo
 end if
 
!--------------------------------------------------------------------------------------------------
! Mechanical twinning part
 gdot_twin = 0.0_pReal
 dgdot_dtautwin = 0.0_pReal
 j = 0_pInt
 twinFamiliesLoop: do f = 1_pInt,size(prm%Ntwin,1)
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph)) ! at which index starts my family
   twinSystemsLoop: do i = 1_pInt,prm%Ntwin(f)
      j = j+1_pInt
 
      !* Calculation of Lp
      !* Resolved shear stress on twin system
      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,ph))

      !* Stress ratios
      if (tau_twin(j) > tol_math_check) then
        StressRatio_r = (state(instance)%threshold_stress_twin(j,of)/tau_twin(j))**prm%r(f)
        !* Shear rates and their derivatives due to twin
        select case(lattice_structure(ph))
          case (LATTICE_fcc_ID)
            s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
            s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
            if (tau_twin(j) < tau_r_twin(j,instance)) then
              Ndot0_twin=(abs(gdot_slip(s1))*(state(instance)%rhoEdge(s2,of)+state(ph)%rhoEdgeDip(s2,of))+&  !!!!! correct?
                     abs(gdot_slip(s2))*(state(instance)%rhoEdge(s1,of)+state(instance)%rhoEdgeDip(s1,of)))/&
                    (prm%L0_twin*prm%burgers_slip(j))*&
                    (1.0_pReal-exp(-prm%VcrossSlip/(kB*Temperature)*&
                    (tau_r_twin(j,instance)-tau_twin(j))))
            else
              Ndot0_twin=0.0_pReal
            end if
          case default
            Ndot0_twin=prm%Ndot0_twin(j)
        end select
        gdot_twin(j) = &
          (1.0_pReal-sumf-sumftr)*lattice_shearTwin(index_myFamily+i,ph)*&
          state(instance)%twinVolume(j,of)*Ndot0_twin*exp(-StressRatio_r)
        dgdot_dtautwin(j) = ((gdot_twin(j)*prm%r(j))/tau_twin(j))*StressRatio_r
      endif
 
      !* Plastic velocity gradient for mechanical twinning
      Lp = Lp + gdot_twin(j)*lattice_Stwin(:,:,index_myFamily+i,ph)
      !* Calculation of the tangent of Lp
      forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtautwin(j)*&
                                  lattice_Stwin(k,l,index_myFamily+i,ph)*&
                                  lattice_Stwin(m,n,index_myFamily+i,ph)
   enddo twinSystemsLoop
 enddo twinFamiliesLoop

 !* Phase transformation part
 gdot_trans = 0.0_pReal
 dgdot_dtautrans = 0.0_pReal
 j = 0_pInt
 transFamiliesLoop: do f = 1_pInt,size(prm%Ntrans,1)
   index_myFamily = sum(lattice_NtransSystem(1:f-1_pInt,ph)) ! at which index starts my family
   transSystemsLoop: do i = 1_pInt,prm%Ntrans(f)
      j = j+1_pInt

      !* Resolved shear stress on transformation system
      tau_trans(j) = dot_product(Tstar_v,lattice_Strans_v(:,index_myFamily+i,ph))

      !* Stress ratios
      if (tau_trans(j) > tol_math_check) then
        StressRatio_s = (state(instance)%threshold_stress_trans(j,of)/tau_trans(j))**prm%s(f)
        !* Shear rates and their derivatives due to transformation
        select case(lattice_structure(ph))
          case (LATTICE_fcc_ID)
            s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
            s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
            if (tau_trans(j) < tau_r_trans(j,instance)) then
              Ndot0_trans=(abs(gdot_slip(s1))*(state(instance)%rhoEdge(s2,of)+state(instance)%rhoEdgeDip(s2,of))+&  !!!!! correct?
                           abs(gdot_slip(s2))*(state(instance)%rhoEdge(s1,of)+state(instance)%rhoEdgeDip(s1,of)))/&
                          (prm%L0_trans*prm%burgers_slip(j))*&
                          (1.0_pReal-exp(-prm%VcrossSlip/(kB*Temperature)*&
                          (tau_r_trans(j,instance)-tau_trans(j))))
            else
              Ndot0_trans=0.0_pReal
            end if
          case default
            Ndot0_trans=prm%Ndot0_trans(j)
        end select
        gdot_trans(j) = &
          (1.0_pReal-sumf-sumftr)*&
          state(instance)%martensiteVolume(j,of)*Ndot0_trans*exp(-StressRatio_s)
        dgdot_dtautrans(j) = ((gdot_trans(j)*prm%s(f))/tau_trans(j))*StressRatio_s
      endif

      !* Plastic velocity gradient for phase transformation
      Lp = Lp + gdot_trans(j)*lattice_Strans(:,:,index_myFamily+i,ph)

      !* Calculation of the tangent of Lp
      forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtautrans(j)*&
                                  lattice_Strans(k,l,index_myFamily+i,ph)*&
                                  lattice_Strans(m,n,index_myFamily+i,ph)

   enddo transSystemsLoop
 enddo transFamiliesLoop
 end associate
 dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)
 
end subroutine plastic_dislotwin_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_dotState(Tstar_v,Temperature,ipc,ip,el)
 use prec, only: &
   tol_math_check, &
   dEq0
 use math, only: &
   pi
 use material, only: &
   material_phase, &
   phase_plasticityInstance, &
   plasticState, &
   phaseAt, phasememberAt
 use lattice,  only: &
   lattice_Sslip_v, &
   lattice_Stwin_v, &
   lattice_Strans_v, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_maxNtransFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_NtransSystem, &   
   lattice_sheartwin, &
   lattice_mu, &
   lattice_structure, &
   lattice_fcc_twinNucleationSlipPair, &
   lattice_fccTobcc_transNucleationTwinPair, &
   lattice_fccTobcc_shearCritTrans, &
   LATTICE_fcc_ID

 implicit none
 real(pReal), dimension(6),  intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 integer(pInt) :: instance,f,i,j,index_myFamily,s1,s2, &
                  ph, &
                  of
 real(pReal) :: sumf,sumftr,StressRatio_p,StressRatio_pminus1,BoltzmannRatio,DotGamma0,&
             EdgeDipMinDistance,AtomicVolume,VacancyDiffusion,StressRatio_r,Ndot0_twin,stressRatio,&
             Ndot0_trans,StressRatio_s,EdgeDipDistance, ClimbVelocity,DotRhoEdgeDipClimb,DotRhoEdgeDipAnnihilation, &
             DotRhoDipFormation,DotRhoMultiplication,DotRhoEdgeEdgeAnnihilation
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%Nslip) :: &
 gdot_slip,tau_slip

 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%Ntwin) :: &
              tau_twin
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%Ntrans) :: &
              tau_trans

 type(tParameters) :: prm

 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance  = phase_plasticityInstance(ph)

 associate(prm => param(instance))
 !* Total twin volume fraction
 sumf = sum(state(instance)%twinFraction(1_pInt:prm%totalNtwin,of)) ! safe for prm%totalNtwin == 0
 plasticState(ph)%dotState(:,of) = 0.0_pReal
 
 !* Total transformed volume fraction
 sumftr = sum(state(instance)%stressTransFraction(1_pInt:prm%totalNtrans,of)) + &
          sum(state(instance)%strainTransFraction(1_pInt:prm%totalNtrans,of))
 
 !* Dislocation density evolution
 gdot_slip = 0.0_pReal
 j = 0_pInt
 do f = 1_pInt,size(prm%Nslip,1)
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))   ! at which index starts my family
   do i = 1_pInt,prm%Nslip(f)          ! process each (active) slip system in family
      j = j+1_pInt
 
      !* Resolved shear stress on slip system
      tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))

      if((abs(tau_slip(j))-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
      !* Stress ratios
        stressRatio =((abs(tau_slip(j))- state(instance)%threshold_stress_slip(j,of))/&
          (prm%SolidSolutionStrength+prm%tau_peierlsPerSlipFamily(f)))
        StressRatio_p       = stressRatio** prm%p(f)
        StressRatio_pminus1 = stressRatio**(prm%p(f)-1.0_pReal)
      !* Boltzmann ratio
        BoltzmannRatio = prm%Qedge(j)/(kB*Temperature)
      !* Initial shear rates
        DotGamma0 = plasticState(ph)%state(j, of)*prm%burgers_slip(j)*prm%v0(j)
 
      !* Shear rates due to slip
        gdot_slip(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)** &
                       prm%q(f))*sign(1.0_pReal,tau_slip(j))
      endif
      !* Multiplication
      DotRhoMultiplication = abs(gdot_slip(j))/(prm%burgers_slip(j)*state(instance)%mfp_slip(j,of))

      !* Dipole formation
      EdgeDipMinDistance = prm%CEdgeDipMinDistance*prm%burgers_slip(j)
      if (dEq0(tau_slip(j))) then
        DotRhoDipFormation = 0.0_pReal
      else
        EdgeDipDistance = &
          (3.0_pReal*lattice_mu(ph)*prm%burgers_slip(j))/&
          (16.0_pReal*pi*abs(tau_slip(j)))
        if (EdgeDipDistance>state(instance)%mfp_slip(j,of)) EdgeDipDistance=state(instance)%mfp_slip(j,of)
        if (EdgeDipDistance<EdgeDipMinDistance)       EdgeDipDistance=EdgeDipMinDistance
        DotRhoDipFormation = &
          ((2.0_pReal*(EdgeDipDistance-EdgeDipMinDistance))/prm%burgers_slip(j))*&
          state(instance)%rhoEdge(j,of)*abs(gdot_slip(j))*prm%dipoleFormationFactor
      endif
 
      !* Spontaneous annihilation of 2 single edge dislocations
      DotRhoEdgeEdgeAnnihilation = &
        ((2.0_pReal*EdgeDipMinDistance)/prm%burgers_slip(j))*&
        state(instance)%rhoEdge(j,of)*abs(gdot_slip(j))
 
      !* Spontaneous annihilation of a single edge dislocation with a dipole constituent
      DotRhoEdgeDipAnnihilation = &
        ((2.0_pReal*EdgeDipMinDistance)/prm%burgers_slip(j))*&
        state(instance)%rhoEdgeDip(j,of)*abs(gdot_slip(j))
 
      !* Dislocation dipole climb
      AtomicVolume = &
        prm%CAtomicVolume*prm%burgers_slip(j)**(3.0_pReal)
      VacancyDiffusion = &
        prm%D0*exp(-prm%Qsd/(kB*Temperature))
      if (dEq0(tau_slip(j))) then
        DotRhoEdgeDipClimb = 0.0_pReal
      else
        if (dEq0(EdgeDipDistance-EdgeDipMinDistance)) then
          DotRhoEdgeDipClimb = 0.0_pReal
        else
          ClimbVelocity = 3.0_pReal*lattice_mu(ph)*VacancyDiffusion*AtomicVolume/ &
                          (2.0_pReal*pi*kB*Temperature*(EdgeDipDistance+EdgeDipMinDistance))
          DotRhoEdgeDipClimb = 4.0_pReal*ClimbVelocity*state(instance)%rhoEdgeDip(j,of)/ &
                          (EdgeDipDistance-EdgeDipMinDistance)
        endif
      endif
      !* Edge dislocation density rate of change
      dotState(instance)%rhoEdge(j,of) = DotRhoMultiplication-DotRhoDipFormation-DotRhoEdgeEdgeAnnihilation
 
      !* Edge dislocation dipole density rate of change
      dotState(instance)%rhoEdgeDip(j,of) = DotRhoDipFormation-DotRhoEdgeDipAnnihilation-DotRhoEdgeDipClimb
 
      !* Dotstate for accumulated shear due to slip
      dotState(instance)%accshear_slip(j,of) = abs(gdot_slip(j))
 
   enddo
 enddo
 
 !* Twin volume fraction evolution
 j = 0_pInt
 do f = 1_pInt,size(prm%Ntwin,1)
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))   ! at which index starts my family
   do i = 1_pInt,prm%Ntwin(f)           ! process each (active) twin system in family
      j = j+1_pInt
 
      !* Resolved shear stress on twin system
      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,ph))
      !* Stress ratios
      if (tau_twin(j) > tol_math_check) then
        StressRatio_r = (state(instance)%threshold_stress_twin(j,of)/&
                         tau_twin(j))**prm%r(j)
        !* Shear rates and their derivatives due to twin
        select case(lattice_structure(ph))
          case (LATTICE_fcc_ID)
            s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
            s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
            if (tau_twin(j) < tau_r_twin(j,instance)) then
              Ndot0_twin=(abs(gdot_slip(s1))*(state(instance)%rhoEdge(s2,of)+state(instance)%rhoEdgeDip(s2,of))+&
                     abs(gdot_slip(s2))*(state(instance)%rhoEdge(s1,of)+state(instance)%rhoEdgeDip(s1,of)))/&
                    (prm%L0_twin*prm%burgers_slip(j))*&
                    (1.0_pReal-exp(-prm%VcrossSlip/(kB*Temperature)*&
                    (tau_r_twin(j,instance)-tau_twin(j))))
            else
              Ndot0_twin=0.0_pReal
            end if
          case default
            Ndot0_twin=prm%Ndot0_twin(j)
        end select
        dotState(instance)%twinFraction(j,of) = &
          (1.0_pReal-sumf-sumftr)*&
          state(instance)%twinVolume(j,of)*Ndot0_twin*exp(-StressRatio_r)
        !* Dotstate for accumulated shear due to twin
        dotState(instance)%accshear_twin(j,of) = dotState(instance)%twinFraction(j,of) * &
                                                          lattice_sheartwin(index_myfamily+i,ph)
      endif
   enddo
 enddo

 !* Transformation volume fraction evolution
 j = 0_pInt
 do f = 1_pInt,size(prm%Ntrans,1)
   index_myFamily = sum(lattice_NtransSystem(1:f-1_pInt,ph))  ! at which index starts my family
   do i = 1_pInt,prm%Ntrans(f)         ! process each (active) trans system in family
      j = j+1_pInt

      !* Resolved shear stress on transformation system
      tau_trans(j) = dot_product(Tstar_v,lattice_Strans_v(:,index_myFamily+i,ph))

      !* Stress ratios
      if (tau_trans(j) > tol_math_check) then
        StressRatio_s = (state(instance)%threshold_stress_trans(j,of)/&
                         tau_trans(j))**prm%s(f)
        !* Shear rates and their derivatives due to transformation
        select case(lattice_structure(ph))
          case (LATTICE_fcc_ID)
            s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
            s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
            if (tau_trans(j) < tau_r_trans(j,instance)) then
              Ndot0_trans=(abs(gdot_slip(s1))*(state(instance)%rhoEdge(s2,of)+state(instance)%rhoEdgeDip(s2,of))+&
                           abs(gdot_slip(s2))*(state(instance)%rhoEdge(s1,of)+state(instance)%rhoEdgeDip(s1,of)))/&
                          (prm%L0_trans*prm%burgers_slip(j))*&
                          (1.0_pReal-exp(-prm%VcrossSlip/(kB*Temperature)*&
                          (tau_r_trans(j,instance)-tau_trans(j))))
            else
              Ndot0_trans=0.0_pReal
            end if
          case default
            Ndot0_trans=prm%Ndot0_trans(j)
        end select
        dotState(instance)%strainTransFraction(j,of) = &
          (1.0_pReal-sumf-sumftr)*&
          state(instance)%martensiteVolume(j,of)*Ndot0_trans*exp(-StressRatio_s)
        !* Dotstate for accumulated shear due to transformation
        !dotState(instance)%accshear_trans(j,of) = dotState(instance)%strainTransFraction(j,of) * &
        !                                                  lattice_sheartrans(index_myfamily+i,ph)
      endif

   enddo
 enddo 
end associate
end subroutine plastic_dislotwin_dotState

 
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_dislotwin_postResults(Tstar_v,Temperature,ipc,ip,el)
 use prec, only: &
   tol_math_check, &
   dEq0
 use math, only: &
   pi, &
   math_Mandel6to33, &
   math_eigenValuesSym33, &
   math_eigenValuesVectorsSym33
 use material, only: &
   material_phase, &
   plasticState, &
   phase_plasticityInstance,& 
   phaseAt, phasememberAt
 use lattice, only: &
   lattice_Sslip_v, &
   lattice_Stwin_v, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_shearTwin, &
   lattice_mu, &
   lattice_structure, &
   lattice_fcc_twinNucleationSlipPair, &
   LATTICE_fcc_ID

 implicit none
 real(pReal), dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%sizePostResults) :: &
                                           plastic_dislotwin_postResults
 integer(pInt) :: &
   instance,&
   f,o,i,c,j,index_myFamily,&
   s1,s2, &
   ph, &
   of
 real(pReal) :: sumf,tau,StressRatio_p,StressRatio_pminus1,BoltzmannRatio,DotGamma0,StressRatio_r,Ndot0_twin,dgdot_dtauslip, &
   stressRatio
 real(preal), dimension(plasticState(material_phase(ipc,ip,el))%Nslip) :: &
   gdot_slip
 real(pReal), dimension(3,3) :: eigVectors
 real(pReal), dimension (3) :: eigValues
 
  type(tParameters) :: prm
 
 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance  = phase_plasticityInstance(ph)


 associate(prm => param(instance))
 !* Total twin volume fraction
 sumf = sum(state(instance)%twinFraction(1_pInt:prm%totalNtwin,of))                          ! safe for prm%totalNtwin == 0
 
 !* Required output
 c = 0_pInt
 plastic_dislotwin_postResults = 0.0_pReal
 do o = 1_pInt,size(prm%outputID)
    select case(prm%outputID(o))
 
      case (edge_density_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+prm%totalNslip) = state(instance)%rhoEdge(1_pInt:prm%totalNslip,of)
        c = c + prm%totalNslip
      case (dipole_density_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+prm%totalNslip) = state(instance)%rhoEdgeDip(1_pInt:prm%totalNslip,of)
        c = c + prm%totalNslip
      case (shear_rate_slip_ID)
        j = 0_pInt
        do f = 1_pInt,size(prm%Nslip,1)                                                        ! loop over all slip families
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
           do i = 1_pInt,prm%Nslip(f)                                        ! process each (active) slip system in family
              j = j + 1_pInt                                                                        ! could be taken from state by now!
 
              !* Resolved shear stress on slip system
              tau = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))
              !* Stress ratios
              if((abs(tau)-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
              !* Stress ratios
                stressRatio = ((abs(tau)-state(ph)%threshold_stress_slip(j,of))/&
                                (prm%SolidSolutionStrength+&
                                 prm%tau_peierlsPerSlipFamily(f)))
                StressRatio_p       = stressRatio** prm%p(f)
                StressRatio_pminus1 = stressRatio**(prm%p(f)-1.0_pReal)
              !* Boltzmann ratio
                BoltzmannRatio = prm%Qedge(j)/(kB*Temperature)
              !* Initial shear rates
                DotGamma0 = &
                  state(instance)%rhoEdge(j,of)*prm%burgers_slip(j)* &
                  prm%v0(j)
 
              !* Shear rates due to slip
                plastic_dislotwin_postResults(c+j) = &
                  DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**&
                               prm%q(f))*sign(1.0_pReal,tau)
              else
                plastic_dislotwin_postResults(c+j) = 0.0_pReal
              endif
              
        enddo ; enddo
        c = c + prm%totalNslip
      case (accumulated_shear_slip_ID)
       plastic_dislotwin_postResults(c+1_pInt:c+prm%totalNslip) = &
                      state(instance)%accshear_slip(1_pInt:prm%totalNslip,of)
        c = c + prm%totalNslip
      case (mfp_slip_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+prm%totalNslip) =&
                      state(instance)%mfp_slip(1_pInt:prm%totalNslip,of)
        c = c + prm%totalNslip
      case (resolved_stress_slip_ID)
        j = 0_pInt
        do f = 1_pInt,size(prm%Nslip,1)                                                        ! loop over all slip families
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
           do i = 1_pInt,prm%Nslip(f)                                        ! process each (active) slip system in family
              j = j + 1_pInt
              plastic_dislotwin_postResults(c+j) =&
                                dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))
        enddo; enddo
        c = c + prm%totalNslip
      case (threshold_stress_slip_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+prm%totalNslip) = &
                                state(instance)%threshold_stress_slip(1_pInt:prm%totalNslip,of)
        c = c + prm%totalNslip
      case (edge_dipole_distance_ID)
        j = 0_pInt
        do f = 1_pInt,size(prm%Nslip,1)                                                        ! loop over all slip families
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
           do i = 1_pInt,prm%Nslip(f)                                         ! process each (active) slip system in family
              j = j + 1_pInt
              plastic_dislotwin_postResults(c+j) = &
                (3.0_pReal*lattice_mu(ph)*prm%burgers_slip(j))/&
                (16.0_pReal*pi*abs(dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))))
              plastic_dislotwin_postResults(c+j)=min(plastic_dislotwin_postResults(c+j),&
                                                            state(instance)%mfp_slip(j,of))
 !            plastic_dislotwin_postResults(c+j)=max(plastic_dislotwin_postResults(c+j),&
 !                                                            plasticState(ph)%state(4*ns+2*nt+2*nr+j, of))
        enddo; enddo
        c = c + prm%totalNslip
       case (resolved_stress_shearband_ID)
         do j = 1_pInt,6_pInt                                                                       ! loop over all shearband families
            plastic_dislotwin_postResults(c+j) = dot_product(Tstar_v, &
                                                       sbSv(1:6,j,ipc,ip,el))
         enddo
         c = c + 6_pInt
       case (shear_rate_shearband_ID)
         do j = 1_pInt,6_pInt                                                                       ! loop over all shearbands
              !* Resolved shear stress on shearband system
              tau = dot_product(Tstar_v,sbSv(1:6,j,ipc,ip,el))
              !* Stress ratios
              if (abs(tau) < tol_math_check) then
                StressRatio_p = 0.0_pReal
                StressRatio_pminus1 = 0.0_pReal
              else
                StressRatio_p = (abs(tau)/prm%sbResistance)**&
                                 prm%pShearBand
                StressRatio_pminus1 = (abs(tau)/prm%sbResistance)**&
                                (prm%pShearBand-1.0_pReal)
              endif
              !* Boltzmann ratio
              BoltzmannRatio = prm%sbQedge/(kB*Temperature)
              !* Initial shear rates
              DotGamma0 = prm%sbVelocity
              ! Shear rate due to shear band
              plastic_dislotwin_postResults(c+j) = &
          DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**prm%qShearBand)*&
                sign(1.0_pReal,tau)
         enddo 
        c = c + 6_pInt
      case (twin_fraction_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+prm%totalNtwin) = state(instance)%twinFraction(1_pInt:prm%totalNtwin,of)
        c = c + prm%totalNtwin
      case (shear_rate_twin_ID)
        if (prm%totalNtwin > 0_pInt) then
        
          j = 0_pInt
          do f = 1_pInt,size(prm%Nslip,1)
             index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                               ! at which index starts my family
             do i = 1_pInt,prm%Nslip(f)                                      ! process each (active) slip system in family
                j = j + 1_pInt
 
               !* Resolved shear stress on slip system
               tau = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))
               !* Stress ratios
               if((abs(tau)-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
               !* Stress ratios
                 StressRatio_p = ((abs(tau)-state(instance)%threshold_stress_slip(j,of))/&
                                 (prm%SolidSolutionStrength+&
                                  prm%tau_peierlsPerSlipFamily(f)))&
                                                    **prm%p(f)
                 StressRatio_pminus1 = ((abs(tau)-state(instance)%threshold_stress_slip(j,of))/&
                                       (prm%SolidSolutionStrength+&
                                        prm%tau_peierlsPerSlipFamily(f)))&
                                                    **(prm%p(f)-1.0_pReal)
               !* Boltzmann ratio
                 BoltzmannRatio = prm%Qedge(j)/(kB*Temperature)
               !* Initial shear rates
                 DotGamma0 = &
                   state(instance)%rhoEdge(j,of)*prm%burgers_slip(j)* &
                   prm%v0(j)
 
               !* Shear rates due to slip
                 gdot_slip(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**&
                                prm%q(f))*sign(1.0_pReal,tau)
               else
                 gdot_slip(j) = 0.0_pReal
               endif
          enddo;enddo

          j = 0_pInt
          do f = 1_pInt,size(prm%Ntwin,1)
            index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                ! at which index starts my family
            do i = 1,prm%Ntwin(f)                                            ! process each (active) twin system in family
              j = j + 1_pInt

              tau = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,ph))


              !* Shear rates due to twin
              if ( tau > 0.0_pReal ) then
                select case(lattice_structure(ph))
                  case (LATTICE_fcc_ID)
                  s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
                  s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
                  if (tau < tau_r_twin(j,instance)) then
                    Ndot0_twin=(abs(gdot_slip(s1))*(state(instance)%rhoEdge(s2,of)+state(instance)%rhoEdgeDip(s2,of))+&
                           abs(gdot_slip(s2))*(state(instance)%rhoEdge(s1,of)+state(instance)%rhoEdgeDip(s1,of)))/&
                          (prm%L0_twin*&
                           prm%burgers_slip(j))*&
                          (1.0_pReal-exp(-prm%VcrossSlip/(kB*Temperature)*&
                          (tau_r_twin(j,instance)-tau)))
                  else
                    Ndot0_twin=0.0_pReal
                  end if
                  case default
                    Ndot0_twin=prm%Ndot0_twin(j)
                end select
                StressRatio_r = (state(instance)%threshold_stress_twin(j,of)/tau) &
                                                   **prm%r(j)
                plastic_dislotwin_postResults(c+j) = &
                  (prm%MaxTwinFraction-sumf)*lattice_shearTwin(index_myFamily+i,ph)*&
                  state(instance)%twinVolume(j,of)*Ndot0_twin*exp(-StressRatio_r)
              endif
 
          enddo ; enddo
        endif
        c = c + prm%totalNtwin
      case (accumulated_shear_twin_ID)
       plastic_dislotwin_postResults(c+1_pInt:c+prm%totalNtwin) = state(instance)%accshear_twin(1_pInt:prm%totalNtwin,of)
        c = c + prm%totalNtwin     
      case (mfp_twin_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+prm%totalNtwin) = state(instance)%mfp_twin(1_pInt:prm%totalNtwin,of)
        c = c + prm%totalNtwin
      case (resolved_stress_twin_ID)
        if (prm%totalNtwin > 0_pInt) then
          j = 0_pInt
          do f = 1_pInt,size(prm%Ntwin,1)
            index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                ! at which index starts my family
            do i = 1_pInt,prm%Ntwin(f)                                  ! process each (active) slip system in family
              j = j + 1_pInt
              plastic_dislotwin_postResults(c+j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,ph))
          enddo; enddo
        endif
        c = c + prm%totalNtwin
      case (threshold_stress_twin_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+prm%totalNtwin) = state(instance)%threshold_stress_twin(1_pInt:prm%totalNtwin,of)
        c = c + prm%totalNtwin
      case (stress_exponent_ID)
        j = 0_pInt
        do f = 1_pInt,size(prm%Nslip,1)
          index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                  ! at which index starts my family
          do i = 1_pInt,prm%Nslip(f)                                    ! process each (active) slip system in family
             j = j + 1_pInt

             !* Resolved shear stress on slip system
             tau = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))
             if((abs(tau)-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
             !* Stress ratios
               StressRatio_p = ((abs(tau)-state(instance)%threshold_stress_slip(j,of))/&
                               (prm%SolidSolutionStrength+&
                                prm%tau_peierlsPerSlipFamily(f)))&
                                                  **prm%p(f)
               StressRatio_pminus1 = ((abs(tau)-state(instance)%threshold_stress_slip(j,of))/&
                                     (prm%SolidSolutionStrength+&
                                      prm%tau_peierlsPerSlipFamily(f)))&
                                                  **(prm%p(f)-1.0_pReal)
             !* Boltzmann ratio
               BoltzmannRatio = prm%Qedge(j)/(kB*Temperature)
             !* Initial shear rates
               DotGamma0 = &
                 state(instance)%rhoEdge(j,of)*prm%burgers_slip(j)* &
                 prm%v0(j)

             !* Shear rates due to slip
               gdot_slip(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**&
                            prm%q(f))*sign(1.0_pReal,tau)

             !* Derivatives of shear rates
               dgdot_dtauslip = &
                 abs(gdot_slip(j))*BoltzmannRatio*prm%p(f)&
                 *prm%q(f)/&
                 (prm%SolidSolutionStrength+&
                  prm%tau_peierlsPerSlipFamily(f))*&
                 StressRatio_pminus1*(1-StressRatio_p)**(prm%q(f)-1.0_pReal)
             
             else
               gdot_slip(j) = 0.0_pReal
               dgdot_dtauslip = 0.0_pReal
             endif

             !* Stress exponent
             plastic_dislotwin_postResults(c+j) = &
               merge(0.0_pReal,(tau/gdot_slip(j))*dgdot_dtauslip,dEq0(gdot_slip(j)))
         enddo ; enddo
         c = c + prm%totalNslip
      case (sb_eigenvalues_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+3_pInt) = math_eigenvaluesSym33(math_Mandel6to33(Tstar_v))
        c = c + 3_pInt
      case (sb_eigenvectors_ID)
        call math_eigenValuesVectorsSym33(math_Mandel6to33(Tstar_v),eigValues,eigVectors)
        plastic_dislotwin_postResults(c+1_pInt:c+9_pInt) = reshape(eigVectors,[9])
        c = c + 9_pInt
      case (stress_trans_fraction_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+prm%totalNtrans) = &
          state(instance)%stressTransFraction(1_pInt:prm%totalNtrans,of)
        c = c + prm%totalNtrans
      case (strain_trans_fraction_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+prm%totalNtrans) = &
          state(instance)%strainTransFraction(1_pInt:prm%totalNtrans,of)
        c = c + prm%totalNtrans
      case (trans_fraction_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+prm%totalNtrans) = &
          state(instance)%stressTransFraction(1_pInt:prm%totalNtrans,of) + &
          state(instance)%strainTransFraction(1_pInt:prm%totalNtrans,of)
        c = c + prm%totalNtrans
    end select
 enddo
 end associate
end function plastic_dislotwin_postResults

end module plastic_dislotwin
