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
   pReal

 implicit none
 private
 integer,                       dimension(:,:),  allocatable, target, public :: &
   plastic_dislotwin_sizePostResult                                                                 !< size of each post result output
 character(len=64),             dimension(:,:),  allocatable, target, public :: &
   plastic_dislotwin_output                                                                         !< name of each post result output

 real(pReal),                                    parameter,           private :: &
   kB = 1.38e-23_pReal                                                                              !< Boltzmann constant in J/Kelvin

 enum, bind(c)
   enumerator :: &
     undefined_ID, &
     rho_mob_ID, &
     rho_dip_ID, &
     gamma_dot_sl_ID, &
     gamma_sl_ID, &
     mfp_slip_ID, &
     resolved_stress_slip_ID, &
     threshold_stress_slip_ID, &
     edge_dipole_distance_ID, &
     f_tw_ID, &
     mfp_twin_ID, &
     resolved_stress_twin_ID, &
     threshold_stress_twin_ID, &
     strain_trans_fraction_ID
 end enum

 type, private :: tParameters
   real(pReal) :: &
     mu, &
     nu, &
     D0, &                                                                                          !< prefactor for self-diffusion coefficient
     Qsd, &                                                                                         !< activation energy for dislocation climb
     GrainSize, &                                                                                   !<grain size
     pShearBand, &                                                                                  !< p-exponent in shear band velocity
     qShearBand, &                                                                                  !< q-exponent in shear band velocity
     CEdgeDipMinDistance, &                                                                         !<
     Cmfptwin, &                                                                                    !<
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
     aTol_rho, &                                                                                    !< absolute tolerance for integration of dislocation density
     aTol_f_tw, &                                                                                   !< absolute tolerance for integration of twin volume fraction
     aTol_f_tr, &                                                                                   !< absolute tolerance for integration of trans volume fraction
     deltaG, &                                                                                      !< Free energy difference between austensite and martensite
     Cmfptrans, &                                                                                   !<
     transStackHeight                                                                               !< Stack height of hex nucleus 
   real(pReal),                  dimension(:),     allocatable :: & 
     rho_mob_0, &                                                                                   !< initial unipolar dislocation density per slip system
     rho_dip_0, &                                                                                   !< initial dipole dislocation density per slip system
     b_sl, &                                                                                        !< absolute length of burgers vector [m] for each slip system
     b_tw, &                                                                                        !< absolute length of burgers vector [m] for each twin system
     burgers_trans, &                                                                               !< absolute length of burgers vector [m] for each transformation system
     Qedge,&                                                                                        !< activation energy for glide [J] for each slip system
     v0, &                                                                                          !< dislocation velocity prefactor [m/s] for each slip system
     tau_peierls,&                                                                                  !< Peierls stress [Pa] for each slip system
     Ndot0_twin, &                                                                                  !< twin nucleation rate [1/m³s] for each twin system
     Ndot0_trans, &                                                                                 !< trans nucleation rate [1/m³s] for each trans system
     twinsize, &                                                                                    !< twin thickness [m] for each twin system
     CLambdaSlip, &                                                                                 !< Adj. parameter for distance between 2 forest dislocations for each slip system
     atomicVolume, &
     lamellarsize, &                                                                                !< martensite lamellar thickness [m] for each trans system and instance
     p, &                                                                                           !< p-exponent in glide velocity
     q, &                                                                                           !< q-exponent in glide velocity
     r, &                                                                                           !< r-exponent in twin nucleation rate
     s, &                                                                                           !< s-exponent in trans nucleation rate
     shear_twin, &                                                                                  !< characteristic shear for twins
     B                                                                                              !< drag coefficient
   real(pReal),                  dimension(:,:),   allocatable :: & 
     h_sl_sl, &                                                                                     !< 
     h_sl_tw, &                                                                                     !<
     h_tw_tw, &                                                                        !< 
     interaction_SlipTrans, &                                                                       !< 
     interaction_TransTrans                                                                         !< 
   integer,                      dimension(:,:),   allocatable :: & 
     fcc_twinNucleationSlipPair                                                                     ! ToDo: Better name? Is also use for trans
   real(pReal),                  dimension(:,:),   allocatable :: & 
     forestProjection, &
     C66
   real(pReal),                  dimension(:,:,:), allocatable :: &
     Schmid_trans, &
     Schmid_slip, &
     Schmid_twin, &
     C66_twin, &
     C66_trans
   integer :: & 
     totalNslip, &                                                                                   !< total number of active slip system
     totalNtwin, &                                                                                   !< total number of active twin system
     totalNtrans                                                                                     !< total number of active transformation system 
   integer,                      dimension(:),     allocatable :: & 
     N_sl, &                                                                                         !< number of active slip systems for each family
     N_tw, &                                                                                         !< number of active twin systems for each family
     Ntrans                                                                                          !< number of active transformation systems for each family
   integer(kind(undefined_ID)),  dimension(:),     allocatable :: &
     outputID                                                                                        !< ID of each post result output
   logical :: &
     fccTwinTransNucleation, &                                                                       !< twinning and transformation models are for fcc
     dipoleFormation                                                                                 !< flag indicating consideration of dipole formation
 end type                                                                                            !< container type for internal constitutive parameters

 type, private :: tDislotwinState
   real(pReal),                  dimension(:,:),   pointer :: &
     rhoEdge, &
     rhoEdgeDip, &
     accshear_slip, &
     twinFraction, &
     strainTransFraction
 end type tDislotwinState

 type, private :: tDislotwinMicrostructure
   real(pReal),                  dimension(:,:),   allocatable :: &
     invLambdaSlip, &
     invLambdaSlipTwin, &
     invLambdaSlipTrans, &
     invLambdaTwin, &
     invLambdaTrans, &
     mfp_slip, &
     mfp_twin, &
     mfp_trans, &
     tau_pass, &
     threshold_stress_twin, &
     threshold_stress_trans, &
     twinVolume, &
     martensiteVolume, &
     tau_r_twin, &                                                                                  !< stress to bring partials close together (twin)
     tau_r_trans                                                                                    !< stress to bring partials close together (trans)
 end type tDislotwinMicrostructure

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
 type(tParameters),              allocatable, dimension(:), private :: param
 type(tDislotwinState),          allocatable, dimension(:), private :: &
   dotState, &
   state
 type(tDislotwinMicrostructure), allocatable, dimension(:), private :: microstructure

 public :: &
   plastic_dislotwin_init, &
   plastic_dislotwin_homogenizedC, &
   plastic_dislotwin_dependentState, &
   plastic_dislotwin_LpAndItsTangent, &
   plastic_dislotwin_dotState, &
   plastic_dislotwin_postResults, &
   plastic_dislotwin_results
 private :: &
   kinetics_slip, &
   kinetics_twin, &
   kinetics_trans

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_init
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
   math_expand,&
   PI
 use IO, only: &
   IO_error
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
   config_phase
 use lattice

 implicit none
 integer :: &
   Ninstance, &
   p, i, &
   NipcMyPhase, outputSize, &
   sizeState, sizeDotState, &
   startIndex, endIndex

 integer,               dimension(0), parameter :: emptyIntArray    = [integer::]
 real(pReal),           dimension(0), parameter :: emptyRealArray   = [real(pReal)::]
 character(len=65536),  dimension(0), parameter :: emptyStringArray = [character(len=65536)::]

 integer(kind(undefined_ID)) :: &
   outputID

 character(len=pStringLen) :: &
   extmsg = ''
 character(len=65536), dimension(:), allocatable :: &
   outputs

 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_DISLOTWIN_label//' init  -+>>>'

 write(6,'(/,a)')   ' Ma and Roters, Acta Materialia 52(12):3603–3612, 2004'
 write(6,'(a)')     ' https://doi.org/10.1016/j.actamat.2004.04.012'

 write(6,'(/,a)')   ' Roters et al., Computational Materials Science 39:91–95, 2007'
 write(6,'(a)')     ' https://doi.org/10.1016/j.commatsci.2006.04.014'

 write(6,'(/,a)')   ' Wong et al., Acta Materialia 118:140–151, 2016'
 write(6,'(a,/)')   ' https://doi.org/10.1016/j.actamat.2016.07.032'

 Ninstance = count(phase_plasticity == PLASTICITY_DISLOTWIN_ID)

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
   write(6,'(a16,1x,i5,/)') '# instances:',Ninstance
 
 allocate(plastic_dislotwin_sizePostResult(maxval(phase_Noutput),Ninstance),source=0)
 allocate(plastic_dislotwin_output(maxval(phase_Noutput),Ninstance))
          plastic_dislotwin_output = ''

 allocate(param(Ninstance))
 allocate(state(Ninstance))
 allocate(dotState(Ninstance))
 allocate(microstructure(Ninstance))

 do p = 1, size(phase_plasticity)
   if (phase_plasticity(p) /= PLASTICITY_DISLOTWIN_ID) cycle
   associate(prm => param(phase_plasticityInstance(p)), &
             dot => dotState(phase_plasticityInstance(p)), &
             stt => state(phase_plasticityInstance(p)), &
             dst => microstructure(phase_plasticityInstance(p)), &
             config   => config_phase(p))

   prm%aTol_rho  = config%getFloat('atol_rho',       defaultVal=0.0_pReal)
   prm%aTol_f_tw = config%getFloat('atol_twinfrac',  defaultVal=0.0_pReal)
   prm%aTol_f_tr = config%getFloat('atol_transfrac', defaultVal=0.0_pReal)

   ! This data is read in already in lattice
   prm%mu = lattice_mu(p)
   prm%nu = lattice_nu(p)
   prm%C66 = lattice_C66(1:6,1:6,p)


!--------------------------------------------------------------------------------------------------
! slip related parameters
   prm%N_sl      = config%getInts('nslip',defaultVal=emptyIntArray)
   prm%totalNslip = sum(prm%N_sl)
   slipActive: if (prm%totalNslip > 0) then
     prm%Schmid_slip          = lattice_SchmidMatrix_slip(prm%N_sl,config%getString('lattice_structure'),&
                                                          config%getFloat('c/a',defaultVal=0.0_pReal))
     prm%h_sl_sl = lattice_interaction_SlipBySlip(prm%N_sl, &
                                                  config%getFloats('interaction_slipslip'), &
                                                  config%getString('lattice_structure'))
     prm%forestProjection     = lattice_forestProjection (prm%N_sl,config%getString('lattice_structure'),&
                                                          config%getFloat('c/a',defaultVal=0.0_pReal))

     prm%fccTwinTransNucleation = merge(.true., .false., lattice_structure(p) == LATTICE_FCC_ID) &
                                .and. (prm%N_sl(1) == 12)
     if(prm%fccTwinTransNucleation) &
       prm%fcc_twinNucleationSlipPair = lattice_fcc_twinNucleationSlipPair

     prm%rho_mob_0            = config%getFloats('rhoedge0',   requiredSize=size(prm%N_sl))
     prm%rho_dip_0            = config%getFloats('rhoedgedip0',requiredSize=size(prm%N_sl))
     prm%v0                   = config%getFloats('v0',         requiredSize=size(prm%N_sl))
     prm%b_sl                 = config%getFloats('slipburgers',requiredSize=size(prm%N_sl))
     prm%Qedge                = config%getFloats('qedge',      requiredSize=size(prm%N_sl))  !ToDo: rename (ask Karo)
     prm%CLambdaSlip          = config%getFloats('clambdaslip',requiredSize=size(prm%N_sl))
     prm%p                    = config%getFloats('p_slip',     requiredSize=size(prm%N_sl))
     prm%q                    = config%getFloats('q_slip',     requiredSize=size(prm%N_sl))
     prm%B                    = config%getFloats('b',          requiredSize=size(prm%N_sl), &
                                                          defaultVal=[(0.0_pReal, i=1,size(prm%N_sl))])
     prm%tau_peierls          = config%getFloats('tau_peierls',requiredSize=size(prm%N_sl), &
                                                          defaultVal=[(0.0_pReal, i=1,size(prm%N_sl))]) ! Deprecated

     prm%CEdgeDipMinDistance  = config%getFloat('cedgedipmindistance')
     prm%D0                   = config%getFloat('d0')
     prm%Qsd                  = config%getFloat('qsd')
     prm%atomicVolume         = config%getFloat('catomicvolume') * prm%b_sl**3.0_pReal

     ! expand: family => system
     prm%rho_mob_0    = math_expand(prm%rho_mob_0,   prm%N_sl)
     prm%rho_dip_0    = math_expand(prm%rho_dip_0,   prm%N_sl)
     prm%v0           = math_expand(prm%v0,          prm%N_sl)
     prm%b_sl         = math_expand(prm%b_sl,prm%N_sl)
     prm%Qedge        = math_expand(prm%Qedge,       prm%N_sl)
     prm%CLambdaSlip  = math_expand(prm%CLambdaSlip, prm%N_sl)
     prm%p            = math_expand(prm%p,           prm%N_sl)
     prm%q            = math_expand(prm%q,           prm%N_sl)
     prm%B            = math_expand(prm%B,           prm%N_sl)
     prm%tau_peierls  = math_expand(prm%tau_peierls, prm%N_sl)
     prm%atomicVolume = math_expand(prm%atomicVolume,prm%N_sl)                                   

     ! sanity checks
     if (    prm%D0           <= 0.0_pReal)          extmsg = trim(extmsg)//' D0'
     if (    prm%Qsd          <= 0.0_pReal)          extmsg = trim(extmsg)//' Qsd'
     if (any(prm%rho_mob_0    <  0.0_pReal))         extmsg = trim(extmsg)//' rho_mob_0'
     if (any(prm%rho_dip_0    <  0.0_pReal))         extmsg = trim(extmsg)//' rho_dip_0'
     if (any(prm%v0           <  0.0_pReal))         extmsg = trim(extmsg)//' v0'
     if (any(prm%b_sl         <= 0.0_pReal))         extmsg = trim(extmsg)//' b_sl'
     if (any(prm%Qedge        <= 0.0_pReal))         extmsg = trim(extmsg)//' Qedge'
     if (any(prm%CLambdaSlip  <= 0.0_pReal))         extmsg = trim(extmsg)//' CLambdaSlip'
     if (any(prm%B            <  0.0_pReal))         extmsg = trim(extmsg)//' B'
     if (any(prm%tau_peierls  <  0.0_pReal))         extmsg = trim(extmsg)//' tau_peierls'
     if (any(prm%p<=0.0_pReal .or. prm%p>1.0_pReal)) extmsg = trim(extmsg)//' p'
     if (any(prm%q< 1.0_pReal .or. prm%q>2.0_pReal)) extmsg = trim(extmsg)//' q'

   else slipActive
     allocate(prm%b_sl(0))
   endif slipActive

!--------------------------------------------------------------------------------------------------
! twin related parameters
   prm%N_tw      = config%getInts('ntwin', defaultVal=emptyIntArray)
   prm%totalNtwin = sum(prm%N_tw)
   if (prm%totalNtwin > 0) then
     prm%Schmid_twin  = lattice_SchmidMatrix_twin(prm%N_tw,config%getString('lattice_structure'),&
                                                  config%getFloat('c/a',defaultVal=0.0_pReal))
     prm%h_tw_tw      = lattice_interaction_TwinByTwin(prm%N_tw,&
                                                       config%getFloats('interaction_twintwin'), &
                                                       config%getString('lattice_structure'))

     prm%b_tw         = config%getFloats('twinburgers',  requiredSize=size(prm%N_tw))
     prm%twinsize     = config%getFloats('twinsize',     requiredSize=size(prm%N_tw))
     prm%r            = config%getFloats('r_twin',       requiredSize=size(prm%N_tw))

     prm%xc_twin         = config%getFloat('xc_twin')
     prm%L0_twin         = config%getFloat('l0_twin')
     prm%Cmfptwin        = config%getFloat('cmfptwin',       defaultVal=0.0_pReal) ! ToDo: How to handle that???

     prm%shear_twin      = lattice_characteristicShear_Twin(prm%N_tw,config%getString('lattice_structure'),&
                                                            config%getFloat('c/a',defaultVal=0.0_pReal))

     prm%C66_twin        = lattice_C66_twin(prm%N_tw,prm%C66,config%getString('lattice_structure'),&
                                            config%getFloat('c/a',defaultVal=0.0_pReal))

     if (.not. prm%fccTwinTransNucleation) then
       prm%Ndot0_twin = config%getFloats('ndot0_twin') 
       prm%Ndot0_twin = math_expand(prm%Ndot0_twin,prm%N_tw)
     endif

     ! expand: family => system
     prm%b_tw         = math_expand(prm%b_tw,prm%N_tw)
     prm%twinsize     = math_expand(prm%twinsize,prm%N_tw)
     prm%r            = math_expand(prm%r,prm%N_tw)
     
   else
     allocate(prm%twinsize(0))
     allocate(prm%b_tw(0))
     allocate(prm%r(0))
   endif
  
!--------------------------------------------------------------------------------------------------
! transformation related parameters
   prm%Ntrans      = config%getInts('ntrans', defaultVal=emptyIntArray)
   prm%totalNtrans = sum(prm%Ntrans)
   if (prm%totalNtrans > 0) then
     prm%burgers_trans = config%getFloats('transburgers')
     prm%burgers_trans = math_expand(prm%burgers_trans,prm%Ntrans)
     
     prm%transStackHeight = config%getFloat('transstackheight', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%Cmfptrans        = config%getFloat('cmfptrans', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%deltaG           = config%getFloat('deltag')
     prm%xc_trans         = config%getFloat('xc_trans', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%L0_trans         = config%getFloat('l0_trans')

     prm%interaction_TransTrans = lattice_interaction_TransByTrans(prm%Ntrans,&
                                                                   config%getFloats('interaction_transtrans'), &
                                                                   config%getString('lattice_structure'))
                                                             
     prm%C66_trans        = lattice_C66_trans(prm%Ntrans,prm%C66, &
                                  config%getString('trans_lattice_structure'), &
                                  0.0_pReal, &
                                  config%getFloat('a_bcc', defaultVal=0.0_pReal), &
                                  config%getFloat('a_fcc', defaultVal=0.0_pReal))
                                  
      prm%Schmid_trans        = lattice_SchmidMatrix_trans(prm%Ntrans, &
                                  config%getString('trans_lattice_structure'), &
                                  0.0_pReal, &
                                  config%getFloat('a_bcc', defaultVal=0.0_pReal), &
                                  config%getFloat('a_fcc', defaultVal=0.0_pReal))
                                                 
     if (lattice_structure(p) /= LATTICE_fcc_ID) then
        prm%Ndot0_trans = config%getFloats('ndot0_trans')
        prm%Ndot0_trans = math_expand(prm%Ndot0_trans,prm%Ntrans)
     endif
     prm%lamellarsize = config%getFloats('lamellarsize')
     prm%lamellarsize = math_expand(prm%lamellarsize,prm%Ntrans)
     prm%s = config%getFloats('s_trans',defaultVal=[0.0_pReal])
     prm%s = math_expand(prm%s,prm%Ntrans)
   else
     allocate(prm%lamellarsize(0))
     allocate(prm%burgers_trans(0))
   endif
   
   if (sum(prm%N_tw) > 0  .or. prm%totalNtrans > 0) then
     prm%SFE_0K     = config%getFloat('sfe_0k')
     prm%dSFE_dT    = config%getFloat('dsfe_dt')
     prm%VcrossSlip = config%getFloat('vcrossslip')
   endif
   
   if (prm%totalNslip > 0 .and. prm%totalNtwin > 0) then
     prm%h_sl_tw = lattice_interaction_SlipByTwin(prm%N_sl,prm%N_tw,&
                                                  config%getFloats('interaction_sliptwin'), &
                                                  config%getString('lattice_structure'))
     if (prm%fccTwinTransNucleation .and. prm%totalNtwin > 12) write(6,*) 'mist' ! ToDo: implement better test. The model will fail also if N_tw is [6,6]
   endif    

   if (prm%totalNslip > 0 .and. prm%totalNtrans > 0) then  
     prm%interaction_SlipTrans = lattice_interaction_SlipByTrans(prm%N_sl,prm%Ntrans,&
                                                                 config%getFloats('interaction_sliptrans'), &
                                                                 config%getString('lattice_structure')) 
     if (prm%fccTwinTransNucleation .and. prm%totalNtrans > 12) write(6,*) 'mist' ! ToDo: implement better test. The model will fail also if ntrans is [6,6]
   endif  
  
!--------------------------------------------------------------------------------------------------
! shearband related parameters
   prm%sbVelocity = config%getFloat('shearbandvelocity',defaultVal=0.0_pReal)
   if (prm%sbVelocity > 0.0_pReal) then  
     prm%sbResistance = config%getFloat('shearbandresistance')
     prm%sbQedge      = config%getFloat('qedgepersbsystem')
     prm%pShearBand   = config%getFloat('p_shearband')
     prm%qShearBand   = config%getFloat('q_shearband')
     
     ! sanity checks
     if (prm%sbResistance  <  0.0_pReal) extmsg = trim(extmsg)//' shearbandresistance'
     if (prm%sbQedge       <  0.0_pReal) extmsg = trim(extmsg)//' qedgepersbsystem'
     if (prm%pShearBand    <= 0.0_pReal) extmsg = trim(extmsg)//' p_shearband'
     if (prm%qShearBand    <= 0.0_pReal) extmsg = trim(extmsg)//' q_shearband'
   endif



   prm%GrainSize             = config%getFloat('grainsize')
   prm%SolidSolutionStrength = config%getFloat('solidsolutionstrength')                              ! Deprecated

   if (config%keyExists('dipoleformationfactor')) call IO_error(1,ext_msg='use /nodipoleformation/')
   prm%dipoleformation = .not. config%keyExists('/nodipoleformation/')
   

       !if (Ndot0PerTwinFamily(f,p) < 0.0_pReal) &
        ! call IO_error(211,el=p,ext_msg='ndot0_twin ('//PLASTICITY_DISLOTWIN_label//')')

   if (any(prm%atomicVolume <= 0.0_pReal)) &
     call IO_error(211,el=p,ext_msg='cAtomicVolume ('//PLASTICITY_DISLOTWIN_label//')')
   if (prm%totalNtwin > 0) then
     if (prm%aTol_rho <= 0.0_pReal) &
       call IO_error(211,el=p,ext_msg='aTol_rho ('//PLASTICITY_DISLOTWIN_label//')')   
     if (prm%aTol_f_tw <= 0.0_pReal) &
       call IO_error(211,el=p,ext_msg='aTol_f_tw ('//PLASTICITY_DISLOTWIN_label//')')
   endif
   if (prm%totalNtrans > 0) then
     if (prm%aTol_f_tr <= 0.0_pReal) &
       call IO_error(211,el=p,ext_msg='aTol_f_tr ('//PLASTICITY_DISLOTWIN_label//')')
   endif
 
   outputs = config%getStrings('(output)', defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   do i= 1, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
       case ('edge_density')
         outputID = merge(rho_mob_ID,undefined_ID,prm%totalNslip > 0)
         outputSize = prm%totalNslip
       case ('dipole_density')
         outputID = merge(rho_dip_ID,undefined_ID,prm%totalNslip > 0)
         outputSize = prm%totalNslip
       case ('shear_rate_slip','shearrate_slip')
         outputID = merge(gamma_dot_sl_ID,undefined_ID,prm%totalNslip > 0)
         outputSize = prm%totalNslip
       case ('accumulated_shear_slip')
         outputID = merge(gamma_sl_ID,undefined_ID,prm%totalNslip > 0)
         outputSize = prm%totalNslip
       case ('mfp_slip')
         outputID = merge(mfp_slip_ID,undefined_ID,prm%totalNslip > 0)
         outputSize = prm%totalNslip
       case ('resolved_stress_slip')
         outputID = merge(resolved_stress_slip_ID,undefined_ID,prm%totalNslip > 0)
         outputSize = prm%totalNslip
       case ('threshold_stress_slip')
         outputID= merge(threshold_stress_slip_ID,undefined_ID,prm%totalNslip > 0)
         outputSize = prm%totalNslip

       case ('twin_fraction')
         outputID = merge(f_tw_ID,undefined_ID,prm%totalNtwin >0)
         outputSize = prm%totalNtwin
       case ('mfp_twin')
         outputID = merge(mfp_twin_ID,undefined_ID,prm%totalNtwin >0)
         outputSize = prm%totalNtwin
       case ('resolved_stress_twin')
         outputID = merge(resolved_stress_twin_ID,undefined_ID,prm%totalNtwin >0)
         outputSize = prm%totalNtwin
       case ('threshold_stress_twin')
         outputID = merge(threshold_stress_twin_ID,undefined_ID,prm%totalNtwin >0)
         outputSize = prm%totalNtwin
         
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
   NipcMyPhase = count(material_phase == p)
   sizeDotState = size(['rho         ','rhoDip      ','accshearslip']) * prm%totalNslip &
                + size(['twinFraction'])                               * prm%totalNtwin &
                + size(['strainTransFraction'])                        * prm%totalNtrans
   sizeState = sizeDotState

   call material_allocatePlasticState(p,NipcMyPhase,sizeState,sizeDotState,0, &
                                      prm%totalNslip,prm%totalNtwin,prm%totalNtrans)
   plasticState(p)%sizePostResults = sum(plastic_dislotwin_sizePostResult(:,phase_plasticityInstance(p)))


!--------------------------------------------------------------------------------------------------
! locally defined state aliases and initialization of state0 and aTolState
   startIndex = 1
   endIndex   = prm%totalNslip
   stt%rhoEdge=>plasticState(p)%state(startIndex:endIndex,:)
   stt%rhoEdge= spread(prm%rho_mob_0,2,NipcMyPhase)
   dot%rhoEdge=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTol_rho

   startIndex = endIndex + 1
   endIndex   = endIndex + prm%totalNslip
   stt%rhoEdgeDip=>plasticState(p)%state(startIndex:endIndex,:)
   stt%rhoEdgeDip= spread(prm%rho_dip_0,2,NipcMyPhase)
   dot%rhoEdgeDip=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTol_rho

   startIndex = endIndex + 1
   endIndex   = endIndex + prm%totalNslip
   stt%accshear_slip=>plasticState(p)%state(startIndex:endIndex,:)
   dot%accshear_slip=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = 1.0e6_pReal  !ToDo: better make optional parameter
   ! global alias
   plasticState(p)%slipRate        => plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%accumulatedSlip => plasticState(p)%state(startIndex:endIndex,:)
   
   startIndex = endIndex + 1
   endIndex   = endIndex + prm%totalNtwin
   stt%twinFraction=>plasticState(p)%state(startIndex:endIndex,:)
   dot%twinFraction=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTol_f_tw
   
   startIndex = endIndex + 1
   endIndex   = endIndex + prm%totalNtrans
   stt%strainTransFraction=>plasticState(p)%state(startIndex:endIndex,:)
   dot%strainTransFraction=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTol_f_tr

   allocate(dst%invLambdaSlip         (prm%totalNslip, NipcMyPhase),source=0.0_pReal)
   allocate(dst%invLambdaSlipTwin     (prm%totalNslip, NipcMyPhase),source=0.0_pReal)
   allocate(dst%invLambdaSlipTrans    (prm%totalNslip, NipcMyPhase),source=0.0_pReal)
   allocate(dst%mfp_slip              (prm%totalNslip, NipcMyPhase),source=0.0_pReal)
   allocate(dst%tau_pass              (prm%totalNslip, NipcMyPhase),source=0.0_pReal)

   allocate(dst%invLambdaTwin         (prm%totalNtwin, NipcMyPhase),source=0.0_pReal)
   allocate(dst%mfp_twin              (prm%totalNtwin, NipcMyPhase),source=0.0_pReal)
   allocate(dst%threshold_stress_twin (prm%totalNtwin, NipcMyPhase),source=0.0_pReal)
   allocate(dst%tau_r_twin            (prm%totalNtwin, NipcMyPhase),source=0.0_pReal)
   allocate(dst%twinVolume            (prm%totalNtwin, NipcMyPhase),source=0.0_pReal)

   allocate(dst%invLambdaTrans        (prm%totalNtrans,NipcMyPhase),source=0.0_pReal)
   allocate(dst%mfp_trans             (prm%totalNtrans,NipcMyPhase),source=0.0_pReal)
   allocate(dst%threshold_stress_trans(prm%totalNtrans,NipcMyPhase),source=0.0_pReal)
   allocate(dst%tau_r_trans           (prm%totalNtrans,NipcMyPhase),source=0.0_pReal)
   allocate(dst%martensiteVolume      (prm%totalNtrans,NipcMyPhase),source=0.0_pReal)


   plasticState(p)%state0 = plasticState(p)%state                                                   ! ToDo: this could be done centrally

   end associate

 enddo

end subroutine plastic_dislotwin_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenized elasticity matrix
!--------------------------------------------------------------------------------------------------
function plastic_dislotwin_homogenizedC(ipc,ip,el) result(homogenizedC)
 use material, only: &
  material_phase, &
  phase_plasticityInstance, &
  phasememberAt
 
 implicit none
 real(pReal), dimension(6,6) :: &
   homogenizedC
 integer,     intent(in) :: &
   ipc, &                                                                                          !< component-ID of integration point
   ip, &                                                                                           !< integration point
   el                                                                                              !< element

 integer :: i, &
                  of
 real(pReal) :: f_unrotated

 of = phasememberAt(ipc,ip,el)
 associate(prm => param(phase_plasticityInstance(material_phase(ipc,ip,el))),&
           stt => state(phase_plasticityInstance(material_phase(ipc,ip,el))))

 f_unrotated = 1.0_pReal &
             - sum(stt%twinFraction(1:prm%totalNtwin,of)) &
             - sum(stt%strainTransFraction(1:prm%totalNtrans,of))

 homogenizedC = f_unrotated * prm%C66
 do i=1,prm%totalNtwin
   homogenizedC = homogenizedC &
                + stt%twinFraction(i,of)*prm%C66_twin(1:6,1:6,i)
 enddo
 do i=1,prm%totalNtrans
   homogenizedC = homogenizedC &
                + stt%strainTransFraction(i,of)*prm%C66_trans(1:6,1:6,i)
 enddo

 end associate
 
end function plastic_dislotwin_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_LpAndItsTangent(Lp,dLp_dMp,Mp,Temperature,instance,of)
 use prec, only: &
   tol_math_check, &
   dNeq0
 use math, only: &
   math_eigenValuesVectorsSym, &
   math_outer, &
   math_symmetric33, &
   math_mul33xx33, &
   math_mul33x3
 
 implicit none
 real(pReal), dimension(3,3),     intent(out) :: Lp
 real(pReal), dimension(3,3,3,3), intent(out) :: dLp_dMp
 real(pReal), dimension(3,3),     intent(in)  :: Mp
 integer,                         intent(in)  :: instance,of
 real(pReal),                     intent(in)  :: Temperature

 integer :: i,k,l,m,n
 real(pReal) :: f_unrotated,StressRatio_p,&
                BoltzmannRatio, &
    dgdot_dtau, &
    tau
 real(pReal), dimension(param(instance)%totalNslip) :: &
    gdot_slip,dgdot_dtau_slip
 real(pReal), dimension(param(instance)%totalNtwin) :: &
    gdot_twin,dgdot_dtau_twin
 real(pReal), dimension(param(instance)%totalNtrans) :: &
    gdot_trans,dgdot_dtau_trans
 real(pReal):: gdot_sb
 real(pReal), dimension(3,3) :: eigVectors, Schmid_shearBand
 real(pReal), dimension(3)   :: eigValues
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

 associate(prm => param(instance), stt => state(instance))

 f_unrotated = 1.0_pReal &
             - sum(stt%twinFraction(1:prm%totalNtwin,of)) &
             - sum(stt%strainTransFraction(1:prm%totalNtrans,of))

 Lp = 0.0_pReal
 dLp_dMp = 0.0_pReal 

 call kinetics_slip(Mp,temperature,instance,of,gdot_slip,dgdot_dtau_slip)
 slipContribution: do i = 1, prm%totalNslip
   Lp = Lp + gdot_slip(i)*prm%Schmid_slip(1:3,1:3,i)
   forall (k=1:3,l=1:3,m=1:3,n=1:3) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + dgdot_dtau_slip(i) * prm%Schmid_slip(k,l,i) * prm%Schmid_slip(m,n,i)
 enddo slipContribution
 
 !ToDo: Why do this before shear banding?
 Lp      = Lp      * f_unrotated
 dLp_dMp = dLp_dMp * f_unrotated
 
 shearBandingContribution: if(dNeq0(prm%sbVelocity)) then

   BoltzmannRatio = prm%sbQedge/(kB*Temperature)
   call math_eigenValuesVectorsSym(Mp,eigValues,eigVectors,error)

   do i = 1,6
     Schmid_shearBand = 0.5_pReal * math_outer(math_mul33x3(eigVectors,sb_sComposition(1:3,i)),&
                                               math_mul33x3(eigVectors,sb_mComposition(1:3,i)))
     tau = math_mul33xx33(Mp,Schmid_shearBand)
   
     significantShearBandStress: if (abs(tau) > tol_math_check) then
       StressRatio_p = (abs(tau)/prm%sbResistance)**prm%pShearBand
       gdot_sb = sign(prm%sbVelocity*exp(-BoltzmannRatio*(1-StressRatio_p)**prm%qShearBand), tau)
       dgdot_dtau = abs(gdot_sb)*BoltzmannRatio* prm%pShearBand*prm%qShearBand/ prm%sbResistance &
                  * (abs(tau)/prm%sbResistance)**(prm%pShearBand-1.0_pReal) &
                  * (1.0_pReal-StressRatio_p)**(prm%qShearBand-1.0_pReal)
 
       Lp = Lp + gdot_sb * Schmid_shearBand
       forall (k=1:3,l=1:3,m=1:3,n=1:3) &
         dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                          + dgdot_dtau * Schmid_shearBand(k,l) * Schmid_shearBand(m,n)
     endif significantShearBandStress
   enddo

 endif shearBandingContribution
 
 call kinetics_twin(Mp,temperature,gdot_slip,instance,of,gdot_twin,dgdot_dtau_twin)
 twinContibution: do i = 1, prm%totalNtwin
   Lp = Lp + gdot_twin(i)*prm%Schmid_twin(1:3,1:3,i) * f_unrotated
   forall (k=1:3,l=1:3,m=1:3,n=1:3) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + dgdot_dtau_twin(i)* prm%Schmid_twin(k,l,i)*prm%Schmid_twin(m,n,i) * f_unrotated
 enddo twinContibution
 
 call kinetics_twin(Mp,temperature,gdot_slip,instance,of,gdot_trans,dgdot_dtau_trans)
 transContibution: do i = 1, prm%totalNtrans
   Lp = Lp + gdot_trans(i)*prm%Schmid_trans(1:3,1:3,i) * f_unrotated
   forall (k=1:3,l=1:3,m=1:3,n=1:3) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + dgdot_dtau_trans(i)* prm%Schmid_trans(k,l,i)*prm%Schmid_trans(m,n,i) * f_unrotated
 enddo transContibution


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
   math_clip, &
   math_mul33xx33, &
   PI 
 use material, only: &
   plasticState

 implicit none
 real(pReal), dimension(3,3),  intent(in):: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                  intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer,                      intent(in) :: &
   instance, &
   of

 integer :: i
 real(pReal) :: f_unrotated,&
             VacancyDiffusion,&
             EdgeDipDistance, ClimbVelocity,DotRhoEdgeDipClimb,DotRhoEdgeDipAnnihilation, &
             DotRhoDipFormation,DotRhoEdgeEdgeAnnihilation, &
            tau
 real(pReal), dimension(param(instance)%totalNslip) :: &
   EdgeDipMinDistance, &
   DotRhoMultiplication, &
   gdot_slip
 real(pReal), dimension(param(instance)%totalNtwin) :: &
   gdot_twin
 real(pReal), dimension(param(instance)%totalNtrans) :: &
   gdot_trans

 associate(prm => param(instance),    stt => state(instance), &
           dot => dotstate(instance), dst => microstructure(instance))

 f_unrotated = 1.0_pReal &
             - sum(stt%twinFraction(1:prm%totalNtwin,of)) &
             - sum(stt%strainTransFraction(1:prm%totalNtrans,of))
 VacancyDiffusion = prm%D0*exp(-prm%Qsd/(kB*Temperature))

 call kinetics_slip(Mp,temperature,instance,of,gdot_slip)
 dot%accshear_slip(:,of) = abs(gdot_slip)
 
 DotRhoMultiplication = abs(gdot_slip)/(prm%b_sl*dst%mfp_slip(:,of))
 EdgeDipMinDistance   = prm%CEdgeDipMinDistance*prm%b_sl
 
 slipState: do i = 1, prm%totalNslip
   tau = math_mul33xx33(Mp,prm%Schmid_slip(1:3,1:3,i))

   significantSlipStress: if (dEq0(tau)) then
     DotRhoDipFormation = 0.0_pReal
     DotRhoEdgeDipClimb = 0.0_pReal
   else significantSlipStress
     EdgeDipDistance = 3.0_pReal*prm%mu*prm%b_sl(i)/(16.0_pReal*PI*abs(tau))
     EdgeDipDistance = math_clip(EdgeDipDistance, right = dst%mfp_slip(i,of))
     EdgeDipDistance = math_clip(EdgeDipDistance, left  = EdgeDipMinDistance(i))

     if (prm%dipoleFormation) then
       DotRhoDipFormation = 2.0_pReal*(EdgeDipDistance-EdgeDipMinDistance(i))/prm%b_sl(i) &
                          * stt%rhoEdge(i,of)*abs(gdot_slip(i))
     else
       DotRhoDipFormation = 0.0_pReal
     endif

     if (dEq0(EdgeDipDistance-EdgeDipMinDistance(i))) then
       DotRhoEdgeDipClimb = 0.0_pReal
     else
       ClimbVelocity = 3.0_pReal*prm%mu*VacancyDiffusion*prm%atomicVolume(i) &
                     / (2.0_pReal*PI*kB*Temperature*(EdgeDipDistance+EdgeDipMinDistance(i)))
       DotRhoEdgeDipClimb = 4.0_pReal*ClimbVelocity*stt%rhoEdgeDip(i,of) &
                          / (EdgeDipDistance-EdgeDipMinDistance(i))
     endif
   endif significantSlipStress
 
   !* Spontaneous annihilation of 2 single edge dislocations
   DotRhoEdgeEdgeAnnihilation = 2.0_pReal*EdgeDipMinDistance(i)/prm%b_sl(i) &
                              * stt%rhoEdge(i,of)*abs(gdot_slip(i))
   !* Spontaneous annihilation of a single edge dislocation with a dipole constituent
   DotRhoEdgeDipAnnihilation = 2.0_pReal*EdgeDipMinDistance(i)/prm%b_sl(i) &
                             * stt%rhoEdgeDip(i,of)*abs(gdot_slip(i))

   dot%rhoEdge(i,of)    = DotRhoMultiplication(i)-DotRhoDipFormation-DotRhoEdgeEdgeAnnihilation
   dot%rhoEdgeDip(i,of) = DotRhoDipFormation-DotRhoEdgeDipAnnihilation-DotRhoEdgeDipClimb
 enddo slipState
 
 call kinetics_twin(Mp,temperature,gdot_slip,instance,of,gdot_twin)
 dot%twinFraction(:,of) = f_unrotated*gdot_twin/prm%shear_twin

 call kinetics_trans(Mp,temperature,gdot_slip,instance,of,gdot_trans)
 dot%twinFraction(:,of) = f_unrotated*gdot_trans

 end associate
 
end subroutine plastic_dislotwin_dotState


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_dependentState(temperature,instance,of)
 use math, only: &
   PI

 implicit none
 integer,       intent(in) :: &
   instance, &
   of
 real(pReal),   intent(in) :: &
   temperature

 integer :: &
   i
 real(pReal) :: &
   sumf_twin,SFE,sumf_trans
 real(pReal), dimension(:), allocatable :: &
   x0, &
   fOverStacksize, &
   ftransOverLamellarSize


 associate(prm => param(instance),&
           stt => state(instance),&
           dst => microstructure(instance))

 sumf_twin  = sum(stt%twinFraction(1:prm%totalNtwin,of))
 sumf_trans = sum(stt%strainTransFraction(1:prm%totalNtrans,of))

 SFE = prm%SFE_0K + prm%dSFE_dT * Temperature
 
 !* rescaled volume fraction for topology
 fOverStacksize         =  stt%twinFraction(1:prm%totalNtwin,of)/prm%twinsize  !ToDo: this is per system
 ftransOverLamellarSize =  sumf_trans/prm%lamellarsize                !ToDo: But this not ...
 !Todo: Physically ok, but naming could be adjusted


 !* 1/mean free distance between 2 forest dislocations seen by a moving dislocation
 forall (i = 1:prm%totalNslip) &
   dst%invLambdaSlip(i,of) = &
     sqrt(dot_product((stt%rhoEdge(1:prm%totalNslip,of)+stt%rhoEdgeDip(1:prm%totalNslip,of)),&
                      prm%forestProjection(1:prm%totalNslip,i)))/prm%CLambdaSlip(i)

 !* 1/mean free distance between 2 twin stacks from different systems seen by a moving dislocation
 if (prm%totalNtwin > 0 .and. prm%totalNslip > 0) &
   dst%invLambdaSlipTwin(1:prm%totalNslip,of) = &
     matmul(transpose(prm%h_sl_tw),fOverStacksize)/(1.0_pReal-sumf_twin)               ! ToDo: Change order and use matmul

 !* 1/mean free distance between 2 twin stacks from different systems seen by a growing twin

  !ToDo: needed? if (prm%totalNtwin > 0) &
 dst%invLambdaTwin(1:prm%totalNtwin,of) = matmul(prm%h_tw_tw,fOverStacksize)/(1.0_pReal-sumf_twin)


 !* 1/mean free distance between 2 martensite lamellar from different systems seen by a moving dislocation
 if (prm%totalNtrans > 0 .and. prm%totalNslip > 0) &
   dst%invLambdaSlipTrans(1:prm%totalNslip,of) = &                                  ! ToDo: does not work if Ntrans is not 12
      matmul(transpose(prm%interaction_SlipTrans),ftransOverLamellarSize)/(1.0_pReal-sumf_trans)    ! ToDo: Transpose needed

 !* 1/mean free distance between 2 martensite stacks from different systems seen by a growing martensite (1/lambda_trans)
 !ToDo: needed? if (prm%totalNtrans > 0) &
 dst%invLambdaTrans(1:prm%totalNtrans,of) = matmul(prm%interaction_TransTrans,ftransOverLamellarSize)/(1.0_pReal-sumf_trans)

 !* mean free path between 2 obstacles seen by a moving dislocation
 do i = 1,prm%totalNslip
    if ((prm%totalNtwin > 0) .or. (prm%totalNtrans > 0)) then              ! ToDo: Change order and use matmul
       dst%mfp_slip(i,of) = &
         prm%GrainSize/(1.0_pReal+prm%GrainSize*&
         (dst%invLambdaSlip(i,of) + dst%invLambdaSlipTwin(i,of) + dst%invLambdaSlipTrans(i,of)))
    else
       dst%mfp_slip(i,of) = prm%GrainSize &
                          / (1.0_pReal+prm%GrainSize*dst%invLambdaSlip(i,of)) !!!!!! correct?
    endif
 enddo

 !* mean free path between 2 obstacles seen by a growing twin/martensite
 dst%mfp_twin(:,of)  = prm%Cmfptwin*prm%GrainSize/ (1.0_pReal+prm%GrainSize*dst%invLambdaTwin(:,of))
 dst%mfp_trans(:,of) = prm%Cmfptrans*prm%GrainSize/(1.0_pReal+prm%GrainSize*dst%invLambdaTrans(:,of))

 !* threshold stress for dislocation motion
 forall (i = 1:prm%totalNslip) dst%tau_pass(i,of) = &
     prm%mu*prm%b_sl(i)*&
     sqrt(dot_product(stt%rhoEdge(1:prm%totalNslip,of)+stt%rhoEdgeDip(1:prm%totalNslip,of),&
                      prm%h_sl_sl(:,i)))

 !* threshold stress for growing twin/martensite
 if(prm%totalNtwin == prm%totalNslip) &
 dst%threshold_stress_twin(:,of) = &
             (SFE/(3.0_pReal*prm%b_tw)+ 3.0_pReal*prm%b_tw*prm%mu/(prm%L0_twin*prm%b_sl)) ! slip burgers here correct?
 if(prm%totalNtrans == prm%totalNslip) &
   dst%threshold_stress_trans(:,of) =  &
         (SFE/(3.0_pReal*prm%burgers_trans) + 3.0_pReal*prm%burgers_trans*prm%mu/&
              (prm%L0_trans*prm%b_sl) + prm%transStackHeight*prm%deltaG/ (3.0_pReal*prm%burgers_trans) )    
 

 dst%twinVolume(:,of)       = (PI/4.0_pReal)*prm%twinsize*dst%mfp_twin(:,of)**2.0_pReal
 dst%martensiteVolume(:,of) = (PI/4.0_pReal)*prm%lamellarsize*dst%mfp_trans(:,of)**2.0_pReal


 x0 = prm%mu*prm%b_tw**2.0_pReal/(SFE*8.0_pReal*PI)*(2.0_pReal+prm%nu)/(1.0_pReal-prm%nu)  ! ToDo: In the paper, this is the burgers vector for slip
 dst%tau_r_twin(:,of) = prm%mu*prm%b_tw/(2.0_pReal*PI)*(1.0_pReal/(x0+prm%xc_twin)+cos(pi/3.0_pReal)/x0)

 x0 = prm%mu*prm%burgers_trans**2.0_pReal/(SFE*8.0_pReal*PI)*(2.0_pReal+prm%nu)/(1.0_pReal-prm%nu) ! ToDo: In the paper, this is the burgers vector for slip
 dst%tau_r_trans(:,of) = prm%mu*prm%burgers_trans/(2.0_pReal*PI)*(1.0_pReal/(x0+prm%xc_trans)+cos(pi/3.0_pReal)/x0)

 end associate

end subroutine plastic_dislotwin_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_dislotwin_postResults(Mp,Temperature,instance,of) result(postResults)
 use prec, only: &
   tol_math_check, &
   dEq0
 use math, only: &
   PI, &
   math_mul33xx33

 implicit none
 real(pReal), dimension(3,3),intent(in) :: &
   Mp                                                                                               !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer,                    intent(in) :: &
   instance, &
   of

 real(pReal), dimension(sum(plastic_dislotwin_sizePostResult(:,instance))) :: &
   postResults

 integer :: &
   o,c,j

 associate(prm => param(instance), stt => state(instance), dst => microstructure(instance))
 
 c = 0

 do o = 1,size(prm%outputID)
   select case(prm%outputID(o))
 
     case (rho_mob_ID)
       postResults(c+1:c+prm%totalNslip) = stt%rhoEdge(1:prm%totalNslip,of)
       c = c + prm%totalNslip
     case (rho_dip_ID)
       postResults(c+1:c+prm%totalNslip) = stt%rhoEdgeDip(1:prm%totalNslip,of)
       c = c + prm%totalNslip
     case (gamma_dot_sl_ID)
       call kinetics_slip(Mp,temperature,instance,of,postResults(c+1:c+prm%totalNslip))
       c = c + prm%totalNslip
     case (gamma_sl_ID)
      postResults(c+1:c+prm%totalNslip)  = stt%accshear_slip(1:prm%totalNslip,of)
       c = c + prm%totalNslip
     case (mfp_slip_ID)
       postResults(c+1:c+prm%totalNslip) = dst%mfp_slip(1:prm%totalNslip,of)
       c = c + prm%totalNslip
     case (resolved_stress_slip_ID)
       do j = 1, prm%totalNslip
         postResults(c+j) = math_mul33xx33(Mp,prm%Schmid_slip(1:3,1:3,j))
       enddo
       c = c + prm%totalNslip
     case (threshold_stress_slip_ID)
       postResults(c+1:c+prm%totalNslip) = dst%tau_pass(1:prm%totalNslip,of)
       c = c + prm%totalNslip

     case (f_tw_ID)
       postResults(c+1:c+prm%totalNtwin) = stt%twinFraction(1:prm%totalNtwin,of)
       c = c + prm%totalNtwin    
     case (mfp_twin_ID)
       postResults(c+1:c+prm%totalNtwin) = dst%mfp_twin(1:prm%totalNtwin,of)
       c = c + prm%totalNtwin
     case (resolved_stress_twin_ID)
       do j = 1, prm%totalNtwin
         postResults(c+j) = math_mul33xx33(Mp,prm%Schmid_twin(1:3,1:3,j))
       enddo
       c = c + prm%totalNtwin
     case (threshold_stress_twin_ID)
       postResults(c+1:c+prm%totalNtwin) = dst%threshold_stress_twin(1:prm%totalNtwin,of)
       c = c + prm%totalNtwin

     case (strain_trans_fraction_ID)
       postResults(c+1:c+prm%totalNtrans) = stt%strainTransFraction(1:prm%totalNtrans,of)
       c = c + prm%totalNtrans
   end select
 enddo
 
 end associate
 
end function plastic_dislotwin_postResults


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_results(instance,group)
#if defined(PETSc) || defined(DAMASKHDF5)
  use results

  implicit none
  integer, intent(in) :: instance
  character(len=*) :: group
  integer :: o

  associate(prm => param(instance), stt => state(instance))
  outputsLoop: do o = 1,size(prm%outputID)
    select case(prm%outputID(o))
    end select
  enddo outputsLoop
  end associate
#else
  integer, intent(in) :: instance
  character(len=*) :: group
#endif

end subroutine plastic_dislotwin_results


!--------------------------------------------------------------------------------------------------
!> @brief Shear rates on slip systems, their derivatives with respect to resolved stress and the
!  resolved stresss
!> @details Derivatives and resolved stress are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_slip(Mp,Temperature,instance,of, &
                              gdot_slip,dgdot_dtau_slip,tau_slip)
 use prec, only: &
  tol_math_check, &
  dNeq0
 use math, only: &
   math_mul33xx33

 implicit none
 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                  intent(in) :: &
   temperature                                                                                      !< temperature
 integer,                      intent(in) :: &
   instance, &
   of
   
 real(pReal), dimension(param(instance)%totalNslip), intent(out) :: &
   gdot_slip
 real(pReal), dimension(param(instance)%totalNslip), optional, intent(out) :: &
   dgdot_dtau_slip, &
   tau_slip
 real(pReal), dimension(param(instance)%totalNslip) :: &
   dgdot_dtau

 real(pReal), dimension(param(instance)%totalNslip) :: &
   tau, &
   stressRatio, &
   StressRatio_p, &
   BoltzmannRatio, &
   v_wait_inverse, &                                                                                !< inverse of the effective velocity of a dislocation waiting at obstacles (unsigned)
   v_run_inverse, &                                                                                 !< inverse of the velocity of a free moving dislocation (unsigned)
   dV_wait_inverse_dTau, &
   dV_run_inverse_dTau, &
   dV_dTau, &
   tau_eff                                                                                          !< effective resolved stress
 integer :: i 
 
 associate(prm => param(instance), stt => state(instance), dst => microstructure(instance))

 do i = 1, prm%totalNslip
   tau(i) = math_mul33xx33(Mp,prm%Schmid_slip(1:3,1:3,i))
 enddo
 
 tau_eff = abs(tau)-dst%tau_pass(:,of)
   
 significantStress: where(tau_eff > tol_math_check)
   stressRatio    = tau_eff/(prm%SolidSolutionStrength+prm%tau_peierls)
   StressRatio_p  = stressRatio** prm%p
   BoltzmannRatio = prm%Qedge/(kB*Temperature)
   v_wait_inverse = prm%v0**(-1.0_pReal) * exp(BoltzmannRatio*(1.0_pReal-StressRatio_p)** prm%q)
   v_run_inverse  = prm%B/(tau_eff*prm%b_sl)

   gdot_slip = sign(stt%rhoEdge(:,of)*prm%b_sl/(v_wait_inverse+v_run_inverse),tau)

   dV_wait_inverse_dTau = v_wait_inverse * prm%p * prm%q * BoltzmannRatio &
                        * (stressRatio**(prm%p-1.0_pReal)) &
                        * (1.0_pReal-StressRatio_p)**(prm%q-1.0_pReal) &
                        / (prm%SolidSolutionStrength+prm%tau_peierls)
   dV_run_inverse_dTau  = v_run_inverse/tau_eff
   dV_dTau              = (dV_wait_inverse_dTau+dV_run_inverse_dTau) &
                        / (v_wait_inverse+v_run_inverse)**2.0_pReal
   dgdot_dtau = dV_dTau*stt%rhoEdge(:,of)*prm%b_sl
 else where significantStress
   gdot_slip  = 0.0_pReal
   dgdot_dtau = 0.0_pReal
 end where significantStress
 
 end associate
 
 if(present(dgdot_dtau_slip)) dgdot_dtau_slip = dgdot_dtau
 if(present(tau_slip))        tau_slip        = tau
 
end subroutine kinetics_slip


!--------------------------------------------------------------------------------------------------
!> @brief calculates shear rates on twin systems
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_twin(Mp,temperature,gdot_slip,instance,of,&
                              gdot_twin,dgdot_dtau_twin)
 use prec, only: &
   tol_math_check, &
   dNeq0
 use math, only: &
   math_mul33xx33

 implicit none
 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                  intent(in) :: &
   temperature                                                                                      !< temperature
 integer,                      intent(in) :: &
   instance, &
   of
 real(pReal), dimension(param(instance)%totalNslip), intent(in) :: &
   gdot_slip
   
 real(pReal), dimension(param(instance)%totalNtwin), intent(out) :: &
   gdot_twin
 real(pReal), dimension(param(instance)%totalNtwin), optional, intent(out) :: &
   dgdot_dtau_twin

 real, dimension(param(instance)%totalNtwin) :: &
   tau, &
   Ndot0, &
   stressRatio_r, &
   dgdot_dtau

 integer :: i,s1,s2
 
 associate(prm => param(instance), stt => state(instance), dst => microstructure(instance))

 do i = 1, prm%totalNtwin
   tau(i) = math_mul33xx33(Mp,prm%Schmid_twin(1:3,1:3,i))
   isFCC: if (prm%fccTwinTransNucleation) then
     s1=prm%fcc_twinNucleationSlipPair(1,i)
     s2=prm%fcc_twinNucleationSlipPair(2,i)
     if (tau(i) < dst%tau_r_twin(i,of)) then
       Ndot0=(abs(gdot_slip(s1))*(stt%rhoEdge(s2,of)+stt%rhoEdgeDip(s2,of))+&
              abs(gdot_slip(s2))*(stt%rhoEdge(s1,of)+stt%rhoEdgeDip(s1,of)))/&                      ! ToDo: MD: it would be more consistent to use shearrates from state
               (prm%L0_twin*prm%b_sl(i))*&
               (1.0_pReal-exp(-prm%VcrossSlip/(kB*Temperature)*&
               (dst%tau_r_twin(i,of)-tau)))
     else
       Ndot0=0.0_pReal
     end if
   else isFCC
     Ndot0=prm%Ndot0_twin(i)
   endif isFCC
 enddo

 significantStress: where(tau > tol_math_check)
   StressRatio_r = (dst%threshold_stress_twin(:,of)/tau)**prm%r
   gdot_twin  = prm%shear_twin * dst%twinVolume(:,of) * Ndot0*exp(-StressRatio_r)
   dgdot_dtau = (gdot_twin*prm%r/tau)*StressRatio_r
 else where significantStress
   gdot_twin  = 0.0_pReal
   dgdot_dtau = 0.0_pReal
 end where significantStress
 
 end associate

 if(present(dgdot_dtau_twin)) dgdot_dtau_twin = dgdot_dtau

end subroutine kinetics_twin


!--------------------------------------------------------------------------------------------------
!> @brief calculates shear rates on twin systems
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_trans(Mp,temperature,gdot_slip,instance,of,&
                              gdot_trans,dgdot_dtau_trans)
 use prec, only: &
   tol_math_check, &
   dNeq0
 use math, only: &
   math_mul33xx33

 implicit none
 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                  intent(in) :: &
   temperature                                                                                      !< temperature
 integer,                      intent(in) :: &
   instance, &
   of
 real(pReal), dimension(param(instance)%totalNslip), intent(in) :: &
   gdot_slip
   
 real(pReal), dimension(param(instance)%totalNtrans), intent(out) :: &
   gdot_trans
 real(pReal), dimension(param(instance)%totalNtrans), optional, intent(out) :: &
   dgdot_dtau_trans

 real, dimension(param(instance)%totalNtrans) :: &
   tau, &
   Ndot0, &
   stressRatio_s, &
   dgdot_dtau

 integer :: i,s1,s2
 
 associate(prm => param(instance), stt => state(instance), dst => microstructure(instance))

 do i = 1, prm%totalNtrans
   tau(i) = math_mul33xx33(Mp,prm%Schmid_trans(1:3,1:3,i))
   isFCC: if (prm%fccTwinTransNucleation) then
     s1=prm%fcc_twinNucleationSlipPair(1,i)
     s2=prm%fcc_twinNucleationSlipPair(2,i)
     if (tau(i) < dst%tau_r_trans(i,of)) then
       Ndot0=(abs(gdot_slip(s1))*(stt%rhoEdge(s2,of)+stt%rhoEdgeDip(s2,of))+&
              abs(gdot_slip(s2))*(stt%rhoEdge(s1,of)+stt%rhoEdgeDip(s1,of)))/&                      ! ToDo: MD: it would be more consistent to use shearrates from state
               (prm%L0_trans*prm%b_sl(i))*&
               (1.0_pReal-exp(-prm%VcrossSlip/(kB*Temperature)*&
               (dst%tau_r_trans(i,of)-tau)))
     else
       Ndot0=0.0_pReal
     end if
   else isFCC
     Ndot0=prm%Ndot0_trans(i)
   endif isFCC
 enddo

 significantStress: where(tau > tol_math_check)
   StressRatio_s = (dst%threshold_stress_trans(:,of)/tau)**prm%s
   gdot_trans  = dst%martensiteVolume(:,of) * Ndot0*exp(-StressRatio_s)
   dgdot_dtau = (gdot_trans*prm%r/tau)*StressRatio_s
 else where significantStress
   gdot_trans  = 0.0_pReal
   dgdot_dtau = 0.0_pReal
 end where significantStress
 
 end associate

 if(present(dgdot_dtau_trans)) dgdot_dtau_trans = dgdot_dtau

end subroutine kinetics_trans

end module plastic_dislotwin
