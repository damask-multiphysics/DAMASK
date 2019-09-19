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
  use prec
  use debug
  use math
  use IO
  use material
  use config
  use lattice
  use discretization
  use results

 implicit none
 private
 
 integer,                       dimension(:,:),  allocatable, target, public :: &
   plastic_dislotwin_sizePostResult                                                                 !< size of each post result output
 character(len=64),             dimension(:,:),  allocatable, target, public :: &
   plastic_dislotwin_output                                                                         !< name of each post result output

 real(pReal),                                    parameter :: &
   kB = 1.38e-23_pReal                                                                              !< Boltzmann constant in J/Kelvin

 enum, bind(c)
   enumerator :: &
     undefined_ID, &
     rho_mob_ID, &
     rho_dip_ID, &
     dot_gamma_sl_ID, &
     gamma_sl_ID, &
     Lambda_sl_ID, &
     resolved_stress_slip_ID, &
     threshold_stress_slip_ID, &
     edge_dipole_distance_ID, &
     f_tw_ID, &
     Lambda_tw_ID, &
     resolved_stress_twin_ID, &
     tau_hat_tw_ID, &
     f_tr_ID
 end enum

 type :: tParameters
   real(pReal) :: &
     mu, &
     nu, &
     D0, &                                                                                          !< prefactor for self-diffusion coefficient
     Qsd, &                                                                                         !< activation energy for dislocation climb
     omega, &                                                                                       !< frequency factor for dislocation climb      
     D, &                                                                                           !< grain size
     p_sb, &                                                                                        !< p-exponent in shear band velocity
     q_sb, &                                                                                        !< q-exponent in shear band velocity
     CEdgeDipMinDistance, &                                                                         !<
     i_tw, &                                                                                        !<
     tau_0, &                                                                                       !< strength due to elements in solid solution
     L_tw, &                                                                                        !< Length of twin nuclei in Burgers vectors
     L_tr, &                                                                                        !< Length of trans nuclei in Burgers vectors
     xc_twin, &                                                                                     !< critical distance for formation of twin nucleus
     xc_trans, &                                                                                    !< critical distance for formation of trans nucleus
     V_cs, &                                                                                        !< cross slip volume
     sbResistance, &                                                                                !< value for shearband resistance (might become an internal state variable at some point)
     sbVelocity, &                                                                                  !< value for shearband velocity_0
     sbQedge, &                                                                                     !< activation energy for shear bands
     SFE_0K, &                                                                                      !< stacking fault energy at zero K
     dSFE_dT, &                                                                                     !< temperature dependance of stacking fault energy
     aTol_rho, &                                                                                    !< absolute tolerance for integration of dislocation density
     aTol_f_tw, &                                                                                   !< absolute tolerance for integration of twin volume fraction
     aTol_f_tr, &                                                                                   !< absolute tolerance for integration of trans volume fraction
     gamma_fcc_hex, &                                                                               !< Free energy difference between austensite and martensite
     i_tr, &                                                                                        !<
     h                                                                                              !< Stack height of hex nucleus 
   real(pReal),                  dimension(:),     allocatable :: & 
     rho_mob_0, &                                                                                   !< initial unipolar dislocation density per slip system
     rho_dip_0, &                                                                                   !< initial dipole dislocation density per slip system
     b_sl, &                                                                                        !< absolute length of burgers vector [m] for each slip system
     b_tw, &                                                                                        !< absolute length of burgers vector [m] for each twin system
     b_tr, &                                                                                        !< absolute length of burgers vector [m] for each transformation system
     Delta_F,&                                                                                      !< activation energy for glide [J] for each slip system
     v0, &                                                                                          !< dislocation velocity prefactor [m/s] for each slip system
     dot_N_0_tw, &                                                                                  !< twin nucleation rate [1/m³s] for each twin system
     dot_N_0_tr, &                                                                                  !< trans nucleation rate [1/m³s] for each trans system
     t_tw, &                                                                                        !< twin thickness [m] for each twin system
     CLambdaSlip, &                                                                                 !< Adj. parameter for distance between 2 forest dislocations for each slip system
     atomicVolume, &
     t_tr, &                                                                                        !< martensite lamellar thickness [m] for each trans system and instance
     p, &                                                                                           !< p-exponent in glide velocity
     q, &                                                                                           !< q-exponent in glide velocity
     r, &                                                                                           !< r-exponent in twin nucleation rate
     s, &                                                                                           !< s-exponent in trans nucleation rate
     gamma_char, &                                                                                  !< characteristic shear for twins
     B                                                                                              !< drag coefficient
   real(pReal),                  dimension(:,:),   allocatable :: & 
     h_sl_sl, &                                                                                     !< 
     h_sl_tw, &                                                                                     !<
     h_tw_tw, &                                                                                     !< 
     h_sl_tr, &                                                                                     !< 
     h_tr_tr                                                                                        !< 
   integer,                      dimension(:,:),   allocatable :: & 
     fcc_twinNucleationSlipPair                                                                     ! ToDo: Better name? Is also use for trans
   real(pReal),                  dimension(:,:),   allocatable :: & 
     n0_sl, &                                                                                       !< slip system normal
     forestProjection, &
     C66
   real(pReal),                  dimension(:,:,:), allocatable :: &
     P_tr, &
     P_sl, &
     P_tw, &
     C66_tw, &
     C66_tr
   integer :: & 
     sum_N_sl, &                                                                                    !< total number of active slip system
     sum_N_tw, &                                                                                    !< total number of active twin system
     sum_N_tr                                                                                       !< total number of active transformation system 
   integer,                      dimension(:),     allocatable :: & 
     N_sl, &                                                                                        !< number of active slip systems for each family
     N_tw, &                                                                                        !< number of active twin systems for each family
     N_tr                                                                                           !< number of active transformation systems for each family
   integer(kind(undefined_ID)),  dimension(:),     allocatable :: &
     outputID                                                                                       !< ID of each post result output
   logical :: &
     ExtendedDislocations, &                                                                        !< consider split into partials for climb calculation
     fccTwinTransNucleation, &                                                                      !< twinning and transformation models are for fcc
     dipoleFormation                                                                                !< flag indicating consideration of dipole formation
 end type                                                                                           !< container type for internal constitutive parameters

 type :: tDislotwinState
   real(pReal),                  dimension(:,:),   pointer :: &
     rho_mob, &
     rho_dip, &
     gamma_sl, &
     f_tw, &
     f_tr
 end type tDislotwinState

 type :: tDislotwinMicrostructure
   real(pReal),                  dimension(:,:),   allocatable :: &
     Lambda_sl, &                                                                                   !< mean free path between 2 obstacles seen by a moving dislocation
     Lambda_tw, &                                                                                   !< mean free path between 2 obstacles seen by a growing twin
     Lambda_tr, &                                                                                   !< mean free path between 2 obstacles seen by a growing martensite
     tau_pass, &
     tau_hat_tw, &
     tau_hat_tr, &
     V_tw, &                                                                                        !< volume of a new twin
     V_tr, &                                                                                        !< volume of a new martensite disc
     tau_r_tw, &                                                                                    !< stress to bring partials close together (twin)
     tau_r_tr                                                                                       !< stress to bring partials close together (trans)
 end type tDislotwinMicrostructure

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
 type(tParameters),              allocatable, dimension(:) :: param
 type(tDislotwinState),          allocatable, dimension(:) :: &
   dotState, &
   state
 type(tDislotwinMicrostructure), allocatable, dimension(:) :: dependentState

 public :: &
   plastic_dislotwin_init, &
   plastic_dislotwin_homogenizedC, &
   plastic_dislotwin_dependentState, &
   plastic_dislotwin_LpAndItsTangent, &
   plastic_dislotwin_dotState, &
   plastic_dislotwin_postResults, &
   plastic_dislotwin_results

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_init

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
 allocate(dependentState(Ninstance))

 do p = 1, size(phase_plasticity)
   if (phase_plasticity(p) /= PLASTICITY_DISLOTWIN_ID) cycle
   associate(prm => param(phase_plasticityInstance(p)), &
             dot => dotState(phase_plasticityInstance(p)), &
             stt => state(phase_plasticityInstance(p)), &
             dst => dependentState(phase_plasticityInstance(p)), &
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
   prm%sum_N_sl = sum(prm%N_sl)
   slipActive: if (prm%sum_N_sl > 0) then
     prm%P_sl          = lattice_SchmidMatrix_slip(prm%N_sl,config%getString('lattice_structure'),&
                                                          config%getFloat('c/a',defaultVal=0.0_pReal))
     prm%h_sl_sl = lattice_interaction_SlipBySlip(prm%N_sl, &
                                                  config%getFloats('interaction_slipslip'), &
                                                  config%getString('lattice_structure'))
     prm%forestProjection     = lattice_forestProjection (prm%N_sl,config%getString('lattice_structure'),&
                                                          config%getFloat('c/a',defaultVal=0.0_pReal))

     prm%n0_sl                = lattice_slip_normal(prm%N_sl,config%getString('lattice_structure'),&
                                                    config%getFloat('c/a',defaultVal=0.0_pReal))  
     prm%fccTwinTransNucleation = merge(.true., .false., lattice_structure(p) == LATTICE_FCC_ID) &
                                .and. (prm%N_sl(1) == 12)
     if(prm%fccTwinTransNucleation) &
       prm%fcc_twinNucleationSlipPair = lattice_fcc_twinNucleationSlipPair

     prm%rho_mob_0            = config%getFloats('rhoedge0',   requiredSize=size(prm%N_sl))
     prm%rho_dip_0            = config%getFloats('rhoedgedip0',requiredSize=size(prm%N_sl))
     prm%v0                   = config%getFloats('v0',         requiredSize=size(prm%N_sl))
     prm%b_sl                 = config%getFloats('slipburgers',requiredSize=size(prm%N_sl))
     prm%Delta_F              = config%getFloats('qedge',      requiredSize=size(prm%N_sl))
     prm%CLambdaSlip          = config%getFloats('clambdaslip',requiredSize=size(prm%N_sl))
     prm%p                    = config%getFloats('p_slip',     requiredSize=size(prm%N_sl))
     prm%q                    = config%getFloats('q_slip',     requiredSize=size(prm%N_sl))
     prm%B                    = config%getFloats('b',          requiredSize=size(prm%N_sl), &
                                                          defaultVal=[(0.0_pReal, i=1,size(prm%N_sl))])

     prm%tau_0                = config%getFloat('solidsolutionstrength')
     prm%CEdgeDipMinDistance  = config%getFloat('cedgedipmindistance')
     prm%D0                   = config%getFloat('d0')
     prm%Qsd                  = config%getFloat('qsd')
     prm%atomicVolume         = config%getFloat('catomicvolume') * prm%b_sl**3.0_pReal
     prm%ExtendedDislocations = config%keyExists('/extend_dislocations/')
     if (prm%ExtendedDislocations) then
       prm%SFE_0K               = config%getFloat('sfe_0k')
       prm%dSFE_dT              = config%getFloat('dsfe_dt')
     endif
  
     ! multiplication factor according to crystal structure (nearest neighbors bcc vs fcc/hex)
     !@details: Refer: Argon & Moffat, Acta Metallurgica, Vol. 29, pg 293 to 299, 1981
     prm%omega                = config%getFloat('omega',  defaultVal = 1000.0_pReal) &
                              * merge(12.0_pReal, &
                                      8.0_pReal, &
                                      lattice_structure(p) == LATTICE_FCC_ID .or. lattice_structure(p) == LATTICE_HEX_ID)


     ! expand: family => system
     prm%rho_mob_0    = math_expand(prm%rho_mob_0,   prm%N_sl)
     prm%rho_dip_0    = math_expand(prm%rho_dip_0,   prm%N_sl)
     prm%v0           = math_expand(prm%v0,          prm%N_sl)
     prm%b_sl         = math_expand(prm%b_sl,        prm%N_sl)
     prm%Delta_F      = math_expand(prm%Delta_F,     prm%N_sl)
     prm%CLambdaSlip  = math_expand(prm%CLambdaSlip, prm%N_sl)
     prm%p            = math_expand(prm%p,           prm%N_sl)
     prm%q            = math_expand(prm%q,           prm%N_sl)
     prm%B            = math_expand(prm%B,           prm%N_sl)
     prm%atomicVolume = math_expand(prm%atomicVolume,prm%N_sl)                                   

     ! sanity checks
     if (    prm%D0           <= 0.0_pReal)          extmsg = trim(extmsg)//' D0'
     if (    prm%Qsd          <= 0.0_pReal)          extmsg = trim(extmsg)//' Qsd'
     if (any(prm%rho_mob_0    <  0.0_pReal))         extmsg = trim(extmsg)//' rho_mob_0'
     if (any(prm%rho_dip_0    <  0.0_pReal))         extmsg = trim(extmsg)//' rho_dip_0'
     if (any(prm%v0           <  0.0_pReal))         extmsg = trim(extmsg)//' v0'
     if (any(prm%b_sl         <= 0.0_pReal))         extmsg = trim(extmsg)//' b_sl'
     if (any(prm%Delta_F      <= 0.0_pReal))         extmsg = trim(extmsg)//' Delta_F'
     if (any(prm%CLambdaSlip  <= 0.0_pReal))         extmsg = trim(extmsg)//' CLambdaSlip'
     if (any(prm%B            <  0.0_pReal))         extmsg = trim(extmsg)//' B'
     if (any(prm%p<=0.0_pReal .or. prm%p>1.0_pReal)) extmsg = trim(extmsg)//' p'
     if (any(prm%q< 1.0_pReal .or. prm%q>2.0_pReal)) extmsg = trim(extmsg)//' q'

   else slipActive
     allocate(prm%b_sl(0))
   endif slipActive

!--------------------------------------------------------------------------------------------------
! twin related parameters
   prm%N_tw      = config%getInts('ntwin', defaultVal=emptyIntArray)
   prm%sum_N_tw = sum(prm%N_tw)
   if (prm%sum_N_tw > 0) then
     prm%P_tw  = lattice_SchmidMatrix_twin(prm%N_tw,config%getString('lattice_structure'),&
                                                  config%getFloat('c/a',defaultVal=0.0_pReal))
     prm%h_tw_tw   = lattice_interaction_TwinByTwin(prm%N_tw,&
                                                    config%getFloats('interaction_twintwin'), &
                                                    config%getString('lattice_structure'))

     prm%b_tw      = config%getFloats('twinburgers',  requiredSize=size(prm%N_tw))
     prm%t_tw      = config%getFloats('twinsize',     requiredSize=size(prm%N_tw))
     prm%r         = config%getFloats('r_twin',       requiredSize=size(prm%N_tw))

     prm%xc_twin   = config%getFloat('xc_twin')
     prm%L_tw      = config%getFloat('l0_twin')
     prm%i_tw      = config%getFloat('cmfptwin')

     prm%gamma_char= lattice_characteristicShear_Twin(prm%N_tw,config%getString('lattice_structure'),&
                                                      config%getFloat('c/a',defaultVal=0.0_pReal))

     prm%C66_tw    = lattice_C66_twin(prm%N_tw,prm%C66,config%getString('lattice_structure'),&
                                      config%getFloat('c/a',defaultVal=0.0_pReal))

     if (.not. prm%fccTwinTransNucleation) then
       prm%dot_N_0_tw = config%getFloats('ndot0_twin') 
       prm%dot_N_0_tw = math_expand(prm%dot_N_0_tw,prm%N_tw)
     endif

     ! expand: family => system
     prm%b_tw         = math_expand(prm%b_tw,prm%N_tw)
     prm%t_tw         = math_expand(prm%t_tw,prm%N_tw)
     prm%r            = math_expand(prm%r,prm%N_tw)
     
   else
     allocate(prm%gamma_char(0))
     allocate(prm%t_tw      (0))
     allocate(prm%b_tw      (0))
     allocate(prm%r         (0))
     allocate(prm%h_tw_tw   (0,0))
   endif
  
!--------------------------------------------------------------------------------------------------
! transformation related parameters
   prm%N_tr     = config%getInts('ntrans', defaultVal=emptyIntArray)
   prm%sum_N_tr = sum(prm%N_tr)
   if (prm%sum_N_tr > 0) then
     prm%b_tr = config%getFloats('transburgers')
     prm%b_tr = math_expand(prm%b_tr,prm%N_tr)
     
     prm%h             = config%getFloat('transstackheight', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%i_tr          = config%getFloat('cmfptrans', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%gamma_fcc_hex = config%getFloat('deltag')
     prm%xc_trans      = config%getFloat('xc_trans', defaultVal=0.0_pReal) ! ToDo: How to handle that???
     prm%L_tr          = config%getFloat('l0_trans')

     prm%h_tr_tr       = lattice_interaction_TransByTrans(prm%N_tr,&
                                                          config%getFloats('interaction_transtrans'), &
                                                          config%getString('lattice_structure'))
                                                             
     prm%C66_tr        = lattice_C66_trans(prm%N_tr,prm%C66, &
                                  config%getString('trans_lattice_structure'), &
                                  0.0_pReal, &
                                  config%getFloat('a_bcc', defaultVal=0.0_pReal), &
                                  config%getFloat('a_fcc', defaultVal=0.0_pReal))
                                  
      prm%P_tr         = lattice_SchmidMatrix_trans(prm%N_tr, &
                                  config%getString('trans_lattice_structure'), &
                                  0.0_pReal, &
                                  config%getFloat('a_bcc', defaultVal=0.0_pReal), &
                                  config%getFloat('a_fcc', defaultVal=0.0_pReal))
                                                 
     if (lattice_structure(p) /= LATTICE_fcc_ID) then
        prm%dot_N_0_tr = config%getFloats('ndot0_trans')
        prm%dot_N_0_tr = math_expand(prm%dot_N_0_tr,prm%N_tr)
     endif
     prm%t_tr = config%getFloats('lamellarsize')
     prm%t_tr = math_expand(prm%t_tr,prm%N_tr)
     prm%s    = config%getFloats('s_trans',defaultVal=[0.0_pReal])
     prm%s    = math_expand(prm%s,prm%N_tr)
   else
     allocate(prm%t_tr   (0))
     allocate(prm%b_tr   (0))
     allocate(prm%s      (0))
     allocate(prm%h_tr_tr(0,0))
   endif
   
   if (sum(prm%N_tw) > 0  .or. prm%sum_N_tr > 0) then
     prm%SFE_0K     = config%getFloat('sfe_0k')
     prm%dSFE_dT    = config%getFloat('dsfe_dt')
     prm%V_cs       = config%getFloat('vcrossslip')
   endif
   
   if (prm%sum_N_sl > 0 .and. prm%sum_N_tw > 0) then
     prm%h_sl_tw = lattice_interaction_SlipByTwin(prm%N_sl,prm%N_tw,&
                                                  config%getFloats('interaction_sliptwin'), &
                                                  config%getString('lattice_structure'))
     if (prm%fccTwinTransNucleation .and. prm%sum_N_tw > 12) write(6,*) 'mist' ! ToDo: implement better test. The model will fail also if N_tw is [6,6]
   endif    

   if (prm%sum_N_sl > 0 .and. prm%sum_N_tr > 0) then  
     prm%h_sl_tr = lattice_interaction_SlipByTrans(prm%N_sl,prm%N_tr,&
                                                                 config%getFloats('interaction_sliptrans'), &
                                                                 config%getString('lattice_structure'))
     if (prm%fccTwinTransNucleation .and. prm%sum_N_tr > 12) write(6,*) 'mist' ! ToDo: implement better test. The model will fail also if N_tr is [6,6]
   endif  
  
!--------------------------------------------------------------------------------------------------
! shearband related parameters
   prm%sbVelocity = config%getFloat('shearbandvelocity',defaultVal=0.0_pReal)
   if (prm%sbVelocity > 0.0_pReal) then  
     prm%sbResistance = config%getFloat('shearbandresistance')
     prm%sbQedge      = config%getFloat('qedgepersbsystem')
     prm%p_sb         = config%getFloat('p_shearband')
     prm%q_sb         = config%getFloat('q_shearband')
     
     ! sanity checks
     if (prm%sbResistance  <  0.0_pReal) extmsg = trim(extmsg)//' shearbandresistance'
     if (prm%sbQedge       <  0.0_pReal) extmsg = trim(extmsg)//' qedgepersbsystem'
     if (prm%p_sb          <= 0.0_pReal) extmsg = trim(extmsg)//' p_shearband'
     if (prm%q_sb          <= 0.0_pReal) extmsg = trim(extmsg)//' q_shearband'
   endif



   prm%D             = config%getFloat('grainsize')

   if (config%keyExists('dipoleformationfactor')) call IO_error(1,ext_msg='use /nodipoleformation/')
   prm%dipoleformation = .not. config%keyExists('/nodipoleformation/')
   

       !if (Ndot0PerTwinFamily(f,p) < 0.0_pReal) &
        ! call IO_error(211,el=p,ext_msg='dot_N_0_tw ('//PLASTICITY_DISLOTWIN_label//')')

   if (any(prm%atomicVolume <= 0.0_pReal)) &
     call IO_error(211,el=p,ext_msg='cAtomicVolume ('//PLASTICITY_DISLOTWIN_label//')')
   if (prm%sum_N_tw > 0) then
     if (prm%aTol_rho <= 0.0_pReal) &
       call IO_error(211,el=p,ext_msg='aTol_rho ('//PLASTICITY_DISLOTWIN_label//')')   
     if (prm%aTol_f_tw <= 0.0_pReal) &
       call IO_error(211,el=p,ext_msg='aTol_f_tw ('//PLASTICITY_DISLOTWIN_label//')')
   endif
   if (prm%sum_N_tr > 0) then
     if (prm%aTol_f_tr <= 0.0_pReal) &
       call IO_error(211,el=p,ext_msg='aTol_f_tr ('//PLASTICITY_DISLOTWIN_label//')')
   endif
 
   outputs = config%getStrings('(output)', defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   do i= 1, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
       case ('rho_mob')
         outputID = merge(rho_mob_ID,undefined_ID,prm%sum_N_sl > 0)
         outputSize = prm%sum_N_sl
       case ('rho_dip')
         outputID = merge(rho_dip_ID,undefined_ID,prm%sum_N_sl > 0)
         outputSize = prm%sum_N_sl
       case ('gamma_sl')
         outputID = merge(gamma_sl_ID,undefined_ID,prm%sum_N_sl > 0)
         outputSize = prm%sum_N_sl
       case ('lambda_sl')
         outputID = merge(Lambda_sl_ID,undefined_ID,prm%sum_N_sl > 0)
         outputSize = prm%sum_N_sl
       case ('tau_pass')
         outputID= merge(threshold_stress_slip_ID,undefined_ID,prm%sum_N_sl > 0)
         outputSize = prm%sum_N_sl

       case ('f_tw')
         outputID = merge(f_tw_ID,undefined_ID,prm%sum_N_tw >0)
         outputSize = prm%sum_N_tw
       case ('lambda_tw')
         outputID = merge(Lambda_tw_ID,undefined_ID,prm%sum_N_tw >0)
         outputSize = prm%sum_N_tw
       case ('tau_hat_tw')
         outputID = merge(tau_hat_tw_ID,undefined_ID,prm%sum_N_tw >0)
         outputSize = prm%sum_N_tw
         
       case ('f_tr')
         outputID = f_tr_ID
         outputSize = prm%sum_N_tr
        
     end select
        
     if (outputID /= undefined_ID) then
       plastic_dislotwin_output(i,phase_plasticityInstance(p)) = outputs(i)
       plastic_dislotwin_sizePostResult(i,phase_plasticityInstance(p)) = outputSize
       prm%outputID = [prm%outputID, outputID]
     endif

   enddo

!--------------------------------------------------------------------------------------------------
! allocate state arrays
   NipcMyPhase  = count(material_phaseAt == p) * discretization_nIP
   sizeDotState = size(['rho_mob ','rho_dip ','gamma_sl']) * prm%sum_N_sl &
                + size(['f_tw'])                           * prm%sum_N_tw &
                + size(['f_tr'])                           * prm%sum_N_tr
   sizeState = sizeDotState

   call material_allocatePlasticState(p,NipcMyPhase,sizeState,sizeDotState,0, &
                                      prm%sum_N_sl,prm%sum_N_tw,prm%sum_N_tr)
   plasticState(p)%sizePostResults = sum(plastic_dislotwin_sizePostResult(:,phase_plasticityInstance(p)))


!--------------------------------------------------------------------------------------------------
! locally defined state aliases and initialization of state0 and aTolState
   startIndex = 1
   endIndex   = prm%sum_N_sl
   stt%rho_mob=>plasticState(p)%state(startIndex:endIndex,:)
   stt%rho_mob= spread(prm%rho_mob_0,2,NipcMyPhase)
   dot%rho_mob=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTol_rho

   startIndex = endIndex + 1
   endIndex   = endIndex + prm%sum_N_sl
   stt%rho_dip=>plasticState(p)%state(startIndex:endIndex,:)
   stt%rho_dip= spread(prm%rho_dip_0,2,NipcMyPhase)
   dot%rho_dip=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTol_rho

   startIndex = endIndex + 1
   endIndex   = endIndex + prm%sum_N_sl
   stt%gamma_sl=>plasticState(p)%state(startIndex:endIndex,:)
   dot%gamma_sl=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = 1.0e6_pReal  !ToDo: better make optional parameter
   ! global alias
   plasticState(p)%slipRate        => plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%accumulatedSlip => plasticState(p)%state(startIndex:endIndex,:)
   
   startIndex = endIndex + 1
   endIndex   = endIndex + prm%sum_N_tw
   stt%f_tw=>plasticState(p)%state(startIndex:endIndex,:)
   dot%f_tw=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTol_f_tw
   
   startIndex = endIndex + 1
   endIndex   = endIndex + prm%sum_N_tr
   stt%f_tr=>plasticState(p)%state(startIndex:endIndex,:)
   dot%f_tr=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTol_f_tr

   allocate(dst%Lambda_sl             (prm%sum_N_sl,NipcMyPhase),source=0.0_pReal)
   allocate(dst%tau_pass              (prm%sum_N_sl,NipcMyPhase),source=0.0_pReal)

   allocate(dst%Lambda_tw             (prm%sum_N_tw,NipcMyPhase),source=0.0_pReal)
   allocate(dst%tau_hat_tw            (prm%sum_N_tw,NipcMyPhase),source=0.0_pReal)
   allocate(dst%tau_r_tw              (prm%sum_N_tw,NipcMyPhase),source=0.0_pReal)
   allocate(dst%V_tw                  (prm%sum_N_tw,NipcMyPhase),source=0.0_pReal)

   allocate(dst%Lambda_tr             (prm%sum_N_tr,NipcMyPhase),source=0.0_pReal)
   allocate(dst%tau_hat_tr            (prm%sum_N_tr,NipcMyPhase),source=0.0_pReal)
   allocate(dst%tau_r_tr              (prm%sum_N_tr,NipcMyPhase),source=0.0_pReal)
   allocate(dst%V_tr                  (prm%sum_N_tr,NipcMyPhase),source=0.0_pReal)


   plasticState(p)%state0 = plasticState(p)%state                                                   ! ToDo: this could be done centrally

   end associate

 enddo

end subroutine plastic_dislotwin_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenized elasticity matrix
!--------------------------------------------------------------------------------------------------
function plastic_dislotwin_homogenizedC(ipc,ip,el) result(homogenizedC)
 
 real(pReal), dimension(6,6) :: &
   homogenizedC
 integer,     intent(in) :: &
   ipc, &                                                                                          !< component-ID of integration point
   ip, &                                                                                           !< integration point
   el                                                                                              !< element

 integer :: i, &
            of
 real(pReal) :: f_unrotated

 of = material_phasememberAt(ipc,ip,el)
 associate(prm => param(phase_plasticityInstance(material_phaseAt(ipc,el))),&
           stt => state(phase_plasticityInstance(material_phaseAT(ipc,el))))

 f_unrotated = 1.0_pReal &
             - sum(stt%f_tw(1:prm%sum_N_tw,of)) &
             - sum(stt%f_tr(1:prm%sum_N_tr,of))

 homogenizedC = f_unrotated * prm%C66
 do i=1,prm%sum_N_tw
   homogenizedC = homogenizedC &
                + stt%f_tw(i,of)*prm%C66_tw(1:6,1:6,i)
 enddo
 do i=1,prm%sum_N_tr
   homogenizedC = homogenizedC &
                + stt%f_tr(i,of)*prm%C66_tr(1:6,1:6,i)
 enddo

 end associate
 
end function plastic_dislotwin_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_LpAndItsTangent(Lp,dLp_dMp,Mp,T,instance,of)
 
 real(pReal), dimension(3,3),     intent(out) :: Lp
 real(pReal), dimension(3,3,3,3), intent(out) :: dLp_dMp
 real(pReal), dimension(3,3),     intent(in)  :: Mp
 integer,                         intent(in)  :: instance,of
 real(pReal),                     intent(in)  :: T

 integer :: i,k,l,m,n
 real(pReal) :: f_unrotated,StressRatio_p,&
                BoltzmannRatio, &
    ddot_gamma_dtau, &
    tau
 real(pReal), dimension(param(instance)%sum_N_sl) :: &
    dot_gamma_sl,ddot_gamma_dtau_slip
 real(pReal), dimension(param(instance)%sum_N_tw) :: &
    dot_gamma_twin,ddot_gamma_dtau_twin
 real(pReal), dimension(param(instance)%sum_N_tr) :: &
    dot_gamma_tr,ddot_gamma_dtau_trans
 real(pReal):: dot_gamma_sb
 real(pReal), dimension(3,3) :: eigVectors, P_sb
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
             - sum(stt%f_tw(1:prm%sum_N_tw,of)) &
             - sum(stt%f_tr(1:prm%sum_N_tr,of))

 Lp = 0.0_pReal
 dLp_dMp = 0.0_pReal 

 call kinetics_slip(Mp,T,instance,of,dot_gamma_sl,ddot_gamma_dtau_slip)
 slipContribution: do i = 1, prm%sum_N_sl
   Lp = Lp + dot_gamma_sl(i)*prm%P_sl(1:3,1:3,i)
   forall (k=1:3,l=1:3,m=1:3,n=1:3) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + ddot_gamma_dtau_slip(i) * prm%P_sl(k,l,i) * prm%P_sl(m,n,i)
 enddo slipContribution
 
 !ToDo: Why do this before shear banding?
 Lp      = Lp      * f_unrotated
 dLp_dMp = dLp_dMp * f_unrotated
 
 shearBandingContribution: if(dNeq0(prm%sbVelocity)) then

   BoltzmannRatio = prm%sbQedge/(kB*T)
   call math_eigenValuesVectorsSym(Mp,eigValues,eigVectors,error)

   do i = 1,6
     P_sb = 0.5_pReal * math_outer(matmul(eigVectors,sb_sComposition(1:3,i)),&
                                   matmul(eigVectors,sb_mComposition(1:3,i)))
     tau = math_mul33xx33(Mp,P_sb)
   
     significantShearBandStress: if (abs(tau) > tol_math_check) then
       StressRatio_p = (abs(tau)/prm%sbResistance)**prm%p_sb
       dot_gamma_sb = sign(prm%sbVelocity*exp(-BoltzmannRatio*(1-StressRatio_p)**prm%q_sb), tau)
       ddot_gamma_dtau = abs(dot_gamma_sb)*BoltzmannRatio* prm%p_sb*prm%q_sb/ prm%sbResistance &
                  * (abs(tau)/prm%sbResistance)**(prm%p_sb-1.0_pReal) &
                  * (1.0_pReal-StressRatio_p)**(prm%q_sb-1.0_pReal)
 
       Lp = Lp + dot_gamma_sb * P_sb
       forall (k=1:3,l=1:3,m=1:3,n=1:3) &
         dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                          + ddot_gamma_dtau * P_sb(k,l) * P_sb(m,n)
     endif significantShearBandStress
   enddo

 endif shearBandingContribution
 
 call kinetics_twin(Mp,T,dot_gamma_sl,instance,of,dot_gamma_twin,ddot_gamma_dtau_twin)
 twinContibution: do i = 1, prm%sum_N_tw
   Lp = Lp + dot_gamma_twin(i)*prm%P_tw(1:3,1:3,i) * f_unrotated
   forall (k=1:3,l=1:3,m=1:3,n=1:3) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + ddot_gamma_dtau_twin(i)* prm%P_tw(k,l,i)*prm%P_tw(m,n,i) * f_unrotated
 enddo twinContibution
 
 call kinetics_trans(Mp,T,dot_gamma_sl,instance,of,dot_gamma_tr,ddot_gamma_dtau_trans)
 transContibution: do i = 1, prm%sum_N_tr
   Lp = Lp + dot_gamma_tr(i)*prm%P_tr(1:3,1:3,i) * f_unrotated
   forall (k=1:3,l=1:3,m=1:3,n=1:3) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + ddot_gamma_dtau_trans(i)* prm%P_tr(k,l,i)*prm%P_tr(m,n,i) * f_unrotated
 enddo transContibution


 end associate
 
end subroutine plastic_dislotwin_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_dotState(Mp,T,instance,of)

 real(pReal), dimension(3,3),  intent(in):: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                  intent(in) :: &
   T                                                                                                !< temperature at integration point
 integer,                      intent(in) :: &
   instance, &
   of

 integer :: i
 real(pReal) :: &
   f_unrotated, &
   VacancyDiffusion, &
   rho_dip_distance, &
   v_cl, &                                                                                          !< climb velocity 
   Gamma, &                                                                                         !< stacking fault energy
   tau, &
   sigma_cl, &                                                                                      !< climb stress 
   b_d                                                                                              !< ratio of burgers vector to stacking fault width
 real(pReal), dimension(param(instance)%sum_N_sl) :: &
   dot_rho_dip_formation, &
   dot_rho_dip_climb, &
   rho_dip_distance_min, &
   dot_gamma_sl
 real(pReal), dimension(param(instance)%sum_N_tw) :: &
   dot_gamma_twin
 real(pReal), dimension(param(instance)%sum_N_tr) :: &
   dot_gamma_tr

 associate(prm => param(instance),    stt => state(instance), &
           dot => dotState(instance), dst => dependentState(instance))

 f_unrotated = 1.0_pReal &
             - sum(stt%f_tw(1:prm%sum_N_tw,of)) &
             - sum(stt%f_tr(1:prm%sum_N_tr,of))
 VacancyDiffusion = prm%D0*exp(-prm%Qsd/(kB*T))

 call kinetics_slip(Mp,T,instance,of,dot_gamma_sl)
 dot%gamma_sl(:,of) = abs(dot_gamma_sl)
 
 rho_dip_distance_min   = prm%CEdgeDipMinDistance*prm%b_sl
 
 slipState: do i = 1, prm%sum_N_sl
   tau = math_mul33xx33(Mp,prm%P_sl(1:3,1:3,i))

   significantSlipStress: if (dEq0(tau)) then
     dot_rho_dip_formation(i) = 0.0_pReal
     dot_rho_dip_climb(i) = 0.0_pReal
   else significantSlipStress
     rho_dip_distance = 3.0_pReal*prm%mu*prm%b_sl(i)/(16.0_pReal*PI*abs(tau))
     rho_dip_distance = math_clip(rho_dip_distance, right = dst%Lambda_sl(i,of))
     rho_dip_distance = math_clip(rho_dip_distance, left  = rho_dip_distance_min(i))

     if (prm%dipoleFormation) then
       dot_rho_dip_formation(i) = 2.0_pReal*(rho_dip_distance-rho_dip_distance_min(i))/prm%b_sl(i) &
                                * stt%rho_mob(i,of)*abs(dot_gamma_sl(i))
     else
       dot_rho_dip_formation(i) = 0.0_pReal
     endif

     if (dEq0(rho_dip_distance-rho_dip_distance_min(i))) then
       dot_rho_dip_climb(i) = 0.0_pReal
     else
     !@details: Refer: Argon & Moffat, Acta Metallurgica, Vol. 29, pg 293 to 299, 1981
       sigma_cl = dot_product(prm%n0_sl(1:3,i),matmul(Mp,prm%n0_sl(1:3,i)))
       if (prm%ExtendedDislocations) then            
         Gamma = prm%SFE_0K + prm%dSFE_dT * T
         b_d = 24.0_pReal*PI*(1.0_pReal - prm%nu)/(2.0_pReal + prm%nu)* Gamma/(prm%mu*prm%b_sl(i))
       else
         b_d = 1.0_pReal      
       endif
       v_cl = 2.0_pReal*prm%omega*b_d**2.0_pReal*exp(-prm%Qsd/(kB*T)) &
            * (exp(abs(sigma_cl)*prm%b_sl(i)**3.0_pReal/(kB*T)) - 1.0_pReal)
              
       dot_rho_dip_climb(i) = 4.0_pReal*v_cl*stt%rho_dip(i,of) &
                            / (rho_dip_distance-rho_dip_distance_min(i))
     endif
   endif significantSlipStress
 enddo slipState

 dot%rho_mob(:,of) = abs(dot_gamma_sl)/(prm%b_sl*dst%Lambda_sl(:,of)) &
                   - dot_rho_dip_formation &
                   - 2.0_pReal*rho_dip_distance_min/prm%b_sl * stt%rho_mob(:,of)*abs(dot_gamma_sl)

 dot%rho_dip(:,of) = dot_rho_dip_formation &
                   - 2.0_pReal*rho_dip_distance_min/prm%b_sl * stt%rho_dip(:,of)*abs(dot_gamma_sl) &
                   - dot_rho_dip_climb

 
 call kinetics_twin(Mp,T,dot_gamma_sl,instance,of,dot_gamma_twin)
 dot%f_tw(:,of) = f_unrotated*dot_gamma_twin/prm%gamma_char

 call kinetics_trans(Mp,T,dot_gamma_sl,instance,of,dot_gamma_tr)
 dot%f_tr(:,of) = f_unrotated*dot_gamma_tr

 end associate
 
end subroutine plastic_dislotwin_dotState


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_dependentState(T,instance,of)

 integer,       intent(in) :: &
   instance, &
   of
 real(pReal),   intent(in) :: &
   T

 integer :: &
   i
 real(pReal) :: &
   sumf_twin,SFE,sumf_trans
 real(pReal), dimension(param(instance)%sum_N_sl) :: &
   inv_lambda_sl_sl, &                                                                              !< 1/mean free distance between 2 forest dislocations seen by a moving dislocation
   inv_lambda_sl_tw, &                                                                              !< 1/mean free distance between 2 twin stacks from different systems seen by a moving dislocation
   inv_lambda_sl_tr                                                                                 !< 1/mean free distance between 2 martensite lamellar from different systems seen by a moving dislocation
 real(pReal), dimension(param(instance)%sum_N_tw) :: &
   inv_lambda_tw_tw, &                                                                              !< 1/mean free distance between 2 twin stacks from different systems seen by a growing twin
   f_over_t_tw
  real(pReal), dimension(param(instance)%sum_N_tr) :: &
   inv_lambda_tr_tr, &                                                                              !< 1/mean free distance between 2 martensite stacks from different systems seen by a growing martensite
   f_over_t_tr
 real(pReal), dimension(:), allocatable :: &
   x0


 associate(prm => param(instance),&
           stt => state(instance),&
           dst => dependentState(instance))

 sumf_twin  = sum(stt%f_tw(1:prm%sum_N_tw,of))
 sumf_trans = sum(stt%f_tr(1:prm%sum_N_tr,of))

 SFE = prm%SFE_0K + prm%dSFE_dT * T
 
 !* rescaled volume fraction for topology
 f_over_t_tw = stt%f_tw(1:prm%sum_N_tw,of)/prm%t_tw  ! this is per system ...
 f_over_t_tr = sumf_trans/prm%t_tr                   ! but this not
                                                     ! ToDo ...Physically correct, but naming could be adjusted


 forall (i = 1:prm%sum_N_sl) &
   inv_lambda_sl_sl(i) = &
     sqrt(dot_product((stt%rho_mob(1:prm%sum_N_sl,of)+stt%rho_dip(1:prm%sum_N_sl,of)),&
                      prm%forestProjection(1:prm%sum_N_sl,i)))/prm%CLambdaSlip(i) ! change order and use matmul

 
 if (prm%sum_N_tw > 0 .and. prm%sum_N_sl > 0) &
   inv_lambda_sl_tw = matmul(prm%h_sl_tw,f_over_t_tw)/(1.0_pReal-sumf_twin)

 inv_lambda_tw_tw = matmul(prm%h_tw_tw,f_over_t_tw)/(1.0_pReal-sumf_twin)
 
 if (prm%sum_N_tr > 0 .and. prm%sum_N_sl > 0) &
   inv_lambda_sl_tr = matmul(prm%h_sl_tr,f_over_t_tr)/(1.0_pReal-sumf_trans)
 
 inv_lambda_tr_tr = matmul(prm%h_tr_tr,f_over_t_tr)/(1.0_pReal-sumf_trans)

 

 if ((prm%sum_N_tw > 0) .or. (prm%sum_N_tr > 0)) then              ! ToDo: better logic needed here
   dst%Lambda_sl(:,of) = prm%D &
                       / (1.0_pReal+prm%D*(inv_lambda_sl_sl + inv_lambda_sl_tw + inv_lambda_sl_tr))
 else
   dst%Lambda_sl(:,of) = prm%D &
                       / (1.0_pReal+prm%D*inv_lambda_sl_sl) !!!!!! correct?
 endif


 dst%Lambda_tw(:,of) = prm%i_tw*prm%D/(1.0_pReal+prm%D*inv_lambda_tw_tw)
 dst%Lambda_tr(:,of) = prm%i_tr*prm%D/(1.0_pReal+prm%D*inv_lambda_tr_tr)

 !* threshold stress for dislocation motion
 dst%tau_pass(:,of) = prm%mu*prm%b_sl* sqrt(matmul(prm%h_sl_sl,stt%rho_mob(:,of)+stt%rho_dip(:,of)))

 !* threshold stress for growing twin/martensite
 if(prm%sum_N_tw == prm%sum_N_sl) &
   dst%tau_hat_tw(:,of) = SFE/(3.0_pReal*prm%b_tw) &
                        + 3.0_pReal*prm%b_tw*prm%mu/(prm%L_tw*prm%b_sl) ! slip burgers here correct?
 if(prm%sum_N_tr == prm%sum_N_sl) &
   dst%tau_hat_tr(:,of) = SFE/(3.0_pReal*prm%b_tr) &
                        + 3.0_pReal*prm%b_tr*prm%mu/(prm%L_tr*prm%b_sl) & ! slip burgers here correct?
                        + prm%h*prm%gamma_fcc_hex/ (3.0_pReal*prm%b_tr)  
 

 dst%V_tw(:,of) = (PI/4.0_pReal)*prm%t_tw*dst%Lambda_tw(:,of)**2.0_pReal
 dst%V_tr(:,of) = (PI/4.0_pReal)*prm%t_tr*dst%Lambda_tr(:,of)**2.0_pReal


 x0 = prm%mu*prm%b_tw**2.0_pReal/(SFE*8.0_pReal*PI)*(2.0_pReal+prm%nu)/(1.0_pReal-prm%nu)  ! ToDo: In the paper, this is the burgers vector for slip and is the same for twin and trans
 dst%tau_r_tw(:,of) = prm%mu*prm%b_tw/(2.0_pReal*PI)*(1.0_pReal/(x0+prm%xc_twin)+cos(pi/3.0_pReal)/x0)

 x0 = prm%mu*prm%b_tr**2.0_pReal/(SFE*8.0_pReal*PI)*(2.0_pReal+prm%nu)/(1.0_pReal-prm%nu) ! ToDo: In the paper, this is the burgers vector for slip
 dst%tau_r_tr(:,of) = prm%mu*prm%b_tr/(2.0_pReal*PI)*(1.0_pReal/(x0+prm%xc_trans)+cos(pi/3.0_pReal)/x0)

 end associate

end subroutine plastic_dislotwin_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_dislotwin_postResults(Mp,T,instance,of) result(postResults)

 real(pReal), dimension(3,3),intent(in) :: &
   Mp                                                                                               !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   T                                                                                                !< temperature at integration point
 integer,                    intent(in) :: &
   instance, &
   of

 real(pReal), dimension(sum(plastic_dislotwin_sizePostResult(:,instance))) :: &
   postResults

 integer :: &
   o,c,j

 associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))
 
 c = 0

 do o = 1,size(prm%outputID)
   select case(prm%outputID(o))
 
     case (rho_mob_ID)
       postResults(c+1:c+prm%sum_N_sl) = stt%rho_mob(1:prm%sum_N_sl,of)
       c = c + prm%sum_N_sl
     case (rho_dip_ID)
       postResults(c+1:c+prm%sum_N_sl) = stt%rho_dip(1:prm%sum_N_sl,of)
       c = c + prm%sum_N_sl
     case (dot_gamma_sl_ID)
       call kinetics_slip(Mp,T,instance,of,postResults(c+1:c+prm%sum_N_sl))
       c = c + prm%sum_N_sl
     case (gamma_sl_ID)
      postResults(c+1:c+prm%sum_N_sl)  = stt%gamma_sl(1:prm%sum_N_sl,of)
       c = c + prm%sum_N_sl
     case (Lambda_sl_ID)
       postResults(c+1:c+prm%sum_N_sl) = dst%Lambda_sl(1:prm%sum_N_sl,of)
       c = c + prm%sum_N_sl
     case (resolved_stress_slip_ID)
       do j = 1, prm%sum_N_sl
         postResults(c+j) = math_mul33xx33(Mp,prm%P_sl(1:3,1:3,j))
       enddo
       c = c + prm%sum_N_sl
     case (threshold_stress_slip_ID)
       postResults(c+1:c+prm%sum_N_sl) = dst%tau_pass(1:prm%sum_N_sl,of)
       c = c + prm%sum_N_sl

     case (f_tw_ID)
       postResults(c+1:c+prm%sum_N_tw) = stt%f_tw(1:prm%sum_N_tw,of)
       c = c + prm%sum_N_tw    
     case (Lambda_tw_ID)
       postResults(c+1:c+prm%sum_N_tw) = dst%Lambda_tw(1:prm%sum_N_tw,of)
       c = c + prm%sum_N_tw
     case (resolved_stress_twin_ID)
       do j = 1, prm%sum_N_tw
         postResults(c+j) = math_mul33xx33(Mp,prm%P_tw(1:3,1:3,j))
       enddo
       c = c + prm%sum_N_tw
     case (tau_hat_tw_ID)
       postResults(c+1:c+prm%sum_N_tw) = dst%tau_hat_tw(1:prm%sum_N_tw,of)
       c = c + prm%sum_N_tw

     case (f_tr_ID)
       postResults(c+1:c+prm%sum_N_tr) = stt%f_tr(1:prm%sum_N_tr,of)
       c = c + prm%sum_N_tr
   end select
 enddo
 
 end associate
 
end function plastic_dislotwin_postResults


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_results(instance,group)
#if defined(PETSc) || defined(DAMASK_HDF5)

  integer, intent(in) :: instance
  character(len=*) :: group
  integer :: o

  associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))
  outputsLoop: do o = 1,size(prm%outputID)
    select case(prm%outputID(o))

      case (rho_mob_ID)
        call results_writeDataset(group,stt%rho_mob,'rho_mob',&
                                  'mobile dislocation density','1/m²')
      case (rho_dip_ID)
        call results_writeDataset(group,stt%rho_dip,'rho_dip',&
                                  'dislocation dipole density''1/m²')
      case (gamma_sl_ID)
        call results_writeDataset(group,stt%gamma_sl,'gamma_sl',&
                                  'plastic shear','1')
      case (Lambda_sl_ID)
        call results_writeDataset(group,dst%Lambda_sl,'Lambda_sl',&
                                  'mean free path for slip','m')
      case (threshold_stress_slip_ID)
        call results_writeDataset(group,dst%tau_pass,'tau_pass',&
                                  'passing stress for slip','Pa')

      case (f_tw_ID)
        call results_writeDataset(group,stt%f_tw,'f_tw',&
                                 'twinned volume fraction','m³/m³')
      case (Lambda_tw_ID)
        call results_writeDataset(group,dst%Lambda_tw,'Lambda_tw',&
                                  'mean free path for twinning','m')   
      case (tau_hat_tw_ID)
        call results_writeDataset(group,dst%tau_hat_tw,'tau_hat_tw',&
                                  'threshold stress for twinning','Pa')

      case (f_tr_ID)
        call results_writeDataset(group,stt%f_tr,'f_tr',&
                                 'martensite volume fraction','m³/m³')
                                 
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
pure subroutine kinetics_slip(Mp,T,instance,of, &
                              dot_gamma_sl,ddot_gamma_dtau_slip,tau_slip)

 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                  intent(in) :: &
   T                                                                                                !< temperature
 integer,                      intent(in) :: &
   instance, &
   of
   
 real(pReal), dimension(param(instance)%sum_N_sl), intent(out) :: &
   dot_gamma_sl
 real(pReal), dimension(param(instance)%sum_N_sl), optional, intent(out) :: &
   ddot_gamma_dtau_slip, &
   tau_slip
 real(pReal), dimension(param(instance)%sum_N_sl) :: &
   ddot_gamma_dtau

 real(pReal), dimension(param(instance)%sum_N_sl) :: &
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
 
 associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))

 do i = 1, prm%sum_N_sl
   tau(i) = math_mul33xx33(Mp,prm%P_sl(1:3,1:3,i))
 enddo
 
 tau_eff = abs(tau)-dst%tau_pass(:,of)
   
 significantStress: where(tau_eff > tol_math_check)
   stressRatio    = tau_eff/prm%tau_0
   StressRatio_p  = stressRatio** prm%p
   BoltzmannRatio = prm%Delta_F/(kB*T)
   v_wait_inverse = prm%v0**(-1.0_pReal) * exp(BoltzmannRatio*(1.0_pReal-StressRatio_p)** prm%q)
   v_run_inverse  = prm%B/(tau_eff*prm%b_sl)

   dot_gamma_sl = sign(stt%rho_mob(:,of)*prm%b_sl/(v_wait_inverse+v_run_inverse),tau)

   dV_wait_inverse_dTau = -1.0_pReal * v_wait_inverse * prm%p * prm%q * BoltzmannRatio &
                        * (stressRatio**(prm%p-1.0_pReal)) &
                        * (1.0_pReal-StressRatio_p)**(prm%q-1.0_pReal) &
                        / prm%tau_0
   dV_run_inverse_dTau  = -1.0_pReal * v_run_inverse/tau_eff
   dV_dTau              = -1.0_pReal * (dV_wait_inverse_dTau+dV_run_inverse_dTau) &
                        / (v_wait_inverse+v_run_inverse)**2.0_pReal
   ddot_gamma_dtau = dV_dTau*stt%rho_mob(:,of)*prm%b_sl
 else where significantStress
   dot_gamma_sl    = 0.0_pReal
   ddot_gamma_dtau = 0.0_pReal
 end where significantStress
 
 end associate
  
 if(present(ddot_gamma_dtau_slip)) ddot_gamma_dtau_slip = ddot_gamma_dtau
 if(present(tau_slip))             tau_slip             = tau
 
end subroutine kinetics_slip


!--------------------------------------------------------------------------------------------------
!> @brief calculates shear rates on twin systems
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_twin(Mp,T,dot_gamma_sl,instance,of,&
                              dot_gamma_twin,ddot_gamma_dtau_twin)

 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                  intent(in) :: &
   T                                                                                                !< temperature
 integer,                      intent(in) :: &
   instance, &
   of
 real(pReal), dimension(param(instance)%sum_N_sl), intent(in) :: &
   dot_gamma_sl
   
 real(pReal), dimension(param(instance)%sum_N_tw), intent(out) :: &
   dot_gamma_twin
 real(pReal), dimension(param(instance)%sum_N_tw), optional, intent(out) :: &
   ddot_gamma_dtau_twin

 real, dimension(param(instance)%sum_N_tw) :: &
   tau, &
   Ndot0, &
   stressRatio_r, &
   ddot_gamma_dtau

 integer :: i,s1,s2
 
 associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))

 do i = 1, prm%sum_N_tw
   tau(i) = math_mul33xx33(Mp,prm%P_tw(1:3,1:3,i))
   isFCC: if (prm%fccTwinTransNucleation) then
     s1=prm%fcc_twinNucleationSlipPair(1,i)
     s2=prm%fcc_twinNucleationSlipPair(2,i)
     if (tau(i) < dst%tau_r_tw(i,of)) then                                                          ! ToDo: correct?
       Ndot0=(abs(dot_gamma_sl(s1))*(stt%rho_mob(s2,of)+stt%rho_dip(s2,of))+&
              abs(dot_gamma_sl(s2))*(stt%rho_mob(s1,of)+stt%rho_dip(s1,of)))/&                      ! ToDo: MD: it would be more consistent to use shearrates from state
               (prm%L_tw*prm%b_sl(i))*&
               (1.0_pReal-exp(-prm%V_cs/(kB*T)*(dst%tau_r_tw(i,of)-tau(i))))                        ! P_ncs
     else
       Ndot0=0.0_pReal
     end if
   else isFCC
     Ndot0=prm%dot_N_0_tw(i)
   endif isFCC
 enddo

 significantStress: where(tau > tol_math_check)
   StressRatio_r   = (dst%tau_hat_tw(:,of)/tau)**prm%r
   dot_gamma_twin  = prm%gamma_char * dst%V_tw(:,of) * Ndot0*exp(-StressRatio_r)
   ddot_gamma_dtau = (dot_gamma_twin*prm%r/tau)*StressRatio_r
 else where significantStress
   dot_gamma_twin  = 0.0_pReal
   ddot_gamma_dtau = 0.0_pReal
 end where significantStress
 
 end associate

 if(present(ddot_gamma_dtau_twin)) ddot_gamma_dtau_twin = ddot_gamma_dtau

end subroutine kinetics_twin


!--------------------------------------------------------------------------------------------------
!> @brief calculates shear rates on twin systems
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_trans(Mp,T,dot_gamma_sl,instance,of,&
                              dot_gamma_tr,ddot_gamma_dtau_trans)

 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                  intent(in) :: &
   T                                                                                                !< temperature
 integer,                      intent(in) :: &
   instance, &
   of
 real(pReal), dimension(param(instance)%sum_N_sl), intent(in) :: &
   dot_gamma_sl
   
 real(pReal), dimension(param(instance)%sum_N_tr), intent(out) :: &
   dot_gamma_tr
 real(pReal), dimension(param(instance)%sum_N_tr), optional, intent(out) :: &
   ddot_gamma_dtau_trans

 real, dimension(param(instance)%sum_N_tr) :: &
   tau, &
   Ndot0, &
   stressRatio_s, &
   ddot_gamma_dtau

 integer :: i,s1,s2
 associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))

 do i = 1, prm%sum_N_tr
   tau(i) = math_mul33xx33(Mp,prm%P_tr(1:3,1:3,i))
   isFCC: if (prm%fccTwinTransNucleation) then
     s1=prm%fcc_twinNucleationSlipPair(1,i)
     s2=prm%fcc_twinNucleationSlipPair(2,i)
     if (tau(i) < dst%tau_r_tr(i,of)) then                                                          ! ToDo: correct?
       Ndot0=(abs(dot_gamma_sl(s1))*(stt%rho_mob(s2,of)+stt%rho_dip(s2,of))+&
              abs(dot_gamma_sl(s2))*(stt%rho_mob(s1,of)+stt%rho_dip(s1,of)))/&                      ! ToDo: MD: it would be more consistent to use shearrates from state
               (prm%L_tr*prm%b_sl(i))*&
               (1.0_pReal-exp(-prm%V_cs/(kB*T)*(dst%tau_r_tr(i,of)-tau(i))))                        ! P_ncs
     else
       Ndot0=0.0_pReal
     end if
   else isFCC
     Ndot0=prm%dot_N_0_tr(i)
   endif isFCC
 enddo

 significantStress: where(tau > tol_math_check)
   StressRatio_s   = (dst%tau_hat_tr(:,of)/tau)**prm%s
   dot_gamma_tr    = dst%V_tr(:,of) * Ndot0*exp(-StressRatio_s)
   ddot_gamma_dtau = (dot_gamma_tr*prm%s/tau)*StressRatio_s
 else where significantStress
   dot_gamma_tr  = 0.0_pReal
   ddot_gamma_dtau = 0.0_pReal
 end where significantStress
 
 end associate

 if(present(ddot_gamma_dtau_trans)) ddot_gamma_dtau_trans = ddot_gamma_dtau

end subroutine kinetics_trans

end module plastic_dislotwin
