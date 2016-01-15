!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Parses material config file, either solverJobName.materialConfig or material.config
!> @details reads the material configuration file, where solverJobName.materialConfig takes
!! precedence over material.config and parses the sections 'homogenization', 'crystallite',
!! 'phase', 'texture', and 'microstucture'
!--------------------------------------------------------------------------------------------------
module material
 use prec, only: &
   pReal, &
   pInt, &
   tState, &
   tPlasticState, &
   tSourceState, &
   tHomogMapping, &
   tPhaseMapping, &
   p_vec, &
   p_intvec

 implicit none
 private
 character(len=*),                         parameter,            public :: &
   ELASTICITY_hooke_label               = 'hooke', &
   PLASTICITY_none_label                = 'none', &
   PLASTICITY_isotropic_label           = 'isotropic', &
   PLASTICITY_j2_label                  = 'j2', &
   PLASTICITY_phenopowerlaw_label       = 'phenopowerlaw', &
   PLASTICITY_phenoplus_label           = 'phenoplus', &
   PLASTICITY_dislotwin_label           = 'dislotwin', &
   PLASTICITY_disloucla_label           = 'disloucla', &
   PLASTICITY_titanmod_label            = 'titanmod', &
   PLASTICITY_nonlocal_label            = 'nonlocal', &
   SOURCE_thermal_dissipation_label     = 'thermal_dissipation', &
   SOURCE_thermal_externalheat_label    = 'thermal_externalheat', &
   SOURCE_damage_isoBrittle_label       = 'damage_isobrittle', &
   SOURCE_damage_isoDuctile_label       = 'damage_isoductile', &
   SOURCE_damage_anisoBrittle_label     = 'damage_anisobrittle', &
   SOURCE_damage_anisoDuctile_label     = 'damage_anisoductile', &
   SOURCE_vacancy_phenoplasticity_label = 'vacancy_phenoplasticity', &
   SOURCE_vacancy_irradiation_label     = 'vacancy_irradiation', &
   SOURCE_vacancy_thermalfluc_label     = 'vacancy_thermalfluctuation', &
   KINEMATICS_thermal_expansion_label   = 'thermal_expansion', &
   KINEMATICS_cleavage_opening_label    = 'cleavage_opening', &
   KINEMATICS_slipplane_opening_label   = 'slipplane_opening', &
   KINEMATICS_vacancy_strain_label      = 'vacancy_strain', &
   KINEMATICS_hydrogen_strain_label     = 'hydrogen_strain', &
   STIFFNESS_DEGRADATION_damage_label   = 'damage', &
   STIFFNESS_DEGRADATION_porosity_label = 'porosity', &
   THERMAL_isothermal_label             = 'isothermal', &
   THERMAL_adiabatic_label              = 'adiabatic', &
   THERMAL_conduction_label             = 'conduction', &
   DAMAGE_none_label                    = 'none', &
   DAMAGE_local_label                   = 'local', &
   DAMAGE_nonlocal_label                = 'nonlocal', &
   VACANCYFLUX_isoconc_label            = 'isoconcentration', &
   VACANCYFLUX_isochempot_label         = 'isochemicalpotential', &
   VACANCYFLUX_cahnhilliard_label       = 'cahnhilliard', &
   POROSITY_none_label                  = 'none', &
   POROSITY_phasefield_label            = 'phasefield', &
   HYDROGENFLUX_isoconc_label           = 'isoconcentration', &
   HYDROGENFLUX_cahnhilliard_label      = 'cahnhilliard', &
   HOMOGENIZATION_none_label            = 'none', &
   HOMOGENIZATION_isostrain_label       = 'isostrain', &
   HOMOGENIZATION_rgc_label             = 'rgc'



 enum, bind(c)
   enumerator :: ELASTICITY_undefined_ID, &
                 ELASTICITY_hooke_ID
 end enum
 enum, bind(c)
   enumerator :: PLASTICITY_undefined_ID, &
                 PLASTICITY_none_ID, &
                 PLASTICITY_isotropic_ID, &
                 PLASTICITY_j2_ID, &
                 PLASTICITY_phenopowerlaw_ID, &
                 PLASTICITY_phenoplus_ID, &
                 PLASTICITY_dislotwin_ID, &
                 PLASTICITY_disloucla_ID, &
                 PLASTICITY_titanmod_ID, &
                 PLASTICITY_nonlocal_ID
 end enum

 enum, bind(c)
   enumerator :: SOURCE_undefined_ID, &
                 SOURCE_thermal_dissipation_ID, &
                 SOURCE_thermal_externalheat_ID, &
                 SOURCE_damage_isoBrittle_ID, &
                 SOURCE_damage_isoDuctile_ID, &
                 SOURCE_damage_anisoBrittle_ID, &
                 SOURCE_damage_anisoDuctile_ID, &
                 SOURCE_vacancy_phenoplasticity_ID, &
                 SOURCE_vacancy_irradiation_ID, &
                 SOURCE_vacancy_thermalfluc_ID
 end enum

 enum, bind(c)
   enumerator :: KINEMATICS_undefined_ID, &
                 KINEMATICS_cleavage_opening_ID, &
                 KINEMATICS_slipplane_opening_ID, &
                 KINEMATICS_thermal_expansion_ID, &
                 KINEMATICS_vacancy_strain_ID, &
                 KINEMATICS_hydrogen_strain_ID
 end enum

 enum, bind(c)
   enumerator :: STIFFNESS_DEGRADATION_undefined_ID, &
                 STIFFNESS_DEGRADATION_damage_ID, &
                 STIFFNESS_DEGRADATION_porosity_ID
 end enum

 enum, bind(c)
   enumerator :: THERMAL_isothermal_ID, &
                 THERMAL_adiabatic_ID, &
                 THERMAL_conduction_ID
 end enum

 enum, bind(c)
   enumerator :: DAMAGE_none_ID, &
                 DAMAGE_local_ID, &
                 DAMAGE_nonlocal_ID
 end enum

 enum, bind(c)
   enumerator :: VACANCYFLUX_isoconc_ID, &
                 VACANCYFLUX_isochempot_ID, &
                 VACANCYFLUX_cahnhilliard_ID
 end enum

 enum, bind(c)
   enumerator :: POROSITY_none_ID, &
                 POROSITY_phasefield_ID
 end enum
 enum, bind(c)
   enumerator :: HYDROGENFLUX_isoconc_ID, &
                 HYDROGENFLUX_cahnhilliard_ID
 end enum

 enum, bind(c)
   enumerator :: HOMOGENIZATION_undefined_ID, &
                 HOMOGENIZATION_none_ID, &
                 HOMOGENIZATION_isostrain_ID, &
                 HOMOGENIZATION_rgc_ID
 end enum

 character(len=*), parameter, public  :: &
   MATERIAL_configFile         = 'material.config', &                                               !< generic name for material configuration file
   MATERIAL_localFileExt       = 'materialConfig'                                                   !< extension of solver job name depending material configuration file

 character(len=*), parameter, public  :: &
   MATERIAL_partHomogenization = 'homogenization', &                                                !< keyword for homogenization part
   MATERIAL_partCrystallite    = 'crystallite', &                                                   !< keyword for crystallite part
   MATERIAL_partPhase          = 'phase'                                                            !< keyword for phase part

 integer(kind(ELASTICITY_undefined_ID)),     dimension(:),   allocatable, public, protected :: &
   phase_elasticity                                                                                 !< elasticity of each phase
 integer(kind(PLASTICITY_undefined_ID)),     dimension(:),   allocatable, public, protected :: &
   phase_plasticity                                                                                 !< plasticity of each phase
 integer(kind(THERMAL_isothermal_ID)),       dimension(:),   allocatable, public, protected :: &
   thermal_type                                                                                     !< thermal transport model
 integer(kind(DAMAGE_none_ID)),              dimension(:),   allocatable, public, protected :: &
   damage_type                                                                                      !< nonlocal damage model
 integer(kind(VACANCYFLUX_isoconc_ID)),      dimension(:),   allocatable, public, protected :: &
   vacancyflux_type                                                                                 !< vacancy transport model
 integer(kind(POROSITY_none_ID)),            dimension(:),   allocatable, public, protected :: &
   porosity_type                                                                                    !< porosity evolution model
 integer(kind(HYDROGENFLUX_isoconc_ID)),     dimension(:),   allocatable, public, protected :: &
   hydrogenflux_type                                                                                !< hydrogen transport model

 integer(kind(SOURCE_undefined_ID)),         dimension(:,:), allocatable, public, protected :: &
   phase_source, &                                                                                  !< active sources mechanisms of each phase
   phase_kinematics, &                                                                              !< active kinematic mechanisms of each phase
   phase_stiffnessDegradation                                                                       !< active stiffness degradation mechanisms of each phase

 integer(kind(HOMOGENIZATION_undefined_ID)), dimension(:),   allocatable, public, protected :: &
   homogenization_type                                                                              !< type of each homogenization

 character(len=64), dimension(:), allocatable, public, protected :: &
   phase_name, &                                                                                    !< name of each phase
   homogenization_name, &                                                                           !< name of each homogenization
   crystallite_name                                                                                 !< name of each crystallite setting

 integer(pInt), public, protected :: &
   homogenization_maxNgrains, &                                                                     !< max number of grains in any USED homogenization
   material_Nphase, &                                                                               !< number of phases
   material_Nhomogenization, &                                                                      !< number of homogenizations
   material_Nmicrostructure, &                                                                      !< number of microstructures
   material_Ncrystallite                                                                            !< number of crystallite settings

 integer(pInt), dimension(:), allocatable, public, protected :: &
   phase_Nsources, &                                                                                !< number of source mechanisms active in each phase
   phase_Nkinematics, &                                                                             !< number of kinematic mechanisms active in each phase
   phase_NstiffnessDegradations, &                                                                  !< number of stiffness degradation mechanisms active in each phase
   phase_Noutput, &                                                                                 !< number of '(output)' items per phase
   phase_elasticityInstance, &                                                                      !< instance of particular elasticity of each phase
   phase_plasticityInstance                                                                         !< instance of particular plasticity of each phase

 integer(pInt), dimension(:), allocatable, public, protected :: &
   crystallite_Noutput                                                                              !< number of '(output)' items per crystallite setting

 integer(pInt), dimension(:), allocatable, public, protected :: &
   homogenization_Ngrains, &                                                                        !< number of grains in each homogenization
   homogenization_Noutput, &                                                                        !< number of '(output)' items per homogenization
   homogenization_typeInstance, &                                                                   !< instance of particular type of each homogenization
   thermal_typeInstance, &                                                                          !< instance of particular type of each thermal transport
   damage_typeInstance, &                                                                           !< instance of particular type of each nonlocal damage
   vacancyflux_typeInstance, &                                                                      !< instance of particular type of each vacancy flux
   porosity_typeInstance, &                                                                         !< instance of particular type of each porosity model
   hydrogenflux_typeInstance, &                                                                     !< instance of particular type of each hydrogen flux
   microstructure_crystallite                                                                       !< crystallite setting ID of each microstructure

 real(pReal), dimension(:), allocatable, public, protected :: &
   thermal_initialT, &                                                                              !< initial temperature per each homogenization
   damage_initialPhi, &                                                                             !< initial damage per each homogenization
   vacancyflux_initialCv, &                                                                         !< initial vacancy concentration per each homogenization
   porosity_initialPhi, &                                                                           !< initial posority per each homogenization
   hydrogenflux_initialCh                                                                           !< initial hydrogen concentration per each homogenization

 integer(pInt), dimension(:,:,:), allocatable, public :: &
   material_phase                                                                                   !< phase (index) of each grain,IP,element
 integer(pInt), dimension(:,:), allocatable, public :: &
   material_homog                                                                                   !< homogenization (index) of each IP,element
 type(tPlasticState), allocatable, dimension(:), public :: &
   plasticState
 type(tSourceState),  allocatable, dimension(:), public :: &
   sourceState
 type(tState),        allocatable, dimension(:), public :: &
   homogState, &
   thermalState, &
   damageState, &
   vacancyfluxState, &
   porosityState, &
   hydrogenfluxState

 integer(pInt), dimension(:,:,:), allocatable, public, protected :: &
   material_texture                                                                                 !< texture (index) of each grain,IP,element

 real(pReal), dimension(:,:,:,:), allocatable, public, protected :: &
   material_EulerAngles                                                                             !< initial orientation of each grain,IP,element

 logical, dimension(:), allocatable, public, protected :: &
   microstructure_active, &
   microstructure_elemhomo, &                                                                       !< flag to indicate homogeneous microstructure distribution over element's IPs
   phase_localPlasticity                                                                            !< flags phases with local constitutive law


 character(len=*), parameter, private :: &
   MATERIAL_partMicrostructure = 'microstructure', &                                                !< keyword for microstructure part
   MATERIAL_partTexture        = 'texture'                                                          !< keyword for texture part

 character(len=64), dimension(:), allocatable, private :: &
   microstructure_name, &                                                                           !< name of each microstructure
   texture_name                                                                                     !< name of each texture

 character(len=256), dimension(:), allocatable, private :: &
   texture_ODFfile                                                                                  !< name of each ODF file

 integer(pInt), private :: &
   material_Ntexture, &                                                                             !< number of textures
   microstructure_maxNconstituents, &                                                               !< max number of constituents in any phase
   texture_maxNgauss, &                                                                             !< max number of Gauss components in any texture
   texture_maxNfiber                                                                                !< max number of Fiber components in any texture

 integer(pInt), dimension(:), allocatable, private :: &
   microstructure_Nconstituents, &                                                                  !< number of constituents in each microstructure
   texture_symmetry, &                                                                              !< number of symmetric orientations per texture
   texture_Ngauss, &                                                                                !< number of Gauss components per texture
   texture_Nfiber                                                                                   !< number of Fiber components per texture

 integer(pInt), dimension(:,:), allocatable, private :: &
   microstructure_phase, &                                                                          !< phase IDs of each microstructure
   microstructure_texture                                                                           !< texture IDs of each microstructure

 real(pReal), dimension(:,:), allocatable, private :: &
   microstructure_fraction                                                                          !< vol fraction of each constituent in microstructure

 real(pReal), dimension(:,:,:), allocatable, private :: &
   material_volume, &                                                                               !< volume of each grain,IP,element
   texture_Gauss, &                                                                                 !< data of each Gauss component
   texture_Fiber, &                                                                                 !< data of each Fiber component
   texture_transformation                                                                           !< transformation for each texture

 logical, dimension(:), allocatable, private :: &
   homogenization_active

 integer(pInt), dimension(:,:,:), allocatable, public :: phaseAt                                    !< phase ID of every material point (ipc,ip,el)
 integer(pInt), dimension(:,:,:), allocatable, public :: phasememberAt                              !< memberID of given phase at every material point (ipc,ip,el)
 integer(pInt), dimension(:,:,:), allocatable, public, target :: mappingCrystallite
 integer(pInt), dimension(:,:,:), allocatable, public, target :: mappingHomogenization              !< mapping from material points to offset in heterogenous state/field
 integer(pInt), dimension(:,:),   allocatable, public, target :: mappingHomogenizationConst         !< mapping from material points to offset in constant state/field

 type(tHomogMapping), allocatable, dimension(:), public :: &
   thermalMapping, &                                                                                !< mapping for thermal state/fields
   damageMapping, &                                                                                 !< mapping for damage state/fields
   vacancyfluxMapping, &                                                                            !< mapping for vacancy conc state/fields
   porosityMapping, &                                                                               !< mapping for porosity state/fields
   hydrogenfluxMapping                                                                              !< mapping for hydrogen conc state/fields

 type(p_vec),         allocatable, dimension(:), public :: &
   temperature, &                                                                                   !< temperature field
   damage, &                                                                                        !< damage field
   vacancyConc, &                                                                                   !< vacancy conc field
   porosity, &                                                                                      !< porosity field
   hydrogenConc, &                                                                                  !< hydrogen conc field
   temperatureRate, &                                                                               !< temperature change rate field
   vacancyConcRate, &                                                                               !< vacancy conc change field
   hydrogenConcRate                                                                                 !< hydrogen conc change field

 public :: &
   material_init, &
   ELASTICITY_hooke_ID ,&
   PLASTICITY_none_ID, &
   PLASTICITY_isotropic_ID, &
   PLASTICITY_J2_ID, &
   PLASTICITY_phenopowerlaw_ID, &
   PLASTICITY_phenoplus_ID, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_titanmod_ID, &
   PLASTICITY_nonlocal_ID, &
   SOURCE_thermal_dissipation_ID, &
   SOURCE_thermal_externalheat_ID, &
   SOURCE_damage_isoBrittle_ID, &
   SOURCE_damage_isoDuctile_ID, &
   SOURCE_damage_anisoBrittle_ID, &
   SOURCE_damage_anisoDuctile_ID, &
   SOURCE_vacancy_phenoplasticity_ID, &
   SOURCE_vacancy_irradiation_ID, &
   SOURCE_vacancy_thermalfluc_ID, &
   KINEMATICS_cleavage_opening_ID, &
   KINEMATICS_slipplane_opening_ID, &
   KINEMATICS_thermal_expansion_ID, &
   KINEMATICS_vacancy_strain_ID, &
   KINEMATICS_hydrogen_strain_ID, &
   STIFFNESS_DEGRADATION_damage_ID, &
   STIFFNESS_DEGRADATION_porosity_ID, &
   THERMAL_isothermal_ID, &
   THERMAL_adiabatic_ID, &
   THERMAL_conduction_ID, &
   DAMAGE_none_ID, &
   DAMAGE_local_ID, &
   DAMAGE_nonlocal_ID, &
   VACANCYFLUX_isoconc_ID, &
   VACANCYFLUX_isochempot_ID, &
   VACANCYFLUX_cahnhilliard_ID, &
   POROSITY_none_ID, &
   POROSITY_phasefield_ID, &
   HYDROGENFLUX_isoconc_ID, &
   HYDROGENFLUX_cahnhilliard_ID, &
   HOMOGENIZATION_none_ID, &
   HOMOGENIZATION_isostrain_ID, &
#ifdef HDF
   material_NconstituentsPhase, &
#endif
   HOMOGENIZATION_RGC_ID

 private :: &
   material_parseHomogenization, &
   material_parseMicrostructure, &
   material_parseCrystallite, &
   material_parsePhase, &
   material_parseTexture, &
   material_populateGrains

contains


!--------------------------------------------------------------------------------------------------
!> @brief parses material configuration file
!> @details figures out if solverJobName.materialConfig is present, if not looks for
!> material.config
!--------------------------------------------------------------------------------------------------
subroutine material_init()
 use, intrinsic :: iso_fortran_env                             ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use IO, only: &
   IO_error, &
   IO_open_file, &
   IO_open_jobFile_stat, &
   IO_timeStamp
 use debug, only: &
   debug_level, &
   debug_material, &
   debug_levelBasic, &
   debug_levelExtensive
 use mesh, only: &
   mesh_maxNips, &
   mesh_NcpElems, &
   mesh_element, &
   FE_Nips, &
   FE_geomtype
 use numerics, only: &
   worldrank

 implicit none
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt)            :: m,c,h, myDebug, myPhase, myHomog
 integer(pInt) :: &
  g, &                                                                                              !< grain number
  i, &                                                                                              !< integration point number
  e, &                                                                                              !< element number
  phase
 integer(pInt), dimension(:), allocatable :: ConstitutivePosition
 integer(pInt), dimension(:), allocatable :: CrystallitePosition
 integer(pInt), dimension(:), allocatable :: HomogenizationPosition

 myDebug = debug_level(debug_material)

 mainProcess: if (worldrank == 0) then
   write(6,'(/,a)') ' <<<+-  material init  -+>>>'
   write(6,'(a)') ' $Id$'
   write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ...open material.config file
 call material_parseHomogenization(FILEUNIT,material_partHomogenization)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Homogenization parsed'; flush(6)
 call material_parseMicrostructure(FILEUNIT,material_partMicrostructure)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Microstructure parsed'; flush(6)
 call material_parseCrystallite(FILEUNIT,material_partCrystallite)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Crystallite parsed'; flush(6)
 call material_parseTexture(FILEUNIT,material_partTexture)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Texture parsed'; flush(6)
 call material_parsePhase(FILEUNIT,material_partPhase)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Phase parsed'; flush(6)
 close(FILEUNIT)

 allocate(plasticState       (material_Nphase))
 allocate(sourceState        (material_Nphase))
 do myPhase = 1,material_Nphase
   allocate(sourceState(myPhase)%p(phase_Nsources(myPhase)))
 enddo

 allocate(homogState         (material_Nhomogenization))
 allocate(thermalState       (material_Nhomogenization))
 allocate(damageState        (material_Nhomogenization))
 allocate(vacancyfluxState   (material_Nhomogenization))
 allocate(porosityState      (material_Nhomogenization))
 allocate(hydrogenfluxState  (material_Nhomogenization))

 allocate(thermalMapping     (material_Nhomogenization))
 allocate(damageMapping      (material_Nhomogenization))
 allocate(vacancyfluxMapping (material_Nhomogenization))
 allocate(porosityMapping    (material_Nhomogenization))
 allocate(hydrogenfluxMapping(material_Nhomogenization))

 allocate(temperature        (material_Nhomogenization))
 allocate(damage             (material_Nhomogenization))
 allocate(vacancyConc        (material_Nhomogenization))
 allocate(porosity           (material_Nhomogenization))
 allocate(hydrogenConc       (material_Nhomogenization))

 allocate(temperatureRate    (material_Nhomogenization))
 allocate(vacancyConcRate    (material_Nhomogenization))
 allocate(hydrogenConcRate   (material_Nhomogenization))

 do m = 1_pInt,material_Nmicrostructure
   if(microstructure_crystallite(m) < 1_pInt .or. &
      microstructure_crystallite(m) > material_Ncrystallite) &
        call IO_error(150_pInt,m,ext_msg='crystallite')
   if(minval(microstructure_phase(1:microstructure_Nconstituents(m),m)) < 1_pInt .or. &
      maxval(microstructure_phase(1:microstructure_Nconstituents(m),m)) > material_Nphase) &
        call IO_error(150_pInt,m,ext_msg='phase')
   if(minval(microstructure_texture(1:microstructure_Nconstituents(m),m)) < 1_pInt .or. &
      maxval(microstructure_texture(1:microstructure_Nconstituents(m),m)) > material_Ntexture) &
        call IO_error(150_pInt,m,ext_msg='texture')
   if(microstructure_Nconstituents(m) < 1_pInt) &
        call IO_error(151_pInt,m)
 enddo

 debugOut: if (iand(myDebug,debug_levelExtensive) /= 0_pInt) then
   write(6,'(/,a,/)') ' MATERIAL configuration'
   write(6,'(a32,1x,a16,1x,a6)') 'homogenization                  ','type            ','grains'
   do h = 1_pInt,material_Nhomogenization
     write(6,'(1x,a32,1x,a16,1x,i6)') homogenization_name(h),homogenization_type(h),homogenization_Ngrains(h)
   enddo
   write(6,'(/,a14,18x,1x,a11,1x,a12,1x,a13)') 'microstructure','crystallite','constituents','homogeneous'
   do m = 1_pInt,material_Nmicrostructure
     write(6,'(1x,a32,1x,i11,1x,i12,1x,l13)') microstructure_name(m), &
                                        microstructure_crystallite(m), &
                                        microstructure_Nconstituents(m), &
                                        microstructure_elemhomo(m)
     if (microstructure_Nconstituents(m) > 0_pInt) then
       do c = 1_pInt,microstructure_Nconstituents(m)
         write(6,'(a1,1x,a32,1x,a32,1x,f7.4)') '>',phase_name(microstructure_phase(c,m)),&
                                                   texture_name(microstructure_texture(c,m)),&
                                                   microstructure_fraction(c,m)
       enddo
       write(6,*)
     endif
   enddo
 endif debugOut

 call material_populateGrains

 allocate(phaseAt                   (  homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems),source=0_pInt)
 allocate(phasememberAt             (  homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems),source=0_pInt)
 allocate(mappingHomogenization     (2,                          mesh_maxNips,mesh_NcpElems),source=0_pInt)
 allocate(mappingCrystallite        (2,homogenization_maxNgrains,             mesh_NcpElems),source=0_pInt)
 allocate(mappingHomogenizationConst(                            mesh_maxNips,mesh_NcpElems),source=1_pInt)

 allocate(ConstitutivePosition  (material_Nphase),         source=0_pInt)
 allocate(HomogenizationPosition(material_Nhomogenization),source=0_pInt)
 allocate(CrystallitePosition   (material_Nphase),         source=0_pInt)

 ElemLoop:do e = 1_pInt,mesh_NcpElems                                                               ! loop over elements
 myHomog = mesh_element(3,e)
   IPloop:do i = 1_pInt,FE_Nips(FE_geomtype(mesh_element(2,e)))                                     ! loop over IPs
     HomogenizationPosition(myHomog) = HomogenizationPosition(myHomog) + 1_pInt
     mappingHomogenization(1:2,i,e) = [HomogenizationPosition(myHomog),myHomog]
     GrainLoop:do g = 1_pInt,homogenization_Ngrains(mesh_element(3,e))                              ! loop over grains
       phase = material_phase(g,i,e)
       ConstitutivePosition(phase) = ConstitutivePosition(phase)+1_pInt                             ! not distinguishing between instances of same phase
       phaseAt(g,i,e)              = phase
       phasememberAt(g,i,e)        = ConstitutivePosition(phase)
     enddo GrainLoop
   enddo IPloop
 enddo ElemLoop

! hack needed to initialize field values used during constitutive and crystallite initializations
 do myHomog = 1,material_Nhomogenization
   thermalMapping     (myHomog)%p => mappingHomogenizationConst
   damageMapping      (myHomog)%p => mappingHomogenizationConst
   vacancyfluxMapping (myHomog)%p => mappingHomogenizationConst
   porosityMapping    (myHomog)%p => mappingHomogenizationConst
   hydrogenfluxMapping(myHomog)%p => mappingHomogenizationConst
   allocate(temperature     (myHomog)%p(1), source=thermal_initialT(myHomog))
   allocate(damage          (myHomog)%p(1), source=damage_initialPhi(myHomog))
   allocate(vacancyConc     (myHomog)%p(1), source=vacancyflux_initialCv(myHomog))
   allocate(porosity        (myHomog)%p(1), source=porosity_initialPhi(myHomog))
   allocate(hydrogenConc    (myHomog)%p(1), source=hydrogenflux_initialCh(myHomog))
   allocate(temperatureRate (myHomog)%p(1), source=0.0_pReal)
   allocate(vacancyConcRate (myHomog)%p(1), source=0.0_pReal)
   allocate(hydrogenConcRate(myHomog)%p(1), source=0.0_pReal)
 enddo

end subroutine material_init


!--------------------------------------------------------------------------------------------------
!> @brief parses the homogenization part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseHomogenization(fileUnit,myPart)
 use IO, only: &
   IO_read, &
   IO_globalTagInPart, &
   IO_countSections, &
   IO_error, &
   IO_countTagInPart, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringValue, &
   IO_intValue, &
   IO_floatValue, &
   IO_stringPos, &
   IO_EOF
 use mesh, only: &
   mesh_element

 implicit none
 character(len=*), intent(in) :: myPart
 integer(pInt),    intent(in) :: fileUnit


 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt)        :: Nsections, section, s, p
 character(len=65536) :: &
   tag, line
 logical              :: echo

 echo = IO_globalTagInPart(fileUnit,myPart,'/echo/')
 Nsections = IO_countSections(fileUnit,myPart)
 material_Nhomogenization = Nsections
 if (Nsections < 1_pInt) call IO_error(160_pInt,ext_msg=myPart)

 allocate(homogenization_name(Nsections));          homogenization_name = ''
 allocate(homogenization_type(Nsections),           source=HOMOGENIZATION_undefined_ID)
 allocate(thermal_type(Nsections),                  source=THERMAL_isothermal_ID)
 allocate(damage_type (Nsections),                  source=DAMAGE_none_ID)
 allocate(vacancyflux_type(Nsections),              source=VACANCYFLUX_isoconc_ID)
 allocate(porosity_type (Nsections),                source=POROSITY_none_ID)
 allocate(hydrogenflux_type(Nsections),             source=HYDROGENFLUX_isoconc_ID)
 allocate(homogenization_typeInstance(Nsections),   source=0_pInt)
 allocate(thermal_typeInstance(Nsections),          source=0_pInt)
 allocate(damage_typeInstance(Nsections),           source=0_pInt)
 allocate(vacancyflux_typeInstance(Nsections),      source=0_pInt)
 allocate(porosity_typeInstance(Nsections),         source=0_pInt)
 allocate(hydrogenflux_typeInstance(Nsections),     source=0_pInt)
 allocate(homogenization_Ngrains(Nsections),        source=0_pInt)
 allocate(homogenization_Noutput(Nsections),        source=0_pInt)
 allocate(homogenization_active(Nsections),         source=.false.)  !!!!!!!!!!!!!!!
 allocate(thermal_initialT(Nsections),              source=300.0_pReal)
 allocate(damage_initialPhi(Nsections),             source=1.0_pReal)
 allocate(vacancyflux_initialCv(Nsections),         source=0.0_pReal)
 allocate(porosity_initialPhi(Nsections),           source=1.0_pReal)
 allocate(hydrogenflux_initialCh(Nsections),        source=0.0_pReal)

 forall (s = 1_pInt:Nsections) homogenization_active(s) = any(mesh_element(3,:) == s)               ! current homogenization used in model? Homogenization view, maximum operations depend on maximum number of homog schemes
   homogenization_Noutput = IO_countTagInPart(fileUnit,myPart,'(output)',Nsections)

 rewind(fileUnit)
 line        = ''                                                                                   ! to have it initialized
 section     = 0_pInt                                                                               !  - " -
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= myPart)                     ! wind forward to <homogenization>
   line = IO_read(fileUnit)
 enddo
 if (echo) write(6,'(/,1x,a)') trim(line)                                                           ! echo part header

 do while (trim(line) /= IO_EOF)                                                                    ! read through sections of material part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit
   endif
   if (echo) write(6,'(2x,a)') trim(line)                                                           ! echo back read lines
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     homogenization_name(section) = IO_getTag(line,'[',']')
   endif
   if (section > 0_pInt) then
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('type')
         select case (IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case(HOMOGENIZATION_NONE_label)
             homogenization_type(section) = HOMOGENIZATION_NONE_ID
             homogenization_Ngrains(section) = 1_pInt
           case(HOMOGENIZATION_ISOSTRAIN_label)
             homogenization_type(section) = HOMOGENIZATION_ISOSTRAIN_ID
           case(HOMOGENIZATION_RGC_label)
             homogenization_type(section) = HOMOGENIZATION_RGC_ID
           case default
             call IO_error(500_pInt,ext_msg=trim(IO_stringValue(line,chunkPos,2_pInt)))
         end select
         homogenization_typeInstance(section) = &
                                          count(homogenization_type==homogenization_type(section))  ! count instances
        case ('thermal')
         select case (IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case(THERMAL_isothermal_label)
             thermal_type(section) = THERMAL_isothermal_ID
           case(THERMAL_adiabatic_label)
             thermal_type(section) = THERMAL_adiabatic_ID
           case(THERMAL_conduction_label)
             thermal_type(section) = THERMAL_conduction_ID
           case default
             call IO_error(500_pInt,ext_msg=trim(IO_stringValue(line,chunkPos,2_pInt)))
         end select

        case ('damage')
         select case (IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case(DAMAGE_NONE_label)
             damage_type(section) = DAMAGE_none_ID
           case(DAMAGE_LOCAL_label)
             damage_type(section) = DAMAGE_local_ID
           case(DAMAGE_NONLOCAL_label)
             damage_type(section) = DAMAGE_nonlocal_ID
           case default
             call IO_error(500_pInt,ext_msg=trim(IO_stringValue(line,chunkPos,2_pInt)))
         end select

        case ('vacancyflux')
         select case (IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case(VACANCYFLUX_isoconc_label)
             vacancyflux_type(section) = VACANCYFLUX_isoconc_ID
           case(VACANCYFLUX_isochempot_label)
             vacancyflux_type(section) = VACANCYFLUX_isochempot_ID
           case(VACANCYFLUX_cahnhilliard_label)
             vacancyflux_type(section) = VACANCYFLUX_cahnhilliard_ID
           case default
             call IO_error(500_pInt,ext_msg=trim(IO_stringValue(line,chunkPos,2_pInt)))
         end select

        case ('porosity')
         select case (IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case(POROSITY_NONE_label)
             porosity_type(section) = POROSITY_none_ID
           case(POROSITY_phasefield_label)
             porosity_type(section) = POROSITY_phasefield_ID
           case default
             call IO_error(500_pInt,ext_msg=trim(IO_stringValue(line,chunkPos,2_pInt)))
         end select

        case ('hydrogenflux')
         select case (IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case(HYDROGENFLUX_isoconc_label)
             hydrogenflux_type(section) = HYDROGENFLUX_isoconc_ID
           case(HYDROGENFLUX_cahnhilliard_label)
             hydrogenflux_type(section) = HYDROGENFLUX_cahnhilliard_ID
           case default
             call IO_error(500_pInt,ext_msg=trim(IO_stringValue(line,chunkPos,2_pInt)))
         end select

       case ('nconstituents','ngrains')
         homogenization_Ngrains(section) = IO_intValue(line,chunkPos,2_pInt)

       case ('initialtemperature','initialt')
         thermal_initialT(section) = IO_floatValue(line,chunkPos,2_pInt)

       case ('initialdamage')
         damage_initialPhi(section) = IO_floatValue(line,chunkPos,2_pInt)

       case ('initialvacancyconc','initialcv')
         vacancyflux_initialCv(section) = IO_floatValue(line,chunkPos,2_pInt)

       case ('initialporosity')
         porosity_initialPhi(section) = IO_floatValue(line,chunkPos,2_pInt)

       case ('initialhydrogenconc','initialch')
         hydrogenflux_initialCh(section) = IO_floatValue(line,chunkPos,2_pInt)

     end select
   endif
 enddo

 do p=1_pInt, Nsections
   homogenization_typeInstance(p)  = count(homogenization_type(1:p)  == homogenization_type(p))
   thermal_typeInstance(p)         = count(thermal_type       (1:p)  == thermal_type       (p))
   damage_typeInstance(p)          = count(damage_type        (1:p)  == damage_type        (p))
   vacancyflux_typeInstance(p)     = count(vacancyflux_type   (1:p)  == vacancyflux_type   (p))
   porosity_typeInstance(p)        = count(porosity_type      (1:p)  == porosity_type      (p))
   hydrogenflux_typeInstance(p)    = count(hydrogenflux_type  (1:p)  == hydrogenflux_type  (p))
 enddo

 homogenization_maxNgrains = maxval(homogenization_Ngrains,homogenization_active)

end subroutine material_parseHomogenization


!--------------------------------------------------------------------------------------------------
!> @brief parses the microstructure part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseMicrostructure(fileUnit,myPart)
 use IO
 use mesh, only: &
   mesh_element, &
   mesh_NcpElems

 implicit none
 character(len=*), intent(in) :: myPart
 integer(pInt),    intent(in) :: fileUnit


 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: Nsections, section, constituent, e, i
 character(len=65536) :: &
   tag, line
 logical              :: echo

 echo = IO_globalTagInPart(fileUnit,myPart,'/echo/')

 Nsections = IO_countSections(fileUnit,myPart)
 material_Nmicrostructure = Nsections
 if (Nsections < 1_pInt) call IO_error(160_pInt,ext_msg=myPart)

 allocate(microstructure_name(Nsections));            microstructure_name = ''
 allocate(microstructure_crystallite(Nsections),          source=0_pInt)
 allocate(microstructure_Nconstituents(Nsections),        source=0_pInt)
 allocate(microstructure_active(Nsections),               source=.false.)
 allocate(microstructure_elemhomo(Nsections),             source=.false.)

 if(any(mesh_element(4,1:mesh_NcpElems) > Nsections)) &
  call IO_error(155_pInt,ext_msg='Microstructure in geometry > Sections in material.config')

 forall (e = 1_pInt:mesh_NcpElems) microstructure_active(mesh_element(4,e)) = .true.                ! current microstructure used in model? Elementwise view, maximum N operations for N elements

 microstructure_Nconstituents = IO_countTagInPart(fileUnit,myPart,'(constituent)',Nsections)
 microstructure_maxNconstituents = maxval(microstructure_Nconstituents)
 microstructure_elemhomo = IO_spotTagInPart(fileUnit,myPart,'/elementhomogeneous/',Nsections)

 allocate(microstructure_phase   (microstructure_maxNconstituents,Nsections),source=0_pInt)
 allocate(microstructure_texture (microstructure_maxNconstituents,Nsections),source=0_pInt)
 allocate(microstructure_fraction(microstructure_maxNconstituents,Nsections),source=0.0_pReal)

 rewind(fileUnit)
 line        = ''                                                                                   ! to have it initialized
 section     = 0_pInt                                                                               !  - " -
 constituent = 0_pInt                                                                               !  - " -

 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= myPart)                     ! wind forward to <microstructure>
   line = IO_read(fileUnit)
 enddo
 if (echo) write(6,'(/,1x,a)') trim(line)                                                           ! echo part header

 do while (trim(line) /= IO_EOF)                                                                    ! read through sections of material part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit
   endif
   if (echo) write(6,'(2x,a)') trim(line)                                                           ! echo back read lines
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     constituent = 0_pInt
     microstructure_name(section) = IO_getTag(line,'[',']')
   endif
   if (section > 0_pInt) then
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('crystallite')
         microstructure_crystallite(section) = IO_intValue(line,chunkPos,2_pInt)
       case ('(constituent)')
         constituent = constituent + 1_pInt
         do i=2_pInt,6_pInt,2_pInt
           tag = IO_lc(IO_stringValue(line,chunkPos,i))
           select case (tag)
             case('phase')
               microstructure_phase(constituent,section) =    IO_intValue(line,chunkPos,i+1_pInt)
             case('texture')
               microstructure_texture(constituent,section) =  IO_intValue(line,chunkPos,i+1_pInt)
             case('fraction')
               microstructure_fraction(constituent,section) = IO_floatValue(line,chunkPos,i+1_pInt)
           end select
         enddo
     end select
   endif
 enddo

end subroutine material_parseMicrostructure


!--------------------------------------------------------------------------------------------------
!> @brief parses the crystallite part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseCrystallite(fileUnit,myPart)
 use IO, only: &
   IO_read, &
   IO_countSections, &
   IO_error, &
   IO_countTagInPart, &
   IO_globalTagInPart, &
   IO_getTag, &
   IO_lc, &
   IO_isBlank, &
   IO_EOF

 implicit none
 character(len=*), intent(in) :: myPart
 integer(pInt),    intent(in) :: fileUnit

 integer(pInt)        :: Nsections, &
                         section
 character(len=65536) :: line
 logical              :: echo

 echo = IO_globalTagInPart(fileUnit,myPart,'/echo/')

 Nsections = IO_countSections(fileUnit,myPart)
 material_Ncrystallite = Nsections
 if (Nsections < 1_pInt) call IO_error(160_pInt,ext_msg=myPart)

 allocate(crystallite_name(Nsections));       crystallite_name = ''
 allocate(crystallite_Noutput(Nsections),           source=0_pInt)

 crystallite_Noutput = IO_countTagInPart(fileUnit,myPart,'(output)',Nsections)

 rewind(fileUnit)
 line        = ''                                                                                   ! to have it initialized
 section     = 0_pInt                                                                               !  - " -
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= myPart)                     ! wind forward to <Crystallite>
   line = IO_read(fileUnit)
 enddo
 if (echo) write(6,'(/,1x,a)') trim(line)                                                           ! echo part header

 do while (trim(line) /= IO_EOF)                                                                    ! read through sections of material part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit
   endif
   if (echo) write(6,'(2x,a)') trim(line)                                                           ! echo back read lines
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     crystallite_name(section) = IO_getTag(line,'[',']')
   endif
 enddo

end subroutine material_parseCrystallite


!--------------------------------------------------------------------------------------------------
!> @brief parses the phase part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parsePhase(fileUnit,myPart)
 use IO, only: &
   IO_read, &
   IO_globalTagInPart, &
   IO_countSections, &
   IO_error, &
   IO_countTagInPart, &
   IO_getTag, &
   IO_spotTagInPart, &
   IO_lc, &
   IO_isBlank, &
   IO_stringValue, &
   IO_stringPos, &
   IO_EOF

 implicit none
 character(len=*), intent(in) :: myPart
 integer(pInt),    intent(in) :: fileUnit


 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: Nsections, section, sourceCtr, kinematicsCtr, stiffDegradationCtr, p
 character(len=65536) :: &
  tag,line
 logical              :: echo

 echo = IO_globalTagInPart(fileUnit,myPart,'/echo/')

 Nsections = IO_countSections(fileUnit,myPart)
 material_Nphase = Nsections
 if (Nsections < 1_pInt) call IO_error(160_pInt,ext_msg=myPart)

 allocate(phase_name(Nsections));                phase_name = ''
 allocate(phase_elasticity(Nsections),           source=ELASTICITY_undefined_ID)
 allocate(phase_elasticityInstance(Nsections),   source=0_pInt)
 allocate(phase_plasticity(Nsections) ,          source=PLASTICITY_undefined_ID)
 allocate(phase_plasticityInstance(Nsections),   source=0_pInt)
 allocate(phase_Nsources(Nsections),             source=0_pInt)
 allocate(phase_Nkinematics(Nsections),          source=0_pInt)
 allocate(phase_NstiffnessDegradations(Nsections),source=0_pInt)
 allocate(phase_Noutput(Nsections),              source=0_pInt)
 allocate(phase_localPlasticity(Nsections),      source=.false.)

 phase_Noutput = IO_countTagInPart(fileUnit,myPart,'(output)',Nsections)
 phase_Nsources = IO_countTagInPart(fileUnit,myPart,'(source)',Nsections)
 phase_Nkinematics = IO_countTagInPart(fileUnit,myPart,'(kinematics)',Nsections)
 phase_NstiffnessDegradations = IO_countTagInPart(fileUnit,myPart,'(stiffness_degradation)',Nsections)
 phase_localPlasticity = .not. IO_spotTagInPart(fileUnit,myPart,'/nonlocal/',Nsections)

 allocate(phase_source(maxval(phase_Nsources),Nsections), source=SOURCE_undefined_ID)
 allocate(phase_kinematics(maxval(phase_Nkinematics),Nsections), source=KINEMATICS_undefined_ID)
 allocate(phase_stiffnessDegradation(maxval(phase_NstiffnessDegradations),Nsections), &
          source=STIFFNESS_DEGRADATION_undefined_ID)

 rewind(fileUnit)
 line        = ''                                                                                   ! to have it initialized
 section     = 0_pInt                                                                               !  - " -
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= myPart)                     ! wind forward to <Phase>
   line = IO_read(fileUnit)
 enddo
 if (echo) write(6,'(/,1x,a)') trim(line)                                                           ! echo part header

 do while (trim(line) /= IO_EOF)                                                                    ! read through sections of material part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit
   endif
   if (echo) write(6,'(2x,a)') trim(line)                                                           ! echo back read lines
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     sourceCtr = 0_pInt
     kinematicsCtr = 0_pInt
     stiffDegradationCtr = 0_pInt
     phase_name(section) = IO_getTag(line,'[',']')
   endif
   if (section > 0_pInt) then
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('elasticity')
         select case (IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case (ELASTICITY_HOOKE_label)
             phase_elasticity(section) = ELASTICITY_HOOKE_ID
           case default
             call IO_error(200_pInt,ext_msg=trim(IO_stringValue(line,chunkPos,2_pInt)))
         end select
       case ('plasticity')
         select case (IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case (PLASTICITY_NONE_label)
             phase_plasticity(section) = PLASTICITY_NONE_ID
           case (PLASTICITY_ISOTROPIC_label)
             phase_plasticity(section) = PLASTICITY_ISOTROPIC_ID
           case (PLASTICITY_J2_label)
             phase_plasticity(section) = PLASTICITY_J2_ID
           case (PLASTICITY_PHENOPOWERLAW_label)
             phase_plasticity(section) = PLASTICITY_PHENOPOWERLAW_ID
           case (PLASTICITY_PHENOPLUS_label)
             phase_plasticity(section) = PLASTICITY_PHENOPLUS_ID
           case (PLASTICITY_DISLOTWIN_label)
             phase_plasticity(section) = PLASTICITY_DISLOTWIN_ID
           case (PLASTICITY_DISLOUCLA_label)
             phase_plasticity(section) = PLASTICITY_DISLOUCLA_ID
           case (PLASTICITY_TITANMOD_label)
             phase_plasticity(section) = PLASTICITY_TITANMOD_ID
           case (PLASTICITY_NONLOCAL_label)
             phase_plasticity(section) = PLASTICITY_NONLOCAL_ID
           case default
             call IO_error(201_pInt,ext_msg=trim(IO_stringValue(line,chunkPos,2_pInt)))
         end select
       case ('(source)')
         sourceCtr = sourceCtr + 1_pInt
         select case (IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case (SOURCE_thermal_dissipation_label)
             phase_source(sourceCtr,section) = SOURCE_thermal_dissipation_ID
           case (SOURCE_thermal_externalheat_label)
             phase_source(sourceCtr,section) = SOURCE_thermal_externalheat_ID
           case (SOURCE_damage_isoBrittle_label)
             phase_source(sourceCtr,section) = SOURCE_damage_isoBrittle_ID
           case (SOURCE_damage_isoDuctile_label)
             phase_source(sourceCtr,section) = SOURCE_damage_isoDuctile_ID
           case (SOURCE_damage_anisoBrittle_label)
             phase_source(sourceCtr,section) = SOURCE_damage_anisoBrittle_ID
           case (SOURCE_damage_anisoDuctile_label)
             phase_source(sourceCtr,section) = SOURCE_damage_anisoDuctile_ID
           case (SOURCE_vacancy_phenoplasticity_label)
             phase_source(sourceCtr,section) = SOURCE_vacancy_phenoplasticity_ID
           case (SOURCE_vacancy_irradiation_label)
             phase_source(sourceCtr,section) = SOURCE_vacancy_irradiation_ID
           case (SOURCE_vacancy_thermalfluc_label)
             phase_source(sourceCtr,section) = SOURCE_vacancy_thermalfluc_ID
         end select
       case ('(kinematics)')
         kinematicsCtr = kinematicsCtr + 1_pInt
         select case (IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case (KINEMATICS_cleavage_opening_label)
             phase_kinematics(kinematicsCtr,section) = KINEMATICS_cleavage_opening_ID
           case (KINEMATICS_slipplane_opening_label)
             phase_kinematics(kinematicsCtr,section) = KINEMATICS_slipplane_opening_ID
           case (KINEMATICS_thermal_expansion_label)
             phase_kinematics(kinematicsCtr,section) = KINEMATICS_thermal_expansion_ID
           case (KINEMATICS_vacancy_strain_label)
             phase_kinematics(kinematicsCtr,section) = KINEMATICS_vacancy_strain_ID
           case (KINEMATICS_hydrogen_strain_label)
             phase_kinematics(kinematicsCtr,section) = KINEMATICS_hydrogen_strain_ID
         end select
       case ('(stiffness_degradation)')
         stiffDegradationCtr = stiffDegradationCtr + 1_pInt
         select case (IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case (STIFFNESS_DEGRADATION_damage_label)
             phase_stiffnessDegradation(stiffDegradationCtr,section) = STIFFNESS_DEGRADATION_damage_ID
           case (STIFFNESS_DEGRADATION_porosity_label)
             phase_stiffnessDegradation(stiffDegradationCtr,section) = STIFFNESS_DEGRADATION_porosity_ID
         end select

     end select
   endif
 enddo

 do p=1_pInt, Nsections
   phase_elasticityInstance(p)  = count(phase_elasticity(1:p)  == phase_elasticity(p))
   phase_plasticityInstance(p)  = count(phase_plasticity(1:p)  == phase_plasticity(p))
 enddo

end subroutine material_parsePhase

!--------------------------------------------------------------------------------------------------
!> @brief parses the texture part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseTexture(fileUnit,myPart)
 use IO, only: &
   IO_read, &
   IO_globalTagInPart, &
   IO_countSections, &
   IO_error, &
   IO_countTagInPart, &
   IO_getTag, &
   IO_spotTagInPart, &
   IO_lc, &
   IO_isBlank, &
   IO_floatValue, &
   IO_stringValue, &
   IO_stringPos, &
   IO_EOF
 use math, only: &
   inRad, &
   math_sampleRandomOri, &
   math_I3, &
   math_inv33

 implicit none
 character(len=*), intent(in) :: myPart
 integer(pInt),    intent(in) :: fileUnit


 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: Nsections, section, gauss, fiber, j
 character(len=65536) :: tag
 character(len=65536) :: line
 logical              :: echo

 echo = IO_globalTagInPart(fileUnit,myPart,'/echo/')

 Nsections = IO_countSections(fileUnit,myPart)
 material_Ntexture = Nsections
 if (Nsections < 1_pInt) call IO_error(160_pInt,ext_msg=myPart)

 allocate(texture_name(Nsections));        texture_name=''
 allocate(texture_ODFfile(Nsections));  texture_ODFfile=''
 allocate(texture_symmetry(Nsections),           source=1_pInt)
 allocate(texture_Ngauss(Nsections),             source=0_pInt)
 allocate(texture_Nfiber(Nsections),             source=0_pInt)

 texture_Ngauss = IO_countTagInPart(fileUnit,myPart,'(gauss)', Nsections) + &
                  IO_countTagInPart(fileUnit,myPart,'(random)',Nsections)
 texture_Nfiber = IO_countTagInPart(fileUnit,myPart,'(fiber)', Nsections)
 texture_maxNgauss = maxval(texture_Ngauss)
 texture_maxNfiber = maxval(texture_Nfiber)
 allocate(texture_Gauss (5,texture_maxNgauss,Nsections), source=0.0_pReal)
 allocate(texture_Fiber (6,texture_maxNfiber,Nsections), source=0.0_pReal)
 allocate(texture_transformation(3,3,Nsections),         source=0.0_pReal)
 texture_transformation = spread(math_I3,3,Nsections)

 rewind(fileUnit)
 line    = ''                                                                                       ! to have in initialized
 section = 0_pInt                                                                                   ! - " -
 gauss   = 0_pInt                                                                                   ! - " -
 fiber   = 0_pInt                                                                                   ! - " -
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= myPart)                     ! wind forward to <texture>
   line = IO_read(fileUnit)
 enddo
 if (echo) write(6,'(/,1x,a)') trim(line)                                                           ! echo part header

 do while (trim(line) /= IO_EOF)
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit
   endif
   if (echo) write(6,'(2x,a)') trim(line)                                                           ! echo back read lines
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     gauss = 0_pInt
     fiber = 0_pInt
     texture_name(section) = IO_getTag(line,'[',']')
   endif
   if (section > 0_pInt) then
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     textureType: select case(tag)

       case ('axes', 'rotation') textureType
         do j = 1_pInt, 3_pInt                                                                      ! look for "x", "y", and "z" entries
           tag = IO_lc(IO_stringValue(line,chunkPos,j+1_pInt))
           select case (tag)
             case('x', '+x')
               texture_transformation(j,1:3,section) = [ 1.0_pReal, 0.0_pReal, 0.0_pReal]           ! original axis is now +x-axis
             case('-x')
               texture_transformation(j,1:3,section) = [-1.0_pReal, 0.0_pReal, 0.0_pReal]           ! original axis is now -x-axis
             case('y', '+y')
               texture_transformation(j,1:3,section) = [ 0.0_pReal, 1.0_pReal, 0.0_pReal]           ! original axis is now +y-axis
             case('-y')
               texture_transformation(j,1:3,section) = [ 0.0_pReal,-1.0_pReal, 0.0_pReal]           ! original axis is now -y-axis
             case('z', '+z')
               texture_transformation(j,1:3,section) = [ 0.0_pReal, 0.0_pReal, 1.0_pReal]           ! original axis is now +z-axis
             case('-z')
               texture_transformation(j,1:3,section) = [ 0.0_pReal, 0.0_pReal,-1.0_pReal]           ! original axis is now -z-axis
             case default
               call IO_error(157_pInt,section)
           end select
         enddo

       case ('hybridia') textureType
         texture_ODFfile(section) = IO_stringValue(line,chunkPos,2_pInt)

       case ('symmetry') textureType
         tag = IO_lc(IO_stringValue(line,chunkPos,2_pInt))
         select case (tag)
           case('orthotropic')
             texture_symmetry(section) = 4_pInt
           case('monoclinic')
             texture_symmetry(section) = 2_pInt
           case default
             texture_symmetry(section) = 1_pInt
         end select

       case ('(random)') textureType
         gauss = gauss + 1_pInt
         texture_Gauss(1:3,gauss,section) = math_sampleRandomOri()
         do j = 2_pInt,4_pInt,2_pInt
           tag = IO_lc(IO_stringValue(line,chunkPos,j))
           select case (tag)
             case('scatter')
                 texture_Gauss(4,gauss,section) = IO_floatValue(line,chunkPos,j+1_pInt)*inRad
             case('fraction')
                 texture_Gauss(5,gauss,section) = IO_floatValue(line,chunkPos,j+1_pInt)
           end select
         enddo

       case ('(gauss)') textureType
         gauss = gauss + 1_pInt
         do j = 2_pInt,10_pInt,2_pInt
           tag = IO_lc(IO_stringValue(line,chunkPos,j))
           select case (tag)
             case('phi1')
                 texture_Gauss(1,gauss,section) = IO_floatValue(line,chunkPos,j+1_pInt)*inRad
             case('phi')
                 texture_Gauss(2,gauss,section) = IO_floatValue(line,chunkPos,j+1_pInt)*inRad
             case('phi2')
                 texture_Gauss(3,gauss,section) = IO_floatValue(line,chunkPos,j+1_pInt)*inRad
             case('scatter')
                 texture_Gauss(4,gauss,section) = IO_floatValue(line,chunkPos,j+1_pInt)*inRad
             case('fraction')
                 texture_Gauss(5,gauss,section) = IO_floatValue(line,chunkPos,j+1_pInt)
           end select
         enddo

       case ('(fiber)') textureType
         fiber = fiber + 1_pInt
         do j = 2_pInt,12_pInt,2_pInt
           tag = IO_lc(IO_stringValue(line,chunkPos,j))
           select case (tag)
             case('alpha1')
                 texture_Fiber(1,fiber,section) = IO_floatValue(line,chunkPos,j+1_pInt)*inRad
             case('alpha2')
                 texture_Fiber(2,fiber,section) = IO_floatValue(line,chunkPos,j+1_pInt)*inRad
             case('beta1')
                 texture_Fiber(3,fiber,section) = IO_floatValue(line,chunkPos,j+1_pInt)*inRad
             case('beta2')
                 texture_Fiber(4,fiber,section) = IO_floatValue(line,chunkPos,j+1_pInt)*inRad
             case('scatter')
                 texture_Fiber(5,fiber,section) = IO_floatValue(line,chunkPos,j+1_pInt)*inRad
             case('fraction')
                 texture_Fiber(6,fiber,section) = IO_floatValue(line,chunkPos,j+1_pInt)
           end select
         enddo

     end select textureType
   endif
 enddo

end subroutine material_parseTexture


!--------------------------------------------------------------------------------------------------
!> @brief populates the grains
!> @details populates the grains by identifying active microstructure/homogenization pairs,
!! calculates the volume of the grains and deals with texture components and hybridIA
!--------------------------------------------------------------------------------------------------
subroutine material_populateGrains
 use math, only: &
   math_RtoEuler, &
   math_EulerToR, &
   math_mul33x33, &
   math_range, &
   math_sampleRandomOri, &
   math_sampleGaussOri, &
   math_sampleFiberOri, &
   math_symmetricEulers
 use mesh, only: &
   mesh_element, &
   mesh_maxNips, &
   mesh_NcpElems, &
   mesh_ipVolume, &
   FE_Nips, &
   FE_geomtype
 use IO, only: &
   IO_error, &
   IO_hybridIA
 use debug, only: &
   debug_level, &
   debug_material, &
   debug_levelBasic

 implicit none
 integer(pInt), dimension (:,:), allocatable :: Ngrains
 integer(pInt), dimension (microstructure_maxNconstituents)  :: &
   NgrainsOfConstituent, &
   currentGrainOfConstituent, &
   randomOrder
 real(pReal), dimension (microstructure_maxNconstituents)  :: &
   rndArray
 real(pReal), dimension (:),     allocatable :: volumeOfGrain
 real(pReal), dimension (:,:),   allocatable :: orientationOfGrain
 real(pReal), dimension (3)                  :: orientation
 real(pReal), dimension (3,3)                :: symOrientation
 integer(pInt), dimension (:),   allocatable :: phaseOfGrain, textureOfGrain
 integer(pInt) :: t,e,i,g,j,m,c,r,homog,micro,sgn,hme, myDebug, &
                  phaseID,textureID,dGrains,myNgrains,myNorientations,myNconstituents, &
                  grain,constituentGrain,ipGrain,symExtension, ip
 real(pReal) :: extreme,rnd
 integer(pInt),  dimension (:,:),   allocatable :: Nelems                                           ! counts number of elements in homog, micro array
 type(p_intvec), dimension (:,:), allocatable :: elemsOfHomogMicro                                  ! lists element number in homog, micro array

 myDebug = debug_level(debug_material)

 allocate(material_volume(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems),       source=0.0_pReal)
 allocate(material_phase(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems),        source=0_pInt)
 allocate(material_homog(mesh_maxNips,mesh_NcpElems),                                  source=0_pInt)
 allocate(material_texture(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems),      source=0_pInt)
 allocate(material_EulerAngles(3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems),source=0.0_pReal)

 allocate(Ngrains(material_Nhomogenization,material_Nmicrostructure),                  source=0_pInt)
 allocate(Nelems(material_Nhomogenization,material_Nmicrostructure),                   source=0_pInt)

! populating homogenization schemes in each
!--------------------------------------------------------------------------------------------------
 do e = 1_pInt, mesh_NcpElems
   material_homog(1_pInt:FE_Nips(FE_geomtype(mesh_element(2,e))),e) = mesh_element(3,e)
 enddo

!--------------------------------------------------------------------------------------------------
! precounting of elements for each homog/micro pair
 do e = 1_pInt, mesh_NcpElems
   homog = mesh_element(3,e)
   micro = mesh_element(4,e)
   Nelems(homog,micro) = Nelems(homog,micro) + 1_pInt
 enddo
 allocate(elemsOfHomogMicro(material_Nhomogenization,material_Nmicrostructure))
 do homog = 1,material_Nhomogenization
   do micro = 1,material_Nmicrostructure
     if (Nelems(homog,micro) > 0_pInt) then
       allocate(elemsOfHomogMicro(homog,micro)%p(Nelems(homog,micro)))
       elemsOfHomogMicro(homog,micro)%p = 0_pInt
    endif
   enddo
 enddo

!--------------------------------------------------------------------------------------------------
! identify maximum grain count per IP (from element) and find grains per homog/micro pair
 Nelems = 0_pInt                                                                                    ! reuse as counter
 elementLooping: do e = 1_pInt,mesh_NcpElems
   t    = FE_geomtype(mesh_element(2,e))
   homog = mesh_element(3,e)
   micro = mesh_element(4,e)
   if (homog < 1_pInt .or. homog > material_Nhomogenization) &                                      ! out of bounds
     call IO_error(154_pInt,e,0_pInt,0_pInt)
   if (micro < 1_pInt .or. micro > material_Nmicrostructure) &                                      ! out of bounds
     call IO_error(155_pInt,e,0_pInt,0_pInt)
   if (microstructure_elemhomo(micro)) then                                                         ! how many grains are needed at this element?
     dGrains = homogenization_Ngrains(homog)                                                        ! only one set of Ngrains (other IPs are plain copies)
   else
     dGrains = homogenization_Ngrains(homog) * FE_Nips(t)                                           ! each IP has Ngrains
   endif
   Ngrains(homog,micro) = Ngrains(homog,micro) + dGrains                                            ! total grain count
   Nelems(homog,micro)  = Nelems(homog,micro) + 1_pInt                                              ! total element count
   elemsOfHomogMicro(homog,micro)%p(Nelems(homog,micro)) = e                                        ! remember elements active in this homog/micro pair
 enddo elementLooping

 allocate(volumeOfGrain(maxval(Ngrains)),       source=0.0_pReal)                                   ! reserve memory for maximum case
 allocate(phaseOfGrain(maxval(Ngrains)),        source=0_pInt)                                      ! reserve memory for maximum case
 allocate(textureOfGrain(maxval(Ngrains)),      source=0_pInt)                                      ! reserve memory for maximum case
 allocate(orientationOfGrain(3,maxval(Ngrains)),source=0.0_pReal)                                   ! reserve memory for maximum case

 if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
     write(6,'(/,a/)') ' MATERIAL grain population'
     write(6,'(a32,1x,a32,1x,a6)') 'homogenization_name','microstructure_name','grain#'
   !$OMP END CRITICAL (write2out)
 endif
 do homog = 1_pInt,material_Nhomogenization                                                         ! loop over homogenizations
   dGrains = homogenization_Ngrains(homog)                                                          ! grain number per material point
   do micro = 1_pInt,material_Nmicrostructure                                                       ! all pairs of homog and micro
     if (Ngrains(homog,micro) > 0_pInt) then                                                        ! an active pair of homog and micro
       myNgrains = Ngrains(homog,micro)                                                             ! assign short name for total number of grains to populate
       myNconstituents = microstructure_Nconstituents(micro)                                        ! assign short name for number of constituents
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
         !$OMP CRITICAL (write2out)
           write(6,'(/,a32,1x,a32,1x,i6)') homogenization_name(homog),microstructure_name(micro),myNgrains
         !$OMP END CRITICAL (write2out)
       endif


!--------------------------------------------------------------------------------------------------
! calculate volume of each grain

       volumeOfGrain = 0.0_pReal
       grain = 0_pInt

       do hme = 1_pInt, Nelems(homog,micro)
         e = elemsOfHomogMicro(homog,micro)%p(hme)                                                  ! my combination of homog and micro, only perform calculations for elements with homog, micro combinations which is indexed in cpElemsindex
         t = FE_geomtype(mesh_element(2,e))
         if (microstructure_elemhomo(micro)) then                                                   ! homogeneous distribution of grains over each element's IPs
           volumeOfGrain(grain+1_pInt:grain+dGrains) = sum(mesh_ipVolume(1:FE_Nips(t),e))/&
                                                                         real(dGrains,pReal)        ! each grain combines size of all IPs in that element
           grain = grain + dGrains                                                                  ! wind forward by Ngrains@IP
         else
           forall (i = 1_pInt:FE_Nips(t)) &                                                         ! loop over IPs
             volumeOfGrain(grain+(i-1)*dGrains+1_pInt:grain+i*dGrains) = &
               mesh_ipVolume(i,e)/dGrains                                                           ! assign IPvolume/Ngrains@IP to all grains of IP
           grain = grain + FE_Nips(t) * dGrains                                                     ! wind forward by Nips*Ngrains@IP
         endif
       enddo

       if (grain /= myNgrains) &
         call IO_error(0,el = homog,ip = micro,ext_msg = 'inconsistent grain count after volume calc')

!--------------------------------------------------------------------------------------------------
! divide myNgrains as best over constituents
!
! example: three constituents with fractions of 0.25, 0.25, and 0.5 distributed over 20 (microstructure) grains
!
!                       ***** ***** **********
! NgrainsOfConstituent: 5,    5,    10
! counters:
!                      |-----> grain (if constituent == 2)
!                            |--> constituentGrain (of constituent 2)
!

       NgrainsOfConstituent = 0_pInt                                                                ! reset counter of grains per constituent
       forall (i = 1_pInt:myNconstituents) &
         NgrainsOfConstituent(i) = nint(microstructure_fraction(i,micro) * myNgrains, pInt)         ! do rounding integer conversion
       do while (sum(NgrainsOfConstituent) /= myNgrains)                                            ! total grain count over constituents wrong?
         sgn = sign(1_pInt, myNgrains - sum(NgrainsOfConstituent))                                  ! direction of required change
         extreme = 0.0_pReal
         t = 0_pInt
         do i = 1_pInt,myNconstituents                                                              ! find largest deviator
           if (real(sgn,pReal)*log(NgrainsOfConstituent(i)/myNgrains/microstructure_fraction(i,micro)) > extreme) then
             extreme = real(sgn,pReal)*log(NgrainsOfConstituent(i)/myNgrains/microstructure_fraction(i,micro))
             t = i
           endif
         enddo
         NgrainsOfConstituent(t) = NgrainsOfConstituent(t) + sgn                                    ! change that by one
       enddo

!--------------------------------------------------------------------------------------------------
! assign phase and texture info

       phaseOfGrain = 0_pInt
       textureOfGrain = 0_pInt
       orientationOfGrain = 0.0_pReal

       texture: do i = 1_pInt,myNconstituents                                                       ! loop over constituents
         grain            = sum(NgrainsOfConstituent(1_pInt:i-1_pInt))                              ! set microstructure grain index of current constituent
                                                                                                    ! "grain" points to start of this constituent's grain population
         constituentGrain = 0_pInt                                                                  ! constituent grain index

         phaseID   = microstructure_phase(i,micro)
         textureID = microstructure_texture(i,micro)
         phaseOfGrain  (grain+1_pInt:grain+NgrainsOfConstituent(i)) = phaseID                       ! assign resp. phase
         textureOfGrain(grain+1_pInt:grain+NgrainsOfConstituent(i)) = textureID                     ! assign resp. texture

         myNorientations = ceiling(real(NgrainsOfConstituent(i),pReal)/&
                                   real(texture_symmetry(textureID),pReal),pInt)                    ! max number of unique orientations (excl. symmetry)

!--------------------------------------------------------------------------------------------------
! ...has texture components
         if (texture_ODFfile(textureID) == '') then
           gauss: do t = 1_pInt,texture_Ngauss(textureID)                                           ! loop over Gauss components
             do g = 1_pInt,int(myNorientations*texture_Gauss(5,t,textureID),pInt)                   ! loop over required grain count
               orientationOfGrain(:,grain+constituentGrain+g) = &
                 math_sampleGaussOri(texture_Gauss(1:3,t,textureID),&
                                     texture_Gauss(  4,t,textureID))
             enddo
             constituentGrain = &
             constituentGrain + int(myNorientations*texture_Gauss(5,t,textureID))                   ! advance counter for grains of current constituent
           enddo gauss

           fiber: do t = 1_pInt,texture_Nfiber(textureID)                                           ! loop over fiber components
             do g = 1_pInt,int(myNorientations*texture_Fiber(6,t,textureID),pInt)                   ! loop over required grain count
               orientationOfGrain(:,grain+constituentGrain+g) = &
                 math_sampleFiberOri(texture_Fiber(1:2,t,textureID),&
                                     texture_Fiber(3:4,t,textureID),&
                                     texture_Fiber(  5,t,textureID))
             enddo
             constituentGrain = &
             constituentGrain + int(myNorientations*texture_fiber(6,t,textureID),pInt)              ! advance counter for grains of current constituent
           enddo fiber

           random: do constituentGrain = constituentGrain+1_pInt,myNorientations                    ! fill remainder with random
              orientationOfGrain(:,grain+constituentGrain) = math_sampleRandomOri()
           enddo random
!--------------------------------------------------------------------------------------------------
! ...has hybrid IA
         else
           orientationOfGrain(1:3,grain+1_pInt:grain+myNorientations) = &
                                            IO_hybridIA(myNorientations,texture_ODFfile(textureID))
           if (all(orientationOfGrain(1:3,grain+1_pInt) == -1.0_pReal)) call IO_error(156_pInt)
         endif

!--------------------------------------------------------------------------------------------------
! ...texture transformation

         do j = 1_pInt,myNorientations                                                              ! loop over each "real" orientation
           orientationOfGrain(1:3,grain+j) = math_RtoEuler( &                                       ! translate back to Euler angles
                                             math_mul33x33( &                                       ! pre-multiply
                                               math_EulertoR(orientationOfGrain(1:3,grain+j)), &    ! face-value orientation
                                               texture_transformation(1:3,1:3,textureID) &          ! and transformation matrix
                                             ) &
                                             )
         enddo

!--------------------------------------------------------------------------------------------------
! ...sample symmetry

         symExtension = texture_symmetry(textureID) - 1_pInt
         if (symExtension > 0_pInt) then                                                            ! sample symmetry (number of additional equivalent orientations)
           constituentGrain = myNorientations                                                       ! start right after "real" orientations
           do j = 1_pInt,myNorientations                                                            ! loop over each "real" orientation
             symOrientation = math_symmetricEulers(texture_symmetry(textureID), &
                                                   orientationOfGrain(1:3,grain+j))                 ! get symmetric equivalents
             e = min(symExtension,NgrainsOfConstituent(i)-constituentGrain)                         ! do not overshoot end of constituent grain array
             if (e > 0_pInt) then
               orientationOfGrain(1:3,grain+constituentGrain+1:   &
                                      grain+constituentGrain+e) = &
                 symOrientation(1:3,1:e)
               constituentGrain = constituentGrain + e                                              ! remainder shrinks by e
             endif
           enddo
         endif

!--------------------------------------------------------------------------------------------------
! shuffle grains within current constituent

         do j = 1_pInt,NgrainsOfConstituent(i)-1_pInt                                               ! walk thru grains of current constituent
           call random_number(rnd)
           t = nint(rnd*(NgrainsOfConstituent(i)-j)+j+0.5_pReal,pInt)                               ! select a grain in remaining list
           m                               = phaseOfGrain(grain+t)                                  ! exchange current with random
           phaseOfGrain(grain+t)           = phaseOfGrain(grain+j)
           phaseOfGrain(grain+j)           = m
           m                               = textureOfGrain(grain+t)                                ! exchange current with random
           textureOfGrain(grain+t)         = textureOfGrain(grain+j)
           textureOfGrain(grain+j)         = m
           orientation                     = orientationOfGrain(1:3,grain+t)                        ! exchange current with random
           orientationOfGrain(1:3,grain+t) = orientationOfGrain(1:3,grain+j)
           orientationOfGrain(1:3,grain+j) = orientation
         enddo

       enddo texture
!< @todo calc fraction after weighing with volumePerGrain, exchange in MC steps to improve result (humbug at the moment)



!--------------------------------------------------------------------------------------------------
! distribute grains of all constituents as accurately as possible to given constituent fractions

       ip = 0_pInt
       currentGrainOfConstituent = 0_pInt

       do hme = 1_pInt, Nelems(homog,micro)
         e = elemsOfHomogMicro(homog,micro)%p(hme)                                                  ! only perform calculations for elements with homog, micro combinations which is indexed in cpElemsindex
         t = FE_geomtype(mesh_element(2,e))
         if (microstructure_elemhomo(micro)) then                                                   ! homogeneous distribution of grains over each element's IPs
           m = 1_pInt                                                                               ! process only first IP
         else
           m = FE_Nips(t)                                                                           ! process all IPs
         endif

         do i = 1_pInt, m                                                                           ! loop over necessary IPs
           ip = ip + 1_pInt                                                                         ! keep track of total ip count
           ipGrain = 0_pInt                                                                         ! count number of grains assigned at this IP
           randomOrder = math_range(microstructure_maxNconstituents)                                ! start out with ordered sequence of constituents
           call random_number(rndArray)                                                             ! as many rnd numbers as (max) constituents
           do j = 1_pInt, myNconstituents - 1_pInt                                                  ! loop over constituents ...
             r = nint(rndArray(j)*(myNconstituents-j)+j+0.5_pReal,pInt)                             ! ... select one in remaining list
             c = randomOrder(r)                                                                     ! ... call it "c"
             randomOrder(r) = randomOrder(j)                                                        ! ... and exchange with present position in constituent list
             grain = sum(NgrainsOfConstituent(1:c-1_pInt))                                          ! figure out actual starting index in overall/consecutive grain population
             do g = 1_pInt, min(dGrains-ipGrain, &                                                  ! leftover number of grains at this IP
                                max(0_pInt, &                                                       ! no negative values
                                    nint(real(ip * dGrains * NgrainsOfConstituent(c)) / &           ! fraction of grains scaled to this constituent...
                                         real(myNgrains),pInt) - &                                  ! ...minus those already distributed
                                         currentGrainOfConstituent(c)))
               ipGrain = ipGrain + 1_pInt                                                           ! advance IP grain counter
               currentGrainOfConstituent(c)  = currentGrainOfConstituent(c) + 1_pInt                ! advance index of grain population for constituent c
               material_volume(ipGrain,i,e)  = volumeOfGrain(grain+currentGrainOfConstituent(c))    ! assign properties
               material_phase(ipGrain,i,e)   = phaseOfGrain(grain+currentGrainOfConstituent(c))
               material_texture(ipGrain,i,e) = textureOfGrain(grain+currentGrainOfConstituent(c))
               material_EulerAngles(1:3,ipGrain,i,e) = orientationOfGrain(1:3,grain+currentGrainOfConstituent(c))
           enddo; enddo

           c = randomOrder(microstructure_Nconstituents(micro))                                     ! look up constituent remaining after random shuffling
           grain = sum(NgrainsOfConstituent(1:c-1_pInt))                                            ! figure out actual starting index in overall/consecutive grain population
           do ipGrain = ipGrain + 1_pInt, dGrains                                                   ! ensure last constituent fills up to dGrains
             currentGrainOfConstituent(c)  = currentGrainOfConstituent(c) + 1_pInt
             material_volume(ipGrain,i,e)  = volumeOfGrain(grain+currentGrainOfConstituent(c))
             material_phase(ipGrain,i,e)   = phaseOfGrain(grain+currentGrainOfConstituent(c))
             material_texture(ipGrain,i,e) = textureOfGrain(grain+currentGrainOfConstituent(c))
             material_EulerAngles(1:3,ipGrain,i,e) = orientationOfGrain(1:3,grain+currentGrainOfConstituent(c))
           enddo

         enddo

         do i = i, FE_Nips(t)                                                                       ! loop over IPs to (possibly) distribute copies from first IP
           material_volume (1_pInt:dGrains,i,e) = material_volume (1_pInt:dGrains,1,e)
           material_phase  (1_pInt:dGrains,i,e) = material_phase  (1_pInt:dGrains,1,e)
           material_texture(1_pInt:dGrains,i,e) = material_texture(1_pInt:dGrains,1,e)
           material_EulerAngles(1:3,1_pInt:dGrains,i,e) = material_EulerAngles(1:3,1_pInt:dGrains,1,e)
         enddo

       enddo
     endif                                                                                          ! active homog,micro pair
   enddo
 enddo

 deallocate(volumeOfGrain)
 deallocate(phaseOfGrain)
 deallocate(textureOfGrain)
 deallocate(orientationOfGrain)
 deallocate(Nelems)
 !> @todo - causing segmentation fault: needs looking into
 !do homog = 1,material_Nhomogenization
 !  do micro = 1,material_Nmicrostructure
 !    if (Nelems(homog,micro) > 0_pInt) deallocate(elemsOfHomogMicro(homog,micro)%p)
 !  enddo
 !enddo
 deallocate(elemsOfHomogMicro)

end subroutine material_populateGrains

#ifdef HDF
integer(pInt) pure function material_NconstituentsPhase(matID)

 implicit none
 integer(pInt), intent(in) :: matID

 material_NconstituentsPhase = count(microstructure_phase == matID)
end function
#endif

end module material
