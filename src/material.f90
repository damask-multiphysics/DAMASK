!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
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
   group_float, &
   group_int

 implicit none
 private
 character(len=*),                         parameter,            public :: &
   ELASTICITY_hooke_label               = 'hooke', &
   PLASTICITY_none_label                = 'none', &
   PLASTICITY_isotropic_label           = 'isotropic', &
   PLASTICITY_phenopowerlaw_label       = 'phenopowerlaw', &
   PLASTICITY_kinehardening_label       = 'kinehardening', &
   PLASTICITY_dislotwin_label           = 'dislotwin', &
   PLASTICITY_disloucla_label           = 'disloucla', &
   PLASTICITY_nonlocal_label            = 'nonlocal', &
   SOURCE_thermal_dissipation_label     = 'thermal_dissipation', &
   SOURCE_thermal_externalheat_label    = 'thermal_externalheat', &
   SOURCE_damage_isoBrittle_label       = 'damage_isobrittle', &
   SOURCE_damage_isoDuctile_label       = 'damage_isoductile', &
   SOURCE_damage_anisoBrittle_label     = 'damage_anisobrittle', &
   SOURCE_damage_anisoDuctile_label     = 'damage_anisoductile', &
   KINEMATICS_thermal_expansion_label   = 'thermal_expansion', &
   KINEMATICS_cleavage_opening_label    = 'cleavage_opening', &
   KINEMATICS_slipplane_opening_label   = 'slipplane_opening', &
   STIFFNESS_DEGRADATION_damage_label   = 'damage', &
   THERMAL_isothermal_label             = 'isothermal', &
   THERMAL_adiabatic_label              = 'adiabatic', &
   THERMAL_conduction_label             = 'conduction', &
   DAMAGE_none_label                    = 'none', &
   DAMAGE_local_label                   = 'local', &
   DAMAGE_nonlocal_label                = 'nonlocal', &
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
                 PLASTICITY_phenopowerlaw_ID, &
                 PLASTICITY_kinehardening_ID, &
                 PLASTICITY_dislotwin_ID, &
                 PLASTICITY_disloucla_ID, &
                 PLASTICITY_nonlocal_ID
 end enum

 enum, bind(c)
   enumerator :: SOURCE_undefined_ID, &
                 SOURCE_thermal_dissipation_ID, &
                 SOURCE_thermal_externalheat_ID, &
                 SOURCE_damage_isoBrittle_ID, &
                 SOURCE_damage_isoDuctile_ID, &
                 SOURCE_damage_anisoBrittle_ID, &
                 SOURCE_damage_anisoDuctile_ID
 end enum

 enum, bind(c)
   enumerator :: KINEMATICS_undefined_ID, &
                 KINEMATICS_cleavage_opening_ID, &
                 KINEMATICS_slipplane_opening_ID, &
                 KINEMATICS_thermal_expansion_ID
 end enum

 enum, bind(c)
   enumerator :: STIFFNESS_DEGRADATION_undefined_ID, &
                 STIFFNESS_DEGRADATION_damage_ID
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
   enumerator :: HOMOGENIZATION_undefined_ID, &
                 HOMOGENIZATION_none_ID, &
                 HOMOGENIZATION_isostrain_ID, &
                 HOMOGENIZATION_rgc_ID
 end enum

 integer(kind(ELASTICITY_undefined_ID)),     dimension(:),   allocatable, public, protected :: &
   phase_elasticity                                                                                 !< elasticity of each phase
 integer(kind(PLASTICITY_undefined_ID)),     dimension(:),   allocatable, public, protected :: &
   phase_plasticity                                                                                 !< plasticity of each phase
 integer(kind(THERMAL_isothermal_ID)),       dimension(:),   allocatable, public, protected :: &
   thermal_type                                                                                     !< thermal transport model
 integer(kind(DAMAGE_none_ID)),              dimension(:),   allocatable, public, protected :: &
   damage_type                                                                                      !< nonlocal damage model

 integer(kind(SOURCE_undefined_ID)),         dimension(:,:), allocatable, public, protected :: &
   phase_source, &                                                                                  !< active sources mechanisms of each phase
   phase_kinematics, &                                                                              !< active kinematic mechanisms of each phase
   phase_stiffnessDegradation                                                                       !< active stiffness degradation mechanisms of each phase

 integer(kind(HOMOGENIZATION_undefined_ID)), dimension(:),   allocatable, public, protected :: &
   homogenization_type                                                                              !< type of each homogenization

 integer(pInt), public, protected :: &
   homogenization_maxNgrains                                                                        !< max number of grains in any USED homogenization

 integer(pInt), dimension(:), allocatable, public, protected :: &
   phase_Nsources, &                                                                                !< number of source mechanisms active in each phase
   phase_Nkinematics, &                                                                             !< number of kinematic mechanisms active in each phase
   phase_NstiffnessDegradations, &                                                                  !< number of stiffness degradation mechanisms active in each phase
   phase_Noutput, &                                                                                 !< number of '(output)' items per phase
   phase_elasticityInstance, &                                                                      !< instance of particular elasticity of each phase
   phase_plasticityInstance, &                                                                      !< instance of particular plasticity of each phase
   crystallite_Noutput, &                                                                           !< number of '(output)' items per crystallite setting
   homogenization_Ngrains, &                                                                        !< number of grains in each homogenization
   homogenization_Noutput, &                                                                        !< number of '(output)' items per homogenization
   homogenization_typeInstance, &                                                                   !< instance of particular type of each homogenization
   thermal_typeInstance, &                                                                          !< instance of particular type of each thermal transport
   damage_typeInstance, &                                                                           !< instance of particular type of each nonlocal damage
   microstructure_crystallite                                                                       !< crystallite setting ID of each microstructure ! DEPRECATED !!!!

 real(pReal), dimension(:), allocatable, public, protected :: &
   thermal_initialT, &                                                                              !< initial temperature per each homogenization
   damage_initialPhi                                                                                !< initial damage per each homogenization

! NEW MAPPINGS 
 integer, dimension(:),     allocatable, public, protected :: &                                     ! (elem)
   material_homogenizationAt                                                                        !< homogenization ID of each element (copy of mesh_homogenizationAt)
 integer, dimension(:,:),   allocatable, public, protected :: &                                     ! (ip,elem)
   material_homogenizationMemberAt                                                                  !< position of the element within its homogenization instance
 integer, dimension(:,:), allocatable, public, protected :: &                                       ! (constituent,elem)
   material_phaseAt                                                                                 !< phase ID of each element
 integer, dimension(:,:,:), allocatable, public, protected :: &                                     ! (constituent,ip,elem)
   material_phaseMemberAt                                                                           !< position of the element within its phase instance
! END NEW MAPPINGS
 
! DEPRECATED: use material_phaseAt
 integer(pInt), dimension(:,:,:), allocatable, public :: &
   material_phase                                                                                   !< phase (index) of each grain,IP,element

 type(tPlasticState), allocatable, dimension(:), public :: &
   plasticState
 type(tSourceState),  allocatable, dimension(:), public :: &
   sourceState
 type(tState),        allocatable, dimension(:), public :: &
   homogState, &
   thermalState, &
   damageState

 integer(pInt), dimension(:,:,:), allocatable, public, protected :: &
   material_texture                                                                                 !< texture (index) of each grain,IP,element

 real(pReal), dimension(:,:,:,:), allocatable, public, protected :: &
   material_EulerAngles                                                                             !< initial orientation of each grain,IP,element

 logical, dimension(:), allocatable, public, protected :: &
   microstructure_active, &
   microstructure_elemhomo, &                                                                       !< flag to indicate homogeneous microstructure distribution over element's IPs
   phase_localPlasticity                                                                            !< flags phases with local constitutive law

 integer(pInt), private :: &
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

! BEGIN DEPRECATED
 integer(pInt), dimension(:,:,:), allocatable, public :: phaseAt                                    !< phase ID of every material point (ipc,ip,el)
 integer(pInt), dimension(:,:,:), allocatable, public :: phasememberAt                              !< memberID of given phase at every material point (ipc,ip,el)

 integer(pInt), dimension(:,:,:), allocatable, public, target :: mappingHomogenization              !< mapping from material points to offset in heterogenous state/field
 integer(pInt), dimension(:,:),   allocatable, private, target :: mappingHomogenizationConst         !< mapping from material points to offset in constant state/field
! END DEPRECATED

 type(tHomogMapping), allocatable, dimension(:), public :: &
   thermalMapping, &                                                                                !< mapping for thermal state/fields
   damageMapping                                                                                    !< mapping for damage state/fields

 type(group_float),  allocatable, dimension(:), public :: &
   temperature, &                                                                                   !< temperature field
   damage, &                                                                                        !< damage field
   temperatureRate                                                                                  !< temperature change rate field

 public :: &
   material_init, &
   material_allocatePlasticState, &
   material_allocateSourceState, &
   ELASTICITY_hooke_ID ,&
   PLASTICITY_none_ID, &
   PLASTICITY_isotropic_ID, &
   PLASTICITY_phenopowerlaw_ID, &
   PLASTICITY_kinehardening_ID, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_nonlocal_ID, &
   SOURCE_thermal_dissipation_ID, &
   SOURCE_thermal_externalheat_ID, &
   SOURCE_damage_isoBrittle_ID, &
   SOURCE_damage_isoDuctile_ID, &
   SOURCE_damage_anisoBrittle_ID, &
   SOURCE_damage_anisoDuctile_ID, &
   KINEMATICS_cleavage_opening_ID, &
   KINEMATICS_slipplane_opening_ID, &
   KINEMATICS_thermal_expansion_ID, &
   STIFFNESS_DEGRADATION_damage_ID, &
   THERMAL_isothermal_ID, &
   THERMAL_adiabatic_ID, &
   THERMAL_conduction_ID, &
   DAMAGE_none_ID, &
   DAMAGE_local_ID, &
   DAMAGE_nonlocal_ID, &
   HOMOGENIZATION_none_ID, &
   HOMOGENIZATION_isostrain_ID, &
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
subroutine material_init
#if defined(PETSc) || defined(DAMASK_HDF5)
 use results
#endif
 use IO, only: &
   IO_error
 use debug, only: &
   debug_level, &
   debug_material, &
   debug_levelBasic, &
   debug_levelExtensive
 use config, only: &
   config_crystallite, &
   config_homogenization, &
   config_microstructure, &
   config_phase, &
   config_texture, &
   homogenization_name, &
   microstructure_name, &
   phase_name, &
   texture_name
 use mesh, only: &
   theMesh

 implicit none
 integer(pInt), parameter :: FILEUNIT = 210_pInt
 integer(pInt)            :: m,c,h, myDebug, myPhase, myHomog
 integer(pInt) :: &
  g, &                                                                                              !< grain number
  i, &                                                                                              !< integration point number
  e                                                                                                 !< element number
 integer(pInt), dimension(:), allocatable :: &
  CounterPhase, &
  CounterHomogenization

 myDebug = debug_level(debug_material)

 write(6,'(/,a)') ' <<<+-  material init  -+>>>'

 call material_parsePhase()
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Phase          parsed'; flush(6)
 
 call material_parseMicrostructure()
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Microstructure parsed'; flush(6)
 
 call material_parseCrystallite()
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Crystallite    parsed'; flush(6)
 
 call material_parseHomogenization()
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Homogenization parsed'; flush(6)
 
 call material_parseTexture()
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Texture        parsed'; flush(6)

 allocate(plasticState       (size(config_phase)))
 allocate(sourceState        (size(config_phase)))
 do myPhase = 1,size(config_phase)
   allocate(sourceState(myPhase)%p(phase_Nsources(myPhase)))
 enddo

 allocate(homogState         (size(config_homogenization)))
 allocate(thermalState       (size(config_homogenization)))
 allocate(damageState        (size(config_homogenization)))

 allocate(thermalMapping     (size(config_homogenization)))
 allocate(damageMapping      (size(config_homogenization)))

 allocate(temperature        (size(config_homogenization)))
 allocate(damage             (size(config_homogenization)))

 allocate(temperatureRate    (size(config_homogenization)))

 do m = 1_pInt,size(config_microstructure)
   if(microstructure_crystallite(m) < 1_pInt .or. &
      microstructure_crystallite(m) > size(config_crystallite)) &
        call IO_error(150_pInt,m,ext_msg='crystallite')
   if(minval(microstructure_phase(1:microstructure_Nconstituents(m),m)) < 1_pInt .or. &
      maxval(microstructure_phase(1:microstructure_Nconstituents(m),m)) > size(config_phase)) &
        call IO_error(150_pInt,m,ext_msg='phase')
   if(minval(microstructure_texture(1:microstructure_Nconstituents(m),m)) < 1_pInt .or. &
      maxval(microstructure_texture(1:microstructure_Nconstituents(m),m)) > size(config_texture)) &
        call IO_error(150_pInt,m,ext_msg='texture')
   if(microstructure_Nconstituents(m) < 1_pInt) &
        call IO_error(151_pInt,m)
 enddo

 debugOut: if (iand(myDebug,debug_levelExtensive) /= 0_pInt) then
   write(6,'(/,a,/)') ' MATERIAL configuration'
   write(6,'(a32,1x,a16,1x,a6)') 'homogenization                  ','type            ','grains'
   do h = 1_pInt,size(config_homogenization)
     write(6,'(1x,a32,1x,a16,1x,i6)') homogenization_name(h),homogenization_type(h),homogenization_Ngrains(h)
   enddo
   write(6,'(/,a14,18x,1x,a11,1x,a12,1x,a13)') 'microstructure','crystallite','constituents','homogeneous'
   do m = 1_pInt,size(config_microstructure)
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
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! new mappings
 allocate(material_homogenizationAt,source=theMesh%homogenizationAt)
 allocate(material_homogenizationMemberAt(theMesh%elem%nIPs,theMesh%Nelems),source=0)

 allocate(CounterHomogenization(size(config_homogenization)),source=0)
 do e = 1, theMesh%Nelems
   do i = 1, theMesh%elem%nIPs
     CounterHomogenization(material_homogenizationAt(e)) = &
     CounterHomogenization(material_homogenizationAt(e)) + 1
     material_homogenizationMemberAt(i,e) = CounterHomogenization(material_homogenizationAt(e))
   enddo
 enddo


 allocate(material_phaseAt(homogenization_maxNgrains,theMesh%Nelems), source=material_phase(:,1,:))
 allocate(material_phaseMemberAt(homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems),source=0)
 
 allocate(CounterPhase(size(config_phase)),source=0)
 do e = 1, theMesh%Nelems
   do i = 1, theMesh%elem%nIPs
     do c = 1, homogenization_maxNgrains
       CounterPhase(material_phaseAt(c,e)) = &
       CounterPhase(material_phaseAt(c,e)) + 1
       material_phaseMemberAt(c,i,e) = CounterPhase(material_phaseAt(c,e))
     enddo
   enddo
 enddo
 
#if defined(PETSc) || defined(DAMASK_HDF5)
 call results_openJobFile
 call results_mapping_constituent(material_phaseAt,material_phaseMemberAt,phase_name)
 call results_mapping_materialpoint(material_homogenizationAt,material_homogenizationMemberAt,homogenization_name)
 call results_closeJobFile
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BEGIN DEPRECATED
 allocate(phaseAt                   (  homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems),source=0_pInt)
 allocate(phasememberAt             (  homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems),source=0_pInt)
 allocate(mappingHomogenization     (2,                          theMesh%elem%nIPs,theMesh%Nelems),source=0_pInt)
 allocate(mappingHomogenizationConst(                            theMesh%elem%nIPs,theMesh%Nelems),source=1_pInt)
 
 CounterHomogenization=0
 CounterPhase         =0


 do e = 1_pInt,theMesh%Nelems
 myHomog = theMesh%homogenizationAt(e)
   do i = 1_pInt, theMesh%elem%nIPs
     CounterHomogenization(myHomog) = CounterHomogenization(myHomog) + 1_pInt
     mappingHomogenization(1:2,i,e) = [CounterHomogenization(myHomog),huge(1)]
     do g = 1_pInt,homogenization_Ngrains(myHomog)
       myPhase = material_phase(g,i,e)
       CounterPhase(myPhase) = CounterPhase(myPhase)+1_pInt                             ! not distinguishing between instances of same phase
       phaseAt(g,i,e)              = myPhase
       phasememberAt(g,i,e)        = CounterPhase(myPhase)
     enddo
   enddo
 enddo
! END DEPRECATED

! REMOVE !!!!!
! hack needed to initialize field values used during constitutive and crystallite initializations
 do myHomog = 1,size(config_homogenization)
   thermalMapping     (myHomog)%p => mappingHomogenizationConst
   damageMapping      (myHomog)%p => mappingHomogenizationConst
   allocate(temperature     (myHomog)%p(1), source=thermal_initialT(myHomog))
   allocate(damage          (myHomog)%p(1), source=damage_initialPhi(myHomog))
   allocate(temperatureRate (myHomog)%p(1), source=0.0_pReal)
 enddo

end subroutine material_init


!--------------------------------------------------------------------------------------------------
!> @brief parses the homogenization part from the material configuration
!--------------------------------------------------------------------------------------------------
subroutine material_parseHomogenization
 use config, only : &
   config_homogenization
 use mesh, only: &
   theMesh
 use IO, only: &
   IO_error

 implicit none
 integer(pInt)        :: h
 character(len=65536) :: tag

 allocate(homogenization_type(size(config_homogenization)),           source=HOMOGENIZATION_undefined_ID)
 allocate(thermal_type(size(config_homogenization)),                  source=THERMAL_isothermal_ID)
 allocate(damage_type (size(config_homogenization)),                  source=DAMAGE_none_ID)
 allocate(homogenization_typeInstance(size(config_homogenization)),   source=0_pInt)
 allocate(thermal_typeInstance(size(config_homogenization)),          source=0_pInt)
 allocate(damage_typeInstance(size(config_homogenization)),           source=0_pInt)
 allocate(homogenization_Ngrains(size(config_homogenization)),        source=0_pInt)
 allocate(homogenization_Noutput(size(config_homogenization)),        source=0_pInt)
 allocate(homogenization_active(size(config_homogenization)),         source=.false.)  !!!!!!!!!!!!!!!
 allocate(thermal_initialT(size(config_homogenization)),              source=300.0_pReal)
 allocate(damage_initialPhi(size(config_homogenization)),             source=1.0_pReal)

 forall (h = 1_pInt:size(config_homogenization)) &
   homogenization_active(h) = any(theMesh%homogenizationAt == h)


 do h=1_pInt, size(config_homogenization)
   homogenization_Noutput(h) = config_homogenization(h)%countKeys('(output)')

   tag = config_homogenization(h)%getString('mech')
   select case (trim(tag))
     case(HOMOGENIZATION_NONE_label)
       homogenization_type(h) = HOMOGENIZATION_NONE_ID
       homogenization_Ngrains(h) = 1_pInt
     case(HOMOGENIZATION_ISOSTRAIN_label)
       homogenization_type(h) = HOMOGENIZATION_ISOSTRAIN_ID
       homogenization_Ngrains(h) = config_homogenization(h)%getInt('nconstituents')
     case(HOMOGENIZATION_RGC_label)
       homogenization_type(h) = HOMOGENIZATION_RGC_ID
       homogenization_Ngrains(h) = config_homogenization(h)%getInt('nconstituents')
     case default
       call IO_error(500_pInt,ext_msg=trim(tag))
   end select
   
   homogenization_typeInstance(h) = count(homogenization_type==homogenization_type(h))

   if (config_homogenization(h)%keyExists('thermal')) then
     thermal_initialT(h) =  config_homogenization(h)%getFloat('t0',defaultVal=300.0_pReal)

     tag = config_homogenization(h)%getString('thermal')
     select case (trim(tag))
       case(THERMAL_isothermal_label)
         thermal_type(h) = THERMAL_isothermal_ID
       case(THERMAL_adiabatic_label)
         thermal_type(h) = THERMAL_adiabatic_ID
       case(THERMAL_conduction_label)
         thermal_type(h) = THERMAL_conduction_ID
       case default
         call IO_error(500_pInt,ext_msg=trim(tag))
     end select

   endif

   if (config_homogenization(h)%keyExists('damage')) then
     damage_initialPhi(h) =  config_homogenization(h)%getFloat('initialdamage',defaultVal=1.0_pReal)

     tag = config_homogenization(h)%getString('damage')
     select case (trim(tag))
       case(DAMAGE_NONE_label)
         damage_type(h) = DAMAGE_none_ID
       case(DAMAGE_LOCAL_label)
         damage_type(h) = DAMAGE_local_ID
       case(DAMAGE_NONLOCAL_label)
         damage_type(h) = DAMAGE_nonlocal_ID
       case default
         call IO_error(500_pInt,ext_msg=trim(tag))
     end select

   endif

 enddo

 do h=1_pInt, size(config_homogenization)
   homogenization_typeInstance(h)  = count(homogenization_type(1:h)  == homogenization_type(h))
   thermal_typeInstance(h)         = count(thermal_type       (1:h)  == thermal_type       (h))
   damage_typeInstance(h)          = count(damage_type        (1:h)  == damage_type        (h))
 enddo

 homogenization_maxNgrains = maxval(homogenization_Ngrains,homogenization_active)

end subroutine material_parseHomogenization


!--------------------------------------------------------------------------------------------------
!> @brief parses the microstructure part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseMicrostructure
 use prec, only: &
  dNeq
 use IO, only: &
   IO_floatValue, &
   IO_intValue, &
   IO_stringValue, &
   IO_stringPos, &
   IO_error
 use config, only: &
   config_microstructure, &
   microstructure_name
 use mesh, only: &
   theMesh

 implicit none
 character(len=65536), dimension(:), allocatable :: &
   strings
 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: e, m, c, i
 character(len=65536) :: &
   tag

 allocate(microstructure_crystallite(size(config_microstructure)),          source=0_pInt)
 allocate(microstructure_Nconstituents(size(config_microstructure)),        source=0_pInt)
 allocate(microstructure_active(size(config_microstructure)),               source=.false.)
 allocate(microstructure_elemhomo(size(config_microstructure)),             source=.false.)

 if(any(theMesh%microstructureAt > size(config_microstructure))) &
  call IO_error(155_pInt,ext_msg='More microstructures in geometry than sections in material.config')

 forall (e = 1_pInt:theMesh%Nelems) &
   microstructure_active(theMesh%microstructureAt(e)) = .true.                                         ! current microstructure used in model? Elementwise view, maximum N operations for N elements

 do m=1_pInt, size(config_microstructure)
   microstructure_Nconstituents(m) =  config_microstructure(m)%countKeys('(constituent)')
   microstructure_crystallite(m)   =  config_microstructure(m)%getInt('crystallite')
   microstructure_elemhomo(m)      =  config_microstructure(m)%keyExists('/elementhomogeneous/')
 enddo

 microstructure_maxNconstituents = maxval(microstructure_Nconstituents)
 allocate(microstructure_phase   (microstructure_maxNconstituents,size(config_microstructure)),source=0_pInt)
 allocate(microstructure_texture (microstructure_maxNconstituents,size(config_microstructure)),source=0_pInt)
 allocate(microstructure_fraction(microstructure_maxNconstituents,size(config_microstructure)),source=0.0_pReal)

 allocate(strings(1))                                                                               ! Intel 16.0 Bug
 do m=1_pInt, size(config_microstructure)
   strings = config_microstructure(m)%getStrings('(constituent)',raw=.true.)
   do c = 1_pInt, size(strings)
     chunkPos = IO_stringPos(strings(c))

     do i = 1_pInt,5_pInt,2_pInt
        tag = IO_stringValue(strings(c),chunkPos,i)

        select case (tag)
          case('phase')
            microstructure_phase(c,m) =    IO_intValue(strings(c),chunkPos,i+1_pInt)
          case('texture')
            microstructure_texture(c,m) =  IO_intValue(strings(c),chunkPos,i+1_pInt)
          case('fraction')
            microstructure_fraction(c,m) =  IO_floatValue(strings(c),chunkPos,i+1_pInt)
        end select
     
     enddo
   enddo
 enddo

 do m = 1_pInt, size(config_microstructure)
   if (dNeq(sum(microstructure_fraction(:,m)),1.0_pReal)) &
     call IO_error(153_pInt,ext_msg=microstructure_name(m))
 enddo
 
end subroutine material_parseMicrostructure


!--------------------------------------------------------------------------------------------------
!> @brief parses the crystallite part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseCrystallite
 use config, only: &
   config_crystallite

 implicit none
 integer(pInt)        :: c

 allocate(crystallite_Noutput(size(config_crystallite)),source=0_pInt)
 do c=1_pInt, size(config_crystallite)
   crystallite_Noutput(c) =  config_crystallite(c)%countKeys('(output)')
 enddo

end subroutine material_parseCrystallite


!--------------------------------------------------------------------------------------------------
!> @brief parses the phase part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parsePhase
 use IO, only: &
   IO_error, &
   IO_getTag, &
   IO_stringValue
 use config, only: &
   config_phase

 implicit none
 integer(pInt) :: sourceCtr, kinematicsCtr, stiffDegradationCtr, p
 character(len=65536), dimension(:), allocatable ::  str 


 allocate(phase_elasticity(size(config_phase)),source=ELASTICITY_undefined_ID)
 allocate(phase_plasticity(size(config_phase)),source=PLASTICITY_undefined_ID)
 allocate(phase_Nsources(size(config_phase)),              source=0_pInt)
 allocate(phase_Nkinematics(size(config_phase)),           source=0_pInt)
 allocate(phase_NstiffnessDegradations(size(config_phase)),source=0_pInt)
 allocate(phase_Noutput(size(config_phase)),               source=0_pInt)
 allocate(phase_localPlasticity(size(config_phase)),       source=.false.)

 do p=1_pInt, size(config_phase)
   phase_Noutput(p) =                 config_phase(p)%countKeys('(output)')
   phase_Nsources(p) =                config_phase(p)%countKeys('(source)')
   phase_Nkinematics(p) =             config_phase(p)%countKeys('(kinematics)')
   phase_NstiffnessDegradations(p) =  config_phase(p)%countKeys('(stiffness_degradation)')
   phase_localPlasticity(p) = .not.   config_phase(p)%KeyExists('/nonlocal/')

   select case (config_phase(p)%getString('elasticity'))
     case (ELASTICITY_HOOKE_label)
       phase_elasticity(p) = ELASTICITY_HOOKE_ID
     case default
       call IO_error(200_pInt,ext_msg=trim(config_phase(p)%getString('elasticity')))
   end select

   select case (config_phase(p)%getString('plasticity'))
     case (PLASTICITY_NONE_label)
       phase_plasticity(p) = PLASTICITY_NONE_ID
     case (PLASTICITY_ISOTROPIC_label)
       phase_plasticity(p) = PLASTICITY_ISOTROPIC_ID
     case (PLASTICITY_PHENOPOWERLAW_label)
       phase_plasticity(p) = PLASTICITY_PHENOPOWERLAW_ID
     case (PLASTICITY_KINEHARDENING_label)
       phase_plasticity(p) = PLASTICITY_KINEHARDENING_ID
     case (PLASTICITY_DISLOTWIN_label)
       phase_plasticity(p) = PLASTICITY_DISLOTWIN_ID
     case (PLASTICITY_DISLOUCLA_label)
       phase_plasticity(p) = PLASTICITY_DISLOUCLA_ID
     case (PLASTICITY_NONLOCAL_label)
       phase_plasticity(p) = PLASTICITY_NONLOCAL_ID
     case default
       call IO_error(201_pInt,ext_msg=trim(config_phase(p)%getString('plasticity')))
   end select

 enddo

 allocate(phase_source(maxval(phase_Nsources),size(config_phase)), source=SOURCE_undefined_ID)
 allocate(phase_kinematics(maxval(phase_Nkinematics),size(config_phase)), source=KINEMATICS_undefined_ID)
 allocate(phase_stiffnessDegradation(maxval(phase_NstiffnessDegradations),size(config_phase)), &
          source=STIFFNESS_DEGRADATION_undefined_ID)
 do p=1_pInt, size(config_phase)
#if defined(__GFORTRAN__) || defined(__PGI)
   str = ['GfortranBug86277']
   str = config_phase(p)%getStrings('(source)',defaultVal=str)
   if (str(1) == 'GfortranBug86277') str = [character(len=65536)::]
#else
   str = config_phase(p)%getStrings('(source)',defaultVal=[character(len=65536)::])
#endif
   do sourceCtr = 1_pInt, size(str)
     select case (trim(str(sourceCtr)))
       case (SOURCE_thermal_dissipation_label)
         phase_source(sourceCtr,p) = SOURCE_thermal_dissipation_ID
       case (SOURCE_thermal_externalheat_label)
         phase_source(sourceCtr,p) = SOURCE_thermal_externalheat_ID
       case (SOURCE_damage_isoBrittle_label)
         phase_source(sourceCtr,p) = SOURCE_damage_isoBrittle_ID
       case (SOURCE_damage_isoDuctile_label)
         phase_source(sourceCtr,p) = SOURCE_damage_isoDuctile_ID
       case (SOURCE_damage_anisoBrittle_label)
         phase_source(sourceCtr,p) = SOURCE_damage_anisoBrittle_ID
       case (SOURCE_damage_anisoDuctile_label)
         phase_source(sourceCtr,p) = SOURCE_damage_anisoDuctile_ID
     end select
   enddo

#if defined(__GFORTRAN__) || defined(__PGI)
   str = ['GfortranBug86277']
   str = config_phase(p)%getStrings('(kinematics)',defaultVal=str)
   if (str(1) == 'GfortranBug86277') str = [character(len=65536)::]
#else
   str = config_phase(p)%getStrings('(kinematics)',defaultVal=[character(len=65536)::])
#endif
   do kinematicsCtr = 1_pInt, size(str)
     select case (trim(str(kinematicsCtr)))
       case (KINEMATICS_cleavage_opening_label)
         phase_kinematics(kinematicsCtr,p) = KINEMATICS_cleavage_opening_ID
       case (KINEMATICS_slipplane_opening_label)
         phase_kinematics(kinematicsCtr,p) = KINEMATICS_slipplane_opening_ID
       case (KINEMATICS_thermal_expansion_label)
         phase_kinematics(kinematicsCtr,p) = KINEMATICS_thermal_expansion_ID
     end select
   enddo
#if defined(__GFORTRAN__) || defined(__PGI)
   str = ['GfortranBug86277']
   str = config_phase(p)%getStrings('(stiffness_degradation)',defaultVal=str)
   if (str(1) == 'GfortranBug86277') str = [character(len=65536)::]
#else
   str = config_phase(p)%getStrings('(stiffness_degradation)',defaultVal=[character(len=65536)::])
#endif
   do stiffDegradationCtr = 1_pInt, size(str)
     select case (trim(str(stiffDegradationCtr)))
       case (STIFFNESS_DEGRADATION_damage_label)
         phase_stiffnessDegradation(stiffDegradationCtr,p) = STIFFNESS_DEGRADATION_damage_ID
    end select
   enddo
 enddo

 allocate(phase_plasticityInstance(size(config_phase)),   source=0_pInt)
 allocate(phase_elasticityInstance(size(config_phase)),   source=0_pInt)

 do p=1_pInt, size(config_phase)
   phase_elasticityInstance(p)  = count(phase_elasticity(1:p)  == phase_elasticity(p))
   phase_plasticityInstance(p)  = count(phase_plasticity(1:p)  == phase_plasticity(p))
 enddo

end subroutine material_parsePhase

!--------------------------------------------------------------------------------------------------
!> @brief parses the texture part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseTexture
 use prec, only: &
   dNeq
 use IO, only: &
   IO_error, &
   IO_stringPos, &
   IO_floatValue, &
   IO_stringValue
 use config, only: &
   config_deallocate, &
   config_texture
 use math, only: &
   inRad, &
   math_sampleRandomOri, &
   math_I3, &
   math_det33

 implicit none
 integer(pInt) :: section, gauss, fiber, j, t, i
 character(len=65536), dimension(:), allocatable ::  strings                                     ! Values for given key in material config 
 integer(pInt), dimension(:), allocatable :: chunkPos

 allocate(texture_symmetry(size(config_texture)), source=1_pInt)
 allocate(texture_Ngauss(size(config_texture)),   source=0_pInt)
 allocate(texture_Nfiber(size(config_texture)),   source=0_pInt)

 do t=1_pInt, size(config_texture)
   texture_Ngauss(t) =  config_texture(t)%countKeys('(gauss)') &
                     +  config_texture(t)%countKeys('(random)')
   texture_Nfiber(t) =  config_texture(t)%countKeys('(fiber)')
 enddo

 texture_maxNgauss = maxval(texture_Ngauss)
 texture_maxNfiber = maxval(texture_Nfiber)
 allocate(texture_Gauss (5,texture_maxNgauss,size(config_texture)), source=0.0_pReal)
 allocate(texture_Fiber (6,texture_maxNfiber,size(config_texture)), source=0.0_pReal)
 allocate(texture_transformation(3,3,size(config_texture)),         source=0.0_pReal)
          texture_transformation = spread(math_I3,3,size(config_texture))

 do t=1_pInt, size(config_texture)
   section = t
   gauss = 0_pInt
   fiber = 0_pInt
   
   if (config_texture(t)%keyExists('axes')) then
     strings = config_texture(t)%getStrings('axes')
     do j = 1_pInt, 3_pInt                                                                          ! look for "x", "y", and "z" entries
       select case (strings(j))
         case('x', '+x')
           texture_transformation(j,1:3,t) = [ 1.0_pReal, 0.0_pReal, 0.0_pReal]                     ! original axis is now +x-axis
         case('-x')
           texture_transformation(j,1:3,t) = [-1.0_pReal, 0.0_pReal, 0.0_pReal]                     ! original axis is now -x-axis
         case('y', '+y')
           texture_transformation(j,1:3,t) = [ 0.0_pReal, 1.0_pReal, 0.0_pReal]                     ! original axis is now +y-axis
         case('-y')
           texture_transformation(j,1:3,t) = [ 0.0_pReal,-1.0_pReal, 0.0_pReal]                     ! original axis is now -y-axis
         case('z', '+z')
           texture_transformation(j,1:3,t) = [ 0.0_pReal, 0.0_pReal, 1.0_pReal]                     ! original axis is now +z-axis
         case('-z')
           texture_transformation(j,1:3,t) = [ 0.0_pReal, 0.0_pReal,-1.0_pReal]                     ! original axis is now -z-axis
         case default
           call IO_error(157_pInt,t)
       end select
     enddo
     if(dNeq(math_det33(texture_transformation(1:3,1:3,t)),1.0_pReal)) call IO_error(157_pInt,t)
   endif

   if (config_texture(t)%keyExists('symmetry')) call IO_error(147,ext_msg='symmetry')
   if (config_texture(t)%keyExists('(random)')) call IO_error(147,ext_msg='(random)')
   if (config_texture(t)%keyExists('(fiber)'))  call IO_error(147,ext_msg='(fiber)')
   
   if (config_texture(t)%keyExists('(gauss)')) then
     gauss = gauss + 1_pInt
     strings = config_texture(t)%getStrings('(gauss)',raw= .true.)
     do i = 1_pInt , size(strings)
       chunkPos = IO_stringPos(strings(i))
       do j = 1_pInt,9_pInt,2_pInt
         select case (IO_stringValue(strings(i),chunkPos,j))
             case('phi1')
                 texture_Gauss(1,gauss,t) = IO_floatValue(strings(i),chunkPos,j+1_pInt)*inRad
             case('phi')
                 texture_Gauss(2,gauss,t) = IO_floatValue(strings(i),chunkPos,j+1_pInt)*inRad
             case('phi2')
                 texture_Gauss(3,gauss,t) = IO_floatValue(strings(i),chunkPos,j+1_pInt)*inRad
             case('scatter')
                 texture_Gauss(4,gauss,t) = IO_floatValue(strings(i),chunkPos,j+1_pInt)*inRad
             case('fraction')
                 texture_Gauss(5,gauss,t) = IO_floatValue(strings(i),chunkPos,j+1_pInt)
          end select
      enddo
     enddo
   endif
 enddo    
 
 call config_deallocate('material.config/texture')

end subroutine material_parseTexture


!--------------------------------------------------------------------------------------------------
!> @brief allocates the plastic state of a phase
!--------------------------------------------------------------------------------------------------
subroutine material_allocatePlasticState(phase,NofMyPhase,&
                                         sizeState,sizeDotState,sizeDeltaState,&
                                         Nslip,Ntwin,Ntrans)
 use numerics, only: &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: &
   phase, &
   NofMyPhase, &
   sizeState, &
   sizeDotState, &
   sizeDeltaState, &
   Nslip, &
   Ntwin, &
   Ntrans

 plasticState(phase)%sizeState        = sizeState
 plasticState(phase)%sizeDotState     = sizeDotState
 plasticState(phase)%sizeDeltaState   = sizeDeltaState
 plasticState(phase)%offsetDeltaState = sizeState-sizeDeltaState                                    ! deltaState occupies latter part of state by definition
 plasticState(phase)%Nslip = Nslip
 plasticState(phase)%Ntwin = Ntwin
 plasticState(phase)%Ntrans= Ntrans

 allocate(plasticState(phase)%aTolState           (sizeState),                source=0.0_pReal)
 allocate(plasticState(phase)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
 allocate(plasticState(phase)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
 allocate(plasticState(phase)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
 allocate(plasticState(phase)%state               (sizeState,NofMyPhase),     source=0.0_pReal)

 allocate(plasticState(phase)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
 if (numerics_integrator == 1_pInt) then
   allocate(plasticState(phase)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
   allocate(plasticState(phase)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
 endif
 if (numerics_integrator == 4_pInt) &
   allocate(plasticState(phase)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
 if (numerics_integrator == 5_pInt) &
   allocate(plasticState(phase)%RKCK45dotState  (6,sizeDotState,NofMyPhase),  source=0.0_pReal)

 allocate(plasticState(phase)%deltaState        (sizeDeltaState,NofMyPhase),  source=0.0_pReal)

end subroutine material_allocatePlasticState


!--------------------------------------------------------------------------------------------------
!> @brief allocates the source state of a phase
!--------------------------------------------------------------------------------------------------
subroutine material_allocateSourceState(phase,of,NofMyPhase,&
                                        sizeState,sizeDotState,sizeDeltaState)
 use numerics, only: &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: &
   phase, &
   of, &
   NofMyPhase, &
   sizeState, sizeDotState,sizeDeltaState

 sourceState(phase)%p(of)%sizeState        = sizeState
 sourceState(phase)%p(of)%sizeDotState     = sizeDotState
 sourceState(phase)%p(of)%sizeDeltaState   = sizeDeltaState
 sourceState(phase)%p(of)%offsetDeltaState = sizeState-sizeDeltaState                               ! deltaState occupies latter part of state by definition

 allocate(sourceState(phase)%p(of)%aTolState           (sizeState),                source=0.0_pReal)
 allocate(sourceState(phase)%p(of)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
 allocate(sourceState(phase)%p(of)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
 allocate(sourceState(phase)%p(of)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
 allocate(sourceState(phase)%p(of)%state               (sizeState,NofMyPhase),     source=0.0_pReal)

 allocate(sourceState(phase)%p(of)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
 if (numerics_integrator == 1_pInt) then
   allocate(sourceState(phase)%p(of)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
   allocate(sourceState(phase)%p(of)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
 endif
 if (numerics_integrator == 4_pInt) &
   allocate(sourceState(phase)%p(of)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
 if (numerics_integrator == 5_pInt) &
   allocate(sourceState(phase)%p(of)%RKCK45dotState  (6,sizeDotState,NofMyPhase),  source=0.0_pReal)

 allocate(sourceState(phase)%p(of)%deltaState        (sizeDeltaState,NofMyPhase),  source=0.0_pReal)

end subroutine material_allocateSourceState


!--------------------------------------------------------------------------------------------------
!> @brief populates the grains
!> @details populates the grains by identifying active microstructure/homogenization pairs,
!! calculates the volume of the grains and deals with texture components
!--------------------------------------------------------------------------------------------------
subroutine material_populateGrains
 use prec, only: &
   dEq
 use math, only: &
   math_RtoEuler, &
   math_EulerToR, &
   math_mul33x33, &
   math_range
 use mesh, only: &
   theMesh, &
   mesh_ipVolume
 use config, only: &
   config_homogenization, &
   config_microstructure, &
   config_deallocate, &
   homogenization_name, &
   microstructure_name 
 use IO, only: &
   IO_error
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
 real(pReal) :: deviation,extreme,rnd
 integer(pInt),         dimension (:,:), allocatable :: Nelems                                      ! counts number of elements in homog, micro array
 type(group_int), dimension (:,:), allocatable :: elemsOfHomogMicro                                 ! lists element number in homog, micro array

 myDebug = debug_level(debug_material)

 allocate(material_volume(homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems),       source=0.0_pReal)
 allocate(material_phase(homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems),        source=0_pInt)
 allocate(material_texture(homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems),      source=0_pInt)
 allocate(material_EulerAngles(3,homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems),source=0.0_pReal)

 allocate(Ngrains(size(config_homogenization),size(config_microstructure)),            source=0_pInt)
 allocate(Nelems (size(config_homogenization),size(config_microstructure)),            source=0_pInt)


!--------------------------------------------------------------------------------------------------
! precounting of elements for each homog/micro pair
 do e = 1_pInt, theMesh%Nelems
   homog = theMesh%homogenizationAt(e)
   micro = theMesh%microstructureAt(e)
   Nelems(homog,micro) = Nelems(homog,micro) + 1_pInt
 enddo
 allocate(elemsOfHomogMicro(size(config_homogenization),size(config_microstructure)))
 do homog = 1,size(config_homogenization)
   do micro = 1,size(config_microstructure)
     if (Nelems(homog,micro) > 0_pInt) then
       allocate(elemsOfHomogMicro(homog,micro)%p(Nelems(homog,micro)))
       elemsOfHomogMicro(homog,micro)%p = 0_pInt
    endif
   enddo
 enddo

!--------------------------------------------------------------------------------------------------
! identify maximum grain count per IP (from element) and find grains per homog/micro pair
 Nelems = 0_pInt                                                                                    ! reuse as counter
 elementLooping: do e = 1_pInt,theMesh%Nelems
   homog = theMesh%homogenizationAt(e)
   micro = theMesh%microstructureAt(e)
   if (homog < 1_pInt .or. homog > size(config_homogenization)) &                                      ! out of bounds
     call IO_error(154_pInt,e,0_pInt,0_pInt)
   if (micro < 1_pInt .or. micro > size(config_microstructure)) &                                      ! out of bounds
     call IO_error(155_pInt,e,0_pInt,0_pInt)
   if (microstructure_elemhomo(micro)) then                                                         ! how many grains are needed at this element?
     dGrains = homogenization_Ngrains(homog)                                                        ! only one set of Ngrains (other IPs are plain copies)
   else
     dGrains = homogenization_Ngrains(homog) * theMesh%elem%nIPs                                     ! each IP has Ngrains
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
     write(6,'(/,a/)') ' MATERIAL grain population'
     write(6,'(a32,1x,a32,1x,a6)') 'homogenization_name','microstructure_name','grain#'
 endif
 homogenizationLoop: do homog = 1_pInt,size(config_homogenization)
   dGrains = homogenization_Ngrains(homog)                                                          ! grain number per material point
   microstructureLoop: do micro = 1_pInt,size(config_microstructure)                                                       ! all pairs of homog and micro
     activePair: if (Ngrains(homog,micro) > 0_pInt) then
       myNgrains = Ngrains(homog,micro)                                                             ! assign short name for total number of grains to populate
       myNconstituents = microstructure_Nconstituents(micro)                                        ! assign short name for number of constituents
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) &
         write(6,'(/,a32,1x,a32,1x,i6)') homogenization_name(homog),microstructure_name(micro),myNgrains


!--------------------------------------------------------------------------------------------------
! calculate volume of each grain

       volumeOfGrain = 0.0_pReal
       grain = 0_pInt

       do hme = 1_pInt, Nelems(homog,micro)
         e = elemsOfHomogMicro(homog,micro)%p(hme)                                                  ! my combination of homog and micro, only perform calculations for elements with homog, micro combinations which is indexed in cpElemsindex
         if (microstructure_elemhomo(micro)) then                                                   ! homogeneous distribution of grains over each element's IPs
           volumeOfGrain(grain+1_pInt:grain+dGrains) = sum(mesh_ipVolume(1:theMesh%elem%nIPs,e))/&
                                                                         real(dGrains,pReal)        ! each grain combines size of all IPs in that element
           grain = grain + dGrains                                                                  ! wind forward by Ngrains@IP
         else
           forall (i = 1_pInt:theMesh%elem%nIPs) &                                                   ! loop over IPs
             volumeOfGrain(grain+(i-1)*dGrains+1_pInt:grain+i*dGrains) = &
               mesh_ipVolume(i,e)/real(dGrains,pReal)                                               ! assign IPvolume/Ngrains@IP to all grains of IP
           grain = grain + theMesh%elem%nIPs * dGrains                                               ! wind forward by Nips*Ngrains@IP
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
         NgrainsOfConstituent(i) = nint(microstructure_fraction(i,micro)*real(myNgrains,pReal),pInt)! do rounding integer conversion
       do while (sum(NgrainsOfConstituent) /= myNgrains)                                            ! total grain count over constituents wrong?
         sgn = sign(1_pInt, myNgrains - sum(NgrainsOfConstituent))                                  ! direction of required change
         extreme = 0.0_pReal
         t = 0_pInt
         do i = 1_pInt,myNconstituents                                                              ! find largest deviator
           deviation = real(sgn,pReal)*log( microstructure_fraction(i,micro) / &
                                           !-------------------------------- &
                                           (real(NgrainsOfConstituent(i),pReal)/real(myNgrains,pReal) ) )
           if (deviation > extreme) then
             extreme = deviation
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
! has texture components
         gauss: do t = 1_pInt,texture_Ngauss(textureID)                                             ! loop over Gauss components
           do g = 1_pInt,int(real(myNorientations,pReal)*texture_Gauss(5,t,textureID),pInt)         ! loop over required grain count
             orientationOfGrain(:,grain+constituentGrain+g) = texture_Gauss(1:3,t,textureID)
           enddo
           constituentGrain = &
           constituentGrain + int(real(myNorientations,pReal)*texture_Gauss(5,t,textureID))         ! advance counter for grains of current constituent
         enddo gauss


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
! shuffle grains within current constituent

         do j = 1_pInt,NgrainsOfConstituent(i)-1_pInt                                               ! walk thru grains of current constituent
           call random_number(rnd)
           t = nint(rnd*real(NgrainsOfConstituent(i)-j,pReal)+real(j,pReal)+0.5_pReal,pInt)       ! select a grain in remaining list
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
         if (microstructure_elemhomo(micro)) then                                                   ! homogeneous distribution of grains over each element's IPs
           m = 1_pInt                                                                               ! process only first IP
         else
           m = theMesh%elem%nIPs
         endif

         do i = 1_pInt, m                                                                           ! loop over necessary IPs
           ip = ip + 1_pInt                                                                         ! keep track of total ip count
           ipGrain = 0_pInt                                                                         ! count number of grains assigned at this IP
           randomOrder = math_range(microstructure_maxNconstituents)                                ! start out with ordered sequence of constituents
           call random_number(rndArray)                                                             ! as many rnd numbers as (max) constituents
           do j = 1_pInt, myNconstituents - 1_pInt                                                  ! loop over constituents ...
             r = nint(rndArray(j)*real(myNconstituents-j,pReal)+real(j,pReal)+0.5_pReal,pInt)       ! ... select one in remaining list
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

         do i = i, theMesh%elem%nIPs                                                                 ! loop over IPs to (possibly) distribute copies from first IP
           material_volume (1_pInt:dGrains,i,e) = material_volume (1_pInt:dGrains,1,e)
           material_phase  (1_pInt:dGrains,i,e) = material_phase  (1_pInt:dGrains,1,e)
           material_texture(1_pInt:dGrains,i,e) = material_texture(1_pInt:dGrains,1,e)
           material_EulerAngles(1:3,1_pInt:dGrains,i,e) = material_EulerAngles(1:3,1_pInt:dGrains,1,e)
         enddo

       enddo
     endif activePair
   enddo microstructureLoop
 enddo homogenizationLoop

 deallocate(texture_transformation)
 deallocate(elemsOfHomogMicro)
 call config_deallocate('material.config/microstructure')

end subroutine material_populateGrains

end module material
