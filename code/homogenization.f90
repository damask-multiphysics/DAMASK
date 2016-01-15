!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Denny Tjahjanto, Max-Planck-Institut für Eisenforschung GmbH
!> @brief homogenization manager, organizing deformation partitioning and stress homogenization
!--------------------------------------------------------------------------------------------------
module homogenization
 use prec, only: &
#ifdef FEM
   tOutputData, &
#endif
   pInt, &
   pReal

!--------------------------------------------------------------------------------------------------
! General variables for the homogenization at a  material point
 implicit none
 private
   real(pReal),   dimension(:,:,:,:),     allocatable, public :: &
   materialpoint_F0, &                                                                              !< def grad of IP at start of FE increment
   materialpoint_F, &                                                                               !< def grad of IP to be reached at end of FE increment
   materialpoint_P                                                                                  !< first P--K stress of IP
 real(pReal),   dimension(:,:,:,:,:,:), allocatable, public ::  &
   materialpoint_dPdF                                                                               !< tangent of first P--K stress at IP
#ifdef FEM
 type(tOutputData), dimension(:),       allocatable, public :: &
   homogOutput
 type(tOutputData), dimension(:,:),     allocatable, public :: &
   crystalliteOutput, &
   phaseOutput
#else
 real(pReal),   dimension(:,:,:),       allocatable, public :: &
   materialpoint_results                                                                            !< results array of material point
#endif
 integer(pInt),                                      public, protected  :: &
   materialpoint_sizeResults, &
   homogenization_maxSizePostResults, &
   thermal_maxSizePostResults, &
   damage_maxSizePostResults, &
   vacancyflux_maxSizePostResults, &
   porosity_maxSizePostResults, &
   hydrogenflux_maxSizePostResults

 real(pReal),   dimension(:,:,:,:),     allocatable, private :: &
   materialpoint_subF0, &                                                                           !< def grad of IP at beginning of homogenization increment
   materialpoint_subF                                                                               !< def grad of IP to be reached at end of homog inc
 real(pReal),   dimension(:,:),         allocatable, private :: &
   materialpoint_subFrac, &
   materialpoint_subStep, &
   materialpoint_subdt
 logical,       dimension(:,:),         allocatable, private :: &
   materialpoint_requested, &
   materialpoint_converged
 logical,       dimension(:,:,:),       allocatable, private :: &
   materialpoint_doneAndHappy

 public ::  &
   homogenization_init, &
   materialpoint_stressAndItsTangent, &
   materialpoint_postResults
 private :: &
   homogenization_partitionDeformation, &
   homogenization_updateState, &
   homogenization_averageStressAndItsTangent, &
   homogenization_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!--------------------------------------------------------------------------------------------------
subroutine homogenization_init
#ifdef HDF
 use hdf5, only: &
   HID_T
 use IO, only : &
   HDF5_mappingHomogenization
#endif
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use math, only: &
   math_I3
 use debug, only: &
   debug_level, &
   debug_homogenization, &
   debug_levelBasic, &
   debug_e, &
   debug_g
 use mesh, only: &
   mesh_maxNips, &
   mesh_NcpElems, &
   mesh_element, &
   FE_Nips, &
   FE_geomtype
#ifdef FEM
 use crystallite, only: &
   crystallite_sizePostResults
#else
 use constitutive, only: &
   constitutive_plasticity_maxSizePostResults, &
   constitutive_source_maxSizePostResults
 use crystallite, only: &
   crystallite_maxSizePostResults
#endif
 use material
 use homogenization_none
 use homogenization_isostrain
 use homogenization_RGC
 use thermal_isothermal
 use thermal_adiabatic
 use thermal_conduction
 use damage_none
 use damage_local
 use damage_nonlocal
 use vacancyflux_isoconc
 use vacancyflux_isochempot
 use vacancyflux_cahnhilliard
 use porosity_none
 use porosity_phasefield
 use hydrogenflux_isoconc
 use hydrogenflux_cahnhilliard
 use IO
 use numerics, only: &
   worldrank

 implicit none
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt) :: e,i,p
 integer(pInt), dimension(:,:), pointer :: thisSize
 integer(pInt), dimension(:)  , pointer :: thisNoutput
 character(len=64), dimension(:,:), pointer :: thisOutput
 character(len=32) :: outputName                                                                    !< name of output, intermediate fix until HDF5 output is ready
 logical :: knownHomogenization, knownThermal, knownDamage, knownVacancyflux, knownPorosity, knownHydrogenflux
#ifdef HDF
 integer(pInt), dimension(:,:), allocatable :: mapping
 integer(pInt), dimension(:), allocatable :: InstancePosition
 allocate(mapping(mesh_ncpelems,4),source=0_pInt)
 allocate(InstancePosition(material_Nhomogenization),source=0_pInt)
#endif

!--------------------------------------------------------------------------------------------------
! parse homogenization from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(homogenization_type == HOMOGENIZATION_NONE_ID)) &
   call homogenization_none_init()
 if (any(homogenization_type == HOMOGENIZATION_ISOSTRAIN_ID)) &
   call homogenization_isostrain_init(FILEUNIT)
 if (any(homogenization_type == HOMOGENIZATION_RGC_ID)) &
   call homogenization_RGC_init(FILEUNIT)


!--------------------------------------------------------------------------------------------------
! parse thermal from config file
 call IO_checkAndRewind(FILEUNIT)
 if (any(thermal_type == THERMAL_isothermal_ID)) &
   call thermal_isothermal_init()
 if (any(thermal_type == THERMAL_adiabatic_ID)) &
   call thermal_adiabatic_init(FILEUNIT)
 if (any(thermal_type == THERMAL_conduction_ID)) &
   call thermal_conduction_init(FILEUNIT)


!--------------------------------------------------------------------------------------------------
! parse damage from config file
 call IO_checkAndRewind(FILEUNIT)
 if (any(damage_type == DAMAGE_none_ID)) &
   call damage_none_init()
 if (any(damage_type == DAMAGE_local_ID)) &
   call damage_local_init(FILEUNIT)
 if (any(damage_type == DAMAGE_nonlocal_ID)) &
   call damage_nonlocal_init(FILEUNIT)


!--------------------------------------------------------------------------------------------------
! parse vacancy transport from config file
 call IO_checkAndRewind(FILEUNIT)
 if (any(vacancyflux_type == VACANCYFLUX_isoconc_ID)) &
   call vacancyflux_isoconc_init()
 if (any(vacancyflux_type == VACANCYFLUX_isochempot_ID)) &
   call vacancyflux_isochempot_init(FILEUNIT)
 if (any(vacancyflux_type == VACANCYFLUX_cahnhilliard_ID)) &
   call vacancyflux_cahnhilliard_init(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! parse porosity from config file
 call IO_checkAndRewind(FILEUNIT)
 if (any(porosity_type == POROSITY_none_ID)) &
   call porosity_none_init()
 if (any(porosity_type == POROSITY_phasefield_ID)) &
   call porosity_phasefield_init(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! parse hydrogen transport from config file
 call IO_checkAndRewind(FILEUNIT)
 if (any(hydrogenflux_type == HYDROGENFLUX_isoconc_ID)) &
   call hydrogenflux_isoconc_init()
 if (any(hydrogenflux_type == HYDROGENFLUX_cahnhilliard_ID)) &
   call hydrogenflux_cahnhilliard_init(FILEUNIT)
 close(FILEUNIT)

!--------------------------------------------------------------------------------------------------
! write description file for homogenization output
 mainProcess2: if (worldrank == 0) then
   call IO_write_jobFile(FILEUNIT,'outputHomogenization')
   do p = 1,material_Nhomogenization
     if (any(material_homog == p)) then
       i = homogenization_typeInstance(p)                                                               ! which instance of this homogenization type
       knownHomogenization = .true.                                                                     ! assume valid
       select case(homogenization_type(p))                                                              ! split per homogenization type
         case (HOMOGENIZATION_NONE_ID)
           outputName = HOMOGENIZATION_NONE_label
           thisNoutput => null()
           thisOutput => null()
           thisSize   => null()
         case (HOMOGENIZATION_ISOSTRAIN_ID)
           outputName = HOMOGENIZATION_ISOSTRAIN_label
           thisNoutput => homogenization_isostrain_Noutput
           thisOutput => homogenization_isostrain_output
           thisSize   => homogenization_isostrain_sizePostResult
         case (HOMOGENIZATION_RGC_ID)
           outputName = HOMOGENIZATION_RGC_label
           thisNoutput => homogenization_RGC_Noutput
           thisOutput => homogenization_RGC_output
           thisSize   => homogenization_RGC_sizePostResult
         case default
           knownHomogenization = .false.
       end select
       write(FILEUNIT,'(/,a,/)')  '['//trim(homogenization_name(p))//']'
       if (knownHomogenization) then
         write(FILEUNIT,'(a)') '(type)'//char(9)//trim(outputName)
         write(FILEUNIT,'(a,i4)') '(ngrains)'//char(9),homogenization_Ngrains(p)
         if (homogenization_type(p) /= HOMOGENIZATION_NONE_ID) then
           do e = 1,thisNoutput(i)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
           enddo
         endif
       endif
       i = thermal_typeInstance(p)                                                                      ! which instance of this thermal type
       knownThermal = .true.                                                                            ! assume valid
       select case(thermal_type(p))                                                                     ! split per thermal type
         case (THERMAL_isothermal_ID)
           outputName = THERMAL_isothermal_label
           thisNoutput => null()
           thisOutput => null()
           thisSize   => null()
         case (THERMAL_adiabatic_ID)
           outputName = THERMAL_adiabatic_label
           thisNoutput => thermal_adiabatic_Noutput
           thisOutput => thermal_adiabatic_output
           thisSize   => thermal_adiabatic_sizePostResult
         case (THERMAL_conduction_ID)
           outputName = THERMAL_conduction_label
           thisNoutput => thermal_conduction_Noutput
           thisOutput => thermal_conduction_output
           thisSize   => thermal_conduction_sizePostResult
         case default
           knownThermal = .false.
       end select
       if (knownThermal) then
         write(FILEUNIT,'(a)') '(thermal)'//char(9)//trim(outputName)
         if (thermal_type(p) /= THERMAL_isothermal_ID) then
           do e = 1,thisNoutput(i)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
           enddo
         endif
       endif
       i = damage_typeInstance(p)                                                                       ! which instance of this damage type
       knownDamage = .true.                                                                             ! assume valid
       select case(damage_type(p))                                                                      ! split per damage type
         case (DAMAGE_none_ID)
           outputName = DAMAGE_none_label
           thisNoutput => null()
           thisOutput => null()
           thisSize   => null()
         case (DAMAGE_local_ID)
           outputName = DAMAGE_local_label
           thisNoutput => damage_local_Noutput
           thisOutput => damage_local_output
           thisSize   => damage_local_sizePostResult
         case (DAMAGE_nonlocal_ID)
           outputName = DAMAGE_nonlocal_label
           thisNoutput => damage_nonlocal_Noutput
           thisOutput => damage_nonlocal_output
           thisSize   => damage_nonlocal_sizePostResult
         case default
           knownDamage = .false.
       end select
       if (knownDamage) then
         write(FILEUNIT,'(a)') '(damage)'//char(9)//trim(outputName)
         if (damage_type(p) /= DAMAGE_none_ID) then
           do e = 1,thisNoutput(i)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
           enddo
         endif
       endif
       i = vacancyflux_typeInstance(p)                                                                  ! which instance of this vacancy flux type
       knownVacancyflux = .true.                                                                        ! assume valid
       select case(vacancyflux_type(p))                                                                 ! split per vacancy flux type
         case (VACANCYFLUX_isoconc_ID)
           outputName = VACANCYFLUX_isoconc_label
           thisNoutput => null()
           thisOutput => null()
           thisSize   => null()
         case (VACANCYFLUX_isochempot_ID)
           outputName = VACANCYFLUX_isochempot_label
           thisNoutput => vacancyflux_isochempot_Noutput
           thisOutput => vacancyflux_isochempot_output
           thisSize   => vacancyflux_isochempot_sizePostResult
         case (VACANCYFLUX_cahnhilliard_ID)
           outputName = VACANCYFLUX_cahnhilliard_label
           thisNoutput => vacancyflux_cahnhilliard_Noutput
           thisOutput => vacancyflux_cahnhilliard_output
           thisSize   => vacancyflux_cahnhilliard_sizePostResult
         case default
           knownVacancyflux = .false.
       end select
       if (knownVacancyflux) then
         write(FILEUNIT,'(a)') '(vacancyflux)'//char(9)//trim(outputName)
         if (vacancyflux_type(p) /= VACANCYFLUX_isoconc_ID) then
           do e = 1,thisNoutput(i)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
           enddo
         endif
       endif
       i = porosity_typeInstance(p)                                                                     ! which instance of this porosity type
       knownPorosity = .true.                                                                           ! assume valid
       select case(porosity_type(p))                                                                    ! split per porosity type
         case (POROSITY_none_ID)
           outputName = POROSITY_none_label
           thisNoutput => null()
           thisOutput => null()
           thisSize   => null()
         case (POROSITY_phasefield_ID)
           outputName = POROSITY_phasefield_label
           thisNoutput => porosity_phasefield_Noutput
           thisOutput => porosity_phasefield_output
           thisSize   => porosity_phasefield_sizePostResult
         case default
           knownPorosity = .false.
       end select
       if (knownPorosity) then
         write(FILEUNIT,'(a)') '(porosity)'//char(9)//trim(outputName)
         if (porosity_type(p) /= POROSITY_none_ID) then
           do e = 1,thisNoutput(i)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
           enddo
         endif
       endif
       i = hydrogenflux_typeInstance(p)                                                                 ! which instance of this hydrogen flux type
       knownHydrogenflux = .true.                                                                       ! assume valid
       select case(hydrogenflux_type(p))                                                                ! split per hydrogen flux type
         case (HYDROGENFLUX_isoconc_ID)
           outputName = HYDROGENFLUX_isoconc_label
           thisNoutput => null()
           thisOutput => null()
           thisSize   => null()
         case (HYDROGENFLUX_cahnhilliard_ID)
           outputName = HYDROGENFLUX_cahnhilliard_label
           thisNoutput => hydrogenflux_cahnhilliard_Noutput
           thisOutput => hydrogenflux_cahnhilliard_output
           thisSize   => hydrogenflux_cahnhilliard_sizePostResult
         case default
           knownHydrogenflux = .false.
       end select
       if (knownHydrogenflux) then
         write(FILEUNIT,'(a)') '(hydrogenflux)'//char(9)//trim(outputName)
         if (hydrogenflux_type(p) /= HYDROGENFLUX_isoconc_ID) then
           do e = 1,thisNoutput(i)
             write(FILEUNIT,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
           enddo
         endif
       endif
     endif
   enddo
   close(FILEUNIT)
 endif mainProcess2

!--------------------------------------------------------------------------------------------------
! allocate and initialize global variables
 allocate(materialpoint_dPdF(3,3,3,3,mesh_maxNips,mesh_NcpElems),       source=0.0_pReal)
 allocate(materialpoint_F0(3,3,mesh_maxNips,mesh_NcpElems),             source=0.0_pReal)
 materialpoint_F0 = spread(spread(math_I3,3,mesh_maxNips),4,mesh_NcpElems)                          ! initialize to identity
 allocate(materialpoint_F(3,3,mesh_maxNips,mesh_NcpElems),              source=0.0_pReal)
 materialpoint_F = materialpoint_F0                                                                 ! initialize to identity
 allocate(materialpoint_subF0(3,3,mesh_maxNips,mesh_NcpElems),          source=0.0_pReal)
 allocate(materialpoint_subF(3,3,mesh_maxNips,mesh_NcpElems),           source=0.0_pReal)
 allocate(materialpoint_P(3,3,mesh_maxNips,mesh_NcpElems),              source=0.0_pReal)
 allocate(materialpoint_subFrac(mesh_maxNips,mesh_NcpElems),            source=0.0_pReal)
 allocate(materialpoint_subStep(mesh_maxNips,mesh_NcpElems),            source=0.0_pReal)
 allocate(materialpoint_subdt(mesh_maxNips,mesh_NcpElems),              source=0.0_pReal)
 allocate(materialpoint_requested(mesh_maxNips,mesh_NcpElems),          source=.false.)
 allocate(materialpoint_converged(mesh_maxNips,mesh_NcpElems),          source=.true.)
 allocate(materialpoint_doneAndHappy(2,mesh_maxNips,mesh_NcpElems),     source=.true.)

!--------------------------------------------------------------------------------------------------
! allocate and initialize global state and postresutls variables
#ifdef HDF
 elementLooping: do e = 1,mesh_NcpElems
   myInstance = homogenization_typeInstance(mesh_element(3,e))
   IpLooping: do i = 1,FE_Nips(FE_geomtype(mesh_element(2,e)))
       InstancePosition(myInstance) = InstancePosition(myInstance)+1_pInt
       mapping(e,1:4) = [instancePosition(myinstance),myinstance,e,i]
   enddo IpLooping
 enddo elementLooping
 call  HDF5_mappingHomogenization(mapping)
#endif

 homogenization_maxSizePostResults = 0_pInt
 thermal_maxSizePostResults        = 0_pInt
 damage_maxSizePostResults         = 0_pInt
 vacancyflux_maxSizePostResults    = 0_pInt
 porosity_maxSizePostResults       = 0_pInt
 hydrogenflux_maxSizePostResults   = 0_pInt
 do p = 1,material_Nhomogenization
   homogenization_maxSizePostResults = max(homogenization_maxSizePostResults,homogState       (p)%sizePostResults)
   thermal_maxSizePostResults        = max(thermal_maxSizePostResults,       thermalState     (p)%sizePostResults)
   damage_maxSizePostResults         = max(damage_maxSizePostResults        ,damageState      (p)%sizePostResults)
   vacancyflux_maxSizePostResults    = max(vacancyflux_maxSizePostResults   ,vacancyfluxState (p)%sizePostResults)
   porosity_maxSizePostResults       = max(porosity_maxSizePostResults      ,porosityState    (p)%sizePostResults)
   hydrogenflux_maxSizePostResults   = max(hydrogenflux_maxSizePostResults  ,hydrogenfluxState(p)%sizePostResults)
 enddo

#ifdef FEM
 allocate(homogOutput      (material_Nhomogenization                         ))
 allocate(crystalliteOutput(material_Ncrystallite,  homogenization_maxNgrains))
 allocate(phaseOutput      (material_Nphase,        homogenization_maxNgrains))
 do p = 1, material_Nhomogenization
   homogOutput(p)%sizeResults = homogState       (p)%sizePostResults + &
                                thermalState     (p)%sizePostResults + &
                                damageState      (p)%sizePostResults + &
                                vacancyfluxState (p)%sizePostResults + &
                                porosityState    (p)%sizePostResults + &
                                hydrogenfluxState(p)%sizePostResults
   homogOutput(p)%sizeIpCells = count(material_homog==p)
   allocate(homogOutput(p)%output(homogOutput(p)%sizeResults,homogOutput(p)%sizeIpCells))
 enddo
 do p = 1, material_Ncrystallite; do e = 1, homogenization_maxNgrains
   crystalliteOutput(p,e)%sizeResults = crystallite_sizePostResults(p)
   crystalliteOutput(p,e)%sizeIpCells = count(microstructure_crystallite(mesh_element(4,:)) == p .and. &
                                              homogenization_Ngrains    (mesh_element(3,:)) >= e)*mesh_maxNips
   allocate(crystalliteOutput(p,e)%output(crystalliteOutput(p,e)%sizeResults,crystalliteOutput(p,e)%sizeIpCells))
 enddo; enddo
 do p = 1, material_Nphase; do e = 1, homogenization_maxNgrains
   phaseOutput(p,e)%sizeResults = plasticState    (p)%sizePostResults + &
                                  sum(sourceState (p)%p(:)%sizePostResults)
   phaseOutput(p,e)%sizeIpCells = count(material_phase(e,:,:) == p)
   allocate(phaseOutput(p,e)%output(phaseOutput(p,e)%sizeResults,phaseOutput(p,e)%sizeIpCells))
 enddo; enddo
#else
 materialpoint_sizeResults = 1 &                                                                    ! grain count
                           + 1 + homogenization_maxSizePostResults &                                ! homogSize & homogResult
                               + thermal_maxSizePostResults        &
                               + damage_maxSizePostResults         &
                               + vacancyflux_maxSizePostResults    &
                               + porosity_maxSizePostResults       &
                               + hydrogenflux_maxSizePostResults   &
                           + homogenization_maxNgrains * (1 + crystallite_maxSizePostResults &      ! crystallite size & crystallite results
                                                        + 1 + constitutive_plasticity_maxSizePostResults &     ! constitutive size & constitutive results
                                                            + constitutive_source_maxSizePostResults)
 allocate(materialpoint_results(materialpoint_sizeResults,mesh_maxNips,mesh_NcpElems))
#endif

 mainProcess: if (worldrank == 0) then
   write(6,'(/,a)')   ' <<<+-  homogenization init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt) then
#ifdef TODO
   write(6,'(a32,1x,7(i8,1x))')   'homogenization_state0:          ', shape(homogenization_state0)
   write(6,'(a32,1x,7(i8,1x))')   'homogenization_subState0:       ', shape(homogenization_subState0)
   write(6,'(a32,1x,7(i8,1x))')   'homogenization_state:           ', shape(homogenization_state)
#endif
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_dPdF:             ', shape(materialpoint_dPdF)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_F0:               ', shape(materialpoint_F0)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_F:                ', shape(materialpoint_F)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subF0:            ', shape(materialpoint_subF0)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subF:             ', shape(materialpoint_subF)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_P:                ', shape(materialpoint_P)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subFrac:          ', shape(materialpoint_subFrac)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subStep:          ', shape(materialpoint_subStep)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subdt:            ', shape(materialpoint_subdt)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_requested:        ', shape(materialpoint_requested)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_converged:        ', shape(materialpoint_converged)
   write(6,'(a32,1x,7(i8,1x),/)') 'materialpoint_doneAndHappy:     ', shape(materialpoint_doneAndHappy)
#ifndef FEM
   write(6,'(a32,1x,7(i8,1x),/)') 'materialpoint_results:          ', shape(materialpoint_results)
#endif
   write(6,'(a32,1x,7(i8,1x))')   'maxSizePostResults: ', homogenization_maxSizePostResults
 endif
 flush(6)

 if (debug_g < 1 .or. debug_g > homogenization_Ngrains(mesh_element(3,debug_e))) &
   call IO_error(602_pInt,ext_msg='component (grain)')

end subroutine homogenization_init


!--------------------------------------------------------------------------------------------------
!> @brief  parallelized calculation of stress and corresponding tangent at material points
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_stressAndItsTangent(updateJaco,dt)
 use numerics, only: &
   subStepMinHomog, &
   subStepSizeHomog, &
   stepIncreaseHomog, &
   nHomog, &
   nMPstate
 use math, only: &
   math_transpose33
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP, &
   terminallyIll
 use mesh, only: &
   mesh_element
 use material, only: &
   plasticState, &
   sourceState, &
   homogState, &
   thermalState, &
   damageState, &
   vacancyfluxState, &
   porosityState, &
   hydrogenfluxState, &
   phase_Nsources, &
   mappingHomogenization, &
   phaseAt, phasememberAt, &
   homogenization_Ngrains
 use crystallite, only: &
   crystallite_F0, &
   crystallite_Fp0, &
   crystallite_Fp, &
   crystallite_Fi0, &
   crystallite_Fi, &
   crystallite_Lp0, &
   crystallite_Lp, &
   crystallite_Li0, &
   crystallite_Li, &
   crystallite_dPdF, &
   crystallite_dPdF0, &
   crystallite_Tstar0_v, &
   crystallite_Tstar_v, &
   crystallite_partionedF0, &
   crystallite_partionedF, &
   crystallite_partionedFp0, &
   crystallite_partionedLp0, &
   crystallite_partionedFi0, &
   crystallite_partionedLi0, &
   crystallite_partioneddPdF0, &
   crystallite_partionedTstar0_v, &
   crystallite_dt, &
   crystallite_requested, &
   crystallite_converged, &
   crystallite_stressAndItsTangent, &
   crystallite_orientations
 use debug, only: &
   debug_level, &
   debug_homogenization, &
   debug_levelBasic, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_MaterialpointLoopDistribution, &
   debug_MaterialpointStateLoopDistribution

 implicit none
 real(pReal), intent(in) :: dt                                                                      !< time increment
 logical,     intent(in) :: updateJaco                                                              !< initiating Jacobian update
 integer(pInt) :: &
   NiterationHomog, &
   NiterationMPstate, &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e, &                                                                                             !< element number
   mySource, &
   myNgrains

!--------------------------------------------------------------------------------------------------
! initialize to starting condition
 if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
     write(6,'(/a,i5,1x,i2)') '<< HOMOG >> Material Point start at el ip ', debug_e, debug_i

     write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< HOMOG >> F0', &
                                     math_transpose33(materialpoint_F0(1:3,1:3,debug_i,debug_e))
     write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< HOMOG >> F', &
                                     math_transpose33(materialpoint_F(1:3,1:3,debug_i,debug_e))
   !$OMP END CRITICAL (write2out)
 endif

!--------------------------------------------------------------------------------------------------
! initialize restoration points of ...
 do e = FEsolving_execElem(1),FEsolving_execElem(2)
   myNgrains = homogenization_Ngrains(mesh_element(3,e))
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e); do g = 1,myNgrains

     plasticState    (phaseAt(g,i,e))%partionedState0(:,phasememberAt(g,i,e)) = &
     plasticState    (phaseAt(g,i,e))%state0(         :,phasememberAt(g,i,e))
     do mySource = 1_pInt, phase_Nsources(phaseAt(g,i,e))
       sourceState(phaseAt(g,i,e))%p(mySource)%partionedState0(:,phasememberAt(g,i,e)) = &
       sourceState(phaseAt(g,i,e))%p(mySource)%state0(         :,phasememberAt(g,i,e))
     enddo

     crystallite_partionedFp0(1:3,1:3,g,i,e) = crystallite_Fp0(1:3,1:3,g,i,e)                       ! ...plastic def grads
     crystallite_partionedLp0(1:3,1:3,g,i,e) = crystallite_Lp0(1:3,1:3,g,i,e)                       ! ...plastic velocity grads
     crystallite_partionedFi0(1:3,1:3,g,i,e) = crystallite_Fi0(1:3,1:3,g,i,e)                       ! ...intermediate def grads
     crystallite_partionedLi0(1:3,1:3,g,i,e) = crystallite_Li0(1:3,1:3,g,i,e)                       ! ...intermediate velocity grads
     crystallite_partioneddPdF0(1:3,1:3,1:3,1:3,g,i,e) = crystallite_dPdF0(1:3,1:3,1:3,1:3,g,i,e)   ! ...stiffness
     crystallite_partionedF0(1:3,1:3,g,i,e) = crystallite_F0(1:3,1:3,g,i,e)                         ! ...def grads
     crystallite_partionedTstar0_v(1:6,g,i,e) = crystallite_Tstar0_v(1:6,g,i,e)                     ! ...2nd PK stress

   enddo; enddo
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e))
     materialpoint_subF0(1:3,1:3,i,e) = materialpoint_F0(1:3,1:3,i,e)                               ! ...def grad
     materialpoint_subFrac(i,e) = 0.0_pReal
     materialpoint_subStep(i,e) = 1.0_pReal/subStepSizeHomog                                        ! <<added to adopt flexibility in cutback size>>
     materialpoint_converged(i,e) = .false.                                                         ! pretend failed step of twice the required size
     materialpoint_requested(i,e) = .true.                                                          ! everybody requires calculation
   endforall
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
     homogState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
       homogState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
       homogState(mappingHomogenization(2,i,e))%State0(   :,mappingHomogenization(1,i,e))           ! ...internal homogenization state
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
     thermalState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
       thermalState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
       thermalState(mappingHomogenization(2,i,e))%State0(   :,mappingHomogenization(1,i,e))         ! ...internal thermal state
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
     damageState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
       damageState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
       damageState(mappingHomogenization(2,i,e))%State0(   :,mappingHomogenization(1,i,e))          ! ...internal damage state
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
     vacancyfluxState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
       vacancyfluxState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
       vacancyfluxState(mappingHomogenization(2,i,e))%State0(   :,mappingHomogenization(1,i,e))     ! ...internal vacancy transport state
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
     porosityState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
       porosityState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
       porosityState(mappingHomogenization(2,i,e))%State0(   :,mappingHomogenization(1,i,e))        ! ...internal porosity state
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
     hydrogenfluxState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
       hydrogenfluxState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
       hydrogenfluxState(mappingHomogenization(2,i,e))%State0(   :,mappingHomogenization(1,i,e))    ! ...internal hydrogen transport state
 enddo
 NiterationHomog = 0_pInt

 cutBackLooping: do while (.not. terminallyIll .and. &
      any(materialpoint_subStep(:,FEsolving_execELem(1):FEsolving_execElem(2)) > subStepMinHomog))

   !$OMP PARALLEL DO PRIVATE(myNgrains)
   elementLooping1: do e = FEsolving_execElem(1),FEsolving_execElem(2)
     myNgrains = homogenization_Ngrains(mesh_element(3,e))
     IpLooping1: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)

       converged: if ( materialpoint_converged(i,e) ) then
#ifndef _OPENMP
         if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt &
            .and. ((e == debug_e .and. i == debug_i) &
                   .or. .not. iand(debug_level(debug_homogenization),debug_levelSelective) /= 0_pInt)) then
           write(6,'(a,1x,f12.8,1x,a,1x,f12.8,1x,a,i8,1x,i2/)') '<< HOMOG >> winding forward from', &
             materialpoint_subFrac(i,e), 'to current materialpoint_subFrac', &
             materialpoint_subFrac(i,e)+materialpoint_subStep(i,e),'in materialpoint_stressAndItsTangent at el ip',e,i
         endif
#endif

!---------------------------------------------------------------------------------------------------
! calculate new subStep and new subFrac
         materialpoint_subFrac(i,e) = materialpoint_subFrac(i,e) + materialpoint_subStep(i,e)
         !$OMP FLUSH(materialpoint_subFrac)
         materialpoint_subStep(i,e) = min(1.0_pReal-materialpoint_subFrac(i,e), &
                                          stepIncreaseHomog*materialpoint_subStep(i,e))             ! introduce flexibility for step increase/acceleration
         !$OMP FLUSH(materialpoint_subStep)

         steppingNeeded: if (materialpoint_subStep(i,e) > subStepMinHomog) then

           ! wind forward grain starting point of...
           crystallite_partionedF0(1:3,1:3,1:myNgrains,i,e) =  &
              crystallite_partionedF(1:3,1:3,1:myNgrains,i,e)                                       ! ...def grads

           crystallite_partionedFp0(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_Fp(1:3,1:3,1:myNgrains,i,e)                                                ! ...plastic def grads

           crystallite_partionedLp0(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_Lp(1:3,1:3,1:myNgrains,i,e)                                                ! ...plastic velocity grads

           crystallite_partionedFi0(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_Fi(1:3,1:3,1:myNgrains,i,e)                                                ! ...intermediate def grads

           crystallite_partionedLi0(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_Li(1:3,1:3,1:myNgrains,i,e)                                                ! ...intermediate velocity grads

           crystallite_partioneddPdF0(1:3,1:3,1:3,1:3,1:myNgrains,i,e) = &
             crystallite_dPdF(1:3,1:3,1:3,1:3,1:myNgrains,i,e)                                      ! ...stiffness

           crystallite_partionedTstar0_v(1:6,1:myNgrains,i,e) = &
             crystallite_Tstar_v(1:6,1:myNgrains,i,e)                                               ! ...2nd PK stress

           do g = 1,myNgrains
             plasticState    (phaseAt(g,i,e))%partionedState0(:,phasememberAt(g,i,e)) = &
             plasticState    (phaseAt(g,i,e))%state(          :,phasememberAt(g,i,e))
             do mySource = 1_pInt, phase_Nsources(phaseAt(g,i,e))
               sourceState(phaseAt(g,i,e))%p(mySource)%partionedState0(:,phasememberAt(g,i,e)) = &
               sourceState(phaseAt(g,i,e))%p(mySource)%state(          :,phasememberAt(g,i,e))
             enddo
           enddo

           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             homogState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               homogState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
               homogState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e))   ! ...internal homogenization state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             thermalState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               thermalState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
               thermalState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e)) ! ...internal thermal state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             damageState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               damageState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
               damageState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e))  ! ...internal damage state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             vacancyfluxState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               vacancyfluxState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
               vacancyfluxState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e))! ...internal vacancy transport state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             porosityState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               porosityState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
               porosityState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e))! ...internal porosity state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             hydrogenfluxState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               hydrogenfluxState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) = &
               hydrogenfluxState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e))! ...internal hydrogen transport state
           materialpoint_subF0(1:3,1:3,i,e) = materialpoint_subF(1:3,1:3,i,e)                       ! ...def grad
           !$OMP FLUSH(materialpoint_subF0)
         elseif (materialpoint_requested(i,e)) then steppingNeeded                                  ! already at final time (??)
           if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt) then
             !$OMP CRITICAL (distributionHomog)
               debug_MaterialpointLoopDistribution(min(nHomog+1,NiterationHomog)) = &
               debug_MaterialpointLoopDistribution(min(nHomog+1,NiterationHomog)) + 1
             !$OMP END CRITICAL (distributionHomog)
           endif
         endif steppingNeeded

       else converged
         if ( (myNgrains == 1_pInt .and. materialpoint_subStep(i,e) <= 1.0 ) .or. &                 ! single grain already tried internal subStepping in crystallite
              subStepSizeHomog * materialpoint_subStep(i,e) <=  subStepMinHomog ) then              ! would require too small subStep
                                                                                                    ! cutback makes no sense
           !$OMP FLUSH(terminallyIll)
           if (.not. terminallyIll) then                                                            ! so first signals terminally ill...
             !$OMP CRITICAL (write2out)
               write(6,*) 'Integration point ', i,' at element ', e, ' terminally ill'
             !$OMP END CRITICAL (write2out)
           endif
           !$OMP CRITICAL (setTerminallyIll)
             terminallyIll = .true.                                                                 ! ...and kills all others
           !$OMP END CRITICAL (setTerminallyIll)
         else                                                                                       ! cutback makes sense
           materialpoint_subStep(i,e) = subStepSizeHomog * materialpoint_subStep(i,e)               ! crystallite had severe trouble, so do a significant cutback
           !$OMP FLUSH(materialpoint_subStep)

#ifndef _OPENMP
           if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt &
              .and. ((e == debug_e .and. i == debug_i) &
                    .or. .not. iand(debug_level(debug_homogenization), debug_levelSelective) /= 0_pInt)) then
             write(6,'(a,1x,f12.8,a,i8,1x,i2/)') &
               '<< HOMOG >> cutback step in materialpoint_stressAndItsTangent with new materialpoint_subStep:',&
               materialpoint_subStep(i,e),' at el ip',e,i
           endif
#endif

!--------------------------------------------------------------------------------------------------
! restore...
           crystallite_Fp(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_partionedFp0(1:3,1:3,1:myNgrains,i,e)                                      ! ...plastic def grads
           crystallite_Lp(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_partionedLp0(1:3,1:3,1:myNgrains,i,e)                                      ! ...plastic velocity grads
           crystallite_Fi(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_partionedFi0(1:3,1:3,1:myNgrains,i,e)                                      ! ...intermediate def grads
           crystallite_Li(1:3,1:3,1:myNgrains,i,e) = &
             crystallite_partionedLi0(1:3,1:3,1:myNgrains,i,e)                                      ! ...intermediate velocity grads
           crystallite_dPdF(1:3,1:3,1:3,1:3,1:myNgrains,i,e) = &
             crystallite_partioneddPdF0(1:3,1:3,1:3,1:3,1:myNgrains,i,e)                            ! ...stiffness
           crystallite_Tstar_v(1:6,1:myNgrains,i,e) = &
              crystallite_partionedTstar0_v(1:6,1:myNgrains,i,e)                                    ! ...2nd PK stress
           do g = 1, myNgrains
             plasticState    (phaseAt(g,i,e))%state(          :,phasememberAt(g,i,e)) = &
             plasticState    (phaseAt(g,i,e))%partionedState0(:,phasememberAt(g,i,e))
             do mySource = 1_pInt, phase_Nsources(phaseAt(g,i,e))
               sourceState(phaseAt(g,i,e))%p(mySource)%state(          :,phasememberAt(g,i,e)) = &
               sourceState(phaseAt(g,i,e))%p(mySource)%partionedState0(:,phasememberAt(g,i,e))
             enddo
           enddo
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             homogState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               homogState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e)) = &
               homogState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e))   ! ...internal homogenization state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             thermalState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               thermalState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e)) = &
               thermalState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e)) ! ...internal thermal state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             damageState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               damageState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e)) = &
               damageState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e))  ! ...internal damage state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             vacancyfluxState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               vacancyfluxState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e)) = &
               vacancyfluxState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e))! ...internal vacancy transport state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             porosityState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               porosityState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e)) = &
               porosityState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e))! ...internal porosity state
           forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), &
             hydrogenfluxState(mappingHomogenization(2,i,e))%sizeState > 0_pInt) &
               hydrogenfluxState(mappingHomogenization(2,i,e))%State(    :,mappingHomogenization(1,i,e)) = &
               hydrogenfluxState(mappingHomogenization(2,i,e))%subState0(:,mappingHomogenization(1,i,e))! ...internal hydrogen transport state
         endif
       endif converged

       if (materialpoint_subStep(i,e) > subStepMinHomog) then
         materialpoint_requested(i,e) = .true.
         materialpoint_subF(1:3,1:3,i,e) = materialpoint_subF0(1:3,1:3,i,e) + &
                                         materialpoint_subStep(i,e) * (materialpoint_F(1:3,1:3,i,e) - materialpoint_F0(1:3,1:3,i,e))
         materialpoint_subdt(i,e) = materialpoint_subStep(i,e) * dt
         materialpoint_doneAndHappy(1:2,i,e) = [.false.,.true.]
       endif
     enddo IpLooping1
   enddo elementLooping1
   !$OMP END PARALLEL DO

   NiterationMPstate = 0_pInt

   convergenceLooping: do while (.not. terminallyIll .and. &
             any(            materialpoint_requested(:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                 .and. .not. materialpoint_doneAndHappy(1,:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                ) .and. &
             NiterationMPstate < nMPstate)
     NiterationMPstate = NiterationMPstate + 1

!--------------------------------------------------------------------------------------------------
! deformation partitioning
! based on materialpoint_subF0,.._subF,crystallite_partionedF0, and homogenization_state,
! results in crystallite_partionedF
     !$OMP PARALLEL DO PRIVATE(myNgrains)
     elementLooping2: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       myNgrains = homogenization_Ngrains(mesh_element(3,e))
       IpLooping2: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
         if (      materialpoint_requested(i,e) .and. &                                             ! process requested but...
             .not. materialpoint_doneAndHappy(1,i,e)) then                                          ! ...not yet done material points
           call homogenization_partitionDeformation(i,e)                                            ! partition deformation onto constituents
           crystallite_dt(1:myNgrains,i,e) = materialpoint_subdt(i,e)                               ! propagate materialpoint dt to grains
           crystallite_requested(1:myNgrains,i,e) = .true.                                          ! request calculation for constituents
         else
           crystallite_requested(1:myNgrains,i,e) = .false.                                         ! calculation for constituents not required anymore
         endif
       enddo IpLooping2
     enddo elementLooping2
     !$OMP END PARALLEL DO

!--------------------------------------------------------------------------------------------------
! crystallite integration
! based on crystallite_partionedF0,.._partionedF
! incrementing by crystallite_dt
     call crystallite_stressAndItsTangent(updateJaco)                                                ! request stress and tangent calculation for constituent grains

!--------------------------------------------------------------------------------------------------
! state update
     !$OMP PARALLEL DO
     elementLooping3: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       IpLooping3: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
         if (      materialpoint_requested(i,e) .and. &
             .not. materialpoint_doneAndHappy(1,i,e)) then
           if (.not. all(crystallite_converged(:,i,e))) then
             materialpoint_doneAndHappy(1:2,i,e) = [.true.,.false.]
             materialpoint_converged(i,e) = .false.
           else
             materialpoint_doneAndHappy(1:2,i,e) = homogenization_updateState(i,e)
             materialpoint_converged(i,e) = all(materialpoint_doneAndHappy(1:2,i,e))                  ! converged if done and happy
           endif
           !$OMP FLUSH(materialpoint_converged)
           if (materialpoint_converged(i,e)) then
             if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt) then
               !$OMP CRITICAL (distributionMPState)
                 debug_MaterialpointStateLoopdistribution(NiterationMPstate) = &
                   debug_MaterialpointStateLoopdistribution(NiterationMPstate) + 1_pInt
               !$OMP END CRITICAL (distributionMPState)
             endif
           endif
         endif
       enddo IpLooping3
     enddo elementLooping3
     !$OMP END PARALLEL DO

   enddo convergenceLooping

   NiterationHomog = NiterationHomog + 1_pInt

 enddo cutBackLooping

 if (.not. terminallyIll ) then
   call crystallite_orientations()                                                                  ! calculate crystal orientations
   !$OMP PARALLEL DO
   elementLooping4: do e = FEsolving_execElem(1),FEsolving_execElem(2)
     IpLooping4: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       call homogenization_averageStressAndItsTangent(i,e)
     enddo IpLooping4
   enddo elementLooping4
   !$OMP END PARALLEL DO
 else
   !$OMP CRITICAL (write2out)
   write(6,'(/,a,/)') '<< HOMOG >> Material Point terminally ill'
   !$OMP END CRITICAL (write2out)
 endif

end subroutine materialpoint_stressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief parallelized calculation of result array at material points
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_postResults
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP
 use mesh, only: &
   mesh_element
 use material, only: &
   mappingHomogenization, &
#ifdef FEM
   phaseAt, phasememberAt, &
   homogenization_maxNgrains, &
   material_Ncrystallite, &
   material_Nphase, &
#else
   homogState, &
   thermalState, &
   damageState, &
   vacancyfluxState, &
   porosityState, &
   hydrogenfluxState, &
#endif
   plasticState, &
   sourceState, &
   material_phase, &
   homogenization_Ngrains, &
   microstructure_crystallite
 use constitutive, only: &
#ifdef FEM
   constitutive_plasticity_maxSizePostResults, &
   constitutive_source_maxSizePostResults, &
#endif
   constitutive_postResults
 use crystallite, only: &
#ifdef FEM
   crystallite_maxSizePostResults, &
#endif
   crystallite_sizePostResults, &
   crystallite_postResults

 implicit none
 integer(pInt) :: &
   thePos, &
   theSize, &
   myNgrains, &
   myCrystallite, &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number
#ifdef FEM
 integer(pInt) :: &
   myHomog, &
   myPhase, &
   crystalliteCtr(material_Ncrystallite,  homogenization_maxNgrains), &
   phaseCtr      (material_Nphase,        homogenization_maxNgrains)
 real(pReal), dimension(1+crystallite_maxSizePostResults + &
                        1+constitutive_plasticity_maxSizePostResults + &
                          constitutive_source_maxSizePostResults) :: &
   crystalliteResults



 crystalliteCtr = 0_pInt; phaseCtr = 0_pInt
 elementLooping: do e = FEsolving_execElem(1),FEsolving_execElem(2)
   myNgrains = homogenization_Ngrains(mesh_element(3,e))
   myCrystallite = microstructure_crystallite(mesh_element(4,e))
   IpLooping: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     myHomog = mappingHomogenization(2,i,e)
     thePos =  mappingHomogenization(1,i,e)
     homogOutput(myHomog)%output(1: &
                                 homogOutput(myHomog)%sizeResults, &
                                 thePos) = homogenization_postResults(i,e)

     grainLooping :do g = 1,myNgrains
       myPhase = phaseAt(g,i,e)
       crystalliteResults(1:1+crystallite_sizePostResults(myCrystallite) + &
                            1+plasticState(myPhase)%sizePostResults + &
                              sum(sourceState(myPhase)%p(:)%sizePostResults)) = crystallite_postResults(g,i,e)
       if (microstructure_crystallite(mesh_element(4,e)) == myCrystallite .and. &
           homogenization_Ngrains    (mesh_element(3,e)) >= g) then
         crystalliteCtr(myCrystallite,g) = crystalliteCtr(myCrystallite,g) + 1_pInt
         crystalliteOutput(myCrystallite,g)% &
           output(1:crystalliteOutput(myCrystallite,g)%sizeResults,crystalliteCtr(myCrystallite,g)) = &
             crystalliteResults(2:1+crystalliteOutput(myCrystallite,g)%sizeResults)
       endif
       if (material_phase(g,i,e) == myPhase) then
         phaseCtr(myPhase,g) = phaseCtr(myPhase,g) + 1_pInt
         phaseOutput(myPhase,g)% &
           output(1:phaseOutput(myPhase,g)%sizeResults,phaseCtr(myPhase,g)) = &
             crystalliteResults(3 + crystalliteOutput(myCrystallite,g)%sizeResults: &
                                1 + crystalliteOutput(myCrystallite,g)%sizeResults + &
                                1 + plasticState    (myphase)%sizePostResults + &
                                    sum(sourceState(myphase)%p(:)%sizePostResults))
       endif
     enddo grainLooping
   enddo IpLooping
 enddo elementLooping
#else

 !$OMP PARALLEL DO PRIVATE(myNgrains,myCrystallite,thePos,theSize)
   elementLooping: do e = FEsolving_execElem(1),FEsolving_execElem(2)
     myNgrains = homogenization_Ngrains(mesh_element(3,e))
     myCrystallite = microstructure_crystallite(mesh_element(4,e))
     IpLooping: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       thePos = 0_pInt

       theSize = homogState       (mappingHomogenization(2,i,e))%sizePostResults &
               + thermalState     (mappingHomogenization(2,i,e))%sizePostResults &
               + damageState      (mappingHomogenization(2,i,e))%sizePostResults &
               + vacancyfluxState (mappingHomogenization(2,i,e))%sizePostResults &
               + porosityState    (mappingHomogenization(2,i,e))%sizePostResults &
               + hydrogenfluxState(mappingHomogenization(2,i,e))%sizePostResults
       materialpoint_results(thePos+1,i,e) = real(theSize,pReal)                                    ! tell size of homogenization results
       thePos = thePos + 1_pInt

       if (theSize > 0_pInt) then                                                                   ! any homogenization results to mention?
         materialpoint_results(thePos+1:thePos+theSize,i,e) = homogenization_postResults(i,e)       ! tell homogenization results
         thePos = thePos + theSize
       endif

       materialpoint_results(thePos+1,i,e) = real(myNgrains,pReal)                                  ! tell number of grains at materialpoint
       thePos = thePos + 1_pInt

       grainLooping :do g = 1,myNgrains
         theSize = 1 + crystallite_sizePostResults(myCrystallite) + &
                   1 + plasticState    (material_phase(g,i,e))%sizePostResults + &                    !ToDo
                       sum(sourceState(material_phase(g,i,e))%p(:)%sizePostResults)
         materialpoint_results(thePos+1:thePos+theSize,i,e) = crystallite_postResults(g,i,e)        ! tell crystallite results
         thePos = thePos + theSize
       enddo grainLooping
     enddo IpLooping
   enddo elementLooping
 !$OMP END PARALLEL DO
#endif

end subroutine materialpoint_postResults


!--------------------------------------------------------------------------------------------------
!> @brief  partition material point def grad onto constituents
!--------------------------------------------------------------------------------------------------
subroutine homogenization_partitionDeformation(ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_type, &
   homogenization_maxNgrains, &
   HOMOGENIZATION_NONE_ID, &
   HOMOGENIZATION_ISOSTRAIN_ID, &
   HOMOGENIZATION_RGC_ID
 use crystallite, only: &
   crystallite_partionedF
 use homogenization_isostrain, only: &
   homogenization_isostrain_partitionDeformation
 use homogenization_RGC, only: &
   homogenization_RGC_partitionDeformation

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number

 chosenHomogenization: select case(homogenization_type(mesh_element(3,el)))

   case (HOMOGENIZATION_NONE_ID) chosenHomogenization
     crystallite_partionedF(1:3,1:3,1:homogenization_maxNgrains,ip,el) = 0.0_pReal
     crystallite_partionedF(1:3,1:3,1:1,ip,el) = &
       spread(materialpoint_subF(1:3,1:3,ip,el),3,1)

   case (HOMOGENIZATION_ISOSTRAIN_ID) chosenHomogenization
     call homogenization_isostrain_partitionDeformation(&
                          crystallite_partionedF(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                          materialpoint_subF(1:3,1:3,ip,el),&
                          el)
   case (HOMOGENIZATION_RGC_ID) chosenHomogenization
     call homogenization_RGC_partitionDeformation(&
                         crystallite_partionedF(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                         materialpoint_subF(1:3,1:3,ip,el),&
                         ip, &
                         el)
 end select chosenHomogenization

end subroutine homogenization_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief update the internal state of the homogenization scheme and tell whether "done" and
!> "happy" with result
!--------------------------------------------------------------------------------------------------
function homogenization_updateState(ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_type, &
   thermal_type, &
   damage_type, &
   vacancyflux_type, &
   homogenization_maxNgrains, &
   HOMOGENIZATION_RGC_ID, &
   THERMAL_adiabatic_ID, &
   DAMAGE_local_ID, &
   VACANCYFLUX_isochempot_ID
 use crystallite, only: &
   crystallite_P, &
   crystallite_dPdF, &
   crystallite_partionedF,&
   crystallite_partionedF0
 use homogenization_RGC, only: &
   homogenization_RGC_updateState
 use thermal_adiabatic, only: &
   thermal_adiabatic_updateState
 use damage_local, only: &
   damage_local_updateState
 use vacancyflux_isochempot, only: &
   vacancyflux_isochempot_updateState

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number
 logical, dimension(2) :: homogenization_updateState

 homogenization_updateState = .true.
 chosenHomogenization: select case(homogenization_type(mesh_element(3,el)))
   case (HOMOGENIZATION_RGC_ID) chosenHomogenization
     homogenization_updateState = &
       homogenization_updateState .and. &
        homogenization_RGC_updateState(crystallite_P(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                                       crystallite_partionedF(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                                       crystallite_partionedF0(1:3,1:3,1:homogenization_maxNgrains,ip,el),&
                                       materialpoint_subF(1:3,1:3,ip,el),&
                                       materialpoint_subdt(ip,el), &
                                       crystallite_dPdF(1:3,1:3,1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                                       ip, &
                                       el)
   case default chosenHomogenization
     homogenization_updateState = &
       homogenization_updateState .and. [.true., .true.]
 end select chosenHomogenization

 chosenThermal: select case (thermal_type(mesh_element(3,el)))
   case (THERMAL_adiabatic_ID) chosenThermal
     homogenization_updateState = &
       homogenization_updateState .and. &
       thermal_adiabatic_updateState(materialpoint_subdt(ip,el), &
                                     ip, &
                                     el)
   case default chosenThermal
     homogenization_updateState = &
       homogenization_updateState .and. [.true., .true.]
 end select chosenThermal

 chosenDamage: select case (damage_type(mesh_element(3,el)))
   case (DAMAGE_local_ID) chosenDamage
     homogenization_updateState = &
       homogenization_updateState .and. &
       damage_local_updateState(materialpoint_subdt(ip,el), &
                                ip, &
                                el)
   case default chosenDamage
     homogenization_updateState = &
       homogenization_updateState .and. [.true., .true.]
 end select chosenDamage

 chosenVacancyflux: select case (vacancyflux_type(mesh_element(3,el)))
   case (VACANCYFLUX_isochempot_ID) chosenVacancyflux
     homogenization_updateState = &
       homogenization_updateState .and. &
       vacancyflux_isochempot_updateState(materialpoint_subdt(ip,el), &
                                          ip, &
                                          el)
   case default chosenVacancyflux
     homogenization_updateState = &
       homogenization_updateState .and. [.true., .true.]
 end select chosenVacancyflux

end function homogenization_updateState


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities
!--------------------------------------------------------------------------------------------------
subroutine homogenization_averageStressAndItsTangent(ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_type, &
   homogenization_maxNgrains, &
   HOMOGENIZATION_NONE_ID, &
   HOMOGENIZATION_ISOSTRAIN_ID, &
   HOMOGENIZATION_RGC_ID
 use crystallite, only: &
   crystallite_P,crystallite_dPdF
 use homogenization_isostrain, only: &
   homogenization_isostrain_averageStressAndItsTangent
 use homogenization_RGC, only: &
   homogenization_RGC_averageStressAndItsTangent

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number

 chosenHomogenization: select case(homogenization_type(mesh_element(3,el)))
   case (HOMOGENIZATION_NONE_ID) chosenHomogenization
       materialpoint_P(1:3,1:3,ip,el) = sum(crystallite_P(1:3,1:3,1:1,ip,el),3)
       materialpoint_dPdF(1:3,1:3,1:3,1:3,ip,el) &
        = sum(crystallite_dPdF(1:3,1:3,1:3,1:3,1:1,ip,el),5)

   case (HOMOGENIZATION_ISOSTRAIN_ID) chosenHomogenization
     call homogenization_isostrain_averageStressAndItsTangent(&
       materialpoint_P(1:3,1:3,ip,el), &
       materialpoint_dPdF(1:3,1:3,1:3,1:3,ip,el),&
       crystallite_P(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
       crystallite_dPdF(1:3,1:3,1:3,1:3,1:homogenization_maxNgrains,ip,el), &
       el)
   case (HOMOGENIZATION_RGC_ID) chosenHomogenization
     call homogenization_RGC_averageStressAndItsTangent(&
       materialpoint_P(1:3,1:3,ip,el), &
       materialpoint_dPdF(1:3,1:3,1:3,1:3,ip,el),&
       crystallite_P(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
       crystallite_dPdF(1:3,1:3,1:3,1:3,1:homogenization_maxNgrains,ip,el), &
       el)
 end select chosenHomogenization

end subroutine homogenization_averageStressAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion. call only,
!> if homogenization_sizePostResults(i,e) > 0 !!
!--------------------------------------------------------------------------------------------------
function homogenization_postResults(ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   mappingHomogenization, &
   homogState, &
   thermalState, &
   damageState, &
   vacancyfluxState, &
   porosityState, &
   hydrogenfluxState, &
   homogenization_type, &
   thermal_type, &
   damage_type, &
   vacancyflux_type, &
   porosity_type, &
   hydrogenflux_type, &
   HOMOGENIZATION_NONE_ID, &
   HOMOGENIZATION_ISOSTRAIN_ID, &
   HOMOGENIZATION_RGC_ID, &
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
   HYDROGENFLUX_cahnhilliard_ID
 use homogenization_isostrain, only: &
   homogenization_isostrain_postResults
 use homogenization_RGC, only: &
   homogenization_RGC_postResults
 use thermal_adiabatic, only: &
   thermal_adiabatic_postResults
 use thermal_conduction, only: &
   thermal_conduction_postResults
 use damage_local, only: &
   damage_local_postResults
 use damage_nonlocal, only: &
   damage_nonlocal_postResults
 use vacancyflux_isochempot, only: &
   vacancyflux_isochempot_postResults
 use vacancyflux_cahnhilliard, only: &
   vacancyflux_cahnhilliard_postResults
 use porosity_phasefield, only: &
   porosity_phasefield_postResults
 use hydrogenflux_cahnhilliard, only: &
   hydrogenflux_cahnhilliard_postResults

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number
 real(pReal), dimension(  homogState       (mappingHomogenization(2,ip,el))%sizePostResults &
                        + thermalState     (mappingHomogenization(2,ip,el))%sizePostResults &
                        + damageState      (mappingHomogenization(2,ip,el))%sizePostResults &
                        + vacancyfluxState (mappingHomogenization(2,ip,el))%sizePostResults &
                        + porosityState    (mappingHomogenization(2,ip,el))%sizePostResults &
                        + hydrogenfluxState(mappingHomogenization(2,ip,el))%sizePostResults) :: &
   homogenization_postResults
 integer(pInt) :: &
   startPos, endPos

 homogenization_postResults = 0.0_pReal

 startPos = 1_pInt
 endPos   = homogState(mappingHomogenization(2,ip,el))%sizePostResults
 chosenHomogenization: select case (homogenization_type(mesh_element(3,el)))
   case (HOMOGENIZATION_NONE_ID) chosenHomogenization

   case (HOMOGENIZATION_ISOSTRAIN_ID) chosenHomogenization
     homogenization_postResults(startPos:endPos) = &
       homogenization_isostrain_postResults(&
                                  ip, &
                                  el, &
                                  materialpoint_P(1:3,1:3,ip,el), &
                                  materialpoint_F(1:3,1:3,ip,el))
   case (HOMOGENIZATION_RGC_ID) chosenHomogenization
     homogenization_postResults(startPos:endPos) = &
       homogenization_RGC_postResults(&
                                  ip, &
                                  el, &
                                  materialpoint_P(1:3,1:3,ip,el), &
                                  materialpoint_F(1:3,1:3,ip,el))
 end select chosenHomogenization

 startPos = endPos + 1_pInt
 endPos   = endPos + thermalState(mappingHomogenization(2,ip,el))%sizePostResults
 chosenThermal: select case (thermal_type(mesh_element(3,el)))
   case (THERMAL_isothermal_ID) chosenThermal

   case (THERMAL_adiabatic_ID) chosenThermal
     homogenization_postResults(startPos:endPos) = &
       thermal_adiabatic_postResults(ip, el)
   case (THERMAL_conduction_ID) chosenThermal
     homogenization_postResults(startPos:endPos) = &
       thermal_conduction_postResults(ip, el)
 end select chosenThermal

 startPos = endPos + 1_pInt
 endPos   = endPos + damageState(mappingHomogenization(2,ip,el))%sizePostResults
 chosenDamage: select case (damage_type(mesh_element(3,el)))
   case (DAMAGE_none_ID) chosenDamage

   case (DAMAGE_local_ID) chosenDamage
     homogenization_postResults(startPos:endPos) = &
       damage_local_postResults(ip, el)

   case (DAMAGE_nonlocal_ID) chosenDamage
     homogenization_postResults(startPos:endPos) = &
       damage_nonlocal_postResults(ip, el)
 end select chosenDamage

 startPos = endPos + 1_pInt
 endPos   = endPos + vacancyfluxState(mappingHomogenization(2,ip,el))%sizePostResults
 chosenVacancyflux: select case (vacancyflux_type(mesh_element(3,el)))
   case (VACANCYFLUX_isoconc_ID) chosenVacancyflux

   case (VACANCYFLUX_isochempot_ID) chosenVacancyflux
     homogenization_postResults(startPos:endPos) = &
       vacancyflux_isochempot_postResults(ip, el)
   case (VACANCYFLUX_cahnhilliard_ID) chosenVacancyflux
     homogenization_postResults(startPos:endPos) = &
       vacancyflux_cahnhilliard_postResults(ip, el)
 end select chosenVacancyflux

 startPos = endPos + 1_pInt
 endPos   = endPos + porosityState(mappingHomogenization(2,ip,el))%sizePostResults
 chosenPorosity: select case (porosity_type(mesh_element(3,el)))
   case (POROSITY_none_ID) chosenPorosity

   case (POROSITY_phasefield_ID) chosenPorosity
     homogenization_postResults(startPos:endPos) = &
       porosity_phasefield_postResults(ip, el)
 end select chosenPorosity

 startPos = endPos + 1_pInt
 endPos   = endPos + hydrogenfluxState(mappingHomogenization(2,ip,el))%sizePostResults
 chosenHydrogenflux: select case (hydrogenflux_type(mesh_element(3,el)))
   case (HYDROGENFLUX_isoconc_ID) chosenHydrogenflux

   case (HYDROGENFLUX_cahnhilliard_ID) chosenHydrogenflux
     homogenization_postResults(startPos:endPos) = &
       hydrogenflux_cahnhilliard_postResults(ip, el)
 end select chosenHydrogenflux

end function homogenization_postResults

end module homogenization
