!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for conservative transport of vacancy concentration field
!> @details to be done
!--------------------------------------------------------------------------------------------------
module vacancyflux_cahnhilliard
 use prec, only: &
   pReal, &
   pInt, &
   p_vec

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   vacancyflux_cahnhilliard_sizePostResults                                                           !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   vacancyflux_cahnhilliard_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   vacancyflux_cahnhilliard_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   vacancyflux_cahnhilliard_Noutput                                                                   !< number of outputs per instance of this damage 

 real(pReal),                         dimension(:),           allocatable,        private :: &
   vacancyflux_cahnhilliard_flucAmplitude

 type(p_vec),                         dimension(:),           allocatable,        private :: &
   vacancyflux_cahnhilliard_thermalFluc

 real(pReal),                                                 parameter,           private :: &
   kB = 1.3806488e-23_pReal                                                                          !< Boltzmann constant in J/Kelvin

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 vacancyConc_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,        private :: & 
   vacancyflux_cahnhilliard_outputID                                                                  !< ID of each post result output


 public :: &
   vacancyflux_cahnhilliard_init, &
   vacancyflux_cahnhilliard_getSourceAndItsTangent, &
   vacancyflux_cahnhilliard_getMobility33, &
   vacancyflux_cahnhilliard_getDiffusion33, &
   vacancyflux_cahnhilliard_getChemPotAndItsTangent, &
   vacancyflux_cahnhilliard_putVacancyConcAndItsRate, &
   vacancyflux_cahnhilliard_postResults
 private :: &
   vacancyflux_cahnhilliard_getFormationEnergy, &
   vacancyflux_cahnhilliard_getEntropicCoeff, &
   vacancyflux_cahnhilliard_KinematicChemPotAndItsTangent  

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine vacancyflux_cahnhilliard_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
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
   vacancyflux_type, &
   vacancyflux_typeInstance, &
   homogenization_Noutput, &
   VACANCYFLUX_cahnhilliard_label, &
   VACANCYFLUX_cahnhilliard_ID, &
   material_homog, &  
   mappingHomogenization, &
   vacancyfluxState, &
   vacancyfluxMapping, &
   vacancyConc, &
   vacancyConcRate, &
   vacancyflux_initialCv, &
   material_partHomogenization, &
   material_partPhase
 use numerics,only: &
   worldrank

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: maxNinstance,mySize=0_pInt,section,instance,o,offset
 integer(pInt) :: sizeState
 integer(pInt) :: NofMyHomog   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  vacancyflux_'//VACANCYFLUX_cahnhilliard_label//' init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(vacancyflux_type == VACANCYFLUX_cahnhilliard_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 allocate(vacancyflux_cahnhilliard_sizePostResults(maxNinstance),                               source=0_pInt)
 allocate(vacancyflux_cahnhilliard_sizePostResult (maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(vacancyflux_cahnhilliard_output         (maxval(homogenization_Noutput),maxNinstance))
          vacancyflux_cahnhilliard_output = ''
 allocate(vacancyflux_cahnhilliard_outputID       (maxval(homogenization_Noutput),maxNinstance),source=undefined_ID)
 allocate(vacancyflux_cahnhilliard_Noutput        (maxNinstance),                               source=0_pInt) 

 allocate(vacancyflux_cahnhilliard_flucAmplitude  (maxNinstance))
 allocate(vacancyflux_cahnhilliard_thermalFluc    (maxNinstance))

 rewind(fileUnit)
 section = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partHomogenization)! wind forward to <homogenization>
   line = IO_read(fileUnit)
 enddo
 
 parsingHomog: do while (trim(line) /= IO_EOF)                                                      ! read through sections of homog part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif   
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next homog section
     section = section + 1_pInt                                                                     ! advance homog section counter
     cycle                                                                                          ! skip to next line
   endif

   if (section > 0_pInt ) then; if (vacancyflux_type(section) == VACANCYFLUX_cahnhilliard_ID) then  ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran

     instance = vacancyflux_typeInstance(section)                                                   ! which instance of my vacancyflux is present homog
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case ('vacancyconc')
             vacancyflux_cahnhilliard_Noutput(instance) = vacancyflux_cahnhilliard_Noutput(instance) + 1_pInt
             vacancyflux_cahnhilliard_outputID(vacancyflux_cahnhilliard_Noutput(instance),instance) = vacancyConc_ID
             vacancyflux_cahnhilliard_output(vacancyflux_cahnhilliard_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
          end select

       case ('vacancyflux_flucamplitude')
         vacancyflux_cahnhilliard_flucAmplitude(instance) = IO_floatValue(line,chunkPos,2_pInt)
         
     end select
   endif; endif
 enddo parsingHomog
 
 initializeInstances: do section = 1_pInt, size(vacancyflux_type)
   if (vacancyflux_type(section) == VACANCYFLUX_cahnhilliard_ID) then
     NofMyHomog=count(material_homog==section)
     instance = vacancyflux_typeInstance(section)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,vacancyflux_cahnhilliard_Noutput(instance)
       select case(vacancyflux_cahnhilliard_outputID(o,instance))
         case(vacancyConc_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          vacancyflux_cahnhilliard_sizePostResult(o,instance) = mySize
          vacancyflux_cahnhilliard_sizePostResults(instance)  = vacancyflux_cahnhilliard_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop

! allocate state arrays
     sizeState = 0_pInt
     vacancyfluxState(section)%sizeState = sizeState
     vacancyfluxState(section)%sizePostResults = vacancyflux_cahnhilliard_sizePostResults(instance)
     allocate(vacancyfluxState(section)%state0   (sizeState,NofMyHomog))
     allocate(vacancyfluxState(section)%subState0(sizeState,NofMyHomog))
     allocate(vacancyfluxState(section)%state    (sizeState,NofMyHomog))
     
     allocate(vacancyflux_cahnhilliard_thermalFluc(instance)%p(NofMyHomog))
     do offset = 1_pInt, NofMyHomog
       call random_number(vacancyflux_cahnhilliard_thermalFluc(instance)%p(offset))
       vacancyflux_cahnhilliard_thermalFluc(instance)%p(offset) = &
         1.0_pReal - &
         vacancyflux_cahnhilliard_flucAmplitude(instance)* &
         (vacancyflux_cahnhilliard_thermalFluc(instance)%p(offset) - 0.5_pReal)
     enddo  

     nullify(vacancyfluxMapping(section)%p)
     vacancyfluxMapping(section)%p => mappingHomogenization(1,:,:)
     deallocate(vacancyConc    (section)%p)
     allocate  (vacancyConc    (section)%p(NofMyHomog), source=vacancyflux_initialCv(section))
     deallocate(vacancyConcRate(section)%p)
     allocate  (vacancyConcRate(section)%p(NofMyHomog), source=0.0_pReal)
     
   endif
 
 enddo initializeInstances
 
end subroutine vacancyflux_cahnhilliard_init

!--------------------------------------------------------------------------------------------------
!> @brief  calculates homogenized vacancy driving forces  
!--------------------------------------------------------------------------------------------------
subroutine vacancyflux_cahnhilliard_getSourceAndItsTangent(CvDot, dCvDot_dCv, Cv, ip, el)
 use material, only: &
   homogenization_Ngrains, &
   mappingHomogenization, &
   phaseAt, phasememberAt, &
   phase_source, &
   phase_Nsources, &
   SOURCE_vacancy_phenoplasticity_ID, &
   SOURCE_vacancy_irradiation_ID, &
   SOURCE_vacancy_thermalfluc_ID
 use source_vacancy_phenoplasticity, only: &
   source_vacancy_phenoplasticity_getRateAndItsTangent
 use source_vacancy_irradiation, only: &
   source_vacancy_irradiation_getRateAndItsTangent
 use source_vacancy_thermalfluc, only: &
   source_vacancy_thermalfluc_getRateAndItsTangent
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   Cv
 integer(pInt) :: &
   phase, &
   grain, &
   source
 real(pReal) :: &
   CvDot, dCvDot_dCv, localCvDot, dLocalCvDot_dCv  

 CvDot = 0.0_pReal
 dCvDot_dCv = 0.0_pReal
 do grain = 1, homogenization_Ngrains(mappingHomogenization(2,ip,el))
   phase = phaseAt(grain,ip,el)
   do source = 1_pInt, phase_Nsources(phase)
     select case(phase_source(source,phase))                                                   
       case (SOURCE_vacancy_phenoplasticity_ID)
        call source_vacancy_phenoplasticity_getRateAndItsTangent  (localCvDot, dLocalCvDot_dCv, grain, ip,  el)

       case (SOURCE_vacancy_irradiation_ID)
        call source_vacancy_irradiation_getRateAndItsTangent  (localCvDot, dLocalCvDot_dCv, grain, ip,  el)

       case (SOURCE_vacancy_thermalfluc_ID)
        call source_vacancy_thermalfluc_getRateAndItsTangent(localCvDot, dLocalCvDot_dCv, grain, ip,  el)

     end select
     CvDot = CvDot + localCvDot
     dCvDot_dCv = dCvDot_dCv + dLocalCvDot_dCv
   enddo  
 enddo
 
 CvDot = CvDot/homogenization_Ngrains(mappingHomogenization(2,ip,el))
 dCvDot_dCv = dCvDot_dCv/homogenization_Ngrains(mappingHomogenization(2,ip,el))
 
end subroutine vacancyflux_cahnhilliard_getSourceAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized vacancy mobility tensor in reference configuration
!--------------------------------------------------------------------------------------------------
function vacancyflux_cahnhilliard_getMobility33(ip,el)
 use lattice, only: &
   lattice_vacancyfluxMobility33
 use material, only: &
   homogenization_Ngrains, &
   material_phase
 use mesh, only: &
   mesh_element
 use crystallite, only: &
   crystallite_push33ToRef

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   vacancyflux_cahnhilliard_getMobility33
 integer(pInt) :: &
   grain
  
 vacancyflux_cahnhilliard_getMobility33 = 0.0_pReal
 do grain = 1, homogenization_Ngrains(mesh_element(3,el))
   vacancyflux_cahnhilliard_getMobility33 = vacancyflux_cahnhilliard_getMobility33 + &
    crystallite_push33ToRef(grain,ip,el,lattice_vacancyfluxMobility33(:,:,material_phase(grain,ip,el)))
 enddo

 vacancyflux_cahnhilliard_getMobility33 = &
   vacancyflux_cahnhilliard_getMobility33/ &
   homogenization_Ngrains(mesh_element(3,el))
 
end function vacancyflux_cahnhilliard_getMobility33
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized vacancy diffusion tensor in reference configuration
!--------------------------------------------------------------------------------------------------
function vacancyflux_cahnhilliard_getDiffusion33(ip,el)
 use lattice, only: &
   lattice_vacancyfluxDiffusion33
 use material, only: &
   homogenization_Ngrains, &
   material_phase
 use mesh, only: &
   mesh_element
 use crystallite, only: &
   crystallite_push33ToRef

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   vacancyflux_cahnhilliard_getDiffusion33
 integer(pInt) :: &
   grain
  
 vacancyflux_cahnhilliard_getDiffusion33 = 0.0_pReal
 do grain = 1, homogenization_Ngrains(mesh_element(3,el))
   vacancyflux_cahnhilliard_getDiffusion33 = vacancyflux_cahnhilliard_getDiffusion33 + &
    crystallite_push33ToRef(grain,ip,el,lattice_vacancyfluxDiffusion33(:,:,material_phase(grain,ip,el)))
 enddo

 vacancyflux_cahnhilliard_getDiffusion33 = &
   vacancyflux_cahnhilliard_getDiffusion33/ &
   homogenization_Ngrains(mesh_element(3,el))
 
end function vacancyflux_cahnhilliard_getDiffusion33
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized vacancy formation energy
!--------------------------------------------------------------------------------------------------
real(pReal) function vacancyflux_cahnhilliard_getFormationEnergy(ip,el)
 use lattice, only: &
   lattice_vacancyFormationEnergy, &
   lattice_vacancyVol, &
   lattice_vacancySurfaceEnergy
 use material, only: &
   homogenization_Ngrains, &
   material_phase
 use mesh, only: &
   mesh_element

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 integer(pInt) :: &
   grain
  
 vacancyflux_cahnhilliard_getFormationEnergy = 0.0_pReal
 do grain = 1, homogenization_Ngrains(mesh_element(3,el))
   vacancyflux_cahnhilliard_getFormationEnergy = vacancyflux_cahnhilliard_getFormationEnergy + &
    lattice_vacancyFormationEnergy(material_phase(grain,ip,el))/ &
    lattice_vacancyVol(material_phase(grain,ip,el))/ &
    lattice_vacancySurfaceEnergy(material_phase(grain,ip,el))
 enddo

 vacancyflux_cahnhilliard_getFormationEnergy = &
   vacancyflux_cahnhilliard_getFormationEnergy/ &
   homogenization_Ngrains(mesh_element(3,el))
 
end function vacancyflux_cahnhilliard_getFormationEnergy
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized vacancy entropy coefficient
!--------------------------------------------------------------------------------------------------
real(pReal) function vacancyflux_cahnhilliard_getEntropicCoeff(ip,el)
 use lattice, only: &
   lattice_vacancyVol, &
   lattice_vacancySurfaceEnergy
 use material, only: &
   homogenization_Ngrains, &
   material_homog, &
   material_phase, &
   temperature, &
   thermalMapping

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 integer(pInt) :: &
   grain
  
 vacancyflux_cahnhilliard_getEntropicCoeff = 0.0_pReal
 do grain = 1, homogenization_Ngrains(material_homog(ip,el))
   vacancyflux_cahnhilliard_getEntropicCoeff = vacancyflux_cahnhilliard_getEntropicCoeff + &
    kB/ &
    lattice_vacancyVol(material_phase(grain,ip,el))/ &
    lattice_vacancySurfaceEnergy(material_phase(grain,ip,el))
 enddo

 vacancyflux_cahnhilliard_getEntropicCoeff = &
   vacancyflux_cahnhilliard_getEntropicCoeff* &
   temperature(material_homog(ip,el))%p(thermalMapping(material_homog(ip,el))%p(ip,el))/ &
   homogenization_Ngrains(material_homog(ip,el))
 
end function vacancyflux_cahnhilliard_getEntropicCoeff
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized kinematic contribution to chemical potential
!--------------------------------------------------------------------------------------------------
subroutine vacancyflux_cahnhilliard_KinematicChemPotAndItsTangent(KPot, dKPot_dCv, Cv, ip, el)
 use lattice, only: &
   lattice_vacancySurfaceEnergy
 use material, only: &
   homogenization_Ngrains, &
   material_homog, &
   phase_kinematics, &
   phase_Nkinematics, &
   material_phase, &
   KINEMATICS_vacancy_strain_ID
 use crystallite, only: &
   crystallite_Tstar_v, &
   crystallite_Fi0, &
   crystallite_Fi
 use kinematics_vacancy_strain, only: &
   kinematics_vacancy_strain_ChemPotAndItsTangent

 implicit none
 integer(pInt), intent(in)  :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in)  :: &
   Cv
 real(pReal),   intent(out) :: &
   KPot, dKPot_dCv
 real(pReal) :: &
   my_KPot, my_dKPot_dCv
 integer(pInt) :: &
   grain, kinematics
  
 KPot = 0.0_pReal
 dKPot_dCv = 0.0_pReal
 do grain = 1_pInt,homogenization_Ngrains(material_homog(ip,el)) 
   do kinematics = 1_pInt, phase_Nkinematics(material_phase(grain,ip,el))
     select case (phase_kinematics(kinematics,material_phase(grain,ip,el)))
       case (KINEMATICS_vacancy_strain_ID)
         call kinematics_vacancy_strain_ChemPotAndItsTangent(my_KPot, my_dKPot_dCv, &
                                                             crystallite_Tstar_v(1:6,grain,ip,el), &
                                                             crystallite_Fi0(1:3,1:3,grain,ip,el), &
                                                             crystallite_Fi (1:3,1:3,grain,ip,el), &
                                                             grain,ip, el)

       case default
         my_KPot = 0.0_pReal
         my_dKPot_dCv = 0.0_pReal
     
     end select
     KPot = KPot + my_KPot/lattice_vacancySurfaceEnergy(material_phase(grain,ip,el))
     dKPot_dCv = dKPot_dCv + my_dKPot_dCv/lattice_vacancySurfaceEnergy(material_phase(grain,ip,el))
   enddo
 enddo 
 
 KPot = KPot/homogenization_Ngrains(material_homog(ip,el))  
 dKPot_dCv = dKPot_dCv/homogenization_Ngrains(material_homog(ip,el))  

end subroutine vacancyflux_cahnhilliard_KinematicChemPotAndItsTangent
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized chemical potential and its tangent
!--------------------------------------------------------------------------------------------------
subroutine vacancyflux_cahnhilliard_getChemPotAndItsTangent(ChemPot,dChemPot_dCv,Cv,ip,el)
 use numerics, only: &
   vacancyBoundPenalty, &
   vacancyPolyOrder
 use material, only: &
   mappingHomogenization, &
   vacancyflux_typeInstance, &
   porosity, &
   porosityMapping

 implicit none
 integer(pInt), intent(in)  :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in)  :: &
   Cv
 real(pReal),   intent(out) :: &
   ChemPot, &
   dChemPot_dCv
 real(pReal) :: &
   VoidPhaseFrac, kBT, KPot, dKPot_dCv  
 integer(pInt) :: &
   homog, o  
  
 homog = mappingHomogenization(2,ip,el)
 VoidPhaseFrac = porosity(homog)%p(porosityMapping(homog)%p(ip,el))
 kBT = vacancyflux_cahnhilliard_getEntropicCoeff(ip,el)
 
 ChemPot = vacancyflux_cahnhilliard_getFormationEnergy(ip,el)    
 dChemPot_dCv = 0.0_pReal
 do o = 1_pInt, vacancyPolyOrder
   ChemPot = ChemPot + kBT*((2.0_pReal*Cv - 1.0_pReal)**real(2_pInt*o-1_pInt,pReal))/ &
                       real(2_pInt*o-1_pInt,pReal)
   dChemPot_dCv = dChemPot_dCv + 2.0_pReal*kBT*(2.0_pReal*Cv - 1.0_pReal)**real(2_pInt*o-2_pInt,pReal) 
 enddo
 
 ChemPot =   VoidPhaseFrac*VoidPhaseFrac*ChemPot &
           - 2.0_pReal*(1.0_pReal - Cv)*(1.0_pReal - VoidPhaseFrac)*(1.0_pReal - VoidPhaseFrac)

 dChemPot_dCv =   VoidPhaseFrac*VoidPhaseFrac*dChemPot_dCv &
                + 2.0_pReal*(1.0_pReal - VoidPhaseFrac)*(1.0_pReal - VoidPhaseFrac)

 call vacancyflux_cahnhilliard_KinematicChemPotAndItsTangent(KPot, dKPot_dCv, Cv, ip, el)
 ChemPot = ChemPot + KPot 
 dChemPot_dCv = dChemPot_dCv + dKPot_dCv
 
 if (Cv < 0.0_pReal) then
   ChemPot = ChemPot - 3.0_pReal*vacancyBoundPenalty*Cv*Cv
   dChemPot_dCv = dChemPot_dCv - 6.0_pReal*vacancyBoundPenalty*Cv
 elseif (Cv > 1.0_pReal) then
   ChemPot = ChemPot + 3.0_pReal*vacancyBoundPenalty*(1.0_pReal - Cv)*(1.0_pReal - Cv)
   dChemPot_dCv = dChemPot_dCv - 6.0_pReal*vacancyBoundPenalty*(1.0_pReal - Cv)
 endif
 
 ChemPot = ChemPot* &
           vacancyflux_cahnhilliard_thermalFluc(vacancyflux_typeInstance(homog))%p(mappingHomogenization(1,ip,el))
 dChemPot_dCv = dChemPot_dCv* &
           vacancyflux_cahnhilliard_thermalFluc(vacancyflux_typeInstance(homog))%p(mappingHomogenization(1,ip,el))                 
 
end subroutine vacancyflux_cahnhilliard_getChemPotAndItsTangent
 
!--------------------------------------------------------------------------------------------------
!> @brief updated vacancy concentration and its rate with solution from transport PDE
!--------------------------------------------------------------------------------------------------
subroutine vacancyflux_cahnhilliard_putVacancyConcAndItsRate(Cv,Cvdot,ip,el)
 use material, only: &
   mappingHomogenization, &
   vacancyConc, &
   vacancyConcRate, &
   vacancyfluxMapping

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   Cv, &
   Cvdot
 integer(pInt) :: &
   homog, &
   offset  
 
 homog  = mappingHomogenization(2,ip,el)
 offset = vacancyfluxMapping(homog)%p(ip,el)
 vacancyConc    (homog)%p(offset) = Cv
 vacancyConcRate(homog)%p(offset) = Cvdot

end subroutine vacancyflux_cahnhilliard_putVacancyConcAndItsRate
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of vacancy transport results
!--------------------------------------------------------------------------------------------------
function vacancyflux_cahnhilliard_postResults(ip,el)
 use material, only: &
   mappingHomogenization, &
   vacancyflux_typeInstance, &
   vacancyConc, &
   vacancyfluxMapping

 implicit none
 integer(pInt),              intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(vacancyflux_cahnhilliard_sizePostResults(vacancyflux_typeInstance(mappingHomogenization(2,ip,el)))) :: &
   vacancyflux_cahnhilliard_postResults

 integer(pInt) :: &
   instance, homog, offset, o, c
   
 homog     = mappingHomogenization(2,ip,el)
 offset    = vacancyfluxMapping(homog)%p(ip,el)
 instance  = vacancyflux_typeInstance(homog)

 c = 0_pInt
 vacancyflux_cahnhilliard_postResults = 0.0_pReal

 do o = 1_pInt,vacancyflux_cahnhilliard_Noutput(instance)
    select case(vacancyflux_cahnhilliard_outputID(o,instance))
 
      case (vacancyConc_ID)
        vacancyflux_cahnhilliard_postResults(c+1_pInt) = vacancyConc(homog)%p(offset)
        c = c + 1
    end select
 enddo
end function vacancyflux_cahnhilliard_postResults

end module vacancyflux_cahnhilliard
