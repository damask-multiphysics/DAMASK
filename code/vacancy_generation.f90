!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for plastically generated vacancy concentrations
!> @details to be done
!--------------------------------------------------------------------------------------------------
module vacancy_generation
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   vacancy_generation_sizePostResults                                                           !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   vacancy_generation_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   vacancy_generation_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   vacancy_generation_Noutput                                                                   !< number of outputs per instance of this damage 

 real(pReal),                         dimension(:),           allocatable,         public :: &
   vacancy_generation_aTol, &
   vacancy_generation_freq, &
   vacancy_generation_formationEnergy, &
   vacancy_generation_diffusionEnergy, &
   vacancy_generation_diffusionCoeff0, &                                                      !< the temperature-independent pre-exponential of diffusion coefficient D_0
   vacancy_generation_stressCoeff, &
   vacancy_generation_jogHeight, &                                                              !< the height of jogs in Burgers vectors
   vacancy_generation_jogSeparation, &                                                          !< the jog seperation
   vacancy_generation_nLatticeSites, &                                                          !< the number of lattice sites per unit volume
   vacancy_generation_burgersVec, &                                                             !< the Burgers vector
   vacancy_generation_dislocationCoeff, &
   vacancy_generation_equilibConcentration                                                      !< the equilibrium concentration of vacancy

 real(pReal),                         dimension(:),           allocatable,         public :: &
   pore_nucleation_surfaceEnergy, &                                                             !< surface energy of metal which controls the necleation of pores
   pore_nucleation_atomVolume, &                                                                !< the volume of atom
   pore_nucleation_shellThickness, &                                                            !< the thickness of spherical shell surrounding the pore
   pore_nucleation_concentrationCoeff0                                                          !< the pre-exponential of equilibrium concentration of critical pore

 real(pReal),                                                 parameter,           private :: &
   kB = 1.38e-23_pReal                                                                          !< Boltzmann constant in J/Kelvin

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 vacancy_concentration_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   vacancy_generation_outputID                                                                  !< ID of each post result output


 public :: &
   vacancy_generation_init, &
   vacancy_generation_stateInit, &
   vacancy_generation_aTolState, &
   vacancy_generation_dotState, &
   vacancy_generation_getConcentration, &
   vacancy_generation_putConcentration, &
   vacancy_generation_getVacancyDiffusion33, &
   vacancy_generation_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine vacancy_generation_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
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
   phase_vacancy, &
   phase_vacancyInstance, &
   phase_Noutput, &
   LOCAL_VACANCY_GENERATION_label, &
   LOCAL_VACANCY_generation_ID, &
   material_phase, &  
   vacancyState, &
   MATERIAL_partPhase
 use numerics,only: &
   worldrank, &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), parameter :: MAXNCHUNKS = 7_pInt
 integer(pInt), dimension(1+2*MAXNCHUNKS) :: positions
 integer(pInt) :: maxNinstance,mySize=0_pInt,phase,instance,o
 integer(pInt) :: sizeState, sizeDotState
 integer(pInt) :: NofMyPhase   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  vacancy_'//LOCAL_VACANCY_GENERATION_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(phase_vacancy == LOCAL_VACANCY_generation_ID),pInt)
 if (maxNinstance == 0_pInt) return
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(vacancy_generation_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(vacancy_generation_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(vacancy_generation_output(maxval(phase_Noutput),maxNinstance))
          vacancy_generation_output = ''
 allocate(vacancy_generation_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(vacancy_generation_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(vacancy_generation_aTol(maxNinstance),                                source=0.0_pReal) 
 allocate(vacancy_generation_freq(maxNinstance),                                source=0.0_pReal) 
 allocate(vacancy_generation_formationEnergy(maxNinstance),                     source=0.0_pReal) 
 allocate(vacancy_generation_diffusionEnergy(maxNinstance),                     source=0.0_pReal) 
 allocate(vacancy_generation_stressCoeff(maxNinstance),                         source=0.0_pReal) 
 allocate(vacancy_generation_jogHeight(maxNinstance),                           source=0.0_pReal)
 allocate(vacancy_generation_jogSeparation(maxNinstance),                       source=0.0_pReal) 
 allocate(vacancy_generation_nLatticeSites(maxNinstance),                       source=0.0_pReal)
 allocate(vacancy_generation_burgersVec(maxNinstance),                          source=0.0_pReal)
 allocate(vacancy_generation_diffusionCoeff0(maxNinstance),                     source=0.0_pReal)
 allocate(vacancy_generation_equilibConcentration(maxNinstance),                source=0.0_pReal)
 
 allocate(vacancy_generation_dislocationCoeff(maxNinstance),                    source=0.0_pReal)
 
 allocate(pore_nucleation_surfaceEnergy(maxNinstance),                          source=0.0_pReal)
 allocate(pore_nucleation_atomVolume(maxNinstance),                             source=0.0_pReal)
 allocate(pore_nucleation_shellThickness(maxNinstance),                         source=0.0_pReal)
 allocate(pore_nucleation_concentrationCoeff0(maxNinstance),                    source=0.0_pReal)

 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= MATERIAL_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo
 
 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif   
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next phase section
     phase = phase + 1_pInt                                                                         ! advance phase section counter
     cycle                                                                                          ! skip to next line
   endif

   if (phase > 0_pInt ) then; if (phase_vacancy(phase) == LOCAL_VACANCY_generation_ID) then               ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran

     instance = phase_vacancyInstance(phase)                                                     ! which instance of my vacancy is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('vacancy_concentration')
             vacancy_generation_Noutput(instance) = vacancy_generation_Noutput(instance) + 1_pInt
             vacancy_generation_outputID(vacancy_generation_Noutput(instance),instance) = vacancy_concentration_ID
             vacancy_generation_output(vacancy_generation_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select

       case ('atol_vacancygeneration')
         vacancy_generation_aTol(instance) = IO_floatValue(line,positions,2_pInt)

       case ('vacancy_frequency')
         vacancy_generation_freq(instance) = IO_floatValue(line,positions,2_pInt)

       case ('vacancy_formationenergy')
         vacancy_generation_formationEnergy(instance) = IO_floatValue(line,positions,2_pInt)

       case ('vacancy_equilibconcentration')
         vacancy_generation_equilibConcentration(instance) = IO_floatValue(line,positions,2_pInt)

       case ('vacancy_diffusionenergy')
         vacancy_generation_diffusionEnergy(instance) = IO_floatValue(line,positions,2_pInt)

       case ('vacancy_diffusioncoeff0')
         vacancy_generation_diffusionCoeff0(instance) = IO_floatValue(line,positions,2_pInt)

       case ('vacancy_stresscoeff')
         vacancy_generation_stressCoeff(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('vacancy_jogheight')
         vacancy_generation_jogHeight(instance) = IO_floatValue(line,positions,2_pInt)

       case ('vacancy_jogseparation')
         vacancy_generation_jogSeparation(instance) = IO_floatValue(line,positions,2_pInt)

       case ('vacancy_nlatticesites')
         vacancy_generation_nLatticeSites(instance) = IO_floatValue(line,positions,2_pInt)

       case ('vacancy_burgersvec')
         vacancy_generation_burgersVec(instance) = IO_floatValue(line,positions,2_pInt)

       case ('pore_surfacefnergy')
         pore_nucleation_surfaceEnergy(instance) = IO_floatValue(line,positions,2_pInt)

       case ('pore_atomvolume')
         pore_nucleation_atomVolume(instance) = IO_floatValue(line,positions,2_pInt)

       case ('pore_shellthickness')
         pore_nucleation_shellThickness(instance) = IO_floatValue(line,positions,2_pInt)

       case ('pore_concentrationcoeff0')
         pore_nucleation_concentrationCoeff0(instance) = IO_floatValue(line,positions,2_pInt)

     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do phase = 1_pInt, size(phase_vacancy)
   if (phase_vacancy(phase) == LOCAL_VACANCY_generation_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_vacancyInstance(phase)

!--------------------------------------------------------------------------------------------------
!  Calculate the coefficient for dislocation motion induced vacancy generation
     vacancy_generation_dislocationCoeff(instance) = vacancy_generation_jogHeight(instance)/     &
                                                     vacancy_generation_jogSeparation(instance)/ &
                                                     vacancy_generation_nLatticeSites(instance)/ &
                                                     vacancy_generation_burgersVec(instance)/    &
                                                     vacancy_generation_burgersVec(instance)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,vacancy_generation_Noutput(instance)
       select case(vacancy_generation_outputID(o,instance))
         case(vacancy_concentration_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          vacancy_generation_sizePostResult(o,instance) = mySize
          vacancy_generation_sizePostResults(instance)  = vacancy_generation_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop
! Determine size of state array
     sizeDotState              =   1_pInt
     sizeState                 =   1_pInt
     vacancyState(phase)%sizeState =       sizeState
     vacancyState(phase)%sizeDotState =    sizeDotState
     vacancyState(phase)%sizePostResults = vacancy_generation_sizePostResults(instance)
     allocate(vacancyState(phase)%aTolState           (sizeState),                source=0.0_pReal)
     allocate(vacancyState(phase)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(vacancyState(phase)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(vacancyState(phase)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(vacancyState(phase)%state               (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(vacancyState(phase)%state_backup        (sizeState,NofMyPhase),     source=0.0_pReal)

     allocate(vacancyState(phase)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(vacancyState(phase)%deltaState          (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(vacancyState(phase)%dotState_backup     (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(vacancyState(phase)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
       allocate(vacancyState(phase)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(vacancyState(phase)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(vacancyState(phase)%RKCK45dotState    (6,sizeDotState,NofMyPhase),source=0.0_pReal)      

     call vacancy_generation_stateInit(phase)
     call vacancy_generation_aTolState(phase,instance)
   endif
 
 enddo initializeInstances
end subroutine vacancy_generation_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant  NEW state values for a given instance of this vacancy model
!--------------------------------------------------------------------------------------------------
subroutine vacancy_generation_stateInit(phase)
 use material, only: &
   vacancyState
 use lattice, only: &
  lattice_equilibriumVacancyConcentration
 
 implicit none
 integer(pInt), intent(in) :: phase                                                    !< number specifying the phase of the vacancy
 real(pReal), dimension(vacancyState(phase)%sizeState) :: tempState

 tempState = lattice_equilibriumVacancyConcentration(phase)
 vacancyState(phase)%state = spread(tempState,2,size(vacancyState(phase)%state(1,:)))
 vacancyState(phase)%state0 = vacancyState(phase)%state
 vacancyState(phase)%partionedState0 = vacancyState(phase)%state
end subroutine vacancy_generation_stateInit

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this vacancy model
!--------------------------------------------------------------------------------------------------
subroutine vacancy_generation_aTolState(phase,instance)
 use material, only: &
  vacancyState

 implicit none
 integer(pInt), intent(in) ::  &
   phase, &
   instance                                                                                         ! number specifying the current instance of the vacancy
 real(pReal), dimension(vacancyState(phase)%sizeState) :: tempTol

 tempTol = vacancy_generation_aTol(instance)
 vacancyState(phase)%aTolState = tempTol
end subroutine vacancy_generation_aTolState
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine vacancy_generation_dotState(nSlip, accumulatedSlip, Tstar_v, Temperature, ipc, ip, el)
 use lattice, only: &
   lattice_massDensity, &
   lattice_specificHeat
 use material, only: &
   mappingConstitutive, &
   phase_vacancyInstance, &
   vacancyState
 use math, only: &
   math_Mandel6to33, &
   math_trace33, &
   pi

 implicit none
 integer(pInt), intent(in) :: &
   nSlip, &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(nSlip), intent(in) :: &
   accumulatedSlip
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal),  intent(in) :: &
   Temperature                                                                                      !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal) :: &
   pressure                                                                                         !< 2nd Piola Kirchhoff stress tensor (Mandel)
 integer(pInt) :: &
   instance, phase, constituent 
 real(pReal) :: &
   vacancyConcentration, &                                                                          !< current vacancy concentration
   vacancyDiffusion, &                                                                              !< the diffusion coefficient D_v
   poleZeldovichCoeff, &                                                                            !< Zeldovich factor of pore nucleation
   vacancyAbsorpRateCoeff, &                                                                        !< vacancy absorption rate
   chemicalPotential, &                                                                             !< the chemical potential due to vacancy concentration
   criticalRadius, &                                                                                !< the critical pore radius
   Gibbs4Pore, &                                                                                    !< the Gibbs free energy for generating a critical pore
   equilibPoreConcentration, &                                                                      !< the equilibrium pore concentration
   nucleationRatePore, &                                                                            !< the nucleation rate of pore
   ratioCvCve                                                                                       !< the ratio of Cv with respect to Cve
 real(pReal) :: &
   threshold4ratioCvCve = 2.0_pReal                                                                 !< the threshold value for Cv/Cve

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_vacancyInstance(phase)
 pressure = math_trace33(math_Mandel6to33(Tstar_v))

!--------------------------------------------------------------------------------------------------
 vacancyConcentration = vacancy_generation_getConcentration(ipc, ip, el)
 ratioCvCve = vacancyConcentration/vacancy_generation_equilibConcentration(instance)

 if(ratioCvCve < threshold4ratioCvCve) then
   nucleationRatePore = 0.0_pReal
 else
!  Calculate nucleation rate of pore
   vacancyDiffusion = vacancy_generation_diffusionCoeff0(instance)* &
                      exp( -vacancy_generation_diffusionEnergy(instance)/(kB*temperature) )
   chemicalPotential = kB*Temperature * log(vacancyConcentration/ &
                       vacancy_generation_equilibConcentration(instance))
   criticalRadius = 2_pReal/chemicalPotential* &
                    pore_nucleation_surfaceEnergy(instance) * pore_nucleation_atomVolume(instance)
   Gibbs4Pore = 4_pReal/3_pReal * pi * pore_nucleation_surfaceEnergy(instance)* &
                criticalRadius * criticalRadius
   equilibPoreConcentration = pore_nucleation_concentrationCoeff0(instance)* &
                              exp( -Gibbs4Pore/(kB*temperature) )

   vacancyAbsorpRateCoeff = 2_pReal/pore_nucleation_shellThickness(instance) * &
                            vacancyDiffusion * vacancyConcentration
   poleZeldovichCoeff = pore_nucleation_atomVolume(instance)* &
                        sqrt( pore_nucleation_surfaceEnergy(instance)/(kB*temperature) )
   nucleationRatePore = poleZeldovichCoeff * vacancyAbsorpRateCoeff* equilibPoreConcentration
 endif

!--------------------------------------------------------------------------------------------------
!  the net generating rate vacancy                            
 vacancyState(phase)%dotState(1,constituent) = &
    vacancy_generation_freq(instance)* &
    exp(-(vacancy_generation_formationEnergy(instance) - vacancy_generation_stressCoeff(instance)*pressure)/ &
         (kB*Temperature)) + &
    sum(accumulatedSlip) * vacancy_generation_dislocationCoeff(instance)- &                         !< Induced by dislocation motion
    nucleationRatePore * (4_pReal/3_pReal * pi * criticalRadius**3_pReal)/ &                        !< Reduced by the formation of pore
    pore_nucleation_atomVolume(instance)   

end subroutine vacancy_generation_dotState

!--------------------------------------------------------------------------------------------------
!> @brief returns vacancy concentration based on state layout 
!--------------------------------------------------------------------------------------------------
function vacancy_generation_getConcentration(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   vacancyState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: vacancy_generation_getConcentration
 
 vacancy_generation_getConcentration = &
   vacancyState(mappingConstitutive(2,ipc,ip,el))%state(1,mappingConstitutive(1,ipc,ip,el))
 
end function vacancy_generation_getConcentration
 
!--------------------------------------------------------------------------------------------------
!> @brief returns temperature based on local damage model state layout 
!--------------------------------------------------------------------------------------------------
subroutine vacancy_generation_putConcentration(ipc, ip, el, localVacancyConcentration)
 use material, only: &
   mappingConstitutive, &
   vacancyState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   localVacancyConcentration
 
 vacancyState(mappingConstitutive(2,ipc,ip,el))%state(1,mappingConstitutive(1,ipc,ip,el))= &
   localVacancyConcentration
 
end subroutine vacancy_generation_putConcentration
 
!--------------------------------------------------------------------------------------------------
!> @brief returns generation vacancy diffusion tensor 
!--------------------------------------------------------------------------------------------------
function vacancy_generation_getVacancyDiffusion33(nSlip,accumulatedSlip,temperature,ipc,ip,el)
 use lattice, only: &
   lattice_VacancyDiffusion33
 use material, only: &
   mappingConstitutive, &
   phase_vacancyInstance, &
   vacancyState

 implicit none
 integer(pInt), intent(in) :: &
   nSlip, &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   vacancy_generation_getVacancyDiffusion33
 real(pReal), dimension(nSlip) :: &
   accumulatedSlip
 real(pReal) :: &
   temperature
 integer(pInt) :: &
   phase, constituent, instance
 
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_vacancyInstance(phase)

 vacancy_generation_getVacancyDiffusion33 = &
   lattice_VacancyDiffusion33(1:3,1:3,phase)* &
   exp(-vacancy_generation_diffusionEnergy(instance)/(kB*temperature))
    
end function vacancy_generation_getVacancyDiffusion33

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function vacancy_generation_postResults(ipc,ip,el)
 use material, only: &
   mappingConstitutive, &
   phase_vacancyInstance, &
   vacancyState

 implicit none
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(vacancy_generation_sizePostResults(phase_vacancyInstance(mappingConstitutive(2,ipc,ip,el)))) :: &
   vacancy_generation_postResults

 integer(pInt) :: &
   instance, phase, constituent, o, c
   
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance  = phase_vacancyInstance(phase)

 c = 0_pInt
 vacancy_generation_postResults = 0.0_pReal

 do o = 1_pInt,vacancy_generation_Noutput(instance)
    select case(vacancy_generation_outputID(o,instance))
 
      case (vacancy_concentration_ID)
        vacancy_generation_postResults(c+1_pInt) = vacancyState(phase)%state(1,constituent)
        c = c + 1
    end select
 enddo
end function vacancy_generation_postResults

end module vacancy_generation
