!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for conservative transport of solute hydrogen
!> @details to be done
!--------------------------------------------------------------------------------------------------
module hydrogenflux_cahnhilliard
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   hydrogenflux_cahnhilliard_sizePostResults                                                           !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   hydrogenflux_cahnhilliard_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   hydrogenflux_cahnhilliard_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   hydrogenflux_cahnhilliard_Noutput                                                                   !< number of outputs per instance of this damage 

 real(pReal),                                                 parameter,           private :: &
   kB = 1.3806488e-23_pReal                                                                          !< Boltzmann constant in J/Kelvin

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 hydrogenConc_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   hydrogenflux_cahnhilliard_outputID                                                                  !< ID of each post result output


 public :: &
   hydrogenflux_cahnhilliard_init, &
   hydrogenflux_cahnhilliard_getMobility33, &
   hydrogenflux_cahnhilliard_getDiffusion33, &
   hydrogenflux_cahnhilliard_getFormationEnergy, &
   hydrogenflux_cahnhilliard_KinematicChemPotAndItsTangent, &
   hydrogenflux_cahnhilliard_getChemPotAndItsTangent, &
   hydrogenflux_cahnhilliard_putHydrogenConcAndItsRate, &
   hydrogenflux_cahnhilliard_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine hydrogenflux_cahnhilliard_init(fileUnit)
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
   hydrogenflux_type, &
   hydrogenflux_typeInstance, &
   homogenization_Noutput, &
   HYDROGENFLUX_cahnhilliard_label, &
   HYDROGENFLUX_cahnhilliard_ID, &
   material_homog, &  
   mappingHomogenization, &
   hydrogenfluxState, &
   hydrogenfluxMapping, &
   hydrogenConc, &
   hydrogenConcRate, &
   hydrogenflux_initialCh, &
   material_partHomogenization, &
   material_partPhase
 use numerics,only: &
   worldrank

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: maxNinstance,mySize=0_pInt,section,instance,o
 integer(pInt) :: sizeState
 integer(pInt) :: NofMyHomog   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  hydrogenflux_'//HYDROGENFLUX_cahnhilliard_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(hydrogenflux_type == HYDROGENFLUX_cahnhilliard_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 allocate(hydrogenflux_cahnhilliard_sizePostResults(maxNinstance),                               source=0_pInt)
 allocate(hydrogenflux_cahnhilliard_sizePostResult (maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(hydrogenflux_cahnhilliard_output         (maxval(homogenization_Noutput),maxNinstance))
          hydrogenflux_cahnhilliard_output = ''
 allocate(hydrogenflux_cahnhilliard_outputID       (maxval(homogenization_Noutput),maxNinstance),source=undefined_ID)
 allocate(hydrogenflux_cahnhilliard_Noutput        (maxNinstance),                               source=0_pInt) 

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

   if (section > 0_pInt ) then; if (hydrogenflux_type(section) == HYDROGENFLUX_cahnhilliard_ID) then  ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran

     instance = hydrogenflux_typeInstance(section)                                                   ! which instance of my hydrogenflux is present homog
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case ('hydrogenconc')
             hydrogenflux_cahnhilliard_Noutput(instance) = hydrogenflux_cahnhilliard_Noutput(instance) + 1_pInt
             hydrogenflux_cahnhilliard_outputID(hydrogenflux_cahnhilliard_Noutput(instance),instance) = hydrogenConc_ID
             hydrogenflux_cahnhilliard_output(hydrogenflux_cahnhilliard_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
          end select

     end select
   endif; endif
 enddo parsingHomog
 
 rewind(fileUnit)
 section = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partPhase)         ! wind forward to <homogenization>
   line = IO_read(fileUnit)
 enddo
 
 initializeInstances: do section = 1_pInt, size(hydrogenflux_type)
   if (hydrogenflux_type(section) == HYDROGENFLUX_cahnhilliard_ID) then
     NofMyHomog=count(material_homog==section)
     instance = hydrogenflux_typeInstance(section)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,hydrogenflux_cahnhilliard_Noutput(instance)
       select case(hydrogenflux_cahnhilliard_outputID(o,instance))
         case(hydrogenConc_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          hydrogenflux_cahnhilliard_sizePostResult(o,instance) = mySize
          hydrogenflux_cahnhilliard_sizePostResults(instance)  = hydrogenflux_cahnhilliard_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop

! allocate state arrays
     sizeState = 0_pInt
     hydrogenfluxState(section)%sizeState = sizeState
     hydrogenfluxState(section)%sizePostResults = hydrogenflux_cahnhilliard_sizePostResults(instance)
     allocate(hydrogenfluxState(section)%state0   (sizeState,NofMyHomog))
     allocate(hydrogenfluxState(section)%subState0(sizeState,NofMyHomog))
     allocate(hydrogenfluxState(section)%state    (sizeState,NofMyHomog))

     nullify(hydrogenfluxMapping(section)%p)
     hydrogenfluxMapping(section)%p => mappingHomogenization(1,:,:)
     deallocate(hydrogenConc    (section)%p)
     deallocate(hydrogenConcRate(section)%p)
     allocate  (hydrogenConc    (section)%p(NofMyHomog), source=hydrogenflux_initialCh(section))
     allocate  (hydrogenConcRate(section)%p(NofMyHomog), source=0.0_pReal)
     
   endif
 
 enddo initializeInstances
 
end subroutine hydrogenflux_cahnhilliard_init

!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized solute mobility tensor in reference configuration
!--------------------------------------------------------------------------------------------------
function hydrogenflux_cahnhilliard_getMobility33(ip,el)
 use lattice, only: &
   lattice_hydrogenfluxMobility33
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
   hydrogenflux_cahnhilliard_getMobility33
 integer(pInt) :: &
   grain
  
 hydrogenflux_cahnhilliard_getMobility33 = 0.0_pReal
 do grain = 1, homogenization_Ngrains(mesh_element(3,el))
   hydrogenflux_cahnhilliard_getMobility33 = hydrogenflux_cahnhilliard_getMobility33 + &
    crystallite_push33ToRef(grain,ip,el,lattice_hydrogenfluxMobility33(:,:,material_phase(grain,ip,el)))
 enddo

 hydrogenflux_cahnhilliard_getMobility33 = &
   hydrogenflux_cahnhilliard_getMobility33/ &
   homogenization_Ngrains(mesh_element(3,el))
 
end function hydrogenflux_cahnhilliard_getMobility33
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized solute nonlocal diffusion tensor in reference configuration
!--------------------------------------------------------------------------------------------------
function hydrogenflux_cahnhilliard_getDiffusion33(ip,el)
 use lattice, only: &
   lattice_hydrogenfluxDiffusion33
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
   hydrogenflux_cahnhilliard_getDiffusion33
 integer(pInt) :: &
   grain
  
 hydrogenflux_cahnhilliard_getDiffusion33 = 0.0_pReal
 do grain = 1, homogenization_Ngrains(mesh_element(3,el))
   hydrogenflux_cahnhilliard_getDiffusion33 = hydrogenflux_cahnhilliard_getDiffusion33 + &
    crystallite_push33ToRef(grain,ip,el,lattice_hydrogenfluxDiffusion33(:,:,material_phase(grain,ip,el)))
 enddo

 hydrogenflux_cahnhilliard_getDiffusion33 = &
   hydrogenflux_cahnhilliard_getDiffusion33/ &
   homogenization_Ngrains(mesh_element(3,el))
 
end function hydrogenflux_cahnhilliard_getDiffusion33
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized solution energy
!--------------------------------------------------------------------------------------------------
function hydrogenflux_cahnhilliard_getFormationEnergy(ip,el)
 use lattice, only: &
   lattice_hydrogenFormationEnergy, &
   lattice_hydrogenVol, &
   lattice_hydrogenSurfaceEnergy
 use material, only: &
   homogenization_Ngrains, &
   material_phase
 use mesh, only: &
   mesh_element

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: &
   hydrogenflux_cahnhilliard_getFormationEnergy
 integer(pInt) :: &
   grain
  
 hydrogenflux_cahnhilliard_getFormationEnergy = 0.0_pReal
 do grain = 1, homogenization_Ngrains(mesh_element(3,el))
   hydrogenflux_cahnhilliard_getFormationEnergy = hydrogenflux_cahnhilliard_getFormationEnergy + &
    lattice_hydrogenFormationEnergy(material_phase(grain,ip,el))/ &
    lattice_hydrogenVol(material_phase(grain,ip,el))/ &
    lattice_hydrogenSurfaceEnergy(material_phase(grain,ip,el))
 enddo

 hydrogenflux_cahnhilliard_getFormationEnergy = &
   hydrogenflux_cahnhilliard_getFormationEnergy/ &
   homogenization_Ngrains(mesh_element(3,el))
 
end function hydrogenflux_cahnhilliard_getFormationEnergy
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized hydrogen entropy coefficient
!--------------------------------------------------------------------------------------------------
function hydrogenflux_cahnhilliard_getEntropicCoeff(ip,el)
 use lattice, only: &
   lattice_hydrogenVol, &
   lattice_hydrogenSurfaceEnergy
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
 real(pReal) :: &
   hydrogenflux_cahnhilliard_getEntropicCoeff
 integer(pInt) :: &
   grain
  
 hydrogenflux_cahnhilliard_getEntropicCoeff = 0.0_pReal
 do grain = 1, homogenization_Ngrains(material_homog(ip,el))
   hydrogenflux_cahnhilliard_getEntropicCoeff = hydrogenflux_cahnhilliard_getEntropicCoeff + &
    kB/ &
    lattice_hydrogenVol(material_phase(grain,ip,el))/ &
    lattice_hydrogenSurfaceEnergy(material_phase(grain,ip,el))
 enddo

 hydrogenflux_cahnhilliard_getEntropicCoeff = &
   hydrogenflux_cahnhilliard_getEntropicCoeff* &
   temperature(material_homog(ip,el))%p(thermalMapping(material_homog(ip,el))%p(ip,el))/ &
   homogenization_Ngrains(material_homog(ip,el))
 
end function hydrogenflux_cahnhilliard_getEntropicCoeff
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized kinematic contribution to chemical potential
!--------------------------------------------------------------------------------------------------
subroutine hydrogenflux_cahnhilliard_KinematicChemPotAndItsTangent(KPot, dKPot_dCh, Ch, ip, el)
 use lattice, only: &
   lattice_hydrogenSurfaceEnergy
 use material, only: &
   homogenization_Ngrains, &
   material_homog, &
   phase_kinematics, &
   phase_Nkinematics, &
   material_phase, &
   KINEMATICS_hydrogen_strain_ID
 use crystallite, only: &
   crystallite_Tstar_v, &
   crystallite_Fi0, &
   crystallite_Fi
 use kinematics_hydrogen_strain, only: &
   kinematics_hydrogen_strain_ChemPotAndItsTangent

 implicit none
 integer(pInt), intent(in)  :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in)  :: &
   Ch
 real(pReal),   intent(out) :: &
   KPot, dKPot_dCh
 real(pReal) :: &
   my_KPot, my_dKPot_dCh
 integer(pInt) :: &
   grain, kinematics
  
 KPot = 0.0_pReal
 dKPot_dCh = 0.0_pReal
 do grain = 1_pInt,homogenization_Ngrains(material_homog(ip,el)) 
   do kinematics = 1_pInt, phase_Nkinematics(material_phase(grain,ip,el))
     select case (phase_kinematics(kinematics,material_phase(grain,ip,el)))
       case (KINEMATICS_hydrogen_strain_ID)
         call kinematics_hydrogen_strain_ChemPotAndItsTangent(my_KPot, my_dKPot_dCh, &
                                                             crystallite_Tstar_v(1:6,grain,ip,el), &
                                                             crystallite_Fi0(1:3,1:3,grain,ip,el), &
                                                             crystallite_Fi (1:3,1:3,grain,ip,el), &
                                                             grain,ip, el)

       case default
         my_KPot = 0.0_pReal
         my_dKPot_dCh = 0.0_pReal
     
     end select
     KPot = KPot + my_KPot/lattice_hydrogenSurfaceEnergy(material_phase(grain,ip,el))
     dKPot_dCh = dKPot_dCh + my_dKPot_dCh/lattice_hydrogenSurfaceEnergy(material_phase(grain,ip,el))
   enddo
 enddo 
 
 KPot = KPot/homogenization_Ngrains(material_homog(ip,el))  
 dKPot_dCh = dKPot_dCh/homogenization_Ngrains(material_homog(ip,el))  

end subroutine hydrogenflux_cahnhilliard_KinematicChemPotAndItsTangent
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized chemical potential
!--------------------------------------------------------------------------------------------------
subroutine hydrogenflux_cahnhilliard_getChemPotAndItsTangent(ChemPot,dChemPot_dCh,Ch,ip,el)
 use numerics, only: &
   hydrogenBoundPenalty, &
   hydrogenPolyOrder

 implicit none
 integer(pInt), intent(in)  :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in)  :: &
   Ch
 real(pReal),   intent(out) :: &
   ChemPot, &
   dChemPot_dCh
 real(pReal) :: &
   kBT, KPot, dKPot_dCh  
 integer(pInt) :: &
   o  
  
 ChemPot = hydrogenflux_cahnhilliard_getFormationEnergy(ip,el)
 dChemPot_dCh = 0.0_pReal
 kBT = hydrogenflux_cahnhilliard_getEntropicCoeff(ip,el)
 do o = 1_pInt, hydrogenPolyOrder
   ChemPot = ChemPot + kBT*((2.0_pReal*Ch - 1.0_pReal)**real(2_pInt*o-1_pInt,pReal))/ &
                       real(2_pInt*o-1_pInt,pReal)
   dChemPot_dCh = dChemPot_dCh + 2.0_pReal*kBT*(2.0_pReal*Ch - 1.0_pReal)**real(2_pInt*o-2_pInt,pReal) 
 enddo
 
 call hydrogenflux_cahnhilliard_KinematicChemPotAndItsTangent(KPot, dKPot_dCh, Ch, ip, el)
 ChemPot = ChemPot + KPot 
 dChemPot_dCh = dChemPot_dCh + dKPot_dCh

 if (Ch < 0.0_pReal) then
   ChemPot = ChemPot - 3.0_pReal*hydrogenBoundPenalty*Ch*Ch
   dChemPot_dCh = dChemPot_dCh - 6.0_pReal*hydrogenBoundPenalty*Ch
 elseif (Ch > 1.0_pReal) then
   ChemPot = ChemPot + 3.0_pReal*hydrogenBoundPenalty*(1.0_pReal - Ch)*(1.0_pReal - Ch)
   dChemPot_dCh = dChemPot_dCh - 6.0_pReal*hydrogenBoundPenalty*(1.0_pReal - Ch)
 endif       
 
end subroutine hydrogenflux_cahnhilliard_getChemPotAndItsTangent
 
!--------------------------------------------------------------------------------------------------
!> @brief updates hydrogen concentration with solution from Cahn-Hilliard PDE for solute transport
!--------------------------------------------------------------------------------------------------
subroutine hydrogenflux_cahnhilliard_putHydrogenConcAndItsRate(Ch,Chdot,ip,el)
 use material, only: &
   mappingHomogenization, &
   hydrogenConc, &
   hydrogenConcRate, &
   hydrogenfluxMapping

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   Ch, &
   Chdot
 integer(pInt) :: &
   homog, &
   offset  
 
 homog  = mappingHomogenization(2,ip,el)
 offset = hydrogenfluxMapping(homog)%p(ip,el)
 hydrogenConc    (homog)%p(offset) = Ch
 hydrogenConcRate(homog)%p(offset) = Chdot

end subroutine hydrogenflux_cahnhilliard_putHydrogenConcAndItsRate
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of hydrogen transport results
!--------------------------------------------------------------------------------------------------
function hydrogenflux_cahnhilliard_postResults(ip,el)
 use material, only: &
   mappingHomogenization, &
   hydrogenflux_typeInstance, &
   hydrogenConc, &
   hydrogenfluxMapping

 implicit none
 integer(pInt),              intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(hydrogenflux_cahnhilliard_sizePostResults(hydrogenflux_typeInstance(mappingHomogenization(2,ip,el)))) :: &
   hydrogenflux_cahnhilliard_postResults

 integer(pInt) :: &
   instance, homog, offset, o, c
   
 homog     = mappingHomogenization(2,ip,el)
 offset    = hydrogenfluxMapping(homog)%p(ip,el)
 instance  = hydrogenflux_typeInstance(homog)

 c = 0_pInt
 hydrogenflux_cahnhilliard_postResults = 0.0_pReal

 do o = 1_pInt,hydrogenflux_cahnhilliard_Noutput(instance)
    select case(hydrogenflux_cahnhilliard_outputID(o,instance))
 
      case (hydrogenConc_ID)
        hydrogenflux_cahnhilliard_postResults(c+1_pInt) = hydrogenConc(homog)%p(offset)
        c = c + 1
    end select
 enddo
end function hydrogenflux_cahnhilliard_postResults

end module hydrogenflux_cahnhilliard
