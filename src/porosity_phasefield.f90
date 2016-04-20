!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for phase field modelling of pore nucleation and growth
!> @details phase field model for pore nucleation and growth based on vacancy clustering
!--------------------------------------------------------------------------------------------------
module porosity_phasefield
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   porosity_phasefield_sizePostResults                                                           !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   porosity_phasefield_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   porosity_phasefield_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   porosity_phasefield_Noutput                                                                   !< number of outputs per instance of this porosity 

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 porosity_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   porosity_phasefield_outputID                                                                  !< ID of each post result output


 public :: &
   porosity_phasefield_init, &
   porosity_phasefield_getFormationEnergy, &
   porosity_phasefield_getSurfaceEnergy, &
   porosity_phasefield_getSourceAndItsTangent, &
   porosity_phasefield_getDiffusion33, &
   porosity_phasefield_getMobility, &
   porosity_phasefield_putPorosity, &
   porosity_phasefield_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine porosity_phasefield_init(fileUnit)
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
   porosity_type, &
   porosity_typeInstance, &
   homogenization_Noutput, &
   POROSITY_phasefield_label, &
   POROSITY_phasefield_ID, &
   material_homog, & 
   mappingHomogenization, & 
   porosityState, &
   porosityMapping, &
   porosity, &
   porosity_initialPhi, &
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
   write(6,'(/,a)')   ' <<<+-  porosity_'//POROSITY_phasefield_label//' init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(porosity_type == POROSITY_phasefield_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 allocate(porosity_phasefield_sizePostResults(maxNinstance),                               source=0_pInt)
 allocate(porosity_phasefield_sizePostResult (maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(porosity_phasefield_output         (maxval(homogenization_Noutput),maxNinstance))
          porosity_phasefield_output = ''
 allocate(porosity_phasefield_outputID       (maxval(homogenization_Noutput),maxNinstance),source=undefined_ID)
 allocate(porosity_phasefield_Noutput        (maxNinstance),                               source=0_pInt) 

 rewind(fileUnit)
 section = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partHomogenization)! wind forward to <homogenization>
   line = IO_read(fileUnit)
 enddo
 
 parsingHomog: do while (trim(line) /= IO_EOF)                                                       ! read through sections of homog part
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

   if (section > 0_pInt ) then; if (porosity_type(section) == POROSITY_phasefield_ID) then          ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran

     instance = porosity_typeInstance(section)                                                      ! which instance of my porosity is present homog
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case ('porosity')
             porosity_phasefield_Noutput(instance) = porosity_phasefield_Noutput(instance) + 1_pInt
             porosity_phasefield_outputID(porosity_phasefield_Noutput(instance),instance) = porosity_ID
             porosity_phasefield_output(porosity_phasefield_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
          end select

     end select
   endif; endif
 enddo parsingHomog
 
 initializeInstances: do section = 1_pInt, size(porosity_type)
   if (porosity_type(section) == POROSITY_phasefield_ID) then
     NofMyHomog=count(material_homog==section)
     instance = porosity_typeInstance(section)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,porosity_phasefield_Noutput(instance)
       select case(porosity_phasefield_outputID(o,instance))
         case(porosity_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          porosity_phasefield_sizePostResult(o,instance) = mySize
          porosity_phasefield_sizePostResults(instance)  = porosity_phasefield_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop

! allocate state arrays
     sizeState = 0_pInt
     porosityState(section)%sizeState = sizeState
     porosityState(section)%sizePostResults = porosity_phasefield_sizePostResults(instance)
     allocate(porosityState(section)%state0   (sizeState,NofMyHomog))
     allocate(porosityState(section)%subState0(sizeState,NofMyHomog))
     allocate(porosityState(section)%state    (sizeState,NofMyHomog))

     nullify(porosityMapping(section)%p)
     porosityMapping(section)%p => mappingHomogenization(1,:,:)
     deallocate(porosity(section)%p)
     allocate(porosity(section)%p(NofMyHomog), source=porosity_initialPhi(section))
     
   endif
 
 enddo initializeInstances
end subroutine porosity_phasefield_init

!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized vacancy formation energy
!--------------------------------------------------------------------------------------------------
function porosity_phasefield_getFormationEnergy(ip,el)
 use lattice, only: &
   lattice_vacancyFormationEnergy, &
   lattice_vacancyVol
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
   porosity_phasefield_getFormationEnergy
 integer(pInt) :: &
   grain
  
 porosity_phasefield_getFormationEnergy = 0.0_pReal
 do grain = 1, homogenization_Ngrains(mesh_element(3,el))
   porosity_phasefield_getFormationEnergy = porosity_phasefield_getFormationEnergy + &
    lattice_vacancyFormationEnergy(material_phase(grain,ip,el))/ &
    lattice_vacancyVol(material_phase(grain,ip,el))
 enddo

 porosity_phasefield_getFormationEnergy = &
   porosity_phasefield_getFormationEnergy/ &
   homogenization_Ngrains(mesh_element(3,el))
 
end function porosity_phasefield_getFormationEnergy
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized pore surface energy (normalized by characteristic length)
!--------------------------------------------------------------------------------------------------
function porosity_phasefield_getSurfaceEnergy(ip,el)
 use lattice, only: &
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
 real(pReal) :: &
   porosity_phasefield_getSurfaceEnergy
 integer(pInt) :: &
   grain
  
 porosity_phasefield_getSurfaceEnergy = 0.0_pReal
 do grain = 1, homogenization_Ngrains(mesh_element(3,el))
   porosity_phasefield_getSurfaceEnergy = porosity_phasefield_getSurfaceEnergy + &
    lattice_vacancySurfaceEnergy(material_phase(grain,ip,el))
 enddo

 porosity_phasefield_getSurfaceEnergy = &
   porosity_phasefield_getSurfaceEnergy/ &
   homogenization_Ngrains(mesh_element(3,el))
 
end function porosity_phasefield_getSurfaceEnergy
 
!--------------------------------------------------------------------------------------------------
!> @brief  calculates homogenized local driving force for pore nucleation and growth  
!--------------------------------------------------------------------------------------------------
subroutine porosity_phasefield_getSourceAndItsTangent(phiDot, dPhiDot_dPhi, phi, ip, el)
 use math, only : &
   math_mul33x33, &
   math_mul66x6, &
   math_Mandel33to6, &
   math_transpose33, &
   math_I3
 use material, only: &
   homogenization_Ngrains, &
   material_homog, &
   material_phase, &
   phase_NstiffnessDegradations, &
   phase_stiffnessDegradation, &
   vacancyConc, &
   vacancyfluxMapping, &
   damage, &
   damageMapping, &
   STIFFNESS_DEGRADATION_damage_ID
 use crystallite, only: &
   crystallite_Fe
 use constitutive, only: &
   constitutive_homogenizedC  
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   phi
 integer(pInt) :: &
   phase, &
   grain, &
   homog, &
   mech
 real(pReal) :: &
   phiDot, dPhiDot_dPhi, Cv, W_e, strain(6), C(6,6)  

 homog = material_homog(ip,el)
 Cv = vacancyConc(homog)%p(vacancyfluxMapping(homog)%p(ip,el))

 W_e = 0.0_pReal
 do grain = 1, homogenization_Ngrains(homog)
   phase = material_phase(grain,ip,el)
   strain = math_Mandel33to6(math_mul33x33(math_transpose33(crystallite_Fe(1:3,1:3,grain,ip,el)), &
                                           crystallite_Fe(1:3,1:3,grain,ip,el)) - math_I3)/2.0_pReal
   C = constitutive_homogenizedC(grain,ip,el)
   do mech = 1_pInt, phase_NstiffnessDegradations(phase)
     select case(phase_stiffnessDegradation(mech,phase))
       case (STIFFNESS_DEGRADATION_damage_ID)
         C = damage(homog)%p(damageMapping(homog)%p(ip,el))* &
             damage(homog)%p(damageMapping(homog)%p(ip,el))* &
             C                                        
     
     end select
   enddo                                                 
   W_e = W_e + sum(abs(strain*math_mul66x6(C,strain)))
 enddo
 W_e = W_e/homogenization_Ngrains(homog)
 
 phiDot = 2.0_pReal*(1.0_pReal - phi)*(1.0_pReal - Cv)*(1.0_pReal - Cv) - &
          2.0_pReal*phi*(W_e + Cv*porosity_phasefield_getFormationEnergy(ip,el))/ &
                                  porosity_phasefield_getSurfaceEnergy  (ip,el)
 dPhiDot_dPhi = - 2.0_pReal*(1.0_pReal - Cv)*(1.0_pReal - Cv) &
                - 2.0_pReal*(W_e + Cv*porosity_phasefield_getFormationEnergy(ip,el))/ &
                                      porosity_phasefield_getSurfaceEnergy  (ip,el)

end subroutine porosity_phasefield_getSourceAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized nonlocal diffusion tensor in reference configuration
!--------------------------------------------------------------------------------------------------
function porosity_phasefield_getDiffusion33(ip,el)
 use lattice, only: &
   lattice_PorosityDiffusion33
 use material, only: &
   homogenization_Ngrains, &
   material_phase, &
   mappingHomogenization
 use crystallite, only: &
   crystallite_push33ToRef

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   porosity_phasefield_getDiffusion33
 integer(pInt) :: &
   homog, &
   grain
   
 homog  = mappingHomogenization(2,ip,el)
 porosity_phasefield_getDiffusion33 = 0.0_pReal  
 do grain = 1, homogenization_Ngrains(homog)
   porosity_phasefield_getDiffusion33 = porosity_phasefield_getDiffusion33 + &
     crystallite_push33ToRef(grain,ip,el,lattice_PorosityDiffusion33(1:3,1:3,material_phase(grain,ip,el)))
 enddo

 porosity_phasefield_getDiffusion33 = &
   porosity_phasefield_getDiffusion33/ &
   homogenization_Ngrains(homog)
 
end function porosity_phasefield_getDiffusion33
 
!--------------------------------------------------------------------------------------------------
!> @brief Returns homogenized phase field mobility 
!--------------------------------------------------------------------------------------------------
real(pReal) function porosity_phasefield_getMobility(ip,el)
 use mesh, only: &
   mesh_element
 use lattice, only: &
   lattice_PorosityMobility
 use material, only: &
   material_phase, &
   homogenization_Ngrains

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 integer(pInt) :: &
   ipc
 
 porosity_phasefield_getMobility = 0.0_pReal
                                                
 do ipc = 1, homogenization_Ngrains(mesh_element(3,el))
   porosity_phasefield_getMobility = porosity_phasefield_getMobility + lattice_PorosityMobility(material_phase(ipc,ip,el))
 enddo

 porosity_phasefield_getMobility = porosity_phasefield_getMobility/homogenization_Ngrains(mesh_element(3,el))

end function porosity_phasefield_getMobility

!--------------------------------------------------------------------------------------------------
!> @brief updates porosity with solution from phasefield PDE
!--------------------------------------------------------------------------------------------------
subroutine porosity_phasefield_putPorosity(phi,ip,el)
 use material, only: &
   material_homog, &
   porosityMapping, &
   porosity

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   phi
 integer(pInt) :: &
   homog, &
   offset
 
 homog  = material_homog(ip,el)
 offset = porosityMapping(homog)%p(ip,el)
 porosity(homog)%p(offset) = phi

end subroutine porosity_phasefield_putPorosity
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of porosity results
!--------------------------------------------------------------------------------------------------
function porosity_phasefield_postResults(ip,el)
 use material, only: &
   mappingHomogenization, &
   porosity_typeInstance, &
   porosity

 implicit none
 integer(pInt),              intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(porosity_phasefield_sizePostResults(porosity_typeInstance(mappingHomogenization(2,ip,el)))) :: &
   porosity_phasefield_postResults

 integer(pInt) :: &
   instance, homog, offset, o, c
   
 homog     = mappingHomogenization(2,ip,el)
 offset    = mappingHomogenization(1,ip,el)
 instance  = porosity_typeInstance(homog)

 c = 0_pInt
 porosity_phasefield_postResults = 0.0_pReal

 do o = 1_pInt,porosity_phasefield_Noutput(instance)
    select case(porosity_phasefield_outputID(o,instance))
 
      case (porosity_ID)
        porosity_phasefield_postResults(c+1_pInt) = porosity(homog)%p(offset)
        c = c + 1
    end select
 enddo
end function porosity_phasefield_postResults

end module porosity_phasefield
