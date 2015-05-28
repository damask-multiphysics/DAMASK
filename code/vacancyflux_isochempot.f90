!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for locally evolving vacancy concentration
!> @details to be done
!--------------------------------------------------------------------------------------------------
module vacancyflux_isochempot
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   vacancyflux_isochempot_sizePostResults                                                           !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   vacancyflux_isochempot_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   vacancyflux_isochempot_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   vacancyflux_isochempot_Noutput                                                                   !< number of outputs per instance of this damage 

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 vacancyconc_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   vacancyflux_isochempot_outputID                                                                  !< ID of each post result output


 public :: &
   vacancyflux_isochempot_init, &
   vacancyflux_isochempot_updateState, &
   vacancyflux_isochempot_getSourceAndItsTangent, &
   vacancyflux_isochempot_postResults 

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine vacancyflux_isochempot_init(fileUnit)
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
   VACANCYFLUX_isochempot_label, &
   VACANCYFLUX_isochempot_ID, &
   material_homog, &  
   mappingHomogenization, &
   vacancyfluxState, &
   vacancyfluxMapping, &
   vacancyConc, &
   vacancyConcRate, &
   material_partHomogenization
 use numerics,only: &
   worldrank

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), parameter :: MAXNCHUNKS = 7_pInt
 integer(pInt), dimension(1+2*MAXNCHUNKS) :: positions
 integer(pInt) :: maxNinstance,mySize=0_pInt,section,instance,o
 integer(pInt) :: sizeState
 integer(pInt) :: NofMyHomog   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  vacancyflux_'//VACANCYFLUX_isochempot_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(vacancyflux_type == VACANCYFLUX_isochempot_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 allocate(vacancyflux_isochempot_sizePostResults(maxNinstance),                               source=0_pInt)
 allocate(vacancyflux_isochempot_sizePostResult (maxval(homogenization_Noutput),maxNinstance),source=0_pInt)
 allocate(vacancyflux_isochempot_output         (maxval(homogenization_Noutput),maxNinstance))
          vacancyflux_isochempot_output = ''
 allocate(vacancyflux_isochempot_outputID       (maxval(homogenization_Noutput),maxNinstance),source=undefined_ID)
 allocate(vacancyflux_isochempot_Noutput        (maxNinstance),                               source=0_pInt) 

 rewind(fileUnit)
 section = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partHomogenization)! wind forward to <homogenization>
   line = IO_read(fileUnit)
 enddo
 
 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of homog part
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

   if (section > 0_pInt ) then; if (vacancyflux_type(section) == VACANCYFLUX_isochempot_ID) then             ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran

     instance = vacancyflux_typeInstance(section)                                                       ! which instance of my vacancyflux is present homog
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('vacancyconc')
             vacancyflux_isochempot_Noutput(instance) = vacancyflux_isochempot_Noutput(instance) + 1_pInt
             vacancyflux_isochempot_outputID(vacancyflux_isochempot_Noutput(instance),instance) = vacancyconc_ID
             vacancyflux_isochempot_output(vacancyflux_isochempot_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select

     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do section = 1_pInt, size(vacancyflux_type)
   if (vacancyflux_type(section) == VACANCYFLUX_isochempot_ID) then
     NofMyHomog=count(material_homog==section)
     instance = vacancyflux_typeInstance(section)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,vacancyflux_isochempot_Noutput(instance)
       select case(vacancyflux_isochempot_outputID(o,instance))
         case(vacancyconc_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          vacancyflux_isochempot_sizePostResult(o,instance) = mySize
          vacancyflux_isochempot_sizePostResults(instance)  = vacancyflux_isochempot_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop

! allocate state arrays
     sizeState = 1_pInt
     vacancyfluxState(section)%sizeState = sizeState
     vacancyfluxState(section)%sizePostResults = vacancyflux_isochempot_sizePostResults(instance)
     allocate(vacancyfluxState(section)%state0   (sizeState,NofMyHomog), source=0.0_pReal)
     allocate(vacancyfluxState(section)%subState0(sizeState,NofMyHomog), source=0.0_pReal)
     allocate(vacancyfluxState(section)%state    (sizeState,NofMyHomog), source=0.0_pReal)

     nullify(vacancyfluxMapping(section)%p)
     vacancyfluxMapping(section)%p => mappingHomogenization(1,:,:)
     deallocate(vacancyConc(section)%p)
     vacancyConc(section)%p => vacancyfluxState(section)%state(1,:)
     deallocate(vacancyConcRate(section)%p)
     allocate(vacancyConcRate(section)%p(NofMyHomog), source=0.0_pReal)
     
   endif
 
 enddo initializeInstances
end subroutine vacancyflux_isochempot_init

!--------------------------------------------------------------------------------------------------
!> @brief  calculates change in vacancy concentration based on local vacancy generation model  
!--------------------------------------------------------------------------------------------------
function vacancyflux_isochempot_updateState(subdt, ip, el)
 use numerics, only: &
   err_vacancyflux_tolAbs, &
   err_vacancyflux_tolRel
 use material, only: &
   mappingHomogenization, &
   vacancyflux_typeInstance, &
   vacancyfluxState, &
   vacancyConc, &
   vacancyConcRate, &
   vacancyfluxMapping

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   subdt
 logical,                    dimension(2)                             :: &
   vacancyflux_isochempot_updateState
 integer(pInt) :: &
   homog, &
   offset, &
   instance
 real(pReal) :: &
   Cv, Cvdot, dCvDot_dCv  

 homog  = mappingHomogenization(2,ip,el)
 offset = mappingHomogenization(1,ip,el)
 instance = vacancyflux_typeInstance(homog)
 
 Cv = vacancyfluxState(homog)%subState0(1,offset)
 call vacancyflux_isochempot_getSourceAndItsTangent(CvDot, dCvDot_dCv, Cv, ip, el)
 Cv = Cv + subdt*Cvdot 
 
 vacancyflux_isochempot_updateState = [     abs(Cv - vacancyfluxState(homog)%state(1,offset)) &
                                         <= err_vacancyflux_tolAbs &
                                       .or. abs(Cv - vacancyfluxState(homog)%state(1,offset)) &
                                         <= err_vacancyflux_tolRel*abs(vacancyfluxState(homog)%state(1,offset)), &
                                      .true.]

 vacancyConc    (homog)%p(vacancyfluxMapping(homog)%p(ip,el)) = Cv
 vacancyConcRate(homog)%p(vacancyfluxMapping(homog)%p(ip,el)) = &
   (vacancyfluxState(homog)%state(1,offset) - vacancyfluxState(homog)%subState0(1,offset))/subdt
 
end function vacancyflux_isochempot_updateState

!--------------------------------------------------------------------------------------------------
!> @brief  calculates homogenized vacancy driving forces  
!--------------------------------------------------------------------------------------------------
subroutine vacancyflux_isochempot_getSourceAndItsTangent(CvDot, dCvDot_dCv, Cv, ip, el)
 use material, only: &
   homogenization_Ngrains, &
   mappingHomogenization, &
   mappingConstitutive, &
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
   phase = mappingConstitutive(2,grain,ip,el)
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
 
end subroutine vacancyflux_isochempot_getSourceAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief return array of vacancy transport results
!--------------------------------------------------------------------------------------------------
function vacancyflux_isochempot_postResults(ip,el)
 use material, only: &
   mappingHomogenization, &
   vacancyflux_typeInstance, &
   vacancyConc, &
   vacancyfluxMapping

 implicit none
 integer(pInt),              intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(vacancyflux_isochempot_sizePostResults(vacancyflux_typeInstance(mappingHomogenization(2,ip,el)))) :: &
   vacancyflux_isochempot_postResults

 integer(pInt) :: &
   instance, homog, offset, o, c
   
 homog     = mappingHomogenization(2,ip,el)
 offset    = vacancyfluxMapping(homog)%p(ip,el)
 instance  = vacancyflux_typeInstance(homog)

 c = 0_pInt
 vacancyflux_isochempot_postResults = 0.0_pReal

 do o = 1_pInt,vacancyflux_isochempot_Noutput(instance)
    select case(vacancyflux_isochempot_outputID(o,instance))
 
      case (vacancyconc_ID)
        vacancyflux_isochempot_postResults(c+1_pInt) = vacancyConc(homog)%p(offset)
        c = c + 1
    end select
 enddo
end function vacancyflux_isochempot_postResults

end module vacancyflux_isochempot
