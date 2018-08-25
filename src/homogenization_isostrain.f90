!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Isostrain (full constraint Taylor assuption) homogenization scheme
!--------------------------------------------------------------------------------------------------
module homogenization_isostrain
 use prec, only: &
   pInt
 
 implicit none
 private
 integer(pInt),               dimension(:),   allocatable,         private :: &
   homogenization_isostrain_Ngrains

 enum, bind(c) 
   enumerator :: parallel_ID, &
                 average_ID
 end enum

 integer(kind(average_ID)),   dimension(:),   allocatable,         private :: &
  homogenization_isostrain_mapping                                                                  !< mapping type


 public :: &
   homogenization_isostrain_init, &
   homogenization_isostrain_partitionDeformation, &
   homogenization_isostrain_averageStressAndItsTangent

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostrain_init()
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   pReal
 use debug, only: &
   debug_HOMOGENIZATION, &
   debug_level, &
   debug_levelBasic
 use IO
 use material
 use config
 
 implicit none
 integer(pInt) :: &
   h
 integer :: &
   maxNinstance, &
   instance
 integer :: &
   NofMyHomog                                                                                       ! no pInt (stores a system dependen value from 'count'
 character(len=65536) :: &
   tag  = ''
 
 write(6,'(/,a)')   ' <<<+-  homogenization_'//HOMOGENIZATION_ISOSTRAIN_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = count(homogenization_type == HOMOGENIZATION_ISOSTRAIN_ID)
 if (maxNinstance == 0) return
 
 if (iand(debug_level(debug_HOMOGENIZATION),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 allocate(homogenization_isostrain_Ngrains(maxNinstance),source=0_pInt)
 allocate(homogenization_isostrain_mapping(maxNinstance),source=average_ID)

 do h = 1_pInt, size(homogenization_type)
   if (homogenization_type(h) /= HOMOGENIZATION_ISOSTRAIN_ID) cycle
   instance = homogenization_typeInstance(h)
  
   homogenization_isostrain_Ngrains(instance) = config_homogenization(h)%getInt('nconstituents')
   tag = 'sum'
   tag = config_homogenization(h)%getString('mapping',defaultVal = tag)
   select case(trim(tag))
     case ('parallel','sum')
       homogenization_isostrain_mapping(instance) = parallel_ID
     case ('average','mean','avg')
       homogenization_isostrain_mapping(instance) = average_ID
     case default
       call IO_error(211_pInt,ext_msg=trim(tag)//' ('//HOMOGENIZATION_isostrain_label//')')
   end select

   NofMyHomog = count(material_homog == h)

   homogState(h)%sizeState       = 0_pInt
   homogState(h)%sizePostResults = 0_pInt
   allocate(homogState(h)%state0   (0_pInt,NofMyHomog), source=0.0_pReal)
   allocate(homogState(h)%subState0(0_pInt,NofMyHomog), source=0.0_pReal)
   allocate(homogState(h)%state    (0_pInt,NofMyHomog), source=0.0_pReal)

 enddo

end subroutine homogenization_isostrain_init


!--------------------------------------------------------------------------------------------------
!> @brief partitions the deformation gradient onto the constituents
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostrain_partitionDeformation(F,avgF,el)
 use prec, only: &
   pReal
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_maxNgrains, &
   homogenization_Ngrains
 
 implicit none
 real(pReal),   dimension (3,3,homogenization_maxNgrains), intent(out) :: F                         !< partioned def grad per grain
 real(pReal),   dimension (3,3),                           intent(in)  :: avgF                      !< my average def grad
 integer(pInt),                                            intent(in)  :: &
   el                                                                                               !< element number
 F = 0.0_pReal
 F(1:3,1:3,1:homogenization_Ngrains(mesh_element(3,el))) = &
                                          spread(avgF,3,homogenization_Ngrains(mesh_element(3,el)))

end subroutine homogenization_isostrain_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities 
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostrain_averageStressAndItsTangent(avgP,dAvgPdAvgF,P,dPdF,el)
 use prec, only: &
   pReal
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_maxNgrains, &
   homogenization_Ngrains, &
   homogenization_typeInstance
 
 implicit none
 real(pReal),   dimension (3,3),                               intent(out) :: avgP                  !< average stress at material point
 real(pReal),   dimension (3,3,3,3),                           intent(out) :: dAvgPdAvgF            !< average stiffness at material point
 real(pReal),   dimension (3,3,homogenization_maxNgrains),     intent(in)  :: P                     !< array of current grain stresses
 real(pReal),   dimension (3,3,3,3,homogenization_maxNgrains), intent(in)  :: dPdF                  !< array of current grain stiffnesses
 integer(pInt),                                                intent(in)  :: el                    !< element number
 integer(pInt) :: &
   homID, & 
   Ngrains

 homID = homogenization_typeInstance(mesh_element(3,el))
 Ngrains = homogenization_Ngrains(mesh_element(3,el))

 select case (homogenization_isostrain_mapping(homID))
   case (parallel_ID)
     avgP       = sum(P,3)
     dAvgPdAvgF = sum(dPdF,5)
   case (average_ID)
     avgP       = sum(P,3)   /real(Ngrains,pReal)
     dAvgPdAvgF = sum(dPdF,5)/real(Ngrains,pReal)
 end select

end subroutine homogenization_isostrain_averageStressAndItsTangent

end module homogenization_isostrain
