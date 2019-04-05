!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Isostrain (full constraint Taylor assuption) homogenization scheme
!--------------------------------------------------------------------------------------------------
module homogenization_isostrain
 use prec, only: &
   pInt

 implicit none
 private
 enum, bind(c) 
   enumerator :: &
     parallel_ID, &
     average_ID
 end enum

 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   integer :: &
     Nconstituents
   integer(kind(average_ID)) :: &
     mapping
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)

 public :: &
   homogenization_isostrain_init, &
   homogenization_isostrain_partitionDeformation, &
   homogenization_isostrain_averageStressAndItsTangent

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostrain_init()
 use debug, only: &
   debug_HOMOGENIZATION, &
   debug_level, &
   debug_levelBasic
 use IO, only: &
   IO_error
 use material, only: &
   homogenization_type, &
   material_homogenizationAt, &
   homogState, &
   HOMOGENIZATION_ISOSTRAIN_ID, &
   HOMOGENIZATION_ISOSTRAIN_LABEL, &
   homogenization_typeInstance
 use config, only: &
   config_homogenization
 
 implicit none
 integer :: &
   Ninstance, &
   h, &
   NofMyHomog
 character(len=65536) :: &
   tag  = ''

 write(6,'(/,a)')   ' <<<+-  homogenization_'//HOMOGENIZATION_ISOSTRAIN_label//' init  -+>>>'

 Ninstance = int(count(homogenization_type == HOMOGENIZATION_ISOSTRAIN_ID),pInt)
 if (iand(debug_level(debug_HOMOGENIZATION),debug_levelBasic) /= 0) &
   write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

 allocate(param(Ninstance))                                                                         ! one container of parameters per instance

 do h = 1, size(homogenization_type)
   if (homogenization_type(h) /= HOMOGENIZATION_ISOSTRAIN_ID) cycle
   
   associate(prm => param(homogenization_typeInstance(h)),&
             config => config_homogenization(h))
  
   prm%Nconstituents = config_homogenization(h)%getInt('nconstituents')
   tag = 'sum'
   select case(trim(config%getString('mapping',defaultVal = tag)))
     case ('sum')
       prm%mapping = parallel_ID
     case ('avg')
       prm%mapping = average_ID
     case default
       call IO_error(211,ext_msg=trim(tag)//' ('//HOMOGENIZATION_isostrain_label//')')
   end select

   NofMyHomog = count(material_homogenizationAt == h)
   homogState(h)%sizeState       = 0
   homogState(h)%sizePostResults = 0
   allocate(homogState(h)%state0   (0,NofMyHomog))
   allocate(homogState(h)%subState0(0,NofMyHomog))
   allocate(homogState(h)%state    (0,NofMyHomog))
   
   end associate

 enddo

end subroutine homogenization_isostrain_init


!--------------------------------------------------------------------------------------------------
!> @brief partitions the deformation gradient onto the constituents
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostrain_partitionDeformation(F,avgF)
 use prec, only: &
   pReal
 
 implicit none
 real(pReal),   dimension (:,:,:), intent(out) :: F                                                 !< partitioned deformation gradient
 
 real(pReal),   dimension (3,3),   intent(in)  :: avgF                                              !< average deformation gradient at material point

 F = spread(avgF,3,size(F,3))

end subroutine homogenization_isostrain_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities 
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostrain_averageStressAndItsTangent(avgP,dAvgPdAvgF,P,dPdF,instance)
 use prec, only: &
   pReal
 
 implicit none
 real(pReal),   dimension (3,3),       intent(out) :: avgP                                          !< average stress at material point
 real(pReal),   dimension (3,3,3,3),   intent(out) :: dAvgPdAvgF                                    !< average stiffness at material point

 real(pReal),   dimension (:,:,:),     intent(in)  :: P                                             !< partitioned stresses
 real(pReal),   dimension (:,:,:,:,:), intent(in)  :: dPdF                                          !< partitioned stiffnesses
 integer,                              intent(in)  :: instance 

 associate(prm => param(instance))
 
 select case (prm%mapping)
   case (parallel_ID)
     avgP       = sum(P,3)
     dAvgPdAvgF = sum(dPdF,5)
   case (average_ID)
     avgP       = sum(P,3)   /real(prm%Nconstituents,pReal)
     dAvgPdAvgF = sum(dPdF,5)/real(prm%Nconstituents,pReal)
 end select
 
 end associate

end subroutine homogenization_isostrain_averageStressAndItsTangent

end module homogenization_isostrain
