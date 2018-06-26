!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Reads in the material configuration from file
!> @details Reads the material configuration file, where solverJobName.materialConfig takes
!! precedence over material.config. Stores the raw strings and the positions of delimiters for the
!! parts 'homogenization', 'crystallite', 'phase', 'texture', and 'microstucture'
!--------------------------------------------------------------------------------------------------
module config
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private 
  type, private :: tPartitionedString
   character(len=:),            allocatable :: val
   integer(pInt), dimension(:), allocatable :: pos
 end type tPartitionedString
 
 type, public :: tPartitionedStringList
   type(tPartitionedString)               :: string
   type(tPartitionedStringList),  pointer :: next => null()

   contains
     procedure :: add            => add
     procedure :: show           => show
     procedure :: free           => free

     procedure :: keyExists      => keyExists
     procedure :: countKeys      => countKeys

     procedure :: getFloat       => getFloat
     procedure :: getInt         => getInt
     procedure :: getString      => getString

     procedure :: getFloats      => getFloats
     procedure :: getInts        => getInts
     procedure :: getStrings     => getStrings

 end type tPartitionedStringList

 type(tPartitionedStringList), public :: emptyList

 type(tPartitionedStringList), public, protected, allocatable, dimension(:) :: &                    ! QUESTION: rename to config_XXX?
   phaseConfig, &
   microstructureConfig, &
   homogenizationConfig, &
   textureConfig, &
   crystalliteConfig
 
 character(len=64), dimension(:), allocatable, public, protected :: &
   phase_name, &                                                                                    !< name of each phase
   homogenization_name, &                                                                           !< name of each homogenization
   crystallite_name, &                                                                              !< name of each crystallite setting
   microstructure_name, &                                                                           !< name of each microstructure
   texture_name                                                                                     !< name of each texture

! ToDo: make private, no one needs to know that
 character(len=*), parameter, public  :: &
   MATERIAL_partHomogenization = 'homogenization', &                                                !< keyword for homogenization part
   MATERIAL_partCrystallite    = 'crystallite', &                                                   !< keyword for crystallite part
   MATERIAL_partPhase          = 'phase', &                                                         !< keyword for phase part
   MATERIAL_partMicrostructure = 'microstructure', &                                                !< keyword for microstructure part
   MATERIAL_partTexture        = 'texture'                                                          !< keyword for texture part

! ToDo: Remove, use size(phaseConfig) etc
 integer(pInt), public, protected :: &
   material_Ntexture, &                                                                             !< number of textures
   material_Nphase, &                                                                               !< number of phases
   material_Nhomogenization, &                                                                      !< number of homogenizations
   material_Nmicrostructure, &                                                                      !< number of microstructures
   material_Ncrystallite                                                                            !< number of crystallite settings

! ToDo: make private, no one needs to know that
 character(len=*), parameter, public  :: &
   MATERIAL_configFile         = 'material.config', &                                               !< generic name for material configuration file
   MATERIAL_localFileExt       = 'materialConfig'                                                   !< extension of solver job name depending material configuration file


public :: config_init

contains

subroutine config_init()
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use IO, only: &
   IO_error, &
   IO_open_file, &
   IO_read, &
   IO_lc, &
   IO_open_jobFile_stat, &
   IO_getTag, &
   IO_timeStamp, &
   IO_EOF
 use debug, only: &
   debug_level, &
   debug_material, &
   debug_levelBasic

 implicit none
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt)            :: myDebug

 character(len=65536) :: &                                                                          
  line, &
  part


 myDebug = debug_level(debug_material)

 write(6,'(/,a)') ' <<<+-  config init  -+>>>'
 write(6,'(a15,a)')   ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ...open material.config file

 rewind(fileUnit)
 line        = ''                                                                                   ! to have it initialized
 do while (trim(line) /= IO_EOF)
   part = IO_lc(IO_getTag(line,'<','>'))

   select case (trim(part))
    
     case (trim(material_partPhase))
       call parseFile(line,phase_name,phaseConfig,FILEUNIT)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Phase parsed'; flush(6)
    
     case (trim(material_partMicrostructure))
       call parseFile(line,microstructure_name,microstructureConfig,FILEUNIT)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Microstructure parsed'; flush(6)
    
     case (trim(material_partCrystallite))
       call parseFile(line,crystallite_name,crystalliteConfig,FILEUNIT)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Crystallite parsed'; flush(6)
    
     case (trim(material_partHomogenization))
       call parseFile(line,homogenization_name,homogenizationConfig,FILEUNIT)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Homogenization parsed'; flush(6)
    
     case (trim(material_partTexture))
       call parseFile(line,texture_name,textureConfig,FILEUNIT)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Texture parsed'; flush(6)

     case default
       line = IO_read(fileUnit)

   end select

 enddo

 material_Nhomogenization = size(homogenizationConfig)
 if (material_Nhomogenization < 1_pInt) call IO_error(160_pInt,ext_msg=material_partHomogenization)
 material_Nmicrostructure = size(microstructureConfig)
 if (material_Nmicrostructure < 1_pInt) call IO_error(160_pInt,ext_msg=material_partMicrostructure)
 material_Ncrystallite = size(crystalliteConfig)
 if (material_Ncrystallite < 1_pInt) call IO_error(160_pInt,ext_msg=material_partCrystallite)
 material_Nphase = size(phaseConfig)
 if (material_Nphase < 1_pInt) call IO_error(160_pInt,ext_msg=material_partPhase)
 material_Ntexture = size(textureConfig)
 if (material_Ntexture < 1_pInt) call IO_error(160_pInt,ext_msg=material_partTexture)


end subroutine config_init


!--------------------------------------------------------------------------------------------------
!> @brief parses the material.config file
!--------------------------------------------------------------------------------------------------
subroutine parseFile(line,&
                     sectionNames,part,fileUnit)
 use IO, only: &
   IO_read, &
   IO_error, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringValue, &
   IO_stringPos, &
   IO_EOF

 implicit none
 integer(pInt),    intent(in) :: fileUnit
 character(len=*),  dimension(:), allocatable, intent(inout)  :: sectionNames
 type(tPartitionedStringList), allocatable, dimension(:), intent(inout) :: part
 character(len=65536),intent(out) :: line

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt)        :: s
 character(len=65536) :: devNull
 character(len=64)    :: tag
 logical              :: echo

 echo = .false. 
 allocate(part(0))

 s = 0_pInt
 do while (trim(line) /= IO_EOF)                                                                    ! read through sections of material part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   foundNextPart: if (IO_getTag(line,'<','>') /= '') then
     devNull = IO_read(fileUnit, .true.)                                                            ! reset IO_read to close any recursively included files
     exit
   endif foundNextPart
   nextSection: if (IO_getTag(line,'[',']') /= '') then
     s = s + 1_pInt
     part = [part, emptyList]
     tag = IO_getTag(line,'[',']')
     GfortranBug86033: if (.not. allocated(sectionNames)) then
       allocate(sectionNames(1),source=tag)
     else GfortranBug86033
       sectionNames  = [sectionNames,tag]
     endif GfortranBug86033
     cycle
   endif nextSection
   chunkPos = IO_stringPos(line)
   tag = IO_lc(IO_stringValue(trim(line),chunkPos,1_pInt))                                          ! extract key
   inSection: if (s > 0_pInt) then
     call part(s)%add(IO_lc(trim(line)))
   else inSection
     echo = (trim(tag) == '/echo/')
   endif inSection
 enddo

 if (echo) then
   do s = 1, size(sectionNames)
     call part(s)%show()
   end do
 end if 

end subroutine parseFile

!--------------------------------------------------------------------------------------------------
!> @brief add element
!> @details Adds a string together with the start/end position of chunks in this string. The new 
!! element is added at the end of the list. Empty strings are not added. All strings are converted
!! to lower case
!--------------------------------------------------------------------------------------------------
subroutine add(this,string)
  use IO, only: &
    IO_isBlank, &
    IO_lc, &
    IO_stringPos

  implicit none
  class(tPartitionedStringList),  target, intent(in) :: this
  character(len=*),                       intent(in) :: string
  type(tPartitionedStringList),   pointer            :: new, item

  if (IO_isBlank(string)) return

  allocate(new)
  new%string%val = IO_lc       (trim(string))
  new%string%pos = IO_stringPos(trim(string))

  item => this
  do while (associated(item%next))
    item => item%next
  enddo
  item%next => new

end subroutine add


!--------------------------------------------------------------------------------------------------
!> @brief prints all elements
!> @details Strings are printed in order of insertion (FIFO)
!--------------------------------------------------------------------------------------------------
subroutine show(this)

 implicit none
 class(tPartitionedStringList) :: this
 type(tPartitionedStringList), pointer :: item

 item => this%next
 do while (associated(item))
   write(6,'(a)') trim(item%string%val)
   item => item%next
 end do

end subroutine show


!--------------------------------------------------------------------------------------------------
!> @brief cleans entire list
!> @details list head is remains alive
!--------------------------------------------------------------------------------------------------
subroutine free(this)

  implicit none
  class(tPartitionedStringList),  target, intent(in) :: this
  type(tPartitionedStringList),   pointer            :: new, item

  item => this%next
  do while (associated(item%next))
    new => item
    deallocate(item)
    item => new%next
  enddo
  deallocate(item)

end subroutine free


!--------------------------------------------------------------------------------------------------
!> @brief reports wether a given key (string value at first position) exists in the list
!--------------------------------------------------------------------------------------------------
logical function keyExists(this,key)
 use IO, only: &
   IO_stringValue

 implicit none
 class(tPartitionedStringList), intent(in) :: this
 character(len=*), intent(in)              :: key
 type(tPartitionedStringList), pointer     :: item

 keyExists = .false.

 item => this%next
 do while (associated(item) .and. .not. keyExists)
   keyExists = trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)
   item => item%next
 end do

end function keyExists


!--------------------------------------------------------------------------------------------------
!> @brief count number of key appearances
!> @details traverses list and counts each occurrence of specified key
!--------------------------------------------------------------------------------------------------
integer(pInt) function countKeys(this,key)
 use IO, only: &
   IO_stringValue

 implicit none

 class(tPartitionedStringList), intent(in) :: this
 character(len=*), intent(in)              :: key
 type(tPartitionedStringList), pointer     :: item

 countKeys = 0_pInt

 item => this%next
 do while (associated(item))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) &
     countKeys = countKeys + 1_pInt
   item => item%next
 end do

end function countKeys


!--------------------------------------------------------------------------------------------------
!> @brief gets float value of for a given key from a linked list
!> @details gets the last value if the key occurs more than once. If key is not found exits with 
!! error unless default is given
!--------------------------------------------------------------------------------------------------
real(pReal) function getFloat(this,key,defaultVal)
 use IO, only : &
   IO_error, &
   IO_stringValue, &
   IO_FloatValue

 implicit none
 class(tPartitionedStringList), intent(in)           :: this
 character(len=*),              intent(in)           :: key
 real(pReal),                   intent(in), optional :: defaultVal
 type(tPartitionedStringList),  pointer              :: item
 logical                                             :: found

 found = present(defaultVal)
 if (found) getFloat = defaultVal
 
 item => this%next
 do while (associated(item))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (item%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     getFloat = IO_FloatValue(item%string%val,item%string%pos,2)
   endif
   item => item%next
 end do

 if (.not. found) call IO_error(140_pInt,ext_msg=key)

end function getFloat


!--------------------------------------------------------------------------------------------------
!> @brief gets integer value of for a given key from a linked list
!> @details gets the last value if the key occurs more than once. If key is not found exits with 
!! error unless default is given
!--------------------------------------------------------------------------------------------------
integer(pInt) function getInt(this,key,defaultVal)
 use IO, only: &
   IO_error, &
   IO_stringValue, &
   IO_IntValue

 implicit none
 class(tPartitionedStringList), intent(in)           :: this
 character(len=*),              intent(in)           :: key
 integer(pInt),                 intent(in), optional :: defaultVal
 type(tPartitionedStringList),  pointer              :: item
 logical                                             :: found

 found = present(defaultVal)
 if (found) getInt = defaultVal
 
 item => this%next
 do while (associated(item))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (item%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     getInt = IO_IntValue(item%string%val,item%string%pos,2)
   endif
   item => item%next
 end do

 if (.not. found) call IO_error(140_pInt,ext_msg=key)

end function getInt


!--------------------------------------------------------------------------------------------------
!> @brief gets string value of for a given key from a linked list
!> @details gets the last value if the key occurs more than once. If key is not found exits with 
!! error unless default is given. If raw is true, the the complete string is returned, otherwise 
!! the individual chunks are returned
!--------------------------------------------------------------------------------------------------
character(len=65536) function getString(this,key,defaultVal,raw)
 use IO, only: &
   IO_error, &
   IO_stringValue

 implicit none
 class(tPartitionedStringList), intent(in)           :: this
 character(len=*),              intent(in)           :: key
 character(len=65536),          intent(in), optional :: defaultVal
 logical,                       intent(in), optional :: raw
 type(tPartitionedStringList),  pointer              :: item
 logical                                             :: found, &
                                                        whole

 whole = merge(raw,.false.,present(raw))                                                            ! whole string or white space splitting
 found = present(defaultVal)
 if (found) getString = trim(defaultVal)

 item => this%next
 do while (associated(item))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (item%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)

     if (whole) then
       getString = trim(item%string%val(item%string%pos(4):))                                       ! raw string starting a second chunk
     else
       getString = IO_StringValue(item%string%val,item%string%pos,2)
     endif
   endif
   item => item%next
 end do

 if (.not. found) call IO_error(140_pInt,ext_msg=key)
 if (present(defaultVal)) then
   if (len_trim(getString) /= len_trim(defaultVal)) call IO_error(0_pInt,ext_msg='getString')
 endif
end function getString


!--------------------------------------------------------------------------------------------------
!> @brief gets array of float values of for a given key from a linked list
!> @details for cumulative keys, "()", values from all occurrences are return. Otherwise only all
!! values from the last occurrence. If key is not found exits with error unless default is given.
!--------------------------------------------------------------------------------------------------
function getFloats(this,key,defaultVal)
 use IO, only: &
   IO_error, &
   IO_stringValue, &
   IO_FloatValue

 implicit none
 real(pReal),     dimension(:), allocatable          :: getFloats
 class(tPartitionedStringList), intent(in)           :: this
 character(len=*),              intent(in)           :: key
 real(pReal),   dimension(:),   intent(in), optional :: defaultVal
 type(tPartitionedStringList),  pointer              :: item
 integer(pInt)                                       :: i
 logical                                             :: found, &
                                                        cumulative

 cumulative = (key(1:1) == '(' .and. key(len_trim(key):len_trim(key)) == ')')
 found = .false.

 allocate(getFloats(0))

 item => this%next
 do while (associated(item))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (.not. cumulative) then
       deallocate(getFloats) ! use here rhs allocation with empty list
       allocate(getFloats(0))
     endif
     if (item%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     do i = 2_pInt, item%string%pos(1)
       getFloats = [getFloats,IO_FloatValue(item%string%val,item%string%pos,i)]
     enddo
   endif
   item => item%next
 end do

 if (.not. found) then
   if (present(defaultVal)) then; getFloats = defaultVal; else; call IO_error(140_pInt,ext_msg=key); endif
 endif

end function getFloats


!--------------------------------------------------------------------------------------------------
!> @brief gets array of integer values of for a given key from a linked list
!> @details for cumulative keys, "()", values from all occurrences are return. Otherwise only all
!! values from the last occurrence. If key is not found exits with error unless default is given.
!--------------------------------------------------------------------------------------------------
function getInts(this,key,defaultVal)
 use IO, only: &
   IO_error, &
   IO_stringValue, &
   IO_IntValue

 implicit none
 integer(pInt), dimension(:), allocatable            :: getInts
 class(tPartitionedStringList), intent(in)           :: this
 character(len=*),              intent(in)           :: key
 integer(pInt), dimension(:),   intent(in), optional :: defaultVal
 type(tPartitionedStringList),  pointer              :: item
 integer(pInt)                                       :: i
 logical                                             :: found, &
                                                        cumulative

 cumulative = (key(1:1) == '(' .and. key(len_trim(key):len_trim(key)) == ')')
 found = .false.

 allocate(getInts(0))

 item => this%next
 do while (associated(item))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (.not. cumulative) then
       deallocate(getInts) ! use here rhs allocation with empty list
       allocate(getInts(0))
     endif
     if (item%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     do i = 2_pInt, item%string%pos(1)
       getInts = [getInts,IO_IntValue(item%string%val,item%string%pos,i)]
     enddo
   endif
   item => item%next
 end do

 if (.not. found) then
   if (present(defaultVal)) then; getInts = defaultVal; else; call IO_error(140_pInt,ext_msg=key); endif
 endif

end function getInts


!--------------------------------------------------------------------------------------------------
!> @brief gets array of string values of for a given key from a linked list
!> @details for cumulative keys, "()", values from all occurrences are return. Otherwise only all
!! values from the last occurrence. If key is not found exits with error unless default is given.
!! If raw is true, the the complete string is returned, otherwise the individual chunks are returned
!--------------------------------------------------------------------------------------------------
function getStrings(this,key,defaultVal,raw)
 use IO, only: &
   IO_error, &
   IO_StringValue

 implicit none
 character(len=65536),dimension(:), allocatable           :: getStrings
 class(tPartitionedStringList),      intent(in)           :: this
 character(len=*),                   intent(in)           :: key
 character(len=65536),dimension(:),  intent(in), optional :: defaultVal
 logical,                            intent(in), optional :: raw
 type(tPartitionedStringList), pointer                    :: item
 character(len=65536)                                     :: str
 integer(pInt)                                            :: i
 logical                                                  :: found, &
                                                             whole, &
                                                             cumulative

 cumulative = (key(1:1) == '(' .and. key(len_trim(key):len_trim(key)) == ')')
 whole = merge(raw,.false.,present(raw))
 found = .false.

 item => this%next
 do while (associated(item))
   if (trim(IO_stringValue(item%string%val,item%string%pos,1)) == trim(key)) then
     found = .true.
     if (allocated(getStrings) .and. .not. cumulative) deallocate(getStrings)
     if (item%string%pos(1) < 2_pInt) call IO_error(143_pInt,ext_msg=key)
     
     notAllocated: if (.not. allocated(getStrings)) then
       if (whole) then
         str = item%string%val(item%string%pos(4):)
         getStrings = [str]
       else
         str = IO_StringValue(item%string%val,item%string%pos,2_pInt)
         allocate(getStrings(1),source=str)
         do i=3_pInt,item%string%pos(1)
           str = IO_StringValue(item%string%val,item%string%pos,i)
           getStrings = [getStrings,str]
         enddo
       endif
     else notAllocated
       if (whole) then
         getStrings = [getStrings,str]
       else
         do i=2_pInt,item%string%pos(1)
           str = IO_StringValue(item%string%val,item%string%pos,i)
           getStrings = [getStrings,str]
         enddo
       endif
     endif notAllocated
   endif
   item => item%next
 end do

 if (.not. found) then
   if (present(defaultVal)) then; getStrings = defaultVal; else; call IO_error(140_pInt,ext_msg=key); endif
 endif

end function getStrings


end module config
