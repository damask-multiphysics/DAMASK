module config_material
 use chained_list
 use prec, only: &
   pReal, &
   pInt
 implicit none
 private 
 type(tPartitionedStringList), private,protected, allocatable, dimension(:) :: &
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
 character(len=*), parameter  :: &
   MATERIAL_partHomogenization = 'homogenization', &                                                !< keyword for homogenization part
   MATERIAL_partCrystallite    = 'crystallite', &                                                   !< keyword for crystallite part
   MATERIAL_partPhase          = 'phase',&                                                            !< keyword for phase part
   MATERIAL_partMicrostructure = 'microstructure', &                                                !< keyword for microstructure part
   MATERIAL_partTexture        = 'texture'                                                          !< keyword for texture part

 integer(pInt), public, protected :: &
   material_Ntexture, &                                                                             !< number of textures
   material_Nphase, &                                                                               !< number of phases
   material_Nhomogenization, &                                                                      !< number of homogenizations
   material_Nmicrostructure, &                                                                      !< number of microstructures
   material_Ncrystallite                                                                            !< number of crystallite settings


contains

subroutine config_material_init()
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
   debug_levelBasic, &
   debug_levelExtensive
 use mesh, only: &
   mesh_maxNips, &
   mesh_NcpElems, &
   mesh_element, &
   FE_Nips, &
   FE_geomtype

 implicit none
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt)            :: m,c,h, myDebug, myPhase, myHomog
 integer(pInt) :: &
  g, &                                                                                              !< grain number
  i, &                                                                                              !< integration point number
  e, &                                                                                              !< element number
  phase
 integer(pInt), dimension(:), allocatable :: ConstitutivePosition
 integer(pInt), dimension(:), allocatable :: CrystallitePosition
 integer(pInt), dimension(:), allocatable :: HomogenizationPosition

 character(len=65536) :: &                                                                          
  line,part

 character(len=*), parameter  :: &
   MATERIAL_configFile         = 'material.config', &                                               !< generic name for material configuration file
   MATERIAL_localFileExt       = 'materialConfig'                                                   !< extension of solver job name depending material configuration file


 myDebug = debug_level(debug_material)

 write(6,'(/,a)') ' <<<+-  material init  -+>>>'
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
       line = material_parsePhase(FILEUNIT)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Phase parsed'; flush(6)
    
     case (trim(material_partMicrostructure))
       line = material_parseMicrostructure(FILEUNIT)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Microstructure parsed'; flush(6)
    
     case (trim(material_partCrystallite))
       line = material_parseCrystallite(FILEUNIT)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Crystallite parsed'; flush(6)
    
     case (trim(material_partHomogenization))
       line = material_parseHomogenization(FILEUNIT)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Homogenization parsed'; flush(6)
    
     case (trim(material_partTexture))
       line = material_parseTexture(FILEUNIT)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Texture parsed'; flush(6)

     case default
       line = IO_read(fileUnit)

   end select

 enddo
end subroutine config_material_init

!--------------------------------------------------------------------------------------------------
!> @brief parses the homogenization part in the material configuration file
!--------------------------------------------------------------------------------------------------
character(len=65536) function material_parseHomogenization(fileUnit)
 use IO, only: &
   IO_read, &
   IO_globalTagInPart, &
   IO_countSections, &
   IO_error, &
   IO_countTagInPart, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringValue, &
   IO_intValue, &
   IO_floatValue, &
   IO_stringPos, &
   IO_EOF
 use mesh, only: &
   mesh_element

 implicit none
 integer(pInt),    intent(in) :: fileUnit


 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt)        :: Nsections,  h
 character(len=65536) :: line, tag,devNull
 character(len=64) ::  tag2
 logical              :: echo
 
 allocate(homogenizationConfig(0))

 h = 0_pInt
 do while (trim(line) /= IO_EOF)                                                                    ! read through sections of material part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   foundNextPart: if (IO_getTag(line,'<','>') /= '') then
     devNull = IO_read(fileUnit, .true.)                                                            ! reset IO_read
     exit
   endif foundNextPart
   nextSection: if (IO_getTag(line,'[',']') /= '') then
     h = h + 1_pInt
     homogenizationConfig = [homogenizationConfig, emptyList]
     tag2 = IO_getTag(line,'[',']')
     GfortranBug86033: if (.not. allocated(homogenization_name)) then
       allocate(homogenization_name(1),source=tag2)
     else GfortranBug86033
       homogenization_name  = [homogenization_name,tag2]
     endif GfortranBug86033
   endif nextSection
   chunkPos = IO_stringPos(line)
   tag = IO_lc(IO_stringValue(trim(line),chunkPos,1_pInt))                                          ! extract key
   inSection: if (h > 0_pInt) then
     chunkPos = IO_stringPos(line)
     call homogenizationConfig(h)%add(IO_lc(trim(line)),chunkPos)
   else inSection
     echo = (trim(tag) == '/echo/')
   endif inSection
 enddo

 material_Nhomogenization = size(homogenizationConfig)
 if (material_Nhomogenization < 1_pInt) call IO_error(160_pInt,ext_msg=material_partHomogenization)
 material_parseHomogenization=line

end function material_parseHomogenization


!--------------------------------------------------------------------------------------------------
!> @brief parses the microstructure part in the material configuration file
!--------------------------------------------------------------------------------------------------
character(len=65536) function material_parseMicrostructure(fileUnit)
 use prec, only: &
  dNeq
 use IO
 use mesh, only: &
   mesh_element, &
   mesh_NcpElems

 implicit none
 integer(pInt),    intent(in) :: fileUnit

 character(len=256), dimension(:), allocatable :: &
   str 
 character(len=64) :: tag2
 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt), allocatable, dimension(:,:) :: chunkPoss
 integer(pInt) :: e, m, constituent, i
 character(len=65536) :: &
   tag,line,devNull
 logical              :: echo

 allocate(MicrostructureConfig(0))
 line    = ''                                                                                       ! to have it initialized
 m       = 0_pInt
 echo    =.false.

 do while (trim(line) /= IO_EOF)                                                                    ! read through sections of material part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   foundNextPart: if (IO_getTag(line,'<','>') /= '') then
     devNull = IO_read(fileUnit, .true.)                                                            ! reset IO_read
     exit
   endif foundNextPart
   nextSection: if (IO_getTag(line,'[',']') /= '') then
     m = m + 1_pInt
     microstructureConfig = [microstructureConfig, emptyList]
     tag2 = IO_getTag(line,'[',']')
     GfortranBug86033: if (.not. allocated(microstructure_name)) then
       allocate(microstructure_name(1),source=tag2)
     else GfortranBug86033
       microstructure_name  = [microstructure_name,tag2]
     endif GfortranBug86033
   endif nextSection
   chunkPos = IO_stringPos(line)
   tag = IO_lc(IO_stringValue(trim(line),chunkPos,1_pInt))                                          ! extract key
   inSection: if (m > 0_pInt) then
     chunkPos = IO_stringPos(line)
     call microstructureConfig(m)%add(IO_lc(trim(line)),chunkPos)
   else inSection
     echo = (trim(tag) == '/echo/')
   endif inSection
 enddo

 material_Nmicrostructure = size(microstructureConfig)
 if (material_Nmicrostructure < 1_pInt) call IO_error(160_pInt,ext_msg=material_partMicrostructure)
 material_parseMicrostructure = line
end function material_parseMicrostructure


!--------------------------------------------------------------------------------------------------
!> @brief parses the crystallite part in the material configuration file
!--------------------------------------------------------------------------------------------------
character(len=65536) function material_parseCrystallite(fileUnit)
 use IO, only: &
   IO_read, &
   IO_error, &
   IO_getTag, &
   IO_lc, &
   IO_stringPos, &
   IO_stringValue, &
   IO_isBlank, &
   IO_EOF

 implicit none
 integer(pInt),    intent(in) :: fileUnit
 integer(pInt), allocatable, dimension(:) :: chunkPos

 character(len=64) ::  tag2
 integer(pInt)        :: c
 character(len=65536) :: line, tag,devNull
 logical              :: echo

 allocate(crystalliteConfig(0))
 c = 0_pInt
 do while (trim(line) /= IO_EOF)                                                                    ! read through sections of material part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   foundNextPart: if (IO_getTag(line,'<','>') /= '') then
     devNull = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit
   endif foundNextPart
   nextSection: if (IO_getTag(line,'[',']') /= '') then
     c = c + 1_pInt
     crystalliteConfig = [crystalliteConfig, emptyList]
     tag2 = IO_getTag(line,'[',']')
     GfortranBug86033: if (.not. allocated(crystallite_name)) then
       allocate(crystallite_name(1),source=tag2)
     else GfortranBug86033
       crystallite_name  = [crystallite_name,tag2]
     endif GfortranBug86033
   endif nextSection
   chunkPos = IO_stringPos(line)
   tag = IO_lc(IO_stringValue(trim(line),chunkPos,1_pInt))                                          ! extract key
   inSection: if (c > 0_pInt) then
     chunkPos = IO_stringPos(line)
     call crystalliteConfig(c)%add(IO_lc(trim(line)),chunkPos)
   else inSection
     echo = (trim(tag) == '/echo/')
   endif inSection
 enddo

 material_Ncrystallite = size(crystalliteConfig)
 if (material_Ncrystallite < 1_pInt) call IO_error(160_pInt,ext_msg=material_partCrystallite)
 material_parseCrystallite = line

end function material_parseCrystallite


!--------------------------------------------------------------------------------------------------
!> @brief parses the phase part in the material configuration file
!--------------------------------------------------------------------------------------------------
character(len=65536) function material_parsePhase(fileUnit)
 use chained_list, only: &
   emptyList 
 use IO, only: &
   IO_read, &
   IO_globalTagInPart, &
   IO_countSections, &
   IO_error, &
   IO_countTagInPart, &
   IO_getTag, &
   IO_spotTagInPart, &
   IO_lc, &
   IO_isBlank, &
   IO_stringValue, &
   IO_stringPos, &
   IO_EOF

 implicit none
 integer(pInt),    intent(in) :: fileUnit


 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: sourceCtr, kinematicsCtr, stiffDegradationCtr, p
 character(len=65536) :: &
  tag,line,devNull
 character(len=64) ::  tag2
 character(len=64), dimension(:), allocatable :: &
   str 
 logical              :: echo

 allocate(phaseConfig(0))
 line    = ''                                                                                       ! to have it initialized
 p = 0_pInt                                                                                         !  - " -
 echo    =.false.

 do while (trim(line) /= IO_EOF)                                                                    ! read through sections of material part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   foundNextPart: if (IO_getTag(line,'<','>') /= '') then
     devNull = IO_read(fileUnit, .true.)                                                            ! reset IO_read
     exit
   endif foundNextPart
   nextSection: if (IO_getTag(line,'[',']') /= '') then
     p = p + 1_pInt
     phaseConfig = [phaseConfig, emptyList]
     tag2 = IO_getTag(line,'[',']')
     GfortranBug86033: if (.not. allocated(phase_name)) then
       allocate(phase_name(1),source=tag2)
     else GfortranBug86033
       phase_name  = [phase_name,tag2]
     endif GfortranBug86033
   endif nextSection
   chunkPos = IO_stringPos(line)
   tag = IO_lc(IO_stringValue(trim(line),chunkPos,1_pInt))                                          ! extract key
   inSection: if (p > 0_pInt) then
     chunkPos = IO_stringPos(line)
     call phaseConfig(p)%add(IO_lc(trim(line)),chunkPos)
   else inSection
     echo = (trim(tag) == '/echo/')
   endif inSection
 enddo

 material_Nphase = size(phaseConfig)
 if (material_Nphase < 1_pInt) call IO_error(160_pInt,ext_msg=material_partPhase)
  material_parsePhase = line
end function material_parsePhase

!--------------------------------------------------------------------------------------------------
!> @brief parses the texture part in the material configuration file
!--------------------------------------------------------------------------------------------------
character(len=65536) function material_parseTexture(fileUnit)
 use prec, only: &
   dNeq
 use IO, only: &
   IO_read, &
   IO_globalTagInPart, &
   IO_countSections, &
   IO_error, &
   IO_countTagInPart, &
   IO_getTag, &
   IO_spotTagInPart, &
   IO_lc, &
   IO_isBlank, &
   IO_floatValue, &
   IO_stringValue, &
   IO_stringPos, &
   IO_EOF
 use math, only: &
   inRad, &
   math_sampleRandomOri, &
   math_I3, &
   math_det33, &
   math_inv33

 implicit none
 integer(pInt),    intent(in) :: fileUnit


 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: Nsections, section, gauss, fiber, j, t, i
 character(len=64) ::  tag2
 character(len=256), dimension(:), allocatable ::  bla
 logical              :: echo

 character(len=65536) :: line, tag,devNull, line2

 allocate(textureConfig(0))

 t = 0_pInt
 do while (trim(line2) /= IO_EOF)                                                                    ! read through sections of material part
   line2 = IO_read(fileUnit)
   if (IO_isBlank(line2)) cycle                                                                      ! skip empty lines
   foundNextPart: if (IO_getTag(line2,'<','>') /= '') then
     devNull = IO_read(fileUnit, .true.)                                                            ! reset IO_read
     exit
   endif foundNextPart
   nextSection: if (IO_getTag(line2,'[',']') /= '') then
     t = t + 1_pInt
     textureConfig = [textureConfig, emptyList]
     tag2 = IO_getTag(line2,'[',']')
     GfortranBug86033: if (.not. allocated(texture_name)) then
       allocate(texture_name(1),source=tag2)
     else GfortranBug86033
       texture_name  = [texture_name,tag2]
     endif GfortranBug86033
   endif nextSection
   chunkPos = IO_stringPos(line2)
   tag = IO_lc(IO_stringValue(trim(line2),chunkPos,1_pInt))                                          ! extract key
   inSection: if (t > 0_pInt) then
     chunkPos = IO_stringPos(line2)
     call textureConfig(t)%add(IO_lc(trim(line2)),chunkPos)
   else inSection
     echo = (trim(tag) == '/echo/')
   endif inSection
 enddo

 material_Ntexture = size(textureConfig)
 if (material_Ntexture < 1_pInt) call IO_error(160_pInt,ext_msg=material_partTexture)

 material_parseTexture = line2
end function material_parseTexture


end module config_material
