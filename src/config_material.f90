!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Reads in the material configuration from file
!> @details Reads the material configuration file, where solverJobName.materialConfig takes
!! precedence over material.config. Stores the raw strings and the positions of delimiters for the
!! parts 'homogenization', 'crystallite', 'phase', 'texture', and 'microstucture'
!--------------------------------------------------------------------------------------------------
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
       call parseFile(phase_name,phaseConfig,FILEUNIT,line)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Phase parsed'; flush(6)
    
     case (trim(material_partMicrostructure))
       call parseFile(microstructure_name,microstructureConfig,FILEUNIT,line)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Microstructure parsed'; flush(6)
    
     case (trim(material_partCrystallite))
       call parseFile(crystallite_name,crystalliteConfig,FILEUNIT,line)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Crystallite parsed'; flush(6)
    
     case (trim(material_partHomogenization))
       call parseFile(homogenization_name,homogenizationConfig,FILEUNIT,line)
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Homogenization parsed'; flush(6)
    
     case (trim(material_partTexture))
       call parseFile(texture_name,textureConfig,FILEUNIT,line)
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


end subroutine config_material_init

!--------------------------------------------------------------------------------------------------
!> @brief parses the homogenization part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine parseFile(partLabel,part,fileUnit,nextLine)
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


 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt)        :: Nsections,  h
 character(len=65536) :: line, tag,devNull
 character(len=65536) :: nextLine
 character(len=64) ::  tag2
 logical              :: echo
 type(tPartitionedStringList), allocatable, dimension(:), intent(inout) :: &
   part
 character(len=*),  dimension(:), allocatable, intent(inout)  :: partLabel
 
 allocate(part(0))

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
     part = [part, emptyList]
     tag2 = IO_getTag(line,'[',']')
     GfortranBug86033: if (.not. allocated(partLabel)) then
       allocate(partLabel(1),source=tag2)
     else GfortranBug86033
       partLabel  = [partLabel,tag2]
     endif GfortranBug86033
   endif nextSection
   chunkPos = IO_stringPos(line)
   tag = IO_lc(IO_stringValue(trim(line),chunkPos,1_pInt))                                          ! extract key
   inSection: if (h > 0_pInt) then
     chunkPos = IO_stringPos(line)
     call part(h)%add(IO_lc(trim(line)))
   else inSection
     echo = (trim(tag) == '/echo/')
   endif inSection
 enddo

end subroutine parseFile

end module config_material
