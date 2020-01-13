!-------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Reads in the material configuration from file
!> @details Reads the material configuration file, where solverJobName.materialConfig takes
!! precedence over material.config. Stores the raw strings and the positions of delimiters for the
!! parts 'homogenization', 'crystallite', 'phase', 'texture', and 'microstucture'
!! Reads numerics.config and debug.config
!--------------------------------------------------------------------------------------------------
module config
  use prec
  use DAMASK_interface
  use IO
  use debug
  use list

  implicit none
  private

  type(tPartitionedStringList), public, protected, allocatable, dimension(:) :: &
    config_phase, &
    config_microstructure, &
    config_homogenization, &
    config_texture, &
    config_crystallite
   
  type(tPartitionedStringList), public, protected :: &
    config_numerics, &
    config_debug
 
  character(len=pStringLen), dimension(:), allocatable, public, protected :: &
    config_name_phase, &                                                                            !< name of each phase
    config_name_homogenization, &                                                                   !< name of each homogenization
    config_name_crystallite, &                                                                      !< name of each crystallite setting
    config_name_microstructure, &                                                                   !< name of each microstructure
    config_name_texture                                                                             !< name of each texture

  public :: &
    config_init, &
    config_deallocate

contains

!--------------------------------------------------------------------------------------------------
!> @brief reads material.config and stores its content per part
!--------------------------------------------------------------------------------------------------
subroutine config_init

  integer :: i
  logical :: verbose

  character(len=pStringLen) :: &
    line, &
    part
  character(len=pStringLen), dimension(:), allocatable :: fileContent
  logical :: fileExists

  write(6,'(/,a)') ' <<<+-  config init  -+>>>'; flush(6)

  verbose = iand(debug_level(debug_material),debug_levelBasic) /= 0

  inquire(file=trim(getSolverJobName())//'.materialConfig',exist=fileExists)
  if(fileExists) then
    write(6,'(/,a)') ' reading '//trim(getSolverJobName())//'.materialConfig'; flush(6)
    fileContent = read_materialConfig(trim(getSolverJobName())//'.materialConfig')
  else
    inquire(file='material.config',exist=fileExists)
    if(.not. fileExists) call IO_error(100,ext_msg='material.config')
    write(6,'(/,a)') ' reading material.config'; flush(6)
    fileContent = read_materialConfig('material.config')
  endif

  do i = 1, size(fileContent)
    line = trim(fileContent(i))
    part = IO_lc(IO_getTag(line,'<','>'))
    select case (trim(part))
    
      case (trim('phase'))
        call parse_materialConfig(config_name_phase,config_phase,line,fileContent(i+1:))
        if (verbose) write(6,'(a)') ' Phase          parsed'; flush(6)
    
      case (trim('microstructure'))
        call parse_materialConfig(config_name_microstructure,config_microstructure,line,fileContent(i+1:))
        if (verbose) write(6,'(a)') ' Microstructure parsed'; flush(6)
    
      case (trim('crystallite'))
        call parse_materialConfig(config_name_crystallite,config_crystallite,line,fileContent(i+1:))
        if (verbose) write(6,'(a)') ' Crystallite    parsed'; flush(6)
        deallocate(config_crystallite)
    
      case (trim('homogenization'))
        call parse_materialConfig(config_name_homogenization,config_homogenization,line,fileContent(i+1:))
        if (verbose) write(6,'(a)') ' Homogenization parsed'; flush(6)
    
      case (trim('texture'))
        call parse_materialConfig(config_name_texture,config_texture,line,fileContent(i+1:))
        if (verbose) write(6,'(a)') ' Texture        parsed'; flush(6)

    end select

  enddo
 
  if (.not. allocated(config_homogenization) .or. size(config_homogenization) < 1) &
    call IO_error(160,ext_msg='<homogenization>')
  if (.not. allocated(config_microstructure) .or. size(config_microstructure) < 1) &
    call IO_error(160,ext_msg='<microstructure>')
  if (.not. allocated(config_phase)          .or. size(config_phase)          < 1) &
    call IO_error(160,ext_msg='<phase>')
  if (.not. allocated(config_texture)        .or. size(config_texture)        < 1) &
    call IO_error(160,ext_msg='<texture>')
 
 
  inquire(file='numerics.config', exist=fileExists)
  if (fileExists) then
    write(6,'(/,a)') ' reading numerics.config'; flush(6)
    fileContent = IO_read_ASCII('numerics.config')
    call parse_debugAndNumericsConfig(config_numerics,fileContent)
  endif
  
  inquire(file='debug.config', exist=fileExists)
  if (fileExists) then
    write(6,'(/,a)') ' reading debug.config'; flush(6)
    fileContent = IO_read_ASCII('debug.config')
    call parse_debugAndNumericsConfig(config_debug,fileContent)
  endif

contains


!--------------------------------------------------------------------------------------------------
!> @brief reads material.config
!!        Recursion is triggered by "{path/to/inputfile}" in a line
!--------------------------------------------------------------------------------------------------
recursive function read_materialConfig(fileName,cnt) result(fileContent)

  character(len=*),          intent(in)                :: fileName
  integer,                   intent(in), optional      :: cnt                                       !< recursion counter
  character(len=pStringLen), dimension(:), allocatable :: fileContent                               !< file content, separated per lines
  character(len=pStringLen), dimension(:), allocatable :: includedContent
  character(len=pStringLen)                            :: line
  character(len=pStringLen), parameter                 :: dummy = 'https://damask.mpie.de'          !< to fill up remaining array
  character(len=:),                        allocatable :: rawData
  integer ::  &
    fileLength, &
    fileUnit, &
    startPos, endPos, &
    myTotalLines, &                                                                                 !< # lines read from file without include statements
    l,i, &
    myStat
  logical :: warned
  
  if (present(cnt)) then
    if (cnt>10) call IO_error(106,ext_msg=trim(fileName))
  endif

!--------------------------------------------------------------------------------------------------
! read data as stream
  inquire(file = fileName, size=fileLength)
  if (fileLength == 0) then
    allocate(fileContent(0))
    return
  endif
  open(newunit=fileUnit, file=fileName, access='stream',&
       status='old', position='rewind', action='read',iostat=myStat)
  if(myStat /= 0) call IO_error(100,ext_msg=trim(fileName))
  allocate(character(len=fileLength)::rawData)
  read(fileUnit) rawData
  close(fileUnit)

!--------------------------------------------------------------------------------------------------
! count lines to allocate string array
  myTotalLines = 1
  do l=1, len(rawData)
    if (rawData(l:l) == new_line('')) myTotalLines = myTotalLines+1
  enddo
  allocate(fileContent(myTotalLines))

!--------------------------------------------------------------------------------------------------
! split raw data at end of line and handle includes
  warned = .false.
  startPos = 1
  l = 1
  do while (l <= myTotalLines)
    endPos = merge(startPos + scan(rawData(startPos:),new_line('')) - 2,len(rawData),l /= myTotalLines)
    if (endPos - startPos > pStringLen -1) then
      line = rawData(startPos:startPos+pStringLen-1)
      if (.not. warned) then
        call IO_warning(207,ext_msg=trim(fileName),el=l)
        warned = .true.
      endif
    else
      line = rawData(startPos:endpos)
    endif
    startPos = endPos + 2                                                                           ! jump to next line start

    recursion: if (scan(trim(adjustl(line)),'{') == 1 .and. scan(trim(line),'}') > 2) then
      includedContent = read_materialConfig(trim(line(scan(line,'{')+1:scan(line,'}')-1)), &
                        merge(cnt,1,present(cnt)))                                                  ! to track recursion depth
      fileContent     = [ fileContent(1:l-1), includedContent, [(dummy,i=1,myTotalLines-l)] ]       ! add content and grow array
      myTotalLines    = myTotalLines - 1 + size(includedContent)
      l               = l            - 1 + size(includedContent)
    else recursion
      fileContent(l) = line
      l = l + 1
    endif recursion

  enddo

end function read_materialConfig


!--------------------------------------------------------------------------------------------------
!> @brief parses the material.config file
!--------------------------------------------------------------------------------------------------
subroutine parse_materialConfig(sectionNames,part,line, &
                                fileContent)

  character(len=pStringLen),    allocatable, dimension(:), intent(out)   :: sectionNames
  type(tPartitionedStringList), allocatable, dimension(:), intent(inout) :: part
  character(len=pStringLen),                               intent(inout) :: line
  character(len=pStringLen),                 dimension(:), intent(in)    :: fileContent

  integer, allocatable, dimension(:) :: partPosition                                                !< position of [] tags + last line in section
  integer :: i, j
  logical :: echo
  character(len=pStringLen) :: sectionName

  echo = .false. 

  if (allocated(part)) call IO_error(161,ext_msg=trim(line))
  allocate(partPosition(0))
 
  do i = 1, size(fileContent)
    line = trim(fileContent(i))
    if (IO_getTag(line,'<','>') /= '') exit
    nextSection: if (IO_getTag(line,'[',']') /= '') then
      partPosition = [partPosition, i]
      cycle
    endif nextSection
    if (size(partPosition) < 1) &
      echo = (trim(IO_getTag(line,'/','/')) == 'echo') .or. echo
  enddo

  allocate(sectionNames(size(partPosition)))
  allocate(part(size(partPosition)))

  partPosition = [partPosition, i]                                                                  ! needed when actually storing content

  do i = 1, size(partPosition) -1
    write(sectionName,'(i0,a,a)') i,'_',trim(IO_getTag(fileContent(partPosition(i)),'[',']'))
    sectionNames(i) = sectionName
    do j = partPosition(i) + 1,  partPosition(i+1) -1
      call part(i)%add(trim(adjustl(fileContent(j))))
    enddo
    if (echo) then
      write(6,*) 'section',i, '"'//trim(sectionNames(i))//'"'
      call part(i)%show()
    endif
  enddo

end subroutine parse_materialConfig


!--------------------------------------------------------------------------------------------------
!> @brief parses the material.config file
!--------------------------------------------------------------------------------------------------
subroutine parse_debugAndNumericsConfig(config_list, &
                                        fileContent)

  type(tPartitionedStringList),              intent(out) :: config_list
  character(len=pStringLen),   dimension(:), intent(in)  :: fileContent
  integer :: i

  do i = 1, size(fileContent)
    call config_list%add(trim(adjustl(fileContent(i))))
  enddo

end subroutine parse_debugAndNumericsConfig

end subroutine config_init


!--------------------------------------------------------------------------------------------------
!> @brief deallocates the linked lists that store the content of the configuration files
!--------------------------------------------------------------------------------------------------
subroutine config_deallocate(what)

  character(len=*), intent(in) :: what

  select case(trim(what))

    case('material.config/phase')
      deallocate(config_phase)

    case('material.config/microstructure')
      deallocate(config_microstructure)

    case('material.config/homogenization')
      deallocate(config_homogenization)

    case('material.config/texture')
      deallocate(config_texture)
     
    case('debug.config')
      call config_debug%free
     
    case('numerics.config')
      call config_numerics%free
     
    case default
      call IO_error(0,ext_msg='config_deallocate')

  end select

end subroutine config_deallocate

end module config
