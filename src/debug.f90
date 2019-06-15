!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Reading in and interpretating the debugging settings for the various modules
!--------------------------------------------------------------------------------------------------
module debug
  use prec
  use IO

  implicit none
  private

  integer, parameter, public :: &
    debug_LEVELSELECTIVE     = 2**0, &
    debug_LEVELBASIC         = 2**1, &
    debug_LEVELEXTENSIVE     = 2**2
  integer, parameter, private :: &
    debug_MAXGENERAL         = debug_LEVELEXTENSIVE                                                 ! must be set to the last bitcode used by (potentially) all debug types
  integer, parameter, public :: &
    debug_SPECTRALRESTART    = debug_MAXGENERAL*2**1, &
    debug_SPECTRALFFTW       = debug_MAXGENERAL*2**2, &
    debug_SPECTRALDIVERGENCE = debug_MAXGENERAL*2**3, &
    debug_SPECTRALROTATION   = debug_MAXGENERAL*2**4, &
    debug_SPECTRALPETSC      = debug_MAXGENERAL*2**5
    
  integer, parameter, public :: &
    debug_DEBUG                   =  1, &
    debug_MATH                    =  2, &
    debug_FESOLVING               =  3, &
    debug_MESH                    =  4, &                                                           !< stores debug level for mesh part of DAMASK bitwise coded
    debug_MATERIAL                =  5, &                                                           !< stores debug level for material part of DAMASK bitwise coded
    debug_LATTICE                 =  6, &                                                           !< stores debug level for lattice part of DAMASK bitwise coded
    debug_CONSTITUTIVE            =  7, &                                                           !< stores debug level for constitutive part of DAMASK bitwise coded
    debug_CRYSTALLITE             =  8, &
    debug_HOMOGENIZATION          =  9, &
    debug_CPFEM                   = 10, &
    debug_SPECTRAL                = 11, &
    debug_MARC                    = 12, &
    debug_ABAQUS                  = 13
  integer, parameter, private :: &
    debug_MAXNTYPE                = debug_ABAQUS                                                    !< must be set to the maximum defined debug type

  integer,protected, dimension(debug_maxNtype+2),  public :: &                                      ! specific ones, and 2 for "all" and "other"
    debug_level                    = 0

  integer, protected, public :: &
    debug_e                       = 1, &
    debug_i                       = 1, &
    debug_g                       = 1

  integer, dimension(2), public :: &
    debug_stressMaxLocation       = 0, &
    debug_stressMinLocation       = 0, &
    debug_jacobianMaxLocation     = 0, &
    debug_jacobianMinLocation     = 0


  real(pReal), public :: &
    debug_stressMax               = -huge(1.0_pReal), &
    debug_stressMin               =  huge(1.0_pReal), &
    debug_jacobianMax             = -huge(1.0_pReal), &
    debug_jacobianMin             =  huge(1.0_pReal)

#ifdef PETSc
  character(len=1024), parameter, public :: &
    PETSCDEBUG = ' -snes_view -snes_monitor '
#endif
  public :: debug_init, &
            debug_reset, &
            debug_info

contains


!--------------------------------------------------------------------------------------------------
!> @brief reads in parameters from debug.config and allocates arrays
!--------------------------------------------------------------------------------------------------
subroutine debug_init

  character(len=pStringLen), dimension(:), allocatable :: fileContent

  integer                            :: i, what, j
  integer, allocatable, dimension(:) :: chunkPos
  character(len=pStringLen)          :: tag, line
  logical :: fexist 

  write(6,'(/,a)')   ' <<<+-  debug init  -+>>>'
#ifdef DEBUG
  write(6,'(a)') achar(27)//'[31m <<<+-  DEBUG version  -+>>>'//achar(27)//'[0m'
#endif


  inquire(file='debug.config', exist=fexist)

  fileExists: if (fexist) then
    fileContent = IO_read_ASCII('debug.config') 
    do j=1, size(fileContent)
      line = fileContent(j)
      if (IO_isBlank(line)) cycle                                                                    ! skip empty lines
      chunkPos = IO_stringPos(line)
      tag = IO_lc(IO_stringValue(line,chunkPos,1))                                                   ! extract key
      select case(tag)
        case ('element','e','el')
          debug_e = IO_intValue(line,chunkPos,2)
        case ('integrationpoint','i','ip')
          debug_i = IO_intValue(line,chunkPos,2)
        case ('grain','g','gr')
          debug_g = IO_intValue(line,chunkPos,2)
      end select

      what = 0
      select case(tag)
        case ('debug')
          what = debug_DEBUG
        case ('math')
          what = debug_MATH
        case ('fesolving', 'fe')
          what = debug_FESOLVING
        case ('mesh')
          what = debug_MESH
        case ('material')
          what = debug_MATERIAL
        case ('lattice')
          what = debug_LATTICE
        case ('constitutive')
          what = debug_CONSTITUTIVE
        case ('crystallite')
          what = debug_CRYSTALLITE
        case ('homogenization')
          what = debug_HOMOGENIZATION
        case ('cpfem')
          what = debug_CPFEM
        case ('spectral')
          what = debug_SPECTRAL
        case ('marc')
          what = debug_MARC
        case ('abaqus')
          what = debug_ABAQUS
        case ('all')
          what = debug_MAXNTYPE + 1
        case ('other')
          what = debug_MAXNTYPE + 2
      end select
      if (what /= 0) then
        do i = 2, chunkPos(1)
          select case(IO_lc(IO_stringValue(line,chunkPos,i)))
            case('basic')
              debug_level(what) = ior(debug_level(what), debug_LEVELBASIC)
            case('extensive')
              debug_level(what) = ior(debug_level(what), debug_LEVELEXTENSIVE)
            case('selective')
              debug_level(what) = ior(debug_level(what), debug_LEVELSELECTIVE)
            case('restart')
              debug_level(what) = ior(debug_level(what), debug_SPECTRALRESTART)
            case('fft','fftw')
              debug_level(what) = ior(debug_level(what), debug_SPECTRALFFTW)
            case('divergence')
              debug_level(what) = ior(debug_level(what), debug_SPECTRALDIVERGENCE)
            case('rotation')
              debug_level(what) = ior(debug_level(what), debug_SPECTRALROTATION)
            case('petsc')
              debug_level(what) = ior(debug_level(what), debug_SPECTRALPETSC)
          end select
        enddo
       endif
    enddo

    do i = 1, debug_maxNtype
      if (debug_level(i) == 0) &
        debug_level(i) = ior(debug_level(i), debug_level(debug_MAXNTYPE + 2))                        ! fill undefined debug types with levels specified by "other"

      debug_level(i) = ior(debug_level(i), debug_level(debug_MAXNTYPE + 1))                          ! fill all debug types with levels specified by "all"
    enddo

    if (iand(debug_level(debug_debug),debug_LEVELBASIC) /= 0) &
      write(6,'(a,/)') ' using values from config file'
  else fileExists
    if (iand(debug_level(debug_debug),debug_LEVELBASIC) /= 0) &
      write(6,'(a,/)') ' using standard values'
  endif fileExists

!--------------------------------------------------------------------------------------------------
! output switched on (debug level for debug must be extensive)
   if (iand(debug_level(debug_debug),debug_LEVELEXTENSIVE) /= 0) then
     do i = 1, debug_MAXNTYPE
       select case(i)
         case (debug_DEBUG)
           tag = ' Debug'
         case (debug_MATH)
           tag = ' Math'
         case (debug_FESOLVING)
           tag = ' FEsolving'
         case (debug_MESH)
           tag = ' Mesh'
         case (debug_MATERIAL)
           tag = ' Material'
         case (debug_LATTICE)
           tag = ' Lattice'
         case (debug_CONSTITUTIVE)
           tag = ' Constitutive'
         case (debug_CRYSTALLITE)
           tag = ' Crystallite'
         case (debug_HOMOGENIZATION)
           tag = ' Homogenizaiton'
         case (debug_CPFEM)
           tag = ' CPFEM'
         case (debug_SPECTRAL)
           tag = ' Spectral solver'
         case (debug_MARC)
           tag = ' MSC.MARC FEM solver'
         case (debug_ABAQUS)
           tag = ' ABAQUS FEM solver'
       end select

       if(debug_level(i) /= 0) then
         write(6,'(3a)') ' debug level for ', trim(tag), ':'
         if(iand(debug_level(i),debug_LEVELBASIC)        /= 0) write(6,'(a)') '  basic'
         if(iand(debug_level(i),debug_LEVELEXTENSIVE)    /= 0) write(6,'(a)') '  extensive'
         if(iand(debug_level(i),debug_LEVELSELECTIVE)    /= 0) then
           write(6,'(a)') ' selective on:'
           write(6,'(a24,1x,i8)') '  element:              ',debug_e
           write(6,'(a24,1x,i8)') '  ip:                   ',debug_i
           write(6,'(a24,1x,i8)') '  grain:                ',debug_g
         endif
         if(iand(debug_level(i),debug_SPECTRALRESTART)   /= 0) write(6,'(a)') '  restart'
         if(iand(debug_level(i),debug_SPECTRALFFTW)      /= 0) write(6,'(a)') '  FFTW'
         if(iand(debug_level(i),debug_SPECTRALDIVERGENCE)/= 0) write(6,'(a)') '  divergence'
         if(iand(debug_level(i),debug_SPECTRALROTATION)  /= 0) write(6,'(a)') '  rotation'
         if(iand(debug_level(i),debug_SPECTRALPETSC)     /= 0) write(6,'(a)') '  PETSc'
       endif
     enddo
   endif

end subroutine debug_init


!--------------------------------------------------------------------------------------------------
!> @brief resets all debug values
!--------------------------------------------------------------------------------------------------
subroutine debug_reset

  debug_stressMaxLocation   = 0
  debug_stressMinLocation   = 0
  debug_jacobianMaxLocation = 0
  debug_jacobianMinLocation = 0
  debug_stressMax           = -huge(1.0_pReal)
  debug_stressMin           =  huge(1.0_pReal)
  debug_jacobianMax         = -huge(1.0_pReal)
  debug_jacobianMin         =  huge(1.0_pReal)

end subroutine debug_reset


!--------------------------------------------------------------------------------------------------
!> @brief writes debug statements to standard out
!--------------------------------------------------------------------------------------------------
subroutine debug_info

  !$OMP CRITICAL (write2out)
    debugOutputCPFEM: if (iand(debug_level(debug_CPFEM),debug_LEVELBASIC) /= 0 &
                       .and. any(debug_stressMinLocation /= 0) &
                       .and. any(debug_stressMaxLocation /= 0) ) then
      write(6,'(2/,a,/)') ' Extreme values of returned stress and Jacobian'
      write(6,'(a39)')                      '                      value     el   ip'
      write(6,'(a14,1x,e12.3,1x,i8,1x,i4)')   ' stress   min :', debug_stressMin, debug_stressMinLocation
      write(6,'(a14,1x,e12.3,1x,i8,1x,i4)')   '          max :', debug_stressMax, debug_stressMaxLocation
      write(6,'(a14,1x,e12.3,1x,i8,1x,i4)')   ' Jacobian min :', debug_jacobianMin, debug_jacobianMinLocation
      write(6,'(a14,1x,e12.3,1x,i8,1x,i4,/)') '          max :', debug_jacobianMax, debug_jacobianMaxLocation
    endif debugOutputCPFEM
  !$OMP END CRITICAL (write2out)

end subroutine debug_info

end module debug
