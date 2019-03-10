!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Reading in and interpretating the debugging settings for the various modules
!--------------------------------------------------------------------------------------------------
module debug
 use prec, only: &
   pInt, &
   pReal

 implicit none
 private
 integer(pInt), parameter, public :: &
   debug_LEVELSELECTIVE     = 2_pInt**0_pInt, &
   debug_LEVELBASIC         = 2_pInt**1_pInt, &
   debug_LEVELEXTENSIVE     = 2_pInt**2_pInt
 integer(pInt), parameter, private :: &
   debug_MAXGENERAL         = debug_LEVELEXTENSIVE                                                  ! must be set to the last bitcode used by (potentially) all debug types
 integer(pInt), parameter, public :: &
   debug_SPECTRALRESTART    = debug_MAXGENERAL*2_pInt**1_pInt, &
   debug_SPECTRALFFTW       = debug_MAXGENERAL*2_pInt**2_pInt, &
   debug_SPECTRALDIVERGENCE = debug_MAXGENERAL*2_pInt**3_pInt, &
   debug_SPECTRALROTATION   = debug_MAXGENERAL*2_pInt**4_pInt, &
   debug_SPECTRALPETSC      = debug_MAXGENERAL*2_pInt**5_pInt
   
 integer(pInt), parameter, public :: &
   debug_DEBUG                   =  1_pInt, &
   debug_MATH                    =  2_pInt, &
   debug_FESOLVING               =  3_pInt, &
   debug_MESH                    =  4_pInt, &                                                       !< stores debug level for mesh part of DAMASK bitwise coded
   debug_MATERIAL                =  5_pInt, &                                                       !< stores debug level for material part of DAMASK bitwise coded
   debug_LATTICE                 =  6_pInt, &                                                       !< stores debug level for lattice part of DAMASK bitwise coded
   debug_CONSTITUTIVE            =  7_pInt, &                                                       !< stores debug level for constitutive part of DAMASK bitwise coded
   debug_CRYSTALLITE             =  8_pInt, &
   debug_HOMOGENIZATION          =  9_pInt, &
   debug_CPFEM                   = 10_pInt, &
   debug_SPECTRAL                = 11_pInt, &
   debug_MARC                    = 12_pInt, &
   debug_ABAQUS                  = 13_pInt
 integer(pInt), parameter, private :: &
   debug_MAXNTYPE                = debug_ABAQUS                                                     !< must be set to the maximum defined debug type

 integer(pInt),protected, dimension(debug_maxNtype+2_pInt),  public :: &                            ! specific ones, and 2 for "all" and "other"
   debug_level                    = 0_pInt

 integer(pInt), protected, public :: &
   debug_e                       = 1_pInt, &
   debug_i                       = 1_pInt, &
   debug_g                       = 1_pInt

 integer(pInt), dimension(2), public :: &
   debug_stressMaxLocation       = 0_pInt, &
   debug_stressMinLocation       = 0_pInt, &
   debug_jacobianMaxLocation     = 0_pInt, &
   debug_jacobianMinLocation     = 0_pInt


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
 use prec, only: &
   pStringLen
 use IO, only: &
   IO_read_ASCII, &
   IO_error, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_lc, &
   IO_floatValue, &
   IO_intValue

 implicit none
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
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                              ! extract key
     select case(tag)
       case ('element','e','el')
         debug_e = IO_intValue(line,chunkPos,2_pInt)
       case ('integrationpoint','i','ip')
         debug_i = IO_intValue(line,chunkPos,2_pInt)
       case ('grain','g','gr')
         debug_g = IO_intValue(line,chunkPos,2_pInt)
     end select

     what = 0_pInt
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
         what = debug_MAXNTYPE + 1_pInt
       case ('other')
         what = debug_MAXNTYPE + 2_pInt
     end select
     if (what /= 0) then
       do i = 2_pInt, chunkPos(1)
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

   do i = 1_pInt, debug_maxNtype
     if (debug_level(i) == 0) &
       debug_level(i) = ior(debug_level(i), debug_level(debug_MAXNTYPE + 2_pInt))                   ! fill undefined debug types with levels specified by "other"

     debug_level(i) = ior(debug_level(i), debug_level(debug_MAXNTYPE + 1_pInt))                     ! fill all debug types with levels specified by "all"
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
     do i = 1_pInt, debug_MAXNTYPE
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

 implicit none

 debug_stressMaxLocation                   = 0_pInt
 debug_stressMinLocation                   = 0_pInt
 debug_jacobianMaxLocation                 = 0_pInt
 debug_jacobianMinLocation                 = 0_pInt
 debug_stressMax                           = -huge(1.0_pReal)
 debug_stressMin                           =  huge(1.0_pReal)
 debug_jacobianMax                         = -huge(1.0_pReal)
 debug_jacobianMin                         =  huge(1.0_pReal)

end subroutine debug_reset


!--------------------------------------------------------------------------------------------------
!> @brief writes debug statements to standard out
!--------------------------------------------------------------------------------------------------
subroutine debug_info

 implicit none
 
 !$OMP CRITICAL (write2out)
   debugOutputCPFEM: if (iand(debug_level(debug_CPFEM),debug_LEVELBASIC) /= 0 &
                      .and. any(debug_stressMinLocation /= 0_pInt) &
                      .and. any(debug_stressMaxLocation /= 0_pInt) ) then
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
