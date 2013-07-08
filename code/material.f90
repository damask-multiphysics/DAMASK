! Copyright 2011-13 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!--------------------------------------------------------------------------------------------------
!* $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Parses material config file, either solverJobName.materialConfig or material.config
!> @details reads the material configuration file, where solverJobName.materialConfig takes
!! precedence over material.config and parses the sections 'homogenization', 'crystallite',
!! 'phase', 'texture', and 'microstucture'
!--------------------------------------------------------------------------------------------------
module material
 use prec, only: &
   pReal, &
   pInt, &
   p_intvec

 implicit none
 private
 character(len=64), parameter, public  :: &
   MATERIAL_configFile         = 'material.config', &                                               !< generic name for material configuration file 
   MATERIAL_localFileExt       = 'materialConfig'                                                   !< extension of solver job name depending material configuration file  
   
 character(len=32), parameter, public  :: &
   MATERIAL_partHomogenization = 'homogenization', &                                                !< keyword for homogenization part
   MATERIAL_partCrystallite    = 'crystallite', &                                                   !< keyword for crystallite part
   MATERIAL_partPhase          = 'phase'                                                            !< keyword for phase part
 
 character(len=64), dimension(:), allocatable, public, protected :: &
   phase_elasticity, &                                                                              !< elasticity of each phase  
   phase_plasticity, &                                                                              !< plasticity of each phase  
   phase_name, &                                                                                    !< name of each phase
   homogenization_name, &                                                                           !< name of each homogenization
   homogenization_type, &                                                                           !< type of each homogenization
   crystallite_name                                                                                 !< name of each crystallite setting

 integer(pInt), public, protected :: &
   homogenization_maxNgrains, &                                                                     !< max number of grains in any USED homogenization
   material_Nphase, &                                                                               !< number of phases
   material_Nhomogenization, &                                                                      !< number of homogenizations
   material_Nmicrostructure, &                                                                      !< number of microstructures
   material_Ncrystallite                                                                            !< number of crystallite settings
 
 integer(pInt), dimension(:), allocatable, public, protected :: &
   homogenization_Ngrains, &                                                                        !< number of grains in each homogenization
   homogenization_Noutput, &                                                                        !< number of '(output)' items per homogenization
   phase_Noutput, &                                                                                 !< number of '(output)' items per phase
   phase_elasticityInstance, &                                                                      !< instance of particular elasticity of each phase
   phase_plasticityInstance, &                                                                      !< instance of particular plasticity of each phase
   crystallite_Noutput, &                                                                           !< number of '(output)' items per crystallite setting
   homogenization_typeInstance, &                                                                   !< instance of particular type of each homogenization
   microstructure_crystallite                                                                       !< crystallite setting ID of each microstructure

 integer(pInt), dimension(:,:,:), allocatable, public:: &
   material_phase                                                                                   !< phase   (index) of each grain,IP,element
 integer(pInt), dimension(:,:,:), allocatable, public, protected :: &
   material_texture                                                                                 !< texture (index) of each grain,IP,element
 
 real(pReal), dimension(:,:,:,:), allocatable, public, protected :: &
   material_EulerAngles                                                                             !< initial orientation of each grain,IP,element
 
 logical, dimension(:), allocatable, public, protected :: &
   microstructure_active, & 
   microstructure_elemhomo, &                                                                       !< flag to indicate homogeneous microstructure distribution over element's IPs
   phase_localPlasticity                                                                            !< flags phases with local constitutive law


 character(len=32), parameter, private :: &
   MATERIAL_partMicrostructure = 'microstructure', &                                                !< keyword for microstructure part
   MATERIAL_partTexture        = 'texture'                                                          !< keyword for texture part
   
 character(len=64), dimension(:), allocatable, private :: &
   microstructure_name, &                                                                           !< name of each microstructure
   texture_name                                                                                     !< name of each texture
     
 character(len=256), dimension(:), allocatable, private :: &
   texture_ODFfile                                                                                  !< name of each ODF file         

 integer(pInt), private :: &
   material_Ntexture, &                                                                             !< number of textures
   microstructure_maxNconstituents, &                                                               !< max number of constituents in any phase
   texture_maxNgauss, &                                                                             !< max number of Gauss components in any texture
   texture_maxNfiber                                                                                !< max number of Fiber components in any texture

 integer(pInt), dimension(:), allocatable, private :: &
   microstructure_Nconstituents, &                                                                  !< number of constituents in each microstructure
   texture_symmetry, &                                                                              !< number of symmetric orientations per texture
   texture_Ngauss, &                                                                                !< number of Gauss components per texture
   texture_Nfiber                                                                                   !< number of Fiber components per texture
 
 integer(pInt), dimension(:,:), allocatable, private :: &
   microstructure_phase, &                                                                          !< phase IDs of each microstructure
   microstructure_texture                                                                           !< texture IDs of each microstructure
 
 real(pReal), dimension(:,:), allocatable, private :: &
   microstructure_fraction                                                                          !< vol fraction of each constituent in microstructure
 
 real(pReal), dimension(:,:,:), allocatable, private :: &
   material_volume, &                                                                               !< volume of each grain,IP,element
   texture_Gauss, &                                                                                 !< data of each Gauss component
   texture_Fiber, &                                                                                 !< data of each Fiber component
   texture_rotation                                                                                 !< rotation of each texture
 
 logical, dimension(:), allocatable, private :: &
   homogenization_active


 public  :: material_init
 
 private :: material_parseHomogenization, &
            material_parseMicrostructure, &
            material_parseCrystallite, &
            material_parsePhase, &
            material_parseTexture, &
            material_populateGrains
contains


!--------------------------------------------------------------------------------------------------
!> @brief parses material configuration file
!> @details figures out if solverJobName.materialConfig is present, if not looks for 
!> material.config
!--------------------------------------------------------------------------------------------------
subroutine material_init
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use IO, only: &
   IO_error, &
   IO_open_file, &
   IO_open_jobFile_stat, &
   IO_timeStamp
 use debug, only: &
   debug_level, &
   debug_material, &
   debug_levelBasic, &
   debug_levelExtensive
 
 implicit none
 integer(pInt), parameter :: fileunit = 200_pInt
 integer(pInt)            :: m,c,h, myDebug
 myDebug = debug_level(debug_material)
 
 write(6,'(/,a)') ' <<<+-  material init  -+>>>'
 write(6,'(a)') ' $Id$'
 write(6,'(a16,a)')   ' Current time : ',IO_timeStamp()
#include "compilation_info.f90"

 if (.not. IO_open_jobFile_stat(fileunit,material_localFileExt)) then                               ! no local material configuration present...
   call IO_open_file(fileunit,material_configFile)                                                  ! ...open material.config file
 endif
 call material_parseHomogenization(fileunit,material_partHomogenization)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
   write(6,'(a)') ' Homogenization parsed'
 endif
 call material_parseMicrostructure(fileunit,material_partMicrostructure)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
   write(6,'(a)') ' Microstructure parsed'
 endif
 call material_parseCrystallite(fileunit,material_partCrystallite)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
   write(6,'(a)') ' Crystallite parsed'
 endif
 call material_parseTexture(fileunit,material_partTexture)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
   write(6,'(a)') ' Texture parsed'
 endif
 call material_parsePhase(fileunit,material_partPhase)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
   write(6,'(a)') ' Phase parsed'
 endif
 close(fileunit)

 do m = 1_pInt,material_Nmicrostructure
   if (microstructure_crystallite(m) < 1_pInt .or. &
       microstructure_crystallite(m) > material_Ncrystallite) & 
         call IO_error(150_pInt,m)
   if (minval(microstructure_phase(1:microstructure_Nconstituents(m),m)) < 1_pInt .or. &
       maxval(microstructure_phase(1:microstructure_Nconstituents(m),m)) > material_Nphase) &
         call IO_error(151_pInt,m)
   if (minval(microstructure_texture(1:microstructure_Nconstituents(m),m)) < 1_pInt .or. &
       maxval(microstructure_texture(1:microstructure_Nconstituents(m),m)) > material_Ntexture) &
         call IO_error(152_pInt,m)
!   if (abs(sum(microstructure_fraction(:,m)) - 1.0_pReal) >= 1.0e-6_pReal) then                     ! have ppm precision in fractions
!     if (iand(myDebug,debug_levelExtensive) /= 0_pInt) then
!         write(6,'(a,1x,f12.9)') ' sum of microstructure fraction = ',sum(microstructure_fraction(:,m))
!     endif
!     call IO_error(153_pInt,m)
!   endif
 enddo
 debugOut: if (iand(myDebug,debug_levelExtensive) /= 0_pInt) then
   write(6,'(/,a,/)') ' MATERIAL configuration'
   write(6,'(a32,1x,a16,1x,a6)') 'homogenization                  ','type            ','grains'
   do h = 1_pInt,material_Nhomogenization
     write(6,'(1x,a32,1x,a16,1x,i6)') homogenization_name(h),homogenization_type(h),homogenization_Ngrains(h)
   enddo
   write(6,'(/,a14,18x,1x,a11,1x,a12,1x,a13)') 'microstructure','crystallite','constituents','homogeneous'
   do m = 1_pInt,material_Nmicrostructure
     write(6,'(1x,a32,1x,i11,1x,i12,1x,l13)') microstructure_name(m), &
                                        microstructure_crystallite(m), &
                                        microstructure_Nconstituents(m), &
                                        microstructure_elemhomo(m)
     if (microstructure_Nconstituents(m) > 0_pInt) then
       do c = 1_pInt,microstructure_Nconstituents(m)
         write(6,'(a1,1x,a32,1x,a32,1x,f7.4)') '>',phase_name(microstructure_phase(c,m)),&
                                                   texture_name(microstructure_texture(c,m)),&
                                                   microstructure_fraction(c,m)
       enddo
       write(6,*)
     endif
   enddo
 endif debugOut
 
 call material_populateGrains

end subroutine material_init


!--------------------------------------------------------------------------------------------------
!> @brief parses the homogenization part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseHomogenization(myFile,myPart)
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
   IO_stringPos
 use mesh, only: &
   mesh_element
 
 implicit none
 character(len=*), intent(in) :: myPart
 integer(pInt),    intent(in) :: myFile
 
 integer(pInt),     parameter :: maxNchunks = 2_pInt
 
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) Nsections, section, s
 character(len=65536) :: tag
 character(len=65536) :: line
 logical              :: echo
 
 echo = IO_globalTagInPart(myFile,myPart,'/echo/')
 
 Nsections = IO_countSections(myFile,myPart)
 material_Nhomogenization = Nsections
 if (Nsections < 1_pInt) call IO_error(160_pInt,ext_msg=myPart)
 
 allocate(homogenization_name(Nsections));          homogenization_name = ''
 allocate(homogenization_type(Nsections));          homogenization_type = ''
 allocate(homogenization_typeInstance(Nsections));  homogenization_typeInstance = 0_pInt
 allocate(homogenization_Ngrains(Nsections));       homogenization_Ngrains = 0_pInt
 allocate(homogenization_Noutput(Nsections));       homogenization_Noutput = 0_pInt
 allocate(homogenization_active(Nsections));        homogenization_active = .false.

 forall (s = 1_pInt:Nsections) homogenization_active(s) = any(mesh_element(3,:) == s)               ! current homogenization used in model? Homogenization view, maximum operations depend on maximum number of homog schemes
   homogenization_Noutput = IO_countTagInPart(myFile,myPart,'(output)',Nsections)
 
 rewind(myFile)
 line = ''
 section = 0_pInt
 
 do while (trim(line) /= '#EOF#' .and. IO_lc(IO_getTag(line,'<','>')) /= myPart)                    ! wind forward to myPart
   line = IO_read(myFile)
 enddo
 if (echo) write(6,'(/,a)') trim(line)                                                              ! echo part header

 do while (trim(line) /= '#EOF#')
   line = IO_read(myFile)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (echo) write(6,'(a)') trim(line)                                                              ! echo back read lines
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     homogenization_name(section) = IO_getTag(line,'[',']')
   endif
   if (section > 0_pInt) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('type')
         homogenization_type(section) = IO_lc(IO_stringValue(line,positions,2_pInt))                ! adding: IO_lc function
         do s = 1_pInt,section
           if (homogenization_type(s) == homogenization_type(section)) &
             homogenization_typeInstance(section) = homogenization_typeInstance(section) + 1_pInt   ! count instances
         enddo
       case ('ngrains')
         homogenization_Ngrains(section) = IO_intValue(line,positions,2_pInt)
     end select
   endif
 enddo

 homogenization_maxNgrains = maxval(homogenization_Ngrains,homogenization_active)

end subroutine material_parseHomogenization


!--------------------------------------------------------------------------------------------------
!> @brief parses the microstructure part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseMicrostructure(myFile,myPart)
 use IO
 use mesh, only: &
   mesh_element, &
   mesh_NcpElems
 
 implicit none
 character(len=*), intent(in) :: myPart
 integer(pInt),    intent(in) :: myFile
 
 integer(pInt), parameter :: maxNchunks = 7_pInt
 
 integer(pInt), dimension(1_pInt+2_pInt*maxNchunks) :: positions
 integer(pInt) :: Nsections, section, constituent, e, i
 character(len=65536) :: tag
 character(len=65536) :: line
 logical              :: echo

 echo = IO_globalTagInPart(myFile,myPart,'/echo/')

 Nsections = IO_countSections(myFile,myPart)
 material_Nmicrostructure = Nsections
 if (Nsections < 1_pInt) call IO_error(160_pInt,ext_msg=myPart)

 allocate(microstructure_name(Nsections));            microstructure_name = ''
 allocate(microstructure_crystallite(Nsections));     microstructure_crystallite = 0_pInt
 allocate(microstructure_Nconstituents(Nsections))
 allocate(microstructure_active(Nsections))
 allocate(microstructure_elemhomo(Nsections))

 forall (e = 1_pInt:mesh_NcpElems) microstructure_active(mesh_element(4,e)) = .true.                ! current microstructure used in model? Elementwise view, maximum N operations for N elements
  
 microstructure_Nconstituents = IO_countTagInPart(myFile,myPart,'(constituent)',Nsections)
 microstructure_maxNconstituents = maxval(microstructure_Nconstituents)
 microstructure_elemhomo = IO_spotTagInPart(myFile,myPart,'/elementhomogeneous/',Nsections)

 allocate(microstructure_phase   (microstructure_maxNconstituents,Nsections))
   microstructure_phase    = 0_pInt
 allocate(microstructure_texture (microstructure_maxNconstituents,Nsections))
   microstructure_texture  = 0_pInt
 allocate(microstructure_fraction(microstructure_maxNconstituents,Nsections))
   microstructure_fraction = 0.0_pReal
 
 rewind(myFile)
 line        = ''                                                                                   ! to have it initialized
 section     = 0_pInt                                                                               !  - " -
 constituent = 0_pInt                                                                               !  - " -
 
 do while (trim(line) /= '#EOF#' .and. IO_lc(IO_getTag(line,'<','>')) /= myPart)                    ! wind forward to myPart
   line = IO_read(myFile)
 enddo
 if (echo) write(6,'(/,a)') trim(line)                                                              ! echo part header

 do while (trim(line) /= '#EOF#')
   line = IO_read(myFile)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (echo) write(6,'(a)') trim(line)                                                              ! echo back read lines
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     constituent = 0_pInt
     microstructure_name(section) = IO_getTag(line,'[',']')
   endif
   if (section > 0_pInt) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('crystallite')
         microstructure_crystallite(section) = IO_intValue(line,positions,2_pInt)
       case ('(constituent)')
         constituent = constituent + 1_pInt
         do i=2_pInt,6_pInt,2_pInt
           tag = IO_lc(IO_stringValue(line,positions,i))
           select case (tag)
             case('phase')
               microstructure_phase(constituent,section) =    IO_intValue(line,positions,i+1_pInt)
             case('texture')
               microstructure_texture(constituent,section) =  IO_intValue(line,positions,i+1_pInt)
             case('fraction')
               microstructure_fraction(constituent,section) = IO_floatValue(line,positions,i+1_pInt)
           end select
         enddo
     end select
   endif
 enddo

end subroutine material_parseMicrostructure


!--------------------------------------------------------------------------------------------------
!> @brief parses the crystallite part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseCrystallite(myFile,myPart)
 use IO, only: &
   IO_read, &
   IO_countSections, &
   IO_error, &
   IO_countTagInPart, &
   IO_globalTagInPart, &
   IO_getTag, &
   IO_lc, &
   IO_isBlank

 implicit none
 character(len=*), intent(in) :: myPart
 integer(pInt),    intent(in) :: myFile
 
 integer(pInt)        :: Nsections, &
                         section
 character(len=65536) :: line
 logical              :: echo

 echo = IO_globalTagInPart(myFile,myPart,'/echo/')
 
 Nsections = IO_countSections(myFile,myPart)
 material_Ncrystallite = Nsections
 if (Nsections < 1_pInt) call IO_error(160_pInt,ext_msg=myPart)

 allocate(crystallite_name(Nsections));       crystallite_name = ''
 allocate(crystallite_Noutput(Nsections));    crystallite_Noutput = 0_pInt

 crystallite_Noutput = IO_countTagInPart(myFile,myPart,'(output)',Nsections)
 
 rewind(myFile)
 line = ''
 section = 0_pInt
 
 do while (trim(line) /= '#EOF#' .and. IO_lc(IO_getTag(line,'<','>')) /= myPart)                    ! wind forward to myPart
   line = IO_read(myFile)
 enddo
 if (echo) write(6,'(/,a)') trim(line)                                                              ! echo part header

 do while (trim(line) /= '#EOF#')
   line = IO_read(myFile)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (echo) write(6,'(a)') trim(line)                                                              ! echo back read lines
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     crystallite_name(section) = IO_getTag(line,'[',']')
   endif
 enddo

end subroutine material_parseCrystallite


!--------------------------------------------------------------------------------------------------
!> @brief parses the phase part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parsePhase(myFile,myPart)
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
   IO_stringPos

 implicit none
 character(len=*), intent(in) :: myPart
 integer(pInt),    intent(in) :: myFile
 
 integer(pInt), parameter :: maxNchunks = 2_pInt
 
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) Nsections, section, s
 character(len=65536) :: tag
 character(len=65536) :: line
 logical              :: echo

 echo = IO_globalTagInPart(myFile,myPart,'/echo/')
 
 Nsections = IO_countSections(myFile,myPart)
 material_Nphase = Nsections
 if (Nsections < 1_pInt) call IO_error(160_pInt,ext_msg=myPart)

 allocate(phase_name(Nsections));                phase_name = ''
 allocate(phase_elasticity(Nsections));          phase_elasticity = ''
 allocate(phase_elasticityInstance(Nsections));  phase_elasticityInstance = 0_pInt
 allocate(phase_plasticity(Nsections));          phase_plasticity = ''
 allocate(phase_plasticityInstance(Nsections));  phase_plasticityInstance = 0_pInt
 allocate(phase_Noutput(Nsections));             phase_Noutput = 0_pInt
 allocate(phase_localPlasticity(Nsections));     phase_localPlasticity = .false.

 phase_Noutput = IO_countTagInPart(myFile,myPart,'(output)',Nsections)
 phase_localPlasticity = .not. IO_spotTagInPart(myFile,myPart,'/nonlocal/',Nsections)
 
 rewind(myFile)
 line = ''
 section = 0_pInt
 
 do while (trim(line) /= '#EOF#' .and. IO_lc(IO_getTag(line,'<','>')) /= myPart)                    ! wind forward to myPart
   line = IO_read(myFile)
 enddo
 if (echo) write(6,'(/,a)') trim(line)                                                              ! echo part header

 do while (trim(line) /= '#EOF#')
   line = IO_read(myFile)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (echo) write(6,'(a)') trim(line)                                                              ! echo back read lines
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     phase_name(section) = IO_getTag(line,'[',']')
   endif
   if (section > 0_pInt) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('elasticity')
         phase_elasticity(section) = IO_lc(IO_stringValue(line,positions,2_pInt))
         do s = 1_pInt,section
           if (phase_elasticity(s) == phase_elasticity(section)) &
             phase_elasticityInstance(section) = phase_elasticityInstance(section) + 1_pInt         ! count instances
         enddo
       case ('plasticity')
         phase_plasticity(section) = IO_lc(IO_stringValue(line,positions,2_pInt))
         do s = 1_pInt,section
           if (phase_plasticity(s) == phase_plasticity(section)) &
             phase_plasticityInstance(section) = phase_plasticityInstance(section) + 1_pInt         ! count instances
         enddo
     end select
   endif
 enddo

end subroutine material_parsePhase


!--------------------------------------------------------------------------------------------------
!> @brief parses the texture part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseTexture(myFile,myPart)
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
   IO_stringPos
 use math, only: &
   inRad, &
   math_sampleRandomOri, &
   math_I3, &
   math_inv33
 
 implicit none
 character(len=*), intent(in) :: myPart
 integer(pInt),    intent(in) :: myFile
 
 integer(pInt), parameter     :: maxNchunks = 13_pInt
 
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) :: Nsections, section, gauss, fiber, j
 character(len=65536) :: tag
 character(len=65536) :: line
 logical              :: echo

 echo = IO_globalTagInPart(myFile,myPart,'/echo/')
 
 Nsections = IO_countSections(myFile,myPart)
 material_Ntexture = Nsections
 if (Nsections < 1_pInt) call IO_error(160_pInt,ext_msg=myPart)

 allocate(texture_name(Nsections));     texture_name = ''
 allocate(texture_ODFfile(Nsections));  texture_ODFfile = ''
 allocate(texture_symmetry(Nsections)); texture_symmetry = 1_pInt
 allocate(texture_Ngauss(Nsections));   texture_Ngauss = 0_pInt
 allocate(texture_Nfiber(Nsections));   texture_Nfiber = 0_pInt

 texture_Ngauss = IO_countTagInPart(myFile,myPart,'(gauss)', Nsections) + &
                  IO_countTagInPart(myFile,myPart,'(random)',Nsections)
 texture_Nfiber = IO_countTagInPart(myFile,myPart,'(fiber)', Nsections)
 texture_maxNgauss = maxval(texture_Ngauss)
 texture_maxNfiber = maxval(texture_Nfiber)
 allocate(texture_Gauss   (5,texture_maxNgauss,Nsections)); texture_Gauss    = 0.0_pReal
 allocate(texture_Fiber   (6,texture_maxNfiber,Nsections)); texture_Fiber    = 0.0_pReal
 allocate(texture_rotation(3,3,Nsections));                 
 do j = 1_pInt, Nsections
   texture_rotation(1:3,1:3,j) = math_I3
 enddo
 
 rewind(myFile)
 line    = ''                                                                                       ! to have in initialized
 section = 0_pInt                                                                                   ! - " -
 gauss   = 0_pInt                                                                                   ! - " - 
 fiber   = 0_pInt                                                                                   ! - " - 
 
 do while (trim(line) /= '#EOF#' .and. IO_lc(IO_getTag(line,'<','>')) /= myPart)                    ! wind forward to myPart
   line = IO_read(myFile)
 enddo
 if (echo) write(6,'(/,a)') trim(line)                                                              ! echo part header

 do while (trim(line) /= '#EOF#')
   line = IO_read(myFile)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (echo) write(6,'(a)') trim(line)                                                              ! echo back read lines
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     gauss = 0_pInt
     fiber = 0_pInt
     texture_name(section) = IO_getTag(line,'[',']')
   endif
   if (section > 0_pInt) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     textureType: select case(tag)

       case ('rotation') textureType
         do j = 1_pInt, 3_pInt                                                                      ! look for "x", "y", and "z" entries
           tag = IO_lc(IO_stringValue(line,positions,j+1_pInt))
           select case (tag)
             case('x', '+x')
               texture_rotation(j,1:3,section) = (/ 1.0_pReal, 0.0_pReal, 0.0_pReal/)               ! original axis is now +x-axis
             case('-x')
               texture_rotation(j,1:3,section) = (/-1.0_pReal, 0.0_pReal, 0.0_pReal/)               ! original axis is now -x-axis
             case('y', '+y')
               texture_rotation(j,1:3,section) = (/ 0.0_pReal, 1.0_pReal, 0.0_pReal/)               ! original axis is now +y-axis
             case('-y')
               texture_rotation(j,1:3,section) = (/ 0.0_pReal,-1.0_pReal, 0.0_pReal/)               ! original axis is now -y-axis
             case('z', '+z')
               texture_rotation(j,1:3,section) = (/ 0.0_pReal, 0.0_pReal, 1.0_pReal/)               ! original axis is now +z-axis
             case('-z')
               texture_rotation(j,1:3,section) = (/ 0.0_pReal, 0.0_pReal,-1.0_pReal/)               ! original axis is now -z-axis
             case default
               call IO_error(157_pInt,section)
           end select
         enddo
       
       case ('hybridia') textureType
         texture_ODFfile(section) = IO_stringValue(line,positions,2_pInt)

       case ('symmetry') textureType
         tag = IO_lc(IO_stringValue(line,positions,2_pInt))
         select case (tag)
           case('orthotropic')
             texture_symmetry(section) = 4_pInt
           case('monoclinic')
             texture_symmetry(section) = 2_pInt
           case default
             texture_symmetry(section) = 1_pInt
         end select
         
       case ('(random)') textureType
         gauss = gauss + 1_pInt
         texture_Gauss(1:3,gauss,section) = math_sampleRandomOri()
         do j = 2_pInt,4_pInt,2_pInt
           tag = IO_lc(IO_stringValue(line,positions,j))
           select case (tag)
             case('scatter')
                 texture_Gauss(4,gauss,section) = IO_floatValue(line,positions,j+1_pInt)*inRad
             case('fraction')
                 texture_Gauss(5,gauss,section) = IO_floatValue(line,positions,j+1_pInt)
           end select
         enddo

       case ('(gauss)') textureType
         gauss = gauss + 1_pInt
         do j = 2_pInt,10_pInt,2_pInt
           tag = IO_lc(IO_stringValue(line,positions,j))
           select case (tag)
             case('phi1')
                 texture_Gauss(1,gauss,section) = IO_floatValue(line,positions,j+1_pInt)*inRad
             case('phi')
                 texture_Gauss(2,gauss,section) = IO_floatValue(line,positions,j+1_pInt)*inRad
             case('phi2')
                 texture_Gauss(3,gauss,section) = IO_floatValue(line,positions,j+1_pInt)*inRad
             case('scatter')
                 texture_Gauss(4,gauss,section) = IO_floatValue(line,positions,j+1_pInt)*inRad
             case('fraction')
                 texture_Gauss(5,gauss,section) = IO_floatValue(line,positions,j+1_pInt)
           end select
         enddo

       case ('(fiber)') textureType
         fiber = fiber + 1_pInt
         do j = 2_pInt,12_pInt,2_pInt
           tag = IO_lc(IO_stringValue(line,positions,j))
           select case (tag)
             case('alpha1')
                 texture_Fiber(1,fiber,section) = IO_floatValue(line,positions,j+1_pInt)*inRad
             case('alpha2')
                 texture_Fiber(2,fiber,section) = IO_floatValue(line,positions,j+1_pInt)*inRad
             case('beta1')
                 texture_Fiber(3,fiber,section) = IO_floatValue(line,positions,j+1_pInt)*inRad
             case('beta2')
                 texture_Fiber(4,fiber,section) = IO_floatValue(line,positions,j+1_pInt)*inRad
             case('scatter')
                 texture_Fiber(5,fiber,section) = IO_floatValue(line,positions,j+1_pInt)*inRad
             case('fraction')
                 texture_Fiber(6,fiber,section) = IO_floatValue(line,positions,j+1_pInt)
           end select
         enddo

     end select textureType
   endif
 enddo

end subroutine material_parseTexture


!--------------------------------------------------------------------------------------------------
!> @brief populates the grains
!> @details populates the grains by identifying active microstructure/homogenization pairs,
!! calculates the volume of the grains and deals with texture components and hybridIA
!--------------------------------------------------------------------------------------------------
subroutine material_populateGrains
 use math, only: &
   math_RtoEuler, &
   math_EulerToR, &
   math_mul33x33, &
   math_range, &
   math_sampleRandomOri, &
   math_sampleGaussOri, &
   math_sampleFiberOri, &
   math_symmetricEulers
 use mesh, only: &
   mesh_element, &
   mesh_maxNips, &
   mesh_NcpElems, &
   mesh_ipVolume, &
   FE_Nips, &
   FE_geomtype
 use IO, only: &
   IO_error, &
   IO_hybridIA
 use FEsolving, only: &
   FEsolving_execIP
 use debug, only: &
   debug_level, &
   debug_material, &
   debug_levelBasic
 
 implicit none
 integer(pInt), dimension (:,:), allocatable :: Ngrains
 integer(pInt), dimension (microstructure_maxNconstituents)  :: &
   NgrainsOfConstituent, &
   currentGrainOfConstituent, &
   randomOrder
 real(pReal), dimension (microstructure_maxNconstituents)  :: &
   rndArray
 real(pReal), dimension (:),     allocatable :: volumeOfGrain
 real(pReal), dimension (:,:),   allocatable :: orientationOfGrain
 real(pReal), dimension (3)                  :: orientation
 real(pReal), dimension (3,3)                :: symOrientation
 integer(pInt), dimension (:),   allocatable :: phaseOfGrain, textureOfGrain
 integer(pInt) :: t,e,i,g,j,m,c,r,homog,micro,sgn,hme, myDebug, &
                  phaseID,textureID,dGrains,myNgrains,myNorientations,myNconstituents, &
                  grain,constituentGrain,ipGrain,symExtension, ip
 real(pReal) :: extreme,rnd
 integer(pInt), dimension (:,:),   allocatable :: Nelems                                            ! counts number of elements in homog, micro array
 type(p_intvec), dimension (:,:), allocatable :: elemsOfHomogMicro                                  ! lists element number in homog, micro array

 myDebug = debug_level(debug_material)
 
 allocate(material_volume(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ;        material_volume      = 0.0_pReal
 allocate(material_phase(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ;         material_phase       = 0_pInt
 allocate(material_texture(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ;       material_texture     = 0_pInt
 allocate(material_EulerAngles(3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; material_EulerAngles = 0.0_pReal
 
 allocate(Ngrains(material_Nhomogenization,material_Nmicrostructure)); Ngrains = 0_pInt
 allocate(Nelems(material_Nhomogenization,material_Nmicrostructure));  Nelems = 0_pInt
 
!--------------------------------------------------------------------------------------------------
! precounting of elements for each homog/micro pair
 do e = 1_pInt, mesh_NcpElems
   homog = mesh_element(3,e)
   micro = mesh_element(4,e)
   Nelems(homog,micro) = Nelems(homog,micro) + 1_pInt
 enddo
 allocate(elemsOfHomogMicro(material_Nhomogenization,material_Nmicrostructure))
 do homog = 1,material_Nhomogenization
   do micro = 1,material_Nmicrostructure
     if (Nelems(homog,micro) > 0_pInt) then
       allocate(elemsOfHomogMicro(homog,micro)%p(Nelems(homog,micro)))
       elemsOfHomogMicro(homog,micro)%p = 0_pInt
    endif
   enddo
 enddo

!--------------------------------------------------------------------------------------------------
! identify maximum grain count per IP (from element) and find grains per homog/micro pair
 Nelems = 0_pInt                                                                                    ! reuse as counter
 elementLooping: do e = 1_pInt,mesh_NcpElems
   t    = FE_geomtype(mesh_element(2,e))
   homog = mesh_element(3,e)
   micro = mesh_element(4,e)
   if (homog < 1_pInt .or. homog > material_Nhomogenization) &                                      ! out of bounds
     call IO_error(154_pInt,e,0_pInt,0_pInt)
   if (micro < 1_pInt .or. micro > material_Nmicrostructure) &                                      ! out of bounds
     call IO_error(155_pInt,e,0_pInt,0_pInt)
   if (microstructure_elemhomo(micro)) then                                                         ! how many grains are needed at this element?
     dGrains = homogenization_Ngrains(homog)                                                        ! only one set of Ngrains (other IPs are plain copies)
   else
     dGrains = homogenization_Ngrains(homog) * FE_Nips(t)                                           ! each IP has Ngrains
   endif
   Ngrains(homog,micro) = Ngrains(homog,micro) + dGrains                                            ! total grain count
   Nelems(homog,micro)  = Nelems(homog,micro) + 1_pInt                                              ! total element count
   elemsOfHomogMicro(homog,micro)%p(Nelems(homog,micro)) = e                                        ! remember elements active in this homog/micro pair
 enddo elementLooping
   
 allocate(volumeOfGrain(maxval(Ngrains)))                                                           ! reserve memory for maximum case
 allocate(phaseOfGrain(maxval(Ngrains)))                                                            ! reserve memory for maximum case
 allocate(textureOfGrain(maxval(Ngrains)))                                                          ! reserve memory for maximum case
 allocate(orientationOfGrain(3,maxval(Ngrains)))                                                    ! reserve memory for maximum case
 
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
     write(6,'(/,a/)') ' MATERIAL grain population'
     write(6,'(a32,1x,a32,1x,a6)') 'homogenization_name','microstructure_name','grain#'
   !$OMP END CRITICAL (write2out)
 endif
 do homog = 1_pInt,material_Nhomogenization                                                         ! loop over homogenizations
   dGrains = homogenization_Ngrains(homog)                                                          ! grain number per material point
   do micro = 1_pInt,material_Nmicrostructure                                                       ! all pairs of homog and micro
     if (Ngrains(homog,micro) > 0_pInt) then                                                        ! an active pair of homog and micro
       myNgrains = Ngrains(homog,micro)                                                             ! assign short name for total number of grains to populate
       myNconstituents = microstructure_Nconstituents(micro)                                        ! assign short name for number of constituents
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
         !$OMP CRITICAL (write2out)
           write(6,'(/,a32,1x,a32,1x,i6)') homogenization_name(homog),microstructure_name(micro),myNgrains
         !$OMP END CRITICAL (write2out)
       endif


!--------------------------------------------------------------------------------------------------
! calculate volume of each grain

       volumeOfGrain = 0.0_pReal
       grain = 0_pInt

       do hme = 1_pInt, Nelems(homog,micro)
         e = elemsOfHomogMicro(homog,micro)%p(hme)                                                  ! my combination of homog and micro, only perform calculations for elements with homog, micro combinations which is indexed in cpElemsindex
         t = FE_geomtype(mesh_element(2,e))
         if (microstructure_elemhomo(micro)) then                                                   ! homogeneous distribution of grains over each element's IPs
           volumeOfGrain(grain+1_pInt:grain+dGrains) = sum(mesh_ipVolume(1:FE_Nips(t),e))/&
                                                                         real(dGrains,pReal)        ! each grain combines size of all IPs in that element
           grain = grain + dGrains                                                                  ! wind forward by Ngrains@IP
         else
           forall (i = 1_pInt:FE_Nips(t)) &                                                         ! loop over IPs
             volumeOfGrain(grain+(i-1)*dGrains+1_pInt:grain+i*dGrains) = &
               mesh_ipVolume(i,e)/dGrains                                                           ! assign IPvolume/Ngrains@IP to all grains of IP
           grain = grain + FE_Nips(t) * dGrains                                                     ! wind forward by Nips*Ngrains@IP
         endif
       enddo

       if (grain /= myNgrains) &
         call IO_error(0,e = homog,i = micro,ext_msg = 'inconsistent grain count after volume calc')

!--------------------------------------------------------------------------------------------------
! divide myNgrains as best over constituents
!
! example: three constituents with fractions of 0.25, 0.25, and 0.5 distributed over 20 (microstructure) grains
!
!                       ***** ***** **********
! NgrainsOfConstituent: 5,    5,    10
! counters:
!                      |-----> grain (if constituent == 2)
!                            |--> constituentGrain (of constituent 2)
!

       NgrainsOfConstituent = 0_pInt                                                                ! reset counter of grains per constituent
       forall (i = 1_pInt:myNconstituents) &
         NgrainsOfConstituent(i) = nint(microstructure_fraction(i,micro) * myNgrains, pInt)         ! do rounding integer conversion
       do while (sum(NgrainsOfConstituent) /= myNgrains)                                            ! total grain count over constituents wrong?
         sgn = sign(1_pInt, myNgrains - sum(NgrainsOfConstituent))                                  ! direction of required change
         extreme = 0.0_pReal
         t = 0_pInt
         do i = 1_pInt,myNconstituents                                                              ! find largest deviator
           if (real(sgn,pReal)*log(NgrainsOfConstituent(i)/myNgrains/microstructure_fraction(i,micro)) > extreme) then
             extreme = real(sgn,pReal)*log(NgrainsOfConstituent(i)/myNgrains/microstructure_fraction(i,micro))
             t = i
           endif
         enddo
         NgrainsOfConstituent(t) = NgrainsOfConstituent(t) + sgn                                    ! change that by one
       enddo

!--------------------------------------------------------------------------------------------------
! assign phase and texture info

       phaseOfGrain = 0_pInt
       textureOfGrain = 0_pInt
       orientationOfGrain = 0.0_pReal

       texture: do i = 1_pInt,myNconstituents                                                       ! loop over constituents
         grain            = sum(NgrainsOfConstituent(1_pInt:i-1_pInt))                              ! set microstructure grain index of current constituent
                                                                                                    ! "grain" points to start of this constituent's grain population
         constituentGrain = 0_pInt                                                                  ! constituent grain index

         phaseID   = microstructure_phase(i,micro)
         textureID = microstructure_texture(i,micro)
         phaseOfGrain  (grain+1_pInt:grain+NgrainsOfConstituent(i)) = phaseID                       ! assign resp. phase
         textureOfGrain(grain+1_pInt:grain+NgrainsOfConstituent(i)) = textureID                     ! assign resp. texture

         myNorientations = ceiling(real(NgrainsOfConstituent(i),pReal)/&
                                   real(texture_symmetry(textureID),pReal),pInt)                    ! max number of unique orientations (excl. symmetry)

!--------------------------------------------------------------------------------------------------
! ...has texture components
         if (texture_ODFfile(textureID) == '') then
           gauss: do t = 1_pInt,texture_Ngauss(textureID)                                           ! loop over Gauss components
             do g = 1_pInt,int(myNorientations*texture_Gauss(5,t,textureID),pInt)                   ! loop over required grain count
               orientationOfGrain(:,grain+constituentGrain+g) = &
                 math_sampleGaussOri(texture_Gauss(1:3,t,textureID),&
                                     texture_Gauss(  4,t,textureID))
             enddo
             constituentGrain = &
             constituentGrain + int(myNorientations*texture_Gauss(5,t,textureID))                   ! advance counter for grains of current constituent
           enddo gauss

           fiber: do t = 1_pInt,texture_Nfiber(textureID)                                           ! loop over fiber components
             do g = 1_pInt,int(myNorientations*texture_Fiber(6,t,textureID),pInt)                   ! loop over required grain count
               orientationOfGrain(:,grain+constituentGrain+g) = &
                 math_sampleFiberOri(texture_Fiber(1:2,t,textureID),&
                                     texture_Fiber(3:4,t,textureID),&
                                     texture_Fiber(  5,t,textureID))
             enddo
             constituentGrain = &
             constituentGrain + int(myNorientations*texture_fiber(6,t,textureID),pInt)              ! advance counter for grains of current constituent
           enddo fiber

           random: do constituentGrain = constituentGrain+1_pInt,myNorientations                    ! fill remainder with random
              orientationOfGrain(:,grain+constituentGrain) = math_sampleRandomOri()
           enddo random
!--------------------------------------------------------------------------------------------------
! ...has hybrid IA
         else
           orientationOfGrain(1:3,grain+1_pInt:grain+myNorientations) = &
                                            IO_hybridIA(myNorientations,texture_ODFfile(textureID))
           if (all(orientationOfGrain(1:3,grain+1_pInt) == -1.0_pReal)) call IO_error(156_pInt)
         endif

!--------------------------------------------------------------------------------------------------
! ...texture rotation

         do j = 1_pInt,myNorientations                                                              ! loop over each "real" orientation
           orientationOfGrain(1:3,grain+j) = math_RtoEuler( &                                       ! translate back to Euler angles
                                             math_mul33x33( &                                       ! pre-multiply
                                               math_EulertoR(orientationOfGrain(1:3,grain+j)), &    ! face-value orientation
                                               texture_rotation(1:3,1:3,textureID) &               ! rotation matrix and
                                             ) &
                                             )
         enddo

!--------------------------------------------------------------------------------------------------
! ...sample symmetry

         symExtension = texture_symmetry(textureID) - 1_pInt
         if (symExtension > 0_pInt) then                                                            ! sample symmetry (number of additional equivalent orientations)
           constituentGrain = myNorientations                                                       ! start right after "real" orientations
           do j = 1_pInt,myNorientations                                                            ! loop over each "real" orientation
             symOrientation = math_symmetricEulers(texture_symmetry(textureID), &
                                                   orientationOfGrain(1:3,grain+j))                 ! get symmetric equivalents
             e = min(symExtension,NgrainsOfConstituent(i)-constituentGrain)                         ! do not overshoot end of constituent grain array
             if (e > 0_pInt) then
               orientationOfGrain(1:3,grain+constituentGrain+1:   &
                                      grain+constituentGrain+e) = &
                 symOrientation(1:3,1:e)
               constituentGrain = constituentGrain + e                                              ! remainder shrinks by e
             endif
           enddo
         endif

!--------------------------------------------------------------------------------------------------
! shuffle grains within current constituent

         do j = 1_pInt,NgrainsOfConstituent(i)-1_pInt                                               ! walk thru grains of current constituent
           call random_number(rnd)
           t = nint(rnd*(NgrainsOfConstituent(i)-j)+j+0.5_pReal,pInt)                               ! select a grain in remaining list
           m                               = phaseOfGrain(grain+t)                                  ! exchange current with random
           phaseOfGrain(grain+t)           = phaseOfGrain(grain+j)
           phaseOfGrain(grain+j)           = m
           m                               = textureOfGrain(grain+t)                                ! exchange current with random
           textureOfGrain(grain+t)         = textureOfGrain(grain+j)
           textureOfGrain(grain+j)         = m
           orientation                     = orientationOfGrain(1:3,grain+t)                        ! exchange current with random
           orientationOfGrain(1:3,grain+t) = orientationOfGrain(1:3,grain+j)
           orientationOfGrain(1:3,grain+j) = orientation
         enddo
         
       enddo texture
!< @todo calc fraction after weighing with volumePerGrain, exchange in MC steps to improve result (humbug at the moment)

 

!--------------------------------------------------------------------------------------------------
! distribute grains of all constituents as accurately as possible to given constituent fractions

       ip = 0_pInt
       currentGrainOfConstituent = 0_pInt

       do hme = 1_pInt, Nelems(homog,micro)
         e = elemsOfHomogMicro(homog,micro)%p(hme)                                                  ! only perform calculations for elements with homog, micro combinations which is indexed in cpElemsindex
         t = FE_geomtype(mesh_element(2,e))
         if (microstructure_elemhomo(micro)) then                                                   ! homogeneous distribution of grains over each element's IPs
           m = 1_pInt                                                                               ! process only first IP
         else
           m = FE_Nips(t)                                                                           ! process all IPs
         endif

         do i = 1_pInt, m                                                                           ! loop over necessary IPs
           ip = ip + 1_pInt                                                                         ! keep track of total ip count
           ipGrain = 0_pInt                                                                         ! count number of grains assigned at this IP
           randomOrder = math_range(microstructure_maxNconstituents)                                ! start out with ordered sequence of constituents
           call random_number(rndArray)                                                             ! as many rnd numbers as (max) constituents
           do j = 1_pInt, myNconstituents - 1_pInt                                                  ! loop over constituents ...
             r = nint(rndArray(j)*(myNconstituents-j)+j+0.5_pReal,pInt)                             ! ... select one in remaining list
             c = randomOrder(r)                                                                     ! ... call it "c"
             randomOrder(r) = randomOrder(j)                                                        ! ... and exchange with present position in constituent list
             grain = sum(NgrainsOfConstituent(1:c-1_pInt))                                          ! figure out actual starting index in overall/consecutive grain population
             do g = 1_pInt, min(dGrains-ipGrain, &                                                  ! leftover number of grains at this IP
                                max(0_pInt, &                                                       ! no negative values
                                    nint(real(ip * dGrains * NgrainsOfConstituent(c)) / &           ! fraction of grains scaled to this constituent...
                                         real(myNgrains),pInt) - &                                  ! ...minus those already distributed
                                         currentGrainOfConstituent(c)))
               ipGrain = ipGrain + 1_pInt                                                           ! advance IP grain counter
               currentGrainOfConstituent(c)  = currentGrainOfConstituent(c) + 1_pInt                ! advance index of grain population for constituent c
               material_volume(ipGrain,i,e)  = volumeOfGrain(grain+currentGrainOfConstituent(c))    ! assign properties
               material_phase(ipGrain,i,e)   = phaseOfGrain(grain+currentGrainOfConstituent(c))
               material_texture(ipGrain,i,e) = textureOfGrain(grain+currentGrainOfConstituent(c))
               material_EulerAngles(1:3,ipGrain,i,e) = orientationOfGrain(1:3,grain+currentGrainOfConstituent(c))
           enddo; enddo

           c = randomOrder(microstructure_Nconstituents(micro))                                     ! look up constituent remaining after random shuffling
           grain = sum(NgrainsOfConstituent(1:c-1_pInt))                                            ! figure out actual starting index in overall/consecutive grain population
           do ipGrain = ipGrain + 1_pInt, dGrains                                                   ! ensure last constituent fills up to dGrains
             currentGrainOfConstituent(c)  = currentGrainOfConstituent(c) + 1_pInt
             material_volume(ipGrain,i,e)  = volumeOfGrain(grain+currentGrainOfConstituent(c))
             material_phase(ipGrain,i,e)   = phaseOfGrain(grain+currentGrainOfConstituent(c))
             material_texture(ipGrain,i,e) = textureOfGrain(grain+currentGrainOfConstituent(c))
             material_EulerAngles(1:3,ipGrain,i,e) = orientationOfGrain(1:3,grain+currentGrainOfConstituent(c))
           enddo

         enddo

         do i = i, FE_Nips(t)                                                                       ! loop over IPs to (possibly) distribute copies from first IP
           material_volume (1_pInt:dGrains,i,e) = material_volume (1_pInt:dGrains,1,e)
           material_phase  (1_pInt:dGrains,i,e) = material_phase  (1_pInt:dGrains,1,e)
           material_texture(1_pInt:dGrains,i,e) = material_texture(1_pInt:dGrains,1,e)
           material_EulerAngles(1:3,1_pInt:dGrains,i,e) = material_EulerAngles(1:3,1_pInt:dGrains,1,e)
         enddo

       enddo
     endif                                                                                          ! active homog,micro pair
   enddo
 enddo
 
 deallocate(volumeOfGrain)
 deallocate(phaseOfGrain)
 deallocate(textureOfGrain)
 deallocate(orientationOfGrain)
 deallocate(Nelems)
 !> @todo - causing segmentation fault: needs looking into
 !do homog = 1,material_Nhomogenization
 !  do micro = 1,material_Nmicrostructure
 !    if (Nelems(homog,micro) > 0_pInt) deallocate(elemsOfHomogMicro(homog,micro)%p)
 !  enddo
 !enddo
 deallocate(elemsOfHomogMicro)

end subroutine material_populateGrains

end module material
