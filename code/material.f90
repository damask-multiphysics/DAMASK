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
!! Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
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
   material_configFile         = 'material.config', &                                               !< generic name for material configuration file 
   material_localFileExt       = 'materialConfig'                                                   !< extension of solver job name depending material configuration file  
   
 character(len=32), parameter, public  :: &
   material_partHomogenization = 'homogenization', &                                                !< keyword for homogenization part
   material_partCrystallite    = 'crystallite', &                                                   !< keyword for crystallite part
   material_partPhase          = 'phase'                                                            !< keyword for phase part
 
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
   material_partMicrostructure = 'microstructure', &                                                !< keyword for microstructure part
   material_partTexture        = 'texture'                                                          !< keyword for texture part
   
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
   texture_Fiber                                                                                    !< data of each Fiber component
 
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
 integer(pInt)            :: i,j, myDebug
 
 myDebug = debug_level(debug_material)
 
 write(6,'(/,a)') ' <<<+-  material init  -+>>>'
 write(6,'(a)') ' $Id$'
 write(6,'(a16,a)')   ' Current time : ',IO_timeStamp()
#include "compilation_info.f90"

 
 if (.not. IO_open_jobFile_stat(fileunit,material_localFileExt)) then        ! no local material configuration present...
   call IO_open_file(fileunit,material_configFile)                           ! ...open material.config file
 endif
 call material_parseHomogenization(fileunit,material_partHomogenization)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
   write(6,*) 'Homogenization parsed'
 endif
 call material_parseMicrostructure(fileunit,material_partMicrostructure)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
   write(6,*) 'Microstructure parsed'
 endif
 call material_parseCrystallite(fileunit,material_partCrystallite)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
   write(6,*) 'Crystallite parsed'
 endif
 call material_parseTexture(fileunit,material_partTexture)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
   write(6,*) 'Texture parsed'
 endif
 call material_parsePhase(fileunit,material_partPhase)
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
   write(6,*) 'Phase parsed'
 endif
 close(fileunit)

 do i = 1_pInt,material_Nmicrostructure
   if (microstructure_crystallite(i) < 1_pInt .or. &
       microstructure_crystallite(i) > material_Ncrystallite) call IO_error(150_pInt,i)
   if (minval(microstructure_phase(1:microstructure_Nconstituents(i),i)) < 1_pInt .or. &
       maxval(microstructure_phase(1:microstructure_Nconstituents(i),i)) > material_Nphase) call IO_error(151_pInt,i)
   if (minval(microstructure_texture(1:microstructure_Nconstituents(i),i)) < 1_pInt .or. &
       maxval(microstructure_texture(1:microstructure_Nconstituents(i),i)) > material_Ntexture) call IO_error(152_pInt,i)
   if (abs(sum(microstructure_fraction(:,i)) - 1.0_pReal) >= 1.0e-10_pReal) then
     if (iand(myDebug,debug_levelExtensive) /= 0_pInt) then
         write(6,*)'sum of microstructure fraction = ',sum(microstructure_fraction(:,i))
     endif
     call IO_error(153_pInt,i)
   endif
 enddo
 if (iand(myDebug,debug_levelExtensive) /= 0_pInt) then
   write(6,'(/,a,/)') ' MATERIAL configuration'
   write(6,'(a32,1x,a16,1x,a6)') 'homogenization                  ','type            ','grains'
   do i = 1_pInt,material_Nhomogenization
     write(6,'(1x,a32,1x,a16,1x,i4)') homogenization_name(i),homogenization_type(i),homogenization_Ngrains(i)
   enddo
   write(6,*)
   write(6,'(a32,1x,a11,1x,a12,1x,a13)') 'microstructure                  ','crystallite','constituents','homogeneous'
   do i = 1_pInt,material_Nmicrostructure
     write(6,'(a32,4x,i4,8x,i4,8x,l1)') microstructure_name(i), &
                                        microstructure_crystallite(i), &
                                        microstructure_Nconstituents(i), &
                                        microstructure_elemhomo(i)
     if (microstructure_Nconstituents(i) > 0_pInt) then
       do j = 1_pInt,microstructure_Nconstituents(i)
         write(6,'(a1,1x,a32,1x,a32,1x,f7.4)') '>',phase_name(microstructure_phase(j,i)),&
                                                   texture_name(microstructure_texture(j,i)),&
                                                   microstructure_fraction(j,i)
       enddo
       write(6,*)
     endif
   enddo
 endif
 
 call material_populateGrains

end subroutine material_init


!--------------------------------------------------------------------------------------------------
!> @brief parses the homogenization part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseHomogenization(myFile,myPart)
 use IO
 use mesh, only: &
   mesh_element
 
 implicit none
 character(len=*), intent(in) :: myPart
 integer(pInt),    intent(in) :: myFile
 
 integer(pInt),     parameter :: maxNchunks = 2_pInt
 
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) Nsections, section, s
 character(len=64)   :: tag
 character(len=1024) :: line
 logical             :: echo
 
 echo = IO_globalTagInPart(myFile,myPart,'/echo/')
 
 Nsections = IO_countSections(myFile,myPart)
 material_Nhomogenization = Nsections
 if (Nsections < 1_pInt) call IO_error(160_pInt,ext_msg=myPart)
 
 allocate(homogenization_name(Nsections));    homogenization_name = ''
 allocate(homogenization_type(Nsections));    homogenization_type = ''
 allocate(homogenization_typeInstance(Nsections));  homogenization_typeInstance = 0_pInt
 allocate(homogenization_Ngrains(Nsections)); homogenization_Ngrains = 0_pInt
 allocate(homogenization_Noutput(Nsections)); homogenization_Noutput = 0_pInt
 allocate(homogenization_active(Nsections));  homogenization_active = .false.

 forall (s = 1_pInt:Nsections) homogenization_active(s) = any(mesh_element(3,:) == s)               ! current homogenization used in model? Homogenization view, maximum operations depend on maximum number of homog schemes
 homogenization_Noutput = IO_countTagInPart(myFile,myPart,'(output)',Nsections)
 
 rewind(myFile)
 line = ''
 section = 0_pInt
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= myPart)                                                ! wind forward to myPart
   read(myFile,'(a1024)',END=100) line
 enddo
 if (echo) write(6,*) trim(line)                                                                    ! echo part header

 do
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (echo) write(6,*) trim(line)                                                                  ! echo back read lines
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

100 homogenization_maxNgrains = maxval(homogenization_Ngrains,homogenization_active)

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
 character(len=64)   :: tag
 character(len=1024) :: line
 logical             :: echo

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
 line = ''
 section = 0_pInt
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= myPart)                                                ! wind forward to myPart
   read(myFile,'(a1024)',END=100) line
 enddo
 if (echo) write(6,*) trim(line)                                                                    ! echo part header

 do
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (echo) write(6,*) trim(line)                                                                  ! echo back read lines
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

100 end subroutine material_parseMicrostructure


!--------------------------------------------------------------------------------------------------
!> @brief parses the crystallite part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseCrystallite(myFile,myPart)
 use IO, only: &
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
 
 integer(pInt)       :: Nsections, &
                        section
 character(len=1024) :: line
 logical             :: echo

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
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= myPart)                                                ! wind forward to myPart
   read(myFile,'(a1024)',END=100) line
 enddo
 if (echo) write(6,*) trim(line)                                                                    ! echo part header

 do
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (echo) write(6,*) trim(line)                                                                  ! echo back read lines
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     crystallite_name(section) = IO_getTag(line,'[',']')
   endif
 enddo

100 end subroutine material_parseCrystallite


!--------------------------------------------------------------------------------------------------
!> @brief parses the phase part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parsePhase(myFile,myPart)
 use IO
 
 implicit none
 character(len=*), intent(in) :: myPart
 integer(pInt),    intent(in) :: myFile
 
 integer(pInt), parameter :: maxNchunks = 2_pInt
 
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) Nsections, section, s
 character(len=64)   :: tag
 character(len=1024) :: line
 logical             :: echo

 echo = IO_globalTagInPart(myFile,myPart,'/echo/')
 
 Nsections = IO_countSections(myFile,myPart)
 material_Nphase = Nsections
 if (Nsections < 1_pInt) call IO_error(160_pInt,ext_msg=myPart)

 allocate(phase_name(Nsections));          phase_name = ''
 allocate(phase_elasticity(Nsections));  phase_elasticity = ''
 allocate(phase_elasticityInstance(Nsections));  phase_elasticityInstance = 0_pInt
 allocate(phase_plasticity(Nsections));  phase_plasticity = ''
 allocate(phase_plasticityInstance(Nsections));  phase_plasticityInstance = 0_pInt
 allocate(phase_Noutput(Nsections))
 allocate(phase_localPlasticity(Nsections))

 phase_Noutput = IO_countTagInPart(myFile,myPart,'(output)',Nsections)
 phase_localPlasticity = .not. IO_spotTagInPart(myFile,myPart,'/nonlocal/',Nsections)
 
 rewind(myFile)
 line = ''
 section = 0_pInt
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= myPart)                                                ! wind forward to myPart
   read(myFile,'(a1024)',END=100) line
 enddo
 if (echo) write(6,*) trim(line)                                                                    ! echo part header

 do
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (echo) write(6,*) trim(line)                                                                  ! echo back read lines
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

100 end subroutine material_parsePhase


!--------------------------------------------------------------------------------------------------
!> @brief parses the texture part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseTexture(myFile,myPart)
 use IO
 use math, only: &
   inRad, &
   math_sampleRandomOri
 
 implicit none
 character(len=*), intent(in) :: myPart
 integer(pInt),    intent(in) :: myFile
 
 integer(pInt), parameter     :: maxNchunks = 13_pInt
 
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) :: Nsections, section, gauss, fiber, i
 character(len=64) :: tag
 character(len=1024) :: line
 logical             :: echo

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
 
 rewind(myFile)
 line = ''
 section = 0_pInt
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= myPart)                                                ! wind forward to myPart
   read(myFile,'(a1024)',END=100) line
 enddo
 if (echo) write(6,*) trim(line)                                                                    ! echo part header

 do
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (echo) write(6,*) trim(line)                                                                  ! echo back read lines
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     gauss = 0_pInt
     fiber = 0_pInt
     texture_name(section) = IO_getTag(line,'[',']')
   endif
   if (section > 0_pInt) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)

       case ('hybridia')
         texture_ODFfile(section) = IO_stringValue(line,positions,2_pInt)

       case ('symmetry')
         tag = IO_lc(IO_stringValue(line,positions,2_pInt))
         select case (tag)
           case('orthotropic')
             texture_symmetry(section) = 4_pInt
           case('monoclinic')
             texture_symmetry(section) = 2_pInt
           case default
             texture_symmetry(section) = 1_pInt
         end select
         
       case ('(random)')
         gauss = gauss + 1_pInt
         texture_Gauss(1:3,gauss,section) = math_sampleRandomOri()
         do i = 2_pInt,4_pInt,2_pInt
           tag = IO_lc(IO_stringValue(line,positions,i))
           select case (tag)
             case('scatter')
                 texture_Gauss(4,gauss,section) = IO_floatValue(line,positions,i+1_pInt)*inRad
             case('fraction')
                 texture_Gauss(5,gauss,section) = IO_floatValue(line,positions,i+1_pInt)
           end select
         enddo

       case ('(gauss)')
         gauss = gauss + 1_pInt
         do i = 2_pInt,10_pInt,2_pInt
           tag = IO_lc(IO_stringValue(line,positions,i))
           select case (tag)
             case('phi1')
                 texture_Gauss(1,gauss,section) = IO_floatValue(line,positions,i+1_pInt)*inRad
             case('phi')
                 texture_Gauss(2,gauss,section) = IO_floatValue(line,positions,i+1_pInt)*inRad
             case('phi2')
                 texture_Gauss(3,gauss,section) = IO_floatValue(line,positions,i+1_pInt)*inRad
             case('scatter')
                 texture_Gauss(4,gauss,section) = IO_floatValue(line,positions,i+1_pInt)*inRad
             case('fraction')
                 texture_Gauss(5,gauss,section) = IO_floatValue(line,positions,i+1_pInt)
           end select
         enddo

       case ('(fiber)')
         fiber = fiber + 1_pInt
         do i = 2_pInt,12_pInt,2_pInt
           tag = IO_lc(IO_stringValue(line,positions,i))
           select case (tag)
             case('alpha1')
                 texture_Fiber(1,fiber,section) = IO_floatValue(line,positions,i+1_pInt)*inRad
             case('alpha2')
                 texture_Fiber(2,fiber,section) = IO_floatValue(line,positions,i+1_pInt)*inRad
             case('beta1')
                 texture_Fiber(3,fiber,section) = IO_floatValue(line,positions,i+1_pInt)*inRad
             case('beta2')
                 texture_Fiber(4,fiber,section) = IO_floatValue(line,positions,i+1_pInt)*inRad
             case('scatter')
                 texture_Fiber(5,fiber,section) = IO_floatValue(line,positions,i+1_pInt)*inRad
             case('fraction')
                 texture_Fiber(6,fiber,section) = IO_floatValue(line,positions,i+1_pInt)
           end select
         enddo

     end select
   endif
 enddo

100 end subroutine material_parseTexture


!--------------------------------------------------------------------------------------------------
!> @brief populates the grains
!> @details populates the grains by identifying active microstructure/homogenization pairs,
!! calculates the volume of the grains and deals with texture components and hybridIA
!--------------------------------------------------------------------------------------------------
subroutine material_populateGrains
 use math, only: &
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
 integer(pInt), dimension (microstructure_maxNconstituents) &
                                             :: NgrainsOfConstituent
 real(pReal), dimension (:),     allocatable :: volumeOfGrain
 real(pReal), dimension (:,:),   allocatable :: orientationOfGrain
 real(pReal), dimension (3)                  :: orientation
 real(pReal), dimension (3,3)                :: symOrientation
 integer(pInt), dimension (:),   allocatable :: phaseOfGrain, textureOfGrain
 integer(pInt) :: t,e,i,g,j,m,homog,micro,sgn,hme, myDebug
 integer(pInt) :: phaseID,textureID,dGrains,myNgrains,myNorientations, &
               grain,constituentGrain,symExtension
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
 do e = 1_pInt,mesh_NcpElems
   t     = FE_geomtype(mesh_element(2,e))
   homog = mesh_element(3,e)
   micro = mesh_element(4,e)
   if (homog < 1_pInt .or. homog > material_Nhomogenization) &                                      ! out of bounds
     call IO_error(154_pInt,e,0_pInt,0_pInt)
   if (micro < 1_pInt .or. micro > material_Nmicrostructure) &                                      ! out of bounds
     call IO_error(155_pInt,e,0_pInt,0_pInt)
   if (microstructure_elemhomo(micro)) then
     dGrains = homogenization_Ngrains(homog)
   else
     dGrains = homogenization_Ngrains(homog) * FE_Nips(t)
   endif
   Ngrains(homog,micro) = Ngrains(homog,micro) + dGrains
   Nelems(homog,micro)  = Nelems(homog,micro) + 1_pInt
   elemsOfHomogMicro(homog,micro)%p(Nelems(homog,micro)) = e                                        ! remember elements active in this homog/micro pair
   
 enddo
 allocate(volumeOfGrain(maxval(Ngrains)))                                                           ! reserve memory for maximum case
 allocate(phaseOfGrain(maxval(Ngrains)))                                                            ! reserve memory for maximum case
 allocate(textureOfGrain(maxval(Ngrains)))                                                          ! reserve memory for maximum case
 allocate(orientationOfGrain(3,maxval(Ngrains)))                                                    ! reserve memory for maximum case
 
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
     write(6,*)
     write(6,*) 'MATERIAL grain population'
     write(6,*)
     write(6,'(a32,1x,a32,1x,a6)') 'homogenization_name','microstructure_name','grain#'
   !$OMP END CRITICAL (write2out)
 endif
 do homog = 1_pInt,material_Nhomogenization                                                         ! loop over homogenizations
   dGrains = homogenization_Ngrains(homog)                                                          ! grain number per material point
   do micro = 1_pInt,material_Nmicrostructure                                                       ! all pairs of homog and micro
     if (Ngrains(homog,micro) > 0_pInt) then                                                        ! an active pair of homog and micro
       myNgrains = Ngrains(homog,micro)                                                             ! assign short name for total number of grains to populate
       if (iand(myDebug,debug_levelBasic) /= 0_pInt) then
         !$OMP CRITICAL (write2out)
           write(6,*)
           write(6,'(a32,1x,a32,1x,i6)') homogenization_name(homog),microstructure_name(micro),myNgrains
         !$OMP END CRITICAL (write2out)
       endif
     
!--------------------------------------------------------------------------------------------------
! calculate volume of each grain
       volumeOfGrain = 0.0_pReal
       grain = 0_pInt
       do hme = 1_pInt, Nelems(homog,micro)
         e = elemsOfHomogMicro(homog,micro)%p(hme)                                                     ! my combination of homog and micro, only perform calculations for elements with homog, micro combinations which is indexed in cpElemsindex
         t = FE_geomtype(mesh_element(2,e))
         if (microstructure_elemhomo(micro)) then                                                   ! homogeneous distribution of grains over each element's IPs
           volumeOfGrain(grain+1_pInt:grain+dGrains) = sum(mesh_ipVolume(1:FE_Nips(t),e))/&
                                                                         real(dGrains,pReal)
           grain = grain + dGrains                                                                  ! wind forward by NgrainsPerIP
         else
           forall (i = 1_pInt:FE_Nips(t)) &                                                         ! loop over IPs
             volumeOfGrain(grain+(i-1)*dGrains+1_pInt:grain+i*dGrains) = &
               mesh_ipVolume(i,e)/dGrains                                                           ! assign IPvolume/Ngrains to all grains of IP
           grain = grain + FE_Nips(t) * dGrains                                                     ! wind forward by Nips*NgrainsPerIP
         endif
       enddo
       
!--------------------------------------------------------------------------------------------------
! divide myNgrains as best over constituents
       NgrainsOfConstituent = 0_pInt
       forall (i = 1_pInt:microstructure_Nconstituents(micro)) &
         NgrainsOfConstituent(i) = nint(microstructure_fraction(i,micro) * myNgrains, pInt)        ! do rounding integer conversion
       do while (sum(NgrainsOfConstituent) /= myNgrains)                                           ! total grain count over constituents wrong?
         sgn = sign(1_pInt, myNgrains - sum(NgrainsOfConstituent))                                 ! direction of required change
         extreme = 0.0_pReal
         t = 0_pInt
         do i = 1_pInt,microstructure_Nconstituents(micro)                                         ! find largest deviator
           if (real(sgn,pReal)*log(NgrainsOfConstituent(i)/myNgrains/microstructure_fraction(i,micro)) > extreme) then
             extreme = real(sgn,pReal)*log(NgrainsOfConstituent(i)/myNgrains/microstructure_fraction(i,micro))
             t = i
           endif
         enddo
         NgrainsOfConstituent(t) = NgrainsOfConstituent(t) + sgn                                   ! change that by one
       enddo

       phaseOfGrain = 0_pInt
       textureOfGrain = 0_pInt
       orientationOfGrain = 0.0_pReal
       grain = 0_pInt                                                                              ! reset microstructure grain index

       constituents: do i = 1_pInt,microstructure_Nconstituents(micro)                             ! loop over constituents
         phaseID   = microstructure_phase(i,micro)
         textureID = microstructure_texture(i,micro)
         phaseOfGrain(grain+1_pInt:grain+NgrainsOfConstituent(i)) = phaseID                        ! assign resp. phase
         textureOfGrain(grain+1_pInt:grain+NgrainsOfConstituent(i)) = textureID                    ! assign resp. texture

         myNorientations = ceiling(real(NgrainsOfConstituent(i),pReal)/&
                                   real(texture_symmetry(textureID),pReal),pInt)                   ! max number of unique orientations (excl. symmetry)

         constituentGrain = 0_pInt                                                                 ! constituent grain index
!--------------------------------------------------------------------------------------------------
! dealing with texture components
         if (texture_ODFfile(textureID) == '') then
           do t = 1_pInt,texture_Ngauss(textureID)                                                 ! loop over Gauss components
             do g = 1_pInt,int(myNorientations*texture_Gauss(5,t,textureID),pInt)                  ! loop over required grain count
               orientationOfGrain(:,grain+constituentGrain+g) = &
                 math_sampleGaussOri(texture_Gauss(1:3,t,textureID),&
                                     texture_Gauss(  4,t,textureID))
             enddo
             constituentGrain = constituentGrain + int(myNorientations*texture_Gauss(5,t,textureID))
           enddo

           do t = 1_pInt,texture_Nfiber(textureID)                                                 ! loop over fiber components
             do g = 1_pInt,int(myNorientations*texture_Fiber(6,t,textureID),pInt)                  ! loop over required grain count
               orientationOfGrain(:,grain+constituentGrain+g) = &
                 math_sampleFiberOri(texture_Fiber(1:2,t,textureID),&
                                     texture_Fiber(3:4,t,textureID),&
                                     texture_Fiber(  5,t,textureID))
             enddo
             constituentGrain = constituentGrain + int(myNorientations*texture_fiber(6,t,textureID),pInt)
           enddo

           do j = constituentGrain+1_pInt,myNorientations                                           ! fill remainder with random
              orientationOfGrain(:,grain+j) = math_sampleRandomOri()
           enddo
         else
!--------------------------------------------------------------------------------------------------
! hybrid IA
           orientationOfGrain(:,grain+1:grain+myNorientations) = IO_hybridIA(myNorientations,texture_ODFfile(textureID))
           if (all(orientationOfGrain(:,grain+1) == -1.0_pReal)) call IO_error(156_pInt)  
           constituentGrain = constituentGrain + myNorientations
         endif

         symExtension = texture_symmetry(textureID) - 1_pInt
         if (symExtension > 0_pInt) then                                                            ! sample symmetry
           constituentGrain = NgrainsOfConstituent(i)-myNorientations                               ! calc remainder of array
           do j = 1_pInt,myNorientations                                                            ! loop over each "real" orientation
             symOrientation = math_symmetricEulers(texture_symmetry(textureID),orientationOfGrain(:,j))  ! get symmetric equivalents
             e = min(symExtension,constituentGrain)                                                 ! are we at end of constituent grain array?
             if (e > 0_pInt) then
               orientationOfGrain(:,grain+myNorientations+1+(j-1_pInt)*symExtension:&
                                    grain+myNorientations+e+(j-1_pInt)*symExtension) = &
                 symOrientation(:,1:e)
               constituentGrain = constituentGrain - e                                              ! remainder shrinks by e
             endif
           enddo
         endif
         
         grain = grain + NgrainsOfConstituent(i)                                                    ! advance microstructure grain index
       enddo constituents


! unless element homogeneous, reshuffle grains
       if (.not. microstructure_elemhomo(micro)) then                                               
         do i=1_pInt,myNgrains-1_pInt                                                               ! walk thru grains
           call random_number(rnd)
           t = nint(rnd*(myNgrains-i)+i+0.5_pReal,pInt)                                             ! select a grain in remaining list
           m                       = phaseOfGrain(t)                                                ! exchange current with random
           phaseOfGrain(t)         = phaseOfGrain(i)
           phaseOfGrain(i)         = m
           m                       = textureOfGrain(t)                                              ! exchange current with random
           textureOfGrain(t)       = textureOfGrain(i)
           textureOfGrain(i)       = m
           orientation             = orientationOfGrain(:,t)
           orientationOfGrain(:,t) = orientationOfGrain(:,i)
           orientationOfGrain(:,i) = orientation
         enddo
       endif
       
!--------------------------------------------------------------------------------------------------
! calc fraction after weighing with volumePerGrain, exchange in MC steps to improve result...
       grain = 0_pInt
       do hme = 1_pInt, Nelems(homog,micro)
         e = elemsOfHomogMicro(homog,micro)%p(hme)                                                  ! only perform calculations for elements with homog, micro combinations which is indexed in cpElemsindex
         t = FE_geomtype(mesh_element(2,e))
         if (microstructure_elemhomo(micro)) then                                                   ! homogeneous distribution of grains over each element's IPs
           forall (i = 1_pInt:FE_Nips(t), g = 1_pInt:dGrains)                                       ! loop over IPs and grains
             material_volume(g,i,e)        = volumeOfGrain(grain+g)
             material_phase(g,i,e)         = phaseOfGrain(grain+g)
             material_texture(g,i,e)       = textureOfGrain(grain+g)
             material_EulerAngles(:,g,i,e) = orientationOfGrain(:,grain+g)
           end forall
           FEsolving_execIP(2,e) = 1_pInt                                                           ! restrict calculation to first IP only, since all other results are to be copied from this
           grain = grain + dGrains                                                                  ! wind forward by NgrainsPerIP
         else
           forall (i = 1_pInt:FE_Nips(t), g = 1_pInt:dGrains)                                       ! loop over IPs and grains
             material_volume(g,i,e)        = volumeOfGrain(grain+(i-1_pInt)*dGrains+g)
             material_phase(g,i,e)         = phaseOfGrain(grain+(i-1_pInt)*dGrains+g)
             material_texture(g,i,e)       = textureOfGrain(grain+(i-1_pInt)*dGrains+g)
             material_EulerAngles(:,g,i,e) = orientationOfGrain(:,grain+(i-1_pInt)*dGrains+g)
           end forall
           grain = grain + FE_Nips(t) * dGrains                                                     ! wind forward by Nips*NgrainsPerIP
         endif
       enddo
     endif                                                                                          ! active homog,micro pair
   enddo
 enddo
 
 deallocate(volumeOfGrain)
 deallocate(phaseOfGrain)
 deallocate(textureOfGrain)
 deallocate(orientationOfGrain)
 deallocate(Nelems)
 !do homog = 1,material_Nhomogenization
 !  do micro = 1,material_Nmicrostructure
 !    if (Nelems(homog,micro) > 0_pInt) deallocate(elemsOfHomogMicro(homog,micro)%p)                ! ToDo - causing segmentation fault: needs looking into
 !  enddo
 !enddo
 deallocate(elemsOfHomogMicro)

end subroutine material_populateGrains

end module material
