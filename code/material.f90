!* $Id$
!************************************
!*      Module: MATERIAL            *
!************************************
!* contains:                        *
!* - parsing of material.config     *
!************************************

MODULE material

!*** Include other modules ***
use prec, only: pReal,pInt
implicit none

character(len=64), parameter :: material_configFile = 'material.config'
character(len=32), parameter :: material_partHomogenization = 'homogenization'
character(len=32), parameter :: material_partMicrostructure = 'microstructure'
character(len=32), parameter :: material_partPhase =          'phase'
character(len=32), parameter :: material_partTexture =        'texture'


!*************************************
!* Definition of material properties *
!*************************************
!* Number of materials
integer(pInt) material_Nhomogenization, &                                              ! number of homogenizations
              material_Nmicrostructure, &                                              ! number of microstructures
              material_Nphase, &                                                       ! number of phases
              material_Ntexture, &                                                     ! number of textures
              microstructure_maxNconstituents, &                                       ! max number of constituents in any phase
              homogenization_maxNgrains, &                                             ! max number of grains in any USED homogenization
              texture_maxNgauss, &                                                     ! max number of Gauss components in any texture
              texture_maxNfiber                                                        ! max number of Fiber components in any texture
character(len=64), dimension(:),       allocatable :: homogenization_name, &           ! name of each homogenization
                                                      homogenization_type, &           ! type of each homogenization
                                                      microstructure_name, &           ! name of each microstructure
                                                      phase_name, &                    ! name of each phase
                                                      phase_constitution, &            ! constitution of each phase
                                                      texture_name                     ! name of each texture
character(len=256),dimension(:),       allocatable :: texture_ODFfile                  ! name of each ODF file
integer(pInt),     dimension(:),       allocatable :: homogenization_Ngrains, &        ! number of grains in each homogenization
                                                      homogenization_typeInstance, &   ! instance of particular type of each homogenization
                                                      homogenization_Noutput, &        ! number of '(output)' items per homogenization
                                                      microstructure_Nconstituents, &  ! number of constituents in each microstructure
                                                      phase_constitutionInstance, &    ! instance of particular constitution of each phase
                                                      phase_Noutput, &                 ! number of '(output)' items per phase
                                                      texture_symmetry, &              ! number of symmetric orientations per texture
                                                      texture_Ngauss, &                ! number of Gauss components per texture
                                                      texture_Nfiber                   ! number of Fiber components per texture
logical,           dimension(:),       allocatable :: homogenization_active, &         !
                                                      microstructure_active, &         ! 
                                                      microstructure_elemhomo, &       ! flag to indicate homogeneous microstructure distribution over element's IPs
                                                      phase_localConstitution          ! flags phases with local constitutive law
integer(pInt),     dimension(:,:),     allocatable :: microstructure_phase, &          ! phase IDs of each microstructure
                                                      microstructure_texture           ! texture IDs of each microstructure
real(pReal),       dimension(:,:),     allocatable :: microstructure_fraction          ! vol fraction of each constituent in microstructure
real(pReal),       dimension(:,:,:),   allocatable :: material_volume                  ! volume of each grain,IP,element
integer(pInt),     dimension(:,:,:),   allocatable :: material_phase                   ! phase of each grain,IP,element
real(pReal),       dimension(:,:,:,:), allocatable :: material_EulerAngles             ! initial orientation of each grain,IP,element
real(pReal),       dimension(:,:,:),   allocatable :: texture_Gauss, &                 ! data of each Gauss component
                                                      texture_Fiber                    ! data of each Fiber component

CONTAINS


!*********************************************************************
subroutine material_init()
!*********************************************************************
!*      Module initialization         *
!**************************************
 use prec, only: pReal,pInt
 use IO, only: IO_error, IO_open_file
 implicit none

!* Definition of variables
 integer(pInt), parameter :: fileunit = 200
 integer(pInt) i,j
 
 write(6,*)
 write(6,*) '<<<+-  material init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)
  
 if(.not. IO_open_file(fileunit,material_configFile)) call IO_error (100) ! corrupt config file
 call material_parseHomogenization(fileunit,material_partHomogenization)
 call material_parseMicrostructure(fileunit,material_partMicrostructure)
 call material_parseTexture(fileunit,material_partTexture)
 call material_parsePhase(fileunit,material_partPhase)
 close(fileunit)

 do i = 1,material_Nmicrostructure
   if (minval(microstructure_phase(1:microstructure_Nconstituents(i),i)) < 1 .or. &
       maxval(microstructure_phase(1:microstructure_Nconstituents(i),i)) > material_Nphase) call IO_error(150,i)
   if (minval(microstructure_texture(1:microstructure_Nconstituents(i),i)) < 1 .or. &
       maxval(microstructure_texture(1:microstructure_Nconstituents(i),i)) > material_Ntexture) call IO_error(160,i)
   if (abs(sum(microstructure_fraction(:,i)) - 1.0_pReal) >= 1.0e-10_pReal) then
     write(6,*)'sum of microstructure fraction = ',sum(microstructure_fraction(:,i))
     call IO_error(170,i)
   endif
 enddo
 write (6,*)
 write (6,*) 'MATERIAL configuration'
 write (6,*)
 write (6,'(a32,x,a16,x,a6)') 'homogenization                  ','type            ','grains'
 do i = 1,material_Nhomogenization
   write (6,'(x,a32,x,a16,x,i4)') homogenization_name(i),homogenization_type(i),homogenization_Ngrains(i)
 enddo
 write (6,*)
 write (6,'(a32,x,a12,x,a13)') 'microstructure                  ','constituents','homogeneous'
 do i = 1,material_Nmicrostructure
   write (6,'(a32,4x,i4,8x,l)') microstructure_name(i),microstructure_Nconstituents(i),microstructure_elemhomo(i)
   if (microstructure_Nconstituents(i) > 0_pInt) then
     do j = 1,microstructure_Nconstituents(i)
       write (6,'(a1,x,a32,x,a32,x,f6.4)') '>',phase_name(microstructure_phase(j,i)),&
                                               texture_name(microstructure_texture(j,i)),&
                                               microstructure_fraction(j,i)
     enddo
     write (6,*)
   endif
 enddo

 call material_populateGrains()

endsubroutine


!*********************************************************************
subroutine material_parseHomogenization(file,myPart)
!*********************************************************************

 use prec, only: pInt
 use IO
 use mesh, only: mesh_element
 implicit none

 character(len=*), intent(in) :: myPart
 integer(pInt), intent(in) :: file
 integer(pInt), parameter :: maxNchunks = 2
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) Nsections, section, s
 character(len=64) tag
 character(len=1024) line
 
 Nsections = IO_countSections(file,myPart)
 material_Nhomogenization = Nsections
 
 allocate(homogenization_name(Nsections));    homogenization_name = ''
 allocate(homogenization_type(Nsections));    homogenization_type = ''
 allocate(homogenization_typeInstance(Nsections));  homogenization_typeInstance = 0_pInt
 allocate(homogenization_Ngrains(Nsections)); homogenization_Ngrains = 0_pInt
 allocate(homogenization_Noutput(Nsections)); homogenization_Noutput = 0_pInt
 allocate(homogenization_active(Nsections));  homogenization_active = .false.

 forall (s = 1:Nsections) homogenization_active(s) = any(mesh_element(3,:) == s)    ! current homogenization used in model?
 homogenization_Noutput = IO_countTagInPart(file,myPart,'(output)',Nsections)
 
 rewind(file)
 line = ''
 section = 0
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= myPart)      ! wind forward to myPart
   read(file,'(a1024)',END=100) line
 enddo

 do
   read(file,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                ! next section
     section = section + 1
     homogenization_name(section) = IO_getTag(line,'[',']')
   endif
   if (section > 0) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     select case(tag)
       case ('type')
         homogenization_type(section) = IO_lc(IO_stringValue(line,positions,2))  ! adding: IO_lc function <<<updated 31.07.2009>>>
         do s = 1,section
           if (homogenization_type(s) == homogenization_type(section)) &
             homogenization_typeInstance(section) = homogenization_typeInstance(section) + 1  ! count instances
         enddo
       case ('ngrains')
         homogenization_Ngrains(section) = IO_intValue(line,positions,2)
     end select
   endif
 enddo

100 homogenization_maxNgrains = maxval(homogenization_Ngrains,homogenization_active)
  return

 endsubroutine


!*********************************************************************
subroutine material_parseMicrostructure(file,myPart)
!*********************************************************************

 use prec, only: pInt
 use IO
 use mesh, only: mesh_element
 implicit none

 character(len=*), intent(in) :: myPart
 integer(pInt), intent(in) :: file
 integer(pInt), parameter :: maxNchunks = 7
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) Nsections, section, constituent, i
 character(len=64) tag
 character(len=1024) line
 
 Nsections = IO_countSections(file,myPart)
 material_Nmicrostructure = Nsections
 allocate(microstructure_name(Nsections));     microstructure_name = ''
 allocate(microstructure_Nconstituents(Nsections))
 allocate(microstructure_active(Nsections))
 allocate(microstructure_elemhomo(Nsections))

 forall (i = 1:Nsections) microstructure_active(i) = any(mesh_element(4,:) == i)    ! current microstructure used in model?

 microstructure_Nconstituents = IO_countTagInPart(file,myPart,'(constituent)',Nsections)
 microstructure_maxNconstituents = maxval(microstructure_Nconstituents)
 microstructure_elemhomo = IO_spotTagInPart(file,myPart,'/elementhomogeneous/',Nsections)

 allocate(microstructure_phase   (microstructure_maxNconstituents,Nsections)); microstructure_phase    = 0_pInt
 allocate(microstructure_texture (microstructure_maxNconstituents,Nsections)); microstructure_texture  = 0_pInt
 allocate(microstructure_fraction(microstructure_maxNconstituents,Nsections)); microstructure_fraction = 0.0_pReal
 
 rewind(file)
 line = ''
 section = 0
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= myPart)      ! wind forward to myPart
   read(file,'(a1024)',END=100) line
 enddo

 do
   read(file,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                ! next section
     section = section + 1
     constituent = 0
     microstructure_name(section) = IO_getTag(line,'[',']')
   endif
   if (section > 0) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     select case(tag)
       case ('(constituent)')
         constituent = constituent + 1
         do i=2,6,2
           tag = IO_lc(IO_stringValue(line,positions,i))
           select case (tag)
             case('phase')
               microstructure_phase(constituent,section) =    IO_intValue(line,positions,i+1)
             case('texture')
               microstructure_texture(constituent,section) =  IO_intValue(line,positions,i+1)
             case('fraction')
               microstructure_fraction(constituent,section) = IO_floatValue(line,positions,i+1)
           end select
         enddo
     end select
   endif
 enddo

100 return

 endsubroutine


!*********************************************************************
subroutine material_parsePhase(file,myPart)
!*********************************************************************

 use prec, only: pInt
 use IO
 implicit none

 character(len=*), intent(in) :: myPart
 integer(pInt), intent(in) :: file
 integer(pInt), parameter :: maxNchunks = 2
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) Nsections, section, s
 character(len=64) tag
 character(len=1024) line
 
 Nsections = IO_countSections(file,myPart)
 material_Nphase = Nsections
 allocate(phase_name(Nsections));          phase_name = ''
 allocate(phase_constitution(Nsections));  phase_constitution = ''
 allocate(phase_constitutionInstance(Nsections));  phase_constitutionInstance = 0_pInt
 allocate(phase_Noutput(Nsections))
 allocate(phase_localConstitution(Nsections))

 phase_Noutput = IO_countTagInPart(file,myPart,'(output)',Nsections)
 phase_localConstitution = .not. IO_spotTagInPart(file,myPart,'/nonlocal/',Nsections)
 
 rewind(file)
 line = ''
 section = 0
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= myPart)      ! wind forward to myPart
   read(file,'(a1024)',END=100) line
 enddo

 do
   read(file,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                ! next section
     section = section + 1
     phase_name(section) = IO_getTag(line,'[',']')
   endif
   if (section > 0) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     select case(tag)
       case ('constitution')
         phase_constitution(section) = IO_lc(IO_stringValue(line,positions,2))
         do s = 1,section
           if (phase_constitution(s) == phase_constitution(section)) &
             phase_constitutionInstance(section) = phase_constitutionInstance(section) + 1  ! count instances
         enddo
     end select
   endif
 enddo

100 return

 endsubroutine


!*********************************************************************
subroutine material_parseTexture(file,myPart)
!*********************************************************************

 use prec, only: pInt, pReal
 use IO
 use math, only: inRad
 implicit none

 character(len=*), intent(in) :: myPart
 integer(pInt), intent(in) :: file
 integer(pInt), parameter :: maxNchunks = 13
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) Nsections, section, gauss, fiber, i
 character(len=64) tag
 character(len=1024) line
 
 
 Nsections = IO_countSections(file,myPart)
 material_Ntexture = Nsections
 allocate(texture_name(Nsections));     texture_name = ''
 allocate(texture_ODFfile(Nsections));  texture_ODFfile = ''
 allocate(texture_symmetry(Nsections)); texture_symmetry = 1_pInt
 allocate(texture_Ngauss(Nsections));   texture_Ngauss = 0_pInt
 allocate(texture_Nfiber(Nsections));   texture_Nfiber = 0_pInt

 texture_Ngauss = IO_countTagInPart(file,myPart,'(gauss)',Nsections)
 texture_Nfiber = IO_countTagInPart(file,myPart,'(fiber)',Nsections)
 texture_maxNgauss = maxval(texture_Ngauss)
 texture_maxNfiber = maxval(texture_Nfiber)
 allocate(texture_Gauss   (5,texture_maxNgauss,Nsections)); texture_Gauss    = 0.0_pReal
 allocate(texture_Fiber   (6,texture_maxNfiber,Nsections)); texture_Fiber    = 0.0_pReal
 
 rewind(file)
 line = ''
 section = 0
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= myPart)      ! wind forward to myPart
   read(file,'(a1024)',END=100) line
 enddo

 do
   read(file,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                ! next section
     section = section + 1
     gauss = 0
     fiber = 0
     texture_name(section) = IO_getTag(line,'[',']')
   endif
   if (section > 0) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     select case(tag)

       case ('hybridia')
         texture_ODFfile(section) = IO_stringValue(line,positions,2)

       case ('symmetry')
         tag = IO_lc(IO_stringValue(line,positions,2))
         select case (tag)
           case('orthotropic')
             texture_symmetry(section) = 4
           case('monoclinic')
             texture_symmetry(section) = 2
           case default
             texture_symmetry(section) = 1
         end select
         
       case ('(gauss)')
         gauss = gauss + 1
         do i = 2,10,2
           tag = IO_lc(IO_stringValue(line,positions,i))
           select case (tag)
             case('phi1')
                 texture_Gauss(1,gauss,section) = IO_floatValue(line,positions,i+1)*inRad
             case('phi')
                 texture_Gauss(2,gauss,section) = IO_floatValue(line,positions,i+1)*inRad
             case('phi2')
                 texture_Gauss(3,gauss,section) = IO_floatValue(line,positions,i+1)*inRad
             case('scatter')
                 texture_Gauss(4,gauss,section) = IO_floatValue(line,positions,i+1)*inRad
             case('fraction')
                 texture_Gauss(5,gauss,section) = IO_floatValue(line,positions,i+1)
           end select
         enddo
       case ('(fiber)')
         fiber = fiber + 1
         do i = 2,12,2
           tag = IO_lc(IO_stringValue(line,positions,i))
           select case (tag)
             case('alpha1')
                 texture_Fiber(1,fiber,section) = IO_floatValue(line,positions,i+1)*inRad
             case('alpha2')
                 texture_Fiber(2,fiber,section) = IO_floatValue(line,positions,i+1)*inRad
             case('beta1')
                 texture_Fiber(3,fiber,section) = IO_floatValue(line,positions,i+1)*inRad
             case('beta2')
                 texture_Fiber(4,fiber,section) = IO_floatValue(line,positions,i+1)*inRad
             case('scatter')
                 texture_Fiber(5,fiber,section) = IO_floatValue(line,positions,i+1)*inRad
             case('fraction')
                 texture_Fiber(6,fiber,section) = IO_floatValue(line,positions,i+1)
           end select
         enddo

     end select
   endif
 enddo

100 return

 endsubroutine


!*********************************************************************
subroutine material_populateGrains()
!*********************************************************************

 use prec, only: pInt, pReal
 use math, only: math_sampleRandomOri, math_sampleGaussOri, math_sampleFiberOri, math_symmetricEulers
 use mesh, only: mesh_element, mesh_maxNips, mesh_NcpElems, mesh_ipVolume, FE_Nips
 use IO,   only: IO_error, IO_hybridIA
 implicit none

 integer(pInt), dimension (:,:), allocatable :: Ngrains
 integer(pInt), dimension (microstructure_maxNconstituents) :: NgrainsOfConstituent
 real(pReal), dimension (:),     allocatable :: volumeOfGrain, phaseOfGrain
 real(pReal), dimension (:,:),   allocatable :: orientationOfGrain
 real(pReal), dimension (3) :: orientation
 real(pReal), dimension (3,3) :: symOrientation
 integer(pInt) t,e,i,g,j,m,homog,micro,sgn
 integer(pInt) phaseID,textureID,dGrains,myNgrains,myNorientations, &
               grain,constituentGrain,symExtension
 real(pReal) extreme,rnd

 allocate(material_volume(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; material_volume = 0.0_pReal
 allocate(material_phase(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; material_phase = 0_pInt
 allocate(material_EulerAngles(3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; material_EulerAngles = 0.0_pReal
 
 allocate(Ngrains(material_Nhomogenization,material_Nmicrostructure)); Ngrains = 0_pInt

! identify maximum grain count per IP (from element) and find grains per homog/micro pair
 do e = 1,mesh_NcpElems
   homog = mesh_element(3,e)
   micro = mesh_element(4,e)
   if (homog < 1 .or. homog > material_Nhomogenization) &   ! out of bounds
     call IO_error(130,e,0,0)
   if (micro < 1 .or. micro > material_Nmicrostructure) &   ! out of bounds
     call IO_error(140,e,0,0)
   if (microstructure_elemhomo(micro)) then
     dGrains = homogenization_Ngrains(homog)
   else
     dGrains = homogenization_Ngrains(homog) * FE_Nips(mesh_element(2,e))
   endif
   Ngrains(homog,micro) = Ngrains(homog,micro) + dGrains
 enddo

 allocate(volumeOfGrain(maxval(Ngrains)))           ! reserve memory for maximum case
 allocate(phaseOfGrain(maxval(Ngrains)))            ! reserve memory for maximum case
 allocate(orientationOfGrain(3,maxval(Ngrains)))    ! reserve memory for maximum case
 
 write (6,*)
 write (6,*) 'MATERIAL grain population'
 write (6,*)
 write (6,'(a32,x,a32,x,a6)') 'homogenization_name','microstructure_name','grain#'
 do homog = 1,material_Nhomogenization              ! loop over homogenizations
   dGrains = homogenization_Ngrains(homog)          ! grain number per material point
   do micro = 1,material_Nmicrostructure            ! all pairs of homog and micro
     if (Ngrains(homog,micro) > 0) then             ! an active pair of homog and micro
       myNgrains = Ngrains(homog,micro)             ! assign short name for total number of grains to populate
       write (6,*)
       write (6,'(a32,x,a32,x,i6)') homogenization_name(homog),microstructure_name(micro),myNgrains
     
! ----------------------------------------------------------------------------  calculate volume of each grain
       volumeOfGrain = 0.0_pReal
       grain = 0_pInt                               ! microstructure grain index
       do e = 1,mesh_NcpElems                       ! check each element
         if (mesh_element(3,e) == homog .and. mesh_element(4,e) == micro) then  ! my combination of homog and micro
           if (microstructure_elemhomo(micro)) then                             ! homogeneous distribution of grains over each element's IPs
             volumeOfGrain(grain+1:grain+dGrains) = sum(mesh_ipVolume(1:FE_Nips(mesh_element(2,e)),e))/dGrains
             grain = grain + dGrains                                            ! wind forward by NgrainsPerIP
           else
             forall (i = 1:FE_Nips(mesh_element(2,e))) &                        ! loop over IPs
               volumeOfGrain(grain+(i-1)*dGrains+1:grain+i*dGrains) = &
                 mesh_ipVolume(i,e)/dGrains                                     ! assign IPvolume/Ngrains to all grains of IP
             grain = grain + FE_Nips(mesh_element(2,e)) * dGrains               ! wind forward by Nips*NgrainsPerIP
           endif
         endif
       enddo
! ----------------------------------------------------------------------------  divide myNgrains as best over constituents
       NgrainsOfConstituent = 0_pInt
       forall (i = 1:microstructure_Nconstituents(micro)) &
         NgrainsOfConstituent(i) = nint(microstructure_fraction(i,micro) * myNgrains, pInt)  ! do rounding integer conversion
       do while (sum(NgrainsOfConstituent) /= myNgrains)                        ! total grain count over constituents wrong?
         sgn = sign(1_pInt, myNgrains - sum(NgrainsOfConstituent))              ! direction of required change
         extreme = 0.0_pReal
         t = 0_pInt
         do i = 1,microstructure_Nconstituents(micro)                           ! find largest deviator
           if (sgn*log(NgrainsOfConstituent(i)/myNgrains/microstructure_fraction(i,micro)) > extreme) then
             extreme = sgn*log(NgrainsOfConstituent(i)/myNgrains/microstructure_fraction(i,micro))
             t = i
           endif
         enddo
         NgrainsOfConstituent(t) = NgrainsOfConstituent(t) + sgn               ! change that by one
       enddo
! ----------------------------------------------------------------------------
       phaseOfGrain = 0_pInt
       orientationOfGrain = 0.0_pReal
       grain = 0_pInt                                                         ! reset microstructure grain index

       do i = 1,microstructure_Nconstituents(micro)                           ! loop over constituents
         phaseID   = microstructure_phase(i,micro)
         textureID = microstructure_texture(i,micro)
         phaseOfGrain(grain+1:grain+NgrainsOfConstituent(i)) = phaseID        ! assign resp. phase

         myNorientations = ceiling(float(NgrainsOfConstituent(i))/texture_symmetry(textureID))   ! max number of unique orientations (excl. symmetry)

         constituentGrain = 0_pInt                                            ! constituent grain index
                                                                              ! ---------
         if (texture_ODFfile(textureID) == '') then                           ! dealing with texture components
                                                                              ! ---------
           do t = 1,texture_Ngauss(textureID)                                 ! loop over Gauss components
             do g = 1,int(myNorientations*texture_Gauss(5,t,textureID))       ! loop over required grain count
               orientationOfGrain(:,grain+constituentGrain+g) = &
                 math_sampleGaussOri(texture_Gauss(1:3,t,textureID),&
                                     texture_Gauss(  4,t,textureID))
             enddo
             constituentGrain = constituentGrain + int(myNorientations*texture_Gauss(5,t,textureID))
           enddo

           do t = 1,texture_Nfiber(textureID)                                 ! loop over fiber components
             do g = 1,int(myNorientations*texture_Fiber(6,t,textureID))       ! loop over required grain count
               orientationOfGrain(:,grain+constituentGrain+g) = &
                 math_sampleFiberOri(texture_Fiber(1:2,t,textureID),&
                                     texture_Fiber(3:4,t,textureID),&
                                     texture_Fiber(  5,t,textureID))
             enddo
             constituentGrain = constituentGrain + int(myNorientations*texture_fiber(6,t,textureID))
           enddo

           do j = constituentGrain+1,myNorientations                          ! fill remainder with random
              orientationOfGrain(:,grain+j) = math_sampleRandomOri()
           enddo
                                                                              ! ---------
         else                                                                 ! hybrid IA
                                                                              ! ---------
           orientationOfGrain(:,grain+1:grain+myNorientations) = IO_hybridIA(myNorientations,texture_ODFfile(textureID))
           if (all(orientationOfGrain(:,grain+1) == -1.0_pReal)) call IO_error(105)  
           constituentGrain = constituentGrain + myNorientations

         endif
! ----------------------------------------------------------------------------
         symExtension = texture_symmetry(textureID) - 1_pInt
         if (symExtension > 0_pInt) then                                      ! sample symmetry
           constituentGrain = NgrainsOfConstituent(i)-myNorientations         ! calc remainder of array
           do j = 1,myNorientations                                           ! loop over each "real" orientation
             symOrientation = math_symmetricEulers(texture_symmetry(textureID),orientationOfGrain(:,j))  ! get symmetric equivalents
             e = min(symExtension,constituentGrain)                           ! are we at end of constituent grain array?
             if (e > 0_pInt) then
               orientationOfGrain(:,grain+myNorientations+1+(j-1)*symExtension:&
                                    grain+myNorientations+e+(j-1)*symExtension) = &
                 symOrientation(:,1:e)
               constituentGrain = constituentGrain - e                        ! remainder shrinks by e
             endif
           enddo
         endif
         
         grain = grain + NgrainsOfConstituent(i)                              ! advance microstructure grain index
       enddo  ! constituent

! ----------------------------------------------------------------------------
       do i=1,myNgrains-1                                                     ! walk thru grains
         call random_number(rnd)
         t = nint(rnd*(myNgrains-i)+i+0.5_pReal,pInt)                         ! select a grain in remaining list
         m                       = phaseOfGrain(t)                            ! exchange current with random
         phaseOfGrain(t)         = phaseOfGrain(i)
         phaseOfGrain(i)         = m
         orientation             = orientationOfGrain(:,t)
         orientationOfGrain(:,t) = orientationOfGrain(:,i)
         orientationOfGrain(:,i) = orientation
       enddo
       !calc fraction after weighing with volumePerGrain
       !exchange in MC steps to improve result...

! ----------------------------------------------------------------------------
       grain = 0_pInt                                                           ! microstructure grain index
       do e = 1,mesh_NcpElems                                                   ! check each element
         if (mesh_element(3,e) == homog .and. mesh_element(4,e) == micro) then  ! my combination of homog and micro
           if (microstructure_elemhomo(micro)) then                             ! homogeneous distribution of grains over each element's IPs
             forall (i = 1:FE_Nips(mesh_element(2,e)), g = 1:dGrains)           ! loop over IPs and grains
               material_volume(g,i,e) = volumeOfGrain(grain+g)
               material_phase(g,i,e) = phaseOfGrain(grain+g)
               material_EulerAngles(:,g,i,e) = orientationOfGrain(:,grain+g)
             end forall
             grain = grain + dGrains                                            ! wind forward by NgrainsPerIP
           else
             forall (i = 1:FE_Nips(mesh_element(2,e)), g = 1:dGrains)           ! loop over IPs and grains
               material_volume(g,i,e) = volumeOfGrain(grain+(i-1)*dGrains+g)
               material_phase(g,i,e) = phaseOfGrain(grain+(i-1)*dGrains+g)
               material_EulerAngles(:,g,i,e) = orientationOfGrain(:,grain+(i-1)*dGrains+g)
             end forall
             grain = grain + FE_Nips(mesh_element(2,e)) * dGrains               ! wind forward by Nips*NgrainsPerIP
           endif
           ! write (6,*) e
           ! write (6,*) material_phase(:,:,e)
           ! write (6,*) material_EulerAngles(:,:,:,e)
         endif
       enddo

     endif   ! active homog,micro pair
   enddo
 enddo

 
 deallocate(volumeOfGrain)
 deallocate(phaseOfGrain)
 deallocate(orientationOfGrain)
 
 return

 endsubroutine


END MODULE
