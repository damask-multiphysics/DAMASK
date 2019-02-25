!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  contains lattice structure definitions including Schmid matrices for slip, twin, trans,
!          and cleavage as well as interaction among the various systems
!--------------------------------------------------------------------------------------------------
module lattice
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private

! BEGIN DEPRECATED
 integer(pInt), parameter, public :: &
   LATTICE_maxNslipFamily     =  13_pInt, &                                                         !< max # of slip system families over lattice structures
   LATTICE_maxNcleavageFamily =  3_pInt                                                             !< max # of transformation system families over lattice structures

 integer(pInt), allocatable, dimension(:,:), protected, public :: &
   lattice_NslipSystem, &                                                                           !< total # of slip systems in each family
   lattice_NcleavageSystem                                                                          !< total # of transformation systems in each family

 real(pReal), allocatable, dimension(:,:,:,:,:), protected, public :: &
   lattice_Scleavage                                                                                !< Schmid matrices for cleavage systems

 real(pReal), allocatable, dimension(:,:,:), protected, public :: &
   lattice_sn, &                                                                                    !< normal direction of slip system
   lattice_st, &                                                                                    !< sd x sn
   lattice_sd                                                                                       !< slip direction of slip system
! END DEPRECATED


!--------------------------------------------------------------------------------------------------
! face centered cubic
 integer(pInt), dimension(LATTICE_maxNslipFamily), parameter, private :: &
   LATTICE_FCC_NSLIPSYSTEM = int([12, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],pInt)                     !< # of slip systems per family for fcc

 integer(pInt), dimension(1), parameter, private :: &
   LATTICE_FCC_NTWINSYSTEM = int([12],pInt)                                                         !< # of twin systems per family for fcc

 integer(pInt), dimension(1), parameter, private :: &
   LATTICE_FCC_NTRANSSYSTEM = int([12],pInt)                                                        !< # of transformation systems per family for fcc

 integer(pInt), dimension(LATTICE_maxNcleavageFamily), parameter, private :: &
   LATTICE_FCC_NCLEAVAGESYSTEM = int([3, 4, 0],pInt)                                                !< # of cleavage systems per family for fcc

 integer(pInt), parameter, private  :: &
   LATTICE_FCC_NSLIP      = sum(LATTICE_FCC_NSLIPSYSTEM), &                                         !< total # of slip systems for fcc
   LATTICE_FCC_NTWIN      = sum(LATTICE_FCC_NTWINSYSTEM), &                                         !< total # of twin systems for fcc
   LATTICE_FCC_NTRANS     = sum(LATTICE_FCC_NTRANSSYSTEM), &                                        !< total # of transformation systems for fcc
   LATTICE_FCC_NCLEAVAGE  = sum(LATTICE_FCC_NCLEAVAGESYSTEM)                                        !< total # of cleavage systems for fcc

 real(pReal), dimension(3+3,LATTICE_FCC_NSLIP), parameter, private :: &
   LATTICE_FCC_SYSTEMSLIP = reshape(real([&
    ! Slip direction     Plane normal                                                               ! SCHMID-BOAS notation
      0, 1,-1,     1, 1, 1, &                                                                       ! B2
     -1, 0, 1,     1, 1, 1, &                                                                       ! B4
      1,-1, 0,     1, 1, 1, &                                                                       ! B5
      0,-1,-1,    -1,-1, 1, &                                                                       ! C1
      1, 0, 1,    -1,-1, 1, &                                                                       ! C3
     -1, 1, 0,    -1,-1, 1, &                                                                       ! C5
      0,-1, 1,     1,-1,-1, &                                                                       ! A2
     -1, 0,-1,     1,-1,-1, &                                                                       ! A3
      1, 1, 0,     1,-1,-1, &                                                                       ! A6
      0, 1, 1,    -1, 1,-1, &                                                                       ! D1
      1, 0,-1,    -1, 1,-1, &                                                                       ! D4
     -1,-1, 0,    -1, 1,-1, &                                                                       ! D6
     ! Slip system <110>{110}
      1, 1, 0,     1,-1, 0, &
      1,-1, 0,     1, 1, 0, &
      1, 0, 1,     1, 0,-1, &
      1, 0,-1,     1, 0, 1, &
      0, 1, 1,     0, 1,-1, &
      0, 1,-1,     0, 1, 1  &
     ],pReal),shape(LATTICE_FCC_SYSTEMSLIP))                                                        !< Slip system <110>{111} directions. Sorted according to Eisenlohr & Hantcherli

 character(len=*), dimension(2), parameter, private :: LATTICE_FCC_SLIPFAMILY_NAME = &
   ['<0 1 -1>{1 1 1}', &
    '<0 1 -1>{0 1 1}']

 real(pReal), dimension(3+3,LATTICE_FCC_NTWIN), parameter, private :: &
   LATTICE_FCC_SYSTEMTWIN = reshape(real( [&
     -2, 1, 1,     1, 1, 1, &
      1,-2, 1,     1, 1, 1, &
      1, 1,-2,     1, 1, 1, &
      2,-1, 1,    -1,-1, 1, &
     -1, 2, 1,    -1,-1, 1, &
     -1,-1,-2,    -1,-1, 1, &
     -2,-1,-1,     1,-1,-1, &
      1, 2,-1,     1,-1,-1, &
      1,-1, 2,     1,-1,-1, &
      2, 1,-1,    -1, 1,-1, &
     -1,-2,-1,    -1, 1,-1, &
     -1, 1, 2,    -1, 1,-1  &
     ],pReal),shape(LATTICE_FCC_SYSTEMTWIN))                                                        !< Twin system <112>{111} directions. Sorted according to Eisenlohr & Hantcherli

 character(len=*), dimension(1), parameter, private :: LATTICE_FCC_TWINFAMILY_NAME = &
   ['<-2 1 1>{1 1 1}']


 integer(pInt), dimension(2_pInt,LATTICE_FCC_NTWIN), parameter, public :: &
   LATTICE_FCC_TWINNUCLEATIONSLIPPAIR = reshape(int( [&
     2,3, &
     1,3, &
     1,2, &
     5,6, &
     4,6, &
     4,5, &
     8,9, &
     7,9, &
     7,8, &
     11,12, &
     10,12, &
     10,11 &
     ],pInt),shape(LATTICE_FCC_TWINNUCLEATIONSLIPPAIR))

 real(pReal), dimension(3+3,LATTICE_fcc_Ncleavage), parameter, private :: &
   LATTICE_fcc_systemCleavage = reshape(real([&
    ! Cleavage direction     Plane normal
      0, 1, 0,     1, 0, 0, &
      0, 0, 1,     0, 1, 0, &
      1, 0, 0,     0, 0, 1, &
      0, 1,-1,     1, 1, 1, &
      0,-1,-1,    -1,-1, 1, &
     -1, 0,-1,     1,-1,-1, &
      0, 1, 1,    -1, 1,-1  &
     ],pReal),shape(LATTICE_FCC_SYSTEMCLEAVAGE))

!--------------------------------------------------------------------------------------------------
! body centered cubic
 integer(pInt), dimension(LATTICE_maxNslipFamily), parameter, private :: &
   LATTICE_BCC_NSLIPSYSTEM = int([ 12, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], pInt)                  !< # of slip systems per family for bcc

 integer(pInt), dimension(1), parameter, private :: &
   LATTICE_BCC_NTWINSYSTEM = int([12], pInt)                                                        !< # of twin systems per family for bcc

 integer(pInt), dimension(LATTICE_maxNcleavageFamily), parameter, private :: &
   LATTICE_bcc_NcleavageSystem = int([3, 6, 0],pInt)                                                !< # of cleavage systems per family for bcc

 integer(pInt), parameter, private  :: &
   LATTICE_BCC_NSLIP      = sum(LATTICE_BCC_NSLIPSYSTEM), &                                         !< total # of slip systems for bcc
   LATTICE_BCC_NTWIN      = sum(LATTICE_BCC_NTWINSYSTEM), &                                         !< total # of twin systems for bcc
   LATTICE_bcc_Ncleavage  = sum(lattice_bcc_NcleavageSystem)                                        !< total # of cleavage systems for bcc

 real(pReal), dimension(3+3,LATTICE_BCC_NSLIP), parameter, private :: &
   LATTICE_bcc_systemSlip = reshape(real([&
    ! Slip direction     Plane normal
    ! Slip system <111>{110}
      1,-1, 1,     0, 1, 1, &
     -1,-1, 1,     0, 1, 1, &
      1, 1, 1,     0,-1, 1, &
     -1, 1, 1,     0,-1, 1, &
     -1, 1, 1,     1, 0, 1, &
     -1,-1, 1,     1, 0, 1, &
      1, 1, 1,    -1, 0, 1, &
      1,-1, 1,    -1, 0, 1, &
     -1, 1, 1,     1, 1, 0, &
     -1, 1,-1,     1, 1, 0, &
      1, 1, 1,    -1, 1, 0, &
      1, 1,-1,    -1, 1, 0, &
    ! Slip system <111>{112}
     -1, 1, 1,     2, 1, 1, &
      1, 1, 1,    -2, 1, 1, &
      1, 1,-1,     2,-1, 1, &
      1,-1, 1,     2, 1,-1, &
      1,-1, 1,     1, 2, 1, &
      1, 1,-1,    -1, 2, 1, &
      1, 1, 1,     1,-2, 1, &
     -1, 1, 1,     1, 2,-1, &
      1, 1,-1,     1, 1, 2, &
      1,-1, 1,    -1, 1, 2, &
     -1, 1, 1,     1,-1, 2, &
      1, 1, 1,     1, 1,-2  &
     ],pReal),shape(LATTICE_BCC_SYSTEMSLIP))

 character(len=*), dimension(2), parameter, private :: LATTICE_BCC_SLIPFAMILY_NAME = &
   ['<1 -1 1>{0 1 1}', &
    '<1 -1 1>{2 1 1}']

 real(pReal), dimension(3+3,LATTICE_BCC_NTWIN), parameter, private :: &
   LATTICE_bcc_systemTwin = reshape(real([&
    ! Twin system <111>{112}
     -1, 1, 1,     2, 1, 1, &
      1, 1, 1,    -2, 1, 1, &
      1, 1,-1,     2,-1, 1, &
      1,-1, 1,     2, 1,-1, &
      1,-1, 1,     1, 2, 1, &
      1, 1,-1,    -1, 2, 1, &
      1, 1, 1,     1,-2, 1, &
     -1, 1, 1,     1, 2,-1, &
      1, 1,-1,     1, 1, 2, &
      1,-1, 1,    -1, 1, 2, &
     -1, 1, 1,     1,-1, 2, &
      1, 1, 1,     1, 1,-2  &
     ],pReal),shape(LATTICE_BCC_SYSTEMTWIN))

 character(len=*), dimension(1), parameter, private :: LATTICE_BCC_TWINFAMILY_NAME = &
   ['<1 1 1>{2 1 1}']

 real(pReal), dimension(3+3,LATTICE_bcc_Ncleavage), parameter, private :: &
   LATTICE_bcc_systemCleavage = reshape(real([&
    ! Cleavage direction     Plane normal
      0, 1, 0,     1, 0, 0, &
      0, 0, 1,     0, 1, 0, &
      1, 0, 0,     0, 0, 1, &
      1,-1, 1,     0, 1, 1, &
      1, 1, 1,     0,-1, 1, &
     -1, 1, 1,     1, 0, 1, &
      1, 1, 1,    -1, 0, 1, &
     -1, 1, 1,     1, 1, 0, &
      1, 1, 1,    -1, 1, 0  &
     ],pReal),shape(LATTICE_BCC_SYSTEMCLEAVAGE))

!--------------------------------------------------------------------------------------------------
! hexagonal
 integer(pInt), dimension(LATTICE_maxNslipFamily), parameter, private :: &
   LATTICE_HEX_NSLIPSYSTEM = int([ 3, 3, 3, 6, 12, 6, 0, 0, 0, 0, 0, 0, 0],pInt)                    !< # of slip systems per family for hex

 integer(pInt), dimension(4), parameter, private :: &
   LATTICE_HEX_NTWINSYSTEM = int([ 6, 6, 6, 6],pInt)                                                !< # of slip systems per family for hex

 integer(pInt), dimension(LATTICE_maxNcleavageFamily), parameter, private :: &
   LATTICE_hex_NcleavageSystem = int([3, 0, 0],pInt)                                                !< # of cleavage systems per family for hex

 integer(pInt), parameter, private  :: &
   LATTICE_HEX_NSLIP      = sum(LATTICE_HEX_NSLIPSystem), &                                         !< total # of slip systems for hex
   LATTICE_HEX_NTWIN      = sum(LATTICE_HEX_NTWINSYSTEM), &                                         !< total # of twin systems for hex
   LATTICE_hex_Ncleavage  = sum(lattice_hex_NcleavageSystem)                                        !< total # of cleavage systems for hex

 real(pReal), dimension(4+4,LATTICE_HEX_NSLIP), parameter, private :: &
   LATTICE_hex_systemSlip = reshape(real([&
    ! Slip direction     Plane normal
    ! Basal systems <11.0>{00.1} (independent of c/a-ratio, Bravais notation (4 coordinate base))
      2, -1, -1,  0,     0,  0,  0,  1, &
     -1,  2, -1,  0,     0,  0,  0,  1, &
     -1, -1,  2,  0,     0,  0,  0,  1, &
    ! 1st type prismatic systems <11.0>{10.0}  (independent of c/a-ratio)
      2, -1, -1,  0,     0,  1, -1,  0, &
     -1,  2, -1,  0,    -1,  0,  1,  0, &
     -1, -1,  2,  0,     1, -1,  0,  0, &
    ! 2nd type prismatic systems <10.0>{11.0} -- a slip; plane normals independent of c/a-ratio
      0,  1,  -1, 0,     2, -1, -1,  0, &
     -1,  0,  1,  0,    -1,  2, -1,  0, &
      1, -1,  0,  0,    -1, -1,  2,  0,  &
    ! 1st type 1st order pyramidal systems <11.0>{-11.1} -- plane normals depend on the c/a-ratio
      2, -1, -1,  0,     0,  1, -1,  1, &
     -1,  2, -1,  0,    -1,  0,  1,  1, &
     -1, -1,  2,  0,     1, -1,  0,  1, &
      1,  1, -2,  0,    -1,  1,  0,  1, &
     -2,  1,  1,  0,     0, -1,  1,  1, &
      1, -2,  1,  0,     1,  0, -1,  1, &
    ! pyramidal system: c+a slip <11.3>{-10.1} -- plane normals depend on the c/a-ratio
      2, -1, -1,  3,    -1,  1,  0,  1, &
      1, -2,  1,  3,    -1,  1,  0,  1, &
     -1, -1,  2,  3,     1,  0, -1,  1, &
     -2,  1,  1,  3,     1,  0, -1,  1, &
     -1,  2, -1,  3,     0, -1,  1,  1, &
      1,  1, -2,  3,     0, -1,  1,  1, &
     -2,  1,  1,  3,     1, -1,  0,  1, &
     -1,  2, -1,  3,     1, -1,  0,  1, &
      1,  1, -2,  3,    -1,  0,  1,  1, &
      2, -1, -1,  3,    -1,  0,  1,  1, &
      1, -2,  1,  3,     0,  1, -1,  1, &
     -1, -1,  2,  3,     0,  1, -1,  1, &
    ! pyramidal system: c+a slip <11.3>{-1-1.2} -- as for hexagonal ice (Castelnau et al. 1996, similar to twin system found below)
      2, -1, -1,  3,    -2,  1,  1,  2, & ! sorted according to similar twin system
     -1,  2, -1,  3,     1, -2,  1,  2, & ! <11.3>{-1-1.2} shear = 2((c/a)^2-2)/(3 c/a)
     -1, -1,  2,  3,     1,  1, -2,  2, &
     -2,  1,  1,  3,     2, -1, -1,  2, &
      1, -2,  1,  3,    -1,  2, -1,  2, &
      1,  1, -2,  3,    -1, -1,  2,  2  &
     ],pReal),shape(LATTICE_HEX_SYSTEMSLIP))                                                        !< slip systems for hex sorted by A. Alankar & P. Eisenlohr

 character(len=*), dimension(6), parameter, private :: LATTICE_HEX_SLIPFAMILY_NAME = &
   ['<1 1 . 1>{0 0 . 1}  ', &
    '<1 1 . 1>{1 0 . 0}  ', &
    '<1 0 . 0>{1 1 . 0}  ', &
    '<1 1 . 0>{-1 1 . 1} ', &
    '<1 1 . 3>{-1 0 . 1} ', &
    '<1 1 . 3>{-1 -1 . 2}']

 real(pReal), dimension(4+4,LATTICE_hex_Ntwin), parameter, private :: &
   LATTICE_hex_systemTwin =  reshape(real([&
    ! Compression or Tension =f(twinning shear=f(c/a)) for each metal ! (according to Yoo 1981)
      1, -1,  0,  1,    -1,  1,  0,  2, & ! <-10.1>{10.2} shear = (3-(c/a)^2)/(sqrt(3) c/a)
     -1,  0,  1,  1,     1,  0, -1,  2, &
      0,  1, -1,  1,     0, -1,  1,  2, &
     -1,  1,  0,  1,     1, -1,  0,  2, &
      1,  0, -1,  1,    -1,  0,  1,  2, &
      0, -1,  1,  1,     0,  1, -1,  2, &
!
      2, -1, -1,  6,    -2,  1,  1,  1, & ! <11.6>{-1-1.1} shear = 1/(c/a)
     -1,  2, -1,  6,     1, -2,  1,  1, &
     -1, -1,  2,  6,     1,  1, -2,  1, &
     -2,  1,  1,  6,     2, -1, -1,  1, &
      1, -2,  1,  6,    -1,  2, -1,  1, &
      1,  1, -2,  6,    -1, -1,  2,  1, &
!
     -1,  1,  0, -2,    -1,  1,  0,  1, & !! <10.-2>{10.1} shear = (4(c/a)^2-9)/(4 sqrt(3) c/a)
      1,  0, -1, -2,     1,  0, -1,  1, &
      0, -1,  1, -2,     0, -1,  1,  1, &
      1, -1,  0, -2,     1, -1,  0,  1, &
     -1,  0,  1, -2,    -1,  0,  1,  1, &
      0,  1, -1, -2,     0,  1, -1,  1, &
!
      2, -1, -1, -3,     2, -1, -1,  2, & ! <11.-3>{11.2} shear = 2((c/a)^2-2)/(3 c/a)
     -1,  2, -1, -3,    -1,  2, -1,  2, &
     -1, -1,  2, -3,    -1, -1,  2,  2, &
     -2,  1,  1, -3,    -2,  1,  1,  2, &
      1, -2,  1, -3,     1, -2,  1,  2, &
      1,  1, -2, -3,     1,  1, -2,  2  &
     ],pReal),shape(LATTICE_HEX_SYSTEMTWIN))                                                        !< twin systems for hex, order follows Prof. Tom Bieler's scheme

 character(len=*), dimension(4), parameter, private :: LATTICE_HEX_TWINFAMILY_NAME = &
   ['<-1 0 . 1>{1 0 . 2} ', &
    '<1 1 . 6>{-1 -1 . 1}', &
    '<1 0 . -2>{1 0 . 1} ', &
    '<1 1 . -3>{1 1 . 2} ']

 real(pReal), dimension(4+4,LATTICE_hex_Ncleavage), parameter, private :: &
   LATTICE_hex_systemCleavage = reshape(real([&
    ! Cleavage direction     Plane normal
      2,-1,-1, 0,     0, 0, 0, 1, &
      0, 0, 0, 1,     2,-1,-1, 0, &
      0, 0, 0, 1,     0, 1,-1, 0  &
     ],pReal),shape(LATTICE_HEX_SYSTEMCLEAVAGE))


!--------------------------------------------------------------------------------------------------
! body centered tetragonal
 integer(pInt), dimension(LATTICE_maxNslipFamily), parameter, private :: &
   LATTICE_bct_NslipSystem = int([2, 2, 2, 4, 2, 4, 2, 2, 4, 8, 4, 8, 8 ],pInt)                     !< # of slip systems per family for bct (Sn) Bieler J. Electr Mater 2009

 integer(pInt), parameter, private  :: &
   LATTICE_bct_Nslip      = sum(lattice_bct_NslipSystem)                                            !< total # of slip systems for bct

 real(pReal), dimension(3+3,LATTICE_bct_Nslip), parameter, private :: &
   LATTICE_bct_systemSlip = reshape(real([&
    ! Slip direction     Plane normal
    ! Slip family 1 {100)<001] (Bravais notation {hkl)<uvw] for bct c/a = 0.5456)
      0, 0, 1,      1, 0, 0, &
      0, 0, 1,      0, 1, 0, &
    ! Slip family 2 {110)<001]
      0, 0, 1,      1, 1, 0, &
      0, 0, 1,     -1, 1, 0, &
    ! slip family 3 {100)<010]
      0,  1, 0,     1, 0, 0, &
      1,  0, 0,     0, 1, 0, &
    ! Slip family 4 {110)<1-11]/2
      1,-1, 1,      1, 1, 0, &
      1,-1,-1,      1, 1, 0, &
     -1,-1,-1,     -1, 1, 0, &
     -1,-1, 1,     -1, 1, 0, &
    ! Slip family 5 {110)<1-10]
      1, -1, 0,     1, 1, 0, &
      1,  1, 0,     1,-1, 0, &
    ! Slip family 6 {100)<011]
      0, 1, 1,      1, 0, 0, &
      0,-1, 1,      1, 0, 0, &
     -1, 0, 1,      0, 1, 0, &
      1, 0, 1,      0, 1, 0, &
    ! Slip family 7 {001)<010]
      0, 1, 0,      0, 0, 1, &
      1, 0, 0,      0, 0, 1, &
    ! Slip family 8 {001)<110]
      1, 1, 0,      0, 0, 1, &
     -1, 1, 0,      0, 0, 1, &
    ! Slip family 9 {011)<01-1]
      0, 1,-1,      0, 1, 1, &
      0,-1,-1,      0,-1, 1, &
     -1, 0,-1,     -1, 0, 1, &
      1, 0,-1,      1, 0, 1, &
    ! Slip family 10 {011)<1-11]/2
      1,-1, 1,      0, 1, 1, &
      1, 1,-1,      0, 1, 1, &
      1, 1, 1,      0, 1,-1, &
     -1, 1, 1,      0, 1,-1, &
      1,-1,-1,      1, 0, 1, &
     -1,-1, 1,      1, 0, 1, &
      1, 1, 1,      1, 0,-1, &
      1,-1, 1,      1, 0,-1,  &
    ! Slip family 11 {011)<100]
      1, 0, 0,      0, 1, 1, &
      1, 0, 0,      0, 1,-1, &
      0, 1, 0,      1, 0, 1, &
      0, 1, 0,      1, 0,-1, &
    ! Slip family 12 {211)<01-1]
      0, 1,-1,      2, 1, 1, &
      0,-1,-1,      2,-1, 1, &
      1, 0,-1,      1, 2, 1, &
     -1, 0,-1,     -1, 2, 1, &
      0, 1,-1,     -2, 1, 1, &
      0,-1,-1,     -2,-1, 1, &
     -1, 0,-1,     -1,-2, 1, &
      1, 0,-1,      1,-2, 1, &
    ! Slip family 13 {211)<-111]/2
     -1, 1, 1,      2, 1, 1, &
     -1,-1, 1,      2,-1, 1, &
      1,-1, 1,      1, 2, 1, &
     -1,-1, 1,     -1, 2, 1, &
      1, 1, 1,     -2, 1, 1, &
      1,-1, 1,     -2,-1, 1, &
     -1, 1, 1,     -1,-2, 1, &
      1, 1, 1,      1,-2, 1  &
      ],pReal),[ 3_pInt + 3_pInt,LATTICE_bct_Nslip])                                                !< slip systems for bct sorted by Bieler

 character(len=*), dimension(13), parameter, private :: LATTICE_BCT_SLIPFAMILY_NAME = &
   ['{1 0 0)<0 0 1] ', &
    '{1 1 0)<0 0 1] ', &
    '{1 0 0)<0 1 0] ', &
    '{1 1 0)<1 -1 1]', &
    '{1 1 0)<1 -1 0]', &
    '{1 0 0)<0 1 1] ', &
    '{0 0 1)<0 1 0] ', &
    '{0 0 1)<1 1 0] ', &
    '{0 1 1)<0 1 -1]', &
    '{0 1 1)<1 -1 1]', &
    '{0 1 1)<1 0 0] ', &
    '{2 1 1)<0 1 -1]', &
    '{2 1 1)<-1 1 1]']


!--------------------------------------------------------------------------------------------------
! isotropic
 integer(pInt), dimension(LATTICE_maxNcleavageFamily), parameter, private :: &
   LATTICE_iso_NcleavageSystem = int([3, 0, 0],pInt)                                                !< # of cleavage systems per family for iso

 integer(pInt), parameter, private  :: &
   LATTICE_iso_Ncleavage  = sum(lattice_iso_NcleavageSystem)                                        !< total # of cleavage systems for iso

 real(pReal), dimension(3+3,LATTICE_iso_Ncleavage), parameter, private :: &
   LATTICE_iso_systemCleavage = reshape(real([&
    ! Cleavage direction     Plane normal
      0, 1, 0,     1, 0, 0, &
      0, 0, 1,     0, 1, 0, &
      1, 0, 0,     0, 0, 1  &
     ],pReal),[ 3_pInt + 3_pInt,LATTICE_iso_Ncleavage])


!--------------------------------------------------------------------------------------------------
! orthorhombic
 integer(pInt), dimension(LATTICE_maxNcleavageFamily), parameter, private :: &
   LATTICE_ort_NcleavageSystem = int([1, 1, 1],pInt)                                                !< # of cleavage systems per family for ortho

 integer(pInt), parameter, private  :: &
   LATTICE_ort_Ncleavage  = sum(lattice_ort_NcleavageSystem)                                        !< total # of cleavage systems for ortho

 real(pReal), dimension(3+3,LATTICE_ort_Ncleavage), parameter, private :: &
   LATTICE_ort_systemCleavage = reshape(real([&
    ! Cleavage direction     Plane normal
      0, 1, 0,     1, 0, 0, &
      0, 0, 1,     0, 1, 0, &
      1, 0, 0,     0, 0, 1  &
     ],pReal),[ 3_pInt + 3_pInt,LATTICE_ort_Ncleavage])

! BEGIN DEPRECATED
 integer(pInt), parameter, public :: &
   LATTICE_maxNslip        = max(LATTICE_FCC_NSLIP,LATTICE_BCC_NSLIP,LATTICE_HEX_NSLIP, &
                                 LATTICE_bct_Nslip), &                                              !< max # of slip systems over lattice structures
   LATTICE_maxNcleavage    = max(LATTICE_fcc_Ncleavage,LATTICE_bcc_Ncleavage, &
                                 LATTICE_hex_Ncleavage, &
                                 LATTICE_iso_Ncleavage,LATTICE_ort_Ncleavage)                    !< max # of cleavage systems over lattice structures
!END DEPRECATED

 real(pReal),                          dimension(:,:,:),     allocatable, public, protected :: &
   lattice_C66
 real(pReal),                              dimension(:,:,:,:,:),   allocatable, public, protected :: &
   lattice_C3333
 real(pReal),                              dimension(:),   allocatable, public, protected :: &
   lattice_mu, lattice_nu

! SHOULD NOT BE PART OF LATTICE BEGIN
 real(pReal),                              dimension(:,:,:,:), allocatable, public, protected :: &  ! with higher-order parameters (e.g. temperature-dependent)
   lattice_thermalExpansion33
 real(pReal),                              dimension(:,:,:),   allocatable, public, protected :: &
   lattice_thermalConductivity33, &
   lattice_damageDiffusion33
 real(pReal),                              dimension(:),       allocatable, public, protected :: &
   lattice_damageMobility, &
   lattice_massDensity, &
   lattice_specificHeat, &
   lattice_referenceTemperature
! SHOULD NOT BE PART OF LATTICE END

 enum, bind(c)
   enumerator :: LATTICE_undefined_ID, &
                 LATTICE_iso_ID, &
                 LATTICE_fcc_ID, &
                 LATTICE_bcc_ID, &
                 LATTICE_hex_ID, &
                 LATTICE_bct_ID, &
                 LATTICE_ort_ID
 end enum

 integer(kind(LATTICE_undefined_ID)),        dimension(:),       allocatable, public, protected :: &
   lattice_structure, trans_lattice_structure

 
 interface lattice_forestProjection ! DEPRECATED, use lattice_forestProjection_edge
   module procedure slipProjection_transverse
 end interface lattice_forestProjection
 
 interface lattice_forestProjection_edge
   module procedure slipProjection_transverse
 end interface lattice_forestProjection_edge
 
 interface lattice_forestProjection_screw
   module procedure slipProjection_direction
 end interface lattice_forestProjection_screw
 
 interface lattice_slipProjection_modeI
   module procedure slipProjection_normal
 end interface lattice_slipProjection_modeI
 
  interface lattice_slipProjection_modeII
   module procedure slipProjection_direction
 end interface lattice_slipProjection_modeII
 
  interface lattice_slipProjection_modeIII
   module procedure slipProjection_transverse
 end interface lattice_slipProjection_modeIII


 public :: &
  lattice_init, &
  lattice_qDisorientation, &
  LATTICE_fcc_ID, &
  LATTICE_bcc_ID, &
  LATTICE_bct_ID, &
  LATTICE_hex_ID, &
  lattice_SchmidMatrix_slip, &
  lattice_SchmidMatrix_twin, &
  lattice_SchmidMatrix_trans, &
  lattice_SchmidMatrix_cleavage, &
  lattice_nonSchmidMatrix, &
  lattice_interaction_SlipSlip, &
  lattice_interaction_TwinTwin, &
  lattice_interaction_TransTrans, &
  lattice_interaction_SlipTwin, &
  lattice_interaction_SlipTrans, &
  lattice_interaction_TwinSlip, &
  lattice_characteristicShear_Twin, &
  lattice_C66_twin, &
  lattice_C66_trans, &
  lattice_forestProjection, &
  lattice_forestProjection_edge, &
  lattice_forestProjection_screw, &
  lattice_slipProjection_modeI, &
  lattice_slipProjection_modeII, &
  lattice_slipProjection_modeIII, &
  lattice_slip_normal, &
  lattice_slip_direction, &
  lattice_slip_transverse


contains

!--------------------------------------------------------------------------------------------------
!> @brief Module initialization
!--------------------------------------------------------------------------------------------------
subroutine lattice_init
 use IO, only: &
   IO_error
 use config, only: &
   config_phase

 implicit none
 integer(pInt) :: Nphases
 character(len=65536) :: &
   tag  = ''
 integer(pInt) :: i,p
 real(pReal),  dimension(:), allocatable :: &
   temp, &
   CoverA                                                                                           !< c/a ratio for low symmetry type lattice


 write(6,'(/,a)') ' <<<+-  lattice init  -+>>>'

 Nphases = size(config_phase)

 allocate(lattice_structure(Nphases),source = LATTICE_undefined_ID)
 allocate(trans_lattice_structure(Nphases),source = LATTICE_undefined_ID)
 allocate(lattice_C66(6,6,Nphases),  source=0.0_pReal)
 allocate(lattice_C3333(3,3,3,3,Nphases),  source=0.0_pReal)
 
 allocate(lattice_thermalExpansion33     (3,3,3,Nphases), source=0.0_pReal)                       ! constant, linear, quadratic coefficients
 allocate(lattice_thermalConductivity33  (3,3,Nphases), source=0.0_pReal)
 allocate(lattice_damageDiffusion33      (3,3,Nphases), source=0.0_pReal)
 allocate(lattice_damageMobility         (    Nphases), source=0.0_pReal)
 allocate(lattice_massDensity            (    Nphases), source=0.0_pReal)
 allocate(lattice_specificHeat           (    Nphases), source=0.0_pReal)
 allocate(lattice_referenceTemperature   (    Nphases), source=300.0_pReal)

 allocate(lattice_mu(Nphases),       source=0.0_pReal)
 allocate(lattice_nu(Nphases),       source=0.0_pReal)

 allocate(lattice_NslipSystem(lattice_maxNslipFamily,Nphases),source=0_pInt)

 allocate(lattice_Scleavage(3,3,3,lattice_maxNslip,Nphases),source=0.0_pReal)
 allocate(lattice_NcleavageSystem(lattice_maxNcleavageFamily,Nphases),source=0_pInt)

 allocate(CoverA(Nphases),source=0.0_pReal)

 allocate(lattice_sd(3,lattice_maxNslip,Nphases),source=0.0_pReal)
 allocate(lattice_st(3,lattice_maxNslip,Nphases),source=0.0_pReal)
 allocate(lattice_sn(3,lattice_maxNslip,Nphases),source=0.0_pReal)


 do p = 1, size(config_phase)
   tag = config_phase(p)%getString('lattice_structure')
   select case(trim(tag))
     case('iso','isotropic')
       lattice_structure(p) = LATTICE_iso_ID
     case('fcc')
       lattice_structure(p) = LATTICE_fcc_ID
     case('bcc')
       lattice_structure(p) = LATTICE_bcc_ID
     case('hex','hexagonal')
       lattice_structure(p) = LATTICE_hex_ID
     case('bct')
       lattice_structure(p) = LATTICE_bct_ID
     case('ort','orthorhombic')
       lattice_structure(p) = LATTICE_ort_ID
   end select

   tag = 'undefined'
   tag = config_phase(p)%getString('trans_lattice_structure',defaultVal=tag)
   select case(trim(tag))
     case('bcc')
        trans_lattice_structure(p) = LATTICE_bcc_ID
     case('hex','hexagonal')
        trans_lattice_structure(p) = LATTICE_hex_ID
   end select

   lattice_C66(1,1,p) = config_phase(p)%getFloat('c11',defaultVal=0.0_pReal)
   lattice_C66(1,2,p) = config_phase(p)%getFloat('c12',defaultVal=0.0_pReal)
   lattice_C66(1,3,p) = config_phase(p)%getFloat('c13',defaultVal=0.0_pReal)
   lattice_C66(2,2,p) = config_phase(p)%getFloat('c22',defaultVal=0.0_pReal)
   lattice_C66(2,3,p) = config_phase(p)%getFloat('c23',defaultVal=0.0_pReal)
   lattice_C66(3,3,p) = config_phase(p)%getFloat('c33',defaultVal=0.0_pReal)
   lattice_C66(4,4,p) = config_phase(p)%getFloat('c44',defaultVal=0.0_pReal)
   lattice_C66(5,5,p) = config_phase(p)%getFloat('c55',defaultVal=0.0_pReal)
   lattice_C66(6,6,p) = config_phase(p)%getFloat('c66',defaultVal=0.0_pReal)


   CoverA(p)       = config_phase(p)%getFloat('c/a',defaultVal=0.0_pReal)

   lattice_thermalConductivity33(1,1,p) = config_phase(p)%getFloat('thermal_conductivity11',defaultVal=0.0_pReal)
   lattice_thermalConductivity33(2,2,p) = config_phase(p)%getFloat('thermal_conductivity22',defaultVal=0.0_pReal)
   lattice_thermalConductivity33(3,3,p) = config_phase(p)%getFloat('thermal_conductivity33',defaultVal=0.0_pReal)

   temp = config_phase(p)%getFloats('thermal_expansion11',defaultVal=[0.0_pReal])                   ! read up to three parameters (constant, linear, quadratic with T)
   lattice_thermalExpansion33(1,1,1:size(temp),p) = temp
   temp = config_phase(p)%getFloats('thermal_expansion22',defaultVal=[0.0_pReal])                   ! read up to three parameters (constant, linear, quadratic with T)
   lattice_thermalExpansion33(2,2,1:size(temp),p) = temp
   temp = config_phase(p)%getFloats('thermal_expansion33',defaultVal=[0.0_pReal])                   ! read up to three parameters (constant, linear, quadratic with T)
   lattice_thermalExpansion33(3,3,1:size(temp),p) = temp

   lattice_specificHeat(p) = config_phase(p)%getFloat( 'specific_heat',defaultVal=0.0_pReal)
   lattice_massDensity(p) = config_phase(p)%getFloat( 'mass_density',defaultVal=0.0_pReal)
   lattice_referenceTemperature(p) = config_phase(p)%getFloat( 'reference_temperature',defaultVal=0.0_pReal)
   lattice_DamageDiffusion33(1,1,p) = config_phase(p)%getFloat( 'damage_diffusion11',defaultVal=0.0_pReal)
   lattice_DamageDiffusion33(2,2,p) = config_phase(p)%getFloat( 'damage_diffusion22',defaultVal=0.0_pReal)
   lattice_DamageDiffusion33(3,3,p) = config_phase(p)%getFloat( 'damage_diffusion33',defaultVal=0.0_pReal)
   lattice_DamageMobility(p) = config_phase(p)%getFloat( 'damage_mobility',defaultVal=0.0_pReal)
 enddo

 do i = 1_pInt,Nphases
   if ((CoverA(i) < 1.0_pReal .or. CoverA(i) > 2.0_pReal) &
       .and. lattice_structure(i) == LATTICE_hex_ID) call IO_error(131_pInt,el=i)                        ! checking physical significance of c/a
   if ((CoverA(i) > 2.0_pReal) &
       .and. lattice_structure(i) == LATTICE_bct_ID) call IO_error(131_pInt,el=i)                        ! checking physical significance of c/a
   call lattice_initializeStructure(i, CoverA(i))
 enddo

end subroutine lattice_init


!--------------------------------------------------------------------------------------------------
!> @brief !!!!!!!DEPRECTATED!!!!!!
!--------------------------------------------------------------------------------------------------
subroutine lattice_initializeStructure(myPhase,CoverA)
 use prec, only: &
  tol_math_check
 use math, only: &
   math_mul33x33, &
   math_sym3333to66, &
   math_Voigt66to3333, &
   math_crossproduct
 use IO, only: &
   IO_error

 implicit none
 integer(pInt), intent(in) :: myPhase
 real(pReal), intent(in) :: &
   CoverA

 real(pReal), dimension(3,lattice_maxNslip) :: &
   sd,  sn
 integer(pInt) :: &
  i,  &
  myNslip, myNcleavage

 lattice_C66(1:6,1:6,myPhase) = lattice_symmetrizeC66(lattice_structure(myPhase),&
                                                      lattice_C66(1:6,1:6,myPhase))

 lattice_mu(myPhase) = 0.2_pReal *(  lattice_C66(1,1,myPhase) &
                                   - lattice_C66(1,2,myPhase) &
                                   + 3.0_pReal*lattice_C66(4,4,myPhase))                             ! (C11iso-C12iso)/2 with C11iso=(3*C11+2*C12+4*C44)/5 and C12iso=(C11+4*C12-2*C44)/5
 lattice_nu(myPhase) = (  lattice_C66(1,1,myPhase) &
                        + 4.0_pReal*lattice_C66(1,2,myPhase) &
                        - 2.0_pReal*lattice_C66(4,4,myPhase)) &
                                                             /(  4.0_pReal*lattice_C66(1,1,myPhase) &
                                                               + 6.0_pReal*lattice_C66(1,2,myPhase) &
                                                               + 2.0_pReal*lattice_C66(4,4,myPhase))! C12iso/(C11iso+C12iso) with C11iso=(3*C11+2*C12+4*C44)/5 and C12iso=(C11+4*C12-2*C44)/5
 lattice_C3333(1:3,1:3,1:3,1:3,myPhase) = math_Voigt66to3333(lattice_C66(1:6,1:6,myPhase))          ! Literature data is Voigt
 lattice_C66(1:6,1:6,myPhase) = math_sym3333to66(lattice_C3333(1:3,1:3,1:3,1:3,myPhase))            ! DAMASK uses Mandel-weighting
 do i = 1_pInt, 6_pInt
   if (abs(lattice_C66(i,i,myPhase))<tol_math_check) &
     call IO_error(135_pInt,el=i,ip=myPhase,ext_msg='matrix diagonal "el"ement of phase "ip"')
 enddo

 forall (i = 1_pInt:3_pInt) &
   lattice_thermalExpansion33 (1:3,1:3,i,myPhase) = lattice_symmetrize33(lattice_structure(myPhase),&
                                                                         lattice_thermalExpansion33   (1:3,1:3,i,myPhase))

 lattice_thermalConductivity33  (1:3,1:3,myPhase) = lattice_symmetrize33(lattice_structure(myPhase),&
                                                                         lattice_thermalConductivity33  (1:3,1:3,myPhase))
 lattice_DamageDiffusion33      (1:3,1:3,myPhase) = lattice_symmetrize33(lattice_structure(myPhase),&
                                                                         lattice_DamageDiffusion33      (1:3,1:3,myPhase))
 myNslip       = 0_pInt
 myNcleavage   = 0_pInt

 select case(lattice_structure(myPhase))
!--------------------------------------------------------------------------------------------------
! fcc
   case (LATTICE_fcc_ID)
     myNslip     = LATTICE_FCC_NSLIP
     myNcleavage = lattice_fcc_Ncleavage
     lattice_NslipSystem    (1:lattice_maxNslipFamily,myPhase)      = lattice_fcc_NslipSystem
     lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,myPhase)  = lattice_fcc_NcleavageSystem

     lattice_Scleavage(1:3,1:3,1:3,1:myNcleavage,myPhase) = &
     lattice_SchmidMatrix_cleavage(lattice_fcc_ncleavageSystem,'fcc',covera)

     do i = 1_pInt,myNslip
       sd(1:3,i) = lattice_fcc_systemSlip(1:3,i)
       sn(1:3,i) = lattice_fcc_systemSlip(4:6,i)
     enddo


!--------------------------------------------------------------------------------------------------
! bcc
   case (LATTICE_bcc_ID)
     myNslip     = LATTICE_BCC_NSLIP
     myNcleavage = lattice_bcc_Ncleavage
     lattice_NslipSystem(1:lattice_maxNslipFamily,myPhase)          = lattice_bcc_NslipSystem
     lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,myPhase)  = lattice_bcc_NcleavageSystem

     lattice_Scleavage(1:3,1:3,1:3,1:myNcleavage,myPhase) = &
     lattice_SchmidMatrix_cleavage(lattice_bcc_ncleavagesystem,'bcc',covera)

     do i = 1_pInt,myNslip
       sd(1:3,i) = lattice_bcc_systemSlip(1:3,i)
       sn(1:3,i) = lattice_bcc_systemSlip(4:6,i)
     enddo

!--------------------------------------------------------------------------------------------------
! hex (including conversion from miller-bravais (a1=a2=a3=c) to miller (a, b, c) indices)
   case (LATTICE_hex_ID)
     myNslip     = LATTICE_HEX_NSLIP
     myNcleavage = lattice_hex_Ncleavage
     lattice_NslipSystem(1:lattice_maxNslipFamily,myPhase)         = LATTICE_HEX_NSLIPSystem
     lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,myPhase) = lattice_hex_NcleavageSystem

     lattice_Scleavage(1:3,1:3,1:3,1:myNcleavage,myPhase) = &
     lattice_SchmidMatrix_cleavage(lattice_hex_ncleavagesystem,'hex',covera)

     do i = 1_pInt,myNslip                                                                          ! assign slip system vectors
       sd(1,i) =  lattice_hex_systemSlip(1,i)*1.5_pReal                                             ! direction [uvtw]->[3u/2 (u+2v)*sqrt(3)/2 w*(c/a)]
       sd(2,i) = (lattice_hex_systemSlip(1,i)+2.0_pReal*lattice_hex_systemSlip(2,i))*&
                                                                      0.5_pReal*sqrt(3.0_pReal)
       sd(3,i) =  lattice_hex_systemSlip(4,i)*CoverA
       sn(1,i) =  lattice_hex_systemSlip(5,i)                                                       ! plane (hkil)->(h (h+2k)/sqrt(3) l/(c/a))
       sn(2,i) = (lattice_hex_systemSlip(5,i)+2.0_pReal*lattice_hex_systemSlip(6,i))/sqrt(3.0_pReal)
       sn(3,i) =  lattice_hex_systemSlip(8,i)/CoverA
     enddo

!--------------------------------------------------------------------------------------------------
! bct
   case (LATTICE_bct_ID)
     myNslip = lattice_bct_Nslip
     lattice_NslipSystem(1:lattice_maxNslipFamily,myPhase)          = lattice_bct_NslipSystem

     do i = 1_pInt,myNslip                                                                          ! assign slip system vectors
       sd(1:2,i) = lattice_bct_systemSlip(1:2,i)
       sd(3,i) = lattice_bct_systemSlip(3,i)*CoverA
       sn(1:2,i) = lattice_bct_systemSlip(4:5,i)
       sn(3,i) =  lattice_bct_systemSlip(6,i)/CoverA
     enddo

!--------------------------------------------------------------------------------------------------
! orthorhombic (no crystal plasticity)
   case (LATTICE_ort_ID)
     myNcleavage   = lattice_ort_Ncleavage
     lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,myPhase)  = lattice_ort_NcleavageSystem

     lattice_Scleavage(1:3,1:3,1:3,1:myNcleavage,myPhase) = &
     lattice_SchmidMatrix_cleavage(lattice_ort_NcleavageSystem,'ort',covera)

!--------------------------------------------------------------------------------------------------
! isotropic (no crystal plasticity)
   case (LATTICE_iso_ID)
     myNcleavage   = lattice_iso_Ncleavage
     lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,myPhase)  = lattice_iso_NcleavageSystem

     lattice_Scleavage(1:3,1:3,1:3,1:myNcleavage,myPhase) = &
     lattice_SchmidMatrix_cleavage(lattice_iso_NcleavageSystem,'iso',covera)

!--------------------------------------------------------------------------------------------------
! something went wrong
   case default
     call IO_error(130_pInt,ext_msg='lattice_initializeStructure')
 end select


 do i = 1_pInt,myNslip                                                                              ! store slip system vectors and Schmid matrix for my structure
   lattice_sd(1:3,i,myPhase) = sd(1:3,i)/norm2(sd(1:3,i))                                           ! make unit vector
   lattice_sn(1:3,i,myPhase) = sn(1:3,i)/norm2(sn(1:3,i))                                           ! make unit vector
   lattice_st(1:3,i,myPhase) = math_crossproduct(lattice_sd(1:3,i,myPhase),lattice_sn(1:3,i,myPhase))
 enddo

end subroutine lattice_initializeStructure


!--------------------------------------------------------------------------------------------------
!> @brief Symmetrizes stiffness matrix according to lattice type
!> @details J. A. Rayne and B. S. Chandrasekhar Phys. Rev. 120, 1658 Erratum Phys. Rev. 122, 1962
!--------------------------------------------------------------------------------------------------
pure function lattice_symmetrizeC66(struct,C66)

 implicit none
 integer(kind(LATTICE_undefined_ID)), intent(in) :: struct
 real(pReal), dimension(6,6), intent(in) :: C66
 real(pReal), dimension(6,6) :: lattice_symmetrizeC66
 integer(pInt) :: j,k

 lattice_symmetrizeC66 = 0.0_pReal

 select case(struct)
   case (LATTICE_iso_ID)
     forall(k=1_pInt:3_pInt)
       forall(j=1_pInt:3_pInt) lattice_symmetrizeC66(k,j) = C66(1,2)
       lattice_symmetrizeC66(k,k) = C66(1,1)
       lattice_symmetrizeC66(k+3,k+3) = 0.5_pReal*(C66(1,1)-C66(1,2))
     end forall
   case (LATTICE_fcc_ID,LATTICE_bcc_ID)
     forall(k=1_pInt:3_pInt)
       forall(j=1_pInt:3_pInt) lattice_symmetrizeC66(k,j) =   C66(1,2)
       lattice_symmetrizeC66(k,k) =     C66(1,1)
       lattice_symmetrizeC66(k+3_pInt,k+3_pInt) = C66(4,4)
     end forall
   case (LATTICE_hex_ID)
     lattice_symmetrizeC66(1,1) = C66(1,1)
     lattice_symmetrizeC66(2,2) = C66(1,1)
     lattice_symmetrizeC66(3,3) = C66(3,3)
     lattice_symmetrizeC66(1,2) = C66(1,2)
     lattice_symmetrizeC66(2,1) = C66(1,2)
     lattice_symmetrizeC66(1,3) = C66(1,3)
     lattice_symmetrizeC66(3,1) = C66(1,3)
     lattice_symmetrizeC66(2,3) = C66(1,3)
     lattice_symmetrizeC66(3,2) = C66(1,3)
     lattice_symmetrizeC66(4,4) = C66(4,4)
     lattice_symmetrizeC66(5,5) = C66(4,4)
     lattice_symmetrizeC66(6,6) = 0.5_pReal*(C66(1,1)-C66(1,2))
   case (LATTICE_ort_ID)
     lattice_symmetrizeC66(1,1) = C66(1,1)
     lattice_symmetrizeC66(2,2) = C66(2,2)
     lattice_symmetrizeC66(3,3) = C66(3,3)
     lattice_symmetrizeC66(1,2) = C66(1,2)
     lattice_symmetrizeC66(2,1) = C66(1,2)
     lattice_symmetrizeC66(1,3) = C66(1,3)
     lattice_symmetrizeC66(3,1) = C66(1,3)
     lattice_symmetrizeC66(2,3) = C66(2,3)
     lattice_symmetrizeC66(3,2) = C66(2,3)
     lattice_symmetrizeC66(4,4) = C66(4,4)
     lattice_symmetrizeC66(5,5) = C66(5,5)
     lattice_symmetrizeC66(6,6) = C66(6,6)
   case (LATTICE_bct_ID)
     lattice_symmetrizeC66(1,1) = C66(1,1)
     lattice_symmetrizeC66(2,2) = C66(1,1)
     lattice_symmetrizeC66(3,3) = C66(3,3)
     lattice_symmetrizeC66(1,2) = C66(1,2)
     lattice_symmetrizeC66(2,1) = C66(1,2)
     lattice_symmetrizeC66(1,3) = C66(1,3)
     lattice_symmetrizeC66(3,1) = C66(1,3)
     lattice_symmetrizeC66(2,3) = C66(1,3)
     lattice_symmetrizeC66(3,2) = C66(1,3)
     lattice_symmetrizeC66(4,4) = C66(4,4)
     lattice_symmetrizeC66(5,5) = C66(4,4)
     lattice_symmetrizeC66(6,6) = C66(6,6)
   case default
     lattice_symmetrizeC66 = C66
  end select

 end function lattice_symmetrizeC66


!--------------------------------------------------------------------------------------------------
!> @brief Symmetrizes 2nd order tensor according to lattice type
!--------------------------------------------------------------------------------------------------
pure function lattice_symmetrize33(struct,T33)

 implicit none
 integer(kind(LATTICE_undefined_ID)), intent(in) :: struct
 real(pReal), dimension(3,3), intent(in) :: T33
 real(pReal), dimension(3,3) :: lattice_symmetrize33
 integer(pInt) :: k

 lattice_symmetrize33 = 0.0_pReal

 select case(struct)
   case (LATTICE_iso_ID,LATTICE_fcc_ID,LATTICE_bcc_ID)
     forall(k=1_pInt:3_pInt) lattice_symmetrize33(k,k) = T33(1,1)
   case (LATTICE_hex_ID)
     lattice_symmetrize33(1,1) = T33(1,1)
     lattice_symmetrize33(2,2) = T33(1,1)
     lattice_symmetrize33(3,3) = T33(3,3)
   case (LATTICE_ort_ID,lattice_bct_ID)
     lattice_symmetrize33(1,1) = T33(1,1)
     lattice_symmetrize33(2,2) = T33(2,2)
     lattice_symmetrize33(3,3) = T33(3,3)
   case default
     lattice_symmetrize33 = T33
  end select

 end function lattice_symmetrize33


!--------------------------------------------------------------------------------------------------
!> @brief figures whether unit quat falls into stereographic standard triangle
!--------------------------------------------------------------------------------------------------
logical pure function lattice_qInSST(Q, struct)
 use, intrinsic :: &
   IEEE_arithmetic
 use math, only: &
   math_qToRodrig

 implicit none
 real(pReal), dimension(4), intent(in) ::      Q                                                    ! orientation
 integer(kind(LATTICE_undefined_ID)), intent(in) :: struct                                          ! lattice structure
 real(pReal), dimension(3) ::                  Rodrig                                               ! Rodrigues vector of Q

 Rodrig = math_qToRodrig(Q)
 if (any(IEEE_is_NaN(Rodrig))) then
   lattice_qInSST = .false.
 else
   select case (struct)
     case (LATTICE_bcc_ID,LATTICE_fcc_ID)
       lattice_qInSST = Rodrig(1) > Rodrig(2) .and. &
                        Rodrig(2) > Rodrig(3) .and. &
                        Rodrig(3) > 0.0_pReal
     case (LATTICE_hex_ID)
       lattice_qInSST = Rodrig(1) > sqrt(3.0_pReal)*Rodrig(2) .and. &
                        Rodrig(2) > 0.0_pReal .and. &
                        Rodrig(3) > 0.0_pReal
     case default
       lattice_qInSST = .true.
   end select
 endif

end function lattice_qInSST


!--------------------------------------------------------------------------------------------------
!> @brief calculates the disorientation for 2 unit quaternions
!--------------------------------------------------------------------------------------------------
pure function lattice_qDisorientation(Q1, Q2, struct)
 use prec, only: &
  tol_math_check
 use math, only: &
   math_qMul, &
   math_qConj

 implicit none
 real(pReal), dimension(4) ::                  lattice_qDisorientation
 real(pReal), dimension(4), intent(in) :: &
   Q1, &                                                                                            ! 1st orientation
   Q2                                                                                               ! 2nd orientation
 integer(kind(LATTICE_undefined_ID)), optional, intent(in) :: &                                     ! if given, symmetries between the two orientation will be considered
   struct

 real(pReal), dimension(4) ::                  dQ,dQsymA,mis
 integer(pInt)    ::                           i,j,k,s,symmetry
 integer(kind(LATTICE_undefined_ID)) :: myStruct

 integer(pInt), dimension(2), parameter :: &
   NsymOperations = [24_pInt,12_pInt]

real(pReal), dimension(4,36), parameter :: &
  symOperations = reshape([&
     1.0_pReal,                  0.0_pReal,                  0.0_pReal,                  0.0_pReal, &                  ! cubic symmetry operations
     0.0_pReal,                  0.0_pReal,                  1.0_pReal/sqrt(2.0_pReal),  1.0_pReal/sqrt(2.0_pReal), &  !     2-fold symmetry
     0.0_pReal,                  1.0_pReal/sqrt(2.0_pReal),  0.0_pReal,                  1.0_pReal/sqrt(2.0_pReal), &
     0.0_pReal,                  1.0_pReal/sqrt(2.0_pReal),  1.0_pReal/sqrt(2.0_pReal),  0.0_pReal, &
     0.0_pReal,                  0.0_pReal,                  1.0_pReal/sqrt(2.0_pReal), -1.0_pReal/sqrt(2.0_pReal), &
     0.0_pReal,                 -1.0_pReal/sqrt(2.0_pReal),  0.0_pReal,                  1.0_pReal/sqrt(2.0_pReal), &
     0.0_pReal,                  1.0_pReal/sqrt(2.0_pReal), -1.0_pReal/sqrt(2.0_pReal),  0.0_pReal, &
     0.5_pReal,                  0.5_pReal,                  0.5_pReal,                  0.5_pReal, &                  !     3-fold symmetry
    -0.5_pReal,                  0.5_pReal,                  0.5_pReal,                  0.5_pReal, &
     0.5_pReal,                 -0.5_pReal,                  0.5_pReal,                  0.5_pReal, &
    -0.5_pReal,                 -0.5_pReal,                  0.5_pReal,                  0.5_pReal, &
     0.5_pReal,                  0.5_pReal,                 -0.5_pReal,                  0.5_pReal, &
    -0.5_pReal,                  0.5_pReal,                 -0.5_pReal,                  0.5_pReal, &
     0.5_pReal,                  0.5_pReal,                  0.5_pReal,                 -0.5_pReal, &
    -0.5_pReal,                  0.5_pReal,                  0.5_pReal,                 -0.5_pReal, &
     1.0_pReal/sqrt(2.0_pReal),  1.0_pReal/sqrt(2.0_pReal),  0.0_pReal,                  0.0_pReal, &                  !     4-fold symmetry
     0.0_pReal,                  1.0_pReal,                  0.0_pReal,                  0.0_pReal, &
    -1.0_pReal/sqrt(2.0_pReal),  1.0_pReal/sqrt(2.0_pReal),  0.0_pReal,                  0.0_pReal, &
     1.0_pReal/sqrt(2.0_pReal),  0.0_pReal,                  1.0_pReal/sqrt(2.0_pReal),  0.0_pReal, &
     0.0_pReal,                  0.0_pReal,                  1.0_pReal,                  0.0_pReal, &
    -1.0_pReal/sqrt(2.0_pReal),  0.0_pReal,                  1.0_pReal/sqrt(2.0_pReal),  0.0_pReal, &
     1.0_pReal/sqrt(2.0_pReal),  0.0_pReal,                  0.0_pReal,                  1.0_pReal/sqrt(2.0_pReal), &
     0.0_pReal,                  0.0_pReal,                  0.0_pReal,                  1.0_pReal, &
    -1.0_pReal/sqrt(2.0_pReal),  0.0_pReal,                  0.0_pReal,                  1.0_pReal/sqrt(2.0_pReal), &
!
     1.0_pReal,                  0.0_pReal,                  0.0_pReal,                  0.0_pReal, &                  ! hexagonal symmetry operations
     0.0_pReal,                  1.0_pReal,                  0.0_pReal,                  0.0_pReal, &                  !     2-fold symmetry
     0.0_pReal,                  0.0_pReal,                  1.0_pReal,                  0.0_pReal, &
     0.0_pReal,                  0.5_pReal,                  2.0_pReal/sqrt(3.0_pReal),  0.0_pReal, &
     0.0_pReal,                 -0.5_pReal,                  2.0_pReal/sqrt(3.0_pReal),  0.0_pReal, &
     0.0_pReal,                  2.0_pReal/sqrt(3.0_pReal),  0.5_pReal,                  0.0_pReal, &
     0.0_pReal,                 -2.0_pReal/sqrt(3.0_pReal),  0.5_pReal,                  0.0_pReal, &
     2.0_pReal/sqrt(3.0_pReal),  0.0_pReal,                  0.0_pReal,                  0.5_pReal, &                  !     6-fold symmetry
    -2.0_pReal/sqrt(3.0_pReal),  0.0_pReal,                  0.0_pReal,                  0.5_pReal, &
     0.5_pReal,                  0.0_pReal,                  0.0_pReal,                  2.0_pReal/sqrt(3.0_pReal), &
    -0.5_pReal,                  0.0_pReal,                  0.0_pReal,                  2.0_pReal/sqrt(3.0_pReal), &
     0.0_pReal,                  0.0_pReal,                  0.0_pReal,                  1.0_pReal &
     ],[4,36])                                                                                      !< Symmetry operations as quaternions 24 for cubic, 12 for hexagonal = 36

!--------------------------------------------------------------------------------------------------
! check if a structure with known symmetries is given
 if (present(struct)) then
   myStruct = struct
   select case (struct)
     case(LATTICE_fcc_ID,LATTICE_bcc_ID)
       symmetry = 1_pInt
     case(LATTICE_hex_ID)
       symmetry = 2_pInt
     case default
       symmetry = 0_pInt
   end select
 else
   symmetry = 0_pInt
   myStruct = LATTICE_undefined_ID
 endif


!--------------------------------------------------------------------------------------------------
! calculate misorientation, for cubic and hexagonal structure find symmetries
 dQ = math_qMul(math_qConj(Q1),Q2)
 lattice_qDisorientation = dQ

 select case(symmetry)

    case (1_pInt,2_pInt)
     s = sum(NsymOperations(1:symmetry-1_pInt))
      do i = 1_pInt,2_pInt
        dQ = math_qConj(dQ)                                                                         ! switch order of "from -- to"
       do j = 1_pInt,NsymOperations(symmetry)                                                       ! run through first crystal's symmetries
          dQsymA = math_qMul(symOperations(1:4,s+j),dQ)                                             ! apply sym
         do k = 1_pInt,NsymOperations(symmetry)                                                     ! run through 2nd crystal's symmetries
            mis = math_qMul(dQsymA,symOperations(1:4,s+k))                                          ! apply sym
            if (mis(1) < 0.0_pReal) &                                                               ! want positive angle
              mis = -mis
            if (mis(1)-lattice_qDisorientation(1) > -tol_math_check &
             .and. lattice_qInSST(mis,LATTICE_undefined_ID)) lattice_qDisorientation = mis          ! found better one
      enddo; enddo; enddo
   case (0_pInt)
     if (lattice_qDisorientation(1) < 0.0_pReal) lattice_qDisorientation = -lattice_qDisorientation ! keep omega within 0 to 180 deg
  end select

end function lattice_qDisorientation


!--------------------------------------------------------------------------------------------------
!> @brief Characteristic shear for twinning
!--------------------------------------------------------------------------------------------------
function lattice_characteristicShear_Twin(Ntwin,structure,CoverA) result(characteristicShear)
 use IO, only: &
   IO_error

 implicit none
 integer(pInt),    dimension(:),            intent(in) :: Ntwin                                     !< number of active twin systems per family
 character(len=3),                          intent(in) :: structure                                 !< lattice structure
 real(pReal),                               intent(in) :: cOverA                                    !< c/a ratio
 real(pReal),     dimension(sum(Ntwin))                :: characteristicShear
 integer(pInt) :: &
   a, &                                                                                             !< index of active system
   c, &                                                                                             !< index in complete system list
   mf, &                                                                                            !< index of my family
   ms                                                                                               !< index of my system in current family

 integer(pInt), dimension(LATTICE_HEX_NTWIN), parameter :: &
   HEX_SHEARTWIN = reshape(int( [&
     1, &  ! <-10.1>{10.2}
     1, &
     1, &
     1, &
     1, &
     1, &
     2, &  ! <11.6>{-1-1.1}
     2, &
     2, &
     2, &
     2, &
     2, &
     3, &  ! <10.-2>{10.1}
     3, &
     3, &
     3, &
     3, &
     3, &
     4, &  ! <11.-3>{11.2}
     4, &
     4, &
     4, &
     4, &
     4  &
     ],pInt),[LATTICE_HEX_NTWIN])                                                                   ! indicator to formulas below

 if (len_trim(structure) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='lattice_characteristicShear_Twin: '//trim(structure))

 a = 0_pInt
 myFamilies: do mf = 1_pInt,size(Ntwin,1)
   mySystems: do ms = 1_pInt,Ntwin(mf)
     a = a + 1_pInt
     select case(structure(1:3))
       case('fcc','bcc')
         characteristicShear(a) = 0.5_pReal*sqrt(2.0_pReal)
       case('hex')
         if (cOverA < 1.0_pReal .or. cOverA > 2.0_pReal) &
           call IO_error(131_pInt,ext_msg='lattice_characteristicShear_Twin')
         c = sum(LATTICE_HEX_NTWINSYSTEM(1:mf-1))+ms
         select case(HEX_SHEARTWIN(c))                                                              ! from Christian & Mahajan 1995 p.29
           case (1_pInt)                                                                            ! <-10.1>{10.2}
             characteristicShear(a) = (3.0_pReal-cOverA**2.0_pReal)/sqrt(3.0_pReal)/CoverA
           case (2_pInt)                                                                            ! <11.6>{-1-1.1}
             characteristicShear(a) = 1.0_pReal/cOverA
           case (3_pInt)                                                                            ! <10.-2>{10.1}
             characteristicShear(a) = (4.0_pReal*cOverA**2.0_pReal-9.0_pReal)/sqrt(48.0_pReal)/cOverA
           case (4_pInt)                                                                            ! <11.-3>{11.2}
             characteristicShear(a) = 2.0_pReal*(cOverA**2.0_pReal-2.0_pReal)/3.0_pReal/cOverA
         end select
       case default
         call IO_error(137_pInt,ext_msg='lattice_characteristicShear_Twin: '//trim(structure))
     end select
   enddo mySystems
 enddo myFamilies

end function lattice_characteristicShear_Twin


!--------------------------------------------------------------------------------------------------
!> @brief Rotated elasticity matrices for twinning in 66-vector notation
!--------------------------------------------------------------------------------------------------
function lattice_C66_twin(Ntwin,C66,structure,CoverA)
 use IO, only: &
   IO_error
 use math, only: &
   INRAD, &
   math_axisAngleToR, &
   math_sym3333to66, &
   math_66toSym3333, &
   math_rotate_forward3333

 implicit none
 integer(pInt),    dimension(:),            intent(in) :: Ntwin                                     !< number of active twin systems per family
 character(len=*),                          intent(in) :: structure                                 !< lattice structure
 real(pReal),      dimension(6,6),          intent(in) :: C66                                       !< unrotated parent stiffness matrix
 real(pReal),                               intent(in) :: cOverA                                    !< c/a ratio
 real(pReal),      dimension(6,6,sum(Ntwin))           :: lattice_C66_twin

 real(pReal),      dimension(3,3,sum(Ntwin))           :: coordinateSystem

 real(pReal),      dimension(3,3)                      :: R
 integer(pInt) :: i
 
 if (len_trim(structure) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='lattice_C66_twin: '//trim(structure))

 select case(structure(1:3))
   case('fcc')
     coordinateSystem = buildCoordinateSystem(Ntwin,LATTICE_FCC_NSLIPSYSTEM,LATTICE_FCC_SYSTEMTWIN,&
                                              trim(structure),0.0_pReal)
   case('bcc')
     coordinateSystem = buildCoordinateSystem(Ntwin,LATTICE_BCC_NSLIPSYSTEM,LATTICE_BCC_SYSTEMTWIN,&
                                              trim(structure),0.0_pReal)
   case('hex')
     coordinateSystem = buildCoordinateSystem(Ntwin,LATTICE_HEX_NSLIPSYSTEM,LATTICE_HEX_SYSTEMTWIN,&
                                              'hex',cOverA)
   case default
     call IO_error(137_pInt,ext_msg='lattice_C66_twin: '//trim(structure))
 end select

 do i = 1, sum(Ntwin)
   R = math_axisAngleToR(coordinateSystem(1:3,2,i), 180.0_pReal * INRAD)                            ! ToDo: Why always 180 deg?
   lattice_C66_twin(1:6,1:6,i) = math_sym3333to66(math_rotate_forward3333(math_66toSym3333(C66),R))
 enddo
end function lattice_C66_twin


!--------------------------------------------------------------------------------------------------
!> @brief Rotated elasticity matrices for transformation in 66-vector notation
!> ToDo: Completely untested and incomplete and undocumented
!--------------------------------------------------------------------------------------------------
function lattice_C66_trans(Ntrans,C_parent66,structure_target, &
                           CoverA_trans,a_bcc,a_fcc)
 use prec, only: &
  tol_math_check
 use IO, only: &
   IO_error
 use math, only: &
   INRAD, &
   MATH_I3, &
   math_axisAngleToR, &
   math_sym3333to66, &
   math_66toSym3333, &
   math_rotate_forward3333, &
   math_mul33x33, &
   math_tensorproduct33, &
  math_crossproduct

 implicit none
 integer(pInt),    dimension(:),  intent(in) :: Ntrans                                              !< number of active twin systems per family
 character(len=*),                intent(in) :: &
   structure_target                                                                                 !< lattice structure
 real(pReal),     dimension(6,6), intent(in) :: C_parent66
 real(pReal),     dimension(6,6)             :: C_bar66, C_target_unrotated66
 real(pReal),     dimension(3,3,3,3)         :: C_target_unrotated
 real(pReal),     dimension(6,6,sum(Ntrans)) :: lattice_C66_trans
 real(pReal),     dimension(3,3,sum(Ntrans)) :: Q,S
 real(pReal)   :: a_bcc, a_fcc, CoverA_trans
 integer(pInt) :: i

 if (len_trim(structure_target) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='lattice_C66_trans (target): '//trim(structure_target))

 !ToDo: add checks for CoverA_trans,a_fcc,a_bcc
   
!--------------------------------------------------------------------------------------------------
! elasticity matrix of the target phase in cube orientation
 if (structure_target(1:3) == 'hex') then
   C_bar66(1,1) = (C_parent66(1,1) + C_parent66(1,2) + 2.0_pReal*C_parent66(4,4))/2.0_pReal
   C_bar66(1,2) = (C_parent66(1,1) + 5.0_pReal*C_parent66(1,2) - 2.0_pReal*C_parent66(4,4))/6.0_pReal
   C_bar66(3,3) = (C_parent66(1,1) + 2.0_pReal*C_parent66(1,2) + 4.0_pReal*C_parent66(4,4))/3.0_pReal
   C_bar66(1,3) = (C_parent66(1,1) + 2.0_pReal*C_parent66(1,2) - 2.0_pReal*C_parent66(4,4))/3.0_pReal
   C_bar66(4,4) = (C_parent66(1,1) - C_parent66(1,2) + C_parent66(4,4))/3.0_pReal
   C_bar66(1,4) = (C_parent66(1,1) - C_parent66(1,2) - 2.0_pReal*C_parent66(4,4)) /(3.0_pReal*sqrt(2.0_pReal))

   C_target_unrotated66 = 0.0_pReal
   C_target_unrotated66(1,1) = C_bar66(1,1) - C_bar66(1,4)**2.0_pReal/C_bar66(4,4)
   C_target_unrotated66(1,2) = C_bar66(1,2) + C_bar66(1,4)**2.0_pReal/C_bar66(4,4)
   C_target_unrotated66(1,3) = C_bar66(1,3)
   C_target_unrotated66(3,3) = C_bar66(3,3)
   C_target_unrotated66(4,4) = C_bar66(4,4) -  C_bar66(1,4)**2.0_pReal/(0.5_pReal*(C_bar66(1,1) - C_bar66(1,2)))
   C_target_unrotated66 = lattice_symmetrizeC66(LATTICE_HEX_ID,C_target_unrotated66)
 elseif (structure_target(1:3)  == 'bcc') then
   C_target_unrotated66 = C_parent66
 else
   call IO_error(137_pInt,ext_msg='lattice_C66_trans (target): '//trim(structure_target))
 endif


 do i = 1_pInt, 6_pInt
   if (abs(C_target_unrotated66(i,i))<tol_math_check) &
   call IO_error(135_pInt,el=i,ext_msg='matrix diagonal "el"ement in transformation')
 enddo
 C_target_unrotated = math_66toSym3333(C_target_unrotated66)
 call buildTransformationSystem(Q,S,Ntrans,CoverA_trans,a_fcc,a_bcc)

 do i = 1, sum(Ntrans)
   lattice_C66_trans(1:6,1:6,i) = math_sym3333to66(math_rotate_forward3333(C_target_unrotated,Q(1:3,1:3,i)))
 enddo
end function lattice_C66_trans


!--------------------------------------------------------------------------------------------------
!> @brief Non-schmid projections for bcc with up to 6 coefficients
! Koester et al. 2012, Acta Materialia 60 (2012) 3894–3901, eq. (17)
! Gröger et al. 2008, Acta Materialia 56 (2008) 5412–5425, table 1
!--------------------------------------------------------------------------------------------------
function lattice_nonSchmidMatrix(Nslip,nonSchmidCoefficients,sense) result(nonSchmidMatrix)
 use IO, only: &
   IO_error
 use math, only: &
   INRAD, &
   math_tensorproduct33, &
   math_crossproduct, &
   math_mul33x3, &
   math_axisAngleToR
 implicit none
 integer(pInt), dimension(:),                   intent(in) :: Nslip                                 !< number of active slip systems per family
 real(pReal),   dimension(:),                   intent(in) :: nonSchmidCoefficients                 !< non-Schmid coefficients for projections
 integer(pInt),                                 intent(in) :: sense                                 !< sense (-1,+1)
 real(pReal),   dimension(1:3,1:3,sum(Nslip))              :: nonSchmidMatrix

 real(pReal),   dimension(1:3,1:3,sum(Nslip))              :: coordinateSystem                      !< coordinate system of slip system
 real(pReal),   dimension(:),                  allocatable :: &
   direction, normal, np
 integer(pInt) :: i

 if (abs(sense) /= 1_pInt) call IO_error(0_pInt,ext_msg='lattice_nonSchmidMatrix')

 coordinateSystem  = buildCoordinateSystem(Nslip,LATTICE_BCC_NSLIPSYSTEM,LATTICE_BCC_SYSTEMSLIP,&
                                           'bcc',0.0_pReal)
 coordinateSystem(1:3,1,1:sum(Nslip)) = coordinateSystem(1:3,1,1:sum(Nslip)) *real(sense,pReal)     ! convert unidirectional coordinate system
 nonSchmidMatrix = lattice_SchmidMatrix_slip(Nslip,'bcc',0.0_pReal)                                 ! Schmid contribution

 do i = 1_pInt,sum(Nslip)
   direction = coordinateSystem(1:3,1,i)
   normal    = coordinateSystem(1:3,2,i)
   np = math_mul33x3(math_axisAngleToR(direction,60.0_pReal*INRAD), normal)
   if (size(nonSchmidCoefficients)>0) nonSchmidMatrix(1:3,1:3,i) = nonSchmidMatrix(1:3,1:3,i) &
     + nonSchmidCoefficients(1) * math_tensorproduct33(direction, np)
   if (size(nonSchmidCoefficients)>1) nonSchmidMatrix(1:3,1:3,i) = nonSchmidMatrix(1:3,1:3,i) &
     + nonSchmidCoefficients(2) * math_tensorproduct33(math_crossproduct(normal, direction), normal)
   if (size(nonSchmidCoefficients)>2) nonSchmidMatrix(1:3,1:3,i) = nonSchmidMatrix(1:3,1:3,i) &
     + nonSchmidCoefficients(3) * math_tensorproduct33(math_crossproduct(np, direction), np)
   if (size(nonSchmidCoefficients)>3) nonSchmidMatrix(1:3,1:3,i) = nonSchmidMatrix(1:3,1:3,i) &
     + nonSchmidCoefficients(4) * math_tensorproduct33(normal, normal)
   if (size(nonSchmidCoefficients)>4) nonSchmidMatrix(1:3,1:3,i) = nonSchmidMatrix(1:3,1:3,i) &
     + nonSchmidCoefficients(5) * math_tensorproduct33(math_crossproduct(normal, direction), &
                                                       math_crossproduct(normal, direction))
   if (size(nonSchmidCoefficients)>5) nonSchmidMatrix(1:3,1:3,i) = nonSchmidMatrix(1:3,1:3,i) &
     + nonSchmidCoefficients(6) * math_tensorproduct33(direction, direction)
 enddo
end function lattice_nonSchmidMatrix


!--------------------------------------------------------------------------------------------------
!> @brief Slip-slip interaction matrix
!> details only active slip systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_interaction_SlipSlip(Nslip,interactionValues,structure) result(interactionMatrix)
 use IO, only: &
   IO_error

 implicit none
 integer(pInt),    dimension(:),                   intent(in) :: Nslip                              !< number of active slip systems per family
 real(pReal),      dimension(:),                   intent(in) :: interactionValues                  !< values for slip-slip interaction
 character(len=*),                                 intent(in) :: structure                          !< lattice structure
 real(pReal),     dimension(sum(Nslip),sum(Nslip))            :: interactionMatrix

 integer(pInt),   dimension(:),                   allocatable :: NslipMax
 integer(pInt),   dimension(:,:),                 allocatable :: interactionTypes
 
 integer(pInt), dimension(LATTICE_FCC_NSLIP,LATTICE_FCC_NSLIP), parameter :: &
   FCC_INTERACTIONSLIPSLIP = reshape(int( [&
      1, 2, 2, 4, 6, 5, 3, 5, 5, 4, 5, 6,  9,10, 9,10,11,12, &  ! ---> slip
      2, 1, 2, 6, 4, 5, 5, 4, 6, 5, 3, 5,  9,10,11,12, 9,10, &  ! |
      2, 2, 1, 5, 5, 3, 5, 6, 4, 6, 5, 4, 11,12, 9,10, 9,10, &  ! |
      4, 6, 5, 1, 2, 2, 4, 5, 6, 3, 5, 5,  9,10,10, 9,12,11, &  ! v slip
      6, 4, 5, 2, 1, 2, 5, 3, 5, 5, 4, 6,  9,10,12,11,10, 9, &
      5, 5, 3, 2, 2, 1, 6, 5, 4, 5, 6, 4, 11,12,10, 9,10, 9, &
      3, 5, 5, 4, 5, 6, 1, 2, 2, 4, 6, 5, 10, 9,10, 9,11,12, &
      5, 4, 6, 5, 3, 5, 2, 1, 2, 6, 4, 5, 10, 9,12,11, 9,10, &
      5, 6, 4, 6, 5, 4, 2, 2, 1, 5, 5, 3, 12,11,10, 9, 9,10, &
      4, 5, 6, 3, 5, 5, 4, 6, 5, 1, 2, 2, 10, 9, 9,10,12,11, &
      5, 3, 5, 5, 4, 6, 6, 4, 5, 2, 1, 2, 10, 9,11,12,10, 9, &
      6, 5, 4, 5, 6, 4, 5, 5, 3, 2, 2, 1, 12,11, 9,10,10, 9, &

      9, 9,11, 9, 9,11,10,10,12,10,10,12,  1, 7, 8, 8, 8, 8, &
     10,10,12,10,10,12, 9, 9,11, 9, 9,11,  7, 1, 8, 8, 8, 8, &
      9,11, 9,10,12,10,10,12,10, 9,11, 9,  8, 8, 1, 7, 8, 8, &
     10,12,10, 9,11, 9, 9,11, 9,10,12,10,  8, 8, 7, 1, 8, 8, &
     11, 9, 9,12,10,10,11, 9, 9,12,10,10,  8, 8, 8, 8, 1, 7, &
     12,10,10,11, 9, 9,12,10,10,11, 9, 9,  8, 8, 8, 8, 7, 1  &
     ],pInt),shape(FCC_INTERACTIONSLIPSLIP),order=[2,1])                                            !< Slip--slip interaction types for fcc
                                                                                                    !< 1: self interaction
                                                                                                    !< 2: coplanar interaction
                                                                                                    !< 3: collinear interaction
                                                                                                    !< 4: Hirth locks
                                                                                                    !< 5: glissile junctions
                                                                                                    !< 6: Lomer locks
                                                                                                    !< 7: crossing (similar to Hirth locks in <110>{111} for two {110} planes)
                                                                                                    !< 8: similar to Lomer locks in <110>{111} for two {110} planes
                                                                                                    !< 9: similar to Lomer locks in <110>{111} btw one {110} and one {111} plane
                                                                                                    !<10: similar to glissile junctions in <110>{111} btw one {110} and one {111} plane
                                                                                                    !<11: crossing btw one {110} and one {111} plane
                                                                                                    !<12: collinear btw one {110} and one {111} plane

 integer(pInt), dimension(LATTICE_BCC_NSLIP,LATTICE_BCC_NSLIP), parameter :: &
   BCC_INTERACTIONSLIPSLIP = reshape(int( [&
     1,2,6,6,5,4,4,3,4,3,5,4, 6,6,4,3,3,4,6,6,4,3,6,6, &  ! ---> slip
     2,1,6,6,4,3,5,4,5,4,4,3, 6,6,3,4,4,3,6,6,3,4,6,6, &  ! |
     6,6,1,2,4,5,3,4,4,5,3,4, 4,3,6,6,6,6,3,4,6,6,4,3, &  ! |
     6,6,2,1,3,4,4,5,3,4,4,5, 3,4,6,6,6,6,4,3,6,6,3,4, &  ! v slip
     5,4,4,3,1,2,6,6,3,4,5,4, 3,6,4,6,6,4,6,3,4,6,3,6, &
     4,3,5,4,2,1,6,6,4,5,4,3, 4,6,3,6,6,3,6,4,3,6,4,6, &
     4,5,3,4,6,6,1,2,5,4,3,4, 6,3,6,4,4,6,3,6,6,4,6,3, &
     3,4,4,5,6,6,2,1,4,3,4,5, 6,4,6,3,3,6,4,6,6,3,6,4, &
     4,5,4,3,3,4,5,4,1,2,6,6, 3,6,6,4,4,6,6,3,6,4,3,6, &
     3,4,5,4,4,5,4,3,2,1,6,6, 4,6,6,3,3,6,6,4,6,3,4,6, &
     5,4,3,4,5,4,3,4,6,6,1,2, 6,3,4,6,6,4,3,6,4,6,6,3, &
     4,3,4,5,4,3,4,5,6,6,2,1, 6,4,3,6,6,3,4,6,3,6,6,4, &
     !
     6,6,4,3,3,4,6,6,3,4,6,6, 1,5,6,6,5,6,6,3,5,6,3,6, &
     6,6,3,4,6,6,3,4,6,6,3,4, 5,1,6,6,6,5,3,6,6,5,6,3, &
     4,3,6,6,4,3,6,6,6,6,4,3, 6,6,1,5,6,3,5,6,3,6,5,6, &
     3,4,6,6,6,6,4,3,4,3,6,6, 6,6,5,1,3,6,6,5,6,3,6,5, &
     3,4,6,6,6,6,4,3,4,3,6,6, 5,6,6,3,1,6,5,6,5,3,6,6, &
     4,3,6,6,4,3,6,6,6,6,4,3, 6,5,3,6,6,1,6,5,3,5,6,6, &
     6,6,3,4,6,6,3,4,6,6,3,4, 6,3,5,6,5,6,1,6,6,6,5,3, &
     6,6,4,3,3,4,6,6,3,4,6,6, 3,6,6,5,6,5,6,1,6,6,3,5, &
     4,3,6,6,4,3,6,6,6,6,4,3, 5,6,3,6,5,3,6,6,1,6,6,5, &
     3,4,6,6,6,6,4,3,4,3,6,6, 6,5,6,3,3,5,6,6,6,1,5,6, &
     6,6,4,3,3,4,6,6,3,4,6,6, 3,6,5,6,6,6,5,3,6,5,1,6, &
     6,6,3,4,6,6,3,4,6,6,3,4, 6,3,6,5,6,6,3,5,5,6,6,1  &
     ],pInt),shape(BCC_INTERACTIONSLIPSLIP),order=[2,1])                                            !< Slip--slip interaction types for bcc from Queyreau et al. Int J Plast 25 (2009) 361–377
                                                                                                    !< 1: self interaction
                                                                                                    !< 2: coplanar interaction
                                                                                                    !< 3: collinear interaction
                                                                                                    !< 4: mixed-asymmetrical junction
                                                                                                    !< 5: mixed-symmetrical junction
                                                                                                    !< 6: edge junction

 integer(pInt), dimension(LATTICE_HEX_NSLIP,LATTICE_HEX_NSLIP), parameter :: &
   HEX_INTERACTIONSLIPSLIP = reshape(int( [&
      1, 2, 2,   3, 3, 3,   7, 7, 7,  13,13,13,13,13,13,  21,21,21,21,21,21,21,21,21,21,21,21,  31,31,31,31,31,31, &  ! ---> slip
      2, 1, 2,   3, 3, 3,   7, 7, 7,  13,13,13,13,13,13,  21,21,21,21,21,21,21,21,21,21,21,21,  31,31,31,31,31,31, &  ! |
      2, 2, 1,   3, 3, 3,   7, 7, 7,  13,13,13,13,13,13,  21,21,21,21,21,21,21,21,21,21,21,21,  31,31,31,31,31,31, &  ! |
    !                                                                                                                   v slip
      6, 6, 6,   4, 5, 5,   8, 8, 8,  14,14,14,14,14,14,  22,22,22,22,22,22,22,22,22,22,22,22,  32,32,32,32,32,32, &
      6, 6, 6,   5, 4, 5,   8, 8, 8,  14,14,14,14,14,14,  22,22,22,22,22,22,22,22,22,22,22,22,  32,32,32,32,32,32, &
      6, 6, 6,   5, 5, 4,   8, 8, 8,  14,14,14,14,14,14,  22,22,22,22,22,22,22,22,22,22,22,22,  32,32,32,32,32,32, &
    !
     12,12,12,  11,11,11,   9,10,10,  15,15,15,15,15,15,  23,23,23,23,23,23,23,23,23,23,23,23,  33,33,33,33,33,33, &
     12,12,12,  11,11,11,  10, 9,10,  15,15,15,15,15,15,  23,23,23,23,23,23,23,23,23,23,23,23,  33,33,33,33,33,33, &
     12,12,12,  11,11,11,  10,10, 9,  15,15,15,15,15,15,  23,23,23,23,23,23,23,23,23,23,23,23,  33,33,33,33,33,33, &
    !
     20,20,20,  19,19,19,  18,18,18,  16,17,17,17,17,17,  24,24,24,24,24,24,24,24,24,24,24,24,  34,34,34,34,34,34, &
     20,20,20,  19,19,19,  18,18,18,  17,16,17,17,17,17,  24,24,24,24,24,24,24,24,24,24,24,24,  34,34,34,34,34,34, &
     20,20,20,  19,19,19,  18,18,18,  17,17,16,17,17,17,  24,24,24,24,24,24,24,24,24,24,24,24,  34,34,34,34,34,34, &
     20,20,20,  19,19,19,  18,18,18,  17,17,17,16,17,17,  24,24,24,24,24,24,24,24,24,24,24,24,  34,34,34,34,34,34, &
     20,20,20,  19,19,19,  18,18,18,  17,17,17,17,16,17,  24,24,24,24,24,24,24,24,24,24,24,24,  34,34,34,34,34,34, &
     20,20,20,  19,19,19,  18,18,18,  17,17,17,17,17,16,  24,24,24,24,24,24,24,24,24,24,24,24,  34,34,34,34,34,34, &
    !
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  25,26,26,26,26,26,26,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,25,26,26,26,26,26,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,25,26,26,26,26,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,25,26,26,26,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,25,26,26,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,25,26,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,26,25,26,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,26,26,25,26,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,26,26,26,25,26,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,26,26,26,26,25,26,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,26,26,26,26,26,25,26,  35,35,35,35,35,35, &
     30,30,30,  29,29,29,  28,28,28,  27,27,27,27,27,27,  26,26,26,26,26,26,26,26,26,26,26,25,  35,35,35,35,35,35, &
    !
     42,42,42,  41,41,41,  40,40,40,  39,39,39,39,39,39,  38,38,38,38,38,38,38,38,38,38,38,38,  36,37,37,37,37,37, &
     42,42,42,  41,41,41,  40,40,40,  39,39,39,39,39,39,  38,38,38,38,38,38,38,38,38,38,38,38,  37,36,37,37,37,37, &
     42,42,42,  41,41,41,  40,40,40,  39,39,39,39,39,39,  38,38,38,38,38,38,38,38,38,38,38,38,  37,37,36,37,37,37, &
     42,42,42,  41,41,41,  40,40,40,  39,39,39,39,39,39,  38,38,38,38,38,38,38,38,38,38,38,38,  37,37,37,36,37,37, &
     42,42,42,  41,41,41,  40,40,40,  39,39,39,39,39,39,  38,38,38,38,38,38,38,38,38,38,38,38,  37,37,37,37,36,37, &
     42,42,42,  41,41,41,  40,40,40,  39,39,39,39,39,39,  38,38,38,38,38,38,38,38,38,38,38,38,  37,37,37,37,37,36  &
    ],pInt),shape(HEX_INTERACTIONSLIPSLIP),order=[2,1])                                             !< Slip--slip interaction types for hex (onion peel naming scheme)
    
 integer(pInt), dimension(LATTICE_BCT_NSLIP,LATTICE_BCT_NSLIP), parameter :: &
   BCT_INTERACTIONSLIPSLIP = reshape(int( [&
     1,  2,   3,  3,   7,  7,  13, 13, 13, 13,  21, 21,  31, 31, 31, 31,  43, 43,  57, 57,  73, 73, 73, 73,  91, 91, 91, 91, 91, 91, 91, 91,  111, 111, 111, 111, 133,133,133,133,133,133,133,133, 157,157,157,157,157,157,157,157, &
     2,  1,   3,  3,   7,  7,  13, 13, 13, 13,  21, 21,  31, 31, 31, 31,  43, 43,  57, 57,  73, 73, 73, 73,  91, 91, 91, 91, 91, 91, 91, 91,  111, 111, 111, 111, 133,133,133,133,133,133,133,133, 157,157,157,157,157,157,157,157, &
    !
     6,  6,   4,  5,   8,  8,  14, 14, 14, 14,  22, 22,  32, 32, 32, 32,  44, 44,  58, 58,  74, 74, 74, 74,  92, 92, 92, 92, 92, 92, 92, 92,  112, 112, 112, 112, 134,134,134,134,134,134,134,134, 158,158,158,158,158,158,158,158, &
     6,  6,   5,  4,   8,  8,  14, 14, 14, 14,  22, 22,  32, 32, 32, 32,  44, 44,  58, 58,  74, 74, 74, 74,  92, 92, 92, 92, 92, 92, 92, 92,  112, 112, 112, 112, 134,134,134,134,134,134,134,134, 158,158,158,158,158,158,158,158, &
    !
    12, 12,  11, 11,   9, 10,  15, 15, 15, 15,  23, 23,  33, 33, 33, 33,  45, 45,  59, 59,  75, 75, 75, 75,  93, 93, 93, 93, 93, 93, 93, 93,  113, 113, 113, 113, 135,135,135,135,135,135,135,135, 159,159,159,159,159,159,159,159, &
    12, 12,  11, 11,  10,  9,  15, 15, 15, 15,  23, 23,  33, 33, 33, 33,  45, 45,  59, 59,  75, 75, 75, 75,  93, 93, 93, 93, 93, 93, 93, 93,  113, 113, 113, 113, 135,135,135,135,135,135,135,135, 159,159,159,159,159,159,159,159, &
    !
    20, 20,  19, 19,  18, 18,  16, 17, 17, 17,  24, 24,  34, 34, 34, 34,  46, 46,  60, 60,  76, 76, 76, 76,  94, 94, 94, 94, 94, 94, 94, 94,  114, 114, 114, 114, 136,136,136,136,136,136,136,136, 160,160,160,160,160,160,160,160, &
    20, 20,  19, 19,  18, 18,  17, 16, 17, 17,  24, 24,  34, 34, 34, 34,  46, 46,  60, 60,  76, 76, 76, 76,  94, 94, 94, 94, 94, 94, 94, 94,  114, 114, 114, 114, 136,136,136,136,136,136,136,136, 160,160,160,160,160,160,160,160, &
    20, 20,  19, 19,  18, 18,  17, 17, 16, 17,  24, 24,  34, 34, 34, 34,  46, 46,  60, 60,  76, 76, 76, 76,  94, 94, 94, 94, 94, 94, 94, 94,  114, 114, 114, 114, 136,136,136,136,136,136,136,136, 160,160,160,160,160,160,160,160, &
    20, 20,  19, 19,  18, 18,  17, 17, 17, 16,  24, 24,  34, 34, 34, 34,  46, 46,  60, 60,  76, 76, 76, 76,  94, 94, 94, 94, 94, 94, 94, 94,  114, 114, 114, 114, 136,136,136,136,136,136,136,136, 160,160,160,160,160,160,160,160, &
    !
    30, 30,  29, 29,  28, 28,  27, 27, 27, 27,  25, 26,  35, 35, 35, 35,  47, 47,  61, 61,  77, 77, 77, 77,  95, 95, 95, 95, 95, 95, 95, 95,  115, 115, 115, 115, 137,137,137,137,137,137,137,137, 161,161,161,161,161,161,161,161, &
    30, 30,  29, 29,  28, 28,  27, 27, 27, 27,  26, 25,  35, 35, 35, 35,  47, 47,  61, 61,  77, 77, 77, 77,  95, 95, 95, 95, 95, 95, 95, 95,  115, 115, 115, 115, 137,137,137,137,137,137,137,137, 161,161,161,161,161,161,161,161, &
    !
    42, 42,  41, 41,  40, 40,  39, 39, 39, 39,  38, 38,  36, 37, 37, 37,  48, 48,  62, 62,  78, 78, 78, 78,  96, 96, 96, 96, 96, 96, 96, 96,  116, 116, 116, 116, 138,138,138,138,138,138,138,138, 162,162,162,162,162,162,162,162, &
    42, 42,  41, 41,  40, 40,  39, 39, 39, 39,  38, 38,  37, 36, 37, 37,  48, 48,  62, 62,  78, 78, 78, 78,  96, 96, 96, 96, 96, 96, 96, 96,  116, 116, 116, 116, 138,138,138,138,138,138,138,138, 162,162,162,162,162,162,162,162, &
    42, 42,  41, 41,  40, 40,  39, 39, 39, 39,  38, 38,  37, 37, 36, 37,  48, 48,  62, 62,  78, 78, 78, 78,  96, 96, 96, 96, 96, 96, 96, 96,  116, 116, 116, 116, 138,138,138,138,138,138,138,138, 162,162,162,162,162,162,162,162, &
    42, 42,  41, 41,  40, 40,  39, 39, 39, 39,  38, 38,  37, 37, 37, 36,  48, 48,  62, 62,  78, 78, 78, 78,  96, 96, 96, 96, 96, 96, 96, 96,  116, 116, 116, 116, 138,138,138,138,138,138,138,138, 162,162,162,162,162,162,162,162, &
    !
    56, 56,  55, 55,  54, 54,  53, 53, 53, 53,  52, 52,  51, 51, 51, 51,  49, 50,  63, 63,  79, 79, 79, 79,  97, 97, 97, 97, 97, 97, 97, 97,  117, 117, 117, 117, 139,139,139,139,139,139,139,139, 163,163,163,163,163,163,163,163, &
    56, 56,  55, 55,  54, 54,  53, 53, 53, 53,  52, 52,  51, 51, 51, 51,  50, 49,  63, 63,  79, 79, 79, 79,  97, 97, 97, 97, 97, 97, 97, 97,  117, 117, 117, 117, 139,139,139,139,139,139,139,139, 163,163,163,163,163,163,163,163, &
    !
    72, 72,  71, 71,  70, 70,  69, 69, 69, 69,  68, 68,  67, 67, 67, 67,  66, 66,  64, 65,  80, 80, 80, 80,  98, 98, 98, 98, 98, 98, 98, 98,  118, 118, 118, 118, 140,140,140,140,140,140,140,140, 164,164,164,164,164,164,164,164, &
    72, 72,  71, 71,  70, 70,  69, 69, 69, 69,  68, 68,  67, 67, 67, 67,  66, 66,  65, 64,  80, 80, 80, 80,  98, 98, 98, 98, 98, 98, 98, 98,  118, 118, 118, 118, 140,140,140,140,140,140,140,140, 164,164,164,164,164,164,164,164, &
    !
    90, 90,  89, 89,  88, 88,  87, 87, 87, 87,  86, 86,  85, 85, 85, 85,  84, 84,  83, 83,  81, 82, 82, 82,  99, 99, 99, 99, 99, 99, 99, 99,  119, 119, 119, 119, 141,141,141,141,141,141,141,141, 165,165,165,165,165,165,165,165, &
    90, 90,  89, 89,  88, 88,  87, 87, 87, 87,  86, 86,  85, 85, 85, 85,  84, 84,  83, 83,  82, 81, 82, 82,  99, 99, 99, 99, 99, 99, 99, 99,  119, 119, 119, 119, 141,141,141,141,141,141,141,141, 165,165,165,165,165,165,165,165, &
    90, 90,  89, 89,  88, 88,  87, 87, 87, 87,  86, 86,  85, 85, 85, 85,  84, 84,  83, 83,  82, 82, 81, 82,  99, 99, 99, 99, 99, 99, 99, 99,  119, 119, 119, 119, 141,141,141,141,141,141,141,141, 165,165,165,165,165,165,165,165, &
    90, 90,  89, 89,  88, 88,  87, 87, 87, 87,  86, 86,  85, 85, 85, 85,  84, 84,  83, 83,  82, 82, 82, 81,  99, 99, 99, 99, 99, 99, 99, 99,  119, 119, 119, 119, 141,141,141,141,141,141,141,141, 165,165,165,165,165,165,165,165, &
    !
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 100,101,101,101,101,101,101,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,100,101,101,101,101,101,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,100,101,101,101,101,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,100,101,101,101,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,101,100,101,101,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,101,101,100,101,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,101,101,101,100,101,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
   110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,101,101,101,101,100,  120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
    !
   132,132, 131,131, 130,130, 129,129,129,129, 128,128, 127,127,127,127, 126,126, 125,125, 124,124,124,124, 123,123,123,123,123,123,123,123,  121, 122, 122, 122, 143,143,143,143,143,143,143,143, 167,167,167,167,167,167,167,167, &
   132,132, 131,131, 130,130, 129,129,129,129, 128,128, 127,127,127,127, 126,126, 125,125, 124,124,124,124, 123,123,123,123,123,123,123,123,  121, 121, 122, 122, 143,143,143,143,143,143,143,143, 167,167,167,167,167,167,167,167, &
   132,132, 131,131, 130,130, 129,129,129,129, 128,128, 127,127,127,127, 126,126, 125,125, 124,124,124,124, 123,123,123,123,123,123,123,123,  121, 122, 121, 122, 143,143,143,143,143,143,143,143, 167,167,167,167,167,167,167,167, &
   132,132, 131,131, 130,130, 129,129,129,129, 128,128, 127,127,127,127, 126,126, 125,125, 124,124,124,124, 123,123,123,123,123,123,123,123,  121, 122, 122, 121, 143,143,143,143,143,143,143,143, 167,167,167,167,167,167,167,167, &
    !
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 144,145,145,145,145,145,145,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,144,145,145,145,145,145,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,145,144,145,145,145,145,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,145,145,144,145,145,145,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,145,145,145,144,145,145,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,145,145,145,145,144,145,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,145,145,145,145,145,144,145, 168,168,168,168,168,168,168,168, &
   156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147,  146, 146, 146, 146, 145,145,145,145,145,145,145,144, 168,168,168,168,168,168,168,168, &
    !
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 169,170,170,170,170,170,170,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 170,169,170,170,170,170,170,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 170,170,169,170,170,170,170,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 170,170,170,169,170,170,170,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 170,170,170,170,169,170,170,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 169,170,170,170,170,169,170,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 169,170,170,170,170,170,169,170, &
   182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173,  172, 172, 172, 172, 171,171,171,171,171,171,171,171, 169,170,170,170,170,170,170,169  &
 ],pInt),shape(BCT_INTERACTIONSLIPSLIP),order=[2,1])


 if (len_trim(structure) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='lattice_interaction_SlipSlip: '//trim(structure))

 select case(structure(1:3))
   case('fcc')
     interactionTypes = FCC_INTERACTIONSLIPSLIP
     NslipMax         = LATTICE_FCC_NSLIPSYSTEM
   case('bcc')
     interactionTypes = BCC_INTERACTIONSLIPSLIP
     NslipMax         = LATTICE_BCC_NSLIPSYSTEM
   case('hex')
     interactionTypes = HEX_INTERACTIONSLIPSLIP
     NslipMax         = LATTICE_HEX_NSLIPSYSTEM
   case('bct')
     interactionTypes = BCT_INTERACTIONSLIPSLIP
     NslipMax         = LATTICE_BCT_NSLIPSYSTEM
   case default
     call IO_error(137_pInt,ext_msg='lattice_interaction_SlipSlip: '//trim(structure))
 end select

 interactionMatrix = buildInteraction(Nslip,Nslip,NslipMax,NslipMax,interactionValues,interactionTypes)

end function lattice_interaction_SlipSlip


!--------------------------------------------------------------------------------------------------
!> @brief Twin-twin interaction matrix
!> details only active twin systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_interaction_TwinTwin(Ntwin,interactionValues,structure) result(interactionMatrix)
 use IO, only: &
   IO_error

 implicit none
 integer(pInt),    dimension(:),                   intent(in) :: Ntwin                              !< number of active twin systems per family
 real(pReal),      dimension(:),                   intent(in) :: interactionValues                  !< values for twin-twin interaction
 character(len=*),                                 intent(in) :: structure                          !< lattice structure
 real(pReal),     dimension(sum(Ntwin),sum(Ntwin))            :: interactionMatrix

 integer(pInt),   dimension(:),                   allocatable :: NtwinMax
 integer(pInt),   dimension(:,:),                 allocatable :: interactionTypes

 integer(pInt), dimension(LATTICE_FCC_NTWIN,LATTICE_FCC_NTWIN), parameter :: &
   FCC_INTERACTIONTWINTWIN = reshape(int( [&
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! ---> twin
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! |
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! |
     2,2,2,1,1,1,2,2,2,2,2,2, &  ! v twin
     2,2,2,1,1,1,2,2,2,2,2,2, &
     2,2,2,1,1,1,2,2,2,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,2,2,2,1,1,1, &
     2,2,2,2,2,2,2,2,2,1,1,1, &
     2,2,2,2,2,2,2,2,2,1,1,1  &
     ],pInt),shape(FCC_INTERACTIONTWINTWIN),order=[2,1])                                            !< Twin-twin interaction types for fcc

 integer(pInt), dimension(LATTICE_BCC_NTWIN,LATTICE_BCC_NTWIN), parameter :: &
   BCC_INTERACTIONTWINTWIN = reshape(int( [&
     1,3,3,3,3,3,3,2,3,3,2,3, &  ! ---> twin
     3,1,3,3,3,3,2,3,3,3,3,2, &  ! |
     3,3,1,3,3,2,3,3,2,3,3,3, &  ! |
     3,3,3,1,2,3,3,3,3,2,3,3, &  ! v twin
     3,3,3,2,1,3,3,3,3,2,3,3, &
     3,3,2,3,3,1,3,3,2,3,3,3, &
     3,2,3,3,3,3,1,3,3,3,3,2, &
     2,3,3,3,3,3,3,1,3,3,2,3, &
     3,3,2,3,3,2,3,3,1,3,3,3, &
     3,3,3,2,2,3,3,3,3,1,3,3, &
     2,3,3,3,3,3,3,2,3,3,1,3, &
     3,2,3,3,3,3,2,3,3,3,3,1  &
     ],pInt),shape(BCC_INTERACTIONTWINTWIN),order=[2,1])                                            !< Twin-twin interaction types for bcc
                                                                                                    !< 1: self interaction
                                                                                                    !< 2: collinear interaction
                                                                                                    !< 3: other interaction
 integer(pInt), dimension(LATTICE_HEX_NTWIN,LATTICE_HEX_NTWIN), parameter :: &
   HEX_INTERACTIONTWINTWIN = reshape(int( [&
      1, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, &  ! ---> twin
      2, 1, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, &  ! |
      2, 2, 1, 2, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, &  ! |
      2, 2, 2, 1, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, &  ! v twin
      2, 2, 2, 2, 1, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, &
      2, 2, 2, 2, 2, 1,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, &
    !
      6, 6, 6, 6, 6, 6,   4, 5, 5, 5, 5, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
      6, 6, 6, 6, 6, 6,   5, 4, 5, 5, 5, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
      6, 6, 6, 6, 6, 6,   5, 5, 4, 5, 5, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
      6, 6, 6, 6, 6, 6,   5, 5, 5, 4, 5, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
      6, 6, 6, 6, 6, 6,   5, 5, 5, 5, 4, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
      6, 6, 6, 6, 6, 6,   5, 5, 5, 5, 5, 4,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
    !
     12,12,12,12,12,12,  11,11,11,11,11,11,   9,10,10,10,10,10,  15,15,15,15,15,15, &
     12,12,12,12,12,12,  11,11,11,11,11,11,  10, 9,10,10,10,10,  15,15,15,15,15,15, &
     12,12,12,12,12,12,  11,11,11,11,11,11,  10,10, 9,10,10,10,  15,15,15,15,15,15, &
     12,12,12,12,12,12,  11,11,11,11,11,11,  10,10,10, 9,10,10,  15,15,15,15,15,15, &
     12,12,12,12,12,12,  11,11,11,11,11,11,  10,10,10,10, 9,10,  15,15,15,15,15,15, &
     12,12,12,12,12,12,  11,11,11,11,11,11,  10,10,10,10,10, 9,  15,15,15,15,15,15, &
    !
     20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  16,17,17,17,17,17, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,16,17,17,17,17, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,17,16,17,17,17, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,17,17,16,17,17, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,17,17,17,16,17, &
     20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,17,17,17,17,16  &
     ],pInt),shape(HEX_INTERACTIONTWINTWIN),order=[2,1])                                            !< Twin-twin interaction types for hex

 if (len_trim(structure) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='lattice_interaction_TwinTwin: '//trim(structure))

 select case(structure(1:3))
   case('fcc')
     interactionTypes = FCC_INTERACTIONTWINTWIN
     NtwinMax         = LATTICE_FCC_NTWINSYSTEM
   case('bcc')
     interactionTypes = BCC_INTERACTIONTWINTWIN
     NtwinMax         = LATTICE_BCC_NTWINSYSTEM
   case('hex')
     interactionTypes = HEX_INTERACTIONTWINTWIN
     NtwinMax         = LATTICE_HEX_NTWINSYSTEM
   case default
     call IO_error(137_pInt,ext_msg='lattice_interaction_TwinTwin: '//trim(structure))
 end select

 interactionMatrix = buildInteraction(Ntwin,Ntwin,NtwinMax,NtwinMax,interactionValues,interactionTypes)

end function lattice_interaction_TwinTwin


!--------------------------------------------------------------------------------------------------
!> @brief Trans-trans interaction matrix
!> details only active trans systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_interaction_TransTrans(Ntrans,interactionValues,structure) result(interactionMatrix)
 use IO, only: &
   IO_error

 implicit none
 integer(pInt),    dimension(:),                   intent(in) :: Ntrans                             !< number of active trans systems per family
 real(pReal),      dimension(:),                   intent(in) :: interactionValues                  !< values for trans-trans interaction
 character(len=*),                                 intent(in) :: structure                          !< lattice structure (parent crystal)
 real(pReal),     dimension(sum(Ntrans),sum(Ntrans))          :: interactionMatrix

 integer(pInt),   dimension(:),                   allocatable :: NtransMax
 integer(pInt),   dimension(:,:),                 allocatable :: interactionTypes

 integer(pInt), dimension(LATTICE_FCC_NTRANS,LATTICE_FCC_NTRANS), parameter :: &
   FCC_INTERACTIONTRANSTRANS = reshape(int( [&
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! ---> trans
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! |
     1,1,1,2,2,2,2,2,2,2,2,2, &  ! |
     2,2,2,1,1,1,2,2,2,2,2,2, &  ! v trans
     2,2,2,1,1,1,2,2,2,2,2,2, &
     2,2,2,1,1,1,2,2,2,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,1,1,1,2,2,2, &
     2,2,2,2,2,2,2,2,2,1,1,1, &
     2,2,2,2,2,2,2,2,2,1,1,1, &
     2,2,2,2,2,2,2,2,2,1,1,1  &
     ],pInt),shape(FCC_INTERACTIONTRANSTRANS),order=[2,1])                                          !< Trans-trans interaction types for fcc

 if (len_trim(structure) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='lattice_interaction_TransTrans: '//trim(structure))

 if(structure(1:3) == 'fcc') then
   interactionTypes = FCC_INTERACTIONTRANSTRANS
   NtransMax        = LATTICE_FCC_NTRANSSYSTEM
 else
   call IO_error(137_pInt,ext_msg='lattice_interaction_TransTrans: '//trim(structure))
 end if

 interactionMatrix = buildInteraction(Ntrans,Ntrans,NtransMax,NtransMax,interactionValues,interactionTypes)

end function lattice_interaction_TransTrans


!--------------------------------------------------------------------------------------------------
!> @brief Slip-twin interaction matrix
!> details only active slip and twin systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_interaction_SlipTwin(Nslip,Ntwin,interactionValues,structure) result(interactionMatrix)
 use IO, only: &
   IO_error

 implicit none
 integer(pInt),    dimension(:),                   intent(in) :: Nslip, &                           !< number of active slip systems per family
                                                                 Ntwin                              !< number of active twin systems per family
 real(pReal),      dimension(:),                   intent(in) :: interactionValues                  !< values for slip-twin interaction
 character(len=*),                                 intent(in) :: structure                          !< lattice structure
 real(pReal),     dimension(sum(Nslip),sum(Ntwin))            :: interactionMatrix

 integer(pInt),   dimension(:),                   allocatable :: NslipMax, &
                                                                 NtwinMax
 integer(pInt),   dimension(:,:),                 allocatable :: interactionTypes

 integer(pInt), dimension(LATTICE_FCC_NSLIP,LATTICE_FCC_NTWIN), parameter :: &
   FCC_INTERACTIONSLIPTWIN = reshape(int( [&
     1,1,1,3,3,3,2,2,2,3,3,3, & ! ---> twin
     1,1,1,3,3,3,3,3,3,2,2,2, & ! |
     1,1,1,2,2,2,3,3,3,3,3,3, & ! |
     3,3,3,1,1,1,3,3,3,2,2,2, & ! v slip
     3,3,3,1,1,1,2,2,2,3,3,3, &
     2,2,2,1,1,1,3,3,3,3,3,3, &
     2,2,2,3,3,3,1,1,1,3,3,3, &
     3,3,3,2,2,2,1,1,1,3,3,3, &
     3,3,3,3,3,3,1,1,1,2,2,2, &
     3,3,3,2,2,2,3,3,3,1,1,1, &
     2,2,2,3,3,3,3,3,3,1,1,1, &
     3,3,3,3,3,3,2,2,2,1,1,1, &

     4,4,4,4,4,4,4,4,4,4,4,4, &
     4,4,4,4,4,4,4,4,4,4,4,4, &
     4,4,4,4,4,4,4,4,4,4,4,4, &
     4,4,4,4,4,4,4,4,4,4,4,4, &
     4,4,4,4,4,4,4,4,4,4,4,4, &
     4,4,4,4,4,4,4,4,4,4,4,4  &
     ],pInt),shape(FCC_INTERACTIONSLIPTWIN),order=[2,1])                                            !< Slip-twin interaction types for fcc
                                                                                                    !< 1: coplanar interaction
                                                                                                    !< 2: screw trace between slip system and twin habit plane (easy cross slip)
                                                                                                    !< 3: other interaction
 integer(pInt), dimension(LATTICE_BCC_NSLIP,LATTICE_BCC_NTWIN), parameter :: &
  BCC_INTERACTIONSLIPTWIN = reshape(int( [&
     3,3,3,2,2,3,3,3,3,2,3,3, &  ! ---> twin
     3,3,2,3,3,2,3,3,2,3,3,3, &  ! |
     3,2,3,3,3,3,2,3,3,3,3,2, &  ! |
     2,3,3,3,3,3,3,2,3,3,2,3, &  ! v slip
     2,3,3,3,3,3,3,2,3,3,2,3, &
     3,3,2,3,3,2,3,3,2,3,3,3, &
     3,2,3,3,3,3,2,3,3,3,3,2, &
     3,3,3,2,2,3,3,3,3,2,3,3, &
     2,3,3,3,3,3,3,2,3,3,2,3, &
     3,3,3,2,2,3,3,3,3,2,3,3, &
     3,2,3,3,3,3,2,3,3,3,3,2, &
     3,3,2,3,3,2,3,3,2,3,3,3, &
     !
     1,3,3,3,3,3,3,2,3,3,2,3, &
     3,1,3,3,3,3,2,3,3,3,3,2, &
     3,3,1,3,3,2,3,3,2,3,3,3, &
     3,3,3,1,2,3,3,3,3,2,3,3, &
     3,3,3,2,1,3,3,3,3,2,3,3, &
     3,3,2,3,3,1,3,3,2,3,3,3, &
     3,2,3,3,3,3,1,3,3,3,3,2, &
     2,3,3,3,3,3,3,1,3,3,2,3, &
     3,3,2,3,3,2,3,3,1,3,3,3, &
     3,3,3,2,2,3,3,3,3,1,3,3, &
     2,3,3,3,3,3,3,2,3,3,1,3, &
     3,2,3,3,3,3,2,3,3,3,3,1  &
     ],pInt),shape(BCC_INTERACTIONSLIPTWIN),order=[2,1])                                            !< Slip-twin interaction types for bcc
                                                                                                    !< 1: coplanar interaction
                                                                                                    !< 2: screw trace between slip system and twin habit plane (easy cross slip)
                                                                                                    !< 3: other interaction
 integer(pInt), dimension(LATTICE_HEX_NSLIP,LATTICE_HEX_NTWIN), parameter :: &
   HEX_INTERACTIONSLIPTWIN = reshape(int( [&
      1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! --> twin
      1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! |
      1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! |
    !                                                                                   v
      5, 5, 5, 5, 5, 5,   6, 6, 6, 6, 6, 6,   7, 7, 7, 7, 7, 7,   8, 8, 8, 8, 8, 8, & ! slip
      5, 5, 5, 5, 5, 5,   6, 6, 6, 6, 6, 6,   7, 7, 7, 7, 7, 7,   8, 8, 8, 8, 8, 8, &
      5, 5, 5, 5, 5, 5,   6, 6, 6, 6, 6, 6,   7, 7, 7, 7, 7, 7,   8, 8, 8, 8, 8, 8, &
    !
      9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
      9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
      9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
    !
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
     13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
    !
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
     17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
    !
     21,21,21,21,21,21,  22,22,22,22,22,22,  23,23,23,23,23,23,  24,24,24,24,24,24, &
     21,21,21,21,21,21,  22,22,22,22,22,22,  23,23,23,23,23,23,  24,24,24,24,24,24, &
     21,21,21,21,21,21,  22,22,22,22,22,22,  23,23,23,23,23,23,  24,24,24,24,24,24, &
     21,21,21,21,21,21,  22,22,22,22,22,22,  23,23,23,23,23,23,  24,24,24,24,24,24, &
     21,21,21,21,21,21,  22,22,22,22,22,22,  23,23,23,23,23,23,  24,24,24,24,24,24, &
     21,21,21,21,21,21,  22,22,22,22,22,22,  23,23,23,23,23,23,  24,24,24,24,24,24  &
    !
     ],pInt),shape(HEX_INTERACTIONSLIPTWIN),order=[2,1])                                            !< Slip-twin interaction types for hex

 if (len_trim(structure) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='lattice_interaction_SlipTwin: '//trim(structure))
   
 select case(structure(1:3))
   case('fcc')
     interactionTypes = FCC_INTERACTIONSLIPTWIN
     NslipMax         = LATTICE_FCC_NSLIPSYSTEM
     NtwinMax         = LATTICE_FCC_NTWINSYSTEM
   case('bcc')
     interactionTypes = BCC_INTERACTIONSLIPTWIN
     NslipMax         = LATTICE_BCC_NSLIPSYSTEM
     NtwinMax         = LATTICE_BCC_NTWINSYSTEM
   case('hex')
     interactionTypes = HEX_INTERACTIONSLIPTWIN
     NslipMax         = LATTICE_HEX_NSLIPSYSTEM
     NtwinMax         = LATTICE_HEX_NTWINSYSTEM
   case default
     call IO_error(137_pInt,ext_msg='lattice_interaction_SlipTwin: '//trim(structure))
 end select

 interactionMatrix = buildInteraction(Nslip,Ntwin,NslipMax,NtwinMax,interactionValues,interactionTypes)

end function lattice_interaction_SlipTwin


!--------------------------------------------------------------------------------------------------
!> @brief Slip-trans interaction matrix
!> details only active slip and trans systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_interaction_SlipTrans(Nslip,Ntrans,interactionValues,structure) result(interactionMatrix)
 use IO, only: &
   IO_error

 implicit none
 integer(pInt),    dimension(:),                   intent(in) :: Nslip, &                           !< number of active slip systems per family
                                                                 Ntrans                             !< number of active trans systems per family
 real(pReal),      dimension(:),                   intent(in) :: interactionValues                  !< values for slip-trans interaction
 character(len=*),                                 intent(in) :: &
   structure                                                                                        !< lattice structure (parent crystal)
 real(pReal),     dimension(sum(Nslip),sum(Ntrans))           :: interactionMatrix

 integer(pInt),   dimension(:),                   allocatable :: NslipMax, &
                                                                 NtransMax
 integer(pInt),   dimension(:,:),                 allocatable :: interactionTypes

 integer(pInt), dimension(LATTICE_FCC_NSLIP,LATTICE_fcc_Ntrans), parameter :: &
   FCC_INTERACTIONSLIPTRANS = reshape(int( [&
     1,1,1,3,3,3,2,2,2,3,3,3, & ! ---> trans
     1,1,1,3,3,3,3,3,3,2,2,2, & ! |
     1,1,1,2,2,2,3,3,3,3,3,3, & ! |
     3,3,3,1,1,1,3,3,3,2,2,2, & ! v slip
     3,3,3,1,1,1,2,2,2,3,3,3, &
     2,2,2,1,1,1,3,3,3,3,3,3, &
     2,2,2,3,3,3,1,1,1,3,3,3, &
     3,3,3,2,2,2,1,1,1,3,3,3, &
     3,3,3,3,3,3,1,1,1,2,2,2, &
     3,3,3,2,2,2,3,3,3,1,1,1, &
     2,2,2,3,3,3,3,3,3,1,1,1, &
     3,3,3,3,3,3,2,2,2,1,1,1, &

     4,4,4,4,4,4,4,4,4,4,4,4, &
     4,4,4,4,4,4,4,4,4,4,4,4, &
     4,4,4,4,4,4,4,4,4,4,4,4, &
     4,4,4,4,4,4,4,4,4,4,4,4, &
     4,4,4,4,4,4,4,4,4,4,4,4, &
     4,4,4,4,4,4,4,4,4,4,4,4  &
     ],pInt),shape(FCC_INTERACTIONSLIPTRANS),order=[2,1])                                           !< Slip-trans interaction types for fcc

 if (len_trim(structure) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='lattice_interaction_SlipTrans: '//trim(structure))

 select case(structure(1:3))
   case('fcc')
     interactionTypes = FCC_INTERACTIONSLIPTRANS
     NslipMax         = LATTICE_FCC_NSLIPSYSTEM
     NtransMax        = LATTICE_FCC_NTRANSSYSTEM
   case default
     call IO_error(137_pInt,ext_msg='lattice_interaction_SlipTrans: '//trim(structure))
 end select

 interactionMatrix = buildInteraction(Nslip,Ntrans,NslipMax,NtransMax,interactionValues,interactionTypes)

end function lattice_interaction_SlipTrans


!--------------------------------------------------------------------------------------------------
!> @brief Twin-slip interaction matrix
!> details only active twin and slip systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_interaction_TwinSlip(Ntwin,Nslip,interactionValues,structure) result(interactionMatrix)
 use IO, only: &
   IO_error

 implicit none
 integer(pInt),    dimension(:),                   intent(in) :: Ntwin, &                           !< number of active twin systems per family
                                                                 Nslip                              !< number of active slip systems per family
 real(pReal),      dimension(:),                   intent(in) :: interactionValues                  !< values for twin-twin interaction
 character(len=*),                                 intent(in) :: structure                          !< lattice structure
 real(pReal),     dimension(sum(Ntwin),sum(Nslip))            :: interactionMatrix

 integer(pInt),   dimension(:),                   allocatable :: NtwinMax, &
                                                                 NslipMax
 integer(pInt),   dimension(:,:),                 allocatable :: interactionTypes

 integer(pInt), dimension(LATTICE_FCC_NTWIN,LATTICE_FCC_NSLIP), parameter :: &
   FCC_INTERACTIONTWINSLIP = 1_pInt                                                                 !< Twin-Slip interaction types for fcc

 integer(pInt), dimension(LATTICE_BCC_NTWIN,LATTICE_BCC_NSLIP), parameter :: &
   BCC_INTERACTIONTWINSLIP = 1_pInt                                                                 !< Twin-slip interaction types for bcc

 integer(pInt), dimension(LATTICE_HEX_NTWIN,LATTICE_HEX_NSLIP), parameter :: &
   HEX_INTERACTIONTWINSLIP = reshape(int( [&
      1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! --> slip
      1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! |
      1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! |
      1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! v
      1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! twin
      1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, &
    !
      2, 2, 2,   6, 6, 6,  10,10,10,  14,14,14,14,14,14,  18,18,18,18,18,18,18,18,18,18,18,18,  22,22,22,22,22,22, &
      2, 2, 2,   6, 6, 6,  10,10,10,  14,14,14,14,14,14,  18,18,18,18,18,18,18,18,18,18,18,18,  22,22,22,22,22,22, &
      2, 2, 2,   6, 6, 6,  10,10,10,  14,14,14,14,14,14,  18,18,18,18,18,18,18,18,18,18,18,18,  22,22,22,22,22,22, &
      2, 2, 2,   6, 6, 6,  10,10,10,  14,14,14,14,14,14,  18,18,18,18,18,18,18,18,18,18,18,18,  22,22,22,22,22,22, &
      2, 2, 2,   6, 6, 6,  10,10,10,  14,14,14,14,14,14,  18,18,18,18,18,18,18,18,18,18,18,18,  22,22,22,22,22,22, &
      2, 2, 2,   6, 6, 6,  10,10,10,  14,14,14,14,14,14,  18,18,18,18,18,18,18,18,18,18,18,18,  22,22,22,22,22,22, &
    !
      3, 3, 3,   7, 7, 7,  11,11,11,  15,15,15,15,15,15,  19,19,19,19,19,19,19,19,19,19,19,19,  23,23,23,23,23,23, &
      3, 3, 3,   7, 7, 7,  11,11,11,  15,15,15,15,15,15,  19,19,19,19,19,19,19,19,19,19,19,19,  23,23,23,23,23,23, &
      3, 3, 3,   7, 7, 7,  11,11,11,  15,15,15,15,15,15,  19,19,19,19,19,19,19,19,19,19,19,19,  23,23,23,23,23,23, &
      3, 3, 3,   7, 7, 7,  11,11,11,  15,15,15,15,15,15,  19,19,19,19,19,19,19,19,19,19,19,19,  23,23,23,23,23,23, &
      3, 3, 3,   7, 7, 7,  11,11,11,  15,15,15,15,15,15,  19,19,19,19,19,19,19,19,19,19,19,19,  23,23,23,23,23,23, &
      3, 3, 3,   7, 7, 7,  11,11,11,  15,15,15,15,15,15,  19,19,19,19,19,19,19,19,19,19,19,19,  23,23,23,23,23,23, &
    !
      4, 4, 4,   8, 8, 8,  12,12,12,  16,16,16,16,16,16,  20,20,20,20,20,20,20,20,20,20,20,20,  24,24,24,24,24,24, &
      4, 4, 4,   8, 8, 8,  12,12,12,  16,16,16,16,16,16,  20,20,20,20,20,20,20,20,20,20,20,20,  24,24,24,24,24,24, &
      4, 4, 4,   8, 8, 8,  12,12,12,  16,16,16,16,16,16,  20,20,20,20,20,20,20,20,20,20,20,20,  24,24,24,24,24,24, &
      4, 4, 4,   8, 8, 8,  12,12,12,  16,16,16,16,16,16,  20,20,20,20,20,20,20,20,20,20,20,20,  24,24,24,24,24,24, &
      4, 4, 4,   8, 8, 8,  12,12,12,  16,16,16,16,16,16,  20,20,20,20,20,20,20,20,20,20,20,20,  24,24,24,24,24,24, &
      4, 4, 4,   8, 8, 8,  12,12,12,  16,16,16,16,16,16,  20,20,20,20,20,20,20,20,20,20,20,20,  24,24,24,24,24,24  &
     ],pInt),shape(HEX_INTERACTIONTWINSLIP),order=[2,1])                                            !< Twin-twin interaction types for hex
     
 if (len_trim(structure) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='lattice_interaction_TwinSlip: '//trim(structure))

 select case(structure(1:3))
   case('fcc')
     interactionTypes = FCC_INTERACTIONTWINSLIP
     NtwinMax         = LATTICE_FCC_NTWINSYSTEM
     NslipMax         = LATTICE_FCC_NSLIPSYSTEM
   case('bcc')
     interactionTypes = BCC_INTERACTIONTWINSLIP
     NtwinMax         = LATTICE_BCC_NTWINSYSTEM
     NslipMax         = LATTICE_BCC_NSLIPSYSTEM
   case('hex')
     interactionTypes = HEX_INTERACTIONTWINSLIP
     NtwinMax         = LATTICE_HEX_NTWINSYSTEM
     NslipMax         = LATTICE_HEX_NSLIPSYSTEM
   case default
     call IO_error(137_pInt,ext_msg='lattice_interaction_TwinSlip: '//trim(structure))
 end select

 interactionMatrix = buildInteraction(Ntwin,Nslip,NtwinMax,NslipMax,interactionValues,interactionTypes)

end function lattice_interaction_TwinSlip


!--------------------------------------------------------------------------------------------------
!> @brief Schmid matrix for slip
!> details only active slip systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_SchmidMatrix_slip(Nslip,structure,cOverA) result(SchmidMatrix)
 use prec, only: &
   tol_math_check
 use IO, only: &
   IO_error
 use math, only: &
   math_trace33, &
   math_tensorproduct33

 implicit none
 integer(pInt),    dimension(:),            intent(in) :: Nslip                                     !< number of active slip systems per family
 character(len=*),                          intent(in) :: structure                                 !< lattice structure
 real(pReal),                               intent(in) :: cOverA
 real(pReal),     dimension(3,3,sum(Nslip))            :: SchmidMatrix

 real(pReal),     dimension(3,3,sum(Nslip))            :: coordinateSystem
 real(pReal),     dimension(:,:), allocatable          :: slipSystems
 integer(pInt),   dimension(:),  allocatable           :: NslipMax
 integer(pInt) :: i
 
 if (len_trim(structure) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='lattice_SchmidMatrix_slip: '//trim(structure))

 select case(structure(1:3))
   case('fcc')
     NslipMax    = LATTICE_FCC_NSLIPSYSTEM
     slipSystems = LATTICE_FCC_SYSTEMSLIP
   case('bcc')
     NslipMax    = LATTICE_BCC_NSLIPSYSTEM
     slipSystems = LATTICE_BCC_SYSTEMSLIP
   case('hex')
     NslipMax    = LATTICE_HEX_NSLIPSYSTEM
     slipSystems = LATTICE_HEX_SYSTEMSLIP
   case('bct')
     NslipMax    = LATTICE_BCT_NSLIPSYSTEM
     slipSystems = LATTICE_BCT_SYSTEMSLIP
   case default
     call IO_error(137_pInt,ext_msg='lattice_SchmidMatrix_slip: '//trim(structure))
 end select

 if (any(NslipMax(1:size(Nslip)) - Nslip < 0_pInt)) &
   call IO_error(145_pInt,ext_msg='Nslip '//trim(structure))
 if (any(Nslip < 0_pInt)) &
   call IO_error(144_pInt,ext_msg='Nslip '//trim(structure))

 coordinateSystem = buildCoordinateSystem(Nslip,NslipMax,slipSystems,structure,cOverA)

 do i = 1, sum(Nslip)
   SchmidMatrix(1:3,1:3,i) = math_tensorproduct33(coordinateSystem(1:3,1,i),coordinateSystem(1:3,2,i))
   if (abs(math_trace33(SchmidMatrix(1:3,1:3,i))) > tol_math_check) &
     call IO_error(0_pInt,i,ext_msg = 'dilatational Schmid matrix for slip')
 enddo

end function lattice_SchmidMatrix_slip


!--------------------------------------------------------------------------------------------------
!> @brief Schmid matrix for twinning
!> details only active twin systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_SchmidMatrix_twin(Ntwin,structure,cOverA) result(SchmidMatrix)
 use prec, only: &
   tol_math_check
 use IO, only: &
   IO_error
 use math, only: &
   math_trace33, &
   math_tensorproduct33

 implicit none
 integer(pInt),    dimension(:),            intent(in) :: Ntwin                                     !< number of active twin systems per family
 character(len=*),                          intent(in) :: structure                                 !< lattice structure
 real(pReal),                               intent(in) :: cOverA                                    !< c/a ratio
 real(pReal),     dimension(3,3,sum(Ntwin))            :: SchmidMatrix

 real(pReal),     dimension(3,3,sum(Ntwin))            :: coordinateSystem
 real(pReal),     dimension(:,:), allocatable          :: twinSystems
 integer(pInt),   dimension(:),   allocatable          :: NtwinMax
 integer(pInt) :: i

 if (len_trim(structure) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='lattice_SchmidMatrix_twin: '//trim(structure))

 select case(structure(1:3))
   case('fcc')
     NtwinMax    = LATTICE_FCC_NTWINSYSTEM
     twinSystems = LATTICE_FCC_SYSTEMTWIN
   case('bcc')
     NtwinMax    = LATTICE_BCC_NTWINSYSTEM
     twinSystems = LATTICE_BCC_SYSTEMTWIN
   case('hex')
     NtwinMax    = LATTICE_HEX_NTWINSYSTEM
     twinSystems = LATTICE_HEX_SYSTEMTWIN
   case default
     call IO_error(137_pInt,ext_msg='lattice_SchmidMatrix_twin: '//trim(structure))
 end select

 if (any(NtwinMax(1:size(Ntwin)) - Ntwin < 0_pInt)) &
   call IO_error(145_pInt,ext_msg='Ntwin '//trim(structure))
 if (any(Ntwin < 0_pInt)) &
   call IO_error(144_pInt,ext_msg='Ntwin '//trim(structure))

 coordinateSystem = buildCoordinateSystem(Ntwin,NtwinMax,twinSystems,structure,cOverA)

 do i = 1, sum(Ntwin)
   SchmidMatrix(1:3,1:3,i) = math_tensorproduct33(coordinateSystem(1:3,1,i),coordinateSystem(1:3,2,i))
   if (abs(math_trace33(SchmidMatrix(1:3,1:3,i))) > tol_math_check) &
     call IO_error(0_pInt,i,ext_msg = 'dilatational Schmid matrix for twin')
 enddo

end function lattice_SchmidMatrix_twin


!--------------------------------------------------------------------------------------------------
!> @brief Schmid matrix for twinning
!> details only active twin systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_SchmidMatrix_trans(Ntrans,structure_target,cOverA,a_bcc,a_fcc) result(SchmidMatrix)
 use prec, only: &
   tol_math_check
 use IO, only: &
   IO_error
 use math, only: &
   math_trace33, &
   math_tensorproduct33

 implicit none
 integer(pInt),    dimension(:),            intent(in) :: Ntrans                                    !< number of active twin systems per family
 real(pReal),                               intent(in) :: cOverA                                    !< c/a ratio
 real(pReal),     dimension(3,3,sum(Ntrans))           :: SchmidMatrix

 character(len=*),                          intent(in) :: &
   structure_target                                                                                 !< lattice structure

 real(pReal),     dimension(3,3,sum(Ntrans))           :: devNull
 real(pReal) :: a_bcc, a_fcc

 if (len_trim(structure_target) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='lattice_SchmidMatrix_trans (target): '//trim(structure_target))
 if (structure_target(1:3) /= 'bcc' .and. structure_target(1:3) /= 'hex') &
   call IO_error(137_pInt,ext_msg='lattice_SchmidMatrix_trans (target): '//trim(structure_target))

 !ToDo: add checks for CoverA_trans,a_fcc,a_bcc

 call buildTransformationSystem(devNull,SchmidMatrix,Ntrans,cOverA,a_fcc,a_bcc)

 end function lattice_SchmidMatrix_trans


!--------------------------------------------------------------------------------------------------
!> @brief Schmid matrix for cleavage
!> details only active cleavage systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_SchmidMatrix_cleavage(Ncleavage,structure,cOverA) result(SchmidMatrix)
 use math, only: &
   math_tensorproduct33
 use IO, only: &
   IO_error

 implicit none
 integer(pInt),    dimension(:),            intent(in) :: Ncleavage                                 !< number of active cleavage systems per family
 character(len=*),                          intent(in) :: structure                                 !< lattice structure
 real(pReal),                               intent(in) :: cOverA                                    !< c/a ratio
 real(pReal),     dimension(3,3,3,sum(Ncleavage))      :: SchmidMatrix

 real(pReal),     dimension(3,3,sum(Ncleavage))        :: coordinateSystem
 real(pReal),     dimension(:,:), allocatable          :: cleavageSystems
 integer(pInt),   dimension(:), allocatable            :: NcleavageMax
 integer(pInt) :: i
 
 if (len_trim(structure) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='lattice_SchmidMatrix_cleavage: '//trim(structure))

 select case(structure(1:3))
   case('iso')
     NcleavageMax    = LATTICE_ISO_NCLEAVAGESYSTEM
     cleavageSystems = LATTICE_ISO_SYSTEMCLEAVAGE
   case('ort')
     NcleavageMax    = LATTICE_ORT_NCLEAVAGESYSTEM
     cleavageSystems = LATTICE_ORT_SYSTEMCLEAVAGE
   case('fcc')
     NcleavageMax    = LATTICE_FCC_NCLEAVAGESYSTEM
     cleavageSystems = LATTICE_FCC_SYSTEMCLEAVAGE
   case('bcc')
     NcleavageMax    = LATTICE_BCC_NCLEAVAGESYSTEM
     cleavageSystems = LATTICE_BCC_SYSTEMCLEAVAGE
   case('hex')
     NcleavageMax    = LATTICE_HEX_NCLEAVAGESYSTEM
     cleavageSystems = LATTICE_HEX_SYSTEMCLEAVAGE
   case default
     call IO_error(137_pInt,ext_msg='lattice_SchmidMatrix_cleavage: '//trim(structure))
 end select

 if (any(NcleavageMax(1:size(Ncleavage)) - Ncleavage < 0_pInt)) &
   call IO_error(145_pInt,ext_msg='Ncleavage '//trim(structure))
 if (any(Ncleavage < 0_pInt)) &
   call IO_error(144_pInt,ext_msg='Ncleavage '//trim(structure))

 coordinateSystem = buildCoordinateSystem(Ncleavage,NcleavageMax,cleavageSystems,structure,cOverA)

 do i = 1, sum(Ncleavage)
   SchmidMatrix(1:3,1:3,1,i) = math_tensorproduct33(coordinateSystem(1:3,1,i),coordinateSystem(1:3,2,i))
   SchmidMatrix(1:3,1:3,2,i) = math_tensorproduct33(coordinateSystem(1:3,3,i),coordinateSystem(1:3,2,i))
   SchmidMatrix(1:3,1:3,3,i) = math_tensorproduct33(coordinateSystem(1:3,2,i),coordinateSystem(1:3,2,i))
 enddo

end function lattice_SchmidMatrix_cleavage


!--------------------------------------------------------------------------------------------------
!> @brief Normal direction of slip systems (n)
!--------------------------------------------------------------------------------------------------
function lattice_slip_normal(Nslip,structure,cOverA) result(n)

 implicit none
 integer(pInt),    dimension(:),            intent(in) :: Nslip                                     !< number of active slip systems per family
 character(len=*),                          intent(in) :: structure                                 !< lattice structure
 real(pReal),                               intent(in) :: cOverA                                    !< c/a ratio
 real(pReal),     dimension(3,sum(Nslip))              :: n

 real(pReal),     dimension(3,3,sum(Nslip))            :: coordinateSystem
 
 coordinateSystem = coordinateSystem_slip(Nslip,structure,cOverA)
 n = coordinateSystem(1:3,2,1:sum(Nslip))

end function lattice_slip_normal


!--------------------------------------------------------------------------------------------------
!> @brief Slip direction of slip systems (|| b)
!> @details: t = b x n
!--------------------------------------------------------------------------------------------------
function lattice_slip_direction(Nslip,structure,cOverA) result(d)

 implicit none
 integer(pInt),    dimension(:),            intent(in) :: Nslip                                     !< number of active slip systems per family
 character(len=*),                          intent(in) :: structure                                 !< lattice structure
 real(pReal),                               intent(in) :: cOverA                                    !< c/a ratio
 real(pReal),     dimension(3,sum(Nslip))              :: d

 real(pReal),     dimension(3,3,sum(Nslip))            :: coordinateSystem
 
 coordinateSystem = coordinateSystem_slip(Nslip,structure,cOverA)
 d = coordinateSystem(1:3,1,1:sum(Nslip))

end function lattice_slip_direction


!--------------------------------------------------------------------------------------------------
!> @brief Transverse direction of slip systems (||t, t = b x n)
!--------------------------------------------------------------------------------------------------
function lattice_slip_transverse(Nslip,structure,cOverA) result(t)

 implicit none
 integer(pInt),    dimension(:),            intent(in) :: Nslip                                     !< number of active slip systems per family
 character(len=*),                          intent(in) :: structure                                 !< lattice structure
 real(pReal),                               intent(in) :: cOverA                                    !< c/a ratio
 real(pReal),     dimension(3,sum(Nslip))              :: t

 real(pReal),     dimension(3,3,sum(Nslip))            :: coordinateSystem
 
 coordinateSystem = coordinateSystem_slip(Nslip,structure,cOverA)
 t = coordinateSystem(1:3,3,1:sum(Nslip))

end function lattice_slip_transverse


!--------------------------------------------------------------------------------------------------
!> @brief Projection of the transverse direction onto the slip plane
!> @details: This projection is used to calculate forest hardening for edge dislocations and for 
! mode III failure (ToDo: MD I am not 100% sure about mode III)
!--------------------------------------------------------------------------------------------------
function slipProjection_transverse(Nslip,structure,cOverA) result(projection)
 use math, only: &
   math_mul3x3

 implicit none
 integer(pInt),    dimension(:),            intent(in) :: Nslip                                     !< number of active slip systems per family
 character(len=*),                          intent(in) :: structure                                 !< lattice structure
 real(pReal),                               intent(in) :: cOverA                                    !< c/a ratio
 real(pReal),     dimension(sum(Nslip),sum(Nslip))     :: projection

 real(pReal),     dimension(3,3,sum(Nslip))            :: coordinateSystem
 integer(pInt) :: i, j
 
 coordinateSystem = coordinateSystem_slip(Nslip,structure,cOverA)

 do i=1_pInt, sum(Nslip); do j=1_pInt, sum(Nslip)
   projection(i,j) = abs(math_mul3x3(coordinateSystem(1:3,2,i),coordinateSystem(1:3,3,j)))
 enddo; enddo

end function slipProjection_transverse


!--------------------------------------------------------------------------------------------------
!> @brief Projection of the slip direction onto the slip plane
!> @details: This projection is used to calculate forest hardening for screw dislocations and for 
! mode II failure (ToDo: MD I am not 100% sure about mode II)
!--------------------------------------------------------------------------------------------------
function slipProjection_direction(Nslip,structure,cOverA) result(projection)
 use math, only: &
   math_mul3x3

 implicit none
 integer(pInt),    dimension(:),            intent(in) :: Nslip                                     !< number of active slip systems per family
 character(len=*),                          intent(in) :: structure                                 !< lattice structure
 real(pReal),                               intent(in) :: cOverA                                    !< c/a ratio
 real(pReal),     dimension(sum(Nslip),sum(Nslip))     :: projection

 real(pReal),     dimension(3,3,sum(Nslip))            :: coordinateSystem
 integer(pInt) :: i, j
 
 coordinateSystem = coordinateSystem_slip(Nslip,structure,cOverA)

 do i=1_pInt, sum(Nslip); do j=1_pInt, sum(Nslip)
   projection(i,j) = abs(math_mul3x3(coordinateSystem(1:3,2,i),coordinateSystem(1:3,1,j)))
 enddo; enddo

end function slipProjection_direction


!--------------------------------------------------------------------------------------------------
!> @brief Projection of the slip plane onto itself
!> @details: This projection is used for mode I failure
!--------------------------------------------------------------------------------------------------
function slipProjection_normal(Nslip,structure,cOverA) result(projection)
 use math, only: &
   math_mul3x3

 implicit none
 integer(pInt),    dimension(:),            intent(in) :: Nslip                                     !< number of active slip systems per family
 character(len=*),                          intent(in) :: structure                                 !< lattice structure
 real(pReal),                               intent(in) :: cOverA                                    !< c/a ratio
 real(pReal),     dimension(sum(Nslip),sum(Nslip))     :: projection

 real(pReal),     dimension(3,3,sum(Nslip))            :: coordinateSystem
 integer(pInt) :: i, j

 coordinateSystem = coordinateSystem_slip(Nslip,structure,cOverA)

 do i=1_pInt, sum(Nslip); do j=1_pInt, sum(Nslip)
   projection(i,j) = abs(math_mul3x3(coordinateSystem(1:3,2,i),coordinateSystem(1:3,2,j)))
 enddo; enddo

end function slipProjection_normal


!--------------------------------------------------------------------------------------------------
!> @brief build a local coordinate system on slip systems
!> @details Order: Direction, plane (normal), and common perpendicular
!--------------------------------------------------------------------------------------------------
function coordinateSystem_slip(Nslip,structure,cOverA) result(coordinateSystem)
 use math, only: &
   math_mul3x3
 use IO, only: &
   IO_error

 implicit none
 integer(pInt),    dimension(:),            intent(in) :: Nslip                                     !< number of active slip systems per family
 character(len=*),                          intent(in) :: structure                                 !< lattice structure
 real(pReal),                               intent(in) :: cOverA                                    !< c/a ratio
 real(pReal),     dimension(3,3,sum(Nslip))            :: coordinateSystem
 real(pReal),     dimension(:,:), allocatable          :: slipSystems
 integer(pInt),   dimension(:), allocatable            :: NslipMax

 if (len_trim(structure) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='coordinateSystem_slip: '//trim(structure))

 select case(structure(1:3))
   case('fcc')
     NslipMax    = LATTICE_FCC_NSLIPSYSTEM
     slipSystems = LATTICE_FCC_SYSTEMSLIP
   case('bcc')
     NslipMax    = LATTICE_BCC_NSLIPSYSTEM
     slipSystems = LATTICE_BCC_SYSTEMSLIP
   case('hex')
     NslipMax    = LATTICE_HEX_NSLIPSYSTEM
     slipSystems = LATTICE_HEX_SYSTEMSLIP
   case('bct')
     NslipMax    = LATTICE_BCT_NSLIPSYSTEM
     slipSystems = LATTICE_BCT_SYSTEMSLIP
   case default
     call IO_error(137_pInt,ext_msg='coordinateSystem_slip: '//trim(structure))
 end select

 if (any(NslipMax(1:size(Nslip)) - Nslip < 0_pInt)) &
   call IO_error(145_pInt,ext_msg='Nslip '//trim(structure))
 if (any(Nslip < 0_pInt)) &
   call IO_error(144_pInt,ext_msg='Nslip '//trim(structure))

 coordinateSystem = buildCoordinateSystem(Nslip,NslipMax,slipSystems,structure,cOverA)

end function coordinateSystem_slip


!--------------------------------------------------------------------------------------------------
!> @brief Populates reduced interaction matrix
!--------------------------------------------------------------------------------------------------
function buildInteraction(activeA,activeB,maxA,maxB,values,matrix)
 use IO, only: &
   IO_error
 implicit none
 integer(pInt), dimension(:),           intent(in) :: &
    activeA, &                                                                                      !< number of active systems as specified in material.config
    activeB, &                                                                                      !< number of active systems as specified in material.config
    maxA, &                                                                                         !< number of maximum available systems
    maxB                                                                                            !< number of maximum available systems
 real(pReal),   dimension(:),           intent(in)   :: values                                      !< interaction values
 integer(pInt), dimension(:,:),         intent(in)   :: matrix                                      !< complete interaction matrix
 real(pReal),   dimension(sum(activeA),sum(activeB)) :: buildInteraction

 integer(pInt) :: &
   index_myFamily, index_otherFamily, &
   mf, ms, of, os

 myFamilies: do mf = 1_pInt,size(activeA,1)
   index_myFamily = sum(activeA(1:mf-1_pInt))
   mySystems: do ms = 1_pInt,activeA(mf)

     otherFamilies: do of = 1_pInt,size(activeB,1)
     index_otherFamily = sum(activeB(1:of-1_pInt))
       otherSystems: do os = 1_pInt,activeB(of)
         if(matrix(sum(maxA(1:mf-1))+ms, sum(maxB(1:of-1))+os) > size(values)) &
           call IO_error(138,ext_msg='buildInteraction')
         buildInteraction(index_myFamily+ms,index_otherFamily+os) = &
                                       values(matrix(sum(maxA(1:mf-1))+ms, sum(maxB(1:of-1))+os))
     enddo otherSystems; enddo otherFamilies;

 enddo mySystems;enddo myFamilies

end function buildInteraction


!--------------------------------------------------------------------------------------------------
!> @brief build a local coordinate system on slip, twin, trans, cleavage systems
!> @details Order: Direction, plane (normal), and common perpendicular
!--------------------------------------------------------------------------------------------------
function buildCoordinateSystem(active,complete,system,structure,cOverA)
 use IO, only: &
   IO_error
 use math, only: &
   math_crossproduct

 implicit none
 integer(pInt), dimension(:), intent(in) :: &
   active, &
   complete
 real(pReal), dimension(:,:), intent(in) :: &
   system
 character(len=*),            intent(in) :: &
   structure                                                                                        !< lattice structure
 real(pReal),                 intent(in) :: &
   cOverA
 real(pReal), dimension(3,3,sum(active)) :: &
  buildCoordinateSystem

 real(pReal), dimension(3) :: &
   direction, normal
 integer(pInt) :: &
   a, &                                                                                             !< index of active system
   c, &                                                                                             !< index in complete system matrix
   f, &                                                                                             !< index of my family
   s                                                                                                !< index of my system in current family

 if (len_trim(structure) /= 3_pInt) &
   call IO_error(137_pInt,ext_msg='buildCoordinateSystem: '//trim(structure))
 if (trim(structure(1:3)) == 'bct' .and. cOverA > 2.0_pReal) &
   call IO_error(131_pInt,ext_msg='buildCoordinateSystem:'//trim(structure))
 if (trim(structure(1:3)) == 'hex' .and. (cOverA < 1.0_pReal .or. cOverA > 2.0_pReal)) &
   call IO_error(131_pInt,ext_msg='buildCoordinateSystem:'//trim(structure))

 a = 0_pInt
 activeFamilies: do f = 1_pInt,size(active,1)
   activeSystems: do s = 1_pInt,active(f)
     a = a + 1_pInt
     c = sum(complete(1:f-1))+s

     select case(trim(structure(1:3)))

       case ('fcc','bcc','iso','ort','bct')
         direction = system(1:3,c)
         normal    = system(4:6,c)

       case ('hex')
         direction = [ system(1,c)*1.5_pReal, &
                      (system(1,c)+2.0_pReal*system(2,c))*sqrt(0.75_pReal), &
                       system(4,c)*cOverA ]                                                         ! direction [uvtw]->[3u/2 (u+2v)*sqrt(3)/2 w*(c/a)])
         normal    = [ system(5,c), &
                      (system(5,c)+2.0_pReal*system(6,c))/sqrt(3.0_pReal), &
                       system(8,c)/cOverA ]                                                         ! plane (hkil)->(h (h+2k)/sqrt(3) l/(c/a))

       case default
         call IO_error(137_pInt,ext_msg='buildCoordinateSystem: '//trim(structure))

     end select

     buildCoordinateSystem(1:3,1,a) = direction/norm2(direction)
     buildCoordinateSystem(1:3,2,a) = normal/norm2(normal)
     buildCoordinateSystem(1:3,3,a) = math_crossproduct(buildCoordinateSystem(1:3,1,a),&
                                                        buildCoordinateSystem(1:3,2,a))

   enddo activeSystems
 enddo activeFamilies

end function buildCoordinateSystem


!--------------------------------------------------------------------------------------------------
!> @brief Helper function to define transformation systems
! Needed to calculate Schmid matrix and rotated stiffness matrices.
! @details: set c/a   = 0.0 for fcc -> bcc transformation
!           set a_bcc = 0.0 for fcc -> hex transformation
!--------------------------------------------------------------------------------------------------
subroutine buildTransformationSystem(Q,S,Ntrans,cOverA,a_fcc,a_bcc)
 use prec, only: &
   dEq0
 use math, only: &
   math_crossproduct, &
   math_tensorproduct33, &
   math_mul33x33, &
   math_mul33x3, &
   math_axisAngleToR, &
   INRAD, &
   MATH_I3
 use IO, only: &
   IO_error

 implicit none
 integer(pInt), dimension(:), intent(in) :: &
   Ntrans
 real(pReal),  dimension(3,3,sum(Ntrans)), intent(out) :: &
   Q, &                                                                                             !< Total rotation: Q = R*B
   S                                                                                                !< Eigendeformation tensor for phase transformation
 real(pReal),                 intent(in) :: &
   cOverA, &                                                                                        !< c/a for target hex structure
   a_bcc, &                                                                                         !< lattice parameter a for target bcc structure
   a_fcc                                                                                            !< lattice parameter a for parent fcc structure

 real(pReal), dimension(3,3) :: &
   R, &                                                                                             !< Pitsch rotation
   U, &                                                                                             !< Bain deformation
   B, &                                                                                             !< Rotation of fcc to Bain coordinate system
   ss, sd
 real(pReal), dimension(3) :: &
   x, y, z
 integer(pInt) :: &
   i
 real(pReal), dimension(3+3,LATTICE_FCC_NTRANS), parameter :: &
   LATTICE_FCCTOHEX_SYSTEMTRANS = reshape(real( [&
     -2, 1, 1,     1, 1, 1, &
      1,-2, 1,     1, 1, 1, &
      1, 1,-2,     1, 1, 1, &
      2,-1, 1,    -1,-1, 1, &
     -1, 2, 1,    -1,-1, 1, &
     -1,-1,-2,    -1,-1, 1, &
     -2,-1,-1,     1,-1,-1, &
      1, 2,-1,     1,-1,-1, &
      1,-1, 2,     1,-1,-1, &
      2, 1,-1,    -1, 1,-1, &
     -1,-2,-1,    -1, 1,-1, &
     -1, 1, 2,    -1, 1,-1  &
     ],pReal),shape(LATTICE_FCCTOHEX_SYSTEMTRANS))
      real(pReal), dimension(4,LATTICE_fcc_Ntrans), parameter :: &
   LATTICE_FCCTOBCC_SYSTEMTRANS = reshape([&
     0.0, 1.0, 0.0,     10.26, &                                                                    ! Pitsch OR (Ma & Hartmaier 2014, Table 3)
     0.0, 1.0, 0.0,    -10.26, &
     0.0, 0.0, 1.0,     10.26, &
     0.0, 0.0, 1.0,    -10.26, &
     1.0, 0.0, 0.0,     10.26, &
     1.0, 0.0, 0.0,    -10.26, &
     0.0, 0.0, 1.0,     10.26, &
     0.0, 0.0, 1.0,    -10.26, &
     1.0, 0.0, 0.0,     10.26, &
     1.0, 0.0, 0.0,    -10.26, &
     0.0, 1.0, 0.0,     10.26, &
     0.0, 1.0, 0.0,    -10.26  &
     ],shape(LATTICE_FCCTOBCC_SYSTEMTRANS))

 integer(pInt), dimension(9,LATTICE_fcc_Ntrans), parameter :: &
   LATTICE_FCCTOBCC_BAINVARIANT = reshape(int( [&
     1, 0, 0, 0, 1, 0, 0, 0, 1, &                                                                   ! Pitsch OR (Ma & Hartmaier 2014, Table 3)
     1, 0, 0, 0, 1, 0, 0, 0, 1, &
     1, 0, 0, 0, 1, 0, 0, 0, 1, &
     1, 0, 0, 0, 1, 0, 0, 0, 1, &
     0, 1, 0, 1, 0, 0, 0, 0, 1, &
     0, 1, 0, 1, 0, 0, 0, 0, 1, &
     0, 1, 0, 1, 0, 0, 0, 0, 1, &
     0, 1, 0, 1, 0, 0, 0, 0, 1, &
     0, 0, 1, 1, 0, 0, 0, 1, 0, &
     0, 0, 1, 1, 0, 0, 0, 1, 0, &
     0, 0, 1, 1, 0, 0, 0, 1, 0, &
     0, 0, 1, 1, 0, 0, 0, 1, 0  &
     ],pInt),shape(LATTICE_FCCTOBCC_BAINVARIANT))

 real(pReal), dimension(4,LATTICE_fcc_Ntrans), parameter :: &
   LATTICE_FCCTOBCC_BAINROT = reshape([&
     1.0, 0.0, 0.0,     45.0, &                                                                     ! Rotate fcc austensite to bain variant
     1.0, 0.0, 0.0,     45.0, &
     1.0, 0.0, 0.0,     45.0, &
     1.0, 0.0, 0.0,     45.0, &
     0.0, 1.0, 0.0,     45.0, &
     0.0, 1.0, 0.0,     45.0, &
     0.0, 1.0, 0.0,     45.0, &
     0.0, 1.0, 0.0,     45.0, &
     0.0, 0.0, 1.0,     45.0, &
     0.0, 0.0, 1.0,     45.0, &
     0.0, 0.0, 1.0,     45.0, &
     0.0, 0.0, 1.0,     45.0  &
     ],shape(LATTICE_FCCTOBCC_BAINROT))

 if (size(Ntrans) < 1_pInt .or. size(Ntrans) > 1_pInt) print*, 'mist'                               ! ToDo

 if (a_bcc > 0.0_pReal .and. dEq0(cOverA)) then                                                     ! fcc -> bcc transformation
   do i = 1_pInt,sum(Ntrans)
     R = math_axisAngleToR(lattice_fccTobcc_systemTrans(1:3,i), &
                           lattice_fccTobcc_systemTrans(4,i)*INRAD)
     B = math_axisAngleToR(lattice_fccTobcc_bainRot(1:3,i), & 
                           lattice_fccTobcc_bainRot(4,i)*INRAD)
     x = real(LATTICE_fccTobcc_bainVariant(1:3,i),pReal)
     y = real(LATTICE_fccTobcc_bainVariant(4:6,i),pReal)
     z = real(LATTICE_fccTobcc_bainVariant(7:9,i),pReal)

     U = (a_bcc/a_fcc)*math_tensorproduct33(x,x) &
       + (a_bcc/a_fcc)*math_tensorproduct33(y,y) * sqrt(2.0_pReal) &
       + (a_bcc/a_fcc)*math_tensorproduct33(z,z) * sqrt(2.0_pReal)
     Q(1:3,1:3,i) = math_mul33x33(R,B)
     S(1:3,1:3,i) = math_mul33x33(R,U) - MATH_I3
   enddo
 elseif (cOverA > 0.0_pReal .and. dEq0(a_bcc)) then                                                 ! fcc -> hex transformation
   ss      = MATH_I3
   sd      = MATH_I3
   ss(1,3) = sqrt(2.0_pReal)/4.0_pReal
   if (cOverA > 1.0_pReal .and. cOverA < 2.0_pReal) &
     sd(3,3) = cOverA/sqrt(8.0_pReal/3.0_pReal)

   do i = 1_pInt,sum(Ntrans)
     x = lattice_fccTohex_systemTrans(1:3,i)/norm2(lattice_fccTohex_systemTrans(1:3,i))
     z = lattice_fccTohex_systemTrans(4:6,i)/norm2(lattice_fccTohex_systemTrans(4:6,i))
     y = -math_crossproduct(x,z)
     Q(1:3,1,i) = x
     Q(1:3,2,i) = y
     Q(1:3,3,i) = z
     S(1:3,1:3,i) = math_mul33x33(Q(1:3,1:3,i), math_mul33x33(math_mul33x33(sd,ss), transpose(Q(1:3,1:3,i)))) - MATH_I3  ! ToDo: This is of interest for the Schmid matrix only
   enddo
 else
   call IO_error(0_pInt) !ToDo: define error
 endif

end subroutine buildTransformationSystem

end module lattice
