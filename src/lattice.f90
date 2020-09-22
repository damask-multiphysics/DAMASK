!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  contains lattice structure definitions including Schmid matrices for slip, twin, trans,
!          and cleavage as well as interaction among the various systems
!--------------------------------------------------------------------------------------------------
module lattice
  use prec
  use IO
  use config
  use math
  use rotations

  implicit none
  private

!--------------------------------------------------------------------------------------------------
! face centered cubic
  integer, dimension(*), parameter :: &
    FCC_NSLIPSYSTEM = [12, 6]                                                                       !< # of slip systems per family for fcc

  integer, dimension(*), parameter :: &
    FCC_NTWINSYSTEM = [12]                                                                          !< # of twin systems per family for fcc

  integer, dimension(*), parameter :: &
    FCC_NTRANSSYSTEM = [12]                                                                         !< # of transformation systems per family for fcc

  integer, dimension(*), parameter :: &
    FCC_NCLEAVAGESYSTEM = [3]                                                                       !< # of cleavage systems per family for fcc

  integer, parameter  :: &
#ifndef __PGI
    FCC_NSLIP     = sum(FCC_NSLIPSYSTEM), &                                                         !< total # of slip systems for fcc
    FCC_NTWIN     = sum(FCC_NTWINSYSTEM), &                                                         !< total # of twin systems for fcc
    FCC_NTRANS    = sum(FCC_NTRANSSYSTEM), &                                                        !< total # of transformation systems for fcc
    FCC_NCLEAVAGE = sum(FCC_NCLEAVAGESYSTEM)                                                        !< total # of cleavage systems for fcc
#else
    FCC_NSLIP     = 18, &
    FCC_NTWIN     = 12, &
    FCC_NTRANS    = 12, &
    FCC_NCLEAVAGE = 3
#endif

  real(pReal), dimension(3+3,FCC_NSLIP), parameter :: &
    FCC_SYSTEMSLIP = reshape(real([&
     ! Slip direction     Plane normal                                                              ! SCHMID-BOAS notation
       0, 1,-1,     1, 1, 1, &                                                                      ! B2
      -1, 0, 1,     1, 1, 1, &                                                                      ! B4
       1,-1, 0,     1, 1, 1, &                                                                      ! B5
       0,-1,-1,    -1,-1, 1, &                                                                      ! C1
       1, 0, 1,    -1,-1, 1, &                                                                      ! C3
      -1, 1, 0,    -1,-1, 1, &                                                                      ! C5
       0,-1, 1,     1,-1,-1, &                                                                      ! A2
      -1, 0,-1,     1,-1,-1, &                                                                      ! A3
       1, 1, 0,     1,-1,-1, &                                                                      ! A6
       0, 1, 1,    -1, 1,-1, &                                                                      ! D1
       1, 0,-1,    -1, 1,-1, &                                                                      ! D4
      -1,-1, 0,    -1, 1,-1, &                                                                      ! D6
      ! Slip system <110>{110}
       1, 1, 0,     1,-1, 0, &
       1,-1, 0,     1, 1, 0, &
       1, 0, 1,     1, 0,-1, &
       1, 0,-1,     1, 0, 1, &
       0, 1, 1,     0, 1,-1, &
       0, 1,-1,     0, 1, 1  &
      ],pReal),shape(FCC_SYSTEMSLIP))                                                               !< Slip system <110>{111} directions. Sorted according to Eisenlohr & Hantcherli

  real(pReal), dimension(3+3,FCC_NTWIN), parameter :: &
    FCC_SYSTEMTWIN = reshape(real( [&
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
      ],pReal),shape(FCC_SYSTEMTWIN))                                                               !< Twin system <112>{111} directions. Sorted according to Eisenlohr & Hantcherli

  integer, dimension(2,FCC_NTWIN), parameter, public :: &
    lattice_FCC_TWINNUCLEATIONSLIPPAIR = reshape( [&
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
      ],shape(lattice_FCC_TWINNUCLEATIONSLIPPAIR))

  real(pReal), dimension(3+3,FCC_NCLEAVAGE), parameter :: &
    FCC_SYSTEMCLEAVAGE = reshape(real([&
     ! Cleavage direction     Plane normal
       0, 1, 0,     1, 0, 0, &
       0, 0, 1,     0, 1, 0, &
       1, 0, 0,     0, 0, 1  &
      ],pReal),shape(FCC_SYSTEMCLEAVAGE))

!--------------------------------------------------------------------------------------------------
! body centered cubic
  integer, dimension(*), parameter :: &
    BCC_NSLIPSYSTEM = [12, 12]                                                                      !< # of slip systems per family for bcc

  integer, dimension(*), parameter :: &
    BCC_NTWINSYSTEM = [12]                                                                          !< # of twin systems per family for bcc

  integer, dimension(*), parameter :: &
    BCC_NCLEAVAGESYSTEM = [3]                                                                       !< # of cleavage systems per family for bcc

  integer, parameter  :: &
#ifndef __PGI
    BCC_NSLIP     = sum(BCC_NSLIPSYSTEM), &                                                         !< total # of slip systems for bcc
    BCC_NTWIN     = sum(BCC_NTWINSYSTEM), &                                                         !< total # of twin systems for bcc
    BCC_NCLEAVAGE = sum(BCC_NCLEAVAGESYSTEM)                                                        !< total # of cleavage systems for bcc
#else
    BCC_NSLIP     = 24, &
    BCC_NTWIN     = 12, &
    BCC_NCLEAVAGE = 3
#endif

  real(pReal), dimension(3+3,BCC_NSLIP), parameter :: &
    BCC_SYSTEMSLIP = reshape(real([&
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
      ],pReal),shape(BCC_SYSTEMSLIP))

  real(pReal), dimension(3+3,BCC_NTWIN), parameter :: &
    BCC_SYSTEMTWIN = reshape(real([&
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
      ],pReal),shape(BCC_SYSTEMTWIN))

  real(pReal), dimension(3+3,BCC_NCLEAVAGE), parameter :: &
    BCC_SYSTEMCLEAVAGE = reshape(real([&
     ! Cleavage direction     Plane normal
       0, 1, 0,     1, 0, 0, &
       0, 0, 1,     0, 1, 0, &
       1, 0, 0,     0, 0, 1  &
      ],pReal),shape(BCC_SYSTEMCLEAVAGE))

!--------------------------------------------------------------------------------------------------
! hexagonal
  integer, dimension(*), parameter :: &
    HEX_NSLIPSYSTEM = [3, 3, 3, 6, 12, 6]                                                           !< # of slip systems per family for hex

  integer, dimension(*), parameter :: &
    HEX_NTWINSYSTEM = [6, 6, 6, 6]                                                                  !< # of slip systems per family for hex

  integer, parameter  :: &
#ifndef __PGI
    HEX_NSLIP     = sum(HEX_NSLIPSYSTEM), &                                                         !< total # of slip systems for hex
    HEX_NTWIN     = sum(HEX_NTWINSYSTEM)                                                            !< total # of twin systems for hex
#else
    HEX_NSLIP     = 33, &
    HEX_NTWIN     = 24
#endif

  real(pReal), dimension(4+4,HEX_NSLIP), parameter :: &
    HEX_SYSTEMSLIP = reshape(real([&
     ! Slip direction     Plane normal
     ! Basal systems <-1-1.0>{00.1} (independent of c/a-ratio, Bravais notation (4 coordinate base))
       2, -1, -1,  0,     0,  0,  0,  1, &
      -1,  2, -1,  0,     0,  0,  0,  1, &
      -1, -1,  2,  0,     0,  0,  0,  1, &
     ! 1st type prismatic systems <-1-1.0>{1-1.0}  (independent of c/a-ratio)
       2, -1, -1,  0,     0,  1, -1,  0, &
      -1,  2, -1,  0,    -1,  0,  1,  0, &
      -1, -1,  2,  0,     1, -1,  0,  0, &
     ! 2nd type prismatic systems <-11.0>{11.0} -- a slip; plane normals independent of c/a-ratio
      -1,  1,  0,  0,     1,  1, -2,  0, &
       0, -1,  1,  0,    -2,  1,  1,  0, &
       1,  0, -1,  0,     1, -2,  1,  0, &
     ! 1st type 1st order pyramidal systems <-1-1.0>{-11.1} -- plane normals depend on the c/a-ratio
      -1,  2, -1,  0,     1,  0, -1,  1, &
      -2,  1,  1,  0,     0,  1, -1,  1, &
      -1, -1,  2,  0,    -1,  1,  0,  1, &
       1, -2,  1,  0,    -1,  0,  1,  1, &
       2, -1, -1,  0,     0, -1,  1,  1, &
       1,  1, -2,  0,     1, -1,  0,  1, &
     ! pyramidal system: c+a slip <11.3>{-10.1} -- plane normals depend on the c/a-ratio
      -2,  1,  1,  3,     1,  0, -1,  1, &
      -1, -1,  2,  3,     1,  0, -1,  1, &
      -1, -1,  2,  3,     0,  1, -1,  1, &
       1, -2,  1,  3,     0,  1, -1,  1, &
       1, -2,  1,  3,    -1,  1,  0,  1, &
       2, -1, -1,  3,    -1,  1,  0,  1, &
       2, -1, -1,  3,    -1,  0,  1,  1, &
       1,  1, -2,  3,    -1,  0,  1,  1, &
       1,  1, -2,  3,     0, -1,  1,  1, &
      -1,  2, -1,  3,     0, -1,  1,  1, &
      -1,  2, -1,  3,     1, -1,  0,  1, &
      -2,  1,  1,  3,     1, -1,  0,  1, &
     ! pyramidal system: c+a slip <11.3>{-1-1.2} -- as for hexagonal ice (Castelnau et al. 1996, similar to twin system found below)
      -1, -1,  2,  3,     1,  1, -2,  2, & ! <11.3>{-1-1.2} shear = 2((c/a)^2-2)/(3 c/a)
       1, -2,  1,  3,    -1,  2, -1,  2, &
       2, -1, -1,  3,    -2,  1,  1,  2, &
       1,  1, -2,  3,    -1, -1,  2,  2, &
      -1,  2, -1,  3,     1, -2,  1,  2, &
      -2,  1,  1,  3,     2, -1, -1,  2  &
      ],pReal),shape(HEX_SYSTEMSLIP))                                                               !< slip systems for hex, sorted by P. Eisenlohr CCW around <c> starting next to a_1 axis

  real(pReal), dimension(4+4,HEX_NTWIN), parameter :: &
    HEX_SYSTEMTWIN =  reshape(real([&
     ! Compression or Tension = f(twinning shear=f(c/a)) for each metal ! (according to Yoo 1981)
      -1,  0,  1,  1,     1,  0, -1,  2, & ! <-10.1>{10.2} shear = (3-(c/a)^2)/(sqrt(3) c/a)
       0, -1,  1,  1,     0,  1, -1,  2, &
       1, -1,  0,  1,    -1,  1,  0,  2, &
       1,  0, -1,  1,    -1,  0,  1,  2, &
       0,  1, -1,  1,     0, -1,  1,  2, &
      -1,  1,  0,  1,     1, -1,  0,  2, &
!
      -1, -1,  2,  6,     1,  1, -2,  1, & ! <11.6>{-1-1.1} shear = 1/(c/a)
       1, -2,  1,  6,    -1,  2, -1,  1, &
       2, -1, -1,  6,    -2,  1,  1,  1, &
       1,  1, -2,  6,    -1, -1,  2,  1, &
      -1,  2, -1,  6,     1, -2,  1,  1, &
      -2,  1,  1,  6,     2, -1, -1,  1, &
!
       1,  0, -1, -2,     1,  0, -1,  1, & ! <10.-2>{10.1} shear = (4(c/a)^2-9)/(4 sqrt(3) c/a)
       0,  1, -1, -2,     0,  1, -1,  1, &
      -1,  1,  0, -2,    -1,  1,  0,  1, &
      -1,  0,  1, -2,    -1,  0,  1,  1, &
       0, -1,  1, -2,     0, -1,  1,  1, &
       1, -1,  0, -2,     1, -1,  0,  1, &
!
       1,  1, -2, -3,     1,  1, -2,  2, & ! <11.-3>{11.2} shear = 2((c/a)^2-2)/(3 c/a)
      -1,  2, -1, -3,    -1,  2, -1,  2, &
      -2,  1,  1, -3,    -2,  1,  1,  2, &
      -1, -1,  2, -3,    -1, -1,  2,  2, &
       1, -2,  1, -3,     1, -2,  1,  2, &
       2, -1, -1, -3,     2, -1, -1,  2  &
      ],pReal),shape(HEX_SYSTEMTWIN))                                                               !< twin systems for hex, sorted by P. Eisenlohr CCW around <c> starting next to a_1 axis

!--------------------------------------------------------------------------------------------------
! body centered tetragonal
  integer, dimension(*), parameter :: &
    BCT_NSLIPSYSTEM = [2, 2, 2, 4, 2, 4, 2, 2, 4, 8, 4, 8, 8 ]                                      !< # of slip systems per family for bct (Sn) Bieler J. Electr Mater 2009

  integer, parameter :: &
#ifndef __PGI
    BCT_NSLIP = sum(BCT_NSLIPSYSTEM)                                                                !< total # of slip systems for bct
#else
    BCT_NSLIP = 52
#endif

  real(pReal), dimension(3+3,BCT_NSLIP), parameter :: &
    BCT_SYSTEMSLIP = reshape(real([&
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
       ],pReal),shape(BCT_SYSTEMSLIP))                                                              !< slip systems for bct sorted by Bieler

!--------------------------------------------------------------------------------------------------
! orthorhombic
  integer, dimension(*), parameter :: &
    ORT_NCLEAVAGESYSTEM = [1, 1, 1]                                                                 !< # of cleavage systems per family for ortho

  integer, parameter  :: &
#ifndef __PGI
    ORT_NCLEAVAGE = sum(ORT_NCLEAVAGESYSTEM)                                                        !< total # of cleavage systems for ortho
#else
    ORT_NCLEAVAGE = 3
#endif

  real(pReal), dimension(3+3,ORT_NCLEAVAGE), parameter :: &
    ORT_SYSTEMCLEAVAGE = reshape(real([&
     ! Cleavage direction     Plane normal
       0, 1, 0,     1, 0, 0, &
       0, 0, 1,     0, 1, 0, &
       1, 0, 0,     0, 0, 1  &
      ],pReal),shape(ORT_SYSTEMCLEAVAGE))


  enum, bind(c); enumerator :: &
    lattice_UNDEFINED_ID, &
    lattice_ISO_ID, &
    lattice_FCC_ID, &
    lattice_BCC_ID, &
    lattice_BCT_ID, &
    lattice_HEX_ID, &
    lattice_ORT_ID
  end enum

! SHOULD NOT BE PART OF LATTICE BEGIN
  real(pReal),                        dimension(:),     allocatable, public, protected :: &
    lattice_mu, lattice_nu, &
    lattice_damageMobility, &
    lattice_massDensity, &
    lattice_specificHeat
   real(pReal),                       dimension(:,:,:), allocatable, public, protected :: &
    lattice_C66, &
    lattice_thermalConductivity, &
    lattice_damageDiffusion
 integer(kind(lattice_UNDEFINED_ID)), dimension(:),     allocatable, public, protected :: &
    lattice_structure
! SHOULD NOT BE PART OF LATTICE END

  interface lattice_forestProjection_edge
    module procedure slipProjection_transverse
  end interface lattice_forestProjection_edge

  interface lattice_forestProjection_screw
    module procedure slipProjection_direction
  end interface lattice_forestProjection_screw

  public :: &
    lattice_init, &
    lattice_ISO_ID, &
    lattice_FCC_ID, &
    lattice_BCC_ID, &
    lattice_BCT_ID, &
    lattice_HEX_ID, &
    lattice_ORT_ID, &
    lattice_applyLatticeSymmetry33, &
    lattice_SchmidMatrix_slip, &
    lattice_SchmidMatrix_twin, &
    lattice_SchmidMatrix_trans, &
    lattice_SchmidMatrix_cleavage, &
    lattice_nonSchmidMatrix, &
    lattice_interaction_SlipBySlip, &
    lattice_interaction_TwinByTwin, &
    lattice_interaction_TransByTrans, &
    lattice_interaction_SlipByTwin, &
    lattice_interaction_SlipByTrans, &
    lattice_interaction_TwinBySlip, &
    lattice_characteristicShear_Twin, &
    lattice_C66_twin, &
    lattice_C66_trans, &
    lattice_forestProjection_edge, &
    lattice_forestProjection_screw, &
    lattice_slip_normal, &
    lattice_slip_direction, &
    lattice_slip_transverse, &
    lattice_labels_slip, &
    lattice_labels_twin

contains

!--------------------------------------------------------------------------------------------------
!> @brief Module initialization
!--------------------------------------------------------------------------------------------------
subroutine lattice_init

  integer :: Nphases, p,i
  class(tNode), pointer :: &
    phases, &
    phase, &
    elasticity
 
  print'(/,a)', ' <<<+-  lattice init  -+>>>'; flush(OUTPUT_UNIT)

  phases => config_material%get('phase')
  Nphases = phases%length

  allocate(lattice_structure(Nphases),source = lattice_UNDEFINED_ID)
  allocate(lattice_C66(6,6,Nphases),  source=0.0_pReal)

  allocate(lattice_thermalConductivity  (3,3,Nphases), source=0.0_pReal)
  allocate(lattice_damageDiffusion      (3,3,Nphases), source=0.0_pReal)

  allocate(lattice_damageMobility,&
           lattice_massDensity,lattice_specificHeat, &
           lattice_mu, lattice_nu,&
           source=[(0.0_pReal,i=1,Nphases)])

  do p = 1, phases%length
    phase => phases%get(p)
    elasticity => phase%get('elasticity')
    lattice_C66(1,1,p) = elasticity%get_asFloat('C_11')
    lattice_C66(1,2,p) = elasticity%get_asFloat('C_12')

    lattice_C66(1,3,p) = elasticity%get_asFloat('C_13',defaultVal=0.0_pReal)
    lattice_C66(2,2,p) = elasticity%get_asFloat('C_22',defaultVal=0.0_pReal)
    lattice_C66(2,3,p) = elasticity%get_asFloat('C_23',defaultVal=0.0_pReal)
    lattice_C66(3,3,p) = elasticity%get_asFloat('C_33',defaultVal=0.0_pReal)
    lattice_C66(4,4,p) = elasticity%get_asFloat('C_44',defaultVal=0.0_pReal)
    lattice_C66(5,5,p) = elasticity%get_asFloat('C_55',defaultVal=0.0_pReal)
    lattice_C66(6,6,p) = elasticity%get_asFloat('C_66',defaultVal=0.0_pReal)

    select case(phase%get_asString('lattice'))
      case('iso')
        lattice_structure(p) = lattice_ISO_ID
      case('fcc')
        lattice_structure(p) = lattice_FCC_ID
      case('bcc')
        lattice_structure(p) = lattice_BCC_ID
      case('hex')
        lattice_structure(p) = lattice_HEX_ID
      case('bct')
        lattice_structure(p) = lattice_BCT_ID
      case('ort')
        lattice_structure(p) = lattice_ORT_ID
      case default
        call IO_error(130,ext_msg='lattice_init: '//phase%get_asString('lattice'))
    end select

    lattice_C66(1:6,1:6,p) = applyLatticeSymmetryC66(lattice_C66(1:6,1:6,p),phase%get_asString('lattice'))

    lattice_mu(p) = equivalent_mu(lattice_C66(1:6,1:6,p),'voigt')
    lattice_nu(p) = equivalent_nu(lattice_C66(1:6,1:6,p),'voigt')

    lattice_C66(1:6,1:6,p) = math_sym3333to66(math_Voigt66to3333(lattice_C66(1:6,1:6,p)))           ! Literature data is in Voigt notation
    do i = 1, 6
      if (abs(lattice_C66(i,i,p))<tol_math_check) &
        call IO_error(135,el=i,ip=p,ext_msg='matrix diagonal "el"ement of phase "ip"')
    enddo


    ! SHOULD NOT BE PART OF LATTICE BEGIN
    lattice_thermalConductivity(1,1,p) = phase%get_asFloat('K_11',defaultVal=0.0_pReal)
    lattice_thermalConductivity(2,2,p) = phase%get_asFloat('K_22',defaultVal=0.0_pReal)
    lattice_thermalConductivity(3,3,p) = phase%get_asFloat('K_33',defaultVal=0.0_pReal)
    lattice_thermalConductivity(1:3,1:3,p) = lattice_applyLatticeSymmetry33(lattice_thermalConductivity(1:3,1:3,p), &
                                                                                 phase%get_asString('lattice'))

    lattice_specificHeat(p) = phase%get_asFloat('c_p',defaultVal=0.0_pReal)
    lattice_massDensity(p)  = phase%get_asFloat('rho', defaultVal=0.0_pReal)

    lattice_DamageDiffusion(1,1,p) = phase%get_asFloat('D_11',defaultVal=0.0_pReal)
    lattice_DamageDiffusion(2,2,p) = phase%get_asFloat('D_22',defaultVal=0.0_pReal)
    lattice_DamageDiffusion(3,3,p) = phase%get_asFloat('D_33',defaultVal=0.0_pReal)
    lattice_DamageDiffusion(1:3,1:3,p) = lattice_applyLatticeSymmetry33(lattice_DamageDiffusion(1:3,1:3,p), &
                                                                     phase%get_asString('lattice'))

    lattice_DamageMobility(p) = phase%get_asFloat('M',defaultVal=0.0_pReal)
    ! SHOULD NOT BE PART OF LATTICE END

    call selfTest

  enddo

end subroutine lattice_init


!--------------------------------------------------------------------------------------------------
!> @brief Characteristic shear for twinning
!--------------------------------------------------------------------------------------------------
function lattice_characteristicShear_Twin(Ntwin,structure,CoverA) result(characteristicShear)

  integer,     dimension(:),            intent(in) :: Ntwin                                         !< number of active twin systems per family
  character(len=*),                     intent(in) :: structure                                     !< lattice structure
  real(pReal),                          intent(in) :: cOverA                                        !< c/a ratio
  real(pReal), dimension(sum(Ntwin))               :: characteristicShear

  integer :: &
    a, &                                                                                            !< index of active system
    p, &                                                                                            !< index in potential system list
    f, &                                                                                            !< index of my family
    s                                                                                               !< index of my system in current family

  integer, dimension(HEX_NTWIN), parameter :: &
    HEX_SHEARTWIN = reshape( [&
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
      ],[HEX_NTWIN])                                                                                ! indicator to formulas below

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_characteristicShear_Twin: '//trim(structure))

  a = 0
  myFamilies: do f = 1,size(Ntwin,1)
    mySystems: do s = 1,Ntwin(f)
      a = a + 1
      select case(structure)
        case('fcc','bcc')
          characteristicShear(a) = 0.5_pReal*sqrt(2.0_pReal)
        case('hex')
          if (cOverA < 1.0_pReal .or. cOverA > 2.0_pReal) &
            call IO_error(131,ext_msg='lattice_characteristicShear_Twin')
          p = sum(HEX_NTWINSYSTEM(1:f-1))+s
          select case(HEX_SHEARTWIN(p))                                                             ! from Christian & Mahajan 1995 p.29
            case (1)                                                                                ! <-10.1>{10.2}
              characteristicShear(a) = (3.0_pReal-cOverA**2.0_pReal)/sqrt(3.0_pReal)/CoverA
            case (2)                                                                                ! <11.6>{-1-1.1}
              characteristicShear(a) = 1.0_pReal/cOverA
            case (3)                                                                                ! <10.-2>{10.1}
              characteristicShear(a) = (4.0_pReal*cOverA**2.0_pReal-9.0_pReal)/sqrt(48.0_pReal)/cOverA
            case (4)                                                                                ! <11.-3>{11.2}
              characteristicShear(a) = 2.0_pReal*(cOverA**2.0_pReal-2.0_pReal)/3.0_pReal/cOverA
          end select
        case default
          call IO_error(137,ext_msg='lattice_characteristicShear_Twin: '//trim(structure))
      end select
    enddo mySystems
  enddo myFamilies

end function lattice_characteristicShear_Twin


!--------------------------------------------------------------------------------------------------
!> @brief Rotated elasticity matrices for twinning in 66-vector notation
!--------------------------------------------------------------------------------------------------
function lattice_C66_twin(Ntwin,C66,structure,CoverA)

  integer,     dimension(:),            intent(in) :: Ntwin                                         !< number of active twin systems per family
  character(len=*),                     intent(in) :: structure                                     !< lattice structure
  real(pReal), dimension(6,6),          intent(in) :: C66                                           !< unrotated parent stiffness matrix
  real(pReal),                          intent(in) :: cOverA                                        !< c/a ratio
  real(pReal), dimension(6,6,sum(Ntwin))           :: lattice_C66_twin

  real(pReal), dimension(3,3,sum(Ntwin)):: coordinateSystem
  type(rotation)                        :: R
  integer                               :: i

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_C66_twin: '//trim(structure))

  select case(structure)
    case('fcc')
      coordinateSystem = buildCoordinateSystem(Ntwin,FCC_NSLIPSYSTEM,FCC_SYSTEMTWIN,&
                                               trim(structure),0.0_pReal)
    case('bcc')
      coordinateSystem = buildCoordinateSystem(Ntwin,BCC_NSLIPSYSTEM,BCC_SYSTEMTWIN,&
                                               trim(structure),0.0_pReal)
    case('hex')
      coordinateSystem = buildCoordinateSystem(Ntwin,HEX_NSLIPSYSTEM,HEX_SYSTEMTWIN,&
                                               'hex',cOverA)
    case default
      call IO_error(137,ext_msg='lattice_C66_twin: '//trim(structure))
  end select

  do i = 1, sum(Ntwin)
    call R%fromAxisAngle([coordinateSystem(1:3,2,i),PI],P=1)                                        ! ToDo: Why always 180 deg?
    lattice_C66_twin(1:6,1:6,i) = R%rotTensor4sym(C66)
  enddo

end function lattice_C66_twin


!--------------------------------------------------------------------------------------------------
!> @brief Rotated elasticity matrices for transformation in 66-vector notation
!--------------------------------------------------------------------------------------------------
function lattice_C66_trans(Ntrans,C_parent66,structure_target, &
                           cOverA_trans,a_bcc,a_fcc)

  integer,     dimension(:),             intent(in) :: Ntrans                                       !< number of active twin systems per family
  character(len=*),                      intent(in) :: structure_target                             !< lattice structure
  real(pReal), dimension(6,6),           intent(in) :: C_parent66
  real(pReal), dimension(6,6,sum(Ntrans))           :: lattice_C66_trans

  real(pReal), dimension(6,6)             :: C_bar66, C_target_unrotated66
  real(pReal), dimension(3,3,sum(Ntrans)) :: Q,S
  type(rotation)                          :: R
  real(pReal)                             :: a_bcc, a_fcc, cOverA_trans
  integer                                 :: i

  if (len_trim(structure_target) /= 3) &
    call IO_error(137,ext_msg='lattice_C66_trans (target): '//trim(structure_target))

 !--------------------------------------------------------------------------------------------------
 ! elasticity matrix of the target phase in cube orientation
  if (structure_target(1:3) == 'hex') then
    if (cOverA_trans < 1.0_pReal .or. cOverA_trans > 2.0_pReal) &
      call IO_error(131,ext_msg='lattice_C66_trans: '//trim(structure_target))
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
    C_target_unrotated66(4,4) = C_bar66(4,4) - C_bar66(1,4)**2.0_pReal/(0.5_pReal*(C_bar66(1,1) - C_bar66(1,2)))
    C_target_unrotated66 = applyLatticeSymmetryC66(C_target_unrotated66,'hex')
  elseif (structure_target(1:3)  == 'bcc') then
    if (a_bcc <= 0.0_pReal .or. a_fcc <= 0.0_pReal) &
      call IO_error(134,ext_msg='lattice_C66_trans: '//trim(structure_target))
    C_target_unrotated66 = C_parent66
  else
    call IO_error(137,ext_msg='lattice_C66_trans : '//trim(structure_target))
  endif

  do i = 1, 6
    if (abs(C_target_unrotated66(i,i))<tol_math_check) &
    call IO_error(135,el=i,ext_msg='matrix diagonal "el"ement in transformation')
  enddo

  call buildTransformationSystem(Q,S,Ntrans,cOverA_trans,a_fcc,a_bcc)

  do i = 1, sum(Ntrans)
    call R%fromMatrix(Q(1:3,1:3,i))
    lattice_C66_trans(1:6,1:6,i) = R%rotTensor4sym(C_target_unrotated66)
  enddo

 end function lattice_C66_trans


!--------------------------------------------------------------------------------------------------
!> @brief Non-schmid projections for bcc with up to 6 coefficients
! Koester et al. 2012, Acta Materialia 60 (2012) 3894–3901, eq. (17)
! Gröger et al. 2008, Acta Materialia 56 (2008) 5412–5425, table 1
!--------------------------------------------------------------------------------------------------
function lattice_nonSchmidMatrix(Nslip,nonSchmidCoefficients,sense) result(nonSchmidMatrix)

  integer,     dimension(:),                intent(in) :: Nslip                                     !< number of active slip systems per family
  real(pReal), dimension(:),                intent(in) :: nonSchmidCoefficients                     !< non-Schmid coefficients for projections
  integer,                                  intent(in) :: sense                                     !< sense (-1,+1)
  real(pReal), dimension(1:3,1:3,sum(Nslip))           :: nonSchmidMatrix

  real(pReal), dimension(1:3,1:3,sum(Nslip))           :: coordinateSystem                          !< coordinate system of slip system
  real(pReal), dimension(3)                            :: direction, normal, np
  type(rotation)                                       :: R
  integer                                              :: i

  if (abs(sense) /= 1) error stop 'Sense in lattice_nonSchmidMatrix'

  coordinateSystem  = buildCoordinateSystem(Nslip,BCC_NSLIPSYSTEM,BCC_SYSTEMSLIP,&
                                            'bcc',0.0_pReal)
  coordinateSystem(1:3,1,1:sum(Nslip)) = coordinateSystem(1:3,1,1:sum(Nslip))*real(sense,pReal)     ! convert unidirectional coordinate system
  nonSchmidMatrix = lattice_SchmidMatrix_slip(Nslip,'bcc',0.0_pReal)                                ! Schmid contribution

  do i = 1,sum(Nslip)
    direction = coordinateSystem(1:3,1,i)
    normal    = coordinateSystem(1:3,2,i)
    call R%fromAxisAngle([direction,60.0_pReal],degrees=.true.,P=1)
    np = R%rotate(normal)

    if (size(nonSchmidCoefficients)>0) nonSchmidMatrix(1:3,1:3,i) = nonSchmidMatrix(1:3,1:3,i) &
      + nonSchmidCoefficients(1) * math_outer(direction, np)
    if (size(nonSchmidCoefficients)>1) nonSchmidMatrix(1:3,1:3,i) = nonSchmidMatrix(1:3,1:3,i) &
      + nonSchmidCoefficients(2) * math_outer(math_cross(normal, direction), normal)
    if (size(nonSchmidCoefficients)>2) nonSchmidMatrix(1:3,1:3,i) = nonSchmidMatrix(1:3,1:3,i) &
      + nonSchmidCoefficients(3) * math_outer(math_cross(np, direction), np)
    if (size(nonSchmidCoefficients)>3) nonSchmidMatrix(1:3,1:3,i) = nonSchmidMatrix(1:3,1:3,i) &
      + nonSchmidCoefficients(4) * math_outer(normal, normal)
    if (size(nonSchmidCoefficients)>4) nonSchmidMatrix(1:3,1:3,i) = nonSchmidMatrix(1:3,1:3,i) &
      + nonSchmidCoefficients(5) * math_outer(math_cross(normal, direction), &
                                              math_cross(normal, direction))
    if (size(nonSchmidCoefficients)>5) nonSchmidMatrix(1:3,1:3,i) = nonSchmidMatrix(1:3,1:3,i) &
      + nonSchmidCoefficients(6) * math_outer(direction, direction)
  enddo

end function lattice_nonSchmidMatrix


!--------------------------------------------------------------------------------------------------
!> @brief Slip-slip interaction matrix
!> details only active slip systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_interaction_SlipBySlip(Nslip,interactionValues,structure) result(interactionMatrix)

  integer,         dimension(:),                   intent(in) :: Nslip                              !< number of active slip systems per family
  real(pReal),     dimension(:),                   intent(in) :: interactionValues                  !< values for slip-slip interaction
  character(len=*),                                intent(in) :: structure                          !< lattice structure
  real(pReal),     dimension(sum(Nslip),sum(Nslip))           :: interactionMatrix

  integer, dimension(:),   allocatable :: NslipMax
  integer, dimension(:,:), allocatable :: interactionTypes

  integer, dimension(FCC_NSLIP,FCC_NSLIP), parameter :: &
    FCC_INTERACTIONSLIPSLIP = reshape( [&
       1, 2, 2, 4, 6, 5, 3, 5, 5, 4, 5, 6,  9,10, 9,10,11,12, & ! -----> acting
       2, 1, 2, 6, 4, 5, 5, 4, 6, 5, 3, 5,  9,10,11,12, 9,10, & ! |
       2, 2, 1, 5, 5, 3, 5, 6, 4, 6, 5, 4, 11,12, 9,10, 9,10, & ! |
       4, 6, 5, 1, 2, 2, 4, 5, 6, 3, 5, 5,  9,10,10, 9,12,11, & ! v
       6, 4, 5, 2, 1, 2, 5, 3, 5, 5, 4, 6,  9,10,12,11,10, 9, & ! reacting
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
      ],shape(FCC_INTERACTIONSLIPSLIP))                                                             !< Slip--slip interaction types for fcc
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

  integer, dimension(BCC_NSLIP,BCC_NSLIP), parameter :: &
    BCC_INTERACTIONSLIPSLIP = reshape( [&
      1,2,6,6,5,4,4,3,4,3,5,4, 6,6,4,3,3,4,6,6,4,3,6,6, & ! -----> acting
      2,1,6,6,4,3,5,4,5,4,4,3, 6,6,3,4,4,3,6,6,3,4,6,6, & ! |
      6,6,1,2,4,5,3,4,4,5,3,4, 4,3,6,6,6,6,3,4,6,6,4,3, & ! |
      6,6,2,1,3,4,4,5,3,4,4,5, 3,4,6,6,6,6,4,3,6,6,3,4, & ! v
      5,4,4,3,1,2,6,6,3,4,5,4, 3,6,4,6,6,4,6,3,4,6,3,6, & ! reacting
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
      ],shape(BCC_INTERACTIONSLIPSLIP))                                                             !< Slip--slip interaction types for bcc from Queyreau et al. Int J Plast 25 (2009) 361–377
                                                                                                    !< 1: self interaction
                                                                                                    !< 2: coplanar interaction
                                                                                                    !< 3: collinear interaction
                                                                                                    !< 4: mixed-asymmetrical junction
                                                                                                    !< 5: mixed-symmetrical junction
                                                                                                    !< 6: edge junction

  integer, dimension(HEX_NSLIP,HEX_NSLIP), parameter :: &
    HEX_INTERACTIONSLIPSLIP = reshape( [&
       1, 2, 2,   3, 3, 3,   7, 7, 7,  13,13,13,13,13,13,  21,21,21,21,21,21,21,21,21,21,21,21,  31,31,31,31,31,31, & ! -----> acting
       2, 1, 2,   3, 3, 3,   7, 7, 7,  13,13,13,13,13,13,  21,21,21,21,21,21,21,21,21,21,21,21,  31,31,31,31,31,31, & ! |
       2, 2, 1,   3, 3, 3,   7, 7, 7,  13,13,13,13,13,13,  21,21,21,21,21,21,21,21,21,21,21,21,  31,31,31,31,31,31, & ! |
     !                                                                                                                ! v
       6, 6, 6,   4, 5, 5,   8, 8, 8,  14,14,14,14,14,14,  22,22,22,22,22,22,22,22,22,22,22,22,  32,32,32,32,32,32, & ! reacting
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
     ],shape(HEX_INTERACTIONSLIPSLIP))                                                              !< Slip--slip interaction types for hex (onion peel naming scheme)

  integer, dimension(BCT_NSLIP,BCT_NSLIP), parameter :: &
    BCT_INTERACTIONSLIPSLIP = reshape( [&
      1,  2,   3,  3,   7,  7,  13, 13, 13, 13,  21, 21,  31, 31, 31, 31,  43, 43,  57, 57,  73, 73, 73, 73,  91, 91, 91, 91, 91, 91, 91, 91,  111, 111, 111, 111, 133,133,133,133,133,133,133,133, 157,157,157,157,157,157,157,157, & ! -----> acting
      2,  1,   3,  3,   7,  7,  13, 13, 13, 13,  21, 21,  31, 31, 31, 31,  43, 43,  57, 57,  73, 73, 73, 73,  91, 91, 91, 91, 91, 91, 91, 91,  111, 111, 111, 111, 133,133,133,133,133,133,133,133, 157,157,157,157,157,157,157,157, & ! |
     !                                                                                                                                                                                                                                   |
      6,  6,   4,  5,   8,  8,  14, 14, 14, 14,  22, 22,  32, 32, 32, 32,  44, 44,  58, 58,  74, 74, 74, 74,  92, 92, 92, 92, 92, 92, 92, 92,  112, 112, 112, 112, 134,134,134,134,134,134,134,134, 158,158,158,158,158,158,158,158, & ! v
      6,  6,   5,  4,   8,  8,  14, 14, 14, 14,  22, 22,  32, 32, 32, 32,  44, 44,  58, 58,  74, 74, 74, 74,  92, 92, 92, 92, 92, 92, 92, 92,  112, 112, 112, 112, 134,134,134,134,134,134,134,134, 158,158,158,158,158,158,158,158, & ! reacting
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
  ],shape(BCT_INTERACTIONSLIPSLIP))


  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_interaction_SlipBySlip: '//trim(structure))

  select case(structure)
    case('fcc')
      interactionTypes = FCC_INTERACTIONSLIPSLIP
      NslipMax         = FCC_NSLIPSYSTEM
    case('bcc')
      interactionTypes = BCC_INTERACTIONSLIPSLIP
      NslipMax         = BCC_NSLIPSYSTEM
    case('hex')
      interactionTypes = HEX_INTERACTIONSLIPSLIP
      NslipMax         = HEX_NSLIPSYSTEM
    case('bct')
      interactionTypes = BCT_INTERACTIONSLIPSLIP
      NslipMax         = BCT_NSLIPSYSTEM
    case default
      call IO_error(137,ext_msg='lattice_interaction_SlipBySlip: '//trim(structure))
  end select

  interactionMatrix = buildInteraction(Nslip,Nslip,NslipMax,NslipMax,interactionValues,interactionTypes)

end function lattice_interaction_SlipBySlip


!--------------------------------------------------------------------------------------------------
!> @brief Twin-twin interaction matrix
!> details only active twin systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_interaction_TwinByTwin(Ntwin,interactionValues,structure) result(interactionMatrix)

  integer,         dimension(:),                   intent(in) :: Ntwin                              !< number of active twin systems per family
  real(pReal),     dimension(:),                   intent(in) :: interactionValues                  !< values for twin-twin interaction
  character(len=*),                                intent(in) :: structure                          !< lattice structure
  real(pReal),     dimension(sum(Ntwin),sum(Ntwin))           :: interactionMatrix

  integer, dimension(:),   allocatable :: NtwinMax
  integer, dimension(:,:), allocatable :: interactionTypes

  integer, dimension(FCC_NTWIN,FCC_NTWIN), parameter :: &
    FCC_INTERACTIONTWINTWIN = reshape( [&
      1,1,1,2,2,2,2,2,2,2,2,2, & ! -----> acting
      1,1,1,2,2,2,2,2,2,2,2,2, & ! |
      1,1,1,2,2,2,2,2,2,2,2,2, & ! |
      2,2,2,1,1,1,2,2,2,2,2,2, & ! v
      2,2,2,1,1,1,2,2,2,2,2,2, & ! reacting
      2,2,2,1,1,1,2,2,2,2,2,2, &
      2,2,2,2,2,2,1,1,1,2,2,2, &
      2,2,2,2,2,2,1,1,1,2,2,2, &
      2,2,2,2,2,2,1,1,1,2,2,2, &
      2,2,2,2,2,2,2,2,2,1,1,1, &
      2,2,2,2,2,2,2,2,2,1,1,1, &
      2,2,2,2,2,2,2,2,2,1,1,1  &
      ],shape(FCC_INTERACTIONTWINTWIN))                                                             !< Twin-twin interaction types for fcc

  integer, dimension(BCC_NTWIN,BCC_NTWIN), parameter :: &
    BCC_INTERACTIONTWINTWIN = reshape( [&
      1,3,3,3,3,3,3,2,3,3,2,3, & ! -----> acting
      3,1,3,3,3,3,2,3,3,3,3,2, & ! |
      3,3,1,3,3,2,3,3,2,3,3,3, & ! |
      3,3,3,1,2,3,3,3,3,2,3,3, & ! v
      3,3,3,2,1,3,3,3,3,2,3,3, & ! reacting
      3,3,2,3,3,1,3,3,2,3,3,3, &
      3,2,3,3,3,3,1,3,3,3,3,2, &
      2,3,3,3,3,3,3,1,3,3,2,3, &
      3,3,2,3,3,2,3,3,1,3,3,3, &
      3,3,3,2,2,3,3,3,3,1,3,3, &
      2,3,3,3,3,3,3,2,3,3,1,3, &
      3,2,3,3,3,3,2,3,3,3,3,1  &
      ],shape(BCC_INTERACTIONTWINTWIN))                                                             !< Twin-twin interaction types for bcc
                                                                                                    !< 1: self interaction
                                                                                                    !< 2: collinear interaction
                                                                                                    !< 3: other interaction
  integer, dimension(HEX_NTWIN,HEX_NTWIN), parameter :: &
    HEX_INTERACTIONTWINTWIN = reshape( [&
       1, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, & ! -----> acting
       2, 1, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, & ! |
       2, 2, 1, 2, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, & ! |
       2, 2, 2, 1, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, & ! v
       2, 2, 2, 2, 1, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, & ! reacting
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
      ],shape(HEX_INTERACTIONTWINTWIN))                                                             !< Twin-twin interaction types for hex

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_interaction_TwinByTwin: '//trim(structure))

  select case(structure)
    case('fcc')
      interactionTypes = FCC_INTERACTIONTWINTWIN
      NtwinMax         = FCC_NTWINSYSTEM
    case('bcc')
      interactionTypes = BCC_INTERACTIONTWINTWIN
      NtwinMax         = BCC_NTWINSYSTEM
    case('hex')
      interactionTypes = HEX_INTERACTIONTWINTWIN
      NtwinMax         = HEX_NTWINSYSTEM
    case default
      call IO_error(137,ext_msg='lattice_interaction_TwinByTwin: '//trim(structure))
  end select

  interactionMatrix = buildInteraction(Ntwin,Ntwin,NtwinMax,NtwinMax,interactionValues,interactionTypes)

end function lattice_interaction_TwinByTwin


!--------------------------------------------------------------------------------------------------
!> @brief Trans-trans interaction matrix
!> details only active trans systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_interaction_TransByTrans(Ntrans,interactionValues,structure) result(interactionMatrix)

  integer,         dimension(:),                     intent(in) :: Ntrans                           !< number of active trans systems per family
  real(pReal),     dimension(:),                     intent(in) :: interactionValues                !< values for trans-trans interaction
  character(len=*),                                  intent(in) :: structure                        !< lattice structure (parent crystal)
  real(pReal),     dimension(sum(Ntrans),sum(Ntrans))           :: interactionMatrix

  integer, dimension(:),   allocatable :: NtransMax
  integer, dimension(:,:), allocatable :: interactionTypes

  integer, dimension(FCC_NTRANS,FCC_NTRANS), parameter :: &
    FCC_INTERACTIONTRANSTRANS = reshape( [&
      1,1,1,2,2,2,2,2,2,2,2,2, & ! -----> acting
      1,1,1,2,2,2,2,2,2,2,2,2, & ! |
      1,1,1,2,2,2,2,2,2,2,2,2, & ! |
      2,2,2,1,1,1,2,2,2,2,2,2, & ! v
      2,2,2,1,1,1,2,2,2,2,2,2, & ! reacting
      2,2,2,1,1,1,2,2,2,2,2,2, &
      2,2,2,2,2,2,1,1,1,2,2,2, &
      2,2,2,2,2,2,1,1,1,2,2,2, &
      2,2,2,2,2,2,1,1,1,2,2,2, &
      2,2,2,2,2,2,2,2,2,1,1,1, &
      2,2,2,2,2,2,2,2,2,1,1,1, &
      2,2,2,2,2,2,2,2,2,1,1,1  &
      ],shape(FCC_INTERACTIONTRANSTRANS))                                                           !< Trans-trans interaction types for fcc

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_interaction_TransByTrans: '//trim(structure))

  if(structure == 'fcc') then
    interactionTypes = FCC_INTERACTIONTRANSTRANS
    NtransMax        = FCC_NTRANSSYSTEM
  else
    call IO_error(137,ext_msg='lattice_interaction_TransByTrans: '//trim(structure))
  end if

  interactionMatrix = buildInteraction(Ntrans,Ntrans,NtransMax,NtransMax,interactionValues,interactionTypes)

end function lattice_interaction_TransByTrans


!--------------------------------------------------------------------------------------------------
!> @brief Slip-twin interaction matrix
!> details only active slip and twin systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_interaction_SlipByTwin(Nslip,Ntwin,interactionValues,structure) result(interactionMatrix)

  integer,         dimension(:),                   intent(in) :: Nslip, &                           !< number of active slip systems per family
                                                                 Ntwin                              !< number of active twin systems per family
  real(pReal),     dimension(:),                   intent(in) :: interactionValues                  !< values for slip-twin interaction
  character(len=*),                                intent(in) :: structure                          !< lattice structure
  real(pReal),     dimension(sum(Nslip),sum(Ntwin))           :: interactionMatrix

  integer, dimension(:),   allocatable :: NslipMax, &
                                          NtwinMax
  integer, dimension(:,:), allocatable :: interactionTypes

  integer, dimension(FCC_NTWIN,FCC_NSLIP), parameter :: &
    FCC_INTERACTIONSLIPTWIN = reshape( [&
      1,1,1,3,3,3,2,2,2,3,3,3, & ! -----> twin (acting)
      1,1,1,3,3,3,3,3,3,2,2,2, & ! |
      1,1,1,2,2,2,3,3,3,3,3,3, & ! |
      3,3,3,1,1,1,3,3,3,2,2,2, & ! v
      3,3,3,1,1,1,2,2,2,3,3,3, & ! slip (reacting)
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
      ],shape(FCC_INTERACTIONSLIPTWIN))                                                             !< Slip-twin interaction types for fcc
                                                                                                    !< 1: coplanar interaction
                                                                                                    !< 2: screw trace between slip system and twin habit plane (easy cross slip)
                                                                                                    !< 3: other interaction
  integer, dimension(BCC_NTWIN,BCC_NSLIP), parameter :: &
   BCC_INTERACTIONSLIPTWIN = reshape( [&
      3,3,3,2,2,3,3,3,3,2,3,3, & ! -----> twin (acting)
      3,3,2,3,3,2,3,3,2,3,3,3, & ! |
      3,2,3,3,3,3,2,3,3,3,3,2, & ! |
      2,3,3,3,3,3,3,2,3,3,2,3, & ! v
      2,3,3,3,3,3,3,2,3,3,2,3, & ! slip (reacting)
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
      ],shape(BCC_INTERACTIONSLIPTWIN))                                                             !< Slip-twin interaction types for bcc
                                                                                                    !< 1: coplanar interaction
                                                                                                    !< 2: screw trace between slip system and twin habit plane (easy cross slip)
                                                                                                    !< 3: other interaction
  integer, dimension(HEX_NTWIN,HEX_NSLIP), parameter :: &
    HEX_INTERACTIONSLIPTWIN = reshape( [&
       1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! ----> twin (acting)
       1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! |
       1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! |
     !                                                                                   v
       5, 5, 5, 5, 5, 5,   6, 6, 6, 6, 6, 6,   7, 7, 7, 7, 7, 7,   8, 8, 8, 8, 8, 8, & ! slip (reacting)
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
      ],shape(HEX_INTERACTIONSLIPTWIN))                                                             !< Slip-twin interaction types for hex

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_interaction_SlipByTwin: '//trim(structure))

  select case(structure)
    case('fcc')
      interactionTypes = FCC_INTERACTIONSLIPTWIN
      NslipMax         = FCC_NSLIPSYSTEM
      NtwinMax         = FCC_NTWINSYSTEM
    case('bcc')
      interactionTypes = BCC_INTERACTIONSLIPTWIN
      NslipMax         = BCC_NSLIPSYSTEM
      NtwinMax         = BCC_NTWINSYSTEM
    case('hex')
      interactionTypes = HEX_INTERACTIONSLIPTWIN
      NslipMax         = HEX_NSLIPSYSTEM
      NtwinMax         = HEX_NTWINSYSTEM
    case default
      call IO_error(137,ext_msg='lattice_interaction_SlipByTwin: '//trim(structure))
  end select

  interactionMatrix = buildInteraction(Nslip,Ntwin,NslipMax,NtwinMax,interactionValues,interactionTypes)

end function lattice_interaction_SlipByTwin


!--------------------------------------------------------------------------------------------------
!> @brief Slip-trans interaction matrix
!> details only active slip and trans systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_interaction_SlipByTrans(Nslip,Ntrans,interactionValues,structure) result(interactionMatrix)

  integer,         dimension(:),                    intent(in) :: Nslip, &                          !< number of active slip systems per family
                                                                  Ntrans                            !< number of active trans systems per family
  real(pReal),     dimension(:),                    intent(in) :: interactionValues                 !< values for slip-trans interaction
  character(len=*),                                 intent(in) :: structure                         !< lattice structure (parent crystal)
  real(pReal),     dimension(sum(Nslip),sum(Ntrans))           :: interactionMatrix

  integer, dimension(:),   allocatable :: NslipMax, &
                                          NtransMax
  integer, dimension(:,:), allocatable :: interactionTypes

  integer, dimension(FCC_NTRANS,FCC_NSLIP), parameter :: &
    FCC_INTERACTIONSLIPTRANS = reshape( [&
      1,1,1,3,3,3,2,2,2,3,3,3, & ! -----> trans (acting)
      1,1,1,3,3,3,3,3,3,2,2,2, & ! |
      1,1,1,2,2,2,3,3,3,3,3,3, & ! |
      3,3,3,1,1,1,3,3,3,2,2,2, & ! v
      3,3,3,1,1,1,2,2,2,3,3,3, & ! slip (reacting)
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
      ],shape(FCC_INTERACTIONSLIPTRANS))                                                            !< Slip-trans interaction types for fcc

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_interaction_SlipByTrans: '//trim(structure))

  select case(structure)
    case('fcc')
      interactionTypes = FCC_INTERACTIONSLIPTRANS
      NslipMax         = FCC_NSLIPSYSTEM
      NtransMax        = FCC_NTRANSSYSTEM
    case default
      call IO_error(137,ext_msg='lattice_interaction_SlipByTrans: '//trim(structure))
  end select

  interactionMatrix = buildInteraction(Nslip,Ntrans,NslipMax,NtransMax,interactionValues,interactionTypes)

 end function lattice_interaction_SlipByTrans


!--------------------------------------------------------------------------------------------------
!> @brief Twin-slip interaction matrix
!> details only active twin and slip systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_interaction_TwinBySlip(Ntwin,Nslip,interactionValues,structure) result(interactionMatrix)

  integer,         dimension(:),                   intent(in) :: Ntwin, &                           !< number of active twin systems per family
                                                                 Nslip                              !< number of active slip systems per family
  real(pReal),     dimension(:),                   intent(in) :: interactionValues                  !< values for twin-twin interaction
  character(len=*),                                intent(in) :: structure                          !< lattice structure
  real(pReal),     dimension(sum(Ntwin),sum(Nslip))           :: interactionMatrix

  integer, dimension(:),   allocatable :: NtwinMax, &
                                          NslipMax
  integer, dimension(:,:), allocatable :: interactionTypes

  integer, dimension(FCC_NSLIP,FCC_NTWIN), parameter :: &
    FCC_INTERACTIONTWINSLIP = 1                                                                     !< Twin-slip interaction types for fcc

  integer, dimension(BCC_NSLIP,BCC_NTWIN), parameter :: &
    BCC_INTERACTIONTWINSLIP = 1                                                                     !< Twin-slip interaction types for bcc

  integer, dimension(HEX_NSLIP,HEX_NTWIN), parameter :: &
    HEX_INTERACTIONTWINSLIP = reshape( [&
       1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! ----> slip (acting)
       1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! |
       1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! |
       1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! v
       1, 1, 1,   5, 5, 5,   9, 9, 9,  13,13,13,13,13,13,  17,17,17,17,17,17,17,17,17,17,17,17,  21,21,21,21,21,21, & ! twin (reacting)
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
      ],shape(HEX_INTERACTIONTWINSLIP))                                                             !< Twin-slip interaction types for hex

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_interaction_TwinBySlip: '//trim(structure))

  select case(structure)
    case('fcc')
      interactionTypes = FCC_INTERACTIONTWINSLIP
      NtwinMax         = FCC_NTWINSYSTEM
      NslipMax         = FCC_NSLIPSYSTEM
    case('bcc')
      interactionTypes = BCC_INTERACTIONTWINSLIP
      NtwinMax         = BCC_NTWINSYSTEM
      NslipMax         = BCC_NSLIPSYSTEM
    case('hex')
      interactionTypes = HEX_INTERACTIONTWINSLIP
      NtwinMax         = HEX_NTWINSYSTEM
      NslipMax         = HEX_NSLIPSYSTEM
    case default
      call IO_error(137,ext_msg='lattice_interaction_TwinBySlip: '//trim(structure))
  end select

  interactionMatrix = buildInteraction(Ntwin,Nslip,NtwinMax,NslipMax,interactionValues,interactionTypes)

end function lattice_interaction_TwinBySlip


!--------------------------------------------------------------------------------------------------
!> @brief Schmid matrix for slip
!> details only active slip systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_SchmidMatrix_slip(Nslip,structure,cOverA) result(SchmidMatrix)

  integer,         dimension(:),            intent(in) :: Nslip                                     !< number of active slip systems per family
  character(len=*),                         intent(in) :: structure                                 !< lattice structure
  real(pReal),                              intent(in) :: cOverA
  real(pReal),     dimension(3,3,sum(Nslip))           :: SchmidMatrix

  real(pReal), dimension(3,3,sum(Nslip))             :: coordinateSystem
  real(pReal), dimension(:,:),           allocatable :: slipSystems
  integer,     dimension(:),             allocatable :: NslipMax
  integer                                            :: i

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_SchmidMatrix_slip: '//trim(structure))

  select case(structure)
    case('fcc')
      NslipMax    = FCC_NSLIPSYSTEM
      slipSystems = FCC_SYSTEMSLIP
    case('bcc')
      NslipMax    = BCC_NSLIPSYSTEM
      slipSystems = BCC_SYSTEMSLIP
    case('hex')
      NslipMax    = HEX_NSLIPSYSTEM
      slipSystems = HEX_SYSTEMSLIP
    case('bct')
      NslipMax    = BCT_NSLIPSYSTEM
      slipSystems = BCT_SYSTEMSLIP
    case default
      allocate(NslipMax(0))
      call IO_error(137,ext_msg='lattice_SchmidMatrix_slip: '//trim(structure))
  end select

  if (any(NslipMax(1:size(Nslip)) - Nslip < 0)) &
    call IO_error(145,ext_msg='Nslip '//trim(structure))
  if (any(Nslip < 0)) &
    call IO_error(144,ext_msg='Nslip '//trim(structure))

  coordinateSystem = buildCoordinateSystem(Nslip,NslipMax,slipSystems,structure,cOverA)

  do i = 1, sum(Nslip)
    SchmidMatrix(1:3,1:3,i) = math_outer(coordinateSystem(1:3,1,i),coordinateSystem(1:3,2,i))
    if (abs(math_trace33(SchmidMatrix(1:3,1:3,i))) > tol_math_check) &
      call IO_error(0,i,ext_msg = 'dilatational Schmid matrix for slip')
  enddo

end function lattice_SchmidMatrix_slip


!--------------------------------------------------------------------------------------------------
!> @brief Schmid matrix for twinning
!> details only active twin systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_SchmidMatrix_twin(Ntwin,structure,cOverA) result(SchmidMatrix)

  integer,         dimension(:),            intent(in) :: Ntwin                                     !< number of active twin systems per family
  character(len=*),                         intent(in) :: structure                                 !< lattice structure
  real(pReal),                              intent(in) :: cOverA                                    !< c/a ratio
  real(pReal),     dimension(3,3,sum(Ntwin))           :: SchmidMatrix

  real(pReal), dimension(3,3,sum(Ntwin))             :: coordinateSystem
  real(pReal), dimension(:,:),           allocatable :: twinSystems
  integer,     dimension(:),             allocatable :: NtwinMax
  integer                                            :: i

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_SchmidMatrix_twin: '//trim(structure))

  select case(structure)
    case('fcc')
      NtwinMax    = FCC_NTWINSYSTEM
      twinSystems = FCC_SYSTEMTWIN
    case('bcc')
      NtwinMax    = BCC_NTWINSYSTEM
      twinSystems = BCC_SYSTEMTWIN
    case('hex')
      NtwinMax    = HEX_NTWINSYSTEM
      twinSystems = HEX_SYSTEMTWIN
    case default
      allocate(NtwinMax(0))
      call IO_error(137,ext_msg='lattice_SchmidMatrix_twin: '//trim(structure))
  end select

  if (any(NtwinMax(1:size(Ntwin)) - Ntwin < 0)) &
    call IO_error(145,ext_msg='Ntwin '//trim(structure))
  if (any(Ntwin < 0)) &
    call IO_error(144,ext_msg='Ntwin '//trim(structure))

  coordinateSystem = buildCoordinateSystem(Ntwin,NtwinMax,twinSystems,structure,cOverA)

  do i = 1, sum(Ntwin)
    SchmidMatrix(1:3,1:3,i) = math_outer(coordinateSystem(1:3,1,i),coordinateSystem(1:3,2,i))
    if (abs(math_trace33(SchmidMatrix(1:3,1:3,i))) > tol_math_check) &
      call IO_error(0,i,ext_msg = 'dilatational Schmid matrix for twin')
  enddo

end function lattice_SchmidMatrix_twin


!--------------------------------------------------------------------------------------------------
!> @brief Schmid matrix for twinning
!> details only active twin systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_SchmidMatrix_trans(Ntrans,structure_target,cOverA,a_bcc,a_fcc) result(SchmidMatrix)

  integer,         dimension(:),             intent(in) :: Ntrans                                   !< number of active twin systems per family
  character(len=*),                          intent(in) :: structure_target                         !< lattice structure
  real(pReal),                               intent(in) :: cOverA                                   !< c/a ratio
  real(pReal),     dimension(3,3,sum(Ntrans))           :: SchmidMatrix

  real(pReal), dimension(3,3,sum(Ntrans)) :: devNull
  real(pReal)                             :: a_bcc, a_fcc

  if (len_trim(structure_target) /= 3) &
    call IO_error(137,ext_msg='lattice_SchmidMatrix_trans: '//trim(structure_target))
  if (structure_target(1:3) /= 'bcc' .and. structure_target(1:3) /= 'hex') &
    call IO_error(137,ext_msg='lattice_SchmidMatrix_trans: '//trim(structure_target))

  if (structure_target(1:3) == 'hex' .and. (cOverA < 1.0_pReal .or. cOverA > 2.0_pReal)) &
    call IO_error(131,ext_msg='lattice_SchmidMatrix_trans: '//trim(structure_target))

  if (structure_target(1:3) == 'bcc' .and. (a_bcc <= 0.0_pReal .or. a_fcc <= 0.0_pReal)) &
    call IO_error(134,ext_msg='lattice_SchmidMatrix_trans: '//trim(structure_target))

  call buildTransformationSystem(devNull,SchmidMatrix,Ntrans,cOverA,a_fcc,a_bcc)

end function lattice_SchmidMatrix_trans


!--------------------------------------------------------------------------------------------------
!> @brief Schmid matrix for cleavage
!> details only active cleavage systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_SchmidMatrix_cleavage(Ncleavage,structure,cOverA) result(SchmidMatrix)

  integer,         dimension(:),                  intent(in) :: Ncleavage                           !< number of active cleavage systems per family
  character(len=*),                               intent(in) :: structure                           !< lattice structure
  real(pReal),                                    intent(in) :: cOverA                              !< c/a ratio
  real(pReal),     dimension(3,3,3,sum(Ncleavage))           :: SchmidMatrix

  real(pReal), dimension(3,3,sum(Ncleavage))             :: coordinateSystem
  real(pReal), dimension(:,:),               allocatable :: cleavageSystems
  integer,     dimension(:),                 allocatable :: NcleavageMax
  integer                                                :: i

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_SchmidMatrix_cleavage: '//trim(structure))

  select case(structure)
    case('ort')
      NcleavageMax    = ORT_NCLEAVAGESYSTEM
      cleavageSystems = ORT_SYSTEMCLEAVAGE
    case('fcc')
      NcleavageMax    = FCC_NCLEAVAGESYSTEM
      cleavageSystems = FCC_SYSTEMCLEAVAGE
    case('bcc')
      NcleavageMax    = BCC_NCLEAVAGESYSTEM
      cleavageSystems = BCC_SYSTEMCLEAVAGE
    case default
      allocate(NcleavageMax(0))
      call IO_error(137,ext_msg='lattice_SchmidMatrix_cleavage: '//trim(structure))
  end select

  if (any(NcleavageMax(1:size(Ncleavage)) - Ncleavage < 0)) &
    call IO_error(145,ext_msg='Ncleavage '//trim(structure))
  if (any(Ncleavage < 0)) &
    call IO_error(144,ext_msg='Ncleavage '//trim(structure))

  coordinateSystem = buildCoordinateSystem(Ncleavage,NcleavageMax,cleavageSystems,structure,cOverA)

  do i = 1, sum(Ncleavage)
    SchmidMatrix(1:3,1:3,1,i) = math_outer(coordinateSystem(1:3,1,i),coordinateSystem(1:3,2,i))
    SchmidMatrix(1:3,1:3,2,i) = math_outer(coordinateSystem(1:3,3,i),coordinateSystem(1:3,2,i))
    SchmidMatrix(1:3,1:3,3,i) = math_outer(coordinateSystem(1:3,2,i),coordinateSystem(1:3,2,i))
  enddo

end function lattice_SchmidMatrix_cleavage


!--------------------------------------------------------------------------------------------------
!> @brief Slip direction of slip systems (|| b)
!--------------------------------------------------------------------------------------------------
function lattice_slip_direction(Nslip,structure,cOverA) result(d)

  integer,         dimension(:),           intent(in) :: Nslip                                      !< number of active slip systems per family
  character(len=*),                        intent(in) :: structure                                  !< lattice structure
  real(pReal),                             intent(in) :: cOverA                                     !< c/a ratio
  real(pReal),     dimension(3,sum(Nslip))            :: d

  real(pReal), dimension(3,3,sum(Nslip)) :: coordinateSystem

  coordinateSystem = coordinateSystem_slip(Nslip,structure,cOverA)
  d = coordinateSystem(1:3,1,1:sum(Nslip))

end function lattice_slip_direction


!--------------------------------------------------------------------------------------------------
!> @brief Normal direction of slip systems (|| n)
!--------------------------------------------------------------------------------------------------
function lattice_slip_normal(Nslip,structure,cOverA) result(n)

  integer,         dimension(:),           intent(in) :: Nslip                                      !< number of active slip systems per family
  character(len=*),                        intent(in) :: structure                                  !< lattice structure
  real(pReal),                             intent(in) :: cOverA                                     !< c/a ratio
  real(pReal),     dimension(3,sum(Nslip))            :: n

  real(pReal), dimension(3,3,sum(Nslip)) :: coordinateSystem

  coordinateSystem = coordinateSystem_slip(Nslip,structure,cOverA)
  n = coordinateSystem(1:3,2,1:sum(Nslip))

end function lattice_slip_normal


!--------------------------------------------------------------------------------------------------
!> @brief Transverse direction of slip systems ( || t = b x n)
!--------------------------------------------------------------------------------------------------
function lattice_slip_transverse(Nslip,structure,cOverA) result(t)

  integer,         dimension(:),           intent(in) :: Nslip                                      !< number of active slip systems per family
  character(len=*),                        intent(in) :: structure                                  !< lattice structure
  real(pReal),                             intent(in) :: cOverA                                     !< c/a ratio
  real(pReal),     dimension(3,sum(Nslip))            :: t

  real(pReal), dimension(3,3,sum(Nslip)) :: coordinateSystem

  coordinateSystem = coordinateSystem_slip(Nslip,structure,cOverA)
  t = coordinateSystem(1:3,3,1:sum(Nslip))

end function lattice_slip_transverse


!--------------------------------------------------------------------------------------------------
!> @brief Labels for slip systems
!> details only active slip systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_labels_slip(Nslip,structure) result(labels)

  integer,         dimension(:),  intent(in)  :: Nslip                                              !< number of active slip systems per family
  character(len=*),               intent(in)  :: structure                                          !< lattice structure

  character(len=:), dimension(:),   allocatable :: labels

  real(pReal),      dimension(:,:), allocatable :: slipSystems
  integer,          dimension(:),   allocatable :: NslipMax

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_labels_slip: '//trim(structure))

  select case(structure)
    case('fcc')
      NslipMax    = FCC_NSLIPSYSTEM
      slipSystems = FCC_SYSTEMSLIP
    case('bcc')
      NslipMax    = BCC_NSLIPSYSTEM
      slipSystems = BCC_SYSTEMSLIP
    case('hex')
      NslipMax    = HEX_NSLIPSYSTEM
      slipSystems = HEX_SYSTEMSLIP
    case('bct')
      NslipMax    = BCT_NSLIPSYSTEM
      slipSystems = BCT_SYSTEMSLIP
    case default
      call IO_error(137,ext_msg='lattice_labels_slip: '//trim(structure))
  end select

  if (any(NslipMax(1:size(Nslip)) - Nslip < 0)) &
    call IO_error(145,ext_msg='Nslip '//trim(structure))
  if (any(Nslip < 0)) &
    call IO_error(144,ext_msg='Nslip '//trim(structure))

  labels = getLabels(Nslip,NslipMax,slipSystems)

end function lattice_labels_slip


!--------------------------------------------------------------------------------------------------
!> @brief Return 3x3 tensor with symmetry according to given crystal structure
!--------------------------------------------------------------------------------------------------
function lattice_applyLatticeSymmetry33(T,structure) result(T_sym)

  real(pReal), dimension(3,3) :: T_sym

  real(pReal), dimension(3,3), intent(in) :: T
  character(len=*),            intent(in) :: structure

  integer :: k

  T_sym = 0.0_pReal

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_applyLatticeSymmetry33: '//trim(structure))

  select case(structure)
    case('iso','fcc','bcc')
      do k=1,3
        T_sym(k,k) = T(1,1)
      enddo
    case('hex')
      T_sym(1,1) = T(1,1)
      T_sym(2,2) = T(1,1)
      T_sym(3,3) = T(3,3)
    case('ort','bct')
      T_sym(1,1) = T(1,1)
      T_sym(2,2) = T(2,2)
      T_sym(3,3) = T(3,3)
    case default
      call IO_error(137,ext_msg='lattice_applyLatticeSymmetry33: '//trim(structure))
   end select

end function lattice_applyLatticeSymmetry33


!--------------------------------------------------------------------------------------------------
!> @brief Return stiffness matrix in 6x6 notation with symmetry according to given crystal structure
!> @details J. A. Rayne and B. S. Chandrasekhar Phys. Rev. 120, 1658 Erratum Phys. Rev. 122, 1962
!--------------------------------------------------------------------------------------------------
function applyLatticeSymmetryC66(C66,structure) result(C66_sym)

  real(pReal), dimension(6,6) :: C66_sym

  real(pReal), dimension(6,6), intent(in) :: C66
  character(len=*),            intent(in) :: structure

  integer :: j,k

  C66_sym = 0.0_pReal

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='applyLatticeSymmetryC66: '//trim(structure))

  select case(structure)
    case ('iso')
      do k=1,3
        do j=1,3
          C66_sym(k,j) = C66(1,2)
        enddo
        C66_sym(k,k)     = C66(1,1)
        C66_sym(k+3,k+3) = 0.5_pReal*(C66(1,1)-C66(1,2))
      enddo
    case ('fcc','bcc')
      do k=1,3
        do j=1,3
          C66_sym(k,j) = C66(1,2)
        enddo
        C66_sym(k,k)     = C66(1,1)
        C66_sym(k+3,k+3) = C66(4,4)
      enddo
    case ('hex')
      C66_sym(1,1) = C66(1,1)
      C66_sym(2,2) = C66(1,1)
      C66_sym(3,3) = C66(3,3)
      C66_sym(1,2) = C66(1,2)
      C66_sym(2,1) = C66(1,2)
      C66_sym(1,3) = C66(1,3)
      C66_sym(3,1) = C66(1,3)
      C66_sym(2,3) = C66(1,3)
      C66_sym(3,2) = C66(1,3)
      C66_sym(4,4) = C66(4,4)
      C66_sym(5,5) = C66(4,4)
      C66_sym(6,6) = 0.5_pReal*(C66(1,1)-C66(1,2))
    case ('ort')
      C66_sym(1,1) = C66(1,1)
      C66_sym(2,2) = C66(2,2)
      C66_sym(3,3) = C66(3,3)
      C66_sym(1,2) = C66(1,2)
      C66_sym(2,1) = C66(1,2)
      C66_sym(1,3) = C66(1,3)
      C66_sym(3,1) = C66(1,3)
      C66_sym(2,3) = C66(2,3)
      C66_sym(3,2) = C66(2,3)
      C66_sym(4,4) = C66(4,4)
      C66_sym(5,5) = C66(5,5)
      C66_sym(6,6) = C66(6,6)
    case ('bct')
      C66_sym(1,1) = C66(1,1)
      C66_sym(2,2) = C66(1,1)
      C66_sym(3,3) = C66(3,3)
      C66_sym(1,2) = C66(1,2)
      C66_sym(2,1) = C66(1,2)
      C66_sym(1,3) = C66(1,3)
      C66_sym(3,1) = C66(1,3)
      C66_sym(2,3) = C66(1,3)
      C66_sym(3,2) = C66(1,3)
      C66_sym(4,4) = C66(4,4)
      C66_sym(5,5) = C66(4,4)
      C66_sym(6,6) = C66(6,6)
    case default
      call IO_error(137,ext_msg='applyLatticeSymmetryC66: '//trim(structure))
   end select

end function applyLatticeSymmetryC66


!--------------------------------------------------------------------------------------------------
!> @brief Labels for twin systems
!> details only active twin systems are considered
!--------------------------------------------------------------------------------------------------
function lattice_labels_twin(Ntwin,structure) result(labels)

  integer,         dimension(:),  intent(in)  :: Ntwin                                              !< number of active slip systems per family
  character(len=*),               intent(in)  :: structure                                          !< lattice structure

  character(len=:), dimension(:),   allocatable :: labels

  real(pReal),      dimension(:,:), allocatable :: twinSystems
  integer,          dimension(:),   allocatable :: NtwinMax

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='lattice_labels_twin: '//trim(structure))

  select case(structure)
    case('fcc')
      NtwinMax    = FCC_NTWINSYSTEM
      twinSystems = FCC_SYSTEMTWIN
    case('bcc')
      NtwinMax    = BCC_NTWINSYSTEM
      twinSystems = BCC_SYSTEMTWIN
    case('hex')
      NtwinMax    = HEX_NTWINSYSTEM
      twinSystems = HEX_SYSTEMTWIN
    case default
      call IO_error(137,ext_msg='lattice_labels_twin: '//trim(structure))
  end select

  if (any(NtwinMax(1:size(Ntwin)) - Ntwin < 0)) &
    call IO_error(145,ext_msg='Ntwin '//trim(structure))
  if (any(Ntwin < 0)) &
    call IO_error(144,ext_msg='Ntwin '//trim(structure))

  labels = getLabels(Ntwin,NtwinMax,twinSystems)

end function lattice_labels_twin


!--------------------------------------------------------------------------------------------------
!> @brief Projection of the transverse direction onto the slip plane
!> @details: This projection is used to calculate forest hardening for edge dislocations
!--------------------------------------------------------------------------------------------------
function slipProjection_transverse(Nslip,structure,cOverA) result(projection)

  integer,         dimension(:),                   intent(in) :: Nslip                              !< number of active slip systems per family
  character(len=*),                                intent(in) :: structure                          !< lattice structure
  real(pReal),                                     intent(in) :: cOverA                             !< c/a ratio
  real(pReal),     dimension(sum(Nslip),sum(Nslip))           :: projection

  real(pReal), dimension(3,sum(Nslip)) :: n, t
  integer                              :: i, j

  n = lattice_slip_normal    (Nslip,structure,cOverA)
  t = lattice_slip_transverse(Nslip,structure,cOverA)

  do i=1, sum(Nslip); do j=1, sum(Nslip)
    projection(i,j) = abs(math_inner(n(:,i),t(:,j)))
  enddo; enddo

end function slipProjection_transverse


!--------------------------------------------------------------------------------------------------
!> @brief Projection of the slip direction onto the slip plane
!> @details: This projection is used to calculate forest hardening for screw dislocations
!--------------------------------------------------------------------------------------------------
function slipProjection_direction(Nslip,structure,cOverA) result(projection)

  integer,         dimension(:),                   intent(in) :: Nslip                              !< number of active slip systems per family
  character(len=*),                                intent(in) :: structure                          !< lattice structure
  real(pReal),                                     intent(in) :: cOverA                             !< c/a ratio
  real(pReal),     dimension(sum(Nslip),sum(Nslip))           :: projection

  real(pReal), dimension(3,sum(Nslip)) :: n, d
  integer                              :: i, j

  n = lattice_slip_normal   (Nslip,structure,cOverA)
  d = lattice_slip_direction(Nslip,structure,cOverA)

  do i=1, sum(Nslip); do j=1, sum(Nslip)
    projection(i,j) = abs(math_inner(n(:,i),d(:,j)))
  enddo; enddo

end function slipProjection_direction


!--------------------------------------------------------------------------------------------------
!> @brief build a local coordinate system on slip systems
!> @details Order: Direction, plane (normal), and common perpendicular
!--------------------------------------------------------------------------------------------------
function coordinateSystem_slip(Nslip,structure,cOverA) result(coordinateSystem)

  integer,          dimension(:),            intent(in) :: Nslip                                    !< number of active slip systems per family
  character(len=*),                          intent(in) :: structure                                !< lattice structure
  real(pReal),                               intent(in) :: cOverA                                   !< c/a ratio
  real(pReal),     dimension(3,3,sum(Nslip))            :: coordinateSystem

  real(pReal), dimension(:,:), allocatable :: slipSystems
  integer,     dimension(:),   allocatable :: NslipMax

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='coordinateSystem_slip: '//trim(structure))

  select case(structure)
    case('fcc')
      NslipMax    = FCC_NSLIPSYSTEM
      slipSystems = FCC_SYSTEMSLIP
    case('bcc')
      NslipMax    = BCC_NSLIPSYSTEM
      slipSystems = BCC_SYSTEMSLIP
    case('hex')
      NslipMax    = HEX_NSLIPSYSTEM
      slipSystems = HEX_SYSTEMSLIP
    case('bct')
      NslipMax    = BCT_NSLIPSYSTEM
      slipSystems = BCT_SYSTEMSLIP
    case default
      allocate(NslipMax(0))
      call IO_error(137,ext_msg='coordinateSystem_slip: '//trim(structure))
  end select

  if (any(NslipMax(1:size(Nslip)) - Nslip < 0)) &
    call IO_error(145,ext_msg='Nslip '//trim(structure))
  if (any(Nslip < 0)) &
    call IO_error(144,ext_msg='Nslip '//trim(structure))

  coordinateSystem = buildCoordinateSystem(Nslip,NslipMax,slipSystems,structure,cOverA)

end function coordinateSystem_slip


!--------------------------------------------------------------------------------------------------
!> @brief Populate reduced interaction matrix
!--------------------------------------------------------------------------------------------------
function buildInteraction(reacting_used,acting_used,reacting_max,acting_max,values,matrix)

  integer,     dimension(:),                                 intent(in) :: &
    reacting_used, &                                                                                !< # of reacting systems per family as specified in material.config
    acting_used, &                                                                                  !< # of   acting systems per family as specified in material.config
    reacting_max, &                                                                                 !< max # of reacting systems per family for given lattice
    acting_max                                                                                      !< max # of   acting systems per family for given lattice
  real(pReal), dimension(:),                                 intent(in) :: values                   !< interaction values
  integer,     dimension(:,:),                               intent(in) :: matrix                   !< interaction types
  real(pReal), dimension(sum(reacting_used),sum(acting_used))           :: buildInteraction

  integer :: &
    acting_family_index,     acting_family,   acting_system, &
    reacting_family_index, reacting_family, reacting_system, &
    i,j,k,l

  do acting_family = 1,size(acting_used,1)
    acting_family_index = sum(acting_used(1:acting_family-1))
    do acting_system = 1,acting_used(acting_family)

      do reacting_family = 1,size(reacting_used,1)
        reacting_family_index = sum(reacting_used(1:reacting_family-1))
        do reacting_system = 1,reacting_used(reacting_family)

          i = sum(  acting_max(1:  acting_family-1)) +   acting_system
          j = sum(reacting_max(1:reacting_family-1)) + reacting_system

          k =   acting_family_index +   acting_system
          l = reacting_family_index + reacting_system

          if (matrix(i,j) > size(values)) call IO_error(138,ext_msg='buildInteraction')

          buildInteraction(l,k) = values(matrix(i,j))

      enddo; enddo
  enddo; enddo

end function buildInteraction


!--------------------------------------------------------------------------------------------------
!> @brief Build a local coordinate system on slip, twin, trans, cleavage systems
!> @details Order: Direction, plane (normal), and common perpendicular
!--------------------------------------------------------------------------------------------------
function buildCoordinateSystem(active,potential,system,structure,cOverA)

  integer, dimension(:), intent(in) :: &
    active, &                                                                                       !< # of active systems per family
    potential                                                                                       !< # of potential systems per family
  real(pReal), dimension(:,:), intent(in) :: &
    system
  character(len=*),            intent(in) :: &
    structure                                                                                       !< lattice structure
  real(pReal),                 intent(in) :: &
    cOverA
  real(pReal), dimension(3,3,sum(active)) :: &
    buildCoordinateSystem

  real(pReal), dimension(3) :: &
    direction, normal
  integer :: &
    a, &                                                                                            !< index of active system
    p, &                                                                                            !< index in potential system matrix
    f, &                                                                                            !< index of my family
    s                                                                                               !< index of my system in current family

  if (len_trim(structure) /= 3) &
    call IO_error(137,ext_msg='buildCoordinateSystem: '//trim(structure))
  if (trim(structure) == 'bct' .and. cOverA > 2.0_pReal) &
    call IO_error(131,ext_msg='buildCoordinateSystem:'//trim(structure))
  if (trim(structure) == 'hex' .and. (cOverA < 1.0_pReal .or. cOverA > 2.0_pReal)) &
    call IO_error(131,ext_msg='buildCoordinateSystem:'//trim(structure))

  a = 0
  activeFamilies: do f = 1,size(active,1)
    activeSystems: do s = 1,active(f)
      a = a + 1
      p = sum(potential(1:f-1))+s

      select case(trim(structure))

        case ('fcc','bcc','iso','ort','bct')
          direction = system(1:3,p)
          normal    = system(4:6,p)

        case ('hex')
          direction = [ system(1,p)*1.5_pReal, &
                       (system(1,p)+2.0_pReal*system(2,p))*sqrt(0.75_pReal), &
                        system(4,p)*cOverA ]                                                        ! direction [uvtw]->[3u/2 (u+2v)*sqrt(3)/2 w*(p/a)])
          normal    = [ system(5,p), &
                       (system(5,p)+2.0_pReal*system(6,p))/sqrt(3.0_pReal), &
                        system(8,p)/cOverA ]                                                        ! plane (hkil)->(h (h+2k)/sqrt(3) l/(p/a))

        case default
          call IO_error(137,ext_msg='buildCoordinateSystem: '//trim(structure))

      end select

      buildCoordinateSystem(1:3,1,a) = direction/norm2(direction)
      buildCoordinateSystem(1:3,2,a) = normal   /norm2(normal)
      buildCoordinateSystem(1:3,3,a) = math_cross(direction/norm2(direction),&
                                                  normal   /norm2(normal))

    enddo activeSystems
  enddo activeFamilies

end function buildCoordinateSystem


!--------------------------------------------------------------------------------------------------
!> @brief Helper function to define transformation systems
! Needed to calculate Schmid matrix and rotated stiffness matrices.
! @details: set c/a   = 0.0 for fcc -> bcc transformation
!           set a_Xcc = 0.0 for fcc -> hex transformation
!--------------------------------------------------------------------------------------------------
subroutine buildTransformationSystem(Q,S,Ntrans,cOverA,a_fcc,a_bcc)

  integer, dimension(:), intent(in) :: &
    Ntrans
  real(pReal),  dimension(3,3,sum(Ntrans)), intent(out) :: &
    Q, &                                                                                            !< Total rotation: Q = R*B
    S                                                                                               !< Eigendeformation tensor for phase transformation
  real(pReal),                 intent(in) :: &
    cOverA, &                                                                                       !< c/a for target hex structure
    a_bcc, &                                                                                        !< lattice parameter a for target bcc structure
    a_fcc                                                                                           !< lattice parameter a for parent fcc structure

  type(rotation) :: &
    R, &                                                                                            !< Pitsch rotation
    B                                                                                               !< Rotation of fcc to Bain coordinate system
  real(pReal), dimension(3,3) :: &
    U, &                                                                                            !< Bain deformation
    ss, sd
  real(pReal), dimension(3) :: &
    x, y, z
  integer :: &
    i
  real(pReal), dimension(3+3,FCC_NTRANS), parameter :: &
    FCCTOHEX_SYSTEMTRANS = reshape(real( [&
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
      ],pReal),shape(FCCTOHEX_SYSTEMTRANS))
       real(pReal), dimension(4,fcc_Ntrans), parameter :: &
    FCCTOBCC_SYSTEMTRANS = reshape([&
      0.0, 1.0, 0.0,     10.26, &                                                                   ! Pitsch OR (Ma & Hartmaier 2014, Table 3)
      0.0,-1.0, 0.0,     10.26, &
      0.0, 0.0, 1.0,     10.26, &
      0.0, 0.0,-1.0,     10.26, &
      1.0, 0.0, 0.0,     10.26, &
     -1.0, 0.0, 0.0,     10.26, &
      0.0, 0.0, 1.0,     10.26, &
      0.0, 0.0,-1.0,     10.26, &
      1.0, 0.0, 0.0,     10.26, &
     -1.0, 0.0, 0.0,     10.26, &
      0.0, 1.0, 0.0,     10.26, &
      0.0,-1.0, 0.0,     10.26  &
      ],shape(FCCTOBCC_SYSTEMTRANS))

  integer, dimension(9,fcc_Ntrans), parameter :: &
    FCCTOBCC_BAINVARIANT = reshape( [&
      1, 0, 0, 0, 1, 0, 0, 0, 1, &                                                                  ! Pitsch OR (Ma & Hartmaier 2014, Table 3)
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
      ],shape(FCCTOBCC_BAINVARIANT))

  real(pReal), dimension(4,fcc_Ntrans), parameter :: &
    FCCTOBCC_BAINROT = reshape([&
      1.0, 0.0, 0.0,     45.0, &                                                                    ! Rotate fcc austensite to bain variant
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
      ],shape(FCCTOBCC_BAINROT))

  if (a_bcc > 0.0_pReal .and. a_fcc > 0.0_pReal .and. dEq0(cOverA)) then                            ! fcc -> bcc transformation
    do i = 1,sum(Ntrans)
      call R%fromAxisAngle(FCCTOBCC_SYSTEMTRANS(:,i),degrees=.true.,P=1)
      call B%fromAxisAngle(FCCTOBCC_BAINROT(:,i),    degrees=.true.,P=1)
      x = real(FCCTOBCC_BAINVARIANT(1:3,i),pReal)
      y = real(FCCTOBCC_BAINVARIANT(4:6,i),pReal)
      z = real(FCCTOBCC_BAINVARIANT(7:9,i),pReal)

      U = (a_bcc/a_fcc)*math_outer(x,x) &
        + (a_bcc/a_fcc)*math_outer(y,y) * sqrt(2.0_pReal) &
        + (a_bcc/a_fcc)*math_outer(z,z) * sqrt(2.0_pReal)
      Q(1:3,1:3,i) = matmul(R%asMatrix(),B%asMatrix())
      S(1:3,1:3,i) = matmul(R%asMatrix(),U) - MATH_I3
    enddo
  elseif (cOverA > 0.0_pReal .and. dEq0(a_bcc)) then                                                ! fcc -> hex transformation
    ss      = MATH_I3
    sd      = MATH_I3
    ss(1,3) = sqrt(2.0_pReal)/4.0_pReal
    sd(3,3) = cOverA/sqrt(8.0_pReal/3.0_pReal)

    do i = 1,sum(Ntrans)
      x = FCCTOHEX_SYSTEMTRANS(1:3,i)/norm2(FCCTOHEX_SYSTEMTRANS(1:3,i))
      z = FCCTOHEX_SYSTEMTRANS(4:6,i)/norm2(FCCTOHEX_SYSTEMTRANS(4:6,i))
      y = -math_cross(x,z)
      Q(1:3,1,i) = x
      Q(1:3,2,i) = y
      Q(1:3,3,i) = z
      S(1:3,1:3,i) = matmul(Q(1:3,1:3,i), matmul(matmul(sd,ss), transpose(Q(1:3,1:3,i)))) - MATH_I3 ! ToDo: This is of interest for the Schmid matrix only
    enddo
  else
    call IO_error(132,ext_msg='buildTransformationSystem')
  endif

end subroutine buildTransformationSystem


!--------------------------------------------------------------------------------------------------
!> @brief select active systems as strings
!--------------------------------------------------------------------------------------------------
function getlabels(active,potential,system) result(labels)

  integer,          dimension(:),   intent(in) :: &
    active, &                                                                                       !< # of active systems per family
    potential                                                                                       !< # of potential systems per family
  real(pReal),      dimension(:,:), intent(in) :: &
    system

  character(len=:), dimension(:), allocatable :: labels
  character(len=:),               allocatable :: label

  integer :: i,j
  integer :: &
    a, &                                                                                            !< index of active system
    p, &                                                                                            !< index in potential system matrix
    f, &                                                                                            !< index of my family
    s                                                                                               !< index of my system in current family

  i = 2*size(system,1) + (size(system,1) - 2) + 4                                                   ! 2 letters per index + spaces + brackets
  allocate(character(len=i) :: labels(sum(active)), label)

  a = 0
  activeFamilies: do f = 1,size(active,1)
    activeSystems: do s = 1,active(f)
      a = a + 1
      p = sum(potential(1:f-1))+s

      i = 1
      label(i:i) = '['
      direction: do j = 1, size(system,1)/2
        write(label(i+1:i+2),'(I2.1)') int(system(j,p))
        label(i+3:i+3) = ' '
        i = i + 3
      enddo direction
      label(i:i) = ']'

      i = i +1
      label(i:i) = '('
      normal: do j = size(system,1)/2+1, size(system,1)
        write(label(i+1:i+2),'(I2.1)') int(system(j,p))
        label(i+3:i+3) = ' '
        i = i + 3
      enddo normal
      label(i:i) = ')'

      labels(s) = label

    enddo activeSystems
  enddo activeFamilies

end function getlabels


!--------------------------------------------------------------------------------------------------
!> @brief Equivalent Poisson's ratio (ν)
!> @details https://doi.org/10.1143/JPSJ.20.635
!--------------------------------------------------------------------------------------------------
function equivalent_nu(C,assumption) result(nu)

  real(pReal), dimension(6,6), intent(in) :: C                                                      !< Stiffness tensor (Voigt notation)
  character(len=*),            intent(in) :: assumption                                             !< Assumption ('Voigt' = isostrain, 'Reuss' = isostress)

  real(pReal)                 :: K, mu, nu
  logical                     :: error
  real(pReal), dimension(6,6) :: S

  if    (IO_lc(assumption) == 'voigt') then
    K = (C(1,1)+C(2,2)+C(3,3) +2.0_pReal*(C(1,2)+C(2,3)+C(1,3))) &
      / 9.0_pReal
  elseif(IO_lc(assumption) == 'reuss') then
    call math_invert(S,error,C)
    if(error) call IO_error(0)
    K = 1.0_pReal &
      / (S(1,1)+S(2,2)+S(3,3) +2.0_pReal*(S(1,2)+S(2,3)+S(1,3)))
  else
    call IO_error(0)
    K = 0.0_pReal
  endif

  mu = equivalent_mu(C,assumption)
  nu = (1.5_pReal*K -mu)/(3.0_pReal*K+mu)

end function equivalent_nu


!--------------------------------------------------------------------------------------------------
!> @brief Equivalent shear modulus (μ)
!> @details https://doi.org/10.1143/JPSJ.20.635
!--------------------------------------------------------------------------------------------------
function equivalent_mu(C,assumption) result(mu)

  real(pReal), dimension(6,6), intent(in) :: C                                                      !< Stiffness tensor (Voigt notation)
  character(len=*),            intent(in) :: assumption                                             !< Assumption ('Voigt' = isostrain, 'Reuss' = isostress)

  real(pReal)                 :: mu
  logical                     :: error
  real(pReal), dimension(6,6) :: S

  if    (IO_lc(assumption) == 'voigt') then
    mu = (1.0_pReal*(C(1,1)+C(2,2)+C(3,3)) -1.0_pReal*(C(1,2)+C(2,3)+C(1,3)) +3.0_pReal*(C(4,4)+C(5,5)+C(6,6))) &
       / 15.0_pReal
  elseif(IO_lc(assumption) == 'reuss') then
    call math_invert(S,error,C)
    if(error) call IO_error(0)
    mu = 15.0_pReal &
       / (4.0_pReal*(S(1,1)+S(2,2)+S(3,3)) -4.0_pReal*(S(1,2)+S(2,3)+S(1,3)) +3.0_pReal*(S(4,4)+S(5,5)+S(6,6)))
  else
    call IO_error(0)
    mu = 0.0_pReal
  endif

end function equivalent_mu


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some lattice functions.
!--------------------------------------------------------------------------------------------------
subroutine selfTest

  real(pReal), dimension(:,:,:), allocatable :: CoSy
  real(pReal), dimension(:,:),   allocatable :: system

  real(pReal), dimension(6,6) :: C
  real(pReal), dimension(2)   :: r
  real(pReal)                 :: lambda

  call random_number(r)

  system = reshape([1.0_pReal+r(1),0.0_pReal,0.0_pReal, 0.0_pReal,1.0_pReal+r(2),0.0_pReal],[6,1])
  CoSy   = buildCoordinateSystem([1],[1],system,'fcc',0.0_pReal)
  if(any(dNeq(CoSy(1:3,1:3,1),math_I3))) error stop 'buildCoordinateSystem'

  call random_number(C)
  C(1,1) = C(1,1) + 1.0_pReal
  C = applyLatticeSymmetryC66(C,'iso')
  if(dNeq(C(6,6),equivalent_mu(C,'voigt'),1.0e-12_pReal)) error stop 'equivalent_mu/voigt'
  if(dNeq(C(6,6),equivalent_mu(C,'voigt'),1.0e-12_pReal)) error stop 'equivalent_mu/reuss'
  
  lambda = C(1,2)
  if(dNeq(lambda*0.5_pReal/(lambda+equivalent_mu(C,'voigt')),equivalent_nu(C,'voigt'),1.0e-12_pReal)) &
    error stop 'equivalent_nu/voigt'
  if(dNeq(lambda*0.5_pReal/(lambda+equivalent_mu(C,'reuss')),equivalent_nu(C,'reuss'),1.0e-12_pReal)) &
    error stop 'equivalent_nu/reuss'

end subroutine selfTest

end module lattice
