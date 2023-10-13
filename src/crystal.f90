!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Contains crystal definitions including Schmid matrices for slip, twin, trans,
!          and cleavage as well as interaction among the various systems.
!--------------------------------------------------------------------------------------------------
module crystal
  use prec
  use misc
  use IO
  use config
  use math
  use rotations

  implicit none(type,external)
  private

!--------------------------------------------------------------------------------------------------
! cF: face centered cubic (fcc)

  integer, dimension(*), parameter :: &
    CF_NSLIPSYSTEM = [12, 6]                                                                        !< # of slip systems per family for cF

  integer, dimension(*), parameter :: &
    CF_NTWINSYSTEM = [12]                                                                           !< # of twin systems per family for cF

  integer, dimension(*), parameter :: &
    CF_NTRANSSYSTEM = [12]                                                                          !< # of transformation systems per family for cF

  integer, dimension(*), parameter :: &
    CF_NCLEAVAGESYSTEM = [3]                                                                        !< # of cleavage systems per family for cF

  integer, parameter  :: &
    CF_NSLIP     = sum(CF_NSLIPSYSTEM), &                                                           !< total # of slip systems for cF
    CF_NTWIN     = sum(CF_NTWINSYSTEM), &                                                           !< total # of twin systems for cF
    CF_NTRANS    = sum(CF_NTRANSSYSTEM), &                                                          !< total # of transformation systems for cF
    CF_NCLEAVAGE = sum(CF_NCLEAVAGESYSTEM)                                                          !< total # of cleavage systems for cF

  real(pREAL), dimension(3+3,CF_NSLIP), parameter :: &
    CF_SYSTEMSLIP = reshape(real([&
    ! <110>{111} systems (Thompson tetrahedron labeling according to Fig. 3 of 10.1016/S1572-4859(05)80003-8)
       0, 1,-1,     1, 1, 1, & ! AC(d)
      -1, 0, 1,     1, 1, 1, & ! CB(d)
       1,-1, 0,     1, 1, 1, & ! BA(d)
       0,-1,-1,    -1,-1, 1, & ! BD(c)
       1, 0, 1,    -1,-1, 1, & ! DA(c)
      -1, 1, 0,    -1,-1, 1, & ! AB(c)
       0,-1, 1,     1,-1,-1, & ! CA(b)
      -1, 0,-1,     1,-1,-1, & ! AD(b)
       1, 1, 0,     1,-1,-1, & ! DC(b)
       0, 1, 1,    -1, 1,-1, & ! DB(a)
       1, 0,-1,    -1, 1,-1, & ! BC(a)
      -1,-1, 0,    -1, 1,-1, & ! CD(a)
     ! <110>{110}/non-octahedral systems
       1, 1, 0,     1,-1, 0, &
       1,-1, 0,     1, 1, 0, &
       1, 0, 1,     1, 0,-1, &
       1, 0,-1,     1, 0, 1, &
       0, 1, 1,     0, 1,-1, &
       0, 1,-1,     0, 1, 1  &
      ],pREAL),shape(CF_SYSTEMSLIP))                                                                !< cF slip systems

  real(pREAL), dimension(3+3,CF_NTWIN), parameter :: &
    CF_SYSTEMTWIN = reshape(real( [&
    ! <112>{111} systems
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
      ],pREAL),shape(CF_SYSTEMTWIN))                                                                !< cF twin systems

  integer, dimension(2,CF_NTWIN), parameter, public :: &
    crystal_CF_TWINNUCLEATIONSLIPPAIR = reshape( [&
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
      ],shape(crystal_CF_TWINNUCLEATIONSLIPPAIR))

  real(pREAL), dimension(3+3,CF_NCLEAVAGE), parameter :: &
    CF_SYSTEMCLEAVAGE = reshape(real([&
    ! <001>{001} systems
       0, 1, 0,     1, 0, 0, &
       0, 0, 1,     0, 1, 0, &
       1, 0, 0,     0, 0, 1  &
      ],pREAL),shape(CF_SYSTEMCLEAVAGE))                                                            !< cF cleavage systems

!--------------------------------------------------------------------------------------------------
! cI: body centered cubic (bcc)

  integer, dimension(*), parameter :: &
    CI_NSLIPSYSTEM = [12, 12, 24]                                                                   !< # of slip systems per family for cI

  integer, dimension(*), parameter :: &
    CI_NTWINSYSTEM = [12]                                                                           !< # of twin systems per family for cI

  integer, dimension(*), parameter :: &
    CI_NCLEAVAGESYSTEM = [3]                                                                        !< # of cleavage systems per family for cI

  integer, parameter  :: &
    CI_NSLIP     = sum(CI_NSLIPSYSTEM), &                                                           !< total # of slip systems for cI
    CI_NTWIN     = sum(CI_NTWINSYSTEM), &                                                           !< total # of twin systems for cI
    CI_NCLEAVAGE = sum(CI_NCLEAVAGESYSTEM)                                                          !< total # of cleavage systems for cI

  real(pREAL), dimension(3+3,CI_NSLIP), parameter :: &
    CI_SYSTEMSLIP = reshape(real([&
    ! <111>{110} systems
    ! Sign convention follows Table 1 of 10.1016/j.ijplas.2020.102733
    ! to allow for universal calculation of non-glide plane normal n1 = Rot(-m,60°) @ n
    ! The choice matters since Rot(-m,60°) @ n ≠ Rot(m,60°) @ -n ..!
       1,-1, 1,     0, 1, 1, &
      -1,-1, 1,     0,-1,-1, &
       1, 1, 1,     0, 1,-1, &
      -1, 1, 1,     0,-1, 1, &
      -1, 1, 1,    -1, 0,-1, &
      -1,-1, 1,     1, 0, 1, &
       1, 1, 1,    -1, 0, 1, &
       1,-1, 1,     1, 0,-1, &
      -1, 1, 1,     1, 1, 0, &
       1,-1, 1,    -1,-1, 0, &
       1, 1, 1,     1,-1, 0, &
      -1,-1, 1,    -1, 1, 0, &
     ! <111>{112} systems
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
       1, 1, 1,     1, 1,-2, &
     ! Slip system <111>{123}
       1, 1,-1,     1, 2, 3, &
       1,-1, 1,    -1, 2, 3, &
      -1, 1, 1,     1,-2, 3, &
       1, 1, 1,     1, 2,-3, &
       1,-1, 1,     1, 3, 2, &
       1, 1,-1,    -1, 3, 2, &
       1, 1, 1,     1,-3, 2, &
      -1, 1, 1,     1, 3,-2, &
       1, 1,-1,     2, 1, 3, &
       1,-1, 1,    -2, 1, 3, &
      -1, 1, 1,     2,-1, 3, &
       1, 1, 1,     2, 1,-3, &
       1,-1, 1,     2, 3, 1, &
       1, 1,-1,    -2, 3, 1, &
       1, 1, 1,     2,-3, 1, &
      -1, 1, 1,     2, 3,-1, &
      -1, 1, 1,     3, 1, 2, &
       1, 1, 1,    -3, 1, 2, &
       1, 1,-1,     3,-1, 2, &
       1,-1, 1,     3, 1,-2, &
      -1, 1, 1,     3, 2, 1, &
       1, 1, 1,    -3, 2, 1, &
       1, 1,-1,     3,-2, 1, &
       1,-1, 1,     3, 2,-1  &
      ],pREAL),shape(CI_SYSTEMSLIP))                                                                !< cI slip systems

  real(pREAL), dimension(3+3,CI_NTWIN), parameter :: &
    CI_SYSTEMTWIN = reshape(real([&
    ! <111>{112} systems
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
      ],pREAL),shape(CI_SYSTEMTWIN))                                                                !< cI twin systems

  real(pREAL), dimension(3+3,CI_NCLEAVAGE), parameter :: &
    CI_SYSTEMCLEAVAGE = reshape(real([&
    ! <001>{001} systems
       0, 1, 0,     1, 0, 0, &
       0, 0, 1,     0, 1, 0, &
       1, 0, 0,     0, 0, 1  &
      ],pREAL),shape(CI_SYSTEMCLEAVAGE))                                                            !< cI cleavage systems

!--------------------------------------------------------------------------------------------------
! hP: hexagonal [close packed] (hex, hcp)

  integer, dimension(*), parameter :: &
    HP_NSLIPSYSTEM = [3, 3, 6, 12, 6]                                                               !< # of slip systems per family for hP

  integer, dimension(*), parameter :: &
    HP_NTWINSYSTEM = [6, 6, 6, 6]                                                                   !< # of slip systems per family for hP

  integer, parameter  :: &
    HP_NSLIP     = sum(HP_NSLIPSYSTEM), &                                                           !< total # of slip systems for hP
    HP_NTWIN     = sum(HP_NTWINSYSTEM)                                                              !< total # of twin systems for hP

  real(pREAL), dimension(4+4,HP_NSLIP), parameter :: &
    HP_SYSTEMSLIP = reshape(real([&
    ! <-1-1.0>{00.1}/basal systems (independent of c/a-ratio)
       2, -1, -1,  0,     0,  0,  0,  1, &
      -1,  2, -1,  0,     0,  0,  0,  1, &
      -1, -1,  2,  0,     0,  0,  0,  1, &
    ! <-1-1.0>{1-1.0}/prismatic systems (independent of c/a-ratio)
       2, -1, -1,  0,     0,  1, -1,  0, &
      -1,  2, -1,  0,    -1,  0,  1,  0, &
      -1, -1,  2,  0,     1, -1,  0,  0, &
    ! <-1-1.0>{-11.1}/1. order pyramidal <a> systems (direction independent of c/a-ratio)
      -1,  2, -1,  0,     1,  0, -1,  1, &
      -2,  1,  1,  0,     0,  1, -1,  1, &
      -1, -1,  2,  0,    -1,  1,  0,  1, &
       1, -2,  1,  0,    -1,  0,  1,  1, &
       2, -1, -1,  0,     0, -1,  1,  1, &
       1,  1, -2,  0,     1, -1,  0,  1, &
    ! <11.3>{-10.1}/1. order pyramidal <c+a> systems (direction independent of c/a-ratio)
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
    ! <11.3>{-1-1.2}/2. order pyramidal <c+a> systems
      -1, -1,  2,  3,     1,  1, -2,  2, &
       1, -2,  1,  3,    -1,  2, -1,  2, &
       2, -1, -1,  3,    -2,  1,  1,  2, &
       1,  1, -2,  3,    -1, -1,  2,  2, &
      -1,  2, -1,  3,     1, -2,  1,  2, &
      -2,  1,  1,  3,     2, -1, -1,  2  &
      ],pREAL),shape(HP_SYSTEMSLIP))                                                                !< hP slip systems, sorted by P. Eisenlohr CCW around <c> starting next to a_1 axis

  real(pREAL), dimension(4+4,HP_NTWIN), parameter :: &
    HP_SYSTEMTWIN =  reshape(real([&
    ! <-10.1>{10.2} systems, shear = (3-(c/a)^2)/(sqrt(3) c/a)
    ! tension in Co, Mg, Zr, Ti, and Be; compression in Cd and Zn
      -1,  0,  1,  1,     1,  0, -1,  2, & !
       0, -1,  1,  1,     0,  1, -1,  2, &
       1, -1,  0,  1,    -1,  1,  0,  2, &
       1,  0, -1,  1,    -1,  0,  1,  2, &
       0,  1, -1,  1,     0, -1,  1,  2, &
      -1,  1,  0,  1,     1, -1,  0,  2, &
    ! <11.6>{-1-1.1} systems, shear = 1/(c/a)
    ! tension in Co, Re, and Zr
      -1, -1,  2,  6,     1,  1, -2,  1, &
       1, -2,  1,  6,    -1,  2, -1,  1, &
       2, -1, -1,  6,    -2,  1,  1,  1, &
       1,  1, -2,  6,    -1, -1,  2,  1, &
      -1,  2, -1,  6,     1, -2,  1,  1, &
      -2,  1,  1,  6,     2, -1, -1,  1, &
    ! <10.-2>{10.1} systems, shear = (4(c/a)^2-9)/(4 sqrt(3) c/a)
    ! compression in Mg
       1,  0, -1, -2,     1,  0, -1,  1, &
       0,  1, -1, -2,     0,  1, -1,  1, &
      -1,  1,  0, -2,    -1,  1,  0,  1, &
      -1,  0,  1, -2,    -1,  0,  1,  1, &
       0, -1,  1, -2,     0, -1,  1,  1, &
       1, -1,  0, -2,     1, -1,  0,  1, &
    ! <11.-3>{11.2} systems, shear = 2((c/a)^2-2)/(3 c/a)
    ! compression in Ti and Zr
       1,  1, -2, -3,     1,  1, -2,  2, &
      -1,  2, -1, -3,    -1,  2, -1,  2, &
      -2,  1,  1, -3,    -2,  1,  1,  2, &
      -1, -1,  2, -3,    -1, -1,  2,  2, &
       1, -2,  1, -3,     1, -2,  1,  2, &
       2, -1, -1, -3,     2, -1, -1,  2  &
      ],pREAL),shape(HP_SYSTEMTWIN))                                                                !< hP twin systems, sorted by P. Eisenlohr CCW around <c> starting next to a_1 axis

!--------------------------------------------------------------------------------------------------
! tI: body centered tetragonal (bct)

  integer, dimension(*), parameter :: &
    TI_NSLIPSYSTEM = [2, 2, 2, 4, 2, 4, 2, 2, 4, 8, 4, 8, 8 ]                                       !< # of slip systems per family for tI

  integer, parameter :: &
    TI_NSLIP = sum(TI_NSLIPSYSTEM)                                                                 !< total # of slip systems for tI

  real(pREAL), dimension(3+3,TI_NSLIP), parameter :: &
    TI_SYSTEMSLIP = reshape(real([&
    ! {100)<001] systems
       0, 0, 1,      1, 0, 0, &
       0, 0, 1,      0, 1, 0, &
    ! {110)<001] systems
       0, 0, 1,      1, 1, 0, &
       0, 0, 1,     -1, 1, 0, &
    ! {100)<010] systems
       0,  1, 0,     1, 0, 0, &
       1,  0, 0,     0, 1, 0, &
    ! {110)<1-11]/2 systems
       1,-1, 1,      1, 1, 0, &
       1,-1,-1,      1, 1, 0, &
      -1,-1,-1,     -1, 1, 0, &
      -1,-1, 1,     -1, 1, 0, &
    ! {110)<1-10] systems
       1, -1, 0,     1, 1, 0, &
       1,  1, 0,     1,-1, 0, &
    ! {100)<011] systems
       0, 1, 1,      1, 0, 0, &
       0,-1, 1,      1, 0, 0, &
      -1, 0, 1,      0, 1, 0, &
       1, 0, 1,      0, 1, 0, &
    ! {001)<010] systems
       0, 1, 0,      0, 0, 1, &
       1, 0, 0,      0, 0, 1, &
    ! {001)<110] systems
       1, 1, 0,      0, 0, 1, &
      -1, 1, 0,      0, 0, 1, &
    ! {011)<01-1] systems
       0, 1,-1,      0, 1, 1, &
       0,-1,-1,      0,-1, 1, &
      -1, 0,-1,     -1, 0, 1, &
       1, 0,-1,      1, 0, 1, &
    ! {011)<1-11]/2 systems
       1,-1, 1,      0, 1, 1, &
       1, 1,-1,      0, 1, 1, &
       1, 1, 1,      0, 1,-1, &
      -1, 1, 1,      0, 1,-1, &
       1,-1,-1,      1, 0, 1, &
      -1,-1, 1,      1, 0, 1, &
       1, 1, 1,      1, 0,-1, &
       1,-1, 1,      1, 0,-1,  &
    ! {011)<100] systems
       1, 0, 0,      0, 1, 1, &
       1, 0, 0,      0, 1,-1, &
       0, 1, 0,      1, 0, 1, &
       0, 1, 0,      1, 0,-1, &
    ! {211)<01-1] systems
       0, 1,-1,      2, 1, 1, &
       0,-1,-1,      2,-1, 1, &
       1, 0,-1,      1, 2, 1, &
      -1, 0,-1,     -1, 2, 1, &
       0, 1,-1,     -2, 1, 1, &
       0,-1,-1,     -2,-1, 1, &
      -1, 0,-1,     -1,-2, 1, &
       1, 0,-1,      1,-2, 1, &
    ! {211)<-111]/2 systems
      -1, 1, 1,      2, 1, 1, &
      -1,-1, 1,      2,-1, 1, &
       1,-1, 1,      1, 2, 1, &
      -1,-1, 1,     -1, 2, 1, &
       1, 1, 1,     -2, 1, 1, &
       1,-1, 1,     -2,-1, 1, &
      -1, 1, 1,     -1,-2, 1, &
       1, 1, 1,      1,-2, 1  &
       ],pREAL),shape(TI_SYSTEMSLIP))                                                               !< tI slip systems for c/a = 0.5456 (Sn), sorted by Bieler 2009 (https://doi.org/10.1007/s11664-009-0909-x)


  interface crystal_forestProjection_edge
    module procedure slipProjection_transverse
  end interface crystal_forestProjection_edge

  interface crystal_forestProjection_screw
    module procedure slipProjection_direction
  end interface crystal_forestProjection_screw

  public :: &
    crystal_init, &
    crystal_selfTest, &
    crystal_isotropic_nu, &
    crystal_isotropic_mu, &
    crystal_symmetrize_33, &
    crystal_symmetrize_C66, &
    crystal_SchmidMatrix_slip, &
    crystal_SchmidMatrix_twin, &
    crystal_SchmidMatrix_trans, &
    crystal_SchmidMatrix_cleavage, &
    crystal_interaction_SlipBySlip, &
    crystal_interaction_TwinByTwin, &
    crystal_interaction_TransByTrans, &
    crystal_interaction_SlipByTwin, &
    crystal_interaction_SlipByTrans, &
    crystal_interaction_TwinBySlip, &
    crystal_characteristicShear_Twin, &
    crystal_C66_twin, &
    crystal_C66_trans, &
    crystal_forestProjection_edge, &
    crystal_forestProjection_screw, &
    crystal_slip_normal, &
    crystal_slip_direction, &
    crystal_slip_transverse, &
    crystal_labels_slip, &
    crystal_labels_twin

contains

!--------------------------------------------------------------------------------------------------
!> @brief Run self test.
!--------------------------------------------------------------------------------------------------
subroutine crystal_init()

  print'(/,1x,a)', '<<<+-  crystal init  -+>>>'; flush(IO_STDOUT)

  call crystal_selfTest()

end subroutine crystal_init


!--------------------------------------------------------------------------------------------------
!> @brief Characteristic shear for twinning
!--------------------------------------------------------------------------------------------------
function crystal_characteristicShear_Twin(Ntwin,lattice,CoverA) result(characteristicShear)

  integer,     dimension(:),            intent(in) :: Ntwin                                         !< number of active twin systems per family
  character(len=*),                     intent(in) :: lattice                                       !< Bravais lattice (Pearson symbol)
  real(pREAL),                          intent(in) :: cOverA                                        !< c/a ratio
  real(pREAL), dimension(sum(Ntwin))               :: characteristicShear

  integer :: &
    a, &                                                                                            !< index of active system
    p, &                                                                                            !< index in potential system list
    f, &                                                                                            !< index of my family
    s                                                                                               !< index of my system in current family

  integer, dimension(HP_NTWIN), parameter :: &
    HP_SHEARTWIN = reshape( [&
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
      ],[HP_NTWIN])                                                                                 !< indicator to formulas below

  a = 0
  myFamilies: do f = 1,size(Ntwin,1)
    mySystems: do s = 1,Ntwin(f)
      a = a + 1
      select case(lattice)
        case('cF','cI')
          characteristicShear(a) = 0.5_pREAL*sqrt(2.0_pREAL)
        case('hP')
          if (cOverA < 1.0_pREAL .or. cOverA > 2.0_pREAL) &
            call IO_error(131,ext_msg='crystal_characteristicShear_Twin')
          p = sum(HP_NTWINSYSTEM(1:f-1))+s
          select case(HP_SHEARTWIN(p))                                                              ! from Christian & Mahajan 1995 p.29
            case (1)                                                                                ! <-10.1>{10.2}
              characteristicShear(a) = (3.0_pREAL-cOverA**2)/sqrt(3.0_pREAL)/CoverA
            case (2)                                                                                ! <11.6>{-1-1.1}
              characteristicShear(a) = 1.0_pREAL/cOverA
            case (3)                                                                                ! <10.-2>{10.1}
              characteristicShear(a) = (4.0_pREAL*cOverA**2-9.0_pREAL)/sqrt(48.0_pREAL)/cOverA
            case (4)                                                                                ! <11.-3>{11.2}
              characteristicShear(a) = 2.0_pREAL*(cOverA**2-2.0_pREAL)/3.0_pREAL/cOverA
          end select
        case default
          call IO_error(137,ext_msg='crystal_characteristicShear_Twin: '//trim(lattice))
      end select
    end do mySystems
  end do myFamilies

end function crystal_characteristicShear_Twin


!--------------------------------------------------------------------------------------------------
!> @brief Rotated elasticity matrices for twinning in 6x6-matrix notation
!--------------------------------------------------------------------------------------------------
function crystal_C66_twin(Ntwin,C66,lattice,CoverA)

  integer,     dimension(:),            intent(in) :: Ntwin                                         !< number of active twin systems per family
  character(len=*),                     intent(in) :: lattice                                       !< Bravais lattice (Pearson symbol)
  real(pREAL), dimension(6,6),          intent(in) :: C66                                           !< unrotated parent stiffness matrix
  real(pREAL),                          intent(in) :: cOverA                                        !< c/a ratio
  real(pREAL), dimension(6,6,sum(Ntwin))           :: crystal_C66_twin

  real(pREAL), dimension(3,3,sum(Ntwin)):: coordinateSystem
  type(tRotation)                       :: R
  integer                               :: i


  select case(lattice)
    case('cF')
      coordinateSystem = buildCoordinateSystem(Ntwin,CF_NSLIPSYSTEM,CF_SYSTEMTWIN,&
                                               lattice,0.0_pREAL)
    case('cI')
      coordinateSystem = buildCoordinateSystem(Ntwin,CI_NSLIPSYSTEM,CI_SYSTEMTWIN,&
                                               lattice,0.0_pREAL)
    case('hP')
      coordinateSystem = buildCoordinateSystem(Ntwin,HP_NSLIPSYSTEM,HP_SYSTEMTWIN,&
                                               lattice,cOverA)
    case default
      call IO_error(137,ext_msg='crystal_C66_twin: '//trim(lattice))
  end select

  do i = 1, sum(Ntwin)
    call R%fromAxisAngle([coordinateSystem(1:3,2,i),PI],P=1)                                        ! ToDo: Why always 180 deg?
    crystal_C66_twin(1:6,1:6,i) = R%rotStiffness(C66)
  end do

end function crystal_C66_twin


!--------------------------------------------------------------------------------------------------
!> @brief Rotated elasticity matrices for transformation in 6x6-matrix notation
!--------------------------------------------------------------------------------------------------
function crystal_C66_trans(Ntrans,C_parent66,crystal_target, &
                           cOverA_trans,a_cF,a_cI)

  integer,     dimension(:),             intent(in) :: Ntrans                                       !< number of active twin systems per family
  character(len=*),                      intent(in) :: crystal_target                               !< Bravais lattice (Pearson symbol)
  real(pREAL), dimension(6,6),           intent(in) :: C_parent66
  real(pREAL),                 optional, intent(in) :: cOverA_trans, a_cF, a_cI
  real(pREAL), dimension(6,6,sum(Ntrans))           :: crystal_C66_trans

  real(pREAL), dimension(6,6)             :: C_bar66, C_target_unrotated66
  real(pREAL), dimension(3,3,sum(Ntrans)) :: Q,S
  type(tRotation)                         :: R
  integer                                 :: i

 !--------------------------------------------------------------------------------------------------
 ! elasticity matrix of the target phase in cube orientation
  if (crystal_target == 'hP' .and. present(cOverA_trans)) then
    ! https://doi.org/10.1063/1.1663858 eq. (16), eq. (18), eq. (19)
    ! https://doi.org/10.1016/j.actamat.2016.07.032 eq. (47), eq. (48)
    if (cOverA_trans < 1.0_pREAL .or. cOverA_trans > 2.0_pREAL) &
      call IO_error(131,ext_msg='crystal_C66_trans: '//trim(crystal_target))
    C_bar66(1,1) = (C_parent66(1,1) + C_parent66(1,2) + 2.0_pREAL*C_parent66(4,4))/2.0_pREAL
    C_bar66(1,2) = (C_parent66(1,1) + 5.0_pREAL*C_parent66(1,2) - 2.0_pREAL*C_parent66(4,4))/6.0_pREAL
    C_bar66(3,3) = (C_parent66(1,1) + 2.0_pREAL*C_parent66(1,2) + 4.0_pREAL*C_parent66(4,4))/3.0_pREAL
    C_bar66(1,3) = (C_parent66(1,1) + 2.0_pREAL*C_parent66(1,2) - 2.0_pREAL*C_parent66(4,4))/3.0_pREAL
    C_bar66(4,4) = (C_parent66(1,1) - C_parent66(1,2) + C_parent66(4,4))/3.0_pREAL
    C_bar66(1,4) = (C_parent66(1,1) - C_parent66(1,2) - 2.0_pREAL*C_parent66(4,4)) /(3.0_pREAL*sqrt(2.0_pREAL))

    C_target_unrotated66 = 0.0_pREAL
    C_target_unrotated66(1,1) = C_bar66(1,1) - C_bar66(1,4)**2/C_bar66(4,4)
    C_target_unrotated66(1,2) = C_bar66(1,2) + C_bar66(1,4)**2/C_bar66(4,4)
    C_target_unrotated66(1,3) = C_bar66(1,3)
    C_target_unrotated66(3,3) = C_bar66(3,3)
    C_target_unrotated66(4,4) = C_bar66(4,4) - C_bar66(1,4)**2/(0.5_pREAL*(C_bar66(1,1) - C_bar66(1,2)))
    C_target_unrotated66 = crystal_symmetrize_C66(C_target_unrotated66,'hP')
  elseif (crystal_target  == 'cI' .and. present(a_cF) .and. present(a_cI)) then
    if (a_cI <= 0.0_pREAL .or. a_cF <= 0.0_pREAL) &
      call IO_error(134,ext_msg='crystal_C66_trans: '//trim(crystal_target))
    C_target_unrotated66 = C_parent66
  else
    call IO_error(137,ext_msg='crystal_C66_trans : '//trim(crystal_target))
  end if

  do i = 1,6
    if (abs(C_target_unrotated66(i,i))<tol_math_check) &
    call IO_error(135,'matrix diagonal in transformation',label1='entry',ID1=i)
  end do

  call buildTransformationSystem(Q,S,Ntrans,cOverA_trans,a_cF,a_cI)

  do i = 1,sum(Ntrans)
    call R%fromMatrix(Q(1:3,1:3,i))
    crystal_C66_trans(1:6,1:6,i) = R%rotStiffness(C_target_unrotated66)
  end do

 end function crystal_C66_trans


!--------------------------------------------------------------------------------------------------
!> @brief Slip-slip interaction matrix
!> @details only active slip systems are considered
!> @details https://doi.org/10.1016/j.actamat.2016.12.040 (cF: Tab S4-1, cI: Tab S5-1)
!> @details https://doi.org/10.1016/j.ijplas.2014.06.010 (hP: Tab 3b)
!--------------------------------------------------------------------------------------------------
function crystal_interaction_SlipBySlip(Nslip,interactionValues,lattice) result(interactionMatrix)

  integer,         dimension(:),                   intent(in) :: Nslip                              !< number of active slip systems per family
  real(pREAL),     dimension(:),                   intent(in) :: interactionValues                  !< values for slip-slip interaction
  character(len=*),                                intent(in) :: lattice                            !< Bravais lattice (Pearson symbol)
  real(pREAL),     dimension(sum(Nslip),sum(Nslip))           :: interactionMatrix

  integer, dimension(:),   allocatable :: NslipMax
  integer, dimension(:,:), allocatable :: interactionTypes


  integer, dimension(CF_NSLIP,CF_NSLIP), parameter :: &
    CF_INTERACTIONSLIPSLIP = reshape( [&
       1, 2, 2, 4, 7, 5, 3, 5, 5, 4, 6, 7,  10,11,10,11,12,13, & ! -----> acting (forest)
       2, 1, 2, 7, 4, 5, 6, 4, 7, 5, 3, 5,  10,11,12,13,10,11, & ! |
       2, 2, 1, 5, 5, 3, 6, 7, 4, 7, 6, 4,  12,13,10,11,10,11, & ! |
       4, 7, 6, 1, 2, 2, 4, 6, 7, 3, 5, 5,  10,11,11,10,13,12, & ! v
       7, 4, 6, 2, 1, 2, 5, 3, 5, 6, 4, 7,  10,11,13,12,11,10, & ! reacting (primary)
       5, 5, 3, 2, 2, 1, 7, 6, 4, 6, 7, 4,  12,13,11,10,11,10, &
       3, 5, 5, 4, 6, 7, 1, 2, 2, 4, 7, 6,  11,10,11,10,12,13, &
       6, 4, 7, 5, 3, 5, 2, 1, 2, 7, 4, 6,  11,10,13,12,10,11, &
       6, 7, 4, 7, 6, 4, 2, 2, 1, 5, 5, 3,  13,12,11,10,10,11, &
       4, 6, 7, 3, 5, 5, 4, 7, 6, 1, 2, 2,  11,10,10,11,13,12, &
       5, 3, 5, 6, 4, 7, 7, 4, 6, 2, 1, 2,  11,10,12,13,11,10, &
       7, 6, 4, 6, 7, 4, 5, 5, 3, 2, 2, 1,  13,12,10,11,11,10, &

      10,10,12,10,10,12,11,11,13,11,11,13,   1, 8, 9, 9, 9, 9, &
      11,11,13,11,11,13,10,10,12,10,10,12,   8, 1, 9, 9, 9, 9, &
      10,12,10,11,13,11,11,13,11,10,12,10,   9, 9, 1, 8, 9, 9, &
      11,13,11,10,12,10,10,12,10,11,13,11,   9, 9, 8, 1, 9, 9, &
      12,10,10,13,11,11,12,10,10,13,11,11,   9, 9, 9, 9, 1, 8, &
      13,11,11,12,10,10,13,11,11,12,10,10,   9, 9, 9, 9, 8, 1  &
      ],shape(CF_INTERACTIONSLIPSLIP))                                                              !< Slip-slip interaction types for cF / Madec 2017 (https://doi.org/10.1016/j.actamat.2016.12.040)
                                                                                                    !< 1: self interaction         --> alpha 0
                                                                                                    !< 2: coplanar interaction     --> alpha copla
                                                                                                    !< 3: collinear interaction    --> alpha coli
                                                                                                    !< 4: Hirth locks              --> alpha 1
                                                                                                    !< 5: glissile junctions I     --> alpha 2
                                                                                                    !< 6: glissile junctions II    --> alpha 2*
                                                                                                    !< 7: Lomer locks              --> alpha 3
                                                                                                    !< 8: crossing (similar to Hirth locks in <110>{111} for two {110} planes)
                                                                                                    !< 9: similar to Lomer locks in <110>{111} for two {110} planes
                                                                                                    !<10: similar to Lomer locks in <110>{111} btw one {110} and one {111} plane
                                                                                                    !<11: similar to glissile junctions in <110>{111} btw one {110} and one {111} plane
                                                                                                    !<12: crossing btw one {110} and one {111} plane
                                                                                                    !<13: collinear btw one {110} and one {111} plane

  integer, dimension(CI_NSLIP,CI_NSLIP), parameter :: &
    CI_INTERACTIONSLIPSLIP = reshape( [&
       1, 3, 6, 6, 7, 5, 4, 2, 4, 2, 7, 5,  18,18,11, 8, 9,13,17,14,13, 9,17,14,  28,25,28,28,25,28,28,28,28,25,28,28,25,28,28,28,28,28,28,25,28,28,28,25, &! -----> acting (forest)
       3, 1, 6, 6, 4, 2, 7, 5, 7, 5, 4, 2,  18,18, 8,11,13, 9,14,17, 9,13,14,17,  25,28,28,28,28,25,28,28,25,28,28,28,28,25,28,28,28,28,25,28,28,28,25,28, &! |
       6, 6, 1, 3, 5, 7, 2, 4, 5, 7, 2, 4,  11, 8,18,18,17,14, 9,13,17,14,13, 9,  28,28,28,25,28,28,25,28,28,28,28,25,28,28,25,28,28,25,28,28,28,25,28,28, &! |
       6, 6, 3, 1, 2, 4, 5, 7, 2, 4, 5, 7,   8,11,18,18,14,17,13, 9,14,17, 9,13,  28,28,25,28,28,28,28,25,28,28,25,28,28,28,28,25,25,28,28,28,25,28,28,28, &! v
       7, 5, 4, 2, 1, 3, 6, 6, 2, 4, 7, 5,   9,17,13,14,18,11,18, 8,13,17, 9,14,  28,28,25,28,28,28,28,25,28,28,25,28,28,28,28,25,25,28,28,28,25,28,28,28, &! reacting (primary)
       4, 2, 7, 5, 3, 1, 6, 6, 5, 7, 4, 2,  13,14, 9,17,18, 8,18,11, 9,14,13,17,  25,28,28,28,28,25,28,28,25,28,28,28,28,25,28,28,28,28,25,28,28,28,25,28, &
       5, 7, 2, 4, 6, 6, 1, 3, 7, 5, 2, 4,  17, 9,14,13,11,18, 8,18,17,13,14, 9,  28,28,28,25,28,28,25,28,28,28,28,25,28,28,25,28,28,25,28,28,28,25,28,28, &
       2, 4, 5, 7, 6, 6, 3, 1, 4, 2, 5, 7,  14,13,17, 9, 8,18,11,18,14, 9,17,13,  28,25,28,28,25,28,28,28,28,25,28,28,25,28,28,28,28,28,28,25,28,28,28,25, &
       5, 7, 4, 2, 2, 4, 7, 5, 1, 3, 6, 6,   9,17,14,13,13,17,14, 9,18,11, 8,18,  28,28,25,28,28,28,28,25,28,28,25,28,28,28,28,25,25,28,28,28,25,28,28,28, &
       2, 4, 7, 5, 5, 7, 4, 2, 3, 1, 6, 6,  13,14,17, 9, 9,14,17,13,18, 8,11,18,  28,25,28,28,25,28,28,28,28,25,28,28,25,28,28,28,28,28,28,25,28,28,28,25, &
       7, 5, 2, 4, 7, 5, 2, 4, 6, 6, 1, 3,  17, 9,13,14,17,13, 9,14,11,18,18, 8,  28,28,28,25,28,28,25,28,28,28,28,25,28,28,25,28,28,25,28,28,28,25,28,28, &
       4, 2, 5, 7, 4, 2, 5, 7, 6, 6, 3, 1,  14,13, 9,17,14, 9,13,17, 8,18,18,11,  25,28,28,28,28,25,28,28,25,28,28,28,28,25,28,28,28,28,25,28,28,28,25,28, &

      19,19,10, 8, 9,12,16,15, 9,12,16,15,   1,20,24,24,23,22,21, 2,23,22, 2,21,  28,28,26,28,28,28,28,26,28,28,26,28,28,28,28,26,26,28,28,28,26,28,28,28, &
      19,19, 8,10,16,15, 9,12,16,15, 9,12,  20, 1,24,24,22,23, 2,21,22,23,21, 2,  28,28,28,26,28,28,26,28,28,28,28,26,28,28,26,28,28,26,28,28,28,26,28,28, &
      10, 8,19,19,12, 9,15,16,15,16,12, 9,  24,24, 1,20,21, 2,23,22, 2,21,23,22,  26,28,28,28,28,26,28,28,26,28,28,28,28,26,28,28,28,28,26,28,28,28,26,28, &
       8,10,19,19,15,16,12, 9,12, 9,15,16,  24,24,20, 1, 2,21,22,23,21, 2,22,23,  28,26,28,28,26,28,28,28,28,26,28,28,26,28,28,28,28,28,28,26,28,28,28,26, &
       9,12,16,15,19,19,10, 8,12, 9,16,15,  23,21,22, 2, 1,24,20,24,23, 2,22,21,  28,26,28,28,26,28,28,28,28,26,28,28,26,28,28,28,28,28,28,26,28,28,28,26, &
      12, 9,15,16,10, 8,19,19,16,15,12, 9,  21,23, 2,21,24, 1,24,20, 2,23,21,22,  26,28,28,28,28,26,28,28,26,28,28,28,28,26,28,28,28,28,26,28,28,28,26,28, &
      16,15, 9,12,19,19, 8,10,15,16, 9,12,  22, 2,23,22,20,24, 1,24,22,21,23, 2,  28,28,28,26,28,28,26,28,28,28,28,26,28,28,26,28,28,26,28,28,28,26,28,28, &
      15,16,12, 9, 8,10,19,19, 9,12,15,16,   2,22,21,23,24,20,24, 1,21,22, 2,23,  28,28,26,28,28,28,28,26,28,28,26,28,28,28,28,26,26,28,28,28,26,28,28,28, &
      12, 9,16,15,12, 9,16,15,19,19,10, 8,  23,21, 2,22,23, 2,21,22, 1,24,24,20,  26,28,28,28,28,26,28,28,26,28,28,28,28,26,28,28,28,28,26,28,28,28,26,28, &
       9,12,15,16,16,15,12, 9,10, 8,19,19,  21,23,22, 2, 2,23,22,21,24, 1,20,24,  28,26,28,28,26,28,28,28,28,26,28,28,26,28,28,28,28,28,28,26,28,28,28,26, &
      16,15,12, 9, 9,12,15,16, 8,10,19,19,   2,22,23,21,21,22,23, 2,24,20, 1,24,  28,28,26,28,28,28,28,26,28,28,26,28,28,28,28,26,26,28,28,28,26,28,28,28, &
      15,16, 9,12,15,16, 9,12,19,19, 8,10,  22, 2,21,23,22,21, 2,23,20,24,24, 1,  28,28,28,26,28,28,26,28,28,28,28,26,28,28,26,28,28,26,28,28,28,26,28,28, &

      28,25,28,28,28,25,28,28,28,28,28,25,  28,28,26,28,28,26,28,28,26,28,28,28,   1,28,28,28,28,27,28,28,27,28,28,28,28,27,28,28,28,28,27,28,28,28,27,28, &
      25,28,28,28,28,28,28,25,28,25,28,28,  28,28,28,26,26,28,28,28,28,26,28,28,  28, 1,28,28,27,28,28,28,28,27,28,28,27,28,28,28,28,28,28,27,28,28,28,27, &
      28,28,28,25,25,28,28,28,25,28,28,28,  26,28,28,28,28,28,28,26,28,28,26,28,  28,28, 1,28,28,28,28,27,28,28,27,28,28,28,28,27,27,28,28,28,27,28,28,28, &
      28,28,25,28,28,28,25,28,28,28,25,28,  28,26,28,28,28,28,26,28,28,28,28,26,  28,28,28, 1,28,28,27,28,28,28,28,27,28,28,27,28,28,27,28,28,28,27,28,28, &
      25,28,28,28,28,28,28,25,28,25,28,28,  28,28,28,26,26,28,28,28,28,26,28,28,  28,27,28,28, 1,28,28,28,28,27,28,28,27,28,28,28,28,28,28,27,28,28,28,27, &
      28,25,28,28,28,25,28,28,28,28,28,25,  28,28,26,28,28,26,28,28,26,28,28,28,  27,28,28,28,28, 1,28,28,27,28,28,28,28,27,28,28,28,28,27,28,28,28,27,28, &
      28,28,25,28,28,28,25,28,28,28,25,28,  28,26,28,28,28,28,26,28,28,28,28,26,  28,28,28,27,28,28, 1,28,28,28,28,27,28,28,27,28,28,27,28,28,28,27,28,28, &
      28,28,28,25,25,28,28,28,25,28,28,28,  26,28,28,28,28,28,28,26,28,28,26,28,  28,28,27,28,28,28,28, 1,28,28,27,28,28,28,28,27,27,28,28,28,27,28,28,28, &
      28,25,28,28,28,25,28,28,28,28,28,25,  28,28,26,28,28,26,28,28,26,28,28,28,  27,28,28,28,28,27,28,28, 1,28,28,28,28,27,28,28,28,28,27,28,28,28,27,28, &
      25,28,28,28,28,28,28,25,28,25,28,28,  28,28,28,26,26,28,28,28,28,26,28,28,  28,27,28,28,27,28,28,28,28, 1,28,28,27,28,28,28,28,28,28,27,28,28,28,27, &
      28,28,28,25,25,28,28,28,25,28,28,28,  26,28,28,28,28,28,28,26,28,28,26,28,  28,28,27,28,28,28,28,27,28,28, 1,28,28,28,28,27,27,28,28,28,27,28,28,28, &
      28,28,25,28,28,28,25,28,28,28,25,28,  28,26,28,28,28,28,26,28,28,28,28,26,  28,28,28,27,28,28,27,28,28,28,28, 1,28,28,27,28,28,27,28,28,28,27,28,28, &
      25,28,28,28,28,28,28,25,28,25,28,28,  28,28,28,26,26,28,28,28,28,26,28,28,  28,27,28,28,27,28,28,28,28,27,28,28, 1,28,28,28,28,28,28,27,28,28,28,27, &
      28,25,28,28,28,25,28,28,28,28,28,25,  28,28,26,28,28,26,28,28,26,28,28,28,  27,28,28,28,28,27,28,28,27,28,28,28,28, 1,28,28,28,28,27,28,28,28,27,28, &
      28,28,25,28,28,28,25,28,28,28,25,28,  28,26,28,28,28,28,26,28,28,28,28,26,  28,28,28,27,28,28,27,28,28,28,28,27,28,28, 1,28,28,27,28,28,28,27,28,28, &
      28,28,28,25,25,28,28,28,25,28,28,28,  26,28,28,28,28,28,28,26,28,28,26,28,  28,28,27,28,28,28,28,27,28,28,27,28,28,28,28, 1,27,28,28,28,27,28,28,28, &
      28,28,28,25,25,28,28,28,25,28,28,28,  26,28,28,28,28,28,28,26,28,28,26,28,  28,28,27,28,28,28,28,27,28,28,27,28,28,28,28,27, 1,28,28,28,27,28,28,28, &
      28,28,25,28,28,28,25,28,28,28,25,28,  28,26,28,28,28,28,26,28,28,28,28,26,  28,28,28,27,28,28,27,28,28,28,28,27,28,28,27,28,28, 1,28,28,28,27,28,28, &
      28,25,28,28,28,25,28,28,28,28,28,25,  28,28,26,28,28,26,28,28,26,28,28,28,  27,28,28,28,28,27,28,28,27,28,28,28,28,27,28,28,28,28, 1,28,28,28,27,28, &
      25,28,28,28,28,28,28,25,28,25,28,28,  28,28,28,26,26,28,28,28,28,26,28,28,  28,27,28,28,27,28,28,28,28,27,28,28,27,28,28,28,28,28,28, 1,28,28,28,27, &
      28,28,28,25,25,28,28,28,25,28,28,28,  26,28,28,28,28,28,28,26,28,28,26,28,  28,28,27,28,28,28,28,27,28,28,27,28,28,28,28,27,27,28,28,28, 1,28,28,28, &
      28,28,25,28,28,28,25,28,28,28,25,28,  28,26,28,28,28,28,26,28,28,28,28,26,  28,28,28,27,28,28,27,28,28,28,28,27,28,28,27,28,28,27,28,28,28, 1,28,28, &
      28,25,28,28,28,25,28,28,28,28,28,25,  28,28,26,28,28,26,28,28,26,28,28,28,  27,28,28,28,28,27,28,28,27,28,28,28,28,27,28,28,28,28,27,28,28,28, 1,28, &
      25,28,28,28,28,28,28,25,28,25,28,28,  28,28,28,26,26,28,28,28,28,26,28,28,  28,27,28,28,27,28,28,28,28,27,28,28,27,28,28,28,28,28,28,27,28,28,28, 1  &
      ],shape(CI_INTERACTIONSLIPSLIP))                                                              !< Slip-slip interaction types for cI / Madec 2017 (https://doi.org/10.1016/j.actamat.2016.12.040)
                                                                                                    !< 1: self interaction      --> alpha 0
                                                                                                    !< 2: collinear interaction --> alpha 1
                                                                                                    !< 3: coplanar interaction  --> alpha 2
                                                                                                    !< 4-7: other coefficients
                                                                                                    !< 8: {110}-{112} collinear and perpendicular planes --> alpha 6
                                                                                                    !< 9: {110}-{112} collinear                          --> alpha 7
                                                                                                    !< 10-24: other coefficients
                                                                                                    !< 25: {110}-{123} collinear
                                                                                                    !< 26: {112}-{123} collinear
                                                                                                    !< 27: {123}-{123} collinear
                                                                                                    !< 28: other interaction

  integer, dimension(HP_NSLIP,HP_NSLIP), parameter :: &
    HP_INTERACTIONSLIPSLIP = reshape( [&
    ! basal      prism      1. pyr<a>           1. pyr<c+a>                           2. pyr<c+a>
       1, 2, 2,   3, 4, 4,   9,10, 9, 9,10, 9,  20,21,22,22,21,20,20,21,22,22,21,20,  47,47,48,47,47,48, & ! -----> acting (forest)
       2, 1, 2,   4, 3, 4,  10, 9, 9,10, 9, 9,  22,22,21,20,20,21,22,22,21,20,20,21,  47,48,47,47,48,47, & ! | basal
       2, 2, 1,   4, 4, 3,   9, 9,10, 9, 9,10,  21,20,20,21,22,22,21,20,20,21,22,22,  48,47,47,48,47,47, & ! |
                                                                                                           ! v
       7, 8, 8,   5, 6, 6,  11,12,11,11,12,11,  23,24,25,25,24,23,23,24,25,25,24,23,  49,49,50,49,49,50, & ! reacting (primary)
       8, 7, 8,   6, 5, 6,  12,11,11,12,11,11,  25,25,24,23,23,24,25,25,24,23,23,24,  49,50,49,49,50,49, & ! prism
       8, 8, 7,   6, 6, 5,  11,11,12,11,11,12,  24,23,23,24,25,25,24,23,23,24,25,25,  50,49,49,50,49,49, &

      18,19,18,  16,17,16,  13,14,14,15,14,14,  26,26,27,28,28,27,29,29,27,28,28,27,  51,52,51,51,52,51, &
      19,18,18,  17,16,16,  14,13,14,14,15,14,  28,27,26,26,27,28,28,27,29,29,27,28,  51,51,52,51,51,52, &
      18,18,19,  16,16,17,  14,14,13,14,14,15,  27,28,28,27,26,26,27,28,28,27,29,29,  52,51,51,52,51,51, &
      18,19,18,  16,17,16,  15,14,14,13,14,14,  29,29,27,28,28,27,26,26,27,28,28,27,  51,52,51,51,52,51, & ! 1. pyr<a>
      19,18,18,  17,16,16,  14,15,14,14,13,14,  28,27,29,29,27,28,28,27,26,26,27,28,  51,51,52,51,51,52, &
      18,18,19,  16,16,17,  14,14,15,14,14,13,  27,28,28,27,29,29,27,28,28,27,26,26,  52,51,51,52,51,51, &

      44,45,46,  41,42,43,  37,38,39,40,38,39,  30,31,32,32,32,33,34,35,32,32,32,36,  53,54,55,53,54,56, &
      46,45,44,  43,42,41,  37,39,38,40,39,38,  31,30,36,32,32,32,35,34,33,32,32,32,  56,54,53,55,54,53, &
      45,46,44,  42,43,41,  39,37,38,39,40,38,  32,36,30,31,32,32,32,33,34,35,32,32,  56,53,54,55,53,54, &
      45,44,46,  42,41,43,  38,37,39,38,40,39,  32,32,31,30,36,32,32,32,35,34,33,32,  53,56,54,53,55,54, &
      46,44,45,  43,41,42,  38,39,37,38,39,40,  32,32,32,36,30,31,32,32,32,33,34,35,  54,56,53,54,55,53, &
      44,46,45,  41,43,42,  39,38,37,39,38,40,  33,32,32,32,31,30,36,32,32,32,35,34,  54,53,56,54,53,55, &
      44,45,46,  41,42,43,  40,38,39,37,38,39,  34,35,32,32,32,36,30,31,32,32,32,33,  53,54,56,53,54,55, & ! 1. pyr<c+a>
      46,45,44,  43,42,41,  40,39,38,37,39,38,  35,34,33,32,32,32,31,30,36,32,32,32,  55,54,53,56,54,53, &
      45,46,44,  42,43,41,  39,40,38,39,37,38,  32,33,34,35,32,32,32,36,30,31,32,32,  55,53,54,56,53,54, &
      45,44,46,  42,41,43,  38,40,39,38,37,39,  32,32,35,34,33,32,32,32,31,30,36,32,  53,55,54,53,56,54, &
      46,44,45,  43,41,42,  38,39,40,38,39,37,  32,32,32,33,34,35,32,32,32,36,30,31,  54,55,53,54,56,53, &
      44,46,45,  41,43,42,  39,38,40,39,38,37,  36,32,32,32,35,34,33,32,32,32,31,30,  54,53,55,54,53,56, &

      68,68,69,  66,66,67,  64,64,65,64,65,65,  60,61,61,60,62,62,60,63,63,60,62,62,  57,58,58,59,58,58, &
      68,69,68,  66,67,66,  65,64,64,65,64,64,  62,62,60,61,61,60,62,62,60,63,63,60,  58,57,58,58,59,58, &
      69,68,68,  67,66,66,  64,65,64,64,65,64,  63,60,62,62,60,61,61,60,62,62,60,63,  58,58,57,58,58,59, &
      68,68,69,  66,66,67,  64,64,65,64,64,65,  60,63,63,60,62,62,60,61,61,60,62,62,  59,58,58,57,58,58, & ! 2. pyr<c+a>
      68,69,68,  66,67,66,  65,64,64,65,64,64,  62,62,60,63,63,60,62,62,60,61,61,60,  58,59,58,58,57,58, &
      69,68,68,  67,66,66,  64,65,64,64,65,64,  61,60,62,62,60,63,63,60,62,62,60,61,  58,58,59,58,58,57  &
      ],shape(HP_INTERACTIONSLIPSLIP))                                                              !< Slip-slip interaction types for hP (onion peel naming scheme)
                                                                                                    !< 10.1016/j.ijplas.2014.06.010 table 3
                                                                                                    !< 10.1080/14786435.2012.699689 table 2 and 3
                                                                                                    !< index & label & description
                                                                                                    !<  1 & S1 & basal self-interaction
                                                                                                    !<  2 &  1 & basal/basal coplanar
                                                                                                    !<  3 &  3 & basal/prismatic collinear
                                                                                                    !<  4 &  4 & basal/prismatic non-collinear
                                                                                                    !<  5 & S2 & prismatic self-interaction
                                                                                                    !<  6 &  2 & prismatic/prismatic
                                                                                                    !<  7 &  5 & prismatic/basal collinear
                                                                                                    !<  8 &  6 & prismatic/basal non-collinear
                                                                                                    !<  9 &  - & basal/pyramidal <a> non-collinear
                                                                                                    !< 10 &  - & basal/pyramidal <a> collinear
                                                                                                    !< 11 &  - & prismatic/pyramidal <a> non-collinear
                                                                                                    !< 12 &  - & prismatic/pyramidal <a> collinear
                                                                                                    !< 13 &  - & pyramidal <a> self-interaction
                                                                                                    !< 14 &  - & pyramidal <a> non-collinear
                                                                                                    !< 15 &  - & pyramidal <a> collinear
                                                                                                    !< 16 &  - & pyramidal <a>/prismatic non-collinear
                                                                                                    !< 17 &  - & pyramidal <a>/prismatic collinear
                                                                                                    !< 18 &  - & pyramidal <a>/basal non-collinear
                                                                                                    !< 19 &  - & pyramidal <a>/basal collinear
                                                                                                    !< 20 &  - & basal/1. order pyramidal <c+a> semi-collinear
                                                                                                    !< 21 &  - & basal/1. order pyramidal <c+a>
                                                                                                    !< 22 &  - & basal/1. order pyramidal <c+a>
                                                                                                    !< 23 &  - & prismatic/1. order pyramidal <c+a> semi-collinear
                                                                                                    !< 24 &  - & prismatic/1. order pyramidal <c+a>
                                                                                                    !< 25 &  - & prismatic/1. order pyramidal <c+a> semi-coplanar?
                                                                                                    !< 26 &  - & pyramidal <a>/1. order pyramidal <c+a> coplanar
                                                                                                    !< 27 &  - & pyramidal <a>/1. order pyramidal <c+a>
                                                                                                    !< 28 &  - & pyramidal <a>/1. order pyramidal <c+a> semi-collinear
                                                                                                    !< 29 &  - & pyramidal <a>/1. order pyramidal <c+a> semi-coplanar
                                                                                                    !< 30 &  - & 1. order pyramidal <c+a> self-interaction
                                                                                                    !< 31 &  - & 1. order pyramidal <c+a> coplanar
                                                                                                    !< 32 &  - & 1. order pyramidal <c+a>
                                                                                                    !< 33 &  - & 1. order pyramidal <c+a>
                                                                                                    !< 34 &  - & 1. order pyramidal <c+a> semi-coplanar
                                                                                                    !< 35 &  - & 1. order pyramidal <c+a> semi-coplanar
                                                                                                    !< 36 &  - & 1. order pyramidal <c+a> collinear
                                                                                                    !< 37 &  - & 1. order pyramidal <c+a>/pyramidal <a> coplanar
                                                                                                    !< 38 &  - & 1. order pyramidal <c+a>/pyramidal <a> semi-collinear
                                                                                                    !< 39 &  - & 1. order pyramidal <c+a>/pyramidal <a>
                                                                                                    !< 40 &  - & 1. order pyramidal <c+a>/pyramidal <a> semi-coplanar
                                                                                                    !< 41 &  - & 1. order pyramidal <c+a>/prismatic semi-collinear
                                                                                                    !< 42 &  - & 1. order pyramidal <c+a>/prismatic semi-coplanar
                                                                                                    !< 43 &  - & 1. order pyramidal <c+a>/prismatic
                                                                                                    !< 44 &  - & 1. order pyramidal <c+a>/basal semi-collinear
                                                                                                    !< 45 &  - & 1. order pyramidal <c+a>/basal
                                                                                                    !< 46 &  - & 1. order pyramidal <c+a>/basal
                                                                                                    !< 47 &  8 & basal/2. order pyramidal <c+a> non-collinear
                                                                                                    !< 48 &  7 & basal/2. order pyramidal <c+a> semi-collinear
                                                                                                    !< 49 & 10 & prismatic/2. order pyramidal <c+a>
                                                                                                    !< 50 &  9 & prismatic/2. order pyramidal <c+a> semi-collinear
                                                                                                    !< 51 &  - & pyramidal <a>/2. order pyramidal <c+a>
                                                                                                    !< 52 &  - & pyramidal <a>/2. order pyramidal <c+a> semi collinear
                                                                                                    !< 53 &  - & 1. order pyramidal <c+a>/2. order pyramidal <c+a>
                                                                                                    !< 54 &  - & 1. order pyramidal <c+a>/2. order pyramidal <c+a>
                                                                                                    !< 55 &  - & 1. order pyramidal <c+a>/2. order pyramidal <c+a>
                                                                                                    !< 56 &  - & 1. order pyramidal <c+a>/2. order pyramidal <c+a> collinear
                                                                                                    !< 57 & S3 & 2. order pyramidal <c+a> self-interaction
                                                                                                    !< 58 & 16 & 2. order pyramidal <c+a> non-collinear
                                                                                                    !< 59 & 15 & 2. order pyramidal <c+a> semi-collinear
                                                                                                    !< 60 &  - & 2. order pyramidal <c+a>/1. order pyramidal <c+a>
                                                                                                    !< 61 &  - & 2. order pyramidal <c+a>/1. order pyramidal <c+a> collinear
                                                                                                    !< 62 &  - & 2. order pyramidal <c+a>/1. order pyramidal <c+a>
                                                                                                    !< 63 &  - & 2. order pyramidal <c+a>/1. order pyramidal <c+a>
                                                                                                    !< 64 &  - & 2. order pyramidal <c+a>/pyramidal <a> non-collinear
                                                                                                    !< 65 &  - & 2. order pyramidal <c+a>/pyramidal <a> semi-collinear
                                                                                                    !< 66 & 14 & 2. order pyramidal <c+a>/prismatic non-collinear
                                                                                                    !< 67 & 13 & 2. order pyramidal <c+a>/prismatic semi-collinear
                                                                                                    !< 68 & 12 & 2. order pyramidal <c+a>/basal non-collinear
                                                                                                    !< 69 & 11 & 2. order pyramidal <c+a>/basal semi-collinear

  integer, dimension(TI_NSLIP,TI_NSLIP), parameter :: &
    TI_INTERACTIONSLIPSLIP = reshape( [&
        1,  2,   3,  3,   7,  7,  13, 13, 13, 13,  21, 21,  31, 31, 31, 31,  43, 43,  57, 57,  73, 73, 73, 73,  91, 91, 91, 91, 91, 91, 91, 91, 111, 111, 111, 111, 133,133,133,133,133,133,133,133, 157,157,157,157,157,157,157,157, & ! -----> acting
        2,  1,   3,  3,   7,  7,  13, 13, 13, 13,  21, 21,  31, 31, 31, 31,  43, 43,  57, 57,  73, 73, 73, 73,  91, 91, 91, 91, 91, 91, 91, 91, 111, 111, 111, 111, 133,133,133,133,133,133,133,133, 157,157,157,157,157,157,157,157, & ! |
                                                                                                                                                                                                                                        ! |
        6,  6,   4,  5,   8,  8,  14, 14, 14, 14,  22, 22,  32, 32, 32, 32,  44, 44,  58, 58,  74, 74, 74, 74,  92, 92, 92, 92, 92, 92, 92, 92, 112, 112, 112, 112, 134,134,134,134,134,134,134,134, 158,158,158,158,158,158,158,158, & ! v
        6,  6,   5,  4,   8,  8,  14, 14, 14, 14,  22, 22,  32, 32, 32, 32,  44, 44,  58, 58,  74, 74, 74, 74,  92, 92, 92, 92, 92, 92, 92, 92, 112, 112, 112, 112, 134,134,134,134,134,134,134,134, 158,158,158,158,158,158,158,158, & ! reacting

       12, 12,  11, 11,   9, 10,  15, 15, 15, 15,  23, 23,  33, 33, 33, 33,  45, 45,  59, 59,  75, 75, 75, 75,  93, 93, 93, 93, 93, 93, 93, 93, 113, 113, 113, 113, 135,135,135,135,135,135,135,135, 159,159,159,159,159,159,159,159, &
       12, 12,  11, 11,  10,  9,  15, 15, 15, 15,  23, 23,  33, 33, 33, 33,  45, 45,  59, 59,  75, 75, 75, 75,  93, 93, 93, 93, 93, 93, 93, 93, 113, 113, 113, 113, 135,135,135,135,135,135,135,135, 159,159,159,159,159,159,159,159, &

       20, 20,  19, 19,  18, 18,  16, 17, 17, 17,  24, 24,  34, 34, 34, 34,  46, 46,  60, 60,  76, 76, 76, 76,  94, 94, 94, 94, 94, 94, 94, 94, 114, 114, 114, 114, 136,136,136,136,136,136,136,136, 160,160,160,160,160,160,160,160, &
       20, 20,  19, 19,  18, 18,  17, 16, 17, 17,  24, 24,  34, 34, 34, 34,  46, 46,  60, 60,  76, 76, 76, 76,  94, 94, 94, 94, 94, 94, 94, 94, 114, 114, 114, 114, 136,136,136,136,136,136,136,136, 160,160,160,160,160,160,160,160, &
       20, 20,  19, 19,  18, 18,  17, 17, 16, 17,  24, 24,  34, 34, 34, 34,  46, 46,  60, 60,  76, 76, 76, 76,  94, 94, 94, 94, 94, 94, 94, 94, 114, 114, 114, 114, 136,136,136,136,136,136,136,136, 160,160,160,160,160,160,160,160, &
       20, 20,  19, 19,  18, 18,  17, 17, 17, 16,  24, 24,  34, 34, 34, 34,  46, 46,  60, 60,  76, 76, 76, 76,  94, 94, 94, 94, 94, 94, 94, 94, 114, 114, 114, 114, 136,136,136,136,136,136,136,136, 160,160,160,160,160,160,160,160, &

       30, 30,  29, 29,  28, 28,  27, 27, 27, 27,  25, 26,  35, 35, 35, 35,  47, 47,  61, 61,  77, 77, 77, 77,  95, 95, 95, 95, 95, 95, 95, 95, 115, 115, 115, 115, 137,137,137,137,137,137,137,137, 161,161,161,161,161,161,161,161, &
       30, 30,  29, 29,  28, 28,  27, 27, 27, 27,  26, 25,  35, 35, 35, 35,  47, 47,  61, 61,  77, 77, 77, 77,  95, 95, 95, 95, 95, 95, 95, 95, 115, 115, 115, 115, 137,137,137,137,137,137,137,137, 161,161,161,161,161,161,161,161, &

       42, 42,  41, 41,  40, 40,  39, 39, 39, 39,  38, 38,  36, 37, 37, 37,  48, 48,  62, 62,  78, 78, 78, 78,  96, 96, 96, 96, 96, 96, 96, 96, 116, 116, 116, 116, 138,138,138,138,138,138,138,138, 162,162,162,162,162,162,162,162, &
       42, 42,  41, 41,  40, 40,  39, 39, 39, 39,  38, 38,  37, 36, 37, 37,  48, 48,  62, 62,  78, 78, 78, 78,  96, 96, 96, 96, 96, 96, 96, 96, 116, 116, 116, 116, 138,138,138,138,138,138,138,138, 162,162,162,162,162,162,162,162, &
       42, 42,  41, 41,  40, 40,  39, 39, 39, 39,  38, 38,  37, 37, 36, 37,  48, 48,  62, 62,  78, 78, 78, 78,  96, 96, 96, 96, 96, 96, 96, 96, 116, 116, 116, 116, 138,138,138,138,138,138,138,138, 162,162,162,162,162,162,162,162, &
       42, 42,  41, 41,  40, 40,  39, 39, 39, 39,  38, 38,  37, 37, 37, 36,  48, 48,  62, 62,  78, 78, 78, 78,  96, 96, 96, 96, 96, 96, 96, 96, 116, 116, 116, 116, 138,138,138,138,138,138,138,138, 162,162,162,162,162,162,162,162, &

       56, 56,  55, 55,  54, 54,  53, 53, 53, 53,  52, 52,  51, 51, 51, 51,  49, 50,  63, 63,  79, 79, 79, 79,  97, 97, 97, 97, 97, 97, 97, 97, 117, 117, 117, 117, 139,139,139,139,139,139,139,139, 163,163,163,163,163,163,163,163, &
       56, 56,  55, 55,  54, 54,  53, 53, 53, 53,  52, 52,  51, 51, 51, 51,  50, 49,  63, 63,  79, 79, 79, 79,  97, 97, 97, 97, 97, 97, 97, 97, 117, 117, 117, 117, 139,139,139,139,139,139,139,139, 163,163,163,163,163,163,163,163, &

       72, 72,  71, 71,  70, 70,  69, 69, 69, 69,  68, 68,  67, 67, 67, 67,  66, 66,  64, 65,  80, 80, 80, 80,  98, 98, 98, 98, 98, 98, 98, 98, 118, 118, 118, 118, 140,140,140,140,140,140,140,140, 164,164,164,164,164,164,164,164, &
       72, 72,  71, 71,  70, 70,  69, 69, 69, 69,  68, 68,  67, 67, 67, 67,  66, 66,  65, 64,  80, 80, 80, 80,  98, 98, 98, 98, 98, 98, 98, 98, 118, 118, 118, 118, 140,140,140,140,140,140,140,140, 164,164,164,164,164,164,164,164, &

       90, 90,  89, 89,  88, 88,  87, 87, 87, 87,  86, 86,  85, 85, 85, 85,  84, 84,  83, 83,  81, 82, 82, 82,  99, 99, 99, 99, 99, 99, 99, 99, 119, 119, 119, 119, 141,141,141,141,141,141,141,141, 165,165,165,165,165,165,165,165, &
       90, 90,  89, 89,  88, 88,  87, 87, 87, 87,  86, 86,  85, 85, 85, 85,  84, 84,  83, 83,  82, 81, 82, 82,  99, 99, 99, 99, 99, 99, 99, 99, 119, 119, 119, 119, 141,141,141,141,141,141,141,141, 165,165,165,165,165,165,165,165, &
       90, 90,  89, 89,  88, 88,  87, 87, 87, 87,  86, 86,  85, 85, 85, 85,  84, 84,  83, 83,  82, 82, 81, 82,  99, 99, 99, 99, 99, 99, 99, 99, 119, 119, 119, 119, 141,141,141,141,141,141,141,141, 165,165,165,165,165,165,165,165, &
       90, 90,  89, 89,  88, 88,  87, 87, 87, 87,  86, 86,  85, 85, 85, 85,  84, 84,  83, 83,  82, 82, 82, 81,  99, 99, 99, 99, 99, 99, 99, 99, 119, 119, 119, 119, 141,141,141,141,141,141,141,141, 165,165,165,165,165,165,165,165, &

      110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 100,101,101,101,101,101,101,101, 120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
      110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,100,101,101,101,101,101,101, 120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
      110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,100,101,101,101,101,101, 120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
      110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,100,101,101,101,101, 120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
      110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,101,100,101,101,101, 120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
      110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,101,101,100,101,101, 120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
      110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,101,101,101,100,101, 120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &
      110,110, 109,109, 108,108, 107,107,107,107, 106,106, 105,105,105,105, 104,104, 103,103, 102,102,102,102, 101,101,101,101,101,101,101,100, 120, 120, 120, 120, 142,142,142,142,142,142,142,142, 166,166,166,166,166,166,166,166, &

      132,132, 131,131, 130,130, 129,129,129,129, 128,128, 127,127,127,127, 126,126, 125,125, 124,124,124,124, 123,123,123,123,123,123,123,123, 121, 122, 122, 122, 143,143,143,143,143,143,143,143, 167,167,167,167,167,167,167,167, &
      132,132, 131,131, 130,130, 129,129,129,129, 128,128, 127,127,127,127, 126,126, 125,125, 124,124,124,124, 123,123,123,123,123,123,123,123, 121, 121, 122, 122, 143,143,143,143,143,143,143,143, 167,167,167,167,167,167,167,167, &
      132,132, 131,131, 130,130, 129,129,129,129, 128,128, 127,127,127,127, 126,126, 125,125, 124,124,124,124, 123,123,123,123,123,123,123,123, 121, 122, 121, 122, 143,143,143,143,143,143,143,143, 167,167,167,167,167,167,167,167, &
      132,132, 131,131, 130,130, 129,129,129,129, 128,128, 127,127,127,127, 126,126, 125,125, 124,124,124,124, 123,123,123,123,123,123,123,123, 121, 122, 122, 121, 143,143,143,143,143,143,143,143, 167,167,167,167,167,167,167,167, &

      156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147, 146, 146, 146, 146, 144,145,145,145,145,145,145,145, 168,168,168,168,168,168,168,168, &
      156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147, 146, 146, 146, 146, 145,144,145,145,145,145,145,145, 168,168,168,168,168,168,168,168, &
      156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147, 146, 146, 146, 146, 145,145,144,145,145,145,145,145, 168,168,168,168,168,168,168,168, &
      156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147, 146, 146, 146, 146, 145,145,145,144,145,145,145,145, 168,168,168,168,168,168,168,168, &
      156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147, 146, 146, 146, 146, 145,145,145,145,144,145,145,145, 168,168,168,168,168,168,168,168, &
      156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147, 146, 146, 146, 146, 145,145,145,145,145,144,145,145, 168,168,168,168,168,168,168,168, &
      156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147, 146, 146, 146, 146, 145,145,145,145,145,145,144,145, 168,168,168,168,168,168,168,168, &
      156,156, 155,155, 154,154, 153,153,153,153, 152,152, 151,151,151,151, 150,150, 149,149, 148,148,148,148, 147,147,147,147,147,147,147,147, 146, 146, 146, 146, 145,145,145,145,145,145,145,144, 168,168,168,168,168,168,168,168, &

      182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173, 172, 172, 172, 172, 171,171,171,171,171,171,171,171, 169,170,170,170,170,170,170,170, &
      182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173, 172, 172, 172, 172, 171,171,171,171,171,171,171,171, 170,169,170,170,170,170,170,170, &
      182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173, 172, 172, 172, 172, 171,171,171,171,171,171,171,171, 170,170,169,170,170,170,170,170, &
      182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173, 172, 172, 172, 172, 171,171,171,171,171,171,171,171, 170,170,170,169,170,170,170,170, &
      182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173, 172, 172, 172, 172, 171,171,171,171,171,171,171,171, 170,170,170,170,169,170,170,170, &
      182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173, 172, 172, 172, 172, 171,171,171,171,171,171,171,171, 169,170,170,170,170,169,170,170, &
      182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173, 172, 172, 172, 172, 171,171,171,171,171,171,171,171, 169,170,170,170,170,170,169,170, &
      182,182, 181,181, 180,180, 179,179,179,179, 178,178, 177,177,177,177, 176,176, 175,175, 174,174,174,174, 173,173,173,173,173,173,173,173, 172, 172, 172, 172, 171,171,171,171,171,171,171,171, 169,170,170,170,170,170,170,169  &
      ],shape(TI_INTERACTIONSLIPSLIP))


  select case(lattice)
    case('cF')
      interactionTypes = CF_INTERACTIONSLIPSLIP
      NslipMax         = CF_NSLIPSYSTEM
    case('cI')
      interactionTypes = CI_INTERACTIONSLIPSLIP
      NslipMax         = CI_NSLIPSYSTEM
    case('hP')
      interactionTypes = HP_INTERACTIONSLIPSLIP
      NslipMax         = HP_NSLIPSYSTEM
    case('tI')
      interactionTypes = TI_INTERACTIONSLIPSLIP
      NslipMax         = TI_NSLIPSYSTEM
    case default
      call IO_error(137,ext_msg='crystal_interaction_SlipBySlip: '//trim(lattice))
  end select

  interactionMatrix = buildInteraction(Nslip,Nslip,NslipMax,NslipMax,interactionValues,interactionTypes)

end function crystal_interaction_SlipBySlip


!--------------------------------------------------------------------------------------------------
!> @brief Twin-twin interaction matrix
!> details only active twin systems are considered
!--------------------------------------------------------------------------------------------------
function crystal_interaction_TwinByTwin(Ntwin,interactionValues,lattice) result(interactionMatrix)

  integer,         dimension(:),                   intent(in) :: Ntwin                              !< number of active twin systems per family
  real(pREAL),     dimension(:),                   intent(in) :: interactionValues                  !< values for twin-twin interaction
  character(len=*),                                intent(in) :: lattice                            !< Bravais lattice (Pearson symbol)
  real(pREAL),     dimension(sum(Ntwin),sum(Ntwin))           :: interactionMatrix

  integer, dimension(:),   allocatable :: NtwinMax
  integer, dimension(:,:), allocatable :: interactionTypes

  integer, dimension(CF_NTWIN,CF_NTWIN), parameter :: &
    CF_INTERACTIONTWINTWIN = reshape( [&
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
      ],shape(CF_INTERACTIONTWINTWIN))                                                              !< Twin-twin interaction types for cF

  integer, dimension(CI_NTWIN,CI_NTWIN), parameter :: &
    CI_INTERACTIONTWINTWIN = reshape( [&
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
      ],shape(CI_INTERACTIONTWINTWIN))                                                              !< Twin-twin interaction types for cI
                                                                                                    !< 1: self interaction
                                                                                                    !< 2: collinear interaction
                                                                                                    !< 3: other interaction
  integer, dimension(HP_NTWIN,HP_NTWIN), parameter :: &
    HP_INTERACTIONTWINTWIN = reshape( [&
    ! <-10.1>{10.2}       <11.6>{-1-1.1}      <10.-2>{10.1}       <11.-3>{11.2}
       1, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, & ! -----> acting
       2, 1, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, & ! |
       2, 2, 1, 2, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, & ! |
       2, 2, 2, 1, 2, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, & ! v <-10.1>{10.2}
       2, 2, 2, 2, 1, 2,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, & ! reacting
       2, 2, 2, 2, 2, 1,   3, 3, 3, 3, 3, 3,   7, 7, 7, 7, 7, 7,  13,13,13,13,13,13, &

       6, 6, 6, 6, 6, 6,   4, 5, 5, 5, 5, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
       6, 6, 6, 6, 6, 6,   5, 4, 5, 5, 5, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
       6, 6, 6, 6, 6, 6,   5, 5, 4, 5, 5, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
       6, 6, 6, 6, 6, 6,   5, 5, 5, 4, 5, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, & ! <11.6>{-1-1.1}
       6, 6, 6, 6, 6, 6,   5, 5, 5, 5, 4, 5,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &
       6, 6, 6, 6, 6, 6,   5, 5, 5, 5, 5, 4,   8, 8, 8, 8, 8, 8,  14,14,14,14,14,14, &

      12,12,12,12,12,12,  11,11,11,11,11,11,   9,10,10,10,10,10,  15,15,15,15,15,15, &
      12,12,12,12,12,12,  11,11,11,11,11,11,  10, 9,10,10,10,10,  15,15,15,15,15,15, &
      12,12,12,12,12,12,  11,11,11,11,11,11,  10,10, 9,10,10,10,  15,15,15,15,15,15, &
      12,12,12,12,12,12,  11,11,11,11,11,11,  10,10,10, 9,10,10,  15,15,15,15,15,15, & ! <10.-2>{10.1}
      12,12,12,12,12,12,  11,11,11,11,11,11,  10,10,10,10, 9,10,  15,15,15,15,15,15, &
      12,12,12,12,12,12,  11,11,11,11,11,11,  10,10,10,10,10, 9,  15,15,15,15,15,15, &

      20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  16,17,17,17,17,17, &
      20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,16,17,17,17,17, &
      20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,17,16,17,17,17, &
      20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,17,17,16,17,17, & ! <11.-3>{11.2}
      20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,17,17,17,16,17, &
      20,20,20,20,20,20,  19,19,19,19,19,19,  18,18,18,18,18,18,  17,17,17,17,17,16  &
      ],shape(HP_INTERACTIONTWINTWIN))                                                              !< Twin-twin interaction types for hP

  select case(lattice)
    case('cF')
      interactionTypes = CF_INTERACTIONTWINTWIN
      NtwinMax         = CF_NTWINSYSTEM
    case('cI')
      interactionTypes = CI_INTERACTIONTWINTWIN
      NtwinMax         = CI_NTWINSYSTEM
    case('hP')
      interactionTypes = HP_INTERACTIONTWINTWIN
      NtwinMax         = HP_NTWINSYSTEM
    case default
      call IO_error(137,ext_msg='crystal_interaction_TwinByTwin: '//trim(lattice))
  end select

  interactionMatrix = buildInteraction(Ntwin,Ntwin,NtwinMax,NtwinMax,interactionValues,interactionTypes)

end function crystal_interaction_TwinByTwin


!--------------------------------------------------------------------------------------------------
!> @brief Trans-trans interaction matrix
!> details only active trans systems are considered
!--------------------------------------------------------------------------------------------------
function crystal_interaction_TransByTrans(Ntrans,interactionValues,lattice) result(interactionMatrix)

  integer,         dimension(:),                     intent(in) :: Ntrans                           !< number of active trans systems per family
  real(pREAL),     dimension(:),                     intent(in) :: interactionValues                !< values for trans-trans interaction
  character(len=*),                                  intent(in) :: lattice                          !<Bravais lattice (Pearson symbol) (parent crystal)
  real(pREAL),     dimension(sum(Ntrans),sum(Ntrans))           :: interactionMatrix

  integer, dimension(:),   allocatable :: NtransMax
  integer, dimension(:,:), allocatable :: interactionTypes

  integer, dimension(CF_NTRANS,CF_NTRANS), parameter :: &
    CF_INTERACTIONTRANSTRANS = reshape( [&
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
      ],shape(CF_INTERACTIONTRANSTRANS))                                                            !< Trans-trans interaction types for cF

  if (lattice == 'cF') then
    interactionTypes = CF_INTERACTIONTRANSTRANS
    NtransMax        = CF_NTRANSSYSTEM
  else
    call IO_error(137,ext_msg='crystal_interaction_TransByTrans: '//trim(lattice))
  end if

  interactionMatrix = buildInteraction(Ntrans,Ntrans,NtransMax,NtransMax,interactionValues,interactionTypes)

end function crystal_interaction_TransByTrans


!--------------------------------------------------------------------------------------------------
!> @brief Slip-twin interaction matrix
!> details only active slip and twin systems are considered
!--------------------------------------------------------------------------------------------------
function crystal_interaction_SlipByTwin(Nslip,Ntwin,interactionValues,lattice) result(interactionMatrix)

  integer,         dimension(:),                   intent(in) :: Nslip, &                           !< number of active slip systems per family
                                                                 Ntwin                              !< number of active twin systems per family
  real(pREAL),     dimension(:),                   intent(in) :: interactionValues                  !< values for slip-twin interaction
  character(len=*),                                intent(in) :: lattice                            !< Bravais lattice (Pearson symbol)
  real(pREAL),     dimension(sum(Nslip),sum(Ntwin))           :: interactionMatrix

  integer, dimension(:),   allocatable :: NslipMax, &
                                          NtwinMax
  integer, dimension(:,:), allocatable :: interactionTypes

  integer, dimension(CF_NTWIN,CF_NSLIP), parameter :: &
    CF_INTERACTIONSLIPTWIN = reshape( [&
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
      ],shape(CF_INTERACTIONSLIPTWIN))                                                              !< Slip-twin interaction types for cF
                                                                                                    !< 1: coplanar interaction
                                                                                                    !< 2: screw trace between slip system and twin habit plane (easy cross slip)
                                                                                                    !< 3: other interaction
  integer, dimension(CI_NTWIN,CI_NSLIP), parameter :: &
   CI_INTERACTIONSLIPTWIN = reshape( [&
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
      3,2,3,3,3,3,2,3,3,3,3,1, &

      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4, &
      4,4,4,4,4,4,4,4,4,4,4,4  &
      ],shape(CI_INTERACTIONSLIPTWIN))                                                              !< Slip-twin interaction types for cI
                                                                                                    !< 1: coplanar interaction
                                                                                                    !< 2: screw trace between slip system and twin habit plane (easy cross slip)
                                                                                                    !< 3: other interaction
                                                                                                    !< 4: other interaction with slip family {123}

  integer, dimension(HP_NTWIN,HP_NSLIP), parameter :: &
    HP_INTERACTIONSLIPTWIN = reshape( [&
    ! <-10.1>{10.2}       <11.6>{-1-1.1}      <10.-2>{10.1}       <11.-3>{11.2}
       1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! ----> twin (acting)
       1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! | basal
       1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2,   3, 3, 3, 3, 3, 3,   4, 4, 4, 4, 4, 4, & ! |
                                                                                       ! v
       5, 5, 5, 5, 5, 5,   6, 6, 6, 6, 6, 6,   7, 7, 7, 7, 7, 7,   8, 8, 8, 8, 8, 8, & ! slip (reacting)
       5, 5, 5, 5, 5, 5,   6, 6, 6, 6, 6, 6,   7, 7, 7, 7, 7, 7,   8, 8, 8, 8, 8, 8, & ! prism
       5, 5, 5, 5, 5, 5,   6, 6, 6, 6, 6, 6,   7, 7, 7, 7, 7, 7,   8, 8, 8, 8, 8, 8, &

       9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
       9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
       9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
       9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, & ! 1. pyr<a>
       9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &
       9, 9, 9, 9, 9, 9,  10,10,10,10,10,10,  11,11,11,11,11,11,  12,12,12,12,12,12, &

      13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
      13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
      13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
      13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
      13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
      13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
      13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, & ! 1. pyr<c+a>
      13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
      13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
      13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
      13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &
      13,13,13,13,13,13,  14,14,14,14,14,14,  15,15,15,15,15,15,  16,16,16,16,16,16, &

      17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
      17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
      17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
      17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, & ! 2. pyr<c+a>
      17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20, &
      17,17,17,17,17,17,  18,18,18,18,18,18,  19,19,19,19,19,19,  20,20,20,20,20,20  &
      ],shape(HP_INTERACTIONSLIPTWIN))                                                              !< Slip-twin interaction types for hP

  select case(lattice)
    case('cF')
      interactionTypes = CF_INTERACTIONSLIPTWIN
      NslipMax         = CF_NSLIPSYSTEM
      NtwinMax         = CF_NTWINSYSTEM
    case('cI')
      interactionTypes = CI_INTERACTIONSLIPTWIN
      NslipMax         = CI_NSLIPSYSTEM
      NtwinMax         = CI_NTWINSYSTEM
    case('hP')
      interactionTypes = HP_INTERACTIONSLIPTWIN
      NslipMax         = HP_NSLIPSYSTEM
      NtwinMax         = HP_NTWINSYSTEM
    case default
      call IO_error(137,ext_msg='crystal_interaction_SlipByTwin: '//trim(lattice))
  end select

  interactionMatrix = buildInteraction(Nslip,Ntwin,NslipMax,NtwinMax,interactionValues,interactionTypes)

end function crystal_interaction_SlipByTwin


!--------------------------------------------------------------------------------------------------
!> @brief Slip-trans interaction matrix
!> details only active slip and trans systems are considered
!--------------------------------------------------------------------------------------------------
function crystal_interaction_SlipByTrans(Nslip,Ntrans,interactionValues,lattice) result(interactionMatrix)

  integer,         dimension(:),                    intent(in) :: Nslip, &                          !< number of active slip systems per family
                                                                  Ntrans                            !< number of active trans systems per family
  real(pREAL),     dimension(:),                    intent(in) :: interactionValues                 !< values for slip-trans interaction
  character(len=*),                                 intent(in) :: lattice                           !< Bravais lattice (Pearson symbol) (parent crystal)
  real(pREAL),     dimension(sum(Nslip),sum(Ntrans))           :: interactionMatrix

  integer, dimension(:),   allocatable :: NslipMax, &
                                          NtransMax
  integer, dimension(:,:), allocatable :: interactionTypes

  integer, dimension(CF_NTRANS,CF_NSLIP), parameter :: &
    CF_INTERACTIONSLIPTRANS = reshape( [&
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
      ],shape(CF_INTERACTIONSLIPTRANS))                                                             !< Slip-trans interaction types for cF

  select case(lattice)
    case('cF')
      interactionTypes = CF_INTERACTIONSLIPTRANS
      NslipMax         = CF_NSLIPSYSTEM
      NtransMax        = CF_NTRANSSYSTEM
    case default
      call IO_error(137,ext_msg='crystal_interaction_SlipByTrans: '//trim(lattice))
  end select

  interactionMatrix = buildInteraction(Nslip,Ntrans,NslipMax,NtransMax,interactionValues,interactionTypes)

 end function crystal_interaction_SlipByTrans


!--------------------------------------------------------------------------------------------------
!> @brief Twin-slip interaction matrix
!> details only active twin and slip systems are considered
!--------------------------------------------------------------------------------------------------
function crystal_interaction_TwinBySlip(Ntwin,Nslip,interactionValues,lattice) result(interactionMatrix)

  integer,         dimension(:),                   intent(in) :: Ntwin, &                           !< number of active twin systems per family
                                                                 Nslip                              !< number of active slip systems per family
  real(pREAL),     dimension(:),                   intent(in) :: interactionValues                  !< values for twin-twin interaction
  character(len=*),                                intent(in) :: lattice                            !< Bravais lattice (Pearson symbol)
  real(pREAL),     dimension(sum(Ntwin),sum(Nslip))           :: interactionMatrix

  integer, dimension(:),   allocatable :: NtwinMax, &
                                          NslipMax
  integer, dimension(:,:), allocatable :: interactionTypes

  integer, dimension(CF_NSLIP,CF_NTWIN), parameter :: &
    CF_INTERACTIONTWINSLIP = 1                                                                      !< Twin-slip interaction types for cF

  integer, dimension(CI_NSLIP,CI_NTWIN), parameter :: &
    CI_INTERACTIONTWINSLIP = 1                                                                      !< Twin-slip interaction types for cI

  integer, dimension(HP_NSLIP,HP_NTWIN), parameter :: &
    HP_INTERACTIONTWINSLIP = reshape( [&
    ! basal      prism     1. pyr<a>           1. pyr<c+a>                           2. pyr<c+a>
      1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17, & ! ----> slip (acting)
      1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17, & ! |
      1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17, & ! |
      1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17, & ! v <-10.1>{10.2}
      1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17, & ! twin (reacting)
      1, 1, 1,   5, 5, 5,   9, 9, 9, 9, 9, 9,  13,13,13,13,13,13,13,13,13,13,13,13,  17,17,17,17,17,17, &

      2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14,  18,18,18,18,18,18, &
      2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14,  18,18,18,18,18,18, &
      2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14,  18,18,18,18,18,18, &
      2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14,  18,18,18,18,18,18, & ! <11.6>{-1-1.1}
      2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14,  18,18,18,18,18,18, &
      2, 2, 2,   6, 6, 6,  10,10,10,10,10,10,  14,14,14,14,14,14,14,14,14,14,14,14,  18,18,18,18,18,18, &

      3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15,  19,19,19,19,19,19, &
      3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15,  19,19,19,19,19,19, &
      3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15,  19,19,19,19,19,19, &
      3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15,  19,19,19,19,19,19, & ! <10.-2>{10.1}
      3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15,  19,19,19,19,19,19, &
      3, 3, 3,   7, 7, 7,  11,11,11,11,11,11,  15,15,15,15,15,15,15,15,15,15,15,15,  19,19,19,19,19,19, &

      4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  20,20,20,20,20,20, &
      4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  20,20,20,20,20,20, &
      4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  20,20,20,20,20,20, &
      4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  20,20,20,20,20,20, & ! <11.-3>{11.2}
      4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  20,20,20,20,20,20, &
      4, 4, 4,   8, 8, 8,  12,12,12,12,12,12,  16,16,16,16,16,16,16,16,16,16,16,16,  20,20,20,20,20,20  &
      ],shape(HP_INTERACTIONTWINSLIP))                                                              !< Twin-slip interaction types for hP

  select case(lattice)
    case('cF')
      interactionTypes = CF_INTERACTIONTWINSLIP
      NtwinMax         = CF_NTWINSYSTEM
      NslipMax         = CF_NSLIPSYSTEM
    case('cI')
      interactionTypes = CI_INTERACTIONTWINSLIP
      NtwinMax         = CI_NTWINSYSTEM
      NslipMax         = CI_NSLIPSYSTEM
    case('hP')
      interactionTypes = HP_INTERACTIONTWINSLIP
      NtwinMax         = HP_NTWINSYSTEM
      NslipMax         = HP_NSLIPSYSTEM
    case default
      call IO_error(137,ext_msg='crystal_interaction_TwinBySlip: '//trim(lattice))
  end select

  interactionMatrix = buildInteraction(Ntwin,Nslip,NtwinMax,NslipMax,interactionValues,interactionTypes)

end function crystal_interaction_TwinBySlip


!--------------------------------------------------------------------------------------------------
!> @brief Schmid matrix for slip
!> details only active slip systems are considered
! Non-schmid projections for cI with up to 6 coefficients
! https://doi.org/10.1016/j.actamat.2012.03.053, eq. (17)
! https://doi.org/10.1016/j.actamat.2008.07.037, table 1
!--------------------------------------------------------------------------------------------------
function crystal_SchmidMatrix_slip(Nslip,lattice,cOverA,nonSchmidCoefficients,sense) result(SchmidMatrix)

  integer,     dimension(:),              intent(in) :: Nslip                                       !< number of active slip systems per family
  character(len=*),                       intent(in) :: lattice                                     !< Bravais lattice (Pearson symbol)
  real(pREAL),                            intent(in) :: cOverA
  real(pREAL), dimension(:,:), optional,  intent(in) :: nonSchmidCoefficients                       !< non-Schmid coefficients for projections
  integer,                     optional,  intent(in) :: sense                                       !< sense (-1,+1)
  real(pREAL), dimension(3,3,sum(Nslip))             :: SchmidMatrix

  real(pREAL), dimension(3,3,sum(Nslip))             :: coordinateSystem
  real(pREAL), dimension(:,:),           allocatable :: slipSystems
  integer,     dimension(:),             allocatable :: NslipMax
  integer,     dimension(:),             allocatable :: slipFamily
  real(pREAL), dimension(3)                          :: direction, normal, np
  real(pREAL), dimension(6)                          :: coeff                                       !< local nonSchmid coefficient variable
  type(tRotation)                                    :: R
  integer                                            :: i

  select case(lattice)
    case('cF')
      NslipMax    = CF_NSLIPSYSTEM
      slipSystems = CF_SYSTEMSLIP
    case('cI')
      NslipMax    = CI_NSLIPSYSTEM
      slipSystems = CI_SYSTEMSLIP
    case('hP')
      NslipMax    = HP_NSLIPSYSTEM
      slipSystems = HP_SYSTEMSLIP
    case('tI')
      NslipMax    = TI_NSLIPSYSTEM
      slipSystems = TI_SYSTEMSLIP
    case default
      allocate(NslipMax(0))
      call IO_error(137,ext_msg='crystal_SchmidMatrix_slip: '//trim(lattice))
  end select

  if (any(NslipMax(1:size(Nslip)) - Nslip < 0)) &
    call IO_error(145,ext_msg='Nslip '//trim(lattice))
  if (any(Nslip < 0)) &
    call IO_error(144,ext_msg='Nslip '//trim(lattice))

  slipFamily = math_expand([(i, i=1,size(Nslip))],Nslip)
  coordinateSystem = buildCoordinateSystem(Nslip,NslipMax,slipSystems,lattice,cOverA)
  if (present(sense)) then
    if (abs(sense) /= 1) error stop 'neither +1 nor -1 sense in crystal_SchmidMatrix_slip'
    coordinateSystem(1:3,1,1:sum(Nslip)) = coordinateSystem(1:3,1,1:sum(Nslip)) * real(sense,pREAL)
  end if

  do i = 1,sum(Nslip)
    direction = coordinateSystem(1:3,1,i)
    normal    = coordinateSystem(1:3,2,i)

    SchmidMatrix(1:3,1:3,i) = math_outer(direction,normal)
    if (abs(math_trace33(SchmidMatrix(1:3,1:3,i))) > tol_math_check) &
      error stop 'dilatational Schmid matrix for slip'

    if (present(nonSchmidCoefficients)) then
      select case(lattice)
        case('cI')
          coeff(:) = 0.0_pREAL
          select case(slipFamily(i))
            case(1)
              coeff(:size(nonSchmidCoefficients(1,:))) = nonSchmidCoefficients(1,:)
              call R%fromAxisAngle([direction,60.0_pREAL],degrees=.true.,P=1)
              np = R%rotate(normal)
              SchmidMatrix(1:3,1:3,i) = SchmidMatrix(1:3,1:3,i) &
                                      + coeff(1) * math_outer(direction, np) &
                                      + coeff(2) * math_outer(math_cross(normal, direction), normal) &
                                      + coeff(3) * math_outer(math_cross(np, direction), np) &
                                      + coeff(4) * math_outer(normal, normal) &
                                      + coeff(5) * math_outer(math_cross(normal, direction), &
                                                              math_cross(normal, direction)) &
                                      + coeff(6) * math_outer(direction, direction)
          end select
      end select
    end if
  end do

end function crystal_SchmidMatrix_slip


!--------------------------------------------------------------------------------------------------
!> @brief Schmid matrix for twinning
!> details only active twin systems are considered
!--------------------------------------------------------------------------------------------------
function crystal_SchmidMatrix_twin(Ntwin,lattice,cOverA) result(SchmidMatrix)

  integer,     dimension(:),              intent(in) :: Ntwin                                       !< number of active twin systems per family
  character(len=*),                       intent(in) :: lattice                                     !< Bravais lattice (Pearson symbol)
  real(pREAL),                            intent(in) :: cOverA                                      !< c/a ratio
  real(pREAL), dimension(3,3,sum(Ntwin))             :: SchmidMatrix

  real(pREAL), dimension(3,3,sum(Ntwin))             :: coordinateSystem
  real(pREAL), dimension(:,:),           allocatable :: twinSystems
  integer,     dimension(:),             allocatable :: NtwinMax
  integer                                            :: i

  select case(lattice)
    case('cF')
      NtwinMax    = CF_NTWINSYSTEM
      twinSystems = CF_SYSTEMTWIN
    case('cI')
      NtwinMax    = CI_NTWINSYSTEM
      twinSystems = CI_SYSTEMTWIN
    case('hP')
      NtwinMax    = HP_NTWINSYSTEM
      twinSystems = HP_SYSTEMTWIN
    case default
      allocate(NtwinMax(0))
      call IO_error(137,ext_msg='crystal_SchmidMatrix_twin: '//trim(lattice))
  end select

  if (any(NtwinMax(1:size(Ntwin)) - Ntwin < 0)) &
    call IO_error(145,ext_msg='Ntwin '//trim(lattice))
  if (any(Ntwin < 0)) &
    call IO_error(144,ext_msg='Ntwin '//trim(lattice))

  coordinateSystem = buildCoordinateSystem(Ntwin,NtwinMax,twinSystems,lattice,cOverA)

  do i = 1, sum(Ntwin)
    SchmidMatrix(1:3,1:3,i) = math_outer(coordinateSystem(1:3,1,i),coordinateSystem(1:3,2,i))
    if (abs(math_trace33(SchmidMatrix(1:3,1:3,i))) > tol_math_check) &
      error stop 'dilatational Schmid matrix for twin'
  end do

end function crystal_SchmidMatrix_twin


!--------------------------------------------------------------------------------------------------
!> @brief Schmid matrix for transformation
!> details only active twin systems are considered
!--------------------------------------------------------------------------------------------------
function crystal_SchmidMatrix_trans(Ntrans,crystal_target,cOverA,a_cF,a_cI) result(SchmidMatrix)

  integer,         dimension(:),             intent(in) :: Ntrans                                   !< number of active twin systems per family
  character(len=*),                          intent(in) :: crystal_target                           !< Bravais lattice (Pearson symbol)
  real(pREAL),                     optional, intent(in) :: cOverA, a_cI, a_cF
  real(pREAL),     dimension(3,3,sum(Ntrans))           :: SchmidMatrix

  real(pREAL), dimension(3,3,sum(Ntrans)) :: devNull


  if (crystal_target == 'hP' .and. present(cOverA)) then
    if (cOverA < 1.0_pREAL .or. cOverA > 2.0_pREAL) &
    call IO_error(131,ext_msg='crystal_SchmidMatrix_trans: '//trim(crystal_target))
    call buildTransformationSystem(devNull,SchmidMatrix,Ntrans,cOverA=cOverA)
  else if (crystal_target == 'cI' .and. present(a_cF) .and. present(a_cI)) then
    if (a_cI <= 0.0_pREAL .or. a_cF <= 0.0_pREAL) &
    call IO_error(134,ext_msg='crystal_SchmidMatrix_trans: '//trim(crystal_target))
    call buildTransformationSystem(devNull,SchmidMatrix,Ntrans,a_cF=a_cF,a_cI=a_cI)
  else
    call IO_error(131,ext_msg='crystal_SchmidMatrix_trans: '//trim(crystal_target))
  end if

end function crystal_SchmidMatrix_trans


!--------------------------------------------------------------------------------------------------
!> @brief Schmid matrix for cleavage
!> details only active cleavage systems are considered
!--------------------------------------------------------------------------------------------------
function crystal_SchmidMatrix_cleavage(Ncleavage,lattice,cOverA) result(SchmidMatrix)

  integer,     dimension(:),                  intent(in) :: Ncleavage                               !< number of active cleavage systems per family
  character(len=*),                           intent(in) :: lattice                                 !< Bravais lattice (Pearson symbol)
  real(pREAL),                                intent(in) :: cOverA                                  !< c/a ratio
  real(pREAL), dimension(3,3,3,sum(Ncleavage))           :: SchmidMatrix

  real(pREAL), dimension(3,3,sum(Ncleavage))             :: coordinateSystem
  real(pREAL), dimension(:,:),               allocatable :: cleavageSystems
  integer,     dimension(:),                 allocatable :: NcleavageMax
  integer                                                :: i

  select case(lattice)
    case('cF')
      NcleavageMax    = CF_NCLEAVAGESYSTEM
      cleavageSystems = CF_SYSTEMCLEAVAGE
    case('cI')
      NcleavageMax    = CI_NCLEAVAGESYSTEM
      cleavageSystems = CI_SYSTEMCLEAVAGE
    case default
      allocate(NcleavageMax(0))
      call IO_error(137,ext_msg='crystal_SchmidMatrix_cleavage: '//trim(lattice))
  end select

  if (any(NcleavageMax(1:size(Ncleavage)) - Ncleavage < 0)) &
    call IO_error(145,ext_msg='Ncleavage '//trim(lattice))
  if (any(Ncleavage < 0)) &
    call IO_error(144,ext_msg='Ncleavage '//trim(lattice))

  coordinateSystem = buildCoordinateSystem(Ncleavage,NcleavageMax,cleavageSystems,lattice,cOverA)

  do i = 1, sum(Ncleavage)
    SchmidMatrix(1:3,1:3,1,i) = math_outer(coordinateSystem(1:3,1,i),coordinateSystem(1:3,2,i))
    SchmidMatrix(1:3,1:3,2,i) = math_outer(coordinateSystem(1:3,3,i),coordinateSystem(1:3,2,i))
    SchmidMatrix(1:3,1:3,3,i) = math_outer(coordinateSystem(1:3,2,i),coordinateSystem(1:3,2,i))
  end do

end function crystal_SchmidMatrix_cleavage


!--------------------------------------------------------------------------------------------------
!> @brief Slip direction of slip systems (|| b)
!--------------------------------------------------------------------------------------------------
function crystal_slip_direction(Nslip,lattice,cOverA) result(d)

  integer,         dimension(:),           intent(in) :: Nslip                                      !< number of active slip systems per family
  character(len=*),                        intent(in) :: lattice                                    !< Bravais lattice (Pearson symbol)
  real(pREAL),                             intent(in) :: cOverA                                     !< c/a ratio
  real(pREAL),     dimension(3,sum(Nslip))            :: d

  real(pREAL), dimension(3,3,sum(Nslip)) :: coordinateSystem

  coordinateSystem = coordinateSystem_slip(Nslip,lattice,cOverA)
  d = coordinateSystem(1:3,1,1:sum(Nslip))

end function crystal_slip_direction


!--------------------------------------------------------------------------------------------------
!> @brief Normal direction of slip systems (|| n)
!--------------------------------------------------------------------------------------------------
function crystal_slip_normal(Nslip,lattice,cOverA) result(n)

  integer,         dimension(:),           intent(in) :: Nslip                                      !< number of active slip systems per family
  character(len=*),                        intent(in) :: lattice                                    !< Bravais lattice (Pearson symbol)
  real(pREAL),                             intent(in) :: cOverA                                     !< c/a ratio
  real(pREAL),     dimension(3,sum(Nslip))            :: n

  real(pREAL), dimension(3,3,sum(Nslip)) :: coordinateSystem

  coordinateSystem = coordinateSystem_slip(Nslip,lattice,cOverA)
  n = coordinateSystem(1:3,2,1:sum(Nslip))

end function crystal_slip_normal


!--------------------------------------------------------------------------------------------------
!> @brief Transverse direction of slip systems (|| t = b x n)
!--------------------------------------------------------------------------------------------------
function crystal_slip_transverse(Nslip,lattice,cOverA) result(t)

  integer,         dimension(:),           intent(in) :: Nslip                                      !< number of active slip systems per family
  character(len=*),                        intent(in) :: lattice                                    !< Bravais lattice (Pearson symbol)
  real(pREAL),                             intent(in) :: cOverA                                     !< c/a ratio
  real(pREAL),     dimension(3,sum(Nslip))            :: t

  real(pREAL), dimension(3,3,sum(Nslip)) :: coordinateSystem

  coordinateSystem = coordinateSystem_slip(Nslip,lattice,cOverA)
  t = coordinateSystem(1:3,3,1:sum(Nslip))

end function crystal_slip_transverse


!--------------------------------------------------------------------------------------------------
!> @brief Labels of slip systems
!> details only active slip systems are considered
!--------------------------------------------------------------------------------------------------
function crystal_labels_slip(Nslip,lattice) result(labels)

  integer,         dimension(:),  intent(in)  :: Nslip                                              !< number of active slip systems per family
  character(len=*),               intent(in)  :: lattice                                            !< Bravais lattice (Pearson symbol)

  character(len=:), dimension(:),   allocatable :: labels

  real(pREAL),      dimension(:,:), allocatable :: slipSystems
  integer,          dimension(:),   allocatable :: NslipMax

  select case(lattice)
    case('cF')
      NslipMax    = CF_NSLIPSYSTEM
      slipSystems = CF_SYSTEMSLIP
    case('cI')
      NslipMax    = CI_NSLIPSYSTEM
      slipSystems = CI_SYSTEMSLIP
    case('hP')
      NslipMax    = HP_NSLIPSYSTEM
      slipSystems = HP_SYSTEMSLIP
    case('tI')
      NslipMax    = TI_NSLIPSYSTEM
      slipSystems = TI_SYSTEMSLIP
    case default
      call IO_error(137,ext_msg='crystal_labels_slip: '//trim(lattice))
  end select

  if (any(NslipMax(1:size(Nslip)) - Nslip < 0)) &
    call IO_error(145,ext_msg='Nslip '//trim(lattice))
  if (any(Nslip < 0)) &
    call IO_error(144,ext_msg='Nslip '//trim(lattice))

  labels = getLabels(Nslip,NslipMax,slipSystems)

end function crystal_labels_slip


!--------------------------------------------------------------------------------------------------
!> @brief Return 3x3 tensor with symmetry according to given Bravais lattice
!--------------------------------------------------------------------------------------------------
pure function crystal_symmetrize_33(T,lattice) result(T_sym)

  real(pREAL), dimension(3,3) :: T_sym

  real(pREAL), dimension(3,3), intent(in) :: T
  character(len=*),            intent(in) :: lattice                                                !< Bravais lattice (Pearson symbol)


  T_sym = 0.0_pREAL

  select case(lattice)
    case('cF','cI')
      T_sym(1,1) = T(1,1)
      T_sym(2,2) = T(1,1)
      T_sym(3,3) = T(1,1)
    case('hP','tI')
      T_sym(1,1) = T(1,1)
      T_sym(2,2) = T(1,1)
      T_sym(3,3) = T(3,3)
   end select

end function crystal_symmetrize_33


!--------------------------------------------------------------------------------------------------
!> @brief Return stiffness matrix in 6x6 notation with symmetry according to given Bravais lattice
!> @details J. A. Rayne and B. S. Chandrasekhar Phys. Rev. 120, 1658 Erratum Phys. Rev. 122, 1962
!--------------------------------------------------------------------------------------------------
pure function crystal_symmetrize_C66(C66,lattice) result(C66_sym)

  real(pREAL), dimension(6,6) :: C66_sym

  real(pREAL), dimension(6,6), intent(in) :: C66
  character(len=*),            intent(in) :: lattice                                                !< Bravais lattice (Pearson symbol)

  integer :: i,j


  C66_sym = 0.0_pREAL

  select case(lattice)
    case ('cF','cI')
      C66_sym(1,1) = C66(1,1); C66_sym(2,2) = C66(1,1); C66_sym(3,3) = C66(1,1)
      C66_sym(1,2) = C66(1,2); C66_sym(1,3) = C66(1,2); C66_sym(2,3) = C66(1,2)
      C66_sym(4,4) = C66(4,4); C66_sym(5,5) = C66(4,4); C66_sym(6,6) = C66(4,4)                     ! isotropic C_44 = (C_11-C_12)/2
    case ('hP')
      C66_sym(1,1) = C66(1,1); C66_sym(2,2) = C66(1,1)
      C66_sym(3,3) = C66(3,3)
      C66_sym(1,2) = C66(1,2)
      C66_sym(1,3) = C66(1,3); C66_sym(2,3) = C66(1,3)
      C66_sym(4,4) = C66(4,4); C66_sym(5,5) = C66(4,4)
      C66_sym(6,6) = 0.5_pREAL*(C66(1,1)-C66(1,2))
    case ('tI')
      C66_sym(1,1) = C66(1,1); C66_sym(2,2) = C66(1,1)
      C66_sym(3,3) = C66(3,3)
      C66_sym(1,2) = C66(1,2)
      C66_sym(1,3) = C66(1,3); C66_sym(2,3) = C66(1,3)
      C66_sym(4,4) = C66(4,4); C66_sym(5,5) = C66(4,4)
      C66_sym(6,6) = C66(6,6)
   end select

   do i = 1, 6
     do j = i+1, 6
       C66_sym(j,i) = C66_sym(i,j)
     end do
   end do

end function crystal_symmetrize_C66


!--------------------------------------------------------------------------------------------------
!> @brief Labels for twin systems
!> details only active twin systems are considered
!--------------------------------------------------------------------------------------------------
function crystal_labels_twin(Ntwin,lattice) result(labels)

  integer,         dimension(:),  intent(in)  :: Ntwin                                              !< number of active slip systems per family
  character(len=*),               intent(in)  :: lattice                                            !< Bravais lattice (Pearson symbol)

  character(len=:), dimension(:),   allocatable :: labels

  real(pREAL),      dimension(:,:), allocatable :: twinSystems
  integer,          dimension(:),   allocatable :: NtwinMax

  select case(lattice)
    case('cF')
      NtwinMax    = CF_NTWINSYSTEM
      twinSystems = CF_SYSTEMTWIN
    case('cI')
      NtwinMax    = CI_NTWINSYSTEM
      twinSystems = CI_SYSTEMTWIN
    case('hP')
      NtwinMax    = HP_NTWINSYSTEM
      twinSystems = HP_SYSTEMTWIN
    case default
      call IO_error(137,ext_msg='crystal_labels_twin: '//trim(lattice))
  end select

  if (any(NtwinMax(1:size(Ntwin)) - Ntwin < 0)) &
    call IO_error(145,ext_msg='Ntwin '//trim(lattice))
  if (any(Ntwin < 0)) &
    call IO_error(144,ext_msg='Ntwin '//trim(lattice))

  labels = getLabels(Ntwin,NtwinMax,twinSystems)

end function crystal_labels_twin


!--------------------------------------------------------------------------------------------------
!> @brief Projection of the transverse direction onto the slip plane
!> @details: This projection is used to calculate forest hardening for edge dislocations
!--------------------------------------------------------------------------------------------------
function slipProjection_transverse(Nslip,lattice,cOverA) result(projection)

  integer,         dimension(:),                   intent(in) :: Nslip                              !< number of active slip systems per family
  character(len=*),                                intent(in) :: lattice                            !< Bravais lattice (Pearson symbol)
  real(pREAL),                                     intent(in) :: cOverA                             !< c/a ratio
  real(pREAL),     dimension(sum(Nslip),sum(Nslip))           :: projection

  real(pREAL), dimension(3,sum(Nslip)) :: n, t
  integer                              :: i, j

  n = crystal_slip_normal    (Nslip,lattice,cOverA)
  t = crystal_slip_transverse(Nslip,lattice,cOverA)

  do i=1, sum(Nslip); do j=1, sum(Nslip)
    projection(i,j) = abs(math_inner(n(:,i),t(:,j)))
  end do; end do

end function slipProjection_transverse


!--------------------------------------------------------------------------------------------------
!> @brief Projection of the slip direction onto the slip plane
!> @details: This projection is used to calculate forest hardening for screw dislocations
!--------------------------------------------------------------------------------------------------
function slipProjection_direction(Nslip,lattice,cOverA) result(projection)

  integer,         dimension(:),                   intent(in) :: Nslip                              !< number of active slip systems per family
  character(len=*),                                intent(in) :: lattice                            !< Bravais lattice (Pearson symbol)
  real(pREAL),                                     intent(in) :: cOverA                             !< c/a ratio
  real(pREAL),     dimension(sum(Nslip),sum(Nslip))           :: projection

  real(pREAL), dimension(3,sum(Nslip)) :: n, d
  integer                              :: i, j

  n = crystal_slip_normal   (Nslip,lattice,cOverA)
  d = crystal_slip_direction(Nslip,lattice,cOverA)

  do i=1, sum(Nslip); do j=1, sum(Nslip)
    projection(i,j) = abs(math_inner(n(:,i),d(:,j)))
  end do; end do

end function slipProjection_direction


!--------------------------------------------------------------------------------------------------
!> @brief build a local coordinate system on slip systems
!> @details Order: Direction, plane (normal), and common perpendicular
!--------------------------------------------------------------------------------------------------
function coordinateSystem_slip(Nslip,lattice,cOverA) result(coordinateSystem)

  integer,          dimension(:),            intent(in) :: Nslip                                    !< number of active slip systems per family
  character(len=*),                          intent(in) :: lattice                                  !< Bravais lattice (Pearson symbol)
  real(pREAL),                               intent(in) :: cOverA                                   !< c/a ratio
  real(pREAL),     dimension(3,3,sum(Nslip))            :: coordinateSystem

  real(pREAL), dimension(:,:), allocatable :: slipSystems
  integer,     dimension(:),   allocatable :: NslipMax

  select case(lattice)
    case('cF')
      NslipMax    = CF_NSLIPSYSTEM
      slipSystems = CF_SYSTEMSLIP
    case('cI')
      NslipMax    = CI_NSLIPSYSTEM
      slipSystems = CI_SYSTEMSLIP
    case('hP')
      NslipMax    = HP_NSLIPSYSTEM
      slipSystems = HP_SYSTEMSLIP
    case('tI')
      NslipMax    = TI_NSLIPSYSTEM
      slipSystems = TI_SYSTEMSLIP
    case default
      allocate(NslipMax(0))
      call IO_error(137,ext_msg='coordinateSystem_slip: '//trim(lattice))
  end select

  if (any(NslipMax(1:size(Nslip)) - Nslip < 0)) &
    call IO_error(145,ext_msg='Nslip '//trim(lattice))
  if (any(Nslip < 0)) &
    call IO_error(144,ext_msg='Nslip '//trim(lattice))

  coordinateSystem = buildCoordinateSystem(Nslip,NslipMax,slipSystems,lattice,cOverA)

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
  real(pREAL), dimension(:),                                 intent(in) :: values                   !< interaction values
  integer,     dimension(:,:),                               intent(in) :: matrix                   !< interaction types
  real(pREAL), dimension(sum(reacting_used),sum(acting_used))           :: buildInteraction

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

      end do; end do
  end do; end do

end function buildInteraction


!--------------------------------------------------------------------------------------------------
!> @brief Build a local coordinate system on slip, twin, trans, cleavage systems
!> @details Order: Direction, plane (normal), and common perpendicular
!--------------------------------------------------------------------------------------------------
function buildCoordinateSystem(active,potential,system,lattice,cOverA) result(coordinateSystem)

  integer,     dimension(:),   intent(in) :: &
    active, &                                                                                       !< # of active systems per family
    potential                                                                                       !< # of potential systems per family
  real(pREAL), dimension(:,:), intent(in) :: &
    system
  character(len=*),            intent(in) :: &
    lattice                                                                                         !< Bravais lattice (Pearson symbol)
  real(pREAL),                 intent(in) :: &
    cOverA
  real(pREAL), dimension(3,3,sum(active)) :: &
    coordinateSystem

  real(pREAL), dimension(3) :: &
    direction, normal
  integer :: &
    a, &                                                                                            !< index of active system
    p, &                                                                                            !< index in potential system matrix
    f, &                                                                                            !< index of my family
    s                                                                                               !< index of my system in current family

  if (lattice == 'tI' .and. cOverA > 2.0_pREAL) &
    call IO_error(131,ext_msg='buildCoordinateSystem:'//trim(lattice))
  if (lattice == 'hP' .and. (cOverA < 1.0_pREAL .or. cOverA > 2.0_pREAL)) &
    call IO_error(131,ext_msg='buildCoordinateSystem:'//trim(lattice))

  a = 0
  activeFamilies: do f = 1,size(active,1)
    activeSystems: do s = 1,active(f)
      a = a + 1
      p = sum(potential(1:f-1))+s

      select case(lattice)

        case ('cF','cI')
          direction = system(1:3,p)
          normal    = system(4:6,p)

        case ('tI')
          direction = [ system(1,p), system(2,p), system(3,p)*cOverA ]
          normal    = [ system(4,p), system(5,p), system(6,p)/cOverA ]

        case ('hP')
          direction = [ system(1,p)*1.5_pREAL, &
                       (system(1,p)+2.0_pREAL*system(2,p))*sqrt(0.75_pREAL), &
                        system(4,p)*cOverA ]                                                        ! direction [uvtw]->[3u/2 (u+2v)*sqrt(3)/2 w*(p/a)])
          normal    = [ system(5,p), &
                       (system(5,p)+2.0_pREAL*system(6,p))/sqrt(3.0_pREAL), &
                        system(8,p)/cOverA ]                                                        ! plane (hkil)->(h (h+2k)/sqrt(3) l/(p/a))

        case default
          call IO_error(137,ext_msg='buildCoordinateSystem: '//trim(lattice))

      end select

      coordinateSystem(1:3,1,a) = direction/norm2(direction)
      coordinateSystem(1:3,2,a) = normal   /norm2(normal)
      coordinateSystem(1:3,3,a) = math_cross(direction/norm2(direction),&
                                             normal   /norm2(normal))

    end do activeSystems
  end do activeFamilies

end function buildCoordinateSystem


!--------------------------------------------------------------------------------------------------
!> @brief Helper function to define transformation systems
! Needed to calculate Schmid matrix and rotated stiffness matrices.
! @details: use c/a for cF -> hP transformation
!           use a_cX for cF -> cI transformation
!--------------------------------------------------------------------------------------------------
subroutine buildTransformationSystem(Q,S,Ntrans,cOverA,a_cF,a_cI)

  integer,      dimension(:),               intent(in) :: &
    Ntrans
  real(pREAL),  dimension(3,3,sum(Ntrans)), intent(out) :: &
    Q, &                                                                                            !< Total rotation: Q = R*B
    S                                                                                               !< Eigendeformation tensor for phase transformation
  real(pREAL),  optional,                   intent(in) :: &
    cOverA, &                                                                                       !< c/a for target hP lattice
    a_cF, &                                                                                         !< lattice parameter a for cF target lattice
    a_cI                                                                                            !< lattice parameter a for cI parent lattice

  type(tRotation) :: &
    R, &                                                                                            !< Pitsch rotation
    B                                                                                               !< Rotation of cF to Bain coordinate system
  real(pREAL), dimension(3,3) :: &
    U, &                                                                                            !< Bain deformation
    ss, sd
  real(pREAL), dimension(3) :: &
    x, y, z
  integer :: &
    i
  real(pREAL), dimension(3+3,CF_NTRANS), parameter :: &
    CFTOHP_SYSTEMTRANS = reshape(real( [&
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
      ],pREAL),shape(CFTOHP_SYSTEMTRANS))

  real(pREAL), dimension(4,cF_Ntrans), parameter :: &
    CFTOCI_SYSTEMTRANS = real(reshape([&
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
      ],shape(CFTOCI_SYSTEMTRANS)),pREAL)

  integer, dimension(9,cF_Ntrans), parameter :: &
    CFTOCI_BAINVARIANT = reshape( [&
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
      ],shape(CFTOCI_BAINVARIANT))

  real(pREAL), dimension(4,cF_Ntrans), parameter :: &
    CFTOCI_BAINROT = real(reshape([&
      1.0, 0.0, 0.0,     45.0, &                                                                    ! Rotate cF austensite to bain variant
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
      ],shape(CFTOCI_BAINROT)),pREAL)

  if (present(a_cI) .and. present(a_cF)) then
    do i = 1,sum(Ntrans)
      call R%fromAxisAngle(CFTOCI_SYSTEMTRANS(:,i),degrees=.true.,P=1)
      call B%fromAxisAngle(CFTOCI_BAINROT(:,i),    degrees=.true.,P=1)
      x = real(CFTOCI_BAINVARIANT(1:3,i),pREAL)
      y = real(CFTOCI_BAINVARIANT(4:6,i),pREAL)
      z = real(CFTOCI_BAINVARIANT(7:9,i),pREAL)

      U = (a_cI/a_cF) * (math_outer(x,x) + (math_outer(y,y)+math_outer(z,z)) * sqrt(2.0_pREAL))
      Q(1:3,1:3,i) = matmul(R%asMatrix(),B%asMatrix())
      S(1:3,1:3,i) = matmul(R%asMatrix(),U) - MATH_I3
    end do
  else if (present(cOverA)) then
    ss      = MATH_I3
    sd      = MATH_I3
    ss(1,3) = sqrt(2.0_pREAL)/4.0_pREAL
    sd(3,3) = cOverA/sqrt(8.0_pREAL/3.0_pREAL)

    do i = 1,sum(Ntrans)
      x = CFTOHP_SYSTEMTRANS(1:3,i)/norm2(CFTOHP_SYSTEMTRANS(1:3,i))
      z = CFTOHP_SYSTEMTRANS(4:6,i)/norm2(CFTOHP_SYSTEMTRANS(4:6,i))
      y = -math_cross(x,z)
      Q(1:3,1,i) = x
      Q(1:3,2,i) = y
      Q(1:3,3,i) = z
      S(1:3,1:3,i) = matmul(Q(1:3,1:3,i), matmul(matmul(sd,ss), transpose(Q(1:3,1:3,i)))) - MATH_I3 ! ToDo: This is of interest for the Schmid matrix only
    end do
  else
    call IO_error(132,ext_msg='buildTransformationSystem')
  end if

end subroutine buildTransformationSystem


!--------------------------------------------------------------------------------------------------
!> @brief select active systems as strings
!--------------------------------------------------------------------------------------------------
function getlabels(active,potential,system) result(labels)

  integer,          dimension(:),   intent(in) :: &
    active, &                                                                                       !< # of active systems per family
    potential                                                                                       !< # of potential systems per family
  real(pREAL),      dimension(:,:), intent(in) :: &
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
      end do direction
      label(i:i) = ']'

      i = i +1
      label(i:i) = '('
      normal: do j = size(system,1)/2+1, size(system,1)
        write(label(i+1:i+2),'(I2.1)') int(system(j,p))
        label(i+3:i+3) = ' '
        i = i + 3
      end do normal
      label(i:i) = ')'

      labels(a) = label

    end do activeSystems
  end do activeFamilies

end function getlabels


!--------------------------------------------------------------------------------------------------
!> @brief Equivalent Poisson's ratio (ν)
!> @details https://doi.org/10.1143/JPSJ.20.635
!--------------------------------------------------------------------------------------------------
pure function crystal_isotropic_nu(C,assumption,lattice) result(nu)

  real(pREAL), dimension(6,6), intent(in) :: C                                                      !< Stiffness tensor (Voigt notation)
  character(len=*),            intent(in) :: assumption                                             !< Assumption (isostrain = 'Voigt', isostress = 'Reuss')
  character(len=*), optional,  intent(in) :: lattice
  real(pREAL) :: nu

  real(pREAL) :: K, mu
  logical                       :: error
  real(pREAL), dimension(6,6)   :: S


  if     (IO_lc(assumption) == 'isostrain') then
    K = sum(C(1:3,1:3)) / 9.0_pREAL
  elseif (IO_lc(assumption) == 'isostress') then
    call math_invert(S,error,C)
    if (error) error stop 'matrix inversion failed'
    K = 1.0_pREAL / sum(S(1:3,1:3))
  else
    error stop 'invalid assumption'
  end if

  mu = crystal_isotropic_mu(C,assumption,lattice)
  nu = (1.5_pREAL*K-mu)/(3.0_pREAL*K+mu)

end function crystal_isotropic_nu


!--------------------------------------------------------------------------------------------------
!> @brief Equivalent shear modulus (μ)
!> @details https://doi.org/10.1143/JPSJ.20.635
!> @details Nonlinear Mechanics of Crystals 10.1007/978-94-007-0350-6, pp 563
!--------------------------------------------------------------------------------------------------
pure function crystal_isotropic_mu(C,assumption,lattice) result(mu)

  real(pREAL), dimension(6,6), intent(in) :: C                                                      !< Stiffness tensor (Voigt notation)
  character(len=*),            intent(in) :: assumption                                             !< Assumption (isostrain = 'Voigt', isostress = 'Reuss')
  character(len=*), optional,  intent(in) :: lattice
  real(pREAL) :: mu

  logical                       :: error
  real(pREAL), dimension(6,6)   :: S


  if     (IO_lc(assumption) == 'isostrain') then
      select case(misc_optional(lattice,''))
        case('cF','cI')
          mu = ( C(1,1) - C(1,2) + C(4,4)*3.0_pREAL) / 5.0_pREAL
        case default
          mu = (  C(1,1)+C(2,2)+C(3,3) &
                - C(1,2)-C(2,3)-C(1,3) &
                +(C(4,4)+C(5,5)+C(6,6)) * 3.0_pREAL &
               ) / 15.0_pREAL
      end select

  elseif (IO_lc(assumption) == 'isostress') then
      select case(misc_optional(lattice,''))
        case('cF','cI')
          mu = 5.0_pREAL &
               / (4.0_pREAL/(C(1,1)-C(1,2)) + 3.0_pREAL/C(4,4))
        case default
          call math_invert(S,error,C)
          if (error) error stop 'matrix inversion failed'
          mu = 15.0_pREAL &
              / (4.0_pREAL*(S(1,1)+S(2,2)+S(3,3)-S(1,2)-S(2,3)-S(1,3)) + 3.0_pREAL*(S(4,4)+S(5,5)+S(6,6)))
      end select
  else
    error stop 'invalid assumption'
  end if

end function crystal_isotropic_mu


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some crystal functions.
!--------------------------------------------------------------------------------------------------
subroutine crystal_selfTest()

  real(pREAL), dimension(:,:,:), allocatable :: CoSy
  real(pREAL), dimension(:,:),   allocatable :: system

  real(pREAL), dimension(6,6) :: C, C_cF, C_cI, C_hP, C_tI
  real(pREAL), dimension(3,3) :: T, T_cF, T_cI, T_hP, T_tI
  real(pREAL), dimension(2)   :: r
  real(pREAL) :: lambda
  integer     :: i


  call random_number(r)

  system = reshape([1.0_pREAL+r(1),0.0_pREAL,0.0_pREAL, 0.0_pREAL,1.0_pREAL+r(2),0.0_pREAL],[6,1])
  CoSy   = buildCoordinateSystem([1],[1],system,'cF',0.0_pREAL)
  if (any(dNeq(CoSy(1:3,1:3,1),math_I3))) error stop 'buildCoordinateSystem'

  if (any(dNeq(buildCoordinateSystem(TI_NSLIPSYSTEM,TI_NSLIPSYSTEM,TI_SYSTEMSLIP,'cI',0.0_pReal), &
               buildCoordinateSystem(TI_NSLIPSYSTEM,TI_NSLIPSYSTEM,TI_SYSTEMSLIP,'tI',1.0_pReal)))) &
                                                           error stop 'cI/tI coordinate system'
  if (all(dEq( buildCoordinateSystem(TI_NSLIPSYSTEM,TI_NSLIPSYSTEM,TI_SYSTEMSLIP,'tI',1.1_pReal + r(1)*0.9_pReal), &
               buildCoordinateSystem(TI_NSLIPSYSTEM,TI_NSLIPSYSTEM,TI_SYSTEMSLIP,'tI',1.0_pReal)))) &
                                                           error stop 'tI coordinate system'
  do i = 1, 10
    call random_number(C)
    C_cF = crystal_symmetrize_C66(C,'cI')
    C_cI = crystal_symmetrize_C66(C,'cF')
    C_hP = crystal_symmetrize_C66(C,'hP')
    C_tI = crystal_symmetrize_C66(C,'tI')

    if (any(dNeq(C_cI,transpose(C_cF))))                   error stop 'SymmetryC66/cI-cF'
    if (any(dNeq(C_cF,transpose(C_cI))))                   error stop 'SymmetryC66/cF-cI'
    if (any(dNeq(C_hP,transpose(C_hP))))                   error stop 'SymmetryC66/hP'
    if (any(dNeq(C_tI,transpose(C_tI))))                   error stop 'SymmetryC66/tI'

    if (any(dNeq(C(1,1),[C_cF(1,1),C_cF(2,2),C_cF(3,3)]))) error stop 'SymmetryC_11-22-33/c'
    if (any(dNeq(C(1,2),[C_cF(1,2),C_cF(1,3),C_cF(2,3)]))) error stop 'SymmetryC_12-13-23/c'
    if (any(dNeq(C(4,4),[C_cF(4,4),C_cF(5,5),C_cF(6,6)]))) error stop 'SymmetryC_44-55-66/c'

    if (any(dNeq(C(1,1),[C_hP(1,1),C_hP(2,2)])))           error stop 'SymmetryC_11-22/hP'
    if (any(dNeq(C(1,3),[C_hP(1,3),C_hP(2,3)])))           error stop 'SymmetryC_13-23/hP'
    if (any(dNeq(C(4,4),[C_hP(4,4),C_hP(5,5)])))           error stop 'SymmetryC_44-55/hP'

    if (any(dNeq(C(1,1),[C_tI(1,1),C_tI(2,2)])))           error stop 'SymmetryC_11-22/tI'
    if (any(dNeq(C(1,3),[C_tI(1,3),C_tI(2,3)])))           error stop 'SymmetryC_13-23/tI'
    if (any(dNeq(C(4,4),[C_tI(4,4),C_tI(5,5)])))           error stop 'SymmetryC_44-55/tI'

    call random_number(T)
    T_cF = crystal_symmetrize_33(T,'cI')
    T_cI = crystal_symmetrize_33(T,'cF')
    T_hP = crystal_symmetrize_33(T,'hP')
    T_tI = crystal_symmetrize_33(T,'tI')

    if (any(dNeq0(T_cF) .and. math_I3<1.0_pREAL))          error stop 'Symmetry33/c'
    if (any(dNeq0(T_hP) .and. math_I3<1.0_pREAL))          error stop 'Symmetry33/hP'
    if (any(dNeq0(T_tI) .and. math_I3<1.0_pREAL))          error stop 'Symmetry33/tI'

    if (any(dNeq(T(1,1),[T_cI(1,1),T_cI(2,2),T_cI(3,3)]))) error stop 'Symmetry33_11-22-33/c'
    if (any(dNeq(T(1,1),[T_hP(1,1),T_hP(2,2)])))           error stop 'Symmetry33_11-22/hP'
    if (any(dNeq(T(1,1),[T_tI(1,1),T_tI(2,2)])))           error stop 'Symmetry33_11-22/tI'

  end do

  call random_number(C)
  C(1,1) = C(1,1) + C(1,2) + 0.1_pREAL
  C(1,3) = C(1,2)
  C(3,3) = C(1,1)
  C(4,4) = 0.5_pREAL * (C(1,1) - C(1,2))
  C(6,6) = C(4,4)

  C_cI = crystal_symmetrize_C66(C,'cI')
  if (dNeq(C_cI(4,4),crystal_isotropic_mu(C_cI,'isostrain','cI'),1.0e-12_pREAL)) error stop 'isotropic_mu/isostrain/cI'
  if (dNeq(C_cI(4,4),crystal_isotropic_mu(C_cI,'isostress','cI'),1.0e-12_pREAL)) error stop 'isotropic_mu/isostress/cI'

  lambda = C_cI(1,2)
  if (dNeq(lambda*0.5_pREAL/(lambda+crystal_isotropic_mu(C_cI,'isostrain','cI')), &
          crystal_isotropic_nu(C_cI,'isostrain','cI'),1.0e-12_pREAL)) error stop 'isotropic_nu/isostrain/cI'
  if (dNeq(lambda*0.5_pREAL/(lambda+crystal_isotropic_mu(C_cI,'isostress','cI')), &
          crystal_isotropic_nu(C_cI,'isostress','cI'),1.0e-12_pREAL)) error stop 'isotropic_nu/isostress/cI'


  C_hP = crystal_symmetrize_C66(C,'hP')
  if (dNeq(C(4,4),crystal_isotropic_mu(C_hP,'isostrain','hP'),1.0e-12_pREAL)) error stop 'isotropic_mu/isostrain/hP'
  if (dNeq(C(4,4),crystal_isotropic_mu(C_hP,'isostress','hP'),1.0e-12_pREAL)) error stop 'isotropic_mu/isostress/hP'

  lambda = C_hP(1,2)
  if (dNeq(lambda*0.5_pREAL/(lambda+crystal_isotropic_mu(C_hP,'isostrain','hP')), &
          crystal_isotropic_nu(C_hP,'isostrain','hP'),1.0e-12_pREAL)) error stop 'isotropic_nu/isostrain/hP'
  if (dNeq(lambda*0.5_pREAL/(lambda+crystal_isotropic_mu(C_hP,'isostress','hP')), &
          crystal_isotropic_nu(C_hP,'isostress','hP'),1.0e-12_pREAL)) error stop 'isotropic_nu/isostress/hP'

  C_tI = crystal_symmetrize_C66(C,'tI')
  if (dNeq(C(6,6),crystal_isotropic_mu(C_tI,'isostrain','tI'),1.0e-12_pREAL)) error stop 'isotropic_mu/isostrain/tI'
  if (dNeq(C(6,6),crystal_isotropic_mu(C_tI,'isostress','tI'),1.0e-12_pREAL)) error stop 'isotropic_mu/isostress/tI'

  lambda = C_tI(1,2)
  if (dNeq(lambda*0.5_pREAL/(lambda+crystal_isotropic_mu(C_tI,'isostrain','tI')), &
          crystal_isotropic_nu(C_tI,'isostrain','tI'),1.0e-12_pREAL)) error stop 'isotropic_nu/isostrain/tI'
  if (dNeq(lambda*0.5_pREAL/(lambda+crystal_isotropic_mu(C_tI,'isostress','tI')), &
          crystal_isotropic_nu(C_tI,'isostress','tI'),1.0e-12_pREAL)) error stop 'isotropic_nu/isostress/tI'

  call random_number(C)
  C = crystal_symmetrize_C66(C+math_eye(6),'cI')
  if (dNeq(crystal_isotropic_mu(C,'isostrain','cI'), crystal_isotropic_mu(C,'isostrain','hP'), 1.0e-12_pREAL)) &
    error stop 'isotropic_mu/isostrain/cI-hP'
  if (dNeq(crystal_isotropic_nu(C,'isostrain','cF'), crystal_isotropic_nu(C,'isostrain','cI'), 1.0e-12_pREAL)) &
    error stop 'isotropic_nu/isostrain/cF-tI'
  if (dNeq(crystal_isotropic_mu(C,'isostress','cI'), crystal_isotropic_mu(C,'isostress'), 1.0e-12_pREAL)) &
    error stop 'isotropic_mu/isostress/cI-hP'
  if (dNeq(crystal_isotropic_nu(C,'isostress','cF'), crystal_isotropic_nu(C,'isostress'), 1.0e-12_pREAL)) &
    error stop 'isotropic_nu/isostress/cF-tI'

end subroutine crystal_selfTest

end module crystal
