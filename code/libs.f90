!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief dummy source for inclusion of Library files
!--------------------------------------------------------------------------------------------------
module libs
!nothing in here
end module libs

#ifdef Spectral
#include "../lib/kdtree2.f90"
#endif
#ifdef FEM
#include "../lib/kdtree2.f90"
#endif
#include "../lib/IR_Precision.f90"
#include "../lib/Lib_Base64.f90"
#include "../lib/Lib_VTK_IO.f90"

