!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Fortran interfaces for LAPACK routines
!> @details https://www.netlib.org/lapack/
!--------------------------------------------------------------------------------------------------
module LAPACK_interface
  interface

    pure subroutine dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
      use prec
      implicit none

      character,   intent(in)                             :: jobvl,jobvr
      integer,     intent(in)                             :: n,lda,ldvl,ldvr,lwork
      real(pReal), intent(inout), dimension(lda,n)        :: a
      real(pReal), intent(out),   dimension(n)            :: wr,wi
      real(pReal), intent(out),   dimension(ldvl,n)       :: vl
      real(pReal), intent(out),   dimension(ldvr,n)       :: vr
      real(pReal), intent(out),   dimension(max(1,lwork)) :: work
      integer,     intent(out)                            :: info
    end subroutine dgeev

    pure subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
      use prec
      implicit none

      integer,     intent(in)                             :: n,nrhs,lda,ldb
      real(pReal), intent(inout), dimension(lda,n)        :: a
      integer,     intent(out),   dimension(n)            :: ipiv
      real(pReal), intent(inout), dimension(ldb,nrhs)     :: b
      integer,     intent(out)                            :: info
    end subroutine dgesv

    pure subroutine dgetrf(m,n,a,lda,ipiv,info)
      use prec
      implicit none

      integer,     intent(in)                             :: m,n,lda
      real(pReal), intent(inout), dimension(lda,n)        :: a
      integer,     intent(out),   dimension(min(m,n))     :: ipiv
      integer,     intent(out)                            :: info
    end subroutine dgetrf

    pure subroutine dgetri(n,a,lda,ipiv,work,lwork,info)
      use prec
      implicit none

      integer,     intent(in)                             :: n,lda,lwork
      real(pReal), intent(inout), dimension(lda,n)        :: a
      integer,     intent(in),    dimension(n)            :: ipiv
      real(pReal), intent(out),   dimension(max(1,lwork)) :: work
      integer,     intent(out)                            :: info
    end subroutine dgetri

    pure subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      use prec
      implicit none

      character,   intent(in)                             :: jobz,uplo
      integer,     intent(in)                             :: n,lda,lwork
      real(pReal), intent(inout), dimension(lda,n)        :: a
      real(pReal), intent(out),   dimension(n)            :: w
      real(pReal), intent(out),   dimension(max(1,lwork)) :: work
      integer,     intent(out)                            :: info
    end subroutine dsyev

  end interface

end module LAPACK_interface
