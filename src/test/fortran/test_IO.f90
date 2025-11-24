module test_IO
  use prec
  use parallelization
  use IO

  implicit none(type,external)

  private
  public :: test_IO_run

  contains

subroutine test_IO_run()

  real, dimension(30) :: rnd_real
  character(len=size(rnd_real)) :: rnd_str
  character(len=pSTRLEN), dimension(1) :: strarray_out
  character(len=:), allocatable :: str_out, fname
  integer :: u,i


  call IO_selfTest()

  call random_number(rnd_real)
  fname = 'test'//IO_intAsStr(worldrank)//'.txt'

  do i = 1, size(rnd_real)
    rnd_str(i:i) = char(32 + int(rnd_real(i)*(127.-32.)))
  end do
  open(newunit=u,file=fname,status='replace',form='formatted')
  write(u,'(a)') rnd_str
  close(u)

  str_out = IO_read(fname)
  if (rnd_str//IO_EOL /= str_out) error stop 'IO_read'

end subroutine test_IO_run

end module test_IO
