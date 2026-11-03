! Driver for subs/mat5pen.f - the block-Thomas Schur-complement
! kernel used by the TRI path in cmf_blkband_v3.f.
!
! Reads one case file and prints A_out:
!   line 1: N NR ENOTZERO DNOTZERO
!   then A (N x NR), B (N x N), C (N x NR), D (N x N),
!        E (N x NR), F (N x N), each row-major, one row per line.
!
! Production semantics, both flags:
!   DNOTZERO: WRKMAT = F.C - E ; A = A + D.WRKMAT
!   always:   A = A - B.C
!
! NOTE: MAT5PEN passes the matrices to DGEMM with LDA=N, so the
! driver must allocate exact-size arrays - a padded (NMAX,NMAX)
! array misaligns every column after the first.
!
program drv_mat5pen
  use set_kind_module
  implicit none
  integer :: lu, ios, n, nr, i, j
  logical :: enotzero, dnotzero
  real(kind=ldp), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:)
  real(kind=ldp), allocatable :: e(:,:), f(:,:), wrk(:,:), vec(:)
  character(len=512) :: cfile

  call get_command_argument(1, cfile)
  open(newunit=lu, file=trim(cfile), status='old', action='read', &
       iostat=ios)
  if (ios /= 0) stop 2
  read(lu,*) n, nr, enotzero, dnotzero
  allocate(a(n,nr), b(n,n), c(n,nr), d(n,n), e(n,nr), f(n,n), &
           wrk(n,nr), vec(n))
  do i = 1, n
     read(lu,*) (a(i,j), j=1,nr)
  end do
  do i = 1, n
     read(lu,*) (b(i,j), j=1,n)
  end do
  do i = 1, n
     read(lu,*) (c(i,j), j=1,nr)
  end do
  do i = 1, n
     read(lu,*) (d(i,j), j=1,n)
  end do
  do i = 1, n
     read(lu,*) (e(i,j), j=1,nr)
  end do
  do i = 1, n
     read(lu,*) (f(i,j), j=1,n)
  end do
  close(lu)

  call mat5pen(a,b,c,d,e,f,wrk,vec,n,nr,enotzero,dnotzero)

  write(*,'(a)') 'A_OUT'
  do i = 1, n
     write(*,'(*(es24.16e3,1x))') (a(i,j), j=1,nr)
  end do

end program drv_mat5pen
