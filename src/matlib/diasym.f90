!---------------------------------------------------------!
!Calls the LAPACK diagonalization subroutine DSYEV        !
!input:  a(n,n) = real symmetric matrix to be diagonalized!
!            n  = size of a                               !
!output: a(n,n) = orthonormal eigenvectors of a           !
!        eig(n) = eigenvalues of a in ascending order     !
!---------------------------------------------------------!
! This code is part of FIOs-EDOs program, licensed with GPLv3 .
! This is too simple to have copyrights and I don't rememeber
! if it was taken from the manual directly
subroutine diasym(a,eig,n)

  implicit none

  integer ( kind = 4 ) :: n,l,inf
  real    ( kind = 8 ) :: a(n,n),eig(n),work(n*(3+n/2))

  external dsyev

  l=n*(3+n/2)
  call dsyev('V','U',n,a,n,eig,work,l,inf)

end subroutine diasym
