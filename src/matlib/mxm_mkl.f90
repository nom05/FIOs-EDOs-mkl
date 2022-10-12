! ** matlib/mxm_mkl.f90 >> Interface for matrices product: m x m
!
!  Copyright (c) 2022  Nicolás Otero Martínez - Marcos Mandado Alonso
!  This file is part of the FIOs-EDOs program available in:
!      https://github.com/nom05/FIOs-EDOs
!
!  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
!  by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!  
!  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License along with this code.  If not, see 
!  <http://www.gnu.org/licenses/>.

subroutine         mxm_mkl(a,b,c,m,n,k)

!
!     Info: https://www.math.utah.edu/software/lapack/lapack-blas/dgemm.html
!       * Original call: SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA,
!                                            A, LDA, B, LDB, BETA, C, LDC )
!

  use omp_lib
  use commonmod ,only                              : debug,debuglabel,ilendbglbl
  use SUFR_kinds,only                              : dbl

  implicit none
  integer ( kind = 4 ),intent(in )                :: m,n,k
  real    ( kind = 8 )                            :: wtime,alpha=1._dbl,beta=0._dbl
  real    ( kind = 8 ),intent(in ),dimension(m,k) :: a  !! Original: a(,).               op-a(m,k)
  real    ( kind = 8 ),intent(in ),dimension(k,n) :: b  !! Original: b(,).               op-b(k,n)
  real    ( kind = 8 ),intent(out),dimension(m,n) :: c  !!                                  c(m,n)

  external dgemm

  if (debug) then
     write( *,'(/,1X,A,A)'        ) debuglabel(:ilendbglbl),'MXM_MKL: Starting ...'
     write( *,'(  8X,A,I8,A,I8)') 'The matrix dimension of A is ',m,' x ',k
     write( *,'(  8X,A,I8,A,I8)') 'The matrix dimension of B is ',k,' x ',n
     wtime = omp_get_wtime()
  endif !! (debug) then

  call dgemm('N','N',m,n,k,alpha,a,m,b,k,beta,c,m)

  if (debug) then
     wtime = omp_get_wtime()-wtime
     write(*,'(8X,A,G14.6)')   'Elapsed seconds = ', wtime
     write(*,'(8X,A,/)') 'MXM_MKL: Normal end of execution.'
  endif !! (debug) then

end subroutine     mxm_mkl

