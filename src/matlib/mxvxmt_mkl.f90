! ** matlib/mxvxmt_mkl.f90 >> Performing matrices product: m x vector x m(transposed)
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

subroutine mxvxmt_mkl(a,b,c,d,l,m,n)

  use omp_lib
  use SUFR_kinds ,only                             : dbl
  use SUFR_system,only                             : quit_program_error
  use commonmod  ,only                             : debug,debuglabel,ilendbglbl,nocolor

  implicit none

  integer ( kind = 4 )                            :: i
  integer ( kind = 4 ),intent(in )                :: l,m,n
  real    ( kind = 8 )                            :: wtime
  real    ( kind = 8 ),intent(in ),dimension(l,m) :: a
  real    ( kind = 8 ),intent(in ),dimension(m  ) :: b
  real    ( kind = 8 ),intent(in ),dimension(n,m) :: c
  real    ( kind = 8 ),intent(out),dimension(l,n) :: d
  real    ( kind = 8 ),allocatable,dimension(:,:) :: tmp

  interface
    subroutine mxmt_mkl(a,b,c,m,k,n)
      integer ( kind = 4 ),intent(in )                :: m,n,k
      real    ( kind = 8 ),intent(in ),dimension(m,k) :: a
      real    ( kind = 8 ),intent(in ),dimension(n,k) :: b
      real    ( kind = 8 ),intent(out),dimension(m,n) :: c
    end subroutine mxmt_mkl
  end interface

  if (debug) then
     write( *,'(/,1X,A,A)'        ) debuglabel(:ilendbglbl),'MXVXMt_MKL (2nd array transposed): Starting ...'
     write( *,'(  8X,A,I8,A,I8)') 'The matrix dimension of A is ', l,' x ', m
     write( *,'(  8X,A,I8,A,I8)') 'The vector dimension of B is ', m
     write( *,'(  8X,A,I8,A,I8)') 'The matrix dimension of C is ', n,' x ', m
     wtime = omp_get_wtime()
  endif !! (debug) then

  if (size(a,2).NE.size(b  ).OR.m.NE.size(c,2)) &
          call quit_program_error('WRONG MATRIX DIMENSIONS',1,nocolor)

  allocate(tmp(l,m))
  do i = 1,l
  !$omp parallel default ( none ) shared ( a,b,tmp,i )
     !$omp workshare
    tmp(i,:) = a(i,:)*b(:)
     !$omp end  workshare
!$omp end parallel
  enddo !! i = 1,l

  call mxmt_mkl(tmp,c,d,l,m,n)

  deallocate(tmp)

  if (debug) then
     wtime = omp_get_wtime()-wtime
     write(*,'(8X,A,G14.6)')   'Elapsed seconds = ', wtime
     write(*,'(8X,A,/)') 'MXVXM_MKL: Normal end of execution.'
  endif !! (debug) then

end subroutine mxvxmt_mkl
