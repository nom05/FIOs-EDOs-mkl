! ** matlib/mtxmxm_mkl.f90 >> Interface to multiply matrices: m(transposed) x m x m
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


subroutine mtxmxm_mkl(a,b,c,d,l,m,n,o)

  use omp_lib
  use SUFR_kinds ,only                             : dbl
  use SUFR_system,only                             : quit_program_error
  use commonmod  ,only                             : debug,debuglabel,ilendbglbl,nocolor

  implicit none

  integer ( kind = 4 ),intent(in )                :: l,m,n,o
  real    ( kind = 8 )                            :: wtime
  real    ( kind = 8 ),intent(in ),dimension(m,l) :: a
  real    ( kind = 8 ),intent(in ),dimension(m,n) :: b
  real    ( kind = 8 ),intent(in ),dimension(n,o) :: c
  real    ( kind = 8 ),intent(out),dimension(l,o) :: d
  real    ( kind = 8 ),allocatable,dimension(:,:) :: tmp

  interface
    subroutine     mtxm_mkl(a,b,c,k,m,n)
      integer ( kind = 4 ),intent(in )                :: m,n,k
      real    ( kind = 8 ),intent(in ),dimension(k,m) :: a
      real    ( kind = 8 ),intent(in ),dimension(k,n) :: b
      real    ( kind = 8 ),intent(out),dimension(m,n) :: c
    end subroutine mtxm_mkl
    subroutine     mxm_mkl(a,b,c,m,n,k)
      integer ( kind = 4 ),intent(in )                :: m,n,k
      real    ( kind = 8 ),intent(in ),dimension(m,k) :: a
      real    ( kind = 8 ),intent(in ),dimension(k,n) :: b
      real    ( kind = 8 ),intent(out),dimension(m,n) :: c
    end subroutine mxm_mkl
  end interface

  if (debug) then
     write( *,'(/,1X,A,A)'        ) debuglabel(:ilendbglbl),'MtXMXM_MKL: Starting ...'
     write( *,'(  8X,A,I8,A,I8)') 'The matrix dimension of A is ', m,' x ', l
     write( *,'(  8X,A,I8,A,I8)') 'The matrix dimension of B is ', m,' x ', n
     write( *,'(  8X,A,I8,A,I8)') 'The matrix dimension of C is ', n,' x ', o
     wtime = omp_get_wtime()
  endif !! (debug) then

  if (size(a,1).NE.size(b,1).OR.size(b,2).NE.size(c,1)) &
          call quit_program_error('WRONG MATRIX DIMENSIONS',1,nocolor)

  allocate(tmp(l,n))

  call mtxm_mkl(  a,b,tmp,m,l,n)
  call  mxm_mkl(tmp,c,  d,l,o,n)

  deallocate(tmp)

  if (debug) then
     wtime = omp_get_wtime()-wtime
     write(*,'(8X,A,G14.6)')   'Elapsed seconds = ', wtime
     write(*,'(8X,A,/)') 'MtXMXM_MKL: Normal end of execution.'
  endif !! (debug) then

end subroutine mtxmxm_mkl
