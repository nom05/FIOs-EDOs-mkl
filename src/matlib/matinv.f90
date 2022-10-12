! ** matlib/matinv.f90 >> Interface for the inversion of matrices
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

subroutine matinv90(A,n) !!,nproc)

  use SUFR_system                           ,only : quit_program_error
  use commonmod                             ,only : stealth,nocolor
  implicit none

  integer( kind = 4 ),intent(in)                 :: n !!,nproc
  integer( kind = 4 ),allocatable,dimension(:)   :: IPIV
  real   ( kind = 8 )            ,dimension(n,n) :: A
  real   ( kind = 8 ),allocatable,dimension(:  ) :: WORK
  integer( kind = 4 )                            :: info,error

  external DGETRF,DGETRI

  !!call MKL_SET_NUM_THREADS(nproc)

  allocate(WORK(n),IPIV(n),stat=error)
  if (error.NE.0) call quit_program_error('Not enough memory',1,nocolor)
  call DGETRF(n,n,A,n,IPIV,info)
  if (info.EQ.0) then
     if (.NOT.stealth) write(*,'(/,8X,A)') "- Factorization succeeded"
  else
     call quit_program_error('Factorization failed',1,nocolor)
  endif !! (info.EQ.0) then
  call DGETRI(n,A,n,IPIV,WORK,n,info)
  if (info.EQ.0) then
     if (.NOT.stealth) write(*,'(8X,A,/)') "- Matrix inversion succeeded"
  else
     call quit_program_error('Matrix inversion failed',1,nocolor)
  endif !! (info.EQ.0) then
  deallocate(IPIV,WORK,stat=error)
  if (error.NE.0) call quit_program_error('Releasing process failed',1,nocolor)
end subroutine
