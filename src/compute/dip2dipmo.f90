! ** compute/dip2dipmo.f90 >> Transforming multipolar matrices from basis functions to molecular orbitals
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


subroutine     dip2dipmo
  use omp_lib
  use SUFR_kinds                                  ,only: double,dbl
  use SUFR_system                                 ,only: quit_program_error
  use SUFR_text                                   ,only: uppercase
  use memory_use
  use writemod                                    ,only: print_asym_mat
  use commonmod                                   ,only: nfb,nmo,c,dip,dipmo,dipmotrunc,print_message,iout      &
                                                       , partial_time,ntime,dipmotrunc,nmotrunc,textxyz,momap   &
                                                       , debuglabel,ilendbglbl &
                                                       , nocolor                      & !! loptions(  6)
                                                       , debug,debuglabel,ilendbglbl  & !! loptions(  7)
                                                       , print_large_matrices         & !! loptions( 12)
                                                       , mem_info                    !& !! loptions( 14)

  implicit none

  integer  (kind =      4)                            :: i !! tmp
  character( len =     80)                            :: label
  real     (kind = double),allocatable,dimension(:,:) :: t

  interface
    subroutine mtxmxm_mkl(a,b,c,d,l,m,n,o)
      integer ( kind = 4 ),intent(in )                :: l,m,n,o
      real    ( kind = 8 ),intent(in ),dimension(m,l) :: a
      real    ( kind = 8 ),intent(in ),dimension(m,n) :: b
      real    ( kind = 8 ),intent(in ),dimension(n,o) :: c
      real    ( kind = 8 ),intent(out),dimension(l,o) :: d
    end subroutine mtxmxm_mkl
  end interface

  if (.NOT.allocated(dip)) call quit_program_error('Multipole matrices array is unexpectly deallocated',1,nocolor)

  call myalloc(dipmo,nmo,nmo,3,'dip2dipmo','dipmo',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:FIOs
  do i = 1,3
     if (debug) write(*,'(1X,A,"<+> PERFORMING PRODUCT (dip->dipmo):",1X,A," <+>")') &
                               debuglabel(:ilendbglbl),textxyz(i)
     call    mtxmxm_mkl(c(:,:,1),dip(:,:,i),c(:,:,1),dipmo(:,:,i),nmo,nfb,nfb,nmo) !! 1st array is transposed
  enddo !! i = 1,3
  call myalloc(dipmotrunc,nmotrunc,nmotrunc,3,'dip2dipmo','dipmotrunc',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:FIOs
  dipmotrunc = dipmo(momap(:),momap(:),:3)

  if (print_large_matrices) then
     call myalloc(t,nmo,nmo,'dip2dipmo','t',nocolor=nocolor,verbose=mem_info) !+D:here
     if (debug) PRINT *,debuglabel(:ilendbglbl),'Printing Multipole matrices(MO basis)...'
     do i = 1,3
        t = dipmo(:,:,i)
        label(:) = ''
        label    = 'Multipole matrix (MO) array for component '//uppercase(textxyz(i))
        call print_asym_mat(iout,nmo,nmo,t,trim(label))
     enddo !! i = 1,3
     if (allocated(t)) call mydealloc(t,'dip2dipmo','t',nocolor=nocolor,verbose=mem_info) !+A:here
     if (debug) PRINT *,debuglabel(:ilendbglbl),'Printed Multipole matrices(MO basis)...'
  endif !! (print_large_matrices) then

  call           partial_time(ntime,'AO to MO (mult.mat) spent')

end subroutine dip2dipmo
