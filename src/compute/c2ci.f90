! ** compute/c2ci.f90 >> Computing inversion of MO coefficients
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


subroutine     inter_c2ci

  use SUFR_kinds                                  ,only: double,dbl
  use SUFR_system                                 ,only: quit_program_error
  use memory_use
  use commonmod                                   ,only: nfb,nmo,c,ci,print_message,partial_time,ntime,iout &
                                                       , edo                          & !! loptions( -4)
                                                       , nocolor                      & !! loptions(  6)
                                                       , debug,debuglabel,ilendbglbl  & !! loptions(  7)
                                                       , mem_info                    !& !! loptions( 14)
  implicit none

  integer :: i !! tmp

  if (.NOT.allocated(c) ) call quit_program_error('MO coefficients array is unexpectly deallocated',1,nocolor)
  if (     allocated(ci)) &
                 call quit_program_error('Inverted MO coefficients array is unexpectly allocated',1,nocolor)
  i = size(c(1,1,:))
  call myalloc(ci,nmo,nfb,i,'inter_c2ci','ci',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:d2dmo OR EDOs
  call         c2ci(nfb,nmo,c(:,:,1),ci(:,:,1),'')
  if (edo) then
     if (i.NE.2) call quit_program_error('2nd MO coefficients array is unexpectly not included',1,nocolor)
     if (i.NE.2) call quit_program_error('2nd Inv. MO coefficients array is unexpectly not included',1,nocolor)
     call      c2ci(nfb,nmo,c(:,:,2),ci(:,:,2),' (perturbed by field)')
  endif !! (edo) then

  contains

     subroutine     c2ci(nfb,nmo,c,ci,text)
       use omp_lib
       use commonmod                                       ,only: &
                                                                  print_large_matrices         & !! loptions( 12)
                                                                , mem_info                    !& !! loptions( 14)
       use writemod                                        ,only: print_asym_mat
     
       implicit none
     
       integer  (kind =      4)                                :: nfb,nmo
       real     (kind = double)                                :: error
       real     (kind = double),intent(in) ,dimension(nfb,nmo) :: c
       real     (kind = double)            ,dimension(nmo,nfb) :: ci
       real     (kind = double),allocatable,dimension(:,:)     :: t
       character( len =   1000)                                :: label
       character( len =      *),intent(in)                     :: text
          
       interface
         subroutine     mtxm_mkl(a,b,c,k,m,n)
           integer ( kind = 4 ),intent(in )                :: m,n,k
           real    ( kind = 8 ),intent(in ),dimension(k,m) :: a
           real    ( kind = 8 ),intent(in ),dimension(k,n) :: b
           real    ( kind = 8 ),intent(out),dimension(m,n) :: c
         end subroutine mtxm_mkl
         subroutine     matinv90(A,n)
           integer( kind = 4 ),intent(in)                  :: n
           real   ( kind = 8 )             ,dimension(n,n) :: A
         end subroutine matinv90
         subroutine     mxmt_mkl(a,b,c,m,k,n)
           integer ( kind = 4 ),intent(in )                :: m,n,k
           real    ( kind = 8 ),intent(in ),dimension(m,k) :: a
           real    ( kind = 8 ),intent(in ),dimension(n,k) :: b
           real    ( kind = 8 ),intent(out),dimension(m,n) :: c
         end subroutine mxmt_mkl
         subroutine     mxm_mkl(a,b,c,m,n,k)
           integer ( kind = 4 ),intent(in )                :: m,n,k
           real    ( kind = 8 ),intent(in ),dimension(m,k) :: a
           real    ( kind = 8 ),intent(in ),dimension(k,n) :: b
           real    ( kind = 8 ),intent(out),dimension(m,n) :: c
         end subroutine mxm_mkl
       end interface
     
       if (debug)  PRINT *,debuglabel(:ilendbglbl),'c2ci:',nfb,nmo,'+',trim(text),'+'
       if (print_large_matrices) call print_asym_mat(iout,nfb,nmo,c,'MO coefficients'//trim(text) &
                                                    ,numcolopt=7,reverse=.TRUE.)
       call myalloc(t,nmo,nmo,'c2ci','t',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:here
       call    mtxm_mkl(c ,c,t,nfb,nmo,nmo)   !! 1st array is transposed
       call print_message('i',"Performing matrix inversion of MOs coefficients array")
       call matinv90(t,nmo)
       call    mxmt_mkl(t,c ,ci,nmo,nmo,nfb)  !! 2nd array is transposed
       !
       !     Checking the inverted matrix ci
       !
       call    mxm_mkl(ci,c,t,nmo,nmo,nfb)
       if (abs(1._dbl-t(1,1)).GT..0001_dbl) then
          label(:) = ''
          label    = 'One or more diagonal ci identity matrix terms too far away from the optimal value 1.0'
          call quit_program_error(trim(label),1,nocolor)
       endif !! (abs(1._dbl-t(1,1)).GT..0001_dbl) then
       
       error = 0._dbl
       !$omp parallel default( none ) shared ( t,nmo,error ) private ( i )
         !$omp do
         do i=1,nmo
            t(i,i) = 1._dbl-t(i,i)
         enddo !! i=1,nmo
         !$omp end do
         !$omp workshare
         error = maxval(abs(t(:,:)))
         !$omp end workshare
       !$omp end parallel
       
       label(:) = ''
       write(label,'("Inversion error = ",1PE13.6)') error
       write(iout,'(/)')
       call print_message('i',trim(label))
       if (abs(error).GT..0001_dbl) then
          label(:) = ''
          label    = 'One or more diagonal ci identity matrix terms too far away from the optimal value 1.0'
          call quit_program_error(trim(label),1,nocolor)
       endif !! (abs(1._dbl-t(1,1)).GT..0001_dbl) then
       call mydealloc(t,'c2ci','t',nocolor=nocolor,verbose=mem_info)  !+A:here
       if (print_large_matrices) call print_asym_mat(iout,nmo,nfb,ci,'Inv. MO coefficients'//trim(text) &
                                                    ,numcolopt=7,reverse=.TRUE.)
       call           partial_time(ntime,'Matrix inv. spent')
       
     end subroutine c2ci
     
end subroutine inter_c2ci
