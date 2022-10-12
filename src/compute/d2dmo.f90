! ** compute/d2dmo.f90 >> Transforming (derivative of) density matrix (wrt a electric field) from basis functions to MO 
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


subroutine     d2dmo

  use omp_lib
  use SUFR_kinds   ,only: double,dbl
  use SUFR_system  ,only: quit_program_error
  use SUFR_text    ,only: uppercase,int2str
  use memory_use
  use writemod     ,only: print_sym_mat,print_asym_mat
  use commonmod    ,only: nfb,nmo,ci,print_message,partial_time,ntime,nop,lop,ncmpeqv,comptext,d,dmo &
                        , iout,comptext &
                        , edo                          & !! loptions( -4)
                        , compute_dmo                  & !! loptions(  1)
                        , dynamic_alpha                & !! loptions(  2)
                        , nocolor                      & !! loptions(  6)
                        , debug,debuglabel,ilendbglbl  & !! loptions(  7)
                        , print_large_matrices         & !! loptions( 12)
                        , mem_info                    !& !! loptions( 14)

  implicit none

  integer  (kind =  4) :: i !! tmp
  character( len = 80) :: label

  interface
    subroutine mxmxmt_mkl(a,b,c,d,l,m,n,o)
      integer ( kind = 4 ),intent(in )                :: l,m,n,o
      real    ( kind = 8 ),intent(in ),dimension(l,m) :: a
      real    ( kind = 8 ),intent(in ),dimension(m,n) :: b
      real    ( kind = 8 ),intent(in ),dimension(o,n) :: c
      real    ( kind = 8 ),intent(out),dimension(l,o) :: d
    end subroutine mxmxmt_mkl
  end interface

  if (compute_dmo) then
     if (allocated(dmo)) call quit_program_error("DMO allocated unexpectly before computation",1,nocolor)
     call myalloc(dmo,nmo,nmo,nop,'d2dmo','dmo',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:EDOs OR FIOs
  endif !! (compute_dmo) then
  do i = 1,nop
     if (edo) then
        if (.NOT.compute_dmo) call quit_program_error("DMO computation is not enabled for EDOs",1,nocolor)
        if (debug) PRINT *,debuglabel(:ilendbglbl),'<+> PERFORMING PRODUCT (EDO) for DBF',i
        call mxmxmt_mkl(ci(:,:,i),d(:,:,i),ci(:,:,i),dmo(:,:,i),nmo,nfb,nfb,nmo)  !! 3rd array is transposed
     else
        if (ncmpeqv(lop(i)).LE.3) then
           if (compute_dmo) then
              if (debug) write(*,'(1X,A,"<+> PERFORMING PRODUCT (alpha)",I3,":",1X,A," <+>")') &
                                 debuglabel(:ilendbglbl),i,comptext(lop(i))
              call mxmxmt_mkl(ci(:,:,1),d(:,:,i),ci(:,:,1),dmo(:,:,i),nmo,nfb,nfb,nmo)  !! 3rd array is transposed
           else 
              if (debug) PRINT *,debuglabel(:ilendbglbl),'DMO "',trim(comptext(lop(i))),'" is already stored'
           endif !! (compute_dmo) then
        else
           if (debug) write(*,'(1X,A,"<+> PERFORMING PRODUCT (beta)",I3,":",1X,A," <+>")') &
                              debuglabel(:ilendbglbl),i,comptext(lop(i))
           call mxmxmt_mkl(ci(:,:,1),d(:,:,i),ci(:,:,1),dmo(:,:,i),nmo,nfb,nfb,nmo)  !! 3rd array is transposed
        endif !! (ncmpeqv(lop(i)).LE.3) then
     endif !! (edo) then
     if (print_large_matrices) then
        label(:) = ''
        if (edo) then
           label = 'Density matrix '//int2str(i)//' built'
        else
           label = trim(uppercase(comptext(lop(i))))//' derivative wrt untruncated density matrix built'
        endif !! (edo) then
        if (dynamic_alpha) then
           call  print_asym_mat(iout,nfb,nfb,  d(:,:,i),trim(label)//' from BFs')
           call  print_asym_mat(iout,nmo,nmo,dmo(:,:,i),trim(label)//' from MOs')
        else
           call  print_sym_mat (iout,nfb,      d(:,:,i),trim(label)//' from BFs')
           call  print_sym_mat (iout,nmo,    dmo(:,:,i),trim(label)//' from MOs')
        endif !! (dynamic_alpha) then
     endif !! (print_large_matrices) then
  enddo !! i = 1,nop
  if (.NOT.edo) call mydealloc(ci,'d2dmo','ci',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:inter_c2ci

  call           partial_time(ntime,'Transformation of mat. spent')

end subroutine d2dmo
