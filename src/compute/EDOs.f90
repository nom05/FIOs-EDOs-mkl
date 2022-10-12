! ** compute/EDOs.f90 >> Main EDOs calculation subroutine. See refs. in manual.
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


subroutine     EDOs

  use omp_lib
  use SUFR_kinds                                  ,only: double,dbl
  use SUFR_system                                 ,only: quit_program_warning
  use SUFR_text                                   ,only: dbl2str,int2str,uppercase
  use memory_use
  use commonmod                                   ,only: nfb,nmo,c,ci,trace2D,print_message,iout,nhomo,dmo,d      &
                                                       , nproc,ordD2Acedo,D2,trace2D,edocomp,momap,nmotrunc,ntime &
                                                       , partial_time,ffchkdef,extension,ifio,fileopen,cv=>cedo   &
                                                       , sqrccutoff &
                                                       , squarec_cutoff              & !! loptions( -6)
                                                       , overwrite_output            & !! loptions(  0)
                                                       , nocolor                     & !! loptions(  6)
                                                       , debug,debuglabel,ilendbglbl & !! loptions(  7)
                                                       , skip_fio_fchk               & !! loptions(  9)
                                                       , print_large_matrices        & !! loptions( 12)
                                                       , print_pops_only             & !! loptions( 13)
                                                       , print_EF_pops               & !! loptions( 14)
                                                       , mem_info                   !& !! loptions( 15)
  use writemod                                    ,only: print_sym_mat,print_asym_mat,mk_fio_edo_fchk

  implicit none

  real     (kind = double),allocatable,dimension(:  ) :: arr
  real     (kind = double),allocatable,dimension(:,:) :: dmov,A
  real     (kind = double)                            :: ctrans,polpos,polneg,poptrans,coef
  integer  (kind =     4 )                            :: i,iend,j,k !! tmp
  integer  (kind =     4 )                            :: nhomotrunc
  integer  (kind =     4 ),allocatable,dimension(:  ) :: indx
  character( len =  1000 )                            :: label,text
  character( len =     3 ),dimension(-3:3),parameter  :: txtfldcmp=(/ 'f-z' &
                                                                    , 'f-y' &
                                                                    , 'f-x' &
                                                                    , 'f0d' &
                                                                    , 'f+x' &
                                                                    , 'f+y' &
                                                                    , 'f+z' /)

  interface
    subroutine     mxvxmt_mkl(a,b,c,d,l,m,n)
      integer ( kind = 4 ),intent(in )                :: l,m,n
      real    ( kind = 8 ),intent(in ),dimension(l,m) :: a
      real    ( kind = 8 ),intent(in ),dimension(m  ) :: b
      real    ( kind = 8 ),intent(in ),dimension(n,m) :: c
      real    ( kind = 8 ),intent(out),dimension(l,n) :: d
    end subroutine mxvxmt_mkl
    subroutine     mxm_mkl(a,b,c,m,n,k)
      integer ( kind = 4 ),intent(in )                :: m,n,k
      real    ( kind = 8 ),intent(in ),dimension(m,k) :: a
      real    ( kind = 8 ),intent(in ),dimension(k,n) :: b
      real    ( kind = 8 ),intent(out),dimension(m,n) :: c
    end subroutine mxm_mkl
    subroutine     indexxabs(B,matsort,nA,labs,nproc,debug)
      integer ( kind = 4 ),intent( in)                :: nA
      integer ( kind = 4 ),intent(out),dimension(nA)  :: matsort
      integer ( kind = 4 ),intent( in),optional       :: nproc
      real    ( kind = 8 ),intent( in),dimension(nA)  :: B
      logical             ,intent( in),optional       :: debug,labs
    end subroutine indexxabs
    subroutine     diasym(a,eig,n)
      integer ( kind = 4 )                            :: n
      real    ( kind = 8 )                            :: a(n,n),eig(n)
    end subroutine diasym
    subroutine     mxmxmt_mkl(a,b,c,d,l,m,n,o)
      integer ( kind = 4 ),intent(in )                :: l,m,n,o
      real    ( kind = 8 ),intent(in ),dimension(l,m) :: a
      real    ( kind = 8 ),intent(in ),dimension(m,n) :: b
      real    ( kind = 8 ),intent(in ),dimension(o,n) :: c
      real    ( kind = 8 ),intent(out),dimension(l,o) :: d
    end subroutine mxmxmt_mkl
    subroutine     mtxmxm_mkl(a,b,c,d,l,m,n,o)
      integer ( kind = 4 ),intent(in )                :: l,m,n,o
      real    ( kind = 8 ),intent(in ),dimension(m,l) :: a
      real    ( kind = 8 ),intent(in ),dimension(m,n) :: b
      real    ( kind = 8 ),intent(in ),dimension(n,o) :: c
      real    ( kind = 8 ),intent(out),dimension(l,o) :: d
    end subroutine mtxmxm_mkl
  end interface

  call print_message('i','Computing EDOs...')
  label(:) = ''
  label    = '(2(/),3X,68("*"),/,3X,&
           &"** ANALYSING ",A," FIELD COMPONENT OF ELECTRON DEFORMATION ORBITALS **",/,&
           &3X,68("*"),/)'
  write(iout,label) trim(uppercase(txtfldcmp(edocomp)))
  
  call myalloc(dmov,nmo,nmo,'EDOs','dmov',nocolor=nocolor,verbose=mem_info) !+D:here
  call myalloc(A   ,nmo,nmo,'EDOs','A'   ,nocolor=nocolor,verbose=mem_info) !+D:here
  call    mxm_mkl(ci(:,:,1),c(:,:,2),A,nmo,nmo,nfb)
  call mydealloc(ci,'EDOs','ci',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:c2ci
  call    mxmxmt_mkl(A(:,:),dmo(:,:,2),A,dmov,nmo,nmo,nmo,nmo)

!$omp parallel default ( none ) shared ( momap,nhomotrunc,nhomo )
  !$omp workshare
  nhomotrunc = count(momap<=nhomo)
  !$omp end workshare
!$omp end parallel
  if (debug) PRINT *,debuglabel(:ilendbglbl),'Truncated ** ',int2str(nmo),'/',int2str(nmotrunc),'=nmo/nmotrunc &
  &nhomo/nhomotrunc=',int2str(nhomo),'/',int2str(nhomotrunc)

  if (print_large_matrices) &
          call  print_sym_mat (iout,nmo,dmov(:,:),'Density matrix 3 built from MOs of 1 and 2')
  write(iout,'(/,1X,"Trace of DMs built from MOs: ",A,1X,A,1X,A)') dbl2str(trace2D(nmo, dmo(:,:,1)),4) &
                                                                 , dbl2str(trace2D(nmo, dmo(:,:,2)),4) &
                                                                 , dbl2str(trace2D(nmo,dmov(:,:  )),4)
  write(iout,'(/)')

  if (print_EF_pops) then
     write(iout,'(/,3X,">> Occupations:",/)')
     if (nmotrunc.NE.nmo) write(iout,'(10X,"(REMEMBER: MOs set was truncated [",A,"/",A,"])",/)') &
                                     int2str(nmotrunc),int2str(nmo     )
     write(iout,'(6X,"MO",1X,I4,1X,"occ = ",F9.6)') (momap(i),dmov(momap(i),momap(i)),i=1,nmotrunc)
  endif !! (print_EF_pops) then
  if (print_pops_only) call quit_program_warning('Stopping under user request ...',1,nocolor)

!
!     Calculating the occupied-->virtual electron transfer
!
  ctrans = trace2D(nmotrunc-nhomotrunc,dmov(momap(nhomotrunc+1:),momap(nhomotrunc+1:)))

!
!     Diagonalizing deformation density matrix, calculating the eigenvectors
!     and writing the corresponding "deformation orbitals", occupation numbers
!     and deformation density matrix in a formatted checkpoint file
!       ** Reuse A to save the deformation DM. **
!
  call mydealloc(A,'EDOs','A',nocolor=nocolor,verbose=mem_info) !+A:here
  call   myalloc(A         ,nmotrunc,nmotrunc,'EDOs','A'         ,nocolor=nocolor,verbose=mem_info) !+D:here
  call   myalloc(D2        ,nmotrunc         ,'EDOs','D2'        ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.)!+D:conductance_calc
  call   myalloc(ordD2Acedo,nmotrunc         ,'EDOs','ordD2Acedo',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.)!+D:conductance_calc
  call   myalloc(cv        ,nfb     ,nmotrunc,'EDOs','cv'        ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.)!+D:conductance_calc
  call   myalloc(arr       ,nmotrunc         ,'EDOs','arr'       ,nocolor=nocolor,verbose=mem_info) !+D:here
  call   myalloc(indx      ,nmotrunc         ,'EDOs','indx'      ,nocolor=nocolor,verbose=mem_info) !+D:here
  if (debug) PRINT *,debuglabel(:ilendbglbl),'nmotrunc from momap=',size(momap)

!$omp  parallel default ( none ) &
!$omp& shared           ( dmo,dmov,A,momap,nmotrunc )
  !$omp workshare
  A(:,:) = dmov(momap(:nmotrunc),momap(:nmotrunc))-dmo(momap(:nmotrunc),momap(:nmotrunc),1)
  !$omp end workshare
!$omp end parallel
  call mydealloc(dmo,'EDOs','dmo',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:d2dmo OR readg09out
  if (print_large_matrices) &
          call  print_sym_mat (iout,nmotrunc,A(:,:),'Truncated density matrix 3 built from MOs of 1 and 2')
  call  partial_time(ntime,'Preliminar steps spent')

  call diasym(A,D2,nmotrunc)
  call partial_time(ntime,'Diagonalization spent')
! 
!     Ordering EDOs by their occupation numbers
! 
  call indexxabs(D2,ordD2Acedo,nmotrunc,.TRUE.,nproc,debug)
!                      if true, abs <---^^^^^^

  if (debug) then
     i = 30
     if (nmotrunc.LT.30) i = nmo
     PRINT*,debuglabel(:ilendbglbl),'ordD2Acedo( 1,:',int2str(i),')=',      ordD2Acedo(:i)
     PRINT*,debuglabel(:ilendbglbl),'      D2  (   :',int2str(i),')=',D2(              :i )
     PRINT*,debuglabel(:ilendbglbl),'      Dord(   :',int2str(i),')=',D2(   ordD2Acedo(:i))
     PRINT*,debuglabel(:ilendbglbl),'      Aord( 1,:',int2str(i),')=', A( 1,ordD2Acedo(:i))
  endif !! (debug) then

  write(iout,'(/,3X,">>",1X,"Deformation orbitals on the basis of unperturbed MOs",/)')
  if (nmotrunc.NE.nmo) write(iout,'(10X,"(REMEMBER: MOs set was truncated)",/)')
  do i = 1,nmotrunc
  !$omp  parallel default ( none ) &
  !$omp& shared           ( coef,poptrans,A,D2,ordD2Acedo,nhomotrunc,i,arr )
     !$omp workshare
     coef     = sum(A(:,i)*A(:,i))
     poptrans = sum(D2(ordD2Acedo(i))*(A(nhomotrunc+1:,ordD2Acedo(i))*A(nhomotrunc+1:,ordD2Acedo(i))))
     arr(:)   = A(:,ordD2Acedo(i))*A(:,ordD2Acedo(i))
     !$omp end workshare
  !$omp end parallel
     call indexxabs(arr,indx,nmotrunc,.FALSE.,nproc,debug)
!                    if true, abs <---^^^^^^^
     write(iout,'(1X,"EDO ",I4,1X,"occ = ",F9.6,1X,"poptrans = ",F9.6,1X,"Norm = ",F9.6)') &
                                            i,D2(ordD2Acedo(i)),poptrans,coef
     write(iout,'("  # ",        &
                 &"position"  ,8X, &
                 &"value(c)"  ,3X, &
                 &"square(c)" ,6X, &
                 &"(HOMO:",3X,A,")")',advance='NO') int2str(nhomo)
     if (squarec_cutoff) then
!$omp parallel default ( none ) shared ( iend,arr,sqrccutoff )
  !$omp workshare
        iend = count(arr>=sqrccutoff)
  !$omp end workshare
!$omp end parallel
        write(iout,'(2X,"Printing",1X,"squared MO coefficients>=",A,"%")') dbl2str(sqrccutoff*100,2)
     else
        iend = 10
        write(iout,'(2X,"Printing",1X,A,1X,"largest squared MO coefficients")') int2str(iend)
     endif !! (squarec_cutoff) then
     write(iout,'(I3,I6,5X,ES14.6,3X,F9.6)') &
          ((                                 &
                 momap(j)                  , &
                  indx(j)                  , &
                A(indx(j),ordD2Acedo(i))   , &
              arr(indx(j))                   &
                                                ),j=1,iend)
     if (nmotrunc.LE.100.OR.print_large_matrices) then
        write(iout,'(/,"MO Coefficients:")',advance='NO')
        if (nmotrunc.NE.nmo) then
            write(iout,'(2X,"(REMEMBER: MOs set was truncated)")')
        else
            write(iout,'(1X)')
        endif !! (nmotrunc.NE.nmo) then
        write(iout,'(5ES14.6)') (A(j,ordD2Acedo(i)),j=1,nmotrunc)
        write(iout,'(/)')
     endif !! (nmotrunc.LE.100) then
  enddo !! i = 1,nmotrunc

  call    mxm_mkl(c(:,momap(:),1),A(:,ordD2Acedo(:)),cv,nfb,nmotrunc,nmotrunc)

  call mydealloc(c    ,'EDOs','c'    ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:readfchk
  call mydealloc(momap,'EDOs','momap',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:MAIN

  if (print_large_matrices) call print_asym_mat(iout,nfb,nmotrunc,cv,'Sorted EDOs coefficients' &
                                                                    ,numcolopt=7,reverse=.TRUE.)

!
!     Writing FCHK file containing electron deformation density and EDOs for latter
!     visualization and use
!
  if (.NOT.skip_fio_fchk) then
     k        = len(trim(ffchkdef))
     ffchkdef = trim(ffchkdef)//'_EDO-'//trim(txtfldcmp(edocomp))
     label(:) = ''
     write(label,'(2(A),".",A)') "Generating fchk file with EDOs: ",trim(ffchkdef),trim(extension(1))
     call print_message('i',trim(label))
     call fileopen(trim(ffchkdef),trim(extension(1)),text,ifio ,&
                   lexist=.FALSE.,sttsfl='UNKNOWN',overwrite=overwrite_output)
     ffchkdef(k+1:len(ffchkdef)) = '' !! Cleaning variable ffchkdef
     call mxvxmt_mkl(cv,D2(ordD2Acedo(:)),cv,d(:,:,3),nfb,nmotrunc,nfb)
     call mk_fio_edo_fchk(  3 ,ledo=.TRUE.,debugval=debug,nocolorval=nocolor)       !! Write fchk
     ffchkdef(k+1:len(ffchkdef)) = '' !! Cleaning variable ffchkdef
     close(ifio)
  else
     call print_message('i',"Skipping fchk file with EDOs")
  endif !! (.NOT.skip_fio_fchk) then

  !$omp  parallel default ( none ) &
  !$omp& shared           ( polpos,D2,polneg )
     !$omp workshare
  polpos = sum(D2,MASK=D2>0._dbl)
  polneg = sum(D2,MASK=D2<0._dbl)
     !$omp end workshare
  !$omp  end parallel

  write(iout,'(/,"N electrons + =", F8.4,&
              &/,"N electrons - =", F8.4,&
              &/," Occ --> Vir  =", F8.4,&
              &/,"Def pop       =",G14.8)') polpos,polneg,ctrans,polpos+polneg

  call mydealloc(indx ,'EDOs','indx',nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(arr  ,'EDOs','arr' ,nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(A    ,'EDOs','A',nocolor=nocolor,verbose=mem_info)    !+A:here
  call mydealloc(dmov ,'EDOs','dmov',nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(d    ,'EDOs','dbf' ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:readfchk

end subroutine EDOs
