! ** compute/FIOs.f90 >> Main subroutine to perform FIOs decomposition analysis
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

subroutine     FIOs

  use omp_lib
  use SUFR_kinds                                  ,only: double,dbl
  use SUFR_system                                 ,only: quit_program_error
  use SUFR_text                                   ,only: dbl2str,int2str,uppercase
  use memory_use
  use commonmod                                   ,only: nfb,nmo,nel,c,print_message,partial_time,ntime,d,dmo &
                                                       , nop,lop,ncmpeqv,comptext,ilendbglbl,debuglabel,dip   &
                                                       , dipmo,dipmotrunc,nmotrunc,textxyz,momap,comptext     &
                                                       , iout,nproc,ffiodat,fgnuplot,extension,idat,trace2D   &
                                                       , ignu,fileopen,ifio,ffchkdef,D2,ordD2Acedo,sqrccutoff &
                                                       , cedo,nhomo &
                                                       , squarec_cutoff         & !! loptions( -6)
                                                       , overwrite_output       & !! loptions(  0)
                                                       , dynamic_alpha          & !! loptions(  2)
                                                       , test_deriv             & !! loptions(  3)
                                                       , stealth                & !! loptions(  5)
                                                       , nocolor                & !! loptions(  6)
                                                       , debug                  & !! loptions(  7)
                                                       , skip_fio_fchk          & !! loptions(  9)
                                                       , mk_gnuplot_scr         & !! loptions( 10)
                                                       , print_large_matrices   & !! loptions( 12)
                                                       , print_EF_pops          & !! loptions( 14)
                                                       , mem_info              !& !! loptions( 15)
  use writemod                                    ,only: mk_fio_edo_fchk,print_sym_mat,print_asym_mat,mk_gnuplot

  implicit none

  integer  (kind =      4)                            :: i,j,k,ii,jj,iend !! tmp
  integer  (kind =      4)                            :: nhomotrunc,iprop
  integer  (kind =      4),allocatable,dimension(:  ) :: indx
  real     (kind = double)                            :: ctrans,poptrans,coef,polx,poly,polz,polpos,polneg,popdef
  real     (kind = double)                            :: poltensor(12,3),sumpol(3)
  real     (kind = double),allocatable,dimension(:  ) :: wMO,arr
  real     (kind = double),allocatable,dimension(:,:) :: A,B,polweight
  character(len  =      2)                            :: comp
  character(len  =      3),allocatable,dimension(:  ) :: MOstat
  character(len  =      5),parameter  ,dimension(2  ) :: prop = (/'Alpha','Beta '/)
  character(len  =   1000)                            :: label,label2,text
  logical                                             :: lbeta=.FALSE.,lalpha=.FALSE.

  interface
    subroutine     diasym(a,eig,n)
      integer ( kind = 4 ) :: n
      real    ( kind = 8 ) :: a(n,n),eig(n)
    end subroutine diasym
    subroutine     mxm_mkl(a,b,c,m,n,k)
      integer ( kind = 4 )                           ,intent(in ) :: m,n,k
      real    ( kind = 8 ),            dimension(m,k),intent(in ) :: a
      real    ( kind = 8 ),            dimension(k,n),intent(in ) :: b
      real    ( kind = 8 ),            dimension(m,n),intent(out) :: c
    end subroutine mxm_mkl
    subroutine     indexxabs(B,matsort,nA,l,n,d)
      integer ( kind = 4 )                           ,intent( in) :: nA
      integer ( kind = 4 ),            dimension(nA) ,intent(out) :: matsort
      integer ( kind = 4 ),optional                  ,intent( in) :: n
      real   (kind=8),                 dimension(nA) ,intent( in) :: B
      logical        ,optional                       ,intent( in) :: d,l
    end subroutine indexxabs
  end interface
! 
!    Loop that performs the FIOs analysis of (hyper)polarizability
! 
  call myalloc(         A,nmotrunc,nmotrunc,'FIOs','A'         ,nocolor=nocolor,verbose=mem_info) !+D:here
  call myalloc(        D2,nmotrunc         ,'FIOs','D2'        ,nocolor=nocolor,verbose=mem_info) !+D:here
  call myalloc(ordD2Acedo,nmotrunc         ,'FIOs','ordD2Acedo',nocolor=nocolor,verbose=mem_info) !+D:here
  call myalloc(      cedo,     nfb,nmotrunc,'FIOs','cedo'      ,nocolor=nocolor,verbose=mem_info) !+D:here
  call myalloc(       wMO,nmotrunc         ,'FIOs','wMO'       ,nocolor=nocolor,verbose=mem_info) !+D:here
  call myalloc(    MOstat,nmotrunc         ,'FIOs','MOstat'    ,nocolor=nocolor,verbose=mem_info) !+D:here
  call myalloc(      indx,nmotrunc         ,'FIOs','indx'      ,nocolor=nocolor,verbose=mem_info) !+D:here
  call myalloc(       arr,nmotrunc         ,'FIOs','arr'       ,nocolor=nocolor,verbose=mem_info) !+D:here
  call myalloc( polweight,nmotrunc,3       ,'FIOs','polweight' ,nocolor=nocolor,verbose=mem_info) !+D:here
  poltensor  = 0._dbl

!$omp parallel default ( none ) shared ( momap,nhomotrunc,nhomo )
  !$omp workshare
  nhomotrunc = count(momap<=nhomo)
  !$omp end workshare
!$omp end parallel

  if (debug) PRINT *,debuglabel(:ilendbglbl),'Truncated ** ',int2str(nmo),'/',int2str(nmotrunc),'=nmo/nmotrunc &
  &nhomo/nhomotrunc=',int2str(nhomo),'/',int2str(nhomotrunc)
  lbeta  = any(ncmpeqv(lop(:nop)).GT.3)
  lalpha = any(ncmpeqv(lop(:nop)).LE.3)
  if (debug) PRINT *,debuglabel(:ilendbglbl),'lalpha=',lalpha,lbeta,'=lbeta'

 lii: do ii=1,nop
           jj = ncmpeqv(lop(ii))
           iprop = 2
           if (jj.LE.3) iprop = 1
           comp   = trim(uppercase(comptext(lop(ii))))
           if (debug) PRINT *,debuglabel(:ilendbglbl),'|=| Performing calc ',ii,'for: ',trim(comp),'(',jj,')'
           label(:) = ''
           i = len_trim(comp)
           j = 0
           label2(:) = ''
           if (iprop.EQ.2) then
              j = 5
              label2(:) = 'HYPER'
           endif !! (iprop.EQ.2) then
           label    = '(2/,3X,'//int2str(i+50+j)//'("*"),/,3X, &
                    &"** ANALYSING ",A," FIELD COMPONENT OF ",A,"POLARIZABILITY **",/,&
                    &3X,'//int2str(i+50+j)//'("*"),/)'
           write(iout,label) trim(comp),trim(label2)

           if (dynamic_alpha) then    !! DYNAMIC ALPHA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          A(i,j)=2*dmo(i,j,ii)  !! if DMO is read, this operation must be carried out.
!                                !!          (DONE in readg09.f90)
! asymm dmo=
!   * sym  dmo = 1/2 * (dmo + dmot) (Real part) (t=tranposed)
!   * skew dmo = 1/2 * (dmo - dmot) (Im.  part) (t=tranposed)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE: skew dmo calculated if tests are enabled only.
  
               if (test_deriv) then   !! TEST. Calculation of skew dmo is useless without testing
                  if (debug) PRINT *,debuglabel(:ilendbglbl),'Running some tests ...'

!$omp parallel default ( none ) shared ( dmo,ii,A,nmotrunc,momap )
  !$omp workshare
                  A(:,:) = (            dmo(momap(:nmotrunc),momap(:nmotrunc),ii)  &
                            - transpose(dmo(momap(:nmotrunc),momap(:nmotrunc),ii)) &
                           )*.5_dbl
  !$omp end workshare
!$omp end parallel

                  write(iout,'(1X,"<+> IMAGINARY PART (only DMO):")')
                  call  sometests(ii,A,dipmotrunc,nmotrunc,trunc_par=500,skip_dbf=.TRUE.)
                  if (nmotrunc.NE.nmo) then
                     call myalloc(B,nmo,nmo,'FIOs','B',nocolor=nocolor,verbose=mem_info) !+D:here
!$omp parallel default ( none ) shared ( dmo,ii,B,nmo )
  !$omp workshare
                  B(:,:) = (            dmo(      :nmo      ,      :nmo      ,ii)  &
                            - transpose(dmo(      :nmo      ,      :nmo      ,ii)) &
                           )*.5_dbl
  !$omp end workshare
!$omp end parallel
                     write(iout,'(/,1X,"<+> IMAGINARY PART (full DMO):")')
                     call sometests(ii,B,dipmo     ,nmo     ,trunc_par=500,skip_dbf=.TRUE.)
                     if (allocated(B)) call mydealloc(B,'FIOs','B',nocolor=nocolor,verbose=mem_info) !+A:here
                  endif !! (nmotrunc.NE.nmo) then
               endif !! (test_deriv) then

!$omp parallel default ( none ) shared ( dmo,ii,A,nmotrunc,momap )
  !$omp workshare
                  A(:,:) = (            dmo(momap(:nmotrunc),momap(:nmotrunc),ii)  &
                            + transpose(dmo(momap(:nmotrunc),momap(:nmotrunc),ii)) &
                           )*.5_dbl
  !$omp end workshare
!$omp end parallel

               if (test_deriv) then
                  write(iout,'(1X,/,"<+> REAL      PART:")')
                  call  sometests(ii,A,dipmotrunc,nmotrunc,trunc_par=500,skip_dbf=.FALSE.)
                  if (nmotrunc.NE.nmo) then
                     call myalloc(B,nmo,nmo,'FIOs','B',nocolor=nocolor,verbose=mem_info) !! +D:here
!$omp parallel default ( none ) shared ( dmo,ii,B,nmo )
  !$omp workshare
                  B(:,:) = (            dmo(      :nmo      ,      :nmo      ,ii)  &
                            + transpose(dmo(      :nmo      ,      :nmo      ,ii)) &
                           )*.5_dbl
  !$omp end workshare
!$omp end parallel
                     write(iout,'(/,1X,"<+> REAL      PART (full DMO):")')
                     call sometests(ii,B,dipmo     ,nmo     ,trunc_par=500,skip_dbf=.TRUE.)
                     if (allocated(B)) call mydealloc(B,'FIOs','B',nocolor=nocolor,verbose=mem_info) !+A:here
                  endif !! (nmotrunc.NE.nmo) then
               endif !! (test_deriv) then
           else

!$omp parallel default ( none ) shared ( dmo,ii,A,momap )
  !$omp workshare
               A(:,:) = dmo(momap(:),momap(:),ii)
  !$omp end workshare
!$omp end parallel

               if (test_deriv) then
                  if (debug) PRINT *,debuglabel(:ilendbglbl),'Running some tests ...'
                  call  sometests(ii,A,dipmotrunc,nmotrunc,trunc_par=500)
                  if (nmotrunc.NE.nmo) then
                     call myalloc(B,nmo,nmo,'FIOs','B',nocolor=nocolor,verbose=mem_info) !+D:here
!$omp parallel default ( none ) shared ( dmo,ii,B )
  !$omp workshare
                     B(:,:) = dmo(:,:,ii)
  !$omp end workshare
!$omp end parallel
                     write(iout,'(/,11X,"--- full DMO ---")')
                     call  sometests(ii,B,dipmo     ,nmo     ,trunc_par=500,skip_dbf=.TRUE.)
                     if (allocated(B)) call mydealloc(B,'FIOs','B',nocolor=nocolor,verbose=mem_info) !+D:here
                  endif !! (nmotrunc.NE.nmo) then
               endif !! (test_deriv) then

           endif !! (dynamic_alpha) then
           if (print_large_matrices) call print_sym_mat(iout,nmotrunc,A, &
                      'Truncated and symmetrized derivative matrix (MOs basis)',numcolopt=7)
           label(:) = ''
           label    = 'Arrays alloc.+tests '//int2str(ii)//' spent'
           call  partial_time(ntime,trim(label))

! 
!     Diagonalizing perturbed density matrix, calculating the eigenvectors
!     and writing the corresponding "perturbation orbitals", occupation numbers
!     and perturbed density matrix in a formatted checkpoint file
! 
           if (debug) then
              i = 30
              if (nmotrunc.LT.30) i = nmotrunc
              PRINT *,debuglabel(:ilendbglbl),'before diag, A( 1,:',int2str(i),')=',A( 1,:i)
           endif !! (debug) then
           call diasym(A,D2,nmotrunc)

           if (print_large_matrices) call print_asym_mat(iout,nmotrunc,nmotrunc,A, &
                      'Unsorted derivative eigenvectors (MOs basis)',numcolopt=7)
           if (debug) then
              PRINT *,debuglabel(:ilendbglbl),'Eigenvector  A( 1,:',int2str(i),')=', A( 1,:i)
              PRINT *,debuglabel(:ilendbglbl),'Eigenvalues D2(   :',int2str(i),')=',D2(   :i)
           endif !! (debug) then
           label(:) = ''
           label = 'Diagonalization '//int2str(ii)//' spent'
           call  partial_time(ntime,trim(label))

! 
!     Calculating FIOs in the basis of AOs
! 
           call mxm_mkl(c(:,momap(:),1),A,cedo,nfb,nmotrunc,nmotrunc)

! 
!     Obtaining the weight of each MO in the whole set of FIOs
! 
           if (print_EF_pops) then
              write(iout,'(/,"Contribution of each MO to the whole set of FIOs:")')
              wMO = 0._dbl

!$omp  parallel default ( none )                          &
!$omp&          shared  ( nmotrunc,D2,A,nel,wMO,MOstat,nhomo,momap ) &
!$omp&          private (  i,j )
  !$omp do
              do i = 1,nmotrunc
                 do j = 1,nmotrunc
                    wMO(i) = wMO(i)+abs(D2(j))*A(i,j)*A(i,j)
                 enddo !! j = 1,nmo
                 if (momap(i).LE.nhomo) then   !! Unrestricted not implemented. Check
                    MOstat(i)='Occ'
                 else
                    MOstat(i)='Vir'
                 endif !! (momap(i).LE.nhomo) then
              enddo !! i = 1,nmotrunc
  !$omp end do
!$omp end parallel

              if (nmotrunc.NE.nmo) write(iout,'(10X,"(REMEMBER: MOs set was truncated)",/)')
              write(iout,'("MO ",I4," = ",ES13.6,3x,A)') (momap(i),wMO(i),trim(MOstat(i)),i=1,nmotrunc)
           endif !! (print_EF_pops) then

! 
!     Ordering FIOs by their occupation numbers
! 
           call indexxabs(D2,ordD2Acedo,nmotrunc,.TRUE.,nproc,debug)
!                          if true, abs <--------^^^^^^

           if (debug) then
              i = 30
              if (nmotrunc.LT.30) i = nmotrunc
              PRINT*,debuglabel(:ilendbglbl),'ordD2Acedo( 1,:',int2str(i),')=',      ordD2Acedo(:i)
              PRINT*,debuglabel(:ilendbglbl),'      Dord(   :',int2str(i),')=',D2(   ordD2Acedo(:i))
              PRINT*,debuglabel(:ilendbglbl),'      Aord( 1,:',int2str(i),')=', A( 1,ordD2Acedo(:i))
           endif !! (debug) then

! 
!     Calculating values
! 
           write(iout,'(/,"Perturbation orbitals on the basis of unperturbed MOs",2(/)   &
                        & ,"Eigenvalue, Electron transfer from occupied to virtual MOs, "&
                        & ,"Normalization (must be 1) ",/                                &
                        & ,"and X,Y,Z contributions to the polarizability tensor",/)') 
           ctrans    = 0._dbl
           polweight = 0._dbl

           do i=1,nmotrunc
              poptrans = 0._dbl
              polx     = 0._dbl
              poly     = 0._dbl
              polz     = 0._dbl

!$omp  parallel default ( none ) &
!$omp&          shared  ( i,nmotrunc,polx,poly,polz,D2,ordD2Acedo,dipmotrunc,poptrans,coef,nhomo,A,nhomotrunc )&
!$omp&          private ( j,k )

  !$omp workshare
              coef     = sum(A(:,i)*A(:,i))
              poptrans = sum(D2(ordD2Acedo(i))*(A(nhomo+1:,ordD2Acedo(i))*A(nhomo+1:,ordD2Acedo(i))))
  !$omp end workshare

  !$omp do reduction(+:polx,poly,polz)
              do j=1,nmotrunc
                 do k=1,nmotrunc
                    polx = polx +         D2(  ordD2Acedo(i))     &
                                * dipmotrunc(              j,k,1) &
                                *          A(j,ordD2Acedo(i))     &
                                *          A(k,ordD2Acedo(i))
                    poly = poly +         D2(  ordD2Acedo(i))     &
                                * dipmotrunc(              j,k,2) &
                                *          A(j,ordD2Acedo(i))     &
                                *          A(k,ordD2Acedo(i))
                    polz = polz +         D2(  ordD2Acedo(i))     &
                                * dipmotrunc(              j,k,3) &
                                *          A(j,ordD2Acedo(i))     &
                                *          A(k,ordD2Acedo(i))
                 enddo !! k=1,nmotrunc
              enddo !! j=1,nmotrunc
  !$omp end do

!$omp end parallel


              write(iout, '(1X,"FIO",1X,I4," occ = ",ES13.6," poptrans = ",ES13.6, &
                           &" Norm = ",ES13.6," PolX = ",ES13.6," PolY = ",ES13.6, &
                           &" PolZ = ",ES13.6)') &
                                 i,D2(ordD2Acedo(i)),poptrans,coef,polx,poly,polz



!$omp  parallel default ( none ) &
!$omp&          shared  ( i,nmotrunc,ordD2Acedo,A,arr,polweight,polx,poly,polz ) &
!$omp&          private ( j )
  !$omp do
              do j = 1,nmotrunc
                       arr(j  ) = A(j,ordD2Acedo(i))*A(j,ordD2Acedo(i))
                 polweight(j,1) = polweight(j,1)+arr(j)*polx
                 polweight(j,2) = polweight(j,2)+arr(j)*poly
                 polweight(j,3) = polweight(j,3)+arr(j)*polz
              enddo !! j = 1,nmo
  !$omp end do
!$omp end parallel

              call indexxabs(arr,indx,nmotrunc,.FALSE.,nproc,debug)
!                        if true, abs <--------^^^^^^^
              write(iout,'("  # ",         &
            &           "position" ,8X,    &
            &           "value(c)" ,3X,    &
            &           "square(c)",6X,    &
            &           "c*c*PolX" ,6X,    &
            &           "c*c*PolY" ,6X,    &
            &           "c*c*PolZ" ,2X,    &
            &           "(HOMO:",3X,A,")")',advance='NO') int2str(nhomo)
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
              write(iout,'(I3,I6,5X,ES14.6,3X,F9.6,ES14.6,ES14.6,ES14.6)') &
                        ((                              &
                                j                     , &
                          momap(indx(j))              , &
                              A(indx(j),ordD2Acedo(i)), &
                            arr(indx(j))              , &
                            arr(indx(j))*polx         , &
                            arr(indx(j))*poly         , &
                            arr(indx(j))*polz           &
                                                          ),j=1,iend)

              if (nmotrunc.LE.100) then
                 write(iout,'(/,"MO Coefficients:")',advance='no')
                 if (nmotrunc.NE.nmo) then
                     write(iout,'(2X,"(REMEMBER: MO set was truncated)")')
                 else
                     write(iout,'(1X)')
                 endif !! (nmotrunc.NE.nmo) then
                 write(iout,'(5ES14.6)') (A(j,ordD2Acedo(i)),j=1,nmotrunc)
                 write(iout,'(/)')
              endif !! (nmo.LE.100) then

              poltensor(jj,1) = poltensor(jj,1) + polx
              poltensor(jj,2) = poltensor(jj,2) + poly
              poltensor(jj,3) = poltensor(jj,3) + polz
              ctrans = ctrans+poptrans
           enddo !! i=1,nmotrunc

           if (mk_gnuplot_scr) then
              k        = len_trim(ffiodat)
              ffiodat  = trim(ffiodat )//'_FIO'//trim(comp)
              fgnuplot = trim(fgnuplot)//'_FIO'//trim(comp)
              label(:) = ''
              write(label,'(A,A,".{",A,",",A,"}")') &
                          "Generating dat and Gnuplot files with FIOs contributions: ",trim(ffiodat) &
                          ,trim(extension(4)),trim(extension(5))
              call print_message('i',trim(label))
              call fileopen(trim(ffiodat),trim(extension(4)),text,idat ,&
                         lexist=.FALSE.,sttsfl='UNKNOWN',overwrite=overwrite_output)
           endif !! (mk_gnuplot_scr) then
           write(iout,'(/," >> Value of pol wrt each MO:")')
           do j=1,nmotrunc
              write(iout,'("Orb.",I6,2X,ES14.6,ES14.6,ES14.6)') momap(j),(polweight(j,i),i=1,3)
              do i = 1,3
                 sumpol(i) = sum(polweight(:j,i))
              enddo !! i = 1,3
              if (mk_gnuplot_scr) write(idat,'(I6,2X,ES14.6,ES14.6,ES14.6)') momap(j),(sumpol(i),i=1,3)
           enddo !! j=1,nmotrunc
           write(iout,'(54("-"),/,"SUM(X,Y,Z)",2X,3(ES14.6))') (sumpol(i),i=1,3)
           if (mk_gnuplot_scr) then
              call fileopen(trim(fgnuplot),trim(extension(5)),text,ignu ,&
                         lexist=.FALSE.,sttsfl='UNKNOWN',overwrite=overwrite_output)
              call mk_gnuplot(ignu,fgnuplot,jj.GT.3,lop(ii),nhomo,nmo,ffiodat,sumpol)
               ffiodat(k+1:) = '' !! Cleaning variable
              fgnuplot(k+1:) = '' !! Cleaning variable
              close(ignu)
              close(idat)
           endif !! (mk_gnuplot_scr) then

! 
!     Sum of polarized electrons
! 
  !$omp  parallel default ( none ) &
  !$omp& shared           ( polpos,D2,polneg )
     !$omp workshare
           polpos = sum(D2,MASK=D2>0._dbl)
           polneg = sum(D2,MASK=D2<0._dbl)
     !$omp end workshare
  !$omp  end parallel

           write(iout,'(/,"Sum of positive and negative eigenvalues of FIOs and total electron&
                     & transfer from occupied to virtual MOs")')
           write(iout,'("N electrons + =",ES14.6,/,"N electrons - =",ES14.6,/," Occ -->&
                     & Vir  =",ES14.6,1(/))') polpos,polneg,ctrans
           if (lbeta) then
              write(iout,'(/,&
                             & "X,Y,Z components of the total hyperpolarizability",/,&
                             & A,A,"X =",ES14.6,/,               &
                             & A,A,"Y =",ES14.6,/,               &
                             & A,A,"Z =",ES14.6                  &
                             & )') trim(prop(iprop)),trim(comp),poltensor(jj,1) &
                                  ,trim(prop(iprop)),trim(comp),poltensor(jj,2) &
                                  ,trim(prop(iprop)),trim(comp),poltensor(jj,3)
           else
              write(iout,'(/,&
                             & "X,Y,Z components of the total polarizability",/,&
                             & A,A,"X =",ES14.6,/,               &
                             & A,A,"Y =",ES14.6,/,               &
                             & A,A,"Z =",ES14.6                  &
                             & )') trim(prop(iprop)),trim(comp),poltensor(jj,1) &
                                  ,trim(prop(iprop)),trim(comp),poltensor(jj,2) &
                                  ,trim(prop(iprop)),trim(comp),poltensor(jj,3)
           endif !! (lbeta) then

  !$omp parallel default ( none ) shared ( popdef,D2 )
     !$omp workshare
           popdef = sum(D2)
     !$omp end workshare
  !$omp end parallel

           label = 'Sum of eigenvalues from FIOs for field component '//trim(comp)//' (must be zero) ='
           write(iout,'(A,F10.6)') trim(label),popdef

           call  partial_time(ntime,'Generation of FIOs spent')
!
!     Writing FCHK file containing field-induced density and FIOs for latter
!     visualization and use
!
           if (.NOT.skip_fio_fchk) then
              k = len(trim(ffchkdef))
              ffchkdef = trim(ffchkdef)//'_FIO'//trim(comp)
              label(:) = ''
              write(label,'(2(A),".",A)') "Generating fchk file with FIOs: " &
                                        , trim(ffchkdef),trim(extension(1))
              call print_message('i',trim(label))
              call fileopen(trim(ffchkdef),trim(extension(1)),text,ifio ,&
                            lexist=.FALSE.,sttsfl='UNKNOWN',overwrite=overwrite_output)
       
              call mk_fio_edo_fchk(ii  ,lfio=.TRUE.,debugval=debug,meminfo=mem_info)       !! Write fchk
              ffchkdef(k+1:len(ffchkdef)) = '' !! Cleaning variable ffchkdef
              close(ifio)
              call  partial_time(ntime,'FIOs fchk file gen. spent')
           else
              call print_message('i',"Skipping fchk file with FIOs")
           endif !! (.NOT.skip_fio_fchk) then

     end do lii

     if (lalpha) then
        call print_message('i',"RECOMPUTED POL  TENS.:")
        write(*,*)
        write(iout,'(/,"RECOMPUTED POL  TENS.:")')
        label(:) = ''
        do i = 1,nop
           if (ncmpeqv(lop(i)).LE.3) then
              label = '('//int2str(ncmpeqv(lop(i)))//'(F10.2,"='//uppercase(trim(comptext(lop(i))))//'",A))'
              if (.NOT.stealth) &
                 write(*   ,label) (poltensor(ncmpeqv(lop(i)),j),uppercase(textxyz(j)),j=1,ncmpeqv(lop(i)))
              write(iout,label) (poltensor(ncmpeqv(lop(i)),j),uppercase(textxyz(j)),j=1,ncmpeqv(lop(i)))
           endif !! (ncmpeqv(lop(i).LE.3) then
        enddo !! i = 1,3
     endif !! (lalpha) then

     if (lbeta) then
        label(:) = ''
        if (.NOT.stealth) then
           write(*,*)
           call print_message('i',"RECOMPUTED BETA TENS.:")
           write(*,*)
        endif !! (.NOT.stealth) then
        write(iout,'(/,"RECOMPUTED BETA TENS.:")')
        do i = 1,nop
           if (ncmpeqv(lop(i)).GT.3) then
              label = '(3(F10.2,"='//uppercase(trim(comptext(lop(i))))//'",A))'
              if (.NOT.stealth) &
                 write(*   ,label) (poltensor(ncmpeqv(lop(i)),j),uppercase(textxyz(j)),j=1,3)
              write(iout,label) (poltensor(ncmpeqv(lop(i)),j),uppercase(textxyz(j)),j=1,3)
           endif !! (ncmpeqv(lop(i)).GT.3) then
        enddo !! i = 1,nop
     endif !! (lbeta) then
  call mydealloc(polweight ,'FIOs','polweight' ,nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(arr       ,'FIOs','arr'       ,nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(indx      ,'FIOs','indx'      ,nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(A         ,'FIOs','A'         ,nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(D2        ,'FIOs','D2'        ,nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(ordD2Acedo,'FIOs','ordD2Acedo',nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(cedo      ,'FIOs','cedo'      ,nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(wMO       ,'FIOs','wMO'       ,nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(MOstat    ,'FIOs','MOstat'    ,nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(dipmotrunc,'FIOs','dipmotrunc',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:dip2dipmo
  call mydealloc(dipmo     ,'FIOs','dipmo'     ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:dip2dipmo
  call mydealloc(dip       ,'FIOs','dip'       ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:proc_dip OR readg09out
  call mydealloc(dmo       ,'FIOs','dmo'       ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:d2dmo OR readg09out
  call mydealloc(d         ,'FIOs','dbf'       ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:readg09out
  call mydealloc(c         ,'FIOs','c'         ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:readfchk
  call mydealloc(momap     ,'FIOs','momap'     ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:MAIN

  contains

  subroutine     sometests(jj,A,M,nom,trunc_par,skip_dbf)

     use SUFR_text                            ,only: int2str,uppercase
     use memory_use
     use writemod                             ,only: print_sym_mat,print_asym_mat

     implicit none
     real(kind=double),allocatable,dimension(:,:) :: test
     real(kind=double),intent(in)                 :: A(:,:),M(:,:,:)
     integer          ,intent(in) ,optional       :: trunc_par
     integer                                      :: trunc_pardef=150
     integer                                      :: i
     integer          ,intent(in)                 :: jj,nom
     logical          ,intent(in) ,optional       :: skip_dbf
     logical                                      :: skip_dbfdef=.FALSE.
     character(len=100)                           :: label

     if (present(trunc_par)) trunc_pardef = trunc_par
     if (present(skip_dbf) ) skip_dbfdef  = skip_dbf
     if (debug) PRINT *,debuglabel(:ilendbglbl),trunc_pardef,'=trunc_par skip_dbf=',skip_dbfdef,'for ',jj

     if (.NOT.skip_dbfdef) then
        if (nfb*(nfb+1)/2.GT.trunc_pardef) then
           write(iout,'(/,&
               &"Skipping Perturbed DBF on the basis of&
               & unperturbed BFs (",A,">",A," values)")') trim(int2str(nfb*(nfb+1)/2)),trim(int2str(trunc_pardef))
        elseif (print_large_matrices) then
           write(iout,'(/,"Skipping Perturbed DBF on the basis of unperturbed BFs (already printed)")')
        else
           label = trim(uppercase(comptext(lop(jj))))//' derivative wrt untruncated density matrix built'
           call  print_sym_mat (iout,nfb,      d(:,:,jj),trim(label)//' from BFs')
        endif !! (nfb*(nfb+1)/2.LE.trunc_pardef) then
     endif !! (.NOT.skip_dbfdef) then

     if (nom*nom.GT.trunc_pardef) then
        write(iout,'(/,&
             &"Skipping Perturbed DMO on the basis of&
             & unperturbed MOs (",A,">",A," values)")') trim(int2str(nmo*nmo)),trim(int2str(trunc_pardef))
     elseif (print_large_matrices) then
        write(iout,'(/,"Skipping Perturbed DMO on the basis of unperturbed MOs (already printed)")')
     else
        call print_asym_mat(iout,nom,nom,A,'Perturbed DMO on the basis of unperturbed MOs',numcolopt=7)
     endif !! (nom*(nom+1)/2.LE.trunc_pardef) then

     write(iout,'(/,"TESTS:")')
     if (.NOT.skip_dbfdef) then
        write(iout,'(3X,">> Trace of the DBF =",F12.6)') trace2D(nfb,d(:,:,jj))
        if (allocated(test)) call mydealloc(test,'sometests','test',nocolor=nocolor,verbose=mem_info) !+A:here
        call myalloc(test,nfb,nfb,'sometests','test',nocolor=nocolor,verbose=mem_info) !+D:here
        do i = 1,3
           call mxm_mkl(  d(:,:,jj),  dip(:,:,i),test,nfb,nfb,nfb)
           write(iout,'(3X,">> tr(DBF*DIPBF",A,")   = ",G12.4)') trim(textxyz(i)),trace2D(nfb,test)
        enddo !! i = 1,3
     endif !! (.NOT.skip_dbfdef) then

     write(iout,'(3X,">> Trace of the DMO =",F12.6)') trace2D(nom,A)
     if(allocated(test)) call mydealloc(test,'sometests','test',nocolor=nocolor,verbose=mem_info) !+A:here
     call myalloc(test,nom,nom,'sometests','test',nocolor=nocolor,verbose=mem_info) !+D:here
     do i = 1,3
        call mxm_mkl(A(:,:),M(:,:,i),test,nom,nom,nom)
        write(iout,'(3X,">> tr(DMO*DIPMO",A,")   = ",G12.4)') trim(textxyz(i)),trace2D(nom,test)
     enddo !! i = 1,3
     if(allocated(test)) call mydealloc(test,'sometests','test',nocolor=nocolor,verbose=mem_info) !+A:here
  end subroutine sometests

end subroutine FIOs
