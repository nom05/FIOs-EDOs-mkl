! ** modules/writemod.f90 >> Main module with tools to write output info on screen or a file
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


module     writemod

  use SUFR_kinds,only : dbl
  use SUFR_system,only: quit_program_error
  use memory_use
  use commonmod  ,only: ndimtriang,debug,debuglabel,ilendbglbl,nocolor,mem_info

  implicit none

  contains

    subroutine write_p1fromgauss(inp,ncf,S,text,numcolopt)
    
      use SUFR_kinds                        ,only : double !!,dbl
      use SUFR_text                         ,only : int2str

      implicit none
      integer               ,intent(in)          :: inp,ncf
      integer               ,intent(in),optional :: numcolopt
      character( len=     *),intent(in),optional :: text
      character( len=   100)                     :: label,label2
      integer                                    :: istart,irow,l,ir,ierr,i
      integer                                    :: numcol
      real     (kind=double)                     :: S(:)
    
      numcol = 5
      if (present(numcolopt)) numcol = numcolopt
      if (present(text)) then
         write(inp,'(/,3X,">> ",A,1X,"(",A,1X,"elements):",/)') trim(text),trim(int2str(ndimtriang(ncf,.FALSE.)))
      else
         write(inp,'(/,3X,">> ARRAY ("  ,A,1X,"elements):",/)')            trim(int2str(ndimtriang(ncf,.FALSE.)))
      endif !! (present(text)) then
       label(:) = ''
      label2(:) = ''
      label     = '(3X,'//trim(int2str(numcol))//'(I14))'
      label2    = '(3X,I5,'//trim(int2str(numcol))//'(1PE14.6))'
      do istart = 1,ncf,numcol
         write(inp,label) (i,i=istart,min(ncf,istart+numcol-1))
         do irow = istart,ncf
            l  = min(irow,istart+numcol-1)
            ir = ndimtriang(irow,.TRUE.)
            write(inp,label2,iostat=ierr) irow,S(ir+istart:ir+l)
            if (ierr.NE.0) call quit_program_error('Error 1 in write_p1fromgauss',1,nocolor)
         enddo !! irow = istart,ncf
      enddo !! istart = 1,ncf,numcol

    end subroutine write_p1fromgauss

    subroutine     print_sym_mat(inp,ncf1,S,text,numcolopt)
      use memory_use
      implicit none
      integer  (kind=     4),intent(in)                 :: ncf1,inp
      integer  (kind=     4),intent(in) ,optional       :: numcolopt
      integer  (kind=     4)                            :: i,ii,j,numcol
      real     (kind=double),intent(in) ,dimension(:,:) :: S
      real     (kind=double),allocatable,dimension(:  ) :: t
      character( len=     *)            ,optional       :: text
      numcol     = 5
      if (present(numcolopt)) numcol     = numcolopt
      if (.NOT.present(text)) text       = 'ARRAY'
      call myalloc(t,ndimtriang(ncf1,.FALSE.),'print_sym_mat','t',nocolor=nocolor,verbose=mem_info) !+D:here

!$omp parallel default(none) shared ( t,ncf1,S ) private ( i,j,ii )
  !$omp do
      do i = 1,ncf1
         ii = ndimtriang(i,.TRUE.)
         do j = 1,i
            t(ii+j) = S(i,j)
         enddo !! j = 1,ncf1
      enddo !! i = 1,ncf1
  !$omp end do
!$omp end parallel
      call write_p1fromgauss(inp,ncf1,t,trim(text),numcol)
      call mydealloc(t,'print_sym_mat','t',nocolor=nocolor,verbose=mem_info) !+A:here
    end subroutine print_sym_mat

    subroutine write_pmofromgauss(inp,nbf,nom,S,text,numcol)

      use SUFR_kinds,only : double !!,dbl
      use SUFR_text ,only : int2str
   
      implicit none
      integer               ,intent(in) :: inp,nom,nbf
      integer               ,intent(in) :: numcol
      character( len=     *),intent(in) :: text
      character( len=   100)            :: label,label2
      integer                           :: istart,irow,l,ir,ierr,i
      real     (kind=double)            :: S(:)
   
      write(inp,'(/,3X,">> ",A,1X,"(",A,1X,"elements):",/)') trim(text),trim(int2str(nbf*nom))
       label(:) = ''
      label2(:) = ''
      label     = '(3X,'//trim(int2str(numcol))//'(I14))'
      label2    = '(3X,I5,'//trim(int2str(numcol))//'(1PE14.6))'
      do istart = 1,nom,numcol
         l  = min(nom,istart+numcol-1)
         write(inp,label) (i,i=istart,l)
         do irow = 1,nbf
            ir = (irow-1)*nom
            write(inp,label2,iostat=ierr) irow,S(ir+istart:ir+l)
            if (ierr.NE.0) call quit_program_error('Error 1 in write_pmofromgauss',1,nocolor)
         enddo !! irow = istart,ncf
      enddo !! istart = 1,nmo,numcol

    end subroutine write_pmofromgauss

    subroutine     print_asym_mat(inp,ncf1,ncf2,S,text,numcolopt,reverse)
      use memory_use
      implicit none
      integer  (kind=     4),intent(in)                 :: ncf1,ncf2,inp
      integer  (kind=     4),intent(in) ,optional       :: numcolopt
      integer  (kind=     4)                            :: i,j,numcol
      real     (kind=double),intent(in) ,dimension(:,:) :: S
      real     (kind=double),allocatable,dimension(:  ) :: t
      character( len=     *)            ,optional       :: text
      logical               ,intent(in) ,optional       :: reverse
      logical                                           :: reversedef
      reversedef = .FALSE.
      numcol     = 5
      if (present(numcolopt)) numcol     = numcolopt
      if (present(reverse  )) reversedef = reverse
      if (.NOT.present(text)) text       = 'ARRAY'
      call myalloc(t,ncf1*ncf2,'print_asym_mat','t',nocolor=nocolor,verbose=mem_info) !+D:here
   
      if (reversedef) then
!$omp parallel default(none) shared ( t,ncf1,ncf2,S ) private ( i,j )
  !$omp do
         do i = 1,ncf2
            do j = 1,ncf1
               t(ncf1*(i-1)+j) = S(j,i)
            enddo !! j = 1,ncf2
         enddo !! i = 1,ncf1
  !$omp end do
!$omp end parallel
      else
!$omp parallel default(none) shared ( t,ncf1,ncf2,S ) private ( i,j )
  !$omp do
         do i = 1,ncf1
            do j = 1,ncf2
               t(ncf2*(i-1)+j) = S(i,j)
            enddo !! j = 1,ncf2
         enddo !! i = 1,ncf1
  !$omp end do
!$omp end parallel
      endif !! (reversedef) then
      call write_pmofromgauss(inp,ncf1,ncf2,t,trim(text),numcol)
      call mydealloc(t,'print_asym_mat','t',nocolor=nocolor,verbose=mem_info) !+A:here
    end subroutine print_asym_mat

    subroutine     mk_fio_edo_fchk(ipos,ledo,lfio,debugval,nocolorval,meminfo)
    
        use SUFR_kinds                              ,only: double,dbl
        use commonmod                               ,only: ifchk,nmotrunc,ifio,D2,ordD2Acedo,cedo,nmo,nfb,d    &
                                                           , ndimtriang,cedo

        implicit none

        integer(kind =      4),intent(in)               :: ipos
        integer(kind =      4)                          :: icde,jjj,i,j
        logical,optional      ,intent(in)               :: debugval,nocolorval,meminfo,ledo,lfio
        logical                                         :: debugvaldef,nocolorvaldef,meminfo_def,ledodef,lfiodef
        real   (kind = double),allocatable,dimension(:) :: tmp1

        meminfo_def   = .FALSE.
        debugvaldef   = .FALSE.
        nocolorvaldef = .FALSE.
        ledodef       = .FALSE.
        lfiodef       = .FALSE.
        if (present(debugval)  ) debugvaldef   = debugval
        if (present(nocolorval)) nocolorvaldef = nocolorval
        if (present(meminfo)   ) meminfo_def   = meminfo
        if (present(ledo)      ) ledodef       = ledo
        if (present(lfio)      ) lfiodef       = lfio
        call copyuptofchk('Number of independent functions',ifchk,nmotrunc    ,icde,0  ,ifio,jjj,debugval=debugval)
!!                  rewind -- yes = 0 / 1 = no -------------------------------------^^^
!!                  |||||| ---------------------------------------------------------VVV
        call copyuptofchk('Alpha Orbital Energi'           ,ifchk,nmotrunc    ,icde,1  ,ifio,jjj,debugval=debugval)
        write(ifio,'(5(1PE16.8))') D2(ordD2Acedo(:))
        if (allocated(tmp1)) call mydealloc(tmp1,'mk_fio_edo_fchk','tmp1',nocolor=nocolorvaldef,verbose=mem_info) !+A:here
        call myalloc(tmp1,nmo,'mk_fio_edo_fchk','tmp1',nocolor=nocolorvaldef,verbose=meminfo_def) !+D:here
        read (ifchk,*) tmp1(:nmo) !! Skipping values from fchk
        call copyuptofchk('Alpha MO coefficient'           ,ifchk,nfb*nmotrunc,icde,1  ,ifio,jjj,debugval=debugval)
        if (ledodef)      then !! In EDOs, cedo (cv) is already sorted
           write(ifio,'(5(1PE16.8))') ((cedo(i,           j ),i=1,nfb),j=1,nmotrunc)
        else if (lfiodef) then
           write(ifio,'(5(1PE16.8))') ((cedo(i,ordD2Acedo(j)),i=1,nfb),j=1,nmotrunc)
        else

                STOP 'Not implemented'


        endif !!(ledodef) then
        if (allocated(tmp1)) call mydealloc(tmp1,'mk_fio_edo_fchk','tmp1',nocolor=nocolorvaldef,verbose=meminfo_def) !+A:here
        call myalloc(tmp1,nfb*nmo,'mk_fio_edo_fchk','tmp1',nocolor=nocolorvaldef,verbose=meminfo_def) !+D:here
        read (ifchk,*) tmp1(:nfb*nmo)
        call copyfchk('Total SCF Density   ',ifchk,icde,1  ,ifio,jjj,debugval=debugval)
        write(ifio,'(5(1PE16.8))') ((d(i,j,ipos),j=1,i),i=1,nfb)
        if (allocated(tmp1)) call mydealloc(tmp1,'mk_fio_edo_fchk','tmp1',nocolor=nocolorvaldef,verbose=meminfo_def)  !+A:here
        call myalloc(tmp1,ndimtriang(nfb,.FALSE.),'mk_fio_edo_fchk','tmp1',nocolor=nocolorvaldef,verbose=meminfo_def) !+D:here
        read (ifchk,*) tmp1(:ndimtriang(nfb,.FALSE.))
        call mydealloc(tmp1,'mk_fio_edo_fchk','tmp1',nocolor=nocolorvaldef,verbose=meminfo_def) !+A:here
        call copyfchk('ve al final del fich',ifchk,icde,1  ,ifio,jjj,debugval=debugval)
                                                                      !! TRICK: put any inexistent text
                                                                      !!        to copy the rest of the file.

  end subroutine mk_fio_edo_fchk

  subroutine copyfchk(text,inp,icde,irw,iout,jjj,debugval)

    implicit none

    integer          ,intent(in) :: inp,irw,iout
    integer                      :: icde,jjj,iii
    logical,optional ,intent(in) :: debugval
    logical                      :: debugvaldef
    character(len=20),intent(in) :: text
    character(len=85)            :: line

    debugvaldef = .FALSE.
    if (present(debugval)) debugvaldef = debugval
    if (debugval) PRINT *,debuglabel(:ilendbglbl),'In copyfchk: +',trim(text),'+',inp,icde,irw,iout,jjj
    icde = 0
    jjj  = 0
    if (irw.EQ.0) rewind(inp)
    do
      read (inp,'(A)',iostat=iii) line
      if (index(trim(line),trim(text)).GT.0) then
         write(iout,'(A)') trim(line)
         if (debugval) PRINT *,debuglabel(:ilendbglbl),'+',line,'+',text,'+'
         exit
      else if (iii .NE.0) then
         jjj = 1
         exit
      else !! if (iout.NE.0) then
         write(iout,'(A)') trim(line)
      endif
    enddo
    icde = 1

  end subroutine copyfchk

  subroutine copyuptofchk(text                  ,inp  ,ivalue      ,icde,irw,iout,jjj,debugval)

    implicit none

    integer          ,intent(in) :: inp,irw,iout,ivalue
    logical,optional ,intent(in) :: debugval
    logical                      :: debugvaldef
    integer                      :: icde,jjj,iii
    character(len= *),intent(in) :: text
    character(len=85)            :: line

    debugvaldef = .FALSE.
    if (present(debugval)) debugvaldef = debugval
    if (debugvaldef) PRINT *,debuglabel(:ilendbglbl) &
                            ,'In copyuptofchk: +',trim(text),'+',inp,ivalue,icde,irw,iout,jjj
    icde = 0
    jjj  = 0
    if (irw.EQ.0) rewind(inp)
    do
      read (inp,'(A)',iostat=iii) line
      if (index(trim(line),trim(text)).GT.0) then
         write(iout,'(A49,I12)') trim(line),ivalue
         if (debugvaldef) PRINT *,debuglabel(:ilendbglbl),'+',line,'+',text,'+'
         exit
      else if (iii .NE.0) then
         jjj = 1
         exit
      else !! if (iout.NE.0) then
         write(iout,'(A)') trim(line)
      endif
    enddo
    icde = 1

  end subroutine copyuptofchk

  subroutine mk_gnuplot(ignu,fpdf,lbeta,tipo,nhomo,nmo,fdat,sumpol)

     use SUFR_text                             ,only: uppercase,int2str,dbl2str
     use SUFR_system                           ,only: warn
     use commonmod                             ,only: comptext,textxyz,extension

     implicit none
     character( len  = * ),intent(in)              :: fpdf,fdat
     logical              ,intent(in)              :: lbeta
     integer  ( kind = 4 ),intent(in)              :: tipo,nhomo,ignu,nmo
     real     ( kind = 8 ),intent(in),dimension(3) :: sumpol
     real     ( kind = 8 )                         :: cutoff
     integer  ( kind = 4 )                         :: ibeta=1,maxcomp,i,nrepr,ii
     character( len  = 1 )           ,dimension(3) :: abc=(/'a','b','c'/)
     logical                                       :: lpos,lrepresent(3)=.FALSE.

     cutoff  = maxval(abs(sumpol),1)/50._dbl
     maxcomp = maxloc(abs(sumpol),1)
     lpos    = sumpol(maxloc(abs(sumpol),1)).GT.0
     if (lbeta) ibeta = 2
     do i = 1,3
        lrepresent(i) = sumpol(i).LT.-cutoff.OR.sumpol(i).GT.cutoff
     enddo !! i = 1,3

     if (.NOT.any(lrepresent)) then
        call warn('All components are almost zero!!',1,nocolor=nocolor)
        lrepresent = .TRUE.
     endif !! (.NOT.any(lrepresent)) then
     nrepr = count(lrepresent)

     if (debug) PRINT *,debuglabel(:ilendbglbl),'mk_gnuplot(lbeta,tipo,nhomo,maxcomp,lpos,nrepr):' &
                                               ,lbeta,tipo,nhomo,maxcomp,lpos,nrepr

     write(ignu,'(   A )') 'set terminal pdf size 11cm, 7cm dashed font "Helvetica,14"'
     write(ignu,'(   A )') 'set termoption enhanced'
     write(ignu,'( 3(A))') 'set output "',trim(fpdf),'.pdf"'
     write(ignu,'(   A )') 'set autoscale'
     write(ignu,'(   A )') 'set format y "%.2f"'
     if (nrepr.GT.1) then
        if (lpos) then
           write(ignu,'(   A )') 'set key bottom right box opaque'
        else
           write(ignu,'(   A )') 'set key top right box opaque'
        endif !! (lpos) then
        write(ignu,'(   A )') 'set border back'
        write(ignu,'(   A )') 'set key spacing 1.7 font ",12"'
        write(ignu,'( 5(A))') 'set ylabel "{/Symbol ',abc(ibeta),'}_{',trim(uppercase(comptext(tipo))),'} (au)"'
     else
        do i = 1,3
           if (lrepresent(i)) ii = i
        enddo !! i = 1,3
        write(ignu,'(   A )') 'set nokey'
        write(ignu,'( 6(A))') 'set ylabel "{/Symbol ',abc(ibeta),'}_{',trim(uppercase(comptext(tipo))) &
                   ,trim(uppercase(textxyz(ii))),'} (au)"'
     endif !! (nrepr.GT.1) then
     write(ignu,'(   A )') 'set xlabel "MO number"'
     write(ignu,'( 5(A))') 'set arrow from ',int2str(nhomo),', graph 0 to ' &
                                             ,int2str(nhomo),', graph 1 nohead ls 3 lw 2'
     write(ignu,'(   A )') 'set grid'
     write(ignu,'(   A )') '#set yrange[120:-80]'
     write(ignu,'(   A )') '#set xrange[0:80]'
     write(ignu,'(   A )') '#set ytics 40'
     write(ignu,'(   A )') 'set style line 1 \'
     write(ignu,'(   A )') '    linecolor rgb "#0060ad" \'
     write(ignu,'(   A )') '    linetype 1 linewidth 2 \'
     write(ignu,'(   A )') '    pointtype 6 pointsize .5'
     write(ignu,'(   A )') 'set style line 2 \'
     write(ignu,'(   A )') '    linecolor rgb "#aa60ad" \'
     write(ignu,'(   A )') '    linetype 1 linewidth 2 \'
     write(ignu,'(   A )') '    pointtype 4 pointsize .5'
     write(ignu,'(   A )') 'set style line 3 \'
     write(ignu,'(   A )') '    linecolor rgb "#aaad60" \'
     write(ignu,'(   A )') '    linetype 1 linewidth 2 \'
     write(ignu,'(   A )') '    pointtype 2 pointsize .5'
     do i = 1,3
        if (lrepresent(i)) &
           write(ignu,'( 9(A))') 'set label ',int2str(i),' sprintf("%.2f", ',trim(dbl2str(sumpol(i),4)) &
                                 ,') center at first ',int2str(nmo),','     ,trim(dbl2str(sumpol(i),4)) &
                                 ,' font ",9" point pt .7 ps 1 offset -1.5,-1.'
     enddo !! i = 1,3
     write(ignu,'(   A )',advance='no') 'plot '
     do i = 1,3
        if (lrepresent(i)) then
           write(ignu,'(   A )') '\'
           write(ignu,'(14(A))',advance='no') '"',trim(fdat),'.',trim(extension(4)),'" u 1:($',int2str(i+1) &
                                              ,'*1) w linespoints title "{/Symbol ',abc(ibeta),'}_{/=7 ' &
                                              ,trim(uppercase(comptext(tipo))),trim(uppercase(textxyz(i))) &
                                              ,'}" linestyle ',int2str(i),','
        endif !! (lrepresent(i)) then
     enddo !! i = 1,3

  end subroutine mk_gnuplot


end module writemod
