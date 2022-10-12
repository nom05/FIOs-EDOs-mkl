! ** modules/commonmod.f90 >> FIOs-EDOs Main module for arrays, options, subroutines and functions
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

 module commonmod

  use SUFR_kinds,only : double,dbl

  implicit none

  character(len=6)                  ,public  :: version='2022.1'
  integer( kind = 4      )          ,public  :: nproc=1,nprocmkl=1,nato=0,nel=0,nfb=0,nmo=0,nsh=0,nhomo=0 &
                                              , nalpha=0,nbeta=0,ntime=0,nheavy=0,nop=0,nprimsh=0,ilendbglbl=19 &
                                              , edocomp=0,nmotrunc=0
  integer( kind = 4      ),parameter,public  :: ifchk=10,igau=11,iout=12,ifio=13,idat=14,ignu=15
  real   ( kind = double )          ,public  :: Field = 0._dbl,cptim0,rmempart,sqrccutoff=0._dbl
  real   ( kind = double ),dimension(100)    :: ptime = 0._dbl

  integer,dimension(-4:4),parameter          :: nbfsize=(/  9 &  !! -4 (g9 )
                                                         ,  7 &  !! -3 (f7 )
                                                         ,  5 &  !! -2 (d5 )
                                                         ,  4 &  !! -1 (sp )
                                                         ,  1 &  !!  0 ( s )
                                                         ,  3 &  !!  1 ( p )
                                                         ,  6 &  !!  2 (d6 )
                                                         , 10 &  !!  3 (f10)
                                                         , 15 /) !!  4 (g15)
  integer,                parameter          :: maxnbf   =  (size(nbfsize)-1)/2 !! Up to g functions 
                                                                                !! are currently supported

  integer,dimension( 0:9),parameter          :: ncmpeqv=(/  0 &  !!  EDO
                                                         ,  1 &  !!  x -> 0,1
                                                         ,  2 &  !!  y -> 0,2
                                                         ,  3 &  !!  z -> 0,3
                                                         ,  4 &  !! xx -> 1,1
                                                         ,  5 &  !! xy -> 1,2
                                                         ,  6 &  !! xz -> 1,3
                                                         ,  8 &  !! yy -> 2,2
                                                         ,  9 &  !! yz -> 2,3
                                                         , 12 /) !! zz -> 3,3

  integer,dimension(   9),parameter          :: neqread=(/  1 &  !!  x -> 0,1 f0
                                                         ,  1 &  !!  y -> 0,2 f0
                                                         ,  1 &  !!  z -> 0,3 f0
                                                         ,  2 &  !! xx -> 1,1 fx
                                                         ,  2 &  !! xy -> 1,2 fx
                                                         ,  2 &  !! xz -> 1,3 fx
                                                         ,  3 &  !! yy -> 2,2 fy
                                                         ,  3 &  !! yz -> 2,3 fy
                                                         ,  4 /) !! zz -> 3,3 fz

  integer,dimension(   9),parameter          :: nxyzcom=(/  1 &  !!  x -> 0,1 f0
                                                         ,  2 &  !!  y -> 0,2 f0
                                                         ,  3 &  !!  z -> 0,3 f0
                                                         ,  1 &  !! xx -> 1,1 fx
                                                         ,  2 &  !! xy -> 1,2 fx
                                                         ,  3 &  !! xz -> 1,3 fx
                                                         ,  2 &  !! yy -> 2,2 fy
                                                         ,  3 &  !! yz -> 2,3 fy
                                                         ,  3 /) !! zz -> 3,3 fz

  integer                                    :: nclas(4)

  character(len=19)                          :: debuglabel='  ('//achar(27)//'[90mDG'//achar(27)//'[0m'//') '
  character(len= 1),dimension(0:3),parameter :: textxyz  =(/' ','x','y','z'/)
  character(len= 4),dimension(  5),parameter :: extension=(/'fchk','out ','mat ','dat ','gnu '/)
  character(len= 3),dimension(0:9),parameter :: comptext =(/ 'EDO', 'x  ' ,'y  ' ,'z  ' &
                                                           , 'xx ' &
                                                           , 'xy ' &
                                                           , 'xz ' &
                                                           , 'yy ' &
                                                           , 'yz ' &
                                                           , 'zz ' /)
  logical,public           ,dimension(-6:16) :: loptions =(/        &  !! if TRUE:
                                                            .FALSE. &  !!  - -6 - Square(c) cut-off.
                                                          , .FALSE. &  !!  - -5 - EDOs + conductance.
                                                          , .FALSE. &  !!  - -4 - EDOs.
                                                          , .FALSE. &  !!  - -3 - Truncate DMO.
                                                          , .FALSE. &  !!  - -2 - Spec. comps. will be computed
                                                          , .FALSE. &  !!  - -1 - Field introduced by hand.
                                                          , .FALSE. &  !!  -  0 - Overwrite output files.
                                                          , .FALSE. &  !!  -  1 - DMO will be computed.
                                                          , .FALSE. &  !!  -  2 - Dynamic alpha FIOs.
                                                          , .FALSE. &  !!  -  3 - Test w/ DM deriv and dip.mats.
                                                          , .FALSE. &  !!  -  4 - Read multipole matrices.
                                                          , .FALSE. &  !!  -  5 - Stealth mode.
                                                          , .FALSE. &  !!  -  6 - No colors on screen
                                                          , .FALSE. &  !!  -  7 - Debug
                                                          , .FALSE. &  !!  -  8 - Output file name
                                                          , .FALSE. &  !!  -  9 - Skip FIOs fchk files writing
                                                          , .FALSE. &  !!  - 10 - Make Gnuplot scripts
                                                          , .FALSE. &  !!  - 11 - Estimate memory and stop
                                                          , .FALSE. &  !!  - 12 - Print large matrices
                                                          , .FALSE. &  !!  - 13 - Print deformation pops. and stop
                                                          , .FALSE. &  !!  - 14 - Print def. pop. (EDOs) or MO contrib (FIOs)
                                                          , .FALSE. &  !!  - 15 - Print detailed memory info
                                                          , .FALSE. &  !!  - 16 - Print conductance matrices
                                                         /)
  logical,public                            :: truncate_dmo,spec_comps,lfield,overwrite_output            &
                                             , compute_dmo,dynamic_alpha,test_deriv,read_mult_mat,stealth &
                                             , nocolor,debug,output_file,skip_fio_fchk,mk_gnuplot_scr,edo &
                                             , lconduc,estimate_mem,print_large_matrices,print_pops_only  &
                                             , mem_info,print_cond_mat,squarec_cutoff,print_EF_pops
  equivalence                                  (loptions( -6),squarec_cutoff      ) &
                                             , (loptions( -5),lconduc             ) &
                                             , (loptions( -4),edo                 ) &
                                             , (loptions( -3),truncate_dmo        ) &
                                             , (loptions( -2),spec_comps          ) &
                                             , (loptions( -1),lfield              ) &
                                             , (loptions(  0),overwrite_output    ) &
                                             , (loptions(  1),compute_dmo         ) &
                                             , (loptions(  2),dynamic_alpha       ) &
                                             , (loptions(  3),test_deriv          ) &
                                             , (loptions(  4),read_mult_mat       ) &
                                             , (loptions(  5),stealth             ) &
                                             , (loptions(  6),nocolor             ) &
                                             , (loptions(  7),debug               ) &
                                             , (loptions(  8),output_file         ) &
                                             , (loptions(  9),skip_fio_fchk       ) &
                                             , (loptions( 10),mk_gnuplot_scr      ) &
                                             , (loptions( 11),estimate_mem        ) &
                                             , (loptions( 12),print_large_matrices) &
                                             , (loptions( 13),print_pops_only     ) &
                                             , (loptions( 14),print_EF_pops       ) &
                                             , (loptions( 15),mem_info            ) &
                                             , (loptions( 16),print_cond_mat      )!&
  character(len= 36),public,dimension(-6:16):: textopts  =(/                                        &
                                                             "Square(c) cut-off in percentage     " &  !!  - -6 - !  1
                                                           , "Computation of EDOs and conductance " &  !!  - -5 - !  2
                                                           , "Computation of EDOs                 " &  !!  - -4 - !  3
                                                           , "Truncated set of MOs to build DMO   " &  !!  - -3 - !  4
                                                           , "Specific components to be computed  " &  !!  - -2 - !  5
                                                           , "Field introduced by hand            " &  !!  - -1 - !  6
                                                           , "Overwrite output file               " &  !!  -  0 - !  7
                                                           , "DMO will be completely computed     " &  !!  -  1 - !  8
                                                           , "Dynamic alpha FIOs                  " &  !!  -  2 - !  9
                                                           , "Test using DM deriv. and dip. mats. " &  !!  -  3 - ! 10
                                                           , "Multipole matrices read from G09 log" &  !!  -  4 - ! 11
                                                           , "Print on screen warnings and errors " &  !!  -  5 - ! 12
                                                           , "No colors on screen                 " &  !!  -  6 - ! 13
                                                           , "Debug messages on screen            " &  !!  -  7 - ! 14
                                                           , "Change name of the output file      " &  !!  -  8 - ! 15
                                                           , "Skip FIOs fchk files writing        " &  !!  -  9 - ! 16
                                                           , "Make Gnuplot scripts                " &  !!  - 10 - ! 17
                                                           , "Estimate memory use                 " &  !!  - 11 - ! 18
                                                           , "Print large matrices                " &  !!  - 12 - ! 19
                                                           , "Print deform. occup. and stop (EDO) " &  !!  - 13 - ! 20
                                                           , "Print def occ.(EDO)/MO contrb.(FIO) " &  !!  - 14 - ! 21
                                                           , "Print detailed memory info          " &  !!  - 15 - ! 22
                                                           , "Print conductance matrices          " &  !!  - 16 - ! 23
                                                          /)
  character(len  = 105)                     :: ffiodat,fgnuplot,ffchkdef
  integer,public  ,dimension(  9)           :: lop
  logical,public  ,dimension(  9)           :: lcomponents=.FALSE.

  integer( kind = 4      ),dimension(:    ),allocatable        :: itsh,icsh,inpps,momap,ordD2Acedo
  real   ( kind = double ),dimension(:    ),allocatable        :: D2,wMO,concoeffs,primexps,pspconcoeffs
  real   ( kind = double ),dimension(:    ),allocatable,target :: x,y,z
  real   ( kind = double ),dimension(:,:  ),allocatable        :: cedo,over
  real   ( kind = double ),dimension(:,:,:),allocatable        :: c,ci
  real   ( kind = double ),dimension(:,:,:),allocatable        :: dip,d,dmo,dipmo,dipmotrunc

  type,public                                                  :: nat
    integer                                                    :: ifb,ifbend,ifbstart
  end type                                                        nat
  type(nat)        ,allocatable,dimension(:)                   :: ato2basis

contains

  integer function ndimtriang(n,lminus)
 !  Dimension of triangular array w/ (FALSE) or w/o (TRUE) the diagonal
        implicit none
        integer,intent(in) :: n
        logical,intent(in) :: lminus
        if (.NOT.lminus) then
           ndimtriang = n*(n+1)/2
        else
           ndimtriang = n*(n-1)/2
        endif !! (.NOT.lminus) then
  end function ndimtriang

  real ( kind = double ) function trace2D(ndim,mat)
!   Trace of 2D array
        use omp_lib
        implicit none
        integer,intent(in)                               :: ndim
        integer                                          :: i
        real ( kind = double ),intent(in),dimension(:,:) :: mat
        trace2D = 0._dbl

!$omp parallel default( none ) shared( ndim,trace2D,mat ) private( i )
  !$omp do reduction(+:trace2D)
           do i = 1,ndim
              trace2D = trace2D + mat(i,i)
           enddo !! i = 1,ndim
  !$omp end do
!$omp end parallel

  end function trace2D

  subroutine     fileopen(filga,exten,text,inp,lexist,sttsfl    &
                                                     ,accsfl    &
                                                     ,formfl    &
                                                     ,overwrite )
 
      use SUFR_system   ,only: file_open_error_quit,quit_program_error
      implicit none

      integer                                 :: ios
      integer            ,intent(in)          :: inp
      character(len=*   ),intent(in)          :: filga,exten
      character(len= 10 )                     :: extennew
      character(len=*   )                     :: text
      character(len=*   ),intent(in),optional :: sttsfl,accsfl,formfl
      character(len=100 )                     :: sttsfldef,accsfldef,formfldef
      logical            ,intent(in),optional :: lexist,overwrite
      logical                                 :: lchk,llexist,overwritedef

      sttsfldef    ='OLD'
      accsfldef    ='SEQUENTIAL'
      formfldef    ='FORMATTED'
      llexist      = .TRUE.
      overwritedef = .FALSE.
      if (present(sttsfl   )) sttsfldef    = sttsfl
      if (present(accsfl   )) accsfldef    = accsfl
      if (present(formfl   )) formfldef    = formfl
      if (present(lexist   )) llexist      = lexist
      if (present(overwrite)) overwritedef = overwrite

      if (debug) PRINT *,debuglabel(:ilendbglbl),'Entering fileopen: filga,exten,inp: +',trim(filga) &
                        ,'+',trim(exten),'+',inp,llexist,overwritedef

      extennew(:) = ''
      extennew    = trim(exten)
      if (llexist) then
         inquire(file=trim(filga),exist=lchk)
         if (.NOT.lchk) then
            inquire(file=trim(filga)//'.'//trim(exten),exist=lchk)
            if (.NOT.lchk) call file_open_error_quit("Requested file does not exist!!",0,1)
            text(:)     = trim(filga)
         else
            text     = filga(1:index(trim(filga),'.',.TRUE.)-1)
            extennew = filga(index(trim(filga),'.',.TRUE.)+1:len_trim(filga))
            if (trim(extennew).EQ.trim(exten)) then
               extennew(:) = ''
               extennew    = trim(exten)
            endif !! (trim(extennew).EQ.trim(exten)) then
         endif !! (.NOT.lchk) then (including extension)
      else if (.NOT.overwrite) then
         text(:)     = trim(filga)
         inquire(file=trim(text)//'.'//trim(extennew),exist=lchk)
         if (lchk) call quit_program_error('Specified output file will not be overwritten',1,nocolor)
      else
         text(:)     = trim(filga)
      endif !! (lexist) then
      if (debug) PRINT *,debuglabel(:ilendbglbl),'final text: +',trim(text),'+'
      if (debug) PRINT *,debuglabel(:ilendbglbl),'final file: +',trim(text),'.',trim(extennew),'+'
      open (unit=inp,file=trim(text)//'.'//trim(extennew)        &
                                     ,iostat   = ios             &
                                     ,status   = trim(sttsfldef) &
                                     ,access   = trim(accsfldef) &
                                     ,form     = trim(formfldef) )
      if (ios.NE.0) call file_open_error_quit(trim(filga),1,1)

  end subroutine fileopen

  subroutine     partial_time(n    ,text,lfinal,iprint)

      use omp_lib
      use SUFR_system     ,only         : quit_program_warning,warn
      use SUFR_text       ,only         : dbl2str,int2str

      implicit none

      integer  ( kind=4   )            :: n,nlength
      integer  ( kind=4   ),parameter  :: ntotchar=30
      character(  len=*   ),intent(in) :: text
      character(  len=100 )            :: label
      integer     ,optional,intent(in) :: iprint
      integer                          :: iprintdef=0
      logical     ,optional,intent(in) :: lfinal
      logical                          :: lfinaldef=.FALSE.

      if (present(lfinal)) lfinaldef = lfinal
      if (present(iprint)) then
         if (iprint.NE.0.OR.iprint.NE.5.OR.iprint.NE.6) iprintdef = iprint
      endif !! (present(iprint)) then
      if (.NOT.lfinaldef) then
         if (n.GT.size(ptime)) call quit_program_warning('Increase dimension of ptime',1,nocolor)
         if (n.EQ.0) then
            cptim0 = omp_get_wtime()
         else
            ptime(n) = omp_get_wtime()-cptim0-sum(ptime)
            nlength  = len_trim(text)
            if (debug) PRINT *,debuglabel(:ilendbglbl),'+',trim(text),'+ (',len_trim(text),') ntotchar set to=',ntotchar
            if (nlength.GT.ntotchar) then
               call warn('Text in partial_time too long')
               write(label, &
                   '("This step spent",1X,18("."),1X,A," seconds")') &
                              dbl2str(ptime(n),2)
            else
                label(:) = ''
                label    = '(A,1X,'//int2str(33-nlength)//'("."),1X,A," seconds")'
                write(label,label) trim(text),dbl2str(ptime(n),2)
            endif !! (nlength.GT.ntotchar) then
         endif !! (n.EQ.0) then
         n = n+1
      else
               write(label, &
                   '("TOTAL ELAPSED TIME",1X,15("."),1X,A," seconds")') &
                              dbl2str(                        &
                                       omp_get_wtime()-cptim0 &
                                     ,2)
      endif !! (.NOT.lfinaldef) then
      if (n.GT.1) then
         if (iprintdef.NE.0) then
            call print_message('t',trim(label),iprintdef)
         else
            call print_message('t',trim(label))
         endif !! (iprintdef.NE.0) then
      endif !! (n.GT.1) then
      if (present(iprint)) iprintdef=0

  end subroutine partial_time

  subroutine intcomparegt(left,nvall,nvalr,right,label)
     implicit none
     integer  (kind =   4),intent(in) :: nvall,nvalr
     character(len  =   *),intent(in) :: left,right
     character(len  =   6)            :: charint
     character(len  = 100)            :: label
     label(:) = ''
     write(charint,'(I6)') nvall
     label = trim(left)//' ('//trim(charint(verify(charint,' '):6))//') > '//trim(right)//' ('
     write(charint,'(I6)') nvalr
     label = trim(label)//trim(charint(verify(charint,' '):6))//')!!'
  end subroutine intcomparegt

  subroutine     process_text(text)
     use SUFR_system,only : quit_program_error
     use SUFR_text  ,only : lowercase
     implicit none
     integer             :: i,j,l
     character(len=100)  :: text,word,label
     logical             :: lalpha=.FALSE.,lbeta=.FALSE.
     l = len(trim(text))
     i = 1
     if (debug) PRINT *,debuglabel(:ilendbglbl),'---',trim(text),'--- (initial)'
     text = lowercase(text)
     if (debug) PRINT *,debuglabel(:ilendbglbl),'---',trim(text),'--- (lowercase)'
     do while (i.LE.len(trim(text)))
        word(:) = ''
        j = index(trim(text(i:l)),';')-2+i
        if (j.EQ.i-2) j = l
        word = trim(text(i:j))
        if (debug) PRINT *,debuglabel(:ilendbglbl),'+',trim(word),'+',i,j,j+2
        if (lalpha.OR.index(trim(word),'a').GT.0) then
           if (index(trim(word),'a').GT.0) word = trim(word(2:))
           if (index(trim(word),'('  ).GT.0) then
              lalpha = .TRUE.
              word   = trim(word(index(trim(word),'(')+1:))
           endif
           if (index(trim(word),')'  ).GT.0) then
              lalpha = .FALSE.
              word   = trim(word(:index(trim(word),')')-1))
           endif
           if (debug) PRINT *,debuglabel(:ilendbglbl),trim(word)
           if (len(trim(word)).EQ.0) then
                     lcomponents(:3) = .TRUE.
           else
              select case(trim(word))
                 case('x')
                     lcomponents( 1) = .TRUE.
                 case('y')
                     lcomponents( 2) = .TRUE.
                 case('z')
                     lcomponents( 3) = .TRUE.
                 case default
                     label = 'Syntax error in alpha (**'//trim(word)//'**)'
                     call quit_program_error(trim(label),1,nocolor)
              end select !! case(trim(word))
           endif !! (  len(trim(word)).LT.2)       then
        else if (lbeta .OR.index(trim(word),'b').GT.0) then
           if (index(trim(word),'b'  ).GT.0) word = trim(word(2:))
           if (index(trim(word),'('  ).GT.0) then
              lbeta = .TRUE.
              word  = trim(word(index(trim(word),'(')+1:))
           endif
           if (index(trim(word),')'  ).GT.0) then
              lbeta = .FALSE.
              word  = trim(word(:index(trim(word),')')-1))
           endif
           if (debug) PRINT *,debuglabel(:ilendbglbl),trim(word)
           if (len(trim(word)).EQ.0) then
              lcomponents(4:)   = .TRUE.
           else
              select case(trim(word))
                 case('xx')
                     lcomponents( 4) = .TRUE.
                 case('xy')
                     lcomponents( 5) = .TRUE.
                 case('xz')
                     lcomponents( 6) = .TRUE.
                 case('yy')
                     lcomponents( 7) = .TRUE.
                 case('yz')
                     lcomponents( 8) = .TRUE.
                 case('zz')
                     lcomponents( 9) = .TRUE.
                 case default
                     label = 'Syntax error in beta  (**'//trim(word)//'**)'
                     call quit_program_error(trim(label),1,nocolor)
              end select !! case(trim(word))
           endif !! (  len(trim(word)).LT.2)       then
        else
           label = 'Syntax error (**'//trim(word)//'**)'
           call quit_program_error(trim(label),1,nocolor)
        endif
        i = j+2
     enddo
     if (debug) PRINT *,debuglabel(:ilendbglbl),'lcomponents=',lcomponents
  end subroutine process_text

  subroutine     hmnumb(text,iato,ifrag)
     use SUFR_system  ,only         : quit_program_error

     implicit none
     integer          ,intent(out) :: iato,ifrag
     integer                       :: i,j,k,n,iii
     character(len=20)             :: charint
     character(len=* ),intent(in)  :: text
     logical                       :: filex

     if (debug) PRINT *,debuglabel(:ilendbglbl),'hmnumb: +',trim(text),'+','(',len_trim(text),')'
     iato       = 0
     ifrag      = 0
     j          = 1
     charint(:) = ''
     do i = 1,len_trim(text)
        if (text(i:i).EQ.';') then
            iato  = iato+1
            ifrag = ifrag+1
            filex = .TRUE.
        else
            filex = .FALSE.
            charint(j:j) = text(i:i)
            j = j+1
        endif !! (text(i:i).EQ.';') then
        if (filex.OR.i.EQ.len_trim(text)) then
           j = 1
           if (debug) PRINT *,debuglabel(:ilendbglbl),'Processing numbers... ',trim(charint)
           if (index(trim(charint),'-').GT.0) then
              read (charint(:index(trim(charint),'-')-1),*,iostat=iii) k
              if (iii.NE.0) call quit_program_error('PROBLEM while a number was read',1,nocolor)
              read (charint(index(trim(charint),'-')+1:len(charint)),*,iostat=iii) n
              if (iii.NE.0) call quit_program_error('PROBLEM while a number was read -2-',1,nocolor)
              if (debug) PRINT *,debuglabel(:ilendbglbl),'range: ',k,n
              iato = iato+n-k
           else
              read (charint(1:7),*,iostat=iii) k
              if (iii.NE.0) call quit_program_error('PROBLEM while a number was read -3-',1,nocolor)
              if (debug) PRINT *,debuglabel(:ilendbglbl),k
           endif !! (index(trim(charint),'-').GT.0) then
           charint(:) = ''
        endif !! (filex.OR.i.EQ.len_trim(text)) then
     enddo !! i = 1,len_trim(text)
     ifrag = ifrag+1
     iato  =  iato+1
  end subroutine hmnumb

  subroutine procnumb(text,kk,imos,ifrmo,moc)
     use SUFR_system  ,only         : quit_program_error
     implicit none
     integer          ,intent(in ) :: imos,ifrmo,kk
     integer          ,intent(out) :: moc(imos)
     integer                       :: i,j,imo,iii,ntemp,nn,l
     character(len=20)             :: charint
     character(len=kk)             :: text

     text(kk:kk) = ';'
     imo = 0
     j   = 1
     do i = 1,ifrmo
        if (debug) PRINT *,debuglabel(:ilendbglbl),'j=',j,',j+1=',index(text(j:kk),';')-2+j
        read (text(j:index(text(j:kk),';')-2+j),*,iostat=iii) charint
        if (debug) PRINT *,debuglabel(:ilendbglbl),charint
        if (iii.NE.0) call quit_program_error("PROBLEM1 while a number was read",1,nocolor)
        if (index(trim(charint),'-').GT.0) then
           imo = imo+1
           read (charint(1:index(trim(charint),'-')-1),*,iostat=iii) moc(imo)
           if (debug) PRINT *,debuglabel(:ilendbglbl),imo,moc(imo)
           if (iii.NE.0) call quit_program_error("PROBLEM2 while a number was read",1,nocolor)
           read (charint(index(trim(charint),'-')+1:20),*,iostat=iii) ntemp
           if (iii.NE.0) call quit_program_error("PROBLEM3 while a number was read",1,nocolor)
           nn = imo
           do l = 1,ntemp-moc(nn)
              imo = imo+1
              moc(imo) = moc(nn)+l
              if (debug) PRINT *,debuglabel(:ilendbglbl),imo,moc(imo),' gen'
           enddo !! l = 1,ntemp-iatomat(n)
        else
           imo = imo+1
           read (charint(1:20),*,iostat=iii)     moc(imo)
           if (debug) PRINT *,debuglabel(:ilendbglbl),imo,moc(imo)
        endif !! (index(trim(charint),'-').GT.0) then
        j = index(text(j:kk),';')+j
     enddo !! i = 1,ifrmo
     if (debug) PRINT *,debuglabel(:ilendbglbl),imo,'=imo imos=',imos
     if (imo.NE.imos) call quit_program_error("Wrong number of values",1,nocolor)
  end subroutine procnumb

  subroutine     print_presentation(fchkname,g09name,outfile)

     use SUFR_text,only               : int2str,uppercase,dbl2str

     implicit none
     character(len  =*  ),intent(in) :: fchkname,g09name,outfile
     character(len  =100)            :: label,label2,label3
     integer                         :: i
      label(:) = ''
     label2(:) = ''
     label3(:) = ''
     label     = "* Computation of FIELD-INDUCED ORBITALS - Intel(R) MKL version *"
     i         = len_trim(label)
     label2    = '(6X,'//int2str(i)//'("*"))'
     label3    = '(6X,"*",'//int2str(i-2)//'X,"*")'
     write(iout,'(/)')
     write(iout,label2)
     write(iout,label3)
     write(iout,'(6X,A)') trim(label)
     label     = "*            and ELECTRON-DEFORMATION ORBITALS                 *"
     write(iout,'(6X,A)') trim(label)
     label3(:) = ''
     label3    = '(6X,"*",25X,"v.",A,'//int2str(i-35)//'X,"*",/,6X,"*",62X,"*")'
     write(iout,label3) version
     label3(:) = ''
      label(:) = ''
      label(:) = 'AUTHORS: Nicolas Otero Martinez'
     label3    = '(6X,"*",1X,A,'//int2str(i-3-len_trim(label))//'X,"*")'
     write(iout,label3) trim(label)
     label3(:) = ''
      label(:) = ''
      label(:) = 'Marcos Mandado Alonso'
     label3    = '(6X,"*",10X,A,'//int2str(i-12-len_trim(label))//'X,"*")'
     write(iout,label3) trim(label)
     label3(:) = ''
      label(:) = ''
      label(:) = 'LICENSE: GPLv3'
     label3    = '(6X,"*",1X,A,'//int2str(i-3-len_trim(label))//'X,"*")'
     write(iout,label3) trim(label)
     label3(:) = ''
      label(:) = ''
      label(:) = 'REFs   : FIOs - Phys. Chem. Chem. Phys., 2019,21, 6274-6286.'
     label3    = '(6X,"*",1X,A,'//int2str(i-3-len_trim(label))//'X,"*")'
     write(iout,label3) trim(label)
     label3(:) = ''
      label(:) = ''
      label(:) = '         EDOs - J. Comput. Chem., 2014, 35, 1261-1269.'
     label3    = '(6X,"*",1X,A,'//int2str(i-3-len_trim(label))//'X,"*")'
     write(iout,label3) trim(label)
     label3(:) = ''
      label(:) = ''
      label(:) = '         Cond.- Int. J. Quantum Chem., 2018, 118, e25651.'
     label3    = '(6X,"*",1X,A,'//int2str(i-3-len_trim(label))//'X,"*")'
     write(iout,label3) trim(label)
     label3(:) = ''
     label3    = '(6X,"*",'//int2str(i-2)//'X,"*")'
     write(iout,label3)
     label3(:) = ''
      label(:) = ''
      label(:) = 'Parts of the program BASED ON/USING:'
     label3    = '(6X,"*",1X,A,'//int2str(i-3-len_trim(label))//'X,"*")'
     write(iout,label3) trim(label)
     label3(:) = ''
      label(:) = ''
      label(:) = 'Computation of Multipole Matrices (Multiwfn, MIT):'
     label3    = '(6X,"*",3X,"(o)",1X,A,'//int2str(i-9-len_trim(label))//'X,"*")'
     write(iout,label3) trim(label)
     label3(:) = ''
      label(:) = ''
      label(:) = 'Tian Lu, Feiwu Chen, JCC, 2012, 33, 580-592'
     label3    = '(6X,"*",10X,A,'//int2str(i-12-len_trim(label))//'X,"*")'
     write(iout,label3) trim(label)
     label3(:) = ''
      label(:) = ''
      label(:) = 'libSUFR: "Some Useful FORTRAN Routines" (GPL3)'
     label3    = '(6X,"*",3X,"(o)",1X,A,'//int2str(i-9-len_trim(label))//'X,"*")'
     write(iout,label3) trim(label)
     label3(:) = ''
      label(:) = ''
      label(:) = 'http://libsufr.sourceforge.net/'
     label3    = '(6X,"*",10X,A,'//int2str(i-12-len_trim(label))//'X,"*")'
     write(iout,label3) trim(label)
     label3(:) = ''
      label(:) = ''
      label(:) = 'Parallel Quicksort routines: JAMS Fortran library'
     label3    = '(6X,"*",3X,"(o)",1X,A,'//int2str(i-9-len_trim(label))//'X,"*")'
     write(iout,label3) trim(label)
     label3(:) = ''
      label(:) = ''
      label(:) = 'https://github.com/mcuntz/jams_fortran (MIT)'
     label3    = '(6X,"*",10X,A,'//int2str(i-12-len_trim(label))//'X,"*")'
     write(iout,label3) trim(label)
     label3(:) = ''
     label3    = '(6X,"*",'//int2str(i-2)//'X,"*")'
     write(iout,label3)
     write(iout,label2)
     write(iout,'(/,3X,">> COMPUTATIONAL DETAILS:")')
     write(iout,'(6X,"o",1X,2(2(A),3X)                   )') &
                       "# threads=",int2str(nproc),"# threads for MKL=",int2str(nprocmkl)
     write(iout,'(6X,"o",1X,A,":",1X,A,1X,"(",A,")")') "fchk file name"   ,trim(fchkname),trim(extension(1))
     if (loptions(-4)) then
        write(iout,'(6X,"o",1X,A,":",1X,A,1X,"(",A,")")') "2nd fchk file name (perturbed wf)" &
                                                         ,trim(g09name) ,trim(extension(1))
     else
        write(iout,'(6X,"o",1X,A,":",1X,A,1X,"(",A,")")') "G09 log file name",trim(g09name) ,trim(extension(2))
     endif !! (loptions(-4)) then
     write(iout,'(6X,"o",1X,A,":",1X,A,1X,"(",A,")")') "Output file name" ,trim(outfile) ,trim(extension(3))
     write(iout,'(6X,"o",1X,2(2(A),3X),2(/,8X,2(2(A),3X)))') &
                       "# atoms=",int2str(nato),"# electrons=",int2str(nel),&
                       "# basis functions=",int2str(nfb),"# MOs=",int2str(nmo), &
                       "# heavy atoms=",int2str(nheavy),"# shells=",int2str(nsh)
     if (debug) PRINT *,debuglabel(:ilendbglbl),'loptions=',loptions
     if    (loptions(-1)) write(iout,'(6X,"o",1X,A,":",1X,A)'      ) trim(textopts(-1)),dbl2str(Field,4)
     if    (loptions(-2)) write(iout,'(6X,"o",1X,A,":",1X,9(A,1X))') trim(textopts(-2)),comptext(lop(:nop))
     if    (loptions(-3)) then
           call prepnumb(size(momap),momap,label)
           write(iout,'(6X,"o",1X,A,":",1X,A)'      ) trim(textopts(-3)),trim(label)
     endif    !! (loptions(-3)) then
     if    (loptions(-4)) then
           if (loptions(-5)) then
              write(iout,'(6X,"o",1X,A         )'      ) trim(textopts(-5))
           else
              write(iout,'(6X,"o",1X,A         )'      ) trim(textopts(-4))
           endif !! (loptions(-5)) then
     else
           write(iout,'(6X,"o",1X,A         )'      ) "Computation of FIOs"
     endif !! (loptions(-4)) then
     if    (loptions(-6)) write(iout,'(6X,"o",1X,A,":",1X,A,"%")'      ) &
                                trim(textopts(-6)),dbl2str(sqrccutoff*100._dbl,2)
     do i = 0,size(textopts(0:))
        if (loptions( i)) write(iout,'(6X,"o",1X,A)') trim(textopts(i))
     enddo !! i = 0,3
  end subroutine print_presentation

  subroutine     print_message(code,text,iprint)

     use SUFR_system        ,only  : quit_program_error

     implicit none
     character(len=* ),intent(in) :: text
     character(len= 1),intent(in) :: code
     character(len= 2)            :: bicode,color
     character(len=11)            :: colcode
     integer,optional,intent(in)  :: iprint
     integer                      :: iprintdef

     iprintdef = 0
     if (present(iprint)) then
        if (iprint.NE.0.OR.iprint.NE.5.OR.iprint.NE.6) iprintdef = iprint
     endif !! (present(iprint)) then
     bicode(:) = ''
     select case(code)
       case('d')
            bicode = 'DG'
            color  = '90'
       case('i')
            bicode = 'II'
            color  = '32'
       case('t')
            bicode = 'TT'
            color  = '34'
       case default
          call quit_program_error('Code to be printed is UNKNOWN',1,nocolor)
     end select
     if (.NOT.loptions(5)) then
        if (.NOT.loptions(6)) then
           colcode = achar(27)//'['//trim(color)//'m'//trim(bicode)//achar(27)//'[0m'
           write(*,'(3X,"(",A ,")",1X,A)') trim(colcode),trim(text)
        else
           write(*,'(3X,"(",A2,")",1X,A)') trim(bicode ),trim(text)
        endif !! (.NOT.loptions(6)) then
     endif !! (.NOT.loptions(5)) then
     if (iprintdef.NE.0) write(iprintdef,'(3X,"(",A2,")",1X,A)') trim(bicode),trim(text)

  end subroutine print_message

  subroutine prepnumb(imos,moc,nome)

     implicit none
     integer         ,intent(in ) :: imos
     character(len=7)             :: charint,charintd
     character(len=*),intent(out) :: nome
     integer         ,intent(in ) :: moc(imos)
     integer                      :: i,j,n

     if (debug) PRINT *,debuglabel(:ilendbglbl),'In prepnumb:',imos,moc

     nome(:) = ''
     j = 0
     n = 0
     do i = 1,imos
        if (moc(i)-j.NE.n) then
           if (debug) PRINT *,debuglabel(:ilendbglbl),moc(i),moc(i)-j,n
           n = moc(i)
           write(charint ,'(I7)') n
           if (i.NE.1) then
              if (j.EQ.1) then
                 nome = nome(1:len(trim(nome)))//';'//trim(charint(verify(charint ,' '):7))
              else
                 write(charintd,'(I7)') moc(i-1)
                 nome = nome(1:len(trim(nome)))//'-'//trim(charintd(verify(charintd,' '):7)) &
                                               //';'//trim(charint (verify(charint ,' '):7))
              endif !! (j.EQ.1) then
           else
                 nome = nome(1:len(trim(nome)))//trim(charint(verify(charint,' '):7))                      
           endif !! (i.NE.1) then
                       j = 1
           if (debug) print *,trim(nome)
        else
           j = j+1
           if (i.EQ.imos) then
                 write(charint ,'(I7)') moc(i)
                 nome = nome(1:len(trim(nome)))//'-'//trim(charint (verify(charint ,' '):7))
           endif !! (i.EQ.imos) then
        endif !! (moc(i)-j.NE.n) then
     enddo !! i = 1,imos
     if (debug) print *,'numbers= "',trim(nome),'"'
  end subroutine prepnumb

  subroutine     mymemunits(val    ,iunits,unitchgdec,unitchgbin)
     use SUFR_system         ,only: warn
     implicit none
     integer(kind=     4),intent(out) :: iunits
     real   (kind=double),intent(out) :: unitchgdec,unitchgbin
     real   (kind=double),intent(in ) :: val
     select case (floor(log10(val)))
        case ( :2)
           iunits      = 1
           unitchgdec  = 1._dbl
           unitchgbin  = 1._dbl
        case (3:5)
           iunits      = 2
           unitchgdec  = 1.D03
           unitchgbin  = 1024._dbl**(1._dbl)
        case (6:8)
           iunits      = 3
           unitchgdec  = 1.D06
           unitchgbin  = 1024._dbl**(2._dbl)
        case (9:11)
           iunits      = 4
           unitchgdec  = 1.D09
           unitchgbin  = 1024._dbl**(3._dbl)
        case (12:14)
           iunits      = 5
           unitchgdec  = 1.D12
           unitchgbin  = 1024._dbl**(4._dbl)
        case default
           iunits      = 1
           unitchgdec  = 1._dbl
           unitchgbin  = 1._dbl
           call warn('Suitable units not detected to print memory estimation',nocolor=nocolor)
     end select
  end subroutine mymemunits

  subroutine     mymemory
     use SUFR_text                           ,only: dbl2str
     implicit none
     integer  (kind=4)                           :: iunits=1
     real     (kind=8)                           :: unitchgdec,unitchgbin &  !! Units conversion.
                                                  , rmemwfn  = 0._dbl  &  !! wfn(basis set, nuclear coordinates).
                                                  , rmempsp  = 0._dbl  &  !! P=SP basis set coeffs if available.
                                                  , rmemmoc  = 0._dbl  &  !! MOs coeffs. + Inv. mat. of MOs c.
                                                  , rmemmo2  = 0._dbl  &  !! (EDO) MOs coeffs. + Inv. mat. of MOs c.
                                                  , rmemcond = 0._dbl  &  !! Conductance.
                                                  , rmemdip  = 0._dbl  &  !! Dipole matrices.
                                                  , rmemdbf  = 0._dbl  &  !! BF DM derivatives array.
                                                  , rmemdmo  = 0._dbl  &  !! MO DM derivatives array.
                                                  , rmemfio  = 0._dbl  &  !! FIOs computation.
                                                  , rmemedo  = 0._dbl  &  !! EDOs computation.
                                                  , rmemcomm           &  !! Common memory.
                                                  , rmemfioedo         &  !! FIOs or EDOs, DBF, DMO
                                                  , rmemtot           !&  !! Total memory.
     character( len=   2),dimension(5),parameter :: unitsdec=(/'B ' ,'kB' ,'MB' ,'GB' ,'TB' /)
     character( len=   3),dimension(5),parameter :: unitsbin=(/'B  ','KiB','MiB','GiB','TiB'/)
     character( len=1000)                        :: label
     if (.NOT.read_mult_mat) then  !! Compute Multipole matrices
        rmemwfn = rmemwfn + 8._dbl*dble( &
                                            +nato                 & !!            x(nato                 )
                                            +nato                 & !!            y(nato                 )
                                            +nato                 & !!            z(nato                 )
                                       )
     endif !! (.NOT.read_mult_mat) then  !! Compute Multipole matrices
        rmemwfn = rmemwfn + 4._dbl*dble( &
                                            +nsh                  & !!         itsh(nsh                  )
                                            +nsh                  & !!         icsh(nsh                  )
                                            +nato*3               & !!    ato2basis(nato                 )*3
                                       )
        rmemmoc = rmemmoc + 8._dbl*dble( &
                                            +nfb*nmo              & !!            c( nfb    ,nmo    ,1   )
                                            +nmo*nfb              & !!           ci(nmo     ,nfb    ,1   )
                                       )
     if (.NOT.read_mult_mat) then !! Compute Multipole matrices
        rmemwfn = rmemwfn + 4._dbl*dble( &
                                            +nsh                 & !!        inpps(nsh                  )
                                       )
        rmemwfn = rmemwfn + 8._dbl*dble( &
                                            +nprimsh             & !!    concoeffs(nprimsh              )
                                            +nprimsh             & !!     primexps(nprimsh              )
                            )
        rmempsp = rmempsp + 8._dbl*dble( &
                                            +nprimsh             & !! pspconcoeffs(nprimsh              )
                            )
     endif !! (.NOT.read_mult_mat) then !! Compute Multipole matrices
        rmemdip = rmemdip + 8._dbl*dble( &
                                            +nfb*nfb*3           & !!          dip(nfb     ,nfb     ,3  )
                                            +nmo*nmo*3           & !!        dipmo(nmo     ,nmo     ,3  )
                                            +nmotrunc*nmotrunc*3 & !!   dipmotrunc(nmotrunc,nmotrunc,3  )
                            )
     if (spec_comps.OR.edo) then !! After specifying the components to be calculated
        rmemdbf = rmemdbf + 8._dbl*dble( &
                                            +nfb*nfb*nop         & !!            d(nfb     ,nfb     ,nop)
                            )
        if (.NOT.skip_fio_fchk.AND.edo) then !! Keep in mind we need to calculate the DM for the fchk file
        rmemdbf = rmemdbf + 8._dbl*dble( &
                                            +nfb*nfb             & !!            d(nfb     ,nfb     ,3  )
                            )
        endif !! (.NOT.skip_fio_fchk.AND.edo) then
        rmemdmo = rmemdmo + 8._dbl*dble( &
                                            +nmo*nmo*nop         & !!          dmo(nmo     ,nmo     ,nop)
                            )
     else
        rmemdbf = rmemdbf + 8._dbl*dble( &
                                            +nfb*nfb*1           & !!            d(nfb     ,nfb     ,nop)
                            )
        rmemdmo = rmemdmo + 8._dbl*dble( &
                                            +nmo*nmo*1           & !!          dmo(nmo     ,nmo     ,nop)
                            )
     endif !! (spec_comps) then !! After specifying the components to be calculated
     if (.NOT.edo) then !! Allocation for FIOs
        rmemfio = rmemfio + 8._dbl*dble( &
                                            +nmotrunc*nmotrunc   & !!            A(nmotrunc,nmotrunc    )
                                            +nmotrunc*nmotrunc   & !!            B(nmotrunc,nmotrunc    )
                                            +nmotrunc            & !!           D2(nmotrunc             )
                                            +nfb*nmotrunc        & !!         cedo(nfb     ,nmotrunc    )
                                            +nmotrunc            & !!          wMO(nmotrunc             )
                                            +nmotrunc            & !!          arr(nmotrunc             )
                                            +nmotrunc*3          & !!    polweight(nmotrunc,3           )
                            )
        rmemfio = rmemfio + dble(selected_char_kind ('ISO_10646'))*( &
                                            +nmotrunc            & !!       MOstat(nmotrunc             )
                                                                 )
        rmemfio = rmemfio + 4._dbl*dble( &
                                            +nmotrunc            & !!   ordD2Acedo(nmotrunc             )
                                            +nmotrunc            & !!         indx(nmotrunc             )
                            )
     else
        rmemmo2 = rmemmo2 + 8._dbl*dble( &
                                            +nfb*nmo              & !!           c( nfb    ,nmo     ,2  )
                                            +nmo*nfb              & !!          ci(nmo     ,nfb     ,2  )
                                       )
        rmemedo = rmemmo2
     endif !! (.NOT.edo) then !! Allocation for FIOs
     if (lconduc) then
        rmemcond = rmemcond + 8._dbl*dble( &
                                            +nfb*nfb             & !!         over(nfb     ,nfb         )
                              )
     endif !! (lconduc) then
     rmemtot      = rmemwfn+rmempsp+rmemmoc+rmemcond+rmemdip+rmemdbf+rmemdmo+rmemfio+rmemedo
     rmemcomm     = rmemwfn+rmempsp+rmemmoc         +rmemdip
     rmemfioedo   =                                          rmemdbf+rmemdmo+rmemfio+rmemedo
     write(iout,'(/,3X,">> MEMORY ESTIMATION:")')
     call        mymemunits(rmemtot,iunits,unitchgdec,unitchgbin)
     write(iout,'(6X,"|->",1X,A                   )',advance='no') "Total (Common"
     if (edo) then
        write(iout,'(          A                   )',advance='no') "+EDOs"
        if (lconduc) write(iout,'(A                )',advance='no') "+Conductance"
     else
        write(iout,'(          A                   )',advance='no') "+FIOs"
     endif !! (edo) then
     write(iout,'(          A,":",1X,A,A,"/",A,A)') ")" &
                                                           ,dbl2str(rmemtot/unitchgdec,2),unitsdec(iunits) &
                                                           ,dbl2str(rmemtot/unitchgbin,2),unitsbin(iunits)
     label(:) = ''
     label    = 'Memory estimation: '//dbl2str(rmemtot/unitchgdec,2)//unitsdec(iunits)//"/"// &
                dbl2str(rmemtot/unitchgbin,2)//unitsbin(iunits)
     call print_message('i',trim(label))

     call        mymemunits(rmemfioedo,iunits,unitchgdec,unitchgbin)
     if (edo) then
        write(iout,'(6X,"|--->",1X,A,":",1X,A,A,"/",A,A,1X,"(",A,"%)")') &
                                               "EDOs (Computation+DBF+DMO)" &
                                                           ,dbl2str(rmemfioedo/unitchgdec,2),unitsdec(iunits) &
                                                           ,dbl2str(rmemfioedo/unitchgbin,2),unitsbin(iunits) &
                                                           ,dbl2str(100._dbl*rmemfioedo/rmemtot,2)
        call     mymemunits(rmemedo,iunits,unitchgdec,unitchgbin)
        write(iout,'(6X,"|  |-->",1X,A,":",1X,A,A,"/",A,A,1X,"(",A,"%)")') &
                                               "EDOs (Computation)" &
                                                           ,dbl2str(rmemedo/unitchgdec,2),unitsdec(iunits) &
                                                           ,dbl2str(rmemedo/unitchgbin,2),unitsbin(iunits) &
                                                           ,dbl2str(100._dbl*rmemedo/rmemtot,2)
        call     mymemunits(rmemdbf,iunits,unitchgdec,unitchgbin)
        write(iout,'(6X,"|  |-->",1X,A,":",1X,A,A,"/",A,A,1X,"(",A,"%)")') &
                                               "Def.mat. expressed as BF (DBF)" &
                                                           ,dbl2str(rmemdbf/unitchgdec,2),unitsdec(iunits) &
                                                           ,dbl2str(rmemdbf/unitchgbin,2),unitsbin(iunits) &
                                                           ,dbl2str(100._dbl*rmemdbf/rmemtot,2)
        call     mymemunits(rmemdmo,iunits,unitchgdec,unitchgbin)
        write(iout,'(6X,"|  |-->",1X,A,":",1X,A,A,"/",A,A,1X,"(",A,"%)")') &
                                               "Def.mat. expressed as MO (DMO)" &
                                                           ,dbl2str(rmemdmo/unitchgdec,2),unitsdec(iunits) &
                                                           ,dbl2str(rmemdmo/unitchgbin,2),unitsbin(iunits) &
                                                           ,dbl2str(100._dbl*rmemdmo/rmemtot,2)
     else
        write(iout,'(6X,"|--->",1X,A,":",1X,A,A,"/",A,A,1X,"(",A,"%)")') &
                                               "FIOs (Computation+DBF+DMO)" &
                                                           ,dbl2str(rmemfioedo/unitchgdec,2),unitsdec(iunits) &
                                                           ,dbl2str(rmemfioedo/unitchgbin,2),unitsbin(iunits) &
                                                           ,dbl2str(100._dbl*rmemfioedo/rmemtot,2)
        call     mymemunits(rmemfio,iunits,unitchgdec,unitchgbin)
        write(iout,'(6X,"|  |-->",1X,A,":",1X,A,A,"/",A,A,1X,"(",A,"%)")') &
                                               "FIOs (Computation)" &
                                                           ,dbl2str(rmemfio/unitchgdec,2),unitsdec(iunits) &
                                                           ,dbl2str(rmemfio/unitchgbin,2),unitsbin(iunits) &
                                                           ,dbl2str(100._dbl*rmemfio/rmemtot,2)
        call     mymemunits(rmemdbf,iunits,unitchgdec,unitchgbin)
        write(iout,'(6X,"|  |-->",1X,A,":",1X,A,A,"/",A,A,1X,"(",A,"%)")') &
                                               "Der.mat. expressed as BF (DBF)" &
                                                           ,dbl2str(rmemdbf/unitchgdec,2),unitsdec(iunits) &
                                                           ,dbl2str(rmemdbf/unitchgbin,2),unitsbin(iunits) &
                                                           ,dbl2str(100._dbl*rmemdbf/rmemtot,2)
        call     mymemunits(rmemdmo,iunits,unitchgdec,unitchgbin)
        write(iout,'(6X,"|  \-->",1X,A,":",1X,A,A,"/",A,A,1X,"(",A,"%)")') &
                                               "Der.mat. expressed as MO (DMO)" &
                                                           ,dbl2str(rmemdmo/unitchgdec,2),unitsdec(iunits) &
                                                           ,dbl2str(rmemdmo/unitchgbin,2),unitsbin(iunits) &
                                                           ,dbl2str(100._dbl*rmemdmo/rmemtot,2)
     endif !! then (edo)
     if (.NOT.spec_comps.AND..NOT.edo) write(iout,'(6X,"|",/,&
                                      &6X,"|  (NOTE: the size of derivatives matrices array is impossible",/ &
                                      &6X,"|  to know beforehand without reading Gaussian log file.",/        &
                                      &6X,"|  Real memory use will be printed later)",/,6X,"|")')
     call        mymemunits(rmemcomm,iunits,unitchgdec,unitchgbin)
     write(iout,'(6X,"|--->",1X,A,":",1X,A,A,"/",A,A,1X,"(",A,"%)")') &
                                               "Common (Basis set+dipole mats.+MOs coeffs+Inv MOs coeffs)" &
                                                          ,dbl2str(rmemcomm/unitchgdec,2),unitsdec(iunits) &
                                                          ,dbl2str(rmemcomm/unitchgbin,2),unitsbin(iunits) &
                                                          ,dbl2str(100._dbl*rmemcomm/rmemtot,2)
     call        mymemunits(rmemwfn+rmempsp,iunits,unitchgdec,unitchgbin)
     write   (iout,'(9X,   "|-->",1X,A,":",1X,A,A,"/",A,A,1X,"(",A,"%)")') &
                                               "Basis set" &
                                                     ,dbl2str((rmemwfn+rmempsp)/unitchgdec,2),unitsdec(iunits) &
                                                     ,dbl2str((rmemwfn+rmempsp)/unitchgbin,2),unitsbin(iunits) &
                                                     ,dbl2str(100._dbl*(rmemwfn+rmempsp)/rmemtot,2)
     call        mymemunits(rmemdip,iunits,unitchgdec,unitchgbin)
     write   (iout,'(6X,"   |-->",1X,A,":",1X,A,A,"/",A,A,1X,"(",A,"%)")') &
                                               "Dip.mat. arrays" &
                                                           ,dbl2str(rmemdip/unitchgdec,2),unitsdec(iunits) &
                                                           ,dbl2str(rmemdip/unitchgbin,2),unitsbin(iunits) &
                                                           ,dbl2str(100._dbl*rmemdip/rmemtot,2)
     call        mymemunits(rmemmoc,iunits,unitchgdec,unitchgbin)
     write   (iout,'(6X,"   \-->",1X,A,":",1X,A,A,"/",A,A,1X,"(",A,"%)")') &
                                               "MOs coeffs. & Inv. MOs coeffs. arrays" &
                                                           ,dbl2str(rmemmoc/unitchgdec,2),unitsdec(iunits) &
                                                           ,dbl2str(rmemmoc/unitchgbin,2),unitsbin(iunits) &
                                                           ,dbl2str(100._dbl*rmemmoc/rmemtot,2)

     rmempart = rmemtot-rmemdbf-rmemdmo   !! To use after reading g09 log file (real dimension of DBF and DMO)

  end subroutine mymemory

  subroutine     myrealmemory
     use SUFR_text                             ,only: dbl2str
     implicit none
     integer  (kind=     4)                        :: iunits
     real     (kind=double)                        :: rmemdbf,rmemdmo,rmemtot,unitchgdec,unitchgbin
     character( len=     2),dimension(5),parameter :: unitsdec=(/'B ' ,'kB' ,'MB' ,'GB' ,'TB' /)
     character( len=     3),dimension(5),parameter :: unitsbin=(/'B  ','KiB','MiB','GiB','TiB'/)
     character( len=  1000)                        :: label
     rmemdbf = rmemdbf + 8._dbl*dble( &
                                         +nfb*nfb*nop         & !!            d(nfb     ,nfb     ,nop)
                         )
     if (.NOT.skip_fio_fchk.AND.edo) rmemdbf = rmemdbf+nfb*nfb
     rmemdmo = rmemdmo + 8._dbl*dble( &
                                         +nmo*nmo*nop         & !!          dmo(nmo     ,nmo     ,nop)
                         )
     rmemtot = rmempart + rmemdbf + rmemdmo
     write(iout,'(/,3X,">> MEMORY USE ESTIMATION (after reading input files):")')
     call        mymemunits(rmemtot,iunits,unitchgdec,unitchgbin)
     write(iout,'(6X,"o",1X,"Total:"   ,1X,A,A,"/",A,A, &
                 &1X,"-",1X,"Real DBF:",1X,A,A,"/",A,A, &
                 &1X,"-",1X,"Real DMO:",1X,A,A,"/",A,A  &
                &)') dbl2str(rmemtot/unitchgdec,2),unitsdec(iunits) &
                    ,dbl2str(rmemtot/unitchgbin,2),unitsbin(iunits) &
                    ,dbl2str(rmemdbf/unitchgdec,2),unitsdec(iunits) &
                    ,dbl2str(rmemdbf/unitchgbin,2),unitsbin(iunits) &
                    ,dbl2str(rmemdmo/unitchgdec,2),unitsdec(iunits) &
                    ,dbl2str(rmemdmo/unitchgbin,2),unitsbin(iunits)
     label(:) = ''
     label    = 'Memory use: '//dbl2str(rmemtot/unitchgdec,2)//unitsdec(iunits)//"/"// &
                                dbl2str(rmemtot/unitchgbin,2)//unitsbin(iunits)
     call print_message('i',trim(label))
  end subroutine myrealmemory

  subroutine     help_legend
    implicit none
    write(*,'(/,A,/,5X,A,/)') 'Legend for opts.:','[CM] = common; [ED] = EDOs-only; [CO] = EDOs+Cond.; [FI] = FIOs-only; &
    &** = Mandatory'
  end subroutine help_legend

end module commonmod
