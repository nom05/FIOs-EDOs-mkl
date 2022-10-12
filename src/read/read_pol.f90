! ** read/read_pol.f90 >> Reading information from the problematic results lines from Gaussian log file.
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


subroutine read_pol(inp,debug,iout)

  use SUFR_system           ,only: quit_program_error
  use commonmod             ,only: nocolor  !& !! loptions(  6)
  implicit none
  integer  (kind= 4),intent(in) :: inp,iout
  integer  (kind= 4)            :: itext,itline,iok,iline,nline,ifline
  character( len=80)            :: text,lline
  logical           ,intent(in) :: debug

  interface
    subroutine       goon(ia,ib,ic,l,id)
      integer  (kind=        4),intent(in) :: ia,ib,ic,id
      logical                  ,intent(in) :: l
    end subroutine   goon
    subroutine       oldsudgfchk(ie,c,ig,ih,ij,l)
      integer          ,intent(out) :: ih,ij
      integer          ,intent(in ) :: ig,ie
      character(len= *),intent(in ) :: c
      logical          ,intent(in ) :: l
    end subroutine   oldsudgfchk
    subroutine       crabmode(iout,iext,text,icde,iline,line,debug)
      character(len= *),intent(in ) :: text
      character(len=80)             :: line
      logical          ,intent(in ) :: debug
      integer          ,intent(in ) :: iout,iext
      integer          ,intent(out) :: icde,iline
    end subroutine   crabmode
    subroutine       extprop(iext,text,prop,iprop,iline,ldline)
      integer          ,intent(out) :: iline,ldline
      integer          ,intent(in ) :: iext,iprop
      character(len= *),intent(in ) :: text,prop
    end subroutine   extprop
    integer function nlines(iout)
      integer,intent(in )           :: iout
    end function     nlines
  end interface
!
! Checking if the calculation has finished ...
!
  text = 'Normal termination o'
  itext= 20
  call oldsudgfchk(itext,text,inp,iok,nline,debug)
  if (iok.EQ.0) call quit_program_error('This calculation has not finished!!!',1,nocolor)
!
! Searching magical line using crab-mode ...
!
  text   = '@'
  itext  = 1
  call crabmode(inp,itext,text,iok,iline,lline,debug)
  if (debug) print *,'DEBUG ... N_LINE1,N_LINE1-N_LINE2: ',nline,iline  !! DEBUG
  itline = iline
!
! Searching the start of the magical line ...
!
  text   = '1\1\GINC-'
  itext  = 9
  call crabmode(inp,itext,text,iok,iline,lline,debug)
  itline = itline+iline
  ifline = iline+1
  if (debug) write(*,*) 'DEBUG ... NLINE,ILINE,NLINE-ITLINE,LLINE: ',nline,iline,nline-itline &
                        ,'+',trim(lline),'+'  !! DEBUG
!
! Rewinding and reading the magical line ...
!
  itext  = 80
  call oldsudgfchk(itext,lline,inp,iok,iline,debug)
  backspace(inp)
  if (debug) print *,'ILINE,NLINE-ITLINE,ITLINE,iok: ',iline,nline-itline,itline,iok  !! DEBUG
  if (iok.EQ.0) call quit_program_error('A problem has been found!!!',1,nocolor)
  if (nline-itline.NE.iline) call quit_program_error(' **A problem has been found!!!',1,nocolor)
  call goon(ifline,nline-itline,inp,debug,iout)

end subroutine read_pol

subroutine goon(itline,nline,inp,debug,iout)

  use commonmod              ,only: loptions
  use SUFR_kinds             ,only: double,dbl

  implicit none

  integer  (kind=        4),intent(in) :: itline,nline,inp,iout
  integer  (kind=        4)            :: icum,i,iprop,itext,iok,iline,ldline,itemp
  character( len=       70)            :: prop
  character( len=70*itline)            :: text
  character( len=       80)            :: line
  logical                  ,intent(in) :: debug
  logical                              :: field,dyn
  real     (kind=   double),parameter  :: debye2ea0=.393430307_dbl
  real     (kind=   double)            :: fieldm(4),polar(6),hyperpolar(10),dip(3),freq(2),dynpolar(6) &
                                        , dynhypolar(18), energy,atemp,betax,betay,betaz,betatot

  interface
    subroutine       oldsudgfchk(ie,c,ig,ih,ij,l)
      integer          ,intent(out) :: ih,ij
      integer          ,intent(in ) :: ig,ie
      character(len= *),intent(in ) :: c
      logical          ,intent(in ) :: l
    end subroutine   oldsudgfchk
    subroutine       extprop(iext,text,prop,iprop,iline,ldline)
      integer          ,intent(out) :: iline,ldline
      integer          ,intent(in ) :: iext,iprop
      character(len= *),intent(in ) :: text,prop
    end subroutine   extprop
  end interface

  if (debug) PRINT *,'NEW ITLINE,NLINE:',itline,nline
  text(:) = ''
  icum = 1
  do i = 1,itline
     read (inp,'(1X,A70)') line
     if (debug) write(*,'("DEBUG: ",1X,2(I8),1X,A)') icum,icum+69,line  !! DEBUG
     text(icum:icum+69) = line(1:70)
     icum=icum+70
  enddo !! i = 1,itline
  if (debug) print *,'DEBUG.CONCATENATED LINE: +',trim(text),'+'  !! DEBUG
! >>>>>>>>> FIELD <<<<<<<<<<<
  field=.FALSE.
  prop(:) = ''
  prop  = 'field=read'
  iprop = 10
  if (index(text(1:itline*70),prop(1:iprop)).GT.0) then
     if (.NOT.loptions(5)) write(*,'(9X,"** Field detected:")',advance='no')
     write(iout,'(9X,"** Field detected:")',advance='no')
     field=.TRUE.
     line = 'finite field(s) will'
     itext= 20
     call oldsudgfchk(itext,line,inp,iok,itemp,debug)
     if (iok.EQ.0) stop ' **Field detected but not found** '
     read (inp,'(33X,3(F8.4))') (fieldm(i),i=1,3)
     if (.NOT.loptions(5)) write(*,'(3(F8.4),1X)') (fieldm(i),i=1,3)
     write(iout,'(3(F8.4),1X)') (fieldm(i),i=1,3)
     fieldm(4) = sum(fieldm(:3))/3
  endif
! >>>>>>>> DYNAMIC <<<<<<<<<<
  dyn=.FALSE.
  prop(:) = ''
  prop  = 'cphf=rdfreq'
  iprop = 11
  if (index(text(1:itline*70),prop(1:iprop)).gt.0) then
     if (.NOT.loptions(5)) write(*,'(9X,"** Dynamic frequencies detected:")',advance='no')
     write(iout,'(9X,"** Dynamic frequencies detected:")',advance='no')
     dyn=.TRUE.
     line = 'urbation frequencies'
     itext= 20
     call oldsudgfchk(itext,line,inp,iok,itemp,debug)
     if (iok.EQ.0) stop ' **Frequency detected but not found** '
     backspace(inp)
     read (inp,'(32X,2(F12.6))') (freq(i),i=1,2)
     if (.NOT.loptions(5)) write(*,'(2(F12.5),1X)') (freq(i),i=1,2)
     write(iout,'(2(F12.5),1X)') (freq(i),i=1,2)
  endif
! >>>>>>>>> ENERGY <<<<<<<<<<<
  prop(:) = ''
  prop  = 'HF='
  iprop = 3
  call extprop(70*itline,text,prop,iprop,iline,ldline)
  if (debug) print *,'DEBUG.from g09 output: ',text(iline:ldline)  !! DEBUG
  read (text(iline:ldline),'(F18.5)') energy
  if (.NOT.loptions(5)) write(*,'(9X,"E(HF/KS)",3X,"= ",F15.7)') energy
  write(iout,'(9X,"E(HF/KS)",3X,"= ",F15.7)') energy
! >>>>>>>>> DIP. MOM. <<<<<<<
  prop(:) = ''
  prop    = 'Dipole='
  iprop   = 7
  call extprop(70*itline,text,prop,iprop,iline,ldline)
  if (debug) print *,'DEBUG.from g09 output: ',text(iline:ldline)  !! DEBUG
  itemp = iline
  do i = 1,2    !! <- 2 + 1 without ','
     prop  = text(itemp:itemp+index(text(itemp:ldline),',')-2)
     itemp = itemp+index(text(itemp:ldline),',')
     if (debug) print *,'DEBUG.from g09 output: ',prop  !! DEBUG
     read (prop,'(F20.7)') dip(i)
  enddo !! i = 1,2
  prop  = text(itemp:ldline)
  read (prop,'(F20.7)') dip(3)
  if (debug) print *,'DEBUG: ',prop,dip  !! DEBUG
  atemp = sqrt(dot_product(dip,dip))
  if (.NOT.loptions(5)) write(*,'(9X,"Dip. Mom.",2X,"= ",3(F10.5,1X),"(tot=",F10.5," au,",F10.5," D)")') &
                                (dip(i),i=1,3),atemp,atemp/debye2ea0
  write(iout,'(9X,"Dip. Mom.",2X,"= ",3(F10.5,1X),"(tot=",F10.5," au,",F10.5," D)")') &
                                (dip(i),i=1,3),atemp,atemp/debye2ea0
! >>>>>>>>> POLAR <<<<<<<<<<<
  prop(:) = ''
  prop    = 'Polar='
  iprop   = 6
  call extprop(70*itline,text,prop,iprop,iline,ldline)
  if (debug) print *,'DEBUG.from g09 output: ',text(iline:ldline)  !! DEBUG
  itemp = iline
  do i = 1,5    !! <- 5 + 1 without ','
     prop  = text(itemp:itemp+index(text(itemp:ldline),',')-2)
     itemp = itemp+index(text(itemp:ldline),',')
     if (debug) print *,'DEBUG.from g09 output: ',prop  !! DEBUG
     read (prop,'(F20.7)') polar(i)
  enddo !! i = 1,5
  prop  = text(itemp:ldline)
  read (prop,'(F20.7)') polar(6)
  if (debug) print *,'DEBUG: ',prop,polar  !! DEBUG
  if (dyn) then
     line = 'Alpha (input orient'
     itext= 20
     call oldsudgfchk(itext,line,inp,iok,itemp,debug)
     line = 'Alpha(-w;w) w='
     itext= 14
     call oldsudgfchk(itext,line,inp,iok,itemp,debug)
     do i = 1,3
        read (inp,*)
     enddo !! i = 1,3
     do i = 1,6
        read (inp,'(8X,D18.6)') dynpolar(i)
     enddo !! i = 1,6
     if (debug) print *,'DEBUG.dynpolar:',dynpolar
     if (.NOT.loptions(5)) write(*,'(2X,19("-"),6(A11),1X)') 'XX','XY','YY','XZ','YZ','ZZ'
     write(iout,'(2X,19("-"),6(A11),1X)') 'XX','XY','YY','XZ','YZ','ZZ'
     if (.NOT.loptions(5)) write(*,'(2X,"Static Polar",6X,"= ",6(F10.2,1X),"(tr/3=",F10.2,")")') &
                                    (polar(i),i=1,6),(polar(1)+polar(3)+polar(6))/3
     write(iout,'(2X,"Static Polar",6X,"= ",6(F10.2,1X),"(tr/3=",F10.2,")")') &
                                    (polar(i),i=1,6),(polar(1)+polar(3)+polar(6))/3
     if (.NOT.loptions(5)) write(*,'(2X,"Dynam. Polar",6X,"= ",6(F10.2,1X),"(tr/3=",F10.2,")")') &
                                    (dynpolar(i),i=1,6),(dynpolar(1)+dynpolar(3)+dynpolar(6))/3
     write(iout,'(2X,"Dynam. Polar",6X,"= ",6(F10.2,1X),"(tr/3=",F10.2,")")') &
                                    (dynpolar(i),i=1,6),(dynpolar(1)+dynpolar(3)+dynpolar(6))/3
  else
     if (.NOT.loptions(5)) write(*,'(2X,19("-"),6(A11),1X)') 'XX','XY','YY','XZ','YZ','ZZ'
     write(iout,'(2X,19("-"),6(A11),1X)') 'XX','XY','YY','XZ','YZ','ZZ'
     if (.NOT.loptions(5)) write(*,'(9X,"Polar",6X,"= ",6(F10.2,1X),"(tr/3=",F10.2,")")') &
                                    (polar(i),i=1,6),(polar(1)+polar(3)+polar(6))/3
     write(iout,'(9X,"Polar",6X,"= ",6(F10.2,1X),"(tr/3=",F10.2,")")') &
                                    (polar(i),i=1,6),(polar(1)+polar(3)+polar(6))/3
  endif !! (dyn) then
! >>>>>>>>> HYPERPOLAR <<<<<<<<<<<
  prop(:) = ''
  prop      = 'HyperPolar='
  iprop     = 11
  call extprop(70*itline,text,prop,iprop,iline,ldline)
  if (debug) print *,'DEBUG.from g09 output: ',text(iline:ldline)  !! DEBUG
  itemp = iline
  do i = 1,9    !! <- 9 + 1 without ','
     prop  = text(itemp:itemp+index(text(itemp:ldline),',')-2)
     itemp = itemp+index(text(itemp:ldline),',')
     if (debug) print *,'DEBUG.from g09 output: ',prop  !! DEBUG
     read (prop,'(F20.7)') hyperpolar(i)
  enddo !! i = 1,9
  prop  = text(itemp:ldline)
  read (prop,'(F20.7)') hyperpolar(10)
  betax   = hyperpolar(1)+hyperpolar(3)+hyperpolar( 8)
  betay   = hyperpolar(2)+hyperpolar(4)+hyperpolar( 9)
  betaz   = hyperpolar(5)+hyperpolar(7)+hyperpolar(10)
  betatot = sqrt(betax*betax + betay*betay + betaz*betaz)
  if (debug) print *,'DEBUG: ',prop,hyperpolar  !! DEBUG
  if (dyn) then
     line  = 'Beta (input orientat'
     itext = 20
     call oldsudgfchk(itext,line,inp,iok,itemp,debug)
     line  = 'Beta(-w;w,0) w='
     itext = 15
     call oldsudgfchk(itext,line,inp,iok,itemp,debug)
     do i = 1,7
        read (inp,*)
     enddo !! i = 1,3
     do i = 1,18
        read (inp,'(8X,D18.6)') dynhypolar(i)
     enddo !! i = 1,18
     if (debug) print *,'DEBUG.dynhypolar:',dynhypolar
     if (.NOT.loptions(5)) write(*,'(2X,19("-"),1X,10(A10,1X))') &
                       'XXX','XXY','XYY','YYY','XXZ','XYZ','YYZ','XZZ','YZZ','ZZZ'
     write(iout,'(2X,19("-"),1X,10(A10,1X))') 'XXX','XXY','XYY','YYY','XXZ','XYZ','YYZ','XZZ','YZZ','ZZZ'
     if (.NOT.loptions(5)) write(*,'(2X,"Static HyperPolar = ",10(F10.2,1X),1X,"(tot=",F10.2,")")') &
                                          (hyperpolar(i),i=1,10),betatot
     write(iout,'(2X,"Static HyperPolar = ",10(F10.2,1X),1X,"(tot=",F10.2,")")') &
                                          (hyperpolar(i),i=1,10),betatot
     if (.NOT.loptions(5)) write(*,'(22X,9(A10,1X))') 'XXX','YXX','YYX','ZXX','ZYX','ZZX','XXY','YXY','YYY'
     write(iout,'(22X,9(A10,1X))') 'XXX','YXX','YYX','ZXX','ZYX','ZZX','XXY','YXY','YYY'
     if (.NOT.loptions(5)) write(*,'(2X,"Dynam. HyperPolar = ",9(F10.2,1X))') (dynhypolar(i),i=1,9)
     if (.NOT.loptions(5)) write(*,'(22X,9(F10.2,1X))') (dynhypolar(i),i=10,18)
     write(iout,'(2X,"Dynam. HyperPolar = ",9(F10.2,1X))') (dynhypolar(i),i=1,9)
     write(iout,'(22X,9(F10.2,1X))') (dynhypolar(i),i=10,18)
     if (.NOT.loptions(5)) write(*,'(22X, 9(A10,1X))') 'ZXY','ZYY','ZZY','XXZ','YXZ','YYZ','ZXZ','ZYZ','ZZZ'
     write(iout,'(22X, 9(A10,1X))') 'ZXY','ZYY','ZZY','XXZ','YXZ','YYZ','ZXZ','ZYZ','ZZZ'
     if (field) then
        if (.NOT.loptions(5)) write(*,'(2X,"Static HyPol/field= ",10(F10.2,1X))') (hyperpolar(i)/fieldm(4),i=1,10)
        write(iout,'(2X,"Static HyPol/field= ",10(F10.2,1X))') (hyperpolar(i)/fieldm(4),i=1,10)
     endif !! (field) then
  else
     if (.NOT.loptions(5)) write(*,'(2X,19X,1X,10(A10,1X))') &
                       'XXX','XXY','XYY','YYY','XXZ','XYZ','YYZ','XZZ','YZZ','ZZZ'
     write(iout,'(2X,19("-"),1X,10(A10,1X))') 'XXX','XXY','XYY','YYY','XXZ','XYZ','YYZ','XZZ','YZZ','ZZZ'
     if (.NOT.loptions(5)) write(*,'(9X,"HyperPolar = ",10(F10.2,1X),1X,"(tot=",F10.2,")")') &
                                          (hyperpolar(i),i=1,10),betatot
     write(iout,'(9X,"HyperPolar = ",10(F10.2,1X),1X,"(tot=",F10.2,")")') (hyperpolar(i),i=1,10),betatot
     if (field) then
        if (.NOT.loptions(5)) write(*,'(9X,"HyPol/field= ",10(F10.2,1X))') (hyperpolar(i)/fieldm(4),i=1,10)
        write(iout,'(9X,"HyPol/field= ",10(F10.2,1X))') (hyperpolar(i)/fieldm(4),i=1,10)
      endif !! (field) then
  endif !! (dyn) then

end subroutine goon

subroutine oldsudgfchk(iext,text,iout,icde,iline,debug)
  implicit none
  integer          ,intent(out) :: icde,iline
  integer          ,intent(in ) :: iout,iext
  character(len= *),intent(in ) :: text
  character(len=80)             :: line
  logical          ,intent(in ) :: debug

         icde  = 0
         iline = 1
         rewind(iout)
         if (debug) PRINT *,'text ... ','"',trim(text(:iext)),'"'
 2       continue
         read (iout,'(A)',end=999) line
         if (index(line(:80),text(1:iext)).gt.0) goto 3
         iline = iline+1
         goto 2
 3       icde = 1
         if (debug) print *,'iline,icde ',iline,icde
         return

 999     continue

end subroutine oldsudgfchk

subroutine crabmode(iout,iext,text,icde,iline,line,debug)
  implicit none
  character(len= *),intent(in ) :: text
  character(len=80)             :: line
  logical          ,intent(in ) :: debug
  integer          ,intent(in ) :: iout,iext
  integer          ,intent(out) :: icde,iline

      if (debug) print *,'** CRABMODE SUBROUTINE **'
      icde  = 0
      iline = 1
      backspace(iout)
 1    continue
      backspace(iout)
      read(iout,'(a)',end=999) line
      if (debug) print *,line !! DEBUG
      if (index(line(:80),text(:iext)).gt.0) goto 2
      iline = iline+1
      backspace(iout)
      goto 1
 2    icde = 1
      return

 999  continue

end subroutine crabmode

subroutine extprop(iext,text,prop,iprop,iline,ldline)
      implicit none
      integer         ,intent(out) :: iline,ldline
      integer         ,intent(in ) :: iext,iprop
      character(len=*),intent(in ) :: text,prop
!
! Determining the start:
!
      iline  = index(text(:iext),prop(:iprop))+iprop
!
! Determining the end:
!
      ldline = iline+index(text(iline:iext),'\')-2

end subroutine extprop

function nlines(iout)
      implicit none
      integer,intent(in ) :: iout
      integer             :: nlines
      character(len=1)    :: text

      nlines = 0
      rewind(iout)
 1    continue
      read(iout,'(a)',end=999,err=2) text
      nlines = nlines+1
      goto 1
 2    stop '** An error has happened while the file was read'
 999  continue
      rewind(iout)

end function nlines

