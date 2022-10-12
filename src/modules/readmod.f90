! ** modules/readmod.f90 >> Main module with tools to read input files
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


module readmod

  use SUFR_kinds ,only     : double !!,dbl
  use SUFR_system,only     : quit_program_error
  use SUFR_text  ,only     : int2str
  use commonmod  ,only     : debug,debuglabel,ilendbglbl,nocolor

  implicit none

  contains

  subroutine sudgfchk(inp,text,icde,irw,jjj)
  
    implicit none
    character(len=80)             :: line
    character        ,intent(in ) :: text*(*)
    integer          ,intent(in ) :: inp,irw
    integer          ,intent(out) :: icde,jjj
    integer                       :: iii
  
    if (len(text).GT.80) &
        call quit_program_error('Length of searched text is greater than the expected line length'&
                                ,1,nocolor)
    icde = 0
    jjj  = 0
    if (irw.EQ.0) rewind(inp)
    do
      read (inp,'(A)',iostat=iii) line
!     PRINT *,trim(line)
      if (iii.NE.0) then
        jjj = 1
        return
      endif
      if (index(line,text).GT.0) exit
    enddo
    icde = 1
  
  end subroutine sudgfchk
  
  integer function intread(inp,ncolumn,text)
  
    implicit none
    integer  ,intent(in) :: inp,ncolumn
    character,intent(in) :: text*(*)
    character(len=1000)  :: line,word
    integer              :: icde,irw,jjj,i
  
    irw = 1 !! rewind if = 0
    call     sudgfchk(inp,text,icde,irw,jjj)
    if (icde.NE.1.OR.jjj.EQ.1) &
       call quit_program_error('Bad news, text not found',1,nocolor)
    backspace(inp)
    read (inp,'(A)',iostat=jjj) line
    if (jjj.NE.0) call quit_program_error('Error 1 in intread',1,nocolor)
    do i = 1,ncolumn-1
       word(:) = ' '
       read (line,*,iostat=jjj) word
       if (jjj.NE.0) call quit_program_error('Error 2 in intread',1,nocolor)
       call    remove_substring(line,trim(word),.TRUE.)
    enddo !! i = 1,ncolumn-1
    read (line,*,iostat=jjj) intread
    if (jjj.NE.0) call quit_program_error('Error 3 in intread',1,nocolor)
  
  end function intread

  integer function intrdrmlch(inp ,ncolumn,nrm,text)
  
    implicit none
    integer  ,intent(in) :: inp,ncolumn,nrm
    integer              :: i
    character,intent(in) :: text*(*)
    character(len=1000)  :: extrtxt
    character(len=100)   :: label
  
    extrtxt = trim(charread(inp,ncolumn,text))
    i       = len(trim(extrtxt))-nrm
    write(label       ,'("(I",I10,")")') i
    read (extrtxt(1:i),label           ) intrdrmlch

end function intrdrmlch

function charread(inp,ncolumn,text)
  
    implicit none
    integer  ,intent(in) :: inp,ncolumn
    character,intent(in) :: text*(*)
    character(len=1000)  :: line,word,charread
    integer              :: icde,irw,jjj,i,j
  
    irw = 1 !! rewind if = 0
    call     sudgfchk(inp,text,icde,irw,jjj)
    if (icde.NE.1.OR.jjj.EQ.1) &
&      call quit_program_error('Bad news, text not found',1,nocolor)
    backspace(inp)
    read (inp,'(A)',iostat=jjj) line
    if (jjj.NE.0) call quit_program_error('Error 1 in charread',1,nocolor)
    do i = 1,ncolumn-1
       word(:) = ''
       read (line,*,iostat=jjj) word
       if (len(trim(word)).EQ.0) then
          j = 1
          do while (trim(word(:)).EQ.''.AND.j.LE.10)
             word(j:j) = trim(line(j:j))
             j = j+1
          enddo !! while ()
       endif !! (len(trim(word)).EQ.0) then
       if (jjj.NE.0) call quit_program_error('Error 2 in charread',1,nocolor)
       call    remove_substring(line,trim(word),.TRUE.)
    enddo !! i = 1,ncolumn-1
    read (line,*,iostat=jjj) charread
    if (jjj.NE.0) call quit_program_error('Error 3 in charread',1,nocolor)
  
end function charread

  real(double) function realread(ncolumn,text)
  
    implicit none
    integer  ,intent(in) :: ncolumn
    character,intent(in) :: text*(*)
    character            :: copy*(len(text))
    character(len=1000)  :: word
    integer              :: i,ierr !! tmp
  
    copy(:) = text(:)
    do i = 1,ncolumn-1
       word(:) = ''
       read (copy,*,iostat=ierr) word
       if (ierr.NE.0) call quit_program_error('Error 1 in realread',1,nocolor)
       call  remove_substring(copy,trim(word),.TRUE.)
!!  Remove first match only ------------------^^^^^^
    enddo !! i = 1,ncolumn-1
    read (copy,*,iostat=ierr) realread
    if (ierr.NE.0) call quit_program_error('Error 2 in realread',1,nocolor)
  
  end function realread

  subroutine realvec(inp,ncolumn,text,ndim,vec)
  
    implicit none
    integer  ,intent(in) :: inp,ncolumn
    character,intent(in) :: text*(*)
    character(len=1000)  :: line
    integer              :: icde,irw,jjj,i
    integer,intent(in)   :: ndim
    real(double)         :: vec(ndim)
  
    irw = 1 !! rewind if = 0
    call     sudgfchk(inp,text,icde,irw,jjj)
    if (icde.NE.1.OR.jjj.EQ.1) &
       call quit_program_error('Bad news, text not found',1,nocolor)
    backspace(inp)
    read (inp,'(A)',iostat=jjj) line
    if (jjj.NE.0) call quit_program_error('Error 1 in realvec',1,nocolor)
    do i = 1,ndim
       vec(i) = realread(ncolumn-1+i,line)
    enddo !! i = 1,ndim
  
  end subroutine realvec
  
  subroutine pmofromgauss(inp,nom,S)

    use SUFR_kinds,only : double !!,dbl

    implicit none
    integer,intent(in) :: inp
    integer,intent(in) :: nom
    integer            :: istart,irow,l,ir,ierr
    integer,parameter  :: numcol=5
    real(double)       :: S(:)

    do istart = 1,nom,numcol
       read (inp,*)
       l  = min(nom,istart+numcol-1)
       do irow = 1,nom
          ir = (irow-1)*nom
          read (inp,'(8X,5(D14.6))',iostat=ierr) S(ir+istart:ir+l)
          if (ierr.NE.0) call quit_program_error('Error 1 in pmofromgauss',1,nocolor)
       enddo !! irow = istart,nom
    enddo !! istart = 1,nom,numcol

  end subroutine pmofromgauss

  subroutine p1fromgauss(inp,ncf,S)
  
    use SUFR_kinds,only : double !!,dbl
    use commonmod ,only : ndimtriang

    implicit none
    integer,intent(in) :: inp
    integer,intent(in) :: ncf
    integer            :: istart,irow,l,ir,ierr
    integer,parameter  :: numcol=5
    real(double)       :: S(:)
  
    do istart = 1,ncf,numcol
       read (inp,*)
       do irow = istart,ncf
          l  = min(irow,istart+numcol-1)
          ir = ndimtriang(irow,.TRUE.)
          read (inp,'(7X,5(D14.6))',iostat=ierr) S(ir+istart:ir+l)
          if (ierr.NE.0) call quit_program_error('Error 1 in p1fromgauss',1,nocolor)
       enddo !! irow = istart,ncf
    enddo !! istart = 1,ncf

  end subroutine p1fromgauss

  subroutine remove_substring(string,substr,first,debug)
!  Remove a substring from a string, if present
    implicit none
    character,intent(inout)       :: string*(*)
    character,intent(in)          :: substr*(*)
    logical  ,intent(in),optional :: first,debug
    
    integer :: l,ls, i1, il,maxloop
    character :: tstr*(len(string))
    logical :: print_debug,first_loop
    
    print_debug = .FALSE.
    first_loop  = .FALSE.
    if(present(debug)) print_debug = debug
    if(present(first)) first_loop  = first
    
    ls = len(substr)     ! Length of the substring to remove
    if(ls.lt.1) return   ! Zero-length string
    
    i1 = -1
    if (first_loop) then
       maxloop = 1
    else
       maxloop = ceiling( real(len(string))/real(ls) )  ! Prevent infinite loops
    endif !! (first_loop) then
    do il = 1,maxloop
       l = len_trim(string)
       
       i1 = index(string,substr,back=.false.)
       if(i1.le.0) exit
       
       tstr = string(1:i1-1)//string(i1+ls:l)  ! String gets shorter by ls
       
       if(print_debug) then
          PRINT *,debuglabel(:ilendbglbl),string(1:i1-1)
          PRINT *,debuglabel(:ilendbglbl),string(i1+ls:l)
          PRINT *,debuglabel(:ilendbglbl),string(i1:i1+ls),i1,l
          PRINT *,debuglabel(:ilendbglbl),trim(tstr)
       end if
       
       string = tstr
    end do
    
  end subroutine remove_substring

  subroutine     read_check_value(inp,ncolumn,initial,next,text,debugval)

    implicit none
    integer,intent(in) :: inp,ncolumn
    integer            :: initial,next
    character(len=*)   :: text
    character(len=1000):: message
    logical,optional   :: debugval
    logical            :: debugvaldef

    debugvaldef = .FALSE.
    if (present(debugval)) debugvaldef = debugval
    if (initial.EQ.0) then
       initial = intread(inp,ncolumn,trim(text))
       next    = initial
       if (debugvaldef) PRINT *,debuglabel(:ilendbglbl),trim(text),initial
    else
       next    = intread(inp,ncolumn,trim(text))
    endif !! (natoms.EQ.0) then
    if (initial.NE.next) then
       message(:) = ''
       message = "Value is not kept unchanged: " &
              // trim(text)                      &
              // trim(int2str(next   ))          &
              // " != "                          &
              // trim(int2str(initial))          &
              // " (Initial value)"
       call quit_program_error(trim(message),1,nocolor)
    endif !! (initial.NE.next) then

  end subroutine read_check_value

  subroutine     read_int_array_fchk(inp,ncolumn,text,iref,array)

    implicit none
    integer            ,intent(in)      :: inp,ncolumn,iref
    integer                             :: i !! tmp
    integer            ,dimension(iref) :: array
    character(len=*   ),intent(in)      :: text
    character(len=1000)                 :: message

    i = intread(inp,ncolumn,trim(text))
    message(:) = ''
    message    = trim(text)               &
              // ' ('                     &
              // int2str(i)               &
              // ') != Reference value (' &
              // int2str(iref)            &
              // ')'
    if (i.NE.iref) call quit_program_error(trim(message),1,nocolor)
    read (inp,'(6(I12))') array(1:iref)
  end subroutine read_int_array_fchk

  subroutine     read_real_array_fchk(inp,ncolumn,text,iref,array)

    implicit none
    integer            ,intent(in)      :: inp,ncolumn,iref
    integer                             :: i !! tmp
    real(double)       ,dimension(iref) :: array
    character(len=*   ),intent(in)      :: text
    character(len=1000)                 :: message

    i = intread(inp,ncolumn,trim(text))
    message(:) = ''
    message    = trim(text)               &
              // ' ('                     &
              // int2str(i)               &
              // ') != Reference value (' &
              // int2str(iref)            &
              // ')'
    if (i.NE.iref) call quit_program_error(trim(message),1,nocolor)
    read (inp,'(5(E16.8))') array(1:iref)
  end subroutine read_real_array_fchk

end module readmod
