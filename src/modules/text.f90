!> \file text.f90  Procedures to manipulate text/strings


!  Copyright (c) 2002-2017  AstroFloyd - astrofloyd.org
!   
!  This file is part of the libSUFR package, 
!  see: http://libsufr.sourceforge.net/
!   
!  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
!  by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!  
!  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License along with this code.  If not, see 
!  <http://www.gnu.org/licenses/>.




!***********************************************************************************************************************************
!> \brief  Procedures to manipulate text/strings

module SUFR_text
  implicit none
  save
  
contains
  
  
  !*********************************************************************************************************************************
  !> \brief  Make a string lower case
  !!
  !! \param str  String
  
  function lowercase(str)
    implicit none
    character, intent(in) :: str*(*)
    character :: lowercase*(len(str))
    integer :: i,ch
    
    lowercase = str
    do i=1,len_trim(lowercase)
       ch = ichar(lowercase(i:i))
       if(ch.ge.65.and.ch.le.91) ch = ch + 32
       lowercase(i:i) = char(ch)
    end do
    
  end function lowercase
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Make a string upper case
  !!
  !! \param str  String
  
  function uppercase(str)
    implicit none
    character, intent(in) :: str*(*)
    character :: uppercase*(len(str))
    integer :: i,ch
    
    uppercase = str
    do i=1,len_trim(uppercase)
       ch = ichar(uppercase(i:i))
       if(ch.ge.97.and.ch.le.123) ch = ch - 32
       uppercase(i:i) = char(ch)
    end do
    
  end function uppercase
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Make a string lower case with an upper-case initial
  !!
  !! \param str  String
  
  function uppercaseinitial(str)
    implicit none
    character, intent(in) :: str*(*)
    character :: uppercaseinitial*(len(str))
    integer :: i,ic
    
    uppercaseinitial = str
    
    ! Capitalise first letter:
    ic = ichar(uppercaseinitial(1:1))
    if(ic.ge.97.and.ic.le.122) uppercaseinitial(1:1) = char(ic-32)
    
    ! Make the rest of the letters lower case:
    do i=2,len_trim(uppercaseinitial)
       ic = ichar(uppercaseinitial(i:i))
       if(ic.ge.65.and.ic.le.90) uppercaseinitial(i:i) = char(ic+32)
    end do
    
  end function uppercaseinitial
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert a UTF-16 string to UTF-8
  !!
  !! \param  str16     UTF-16 string
  !! \retval UTF16to8  UTF-8 string (about half the length of str16)
  
  function UTF16to8(str16)
    implicit none
    character, intent(in) :: str16*(*)
    character :: UTF16to8*((len(str16)+1)/2)
    integer :: ic16, ic8
    
    UTF16to8(1:len(UTF16to8)) = ' '
    do ic16=1,len(str16),2
       ic8 = (ic16+1)/2
       UTF16to8(ic8:ic8) = str16(ic16:ic16)
    end do
    
  end function UTF16to8
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Search and replace occurences of a substring in a string
  !!
  !! \param string    Original string to replace in.  Trailing spaces are retained, call with string(1:len_trim(string))
  !!                  to ignore them and speed things up.
  !! \param str_srch  Search string
  !! \param str_repl  Replacement string
  
  subroutine replace_substring(string, str_srch, str_repl)
    implicit none
    character, intent(inout) :: string*(*)
    character, intent(in) :: str_srch*(*),str_repl*(*)
    integer :: lstr,is,lsrch, il, maxloop
    
    lstr  = len(string)
    lsrch = len(str_srch)
    if(lsrch.gt.lstr) return  ! Search string is longer than string
    
    is = huge(is)
    maxloop = lstr+1-lsrch  ! Prevent infinite loops
    do il = 1,maxloop
       is = index(string, str_srch, back=.false.)
       if(is.le.0) exit
       if(is.gt.maxloop) exit
       !print*,il,maxloop,lstr,is,'###'//string(max(is-5,1):min(is+5,lstr))//'###'  ! Debug output
       string = string(1:is-1)//str_repl//trim(string(is+lsrch:))
    end do
    
  end subroutine replace_substring
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Remove a substring from a string, if present
  !!
  !! \param string  String to remove the substring from.  Trailing spaces are retained, call with string(1:len_trim(string))
  !!                to ignore them and speed things up.
  !! \param substr  Substring to remove
  !! \param debug   Print debug info (T/F, optional)
  
  subroutine remove_substring(string,substr, debug)
    implicit none
    character, intent(inout) :: string*(*)
    character, intent(in) :: substr*(*)
    logical, intent(in), optional :: debug
    
    integer :: l,ls, i1, il,maxloop
    character :: tstr*(len(string))
    logical :: print_debug
    
    print_debug = .false.
    if(present(debug)) print_debug = debug
    
    ls = len(substr)     ! Length of the substring to remove
    if(ls.lt.1) return   ! Zero-length string
    
    i1 = -1
    maxloop = ceiling( real(len(string))/real(ls) )  ! Prevent infinite loops
    do il = 1,maxloop
       l = len_trim(string)
       
       i1 = index(string,substr,back=.false.)
       if(i1.le.0) exit
       
       tstr = string(1:i1-1)//string(i1+ls:l)  ! String gets shorter by ls
       
       if(print_debug) then
          print*,string(1:i1-1)
          print*,string(i1+ls:l)
          print*,string(i1:i1+ls),i1,l
          print*,trim(tstr)
       end if
       
       string = tstr
    end do
    
  end subroutine remove_substring
  !*********************************************************************************************************************************
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Verify whether a string contains any of a given list of characters
  !!
  !! \param string      String to verify
  !! \param characters  List of characters
  !!
  !! \retval string_contains_one_of  True if the string contains one of the specified characters, otherwise false
  
  function string_contains_one_of(string, characters)
    implicit none
    character, intent(in) :: string*(*), characters*(*)
    logical :: string_contains_one_of
    integer :: ci
    
    string_contains_one_of = .true.
    do ci=1,len_trim(characters)
       if(index(trim(string),characters(ci:ci)).gt.0) return  ! Match found
    end do
    
    string_contains_one_of = .false.
    
  end function string_contains_one_of
  !*********************************************************************************************************************************
  
  
  
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Print multiple tab characters
  !!
  !! \param number  Desired number of tab characters
  
  function tabs(number)
    implicit none
    integer, intent(in) :: number
    character :: tabs*(max(number,1))
    integer :: count
    
    tabs = ''
    if(number.gt.0) then
       do count=1,number
          tabs(count:count) = char(9)
       end do
    end if
    
  end function tabs
  !*********************************************************************************************************************************
  
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert an integer to a character string
  !!
  !! \param number  Integer number to convert
  
  function int2str(number)
    implicit none
    integer, intent(in) :: number
    character :: int2str*(max(ceiling(log10(dble(abs(number)+1))),1) - (sign(1,number)-1)/2)  ! 0-9 -> 1; 10-99 -> 2; +1 if <0
    
    write(int2str,'(I0)') number
    
  end function int2str
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert a double-precision real to a nice character string. Difference with the F0 format descriptor: 
  !!         replace leading . or -. with 0. and -0. respectively (0.1 iso .1; -0.1 iso -.1).
  !!
  !! \param number  Value to convert
  !! \param decim   Number of decimals to use
  !! \param mark    Decimal mark to separate the integer and fractional parts; single character, e.g. "," (optional; default: ".")
  
  function dbl2str(number, decim, mark)
    use SUFR_kinds, only: double, long
    implicit none
    real(double), intent(in) :: number
    integer, intent(in) :: decim
    character, intent(in), optional :: mark*(*)
    real(double), parameter :: eps = sqrt(epsilon(number))  ! sqrt of epsilon for a double real
    character :: fmt*(9)
    
    !> -  ceiling(log10(abs((number)))): 99 gives 2, 999 3, etc.
    !! -  + 10.d0**(-decim)/2.d0: to catch rounding up.  E.g. 99.97 with decim=1 gives ceiling(log10(abs((number)))) = 2,
    !!      but we need 3 since 100.0 must be printed - eps no longer needed?
    !! -  - (sign(1_long,floor(number,long))-1)/2: space for negative sign
    !! -  + decim: add the decimals to the total string length
    !! -  + 1:     add the decimal separator
    character :: dbl2str*(max(ceiling(log10(abs((number + 10.d0**(-decim)/2.d0) * (1.d0+eps)))),1) -  &
         (sign(1_long,floor(number,long))-1)/2 + decim + 1)
    
    write(fmt,'(A,I0,A)') '(F0.',max(decim,0),')'
    write(dbl2str, trim(fmt)) number
    
    ! Remove leading decimal points:
    if(dbl2str(1:1).eq.'.') then
       dbl2str = '0'//trim(dbl2str)
    else if(dbl2str(1:2).eq.'-.') then
       call replace_substring(dbl2str, '-.', '-0.')
    end if
    
    if(present(mark)) call replace_substring(dbl2str, '.', mark(1:1))  ! Replace default decimal point with a specified mark
    
  end function dbl2str
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert a double-precision real to a nice character string. Short alias for dbl2str().
  !!
  !! \param number  Value to convert
  !! \param decim   Number of decimals to use
  !! \param mark    Decimal mark to separate the integer and fractional parts; single character, e.g. "," (optional; default: ".")
  
  function d2s(number, decim, mark)
    use SUFR_kinds, only: double, long
    implicit none
    real(double), intent(in) :: number
    integer, intent(in) :: decim
    character, intent(in), optional :: mark*(*)
    
    real(double), parameter :: eps = sqrt(epsilon(number))  ! sqrt of epsilon for a double real
    character :: d2s*(max(ceiling(log10(abs((number + 10.d0**(-decim)/2.d0) * (1.d0+eps)))),1) -  &
         (sign(1_long,floor(number,long))-1)/2 + decim + 1)
    
    if(present(mark)) then
       d2s = dbl2str(number, decim, mark)
    else
       d2s = dbl2str(number, decim)
    end if
    
  end function d2s
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert a double-precision real to a nice character string using a comma as decimal mark.
  !!         Alias for dbl2str(number, decim, ',').
  !!
  !! \param number  Value to convert
  !! \param decim   Number of decimals to use
  
  function d2sc(number, decim)
    use SUFR_kinds, only: double, long
    implicit none
    real(double), intent(in) :: number
    integer, intent(in) :: decim
    
    real(double), parameter :: eps = sqrt(epsilon(number))  ! sqrt of epsilon for a double real
    character :: d2sc*(max(ceiling(log10(abs((number + 10.d0**(-decim)/2.d0) * (1.d0+eps)))),1) -  &
         (sign(1_long,floor(number,long))-1)/2 + decim + 1)
    
    d2sc = dbl2str(number, decim, ',')
    
  end function d2sc
  !*********************************************************************************************************************************
  
  
  !*********************************************************************************************************************************
  !> \brief  Convert a single-precision real to a nice character string.  Single-precision wrapper for dbl2str.
  !!
  !! \param number  Value to convert
  !! \param decim   Number of decimals to use
  !! \param mark    Decimal mark to separate the integer and fractional parts; single character, e.g. "," (optional; default: ".")
  
  function real2str(number, decim, mark)
    implicit none
    real, intent(in) :: number
    integer, intent(in) :: decim
    character, intent(in), optional :: mark*(*)
    character :: real2str*(max(ceiling(log10(abs(number)+sqrt(epsilon(number)))),1) - (sign(1,floor(number))-1)/2 + decim + 1)
    
    if(present(mark)) then
       real2str = dbl2str(dble(number), decim, mark)
    else
       real2str = dbl2str(dble(number), decim)
    end if
    
  end function real2str
  !*********************************************************************************************************************************
  
  
  
end module SUFR_text
!***********************************************************************************************************************************

