! ** compute/basis2atom_map.f90 >> Mapping basis functions to atoms
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

subroutine     basis2atom_map

  use SUFR_system    ,only: warn,quit_program_error
  use SUFR_text      ,only: int2str
  use commonmod      ,only: nato,nsh,ato2basis,icsh,itsh,nbfsize,print_message &
                          , edo       & !! loptions( -4)
                          , stealth   & !! loptions(  5)
                          , nocolor  !& !! loptions(  6)

  implicit none

  integer  (kind =    4) :: i,j,k !! tmp
  character( len = 1000) :: label

!
!     Basis to Atom map
!
  if (allocated(ato2basis)) call quit_program_error('ato2basis is unexpectedly allocated',1,nocolor)
  allocate(ato2basis(nato))  !+D:basis2atom_map OR conductance_calc
                             !! This is a FORTRAN class. I don't know how to create a procedure in "myalloc" to
                             !! to allocate it there. However, the array is relatively small to represent a 
                             !! remarkable fingerprint in the memory.
  ato2basis%ifb = 0
  do i = 1,nsh
     ato2basis(icsh(i))%ifb  = ato2basis(icsh(i))%ifb &
                             +   nbfsize(itsh(i))
  enddo !! i = 1,nsh
  ato2basis%ifbstart         = 1
  ato2basis%ifbend           = ato2basis(1)%ifb
  do i = 2,nato
     ato2basis(i)%ifbstart   = ato2basis(i-1)%ifbstart &
                             + ato2basis(i-1)%ifb
     ato2basis(i)%ifbend     = ato2basis(i-1)%ifbend   &
                             + ato2basis(i  )%ifb
  enddo !! i = 2,nato
  if (.NOT.stealth) then
     if (nato.LT.30) then
        i = max(ceiling(log10(dble(abs(nato)))),1)
        j = max(ceiling(log10(dble(abs(ato2basis(nato)%ifbstart)))),1)
        k = max(ceiling(log10(dble(abs(ato2basis(nato)%ifbend  )))),1)
        label = '(8X,"- First and last basis functions of Atom ",A'// &
                int2str(i)//'," from ",A'//int2str(j)//'" to ",A'//int2str(k)//')'
        call print_message('i',"Assignation of basis functions to atoms:")
        write(*,*)
        write(*,label) &
           (int2str(i),int2str(ato2basis(i)%ifbstart),int2str(ato2basis(i)%ifbend),i=1,nato)
        write(*,*)
     else
        call warn('Number of atoms too large (>30) to be printed on screen',nocolor=nocolor)
     endif !! (nato.LT.30) then
  endif !! (.NOT.stealth) then
  if (.NOT.edo) deallocate(ato2basis) !+A:basis2atom_map

end subroutine basis2atom_map
