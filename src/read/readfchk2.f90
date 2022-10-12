! ** read/readfchk2.f90 >> Subroutine to read the 2nd fchk file (wf under an electric field)
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

subroutine     readfchk2
  use SUFR_kinds ,only                             : double,dbl
  use SUFR_system,only                             : quit_program_error
  use SUFR_text  ,only                             : int2str
  use memory_use
  use commonmod  ,only                             : igau ,nato,nel,nfb,nmo,nsh,itsh,icsh,d &
                                                   , debuglabel,ilendbglbl,ndimtriang,Field &
                                                   , edocomp,print_message,c &
                                                   , lconduc          &       !! loptions( -5)
                                                   , nocolor          &       !! loptions(  6)
                                                   , debug            &       !! loptions(  7)
                                                   , mem_info        !&       !! loptions( 14)
  use readmod    ,only                             : read_check_value,read_int_array_fchk &
                                                   , read_real_array_fchk
  use conductance,only                             : electrodes

  implicit none

  integer  ( kind =  4 )                          :: natfchk,nelfchk,nfbfchk,nmofchk,nshfchk,nfbddim &
                                                   , ifieldp,ifieldn,nfbnmo
  integer  ( kind =  4 )                          :: i,j !! tmp
  integer  ( kind =  4 ),allocatable,dimension(:) :: itshfchk2,icshfchk2
  real     ( kind =  8 ),allocatable,dimension(:) :: arrabf
  real     ( kind =  8 )                          :: Fieldvec(35)
  character(  len=   3 ),dimension(-3:3),parameter:: txtfldcmp=(/ 'f-z' &
                                                                , 'f-y' &
                                                                , 'f-x' &
                                                                , 'f0d' &
                                                                , 'f+x' &
                                                                , 'f+y' &
                                                                , 'f+z' /)
  character(  len = 80 )                          :: label
!
!                  << Checking some variables previously obtained from fchk file >>
!
  call read_check_value(igau ,5,nato,natfchk,'Number of atoms'                )
  call read_check_value(igau ,5,nel ,nelfchk,'Number of electrons'            )
  call read_check_value(igau ,6,nfb ,nfbfchk,'Number of basis functions'      )
  call read_check_value(igau ,6,nmo ,nmofchk,'Number of independent functions')
  call read_check_value(igau ,6,nsh ,nshfchk,'Number of contracted shells'    )
  call myalloc(itshfchk2,nsh,'readfchk2','itshfchk2',nocolor=nocolor,verbose=mem_info) !+D:here
  call myalloc(icshfchk2,nsh,'readfchk2','icshfchk2',nocolor=nocolor,verbose=mem_info) !+D:here
  call read_int_array_fchk(igau ,5,'Shell types'      ,nshfchk,itshfchk2(1:nshfchk))
  if (debug) PRINT *,debuglabel(:ilendbglbl),itshfchk2(1),'=itsh(1)',itshfchk2(nsh),'=itsh(nsh)'
  call read_int_array_fchk(igau ,7,'Shell to atom map',nshfchk,icshfchk2(1:nshfchk))
  if (debug) PRINT *,debuglabel(:ilendbglbl),icshfchk2(1),'=icsh(1)',icshfchk2(nsh),'=icsh(nsh)'
  if ((.NOT.all(itshfchk2.EQ.itsh)).OR.(.NOT.all(icshfchk2.EQ.icsh))) &
          call quit_program_error('The 2nd fchk file contains different basis set info',1,nocolor)
  if (.NOT.allocated(d)) call quit_program_error('DBF is NOT allocated before reading the 2nd fchk file',1,nocolor)
  call read_real_array_fchk(igau,5,'External E-field',35,Fieldvec(:35))
  if (all(Fieldvec(2:4).EQ.0._dbl)) call quit_program_error('Null electric field vector?',1,nocolor)
  Field   = sqrt(dot_product(Fieldvec(2:4),Fieldvec(2:4)))
  ifieldp = maxloc(Fieldvec(1:4),1)
  ifieldn = minloc(Fieldvec(1:4),1)
  if (Fieldvec(ifieldp).GT.0._dbl) then
     edocomp =   ifieldp-1
  else if (Fieldvec(ifieldn).LT.0._dbl) then
     edocomp = -(ifieldn-1)
  endif !! (Fieldvec(ifieldp).GT.0._dbl) then
  if (debug) PRINT *,debuglabel(:ilendbglbl),'edocomp=',edocomp
  label(:) = ''
  write(label,'(A,A)') "Electric field type: ",trim(txtfldcmp(edocomp))
  call print_message('i',trim(label))
  call read_check_value(igau,6,nfb*nmo,nfbnmo,'Alpha MO coefficients')
  read (igau,'(5E16.8)') ((c(i,j,2),i=1,nfb),j=1,nmo)
  if (debug) PRINT *,debuglabel(:ilendbglbl),c(1  ,1,2),'=cv(1,1)',c(1,nmo,2),'=cv(1,nmo)' &
                    ,c(nfb,1,2),'=cv(nfb,1)',c(nfb,nmo,2),'=cv(nfb,nmo)'
  if (allocated(arrabf)) call quit_program_error('contr.DBF is allocated unexpectly',1,nocolor)
  nfbddim = ndimtriang(nfb,.FALSE.)
  call myalloc(arrabf,nfbddim,'readfchk2','arrabf',nocolor=nocolor,verbose=mem_info) !+D:here
  call read_real_array_fchk(igau,6,'Total SCF Density',nfbddim,arrabf(:nfbddim))
  if (debug) PRINT *,debuglabel(:ilendbglbl),arrabf(1      ),'=cDBF(1,1) cDBF(1,nbf) =' &
                                            ,arrabf(ndimtriang(nfb,.TRUE.)+1)  &
                                            ,arrabf(nfbddim),'=cDBF(nfb,nfb)'

!$omp parallel shared ( nfb,d,arrabf ) private ( i,j )
  !$omp do
     do i = 1,nfb
        d(i,i,2)    = arrabf(ndimtriang(i,.FALSE.)  )
        do j = 1,i-1
           d(i,j,2) = arrabf(ndimtriang(i,.TRUE. )+j)
           d(j,i,2) = d(i,j,2)
        enddo !! j = 1,i-1
     enddo !! i = 1,nfb
  !$omp end do
!$omp end parallel

  if (debug) write(*,'(1X,A,      E13.6,A,E13.6,/,8X,E13.6,A,E13.6)')   &
                                              debuglabel(:ilendbglbl)   &
                         ,d(1  ,1,2),'=dv(1,1) dv(1,n)=',d(1  ,nfb,2) &
                         ,d(nfb,1,2),'=dv(n,1) dv(n,n)=',d(nfb,nfb,2)
  call mydealloc(arrabf,'readfchk2','arrabf',nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(itshfchk2,'readfchk2','itshfchk2',nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(icshfchk2,'readfchk2','icshfchk2',nocolor=nocolor,verbose=mem_info) !+A:here
  if (lconduc) then
     do i = 1,2
        if (any(electrodes(:electrodes(size(electrodes(:,i)),i),i).GT.nato)) &
               call quit_program_error('Wrong atom specification in electrodes array',1,nocolor)
     enddo !! i = 1,2
  endif !! (lconduc) then

end subroutine readfchk2
