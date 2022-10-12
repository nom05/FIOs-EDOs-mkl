! ** read/readfchk.f90 >> Subroutine to read a generic fchk file.
!                         Despite the fact it is a generic subroutine, only necessary information is read.
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


subroutine     readfchk(alloc)

  use SUFR_kinds ,only                                : double !!,dbl
  use SUFR_system,only                                : quit_program_error
  use SUFR_text  ,only                                : int2str
  use memory_use
  use commonmod  ,only                                : ifchk,nato,nel,nfb,nmo,nsh,itsh,icsh,c,nheavy,x,y,z    &
                                                      , nprimsh,concoeffs,maxnbf,primexps,inpps,pspconcoeffs   &
                                                      , ndimtriang,d,print_message,nhomo &
                                                      , read_mult_mat                &       !! loptions(  4)
                                                      , nocolor                      &       !! loptions(  6)
                                                      , debug,debuglabel,ilendbglbl  &       !! loptions(  7)
                                                      , edo                          &       !! loptions( 10)
                                                      , mem_info                    !&       !! loptions( 14)
  use readmod    ,only                                : read_check_value,intread,read_int_array_fchk  &
                                                      , read_real_array_fchk

  implicit none

  integer                                            :: natfchk,nelfchk,nfbfchk,nmofchk,nshfchk &
                                                      , nfbnmo,ncoord,nprimshfchk
  real     ( kind = double),allocatable,dimension(:) :: arrabf
  integer  ( kind =     4 )                          :: i,j,nfbddim
  integer  ( kind =     4 ),allocatable,dimension(:) :: iznumb
  character( len  =   100 )                          :: label
  logical                  ,intent(in),optional      :: alloc
  logical                                            :: allocdef=.FALSE.

  if (present(alloc)) allocdef    = alloc

  if (debug) PRINT *,debuglabel(:ilendbglbl),'fchk file full reading?',.NOT.allocdef
  rewind(ifchk)

  if (allocdef) then
     if (debug) call print_message('d','1st execution to read some variables from fchk file')
     call    read_check_value(ifchk,5,nato   ,natfchk    ,'Number of atoms'                ,debugval=debug)
     call    read_check_value(ifchk,5,nel    ,nelfchk    ,'Number of electrons'            ,debugval=debug)
     call    read_check_value(ifchk,6,nfb    ,nfbfchk    ,'Number of basis functions'      ,debugval=debug)
     call    read_check_value(ifchk,6,nmo    ,nmofchk    ,'Number of independent functions',debugval=debug)
     nhomo   = nel/2   !! Unrestricted calc to be considered ...
     if (allocated(iznumb)) deallocate(iznumb)
     call myalloc(iznumb,nato,'readfchk','iznumb',nocolor=nocolor,verbose=mem_info) !+D:here
     call read_int_array_fchk(ifchk,5,'Atomic numbers'   ,nato,iznumb(:nato))
     do i = 1,nato
        if (iznumb(i).GT.1) nheavy = nheavy+1
     enddo !! i = 1,nato
     call mydealloc(iznumb,'readfchk','iznumb',nocolor=nocolor,verbose=mem_info) !+A:here
     call    read_check_value(ifchk,6,nsh    ,nshfchk    ,'Number of contracted shells'    ,debugval=debug)
     if (.NOT.read_mult_mat)  &  !! Compute Multipole matrices
        call read_check_value(ifchk,6,nprimsh,nprimshfchk,'Number of primitive shells'     ,debugval=debug)
     if (debug) PRINT *,debuglabel(:ilendbglbl),nato,'=nato nel=',nel,nfb,'=nfb nmo=',nmo,nsh,'=nsh'
     allocdef = .FALSE.
     return
  else
     if (nato.EQ.0.OR.nel.EQ.0.OR.nfb.EQ.0.OR.nmo.EQ.0.OR.nsh.EQ.0) &
           call quit_program_error('# electrons, # basis functions, # MOs are unread',1,nocolor)
  endif !! (alloc) then

  if (.NOT.read_mult_mat) then  !! Compute Multipole matrices
     call read_check_value(ifchk,6,nato*3,ncoord,'Current cartesian coordinates')
     call myalloc(x,nato,'readfchk','x',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:calcDbas
     call myalloc(y,nato,'readfchk','y',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:calcDbas
     call myalloc(z,nato,'readfchk','z',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:calcDbas
     read (ifchk,'(5E16.8)') (x(i),y(i),z(i),i=1,nato)
     if (debug) then
        j = 4
        if (nato.LT.4) j = nato
        PRINT *,debuglabel(:ilendbglbl),j,' 1st coordinates=',(i,x(i),y(i),z(i),i=1,j)
     endif !! (debug) then
  endif !! (.NOT.read_mult_mat) then

  call read_check_value(ifchk,6,nsh ,nshfchk,'Number of contracted shells'    )
  if (.NOT.read_mult_mat) then  !! Compute Multipole matrices
      call read_check_value(ifchk,6,nprimsh,nprimshfchk,'Number of primitive shells')
      if (debug) PRINT *,debuglabel(:ilendbglbl),nsh,'=nsh nprimsh=',nprimsh
  endif !! (.NOT.read_mult_mat) then  !! Compute Multipole matrices
  call myalloc(itsh,nsh,'readfchk','itsh',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:calcDbas
  call read_int_array_fchk(ifchk,5,'Shell types'      ,nsh ,itsh(:nsh))
  if (debug) PRINT *,debuglabel(:ilendbglbl),itsh(1),'=itsh(1)',itsh(nsh),'=itsh(nsh)'
  if (any(abs(itsh)>maxnbf)) then
     label(:) = ''
     label    = 'GTFs with angular moment > '//int2str(maxnbf)//' are unsupported'
     call quit_program_error(trim(label),1,nocolor)
  endif !! (any(abs(itsh)>maxnbf)) then
  if (.NOT.read_mult_mat) then !! Compute Multipole matrices
     call myalloc(inpps,nsh,'readfchk','inpps',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:proc_dip
     call read_int_array_fchk(ifchk,8,'Number of primitives per shell',nsh,inpps(:nsh))
     i = sum(inpps)
     if (i.NE.nprimsh) then
        label(:) = ''
        label    = 'Sum of # prims per shell('//int2str(i)//') != # primitive shells('//int2str(nprimsh)//')'
        call quit_program_error(trim(label),1,nocolor)
     endif !! (sum(inpps).NE.nprimsh) then
  endif !! (.NOT.read_mult_mat) then !! Compute Multipole matrices
  call myalloc(icsh,nsh,'readfchk','icsh',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:proc_dip
  call read_int_array_fchk(ifchk,7,'Shell to atom map',nsh ,icsh(:nsh))
  if (debug) PRINT *,debuglabel(:ilendbglbl),icsh(1),'=icsh(1)',icsh(nsh),'=icsh(nsh)'

  if (.NOT.read_mult_mat) then  !! Compute Multipole matrices
     call myalloc(primexps,nprimsh,'readfchk','primexps',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:proc_dip
     call             read_real_array_fchk(ifchk,5,'Primitive exponents',nprimsh,primexps(:nprimsh))
     if (debug) PRINT *,debuglabel(:ilendbglbl),primexps(1),'=primexps(1)',primexps(nprimsh),'=primexps(nprimsh)'
     call myalloc(concoeffs,nprimsh,'readfchk','concoeffs',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:proc_dip
     call             read_real_array_fchk(ifchk,5,'Contraction coefficients',nprimsh,concoeffs(:nprimsh))
     if (debug)PRINT *,debuglabel(:ilendbglbl),concoeffs(1),'=concoeffs(1)',concoeffs(nprimsh),'=concoeffs(nprimsh)'
     read(ifchk,"(A)") label(:80)
     if (index(trim(label),"P(S=P) Contraction coefficients").NE.0) then
        if (debug) PRINT *,debuglabel(:ilendbglbl),'"P(S=P) Contraction coeffs." detected'
        backspace(ifchk)
        call myalloc(pspconcoeffs,nprimsh,'readfchk','pspconcoeffs',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:proc_dip
        call read_real_array_fchk(ifchk,6,'P(S=P) Contraction coefficients',nprimsh,pspconcoeffs(:nprimsh))
        if (debug) PRINT *,debuglabel(:ilendbglbl),pspconcoeffs(1),'=pspconcoeffs(1)' &
                          ,pspconcoeffs(nprimsh),'=pspconcoeffs(nprimsh)'
     end if
  endif !! (.NOT.read_mult_mat) then  !! Compute Multipole matrices

  call read_check_value(ifchk,6,nfb*nmo,nfbnmo,'Alpha MO coefficients')
  i = 1
  if (edo) i = 2
  call myalloc(c,nfb,nmo,i,'readfchk','c',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:EDOs&FIOs
  read (ifchk,'(5E16.8)') ((c(i,j,1),i=1,nfb),j=1,nmo)
  if (debug) PRINT *,debuglabel(:ilendbglbl),c(1  ,1,1),'=c(1,1)',c(1,nmo,1),'=c(1,nmo)' &
                    ,c(nfb,1,1),'=c(nfb,1)',c(nfb,nmo,1),'=c(nfb,nmo)'
  if (edo) then
     if (allocated(d).OR.allocated(arrabf)) call quit_program_error('DBF is allocated unexpectly',1,nocolor)
     nfbddim = ndimtriang(nfb,.FALSE.)
     call myalloc(d,nfb,nfb,3,'readfchk','dbf',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:EDOs
     call myalloc(arrabf,nfbddim,'readfchk','arrabf',nocolor=nocolor,verbose=mem_info) !+D:here
     call read_real_array_fchk(ifchk,6,'Total SCF Density',nfbddim,arrabf(:nfbddim))
     if (debug) PRINT *,debuglabel(:ilendbglbl),arrabf(1      ),'=contr.DBF(1,1) contr.DBF(1,nbf) =' &
                                               ,arrabf(nfb    )  &
                                               ,arrabf(nfbddim),'=contr.DBF(nfb,nfb)'

!$omp parallel shared ( nfb,d,arrabf ) private ( i,j )
  !$omp do
     do i = 1,nfb
        d(i,i,1)    = arrabf(ndimtriang(i,.FALSE.)  )
        do j = 1,i-1
           d(i,j,1) = arrabf(ndimtriang(i,.TRUE. )+j)
           d(j,i,1) = d(i,j,1)
        enddo !! j = 1,nfb
     enddo !! i = 1,nfb
  !$omp end do
!$omp end parallel

     if (debug) write(*,'(1X,A,I1,1X,E13.6,A,E13.6,/,8X,E13.6,A,E13.6)') &
                                                 debuglabel(:ilendbglbl) &
                          ,i,d(1  ,1,1),'=d(1,1) d(1,n)=',d(1  ,nfb,1)   &
                            ,d(nfb,1,1),'=d(n,1) d(n,n)=',d(nfb,nfb,1)
     call mydealloc(arrabf,'readfchk','arrabf',nocolor=nocolor,verbose=mem_info) !+A:here
  endif !! (edo) then

end subroutine readfchk
