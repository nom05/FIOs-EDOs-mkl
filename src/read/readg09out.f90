! ** read/readg09out.f90 >> Reading a Gaussian log file
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


subroutine readg09out(irewind,fieldinit)

  use omp_lib
  use SUFR_kinds ,only : double,dbl
  use SUFR_system,only : quit_program_error,warn
  use SUFR_text  ,only : int2str,lowercase,dbl2str
  use memory_use
  use commonmod  ,only : igau,fileopen,nato,nfb,nalpha,nbeta,nel,ndimtriang,dip,lop,nop,textxyz,d,Field &
                       , nmo,dmo,debuglabel,ilendbglbl,neqread,comptext,print_message,nclas             &
                       , ncmpeqv,nxyzcom &
                       , spec_comps           &       !! loptions( -2)
                       , lfieldbyhand=>lfield &       !! loptions( -1)
                       , compute_dmo          &       !! loptions(  1)
                       , dynamic_alpha        &       !! loptions(  2)
                       , read_mult_mat        &       !! loptions(  4)
                       , nocolor              &       !! loptions(  6)
                       , debug                &       !! loptions(  7)
                       , mem_info            !&       !! loptions( 14)
  use readmod    ,only : read_check_value,intread,p1fromgauss,intrdrmlch,realvec,pmofromgauss

  implicit none
  integer  ( kind=   4 )                          :: NAtoms,NBasis,nfbdim,nmodim,ndim
  integer  ( kind=   4 )                          :: i,ii,j,k,l,iii,jjj !! temp
  integer  ( kind=   4 ),allocatable,dimension(:) :: itmp,itmp2 !! temp
  integer  ( kind=   4 ),optional,intent(in)      :: irewind
  integer  ( kind=   4 )                          :: irewinddef=0,noop,noopmo,nroute,nbldm &
                                                   , itype(size(neqread))=1000,iclasp(maxval(neqread))=1000 &
                                                   , iclasn(maxval(neqread))=1000
  real     ( kind=   8 )                          :: fieldvec(3),fieldinic,temp
  real     ( kind=   8 ),optional,intent(in)      :: fieldinit
  real     ( kind=   8 ),allocatable,dimension(:) :: arrabf,arramo,arratmp,arratmp2
  character(  len=  80 )                          :: label,line
  character(  len=1000 )                          :: arg
  character(  len=   3 ),dimension(-4:4),parameter:: txtfldcmp=(/ 'f-z' &
                                                                , 'f-y' &
                                                                , 'f-x' &
                                                                , 'f0d' &
                                                                , 'edo' &
                                                                , 'f0s' &
                                                                , 'f+x' &
                                                                , 'f+y' &
                                                                , 'f+z' /)
  logical,target                    ,dimension(5) :: lgauopts = .FALSE.
  logical,pointer                                 :: lpolar    &  !! Option: polar      , lgauopts(1)
                                                   , ldyna     &  !! Option: cphf=rdfreq, lgauopts(2)
                                                   , lfield    &  !! Option: field      , lgauopts(3)
                                                   , liop10eq2 &  !! Option: iop33(10=2), lgauopts(4)
                                                   , liop3eq1     !! Option: iop33(3=1) , lgauopts(5)
  logical                                         :: lfirst = .TRUE.
   lpolar    => lgauopts(1)
   ldyna     => lgauopts(2)
   lfield    => lgauopts(3)
   liop10eq2 => lgauopts(4)
   liop3eq1  => lgauopts(5)
   nfbdim    =  ndimtriang(nfb,.FALSE.)
   nmodim    =  nmo*nmo

   if (present(irewind)) irewinddef = irewind
   if (present(fieldinit)) then
      if (fieldinit.NE.0._dbl.AND..NOT.lfieldbyhand) &
          call quit_program_error('A field was specified but the option "field" is unset!',1,nocolor)
      if (fieldinit.EQ.0._dbl.AND.     lfieldbyhand) &
          call quit_program_error('A null field was specified !',1,nocolor)
   endif !! (present(fieldinit) then
   if (allocated(d).OR.allocated(dmo)) call quit_program_error('DBF or DMO are unexpectedly allocated',1,nocolor)
   if (debug) PRINT *,debuglabel(:ilendbglbl),'Initial requested calcs = ',lop(:nop)
   if (debug) PRINT *,debuglabel(:ilendbglbl),'Initial equivalent types = ',neqread(lop(:nop))
!
! NOTE: Allocating here the largest arrays, memory fragmentation is avoided if the program knows the requested
!       components.
!
   if (spec_comps) then
      call myalloc(d,nfb,nfb,nop,'readg09out','dbf',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:FIOs
      if (.NOT.compute_dmo) call myalloc(dmo,nmo,nmo,nop,'readg09out','dmo',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:FIOs
   endif !! (spec_comps) then
! ---------------------------------------------------------
   if (irewinddef.EQ.0 ) rewind(igau)
   noop   = 1
   noopmo = 0
   nroute = 0
 reading: do
      read (igau,'(A)',iostat=iii) line
      if (iii.NE.0) then
         jjj = 1
         exit
      endif
!
!                  << Find 'Normal termination' or STOP >>
!
      if (index(line(1:35),'Normal termination of Gaussian').GT.0) noop = noop+1
!
!                  << Find 'Gaussian route line' >>
!
      if (index(line(1:4),'-').GT.0) then      !! Finding text blocks including '-'
         read (igau,'(A)',iostat=iii) line
         if (iii.NE.0) then
           jjj = 1
           call quit_program_error("A problem was found while G09 output reading (-)",1,nocolor)
         endif
         if (index(line(1:4),'#').GT.0) then   !! Finding text blocks including '#' (route section)
            nroute = nroute+1
            if (debug) PRINT *,debuglabel(:ilendbglbl),nroute,'=nroute noop=',noop
!                         vvvvvv---- Task is not finished 
            if (nroute.NE.noop) call quit_program_error('At least one task of Gaussian calculation &
                                                        &is unfinished correctly',1,nocolor)
            arg (1:) = trim(line(2:))
    mkarg:  do                                 !! Concatenating route section into one line
               read (igau,'(A)',iostat=iii) line
               if (iii.NE.0) then
                 jjj = 1
                 call quit_program_error("A problem was found while G09 output reading (#)",1,nocolor)
               endif
               if (index(line(1:4),'-').GT.0) then
                  exit
               else
                  arg (1+len(trim(arg)):) = trim(line(2:))
               endif !! (index(line(1:4),'-').GT.0) then
            enddo mkarg
            if (debug) PRINT *,debuglabel(:ilendbglbl),'LOOP = "mkarg" finished'
            lpolar    = index(lowercase(trim(arg)),'polar'      ).GT.0    !! Searching options
            ldyna     = index(lowercase(trim(arg)),'cphf=rdfreq').GT.0
            lfield    = index(lowercase(trim(arg)),'field='     ).GT.0
            liop10eq2 = index(lowercase(trim(arg)),'iop33(10=2)').GT.0 &
                    .OR.index(lowercase(trim(arg)),'10/33=2'    ).GT.0
            liop3eq1  = index(lowercase(trim(arg)),    '3=1').GT.0 &
                    .OR.index(lowercase(trim(arg)), '3/33=1').GT.0
     iopif: if (index(lowercase(trim(arg)),'iop33(10=2,3=1)'    ).GT.0 &
            .OR.index(lowercase(trim(arg)),'iop33(3=1,10=2)'    ).GT.0 &
            .OR.index(lowercase(trim(arg)),'iop(3/33=1,10/33=2)').GT.0 &
            .OR.index(lowercase(trim(arg)),'iop(10/33=2,3/33=1)').GT.0 ) then
               liop10eq2 = .TRUE.
               liop3eq1  = .TRUE.
            endif iopif
            if (debug) then
               PRINT *,debuglabel(:ilendbglbl),'COND = "iopif" finished'
               PRINT *,debuglabel(:ilendbglbl),'>>',trim(arg),'<<','(',noop,')'
               PRINT *,'        ','lpolar    = ',lpolar
               PRINT *,'        ','ldyna     = ',ldyna
               PRINT *,'        ','lfield    = ',lfield
               PRINT *,'        ','liop10eq2 = ',liop10eq2
               PRINT *,'        ','liop3eq1  = ',liop3eq1                      !! Multipole matrices
            endif !! (debug) then
!
!                  << Classifying types of Gaussian calculations >>
!
            if (lpolar.AND.liop10eq2.AND..NOT.ldyna.AND..NOT.lfield) itype(noop) =  1 !! Static  pol.
            if (lpolar.AND.liop10eq2.AND.     ldyna.AND..NOT.lfield) itype(noop) = -1 !! Dynamic pol.
            if (lpolar.AND.                   ldyna.AND.     lfield) &
               call quit_program_error('Dynamic polarizability with an electric field is incompatible',1,nocolor)
            if (lpolar.AND.liop10eq2.AND.                    lfield) then       !! Determine beta comp:
               if (lfieldbyhand) call warn('The field vector will be used to consider the component&
                                           & and sign. The strength will be replaced at the end of&
                                           & the reading process by the specified one by hand.')
               if (debug) PRINT *,debuglabel(:ilendbglbl),'Searching field in G09 log ...'
               call realvec(igau,5,'An electric field of',3,fieldvec)
               itype(noop) =  1
               j           =  0
               if (debug) PRINT *,debuglabel(:ilendbglbl),'Field vector:',fieldvec
               fieldinic = sqrt(dot_product(fieldvec,fieldvec))
               if (debug) PRINT *,debuglabel(:ilendbglbl),'Field strength:',fieldinic
               if (Field.NE.0._dbl) then
                  if (fieldinic.NE.Field) then
                     label(:) = ''
                     label    = 'Previous field ('//dbl2str(Field,5)//') != Current field found ('&
                                //dbl2str(fieldinic,5)//')'
                     call quit_program_error(trim(label),1,nocolor)
                  endif !! (fieldinic.NE.Field) then
                  Field = fieldinic
               else
                  if (fieldinic.EQ.0._dbl) then
                     call quit_program_error('Current field purposely specified is zero!! &
                                             &Remove option',1,nocolor)
                  else
                     Field = fieldinic
                  endif !! (fieldinic.EQ.0._dbl) then
               endif !! (Field.NE.0._dbl) then
               do i = 1,3
                  if (fieldvec(i).GT.0._dbl) then
                     itype(noop) = itype(noop)+i
                     j           = j+1         !! Check repetitions
                  endif !! (fieldvec(i).GT.0._dbl) then
               enddo !! i = 1,3
               if (itype(noop).EQ.1) then
                  itype(noop) = -1
                  j           =  0
                  do i = 1,3
                     if (fieldvec(i).LT.0._dbl) then
                        itype(noop) = itype(noop)-i
                        j           = j+1      !! Check repetitions
                     endif !! (fieldvec(i).LT.0._dbl) then
                  enddo !! i = 1,3
               endif !! (itype.EQ.1) then
               if (j.GT.1.OR.j.LE.0) call quit_program_error('Field component is repeated?',1,nocolor)
               if (abs(itype(noop)).GT.4 ) &
                       call quit_program_error('Wrong field was specified in the G09 input',1,nocolor)
            endif !! (lpolar.AND.liop10eq2.AND.                    lfield) then
            if (debug) PRINT *,debuglabel(:ilendbglbl),'Processing itype=',itype(noop)
!
!                  << Checking some variables previously obtained from fchk file >>
!
            call read_check_value(igau,1,nfb   ,NBasis,'basis functions,')
            call read_check_value(igau,1,nalpha,iii   ,'alpha electrons' )
            backspace(igau)
            call read_check_value(igau,4,nbeta ,iii   ,'beta electrons'  )
            if (nel.NE.nalpha+nbeta) call quit_program_error("# electrons value is not kept unchanged"&
                                                             ,1,nocolor)
            call read_check_value(igau,2,nato  ,NAtoms,'NAtoms='         )
            if (debug) PRINT *,debuglabel(:ilendbglbl),NBasis,'=NBasis nalpha=',nalpha,&
                                                       nbeta,'=nbeta NAtoms=',NAtoms
!
!                  << Read multipole matrices if requested >>
!
            if (debug) PRINT *,debuglabel(:ilendbglbl),'dip?->',read_mult_mat.AND.liop3eq1.AND.lfirst,'(',&
                                                               read_mult_mat,liop3eq1,lfirst,')'
            if (read_mult_mat.AND.liop3eq1.AND.lfirst) then
               if (     allocated(arratmp)) call mydealloc(arratmp,'readg09out','arratmp',nocolor=nocolor,verbose=mem_info) !+A:here
               call myalloc(arratmp,nfbdim,'readg09out','arratmp',nocolor=nocolor,verbose=mem_info) !+D:here
               call myalloc(dip,nfb,nfb,3,'readg09out','dip',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.)!+D:conductance_calc&FIOs
               do i = 1,3
                  iii  = intread(igau,6,'Multipole matrices I')
                  if (iii.NE.i) call quit_program_error('One of the multipole matrices was not &
                                                        &found ?!?!',1,nocolor)
                  call p1fromgauss(igau,nfb,arratmp(:nfbdim))

!$omp parallel shared ( nfb,dip,arratmp,i ) private ( j,k )
  !$omp do
                  do j = 1,nfb
                     dip(j,j,i)    = arratmp(ndimtriang(j,.FALSE.)  )
                     do k = 1,j-1
                        dip(j,k,i) = arratmp(ndimtriang(j,.TRUE. )+k)
                        dip(k,j,i) = dip(j,k,i)
                     enddo !! k = 1,nfb
                  enddo !! j = 1,nfb
  !$omp end do
!$omp end parallel

                  if (debug) write(*,'(1X,A,I1,1X,E13.6,A,E13.6,/,8X,E13.6,A,E13.6)') &
                                                              debuglabel(:ilendbglbl) &
                           ,i,dip(1  ,1,i),'=dip(1,1,i) dip(1,n,i)=',dip(1  ,nfb,i)   &
                             ,dip(nfb,1,i),'=dip(n,1,i) dip(n,n,i)=',dip(nfb,nfb,i)
               enddo !! i = 1,3
               if (allocated(arratmp)) call mydealloc(arratmp,'readg09out','arratmp',nocolor=nocolor,verbose=mem_info) !+A:here
            else if (read_mult_mat.AND.liop3eq1.AND.     allocated(dip)) then
               call warn('Mult. matrices repeated several times will be skipped',1,nocolor=nocolor)
            endif !! (read_mult_mat.AND.liop3eq1.AND.lfirst) then
!
!                  << Read DBF and DMO (the latter if requested) >>
!
!                       * Reading or skipping if not requested:
            if (any(neqread(lop(:nop)).EQ.abs(itype(noop)))) then
!                       * Expanding temporary BF/MO 'arrays' to read matrices:

               if (dynamic_alpha.AND.itype(noop).EQ.-1) then
                  ii = 3
               elseif (dynamic_alpha.AND.itype(noop).NE.-1) then
                       call quit_program_error('Incompatible Gaussian task detected. If static calculations &
                                               &are required, specify independently in another calc.',1,nocolor)
               else
                  ii = 0
               endif !! (dynamic_alpha.AND.itype(noop).EQ.-1) then

               if (allocated(arrabf)) then
                  if (3*nfbdim*(noop-1).NE.size(arrabf)) call quit_program_error('BF array size and &
                                                             &current Gaussian task do not correspond &
                                                             &consistently',1,nocolor)
                  call myalloc(arratmp,size(arrabf) + 3*nfbdim,'readg09out','arratmp',nocolor=nocolor,verbose=mem_info) !+D:here
                  arratmp(:size(arrabf)) = arrabf
                  call mydealloc(arrabf,'readg09out','arrabf',nocolor=nocolor,verbose=mem_info) !+A:here
               else
                  if (noop.GT.1) call quit_program_error('BF array size and current Gaussian task &
                                                         &do not correspond consistently',1,nocolor)
                  if (debug) call print_message('d',">>>> BF P1 First allocation")
                  call myalloc(arratmp,3*nfbdim,'readg09out','arratmp',nocolor=nocolor,verbose=mem_info) !+D:here
               endif !! (allocated(arrabf)) then

               if (debug) PRINT *,debuglabel(:ilendbglbl),'Expand DMO?->',.NOT.compute_dmo.AND.&
                                                           (itype(noop).EQ.-1.OR.itype(noop).EQ.1),&
                            '(',.NOT.compute_dmo,(itype(noop).EQ.-1.OR.itype(noop).EQ.1),')'
               if (allocated(arramo).AND..NOT.compute_dmo.AND.(itype(noop).EQ.-1.OR.itype(noop).EQ.1)) then
                  if (3*nmodim*(noopmo).NE.size(arramo)) call quit_program_error('MO array size and &
                                                             &current Gaussian task do not correspond &
                                                             &consistently',1,nocolor)
                  noopmo = noopmo+1
                  call myalloc(arratmp2,size(arramo) + 3*nmodim,'readg09out','arratmp2',nocolor=nocolor,verbose=mem_info) !+D:here
                  arratmp2(:size(arramo)) = arramo
                  call mydealloc(arramo,'readg09out','arramo',nocolor=nocolor,verbose=mem_info) !+A:here
               else if (.NOT.compute_dmo.AND.(itype(noop).EQ.-1.OR.itype(noop).EQ.1)) then
                  if (noopmo.GT.1) call quit_program_error('MO array size and current Gaussian task &
                                                         &do not correspond consistently',1,nocolor)
                  if (debug) call print_message('d',">>>> MO P1 First allocation")
                  noopmo = noopmo+1
                  call myalloc(arratmp2,3*nmodim,'readg09out','arratmp2',nocolor=nocolor,verbose=mem_info) !+D:here
               endif !! (allocated(arramo).AND..NOT.compute_dmo) then

               if (debug) PRINT *,debuglabel(:ilendbglbl),'Incremental BF Array size=',size(arratmp )
               if (debug.AND..NOT.compute_dmo) PRINT *,debuglabel(:ilendbglbl), &
                                                          'Incremental MO Array size=',size(arratmp2)
               do j = 1,3
                  label(:) = ''
                  label = 'P1 alpha (symm , AO basis) for IC =  '//int2str(j+ii)
                  iii = intrdrmlch(igau,10     ,1  ,trim(label))
                  if (iii.NE.j+ii) then
                     label(:) = ''
                     label = 'DBF derivative ('//int2str(i-1)//','//trim(textxyz(j))//') was not found'
                     call quit_program_error(trim(label),1,nocolor)
                  endif !! (iii.NE.j+ii) then
                  if (debug) then
                     PRINT *,debuglabel(:ilendbglbl),'P1 alpha symm with AO basis ',iii,' FOUND->',j+ii
                     PRINT *,debuglabel(:ilendbglbl),'Limits of BF array for new values:' &
                                                    , nfbdim*(3*noop+j-4)+1,nfbdim*(3*noop+j-3)
                  endif !! (debug) then
                  call p1fromgauss(igau,nfb,arratmp(nfbdim*(3*noop+j-4)+1:nfbdim*(3*noop+j-3)))
                  if (debug) then
                     PRINT *,debuglabel(:ilendbglbl),arratmp(nfbdim*(3*noop+j-4)+1) &
                            ,"=(1,1) (nbf,nbf)="    ,arratmp(nfbdim*(3*noop+j-3))
                     PRINT *,debuglabel(:ilendbglbl),arratmp(nfbdim*(3*noop+j-3)-nfb+1),"=(nbf,1)"
                  endif !! (debug) then
               enddo !! j = 1,3
               if (allocated(arrabf)) call quit_program_error('BF Array is uncorrectly allocated in this&
                                                             & program step',1,nocolor)
               j = size(arratmp)*bit_size(j)*size_double_over_size_real/8
               call ProcDEAlloc(j,nocolor,mem_info,'readg09out','arratmp',perform_sub=.TRUE.)
               call move_alloc(arratmp ,arrabf)
               call   ProcAlloc(j,0,nocolor,mem_info,.FALSE.,'readg09out','arrabf',perform_sum=.TRUE.)
               if (allocated(arratmp )) call mydealloc(arratmp,'readg09out','arratmp',nocolor=nocolor,verbose=mem_info) !+A:here
               if (debug) call print_message('d','"move_alloc" correctly performed from tmp to BF array')
               if (.NOT.compute_dmo.AND.(itype(noop).EQ.-1.OR.itype(noop).EQ.1)) then
                  do j = 1,3
                     label(:) = ''
                     label = 'P1 (MO basis) (alpha) for IMat=    '//int2str(j+ii)
                     iii = intrdrmlch(igau,7     ,1  ,trim(label))
                     if (iii.NE.j+ii) then
                        label(:) = ''
                        label = 'DMO derivative ('//int2str(i-1)//','//trim(textxyz(j))//') was not found'
                        call quit_program_error(trim(label),1,nocolor)
                     endif !! (iii.NE.j+ii) then
                     if (debug) then
                        PRINT *,debuglabel(:ilendbglbl),'P1 alpha symm with MO basis ',iii,' FOUND->',j+ii
                        PRINT *,debuglabel(:ilendbglbl),'Limits of MO array for new values:' &
                                                       , nmodim*(3*noopmo+j-4)+1,nmodim*(3*noopmo+j-3) &
                                                       , noopmo
                     endif !! (debug) then
                     call pmofromgauss(igau,nmo,arratmp2(nmodim*(3*noopmo+j-4)+1:nmodim*(3*noopmo+j-3)))
                     if (debug) then
                        PRINT *,debuglabel(:ilendbglbl),arratmp2(nmodim*(3*noopmo+j-4)+1) &
                               ,"=(1,1) (nmo,nmo)="    ,arratmp2(nmodim*(3*noopmo+j-3))
                        PRINT *,debuglabel(:ilendbglbl),arratmp2(nmodim*(3*noopmo+j-3)-nmo+1),"=(nmo,1)"
                     endif !! (debug) then
                  enddo !! j = 1,3
                  if (allocated(arramo)) call quit_program_error('MO array is uncorrectly allocated in this&
                                                                & program step',1,nocolor)
                  j = size(arratmp2)*bit_size(j)*size_double_over_size_real/8
                  call ProcDEAlloc(j,nocolor,mem_info,'readg09out','arratmp2',perform_sub=.TRUE.)
                  call move_alloc(arratmp2,arramo)
                  call   ProcAlloc(j,0,nocolor,mem_info,.FALSE.,'readg09out','arramo',perform_sum=.TRUE.)
                  if (debug) call print_message('d','"move_alloc" correctly performed from tmp to MO array')
                  if (allocated(arratmp2)) &
                          call mydealloc(arratmp2,'readg09out','arratmp2',nocolor=nocolor,verbose=mem_info) !+A:here
               elseif (.NOT.compute_dmo) then
                  if (debug) call print_message('d','DMO is computed for beta')
               endif !! (.NOT.compute_dmo) then
            else
               if (debug) PRINT *,debuglabel(:ilendbglbl),'Skipping unrequested itype=',itype(noop)
               noop   = noop  -1
               noopmo = noopmo-1
               nroute = nroute-1
            endif !! (any(neqread(lop(:nop)).EQ.abs(itype(noop)))) then
            if (debug) PRINT *,debuglabel(:ilendbglbl),'LOOP = DBF finished'
!
!                  << Cleaning >>
!
            lgauopts = .FALSE.
            if (lfirst) lfirst = .FALSE.
            if (debug) PRINT *,debuglabel(:ilendbglbl),'COND = ROUTE section finished'
         endif !! (index(line(1:4),'#').GT.0) then
         if (debug) PRINT *,debuglabel(:ilendbglbl),'COND = "---" finished'
      endif !! (index(line(1:4),'-').GT.0) then
   enddo reading
   noop = noop-1
   if (debug) PRINT *,debuglabel(:ilendbglbl),'LOOP = "reading" finished, itype=',itype(:noop)
!
!                  << Some checks. Reduce 'lop' if necessary according to Gaussian output >>
!
   if (any(itype(:noop).EQ.1000)) call quit_program_error('All tasks of Gaussian output are not &
                                                 &classified or incompatible with FIOs computation. &
                                                 &Check carefully',1,nocolor)
   if (any(abs(itype(:noop)).GT.4)) call quit_program_error('Error assigning itype',1,nocolor)
   do i = 1,noop
      if (abs(itype(i)).GT.1) then
         if (.NOT.any(itype(:noop).EQ.-itype(i))) then
            label(:) = ''
            label    = 'One of the fields ('// trim(txtfldcmp(-itype(i)))//') was not considered'
            call quit_program_error(trim(label),1,nocolor)
         endif !! (.NOT.any(itype(:noop).EQ.-itype(i))) then
      endif !! (abs(itype).GT.1) then
   enddo !! i = 1,noop
   if (spec_comps) then
      do i = 1,nop
         if (.NOT.any(abs(itype(:noop)).EQ.neqread(lop(i)))) then
            label(:) = ''
            label    = 'Requested component ('//comptext(lop(i))//') is not completely included'
            call quit_program_error(trim(label),1,nocolor)
         endif !! (.NOT.any(itype(:noop).EQ.neqread(lop(i)))) then
      enddo !! i = 1,nop
   else
      if (debug) PRINT *,debuglabel(:ilendbglbl),'Initial lop:',lop(:nop)
      if (allocated(itmp )) call mydealloc(itmp,'readg09out','itmp',nocolor=nocolor,verbose=mem_info) !+A:here
      if (allocated(itmp2)) call mydealloc(itmp2,'readg09out','itmp2',nocolor=nocolor,verbose=mem_info) !+A:here
      call myalloc(itmp,nop,'readg09out','itmp',nocolor=nocolor,verbose=mem_info) !+D:here
      call myalloc(itmp2,nop,'readg09out','itmp2',nocolor=nocolor,verbose=mem_info) !+D:here
      itmp = lop(:nop)
      iii  = nop
      do i = nop,1,-1
         if (.NOT.any(abs(itype(:noop)).EQ.neqread(lop(i)))) then
            itmp2 = itmp
            if (i-1.GT.0) itmp(:i-1) = itmp2(:i-1)
            iii   = iii-1
            itmp(i:iii) = itmp2(i+1:iii+1)
         endif !! (.NOT.any(itype(:noop).EQ.neqread(lop(i)))) then
      enddo !! i = 1,nop
      nop = iii
      lop(:nop) = itmp(:nop)
      if (allocated(itmp )) call mydealloc(itmp,'readg09out','itmp',nocolor=nocolor,verbose=mem_info) !+A:here
      if (allocated(itmp2)) call mydealloc(itmp2,'readg09out','itmp2',nocolor=nocolor,verbose=mem_info) !+A:here
      if (debug) PRINT *,debuglabel(:ilendbglbl),'Final lop:',lop(:nop)
   endif !! (spec_comps) then

   if (lfieldbyhand) Field = fieldinit
   if (debug) then
      if (jjj.EQ.1) PRINT *,debuglabel(:ilendbglbl),noop,'=noop, END OF FILE. Processing arrays...'
   endif !! (debug) then
   if (associated(lpolar)   ) nullify(lpolar   )
   if (associated(ldyna)    ) nullify(ldyna    )
   if (associated(lfield)   ) nullify(lfield   )
   if (associated(liop10eq2)) nullify(liop10eq2)
   if (associated(liop3eq1) ) nullify(liop3eq1 )
   
   nbldm = 0
   do i = 1,noop
      select case (itype(i))
         case (-1,1)
              nbldm           = nbldm+1
              iclasp(1)       = nbldm
              iclasn(1)       = nbldm
              nclas(nbldm)    = 1
         case (-2,2)
              if (itype(i).GT.0) then
                 nbldm        = nbldm+1
                 nclas(nbldm) = 2
                 iclasp(2)    = i
              else
                 iclasn(2)    = i
              endif !! (itype(i).GT.0) then
         case (-3,3)
              if (itype(i).GT.0) then
                 nbldm        = nbldm+1
                 nclas(nbldm) = 3
                 iclasp(3)    = i
              else
                 iclasn(3)    = i
              endif !! (itype(i).GT.0) then
         case (-4,4)
              if (itype(i).GT.0) then
                 nbldm        = nbldm+1
                 nclas(nbldm) = 4
                 iclasp(4)    = i
              else
                 iclasn(4)    = i
              endif !! (itype(i).GT.0) then
         case default
              call quit_program_error('Unclassified calculation type',1,nocolor)
      end select
   enddo !! i = 1,noop
   if (debug) then
      PRINT *,debuglabel(:ilendbglbl),'nbldm=' ,nbldm
      PRINT *,debuglabel(:ilendbglbl),'iclasp=',iclasp(:maxval(neqread))
      PRINT *,debuglabel(:ilendbglbl),'iclasn=',iclasn(:maxval(neqread))
      PRINT *,debuglabel(:ilendbglbl),'nclas =',nclas (:nbldm)
   endif !! (debug) then
   if (nbldm.GT.maxval(neqread)) call quit_program_error('Number of DMs is wrong',1,nocolor)
   ndim=0
   do i = 1,nbldm
      ndim = ndim +count(neqread==nclas(i))
   enddo !! i = 1,nbldm
   if (debug.AND..NOT.spec_comps) PRINT *,debuglabel(:ilendbglbl),'# arrays to be allocated in DBF & DMO=',ndim
   if (nop.NE.ndim) then
      label(:) = ''
      label    = 'Number of operations ('//int2str(nop)//') != dimension of operations ('//int2str(ndim)//')'
      call quit_program_error(trim(label),1,nocolor)
   endif !! (nop.EQ.ndim) then
   if (.NOT.allocated(d  )) &
           call myalloc(d  ,nfb,nfb,nop,'readg09out','dbf',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:FIOs
   if (.NOT.compute_dmo.AND..NOT.allocated(dmo)) &
           call myalloc(dmo,nmo,nmo,nop,'readg09out','dmo',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:FIOs

   do i = 1,nop
      if (debug) then
         PRINT *,debuglabel(:ilendbglbl),'>>>> Operation number=',i,neqread(lop(i)) &
                ,'(',trim(comptext(lop(i))),')=type  coord=',trim(textxyz(nxyzcom(lop(i))))
         PRINT *,debuglabel(:ilendbglbl),'     Initial position DBF (>0, <0):' &
                ,3*nfbdim*(iclasp(neqread(lop(i)))-1) + nfbdim*(nxyzcom(lop(i))-1) + 1 &
                ,3*nfbdim*(iclasn(neqread(lop(i)))-1) + nfbdim*(nxyzcom(lop(i))-1) + 1
         if (.NOT.compute_dmo.AND.neqread(lop(i)).EQ.1) PRINT *,debuglabel(:ilendbglbl) &
                ,'     Initial position DMO (>0, <0):' &
                ,3*nmodim*(iclasp(neqread(lop(i)))-1) + nmodim*(nxyzcom(lop(i))-1) + 1 &
                ,3*nmodim*(iclasn(neqread(lop(i)))-1) + nmodim*(nxyzcom(lop(i))-1) + 1
      endif !! (debug) then
      select case (neqread(lop(i)))
         case (1)
                    iii  = 3*nfbdim*(iclasp(neqread(lop(i)))-1) + nfbdim*(nxyzcom(lop(i))-1)

!$omp parallel default(none) shared( i,iii,nfb,d,arrabf ) private( k,l )
  !$omp do
                    do k = 1,nfb
                       d(k,k,i)    = arrabf( iii + ndimtriang(k,.FALSE.)   )
                       do l = 1,k-1
                          d(k,l,i) = arrabf( iii + ndimtriang(k,.TRUE. )+l )
                          d(l,k,i) = d(k,l,i)
                       enddo !! l = 1,nfb
                    enddo !! k = 1,nfb
  !$omp end do
!$omp end parallel

                 if (debug) call printdebug( &
                                        "arrabf"  , &
                                        "DBF", &
                                        textxyz(nxyzcom(lop(i))), &
                                        iii+1, &
                                        iii+ndimtriang(nfb,.FALSE.), &
                                        i, &
                                        d(1  ,1  ,i), &
                                        d(1  ,nfb,i), &
                                        d(nfb,1  ,i), &
                                        d(nfb,nfb,i)  &
                                           )

                 if (.NOT.compute_dmo) then
                    iii  = 3*nmodim*(iclasp(neqread(lop(i)))-1) + nmodim*(nxyzcom(lop(i))-1)

!$omp parallel default(none) shared( i,j,iii,nmo,dmo,arramo,nmodim,iclasp,jjj,nclas ) private( k,l )
  !$omp do
                    do k = 1,nmo
                       do l = 1,nmo
                          dmo(k,l,i) = 2._dbl*arramo( iii + nmo*(k-1) + l )
                       enddo !! l = 1,nmo
                    enddo !! k = 1,nmo
  !$omp end do
!$omp end parallel

                    if (debug) call printdebug( &
                                        "arramo"  , &
                                        "DMO (x2)", &
                                        textxyz(nxyzcom(lop(i))), &
                                        iii+1, &
                                        iii+nmo*nmo, &
                                        i, &
                                        dmo(1  ,1  ,i), &
                                        dmo(1  ,nmo,i), &
                                        dmo(nmo,1  ,i), &
                                        dmo(nmo,nmo,i)  &
                                             )

                 endif !! (.NOT.compute_dmo) then

         case default
                    temp = .5_dbl/Field
                    iii  = 3*nfbdim*(iclasp(neqread(lop(i)))-1) + nfbdim*(nxyzcom(lop(i))-1)
                    jjj  = 3*nfbdim*(iclasn(neqread(lop(i)))-1) + nfbdim*(nxyzcom(lop(i))-1)

!$omp parallel default(none) shared( i,iii,jjj,nfb,d,arrabf,temp ) private( k,l )
  !$omp do
                    do k = 1,nfb
                       d(k,k,i)    = temp*( arrabf( iii + ndimtriang(k,.FALSE.)     ) &
                                           -arrabf( jjj + ndimtriang(k,.FALSE.)     ) )
                       do l = 1,k-1
                          d(k,l,i) = temp*( arrabf( iii + ndimtriang(k,.TRUE. ) + l ) &
                                           -arrabf( jjj + ndimtriang(k,.TRUE. ) + l ) )
                          d(l,k,i)= d(k,l,i)
                       enddo !! l = 1,k-1
                    enddo !! k = 1,nfb
  !$omp end do
!$omp end parallel

                 if (debug) call printdebugdbl(           &
                                        "arrabf","DBF",textxyz(nxyzcom(lop(i))), &
                                        iii+1, &
                                        iii+ndimtriang(nfb,.FALSE.), &
                                        jjj+1, &
                                        jjj+ndimtriang(nfb,.FALSE.), &
                                        i, &
                                        d(1  ,1  ,i), &
                                        d(1  ,nfb,i), &
                                        d(nfb,1  ,i), &
                                        d(nfb,nfb,i)  &
                                              )

      end select !! case (neqread(lop(i)))
   enddo !! i = 1,nop

   if (allocated(arrabf)) call mydealloc(arrabf,'readg09out','arrabf',nocolor=nocolor,verbose=mem_info) !+A:here
   if (allocated(arramo)) call mydealloc(arramo,'readg09out','arramo',nocolor=nocolor,verbose=mem_info) !+A:here
     arg(:) = ''
   label(:) = ''
   label    = '(A,1X,'//int2str(noop)//'(A,1X))'
   write(arg,label) "Types of DMs detected in the G09 log file:",txtfldcmp(itype(:noop))
   call print_message('i',trim(arg))
     arg(:) = ''
   label(:) = ''
   label    = '(A,1X,'//int2str(nop)//'(A,1X))'
   write(arg,label) "Components to be calculated:",comptext(lop(:nop))
   call print_message('i',trim(arg))
   if (   read_mult_mat.AND..NOT.allocated(dip)) call quit_program_error('Multpol. matrices NOT READ!!',1,nocolor)
   if (.NOT.compute_dmo.AND..NOT.allocated(dmo)) call quit_program_error('DMO NOT READ!!',1,nocolor)
   if (     compute_dmo.AND.     allocated(dmo)) call quit_program_error('DMO is unexpectedly allocated!!',1,nocolor)
   if (                     .NOT.allocated(d  )) call quit_program_error('DBF NOT READ!!',1,nocolor)

 

!!!!!!!!!!!!!!!!!!
!!!! OLD CODE !!!!
!!!!!!!!!!!!!!!!!!
!!nfbdim  = ndimtriang(nfb,.FALSE.)
!!nmodim  = nmo*nmo
!!if (read_mult_mat) allocate(dip(nfb,nfb,3))
!!allocate(d(nfb,nfb,12),dmo(nmo,nmo,12),array(21*nfbdim))
!!if (.NOT.compute_dmo) allocate(arra2(3*nmodim))
!!call read_check_value(igau,2,nato  ,NAtoms,'NAtoms='         )
!!call read_check_value(igau,1,nfb   ,NBasis,'basis functions,')
!!call read_check_value(igau,1,nalpha,iii   ,'alpha electrons' )
!!backspace(igau)
!!call read_check_value(igau,4,nbeta ,iii   ,'beta electrons'  )
!!if (nel.NE.nalpha+nbeta) call quit_program_error("# electrons value is not kept unchanged",1,nocolor)
!
!     I will read the electric field strength and rewind the file. It is not
!     the best way but I will avoid to suppose this value is after the 
!     multipole matrices. They could be printed in any step of the Gaussian
!     calculation. Only for the static approach
!
!!if (.NOT.dynamic_alpha) then
!!   if (.NOT.lfieldbyhand) then
!!      if (debug) PRINT *,debuglabel(:ilendbglbl),'Searching field in G09 log ...'
!!      call realvec(igau,5,'An electric field of',3,fieldvec)
!!      if (debug) PRINT *,debuglabel(:ilendbglbl),'Field vector:',fieldvec
!!      Field = sqrt(dot_product(fieldvec,fieldvec))
!!      if (debug) PRINT *,debuglabel(:ilendbglbl),'Field strength:',Field
!!      rewind(igau)
!!   endif !! (.NOT.lfieldbyhand) then
!!endif !! (.NOT.dynamic_alpha) then

!!if (read_mult_mat) then
!!   do i = 1,3
!!      iii  = intread(igau,6,'Multipole matrices I')
!!      if (iii.NE.i) call quit_program_error('One of the multipole matrices was not found ?!?!',1,nocolor)
!!      call p1fromgauss(igau,nfb,array(1:nfbdim))

!!omp parallel shared ( nfb,dip,array,i ) private ( j,k )

!!!$omp do
!!      do j = 1,nfb
!!         dip(j,j,i)    = array(ndimtriang(j,.FALSE.)  )
!!         do k = 1,j-1
!!            dip(j,k,i) = array(ndimtriang(j,.TRUE. )+k)
!!            dip(k,j,i) = dip(j,k,i)
!!         enddo !! k = 1,nfb
!!      enddo !! j = 1,nfb
!!!$omp end do

!!omp end parallel

!!      if (debug) write(*,'(1X,A,I1,1X,E13.6,A,E13.6,/,8X,E13.6,A,E13.6)') debuglabel(:ilendbglbl)&
!!                      ,i,dip(1  ,1,i),'=dip(1,1,i) dip(1,n,i)=',dip(1  ,nfb,i) &
!!                        ,dip(nfb,1,i),'=dip(n,1,i) dip(n,n,i)=',dip(nfb,nfb,i)
!!   enddo !! i = 1,3
!!endif !! (read_mult_mat) then
!
!!    Loop that reads perturbed density matrices on the basis of AOs
!!    for w/o and w/ electric field applied in different directions
!
!!    Structure of array:
!!--- A L P H A ---------------------------------
!!         *  1 -> x( 0) (no field) i=1
!!         *  2 -> y( 0)
!          *  3 -> z( 0)
! --- B E T A   ----------------------------------
!          *  4 -> x(+x) (+x field) i=2
!          *  5 -> y(+x)
!          *  6 -> z(+x)
!          *  7 -> x(+y) (+y field) i=3
!          *  8 -> y(+y)
!          *  9 -> z(+y)
!          * 10 -> x(+z) (+z field) i=4
!          * 11 -> y(+z)
!          * 12 -> z(+z)
!          * 13 -> x(-x) (-x field) i=5
!          * 14 -> y(-x)
!          * 15 -> z(-x)
!          * 16 -> x(-y) (-y field) i=6
!          * 17 -> y(-y)
!          * 18 -> z(-y)
!          * 19 -> x(-z) (-z field) i=7
!          * 20 -> y(-z)
!          * 21 -> z(-z)
!!i = 1
!!   ii = 0
!!   if (dynamic_alpha) ii = 3
!!   if (debug) PRINT *,debuglabel(:ilendbglbl),'ii=',ii
!!   do j = 1,3
!!      label(:) = ''
!!      label = 'P1 alpha (symm , AO basis) for IC =  '//int2str(j+ii)
!!      iii = intrdrmlch(igau,10     ,1  ,trim(label))
!!      if (iii.NE.j+ii) then
!!         label(:) = ''
!!         label = 'DBF derivative ('//int2str(i-1)//','//trim(textxyz(j))//') was not found'
!!         call quit_program_error(trim(label),1,nocolor)
!!      endif !! (iii.NE.j+ii) then
!!      if (debug) then
!!         PRINT *,debuglabel(:ilendbglbl),'P1 alpha symm with AO basis ',iii,' FOUND'
!!         PRINT *,debuglabel(:ilendbglbl),'Limits of reduced array(DBF):',nfbdim*(3*i+j-4)+1,nfbdim*(3*i+j-3)
!!      endif !! (debug) then
!!      call p1fromgauss(igau,nfb,array(nfbdim*(3*i+j-4)+1:nfbdim*(3*i+j-3)))
!!! n=nfbdim         3n(i-1) + n(j-1) <-|||||||||||||||    |||||||||||||||
!!!                                 3n(i-1) + n(j-1) + n <-|||||||||||||||
!!      if (debug) then
!!         PRINT *,debuglabel(:ilendbglbl),array(nfbdim*(3*i+j-4)+1),"=(1,1) (nbf,nbf)=",array(nfbdim*(3*i+j-3))
!!         PRINT *,debuglabel(:ilendbglbl),array(nfbdim*(3*i+j-3)-nfb+1),"=(nbf,1)"
!!      endif !! (debug) then
!!   enddo !! j = 1+i-1,3+i-1

!!   if (.NOT.compute_dmo) then    !! Read DMO. ii comes from previous steps (static or dynamic)
!!      if (debug) PRINT *,debuglabel(:ilendbglbl),'Reading of DMO requested'
!!      do j = 1,3
!!         label(:) = ''
!!         label = 'P1 (MO basis) (alpha) for IMat=    '//int2str(j+ii)  !! Check this when UHF was implemented
!!         iii = intrdrmlch(igau, 7     ,1  ,trim(label))
!!         if (iii.NE.j+ii) then
!!            label(:) = ''
!!            label = 'DMO derivative ('//int2str(i-1)//','//trim(textxyz(j))//') was not found'
!!            call quit_program_error(trim(label),1,nocolor)
!!         endif !! (iii.NE.j+ii) then
!!         if (debug) then
!!            PRINT *,debuglabel(:ilendbglbl),'P1 alpha with MO basis ',iii,' FOUND'
!!            PRINT *,debuglabel(:ilendbglbl),'Limits of reduced array(DMO):',nmodim*(3*i+j-4)+1,nmodim*(3*i+j-3)
!!         endif !! (debug) then
!!         call pmofromgauss(igau,nmo,arra2(nmodim*(3*i+j-4)+1:nmodim*(3*i+j-3)))
!!!    n=nmodim          3n(i-1) + n(j-1) <-||||||||||||||||   ||||||||||||||||
!!!                                     3n(i-1) + n(j-1) + n <-||||||||||||||||
!!         if (debug) then
!!            PRINT *,debuglabel(:ilendbglbl),arra2(nmodim*(3*i+j-4)+1  ),"=(1,1) (nmo,nmo)=",arra2(nmodim*(3*i+j-3)      )
!!            PRINT *,debuglabel(:ilendbglbl),arra2(nmodim*(3*i+j-4)+nmo),"=(1,nmo) (nmo,1)=",arra2(nmodim*(3*i+j-3)-nmo+1)
!!         endif !! (debug) then
!!      enddo !! j = 1,3
!!   endif !! (.NOT.compute_dmo) then

!!if (.NOT.dynamic_alpha) then  !! Static case (beta)
!!   do i = 2,7
!!      do j = 1,3
!!         label(:) = ''
!!         label = 'P1 alpha (symm , AO basis) for IC =  '//int2str(j)
!!         iii = intrdrmlch(igau,10     ,1  ,trim(label))
!!         if (iii.NE.j) then
!!            label(:) = ''
!!            label = 'DBF derivative ('//int2str(i-1)//','//trim(textxyz(j))//') was not found'
!!            call quit_program_error(trim(label),1,nocolor)
!!         endif !! (iii.NE.j) then
!!         call p1fromgauss(igau,nfb,array(nfbdim*(3*i+j-4)+1:nfbdim*(3*i+j-3)))
!!  !! n=nfbdim         3n(i-1) + n(j-1) <-|||||||||||||||    |||||||||||||||
!!  !!                                 3n(i-1) + n(j-1) + n <-|||||||||||||||
!!      enddo !! j = 1,3
!!   enddo !! i = 2,7

!
!!    Computing derivatives for beta
!

!!   temp = Field*2._dbl
!!   do i = 1,3
!!      do j = 1,3

!!omp parallel default(none) shared( i,j,nfbdim,array,temp ) private( k )
!!!$omp do
!!         do k = 1,nfbdim
!!           array(nfbdim*(3*i+j-1)+k) = (array(nfbdim*(3*i+j-1)+k) - array(nfbdim*(3*i+j+8)+k)) / temp
!!! n = nfbdim           |||||||-> 3n(no field) + 3n(i-1) + n(j-1)        ||||||||||||||||
!!!                        3n(no field) + 9n(+field) + 3n(i-1) + n(j-1) <-||||||||||||||||
!!          enddo !! k = 1,nfbdim
!!!$omp end do
!!omp end parallel

!!      enddo !! j = 1,3
!!   enddo !! i = 1,3
!!endif !! (.NOT.dynamic_alpha) then
!! 
!!i = 0
!!      do j = 1,3
!!omp parallel default(none) shared( i,j,nfb,d,array,nfbdim ) private( k,l )
!!!$omp do
!!         do k = 1,nfb
!!            d(k,k,3*i+j)      = array(ndimtriang(k,.FALSE.)    + nfbdim*(j-1) + 3*nfbdim*i)
!!            do l = 1,k-1
!!               d(k,l,3*i+j)   = array(ndimtriang(k,.TRUE. )+ l + nfbdim*(j-1) + 3*nfbdim*i)
!!               d(l,k,3*i+j)   = d(k,l,3*i+j)
!!            enddo !! l = 1,nfb
!!         enddo !! k = 1,nfb
!!!$omp end do
!!omp end parallel
!!         if (debug) write(*,'(1X,A,"Expanding DBF ...",/,9X,I2,A,1X,E13.6,A,E13.6,/,13X,E13.6,A,E13.6)') &
!!                debuglabel(:ilendbglbl),i,textxyz(j),d(1  ,1,3*i+j),'=(1,1,i) (1,n,i)=',d(1  ,nfb,3*i+j) &
!!                                       ,d(nfb,1,3*i+j),'=(n,1,i) (n,n,i)=',d(nfb,nfb,3*i+j)
!!      enddo !! j = 1,3

!!if (.NOT.dynamic_alpha) then
!!   do i = 1,3
!!      do j = 1,3
!!omp parallel default(none) shared( i,j,nfb,d,array,nfbdim ) private( k,l )
!!!$omp do
!!         do k = 1,nfb
!!            d(k,k,3*i+j)      = array(ndimtriang(k,.FALSE.)    + nfbdim*(j-1) + 3*nfbdim*i)
!!            do l = 1,k-1
!!               d(k,l,3*i+j)   = array(ndimtriang(k,.TRUE. )+ l + nfbdim*(j-1) + 3*nfbdim*i)
!!               d(l,k,3*i+j)   = d(k,l,3*i+j)
!!            enddo !! l = 1,nfb
!!         enddo !! k = 1,nfb
!!!$omp end do
!!omp end parallel
!!         if (debug) write(*,'(1X,A,"Expanding DBF ...",/,9X,I2,A,1X,E13.6,A,E13.6,/,13X,E13.6,A,E13.6)') &
!!             debuglabel(:ilendbglbl),i,textxyz(j),d(1  ,1,3*i+j),'=(1,1,i) (1,n,i)=',d(1  ,nfb,3*i+j) &
!!                             ,d(nfb,1,3*i+j),'=(n,1,i) (n,n,i)=',d(nfb,nfb,3*i+j)
!!      enddo !! j = 1,3
!!   enddo !! i = 1,3
!!endif !! (.NOT.dynamic_alpha) then

!
!!Creating dmo if requested. DMO is read from G09 output and multiplied by 2.
!
!!if (.NOT.compute_dmo) then
!!   do i = 1,3
!!omp parallel default(none) shared( i,j,nmo,dmo,arra2,nmodim ) private( k,l )
!!!$omp do
!!      do k = 1,nmo
!!         do l = 1,nmo
!!            dmo(k,l,i) = 2._dbl*arra2(nmodim*(i-1) + nmo*(k-1) + l)
!!         enddo !! l = 1,nmo
!!      enddo !! k = 1,nmo
!!!$omp end do
!!omp end parallel
!!      if (debug) write(*,'(1X,A,"Expanding DMO (after multiplying by 2)..."&
!!                        &,/,9X,I2,A,1X,E13.6,A,E13.6,/,13X,E13.6,A,E13.6)') &
!!               debuglabel(:ilendbglbl),0,textxyz(i),dmo(1  ,1,i),'=(1,1,i) (1,n,i)=',dmo(1  ,nmo,i) &
!!                     ,dmo(nmo,1,i),'=(n,1,i) (n,n,i)=',dmo(nmo,nmo,i)
!!   enddo !! i = 1,3
!!endif !! (.NOT.compute_dmo) then

!!deallocate(array)
!!if (allocated(arra2)) deallocate(arra2)
!!!!!!!!!!!!!!!!!!!!!!!!
!!!! END - OLD CODE !!!!
!!!!!!!!!!!!!!!!!!!!!!!!

   contains

      subroutine     printdebug(text1,text2,text3,inicp,endp,ini,a,b,c,d)
         implicit none
         character( len=*) :: text1,text2,text3
         integer  (kind=4) :: inicp,endp,ini
         real     (kind=8) :: a,b,c,d
         write(*,'(1X,A,"Taken ",A," from ...",I12,1X,"to",1X,I12)')     &
                      debuglabel(:ilendbglbl),trim(text1),inicp,endp
         write(*,'(1X,A,"Expanded ",A,1X,"...",/,9X,I2,A,1X,E13.6,A,E13.6,/,13X,E13.6,A,E13.6)') &
                      debuglabel(:ilendbglbl),trim(text2),ini,trim(text3), &
                      a,'=(1,1,i) (1,n,i)=',b,c,'=(n,1,i) (n,n,i)=',d
      end subroutine printdebug

      subroutine     printdebugdbl(text1,text2,text3,inicp,endp,inicn,endn,ini,a,b,c,d)
         implicit none
         character( len=*) :: text1,text2,text3
         integer  (kind=4) :: inicp,endp,inicn,endn,ini
         real     (kind=8) :: a,b,c,d
         write(*,'(1X,A,"Taken ",A,"+ from ...",I12,1X,"to",1X,I12)')     &
                      debuglabel(:ilendbglbl),trim(text1),inicp,endp
         write(*,'(1X,A,"Taken ",A,"- from ...",I12,1X,"to",1X,I12)')     &
                      debuglabel(:ilendbglbl),trim(text1),inicn,endn
         write(*,'(1X,A,"Expanded ",A,1X,"...",/,9X,I2,A,1X,E13.6,A,E13.6,/,13X,E13.6,A,E13.6)') &
                      debuglabel(:ilendbglbl),trim(text2),ini,trim(text3), &
                      a,'=(1,1,i) (1,n,i)=',b,c,'=(n,1,i) (n,n,i)=',d
      end subroutine printdebugdbl

end subroutine readg09out
