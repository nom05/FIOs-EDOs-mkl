! ** modules/conductance.f90 >> Module for conductance calculations. See references in manual.
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


module     conductance

  use SUFR_kinds                                  ,only: double,dbl
  use SUFR_system                                 ,only: quit_program_error
  use memory_use

  implicit none

  real   ( kind = double )                            :: dd = 0._dbl,GEDO
  real   ( kind = double ),parameter                  :: a2bohr = 1.8897259885789_dbl &
                                                       ,      e = 1.60217662E-19_dbl  &
                                                       ,      h = 1.0545718E-34_dbl
  integer( kind =      4 ),dimension(:  ),allocatable :: e1,e2
  integer( kind =      4 ),dimension(:,:),allocatable :: electrodes
  contains

    subroutine process_cond_par(label,lend,input2)

       use commonmod               ,only: nocolor,procnumb,hmnumb,mem_info &
                                        , debuglabel,ilendbglbl,debug,print_cond_mat
       use SUFR_system             ,only: quit_program_error
       use SUFR_text               ,only: int2str,lowercase

       implicit none

       character( len = * ),intent(in) :: label
       character( len = * )            :: input2
       integer                         :: lstatus,ii,jj
       logical             ,intent(in) :: lend

       if (index(lowercase(label(1:2)),'f=').GT.0) then
          if (debug) PRINT *,debuglabel(:ilendbglbl),'           * File: +'      ,trim(label(3:)),'+'
          input2(:) = ''
          read (label(3:),*,iostat=lstatus) input2
          if (debug) PRINT *,debuglabel(:ilendbglbl),'               input(2fch): +',trim(input2),'+'
          if (lstatus.NE.0) call quit_program_error("Problem reading 2nd fchk file name",1,nocolor)
       elseif (index(lowercase(label(1:2)),'d=').GT.0) then
          if (debug) PRINT *,debuglabel(:ilendbglbl),'           * Distance: +'  ,trim(label(3:)),'+'
          read (label(3:),*,iostat=lstatus) dd
          if (lstatus.NE.0) call quit_program_error("Problem reading electrode distance",1,nocolor)
          if (debug) PRINT *,debuglabel(:ilendbglbl),'               d=',dd
       elseif (index(lowercase(label(1:3)),'e1=').GT.0) then
          if (debug) PRINT *,debuglabel(:ilendbglbl),'           * Electrode1: +',trim(label(4:)),'+'
          call hmnumb(trim(label(4:)),ii,jj)
          if (debug) PRINT *,debuglabel(:ilendbglbl),'# Group of values (e1): ',trim(int2str(jj)),&
                       '; # atoms (e1): ',trim(int2str(ii))
          if (allocated(e1)) &
             call quit_program_error("Electrode-1 atoms array unexpectedly allocated",1,nocolor)
          call myalloc(e1,ii,'process_cond_par','e1',nocolor=nocolor,verbose=mem_info) !+D:here
          call procnumb(trim(label(4:)),len_trim(label(4:))+1,ii,jj,e1)
          if (debug) PRINT *,debuglabel(:ilendbglbl),'               e1: ',e1
       elseif (index(lowercase(label(1:3)),'e2=').GT.0) then
          if (debug) PRINT *,debuglabel(:ilendbglbl),'           * Electrode2: +',trim(label(4:)),'+'
          call hmnumb(trim(label(4:)),ii,jj)
          if (debug) PRINT *,debuglabel(:ilendbglbl),'# Group of values (e2): ',trim(int2str(jj)),&
                       '; # atoms (e2): ',trim(int2str(ii))
          if (allocated(e2)) &
             call quit_program_error("Electrode-2 atoms array unexpectedly allocated",1,nocolor)
             call myalloc(e2,ii,'process_cond_par','e2',nocolor=nocolor,verbose=mem_info) !+D:here
          call procnumb(trim(label(4:)),len_trim(label(4:))+1,ii,jj,e2)
          if (debug) PRINT *,debuglabel(:ilendbglbl),'               e2: ',e2
       elseif (index(lowercase(label(1:5)),'print').GT.0) then
          if (debug) PRINT *,debuglabel(:ilendbglbl),'           * Print conductance matrices'
          print_cond_mat = .TRUE.
       endif !! (index(label(1:2),'f=').GT.0) then

       if (lend) then
          if (dd.EQ.0._dbl) call quit_program_error("Null or unfilled distance",1,nocolor)
          if (.NOT.allocated(e1)) call quit_program_error("Atoms for Electrode-1 are not included",1,nocolor)
          if (.NOT.allocated(e2)) call quit_program_error("Atoms for Electrode-2 are not included",1,nocolor)
          if (allocated(electrodes)) then
             call quit_program_error("Array for electrodes unexpectedly allocated",1,nocolor)
          else
             call myalloc(electrodes,max(size(e1),size(e2))+1,2,'process_cond_par','electrodes' & !+D:conductance_calc
                                                             ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.)
             electrodes(                     :size(e1),1) = e1(:)
             electrodes(                     :size(e2),2) = e2(:)
             electrodes(size(electrodes(:,1))         ,1) = size(e1)
             electrodes(size(electrodes(:,2))         ,2) = size(e2)
             call mydealloc(e1,'process_cond_par','e1',nocolor=nocolor,verbose=mem_info) !+A:here
             call mydealloc(e2,'process_cond_par','e2',nocolor=nocolor,verbose=mem_info) !+A:here
          endif !! (allocated(electrode)) then
          dd = dd*a2bohr
       endif !! (lend) then

    end subroutine process_cond_par

    subroutine     conductance_calc

      use omp_lib
      use SUFR_kinds                                  ,only: double,dbl
      use SUFR_system                                 ,only: quit_program_error
      use SUFR_text                                   ,only: int2str
      use commonmod                                   ,only: nfb,edocomp,dip,over,ato2basis,ordD2Acedo,D2,iout  &
                                                           , print_message,nmotrunc,partial_time,ntime,cv=>cedo &
                                                           , nocolor                      & !! loptions(  6)
                                                           , debug,debuglabel,ilendbglbl  & !! loptions(  7)
                                                           , print_large_matrices         & !! loptions( 12)
                                                           , mem_info                     & !! loptions( 14)
                                                           , print_cond_mat              !& !! loptions( 15)

      use writemod                                    ,only: print_sym_mat,print_asym_mat

      implicit none

      integer  (kind =     4 )                            :: nedotrunc
      integer  (kind =     4 )                            :: i,j,k,l,ilec !! tmp
      real     (kind = double)                            :: CTE,GTotal
      real     (kind = double),allocatable,dimension(:  ) :: sEDO,rEDO,SSEDO,RREDO,Dr,S,R
      real     (kind = double),allocatable,dimension(:,:) :: rr,G,GS
      character( len =   200 )                            :: label

      call myalloc(sEDO,nmotrunc    ,'conductance_calc','sEDO',nocolor=nocolor,verbose=mem_info) !+D:here
      call myalloc(rEDO,nmotrunc    ,'conductance_calc','rEDO',nocolor=nocolor,verbose=mem_info) !+D:here
      call myalloc(  rr,nfb     ,nfb,'conductance_calc','rr'  ,nocolor=nocolor,verbose=mem_info) !+D:here

      if    (edocomp/abs(edocomp).GT.0) then
!$omp parallel default ( none ) shared ( rr,dip,edocomp )
  !$omp workshare
         rr(:,:) = -dip(:,:,abs(edocomp))
  !$omp end workshare
!$omp end parallel
         if (debug) call print_message('d','rr = -dip')
      elseif(edocomp/abs(edocomp).LT.0) then
!$omp parallel default ( none ) shared ( rr,dip,edocomp )
  !$omp workshare
         rr(:,:) =  dip(:,:,abs(edocomp))
  !$omp end workshare
!$omp end parallel
         if (debug) call print_message('d','rr =  dip')
      else
         call quit_program_error("Wrong field component detected",1,nocolor)
      endif !! (edocomp/abs(edocomp).GT.0) then

      call mydealloc(dip,'conductance_calc','dip',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:proc_dip&readg09out

!
!         Building and writing the Gmatrix
! 
      if (debug) then
         PRINT *,debuglabel(:ilendbglbl),'nmotrunc,nfb,electrodes:',nmotrunc,nfb &
                ,electrodes(:electrodes(size(electrodes(:,1)),1),1)
         PRINT *,'over(:30)=',over(:30,1)
         PRINT *,  'rr(:30)=',  rr(:30,1)
         PRINT *,debuglabel(:ilendbglbl),'Basis start for electrode=' &
                ,ato2basis(electrodes(:electrodes(size(electrodes(:,1)),1),1))%ifbstart
         PRINT *,debuglabel(:ilendbglbl),'Basis end for electrode='   &
                ,ato2basis(electrodes(:electrodes(size(electrodes(:,1)),1),1))%ifbend
         PRINT *,debuglabel(:ilendbglbl),'ordD2Acedo(30:)=',ordD2Acedo(nmotrunc:nmotrunc-29:-1)
         PRINT *,debuglabel(:ilendbglbl),'ordD2Acedo(:30)=',ordD2Acedo(:30)
         PRINT *,debuglabel(:ilendbglbl),'unsorted cedo(:30,1)=',cv(:30,1)
         PRINT *,debuglabel(:ilendbglbl),'unsorted cedo(:30,2)=',cv(:30,2)
         PRINT *,debuglabel(:ilendbglbl),'cedo(1,30:)=',cv(1,ordD2Acedo(nmotrunc:nmotrunc-29:-1))
      endif !! (debug) then

elec: do ilec = 1,2
         label    = '(2/,3X,42("*"),/,3X, &
                  &"** ANALYSING CONDUCTANCE OF ELECTRODE ",A," **",/,&
                  &3X,42("*"),/)'
         write(iout,label) int2str(ilec)
!$omp parallel default( none ) &
!$omp&         shared ( sEDO,rEDO,nmotrunc,nfb,electrodes,ilec,ato2basis,cv,over,rr ) &
!$omp&        private ( i,j,k,l )
  !$omp do
         do i = 1,nmotrunc
            sEDO(i) = 0._dbl
            do j = 1,nfb
               do k = 1,electrodes(size(electrodes(:,ilec)),ilec)
                  do l = ato2basis(electrodes(k,ilec))%ifbstart,ato2basis(electrodes(k,ilec))%ifbend
                     sEDO(i) = sEDO(i) + cv(j,i)*cv(l,i)*over(j,l)
                  enddo !! l = ato2basis(k)%ifbstart,ato2basis(k)%ifbend
               enddo !! k = 1,electrodes(size(electrodes(:,1)),1)
            enddo !! j = 1,nfb
            rEDO(i) = 0._dbl
            do j = 1,nfb
               do l = 1,nfb
                  rEDO(i) = rEDO(i) + cv(j,i)*cv(l,i)*rr(j,l)
               enddo !! l = 1,nfb
            enddo !! j = 1,nfb
         enddo !! i=1,nmotrunc
  !$omp end do
!$omp end parallel

         if (debug) then
            PRINT *,debuglabel(:ilendbglbl),'sEDO(:30)=',sEDO(:30)
            PRINT *,debuglabel(:ilendbglbl),'rEDO(:30)=',rEDO(:30)
         endif !! (debug) then
    
!$omp parallel default ( none ) shared ( nedotrunc,D2 )
  !$omp workshare
         nedotrunc = count(D2>=.00001_dbl)
  !$omp end workshare
!$omp end parallel

         if (debug) then
            PRINT *,debuglabel(:ilendbglbl),'EDO trunc=',nedotrunc
            PRINT *,debuglabel(:ilendbglbl),'15 last + 1 first discarded EDOs pairs (<1E-5):' &
                                           ,D2(ordD2Acedo(2*nedotrunc-29:2*nedotrunc+2))
         endif !! (debug) then

         call myalloc(SSEDO,nedotrunc          ,'conductance_calc','SSEDO',nocolor=nocolor,verbose=mem_info) !+D:here
         call myalloc(RREDO,nedotrunc          ,'conductance_calc','RREDO',nocolor=nocolor,verbose=mem_info) !+D:here
         call myalloc(   Dr,nedotrunc          ,'conductance_calc','Dr'   ,nocolor=nocolor,verbose=mem_info) !+D:here
         call myalloc(    S,nedotrunc          ,'conductance_calc','S'    ,nocolor=nocolor,verbose=mem_info) !+D:here
         call myalloc(    R,nedotrunc          ,'conductance_calc','R'    ,nocolor=nocolor,verbose=mem_info) !+D:here
         call myalloc(    G,nedotrunc,nedotrunc,'conductance_calc','G'    ,nocolor=nocolor,verbose=mem_info) !+D:here
         call myalloc(   GS,nedotrunc,nedotrunc,'conductance_calc','GS'   ,nocolor=nocolor,verbose=mem_info) !+D:here
    
!$omp parallel default( none ) &
!$omp&         shared ( nedotrunc,D2,ordD2Acedo,SSEDO,RREDO,Dr,sEDO,rEDO,S,R,dd ) &
!$omp&        private ( i )
  !$omp do
         do i=1,nedotrunc
            if (D2(ordD2Acedo(2*(i-1)+1)).GT.0) then
               SSEDO(i) = D2(ordD2Acedo(2*(i-1)+1)) * (sEDO(2*(i  )  ) - sEDO(2*(i-1)+1))
               RREDO(i) = D2(ordD2Acedo(2*(i-1)+1)) * (rEDO(2*(i  )  ) - rEDO(2*(i-1)+1))
                  Dr(i) = D2(ordD2Acedo(2*(i-1)+1))
            else
               SSEDO(i) = D2(ordD2Acedo(2*(i  )  )) * (sEDO(2*(i-1)+1) - sEDO(2*(i  )  ))
               RREDO(i) = D2(ordD2Acedo(2*(i  )  )) * (rEDO(2*(i-1)+1) - rEDO(2*(i  )  ))
                  Dr(i) = D2(ordD2Acedo(2*(i  )  ))
            endif
            S(i)        = SSEDO(i)/Dr(i)
            R(i)        = RREDO(i)/Dr(i)/dd
         enddo !! i=1,nedotrunc
  !$omp end do
!$omp end parallel
         if (debug) then
            PRINT *,debuglabel(:ilendbglbl),'SSEDO(:30)=',SSEDO(:30)
            PRINT *,debuglabel(:ilendbglbl),'RREDO(:30)=',RREDO(:30)
            PRINT *,debuglabel(:ilendbglbl),   'Dr(:30)=',   Dr(:30)
            PRINT *,debuglabel(:ilendbglbl),    'S(:30)=',    S(:30)
            PRINT *,debuglabel(:ilendbglbl),    'R(:30)=',    R(:30)
         endif !! (debug) then
    
         CTE    = (e*e)/(h*dd)
!$omp parallel default( none ) &
!$omp&         shared ( nedotrunc,SSEDO,RREDO,GTotal,G,ilendbglbl,GS,CTE ) &
!$omp&        private ( i )
  !$omp do
         do i=1,nedotrunc
            do j=1,nedotrunc
               G(i,j) = CTE*(1000000._dbl)*SSEDO(i)*RREDO(j)
            enddo !! j=1,nedotrunc
         enddo !! i=1,nedotrunc
  !$omp end do
  !$omp workshare
         GTotal = sum(G)
  !$omp end workshare
!$omp end parallel
         if (print_large_matrices.OR.print_cond_mat) call print_asym_mat(iout,nedotrunc,nedotrunc,G  &
                                                      ,'Conductance matrix on the basis of EDOs',numcolopt=5)
!$omp parallel default( none ) &
!$omp&         shared ( G,GS )
  !$omp workshare
         GS = (G+transpose(G))*.5_dbl
  !$omp end workshare
!$omp end parallel
         if (print_large_matrices.OR.print_cond_mat) call print_sym_mat(iout,nedotrunc,GS &
                                    ,'Symmetric part of conductance matrix on the basis of EDOs',numcolopt=5)
         write(iout,'(/,3X,">> EDOs CONDUCTANCE ",/)')
         do i=1,nedotrunc
            j     = max(ceiling(log10(dble(abs(nedotrunc)))),1)
            label = '(5X,&
                         &  " G(",I'//int2str(j)//',") = ",E13.6 &
                         &, " n(",I'//int2str(j)//',") = ",E13.6 &
                         &," Dq(",I'//int2str(j)//',") = ",E13.6 &
                         &," Dr(",I'//int2str(j)//',") = ",E13.6 &
                         & )'
            write(iout,label) i,sum(GS(i,:)),i,Dr(i),i,S(i),i,R(i)
         enddo !! i=1,nedotrunc
         write(iout,'(2(/),3X,">> Total G (in microsiemens) =",E13.6,/)') GTotal
         call mydealloc(   GS,'conductance_calc','GS'   ,nocolor=nocolor,verbose=mem_info) !+A:here
         call mydealloc(    G,'conductance_calc','G'    ,nocolor=nocolor,verbose=mem_info) !+A:here
         call mydealloc(    R,'conductance_calc','R'    ,nocolor=nocolor,verbose=mem_info) !+A:here
         call mydealloc(    S,'conductance_calc','S'    ,nocolor=nocolor,verbose=mem_info) !+A:here
         call mydealloc(   Dr,'conductance_calc','Dr'   ,nocolor=nocolor,verbose=mem_info) !+A:here
         call mydealloc(RREDO,'conductance_calc','RREDO',nocolor=nocolor,verbose=mem_info) !+A:here
         call mydealloc(SSEDO,'conductance_calc','SSEDO',nocolor=nocolor,verbose=mem_info) !+A:here
      enddo elec

      call mydealloc(  rr      ,'conductance_calc','rr'  ,nocolor=nocolor,verbose=mem_info) !+A:here
      call mydealloc(rEDO      ,'conductance_calc','rEDO',nocolor=nocolor,verbose=mem_info) !+A:here
      call mydealloc(sEDO      ,'conductance_calc','sEDO',nocolor=nocolor,verbose=mem_info) !+A:here
      call mydealloc(electrodes,'conductance_calc','electrodes',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:process_cond_par
      call mydealloc(cv,'conductance_calc','cv',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:EDOs
      call mydealloc(ordD2Acedo,'conductance_calc','ordD2Acedo',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:EDOs
      call mydealloc(D2,'conductance_calc','D2',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:EDOs
      call mydealloc(over,'conductance_calc','over',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:proc_dip
      deallocate(ato2basis) !+A:basis2atom_map

      call partial_time(ntime,'Conductance computation spent')

    end subroutine conductance_calc

end module conductance
