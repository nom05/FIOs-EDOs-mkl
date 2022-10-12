! ** ./FIOs-EDOs-mkl.f90 >> Computes FIOs/EDOs decomposition analysis. See references and the theoretical background in the manual.
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program FIOsEDOs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use omp_lib
  use memory_use
  use commonmod      ,only: nproc,nprocmkl,Field,ifchk,nmo,textxyz,iout,fileopen,igau,partial_time,ntime &
                          , intcomparegt,extension,ifio,lcomponents,process_text,comptext,ncmpeqv,hmnumb &
                          , procnumb,nheavy,print_presentation,print_message,debuglabel,nop,lop,momap    &
                          , ilendbglbl,idat,ignu,nmotrunc,ffiodat,fgnuplot,ffchkdef,mymemory,version  &
                          , myrealmemory,sqrccutoff,help_legend &
                          , squarec_cutoff              & !! loptions( -6) !  1
                          , lconduc                     & !! loptions( -5) !  2
                          , edo                         & !! loptions( -4) !  3
                          , truncate_dmo                & !! loptions( -3) !  4
                          , spec_comps                  & !! loptions( -2) !  5
                          , lfield                      & !! loptions( -1) !  6
                          , overwrite_output            & !! loptions(  0) !  7
                          , compute_dmo                 & !! loptions(  1) !  8
                          , dynamic_alpha               & !! loptions(  2) !  9
                          , test_deriv                  & !! loptions(  3) ! 10
                          , read_mult_mat               & !! loptions(  4) ! 11
                          , stealth                     & !! loptions(  5) ! 12
                          , nocolor                     & !! loptions(  6) ! 13
                          , debug,debuglabel,ilendbglbl & !! loptions(  7) ! 14
                          , output_file                 & !! loptions(  8) ! 15
                          , skip_fio_fchk               & !! loptions(  9) ! 16
                          , mk_gnuplot_scr              & !! loptions( 10) ! 17
                          , estimate_mem                & !! loptions( 11) ! 18
                          , print_large_matrices        & !! loptions( 12) ! 19
                          , print_pops_only             & !! loptions( 13) ! 20
                          , print_EF_pops               & !! loptions( 14) ! 21
                          , mem_info                    & !! loptions( 15) ! 22
                          , print_cond_mat             !& !! loptions( 16) ! 23

  use gen_multpol_mat,only: proc_dip
  use conductance    ,only: process_cond_par,conductance_calc

!------------------------------------------------------------------------!
!! libSUFR !! libSUFR !! libSUFR !! libSUFR !! libSUFR !! libSUFR !!!!!
  use SUFR_constants ,only: set_SUFR_constants_environment
  use SUFR_getopt    ,only: getopt_t, getopt_long, longoption,&
&                           optarg, getopt_long_help
  use SUFR_text      ,only: int2str,uppercase,dbl2str
  use SUFR_kinds     ,only: double,dbl
  use SUFR_system    ,only: warn,quit_program_error,quit_program_warning,quit_program
!------------------------------------------------------------------------!

  implicit none

  integer  (kind =   4)                              :: nargs,lstatus,jproc
  integer  (kind =   4)                              :: i,j !! tmp
  integer  (kind =   4),allocatable,dimension(:  )   :: motemp,motemp2
  real     (kind =   8)                              :: fieldinit=0._dbl
  character(len  =   1)                              :: option
  character(len  =1000)                              :: label,text
  character(len  = 105)                              :: fout,input,input2,input3
  character(len  =   5),parameter  ,dimension(2  )   :: prop = (/'Alpha','Beta '/)
  character(len  =   1),parameter  ,dimension(-1:1)  :: signo = (/'-','0','+'/)
  logical                                            :: lheavy=.FALSE. &
                                                      , set_OMP=.FALSE.,set_MKL=.FALSE.,lchk=.FALSE.,fio=.FALSE.

  interface
     integer function mkl_get_max_threads()
     end     function mkl_get_max_threads
     subroutine       mkl_set_num_threads(i)
       integer(kind=4)                                       :: i
     end subroutine   mkl_set_num_threads
  end interface

  interface
     subroutine       FIOs
     end subroutine   FIOs
     subroutine       EDOs
     end subroutine   EDOs
     subroutine       d2dmo
     end subroutine   d2dmo
     subroutine       dip2dipmo
     end subroutine   dip2dipmo
     subroutine       inter_c2ci
     end subroutine   inter_c2ci
     subroutine       readfchk2
     end subroutine   readfchk2
     subroutine       basis2atom_map
     end subroutine   basis2atom_map
     subroutine       readfchk(l)
       logical                                               :: l
     end subroutine   readfchk
     subroutine       read_pol(i,l,iou)
       integer(kind=4)                                       :: i,iou
       logical                                               :: l
     end subroutine   read_pol
     subroutine       readg09out(i,f)
       integer(kind=4)                                       :: i
       real   (kind=8)                                       :: f
     end subroutine   readg09out
     subroutine       indexxabs(B,matsort,nA,l,n,d)
       integer(kind=4)                          ,intent( in) :: nA
       integer(kind=4),            dimension(nA),intent(out) :: matsort
       integer(kind=4),optional                 ,intent( in) :: n
       real   (kind=8),            dimension(nA),intent( in) :: B
       logical        ,optional                 ,intent( in) :: d,l
     end subroutine   indexxabs
  end interface

!! * En el despiece por MO (incluso FIOs) poner un porcentaje. Como hay pos. y neg., separarlos.
!! * Determinar uso de memoria. Administración correcta de memoria.
!! * Cálculos unrestricted.
!! * Pseudopotenciales.

! Set up the longopts struct to define the valid options: short option, long option, argument (0/1), short description:
  type(getopt_t) :: longopts(28) = [ &
       getopt_t('h', 'help'     ,0, 'Print help.                                         [CM]')   & !!  1
      ,getopt_t('a', 'calc-dmo' ,0, 'alpha DMO will be computed instead of being read.   [FI]')   & !!  2
      ,getopt_t('C', 'conduc'   ,1, 'Computation of EDOs + conductance.                  [CO]**') & !!  3
      ,getopt_t('c', 'calc'     ,1, 'Choose components to be calculated (FIOs).          [FI]')   & !!  4
      ,getopt_t('D', 'dynamic'  ,0, 'Calculation of Dynamic alpha FIOs.                  [FI]')   & !!  5
      ,getopt_t('d', 'debug'    ,0, 'Print DEBUG messages.                               [CM]')   & !!  6
      ,getopt_t('E', 'edo'      ,1, 'Computation of EDOs. Include 2nd fchk file.         [ED]**') & !!  7
      ,getopt_t('e', 'estimate' ,0, 'Estimate memory only. Stop before computation.      [CM]')   & !!  8
      ,getopt_t('F', 'fchk'     ,1, 'fchk file to be employed w/ or w/o extension.       [CM]**') & !!  9
      ,getopt_t('f', 'field'    ,1, 'Set field by hand.                                  [FI]')   & !! 10
      ,getopt_t('G', 'gau'      ,1, 'G09 log file to be employed w/ or w/o extension.    [FI]**') & !! 11
      ,getopt_t('g', 'mkgnuplot',0, 'Make Gnuplot scripts to represent props vs # MO.    [FI]')   & !! 12
      ,getopt_t('i', 'print-mat',0, 'Print conductance matrices.                         [CO]')   & !! 13
      ,getopt_t('m', 'mem-info' ,0, 'Print memory info for each array (de)allocation.    [CM]')   & !! 14
      ,getopt_t('N', 'nmkl'     ,1, 'Number of threads for Intel(R) MKL libraries.       [CM]')   & !! 15
      ,getopt_t('n', 'no-color' ,0, 'Do not use colors on screen.                        [CM]')   & !! 16
      ,getopt_t('O', 'print-occ',0, 'Print def occup./MO contrib. to EDOs/FIOs,respectiv.[CM]')   & !! 17
      ,getopt_t('o', 'output'   ,1, 'Change name of the output file (w/o extension)      [CM]')   & !! 18
      ,getopt_t('P', 'proc'     ,1, 'Number of threads for OpenMP.                       [CM]')   & !! 19
      ,getopt_t('p', 'print-mat',0, 'Print large matrices.                               [CM]')   & !! 20
      ,getopt_t('r', 'read-mulm',0, 'Read multipole matrices from G09 log file.          [FI]')   & !! 21
      ,getopt_t('S', 'skip-ffch',0, 'Skip FIOs/EDOs fchk files writing                   [CM]')   & !! 22
      ,getopt_t('s', 'stealth'  ,0, 'Print on screen errors and warnings only.           [CM]')   & !! 23
      ,getopt_t('T', 'trunc-dmo',1, 'Truncate DMO with a set of MOs.                     [CM]')   & !! 24
      ,getopt_t('t', 'test'     ,0, 'Perform some tests with DM derivatives and dip. mat.[FI]')   & !! 25
      ,getopt_t('u', 'cutoff'   ,1, 'Square(c) cut-off in percentage.                    [CM]')   & !! 26
      ,getopt_t('W', 'write-pop',0, 'Write deformation orbital occupations and stop.     [FI]')   & !! 27
      ,getopt_t('w', 'overwrite',0, 'Overwrite output files.                             [CM]')   & !! 28
                                   ]
      fout(:) = ''
  ffchkdef(:) = ''
   ffiodat(:) = ''
  fgnuplot(:) = ''
     input(:) = ''
    input2(:) = ''

  call init_memory_counter

  call set_SUFR_constants_environment()
  nargs = command_argument_count()
  do   ! scan all the command-line parameters
       ! getopt_long() returns a single character" ">","!",".", or the short-option character (e.g. "a" for -a).
       !   It also sets two 'global' variables through the SUFR_getopt module:
       !   - longOption:  the full option (e.g. "-a" or "--all") including the dashes
       !   - optArg:      the argument following the option (if required and present)
     option = getopt_long(longopts)
       ! Do different things depending on the option returned:
     select case(option)
       case('>')  ! Last parameter
          if (nargs.EQ.0) then
             call warn('No arguments found!!',1,nocolor=nocolor)
             call help_legend
             call getopt_long_help(longopts)  ! No parameters found - print short help
             call quit_program_warning("Exiting ...",0,nocolor)
          endif !! (nargs.EQ.0) then
          exit
       case('!')  ! Unknown option (starting with "-" or "--")
          label = 'unknown option:  '//trim(optarg)//'. Use --help for a list of valid options'
          call quit_program_error(trim(label),1,nocolor)
       case('n')
          debuglabel  = '  (DG) '
          ilendbglbl  = 8
          nocolor     = .TRUE.
       case('s')
          stealth = .TRUE.
          debug   = .FALSE.
       case('c')
          read (optarg,*,iostat=lstatus) text
          if(len_trim(text).EQ.len(text)-5) call quit_program_error("Components specification too long!",1,nocolor)
          if(lstatus.NE.0)call quit_program_error("Problem with specification of comp. to be calculated!",1,nocolor)
          call   process_text(text)
          nop = 0
          do i = 1,size(lcomponents)
             if (lcomponents(i)) then
                nop      = nop+1
                lop(nop) = i
             endif !! (lcomponents(i)) then
          enddo !! i = 1,size(lcomponents)
          if (nop.EQ.0) then
             call quit_program_error("Problem with specification of comp.!",1,nocolor)
          else
             spec_comps = .TRUE.
          endif !! (nop.EQ.0) then
          if (debug) PRINT *,debuglabel(:ilendbglbl),nop,'=# comps - Components to be calculated: ' &
                                                        ,comptext(lop(:nop))
       case('D')
          dynamic_alpha = .TRUE.
       case('f')
          read (optarg,*,iostat=lstatus) fieldinit
          if(lstatus.NE.0) call quit_program_error("Wrong electric field specificied!",1,nocolor)
          lfield      = .TRUE.
       case('F')
          input(:) = ''
          read (optarg,*,iostat=lstatus) input
          if (len_trim(input ).EQ.len(input )-5) call quit_program_error("Fchk file name too long!",1,nocolor)
          if (debug) PRINT *,debuglabel(:ilendbglbl),'input(fchk): +',trim(input),'+'
          if(lstatus.NE.0) call quit_program_error("Unexpected problem while reading fchk file name",1,nocolor)
       case('G')
          if (edo) call quit_program_error("Simultaneous FIOs and EDOs computations are incompatible",1,nocolor)
          fio = .TRUE.
          input2(:) = ''
          read (optarg,*,iostat=lstatus) input2
          if (len_trim(input2).EQ.len(input2)-5) call quit_program_error("G09 log file name too long!",1,nocolor)
          if (debug) PRINT *,debuglabel(:ilendbglbl),'input(gau ): +',trim(input2),'+'
          if(lstatus.NE.0) call quit_program_error("Unexpected problem while reading G09 log file name",1,nocolor)
       case('r')
          read_mult_mat = .TRUE.
       case('a')
          compute_dmo = .TRUE.
       case('N')
          read (optarg,*,iostat=lstatus) nprocmkl
          if(lstatus.NE.0) call quit_program_error("Wrong # processors for Intel(R) MKL libraries!",1,nocolor)
          set_MKL=.TRUE.
       case('P')
          read (optarg,*,iostat=lstatus) nproc
          if(lstatus.NE.0) call quit_program_error("Wrong # processors for OpenMP!",1,nocolor)
          set_OMP=.TRUE.
       case('T')
          truncate_dmo = .TRUE.
          read (optarg,*,iostat=lstatus) label
          lheavy = index(uppercase(trim(label)),'HEAVY').GT.0
          if (debug) PRINT *,debuglabel(:ilendbglbl),'lheavy=',lheavy
          if (.NOT.lheavy) then
             call    hmnumb(trim(label),nmotrunc,jproc)
             if (debug) PRINT *,debuglabel(:ilendbglbl),'# Group of values: ',trim(int2str(jproc)),&
                                '; # orbs specified: ',trim(int2str(nmotrunc))
             if(allocated(motemp))call quit_program_error("New MO mapping array unexpectedly allocated",1,nocolor)
             call myalloc(motemp,nmotrunc,'main','motemp',nocolor=nocolor,verbose=mem_info) !+D:here
             call    procnumb(trim(label),len_trim(label)+1,nmotrunc,jproc,motemp)
          endif !! (.NOT.lheavy) then
       case('t')
          test_deriv = .TRUE.
       case('w')
          overwrite_output = .TRUE.
       case('d')
          debug = .TRUE.
       case('o')
          input3(:) = ''
          read (optarg,*,iostat=lstatus) input3
          if (len_trim(input3).EQ.len(input3)-5) call quit_program_error("Output file name too long!",1,nocolor)
          if (debug) PRINT *,debuglabel(:ilendbglbl),'input(mat ): +',trim(input3),'+'
          if(lstatus.NE.0) call quit_program_error("Unexpected problem while reading output file name",1,nocolor)
          output_file    = .TRUE.
       case('S')
          skip_fio_fchk  = .TRUE.
       case('g')
          mk_gnuplot_scr = .TRUE.
       case('E')
          if (fio) call quit_program_error("Simultaneous FIOs and EDOs computations are incompatible",1,nocolor)
          edo = .TRUE.
          input2(:) = ''
          read (optarg,*,iostat=lstatus) input2
          if (len_trim(input2).EQ.len(input2)-5) call quit_program_error("2nd fchk file name too long!",1,nocolor)
          if (debug) PRINT *,debuglabel(:ilendbglbl),'input(2fch): +',trim(input2),'+'
          if(lstatus.NE.0)call quit_program_error("Unexpected problem while reading 2nd fchk file name",1,nocolor)
       case('C')
          if (fio) call quit_program_error("Simultaneous FIOs and EDOs computations are incompatible",1,nocolor)
          edo     = .TRUE.
          lconduc = .TRUE.
          input2(:) = ''
            text(:) = ''
          read (optarg,*,iostat=lstatus) text
          if (len_trim(text).EQ.len(text)-5) call quit_program_error("Argument text too long",1,nocolor)
          if(lstatus.NE.0)call quit_program_error("Unexpected problem while reading opts. for conductance",1,nocolor)
          if (debug) PRINT *,debuglabel(:ilendbglbl),'arg text for conductance: +',trim(text),'+'
          label(:) = ''
          j = 0
          do i = 1,len_trim(text)
             if (trim(text(i:i)).EQ.':') then
                j = 0
                if (debug) PRINT *,debuglabel(:ilendbglbl),'         variable: +',trim(label),'+'
                call process_cond_par(label,.FALSE.,input2)
                label(:) = ''
                cycle
             elseif (i.EQ.len_trim(text)) then
                j = j+1
                label(j:j) = trim(text(i:i))
                if (debug) PRINT *,debuglabel(:ilendbglbl),'         variable: +',trim(label),'+'
                call process_cond_par(label,.TRUE.,input2)
                label(:) = ''
             endif !! (trim(text(i:i)).EQ.':') then
             j = j+1
             label(j:j) = trim(text(i:i))
          enddo
       case('e')
          estimate_mem         = .TRUE.
       case('p')
          print_large_matrices = .TRUE.
       case('O')
          print_EF_pops        = .TRUE.
       case('W')
          print_pops_only      = .TRUE.
       case('m')
          mem_info             = .TRUE.
       case('i')
          print_cond_mat       = .TRUE.
       case('u')
          read (optarg,*,iostat=lstatus) sqrccutoff
          sqrccutoff = sqrccutoff/100._dbl
          if(lstatus.NE.0) call quit_program_error("Wrong square c cut-off!",1,nocolor)
          squarec_cutoff       = .TRUE.
       case('h')
          write(*,'(/,6X,"HELP of FIOs-EDOs v.",A," - Intel(R) MKL version",2(/))') trim(version)
          call help_legend
          call getopt_long_help(longopts)
          write(*,'(/,"Conductance input format:")')
          PRINT *,'  * "f=file"       : 2nd fchk file.'
          PRINT *,'  * "d=real_number": electrode distance in A.'
          PRINT *,'  * "e1=1-3,5"     : set of atoms for electrode 1.'
          PRINT *,'  * "e2=12-20,22"  : set of atoms for electrode 2.'
          PRINT *,'  * "print"        : print conductance matrices.'
          PRINT *,'  * Arguments separated by colon(:) and all the text in quotes.'
          write(*,'(/,"Computational details:")')
          PRINT *,'  * FIOs: fchk and G09 log files are mandatory.'
          PRINT *,'  * EDOs: both fchk files are mandatory.'
          PRINT *,'  * "#p"        -> to obtain correct calculation type (mandatory).'
          PRINT *,'  * iop33(3=1)  -> to read multipole matrices from G09 log file (optional).'
          PRINT *,'  * iop33(10=2) -> mandatory writing of density matrices derivatives.'
          PRINT *,'  * field=read  -> mandatory field to perform numerical derivative of alpha.'
          PRINT *,'  * We recommend: "integral=ultrafinegrid" for DFT calculations.'
          PRINT *,'                  "scf=(conver=11,xqc)" for HF/DFT calculations.'
          write(*,'(/,"Example: G09 input file - Dynamic/static Alpha + static Beta_ZZZ")')
          write(*,'(/,   10("+"),1X,"START of input file")')
          write(*,'(2X,A)') '%oldchk=HF.chk'
          write(*,'(2X,A)') '%chk=HF_f0.chk'
          write(*,'(2X,A)') '#p hf chkbas geom=check guess=tcheck scf=(conver=11,xqc) nosym'
          write(*,'(2X,A)') 'polar IOp33(10=2,3=1) cphf=rdfreq'
          write(*,*       ) 
          write(*,'(2X,A)') 'comment'
          write(*,*       ) 
          write(*,'(2X,A)') '0 1'
          write(*,*       ) 
          write(*,'(2X,A)') '0.78992'
          write(*,*       ) 
          write(*,'(2X,A)') '--Link1--'
          write(*,'(2X,A)') '%oldchk=HF.chk'
          write(*,'(2X,A)') '%chk=HF_fz+.chk'
          write(*,'(2X,A)') '#p hf geom=allcheck guess=tcheck chkbasis scf=(conver=11,xqc)'
          write(*,'(2X,A)') 'polar IOp33(10=2) field=read'
          write(*,*       ) 
          write(*,'(2X,A)') '0. 0. 0.001'
          write(*,*       ) 
          write(*,'(2X,A)') '--Link1--'
          write(*,'(2X,A)') '%oldchk=HF.chk'
          write(*,'(2X,A)') '%chk=HF_fz-.chk'
          write(*,'(2X,A)') '#p hf geom=allcheck guess=tcheck chkbasis scf=(conver=11,xqc)'
          write(*,'(2X,A)') 'polar IOp33(10=2) field=read'
          write(*,*       ) 
          write(*,'(2X,A)') '0. 0. -.001'
          write(*,'(/,10("+"),1X,"END   of input file")')
          call quit_program("The program will stop. Rerun it without help")
 !      case('.')  ! Parameter is not an option (i.e., it doesn't start with "-" or "--")
 !         np = np+1
 !         files(np)%filename = trim(optarg)
 !         files(np)%fileunit = input+np
 !         write(charint,'(I6)') np
 !         write(*      ,'(A)') ' o File name'//        &
 !&                trim(charint(verify(charint,' '):6)) &
 !&                                          //' specified: '//trim(files(np)%filename)
       case default
          label(:) = ''
          label    = 'Valid option unhandled: '//trim(longoption)
          call quit_program_error(trim(label),1,nocolor)
     end select
  enddo ! scan all the command-line parameters
  call init_memory_counter
  if (.NOT.stealth) then
     if (.NOT.nocolor) then
        write(*,'(/,6X,"Computation of FIELD-INDUCED/ELECTRON-DEFORMATION ORBITALS - Intel(R) MKL version",&
                 &2(/),3X,"** Legend: (",A,") - Informative msg"   ,/,   &
                               &13X," (",A,") - Elapsed time msg"    )') &
                                achar(27)//'[32mII'//achar(27)//'[0m', &
                                achar(27)//'[34mTT'//achar(27)//'[0m'
        if (debug   ) write(*,'(13X," (",A,") - Debug msg"           )') achar(27)//'[90mDG'//achar(27)//'[0m'
        if (mem_info) write(*,'(13X," (",A,") - Verbose memory msg"  )') achar(27)//'[35mME'//achar(27)//'[0m'
        write(*,'(/)')
     else
        write(*,'(/,6X,"Computation of FIELD-INDUCED/ELECTRON-DEFORMATION ORBITALS - Intel(R) MKL version",&
                 &2(/),3X,"** Legend: (",A,") - Informative msg"    ,/,   &
                               &13X," (",A,") - Elapsed time msg"     )') &
                                               'II'                  , &
                                               'TT'                  
        if (debug   ) write(*,'(13X," (",A,") - Debug msg"            )')                'DG'
        if (mem_info) write(*,'(13X," (",A,") - Verbose memory msg"   )')                'ME'
        write(*,'(/)')
     endif !! (.NOT.nocolor) then
  endif !! (.NOT.stealth then

  if (debug.AND.mem_info) call print_message('d','!! Verbose memory info will be printed !!')
  if (debug.AND.stealth) then
     call warn('"DEBUG msgs." and "stealth mode" are both enabled. Disabling the 1st one',nocolor=nocolor)
     debug = .FALSE.
  else if (debug.AND..NOT.stealth)  then
     call print_message('d','!! DEBUG mode is ON !!')
  endif !! (debug) then
  if (edo) then
     if (fio) call quit_program_error('Simultaneous EDOs and FIOs computations are incompatible',1,nocolor)
     if (dynamic_alpha) call quit_program_error('Dynamic alpha is incompatible with EDOs',1,nocolor)
     call print_message('i','Computation of ELECTRON DEFORMATION ORBITALs (EDOs) is enabled')
     call warn('"--calc-dmo" is automatically set',1,nocolor=nocolor)
     call warn('"--read-mulm" is automatically unset',1,nocolor=nocolor)
     compute_dmo   = .TRUE.
     read_mult_mat = .FALSE.
     nop           = 2 !! We will save the (un)perturbed DMs wrt the electric field
  elseif (fio) then
     call        print_message('i','FIOs computations will be performed')
  endif !! (edo) then
  if (debug.AND.read_mult_mat) call print_message('d','Multipole matrices read from G09 log file')
  if (          dynamic_alpha) &
     call warn('Dynamic alpha FIOs will be computed. "--calc-dmo" is automatically unset',1,nocolor=nocolor)
  if (.NOT.edo.AND..NOT.spec_comps) &
     call warn('TIP: to avoid mem. fragmentation and a correct memory estimation, request the components &
               &you need with "--calc"',1,nocolor=nocolor)
  if (          overwrite_output) call warn('Previous output files will be overwritten',1,nocolor=nocolor)

  if (len(trim(input )).LT.1)         call quit_program_error('fchk file     is kept unset' ,1,nocolor)
  if (edo) then
     if (len(trim(input2)).LT.1)      call quit_program_error('2nd fchk file is kept unset' ,1,nocolor)
  else
     if (len(trim(input2)).LT.1)      call quit_program_error('G09  log file is kept unset' ,1,nocolor)
  endif !! (edo) then
  if (lfield.AND.fieldinit.EQ.0._dbl) call quit_program_error('Electric Field is kept null or unset',1,nocolor)
  if (dynamic_alpha) then
     compute_dmo = .FALSE.
     if (.NOT.spec_comps) then
        nop = 3 !! x,y,z
        lop = (/ (i, i = 1,nop) /)
     else
        if (nop.GT.3) then
           call warn('Dynamic polarizabilities are uniquely available for alpha',1,nocolor=nocolor)
           nop = 3 !! x,y,z
        endif !!(nop.GT.3) then
        do i = nop,1,-1
           if (lop(i).GT.3) then
              call warn('Dynamic polarizabilities are uniquely available for alpha',1,nocolor=nocolor)
              nop = nop-1
           endif !! (lop(i).GT.3) then
        enddo !! i = 1,nop
     endif !! (.NOT.spec_comps) then
  elseif (edo) then
     if (spec_comps) call warn('"--calc" is unnecessary with EDOs',1,nocolor=nocolor)
     spec_comps = .FALSE.
     nop    = 2
     lop(1) = 0
  else
     if (.NOT.spec_comps) then
        nop = size(lcomponents)
        lop = (/ (i, i = 1,nop) /)
     endif !! (.NOT.spec_comps) then
  endif !! (compute_dmo.AND.dynamic_alpha) then
  if (estimate_mem) call print_message('i','Memory allocation will be estimated only')
  if (print_pops_only.AND.edo) then
     print_EF_pops = .TRUE.
     call warn('The program will stop after printing deformation orbital occupations' &
               ,1,nocolor=nocolor)
  elseif (print_pops_only.AND.fio) then
     call warn('FIOs and "--write-pop" are incompatible. The latter will be obviated' &
               ,1,nocolor=nocolor)
     print_pops_only = .FALSE.
  endif !! (print_pops_only.AND.edo) then
  if (print_large_matrices) then
     call warn('Large matrices will be printed. Use with debug purposes only to avoid huge output files' &
               ,1,nocolor=nocolor)
     write(*,*)
  endif !! (print_large_matrices) then
  if (squarec_cutoff) then
     label(:) = ''
     label    = "Squared c coefficients will be printed using the cut-off "//dbl2str(sqrccutoff*100._dbl,2)//"%"
     call print_message('i',trim(label))
  endif !! (squarec_cutoff) then
  call            partial_time(ntime,'') !! First time to take wall time

  if (set_OMP.AND..NOT.set_MKL) then
     call warn('Intel(R) MKL threads not set, employing the same number of OpenMP threads',1,nocolor=nocolor)
     nprocmkl = nproc
  endif !! (set_OMP) then
  if (nproc   .GT.OMP_GET_MAX_THREADS()) then
     call    intcomparegt('nproc',nproc,OMP_GET_MAX_THREADS(),'OMP_GET_MAX_THREADS',label)
     call warn(trim(label),1,nocolor=nocolor)
     nproc    = OMP_GET_MAX_THREADS()
  endif !! (nproc.GT.OMP_GET_MAX_THREADS()) then
  if (nprocmkl.GT.MKL_GET_MAX_THREADS()) then
     call    intcomparegt('nmkl',nprocmkl,MKL_GET_MAX_THREADS(),'MKL_GET_MAX_THREADS',label)
     call warn(trim(label),1,nocolor=nocolor)
     nprocmkl = MKL_GET_MAX_THREADS()
  endif !! (nproc.GT.OMP_GET_MAX_THREADS()) then
  label(:) = ''
  label    = "Running      program with "//int2str(nproc   )//" threads"
  call print_message('i',trim(label))
  label    = "Running Intel(R) MKL with "//int2str(nprocmkl)//" threads"
  call print_message('i',trim(label))
  call OMP_SET_NUM_THREADS(nproc   )
  call MKL_SET_NUM_THREADS(nprocmkl)

  text(:) = ''
  call fileopen(trim(input ),trim(extension(1)),text,ifchk)
  label(:) = ''
  label    = "fchk file name (w/o extension): "//trim(text)
  call print_message('i',trim(label))
  if (output_file) then
     fout(:)     = trim(input3(:))
     ffchkdef(:) = trim(input3(:))
      ffiodat(:) = trim(input3(:))
     fgnuplot(:) = trim(input3(:))
  else
     fout(:)     = trim(text(:))
     ffchkdef(:) = trim(text(:))
      ffiodat(:) = trim(text(:))
     fgnuplot(:) = trim(text(:))
  endif !! (output_file) then
  label    = "Output file                   : "//trim(fout)//'.'//trim(extension(3))
  if (.NOT.overwrite_output) then
     inquire(file=trim(fout)//'.'//trim(extension(3)),exist=lchk)
     if (lchk) call quit_program_error('Specified output file will not be overwritten',1,nocolor)
  endif !! (.NOT.overwrite_output) then
  
  call print_message('i',trim(label))
  call readfchk(.TRUE.) ! Read nmo,nfb,nel,nato
!
! Create (if not exist) and check MO map
!
  if (allocated(momap)) call quit_program_error('The map of MOs is allocated improperly',1,nocolor)
  if (truncate_dmo) then
     if (.NOT.lheavy) then
        if (size(motemp).GT.nmo) then
           label = 'Dimension of MO map ('//trim(int2str(size(motemp)))//') > # MO ('//trim(int2str(nmo))//')'
           call quit_program_error(trim(label),1,nocolor)
        endif !! (size(motemp).GT.nmo) then
        do i = 1,nmotrunc
           if (motemp(i).GT.nmo.OR.motemp(i).LE.0) then
              label = 'Value of MO map ('//trim(int2str(motemp(i)))//') wrong --> # MO ('//trim(int2str(nmo))//')'
              call quit_program_error(trim(label),1,nocolor)
           endif !! (motemp.GT.nmo)
        enddo !! i = 1,nmotrunc
        call myalloc(momap,nmotrunc,'main','momap',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:EDOs&FIOs
        call myalloc(motemp2,nmotrunc,'main','motemp2',nocolor=nocolor,verbose=mem_info) !+D:here
        call       indexxabs(dble(motemp),motemp2,nmotrunc,.FALSE.,nproc,debug)
        do concurrent (i = 1:nmotrunc)
           j        = nmotrunc+1-i
           momap(j) = motemp(motemp2(i))
        enddo
        if (allocated(motemp )) call mydealloc(motemp ,'main','motemp' ,nocolor=nocolor,verbose=mem_info) !+A:here
        if (allocated(motemp2)) call mydealloc(motemp2,'main','motemp2',nocolor=nocolor,verbose=mem_info) !+A:here
     else  !! Skip "1s" MOs from heavy atoms.
        nmotrunc = nmo-2*nheavy
        call myalloc(momap,nmotrunc,'main','momap',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:EDOs&FIOs
        momap = (/ (i, i = nheavy+1,nmo-nheavy) /)
     endif !! (.NOT.lheavy) then
  else
     call myalloc(momap,nmo,'main','momap',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:EDOs&FIOs
     nmotrunc = nmo
     momap = (/ (i, i = 1,nmotrunc) /)
  endif !! (.NOT.truncate_dmo) then
  if (debug)  PRINT *,debuglabel(:ilendbglbl),'Final map of MOs:',momap,'(',size(momap),')'

  call fileopen(trim(fout  ),trim(extension(3)),text,iout ,lexist=.FALSE.,sttsfl='UNKNOWN' &
                            ,overwrite=overwrite_output)

  call print_presentation(trim(input),trim(input2),trim(fout)//'.'//trim(extension(3)))

  call mymemory                !! Memory estimation. Without 'spec_comps' d and dmo are estimated for one component
  if (estimate_mem) call quit_program("The program will stop now")

  call readfchk(.FALSE.)       !! Full reading
  call partial_time(ntime,'Fchk file read in')

  call basis2atom_map          !! Atom map of basis functions

  text(:) = ''
  label(:) = ''
  if (edo) then
     call fileopen(trim(input2),trim(extension(1)),text,igau )
     label    = '2nd fchk file name (w/o extension): '//trim(text)
     call print_message('i',trim(label))
     call readfchk2
     call        partial_time(ntime,'2nd fchk file read in')
  else
     call fileopen(trim(input2),trim(extension(2)),text,igau )
     label    = 'G09  file name (w/o extension): '//trim(text)
     call print_message('i',trim(label))
     call     readg09out(0,fieldinit)
     !! Rewind if zero---^
     call        partial_time(ntime,'G09 log file read in')
  endif !! (edo) then

  if (lfield) then
     label(:) = ''
     label    = 'Electric field strength: '//dbl2str(Field,4)
     call print_message('i',trim(label))
     if (dynamic_alpha.OR.maxval(ncmpeqv).LE.3) &
             call warn('Electric field is not used with alpha',1,nocolor=nocolor)
  endif !! (.NOT.lfield) then

  if (.NOT.spec_comps.OR..NOT.edo) call myrealmemory   !! Memory use

  call inter_c2ci                   !! Invert MO coefficients matrix

                                    !! Read or compute Multipole/overlap array (built from BFs):
  call proc_dip(read_mult_mat,debug,print_large_matrices,lconduc,nocolor,mem_info)
  if (.NOT.edo) call dip2dipmo      !! Compute Multipole matrices array (built from MOs), including truncated one

  call d2dmo                        !! Compute DMO if necessary, trunc. is not included because DMO is cp to A

  if (.NOT.edo) then
     write(iout,'(/,3X,       ">> MOL. PROPERTIES (from Gaussian log file):")')
     call print_message('i',"MOL. PROPERTIES (from Gaussian log file):")
     if (.NOT.stealth) write(*,*)
     call read_pol(igau,.FALSE.,iout)
     if (.NOT.stealth) write(*,*)
     call           partial_time(ntime,'MOL. PROPS spent')
  endif !! (.NOT.edo) then
  close(igau)

  if (edo) call EDOs

  if (edo.AND.lconduc) call conductance_calc

  if (fio) call FIOs

  write(*,*)

  call  partial_time(ntime,'Arrays deallocatation spent')
  if (mem_info) call mem_report(nocolor)

  close(ifchk)
  write(iout,*)
  call  partial_time(ntime,'',.TRUE.,iout)
  close(iout )

end program
