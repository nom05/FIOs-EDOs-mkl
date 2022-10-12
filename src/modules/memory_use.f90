! ** modules/memory_use.f90 >> Main module to allocate and save information about array allocation, etc.
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


module     memory_use
!
! Module based on "memory_use.f90" from CRYSTAL source code
!
  use SUFR_kinds         ,only: double,dbl
  implicit none
  integer( kind = 4 )        :: memory_used,max_mem,max_avail_mem,ncounta,ncountd
  interface     myalloc
    module procedure   IALLOC1
    module procedure   IALLOC2
    module procedure   RALLOC1
    module procedure   RALLOC2
    module procedure   RALLOC3
    module procedure   TALLOC1 !! For text or character array
  end interface myalloc
  interface     mydealloc
    module procedure IDEALLOC1
    module procedure IDEALLOC2
    module procedure RDEALLOC1
    module procedure RDEALLOC2
    module procedure RDEALLOC3
    module procedure TDEALLOC1 !! For text or character array
  end interface mydealloc
!
! 8 bit = 1 word = 1 byte
! 1 char = 1 byte
! Calculations use "double precision". To use "single precision", use
! SIZE_DOUBLE_OVER_SIZE_REAL=1
!
  integer( kind = 4 ),parameter,public  :: size_double_over_size_real = 2
  integer( kind = 4 ),parameter,private :: size_int                   = bit_size(size_int) / 8
  integer( kind = 4 ),parameter,private :: size_logic                 = size_int
  integer( kind = 4 ),parameter,private :: size_double                = size_int * size_double_over_size_real
  integer( kind = 4 ),parameter,private :: size_char                  = 1
  integer  ( kind =   4 ),allocatable,dimension(:)   :: iallocmem,ialloctype
  character(  len = 100 ),allocatable,dimension(:)   :: allocname,allocwhere,deallocwhere,allocmem
  logical                ,allocatable,dimension(:)   :: lallocsaved
  character(  len =   1 ),parameter  ,dimension(0:3) :: allchartype=(/ 'U' &  !! 0
                                                                      ,'I' &  !! 1
                                                                      ,'R' &  !! 2
                                                                      ,'T' &  !! 3
                                                                    /)
  character (  len =   2 ),dimension(5),parameter    :: unitsdec=(/'B ' ,'kB' ,'MB' ,'GB' ,'TB' /)
  character (  len =   3 ),dimension(5),parameter    :: unitsbin=(/'B  ','KiB','MiB','GiB','TiB'/)

!! Max dimensions for mem_report
  integer  ( kind =   4 )                            :: maxn=1,maxname=4,maxsize=4,maxtype=4,maxallwh=8,mxalldewh=10

  contains

  subroutine     init_memory_counter
    memory_used = 0
    max_mem     = 0
    ncounta     = 0
    ncountd     = 0
  end subroutine init_memory_counter

  subroutine     IALLOC1(array,ndim,routine_name,array_name,verbose,nocolor,lsave)
    implicit none
    integer   ( kind = 4 ),allocatable,dimension(:) :: array
    integer   ( kind = 4 ),intent(in)               :: ndim
    character             ,intent(in)               :: routine_name*(*),array_name*(*)
    integer   ( kind = 4 )                          :: ierr,isize
    logical,optional      ,intent(in)               :: verbose   ,nocolor   ,lsave
    logical                                         :: verbosedef,nocolordef,lsavedef

    verbosedef = .TRUE.
    nocolordef = .FALSE.
    lsavedef   = .FALSE.
    if (present(verbose)) verbosedef = verbose
    if (present(nocolor)) nocolordef = nocolor
    if (present(lsave)  ) lsavedef   = lsave
    call print_mem_text('IALLOC1 called by '//trim(array_name)//' in '//trim(routine_name) &
                       ,verbosedef,nocolordef)
    isize       = ndim*size_int
    memory_used = memory_used+isize
    allocate(array(ndim),STAT=ierr)
    call ProcAlloc(isize,ierr,nocolordef,verbosedef,lsavedef,routine_name,array_name,iarrtype=1)
    return
  end subroutine IALLOC1

  subroutine     IALLOC2(array,ndim1,ndim2,routine_name,array_name,verbose,nocolor,lsave)
    implicit none
    integer   ( kind = 4 ),allocatable,dimension(:,:) :: array
    integer   ( kind = 4 ),intent(in)                 :: ndim1,ndim2
    character             ,intent(in)                 :: routine_name*(*),array_name*(*)
    integer   ( kind = 4 )                            :: ierr,isize
    logical,optional      ,intent(in)                 :: verbose   ,nocolor   ,lsave
    logical                                           :: verbosedef,nocolordef,lsavedef

    verbosedef = .TRUE.
    nocolordef = .FALSE.
    lsavedef   = .FALSE.
    if (present(verbose)) verbosedef = verbose
    if (present(nocolor)) nocolordef = nocolor
    if (present(lsave)  ) lsavedef   = lsave
    call print_mem_text('IALLOC2 called by '//trim(array_name)//' in '//trim(routine_name) &
                       ,verbosedef,nocolordef)
    isize       = ndim1*ndim2*size_int
    memory_used = memory_used+isize
    allocate(array(ndim1,ndim2),STAT=ierr)
    call ProcAlloc(isize,ierr,nocolordef,verbosedef,lsavedef,routine_name,array_name,iarrtype=1)
    return
  end subroutine IALLOC2

  subroutine     RALLOC1(array,ndim,routine_name,array_name,verbose,nocolor,lsave)
    implicit none
    real      ( kind = double ),allocatable,dimension(:) :: array
    integer   ( kind =      4 ),intent(in)               :: ndim
    character                  ,intent(in)               :: routine_name*(*),array_name*(*)
    integer   ( kind =      4 )                          :: ierr,isize
    logical,optional      ,intent(in)                    :: verbose   ,nocolor   ,lsave
    logical                                              :: verbosedef,nocolordef,lsavedef

    verbosedef = .TRUE.
    nocolordef = .FALSE.
    lsavedef   = .FALSE.
    if (present(verbose)) verbosedef = verbose
    if (present(nocolor)) nocolordef = nocolor
    if (present(lsave)  ) lsavedef   = lsave
    call print_mem_text('RALLOC1 called by '//trim(array_name)//' in '//trim(routine_name) &
                       ,verbosedef,nocolordef)
    isize       = ndim*size_double
    memory_used = memory_used+isize
    allocate(array(ndim),STAT=ierr)
    call ProcAlloc(isize,ierr,nocolordef,verbosedef,lsavedef,routine_name,array_name,iarrtype=2)
    return
  end subroutine RALLOC1

  subroutine     RALLOC2(array,ndim1,ndim2,routine_name,array_name,verbose,nocolor,lsave)
    implicit none
    real      ( kind = double ),allocatable,dimension(:,:) :: array
    integer   ( kind =      4 ),intent(in)                 :: ndim1,ndim2
    character                  ,intent(in)                 :: routine_name*(*),array_name*(*)
    integer   ( kind =      4 )                            :: ierr,isize
    logical,optional      ,intent(in)                      :: verbose   ,nocolor   ,lsave
    logical                                                :: verbosedef,nocolordef,lsavedef

    verbosedef = .TRUE.
    nocolordef = .FALSE.
    lsavedef   = .FALSE.
    if (present(verbose)) verbosedef = verbose
    if (present(nocolor)) nocolordef = nocolor
    if (present(lsave)  ) lsavedef   = lsave
    call print_mem_text('RALLOC2 called by '//trim(array_name)//' in '//trim(routine_name) &
                       ,verbosedef,nocolordef)
    isize       = ndim1*ndim2*size_double
    memory_used = memory_used+isize
    allocate(array(ndim1,ndim2),STAT=ierr)
    call ProcAlloc(isize,ierr,nocolordef,verbosedef,lsavedef,routine_name,array_name,iarrtype=2)
    return
  end subroutine RALLOC2

  subroutine     RALLOC3(array,ndim1,ndim2,ndim3,routine_name,array_name,verbose,nocolor,lsave)
    implicit none
    real      ( kind = double ),allocatable,dimension(:,:,:) :: array
    integer   ( kind =      4 ),intent(in)                   :: ndim1,ndim2,ndim3
    character                  ,intent(in)                   :: routine_name*(*),array_name*(*)
    integer   ( kind =      4 )                              :: ierr,isize
    logical,optional      ,intent(in)                        :: verbose   ,nocolor   ,lsave
    logical                                                  :: verbosedef,nocolordef,lsavedef

    verbosedef = .TRUE.
    nocolordef = .FALSE.
    lsavedef   = .FALSE.
    if (present(verbose)) verbosedef = verbose
    if (present(nocolor)) nocolordef = nocolor
    if (present(lsave)  ) lsavedef   = lsave
    call print_mem_text('RALLOC3 called by '//trim(array_name)//' in '//trim(routine_name) &
                       ,verbosedef,nocolordef)
    isize       = ndim1*ndim2*ndim3*size_double
    memory_used = memory_used+isize
    allocate(array(ndim1,ndim2,ndim3),STAT=ierr)
    call ProcAlloc(isize,ierr,nocolordef,verbosedef,lsavedef,routine_name,array_name,iarrtype=2)
    return
  end subroutine RALLOC3

  subroutine     TALLOC1(array,ndim,routine_name,array_name,verbose,nocolor,lsave)
    implicit none
    character (  len =      * ),allocatable,dimension(:) :: array
    integer   ( kind =      4 ),intent(in)               :: ndim
    character                  ,intent(in)               :: routine_name*(*),array_name*(*)
    integer   ( kind =      4 )                          :: ierr,isize
    logical,optional      ,intent(in)                    :: verbose   ,nocolor   ,lsave
    logical                                              :: verbosedef,nocolordef,lsavedef

    verbosedef = .TRUE.
    nocolordef = .FALSE.
    lsavedef   = .FALSE.
    if (present(verbose)) verbosedef = verbose
    if (present(nocolor)) nocolordef = nocolor
    if (present(lsave)  ) lsavedef   = lsave
    call print_mem_text('TALLOC1 called by '//trim(array_name)//' in '//trim(routine_name) &
                       ,verbosedef,nocolordef)
    isize       = ndim*size_char
    memory_used = memory_used+isize
    allocate(array(ndim),STAT=ierr)
    call ProcAlloc(isize,ierr,nocolordef,verbosedef,lsavedef,routine_name,array_name,iarrtype=3)
    return
  end subroutine TALLOC1

  subroutine     IDEALLOC1(array,routine_name,array_name,verbose,nocolor,lsave)
    implicit none
    integer   ( kind = 4 ),allocatable,dimension(:) :: array
    character             ,intent(in)               :: routine_name*(*),array_name*(*)
    integer   ( kind = 4 )                          :: isize
    logical,optional      ,intent(in)               :: verbose   ,nocolor   ,lsave
    logical                                         :: verbosedef,nocolordef,lsavedef

    verbosedef = .TRUE.
    nocolordef = .FALSE.
    lsavedef   = .FALSE.
    if (present(verbose)) verbosedef = verbose
    if (present(nocolor)) nocolordef = nocolor
    if (present(lsave)  ) lsavedef   = lsave
    call print_mem_text('IDEALLOC1 called by '//trim(array_name)//' in '//trim(routine_name) &
                       ,verbosedef,nocolordef)
    isize       = size(array)*size_int
    memory_used = memory_used-isize
    deallocate(array)
    call ProcDEAlloc(isize,nocolor,verbose,routine_name,array_name,lsearch=lsavedef)
  end subroutine IDEALLOC1

  subroutine     IDEALLOC2(array,routine_name,array_name,verbose,nocolor,lsave)
    implicit none
    integer   ( kind = 4 ),allocatable,dimension(:,:) :: array
    character             ,intent(in)                 :: routine_name*(*),array_name*(*)
    integer   ( kind = 4 )                            :: isize
    logical,optional      ,intent(in)                 :: verbose   ,nocolor   ,lsave
    logical                                           :: verbosedef,nocolordef,lsavedef

    verbosedef = .TRUE.
    nocolordef = .FALSE.
    lsavedef   = .FALSE.
    if (present(verbose)) verbosedef = verbose
    if (present(nocolor)) nocolordef = nocolor
    if (present(lsave)  ) lsavedef   = lsave
    call print_mem_text('IDEALLOC2 called by '//trim(array_name)//' in '//trim(routine_name) &
                       ,verbosedef,nocolordef)
    isize       = size(array)*size_int
    memory_used = memory_used-isize
    deallocate(array)
    call ProcDEAlloc(isize,nocolor,verbose,routine_name,array_name,lsearch=lsavedef)
  end subroutine IDEALLOC2

  subroutine     RDEALLOC1(array,routine_name,array_name,verbose,nocolor,lsave)
    implicit none
    real      ( kind = double),allocatable,dimension(:) :: array
    character                 ,intent(in)               :: routine_name*(*),array_name*(*)
    integer   ( kind =     4 )                          :: isize
    logical,optional      ,intent(in)                   :: verbose   ,nocolor   ,lsave
    logical                                             :: verbosedef,nocolordef,lsavedef

    verbosedef = .TRUE.
    nocolordef = .FALSE.
    lsavedef   = .FALSE.
    if (present(verbose)) verbosedef = verbose
    if (present(nocolor)) nocolordef = nocolor
    if (present(lsave)  ) lsavedef   = lsave
    call print_mem_text('RDEALLOC1 called by '//trim(array_name)//' in '//trim(routine_name) &
                       ,verbosedef,nocolordef)
    isize       = size(array)*size_double
    memory_used = memory_used-isize
    deallocate(array)
    if (lsavedef) call print_mem_text('    (Saving array data ...)',verbose,nocolor)
    call ProcDEAlloc(isize,nocolor,verbose,routine_name,array_name,lsearch=lsavedef)
  end subroutine RDEALLOC1

  subroutine     RDEALLOC2(array,routine_name,array_name,verbose,nocolor,lsave)
    implicit none
    real      ( kind = double),allocatable,dimension(:,:) :: array
    character                 ,intent(in)                 :: routine_name*(*),array_name*(*)
    integer   ( kind =     4 )                            :: isize
    logical,optional      ,intent(in)                     :: verbose   ,nocolor   ,lsave
    logical                                               :: verbosedef,nocolordef,lsavedef

    verbosedef = .TRUE.
    nocolordef = .FALSE.
    lsavedef   = .FALSE.
    if (present(verbose)) verbosedef = verbose
    if (present(nocolor)) nocolordef = nocolor
    if (present(lsave)  ) lsavedef   = lsave
    call print_mem_text('RDEALLOC2 called by '//trim(array_name)//' in '//trim(routine_name) &
                       ,verbosedef,nocolordef)
    isize       = size(array)*size_double
    memory_used = memory_used-isize
    deallocate(array)
    call ProcDEAlloc(isize,nocolor,verbose,routine_name,array_name,lsearch=lsavedef)
  end subroutine RDEALLOC2

  subroutine     RDEALLOC3(array,routine_name,array_name,verbose,nocolor,lsave)
    implicit none
    real      ( kind = double),allocatable,dimension(:,:,:) :: array
    character                 ,intent(in)                   :: routine_name*(*),array_name*(*)
    integer   ( kind =     4 )                              :: isize
    logical,optional      ,intent(in)                       :: verbose   ,nocolor   ,lsave
    logical                                                 :: verbosedef,nocolordef,lsavedef

    verbosedef = .TRUE.
    nocolordef = .FALSE.
    lsavedef   = .FALSE.
    if (present(verbose)) verbosedef = verbose
    if (present(nocolor)) nocolordef = nocolor
    if (present(lsave)  ) lsavedef   = lsave
    call print_mem_text('RDEALLOC3 called by '//trim(array_name)//' in '//trim(routine_name) &
                       ,verbosedef,nocolordef)
    isize       = size(array)*size_double
    memory_used = memory_used-isize
    deallocate(array)
    call ProcDEAlloc(isize,nocolor,verbose,routine_name,array_name,lsearch=lsavedef)
  end subroutine RDEALLOC3

  subroutine     TDEALLOC1(array,routine_name,array_name,verbose,nocolor,lsave)
    implicit none
    character (  len =     * ),allocatable,dimension(:) :: array
    character                 ,intent(in)               :: routine_name*(*),array_name*(*)
    integer   ( kind = 4 )                              :: isize
    logical,optional      ,intent(in)                   :: verbose   ,nocolor   ,lsave
    logical                                             :: verbosedef,nocolordef,lsavedef

    verbosedef = .TRUE.
    nocolordef = .FALSE.
    lsavedef   = .FALSE.
    if (present(verbose)) verbosedef = verbose
    if (present(nocolor)) nocolordef = nocolor
    if (present(lsave)  ) lsavedef   = lsave
    call print_mem_text('TDEALLOC1 called by '//trim(array_name)//' in '//trim(routine_name) &
                       ,verbosedef,nocolordef)
    isize       = size(array)*size_char
    memory_used = memory_used-isize
    deallocate(array)
    call ProcDEAlloc(isize,nocolor,verbose,routine_name,array_name,lsearch=lsavedef)
  end subroutine TDEALLOC1

  subroutine     ProcAlloc(isize,ierr,nocolor,verbose,lsave,routine_name,array_name,iarrtype,perform_sum)
    use SUFR_system        ,only                       : quit_program_error,warn
    use SUFR_text          ,only                       : int2str,dbl2str
    implicit none
    integer   ( kind =   4 ),intent(in)               :: isize,ierr
    integer   ( kind =   4 )                          :: nsize,iunits
    integer   ( kind =   4 ),intent(in),optional      :: iarrtype
    integer   ( kind =   4 )                          :: iarrtype_def
    character               ,intent(in),optional      :: routine_name*(*),array_name*(*)
    character (  len =  30 )                          :: rouname_def,arrayname_def
    logical                 ,intent(in)               :: nocolor,verbose
    logical                 ,intent(in),optional      :: lsave,perform_sum
    logical                                           :: lsave_def,perform_sum_def
    real      ( kind =   8 )                          :: unitchgdec,unitchgbin
    integer   ( kind =   4 ),allocatable,dimension(:) :: itmparr
    character (  len = 100 ),allocatable,dimension(:) :: temptxt
    logical                 ,allocatable,dimension(:) :: ltemp

    lsave_def       = .FALSE.
    perform_sum_def = .FALSE.
    iarrtype_def    = 0
    if (present(routine_name)) then
       if (len_trim(routine_name).GE.len(rouname_def)) &
               call warn('Too long subroutine name will be truncated',nocolor=nocolor)
       rouname_def   = routine_name(:)
    else
       rouname_def   = 'UNKNOWN ROUTINE'
    endif !! (present(routine_name)) then
    if (present(array_name)  ) then
       if (len_trim(array_name).GE.len(arrayname_def)) &
               call warn('Too long array name will be truncated',nocolor=nocolor)
       arrayname_def =   array_name(:)
    else
       arrayname_def = 'UNKNOWN ARRAY'
    endif !! (present(array_name)  ) then
    if (present(iarrtype)   ) iarrtype_def    = iarrtype
    if (present(lsave)      ) lsave_def       = lsave
    if (present(perform_sum)) perform_sum_def = perform_sum
    if (perform_sum_def) memory_used = memory_used+isize
    max_mem = max(max_mem,memory_used)
    if (ierr.NE.0) then                                                                            
       ! Error codes are compiler dependent
 iallc:if (ierr.EQ.672.OR.ierr.EQ.718.OR.ierr.EQ.722.OR.ierr.EQ.724.OR.ierr.EQ.727.OR.ierr.EQ.728) then
          call quit_program_error('ALLOCATION ERROR: OUT OF MEMORY',1,nocolor)
       elseif(ierr.EQ.151.OR.ierr.EQ.582) then
          call quit_program_error('ALLOCATION ERROR: ALREADY ALLOCATED',1,nocolor)
       else
!         STOP 'DURING ALLOCATION, FORTRAN ERROR ',ierr,' OCCURRED'
          call quit_program_error('DURING ALLOCATION, FORTRAN ERROR OCCURRED',1,nocolor)
       endif iallc
    endif !! (ierr.NE.0) then                                                                            
    ncounta = ncounta+1
    call mymemunits(dble(memory_used),iunits,unitchgdec,unitchgbin,nocolor)
    call print_mem_text(' * ProcAlloc invoked by '//trim(array_name)//' with type (' &
                       &//trim(allchartype(iarrtype_def))//')'//'. Total calls: ' &
                       &//int2str(ncounta)//' times. Use: '//dbl2str(memory_used/unitchgdec,2)&
                       &//unitsdec(iunits),verbose,nocolor)
    if (lsave_def) then
       call print_mem_text('    (Saving array data ...)',verbose,nocolor)
       nsize = 1
       if (allocated(iallocmem)) then
          nsize = size(iallocmem(:))+1
          if (allocated(temptxt)) call quit_program_error('temp array allocated unexpectly',1,nocolor)
          if (allocated(itmparr)) call quit_program_error('int temp array allocated unexpectly',1,nocolor)
          if (allocated(ltemp)  ) call quit_program_error('log temp array allocated unexpectly',1,nocolor)
          !! Expanding allocname:
             call move_alloc(allocname,temptxt)
             allocate(allocname(nsize))
             allocname(:nsize-1)(:) = temptxt(:)(:)
             deallocate(temptxt)
          !! Expanding lallocsaved:
             call move_alloc(lallocsaved,ltemp)
             allocate(lallocsaved(nsize))
             lallocsaved(:nsize-1) = ltemp(:)
             deallocate(ltemp)
          !! Expanding iallocmem:
             call move_alloc(iallocmem,itmparr)
             allocate(iallocmem(nsize))
             iallocmem(:nsize-1)    = itmparr(:)
             deallocate(itmparr)
          !! Expanding ialloctype:
             call move_alloc(ialloctype,itmparr)
             allocate(ialloctype(nsize))
             ialloctype(:nsize-1)    = itmparr(:)
             deallocate(itmparr)
          !! Expanding  allocmem:
             call move_alloc( allocmem,temptxt)
             allocate( allocmem(nsize))
              allocmem(:nsize-1)(:) = temptxt(:)
             deallocate(temptxt)
          !! Expanding   allocwhere:
             call move_alloc(allocwhere,temptxt)
             allocate(allocwhere(nsize))
             allocwhere(:nsize-1)(:) = temptxt(:)(:)
             deallocate(temptxt)
          !! Expanding deallocwhere:
             call move_alloc(deallocwhere,temptxt)
             allocate(deallocwhere(nsize))
             deallocwhere(:nsize-1)(:) = temptxt(:)(:)
             deallocate(temptxt)
       else
          if(allocated(  allocname) ) call quit_program_error('array of names allocated unexpectly'     ,1,nocolor)
          if(allocated( ialloctype) ) call quit_program_error('array of types allocated unexpectly'     ,1,nocolor)
          if(allocated(deallocwhere)) call quit_program_error('dealloc state array allocated unexpectly',1,nocolor)
          if(allocated(  allocwhere)) call quit_program_error('alloc state array allocated unexpectly'  ,1,nocolor)
          if(allocated( lallocsaved)) call quit_program_error('state array allocated unexpectly'        ,1,nocolor)
          allocate(iallocmem(nsize),allocmem(nsize),allocname(nsize),deallocwhere(nsize),allocwhere(nsize) &
                  ,ialloctype(nsize),lallocsaved(nsize))
       endif !! (allocated(alloctxt)) then
       call mymemunits(dble(isize),iunits,unitchgdec,unitchgbin,nocolor)
       allocname(nsize)(:)    = '+'//trim(arrayname_def(:))//'+'
       ialloctype(nsize)      = iarrtype_def
       iallocmem(nsize)       = isize
       allocmem(nsize)        = dbl2str(dble(isize)/unitchgbin,2)//unitsbin(iunits)
       allocwhere(nsize)(:)   = trim(routine_name(:))
       deallocwhere(nsize)(:) = 'NOT YET'
       lallocsaved(nsize)     = .TRUE.
       if (verbose.AND.allocated(iallocmem)) then
          call   mem_report(nocolor)
       else if (verbose) then
           call quit_program_error('Something is wrong while saving array data',1,nocolor)
       endif !! (verbose.AND.allocated(alloctxt)) then
    else
       call print_mem_text('    (Array info is kept unsaved ...)',verbose,nocolor)
    endif !! (lsave_def) then
  end subroutine ProcAlloc

  subroutine     ProcDEAlloc(isize,nocolor,verbose,routine_name,array_name,lsearch,perform_sub)

    use SUFR_system        ,only                 : quit_program_error
    use SUFR_text          ,only                 : int2str,dbl2str
    integer   ( kind =  4 )                     :: i,nsize,iunits
    integer   ( kind =  4 ),intent(in)          :: isize
    logical                ,intent(in)          :: nocolor,verbose
    logical                ,intent(in),optional :: lsearch,perform_sub
    logical                                     :: lsearchdef,perform_sub_def
    real      (kind=8)                          :: unitchgdec,unitchgbin
    character              ,intent(in),optional :: routine_name*(*),array_name*(*)
    character (  len = 30 )                     :: rouname_def,arrayname_def

    lsearchdef      = .FALSE.
    perform_sub_def = .FALSE.
    if (present(routine_name)) then
       rouname_def   = routine_name(:)
    else
       rouname_def   = 'UNKNOWN ROUTINE'
    endif !! (present(routine_name)) then
    if (present(array_name)  ) then
       arrayname_def =   array_name(:)
    else
       arrayname_def = 'UNKNOWN ARRAY'
    endif !! (present(array_name)  ) then
    if (present(lsearch)    ) lsearchdef      = lsearch
    if (present(perform_sub)) perform_sub_def = perform_sub
    if (perform_sub_def) memory_used = memory_used-isize
    nsize = size(iallocmem(:))
    max_mem = max(max_mem,memory_used)
    ncountd = ncountd+1
    call mymemunits(dble(memory_used),iunits,unitchgdec,unitchgbin,nocolor)
    call print_mem_text(' * ProcDEAlloc invoked by '//trim(array_name)//'. Total calls: '&
         &//int2str(ncountd)//' times. Use: '//dbl2str(memory_used/unitchgdec,2)&
         &//unitsdec(iunits),verbose,nocolor)
    if (lsearchdef) then
       i = 1
       do while (.NOT.index('+'//trim(array_name)//'+',trim(allocname(i))).GT.0)
          i = i+1
          if (i.GT.nsize) call quit_program_error('Requested array not found',1,nocolor)
       enddo !! while (if (index(trim(array_name),trim(allocname(i))).GT.0) then
       deallocwhere(i)(:) = ''
       deallocwhere(i)(:) = routine_name(:)
       lallocsaved(i)     = .FALSE.
       if (verbose.AND.allocated(iallocmem)) call mem_report(nocolor)
    endif !! (lsearchdef) then

  end subroutine ProcDEAlloc

  subroutine     print_mem_text(text,verbose,nocolor,iprint)

     use SUFR_system        ,only  : quit_program_error

     implicit none
     character(len=* ),intent(in) :: text
     character(len= 2)            :: bicode,color
     character(len=11)            :: colcode
     integer,optional,intent(in)  :: iprint
     integer                      :: iprintdef
     logical         ,intent(in)  :: nocolor,verbose

     iprintdef = 0
     if (present(iprint)) then
        if (iprint.NE.0.OR.iprint.NE.5.OR.iprint.NE.6) iprintdef = iprint
     endif !! (present(iprint)) then
     bicode(:) = ''
     bicode    = 'ME'
     color     = '35'
     if (verbose) then
        if (.NOT.nocolor) then
           colcode = achar(27)//'['//trim(color)//'m'//trim(bicode)//achar(27)//'[0m'
           write(*,'(3X,"(",A ,")",1X,A)') trim(colcode),trim(text)
        else
           write(*,'(3X,"(",A2,")",1X,A)') trim(bicode ),trim(text)
        endif !! (.NOT.nocolor) then
     endif !! (verbose) then
     if (iprintdef.NE.0) write(iprintdef,'(3X,"(",A2,")",1X,A)') trim(bicode),trim(text)

  end subroutine     print_mem_text

  subroutine     mymemunits(val    ,iunits,unitchgdec,unitchgbin,nocolor)
     use SUFR_system         ,only: warn
     implicit none
     integer(kind=     4),intent(out) :: iunits
     real   (kind=double),intent(out) :: unitchgdec,unitchgbin
     real   (kind=double),intent(in ) :: val
     logical             ,intent(in ) :: nocolor
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

  subroutine     mem_report(nocolor)

     use SUFR_text  ,only                     : int2str,dbl2str
     implicit none
     integer  (kind=4)                       :: i,iunits
     real     (kind=8)                       :: unitchgdec,unitchgbin
     character(len=100)                      :: label,label2
     logical           ,intent(in ),optional :: nocolor
     logical                                 :: nocolor_def

     nocolor_def = .FALSE.
     if (present(nocolor)) nocolor_def = nocolor
     write(*,*)
     call print_mem_text(' >>> CURRENT MEMORY USE SUMMARY <<<',.TRUE.,nocolor=nocolor)
     write(*,'(8X,A)') "(Table contains saved arrays only)"
     i         = size(iallocmem)
     maxn      = max(maxn     ,len_trim(int2str(i)))
     maxname   = max(maxname  ,len_trim(allocname(i)(:))-2)
     maxsize   = max(maxsize  ,len_trim(allocmem(i)))
     maxtype   = max(maxtype  ,len_trim(allchartype(ialloctype(i))))
     maxallwh  = max(maxallwh ,len_trim(allocwhere(i)))
     mxalldewh = max(mxalldewh,len_trim(deallocwhere(i)),10)
     i         = maxn+maxname+maxsize+maxtype+maxallwh+mxalldewh+6
     label     = '(3X,A'//int2str(maxn     )//',1X,A'//&
                         &int2str(maxname  )//',1X,A'//&
                         &int2str(maxsize  )//',1X,A'//&
                         &int2str(maxtype  )//',1X,A'//&
                         &int2str(maxallwh )//',2X,A'//&
                         &int2str(mxalldewh)//')'
     label2    = '(3X,'//int2str(i)//'("-"))'
     write(*,label2)
     write(*,label) '#','NAME','SIZE','TYPE','ALLOC in','DEALLOC in'
     write(*,label2)
     do i = 1,size(iallocmem)
        write(*,label) int2str(i)                                  &
                     , allocname(i)(2:len_trim(allocname(i)(:))-1) &
                     , trim(allocmem(i))                           &
                     , trim(allchartype(ialloctype(i)))            &
                     , trim(allocwhere(i)),trim(deallocwhere(i))
     enddo !! i = 1,size(iallocmem)
     write(*,label2)
     label2(:) = ''
     i         = maxn+maxname-13
     if (i.LT.1) i = 1
     label2    = '(3X,A,'//int2str(i)//'X,A)'
     call mymemunits(dble(memory_used)  ,iunits,unitchgdec,unitchgbin,nocolor_def)
     write(*,label2) "TOTAL          ",dbl2str(memory_used/unitchgbin,2)//unitsbin(iunits)
     i = sum(iallocmem,lallocsaved)
     call mymemunits(dble(i)            ,iunits,unitchgdec,unitchgbin,nocolor_def)
     write(*,label2) "        (saved)",dbl2str(i/unitchgbin,2)//unitsbin(iunits)
     call mymemunits(dble(memory_used-i),iunits,unitchgdec,unitchgbin,nocolor_def)
     write(*,label2) "      (unsaved)",dbl2str((memory_used-i)/unitchgbin,2)//unitsbin(iunits)
     call mymemunits(dble(max_mem)      ,iunits,unitchgdec,unitchgbin,nocolor_def)
     write(*,label2) "          (max)",dbl2str(max_mem/unitchgbin,2)//unitsbin(iunits)
     write(*,*)

  end subroutine mem_report

end module memory_use
