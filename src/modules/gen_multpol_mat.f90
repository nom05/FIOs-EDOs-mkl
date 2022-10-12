module     gen_multpol_mat
!----------------------------------------------------------------------------------------!
! Module written from the source code of Multiwfn (MIT-based license)                    !
!                                                                                        !
! LICENSE INFORMATION from the author:                                                   !
! To download and use Multiwfn, you are required to read and agree the following terms:  !
! (a) Currently Multiwfn is free of charge and open-source for both academic and         !
!     commercial usages, anyone is allowed to freely distribute the original or their    !
!     modified Multiwfn codes to others.                                                 !
! (b) Multiwfn can be distributed as a free component of commercial code. Selling        !
!     modified version of Multiwfn may also be granted, however, obtaining prior consent !
!     from the original author of Multiwfn (Tian Lu) is needed.                          !
! (c) If Multiwfn is utilized in your work, or your own code incorporated any part of    !
!     Multiwfn code, at least the original paper of Multiwfn MUST BE cited in main text  !
!     of your work or code:                                                              !
!                                                                                        !
!                Tian Lu, Feiwu Chen, J. Comput. Chem., 33, 580-592 (2012).              !
!                                                                                        !
! (d) There is no warranty of correctness of the results produced by Multiwfn, the author!
!     of Multiwfn does not hold responsibility in any way for any consequences arising   !
!     from the use of the Multiwfn.                                                      !
!----------------------------------------------------------------------------------------!
        
  use SUFR_kinds     ,only: double,dbl
  use SUFR_system    ,only: quit_program_error
  use commonmod      ,only: maxnbf,nfb,ncenter=>nato,x,y,z,nshell=>nsh,itsh &
                          , shell2atom=>icsh,shellcon=>inpps,primexp=>primexps,concoeff=>concoeffs &
                          , SPconcoeff=>pspconcoeffs,print_message,debuglabel,nthreads=>nproc &
                          , ilendbglbl,lconduc,dip,over
          
  implicit none

  type primtype
    integer              :: center,type !# nuclei that the basis function centered on and its function type
    real(double)         :: exp         !Exponent
  end type

  real(double),parameter :: pi=3.141592653589793_dbl
  integer                :: isphergau=0 !By default, all basis functions are cartesian type, =1 means spherical 
!Here s,p,d sequences are identical to .wfn, .wfx, .fch, .molden  !Note: Sequence in .fch = sequence in Gaussian
!Here f sequence is identical to .wfn, .wfx, but not identical to .fch and .molden
!Here g sequence is identical to .fch, .wfn does not support higher than f function, not identical to .wfx and .molden
!here h sequence is identical to .wfx and .fch, .molden doesn't support h
!Notice: The .wfn produced by G09 B.01 and later supports g and h, the definition is identical to here, and
!thus can be normally loaded
!Overall, spd: Multiwfn=wfn=wfx=fch=molden   f: Multiwfn=wfn=wfx!=fch=molden   
!           g: Multiwfn=fch!=wfx=molden=Molden2AIM   h: Multiwfn=wfx=fch
  integer                :: type2ix(56)=(/ 0, 1,0,0, 2,0,0,1,1,0, 3,0,0,2,2,0,1,1,0,1, &
                                           0,0,0,0,0,1,1,1,1,2,2,2,3,3,4, &
                                           0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5 /)
  integer                :: type2iy(56)=(/ 0, 0,1,0, 0,2,0,1,0,1, 0,3,0,1,0,2,2,0,1,1, &
                                           0,1,2,3,4,0,1,2,3,0,1,2,0,1,0, &
                                           0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0 /)
  integer                :: type2iz(56)=(/ 0, 0,0,1, 0,0,2,0,1,1, 0,0,3,0,1,1,0,2,2,1, &
                                           4,3,2,1,0,3,2,1,0,2,1,0,1,0,0, &
                                           5,4,3,2,1,0,4,3,2,1,0,3,2,1,0,2,1,0,1,0,0 /)
!Convert shell type to the number of basis functions in the shell: 0=s,1=p,-1=sp,2=6d,-2=5d,3=10f,-3=7f,4=15g,-4=9g,5=21h,-5=11h
  integer                :: shtype2nbas(-5:5)=(/ 11,9,7,5,4,1,3,6,10,15,21 /) 

!-------- Variables for wfn information
  integer :: nprims=0 !Number of primitive functions
!-------- Variables for nuclei & GTF & Orbitals
  type(primtype),allocatable :: b(:)
!-------- Variables when basis functions are basis rather than primitive function as basis
! integer :: nbasis=0,nshell=0,nprimshell=0 !The number of basis, basis shell and primitive shell. SP shell is 
  integer :: nbasis=0
!counted as S and P shell separately
  integer,allocatable :: primstart(:),primend(:) !The ith element means the GTF from where to where is attributed
!to the ith Cartesian basis function (which may be yielded during reading spherical wavefunction)
  real(double),allocatable :: primconnorm(:) !element i means the contract. coeff. * normalization coeff. of GTF i,
! can be used for e.g. constructing basis integral from GTF integral
  real(double),allocatable :: Dbas(:,:,:) !Dipole moment integral matrix, the first index 1,2,3=X,Y,Z,
  real(double),allocatable :: Sbas(:,:) !Overlap matrix
!the last two indices are basis index

!!! Other external parameter !!!
! integer :: nthreads=4

!!------------------- Root and weight of hermite polynomial
  real(double) :: Rhm(10,10),Whm(10,10)
  data Rhm( 1, 1) / 0.0D0                      /
  data Rhm( 2, 1) /-0.70710678118654752440D+00 /
  data Rhm( 2, 2) / 0.70710678118654752440D+00 /
  data Rhm( 3, 1) /-1.22474487139158904910D+00 /
  data Rhm( 3, 2) / 0.0D0                      /
  data Rhm( 3, 3) / 1.22474487139158904910D+00 /
  data Rhm( 4, 1) /-1.65068012388578455588D+00 /
  data Rhm( 4, 2) /-0.52464762327529031788D+00 /
  data Rhm( 4, 3) / 0.52464762327529031788D+00 /
  data Rhm( 4, 4) / 1.65068012388578455588D+00 /
  data Rhm( 5, 1) /-2.02018287045608563293D+00 /
  data Rhm( 5, 2) /-0.95857246461381850711D+00 /
  data Rhm( 5, 3) / 0.0D0                      /
  data Rhm( 5, 4) / 0.95857246461381850711D+00 /
  data Rhm( 5, 5) / 2.02018287045608563293D+00 /
  data Rhm( 6, 1) /-2.35060497367449222283D+00 /
  data Rhm( 6, 2) /-1.33584907401369694971D+00 /
  data Rhm( 6, 3) /-0.43607741192761650868D+00 /
  data Rhm( 6, 4) / 0.43607741192761650868D+00 /
  data Rhm( 6, 5) / 1.33584907401369694971D+00 /
  data Rhm( 6, 6) / 2.35060497367449222283D+00 /
  data Rhm( 7, 1) /-2.65196135683523349245D+00 /
  data Rhm( 7, 2) /-1.67355162876747144503D+00 /
  data Rhm( 7, 3) /-0.81628788285896466304D+00 /
  data Rhm( 7, 4) / 0.0D0                      /
  data Rhm( 7, 5) / 0.81628788285896466304D+00 /
  data Rhm( 7, 6) / 1.67355162876747144503D+00 /
  data Rhm( 7, 7) / 2.65196135683523349245D+00 /
  data Rhm( 8, 1) /-2.93063742025724401922D+00 /
  data Rhm( 8, 2) /-1.98165675669584292585D+00 /
  data Rhm( 8, 3) /-1.15719371244678019472D+00 /
  data Rhm( 8, 4) /-0.38118699020732211685D+00 /
  data Rhm( 8, 5) / 0.38118699020732211685D+00 /
  data Rhm( 8, 6) / 1.15719371244678019472D+00 /
  data Rhm( 8, 7) / 1.98165675669584292585D+00 /
  data Rhm( 8, 8) / 2.93063742025724401922D+00 /
  data Rhm( 9, 1) /-3.19099320178152760723D+00 /
  data Rhm( 9, 2) /-2.26658058453184311180D+00 /
  data Rhm( 9, 3) /-1.46855328921666793167D+00 /
  data Rhm( 9, 4) /-0.72355101875283757332D+00 /
  data Rhm( 9, 5) / 0.0D0                      /
  data Rhm( 9, 6) / 0.72355101875283757332D+00 /
  data Rhm( 9, 7) / 1.46855328921666793167D+00 /
  data Rhm( 9, 8) / 2.26658058453184311180D+00 /
  data Rhm( 9, 9) / 3.19099320178152760723D+00 /
  data Rhm(10, 1) /-3.43615911883773760333D+00 /
  data Rhm(10, 2) /-2.53273167423278979641D+00 /
  data Rhm(10, 3) /-1.75668364929988177345D+00 /
  data Rhm(10, 4) /-1.03661082978951365418D+00 /
  data Rhm(10, 5) /-0.34290132722370460879D+00 /
  data Rhm(10, 6) / 0.34290132722370460879D+00 /
  data Rhm(10, 7) / 1.03661082978951365418D+00 /
  data Rhm(10, 8) / 1.75668364929988177345D+00 /
  data Rhm(10, 9) / 2.53273167423278979641D+00 /
  data Rhm(10,10) / 3.43615911883773760333D+00 /

  data Whm( 1, 1) / 1.77245385090551602730D+00 / ! SQRT(PI)
  data Whm( 2, 1) / 8.86226925452758013649D-01 /
  data Whm( 2, 2) / 8.86226925452758013649D-01 /
  data Whm( 3, 1) / 2.95408975150919337883D-01 /
  data Whm( 3, 2) / 1.18163590060367735153D+00 /
  data Whm( 3, 3) / 2.95408975150919337883D-01 /
  data Whm( 4, 1) / 8.13128354472451771430D-02 /
  data Whm( 4, 2) / 8.04914090005512836506D-01 /
  data Whm( 4, 3) / 8.04914090005512836506D-01 /
  data Whm( 4, 4) / 8.13128354472451771430D-02 /
  data Whm( 5, 1) / 1.99532420590459132077D-02 /
  data Whm( 5, 2) / 3.93619323152241159828D-01 /
  data Whm( 5, 3) / 9.45308720482941881226D-01 /
  data Whm( 5, 4) / 3.93619323152241159828D-01 /
  data Whm( 5, 5) / 1.99532420590459132077D-02 /
  data Whm( 6, 1) / 4.53000990550884564086D-03 /
  data Whm( 6, 2) / 1.57067320322856643916D-01 /
  data Whm( 6, 3) / 7.24629595224392524092D-01 /
  data Whm( 6, 4) / 7.24629595224392524092D-01 /
  data Whm( 6, 5) / 1.57067320322856643916D-01 /
  data Whm( 6, 6) / 4.53000990550884564086D-03 /
  data Whm( 7, 1) / 9.71781245099519154149D-04 /
  data Whm( 7, 2) / 5.45155828191270305922D-02 /
  data Whm( 7, 3) / 4.25607252610127800520D-01 /
  data Whm( 7, 4) / 8.10264617556807326765D-01 /
  data Whm( 7, 5) / 4.25607252610127800520D-01 /
  data Whm( 7, 6) / 5.45155828191270305922D-02 /
  data Whm( 7, 7) / 9.71781245099519154149D-04 /
  data Whm( 8, 1) / 1.99604072211367619206D-04 /
  data Whm( 8, 2) / 1.70779830074134754562D-02 /
  data Whm( 8, 3) / 2.07802325814891879543D-01 /
  data Whm( 8, 4) / 6.61147012558241291030D-01 /
  data Whm( 8, 5) / 6.61147012558241291030D-01 /
  data Whm( 8, 6) / 2.07802325814891879543D-01 /
  data Whm( 8, 7) / 1.70779830074134754562D-02 /
  data Whm( 8, 8) / 1.99604072211367619206D-04 /
  data Whm( 9, 1) / 3.96069772632643819046D-05 /
  data Whm( 9, 2) / 4.94362427553694721722D-03 /
  data Whm( 9, 3) / 8.84745273943765732880D-02 /
  data Whm( 9, 4) / 4.32651559002555750200D-01 /
  data Whm( 9, 5) / 7.20235215606050957124D-01 /
  data Whm( 9, 6) / 4.32651559002555750200D-01 /
  data Whm( 9, 7) / 8.84745273943765732880D-02 /
  data Whm( 9, 8) / 4.94362427553694721722D-03 /
  data Whm( 9, 9) / 3.96069772632643819046D-05 /
  data Whm(10, 1) / 7.64043285523262062916D-06 /
  data Whm(10, 2) / 1.34364574678123269220D-03 /
  data Whm(10, 3) / 3.38743944554810631362D-02 /
  data Whm(10, 4) / 2.40138611082314686417D-01 /
  data Whm(10, 5) / 6.10862633735325798784D-01 /
  data Whm(10, 6) / 6.10862633735325798784D-01 /
  data Whm(10, 7) / 2.40138611082314686417D-01 /
  data Whm(10, 8) / 3.38743944554810631362D-02 /
  data Whm(10, 9) / 1.34364574678123269220D-03 /
  data Whm(10,10) / 7.64043285523262062916D-06 /

  contains

subroutine     proc_dip(dont_calc,debug,prnt_big_mats,lconduc,nocolor,mem_info)

  use SUFR_text      ,only: uppercase
  use memory_use
  use commonmod      ,only: partial_time,ntime,iout,textxyz
  use writemod       ,only: print_sym_mat
  implicit none

  integer                :: i
  character( len=    80) :: label
  logical,intent(in)     :: dont_calc,debug,prnt_big_mats,lconduc,mem_info,nocolor

  if (dont_calc) then
     if (.NOT.allocated(dip)) call quit_program_error('Multipole matrices array is unexpectly deallocated' &
                                                      ,1,nocolor)
!$omp parallel
  !$omp workshare
     dip = -dip
  !$omp end workshare
!$omp end parallel

     call mydealloc(itsh,'proc_dip','itsh',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:readfchk

  else
     call myalloc(dip,nfb,nfb,3,'proc_dip','dip',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:conductance_calc&FIOs
     if (lconduc) call myalloc(over,nfb,nfb,'proc_dip','over',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+D:conductance_calc
     call calcDbas(debug,mem_info,nocolor,lconduc)
     if (allocated(SPconcoeff)) &
        call mydealloc(SPconcoeff,'proc_dip','pspconcoeffs',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:readfchk
     call mydealloc(concoeff,'proc_dip','concoeffs',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:readfchk
     call mydealloc(primexp ,'proc_dip','primexps' ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:readfchk
  endif !! (read_mult_mat) then
  if (prnt_big_mats) then
     if (debug) PRINT *,debuglabel(:ilendbglbl),'Printing Multipole matrices...'
     if (debug.AND.lconduc) PRINT *,debuglabel(:ilendbglbl),'...and overlap matrix...'
     do i = 1,3
        label(:) = ''
        label    = 'Multipole matrix (BF) array for '//uppercase(textxyz(i))//' component'
        call print_sym_mat(iout,nfb,dip(:,:,i),trim(label))
     enddo !! i = 1,3
     if (lconduc) then
        label(:) = ''
        label    = 'Overlap matrix (BF) array'
        call print_sym_mat(iout,nfb,over(:,:),trim(label))
        if (debug) PRINT *,debuglabel(:ilendbglbl),'Printed Overlap matrix...'
     endif !! (lconduc) then
  endif !! (prnt_big_mats) then
  call mydealloc(shell2atom,'proc_dip','icsh',nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:readfchk
  if (.NOT.dont_calc) call partial_time(ntime,'Mult. matrices computation')

end subroutine proc_dip

subroutine calcDbas(debug,mem_info,nocolor,lconduc)

  use SUFR_text        ,only: int2str
  use memory_use

  implicit none

  character(len=100)       :: label
! integer,allocatable      :: shelltype(:),shell2atom(:),shellcon(:) !Degree of shell contraction
  integer,allocatable      :: shelltype(:)
! real(double),allocatable :: primexp(:),concoeff(:),SPconcoeff(:)
  integer                  :: i,ibasis,nbasis5D,ipos5D,ipos6D,ish,ishtyp5D,ishtyp6D &
                            , numshorb5D,numshorb6D,k,iexp,j,l,idir
  integer                  :: s2f(-5:5,21)=0 !Give shell type & orbital index to get functype
  real(double)             :: tnormgau,temp
  real(double)             :: conv5d6d(6,5),conv7f10f(10,7),conv9g15g(15,9),conv11h21h(21,11)
  real(double)             :: conv5d6dtr(5,6),conv7f10ftr(7,10),conv9g15gtr(9,15),conv11h21htr(11,21)
  !For backing up spherical basis functions
  integer,allocatable      :: shelltype5D(:)
  real(double),allocatable :: Dbas5D(:,:,:),Sbas5D(:,:)
  logical,intent(in)       :: debug,mem_info,nocolor,lconduc

  nbasis = nfb
  i = size(itsh)*bit_size(i)/8
  call move_alloc(itsh,shelltype) !! The original itsh is deallocated from here
  call ProcDEAlloc(i,nocolor,mem_info,'calcDbas','itsh',lsearch=.TRUE.,perform_sub=.TRUE.) !+A:readfchk
  call ProcAlloc(i,0,nocolor,mem_info,.FALSE.,'calcDbas','shelltype',perform_sum=.TRUE.) !+D:here
  s2f(-5,1:11)=(/ -32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22 /)
  s2f(-4,1:9)=(/ -21,-20,-19,-18,-17,-16,-15,-14,-13 /)
  s2f(-3,1:7)=(/ -12,-11,-10,-9,-8,-7,-6 /)
  s2f(-2,1:5)=(/ -5,-4,-3,-2,-1 /)
  s2f(-1,1:4)=(/ 1,2,3,4 /)
  s2f(0,1)=1
  s2f(1,1:3)=(/ 2,3,4 /)
  s2f(2,1:6)=(/ 5,6,7,8,9,10 /)
  s2f(3,1:10)=(/ 11,12,13,17,14,15,18,19,16,20 /) !Note: The sequence of f functions in Multiwfn is not identical 
  !to .fch, so convert here. While spdgh are identical
  s2f(4,1:15)=(/ 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 /)
  s2f(5,1:21)=(/ 36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56 /)

  call gensphcartab(conv5d6d,conv7f10f,conv9g15g,conv11h21h)

  conv5d6dtr=transpose(conv5d6d)
  conv7f10ftr=transpose(conv7f10f)
  conv9g15gtr=transpose(conv9g15g)
  conv11h21htr=transpose(conv11h21h)

!Note that Multiwfn allows cartesian and spherical harmonic basis functions mixed together.
!    If any basis function is spherical harmonic type, then isphergau=1.
!Only the spherical harmonic ones will be treated specially
  isphergau = 0
  if (any(shelltype<=-2)) isphergau=1
  if (any(abs(shelltype)>maxnbf)) then
     label(:) = ''
     label    = 'GTFs with angular moment > '//int2str(maxnbf)//' are unsupported'
     call quit_program_error(trim(label),1,nocolor)
  end if


!Backup spherical basis information (some of them may be Cartesian ones) with 5D suffix (of course, may be actually 7f, 9g, 11h...),
!convert them to cartesian type temporarily, at final stage recover them back, namely get Sbas, Ptot... in spherical basis
  if (isphergau==1) then
     call myalloc(shelltype5D,nshell,'calcDbas','shelltype5D',nocolor=nocolor,verbose=mem_info) !+D:here
     shelltype5D=shelltype
     where (shelltype<=-2) shelltype=-shelltype !Convert to cartesian type
     nbasis5D=nbasis
     nbasis=0
     do i=1,nshell
        nbasis = nbasis+shtype2nbas(shelltype(i))
     end do
  end if

!Allocate space for arrays
  nprims=0
  do i=1,nshell
     nprims=nprims+shtype2nbas(shelltype(i))*shellcon(i)
  end do
  allocate(b(nprims))
  i = nprims*4+nprims*8
  call ProcAlloc(i,0,nocolor,mem_info,.FALSE.,'calcDbas','b',perform_sum=.TRUE.) !+D:here
  call myalloc(primstart  ,nbasis,'calcDbas','primstart'  ,nocolor=nocolor,verbose=mem_info) !+D:here
  call myalloc(primend    ,nbasis,'calcDbas','primend'    ,nocolor=nocolor,verbose=mem_info) !+D:here
  call myalloc(primconnorm,nprims,'calcDbas','primconnorm',nocolor=nocolor,verbose=mem_info) !+D:here
  !Fill Cobasa and CObasb
  if (isphergau==1) then
     !Map 5D coefficient to 6D coefficient
     ipos5D=1
     ipos6D=1
     do ish=1,nshell
        ishtyp5D   = shelltype5D(ish     )
        ishtyp6D   =   shelltype(ish     )
        numshorb5D = shtype2nbas(ishtyp5D)
        numshorb6D = shtype2nbas(ishtyp6D)
        ipos5D     = ipos5D + numshorb5D
        ipos6D     = ipos6D + numshorb6D
     end do
  end if

  call print_message('i',"Converting basis function information to GTF information")
!Fill: b,primstart,primend,primconnorm
  k      = 1 !Current position of GTF
  iexp   = 1
  ibasis = 1 !Current position of basis
!Note: Below commented with !!! means the line associated to setting basis information
  do i=1,nshell !cycle each shell
     b(k:k+shellcon(i)*shtype2nbas(shelltype(i))-1)%center=shell2atom(i)
     do j=1,shtype2nbas(shelltype(i)) !cycle each basis(orbital) in each shell
        b(k:k+shellcon(i)-1)%type=s2f(shelltype(i),j)
        primstart(ibasis)=k !!! From where the GTFs attributed to ibasis'th basis
        primend(ibasis)=k+shellcon(i)-1 !!! To where the GTFs attributed to ibasis'th basis
        !write(*,*) i,j,ibasis,primstart(ibasis),primend(ibasis)
        do l=1,shellcon(i) !cycle each GTF in each basis in each shell
           b(k)%exp=primexp(iexp+l-1)
           tnormgau=normgau(b(k)%type,b(k)%exp)  !!!Normalization coefficient of GTFs
           temp=concoeff(iexp+l-1)  !!!Contraction coefficient of GTFs
           if (shelltype(i)==-1.and.j/=1) temp=SPconcoeff(iexp+l-1)
           primconnorm(k)=temp*tnormgau !Combines contraction and normalization coefficient
           k=k+1
        end do
        ibasis=ibasis+1
     end do
     iexp=iexp+shellcon(i)
  end do

  if (lconduc) then
     call print_message('i',"Generating overlap and multipole matrices")
     call myalloc(Sbas,nbasis,nbasis,'calcDbas','Sbas',nocolor=nocolor,verbose=mem_info) !+D:here
     call genSbas
  else
     call print_message('i',"Generating multipole matrices")
  endif !! (lconduc) then
  call myalloc(Dbas,nbasis,nbasis,3,'calcDbas','Dbas',nocolor=nocolor,verbose=mem_info) !+D:here
  call genDbas

  deallocate(b)
  i = nprims*4+nprims*8
  call ProcDEAlloc(i,nocolor,mem_info,'calcDbas','b',perform_sub=.TRUE.)                !+A:here
  call mydealloc(primconnorm,'calcDbas','primconnorm',nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(primend    ,'calcDbas','primend'    ,nocolor=nocolor,verbose=mem_info) !+A:here
  call mydealloc(primstart  ,'calcDbas','primstart'  ,nocolor=nocolor,verbose=mem_info) !+A:here

  call mydealloc(shellcon   ,'calcDbas','inpps'      ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:readfchk
  call mydealloc(z          ,'calcDbas','z'          ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:readfchk
  call mydealloc(y          ,'calcDbas','y'          ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:readfchk
  call mydealloc(x          ,'calcDbas','x'          ,nocolor=nocolor,verbose=mem_info,lsave=.TRUE.) !+A:readfchk

!Convert information from Cartesian basis to spherical basis
  if (isphergau==1) then
     call print_message('i',"Back converting BF information from Cartesian to spher. type")
     !Map cartesian overlap matrix to spherical
     call myalloc(Dbas5D,nbasis5D,nbasis5D,3,'calcDbas','Dbas5D',nocolor=nocolor,verbose=mem_info) !+D:here
     if (lconduc) call myalloc(Sbas5D,nbasis5D,nbasis5D,'calcDbas','Sbas5D',nocolor=nocolor,verbose=mem_info) !+A:here
     ipos5D=1
     ipos6D=1
     do ish=1,nshell
        ishtyp5D=shelltype5D(ish)
        ishtyp6D=shelltype(ish)
        numshorb5D=shtype2nbas(ishtyp5D)
        numshorb6D=shtype2nbas(ishtyp6D)
        !Now contract columns of Dbas
        if (ishtyp5D>=-1) Sbas(:,ipos5D:ipos5D+numshorb5D-1) = &
                          Sbas(:,ipos6D:ipos6D+numshorb6D-1) !S, P, SP or other Cartesian shells
        if (ishtyp5D==-2) Sbas(:,ipos5D:ipos5D+numshorb5D-1) = &
                   matmul(Sbas(:,ipos6D:ipos6D+numshorb6D-1),conv5d6d  )
        if (ishtyp5D==-3) Sbas(:,ipos5D:ipos5D+numshorb5D-1) = &
                   matmul(Sbas(:,ipos6D:ipos6D+numshorb6D-1),conv7f10f )
        if (ishtyp5D==-4) Sbas(:,ipos5D:ipos5D+numshorb5D-1) = &
                   matmul(Sbas(:,ipos6D:ipos6D+numshorb6D-1),conv9g15g )
        if (ishtyp5D==-5) Sbas(:,ipos5D:ipos5D+numshorb5D-1) = &
                   matmul(Sbas(:,ipos6D:ipos6D+numshorb6D-1),conv11h21h)
        !Now contract rows of Dbas
        if (ishtyp5D>=-1)              Sbas(ipos5D:ipos5D+numshorb5D-1,:) = &
                                       Sbas(ipos6D:ipos6D+numshorb6D-1,:) !S, P, SP or other Cart. shells
        if (ishtyp5D==-2)              Sbas(ipos5D:ipos5D+numshorb5D-1,:) = &
                   matmul(conv5d6dtr  ,Sbas(ipos6D:ipos6D+numshorb6D-1,:))
        if (ishtyp5D==-3)              Sbas(ipos5D:ipos5D+numshorb5D-1,:) = &
                   matmul(conv7f10ftr ,Sbas(ipos6D:ipos6D+numshorb6D-1,:))
        if (ishtyp5D==-4)              Sbas(ipos5D:ipos5D+numshorb5D-1,:) = &
                   matmul(conv9g15gtr ,Sbas(ipos6D:ipos6D+numshorb6D-1,:))
        if (ishtyp5D==-5)              Sbas(ipos5D:ipos5D+numshorb5D-1,:) = &
                   matmul(conv11h21htr,Sbas(ipos6D:ipos6D+numshorb6D-1,:))
        do idir=1,3
           !Now contract columns of Dbas
           if (ishtyp5D>=-1) Dbas(:,ipos5D:ipos5D+numshorb5D-1,idir) = &
                             Dbas(:,ipos6D:ipos6D+numshorb6D-1,idir) !S, P, SP or other Cartesian shells
           if (ishtyp5D==-2) Dbas(:,ipos5D:ipos5D+numshorb5D-1,idir) = &
                      matmul(Dbas(:,ipos6D:ipos6D+numshorb6D-1,idir),conv5d6d  )
           if (ishtyp5D==-3) Dbas(:,ipos5D:ipos5D+numshorb5D-1,idir) = &
                      matmul(Dbas(:,ipos6D:ipos6D+numshorb6D-1,idir),conv7f10f )
           if (ishtyp5D==-4) Dbas(:,ipos5D:ipos5D+numshorb5D-1,idir) = &
                      matmul(Dbas(:,ipos6D:ipos6D+numshorb6D-1,idir),conv9g15g )
           if (ishtyp5D==-5) Dbas(:,ipos5D:ipos5D+numshorb5D-1,idir) = &
                      matmul(Dbas(:,ipos6D:ipos6D+numshorb6D-1,idir),conv11h21h)
           !Now contract rows of Dbas
           if (ishtyp5D>=-1)              Dbas(ipos5D:ipos5D+numshorb5D-1,:,idir) = &
                                          Dbas(ipos6D:ipos6D+numshorb6D-1,:,idir) !S, P, SP or other Cart. shells
           if (ishtyp5D==-2)              Dbas(ipos5D:ipos5D+numshorb5D-1,:,idir) = &
                      matmul(conv5d6dtr  ,Dbas(ipos6D:ipos6D+numshorb6D-1,:,idir))
           if (ishtyp5D==-3)              Dbas(ipos5D:ipos5D+numshorb5D-1,:,idir) = &
                      matmul(conv7f10ftr ,Dbas(ipos6D:ipos6D+numshorb6D-1,:,idir))
           if (ishtyp5D==-4)              Dbas(ipos5D:ipos5D+numshorb5D-1,:,idir) = &
                      matmul(conv9g15gtr ,Dbas(ipos6D:ipos6D+numshorb6D-1,:,idir))
           if (ishtyp5D==-5)              Dbas(ipos5D:ipos5D+numshorb5D-1,:,idir) = &
                      matmul(conv11h21htr,Dbas(ipos6D:ipos6D+numshorb6D-1,:,idir))
        end do
        ipos5D=ipos5D+numshorb5D
        ipos6D=ipos6D+numshorb6D
     end do
     Dbas5D=Dbas(1:nbasis5D,1:nbasis5D,:)
     if (lconduc) Sbas5D=Sbas(1:nbasis5D,1:nbasis5D)

     !Recover spherical Gaussian basis function information
     nbasis    = nbasis5D
     shelltype = shelltype5D
     call mydealloc(shelltype5D,'calcDbas','shelltype5D',nocolor=nocolor,verbose=mem_info) !+A:here
     ibasis    = 1
     do i=1,nshell
        do j=1,shtype2nbas(shelltype(i))
           ibasis=ibasis+1
        end do
     end do
     call mydealloc(Dbas,'calcDbas','Dbas',nocolor=nocolor,verbose=mem_info) !+A:here
     i = size(Dbas5D)*bit_size(i)*size_double_over_size_real/8
     call ProcDEAlloc(i,nocolor,mem_info,'calcDbas','Dbas5D',perform_sub=.TRUE.) !+A:here
     call move_alloc(Dbas5D,dip)
     if (lconduc) then
        call mydealloc(Sbas,'calcDbas','Sbas',nocolor=nocolor,verbose=mem_info) !+A:here
        i = size(Sbas5D)*bit_size(i)*size_double_over_size_real/8
        call ProcDEAlloc(i,nocolor,mem_info,'calcDbas','Sbas5D',perform_sub=.TRUE.) !+A:here
        call move_alloc(Sbas5D,over)
     endif !! (lconduc) then
  else
     i = size(Dbas)*bit_size(i)*size_double_over_size_real/8
     call ProcDEAlloc(i,nocolor,mem_info,'calcDbas','Dbas',perform_sub=.TRUE.) !+A:here
     call move_alloc(Dbas,dip)
     if (lconduc) then
        i = size(Sbas)*bit_size(i)*size_double_over_size_real/8
        call ProcDEAlloc(i,nocolor,mem_info,'calcDbas','Sbas',perform_sub=.TRUE.) !+A:here
        call move_alloc(Sbas,over)
     endif !! (lconduc) then
  end if
  call ProcDEAlloc(i,nocolor,mem_info,'calcDbas','shelltype',perform_sub=.TRUE.) !+A:readfchk

 if (debug) then
    do i = 1,3
       write(*,'(1X,A,I1,1X,E13.6,A,E13.6,/,8X,E13.6,A,E13.6)') debuglabel(:ilendbglbl),&
                                i,dip(1  ,1,i),'=dip(1,1,i) dip(1,n,i)=',dip(1  ,nfb,i) &
                                 ,dip(nfb,1,i),'=dip(n,1,i) dip(n,n,i)=',dip(nfb,nfb,i)
    enddo !! i = 1,3
    if (lconduc) write(*,'(1X,A,E13.6,A,E13.6,/,8X,E13.6,A,E13.6)') debuglabel(:ilendbglbl) &
                                         ,over(1  ,1),'=over(1,1) over(1,n)=',over(1  ,nfb) &
                                         ,over(nfb,1),'=over(n,1) over(n,n)=',over(nfb,nfb)
 endif !! (debug) then

end subroutine calcDbas


!!!-------- Evaluate dipole moment integral for two unnormalized GTFs, <GTF|-r|GTF>,
!the negative charge of electron has been considered!
!~p arguments are the shifts of GTF index as doSintactual
!xint/yint/zint correspond to dipole moment integral in X/Y/Z
subroutine dodipoleint(iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p,xint,yint,zint)

   implicit none

   real(double) :: xint,yint,zint,x1,y1,z1,x2,y2,z2,ee1,ee2,ep,sqrtep,px,py,pz &
                 , expterm,sx,sy,sz,sxx,syy,szz,tmp,term1,term2
   integer      :: iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p,ix1,iy1,iz1,ix2,iy2,iz2,numx,numy,numz,i

   x1=x(b(iGTF)%center)
   y1=y(b(iGTF)%center)
   z1=z(b(iGTF)%center)
   x2=x(b(jGTF)%center)
   y2=y(b(jGTF)%center)
   z2=z(b(jGTF)%center)
   ee1=b(iGTF)%exp
   ee2=b(jGTF)%exp
   ep=ee1+ee2
   sqrtep=dsqrt(ep)
   px=(ee1*x1+ee2*x2)/ep
   py=(ee1*y1+ee2*y2)/ep
   pz=(ee1*z1+ee2*z2)/ep
   expterm=dexp( -ee1*ee2*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)/ep )
   ix1=type2ix(b(iGTF)%type)+ix1p
   iy1=type2iy(b(iGTF)%type)+iy1p
   iz1=type2iz(b(iGTF)%type)+iz1p
   ix2=type2ix(b(jGTF)%type)+ix2p
   iy2=type2iy(b(jGTF)%type)+iy2p
   iz2=type2iz(b(jGTF)%type)+iz2p
!First, calculate sx,sy,sz as usual as doSint
   numx=ceiling( (ix1+ix2+1)/2D0 ) !Need to calculate n points
   sx=0.0D0
   do i=1,numx
      tmp=Rhm(numx,i)/sqrtep+px
      term1=(tmp-x1)**ix1
      term2=(tmp-x2)**ix2
      sx=sx+Whm(numx,i)*term1*term2
   end do
   sx=sx/sqrtep
   numy=ceiling( (iy1+iy2+1)/2D0 )
   sy=0.0D0
   do i=1,numy
      tmp=Rhm(numy,i)/sqrtep+py
      term1=(tmp-y1)**iy1
      term2=(tmp-y2)**iy2
      sy=sy+Whm(numy,i)*term1*term2
   end do
   sy=sy/sqrtep
   numz=ceiling( (iz1+iz2+1)/2D0 )
   sz=0.0D0
   do i=1,numz
      tmp=Rhm(numz,i)/sqrtep+pz
      term1=(tmp-z1)**iz1
      term2=(tmp-z2)**iz2
      sz=sz+Whm(numz,i)*term1*term2
   end do
   sz=sz/sqrtep
!Second, calculate overlap integral in X,Y,Z directions but with X,Y,Z coordinate variables
!(relative to the original point of the whole system) to produce sxx,syy,szz
   numx=ceiling( (ix1+ix2+2)/2D0 ) !Because X variable is introduced, ix1+ix2+2 is used instead of ix1+ix2+1
   sxx=0.0D0
   do i=1,numx
      tmp=Rhm(numx,i)/sqrtep+px
      term1=(tmp-x1)**ix1
      term2=(tmp-x2)**ix2
      sxx=sxx+Whm(numx,i)*term1*term2*tmp
   end do
   sxx=sxx/sqrtep
   numy=ceiling( (iy1+iy2+2)/2D0 )
   syy=0.0D0
   do i=1,numy
      tmp=Rhm(numy,i)/sqrtep+py
      term1=(tmp-y1)**iy1
      term2=(tmp-y2)**iy2
      syy=syy+Whm(numy,i)*term1*term2*tmp
   end do
   syy=syy/sqrtep
   numz=ceiling( (iz1+iz2+2)/2D0 )
   szz=0.0D0
   do i=1,numz
      tmp=Rhm(numz,i)/sqrtep+pz
      term1=(tmp-z1)**iz1
      term2=(tmp-z2)**iz2
      szz=szz+Whm(numz,i)*term1*term2*tmp
   end do
   szz=szz/sqrtep

   xint=-sxx*sy*sz*expterm
   yint=-sx*syy*sz*expterm
   zint=-sx*sy*szz*expterm
end subroutine

!!!----------- Generate dipole moment integral matrix between all Cartesian basis functions <{basis}|-r|{basis}>
!Dbas should be allocated first. The resultant matrix is for Cartesian basis functions, may be converted to spherical-harmonic later
subroutine genDbas
!  use commonmod ,only: debug,debuglabel,ilendbglbl
   implicit none
   integer :: i,j,ii,jj
   real(double) :: xdiptmp,ydiptmp,zdiptmp
   Dbas=0E0_dbl

!$OMP PARALLEL DO SHARED(Dbas) PRIVATE(i,ii,j,jj,xdiptmp,ydiptmp,zdiptmp) schedule(dynamic) NUM_THREADS(nthreads)
   do i=1,nbasis
      do j=i,nbasis
         do ii=primstart(i),primend(i)
            do jj=primstart(j),primend(j)
               call dodipoleint(ii,jj,0,0,0,0,0,0,xdiptmp,ydiptmp,zdiptmp)
               Dbas(i,j,1)=Dbas(i,j,1)+primconnorm(ii)*primconnorm(jj)*xdiptmp
               Dbas(i,j,2)=Dbas(i,j,2)+primconnorm(ii)*primconnorm(jj)*ydiptmp
               Dbas(i,j,3)=Dbas(i,j,3)+primconnorm(ii)*primconnorm(jj)*zdiptmp
            end do
         end do
         Dbas(j,i,:)=Dbas(i,j,:)
      end do
   end do
!$OMP END PARALLEL DO

end subroutine


!!----------- Generate spherical harmonic -> Cartesian basis function conversion table for d,f,g,h.
!The table comes from IJQC,54,83, which is used by Gaussian
!The sequence of d and f shell is also identical to .molden convention, but for g, another conversion table is used, &
!since in Multiwfn g cartesian shell starts from ZZZZ, but that of .molden starts from xxxx
subroutine gensphcartab(matd,matf,matg,math)
   real(8) :: matd(6,5),matf(10,7),matg(15,9),math(21,11)
   matd = 0E0_dbl
   matf = 0E0_dbl
   matg = 0E0_dbl
   math = 0E0_dbl
! From 5D: D 0,D+1,D-1,D+2,D-2
! To 6D:  1  2  3  4  5  6
!        XX,YY,ZZ,XY,XZ,YZ
!
! D0=-0.5*XX-0.5*YY+ZZ
   matd(1:3, 1) = (/ -0.5D0,-0.5D0,1D0 /)
! D+1=XZ
   matd(  5, 2) = 1D0
! D-1=YZ
   matd(  6, 3) = 1D0
! D+2=SQRT(3)/2*(XX-YY)
   matd(1:2, 4) = (/ sqrt(3D0)/2D0,-sqrt(3D0)/2D0 /)
! D-2=XY
   matd(  4, 5) = 1D0

! From 7F: F 0,F+1,F-1,F+2,F-2,F+3,F-3
! To 10F:  1   2   3   4   5   6   7   8   9  10      
!         XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ (Gaussian sequence, not identical to Multiwfn)
!
   matf(  3, 1) =  1D0
   matf(  6, 1) = -1.5D0/sqrt(5D0)
   matf(  9, 1) = -1.5D0/sqrt(5D0)
   matf(  1, 2) = -sqrt(3D0/8D0)
   matf(  4, 2) = -sqrt(3D0/40D0)
   matf(  7, 2) =  sqrt(6D0/5D0)
   matf(  2, 3) = -sqrt(3D0/8D0)
   matf(  5, 3) = -sqrt(3D0/40D0)
   matf(  8, 3) =  sqrt(6D0/5D0)
   matf(  6, 4) =  sqrt(3D0)/2D0
   matf(  9, 4) = -sqrt(3D0)/2D0
! F-2=XYZ
   matf( 10, 5) =  1D0
   matf(  1, 6) =  sqrt(5D0/8D0)
   matf(  4, 6) = -3D0/sqrt(8D0)
   matf(  2, 7) = -sqrt(5D0/8D0)
   matf(  5, 7) =  3D0/sqrt(8D0)

! From 9G: G 0,G+1,G-1,G+2,G-2,G+3,G-3,G+4,G-4
! To 15G:   1    2    3    4    5    6    7    8
!         ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ
!           9   10   11   12   13   14   15
!         XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
!
   matg(  1, 1) =  1D0
   matg(  3, 1) = -3D0*sqrt(3D0/35D0)
   matg(  5, 1) =  3D0/8D0
   matg( 10, 1) = -3D0*sqrt(3D0/35D0)
   matg( 12, 1) =  3D0/4D0*sqrt(3D0/35D0)
   matg( 15, 1) =  3D0/8D0
   matg(  6, 2) =  2D0*sqrt(5D0/14D0)
   matg(  8, 2) = -1.5D0/sqrt(14D0)
   matg( 13, 2) = -1.5D0*sqrt(5D0/14D0)
   matg(  2, 3) =  2D0*sqrt(5D0/14D0)
   matg(  4, 3) = -1.5D0*sqrt(5D0/14D0)
   matg( 11, 3) = -1.5D0/sqrt(14D0)
   matg(  3, 4) = -3D0*sqrt(3D0/28D0)
   matg(  5, 4) =  sqrt(5D0)/4D0
   matg( 10, 4) =  3D0*sqrt(3D0/28D0)
   matg( 15, 4) = -sqrt(5D0)/4D0
   matg(  7, 5) =  3D0/sqrt(7D0)
   matg(  9, 5) = -sqrt(5D0/28D0)
   matg( 14, 5) = -sqrt(5D0/28D0)
   matg(  8, 6) = -3D0/sqrt(8D0)
   matg( 13, 6) =  sqrt(5D0/8D0)
   matg(  4, 7) = -sqrt(5D0/8D0)
   matg( 11, 7) =  3D0/sqrt(8D0)
   matg(  5, 8) =  sqrt(35D0)/8D0
   matg( 12, 8) = -3D0/4D0*sqrt(3D0)
   matg( 15, 8) =  sqrt(35D0)/8D0
   matg(  9, 9) = -sqrt(5D0)/2D0
   matg( 14, 9) =  sqrt(5D0)/2D0

! From 11H: H 0,H+1,H-1,H+2,H-2,H+3,H-3,H+4,H-4,H+5,H-5
! To 21H:   1     2     3     4     5     6     7     8     9    10
!         ZZZZZ YZZZZ YYZZZ YYYZZ YYYYZ YYYYY XZZZZ XYZZZ XYYZZ XYYYZ 
!          11    12    13    14    15    16    17    18    19    20    21
!         XYYYY XXZZZ XXYZZ XXYYZ XXYYY XXXZZ XXXYZ XXXYY XXXXZ XXXXY XXXXX
!
   math(  1, 1) =  1D0
   math( 12, 1) = -5D0/sqrt(21D0)
   math(  3, 1) = -5D0/sqrt(21D0)
   math( 19, 1) =  5D0/8D0
   math(  5, 1) =  5D0/8D0
   math( 14, 1) =  sqrt(15D0/7D0)/4D0
   math(  7, 2) =  sqrt(5D0/3D0)
   math( 16, 2) = -3D0*sqrt(5D0/28D0)
   math(  9, 2) = -3D0/sqrt(28D0)
   math( 21, 2) =  sqrt(15D0)/8D0
   math( 11, 2) =  sqrt(5D0/3D0)/8D0
   math( 18, 2) =  sqrt(5D0/7D0)/4D0
   math(  2, 3) =  sqrt(5D0/3D0)
   math(  4, 3) = -3D0*sqrt(5D0/28D0)
   math( 13, 3) = -3D0/sqrt(28D0)
   math(  6, 3) =  sqrt(15D0)/8D0
   math( 20, 3) =  sqrt(5D0/3D0)/8D0
   math( 15, 3) =  sqrt(5D0/7D0)/4D0
   math( 12, 4) =  sqrt(5D0)/2D0
   math(  3, 4) = -sqrt(5D0)/2D0
   math( 19, 4) = -sqrt(35D0/3D0)/4D0
   math(  5, 4) =  sqrt(35D0/3D0)/4D0
   math(  8, 5) =  sqrt(5D0/3D0)
   math( 17, 5) = -sqrt(5D0/12D0)
   math( 10, 5) = -sqrt(5D0/12D0)
   math( 16, 6) =  sqrt(5D0/6D0)
   math(  9, 6) = -sqrt(1.5D0)
   math( 21, 6) = -sqrt(17.5D0)/8D0
   math( 11, 6) =  sqrt(17.5D0)/8D0
   math( 18, 6) =  sqrt(5D0/6D0)/4D0
   math(  4, 7) = -sqrt(5D0/6D0)
   math( 13, 7) =  sqrt(1.5D0)
   math( 20, 7) = -sqrt(17.5D0)/8D0
   math(  6, 7) =  sqrt(17.5D0)/8D0
   math( 15, 7) = -sqrt(5D0/6D0)/4D0
   math( 19, 8) =  sqrt(35D0)/8D0
   math(  5, 8) =  sqrt(35D0)/8D0
   math( 14, 8) = -0.75D0*sqrt(3D0)
   math( 17, 9) =  sqrt(5D0)/2D0
   math( 10, 9) = -sqrt(5D0)/2D0
   math( 21,10) =  3D0/8D0*sqrt(3.5D0)
   math( 11,10) =  5D0/8D0*sqrt(3.5D0)
   math( 18,10) = -1.25D0*sqrt(1.5D0)
   math(  6,11) =  3D0/8D0*sqrt(3.5D0)
   math( 20,11) =  5D0/8D0*sqrt(3.5D0)
   math( 15,11) = -1.25D0*sqrt(1.5D0)

end subroutine


!!---------------- Calculate factorial
integer function ft(i)
   integer i,j
   ft=i
   if (i==0) ft=1
   do j=i-1,1,-1
      ft=ft*j
   end do
end function

real(8) function normgau(itype,exp)

   implicit none
   integer :: ix,iy,iz,itype
   real(8) :: exp
   ix=type2ix(itype)
   iy=type2iy(itype)
   iz=type2iz(itype)
   normgau=(2*exp/pi)**0.75D0*dsqrt( (8*exp)**(ix+iy+iz)*ft(ix)*ft(iy)*ft(iz)/(ft(2*ix)*ft(2*iy)*ft(2*iz)) )

end function normgau

!!!------------------ Generate overlap matrix between all basis functions
!Sbas should be allocated first. The resultant matrix is for Cartesian basis functions, may be converted to spherical-harmonic later
subroutine genSbas
   implicit none
   integer :: i,j,ii,jj
   real(8) :: doSint

   Sbas = 0D0

!$OMP PARALLEL DO SHARED(Sbas) PRIVATE(i,ii,j,jj) schedule(dynamic) NUM_THREADS(nthreads)
   do i=1,nbasis
      do j=i,nbasis
         do ii=primstart(i),primend(i)
            do jj=primstart(j),primend(j)
               call doSintactual( ii , jj , 0  , 0  , 0  , 0  , 0  , 0  ,doSint)
               Sbas(i,j)=Sbas(i,j)+primconnorm(ii)*primconnorm(jj)*doSint
            end do
         end do
         Sbas(j,i)=Sbas(i,j)
      end do
   end do
!$OMP END PARALLEL DO

end subroutine

!!!------------------ Evaluate overlap integral for two unnormalized GTFs
!~p arguments are the shifts of GTF index
subroutine doSintactual(iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p,doSint)

   implicit none

   integer :: iGTF,jGTF,ix1p,iy1p,iz1p,ix2p,iy2p,iz2p,ix1,iy1,iz1,ix2,iy2,iz2,numx,numy,numz,i
   real(8) :: x1,y1,z1,x2,y2,z2,ee1,ee2,ep,sqrtep,px,py,pz,expterm,sx,sy,sz,tmp,term1,term2 &
            , doSint

   x1=x(b(iGTF)%center)
   y1=y(b(iGTF)%center)
   z1=z(b(iGTF)%center)
   x2=x(b(jGTF)%center)
   y2=y(b(jGTF)%center)
   z2=z(b(jGTF)%center)
   ee1=b(iGTF)%exp
   ee2=b(jGTF)%exp
   ep=ee1+ee2
   sqrtep=dsqrt(ep)
   px=(ee1*x1+ee2*x2)/ep
   py=(ee1*y1+ee2*y2)/ep
   pz=(ee1*z1+ee2*z2)/ep
   expterm=dexp( -ee1*ee2*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)/ep )
   ix1=type2ix(b(iGTF)%type)+ix1p
   iy1=type2iy(b(iGTF)%type)+iy1p
   iz1=type2iz(b(iGTF)%type)+iz1p
   ix2=type2ix(b(jGTF)%type)+ix2p
   iy2=type2iy(b(jGTF)%type)+iy2p
   iz2=type2iz(b(jGTF)%type)+iz2p
   !chen book,P103
   numx=ceiling( (ix1+ix2+1)/2D0 ) !Need to calculate n points
   sx=0.0D0
   do i=1,numx
      tmp=Rhm(numx,i)/sqrtep+px
      term1=(tmp-x1)**ix1
      term2=(tmp-x2)**ix2
      sx=sx+Whm(numx,i)*term1*term2
   end do
   sx=sx/sqrtep

   numy=ceiling( (iy1+iy2+1)/2D0 )
   sy=0.0D0
   do i=1,numy
      tmp=Rhm(numy,i)/sqrtep+py
      term1=(tmp-y1)**iy1
      term2=(tmp-y2)**iy2
      sy=sy+Whm(numy,i)*term1*term2
   end do
   sy=sy/sqrtep

   numz=ceiling( (iz1+iz2+1)/2D0 )
   sz=0.0D0
   do i=1,numz
      tmp=Rhm(numz,i)/sqrtep+pz
      term1=(tmp-z1)**iz1
      term2=(tmp-z2)**iz2
      sz=sz+Whm(numz,i)*term1*term2
   end do
   sz=sz/sqrtep

   doSint=sx*sy*sz*expterm

end subroutine doSintactual

end module gen_multpol_mat
