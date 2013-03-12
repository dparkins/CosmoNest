! posterior samples code
! place in source directory
! compile using
!ifc -W0 -WB -Vaxlib pos_samp.f90 utils.o ParamNames.o -L../camb -lcamb -I../camb -o ../pos_samples
! run by './pos_samples params.ini', or
! './pos_samples params.ini nchains', if running multiple chains 
! where nchains>0 is the number of chains

program pos_samp

  use AMLutils
  use IniFile
  implicit none
  integer, parameter :: maxnsamp=100000
  integer, parameter :: num_tols=10
  integer, parameter :: max_num_params = 13
  real*8, dimension(maxnsamp) :: like,vnow,wt
  real*8, dimension(:,:), allocatable :: par
  real*8  av,inc,logZ,H,dum
  real*8, dimension(num_tols) :: Ev,N,X,incs
  real, dimension(400) :: liketosort
  real :: rdum, xend, lowwt
  integer :: idum,ichain,nchain
  logical ldum
  real, dimension(max_num_params) :: params_out
  real, dimension(12) :: derived_params
  integer, dimension(max_num_params) :: params_used
  integer i,nsamp,np,numlike,j,nfile,iargc
  character(len=80) evroot,evfile,outfile,livefile
  character(len=1024) InputFile
  character(len=32) fmt
  character(len=1) char_num_chains
  logical bad, multiple_chains
!  logical :: Ini_fail_on_not_found=.false.
! --- initialize tolerances ----  
  incs(1) = 50.0
  incs(2) = 30.
  incs(3) = 25.
  incs(4) = 20.
  incs(5) = 15.
  incs(6) = 10.
  incs(7) = 5.
  incs(8) = 3.
  incs(9) = 1.
  incs(10) = -0.69

  nfile = iargc()
  if(nfile>2) then
     write(*,*) "Too many input files"
     stop
  endif
  if(nfile.eq.1) then
     call getarg(nfile,InputFile)
     nchain=1
     multiple_chains=.false.
  else
     call getarg(1,InputFile)
     call getarg(2,char_num_chains)
     multiple_chains = .true.
     read(char_num_chains,*) nchain
  endif
  call Ini_Open(InputFile, 1, bad, .false.)
  if (bad) then
     write(*,*) 'Error opening parameter file'
     stop
  endif
  Ini_fail_on_not_found = .false.

! --- get number of parameters ---
  call get_num_params(max_num_params,params_used,np,params_out)
  allocate(par(1:maxnsamp,1:np))
  evroot = Ini_Read_String('file_root')


  do ichain=1,nchain
     if(multiple_chains) then
        write(char_num_chains,'(i1)') ichain
        outfile=TRIM(evroot)//'_postsamp_'//TRIM(char_num_chains)//'.txt'
        evfile=TRIM(evroot)//'_'//TRIM(char_num_chains)//'ev.dat'
        livefile = TRIM(evroot)//'_'//TRIM(char_num_chains)//'.txt'
     else
        outfile=TRIM(evroot)//'postsamp.txt'
        evfile=TRIM(evroot)//'ev.dat'
        livefile = TRIM(evroot)//'.txt'
     endif
     ! read in ev file
     open(1,file=evfile,form='formatted',status='old')
     i=0
     write(fmt,'(a,i2.2,a)')  '(',np+4,'E19.8,i8)' 
     do 
        i=i+1
        read(1,fmt,end=20) par(i,1:np),like(i),vnow(i),logZ,inc,numlike  
        ! get results at the jth tolerance level
        do j=1,10
           if(abs(inc-incs(j)).lt.0.1) then
              Ev(j)=logZ
              N(j)=numlike
              X(j)=vnow(i)
           end if
        end do
     end do
20   close(1)
     nsamp=i-2
     par(nsamp+1:maxnsamp,1:np)=0.
     
     ! read in final remaining points
     open(1,file=livefile,form='formatted',status='old')
     read(1,*) rdum
     read(1,*) xend
     read(1,*) idum
     read(1,*) ldum
     idum=400
     do i=1,idum
        read(1,*) par(i,1:np),liketosort(i)
     end do
     close(1)
     
     call piksrt(idum,liketosort)
     
     do i=nsamp+1,nsamp+idum
        like(i)=liketosort(i-nsamp)
        vnow(i)=xend*(idum-(i-nsamp))/real(idum)
        par(i,1:np)=par(i-nsamp,1:np)
     end do
     nsamp=nsamp+idum
     H = 0.
     ! find normalization and information
     do i=1,nsamp
        wt(i)=exp(like(i)-Ev(10))*(vnow(i-1)-vnow(i+1))/2.
        H = H+((like(i)-Ev(10)))*wt(i)
     end do
     av=sum(wt)/real(nsamp-1)
     wt=wt/av
     lowwt=1.e-30
     print*, "Information content", H
     print*, "ln Evidence", Ev(10)
     derived_params(:) = 0.0
     open(2,file=outfile,form='formatted',status='new')
     write(fmt,'(a)')  '(27E16.7)'
     do i=1,nsamp
        do j=1,np
           params_out(params_used(j)) = par(i,j)
        enddo
        if(wt(i).lt.lowwt) wt(i)=0.
        write(2,fmt) wt(i),-like(i),params_out, derived_params
     end do
     close(2)
     write(*,*) 'written posterior samples file ', outfile
  enddo
  deallocate(par)

end program pos_samp

subroutine piksrt(n,arr)
  integer n
  real arr(n)
  integer i,j
  real a
  do j=2,n
     a=arr(j)
     do i=j-1,1,-1
        if(arr(i).le.a) goto 10
        arr(i+1)=arr(i)
     end do
     i=0
10     arr(i+1)=a
  end do
  return
end subroutine piksrt

subroutine get_num_params(max_num_params,params_used,np,params_out)
  use AMLutils
  use IniFile
  use ParamNames
  implicit none
  integer, intent(in) :: max_num_params
  integer, intent(out) :: np,params_used(max_num_params)
  real, intent(out) :: params_out(max_num_params)
  integer i, ncols
  character(len=5000) InLine
  real center,PMin,PMax,wid,PWidth
  character(LEN=120) fname, parameter_names_file, pname
  Type(TParamNames) :: NameMapping

  np = 0
  params_used(:) = 0
  params_out(:) = 0.0
  parameter_names_file = 'params_CMB.paramnames'
  call ParamNames_Init(NameMapping,parameter_names_file) 
  
  ncols = NameMapping%nnames+2
  


  do i=1,max_num_params
     write(pname,'(a,a,a)') 'param[',TRIM(NameMapping%name(i)),']'
     InLine = Ini_Read_String(pname, .true.)
     read(InLine, *, err = 100) center, PMin, PMax, wid, PWidth
      if (PMax < PMin) then
        write(*,*) 'You have param Max  < Min'
        stop
     endif
     if (PWidth /= 0) then
        np = np + 1
        params_used(np) = i
     else
        params_out(i) = center
     end if
  enddo

100  return
end subroutine get_num_params


