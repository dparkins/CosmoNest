! Do nested sampling algorithm to calculate Bayesian evidence
! March 2006

module Nested
  use ParamDef
  use CalcLike
  use Random
  use propose
  implicit none

  integer, parameter :: nlive = 500 ! Number of live points

  real*8, parameter :: ef=1.7 ! enlargement factor
  real, parameter :: nest_tol=5.e-1 ! tolerance at end

  Type(ParamSet) CurParams

contains
  
  subroutine Nestsample
    implicit none
    integer, parameter :: nstep = 10000 ! maximum # of replacements
 
    integer np
    integer i,indx(1),j,ix,ndata,k,n,numlike
    real*8 logZ, vold,vnow,vnext,shrink,h,inc, logEv, X
    real cosparams(num_params),lowlike,likeact,lnew
    real, dimension(num_params_used, num_params_used) :: covmat
    real, dimension(num_params_used, nlive) :: p
    real, dimension(num_params_used) :: covmat_diag
    real, dimension(nlive) :: l
    real, dimension(num_params_used) :: pact,radius,prior_min,prior_max
    Type(ParamSet) Trial, Cur
    logical live_points, mask(num_params_used), accept
    character(len=20) fmt

 
    np=num_params_used
    ! use num_params_used and their ranges as read in from the input file
    accept = .false.
    cosparams = Scales%center
    Trial%P = cosparams
    Cur%P = cosparams
    lnew = -GetLogLike(Trial)
    call AcceptReject(accept,Trial%info,Cur%info)
    if(feedback>1) write(*,*) 'LNEW TEST', lnew, Trial%P, Trial%Info%LastParams%InitPower, &
         Trial%Info%LastParams%omb, Trial%Info%LastParams%omc, Trial%Info%LastParams%omv
    lnew = -GetLogLike(Trial)
    if(feedback>1) write(*,*) 'LNEW TEST', lnew, cosparams

    do i=1,np
       prior_min(i)=Scales%Pmin(params_used(i))
       prior_max(i)=Scales%Pmax(params_used(i))
    enddo
    mask(:)= .true.
    
!! get total number of data points being used
!ndata=0
!do i=1, num_datasets
!   ndata=ndata+datasets(i)%num_points
!end do
!ndata=ndata+899   ! for WMAP
!write(*,*) 'TOTAL NUMBER OF DATA POINTS BEING USED: ', ndata

    if (Feedback > 0) then
       do i=1,np
          write(*,'(a,i2,a,f8.3,a,f8.3)') "Param ",i," min = ",prior_min(i),&
               " max = ",prior_max(i)
          radius(i) = (prior_max(i)-prior_min(i))/2.
       enddo
    endif

    numlike=0
    write(fmt,'(a,i2.2,a)')  '(',np+1,'E19.8)'
    read(outfile_handle,'(E19.8)') logEv
    read(outfile_handle,'(E19.8)') X
    read(outfile_handle,'(i8)') numlike
    read(outfile_handle,'(l)') live_points
    if(feedback.ge.0) print*, logEv, X, numlike, live_points
    if(live_points) then
       do i=1,nlive
          read(outfile_handle,fmt) p(1:np,i), l(i)
       end do
       if(feedback >0) then
          write(*,*) logEv
          write(*,*) X
          write(*,*) numlike
          write(*,*) live_points
          write(*,*) 'got live points from file'
       endif
    else
       call gen_initial_live(np,p,l,prior_min,radius,mask,cosparams)
       numlike=nlive
       live_points = .true.
       rewind(outfile_handle)
       write(outfile_handle,'(E19.8)') logEv
       write(outfile_handle,'(E19.8)') X
       write(outfile_handle,'(i8)') numlike
       write(outfile_handle,'(l)') live_points
       do i=1,nlive
          write(outfile_handle,fmt) p(1:np,i), l(i)
       end do
       if (flush_write) call FlushFile(outfile_handle)
       if(feedback > 0) write(*,*) 'generated live points'
    endif

    vold   = X ! start volume
    vnow   = vold
    logZ = logEv
    i = int(log(X)/log(real(nlive)/real(nlive+1)))
    do while (i <= nstep)
       if(Feedback > 1) write(*,*) "(Re)starting loop"
       i = i+1 
       indx   = minloc(l)  ! lowest like
       pact(:) = p(:,indx(1))
       lowlike = l(indx(1))
       likeact = lowlike

!       shrink=exp(-1.0/real(nlive))

! option 1
! sample log(t) with distribution <log t> = -1/N, dev(log t)=1/N
!        v=(Gaussian1()-1./real(nlive))*1./real(nlive)
!        shrink=exp(v)

       shrink = real(nlive)/real(nlive+1)
       vold = vnow        ! old volume is last steps volume
       vnow = vnow*shrink ! new volume shrinks
       if(mod((i-1),1)==0) then
          call calc_covmat(np,p,covmat)
          covmat_diag=0.
          call Matrix_diagonalize(covmat,covmat_diag,np)
       end if
       call getnewp_box(np,p,l,indx(1),covmat_diag,covmat,prior_min,&
            radius,mask,numlike)

!       p(:,indx(1)) = pact(:)
!       l(indx(1)) = likeact
       vnext = shrink*vnow
!       h    = vold - vnow ! delta-volume = dx
       h = 0.5*(vold-vnext)
       logZ = PLUS(logZ,(lowlike+log(h)))
       if(logZ.gt.-logZero) then
          inc=maxval(l)+log(vnow)-logZ   ! another possibility
          if(inc.lt.(log(nest_tol))) exit
       endif
       write(fmt,'(a,i2.2,a)')  '(',np+1,'E19.8)'
       if (mod(i,1)==0) then
          if(Feedback > 1) then
             write(*,'(a,3E19.8)') 'log evidence so far', logZ,inc, lowlike
             write(*,'(a,2i8)') 'number of replacements', i, numlike
             write(*,'(a,2E19.8)') 'prior mass remaining', vnow, log(vnow)
          endif
          rewind(outfile_handle)
          write(outfile_handle,'(E19.8)') logZ
          write(outfile_handle,'(E19.8)') vnow
          write(outfile_handle,'(i8)') numlike
          write(outfile_handle,'(l)') live_points
          do ix=1,nlive
             write(outfile_handle,fmt) p(1:np,ix), l(ix)
          end do
          if (flush_write) call FlushFile(outfile_handle)
          write(fmt,'(a,i2.2,a)')  '(',np+4,'E19.8,i8)'
          write(indepfile_handle,fmt) pact,likeact,vnow,logZ,inc,numlike
          if (flush_write) call FlushFile(indepfile_handle)
       end if
       if(Feedback > 1) write(*,*) "about to restart loop"
    enddo
    do i = 1,nlive
       logZ = PLUS(logZ,(l(i)+log(vnow/real(nlive))))
    enddo
    write(indepfile_handle,'(e16.7)') logZ
    close(indepfile_handle)


  end subroutine Nestsample

  function PLUS(x,y)

  implicit none
  real*8 x,y,PLUS

  if(x.gt.y) then
     PLUS = x + log(1.+exp(y-x))
  else
     PLUS = y + log(1.+exp(x-y))
  endif

  end function

  subroutine getnewp_box(np,p,l,index,covmat_diag,covmat,prior_min,&
       radius,mask,numlike)
    implicit none
    
    integer, intent(in)   :: np,index
    integer, intent(inout) :: numlike
    real, intent(inout) :: covmat_diag(np),covmat(np,np)
    real, intent(in) :: prior_min(np),radius(np)
    real, intent(inout) :: p(np,nlive),l(nlive)
    logical, intent(in) :: mask(np)
    real trans_p(np,nlive), ulimit(np), llimit(np),try_p(np)
    real*8 u(np),lnew,like,lboundary
    real pnew(np),cosparams(num_params)
    Type(ParamSet) Trial
    logical accept
    integer i,j,num
    
    accept = .false.
    

    cosparams = Scales%center
    Trial = CurParams
    Trial%P = cosparams
    lnew = -GetLogLike(Trial)
    call AcceptReject(accept,CurParams%Info,Trial%Info)
    if(Feedback > 1) write(*,*) 'LNEW TEST', lnew, cosparams


    if(Feedback > 1) write(*,*) "Trying from prior"
    lboundary = l(index)
    call getrandom(np,pnew,radius,prior_min,mask)
! ... check to see if in prior ...
    if(inprior(np,pnew,prior_min,radius,mask)) then
       cosparams = Scales%center
       do i=1,np
          cosparams(params_used(i)) = pnew(i) 
       enddo
       Trial%P = cosparams
       lnew = -GetLogLike(Trial)
       if(min(lnew,0.).eq.0.) then
          !     if(ISNAN(lnew)) then
          !        write(*,*) 'a NAN here' 
          lnew=-LogZero
          !        cycle  
       end if
       
       if(abs(lnew).eq.LogZero) then
          lnew=-LogZero*(1.+0.0001*ranmar())
       end if

       if(Feedback > 1) write(*,'(2(a,e19.8))') "lnew ", lnew, " lboundary ", lboundary
       
       call AcceptReject(accept,CurParams%Info,Trial%Info)
       numlike=numlike+1

    
       if (lnew>lboundary) then
          p(:,index) = pnew(:)
          l(index) = lnew
          return
          !      exit ! if yes, stop
       end if
    endif

    if(Feedback > 1) write(*,*) "Creating rotated population"
    continue
    if(Feedback > 1) write(*,*) "No, honestly!"
    

11  do i=1,nlive
       !     trans_p(:,i) = MATMUL(covmat,p(:,i))
       trans_p(:,i) = MATMUL(transpose(covmat),p(:,i))
    enddo
   
    do i=1,np
       ulimit(i) = maxval(trans_p(i,:))
       llimit(i) = minval(trans_p(i,:))
    enddo
    num=0
    
    if(Feedback > 1) write(*,*) "Trying from spheroid"
    do

! pick point in cube
!     do i=1,np
!        u(1) = ranmar()
!        try_p(i) = ((ulimit(i)-llimit(i))*(ef*(u(1)-0.5)+0.5))+llimit(i)
!     enddo
!     pnew= MATMUL(try_p,covmat)

! or in spheroid
       call spheroid(np,u,ef)
       try_p(:)=u(:)*(ulimit(:)-llimit(:))+llimit(:)
       
       !     pnew= MATMUL(try_p,covmat)
       pnew=MATMUL(covmat,try_p)
       if(.not.inprior(np,pnew,radius,prior_min,mask)) cycle
       do i=1,np
          cosparams(params_used(i)) = pnew(i) 
       enddo
       if(Feedback>1) then
          do i=1,np
             write(*,'(a,i2.2,a,f8.4)') "Param", i, " = ", pnew(i)
          enddo
       endif
       Trial%P = cosparams
       lnew = -GetLogLike(Trial)
       call AcceptReject(accept,CurParams%Info,Trial%Info)
       numlike=numlike+1
       num=num+1
       if(min(lnew,0.).eq.0.) then
          !     if(ISNAN(lnew)) then
          !        write(*,*) 'a NAN here' 
          lnew=-LogZero
       end if
       if(abs(lnew).eq.LogZero) then
          lnew=-LogZero*(1.+0.0001*ranmar())
       end if
       if(Feedback.ge.1) write(*,*) "LNEW =", lnew
       
       if(lnew>lboundary) then
          p(:,index) = pnew(:)
          l(index) = lnew
          !        write(*,*) '**point accepted**',pnew, lnew, lboundary,num
          return
          !      exit ! if yes, stop
       end if
       
       !     write(*,*) 'point not accepted:', pnew
       !     write(*,*) 'likelihood and numlike', lnew, numlike,lboundary
       
       if(num.gt.50) then
          covmat(:,:) = 0.0
          do j=1,np
             covmat(j,j) = 1.0
          enddo
          goto 11
       endif
       
       cosparams = Scales%center

    enddo
    STOP "should never reach this point"
    
  
  end subroutine getnewp_box

  subroutine spheroid(np,u,ef)
   
    implicit none
    integer np,i
    real*8 u(np), mod, ef
    
    do i=1,np
       u(i)=ranmar()
    end do
    u = (u-0.5)
    mod=sum(u(:)**2)
    do while(mod.gt.0.25)
       do i=1,np
          u(i)=ranmar()
       end do
       u = (u-0.5)
       mod=sum(u(:)**2)
    end do
    u = u*ef
    u=u+0.5
    
    return
  end subroutine spheroid
  

  function inprior(np,p,radius,prior_min,mask)
    
    implicit none
    logical inprior
    integer np, i
    real p(np), radius(np), prior_min(np), prob, rnumber
    logical mask(np)
    
    inprior = .true.
    do i=1,np
       if(mask(i)) then
          if(p(i).lt.prior_min(i)) inprior = .false.
          if(p(i).gt.prior_min(i)+2.*radius(i)) inprior = .false.
       else
          prob= exp(-(p(i)-prior_min(i)-2.0*radius(i))/0.05)
          rnumber = ranmar()
          if(prob.le.rnumber) then
             ! print*, p(i), prior_min(i)+2.0*radius(i),  prob, rnumber
             inprior = .false.
          endif
       endif
    enddo
    
  end function inprior

  subroutine gen_initial_live(np,p,l,prior_min,radius,mask,cosparams)
    implicit none
    integer np,i,j
    Type(ParamSet) Trial
    real*8 like_test
    real p(np,nlive), l(nlive),radius(np),prior_min(np),cosparams(num_params)
    logical accept, mask(np)
    accept = .false.

    cosparams = Scales%center
    Trial = CurParams
    Trial%P = cosparams
    like_test = -GetLogLike(Trial)
    call AcceptReject(accept,CurParams%Info,Trial%Info)
    if(feedback > 1) write(*,*) 'LNEW TEST', like_test, Trial%P, Trial%Info%LastParams%InitPower, &
         Trial%Info%LastParams%omb, Trial%Info%LastParams%omc, Trial%Info%LastParams%omv

    
    do i = 1,nlive

       call getrandom(np,p(:,i),prior_min,radius,mask)  ! start points
       cosparams = Scales%center
       do j=1,np
          cosparams(params_used(j)) = p(j,i)
       enddo
       Trial%P = cosparams
       if(feedback.ge.2) print*, "Params = ", Trial%P, Trial%Info%LastParams%InitPower, &
         Trial%Info%LastParams%omb, Trial%Info%LastParams%omc, Trial%Info%LastParams%omv
       l(i) = -GetLogLike(Trial)
       call AcceptReject(accept,CurParams%Info,Trial%Info)
       if(feedback.ge.1) write(*,*) "Like =", l(i)
       if(abs(l(i)).eq.LogZero) then
          l(i)=-LogZero*(1.+0.0001*ranmar())
       end if
       
       !     l(i) = lfunc(np,p(:,i))

       !     if(ISNAN(l(i))) then
       if(min(l(i),0.).eq.0.) then
          !        write(*,*) 'a NAN here', i 
          l(i)=-LogZero*(1.+0.0001*ranmar())
          !        cycle  ! will cause sampling bias
       end if
     
    enddo

    return
  end subroutine gen_initial_live

  subroutine calc_covmat(np,p,covmat)
    implicit none
    integer, intent(in) :: np
    real, intent(in) :: p(np,nlive)
    real, intent(out) :: covmat(np,np)
    real mean(np)
    integer i2,j2,ip
    mean=0.
    do i2=1,np
       do j2=1,nlive
          mean(i2)=mean(i2)+p(i2,j2)/real(nlive)
       end do
    end do
    covmat=0.
    do i2=1,np
       do j2=1,np
          do ip=1,nlive
             covmat(i2,j2)=covmat(i2,j2)+((p(i2,ip)-mean(i2))*&
                  (p(j2,ip)-mean(j2)))/real(nlive)
          end do
       end do
    end do
    return
    
  end subroutine calc_covmat
  
  subroutine getrandom(n,x,prior_min,radius,mask)
    implicit none
    integer, intent(in) :: n
    real, intent(out) :: x(n)
    real, intent(in) :: prior_min(n),radius(n)
    real u(n), g(n)
    logical, intent(in) :: mask(n)
    integer i
    
    ! --- uniform prior ----
    do i = 1,n
       u(i) = ranmar()
    enddo
    
    ! --- gaussian prior ---
    do i = 1,n
       g(i) = Gaussian1() 
    enddo
   
    do i=1,n
       if(mask(i)) then
          x(i) = 2.*radius(i)*u(i)+prior_min(i)
       else
          x(i) = radius(i)*g(i)+(prior_min(i)+radius(i))
       endif
    enddo
    return
  end subroutine getrandom
  
end module Nested
