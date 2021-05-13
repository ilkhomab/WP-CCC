!module precision_type
!  implicit none
!  integer, parameter :: p4 = selected_real_kind(4)
!  integer, parameter :: p8 = selected_real_kind(8)
!  integer, parameter :: p16 = selected_real_kind(16)
!  !integer, parameter :: id = selected_real_kind(16)
!  integer, parameter :: id = selected_real_kind(8)
!  integer, parameter :: mpoud=2*id
!  !integer, parameter :: mpoud_p16=16
!  integer, parameter :: mpoud_p16=32
!  real(kind=id) :: mppic
!  real(kind=id),parameter :: pi=3.141592653589793238462643383279502884197_id
!  real(kind=id),parameter :: hr=27.2107_id
!  real(kind=id),parameter :: pifc=0.797884560802865355879892119869_id
!  real(kind=p16) :: mppic16
!end module precision_type
!
!module data_base
!  use precision_type
!  implicit none
!  real(kind=id)::scal,pmin,pmax
!  integer,allocatable,dimension(:,:,:)::nlm
!  real(kind=id),allocatable,dimension(:,:,:)::pgrid
!  real(kind=id),allocatable,dimension(:,:,:)::wf
!  real(kind=id),allocatable,dimension(:,:)::enc
!  real(kind=id),allocatable,dimension(:)::bin,dis
!  complex(kind=id),parameter::cu=cmplx(0._id,1._id,id)
!  complex(kind=id),parameter::zero=cmplx(0._id,0._id,id)
!  integer::nbin,np,keyc
!end module data_base


 SUBROUTINE locate(xx,n,x,j)
 use precision_type
 implicit none
 real(kind=id),dimension(n):: xx
 real(kind=id)::x
! Given an array xx(1:n), and given a value x,
!returns a value j such that x is between
!xx(j) and xx(j+1). xx(1:n) must be monotonic,
!either increasing or decreasing.
!j=0 or j=n is returned to indicate that x is out of range.
 integer:: jl,jm,ju,n,j
 jl=0                      ! Initialize lower
 ju=n+1                    ! and upper limits.
 10   if(ju-jl.gt.1)then        ! If we are not yet done,
    jm=(ju+jl)/2           ! compute a midpoint,
    if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
       jl=jm               ! and replace either the lower limit
    else
       ju=jm               ! or the upper limit, as appropriate.
    endif
    goto 10                ! Repeat until
 endif                     ! the test condition 10 is satisfied.
 if(x.eq.xx(1))then        ! Then set the output
    j=1
 else if(x.eq.xx(n))then
    j=n-1
 else
    j=jl
 endif
 return                    ! and return.
 END



subroutine get_cwf(charge)
  use precision_type
  use mpi,only:myid
  use data_base
  implicit none
  integer::i,j,ir,nl,kfn,ifail,mode,l,ncust,ic,jc,ntemp
  real(kind=id)::r,p,wp,arg,pp,binmin,charge
  integer,dimension(nbin)::npbin
  real(kind=id),dimension(nbin)::sum,en!,dis
  real(kind=id),dimension(nbin+1)::enbin,binlim
  real(kind=id),dimension(0:lmax_t)::radial
  real(kind=id)::xx,eta,ki,bess_wp
  real(kind=id),dimension(0:lmax_t+1)::fc,gc,fcp,gcp,sig
  real(kind=id),allocatable,dimension(:)::custen



!  call vec_linspace(nbin+1,k_min,k_max,binlim(1:nbin+1))

!  call tcheb_space ( nbin+1,scal, bin )
!  call cheb_grid ( nbin+1,8._id, bin )

!  call set_expgrid(nbin,scal0,scal,bin)
!  call set_customgrid(nbin,scal0,scal,bin)
!  if(key_capt==1)then
    call vec_linspace(nbin+1,k_min,k_max,binlim(1:nbin+1))
!    call set_expgrid(nbin+1,k_min,k_max,binlim(1:nbin+1))

    if(ejenfile)then
      open(666,file='ejenfile',status='old')
      read(666,*)ncust
      allocate(custen(1:ncust))
      read(666,*)custen(1:ncust)
      close(666)
      do ic=1,ncust
        custen(ic)=sqrt(custen(ic)/13.6058_id)
        call locate(binlim(1:nbin+1),nbin+1,custen(ic),jc)
        binlim(jc+1)=0.5_id*(sqrt(3._id*(4*custen(ic)*custen(ic)-binlim(jc)*binlim(jc)))-binlim(jc))
      enddo
    endif

    

!  else
!    call set_customgrid2(binlim(1:nbin+1))
!  endif

!  binlim(1)=0.5d0*bin(1)!+0.5d0*k_min

!  dis(1)=bin(2)-bin(1)

if(myid==0)then
  print*,binlim(1)*binlim(1)*13.6058
  print*,binlim(nbin+1)*binlim(nbin+1)*13.6058
!  stop

  print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print*,'In cwf8'
  print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
endif

  do i=1,nbin
     bin(i)=sqrt((binlim(i)*binlim(i)+binlim(i+1)*binlim(i+1)+binlim(i)*binlim(i+1))/3.e0_id)!0.5d0*(binlim(i)+binlim(i+1))
if(myid==0)then
     print*,i,bin(i)
endif
     dis(i)=binlim(i+1)-binlim(i)
!     binlim(i)=binlim(i-1)+dis(i)
  enddo

!  dis(nbin)=bin(nbin)-bin(nbin-1)
!  binlim(nbin)=0.5d0*(bin(nbin-1)+bin(nbin))
!  binlim(nbin+1)=binlim(nbin)+dis(nbin)


  binmin=minval(dis(1:nbin))

  do i=1,nbin
    ntemp=int(np*dis(i)/binmin)
    if(mod(ntemp,2)/=0)ntemp=ntemp+1
    if(ntemp.gt.npmax)ntemp=npmax
    npbin(i)=ntemp  
  enddo 


!  stop



  do i=1,nbin
    call gauleg(binlim(i),binlim(i+1),pgrid(i,1:npbin(i),1),pgrid(i,1:npbin(i),2),npbin(i))
    ki=bin(i)
    enc(i,0:lmax_t)=ki*ki
!    enbin(i)=bin(i)*bin(i)
if(myid==0)then
    print*,i,enc(i,0)*13.6058e0_id
endif
  enddo
!  enbin(nbin+1)=bin(nbin+1)*bin(nbin+1)
if(probfile)then
  goto 20 
endif

  nl=1
  mode=1
  kfn=0
  ifail=0

!  print*,'Begin'
!
!  do ir=1,meshr
!    print*,ir
!    r=rmesh(ir,1)
!!    if(r.gt.1.e-2_id)then
!!      call coul90(1.e-2_id*r,-100._id,0._id,4,fc,gc,fcp,gcp,kfn,ifail)
!!!      CALL Hydrogenic_Continuum_Wave(1.e-2_id, r, 1._id, 1, Radial)  
!!      print*,'ifail',ifail
!!      write(11,30)rmesh(ir,1),(fc(i)/r,i=0,4)
!!    endif
!!        fc=0d0
!!        gc=0d0
!!        fcp=0d0
!!        gcp=0d0
!!        sig=0d0
!!        ifail=0
!!        pp=0.01_id
!!        eta=-1._id/pp
!!!        if(r.gt.2._id)then
!!          call coul90(pp*r,eta,0._id,lmax,fc,gc,fcp,gcp,kfn,ifail)
!!!        else
!!          CALL Hydrogenic_Continuum_Wave(pp, r, 1._id, lmax, radial)  
!!!        endif
!!!        if(ifail==0)then
!!          write(31,30)r,(fc(l)/r,l=0,lmax)
!!          write(32,30)r,(radial(l)/r,l=0,lmax)
!!!        endif
!    do i=1,nbin
!      wf(i,0:lmax,ir)=0._id
!      do j=1,npbin(i)
!        p=pgrid(i,j,1)
!        arg=p*r
!        wp=pgrid(i,j,2)
!        fc=0d0
!        gc=0d0
!        fcp=0d0
!        gcp=0d0
!        sig=0d0
!        ifail=0
!!        if(arg.gt.1.e-1_id.and.p.gt.1.e-1_id)then
!        if(arg.gt.1.e-3_id)then
!!            print*,arg,p
!!            print*,'coul90'
!            call coul90(arg,-1._id/p,0._id,lmax,fc,gc,fcp,gcp,kfn,ifail)
!        else
!!            print*,'Peng'
!            CALL Hydrogenic_Continuum_Wave(p, r, 1._id, lmax, fc)  
!        endif
!        wf(i,0:lmax,ir)=wf(i,0:lmax,ir)+wp*fc(0:lmax)
!!        if(ifail==0)wf(i,0:lmax,ir)=wf(i,0:lmax,ir)+wp*fc(0:lmax)
!!        xx=p*r
!!        eta=-1._id/p
!!        if(r.gt.1.e-2_id.and.p.gt.1.e-2_id)then
!!          call coulcc(xx,eta,cmplx(0._id,0._id,id),nl,fc,gc,fcp,gcp,sig,mode,kfn,ifail)
!!          CALL Hydrogenic_Continuum_Wave(p, r, 1._id, 1, Radial)  
!!          wf(i,0,ir)=wf(i,0,ir)+wp*radial(0)
!!        endif
!      enddo
!      wf(i,0:lmax,ir)=pifc*wf(i,0:lmax,ir)/sqrt(bin(i+1)-bin(i))!/r
!    enddo
!  enddo
!
!  print*,'done'

if(myid==0)then
  print*,'Begin'
endif

!$omp parallel do default(none)&
!$omp private(r,p,arg,wp,fc,gc,fcp,gcp,sig,ifail,kfn)&
!$omp shared(meshr,rmesh,nbin,wf,wfres,npbin,pgrid,lmax_t,bin,dis,binlim,charge)&
!$omp schedule(dynamic)collapse(2)
  do ir=1,meshr
    do i=1,nbin
      r=rmesh(ir,1)
      wf(i,0:lmax_t,ir)=0._id
      wfres(i,0:lmax_t,ir)=0._id
      do j=1,npbin(i)
        p=pgrid(i,j,1)
        arg=p*r
        wp=pgrid(i,j,2)
        fc=0d0
        gc=0d0
        fcp=0d0
        gcp=0d0
        sig=0d0
        ifail=0
        kfn=0
        if(arg.gt.1.e-3_id.and.arg.lt.1.e4_id)then
           call coul90(arg,-charge/p,0.d0,lmax_t+1,fc,gc,fcp,gcp,kfn,ifail)
        else
           CALL Hydrogenic_Continuum_Wave(p, r, charge, lmax_t+1, fc(0:lmax_t))  
        endif
        wf(i,0:lmax_t,ir)=wf(i,0:lmax_t,ir)+wp*fc(0:lmax_t)
        wfres(i,0:lmax_t,ir)=wfres(i,0:lmax_t,ir)+wp*fc(0:lmax_t)*p*p/2d0
      enddo
      wf(i,0:lmax_t,ir)=pifc*wf(i,0:lmax_t,ir)/sqrt(dis(i))!/r
      wfres(i,0:lmax_t,ir)=pifc*wfres(i,0:lmax_t,ir)/sqrt(dis(i))!/r
    enddo
  enddo
!$omp end parallel do

!  do ir=1,meshr
!    write(666,*)rmesh(ir,1),wf(5,1,ir)
!  enddo

if(myid==0)then
  print*,'done'

!  stop

  print*,bin
endif


  do l=0,lmax_t
      sum(1:nbin)=0._id
      do ir=1,meshr
          r=sqrt(rmesh(ir,3))
!          write(30,30)rmesh(ir,1),(wf(i,l,ir),i=1,nbin)
          sum(1:nbin)=sum(1:nbin)+r*r*wf(1:nbin,l,ir)*wf(1:nbin,l,ir)
      enddo
if(myid==0)then
      print*,'l=',l
      print*,sum(:)
      write(30,*)
      write(30,*)
endif
  enddo

if(myid==0)then
  open(32,file='energies')
  print*,'Energies:'
endif

if(myid==0)then
  do i=1,nbin
    print*,i,enc(i,0),enc(i,0)*hr/2.d0
    write(32,30)enc(i,0),enc(i,0)*hr/2.d0
  enddo
  close(32)
endif
!stop

if(myid==0)then
  open(32,file='bins')
  do i=1,nbin+1
    write(32,30)binlim(i)*binlim(i)
  enddo
  close(32)
endif

20 continue

  if(ejenfile)then
    deallocate(custen)
  endif 

 30   format(2000g14.6)
end subroutine get_cwf

subroutine tcheb_space ( n,s, a )
  use precision_type
  implicit none
  integer(kind=4),intent(in):: n
  real(kind=id),intent(out):: a(n)
  real(kind=id),intent(in):: s
  integer ( kind = 4 ):: i

  do i=1,n
    a(i)=s*tan(real(2*i-1,id)/real(4*n,id)*pi)
  enddo

  return
end

subroutine cheb_grid ( n,amax, a )
  use precision_type
  implicit none
  integer(kind=4),intent(in):: n
  real(kind=id),intent(out):: a(n)
  real(kind=id):: s
  real(kind=id),intent(in):: amax
  integer ( kind = 4 ):: i

  s=amax/tan(real(2*n-1,id)/real(4*n,id)*pi)

  do i=1,n
    a(i)=s*tan(real(2*i-1,id)/real(4*n,id)*pi)
  enddo

  return
end

subroutine set_expgrid(n,a,b,tz)
use precision_type
implicit none
integer::i,n
real*8::a,b
real*8,dimension(1:n)::tz

      
do i=1,n
  tz(i)=a*(b/a)**(dble(i-1)/dble(n-1))
enddo
      
return
end subroutine set_expgrid

subroutine set_customgrid(n,a,b,tz)
use precision_type
use data_base,only:v
implicit none
integer::i,n
real*8::a,b
real*8,dimension(1:n)::tz,wtz

call gauleg(a,v,tz(1:n/2),wtz(1:n/2),n/2)
      
do i=1,n/2
  tz(n/2+i)=v*(b/v)**(dble(i)/dble(n/2))
enddo
      
return
end subroutine set_customgrid

subroutine set_customgrid2(tz)
use precision_type
use data_base,only:v,n1,neps,n2,k_min,k_max,eps,nbin
implicit none
integer::i
real*8,dimension(1:nbin+1)::tz,wtz

do i=1,n1+1
  tz(i)=k_min*((v-0.5d0*eps)/k_min)**(dble(i)/dble(n1+1))
enddo

tz(n1+1)=0d0

call gauleg(v-0.5d0*eps,v+0.5d0*eps,tz(n1+1:n1+neps),wtz(n1+1:n1+neps),neps)
      
!do i=n1+neps+1,n1+neps+n2
do i=1,n2+1
  tz(n1+neps+i)=(v+0.5d0*eps)*(k_max/(v+0.5d0*eps))**(dble(i)/dble(n2+1))
enddo
      
return
end subroutine set_customgrid2

subroutine get_eigen(charge,n,l,lag)
  use data_base
  use precision_type
  use flogs
  
  implicit none 
  integer :: l,i,n
  real(kind=id):: beta,dnorm ,charge
  real(kind=id),dimension(1:meshr):: lag,x,fact2
  real(kind=id),dimension(0:n-l-1,1:meshr):: lagg
!!!! Laguerre/r^l
!!!! for angular z integration

  x(1:meshr)=rmesh(1:meshr,1)/real(n,id)*charge
  fact2(1:meshr)=(2*x(1:meshr))**l*2.d0/real(n*n,id)*exp(-x(1:meshr))

  beta=real(2*l+1,id)
  lagg(0,1:meshr)=1._id
  if(n-l-1>=1)then
    lagg(1,1:meshr)=(1+beta-2*x(1:meshr))
    do i=2,n-l-1
        lagg(i,1:meshr)=((beta+2*(i-1)+1-2*x(1:meshr))*lagg(i-1,1:meshr)-(beta+i-1)*lagg(i-2,1:meshr))/real(i,id)
    enddo
  endif

  dnorm=sqrt(real(faccoef(n,l)/real(n-l,id),id))

  lag(1:meshr)=lagg(n-l-1,1:meshr)*fact2(1:meshr)*dnorm*charge**1.5_id

  do i=1,meshr
    if(abs(lag(i)).lt.1.e-16_id)then
      lag(i)=0._id
    endif
  enddo


return
end subroutine get_eigen

