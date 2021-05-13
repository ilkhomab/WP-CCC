!subroutine gaulag(xx,ww,n,alff)
!
!!use data_base,only:key_num
!implicit none
!integer, parameter :: p16 = selected_real_kind(16)
!integer, parameter :: id = selected_real_kind(8)
!integer n
!real(kind=id),intent(out),dimension(n):: ww,xx
!integer,parameter:: maxit=100
!integer i,its,j
!real(kind=id),intent(in)::alff
!real(kind=p16),dimension(n):: w,x
!real(kind=p16)::eps,z,z1,p1,p2,p3,pp,ai,t1,t2,t3,tmp,arg,alf,tw
!!real(kind=p16) :: qgamma
!
!alf=real(alff,p16)
!
!eps=exp((-16+4)*log(10._p16))
!
!do 13 i=1,n
!  if(i.eq.1)then
!    z=((1._p16)+alf)*((3._p16)+0.92_p16*alf)/((1._p16)+2.4_p16*n+1.8_p16*alf)
!  else if(i.eq.2)then
!    z=z+((15._p16)+6.25_p16*alf)/((1._p16)+0.9_p16*alf+2.5_p16*n)
!  else
!    ai=i-2
!    z=z+(((1._p16)+2.55_p16*ai)/(1.9_p16*ai)+1.26_p16*&
!    &ai*alf/((1._p16)+3.5_p16*ai))*(z-x(i-2))/((1._p16)+0.3_p16*alf)
!  endif
!  do 12 its=1,maxit
!    p1=(1._p16)
!    p2=(0._p16)
!    do 11 j=1,n
!      p3=p2
!      p2=p1
!      p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
!11 continue
!    pp=(n*p1-(n+alf)*p2)/z
!    z1=z
!    z=z1-p1/pp
!    if(abs(z-z1).le.eps)goto 1
!12 continue
!!  pause 'too many iterations in gaulag'
!1 x(i)=z
!t3=alf+real(n,p16)
!
!!call mpgamma(t3,t1)
!t1=gamma(t3)
!t2=gamma(real(n,p16))
!
!!w(i)=-exp(gammln(alf+n)-gammln(float(n)))/(pp*n*p2)
!w(i)=-t1/t2/(pp*real(n,p16)*p2)
!!call mpwrite(6,x(i),w(i))
!13    continue
!
!!  print*,'weights are ready'
!!  pause
!
!  tmp=(0._p16)
!  do i=1,n
!
!    arg=sqrt((x(i)+(1._p16))*(x(i)+(1._p16))-0.5_p16)
!    tmp=tmp+sin(arg)/arg*w(i)!*exp(-x(i))
!
!  enddo
!
!!      print*,'test of gaulag integration of x exp(-x)'
!!      print*,'anaytical integration=1'
!!      print*,'numerical integration='
!
!do i=1,n
!  xx(i)=real(x(i),id)
!!  if(key_num==0)then
!    ww(i)=exp(log(w(i))+x(i))
!!    tw=w(i)*exp(x(i))
!!    ww(i)=real(tw,id)
!!  else
!!    ww(i)=real(w(i),id)
!!  endif
!!  print*,xx(i),ww(i)
!enddo
!
!return
!end


subroutine gaulag(xx,ww,n,alff)

!use precision_type
!use data_base,only:key_num
implicit none
integer, parameter :: p16 = selected_real_kind(16)
integer, parameter :: id = selected_real_kind(8)
integer n
real(kind=id),intent(out),dimension(n):: ww,xx
real(kind=id),dimension(2*n):: wf
integer(kind=id),dimension(2*n)::iwf
integer,parameter:: maxit=100
integer i,its,j,ier
real(kind=id),intent(in)::alff
real(kind=p16),dimension(n):: w,x
real(kind=p16)::eps,z,z1,p1,p2,p3,pp,ai,t1,t2,t3,tmp,arg,alf,tw
!real(kind=p16) :: qgamma

call cgqf(n,xx,ww,5,0d0,0d0,0._id,1._id,0,2*n,wf(1:2*n),2*n,iwf(1:2*n),ier)
do i=1,n
!  ww(i)=exp(xx(i))*ww(i)
  ww(i)=exp(xx(i)+log(ww(i)))
enddo
print*,'ier=',ier
if (ier.eq.0) return

alf=real(alff,p16)

eps=exp((-16+4)*log(10._p16))

do 13 i=1,n
  if(i.eq.1)then
    z=((1._p16)+alf)*((3._p16)+0.92_p16*alf)/((1._p16)+2.4_p16*n+1.8_p16*alf)
  else if(i.eq.2)then
    z=z+((15._p16)+6.25_p16*alf)/((1._p16)+0.9_p16*alf+2.5_p16*n)
  else
    ai=i-2
    z=z+(((1._p16)+2.55_p16*ai)/(1.9_p16*ai)+1.26_p16*&
    &ai*alf/((1._p16)+3.5_p16*ai))*(z-x(i-2))/((1._p16)+0.3_p16*alf)
  endif
  do 12 its=1,maxit
    p1=(1._p16)
    p2=(0._p16)
    do 11 j=1,n
      p3=p2
      p2=p1
      p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
11 continue
    pp=(n*p1-(n+alf)*p2)/z
    z1=z
    z=z1-p1/pp
    if(abs(z-z1).le.eps)goto 1
12 continue
!  pause 'too many iterations in gaulag'
1 x(i)=z
t3=alf+real(n,p16)

!call mpgamma(t3,t1)
t1=gamma(t3)
!call mpgamma(mpreal(n),t2)
t2=gamma(real(n,p16))

!w(i)=-exp(gammln(alf+n)-gammln(float(n)))/(pp*n*p2)
w(i)=-t1/t2/(pp*real(n,p16)*p2)
!call mpwrite(6,x(i),w(i))
13    continue

!  print*,'weights are ready'
!  pause

  tmp=(0._p16)
  do i=1,n

!    tmp=tmp+sin(x(i))*sqrt(x(i)*x(i)+mpreal(1))*w(i)
!    tmp=tmp+sin(x(i))*exp(mpreal('1.5')*log(x(i)))*w(i)
!    tmp=tmp+sin(x(i))*sqrt(mpreal('1.')+mpreal('1.')/x(i)/x(i))*w(i)
    arg=sqrt((x(i)+(1._p16))*(x(i)+(1._p16))-0.5_p16)
    tmp=tmp+sin(arg)/arg*w(i)!*exp(-x(i))

  enddo

!      print*,'test of gaulag integration of x exp(-x)'
!      print*,'anaytical integration=1'
!      print*,'numerical integration='

do i=1,n
  xx(i)=real(x(i),id)
!  if(key_num==0)then
    tw=w(i)*exp(x(i))
    ww(i)=real(tw,id)
!  else
!    ww(i)=real(w(i),id)
!  endif
!  print*,xx(i),ww(i)
enddo

return
end
