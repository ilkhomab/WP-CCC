module precision_type

integer, parameter :: p4 = selected_real_kind(4)
integer, parameter :: p8 = selected_real_kind(8)
integer, parameter :: p16 = selected_real_kind(8)
!integer, parameter :: id = selected_real_kind(16)
integer, parameter :: id = selected_real_kind(8)
integer, parameter :: mpoud=2*id
!integer, parameter :: mpoud_p16=16
integer, parameter :: mpoud_p16=32
real(kind=p16) :: mppic16
real(kind=id),parameter :: hr=27.2116_id
real(kind=id),parameter ::pi=3.1415926535897932384626433832795028841971693993751_id
real(kind=id),parameter :: pifc=0.797884560802865355879892119869_id
real(kind=id),parameter :: mppic=3.1415926535897932384626433832795028841971693993751_id

end module precision_type



program test
use precision_type
implicit none
integer::n,i,ier
integer,allocatable::iwf(:)
real(kind=id),allocatable::x(:),w(:),x1(:),w1(:),wf(:)

print*,'Enter number of gauss laguerre points:'
read(*,*)n
allocate(x(1:n),w(1:n))
allocate(x1(-n/2:n/2),w1(-n/2:n/2))
allocate(wf(1:2*n))
allocate(iwf(1:2*n))
!call gaulag(x,w,n,0.d0)
call gauleg(-1._id,1._id,x,w,n)

call tanhsinh(x1,w1,n)


!call cgqf(n,x1,w1,1,0d0,0d0,0._id,40._id,0,2*n,wf(1:2*n),2*n,iwf(1:2*n),ier)

!call cgqf(n,x1,w1,5,0d0,0d0,0._id,1._id,0,2*n,wf(1:2*n),2*n,iwf(1:2*n),ier)
!print*,'ierrrr=',ier

do i=1,n/2
!  if(x(i).gt.60d0)then
!    exit
!  endif
  write(77,*)x(i),0d0,w(i)
  write(88,*)x1(i),0d0,w1(i)
enddo
!n=i

print*,n,sum(exp(-x(1:n))*w(1:n))
print*,n,sum(exp(-x1(-n/2:n/2))*w1(-n/2:n/2))
!print*,n,sum(w1(1:n))

end program test

!     SUBROUTINE CGQF(NT,T,WTS,KIND,ALPHA,BETA,A,B,LO,
!     1NWF,WF,NIWF,IWF,IER)

subroutine tanhsinh(x,w,n)
  use precision_type
  implicit none
  integer::i,n
  real(kind=id),dimension(-n/2:n/2)::x,w
  real(kind=id)::h

  h=1._id/real(n-1,id)

  do i=-n/2,n/2
    x(i)=tanh(0.5_id*pi*sinh(i*h))
    w(i)=0.5_id*h*pi*cosh(i*h)/cosh(0.5_id*pi*sinh(i*h))/cosh(0.5_id*pi*sinh(i*h))
  enddo

end subroutine tanhsinh

subroutine gauleg(x1,x2,x,w,n)

  use precision_type
!  use data_base,only:num_dig
  implicit none
  real(kind=id) x1,x2
  real(kind=id),dimension(n):: x,w
  real(kind=id),dimension(2*n) :: wf
  integer,dimension(2*n):: iwf
  integer::m,i,j,n
  real(kind=id)::eps,z,z1,p1,p2,p3,pp,xl,xm,tmp

  eps=exp((-mpoud+4)*log(real(10,id)))

  m=(n+1)/2
  xm=0.5_id*(x2+x1)
  xl=0.5_id*(x2-x1)
  do 12 i=1,m
         z=cos(mppic*(i-0.25_id)/(n+0.5_id))
 1       continue
         p1=real(1,id)
         p2=real(0,id)
         do 11 j=1,n
            p3=p2
            p2=p1
            p1=((real(2,id)*j-real(1,id))*z*p2-(j-real(1,id))*p3)/real(j,id)
 11      continue
         pp=n*(z*p1-p2)/(z*z-real(1,id))
         z1=z
         z=z1-p1/pp
         if(abs(z-z1).gt.eps)go to 1
         x(i)=xm-xl*z
         x(n+1-i)=xm+xl*z
         w(i)=real(2,id)*xl/((real(1,id)-z*z)*pp*pp)
         w(n+1-i)=w(i)
 12   continue
      tmp = real(0,id)
      do i = 1, n
!         call mpwrite(6,x(i),w(i))
         tmp = tmp + cos(x(i)) * w(i)
      enddo
!      print*,'Integral from GAULEG:',(sin(x2)-sin(x1)),tmp
!      print*,'test of gauleg integration of cos(x)'
!      print*,'anaytical integration='
!      write(6,*)(sin(x2)-sin(x1))
!      print*,'numerical integration='
!      write(6,*)tmp
return
end

