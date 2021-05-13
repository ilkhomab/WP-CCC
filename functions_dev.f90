module flogs

  use precision_type
  implicit none

  real(kind=p16),allocatable,dimension(:)::fac,fac10,hat,faclog
  real(kind=p16),allocatable,dimension(:,:)::faccoef,lm_fac
  real(kind=p16),dimension(0:100)::fac2,pow2
  real(kind=id),dimension(0:100,0:100)::binom
  real(kind=id),dimension(0:15,-15:15,-15:15)::c7,c54
  real(kind=id),dimension(0:60,0:60,0:60)::g
  real(kind=id),dimension(0:15,0:20,0:20)::c21
!$acc declare create(fac10,faccoef,fac,hat,faclog,lm_fac) 

end module flogs

module data_base

  use precision_type
  implicit none


  real(kind=id)::scal,scal0,pmin,pmax,dmin,step_max,k_min,k_max,eps,unitarity2
  integer,allocatable,dimension(:,:,:)::nlm
  real(kind=id),allocatable,dimension(:,:,:)::pgrid
  real(kind=id),allocatable,dimension(:,:,:)::wf,wfres
  real(kind=id),allocatable,dimension(:,:)::enc
  real(kind=id),allocatable,dimension(:)::wbin,bin,dis,pt,wpt
  complex(kind=id),parameter::cu=cmplx(0._id,1._id,id)
  complex(kind=id),parameter::zero=cmplx(0._id,0._id,id)
  complex(kind=id),allocatable,dimension(:,:) :: zap
  integer::nbin,np,keyc,nneg,n1,neps,n2,key_bas,npt
  complex(kind=id),allocatable,dimension(:) :: ampb
  complex(kind=id),allocatable,dimension(:,:) :: amp

  integer,allocatable,dimension(:)::idir,fdir  
  integer,allocatable,dimension(:)::idirnl,fdirnl,larr,marr1,marr2  
  integer::jdir,jdirnl,wigarr


  logical :: probfile,wffile,ejenfile

  character*10:: target,projectile        ! target, projectile label

  integer::key,nstep,key_c,key_b,key_capt

  type rad_wf
    real(kind=id), allocatable,dimension(:)::rwf
    integer::jmin,jmax
  end type rad_wf

  type(rad_wf),allocatable,dimension(:)::st_t,st_p

  integer::nlnum_t,nlnum_p,imax,int_order,nk1,nk2,nk,npmax,lmar
  integer,allocatable,dimension(:)::nst_t,lst_t,nst_p,lst_p,lar,mar
  real(kind=id),allocatable,dimension(:)::k_t,k_p,alf_t,alf_p,enn_t,enn_p
  real(kind=id),allocatable,dimension(:)::offk,offw
  real(kind=p16),allocatable,dimension(:)::alf_t_p16,alf_p_p16
  integer::num_t,num_p,lmax_t,mmax_t,lmax_p,mmax_p
  integer,allocatable,dimension(:)::nps_t,nps_p,nmax_p,nmax_t
  integer,allocatable,dimension(:)::n_t,l_t,m_t,n_p,l_p,m_p,inl_t,inl_p,lff,mff,lii,mii
  integer,allocatable,dimension(:)::lob,upb

  type ci_coef
     real(kind=id),allocatable,dimension(:,:) ::  cl_t,cl_p
     real(kind=p16),allocatable,dimension(:,:) ::  cl_t_p16,cl_p_p16
  end type ci_coef

!!!Formfactor module

  real(kind=id), allocatable, dimension(:,:,:)::formf
  real(kind=id), allocatable, dimension(:,:,:)::formf_tr
  real(kind=id), allocatable, dimension(:,:,:)::formf_pr
  integer(kind=4),allocatable, dimension(:,:)::istart,istop
  real(kind=id), allocatable, dimension(:,:)::arylm
  real(kind=id), allocatable, dimension(:,:,:,:,:,:)::cleba
  integer(kind=4) :: ilowest


!!!


  integer::nzint,nomega
  real(kind=id) :: zintmax,dzint,omega_max,homega,dk1,dk2,projectile_charge,target_charge
  real(kind=id),allocatable,dimension(:) :: zint
  real(kind=id),allocatable,dimension(:) :: omega
  real(kind=id),allocatable,dimension(:,:,:) :: wig_d
  integer,allocatable,dimension(:,:,:) :: istar
  real(kind=id),allocatable,dimension(:,:) :: vmat_dir
  complex(kind=id),allocatable,dimension(:,:) :: vmat_ex
  complex(kind=id),allocatable,dimension(:,:) :: vmat_exres
  complex(kind=id),allocatable,dimension(:,:) :: vmat_obk
  real(kind=id),allocatable,dimension(:,:,:) :: vmat_dir_ar
  real(kind=id),allocatable,dimension(:,:,:) :: vmat_dir_om
  complex(kind=id),allocatable,dimension(:,:,:) :: vmat_ex_ar
  real(kind=id),allocatable,dimension(:,:,:) :: vmat_ex_om
  real(kind=id),allocatable,dimension(:,:,:) :: vmat_ex0_om
  complex(kind=id),allocatable,dimension(:,:,:) :: vmat_obk_ar

  type(ci_coef),allocatable,dimension(:)::c_t,c_p

  type energies
     real(kind=id),allocatable,dimension(:) ::  enl_t,enl_p
     real(kind=p16),allocatable,dimension(:) ::  enl_t_p16,enl_p_p16
  end type energies

  type(energies),allocatable,dimension(:)::en_t,en_p


  type formf_tt
     real(kind=id),allocatable,dimension(:) ::  ff_ttl
  end type formf_tt

  type formf_t
     type(formf_tt),allocatable,dimension(:) ::  ff_tt
  end type formf_t

  type(formf_t),allocatable,dimension(:)::ff_t



  type formf_pp
     real(kind=id),allocatable,dimension(:) ::  ff_ppl
  end type formf_pp

  type formf_p
     type(formf_pp),allocatable,dimension(:) ::  ff_pp
  end type formf_p

  type(formf_p),allocatable,dimension(:)::ff_p

  integer(kind=4)::ndouble,meshr,njdouble
  integer(kind=4),dimension(1:20)::jdouble
  integer(kind=4),parameter::maxr=100000
  real(kind=p16)::qmax,rmax
  real(kind=p16),dimension(1:20)::rjdouble
  real(kind=p16),dimension(1:maxr,3)::rmesh

  real(kind=id),allocatable,dimension(:,:,:)::radwf_t,radwf_p,radwf_t2,radwf_p2

  real(kind=p16),allocatable,dimension(:,:,:,:,:)::dir_ff_tt


  real(kind=id),allocatable,dimension(:)::u,glag,wglag,gleg,wgleg
  real(kind=id),allocatable,dimension(:)::step,b,wb

  real(kind=id) ::mp,mut,mup!mu=(mp+1)*mp/(2*mp+1)
  real(kind=id) :: a02!d-16         ! cm^2
  real(kind=id) :: ft,ryd,exp_fac,fmax !ft=2pi a02
  complex(kind=id) :: hi
  real(kind=id) :: ein,v,k0
!  real(kind=id) :: bmax
  real(kind=id),allocatable,dimension(:,:,:,:)::plma,plmb
  real(kind=id),allocatable,dimension(:,:,:)::jint
  real(kind=id),allocatable,dimension(:,:,:,:):: basis_a,basis_b,basis_bres
  integer(kind=4) :: nb,i0,nglag,ngleg,key_num,nbf,nbl!,num_dig

! Arrays needed for exchange matrix elements

  real(kind=p16),allocatable,dimension(:,:,:)::bnl
  real(kind=p16),allocatable,dimension(:,:)::gnl

  complex(kind=id),allocatable,dimension(:,:) :: vmatl,vright,vmid,vleft,vltmp,vrtmp,vover
  complex(kind=id),allocatable,dimension(:,:) :: v7m,v23,v7p
  real(kind=id) :: zmin,zmax,hinit,eps_min,eps_max,z,znew,h,zmid,bigr,cosang,ang,hh,h6
  complex(kind=id),allocatable,dimension(:) :: zam,zout,zdadz,dyt,dym,dymm,un,unn
  real(kind=id),allocatable,dimension(:,:)::rpow1,rpow2
  complex(kind=id),allocatable,dimension(:,:)::zij1,zij2

!$acc declare create(rpow1,rpow2) 

!$acc declare create(cleba,meshr,rmesh,nmax_t,nmax_p,lmax_t,lmax_p,mmax_t,mmax_p,jdirnl,fdirnl,idirnl,nst_t,lst_t,radwf_t,formf_tr,radwf_p,formf_pr,nst_p,lst_p)
!$acc declare create(nlnum_t,nlnum_p,njdouble,jdouble) 
!$acc declare create(num_t,num_p,k_t,k_p,enn_t,enn_p,inl_t,inl_p,fdir,idir,jdir,i0,step,nstep,l_t,l_p,m_t,m_p,v,projectile_charge,target_charge)
!$acc declare create(gleg,wgleg,glag,wglag)
!$acc declare create(ngleg,nglag,radwf_t2,radwf_p2)
!$acc declare create(n_t,n_p,en_t,en_p)
!$acc declare create(lar,mar,lmar,istar)
!$acc declare create(larr,marr1,marr2)

!$acc declare create(zam,zout,zdadz,dyt,dym,dymm,vleft,vmid,vright,vmatl,unitarity2,vltmp,vrtmp,un,unn)
!$acc declare create(vover,u,plma,plmb,jint,basis_a,basis_b,basis_bres,wig_d,z,znew,h,zmid,bigr,cosang,ang,hh,h6,zij1,zij2)

end module data_base

module acc_routines

  use precision_type,only:id
  implicit none

contains

function clebsch_gordan (j1, m1, j2, m2, j, m)
!$acc routine
  use flogs,only:fac10
! -------------------------------------------------------------------
! Calculates the Clebsch-Gordan coefficient <j1,m1,j2,m2|j,m>.
! Ref: D.M. Brink and G.R. Satchler, Angular Momentum, second edition,
!      Oxford University Press, p.34
! -------------------------------------------------------------------
! To avoid overflows for very large factorials, we use the
! external function fac10(n) = factorial(n)/10**n
! -------------------------------------------------------------------
!--------------------------------------------------------------------
! formal arguments
!--------------------------------------------------------------------
  integer, intent(in) :: j1, m1, j2, m2, j, m
!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------
  integer :: numin1, numin2, numin3, numin, numax1, numax2, numax3,numax, &
             nu, nuphase,jmin, jmax, a, b, c, alpha, beta,gamma1
  real(kind=id) :: clebsch_gordan, epsil, &
              Delta, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10,f11, &
              f12, f13, f14, f15, f16, p1, p2, sumnu, denom, term
!--------------------------------------------------------------------
! external functions
!--------------------------------------------------------------------
!  real(kind=id), external :: fac10
!--------------------------------------------------------------------
!     *****     P R O G R A M   S T A R T S   H E R E       ****
!--------------------------------------------------------------------
! check for incorrect input (magnetic quantum numbers)
!--------------------------------------------------------------------
!  if ( (iabs(m1) > j1) .or. (iabs(m2) > j2) .or. (iabs(m) > j) ) then
!    write (6, '(//A//)') ' input error: |m| larger than j'
!    return
!  endif
!--------------------------------------------------------------------
! test CG selection rules for magnetic quantum numbers
!    if (m1 + m2) .ne. m  --> CG = 0
!--------------------------------------------------------------------
  epsil = 1.e-3_id
  if ( iabs(m1+m2-m) > epsil ) then
     clebsch_gordan = real(0,id)
     return
  endif
!--------------------------------------------------------------------
! test triangle relations for angular momenta j1, j2, j
! if not fulfilled  -->  CG = 0
!--------------------------------------------------------------------
  jmin = iabs(j1 - j2)
  jmax = j1 + j2
  if (j < jmin .or. j > jmax) then
     clebsch_gordan = real(0,id)
     return
  endif
!--------------------------------------------------------------------
  a=j1; alpha=m1; b=j2; beta=m2; c=j; gamma1=m
!
  f1 = fac10(a+b-c)
  f2 = fac10(a+c-b)
  f3 = fac10(b+c-a)
  f4 = fac10(a+b+c+1)
  Delta = sqrt(f1*f2*f3/f4)
!
  p1 = sqrt(real(int(2 *c + 1),id))
!
  f5 = fac10(a+alpha)
  f6 = fac10(a-alpha)
  f7 = fac10(b+beta)
  f8 = fac10(b-beta)
  f9 = fac10(c+gamma1)
  f10= fac10(c-gamma1)
  p2 = sqrt(f5*f6*f7*f8*f9*f10)
!--------------------------------------------------------------------
! determine lower and upper limits for summation index nu; these
! are derived from the requirement that all factorials n! in the
! denominator are restricted to values with n >=0.
!--------------------------------------------------------------------
  numin1 = b-c-alpha
  numin2 = -c+a+beta
  numin3 = 0
  numin  = max(numin1,numin2,numin3)
!
  numax1 = a-alpha
  numax2 = b+beta
  numax3 = a+b-c
  numax = min(numax1,numax2,numax3)
!
  sumnu = real(0,id)
  do nu = numin, numax
     nuphase = (-1)**nu
!     if (mod(nu,2)==0)then
!        nuphase=1
!     else
!        nuphase=-1
!     endif
     f11 = fac10(a-alpha-nu)
     f12 = fac10(c-b+alpha+nu)
     f13 = fac10(b+beta-nu)
     f14 = fac10(c-a-beta+nu)
     f15 = fac10(nu)
     f16 = fac10(a+b-c-nu)
     denom = f11*f12*f13*f14*f15*f16
     term = nuphase / denom
     sumnu = sumnu + term
  end do
! factor sqrt(10) arises from replacing n! by fac10(n)=(n!/10**n)
  clebsch_gordan = Delta/sqrt(real(10,id)) *p1 *p2 * sumnu
!
  return
  end function clebsch_gordan

subroutine form0(l,fun,minfun,maxfun,f,g,temp)
!$acc routine seq
  use data_base,only:lmax_t,meshr
  implicit none
  integer, intent(in):: l
  real(kind=id), dimension(1:meshr), intent(in):: fun
  integer, intent(in):: minfun, maxfun
  real(kind=id), dimension(1:meshr), intent(out):: temp
  real(kind=id), dimension(0:meshr):: t1, t2
  real(kind=id),dimension(1:meshr,0:2*lmax_t)::f,g
  integer:: i
  integer:: if1, if2

! Set correct limits for integration: \int dr' fun(r') f(r')
  if1 = minfun
  if2 = maxfun



! do not initialise array temp to zero for faster calculation:
!  temp = 0.0


!  Find the integral for fun(r') * f(r') from 0 to istop. The function fun
!  already contains the Simson's integration weights.
! limits for which integral is required are:  if1:istop
         t1(if1 - 1) = 0._id
         do i=if1,if2
            t1(i) = t1(i-1) + fun(i)*f(i,l)
         end do

!  Find the integral of fun(r') * g(r') from infinity to if1
! limits for which integral should be taken:  if1:if2
         t2(if2) = 0.0
         do i=if2,if1,-1
            t2(i-1) = t2(i) + fun(i)*g(i,l)
         end do

!  Make the form factor by summing two parts

         temp(if1:if2) = t1(if1:if2)*g(if1:if2,l)

! if array tmp is initialised to zero then:
!         temp(it2min:istop) = temp(it2min:istop) +  t2(it2min:istop)*rpow_f(it2min:istop,l)
! if not:
         temp(if1:if2) = temp(if1:if2) +  t2(if1:if2)*f(if1:if2,l)





end subroutine form0

subroutine formAll(ni,li,nj,lj,temp)
!$acc routine seq 
  use data_base,only:lmax_t,meshr,rpow1,rpow2,nst_t,lst_t,fdirnl,idirnl,rmesh,radwf_t
  implicit none
  real(kind=id), dimension(1:meshr):: fun
  real(kind=id), dimension(0:(li+lj-iabs(li-lj))/2,1:meshr), intent(out):: temp
  real(kind=id), dimension(0:meshr):: t1, t2
  integer,intent(in):: ni,li,nj,lj
  integer:: i,l,ilmd

  do i=1,meshr
     fun(i)=radwf_t(ni,li,i)*radwf_t(nj,lj,i)*rmesh(i,1)*rmesh(i,1)*rmesh(i,3)
  enddo

ilmd=0
do l=iabs(li-lj),li+lj,2

  t1(0) = 0._id
  do i=1,meshr
     t1(i) = t1(i-1) + fun(i)*rpow1(i,l)
  end do

  t2(meshr) = 0._id
  do i=meshr,1,-1
     t2(i-1) = t2(i) + fun(i)*rpow2(i,l)
  end do

  temp(ilmd,1:meshr) = t1(1:meshr)*rpow2(1:meshr,l)

  temp(ilmd,1:meshr) = temp(ilmd,1:meshr) +  t2(1:meshr)*rpow1(1:meshr,l)

  ilmd=ilmd+1

  do i=1,meshr
    if(abs(temp(ilmd,i)).lt.1e-16_id)then
      temp(ilmd,i)=0._id  
    endif
  enddo

enddo





end subroutine formAll

subroutine formAllp(ni,li,nj,lj,temp)
!$acc routine 
  use data_base,only:lmax_t,meshr,rpow1,rpow2,nst_t,lst_t,fdirnl,idirnl,rmesh,radwf_p
  implicit none
  real(kind=id), dimension(1:meshr):: fun
  real(kind=id), dimension(0:(li+lj-iabs(li-lj))/2,1:meshr), intent(out):: temp
  real(kind=id), dimension(0:meshr):: t1, t2
  integer,intent(in):: ni,li,nj,lj
  integer:: i,l,ilmd

  do i=1,meshr
     fun(i)=radwf_p(ni,li,i)*radwf_p(nj,lj,i)*rmesh(i,1)*rmesh(i,1)*rmesh(i,3)
  enddo

ilmd=0
do l=iabs(li-lj),li+lj,2



! do not initialise array temp to zero for faster calculation:
!  temp = 0.0


!  Find the integral for fun(r') * f(r') from 0 to istop. The function fun
!  already contains the Simson's integration weights.
! limits for which integral is required are:  1:istop
  t1(0) = 0._id
  do i=1,meshr
     t1(i) = t1(i-1) + fun(i)*rpow1(i,l)
  end do

!  Find the integral of fun(r') * g(r') from infinity to 1 
! limits for which integral should be taken:  1:meshr
  t2(meshr) = 0._id
  do i=meshr,1,-1
     t2(i-1) = t2(i) + fun(i)*rpow2(i,l)
  end do

!  Make the form factor by summing two parts

  temp(ilmd,1:meshr) = t1(1:meshr)*rpow2(1:meshr,l)

! if array tmp is initialised to zero then:
!         temp(it2min:istop) = temp(it2min:istop) +  t2(it2min:istop)*rpow_f(it2min:istop,l)
! if not:
  temp(ilmd,1:meshr) = temp(ilmd,1:meshr) +  t2(1:meshr)*rpow1(1:meshr,l)

  ilmd=ilmd+1

  do i=1,meshr
    if(abs(temp(ilmd,i)).lt.1e-16_id)then
      temp(ilmd,i)=0._id  
    endif
  enddo

enddo





end subroutine formAllp


subroutine form01(l,fun,temp)
!$acc routine 
  use data_base,only:lmax_t,meshr,rpow1,rpow2
  implicit none
  real(kind=id), dimension(1:meshr), intent(in):: fun
  real(kind=id), dimension(1:meshr), intent(out):: temp
  real(kind=id), dimension(0:meshr):: t1, t2
  integer,intent(in):: l
  integer:: i
  integer:: if1, if2

! Set correct limits for integration: \int dr' fun(r') f(r')
  if1 = 1 
  if2 = meshr



! do not initialise array temp to zero for faster calculation:
!  temp = 0.0


!  Find the integral for fun(r') * f(r') from 0 to istop. The function fun
!  already contains the Simson's integration weights.
! limits for which integral is required are:  if1:istop
         t1(if1 - 1) = 0._id
         do i=if1,if2
            t1(i) = t1(i-1) + fun(i)*rpow1(i,l)
         end do

!  Find the integral of fun(r') * g(r') from infinity to if1
! limits for which integral should be taken:  if1:if2
         t2(if2) = 0._id
         do i=if2,if1,-1
            t2(i-1) = t2(i) + fun(i)*rpow2(i,l)
         end do

!  Make the form factor by summing two parts

         temp(if1:if2) = t1(if1:if2)*rpow2(if1:if2,l)

! if array tmp is initialised to zero then:
!         temp(it2min:istop) = temp(it2min:istop) +  t2(it2min:istop)*rpow_f(it2min:istop,l)
! if not:
         temp(if1:if2) = temp(if1:if2) +  t2(if1:if2)*rpow1(if1:if2,l)





end subroutine form01

function direct_num_t(ia,ib,costh,arg)
!$acc routine
use data_base,only:inl_t,l_t,m_t,cleba,v,projectile_charge,target_charge
use flogs,only:hat
implicit none
integer(kind=4)::ila,ilb,ia,ib,la,ma,lb,mb,myu,lambda,mini,maxi,ia0,ib0
real(kind=id)::rylm,costh,arg,sum,tmp,direct_num_t,funcn2_t

ila=inl_t(ia)
ilb=inl_t(ib)

la=l_t(ia)
ma=m_t(ia)
lb=l_t(ib)
mb=m_t(ib)
myu=ma-mb

sum=0._id
mini=iabs(la-lb)
maxi=la+lb

do lambda=mini,maxi,2

   if(iabs(myu).gt.lambda)cycle
   tmp=cleba(lambda,0,lb,0,la,0)*cleba(lambda,myu,lb,mb,la,ma)*funcn2_t(arg,ila,ilb,lambda)*rylm(lambda,myu,costh)
   sum=sum+tmp
enddo

direct_num_t=-hat(lb)/hat(la)*sum/v*projectile_charge
if(ia.eq.ib)direct_num_t=direct_num_t+1._id/arg/v*target_charge*projectile_charge   ! was added to add 1/R ! also +35 lines, target only



return
end

function rylm(l,m,x)
!$acc routine

  use flogs,only:faccoef
  implicit none
  integer::l,m,mabs
  real(kind=id):: plgndr,sp,x,rylm


!
! ****  spherical harmonic function (multipled by sqrt(4*pi/(2*l+1))) of order l,
! ****      magnetic projection m
! ****  with arguments theta (in radians). the azimuthal angle, phi
! ****  is zero. the spherical harmonic is a real function in this case.
!
  sp = 0._id 
!      x = cos(theta)
  mabs = iabs(m)
  sp = plgndr(l,mabs,x)
!      rl = dble(2*l+1)
  sp = sp*sqrt(faccoef(l,mabs))
!      sp = sp*dsqrt(factrl(l-mabs)/factrl(l+mabs))
  rylm = sp
  if(m.lt.0) then
     if(mod(mabs,2)/=0)rylm=-rylm
  endif
  return
end

function plgndr(l,mi,x)
!$acc routine

  use flogs,only:faclog,fac
  implicit none
  integer,intent(in)::l,mi
  integer::i,ll,m
  real(kind=id),intent(in):: x
  real(kind=id):: plgndr,pmm,fact,somx2,pll,pmmp1

!  mi=m
!  if(m.lt.0) then
!     m=iabs(m)
!  endif
  m=iabs(mi)

!  if(m.lt.0.or.m.gt.l.or.abs(x).gt.1)then
!     print*,l,m,x
!     stop 'bad arguments in plgndr'
!  endif
  pmm=1._id
  if(m.gt.0) then
    somx2=sqrt((1-x)*(1+x))
    fact=1._id
    do 11 i=1,m
      pmm=-pmm*fact*somx2
      fact=fact+2._id
11    continue
  endif
  if(l.eq.m) then
    plgndr=pmm
  else
    pmmp1=x*(2*m+1)*pmm
    if(l.eq.m+1) then
      plgndr=pmmp1
    else
      do 12 ll=m+2,l
        pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
        pmm=pmmp1
        pmmp1=pll
12      continue
      plgndr=pll
    endif
  endif

  if(mi.lt.0) then
    if(mod(m,2).eq.0)then
       plgndr=plgndr
    else
       plgndr=-plgndr
    endif

!    plgndr=plgndr*fac(l-m)/fac(l+m)
    plgndr=plgndr*exp(faclog(l-m)-faclog(l+m))

  endif


return
end

function funcn2_t(arg,ila,ilb,lambda)
!$acc routine
use data_base,only:rmesh,meshr,lst_t,nlnum_t,formf_tr
implicit none
real(kind=id),dimension(2):: ytmp,xtmp !2th order intrpl is used
integer::ila,ilb,lambda,ilmd,k,jres,la,lb
real(kind=id)::funcn2_t,arg,tmp


if(arg.gt.rmesh(meshr,1)) then
   funcn2_t=0._id
   goto 99
endif

la=lst_t(ila)
lb=lst_t(ilb)
ilmd=(lambda-iabs(la-lb))/2

k=(2*nlnum_t-ilb)*(ilb-1)/2+ila

if(arg.lt.rmesh(1,1)) then
   xtmp(1:2)= rmesh(1:2,1)
   ytmp(1:2)=formf_tr(k,ilmd,1:2)
   call intrpl2(xtmp,ytmp,arg,tmp)
   funcn2_t=tmp
   goto 99
endif

call locatexxr(arg,jres)
!print*,arg,jres
xtmp(1:2)=rmesh(jres:jres+1,1)
ytmp(1:2)=formf_tr(k,ilmd,jres:jres+1)
call intrpl2(xtmp,ytmp,arg,tmp)
funcn2_t=tmp

99 continue
return
end

subroutine intrpl2(xtmp,ytmp,arg,tmp)
!$acc routine
implicit none
real(kind=id),dimension(2):: ytmp,xtmp !2th order intrpl is used
real(kind=id)::denom,a,b,arg,tmp

denom=xtmp(1)-xtmp(2)
a=xtmp(1)*ytmp(2)-xtmp(2)*ytmp(1)
a=a/denom
b=ytmp(1)-ytmp(2)
b=b/denom
tmp=b*arg+a
return
end

subroutine locatexxr(x,jres)
!$acc routine
use data_base,only:rmesh,njdouble,jdouble
implicit none
real(kind=id)::dr,d,x
integer::jres,jd,jdp,idd

dr=rmesh(1,2)
d=1.442695_id*log(x/dr/32._id+1)
idd=int(1+d)
if(idd.gt.njdouble-1)idd=njdouble-1
jd=jdouble(idd)
jdp=int((x-rmesh(jd,1))/rmesh(jd+1,2))
jres=jd+jdp

return
end

function direct_num_p(ia,ib,costh,arg)
!$acc routine
use data_base,only:inl_p,l_p,m_p,cleba,v,projectile_charge,target_charge
use flogs,only:hat
implicit none
integer(kind=4)::ila,ilb,ia,ib,la,ma,lb,mb,myu,lambda,mini,maxi,ia0,ib0
real(kind=id)::costh,arg,sum,tmp,direct_num_p,funcn2_p

ila=inl_p(ia)
ilb=inl_p(ib)

la=l_p(ia)
ma=m_p(ia)
lb=l_p(ib)
mb=m_p(ib)
myu=ma-mb

sum=0._id
mini=iabs(la-lb)
maxi=la+lb

do lambda=mini,maxi,2

   if(iabs(myu).gt.lambda)cycle
   tmp=cleba(lambda,0,lb,0,la,0)*cleba(lambda,myu,lb,mb,la,ma)*funcn2_p(arg,ila,ilb,lambda)*rylm(lambda,myu,costh)
   sum=sum+tmp
enddo

direct_num_p=-hat(lb)/hat(la)*sum/v*target_charge
if(ia.eq.ib)direct_num_p=direct_num_p+1._id/arg/v*target_charge*projectile_charge   ! was added to add 1/R ! also -35 lines, target only



return
end



function funcn2_p(arg,ila,ilb,lambda)
!$acc routine
use data_base,only:rmesh,meshr,lst_p,nlnum_p,formf_pr
implicit none
real(kind=id),dimension(2):: ytmp,xtmp !2th order intrpl is used
integer::ila,ilb,lambda,ilmd,k,jres,la,lb
real(kind=id)::funcn2_p,arg,tmp

la=lst_p(ila)
lb=lst_p(ilb)
ilmd=(lambda-iabs(la-lb))/2

k=(2*nlnum_p-ilb)*(ilb-1)/2+ila

if(arg.gt.rmesh(meshr,1)) then
   funcn2_p=0._id
   goto 99
endif

if(arg.lt.rmesh(1,1)) then
   xtmp(1:2)= rmesh(1:2,1)
   ytmp(1:2)=formf_pr(k,ilmd,1:2)
   call intrpl2(xtmp,ytmp,arg,tmp)
   funcn2_p=tmp
   goto 99
endif

call locatexxr(arg,jres)

xtmp(1:2)=rmesh(jres:jres+1,1)
ytmp(1:2)=formf_pr(k,ilmd,jres:jres+1)
call intrpl2(xtmp,ytmp,arg,tmp)
funcn2_p=tmp

99 continue

return
end

subroutine get_plmijlm(i,j,lm,rho,plma,plmb)
!$acc routine
use data_base
!use ifport
implicit none
integer::i,j,l,m,lm
real(kind=id)::arg,rho,plma,plmb

l=lar(lm)
m=mar(lm)


arg=u(i)*gleg(j)+1._id
arg=arg/(u(i)+gleg(j))
plma=plgndr(l,m,arg)
arg=u(i)*gleg(j)-1._id
arg=arg/(u(i)-gleg(j))
plmb=plgndr(l,m,arg)

end subroutine get_plmijlm

subroutine wf_at_r(k_block,n,l,arg1,arg2,wf_a,wf_b,wf_b_res)
!$acc routine
  use precision_type
  use data_base
  implicit none
  integer,intent(in)::k_block,n,l
  real(kind=id),intent(in)::arg1,arg2
  real(kind=id),intent(out)::wf_a,wf_b,wf_b_res
!  real(kind=id)::wf_inter_t,wf_inter_t2,wf_inter_p,wf_inter_p2

  
if(k_block==0)then
  wf_a=wf_inter_p(n,l,arg1)
  wf_b=wf_inter_t(n,l,arg2)
  wf_b_res=wf_inter_t2(n,l,arg2)
else
  wf_a=wf_inter_t(n,l,arg1)
  wf_b=wf_inter_p(n,l,arg2)
  wf_b_res=wf_inter_p2(n,l,arg2)
endif

end subroutine wf_at_r

function wf_inter_p(n,l,arg)
!$acc routine
use data_base
implicit none
real(kind=id),dimension(2):: ytmp,xtmp !2th order intrpl is used
integer::n,l,jres
real(kind=id)::wf_inter_p,arg,tmp

if(arg.gt.rmesh(meshr,1)) then
   wf_inter_p=0._id
   goto 99
endif

if(arg.lt.rmesh(1,1)) then
   xtmp(1:2)= rmesh(1:2,1)
   ytmp(1:2)=radwf_p(n,l,1:2)
   call intrpl2(xtmp,ytmp,arg,tmp)
   wf_inter_p=tmp
   goto 99
endif

call locatexxr(arg,jres)

xtmp(1:2)=rmesh(jres:jres+1,1)
ytmp(1:2)=radwf_p(n,l,jres:jres+1)
call intrpl2(xtmp,ytmp,arg,tmp)
wf_inter_p=tmp

99 continue

return
end

function wf_inter_p2(n,l,arg)
!$acc routine
use data_base
implicit none
real(kind=id),dimension(2):: ytmp,xtmp !2th order intrpl is used
integer::n,l,jres
real(kind=id)::wf_inter_p2,arg,tmp

if(arg.gt.rmesh(meshr,1)) then
   wf_inter_p2=0._id
   goto 99
endif

if(arg.lt.rmesh(1,1)) then
   xtmp(1:2)= rmesh(1:2,1)
   ytmp(1:2)=radwf_p2(n,l,1:2)
   call intrpl2(xtmp,ytmp,arg,tmp)
   wf_inter_p2=tmp
   goto 99
endif

call locatexxr(arg,jres)

xtmp(1:2)=rmesh(jres:jres+1,1)
ytmp(1:2)=radwf_p2(n,l,jres:jres+1)
call intrpl2(xtmp,ytmp,arg,tmp)
wf_inter_p2=tmp

99 continue

return
end

function wf_inter_t(n,l,arg)
!$acc routine
use data_base
implicit none
real(kind=id),dimension(2):: ytmp,xtmp !2th order intrpl is used
integer::n,l,jres
real(kind=id)::wf_inter_t,arg,tmp

if(arg.gt.rmesh(meshr,1)) then
   wf_inter_t=0._id
   goto 99
endif

if(arg.lt.rmesh(1,1)) then
   xtmp(1:2)= rmesh(1:2,1)
   ytmp(1:2)=radwf_t(n,l,1:2)
   call intrpl2(xtmp,ytmp,arg,tmp)
   wf_inter_t=tmp
   goto 99
endif

call locatexxr(arg,jres)

xtmp(1:2)=rmesh(jres:jres+1,1)
ytmp(1:2)=radwf_t(n,l,jres:jres+1)
call intrpl2(xtmp,ytmp,arg,tmp)
wf_inter_t=tmp

99 continue

return
end

function wf_inter_t2(n,l,arg)
!$acc routine
use data_base
implicit none
real(kind=id),dimension(2):: ytmp,xtmp !2th order intrpl is used
integer::n,l,jres
real(kind=id)::wf_inter_t2,arg,tmp

if(arg.gt.rmesh(meshr,1)) then
   wf_inter_t2=0._id
   goto 99
endif

if(arg.lt.rmesh(1,1)) then
   xtmp(1:2)= rmesh(1:2,1)
   ytmp(1:2)=radwf_t2(n,l,1:2)
   call intrpl2(xtmp,ytmp,arg,tmp)
   wf_inter_t2=tmp
   goto 99
endif

call locatexxr(arg,jres)

xtmp(1:2)=rmesh(jres:jres+1,1)
ytmp(1:2)=radwf_t2(n,l,jres:jres+1)
call intrpl2(xtmp,ytmp,arg,tmp)
wf_inter_t2=tmp

99 continue

return
end

subroutine get_LR_mat_dev_pt2(il,jl,bigr,zz,vmatL,vmatR)
!$acc routine
  use data_base
  use flogs
  implicit none
  integer::mf,mi,lf,li,kf,ki,i,j,mmf,mmi
  integer,intent(in)::il,jl
  real(kind=id)::bigr,zz
  real(kind=id)::tmp,fact,func,com,comres,arg,fact_obk,com_obk,com_dir,var_lam,var_mu,carg,sarg,wig_d_f,wig_d_i
  complex(kind=id),intent(out) ::vmatL,vmatR
!  complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t) ::res,res_obk,resres,res_tmp,res_obk_tmp,resres_tmp
  complex(kind=id) :: sumz,sumresz,sumz_obk,sumzj,sumreszj,sumz_obkj,phase

  kf=n_t(il)
  lf=l_t(il)
  mf=m_t(il)
  ki=n_t(jl)
  li=l_t(jl)
  mi=m_t(jl)

 
  sumz=cmplx(0,0,id)
  sumresz=cmplx(0,0,id)
  sumz_obk=cmplx(0,0,id)

!!$acc loop independent reduction(+:sumz,sumresz,sumz_obk)
  do i=1,nglag
    var_lam=u(i)
    sumzj=cmplx(0,0,id)
    sumreszj=cmplx(0,0,id)
    sumz_obkj=cmplx(0,0,id)
!!$acc loop independent reduction(+:sumzj,sumreszj,sumz_obkj)
    do j=1,ngleg
      var_mu=gleg(j)
      arg=0.5_id*v*zz*var_lam*var_mu
      carg=cos(arg)
      sarg=sin(arg)
      func=var_lam-var_mu
      com_obk=func*plma(lf,mf,i,j)*plmb(li,mi,i,j)*jint(mi-mf,i,j)*basis_a(kf,lf,i,j)*basis_b(ki,li,i,j)
      com=com_obk*(var_lam+var_mu)
      comres=func*plma(lf,mf,i,j)*plmb(li,mi,i,j)*jint(mi-mf,i,j)*basis_a(kf,lf,i,j)*basis_bres(ki,li,i,j)*(var_lam+var_mu)
      sumzj=sumzj+com*cmplx(carg,sarg,id)*wgleg(j)
      sumreszj=sumreszj+comres*cmplx(carg,sarg,id)*wgleg(j)
      sumz_obkj=sumz_obkj+com_obk*cmplx(carg,sarg,id)*wgleg(j)
    enddo
    sumz=sumz+sumzj*wglag(i)
    sumresz=sumresz+sumreszj*wglag(i)
    sumz_obk=sumz_obk+sumz_obkj*wglag(i)
  enddo

  fact_obk=bigr*lm_fac(lf,mf)*lm_fac(li,mi)/8._id
  fact=fact_obk*bigr/2._id

  sumz=sumz*fact*cu**(mi-mf)
  sumresz=sumresz*fact*cu**(mi-mf)
  sumz_obk=sumz_obk*fact_obk*cu**(mi-mf)

  phase=cmplx(cos(zz*(k_t(inl_t(jl))-k_p(inl_t(il)))),sin(zz*(k_t(inl_t(jl))-k_p(inl_t(il)))),id)

  vmatL=sumz*phase
  vmatR=(sumz/bigr*target_charge*projectile_charge-sumz_obk*projectile_charge+sumresz-enn_t(inl_t(jl))*sumz)/v*phase

end subroutine get_LR_mat_dev_pt2

subroutine get_LR_mat_sum_pt(il,jl,bigr,zz,vmatL,vmatR)
!$acc routine
  use data_base
  use flogs
  implicit none
  integer::mf,mi,lf,li,kf,ki,i,j,mmf,mmi
  integer,intent(in)::il,jl
  real(kind=id)::bigr,zz,difc
  real(kind=id)::tmp,fact,func,com,comres,arg,fact_obk,com_obk,com_dir,var_lam,var_mu,carg,sarg,wig_d_f,wig_d_i
  complex(kind=id),intent(out) ::vmatL,vmatR
!  complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t) ::res,res_obk,resres,res_tmp,res_obk_tmp,resres_tmp
  complex(kind=id) :: sumz,sumresz,sumz_obk,sumzj,sumreszj,sumz_obkj,phase

  kf=n_t(il)
  lf=l_t(il)
  mf=m_t(il)
  ki=n_t(jl)
  li=l_t(jl)
  mi=m_t(jl)

!  difc=(projectile_charge-target_charge)/v
 
sumz=sum(zij1(1:nglag,1:ngleg)*plma(lf,mf,1:nglag,1:ngleg)*plmb(li,mi,1:nglag,1:ngleg)*jint(mi-mf,1:nglag,1:ngleg)*basis_a(kf,lf,1:nglag,1:ngleg)*basis_b(ki,li,1:nglag,1:ngleg))
sumresz=sum(zij1(1:nglag,1:ngleg)*plma(lf,mf,1:nglag,1:ngleg)*plmb(li,mi,1:nglag,1:ngleg)*jint(mi-mf,1:nglag,1:ngleg)*basis_a(kf,lf,1:nglag,1:ngleg)*basis_bres(ki,li,1:nglag,1:ngleg))
sumz_obk=sum(zij2(1:nglag,1:ngleg)*plma(lf,mf,1:nglag,1:ngleg)*plmb(li,mi,1:nglag,1:ngleg)*jint(mi-mf,1:nglag,1:ngleg)*basis_a(kf,lf,1:nglag,1:ngleg)*basis_b(ki,li,1:nglag,1:ngleg))

  fact_obk=bigr*lm_fac(lf,mf)*lm_fac(li,mi)/8._id
  fact=fact_obk*bigr/2._id

  sumz=sumz*fact*cu**(mi-mf)
  sumresz=sumresz*fact*cu**(mi-mf)
  sumz_obk=sumz_obk*fact_obk*cu**(mi-mf)

  phase=cmplx(cos(zz*(k_t(inl_t(jl))-k_p(inl_t(il)))),sin(zz*(k_t(inl_t(jl))-k_p(inl_t(il)))),id)

!  phase=phase*cmplx(cos(difc*log(v*(zz+bigr))),-sin(difc*log(v*(zz+bigr))))

  vmatL=sumz*phase
!  vmatR=(sumz/bigr*target_charge*projectile_charge-sumz_obk*projectile_charge)/v*phase  !was added to add 1/R
  vmatR=(sumz/bigr*target_charge*projectile_charge-sumz_obk*projectile_charge+sumresz-enn_t(inl_t(jl))*sumz)/v*phase  !was added to add 1/R
!  vmatR=(-sumz_obk*projectile_charge+sumresz-enn_t(inl_t(jl))*sumz)/v*phase

end subroutine get_LR_mat_sum_pt

subroutine get_LR_mat_sum_all_pt    !(bigr,zz)
  use data_base
  use flogs
  implicit none
  integer::mf,mi,lf,li,kf,ki,i,j,mmf,mmi
  integer::il,jl
  real(kind=id)::tmp,fact,func,com,comres,arg,fact_obk,com_obk,com_dir,var_lam,var_mu,carg,sarg
!  complex(kind=id),dimension(1:num_t,1:num_t),intent(out) ::vmatL,vmatR
!  complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t) ::res,res_obk,resres,res_tmp,res_obk_tmp,resres_tmp
  complex(kind=id) :: sumz,sumresz,sumz_obk,sumzj,sumreszj,sumz_obkj,phase
  complex(kind=id),dimension(1:nglag,1:ngleg) :: tmp1,tmp2,tmp3 

!$acc kernels 
do il=1,num_t
  do jl=1,num_t
     kf=n_t(il)
     lf=l_t(il)
     mf=m_t(il)
     ki=n_t(jl)
     li=l_t(jl)
     mi=m_t(jl)
    
     sumz=sum(zij1(1:nglag,1:ngleg)*plma(lf,mf,1:nglag,1:ngleg)*plmb(li,mi,1:nglag,1:ngleg)*jint(mi-mf,1:nglag,1:ngleg)*basis_a(kf,lf,1:nglag,1:ngleg)*basis_b(ki,li,1:nglag,1:ngleg))
     sumresz=sum(zij1(1:nglag,1:ngleg)*plma(lf,mf,1:nglag,1:ngleg)*plmb(li,mi,1:nglag,1:ngleg)*jint(mi-mf,1:nglag,1:ngleg)*basis_a(kf,lf,1:nglag,1:ngleg)*basis_bres(ki,li,1:nglag,1:ngleg))
     sumz_obk=sum(zij2(1:nglag,1:ngleg)*plma(lf,mf,1:nglag,1:ngleg)*plmb(li,mi,1:nglag,1:ngleg)*jint(mi-mf,1:nglag,1:ngleg)*basis_a(kf,lf,1:nglag,1:ngleg)*basis_b(ki,li,1:nglag,1:ngleg))

     fact_obk=bigr*lm_fac(lf,mf)*lm_fac(li,mi)/8._id
     fact=fact_obk*bigr/2._id
     phase=cmplx(cos(z*(k_t(inl_t(jl))-k_p(inl_t(il)))),sin(z*(k_t(inl_t(jl))-k_p(inl_t(il)))),id)

     sumz=sumz*fact*cu**(mi-mf)
     sumresz=sumresz*fact*cu**(mi-mf)
     sumz_obk=sumz_obk*fact_obk*cu**(mi-mf)

     vLtmp(il,jl)=sumz*phase
     vRtmp(il,jl)=(sumz/bigr*target_charge*projectile_charge-sumz_obk*projectile_charge+sumresz-enn_t(inl_t(jl))*sumz)/v*phase

  enddo
enddo
!$acc end kernels

end subroutine get_LR_mat_sum_all_pt


subroutine get_LR_mat_dev_pt(il,jl,bigr,zz,vmatL,vmatR)
!$acc routine
  use data_base
  use flogs
  implicit none
  integer::mf,mi,lf,li,kf,ki,i,j,mmf,mmi
  integer,intent(in)::il,jl
  real(kind=id)::bigr,zz
  real(kind=id)::tmp,fact,func,com,comres,arg,fact_obk,com_obk,com_dir,var_lam,var_mu,carg,sarg,wig_d_f,wig_d_i
  complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t),intent(out) ::vmatL,vmatR
  complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t) ::res,res_obk,resres,res_tmp,res_obk_tmp,resres_tmp
  complex(kind=id) :: sumz,sumresz,sumz_obk,phase

  kf=nst_t(il)
  lf=lst_t(il)
  ki=nst_t(jl)
  li=lst_t(jl)
  phase=cmplx(cos(zz*(k_t(jl)-k_p(il))),sin(zz*(k_t(jl)-k_p(il))),id)
 
  do mf=-min(lf,mmax_t),min(lf,mmax_t)
    do mi=-min(li,mmax_t),min(li,mmax_t)

       sumz=cmplx(0,0,id)
       sumresz=cmplx(0,0,id)
       sumz_obk=cmplx(0,0,id)

       fact_obk=bigr*lm_fac(lf,mf)*lm_fac(li,mi)/8._id
       fact=fact_obk*bigr/2._id
       do i=1,nglag
         var_lam=u(i)
         do j=1,ngleg
           var_mu=gleg(j)
           arg=0.5_id*v*zz*var_lam*var_mu
           carg=cos(arg)
           sarg=sin(arg)
           func=var_lam-var_mu
           com_obk=func*plma(lf,mf,i,j)*plmb(li,mi,i,j)*jint(mi-mf,i,j)*basis_a(kf,lf,i,j)*basis_b(ki,li,i,j)
           com=com_obk*(var_lam+var_mu)
           comres=func*plma(lf,mf,i,j)*plmb(li,mi,i,j)*jint(mi-mf,i,j)*basis_a(kf,lf,i,j)*basis_bres(ki,li,i,j)*(var_lam+var_mu)
           sumz=sumz+com*cmplx(carg,sarg,id)*wglag(i)*wgleg(j)
           sumresz=sumresz+comres*cmplx(carg,sarg,id)*wglag(i)*wgleg(j)
           sumz_obk=sumz_obk+com_obk*cmplx(carg,sarg,id)*wglag(i)*wgleg(j)
         enddo
       enddo

       res_tmp(mf,mi)=sumz*fact*cu**(mi-mf)
       resres_tmp(mf,mi)=sumresz*fact*cu**(mi-mf)
       res_obk_tmp(mf,mi)=sumz_obk*fact_obk*cu**(mi-mf)

    enddo
  enddo


  do mf=-min(lf,mmax_t),min(lf,mmax_t)
    do mi=-min(li,mmax_t),min(li,mmax_t)

       sumz=cmplx(0,0,id)
       sumresz=cmplx(0,0,id)
       sumz_obk=cmplx(0,0,id)

       do mmf=-min(lf,mmax_t),min(lf,mmax_t)
         wig_d_f=wig_d(lf,mmf,mf)
         do mmi=-min(li,mmax_t),min(li,mmax_t)
           wig_d_i=wig_d(li,mmi,mi)
           sumz=sumz+wig_d_f*wig_d_i*res_tmp(mmf,mmi)
           sumresz=sumresz+wig_d_f*wig_d_i*resres_tmp(mmf,mmi)
           sumz_obk=sumz_obk+wig_d_f*wig_d_i*res_obk_tmp(mmf,mmi)
         enddo
       enddo

       if(mod(lf+li+mf+mi,2)==0)then

          vmatL(mf,mi)=sumz*phase
          vmatR(mf,mi)=(sumz/bigr*target_charge*projectile_charge-sumz_obk*projectile_charge+sumresz-en_t(li)%enl_t(ki)*sumz)/v*phase

       else

          vmatL(mf,mi)=-sumz*phase
          vmatR(mf,mi)=-(sumz/bigr*target_charge*projectile_charge-sumz_obk*projectile_charge+sumresz-en_t(li)%enl_t(ki)*sumz)/v*phase

       endif

    enddo
  enddo

end subroutine get_LR_mat_dev_pt

subroutine get_LR_mat_dev_tp2(il,jl,bigr,zz,vmatL,vmatR)
!$acc routine
  use data_base
  use flogs
  implicit none
  integer::mf,mi,lf,li,kf,ki,i,j,mmf,mmi
  integer,intent(in)::il,jl
  real(kind=id)::bigr,zz
  real(kind=id)::tmp,fact,func,com,comres,arg,fact_obk,com_obk,com_dir,var_lam,var_mu,carg,sarg,wig_d_f,wig_d_i
  complex(kind=id),intent(out) ::vmatL,vmatR
!  complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t) ::res,res_obk,resres,res_tmp,res_obk_tmp,resres_tmp
  complex(kind=id) :: sumz,sumresz,sumz_obk,phase

  kf=n_t(il)
  lf=l_t(il)
  mf=m_t(il)
  ki=n_t(jl)
  li=l_t(jl)
  mi=m_t(jl)
  phase=cmplx(cos(zz*(k_p(inl_t(jl))-k_t(inl_t(il)))),sin(zz*(k_p(inl_t(jl))-k_t(inl_t(il)))),id)
 
  sumz=cmplx(0,0,id)
  sumresz=cmplx(0,0,id)
  sumz_obk=cmplx(0,0,id)

  fact_obk=bigr*lm_fac(lf,mf)*lm_fac(li,mi)/8._id
  fact=fact_obk*bigr/2._id
  do i=1,nglag
    var_lam=u(i)
    do j=1,ngleg
      var_mu=gleg(j)
      arg=0.5_id*v*zz*var_lam*var_mu
      carg=cos(arg)
      sarg=sin(arg)
      func=var_lam-var_mu
      com_obk=func*plma(lf,mf,i,j)*plmb(li,mi,i,j)*jint(mi-mf,i,j)*basis_a(kf,lf,i,j)*basis_b(ki,li,i,j)
      com=com_obk*(var_lam+var_mu)
      comres=func*plma(lf,mf,i,j)*plmb(li,mi,i,j)*jint(mi-mf,i,j)*basis_a(kf,lf,i,j)*basis_bres(ki,li,i,j)*(var_lam+var_mu)
      sumz=sumz+com*cmplx(carg,sarg,id)*wglag(i)*wgleg(j)
      sumresz=sumresz+comres*cmplx(carg,sarg,id)*wglag(i)*wgleg(j)
      sumz_obk=sumz_obk+com_obk*cmplx(carg,sarg,id)*wglag(i)*wgleg(j)
    enddo
  enddo

  sumz=sumz*fact*cu**(mi-mf)
  sumresz=sumresz*fact*cu**(mi-mf)
  sumz_obk=sumz_obk*fact_obk*cu**(mi-mf)



  vmatL=sumz*phase
  vmatR=(sumz/bigr*target_charge*projectile_charge-sumz_obk*target_charge+sumresz-enn_p(inl_p(jl))*sumz)/v*phase


end subroutine get_LR_mat_dev_tp2

subroutine get_LR_mat_sum_tp(il,jl,bigr,zz,vmatL,vmatR)
!$acc routine
  use data_base
  use flogs
  implicit none
  integer::mf,mi,lf,li,kf,ki,i,j,mmf,mmi
  integer,intent(in)::il,jl
  real(kind=id)::bigr,zz,difc
  real(kind=id)::tmp,fact,func,com,comres,arg,fact_obk,com_obk,com_dir,var_lam,var_mu,carg,sarg,wig_d_f,wig_d_i
  complex(kind=id),intent(out) ::vmatL,vmatR
!  complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t) ::res,res_obk,resres,res_tmp,res_obk_tmp,resres_tmp
  complex(kind=id) :: sumz,sumresz,sumz_obk,sumzj,sumreszj,sumz_obkj,phase

  kf=n_t(il)
  lf=l_t(il)
  mf=m_t(il)
  ki=n_t(jl)
  li=l_t(jl)
  mi=m_t(jl)

!  difc=(projectile_charge-target_charge)/v
 
sumz=sum(zij1(1:nglag,1:ngleg)*plma(lf,mf,1:nglag,1:ngleg)*plmb(li,mi,1:nglag,1:ngleg)*jint(mi-mf,1:nglag,1:ngleg)*basis_a(kf,lf,1:nglag,1:ngleg)*basis_b(ki,li,1:nglag,1:ngleg))
sumresz=sum(zij1(1:nglag,1:ngleg)*plma(lf,mf,1:nglag,1:ngleg)*plmb(li,mi,1:nglag,1:ngleg)*jint(mi-mf,1:nglag,1:ngleg)*basis_a(kf,lf,1:nglag,1:ngleg)*basis_bres(ki,li,1:nglag,1:ngleg))
sumz_obk=sum(zij2(1:nglag,1:ngleg)*plma(lf,mf,1:nglag,1:ngleg)*plmb(li,mi,1:nglag,1:ngleg)*jint(mi-mf,1:nglag,1:ngleg)*basis_a(kf,lf,1:nglag,1:ngleg)*basis_b(ki,li,1:nglag,1:ngleg))

  fact_obk=bigr*lm_fac(lf,mf)*lm_fac(li,mi)/8._id
  fact=fact_obk*bigr/2._id

  sumz=sumz*fact*cu**(mi-mf)
  sumresz=sumresz*fact*cu**(mi-mf)
  sumz_obk=sumz_obk*fact_obk*cu**(mi-mf)

  phase=cmplx(cos(zz*(k_p(inl_p(jl))-k_t(inl_t(il)))),sin(zz*(k_p(inl_p(jl))-k_t(inl_t(il)))),id)

!  phase=phase*cmplx(cos(difc*log(v*(zz+bigr))),-sin(difc*log(v*(zz+bigr))))

  vmatL=sumz*phase
  vmatR=(sumz/bigr*target_charge*projectile_charge-sumz_obk*target_charge+sumresz-enn_p(inl_p(jl))*sumz)/v*phase !was added to add 1/R 
!  vmatR=(-sumz_obk*target_charge+sumresz-enn_p(inl_p(jl))*sumz)/v*phase

end subroutine get_LR_mat_sum_tp


subroutine get_LR_mat_dev_tp(il,jl,bigr,zz,vmatL,vmatR)
!$acc routine
  use data_base
  use flogs
  implicit none
  integer::mf,mi,lf,li,kf,ki,i,j,mmf,mmi
  integer,intent(in)::il,jl
  real(kind=id)::bigr,zz
  real(kind=id)::tmp,fact,func,com,comres,arg,fact_obk,com_obk,com_dir,var_lam,var_mu,carg,sarg,wig_d_f,wig_d_i
  complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t),intent(out) ::vmatL,vmatR
  complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t) ::res,res_obk,resres,res_tmp,res_obk_tmp,resres_tmp
  complex(kind=id) :: sumz,sumresz,sumz_obk,phase

  kf=nst_t(il)
  lf=lst_t(il)
  ki=nst_t(jl)
  li=lst_t(jl)
  phase=cmplx(cos(zz*(k_p(jl)-k_t(il))),sin(zz*(k_p(jl)-k_t(il))),id)
 
  do mf=-min(lf,mmax_t),min(lf,mmax_t)
    do mi=-min(li,mmax_t),min(li,mmax_t)

       sumz=cmplx(0,0,id)
       sumresz=cmplx(0,0,id)
       sumz_obk=cmplx(0,0,id)

       fact_obk=bigr*lm_fac(lf,mf)*lm_fac(li,mi)/8._id
       fact=fact_obk*bigr/2._id
       do i=1,nglag
         var_lam=u(i)
         do j=1,ngleg
           var_mu=gleg(j)
           arg=0.5_id*v*zz*var_lam*var_mu
           carg=cos(arg)
           sarg=sin(arg)
           func=var_lam-var_mu
           com_obk=func*plma(lf,mf,i,j)*plmb(li,mi,i,j)*jint(mi-mf,i,j)*basis_a(kf,lf,i,j)*basis_b(ki,li,i,j)
           com=com_obk*(var_lam+var_mu)
           comres=func*plma(lf,mf,i,j)*plmb(li,mi,i,j)*jint(mi-mf,i,j)*basis_a(kf,lf,i,j)*basis_bres(ki,li,i,j)*(var_lam+var_mu)
           sumz=sumz+com*cmplx(carg,sarg,id)*wglag(i)*wgleg(j)
           sumresz=sumresz+comres*cmplx(carg,sarg,id)*wglag(i)*wgleg(j)
           sumz_obk=sumz_obk+com_obk*cmplx(carg,sarg,id)*wglag(i)*wgleg(j)
         enddo
       enddo

       res_tmp(mf,mi)=sumz*fact*cu**(mi-mf)
       resres_tmp(mf,mi)=sumresz*fact*cu**(mi-mf)
       res_obk_tmp(mf,mi)=sumz_obk*fact_obk*cu**(mi-mf)

    enddo
  enddo


  do mf=-min(lf,mmax_t),min(lf,mmax_t)
    do mi=-min(li,mmax_t),min(li,mmax_t)

       sumz=cmplx(0,0,id)
       sumresz=cmplx(0,0,id)
       sumz_obk=cmplx(0,0,id)

       do mmf=-min(lf,mmax_t),min(lf,mmax_t)
         wig_d_f=wig_d(lf,mmf,mf)
         do mmi=-min(li,mmax_t),min(li,mmax_t)
           wig_d_i=wig_d(li,mmi,mi)
           sumz=sumz+wig_d_f*wig_d_i*res_tmp(mmf,mmi)
           sumresz=sumresz+wig_d_f*wig_d_i*resres_tmp(mmf,mmi)
           sumz_obk=sumz_obk+wig_d_f*wig_d_i*res_obk_tmp(mmf,mmi)
         enddo
       enddo

       if(mod(lf+li+mf+mi,2)==0)then

          vmatL(mf,mi)=sumz*phase
          vmatR(mf,mi)=(sumz/bigr*target_charge*projectile_charge-sumz_obk*target_charge+sumresz-en_p(li)%enl_p(ki)*sumz)/v*phase

       else

          vmatL(mf,mi)=-sumz*phase
          vmatR(mf,mi)=-(sumz/bigr*target_charge*projectile_charge-sumz_obk*target_charge+sumresz-en_p(li)%enl_p(ki)*sumz)/v*phase

       endif

    enddo
  enddo

end subroutine get_LR_mat_dev_tp

subroutine rotation_pt(nf,lf,mf,ni,li,mi,zlm1,zrm1)
!$acc routine
use data_base,only:istar,mmax_t,wig_d,vLtmp,vRtmp,zero
complex(kind=id),intent(out)::zlm1,zrm1
complex(kind=id)::zlm2,zrm2
real(kind=id)::wigd1,wigd2
integer,intent(in)::nf,lf,mf,ni,li,mi
integer::m1,m2

    zlm1=zero
    zrm1=zero
    do m1=-min(lf,mmax_t),min(lf,mmax_t)
      wigd1=wig_d(lf,m1,mf)
      zlm2=zero
      zrm2=zero
      do m2=-min(li,mmax_t),min(li,mmax_t)
        wigd2=wig_d(li,m2,mi)
        zlm2=zlm2+wigd2*vLtmp(istar(nf,lf,m1),istar(ni,li,m2))
        zrm2=zrm2+wigd2*vRtmp(istar(nf,lf,m1),istar(ni,li,m2))
      enddo
      zlm1=zlm1+wigd1*zlm2
      zrm1=zrm1+wigd1*zrm2
    enddo
    zlm1=zlm1*(-1._id)**(lf+li+mf+mi)
    zrm1=zrm1*(-1._id)**(lf+li+mf+mi)

end subroutine rotation_pt

subroutine rotation_tp(nf,lf,mf,ni,li,mi,zlm1,zrm1)
!$acc routine
use data_base,only:istar,mmax_t,wig_d,vLtmp,vRtmp,zero
complex(kind=id),intent(out)::zlm1,zrm1
complex(kind=id)::zlm2,zrm2
real(kind=id)::wigd1,wigd2
integer,intent(in)::nf,lf,mf,ni,li,mi
integer::m1,m2

    zlm1=zero
    zrm1=zero
    do m1=-min(lf,mmax_t),min(lf,mmax_t)
      wigd1=wig_d(lf,m1,mf)
      zlm2=zero
      zrm2=zero
      do m2=-min(li,mmax_t),min(li,mmax_t)
        wigd2=wig_d(li,m2,mi)
        zlm2=zlm2+wigd2*vLtmp(istar(nf,lf,m1),istar(ni,li,m2))
        zrm2=zrm2+wigd2*vRtmp(istar(nf,lf,m1),istar(ni,li,m2))
      enddo
      zlm1=zlm1+wigd1*zlm2
      zrm1=zrm1+wigd1*zrm2
    enddo
    zlm1=zlm1*(-1._id)**(mf+mi)
    zrm1=zrm1*(-1._id)**(mf+mi)

end subroutine rotation_tp


end module acc_routines

subroutine readin
  use precision_type
  use data_base
  use flogs
  use mpi,only:myid

  implicit none
  logical:: ex,offshell
  integer:: iostat_data,i,ii,j,n,l,m,jmax,k,h_pow,izmax,izmin,i_eps_min,i_eps_max,lf,li,mf,mi,nchunk,i_on,inl,iread,io,m1,m2
  integer(kind=4),parameter :: nfile=10
  real*8::ein_d,bmin,bmax,bmax_d,tmp,projectile_mass,target_mass
  real*8,allocatable,dimension(:)::dalf_t,dalf_p
  real(kind=id),allocatable,dimension(:) :: enrg,wfeig
  real(kind=id) :: k2,r,f,tmp_mp,arg
  real(kind=id) :: fd,difs,dmax,xtmp
  real(kind=id) :: direct_num,dd
  REAL(KIND=ID),ALLOCATABLE::TMP_STP(:)

  print*,'Entering readin'


!  num_dig=50

  mppic=3.1415926535897932384626433832795028841971693993751_id
  mppic16=3.1415926535897932384626433832795028841971693993751_p16
  mp=1836.15267261_id
  a02=0.2800285204_id; ft=2*mppic*a02;ryd=13.6058_id
  hi=cmplx(0,0.5_id,id)

inquire(file='prob',exist=probfile)
inquire(file='ejenfile',exist=ejenfile)

  inquire(file='data.in',exist=ex)
  if(.not. ex) then
     if(myid==0)then
     print*,'file data.in does not exists. stop'
     endif
     stop
  end if
  open(nfile,file='data.in',iostat=iostat_data)
  if(iostat_data.ne.0) then
     print*, '********  file  canot open file data.in'
     stop
  end if

!  read(nfile,'(A10)') projectile
!    write(*,*) 'projectile: ', projectile

!  read(nfile,'(A10)') target
!    write(*,*) 'target: ', target

  read(nfile,*)projectile_charge,projectile_mass
  read(nfile,*)target_charge,target_mass

if(myid==0)then
  if(projectile_charge==-1d0)then
    write(*,*) 'projectile: antiproton'
  else
    write(*,*) 'projectile: proton'
  endif
endif

if(myid==0)then
  if(target_charge==1d0)then
    write(*,*) 'target: Hydrogen'
  else
    write(*,*) 'target: He+'
  endif
endif

  mut=(target_mass*mp+1)*projectile_mass*mp/((projectile_mass+target_mass)*mp+1)
  mup=(projectile_mass*mp+1)*target_mass*mp/((projectile_mass+target_mass)*mp+1)

  read(nfile,*)key
  read(nfile,*)key_capt

if(myid==0)then
  if(key==0)then
      write(*,*)'First Born Approximation'
  else
      write(*,*)'Full calculations'
  endif
endif



  read(nfile,*) ein,i0
if(myid==0)then
  print*,'energy',ein
endif

if(myid==0)then
    write(*,'(" Ein: ",F10.1," KeV")') ein
    write(*,'(" Initial state is: ",I3)') i0
endif

  v=Sqrt(ein*1000._id/(ryd*projectile_mass*mp))
  k0=mut*v
if(myid==0)then
  write(unit=6,fmt='(" Velocity:",F10.4," au")') v
  write(unit=6,fmt='(" Momentum:",F10.4," au")') k0
endif

  read(nfile,*) qmax,ndouble,rmax

  call grids(qmax,ndouble,rmax,rmesh,maxr,meshr,jdouble,ndouble+2,njdouble)
!$acc update device(meshr,rmesh)

  do i=1,njdouble
    rjdouble(i)=rmesh(jdouble(i),1)
if(myid==0)then
    print*,rjdouble(i)
endif
  enddo

  read(nfile,*)key_bas
if(myid==0)then
  print*,'hihhhhhhhhh'
endif

if(key_bas==0)then

  allocate(nmax_t(0:0))

if(myid==0)then
  write(unit=6,fmt=*) 'Wave-packet pseudostates have been used'
endif

  open(11,file='WPS',iostat=iostat_data)

  read(11,*)nmax_t(0),lmax_t,mmax_t



!$acc update device(nmax_t,lmax_t,mmax_t)
if(myid==0)then
    write(unit=6,fmt='(" Lmax_t=",I2," Mmax_t=",I2)') lmax_t,mmax_t
endif
!  read(unit=nfile,fmt=*) nbin,np,scal
  read(unit=11,fmt=*) k_min,k_max
  read(unit=11,fmt=*) np,n1
  read(unit=11,fmt=*) neps,eps,n2
  nbin=n1+neps+n2
  allocate(radwf_t(1:nmax_t(0),0:lmax_t,1:meshr))
  allocate(radwf_t2(1:nmax_t(0),0:lmax_t,1:meshr))
  allocate(radwf_p(1:nmax_t(0),0:lmax_t,1:meshr))
  allocate(radwf_p2(1:nmax_t(0),0:lmax_t,1:meshr))


  allocate(en_t(0:lmax_t))
  allocate(en_p(0:lmax_t))
  do l=0,lmax_t
    allocate(en_t(l)%enl_t(1:nmax_t(0)))
    allocate(en_p(l)%enl_p(1:nmax_t(0)))
  enddo

  allocate(wfeig(1:meshr))


  nneg=nmax_t(0)-nbin


  call get_faclog
!$acc update device(fac10,faccoef,fac,hat,faclog,lm_fac)
if(myid==0)then
      print*,'hiha'
endif



      do l=0,lmax_t
        do n=1,nneg-l
          call get_eigen(target_charge,n+l,l,wfeig(1:meshr))
          radwf_t(n,l,1:meshr)=wfeig(1:meshr)!*rmesh(1:meshr,1)
          en_t(l)%enl_t(n)=-target_charge*target_charge/2d0/(n+l)**2
          en_p(l)%enl_p(n)=-target_charge*target_charge/2d0/(n+l)**2
          radwf_t2(n,l,1:meshr)=en_t(l)%enl_t(n)*wfeig(1:meshr)!*rmesh(1:meshr,1)
          if(projectile_charge.gt.0d0.and.key_capt==0)then
            call get_eigen(projectile_charge,n+l,l,wfeig(1:meshr))
            radwf_p(n,l,1:meshr)=wfeig(1:meshr)!*rmesh(1:meshr,1)
            en_p(l)%enl_p(n)=-projectile_charge*projectile_charge/2d0/(n+l)**2
            radwf_p2(n,l,1:meshr)=en_p(l)%enl_p(n)*wfeig(1:meshr)!*rmesh(1:meshr,1)
          endif
        enddo
      enddo

  allocate(dis(1:nbin))
  allocate(wf(nbin,0:lmax_t,1:meshr))
  allocate(wfres(nbin,0:lmax_t,1:meshr))
  allocate(enc(nbin,0:lmax_t))

if(myid==0)then
  print*,'eigen done'
endif

  npmax=10000
  allocate(bin(1:nbin+1))
  allocate(wbin(1:nbin+1))
  allocate(pgrid(nbin,npmax,1:2))
  call get_cwf(target_charge)



      do l=0,lmax_t
        do n=nneg-l+1,nmax_t(0)-l
          if(.not.probfile)then
            radwf_t(n,l,1:meshr)=wf(n-nneg+l,l,1:meshr)/rmesh(1:meshr,1)
            radwf_t2(n,l,1:meshr)=wfres(n-nneg+l,l,1:meshr)/rmesh(1:meshr,1)
          endif
          en_t(l)%enl_t(n)=enc(n-nneg+l,l)/2d0
          en_p(l)%enl_p(n)=enc(n-nneg+l,l)/2d0
        enddo
      enddo


if(myid==0)then
  print*,'target + done'
endif
  if(projectile_charge.gt.0d0.and.key_capt==0)then
          call get_cwf(projectile_charge)

        do l=0,lmax_t
          do n=nneg-l+1,nmax_t(0)-l
            if(.not.probfile)then
              radwf_p(n,l,1:meshr)=wf(n-nneg+l,l,1:meshr)/rmesh(1:meshr,1)
              radwf_p2(n,l,1:meshr)=wfres(n-nneg+l,l,1:meshr)/rmesh(1:meshr,1)
            endif
            en_p(l)%enl_p(n)=enc(n-nneg+l,l)/2d0
          enddo
        enddo

        do i=1,meshr


        enddo

if(myid==0)then
    print*,'projectile + done'
endif

  endif


  num_t=0
  nlnum_t=0
  do l=0,lmax_t
     nlnum_t=nlnum_t+nmax_t(0)-l
     m=min(2*l+1,2*mmax_t+1)
     num_t=num_t+(nmax_t(0)-l)*m
  enddo
  num_p=num_t
  nlnum_p=nlnum_t

jdir=num_t*(num_t-1)/2

allocate(fdir(1:jdir))
allocate(idir(1:jdir))


j=0
do f=1,num_t
  do i=1,f-1
    j=j+1
    fdir(j)=f
    idir(j)=i
  enddo
enddo

wigarr=0
do l=0,lmax_t
  do m1=-min(l,mmax_t),min(l,mmax_t)
    do m2=-min(l,mmax_t),min(l,mmax_t)
      wigarr=wigarr+1
    enddo
  enddo
enddo

allocate(larr(1:wigarr))
allocate(marr1(1:wigarr))
allocate(marr2(1:wigarr))
wigarr=0
do l=0,lmax_t
  do m1=-min(l,mmax_t),min(l,mmax_t)
    do m2=-min(l,mmax_t),min(l,mmax_t)
      wigarr=wigarr+1
      larr(wigarr)=l;marr1(wigarr)=m1;marr2(wigarr)=m2
    enddo
  enddo
enddo

!$acc update device(larr,marr1,marr2)


lmar=0
do l=0,lmax_t
  do m1=-min(l,mmax_t),min(l,mmax_t)
      lmar=lmar+1
  enddo
enddo

allocate(lar(1:lmar))
allocate(mar(1:lmar))

lmar=0
do l=0,lmax_t
  do m1=-min(l,mmax_t),min(l,mmax_t)
      lmar=lmar+1
      lar(lmar)=l;mar(lmar)=m1
  enddo
enddo

  allocate(st_t(1:nlnum_t),nst_t(1:nlnum_t),lst_t(1:nlnum_t))
  allocate(k_t(1:nlnum_t))
  allocate(enn_t(1:nlnum_t))
  allocate(n_t(1:num_t),l_t(1:num_t),m_t(1:num_t),inl_t(1:num_t))
  allocate(istar(1:nmax_t(0),0:lmax_t,-mmax_t:mmax_t))


  allocate(st_p(1:nlnum_p),nst_p(1:nlnum_p),lst_p(1:nlnum_p))
  allocate(k_p(1:nlnum_p))
  allocate(enn_p(1:nlnum_p))
  allocate(lob(1:nlnum_p))
  allocate(upb(1:nlnum_p))
  allocate(n_p(1:num_p),l_p(1:num_p),m_p(1:num_p),inl_p(1:num_p))

if(myid==0)then
    write(unit=6,fmt=*) nlnum_t, " Target state energies"
    write(unit=6,fmt=*) "      Exact     n=1       n=2    ..."
endif

  i=0;j=0
  do l=0, lmax_t, 1
     do n=1, nmax_t(0)-l, 1
        j=j+1

        k2=k0*k0-2*mut*(en_t(l)%enl_t(n)-en_t(0)%enl_t(1))
        if( k2.ge.0d0 ) then
           k_t(j)=sqrt(k2)
        else
           k_t(j)=0._id
        end if
        enn_t(j)=en_t(l)%enl_t(n)
!        k2=k0*k0*mup/mut-2*mup*(en_p(l)%enl_p(n)-en_t(0)%enl_t(1))
        k2=k0*k0-2*mut*(en_p(l)%enl_p(n)-en_t(0)%enl_t(1))
        if( k2.ge.0d0 ) then
           k_p(j)=sqrt(k2)
        else
           k_p(j)=0._id
        end if
        enn_p(j)=en_p(l)%enl_p(n)

        nst_t(j)=n
        lst_t(j)=l

        lob(j)=i+1
        do m=-min(l,mmax_t), min(l,mmax_t), 1
           i=i+1
           n_t(i)=n;   l_t(i)=l;   m_t(i)=m
           inl_t(i)=j
           istar(n,l,m)=i
        end do
        upb(j)=i
     end do
if(myid==0)then
        write(unit=6,fmt='( I3,20000F10.4)') l,(dble(en_t(l)%enl_t(n)),n=1,nmax_t(0)-l)
        write(unit=6,fmt='( I3,20000F10.4)') l,(dble(en_p(l)%enl_p(n)),n=1,nmax_t(0)-l)
endif
  end do

  st_p=st_t;nst_p=nst_t;lst_p=lst_t;n_p=n_t;l_p=l_t;m_p=m_t;inl_p=inl_t


  if(i/=num_t.or.j/=nlnum_t)then
    print*,i,num_t,j,nlnum_t
    stop'problem in target basis'
  endif

if(myid==0)then
    write(unit=6,fmt=*) "Altogether, nlm target+projectile states are",num_t+num_p
endif


endif   !   if(key_bas==0)then



if(key_bas==1)then

if(myid==0)then
  write(unit=6,fmt=*) 'Laguerre pseudostates have been used'
endif

!  read(nfile,*)
!  read(nfile,*)
!  read(nfile,*)
!  read(nfile,*)

  open(11,file='LPS',iostat=iostat_data)

  read(11,*)lmax_t,mmax_t
if(myid==0)then
    write(unit=6,fmt='(" Lmax_t=",I2," Mmax_t=",I2)') lmax_t,mmax_t
endif
  allocate(alf_t(0:lmax_t))
  allocate(alf_t_p16(0:lmax_t))
  allocate(dalf_t(0:lmax_t))
  allocate(nmax_t(0:lmax_t))
  allocate(nps_t(0:lmax_t))
  read(unit=11,fmt=*) nmax_t(0:lmax_t)

!  allocate(radwf_t(1:nmax_t(0),0:lmax_t,1:meshr))
  read(unit=11,fmt=*) nps_t(0:lmax_t)
  allocate(radwf_t(1:nps_t(0),0:lmax_t,1:meshr))
  allocate(radwf_t2(1:nps_t(0),0:lmax_t,1:meshr))
if(myid==0)then
  write(unit=6,fmt='(" nps_t(l)=",20I4)') nps_t(0:lmax_t)
endif
  read(unit=11,fmt=*) alf_t_p16(0:lmax_t)
if(myid==0)then
  write(unit=6,fmt='(" alf_t(l)=",20es14.6)') (alf_t_p16(l),l=0,lmax_t)
  print*,'alf_t',alf_t_p16
endif

  do l=0,lmax_t
    alf_t(l)=real(alf_t_p16(l),id)
  enddo

!  allocate(ar_exp_integ(0:2*nps_t(0)+2*lmax_t+4,0:2*lmax_t))

  allocate(c_t(0:lmax_t))
  allocate(en_t(0:lmax_t))
  do l=0,lmax_t
if(myid==0)then
    print*,'l',l,alf_t(l)
endif
    allocate(c_t(l)%cl_t(1:nps_t(l),1:nps_t(l)))
    allocate(en_t(l)%enl_t(1:nps_t(l)))
    allocate(c_t(l)%cl_t_p16(1:nps_t(l),1:nps_t(l)))
    allocate(en_t(l)%enl_t_p16(1:nps_t(l)))
!    call get_ci_coef_and_energy(nps_t(l),l,real(alf_t(l),id),c_t(l)%cl_t(1:nps_t(l),1:nps_t(l)),en_t(l)%enl_t(1:nps_t(l)))
    call get_ci_coef_and_energy_p16(target_charge,nps_t(l),l,alf_t_p16(l),c_t(l)%cl_t_p16(1:nps_t(l),1:nps_t(l)),en_t(l)%enl_t_p16(1:nps_t(l)),radwf_t(1:nps_t(l),l,1:meshr),radwf_t2(1:nps_t(l),l,1:meshr))
    c_t(l)%cl_t(1:nps_t(l),1:nps_t(l))=real(c_t(l)%cl_t_p16(1:nps_t(l),1:nps_t(l)),id)
    en_t(l)%enl_t(1:nps_t(l))=real(en_t(l)%enl_t_p16(1:nps_t(l)),id)
if(myid==0)then
    write(6,fmt='(I4,20es14.6)')l,(en_t(l)%enl_t(i),i=1,nps_t(l))
endif
  enddo

!  write(6,*)sum(radwf_t(13,2,1:meshr)*radwf_t2(13,2,1:meshr)*rmesh(1:meshr,1)*rmesh(1:meshr,1)*rmesh(1:meshr,3))
!  print*,'boldi'
!  stop

  if(projectile_charge.gt.0d0.and.key_capt==0)then

    allocate(radwf_p(1:nps_t(0),0:lmax_t,1:meshr))
    allocate(radwf_p2(1:nps_t(0),0:lmax_t,1:meshr))
    allocate(c_p(0:lmax_t))
    allocate(en_p(0:lmax_t))
    do l=0,lmax_t
if(myid==0)then
      print*,'l',l,alf_t(l)
endif
      allocate(c_p(l)%cl_p(1:nps_t(l),1:nps_t(l)))
      allocate(en_p(l)%enl_p(1:nps_t(l)))
      allocate(c_p(l)%cl_p_p16(1:nps_t(l),1:nps_t(l)))
      allocate(en_p(l)%enl_p_p16(1:nps_t(l)))
!      call get_ci_coef_and_energy(nps_t(l),l,real(alf_t(l),id),c_t(l)%cl_t(1:nps_t(l),1:nps_t(l)),en_t(l)%enl_t(1:nps_t(l)))
      call get_ci_coef_and_energy_p16(projectile_charge,nps_t(l),l,alf_t_p16(l),c_p(l)%cl_p_p16(1:nps_t(l),1:nps_t(l)),en_p(l)%enl_p_p16(1:nps_t(l)),radwf_p(1:nps_t(l),l,1:meshr),radwf_p2(1:nps_t(l),l,1:meshr))
      c_p(l)%cl_p(1:nps_t(l),1:nps_t(l))=real(c_p(l)%cl_p_p16(1:nps_t(l),1:nps_t(l)),id)
      en_p(l)%enl_p(1:nps_t(l))=real(en_p(l)%enl_p_p16(1:nps_t(l)),id)
if(myid==0)then
      write(6,fmt='(I4,20es14.6)')l,(en_p(l)%enl_p(i),i=1,nps_t(l))
endif
    enddo

  do i=1,meshr
if(myid==0)then
    print*,radwf_p(1,0,i)
endif
  enddo

  endif

  num_t=0
  nlnum_t=0
  do l=0,lmax_t
     nlnum_t=nlnum_t+nmax_t(l)
     m=min(2*l+1,2*mmax_t+1)
     num_t=num_t+(nmax_t(l))*m
  enddo
  num_p=num_t
  nlnum_p=nlnum_t

  allocate(st_t(1:nlnum_t),nst_t(1:nlnum_t),lst_t(1:nlnum_t))
  allocate(k_t(1:nlnum_t))
  allocate(n_t(1:num_t),l_t(1:num_t),m_t(1:num_t),inl_t(1:num_t))
  allocate(istar(1:nmax_t(0),0:lmax_t,-mmax_t:mmax_t))


  allocate(st_p(1:nlnum_p),nst_p(1:nlnum_p),lst_p(1:nlnum_p))
  allocate(k_p(1:nlnum_p))
  allocate(n_p(1:num_p),l_p(1:num_p),m_p(1:num_p),inl_p(1:num_p))

if(myid==0)then
    write(unit=6,fmt=*) nlnum_t, " Target state energies"
    write(unit=6,fmt=*) "      Exact     n=1       n=2    ..."
endif



  i=0;j=0
  do l=0, lmax_t, 1
     do n=1, nmax_t(l), 1
        j=j+1

        k2=k0*k0-2*mut*(en_t(l)%enl_t(n)-en_t(0)%enl_t(1))
        if( k2.ge.0d0 ) then
           k_t(j)=sqrt(k2)
        else
           k_t(j)=0._id
        end if
  if(projectile_charge.gt.0d0.and.key_capt==0)then
!        k2=k0*k0*mup/mut-2*mup*(en_p(l)%enl_p(n)-en_t(0)%enl_t(1))
        k2=k0*k0-2*mut*(en_p(l)%enl_p(n)-en_t(0)%enl_t(1))
        if( k2.ge.0d0 ) then
           k_p(j)=sqrt(k2)
        else
           k_p(j)=0._id
        end if
  endif

        nst_t(j)=n
        lst_t(j)=l

        do m=-min(l,mmax_t), min(l,mmax_t), 1
           i=i+1
           n_t(i)=n;   l_t(i)=l;   m_t(i)=m
           inl_t(i)=j
           istar(n,l,m)=i
        end do
     end do


if(myid==0)then
        write(unit=6,fmt='( I3,200F10.4)') l,(dble(en_t(l)%enl_t(n)),n=1,nmax_t(0)-l)
      if(projectile_charge.gt.0d0.and.key_capt==0)then
        write(unit=6,fmt='( I3,200F10.4)') l,(dble(en_p(l)%enl_p(n)),n=1,nmax_t(0)-l)
      endif
endif
  end do

  st_p=st_t;nst_p=nst_t;lst_p=lst_t;n_p=n_t;l_p=l_t;m_p=m_t;inl_p=inl_t


endif   !   if(key_bas==0)then

!  do ii=1,num_t

!    write(6,fmt='(4I4)')ii,n_t(ii),l_t(ii),m_t(ii)
!    print*,ii,inl_t(ii),dble(k_t(inl_t(ii)))

!  enddo

!  stop

  allocate(ff_t(1:nlnum_t))

  do i=1,nlnum_t

    if(key.eq.0)then
      jmax=1
    else
      jmax=i
    endif

    allocate(ff_t(i)%ff_tt(1:jmax))

    do j=1,jmax

      allocate(ff_t(i)%ff_tt(j)%ff_ttl(iabs(lst_t(i)-lst_t(j)):lst_t(i)+lst_t(j)))

    enddo

  enddo

    allocate(arylm(0:2*lmax_t,-2*lmax_t:2*lmax_t))
    allocate(formf((nlnum_t+1)*nlnum_t/2,0:lmax_t,1:meshr))
    allocate(formf_tr((nlnum_t+1)*nlnum_t/2,0:lmax_t,1:meshr))
    allocate(formf_pr((nlnum_t+1)*nlnum_t/2,0:lmax_t,1:meshr))


jdirnl=(nlnum_t+1)*nlnum_t/2

allocate(fdirnl(1:jdirnl))
allocate(idirnl(1:jdirnl))

!            k=(j-1)*(2*nlnum_t-j)/2+i

do f=1,nlnum_t
  do i=1,f
    idirnl(int((i-1)*(2*nlnum_t-i)/2+f))=i
    fdirnl(int((i-1)*(2*nlnum_t-i)/2+f))=f
  enddo
enddo


!$acc update device(jdirnl,fdirnl,idirnl,nst_t,lst_t,radwf_t,nst_p,lst_p,radwf_p)
!$acc update device(nmax_p,lmax_p,mmax_p)

!if(j/=jdirnl)then
!  print*,j,jdirnl
!  stop 'oooooooooooooooooopppppppppppppppppssssssssssssssss'
!endif


    call precalc

    if(.not.probfile)then
      call formfactor_t
      if(key_capt==0)then
        call formfactor_p(formf_pr)
      endif
if(myid==0)then
      print*,'form factors have been saved'
endif
    endif



!    do i=1,meshr
!
!      write(694,fmt='(2es14.6)')rmesh(i,1),radwf_t(10,0,i)
!
!    enddo
!!
!    stop


  if(probfile)then

      read(nfile,*)
!      open(701,file='input.dat')
!      read(701,*)bmin,bmax,nb

!      bmax=real(bmax_d,id)

!open(70,file='points.dat')
open(70,file='prob')
nb=0
io=0
do while(io==0)
  read(70,*,iostat=io)xtmp
  nb=nb+1
enddo
close(70)
nb=nb-1
print*, 'nb',nb


      allocate(b(1:nb),wb(1:nb))

      open(70,file='points.dat')
      open(71,file='weights.dat')

      do iread=1,nb
        read(70,*)b(iread)
        read(71,*)wb(iread)
      enddo
      close(70) 
      close(71) 

  else

      read(unit=nfile,fmt=*)nbf,nbl

      nb=nbl-nbf+1

!      bmax=real(bmax_d,id)

      allocate(b(1:nb),wb(1:nb))

      open(70,file='points.dat')
      open(71,file='weights.dat')

      do iread=1,nbf-1
        read(70,*)
        read(71,*)
      enddo

      do iread=1,nb
        read(70,*)b(iread)
        read(71,*)wb(iread)
      enddo

  endif


  read(unit=nfile,fmt=*) nglag,ngleg,key_num

if(key_capt==0)then
  allocate(glag(nglag))
  allocate(u(nglag))
  allocate(wglag(nglag))

call gaulag(glag,wglag,nglag,0._id)
do i=1,nglag
!    if(glag(i).gt.rmesh(meshr,1))then
    if(glag(i).gt.300d0)then
      nglag=i
      exit
    endif
enddo

if(myid==0)then
  do i=1,nglag
      write(77,*)glag(i),0.
  enddo

  write(unit=6,fmt='(" Number of truncated Gauss-Laguerre points ",I3,es15.6)') nglag,glag(nglag)

endif
!  print*,'___________________________________________________________'
!  print*,glag(:)
!  print*,'___________________________________________________________'
!  print*,wglag(:)
!  print*,'___________________________________________________________'
!  stop


  if(key_num==0)then

  allocate(gleg(ngleg))
  allocate(wgleg(ngleg))

!if(myid==0)then
!endif

  call gauleg(-1._id,1._id,gleg(1:ngleg),wgleg(1:ngleg),ngleg)
if(myid==0)then
    write(unit=6,fmt='(" Number of Gauss-Legendre points ",I3)') ngleg
endif


allocate(plma(0:lmax_t,-lmax_t:lmax_t,1:nglag,1:ngleg))
allocate(plmb(0:lmax_t,-lmax_t:lmax_t,1:nglag,1:ngleg))
plma=0._id
plmb=0._id
allocate(jint(-2*mmax_t:2*mmax_t,1:nglag,1:ngleg))

!allocate(nps_t(0:0))

!nps_t(0)=nmax_t(0)

allocate(basis_a(1:nmax_t(0),0:lmax_t,1:nglag,1:ngleg))
allocate(basis_b(1:nmax_t(0),0:lmax_t,1:nglag,1:ngleg))
allocate(basis_bres(1:nmax_t(0),0:lmax_t,1:nglag,1:ngleg))

  endif   !if(key_num==0)
endif


  read(unit=nfile,fmt=*) dmin,step_max,nstep

ALLOCATE(step(-nstep:nstep),TMP_STP(1:nstep))


DO I=1,NSTEP
  TMP_STP(I)=dmin*(step_max/dmin)**(REAL(I,ID)/REAL(nstep,ID))
ENDDO

DO I=-NSTEP,-1
  STEP(I)=-TMP_STP(-I)
ENDDO

STEP(0)=0D0

DO I=1,NSTEP
  STEP(I)=TMP_STP(I)
ENDDO

!$acc update device(step,nstep)

!  call expgrid

  close(unit=nfile)
if(myid==0)then
    write(unit=6,fmt=*) "All input parameters have been read"
    write(unit=6,fmt=*) "-----------------------------------"
    write(unit=6,fmt=*) "-----------------------------------"
endif



allocate(vmat_dir(1:num_t,1:num_t))
allocate(vmat_ex(1:num_t,1:num_t))
allocate(vmat_exres(1:num_t,1:num_t))
allocate(vmat_obk(1:num_t,1:num_t))
allocate(wig_d(0:lmax_t,-lmax_t:lmax_t,-lmax_t:lmax_t))
wig_d(0:lmax_t,-lmax_t:lmax_t,-lmax_t:lmax_t)=0._id


if(key_capt==0)then
   allocate(vleft(1:2*num_t,1:2*num_t))
   allocate(vright(1:2*num_t,1:2*num_t))
   allocate(vmid(1:2*num_t,1:2*num_t))
   allocate(vmatl(1:2*num_t,1:2*num_t))
   allocate(vover(1:2*num_t,1:2*num_t))
   allocate(vltmp(1:num_t,1:num_t))
   allocate(vrtmp(1:num_t,1:num_t))
else
   allocate(vleft(1:num_t,1:num_t))
   allocate(vright(1:num_t,1:num_t))
   allocate(vmid(1:num_t,1:num_t))
endif
allocate(v7m(1:2*num_t,1:2*num_t))
allocate(v23(1:2*num_t,1:2*num_t))
allocate(v7p(1:2*num_t,1:2*num_t))

if(myid==0)then
print*,'Getting out of readin'
endif

!$acc update device(num_t,num_p,k_t,k_p,enn_t,enn_p,inl_t,inl_p,fdir,idir,jdir,i0,l_t,l_p,m_t,m_p,v,projectile_charge,target_charge,nlnum_t,nlnum_p,njdouble,jdouble)

!$acc update device(gleg,wgleg,glag,wglag,u,plma,plmb,jint)
!$acc update device(ngleg,nglag,radwf_t2,radwf_p2)
!$acc update device(n_t,n_p,en_t,en_p)
!$acc update device(lar,mar,lmar,istar)

!stop


end subroutine readin



subroutine get_plm(bigr,rho)
use precision_type,only:id
use data_base
use acc_routines,only:plgndr
!use ifport
implicit none
integer::i,j,l,m
real(kind=id)::arg,bigr,rho,arg_bess
!real(kind=id),dimension(1:nglag)::u

!u(:)=glag(:)/alf_t(0)/bigr+1._id
u(:)=glag(:)/bigr+1._id



!$omp parallel&
!$omp default(none)&
!$omp private(i,j,l,m,arg,arg_bess)&
!$omp shared(nglag,ngleg,rho,v,gleg,jint,lmax_t,mmax_t,u,plma,plmb)
!$omp do&
!$omp schedule(dynamic)
do i=1,nglag
  do j=1,ngleg
    do l=0,lmax_t
      do m=-l,l,1
        arg=u(i)*gleg(j)+1._id
        arg=arg/(u(i)+gleg(j))
        plma(l,m,i,j)=plgndr(l,m,arg)
        arg=u(i)*gleg(j)-1._id
        arg=arg/(u(i)-gleg(j))
        plmb(l,m,i,j)=plgndr(l,m,arg)
      enddo
    enddo
    arg_bess=0.5_id*v*rho*sqrt((u(i)*u(i)-1._id)*(1._id-gleg(j)*gleg(j)))
    do m=0,2*mmax_t
       jint(m,i,j)=bessel_jn(m,arg_bess)
       if(mod(m,2)==0)then
         jint(-m,i,j)=jint(m,i,j)
       else
         jint(-m,i,j)=-jint(m,i,j)
       endif
    enddo
  enddo
enddo
!$omp end do
!$omp end parallel

!stop

end

subroutine precalc
use precision_type,only:id
use data_base,only:lmax_t,cleba
use acc_routines
implicit none
integer::lmd,myu,la,ma,lap,map,i,j
integer,allocatable,dimension(:)::l1,m1,l2,m2,l3,m3
real(kind=id)::tmp

i=0
do lmd=0,2*lmax_t
   do la=0,lmax_t
      do lap=0,lmax_t
         do myu=-lmd,lmd
            do ma=-la,la
               do map=-lap,lap
                  i=i+1
               enddo
            enddo
         enddo
      enddo
   enddo
enddo

allocate(l1(1:i))
allocate(l2(1:i))
allocate(l3(1:i))
allocate(m1(1:i))
allocate(m2(1:i))
allocate(m3(1:i))

i=0
do lmd=0,2*lmax_t
   do la=0,lmax_t
      do lap=0,lmax_t
         do myu=-lmd,lmd
            do ma=-la,la
               do map=-lap,lap
                  i=i+1
                  l1(i)=lmd
                  m1(i)=myu
                  l2(i)=la
                  m2(i)=ma
                  l3(i)=lap
                  m3(i)=map
               enddo
            enddo
         enddo
      enddo
   enddo
enddo

allocate(cleba(0:2*lmax_t,-2*lmax_t:2*lmax_t,0:lmax_t,-lmax_t:lmax_t,0:lmax_t,-lmax_t:lmax_t))

!$acc data copyin(l1,m1,l2,m2,l3,m3)
!$acc parallel loop 
!!$omp parallel do default(none) shared(i,l1,m1,l2,m2,l3,m3,cleba) schedule(guided)
do j=1,i
  cleba(l1(j),m1(j),l2(j),m2(j),l3(j),m3(j))=clebsch_gordan(l1(j),m1(j),l2(j),m2(j),l3(j),m3(j))
enddo
!!$omp end parallel do
!$acc end parallel loop 
!$acc end data

!$acc update host(cleba)

end subroutine precalc


!======================================================================
!  This subroutine sets up two grids: gridx and rmesh. Gridx will be
!  used to solve the radial Schrodinger equation. Gridr, which is
!  typically every second point of gridx, is used to perform integrations
!  from zero to infinity using Simpsons rule. The structure is such
!  that a wave with momentum QMAX will have NPWAVE points per oscillation. This
!  determines HMAX, the largest step. At some point, XDBLE, the intervals,
!  dX, are progressively halved. This is done to accomodate irregular
!  solutions.
!  INPUT:
!   qmax  - the biggest momentum (a.u.) that can be calculated by this mesh.
!   rmax  - the largest "r=x" in the meshes.
!   nmaxx - array declaration parameter of first declaration of GRIDX
!   nmaxr - array declaration parameter of first declaration of rmesh
!   nmaxj - array declaration parameter of first declaration of JDOUBLE
!  OUTPUT:
!   gridx - X grid
!   nx    - Number of X points
!   rmesh - R grid
!   meshr    - number of R points
!   jdouble - j points where dx doubles
!   njdouble - Number of doublings + 2
!======================================================================
subroutine grids(qmax,ndouble,rmax,rmesh,nmaxr,meshr,jdouble,nmaxj,njdouble)
use precision_type,only:p16
use mpi,only:myid
implicit none
integer::npwave,ndouble,nmaxr,meshr,nmaxj,njdouble,npdbl,nleft,ntot,nd,nj,jj,n,j
real(kind=p16)::hmax,rdble,rleft,hmin,dr,r,s,w,test,qmax,rmax
real(kind=p16),dimension(1:nmaxr,3):: rmesh(nmaxr,3)
integer(kind=4):: jdouble(nmaxj)
!  Set up rmesh
!  The number of points per oscillation (half-wavelength): NPWAVE
npwave = 6
!  The number of doubling of intervals is NDOUBLE
njdouble = ndouble + 2
if (njdouble.gt.nmaxj) then
if(myid==0)then
   print*,'Increase NMAXJ in call to GRIDS to at least:',njdouble
endif
   stop 'Increase NMAXJ in call to GRIDS'
end if
!  let NPDBL be the number of points with the same dx per interval
npdbl = 40
npdbl = 32

!  Make sure NPDBL is even
npdbl=(npdbl/2) * 2
if (npdbl.lt.4) then
if(myid==0)then
   print*,'Need to have at least 4 equally spaced points in rmesh'
endif
   stop 'Need to have at least 4 equally spaced points in rmesh'
end if
hmax = real(3.14,p16)/real(npwave,p16)/qmax
rdble = real(npdbl,p16) * hmax * (2**ndouble - 1) / real(2**ndouble,p16)
rleft = rmax - rdble
nleft = int(rleft / hmax) / 2 * 2
ntot = nleft + npdbl * ndouble
if(myid==0)then
print*,'Estimated R max:',rdble + nleft * hmax,ntot, hmax * (real(npdbl,p16) * (2**ndouble - 1) / real(2**ndouble,p16) + nleft)
endif
hmin = hmax/real(2**ndouble,p16)
!  The value of the R point from which dR is constant, = HMAX, is RDBLE
rdble = real(npdbl,p16) * hmin * (2**ndouble-1)

if(myid==0)then
print*,'Grid R parameters:'
print*,'NDOUBLE:',ndouble
print*,'HMIN:',hmin
print*,'HMAX:',hmax
print*,'NPDBL:',npdbl
print*,'RDBLE:',rdble
endif

dr = hmin
jdouble(1) = 1
do nd = 2, ndouble + 1
   jdouble(nd) = npdbl * (nd - 1)
end do

r=0._p16
j = 0
do nd = 1, ndouble
   do nj = 1, npdbl
      j = j + 1
      rmesh(j,1) = r + real(nj,p16) * dr
      rmesh(j,2) = dr
!  Simpson's rule weights
      rmesh(j,3) = real(mod(j,2) * 2 + 2,p16) * dr / 3._p16
   end do
   rmesh(j,3) = dr
   r = rmesh(j,1)
   dr = dr * 2._p16
end do

if (abs((dr - hmax)/hmax).gt.1e-4) then
if(myid==0)then
   print*,'DR and HMAX should be equal'
endif
   stop 'DR and HMAX should be equal'
end if
nj = nint((rmax-r)/hmax) / 2 * 2
if (nj+j.gt.nmaxr) then
   nj = nmaxr - j
if(myid==0)then
   print*,'Warning: meshr had to be reduced to NMAXR'
endif
end if
do jj = 1, nj
   j = j + 1
   rmesh(j,1) = r + real(jj,p16) * dr
   rmesh(j,2) = dr
!  Simpson's rule weights
   rmesh(j,3) = real(mod(j,2) * 2 + 2,p16) * dr / 3._p16
end do
!  Set the last weight so that the integrals ended exactly at rmax.
!  This is done so that the tail integrals could be done correctly.
rmesh(j,3) = dr / 3._p16
meshr = j
jdouble(ndouble+2) = j
if (mod(j,2).ne.0)print*,'Warning expected J to be even in grids'

do n = 0,0
   if (n.eq.0) then
      s = 2._p16
   else
      s = 3._p16
   endif
   w = 1._p16
   test = 0._p16
   do j = 2**n, meshr, 2**n
      w = s - w
      test = test + rmesh(j,1) ** 2 * rmesh(j,3) * w
   enddo
   if (j-2**n.ne.meshr)print*,'Warning last J should be meshr', j-2**n,meshr
   test = test * 2**n
if(myid==0)then
   print*,2**n,test*3._p16/rmesh(meshr,1)**3
endif
enddo
if(myid==0)then
print*,'Last R and meshr:', rmesh(meshr,1),meshr
endif

return
end







subroutine vec_linspace ( n, a_first, a_last, a )

  use precision_type
  implicit none
  integer(kind=4),intent(in):: n
  real(kind=id),intent(out):: a(n)
  real(kind=id),intent(in):: a_first,a_last
  integer ( kind = 4 ):: i

  if ( n == 1 ) then

    a(1) = ( a_first + a_last ) /2._id 

  else

    do i = 1, n
      a(i) = ( real( n - i,id ) * a_first + real(i - 1,id) * a_last ) / real(n - 1,id)
    end do

  end if

  return
end



subroutine get_cross_sections
   use precision_type
   use data_base
   use acc_routines,only:rylm
!   use omp_lib
   use openacc
   use mpi,only:min_rho,max_rho,nodeid,nodes,myid,ierr
   implicit none
include 'mpif.h'
   real(kind=id),dimension(1:num_t+num_p) :: prob,sig,sigb
   real(kind=id),allocatable,dimension(:) :: tdcs_tot,ddcs_tot,ddcs_tot2,sdcs_tot,sdcs_e_tot,sdcs_e_tot2
   real(kind=id) :: btmp,wtmp,btmp1,kappa,kappa_p,kf,ki,ovrlp,sp,spv,sp1,sp2,sp3,sp4
   real(kind=id) :: tot_prob_dir,tot_ion_prob_dir,gtot_ion_prob_dir,tot_sig_dir,tot_ion_sig_dir,gtot_ion_sig_dir
   real(kind=id) :: tot_prob_dir1,tot_ion_prob_dir1
   real(kind=id) ::tot_prob_ex,tot_ion_prob_ex,gtot_ion_prob_ex,tot_sig_ex,tot_ion_sig_ex,gtot_ion_sig_ex
   real(kind=id) ::tot_prob_ex1,tot_ion_prob_ex1
   real(kind=id) :: btot_prob_dir,btot_ion_prob_dir,btot_sig_dir,btot_ion_sig_dir
   real(kind=id) :: btot_prob_ex,btot_ion_prob_ex,btot_sig_ex,btot_ion_sig_ex
   real(kind=id) :: cost,costL,dsc,qt,q2,rev,test
   real(kind=id) :: sion,sionx,siond,summax,summad,summa,factor,costh,costh_d,rylml,rylmlp!,rylm
!   real(kind=id) :: xmin,xmax,dx,z,dz,y,tt,ttt,temp,ss,radwf,x,zarg,ticsnum,ticsnum2
   real(kind=id) :: tt,ttt,temp,ss,radwf,x,zarg,ticsnum,ticsnum2,test,test1,test2,dir,exc
   real(kind=id),allocatable,dimension(:) ::ddcs_scat_x,ddcs_scat_d,ddcs_scat_coh,sdcs_ej_x,sdcs_ej_d,summal
   real(kind=id),allocatable,dimension(:,:) ::ddcs_dd,ddcs_xx
   real(kind=id),allocatable,dimension(:) :: bnew,wbnew
   complex(kind=id),allocatable,dimension(:) :: ampnew
   integer(kind=4) :: ixmax,izmax,ix,iz,keyp,nbnew
   complex(kind=id),dimension(0:360) ::amp_x,amp_d,tdcs_av_x,tdcs_av_d,tdcs_av_coh
   complex(kind=id),dimension(1:num_t+num_p,1:num_t+num_p) :: ampfi
   complex(kind=id) :: ztdcs,ztmp,coulphase,zsum,zsumm,zddcs_e,zintr,zsdcs_e,zsdcs_e2,c_phase,ztmp2
   complex(kind=id) :: ztdcs_d,ztmp_d,zsum_d,zsumm_d,zddcs,zddcs_e2
   complex(kind=id),allocatable,dimension(:,:,:,:) :: ampfi_int
   complex(kind=id),allocatable,dimension(:) :: ampfi_int_ar
   complex(kind=id),allocatable,dimension(:,:) :: sigmat
   integer(kind=4) :: ir,ir1,j,j_d,i,error,n,l,m,ip,mdif,mdifa,it,lp,jp,nen,ii,nn,ll,mm,jj
   character*3:: extb
   character*2:: extn


    open(12, file='tot_prob', status='unknown')
    open(30, file='gtot_prob', status='unknown')

    allocate(amp(1:num_t+num_p,1:nb))
    allocate(ampb(1:num_t+num_p))


    tot_sig_dir=0._id
    tot_ion_sig_dir=0._id
    gtot_ion_sig_dir=0._id
    tot_sig_ex=0._id
    tot_ion_sig_ex=0._id
    gtot_ion_sig_ex=0._id

    btot_sig_dir=0._id
    btot_ion_sig_dir=0._id
    btot_sig_ex=0._id
    btot_ion_sig_ex=0._id

    sig(1:num_t+num_p)=0._id
    sigb(1:num_t+num_p)=0._id

if(myid==0)then
  print*,'bbb',b(1:nb)
endif
!  allocate(bnl(1:meshr,nps_t(0),0:lmax_t))
!  allocate(gnl(nps_t(0),0:lmax_t))
!    call getoverlap
!  print*,'bbb',b(1:nb)
!
!    print*,'gnl',gnl(7,0),gnl(6,1),gnl(5,2),gnl(4,3),'end gnl'
!    print*,'bnl',bnl(1,7,0),bnl(1,6,1),bnl(1,5,2),bnl(1,4,3),'end bnl'
!


!  stop


!    call get_c8

!inquire(file='prob',exist=probfile)

if(.not.probfile)then

if(myid==0)then
   write(6,fmt='(A)')'Total cs, tics, norm'
endif




if(myid==0)then
   write(extb,fmt="(i3.3)")1
   open(26, file='prob'//'_'//extb(1:3), status='unknown')
   open(28, file='amps_ex'//'_'//extb(1:3), status='unknown')
   open(29, file='amps_dir'//'_'//extb(1:3), status='unknown')
   open(1298, file='norm'//'_'//extb(1:3), status='unknown')
!   if(b(1).lt.1d-10)then
      write(unit=26,fmt='(i4,5es18.8)') 1,0._id,0._id,0._id,0._id,0._id
      write(unit=28,fmt='(i4,5000000es18.8)') 1,0._id,(0._id,j=1,num_p)
      write(unit=29,fmt='(i4,5000000es18.8)') 1,0._id,(0._id,j=1,num_t)
      write(unit=1298,fmt='(2es18.8)') 0._id,1._id
      print*,1,0._id
!   endif
endif

   min_rho=(nodeid-1)*int(nb/nodes)+1
   max_rho=nodeid*int(nb/nodes)

!   do ir=1,nb
   do ir=min_rho,max_rho

      write(extb,fmt="(i3.3)")ir+nbf-1
      open(26, file='prob'//'_'//extb(1:3), status='unknown')
      open(28, file='amps_ex'//'_'//extb(1:3), status='unknown')
      open(29, file='amps_dir'//'_'//extb(1:3), status='unknown')
      open(1298, file='norm'//'_'//extb(1:3), status='unknown')
      open(699, file='normZ'//'_'//extb(1:3), status='unknown')



      btmp=b(ir)
      wtmp=wb(ir)
      prob(:)=0d0
      amp(:,:)=dcmplx(0d0,0d0)
      if(btmp==0d0)then
          goto 10
      else
          if(key==0)then
!              call solveRK_num_born( btmp,prob(1:num_t+num_p),amp(1:num_t+num_p,ir) )
               call RK_1c_FBA_gpu( btmp,prob(1:num_t+num_p),amp(1:num_t+num_p,ir) )
          else
              if(key_capt==0)then
!                  call solveRK_num_2c( btmp,prob(1:num_t+num_p),amp(1:num_t+num_p,ir) )

        call RK_2c_gpu( btmp,prob(1:num_t+num_p),amp(1:num_t+num_p,ir) )

!                  call solveRK_num_1c_prj( btmp,prob(1:num_t+num_p),amp(1:num_t+num_p,ir) )

!                  call solveRK_num_2c_gpu_parallel( btmp )
!                  call solveRK_num_2c_cpu( btmp )
!                  call mat_2c_gpu_parallel( btmp )
!                  call mat_2c_cpu( btmp )
              else
                  call solveRK_num_1c_dev( btmp,prob(num_t+1:num_t+num_p),amp(num_t+1:num_t+num_p,ir) )
!                  call RK_1c_FBA_gpu( btmp,prob(1:num_t+num_p),amp(1:num_t+num_p,ir) )
              endif
          endif
      endif

10    continue
      write(unit=28,fmt='(i4,5000000es18.8)') ir+nbf-1,dble(btmp),(amp(j,ir),j=1,num_t)
      write(unit=29,fmt='(i4,5000000es18.8)') ir+nbf-1,dble(btmp),(amp(j,ir),j=num_t+1,num_t+num_p)
      tot_prob_ex=0._id
      do i=1,num_t
        tot_prob_ex=tot_prob_ex+prob(i)
      enddo
      tot_prob_dir=0._id
      do i=num_t+1,num_t+num_p
        tot_prob_dir=tot_prob_dir+prob(i)
      enddo
      tot_ion_prob_ex=0._id
      tot_ion_prob_dir=0._id
      do j=1,num_t
         if(dble(en_t(l_t(j))%enl_t(n_t(j))).gt.0d0) then
            tot_ion_prob_ex=tot_ion_prob_ex+prob(j)
            tot_ion_prob_dir=tot_ion_prob_dir+prob(num_t+j)
         endif
      enddo
      write(unit=26,fmt='(i4,5000000es18.8)') ir+nbf-1,dble(btmp),dble(tot_prob_ex),dble(tot_ion_prob_ex),dble(tot_prob_dir),dble(tot_ion_prob_dir)

      write(6,fmt='(A,i4,5000000es18.8)')'Impact parameter=',ir+nbf-1,dble(btmp),dble(tot_prob_dir+tot_prob_ex),dble(tot_ion_prob_dir+tot_ion_prob_ex),unitarity2

   enddo

else


!   allocate(tdcs_tot(0:360),ddcs_tot(0:180),sdcs_e_tot(0:180),sdcs_tot(nneg+1:nmax_t(0)))
!   tdcs_tot(0:180)=0d0;ddcs_tot(0:180)=0d0;sdcs_e_tot(0:180)=0d0;sdcs_tot(nneg+1:nmax_t(0))=0d0

  print*,'I am here second run'

   allocate(tdcs_tot(0:360),ddcs_tot(0:180),ddcs_tot2(0:180),sdcs_e_tot(0:180),sdcs_tot(nneg+1:nmax_t(0)),sdcs_e_tot2(0:180))
   tdcs_tot(0:180)=0d0;ddcs_tot(0:180)=0d0;ddcs_tot2(0:180)=0d0;
!   sdcs_e_tot(0:180)=0d0;sdcs_tot(nneg+1:nmax_t(0))=0d0;sdcs_e_tot2(nneg+1:nmax_t(0))=0d0
   sdcs_e_tot(0:180)=0d0;sdcs_tot(nneg+1:nmax_t(0))=0d0;sdcs_e_tot2(0:180)=0d0



   open(26, file='prob')
   open(28, file='amps_ex')
   open(280, file='all_prob_ex')
   open(29, file='amps_dir')
   open(290, file='all_prob_dir')

   npt=200
   allocate(zap(1:num_t+num_p,0:npt))
   allocate(pt(0:npt),wpt(0:npt))
   allocate(ddcs_scat_d(0:npt),ddcs_scat_x(0:npt),ddcs_scat_coh(0:npt))
   call gauleg(0d0,5d0,pt(1:npt),wpt(1:npt),npt)

   zap(1:num_t+num_p,0:npt)=dcmplx(0d0,0d0)

!   pt(0)=0.197839d0     !Corresponds to q=0.25 a.u.
!   pt(0)=0.0580581284021!10keV_1eV_0.1 mrad
!   pt(0)=0.1590995373947424!75keV_1eV_0.1mrad
!   pt(0)=0.2598242683237904_id!200keV_1eV_0.1mrad
!   pt(0)=0._id
   pt(0)=0.5d0
!   pt(0)=0.d0
   wpt(0)=0.d0

   ampfi(1:2*num_t,1:2*num_t)=dcmplx(0d0,0d0)

   allocate(sigmat(1:2*num_t,1:2*num_t))

   sigmat(1:2*num_t,1:2*num_t)=zero 

   do ir=1,nb


      btmp=b(ir)
      wtmp=wb(ir)

      read(unit=28,fmt='(i4,5000000es18.8)') ir1,btmp1,(amp(j,ir),j=1,num_t)
      write(unit=280,fmt='(i4,5000000es18.8)') ir1,btmp1,(btmp1*abs(amp(j,ir))**2,j=1,num_t)
      read(unit=29,fmt='(i4,5000000es18.8)') ir1,btmp1,(amp(j,ir),j=num_t+1,num_t+num_p)
      write(unit=290,fmt='(i4,5000000es18.8)') ir1,btmp1,(btmp1*abs(amp(j,ir))**2,j=num_t+1,num_t+num_p)

      do i=1,2*num_t
        do j=1,2*num_t
          ampfi(i,j)=ampfi(i,j)+amp(i,ir)*conjg(amp(j,ir))*btmp*wtmp
        enddo
      enddo

      tot_prob_ex=0._id
      do i=1,num_t
        tot_prob_ex=tot_prob_ex+btmp*amp(i,ir)*conjg(amp(i,ir))
      enddo
      tot_prob_dir=0._id
      do i=num_t+1,num_t+num_p
        tot_prob_dir=tot_prob_dir+btmp*amp(i,ir)*conjg(amp(i,ir))
      enddo
      tot_ion_prob_ex=0._id
      tot_ion_prob_dir=0._id
      do j=1,num_t
         if(dble(en_t(l_t(j))%enl_t(n_t(j))).gt.0d0) then
            tot_ion_prob_ex=tot_ion_prob_ex+btmp*amp(j,ir)*conjg(amp(j,ir))
            tot_ion_prob_dir=tot_ion_prob_dir+btmp*amp(num_t+j,ir)*conjg(amp(num_t+j,ir))
         endif
      enddo
      read(unit=26,fmt='(i4,5000000es18.8)') ir1,btmp1,tot_prob_ex1,tot_ion_prob_ex1,tot_prob_dir1,tot_ion_prob_dir1

      if(abs(tot_prob_ex1-tot_prob_ex).gt.1e-3_id.or.abs(tot_ion_prob_ex1-tot_ion_prob_ex).gt.1e-3_id.or.abs(tot_prob_dir1-tot_prob_dir).gt.1e-3_id.or.abs(tot_ion_prob_dir1-tot_ion_prob_dir).gt.1e-3_id)then
        print*,'btmp',btmp1,btmp
        print*,'tot_prob_ex',tot_prob_ex1,tot_prob_ex
        print*,'tot_ion_prob_ex1',tot_ion_prob_ex1,tot_ion_prob_ex
        print*,'tot_prob_dir1',tot_prob_dir1,tot_prob_dir
        print*,'tot_ion_prob_dir',tot_ion_prob_dir1,tot_ion_prob_dir
        stop'Read data does not correspond to calculated ones'
      endif


      tot_sig_dir=tot_sig_dir+wtmp*tot_prob_dir
      tot_ion_sig_dir=tot_ion_sig_dir+wtmp*tot_ion_prob_dir

      tot_sig_ex=tot_sig_ex+wtmp*tot_prob_ex
      tot_ion_sig_ex=tot_ion_sig_ex+wtmp*tot_ion_prob_ex


      do i=1,num_t+num_p
        sig(i)=sig(i)+wtmp*btmp*amp(i,ir)*conjg(amp(i,ir))
      enddo

      do i=1,num_t+num_p
        do ii=1,num_t+num_p
          sigmat(i,ii)=sigmat(i,ii)+wtmp*btmp*amp(i,ir)*conjg(amp(ii,ir))
        enddo
      enddo

      print*,'Impact parameter collection=',ir,dble(btmp),dble(tot_ion_prob_dir),tot_ion_prob_dir1

     do i=1,num_t
          mdif=m_t(i)-m_t(i0)
          mdifa=iabs(mdif)
          if(mdif.lt.0)then
            factor=(-1d0)**mdif
          else
            factor=1d0
          endif
          do ip=0,npt
            if(pt(ip).lt.1d-4)then
              if(mdifa==0)then
                zap(i,ip)=zap(i,ip)+factor*wtmp*btmp*amp(i,ir)
              else
                zap(i,ip)=dcmplx(0d0,0d0)
              endif
            else
              zap(i,ip)=zap(i,ip)+factor*wtmp*btmp*bessel_jn(mdifa,btmp*pt(ip))*amp(i,ir)
            endif
          enddo
     enddo

     do i=1+num_t,num_t+num_p
          mdif=m_t(i-num_t)-m_t(i0)
          mdifa=iabs(mdif)
          if(mdif.lt.0)then
            factor=(-1d0)**mdif
          else
            factor=1d0
          endif
          do ip=0,npt
            if(pt(ip).lt.1d-4)then
              if(mdifa==0)then
                zap(i,ip)=zap(i,ip)+factor*wtmp*btmp*amp(i,ir)
              else
                zap(i,ip)=dcmplx(0d0,0d0)
              endif
            else
              zap(i,ip)=zap(i,ip)+factor*wtmp*btmp*bessel_jn(mdifa,btmp*pt(ip))*amp(i,ir)
            endif
          enddo
     enddo


   enddo


   nbnew=5000
   allocate(bnew(1:nbnew),wbnew(1:nbnew),ampnew(1:nbnew))
   call gauleg(0d0,b(nb),bnew(1:nbnew),wbnew(1:nbnew),nbnew)


   do i=1,num_t

      do ir=1,nbnew
        call get_amp_new(nb,b(1:nb),amp(i,1:nb),bnew(ir),ampnew(ir))
      enddo

      mdif=m_t(i)-m_t(i0)
      mdifa=iabs(mdif)
      if(mdif.lt.0)then
        factor=(-1d0)**mdif
      else
        factor=1d0
      endif
      do ip=0,npt
        zap(i,ip)=zero
        do ir=1,nbnew
           if(pt(ip).lt.1d-4)then
             if(mdifa==0)then
               zap(i,ip)=zap(i,ip)+factor*wbnew(ir)*bnew(ir)*ampnew(ir)
             else
               zap(i,ip)=dcmplx(0d0,0d0)
             endif
           else
             zap(i,ip)=zap(i,ip)+factor*wbnew(ir)*bnew(ir)*bessel_jn(mdifa,bnew(ir)*pt(ip))*ampnew(ir)
           endif
        enddo
      enddo

   enddo


   do i=num_t+1,num_t+num_p

      do ir=1,nbnew
        call get_amp_new(nb,b(1:nb),amp(i,1:nb),bnew(ir),ampnew(ir))
      enddo

      mdif=m_t(i-num_t)-m_t(i0)
      mdifa=iabs(mdif)
      if(mdif.lt.0)then
        factor=(-1d0)**mdif
      else
        factor=1d0
      endif
      do ip=0,npt
        zap(i,ip)=zero
        do ir=1,nbnew
           if(pt(ip).lt.1d-4)then
             if(mdifa==0)then
               zap(i,ip)=zap(i,ip)+factor*wbnew(ir)*bnew(ir)*ampnew(ir)
             else
               zap(i,ip)=dcmplx(0d0,0d0)
             endif
           else
             zap(i,ip)=zap(i,ip)+factor*wbnew(ir)*bnew(ir)*bessel_jn(mdifa,bnew(ir)*pt(ip))*ampnew(ir)
           endif
        enddo
      enddo


   enddo


do i=1,num_t+num_p
  do ii=1,num_t+num_p
    sigmat(i,ii)=ft*sigmat(i,ii)
  enddo
enddo

open(310,file='density_matrix')
open(311,file='density_matrix_per_energy')
    


write(unit=310, fmt='(A18,5x)', Advance = 'No')'#(001) Energy, keV'
ii=1
do n=2,4
  do l=0,min(lmax_t,n-1)
    do m=-l,l
      do nn=2,4
        do ll=0,min(lmax_t,nn-1)
          do mm=-ll,ll
             ii=ii+2
             j=istar(n-l,l,m)+num_t
             jj=istar(nn-ll,ll,mm)+num_t
             write(unit=310,fmt='(A1,i4.4,A1,i4.4,A1)', Advance = 'No')'(',ii-1,',',ii,')'           
             write(unit=310,fmt='(i2.1)', Advance = 'No')n,l,m           
             write(unit=310,fmt='(A1)', Advance = 'No')','           
             write(unit=310,fmt='(i2.1)', Advance = 'No')nn,ll,mm           
             write(unit=310,fmt='(5x)', Advance = 'No')           
          enddo
        enddo
      enddo  
    enddo
  enddo
enddo
write(unit=310, *, Advance = 'Yes')


write(unit=310, fmt='(7es15.5,5x)', Advance = 'No')ein
write(unit=311, *)'Index,   n1,l1,m1,n2,l2,m2,    Re[DM], Im[DM] in units of 10^-16 cm^2'
ii=0
do n=2,4
  do l=0,min(lmax_t,n-1)
    do m=-l,l
      do nn=2,4
        do ll=0,min(lmax_t,nn-1)
          do mm=-ll,ll
             ii=ii+1
             j=istar(n-l,l,m)+num_t
             jj=istar(nn-ll,ll,mm)+num_t
             write(unit=310,fmt='(7es15.5,5x)', Advance = 'No')sigmat(j,jj)           
             write(unit=311,fmt='(A1,i4.4,A1)', Advance = 'No')'(',ii,')'           
             write(unit=311,fmt='(5x)', Advance = 'No')           
             write(unit=311,fmt='(i2.1)', Advance = 'No')n,l,m           
             write(unit=311,fmt='(A1)', Advance = 'No')','           
             write(unit=311,fmt='(i2.1)', Advance = 'No')nn,ll,mm           
             write(unit=311,fmt='(5x)', Advance = 'No')      
             write(unit=311,fmt='(7es15.5)', Advance = 'Yes')sigmat(j,jj)
          enddo
        enddo
      enddo  
    enddo
  enddo
enddo


close(310)
close(311)



!      amp(num_t+1)=amp(num_t+1)+real(1,id)
!      xmin=-1.500001d0
!      xmax=1.500001d0
!      ixmax=200
!      dx=(xmax-xmin)/dble(ixmax)
!      zmin=2.500001d0
!      zmax=5.500001d0
!      izmax=200
!      dz=(zmax-zmin)/dble(izmax)
!      open(666,file='distrn.dat')
!      y=0d0
!      x=xmin
!
!      do i=0,nstep
!        print*,i,step(i)
!      enddo
!
!      tt=step(nstep)
!
!      do ix=1,ixmax
!        x=x+dx
!        z=zmin
!        do iz=1,izmax
!!        print*,ix,iz
!          z=z+dz
!          zsum=dcmplx(0d0,0d0)
!!$omp parallel&
!!$omp default(none)&
!!$omp private(j,temp,ss,zarg)&
!!$omp shared(num_t,amp,en_t,n_t,l_t,m_t,x,y,z,pi,v,tt,inl_t)&
!!$omp reduction(+:zsum)
!!$omp do&
!!$omp schedule(dynamic)
!          do j=1,num_t
!            temp=dble(en_t(l_t(j))%enl_t(n_t(j)))*tt/v
!            if(x.ge.0d0)then
!              ss=1d0
!            else
!              ss=(-1d0)**m_t(j)
!            endif
!            zarg=z/dsqrt(x*x+y*y+z*z)
!            if(dabs(zarg).gt.1d0)zarg=zarg/dabs(zarg)
!            zsum=zsum+amp(j+num_t)*dsqrt(2d0*l_t(j)+1d0)/dsqrt(4d0*pi)*rylm(l_t(j),m_t(j),zarg)*radwf(inl_t(j),dsqrt(x*x+y*y+z*z))*ss*dcmplx(dcos(temp),-dsin(temp))!*c_phase(m_t(j),x,y)
!          enddo
!!$omp end do
!!$omp end parallel
!          ttt=zsum*conjg(zsum)
!          write(666,fmt='(4es15.6)')x,z,ttt
!        enddo
!      write(666,*)
!      enddo
!      close(666)
!      print*,'tt=',tt
!      stop


   do i=1,num_t+num_p
      sig(i)=ft*sig(i)
   enddo

   tot_sig_dir=ft*tot_sig_dir
   tot_ion_sig_dir=ft*tot_ion_sig_dir

   tot_sig_ex=ft*tot_sig_ex
   tot_ion_sig_ex=ft*tot_ion_sig_ex

   open(unit=14, file='cross.dat', status='unknown')
   open(unit=31, file='cross_sections', status='unknown')

   write(unit=6,fmt='("total grand and total ionization cross sections in direct channel are",2es15.5)') dble(tot_sig_dir),dble(tot_ion_sig_dir)

   write(unit=6,fmt='("total grand and total ionization cross sections in rearrangement channel are",2es15.5)') dble(tot_sig_ex),dble(tot_ion_sig_ex)

   write(unit=14,fmt='("total grand and total ionization cross sections in direct channel are",2es15.5)') tot_sig_dir,tot_ion_sig_dir

   write(unit=14,fmt='("total grand and total ionization cross sections in rearrangement channel are",2es15.5)') tot_sig_ex,tot_ion_sig_ex

   write(unit=31,fmt='(8es15.5)') ein,tot_sig_ex,tot_sig_ex-tot_ion_sig_ex,tot_ion_sig_ex,tot_sig_dir,tot_sig_dir-tot_ion_sig_dir,tot_ion_sig_dir,tot_sig_dir-tot_ion_sig_dir-sig(1+num_t)


   open(15, file='part_cs_ex', status='unknown')
   sp=0._id
   spv=0._id
   sp1=0._id
   sp2=0._id
   do j=1,num_t
      write(unit=15,fmt='(4i4,2es15.5)')j,n_t(j),l_t(j),m_t(j),sig(j),dble(en_t(l_t(j))%enl_t(n_t(j)))
      sp=sp+sig(j)*(en_t(l_t(j))%enl_t(n_t(j))-en_t(l_t(1))%enl_t(n_t(1)))
      spv=spv+sig(j)*(en_t(l_t(j))%enl_t(n_t(j))-en_t(l_t(1))%enl_t(n_t(1))+v*v/2._id)
      if(en_t(l_t(j))%enl_t(n_t(j)).lt.0._id)then
        sp1=sp1+sig(j)*(en_t(l_t(j))%enl_t(n_t(j))-en_t(l_t(1))%enl_t(n_t(1))+v*v/2._id)
      else
        sp2=sp2+sig(j)*(en_t(l_t(j))%enl_t(n_t(j))-en_t(l_t(1))%enl_t(n_t(1))+v*v/2._id)
      endif
   enddo
   close(15)

   open(16, file='part_cs_dir', status='unknown')
   sp3=0._id
   sp4=0._id
   do j=num_t+1,num_t+num_p
      write(unit=16,fmt='(4i4,2es15.5)') j-num_t,n_t(j-num_t),l_t(j-num_t),m_t(j-num_t),sig(j),dble(en_t(l_t(j-num_t))%enl_t(n_t(j-num_t)))
      sp=sp+sig(j)*(en_t(l_t(j-num_t))%enl_t(n_t(j-num_t))-en_t(l_t(1))%enl_t(n_t(1)))
      spv=spv+sig(j)*(en_t(l_t(j-num_t))%enl_t(n_t(j-num_t))-en_t(l_t(1))%enl_t(n_t(1)))
      if(en_t(l_t(j-num_t))%enl_t(n_t(j-num_t)).lt.0._id)then
        sp3=sp3+sig(j)*(en_t(l_t(j-num_t))%enl_t(n_t(j-num_t))-en_t(l_t(1))%enl_t(n_t(1)))
      else
        sp4=sp4+sig(j)*(en_t(l_t(j-num_t))%enl_t(n_t(j-num_t))-en_t(l_t(1))%enl_t(n_t(1)))
      endif
   enddo
   close(16)
   
   open(16, file='sp', status='unknown')
   write(unit=16,fmt='(3es15.5)') ein, sp*hr,spv*hr
   close(16) 

   open(16, file='spall', status='unknown')
   write(unit=16,fmt='(5es15.5)') ein, sp1*hr,sp2*hr,sp3*hr,sp4*hr
   close(16) 

   close(31)
   open(unit=31, file='pcs_sections', status='unknown')
      write(unit=31,fmt='(5000es15.5)') ein,(sig(j),j=1,num_t+num_p)
   close(31)

  allocate(summal(0:lmax_t))

  do n=1,nneg

    write(extn,fmt="(i2.2)")n
    open(31, file='EC_n='//extn(1:2), status='unknown')

    summa=0d0
    do l=0,min(lmax_t,n-1)
      summal(l)=0d0
      do m=-l,l

          i=istar(n-l,l,m)
          summa=summa+sig(i)
          summal(l)=summal(l)+sig(i)

      enddo

    enddo

    write(unit=31,fmt='(5000es15.5)') ein,summa,(summal(l),l=0,min(lmax_t,n-1))
    
    close(31)

  enddo

  do n=1,nneg

    write(extn,fmt="(i2.2)")n
    open(31, file='DS_n='//extn(1:2), status='unknown')

    summa=0d0
    do l=0,min(lmax_t,n-1)
      summal(l)=0d0
      do m=-l,l

          i=istar(n-l,l,m)+num_t
          summa=summa+sig(i)
          summal(l)=summal(l)+sig(i)

      enddo

    enddo

    write(unit=31,fmt='(5000es15.5)') ein,summa,(summal(l),l=0,min(lmax_t,n-1))
    
    close(31)

  enddo

   open(31,file='SDCS_scattering')

test1=0d0
test2=0d0
do ip=1,npt
   dir=0d0
   exc=0d0
   do n=1,nneg
     summax=0d0
     summad=0d0
     do l=0,min(lmax_t,n-1)
       do m=-l,l
          i=istar(n-l,l,m)
          j=istar(n-l,l,m)+num_t
          summax=summax+abs(zap(i,ip))**2
          summad=summad+abs(zap(j,ip))**2
       enddo
     enddo
     exc=exc+summax
     dir=dir+summad
   enddo
   write(31,11)2d0*asin(pt(ip)/2d0/k0),exc*k0*k0*a02,dir*k0*k0*a02,(abs(zap(i,ip))**2*k0*k0*a02,i=1,nneg) 
   test1=test1+pt(ip)*wpt(ip)*dir*2._id*pi*a02
   test2=test2+pt(ip)*wpt(ip)*exc*2._id*pi*a02
enddo


write(31,*)'#dir',test1
write(31,*)'#exc',test2
close(31)

   open(31,file='SDCS_ion_scattering')

test1=0d0
test2=0d0
do ip=1,npt
   dir=0d0
   exc=0d0
   do n=nneg+1,nmax_t(0)
     summax=0d0
     summad=0d0
     do l=0,min(lmax_t,n-1)
       do m=-l,l
          i=istar(n-l,l,m)
          j=istar(n-l,l,m)+num_t
          summax=summax+abs(zap(i,ip))**2
          summad=summad+abs(zap(j,ip))**2
       enddo
     enddo
     exc=exc+summax
     dir=dir+summad
   enddo
   write(31,11)2d0*asin(pt(ip)/2d0/k0),exc*k0*k0*a02,dir*k0*k0*a02,(abs(zap(i,ip))**2*k0*k0*a02,i=1,nneg) 
   test1=test1+pt(ip)*wpt(ip)*dir*2._id*pi*a02
   test2=test2+pt(ip)*wpt(ip)*exc*2._id*pi*a02
enddo


write(31,*)'#DI',test1
write(31,*)'#ECC',test2
close(31)

!   stop

   nneg=nmax_t(0)-nbin
   open(31,file='SDCS_X')

   sion=0d0
   do n=nneg+1,nmax_t(0)
     summa=0d0
     do l=0,min(lmax_t,n-1)
       do m=-l,l
          i=istar(n-l,l,m)
          summa=summa+sig(i)
       enddo
     enddo
     write(31,fmt='(4es15.5)')dble(en_t(0)%enl_t(n)),summa/dis(n-nmax_t(0)+nbin)/sqrt(2d0*en_t(0)%enl_t(n)),summa/dis(n-nmax_t(0)+nbin),summa
     sdcs_tot(n)=summa/dis(n-nmax_t(0)+nbin)/sqrt(2d0*en_t(0)%enl_t(n))
     sion=sion+summa
   enddo

   print*,'sionx=',sion

   close (31)

   open(31,file='SDCS_D')
   open(32,file='SDCS_TOT')

   sion=0d0
   do n=nneg+1,nmax_t(0)
     summa=0d0
     do l=0,min(lmax_t,n-1)
       do m=-l,l
          i=istar(n-l,l,m)+num_t
          summa=summa+sig(i)
       enddo
     enddo
     write(31,fmt='(5es15.5)')dble(en_t(0)%enl_t(n)),summa/dis(n-nmax_t(0)+nbin)/sqrt(2d0*en_t(0)%enl_t(n)),summa/dis(n-nmax_t(0)+nbin),summa,dble(en_t(0)%enl_t(n))*hr
     write(32,fmt='(3es15.5)')dble(en_t(0)%enl_t(n)),summa/dis(n-nmax_t(0)+nbin)/sqrt(2d0*en_t(0)%enl_t(n))+sdcs_tot(n)
     sion=sion+summa
   enddo


   print*,'siond=',sion

   close(31)
   close(32)

call get_sdcs_eject(1)

call mpi_finalize(ierr)
   stop

   it=30

   open(310,file='angle',status='old')
   read(310,*)it 

   write(extb,fmt="(i3.3)")it
   open(31, file='ddcs_ej_x'//'_'//extb(1:3), status='replace')
   open(32, file='ddcs_ej_tot'//'_'//extb(1:3), status='replace')


   do n=nneg+1,nmax_t(0)
     kappa=dsqrt(2.d0*en_t(0)%enl_t(n))
     costh=kappa*dcos(it*pi/180d0)-v
     costh=costh/dsqrt(kappa*kappa+v*v-2d0*kappa*v*dcos(it*pi/180d0))   !ejected electron in projectile centre
!     costh=dcos(it*pi/180d0)
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     zddcs_e=dcmplx(0d0,0d0)
     zddcs_e2=dcmplx(0d0,0d0)
     do l=0,min(lmax_t,n-1)
       do lp=0,l
         do m=-lp,lp!-lp,lp
           j=istar(n-l,l,m)
           jp=istar(n-lp,lp,m)
           ztmp=dcmplx(0d0,1d0)
           ztmp=ztmp**(l-lp)*(-1d0)**(l+lp)
           ztmp=ztmp*coulphase(-1d0/kappa,l)
           ztmp=ztmp*conjg(coulphase(-1d0/kappa,lp))
           ztmp=ztmp*dsqrt(2*l+1d0)*dsqrt(2*lp+1d0)
           rylml=rylm(l,m,costh)
           rylmlp=rylm(lp,m,costh)
           ztmp=ztmp*rylml*rylmlp
           zintr=dcmplx(0d0,0d0)
           do ip=1,npt
             zintr=zintr+zap(j,ip)*conjg(zap(jp,ip))*pt(ip)*wpt(ip)
           enddo
           if(l==lp)then
             zddcs_e=zddcs_e+ztmp*zintr
             zddcs_e2=zddcs_e2+ztmp*ampfi(j,jp)
           else
             zddcs_e=zddcs_e+2d0*ztmp*zintr
             zddcs_e2=zddcs_e2+2d0*ztmp*ampfi(j,jp)
           endif
         enddo
       enddo
     enddo
     zddcs_e=zddcs_e/2d0/dis(n-nmax_t(0)+nbin)/kappa
     zddcs_e2=zddcs_e2/2d0/dis(n-nmax_t(0)+nbin)/kappa
     write(31,11)dble(en_t(0)%enl_t(n)),zddcs_e,zddcs_e2,costh,en_t(0)%enl_t(n)*hr
     ddcs_tot(it)=dble(zddcs_e)
     ddcs_tot2(it)=dble(zddcs_e2)
   enddo

   close(31)

   open(31, file='ddcs_ej_d'//'_'//extb(1:3), status='replace')

   sion=0d0

   do n=nneg+1,nmax_t(0)
     costh=dcos(it*pi/180d0)
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     zddcs_e=dcmplx(0d0,0d0)
     zddcs_e2=dcmplx(0d0,0d0)
     do l=0,min(lmax_t,n-1)
       do lp=0,l
         do m=-lp,lp!-lp,lp
           j=istar(n-l,l,m)+num_t
           jp=istar(n-lp,lp,m)+num_t
           ztmp=dcmplx(0d0,1d0)
           ztmp=ztmp**(l-lp)*(-1d0)**(l+lp)
            ztmp=ztmp*coulphase(-1d0/dsqrt(2.d0*en_t(0)%enl_t(n)),l)
            ztmp=ztmp*conjg(coulphase(-1d0/dsqrt(2.d0*en_t(0)%enl_t(n)),lp))
           ztmp=ztmp*dsqrt(2*l+1d0)*dsqrt(2*lp+1d0)
           ztmp=ztmp*rylm(l,m,costh)*rylm(lp,m,costh)
           zintr=dcmplx(0d0,0d0)
           do ip=1,npt
             zintr=zintr+zap(j,ip)*conjg(zap(jp,ip))*pt(ip)*wpt(ip)
           enddo
           if(l==lp)then
             zddcs_e=zddcs_e+ztmp*zintr
             zddcs_e2=zddcs_e2+ztmp*ampfi(j,jp)
           else
             zddcs_e=zddcs_e+2d0*ztmp*zintr
             zddcs_e2=zddcs_e2+2d0*ztmp*ampfi(j,jp)
           endif

!           ztmp=ztmp*zintr
!           ztmp2=ztmp*ampfi(j,jp)
!           if(l==lp)then
!             ztmp=ztmp*1d0
!             ztmp2=ztmp2*1d0
!           else
!             ztmp=ztmp*2d0
!             ztmp2=ztmp2*2d0
!           endif
!           zddcs_e=zddcs_e+ztmp
!           zddcs_e2=zddcs_e2+ztmp2
         enddo
       enddo
     enddo
     zddcs_e=zddcs_e/2d0/dis(n-nmax_t(0)+nbin)/dsqrt(2.d0*en_t(0)%enl_t(n))
     zddcs_e2=zddcs_e2/2d0/dis(n-nmax_t(0)+nbin)/dsqrt(2.d0*en_t(0)%enl_t(n))
     write(31,11)dble(en_t(0)%enl_t(n)),zddcs_e,zddcs_e2,costh,en_t(0)%enl_t(n)*hr
     write(32,11)dble(en_t(0)%enl_t(n)),dble(zddcs_e)+ddcs_tot(it),dble(zddcs_e2)+ddcs_tot2(it)
     sion=sion+dble(zddcs_e2/2d0)
   enddo

   write(31,*)it,sion

   close(31)
   close(32)



   n=nneg+3 


   open(31, file='ddcs_ejang_x'//'_'//extb(1:3), status='replace')
   open(32, file='ddcs_ejang_tot'//'_'//extb(1:3), status='replace')


   do it=0,180
     kappa=dsqrt(2.d0*en_t(0)%enl_t(n))
     costh=kappa*dcos(it*pi/180d0)-v
     costh=costh/dsqrt(kappa*kappa+v*v-2d0*kappa*v*dcos(it*pi/180d0))   !ejected electron in projectile centre
!     costh=dcos(it*pi/180d0)
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     zddcs_e=dcmplx(0d0,0d0)
     zddcs_e2=dcmplx(0d0,0d0)
     do l=0,min(lmax_t,n-1)
       do lp=0,l
         do m=-lp,lp!-lp,lp
           j=istar(n-l,l,m)
           jp=istar(n-lp,lp,m)
           ztmp=dcmplx(0d0,1d0)
           ztmp=ztmp**(l-lp)*(-1d0)**(l+lp)
           ztmp=ztmp*coulphase(-1d0/kappa,l)
           ztmp=ztmp*conjg(coulphase(-1d0/kappa,lp))
           ztmp=ztmp*dsqrt(2*l+1d0)*dsqrt(2*lp+1d0)
           rylml=rylm(l,m,costh)
           rylmlp=rylm(lp,m,costh)
           ztmp=ztmp*rylml*rylmlp
           zintr=dcmplx(0d0,0d0)
           do ip=1,npt
             zintr=zintr+zap(j,ip)*conjg(zap(jp,ip))*pt(ip)*wpt(ip)
           enddo
           if(l==lp)then
             zddcs_e=zddcs_e+ztmp*zintr
             zddcs_e2=zddcs_e2+ztmp*ampfi(j,jp)
           else
             zddcs_e=zddcs_e+2d0*ztmp*zintr
             zddcs_e2=zddcs_e2+2d0*ztmp*ampfi(j,jp)
           endif
         enddo
       enddo
     enddo
     zddcs_e=zddcs_e/2d0/dis(n-nmax_t(0)+nbin)/kappa
     zddcs_e2=zddcs_e2/2d0/dis(n-nmax_t(0)+nbin)/kappa
     write(31,11)dble(it),zddcs_e,zddcs_e2,costh,en_t(0)%enl_t(n)*hr
     ddcs_tot(it)=dble(zddcs_e)
     ddcs_tot2(it)=dble(zddcs_e2)
   enddo

   close(31)

   open(31, file='ddcs_ejang_d'//'_'//extb(1:3), status='replace')

   sion=0d0

   do it=0,180
     costh=dcos(it*pi/180d0)
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     zddcs_e=dcmplx(0d0,0d0)
     zddcs_e2=dcmplx(0d0,0d0)
     do l=0,min(lmax_t,n-1)
       do lp=0,l
         do m=-lp,lp!-lp,lp
           j=istar(n-l,l,m)+num_t
           jp=istar(n-lp,lp,m)+num_t
           ztmp=dcmplx(0d0,1d0)
           ztmp=ztmp**(l-lp)*(-1d0)**(l+lp)
            ztmp=ztmp*coulphase(-1d0/dsqrt(2.d0*en_t(0)%enl_t(n)),l)
            ztmp=ztmp*conjg(coulphase(-1d0/dsqrt(2.d0*en_t(0)%enl_t(n)),lp))
           ztmp=ztmp*dsqrt(2*l+1d0)*dsqrt(2*lp+1d0)
           ztmp=ztmp*rylm(l,m,costh)*rylm(lp,m,costh)
           zintr=dcmplx(0d0,0d0)
           do ip=1,npt
             zintr=zintr+zap(j,ip)*conjg(zap(jp,ip))*pt(ip)*wpt(ip)
           enddo
           if(l==lp)then
             zddcs_e=zddcs_e+ztmp*zintr
             zddcs_e2=zddcs_e2+ztmp*ampfi(j,jp)
           else
             zddcs_e=zddcs_e+2d0*ztmp*zintr
             zddcs_e2=zddcs_e2+2d0*ztmp*ampfi(j,jp)
           endif

!           ztmp=ztmp*zintr
!           ztmp2=ztmp*ampfi(j,jp)
!           if(l==lp)then
!             ztmp=ztmp*1d0
!             ztmp2=ztmp2*1d0
!           else
!             ztmp=ztmp*2d0
!             ztmp2=ztmp2*2d0
!           endif
!           zddcs_e=zddcs_e+ztmp
!           zddcs_e2=zddcs_e2+ztmp2
         enddo
       enddo
     enddo
     zddcs_e=zddcs_e/2d0/dis(n-nmax_t(0)+nbin)/dsqrt(2.d0*en_t(0)%enl_t(n))
     zddcs_e2=zddcs_e2/2d0/dis(n-nmax_t(0)+nbin)/dsqrt(2.d0*en_t(0)%enl_t(n))
     write(31,11)dble(it),zddcs_e,zddcs_e2,costh,en_t(0)%enl_t(n)*hr
     write(32,11)dble(it),dble(zddcs_e)+ddcs_tot(it),dble(zddcs_e2)+ddcs_tot2(it)
     sion=sion+dble(zddcs_e2/2d0)
   enddo

   write(31,*)n,sion

   close(31)
   close(32)




   open(31, file='sdcs_ej_xx', status='replace')


test=0d0
do it=0,180
   sion=0d0
   rev=0d0
   do n=nneg+1,nmax_t(0)
     kappa=dsqrt(2.d0*en_t(0)%enl_t(n))
     costh=dcos(it*pi/180d0)
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     zddcs_e=dcmplx(0d0,0d0)
     zddcs_e2=dcmplx(0d0,0d0)
     do l=0,min(lmax_t,n-1)
       do lp=0,l
         do m=-lp,lp!-lp,lp
           j=istar(n-l,l,m)
           jp=istar(n-lp,lp,m)
           ztmp=dcmplx(0d0,1d0)
           ztmp=ztmp**(l-lp)*(-1d0)**(l+lp)
           ztmp=ztmp*coulphase(-1d0/kappa,l)
           ztmp=ztmp*conjg(coulphase(-1d0/kappa,lp))
           ztmp=ztmp*dsqrt(2*l+1d0)*dsqrt(2*lp+1d0)
           rylml=rylm(l,m,costh)
           rylmlp=rylm(lp,m,costh)
           ztmp=ztmp*rylml*rylmlp
           zintr=dcmplx(0d0,0d0)
           do ip=1,npt
             zintr=zintr+zap(j,ip)*conjg(zap(jp,ip))*pt(ip)*wpt(ip)
           enddo
           if(l==lp)then
             zddcs_e=zddcs_e+ztmp*zintr
             zddcs_e2=zddcs_e2+ztmp*ampfi(j,jp)
           else
             zddcs_e=zddcs_e+2d0*ztmp*zintr
             zddcs_e2=zddcs_e2+2d0*ztmp*ampfi(j,jp)
           endif
         enddo
       enddo
     enddo
     sion=sion+dble(zddcs_e2/2d0)
     rev=rev+dble(zddcs_e/2d0)
     zddcs_e=zddcs_e/2d0/dis(n-nmax_t(0)+nbin)/kappa
     zddcs_e2=zddcs_e2/2d0/dis(n-nmax_t(0)+nbin)/kappa
     ddcs_tot(it)=dble(zddcs_e)
     ddcs_tot2(it)=dble(zddcs_e2)
   enddo
   write(31,*)it,sion,rev
   test=test+sion*sin(it*pi/180._id)
enddo
write(31,*)test/180._id*ft*pi

   close(31)


allocate(ddcs_dd(1:nmax_t(0),0:180))
allocate(ddcs_xx(1:nmax_t(0),0:180))
allocate(sdcs_ej_d(0:180))

   open(31, file='sdcs_ej_d', status='replace')

test=0d0
do it=0,180,1
   sion=0d0
   rev=0d0
   do n=nneg+1,nmax_t(0)
     costh=dcos(it*pi/180d0)
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     zddcs_e=dcmplx(0d0,0d0)
     zddcs_e2=dcmplx(0d0,0d0)
     do l=0,min(lmax_t,n-1)
       do lp=0,l
         do m=-lp,lp!-lp,lp
           j=istar(n-l,l,m)+num_t
           jp=istar(n-lp,lp,m)+num_t
           ztmp=dcmplx(0d0,1d0)
           ztmp=ztmp**(l-lp)*(-1d0)**(l+lp)
            ztmp=ztmp*coulphase(-1d0/dsqrt(2.d0*en_t(0)%enl_t(n)),l)
            ztmp=ztmp*conjg(coulphase(-1d0/dsqrt(2.d0*en_t(0)%enl_t(n)),lp))
           ztmp=ztmp*dsqrt(2*l+1d0)*dsqrt(2*lp+1d0)
           ztmp=ztmp*rylm(l,m,costh)*rylm(lp,m,costh)
           zintr=dcmplx(0d0,0d0)
           do ip=1,npt
             zintr=zintr+zap(j,ip)*conjg(zap(jp,ip))*pt(ip)*wpt(ip)
           enddo
           if(l==lp)then
             zddcs_e=zddcs_e+ztmp*zintr
             zddcs_e2=zddcs_e2+ztmp*ampfi(j,jp)
           else
             zddcs_e=zddcs_e+2d0*ztmp*zintr
             zddcs_e2=zddcs_e2+2d0*ztmp*ampfi(j,jp)
           endif

!           ztmp=ztmp*zintr
!           ztmp2=ztmp*ampfi(j,jp)
!           if(l==lp)then
!             ztmp=ztmp*1d0
!             ztmp2=ztmp2*1d0
!           else
!             ztmp=ztmp*2d0
!             ztmp2=ztmp2*2d0
!           endif
!           zddcs_e=zddcs_e+ztmp
!           zddcs_e2=zddcs_e2+ztmp2
         enddo
       enddo
     enddo
     sion=sion+dble(zddcs_e2/2d0)
     rev=rev+dble(zddcs_e/2d0)
     zddcs_e=zddcs_e/2d0/dis(n-nmax_t(0)+nbin)/dsqrt(2.d0*en_t(0)%enl_t(n))
     zddcs_e2=zddcs_e2/2d0/dis(n-nmax_t(0)+nbin)/dsqrt(2.d0*en_t(0)%enl_t(n))
     ddcs_dd(n,it)=dble(zddcs_e2)
   enddo
   write(31,*)it,sion*ft
   sdcs_ej_d(it)=sion*ft
   test=test+sion*sin(it*pi/180._id)
enddo

write(31,*)test/180._id*ft*pi

close(31)




open(31, file='sdcs_ej_x', status='replace')

allocate(ampfi_int(1:nbin,0:lmax_t,0:lmax_t,-lmax_t:lmax_t),ampfi_int_ar(1:nbin))

test=0d0
do it=1,180
   do l=0,lmax_t
     do lp=0,l
       do m=-lp,lp!-lp,lp
          do n=nneg+1,nmax_t(0)
             j=istar(n-l,l,m)
             jp=istar(n-lp,lp,m)
             ampfi_int_ar(n-nmax_t(0)+nbin)=ampfi(j,jp)/dis(n-nmax_t(0)+nbin)       
          enddo
          do n=nneg+1,nmax_t(0)
             kappa=dsqrt(2.d0*en_t(0)%enl_t(n))
             kappa_p=dsqrt(kappa*kappa+v*v-2d0*kappa*v*dcos(it*pi/180d0))
             if(kappa_p.gt.bin(nbin))then
                ampfi_int(n-nmax_t(0)+nbin,l,lp,m)=zero
             else
                call zintrpl(nbin,bin(1:nbin),ampfi_int_ar(1:nbin),kappa_p,ampfi_int(n-nmax_t(0)+nbin,l,lp,m))
             endif
          enddo
       enddo
     enddo
   enddo 

   sion=0d0
   do n=nneg+1,nmax_t(0)
     kappa=dsqrt(2.d0*en_t(0)%enl_t(n))
     kappa_p=dsqrt(kappa*kappa+v*v-2d0*kappa*v*dcos(it*pi/180d0))
     costh=kappa*dcos(it*pi/180d0)-v
     costh=costh/kappa_p   !ejected electron in projectile centre
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     zddcs_e2=dcmplx(0d0,0d0)
     do l=0,lmax_t
       do lp=0,l
         do m=-lp,lp
           ztmp=dcmplx(0d0,1d0)
           ztmp=ztmp**(l-lp)*(-1d0)**(l+lp)
           ztmp=ztmp*coulphase(-1d0/kappa_p,l)
           ztmp=ztmp*conjg(coulphase(-1d0/kappa_p,lp))
           ztmp=ztmp*dsqrt(2*l+1d0)*dsqrt(2*lp+1d0)
           rylml=rylm(l,m,costh)
           rylmlp=rylm(lp,m,costh)
           ztmp=ztmp*rylml*rylmlp
           if(l==lp)then
             zddcs_e2=zddcs_e2+ztmp*ampfi_int(n-nmax_t(0)+nbin,l,lp,m)
           else
             zddcs_e2=zddcs_e2+2d0*ztmp*ampfi_int(n-nmax_t(0)+nbin,l,lp,m)
           endif
         enddo
       enddo
     enddo
     sion=sion+dble(zddcs_e2/2d0)/kappa_p/kappa_p*kappa*kappa*dis(n-nmax_t(0)+nbin)
     zddcs_e2=zddcs_e2/2d0/kappa_p/kappa_p*kappa
     ddcs_xx(n,it)=dble(zddcs_e2)
   enddo
   write(31,*)it,sion*ft,sdcs_ej_d(it)
   test=test+sion*sin(it*pi/180._id)
enddo
write(31,*)test/180._id*ft*pi

close(31)



do n=nneg+1,nmax_t(0)
   write(extb,fmt="(i3.3)")n
   open(32, file='ddcs_ej_tot'//'_'//extb(1:3), status='replace')
   do it=0,180
      write(32,11)dble(it),ddcs_xx(n,it),ddcs_dd(n,it),en_t(0)%enl_t(n)*hr
   enddo
   close(32)
enddo

n=5
ip=0
      open(31,file='TDCS_AV2',status='replace')
do it=0,360
    call get_tdcs(1,sqrt(k0*k0-2d0*mut*(en_t(0)%enl_t(n+nneg)-en_t(0)%enl_t(1))),k0,n+nneg,it,ip,tdcs_av_x,tdcs_av_d,tdcs_av_coh)
enddo
      close(31)

call get_sdcs_eject(1)



endif       !if(.not.probfile)

call mpi_finalize(ierr)
   stop
11 FORMAT(5000000ES18.8)

end subroutine get_cross_sections

subroutine intrpl5(xtmp,ytmp,x,tmp)
use precision_type,only:id
implicit none
real(kind=id),dimension(5):: ytmp,xtmp !2th order intrpl is used
real(kind=id)::a1,a2,a3,a4,a5,x,tmp

a1=(x-xtmp(2))*(x-xtmp(3))*(x-xtmp(4))*(x-xtmp(5))
a1=a1/(xtmp(1)-xtmp(2))/(xtmp(1)-xtmp(3))/(xtmp(1)-xtmp(4))/(xtmp(1)-xtmp(5))

a2=(x-xtmp(1))*(x-xtmp(3))*(x-xtmp(4))*(x-xtmp(5))
a2=a2/(xtmp(2)-xtmp(1))/(xtmp(2)-xtmp(3))/(xtmp(2)-xtmp(4))/(xtmp(2)-xtmp(5))

a3=(x-xtmp(1))*(x-xtmp(2))*(x-xtmp(4))*(x-xtmp(5))
a3=a3/(xtmp(3)-xtmp(1))/(xtmp(3)-xtmp(2))/(xtmp(3)-xtmp(4))/(xtmp(3)-xtmp(5))

a4=(x-xtmp(1))*(x-xtmp(2))*(x-xtmp(3))*(x-xtmp(5))
a4=a4/(xtmp(4)-xtmp(1))/(xtmp(4)-xtmp(2))/(xtmp(4)-xtmp(3))/(xtmp(4)-xtmp(5))

a5=(x-xtmp(1))*(x-xtmp(2))*(x-xtmp(3))*(x-xtmp(4))
a5=a5/(xtmp(5)-xtmp(1))/(xtmp(5)-xtmp(2))/(xtmp(5)-xtmp(3))/(xtmp(5)-xtmp(4))

tmp=a1*ytmp(1)+a2*ytmp(2)+a3*ytmp(3)+a4*ytmp(4)+a5*ytmp(5)

return
end

subroutine zintrpl5(x,y,arg,res)
use precision_type,only:id
   implicit none
   real(kind=id)::arg,rev,imv
   real(kind=id),dimension(5)::x,ry,iy
   complex(kind=id),dimension(5)::y
   complex(kind=id)::res

   ry(1:5)=dble(y(1:5))
   iy(1:5)=dimag(y(1:5))

   call intrpl5(x,ry,arg,rev)
   call intrpl5(x,iy,arg,imv)

   res=dcmplx(rev,imv)

   return
end subroutine zintrpl5

subroutine get_amp_new(nb,b,amp,arg,z)
use precision_type,only:id
implicit none
integer::       nb
real(kind=id),dimension(1:nb):: b
complex(kind=id),dimension(1:nb):: amp
real(kind=id),dimension(1:5):: xtmp
complex(kind=id),dimension(1:5):: ztmp
integer::jres
real(kind=id)::arg
complex(kind=id)::z


if(arg.gt.b(nb-2)) then
!   z=cmplx(0,0,id)
   xtmp(1:5)=b(nb-4:nb)
   ztmp(1:5)=amp(nb-4:nb)
   call zintrpl5(xtmp,ztmp,arg,z)
   goto 99
endif

if(arg.lt.b(3)) then
   xtmp(1:5)=b(1:5)
   ztmp(1:5)=amp(1:5)
   call zintrpl5(xtmp,ztmp,arg,z)
   goto 99
endif

call locate(b(1:nb),nb,arg,jres)
xtmp(1:5)=b(jres-2:jres+2)
ztmp(1:5)=amp(jres-2:jres+2)
call zintrpl5(xtmp,ztmp,arg,z)

99 continue
return
end

subroutine get_sdcs_eject(keyp)
   use precision_type
   use data_base
   implicit none
   integer::keyp,n,ip,it
   real(kind=id),dimension(0:360):: ddcs_scat_x,ddcs_scat_d,ddcs_scat_coh 
   real(kind=id),dimension(0:360):: sdcs_scat_x,sdcs_scat_d,sdcs_scat_coh 
   real(kind=id):: kf,ki, test_x,test_d,test_coh


sdcs_scat_x(0:360)=0._id
sdcs_scat_d(0:360)=0._id
sdcs_scat_coh(0:360)=0._id

  do n=1,nbin
    kf=sqrt(k0*k0-2d0*mut*(en_t(0)%enl_t(n+nneg)-en_t(0)%enl_t(1)))
    ki=k0


    call get_ddcs_eject(0,kf,ki,n+nneg,ddcs_scat_x,ddcs_scat_d,ddcs_scat_coh)

!    do it=0,180,5
    do it=0,15,1
      sdcs_scat_x(it)=sdcs_scat_x(it)+ddcs_scat_x(it)*dis(n)*dsqrt(2.d0*en_t(0)%enl_t(n+nneg))
      sdcs_scat_d(it)=sdcs_scat_d(it)+ddcs_scat_d(it)*dis(n)*dsqrt(2.d0*en_t(0)%enl_t(n+nneg))
      sdcs_scat_coh(it)=sdcs_scat_coh(it)+ddcs_scat_coh(it)*dis(n)*dsqrt(2.d0*en_t(0)%enl_t(n+nneg))
    enddo


  enddo



   if(keyp==1)then
      open(31,file='SDCS_EJECT_ALL',status='replace')
      test_x=0._id
      test_d=0._id
      test_coh=0._id
      do it=0,180,5
        write(31,11)dble(it),sdcs_scat_x(it),sdcs_scat_d(it),sdcs_scat_coh(it)
        test_x=test_x+sdcs_scat_x(it)*sin(it*pi/180._id)
        test_d=test_d+sdcs_scat_d(it)*sin(it*pi/180._id)
        test_coh=test_coh+sdcs_scat_coh(it)*sin(it*pi/180._id)
      enddo
      write(31,*)'#',test_x*pi/180d0*5*ft/2._id/pi,test_d*pi/180d0*5*ft/2._id/pi
      close(31)
   endif

11 format(5000000es18.8)


end subroutine get_sdcs_eject


subroutine get_ddcs_eject(keyp,kf,ki,n,ddcs_scat_x,ddcs_scat_d,ddcs_scat_coh)
   use precision_type
   use data_base
   implicit none
   integer::keyp,n,ip,it
   complex(kind=id),dimension(0:npt):: tdcs_av_x,tdcs_av_d,tdcs_av_coh 
   real(kind=id),dimension(0:360):: ddcs_scat_x,ddcs_scat_d,ddcs_scat_coh 
   real(kind=id):: kf,ki 

tdcs_av_x(:)=zero;tdcs_av_d(:)=zero;tdcs_av_coh(:)=zero

!$omp parallel&
!$omp default(none)&
!$omp private(tdcs_av_x,tdcs_av_d,tdcs_av_coh)&
!$omp shared(pt,wpt,npt,kf,ki,n,ddcs_scat_x,ddcs_scat_d,ddcs_scat_coh)
!$omp do&         
!$omp schedule(dynamic)


!  do it=0,180,5
  do it=0,15,1

    call get_tdcs_phi_av3(0,kf,ki,n,it,tdcs_av_x,tdcs_av_d,tdcs_av_coh)

    ddcs_scat_x(it)=0._id
    ddcs_scat_d(it)=0._id
    ddcs_scat_coh(it)=0._id

    do ip=0,npt

      ddcs_scat_x(it)=ddcs_scat_x(it)+pt(ip)*wpt(ip)*tdcs_av_x(ip)
      ddcs_scat_d(it)=ddcs_scat_d(it)+pt(ip)*wpt(ip)*tdcs_av_d(ip)
      ddcs_scat_coh(it)=ddcs_scat_coh(it)+pt(ip)*wpt(ip)*tdcs_av_coh(ip)
      

    enddo

    ddcs_scat_x(it)=ddcs_scat_x(it)*2._id*pi/kf/ki
    ddcs_scat_d(it)=ddcs_scat_d(it)*2._id*pi/kf/ki
    ddcs_scat_coh(it)=ddcs_scat_coh(it)*2._id*pi/kf/ki

  enddo

!$omp end do         
!$omp end parallel

   if(keyp==1)then
      open(31,file='DDCS_SCAT',status='replace')
      do it=0,360
        write(31,11)dble(it),ddcs_scat_x(it),ddcs_scat_d(it),ddcs_scat_coh(it)
      enddo
      close(31)
   endif

11 format(5000000es18.8)


end subroutine get_ddcs_eject

subroutine get_tdcs(keyp,kf,ki,n,it,ip,tdcs_x_av,tdcs_d_av,tdcs_coh_av)
   use precision_type
   use data_base
   use flogs,only:fac
   use acc_routines,only:rylm
   implicit none
   integer::  ip,keyp,it,l,m,i,j,n,mdif,mdifa,ir,l1,l2,jc,jl,jr,iphi,phi_fix
   complex(kind=id):: ztmp,coulphase,zsumm,zex,zamp 
   complex(kind=id),dimension(1:nbin,0:lmax_t,-lmax_t:lmax_t):: zp
   complex(kind=id),dimension(0:npt):: amp_x,tdcs_x_av,tdcs_d_av,tdcs_coh_av 
   complex(kind=id),dimension(0:npt,0:360):: amp_d
   complex(kind=id),dimension(1:nbin):: zintr
   real(kind=id):: kappa,ovrlp,kf,ki,costh,kappa_p,btmp,wtmp,factor,bessj,rylmp,ptr,ktr,phi,res

   kappa=sqrt(2._id*en_t(0)%enl_t(n))
   ovrlp=1d0/dis(n-nmax_t(0)+nbin)

phi_fix=0
!The direct scattering amplitude:

tdcs_x_av=0._id;tdcs_d_av=0._id;tdcs_coh_av=0._id

if(it.gt.180)then

     tdcs_d_av(ip)=zero  
     do iphi=phi_fix,phi_fix!360
        phi=dble(iphi)*pi/180d0 
        costh=dcos((it-180)*pi/180d0)
        if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
        amp_d(ip,iphi)=dcmplx(0d0,0d0)
        do l=0,lmax_t
          ztmp=dcmplx(0d0,-1d0)
          ztmp=ztmp**l*coulphase(-1d0/kappa,l)
          ztmp=ztmp*dsqrt(2*l+1d0)
          zsumm=dcmplx(0d0,0d0)
          do m=-l,l
            j=istar(n-l,l,m)+num_t
            zsumm=zsumm+rylm(l,m,costh)*zap(j,ip)*dcmplx(0d0,1d0)**(m)*exp(cu*m*phi)!*dcmplx(dcos(m*phi),dsin(m*phi))!*(-1d0)**m
          enddo
          amp_d(ip,iphi)=amp_d(ip,iphi)+ztmp*zsumm*(-1d0)**l
        enddo
        amp_d(ip,iphi)=amp_d(ip,iphi)*sqrt(kf*ki*ovrlp/4d0/pi/kappa)
        tdcs_d_av(ip)=tdcs_d_av(ip)+amp_d(ip,iphi)*conjg(amp_d(ip,iphi))
     enddo
     tdcs_d_av(ip)=tdcs_d_av(ip)!*pi/180d0*5
endif
if(it.le.180)then
     tdcs_d_av(ip)=zero  
     do iphi=phi_fix,phi_fix!360
        phi=dble(iphi)*pi/180d0 
        costh=dcos((it)*pi/180d0)
        if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
        amp_d(ip,iphi)=dcmplx(0d0,0d0)
        do l=0,lmax_t
          ztmp=dcmplx(0d0,-1d0)
          ztmp=ztmp**l*coulphase(-1d0/kappa,l)
          ztmp=ztmp*dsqrt(2*l+1d0)
          zsumm=dcmplx(0d0,0d0)
          do m=-l,l
            j=istar(n-l,l,m)+num_t
            zsumm=zsumm+rylm(l,m,costh)*zap(j,ip)*dcmplx(0d0,1d0)**(m)*exp(cu*m*phi)!*dcmplx(dcos(m*phi),dsin(m*phi))
          enddo
          amp_d(ip,iphi)=amp_d(ip,iphi)+ztmp*zsumm
        enddo
        amp_d(ip,iphi)=amp_d(ip,iphi)*sqrt(kf*ki*ovrlp/4d0/pi/kappa)
        tdcs_d_av(ip)=tdcs_d_av(ip)+amp_d(ip,iphi)*conjg(amp_d(ip,iphi))
     enddo
     tdcs_d_av(ip)=tdcs_d_av(ip)!*pi/180d0*5
endif

!goto 10
!return

!Rearrangement:

if(it.le.180)then
     ktr=kappa*dsin((it)*pi/180d0)
!     ptr=pt(ip)
     kappa_p=sqrt(kappa**2+v**2-2*kappa*v*dcos((it)*pi/180d0))

     costh=dcos((it)*pi/180d0)
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     tdcs_x_av(ip)=zero  
     tdcs_coh_av(ip)=zero  
     do iphi=phi_fix,phi_fix !360
        phi=dble(iphi)*pi/180d0 
        ptr=dsqrt(pt(ip)**2-2d0*pt(ip)*ktr*dcos(phi)+ktr**2)!dble(pt(ip)-ktr)
        amp_x(ip)=dcmplx(0d0,0d0)
        do l=0,lmax_t
          ztmp=dcmplx(0d0,-1d0)
          ztmp=ztmp**l*coulphase(-1d0/kappa_p,l)
          ztmp=ztmp*dsqrt(2*l+1d0)
          zsumm=dcmplx(0d0,0d0)
          do m=-l,l
             mdif=m-m_t(i0)
             mdifa=iabs(mdif)
             if(mdif.lt.0)then
               factor=(-1d0)**mdif
             else
               factor=1d0
             endif

             zex=dcmplx(0d0,0d0)
             do ir=1,nb
                btmp=b(ir)
                wtmp=wb(ir)

                do i=nneg+1,nmax_t(0)
                   
                   zintr(i-nneg)=amp(istar(i-l,l,m),ir)/dsqrt(dis(i-nneg)) 
!                   zintr(i-nneg)=amp(1,istar(i-l,l,m),ir)/dsqrt(dis(i-nneg)) 
          
                enddo

                call zintrpl(nbin,bin(1:nbin),zintr(1:nbin),kappa_p,zamp)

                if(dabs(ptr).lt.1d-4)then
                  if(mdifa==0)then
                    zex=zex+factor*wtmp*btmp*zamp
                  else
                    zex=dcmplx(0d0,0d0)
                  endif
                else
                    zex=zex+factor*wtmp*btmp*bessel_jn(mdifa,btmp*ptr)*zamp
                endif
  
             enddo 

            rylmp=0d0
            do l1=iabs(m),l
               l2=l-l1
               rylmp=rylmp+(-1d0)**l2*kappa**l1*v**l2/sqrt(fac(2*l1)*fac(2*l2))*cleba(l1,m,l2,0,l,m)*rylm(l1,m,costh)
            enddo
            rylmp=rylmp/kappa_p**l*sqrt(fac(2*l)) 

            res=atan2(-ktr*dsin(phi),pt(ip)-ktr*dcos(phi))
            zsumm=zsumm+rylmp*zex*dcmplx(0d0,1d0)**(m)*dcmplx(dcos(m*phi),dsin(m*phi))*dcmplx(dcos(m*res),-dsin(m*res))

          enddo
          amp_x(ip)=amp_x(ip)+ztmp*zsumm
        enddo
        amp_x(ip)=amp_x(ip)*sqrt(kf*ki/4d0/pi*kappa)/kappa_p
        tdcs_x_av(ip)=tdcs_x_av(ip)+amp_x(ip)*conjg(amp_x(ip))
        tdcs_coh_av(ip)=tdcs_coh_av(ip)+(amp_x(ip)+amp_d(ip,iphi))*conjg(amp_x(ip)+amp_d(ip,iphi))
     enddo
     tdcs_x_av(ip)=tdcs_x_av(ip)
     tdcs_coh_av(ip)=tdcs_coh_av(ip)
!   enddo
endif



if(it.gt.180)then
!   do ip=0,npt
     ktr=kappa*dsin((it)*pi/180d0)
!     ptr=dabs(pt(ip)-ktr)
     kappa_p=sqrt(kappa**2+v**2-2*kappa*v*dcos((it)*pi/180d0))

     costh=dcos((it-180)*pi/180d0)
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     tdcs_x_av(ip)=zero  
     tdcs_coh_av(ip)=zero  
     do iphi=phi_fix,phi_fix!0,360
        phi=dble(iphi)*pi/180d0 
        ptr=dsqrt(pt(ip)**2-2d0*pt(ip)*ktr*dcos(phi)+ktr**2)!dble(pt(ip)-ktr)
        amp_x(ip)=dcmplx(0d0,0d0)
        do l=0,lmax_t
          ztmp=dcmplx(0d0,-1d0)
          ztmp=ztmp**l*coulphase(-1d0/kappa_p,l)
          ztmp=ztmp*dsqrt(2*l+1d0)
          zsumm=dcmplx(0d0,0d0)
          do m=-l,l
             mdif=m-m_t(i0)
             mdifa=iabs(mdif)
             if(mdif.lt.0)then
               factor=(-1d0)**mdif
             else
               factor=1d0
             endif

             zex=dcmplx(0d0,0d0)
             do ir=1,nb
                btmp=b(ir)
                wtmp=wb(ir)

                do i=nneg+1,nmax_t(0)
                   
                   zintr(i-nneg)=amp(istar(i-l,l,m),ir)/dsqrt(dis(i-nneg)) 
!                   zintr(i-nneg)=amp(1,istar(i-l,l,m),ir)/dsqrt(dis(i-nneg)) 
          
                enddo

                call zintrpl(nbin,bin(1:nbin),zintr(1:nbin),kappa_p,zamp)

                if(dabs(ptr).lt.1d-4)then
                  if(mdifa==0)then
                    zex=zex+factor*wtmp*btmp*zamp
                  else
                    zex=dcmplx(0d0,0d0)
                  endif
                else
                    zex=zex+factor*wtmp*btmp*bessel_jn(mdifa,btmp*ptr)*zamp
                endif
  
             enddo 

            rylmp=0d0
            do l1=iabs(m),l
               l2=l-l1
               rylmp=rylmp+(-1d0)**l2*kappa**l1*v**l2/sqrt(fac(2*l1)*fac(2*l2))*cleba(l1,m,l2,0,l,m)*rylm(l1,m,costh)*(-1d0)**l1
            enddo
            rylmp=rylmp/kappa_p**l*sqrt(fac(2*l)) 

            res=atan2(-ktr*dsin(phi),pt(ip)-ktr*dcos(phi))
            zsumm=zsumm+rylmp*zex*dcmplx(0d0,1d0)**(m)*dcmplx(dcos(m*phi),dsin(m*phi))*dcmplx(dcos(m*res),-dsin(m*res))

          enddo
          amp_x(ip)=amp_x(ip)+ztmp*zsumm
        enddo
        amp_x(ip)=amp_x(ip)*sqrt(kf*ki/4d0/pi*kappa)/kappa_p
        tdcs_x_av(ip)=tdcs_x_av(ip)+amp_x(ip)*conjg(amp_x(ip))
        tdcs_coh_av(ip)=tdcs_coh_av(ip)+(amp_x(ip)+amp_d(ip,iphi))*conjg(amp_x(ip)+amp_d(ip,iphi))
     enddo
     tdcs_x_av(ip)=tdcs_x_av(ip)
     tdcs_coh_av(ip)=tdcs_coh_av(ip)
!   enddo
endif

10 continue

!Printing:

   if(keyp==1)then
        write(31,11)dble(it),tdcs_x_av(ip),tdcs_d_av(ip),tdcs_coh_av(ip),en_t(0)%enl_t(n)*hr
   endif

11 format(5000000es18.8)

end subroutine get_tdcs 



subroutine get_tdcs_phi_av3(keyp,kf,ki,n,it,tdcs_x_av,tdcs_d_av,tdcs_coh_av)
   use precision_type
   use data_base
   use flogs,only:fac
   use acc_routines,only:rylm
   implicit none
   integer::  ip,keyp,it,l,m,i,j,n,mdif,mdifa,ir,l1,l2,jc,jl,jr,iphi,phi_fix
   complex(kind=id):: ztmp,coulphase,zsumm,zex,zamp 
   complex(kind=id),dimension(1:nbin,0:lmax_t,-lmax_t:lmax_t):: zp
   complex(kind=id),dimension(0:npt):: amp_x,tdcs_x_av,tdcs_d_av,tdcs_coh_av 
   complex(kind=id),dimension(0:npt,0:360):: amp_d
   complex(kind=id),dimension(1:nbin):: zintr
   real(kind=id):: kappa,ovrlp,kf,ki,costh,kappa_p,btmp,wtmp,factor,bessj,rylmp,ptr,ktr,phi,res

   kappa=sqrt(2._id*en_t(0)%enl_t(n))
   ovrlp=1d0/dis(n-nmax_t(0)+nbin)

!phi_fix=270
!The direct scattering amplitude:

tdcs_x_av=0._id;tdcs_d_av=0._id;tdcs_coh_av=0._id

if(it.gt.180)then

   do ip=0,npt
     tdcs_d_av(ip)=zero  
     do iphi=0,360,5!phi_fix,phi_fix!360
        phi=dble(iphi)*pi/180d0 
        costh=dcos((it-180)*pi/180d0)
        if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
        amp_d(ip,iphi)=dcmplx(0d0,0d0)
        do l=0,lmax_t
          ztmp=dcmplx(0d0,-1d0)
          ztmp=ztmp**l*coulphase(-1d0/kappa,l)
          ztmp=ztmp*dsqrt(2*l+1d0)
          zsumm=dcmplx(0d0,0d0)
          do m=-l,l
            j=istar(n-l,l,m)+num_t
            zsumm=zsumm+rylm(l,m,costh)*zap(j,ip)*dcmplx(0d0,1d0)**(m)*exp(cu*m*phi)!*dcmplx(dcos(m*phi),dsin(m*phi))!*(-1d0)**m
          enddo
          amp_d(ip,iphi)=amp_d(ip,iphi)+ztmp*zsumm*(-1d0)**l
        enddo
        amp_d(ip,iphi)=amp_d(ip,iphi)*sqrt(kf*ki*ovrlp/4d0/pi/kappa)
        tdcs_d_av(ip)=tdcs_d_av(ip)+amp_d(ip,iphi)*conjg(amp_d(ip,iphi))
     enddo
     tdcs_d_av(ip)=tdcs_d_av(ip)*pi/180d0*5
   enddo
endif
if(it.le.180)then
   do ip=0,npt
     tdcs_d_av(ip)=zero  
     do iphi=0,360,5!phi_fix,phi_fix!360
        phi=dble(iphi)*pi/180d0 
        costh=dcos((it)*pi/180d0)
        if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
        amp_d(ip,iphi)=dcmplx(0d0,0d0)
        do l=0,lmax_t
          ztmp=dcmplx(0d0,-1d0)
          ztmp=ztmp**l*coulphase(-1d0/kappa,l)
          ztmp=ztmp*dsqrt(2*l+1d0)
          zsumm=dcmplx(0d0,0d0)
          do m=-l,l
            j=istar(n-l,l,m)+num_t
            zsumm=zsumm+rylm(l,m,costh)*zap(j,ip)*dcmplx(0d0,1d0)**(m)*exp(cu*m*phi)!*dcmplx(dcos(m*phi),dsin(m*phi))
          enddo
          amp_d(ip,iphi)=amp_d(ip,iphi)+ztmp*zsumm
        enddo
        amp_d(ip,iphi)=amp_d(ip,iphi)*sqrt(kf*ki*ovrlp/4d0/pi/kappa)
        tdcs_d_av(ip)=tdcs_d_av(ip)+amp_d(ip,iphi)*conjg(amp_d(ip,iphi))
     enddo
     tdcs_d_av(ip)=tdcs_d_av(ip)*pi/180d0*5
   enddo
endif

!return

!Rearrangement:

if(it.le.180)then
   do ip=0,npt
     ktr=kappa*dsin((it)*pi/180d0)
!     ptr=pt(ip)
     kappa_p=sqrt(kappa**2+v**2-2*kappa*v*dcos((it)*pi/180d0))

     costh=dcos((it)*pi/180d0)
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     tdcs_x_av(ip)=zero  
     tdcs_coh_av(ip)=zero  
     do iphi=0,360!phi_fix,phi_fix !360
        phi=dble(iphi)*pi/180d0 
        ptr=dsqrt(pt(ip)**2-2d0*pt(ip)*ktr*dcos(phi)+ktr**2)!dble(pt(ip)-ktr)
        amp_x(ip)=dcmplx(0d0,0d0)
        do l=0,lmax_t
          ztmp=dcmplx(0d0,-1d0)
          ztmp=ztmp**l*coulphase(-1d0/kappa_p,l)
          ztmp=ztmp*dsqrt(2*l+1d0)
          zsumm=dcmplx(0d0,0d0)
          do m=-l,l
             mdif=m-m_t(i0)
             mdifa=iabs(mdif)
             if(mdif.lt.0)then
               factor=(-1d0)**mdif
             else
               factor=1d0
             endif

             zex=dcmplx(0d0,0d0)
             do ir=1,nb
                btmp=b(ir)
                wtmp=wb(ir)

                do i=nneg+1,nmax_t(0)
                   
                   zintr(i-nneg)=amp(istar(i-l,l,m),ir)/dsqrt(dis(i-nneg)) 
!                   zintr(i-nneg)=amp(1,istar(i-l,l,m),ir)/dsqrt(dis(i-nneg)) 
          
                enddo

                call zintrpl(nbin,bin(1:nbin),zintr(1:nbin),kappa_p,zamp)

                if(dabs(ptr).lt.1d-4)then
                  if(mdifa==0)then
                    zex=zex+factor*wtmp*btmp*zamp
                  else
                    zex=dcmplx(0d0,0d0)
                  endif
                else
                    zex=zex+factor*wtmp*btmp*bessel_jn(mdifa,btmp*ptr)*zamp
                endif
  
             enddo 

            rylmp=0d0
            do l1=iabs(m),l
               l2=l-l1
               rylmp=rylmp+(-1d0)**l2*kappa**l1*v**l2/sqrt(fac(2*l1)*fac(2*l2))*cleba(l1,m,l2,0,l,m)*rylm(l1,m,costh)
            enddo
            rylmp=rylmp/kappa_p**l*sqrt(fac(2*l)) 

            res=atan2(-ktr*dsin(phi),pt(ip)-ktr*dcos(phi))
            zsumm=zsumm+rylmp*zex*dcmplx(0d0,1d0)**(m)*dcmplx(dcos(m*phi),dsin(m*phi))*dcmplx(dcos(m*res),-dsin(m*res))

          enddo
          amp_x(ip)=amp_x(ip)+ztmp*zsumm
        enddo
        amp_x(ip)=amp_x(ip)*sqrt(kf*ki/4d0/pi*kappa)/kappa_p
        tdcs_x_av(ip)=tdcs_x_av(ip)+amp_x(ip)*conjg(amp_x(ip))
        tdcs_coh_av(ip)=tdcs_coh_av(ip)+(amp_x(ip)+amp_d(ip,iphi))*conjg(amp_x(ip)+amp_d(ip,iphi))
     enddo
     tdcs_x_av(ip)=tdcs_x_av(ip)*pi/180d0
     tdcs_coh_av(ip)=tdcs_coh_av(ip)*pi/180d0
   enddo
endif



if(it.gt.180)then
   do ip=0,npt
     ktr=kappa*dsin((it)*pi/180d0)
!     ptr=dabs(pt(ip)-ktr)
     kappa_p=sqrt(kappa**2+v**2-2*kappa*v*dcos((it)*pi/180d0))

     costh=dcos((it-180)*pi/180d0)
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     tdcs_x_av(ip)=zero  
     tdcs_coh_av(ip)=zero  
     do iphi=0,360!phi_fix,phi_fix!0,360
        phi=dble(iphi)*pi/180d0 
        ptr=dsqrt(pt(ip)**2-2d0*pt(ip)*ktr*dcos(phi)+ktr**2)!dble(pt(ip)-ktr)
        amp_x(ip)=dcmplx(0d0,0d0)
        do l=0,lmax_t
          ztmp=dcmplx(0d0,-1d0)
          ztmp=ztmp**l*coulphase(-1d0/kappa_p,l)
          ztmp=ztmp*dsqrt(2*l+1d0)
          zsumm=dcmplx(0d0,0d0)
          do m=-l,l
             mdif=m-m_t(i0)
             mdifa=iabs(mdif)
             if(mdif.lt.0)then
               factor=(-1d0)**mdif
             else
               factor=1d0
             endif

             zex=dcmplx(0d0,0d0)
             do ir=1,nb
                btmp=b(ir)
                wtmp=wb(ir)

                do i=nneg+1,nmax_t(0)
                   
                   zintr(i-nneg)=amp(istar(i-l,l,m),ir)/dsqrt(dis(i-nneg)) 
!                   zintr(i-nneg)=amp(1,istar(i-l,l,m),ir)/dsqrt(dis(i-nneg)) 
          
                enddo

                call zintrpl(nbin,bin(1:nbin),zintr(1:nbin),kappa_p,zamp)

                if(dabs(ptr).lt.1d-4)then
                  if(mdifa==0)then
                    zex=zex+factor*wtmp*btmp*zamp
                  else
                    zex=dcmplx(0d0,0d0)
                  endif
                else
                    zex=zex+factor*wtmp*btmp*bessel_jn(mdifa,btmp*ptr)*zamp
                endif
  
             enddo 

            rylmp=0d0
            do l1=iabs(m),l
               l2=l-l1
               rylmp=rylmp+(-1d0)**l2*kappa**l1*v**l2/sqrt(fac(2*l1)*fac(2*l2))*cleba(l1,m,l2,0,l,m)*rylm(l1,m,costh)*(-1d0)**l1
            enddo
            rylmp=rylmp/kappa_p**l*sqrt(fac(2*l)) 

            res=atan2(-ktr*dsin(phi),pt(ip)-ktr*dcos(phi))
            zsumm=zsumm+rylmp*zex*dcmplx(0d0,1d0)**(m)*dcmplx(dcos(m*phi),dsin(m*phi))*dcmplx(dcos(m*res),-dsin(m*res))

          enddo
          amp_x(ip)=amp_x(ip)+ztmp*zsumm
        enddo
        amp_x(ip)=amp_x(ip)*sqrt(kf*ki/4d0/pi*kappa)/kappa_p
        tdcs_x_av(ip)=tdcs_x_av(ip)+amp_x(ip)*conjg(amp_x(ip))
        tdcs_coh_av(ip)=tdcs_coh_av(ip)+(amp_x(ip)+amp_d(ip,iphi))*conjg(amp_x(ip)+amp_d(ip,iphi))
     enddo
     tdcs_x_av(ip)=tdcs_x_av(ip)*pi/180d0
     tdcs_coh_av(ip)=tdcs_coh_av(ip)*pi/180d0
   enddo
endif


!Printing:

   if(keyp==1)then
      open(31,file='TDCS_AV2',status='replace')
      do ip=0,npt
        write(31,11)dble(ip),tdcs_x_av(ip),tdcs_d_av(ip),tdcs_coh_av(ip),en_t(0)%enl_t(n)*hr
      enddo
      close(31)
   endif

11 format(5000000es18.8)

end subroutine get_tdcs_phi_av3 


subroutine zintrpl(n,x,y,arg,res)
   use precision_type
   use data_base
   implicit none
   integer::n
   real(kind=id)::arg,rev,imv
   real(kind=id),dimension(n)::x,ry,iy
   complex(kind=id),dimension(n)::y
   complex(kind=id)::res

   ry(1:n)=dble(y(1:n))
   iy(1:n)=dimag(y(1:n))

   call intrpl(n,x,ry,1,arg,rev)
   call intrpl(n,x,iy,1,arg,imv)

   res=dcmplx(rev,imv)

   return
end subroutine zintrpl



subroutine solveRK_num_1c_dev( rho,prob,zamo )
   use precision_type
   use data_base,only:num_t,vleft,vmid,vright,nstep,step,i0,unitarity2,cu,nglag,ngleg,u,glag,wig_d,lmax_t,mmax_t,fdir,idir,jdir,k_t,inl_t,zam,zout,zdadz,dyt,dym,dymm
   use acc_routines,only:direct_num_t
   implicit none
   integer::ist_f,i,j,ist_i,k
   integer:: fin,iin,finl,iinl
   real(kind=id),intent(in) :: rho
   real(kind=id),intent(inout),dimension(1:num_t) :: prob
!   complex(kind=id),dimension(1:num_t) :: zam,zout,zdadz,dyt,dym,dymm
!   complex(kind=id) :: tmp 
   complex(kind=id),intent(out),dimension(1:num_t) :: zamo
   real(kind=id) :: z,znew,h,zmid,bigr,cosang,hh,h6
   real(kind=id) :: kf,ki
   integer :: t1, t2, dt, count_rate, count_max
   real(kind=id) ::  secs

!    open(220,file='solution',status='unknown')
!    write(220,*)
!    write(220,*)

allocate(zam(1:num_t))
allocate(zout(1:num_t))
allocate(zdadz(1:num_t))
allocate(dyt(1:num_t))
allocate(dym(1:num_t))
allocate(dymm(1:num_t))

print*,'I am here'
   call system_clock(count_max=count_max, count_rate=count_rate)
      call system_clock(t1)

!$acc data present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t) copyout(cosang,bigr,z,znew,zmid,prob) copyin(rho)  !copyout(prob)
!$acc kernels   !present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t) copyin(rho) copyout(prob,cosang,bigr,z) 
!$acc loop independent 
do ist_f=1,num_t
   zam(ist_f)=cmplx(0,0,id)
   zdadz(ist_f)=cmplx(0,0,id)
   zout(ist_f)=cmplx(0,0,id)
enddo
!$acc end loop 
zam(i0)=cmplx(1,0,id)
z=step(-nstep)
bigr=sqrt(rho*rho+z*z)
cosang=z/bigr
!$acc loop independent
do j=1,jdir
  vright(fdir(j),idir(j))=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(z*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(z*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
  vright(idir(j),fdir(j))=conjg(vright(fdir(j),idir(j)))
enddo
!$acc end loop 


!$acc loop
do j=1,num_t
  vright(j,j)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
enddo
!$acc end loop 

!$acc end kernels




do i=-nstep,nstep-1

!$acc kernels

   z=step(i)



   znew=step(i+1)
   zmid=(z+znew)/real(2,id)
   h=znew-z
   hh=h/2._id
   h6=h/6._id


!$acc loop 
   do k=1,num_t 
     vleft(k,1:num_t)=vright(k,1:num_t)
     zdadz(k)=-cu*sum(vleft(k,1:num_t)*zam(1:num_t))
   enddo
!$acc end loop 


   bigr=sqrt(rho*rho+zmid*zmid)
   cosang=zmid/bigr

!$acc loop independent
   do j=1,jdir
     vmid(fdir(j),idir(j))=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(zmid*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(zmid*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
     vmid(idir(j),fdir(j))=conjg(vmid(fdir(j),idir(j)))
   enddo
!$acc end loop 
!$acc loop 
   do j=1,num_t
     vmid(j,j)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
   enddo
!$acc end loop 



   bigr=sqrt(rho*rho+znew*znew)
   cosang=znew/bigr

!$acc loop independent 
   do j=1,jdir
     vright(fdir(j),idir(j))=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(znew*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(znew*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
     vright(idir(j),fdir(j))=conjg(vright(fdir(j),idir(j)))
   enddo
!$acc end loop 
!$acc loop 
   do j=1,num_t
     vright(j,j)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
   enddo
!$acc end loop 



!$acc loop 
   do j=1,num_t 
     dyt(j)=-cu*sum(vmid(j,1:num_t)*(zam(1:num_t)+hh*zdadz(1:num_t)))
   enddo
!$acc end loop
!$acc loop 
   do j=1,num_t 
     dym(j)=-cu*sum(vmid(j,1:num_t)*(zam(1:num_t)+hh*dyt(1:num_t)))
   enddo
!$acc end loop 
!$acc loop 
   do j=1,num_t 
     dymm(j)=dyt(j)+dym(j)
     dyt(j)=-cu*sum(vright(j,1:num_t)*(zam(1:num_t)+h*dym(1:num_t)))
     zout(j)=zam(j)+h6*(zdadz(j)+dyt(j)+2._id*dymm(j))
   enddo
!$acc end loop 
!$acc loop 
   do j=1,num_t 
     zam(j)=zout(j)
   enddo
!$acc end loop 

!$acc end kernels 
!!$acc update host(zam,z,dymm,bigr)
!print*,z,bigr,zam(1:3),dymm(1:3)
!stop

enddo




!$acc kernels 

!
!!$omp do schedule(guided)
!    do j=1,num_t 
!      unitarity2=unitarity2+abs(zam(j))**2
!    enddo
!!$omp end do

unitarity2=sum(abs(zam(1:num_t))**2)

    zam(i0)=zam(i0)-1._id


!$acc loop 
    do ist_f=1,num_t
       prob(ist_f)=zam(ist_f)*conjg(zam(ist_f))*rho
    enddo
!$acc end loop 
!!$acc end parallel

!$acc end kernels 



!$acc end data
!$acc update host(zam,unitarity2)


zamo(:)=zam(:)

    write(1298,*)rho,unitarity2

      call system_clock(t2)
      dt = t2-t1
      secs = real(dt)/real(count_rate)

      open(1299,file='time',status='unknown')
      write(1299,"('number of channels is ',i5,' wall clock time is ',f12.2,' seconds')") num_t,secs
      open(1300,file='timing',status='unknown')
      write(1300,"(i5,f12.2)") num_t,secs

 30   format(i20,10000g14.6)

deallocate(zam)
deallocate(zout)
deallocate(zdadz)
deallocate(dyt)
deallocate(dym)
deallocate(dymm)
end subroutine solveRK_num_1c_dev


subroutine solveRK_num_1c( rho,prob,zamo )
   use precision_type
   use data_base,only:num_t,vleft,vmid,vright,nstep,step,i0,unitarity2,cu,nglag,ngleg,u,glag,wig_d,lmax_t,mmax_t,fdir,idir,jdir,k_t,inl_t,zam,zout,zdadz,dyt,dym,dymm
   use acc_routines,only:direct_num_t
   implicit none
   integer::ist_f,i,j,ist_i,k
   integer:: fin,iin,finl,iinl
   real(kind=id),intent(in) :: rho
   real(kind=id),intent(inout),dimension(1:num_t) :: prob
!   complex(kind=id),dimension(1:num_t) :: zam,zout,zdadz,dyt,dym,dymm
!   complex(kind=id) :: tmp 
   complex(kind=id),intent(out),dimension(1:num_t) :: zamo
   real(kind=id) :: z,znew,h,zmid,bigr,cosang,hh,h6
   real(kind=id) :: kf,ki
   integer :: t1, t2, dt, count_rate, count_max
   real(kind=id) ::  secs

!    open(220,file='solution',status='unknown')
!    write(220,*)
!    write(220,*)

allocate(zam(1:num_t))
allocate(zout(1:num_t))
allocate(zdadz(1:num_t))
allocate(dyt(1:num_t))
allocate(dym(1:num_t))
allocate(dymm(1:num_t))

print*,'I am here'


   call system_clock(count_max=count_max, count_rate=count_rate)
      call system_clock(t1)

!!$acc data present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t) copyout(cosang,bigr,z) copyin(rho) copyout(prob)
!!$acc parallel  
!$acc kernels present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t) copyin(rho) copyout(prob,cosang,bigr,z) 
!!$omp parallel default(none) shared(num_t,dyt,zdadz,vleft,vmid,vright,h,hh,h6,dym,dymm,zam,zout,bigr,znew,zmid,z,k_t,inl_t,cosang,fdir,idir,jdir,prob,rho, i0,step,nstep) reduction(+:unitarity2) 
!!$omp do schedule(guided)
!!$acc loop 
do ist_f=1,num_t
   zam(ist_f)=cmplx(0,0,id)
   zdadz(ist_f)=cmplx(0,0,id)
   zout(ist_f)=cmplx(0,0,id)
enddo
!!$omp end do
!!$omp single
zam(i0)=cmplx(1,0,id)
!!$acc update device(zam)
z=step(-nstep)
bigr=sqrt(rho*rho+z*z)
cosang=z/bigr
!!$omp end single

!!$omp do schedule(guided)
!$acc loop independent
do j=1,jdir
  vright(fdir(j),idir(j))=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(z*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(z*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
  vright(idir(j),fdir(j))=conjg(vright(fdir(j),idir(j)))
enddo
!!$omp end do


!!$omp do schedule(guided)
!$acc loop
do j=1,num_t
  vright(j,j)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
enddo
!!$omp end do



do i=-nstep,nstep-1

!i=0

!!$omp single
   z=step(i)



   znew=step(i+1)
   zmid=(z+znew)/real(2,id)
   h=znew-z
   hh=h/2._id
   h6=h/6._id
!!$omp end single


!!$omp do schedule(guided)
!$acc loop 
   do k=1,num_t 
     vleft(k,1:num_t)=vright(k,1:num_t)
     zdadz(k)=-cu*sum(vleft(k,1:num_t)*zam(1:num_t))
   enddo
!$acc end loop 
!!$omp end do

!!$acc update host(vleft,zdadz,zam)

!!$omp single
   bigr=sqrt(rho*rho+zmid*zmid)
   cosang=zmid/bigr
!!$omp end single

!!$omp do schedule(guided)
!$acc loop independent 
   do j=1,jdir
     vmid(fdir(j),idir(j))=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(zmid*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(zmid*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
     vmid(idir(j),fdir(j))=conjg(vmid(fdir(j),idir(j)))
   enddo
!$acc end loop 
!!$omp end do
!!$omp do schedule(guided)
!$acc loop 
   do j=1,num_t
     vmid(j,j)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
   enddo
!$acc end loop 
!!$omp end do



!!$omp single
   bigr=sqrt(rho*rho+znew*znew)
   cosang=znew/bigr
!!$omp end single 

!!$omp do schedule(guided)
!$acc loop independent 
   do j=1,jdir
     vright(fdir(j),idir(j))=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(znew*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(znew*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
     vright(idir(j),fdir(j))=conjg(vright(fdir(j),idir(j)))
   enddo
!$acc end loop 
!!$omp end do
!!$omp do schedule(guided)
!$acc loop 
   do j=1,num_t
     vright(j,j)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
   enddo
!$acc end loop 
!!$omp end do
!!$omp do schedule(guided)
!$acc loop 
   do j=1,num_t 
     dyt(j)=-cu*sum(vmid(j,1:num_t)*(zam(1:num_t)+hh*zdadz(1:num_t)))
   enddo
!$acc end loop
!!$omp end do
!!$omp do schedule(guided)
!$acc loop 
   do j=1,num_t 
     dym(j)=-cu*sum(vmid(j,1:num_t)*(zam(1:num_t)+hh*dyt(1:num_t)))
   enddo
!$acc end loop 
!!$omp end do
!!$omp do schedule(guided)
!$acc loop 
   do j=1,num_t 
     dymm(j)=dyt(j)+dym(j)
     dyt(j)=-cu*sum(vright(j,1:num_t)*(zam(1:num_t)+h*dym(1:num_t)))
     zout(j)=zam(j)+h6*(zdadz(j)+dyt(j)+2._id*dymm(j))
   enddo
!$acc end loop 
!!$omp end do 
!!$omp do schedule(guided)
!$acc loop 
   do j=1,num_t 
     zam(j)=zout(j)
   enddo
!$acc end loop 
!!$omp end do 



enddo


!!$omp single 
!    unitarity2=0.0_id
!!$omp end single
!
!!$omp do schedule(guided)
!    do j=1,num_t 
!      unitarity2=unitarity2+abs(zam(j))**2
!    enddo
!!$omp end do

!$omp single 
unitarity2=sum(abs(zam(1:num_t))**2)
!$omp end single

!!$omp single 
    zam(i0)=zam(i0)-1._id
!!$omp end single


!!$omp do schedule(guided)
!$acc loop 
    do ist_f=1,num_t
       prob(ist_f)=zam(ist_f)*conjg(zam(ist_f))*rho
    enddo
!$acc end loop 
!!$omp end do
!!$omp end parallel
!!$acc end parallel

!$acc end kernels 
!!$acc end data
!$acc update host(zam,unitarity2)

!!$acc update host(zam,unitarity2)
!!$acc update host(prob)

zamo(:)=zam(:)

    write(1298,*)rho,unitarity2

      call system_clock(t2)
      dt = t2-t1
      secs = real(dt)/real(count_rate)
      write(1299,"(' wall clock time is ',f12.2,' seconds')") &
              secs

 30   format(i20,10000g14.6)

end subroutine solveRK_num_1c

subroutine solveRK_num_2c_cpu( rho )
   use precision_type
   use data_base,only:num_t,vleft,vmid,vright,vmatl,nstep,step,key,i0,unitarity2,cu,nglag,ngleg,u,glag,gleg,wig_d,lmax_t,mmax_t,wigarr,larr,marr1,marr2,&
&   nlnum_t,inl_t,inl_p,nst_t,lst_t,basis_a,basis_b,basis_bres,jdir,idir,fdir,istar,k_t,k_p,l_p,lob,upb,m_t,zero,zam,zout,zdadz,dyt,dym,dymm,plma,plmb,&
&   lar,mar,lmar,v,jint,n_t,l_t,m_t,vLtmp,vRtmp,vover,z,znew,h,zmid,bigr,cosang,ang,hh,h6
!   use acc_routines,only:direct_num_t,direct_num_p,get_plmij
   use acc_routines
   use acc_wig
   implicit none
   integer::ist_f,i,j,ist_i,info,k,l,m1,m2,n,nf,lf,mf,ni,li,mi,lm
   real(kind=id),intent(in) :: rho
   real(kind=id) :: fct,wigd1,wigd2,arg_bess
!   real(kind=id) :: z,znew,h,zmid,bigr,cosang,ang,hh,h6
   integer :: t1, t2, dt, count_rate, count_max
   real(kind=id) ::  secs

!    open(220,file='solution',status='unknown')
!    write(220,*)
!    write(220,*)

print*,'I am here in cpu'

   call system_clock(count_max=count_max, count_rate=count_rate)
      call system_clock(t1)

!$omp parallel default(shared)

!do i=-nstep,nstep
   z=step(1)
   bigr=sqrt(rho*rho+z*z)
   cosang=z/bigr
   ang=acos(z/bigr)


!$omp do schedule(guided)
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$omp end do

!$omp do collapse(3) schedule(guided)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$omp end do

!$omp do collapse(3) schedule(guided)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$omp end do

!$omp do schedule(guided)
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$omp end do

!$omp do collapse(3) schedule(guided)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$omp end do

!$omp do collapse(2) schedule(guided)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_dev_pt2(j,k,bigr,z,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$omp end do

!$omp do collapse(2) schedule(guided)
do j=1,num_t
  do k=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j, k+num_t),vright(j,k+num_t))
  enddo
enddo
!$omp end do
!enddo
!$omp end parallel


      call system_clock(t2)
      dt = t2-t1
      secs = real(dt)/real(count_rate)

      open(1299,file='time',status='unknown')
      write(1299,"('number of channels is ',i5,' wall clock time is ',f12.2,' seconds')") num_t,secs
      open(1300,file='timing',status='unknown')
      write(1300,"(i5,3es15.6)") num_t,secs,sum(vright(:,:))

print*, 'cpu',z,sum(vright(:,:))

 30   format(i20,10000g14.6)

end subroutine solveRK_num_2c_cpu

subroutine mat_2c_cpu( rho )
   use precision_type
   use data_base,only:num_t,vleft,vmid,vright,vmatl,nstep,step,key,i0,unitarity2,cu,nglag,ngleg,u,glag,gleg,wig_d,lmax_t,mmax_t,wigarr,larr,marr1,marr2,&
&   nlnum_t,inl_t,inl_p,nst_t,lst_t,basis_a,basis_b,basis_bres,jdir,idir,fdir,istar,k_t,k_p,l_p,lob,upb,m_t,zero,zam,zout,zdadz,dyt,dym,dymm,plma,plmb,&
&   lar,mar,lmar,v,jint,n_t,l_t,m_t,vLtmp,vRtmp,vover,z,znew,h,zmid,bigr,cosang,ang,hh,h6,zij1,zij2,wglag,wgleg
!   use acc_routines,only:direct_num_t,direct_num_p,get_plmij
   use acc_routines
   use acc_wig
   implicit none
   integer::ist_f,i,j,ist_i,info,k,l,m1,m2,n,nf,lf,mf,ni,li,mi,lm
   real(kind=id),intent(in) :: rho
   real(kind=id) :: fct,wigd1,wigd2,arg_bess,arg
!   real(kind=id) :: z,znew,h,zmid,bigr,cosang,ang,hh,h6
   integer :: t1, t2, dt, count_rate, count_max
   real(kind=id) ::  secs

!    open(220,file='solution',status='unknown')
!    write(220,*)
!    write(220,*)
  allocate(zij1(1:nglag,1:ngleg))
  allocate(zij2(1:nglag,1:ngleg))

print*,'I am here in cpu'

   call system_clock(count_max=count_max, count_rate=count_rate)
      call system_clock(t1)
!$omp parallel default(shared) private(arg)

!do i=-nstep,nstep
   z=step(1)
   bigr=sqrt(rho*rho+z*z)
   cosang=z/bigr
   ang=acos(z/bigr)


!$omp do schedule(guided)
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$omp end do

!$omp do collapse(2) schedule(guided)
    do j=1,nglag
      do k=1,ngleg
        arg=0.5_id*v*z*u(j)*gleg(k)
        zij1(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)*u(j)-gleg(k)*gleg(k))
        zij2(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)-gleg(k))
      enddo
    enddo
!$omp end do

!$omp do collapse(3) schedule(guided)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$omp end do

!$omp do collapse(3) schedule(guided)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$omp end do


!$omp do collapse(3) schedule(guided)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$omp end do


!$omp do collapse(2) schedule(guided)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_sum_pt(j,k,bigr,z,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$omp end do


!$omp do schedule(guided)
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$omp end do

!$omp do collapse(2) schedule(guided)
do j=1,num_t
  do k=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j, k+num_t),vright(j,k+num_t))
  enddo
enddo
!$omp end do
!enddo
!$omp end parallel

      call system_clock(t2)
      dt = t2-t1
      secs = real(dt)/real(count_rate)

      open(1299,file='time',status='unknown')
      write(1299,"('number of channels is ',i5,' wall clock time is ',f12.2,' seconds')") num_t,secs
      open(1300,file='timing',status='unknown')
      write(1300,"(i5,5es15.6)") num_t,secs,sum(vright(:,:)),sum(vmatL(:,:))

print*, 'cpu',z,sum(vright(:,:))

deallocate(zij1,zij2)

 30   format(i20,10000g14.6)

end subroutine mat_2c_cpu

subroutine mat_2c_gpu_parallel( rho )
   use precision_type
   use data_base,only:num_t,vleft,vmid,vright,vmatl,nstep,step,key,i0,unitarity2,cu,nglag,ngleg,u,glag,gleg,wig_d,lmax_t,mmax_t,wigarr,larr,marr1,marr2,&
&   nlnum_t,inl_t,inl_p,nst_t,lst_t,basis_a,basis_b,basis_bres,jdir,idir,fdir,istar,k_t,k_p,l_p,lob,upb,m_t,zero,zam,zout,zdadz,dyt,dym,dymm,plma,plmb,&
&   lar,mar,lmar,v,jint,n_t,l_t,m_t,vLtmp,vRtmp,vover,z,znew,h,zmid,bigr,cosang,ang,hh,h6,zij1,zij2,wglag,wgleg
   use acc_routines
   use acc_wig
   implicit none
   integer::ist_f,i,j,ist_i,info,k,l,m1,m2,n,nf,lf,mf,ni,li,mi,lm
   real(kind=id),intent(in) :: rho
   real(kind=id) :: fct,wigd1,wigd2,arg_bess,arg
!   real(kind=id) :: z,znew,h,zmid,bigr,cosang,ang,hh,h6
   integer :: t1, t2, dt, count_rate, count_max
   real(kind=id) ::  secs

!    open(220,file='solution',status='unknown')
!    write(220,*)
!    write(220,*)

  allocate(zij1(1:nglag,1:ngleg))
  allocate(zij2(1:nglag,1:ngleg))

print*,'I am here in gpu mat'


   call system_clock(count_max=count_max, count_rate=count_rate)
      call system_clock(t1)

!$acc data present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t,lmar,lar,mar,i0,nglag,ngleg,wglag,wgleg,gleg) copyin(rho,v)
!$acc serial
   z=step(1)
   bigr=sqrt(rho*rho+z*z)
   cosang=z/bigr
   ang=acos(z/bigr)
!$acc end serial

!$acc parallel loop
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
    do j=1,nglag
      do k=1,ngleg
        arg=0.5_id*v*z*u(j)*gleg(k)
        zij1(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)*u(j)-gleg(k)*gleg(k))
        zij2(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)-gleg(k))
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop


!$acc parallel loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop


!!$acc parallel loop gang vector_length(128) independent collapse(2)
!do j=1,num_t
!  do k=1,num_t
!    call get_LR_mat_sum_pt(j,k,bigr,z,vLtmp(j,k),vRtmp(j,k))
!  enddo
!enddo
!!$acc end loop


call get_LR_mat_sum_all_pt  !(bigr,z)

!$acc parallel loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop


!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j, k+num_t),vright(j,k+num_t))
  enddo
enddo
!$acc end loop
!$acc end data

      call system_clock(t2)
      dt = t2-t1
      secs = real(dt)/real(count_rate)

!$acc update host(z,vright,vmatL)


      open(1299,file='time',status='unknown')
      write(1299,"('number of channels is ',i5,' wall clock time is ',f12.2,' seconds')") num_t,secs
      open(1300,file='timing',status='unknown')
      write(1300,"(i5,5es15.6)") num_t,secs,sum(vright(:,:)),sum(vmatL(:,:))

print*, 'gpu', z,sum(vright(:,:))

deallocate(zij1,zij2)

 30   format(i20,10000g14.6)

end subroutine mat_2c_gpu_parallel



subroutine solveRK_num_2c_gpu_parallel( rho )
   use precision_type
   use data_base,only:num_t,vleft,vmid,vright,vmatl,nstep,step,key,i0,unitarity2,cu,nglag,ngleg,u,glag,gleg,wig_d,lmax_t,mmax_t,wigarr,larr,marr1,marr2,&
&   nlnum_t,inl_t,inl_p,nst_t,lst_t,basis_a,basis_b,basis_bres,jdir,idir,fdir,istar,k_t,k_p,l_p,lob,upb,m_t,zero,zam,zout,zdadz,dyt,dym,dymm,plma,plmb,&
&   lar,mar,lmar,v,jint,n_t,l_t,m_t,vLtmp,vRtmp,vover,z,znew,h,zmid,bigr,cosang,ang,hh,h6
!   use acc_routines,only:direct_num_t,direct_num_p,get_plmij
   use acc_routines
   use acc_wig
   implicit none
   integer::ist_f,i,j,ist_i,info,k,l,m1,m2,n,nf,lf,mf,ni,li,mi,lm
   real(kind=id),intent(in) :: rho
   real(kind=id) :: fct,wigd1,wigd2,arg_bess
!   real(kind=id) :: z,znew,h,zmid,bigr,cosang,ang,hh,h6
   integer :: t1, t2, dt, count_rate, count_max
   real(kind=id) ::  secs

!    open(220,file='solution',status='unknown')
!    write(220,*)
!    write(220,*)

print*,'I am here in gpu'

   call system_clock(count_max=count_max, count_rate=count_rate)
      call system_clock(t1)


!$acc data present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t,lmar,lar,mar,i0) copyin(rho,v)
!$acc serial
   z=step(1)
   bigr=sqrt(rho*rho+z*z)
   cosang=z/bigr
   ang=acos(z/bigr)
!$acc end serial

!$acc parallel loop
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_dev_pt2(j,k,bigr,z,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j, k+num_t),vright(j,k+num_t))
  enddo
enddo
!$acc end loop
!$acc end data

!$acc update host(z,vright)

      call system_clock(t2)
      dt = t2-t1
      secs = real(dt)/real(count_rate)

      open(1299,file='time',status='unknown')
      write(1299,"('number of channels is ',i5,' wall clock time is ',f12.2,' seconds')") num_t,secs
      open(1300,file='timing',status='unknown')
      write(1300,"(i5,3es15.6)") num_t,secs,sum(vright(:,:))

print*, 'gpu', z,sum(vright(:,:))

 30   format(i20,10000g14.6)

end subroutine solveRK_num_2c_gpu_parallel

subroutine solveRK_num_2c_gpu( rho )
   use precision_type
   use data_base,only:num_t,vleft,vmid,vright,vmatl,nstep,step,key,i0,unitarity2,cu,nglag,ngleg,u,glag,gleg,wig_d,lmax_t,mmax_t,wigarr,larr,marr1,marr2,&
&   nlnum_t,inl_t,inl_p,nst_t,lst_t,basis_a,basis_b,basis_bres,jdir,idir,fdir,istar,k_t,k_p,l_p,lob,upb,m_t,zero,zam,zout,zdadz,dyt,dym,dymm,plma,plmb,&
&   lar,mar,lmar,v,jint,n_t,l_t,m_t,vLtmp,vRtmp,vover,z,znew,h,zmid,bigr,cosang,ang,hh,h6
!   use acc_routines,only:direct_num_t,direct_num_p,get_plmij
   use acc_routines
   use acc_wig
   implicit none
   integer::ist_f,i,j,ist_i,info,k,l,m1,m2,n,nf,lf,mf,ni,li,mi,lm
   real(kind=id),intent(in) :: rho
   real(kind=id) :: fct,wigd1,wigd2,arg_bess
!   real(kind=id) :: z,znew,h,zmid,bigr,cosang,ang,hh,h6
   integer :: t1, t2, dt, count_rate, count_max
   real(kind=id) ::  secs

!    open(220,file='solution',status='unknown')
!    write(220,*)
!    write(220,*)

print*,'I am here in gpu'

   call system_clock(count_max=count_max, count_rate=count_rate)
      call system_clock(t1)


!$acc data present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t,lmar,lar,mar,i0) copyin(rho,v)
!$acc kernels
   z=step(1)
   bigr=sqrt(rho*rho+z*z)
   cosang=z/bigr
   ang=acos(z/bigr)

!$acc loop independent 
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop

!$acc loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop

!$acc loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc loop independent !collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_dev_pt2(j,k,bigr,z,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc loop independent !collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j, k+num_t),vright(j,k+num_t))
  enddo
enddo
!$acc end loop
!$acc end kernels

!$acc end data

!$acc update host(z,vright)

      call system_clock(t2)
      dt = t2-t1
      secs = real(dt)/real(count_rate)

      open(1299,file='time',status='unknown')
      write(1299,"('number of channels is ',i5,' wall clock time is ',f12.2,' seconds')") num_t,secs
      open(1300,file='timing',status='unknown')
      write(1300,"(i5,f12.2)") num_t,secs

print*, 'gpu', z,sum(vright(:,:))

 30   format(i20,10000g14.6)

end subroutine solveRK_num_2c_gpu

subroutine RK_1c_FBA_gpu( rho,prob,zamo )
   use openacc
   use precision_type
   use data_base,only:num_t,vleft,vmid,vright,vmatl,nstep,step,key,i0,unitarity2,cu,nglag,ngleg,u,glag,gleg,wig_d,lmax_t,mmax_t,wigarr,larr,marr1,marr2,&
&   nlnum_t,inl_t,inl_p,nst_t,lst_t,basis_a,basis_b,basis_bres,jdir,idir,fdir,istar,k_t,k_p,l_p,lob,upb,m_t,zero,zam,zout,zdadz,dyt,dym,dymm,plma,plmb,&
&   lar,mar,lmar,v,jint,n_t,l_t,m_t,vLtmp,vRtmp,vover,z,znew,h,zmid,bigr,cosang,ang,hh,h6,zij1,zij2,wglag,wgleg
!   use acc_routines,only:direct_num_t,direct_num_p,get_plmij
   use acc_routines
   use acc_wig
   implicit none
   integer::ist_f,i,j,ist_i,info,k,l,m1,m2,n,nf,lf,mf,ni,li,mi,lm
   real(kind=id),intent(in) :: rho
   real(kind=id),intent(inout),dimension(1:2*num_t) :: prob
   complex(kind=id),intent(out),dimension(1:2*num_t) :: zamo
!   complex(kind=id),dimension(1:num_t,1:num_t) :: vLtmp,vRtmp
!   complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t) :: vmatLtmp,vmatRtmp
   complex(kind=id),dimension(1:2*num_t*2*num_t) ::work
   complex(kind=id)::zrm1,zrm2,zlm1,zlm2
   real(kind=id) :: fct,wigd1,wigd2,arg_bess,arg
!   real(kind=id) :: z,znew,h,zmid,bigr,cosang,ang,hh,h6
   integer,dimension(1:2*num_t)::ipiv
   integer :: t1, t2, dt, count_rate, count_max
   real(kind=id) ::  secs

!    open(220,file='solution',status='unknown')
!    write(220,*)
!    write(220,*)

allocate(zam(1:2*num_t))
allocate(zout(1:2*num_t))
allocate(zdadz(1:2*num_t))
allocate(dyt(1:2*num_t))
allocate(dym(1:2*num_t))
allocate(dymm(1:2*num_t))
allocate(zij1(1:nglag,1:ngleg))
allocate(zij2(1:nglag,1:ngleg))

   call system_clock(count_max=count_max, count_rate=count_rate)
      call system_clock(t1)


!$acc data present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t,lmar,lar,mar,i0) copyout(prob) copyin(rho,v)
!$acc parallel loop independent 
do ist_f=1,2*num_t
   zam(ist_f)=cmplx(0,0,id)
   zdadz(ist_f)=cmplx(0,0,id)
   zout(ist_f)=cmplx(0,0,id)
   vright(ist_f,1)=cmplx(0,0,id)
   vmid(ist_f,1)=cmplx(0,0,id)
   vleft(ist_f,1)=cmplx(0,0,id)
enddo
!$acc end loop

!$acc serial
   zam(num_t+i0)=cmplx(1,0,id)
   z=step(-nstep)
   bigr=sqrt(rho*rho+z*z)
   cosang=z/bigr
   ang=acos(z/bigr)
!$acc end serial

!$acc parallel loop independent
do j=1,num_t
  vright(j+num_t,1)=direct_num_t(j,1,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(1))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(1))-k_t(inl_t(j)))),id)
enddo
!$acc end loop


do i=-nstep,nstep-1

!$acc serial

    z=step(i)
    znew=step(i+1)
    zmid=(z+znew)/real(2,id)
    h=znew-z
    hh=h/2._id
    h6=h/6._id

!$acc end serial

!$acc parallel loop
   do k=1,2*num_t
     vleft(k,1)=vright(k,1)
     zdadz(k)=-cu*vleft(k,1)
   enddo
!$acc end loop

!$acc serial
    bigr=sqrt(rho*rho+zmid*zmid)
    cosang=zmid/bigr
    ang=acos(zmid/bigr)
!$acc end serial

!$acc parallel loop independent
do j=1,num_t
  vmid(j+num_t,1)=direct_num_t(j,1,cosang,bigr)*cmplx(cos(zmid*(k_t(inl_t(1))-k_t(inl_t(j)))),sin(zmid*(k_t(inl_t(1))-k_t(inl_t(j)))),id)
enddo
!$acc end loop

!$acc serial
    bigr=sqrt(rho*rho+znew*znew)
    cosang=znew/bigr
    ang=acos(znew/bigr)
!$acc end serial

!$acc parallel loop independent
do j=1,num_t
  vright(j+num_t,1)=direct_num_t(j,1,cosang,bigr)*cmplx(cos(znew*(k_t(inl_t(1))-k_t(inl_t(j)))),sin(znew*(k_t(inl_t(1))-k_t(inl_t(j)))),id)
enddo
!$acc end loop

!$acc parallel loop
   do j=1,2*num_t 
     dymm(j)=-2._id*cu*vmid(j,1)
     dyt(j)=-cu*vright(j,1)
     zout(j)=zam(j)+h6*(zdadz(j)+dyt(j)+2._id*dymm(j))
   enddo
!$acc end loop

!$acc parallel loop
   do j=1,2*num_t 
     zam(j)=zout(j)
   enddo
!$acc end loop

enddo

print*,'FBA direct'

!$acc serial
zam(num_t+i0)=zam(num_t+i0)-1._id
!$acc end serial
!$acc parallel loop
    do ist_f=1,2*num_t
       prob(ist_f)=zam(ist_f)*conjg(zam(ist_f))*rho
    enddo
!$acc end loop

!$acc end data

!$acc update host(zam,vover)

zamo(:)=zam(:)

zam(num_t+i0)=zam(num_t+i0)+1._id

unitarity2=0.0_id

write(1298,*)rho,dble(unitarity2)

      call system_clock(t2)
      dt = t2-t1
      secs = real(dt)/real(count_rate)

      open(1299,file='time',status='unknown')
      write(1299,"('number of channels is ',i5,' wall clock time is ',f12.2,' seconds')") num_t,secs
      open(1300,file='timing',status='unknown')
      write(1300,"(i5,f12.2)") num_t,secs

 30   format(i20,10000g14.6)

deallocate(zam)
deallocate(zout)
deallocate(zdadz)
deallocate(dyt)
deallocate(dym)
deallocate(dymm)
deallocate(zij1)
deallocate(zij2)

end subroutine RK_1c_FBA_gpu


subroutine RK_2c_FBA_gpu( rho,prob,zamo )
   use openacc
   use precision_type
   use data_base,only:num_t,vleft,vmid,vright,vmatl,nstep,step,key,i0,unitarity2,cu,nglag,ngleg,u,glag,gleg,wig_d,lmax_t,mmax_t,wigarr,larr,marr1,marr2,&
&   nlnum_t,inl_t,inl_p,nst_t,lst_t,basis_a,basis_b,basis_bres,jdir,idir,fdir,istar,k_t,k_p,l_p,lob,upb,m_t,zero,zam,zout,zdadz,dyt,dym,dymm,plma,plmb,&
&   lar,mar,lmar,v,jint,n_t,l_t,m_t,vLtmp,vRtmp,vover,z,znew,h,zmid,bigr,cosang,ang,hh,h6,zij1,zij2,wglag,wgleg
!   use acc_routines,only:direct_num_t,direct_num_p,get_plmij
   use acc_routines
   use acc_wig
   implicit none
   integer::ist_f,i,j,ist_i,info,k,l,m1,m2,n,nf,lf,mf,ni,li,mi,lm
   real(kind=id),intent(in) :: rho
   real(kind=id),intent(inout),dimension(1:2*num_t) :: prob
   complex(kind=id),intent(out),dimension(1:2*num_t) :: zamo
!   complex(kind=id),dimension(1:num_t,1:num_t) :: vLtmp,vRtmp
!   complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t) :: vmatLtmp,vmatRtmp
   complex(kind=id),dimension(1:2*num_t*2*num_t) ::work
   complex(kind=id)::zrm1,zrm2,zlm1,zlm2
   real(kind=id) :: fct,wigd1,wigd2,arg_bess,arg
!   real(kind=id) :: z,znew,h,zmid,bigr,cosang,ang,hh,h6
   integer,dimension(1:2*num_t)::ipiv
   integer :: t1, t2, dt, count_rate, count_max
   real(kind=id) ::  secs

!    open(220,file='solution',status='unknown')
!    write(220,*)
!    write(220,*)

allocate(zam(1:2*num_t))
allocate(zout(1:2*num_t))
allocate(zdadz(1:2*num_t))
allocate(dyt(1:2*num_t))
allocate(dym(1:2*num_t))
allocate(dymm(1:2*num_t))
allocate(zij1(1:nglag,1:ngleg))
allocate(zij2(1:nglag,1:ngleg))

   call system_clock(count_max=count_max, count_rate=count_rate)
      call system_clock(t1)


!$acc data present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t,lmar,lar,mar,i0) copyout(prob) copyin(rho,v)
!$acc parallel loop independent 
do ist_f=1,2*num_t
   zam(ist_f)=cmplx(0,0,id)
   zdadz(ist_f)=cmplx(0,0,id)
   zout(ist_f)=cmplx(0,0,id)
enddo
!$acc end loop

!$acc serial
   zam(num_t+i0)=cmplx(1,0,id)
   z=step(-nstep)
   bigr=sqrt(rho*rho+z*z)
   cosang=z/bigr
   ang=acos(z/bigr)
!$acc end serial

!$acc parallel loop independent 
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
    do j=1,nglag
      do k=1,ngleg
        arg=0.5_id*v*z*u(j)*gleg(k)
        zij1(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)*u(j)-gleg(k)*gleg(k))
        zij2(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)-gleg(k))
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop


!$acc parallel loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc parallel loop independent
do j=1,num_t
    call get_LR_mat_sum_pt(j,1,bigr,z,vLtmp(j,1),vRtmp(j,1))
enddo
!$acc end loop

!$acc parallel loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop

!$acc parallel loop independent
do j=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(1),l_t(1),m_t(1),vmatL(j, 1),vright(j,1))
enddo
!$acc end loop



!$acc parallel loop independent
do j=1,num_t
  vright(j+num_t,1)=direct_num_t(j,1,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(1))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(1))-k_t(inl_t(j)))),id)
enddo
!$acc end loop


do i=-nstep,nstep-1

!$acc serial

    z=step(i)
    znew=step(i+1)
    zmid=(z+znew)/real(2,id)
    h=znew-z
    hh=h/2._id
    h6=h/6._id

!$acc end serial

!$acc parallel loop
   do k=1,2*num_t
     vleft(1:2*num_t,1)=vright(1:2*num_t,1)
     zdadz(k)=-cu*vleft(k,1)
   enddo
!$acc end loop

!$acc serial
    bigr=sqrt(rho*rho+zmid*zmid)
    cosang=zmid/bigr
    ang=acos(zmid/bigr)
!$acc end serial

!$acc parallel loop independent 
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
    do j=1,nglag
      do k=1,ngleg
        arg=0.5_id*v*zmid*u(j)*gleg(k)
        zij1(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)*u(j)-gleg(k)*gleg(k))
        zij2(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)-gleg(k))
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop


!$acc parallel loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc parallel loop independent
do j=1,num_t
    call get_LR_mat_sum_pt(j,1,bigr,zmid,vLtmp(j,1),vRtmp(j,1))
enddo
!$acc end loop

!$acc parallel loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop

!$acc parallel loop independent
do j=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(1),l_t(1),m_t(1),vmatL(j, 1),vmid(j,1))
enddo
!$acc end loop

!$acc parallel loop independent
do j=1,num_t
  vmid(j+num_t,1)=direct_num_t(j,1,cosang,bigr)*cmplx(cos(zmid*(k_t(inl_t(1))-k_t(inl_t(j)))),sin(zmid*(k_t(inl_t(1))-k_t(inl_t(j)))),id)
enddo
!$acc end loop

!$acc serial
    bigr=sqrt(rho*rho+znew*znew)
    cosang=znew/bigr
    ang=acos(znew/bigr)
!$acc end serial

!$acc parallel loop independent 
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
    do j=1,nglag
      do k=1,ngleg
        arg=0.5_id*v*znew*u(j)*gleg(k)
        zij1(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)*u(j)-gleg(k)*gleg(k))
        zij2(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)-gleg(k))
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop


!$acc parallel loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc parallel loop independent
do j=1,num_t
    call get_LR_mat_sum_pt(j,1,bigr,znew,vLtmp(j,1),vRtmp(j,1))
enddo
!$acc end loop

!$acc parallel loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop

!$acc parallel loop independent
do j=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(1),l_t(1),m_t(1),vmatL(j, 1),vright(j,1))
enddo
!$acc end loop


!$acc parallel loop independent
do j=1,num_t
  vright(j+num_t,1)=direct_num_t(j,1,cosang,bigr)*cmplx(cos(znew*(k_t(inl_t(1))-k_t(inl_t(j)))),sin(znew*(k_t(inl_t(1))-k_t(inl_t(j)))),id)
enddo
!$acc end loop

!$acc parallel loop
   do j=1,2*num_t 
     dymm(j)=-2._id*cu*vmid(j,1)
     dyt(j)=-cu*vright(j,1)
     zout(j)=zam(j)+h6*(zdadz(j)+dyt(j)+2._id*dymm(j))
   enddo
!$acc end loop

!$acc parallel loop
   do j=1,2*num_t 
     zam(j)=zout(j)
   enddo
!$acc end loop

enddo

!$acc serial
zam(num_t+i0)=zam(num_t+i0)-1._id
!$acc end serial
!$acc parallel loop
    do ist_f=1,2*num_t
       prob(ist_f)=zam(ist_f)*conjg(zam(ist_f))*rho
    enddo
!$acc end loop

!$acc end data

!$acc update host(zam,vover)

zamo(:)=zam(:)

zam(num_t+i0)=zam(num_t+i0)+1._id

    unitarity2=0.0_id
!!$omp parallel default(shared) reduction(+:unitarity2) 
!!$omp do schedule(guided)collapse(2)
    do k=1,2*num_t 
      do j=1,2*num_t 
        unitarity2=unitarity2+conjg(zam(k))*zam(j)*vover(k,j)
      enddo
    enddo
!!$omp end do
!!$omp end parallel

write(1298,*)rho,dble(unitarity2)

      call system_clock(t2)
      dt = t2-t1
      secs = real(dt)/real(count_rate)

      open(1299,file='time',status='unknown')
      write(1299,"('number of channels is ',i5,' wall clock time is ',f12.2,' seconds')") num_t,secs
      open(1300,file='timing',status='unknown')
      write(1300,"(i5,f12.2)") num_t,secs

 30   format(i20,10000g14.6)

deallocate(zam)
deallocate(zout)
deallocate(zdadz)
deallocate(dyt)
deallocate(dym)
deallocate(dymm)
deallocate(zij1)
deallocate(zij2)

end subroutine RK_2c_FBA_gpu


subroutine RK_2c_gpu( rho,prob,zamo )
   use openacc
   use cublas_v2
   use cusolverDn
   use cudafor
   use precision_type
   use data_base,only:num_t,vleft,vmid,vright,vmatl,nstep,step,key,i0,unitarity2,cu,nglag,ngleg,u,glag,gleg,wig_d,lmax_t,mmax_t,wigarr,larr,marr1,marr2,&
&   nlnum_t,inl_t,inl_p,nst_t,lst_t,basis_a,basis_b,basis_bres,jdir,idir,fdir,istar,k_t,k_p,l_p,lob,upb,m_t,zero,zam,zout,zdadz,dyt,dym,dymm,plma,plmb,&
&   lar,mar,lmar,v,jint,n_t,l_t,m_t,vLtmp,vRtmp,vover,z,znew,h,zmid,bigr,cosang,ang,hh,h6,zij1,zij2,wglag,wgleg,un,unn
!   use acc_routines,only:direct_num_t,direct_num_p,get_plmij
   use acc_routines
   use acc_wig
   implicit none
   integer::ist_f,i,j,ist_i,info,k,l,m1,m2,n,nf,lf,mf,ni,li,mi,lm
   real(kind=id),intent(in) :: rho
   real(kind=id),intent(inout),dimension(1:2*num_t) :: prob
   complex(kind=id),intent(out),dimension(1:2*num_t) :: zamo
!   complex(kind=id),dimension(1:num_t,1:num_t) :: vLtmp,vRtmp
!   complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t) :: vmatLtmp,vmatRtmp
   complex(kind=id),dimension(1:2*num_t*2*num_t) ::work
   complex(kind=id)::zrm1,zrm2,zlm1,zlm2
   real(kind=id) :: fct,wigd1,wigd2,arg_bess,arg
!   real(kind=id) :: z,znew,h,zmid,bigr,cosang,ang,hh,h6
   integer,dimension(1:2*num_t)::ipiv
   integer :: t1, t2, dt, count_rate, count_max
   real(kind=id) ::  secs
   integer::Lwork,status
   complex(kind=id),allocatable,dimension(:) ::work_gpu
   type(cusolverDnHandle) :: handle

!    open(220,file='solution',status='unknown')
!    write(220,*)
!    write(220,*)

allocate(un(-nstep:nstep))
allocate(unn(1:2*num_t))
allocate(zam(1:2*num_t))
allocate(zout(1:2*num_t))
allocate(zdadz(1:2*num_t))
allocate(dyt(1:2*num_t))
allocate(dym(1:2*num_t))
allocate(dymm(1:2*num_t))
allocate(zij1(1:nglag,1:ngleg))
allocate(zij2(1:nglag,1:ngleg))

  status = cusolverDnCreate(handle)
  if (status /= CUSOLVER_STATUS_SUCCESS) &
       write(*,*) 'cusolverDnCreate error: ', status

!$acc enter data copyin(work)
  !$acc host_data use_device(work)
  status=cusolverDnZgetrf_bufferSize(handle,2*num_t,2*num_t,work,2*num_t,Lwork)
    if (status /= CUSOLVER_STATUS_SUCCESS) &
         write(*,*) 'cusolverDnZgetrf_buffersize failed ', status
  !$acc end host_data

allocate(work_gpu(1:Lwork))

!$acc enter data copyin(work_gpu,ipiv,info)

   call system_clock(count_max=count_max, count_rate=count_rate)
      call system_clock(t1)


!$acc data present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t,lmar,lar,mar,i0) copyout(prob) copyin(rho,v)
!$acc parallel loop independent 
do ist_f=1,2*num_t
   zam(ist_f)=cmplx(0,0,id)
   zdadz(ist_f)=cmplx(0,0,id)
   zout(ist_f)=cmplx(0,0,id)
enddo
!$acc end loop

!$acc serial
   zam(num_t+i0)=cmplx(1,0,id)
   z=step(-nstep)
   bigr=sqrt(rho*rho+z*z)
   cosang=z/bigr
   ang=acos(z/bigr)
!$acc end serial

!$acc parallel loop independent 
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
    do j=1,nglag
      do k=1,ngleg
        arg=0.5_id*v*z*u(j)*gleg(k)
        zij1(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)*u(j)-gleg(k)*gleg(k))
        zij2(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)-gleg(k))
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop


!$acc parallel loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop


!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_sum_pt(j,k,bigr,z,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc parallel loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j, k+num_t),vright(j,k+num_t))
  enddo
enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(1,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_sum_tp(j,k,bigr,z,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_tp(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j+num_t, k),vright(j+num_t,k))
  enddo
enddo
!$acc end loop

!$acc parallel loop independent
do j=1,jdir
  if(mod(l_p(fdir(j))+l_p(idir(j)),2)==0)then
    vright(fdir(j),idir(j))=direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(z*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(z*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
  else
    vright(fdir(j),idir(j))=-direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(z*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(z*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
  endif
  vright(idir(j),fdir(j))=conjg(vright(fdir(j),idir(j)))
  vmatL(fdir(j),idir(j))=cmplx(0,0,id)
  vmatL(idir(j),fdir(j))=cmplx(0,0,id)
enddo
!$acc end loop
!$acc parallel loop independent
do j=1,num_t
  vright(j,j)=direct_num_p(j,j,cosang,bigr)*cmplx(cos(z*(k_p(inl_p(j))-k_p(inl_p(j)))),sin(z*(k_p(inl_p(j))-k_p(inl_p(j)))),id)
  vright(j+num_t,j+num_t)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
  vmatL(j,j)=cmplx(1,0,id)
  vmatL(j+num_t,j+num_t)=cmplx(1,0,id)
enddo
!$acc end loop
!$acc parallel loop independent
do j=1,jdir
  vright(fdir(j)+num_t,idir(j)+num_t)=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(z*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(z*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
  vright(idir(j)+num_t,fdir(j)+num_t)=conjg(vright(fdir(j)+num_t,idir(j)+num_t))
  vmatL(idir(j)+num_t,fdir(j)+num_t)=cmplx(0,0,id)
  vmatL(fdir(j)+num_t,idir(j)+num_t)=cmplx(0,0,id)
enddo
!$acc end loop

!$acc host_data use_device(vmatL,vright,work_gpu,info,ipiv)
status=cusolverDnZgetrf(handle,2*num_t,2*num_t,vmatL,2*num_t,work_gpu,ipiv,INFO)
status=cusolverDnZgetrs(handle,cublas_op_n,2*num_t,2*num_t,vmatL,2*num_t,ipiv,vright,2*num_t,INFO)
!$acc end host_data

!!$acc update host(vmatL,vright)
!call zhesv('U',2*num_t,2*num_t,vmatL(1:2*num_t,1:2*num_t),2*num_t,ipiv(1:2*num_t),vright(1:2*num_t,1:2*num_t),2*num_t,work,2*num_t*2*num_t,info)
!call Zgetrf(2*num_t,2*num_t,vmatL,2*num_t,ipiv,INFO)
!call Zgetrs('N',2*num_t,2*num_t,vmatL,2*num_t,ipiv,vright,2*num_t,INFO)
!!$acc update device(vright)


do i=-nstep,nstep-1

!$acc serial

    z=step(i)
    znew=step(i+1)
    zmid=(z+znew)/real(2,id)
    h=znew-z
    hh=h/2._id
    h6=h/6._id

!$acc end serial

!$acc parallel loop
   do k=1,2*num_t
     vleft(k,1:2*num_t)=vright(k,1:2*num_t)
     zdadz(k)=-cu*sum(vleft(k,1:2*num_t)*zam(1:2*num_t))
   enddo
!$acc end loop

!$acc serial
    bigr=sqrt(rho*rho+zmid*zmid)
    cosang=zmid/bigr
    ang=acos(zmid/bigr)
!$acc end serial

!$acc parallel loop independent 
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
    do j=1,nglag
      do k=1,ngleg
        arg=0.5_id*v*zmid*u(j)*gleg(k)
        zij1(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)*u(j)-gleg(k)*gleg(k))
        zij2(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)-gleg(k))
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop


!$acc parallel loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_sum_pt(j,k,bigr,zmid,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc parallel loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j, k+num_t),vmid(j,k+num_t))
  enddo
enddo
!$acc end loop



!$acc parallel loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(1,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_sum_tp(j,k,bigr,zmid,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_tp(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j+num_t, k),vmid(j+num_t,k))
  enddo
enddo
!$acc end loop

!$acc parallel loop independent
do j=1,jdir
  if(mod(l_p(fdir(j))+l_p(idir(j)),2)==0)then
    vmid(fdir(j),idir(j))=direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(zmid*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(zmid*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
  else
    vmid(fdir(j),idir(j))=-direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(zmid*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(zmid*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
  endif
  vmid(idir(j),fdir(j))=conjg(vmid(fdir(j),idir(j)))
  vmatL(fdir(j),idir(j))=cmplx(0,0,id)
  vmatL(idir(j),fdir(j))=cmplx(0,0,id)
enddo
!$acc end loop
!$acc parallel loop independent
do j=1,num_t
  vmid(j,j)=direct_num_p(j,j,cosang,bigr)*cmplx(cos(zmid*(k_p(inl_p(j))-k_p(inl_p(j)))),sin(zmid*(k_p(inl_p(j))-k_p(inl_p(j)))),id)
  vmid(j+num_t,j+num_t)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(zmid*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(zmid*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
  vmatL(j,j)=cmplx(1,0,id)
  vmatL(j+num_t,j+num_t)=cmplx(1,0,id)
enddo
!$acc end loop
!$acc parallel loop independent
do j=1,jdir
  vmid(fdir(j)+num_t,idir(j)+num_t)=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(zmid*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(zmid*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
  vmid(idir(j)+num_t,fdir(j)+num_t)=conjg(vmid(fdir(j)+num_t,idir(j)+num_t))
  vmatL(idir(j)+num_t,fdir(j)+num_t)=cmplx(0,0,id)
  vmatL(fdir(j)+num_t,idir(j)+num_t)=cmplx(0,0,id)
enddo
!$acc end loop

!!$acc update host(vmatL,vmid)
!call zhesv('U',2*num_t,2*num_t,vmatL(1:2*num_t,1:2*num_t),2*num_t,ipiv(1:2*num_t),vmid(1:2*num_t,1:2*num_t),2*num_t,work,2*num_t*2*num_t,info)
!call Zgetrf(2*num_t,2*num_t,vmatL,2*num_t,ipiv,INFO)
!call Zgetrs('N',2*num_t,2*num_t,vmatL,2*num_t,ipiv,vmid,2*num_t,INFO)
!!$acc update device(vmid)

!$acc host_data use_device(vmatL,vmid,work_gpu,info,ipiv)
status=cusolverDnZgetrf(handle,2*num_t,2*num_t,vmatL,2*num_t,work_gpu,ipiv,INFO)
status=cusolverDnZgetrs(handle,cublas_op_n,2*num_t,2*num_t,vmatL,2*num_t,ipiv,vmid,2*num_t,INFO)
!$acc end host_data


!$acc serial
    bigr=sqrt(rho*rho+znew*znew)
    cosang=znew/bigr
    ang=acos(znew/bigr)
!$acc end serial

!$acc parallel loop independent 
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
    do j=1,nglag
      do k=1,ngleg
        arg=0.5_id*v*znew*u(j)*gleg(k)
        zij1(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)*u(j)-gleg(k)*gleg(k))
        zij2(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)-gleg(k))
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc parallel loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop


!$acc parallel loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_sum_pt(j,k,bigr,znew,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc parallel loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j, k+num_t),vright(j,k+num_t))
  enddo
enddo
!$acc end loop



!$acc parallel loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(1,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_sum_tp(j,k,bigr,znew,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc parallel loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_tp(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j+num_t, k),vright(j+num_t,k))
  enddo
enddo
!$acc end loop

!$acc parallel loop independent
do j=1,jdir
  if(mod(l_p(fdir(j))+l_p(idir(j)),2)==0)then
    vright(fdir(j),idir(j))=direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(znew*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(znew*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
  else
    vright(fdir(j),idir(j))=-direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(znew*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(znew*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
  endif
  vright(idir(j),fdir(j))=conjg(vright(fdir(j),idir(j)))
  vmatL(fdir(j),idir(j))=cmplx(0,0,id)
  vmatL(idir(j),fdir(j))=cmplx(0,0,id)
enddo
!$acc end loop
!$acc parallel loop independent
do j=1,num_t
  vright(j,j)=direct_num_p(j,j,cosang,bigr)*cmplx(cos(znew*(k_p(inl_p(j))-k_p(inl_p(j)))),sin(znew*(k_p(inl_p(j))-k_p(inl_p(j)))),id)
  vright(j+num_t,j+num_t)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(znew*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(znew*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
  vmatL(j,j)=cmplx(1,0,id)
  vmatL(j+num_t,j+num_t)=cmplx(1,0,id)
enddo
!$acc end loop
!$acc parallel loop independent
do j=1,jdir
  vright(fdir(j)+num_t,idir(j)+num_t)=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(znew*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(znew*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
  vright(idir(j)+num_t,fdir(j)+num_t)=conjg(vright(fdir(j)+num_t,idir(j)+num_t))
  vmatL(idir(j)+num_t,fdir(j)+num_t)=cmplx(0,0,id)
  vmatL(fdir(j)+num_t,idir(j)+num_t)=cmplx(0,0,id)
enddo
!$acc end loop
!$acc parallel loop independent collapse(2)
do j=1,2*num_t
  do k=1,2*num_t
    vover(k,j)=vmatL(k,j)
  enddo
enddo
!$acc end loop

!!$acc update host(vmatL,vright)
!call zhesv('U',2*num_t,2*num_t,vmatL(1:2*num_t,1:2*num_t),2*num_t,ipiv(1:2*num_t),vright(1:2*num_t,1:2*num_t),2*num_t,work,2*num_t*2*num_t,info)
!call Zgetrf(2*num_t,2*num_t,vmatL,2*num_t,ipiv,INFO)
!call Zgetrs('N',2*num_t,2*num_t,vmatL,2*num_t,ipiv,vright,2*num_t,INFO)
!!$acc update device(vright)

!$acc host_data use_device(vmatL,vright,work_gpu,info,ipiv)
status=cusolverDnZgetrf(handle,2*num_t,2*num_t,vmatL,2*num_t,work_gpu,ipiv,INFO)
status=cusolverDnZgetrs(handle,cublas_op_n,2*num_t,2*num_t,vmatL,2*num_t,ipiv,vright,2*num_t,INFO)
!$acc end host_data

!$acc parallel loop
   do j=1,2*num_t 
     dyt(j)=-cu*sum(vmid(j,1:2*num_t)*(zam(1:2*num_t)+hh*zdadz(1:2*num_t)))
   enddo
!$acc end loop

!$acc parallel loop
   do j=1,2*num_t 
     dym(j)=-cu*sum(vmid(j,1:2*num_t)*(zam(1:2*num_t)+hh*dyt(1:2*num_t)))
   enddo
!$acc end loop

!$acc parallel loop
   do j=1,2*num_t 
     dymm(j)=dyt(j)+dym(j)
     dyt(j)=-cu*sum(vright(j,1:2*num_t)*(zam(1:2*num_t)+h*dym(1:2*num_t)))
     zout(j)=zam(j)+h6*(zdadz(j)+dyt(j)+2._id*dymm(j))
   enddo
!$acc end loop

!$acc parallel loop
   do j=1,2*num_t 
     zam(j)=zout(j)
   enddo
!$acc end loop

!$acc serial !parallel loop 
    do k=1,2*num_t 
      unn(k)=sum(zam(1:2*num_t)*vover(k,1:2*num_t))
    enddo
!$acc end serial  !loop

!$acc serial
!   un(i)=sum(conjg(zam(1:2*num_t))*unn(1:2*num_t)) 
   un(i)=sum(conjg(zam(1:2*num_t))*unn(1:2*num_t)) 
!$acc end serial

enddo

!$acc serial
zam(num_t+i0)=zam(num_t+i0)-1._id
!$acc end serial
!$acc parallel loop
    do ist_f=1,2*num_t
       prob(ist_f)=zam(ist_f)*conjg(zam(ist_f))*rho
    enddo
!$acc end loop

!$acc end data

!$acc update host(zam,vover,un,vright)

zamo(:)=zam(:)

zam(num_t+i0)=zam(num_t+i0)+1._id

    unitarity2=0.0_id
!!$omp parallel default(shared) reduction(+:unitarity2) 
!!$omp do schedule(guided)collapse(2)
    do k=1,2*num_t 
      do j=1,2*num_t 
        unitarity2=unitarity2+conjg(zam(k))*zam(j)*vover(k,j)
      enddo
    enddo
!!$omp end do
!!$omp end parallel

do i=-nstep,nstep-1
  write(699,*)i,un(i)
enddo

print*, 'gpu', z,sum(vright(:,:))

write(1298,*)rho,dble(unitarity2)

      call system_clock(t2)
      dt = t2-t1
      secs = real(dt)/real(count_rate)

      open(1299,file='time',status='unknown')
      write(1299,"('number of channels is ',i5,' wall clock time is ',f12.2,' seconds')") num_t,secs
      open(1300,file='timing',status='unknown')
      write(1300,"(i5,f12.2)") num_t,secs

 30   format(i20,10000g14.6)

!$acc exit data delete(work,work_gpu,ipiv,info)
deallocate(zam)
deallocate(un)
deallocate(unn)
deallocate(zout)
deallocate(zdadz)
deallocate(dyt)
deallocate(dym)
deallocate(dymm)
deallocate(zij1)
deallocate(zij2)

end subroutine RK_2c_gpu

subroutine solveRK_num_1c_prj( rho,prob,zamo )
   use openacc
   use precision_type
!   use data_base,only:num_t,vleft,vmid,vright,nstep,step,i0,unitarity2,cu,nglag,ngleg,u,glag,wig_d,lmax_t,mmax_t,fdir,idir,jdir,k_t,inl_t,zam,zout,zdadz,dyt,dym,dymm
   use data_base,only:num_t,vleft,vmid,vright,vmatl,nstep,step,key,i0,unitarity2,cu,nglag,ngleg,u,glag,gleg,wig_d,lmax_t,mmax_t,wigarr,larr,marr1,marr2,&
&   nlnum_t,inl_t,inl_p,nst_t,lst_t,basis_a,basis_b,basis_bres,jdir,idir,fdir,istar,k_t,k_p,l_p,lob,upb,m_t,zero,zam,zout,zdadz,dyt,dym,dymm,plma,plmb,&
&   lar,mar,lmar,v,jint,n_t,l_t,m_t,vLtmp,vRtmp,vover,z,znew,h,zmid,bigr,cosang,ang,hh,h6,zij1,zij2,wglag,wgleg
!   use acc_routines,only:direct_num_t
   use acc_routines
   use acc_wig
   implicit none
!   integer::ist_f,i,j,ist_i,k
   integer:: fin,iin,finl,iinl
   integer::ist_f,i,j,ist_i,info,k,l,m1,m2,n,nf,lf,mf,ni,li,mi,lm
   real(kind=id),intent(in) :: rho
   real(kind=id),intent(inout),dimension(1:2*num_t) :: prob
!   complex(kind=id),dimension(1:num_t) :: zam,zout,zdadz,dyt,dym,dymm
!   complex(kind=id) :: tmp 
   complex(kind=id),intent(out),dimension(1:2*num_t) :: zamo
!   real(kind=id) :: z,znew,h,zmid,bigr,cosang,hh,h6
   real(kind=id) :: kf,ki
   integer :: t1, t2, dt, count_rate, count_max
   real(kind=id) ::  secs
   complex(kind=id)::zrm1,zrm2,zlm1,zlm2
   real(kind=id) :: fct,wigd1,wigd2,arg_bess,arg

!    open(220,file='solution',status='unknown')
!    write(220,*)
!    write(220,*)

allocate(zij1(1:nglag,1:ngleg))
allocate(zij2(1:nglag,1:ngleg))
allocate(zam(1:2*num_t))
allocate(zout(1:2*num_t))
allocate(zdadz(1:num_t))
allocate(dyt(1:num_t))
allocate(dym(1:num_t))
allocate(dymm(1:num_t))

print*,'I am here'
   call system_clock(count_max=count_max, count_rate=count_rate)
      call system_clock(t1)

!!$acc data present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t) copyout(cosang,bigr,z,znew,zmid,prob) copyin(rho)  !copyout(prob)
!!$acc data present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t,lmar,lar,mar,i0) copyout(prob,cosang,bigr,z,zmid,znew) copyin(rho,v)
!$acc data present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t,lmar,lar,mar,i0) copyout(cosang,bigr,z,znew,zmid,prob) copyin(rho,v)
!$acc kernels   !present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t) copyin(rho) copyout(prob,cosang,bigr,z) 
!$acc loop independent 
do ist_f=1,num_t
   zam(ist_f)=cmplx(0,0,id)
   zam(ist_f+num_t)=cmplx(0,0,id)
   zdadz(ist_f)=cmplx(0,0,id)
   zout(ist_f)=cmplx(0,0,id)
enddo
!$acc end loop 
!!$acc serial
zam(i0)=cmplx(1,0,id)
z=step(-nstep)
bigr=sqrt(rho*rho+z*z)
cosang=z/bigr
ang=acos(z/bigr)
!!$acc end serial
!$acc loop independent
do j=1,jdir
  vright(fdir(j),idir(j))=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(z*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(z*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
  vright(idir(j),fdir(j))=conjg(vright(fdir(j),idir(j)))
enddo
!$acc end loop 


!$acc loop
do j=1,num_t
  vright(j,j)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
enddo
!$acc end loop 


!$acc loop independent 
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc loop independent collapse(2)
    do j=1,nglag
      do k=1,ngleg
        arg=0.5_id*v*z*u(j)*gleg(k)
        zij1(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)*u(j)-gleg(k)*gleg(k))
        zij2(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)-gleg(k))
      enddo
    enddo
!$acc end loop

!$acc loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop


!$acc loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(1,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_sum_tp(j,k,bigr,z,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_tp(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j, k+num_t),vright(j,k+num_t))
  enddo
enddo
!$acc end loop


!$acc end kernels



do i=-nstep,nstep-1

!$acc kernels

   z=step(i)



   znew=step(i+1)
   zmid=(z+znew)/real(2,id)
   h=znew-z
   hh=h/2._id
   h6=h/6._id


!$acc loop 
   do k=1,num_t 
     vleft(k,1:2*num_t)=vright(k,1:2*num_t)
!     vleft(k,num_t+1:2*num_t)=vmatL(k,num_t+1:2*num_t)
     zdadz(k)=-cu*sum(vleft(k,1:num_t)*zam(1:num_t))
   enddo
!$acc end loop 


   bigr=sqrt(rho*rho+zmid*zmid)
   cosang=zmid/bigr

!$acc loop independent
   do j=1,jdir
     vmid(fdir(j),idir(j))=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(zmid*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(zmid*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
     vmid(idir(j),fdir(j))=conjg(vmid(fdir(j),idir(j)))
   enddo
!$acc end loop 
!$acc loop 
   do j=1,num_t
     vmid(j,j)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
   enddo
!$acc end loop 



   bigr=sqrt(rho*rho+znew*znew)
   cosang=znew/bigr
   ang=acos(znew/bigr)

!$acc loop independent 
   do j=1,jdir
     vright(fdir(j),idir(j))=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(znew*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(znew*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
     vright(idir(j),fdir(j))=conjg(vright(fdir(j),idir(j)))
   enddo
!$acc end loop 
!$acc loop 
   do j=1,num_t
     vright(j,j)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
   enddo
!$acc end loop 








!$acc loop independent 
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc loop independent collapse(2)
    do j=1,nglag
      do k=1,ngleg
        arg=0.5_id*v*znew*u(j)*gleg(k)
        zij1(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)*u(j)-gleg(k)*gleg(k))
        zij2(j,k)=wglag(j)*wgleg(k)*cmplx(cos(arg),sin(arg),id)*(u(j)-gleg(k))
      enddo
    enddo
!$acc end loop

!$acc loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop


!$acc loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(1,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_sum_tp(j,k,bigr,znew,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_tp(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j, k+num_t),vright(j,k+num_t))
  enddo
enddo
!$acc end loop










!$acc loop 
   do j=1,num_t 
     dyt(j)=-cu*sum(vmid(j,1:num_t)*(zam(1:num_t)+hh*zdadz(1:num_t)))
   enddo
!$acc end loop
!$acc loop 
   do j=1,num_t 
     dym(j)=-cu*sum(vmid(j,1:num_t)*(zam(1:num_t)+hh*dyt(1:num_t)))
   enddo
!$acc end loop 
!$acc loop 
   do j=1,num_t 
     dymm(j)=dyt(j)+dym(j)
     dyt(j)=-cu*sum(vright(j,1:num_t)*(zam(1:num_t)+h*dym(1:num_t)))
     zout(j)=zam(j)+h6*(zdadz(j)+dyt(j)+2._id*dymm(j))
   enddo
!$acc end loop 

!$acc loop 
   do j=1,num_t 
     zam(j+num_t)=zam(j+num_t)+(sum(vleft(j,num_t+1:2*num_t)*zam(1:num_t))+sum(vright(j,num_t+1:2*num_t)*zout(1:num_t)))/2._id*h
!     zam(j+num_t)=zam(j+num_t)+(sum(vleft(j,1:num_t)*zam(1:num_t))+sum(vright(j,1:num_t)*zout(1:num_t)))/2._id*h
!     zam(j+num_t)=zam(j+num_t)+(sum(vleft(j,num_t+1:2*num_t)*zam(1:num_t))+sum(vmatL(j,num_t+1:2*num_t)*zout(1:num_t)))/2._id*h
   enddo
!$acc end loop 

!$acc loop 
   do j=1,num_t 
     zam(j)=zout(j)
   enddo
!$acc end loop 



!$acc end kernels 
!!$acc update host(zam,z,dymm,bigr)
!print*,z,bigr,zam(1:3),dymm(1:3)
!stop

enddo




!$acc kernels 

!
!!$omp do schedule(guided)
!    do j=1,num_t 
!      unitarity2=unitarity2+abs(zam(j))**2
!    enddo
!!$omp end do

unitarity2=sum(abs(zam(1:num_t))**2)

    zam(i0)=zam(i0)-1._id


!$acc loop 
    do ist_f=1,2*num_t
       prob(ist_f)=zam(ist_f)*conjg(zam(ist_f))*rho
    enddo
!$acc end loop 
!!$acc end parallel

!$acc end kernels 



!$acc end data
!$acc update host(zam,unitarity2)


zamo(:)=zam(:)

    write(1298,*)rho,unitarity2

      call system_clock(t2)
      dt = t2-t1
      secs = real(dt)/real(count_rate)

      open(1299,file='time',status='unknown')
      write(1299,"('number of channels is ',i5,' wall clock time is ',f12.2,' seconds')") num_t,secs
      open(1300,file='timing',status='unknown')
      write(1300,"(i5,f12.2)") num_t,secs

 30   format(i20,10000g14.6)

deallocate(zam)
deallocate(zout)
deallocate(zdadz)
deallocate(dyt)
deallocate(dym)
deallocate(dymm)
deallocate(zij1)
deallocate(zij2)

end subroutine solveRK_num_1c_prj




subroutine solveRK_num_2c_dev( rho,prob,zamo )
   use precision_type
   use data_base,only:num_t,vleft,vmid,vright,vmatl,nstep,step,key,i0,unitarity2,cu,nglag,ngleg,u,glag,gleg,wig_d,lmax_t,mmax_t,wigarr,larr,marr1,marr2,&
&   nlnum_t,inl_t,inl_p,nst_t,lst_t,basis_a,basis_b,basis_bres,jdir,idir,fdir,istar,k_t,k_p,l_p,lob,upb,m_t,zero,zam,zout,zdadz,dyt,dym,dymm,plma,plmb,&
&   lar,mar,lmar,v,jint,n_t,l_t,m_t,vLtmp,vRtmp,vover,z,znew,h,zmid,bigr,cosang,ang,hh,h6
!   use acc_routines,only:direct_num_t,direct_num_p,get_plmij
   use acc_routines
   use acc_wig
   implicit none
   integer::ist_f,i,j,ist_i,info,k,l,m1,m2,n,nf,lf,mf,ni,li,mi,lm
   real(kind=id),intent(in) :: rho
   real(kind=id),intent(inout),dimension(1:2*num_t) :: prob
   complex(kind=id),intent(out),dimension(1:2*num_t) :: zamo
!   complex(kind=id),dimension(1:num_t,1:num_t) :: vLtmp,vRtmp
!   complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t) :: vmatLtmp,vmatRtmp
   complex(kind=id),dimension(1:2*num_t*2*num_t) ::work
   complex(kind=id)::zrm1,zrm2,zlm1,zlm2
   real(kind=id) :: fct,wigd1,wigd2,arg_bess
!   real(kind=id) :: z,znew,h,zmid,bigr,cosang,ang,hh,h6
   integer,dimension(1:2*num_t)::ipiv
   integer :: t1, t2, dt, count_rate, count_max
   real(kind=id) ::  secs

!    open(220,file='solution',status='unknown')
!    write(220,*)
!    write(220,*)

print*,'I am here in 2c'
allocate(zam(1:2*num_t))
allocate(zout(1:2*num_t))
allocate(zdadz(1:2*num_t))
allocate(dyt(1:2*num_t))
allocate(dym(1:2*num_t))
allocate(dymm(1:2*num_t))

   call system_clock(count_max=count_max, count_rate=count_rate)
      call system_clock(t1)


!$acc data present(step,nstep,num_t,jdir,fdir,idir,k_t,inl_t,lmar,lar,mar,i0) copyout(prob) copyin(rho,v)
!$acc kernels
!$acc loop independent 
do ist_f=1,2*num_t
   zam(ist_f)=cmplx(0,0,id)
   zdadz(ist_f)=cmplx(0,0,id)
   zout(ist_f)=cmplx(0,0,id)
enddo
!$acc end loop

   zam(num_t+i0)=cmplx(1,0,id)
   z=step(-nstep)
   bigr=sqrt(rho*rho+z*z)
   cosang=z/bigr
   ang=acos(z/bigr)

!$acc loop independent 
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop

!$acc loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop

!$acc loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_dev_pt2(j,k,bigr,z,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j, k+num_t),vright(j,k+num_t))
  enddo
enddo
!$acc end loop



!$acc loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(1,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_dev_tp2(j,k,bigr,z,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t




    call rotation_tp(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j+num_t, k),vright(j+num_t,k))
        

  enddo
enddo
!$acc end loop

!$acc loop independent
do j=1,jdir
  if(mod(l_p(fdir(j))+l_p(idir(j)),2)==0)then
    vright(fdir(j),idir(j))=direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(z*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(z*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
  else
    vright(fdir(j),idir(j))=-direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(z*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(z*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
  endif
  vright(idir(j),fdir(j))=conjg(vright(fdir(j),idir(j)))
  vmatL(fdir(j),idir(j))=cmplx(0,0,id)
  vmatL(idir(j),fdir(j))=cmplx(0,0,id)
enddo
!$acc end loop
!$acc loop independent
do j=1,num_t
  vright(j,j)=direct_num_p(j,j,cosang,bigr)*cmplx(cos(z*(k_p(inl_p(j))-k_p(inl_p(j)))),sin(z*(k_p(inl_p(j))-k_p(inl_p(j)))),id)
  vright(j+num_t,j+num_t)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
  vmatL(j,j)=cmplx(1,0,id)
  vmatL(j+num_t,j+num_t)=cmplx(1,0,id)
enddo
!$acc end loop
!$acc loop independent
do j=1,jdir
  vright(fdir(j)+num_t,idir(j)+num_t)=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(z*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(z*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
  vright(idir(j)+num_t,fdir(j)+num_t)=conjg(vright(fdir(j)+num_t,idir(j)+num_t))
  vmatL(idir(j)+num_t,fdir(j)+num_t)=cmplx(0,0,id)
  vmatL(fdir(j)+num_t,idir(j)+num_t)=cmplx(0,0,id)
enddo
!$acc end loop

!$acc end kernels

!!$acc update host(vmatL,vright)
!
!call zhesv('U',2*num_t,2*num_t,vmatL(1:2*num_t,1:2*num_t),2*num_t,ipiv(1:2*num_t),vright(1:2*num_t,1:2*num_t),2*num_t,work,2*num_t*2*num_t,info)
!
!!$acc update device(vright)


do i=-nstep,nstep-1

!$acc kernels

    z=step(i)
    znew=step(i+1)
    zmid=(z+znew)/real(2,id)
    h=znew-z
    hh=h/2._id
    h6=h/6._id

!$acc loop
   do k=1,2*num_t
     vleft(k,1:2*num_t)=vright(k,1:2*num_t)
     zdadz(k)=-cu*sum(vleft(k,1:2*num_t)*zam(1:2*num_t))
   enddo
!$acc end loop

    bigr=sqrt(rho*rho+zmid*zmid)
    cosang=zmid/bigr
    ang=acos(zmid/bigr)

!$acc loop independent
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop

!$acc loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop

!$acc loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_dev_pt2(j,k,bigr,zmid,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop


!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j, k+num_t),vmid(j,k+num_t))
  enddo
enddo
!$acc end loop



!$acc loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(1,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_dev_tp2(j,k,bigr,zmid,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t




    call rotation_tp(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j+num_t, k),vmid(j+num_t,k))
        

  enddo
enddo
!$acc end loop

!$acc loop independent
do j=1,jdir
  if(mod(l_p(fdir(j))+l_p(idir(j)),2)==0)then
    vmid(fdir(j),idir(j))=direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(zmid*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(zmid*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
  else
    vmid(fdir(j),idir(j))=-direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(zmid*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(zmid*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
  endif
  vmid(idir(j),fdir(j))=conjg(vmid(fdir(j),idir(j)))
  vmatL(fdir(j),idir(j))=cmplx(0,0,id)
  vmatL(idir(j),fdir(j))=cmplx(0,0,id)
enddo
!$acc end loop
!$acc loop independent
do j=1,num_t
  vmid(j,j)=direct_num_p(j,j,cosang,bigr)*cmplx(cos(zmid*(k_p(inl_p(j))-k_p(inl_p(j)))),sin(zmid*(k_p(inl_p(j))-k_p(inl_p(j)))),id)
  vmid(j+num_t,j+num_t)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(zmid*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(zmid*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
  vmatL(j,j)=cmplx(1,0,id)
  vmatL(j+num_t,j+num_t)=cmplx(1,0,id)
enddo
!$acc end loop
!$acc loop independent
do j=1,jdir
  vmid(fdir(j)+num_t,idir(j)+num_t)=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(zmid*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(zmid*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
  vmid(idir(j)+num_t,fdir(j)+num_t)=conjg(vmid(fdir(j)+num_t,idir(j)+num_t))
  vmatL(idir(j)+num_t,fdir(j)+num_t)=cmplx(0,0,id)
  vmatL(fdir(j)+num_t,idir(j)+num_t)=cmplx(0,0,id)
enddo
!$acc end loop

!!$acc end kernels
!
!!$acc update host(vmatL,vmid)
!
!
!call zhesv('U',2*num_t,2*num_t,vmatL(1:2*num_t,1:2*num_t),2*num_t,ipiv(1:2*num_t),vmid(1:2*num_t,1:2*num_t),2*num_t,work,2*num_t*2*num_t,info)
!
!!$acc update device(vmid)
!
!!$acc kernels

    bigr=sqrt(rho*rho+znew*znew)
    cosang=znew/bigr
    ang=acos(znew/bigr)

!$acc loop independent
    do j=1,nglag
       u(j)=glag(j)/bigr+1._id
    enddo
!$acc end loop

!$acc loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do lm=1,lmar
          call get_plmijlm(j,k,lm,rho,plma(lar(lm),mar(lm),j,k),plmb(lar(lm),mar(lm),j,k))
        enddo
      enddo
    enddo
!$acc end loop

!$acc loop independent collapse(3)
    do j=1,nglag
      do k=1,ngleg
        do mi=0,2*mmax_t
          arg_bess=0.5_id*v*rho*sqrt((u(j)*u(j)-1._id)*(1._id-gleg(k)*gleg(k)))
          jint(mi,j,k)=bessel_jn(mi,arg_bess)
          if(mod(mi,2)==0)then
            jint(-mi,j,k)=jint(mi,j,k)
          else
            jint(-mi,j,k)=-jint(mi,j,k)
          endif
        enddo
      enddo
    enddo
!$acc end loop

!$acc loop independent
    do l=1,wigarr
      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
    enddo
!$acc end loop

!$acc loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_dev_pt2(j,k,bigr,znew,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call rotation_pt(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j, k+num_t),vright(j,k+num_t))
  enddo
enddo
!$acc end loop



!$acc loop independent collapse(3)
  do n=1,nlnum_t
    do j=1,nglag
      do k=1,ngleg
        call wf_at_r(1,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
      enddo
    enddo
  enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t
    call get_LR_mat_dev_tp2(j,k,bigr,znew,vLtmp(j,k),vRtmp(j,k))
  enddo
enddo
!$acc end loop

!$acc loop independent collapse(2)
do j=1,num_t
  do k=1,num_t




    call rotation_tp(n_t(j),l_t(j),m_t(j),n_t(k),l_t(k),m_t(k),vmatL(j+num_t, k),vright(j+num_t,k))
        

  enddo
enddo
!$acc end loop

!$acc loop independent
do j=1,jdir
  if(mod(l_p(fdir(j))+l_p(idir(j)),2)==0)then
    vright(fdir(j),idir(j))=direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(znew*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(znew*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
  else
    vright(fdir(j),idir(j))=-direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(znew*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(znew*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
  endif
  vright(idir(j),fdir(j))=conjg(vright(fdir(j),idir(j)))
  vmatL(fdir(j),idir(j))=cmplx(0,0,id)
  vmatL(idir(j),fdir(j))=cmplx(0,0,id)
enddo
!$acc end loop
!$acc loop independent
do j=1,num_t
  vright(j,j)=direct_num_p(j,j,cosang,bigr)*cmplx(cos(znew*(k_p(inl_p(j))-k_p(inl_p(j)))),sin(znew*(k_p(inl_p(j))-k_p(inl_p(j)))),id)
  vright(j+num_t,j+num_t)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(znew*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(znew*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
  vmatL(j,j)=cmplx(1,0,id)
  vmatL(j+num_t,j+num_t)=cmplx(1,0,id)
enddo
!$acc end loop
!$acc loop independent
do j=1,jdir
  vright(fdir(j)+num_t,idir(j)+num_t)=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(znew*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(znew*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
  vright(idir(j)+num_t,fdir(j)+num_t)=conjg(vright(fdir(j)+num_t,idir(j)+num_t))
  vmatL(idir(j)+num_t,fdir(j)+num_t)=cmplx(0,0,id)
  vmatL(fdir(j)+num_t,idir(j)+num_t)=cmplx(0,0,id)
enddo
!$acc end loop
!$acc loop independent collapse(2)
    do k=1,2*num_t 
      do j=1,2*num_t 
        vover(k,j)=vmatL(k,j)
      enddo
    enddo
!$acc end loop

!!$acc end kernels
!
!!$acc update host(vmatL,vright)
!
!
!call zhesv('U',2*num_t,2*num_t,vmatL(1:2*num_t,1:2*num_t),2*num_t,ipiv(1:2*num_t),vright(1:2*num_t,1:2*num_t),2*num_t,work,2*num_t*2*num_t,info)
!
!
!
!!$acc update device(vright)
!
!!$acc kernels
!$acc loop
   do j=1,2*num_t 
     dyt(j)=-cu*sum(vmid(j,1:2*num_t)*(zam(1:2*num_t)+hh*zdadz(1:2*num_t)))
   enddo
!$acc end loop

!$acc loop
   do j=1,2*num_t 
     dym(j)=-cu*sum(vmid(j,1:2*num_t)*(zam(1:2*num_t)+hh*dyt(1:2*num_t)))
   enddo
!$acc end loop

!$acc loop
   do j=1,2*num_t 
     dymm(j)=dyt(j)+dym(j)
     dyt(j)=-cu*sum(vright(j,1:2*num_t)*(zam(1:2*num_t)+h*dym(1:2*num_t)))
     zout(j)=zam(j)+h6*(zdadz(j)+dyt(j)+2._id*dymm(j))
   enddo
!$acc end loop

!$acc loop
   do j=1,2*num_t 
     zam(j)=zout(j)
   enddo
!$acc end loop
!$acc end kernels

enddo

!$acc kernels
zam(num_t+i0)=zam(num_t+i0)-1._id
!$acc loop
    do ist_f=1,2*num_t
       prob(ist_f)=zam(ist_f)*conjg(zam(ist_f))*rho
    enddo
!$acc end loop
!$acc end kernels

!$acc end data

!$acc update host(zam,vover)

zamo(:)=zam(:)

zam(num_t+i0)=zam(num_t+i0)+1._id

    unitarity2=0.0_id
!!$omp parallel default(shared) reduction(+:unitarity2) 
!!$omp do schedule(guided)collapse(2)
    do k=1,2*num_t 
      do j=1,2*num_t 
        unitarity2=unitarity2+conjg(zam(k))*zam(j)*vover(k,j)
      enddo
    enddo
!!$omp end do
!!$omp end parallel

write(1298,*)rho,dble(unitarity2)

      call system_clock(t2)
      dt = t2-t1
      secs = real(dt)/real(count_rate)

      open(1299,file='time',status='unknown')
      write(1299,"('number of channels is ',i5,' wall clock time is ',f12.2,' seconds')") num_t,secs
      open(1300,file='timing',status='unknown')
      write(1300,"(i5,f12.2)") num_t,secs

 30   format(i20,10000g14.6)

end subroutine solveRK_num_2c_dev



!subroutine solveRK_num_2c( rho,prob,zam )
!   use precision_type
!   use data_base,only:num_t,vleft,vmid,vright,nstep,step,key,i0,unitarity2,cu,nglag,ngleg,u,glag,gleg,wig_d,lmax_t,mmax_t,wigarr,larr,marr1,marr2,nlnum_t,inl_t,inl_p,nst_t,lst_t,basis_a,basis_b,basis_bres,jdir,idir,fdir,istar,k_t,k_p,l_p,lob,upb,m_t,zero
!   use acc_routines
!   use acc_wig
!   implicit none
!   integer::ist_f,i,j,ist_i,info,k,l,m1,m2,n,nf,lf,mf,ni,li,mi
!   real(kind=id),intent(in) :: rho
!   real(kind=id),intent(inout),dimension(1:2*num_t) :: prob
!   complex(kind=id),dimension(1:2*num_t,1:2*num_t) :: vmatL,vmatR,vover
!   complex(kind=id),dimension(1:num_t,1:num_t) :: vLtmp,vRtmp
!   complex(kind=id),dimension(-lmax_t:lmax_t,-lmax_t:lmax_t) :: vmatLtmp,vmatRtmp
!   complex(kind=id),dimension(1:2*num_t) :: zam,zout,zdadz,dyt,dym,dymm
!   complex(kind=id),dimension(1:2*num_t*2*num_t) ::work
!   real(kind=id) :: z,znew,h,zmid,bigr,cosang,fct,ang,hh,h6,wigd1,wigd2
!   integer,dimension(1:2*num_t)::ipiv
!   integer :: t1, t2, dt, count_rate, count_max
!   real(kind=id) ::  secs
!
!!    open(220,file='solution',status='unknown')
!!    write(220,*)
!!    write(220,*)
!
!   call system_clock(count_max=count_max, count_rate=count_rate)
!      call system_clock(t1)
!
!!$omp parallel default(shared) private(nf,lf,mf,ni,li,mi,wigd1,wigd2,vmatLtmp,vmatRtmp,fct,ipiv,work,info) reduction(+:unitarity2) 
!
!!$omp do schedule(guided)
!do ist_f=1,2*num_t
!   zam(ist_f)=cmplx(0,0,id)
!   zdadz(ist_f)=cmplx(0,0,id)
!   zout(ist_f)=cmplx(0,0,id)
!enddo
!!$omp end do
!
!!$omp single
!   zam(num_t+i0)=cmplx(1,0,id)
!   z=step(-nstep)
!   bigr=sqrt(rho*rho+z*z)
!   cosang=z/bigr
!   ang=acos(z/bigr)
!!$omp end single
!
!!$omp do schedule(guided)
!    do j=1,nglag
!       u(j)=glag(j)/bigr+1._id
!    enddo
!!$omp end do
!
!!$omp do schedule(guided)collapse(2)
!    do j=1,nglag
!      do k=1,ngleg
!        call get_plmij(j,k,rho)
!      enddo
!    enddo
!!$omp end do
!
!!$omp do schedule(guided)
!    do l=1,wigarr
!      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
!    enddo
!!$omp end do
!
!
!!$omp do schedule(guided)collapse(3)
!  do n=1,nlnum_t
!    do j=1,nglag
!      do k=1,ngleg
!        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
!      enddo
!    enddo
!  enddo
!!$omp end do
!
!!$omp do schedule(guided)collapse(2)
!do j=1,num_t
!  do k=1,num_t
!    call get_LR_mat_dev_pt2(j,k,bigr,z,vLtmp(j,k),vRtmp(j,k))
!  enddo
!enddo
!!$omp end do
!
!!$omp do schedule(guided)collapse(2)
!do j=1,nlnum_t
!  do k=1,nlnum_t
!
!
!    nf=nst_t(j)
!    lf=lst_t(j)
!    ni=nst_t(k)
!    li=lst_t(k)
!    do ist_f=lob(j),upb(j)
!      mf=m_t(ist_f)
!      do ist_i=lob(k),upb(k)
!        mi=m_t(ist_i)
!
!        vmatL(ist_f,ist_i+num_t)=zero   
!        vright(ist_f,ist_i+num_t)=zero   
!        do m1=-min(lf,mmax_t),min(lf,mmax_t)
!          wigd1=wig_d(lf,m1,mf)
!          do m2=-min(li,mmax_t),min(li,mmax_t)
!            wigd2=wig_d(li,m2,mi)
!            vmatL(ist_f,ist_i+num_t)=vmatL(ist_f,ist_i+num_t)+wigd1*wigd2*vLtmp(istar(nf,lf,m1),istar(ni,li,m2))
!            vright(ist_f,ist_i+num_t)=vright(ist_f,ist_i+num_t)+wigd1*wigd2*vRtmp(istar(nf,lf,m1),istar(ni,li,m2))
!          enddo
!        enddo
!        vmatL(ist_f, ist_i+num_t)=vmatL(ist_f ,ist_i+num_t)*(-1._id)**(lf+li+mf+mi)
!        vright(ist_f,ist_i+num_t)=vright(ist_f,ist_i+num_t)*(-1._id)**(lf+li+mf+mi)   
!        
!      enddo
!    enddo
!
!  enddo
!enddo
!!$omp end do
!
!!    vmatL(lob(j):upb(j),lob(k):upb(k))=matmul(,matmul(transpose(wig_d(lst_t(k),-min(lst_t(k),mmax_t):min(lst_t(k),mmax_t),-min(lst_t(k),mmax_t):min(lst_t(k),mmax_t))),transpose(vLtmp(lob(j):upb(j),lob(k):upb(k)))))
!
!!$omp do schedule(guided)collapse(3)
!  do n=1,nlnum_t
!    do j=1,nglag
!      do k=1,ngleg
!        call wf_at_r(1,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
!      enddo
!    enddo
!  enddo
!!$omp end do
!
!
!!$omp do schedule(guided)collapse(2)
!do j=1,num_t
!  do k=1,num_t
!    call get_LR_mat_dev_tp2(j,k,bigr,z,vLtmp(j,k),vRtmp(j,k))
!  enddo
!enddo
!!$omp end do
!
!!$omp do schedule(guided)collapse(2)
!do j=1,nlnum_t
!  do k=1,nlnum_t
!
!
!    nf=nst_t(j)
!    lf=lst_t(j)
!    ni=nst_t(k)
!    li=lst_t(k)
!    do ist_f=lob(j),upb(j)
!      mf=m_t(ist_f)
!      do ist_i=lob(k),upb(k)
!        mi=m_t(ist_i)
!        vmatL(ist_f+num_t,ist_i)=zero   
!        vright(ist_f+num_t,ist_i)=zero   
!        do m1=-min(lf,mmax_t),min(lf,mmax_t)
!          wigd1=wig_d(lf,m1,mf)
!          do m2=-min(li,mmax_t),min(li,mmax_t)
!            wigd2=wig_d(li,m2,mi)
!            vmatL(ist_f+num_t,ist_i)=vmatL(ist_f+num_t,ist_i)+wigd1*wigd2*vLtmp(istar(nf,lf,m1),istar(ni,li,m2))
!            vright(ist_f+num_t,ist_i)=vright(ist_f+num_t,ist_i)+wigd1*wigd2*vRtmp(istar(nf,lf,m1),istar(ni,li,m2))
!          enddo
!        enddo
!        vmatL(ist_f+num_t,ist_i)=vmatL(ist_f+num_t,ist_i)*(-1._id)**(mf+mi)
!        vright(ist_f+num_t,ist_i)=vright(ist_f+num_t,ist_i)*(-1._id)**(mf+mi)   
!
!        
!      enddo
!    enddo
!
!  enddo
!enddo
!!$omp end do
!
!
!!$omp do schedule(guided) 
!do j=1,jdir
!  if(mod(l_p(fdir(j))+l_p(idir(j)),2)==0)then
!    vright(fdir(j),idir(j))=direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(z*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(z*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
!  else
!    vright(fdir(j),idir(j))=-direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(z*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(z*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
!  endif
!  vright(idir(j),fdir(j))=conjg(vright(fdir(j),idir(j)))
!  vmatL(fdir(j),idir(j))=cmplx(0,0,id)
!  vmatL(idir(j),fdir(j))=cmplx(0,0,id)
!enddo
!!$omp end do
!!$omp do schedule(guided) 
!do j=1,num_t
!  vright(j,j)=direct_num_p(j,j,cosang,bigr)*cmplx(cos(z*(k_p(inl_p(j))-k_p(inl_p(j)))),sin(z*(k_p(inl_p(j))-k_p(inl_p(j)))),id)
!  vright(j+num_t,j+num_t)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(z*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(z*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
!  vmatL(j,j)=cmplx(1,0,id)
!  vmatL(j+num_t,j+num_t)=cmplx(1,0,id)
!enddo
!!$omp end do
!!$omp do schedule(guided) 
!do j=1,jdir
!  vright(fdir(j)+num_t,idir(j)+num_t)=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(z*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(z*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
!  vright(idir(j)+num_t,fdir(j)+num_t)=conjg(vright(fdir(j)+num_t,idir(j)+num_t))
!  vmatL(idir(j)+num_t,fdir(j)+num_t)=cmplx(0,0,id)
!  vmatL(fdir(j)+num_t,idir(j)+num_t)=cmplx(0,0,id)
!enddo
!!$omp end do
!
!!!$omp single
!
!!$omp master
!call zhesv('U',2*num_t,2*num_t,vmatL(1:2*num_t,1:2*num_t),2*num_t,ipiv(1:2*num_t),vright(1:2*num_t,1:2*num_t),2*num_t,work,2*num_t*2*num_t,info)
!!$omp end master
!!$omp barrier
!
!!!$omp end single
!
!!!$omp end parallel
!
!
!do i=-nstep,nstep-1
!
!!$omp single
!    z=step(i)
!    znew=step(i+1)
!    zmid=(z+znew)/real(2,id)
!    h=znew-z
!    hh=h/2._id
!    h6=h/6._id
!!$omp end single
!
!
!!$omp do schedule(guided)collapse(2)
!    do k=1,2*num_t 
!      do j=1,2*num_t 
!        vleft(k,j)=vright(k,j)
!      enddo
!    enddo
!!$omp end do
!!$omp do schedule(guided)
!    do j=1,2*num_t 
!      zdadz(j)=-cu*sum(vleft(j,1:2*num_t)*zam(1:2*num_t))
!    enddo
!!$omp end do
!
!!$omp single
!    bigr=sqrt(rho*rho+zmid*zmid)
!    cosang=zmid/bigr
!    ang=acos(zmid/bigr)
!!$omp end single
!
!!$omp do schedule(guided)
!    do j=1,nglag
!       u(j)=glag(j)/bigr+1._id
!    enddo
!!$omp end do
!!$omp do schedule(guided)collapse(2)
!    do j=1,nglag
!      do k=1,ngleg
!        call get_plmij(j,k,rho)
!      enddo
!    enddo
!!$omp end do
!
!!$omp do schedule(guided)
!    do l=1,wigarr
!      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
!    enddo
!!$omp end do
!
!!$omp do schedule(guided)collapse(3)
!  do n=1,nlnum_t
!    do j=1,nglag
!      do k=1,ngleg
!        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
!      enddo
!    enddo
!  enddo
!!$omp end do
!
!!$omp do schedule(guided)collapse(2)
!do j=1,num_t
!  do k=1,num_t
!    call get_LR_mat_dev_pt2(j,k,bigr,zmid,vLtmp(j,k),vRtmp(j,k))
!  enddo
!enddo
!!$omp end do
!
!!$omp do schedule(guided)collapse(2)
!do j=1,nlnum_t
!  do k=1,nlnum_t
!
!
!    nf=nst_t(j)
!    lf=lst_t(j)
!    ni=nst_t(k)
!    li=lst_t(k)
!    do ist_f=lob(j),upb(j)
!      mf=m_t(ist_f)
!      do ist_i=lob(k),upb(k)
!        mi=m_t(ist_i)
!
!        vmatL(ist_f,ist_i+num_t)=zero   
!        vmid(ist_f,ist_i+num_t)=zero   
!        do m1=-min(lf,mmax_t),min(lf,mmax_t)
!          wigd1=wig_d(lf,m1,mf)
!          do m2=-min(li,mmax_t),min(li,mmax_t)
!            wigd2=wig_d(li,m2,mi)
!            vmatL(ist_f,ist_i+num_t)=vmatL(ist_f,ist_i+num_t)+wigd1*wigd2*vLtmp(istar(nf,lf,m1),istar(ni,li,m2))
!            vmid(ist_f,ist_i+num_t)=vmid(ist_f,ist_i+num_t)+wigd1*wigd2*vRtmp(istar(nf,lf,m1),istar(ni,li,m2))
!          enddo
!        enddo
!        vmatL(ist_f,ist_i+num_t)=vmatL(ist_f,ist_i+num_t)*(-1._id)**(lf+li+mf+mi)
!        vmid(ist_f, ist_i+num_t)=vmid(ist_f, ist_i+num_t)*(-1._id)**(lf+li+mf+mi)   
!        
!      enddo
!    enddo
!
!  enddo
!enddo
!!$omp end do
!
!
!
!!$omp do schedule(guided)collapse(3)
!  do n=1,nlnum_t
!    do j=1,nglag
!      do k=1,ngleg
!        call wf_at_r(1,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
!      enddo
!    enddo
!  enddo
!!$omp end do
!
!!$omp do schedule(guided)collapse(2)
!do j=1,num_t
!  do k=1,num_t
!    call get_LR_mat_dev_tp2(j,k,bigr,zmid,vLtmp(j,k),vRtmp(j,k))
!  enddo
!enddo
!!$omp end do
!
!!$omp do schedule(guided)collapse(2)
!do j=1,nlnum_t
!  do k=1,nlnum_t
!
!
!    nf=nst_t(j)
!    lf=lst_t(j)
!    ni=nst_t(k)
!    li=lst_t(k)
!    do ist_f=lob(j),upb(j)
!      mf=m_t(ist_f)
!      do ist_i=lob(k),upb(k)
!        mi=m_t(ist_i)
!        vmatL(ist_f+num_t,ist_i)=zero   
!        vmid(ist_f+num_t,ist_i)=zero   
!        do m1=-min(lf,mmax_t),min(lf,mmax_t)
!          wigd1=wig_d(lf,m1,mf)
!          do m2=-min(li,mmax_t),min(li,mmax_t)
!            wigd2=wig_d(li,m2,mi)
!            vmatL(ist_f+num_t,ist_i)=vmatL(ist_f+num_t,ist_i)+wigd1*wigd2*vLtmp(istar(nf,lf,m1),istar(ni,li,m2))
!            vmid(ist_f+num_t,ist_i)=vmid(ist_f+num_t,ist_i)+wigd1*wigd2*vRtmp(istar(nf,lf,m1),istar(ni,li,m2))
!          enddo
!        enddo
!        vmatL(ist_f+num_t,ist_i)=vmatL(ist_f+num_t,ist_i)*(-1._id)**(mf+mi)
!        vmid(ist_f+num_t,ist_i)=vmid(ist_f+num_t,ist_i)*(-1._id)**(mf+mi)   
!        
!      enddo
!    enddo
!
!  enddo
!enddo
!!$omp end do
!
!
!
!!$omp do schedule(guided) 
!do j=1,jdir
!  if(mod(l_p(fdir(j))+l_p(idir(j)),2)==0)then
!    vmid(fdir(j),idir(j))=direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(zmid*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(zmid*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
!  else
!    vmid(fdir(j),idir(j))=-direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(zmid*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(zmid*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
!  endif
!  vmid(idir(j),fdir(j))=conjg(vmid(fdir(j),idir(j)))
!  vmatL(fdir(j),idir(j))=cmplx(0,0,id)
!  vmatL(idir(j),fdir(j))=cmplx(0,0,id)
!enddo
!!$omp end do
!!$omp do schedule(guided) 
!do j=1,num_t
!  vmid(j,j)=direct_num_p(j,j,cosang,bigr)*cmplx(cos(zmid*(k_p(inl_p(j))-k_p(inl_p(j)))),sin(zmid*(k_p(inl_p(j))-k_p(inl_p(j)))),id)
!  vmid(j+num_t,j+num_t)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(zmid*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(zmid*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
!  vmatL(j,j)=cmplx(1,0,id)
!  vmatL(j+num_t,j+num_t)=cmplx(1,0,id)
!enddo
!!$omp end do
!!$omp do schedule(guided) 
!do j=1,jdir
!  vmid(fdir(j)+num_t,idir(j)+num_t)=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(zmid*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(zmid*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
!  vmid(idir(j)+num_t,fdir(j)+num_t)=conjg(vmid(fdir(j)+num_t,idir(j)+num_t))
!  vmatL(idir(j)+num_t,fdir(j)+num_t)=cmplx(0,0,id)
!  vmatL(fdir(j)+num_t,idir(j)+num_t)=cmplx(0,0,id)
!enddo
!!$omp end do
!!!$omp end parallel 
!
!
!
!!!$omp single
!!$omp master
!
!call zhesv('U',2*num_t,2*num_t,vmatL(1:2*num_t,1:2*num_t),2*num_t,ipiv(1:2*num_t),vmid(1:2*num_t,1:2*num_t),2*num_t,work,2*num_t*2*num_t,info)
!
!    bigr=sqrt(rho*rho+znew*znew)
!    cosang=znew/bigr
!    ang=acos(znew/bigr)
!
!!$omp end master
!!$omp barrier
!!!$omp end single
!
!
!
!!$omp do schedule(guided)
!    do j=1,nglag
!       u(j)=glag(j)/bigr+1._id
!    enddo
!!$omp end do
!!$omp do schedule(guided)collapse(2)
!    do j=1,nglag
!      do k=1,ngleg
!        call get_plmij(j,k,rho)
!      enddo
!    enddo
!!$omp end do
!
!!$omp do schedule(guided)
!    do l=1,wigarr
!      wig_d(larr(l),marr1(l),marr2(l))=wigd(larr(l),marr1(l),marr2(l),-ang)
!    enddo
!!$omp end do
!
!!$omp do schedule(guided)collapse(3)
!  do n=1,nlnum_t
!    do j=1,nglag
!      do k=1,ngleg
!        call wf_at_r(0,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
!      enddo
!    enddo
!  enddo
!!$omp end do
!!$omp do schedule(guided)collapse(2)
!do j=1,num_t
!  do k=1,num_t
!    call get_LR_mat_dev_pt2(j,k,bigr,znew,vLtmp(j,k),vRtmp(j,k))
!  enddo
!enddo
!!$omp end do
!
!!$omp do schedule(guided)collapse(2)
!do j=1,nlnum_t
!  do k=1,nlnum_t
!
!
!    nf=nst_t(j)
!    lf=lst_t(j)
!    ni=nst_t(k)
!    li=lst_t(k)
!    do ist_f=lob(j),upb(j)
!      mf=m_t(ist_f)
!      do ist_i=lob(k),upb(k)
!        mi=m_t(ist_i)
!
!        vmatL(ist_f,ist_i+num_t)=zero   
!        vright(ist_f,ist_i+num_t)=zero   
!        do m1=-min(lf,mmax_t),min(lf,mmax_t)
!          wigd1=wig_d(lf,m1,mf)
!          do m2=-min(li,mmax_t),min(li,mmax_t)
!            wigd2=wig_d(li,m2,mi)
!            vmatL(ist_f,ist_i+num_t)=vmatL(ist_f,ist_i+num_t)+wigd1*wigd2*vLtmp(istar(nf,lf,m1),istar(ni,li,m2))
!            vright(ist_f,ist_i+num_t)=vright(ist_f,ist_i+num_t)+wigd1*wigd2*vRtmp(istar(nf,lf,m1),istar(ni,li,m2))
!          enddo
!        enddo
!        vmatL(ist_f, ist_i+num_t)=vmatL(ist_f ,ist_i+num_t)*(-1._id)**(lf+li+mf+mi)
!        vright(ist_f,ist_i+num_t)=vright(ist_f,ist_i+num_t)*(-1._id)**(lf+li+mf+mi)   
!        
!      enddo
!    enddo
!
!  enddo
!enddo
!!$omp end do
!!$omp do schedule(guided)collapse(3)
!  do n=1,nlnum_t
!    do j=1,nglag
!      do k=1,ngleg
!        call wf_at_r(1,nst_t(n),lst_t(n),bigr*(u(j)+gleg(k))/2._id,bigr*(u(j)-gleg(k))/2._id,basis_a(nst_t(n),lst_t(n),j,k),basis_b(nst_t(n),lst_t(n),j,k),basis_bres(nst_t(n),lst_t(n),j,k))
!      enddo
!    enddo
!  enddo
!!$omp end do
!!$omp do schedule(guided)collapse(2)
!do j=1,num_t
!  do k=1,num_t
!    call get_LR_mat_dev_tp2(j,k,bigr,znew,vLtmp(j,k),vRtmp(j,k))
!  enddo
!enddo
!!$omp end do
!
!!$omp do schedule(guided)collapse(2)
!do j=1,nlnum_t
!  do k=1,nlnum_t
!
!
!    nf=nst_t(j)
!    lf=lst_t(j)
!    ni=nst_t(k)
!    li=lst_t(k)
!    do ist_f=lob(j),upb(j)
!      mf=m_t(ist_f)
!      do ist_i=lob(k),upb(k)
!        mi=m_t(ist_i)
!        vmatL(ist_f+num_t,ist_i)=zero   
!        vright(ist_f+num_t,ist_i)=zero   
!        do m1=-min(lf,mmax_t),min(lf,mmax_t)
!          wigd1=wig_d(lf,m1,mf)
!          do m2=-min(li,mmax_t),min(li,mmax_t)
!            wigd2=wig_d(li,m2,mi)
!            vmatL(ist_f+num_t,ist_i)=vmatL(ist_f+num_t,ist_i)+wigd1*wigd2*vLtmp(istar(nf,lf,m1),istar(ni,li,m2))
!            vright(ist_f+num_t,ist_i)=vright(ist_f+num_t,ist_i)+wigd1*wigd2*vRtmp(istar(nf,lf,m1),istar(ni,li,m2))
!          enddo
!        enddo
!        vmatL(ist_f+num_t,ist_i)=vmatL(ist_f+num_t,ist_i)*(-1._id)**(mf+mi)
!        vright(ist_f+num_t,ist_i)=vright(ist_f+num_t,ist_i)*(-1._id)**(mf+mi)   
!        
!      enddo
!    enddo
!
!  enddo
!enddo
!!$omp end do
!!$omp do schedule(guided) 
!do j=1,jdir
!  if(mod(l_p(fdir(j))+l_p(idir(j)),2)==0)then
!    vright(fdir(j),idir(j))=direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(znew*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(znew*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
!  else
!    vright(fdir(j),idir(j))=-direct_num_p(fdir(j),idir(j),cosang,bigr)*cmplx(cos(znew*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),sin(znew*(k_p(inl_p(idir(j)))-k_p(inl_p(fdir(j))))),id)
!  endif
!  vright(idir(j),fdir(j))=conjg(vright(fdir(j),idir(j)))
!  vmatL(fdir(j),idir(j))=cmplx(0,0,id)
!  vmatL(idir(j),fdir(j))=cmplx(0,0,id)
!enddo
!!$omp end do
!!$omp do schedule(guided) 
!do j=1,num_t
!  vright(j,j)=direct_num_p(j,j,cosang,bigr)*cmplx(cos(znew*(k_p(inl_p(j))-k_p(inl_p(j)))),sin(znew*(k_p(inl_p(j))-k_p(inl_p(j)))),id)
!  vright(j+num_t,j+num_t)=direct_num_t(j,j,cosang,bigr)*cmplx(cos(znew*(k_t(inl_t(j))-k_t(inl_t(j)))),sin(znew*(k_t(inl_t(j))-k_t(inl_t(j)))),id)
!  vmatL(j,j)=cmplx(1,0,id)
!  vmatL(j+num_t,j+num_t)=cmplx(1,0,id)
!enddo
!!$omp end do
!!$omp do schedule(guided) 
!do j=1,jdir
!  vright(fdir(j)+num_t,idir(j)+num_t)=direct_num_t(fdir(j),idir(j),cosang,bigr)*cmplx(cos(znew*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),sin(znew*(k_t(inl_t(idir(j)))-k_t(inl_t(fdir(j))))),id)
!  vright(idir(j)+num_t,fdir(j)+num_t)=conjg(vright(fdir(j)+num_t,idir(j)+num_t))
!  vmatL(idir(j)+num_t,fdir(j)+num_t)=cmplx(0,0,id)
!  vmatL(fdir(j)+num_t,idir(j)+num_t)=cmplx(0,0,id)
!enddo
!!$omp end do
!
!
!
!!!$omp end parallel 
!
!!$omp do schedule(guided)collapse(2)
!    do k=1,2*num_t 
!      do j=1,2*num_t 
!        vover(k,j)=vmatL(k,j)
!      enddo
!    enddo
!!$omp end do
!
!!!$omp single
!
!!$omp master
!    call zhesv('U',2*num_t,2*num_t,vmatL(1:2*num_t,1:2*num_t),2*num_t,ipiv(1:2*num_t),vright(1:2*num_t,1:2*num_t),2*num_t,work,2*num_t*2*num_t,info)
!!$omp end master
!!$omp barrier
!
!!!$omp end single
!
!!$omp do schedule(guided)
!   do j=1,2*num_t 
!     dyt(j)=-cu*sum(vmid(j,1:2*num_t)*(zam(1:2*num_t)+hh*zdadz(1:2*num_t)))
!   enddo
!!$omp end do
!!$omp do schedule(guided)
!   do j=1,2*num_t 
!     dym(j)=-cu*sum(vmid(j,1:2*num_t)*(zam(1:2*num_t)+hh*dyt(1:2*num_t)))
!   enddo
!!$omp end do
!!$omp do schedule(guided)
!   do j=1,2*num_t 
!     dymm(j)=dyt(j)+dym(j)
!     dyt(j)=-cu*sum(vright(j,1:2*num_t)*(zam(1:2*num_t)+h*dym(1:2*num_t)))
!     zout(j)=zam(j)+h6*(zdadz(j)+dyt(j)+2._id*dymm(j))
!   enddo
!!$omp end do 
!!$omp do schedule(guided)
!   do j=1,2*num_t 
!     zam(j)=zout(j)
!   enddo
!!$omp end do 
!
!
!enddo
!
!
!!$omp single 
!    unitarity2=0.0_id
!!$omp end single
!
!
!!$omp do schedule(guided)collapse(2)
!    do k=1,2*num_t 
!      do j=1,2*num_t 
!        unitarity2=unitarity2+conjg(zam(k))*zam(j)*vover(k,j)
!      enddo
!    enddo
!!$omp end do
!
!!$omp single 
!    zam(num_t+i0)=zam(num_t+i0)-1._id
!!$omp end single
!
!
!!$omp do schedule(guided)
!    do ist_f=1,2*num_t
!       prob(ist_f)=zam(ist_f)*conjg(zam(ist_f))*rho
!    enddo
!!$omp end do
!!$omp end parallel
!
!    write(1298,*)rho,dble(unitarity2)
!
!      call system_clock(t2)
!      dt = t2-t1
!      secs = real(dt)/real(count_rate)
!
!      open(1299,file='time',status='unknown')
!      write(1299,"('number of channels is ',i5,' wall clock time is ',f12.2,' seconds')") num_t,secs
!      open(1300,file='timing',status='unknown')
!      write(1300,"(i5,es16.6)") num_t,secs
!
! 30   format(i20,10000g14.6)
!
!end subroutine solveRK_num_2c

subroutine get_jint(rho)
   use precision_type
   use data_base
   implicit none
   real(kind=id)::rho,arg,rj,ry,rjp,ryp
   integer::m,i,j

!$omp parallel&
!$omp default(none)&
!$omp private(i,j,m,arg,rj,ry,rjp,ryp)&
!$omp shared(nglag,ngleg,rho,v,glag,gleg,jint,mmax_t)
!$omp do&
!$omp schedule(dynamic)
   do i=1,nglag
      do j=1,ngleg
         arg=rho*v/2._id
         arg=arg*sqrt(((glag(i)+1._id)*(glag(i)+1._id)-1._id)*(1._id-gleg(j)*gleg(j)))
         do m=0,2*mmax_t
            call bessjy(arg,real(m,id),rj,ry,rjp,ryp)

            jint(m,i,j)=rj
            if(mod(m,2)==0)then
              jint(-m,i,j)=jint(m,i,j)
            else
              jint(-m,i,j)=-jint(m,i,j)
            endif
         enddo
      enddo
   enddo
!$omp end do
!$omp end parallel

!  do i=1,nglag
!    do j=1,ngleg
!      do m=-2*mmax_t,2*mmax_t
!         write(606,*)i,j,m,dble(jint(m,i,j))
!      enddo
!    enddo
!  enddo

end subroutine get_jint

subroutine bessjy(x,xnu,rj,ry,rjp,ryp)
  use precision_type,only:id,mppic
  implicit none
  real(kind=id):: eps,fpmin,rj,rjp,ry,ryp,x,xnu,xmin
  !u    uses beschb
  integer(kind=4):: i,isign,l,nl,maxit
  real(kind=id):: a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,temp,w,x2,xi,xi2,xmu,xmu2



  if(id==16)then
    maxit=1000000000
    eps=1.e-50_id
    fpmin=1.e-100_id
    xmin=2._id
  else
    maxit=1000000
    eps=1.e-16_id
    fpmin=1.e-16_id
    xmin=2._id
  endif


  if(x.le.0._id.or.xnu.lt.0._id) pause 'bad arguments in bessjy'
  if(x.lt.xmin)then
    nl=int(xnu+0.5_id)
  else
    nl=max(0,int(xnu-x+1.5_id))
  endif
  xmu=xnu-nl
  xmu2=xmu*xmu
  xi=1._id/x
  xi2=2._id*xi
  w=xi2/mppic
  isign=1
  h=xnu*xi
  if(h.lt.fpmin)h=fpmin
  b=xi2*xnu
  d=0._id
  c=h
  do 11 i=1,maxit
    b=b+xi2
    d=b-d
    if(abs(d).lt.fpmin)d=fpmin
    c=b-1._id/c
    if(abs(c).lt.fpmin)c=fpmin
    d=1._id/d
    del=c*d
    h=del*h
    if(d.lt.0._id)isign=-isign
    if(abs(del-1._id).lt.eps)goto 1
11    continue
  pause 'x too large in bessjy; try asymptotic expansion'
1     continue
  rjl=isign*fpmin
  rjpl=h*rjl
  rjl1=rjl
  rjp1=rjpl
  fact=xnu*xi
  do 12 l=nl,1,-1
    rjtemp=fact*rjl+rjpl
    fact=fact-xi
    rjpl=fact*rjtemp-rjl
    rjl=rjtemp
12    continue
  if(rjl.eq.0._id)rjl=eps
  f=rjpl/rjl
  if(x.lt.xmin) then
    x2=0.5_id*x
    pimu=mppic*xmu
    if(abs(pimu).lt.eps)then
      fact=1._id
    else
      fact=pimu/sin(pimu)
    endif
    d=-log(x2)
    e=xmu*d
    if(abs(e).lt.eps)then
      fact2=1._id
    else
      fact2=sinh(e)/e
    endif
    call beschb(xmu,gam1,gam2,gampl,gammi)
    ff=2._id/mppic*fact*(gam1*cosh(e)+gam2*fact2*d)
    e=exp(e)
    p=e/(gampl*mppic)
    q=1._id/(e*mppic*gammi)
    pimu2=0.5_id*pimu
    if(abs(pimu2).lt.eps)then
      fact3=1._id
    else
      fact3=sin(pimu2)/pimu2
    endif
    r=mppic*pimu2*fact3*fact3
    c=1._id
    d=-x2*x2
    sum=ff+r*q
    sum1=p
    do 13 i=1,maxit
      ff=(i*ff+p+q)/(i*i-xmu2)
      c=c*d/i
      p=p/(i-xmu)
      q=q/(i+xmu)
      del=c*(ff+r*q)
      sum=sum+del
      del1=c*p-i*del
      sum1=sum1+del1
      if(abs(del).lt.(1._id+abs(sum))*eps)goto 2
13      continue
    pause 'bessy series failed to converge'
2       continue
    rymu=-sum
    ry1=-sum1*xi2
    rymup=xmu*xi*rymu-ry1
    rjmu=w/(rymup-f*rymu)
  else
    a=0.25_id-xmu2
    p=-xi/2._id
    q=1._id
    br=2._id*x
    bi=2._id
    fact=a*xi/(p*p+q*q)
    cr=br+q*fact
    ci=bi+p*fact
    den=br*br+bi*bi
    dr=br/den
    di=-bi/den
    dlr=cr*dr-ci*di
    dli=cr*di+ci*dr
    temp=p*dlr-q*dli
    q=p*dli+q*dlr
    p=temp
    do 14 i=2,maxit
      a=a+2*(i-1)
      bi=bi+2._id
      dr=a*dr+br
      di=a*di+bi
      if(abs(dr)+abs(di).lt.fpmin)dr=fpmin
      fact=a/(cr*cr+ci*ci)
      cr=br+cr*fact
      ci=bi-ci*fact
      if(abs(cr)+abs(ci).lt.fpmin)cr=fpmin
      den=dr*dr+di*di
      dr=dr/den
      di=-di/den
      dlr=cr*dr-ci*di
      dli=cr*di+ci*dr
      temp=p*dlr-q*dli
      q=p*dli+q*dlr
      p=temp
      if(abs(dlr-1._id)+abs(dli).lt.eps)goto 3
14  continue
       pause 'cf2 failed in bessjy'
3   continue
      gam=(p-f)/q
      rjmu=sqrt(w/((p-f)*gam+q))
      rjmu=sign(rjmu,rjl)
      rymu=rjmu*gam
      rymup=rymu*(p+q/gam)
      ry1=xmu*xi*rymu-rymup
    endif
    fact=rjmu/rjl
    rj=rjl1*fact
    rjp=rjp1*fact
    do 15 i=1,nl
      rytemp=(xmu+i)*xi2*ry1-rymu
      rymu=ry1
      ry1=rytemp
15  continue
    ry=rymu
    ryp=xnu*xi*rymu-ry1
return
end

subroutine beschb(x,gam1,gam2,gampl,gammi)
  use precision_type
  implicit none
 ! integer(kind=4):: nuse1,nuse2
  real(kind=id):: gam1,gam2,gammi,gampl,x!,xx,chebev

  gam1=-0.577215664901532860606512090082402431042159335939923598805767234884867726777664670936947063291746750_id
  gam2=1._id
  gampl=gam2-x*gam1
  gammi=gam2+x*gam1
return
end






subroutine formfactor_t
use precision_type,only:id
use mpi,only:myid
use data_base
use acc_routines
implicit none
integer(kind=4)::ir,lmd,ni,li,nj,lj,jmax,ilmd,ilmd0,i,j,k,ub
!real(kind=id),dimension(1:meshr,0:2*lmax_t)::rpow1,rpow2
real(kind=id),dimension(1:meshr)::ftmp
real(kind=id),dimension(1:meshr)::fun
real(kind=id),dimension(1:meshr)::tmp
integer(kind=4),allocatable,dimension(:)::karr, lmdarr

allocate(rpow1(1:meshr,0:2*lmax_t))
allocate(rpow2(1:meshr,0:2*lmax_t))

do ir=1,meshr
   rpow1(ir,0)=1._id
   rpow2(ir,0)=1._id/rmesh(ir,1)
enddo

do lmd=1,2*lmax_t
   do ir=1,meshr
      rpow1(ir,lmd)=rpow1(ir,lmd-1)*rmesh(ir,1)
      rpow2(ir,lmd)=rpow2(ir,lmd-1)/rmesh(ir,1)
   enddo
enddo

!$omp parallel&
!$omp default(none)&
!$omp private(ni,li,nj,lj,ub)&
!$omp shared(jdirnl,fdirnl,idirnl,nst_t,lst_t,meshr,formf_tr)
!$omp do&
!$omp schedule(dynamic)
      do k=1,jdirnl
         ni=nst_t(fdirnl(k))
         li=lst_t(fdirnl(k))
         nj=nst_t(idirnl(k))
         lj=lst_t(idirnl(k))
         ub=(li+lj-iabs(li-lj))/2 
         call formAll(ni,li,nj,lj,formf_tr(k,0:ub,1:meshr))
      enddo
!$omp end do
!$omp end parallel


if(myid==0)then
print*,'size of formf_tr:',sizeof(formf_tr)*1.e-6_id, 'Mb'
print*,'st size of formf_tr:',storage_size(formf_tr)*1.e-6_id/8._id*(jdirnl)*(lmax_t+1)*(meshr), 'Mb'
print*,'st size of formf_tr:',storage_size(formf_tr),(jdirnl)*(lmax_t+1)*(meshr)
endif

!$acc update device(formf_tr)

return
end

subroutine formfactor_p
use precision_type,only:id
use mpi,only:myid
use data_base
use acc_routines
implicit none
integer(kind=4)::ir,lmd,ni,li,nj,lj,jmax,ilmd,ilmd0,i,j,k,ub
!real(kind=id),dimension(1:meshr,0:2*lmax_t)::rpow1,rpow2
real(kind=id),dimension(1:meshr)::ftmp
real(kind=id),dimension(1:jdirnl,1:meshr)::funint
real(kind=id),dimension(0:lmax_t,1:meshr)::tmp
integer(kind=4),allocatable,dimension(:)::karr, lmdarr

!$omp parallel&
!$omp default(none)&
!$omp private(ni,li,nj,lj,ub)&
!$omp shared(jdirnl,fdirnl,idirnl,nst_t,lst_t,meshr,formf_pr)
!$omp do&
!$omp schedule(dynamic)
      do k=1,jdirnl
         ni=nst_t(fdirnl(k))
         li=lst_t(fdirnl(k))
         nj=nst_t(idirnl(k))
         lj=lst_t(idirnl(k))
         ub=(li+lj-iabs(li-lj))/2 
         call formAllp(ni,li,nj,lj,formf_pr(k,0:ub,1:meshr))
      enddo
!$omp end do
!$omp end parallel

if(myid==0)then
print*,'size of formf_pr:',sizeof(formf_pr)*1.e-6_id, 'Mb'
print*,'st size of formf_pr:',storage_size(formf_pr)*1.e-6_id/8._id*(jdirnl)*(lmax_t+1)*(meshr), 'Mb'
print*,'st size of formf_pr:',storage_size(formf_pr),(jdirnl)*(lmax_t+1)*(meshr)
endif

!$acc update device(formf_pr)

return
end

subroutine gauleg(x1,x2,x,w,n)

  use precision_type
!  use data_base,only:num_dig
  implicit none
  real(kind=id) x1,x2
  real(kind=id),dimension(n):: x,w
  real(kind=id),dimension(2*n):: wf
  integer,dimension(2*n):: iwf
  integer::m,i,j,n,ier
  integer,dimension(2*n)::iwf
  real(kind=id)::eps,z,z1,p1,p2,p3,pp,xl,xm,tmp
  

  call cgqf(n,x,w,1,0._id,0._id,x1,x2,0,2*n,wf,2*n,iwf,ier)

  print*,'ier gauleg',ier

  if (ier.eq.0) return

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









subroutine get_faclog

  use precision_type
  use flogs,only:fac10,faccoef,hat,fac,faclog,lm_fac
  use data_base,only:lmax_t,nneg
  implicit none
  integer::i,n,mf,mi,ib,i1,i2,i3,li,j,si,ki
  real(kind=id):: q,gtmp

  allocate(fac10(0:4*lmax_t+1))
  allocate(faclog(0:4*lmax_t))
  allocate(faccoef(0:max(2*lmax_t,nneg),0:2*lmax_t))
  allocate(hat(0:lmax_t))
  allocate(fac(0:4*lmax_t))
  allocate(lm_fac(0:lmax_t,-lmax_t:lmax_t))

  print*,'0'

  faccoef(0,0)=1._id

  do i=1,max(2*lmax_t,nneg)
     gtmp=1._id
     faccoef(i,0)=gtmp
     do j=1,min(i,2*lmax_t)
        gtmp=gtmp/real(i+1-j,id)/real(i+j,id)
        faccoef(i,j)=gtmp
     enddo
  enddo

  print*,'1'

  fac(0)=1._p16
  fac(1)=1._p16
  q = 2._p16
  do n=2,4*lmax_t
    fac(n)=fac(n-1)*n
    q = q + real(1,id)
  enddo

  fac10(0)=1._p16
  fac10(1)=0.1_p16
  q = 2._p16
  do n=2,4*lmax_t+1
    fac10(n) = fac10(n-1) * q / real(10,id)
    q = q + real(1,id)
  enddo


  print*,'2'

  do li=0,lmax_t
    hat(li)=sqrt(real(2*li+1,id))
  enddo
  print*,'3'

  faclog(0)=0._p16
  faclog(1)=0._p16
  do 10 n=2,4*lmax_t
10 faclog(n)=faclog(n-1)+log(real(n,p16))

  do li=0,lmax_t
    do mi=-li,li
      lm_fac(li,mi)=sqrt(real(2*li+1,id))*sqrt(fac(li-mi)/fac(li+mi))
    enddo
  enddo

end subroutine get_faclog



