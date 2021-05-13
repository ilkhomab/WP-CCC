subroutine get_ddcs_fixed(keyp,kf,ki,n,it)
   use precision_type
   use data_base
   implicit none
   integer::keyp,n,ip,it
   complex(kind=id),dimension(0:npt):: tdcs_av_x,tdcs_av_d,tdcs_av_coh 
   real(kind=id):: ddcs_scat_x,ddcs_scat_d,ddcs_scat_coh 
   real(kind=id):: kf,ki 

tdcs_av_x=zero;tdcs_av_d=zero;tdcs_av_coh=zero

!!$omp parallel&
!!$omp default(none)&
!!$omp private(tdcs_av_x,tdcs_av_d,tdcs_av_coh)&
!!$omp shared(npt,kf,ki,n,ddcs_scat_x,ddcs_scat_d,ddcs_scat_coh,it)
!!$omp do&         
!!$omp schedule(dynamic)
do ip=1,npt
      print*,'ip=',ip  
      call get_tdcs_phi_fixed(0,kf,ki,n,ip,it,tdcs_av_x(ip),tdcs_av_d(ip),tdcs_av_coh(ip))
enddo 
!!$omp end do         
!!$omp end parallel

ddcs_scat_x=0d0
ddcs_scat_d=0d0
ddcs_scat_coh=0d0

do ip=1,npt

  ddcs_scat_x=ddcs_scat_x+tdcs_av_x(ip)*pt(ip)*wpt(ip)
  ddcs_scat_d=ddcs_scat_d+tdcs_av_d(ip)*pt(ip)*wpt(ip)
  ddcs_scat_coh=ddcs_scat_coh+tdcs_av_coh(ip)*pt(ip)*wpt(ip)

enddo

print*, 'ddcs_scat_x=',ddcs_scat_x
print*, 'ddcs_scat_d=',ddcs_scat_d
print*, 'ddcs_scat_coh=',ddcs_scat_coh

ddcs_scat_x=2d0*pi*ddcs_scat_x/kf/ki
ddcs_scat_d=2d0*pi*ddcs_scat_d/kf/ki
ddcs_scat_coh=2d0*pi*ddcs_scat_coh/kf/ki

   if(keyp==1)then
      write(31,11)en_t(0)%enl_t(n)*hr,ddcs_scat_x,ddcs_scat_d,ddcs_scat_coh
!      write(31,11)en_t(0)%enl_t(n)*hr,dble(tdcs_av_x(ip)),dble(tdcs_av_d(ip)),dble(tdcs_av_coh(ip))
   endif

11 format(5000000es18.8)


end subroutine get_ddcs_fixed


subroutine get_ddcs_scat3(keyp,kf,ki,n,ddcs_scat_x,ddcs_scat_d,ddcs_scat_coh)
   use precision_type
   use data_base
   implicit none
   integer::keyp,n,ip,it
   complex(kind=id),dimension(0:360):: tdcs_av_x,tdcs_av_d,tdcs_av_coh 
   real(kind=id),dimension(0:npt):: ddcs_scat_x,ddcs_scat_d,ddcs_scat_coh 
   real(kind=id):: kf,ki 

tdcs_av_x(:)=zero;tdcs_av_d(:)=zero;tdcs_av_coh(:)=zero

!$omp parallel&
!$omp default(none)&
!$omp private(tdcs_av_x,tdcs_av_d,tdcs_av_coh)&
!$omp shared(npt,kf,ki,n,ddcs_scat_x,ddcs_scat_d,ddcs_scat_coh)
!$omp do&         
!$omp schedule(dynamic)
   do ip=0,npt
      print*,'ip=',ip  
      call get_tdcs_phi_av3(0,kf,ki,n,ip,tdcs_av_x,tdcs_av_d,tdcs_av_coh)
      ddcs_scat_x(ip)=0d0
      do it=0,180
         ddcs_scat_x(ip)=ddcs_scat_x(ip)+dble(tdcs_av_x(it))*dsin(it*pi/180d0)*pi/180d0
      enddo
      ddcs_scat_d(ip)=0d0
      do it=0,180
         ddcs_scat_d(ip)=ddcs_scat_d(ip)+dble(tdcs_av_d(it))*dsin(it*pi/180d0)*pi/180d0
      enddo
      ddcs_scat_coh(ip)=0d0
      do it=0,180
         ddcs_scat_coh(ip)=ddcs_scat_coh(ip)+dble(tdcs_av_coh(it))*dsin(it*pi/180d0)*pi/180d0
      enddo
   enddo 
!$omp end do         
!$omp end parallel

   if(keyp==1)then
      open(31,file='DDCS_SCAT',status='replace')
      do ip=0,npt
        write(31,11)2d0*asin(pt(ip)/2d0/k0),ddcs_scat_x(ip),ddcs_scat_d(ip),ddcs_scat_coh(ip)
      enddo
      close(31)
   endif

11 format(5000000es18.8)


end subroutine get_ddcs_scat3


subroutine get_tdcs_phi_av3(keyp,kf,ki,n,ip,tdcs_x_av,tdcs_d_av,tdcs_coh_av)
   use precision_type
   use data_base
   use flogs,only:fac
   implicit none
   integer::  ip,keyp,it,l,m,i,j,n,mdif,mdifa,ir,l1,l2,jc,jl,jr,iphi,phi_fix
   complex(kind=id):: ztmp,coulphase,zsumm,zex,zamp 
   complex(kind=id),dimension(1:nbin,0:lmax_t,-lmax_t:lmax_t):: zp
   complex(kind=id),dimension(0:360):: amp_x,tdcs_x_av,tdcs_d_av,tdcs_coh_av 
   complex(kind=id),dimension(0:360,0:360):: amp_d
   complex(kind=id),dimension(1:nbin):: zintr
   real(kind=id):: kappa,ovrlp,kf,ki,costh,rylm,kappa_p,btmp,wtmp,factor,bessj,rylmp,ptr,ktr,phi,res

   kappa=sqrt(2._id*en_t(0)%enl_t(n))
   ovrlp=1d0/dis(n-nmax_t(0)+nbin)

!phi_fix=270
!The direct scattering amplitude:

   do it=180,360
     tdcs_d_av(it)=zero  
     do iphi=0,360!phi_fix,phi_fix!360
        phi=dble(iphi)*pi/180d0 
        costh=dcos((it-180)*pi/180d0)
        if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
        amp_d(it,iphi)=dcmplx(0d0,0d0)
        do l=0,lmax_t
          ztmp=dcmplx(0d0,-1d0)
          ztmp=ztmp**l*coulphase(-1d0/kappa,l)
          ztmp=ztmp*dsqrt(2*l+1d0)
          zsumm=dcmplx(0d0,0d0)
          do m=-l,l
            j=istar(n-l,l,m)+num_t
            zsumm=zsumm+rylm(l,m,costh)*zap(j,ip)*dcmplx(0d0,1d0)**(m)*exp(cu*m*phi)!*dcmplx(dcos(m*phi),dsin(m*phi))!*(-1d0)**m
          enddo
          amp_d(it,iphi)=amp_d(it,iphi)+ztmp*zsumm*(-1d0)**l
        enddo
        amp_d(it,iphi)=amp_d(it,iphi)*sqrt(kf*ki*ovrlp/4d0/pi/kappa)
        tdcs_d_av(it)=tdcs_d_av(it)+amp_d(it,iphi)*conjg(amp_d(it,iphi))
     enddo
     tdcs_d_av(it)=tdcs_d_av(it)*pi/180d0
   enddo
   do it=0,180
     tdcs_d_av(it)=zero  
     do iphi=0,360!phi_fix,phi_fix!360
        phi=dble(iphi)*pi/180d0 
        costh=dcos((it)*pi/180d0)
        if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
        amp_d(it,iphi)=dcmplx(0d0,0d0)
        do l=0,lmax_t
          ztmp=dcmplx(0d0,-1d0)
          ztmp=ztmp**l*coulphase(-1d0/kappa,l)
          ztmp=ztmp*dsqrt(2*l+1d0)
          zsumm=dcmplx(0d0,0d0)
          do m=-l,l
            j=istar(n-l,l,m)+num_t
            zsumm=zsumm+rylm(l,m,costh)*zap(j,ip)*dcmplx(0d0,1d0)**(m)*exp(cu*m*phi)!*dcmplx(dcos(m*phi),dsin(m*phi))
          enddo
          amp_d(it,iphi)=amp_d(it,iphi)+ztmp*zsumm
        enddo
        amp_d(it,iphi)=amp_d(it,iphi)*sqrt(kf*ki*ovrlp/4d0/pi/kappa)
        tdcs_d_av(it)=tdcs_d_av(it)+amp_d(it,iphi)*conjg(amp_d(it,iphi))
     enddo
     tdcs_d_av(it)=tdcs_d_av(it)*pi/180d0
   enddo

!Rearrangement:

   do it=0,180
     ktr=kappa*dsin((it)*pi/180d0)
!     ptr=pt(ip)
     kappa_p=sqrt(kappa**2+v**2-2*kappa*v*dcos((it)*pi/180d0))

     costh=dcos((it)*pi/180d0)
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     tdcs_x_av(it)=zero  
     tdcs_coh_av(it)=zero  
     do iphi=0,360!phi_fix,phi_fix !360
        phi=dble(iphi)*pi/180d0 
        ptr=dsqrt(pt(ip)**2-2d0*pt(ip)*ktr*dcos(phi)+ktr**2)!dble(pt(ip)-ktr)
        amp_x(it)=dcmplx(0d0,0d0)
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
                    zex=zex+factor*wtmp*btmp*bessj(mdifa,btmp*ptr)*zamp
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
          amp_x(it)=amp_x(it)+ztmp*zsumm
        enddo
        amp_x(it)=amp_x(it)*sqrt(kf*ki/4d0/pi*kappa)/kappa_p
        tdcs_x_av(it)=tdcs_x_av(it)+amp_x(it)*conjg(amp_x(it))
        tdcs_coh_av(it)=tdcs_coh_av(it)+(amp_x(it)+amp_d(it,iphi))*conjg(amp_x(it)+amp_d(it,iphi))
     enddo
     tdcs_x_av(it)=tdcs_x_av(it)*pi/180d0
     tdcs_coh_av(it)=tdcs_coh_av(it)*pi/180d0
   enddo




   do it=180,360
     ktr=kappa*dsin((it)*pi/180d0)
!     ptr=dabs(pt(ip)-ktr)
     kappa_p=sqrt(kappa**2+v**2-2*kappa*v*dcos((it)*pi/180d0))

     costh=dcos((it-180)*pi/180d0)
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     tdcs_x_av(it)=zero  
     tdcs_coh_av(it)=zero  
     do iphi=0,360!phi_fix,phi_fix!0,360
        phi=dble(iphi)*pi/180d0 
        ptr=dsqrt(pt(ip)**2-2d0*pt(ip)*ktr*dcos(phi)+ktr**2)!dble(pt(ip)-ktr)
        amp_x(it)=dcmplx(0d0,0d0)
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
                    zex=zex+factor*wtmp*btmp*bessj(mdifa,btmp*ptr)*zamp
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
          amp_x(it)=amp_x(it)+ztmp*zsumm
        enddo
        amp_x(it)=amp_x(it)*sqrt(kf*ki/4d0/pi*kappa)/kappa_p
        tdcs_x_av(it)=tdcs_x_av(it)+amp_x(it)*conjg(amp_x(it))
        tdcs_coh_av(it)=tdcs_coh_av(it)+(amp_x(it)+amp_d(it,iphi))*conjg(amp_x(it)+amp_d(it,iphi))
     enddo
     tdcs_x_av(it)=tdcs_x_av(it)*pi/180d0
     tdcs_coh_av(it)=tdcs_coh_av(it)*pi/180d0
   enddo


!Printing:

   if(keyp==1)then
      open(31,file='TDCS_AV2',status='replace')
      do it=0,360
        write(31,11)dble(it),tdcs_x_av(it),tdcs_d_av(it),tdcs_coh_av(it),en_t(0)%enl_t(n)*hr
      enddo
      close(31)
   endif

11 format(5000000es18.8)

end subroutine get_tdcs_phi_av3 



subroutine get_tdcs_phi_fixed(keyp,kf,ki,n,ip,it,tdcs_x_av,tdcs_d_av,tdcs_coh_av)
   use precision_type
   use data_base
   use flogs,only:fac
   implicit none
   integer::  ip,keyp,it,l,m,i,j,n,mdif,mdifa,ir,l1,l2,jc,jl,jr,iphi,phi_fix
   complex(kind=id):: ztmp,coulphase,zsumm,zex,zamp 
   complex(kind=id),dimension(1:nbin,0:lmax_t,-lmax_t:lmax_t):: zp
   complex(kind=id):: amp_x,tdcs_x_av,tdcs_d_av,tdcs_coh_av 
   complex(kind=id),dimension(0:360):: amp_d
   complex(kind=id),dimension(1:nbin):: zintr
   real(kind=id):: kappa,ovrlp,kf,ki,costh,rylm,kappa_p,btmp,wtmp,factor,bessj,rylmp,ptr,ktr,phi,res

   kappa=sqrt(2._id*en_t(0)%enl_t(n))
   ovrlp=1d0/dis(n-nmax_t(0)+nbin)

!phi_fix=270
!The direct scattering amplitude:

     tdcs_d_av=zero  
     do iphi=0,360!phi_fix,phi_fix!360
        phi=dble(iphi)*pi/180d0 
        costh=dcos((it-180)*pi/180d0)
        if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
        amp_d(iphi)=dcmplx(0d0,0d0)
        do l=0,lmax_t
          ztmp=dcmplx(0d0,-1d0)
          ztmp=ztmp**l*coulphase(-1d0/kappa,l)
          ztmp=ztmp*dsqrt(2*l+1d0)
          zsumm=dcmplx(0d0,0d0)
          do m=-l,l
            j=istar(n-l,l,m)+num_t
            zsumm=zsumm+rylm(l,m,costh)*zap(j,ip)*dcmplx(0d0,1d0)**(m)*exp(cu*m*phi)!*dcmplx(dcos(m*phi),dsin(m*phi))!*(-1d0)**m
          enddo
          amp_d(iphi)=amp_d(iphi)+ztmp*zsumm*(-1d0)**l
        enddo
        amp_d(iphi)=amp_d(iphi)*sqrt(kf*ki*ovrlp/4d0/pi/kappa)
        tdcs_d_av=tdcs_d_av+amp_d(iphi)*conjg(amp_d(iphi))
     enddo
     tdcs_d_av=tdcs_d_av*pi/180d0

!Rearrangement:

     ktr=kappa*dsin((it)*pi/180d0)
!     ptr=dabs(pt(ip)-ktr)
     kappa_p=sqrt(kappa**2+v**2-2*kappa*v*dcos((it)*pi/180d0))

     costh=dcos((it-180)*pi/180d0)
     if(dabs(costh).gt.1d0)costh=costh/dabs(costh)
     tdcs_x_av=zero  
     tdcs_coh_av=zero  
     do iphi=0,360!phi_fix,phi_fix!0,360
        phi=dble(iphi)*pi/180d0 
        ptr=dsqrt(pt(ip)**2-2d0*pt(ip)*ktr*dcos(phi)+ktr**2)!dble(pt(ip)-ktr)
        amp_x=dcmplx(0d0,0d0)
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

                if(kappa_p.lt.bin(nbin))then
                  call zintrpl(nbin,bin(1:nbin),zintr(1:nbin),kappa_p,zamp)
                else
                  zamp=zero
                endif

                if(dabs(ptr).lt.1d-4)then
                  if(mdifa==0)then
                    zex=zex+factor*wtmp*btmp*zamp
                  else
                    zex=dcmplx(0d0,0d0)
                  endif
                else
                    zex=zex+factor*wtmp*btmp*bessj(mdifa,btmp*ptr)*zamp
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
          amp_x=amp_x+ztmp*zsumm
        enddo
        amp_x=amp_x*sqrt(kf*ki/4d0/pi*kappa)/kappa_p
        tdcs_x_av=tdcs_x_av+amp_x*conjg(amp_x)
        tdcs_coh_av=tdcs_coh_av+(amp_x+amp_d(iphi))*conjg(amp_x+amp_d(iphi))
     enddo
     tdcs_x_av=tdcs_x_av*pi/180d0
     tdcs_coh_av=tdcs_coh_av*pi/180d0


!Printing:

   if(keyp==1)then
      open(31,file='TDCS_AV_FIXED',status='replace')
        write(31,11)dble(it),tdcs_x_av,tdcs_d_av,tdcs_coh_av,en_t(0)%enl_t(n)*hr
      close(31)
   endif

11 format(5000000es18.8)

end subroutine get_tdcs_phi_fixed


