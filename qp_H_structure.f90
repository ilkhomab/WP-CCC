subroutine get_ci_coef_and_energy_p16(z0,nmax,l,alf,c,enrg,wfar,hwfar)

  use precision_type,only:p16,id
  use flogs
  use data_base,only:rmesh,meshr
  implicit none
  integer::l,nmax,nfac,i,n,k
  real(kind=p16):: alf,summa,wf
  real(kind=id):: z0
  real(kind=p16), dimension(1:nmax):: enrg,dnorm,dnorm2
  real(kind=p16), dimension(1:nmax,1:nmax):: c, h
  real(kind=p16), dimension(1:8*nmax):: work
  real(kind=id), dimension(1:nmax,1:meshr):: wfar,dwfar,hwfar
  integer,dimension(1:5*nmax):: iwork
  real(kind=p16), dimension(0:2*nmax+2*l):: facl
!  z0=real(1,p16)

  real*8:: r, r2
  integer:: ir
  real*8, allocatable:: sum(:), norm_sum(:)

  nfac=2*nmax+2*l

  call basisl_p16(real(z0,p16),l,alf,nmax,nmax,enrg,c,h,dnorm,dnorm2,work,iwork,facl,nfac)


  do n=1,nmax

    summa=real(0,p16)
    do i=1,nmax

!      print*,i
      c(n,i)=c(n,i)*dnorm(n)

      summa=summa+c(n,i)*c(n,i)*fac(2*l+1+n)/real(2,p16)/fac(n-1)


    enddo


      write(66,*)'sum c for n=',n,l

      write(66,*)summa


  enddo

  do n=1,nmax


    call comp_wf_b_p16(n,l,1e-4_p16,nmax,c,alf,wf)


    if(wf.lt.0._p16)then
        print*,'Change of sign occured for n,l',n,l

      do k=1,nmax
        c(k,n)=-c(k,n)
      enddo

    endif

  enddo

  do i=1,meshr
    call comp_wf_b_ar_p16(l,rmesh(i,1),nmax,c(1:nmax,1:nmax),alf,wfar(1:nmax,i))
    call comp_dwf_ar_p16(l,rmesh(i,1),nmax,c(1:nmax,1:nmax),alf,real(z0,p16),real(wfar(1:nmax,i),p16),hwfar(1:nmax,i))
  enddo

  allocate( norm_sum(1:nmax), sum(1:nmax) )


  print*, ' Charge, Z = ', dble(z0), 'alf = ', dble(alf)

  norm_sum(1:nmax)=0._id
  do ir=1,meshr
     r=sqrt(rmesh(ir,3))
     r2=rmesh(ir,1)**2
     norm_sum(1:nmax)=norm_sum(1:nmax)+r2*r*r*wfar(1:nmax,ir)*wfar(1:nmax,ir)
  enddo

  sum(1:nmax)=0._id
  do ir=1,meshr
     r=sqrt(rmesh(ir,3))
     r2=rmesh(ir,1)**2
     !sum(1:nbin)=sum(1:nbin)+r*r*wf(1:nbin,l,ir)*wfres(1:nbin,l,ir)
     sum(1:nmax)=sum(1:nmax)+r2*r*r*wfar(1:nmax,ir)*hwfar(1:nmax,ir)
  enddo

  print*,'l=',l
  print*, 'Exact energy, En check, En check/norm, norm'

  do i = 1, nmax
     write(*,'(i3,1000g14.6)') i, enrg(i), sum(i), sum(i)/norm_sum(i), norm_sum(i)
  enddo

  deallocate( norm_sum, sum )

end subroutine get_ci_coef_and_energy_p16

!==============================================================
!          Define basis functions in local potential -Z0/r + potmtrx(m,n)
!==============================================================
subroutine basisl_p16 (Z0, L, alf, nmax, nload, enrg, c, h,dnorm, dnorm2, work, iwork, fac, nfac)

  use precision_type,only:p16
  implicit none
  integer::l,nmax,nload,l2,ntmp,n ,i,j,jj,jjj,lwork,matc,ierr,nfac
  real(kind=p16), dimension(1:nload,1:nmax):: c, h
  real(kind=p16), dimension(0:nfac):: fac
  real(kind=p16), dimension(1:nmax):: enrg,dnorm,dnorm2
  real(kind=p16), dimension(1:8*nmax):: work
  real(kind=p16):: z0,alf,alfl,a2,tmp,c1,c2,sm1,sm2,diag,potlz,res,fact
  integer,dimension(1:5*nmax):: iwork
  integer,dimension(1:1000):: ifail
!
!      f(n,l,r) = dnorm(n,l) * (2*alfl*r)**(l+1) * exp(-alfl*r) *
!                 * Laguerre(2*l+2;n-1;2*alfl*r)
!
! INPUT:
! -----
!  L     - orbital momentum.
!  alf   - parameter of basis.              alfl = alf
!  nmax  - number of pseudostates.          ------------------
!  nload - dimension of arrays
!  Z0    - The charge of the coulomb potential.
!
! OUTPUT:
! ------
!  enrg - eigenvalues
!  H    - Hamiltonian matrix.
!  C    - matrix of eigenvectors.
!  Work, dnorm2 - work space for the program
!
  L2 = 2 * L
  alfl = alf
  a2 = 2._p16 * alfl
!
!     Do some checks
!
  if (Z0 .lt. 0d0) then
     print*, 'Z0 is out of range,  Z0=', Z0
     stop    'Stop in BASISL'
  end if
!
!     Define factorials as n!=dexp(fac(n))
!
  fac(0) = real(0,p16)
  fact = real(0,p16)
  tmp = real(1,p16)
  ntmp = 2 * nmax + L2
  if (nfac .lt. ntmp) then
     print*, 'nfac is not enough to store array of factorials,'
     print*, 'nfac has to be more than   2 * NMAX + 2 * L =', ntmp,' but  nfac =', nfac
     stop 'nfac has to be more than   2 * NMAX + 2 * L'
  end if
  do n = 1, ntmp
     fact   = fact + log(real(n,p16))
     fac(n) = fact
  end do
!C
!C     Define normalization coeff. Dnorm
!C
  c1 = sqrt(a2)
  c2 = real(1,p16) / a2
  do i = 1, nmax
     dnorm(i)  = c1 * exp(0.5_p16* (fac(i - 1)  -  fac(L2 + 1 + i)))
     dnorm2(i) = c2 * exp(fac(L2 + i)  -  fac(i - 1))
  end do
!C
!C     Define Hamiltonian matrix
!C
  c2 = -a2 * a2 * 0.5_p16
  do i = 1, nmax
     do j = 1, i
        sm2 = real(0,p16)
        do jj = 1, j - 1
           do jjj = 1, min(i, jj)
              sm2 = sm2 + dnorm2(jjj)
           end do
        end do
        sm1 = real(0,p16)
        do jj = 1, min(i, j)
           sm1 = sm1 + dnorm2(jj)
        end do
        diag = real(0,p16)
        if (i .eq. j)  diag = 0.25_p16
        potlz = c2 * Z0 / alfl * sm1
        res = (dnorm(i) * dnorm(j) * (-real(l+j,p16) * sm1 + sm2) + diag) * c2 + dnorm(i) * dnorm(j) * ( potlz)
        h(i, j) = res
        h(j, i) = res
     end do
  end do
!
!  if matc = 0, then  only eigenvalues.
  matc = 1
  call rs_p16(nload, nmax, h, enrg, matc, c, dnorm2, work, ierr)
  lwork = 8 * nmax
!      call dsyevx('V','A','U',nmax,h,nload,0d0,0d0,0,0,0d0,nfound,
!     >   enrg,c,nload,work,lwork,iwork,ifail,info)
!      if (info.ne.0) then
!         print*,'INFO:',info
!         stop
!      endif
  if (ierr .ne. 0)  then
     print*, 'Program "RS" finished abnormaly, ierr =', ierr
     stop    'Program "RS" finished abnormaly'
  end if
return
end



subroutine rs_p16(nm,n,a,w,matz,z,fv1,fv2,ierr)

  use precision_type,only:p16
  implicit none
  real(kind=p16), dimension(1:nm,1:n):: a,z
  real(kind=p16), dimension(1:n):: w,fv1,fv2
  integer:: n,nm,ierr,matz

!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real symmetric matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real symmetric matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1  and  fv2  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
  if (n .le. nm) go to 10
  ierr = 10 * n
  go to 50
!
10 if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
  call  tred1_p16(nm,n,a,w,fv1,fv2)
!  tqlrat encounters catastrophic underflow on the Vax
!     call  tqlrat(n,w,fv2,ierr)
  call  tql1_p16(n,w,fv1,ierr)
  go to 50
!     .......... find both eigenvalues and eigenvectors ..........
20 call  tred2_p16(nm,n,a,w,fv1,z)
  call  tql2_p16(nm,n,w,fv1,z,ierr)
50 return
end

subroutine tql1_p16(n,d,e,ierr)

  use precision_type,only:p16
  implicit none
  real(kind=p16), dimension(1:n):: d,e
  integer:: i,j,l,m,n,ii,l1,l2,mml,ierr
  real(kind=p16):: c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag_p16,one
!
!     this subroutine is a translation of the algol procedure tql1,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues of a symmetric
!     tridiagonal matrix by the ql method.
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        e has been destroyed.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
  one=real(1,p16)
  ierr = 0
  if (n .eq. 1) go to 1001
!
  do 100 i = 2, n
  100 e(i-1) = e(i)
!
  f = real(0,p16)
  tst1 = real(0,p16)
  e(n) = real(0,p16)
!
  do 290 l = 1, n
     j = 0
     h = abs(d(l)) + abs(e(l))
     if (tst1 .lt. h) tst1 = h
!     .......... look for small sub-diagonal element ..........
     do 110 m = l, n
        tst2 = tst1 + abs(e(m))
        if (tst2 .eq. tst1) go to 120
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
110     continue
!
120    if (m .eq. l) go to 210
130    if (j .eq. 30) go to 1000
       j = j + 1
!     .......... form shift ..........
       l1 = l + 1
       l2 = l1 + 1
       g = d(l)
       p = (d(l1) - g) / (real(2,p16) * e(l))
       r = pythag_p16(p,one)
       d(l) = e(l) / (p + sign(r,p))
       d(l1) = e(l) * (p + sign(r,p))
       dl1 = d(l1)
       h = g - d(l)
       if (l2 .gt. n) go to 145
!
       do 140 i = l2, n
140    d(i) = d(i) - h
!
145    f = f + h
!     .......... ql transformation ..........
       p = d(m)
       c = real(1,p16)
       c2 = c
       el1 = e(l1)
       s = real(0,p16)
       mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
      do 200 ii = 1, mml
         c3 = c2
         c2 = c
         s2 = s
         i = m - ii
         g = c * e(i)
         h = c * p
         r = pythag_p16(p,e(i))
         e(i+1) = s * r
         s = e(i) / r
         c = p / r
         p = c * d(i) - s * g
         d(i+1) = h + s * (c * g + s * d(i))
200    continue
!
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + abs(e(l))
         if (tst2 .gt. tst1) go to 130
  210    p = d(l) + f
!     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
!     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
!
  250    i = 1
  270    d(i) = p
  290 continue
!
  go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
end
subroutine tql2_p16(nm,n,d,e,z,ierr)
!

  use precision_type,only:p16
  implicit none
  integer:: i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
  real(kind=p16),dimension(1:n):: d,e
  real(kind=p16),dimension(nm,n):: z
  real(kind=p16):: c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag_p16,one
!
!     this subroutine is a translation of the algol procedure tql2,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      one=real(1,p16)
      ierr = 0
      if (n .eq. 1) go to 1001
!
      do 100 i = 2, n
  100 e(i-1) = e(i)
!
      f = real(0,p16)
      tst1 = real(0,p16)
      e(n) = real(0,p16)
!
      do 240 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (tst1 .lt. h) tst1 = h
!     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + abs(e(m))
            if (tst2 .eq. tst1) go to 120
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
  110    continue
!
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (real(2,p16) * e(l))
         r = pythag_p16(p,one)
         d(l) = e(l) / (p + sign(r,p))
         d(l1) = e(l) * (p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
!
         do 140 i = l2, n
  140    d(i) = d(i) - h
!
  145    f = f + h
!     .......... ql transformation ..........
         p = d(m)
         c = real(1,p16)
         c2 = c
         el1 = e(l1)
         s = real(0,p16)
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag_p16(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
!     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
!
  200    continue
!
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + abs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
!     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
!
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
!
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
!
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
!
  300 continue
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
end

subroutine tred1_p16(nm,n,a,d,e,e2)

  use precision_type,only:p16
  implicit none
  integer:: i,j,k,l,n,ii,nm,jp1
  real(kind=p16),dimension(1:n):: d,e,e2
  real(kind=p16),dimension(nm,n):: a
  real(kind=p16):: f,g,h,scale
!
!     this subroutine is a translation of the algol procedure tred1,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix
!     to a symmetric tridiagonal matrix using
!     orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction in its strict lower
!          triangle.  the full upper triangle of a is unaltered.
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
!     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = real(0,p16)
         scale = real(0,p16)
         if (l .lt. 1) go to 130
!     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + abs(d(k))
!
         if (scale .ne. real(0,p16)) go to 140
!
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = real(0,p16)
  125    continue
!
  130    e(i) = real(0,p16)
         e2(i) = real(0,p16)
         go to 300
!
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
!
         e2(i) = scale * scale * h
         f = d(l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
!     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = real(0,p16)
!
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
!
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       continue
!
  220       e(j) = g
  240    continue
!     .......... form p ..........
         f = real(0,p16)
!
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
!
         h = f / (h + h)
!     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
!     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
!
            do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
!
  280    continue
!
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
!
  300 continue
!
return
end

subroutine tred2_p16(nm,n,a,d,e,z)
!

  use precision_type,only:p16
  implicit none
  integer:: i,j,k,l,n,ii,nm,jp1
  real(kind=p16),dimension(1:n):: d,e
  real(kind=p16),dimension(nm,n):: a,z
  real(kind=p16):: f,g,h,hh,scale
!
!     this subroutine is a translation of the algol procedure tred2,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix to a
!     symmetric tridiagonal matrix using and accumulating
!     orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        z contains the orthogonal transformation matrix
!          produced in the reduction.
!
!        a and z may coincide.  if distinct, a is unaltered.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      do 100 i = 1, n
!
         do 80 j = i, n
   80    z(j,i) = a(j,i)
!
         d(i) = a(n,i)
  100 continue
!
      if (n .eq. 1) go to 510
!     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = real(0,p16)
         scale = real(0,p16)
         if (l .lt. 2) go to 130
!     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + abs(d(k))
!
         if (scale .ne. real(0,p16)) go to 140
  130    e(i) = d(l)
!
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = real(0,p16)
            z(j,i) = real(0,p16)
  135    continue
!
         go to 290
!
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
!
         f = d(l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
!     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = real(0,p16)
!
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
!
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
!
  220       e(j) = g
  240    continue
!     .......... form p ..........
         f = real(0,p16)
!
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
!
         hh = f / (h + h)
!     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
!     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
!
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
!
            d(j) = z(l,j)
            z(i,j) = real(0,p16)
  280    continue
!
  290    d(i) = h
  300 continue
!     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = real(1,p16)
         h = d(i)
         if (h .eq. real(0,p16)) go to 380
!
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
!
         do 360 j = 1, l
            g = real(0,p16)
!
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
!
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
!
  380    do 400 k = 1, l
  400    z(k,i) = real(0,p16)
!
  500 continue
!
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = real(0,p16)
  520 continue
!
      z(n,n) = real(1,p16)
      e(1) = real(0,p16)
return
end

function pythag_p16(a,b)

  use precision_type,only:p16
  implicit none
  real(kind=p16)::a,b,p,r,s,t,u,pythag_p16
!
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!
  p = max(abs(a),abs(b))
  if (p .eq. real(0,p16)) go to 20
  r = (min(abs(a),abs(b))/p)**2
10 continue
  t = real(4,p16) + r
  if (t .eq. real(4,p16)) go to 20
  s = r/t
  u = real(1,p16) + real(2,p16)*s
  p = u*p
  r = (s/u)**2 * r
  go to 10
20 pythag_p16 = p

return
end

subroutine comp_wf_b_ar_p16(l,r,nmax,c,alf,wfid)

  use precision_type,only:p16,id
  implicit none
  integer :: n,l,k,nmax
  real(kind=p16),dimension(0:nmax):: basis
  real(kind=p16),dimension(1:nmax,1:nmax):: c
  real(kind=p16):: r,alf
  real(kind=p16),dimension(1:nmax)::wf
  real(kind=id),dimension(1:nmax)::wfid

  call build_laguerre6_p16(l,r*alf,nmax,alf,basis(0:nmax))

  do n=1,nmax
    wf(n)=0._p16
    do k=1,nmax
      wf(n)=wf(n)+c(k,n)*basis(k-1)
    enddo ! k
    wfid(n)=real(wf(n),id)
  enddo

return
end subroutine comp_wf_b_ar_p16



subroutine comp_wf_b_p16(n,l,r,nmax,c,alf,wf)

  use precision_type,only:p16
  implicit none
  integer :: n,l,k,nmax
  real(kind=p16),dimension(0:nmax):: basis
  real(kind=p16),dimension(1:nmax,1:nmax):: c
  real(kind=p16):: wf,r,alf

  call build_laguerre6_p16(l,r*alf,nmax,alf,basis(0:nmax))

  wf=0._p16
  do k=1,nmax
    wf=wf+c(k,n)*basis(k-1)
  enddo ! k

return
end

subroutine build_laguerre6_p16(l,x,n,alf2,laguerre)

  use precision_type,only:p16
  use data_base
  implicit none
  integer :: l,i,n
  real(kind=p16):: beta,lag,alf2,x,fact2
  real(kind=p16),dimension(0:n):: laguerre
!!!! Laguerre/r^l
!!!! for angular z integration

  beta=real(2*l+2,p16)
!!fact2(:)=2*alf2*exp(l*log(2*x(:))-x(:))
!fact2(:)=(2.0_dp**(l+1))*(alf2**(l+1))*exp(-x(:)) !*x(:)
!  fact2=((2*alf2)**(l+1))*exp(-x) !*x(:)
  fact2=((2*x)**(l+1))*alf2/x*exp(-x) !*x(:)
  laguerre(0)=1._p16
  if(n>0)then
    laguerre(1)=(1+beta-2*x)
    do i=1,n-1
      laguerre(i+1)=((beta+2*i+1-2*x)*laguerre(i)-(beta+i)*laguerre(i-1))/real(i+1,p16)
    enddo
  endif

  do i=0,n
    laguerre(i)=laguerre(i)*fact2
  enddo

return
end

subroutine comp_dwf_ar_p16(l,r,nmax,c,alf,z0,wff,wfid)

  use precision_type,only:p16,id
  implicit none
  integer :: n,l,k,nmax
  real(kind=p16),dimension(0:nmax):: basis,bas
  real(kind=p16),dimension(1:nmax,1:nmax):: c
  real(kind=p16):: r,alf,tmp,z0
  real(kind=p16),dimension(1:nmax)::wf,wff
  real(kind=id),dimension(1:nmax)::wfid

  call build_laguerre6_p16(l,r*alf,nmax,alf,bas(0:nmax))
  call build_laguerred_p16(l,r*alf,nmax,alf,basis(0:nmax))

  do n=1,nmax
    wf(n)=0._p16
    !tmp=-0.5_p16*(2._p16+alf*alf*r-2._p16*alf*(1+l))/r*bas(0)

    tmp=-0.5_p16*(2._p16*z0+alf*alf*r-2._p16*alf*(1+l))/r*bas(0)

    wf(n)=wf(n)+c(1,n)*tmp
    do k=2,nmax
      !tmp=-0.5_p16*(basis(k-2)+(2._p16+alf*alf*r-2._p16*alf*(k+l))/r*bas(k-1))

      tmp=-0.5_p16*(basis(k-2)+(2._p16*z0+alf*alf*r-2._p16*alf*(k+l))/r*bas(k-1))

      wf(n)=wf(n)+c(k,n)*tmp
    enddo ! k
    wfid(n)=real(wf(n),id)
  enddo

return
end subroutine comp_dwf_ar_p16

subroutine build_laguerred_p16(l,x,n,alf2,laguerre)

  use precision_type,only:p16
  use data_base
  implicit none
  integer :: l,i,n
  real(kind=p16):: beta,lag,alf2,x,fact2
  real(kind=p16),dimension(0:n):: laguerre
!!!! Laguerre/r^l
!!!! for angular z integration

  beta=real(2*l+3,p16)
!!fact2(:)=2*alf2*exp(l*log(2*x(:))-x(:))
!fact2(:)=(2.0_dp**(l+1))*(alf2**(l+1))*exp(-x(:)) !*x(:)
!  fact2=((2*alf2)**(l+1))*exp(-x) !*x(:)
  fact2=((2*x)**(l+2))*(alf2/x)**3*exp(-x) !*x(:)
  laguerre(0)=1._p16
  if(n>0)then
    laguerre(1)=(1+beta-2*x)
    do i=1,n-1
      laguerre(i+1)=((beta+2*i+1-2*x)*laguerre(i)-(beta+i)*laguerre(i-1))/real(i+1,p16)
    enddo
  endif

  do i=0,n
    laguerre(i)=laguerre(i)*fact2
  enddo

return
end subroutine build_laguerred_p16


subroutine comp_wf_b(n,l,r,wf)

  use precision_type
  use data_base,only:c_t,alf_t,nps_t
  implicit none
  integer :: n,l,k
  real(kind=id),dimension(0:nps_t(l)):: basis
  real(kind=id):: wf,r

  call build_laguerre6(l,r*real(alf_t(l),id),nps_t(l),real(alf_t(l),id),basis)

  wf=0._id
  do k=1,nps_t(l)
    wf=wf+c_t(l)%cl_t(k,n)*basis(k-1)
  enddo ! k

return
end

subroutine build_laguerre6(l,x,n,alf2,laguerre)

  use precision_type
  use data_base
  implicit none
  integer :: l,i,n
  real(kind=id):: beta,lag,alf2,x,fact2
  real(kind=id),dimension(0:n):: laguerre
!!!! Laguerre/r^l
!!!! for angular z integration

  beta=real(2*l+2,id)
!!fact2(:)=2*alf2*exp(l*log(2*x(:))-x(:))
!fact2(:)=(2.0_dp**(l+1))*(alf2**(l+1))*exp(-x(:)) !*x(:)
!  fact2=((2*alf2)**(l+1))*exp(-x) !*x(:)
  fact2=((2*x)**(l+1))*alf2/x*exp(-x) !*x(:)
  laguerre(0)=real(1,id)
  if(n>0)then
    laguerre(1)=(1+beta-2*x)
    do i=1,n-1
      laguerre(i+1)=((beta+2*i+1-2*x)*laguerre(i)-(beta+i)*laguerre(i-1))/(i+1)
    enddo
  endif

  do i=0,n
    laguerre(i)=laguerre(i)*fact2
  enddo

return
end


