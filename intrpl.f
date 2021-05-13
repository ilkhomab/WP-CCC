C     FILE 'INTRPL.FTN77'
C
      SUBROUTINE INTRPL(L,X,Y,N,U,V)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      real*4 FINT
      integer*8 ideriv
      PARAMETER ( NR=1000 )
C
C     DOUBLE PRECISION INTERPOLATION OF A SINGLE VALUED FUNCTION
C     THIS SUBROUTINE INTERPOLATES, FROM VALUES OF THE FUNCTION
C     GIVEN  AS ORDINATES OF INPUT DATA POINTS IN AN X-Y PLANE
C     AND FOR A GIVEN SET OF X VALUES(ABSCISSAE),THE VALUES OF
C     A SINGLE VALUED FUNCTION Y=Y(X).
C
C     THE SUBROUTINE ALSO CALCULATES FIRST DERIVATIVES DV(X)
C     AND THE INTEGRAL OVER THE WHOLE RANGE OF X FINT.
C     DV(X) IS IN COMMON/DV/, CALCULATED IF IDERIV.NE.0
C     FINT IS IN COMMON/INTEG/, CALCULATED IF INT.NE.0
C     QQ(4,L) ARE THE SPLINE COEFFICIENTS. THEY ARE STORED IN COMMON/QQ/
C
C     THE INPUT PARAMETERS ARE;
C
C     L=NUMBER OF DATA POINTS
C     (MUST BE TWO OR GREATER)
C     X=ARRAY OF DIMENSION L STORING THE X VALUES
C     OF INPUT DATA POINTS (IN ASCENDING ORDER)
C     Y=ARRAY OF DIMENSION L STORING THE Y VALUES OF INPUT DATA POINTS
C     N=NUMBER OF POINTS AT WHICH INTERPOLATION OF THE Y-VALUES
C     IS REQUIRED (MUST BE 1 OR GREATER)
C     U=ARRAY OF DIMENSION N STORING THE X VALUES
C     OF THE DESIRED POINTS
C
C     THE OUTPUT PARAMETER IS V=ARRAY OF DIMENSION N WHERE THE
C     INTERPOLATED Y VALUES ARE TO BE DISPLAYED
C
      COMMON/DV/ IDERIV,DV(NR)
      COMMON/INTEG/ INT,FINT
      COMMON/QQ/ QQ(4,NR)
      DIMENSION X(L),Y(L),U(N),V(N)
      EQUIVALENCE (P0,X3),(Q0,Y3),(Q1,T3)
      REAL*8 M1,M2,M3,M4,M5
      EQUIVALENCE (UK,DX),(IMN,X2,A1,M1),(IMX,X5,A5,M5),
     1         (J,SW,SA),(Y2,W2,W4,Q2),(Y5,W3,Q3)
C
C     PRELIMINARY PROCESSING
C
    
    
      L0 = L
      LM1 = L0 - 1
      LM2 = LM1 - 1
      LP1 = L0 + 1
      N0 = N
      NQQ = NR
        IF(N0.GT.NQQ) THEN
        NQQV = NQQ
        WRITE(6,2089) NQQV,N0
        STOP
        END IF
      IF(LM2.LT.0) GO TO 90
      IF(N0.LE.0) GO TO 91
      DO 11 I = 2,L0
      IF(X(I-1)-X(I))11,95,96
   11 CONTINUE
      IPV = 0
C
C     MAIN LOOP
C
      FINT = 0.0D0
      DO 80 K = 1,N0
      UK = U(K)
C
C     BINARY SEARCH TO LOCATE THE DESIRED POINT
C
      IF(LM2.EQ.0) GO TO 27
      IF(UK.GE.X(L0)) GO TO 26
      IF(UK.LT.X(1)) GO TO 25
      IMN = 2
      IMX = L0
   21 I = (IMN+IMX)/2
      IF(UK.GE.X(I)) GO TO 23
      IMX = I
      GO TO 24
   23 IMN = I + 1
   24 IF(IMX.GT.IMN) GO TO 21
      I = IMX
      GO TO 30
   25 I = 1
      GO TO 30
   26 I = LP1
      GO TO 30
   27 I = 2
C
C     CHECK IF I = IPV
C
   30 IF(I.EQ.IPV) GO TO 70
      IPV = I
C
C     ROUTINES TO PICK UP NECESSARY X AND Y VALUES AND TO
C     ESTIMATE THEM IF NECESSARY
C
      J = I
      IF(J.EQ.1) J = 2
      IF(J.EQ.LP1) J = L0
      X3 = X(J-1)
      Y3 = Y(J-1)
      X4 = X(J)
      Y4 = Y(J)
      A3 = X4 - X3
      M3 = (Y4 - Y3)/A3
      IF(LM2.EQ.0) GO TO 43
      IF(J.EQ.2) GO TO 41
      X2 = X(J-2)
      Y2 = Y(J-2)
      A2 = X3 - X2
      M2 = (Y3 - Y2)/A2
      IF(J.EQ.L0) GO TO 42
   41 X5 = X(J+1)
      Y5 = Y(J+1)
      A4 = X5 - X4
      M4 = (Y5 - Y4)/A4
      IF(J.EQ.2) M2 = M3 + M3 - M4
      GO TO 45
   42 M4 = M3 + M3 - M2
      GO TO 45
   43 M2 = M3
   45 IF(J.LE.3) GO TO 46
      A1 = X2 - X(J-3)
      M1 = (Y2 - Y(J-3))/A1
      GO TO 47
   46 M1 = M2 + M2 - M3
   47 IF(J.GE.LM1) GO TO 48
      A5 = X(J+2) - X5
      M5 = (Y(J+2)-Y5)/A5
      GO TO 50
   48 M5 = M4 + M4 - M3
C
C     NUMERICAL DIFFERENTIATION
C
   50 IF(I.EQ.LP1) GO TO 52
      W2 = ABS(M4-M3)
      W3 = ABS(M2-M1)
      SW = W2 + W3
      IF(SW.NE.0.0D0) GO TO 51
      W2 = 0.5D0
      W3 = 0.5D0
      SW = 1.0D0
   51 T3 = (W2*M2 + W3*M3)/SW
      IF(I.EQ.1) GO TO 54
   52 W3 = ABS(M5-M4)
      W4 = ABS(M3-M2)
      SW = W3 + W4
      IF(SW.NE.0.0D0) GO TO 53
      W3 = 0.5D0
      W4 = 0.5D0
      SW = 1.0D0
   53 T4 = (W3*M3 + W4*M4)/SW
      IF(I.NE.LP1) GO TO 60
      T3 = T4
      SA = A2 + A3
      T4 = 0.5D0*(M4 + M5 - A2*(A2-A3)*(M2-M3)/(SA*SA))
      X3 = X4
      Y3 = Y4
      A3 = A2
      M3 = M4
      GO TO 60
   54 T4 = T3
      SA = A3 + A4
      T3 = 0.5D0*(M1 + M2 - A4*(A3-A4)*(M3-M4)/(SA*SA))
      X3 = X3 - A4
      Y3 = Y3 - M2*A4
      A3 = A4
      M3 = M2
C
C     DETERMINATION OF THE COEFFICIENTS
C
   60 Q2 = (2.0D0*(M3-T3) + M3 - T4)/A3
      Q3 = (-M3 - M3 + T3 + T4)/(A3*A3)
C
C     COMPUTATION OF THE POLYNOMIAL
C
   70 DX = UK - P0
      V(K) = Q0 + DX*(Q1 + DX*(Q2 + DX*Q3))
      IF(IDERIV.EQ.0) GO TO 72
      DV(K) = Q1 + DX*(2.0D0*Q2 + DX*3.0D0*Q3)
      QQ(1,K) = Q0
      QQ(2,K) = Q1
      QQ(3,K) = Q2
      QQ(4,K) = Q3
   72 IF(INT.EQ.0) GO TO 80
      IF(K.LT.N0) DD = U(K+1) - U(K)
      IF(K.LT.N0) FINT = FINT
     1 + DD*(Q0 + DD*(Q1/2.0D0+DD*(Q2/3.0D0+DD*Q3/4.0D0)))
   80 CONTINUE
      RETURN
C
C     ERROR EXITS
C
   90 WRITE(6,2090)
      GO TO 99
   91 WRITE(6,2091)
      GO TO 99
   95 WRITE(6,2095)
      GO TO 97
   96 WRITE(6,2096)
   97 WRITE(6,2097)I,X(I)
   99 WRITE(6,2099) L0,N0
      RETURN
C
C     FORMAT STATEMENTS
C
 2089 FORMAT( 'WARNING ERROR IN INTRPL. MAX ALLOWED VALUE OF N0 IS',
     *I3,' HERE N0 IS',I3)
 2090 FORMAT(1X/' N  =  1 OR LESS.'/)
 2091 FORMAT(1X/' N  =  0 OR LESS.'/)
 2095 FORMAT(1X/' IDENTICAL X VALUES.'/)
 2096 FORMAT(1X/' X VALUES OUT OF SEQUENCE.'/)
 2097 FORMAT(6H   I =,I7,10X,6HX(I) =,E12.3)
 2099 FORMAT(6H   L =,I7,10X,3HN =,I7/
     1 ' ERRROR DETECTED IN ROUTINE INTRPL')
      END
      BLOCK DATA
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      real*4 FINT
      integer*8 ideriv
      PARAMETER ( NR=1000 )
      COMMON/DV/ IDERIV,DV(NR)
      COMMON/INTEG/ INT,FINT
      DATA IDERIV,INT/1,0/
      END
      subroutine intrpl1(n,x,y,nn,u,v)
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 x(n), y(n), b(n), c(n), d(n), u(nn), v(nn), seval, vold(nn)
      ii = 1

c$$$      call intrplold(n,x,y,nn,u,vold)
      call spline(n,x,y,b,c,d)
      do i = 1, nn
         v(i) = seval(ii, n, u(i), x, y, b, c, d)
c$$$         if (mod(i,100).eq.0) print'(i5,1p,5e12.4)',
c$$$     >      i,x(1),u(i),x(n),v(i),vold(i)
      enddo
      return
      
      end
      
      subroutine spline (n, x, y, b, c, d)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer n
      real*8 x(n), y(n), b(n), c(n), d(n)
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
      integer nm1, ib, i
      double precision t
c
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.d0*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.d0
      c(n) = 0.d0
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.d0*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.d0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.d0*c(i)
   40 continue
      c(n) = 3.d0*c(n)
      d(n) = d(n-1)
      return
c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.d0
      d(1) = 0.d0
      b(2) = b(1)
      c(2) = 0.d0
      d(2) = 0.d0
      return
      end

      double precision function seval(i, n, u, x, y, b, c, d)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer n
      double precision  u, x(n), y(n), b(n), c(n), d(n)
c
c  this subroutine evaluates the cubic spline function
c
c    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
c
c    where  x(i) .lt. u .lt. x(i+1), using horner's rule
c
c  if  u .lt. x(1) then  i = 1  is used.
c  if  u .ge. x(n) then  i = n  is used.
c
c  input..
c
c    n = the number of data points
c    u = the abscissa at which the spline is to be evaluated
c    x,y = the arrays of data abscissas and ordinates
c    b,c,d = arrays of spline coefficients computed by spline
c
c  if  u  is not in the same interval as the previous call, then a
c  binary search is performed to determine the proper interval.
c
      integer i, j, k
      double precision dx
c$$$      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
c
c  binary search
c
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
c
c  evaluate spline
c
   30 dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
      end      
 
