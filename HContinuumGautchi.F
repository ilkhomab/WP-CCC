
!  f90_kind.f90
!
!  Module defining the KIND numbers for the NAGWare f90 Compiler.
!
      MODULE f90_kind
!
      INTRINSIC kind,selected_int_kind,selected_real_kind  ! we use these,
!
! Indicator that the KIND= is not available for this compiler/host
      INTEGER, PARAMETER :: not_available = -1  
!
! Real and Complex numbers
!   Single precision
      INTEGER, PARAMETER :: single  = kind(0.0)
!   Double precision
      INTEGER, PARAMETER :: double  = kind(0.0d0)
!   Quadruple precision
      INTEGER, PARAMETER :: quad    = selected_real_kind(p=30)
!
! INTEGERs numbers:

!   Single byte INTEGER
      INTEGER, PARAMETER :: int8    = selected_int_kind(2)
!   Two byte INTEGER   
      INTEGER, PARAMETER :: int16   = selected_int_kind(4)
!   Four byte INTEGER
      INTEGER, PARAMETER :: int32   = selected_int_kind(9)
!   Eight byte INTEGER
      INTEGER, PARAMETER :: int64   = selected_int_kind(18)
!
! Logical values
!   Single byte logical
      INTEGER, PARAMETER :: byte    = 1
!   Two byte logical
      INTEGER, PARAMETER :: twobyte = 2
!   Four byte logical
      INTEGER, PARAMETER :: word    = kind(.TRUE.)
!
! Character type
!   Normal single byte character (ASCII sequence)
      INTEGER, PARAMETER :: ascii   = kind('x')

      INTEGER, PARAMETER :: UB = 40000

      real(double),parameter::pi=3.14159265358979323846264338328d0
      real(double),parameter::pifct=0.797884560802865355879892119869d0
      
      END module f90_kind

         
  
!      PROGRAM TEST
!      
!      USE f90_kind
!      
!      IMPLICIT NONE 
!      
! 
!! The largest angular momentum needed:
!           
!      INTEGER, PARAMETER :: Lmax = 1000 
!
!! The momentum kmom   where radial wave function needed:
!      
!      REAL(double)  :: kmom 
!      
!! The nuclear charge of the Hydrogenic atom:
!      
!      REAL(double)  :: ZE
!   
!! The returned radial wave function at rpt for all L's
!! and the Coulomb phase 
!     
!      REAL(double) ::  Radial(0:Lmax), DeltaL(0:Lmax)
!
!! how many r points you need:
!      
!      INTEGER, PARAMETER :: Nrpt = 1250
!
!! We assume equal spacing, but you can have any points, 
!! the caculations are totally indepdent for each r points:
!
!      REAL(double), PARAMETER :: deltaR = 0.2D0 
!      
!      INTEGER :: i, iL, ir, ik, Lmax0
!      
!      REAL(double)  :: rpt(1:Nrpt), radial_wave(1:Nrpt, 0:Lmax), r
!      
!      REAL(double) :: mom(1:4), rp(1:5)  
!      
!      CHARACTER*3 :: CiL
!      
!      
!
!
!! The nuclear charge, we set it for H atom:
!      
!      ZE = 1.D0 
!      
!     
!! Read in the iL, momentum and radial points need to be calculated: 
!      
!      OPEN(20,file='test_input.dat') 
!      
!      READ(20, *) iL 
!      READ(20, *) mom
!      READ(20, *) rp 
!      
!      IF(iL>Lmax) THEN       
!         WRITE(*, *) 'Please increase Lmax in main program to ', iL 
!	 WRITE(*, *) ' and try again.'  
!         STOP       
!      ENDIF
!      
!       
! 
!      OPEN(10, file='test_output_table_1.dat') 
!     
!      WRITE(10, *)       
!      WRITE(10, *) ' L =', iL, '     Results: ' 
!      WRITE(10, *) 
!      WRITE(10, *) 
!    
!      DO ik = 1, 4 
!      
!         kmom = mom(ik) 
!
!         WRITE(10, 100)   kmom   
! 
!         DO ir = 1, 5 
!      
!            r = rp(ir) 
!	  
!            CALL Hydrogenic_Continuum_Wave(kmom, r, ZE, iL,
!     &                               Radial)     
!
!	    WRITE(10, 200)   rp(ir),  Radial(iL) 
!        
!         ENDDO 
!       
!         WRITE(10,*) 
!       
!      ENDDO 
!       
!       CLOSE(10)
! 
! 100   FORMAT('k = ', F5.2, ' :' )      
! 200   FORMAT('  r = ',  F10.2, 2X, ', RadialWave =',  E24.16) 
! 
!      
!       WRITE(*, *) '=================================================='     
!       WRITE(*, *)
!       WRITE(*, *) '(1) Finished reproducing results in Table 1.'
!       WRITE(*, *) '    Please check the file test_output_table1.dat.'
!       
!       WRITE(*, *)
! 
! 
!      READ(20, *) iL 
!      READ(20, *) mom
!      READ(20, *) rp 
!      
!      IF(iL>Lmax) THEN       
!         WRITE(*, *) 'Please increase Lmax in main program to ', iL 
!	 WRITE(*, *) ' and try again.'  
!         STOP       
!      ENDIF
!      
!       
! 
!      OPEN(10, file='test_output_table_2.dat') 
!     
!      WRITE(10, *)       
!      WRITE(10, *) ' L =', iL, '     Results: ' 
!      WRITE(10, *) 
!      WRITE(10, *) 
!    
!      DO ik = 1, 4 
!      
!         kmom = mom(ik) 
!
!         WRITE(10, 100)   kmom   
! 
!         DO ir = 1, 5 
!      
!            r = rp(ir) 
!	  
!            CALL Hydrogenic_Continuum_Wave(kmom, r, ZE, iL,
!     &                               Radial)     
!
!	    WRITE(10, 200)   rp(ir),  Radial(iL) 
!        
!         ENDDO 
!       
!         WRITE(10,*) 
!       
!      ENDDO 
!       
!       CLOSE(10)
! 
!  
!       WRITE(*, *)
!  
!       WRITE(*, *) '(2) Finished reproducing results in Table 2.'
!       WRITE(*, *) '    Please check the file test_output_table_2.dat.'
!       
!       WRITE(*, *)
!
!
!      READ(20, *) iL 
!      READ(20, *) mom
!      READ(20, *) rp 
!      
!      IF(iL>Lmax) THEN       
!         WRITE(*, *) 'Please increase Lmax in main program to ', iL 
!	 WRITE(*, *) ' and try again.'  
!         STOP       
!      ENDIF
!      
!       
! 
!      OPEN(10, file='test_output_table_3.dat') 
!     
!      WRITE(10, *)       
!      WRITE(10, *) ' L =', iL, '     Results: ' 
!      WRITE(10, *) 
!      WRITE(10, *) 
!    
!      DO ik = 1, 4 
!      
!         kmom = mom(ik) 
!
!         WRITE(10, 100)   kmom   
! 
!         DO ir = 1, 5 
!      
!            r = rp(ir) 
!	  
!            CALL Hydrogenic_Continuum_Wave(kmom, r, ZE, iL,
!     &                               Radial)     
!
!	    WRITE(10, 200)   rp(ir),  Radial(iL) 
!        
!         ENDDO 
!       
!         WRITE(10,*) 
!       
!      ENDDO 
!       
!       CLOSE(10)
! 
!   
!       WRITE(*, *) '(3) Finished reproducing results in Table 3.'
!       WRITE(*, *) '    Please check the file test_output_table_3.dat.'
!       
!       WRITE(*, *)
!
!
!
!      READ(20, *) iL 
!      READ(20, *) mom
!      READ(20, *) rp 
!      
!      IF(iL>Lmax) THEN       
!         WRITE(*, *) 'Please increase Lmax in main program to ', iL 
!	 WRITE(*, *) ' and try again.'  
!         STOP       
!      ENDIF
!      
!      CLOSE(20) 
!      
! 
!      OPEN(10, file='test_output_table_4.dat') 
!     
!      WRITE(10, *)       
!      WRITE(10, *) ' L =', iL, '     Results: ' 
!      WRITE(10, *) 
!      WRITE(10, *) 
!    
!      DO ik = 1, 4 
!      
!         kmom = mom(ik) 
!
!         WRITE(10, 100)   kmom   
! 
!         DO ir = 1, 5 
!      
!            r = rp(ir) 
!	  
!            CALL Hydrogenic_Continuum_Wave(kmom, r, ZE, iL,
!     &                               Radial)     
!
!	    WRITE(10, 200)   rp(ir),  Radial(iL) 
!        
!         ENDDO 
!       
!         WRITE(10,*) 
!       
!      ENDDO 
!       
!       CLOSE(10)
! 
!  
!       WRITE(*, *)
!  
!       WRITE(*, *) '(4) Finished reproducing results in Table 4.'
!       WRITE(*, *) '    Please check the file test_output_table_4.dat.'
!       
!       WRITE(*, *)
!
!
! 
!       kmom = 1.d0 
!       
!       WRITE(*, *) '=================================================='     
!       WRITE(*, *) '(5) Compute the H continuum wave  for k=1.0,'
!
!       WRITE(*, *) '      for L = 0 to Lmax= 200, and r = 0.2 to 250. '
!       
!       WRITE(*, *) ' \phi_kl for each l written into a different file.'                           
!       WRITE(*, *) '=================================================='     
!       WRITE(*, *)
!
!
!
!! For each r point, call the subrouinte :
!
!      Lmax0 = 200
!	
!      DO i =  1, Nrpt
!      
!         rpt(i) =  dble(i) * deltaR
!     
!         CALL Hydrogenic_Continuum_Wave(kmom, rpt(i), ZE, Lmax0,
!     &                               Radial)     
!                 
!	   radial_wave(i, 0:Lmax0) = Radial(0:Lmax0)
!        
!      ENDDO    
!      
!
!! Output the  radial wave functions:
!      
!      DO iL = 0, Lmax0
!      
!         WRITE(CiL, '(I3.3)') iL 
!	 
!	 OPEN(19, file='radial_wave_'//CiL//'.dat')
!	 
!	 DO i = 1, Nrpt
!	 
!	   WRITE(19, '( F7.2, 1X, E23.15)') rpt(i),  radial_wave(i, iL)
!    
!	 ENDDO 
!	 
!	 CLOSE(19) 
!	 
!      ENDDO
!
!       WRITE(*, *)
!       WRITE(*, *) ' Finished! Check 200 wave function files written.'
!       WRITE(*, *)
!
!
!       WRITE(*, *) '=================================================='     
!       WRITE(*, *) '(6) Coulomb phase is calculated but w/t ouput.'
!       WRITE(*, *) '=================================================='     
!       WRITE(*, *)
!
!       
!! Then we calculate the Coullomb phase 
!! for all L's from L = 0, ..., Lmax : 
!        
!      CALL  sigma_arb(kmom, Lmax0, ZE, DeltaL)
!
! 
!      WRITE(*, *) ' All test runs are finished!'
!           
!      WRITE(*, *)
!      
!      END  PROGRAM TEST
      

!=============================================================
!
! SUBROUTINE Hydrogenic_Continuum_Wave(kmom, rpt, ZE, Lmax,
!     &                                      Radial, DeltaL)! 
! 
! computes the radial wave function  Radial and the Coulomb
! phase sigma for given momentum kmom, grid point rpt,
! nuclear charge ZE, and the largest angurlar momentum Lmax,
! for all the L's: 0, 1, 2, 3, ...., Lmax 
! at a single call
!
! Authors: Liang-You Peng, Peking University,
! Dates: June, 2010 
!
!=============================================================            

      SUBROUTINE Hydrogenic_Continuum_Wave(kmom, rpt, ZE, Lmax,Radial)
      USE f90_kind
      
      IMPLICIT NONE 

! The momentum kmom and grid point rpt, where radial wave function needed:
      
      REAL(double)  :: kmom, rpt
      
! The nuclear charge of the Hydrogenic atom:
      
      REAL(double), intent(in)  :: ZE

! The largest angular momentum needed:
      
      INTEGER, intent(in) :: Lmax
  
! The returned radial wave function at rpt for all L's
! and the Coulomb phase 
     
      REAL(double) ::  Radial(0:Lmax)
      
!--------------------------------------------------------      
! Local variables:

      REAL(double) :: rho, eta 
      
      
      rho = rpt * kmom
      
      eta = - ZE / kmom 

      
! First calculate the radial wave function for all L's
! from L = 0, ..., Lmax :
      
       CALL Regular_Coulomb_Fun(eta, rho, Lmax, Radial)

! Due to the different convention, there is a factor of 2,
! which should be multipled in order to be consistent with 
! Landau's momentum normalized continuum states:

!       Radial  =  2.D0 * Radial   
!       Radial  =  pifct * Radial   

 
      END  SUBROUTINE Hydrogenic_Continuum_Wave   
      
      

      
      
!===============================================================
! This code is very accurate, the results can be 
! checked against those from Mathematica by using 
!
! 2^L Exp[ - Pi yita / 2.0] Abs[Gamma[L + 1.0 + I yita ]] / 
!  Gamma[ 2 L + 2.0] rho^(L + 1) Exp[-I rho]
!  Hypergeometric1F1[L + 1 - I yita, 2.0 * L + 2.0, 2 I rho]
!
!===============================================================
!   SUBROUTINE Regular_Coulomb_Fun(eta, ro, Lmax,  F) 
! returns an array F[0:Lmax] for the values of 
! regular Coulomb wave function at given arguments
! eta and ro. Accuracy: around 14 digits.
! [ALGOL  is a short for ALGOrithmic Language.]
!
! Author: Liang-You Peng, May 10, 2010
!
! This code is essentially a Fortran coding of the  
! original algorithm (in ALGOL 60) by  Walter Gautschi 
! " Algorithm 292: Regular Coulomb Wave Funcitons."
! References:
! (1) Communications of the ACM, V9, (1966) p793;
! (2) Communications of the ACM, V12, (1969) p280;
! (3) Communications of the ACM, V13, (1970) p573;
! (4) SIAM Review, V9, (1967) p24,  "Computational  
! Aspects of Three-Term Recurrence Relations", See Sect7   
! (5) Carl-Erik Fr\"oberg, Rev. Mod. Phys, V27 (1955), p399-411;
! "Numerical Treatment of Coulomb Wave Functions".
!================================================================
     
      SUBROUTINE Regular_Coulomb_Fun(eta, ro, Lmax,  F)

      USE f90_kind
      
      IMPLICIT NONE 
                      
      INTEGER, intent(in) :: Lmax
     
      REAL(double)  :: eta, ro, t 
     
      REAL(double)  :: F(0:Lmax)
      
      REAL(double), PARAMETER  :: Ebase = 2.718281828459045235D0
      
      REAL(double), PARAMETER  :: EInv = 0.36787944117144232159D0

!      REAL(double), PARAMETER  :: PI    = 3.141592653589793238D0
      
      INTEGER :: L, nu, nu1, mu, mu1, i,   d 
      
      REAL(double) :: scaling
      
      REAL(double)  :: epsilon, ro1, eta2, omega, d1, sum,  
     &           r, r1, s, t1, t2
      
      
      REAL(double) :: lmin(0:UB),Rr(0:Lmax),Rra(0:UB) 
    
      REAL(double) :: lam(0:UB),lambda(0:UB),Fapprox(0:UB)
      
      
            
      IF (ro<0) STOP 'Ro is negative!'
      
      IF (ro<1.D-16) THEN
         
         F = 0.D0
	 
	    RETURN
	 	 
      ENDIF 

! Scale for prevention of overflow:
      
      scaling = 1.D0**(-300)


! The precision for the calculation:
      
      d = 14
      
      epsilon =  1.0D-14
      
       

! The precision for the calculation.
      
      
!      ro1 = 1.D0/ro
      

      eta2 = eta**2
                        
         
      IF (eta>0) THEN
      
         t1 = 0.5D0 * ro / eta
      
      ELSE 
       
         t1 = 0.D0 
	 
      ENDIF 
              
	            
      IF (eta<1.D0) THEN
      
         omega = 0.0D0 
      
      ELSE IF(t1>=1.D0) THEN 
       
         omega = 0.5D0 * PI / t1
	 
      ELSE
      
         omega = ( 0.5D0*PI -  DATAN(DSQRT(1.D0/t1-1.D0)) +
     &           DSQRT(t1*(1.D0-t1)) ) / t1
      
      ENDIF


      
      lambda(0) = scaling  
  
      lmin(0)   = 1.D0 
       
      lambda(1) = scaling* ( omega - eta ) 
      
      sum       = ro * DEXP(omega*ro) * scaling


! Initialization:
    
      Fapprox = 0.D0 
      
!======================================
! log_e 10 =  2.30258509299404568D0
! 2 log_e 2 = 1.3862943611198906188D0
! e /2 =  2.7182818284590452353602875D0 /2  
!----------------------------                          
      
      d1 = 2.30258509299404568D0 * DBLE(d) +1.3862943611198906188D0
      
      t1 = 1.35914091422952261768D0 * ro

      IF ( Lmax<NINT(t1) ) THEN
      
         L = 1 + NINT(t1)
	 
      ELSE
      
         L = Lmax
	 
      ENDIF
      
      
      t1 = DEXP(PI*0.5D0*eta)
      
      s = DSQRT(1.D0+omega**2)
      
     
      IF (ABS(omega)<1.D-16) THEN
      
         t1 = t1 + 1/ t1
	 
      ELSE
      
         t1 = DEXP( -eta * DATAN(1.D0/omega) )
	 
      ENDIF 

      
      t2 = omega  + s 
      
      r  = 1.35914091422952261768D0 * ro * t2
      
      s  = ( d1 + DLOG(t1*DSQRT(t2/s)) -omega * ro ) / r

! The subroutine  inver_xlogx(y, x) calculates the value  x of 
! inverse  function of x Log x, for given value of y.
! The method is by numerical searching the root for x Log (x) - y = 0.      
    
      IF (s >=-EInv) THEN
      
        CALL inver_xlogx(s,  t) 
      
        nu = NINT( r * t )
      
      ELSE
      
        nu = NINT(r/Ebase) 
	
      ENDIF 

      
      CALL inver_xlogx(0.5D0*d1/DBLE(L),  t) 
                
     
      nu1 = NINT( DBLE(L) * t )
                
      IF (nu<nu1) nu = nu1
                            
      
      nu1 = 1 
      
      IF (ABS(omega)<1.D-14) THEN
      
         i = 1 ! L2
	 
      ELSE
      
         i = 2 ! L1;  i =3 ==> M1
	 
      ENDIF 

      
10    IF (i==1) THEN  ! L0, 10  
      
         GOTO 200 ! i = 1, L2
	 
      ELSE IF(i==2) THEN
      
         GOTO 100 ! i = 2, L1
	 
      ELSE  IF (i==3) THEN
      
         GOTO 300 ! i = 3, M1
	 
      ENDIF 
               
                                            
100   d1 = 2.D0 * eta / ( DEXP(2.D0*eta*DATAN(1.D0/omega)) - 1.D0 )
      

      
      i = 3 
      
300   mu1 = 0 ! 300, M1, i = 3
                 
      mu = NINT(1.25D0*nu)
      
      lam(mu1:nu) = 0.D0 
      

220   r = 0.D0  ! 220, M2
      
                     
      DO L = mu, mu1+ 1, -1 
      
         r = - ( DBLE(L**2)+ eta2 ) /  
     &     ( DBLE(L) * (  DBLE(2*L+1) * omega - DBLE(L+1) * r ) )
     
         IF (L<=nu) Rra(L-1) = r 
	 
      ENDDO
      
             
      
      DO L = mu1 + 1, nu
               
          lmin(L) = Rra(L-1) * lmin(L-1) 
	
	  IF (ABS(t1)<1.D-15) lmin(L) = 0.D0  	 	 	 	  
      
      ENDDO 
             
              
      
      DO L = mu1, nu
       
         IF ( ABS(lmin(L)-lam(L)) > epsilon * ABS(lmin(L)) ) THEN 
      
            lam(mu1:nu) = lmin(mu1:nu)
	    
	    mu = mu + 5 
	    
	    IF (mu<5*nu) THEN
	    
	       GOTO 220 ! M2
	       
	    ELSE
	    
	       STOP ' Convergence difficulty in Lambda[L]'
	       
	    ENDIF 
	    
	 
	 ENDIF 
      
      ENDDO 
              
	 
      r1     = 1.D-10	      
      
      lam(0) = - r1 * scaling 
      
      lam(1) = scaling  
      
      t1     = d1 / ( 1.D0 + r1**2 )
      
      
      DO L = 2, nu
       
        lam(L)   = (2.D0 - 1.d0/DBLE(L)) * omega * lam(L-1) 
     & + ( 1.d0-1.D0/DBLE(L) + eta2/DBLE(L-1)/DBLE(L) ) * lam(L-2)  
     
     
        lambda(L)= lmin(L) * scaling + t1 * (lam(L) 
     &            + scaling * r1 * lmin(L))        
     
      ENDDO 
      
     
      GOTO 330 ! L3
                       
    
200   DO L = nu1, nu-1
            
         lambda(L+1) = (2.D0+1.D0/DBLE(L)) / (1.D0+1.D0/DBLE(L))
     &   * omega * lambda(L) + ( 1.d0/(1.d0 + 1.d0/DBLE(L)) + 	 
     &   eta2/ DBLE(L)  / DBLE(L+1) ) * lambda(L-1)

      ENDDO
      
  
330   r = 0.D0 
      
      s = 0.D0 
      
       DO L = nu, 1, -1
       
         t1 = eta / DBLE(L+1)
	 
	 r  = 1.D0 / ( DBLE(2*L-1) * ( t1 / DBLE(L) + 1d0/ro   !ro1 
     &       - (1.D0 + t1**2) * r / DBLE(2*L+3) ) ) 
     
         s  = r * ( Lambda(L) + s )  
	 
	 IF (L<=Lmax) Rr(L-1) = r 
	 
      ENDDO
      
 
      F(0) = sum / (scaling + s)
      
      DO L = 1, Lmax 
      
         F(L) = Rr(L-1) * F(L-1) 
	 
      ENDDO 
      
      
      DO L = 0, Lmax 
      
         IF ( ABS(F(L)-Fapprox(L)) > epsilon * ABS(F(L)) ) THEN 
	 
	    Fapprox(0:Lmax) = F(0:Lmax)
	  	  
	    mu1 = nu
	  
	    nu1 = mu1 
	  
	    nu = nu + 10
	  	  	  
	    IF (nu<UB) THEN
	  
	       GOTO 10 !L0
	      
	    ELSE
	   
	       STOP 'convergence difficulty in Coulomb!' 
	      
	    ENDIF
	  
	 ENDIF 
	   
      ENDDO 
      
    
      t1 = 2.D0*PI * eta
      
      IF (ABS(t1)<1) THEN 
       
         s  = 1.d0
	 
	 t2 = s
	 
	 L  = 1

400	 L = L + 1

         t2 = t1 * t2 / DBLE(L)
	 
	 s = s + t2
	 
	 IF ( ABS(t2)>epsilon*ABS(s) ) GOTO 400 ! L4
	 
	 s = DSQRT(1.D0/s)
      
      
      ELSE
      
         s = DEXP(-t1/4.D0) / DSQRT( ( DEXP(t1/2.D0)
     &      - DEXP(-t1/2.D0) ) / t1 )
	       
      ENDIF 
           
     
      F(0) = s * F(0) 
      
      
      DO L = 1, Lmax 
            
       s = (1.d0 - 1.d0/DBLE(2*L)) / (1.d0 + 1.d0/DBLE(2*L) ) 
     &         * DSQRT( 1.d0 + (eta/DBLE(L))**2 ) * s     
 
       F(L) = s * F(L) 
       
      ENDDO 
      
                                  
      END SUBROUTINE Regular_Coulomb_Fun ! L5








! The inverse of function y = t * log (t), for given y:

! Note, when x<1, there are two solutions, we take
! the larger one in this program

	
      SUBROUTINE inver_xlogx(y,  t) 
      
      USE f90_kind
       
      IMPLICIT NONE 
      
      external fun 
      
      REAL(double) :: fun 
      
      REAL(double) :: t, y 
      
      REAL(double)::   EInv
      
      REAL(double) :: x0, x2, x1, eps 
      
      INTEGER :: imax
      
      REAL(double) ::  amybshec
      
!      common/aPARAMETER/ amybshec  
      
      

    
      EInv = 0.36787944117144232159D0

! Set precision needed for the root finder:
     
      eps = 1.D-15
     
     
      amybshec = y 
      
 
      IF (amybshec < -EInv ) STOP 'y<-1/e in FUNCTION t(y)'
      

! Only one zero for the function for x*ln (x), at x = 1.d0
      
      IF ( ABS(amybshec) < 1.D-14) THEN
      
          t = 1.D0 	  	
	  
      ENDIF 
      
! If amybshec> 0, then the root is unique (x>1.0):
      
      IF ( amybshec > 0.D0 ) THEN
      
          x0 = 1.D0 + 1.D-10
	  
	  x2 = 1.D20

! Max iteration permitted:
	  
	  imax = UB
	  	              
          CALL iterate(x0, x2, fun,amybshec, eps, imax, x1)
	  
	  t = x1 
	  
      ENDIF 
      
! IF  amybshec < 0.D0  .AND.  amybshec >= -1/e, there are two solutions:
       
    
       IF ( amybshec < 0.D0  .AND.  amybshec >= -EInv ) THEN

! The smaller soluttion:
      
!!          x0 =   1.D-20
	  
!!	  x2 =   EInv
	  
!!	  imax = 8000000
           
!!         CALL iterate(x0, x2, fun, eps, imax, x1)
!!	  
!!	  t = x1 
	  	  	  
! The larger solution
      
          x0 =   EInv
	  
	  x2 =   1.D0 - 1.D-10
	  
	  imax = UB
           
          CALL iterate(x0, x2, fun,amybshec, eps, imax, x1)
	     
	  t = x1
	  
       ENDIF 
    
      END SUBROUTINE inver_xlogx 
	
      
        function fun(x,amybshec)

        USE f90_kind

 		
	implicit none 
	
	REAL(double) :: fun, x 
	
	REAL(double) ::  amybshec
      
!        common/aPARAMETER/ amybshec  
      
	
	fun = x * DLOG(x) - amybshec 
	
	END  function fun



! The world's best root finder alogrithm.
! http://www.embedded.com/story/OEG20030508S0030
! By Jack W. Crenshaw 
! Jack Crenshaw is a senior software engineer at Spectrum-Astro and 
! the author of Math Toolkit for Real-Time Programming, from CMP 
! Books. He holds a PhD in physics from Auburn University. 
! E-mail him at jcrens@earthlink.net. 

! Convert from C to Fortran by Liang-You Peng

! Make sure that there is ONLY one root between x0 and x2 
! Therefore, 
      
        subroutine iterate(x0, x2, fun,amybshec, eps, imax, x1)

        USE f90_kind

        implicit none 
  
        external fun
  
        REAL(double) :: fun 

        REAL(double) :: x1, y0, y1, y2, b, c, temp, y10, y20, y21,
     &            xm, ym, eps, x0,x2
 
        REAL(double) :: xmlast 
  
        INTEGER:: imax, i
      
        REAL(double) ::  amybshec

        
	xmlast = x0 
  
        y0 = fun(x0,amybshec)

! Is x0 a root?
  
        if ( abs(y0)<1.D-15 ) then
         
	   x1 = x0
   
           return
   
        ENDif 
  
        y2 = fun(x2,amybshec)
	
! Is x2 a root?   
 
        if ( abs(y0)<1.D-15 ) then
  
           x1 = x2
  
           return 
  
        ENDif 

! Are y2 and y0 have the same sign?

        if (y2 * y0 > 0.0) THEN
  
           stop 'No zero OR at least two zeros between x0 and x2.'
  
        ENDif

! Now do the iteration:


        do i = 0, imax 
  
           x1 = 0.5d0 * (x2 + x0)
           
	   y1 = fun(x1,amybshec)
   
           if ( abs(y1)<1.D-16 .OR. abs(x1 - x0) < eps) then
	   
              return 
          
	   ENDif 
                 
           if (y1 * y0 > 0.d0) then 
    
               temp = x0
               x0   = x2
               x2   = temp
               temp = y0
               y0   = y2
               y2   = temp
      
           ENDif 
    
           y10 = y1 - y0
           y21 = y2 - y1
           y20 = y2 - y0
    
           
	   if (y2 * y20 <2.d0 * y1 * y10) then 
    
              x2 = x1
              y2 = y1
    
           else 
    
              b = (x1 - x0) / y10   
              c = (y10 -y21) / (y21 * y20) 
              xm = x0 - b * y0 * (1.d0-c * y1)
              ym = fun(xm,amybshec)
      
              if (abs(ym)<1.D-16 .OR. abs (xm - xmlast) < eps) then
                 x1 = xm
             	 return 
              ENDif 
                 
              xmlast = xm
      
              if (ym * y0 < 0.d0) then 
      
                 x2 = xm
                 y2 = ym
       
              else
      
                 x0 = xm
                 y0 = ym
                 x2 = x1
                 y2 = y1
     
              ENDif  ! if (ym * y0 < 0.d0)
   
           ENDif !  if (y2 * y20 <2.d0 * y1 * y10)
    
         ENDdo ! do i = 0, imax 
            
       END subroutine iterate
         

!============================================
!
! Peng, 2010.05.20.
!
! Another way to calculate the Coulomb phase:
! 
! Given momentum xk, L, and   charge ZE,
! the code returns the Coulomb phase for all the
! L's: 0, 1, 2, ..., L,
! in the matrix DeltaL(0:L)
!
!=============================================

      SUBROUTINE sigma_arb(xk, l, Zar, DeltaL)
      
        USE f90_kind
      
      IMPLICIT REAL(double) (a-h,o-z)
      
      REAL(double) :: DeltaL(0:L), sigma0
      
      INTEGER :: LL
      
      
!Computes sigma function appearing in coulomb phase for l=0
!according to Martins, JPB 1,154 (1968):

      c=100.d0-1.d0/(xk*xk)
      IF (c) 10,10,20
 10   xl=0.d0
      GOTO 30
 20   xl=idint(dsqrt(c))
 30   q=xl+1.d0
      z=Zar/xk
      a=-0.5d0*z*dlog(q*q+z*z)+z
      xx=z/q
      b=-(xl+0.5d0)*datan(xx)
      c=z/(12.d0*(q*q+z*z))
      d=z*(z*z-3.d0*q*q)/(360.d0*(q*q+z*z)**3)
      sigma0=a+b+c+d
 35   IF (xl) 50,40,50
 40   GOTO 549
 50   xy=z/xl
      sigma0=sigma0+datan(xy)
      xl=xl-1.d0
      GOTO 35
549   CONTINUE

      LL = NINT(xl) 
       
      DeltaL(0) = sigma0

351   IF (dble(l-LL)) 501,401,501
401   RETURN
501   LL=LL+1.d0
      xy=z/dble(LL)
      sigma0=sigma0-datan(xy)
      
      DeltaL(LL) = sigma0
      
      GOTO 351

      END SUBROUTINE sigma_arb



