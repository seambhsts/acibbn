subroutine CBESK (Z, FNU, KODE, N, CY, NZ, IERR)
!
!! CBESK computes a sequence of the Bessel functions K(a,z) for ...
!            complex argument z and real nonnegative orders a=b,b+1, ...
!            b+2,... where b>0.  A scaling option is available to ...
!            help avoid overflow.
!
!***LIBRARY   SLATEC
!***CATEGORY  C10B4
!***TYPE      COMPLEX (CBESK-C, ZBESK-C)
!***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT, K BESSEL FUNCTIONS,
!             MODIFIED BESSEL FUNCTIONS
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!         On KODE=1, CBESK computes an N member sequence of complex
!         Bessel functions CY(L)=K(FNU+L-1,Z) for real nonnegative
!         orders FNU+L-1, L=1,...,N and complex Z /= 0 in the cut
!         plane -pi<arg(Z)<=pi.  On KODE=2, CBESJ returns the scaled
!         functions
!
!            CY(L) = exp(Z)*K(FNU+L-1,Z),  L=1,...,N
!
!         which remove the exponential growth in both the left and
!         right half planes as Z goes to infinity.  Definitions and
!         notation are found in the NBS Handbook of Mathematical
!         Functions (Ref. 1).
!
!         Input
!           Z      - Nonzero argument of type COMPLEX
!           FNU    - Initial order of type REAL, FNU>=0
!           KODE   - A parameter to indicate the scaling option
!                    KODE=1  returns
!                            CY(L)=K(FNU+L-1,Z), L=1,...,N
!                        =2  returns
!                            CY(L)=K(FNU+L-1,Z)*EXP(Z), L=1,...,N
!           N      - Number of terms in the sequence, N>=1
!
!         Output
!           CY     - Result vector of type COMPLEX
!           NZ     - Number of underflows set to zero
!                    NZ=0    Normal return
!                    NZ>0    CY(L)=0 for NZ values of L (if Re(Z)>0
!                            then CY(L)=0 for L=1,...,NZ; in the
!                            complementary half plane the underflows
!                            may not be in an uninterrupted sequence)
!           IERR   - Error flag
!                    IERR=0  Normal return     - COMPUTATION COMPLETED
!                    IERR=1  Input error       - NO COMPUTATION
!                    IERR=2  Overflow          - NO COMPUTATION
!                            (abs(Z) too small and/or FNU+N-1
!                            too large)
!                    IERR=3  Precision warning - COMPUTATION COMPLETED
!                            (Result has half precision or less
!                            because abs(Z) or FNU+N-1 is large)
!                    IERR=4  Precision error   - NO COMPUTATION
!                            (Result has no precision because
!                            abs(Z) or FNU+N-1 is too large)
!                    IERR=5  Algorithmic error - NO COMPUTATION
!                            (Termination condition not met)
!
! *Long Description:
!
!         Equations of the reference are implemented to compute K(a,z)
!         for small orders a and a+1 in the right half plane Re(z)>=0.
!         Forward recurrence generates higher orders.  The formula
!
!            K(a,z*exp((t)) = exp(-t)*K(a,z) - t*I(a,z),  Re(z)>0
!                         t = i*pi or -i*pi
!
!         continues K to the left half plane.
!
!         For large orders, K(a,z) is computed by means of its uniform
!         asymptotic expansion.
!
!         For negative orders, the formula
!
!            K(-a,z) = K(a,z)
!
!         can be used.
!
!         CBESK assumes that a significant digit sinh function is
!         available.
!
!         In most complex variable computation, one must evaluate ele-
!         mentary functions.  When the magnitude of Z or FNU+N-1 is
!         large, losses of significance by argument reduction occur.
!         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
!         losses exceeding half precision are likely and an error flag
!         IERR=3 is triggered where UR=R1MACH(4)=UNIT ROUNDOFF.  Also,
!         if either is larger than U2=0.5/UR, then all significance is
!         lost and IERR=4.  In order to use the INT function, arguments
!         must be further restricted not to exceed the largest machine
!         integer, U3=I1MACH(9).  Thus, the magnitude of Z and FNU+N-1
!         is restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, and
!         U3 approximate 2.0E+3, 4.2E+6, 2.1E+9 in single precision
!         and 4.7E+7, 2.3E+15 and 2.1E+9 in double precision.  This
!         makes U2 limiting in single precision and U3 limiting in
!         double precision.  This means that one can expect to retain,
!         in the worst cases on IEEE machines, no digits in single pre-
!         cision and only 6 digits in double precision.  Similar con-
!         siderations hold for other machines.
!
!         The approximate relative error in the magnitude of a complex
!         Bessel function can be expressed as P*10**S where P=MAX(UNIT
!         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre-
!         sents the increase in error due to argument reduction in the
!         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))),
!         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF
!         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may
!         have only absolute accuracy.  This is most likely to occur
!         when one component (in magnitude) is larger than the other by
!         several orders of magnitude.  If one component is 10**K larger
!         than the other, then one can expect only MAX(ABS(LOG10(P))-K,
!         0) significant digits; or, stated another way, when K exceeds
!         the exponent of P, no significant digits remain in the smaller
!         component.  However, the phase angle retains absolute accuracy
!         because, in complex arithmetic with precision P, the smaller
!         component will not (as a rule) decrease below P times the
!         magnitude of the larger component.  In these extreme cases,
!         the principal phase angle is on the order of +P, -P, PI/2-P,
!         or -PI/2+P.
!
!***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
!                 matical Functions, National Bureau of Standards
!                 Applied Mathematics Series 55, U. S. Department
!                 of Commerce, Tenth Printing (1972) or later.
!               2. D. E. Amos, Computation of Bessel Functions of
!                 Complex Argument, Report SAND83-0086, Sandia National
!                 Laboratories, Albuquerque, NM, May 1983.
!               3. D. E. Amos, Computation of Bessel Functions of
!                 Complex Argument and Large Order, Report SAND83-0643,
!                 Sandia National Laboratories, Albuquerque, NM, May
!                 1983.
!               4. D. E. Amos, A Subroutine Package for Bessel Functions
!                 of a Complex Argument and Nonnegative Order, Report
!                 SAND85-1018, Sandia National Laboratory, Albuquerque,
!                 NM, May 1985.
!               5. D. E. Amos, A portable package for Bessel functions
!                 of a complex argument and nonnegative order, ACM
!                 Transactions on Mathematical Software, 12 (September
!                 1986), pp. 265-273.
!
!***ROUTINES CALLED  CACON, CBKNU, CBUNK, CUOIK, I1MACH, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   890801  REVISION DATE from Version 3.2
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   920128  Category corrected.  (WRB)
!   920811  Prologue revised.  (DWL)
!***END PROLOGUE  CBESK
!
  COMPLEX CY, Z
  REAL AA, ALIM, ALN, ARG, AZ, DIG, ELIM, FN, FNU, FNUL, RL, R1M5, &
   TOL, UFL, XX, YY, R1MACH, BB
  INTEGER IERR, K, KODE, K1, K2, MR, N, NN, NUF, NW, NZ, I1MACH
  DIMENSION CY(N)
!***FIRST EXECUTABLE STATEMENT  CBESK
  IERR = 0
  NZ=0
  XX = REAL(Z)
  YY = AIMAG(Z)
  if (YY == 0.0E0 .AND. XX == 0.0E0) IERR=1
  if (FNU < 0.0E0) IERR=1
  if (KODE < 1 .OR. KODE > 2) IERR=1
  if (N < 1) IERR=1
  if (IERR /= 0) RETURN
  NN = N
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
!-----------------------------------------------------------------------
  TOL = MAX(R1MACH(4),1.0E-18)
  K1 = I1MACH(12)
  K2 = I1MACH(13)
  R1M5 = R1MACH(5)
  K = MIN(ABS(K1),ABS(K2))
  ELIM = 2.303E0*(K*R1M5-3.0E0)
  K1 = I1MACH(11) - 1
  AA = R1M5*K1
  DIG = MIN(AA,18.0E0)
  AA = AA*2.303E0
  ALIM = ELIM + MAX(-AA,-41.45E0)
  FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
  RL = 1.2E0*DIG + 3.0E0
  AZ = ABS(Z)
  FN = FNU + (NN-1)
!-----------------------------------------------------------------------
!     TEST FOR RANGE
!-----------------------------------------------------------------------
  AA = 0.5E0/TOL
  BB=I1MACH(9)*0.5E0
  AA=MIN(AA,BB)
  if ( AZ > AA) go to 210
  if ( FN > AA) go to 210
  AA=SQRT(AA)
  if ( AZ > AA) IERR=3
  if ( FN > AA) IERR=3
!-----------------------------------------------------------------------
!     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
!-----------------------------------------------------------------------
!     UFL = EXP(-ELIM)
  UFL = R1MACH(1)*1.0E+3
  if (AZ < UFL) go to 180
  if (FNU > FNUL) go to 80
  if (FN <= 1.0E0) go to 60
  if (FN > 2.0E0) go to 50
  if (AZ > TOL) go to 60
  ARG = 0.5E0*AZ
  ALN = -FN*ALOG(ARG)
  if (ALN > ELIM) go to 180
  go to 60
   50 CONTINUE
  call CUOIK(Z, FNU, KODE, 2, NN, CY, NUF, TOL, ELIM, ALIM)
  if (NUF < 0) go to 180
  NZ = NZ + NUF
  NN = NN - NUF
!-----------------------------------------------------------------------
!     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
!     if NUF=NN, THEN CY(I)=CZERO FOR ALL I
!-----------------------------------------------------------------------
  if (NN == 0) go to 100
   60 CONTINUE
  if (XX < 0.0E0) go to 70
!-----------------------------------------------------------------------
!     RIGHT HALF PLANE COMPUTATION, REAL(Z) >= 0.
!-----------------------------------------------------------------------
  call CBKNU(Z, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
  if (NW < 0) go to 200
  NZ=NW
  return
!-----------------------------------------------------------------------
!     LEFT HALF PLANE COMPUTATION
!     PI/2 < ARG(Z) <= PI AND -PI < ARG(Z) < -PI/2.
!-----------------------------------------------------------------------
   70 CONTINUE
  if (NZ /= 0) go to 180
  MR = 1
  if (YY < 0.0E0) MR = -1
  call CACON(Z, FNU, KODE, MR, NN, CY, NW, RL, FNUL, TOL, ELIM, &
   ALIM)
  if (NW < 0) go to 200
  NZ=NW
  return
!-----------------------------------------------------------------------
!     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU > FNUL
!-----------------------------------------------------------------------
   80 CONTINUE
  MR = 0
  if (XX >= 0.0E0) go to 90
  MR = 1
  if (YY < 0.0E0) MR = -1
   90 CONTINUE
  call CBUNK(Z, FNU, KODE, MR, NN, CY, NW, TOL, ELIM, ALIM)
  if (NW < 0) go to 200
  NZ = NZ + NW
  return
  100 CONTINUE
  if (XX < 0.0E0) go to 180
  return
  180 CONTINUE
  NZ = 0
  IERR=2
  return
  200 CONTINUE
  if ( NW == (-1)) go to 180
  NZ=0
  IERR=5
  return
  210 CONTINUE
  NZ=0
  IERR=4
  return
end

subroutine CUNK1 (Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
!
!! CUNK1 is subsidiary to CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNK1-A, ZUNK1-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
!     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
!     UNIFORM ASYMPTOTIC EXPANSION.
!     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
!     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
!
!***SEE ALSO  CBESK
!***ROUTINES CALLED  CS1S2, CUCHK, CUNIK, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUNK1
  COMPLEX CFN, CK, CONE, CRSC, CS, CSCL, CSGN, CSPN, CSR, CSS, &
   CWRK, CY, CZERO, C1, C2, PHI,  RZ, SUM,  S1, S2, Y, Z, &
   ZETA1,  ZETA2,  ZR, PHID, ZETA1D, ZETA2D, SUMD
  REAL ALIM, ANG, APHI, ASC, ASCLE, BRY, CPN, C2I, C2M, C2R, ELIM, &
   FMR, FN, FNF, FNU, PI, RS1, SGN, SPN, TOL, X, R1MACH
  INTEGER I, IB, IFLAG, IFN, IL, INIT, INU, IUF, K, KDFLG, KFLAG, &
   KK, KODE, MR, N, NW, NZ, J, IPARD, INITD, IC, M
  DIMENSION BRY(3), INIT(2), Y(N), SUM(2), PHI(2), ZETA1(2), &
   ZETA2(2), CY(2), CWRK(16,3), CSS(3), CSR(3)
  DATA CZERO, CONE / (0.0E0,0.0E0) , (1.0E0,0.0E0) /
  DATA PI / 3.14159265358979324E0 /
!***FIRST EXECUTABLE STATEMENT  CUNK1
  KDFLG = 1
  NZ = 0
!-----------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
!     THE UNDERFLOW LIMIT
!-----------------------------------------------------------------------
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CRSC = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CRSC
  CSR(1) = CRSC
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = 1.0E+3*R1MACH(1)/TOL
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = R1MACH(2)
  X = REAL(Z)
  ZR = Z
  if (X < 0.0E0) ZR = -Z
  J=2
  DO 70 I=1,N
!-----------------------------------------------------------------------
!     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
!-----------------------------------------------------------------------
    J = 3 - J
    FN = FNU + (I-1)
    INIT(J) = 0
    call CUNIK(ZR, FN, 2, 0, TOL, INIT(J), PHI(J), ZETA1(J), &
     ZETA2(J), SUM(J), CWRK(1,J))
    if (KODE == 1) go to 20
    CFN = CMPLX(FN,0.0E0)
    S1 = ZETA1(J) - CFN*(CFN/(ZR+ZETA2(J)))
    go to 30
   20   CONTINUE
    S1 = ZETA1(J) - ZETA2(J)
   30   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = REAL(S1)
    if (ABS(RS1) > ELIM) go to 60
    if (KDFLG == 1) KFLAG = 2
    if (ABS(RS1) < ALIM) go to 40
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    APHI = ABS(PHI(J))
    RS1 = RS1 + ALOG(APHI)
    if (ABS(RS1) > ELIM) go to 60
    if (KDFLG == 1) KFLAG = 1
    if (RS1 < 0.0E0) go to 40
    if (KDFLG == 1) KFLAG = 3
   40   CONTINUE
!-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!-----------------------------------------------------------------------
    S2 = PHI(J)*SUM(J)
    C2R = REAL(S1)
    C2I = AIMAG(S1)
    C2M = EXP(C2R)*REAL(CSS(KFLAG))
    S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
    S2 = S2*S1
    if (KFLAG /= 1) go to 50
    call CUCHK(S2, NW, BRY(1), TOL)
    if (NW /= 0) go to 60
   50   CONTINUE
    CY(KDFLG) = S2
    Y(I) = S2*CSR(KFLAG)
    if (KDFLG == 2) go to 75
    KDFLG = 2
    go to 70
   60   CONTINUE
    if (RS1 > 0.0E0) go to 290
!-----------------------------------------------------------------------
!     FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
    if (X < 0.0E0) go to 290
    KDFLG = 1
    Y(I) = CZERO
    NZ=NZ+1
    if (I == 1) go to 70
    if (Y(I-1) == CZERO) go to 70
    Y(I-1) = CZERO
    NZ=NZ+1
   70 CONTINUE
  I=N
   75 CONTINUE
  RZ = CMPLX(2.0E0,0.0E0)/ZR
  CK = CMPLX(FN,0.0E0)*RZ
  IB = I+1
  if (N < IB) go to 160
!-----------------------------------------------------------------------
!     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
!     ON UNDERFLOW
!-----------------------------------------------------------------------
  FN = FNU+(N-1)
  IPARD = 1
  if (MR /= 0) IPARD = 0
  INITD = 0
  call CUNIK(ZR,FN,2,IPARD,TOL,INITD,PHID,ZETA1D,ZETA2D,SUMD, &
  CWRK(1,3))
  if (KODE == 1) go to 80
  CFN=CMPLX(FN,0.0E0)
  S1=ZETA1D-CFN*(CFN/(ZR+ZETA2D))
  go to 90
   80 CONTINUE
  S1=ZETA1D-ZETA2D
   90 CONTINUE
  RS1=REAL(S1)
  if (ABS(RS1) > ELIM) go to 95
  if (ABS(RS1) < ALIM) go to 100
!-----------------------------------------------------------------------
!     REFINE ESTIMATE AND TEST
!-----------------------------------------------------------------------
  APHI=ABS(PHID)
  RS1=RS1+ALOG(APHI)
  if (ABS(RS1) < ELIM) go to 100
   95 CONTINUE
  if (RS1 > 0.0E0) go to 290
!-----------------------------------------------------------------------
!     FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
  if (X < 0.0E0) go to 290
  NZ=N
  DO 96 I=1,N
    Y(I) = CZERO
   96 CONTINUE
  return
  100 CONTINUE
!-----------------------------------------------------------------------
!     RECUR FORWARD FOR REMAINDER OF THE SEQUENCE
!-----------------------------------------------------------------------
  S1 = CY(1)
  S2 = CY(2)
  C1 = CSR(KFLAG)
  ASCLE = BRY(KFLAG)
  DO 120 I=IB,N
    C2 = S2
    S2 = CK*S2 + S1
    S1 = C2
    CK = CK + RZ
    C2 = S2*C1
    Y(I) = C2
    if (KFLAG >= 3) go to 120
    C2R = REAL(C2)
    C2I = AIMAG(C2)
    C2R = ABS(C2R)
    C2I = ABS(C2I)
    C2M = MAX(C2R,C2I)
    if (C2M <= ASCLE) go to 120
    KFLAG = KFLAG + 1
    ASCLE = BRY(KFLAG)
    S1 = S1*C1
    S2 = C2
    S1 = S1*CSS(KFLAG)
    S2 = S2*CSS(KFLAG)
    C1 = CSR(KFLAG)
  120 CONTINUE
  160 CONTINUE
  if (MR == 0) RETURN
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION FOR RE(Z) < 0.0E0
!-----------------------------------------------------------------------
  NZ = 0
  FMR = MR
  SGN = -SIGN(PI,FMR)
!-----------------------------------------------------------------------
!     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
!-----------------------------------------------------------------------
  CSGN = CMPLX(0.0E0,SGN)
  INU = FNU
  FNF = FNU - INU
  IFN = INU + N - 1
  ANG = FNF*SGN
  CPN = COS(ANG)
  SPN = SIN(ANG)
  CSPN = CMPLX(CPN,SPN)
  if (MOD(IFN,2) == 1) CSPN = -CSPN
  ASC = BRY(1)
  KK = N
  IUF = 0
  KDFLG = 1
  IB = IB-1
  IC = IB-1
  DO 260 K=1,N
    FN = FNU + (KK-1)
!-----------------------------------------------------------------------
!     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
!     FUNCTION ABOVE
!-----------------------------------------------------------------------
    M=3
    if (N > 2) go to 175
  170   CONTINUE
    INITD = INIT(J)
    PHID = PHI(J)
    ZETA1D = ZETA1(J)
    ZETA2D = ZETA2(J)
    SUMD = SUM(J)
    M = J
    J = 3 - J
    go to 180
  175   CONTINUE
    if ((KK == N).AND.(IB < N)) go to 180
    if ((KK == IB).OR.(KK == IC)) go to 170
    INITD = 0
  180   CONTINUE
    call CUNIK(ZR, FN, 1, 0, TOL, INITD, PHID, ZETA1D, &
     ZETA2D, SUMD, CWRK(1,M))
    if (KODE == 1) go to 190
    CFN = CMPLX(FN,0.0E0)
    S1 = -ZETA1D + CFN*(CFN/(ZR+ZETA2D))
    go to 200
  190   CONTINUE
    S1 = -ZETA1D + ZETA2D
  200   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = REAL(S1)
    if (ABS(RS1) > ELIM) go to 250
    if (KDFLG == 1) IFLAG = 2
    if (ABS(RS1) < ALIM) go to 210
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    APHI = ABS(PHID)
    RS1 = RS1 + ALOG(APHI)
    if (ABS(RS1) > ELIM) go to 250
    if (KDFLG == 1) IFLAG = 1
    if (RS1 < 0.0E0) go to 210
    if (KDFLG == 1) IFLAG = 3
  210   CONTINUE
    S2 = CSGN*PHID*SUMD
    C2R = REAL(S1)
    C2I = AIMAG(S1)
    C2M = EXP(C2R)*REAL(CSS(IFLAG))
    S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
    S2 = S2*S1
    if (IFLAG /= 1) go to 220
    call CUCHK(S2, NW, BRY(1), TOL)
    if (NW /= 0) S2 = CMPLX(0.0E0,0.0E0)
  220   CONTINUE
    CY(KDFLG) = S2
    C2 = S2
    S2 = S2*CSR(IFLAG)
!-----------------------------------------------------------------------
!     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
!-----------------------------------------------------------------------
    S1 = Y(KK)
    if (KODE == 1) go to 240
    call CS1S2(ZR, S1, S2, NW, ASC, ALIM, IUF)
    NZ = NZ + NW
  240   CONTINUE
    Y(KK) = S1*CSPN + S2
    KK = KK - 1
    CSPN = -CSPN
    if (C2 /= CZERO) go to 245
    KDFLG = 1
    go to 260
  245   CONTINUE
    if (KDFLG == 2) go to 265
    KDFLG = 2
    go to 260
  250   CONTINUE
    if (RS1 > 0.0E0) go to 290
    S2 = CZERO
    go to 220
  260 CONTINUE
  K = N
  265 CONTINUE
  IL = N - K
  if (IL == 0) RETURN
!-----------------------------------------------------------------------
!     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
!     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
!     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
!-----------------------------------------------------------------------
  S1 = CY(1)
  S2 = CY(2)
  CS = CSR(IFLAG)
  ASCLE = BRY(IFLAG)
  FN = (INU+IL)
  DO 280 I=1,IL
    C2 = S2
    S2 = S1 + CMPLX(FN+FNF,0.0E0)*RZ*S2
    S1 = C2
    FN = FN - 1.0E0
    C2 = S2*CS
    CK = C2
    C1 = Y(KK)
    if (KODE == 1) go to 270
    call CS1S2(ZR, C1, C2, NW, ASC, ALIM, IUF)
    NZ = NZ + NW
  270   CONTINUE
    Y(KK) = C1*CSPN + C2
    KK = KK - 1
    CSPN = -CSPN
    if (IFLAG >= 3) go to 280
    C2R = REAL(CK)
    C2I = AIMAG(CK)
    C2R = ABS(C2R)
    C2I = ABS(C2I)
    C2M = MAX(C2R,C2I)
    if (C2M <= ASCLE) go to 280
    IFLAG = IFLAG + 1
    ASCLE = BRY(IFLAG)
    S1 = S1*CS
    S2 = CK
    S1 = S1*CSS(IFLAG)
    S2 = S2*CSS(IFLAG)
    CS = CSR(IFLAG)
  280 CONTINUE
  return
  290 CONTINUE
  NZ = -1
  return
end
FUNCTION GAMLN (Z, IERR)
!
!! GAMLN computes the logarithm of the Gamma function.
!
!***LIBRARY   SLATEC
!***CATEGORY  C7A
!***TYPE      SINGLE PRECISION (GAMLN-S, DGAMLN-D)
!***KEYWORDS  LOGARITHM OF GAMMA FUNCTION
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!         GAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
!         Z > 0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
!         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
!         G(Z+1)=Z*G(Z) FOR Z <= ZMIN.  THE FUNCTION WAS MADE AS
!         PORTABLE AS POSSIBLE BY COMPUTING ZMIN FROM THE NUMBER OF BASE
!         10 DIGITS IN A WORD, RLN=MAX(-ALOG10(R1MACH(4)),0.5E-18)
!         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.
!
!         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
!         VALUES IS USED FOR SPEED OF EXECUTION.
!
!     DESCRIPTION OF ARGUMENTS
!
!         INPUT
!           Z      - REAL ARGUMENT, Z > 0.0E0
!
!         OUTPUT
!           GAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z
!           IERR   - ERROR FLAG
!                    IERR=0, NORMAL RETURN, COMPUTATION COMPLETED
!                    IERR=1, Z <= 0.0E0,    NO COMPUTATION
!
!***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
!***ROUTINES CALLED  I1MACH, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   830501  REVISION DATE from Version 3.2
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   920128  Category corrected.  (WRB)
!   921215  GAMLN defined for Z negative.  (WRB)
!***END PROLOGUE  GAMLN
!
  REAL GAMLN
  INTEGER I, I1M, K, MZ, NZ, IERR, I1MACH
  REAL CF, CON, FLN, FZ, GLN, RLN, S, TLG, TRM, TST, T1, WDTOL, Z, &
   ZDMY, ZINC, ZM, ZMIN, ZP, ZSQ
  REAL R1MACH
  DIMENSION CF(22), GLN(100)
!           LNGAMMA(N), N=1,100
  DATA GLN(1), GLN(2), GLN(3), GLN(4), GLN(5), GLN(6), GLN(7), &
       GLN(8), GLN(9), GLN(10), GLN(11), GLN(12), GLN(13), GLN(14), &
       GLN(15), GLN(16), GLN(17), GLN(18), GLN(19), GLN(20), &
       GLN(21), GLN(22)/ &
       0.00000000000000000E+00,     0.00000000000000000E+00, &
       6.93147180559945309E-01,     1.79175946922805500E+00, &
       3.17805383034794562E+00,     4.78749174278204599E+00, &
       6.57925121201010100E+00,     8.52516136106541430E+00, &
       1.06046029027452502E+01,     1.28018274800814696E+01, &
       1.51044125730755153E+01,     1.75023078458738858E+01, &
       1.99872144956618861E+01,     2.25521638531234229E+01, &
       2.51912211827386815E+01,     2.78992713838408916E+01, &
       3.06718601060806728E+01,     3.35050734501368889E+01, &
       3.63954452080330536E+01,     3.93398841871994940E+01, &
       4.23356164607534850E+01,     4.53801388984769080E+01/
  DATA GLN(23), GLN(24), GLN(25), GLN(26), GLN(27), GLN(28), &
       GLN(29), GLN(30), GLN(31), GLN(32), GLN(33), GLN(34), &
       GLN(35), GLN(36), GLN(37), GLN(38), GLN(39), GLN(40), &
       GLN(41), GLN(42), GLN(43), GLN(44)/ &
       4.84711813518352239E+01,     5.16066755677643736E+01, &
       5.47847293981123192E+01,     5.80036052229805199E+01, &
       6.12617017610020020E+01,     6.45575386270063311E+01, &
       6.78897431371815350E+01,     7.12570389671680090E+01, &
       7.46582363488301644E+01,     7.80922235533153106E+01, &
       8.15579594561150372E+01,     8.50544670175815174E+01, &
       8.85808275421976788E+01,     9.21361756036870925E+01, &
       9.57196945421432025E+01,     9.93306124547874269E+01, &
       1.02968198614513813E+02,     1.06631760260643459E+02, &
       1.10320639714757395E+02,     1.14034211781461703E+02, &
       1.17771881399745072E+02,     1.21533081515438634E+02/
  DATA GLN(45), GLN(46), GLN(47), GLN(48), GLN(49), GLN(50), &
       GLN(51), GLN(52), GLN(53), GLN(54), GLN(55), GLN(56), &
       GLN(57), GLN(58), GLN(59), GLN(60), GLN(61), GLN(62), &
       GLN(63), GLN(64), GLN(65), GLN(66)/ &
       1.25317271149356895E+02,     1.29123933639127215E+02, &
       1.32952575035616310E+02,     1.36802722637326368E+02, &
       1.40673923648234259E+02,     1.44565743946344886E+02, &
       1.48477766951773032E+02,     1.52409592584497358E+02, &
       1.56360836303078785E+02,     1.60331128216630907E+02, &
       1.64320112263195181E+02,     1.68327445448427652E+02, &
       1.72352797139162802E+02,     1.76395848406997352E+02, &
       1.80456291417543771E+02,     1.84533828861449491E+02, &
       1.88628173423671591E+02,     1.92739047287844902E+02, &
       1.96866181672889994E+02,     2.01009316399281527E+02, &
       2.05168199482641199E+02,     2.09342586752536836E+02/
  DATA GLN(67), GLN(68), GLN(69), GLN(70), GLN(71), GLN(72), &
       GLN(73), GLN(74), GLN(75), GLN(76), GLN(77), GLN(78), &
       GLN(79), GLN(80), GLN(81), GLN(82), GLN(83), GLN(84), &
       GLN(85), GLN(86), GLN(87), GLN(88)/ &
       2.13532241494563261E+02,     2.17736934113954227E+02, &
       2.21956441819130334E+02,     2.26190548323727593E+02, &
       2.30439043565776952E+02,     2.34701723442818268E+02, &
       2.38978389561834323E+02,     2.43268849002982714E+02, &
       2.47572914096186884E+02,     2.51890402209723194E+02, &
       2.56221135550009525E+02,     2.60564940971863209E+02, &
       2.64921649798552801E+02,     2.69291097651019823E+02, &
       2.73673124285693704E+02,     2.78067573440366143E+02, &
       2.82474292687630396E+02,     2.86893133295426994E+02, &
       2.91323950094270308E+02,     2.95766601350760624E+02, &
       3.00220948647014132E+02,     3.04686856765668715E+02/
  DATA GLN(89), GLN(90), GLN(91), GLN(92), GLN(93), GLN(94), &
       GLN(95), GLN(96), GLN(97), GLN(98), GLN(99), GLN(100)/ &
       3.09164193580146922E+02,     3.13652829949879062E+02, &
       3.18152639620209327E+02,     3.22663499126726177E+02, &
       3.27185287703775217E+02,     3.31717887196928473E+02, &
       3.36261181979198477E+02,     3.40815058870799018E+02, &
       3.45379407062266854E+02,     3.49954118040770237E+02, &
       3.54539085519440809E+02,     3.59134205369575399E+02/
!             COEFFICIENTS OF ASYMPTOTIC EXPANSION
  DATA CF(1), CF(2), CF(3), CF(4), CF(5), CF(6), CF(7), CF(8), &
       CF(9), CF(10), CF(11), CF(12), CF(13), CF(14), CF(15), &
       CF(16), CF(17), CF(18), CF(19), CF(20), CF(21), CF(22)/ &
       8.33333333333333333E-02,    -2.77777777777777778E-03, &
       7.93650793650793651E-04,    -5.95238095238095238E-04, &
       8.41750841750841751E-04,    -1.91752691752691753E-03, &
       6.41025641025641026E-03,    -2.95506535947712418E-02, &
       1.79644372368830573E-01,    -1.39243221690590112E+00, &
       1.34028640441683920E+01,    -1.56848284626002017E+02, &
       2.19310333333333333E+03,    -3.61087712537249894E+04, &
       6.91472268851313067E+05,    -1.52382215394074162E+07, &
       3.82900751391414141E+08,    -1.08822660357843911E+10, &
       3.47320283765002252E+11,    -1.23696021422692745E+13, &
       4.88788064793079335E+14,    -2.13203339609193739E+16/
!
!             LN(2*PI)
  DATA CON                    /     1.83787706640934548E+00/
!
!***FIRST EXECUTABLE STATEMENT  GAMLN
  IERR=0
  if (Z <= 0.0E0) go to 70
  if (Z > 101.0E0) go to 10
  NZ = Z
  FZ = Z - NZ
  if (FZ > 0.0E0) go to 10
  if (NZ > 100) go to 10
  GAMLN = GLN(NZ)
  return
   10 CONTINUE
  WDTOL = R1MACH(4)
  WDTOL = MAX(WDTOL,0.5E-18)
  I1M = I1MACH(11)
  RLN = R1MACH(5)*I1M
  FLN = MIN(RLN,20.0E0)
  FLN = MAX(FLN,3.0E0)
  FLN = FLN - 3.0E0
  ZM = 1.8000E0 + 0.3875E0*FLN
  MZ = ZM + 1
  ZMIN = MZ
  ZDMY = Z
  ZINC = 0.0E0
  if (Z >= ZMIN) go to 20
  ZINC = ZMIN - NZ
  ZDMY = Z + ZINC
   20 CONTINUE
  ZP = 1.0E0/ZDMY
  T1 = CF(1)*ZP
  S = T1
  if (ZP < WDTOL) go to 40
  ZSQ = ZP*ZP
  TST = T1*WDTOL
  DO 30 K=2,22
    ZP = ZP*ZSQ
    TRM = CF(K)*ZP
    if (ABS(TRM) < TST) go to 40
    S = S + TRM
   30 CONTINUE
   40 CONTINUE
  if (ZINC /= 0.0E0) go to 50
  TLG = ALOG(Z)
  GAMLN = Z*(TLG-1.0E0) + 0.5E0*(CON-TLG) + S
  return
   50 CONTINUE
  ZP = 1.0E0
  NZ = ZINC
  DO 60 I=1,NZ
    ZP = ZP*(Z+(I-1))
   60 CONTINUE
  TLG = ALOG(ZDMY)
  GAMLN = ZDMY*(TLG-1.0E0) - ALOG(ZP) + 0.5E0*(CON-TLG) + S
  return
!
!
   70 CONTINUE
  GAMLN = R1MACH(2)
  IERR=1
  return
end
subroutine CRATI (Z, FNU, N, CY, TOL)
!
!! CRATI is subsidiary to CBESH, CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CRATI-A, ZRATI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
!     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
!     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
!     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
!     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
!     BY D. J. SOOKNE.
!
!***SEE ALSO  CBESH, CBESI, CBESK
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CRATI
  COMPLEX CDFNU, CONE, CY, CZERO, PT, P1, P2, RZ, T1, Z
  REAL AK, AMAGZ, AP1, AP2, ARG, AZ, DFNU, FDNU, FLAM, FNU, FNUP, &
   RAP1, RHO, TEST, TEST1, TOL
  INTEGER I, ID, IDNU, INU, ITIME, K, KK, MAGZ, N
  DIMENSION CY(N)
  DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CRATI
  AZ = ABS(Z)
  INU = FNU
  IDNU = INU + N - 1
  FDNU = IDNU
  MAGZ = AZ
  AMAGZ = MAGZ+1
  FNUP = MAX(AMAGZ,FDNU)
  ID = IDNU - MAGZ - 1
  ITIME = 1
  K = 1
  RZ = (CONE+CONE)/Z
  T1 = CMPLX(FNUP,0.0E0)*RZ
  P2 = -T1
  P1 = CONE
  T1 = T1 + RZ
  if (ID > 0) ID = 0
  AP2 = ABS(P2)
  AP1 = ABS(P1)
!-----------------------------------------------------------------------
!     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE call TO CBKNX
!     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
!     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
!     PREMATURELY.
!-----------------------------------------------------------------------
  ARG = (AP2+AP2)/(AP1*TOL)
  TEST1 = SQRT(ARG)
  TEST = TEST1
  RAP1 = 1.0E0/AP1
  P1 = P1*CMPLX(RAP1,0.0E0)
  P2 = P2*CMPLX(RAP1,0.0E0)
  AP2 = AP2*RAP1
   10 CONTINUE
  K = K + 1
  AP1 = AP2
  PT = P2
  P2 = P1 - T1*P2
  P1 = PT
  T1 = T1 + RZ
  AP2 = ABS(P2)
  if (AP1 <= TEST) go to 10
  if (ITIME == 2) go to 20
  AK = ABS(T1)*0.5E0
  FLAM = AK + SQRT(AK*AK-1.0E0)
  RHO = MIN(AP2/AP1,FLAM)
  TEST = TEST1*SQRT(RHO/(RHO*RHO-1.0E0))
  ITIME = 2
  go to 10
   20 CONTINUE
  KK = K + 1 - ID
  AK = KK
  DFNU = FNU + (N-1)
  CDFNU = CMPLX(DFNU,0.0E0)
  T1 = CMPLX(AK,0.0E0)
  P1 = CMPLX(1.0E0/AP2,0.0E0)
  P2 = CZERO
  DO 30 I=1,KK
    PT = P1
    P1 = RZ*(CDFNU+T1)*P1 + P2
    P2 = PT
    T1 = T1 - CONE
   30 CONTINUE
  if (REAL(P1) /= 0.0E0 .OR. AIMAG(P1) /= 0.0E0) go to 40
  P1 = CMPLX(TOL,TOL)
   40 CONTINUE
  CY(N) = P2/P1
  if (N == 1) RETURN
  K = N - 1
  AK = K
  T1 = CMPLX(AK,0.0E0)
  CDFNU = CMPLX(FNU,0.0E0)*RZ
  DO 60 I=2,N
    PT = CDFNU + T1*RZ + CY(K+1)
    if (REAL(PT) /= 0.0E0 .OR. AIMAG(PT) /= 0.0E0) go to 50
    PT = CMPLX(TOL,TOL)
   50   CONTINUE
    CY(K) = CONE/PT
    T1 = T1 - CONE
    K = K - 1
   60 CONTINUE
  return
end
subroutine CUOIK (Z, FNU, KODE, IKFLG, N, Y, NUF, TOL, ELIM, ALIM)
!
!! CUOIK is subsidiary to CBESH, CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUOIK-A, ZUOIK-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
!     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
!     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
!     WHERE ALIM < ELIM. if THE MAGNITUDE, BASED ON THE LEADING
!     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
!     THE RESULT IS ON SCALE. if NOT, THEN A REFINED TEST USING OTHER
!     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
!     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
!     EXP(-ELIM)/TOL
!
!     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
!          =2 MEANS THE K SEQUENCE IS TESTED
!     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
!         =-1 MEANS AN OVERFLOW WOULD OCCUR
!     IKFLG=1 AND NUF > 0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
!             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
!     IKFLG=2 AND NUF == N MEANS ALL Y VALUES WERE SET TO ZERO
!     IKFLG=2 AND 0 < NUF < N NOT CONSIDERED. Y MUST BE SET BY
!             ANOTHER ROUTINE
!
!***SEE ALSO  CBESH, CBESI, CBESK
!***ROUTINES CALLED  CUCHK, CUNHJ, CUNIK, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUOIK
  COMPLEX ARG, ASUM, BSUM, CWRK, CZ, CZERO, PHI, SUM, Y, Z, ZB, &
   ZETA1, ZETA2, ZN, ZR
  REAL AARG, AIC, ALIM, APHI, ASCLE, AX, AY, ELIM, FNN, FNU, GNN, &
   GNU, RCZ, TOL, X, YY, R1MACH
  INTEGER I, IFORM, IKFLG, INIT, KODE, N, NN, NUF, NW
  DIMENSION Y(N), CWRK(16)
  DATA CZERO / (0.0E0,0.0E0) /
  DATA AIC / 1.265512123484645396E+00 /
!***FIRST EXECUTABLE STATEMENT  CUOIK
  NUF = 0
  NN = N
  X = REAL(Z)
  ZR = Z
  if (X < 0.0E0) ZR = -Z
  ZB = ZR
  YY = AIMAG(ZR)
  AX = ABS(X)*1.7321E0
  AY = ABS(YY)
  IFORM = 1
  if (AY > AX) IFORM = 2
  GNU = MAX(FNU,1.0E0)
  if (IKFLG == 1) go to 10
  FNN = NN
  GNN = FNU + FNN - 1.0E0
  GNU = MAX(GNN,FNN)
   10 CONTINUE
!-----------------------------------------------------------------------
!     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
!     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
!     THE SIGN OF THE IMAGINARY PART CORRECT.
!-----------------------------------------------------------------------
  if (IFORM == 2) go to 20
  INIT = 0
  call CUNIK(ZR, GNU, IKFLG, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM, &
   CWRK)
  CZ = -ZETA1 + ZETA2
  go to 40
   20 CONTINUE
  ZN = -ZR*CMPLX(0.0E0,1.0E0)
  if (YY > 0.0E0) go to 30
  ZN = CONJG(-ZN)
   30 CONTINUE
  call CUNHJ(ZN, GNU, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
  CZ = -ZETA1 + ZETA2
  AARG = ABS(ARG)
   40 CONTINUE
  if (KODE == 2) CZ = CZ - ZB
  if (IKFLG == 2) CZ = -CZ
  APHI = ABS(PHI)
  RCZ = REAL(CZ)
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
  if (RCZ > ELIM) go to 170
  if (RCZ < ALIM) go to 50
  RCZ = RCZ + ALOG(APHI)
  if (IFORM == 2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
  if (RCZ > ELIM) go to 170
  go to 100
   50 CONTINUE
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
  if (RCZ < (-ELIM)) go to 60
  if (RCZ > (-ALIM)) go to 100
  RCZ = RCZ + ALOG(APHI)
  if (IFORM == 2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
  if (RCZ > (-ELIM)) go to 80
   60 CONTINUE
  DO 70 I=1,NN
    Y(I) = CZERO
   70 CONTINUE
  NUF = NN
  return
   80 CONTINUE
  ASCLE = 1.0E+3*R1MACH(1)/TOL
  CZ = CZ + CLOG(PHI)
  if (IFORM == 1) go to 90
  CZ = CZ - CMPLX(0.25E0,0.0E0)*CLOG(ARG) - CMPLX(AIC,0.0E0)
   90 CONTINUE
  AX = EXP(RCZ)/TOL
  AY = AIMAG(CZ)
  CZ = CMPLX(AX,0.0E0)*CMPLX(COS(AY),SIN(AY))
  call CUCHK(CZ, NW, ASCLE, TOL)
  if (NW == 1) go to 60
  100 CONTINUE
  if (IKFLG == 2) RETURN
  if (N == 1) RETURN
!-----------------------------------------------------------------------
!     SET UNDERFLOWS ON I SEQUENCE
!-----------------------------------------------------------------------
  110 CONTINUE
  GNU = FNU + (NN-1)
  if (IFORM == 2) go to 120
  INIT = 0
  call CUNIK(ZR, GNU, IKFLG, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM, &
   CWRK)
  CZ = -ZETA1 + ZETA2
  go to 130
  120 CONTINUE
  call CUNHJ(ZN, GNU, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
  CZ = -ZETA1 + ZETA2
  AARG = ABS(ARG)
  130 CONTINUE
  if (KODE == 2) CZ = CZ - ZB
  APHI = ABS(PHI)
  RCZ = REAL(CZ)
  if (RCZ < (-ELIM)) go to 140
  if (RCZ > (-ALIM)) RETURN
  RCZ = RCZ + ALOG(APHI)
  if (IFORM == 2) RCZ = RCZ - 0.25E0*ALOG(AARG) - AIC
  if (RCZ > (-ELIM)) go to 150
  140 CONTINUE
  Y(NN) = CZERO
  NN = NN - 1
  NUF = NUF + 1
  if (NN == 0) RETURN
  go to 110
  150 CONTINUE
  ASCLE = 1.0E+3*R1MACH(1)/TOL
  CZ = CZ + CLOG(PHI)
  if (IFORM == 1) go to 160
  CZ = CZ - CMPLX(0.25E0,0.0E0)*CLOG(ARG) - CMPLX(AIC,0.0E0)
  160 CONTINUE
  AX = EXP(RCZ)/TOL
  AY = AIMAG(CZ)
  CZ = CMPLX(AX,0.0E0)*CMPLX(COS(AY),SIN(AY))
  call CUCHK(CZ, NW, ASCLE, TOL)
  if (NW == 1) go to 140
  return
  170 CONTINUE
  NUF = -1
  return
end
subroutine CWRSK (ZR, FNU, KODE, N, Y, NZ, CW, TOL, ELIM, ALIM)
!
!! CWRSK is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CWRSK-A, ZWRSK-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY
!     NORMALIZING THE I FUNCTION RATIOS FROM CRATI BY THE WRONSKIAN
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  CBKNU, CRATI, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CWRSK
  COMPLEX CINU, CSCL, CT, CW, C1, C2, RCT, ST, Y, ZR
  REAL ACT, ACW, ALIM, ASCLE, ELIM, FNU, S1, S2, TOL, YY, R1MACH
  INTEGER I, KODE, N, NW, NZ
  DIMENSION Y(N), CW(2)
!***FIRST EXECUTABLE STATEMENT  CWRSK
!
!     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
!     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
!     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
!
  NZ = 0
  call CBKNU(ZR, FNU, KODE, 2, CW, NW, TOL, ELIM, ALIM)
  if (NW /= 0) go to 50
  call CRATI(ZR, FNU, N, Y, TOL)
!-----------------------------------------------------------------------
!     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
!     R(FNU+J-1,Z)=Y(J),  J=1,...,N
!-----------------------------------------------------------------------
  CINU = CMPLX(1.0E0,0.0E0)
  if (KODE == 1) go to 10
  YY = AIMAG(ZR)
  S1 = COS(YY)
  S2 = SIN(YY)
  CINU = CMPLX(S1,S2)
   10 CONTINUE
!-----------------------------------------------------------------------
!     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
!     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
!     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
!     THE RESULT IS ON SCALE.
!-----------------------------------------------------------------------
  ACW = ABS(CW(2))
  ASCLE = 1.0E+3*R1MACH(1)/TOL
  CSCL = CMPLX(1.0E0,0.0E0)
  if (ACW > ASCLE) go to 20
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  go to 30
   20 CONTINUE
  ASCLE = 1.0E0/ASCLE
  if (ACW < ASCLE) go to 30
  CSCL = CMPLX(TOL,0.0E0)
   30 CONTINUE
  C1 = CW(1)*CSCL
  C2 = CW(2)*CSCL
  ST = Y(1)
!-----------------------------------------------------------------------
!     CINU=CINU*(CONJG(CT)/ABS(CT))*(1.0E0/ABS(CT) PREVENTS
!     UNDER- OR OVERFLOW PREMATURELY BY SQUARING ABS(CT)
!-----------------------------------------------------------------------
  CT = ZR*(C2+ST*C1)
  ACT = ABS(CT)
  RCT = CMPLX(1.0E0/ACT,0.0E0)
  CT = CONJG(CT)*RCT
  CINU = CINU*RCT*CT
  Y(1) = CINU*CSCL
  if (N == 1) RETURN
  DO 40 I=2,N
    CINU = ST*CINU
    ST = Y(I)
    Y(I) = CINU*CSCL
   40 CONTINUE
  return
   50 CONTINUE
  NZ = -1
  if ( NW == (-2)) NZ=-2
  return
end
subroutine CUNK2 (Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
!
!! CUNK2 is subsidiary to CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNK2-A, ZUNK2-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
!     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
!     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
!     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
!     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z if Z IS IN THE RIGHT
!     HALF PLANE OR ZR=-Z if Z IS IN THE LEFT HALF PLANE. MR INDIC-
!     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
!     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
!
!***SEE ALSO  CBESK
!***ROUTINES CALLED  CAIRY, CS1S2, CUCHK, CUNHJ, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUNK2
  COMPLEX AI, ARG, ASUM, BSUM, CFN, CI, CIP, &
   CK, CONE, CRSC, CR1, CR2, CS, CSCL, CSGN, CSPN, CSR, CSS, CY, &
   CZERO, C1, C2, DAI, PHI,  RZ, S1, S2, Y, Z, ZB, ZETA1, &
   ZETA2, ZN, ZR, PHID, ARGD, ZETA1D, ZETA2D, ASUMD, BSUMD
  REAL AARG, AIC, ALIM, ANG, APHI, ASC, ASCLE, BRY, CAR, CPN, C2I, &
   C2M, C2R, ELIM, FMR, FN, FNF, FNU, HPI, PI, RS1, SAR, SGN, SPN, &
   TOL, X, YY, R1MACH
  INTEGER I, IB, IFLAG, IFN, IL, IN, INU, IUF, K, KDFLG, KFLAG, KK, &
   KODE, MR, N, NAI, NDAI, NW, NZ, IDUM, J, IPARD, IC
  DIMENSION BRY(3), Y(N), ASUM(2), BSUM(2), PHI(2), ARG(2), &
   ZETA1(2), ZETA2(2), CY(2), CIP(4), CSS(3), CSR(3)
  DATA CZERO, CONE, CI, CR1, CR2 / &
           (0.0E0,0.0E0),(1.0E0,0.0E0),(0.0E0,1.0E0), &
  (1.0E0,1.73205080756887729E0),(-0.5E0,-8.66025403784438647E-01)/
  DATA HPI, PI, AIC / &
       1.57079632679489662E+00,     3.14159265358979324E+00, &
       1.26551212348464539E+00/
  DATA CIP(1),CIP(2),CIP(3),CIP(4)/ &
   (1.0E0,0.0E0), (0.0E0,-1.0E0), (-1.0E0,0.0E0), (0.0E0,1.0E0)/
!***FIRST EXECUTABLE STATEMENT  CUNK2
  KDFLG = 1
  NZ = 0
!-----------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
!     THE UNDERFLOW LIMIT
!-----------------------------------------------------------------------
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CRSC = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CRSC
  CSR(1) = CRSC
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = 1.0E+3*R1MACH(1)/TOL
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = R1MACH(2)
  X = REAL(Z)
  ZR = Z
  if (X < 0.0E0) ZR = -Z
  YY = AIMAG(ZR)
  ZN = -ZR*CI
  ZB = ZR
  INU = FNU
  FNF = FNU - INU
  ANG = -HPI*FNF
  CAR = COS(ANG)
  SAR = SIN(ANG)
  CPN = -HPI*CAR
  SPN = -HPI*SAR
  C2 = CMPLX(-SPN,CPN)
  KK = MOD(INU,4) + 1
  CS = CR1*C2*CIP(KK)
  if (YY > 0.0E0) go to 10
  ZN = CONJG(-ZN)
  ZB = CONJG(ZB)
   10 CONTINUE
!-----------------------------------------------------------------------
!     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
!     QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE COMPUTED BY
!     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
!-----------------------------------------------------------------------
  J = 2
  DO 70 I=1,N
!-----------------------------------------------------------------------
!     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
!-----------------------------------------------------------------------
    J = 3 - J
    FN = FNU + (I-1)
    call CUNHJ(ZN, FN, 0, TOL, PHI(J), ARG(J), ZETA1(J), ZETA2(J), &
     ASUM(J), BSUM(J))
    if (KODE == 1) go to 20
    CFN = CMPLX(FN,0.0E0)
    S1 = ZETA1(J) - CFN*(CFN/(ZB+ZETA2(J)))
    go to 30
   20   CONTINUE
    S1 = ZETA1(J) - ZETA2(J)
   30   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = REAL(S1)
    if (ABS(RS1) > ELIM) go to 60
    if (KDFLG == 1) KFLAG = 2
    if (ABS(RS1) < ALIM) go to 40
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    APHI = ABS(PHI(J))
    AARG = ABS(ARG(J))
    RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
    if (ABS(RS1) > ELIM) go to 60
    if (KDFLG == 1) KFLAG = 1
    if (RS1 < 0.0E0) go to 40
    if (KDFLG == 1) KFLAG = 3
   40   CONTINUE
!-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!-----------------------------------------------------------------------
    C2 = ARG(J)*CR2
    call CAIRY(C2, 0, 2, AI, NAI, IDUM)
    call CAIRY(C2, 1, 2, DAI, NDAI, IDUM)
    S2 = CS*PHI(J)*(AI*ASUM(J)+CR2*DAI*BSUM(J))
    C2R = REAL(S1)
    C2I = AIMAG(S1)
    C2M = EXP(C2R)*REAL(CSS(KFLAG))
    S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
    S2 = S2*S1
    if (KFLAG /= 1) go to 50
    call CUCHK(S2, NW, BRY(1), TOL)
    if (NW /= 0) go to 60
   50   CONTINUE
    if (YY <= 0.0E0) S2 = CONJG(S2)
    CY(KDFLG) = S2
    Y(I) = S2*CSR(KFLAG)
    CS = -CI*CS
    if (KDFLG == 2) go to 75
    KDFLG = 2
    go to 70
   60   CONTINUE
    if (RS1 > 0.0E0) go to 300
!-----------------------------------------------------------------------
!     FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
    if (X < 0.0E0) go to 300
    KDFLG = 1
    Y(I) = CZERO
    CS = -CI*CS
    NZ=NZ+1
    if (I == 1) go to 70
    if (Y(I-1) == CZERO) go to 70
    Y(I-1) = CZERO
    NZ=NZ+1
   70 CONTINUE
  I=N
   75 CONTINUE
  RZ = CMPLX(2.0E0,0.0E0)/ZR
  CK = CMPLX(FN,0.0E0)*RZ
  IB = I + 1
  if (N < IB) go to 170
!-----------------------------------------------------------------------
!     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
!     ON UNDERFLOW
!-----------------------------------------------------------------------
  FN = FNU+(N-1)
  IPARD = 1
  if (MR /= 0) IPARD = 0
  call CUNHJ(ZN,FN,IPARD,TOL,PHID,ARGD,ZETA1D,ZETA2D,ASUMD,BSUMD)
  if (KODE == 1) go to 80
  CFN=CMPLX(FN,0.0E0)
  S1=ZETA1D-CFN*(CFN/(ZB+ZETA2D))
  go to 90
   80 CONTINUE
  S1=ZETA1D-ZETA2D
   90 CONTINUE
  RS1=REAL(S1)
  if (ABS(RS1) > ELIM) go to 95
  if (ABS(RS1) < ALIM) go to 100
!-----------------------------------------------------------------------
!     REFINE ESTIMATE AND TEST
!-----------------------------------------------------------------------
  APHI=ABS(PHID)
  AARG = ABS(ARGD)
  RS1=RS1+ALOG(APHI)-0.25E0*ALOG(AARG)-AIC
  if (ABS(RS1) < ELIM) go to 100
   95 CONTINUE
  if (RS1 > 0.0E0) go to 300
!-----------------------------------------------------------------------
!     FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
  if (X < 0.0E0) go to 300
  NZ=N
  DO 96 I=1,N
    Y(I) = CZERO
   96 CONTINUE
  return
  100 CONTINUE
!-----------------------------------------------------------------------
!     SCALED FORWARD RECURRENCE FOR REMAINDER OF THE SEQUENCE
!-----------------------------------------------------------------------
  S1 = CY(1)
  S2 = CY(2)
  C1 = CSR(KFLAG)
  ASCLE = BRY(KFLAG)
  DO 120 I=IB,N
    C2 = S2
    S2 = CK*S2 + S1
    S1 = C2
    CK = CK + RZ
    C2 = S2*C1
    Y(I) = C2
    if (KFLAG >= 3) go to 120
    C2R = REAL(C2)
    C2I = AIMAG(C2)
    C2R = ABS(C2R)
    C2I = ABS(C2I)
    C2M = MAX(C2R,C2I)
    if (C2M <= ASCLE) go to 120
    KFLAG = KFLAG + 1
    ASCLE = BRY(KFLAG)
    S1 = S1*C1
    S2 = C2
    S1 = S1*CSS(KFLAG)
    S2 = S2*CSS(KFLAG)
    C1 = CSR(KFLAG)
  120 CONTINUE
  170 CONTINUE
  if (MR == 0) RETURN
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION FOR RE(Z) < 0.0E0
!-----------------------------------------------------------------------
  NZ = 0
  FMR = MR
  SGN = -SIGN(PI,FMR)
!-----------------------------------------------------------------------
!     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
!-----------------------------------------------------------------------
  CSGN = CMPLX(0.0E0,SGN)
  if (YY <= 0.0E0) CSGN = CONJG(CSGN)
  IFN = INU + N - 1
  ANG = FNF*SGN
  CPN = COS(ANG)
  SPN = SIN(ANG)
  CSPN = CMPLX(CPN,SPN)
  if (MOD(IFN,2) == 1) CSPN = -CSPN
!-----------------------------------------------------------------------
!     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
!     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
!     QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE COMPUTED BY
!     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
!-----------------------------------------------------------------------
  CS = CMPLX(CAR,-SAR)*CSGN
  IN = MOD(IFN,4) + 1
  C2 = CIP(IN)
  CS = CS*CONJG(C2)
  ASC = BRY(1)
  KK = N
  KDFLG = 1
  IB = IB-1
  IC = IB-1
  IUF = 0
  DO 270 K=1,N
!-----------------------------------------------------------------------
!     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
!     FUNCTION ABOVE
!-----------------------------------------------------------------------
    FN = FNU+(KK-1)
    if (N > 2) go to 180
  175   CONTINUE
    PHID = PHI(J)
    ARGD = ARG(J)
    ZETA1D = ZETA1(J)
    ZETA2D = ZETA2(J)
    ASUMD = ASUM(J)
    BSUMD = BSUM(J)
    J = 3 - J
    go to 190
  180   CONTINUE
    if ((KK == N).AND.(IB < N)) go to 190
    if ((KK == IB).OR.(KK == IC)) go to 175
    call CUNHJ(ZN, FN, 0, TOL, PHID, ARGD, ZETA1D, ZETA2D, &
     ASUMD, BSUMD)
  190   CONTINUE
    if (KODE == 1) go to 200
    CFN = CMPLX(FN,0.0E0)
    S1 = -ZETA1D + CFN*(CFN/(ZB+ZETA2D))
    go to 210
  200   CONTINUE
    S1 = -ZETA1D + ZETA2D
  210   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = REAL(S1)
    if (ABS(RS1) > ELIM) go to 260
    if (KDFLG == 1) IFLAG = 2
    if (ABS(RS1) < ALIM) go to 220
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    APHI = ABS(PHID)
    AARG = ABS(ARGD)
    RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
    if (ABS(RS1) > ELIM) go to 260
    if (KDFLG == 1) IFLAG = 1
    if (RS1 < 0.0E0) go to 220
    if (KDFLG == 1) IFLAG = 3
  220   CONTINUE
    call CAIRY(ARGD, 0, 2, AI, NAI, IDUM)
    call CAIRY(ARGD, 1, 2, DAI, NDAI, IDUM)
    S2 = CS*PHID*(AI*ASUMD+DAI*BSUMD)
    C2R = REAL(S1)
    C2I = AIMAG(S1)
    C2M = EXP(C2R)*REAL(CSS(IFLAG))
    S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
    S2 = S2*S1
    if (IFLAG /= 1) go to 230
    call CUCHK(S2, NW, BRY(1), TOL)
    if (NW /= 0) S2 = CMPLX(0.0E0,0.0E0)
  230   CONTINUE
    if (YY <= 0.0E0) S2 = CONJG(S2)
    CY(KDFLG) = S2
    C2 = S2
    S2 = S2*CSR(IFLAG)
!-----------------------------------------------------------------------
!     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
!-----------------------------------------------------------------------
    S1 = Y(KK)
    if (KODE == 1) go to 250
    call CS1S2(ZR, S1, S2, NW, ASC, ALIM, IUF)
    NZ = NZ + NW
  250   CONTINUE
    Y(KK) = S1*CSPN + S2
    KK = KK - 1
    CSPN = -CSPN
    CS = -CS*CI
    if (C2 /= CZERO) go to 255
    KDFLG = 1
    go to 270
  255   CONTINUE
    if (KDFLG == 2) go to 275
    KDFLG = 2
    go to 270
  260   CONTINUE
    if (RS1 > 0.0E0) go to 300
    S2 = CZERO
    go to 230
  270 CONTINUE
  K = N
  275 CONTINUE
  IL = N-K
  if (IL == 0) RETURN
!-----------------------------------------------------------------------
!     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
!     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
!     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
!-----------------------------------------------------------------------
  S1 = CY(1)
  S2 = CY(2)
  CS = CSR(IFLAG)
  ASCLE = BRY(IFLAG)
  FN = INU+IL
  DO 290 I=1,IL
    C2 = S2
    S2 = S1 + CMPLX(FN+FNF,0.0E0)*RZ*S2
    S1 = C2
    FN = FN - 1.0E0
    C2 = S2*CS
    CK = C2
    C1 = Y(KK)
    if (KODE == 1) go to 280
    call CS1S2(ZR, C1, C2, NW, ASC, ALIM, IUF)
    NZ = NZ + NW
  280   CONTINUE
    Y(KK) = C1*CSPN + C2
    KK = KK - 1
    CSPN = -CSPN
    if (IFLAG >= 3) go to 290
    C2R = REAL(CK)
    C2I = AIMAG(CK)
    C2R = ABS(C2R)
    C2I = ABS(C2I)
    C2M = MAX(C2R,C2I)
    if (C2M <= ASCLE) go to 290
    IFLAG = IFLAG + 1
    ASCLE = BRY(IFLAG)
    S1 = S1*CS
    S2 = CK
    S1 = S1*CSS(IFLAG)
    S2 = S2*CSS(IFLAG)
    CS = CSR(IFLAG)
  290 CONTINUE
  return
  300 CONTINUE
  NZ = -1
  return
end
subroutine CUNIK (ZR, FNU, IKFLG, IPMTR, TOL, INIT, PHI, ZETA1, &
     ZETA2, SUM, CWRK)
!
!! CUNIK is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNIK-A, ZUNIK-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!        CUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
!        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
!        RESPECTIVELY BY
!
!        W(FNU,ZR) = PHI*EXP(ZETA)*SUM
!
!        WHERE       ZETA=-ZETA1 + ZETA2       OR
!                          ZETA1 - ZETA2
!
!        THE FIRST call MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
!        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG=
!        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
!        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
!        ZETA1,ZETA2.
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUNIK
  COMPLEX CFN, CON, CONE, CRFN, CWRK, CZERO, PHI, S, SR, SUM, T, &
   T2, ZETA1, ZETA2, ZN, ZR
  REAL AC, C, FNU, RFN, TEST, TOL, TSTR, TSTI, R1MACH
  INTEGER I, IKFLG, INIT, IPMTR, J, K, L
  DIMENSION C(120), CWRK(16), CON(2)
  DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
  DATA CON(1), CON(2)  / &
  (3.98942280401432678E-01,0.0E0),(1.25331413731550025E+00,0.0E0)/
  DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10), &
       C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18), &
       C(19), C(20), C(21), C(22), C(23), C(24)/ &
       1.00000000000000000E+00,    -2.08333333333333333E-01, &
       1.25000000000000000E-01,     3.34201388888888889E-01, &
      -4.01041666666666667E-01,     7.03125000000000000E-02, &
      -1.02581259645061728E+00,     1.84646267361111111E+00, &
      -8.91210937500000000E-01,     7.32421875000000000E-02, &
       4.66958442342624743E+00,    -1.12070026162229938E+01, &
       8.78912353515625000E+00,    -2.36408691406250000E+00, &
       1.12152099609375000E-01,    -2.82120725582002449E+01, &
       8.46362176746007346E+01,    -9.18182415432400174E+01, &
       4.25349987453884549E+01,    -7.36879435947963170E+00, &
       2.27108001708984375E-01,     2.12570130039217123E+02, &
      -7.65252468141181642E+02,     1.05999045252799988E+03/
  DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32), &
       C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40), &
       C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/ &
      -6.99579627376132541E+02,     2.18190511744211590E+02, &
      -2.64914304869515555E+01,     5.72501420974731445E-01, &
      -1.91945766231840700E+03,     8.06172218173730938E+03, &
      -1.35865500064341374E+04,     1.16553933368645332E+04, &
      -5.30564697861340311E+03,     1.20090291321635246E+03, &
      -1.08090919788394656E+02,     1.72772750258445740E+00, &
       2.02042913309661486E+04,    -9.69805983886375135E+04, &
       1.92547001232531532E+05,    -2.03400177280415534E+05, &
       1.22200464983017460E+05,    -4.11926549688975513E+04, &
       7.10951430248936372E+03,    -4.93915304773088012E+02, &
       6.07404200127348304E+00,    -2.42919187900551333E+05, &
       1.31176361466297720E+06,    -2.99801591853810675E+06/
  DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56), &
       C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64), &
       C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/ &
       3.76327129765640400E+06,    -2.81356322658653411E+06, &
       1.26836527332162478E+06,    -3.31645172484563578E+05, &
       4.52187689813627263E+04,    -2.49983048181120962E+03, &
       2.43805296995560639E+01,     3.28446985307203782E+06, &
      -1.97068191184322269E+07,     5.09526024926646422E+07, &
      -7.41051482115326577E+07,     6.63445122747290267E+07, &
      -3.75671766607633513E+07,     1.32887671664218183E+07, &
      -2.78561812808645469E+06,     3.08186404612662398E+05, &
      -1.38860897537170405E+04,     1.10017140269246738E+02, &
      -4.93292536645099620E+07,     3.25573074185765749E+08, &
      -9.39462359681578403E+08,     1.55359689957058006E+09, &
      -1.62108055210833708E+09,     1.10684281682301447E+09/
  DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80), &
       C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88), &
       C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/ &
      -4.95889784275030309E+08,     1.42062907797533095E+08, &
      -2.44740627257387285E+07,     2.24376817792244943E+06, &
      -8.40054336030240853E+04,     5.51335896122020586E+02, &
       8.14789096118312115E+08,    -5.86648149205184723E+09, &
       1.86882075092958249E+10,    -3.46320433881587779E+10, &
       4.12801855797539740E+10,    -3.30265997498007231E+10, &
       1.79542137311556001E+10,    -6.56329379261928433E+09, &
       1.55927986487925751E+09,    -2.25105661889415278E+08, &
       1.73951075539781645E+07,    -5.49842327572288687E+05, &
       3.03809051092238427E+03,    -1.46792612476956167E+10, &
       1.14498237732025810E+11,    -3.99096175224466498E+11, &
       8.19218669548577329E+11,    -1.09837515608122331E+12/
  DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104), &
       C(105), C(106), C(107), C(108), C(109), C(110), C(111), &
       C(112), C(113), C(114), C(115), C(116), C(117), C(118)/ &
       1.00815810686538209E+12,    -6.45364869245376503E+11, &
       2.87900649906150589E+11,    -8.78670721780232657E+10, &
       1.76347306068349694E+10,    -2.16716498322379509E+09, &
       1.43157876718888981E+08,    -3.87183344257261262E+06, &
       1.82577554742931747E+04,     2.86464035717679043E+11, &
      -2.40629790002850396E+12,     9.10934118523989896E+12, &
      -2.05168994109344374E+13,     3.05651255199353206E+13, &
      -3.16670885847851584E+13,     2.33483640445818409E+13, &
      -1.23204913055982872E+13,     4.61272578084913197E+12, &
      -1.19655288019618160E+12,     2.05914503232410016E+11, &
      -2.18229277575292237E+10,     1.24700929351271032E+09/
  DATA C(119), C(120)/ &
      -2.91883881222208134E+07,     1.18838426256783253E+05/
!***FIRST EXECUTABLE STATEMENT  CUNIK
  if (INIT /= 0) go to 40
!-----------------------------------------------------------------------
!     INITIALIZE ALL VARIABLES
!-----------------------------------------------------------------------
  RFN = 1.0E0/FNU
  CRFN = CMPLX(RFN,0.0E0)
!     T = ZR*CRFN
!-----------------------------------------------------------------------
!     OVERFLOW TEST (ZR/FNU TOO SMALL)
!-----------------------------------------------------------------------
  TSTR = REAL(ZR)
  TSTI = AIMAG(ZR)
  TEST = R1MACH(1)*1.0E+3
  AC = FNU*TEST
  if (ABS(TSTR) > AC .OR. ABS(TSTI) > AC) go to 15
  AC = 2.0E0*ABS(ALOG(TEST))+FNU
  ZETA1 = CMPLX(AC,0.0E0)
  ZETA2 = CMPLX(FNU,0.0E0)
  PHI=CONE
  return
   15 CONTINUE
  T=ZR*CRFN
  S = CONE + T*T
  SR = CSQRT(S)
  CFN = CMPLX(FNU,0.0E0)
  ZN = (CONE+SR)/T
  ZETA1 = CFN*CLOG(ZN)
  ZETA2 = CFN*SR
  T = CONE/SR
  SR = T*CRFN
  CWRK(16) = CSQRT(SR)
  PHI = CWRK(16)*CON(IKFLG)
  if (IPMTR /= 0) RETURN
  T2 = CONE/S
  CWRK(1) = CONE
  CRFN = CONE
  AC = 1.0E0
  L = 1
  DO 20 K=2,15
    S = CZERO
    DO 10 J=1,K
      L = L + 1
      S = S*T2 + CMPLX(C(L),0.0E0)
   10   CONTINUE
    CRFN = CRFN*SR
    CWRK(K) = CRFN*S
    AC = AC*RFN
    TSTR = REAL(CWRK(K))
    TSTI = AIMAG(CWRK(K))
    TEST = ABS(TSTR) + ABS(TSTI)
    if (AC < TOL .AND. TEST < TOL) go to 30
   20 CONTINUE
  K = 15
   30 CONTINUE
  INIT = K
   40 CONTINUE
  if (IKFLG == 2) go to 60
!-----------------------------------------------------------------------
!     COMPUTE SUM FOR THE I FUNCTION
!-----------------------------------------------------------------------
  S = CZERO
  DO 50 I=1,INIT
    S = S + CWRK(I)
   50 CONTINUE
  SUM = S
  PHI = CWRK(16)*CON(1)
  return
   60 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE SUM FOR THE K FUNCTION
!-----------------------------------------------------------------------
  S = CZERO
  T = CONE
  DO 70 I=1,INIT
    S = S + T*CWRK(I)
    T = -T
   70 CONTINUE
  SUM = S
  PHI = CWRK(16)*CON(2)
  return
end
subroutine CUNHJ (Z, FNU, IPMTR, TOL, PHI, ARG, ZETA1, ZETA2, &
     ASUM, BSUM)
!
!! CUNHJ is subsidiary to CBESI and CBESK
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNHJ-A, ZUNHJ-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     REFERENCES
!         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
!         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
!
!         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
!         PRESS, N.Y., 1974, PAGE 420
!
!     ABSTRACT
!         CUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
!         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
!         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
!
!         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
!
!         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
!         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
!
!               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
!
!         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
!         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
!
!         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
!         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
!         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUNHJ
  COMPLEX ARG, ASUM, BSUM, CFNU, CONE, CR, CZERO, DR, P, PHI, &
   PRZTH, PTFN, RFN13, RTZTA, RZTH, SUMA, SUMB, TFN, T2, UP, W, W2, &
   Z, ZA, ZB, ZC, ZETA, ZETA1, ZETA2, ZTH
  REAL ALFA, ANG, AP, AR, ATOL, AW2, AZTH, BETA, BR, BTOL, C, EX1, &
   EX2, FNU, FN13, FN23, GAMA, HPI, PI, PP, RFNU, RFNU2, THPI, TOL, &
   WI, WR, ZCI, ZCR, ZETAI, ZETAR, ZTHI, ZTHR, ASUMR, ASUMI, BSUMR, &
   BSUMI, TEST, TSTR, TSTI, AC, R1MACH
  INTEGER IAS, IBS, IPMTR, IS, J, JR, JU, K, KMAX, KP1, KS, L, LR, &
   LRP1, L1, L2, M
  DIMENSION AR(14), BR(14), C(105), ALFA(180), BETA(210), GAMA(30), &
   AP(30), P(30), UP(14), CR(14), DR(14)
  DATA AR(1), AR(2), AR(3), AR(4), AR(5), AR(6), AR(7), AR(8), &
       AR(9), AR(10), AR(11), AR(12), AR(13), AR(14)/ &
       1.00000000000000000E+00,     1.04166666666666667E-01, &
       8.35503472222222222E-02,     1.28226574556327160E-01, &
       2.91849026464140464E-01,     8.81627267443757652E-01, &
       3.32140828186276754E+00,     1.49957629868625547E+01, &
       7.89230130115865181E+01,     4.74451538868264323E+02, &
       3.20749009089066193E+03,     2.40865496408740049E+04, &
       1.98923119169509794E+05,     1.79190200777534383E+06/
  DATA BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8), &
       BR(9), BR(10), BR(11), BR(12), BR(13), BR(14)/ &
       1.00000000000000000E+00,    -1.45833333333333333E-01, &
      -9.87413194444444444E-02,    -1.43312053915895062E-01, &
      -3.17227202678413548E-01,    -9.42429147957120249E-01, &
      -3.51120304082635426E+00,    -1.57272636203680451E+01, &
      -8.22814390971859444E+01,    -4.92355370523670524E+02, &
      -3.31621856854797251E+03,    -2.48276742452085896E+04, &
      -2.04526587315129788E+05,    -1.83844491706820990E+06/
  DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10), &
       C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18), &
       C(19), C(20), C(21), C(22), C(23), C(24)/ &
       1.00000000000000000E+00,    -2.08333333333333333E-01, &
       1.25000000000000000E-01,     3.34201388888888889E-01, &
      -4.01041666666666667E-01,     7.03125000000000000E-02, &
      -1.02581259645061728E+00,     1.84646267361111111E+00, &
      -8.91210937500000000E-01,     7.32421875000000000E-02, &
       4.66958442342624743E+00,    -1.12070026162229938E+01, &
       8.78912353515625000E+00,    -2.36408691406250000E+00, &
       1.12152099609375000E-01,    -2.82120725582002449E+01, &
       8.46362176746007346E+01,    -9.18182415432400174E+01, &
       4.25349987453884549E+01,    -7.36879435947963170E+00, &
       2.27108001708984375E-01,     2.12570130039217123E+02, &
      -7.65252468141181642E+02,     1.05999045252799988E+03/
  DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32), &
       C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40), &
       C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/ &
      -6.99579627376132541E+02,     2.18190511744211590E+02, &
      -2.64914304869515555E+01,     5.72501420974731445E-01, &
      -1.91945766231840700E+03,     8.06172218173730938E+03, &
      -1.35865500064341374E+04,     1.16553933368645332E+04, &
      -5.30564697861340311E+03,     1.20090291321635246E+03, &
      -1.08090919788394656E+02,     1.72772750258445740E+00, &
       2.02042913309661486E+04,    -9.69805983886375135E+04, &
       1.92547001232531532E+05,    -2.03400177280415534E+05, &
       1.22200464983017460E+05,    -4.11926549688975513E+04, &
       7.10951430248936372E+03,    -4.93915304773088012E+02, &
       6.07404200127348304E+00,    -2.42919187900551333E+05, &
       1.31176361466297720E+06,    -2.99801591853810675E+06/
  DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56), &
       C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64), &
       C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/ &
       3.76327129765640400E+06,    -2.81356322658653411E+06, &
       1.26836527332162478E+06,    -3.31645172484563578E+05, &
       4.52187689813627263E+04,    -2.49983048181120962E+03, &
       2.43805296995560639E+01,     3.28446985307203782E+06, &
      -1.97068191184322269E+07,     5.09526024926646422E+07, &
      -7.41051482115326577E+07,     6.63445122747290267E+07, &
      -3.75671766607633513E+07,     1.32887671664218183E+07, &
      -2.78561812808645469E+06,     3.08186404612662398E+05, &
      -1.38860897537170405E+04,     1.10017140269246738E+02, &
      -4.93292536645099620E+07,     3.25573074185765749E+08, &
      -9.39462359681578403E+08,     1.55359689957058006E+09, &
      -1.62108055210833708E+09,     1.10684281682301447E+09/
  DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80), &
       C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88), &
       C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/ &
      -4.95889784275030309E+08,     1.42062907797533095E+08, &
      -2.44740627257387285E+07,     2.24376817792244943E+06, &
      -8.40054336030240853E+04,     5.51335896122020586E+02, &
       8.14789096118312115E+08,    -5.86648149205184723E+09, &
       1.86882075092958249E+10,    -3.46320433881587779E+10, &
       4.12801855797539740E+10,    -3.30265997498007231E+10, &
       1.79542137311556001E+10,    -6.56329379261928433E+09, &
       1.55927986487925751E+09,    -2.25105661889415278E+08, &
       1.73951075539781645E+07,    -5.49842327572288687E+05, &
       3.03809051092238427E+03,    -1.46792612476956167E+10, &
       1.14498237732025810E+11,    -3.99096175224466498E+11, &
       8.19218669548577329E+11,    -1.09837515608122331E+12/
  DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104), &
       C(105)/ &
       1.00815810686538209E+12,    -6.45364869245376503E+11, &
       2.87900649906150589E+11,    -8.78670721780232657E+10, &
       1.76347306068349694E+10,    -2.16716498322379509E+09, &
       1.43157876718888981E+08,    -3.87183344257261262E+06, &
       1.82577554742931747E+04/
  DATA ALFA(1), ALFA(2), ALFA(3), ALFA(4), ALFA(5), ALFA(6), &
       ALFA(7), ALFA(8), ALFA(9), ALFA(10), ALFA(11), ALFA(12), &
       ALFA(13), ALFA(14), ALFA(15), ALFA(16), ALFA(17), ALFA(18), &
       ALFA(19), ALFA(20), ALFA(21), ALFA(22)/ &
      -4.44444444444444444E-03,    -9.22077922077922078E-04, &
      -8.84892884892884893E-05,     1.65927687832449737E-04, &
       2.46691372741792910E-04,     2.65995589346254780E-04, &
       2.61824297061500945E-04,     2.48730437344655609E-04, &
       2.32721040083232098E-04,     2.16362485712365082E-04, &
       2.00738858762752355E-04,     1.86267636637545172E-04, &
       1.73060775917876493E-04,     1.61091705929015752E-04, &
       1.50274774160908134E-04,     1.40503497391269794E-04, &
       1.31668816545922806E-04,     1.23667445598253261E-04, &
       1.16405271474737902E-04,     1.09798298372713369E-04, &
       1.03772410422992823E-04,     9.82626078369363448E-05/
  DATA ALFA(23), ALFA(24), ALFA(25), ALFA(26), ALFA(27), ALFA(28), &
       ALFA(29), ALFA(30), ALFA(31), ALFA(32), ALFA(33), ALFA(34), &
       ALFA(35), ALFA(36), ALFA(37), ALFA(38), ALFA(39), ALFA(40), &
       ALFA(41), ALFA(42), ALFA(43), ALFA(44)/ &
       9.32120517249503256E-05,     8.85710852478711718E-05, &
       8.42963105715700223E-05,     8.03497548407791151E-05, &
       7.66981345359207388E-05,     7.33122157481777809E-05, &
       7.01662625163141333E-05,     6.72375633790160292E-05, &
       6.93735541354588974E-04,     2.32241745182921654E-04, &
      -1.41986273556691197E-05,    -1.16444931672048640E-04, &
      -1.50803558053048762E-04,    -1.55121924918096223E-04, &
      -1.46809756646465549E-04,    -1.33815503867491367E-04, &
      -1.19744975684254051E-04,    -1.06184319207974020E-04, &
      -9.37699549891194492E-05,    -8.26923045588193274E-05, &
      -7.29374348155221211E-05,    -6.44042357721016283E-05/
  DATA ALFA(45), ALFA(46), ALFA(47), ALFA(48), ALFA(49), ALFA(50), &
       ALFA(51), ALFA(52), ALFA(53), ALFA(54), ALFA(55), ALFA(56), &
       ALFA(57), ALFA(58), ALFA(59), ALFA(60), ALFA(61), ALFA(62), &
       ALFA(63), ALFA(64), ALFA(65), ALFA(66)/ &
      -5.69611566009369048E-05,    -5.04731044303561628E-05, &
      -4.48134868008882786E-05,    -3.98688727717598864E-05, &
      -3.55400532972042498E-05,    -3.17414256609022480E-05, &
      -2.83996793904174811E-05,    -2.54522720634870566E-05, &
      -2.28459297164724555E-05,    -2.05352753106480604E-05, &
      -1.84816217627666085E-05,    -1.66519330021393806E-05, &
      -1.50179412980119482E-05,    -1.35554031379040526E-05, &
      -1.22434746473858131E-05,    -1.10641884811308169E-05, &
      -3.54211971457743841E-04,    -1.56161263945159416E-04, &
       3.04465503594936410E-05,     1.30198655773242693E-04, &
       1.67471106699712269E-04,     1.70222587683592569E-04/
  DATA ALFA(67), ALFA(68), ALFA(69), ALFA(70), ALFA(71), ALFA(72), &
       ALFA(73), ALFA(74), ALFA(75), ALFA(76), ALFA(77), ALFA(78), &
       ALFA(79), ALFA(80), ALFA(81), ALFA(82), ALFA(83), ALFA(84), &
       ALFA(85), ALFA(86), ALFA(87), ALFA(88)/ &
       1.56501427608594704E-04,     1.36339170977445120E-04, &
       1.14886692029825128E-04,     9.45869093034688111E-05, &
       7.64498419250898258E-05,     6.07570334965197354E-05, &
       4.74394299290508799E-05,     3.62757512005344297E-05, &
       2.69939714979224901E-05,     1.93210938247939253E-05, &
       1.30056674793963203E-05,     7.82620866744496661E-06, &
       3.59257485819351583E-06,     1.44040049814251817E-07, &
      -2.65396769697939116E-06,    -4.91346867098485910E-06, &
      -6.72739296091248287E-06,    -8.17269379678657923E-06, &
      -9.31304715093561232E-06,    -1.02011418798016441E-05, &
      -1.08805962510592880E-05,    -1.13875481509603555E-05/
  DATA ALFA(89), ALFA(90), ALFA(91), ALFA(92), ALFA(93), ALFA(94), &
       ALFA(95), ALFA(96), ALFA(97), ALFA(98), ALFA(99), ALFA(100), &
       ALFA(101), ALFA(102), ALFA(103), ALFA(104), ALFA(105), &
       ALFA(106), ALFA(107), ALFA(108), ALFA(109), ALFA(110)/ &
      -1.17519675674556414E-05,    -1.19987364870944141E-05, &
       3.78194199201772914E-04,     2.02471952761816167E-04, &
      -6.37938506318862408E-05,    -2.38598230603005903E-04, &
      -3.10916256027361568E-04,    -3.13680115247576316E-04, &
      -2.78950273791323387E-04,    -2.28564082619141374E-04, &
      -1.75245280340846749E-04,    -1.25544063060690348E-04, &
      -8.22982872820208365E-05,    -4.62860730588116458E-05, &
      -1.72334302366962267E-05,     5.60690482304602267E-06, &
       2.31395443148286800E-05,     3.62642745856793957E-05, &
       4.58006124490188752E-05,     5.24595294959114050E-05, &
       5.68396208545815266E-05,     5.94349820393104052E-05/
  DATA ALFA(111), ALFA(112), ALFA(113), ALFA(114), ALFA(115), &
       ALFA(116), ALFA(117), ALFA(118), ALFA(119), ALFA(120), &
       ALFA(121), ALFA(122), ALFA(123), ALFA(124), ALFA(125), &
       ALFA(126), ALFA(127), ALFA(128), ALFA(129), ALFA(130)/ &
       6.06478527578421742E-05,     6.08023907788436497E-05, &
       6.01577894539460388E-05,     5.89199657344698500E-05, &
       5.72515823777593053E-05,     5.52804375585852577E-05, &
       5.31063773802880170E-05,     5.08069302012325706E-05, &
       4.84418647620094842E-05,     4.60568581607475370E-05, &
      -6.91141397288294174E-04,    -4.29976633058871912E-04, &
       1.83067735980039018E-04,     6.60088147542014144E-04, &
       8.75964969951185931E-04,     8.77335235958235514E-04, &
       7.49369585378990637E-04,     5.63832329756980918E-04, &
       3.68059319971443156E-04,     1.88464535514455599E-04/
  DATA ALFA(131), ALFA(132), ALFA(133), ALFA(134), ALFA(135), &
       ALFA(136), ALFA(137), ALFA(138), ALFA(139), ALFA(140), &
       ALFA(141), ALFA(142), ALFA(143), ALFA(144), ALFA(145), &
       ALFA(146), ALFA(147), ALFA(148), ALFA(149), ALFA(150)/ &
       3.70663057664904149E-05,    -8.28520220232137023E-05, &
      -1.72751952869172998E-04,    -2.36314873605872983E-04, &
      -2.77966150694906658E-04,    -3.02079514155456919E-04, &
      -3.12594712643820127E-04,    -3.12872558758067163E-04, &
      -3.05678038466324377E-04,    -2.93226470614557331E-04, &
      -2.77255655582934777E-04,    -2.59103928467031709E-04, &
      -2.39784014396480342E-04,    -2.20048260045422848E-04, &
      -2.00443911094971498E-04,    -1.81358692210970687E-04, &
      -1.63057674478657464E-04,    -1.45712672175205844E-04, &
      -1.29425421983924587E-04,    -1.14245691942445952E-04/
  DATA ALFA(151), ALFA(152), ALFA(153), ALFA(154), ALFA(155), &
       ALFA(156), ALFA(157), ALFA(158), ALFA(159), ALFA(160), &
       ALFA(161), ALFA(162), ALFA(163), ALFA(164), ALFA(165), &
       ALFA(166), ALFA(167), ALFA(168), ALFA(169), ALFA(170)/ &
       1.92821964248775885E-03,     1.35592576302022234E-03, &
      -7.17858090421302995E-04,    -2.58084802575270346E-03, &
      -3.49271130826168475E-03,    -3.46986299340960628E-03, &
      -2.82285233351310182E-03,    -1.88103076404891354E-03, &
      -8.89531718383947600E-04,     3.87912102631035228E-06, &
       7.28688540119691412E-04,     1.26566373053457758E-03, &
       1.62518158372674427E-03,     1.83203153216373172E-03, &
       1.91588388990527909E-03,     1.90588846755546138E-03, &
       1.82798982421825727E-03,     1.70389506421121530E-03, &
       1.55097127171097686E-03,     1.38261421852276159E-03/
  DATA ALFA(171), ALFA(172), ALFA(173), ALFA(174), ALFA(175), &
       ALFA(176), ALFA(177), ALFA(178), ALFA(179), ALFA(180)/ &
       1.20881424230064774E-03,     1.03676532638344962E-03, &
       8.71437918068619115E-04,     7.16080155297701002E-04, &
       5.72637002558129372E-04,     4.42089819465802277E-04, &
       3.24724948503090564E-04,     2.20342042730246599E-04, &
       1.28412898401353882E-04,     4.82005924552095464E-05/
  DATA BETA(1), BETA(2), BETA(3), BETA(4), BETA(5), BETA(6), &
       BETA(7), BETA(8), BETA(9), BETA(10), BETA(11), BETA(12), &
       BETA(13), BETA(14), BETA(15), BETA(16), BETA(17), BETA(18), &
       BETA(19), BETA(20), BETA(21), BETA(22)/ &
       1.79988721413553309E-02,     5.59964911064388073E-03, &
       2.88501402231132779E-03,     1.80096606761053941E-03, &
       1.24753110589199202E-03,     9.22878876572938311E-04, &
       7.14430421727287357E-04,     5.71787281789704872E-04, &
       4.69431007606481533E-04,     3.93232835462916638E-04, &
       3.34818889318297664E-04,     2.88952148495751517E-04, &
       2.52211615549573284E-04,     2.22280580798883327E-04, &
       1.97541838033062524E-04,     1.76836855019718004E-04, &
       1.59316899661821081E-04,     1.44347930197333986E-04, &
       1.31448068119965379E-04,     1.20245444949302884E-04, &
       1.10449144504599392E-04,     1.01828770740567258E-04/
  DATA BETA(23), BETA(24), BETA(25), BETA(26), BETA(27), BETA(28), &
       BETA(29), BETA(30), BETA(31), BETA(32), BETA(33), BETA(34), &
       BETA(35), BETA(36), BETA(37), BETA(38), BETA(39), BETA(40), &
       BETA(41), BETA(42), BETA(43), BETA(44)/ &
       9.41998224204237509E-05,     8.74130545753834437E-05, &
       8.13466262162801467E-05,     7.59002269646219339E-05, &
       7.09906300634153481E-05,     6.65482874842468183E-05, &
       6.25146958969275078E-05,     5.88403394426251749E-05, &
      -1.49282953213429172E-03,    -8.78204709546389328E-04, &
      -5.02916549572034614E-04,    -2.94822138512746025E-04, &
      -1.75463996970782828E-04,    -1.04008550460816434E-04, &
      -5.96141953046457895E-05,    -3.12038929076098340E-05, &
      -1.26089735980230047E-05,    -2.42892608575730389E-07, &
       8.05996165414273571E-06,     1.36507009262147391E-05, &
       1.73964125472926261E-05,     1.98672978842133780E-05/
  DATA BETA(45), BETA(46), BETA(47), BETA(48), BETA(49), BETA(50), &
       BETA(51), BETA(52), BETA(53), BETA(54), BETA(55), BETA(56), &
       BETA(57), BETA(58), BETA(59), BETA(60), BETA(61), BETA(62), &
       BETA(63), BETA(64), BETA(65), BETA(66)/ &
       2.14463263790822639E-05,     2.23954659232456514E-05, &
       2.28967783814712629E-05,     2.30785389811177817E-05, &
       2.30321976080909144E-05,     2.28236073720348722E-05, &
       2.25005881105292418E-05,     2.20981015361991429E-05, &
       2.16418427448103905E-05,     2.11507649256220843E-05, &
       2.06388749782170737E-05,     2.01165241997081666E-05, &
       1.95913450141179244E-05,     1.90689367910436740E-05, &
       1.85533719641636667E-05,     1.80475722259674218E-05, &
       5.52213076721292790E-04,     4.47932581552384646E-04, &
       2.79520653992020589E-04,     1.52468156198446602E-04, &
       6.93271105657043598E-05,     1.76258683069991397E-05/
  DATA BETA(67), BETA(68), BETA(69), BETA(70), BETA(71), BETA(72), &
       BETA(73), BETA(74), BETA(75), BETA(76), BETA(77), BETA(78), &
       BETA(79), BETA(80), BETA(81), BETA(82), BETA(83), BETA(84), &
       BETA(85), BETA(86), BETA(87), BETA(88)/ &
      -1.35744996343269136E-05,    -3.17972413350427135E-05, &
      -4.18861861696693365E-05,    -4.69004889379141029E-05, &
      -4.87665447413787352E-05,    -4.87010031186735069E-05, &
      -4.74755620890086638E-05,    -4.55813058138628452E-05, &
      -4.33309644511266036E-05,    -4.09230193157750364E-05, &
      -3.84822638603221274E-05,    -3.60857167535410501E-05, &
      -3.37793306123367417E-05,    -3.15888560772109621E-05, &
      -2.95269561750807315E-05,    -2.75978914828335759E-05, &
      -2.58006174666883713E-05,    -2.41308356761280200E-05, &
      -2.25823509518346033E-05,    -2.11479656768912971E-05, &
      -1.98200638885294927E-05,    -1.85909870801065077E-05/
  DATA BETA(89), BETA(90), BETA(91), BETA(92), BETA(93), BETA(94), &
       BETA(95), BETA(96), BETA(97), BETA(98), BETA(99), BETA(100), &
       BETA(101), BETA(102), BETA(103), BETA(104), BETA(105), &
       BETA(106), BETA(107), BETA(108), BETA(109), BETA(110)/ &
      -1.74532699844210224E-05,    -1.63997823854497997E-05, &
      -4.74617796559959808E-04,    -4.77864567147321487E-04, &
      -3.20390228067037603E-04,    -1.61105016119962282E-04, &
      -4.25778101285435204E-05,     3.44571294294967503E-05, &
       7.97092684075674924E-05,     1.03138236708272200E-04, &
       1.12466775262204158E-04,     1.13103642108481389E-04, &
       1.08651634848774268E-04,     1.01437951597661973E-04, &
       9.29298396593363896E-05,     8.40293133016089978E-05, &
       7.52727991349134062E-05,     6.69632521975730872E-05, &
       5.92564547323194704E-05,     5.22169308826975567E-05, &
       4.58539485165360646E-05,     4.01445513891486808E-05/
  DATA BETA(111), BETA(112), BETA(113), BETA(114), BETA(115), &
       BETA(116), BETA(117), BETA(118), BETA(119), BETA(120), &
       BETA(121), BETA(122), BETA(123), BETA(124), BETA(125), &
       BETA(126), BETA(127), BETA(128), BETA(129), BETA(130)/ &
       3.50481730031328081E-05,     3.05157995034346659E-05, &
       2.64956119950516039E-05,     2.29363633690998152E-05, &
       1.97893056664021636E-05,     1.70091984636412623E-05, &
       1.45547428261524004E-05,     1.23886640995878413E-05, &
       1.04775876076583236E-05,     8.79179954978479373E-06, &
       7.36465810572578444E-04,     8.72790805146193976E-04, &
       6.22614862573135066E-04,     2.85998154194304147E-04, &
       3.84737672879366102E-06,    -1.87906003636971558E-04, &
      -2.97603646594554535E-04,    -3.45998126832656348E-04, &
      -3.53382470916037712E-04,    -3.35715635775048757E-04/
  DATA BETA(131), BETA(132), BETA(133), BETA(134), BETA(135), &
       BETA(136), BETA(137), BETA(138), BETA(139), BETA(140), &
       BETA(141), BETA(142), BETA(143), BETA(144), BETA(145), &
       BETA(146), BETA(147), BETA(148), BETA(149), BETA(150)/ &
      -3.04321124789039809E-04,    -2.66722723047612821E-04, &
      -2.27654214122819527E-04,    -1.89922611854562356E-04, &
      -1.55058918599093870E-04,    -1.23778240761873630E-04, &
      -9.62926147717644187E-05,    -7.25178327714425337E-05, &
      -5.22070028895633801E-05,    -3.50347750511900522E-05, &
      -2.06489761035551757E-05,    -8.70106096849767054E-06, &
       1.13698686675100290E-06,     9.16426474122778849E-06, &
       1.56477785428872620E-05,     2.08223629482466847E-05, &
       2.48923381004595156E-05,     2.80340509574146325E-05, &
       3.03987774629861915E-05,     3.21156731406700616E-05/
  DATA BETA(151), BETA(152), BETA(153), BETA(154), BETA(155), &
       BETA(156), BETA(157), BETA(158), BETA(159), BETA(160), &
       BETA(161), BETA(162), BETA(163), BETA(164), BETA(165), &
       BETA(166), BETA(167), BETA(168), BETA(169), BETA(170)/ &
      -1.80182191963885708E-03,    -2.43402962938042533E-03, &
      -1.83422663549856802E-03,    -7.62204596354009765E-04, &
       2.39079475256927218E-04,     9.49266117176881141E-04, &
       1.34467449701540359E-03,     1.48457495259449178E-03, &
       1.44732339830617591E-03,     1.30268261285657186E-03, &
       1.10351597375642682E-03,     8.86047440419791759E-04, &
       6.73073208165665473E-04,     4.77603872856582378E-04, &
       3.05991926358789362E-04,     1.60315694594721630E-04, &
       4.00749555270613286E-05,    -5.66607461635251611E-05, &
      -1.32506186772982638E-04,    -1.90296187989614057E-04/
  DATA BETA(171), BETA(172), BETA(173), BETA(174), BETA(175), &
       BETA(176), BETA(177), BETA(178), BETA(179), BETA(180), &
       BETA(181), BETA(182), BETA(183), BETA(184), BETA(185), &
       BETA(186), BETA(187), BETA(188), BETA(189), BETA(190)/ &
      -2.32811450376937408E-04,    -2.62628811464668841E-04, &
      -2.82050469867598672E-04,    -2.93081563192861167E-04, &
      -2.97435962176316616E-04,    -2.96557334239348078E-04, &
      -2.91647363312090861E-04,    -2.83696203837734166E-04, &
      -2.73512317095673346E-04,    -2.61750155806768580E-04, &
       6.38585891212050914E-03,     9.62374215806377941E-03, &
       7.61878061207001043E-03,     2.83219055545628054E-03, &
      -2.09841352012720090E-03,    -5.73826764216626498E-03, &
      -7.70804244495414620E-03,    -8.21011692264844401E-03, &
      -7.65824520346905413E-03,    -6.47209729391045177E-03/
  DATA BETA(191), BETA(192), BETA(193), BETA(194), BETA(195), &
       BETA(196), BETA(197), BETA(198), BETA(199), BETA(200), &
       BETA(201), BETA(202), BETA(203), BETA(204), BETA(205), &
       BETA(206), BETA(207), BETA(208), BETA(209), BETA(210)/ &
      -4.99132412004966473E-03,    -3.45612289713133280E-03, &
      -2.01785580014170775E-03,    -7.59430686781961401E-04, &
       2.84173631523859138E-04,     1.10891667586337403E-03, &
       1.72901493872728771E-03,     2.16812590802684701E-03, &
       2.45357710494539735E-03,     2.61281821058334862E-03, &
       2.67141039656276912E-03,     2.65203073395980430E-03, &
       2.57411652877287315E-03,     2.45389126236094427E-03, &
       2.30460058071795494E-03,     2.13684837686712662E-03, &
       1.95896528478870911E-03,     1.77737008679454412E-03, &
       1.59690280765839059E-03,     1.42111975664438546E-03/
  DATA GAMA(1), GAMA(2), GAMA(3), GAMA(4), GAMA(5), GAMA(6), &
       GAMA(7), GAMA(8), GAMA(9), GAMA(10), GAMA(11), GAMA(12), &
       GAMA(13), GAMA(14), GAMA(15), GAMA(16), GAMA(17), GAMA(18), &
       GAMA(19), GAMA(20), GAMA(21), GAMA(22)/ &
       6.29960524947436582E-01,     2.51984209978974633E-01, &
       1.54790300415655846E-01,     1.10713062416159013E-01, &
       8.57309395527394825E-02,     6.97161316958684292E-02, &
       5.86085671893713576E-02,     5.04698873536310685E-02, &
       4.42600580689154809E-02,     3.93720661543509966E-02, &
       3.54283195924455368E-02,     3.21818857502098231E-02, &
       2.94646240791157679E-02,     2.71581677112934479E-02, &
       2.51768272973861779E-02,     2.34570755306078891E-02, &
       2.19508390134907203E-02,     2.06210828235646240E-02, &
       1.94388240897880846E-02,     1.83810633800683158E-02, &
       1.74293213231963172E-02,     1.65685837786612353E-02/
  DATA GAMA(23), GAMA(24), GAMA(25), GAMA(26), GAMA(27), GAMA(28), &
       GAMA(29), GAMA(30)/ &
       1.57865285987918445E-02,     1.50729501494095594E-02, &
       1.44193250839954639E-02,     1.38184805735341786E-02, &
       1.32643378994276568E-02,     1.27517121970498651E-02, &
       1.22761545318762767E-02,     1.18338262398482403E-02/
  DATA EX1, EX2, HPI, PI, THPI / &
       3.33333333333333333E-01,     6.66666666666666667E-01, &
       1.57079632679489662E+00,     3.14159265358979324E+00, &
       4.71238898038468986E+00/
  DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CUNHJ
  RFNU = 1.0E0/FNU
!     ZB = Z*CMPLX(RFNU,0.0E0)
!-----------------------------------------------------------------------
!     OVERFLOW TEST (Z/FNU TOO SMALL)
!-----------------------------------------------------------------------
  TSTR = REAL(Z)
  TSTI = AIMAG(Z)
  TEST = R1MACH(1)*1.0E+3
  AC = FNU*TEST
  if (ABS(TSTR) > AC .OR. ABS(TSTI) > AC) go to 15
  AC = 2.0E0*ABS(ALOG(TEST))+FNU
  ZETA1 = CMPLX(AC,0.0E0)
  ZETA2 = CMPLX(FNU,0.0E0)
  PHI=CONE
  ARG=CONE
  return
   15 CONTINUE
  ZB = Z*CMPLX(RFNU,0.0E0)
  RFNU2 = RFNU*RFNU
!-----------------------------------------------------------------------
!     COMPUTE IN THE FOURTH QUADRANT
!-----------------------------------------------------------------------
  FN13 = FNU**EX1
  FN23 = FN13*FN13
  RFN13 = CMPLX(1.0E0/FN13,0.0E0)
  W2 = CONE - ZB*ZB
  AW2 = ABS(W2)
  if (AW2 > 0.25E0) go to 130
!-----------------------------------------------------------------------
!     POWER SERIES FOR ABS(W2) <= 0.25E0
!-----------------------------------------------------------------------
  K = 1
  P(1) = CONE
  SUMA = CMPLX(GAMA(1),0.0E0)
  AP(1) = 1.0E0
  if (AW2 < TOL) go to 20
  DO 10 K=2,30
    P(K) = P(K-1)*W2
    SUMA = SUMA + P(K)*CMPLX(GAMA(K),0.0E0)
    AP(K) = AP(K-1)*AW2
    if (AP(K) < TOL) go to 20
   10 CONTINUE
  K = 30
   20 CONTINUE
  KMAX = K
  ZETA = W2*SUMA
  ARG = ZETA*CMPLX(FN23,0.0E0)
  ZA = CSQRT(SUMA)
  ZETA2 = CSQRT(W2)*CMPLX(FNU,0.0E0)
  ZETA1 = ZETA2*(CONE+ZETA*ZA*CMPLX(EX2,0.0E0))
  ZA = ZA + ZA
  PHI = CSQRT(ZA)*RFN13
  if (IPMTR == 1) go to 120
!-----------------------------------------------------------------------
!     SUM SERIES FOR ASUM AND BSUM
!-----------------------------------------------------------------------
  SUMB = CZERO
  DO 30 K=1,KMAX
    SUMB = SUMB + P(K)*CMPLX(BETA(K),0.0E0)
   30 CONTINUE
  ASUM = CZERO
  BSUM = SUMB
  L1 = 0
  L2 = 30
  BTOL = TOL*ABS(BSUM)
  ATOL = TOL
  PP = 1.0E0
  IAS = 0
  IBS = 0
  if (RFNU2 < TOL) go to 110
  DO 100 IS=2,7
    ATOL = ATOL/RFNU2
    PP = PP*RFNU2
    if (IAS == 1) go to 60
    SUMA = CZERO
    DO 40 K=1,KMAX
      M = L1 + K
      SUMA = SUMA + P(K)*CMPLX(ALFA(M),0.0E0)
      if (AP(K) < ATOL) go to 50
   40   CONTINUE
   50   CONTINUE
    ASUM = ASUM + SUMA*CMPLX(PP,0.0E0)
    if (PP < TOL) IAS = 1
   60   CONTINUE
    if (IBS == 1) go to 90
    SUMB = CZERO
    DO 70 K=1,KMAX
      M = L2 + K
      SUMB = SUMB + P(K)*CMPLX(BETA(M),0.0E0)
      if (AP(K) < ATOL) go to 80
   70   CONTINUE
   80   CONTINUE
    BSUM = BSUM + SUMB*CMPLX(PP,0.0E0)
    if (PP < BTOL) IBS = 1
   90   CONTINUE
    if (IAS == 1 .AND. IBS == 1) go to 110
    L1 = L1 + 30
    L2 = L2 + 30
  100 CONTINUE
  110 CONTINUE
  ASUM = ASUM + CONE
  PP = RFNU*REAL(RFN13)
  BSUM = BSUM*CMPLX(PP,0.0E0)
  120 CONTINUE
  return
!-----------------------------------------------------------------------
!     ABS(W2) > 0.25E0
!-----------------------------------------------------------------------
  130 CONTINUE
  W = CSQRT(W2)
  WR = REAL(W)
  WI = AIMAG(W)
  if (WR < 0.0E0) WR = 0.0E0
  if (WI < 0.0E0) WI = 0.0E0
  W = CMPLX(WR,WI)
  ZA = (CONE+W)/ZB
  ZC = CLOG(ZA)
  ZCR = REAL(ZC)
  ZCI = AIMAG(ZC)
  if (ZCI < 0.0E0) ZCI = 0.0E0
  if (ZCI > HPI) ZCI = HPI
  if (ZCR < 0.0E0) ZCR = 0.0E0
  ZC = CMPLX(ZCR,ZCI)
  ZTH = (ZC-W)*CMPLX(1.5E0,0.0E0)
  CFNU = CMPLX(FNU,0.0E0)
  ZETA1 = ZC*CFNU
  ZETA2 = W*CFNU
  AZTH = ABS(ZTH)
  ZTHR = REAL(ZTH)
  ZTHI = AIMAG(ZTH)
  ANG = THPI
  if (ZTHR >= 0.0E0 .AND. ZTHI < 0.0E0) go to 140
  ANG = HPI
  if (ZTHR == 0.0E0) go to 140
  ANG = ATAN(ZTHI/ZTHR)
  if (ZTHR < 0.0E0) ANG = ANG + PI
  140 CONTINUE
  PP = AZTH**EX2
  ANG = ANG*EX2
  ZETAR = PP*COS(ANG)
  ZETAI = PP*SIN(ANG)
  if (ZETAI < 0.0E0) ZETAI = 0.0E0
  ZETA = CMPLX(ZETAR,ZETAI)
  ARG = ZETA*CMPLX(FN23,0.0E0)
  RTZTA = ZTH/ZETA
  ZA = RTZTA/W
  PHI = CSQRT(ZA+ZA)*RFN13
  if (IPMTR == 1) go to 120
  TFN = CMPLX(RFNU,0.0E0)/W
  RZTH = CMPLX(RFNU,0.0E0)/ZTH
  ZC = RZTH*CMPLX(AR(2),0.0E0)
  T2 = CONE/W2
  UP(2) = (T2*CMPLX(C(2),0.0E0)+CMPLX(C(3),0.0E0))*TFN
  BSUM = UP(2) + ZC
  ASUM = CZERO
  if (RFNU < TOL) go to 220
  PRZTH = RZTH
  PTFN = TFN
  UP(1) = CONE
  PP = 1.0E0
  BSUMR = REAL(BSUM)
  BSUMI = AIMAG(BSUM)
  BTOL = TOL*(ABS(BSUMR)+ABS(BSUMI))
  KS = 0
  KP1 = 2
  L = 3
  IAS = 0
  IBS = 0
  DO 210 LR=2,12,2
    LRP1 = LR + 1
!-----------------------------------------------------------------------
!     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
!     NEXT SUMA AND SUMB
!-----------------------------------------------------------------------
    DO 160 K=LR,LRP1
      KS = KS + 1
      KP1 = KP1 + 1
      L = L + 1
      ZA = CMPLX(C(L),0.0E0)
      DO 150 J=2,KP1
        L = L + 1
        ZA = ZA*T2 + CMPLX(C(L),0.0E0)
  150     CONTINUE
      PTFN = PTFN*TFN
      UP(KP1) = PTFN*ZA
      CR(KS) = PRZTH*CMPLX(BR(KS+1),0.0E0)
      PRZTH = PRZTH*RZTH
      DR(KS) = PRZTH*CMPLX(AR(KS+2),0.0E0)
  160   CONTINUE
    PP = PP*RFNU2
    if (IAS == 1) go to 180
    SUMA = UP(LRP1)
    JU = LRP1
    DO 170 JR=1,LR
      JU = JU - 1
      SUMA = SUMA + CR(JR)*UP(JU)
  170   CONTINUE
    ASUM = ASUM + SUMA
    ASUMR = REAL(ASUM)
    ASUMI = AIMAG(ASUM)
    TEST = ABS(ASUMR) + ABS(ASUMI)
    if (PP < TOL .AND. TEST < TOL) IAS = 1
  180   CONTINUE
    if (IBS == 1) go to 200
    SUMB = UP(LR+2) + UP(LRP1)*ZC
    JU = LRP1
    DO 190 JR=1,LR
      JU = JU - 1
      SUMB = SUMB + DR(JR)*UP(JU)
  190   CONTINUE
    BSUM = BSUM + SUMB
    BSUMR = REAL(BSUM)
    BSUMI = AIMAG(BSUM)
    TEST = ABS(BSUMR) + ABS(BSUMI)
    if (PP < BTOL .AND. TEST < TOL) IBS = 1
  200   CONTINUE
    if (IAS == 1 .AND. IBS == 1) go to 220
  210 CONTINUE
  220 CONTINUE
  ASUM = ASUM + CONE
  BSUM = -BSUM*RFN13/RTZTA
  go to 120
end
subroutine CMLRI (Z, FNU, KODE, N, Y, NZ, TOL)
!
!! CMLRI is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CMLRI-A, ZMLRI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY THE
!     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  GAMLN, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CMLRI
  COMPLEX CK, CNORM, CONE, CTWO, CZERO, PT, P1, P2, RZ, SUM, Y, Z
  REAL ACK, AK, AP, AT, AZ, BK, FKAP, FKK, FLAM, FNF, FNU, RHO, &
   RHO2, SCLE, TFNF, TOL, TST, X, GAMLN, R1MACH
  INTEGER I, IAZ, IDUM, IFNU, INU, ITIME, K, KK, KM, KODE, M, N, NZ
  DIMENSION Y(N)
  DATA CZERO,CONE,CTWO /(0.0E0,0.0E0),(1.0E0,0.0E0),(2.0E0,0.0E0)/
  SCLE = 1.0E+3*R1MACH(1)/TOL
!***FIRST EXECUTABLE STATEMENT  CMLRI
  NZ=0
  AZ = ABS(Z)
  X = REAL(Z)
  IAZ = AZ
  IFNU = FNU
  INU = IFNU + N - 1
  AT = IAZ + 1.0E0
  CK = CMPLX(AT,0.0E0)/Z
  RZ = CTWO/Z
  P1 = CZERO
  P2 = CONE
  ACK = (AT+1.0E0)/AZ
  RHO = ACK + SQRT(ACK*ACK-1.0E0)
  RHO2 = RHO*RHO
  TST = (RHO2+RHO2)/((RHO2-1.0E0)*(RHO-1.0E0))
  TST = TST/TOL
!-----------------------------------------------------------------------
!     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
!-----------------------------------------------------------------------
  AK = AT
  DO 10 I=1,80
    PT = P2
    P2 = P1 - CK*P2
    P1 = PT
    CK = CK + RZ
    AP = ABS(P2)
    if (AP > TST*AK*AK) go to 20
    AK = AK + 1.0E0
   10 CONTINUE
  go to 110
   20 CONTINUE
  I = I + 1
  K = 0
  if (INU < IAZ) go to 40
!-----------------------------------------------------------------------
!     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
!-----------------------------------------------------------------------
  P1 = CZERO
  P2 = CONE
  AT = INU + 1.0E0
  CK = CMPLX(AT,0.0E0)/Z
  ACK = AT/AZ
  TST = SQRT(ACK/TOL)
  ITIME = 1
  DO 30 K=1,80
    PT = P2
    P2 = P1 - CK*P2
    P1 = PT
    CK = CK + RZ
    AP = ABS(P2)
    if (AP < TST) go to 30
    if (ITIME == 2) go to 40
    ACK = ABS(CK)
    FLAM = ACK + SQRT(ACK*ACK-1.0E0)
    FKAP = AP/ABS(P1)
    RHO = MIN(FLAM,FKAP)
    TST = TST*SQRT(RHO/(RHO*RHO-1.0E0))
    ITIME = 2
   30 CONTINUE
  go to 110
   40 CONTINUE
!-----------------------------------------------------------------------
!     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
!-----------------------------------------------------------------------
  K = K + 1
  KK = MAX(I+IAZ,K+INU)
  FKK = KK
  P1 = CZERO
!-----------------------------------------------------------------------
!     SCALE P2 AND SUM BY SCLE
!-----------------------------------------------------------------------
  P2 = CMPLX(SCLE,0.0E0)
  FNF = FNU - IFNU
  TFNF = FNF + FNF
  BK = GAMLN(FKK+TFNF+1.0E0,IDUM) - GAMLN(FKK+1.0E0,IDUM) &
       -GAMLN(TFNF+1.0E0,IDUM)
  BK = EXP(BK)
  SUM = CZERO
  KM = KK - INU
  DO 50 I=1,KM
    PT = P2
    P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
    P1 = PT
    AK = 1.0E0 - TFNF/(FKK+TFNF)
    ACK = BK*AK
    SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
    BK = ACK
    FKK = FKK - 1.0E0
   50 CONTINUE
  Y(N) = P2
  if (N == 1) go to 70
  DO 60 I=2,N
    PT = P2
    P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
    P1 = PT
    AK = 1.0E0 - TFNF/(FKK+TFNF)
    ACK = BK*AK
    SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
    BK = ACK
    FKK = FKK - 1.0E0
    M = N - I + 1
    Y(M) = P2
   60 CONTINUE
   70 CONTINUE
  if (IFNU <= 0) go to 90
  DO 80 I=1,IFNU
    PT = P2
    P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
    P1 = PT
    AK = 1.0E0 - TFNF/(FKK+TFNF)
    ACK = BK*AK
    SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
    BK = ACK
    FKK = FKK - 1.0E0
   80 CONTINUE
   90 CONTINUE
  PT = Z
  if (KODE == 2) PT = PT - CMPLX(X,0.0E0)
  P1 = -CMPLX(FNF,0.0E0)*CLOG(RZ) + PT
  AP = GAMLN(1.0E0+FNF,IDUM)
  PT = P1 - CMPLX(AP,0.0E0)
!-----------------------------------------------------------------------
!     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
!     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
!-----------------------------------------------------------------------
  P2 = P2 + SUM
  AP = ABS(P2)
  P1 = CMPLX(1.0E0/AP,0.0E0)
  CK = CEXP(PT)*P1
  PT = CONJG(P2)*P1
  CNORM = CK*PT
  DO 100 I=1,N
    Y(I) = Y(I)*CNORM
  100 CONTINUE
  return
  110 CONTINUE
  NZ=-2
  return
end
subroutine CAIRY (Z, ID, KODE, AI, NZ, IERR)
!
!! CAIRY computes the Airy function Ai(z) or its derivative dAi/dz ...
!            for complex argument z.  A scaling option is available ...
!            to help avoid underflow and overflow.
!
!***LIBRARY   SLATEC
!***CATEGORY  C10D
!***TYPE      COMPLEX (CAIRY-C, ZAIRY-C)
!***KEYWORDS  AIRY FUNCTION, BESSEL FUNCTION OF ORDER ONE THIRD,
!             BESSEL FUNCTION OF ORDER TWO THIRDS
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!         On KODE=1, CAIRY computes the complex Airy function Ai(z)
!         or its derivative dAi/dz on ID=0 or ID=1 respectively. On
!         KODE=2, a scaling option exp(zeta)*Ai(z) or exp(zeta)*dAi/dz
!         is provided to remove the exponential decay in -pi/3<arg(z)
!         <pi/3 and the exponential growth in pi/3<abs(arg(z))<pi where
!         zeta=(2/3)*z**(3/2).
!
!         While the Airy functions Ai(z) and dAi/dz are analytic in
!         the whole z-plane, the corresponding scaled functions defined
!         for KODE=2 have a cut along the negative real axis.
!
!         Input
!           Z      - Argument of type COMPLEX
!           ID     - Order of derivative, ID=0 or ID=1
!           KODE   - A parameter to indicate the scaling option
!                    KODE=1  returns
!                            AI=Ai(z)  on ID=0
!                            AI=dAi/dz on ID=1
!                            at z=Z
!                        =2  returns
!                            AI=exp(zeta)*Ai(z)  on ID=0
!                            AI=exp(zeta)*dAi/dz on ID=1
!                            at z=Z where zeta=(2/3)*z**(3/2)
!
!         Output
!           AI     - Result of type COMPLEX
!           NZ     - Underflow indicator
!                    NZ=0    Normal return
!                    NZ=1    AI=0 due to underflow in
!                            -pi/3<arg(Z)<pi/3 on KODE=1
!           IERR   - Error flag
!                    IERR=0  Normal return     - COMPUTATION COMPLETED
!                    IERR=1  Input error       - NO COMPUTATION
!                    IERR=2  Overflow          - NO COMPUTATION
!                            (Re(Z) too large with KODE=1)
!                    IERR=3  Precision warning - COMPUTATION COMPLETED
!                            (Result has less than half precision)
!                    IERR=4  Precision error   - NO COMPUTATION
!                            (Result has no precision)
!                    IERR=5  Algorithmic error - NO COMPUTATION
!                            (Termination condition not met)
!
! *Long Description:
!
!         Ai(z) and dAi/dz are computed from K Bessel functions by
!
!                Ai(z) =  c*sqrt(z)*K(1/3,zeta)
!               dAi/dz = -c*   z   *K(2/3,zeta)
!                    c =  1/(pi*sqrt(3))
!                 zeta =  (2/3)*z**(3/2)
!
!         when abs(z)>1 and from power series when abs(z)<=1.
!
!         In most complex variable computation, one must evaluate ele-
!         mentary functions.  When the magnitude of Z is large, losses
!         of significance by argument reduction occur.  Consequently, if
!         the magnitude of ZETA=(2/3)*Z**(3/2) exceeds U1=SQRT(0.5/UR),
!         then losses exceeding half precision are likely and an error
!         flag IERR=3 is triggered where UR=R1MACH(4)=UNIT ROUNDOFF.
!         Also, if the magnitude of ZETA is larger than U2=0.5/UR, then
!         all significance is lost and IERR=4.  In order to use the INT
!         function, ZETA must be further restricted not to exceed
!         U3=I1MACH(9)=LARGEST INTEGER.  Thus, the magnitude of ZETA
!         must be restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2,
!         and U3 are approximately 2.0E+3, 4.2E+6, 2.1E+9 in single
!         precision and 4.7E+7, 2.3E+15, 2.1E+9 in double precision.
!         This makes U2 limiting is single precision and U3 limiting
!         in double precision.  This means that the magnitude of Z
!         cannot exceed approximately 3.4E+4 in single precision and
!         2.1E+6 in double precision.  This also means that one can
!         expect to retain, in the worst cases on 32-bit machines,
!         no digits in single precision and only 6 digits in double
!         precision.
!
!         The approximate relative error in the magnitude of a complex
!         Bessel function can be expressed as P*10**S where P=MAX(UNIT
!         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre-
!         sents the increase in error due to argument reduction in the
!         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))),
!         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF
!         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may
!         have only absolute accuracy.  This is most likely to occur
!         when one component (in magnitude) is larger than the other by
!         several orders of magnitude.  If one component is 10**K larger
!         than the other, then one can expect only MAX(ABS(LOG10(P))-K,
!         0) significant digits; or, stated another way, when K exceeds
!         the exponent of P, no significant digits remain in the smaller
!         component.  However, the phase angle retains absolute accuracy
!         because, in complex arithmetic with precision P, the smaller
!         component will not (as a rule) decrease below P times the
!         magnitude of the larger component. In these extreme cases,
!         the principal phase angle is on the order of +P, -P, PI/2-P,
!         or -PI/2+P.
!
!***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
!                 matical Functions, National Bureau of Standards
!                 Applied Mathematics Series 55, U. S. Department
!                 of Commerce, Tenth Printing (1972) or later.
!               2. D. E. Amos, Computation of Bessel Functions of
!                 Complex Argument and Large Order, Report SAND83-0643,
!                 Sandia National Laboratories, Albuquerque, NM, May
!                 1983.
!               3. D. E. Amos, A Subroutine Package for Bessel Functions
!                 of a Complex Argument and Nonnegative Order, Report
!                 SAND85-1018, Sandia National Laboratory, Albuquerque,
!                 NM, May 1985.
!               4. D. E. Amos, A portable package for Bessel functions
!                 of a complex argument and nonnegative order, ACM
!                 Transactions on Mathematical Software, 12 (September
!                 1986), pp. 265-273.
!
!***ROUTINES CALLED  CACAI, CBKNU, I1MACH, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   890801  REVISION DATE from Version 3.2
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   920128  Category corrected.  (WRB)
!   920811  Prologue revised.  (DWL)
!***END PROLOGUE  CAIRY
  COMPLEX AI, CONE, CSQ, CY, S1, S2, TRM1, TRM2, Z, ZTA, Z3
  REAL AA, AD, AK, ALIM, ATRM, AZ, AZ3, BK, CK, COEF, C1, C2, DIG, &
   DK, D1, D2, ELIM, FID, FNU, RL, R1M5, SFAC, TOL, TTH, ZI, ZR, &
   Z3I, Z3R, R1MACH, BB, ALAZ
  INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, I1MACH
  DIMENSION CY(1)
  DATA TTH, C1, C2, COEF /6.66666666666666667E-01, &
   3.55028053887817240E-01,2.58819403792806799E-01, &
   1.83776298473930683E-01/
  DATA  CONE / (1.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CAIRY
  IERR = 0
  NZ=0
  if (ID < 0 .OR. ID > 1) IERR=1
  if (KODE < 1 .OR. KODE > 2) IERR=1
  if (IERR /= 0) RETURN
  AZ = ABS(Z)
  TOL = MAX(R1MACH(4),1.0E-18)
  FID = ID
  if (AZ > 1.0E0) go to 60
!-----------------------------------------------------------------------
!     POWER SERIES FOR ABS(Z) <= 1.
!-----------------------------------------------------------------------
  S1 = CONE
  S2 = CONE
  if (AZ < TOL) go to 160
  AA = AZ*AZ
  if (AA < TOL/AZ) go to 40
  TRM1 = CONE
  TRM2 = CONE
  ATRM = 1.0E0
  Z3 = Z*Z*Z
  AZ3 = AZ*AA
  AK = 2.0E0 + FID
  BK = 3.0E0 - FID - FID
  CK = 4.0E0 - FID
  DK = 3.0E0 + FID + FID
  D1 = AK*DK
  D2 = BK*CK
  AD = MIN(D1,D2)
  AK = 24.0E0 + 9.0E0*FID
  BK = 30.0E0 - 9.0E0*FID
  Z3R = REAL(Z3)
  Z3I = AIMAG(Z3)
  DO 30 K=1,25
    TRM1 = TRM1*CMPLX(Z3R/D1,Z3I/D1)
    S1 = S1 + TRM1
    TRM2 = TRM2*CMPLX(Z3R/D2,Z3I/D2)
    S2 = S2 + TRM2
    ATRM = ATRM*AZ3/AD
    D1 = D1 + AK
    D2 = D2 + BK
    AD = MIN(D1,D2)
    if (ATRM < TOL*AD) go to 40
    AK = AK + 18.0E0
    BK = BK + 18.0E0
   30 CONTINUE
   40 CONTINUE
  if (ID == 1) go to 50
  AI = S1*CMPLX(C1,0.0E0) - Z*S2*CMPLX(C2,0.0E0)
  if (KODE == 1) RETURN
  ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
  AI = AI*CEXP(ZTA)
  return
   50 CONTINUE
  AI = -S2*CMPLX(C2,0.0E0)
  if (AZ > TOL) AI = AI + Z*Z*S1*CMPLX(C1/(1.0E0+FID),0.0E0)
  if (KODE == 1) RETURN
  ZTA = Z*CSQRT(Z)*CMPLX(TTH,0.0E0)
  AI = AI*CEXP(ZTA)
  return
!-----------------------------------------------------------------------
!     CASE FOR ABS(Z) > 1.0
!-----------------------------------------------------------------------
   60 CONTINUE
  FNU = (1.0E0+FID)/3.0E0
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!-----------------------------------------------------------------------
  K1 = I1MACH(12)
  K2 = I1MACH(13)
  R1M5 = R1MACH(5)
  K = MIN(ABS(K1),ABS(K2))
  ELIM = 2.303E0*(K*R1M5-3.0E0)
  K1 = I1MACH(11) - 1
  AA = R1M5*K1
  DIG = MIN(AA,18.0E0)
  AA = AA*2.303E0
  ALIM = ELIM + MAX(-AA,-41.45E0)
  RL = 1.2E0*DIG + 3.0E0
  ALAZ=ALOG(AZ)
!-----------------------------------------------------------------------
!     TEST FOR RANGE
!-----------------------------------------------------------------------
  AA=0.5E0/TOL
  BB=I1MACH(9)*0.5E0
  AA=MIN(AA,BB)
  AA=AA**TTH
  if (AZ > AA) go to 260
  AA=SQRT(AA)
  if (AZ > AA) IERR=3
  CSQ=CSQRT(Z)
  ZTA=Z*CSQ*CMPLX(TTH,0.0E0)
!-----------------------------------------------------------------------
!     RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
!-----------------------------------------------------------------------
  IFLAG = 0
  SFAC = 1.0E0
  ZI = AIMAG(Z)
  ZR = REAL(Z)
  AK = AIMAG(ZTA)
  if (ZR >= 0.0E0) go to 70
  BK = REAL(ZTA)
  CK = -ABS(BK)
  ZTA = CMPLX(CK,AK)
   70 CONTINUE
  if (ZI /= 0.0E0) go to 80
  if (ZR > 0.0E0) go to 80
  ZTA = CMPLX(0.0E0,AK)
   80 CONTINUE
  AA = REAL(ZTA)
  if (AA >= 0.0E0 .AND. ZR > 0.0E0) go to 100
  if (KODE == 2) go to 90
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
  if (AA > (-ALIM)) go to 90
  AA = -AA + 0.25E0*ALAZ
  IFLAG = 1
  SFAC = TOL
  if (AA > ELIM) go to 240
   90 CONTINUE
!-----------------------------------------------------------------------
!     CBKNU AND CACAI RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
!-----------------------------------------------------------------------
  MR = 1
  if (ZI < 0.0E0) MR = -1
  call CACAI(ZTA, FNU, KODE, MR, 1, CY, NN, RL, TOL, ELIM, ALIM)
  if (NN < 0) go to 250
  NZ = NZ + NN
  go to 120
  100 CONTINUE
  if (KODE == 2) go to 110
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
  if (AA < ALIM) go to 110
  AA = -AA - 0.25E0*ALAZ
  IFLAG = 2
  SFAC = 1.0E0/TOL
  if (AA < (-ELIM)) go to 180
  110 CONTINUE
  call CBKNU(ZTA, FNU, KODE, 1, CY, NZ, TOL, ELIM, ALIM)
  120 CONTINUE
  S1 = CY(1)*CMPLX(COEF,0.0E0)
  if (IFLAG /= 0) go to 140
  if (ID == 1) go to 130
  AI = CSQ*S1
  return
  130 AI = -Z*S1
  return
  140 CONTINUE
  S1 = S1*CMPLX(SFAC,0.0E0)
  if (ID == 1) go to 150
  S1 = S1*CSQ
  AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
  return
  150 CONTINUE
  S1 = -S1*Z
  AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
  return
  160 CONTINUE
  AA = 1.0E+3*R1MACH(1)
  S1 = CMPLX(0.0E0,0.0E0)
  if (ID == 1) go to 170
  if (AZ > AA) S1 = CMPLX(C2,0.0E0)*Z
  AI = CMPLX(C1,0.0E0) - S1
  return
  170 CONTINUE
  AI = -CMPLX(C2,0.0E0)
  AA = SQRT(AA)
  if (AZ > AA) S1 = Z*Z*CMPLX(0.5E0,0.0E0)
  AI = AI + S1*CMPLX(C1,0.0E0)
  return
  180 CONTINUE
  NZ = 1
  AI = CMPLX(0.0E0,0.0E0)
  return
  240 CONTINUE
  NZ = 0
  IERR=2
  return
  250 CONTINUE
  if ( NN == (-1)) go to 240
  NZ=0
  IERR=5
  return
  260 CONTINUE
  IERR=4
  NZ=0
  return
end
subroutine CACAI (Z, FNU, KODE, MR, N, Y, NZ, RL, TOL, ELIM, ALIM)
!
!! CACAI is subsidiary to CAIRY.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CACAI-A, ZACAI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
!
!         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
!                 MP=PI*MR*CMPLX(0.0,1.0)
!
!     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
!     HALF Z PLANE FOR USE WITH CAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
!     CACAI IS THE SAME AS CACON WITH THE PARTS FOR LARGER ORDERS AND
!     RECURRENCE REMOVED. A RECURSIVE call TO CACON CAN RESULT if CACON
!     IS CALLED FROM CAIRY.
!
!***SEE ALSO  CAIRY
!***ROUTINES CALLED  CASYI, CBKNU, CMLRI, CS1S2, CSERI, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CACAI
  COMPLEX CSGN, CSPN, C1, C2, Y, Z, ZN, CY
  REAL ALIM, ARG, ASCLE, AZ, CPN, DFNU, ELIM, FMR, FNU, PI, RL, &
   SGN, SPN, TOL, YY, R1MACH
  INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
  DIMENSION Y(N), CY(2)
  DATA PI / 3.14159265358979324E0 /
!***FIRST EXECUTABLE STATEMENT  CACAI
  NZ = 0
  ZN = -Z
  AZ = ABS(Z)
  NN = N
  DFNU = FNU + (N-1)
  if (AZ <= 2.0E0) go to 10
  if (AZ*AZ*0.25E0 > DFNU+1.0E0) go to 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     POWER SERIES FOR THE I FUNCTION
!-----------------------------------------------------------------------
  call CSERI(ZN, FNU, KODE, NN, Y, NW, TOL, ELIM, ALIM)
  go to 40
   20 CONTINUE
  if (AZ < RL) go to 30
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
!-----------------------------------------------------------------------
  call CASYI(ZN, FNU, KODE, NN, Y, NW, RL, TOL, ELIM, ALIM)
  if (NW < 0) go to 70
  go to 40
   30 CONTINUE
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
!-----------------------------------------------------------------------
  call CMLRI(ZN, FNU, KODE, NN, Y, NW, TOL)
  if ( NW < 0) go to 70
   40 CONTINUE
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
!-----------------------------------------------------------------------
  call CBKNU(ZN, FNU, KODE, 1, CY, NW, TOL, ELIM, ALIM)
  if (NW /= 0) go to 70
  FMR = MR
  SGN = -SIGN(PI,FMR)
  CSGN = CMPLX(0.0E0,SGN)
  if (KODE == 1) go to 50
  YY = -AIMAG(ZN)
  CPN = COS(YY)
  SPN = SIN(YY)
  CSGN = CSGN*CMPLX(CPN,SPN)
   50 CONTINUE
!-----------------------------------------------------------------------
!     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
  INU = FNU
  ARG = (FNU-INU)*SGN
  CPN = COS(ARG)
  SPN = SIN(ARG)
  CSPN = CMPLX(CPN,SPN)
  if (MOD(INU,2) == 1) CSPN = -CSPN
  C1 = CY(1)
  C2 = Y(1)
  if (KODE == 1) go to 60
  IUF = 0
  ASCLE = 1.0E+3*R1MACH(1)/TOL
  call CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
  NZ = NZ + NW
   60 CONTINUE
  Y(1) = CSPN*C1 + CSGN*C2
  return
   70 CONTINUE
  NZ = -1
  if ( NW == (-2)) NZ=-2
  return
end
subroutine CUNI2 (Z, FNU, KODE, N, Y, NZ, NLAST, FNUL, TOL, ELIM, &
     ALIM)
!
!! CUNI2 is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNI2-A, ZUNI2-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
!     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
!     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
!
!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST /= 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
!     Y(I)=CZERO FOR I=NLAST+1,N
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  CAIRY, CUCHK, CUNHJ, CUOIK, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUNI2
  COMPLEX AI, ARG, ASUM, BSUM, CFN, CI, CID, CIP, CONE, CRSC, CSCL, &
   CSR, CSS, CY, CZERO, C1, C2, DAI, PHI, RZ, S1, S2, Y, Z, ZB, &
   ZETA1, ZETA2, ZN, ZAR
  REAL AARG, AIC, ALIM, ANG, APHI, ASCLE, AY, BRY, CAR, C2I, C2M, &
   C2R, ELIM, FN, FNU, FNUL, HPI, RS1, SAR, TOL, YY, R1MACH
  INTEGER I, IFLAG, IN, INU, J, K, KODE, N, NAI, ND, NDAI, NLAST, &
   NN, NUF, NW, NZ, IDUM
  DIMENSION BRY(3), Y(N), CIP(4), CSS(3), CSR(3), CY(2)
  DATA CZERO,CONE,CI/(0.0E0,0.0E0),(1.0E0,0.0E0),(0.0E0,1.0E0)/
  DATA CIP(1),CIP(2),CIP(3),CIP(4)/ &
   (1.0E0,0.0E0), (0.0E0,1.0E0), (-1.0E0,0.0E0), (0.0E0,-1.0E0)/
  DATA HPI, AIC  / &
        1.57079632679489662E+00,     1.265512123484645396E+00/
!***FIRST EXECUTABLE STATEMENT  CUNI2
  NZ = 0
  ND = N
  NLAST = 0
!-----------------------------------------------------------------------
!     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
!     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
!     EXP(ALIM)=EXP(ELIM)*TOL
!-----------------------------------------------------------------------
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CRSC = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CRSC
  CSR(1) = CRSC
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = 1.0E+3*R1MACH(1)/TOL
  YY = AIMAG(Z)
!-----------------------------------------------------------------------
!     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
!-----------------------------------------------------------------------
  ZN = -Z*CI
  ZB = Z
  CID = -CI
  INU = FNU
  ANG = HPI*(FNU-INU)
  CAR = COS(ANG)
  SAR = SIN(ANG)
  C2 = CMPLX(CAR,SAR)
  ZAR = C2
  IN = INU + N - 1
  IN = MOD(IN,4)
  C2 = C2*CIP(IN+1)
  if (YY > 0.0E0) go to 10
  ZN = CONJG(-ZN)
  ZB = CONJG(ZB)
  CID = -CID
  C2 = CONJG(C2)
   10 CONTINUE
!-----------------------------------------------------------------------
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!-----------------------------------------------------------------------
  FN = MAX(FNU,1.0E0)
  call CUNHJ(ZN, FN, 1, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
  if (KODE == 1) go to 20
  CFN = CMPLX(FNU,0.0E0)
  S1 = -ZETA1 + CFN*(CFN/(ZB+ZETA2))
  go to 30
   20 CONTINUE
  S1 = -ZETA1 + ZETA2
   30 CONTINUE
  RS1 = REAL(S1)
  if (ABS(RS1) > ELIM) go to 150
   40 CONTINUE
  NN = MIN(2,ND)
  DO 90 I=1,NN
    FN = FNU + (ND-I)
    call CUNHJ(ZN, FN, 0, TOL, PHI, ARG, ZETA1, ZETA2, ASUM, BSUM)
    if (KODE == 1) go to 50
    CFN = CMPLX(FN,0.0E0)
    AY = ABS(YY)
    S1 = -ZETA1 + CFN*(CFN/(ZB+ZETA2)) + CMPLX(0.0E0,AY)
    go to 60
   50   CONTINUE
    S1 = -ZETA1 + ZETA2
   60   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = REAL(S1)
    if (ABS(RS1) > ELIM) go to 120
    if (I == 1) IFLAG = 2
    if (ABS(RS1) < ALIM) go to 70
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    APHI = ABS(PHI)
    AARG = ABS(ARG)
    RS1 = RS1 + ALOG(APHI) - 0.25E0*ALOG(AARG) - AIC
    if (ABS(RS1) > ELIM) go to 120
    if (I == 1) IFLAG = 1
    if (RS1 < 0.0E0) go to 70
    if (I == 1) IFLAG = 3
   70   CONTINUE
!-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!-----------------------------------------------------------------------
    call CAIRY(ARG, 0, 2, AI, NAI, IDUM)
    call CAIRY(ARG, 1, 2, DAI, NDAI, IDUM)
    S2 = PHI*(AI*ASUM+DAI*BSUM)
    C2R = REAL(S1)
    C2I = AIMAG(S1)
    C2M = EXP(C2R)*REAL(CSS(IFLAG))
    S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
    S2 = S2*S1
    if (IFLAG /= 1) go to 80
    call CUCHK(S2, NW, BRY(1), TOL)
    if (NW /= 0) go to 120
   80   CONTINUE
    if (YY <= 0.0E0) S2 = CONJG(S2)
    J = ND - I + 1
    S2 = S2*C2
    CY(I) = S2
    Y(J) = S2*CSR(IFLAG)
    C2 = C2*CID
   90 CONTINUE
  if (ND <= 2) go to 110
  RZ = CMPLX(2.0E0,0.0E0)/Z
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = R1MACH(2)
  S1 = CY(1)
  S2 = CY(2)
  C1 = CSR(IFLAG)
  ASCLE = BRY(IFLAG)
  K = ND - 2
  FN = K
  DO 100 I=3,ND
    C2 = S2
    S2 = S1 + CMPLX(FNU+FN,0.0E0)*RZ*S2
    S1 = C2
    C2 = S2*C1
    Y(K) = C2
    K = K - 1
    FN = FN - 1.0E0
    if (IFLAG >= 3) go to 100
    C2R = REAL(C2)
    C2I = AIMAG(C2)
    C2R = ABS(C2R)
    C2I = ABS(C2I)
    C2M = MAX(C2R,C2I)
    if (C2M <= ASCLE) go to 100
    IFLAG = IFLAG + 1
    ASCLE = BRY(IFLAG)
    S1 = S1*C1
    S2 = C2
    S1 = S1*CSS(IFLAG)
    S2 = S2*CSS(IFLAG)
    C1 = CSR(IFLAG)
  100 CONTINUE
  110 CONTINUE
  return
  120 CONTINUE
  if (RS1 > 0.0E0) go to 140
!-----------------------------------------------------------------------
!     SET UNDERFLOW AND UPDATE PARAMETERS
!-----------------------------------------------------------------------
  Y(ND) = CZERO
  NZ = NZ + 1
  ND = ND - 1
  if (ND == 0) go to 110
  call CUOIK(Z, FNU, KODE, 1, ND, Y, NUF, TOL, ELIM, ALIM)
  if (NUF < 0) go to 140
  ND = ND - NUF
  NZ = NZ + NUF
  if (ND == 0) go to 110
  FN = FNU + (ND-1)
  if (FN < FNUL) go to 130
!      FN = AIMAG(CID)
!      J = NUF + 1
!      K = MOD(J,4) + 1
!      S1 = CIP(K)
!      if (FN < 0.0E0) S1 = CONJG(S1)
!      C2 = C2*S1
  IN = INU + ND - 1
  IN = MOD(IN,4) + 1
  C2 = ZAR*CIP(IN)
  if (YY <= 0.0E0)C2=CONJG(C2)
  go to 40
  130 CONTINUE
  NLAST = ND
  return
  140 CONTINUE
  NZ = -1
  return
  150 CONTINUE
  if (RS1 > 0.0E0) go to 140
  NZ = N
  DO 160 I=1,N
    Y(I) = CZERO
  160 CONTINUE
  return
end
subroutine CUNI1 (Z, FNU, KODE, N, Y, NZ, NLAST, FNUL, TOL, ELIM, &
     ALIM)
!
!! CUNI1 is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUNI1-A, ZUNI1-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
!     EXPANSION FOR I(FNU,Z) IN -PI/3 <= ARG Z <= PI/3.
!
!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST /= 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
!     Y(I)=CZERO FOR I=NLAST+1,N
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  CUCHK, CUNIK, CUOIK, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUNI1
  COMPLEX CFN, CONE, CRSC, CSCL, CSR, CSS, CWRK, CZERO, C1, C2, &
   PHI, RZ, SUM, S1, S2, Y, Z, ZETA1, ZETA2, CY
  REAL ALIM, APHI, ASCLE, BRY, C2I, C2M, C2R, ELIM, FN, FNU, FNUL, &
   RS1, TOL, YY, R1MACH
  INTEGER I, IFLAG, INIT, K, KODE, M, N, ND, NLAST, NN, NUF, NW, NZ
  DIMENSION BRY(3), Y(N), CWRK(16), CSS(3), CSR(3), CY(2)
  DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CUNI1
  NZ = 0
  ND = N
  NLAST = 0
!-----------------------------------------------------------------------
!     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
!     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
!     EXP(ALIM)=EXP(ELIM)*TOL
!-----------------------------------------------------------------------
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CRSC = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CRSC
  CSR(1) = CRSC
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = 1.0E+3*R1MACH(1)/TOL
!-----------------------------------------------------------------------
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!-----------------------------------------------------------------------
  FN = MAX(FNU,1.0E0)
  INIT = 0
  call CUNIK(Z, FN, 1, 1, TOL, INIT, PHI, ZETA1, ZETA2, SUM, CWRK)
  if (KODE == 1) go to 10
  CFN = CMPLX(FN,0.0E0)
  S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2))
  go to 20
   10 CONTINUE
  S1 = -ZETA1 + ZETA2
   20 CONTINUE
  RS1 = REAL(S1)
  if (ABS(RS1) > ELIM) go to 130
   30 CONTINUE
  NN = MIN(2,ND)
  DO 80 I=1,NN
    FN = FNU + (ND-I)
    INIT = 0
    call CUNIK(Z, FN, 1, 0, TOL, INIT, PHI, ZETA1, ZETA2, SUM, CWRK)
    if (KODE == 1) go to 40
    CFN = CMPLX(FN,0.0E0)
    YY = AIMAG(Z)
    S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2)) + CMPLX(0.0E0,YY)
    go to 50
   40   CONTINUE
    S1 = -ZETA1 + ZETA2
   50   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
    RS1 = REAL(S1)
    if (ABS(RS1) > ELIM) go to 110
    if (I == 1) IFLAG = 2
    if (ABS(RS1) < ALIM) go to 60
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    APHI = ABS(PHI)
    RS1 = RS1 + ALOG(APHI)
    if (ABS(RS1) > ELIM) go to 110
    if (I == 1) IFLAG = 1
    if (RS1 < 0.0E0) go to 60
    if (I == 1) IFLAG = 3
   60   CONTINUE
!-----------------------------------------------------------------------
!     SCALE S1 if ABS(S1) < ASCLE
!-----------------------------------------------------------------------
    S2 = PHI*SUM
    C2R = REAL(S1)
    C2I = AIMAG(S1)
    C2M = EXP(C2R)*REAL(CSS(IFLAG))
    S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
    S2 = S2*S1
    if (IFLAG /= 1) go to 70
    call CUCHK(S2, NW, BRY(1), TOL)
    if (NW /= 0) go to 110
   70   CONTINUE
    M = ND - I + 1
    CY(I) = S2
    Y(M) = S2*CSR(IFLAG)
   80 CONTINUE
  if (ND <= 2) go to 100
  RZ = CMPLX(2.0E0,0.0E0)/Z
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = R1MACH(2)
  S1 = CY(1)
  S2 = CY(2)
  C1 = CSR(IFLAG)
  ASCLE = BRY(IFLAG)
  K = ND - 2
  FN = K
  DO 90 I=3,ND
    C2 = S2
    S2 = S1 + CMPLX(FNU+FN,0.0E0)*RZ*S2
    S1 = C2
    C2 = S2*C1
    Y(K) = C2
    K = K - 1
    FN = FN - 1.0E0
    if (IFLAG >= 3) go to 90
    C2R = REAL(C2)
    C2I = AIMAG(C2)
    C2R = ABS(C2R)
    C2I = ABS(C2I)
    C2M = MAX(C2R,C2I)
    if (C2M <= ASCLE) go to 90
    IFLAG = IFLAG + 1
    ASCLE = BRY(IFLAG)
    S1 = S1*C1
    S2 = C2
    S1 = S1*CSS(IFLAG)
    S2 = S2*CSS(IFLAG)
    C1 = CSR(IFLAG)
   90 CONTINUE
  100 CONTINUE
  return
!-----------------------------------------------------------------------
!     SET UNDERFLOW AND UPDATE PARAMETERS
!-----------------------------------------------------------------------
  110 CONTINUE
  if (RS1 > 0.0E0) go to 120
  Y(ND) = CZERO
  NZ = NZ + 1
  ND = ND - 1
  if (ND == 0) go to 100
  call CUOIK(Z, FNU, KODE, 1, ND, Y, NUF, TOL, ELIM, ALIM)
  if (NUF < 0) go to 120
  ND = ND - NUF
  NZ = NZ + NUF
  if (ND == 0) go to 100
  FN = FNU + (ND-1)
  if (FN >= FNUL) go to 30
  NLAST = ND
  return
  120 CONTINUE
  NZ = -1
  return
  130 CONTINUE
  if (RS1 > 0.0E0) go to 120
  NZ = N
  DO 140 I=1,N
    Y(I) = CZERO
  140 CONTINUE
  return
end
subroutine CACON (Z, FNU, KODE, MR, N, Y, NZ, RL, FNUL, TOL, ELIM, &
     ALIM)
!
!! CACON is subsidiary to CBESH and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CACON-A, ZACON-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CACON APPLIES THE ANALYTIC CONTINUATION FORMULA
!
!         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
!                 MP=PI*MR*CMPLX(0.0,1.0)
!
!     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
!     HALF Z PLANE
!
!***SEE ALSO  CBESH, CBESK
!***ROUTINES CALLED  CBINU, CBKNU, CS1S2, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CACON
  COMPLEX CK, CONE, CS, CSCL, CSCR, CSGN, CSPN, CSS, CSR, C1, C2, &
   RZ, SC1, SC2, ST, S1, S2, Y, Z, ZN, CY
  REAL ALIM, ARG, ASCLE, AS2, BSCLE, BRY, CPN, C1I, C1M, C1R, ELIM, &
   FMR, FNU, FNUL, PI, RL, SGN, SPN, TOL, YY, R1MACH
  INTEGER I, INU, IUF, KFLAG, KODE, MR, N, NN, NW, NZ
  DIMENSION Y(N), CY(2), CSS(3), CSR(3), BRY(3)
  DATA PI / 3.14159265358979324E0 /
  DATA CONE / (1.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CACON
  NZ = 0
  ZN = -Z
  NN = N
  call CBINU(ZN, FNU, KODE, NN, Y, NW, RL, FNUL, TOL, ELIM, ALIM)
  if (NW < 0) go to 80
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
!-----------------------------------------------------------------------
  NN = MIN(2,N)
  call CBKNU(ZN, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
  if (NW /= 0) go to 80
  S1 = CY(1)
  FMR = MR
  SGN = -SIGN(PI,FMR)
  CSGN = CMPLX(0.0E0,SGN)
  if (KODE == 1) go to 10
  YY = -AIMAG(ZN)
  CPN = COS(YY)
  SPN = SIN(YY)
  CSGN = CSGN*CMPLX(CPN,SPN)
   10 CONTINUE
!-----------------------------------------------------------------------
!     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
  INU = FNU
  ARG = (FNU-INU)*SGN
  CPN = COS(ARG)
  SPN = SIN(ARG)
  CSPN = CMPLX(CPN,SPN)
  if (MOD(INU,2) == 1) CSPN = -CSPN
  IUF = 0
  C1 = S1
  C2 = Y(1)
  ASCLE = 1.0E+3*R1MACH(1)/TOL
  if (KODE == 1) go to 20
  call CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
  NZ = NZ + NW
  SC1 = C1
   20 CONTINUE
  Y(1) = CSPN*C1 + CSGN*C2
  if (N == 1) RETURN
  CSPN = -CSPN
  S2 = CY(2)
  C1 = S2
  C2 = Y(2)
  if (KODE == 1) go to 30
  call CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
  NZ = NZ + NW
  SC2 = C1
   30 CONTINUE
  Y(2) = CSPN*C1 + CSGN*C2
  if (N == 2) RETURN
  CSPN = -CSPN
  RZ = CMPLX(2.0E0,0.0E0)/ZN
  CK = CMPLX(FNU+1.0E0,0.0E0)*RZ
!-----------------------------------------------------------------------
!     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
!-----------------------------------------------------------------------
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CSCR = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CSCR
  CSR(1) = CSCR
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = ASCLE
  BRY(2) = 1.0E0/ASCLE
  BRY(3) = R1MACH(2)
  AS2 = ABS(S2)
  KFLAG = 2
  if (AS2 > BRY(1)) go to 40
  KFLAG = 1
  go to 50
   40 CONTINUE
  if (AS2 < BRY(2)) go to 50
  KFLAG = 3
   50 CONTINUE
  BSCLE = BRY(KFLAG)
  S1 = S1*CSS(KFLAG)
  S2 = S2*CSS(KFLAG)
  CS = CSR(KFLAG)
  DO 70 I=3,N
    ST = S2
    S2 = CK*S2 + S1
    S1 = ST
    C1 = S2*CS
    ST = C1
    C2 = Y(I)
    if (KODE == 1) go to 60
    if (IUF < 0) go to 60
    call CS1S2(ZN, C1, C2, NW, ASCLE, ALIM, IUF)
    NZ = NZ + NW
    SC1 = SC2
    SC2 = C1
    if (IUF /= 3) go to 60
    IUF = -4
    S1 = SC1*CSS(KFLAG)
    S2 = SC2*CSS(KFLAG)
    ST = SC2
   60   CONTINUE
    Y(I) = CSPN*C1 + CSGN*C2
    CK = CK + RZ
    CSPN = -CSPN
    if (KFLAG >= 3) go to 70
    C1R = REAL(C1)
    C1I = AIMAG(C1)
    C1R = ABS(C1R)
    C1I = ABS(C1I)
    C1M = MAX(C1R,C1I)
    if (C1M <= BSCLE) go to 70
    KFLAG = KFLAG + 1
    BSCLE = BRY(KFLAG)
    S1 = S1*CS
    S2 = ST
    S1 = S1*CSS(KFLAG)
    S2 = S2*CSS(KFLAG)
    CS = CSR(KFLAG)
   70 CONTINUE
  return
   80 CONTINUE
  NZ = -1
  if ( NW == (-2)) NZ=-2
  return
end
subroutine CUCHK (Y, NZ, ASCLE, TOL)
!
!! CUCHK is subsidiary to SERI, CUOIK, CUNK1, CUNK2, CUNI1, CUNI2 and CKSCL.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CUCHK-A, ZUCHK-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
!      EXP(-ALIM)=ASCLE=1.0E+3*R1MACH(1)/TOL. THE TEST IS MADE TO SEE
!      if THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDER FLOW
!      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
!      if THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
!      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
!      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
!
!***SEE ALSO  CKSCL, CUNI1, CUNI2, CUNK1, CUNK2, CUOIK, SERI
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CUCHK
!
  COMPLEX Y
  REAL ASCLE, SS, ST, TOL, YR, YI
  INTEGER NZ
!***FIRST EXECUTABLE STATEMENT  CUCHK
  NZ = 0
  YR = REAL(Y)
  YI = AIMAG(Y)
  YR = ABS(YR)
  YI = ABS(YI)
  ST = MIN(YR,YI)
  if (ST > ASCLE) RETURN
  SS = MAX(YR,YI)
  ST=ST/TOL
  if (SS < ST) NZ = 1
  return
end
subroutine CS1S2 (ZR, S1, S2, NZ, ASCLE, ALIM, IUF)
!
!! CS1S2 is subsidiary to CAIRY and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CS1S2-A, ZS1S2-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
!     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
!     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
!     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
!     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
!     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
!     PRECISION ABOVE THE UNDERFLOW LIMIT.
!
!***SEE ALSO  CAIRY, CBESK
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CS1S2
  COMPLEX CZERO, C1, S1, S1D, S2, ZR
  REAL AA, ALIM, ALN, ASCLE, AS1, AS2, XX
  INTEGER IUF, NZ
  DATA CZERO / (0.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CS1S2
  NZ = 0
  AS1 = ABS(S1)
  AS2 = ABS(S2)
  AA = REAL(S1)
  ALN = AIMAG(S1)
  if (AA == 0.0E0 .AND. ALN == 0.0E0) go to 10
  if (AS1 == 0.0E0) go to 10
  XX = REAL(ZR)
  ALN = -XX - XX + ALOG(AS1)
  S1D = S1
  S1 = CZERO
  AS1 = 0.0E0
  if (ALN < (-ALIM)) go to 10
  C1 = CLOG(S1D) - ZR - ZR
  S1 = CEXP(C1)
  AS1 = ABS(S1)
  IUF = IUF + 1
   10 CONTINUE
  AA = MAX(AS1,AS2)
  if (AA > ASCLE) RETURN
  S1 = CZERO
  S2 = CZERO
  NZ = 1
  IUF = 0
  return
end
subroutine CKSCL (ZR, FNU, N, Y, NZ, RZ, ASCLE, TOL, ELIM)
!
!! CKSCL is subsidiary to CBKNU, CUNK1 and CUNK2.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CKSCL-A, ZKSCL-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
!     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
!     return WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
!
!***SEE ALSO  CBKNU, CUNK1, CUNK2
!***ROUTINES CALLED  CUCHK
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CKSCL
  COMPLEX CK, CS, CY, CZERO, RZ, S1, S2, Y, ZR, ZD, CELM
  REAL AA, ASCLE, ACS, AS, CSI, CSR, ELIM, FN, FNU, TOL, XX, ZRI, &
   ELM, ALAS, HELIM
  INTEGER I, IC, K, KK, N, NN, NW, NZ
  DIMENSION Y(N), CY(2)
  DATA CZERO / (0.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CUCHK
  NZ = 0
  IC = 0
  XX = REAL(ZR)
  NN = MIN(2,N)
  DO 10 I=1,NN
    S1 = Y(I)
    CY(I) = S1
    AS = ABS(S1)
    ACS = -XX + ALOG(AS)
    NZ = NZ + 1
    Y(I) = CZERO
    if (ACS < (-ELIM)) go to 10
    CS = -ZR + CLOG(S1)
    CSR = REAL(CS)
    CSI = AIMAG(CS)
    AA = EXP(CSR)/TOL
    CS = CMPLX(AA,0.0E0)*CMPLX(COS(CSI),SIN(CSI))
    call CUCHK(CS, NW, ASCLE, TOL)
    if (NW /= 0) go to 10
    Y(I) = CS
    NZ = NZ - 1
    IC = I
   10 CONTINUE
  if (N == 1) RETURN
  if (IC > 1) go to 20
  Y(1) = CZERO
  NZ = 2
   20 CONTINUE
  if (N == 2) RETURN
  if (NZ == 0) RETURN
  FN = FNU + 1.0E0
  CK = CMPLX(FN,0.0E0)*RZ
  S1 = CY(1)
  S2 = CY(2)
  HELIM = 0.5E0*ELIM
  ELM = EXP(-ELIM)
  CELM = CMPLX(ELM,0.0E0)
  ZRI =AIMAG(ZR)
  ZD = ZR
!
!     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
!     S2 GETS LARGER THAN EXP(ELIM/2)
!
  DO 30 I=3,N
    KK = I
    CS = S2
    S2 = CK*S2 + S1
    S1 = CS
    CK = CK + RZ
    AS = ABS(S2)
    ALAS = ALOG(AS)
    ACS = -XX + ALAS
    NZ = NZ + 1
    Y(I) = CZERO
    if (ACS < (-ELIM)) go to 25
    CS = -ZD + CLOG(S2)
    CSR = REAL(CS)
    CSI = AIMAG(CS)
    AA = EXP(CSR)/TOL
    CS = CMPLX(AA,0.0E0)*CMPLX(COS(CSI),SIN(CSI))
    call CUCHK(CS, NW, ASCLE, TOL)
    if (NW /= 0) go to 25
    Y(I) = CS
    NZ = NZ - 1
    if (IC == (KK-1)) go to 40
    IC = KK
    go to 30
   25   CONTINUE
    if ( ALAS < HELIM) go to 30
    XX = XX-ELIM
    S1 = S1*CELM
    S2 = S2*CELM
    ZD = CMPLX(XX,ZRI)
   30 CONTINUE
  NZ = N
  if ( IC == N) NZ=N-1
  go to 45
   40 CONTINUE
  NZ = KK - 2
   45 CONTINUE
  DO 50 K=1,NZ
    Y(K) = CZERO
   50 CONTINUE
  return
end
subroutine CSHCH (Z, CSH, CCH)
!
!! CSHCH is subsidiary to CBESH and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CSHCH-A, ZSHCH-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
!     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
!
!***SEE ALSO  CBESH, CBESK
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CSHCH
  COMPLEX CCH, CSH, Z
  REAL CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, X, Y
!***FIRST EXECUTABLE STATEMENT  CSHCH
  X = REAL(Z)
  Y = AIMAG(Z)
  SH = SINH(X)
  CH = COSH(X)
  SN = SIN(Y)
  CN = COS(Y)
  CSHR = SH*CN
  CSHI = CH*SN
  CSH = CMPLX(CSHR,CSHI)
  CCHR = CH*CN
  CCHI = SH*SN
  CCH = CMPLX(CCHR,CCHI)
  return
end
subroutine CSERI (Z, FNU, KODE, N, Y, NZ, TOL, ELIM, ALIM)
!
!! CSERI is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CSERI-A, ZSERI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z) >= 0.0 BY
!     MEANS OF THE POWER SERIES FOR LARGE ABS(Z) IN THE
!     REGION ABS(Z) <= 2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
!     NZ > 0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
!     DUE TO UNDERFLOW. NZ < 0 MEANS UNDERFLOW OCCURRED, BUT THE
!     CONDITION ABS(Z) <= 2*SQRT(FNU+1) WAS VIOLATED AND THE
!     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  CUCHK, GAMLN, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CSERI
  COMPLEX AK1, CK, COEF, CONE, CRSC, CZ, CZERO, HZ, RZ, S1, S2, W, &
   Y, Z
  REAL AA, ACZ, AK, ALIM, ARM, ASCLE, ATOL, AZ, DFNU, ELIM, FNU, &
   FNUP, RAK1, RS, RTR1, S, SS, TOL, X, GAMLN, R1MACH
  INTEGER I, IB, IDUM, IFLAG, IL, K, KODE, L, M, N, NN, NW, NZ
  DIMENSION Y(N), W(2)
  DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CSERI
  NZ = 0
  AZ = ABS(Z)
  if (AZ == 0.0E0) go to 150
  X = REAL(Z)
  ARM = 1.0E+3*R1MACH(1)
  RTR1 = SQRT(ARM)
  CRSC = CMPLX(1.0E0,0.0E0)
  IFLAG = 0
  if (AZ < ARM) go to 140
  HZ = Z*CMPLX(0.5E0,0.0E0)
  CZ = CZERO
  if (AZ > RTR1) CZ = HZ*HZ
  ACZ = ABS(CZ)
  NN = N
  CK = CLOG(HZ)
   10 CONTINUE
  DFNU = FNU + (NN-1)
  FNUP = DFNU + 1.0E0
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
  AK1 = CK*CMPLX(DFNU,0.0E0)
  AK = GAMLN(FNUP,IDUM)
  AK1 = AK1 - CMPLX(AK,0.0E0)
  if (KODE == 2) AK1 = AK1 - CMPLX(X,0.0E0)
  RAK1 = REAL(AK1)
  if (RAK1 > (-ELIM)) go to 30
   20 CONTINUE
  NZ = NZ + 1
  Y(NN) = CZERO
  if (ACZ > DFNU) go to 170
  NN = NN - 1
  if (NN == 0) RETURN
  go to 10
   30 CONTINUE
  if (RAK1 > (-ALIM)) go to 40
  IFLAG = 1
  SS = 1.0E0/TOL
  CRSC = CMPLX(TOL,0.0E0)
  ASCLE = ARM*SS
   40 CONTINUE
  AK = AIMAG(AK1)
  AA = EXP(RAK1)
  if (IFLAG == 1) AA = AA*SS
  COEF = CMPLX(AA,0.0E0)*CMPLX(COS(AK),SIN(AK))
  ATOL = TOL*ACZ/FNUP
  IL = MIN(2,NN)
  DO 80 I=1,IL
    DFNU = FNU + (NN-I)
    FNUP = DFNU + 1.0E0
    S1 = CONE
    if (ACZ < TOL*FNUP) go to 60
    AK1 = CONE
    AK = FNUP + 2.0E0
    S = FNUP
    AA = 2.0E0
   50   CONTINUE
    RS = 1.0E0/S
    AK1 = AK1*CZ*CMPLX(RS,0.0E0)
    S1 = S1 + AK1
    S = S + AK
    AK = AK + 2.0E0
    AA = AA*ACZ*RS
    if (AA > ATOL) go to 50
   60   CONTINUE
    M = NN - I + 1
    S2 = S1*COEF
    W(I) = S2
    if (IFLAG == 0) go to 70
    call CUCHK(S2, NW, ASCLE, TOL)
    if (NW /= 0) go to 20
   70   CONTINUE
    Y(M) = S2*CRSC
    if (I /= IL) COEF = COEF*CMPLX(DFNU,0.0E0)/HZ
   80 CONTINUE
  if (NN <= 2) RETURN
  K = NN - 2
  AK = K
  RZ = (CONE+CONE)/Z
  if (IFLAG == 1) go to 110
  IB = 3
   90 CONTINUE
  DO 100 I=IB,NN
    Y(K) = CMPLX(AK+FNU,0.0E0)*RZ*Y(K+1) + Y(K+2)
    AK = AK - 1.0E0
    K = K - 1
  100 CONTINUE
  return
!-----------------------------------------------------------------------
!     RECUR BACKWARD WITH SCALED VALUES
!-----------------------------------------------------------------------
  110 CONTINUE
!-----------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
!     UNDERFLOW LIMIT = ASCLE = R1MACH(1)*CSCL*1.0E+3
!-----------------------------------------------------------------------
  S1 = W(1)
  S2 = W(2)
  DO 120 L=3,NN
    CK = S2
    S2 = S1 + CMPLX(AK+FNU,0.0E0)*RZ*S2
    S1 = CK
    CK = S2*CRSC
    Y(K) = CK
    AK = AK - 1.0E0
    K = K - 1
    if (ABS(CK) > ASCLE) go to 130
  120 CONTINUE
  return
  130 CONTINUE
  IB = L + 1
  if (IB > NN) RETURN
  go to 90
  140 CONTINUE
  NZ = N
  if (FNU == 0.0E0) NZ = NZ - 1
  150 CONTINUE
  Y(1) = CZERO
  if (FNU == 0.0E0) Y(1) = CONE
  if (N == 1) RETURN
  DO 160 I=2,N
    Y(I) = CZERO
  160 CONTINUE
  return
!-----------------------------------------------------------------------
!     return WITH NZ < 0 if ABS(Z*Z/4) > FNU+N-NZ-1 COMPLETE
!     THE CALCULATION IN CBINU WITH N=N-ABS(NZ)
!-----------------------------------------------------------------------
  170 CONTINUE
  NZ = -NZ
  return
end
subroutine CBUNK (Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
!
!! CBUNK is subsidiary to CBESH and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CBUNK-A, ZBUNK-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU > FNUL.
!     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
!     IN CUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN CUNK2
!
!***SEE ALSO  CBESH, CBESK
!***ROUTINES CALLED  CUNK1, CUNK2
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CBUNK
  COMPLEX Y, Z
  REAL ALIM, AX, AY, ELIM, FNU, TOL, XX, YY
  INTEGER KODE, MR, N, NZ
  DIMENSION Y(N)
!***FIRST EXECUTABLE STATEMENT  CBUNK
  NZ = 0
  XX = REAL(Z)
  YY = AIMAG(Z)
  AX = ABS(XX)*1.7321E0
  AY = ABS(YY)
  if (AY > AX) go to 10
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3 <= ARG(Z) <= PI/3
!-----------------------------------------------------------------------
  call CUNK1(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
  go to 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
!     AND HPI=PI/2
!-----------------------------------------------------------------------
  call CUNK2(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
   20 CONTINUE
  return
end
subroutine CBINU (Z, FNU, KODE, N, CY, NZ, RL, FNUL, TOL, ELIM, ALIM)
!
!! CBINU is subsidiary to CAIRY, CBESH, CBESI, CBESJ, CBESK and CBIRY.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CBINU-A, ZBINU-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
!
!***SEE ALSO  CAIRY, CBESH, CBESI, CBESJ, CBESK, CBIRY
!***ROUTINES CALLED  CASYI, CBUNI, CMLRI, CSERI, CUOIK, CWRSK
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CBINU
  COMPLEX CW, CY, CZERO, Z
  REAL ALIM, AZ, DFNU, ELIM, FNU, FNUL, RL, TOL
  INTEGER I, INW, KODE, N, NLAST, NN, NUI, NW, NZ
  DIMENSION CY(N), CW(2)
  DATA CZERO / (0.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CBINU
  NZ = 0
  AZ = ABS(Z)
  NN = N
  DFNU = FNU + (N-1)
  if (AZ <= 2.0E0) go to 10
  if (AZ*AZ*0.25E0 > DFNU+1.0E0) go to 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     POWER SERIES
!-----------------------------------------------------------------------
  call CSERI(Z, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
  INW = ABS(NW)
  NZ = NZ + INW
  NN = NN - INW
  if (NN == 0) RETURN
  if (NW >= 0) go to 120
  DFNU = FNU + (NN-1)
   20 CONTINUE
  if (AZ < RL) go to 40
  if (DFNU <= 1.0E0) go to 30
  if (AZ+AZ < DFNU*DFNU) go to 50
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR LARGE Z
!-----------------------------------------------------------------------
   30 CONTINUE
  call CASYI(Z, FNU, KODE, NN, CY, NW, RL, TOL, ELIM, ALIM)
  if (NW < 0) go to 130
  go to 120
   40 CONTINUE
  if (DFNU <= 1.0E0) go to 70
   50 CONTINUE
!-----------------------------------------------------------------------
!     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
!-----------------------------------------------------------------------
  call CUOIK(Z, FNU, KODE, 1, NN, CY, NW, TOL, ELIM, ALIM)
  if (NW < 0) go to 130
  NZ = NZ + NW
  NN = NN - NW
  if (NN == 0) RETURN
  DFNU = FNU+(NN-1)
  if (DFNU > FNUL) go to 110
  if (AZ > FNUL) go to 110
   60 CONTINUE
  if (AZ > RL) go to 80
   70 CONTINUE
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE SERIES
!-----------------------------------------------------------------------
  call CMLRI(Z, FNU, KODE, NN, CY, NW, TOL)
  if ( NW < 0) go to 130
  go to 120
   80 CONTINUE
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
!-----------------------------------------------------------------------
  call CUOIK(Z, FNU, KODE, 2, 2, CW, NW, TOL, ELIM, ALIM)
  if (NW >= 0) go to 100
  NZ = NN
  DO 90 I=1,NN
    CY(I) = CZERO
   90 CONTINUE
  return
  100 CONTINUE
  if (NW > 0) go to 130
  call CWRSK(Z, FNU, KODE, NN, CY, NW, CW, TOL, ELIM, ALIM)
  if (NW < 0) go to 130
  go to 120
  110 CONTINUE
!-----------------------------------------------------------------------
!     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
!-----------------------------------------------------------------------
  NUI = FNUL-DFNU + 1
  NUI = MAX(NUI,0)
  call CBUNI(Z, FNU, KODE, NN, CY, NW, NUI, NLAST, FNUL, TOL, ELIM, &
   ALIM)
  if (NW < 0) go to 130
  NZ = NZ + NW
  if (NLAST == 0) go to 120
  NN = NLAST
  go to 60
  120 CONTINUE
  return
  130 CONTINUE
  NZ = -1
  if ( NW == (-2)) NZ=-2
  return
end
subroutine CBUNI (Z, FNU, KODE, N, Y, NZ, NUI, NLAST, FNUL, TOL, &
     ELIM, ALIM)
!
!! CBUNI is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CBUNI-A, ZBUNI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE ABS(Z) >
!     FNUL AND FNU+N-1 < FNUL. THE ORDER IS INCREASED FROM
!     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
!     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
!     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  CUNI1, CUNI2, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CBUNI
  COMPLEX CSCL, CSCR, CY, RZ, ST, S1, S2, Y, Z
  REAL ALIM, AX, AY, DFNU, ELIM, FNU, FNUI, FNUL, GNU, TOL, XX, YY, &
   ASCLE, BRY, STR, STI, STM, R1MACH
  INTEGER I, IFLAG, IFORM, K, KODE, N, NL, NLAST, NUI, NW, NZ
  DIMENSION Y(N), CY(2), BRY(3)
!***FIRST EXECUTABLE STATEMENT  CBUNI
  NZ = 0
  XX = REAL(Z)
  YY = AIMAG(Z)
  AX = ABS(XX)*1.7321E0
  AY = ABS(YY)
  IFORM = 1
  if (AY > AX) IFORM = 2
  if (NUI == 0) go to 60
  FNUI = NUI
  DFNU = FNU + (N-1)
  GNU = DFNU + FNUI
  if (IFORM == 2) go to 10
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3 <= ARG(Z) <= PI/3
!-----------------------------------------------------------------------
  call CUNI1(Z, GNU, KODE, 2, CY, NW, NLAST, FNUL, TOL, ELIM, ALIM)
  go to 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
!     AND HPI=PI/2
!-----------------------------------------------------------------------
  call CUNI2(Z, GNU, KODE, 2, CY, NW, NLAST, FNUL, TOL, ELIM, ALIM)
   20 CONTINUE
  if (NW < 0) go to 50
  if (NW /= 0) go to 90
  AY = ABS(CY(1))
!----------------------------------------------------------------------
!     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
!----------------------------------------------------------------------
  BRY(1) = 1.0E+3*R1MACH(1)/TOL
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = BRY(2)
  IFLAG = 2
  ASCLE = BRY(2)
  AX = 1.0E0
  CSCL = CMPLX(AX,0.0E0)
  if (AY > BRY(1)) go to 21
  IFLAG = 1
  ASCLE = BRY(1)
  AX = 1.0E0/TOL
  CSCL = CMPLX(AX,0.0E0)
  go to 25
   21 CONTINUE
  if (AY < BRY(2)) go to 25
  IFLAG = 3
  ASCLE = BRY(3)
  AX = TOL
  CSCL = CMPLX(AX,0.0E0)
   25 CONTINUE
  AY = 1.0E0/AX
  CSCR = CMPLX(AY,0.0E0)
  S1 = CY(2)*CSCL
  S2 = CY(1)*CSCL
  RZ = CMPLX(2.0E0,0.0E0)/Z
  DO 30 I=1,NUI
    ST = S2
    S2 = CMPLX(DFNU+FNUI,0.0E0)*RZ*S2 + S1
    S1 = ST
    FNUI = FNUI - 1.0E0
    if (IFLAG >= 3) go to 30
    ST = S2*CSCR
    STR = REAL(ST)
    STI = AIMAG(ST)
    STR = ABS(STR)
    STI = ABS(STI)
    STM = MAX(STR,STI)
    if (STM <= ASCLE) go to 30
    IFLAG = IFLAG+1
    ASCLE = BRY(IFLAG)
    S1 = S1*CSCR
    S2 = ST
    AX = AX*TOL
    AY = 1.0E0/AX
    CSCL = CMPLX(AX,0.0E0)
    CSCR = CMPLX(AY,0.0E0)
    S1 = S1*CSCL
    S2 = S2*CSCL
   30 CONTINUE
  Y(N) = S2*CSCR
  if (N == 1) RETURN
  NL = N - 1
  FNUI = NL
  K = NL
  DO 40 I=1,NL
    ST = S2
    S2 = CMPLX(FNU+FNUI,0.0E0)*RZ*S2 + S1
    S1 = ST
    ST = S2*CSCR
    Y(K) = ST
    FNUI = FNUI - 1.0E0
    K = K - 1
    if (IFLAG >= 3) go to 40
    STR = REAL(ST)
    STI = AIMAG(ST)
    STR = ABS(STR)
    STI = ABS(STI)
    STM = MAX(STR,STI)
    if (STM <= ASCLE) go to 40
    IFLAG = IFLAG+1
    ASCLE = BRY(IFLAG)
    S1 = S1*CSCR
    S2 = ST
    AX = AX*TOL
    AY = 1.0E0/AX
    CSCL = CMPLX(AX,0.0E0)
    CSCR = CMPLX(AY,0.0E0)
    S1 = S1*CSCL
    S2 = S2*CSCL
   40 CONTINUE
  return
   50 CONTINUE
  NZ = -1
  if ( NW == (-2)) NZ=-2
  return
   60 CONTINUE
  if (IFORM == 2) go to 70
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3 <= ARG(Z) <= PI/3
!-----------------------------------------------------------------------
  call CUNI1(Z, FNU, KODE, N, Y, NW, NLAST, FNUL, TOL, ELIM, ALIM)
  go to 80
   70 CONTINUE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
!     AND HPI=PI/2
!-----------------------------------------------------------------------
  call CUNI2(Z, FNU, KODE, N, Y, NW, NLAST, FNUL, TOL, ELIM, ALIM)
   80 CONTINUE
  if (NW < 0) go to 50
  NZ = NW
  return
   90 CONTINUE
  NLAST = N
  return
end
subroutine CASYI (Z, FNU, KODE, N, Y, NZ, RL, TOL, ELIM, ALIM)
!
!! CASYI is subsidiary to CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CASYI-A, ZASYI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z) >= 0.0 BY
!     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE ABS(Z) IN THE
!     REGION ABS(Z) > MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
!     NZ < 0 INDICATES AN OVERFLOW ON KODE=1.
!
!***SEE ALSO  CBESI, CBESK
!***ROUTINES CALLED  R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CASYI
  COMPLEX AK1, CK, CONE, CS1, CS2, CZ, CZERO, DK, EZ, P1, RZ, S2, &
   Y, Z
  REAL AA, ACZ, AEZ, AK, ALIM, ARG, ARM, ATOL, AZ, BB, BK, DFNU, &
   DNU2, ELIM, FDN, FNU, PI, RL, RTPI, RTR1, S, SGN, SQK, TOL, X, &
   YY, R1MACH
  INTEGER I, IB, IL, INU, J, JL, K, KODE, KODED, M, N, NN, NZ
  DIMENSION Y(N)
  DATA PI, RTPI  /3.14159265358979324E0 , 0.159154943091895336E0 /
  DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CASYI
  NZ = 0
  AZ = ABS(Z)
  X = REAL(Z)
  ARM = 1.0E+3*R1MACH(1)
  RTR1 = SQRT(ARM)
  IL = MIN(2,N)
  DFNU = FNU + (N-IL)
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
  AK1 = CMPLX(RTPI,0.0E0)/Z
  AK1 = CSQRT(AK1)
  CZ = Z
  if (KODE == 2) CZ = Z - CMPLX(X,0.0E0)
  ACZ = REAL(CZ)
  if (ABS(ACZ) > ELIM) go to 80
  DNU2 = DFNU + DFNU
  KODED = 1
  if ((ABS(ACZ) > ALIM) .AND. (N > 2)) go to 10
  KODED = 0
  AK1 = AK1*CEXP(CZ)
   10 CONTINUE
  FDN = 0.0E0
  if (DNU2 > RTR1) FDN = DNU2*DNU2
  EZ = Z*CMPLX(8.0E0,0.0E0)
!-----------------------------------------------------------------------
!     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
!     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
!     EXPANSION FOR THE IMAGINARY PART.
!-----------------------------------------------------------------------
  AEZ = 8.0E0*AZ
  S = TOL/AEZ
  JL = RL+RL + 2
  YY = AIMAG(Z)
  P1 = CZERO
  if (YY == 0.0E0) go to 20
!-----------------------------------------------------------------------
!     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
!     SIGNIFICANCE WHEN FNU OR N IS LARGE
!-----------------------------------------------------------------------
  INU = FNU
  ARG = (FNU-INU)*PI
  INU = INU + N - IL
  AK = -SIN(ARG)
  BK = COS(ARG)
  if (YY < 0.0E0) BK = -BK
  P1 = CMPLX(AK,BK)
  if (MOD(INU,2) == 1) P1 = -P1
   20 CONTINUE
  DO 50 K=1,IL
    SQK = FDN - 1.0E0
    ATOL = S*ABS(SQK)
    SGN = 1.0E0
    CS1 = CONE
    CS2 = CONE
    CK = CONE
    AK = 0.0E0
    AA = 1.0E0
    BB = AEZ
    DK = EZ
    DO 30 J=1,JL
      CK = CK*CMPLX(SQK,0.0E0)/DK
      CS2 = CS2 + CK
      SGN = -SGN
      CS1 = CS1 + CK*CMPLX(SGN,0.0E0)
      DK = DK + EZ
      AA = AA*ABS(SQK)/BB
      BB = BB + AEZ
      AK = AK + 8.0E0
      SQK = SQK - AK
      if (AA <= ATOL) go to 40
   30   CONTINUE
    go to 90
   40   CONTINUE
    S2 = CS1
    if (X+X < ELIM) S2 = S2 + P1*CS2*CEXP(-Z-Z)
    FDN = FDN + 8.0E0*DFNU + 4.0E0
    P1 = -P1
    M = N - IL + K
    Y(M) = S2*AK1
   50 CONTINUE
  if (N <= 2) RETURN
  NN = N
  K = NN - 2
  AK = K
  RZ = (CONE+CONE)/Z
  IB = 3
  DO 60 I=IB,NN
    Y(K) = CMPLX(AK+FNU,0.0E0)*RZ*Y(K+1) + Y(K+2)
    AK = AK - 1.0E0
    K = K - 1
   60 CONTINUE
  if (KODED == 0) RETURN
  CK = CEXP(CZ)
  DO 70 I=1,NN
    Y(I) = Y(I)*CK
   70 CONTINUE
  return
   80 CONTINUE
  NZ = -1
  return
   90 CONTINUE
  NZ=-2
  return
end
subroutine CBKNU (Z, FNU, KODE, N, Y, NZ, TOL, ELIM, ALIM)
!
!! CBKNU is subsidiary to CAIRY, CBESH, CBESI and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CBKNU-A, ZBKNU-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE
!
!***SEE ALSO  CAIRY, CBESH, CBESI, CBESK
!***ROUTINES CALLED  CKSCL, CSHCH, CUCHK, GAMLN, I1MACH, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CBKNU
!
  COMPLEX CCH, CK, COEF, CONE, CRSC, CS, CSCL, CSH, CSR, CSS, CTWO, &
   CZ, CZERO, F, FMU, P, PT, P1, P2, Q, RZ, SMU, ST, S1, S2, Y, Z, &
   ZD, CELM, CY
  REAL AA, AK, ALIM, ASCLE, A1, A2, BB, BK, BRY, CAZ, CC, DNU, &
   DNU2, ELIM, ETEST, FC, FHS, FK, FKS, FNU, FPI, G1, G2, HPI, PI, &
   P2I, P2M, P2R, RK, RTHPI, R1, S, SPI, TM, TOL, TTH, T1, T2, XX, &
   YY, GAMLN, R1MACH, HELIM, ELM, XD, YD, ALAS, AS
  INTEGER I, IDUM, IFLAG, INU, K, KFLAG, KK, KMAX, KODE, KODED, N, &
   NZ, I1MACH, NW, J, IC, INUB
  DIMENSION BRY(3), CC(8), CSS(3), CSR(3), Y(N), CY(2)
!
  DATA KMAX / 30 /
  DATA R1 / 2.0E0 /
  DATA CZERO,CONE,CTWO /(0.0E0,0.0E0),(1.0E0,0.0E0),(2.0E0,0.0E0)/
!
  DATA PI, RTHPI, SPI ,HPI, FPI, TTH / &
       3.14159265358979324E0,       1.25331413731550025E0, &
       1.90985931710274403E0,       1.57079632679489662E0, &
       1.89769999331517738E0,       6.66666666666666666E-01/
!
  DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)/ &
       5.77215664901532861E-01,    -4.20026350340952355E-02, &
      -4.21977345555443367E-02,     7.21894324666309954E-03, &
      -2.15241674114950973E-04,    -2.01348547807882387E-05, &
       1.13302723198169588E-06,     6.11609510448141582E-09/
!
!***FIRST EXECUTABLE STATEMENT  CBKNU
  XX = REAL(Z)
  YY = AIMAG(Z)
  CAZ = ABS(Z)
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CRSC = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CRSC
  CSR(1) = CRSC
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = 1.0E+3*R1MACH(1)/TOL
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = R1MACH(2)
  NZ = 0
  IFLAG = 0
  KODED = KODE
  RZ = CTWO/Z
  INU = FNU+0.5E0
  DNU = FNU - INU
  if (ABS(DNU) == 0.5E0) go to 110
  DNU2 = 0.0E0
  if (ABS(DNU) > TOL) DNU2 = DNU*DNU
  if (CAZ > R1) go to 110
!-----------------------------------------------------------------------
!     SERIES FOR ABS(Z) <= R1
!-----------------------------------------------------------------------
  FC = 1.0E0
  SMU = CLOG(RZ)
  FMU = SMU*CMPLX(DNU,0.0E0)
  call CSHCH(FMU, CSH, CCH)
  if (DNU == 0.0E0) go to 10
  FC = DNU*PI
  FC = FC/SIN(FC)
  SMU = CSH*CMPLX(1.0E0/DNU,0.0E0)
   10 CONTINUE
  A2 = 1.0E0 + DNU
!-----------------------------------------------------------------------
!     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
!-----------------------------------------------------------------------
  T2 = EXP(-GAMLN(A2,IDUM))
  T1 = 1.0E0/(T2*FC)
  if (ABS(DNU) > 0.1E0) go to 40
!-----------------------------------------------------------------------
!     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
!-----------------------------------------------------------------------
  AK = 1.0E0
  S = CC(1)
  DO 20 K=2,8
    AK = AK*DNU2
    TM = CC(K)*AK
    S = S + TM
    if (ABS(TM) < TOL) go to 30
   20 CONTINUE
   30 G1 = -S
  go to 50
   40 CONTINUE
  G1 = (T1-T2)/(DNU+DNU)
   50 CONTINUE
  G2 = 0.5E0*(T1+T2)*FC
  G1 = G1*FC
  F = CMPLX(G1,0.0E0)*CCH + SMU*CMPLX(G2,0.0E0)
  PT = CEXP(FMU)
  P = CMPLX(0.5E0/T2,0.0E0)*PT
  Q = CMPLX(0.5E0/T1,0.0E0)/PT
  S1 = F
  S2 = P
  AK = 1.0E0
  A1 = 1.0E0
  CK = CONE
  BK = 1.0E0 - DNU2
  if (INU > 0 .OR. N > 1) go to 80
!-----------------------------------------------------------------------
!     GENERATE K(FNU,Z), 0.0D0  <=  FNU  <  0.5D0 AND N=1
!-----------------------------------------------------------------------
  if (CAZ < TOL) go to 70
  CZ = Z*Z*CMPLX(0.25E0,0.0E0)
  T1 = 0.25E0*CAZ*CAZ
   60 CONTINUE
  F = (F*CMPLX(AK,0.0E0)+P+Q)*CMPLX(1.0E0/BK,0.0E0)
  P = P*CMPLX(1.0E0/(AK-DNU),0.0E0)
  Q = Q*CMPLX(1.0E0/(AK+DNU),0.0E0)
  RK = 1.0E0/AK
  CK = CK*CZ*CMPLX(RK,0.0)
  S1 = S1 + CK*F
  A1 = A1*T1*RK
  BK = BK + AK + AK + 1.0E0
  AK = AK + 1.0E0
  if (A1 > TOL) go to 60
   70 CONTINUE
  Y(1) = S1
  if (KODED == 1) RETURN
  Y(1) = S1*CEXP(Z)
  return
!-----------------------------------------------------------------------
!     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
!-----------------------------------------------------------------------
   80 CONTINUE
  if (CAZ < TOL) go to 100
  CZ = Z*Z*CMPLX(0.25E0,0.0E0)
  T1 = 0.25E0*CAZ*CAZ
   90 CONTINUE
  F = (F*CMPLX(AK,0.0E0)+P+Q)*CMPLX(1.0E0/BK,0.0E0)
  P = P*CMPLX(1.0E0/(AK-DNU),0.0E0)
  Q = Q*CMPLX(1.0E0/(AK+DNU),0.0E0)
  RK = 1.0E0/AK
  CK = CK*CZ*CMPLX(RK,0.0E0)
  S1 = S1 + CK*F
  S2 = S2 + CK*(P-F*CMPLX(AK,0.0E0))
  A1 = A1*T1*RK
  BK = BK + AK + AK + 1.0E0
  AK = AK + 1.0E0
  if (A1 > TOL) go to 90
  100 CONTINUE
  KFLAG = 2
  BK = REAL(SMU)
  A1 = FNU + 1.0E0
  AK = A1*ABS(BK)
  if (AK > ALIM) KFLAG = 3
  P2 = S2*CSS(KFLAG)
  S2 = P2*RZ
  S1 = S1*CSS(KFLAG)
  if (KODED == 1) go to 210
  F = CEXP(Z)
  S1 = S1*F
  S2 = S2*F
  go to 210
!-----------------------------------------------------------------------
!     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
!     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
!     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
!     RECURSION
!-----------------------------------------------------------------------
  110 CONTINUE
  COEF = CMPLX(RTHPI,0.0E0)/CSQRT(Z)
  KFLAG = 2
  if (KODED == 2) go to 120
  if (XX > ALIM) go to 290
!     BLANK LINE
  A1 = EXP(-XX)*REAL(CSS(KFLAG))
  PT = CMPLX(A1,0.0E0)*CMPLX(COS(YY),-SIN(YY))
  COEF = COEF*PT
  120 CONTINUE
  if (ABS(DNU) == 0.5E0) go to 300
!-----------------------------------------------------------------------
!     MILLER ALGORITHM FOR ABS(Z) > R1
!-----------------------------------------------------------------------
  AK = COS(PI*DNU)
  AK = ABS(AK)
  if (AK == 0.0E0) go to 300
  FHS = ABS(0.25E0-DNU2)
  if (FHS == 0.0E0) go to 300
!-----------------------------------------------------------------------
!     COMPUTE R2=F(E). if ABS(Z) >= R2, USE FORWARD RECURRENCE TO
!     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
!     12 <= E <= 60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(11))=
!     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
!-----------------------------------------------------------------------
  T1 = (I1MACH(11)-1)*R1MACH(5)*3.321928094E0
  T1 = MAX(T1,12.0E0)
  T1 = MIN(T1,60.0E0)
  T2 = TTH*T1 - 6.0E0
  if (XX /= 0.0E0) go to 130
  T1 = HPI
  go to 140
  130 CONTINUE
  T1 = ATAN(YY/XX)
  T1 = ABS(T1)
  140 CONTINUE
  if (T2 > CAZ) go to 170
!-----------------------------------------------------------------------
!     FORWARD RECURRENCE LOOP WHEN ABS(Z) >= R2
!-----------------------------------------------------------------------
  ETEST = AK/(PI*CAZ*TOL)
  FK = 1.0E0
  if (ETEST < 1.0E0) go to 180
  FKS = 2.0E0
  RK = CAZ + CAZ + 2.0E0
  A1 = 0.0E0
  A2 = 1.0E0
  DO 150 I=1,KMAX
    AK = FHS/FKS
    BK = RK/(FK+1.0E0)
    TM = A2
    A2 = BK*A2 - AK*A1
    A1 = TM
    RK = RK + 2.0E0
    FKS = FKS + FK + FK + 2.0E0
    FHS = FHS + FK + FK
    FK = FK + 1.0E0
    TM = ABS(A2)*FK
    if (ETEST < TM) go to 160
  150 CONTINUE
  go to 310
  160 CONTINUE
  FK = FK + SPI*T1*SQRT(T2/CAZ)
  FHS = ABS(0.25E0-DNU2)
  go to 180
  170 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE BACKWARD INDEX K FOR ABS(Z) < R2
!-----------------------------------------------------------------------
  A2 = SQRT(CAZ)
  AK = FPI*AK/(TOL*SQRT(A2))
  AA = 3.0E0*T1/(1.0E0+CAZ)
  BB = 14.7E0*T1/(28.0E0+CAZ)
  AK = (ALOG(AK)+CAZ*COS(AA)/(1.0E0+0.008E0*CAZ))/COS(BB)
  FK = 0.12125E0*AK*AK/CAZ + 1.5E0
  180 CONTINUE
  K = FK
!-----------------------------------------------------------------------
!     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
!-----------------------------------------------------------------------
  FK = K
  FKS = FK*FK
  P1 = CZERO
  P2 = CMPLX(TOL,0.0E0)
  CS = P2
  DO 190 I=1,K
    A1 = FKS - FK
    A2 = (FKS+FK)/(A1+FHS)
    RK = 2.0E0/(FK+1.0E0)
    T1 = (FK+XX)*RK
    T2 = YY*RK
    PT = P2
    P2 = (P2*CMPLX(T1,T2)-P1)*CMPLX(A2,0.0E0)
    P1 = PT
    CS = CS + P2
    FKS = A1 - FK + 1.0E0
    FK = FK - 1.0E0
  190 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE (P2/CS)=(P2/ABS(CS))*(CONJG(CS)/ABS(CS)) FOR BETTER
!     SCALING
!-----------------------------------------------------------------------
  TM = ABS(CS)
  PT = CMPLX(1.0E0/TM,0.0E0)
  S1 = PT*P2
  CS = CONJG(CS)*PT
  S1 = COEF*S1*CS
  if (INU > 0 .OR. N > 1) go to 200
  ZD = Z
  if ( IFLAG == 1) go to 270
  go to 240
  200 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE P1/P2=(P1/ABS(P2)*CONJG(P2)/ABS(P2) FOR SCALING
!-----------------------------------------------------------------------
  TM = ABS(P2)
  PT = CMPLX(1.0E0/TM,0.0E0)
  P1 = PT*P1
  P2 = CONJG(P2)*PT
  PT = P1*P2
  S2 = S1*(CONE+(CMPLX(DNU+0.5E0,0.0E0)-PT)/Z)
!-----------------------------------------------------------------------
!     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION WITH
!     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
!-----------------------------------------------------------------------
  210 CONTINUE
  CK = CMPLX(DNU+1.0E0,0.0E0)*RZ
  if (N == 1) INU = INU - 1
  if (INU > 0) go to 220
  if (N == 1) S1=S2
  ZD = Z
  if ( IFLAG == 1) go to 270
  go to 240
  220 CONTINUE
  INUB = 1
  if (IFLAG == 1) go to 261
  225 CONTINUE
  P1 = CSR(KFLAG)
  ASCLE = BRY(KFLAG)
  DO 230 I=INUB,INU
    ST = S2
    S2 = CK*S2 + S1
    S1 = ST
    CK = CK + RZ
    if (KFLAG >= 3) go to 230
    P2 = S2*P1
    P2R = REAL(P2)
    P2I = AIMAG(P2)
    P2R = ABS(P2R)
    P2I = ABS(P2I)
    P2M = MAX(P2R,P2I)
    if (P2M <= ASCLE) go to 230
    KFLAG = KFLAG + 1
    ASCLE = BRY(KFLAG)
    S1 = S1*P1
    S2 = P2
    S1 = S1*CSS(KFLAG)
    S2 = S2*CSS(KFLAG)
    P1 = CSR(KFLAG)
  230 CONTINUE
  if (N == 1) S1 = S2
  240 CONTINUE
  Y(1) = S1*CSR(KFLAG)
  if (N == 1) RETURN
  Y(2) = S2*CSR(KFLAG)
  if (N == 2) RETURN
  KK = 2
  250 CONTINUE
  KK = KK + 1
  if (KK > N) RETURN
  P1 = CSR(KFLAG)
  ASCLE = BRY(KFLAG)
  DO 260 I=KK,N
    P2 = S2
    S2 = CK*S2 + S1
    S1 = P2
    CK = CK + RZ
    P2 = S2*P1
    Y(I) = P2
    if (KFLAG >= 3) go to 260
    P2R = REAL(P2)
    P2I = AIMAG(P2)
    P2R = ABS(P2R)
    P2I = ABS(P2I)
    P2M = MAX(P2R,P2I)
    if (P2M <= ASCLE) go to 260
    KFLAG = KFLAG + 1
    ASCLE = BRY(KFLAG)
    S1 = S1*P1
    S2 = P2
    S1 = S1*CSS(KFLAG)
    S2 = S2*CSS(KFLAG)
    P1 = CSR(KFLAG)
  260 CONTINUE
  return
!-----------------------------------------------------------------------
!     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
!-----------------------------------------------------------------------
  261 CONTINUE
  HELIM = 0.5E0*ELIM
  ELM = EXP(-ELIM)
  CELM = CMPLX(ELM,0.0)
  ASCLE = BRY(1)
  ZD = Z
  XD = XX
  YD = YY
  IC = -1
  J = 2
  DO 262 I=1,INU
    ST = S2
    S2 = CK*S2+S1
    S1 = ST
    CK = CK+RZ
    AS = ABS(S2)
    ALAS = ALOG(AS)
    P2R = -XD+ALAS
    if ( P2R < (-ELIM)) go to 263
    P2 = -ZD+CLOG(S2)
    P2R = REAL(P2)
    P2I = AIMAG(P2)
    P2M = EXP(P2R)/TOL
    P1 = CMPLX(P2M,0.0E0)*CMPLX(COS(P2I),SIN(P2I))
    call CUCHK(P1,NW,ASCLE,TOL)
    if ( NW /= 0) go to 263
    J=3-J
    CY(J) = P1
    if ( IC == (I-1)) go to 264
    IC = I
    go to 262
  263   CONTINUE
    if ( ALAS < HELIM) go to 262
    XD = XD-ELIM
    S1 = S1*CELM
    S2 = S2*CELM
    ZD = CMPLX(XD,YD)
  262 CONTINUE
  if ( N == 1) S1 = S2
  go to 270
  264 CONTINUE
  KFLAG = 1
  INUB = I+1
  S2 = CY(J)
  J = 3 - J
  S1 = CY(J)
  if ( INUB <= INU) go to 225
  if ( N == 1) S1 = S2
  go to 240
  270 CONTINUE
  Y(1) = S1
  if (N == 1) go to 280
  Y(2) = S2
  280 CONTINUE
  ASCLE = BRY(1)
  call CKSCL(ZD, FNU, N, Y, NZ, RZ, ASCLE, TOL, ELIM)
  INU = N - NZ
  if (INU <= 0) RETURN
  KK = NZ + 1
  S1 = Y(KK)
  Y(KK) = S1*CSR(1)
  if (INU == 1) RETURN
  KK = NZ + 2
  S2 = Y(KK)
  Y(KK) = S2*CSR(1)
  if (INU == 2) RETURN
  T2 = FNU + (KK-1)
  CK = CMPLX(T2,0.0E0)*RZ
  KFLAG = 1
  go to 250
  290 CONTINUE
!-----------------------------------------------------------------------
!     SCALE BY EXP(Z), IFLAG = 1 CASES
!-----------------------------------------------------------------------
  KODED = 2
  IFLAG = 1
  KFLAG = 2
  go to 120
!-----------------------------------------------------------------------
!     FNU=HALF ODD INTEGER CASE, DNU=-0.5
!-----------------------------------------------------------------------
  300 CONTINUE
  S1 = COEF
  S2 = COEF
  go to 210
  310 CONTINUE
  NZ=-2
  return
end
FUNCTION I1MACH (I)
!
!! I1MACH returns integer machine dependent constants.
!
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      INTEGER (I1MACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   I1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument and can be referenced as follows:
!
!        K = I1MACH(I)
!
!   where I=1,...,16.  The (output) value of K above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   I/O unit numbers:
!     I1MACH( 1) = the standard input unit.
!     I1MACH( 2) = the standard output unit.
!     I1MACH( 3) = the standard punch unit.
!     I1MACH( 4) = the standard error message unit.
!
!   Words:
!     I1MACH( 5) = the number of bits per integer storage unit.
!     I1MACH( 6) = the number of characters per integer storage unit.
!
!   Integers:
!     assume integers are represented in the S-digit, base-A form
!
!                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!                where 0  <=  X(I)  <  A for I=0,...,S-1.
!     I1MACH( 7) = A, the base.
!     I1MACH( 8) = S, the number of base-A digits.
!     I1MACH( 9) = A**S - 1, the largest magnitude.
!
!   Floating-Point Numbers:
!     Assume floating-point numbers are represented in the T-digit,
!     base-B form
!                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!                where 0  <=  X(I)  <  B for I=1,...,T,
!                0  <  X(1), and EMIN  <=  E  <=  EMAX.
!     I1MACH(10) = B, the base.
!
!   Single-Precision:
!     I1MACH(11) = T, the number of base-B digits.
!     I1MACH(12) = EMIN, the smallest exponent E.
!     I1MACH(13) = EMAX, the largest exponent E.
!
!   Double-Precision:
!     I1MACH(14) = T, the number of base-B digits.
!     I1MACH(15) = EMIN, the smallest exponent E.
!     I1MACH(16) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   891012  Added VAX G-floating constants.  (WRB)
!   891012  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900618  Added DEC RISC constants.  (WRB)
!   900723  Added IBM RS 6000 constants.  (WRB)
!   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
!           (RWC)
!   910710  Added HP 730 constants.  (SMR)
!   911114  Added Convex IEEE constants.  (WRB)
!   920121  Added SUN -r8 compiler option constants.  (WRB)
!   920229  Added Touchstone Delta i860 constants.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920625  Added Convex -p8 and -pd8 compiler option constants.
!           (BKS, WRB)
!   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
!           options.  (DWL, RWC and WRB).
!***END PROLOGUE  I1MACH
!
  integer i1mach
  INTEGER IMACH(16),OUTPUT
  SAVE IMACH
  EQUIVALENCE (IMACH(4),OUTPUT)
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT COMPILER
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1022 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE APOLLO
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        129 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1025 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA IMACH( 1) /          7 /
!     DATA IMACH( 2) /          2 /
!     DATA IMACH( 3) /          2 /
!     DATA IMACH( 4) /          2 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         33 /
!     DATA IMACH( 9) / Z1FFFFFFFF /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -256 /
!     DATA IMACH(13) /        255 /
!     DATA IMACH(14) /         60 /
!     DATA IMACH(15) /       -256 /
!     DATA IMACH(16) /        255 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         48 /
!     DATA IMACH( 6) /          6 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         39 /
!     DATA IMACH( 9) / O0007777777777777 /
!     DATA IMACH(10) /          8 /
!     DATA IMACH(11) /         13 /
!     DATA IMACH(12) /        -50 /
!     DATA IMACH(13) /         76 /
!     DATA IMACH(14) /         26 /
!     DATA IMACH(15) /        -50 /
!     DATA IMACH(16) /         76 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         48 /
!     DATA IMACH( 6) /          6 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         39 /
!     DATA IMACH( 9) / O0007777777777777 /
!     DATA IMACH(10) /          8 /
!     DATA IMACH(11) /         13 /
!     DATA IMACH(12) /        -50 /
!     DATA IMACH(13) /         76 /
!     DATA IMACH(14) /         26 /
!     DATA IMACH(15) /     -32754 /
!     DATA IMACH(16) /      32780 /
!
!     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          8 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /      -4095 /
!     DATA IMACH(13) /       4094 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /      -4095 /
!     DATA IMACH(16) /       4094 /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /    6LOUTPUT/
!     DATA IMACH( 5) /         60 /
!     DATA IMACH( 6) /         10 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         48 /
!     DATA IMACH( 9) / 00007777777777777777B /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /       -929 /
!     DATA IMACH(13) /       1070 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /       -929 /
!     DATA IMACH(16) /       1069 /
!
!     MACHINE CONSTANTS FOR THE CELERITY C1260
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          0 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / Z'7FFFFFFF' /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1022 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fn COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fi COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -p8 COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         53 /
!     DATA IMACH(12) /      -1023 /
!     DATA IMACH(13) /       1023 /
!     DATA IMACH(14) /        113 /
!     DATA IMACH(15) /     -16383 /
!     DATA IMACH(16) /      16383 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -pd8 COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         53 /
!     DATA IMACH(12) /      -1023 /
!     DATA IMACH(13) /       1023 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE CRAY
!     USING THE 46 BIT INTEGER COMPILER OPTION
!
!     DATA IMACH( 1) /        100 /
!     DATA IMACH( 2) /        101 /
!     DATA IMACH( 3) /        102 /
!     DATA IMACH( 4) /        101 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          8 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         46 /
!     DATA IMACH( 9) / 1777777777777777B /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /      -8189 /
!     DATA IMACH(13) /       8190 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /      -8099 /
!     DATA IMACH(16) /       8190 /
!
!     MACHINE CONSTANTS FOR THE CRAY
!     USING THE 64 BIT INTEGER COMPILER OPTION
!
!     DATA IMACH( 1) /        100 /
!     DATA IMACH( 2) /        101 /
!     DATA IMACH( 3) /        102 /
!     DATA IMACH( 4) /        101 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          8 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 777777777777777777777B /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /      -8189 /
!     DATA IMACH(13) /       8190 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /      -8099 /
!     DATA IMACH(16) /       8190 /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!
!     DATA IMACH( 1) /         11 /
!     DATA IMACH( 2) /         12 /
!     DATA IMACH( 3) /          8 /
!     DATA IMACH( 4) /         10 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /         16 /
!     DATA IMACH(11) /          6 /
!     DATA IMACH(12) /        -64 /
!     DATA IMACH(13) /         63 /
!     DATA IMACH(14) /         14 /
!     DATA IMACH(15) /        -64 /
!     DATA IMACH(16) /         63 /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING G_FLOAT
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING IEEE_FLOAT
!
      DATA IMACH( 1) /          5 /
      DATA IMACH( 2) /          6 /
      DATA IMACH( 3) /          6 /
      DATA IMACH( 4) /          6 /
      DATA IMACH( 5) /         32 /
      DATA IMACH( 6) /          4 /
      DATA IMACH( 7) /          2 /
      DATA IMACH( 8) /         31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /          2 /
      DATA IMACH(11) /         24 /
      DATA IMACH(12) /       -125 /
      DATA IMACH(13) /        128 /
      DATA IMACH(14) /         53 /
      DATA IMACH(15) /      -1021 /
      DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE DEC RISC
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING D_FLOATING
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING G_FLOATING
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE ELXSI 6400
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         32 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1022 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          0 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         24 /
!     DATA IMACH( 6) /          3 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         23 /
!     DATA IMACH( 9) /    8388607 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         23 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         38 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /         43 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          6 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / O377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         63 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HP 730
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     3 WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          4 /
!     DATA IMACH( 4) /          1 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         23 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         39 /
!     DATA IMACH(15) /       -128 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     4 WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          4 /
!     DATA IMACH( 4) /          1 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         23 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         55 /
!     DATA IMACH(15) /       -128 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HP 9000
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          7 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         32 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1015 /
!     DATA IMACH(16) /       1017 /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) /  Z7FFFFFFF /
!     DATA IMACH(10) /         16 /
!     DATA IMACH(11) /          6 /
!     DATA IMACH(12) /        -64 /
!     DATA IMACH(13) /         63 /
!     DATA IMACH(14) /         14 /
!     DATA IMACH(15) /        -64 /
!     DATA IMACH(16) /         63 /
!
!     MACHINE CONSTANTS FOR THE IBM PC
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          0 /
!     DATA IMACH( 4) /          0 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE IBM RS 6000
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          0 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE INTEL i860
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          5 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         54 /
!     DATA IMACH(15) /       -101 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          5 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         62 /
!     DATA IMACH(15) /       -128 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGER ARITHMETIC.
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGER ARITHMETIC.
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE SUN
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE SUN
!     USING THE -r8 COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         53 /
!     DATA IMACH(12) /      -1021 /
!     DATA IMACH(13) /       1024 /
!     DATA IMACH(14) /        113 /
!     DATA IMACH(15) /     -16381 /
!     DATA IMACH(16) /      16384 /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          1 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / O377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         60 /
!     DATA IMACH(15) /      -1024 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
!
!     DATA IMACH( 1) /          1 /
!     DATA IMACH( 2) /          1 /
!     DATA IMACH( 3) /          0 /
!     DATA IMACH( 4) /          1 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!***FIRST EXECUTABLE STATEMENT  I1MACH
!
  if ( I < 1 .OR. I > 16 ) then
    WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
    STOP
  end if

  I1MACH = IMACH(I)

  return
end
FUNCTION R1MACH (I)
!
!! R1MACH returns floating point machine dependent constants.
!
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      SINGLE PRECISION (R1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   R1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        A = R1MACH(I)
!
!   where I=1,...,5.  The (output) value of A above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   R1MACH(3) = B**(-T), the smallest relative spacing.
!   R1MACH(4) = B**(1-T), the largest relative spacing.
!   R1MACH(5) = LOG10(B)
!
!   Assume single precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0  <=  X(I)  <  B for I=1,...,T, 0  <  X(1), and
!   EMIN  <=  E  <=  EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of R1MACH(1) - R1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   890213  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900618  Added DEC RISC constants.  (WRB)
!   900723  Added IBM RS 6000 constants.  (WRB)
!   910710  Added HP 730 constants.  (SMR)
!   911114  Added Convex IEEE constants.  (WRB)
!   920121  Added SUN -r8 compiler option constants.  (WRB)
!   920229  Added Touchstone Delta i860 constants.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920625  Added CONVEX -p8 and -pd8 compiler option constants.
!           (BKS, WRB)
!   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!***END PROLOGUE  R1MACH
!
  real r1mach
  INTEGER SMALL(2)
  INTEGER LARGE(2)
  INTEGER RIGHT(2)
  INTEGER DIVER(2)
  INTEGER LOG10(2)
!
  REAL RMACH(5)
  SAVE RMACH
!
  EQUIVALENCE (RMACH(1),SMALL(1))
  EQUIVALENCE (RMACH(2),LARGE(1))
  EQUIVALENCE (RMACH(3),RIGHT(1))
  EQUIVALENCE (RMACH(4),DIVER(1))
  EQUIVALENCE (RMACH(5),LOG10(1))
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
!
!     DATA SMALL(1) / Z'00800000' /
!     DATA LARGE(1) / Z'7F7FFFFF' /
!     DATA RIGHT(1) / Z'33800000' /
!     DATA DIVER(1) / Z'34000000' /
!     DATA LOG10(1) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
!
!     DATA SMALL(1) / Z'00800000' /
!     DATA LARGE(1) / Z'7EFFFFFF' /
!     DATA RIGHT(1) / Z'33800000' /
!     DATA DIVER(1) / Z'34000000' /
!     DATA LOG10(1) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE APOLLO
!
!     DATA SMALL(1) / 16#00800000 /
!     DATA LARGE(1) / 16#7FFFFFFF /
!     DATA RIGHT(1) / 16#33800000 /
!     DATA DIVER(1) / 16#34000000 /
!     DATA LOG10(1) / 16#3E9A209B /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA RMACH(1) / Z400800000 /
!     DATA RMACH(2) / Z5FFFFFFFF /
!     DATA RMACH(3) / Z4E9800000 /
!     DATA RMACH(4) / Z4EA800000 /
!     DATA RMACH(5) / Z500E730E8 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS
!
!     DATA RMACH(1) / O1771000000000000 /
!     DATA RMACH(2) / O0777777777777777 /
!     DATA RMACH(3) / O1311000000000000 /
!     DATA RMACH(4) / O1301000000000000 /
!     DATA RMACH(5) / O1157163034761675 /
!
!     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
!
!     DATA RMACH(1) / Z"3001800000000000" /
!     DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" /
!     DATA RMACH(3) / Z"3FD2800000000000" /
!     DATA RMACH(4) / Z"3FD3800000000000" /
!     DATA RMACH(5) / Z"3FFF9A209A84FBCF" /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!
!     DATA RMACH(1) / 00564000000000000000B /
!     DATA RMACH(2) / 37767777777777777776B /
!     DATA RMACH(3) / 16414000000000000000B /
!     DATA RMACH(4) / 16424000000000000000B /
!     DATA RMACH(5) / 17164642023241175720B /
!
!     MACHINE CONSTANTS FOR THE CELERITY C1260
!
!     DATA SMALL(1) / Z'00800000' /
!     DATA LARGE(1) / Z'7F7FFFFF' /
!     DATA RIGHT(1) / Z'33800000' /
!     DATA DIVER(1) / Z'34000000' /
!     DATA LOG10(1) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fn COMPILER OPTION
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7FFFFFFF' /
!     DATA RMACH(3) / Z'34800000' /
!     DATA RMACH(4) / Z'35000000' /
!     DATA RMACH(5) / Z'3F9A209B' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fi COMPILER OPTION
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -p8 OR -pd8 COMPILER OPTION
!
!     DATA RMACH(1) / Z'0010000000000000' /
!     DATA RMACH(2) / Z'7FFFFFFFFFFFFFFF' /
!     DATA RMACH(3) / Z'3CC0000000000000' /
!     DATA RMACH(4) / Z'3CD0000000000000' /
!     DATA RMACH(5) / Z'3FF34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE CRAY
!
!     DATA RMACH(1) / 200034000000000000000B /
!     DATA RMACH(2) / 577767777777777777776B /
!     DATA RMACH(3) / 377224000000000000000B /
!     DATA RMACH(4) / 377234000000000000000B /
!     DATA RMACH(5) / 377774642023241175720B /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
!     STATIC RMACH(5)
!
!     DATA SMALL /    20K,       0 /
!     DATA LARGE / 77777K, 177777K /
!     DATA RIGHT / 35420K,       0 /
!     DATA DIVER / 36020K,       0 /
!     DATA LOG10 / 40423K,  42023K /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING G_FLOAT
!
!     DATA RMACH(1) / '00000080'X /
!     DATA RMACH(2) / 'FFFF7FFF'X /
!     DATA RMACH(3) / '00003480'X /
!     DATA RMACH(4) / '00003500'X /
!     DATA RMACH(5) / '209B3F9A'X /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING IEEE_FLOAT
!
      DATA RMACH(1) / Z'00800000' /
      DATA RMACH(2) / Z'7F7FFFFF' /
      DATA RMACH(3) / Z'33800000' /
      DATA RMACH(4) / Z'34000000' /
      DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE DEC RISC
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
!     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
!
!     DATA SMALL(1) /       128 /
!     DATA LARGE(1) /    -32769 /
!     DATA RIGHT(1) /     13440 /
!     DATA DIVER(1) /     13568 /
!     DATA LOG10(1) / 547045274 /
!
!     DATA SMALL(1) / Z00000080 /
!     DATA LARGE(1) / ZFFFF7FFF /
!     DATA RIGHT(1) / Z00003480 /
!     DATA DIVER(1) / Z00003500 /
!     DATA LOG10(1) / Z209B3F9A /
!
!     MACHINE CONSTANTS FOR THE ELXSI 6400
!     (ASSUMING REAL*4 IS THE DEFAULT REAL)
!
!     DATA SMALL(1) / '00800000'X /
!     DATA LARGE(1) / '7F7FFFFF'X /
!     DATA RIGHT(1) / '33800000'X /
!     DATA DIVER(1) / '34000000'X /
!     DATA LOG10(1) / '3E9A209B'X /
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
!     DATA LARGE(1), LARGE(2) / '37777777, '00000177 /
!     DATA RIGHT(1), RIGHT(2) / '20000000, '00000352 /
!     DATA DIVER(1), DIVER(2) / '20000000, '00000353 /
!     DATA LOG10(1), LOG10(2) / '23210115, '00000377 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
!
!     DATA RMACH(1) / O402400000000 /
!     DATA RMACH(2) / O376777777777 /
!     DATA RMACH(3) / O714400000000 /
!     DATA RMACH(4) / O716400000000 /
!     DATA RMACH(5) / O776464202324 /
!
!     MACHINE CONSTANTS FOR THE HP 730
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     3 WORD DOUBLE PRECISION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) / 40000B,       1 /
!     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
!     DATA DIVER(1), DIVER(2) / 40000B,    327B /
!     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     4 WORD DOUBLE PRECISION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) / 40000B,       1 /
!     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
!     DATA DIVER(1), DIVER(2) / 40000B,    327B /
!     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
!
!     MACHINE CONSTANTS FOR THE HP 9000
!
!     DATA SMALL(1) / 00004000000B /
!     DATA LARGE(1) / 17677777777B /
!     DATA RIGHT(1) / 06340000000B /
!     DATA DIVER(1) / 06400000000B /
!     DATA LOG10(1) / 07646420233B /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA RMACH(1) / Z00100000 /
!     DATA RMACH(2) / Z7FFFFFFF /
!     DATA RMACH(3) / Z3B100000 /
!     DATA RMACH(4) / Z3C100000 /
!     DATA RMACH(5) / Z41134413 /
!
!     MACHINE CONSTANTS FOR THE IBM PC
!
!     DATA SMALL(1) / 1.18E-38      /
!     DATA LARGE(1) / 3.40E+38      /
!     DATA RIGHT(1) / 0.595E-07     /
!     DATA DIVER(1) / 1.19E-07      /
!     DATA LOG10(1) / 0.30102999566 /
!
!     MACHINE CONSTANTS FOR THE IBM RS 6000
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE INTEL i860
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR)
!
!     DATA RMACH(1) / "000400000000 /
!     DATA RMACH(2) / "377777777777 /
!     DATA RMACH(3) / "146400000000 /
!     DATA RMACH(4) / "147400000000 /
!     DATA RMACH(5) / "177464202324 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1) /    8388608 /
!     DATA LARGE(1) / 2147483647 /
!     DATA RIGHT(1) /  880803840 /
!     DATA DIVER(1) /  889192448 /
!     DATA LOG10(1) / 1067065499 /
!
!     DATA RMACH(1) / O00040000000 /
!     DATA RMACH(2) / O17777777777 /
!     DATA RMACH(3) / O06440000000 /
!     DATA RMACH(4) / O06500000000 /
!     DATA RMACH(5) / O07746420233 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1), SMALL(2) /   128,     0 /
!     DATA LARGE(1), LARGE(2) / 32767,    -1 /
!     DATA RIGHT(1), RIGHT(2) / 13440,     0 /
!     DATA DIVER(1), DIVER(2) / 13568,     0 /
!     DATA LOG10(1), LOG10(2) / 16282,  8347 /
!
!     DATA SMALL(1), SMALL(2) / O000200, O000000 /
!     DATA LARGE(1), LARGE(2) / O077777, O177777 /
!     DATA RIGHT(1), RIGHT(2) / O032200, O000000 /
!     DATA DIVER(1), DIVER(2) / O032400, O000000 /
!     DATA LOG10(1), LOG10(2) / O037632, O020233 /
!
!     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE SUN
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE SUN
!     USING THE -r8 COMPILER OPTION
!
!     DATA RMACH(1) / Z'0010000000000000' /
!     DATA RMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA RMACH(3) / Z'3CA0000000000000' /
!     DATA RMACH(4) / Z'3CB0000000000000' /
!     DATA RMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES
!
!     DATA RMACH(1) / O000400000000 /
!     DATA RMACH(2) / O377777777777 /
!     DATA RMACH(3) / O146400000000 /
!     DATA RMACH(4) / O147400000000 /
!     DATA RMACH(5) / O177464202324 /
!
!     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
!
!     DATA SMALL(1), SMALL(2) /     0,    256/
!     DATA LARGE(1), LARGE(2) /    -1,   -129/
!     DATA RIGHT(1), RIGHT(2) /     0,  26880/
!     DATA DIVER(1), DIVER(2) /     0,  27136/
!     DATA LOG10(1), LOG10(2) /  8347,  32538/
!
!***FIRST EXECUTABLE STATEMENT  R1MACH
!
  if ( I < 1 .OR. I > 5 ) then
    call XERMSG ('SLATEC', 'R1MACH', 'I OUT OF BOUNDS', 1, 2)
  end if

  R1MACH = RMACH(I)

  return
end
subroutine XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)

!! XERMSG processes XERROR messages.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERMSG-A)
!***KEYWORDS  ERROR MESSAGE, XERROR
!***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!***DESCRIPTION
!
!   XERMSG processes a diagnostic message in a manner determined by the
!   value of LEVEL and the current value of the library error control
!   flag, KONTRL.  See subroutine XSETF for details.
!
!    LIBRAR   A character constant (or character variable) with the name
!             of the library.  This will be 'SLATEC' for the SLATEC
!             Common Math Library.  The error handling package is
!             general enough to be used by many libraries
!             simultaneously, so it is desirable for the routine that
!             detects and reports an error to identify the library name
!             as well as the routine name.
!
!    SUBROU   A character constant (or character variable) with the name
!             of the routine that detected the error.  Usually it is the
!             name of the routine that is calling XERMSG.  There are
!             some instances where a user callable library routine calls
!             lower level subsidiary routines where the error is
!             detected.  In such cases it may be more informative to
!             supply the name of the routine the user called rather than
!             the name of the subsidiary routine that detected the
!             error.
!
!    MESSG    A character constant (or character variable) with the text
!             of the error or warning message.  In the example below,
!             the message is a character constant that contains a
!             generic message.
!
!                   call XERMSG ('SLATEC', 'MMPY',
!                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
!                  *3, 1)
!
!             It is possible (and is sometimes desirable) to generate a
!             specific message--e.g., one that contains actual numeric
!             values.  Specific numeric values can be converted into
!             character strings using formatted WRITE statements into
!             character variables.  This is called standard Fortran
!             internal file I/O and is exemplified in the first three
!             lines of the following example.  You can also catenate
!             substrings of characters to construct the error message.
!             Here is an example showing the use of both writing to
!             an internal file and catenating character strings.
!
!                   CHARACTER*5 CHARN, CHARL
!                   WRITE (CHARN,10) N
!                   WRITE (CHARL,10) LDA
!                10 FORMAT(I5)
!                   call XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
!                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
!                  *   CHARL, 3, 1)
!
!             There are two subtleties worth mentioning.  One is that
!             the // for character catenation is used to construct the
!             error message so that no single character constant is
!             continued to the next line.  This avoids confusion as to
!             whether there are trailing blanks at the end of the line.
!             The second is that by catenating the parts of the message
!             as an actual argument rather than encoding the entire
!             message into one large character variable, we avoid
!             having to know how long the message will be in order to
!             declare an adequate length for that large character
!             variable.  XERMSG calls XERPRN to print the message using
!             multiple lines if necessary.  If the message is very long,
!             XERPRN will break it into pieces of 72 characters (as
!             requested by XERMSG) for printing on multiple lines.
!             Also, XERMSG asks XERPRN to prefix each line with ' *  '
!             so that the total line length could be 76 characters.
!             Note also that XERPRN scans the error message backwards
!             to ignore trailing blanks.  Another feature is that
!             the substring '$$' is treated as a new line sentinel
!             by XERPRN.  If you want to construct a multiline
!             message without having to count out multiples of 72
!             characters, just use '$$' as a separator.  '$$'
!             obviously must occur within 72 characters of the
!             start of each line to have its intended effect since
!             XERPRN is asked to wrap around at 72 characters in
!             addition to looking for '$$'.
!
!    NERR     An integer value that is chosen by the library routine's
!             author.  It must be in the range -99 to 999 (three
!             printable digits).  Each distinct error should have its
!             own error number.  These error numbers should be described
!             in the machine readable documentation for the routine.
!             The error numbers need be unique only within each routine,
!             so it is reasonable for each routine to start enumerating
!             errors from 1 and proceeding to the next integer.
!
!    LEVEL    An integer value in the range 0 to 2 that indicates the
!             level (severity) of the error.  Their meanings are
!
!            -1  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.  An attempt is made to only print this
!                message once.
!
!             0  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.
!
!             1  A recoverable error.  This is used even if the error is
!                so serious that the routine cannot return any useful
!                answer.  If the user has told the error package to
!                return after recoverable errors, then XERMSG will
!                return to the Library routine which can then return to
!                the user's routine.  The user may also permit the error
!                package to terminate the program upon encountering a
!                recoverable error.
!
!             2  A fatal error.  XERMSG will not return to its caller
!                after it receives a fatal error.  This level should
!                hardly ever be used; it is much better to allow the
!                user a chance to recover.  An example of one of the few
!                cases in which it is permissible to declare a level 2
!                error is a reverse communication Library routine that
!                is likely to be called repeatedly until it integrates
!                across some interval.  If there is a serious error in
!                the input such that another step cannot be taken and
!                the Library routine is called again without the input
!                error having been corrected by the caller, the Library
!                routine will probably be called forever with improper
!                input.  In this case, it is reasonable to declare the
!                error to be fatal.
!
!    Each of the arguments to XERMSG is input; none will be modified by
!    XERMSG.  A routine may make multiple calls to XERMSG with warning
!    level messages; however, after a call to XERMSG with a recoverable
!    error, the routine should return to the user.  Do not try to call
!    XERMSG with a second recoverable error after the first recoverable
!    error because the error package saves the error number.  The user
!    can retrieve this error number by calling another entry point in
!    the error handling package and then clear the error number when
!    recovering from the error.  Calling XERMSG in succession causes the
!    old error number to be overwritten by the latest error number.
!    This is considered harmless for error numbers associated with
!    warning messages but must not be done for error numbers of serious
!    errors.  After a call to XERMSG with a recoverable error, the user
!    must be given a chance to call NUMXER or XERCLR to retrieve or
!    clear the error number.
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
!***REVISION HISTORY  (YYMMDD)
!   880101  DATE WRITTEN
!   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
!           THERE ARE TWO BASIC CHANGES.
!           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
!               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
!               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
!               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
!               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
!               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
!               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
!               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
!           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
!               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
!               OF LOWER CASE.
!   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
!           THE PRINCIPAL CHANGES ARE
!           1.  CLARIFY COMMENTS IN THE PROLOGUES
!           2.  RENAME XRPRNT TO XERPRN
!           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
!               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
!               CHARACTER FOR NEW RECORDS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           CLEAN UP THE CODING.
!   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
!           PREFIX.
!   891013  REVISED TO CORRECT COMMENTS.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
!           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
!           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
!           XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERMSG
  CHARACTER*(*) LIBRAR, SUBROU, MESSG
  CHARACTER*8 XLIBR, XSUBR
  CHARACTER*72  TEMP
  CHARACTER*20  LFIRST
!***FIRST EXECUTABLE STATEMENT  XERMSG
  LKNTRL = J4SAVE (2, 0, .FALSE.)
  MAXMES = J4SAVE (4, 0, .FALSE.)
!
!       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
!       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
!          SHOULD BE PRINTED.
!
!       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
!          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
!          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
!
  if (NERR < -9999999 .OR. NERR > 99999999 .OR. NERR == 0 .OR. &
     LEVEL < -1 .OR. LEVEL > 2) THEN
     call XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' // &
        'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '// &
        'JOB ABORT DUE TO FATAL ERROR.', 72)
     call XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
     call XERHLT (' ***XERMSG -- INVALID INPUT')
     return
  end if
!
!       RECORD THE MESSAGE.
!
  I = J4SAVE (1, NERR, .TRUE.)
  call XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
!
!       HANDLE PRINT-ONCE WARNING MESSAGES.
!
  if (LEVEL == -1 .AND. KOUNT > 1) RETURN
!
!       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
!
  XLIBR  = LIBRAR
  XSUBR  = SUBROU
  LFIRST = MESSG
  LERR   = NERR
  LLEVEL = LEVEL
  call XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
!
  LKNTRL = MAX(-2, MIN(2,LKNTRL))
  MKNTRL = ABS(LKNTRL)
!
!       SKIP PRINTING if THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
!       ZERO AND THE ERROR IS NOT FATAL.
!
  if (LEVEL < 2 .AND. LKNTRL == 0) go to 30
  if (LEVEL == 0 .AND. KOUNT > MAXMES) go to 30
  if (LEVEL == 1 .AND. KOUNT > MAXMES .AND. MKNTRL == 1) go to 30
  if (LEVEL == 2 .AND. KOUNT > MAX(1,MAXMES)) go to 30
!
!       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
!       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
!       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY if CONTROL FLAG
!       IS NOT ZERO.
!
  if (LKNTRL  /=  0) THEN
     TEMP(1:21) = 'MESSAGE FROM ROUTINE '
     I = MIN(LEN(SUBROU), 16)
     TEMP(22:21+I) = SUBROU(1:I)
     TEMP(22+I:33+I) = ' IN LIBRARY '
     LTEMP = 33 + I
     I = MIN(LEN(LIBRAR), 16)
     TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
     TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
     LTEMP = LTEMP + I + 1
     call XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
  end if
!
!       if LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
!       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
!       FROM EACH OF THE FOLLOWING THREE OPTIONS.
!       1.  LEVEL OF THE MESSAGE
!              'INFORMATIVE MESSAGE'
!              'POTENTIALLY RECOVERABLE ERROR'
!              'FATAL ERROR'
!       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
!              'PROG CONTINUES'
!              'PROG ABORTED'
!       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
!           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
!           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
!              'TRACEBACK REQUESTED'
!              'TRACEBACK NOT REQUESTED'
!       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
!       EXCEED 74 CHARACTERS.
!       WE SKIP THE NEXT BLOCK if THE INTRODUCTORY LINE IS NOT NEEDED.
!
  if (LKNTRL  >  0) THEN
!
!       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
!
     if (LEVEL  <=  0) THEN
        TEMP(1:20) = 'INFORMATIVE MESSAGE,'
        LTEMP = 20
     ELSEIF (LEVEL  ==  1) THEN
        TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
        LTEMP = 30
     ELSE
        TEMP(1:12) = 'FATAL ERROR,'
        LTEMP = 12
     ENDIF
!
!       THEN WHETHER THE PROGRAM WILL CONTINUE.
!
     if ((MKNTRL == 2 .AND. LEVEL >= 1) .OR. &
         (MKNTRL == 1 .AND. LEVEL == 2)) THEN
        TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
        LTEMP = LTEMP + 14
     ELSE
        TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
        LTEMP = LTEMP + 16
     ENDIF
!
!       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
!
     if (LKNTRL  >  0) THEN
        TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
        LTEMP = LTEMP + 20
     ELSE
        TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
        LTEMP = LTEMP + 24
     ENDIF
     call XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
  end if
!
!       NOW SEND OUT THE MESSAGE.
!
  call XERPRN (' *  ', -1, MESSG, 72)
!
!       if LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
!          TRACEBACK.
!
  if (LKNTRL  >  0) THEN
     WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
     DO 10 I=16,22
        if (TEMP(I:I)  /=  ' ') go to 20
   10    CONTINUE
!
   20    call XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
     call FDUMP
  end if
!
!       if LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
!
  if (LKNTRL  /=  0) THEN
     call XERPRN (' *  ', -1, ' ', 72)
     call XERPRN (' ***', -1, 'END OF MESSAGE', 72)
     call XERPRN ('    ',  0, ' ', 72)
  end if
!
!       if THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
!       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
!
   30 if (LEVEL <= 0 .OR. (LEVEL == 1 .AND. MKNTRL <= 1)) RETURN
!
!       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
!       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
!       SUMMARY if THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
!
  if (LKNTRL > 0 .AND. KOUNT < MAX(1,MAXMES)) THEN
     if (LEVEL  ==  1) THEN
        call XERPRN &
           (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
     ELSE
        call XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
     ENDIF
     call XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
     call XERHLT (' ')
  ELSE
     call XERHLT (MESSG)
  end if
  return
end
subroutine fdump ( )

!*****************************************************************************80
!
!! FDUMP produces a symbolic dump.
!
!  Discussion:
!
!    This routine is intended to be replaced by a locally written
!    version which produces a symbolic dump.  Failing this,
!    it should be replaced by a version which prints the
!    subprogram nesting list.
!
!    Normally, the dump information should be printed to all the
!    active error output units.  The number and value of these
!    units can be determined by calling XGETUA.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    None
!
  implicit none

  return
end
function j4save ( which, value, set )

!*****************************************************************************80
!
!! J4SAVE saves variables needed by the library error handling routines.
!
!  Discussion:
!
!    The internal parameters are initialized to the following values:
!
!    #1 =  0, NERR, the index of the most recent error;
!    #2 =  0, KONTRL, error control flag (0 means only level 2 errors are fatal,
!             and get a printout, while lower level errors get no printout.)
!    #3 =  0, IUNIT, the main error output unit (0 means use standard output).
!    #4 = 10, MAXMES, the maximum number of times any message is printed.
!    #5 =  1, NUNIT, total number of error output units in use.
!    #6 = -1, second error output unit (-1 means not being used).
!    #7 = -1, third error output unit (-1 means not being used).
!    #8 = -1, fourth error output unit (-1 means not being used).
!    #9 = -1, fifth error output unit (-1 means not being used).
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) WHICH, the index of the item desired.
!    1, NERR, the current error number.
!    2, KONTRL, the current error control flag.
!    3, IUNIT, the current unit number to which error messages are sent.
!       (0 means use standard.)
!    4, MAXMES, the maximum times any message is printed (as set by xermax).
!    5, NUNIT, the number of units to which each error message is written.
!    6, the 2nd unit for error messages.
!    7, the 3rd unit for error messages.
!    8, the 4th unit for error messages.
!    9, the 5th unit for error messages.
!
!    Input, integer ( kind = 4 ) VALUE, the value to be set for the WHICH-th 
!    parameter, if SET is TRUE.
!
!    Input, logical SET.
!    TRUE: the WHICH-th parameter will be given the value, VALUE.
!
!    Output, integer ( kind = 4 ) J4SAVE, the old value of the WHICH-th
!    parameter.
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ), save, dimension ( 9 ) :: param = (/  &
    0, 2, 0, 10, 1, -1, -1, -1, -1 /)
  logical set
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which

  if ( which < 1 .or. 9 < which ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'J4SAVE - Fatal error!'
    write ( *, '(a,i10)' ) '  Illegal input value of WHICH = ', which
    stop
  end if

  j4save = param(which)

  if ( set ) then
    param(which) = value
  end if

  return
end
subroutine xercnt ( librar, subrou, messg, nerr, level, kontrl )

!*****************************************************************************80
!
!! XERCNT allows user control over the handling of errors.
!
!  Description:
!
!    This routine allows user control over handling of individual errors.
!
!    This routine is to be used when the error message routine XERMSG
!    is employed.  The similar routine XERCTL is to be used for the
!    older error message routines XERROR and XERRWV.
!
!    Just after each message is recorded, but before it is
!    processed any further (i.e., before it is printed or
!    a decision to abort is made), a call is made to XERCNT.
!
!    If the user has replaced this default, dummy version of XERCNT
!    with a customized routine, it can then be used to override the
!    value of KONTROL used in processing this message by redefining its value.
!
!    KONTRL may be set to any value from -2 to 2.
!
!    The meanings for KONTRL are the same as in XSETF, except
!    that the value of KONTRL changes only for this message.
!
!    If KONTRL is set to a value outside the range from -2 to 2,
!    it will be moved back into that range.
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) LIBRAR, the library or software package
!    from which the error message is coming.
!
!    Input, character ( len = * ) SUBROU, the subroutine or function within
!    the library, from which the error message is coming.
!
!    Input, character ( len = * ) MESSG, the error message.
!
!    Input, integer ( kind = 4 ) NERR, the error number.
!
!    Input, integer ( kind = 4 ) LEVEL, the error severity level.
!    * 2, this is an unconditionally fatal error.
!    * 1, this is a recoverable error.  It is normally non-fatal, unless
!         KONTRL has been reset by XSETF.
!    * 0, this is a warning message only.
!    *-1, this is a warning message which is to be printed at most once,
!         regardless of how many times this call is executed.
!
!    Input/output, integer ( kind = 4 ) KONTRL.  This routine receives the
!    current value of KONTRL, and may reset it.  The change is effective only
!    for the current error message.  This allows the user to suppress
!    or force printing of certain messages, for instance.
!
  implicit none

  integer ( kind = 4 ) kontrl
  integer ( kind = 4 ) level
  character ( len = * ) librar
  character ( len = * ) messg
  integer ( kind = 4 ) nerr
  character ( len = * ) subrou

  return
end
subroutine xerhlt ( messg )

!*****************************************************************************80
!
!! XERHLT aborts program execution.
!
!  Discussion:
!
!    This routine aborts the execution of the program.
!
!    The error message causing the abort is given in the calling
!    sequence.
!
!    This routine is used when the error message handler XERMSG is
!    employed.  The similar routine XERABT is to be used when the
!    older error message handlers XERROR and XERRWV are used.
!
!  Modified:
!
!    06 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the error message associated
!    with the halt in execution.
!
  implicit none

  character ( len = * ) messg

  stop
end
subroutine xerprn ( prefix, npref, messg, nwrap )

!*****************************************************************************80
!
!! XERPRN prints error messages processed by XERMSG.
!
!  Description:
!
!  Discussion:
!
!    This routine is used by the error handling routine XERMSG.  A related
!    routine, XERPRT, is used by the older error handling routines
!    XERROR and XERRWV.
!
!    This routine sends one or more lines to each of the (up to five)
!    logical units to which error messages are to be sent.  This routine
!    is called several times by XERMSG, sometimes with a single line to
!    print and sometimes with a (potentially very long) message that may
!    wrap around into multiple lines.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Kirby Fong
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) PREFIX, a string to be put at the beginning of
!    each line before the body of the message.  No more than 16 characters
!    of PREFIX will be used.
!
!    Input, integer ( kind = 4 ) NPREF, the number of characters to use from 
!    PREFIX.  If it is negative, the intrinsic function LEN is used to 
!    determine its length.  If it is zero, PREFIX is not used.  If it exceeds 
!    16 or if LEN(PREFIX) exceeds 16, only the first 16 characters will be
!    used.  If NPREF is positive and the length of PREFIX is less
!    than NPREF, a copy of PREFIX extended with blanks to length
!    NPREF will be used.
!
!    Input, character ( len = * ) MESSG, the error message.  If it is a long
!    message, it will be broken into pieces for printing on multiple lines.
!    Each line starts with the appropriate prefix and be followed by a piece
!    of the message.  NWRAP is the number of characters per piece; that is,
!    after each NWRAP characters, we break and start a new line.  In addition,
!    the characters '$$' embedded in MESSG are a sentinel for a new line.
!    The counting of characters up to NWRAP starts over for each new line.
!    The value of NWRAP typically used by XERMSG is 72 since many
!    older error messages in the SLATEC Library are laid out to rely on
!    wrap-around every 72 characters.
!
!    Input, integer ( kind = 4 ) NWRAP, the maximum size piece into which to 
!    break MESSG for printing on multiple lines.  An embedded '$$' ends a line,
!    and the count restarts at the following character.  If a line break does 
!    not occur on a blank (it would split a word) that word is moved to the 
!    next line.  Values of NWRAP less than 16 will be treated as 16.  Values of
!    NWRAP greater than 132 will be treated as 132.  The actual line length will
!    be NPREF + NWRAP after NPREF has been adjusted to fall between 0 and 16
!    and NWRAP has been adjusted to fall between 16 and 132.
!
  implicit none

  character ( len = 148 ) cbuff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) idelta
  integer ( kind = 4 ) iu(5)
  integer ( kind = 4 ) lenmsg
  integer ( kind = 4 ) lpiece
  integer ( kind = 4 ) lpref
  integer ( kind = 4 ) lwrap
  character ( len = * ) messg
  integer ( kind = 4 ) n
  character ( len = 2 ), parameter :: newlin = '$$'
  integer ( kind = 4 ) nextc
  integer ( kind = 4 ) npref
  integer ( kind = 4 ) nunit
  integer ( kind = 4 ) nwrap
  character ( len = * ) prefix

  call xgetua ( iu, nunit )
!
!  A zero value for a logical unit number means to use the standard
!  error message unit instead.  I1MACH(4) retrieves the standard
!  error message unit.
!
  n = i1mach(4)
  do i = 1, nunit
    if ( iu(i) == 0 ) then
      iu(i) = n
    end if
  end do
!
!  LPREF is the length of the prefix.  The prefix is placed at the
!  beginning of CBUFF, the character buffer, and kept there during
!  the rest of this routine.
!
  if ( npref < 0 ) then
    lpref = len ( prefix )
  else
    lpref = npref
  end if

  lpref = min ( 16, lpref )

  if ( lpref /= 0 ) then
    cbuff(1:lpref) = prefix
  end if
!
!  LWRAP is the maximum number of characters we want to take at one
!  time from MESSG to print on one line.
!
  lwrap = max ( 16, min ( 132, nwrap ) )
!
!  Set LENMSG to the length of MESSG, ignore any trailing blanks.
!
  lenmsg = len ( messg )
  n = lenmsg
  do i = 1, n
    if ( messg(lenmsg:lenmsg) /= ' ' ) then
      exit
    end if
    lenmsg = lenmsg - 1
  end do
!
!  If the message is all blanks, then print one blank line.
!
  if ( lenmsg == 0 ) then
    cbuff(lpref+1:lpref+1) = ' '
    do i = 1, nunit
      write ( iu(i), '(a)' ) cbuff(1:lpref+1)
    end do
    return
  end if
!
!  Set NEXTC to the position in MESSG where the next substring
!  starts.  From this position we scan for the new line sentinel.
!  When NEXTC exceeds LENMSG, there is no more to print.
!  We loop back to label 50 until all pieces have been printed.
!
!  We look for the next occurrence of the new line sentinel.  The
!  INDEX intrinsic function returns zero if there is no occurrence
!  or if the length of the first argument is less than the length
!  of the second argument.
!
!  There are several cases which should be checked for in the
!  following order.  We are attempting to set LPIECE to the number
!  of characters that should be taken from MESSG starting at
!  position NEXTC.
!
!  * LPIECE == 0
!  The new line sentinel does not occur in the remainder of the
!  character string.  LPIECE should be set to LWRAP or LENMSG+1-NEXTC,
!  whichever is less.
!
!  * LPIECE == 1
!  The new line sentinel starts at MESSG(NEXTC:NEXTC).  LPIECE is effectively
!  zero, and we print nothing to avoid producing unnecessary blank lines.
!  This takes care of the situation where the library routine has a message of
!  exactly 72 characters followed by a new line sentinel followed by more
!  characters.  NEXTC should be incremented by 2.
!
!  * LWRAP + 1 < LPIECE
!  Reduce LPIECE to LWRAP.
!
!  * Otherwise
!  This last case means 2 <= LPIECE <= LWRAP+1.  Reset LPIECE = LPIECE-1.
!  Note that this properly handles the end case where LPIECE = LWRAP+1.
!  That is, the sentinel falls exactly at the end of a line.
!
  nextc = 1

  do

    lpiece = index ( messg(nextc:lenmsg), newlin )
!
!  There was no new line sentinel found.
!
    if ( lpiece == 0 ) then

      idelta = 0
      lpiece = min ( lwrap, lenmsg + 1 - nextc )

      if ( lpiece < lenmsg + 1 - nextc ) then
        do i = lpiece+1, 2, -1
          if ( messg(nextc+i-1:nextc+i-1) == ' ' ) then
            lpiece = i - 1
            idelta = 1
            exit
          end if
        end do
      end if

      cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
      nextc = nextc + lpiece + idelta
!
!  We have a new line sentinel at MESSG(NEXTC:NEXTC+1).
!  Don't print a blank line.
!
    else if ( lpiece == 1 ) then

      nextc = nextc + 2
      cycle
!
!  LPIECE should be set down to LWRAP.
!
    else if ( lwrap + 1 < lpiece ) then

      idelta = 0
      lpiece = lwrap

      do i = lpiece + 1, 2, -1
        if ( messg(nextc+i-1:nextc+i-1) == ' ' ) then
          lpiece = i - 1
          idelta = 1
          exit
        end if
      end do

      cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
      nextc = nextc + lpiece + idelta
!
!  If we arrive here, it means 2 <= LPIECE <= LWRAP+1.
!  We should decrement LPIECE by one.
!
    else

      lpiece = lpiece - 1
      cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
      nextc = nextc + lpiece + 2

    end if
!
!  Print.
!
    do i = 1, nunit
      write ( iu(i), '(a)' ) cbuff(1:lpref+lpiece)
    end do

    if ( lenmsg < nextc ) then
      exit
    end if

  end do

  return
end
subroutine xersve ( librar, subrou, messg, kflag, nerr, level, icount )

!*****************************************************************************80
!
!! XERSVE records that an error has occurred.
!
!  Discussion:
!
!    This routine is used by the error handling routines associated
!    with XERMSG.  It is a revised version of the routine XERSAV, which
!    was used with the older pair of error handling routines XERROR
!    and XERRWV.
!
!  Modified:
!
!    06 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character (len = * ) LIBRAR, the name of the library or software
!    package from which the error message comes.
!
!    Input, character (len = * ) SUBROU, the name of the subroutine or function
!    from which the error message comes.
!
!    Input, character (len = * ) MESSG, the error message.
!
!    Input, integer ( kind = 4 ) KFLAG, indicates the action to be performed.
!    0 < KFLAG, the message in MESSG is saved.
!    KFLAG=0 the tables will be dumped and cleared.
!    KFLAG < 0, the tables will be dumped and not cleared.
!
!    Input, integer ( kind = 4 ) NERR, the error number.
!
!    Input, integer ( kind = 4 ) LEVEL, the error severity.
!
!    Output, integer ( kind = 4 ) ICOUNT, the number of times this message has 
!    been seen, or zero if the table has overflowed and does not contain this
!    message specifically.  When KFLAG=0, ICOUNT will not be altered from its
!    input value.
!
  implicit none

  integer ( kind = 4 ), parameter :: lentab = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) kflag
  integer ( kind = 4 ), save, dimension ( lentab ) :: kount
  integer ( kind = 4 ), save :: kountx = 0
  integer ( kind = 4 ) kunit
  integer ( kind = 4 ) level
  integer ( kind = 4 ), save, dimension ( lentab ) :: levtab
  character ( len = 8 ) lib
  character (len = * ) librar
  character ( len = 8 ), save, dimension ( lentab ) :: libtab
  integer ( kind = 4 ) lun(5)
  character ( len = 20 ) mes
  character (len = * ) messg
  character ( len = 20 ), save, dimension ( lentab ) :: mestab
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ), save, dimension ( lentab ) :: nertab
  integer ( kind = 4 ), save :: nmsg = 0
  integer ( kind = 4 ) nunit
  character ( len = 8 ) sub
  character (len = * ) subrou
  character ( len = 8 ), save, dimension ( lentab ) :: subtab

  if ( kflag <= 0 ) then
!
!  Dump the table.
!
    if ( nmsg == 0 ) then
      return
    end if
!
!  Print to each unit.
!
    call xgetua ( lun, nunit )

    do kunit = 1, nunit

      iunit = lun(kunit)
      if ( iunit == 0 ) then
        iunit = i1mach(4)
      end if
!
!  Print the table header.
!
      write ( iunit, '(a)' ) ' '
      write ( iunit, '(a)' ) '          Error message summary'
      write ( iunit, '(a,a)' ) &
        'Library    Subroutine Message start             NERR', &
        '     Level     Count'

!
!  Print body of table.
!
      do i = 1, nmsg
        write ( iunit, '(a,3x,a,3x,a,3i10)' ) &
          libtab(i), subtab(i), mestab(i), nertab(i), levtab(i), kount(i)
      end do
!
!  Print the number of other errors.
!
      if ( kountx /= 0 ) then
        write ( iunit, '(a)' ) ' '
        write ( iunit, '(a,i10)' ) &
          'Other errors not individually tabulated = ', kountx
      end if

      write ( iunit, '(1x)' )

    end do
!
!  Clear the error tables.
!
    if ( kflag == 0 ) then
      nmsg = 0
      kountx = 0
    end if

  else
!
!  Process a message.
!
!  Search for this message, or else an empty slot for this message,
!  or else determine that the error table is full.
!
    lib = librar
    sub = subrou
    mes = messg

    do i = 1, nmsg

      if ( &
        lib == libtab(i) .and.  &
        sub == subtab(i) .and.  &
        mes == mestab(i) .and.  &
        nerr == nertab(i) .and. &
        level == levtab(i) ) then
        kount(i) = kount(i) + 1
        icount = kount(i)
        return
      end if

    end do
!
!  Empty slot found for new message.
!
    if ( nmsg < lentab ) then

      nmsg = nmsg + 1
      libtab(i) = lib
      subtab(i) = sub
      mestab(i) = mes
      nertab(i) = nerr
      levtab(i) = level
      kount(i)  = 1
      icount    = 1
!
!  Table is full.
!
    else

      kountx = kountx + 1
      icount = 0

    end if

  end if

  return
end
subroutine xgetua ( iunit, nunit )

!*****************************************************************************80
!
!! XGETUA returns the unit numbers to which error messages are being sent.
!
!  Discussion:
!
!    These unit numbers may have been set by a call to XSETUN,
!    or a call to XSETUA, or may be default values.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT(5), an array into which the routine will
!    store the values of the NUNIT units to which the error messages
!    are being sent.  The value of NUNIT is never more than 5, so
!    using an array of dimension 5 will be sufficient.
!
!    Output, integer ( kind = 4 ) NUNIT, the number of units to which the
!    error messages are being sent.  NUNIT will be in the
!    range from 1 to 5.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit(5)
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) nunit
  logical set
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which

  which = 5
  value = 0
  set = .false.

  nunit = j4save ( which, value, set)

  if ( nunit < 1 ) then
    return
  end if

  which = 3
  value = 0
  set = .false.

  iunit(1) = j4save ( which, value, set )

  do i = 2, nunit

    which = i + 4
    value = 0
    set = .false.

    iunit(i) = j4save ( which, value, set )

  end do

  return
end
