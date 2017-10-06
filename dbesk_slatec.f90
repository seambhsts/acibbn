subroutine DBESK (X, FNU, KODE, N, Y, NZ)
!
!! DBESK implements forward recursion on the three term recursion ...
!            relation for a sequence of non-negative order Bessel
!            functions K/SUB(FNU+I-1)/(X), or scaled Bessel functions
!            EXP(X)*K/SUB(FNU+I-1)/(X), I=1,...,N for real, positive
!            X and non-negative orders FNU.
!
!***LIBRARY   SLATEC
!***CATEGORY  C10B3
!***TYPE      DOUBLE PRECISION (BESK-S, DBESK-D)
!***KEYWORDS  K BESSEL FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract  **** a double precision routine ****
!         DBESK implements forward recursion on the three term
!         recursion relation for a sequence of non-negative order Bessel
!         functions K/sub(FNU+I-1)/(X), or scaled Bessel functions
!         EXP(X)*K/sub(FNU+I-1)/(X), I=1,..,N for real X  >  0.0D0 and
!         non-negative orders FNU.  If FNU  <  NULIM, orders FNU and
!         FNU+1 are obtained from DBSKNU to start the recursion.  If
!         FNU  >=  NULIM, the uniform asymptotic expansion is used for
!         orders FNU and FNU+1 to start the recursion.  NULIM is 35 or
!         70 depending on whether N=1 or N  >=  2.  Under and overflow
!         tests are made on the leading term of the asymptotic expansion
!         before any extensive computation is done.
!
!         The maximum number of significant digits obtainable
!         is the smaller of 14 and the number of digits carried in
!         double precision arithmetic.
!
!     Description of Arguments
!
!         Input      X,FNU are double precision
!           X      - X  >  0.0D0
!           FNU    - order of the initial K function, FNU  >=  0.0D0
!           KODE   - a parameter to indicate the scaling option
!                    KODE=1 returns Y(I)=       K/sub(FNU+I-1)/(X),
!                                        I=1,...,N
!                    KODE=2 returns Y(I)=EXP(X)*K/sub(FNU+I-1)/(X),
!                                        I=1,...,N
!           N      - number of members in the sequence, N  >=  1
!
!         Output     Y is double precision
!           Y      - a vector whose first N components contain values
!                    for the sequence
!                    Y(I)=       k/sub(FNU+I-1)/(X), I=1,...,N  or
!                    Y(I)=EXP(X)*K/sub(FNU+I-1)/(X), I=1,...,N
!                    depending on KODE
!           NZ     - number of components of Y set to zero due to
!                    underflow with KODE=1,
!                    NZ=0   , normal return, computation completed
!                    NZ  /=  0, first NZ components of Y set to zero
!                             due to underflow, Y(I)=0.0D0, I=1,...,NZ
!
!     Error Conditions
!         Improper input arguments - a fatal error
!         Overflow - a fatal error
!         Underflow with KODE=1 -  a non-fatal error (NZ  /=  0)
!
!***REFERENCES  F. W. J. Olver, Tables of Bessel Functions of Moderate
!                 or Large Orders, NPL Mathematical Tables 6, Her
!                 Majesty's Stationery Office, London, 1962.
!               N. M. Temme, On the numerical evaluation of the modified
!                 Bessel function of the third kind, Journal of
!                 Computational Physics 19, (1975), pp. 324-337.
!***ROUTINES CALLED  D1MACH, DASYIK, DBESK0, DBESK1, DBSK0E, DBSK1E,
!                    DBSKNU, I1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790201  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBESK
!
  INTEGER I, J, K, KODE, MZ, N, NB, ND, NN, NUD, NULIM, NZ
  INTEGER I1MACH
  DOUBLE PRECISION CN,DNU,ELIM,ETX,FLGIK,FN,FNN,FNU,GLN,GNU,RTZ, &
   S, S1, S2, T, TM, TRX, W, X, XLIM, Y, ZN
  DOUBLE PRECISION DBESK0, DBESK1, DBSK1E, DBSK0E, D1MACH
  DIMENSION W(2), NULIM(2), Y(*)
  SAVE NULIM
  DATA NULIM(1),NULIM(2) / 35 , 70 /
!***FIRST EXECUTABLE STATEMENT  DBESK
  NN = -I1MACH(15)
  ELIM = 2.303D0*(NN*D1MACH(5)-3.0D0)
  XLIM = D1MACH(1)*1.0D+3
  if (KODE < 1 .OR. KODE > 2) go to 280
  if (FNU < 0.0D0) go to 290
  if (X <= 0.0D0) go to 300
  if (X < XLIM) go to 320
  if (N < 1) go to 310
  ETX = KODE - 1
!
!     ND IS A DUMMY VARIABLE FOR N
!     GNU IS A DUMMY VARIABLE FOR FNU
!     NZ = NUMBER OF UNDERFLOWS ON KODE=1
!
  ND = N
  NZ = 0
  NUD = INT(FNU)
  DNU = FNU - NUD
  GNU = FNU
  NN = MIN(2,ND)
  FN = FNU + N - 1
  FNN = FN
  if (FN < 2.0D0) go to 150
!
!     OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
!     FOR THE LAST ORDER, FNU+N-1 >= NULIM
!
  ZN = X/FN
  if (ZN == 0.0D0) go to 320
  RTZ = SQRT(1.0D0+ZN*ZN)
  GLN = LOG((1.0D0+RTZ)/ZN)
  T = RTZ*(1.0D0-ETX) + ETX/(ZN+RTZ)
  CN = -FN*(T-GLN)
  if (CN > ELIM) go to 320
  if (NUD < NULIM(NN)) go to 30
  if (NN == 1) go to 20
   10 CONTINUE
!
!     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
!     FOR THE FIRST ORDER, FNU >= NULIM
!
  FN = GNU
  ZN = X/FN
  RTZ = SQRT(1.0D0+ZN*ZN)
  GLN = LOG((1.0D0+RTZ)/ZN)
  T = RTZ*(1.0D0-ETX) + ETX/(ZN+RTZ)
  CN = -FN*(T-GLN)
   20 CONTINUE
  if (CN < -ELIM) go to 230
!
!     ASYMPTOTIC EXPANSION FOR ORDERS FNU AND FNU+1 >= NULIM
!
  FLGIK = -1.0D0
  call DASYIK(X,GNU,KODE,FLGIK,RTZ,CN,NN,Y)
  if (NN == 1) go to 240
  TRX = 2.0D0/X
  TM = (GNU+GNU+2.0D0)/X
  go to 130
!
   30 CONTINUE
  if (KODE == 2) go to 40
!
!     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION IN X)
!     FOR ORDER DNU
!
  if (X > ELIM) go to 230
   40 CONTINUE
  if (DNU /= 0.0D0) go to 80
  if (KODE == 2) go to 50
  S1 = DBESK0(X)
  go to 60
   50 S1 = DBSK0E(X)
   60 CONTINUE
  if (NUD == 0 .AND. ND == 1) go to 120
  if (KODE == 2) go to 70
  S2 = DBESK1(X)
  go to 90
   70 S2 = DBSK1E(X)
  go to 90
   80 CONTINUE
  NB = 2
  if (NUD == 0 .AND. ND == 1) NB = 1
  call DBSKNU(X, DNU, KODE, NB, W, NZ)
  S1 = W(1)
  if (NB == 1) go to 120
  S2 = W(2)
   90 CONTINUE
  TRX = 2.0D0/X
  TM = (DNU+DNU+2.0D0)/X
!     FORWARD RECUR FROM DNU TO FNU+1 TO GET Y(1) AND Y(2)
  if (ND == 1) NUD = NUD - 1
  if (NUD > 0) go to 100
  if (ND > 1) go to 120
  S1 = S2
  go to 120
  100 CONTINUE
  DO 110 I=1,NUD
    S = S2
    S2 = TM*S2 + S1
    S1 = S
    TM = TM + TRX
  110 CONTINUE
  if (ND == 1) S1 = S2
  120 CONTINUE
  Y(1) = S1
  if (ND == 1) go to 240
  Y(2) = S2
  130 CONTINUE
  if (ND == 2) go to 240
!     FORWARD RECUR FROM FNU+2 TO FNU+N-1
  DO 140 I=3,ND
    Y(I) = TM*Y(I-1) + Y(I-2)
    TM = TM + TRX
  140 CONTINUE
  go to 240
!
  150 CONTINUE
!     UNDERFLOW TEST FOR KODE=1
  if (KODE == 2) go to 160
  if (X > ELIM) go to 230
  160 CONTINUE
!     OVERFLOW TEST
  if (FN <= 1.0D0) go to 170
  if (-FN*(LOG(X)-0.693D0) > ELIM) go to 320
  170 CONTINUE
  if (DNU == 0.0D0) go to 180
  call DBSKNU(X, FNU, KODE, ND, Y, MZ)
  go to 240
  180 CONTINUE
  J = NUD
  if (J == 1) go to 210
  J = J + 1
  if (KODE == 2) go to 190
  Y(J) = DBESK0(X)
  go to 200
  190 Y(J) = DBSK0E(X)
  200 if (ND == 1) go to 240
  J = J + 1
  210 if (KODE == 2) go to 220
  Y(J) = DBESK1(X)
  go to 240
  220 Y(J) = DBSK1E(X)
  go to 240
!
!     UPDATE PARAMETERS ON UNDERFLOW
!
  230 CONTINUE
  NUD = NUD + 1
  ND = ND - 1
  if (ND == 0) go to 240
  NN = MIN(2,ND)
  GNU = GNU + 1.0D0
  if (FNN < 2.0D0) go to 230
  if (NUD < NULIM(NN)) go to 230
  go to 10
  240 CONTINUE
  NZ = N - ND
  if (NZ == 0) RETURN
  if (ND == 0) go to 260
  DO 250 I=1,ND
    J = N - I + 1
    K = ND - I + 1
    Y(J) = Y(K)
  250 CONTINUE
  260 CONTINUE
  DO 270 I=1,NZ
    Y(I) = 0.0D0
  270 CONTINUE
  return
!
!
!
  280 CONTINUE
  call XERMSG ('SLATEC', 'DBESK', &
     'SCALING OPTION, KODE, NOT 1 OR 2', 2, 1)
  return
  290 CONTINUE
  call XERMSG ('SLATEC', 'DBESK', 'ORDER, FNU, LESS THAN ZERO', 2, &
     1)
  return
  300 CONTINUE
  call XERMSG ('SLATEC', 'DBESK', 'X LESS THAN OR EQUAL TO ZERO', &
     2, 1)
  return
  310 CONTINUE
  call XERMSG ('SLATEC', 'DBESK', 'N LESS THAN ONE', 2, 1)
  return
  320 CONTINUE
  call XERMSG ('SLATEC', 'DBESK', &
     'OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL', 6, 1)
  return
end
DOUBLE PRECISION FUNCTION DBESK0 (X)
!
!! DBESK0 computes the modified (hyperbolic) Bessel function of the ...
!            third kind of order zero.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      DOUBLE PRECISION (BESK0-S, DBESK0-D)
!***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
!             THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBESK0(X) calculates the double precision modified (hyperbolic)
! Bessel function of the third kind of order zero for double
! precision argument X.  The argument must be greater than zero
! but not so large that the result underflows.
!
! Series for BK0        on the interval  0.          to  4.00000E+00
!                                        with weighted error   3.08E-33
!                                         log weighted error  32.51
!                               significant figures required  32.05
!                                    decimal places required  33.11
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DBESI0, DBSK0E, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBESK0
  DOUBLE PRECISION X, BK0CS(16), XMAX, XMAXT, XSML, Y, &
    D1MACH, DCSEVL, DBESI0, DBSK0E
  LOGICAL FIRST
  SAVE BK0CS, NTK0, XSML, XMAX, FIRST
  DATA BK0CS(  1) / -.353273932339027687201140060063153D-1    /
  DATA BK0CS(  2) / +.344289899924628486886344927529213D+0    /
  DATA BK0CS(  3) / +.359799365153615016265721303687231D-1    /
  DATA BK0CS(  4) / +.126461541144692592338479508673447D-2    /
  DATA BK0CS(  5) / +.228621210311945178608269830297585D-4    /
  DATA BK0CS(  6) / +.253479107902614945730790013428354D-6    /
  DATA BK0CS(  7) / +.190451637722020885897214059381366D-8    /
  DATA BK0CS(  8) / +.103496952576336245851008317853089D-10   /
  DATA BK0CS(  9) / +.425981614279108257652445327170133D-13   /
  DATA BK0CS( 10) / +.137446543588075089694238325440000D-15   /
  DATA BK0CS( 11) / +.357089652850837359099688597333333D-18   /
  DATA BK0CS( 12) / +.763164366011643737667498666666666D-21   /
  DATA BK0CS( 13) / +.136542498844078185908053333333333D-23   /
  DATA BK0CS( 14) / +.207527526690666808319999999999999D-26   /
  DATA BK0CS( 15) / +.271281421807298560000000000000000D-29   /
  DATA BK0CS( 16) / +.308259388791466666666666666666666D-32   /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBESK0
  if (FIRST) THEN
     NTK0 = INITDS (BK0CS, 16, 0.1*REAL(D1MACH(3)))
     XSML = SQRT(4.0D0*D1MACH(3))
     XMAXT = -LOG(D1MACH(1))
     XMAX = XMAXT - 0.5D0*XMAXT*LOG(XMAXT)/(XMAXT+0.5D0)
  end if
  FIRST = .FALSE.
!
  if (X  <=  0.D0) call XERMSG ('SLATEC', 'DBESK0', &
     'X IS ZERO OR NEGATIVE', 2, 2)
  if (X > 2.0D0) go to 20
!
  Y = 0.D0
  if (X > XSML) Y = X*X
  DBESK0 = -LOG(0.5D0*X)*DBESI0(X) - 0.25D0 + DCSEVL (.5D0*Y-1.D0, &
    BK0CS, NTK0)
  return
!
 20   DBESK0 = 0.D0
  if (X  >  XMAX) call XERMSG ('SLATEC', 'DBESK0', &
     'X SO BIG K0 UNDERFLOWS', 1, 1)
  if (X > XMAX) RETURN
!
  DBESK0 = EXP(-X) * DBSK0E(X)
!
  return
end
  DOUBLE PRECISION FUNCTION DBESK1 (X)
!
!! DBESK1 computes the modified (hyperbolic) Bessel function of the ...
!            third kind of order one.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      DOUBLE PRECISION (BESK1-S, DBESK1-D)
!***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
!             THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBESK1(X) calculates the double precision modified (hyperbolic)
! Bessel function of the third kind of order one for double precision
! argument X.  The argument must be large enough that the result does
! not overflow and small enough that the result does not underflow.
!
! Series for BK1        on the interval  0.          to  4.00000E+00
!                                        with weighted error   9.16E-32
!                                         log weighted error  31.04
!                               significant figures required  30.61
!                                    decimal places required  31.64
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DBESI1, DBSK1E, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBESK1
  DOUBLE PRECISION X, BK1CS(16), XMAX, XMAXT, XMIN, XSML, Y, &
    D1MACH, DCSEVL, DBESI1, DBSK1E
  LOGICAL FIRST
  SAVE BK1CS, NTK1, XMIN, XSML, XMAX, FIRST
  DATA BK1CS(  1) / +.25300227338947770532531120868533D-1     /
  DATA BK1CS(  2) / -.35315596077654487566723831691801D+0     /
  DATA BK1CS(  3) / -.12261118082265714823479067930042D+0     /
  DATA BK1CS(  4) / -.69757238596398643501812920296083D-2     /
  DATA BK1CS(  5) / -.17302889575130520630176507368979D-3     /
  DATA BK1CS(  6) / -.24334061415659682349600735030164D-5     /
  DATA BK1CS(  7) / -.22133876307347258558315252545126D-7     /
  DATA BK1CS(  8) / -.14114883926335277610958330212608D-9     /
  DATA BK1CS(  9) / -.66669016941993290060853751264373D-12    /
  DATA BK1CS( 10) / -.24274498505193659339263196864853D-14    /
  DATA BK1CS( 11) / -.70238634793862875971783797120000D-17    /
  DATA BK1CS( 12) / -.16543275155100994675491029333333D-19    /
  DATA BK1CS( 13) / -.32338347459944491991893333333333D-22    /
  DATA BK1CS( 14) / -.53312750529265274999466666666666D-25    /
  DATA BK1CS( 15) / -.75130407162157226666666666666666D-28    /
  DATA BK1CS( 16) / -.91550857176541866666666666666666D-31    /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBESK1
  if (FIRST) THEN
     NTK1 = INITDS (BK1CS, 16, 0.1*REAL(D1MACH(3)))
     XMIN = EXP(MAX(LOG(D1MACH(1)), -LOG(D1MACH(2))) + 0.01D0)
     XSML = SQRT(4.0D0*D1MACH(3))
     XMAXT = -LOG(D1MACH(1))
     XMAX = XMAXT - 0.5D0*XMAXT*LOG(XMAXT)/(XMAXT+0.5D0)
  end if
  FIRST = .FALSE.
!
  if (X  <=  0.D0) call XERMSG ('SLATEC', 'DBESK1', &
     'X IS ZERO OR NEGATIVE', 2, 2)
  if (X > 2.0D0) go to 20
!
  if (X  <  XMIN) call XERMSG ('SLATEC', 'DBESK1', &
     'X SO SMALL K1 OVERFLOWS', 3, 2)
  Y = 0.D0
  if (X > XSML) Y = X*X
  DBESK1 = LOG(0.5D0*X)*DBESI1(X) + (0.75D0 + DCSEVL (.5D0*Y-1.D0, &
    BK1CS, NTK1))/X
  return
!
 20   DBESK1 = 0.D0
  if (X  >  XMAX) call XERMSG ('SLATEC', 'DBESK1', &
     'X SO BIG K1 UNDERFLOWS', 1, 1)
  if (X > XMAX) RETURN
!
  DBESK1 = EXP(-X) * DBSK1E(X)
!
  return
end
subroutine DBESKS (XNU, X, NIN, BK)
!
!! DBESKS computes a sequence of modified Bessel functions of the ...
!            third kind of fractional order.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B3
!***TYPE      DOUBLE PRECISION (BESKS-S, DBESKS-D)
!***KEYWORDS  FNLIB, FRACTIONAL ORDER, MODIFIED BESSEL FUNCTION,
!             SEQUENCE OF BESSEL FUNCTIONS, SPECIAL FUNCTIONS,
!             THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBESKS computes a sequence of modified Bessel functions of the third
! kind of order XNU + I at X, where X  >  0, XNU lies in (-1,1),
! and I = 0, 1, ... , NIN - 1, if NIN is positive and I = 0, 1, ... ,
! NIN + 1, if NIN is negative.  On return, the vector BK(.) contains
! the results at X for order starting at XNU.  XNU, X, and BK are
! double precision.  NIN is an integer.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DBSKES, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBESKS
  DOUBLE PRECISION XNU, X, BK(*), EXPXI, XMAX, D1MACH
  SAVE XMAX
  DATA XMAX / 0.D0 /
!***FIRST EXECUTABLE STATEMENT  DBESKS
  if (XMAX == 0.D0) XMAX = -LOG (D1MACH(1))
!
  if (X  >  XMAX) call XERMSG ('SLATEC', 'DBESKS', &
     'X SO BIG BESSEL K UNDERFLOWS', 1, 2)
!
  call DBSKES (XNU, X, NIN, BK)
!
  EXPXI = EXP (-X)
  N = ABS (NIN)
  DO 20 I=1,N
    BK(I) = EXPXI * BK(I)
 20   CONTINUE
!
  return
end
  DOUBLE PRECISION FUNCTION DBSK0E (X)
!
!! DBSK0E computes the exponentially scaled modified (hyperbolic) Bessel ...
!  function of the third kind of order zero.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      DOUBLE PRECISION (BESK0E-S, DBSK0E-D)
!***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
!             THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBSK0E(X) computes the double precision exponentially scaled
! modified (hyperbolic) Bessel function of the third kind of
! order zero for positive double precision argument X.
!
! Series for BK0        on the interval  0.          to  4.00000E+00
!                                        with weighted error   3.08E-33
!                                         log weighted error  32.51
!                               significant figures required  32.05
!                                    decimal places required  33.11
!
! Series for AK0        on the interval  1.25000E-01 to  5.00000E-01
!                                        with weighted error   2.85E-32
!                                         log weighted error  31.54
!                               significant figures required  30.19
!                                    decimal places required  32.33
!
! Series for AK02       on the interval  0.          to  1.25000E-01
!                                        with weighted error   2.30E-32
!                                         log weighted error  31.64
!                               significant figures required  29.68
!                                    decimal places required  32.40
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DBESI0, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBSK0E
  DOUBLE PRECISION X, BK0CS(16), AK0CS(38), AK02CS(33), &
    XSML, Y, D1MACH, DCSEVL, DBESI0
  LOGICAL FIRST
  SAVE BK0CS, AK0CS, AK02CS, NTK0, NTAK0, NTAK02, XSML, FIRST
  DATA BK0CS(  1) / -.353273932339027687201140060063153D-1    /
  DATA BK0CS(  2) / +.344289899924628486886344927529213D+0    /
  DATA BK0CS(  3) / +.359799365153615016265721303687231D-1    /
  DATA BK0CS(  4) / +.126461541144692592338479508673447D-2    /
  DATA BK0CS(  5) / +.228621210311945178608269830297585D-4    /
  DATA BK0CS(  6) / +.253479107902614945730790013428354D-6    /
  DATA BK0CS(  7) / +.190451637722020885897214059381366D-8    /
  DATA BK0CS(  8) / +.103496952576336245851008317853089D-10   /
  DATA BK0CS(  9) / +.425981614279108257652445327170133D-13   /
  DATA BK0CS( 10) / +.137446543588075089694238325440000D-15   /
  DATA BK0CS( 11) / +.357089652850837359099688597333333D-18   /
  DATA BK0CS( 12) / +.763164366011643737667498666666666D-21   /
  DATA BK0CS( 13) / +.136542498844078185908053333333333D-23   /
  DATA BK0CS( 14) / +.207527526690666808319999999999999D-26   /
  DATA BK0CS( 15) / +.271281421807298560000000000000000D-29   /
  DATA BK0CS( 16) / +.308259388791466666666666666666666D-32   /
  DATA AK0CS(  1) / -.7643947903327941424082978270088D-1      /
  DATA AK0CS(  2) / -.2235652605699819052023095550791D-1      /
  DATA AK0CS(  3) / +.7734181154693858235300618174047D-3      /
  DATA AK0CS(  4) / -.4281006688886099464452146435416D-4      /
  DATA AK0CS(  5) / +.3081700173862974743650014826660D-5      /
  DATA AK0CS(  6) / -.2639367222009664974067448892723D-6      /
  DATA AK0CS(  7) / +.2563713036403469206294088265742D-7      /
  DATA AK0CS(  8) / -.2742705549900201263857211915244D-8      /
  DATA AK0CS(  9) / +.3169429658097499592080832873403D-9      /
  DATA AK0CS( 10) / -.3902353286962184141601065717962D-10     /
  DATA AK0CS( 11) / +.5068040698188575402050092127286D-11     /
  DATA AK0CS( 12) / -.6889574741007870679541713557984D-12     /
  DATA AK0CS( 13) / +.9744978497825917691388201336831D-13     /
  DATA AK0CS( 14) / -.1427332841884548505389855340122D-13     /
  DATA AK0CS( 15) / +.2156412571021463039558062976527D-14     /
  DATA AK0CS( 16) / -.3349654255149562772188782058530D-15     /
  DATA AK0CS( 17) / +.5335260216952911692145280392601D-16     /
  DATA AK0CS( 18) / -.8693669980890753807639622378837D-17     /
  DATA AK0CS( 19) / +.1446404347862212227887763442346D-17     /
  DATA AK0CS( 20) / -.2452889825500129682404678751573D-18     /
  DATA AK0CS( 21) / +.4233754526232171572821706342400D-19     /
  DATA AK0CS( 22) / -.7427946526454464195695341294933D-20     /
  DATA AK0CS( 23) / +.1323150529392666866277967462400D-20     /
  DATA AK0CS( 24) / -.2390587164739649451335981465599D-21     /
  DATA AK0CS( 25) / +.4376827585923226140165712554666D-22     /
  DATA AK0CS( 26) / -.8113700607345118059339011413333D-23     /
  DATA AK0CS( 27) / +.1521819913832172958310378154666D-23     /
  DATA AK0CS( 28) / -.2886041941483397770235958613333D-24     /
  DATA AK0CS( 29) / +.5530620667054717979992610133333D-25     /
  DATA AK0CS( 30) / -.1070377329249898728591633066666D-25     /
  DATA AK0CS( 31) / +.2091086893142384300296328533333D-26     /
  DATA AK0CS( 32) / -.4121713723646203827410261333333D-27     /
  DATA AK0CS( 33) / +.8193483971121307640135680000000D-28     /
  DATA AK0CS( 34) / -.1642000275459297726780757333333D-28     /
  DATA AK0CS( 35) / +.3316143281480227195890346666666D-29     /
  DATA AK0CS( 36) / -.6746863644145295941085866666666D-30     /
  DATA AK0CS( 37) / +.1382429146318424677635413333333D-30     /
  DATA AK0CS( 38) / -.2851874167359832570811733333333D-31     /
  DATA AK02CS(  1) / -.1201869826307592239839346212452D-1      /
  DATA AK02CS(  2) / -.9174852691025695310652561075713D-2      /
  DATA AK02CS(  3) / +.1444550931775005821048843878057D-3      /
  DATA AK02CS(  4) / -.4013614175435709728671021077879D-5      /
  DATA AK02CS(  5) / +.1567831810852310672590348990333D-6      /
  DATA AK02CS(  6) / -.7770110438521737710315799754460D-8      /
  DATA AK02CS(  7) / +.4611182576179717882533130529586D-9      /
  DATA AK02CS(  8) / -.3158592997860565770526665803309D-10     /
  DATA AK02CS(  9) / +.2435018039365041127835887814329D-11     /
  DATA AK02CS( 10) / -.2074331387398347897709853373506D-12     /
  DATA AK02CS( 11) / +.1925787280589917084742736504693D-13     /
  DATA AK02CS( 12) / -.1927554805838956103600347182218D-14     /
  DATA AK02CS( 13) / +.2062198029197818278285237869644D-15     /
  DATA AK02CS( 14) / -.2341685117579242402603640195071D-16     /
  DATA AK02CS( 15) / +.2805902810643042246815178828458D-17     /
  DATA AK02CS( 16) / -.3530507631161807945815482463573D-18     /
  DATA AK02CS( 17) / +.4645295422935108267424216337066D-19     /
  DATA AK02CS( 18) / -.6368625941344266473922053461333D-20     /
  DATA AK02CS( 19) / +.9069521310986515567622348800000D-21     /
  DATA AK02CS( 20) / -.1337974785423690739845005311999D-21     /
  DATA AK02CS( 21) / +.2039836021859952315522088960000D-22     /
  DATA AK02CS( 22) / -.3207027481367840500060869973333D-23     /
  DATA AK02CS( 23) / +.5189744413662309963626359466666D-24     /
  DATA AK02CS( 24) / -.8629501497540572192964607999999D-25     /
  DATA AK02CS( 25) / +.1472161183102559855208038400000D-25     /
  DATA AK02CS( 26) / -.2573069023867011283812351999999D-26     /
  DATA AK02CS( 27) / +.4601774086643516587376640000000D-27     /
  DATA AK02CS( 28) / -.8411555324201093737130666666666D-28     /
  DATA AK02CS( 29) / +.1569806306635368939301546666666D-28     /
  DATA AK02CS( 30) / -.2988226453005757788979199999999D-29     /
  DATA AK02CS( 31) / +.5796831375216836520618666666666D-30     /
  DATA AK02CS( 32) / -.1145035994347681332155733333333D-30     /
  DATA AK02CS( 33) / +.2301266594249682802005333333333D-31     /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBSK0E
  if (FIRST) THEN
     ETA = 0.1*REAL(D1MACH(3))
     NTK0 = INITDS (BK0CS, 16, ETA)
     NTAK0 = INITDS (AK0CS, 38, ETA)
     NTAK02 = INITDS (AK02CS, 33, ETA)
     XSML = SQRT(4.0D0*D1MACH(3))
  end if
  FIRST = .FALSE.
!
  if (X  <=  0.D0) call XERMSG ('SLATEC', 'DBSK0E', &
     'X IS ZERO OR NEGATIVE', 2, 2)
  if (X > 2.0D0) go to 20
!
  Y = 0.D0
  if (X > XSML) Y = X*X
  DBSK0E = EXP(X)*(-LOG(0.5D0*X)*DBESI0(X) - 0.25D0 + &
    DCSEVL (.5D0*Y-1.D0, BK0CS, NTK0))
  return
!
 20   if (X <= 8.D0) DBSK0E = (1.25D0 + DCSEVL ((16.D0/X-5.D0)/3.D0, &
    AK0CS, NTAK0))/SQRT(X)
  if (X > 8.D0) DBSK0E = (1.25D0 + &
    DCSEVL (16.D0/X-1.D0, AK02CS, NTAK02))/SQRT(X)
!
  return
end
  DOUBLE PRECISION FUNCTION DBSK1E (X)
!
!! DBSK1E computes the exponentially scaled modified (hyperbolic) Bessel ...
!  function of the third kind of order one.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      DOUBLE PRECISION (BESK1E-S, DBSK1E-D)
!***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
!             THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBSK1E(S) computes the double precision exponentially scaled
! modified (hyperbolic) Bessel function of the third kind of order
! one for positive double precision argument X.
!
! Series for BK1        on the interval  0.          to  4.00000E+00
!                                        with weighted error   9.16E-32
!                                         log weighted error  31.04
!                               significant figures required  30.61
!                                    decimal places required  31.64
!
! Series for AK1        on the interval  1.25000E-01 to  5.00000E-01
!                                        with weighted error   3.07E-32
!                                         log weighted error  31.51
!                               significant figures required  30.71
!                                    decimal places required  32.30
!
! Series for AK12       on the interval  0.          to  1.25000E-01
!                                        with weighted error   2.41E-32
!                                         log weighted error  31.62
!                               significant figures required  30.25
!                                    decimal places required  32.38
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DBESI1, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBSK1E
  DOUBLE PRECISION X, BK1CS(16), AK1CS(38), AK12CS(33), XMIN, &
    XSML, Y, D1MACH, DCSEVL, DBESI1
  LOGICAL FIRST
  SAVE BK1CS, AK1CS, AK12CS, NTK1, NTAK1, NTAK12, XMIN, XSML, &
    FIRST
  DATA BK1CS(  1) / +.25300227338947770532531120868533D-1     /
  DATA BK1CS(  2) / -.35315596077654487566723831691801D+0     /
  DATA BK1CS(  3) / -.12261118082265714823479067930042D+0     /
  DATA BK1CS(  4) / -.69757238596398643501812920296083D-2     /
  DATA BK1CS(  5) / -.17302889575130520630176507368979D-3     /
  DATA BK1CS(  6) / -.24334061415659682349600735030164D-5     /
  DATA BK1CS(  7) / -.22133876307347258558315252545126D-7     /
  DATA BK1CS(  8) / -.14114883926335277610958330212608D-9     /
  DATA BK1CS(  9) / -.66669016941993290060853751264373D-12    /
  DATA BK1CS( 10) / -.24274498505193659339263196864853D-14    /
  DATA BK1CS( 11) / -.70238634793862875971783797120000D-17    /
  DATA BK1CS( 12) / -.16543275155100994675491029333333D-19    /
  DATA BK1CS( 13) / -.32338347459944491991893333333333D-22    /
  DATA BK1CS( 14) / -.53312750529265274999466666666666D-25    /
  DATA BK1CS( 15) / -.75130407162157226666666666666666D-28    /
  DATA BK1CS( 16) / -.91550857176541866666666666666666D-31    /
  DATA AK1CS(  1) / +.27443134069738829695257666227266D+0     /
  DATA AK1CS(  2) / +.75719899531993678170892378149290D-1     /
  DATA AK1CS(  3) / -.14410515564754061229853116175625D-2     /
  DATA AK1CS(  4) / +.66501169551257479394251385477036D-4     /
  DATA AK1CS(  5) / -.43699847095201407660580845089167D-5     /
  DATA AK1CS(  6) / +.35402774997630526799417139008534D-6     /
  DATA AK1CS(  7) / -.33111637792932920208982688245704D-7     /
  DATA AK1CS(  8) / +.34459775819010534532311499770992D-8     /
  DATA AK1CS(  9) / -.38989323474754271048981937492758D-9     /
  DATA AK1CS( 10) / +.47208197504658356400947449339005D-10    /
  DATA AK1CS( 11) / -.60478356628753562345373591562890D-11    /
  DATA AK1CS( 12) / +.81284948748658747888193837985663D-12    /
  DATA AK1CS( 13) / -.11386945747147891428923915951042D-12    /
  DATA AK1CS( 14) / +.16540358408462282325972948205090D-13    /
  DATA AK1CS( 15) / -.24809025677068848221516010440533D-14    /
  DATA AK1CS( 16) / +.38292378907024096948429227299157D-15    /
  DATA AK1CS( 17) / -.60647341040012418187768210377386D-16    /
  DATA AK1CS( 18) / +.98324256232648616038194004650666D-17    /
  DATA AK1CS( 19) / -.16284168738284380035666620115626D-17    /
  DATA AK1CS( 20) / +.27501536496752623718284120337066D-18    /
  DATA AK1CS( 21) / -.47289666463953250924281069568000D-19    /
  DATA AK1CS( 22) / +.82681500028109932722392050346666D-20    /
  DATA AK1CS( 23) / -.14681405136624956337193964885333D-20    /
  DATA AK1CS( 24) / +.26447639269208245978085894826666D-21    /
  DATA AK1CS( 25) / -.48290157564856387897969868800000D-22    /
  DATA AK1CS( 26) / +.89293020743610130180656332799999D-23    /
  DATA AK1CS( 27) / -.16708397168972517176997751466666D-23    /
  DATA AK1CS( 28) / +.31616456034040694931368618666666D-24    /
  DATA AK1CS( 29) / -.60462055312274989106506410666666D-25    /
  DATA AK1CS( 30) / +.11678798942042732700718421333333D-25    /
  DATA AK1CS( 31) / -.22773741582653996232867840000000D-26    /
  DATA AK1CS( 32) / +.44811097300773675795305813333333D-27    /
  DATA AK1CS( 33) / -.88932884769020194062336000000000D-28    /
  DATA AK1CS( 34) / +.17794680018850275131392000000000D-28    /
  DATA AK1CS( 35) / -.35884555967329095821994666666666D-29    /
  DATA AK1CS( 36) / +.72906290492694257991679999999999D-30    /
  DATA AK1CS( 37) / -.14918449845546227073024000000000D-30    /
  DATA AK1CS( 38) / +.30736573872934276300799999999999D-31    /
  DATA AK12CS(  1) / +.6379308343739001036600488534102D-1      /
  DATA AK12CS(  2) / +.2832887813049720935835030284708D-1      /
  DATA AK12CS(  3) / -.2475370673905250345414545566732D-3      /
  DATA AK12CS(  4) / +.5771972451607248820470976625763D-5      /
  DATA AK12CS(  5) / -.2068939219536548302745533196552D-6      /
  DATA AK12CS(  6) / +.9739983441381804180309213097887D-8      /
  DATA AK12CS(  7) / -.5585336140380624984688895511129D-9      /
  DATA AK12CS(  8) / +.3732996634046185240221212854731D-10     /
  DATA AK12CS(  9) / -.2825051961023225445135065754928D-11     /
  DATA AK12CS( 10) / +.2372019002484144173643496955486D-12     /
  DATA AK12CS( 11) / -.2176677387991753979268301667938D-13     /
  DATA AK12CS( 12) / +.2157914161616032453939562689706D-14     /
  DATA AK12CS( 13) / -.2290196930718269275991551338154D-15     /
  DATA AK12CS( 14) / +.2582885729823274961919939565226D-16     /
  DATA AK12CS( 15) / -.3076752641268463187621098173440D-17     /
  DATA AK12CS( 16) / +.3851487721280491597094896844799D-18     /
  DATA AK12CS( 17) / -.5044794897641528977117282508800D-19     /
  DATA AK12CS( 18) / +.6888673850418544237018292223999D-20     /
  DATA AK12CS( 19) / -.9775041541950118303002132480000D-21     /
  DATA AK12CS( 20) / +.1437416218523836461001659733333D-21     /
  DATA AK12CS( 21) / -.2185059497344347373499733333333D-22     /
  DATA AK12CS( 22) / +.3426245621809220631645388800000D-23     /
  DATA AK12CS( 23) / -.5531064394246408232501248000000D-24     /
  DATA AK12CS( 24) / +.9176601505685995403782826666666D-25     /
  DATA AK12CS( 25) / -.1562287203618024911448746666666D-25     /
  DATA AK12CS( 26) / +.2725419375484333132349439999999D-26     /
  DATA AK12CS( 27) / -.4865674910074827992378026666666D-27     /
  DATA AK12CS( 28) / +.8879388552723502587357866666666D-28     /
  DATA AK12CS( 29) / -.1654585918039257548936533333333D-28     /
  DATA AK12CS( 30) / +.3145111321357848674303999999999D-29     /
  DATA AK12CS( 31) / -.6092998312193127612416000000000D-30     /
  DATA AK12CS( 32) / +.1202021939369815834623999999999D-30     /
  DATA AK12CS( 33) / -.2412930801459408841386666666666D-31     /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBSK1E
  if (FIRST) THEN
     ETA = 0.1*REAL(D1MACH(3))
     NTK1 = INITDS (BK1CS, 16, ETA)
     NTAK1 = INITDS (AK1CS, 38, ETA)
     NTAK12 = INITDS (AK12CS, 33, ETA)
!
     XMIN = EXP (MAX(LOG(D1MACH(1)), -LOG(D1MACH(2))) + 0.01D0)
     XSML = SQRT(4.0D0*D1MACH(3))
  end if
  FIRST = .FALSE.
!
  if (X  <=  0.D0) call XERMSG ('SLATEC', 'DBSK1E', &
     'X IS ZERO OR NEGATIVE', 2, 2)
  if (X > 2.0D0) go to 20
!
  if (X  <  XMIN) call XERMSG ('SLATEC', 'DBSK1E', &
     'X SO SMALL K1 OVERFLOWS', 3, 2)
  Y = 0.D0
  if (X > XSML) Y = X*X
  DBSK1E = EXP(X)*(LOG(0.5D0*X)*DBESI1(X) + (0.75D0 + &
    DCSEVL (0.5D0*Y-1.D0, BK1CS, NTK1))/X )
  return
!
 20   if (X <= 8.D0) DBSK1E = (1.25D0 + DCSEVL ((16.D0/X-5.D0)/3.D0, &
    AK1CS, NTAK1))/SQRT(X)
  if (X > 8.D0) DBSK1E = (1.25D0 + &
    DCSEVL (16.D0/X-1.D0, AK12CS, NTAK12))/SQRT(X)
!
  return
end
subroutine DBSKNU (X, FNU, KODE, N, Y, NZ)
!
!! DBSKNU is subsidiary to DBESK.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (BESKNU-S, DBSKNU-D)
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract  **** A DOUBLE PRECISION routine ****
!         DBSKNU computes N member sequences of K Bessel functions
!         K/SUB(FNU+I-1)/(X), I=1,N for non-negative orders FNU and
!         positive X. Equations of the references are implemented on
!         small orders DNU for K/SUB(DNU)/(X) and K/SUB(DNU+1)/(X).
!         Forward recursion with the three term recursion relation
!         generates higher orders FNU+I-1, I=1,...,N. The parameter
!         KODE permits K/SUB(FNU+I-1)/(X) values or scaled values
!         EXP(X)*K/SUB(FNU+I-1)/(X), I=1,N to be returned.
!
!         To start the recursion FNU is normalized to the interval
!         -0.5 <= DNU < 0.5. A special form of the power series is
!         implemented on 0 < X <= X1 while the Miller algorithm for the
!         K Bessel function in terms of the confluent hypergeometric
!         function U(FNU+0.5,2*FNU+1,X) is implemented on X1 < X <= X2.
!         For X > X2, the asymptotic expansion for large X is used.
!         When FNU is a half odd integer, a special formula for
!         DNU=-0.5 and DNU+1.0=0.5 is used to start the recursion.
!
!         The maximum number of significant digits obtainable
!         is the smaller of 14 and the number of digits carried in
!         DOUBLE PRECISION arithmetic.
!
!         DBSKNU assumes that a significant digit SINH function is
!         available.
!
!     Description of Arguments
!
!         INPUT      X,FNU are DOUBLE PRECISION
!           X      - X > 0.0D0
!           FNU    - Order of initial K function, FNU >= 0.0D0
!           N      - Number of members of the sequence, N >= 1
!           KODE   - A parameter to indicate the scaling option
!                    KODE= 1  returns
!                             Y(I)=       K/SUB(FNU+I-1)/(X)
!                                  I=1,...,N
!                        = 2  returns
!                             Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X)
!                                  I=1,...,N
!
!         OUTPUT     Y is DOUBLE PRECISION
!           Y      - A vector whose first N components contain values
!                    for the sequence
!                    Y(I)=       K/SUB(FNU+I-1)/(X), I=1,...,N or
!                    Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X), I=1,...,N
!                    depending on KODE
!           NZ     - Number of components set to zero due to
!                    underflow,
!                    NZ= 0   , normal return
!                    NZ /= 0 , first NZ components of Y set to zero
!                              due to underflow, Y(I)=0.0D0,I=1,...,NZ
!
!     Error Conditions
!         Improper input arguments - a fatal error
!         Overflow - a fatal error
!         Underflow with KODE=1 - a non-fatal error (NZ /= 0)
!
!***SEE ALSO  DBESK
!***REFERENCES  N. M. Temme, On the numerical evaluation of the modified
!                 Bessel function of the third kind, Journal of
!                 Computational Physics 19, (1975), pp. 324-337.
!***ROUTINES CALLED  D1MACH, DGAMMA, I1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790201  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!   900727  Added EXTERNAL statement.  (WRB)
!   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBSKNU
!
  INTEGER I, IFLAG, INU, J, K, KK, KODE, KODED, N, NN, NZ
  INTEGER I1MACH
  DOUBLE PRECISION A,AK,A1,A2,B,BK,CC,CK,COEF,CX,DK,DNU,DNU2,ELIM, &
   ETEST, EX, F, FC, FHS, FK, FKS, FLRX, FMU, FNU, G1, G2, P, PI, &
   PT, P1, P2, Q, RTHPI, RX, S, SMU, SQK, ST, S1, S2, TM, TOL, T1, &
   T2, X, X1, X2, Y
  DIMENSION A(160), B(160), Y(*), CC(8)
  DOUBLE PRECISION DGAMMA, D1MACH
  EXTERNAL DGAMMA
  SAVE X1, X2, PI, RTHPI, CC
  DATA X1, X2 / 2.0D0, 17.0D0 /
  DATA PI,RTHPI        / 3.14159265358979D+00, 1.25331413731550D+00/
  DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8) &
                       / 5.77215664901533D-01,-4.20026350340952D-02, &
  -4.21977345555443D-02, 7.21894324666300D-03,-2.15241674114900D-04, &
  -2.01348547807000D-05, 1.13302723200000D-06, 6.11609500000000D-09/
!***FIRST EXECUTABLE STATEMENT  DBSKNU
  KK = -I1MACH(15)
  ELIM = 2.303D0*(KK*D1MACH(5)-3.0D0)
  AK = D1MACH(3)
  TOL = MAX(AK,1.0D-15)
  if (X <= 0.0D0) go to 350
  if (FNU < 0.0D0) go to 360
  if (KODE < 1 .OR. KODE > 2) go to 370
  if (N < 1) go to 380
  NZ = 0
  IFLAG = 0
  KODED = KODE
  RX = 2.0D0/X
  INU = INT(FNU+0.5D0)
  DNU = FNU - INU
  if (ABS(DNU) == 0.5D0) go to 120
  DNU2 = 0.0D0
  if (ABS(DNU) < TOL) go to 10
  DNU2 = DNU*DNU
   10 CONTINUE
  if (X > X1) go to 120
!
!     SERIES FOR X <= X1
!
  A1 = 1.0D0 - DNU
  A2 = 1.0D0 + DNU
  T1 = 1.0D0/DGAMMA(A1)
  T2 = 1.0D0/DGAMMA(A2)
  if (ABS(DNU) > 0.1D0) go to 40
!     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
  S = CC(1)
  AK = 1.0D0
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
  G2 = (T1+T2)*0.5D0
  SMU = 1.0D0
  FC = 1.0D0
  FLRX = LOG(RX)
  FMU = DNU*FLRX
  if (DNU == 0.0D0) go to 60
  FC = DNU*PI
  FC = FC/SIN(FC)
  if (FMU /= 0.0D0) SMU = SINH(FMU)/FMU
   60 CONTINUE
  F = FC*(G1*COSH(FMU)+G2*FLRX*SMU)
  FC = EXP(FMU)
  P = 0.5D0*FC/T2
  Q = 0.5D0/(FC*T1)
  AK = 1.0D0
  CK = 1.0D0
  BK = 1.0D0
  S1 = F
  S2 = P
  if (INU > 0 .OR. N > 1) go to 90
  if (X < TOL) go to 80
  CX = X*X*0.25D0
   70 CONTINUE
  F = (AK*F+P+Q)/(BK-DNU2)
  P = P/(AK-DNU)
  Q = Q/(AK+DNU)
  CK = CK*CX/AK
  T1 = CK*F
  S1 = S1 + T1
  BK = BK + AK + AK + 1.0D0
  AK = AK + 1.0D0
  S = ABS(T1)/(1.0D0+ABS(S1))
  if (S > TOL) go to 70
   80 CONTINUE
  Y(1) = S1
  if (KODED == 1) RETURN
  Y(1) = S1*EXP(X)
  return
   90 CONTINUE
  if (X < TOL) go to 110
  CX = X*X*0.25D0
  100 CONTINUE
  F = (AK*F+P+Q)/(BK-DNU2)
  P = P/(AK-DNU)
  Q = Q/(AK+DNU)
  CK = CK*CX/AK
  T1 = CK*F
  S1 = S1 + T1
  T2 = CK*(P-AK*F)
  S2 = S2 + T2
  BK = BK + AK + AK + 1.0D0
  AK = AK + 1.0D0
  S = ABS(T1)/(1.0D0+ABS(S1)) + ABS(T2)/(1.0D0+ABS(S2))
  if (S > TOL) go to 100
  110 CONTINUE
  S2 = S2*RX
  if (KODED == 1) go to 170
  F = EXP(X)
  S1 = S1*F
  S2 = S2*F
  go to 170
  120 CONTINUE
  COEF = RTHPI/SQRT(X)
  if (KODED == 2) go to 130
  if (X > ELIM) go to 330
  COEF = COEF*EXP(-X)
  130 CONTINUE
  if (ABS(DNU) == 0.5D0) go to 340
  if (X > X2) go to 280
!
!     MILLER ALGORITHM FOR X1 < X <= X2
!
  ETEST = COS(PI*DNU)/(PI*X*TOL)
  FKS = 1.0D0
  FHS = 0.25D0
  FK = 0.0D0
  CK = X + X + 2.0D0
  P1 = 0.0D0
  P2 = 1.0D0
  K = 0
  140 CONTINUE
  K = K + 1
  FK = FK + 1.0D0
  AK = (FHS-DNU2)/(FKS+FK)
  BK = CK/(FK+1.0D0)
  PT = P2
  P2 = BK*P2 - AK*P1
  P1 = PT
  A(K) = AK
  B(K) = BK
  CK = CK + 2.0D0
  FKS = FKS + FK + FK + 1.0D0
  FHS = FHS + FK + FK
  if (ETEST > FK*P1) go to 140
  KK = K
  S = 1.0D0
  P1 = 0.0D0
  P2 = 1.0D0
  DO 150 I=1,K
    PT = P2
    P2 = (B(KK)*P2-P1)/A(KK)
    P1 = PT
    S = S + P2
    KK = KK - 1
  150 CONTINUE
  S1 = COEF*(P2/S)
  if (INU > 0 .OR. N > 1) go to 160
  go to 200
  160 CONTINUE
  S2 = S1*(X+DNU+0.5D0-P1/P2)/X
!
!     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION
!
  170 CONTINUE
  CK = (DNU+DNU+2.0D0)/X
  if (N == 1) INU = INU - 1
  if (INU > 0) go to 180
  if (N > 1) go to 200
  S1 = S2
  go to 200
  180 CONTINUE
  DO 190 I=1,INU
    ST = S2
    S2 = CK*S2 + S1
    S1 = ST
    CK = CK + RX
  190 CONTINUE
  if (N == 1) S1 = S2
  200 CONTINUE
  if (IFLAG == 1) go to 220
  Y(1) = S1
  if (N == 1) RETURN
  Y(2) = S2
  if (N == 2) RETURN
  DO 210 I=3,N
    Y(I) = CK*Y(I-1) + Y(I-2)
    CK = CK + RX
  210 CONTINUE
  return
!     IFLAG=1 CASES
  220 CONTINUE
  S = -X + LOG(S1)
  Y(1) = 0.0D0
  NZ = 1
  if (S < -ELIM) go to 230
  Y(1) = EXP(S)
  NZ = 0
  230 CONTINUE
  if (N == 1) RETURN
  S = -X + LOG(S2)
  Y(2) = 0.0D0
  NZ = NZ + 1
  if (S < -ELIM) go to 240
  NZ = NZ - 1
  Y(2) = EXP(S)
  240 CONTINUE
  if (N == 2) RETURN
  KK = 2
  if (NZ < 2) go to 260
  DO 250 I=3,N
    KK = I
    ST = S2
    S2 = CK*S2 + S1
    S1 = ST
    CK = CK + RX
    S = -X + LOG(S2)
    NZ = NZ + 1
    Y(I) = 0.0D0
    if (S < -ELIM) go to 250
    Y(I) = EXP(S)
    NZ = NZ - 1
    go to 260
  250 CONTINUE
  return
  260 CONTINUE
  if (KK == N) RETURN
  S2 = S2*CK + S1
  CK = CK + RX
  KK = KK + 1
  Y(KK) = EXP(-X+LOG(S2))
  if (KK == N) RETURN
  KK = KK + 1
  DO 270 I=KK,N
    Y(I) = CK*Y(I-1) + Y(I-2)
    CK = CK + RX
  270 CONTINUE
  return
!
!     ASYMPTOTIC EXPANSION FOR LARGE X, X > X2
!
!     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
!     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
!     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
!     RECURSION
  280 CONTINUE
  NN = 2
  if (INU == 0 .AND. N == 1) NN = 1
  DNU2 = DNU + DNU
  FMU = 0.0D0
  if (ABS(DNU2) < TOL) go to 290
  FMU = DNU2*DNU2
  290 CONTINUE
  EX = X*8.0D0
  S2 = 0.0D0
  DO 320 K=1,NN
    S1 = S2
    S = 1.0D0
    AK = 0.0D0
    CK = 1.0D0
    SQK = 1.0D0
    DK = EX
    DO 300 J=1,30
      CK = CK*(FMU-SQK)/DK
      S = S + CK
      DK = DK + EX
      AK = AK + 8.0D0
      SQK = SQK + AK
      if (ABS(CK) < TOL) go to 310
  300   CONTINUE
  310   S2 = S*COEF
    FMU = FMU + 8.0D0*DNU + 4.0D0
  320 CONTINUE
  if (NN > 1) go to 170
  S1 = S2
  go to 200
  330 CONTINUE
  KODED = 2
  IFLAG = 1
  go to 120
!
!     FNU=HALF ODD INTEGER CASE
!
  340 CONTINUE
  S1 = COEF
  S2 = COEF
  go to 170
!
!
  350 call XERMSG ('SLATEC', 'DBSKNU', 'X NOT GREATER THAN ZERO', 2, 1)
  return
  360 call XERMSG ('SLATEC', 'DBSKNU', 'FNU NOT ZERO OR POSITIVE', 2, &
     1)
  return
  370 call XERMSG ('SLATEC', 'DBSKNU', 'KODE NOT 1 OR 2', 2, 1)
  return
  380 call XERMSG ('SLATEC', 'DBSKNU', 'N NOT GREATER THAN 0', 2, 1)
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
!    NERR     An integer value that is chosen by the library routine''s
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
!                the user''s routine.  The user may also permit the error
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
subroutine DASYIK (X, FNU, KODE, FLGIK, RA, ARG, IN, Y)
!
!! DASYIK is subsidiary to DBESI and DBESK.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (ASYIK-S, DASYIK-D)
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!                    DASYIK computes Bessel functions I and K
!                  for arguments X > 0.0 and orders FNU >= 35
!                  on FLGIK = 1 and FLGIK = -1 respectively.
!
!                                    INPUT
!
!      X    - Argument, X > 0.0D0
!      FNU  - Order of first Bessel function
!      KODE - A parameter to indicate the scaling option
!             KODE=1 returns Y(I)=        I/SUB(FNU+I-1)/(X), I=1,IN
!                    or      Y(I)=        K/SUB(FNU+I-1)/(X), I=1,IN
!                    on FLGIK = 1.0D0 or FLGIK = -1.0D0
!             KODE=2 returns Y(I)=EXP(-X)*I/SUB(FNU+I-1)/(X), I=1,IN
!                    or      Y(I)=EXP( X)*K/SUB(FNU+I-1)/(X), I=1,IN
!                    on FLGIK = 1.0D0 or FLGIK = -1.0D0
!     FLGIK - Selection parameter for I or K FUNCTION
!             FLGIK =  1.0D0 gives the I function
!             FLGIK = -1.0D0 gives the K function
!        RA - SQRT(1.+Z*Z), Z=X/FNU
!       ARG - Argument of the leading exponential
!        IN - Number of functions desired, IN=1 or 2
!
!                                    OUTPUT
!
!         Y - A vector whose first IN components contain the sequence
!
!     Abstract  **** A double precision routine ****
!         DASYIK implements the uniform asymptotic expansion of
!         the I and K Bessel functions for FNU >= 35 and real
!         X > 0.0D0. The forms are identical except for a change
!         in sign of some of the terms. This change in sign is
!         accomplished by means of the FLAG FLGIK = 1 or -1.
!
!***SEE ALSO  DBESI, DBESK
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated the AUTHOR section.  (WRB)
!***END PROLOGUE  DASYIK
!
  INTEGER IN, J, JN, K, KK, KODE, L
  DOUBLE PRECISION AK,AP,ARG,C,COEF,CON,ETX,FLGIK,FN,FNU,GLN,RA, &
   S1, S2, T, TOL, T2, X, Y, Z
  DOUBLE PRECISION D1MACH
  DIMENSION Y(*), C(65), CON(2)
  SAVE CON, C
  DATA CON(1), CON(2)  / &
          3.98942280401432678D-01,    1.25331413731550025D+00/
  DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10), &
       C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18), &
       C(19), C(20), C(21), C(22), C(23), C(24)/ &
         -2.08333333333333D-01,        1.25000000000000D-01, &
          3.34201388888889D-01,       -4.01041666666667D-01, &
          7.03125000000000D-02,       -1.02581259645062D+00, &
          1.84646267361111D+00,       -8.91210937500000D-01, &
          7.32421875000000D-02,        4.66958442342625D+00, &
         -1.12070026162230D+01,        8.78912353515625D+00, &
         -2.36408691406250D+00,        1.12152099609375D-01, &
         -2.82120725582002D+01,        8.46362176746007D+01, &
         -9.18182415432400D+01,        4.25349987453885D+01, &
         -7.36879435947963D+00,        2.27108001708984D-01, &
          2.12570130039217D+02,       -7.65252468141182D+02, &
          1.05999045252800D+03,       -6.99579627376133D+02/
  DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32), &
       C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40), &
       C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/ &
          2.18190511744212D+02,       -2.64914304869516D+01, &
          5.72501420974731D-01,       -1.91945766231841D+03, &
          8.06172218173731D+03,       -1.35865500064341D+04, &
          1.16553933368645D+04,       -5.30564697861340D+03, &
          1.20090291321635D+03,       -1.08090919788395D+02, &
          1.72772750258446D+00,        2.02042913309661D+04, &
         -9.69805983886375D+04,        1.92547001232532D+05, &
         -2.03400177280416D+05,        1.22200464983017D+05, &
         -4.11926549688976D+04,        7.10951430248936D+03, &
         -4.93915304773088D+02,        6.07404200127348D+00, &
         -2.42919187900551D+05,        1.31176361466298D+06, &
         -2.99801591853811D+06,        3.76327129765640D+06/
  DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56), &
       C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64), &
       C(65)/ &
         -2.81356322658653D+06,        1.26836527332162D+06, &
         -3.31645172484564D+05,        4.52187689813627D+04, &
         -2.49983048181121D+03,        2.43805296995561D+01, &
          3.28446985307204D+06,       -1.97068191184322D+07, &
          5.09526024926646D+07,       -7.41051482115327D+07, &
          6.63445122747290D+07,       -3.75671766607634D+07, &
          1.32887671664218D+07,       -2.78561812808645D+06, &
          3.08186404612662D+05,       -1.38860897537170D+04, &
          1.10017140269247D+02/
!***FIRST EXECUTABLE STATEMENT  DASYIK
  TOL = D1MACH(3)
  TOL = MAX(TOL,1.0D-15)
  FN = FNU
  Z  = (3.0D0-FLGIK)/2.0D0
  KK = INT(Z)
  DO 50 JN=1,IN
    if (JN == 1) go to 10
    FN = FN - FLGIK
    Z = X/FN
    RA = SQRT(1.0D0+Z*Z)
    GLN = LOG((1.0D0+RA)/Z)
    ETX = KODE - 1
    T = RA*(1.0D0-ETX) + ETX/(Z+RA)
    ARG = FN*(T-GLN)*FLGIK
   10   COEF = EXP(ARG)
    T = 1.0D0/RA
    T2 = T*T
    T = T/FN
    T = SIGN(T,FLGIK)
    S2 = 1.0D0
    AP = 1.0D0
    L = 0
    DO 30 K=2,11
      L = L + 1
      S1 = C(L)
      DO 20 J=2,K
        L = L + 1
        S1 = S1*T2 + C(L)
   20     CONTINUE
      AP = AP*T
      AK = AP*S1
      S2 = S2 + AK
      if (MAX(ABS(AK),ABS(AP))  < TOL) go to 40
   30   CONTINUE
   40   CONTINUE
  T = ABS(T)
  Y(JN) = S2*COEF*SQRT(T)*CON(KK)
   50 CONTINUE
  return
end
DOUBLE PRECISION FUNCTION DBESI0 (X)
!
!! DBESI0 computes the hyperbolic Bessel function of first kind of order zero.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      DOUBLE PRECISION (BESI0-S, DBESI0-D)
!***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBESI0(X) calculates the double precision modified (hyperbolic)
! Bessel function of the first kind of order zero and double
! precision argument X.
!
! Series for BI0        on the interval  0.          to  9.00000E+00
!                                        with weighted error   9.51E-34
!                                         log weighted error  33.02
!                               significant figures required  33.31
!                                    decimal places required  33.65
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DBSI0E, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBESI0
  DOUBLE PRECISION X, BI0CS(18), XMAX, XSML, Y, D1MACH, &
    DCSEVL, DBSI0E
  LOGICAL FIRST
  SAVE BI0CS, NTI0, XSML, XMAX, FIRST
  DATA BI0CS(  1) / -.7660547252839144951081894976243285D-1   /
  DATA BI0CS(  2) / +.1927337953993808269952408750881196D+1   /
  DATA BI0CS(  3) / +.2282644586920301338937029292330415D+0   /
  DATA BI0CS(  4) / +.1304891466707290428079334210691888D-1   /
  DATA BI0CS(  5) / +.4344270900816487451378682681026107D-3   /
  DATA BI0CS(  6) / +.9422657686001934663923171744118766D-5   /
  DATA BI0CS(  7) / +.1434006289510691079962091878179957D-6   /
  DATA BI0CS(  8) / +.1613849069661749069915419719994611D-8   /
  DATA BI0CS(  9) / +.1396650044535669699495092708142522D-10  /
  DATA BI0CS( 10) / +.9579451725505445344627523171893333D-13  /
  DATA BI0CS( 11) / +.5333981859862502131015107744000000D-15  /
  DATA BI0CS( 12) / +.2458716088437470774696785919999999D-17  /
  DATA BI0CS( 13) / +.9535680890248770026944341333333333D-20  /
  DATA BI0CS( 14) / +.3154382039721427336789333333333333D-22  /
  DATA BI0CS( 15) / +.9004564101094637431466666666666666D-25  /
  DATA BI0CS( 16) / +.2240647369123670016000000000000000D-27  /
  DATA BI0CS( 17) / +.4903034603242837333333333333333333D-30  /
  DATA BI0CS( 18) / +.9508172606122666666666666666666666D-33  /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBESI0
  if (FIRST) THEN
     NTI0 = INITDS (BI0CS, 18, 0.1*REAL(D1MACH(3)))
     XSML = SQRT(4.5D0*D1MACH(3))
     XMAX = LOG (D1MACH(2))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 3.0D0) go to 20
!
  DBESI0 = 1.0D0
  if (Y > XSML) DBESI0 = 2.75D0 + DCSEVL (Y*Y/4.5D0-1.D0, BI0CS, &
    NTI0)
  return
!
 20   if (Y  >  XMAX) call XERMSG ('SLATEC', 'DBESI0', &
     'ABS(X) SO BIG I0 OVERFLOWS', 2, 2)
!
  DBESI0 = EXP(Y) * DBSI0E(X)
!
  return
end
  DOUBLE PRECISION FUNCTION DBESI1 (X)
!
!! DBESI1 computes the modified (hyperbolic) Bessel function of the first ...
!  kind of order one.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      DOUBLE PRECISION (BESI1-S, DBESI1-D)
!***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBESI1(X) calculates the double precision modified (hyperbolic)
! Bessel function of the first kind of order one and double precision
! argument X.
!
! Series for BI1        on the interval  0.          to  9.00000E+00
!                                        with weighted error   1.44E-32
!                                         log weighted error  31.84
!                               significant figures required  31.45
!                                    decimal places required  32.46
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DBSI1E, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBESI1
  DOUBLE PRECISION X, BI1CS(17), XMAX, XMIN, XSML, Y, D1MACH, &
    DCSEVL, DBSI1E
  LOGICAL FIRST
  SAVE BI1CS, NTI1, XMIN, XSML, XMAX, FIRST
  DATA BI1CS(  1) / -.19717132610998597316138503218149D-2     /
  DATA BI1CS(  2) / +.40734887667546480608155393652014D+0     /
  DATA BI1CS(  3) / +.34838994299959455866245037783787D-1     /
  DATA BI1CS(  4) / +.15453945563001236038598401058489D-2     /
  DATA BI1CS(  5) / +.41888521098377784129458832004120D-4     /
  DATA BI1CS(  6) / +.76490267648362114741959703966069D-6     /
  DATA BI1CS(  7) / +.10042493924741178689179808037238D-7     /
  DATA BI1CS(  8) / +.99322077919238106481371298054863D-10    /
  DATA BI1CS(  9) / +.76638017918447637275200171681349D-12    /
  DATA BI1CS( 10) / +.47414189238167394980388091948160D-14    /
  DATA BI1CS( 11) / +.24041144040745181799863172032000D-16    /
  DATA BI1CS( 12) / +.10171505007093713649121100799999D-18    /
  DATA BI1CS( 13) / +.36450935657866949458491733333333D-21    /
  DATA BI1CS( 14) / +.11205749502562039344810666666666D-23    /
  DATA BI1CS( 15) / +.29875441934468088832000000000000D-26    /
  DATA BI1CS( 16) / +.69732310939194709333333333333333D-29    /
  DATA BI1CS( 17) / +.14367948220620800000000000000000D-31    /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBESI1
  if (FIRST) THEN
     NTI1 = INITDS (BI1CS, 17, 0.1*REAL(D1MACH(3)))
     XMIN = 2.0D0*D1MACH(1)
     XSML = SQRT(4.5D0*D1MACH(3))
     XMAX = LOG (D1MACH(2))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 3.0D0) go to 20
!
  DBESI1 = 0.D0
  if (Y == 0.D0)  return
!
  if (Y  <=  XMIN) call XERMSG ('SLATEC', 'DBESI1', &
     'ABS(X) SO SMALL I1 UNDERFLOWS', 1, 1)
  if (Y > XMIN) DBESI1 = 0.5D0*X
  if (Y > XSML) DBESI1 = X*(0.875D0 + DCSEVL (Y*Y/4.5D0-1.D0, &
    BI1CS, NTI1))
  return
!
 20   if (Y  >  XMAX) call XERMSG ('SLATEC', 'DBESI1', &
     'ABS(X) SO BIG I1 OVERFLOWS', 2, 2)
!
  DBESI1 = EXP(Y) * DBSI1E(X)
!
  return
end
  DOUBLE PRECISION FUNCTION DBSI0E (X)
!
!! DBSI0E computes the exponentially scaled modified (hyperbolic) Bessel...
!  function of the first kind of order zero.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      DOUBLE PRECISION (BESI0E-S, DBSI0E-D)
!***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
!             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
!             ORDER ZERO, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBSI0E(X) calculates the double precision exponentially scaled
! modified (hyperbolic) Bessel function of the first kind of order
! zero for double precision argument X.  The result is the Bessel
! function I0(X) multiplied by EXP(-ABS(X)).
!
! Series for BI0        on the interval  0.          to  9.00000E+00
!                                        with weighted error   9.51E-34
!                                         log weighted error  33.02
!                               significant figures required  33.31
!                                    decimal places required  33.65
!
! Series for AI0        on the interval  1.25000E-01 to  3.33333E-01
!                                        with weighted error   2.74E-32
!                                         log weighted error  31.56
!                               significant figures required  30.15
!                                    decimal places required  32.39
!
! Series for AI02       on the interval  0.          to  1.25000E-01
!                                        with weighted error   1.97E-32
!                                         log weighted error  31.71
!                               significant figures required  30.15
!                                    decimal places required  32.63
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCSEVL, INITDS
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DBSI0E
  DOUBLE PRECISION X, BI0CS(18), AI0CS(46), AI02CS(69), &
    XSML, Y, D1MACH, DCSEVL
  LOGICAL FIRST
  SAVE BI0CS, AI0CS, AI02CS, NTI0, NTAI0, NTAI02, XSML, FIRST
  DATA BI0CS(  1) / -.7660547252839144951081894976243285D-1   /
  DATA BI0CS(  2) / +.1927337953993808269952408750881196D+1   /
  DATA BI0CS(  3) / +.2282644586920301338937029292330415D+0   /
  DATA BI0CS(  4) / +.1304891466707290428079334210691888D-1   /
  DATA BI0CS(  5) / +.4344270900816487451378682681026107D-3   /
  DATA BI0CS(  6) / +.9422657686001934663923171744118766D-5   /
  DATA BI0CS(  7) / +.1434006289510691079962091878179957D-6   /
  DATA BI0CS(  8) / +.1613849069661749069915419719994611D-8   /
  DATA BI0CS(  9) / +.1396650044535669699495092708142522D-10  /
  DATA BI0CS( 10) / +.9579451725505445344627523171893333D-13  /
  DATA BI0CS( 11) / +.5333981859862502131015107744000000D-15  /
  DATA BI0CS( 12) / +.2458716088437470774696785919999999D-17  /
  DATA BI0CS( 13) / +.9535680890248770026944341333333333D-20  /
  DATA BI0CS( 14) / +.3154382039721427336789333333333333D-22  /
  DATA BI0CS( 15) / +.9004564101094637431466666666666666D-25  /
  DATA BI0CS( 16) / +.2240647369123670016000000000000000D-27  /
  DATA BI0CS( 17) / +.4903034603242837333333333333333333D-30  /
  DATA BI0CS( 18) / +.9508172606122666666666666666666666D-33  /
  DATA AI0CS(  1) / +.7575994494023795942729872037438D-1      /
  DATA AI0CS(  2) / +.7591380810823345507292978733204D-2      /
  DATA AI0CS(  3) / +.4153131338923750501863197491382D-3      /
  DATA AI0CS(  4) / +.1070076463439073073582429702170D-4      /
  DATA AI0CS(  5) / -.7901179979212894660750319485730D-5      /
  DATA AI0CS(  6) / -.7826143501438752269788989806909D-6      /
  DATA AI0CS(  7) / +.2783849942948870806381185389857D-6      /
  DATA AI0CS(  8) / +.8252472600612027191966829133198D-8      /
  DATA AI0CS(  9) / -.1204463945520199179054960891103D-7      /
  DATA AI0CS( 10) / +.1559648598506076443612287527928D-8      /
  DATA AI0CS( 11) / +.2292556367103316543477254802857D-9      /
  DATA AI0CS( 12) / -.1191622884279064603677774234478D-9      /
  DATA AI0CS( 13) / +.1757854916032409830218331247743D-10     /
  DATA AI0CS( 14) / +.1128224463218900517144411356824D-11     /
  DATA AI0CS( 15) / -.1146848625927298877729633876982D-11     /
  DATA AI0CS( 16) / +.2715592054803662872643651921606D-12     /
  DATA AI0CS( 17) / -.2415874666562687838442475720281D-13     /
  DATA AI0CS( 18) / -.6084469888255125064606099639224D-14     /
  DATA AI0CS( 19) / +.3145705077175477293708360267303D-14     /
  DATA AI0CS( 20) / -.7172212924871187717962175059176D-15     /
  DATA AI0CS( 21) / +.7874493403454103396083909603327D-16     /
  DATA AI0CS( 22) / +.1004802753009462402345244571839D-16     /
  DATA AI0CS( 23) / -.7566895365350534853428435888810D-17     /
  DATA AI0CS( 24) / +.2150380106876119887812051287845D-17     /
  DATA AI0CS( 25) / -.3754858341830874429151584452608D-18     /
  DATA AI0CS( 26) / +.2354065842226992576900757105322D-19     /
  DATA AI0CS( 27) / +.1114667612047928530226373355110D-19     /
  DATA AI0CS( 28) / -.5398891884396990378696779322709D-20     /
  DATA AI0CS( 29) / +.1439598792240752677042858404522D-20     /
  DATA AI0CS( 30) / -.2591916360111093406460818401962D-21     /
  DATA AI0CS( 31) / +.2238133183998583907434092298240D-22     /
  DATA AI0CS( 32) / +.5250672575364771172772216831999D-23     /
  DATA AI0CS( 33) / -.3249904138533230784173432285866D-23     /
  DATA AI0CS( 34) / +.9924214103205037927857284710400D-24     /
  DATA AI0CS( 35) / -.2164992254244669523146554299733D-24     /
  DATA AI0CS( 36) / +.3233609471943594083973332991999D-25     /
  DATA AI0CS( 37) / -.1184620207396742489824733866666D-26     /
  DATA AI0CS( 38) / -.1281671853950498650548338687999D-26     /
  DATA AI0CS( 39) / +.5827015182279390511605568853333D-27     /
  DATA AI0CS( 40) / -.1668222326026109719364501503999D-27     /
  DATA AI0CS( 41) / +.3625309510541569975700684800000D-28     /
  DATA AI0CS( 42) / -.5733627999055713589945958399999D-29     /
  DATA AI0CS( 43) / +.3736796722063098229642581333333D-30     /
  DATA AI0CS( 44) / +.1602073983156851963365512533333D-30     /
  DATA AI0CS( 45) / -.8700424864057229884522495999999D-31     /
  DATA AI0CS( 46) / +.2741320937937481145603413333333D-31     /
  DATA AI02CS(  1) / +.5449041101410883160789609622680D-1      /
  DATA AI02CS(  2) / +.3369116478255694089897856629799D-2      /
  DATA AI02CS(  3) / +.6889758346916823984262639143011D-4      /
  DATA AI02CS(  4) / +.2891370520834756482966924023232D-5      /
  DATA AI02CS(  5) / +.2048918589469063741827605340931D-6      /
  DATA AI02CS(  6) / +.2266668990498178064593277431361D-7      /
  DATA AI02CS(  7) / +.3396232025708386345150843969523D-8      /
  DATA AI02CS(  8) / +.4940602388224969589104824497835D-9      /
  DATA AI02CS(  9) / +.1188914710784643834240845251963D-10     /
  DATA AI02CS( 10) / -.3149916527963241364538648629619D-10     /
  DATA AI02CS( 11) / -.1321581184044771311875407399267D-10     /
  DATA AI02CS( 12) / -.1794178531506806117779435740269D-11     /
  DATA AI02CS( 13) / +.7180124451383666233671064293469D-12     /
  DATA AI02CS( 14) / +.3852778382742142701140898017776D-12     /
  DATA AI02CS( 15) / +.1540086217521409826913258233397D-13     /
  DATA AI02CS( 16) / -.4150569347287222086626899720156D-13     /
  DATA AI02CS( 17) / -.9554846698828307648702144943125D-14     /
  DATA AI02CS( 18) / +.3811680669352622420746055355118D-14     /
  DATA AI02CS( 19) / +.1772560133056526383604932666758D-14     /
  DATA AI02CS( 20) / -.3425485619677219134619247903282D-15     /
  DATA AI02CS( 21) / -.2827623980516583484942055937594D-15     /
  DATA AI02CS( 22) / +.3461222867697461093097062508134D-16     /
  DATA AI02CS( 23) / +.4465621420296759999010420542843D-16     /
  DATA AI02CS( 24) / -.4830504485944182071255254037954D-17     /
  DATA AI02CS( 25) / -.7233180487874753954562272409245D-17     /
  DATA AI02CS( 26) / +.9921475412173698598880460939810D-18     /
  DATA AI02CS( 27) / +.1193650890845982085504399499242D-17     /
  DATA AI02CS( 28) / -.2488709837150807235720544916602D-18     /
  DATA AI02CS( 29) / -.1938426454160905928984697811326D-18     /
  DATA AI02CS( 30) / +.6444656697373443868783019493949D-19     /
  DATA AI02CS( 31) / +.2886051596289224326481713830734D-19     /
  DATA AI02CS( 32) / -.1601954907174971807061671562007D-19     /
  DATA AI02CS( 33) / -.3270815010592314720891935674859D-20     /
  DATA AI02CS( 34) / +.3686932283826409181146007239393D-20     /
  DATA AI02CS( 35) / +.1268297648030950153013595297109D-22     /
  DATA AI02CS( 36) / -.7549825019377273907696366644101D-21     /
  DATA AI02CS( 37) / +.1502133571377835349637127890534D-21     /
  DATA AI02CS( 38) / +.1265195883509648534932087992483D-21     /
  DATA AI02CS( 39) / -.6100998370083680708629408916002D-22     /
  DATA AI02CS( 40) / -.1268809629260128264368720959242D-22     /
  DATA AI02CS( 41) / +.1661016099890741457840384874905D-22     /
  DATA AI02CS( 42) / -.1585194335765885579379705048814D-23     /
  DATA AI02CS( 43) / -.3302645405968217800953817667556D-23     /
  DATA AI02CS( 44) / +.1313580902839239781740396231174D-23     /
  DATA AI02CS( 45) / +.3689040246671156793314256372804D-24     /
  DATA AI02CS( 46) / -.4210141910461689149219782472499D-24     /
  DATA AI02CS( 47) / +.4791954591082865780631714013730D-25     /
  DATA AI02CS( 48) / +.8459470390221821795299717074124D-25     /
  DATA AI02CS( 49) / -.4039800940872832493146079371810D-25     /
  DATA AI02CS( 50) / -.6434714653650431347301008504695D-26     /
  DATA AI02CS( 51) / +.1225743398875665990344647369905D-25     /
  DATA AI02CS( 52) / -.2934391316025708923198798211754D-26     /
  DATA AI02CS( 53) / -.1961311309194982926203712057289D-26     /
  DATA AI02CS( 54) / +.1503520374822193424162299003098D-26     /
  DATA AI02CS( 55) / -.9588720515744826552033863882069D-28     /
  DATA AI02CS( 56) / -.3483339380817045486394411085114D-27     /
  DATA AI02CS( 57) / +.1690903610263043673062449607256D-27     /
  DATA AI02CS( 58) / +.1982866538735603043894001157188D-28     /
  DATA AI02CS( 59) / -.5317498081491816214575830025284D-28     /
  DATA AI02CS( 60) / +.1803306629888392946235014503901D-28     /
  DATA AI02CS( 61) / +.6213093341454893175884053112422D-29     /
  DATA AI02CS( 62) / -.7692189292772161863200728066730D-29     /
  DATA AI02CS( 63) / +.1858252826111702542625560165963D-29     /
  DATA AI02CS( 64) / +.1237585142281395724899271545541D-29     /
  DATA AI02CS( 65) / -.1102259120409223803217794787792D-29     /
  DATA AI02CS( 66) / +.1886287118039704490077874479431D-30     /
  DATA AI02CS( 67) / +.2160196872243658913149031414060D-30     /
  DATA AI02CS( 68) / -.1605454124919743200584465949655D-30     /
  DATA AI02CS( 69) / +.1965352984594290603938848073318D-31     /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBSI0E
  if (FIRST) THEN
     ETA = 0.1*REAL(D1MACH(3))
     NTI0 = INITDS (BI0CS, 18, ETA)
     NTAI0 = INITDS (AI0CS, 46, ETA)
     NTAI02 = INITDS (AI02CS, 69, ETA)
     XSML = SQRT(4.5D0*D1MACH(3))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 3.0D0) go to 20
!
  DBSI0E = 1.0D0 - X
  if (Y > XSML) DBSI0E = EXP(-Y) * (2.75D0 + &
    DCSEVL (Y*Y/4.5D0-1.D0, BI0CS, NTI0) )
  return
!
 20   if (Y <= 8.D0) DBSI0E = (0.375D0 + DCSEVL ((48.D0/Y-11.D0)/5.D0, &
    AI0CS, NTAI0))/SQRT(Y)
  if (Y > 8.D0) DBSI0E = (0.375D0 + DCSEVL (16.D0/Y-1.D0, AI02CS, &
    NTAI02))/SQRT(Y)
!
  return
end
  DOUBLE PRECISION FUNCTION DBSI1E (X)
!
!! DBSI1E computes the exponentially scaled modified (hyperbolic) Bessel...
!  function of the first kind of order one.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      DOUBLE PRECISION (BESI1E-S, DBSI1E-D)
!***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
!             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
!             ORDER ONE, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBSI1E(X) calculates the double precision exponentially scaled
! modified (hyperbolic) Bessel function of the first kind of order
! one for double precision argument X.  The result is I1(X)
! multiplied by EXP(-ABS(X)).
!
! Series for BI1        on the interval  0.          to  9.00000E+00
!                                        with weighted error   1.44E-32
!                                         log weighted error  31.84
!                               significant figures required  31.45
!                                    decimal places required  32.46
!
! Series for AI1        on the interval  1.25000E-01 to  3.33333E-01
!                                        with weighted error   2.81E-32
!                                         log weighted error  31.55
!                               significant figures required  29.93
!                                    decimal places required  32.38
!
! Series for AI12       on the interval  0.          to  1.25000E-01
!                                        with weighted error   1.83E-32
!                                         log weighted error  31.74
!                               significant figures required  29.97
!                                    decimal places required  32.66
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBSI1E
  DOUBLE PRECISION X, BI1CS(17), AI1CS(46), AI12CS(69), XMIN, &
    XSML, Y, D1MACH, DCSEVL
  LOGICAL FIRST
  SAVE BI1CS, AI1CS, AI12CS, NTI1, NTAI1, NTAI12, XMIN, XSML, &
    FIRST
  DATA BI1CS(  1) / -.19717132610998597316138503218149D-2     /
  DATA BI1CS(  2) / +.40734887667546480608155393652014D+0     /
  DATA BI1CS(  3) / +.34838994299959455866245037783787D-1     /
  DATA BI1CS(  4) / +.15453945563001236038598401058489D-2     /
  DATA BI1CS(  5) / +.41888521098377784129458832004120D-4     /
  DATA BI1CS(  6) / +.76490267648362114741959703966069D-6     /
  DATA BI1CS(  7) / +.10042493924741178689179808037238D-7     /
  DATA BI1CS(  8) / +.99322077919238106481371298054863D-10    /
  DATA BI1CS(  9) / +.76638017918447637275200171681349D-12    /
  DATA BI1CS( 10) / +.47414189238167394980388091948160D-14    /
  DATA BI1CS( 11) / +.24041144040745181799863172032000D-16    /
  DATA BI1CS( 12) / +.10171505007093713649121100799999D-18    /
  DATA BI1CS( 13) / +.36450935657866949458491733333333D-21    /
  DATA BI1CS( 14) / +.11205749502562039344810666666666D-23    /
  DATA BI1CS( 15) / +.29875441934468088832000000000000D-26    /
  DATA BI1CS( 16) / +.69732310939194709333333333333333D-29    /
  DATA BI1CS( 17) / +.14367948220620800000000000000000D-31    /
  DATA AI1CS(  1) / -.2846744181881478674100372468307D-1      /
  DATA AI1CS(  2) / -.1922953231443220651044448774979D-1      /
  DATA AI1CS(  3) / -.6115185857943788982256249917785D-3      /
  DATA AI1CS(  4) / -.2069971253350227708882823777979D-4      /
  DATA AI1CS(  5) / +.8585619145810725565536944673138D-5      /
  DATA AI1CS(  6) / +.1049498246711590862517453997860D-5      /
  DATA AI1CS(  7) / -.2918338918447902202093432326697D-6      /
  DATA AI1CS(  8) / -.1559378146631739000160680969077D-7      /
  DATA AI1CS(  9) / +.1318012367144944705525302873909D-7      /
  DATA AI1CS( 10) / -.1448423418183078317639134467815D-8      /
  DATA AI1CS( 11) / -.2908512243993142094825040993010D-9      /
  DATA AI1CS( 12) / +.1266388917875382387311159690403D-9      /
  DATA AI1CS( 13) / -.1664947772919220670624178398580D-10     /
  DATA AI1CS( 14) / -.1666653644609432976095937154999D-11     /
  DATA AI1CS( 15) / +.1242602414290768265232168472017D-11     /
  DATA AI1CS( 16) / -.2731549379672432397251461428633D-12     /
  DATA AI1CS( 17) / +.2023947881645803780700262688981D-13     /
  DATA AI1CS( 18) / +.7307950018116883636198698126123D-14     /
  DATA AI1CS( 19) / -.3332905634404674943813778617133D-14     /
  DATA AI1CS( 20) / +.7175346558512953743542254665670D-15     /
  DATA AI1CS( 21) / -.6982530324796256355850629223656D-16     /
  DATA AI1CS( 22) / -.1299944201562760760060446080587D-16     /
  DATA AI1CS( 23) / +.8120942864242798892054678342860D-17     /
  DATA AI1CS( 24) / -.2194016207410736898156266643783D-17     /
  DATA AI1CS( 25) / +.3630516170029654848279860932334D-18     /
  DATA AI1CS( 26) / -.1695139772439104166306866790399D-19     /
  DATA AI1CS( 27) / -.1288184829897907807116882538222D-19     /
  DATA AI1CS( 28) / +.5694428604967052780109991073109D-20     /
  DATA AI1CS( 29) / -.1459597009090480056545509900287D-20     /
  DATA AI1CS( 30) / +.2514546010675717314084691334485D-21     /
  DATA AI1CS( 31) / -.1844758883139124818160400029013D-22     /
  DATA AI1CS( 32) / -.6339760596227948641928609791999D-23     /
  DATA AI1CS( 33) / +.3461441102031011111108146626560D-23     /
  DATA AI1CS( 34) / -.1017062335371393547596541023573D-23     /
  DATA AI1CS( 35) / +.2149877147090431445962500778666D-24     /
  DATA AI1CS( 36) / -.3045252425238676401746206173866D-25     /
  DATA AI1CS( 37) / +.5238082144721285982177634986666D-27     /
  DATA AI1CS( 38) / +.1443583107089382446416789503999D-26     /
  DATA AI1CS( 39) / -.6121302074890042733200670719999D-27     /
  DATA AI1CS( 40) / +.1700011117467818418349189802666D-27     /
  DATA AI1CS( 41) / -.3596589107984244158535215786666D-28     /
  DATA AI1CS( 42) / +.5448178578948418576650513066666D-29     /
  DATA AI1CS( 43) / -.2731831789689084989162564266666D-30     /
  DATA AI1CS( 44) / -.1858905021708600715771903999999D-30     /
  DATA AI1CS( 45) / +.9212682974513933441127765333333D-31     /
  DATA AI1CS( 46) / -.2813835155653561106370833066666D-31     /
  DATA AI12CS(  1) / +.2857623501828012047449845948469D-1      /
  DATA AI12CS(  2) / -.9761097491361468407765164457302D-2      /
  DATA AI12CS(  3) / -.1105889387626237162912569212775D-3      /
  DATA AI12CS(  4) / -.3882564808877690393456544776274D-5      /
  DATA AI12CS(  5) / -.2512236237870208925294520022121D-6      /
  DATA AI12CS(  6) / -.2631468846889519506837052365232D-7      /
  DATA AI12CS(  7) / -.3835380385964237022045006787968D-8      /
  DATA AI12CS(  8) / -.5589743462196583806868112522229D-9      /
  DATA AI12CS(  9) / -.1897495812350541234498925033238D-10     /
  DATA AI12CS( 10) / +.3252603583015488238555080679949D-10     /
  DATA AI12CS( 11) / +.1412580743661378133163366332846D-10     /
  DATA AI12CS( 12) / +.2035628544147089507224526136840D-11     /
  DATA AI12CS( 13) / -.7198551776245908512092589890446D-12     /
  DATA AI12CS( 14) / -.4083551111092197318228499639691D-12     /
  DATA AI12CS( 15) / -.2101541842772664313019845727462D-13     /
  DATA AI12CS( 16) / +.4272440016711951354297788336997D-13     /
  DATA AI12CS( 17) / +.1042027698412880276417414499948D-13     /
  DATA AI12CS( 18) / -.3814403072437007804767072535396D-14     /
  DATA AI12CS( 19) / -.1880354775510782448512734533963D-14     /
  DATA AI12CS( 20) / +.3308202310920928282731903352405D-15     /
  DATA AI12CS( 21) / +.2962628997645950139068546542052D-15     /
  DATA AI12CS( 22) / -.3209525921993423958778373532887D-16     /
  DATA AI12CS( 23) / -.4650305368489358325571282818979D-16     /
  DATA AI12CS( 24) / +.4414348323071707949946113759641D-17     /
  DATA AI12CS( 25) / +.7517296310842104805425458080295D-17     /
  DATA AI12CS( 26) / -.9314178867326883375684847845157D-18     /
  DATA AI12CS( 27) / -.1242193275194890956116784488697D-17     /
  DATA AI12CS( 28) / +.2414276719454848469005153902176D-18     /
  DATA AI12CS( 29) / +.2026944384053285178971922860692D-18     /
  DATA AI12CS( 30) / -.6394267188269097787043919886811D-19     /
  DATA AI12CS( 31) / -.3049812452373095896084884503571D-19     /
  DATA AI12CS( 32) / +.1612841851651480225134622307691D-19     /
  DATA AI12CS( 33) / +.3560913964309925054510270904620D-20     /
  DATA AI12CS( 34) / -.3752017947936439079666828003246D-20     /
  DATA AI12CS( 35) / -.5787037427074799345951982310741D-22     /
  DATA AI12CS( 36) / +.7759997511648161961982369632092D-21     /
  DATA AI12CS( 37) / -.1452790897202233394064459874085D-21     /
  DATA AI12CS( 38) / -.1318225286739036702121922753374D-21     /
  DATA AI12CS( 39) / +.6116654862903070701879991331717D-22     /
  DATA AI12CS( 40) / +.1376279762427126427730243383634D-22     /
  DATA AI12CS( 41) / -.1690837689959347884919839382306D-22     /
  DATA AI12CS( 42) / +.1430596088595433153987201085385D-23     /
  DATA AI12CS( 43) / +.3409557828090594020405367729902D-23     /
  DATA AI12CS( 44) / -.1309457666270760227845738726424D-23     /
  DATA AI12CS( 45) / -.3940706411240257436093521417557D-24     /
  DATA AI12CS( 46) / +.4277137426980876580806166797352D-24     /
  DATA AI12CS( 47) / -.4424634830982606881900283123029D-25     /
  DATA AI12CS( 48) / -.8734113196230714972115309788747D-25     /
  DATA AI12CS( 49) / +.4045401335683533392143404142428D-25     /
  DATA AI12CS( 50) / +.7067100658094689465651607717806D-26     /
  DATA AI12CS( 51) / -.1249463344565105223002864518605D-25     /
  DATA AI12CS( 52) / +.2867392244403437032979483391426D-26     /
  DATA AI12CS( 53) / +.2044292892504292670281779574210D-26     /
  DATA AI12CS( 54) / -.1518636633820462568371346802911D-26     /
  DATA AI12CS( 55) / +.8110181098187575886132279107037D-28     /
  DATA AI12CS( 56) / +.3580379354773586091127173703270D-27     /
  DATA AI12CS( 57) / -.1692929018927902509593057175448D-27     /
  DATA AI12CS( 58) / -.2222902499702427639067758527774D-28     /
  DATA AI12CS( 59) / +.5424535127145969655048600401128D-28     /
  DATA AI12CS( 60) / -.1787068401578018688764912993304D-28     /
  DATA AI12CS( 61) / -.6565479068722814938823929437880D-29     /
  DATA AI12CS( 62) / +.7807013165061145280922067706839D-29     /
  DATA AI12CS( 63) / -.1816595260668979717379333152221D-29     /
  DATA AI12CS( 64) / -.1287704952660084820376875598959D-29     /
  DATA AI12CS( 65) / +.1114548172988164547413709273694D-29     /
  DATA AI12CS( 66) / -.1808343145039336939159368876687D-30     /
  DATA AI12CS( 67) / -.2231677718203771952232448228939D-30     /
  DATA AI12CS( 68) / +.1619029596080341510617909803614D-30     /
  DATA AI12CS( 69) / -.1834079908804941413901308439210D-31     /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBSI1E
  if (FIRST) THEN
     ETA = 0.1*REAL(D1MACH(3))
     NTI1 = INITDS (BI1CS, 17, ETA)
     NTAI1 = INITDS (AI1CS, 46, ETA)
     NTAI12 = INITDS (AI12CS, 69, ETA)
!
     XMIN = 2.0D0*D1MACH(1)
     XSML = SQRT(4.5D0*D1MACH(3))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 3.0D0) go to 20
!
  DBSI1E = 0.0D0
  if (Y == 0.D0)  return
!
  if (Y  <=  XMIN) call XERMSG ('SLATEC', 'DBSI1E', &
     'ABS(X) SO SMALL I1 UNDERFLOWS', 1, 1)
  if (Y > XMIN) DBSI1E = 0.5D0*X
  if (Y > XSML) DBSI1E = X*(0.875D0 + DCSEVL (Y*Y/4.5D0-1.D0, &
    BI1CS, NTI1) )
  DBSI1E = EXP(-Y) * DBSI1E
  return
!
 20   if (Y <= 8.D0) DBSI1E = (0.375D0 + DCSEVL ((48.D0/Y-11.D0)/5.D0, &
    AI1CS, NTAI1))/SQRT(Y)
  if (Y > 8.D0) DBSI1E = (0.375D0 + DCSEVL (16.D0/Y-1.D0, AI12CS, &
    NTAI12))/SQRT(Y)
  DBSI1E = SIGN (DBSI1E, X)
!
  return
end
subroutine DBSKES (XNU, X, NIN, BKE)
!
!! DBSKES computes a sequence of exponentially scaled modified Bessel ...
!  functions of the third kind of fractional order.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B3
!***TYPE      DOUBLE PRECISION (BESKES-S, DBSKES-D)
!***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, FRACTIONAL ORDER,
!             MODIFIED BESSEL FUNCTION, SEQUENCE OF BESSEL FUNCTIONS,
!             SPECIAL FUNCTIONS, THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBSKES(XNU,X,NIN,BKE) computes a double precision sequence
! of exponentially scaled modified Bessel functions
! of the third kind of order XNU + I at X, where X  >  0,
! XNU lies in (-1,1), and I = 0, 1, ... , NIN - 1, if NIN is positive
! and I = 0, -1, ... , NIN + 1, if NIN is negative.  On return, the
! vector BKE(.) contains the results at X for order starting at XNU.
! XNU, X, and BKE are double precision.  NIN is integer.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, D9KNUS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBSKES
  DOUBLE PRECISION XNU, X, BKE(*), BKNU1, V, VINCR, VEND, ALNBIG, &
    D1MACH, DIRECT
  SAVE ALNBIG
  DATA ALNBIG / 0.D0 /
!***FIRST EXECUTABLE STATEMENT  DBSKES
  if (ALNBIG == 0.D0) ALNBIG = LOG (D1MACH(2))
!
  V = ABS(XNU)
  N = ABS(NIN)
!
  if (V  >=  1.D0) call XERMSG ('SLATEC', 'DBSKES', &
     'ABS(XNU) MUST BE LT 1', 2, 2)
  if (X  <=  0.D0) call XERMSG ('SLATEC', 'DBSKES', 'X IS LE 0', 3, &
     2)
  if (N  ==  0) call XERMSG ('SLATEC', 'DBSKES', &
     'N THE NUMBER IN THE SEQUENCE IS 0', 4, 2)
!
  call D9KNUS (V, X, BKE(1), BKNU1, ISWTCH)
  if (N == 1) RETURN
!
  VINCR = SIGN (1.0, REAL(NIN))
  DIRECT = VINCR
  if (XNU /= 0.D0) DIRECT = VINCR*SIGN(1.D0, XNU)
  if (ISWTCH  ==  1 .AND. DIRECT  >  0.) call XERMSG ('SLATEC', &
     'DBSKES', 'X SO SMALL BESSEL K-SUB-XNU+1 OVERFLOWS', 5, 2)
  BKE(2) = BKNU1
!
  if (DIRECT < 0.) call D9KNUS (ABS(XNU+VINCR), X, BKE(2), BKNU1, &
    ISWTCH)
  if (N == 2) RETURN
!
  VEND = ABS (XNU+NIN) - 1.0D0
  if ((VEND-.5D0)*LOG(VEND)+0.27D0-VEND*(LOG(X)-.694D0)  >  &
     ALNBIG) call XERMSG ('SLATEC', 'DBSKES', &
        'X SO SMALL OR ABS(NU) SO BIG THAT BESSEL K-SUB-NU ' // &
        'OVERFLOWS', 5, 2)
!
  V = XNU
  DO 10 I=3,N
    V = V + VINCR
    BKE(I) = 2.0D0*V*BKE(I-1)/X + BKE(I-2)
 10   CONTINUE
!
  return
end
subroutine D9KNUS (XNU, X, BKNU, BKNU1, ISWTCH)
!
!! D9KNUS computes Bessel functions EXP(X)*K-SUB-XNU(X) and ...
!  EXP(X)*K-SUB-XNU+1(X) for 0.0  <=  XNU  <  1.0.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B3
!***TYPE      DOUBLE PRECISION (R9KNUS-S, D9KNUS-D)
!***KEYWORDS  BESSEL FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute Bessel functions EXP(X) * K-sub-XNU (X)  and
! EXP(X) * K-sub-XNU+1 (X) for 0.0  <=  XNU  <  1.0 .
!
! Series for C0K        on the interval  0.          to  2.50000E-01
!                                        with weighted error   2.16E-32
!                                         log weighted error  31.67
!                               significant figures required  30.86
!                                    decimal places required  32.40
!
! Series for ZNU1       on the interval -7.00000E-01 to  0.
!                                        with weighted error   2.45E-33
!                                         log weighted error  32.61
!                               significant figures required  31.85
!                                    decimal places required  33.26
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCSEVL, DGAMMA, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!   900727  Added EXTERNAL statement.  (WRB)
!   920618  Removed space from variable names.  (RWC, WRB)
!***END PROLOGUE  D9KNUS
  DOUBLE PRECISION XNU, X, BKNU, BKNU1, ALPHA(32), BETA(32), A(32), &
    C0KCS(29), ZNU1CS(20), ALNZ, ALN2, A0, BKNUD, BKNU0, &
    B0, C0, EULER, EXPX, P1, P2, P3, QQ, RESULT, SQPI2, SQRTX, V, &
    VLNZ, XI, XMU, XNUSML, XSML, X2N, X2TOV, Z, ZTOV, ALNSML, &
    ALNBIG
  REAL ALNEPS
  DOUBLE PRECISION D1MACH, DCSEVL, DGAMMA
  LOGICAL FIRST
  EXTERNAL DGAMMA
  SAVE C0KCS, ZNU1CS, EULER, SQPI2, ALN2, NTC0K, &
   NTZNU1, XNUSML, XSML, ALNSML, ALNBIG, ALNEPS, FIRST
  DATA C0KCS(  1) / +.60183057242626108387577445180329D-1     /
  DATA C0KCS(  2) / -.15364871433017286092959755943124D+0     /
  DATA C0KCS(  3) / -.11751176008210492040068229226213D-1     /
  DATA C0KCS(  4) / -.85248788891979509827048401550987D-3     /
  DATA C0KCS(  5) / -.61329838767496791874098176922111D-4     /
  DATA C0KCS(  6) / -.44052281245510444562679889548505D-5     /
  DATA C0KCS(  7) / -.31631246728384488192915445892199D-6     /
  DATA C0KCS(  8) / -.22710719382899588330673771793396D-7     /
  DATA C0KCS(  9) / -.16305644608077609552274620515360D-8     /
  DATA C0KCS( 10) / -.11706939299414776568756044043130D-9     /
  DATA C0KCS( 11) / -.84052063786464437174546593413792D-11    /
  DATA C0KCS( 12) / -.60346670118979991487096050737198D-12    /
  DATA C0KCS( 13) / -.43326960335681371952045997366903D-13    /
  DATA C0KCS( 14) / -.31107358030203546214634697772237D-14    /
  DATA C0KCS( 15) / -.22334078226736982254486133409840D-15    /
  DATA C0KCS( 16) / -.16035146716864226300635791528610D-16    /
  DATA C0KCS( 17) / -.11512717363666556196035697705305D-17    /
  DATA C0KCS( 18) / -.82657591746836959105169479089258D-19    /
  DATA C0KCS( 19) / -.59345480806383948172333436695984D-20    /
  DATA C0KCS( 20) / -.42608138196467143926499613023976D-21    /
  DATA C0KCS( 21) / -.30591266864812876299263698370542D-22    /
  DATA C0KCS( 22) / -.21963541426734575224975501815516D-23    /
  DATA C0KCS( 23) / -.15769113261495836071105750684760D-24    /
  DATA C0KCS( 24) / -.11321713935950320948757731048056D-25    /
  DATA C0KCS( 25) / -.81286248834598404082792349714433D-27    /
  DATA C0KCS( 26) / -.58360900893453226552829349315949D-28    /
  DATA C0KCS( 27) / -.41901241623610922519452337780905D-29    /
  DATA C0KCS( 28) / -.30083737960206435069530504212862D-30    /
  DATA C0KCS( 29) / -.21599152067808647728342168089832D-31    /
  DATA ZNU1CS(  1) / +.203306756994191729674444001216911D+0    /
  DATA ZNU1CS(  2) / +.140077933413219771062943670790563D+0    /
  DATA ZNU1CS(  3) / +.791679696100161352840972241972320D-2    /
  DATA ZNU1CS(  4) / +.339801182532104045352930092205750D-3    /
  DATA ZNU1CS(  5) / +.117419756889893366664507228352690D-4    /
  DATA ZNU1CS(  6) / +.339357570612261680333825865475121D-6    /
  DATA ZNU1CS(  7) / +.842594176976219910194629891264803D-8    /
  DATA ZNU1CS(  8) / +.183336677024850089184748150900090D-9    /
  DATA ZNU1CS(  9) / +.354969844704416310863007064469557D-11   /
  DATA ZNU1CS( 10) / +.619032496469887332205244342078407D-13   /
  DATA ZNU1CS( 11) / +.981964535680439424960346115456527D-15   /
  DATA ZNU1CS( 12) / +.142851314396490474211473563005985D-16   /
  DATA ZNU1CS( 13) / +.191894921887825298966162467488436D-18   /
  DATA ZNU1CS( 14) / +.239430979739498914162313140597128D-20   /
  DATA ZNU1CS( 15) / +.278890246815347354835870465474995D-22   /
  DATA ZNU1CS( 16) / +.304606650633033442582845214092865D-24   /
  DATA ZNU1CS( 17) / +.313173237042191815771564260932089D-26   /
  DATA ZNU1CS( 18) / +.304133098987854951645174908005034D-28   /
  DATA ZNU1CS( 19) / +.279840384636833084343185097659733D-30   /
  DATA ZNU1CS( 20) / +.244637186274497596485238794922666D-32   /
  DATA EULER / 0.57721566490153286060651209008240D0 /
  DATA SQPI2 / +1.2533141373155002512078826424055D0      /
  DATA ALN2 / 0.69314718055994530941723212145818D0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  D9KNUS
  if (FIRST) THEN
     ETA = 0.1D0*D1MACH(3)
     NTC0K = INITDS (C0KCS, 29, ETA)
     NTZNU1 = INITDS (ZNU1CS, 20, ETA)
!
     XNUSML = SQRT(D1MACH(3)/8.D0)
     XSML = 0.1D0*D1MACH(3)
     ALNSML = LOG (D1MACH(1))
     ALNBIG = LOG (D1MACH(2))
     ALNEPS = LOG (0.1D0*D1MACH(3))
  end if
  FIRST = .FALSE.
!
  if (XNU  <  0.D0 .OR. XNU  >=  1.D0) call XERMSG ('SLATEC', &
     'D9KNUS', 'XNU MUST BE GE 0 AND LT 1', 1, 2)
  if (X  <=  0.) call XERMSG ('SLATEC', 'D9KNUS', 'X MUST BE GT 0', &
     2, 2)
!
  ISWTCH = 0
  if (X > 2.0D0) go to 50
!
! X IS SMALL.  COMPUTE K-SUB-XNU (X) AND THE DERIVATIVE OF K-SUB-XNU (X)
! THEN FIND K-SUB-XNU+1 (X).  XNU IS REDUCED TO THE INTERVAL (-.5,+.5)
! THEN TO (0., .5), BECAUSE K OF NEGATIVE ORDER (-NU) = K OF POSITIVE
! ORDER (+NU).
!
  V = XNU
  if (XNU > 0.5D0) V = 1.0D0 - XNU
!
! CAREFULLY FIND (X/2)**XNU AND Z**XNU WHERE Z = X*X/4.
  ALNZ = 2.D0 * (LOG(X) - ALN2)
!
  if (X > XNU) go to 20
  if (-0.5D0*XNU*ALNZ-ALN2-LOG(XNU)  >  ALNBIG) call XERMSG &
     ('SLATEC', 'D9KNUS', 'X SO SMALL BESSEL K-SUB-XNU OVERFLOWS', &
     3, 2)
!
 20   VLNZ = V*ALNZ
  X2TOV = EXP (0.5D0*VLNZ)
  ZTOV = 0.0D0
  if (VLNZ > ALNSML) ZTOV = X2TOV**2
!
  A0 = 0.5D0*DGAMMA(1.0D0+V)
  B0 = 0.5D0*DGAMMA(1.0D0-V)
  C0 = -EULER
  if (ZTOV > 0.5D0 .AND. V > XNUSML) C0 = -0.75D0 + &
    DCSEVL ((8.0D0*V)*V-1.0D0, C0KCS, NTC0K)
!
  if (ZTOV <= 0.5D0) ALPHA(1) = (A0-ZTOV*B0)/V
  if (ZTOV > 0.5D0) ALPHA(1) = C0 - ALNZ*(0.75D0 + &
    DCSEVL (VLNZ/0.35D0+1.0D0, ZNU1CS, NTZNU1))*B0
  BETA(1) = -0.5D0*(A0+ZTOV*B0)
!
  Z = 0.0D0
  if (X > XSML) Z = 0.25D0*X*X
  NTERMS = MAX (2.0, 11.0+(8.*REAL(ALNZ)-25.19-ALNEPS) &
    /(4.28-REAL(ALNZ)))
  DO 30 I=2,NTERMS
    XI = I - 1
    A0 = A0/(XI*(XI-V))
    B0 = B0/(XI*(XI+V))
    ALPHA(I) = (ALPHA(I-1)+2.0D0*XI*A0)/(XI*(XI+V))
    BETA(I) = (XI-0.5D0*V)*ALPHA(I) - ZTOV*B0
 30   CONTINUE
!
  BKNU = ALPHA(NTERMS)
  BKNUD = BETA(NTERMS)
  DO 40 II=2,NTERMS
    I = NTERMS + 1 - II
    BKNU = ALPHA(I) + BKNU*Z
    BKNUD = BETA(I) + BKNUD*Z
 40   CONTINUE
!
  EXPX = EXP(X)
  BKNU = EXPX*BKNU/X2TOV
!
  if (-0.5D0*(XNU+1.D0)*ALNZ-2.0D0*ALN2 > ALNBIG) ISWTCH = 1
  if (ISWTCH == 1) RETURN
  BKNUD = EXPX*BKNUD*2.0D0/(X2TOV*X)
!
  if (XNU <= 0.5D0) BKNU1 = V*BKNU/X - BKNUD
  if (XNU <= 0.5D0) RETURN
!
  BKNU0 = BKNU
  BKNU = -V*BKNU/X - BKNUD
  BKNU1 = 2.0D0*XNU*BKNU/X + BKNU0
  return
!
! X IS LARGE.  FIND K-SUB-XNU (X) AND K-SUB-XNU+1 (X) WITH Y. L. LUKE-S
! RATIONAL EXPANSION.
!
 50   SQRTX = SQRT(X)
  if (X > 1.0D0/XSML) go to 90
  AN = -0.60 - 1.02/REAL(X)
  BN = -0.27 - 0.53/REAL(X)
  NTERMS = MIN (32, MAX1 (3.0, AN+BN*ALNEPS))
!
  DO 80 INU=1,2
    XMU = 0.D0
    if (INU == 1 .AND. XNU > XNUSML) XMU = (4.0D0*XNU)*XNU
    if (INU == 2) XMU = 4.0D0*(ABS(XNU)+1.D0)**2
!
    A(1) = 1.0D0 - XMU
    A(2) = 9.0D0 - XMU
    A(3) = 25.0D0 - XMU
    if (A(2) == 0.D0) RESULT = SQPI2*(16.D0*X+XMU+7.D0) / &
      (16.D0*X*SQRTX)
    if (A(2) == 0.D0) go to 70
!
    ALPHA(1) = 1.0D0
    ALPHA(2) = (16.D0*X+A(2))/A(2)
    ALPHA(3) = ((768.D0*X+48.D0*A(3))*X + A(2)*A(3))/(A(2)*A(3))
!
    BETA(1) = 1.0D0
    BETA(2) = (16.D0*X+(XMU+7.D0))/A(2)
    BETA(3) = ((768.D0*X+48.D0*(XMU+23.D0))*X + &
      ((XMU+62.D0)*XMU+129.D0))/(A(2)*A(3))
!
    if (NTERMS < 4) go to 65
    DO 60 I=4,NTERMS
      N = I - 1
      X2N = 2*N - 1
!
      A(I) = (X2N+2.D0)**2 - XMU
      QQ = 16.D0*X2N/A(I)
      P1 = -X2N*((12*N*N-20*N)-A(1))/((X2N-2.D0)*A(I)) &
        - QQ*X
      P2 = ((12*N*N-28*N+8)-A(1))/A(I) - QQ*X
      P3 = -X2N*A(I-3)/((X2N-2.D0)*A(I))
!
      ALPHA(I) = -P1*ALPHA(I-1) - P2*ALPHA(I-2) - P3*ALPHA(I-3)
      BETA(I) = -P1*BETA(I-1) - P2*BETA(I-2) - P3*BETA(I-3)
 60     CONTINUE
!
 65     RESULT = SQPI2*BETA(NTERMS)/(SQRTX*ALPHA(NTERMS))
!
 70     if (INU == 1) BKNU = RESULT
    if (INU == 2) BKNU1 = RESULT
 80   CONTINUE
  return
!
 90   BKNU = SQPI2/SQRTX
  BKNU1 = BKNU
  return
!
end
  DOUBLE PRECISION FUNCTION DCSEVL (X, CS, N)
!
!! DCSEVL evaluates a Chebyshev series.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C3A2
!***TYPE      DOUBLE PRECISION (CSEVL-S, DCSEVL-D)
!***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!  Evaluate the N-term Chebyshev series CS at X.  Adapted from
!  a method presented in the paper by Broucke referenced below.
!
!       Input Arguments --
!  X    value at which the series is to be evaluated.
!  CS   array of N terms of a Chebyshev series.  In evaluating
!       CS, only half the first coefficient is summed.
!  N    number of terms in array CS.
!
!***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
!                 Chebyshev series, Algorithm 446, Communications of
!                 the A.C.M. 16, (1973) pp. 254-256.
!               L. Fox and I. B. Parker, Chebyshev Polynomials in
!                 Numerical Analysis, Oxford University Press, 1968,
!                 page 56.
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900329  Prologued revised extensively and code rewritten to allow
!           X to be slightly outside interval (-1,+1).  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DCSEVL
  DOUBLE PRECISION B0, B1, B2, CS(*), ONEPL, TWOX, X, D1MACH
  LOGICAL FIRST
  SAVE FIRST, ONEPL
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DCSEVL
  if (FIRST) ONEPL = 1.0D0 + D1MACH(4)
  FIRST = .FALSE.
  if (N  <  1) call XERMSG ('SLATEC', 'DCSEVL', &
     'NUMBER OF TERMS  <=  0', 2, 2)
  if (N  >  1000) call XERMSG ('SLATEC', 'DCSEVL', &
     'NUMBER OF TERMS  >  1000', 3, 2)
  if (ABS(X)  >  ONEPL) call XERMSG ('SLATEC', 'DCSEVL', &
     'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1)
!
  B1 = 0.0D0
  B0 = 0.0D0
  TWOX = 2.0D0*X
  DO 10 I = 1,N
     B2 = B1
     B1 = B0
     NI = N + 1 - I
     B0 = TWOX*B1 - B2 + CS(NI)
   10 CONTINUE
!
  DCSEVL = 0.5D0*(B0-B2)
!
  return
end
  DOUBLE PRECISION FUNCTION DGAMMA (X)
!
!! DGAMMA computes the complete Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      DOUBLE PRECISION (GAMMA-S, DGAMMA-D, CGAMMA-C)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DGAMMA(X) calculates the double precision complete Gamma function
! for double precision argument X.
!
! Series for GAM        on the interval  0.          to  1.00000E+00
!                                        with weighted error   5.79E-32
!                                         log weighted error  31.24
!                               significant figures required  30.00
!                                    decimal places required  32.05
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, D9LGMC, DCSEVL, DGAMLM, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920618  Removed space from variable name.  (RWC, WRB)
!***END PROLOGUE  DGAMMA
  DOUBLE PRECISION X, GAMCS(42), DXREL, PI, SINPIY, SQ2PIL, XMAX, &
    XMIN, Y, D9LGMC, DCSEVL, D1MACH
  LOGICAL FIRST
!
  SAVE GAMCS, PI, SQ2PIL, NGAM, XMIN, XMAX, DXREL, FIRST
  DATA GAMCS(  1) / +.8571195590989331421920062399942D-2      /
  DATA GAMCS(  2) / +.4415381324841006757191315771652D-2      /
  DATA GAMCS(  3) / +.5685043681599363378632664588789D-1      /
  DATA GAMCS(  4) / -.4219835396418560501012500186624D-2      /
  DATA GAMCS(  5) / +.1326808181212460220584006796352D-2      /
  DATA GAMCS(  6) / -.1893024529798880432523947023886D-3      /
  DATA GAMCS(  7) / +.3606925327441245256578082217225D-4      /
  DATA GAMCS(  8) / -.6056761904460864218485548290365D-5      /
  DATA GAMCS(  9) / +.1055829546302283344731823509093D-5      /
  DATA GAMCS( 10) / -.1811967365542384048291855891166D-6      /
  DATA GAMCS( 11) / +.3117724964715322277790254593169D-7      /
  DATA GAMCS( 12) / -.5354219639019687140874081024347D-8      /
  DATA GAMCS( 13) / +.9193275519859588946887786825940D-9      /
  DATA GAMCS( 14) / -.1577941280288339761767423273953D-9      /
  DATA GAMCS( 15) / +.2707980622934954543266540433089D-10     /
  DATA GAMCS( 16) / -.4646818653825730144081661058933D-11     /
  DATA GAMCS( 17) / +.7973350192007419656460767175359D-12     /
  DATA GAMCS( 18) / -.1368078209830916025799499172309D-12     /
  DATA GAMCS( 19) / +.2347319486563800657233471771688D-13     /
  DATA GAMCS( 20) / -.4027432614949066932766570534699D-14     /
  DATA GAMCS( 21) / +.6910051747372100912138336975257D-15     /
  DATA GAMCS( 22) / -.1185584500221992907052387126192D-15     /
  DATA GAMCS( 23) / +.2034148542496373955201026051932D-16     /
  DATA GAMCS( 24) / -.3490054341717405849274012949108D-17     /
  DATA GAMCS( 25) / +.5987993856485305567135051066026D-18     /
  DATA GAMCS( 26) / -.1027378057872228074490069778431D-18     /
  DATA GAMCS( 27) / +.1762702816060529824942759660748D-19     /
  DATA GAMCS( 28) / -.3024320653735306260958772112042D-20     /
  DATA GAMCS( 29) / +.5188914660218397839717833550506D-21     /
  DATA GAMCS( 30) / -.8902770842456576692449251601066D-22     /
  DATA GAMCS( 31) / +.1527474068493342602274596891306D-22     /
  DATA GAMCS( 32) / -.2620731256187362900257328332799D-23     /
  DATA GAMCS( 33) / +.4496464047830538670331046570666D-24     /
  DATA GAMCS( 34) / -.7714712731336877911703901525333D-25     /
  DATA GAMCS( 35) / +.1323635453126044036486572714666D-25     /
  DATA GAMCS( 36) / -.2270999412942928816702313813333D-26     /
  DATA GAMCS( 37) / +.3896418998003991449320816639999D-27     /
  DATA GAMCS( 38) / -.6685198115125953327792127999999D-28     /
  DATA GAMCS( 39) / +.1146998663140024384347613866666D-28     /
  DATA GAMCS( 40) / -.1967938586345134677295103999999D-29     /
  DATA GAMCS( 41) / +.3376448816585338090334890666666D-30     /
  DATA GAMCS( 42) / -.5793070335782135784625493333333D-31     /
  DATA PI / 3.14159265358979323846264338327950D0 /
  DATA SQ2PIL / 0.91893853320467274178032973640562D0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DGAMMA
  if (FIRST) THEN
     NGAM = INITDS (GAMCS, 42, 0.1*REAL(D1MACH(3)) )
!
     call DGAMLM (XMIN, XMAX)
     DXREL = SQRT(D1MACH(4))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 10.D0) go to 50
!
! COMPUTE GAMMA(X) FOR -XBND  <=  X  <=  XBND.  REDUCE INTERVAL AND FIND
! GAMMA(1+Y) FOR 0.0  <=  Y  <  1.0 FIRST OF ALL.
!
  N = X
  if (X < 0.D0) N = N - 1
  Y = X - N
  N = N - 1
  DGAMMA = 0.9375D0 + DCSEVL (2.D0*Y-1.D0, GAMCS, NGAM)
  if (N == 0) RETURN
!
  if (N > 0) go to 30
!
! COMPUTE GAMMA(X) FOR X  <  1.0
!
  N = -N
  if (X  ==  0.D0) call XERMSG ('SLATEC', 'DGAMMA', 'X IS 0', 4, 2)
  if (X  <  0.0 .AND. X+N-2  ==  0.D0) call XERMSG ('SLATEC', &
     'DGAMMA', 'X IS A NEGATIVE INTEGER', 4, 2)
  if (X  <  (-0.5D0) .AND. ABS((X-AINT(X-0.5D0))/X)  <  DXREL) &
     call XERMSG ('SLATEC', 'DGAMMA', &
     'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', &
     1, 1)
!
  DO 20 I=1,N
    DGAMMA = DGAMMA/(X+I-1 )
 20   CONTINUE
  return
!
! GAMMA(X) FOR X  >=  2.0 AND X  <=  10.0
!
 30   DO 40 I=1,N
    DGAMMA = (Y+I) * DGAMMA
 40   CONTINUE
  return
!
! GAMMA(X) FOR ABS(X)  >  10.0.  RECALL Y = ABS(X).
!
 50   if (X  >  XMAX) call XERMSG ('SLATEC', 'DGAMMA', &
     'X SO BIG GAMMA OVERFLOWS', 3, 2)
!
  DGAMMA = 0.D0
  if (X  <  XMIN) call XERMSG ('SLATEC', 'DGAMMA', &
     'X SO SMALL GAMMA UNDERFLOWS', 2, 1)
  if (X < XMIN) RETURN
!
  DGAMMA = EXP ((Y-0.5D0)*LOG(Y) - Y + SQ2PIL + D9LGMC(Y) )
  if (X > 0.D0) RETURN
!
  if (ABS((X-AINT(X-0.5D0))/X)  <  DXREL) call XERMSG ('SLATEC', &
     'DGAMMA', &
     'ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER', 1, 1)
!
  SINPIY = SIN (PI*Y)
  if (SINPIY  ==  0.D0) call XERMSG ('SLATEC', 'DGAMMA', &
     'X IS A NEGATIVE INTEGER', 4, 2)
!
  DGAMMA = -PI/(Y*SINPIY*DGAMMA)
!
  return
end
  DOUBLE PRECISION FUNCTION D9LGMC (X)
!
!! D9LGMC computes the log Gamma correction factor so that ...
!  LOG(DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-5.)*LOG(X) - X + D9LGMC(X).
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-C)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
!             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the log gamma correction factor for X  >=  10. so that
! LOG (DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + D9lGMC(X)
!
! Series for ALGM       on the interval  0.          to  1.00000E-02
!                                        with weighted error   1.28E-31
!                                         log weighted error  30.89
!                               significant figures required  29.81
!                                    decimal places required  31.48
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  D9LGMC
  DOUBLE PRECISION X, ALGMCS(15), XBIG, XMAX, DCSEVL, D1MACH
  LOGICAL FIRST
  SAVE ALGMCS, NALGM, XBIG, XMAX, FIRST
  DATA ALGMCS(  1) / +.1666389480451863247205729650822D+0      /
  DATA ALGMCS(  2) / -.1384948176067563840732986059135D-4      /
  DATA ALGMCS(  3) / +.9810825646924729426157171547487D-8      /
  DATA ALGMCS(  4) / -.1809129475572494194263306266719D-10     /
  DATA ALGMCS(  5) / +.6221098041892605227126015543416D-13     /
  DATA ALGMCS(  6) / -.3399615005417721944303330599666D-15     /
  DATA ALGMCS(  7) / +.2683181998482698748957538846666D-17     /
  DATA ALGMCS(  8) / -.2868042435334643284144622399999D-19     /
  DATA ALGMCS(  9) / +.3962837061046434803679306666666D-21     /
  DATA ALGMCS( 10) / -.6831888753985766870111999999999D-23     /
  DATA ALGMCS( 11) / +.1429227355942498147573333333333D-24     /
  DATA ALGMCS( 12) / -.3547598158101070547199999999999D-26     /
  DATA ALGMCS( 13) / +.1025680058010470912000000000000D-27     /
  DATA ALGMCS( 14) / -.3401102254316748799999999999999D-29     /
  DATA ALGMCS( 15) / +.1276642195630062933333333333333D-30     /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  D9LGMC
  if (FIRST) THEN
     NALGM = INITDS (ALGMCS, 15, REAL(D1MACH(3)) )
     XBIG = 1.0D0/SQRT(D1MACH(3))
     XMAX = EXP (MIN(LOG(D1MACH(2)/12.D0), -LOG(12.D0*D1MACH(1))))
  end if
  FIRST = .FALSE.
!
  if (X  <  10.D0) call XERMSG ('SLATEC', 'D9LGMC', &
     'X MUST BE GE 10', 1, 2)
  if (X >= XMAX) go to 20
!
  D9LGMC = 1.D0/(12.D0*X)
  if (X < XBIG) D9LGMC = DCSEVL (2.0D0*(10.D0/X)**2-1.D0, ALGMCS, &
    NALGM) / X
  return
!
 20   D9LGMC = 0.D0
  call XERMSG ('SLATEC', 'D9LGMC', 'X SO BIG D9LGMC UNDERFLOWS', 2, &
     1)
  return
!
end
subroutine DGAMLM (XMIN, XMAX)
!
!! DGAMLM computes bounds for the argument in the Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A, R2
!***TYPE      DOUBLE PRECISION (GAMLIM-S, DGAMLM-D)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Calculate the minimum and maximum legal bounds for X in gamma(X).
! XMIN and XMAX are not the only bounds, but they are the only non-
! trivial ones to calculate.
!
!             Output Arguments --
! XMIN   double precision minimum legal value of X in gamma(X).  Any
!        smaller value of X might result in underflow.
! XMAX   double precision maximum legal value of X in gamma(X).  Any
!        larger value of X might cause overflow.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DGAMLM
  DOUBLE PRECISION XMIN, XMAX, ALNBIG, ALNSML, XLN, XOLD, D1MACH
!***FIRST EXECUTABLE STATEMENT  DGAMLM
  ALNSML = LOG(D1MACH(1))
  XMIN = -ALNSML
  DO 10 I=1,10
    XOLD = XMIN
    XLN = LOG(XMIN)
    XMIN = XMIN - XMIN*((XMIN+0.5D0)*XLN - XMIN - 0.2258D0 + ALNSML) &
      / (XMIN*XLN+0.5D0)
    if (ABS(XMIN-XOLD) < 0.005D0) go to 20
 10   CONTINUE
  call XERMSG ('SLATEC', 'DGAMLM', 'UNABLE TO FIND XMIN', 1, 2)
!
 20   XMIN = -XMIN + 0.01D0
!
  ALNBIG = LOG (D1MACH(2))
  XMAX = ALNBIG
  DO 30 I=1,10
    XOLD = XMAX
    XLN = LOG(XMAX)
    XMAX = XMAX - XMAX*((XMAX-0.5D0)*XLN - XMAX + 0.9189D0 - ALNBIG) &
      / (XMAX*XLN-0.5D0)
    if (ABS(XMAX-XOLD) < 0.005D0) go to 40
 30   CONTINUE
  call XERMSG ('SLATEC', 'DGAMLM', 'UNABLE TO FIND XMAX', 2, 2)
!
 40   XMAX = XMAX - 0.01D0
  XMIN = MAX (XMIN, -XMAX+1.D0)
!
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
!  Don''t print a blank line.
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
function INITDS (OS, NOS, ETA)
!
!! INITDS determines the number of terms needed in an orthogonal ...
!            polynomial series so that it meets a specified accuracy.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C3A2
!***TYPE      DOUBLE PRECISION (INITS-S, INITDS-D)
!***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
!             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!  Initialize the orthogonal series, represented by the array OS, so
!  that INITDS is the number of terms needed to insure the error is no
!  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
!  machine precision.
!
!             Input Arguments --
!   OS     double precision array of NOS coefficients in an orthogonal
!          series.
!   NOS    number of coefficients in OS.
!   ETA    single precision scalar containing requested accuracy of
!          series.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891115  Modified error message.  (WRB)
!   891115  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  INITDS
  DOUBLE PRECISION OS(*)
!***FIRST EXECUTABLE STATEMENT  INITDS
  if (NOS  <  1) call XERMSG ('SLATEC', 'INITDS', &
     'Number of coefficients is less than 1', 2, 1)
!
  ERR = 0.
  DO 10 II = 1,NOS
    I = NOS + 1 - II
    ERR = ERR + ABS(REAL(OS(I)))
    if (ERR > ETA) go to 20
   10 CONTINUE
!
   20 if (I  ==  NOS) call XERMSG ('SLATEC', 'INITDS', &
     'Chebyshev series too short for specified accuracy', 1, 1)
  INITDS = I
!
  return
end
