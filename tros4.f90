!****************************************************************
!* NUMERICAL SOLUTION OF A STIFF SYSTEM OF FIRST 0RDER ORDINARY *
!* DIFFERENTIAL EQUATIONS Y'=F(X,Y) BY ROSENBROCK METHOD.       *
!* ------------------------------------------------------------ *
!* SAMPLE RUN:                                                  *
!* Example #1:                                                  *
!* (Solve set of differential equations (N=2):                  *
!*  F(1) = Y(1) * Y(2) + COS(X) - HALF * SIN(TWO * X)           *
!*  F(2) = Y(1) * Y(1) + Y(2) * Y(2) - (ONE + SIN(X))           *
!*  Find values of F(1), F(2) at X=1.5).                        *
!*                                                              *
!*  SOLUTION AT X=    1.50000000000000                          *
!*  Y(1) =  0.12359935E+01                                      *
!*  Y(2) = -0.10494372E+00                                      *
!*                                                              *
!*  LAST STEP SIZE =  4.150113101356574E-002                    *
!*  ERROR CODE =           1                                    *
!*                                                              *
!* Example #2:                                                  *
!* (Solve set of differential equations (N=5):                  *
!*  F(1) = Y(2)                                                 *
!*  F(2) = Y(3)                                                 *
!*  F(3) = Y(4)                                                 *
!*  F(4) = Y(5)                                                 *
!*  F(5) = (45.d0 * Y(3) * Y(4) * Y(5) -                        *
!*          40.d0 * Y(4) * Y(4) * Y(4)) / (NINE * Y(3) * Y(3))  *
!*  Find values of F(1), F(2), ..., F(5) at X=1.5).             *
!*                                                              *
!*  SOLUTION AT X=    1.50000000000000                          *
!*  Y(1) =  0.43639610E+01                                      *
!*  Y(2) =  0.40000000E+01                                      *
!*  Y(3) =  0.28284271E+01                                      *
!*  Y(4) =  0.14790900E-10                                      *
!*  Y(5) = -0.37712362E+01                                      * 
!*                                                              *
!*  LAST STEP SIZE =  3.825256949194526E-003                    *
!*  ERROR CODE =           1                                    *
!* ------------------------------------------------------------ *
!* Ref.: From Numath Library By Tuan Dang Trong in Fortran 77   *
!*       [BIBLI 18].                                            * 
!*                                                              *
!*                       F90 Release 1.0 By J-P Moreau, Paris   *
!*                                (www.jpmoreau.fr)             *
!****************************************************************
! LIST OF USED SUBROUTINES (HERE INCLUDED):
! ========================================
! ROS4, RO4COR, SHAMP, GRK4A, GRK4T, VELDD, VELDS, LSTAB,
! DEC, DECB, SOL, SOLB.
! (Other used subroutines/functions are standard Fortran Library).
! ---------------------------------------------------------------
!PROGRAM TROS4

!REAL*8 X,XEND, H,RTOL,ATOL
!REAL*8,POINTER :: Y(:)
!REAL*8,POINTER :: WORK(:)
!INTEGER,POINTER :: IWORK(:)

!Initialize parameters (see ROS4) 
!N=2            !DIMENSION OF THE SYSTEM (N=5 for #2)
!IFCN=1         !FCN(N,X,Y,F) MAY DEPEND ON X
!X=0.d0         !INITIAL X-VALUE
!XEND=1.5d0     !FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
!H=0.001d0      !INITIAL STEP SIZE GUESS
!RTOL=1.d-6     !RELATIVE ERROR TOLERANCE (HERE SCALAR)
!ATOL=1.d-8     !ABSOLUTE ERROR TOLERANCE (HERE SCALAR)
               !1.d-10 for both in #2.
!ITOL=0         !BOTH RTOL AND ATOL ARE SCALARS
!IJAC=0         !JACOBIAN IS COMPUTED INTERNALLY BY FINITE
               !DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED
!MLJAC=N        !JACOBIAN IS A FULL MATRIX. THE LINEAR ALGEBRA
               !IS DONE BY FULL-MATRIX GAUSS-ELIMINATION
!IDFX=0         !DF/DX IS COMPUTED INTERNALLY BY FINITE
               !DIFFERENCES, SUBROUTINE "DFX" IS NEVER CALLED
!IMAS=0         !M IS SUPPOSED TO BE THE IDENTITY
               !MATRIX, MAS IS NEVER CALLED
!MLMAS=N        !MLMAS=N: THE FULL MATRIX CASE. THE LINEAR ALGEBRA
               !IS DONE BY FULL-MATRIX GAUSS-ELIMINATION
!IOUT=0         !SUBROUTINE SOLOUT IS NEVER CALLED
!LE1=N          !IF MLJAC=N (FULL JACOBIAN)
!LJAC=N         !IF MLJAC=N (FULL JACOBIAN)
!LMAS=0         !IF IMAS=0

!LIWORK= N+2                     !DECLARED LENGTH OF ARRAY "IWORK"
!LWORK = N*(LJAC+LMAS+LE1+8)+5   !DECLARED LENGTH OF ARRAY "LWORK"

!dynamic allocations 
!allocate(Y(1:N),stat=ialloc)
!allocate(WORK(1:LWORK),stat=ialloc)
!allocate(IWORK(1:LIWORK),stat=ialloc)

!WORK=0.d0      !This triggers default values (see ROS4)
!IWORK=0

!Y(1)=0.5d0     !INITIAL VALUES FOR Y 
!Y(2)=0.5d0     !In #2, Y(1) = Y(2) = ... = Y(5) = 1.d0

!call Rosenbrock subroutine with appropriate parameters
!(here, FCN has been removed from input parameters).
!CALL ROS4(N,IFCN,X,Y,XEND,H,               &
!          RTOL,ATOL,ITOL,                  &
!          JAC ,IJAC,MLJAC,MUJAC,DFX,IDFX,  &
!          MAS ,IMAS,MLMAS,MUMAS,           &
!          SOLOUT,IOUT,                     &
!          WORK,LWORK,IWORK,LIWORK,IDID)

!print results
!print *,' '
!print *,' SOLUTION AT X=', X
!do I=1,N
!  write(*,10) I, Y(I)
!end do
!print *,' '
!print *,' LAST STEP SIZE =', H
!print *,' ERROR CODE =', IDID
!print *,' '

!10 format('  Y(',I1,') = ',E15.8)

!END  !of main program


!define example #1
!SUBROUTINE FCN(N,X,Y,F)
!parameter(HALF=0.5d0,ONE=1.d0,TWO=2.d0)
!REAL*8 X,Y(N),F(N)
!  F(1) = Y(1) * Y(2) + DCOS(X) - HALF * DSIN(TWO * X)
!  F(2) = Y(1) * Y(1) + Y(2) * Y(2) - (ONE + DSIN(X))
!RETURN
!END

!define example #2
!SUBROUTINE FCN(N,X,Y,F)
!parameter(NINE=9.d0)
!REAL*8 X,Y(N),F(N)
!  F(1) = Y(2)
!  F(2) = Y(3)
!  F(3) = Y(4)
!  F(4) = Y(5)
!  F(5) = (45.d0 * Y(3) * Y(4) * Y(5) -  &
!          40.d0 * Y(4) * Y(4) * Y(4)) / (NINE * Y(3) * Y(3)) 
!RETURN
!END

!**********************************************************************
SUBROUTINE ROS4(N,IFCN,X,Y,XEND,H,               &
                RTOL,ATOL,ITOL,                  &
                JAC ,IJAC,MLJAC,MUJAC,DFX,IDFX,  &
                MAS ,IMAS,MLMAS,MUMAS,           &
                SOLOUT,IOUT,                     &
                WORK,LWORK,IWORK,LIWORK,IDID)
! ---------------------------------------------------------------------
!     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
!     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS  MY'=F(X,Y).
!     THIS IS AN EMBEDDED ROSENBROCK METHOD OF ORDER (3)4
!     (WITH STEP SIZE CONTROL).
!     C.F. SECTION IV.7
!
!     AUTHORS: E. HAIRER AND G. WANNER
!              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
!              CH-1211 GENEVE 24, SWITZERLAND
!              E-MAIL:  HAIRER@CGEUGE51.BITNET,  WANNER@CGEUGE51.BITNET
!
!     THIS CODE IS PART OF THE BOOK:
!         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!         SPRINGER-VERLAG (1990)
!
!     VERSION OF OCTOBER 12, 1990
!
!     INPUT PARAMETERS
!     ----------------
!     N           DIMENSION OF THE SYSTEM
!
!     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
!                 VALUE OF F(X,Y):
!                    SUBROUTINE FCN(N,X,Y,F)
!                    REAL*8 X,Y(N),F(N)
!                    F(1)=...   ETC.
!
!     IFCN        GIVES INFORMATION ON FCN:
!                    IFCN=0: F(X,Y) INDEPENDENT OF X (AUTONOMOUS)
!                    IFCN=1: F(X,Y) MAY DEPEND ON X (NON-AUTONOMOUS)
!
!     X           INITIAL X-VALUE
!
!     Y(N)        INITIAL VALUES FOR Y
!
!     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
!
!     H           INITIAL STEP SIZE GUESS;
!                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,
!                 H=1.D0/(NORM OF F'), USUALLY 1.D-2 OR 1.D-3, IS GOOD.
!                 THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY
!                 ADAPTS ITS STEP SIZE. STUDY THE CHOSEN VALUES FOR A FEW
!                 STEPS IN SUBROUTINE "SOLOUT", WHEN YOU ARE NOT SURE.
!                 (IF H=0.D0, THE CODE PUTS H=1.D-6).
!
!     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
!                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
!
!     ITOL        SWITCH FOR RTOL AND ATOL:
!                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
!                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
!                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
!                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
!                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
!                     RTOL(I)*ABS(Y(I))+ATOL(I).
!
!     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
!                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y
!                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
!                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
!                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM:
!                    SUBROUTINE JAC(N,X,Y,DFY,LDFY)
!                    REAL*8 X,Y(N),DFY(LDFY,N)
!                    DFY(1,1)= ...
!                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS
!                 FURNISHED BY THE CALLING PROGRAM.
!                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO
!                    BE FULL AND THE PARTIAL DERIVATIVES ARE
!                    STORED IN DFY AS
!                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
!                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND
!                    THE PARTIAL DERIVATIVES ARE STORED
!                    DIAGONAL-WISE AS
!                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J).
!
!     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
!                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
!                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED.
!                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC.
!
!     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN:
!                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR
!                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
!                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN
!                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
!                       THE MAIN DIAGONAL).
!
!     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON-
!                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
!                 NEED NOT BE DEFINED IF MLJAC=N.
!
!     DFX         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
!                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO X
!                 (THIS ROUTINE IS ONLY CALLED IF IDFX=1 AND IFCN=1;
!                 SUPPLY A DUMMY SUBROUTINE IN THE CASE IDFX=0 OR IFCN=0).
!                 OTHERWISE, THIS SUBROUTINE MUST HAVE THE FORM
!                    SUBROUTINE DFX(N,X,Y,FX)
!                    REAL*8 X,Y(N),FX(N)
!                    FX(1)= ...
!
!     IDFX        SWITCH FOR THE COMPUTATION OF THE DF/DX:
!                    IDFX=0: DF/DX IS COMPUTED INTERNALLY BY FINITE
!                       DIFFERENCES, SUBROUTINE "DFX" IS NEVER CALLED.
!                    IDFX=1: DF/DX IS SUPPLIED BY SUBROUTINE DFX.
!
!     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      -----
!     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): -
!
!     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
!                 MATRIX M.
!                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY
!                 MATRIX AND NEEDS NOT TO BE DEFINED;
!                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
!                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM
!                    SUBROUTINE MAS(N,AM,LMAS)
!                    REAL*8 AM(LMAS,N)
!                    AM(1,1)= ....
!                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED
!                    AS FULL MATRIX LIKE
!                         AM(I,J) = M(I,J)
!                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED
!                    DIAGONAL-WISE AS
!                         AM(I-J+MUMAS+1,J) = M(I,J).
!
!     IMAS       GIVES INFORMATION ON THE MASS-MATRIX:
!                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY
!                       MATRIX, MAS IS NEVER CALLED.
!                    IMAS=1: MASS-MATRIX  IS SUPPLIED.
!
!     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX:
!                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR
!                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
!                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE
!                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
!                       THE MAIN DIAGONAL).
!                 MLMAS IS SUPPOSED TO BE <= MLJAC.
!
!     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON-
!                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
!                 NEED NOT BE DEFINED IF MLMAS=N.
!                 MUMAS IS SUPPOSED TO BE <= MUJAC.
!
!     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
!                 NUMERICAL SOLUTION DURING INTEGRATION.
!                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
!                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0.
!                 IT MUST HAVE THE FORM
!                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,IRTRN)
!                    REAL*8 X,Y(N)
!                    ....
!                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
!                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
!                    THE FIRST GRID-POINT).
!                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
!                    IS SET <0, ROS4 RETURNS TO THE CALLING PROGRAM.
!
!     IOUT        GIVES INFORMATION ON THE SUBROUTINE SOLOUT:
!                    IOUT=0: SUBROUTINE IS NEVER CALLED
!                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT
!
!     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
!                 SERVES AS WORKING SPACE FOR ALL VECTORS AND MATRICES.
!                 "LWORK" MUST BE AT LEAST
!                             N*(LJAC+LMAS+LE1+8)+5
!                 WHERE
!                    LJAC=N              IF MLJAC=N (FULL JACOBIAN)
!                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JAC.)
!                 AND
!                    LMAS=0              IF IMAS=0
!                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL)
!                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M.)
!                 AND
!                    LE1=N               IF MLJAC=N (FULL JACOBIAN)
!                    LE1=2*MLJAC+MUJAC+1 IF MLJAC<N (BANDED JAC.).
!
!                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE
!                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM
!                 STORAGE REQUIREMENT IS
!                             LWORK = 2*N*N+8*N+5.
!
!     LWORK       DECLARED LENGHT OF ARRAY "WORK".
!
!     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
!                 "LIWORK" MUST BE AT LEAST N+2.
!
!     LIWORK      DECLARED LENGHT OF ARRAY "IWORK".
!
! ----------------------------------------------------------------------
!
!     SOPHISTICATED SETTING OF PARAMETERS
!     -----------------------------------
!              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK
!              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),..,WORK(5)
!              AS WELL AS IWORK(1),IWORK(2) DIFFERENT FROM ZERO.
!              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
!
!    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
!              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 100000.
!
!    IWORK(2)  SWITCH FOR THE CHOICE OF THE COEFFICIENTS
!              IF IWORK(2).EQ.1  METHOD OF SHAMPINE
!              IF IWORK(2).EQ.2  METHOD GRK4T OF KAPS-RENTROP
!              IF IWORK(2).EQ.3  METHOD GRK4A OF KAPS-RENTROP
!              IF IWORK(2).EQ.4  METHOD OF VAN VELDHUIZEN (GAMMA=1/2)
!              IF IWORK(2).EQ.5  METHOD OF VAN VELDHUIZEN ("D-STABLE")
!              IF IWORK(2).EQ.6  AN L-STABLE METHOD
!              THE DEFAULT VALUE (FOR IWORK(2)=0) IS IWORK(2)=2.
!
!    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
!
!    WORK(2)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
!
!    WORK(3), WORK(4)   PARAMETERS FOR STEP SIZE SELECTION
!              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
!                 WORK(3) <= HNEW/HOLD <= WORK(4)
!              DEFAULT VALUES: WORK(3)=0.2D0, WORK(4)=6.D0
!
!    WORK(5)   AVOID THE HUMP: AFTER TWO CONSECUTIVE STEP REJECTIONS
!              THE STEP SIZE IS MULTIPLIED BY WORK(5)
!              DEFAULT VALUES: WORK(5)=0.1D0
!
!-----------------------------------------------------------------------
!
!     OUTPUT PARAMETERS
!     -----------------
!     X           X-VALUE WHERE THE SOLUTION IS COMPUTED
!                 (AFTER SUCCESSFUL RETURN X=XEND)
!
!     Y(N)        SOLUTION AT X
!
!     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
!
!     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
!                   IDID=1  COMPUTATION SUCCESSFUL,
!                   IDID=-1 COMPUTATION UNSUCCESSFUL.
!
! ---------------------------------------------------------
! *** *** *** *** *** *** *** *** *** *** *** *** ***
!          DECLARATIONS
! *** *** *** *** *** *** *** *** *** *** *** *** ***
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(N),ATOL(1),RTOL(1),WORK(LWORK),IWORK(LIWORK)
      LOGICAL AUTNMS,IMPLCT,JBAND,ARRET
      EXTERNAL FCN,JAC,DFX,MAS,SOLOUT
      COMMON/STAT/NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL
! --------------------------------------------------------------------
! --- COMMON STAT CAN BE USED FOR STATISTICS
! ---    NFCN      NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
!                  EVALUATION OF THE JACOBIAN ARE NOT COUNTED)
! ---    NJAC      NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
!                  OR NUMERICALLY)
! ---    NSTEP     NUMBER OF COMPUTED STEPS
! ---    NACCPT    NUMBER OF ACCEPTED STEPS
! ---    NREJCT    NUMBER OF REJECTED STEPS (AFTER AT LEAST ONE STEP
!                  HAS BEEN ACCEPTED)
! ---    NDEC      NUMBER OF LU-DECOMPOSITIONS (N-DIMENSIONAL MATRIX)
! ---    NSOL      NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS
! --------------------------------------------------------------------
! *** *** *** *** *** *** ***
!    SETTING THE PARAMETERS
! *** *** *** *** *** *** ***
      NFCN=0
      NJAC=0
      NSTEP=0
      NACCPT=0
      NREJCT=0
      NDEC=0
      NSOL=0
      ARRET=.FALSE.
! -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
      IF(IWORK(1).EQ.0)THEN
         NMAX=100000
      ELSE
         NMAX=IWORK(1)
         IF(NMAX.LE.0)THEN
            WRITE(6,*)' WRONG INPUT IWORK(1)=',IWORK(1)
            ARRET=.TRUE.
         END IF
      END IF
! -------- METH   COEFFICIENTS OF THE METHOD
      IF(IWORK(2).EQ.0)THEN
         METH=2
      ELSE
         METH=IWORK(2)
         IF(METH.LE.0.OR.METH.GE.7)THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(2)=',IWORK(2)
            ARRET=.TRUE.
         END IF
      END IF
! -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0
      IF(WORK(1).EQ.0.D0)THEN
         UROUND=1.D-16
      ELSE
         UROUND=WORK(1)
         IF(UROUND.LE.1.D-14.OR.UROUND.GE.1.D0)THEN
            WRITE(6,*)' COEFFICIENTS HAVE 16 DIGITS, UROUND=',WORK(1)
            ARRET=.TRUE.
         END IF
      END IF
! -------- MAXIMAL STEP SIZE
      IF(WORK(2).EQ.0.D0)THEN
         HMAX=XEND-X
      ELSE
         HMAX=WORK(2)
      END IF
! -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION
      IF(WORK(3).EQ.0.D0)THEN
         FAC1=5.D0
      ELSE
         FAC1=1.D0/WORK(3)
      END IF
      IF(WORK(4).EQ.0.D0)THEN
         FAC2=1.D0/6.0D0
      ELSE
         FAC2=1.D0/WORK(4)
      END IF
! -------  FACREJ    FOR THE HUMP
      IF(WORK(5).EQ.0.D0)THEN
         FACREJ=0.1D0
      ELSE
         FACREJ=WORK(5)
      END IF
! --------- CHECK IF TOLERANCES ARE O.K.
      IF (ITOL.EQ.0) THEN
          IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES ARE TOO SMALL'
              ARRET=.TRUE.
          END IF
      ELSE
          DO 15 I=1,N
          IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES(',I,') ARE TOO SMALL'
              ARRET=.TRUE.
          END IF
  15      CONTINUE
      END IF
! *** *** *** *** *** *** *** *** *** *** *** *** ***
!         COMPUTATION OF ARRAY ENTRIES
! *** *** *** *** *** *** *** *** *** *** *** *** ***
! ---- AUTONOMOUS, IMPLICIT, BANDED OR NOT ?
      AUTNMS=IFCN.EQ.0
      IMPLCT=IMAS.NE.0
      JBAND=MLJAC.NE.N
      ARRET=.FALSE.
! -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS ---
! -- JACOBIAN
      IF(JBAND)THEN
         LDJAC=MLJAC+MUJAC+1
      ELSE
         LDJAC=N
      END IF
! -- MATRIX E FOR LINEAR ALGEBRA
      IF(JBAND)THEN
         LDE=2*MLJAC+MUJAC+1
      ELSE
         LDE=N
      END IF
! -- MASS MATRIX
      IF (IMPLCT) THEN
          IF (MLMAS.NE.N) THEN
              LDMAS=MLMAS+MUMAS+1
          ELSE
              LDMAS=N
          END IF
! ------ BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF "JAC"
          IF (MLMAS.GT.MLJAC.OR.MUMAS.GT.MUJAC) THEN
            WRITE (6,*) 'BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF "JAC"'
            ARRET=.TRUE.
          END IF
      ELSE
          LDMAS=0
      END IF
      LDMAS2=MAX(1,LDMAS)
! ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK -----
      IEYNEW=6
      IEDY1=IEYNEW+N
      IEDY=IEDY1+N
      IEAK1=IEDY+N
      IEAK2=IEAK1+N
      IEAK3=IEAK2+N
      IEAK4=IEAK3+N
      IEFX =IEAK4+N
      IEJAC=IEFX +N
      IEMAS=IEJAC+N*LDJAC
      IEE  =IEMAS+N*LDMAS
! ------ TOTAL STORAGE REQUIREMENT -----------
      ISTORE=IEE+N*LDE-1
      IF(ISTORE.GT.LWORK)THEN
         WRITE(6,*)' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE
         ARRET=.TRUE.
      END IF
! ------- ENTRY POINTS FOR INTEGER WORKSPACE -----
      IEIP=3
! --------- TOTAL REQUIREMENT ---------------
      ISTORE=IEIP+N-1
      IF(ISTORE.GT.LIWORK)THEN
         WRITE(6,*)' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',ISTORE
         ARRET=.TRUE.
      END IF
! ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
      IF (ARRET) THEN
         IDID=-1
         RETURN
      END IF
! -------- CALL TO CORE INTEGRATOR ------------
      CALL RO4COR(N,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,JAC,IJAC,        &
         MLJAC,MUJAC,DFX,IDFX,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,     &
         NMAX,UROUND,METH,FAC1,FAC2,FACREJ,AUTNMS,IMPLCT,JBAND,     &
         LDJAC,LDE,LDMAS2,WORK(IEYNEW),WORK(IEDY1),WORK(IEDY),      &
         WORK(IEAK1),WORK(IEAK2),WORK(IEAK3),WORK(IEAK4),           &
         WORK(IEFX),WORK(IEJAC),WORK(IEE),WORK(IEMAS),IWORK(IEIP))
! ----------- RETURN -----------
      RETURN
      END
!
!
!  --------- ... AND HERE IS THE CORE INTEGRATOR  ----------
!
      SUBROUTINE RO4COR(N,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,JAC,         &
        IJAC,MLJAC,MUJAC,DFX,IDFX,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,   &
        NMAX,UROUND,METH,FAC1,FAC2,FACREJ,AUTNMS,IMPLCT,BANDED,       &
        LFJAC,LE,LDMAS,YNEW,DY1,DY,AK1,AK2,AK3,AK4,FX,FJAC,E,FMAS,IP)
! ----------------------------------------------------------
!     CORE INTEGRATOR FOR ROS4
!     PARAMETERS SAME AS IN ROS4 WITH WORKSPACE ADDED
! ----------------------------------------------------------
!         DECLARATIONS
! ----------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(N),YNEW(N),DY1(N),DY(N),AK1(N), &
	    AK2(N),AK3(N),AK4(N),FX(N),            &
        FJAC(LFJAC,N),E(LE,N),FMAS(LDMAS,N),ATOL(1),RTOL(1)
      INTEGER IP(N)
      LOGICAL REJECT,RJECT2,AUTNMS,IMPLCT,BANDED
      COMMON/STAT/NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL

! ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ----------
      IF (IMPLCT) CALL MAS(N,FMAS,LDMAS)
! ---- PREPARE BANDWIDTHS -----
      IF (BANDED) THEN
          MLE=MLJAC
          MUE=MUJAC
          MBJAC=MLJAC+MUJAC+1
          MBB=MLMAS+MUMAS+1
          MDIAG=MLE+MUE+1
          MBDIAG=MUMAS+1
          MDIFF=MLE+MUE-MUMAS
      END IF
! *** *** *** *** *** *** ***
!  INITIALISATIONS
! *** *** *** *** *** *** ***
      POSNEG=DSIGN(1.D0,XEND-X)
      IF (METH.EQ.1) CALL SHAMP (A21,A31,A32,C21,C31,C32,C41,C42,C43,  &
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IF (METH.EQ.2) CALL GRK4T (A21,A31,A32,C21,C31,C32,C41,C42,C43,  &
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IF (METH.EQ.3) CALL GRK4A (A21,A31,A32,C21,C31,C32,C41,C42,C43,  &
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IF (METH.EQ.4) CALL VELDS (A21,A31,A32,C21,C31,C32,C41,C42,C43,  &
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IF (METH.EQ.5) CALL VELDD (A21,A31,A32,C21,C31,C32,C41,C42,C43,  &
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IF (METH.EQ.6) CALL LSTAB (A21,A31,A32,C21,C31,C32,C41,C42,C43,  &
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)

! --- INITIAL PREPARATIONS
      HMAXN=DMIN1(DABS(HMAX),DABS(XEND-X))
      H=DMIN1(DMAX1(1.D-10,DABS(H)),HMAXN)
      H=DSIGN(H,POSNEG)
      REJECT=.FALSE.
      NSING=0
      IRTRN=1
      XXOLD=X

      IF (IOUT.NE.0) CALL SOLOUT(NACCPT+1,XXOLD,X,Y,N,IRTRN)
      IF (IRTRN.LT.0) GOTO 79
! --- BASIC INTEGRATION STEP
   1  IF (NSTEP.GT.NMAX.OR.X+.1D0*H.EQ.X.OR.ABS(H).LE.UROUND) GOTO 79
      IF ((X-XEND)*POSNEG+UROUND.GT.0.D0) THEN
          H=HOPT
          IDID=1
          RETURN
      END IF
      HOPT=H
      IF ((X+H-XEND)*POSNEG.GT.0.D0) H=XEND-X

      CALL FCN(N,X,Y,DY1)

      NFCN=NFCN+1
! *** *** *** *** *** *** ***
!  COMPUTATION OF THE JACOBIAN
! *** *** *** *** *** *** ***
      NJAC=NJAC+1
      IF (IJAC.EQ.0) THEN
! --- COMPUTE JACOBIAN MATRIX NUMERICALLY
          IF (BANDED) THEN
! --- JACOBIAN IS BANDED
              MUJACP=MUJAC+1
              MD=MIN(MBJAC,N)
              DO 16 K=1,MD
              J=K
 12           AK2(J)=Y(J)
              AK3(J)=DSQRT(UROUND*MAX(1.D-5,ABS(Y(J))))
              Y(J)=Y(J)+AK3(J)
              J=J+MD
              IF (J.LE.N) GOTO 12
              CALL FCN(N,X,Y,AK1)
              J=K
              LBEG=MAX(1,J-MUJAC)
 14           LEND=MIN(N,J+MLJAC)
              Y(J)=AK2(J)
              MUJACJ=MUJACP-J
              DO 15 L=LBEG,LEND
 15           FJAC(L+MUJACJ,J)=(AK1(L)-DY1(L))/AK3(J)
              J=J+MD
              LBEG=LEND+1
              IF (J.LE.N) GOTO 14
 16           CONTINUE
          ELSE
! --- JACOBIAN IS FULL
              DO 18 I=1,N
              YSAFE=Y(I)
              DELT=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE)))
              Y(I)=YSAFE+DELT
              CALL FCN(N,X,Y,AK1)
              DO 17 J=1,N
  17          FJAC(J,I)=(AK1(J)-DY1(J))/DELT
  18          Y(I)=YSAFE
              MLJAC=N
          END IF
      ELSE
! --- COMPUTE JACOBIAN MATRIX ANALYTICALLY
          CALL JAC(N,X,Y,FJAC,LFJAC)
      END IF
      IF (.NOT.AUTNMS) THEN
          IF (IDFX.EQ.0) THEN
! --- COMPUTE NUMERICALLY THE DERIVATIVE WITH RESPECT TO X
              DELT=DSQRT(UROUND*MAX(1.D-5,ABS(X)))
              XDELT=X+DELT
              CALL FCN(N,XDELT,Y,AK1)
              DO 19 J=1,N
  19          FX(J)=(AK1(J)-DY1(J))/DELT
          ELSE
! --- COMPUTE ANALYTICALLY THE DERIVATIVE WITH RESPECT TO X
              CALL DFX(N,X,Y,FX)
          END IF
      END IF
   2  CONTINUE
! *** *** *** *** *** *** ***
!  COMPUTE THE STAGES
! *** *** *** *** *** *** ***
      NDEC=NDEC+1
      HC21=C21/H
      HC31=C31/H
      HC32=C32/H
      HC41=C41/H
      HC42=C42/H
      HC43=C43/H
      FAC=1.D0/(H*GAMMA)
      IF (IMPLCT) THEN
          IF (BANDED) THEN
! --- THE MATRIX E (B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX)
              DO 101 J=1,N
              I1=MAX0(1,MUJAC+2-J)
              I2=MIN0(MBJAC,N+MUJAC+1-J)
              DO 101 I=I1,I2
  101         E(I+MLE,J)=-FJAC(I,J)
              DO 102 J=1,N
              I1=MAX0(1,MUMAS+2-J)
              I2=MIN0(MBB,N+MUMAS+1-J)
              DO 102 I=I1,I2
              IB=I+MDIFF
  102         E(IB,J)=E(IB,J)+FAC*FMAS(I,J)
              CALL DECB(N,LE,E,MLE,MUE,IP,INFO)
              IF (INFO.NE.0) GOTO 80
              IF (AUTNMS) THEN
! --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
! ---   1) THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
! ---   2) THE MATRIX B AND THE JACOBIAN OF F ARE BANDED
! ---   3) THE DIFFERENTIAL EQUATION IS AUTONOMOUS
                  DO 103 I=1,N
  103             AK1(I)=DY1(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK1,IP)
                  DO 110 I=1,N
  110             YNEW(I)=Y(I)+A21*AK1(I)
                      CALL FCN(N,X,YNEW,DY)
                  DO 111 I=1,N
  111             YNEW(I)=HC21*AK1(I)
                  DO 114 I=1,N
                  SUM=0.D0
                  J1=MAX0(1,I-MLMAS)
                  J2=MIN0(N,I+MUMAS)
                  DO 113 J=J1,J2
  113             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
  114             AK2(I)=SUM+DY(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK2,IP)
                  DO 120 I=1,N
  120             YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)
                      CALL FCN(N,X,YNEW,DY)
                  DO 121 I=1,N
  121             YNEW(I)=HC31*AK1(I)+HC32*AK2(I)
                  DO 124 I=1,N
                  SUM=0.D0
                  J1=MAX0(1,I-MLMAS)
                  J2=MIN0(N,I+MUMAS)
                  DO 123 J=J1,J2
  123             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
  124             AK3(I)=SUM+DY(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK3,IP)
                  DO 131 I=1,N
  131             YNEW(I)=HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I)
                  DO 134 I=1,N
                  SUM=0.D0
                  J1=MAX0(1,I-MLMAS)
                  J2=MIN0(N,I+MUMAS)
                  DO 133 J=J1,J2
  133             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
  134             AK4(I)=SUM+DY(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK4,IP)
              ELSE
! --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
! ---   1) THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
! ---   2) THE MATRIX B AND THE JACOBIAN OF F ARE BANDED
! ---   3) THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS
                  HD1=H*D1
                  HD2=H*D2
                  HD3=H*D3
                  HD4=H*D4
                  DO 203 I=1,N
  203             AK1(I)=DY1(I)+HD1*FX(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK1,IP)
                  DO 210 I=1,N
  210             YNEW(I)=Y(I)+A21*AK1(I)
                      CALL FCN(N,X+C2*H,YNEW,DY)
                  DO 211 I=1,N
  211             YNEW(I)=HC21*AK1(I)
                  DO 214 I=1,N
                  SUM=0.D0
                  J1=MAX0(1,I-MLMAS)
                  J2=MIN0(N,I+MUMAS)
                  DO 213 J=J1,J2
  213             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
  214             AK2(I)=SUM+DY(I)+HD2*FX(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK2,IP)
                  DO 220 I=1,N
  220             YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)
                      CALL FCN(N,X+C3*H,YNEW,DY)
                  DO 221 I=1,N
  221             YNEW(I)=HC31*AK1(I)+HC32*AK2(I)
                  DO 224 I=1,N
                  SUM=0.D0
                  J1=MAX0(1,I-MLMAS)
                  J2=MIN0(N,I+MUMAS)
                  DO 223 J=J1,J2
  223             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
  224             AK3(I)=SUM+DY(I)+HD3*FX(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK3,IP)
                  DO 231 I=1,N
  231             YNEW(I)=HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I)
                  DO 234 I=1,N
                  SUM=0.D0
                  J1=MAX0(1,I-MLMAS)
                  J2=MIN0(N,I+MUMAS)
                  DO 233 J=J1,J2
  233             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
  234             AK4(I)=SUM+DY(I)+HD4*FX(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK4,IP)
              END IF
          ELSE
              IF (MLMAS.NE.N) THEN
! --- THE MATRIX E (B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX)
                  MADD=MUMAS+1
                  DO 301 J=1,N
                  DO 301 I=1,N
  301             E(I,J)=-FJAC(I,J)
                  DO 302 J=1,N
                  I1=MAX0(1,J-MUMAS)
                  I2=MIN0(N,J+MLMAS)
                  DO 302 I=I1,I2
  302             E(I,J)=E(I,J)+FAC*FMAS(I-J+MADD,J)
                  CALL DEC(N,LE,E,IP,INFO)
                  IF (INFO.NE.0) GOTO 80
                  IF (AUTNMS) THEN
! --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
! ---   1) THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
! ---   2) THE MATRIX B IS BANDED BUT THE JACOBIAN OF F IS NOT
! ---   3) THE DIFFERENTIAL EQUATION IS AUTONOMOUS
                      DO 303 I=1,N
  303                 AK1(I)=DY1(I)
                          CALL SOL(N,LE,E,AK1,IP)
                      DO 310 I=1,N
  310                 YNEW(I)=Y(I)+A21*AK1(I)
                          CALL FCN(N,X,YNEW,DY)
                      DO 311 I=1,N
  311                 YNEW(I)=HC21*AK1(I)
                      DO 314 I=1,N
                      SUM=DY(I)
                      J1=MAX0(1,I-MLMAS)
                      J2=MIN0(N,I+MUMAS)
                      DO 313 J=J1,J2
  313                 SUM=SUM+FMAS(I-J+MADD,J)*YNEW(J)
  314                 AK2(I)=SUM
                          CALL SOL(N,LE,E,AK2,IP)
                      DO 320 I=1,N
  320                 YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)
                           CALL FCN(N,X,YNEW,DY)
                      DO 321 I=1,N
  321                 YNEW(I)=HC31*AK1(I)+HC32*AK2(I)
                      DO 324 I=1,N
                      SUM=DY(I)
                      J1=MAX0(1,I-MLMAS)
                      J2=MIN0(N,I+MUMAS)
                      DO 323 J=J1,J2
  323                 SUM=SUM+FMAS(I-J+MADD,J)*YNEW(J)
  324                 AK3(I)=SUM
                          CALL SOL(N,LE,E,AK3,IP)
                      DO 331 I=1,N
  331                 YNEW(I)=HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I)
                      DO 334 I=1,N
                      SUM=DY(I)
                      J1=MAX0(1,I-MLMAS)
                      J2=MIN0(N,I+MUMAS)
                      DO 333 J=J1,J2
  333                 SUM=SUM+FMAS(I-J+MADD,J)*YNEW(J)
  334                 AK4(I)=SUM
                          CALL SOL(N,LE,E,AK4,IP)
                  ELSE
! --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
! ---   1) THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
! ---   2) THE MATRIX B IS BANDED BUT THE JACOBIAN OF F IS NOT
! ---   3) THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS
                      HD1=H*D1
                      HD2=H*D2
                      HD3=H*D3
                      HD4=H*D4
                      DO 353 I=1,N
  353                 AK1(I)=DY1(I)+HD1*FX(I)
                          CALL SOL(N,LE,E,AK1,IP)
                      DO 360 I=1,N
  360                 YNEW(I)=Y(I)+A21*AK1(I)
                          CALL FCN(N,X+C2*H,YNEW,DY)
                      DO 361 I=1,N
  361                 YNEW(I)=HC21*AK1(I)
                      DO 364 I=1,N
                      SUM=DY(I)+HD2*FX(I)
                      J1=MAX0(1,I-MLMAS)
                      J2=MIN0(N,I+MUMAS)
                      DO 363 J=J1,J2
  363                 SUM=SUM+FMAS(I-J+MADD,J)*YNEW(J)
  364                 AK2(I)=SUM
                          CALL SOL(N,LE,E,AK2,IP)
                      DO 370 I=1,N
  370                 YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)
                           CALL FCN(N,X+C3*H,YNEW,DY)
                      DO 371 I=1,N
  371                 YNEW(I)=HC31*AK1(I)+HC32*AK2(I)
                      DO 374 I=1,N
                      SUM=DY(I)+HD3*FX(I)
                      J1=MAX0(1,I-MLMAS)
                      J2=MIN0(N,I+MUMAS)
                      DO 373 J=J1,J2
  373                 SUM=SUM+FMAS(I-J+MADD,J)*YNEW(J)
  374                 AK3(I)=SUM
                          CALL SOL(N,LE,E,AK3,IP)
                      DO 381 I=1,N
  381                 YNEW(I)=HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I)
                      DO 384 I=1,N
                      SUM=DY(I)+HD4*FX(I)
                      J1=MAX0(1,I-MLMAS)
                      J2=MIN0(N,I+MUMAS)
                      DO 383 J=J1,J2
  383                 SUM=SUM+FMAS(I-J+MADD,J)*YNEW(J)
  384                 AK4(I)=SUM
                          CALL SOL(N,LE,E,AK4,IP)
                  END IF
              ELSE
! --- THE MATRIX E (B IS A FULL MATRIX, JACOBIAN A FULL OR BANDED MATRIX)
                  IF (MLJAC.EQ.N) THEN
                      DO 401 J=1,N
                      DO 401 I=1,N
  401                 E(I,J)=FMAS(I,J)*FAC-FJAC(I,J)
                  ELSE
                      MADD=MUJAC+1
                      DO 405 J=1,N
                      DO 405 I=1,N
  405                 E(I,J)=FMAS(I,J)*FAC
                      DO 406 J=1,N
                      I1=MAX0(1,J-MUJAC)
                      I2=MIN0(N,J+MLJAC)
                      DO 406 I=I1,I2
  406                 E(I,J)=E(I,J)-FJAC(I-J+MADD,J)
                  END IF
                  CALL DEC(N,LE,E,IP,INFO)
                  IF (INFO.NE.0) GOTO 80
                  IF (AUTNMS) THEN
! --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
! ---   1) THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
! ---   2) THE MATRIX B IS NOT BANDED
! ---   3) THE DIFFERENTIAL EQUATION IS AUTONOMOUS
                      DO 403 I=1,N
  403                 AK1(I)=DY1(I)
                          CALL SOL(N,LE,E,AK1,IP)
                      DO 410 I=1,N
  410                 YNEW(I)=Y(I)+A21*AK1(I)
                          CALL FCN(N,X,YNEW,DY)
                      DO 411 I=1,N
  411                 YNEW(I)=HC21*AK1(I)
                      DO 414 I=1,N
                      SUM=DY(I)
                      DO 413 J=1,N
  413                 SUM=SUM+FMAS(I,J)*YNEW(J)
  414                 AK2(I)=SUM
                          CALL SOL(N,LE,E,AK2,IP)
                      DO 420 I=1,N
  420                 YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)
                           CALL FCN(N,X,YNEW,DY)
                      DO 421 I=1,N
  421                 YNEW(I)=HC31*AK1(I)+HC32*AK2(I)
                      DO 424 I=1,N
                      SUM=DY(I)
                      DO 423 J=1,N
  423                 SUM=SUM+FMAS(I,J)*YNEW(J)
  424                 AK3(I)=SUM
                          CALL SOL(N,LE,E,AK3,IP)
                      DO 431 I=1,N
  431                 YNEW(I)=HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I)
                      DO 434 I=1,N
                      SUM=DY(I)
                      DO 433 J=1,N
  433                 SUM=SUM+FMAS(I,J)*YNEW(J)
  434                 AK4(I)=SUM
                          CALL SOL(N,LE,E,AK4,IP)
                  ELSE
! --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
! ---   1) THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
! ---   2) THE MATRIX B IS NOT BANDED
! ---   3) THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS
                      HD1=H*D1
                      HD2=H*D2
                      HD3=H*D3
                      HD4=H*D4
                      DO 503 I=1,N
  503                 AK1(I)=DY1(I)+HD1*FX(I)
                          CALL SOL(N,LE,E,AK1,IP)
                      DO 510 I=1,N
  510                 YNEW(I)=Y(I)+A21*AK1(I)
                          CALL FCN(N,X+C2*H,YNEW,DY)
                      DO 511 I=1,N
  511                 YNEW(I)=HC21*AK1(I)
                      DO 514 I=1,N
                      SUM=DY(I)+HD2*FX(I)
                      DO 513 J=1,N
  513                 SUM=SUM+FMAS(I,J)*YNEW(J)
  514                 AK2(I)=SUM
                          CALL SOL(N,LE,E,AK2,IP)
                      DO 520 I=1,N
  520                 YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)
                           CALL FCN(N,X+C3*H,YNEW,DY)
                      DO 521 I=1,N
  521                 YNEW(I)=HC31*AK1(I)+HC32*AK2(I)
                      DO 524 I=1,N
                      SUM=DY(I)+HD3*FX(I)
                      DO 523 J=1,N
  523                 SUM=SUM+FMAS(I,J)*YNEW(J)
  524                 AK3(I)=SUM
                          CALL SOL(N,LE,E,AK3,IP)
                      DO 531 I=1,N
  531                 YNEW(I)=HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I)
                      DO 534 I=1,N
                      SUM=DY(I)+HD4*FX(I)
                      DO 533 J=1,N
  533                 SUM=SUM+FMAS(I,J)*YNEW(J)
  534                 AK4(I)=SUM
                          CALL SOL(N,LE,E,AK4,IP)
                  END IF
              END IF
          END IF
      ELSE
          IF (BANDED) THEN
! --- THE MATRIX E (B=IDENTITY, JACOBIAN A BANDED MATRIX)
              DO 601 J=1,N
              I1=MAX0(1,MUJAC+2-J)
              I2=MIN0(MBJAC,N+MUJAC+1-J)
              DO 600 I=I1,I2
  600         E(I+MLE,J)=-FJAC(I,J)
  601         E(MDIAG,J)=E(MDIAG,J)+FAC
              CALL DECB(N,LE,E,MLE,MUE,IP,INFO)
              IF (INFO.NE.0) GOTO 80
              IF (AUTNMS) THEN
! --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
! ---   1) THE DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
! ---   2) THE JACOBIAN OF THE PROBLEM IS A BANDED MATRIX
! ---   3) THE DIFFERENTIAL EQUATION IS AUTONOMOUS
                  DO 603 I=1,N
  603             AK1(I)=DY1(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK1,IP)
                  DO 610 I=1,N
  610             YNEW(I)=Y(I)+A21*AK1(I)
                  CALL FCN(N,X,YNEW,DY)
                  DO 611 I=1,N
  611             AK2(I)=DY(I)+HC21*AK1(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK2,IP)
                  DO 620 I=1,N
  620             YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)
                  CALL FCN(N,X,YNEW,DY)
                  DO 621 I=1,N
  621             AK3(I)=DY(I)+HC31*AK1(I)+HC32*AK2(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK3,IP)
                  DO 631 I=1,N
  631             AK4(I)=DY(I)+HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK4,IP)
              ELSE
! --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
! ---   1) THE DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
! ---   2) THE JACOBIAN OF THE PROBLEM IS A BANDED MATRIX
! ---   3) THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS
                  HD1=H*D1
                  HD2=H*D2
                  HD3=H*D3
                  HD4=H*D4
                  DO 703 I=1,N
  703             AK1(I)=DY1(I)+HD1*FX(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK1,IP)
                  DO 710 I=1,N
  710             YNEW(I)=Y(I)+A21*AK1(I)
                  CALL FCN(N,X+C2*H,YNEW,DY)
                  DO 711 I=1,N
  711             AK2(I)=DY(I)+HD2*FX(I)+HC21*AK1(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK2,IP)
                  DO 720 I=1,N
  720             YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)
                  CALL FCN(N,X+C3*H,YNEW,DY)
                  DO 721 I=1,N
  721             AK3(I)=DY(I)+HD3*FX(I)+HC31*AK1(I)+HC32*AK2(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK3,IP)
                  DO 731 I=1,N
  731             AK4(I)=DY(I)+HD4*FX(I)+HC41*AK1(I)+HC42*AK2(I)  &
                        +HC43*AK3(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK4,IP)
              END IF
          ELSE
! --- THE MATRIX E (B=IDENTITY, JACOBIAN A FULL MATRIX)
              DO 801 J=1,N
              DO 800 I=1,N
  800         E(I,J)=-FJAC(I,J)
  801         E(J,J)=E(J,J)+FAC
              CALL DEC(N,LE,E,IP,INFO)
              IF (INFO.NE.0) GOTO 80
              IF (AUTNMS) THEN
! --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
! ---   1) THE DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
! ---   2) THE JACOBIAN OF THE PROBLEM IS A FULL MATRIX
! ---   3) THE DIFFERENTIAL EQUATION IS AUTONOMOUS
                  DO 803 I=1,N
  803             AK1(I)=DY1(I)
                  CALL SOL(N,LE,E,AK1,IP)
                  DO 810 I=1,N
  810             YNEW(I)=Y(I)+A21*AK1(I)
                  CALL FCN(N,X,YNEW,DY)
                  DO 811 I=1,N
  811             AK2(I)=DY(I)+HC21*AK1(I)
                  CALL SOL(N,LE,E,AK2,IP)
                  DO 820 I=1,N
  820             YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)
                  CALL FCN(N,X,YNEW,DY)
                  DO 821 I=1,N
  821             AK3(I)=DY(I)+HC31*AK1(I)+HC32*AK2(I)
                  CALL SOL(N,LE,E,AK3,IP)
                  DO 831 I=1,N
  831             AK4(I)=DY(I)+HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I)
                  CALL SOL(N,LE,E,AK4,IP)
              ELSE
! --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
! ---   1) THE DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
! ---   2) THE JACOBIAN OF THE PROBLEM IS A FULL MATRIX
! ---   3) THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS
                  HD1=H*D1
                  HD2=H*D2
                  HD3=H*D3
                  HD4=H*D4
                  DO 903 I=1,N
  903             AK1(I)=DY1(I)+HD1*FX(I)
                  CALL SOL(N,LE,E,AK1,IP)
                  DO 910 I=1,N
  910             YNEW(I)=Y(I)+A21*AK1(I)
                  CALL FCN(N,X+C2*H,YNEW,DY)
                  DO 911 I=1,N
  911             AK2(I)=DY(I)+HD2*FX(I)+HC21*AK1(I)
                  CALL SOL(N,LE,E,AK2,IP)
                  DO 920 I=1,N
  920             YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)
                  CALL FCN(N,X+C3*H,YNEW,DY)
                  DO 921 I=1,N
  921             AK3(I)=DY(I)+HD3*FX(I)+HC31*AK1(I)+HC32*AK2(I)
                  CALL SOL(N,LE,E,AK3,IP)
                  DO 931 I=1,N
  931             AK4(I)=DY(I)+HD4*FX(I)+HC41*AK1(I)+HC42*AK2(I)  &
                        +HC43*AK3(I)
                  CALL SOL(N,LE,E,AK4,IP)
              END IF
          END IF
      END IF
      NSOL=NSOL+4
      NFCN=NFCN+2
! *** *** *** *** *** *** ***
!  ERROR ESTIMATION
! *** *** *** *** *** *** ***
      NSTEP=NSTEP+1
! ------------ NEW SOLUTION ---------------
      DO 240 I=1,N
  240 YNEW(I)=Y(I)+B1*AK1(I)+B2*AK2(I)+B3*AK3(I)+B4*AK4(I)
! ------------ COMPUTE ERROR ESTIMATION ----------------
      ERR=0.D0
      DO 300 I=1,N
      S=E1*AK1(I)+E2*AK2(I)+E3*AK3(I)+E4*AK4(I)
      IF (ITOL.EQ.0) THEN
         SK=ATOL(1)+RTOL(1)*DMAX1(DABS(Y(I)),DABS(YNEW(I)))
      ELSE
         SK=ATOL(I)+RTOL(I)*DMAX1(DABS(Y(I)),DABS(YNEW(I)))
      END IF
  300 ERR=ERR+(S/SK)**2
      ERR=DSQRT(ERR/N)
! --- COMPUTATION OF HNEW
! --- WE REQUIRE .2<=HNEW/H<=6.
      FAC=DMAX1(FAC2,DMIN1(FAC1,(ERR)**.25D0/.9D0))
      HNEW=H/FAC
! *** *** *** *** *** *** ***
!  IS THE ERROR SMALL ENOUGH ?
! *** *** *** *** *** *** ***
      IF (ERR.LE.1.D0) THEN
! --- STEP IS ACCEPTED
         NACCPT=NACCPT+1
         DO 44 I=1,N
  44     Y(I)=YNEW(I)
         XXOLD=X
         X=X+H
         IF (IOUT.NE.0) CALL SOLOUT(NACCPT+1,XXOLD,X,Y,N,IRTRN)
         IF (IRTRN.LT.0) GOTO 79
         IF (DABS(HNEW).GT.HMAXN) HNEW=POSNEG*HMAXN
         IF (REJECT) HNEW=POSNEG*DMIN1(DABS(HNEW),DABS(H))
         REJECT=.FALSE.
         RJECT2=.FALSE.
         H=HNEW
         GOTO 1
      ELSE
! --- STEP IS REJECTED
         IF (RJECT2) HNEW=H*FACREJ
         IF (REJECT) RJECT2=.TRUE.
         REJECT=.TRUE.
         H=HNEW
         IF (NACCPT.GE.1) NREJCT=NREJCT+1
         GOTO 2
      END IF
! --- EXIT
  80  WRITE (6,*) ' MATRIX E IS SINGULAR, INFO = ',INFO
      NSING=NSING+1
      IF (NSING.GE.5) GOTO 79
      H=H*0.5D0
      GOTO 2
  79  WRITE(6,979)X,H
 979  FORMAT(' EXIT OF ROS4 AT X=',D16.7,'   H=',D16.7)
      IDID=-1
      RETURN
      END

      SUBROUTINE SHAMP (A21,A31,A32,C21,C31,C32,C41,C42,C43,      &
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IMPLICIT REAL*8 (A-H,O-Z)
         A21=2.D0
         A31=48.D0/25.D0
         A32=6.D0/25.D0
         C21=-8.D0
         C31=372.D0/25.D0
         C32=12.D0/5.D0
         C41=-112.D0/125.D0
         C42=-54.D0/125.D0
         C43=-2.D0/5.D0
         B1=19.D0/9.D0
         B2=1.D0/2.D0
         B3=25.D0/108.D0
         B4=125.D0/108.D0
         E1=17.D0/54.D0
         E2=7.D0/36.D0
         E3=0.D0
         E4=125.D0/108.D0
         GAMMA=.5D0
         C2= 0.1000000000000000D+01
         C3= 0.6000000000000000D+00
         D1= 0.5000000000000000D+00
         D2=-0.1500000000000000D+01
         D3= 0.2420000000000000D+01
         D4= 0.1160000000000000D+00
      RETURN
      END
!
      SUBROUTINE GRK4A (A21,A31,A32,C21,C31,C32,C41,C42,C43,     &
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IMPLICIT REAL*8 (A-H,O-Z)
       A21= 0.1108860759493671D+01
       A31= 0.2377085261983360D+01
       A32= 0.1850114988899692D+00
       C21=-0.4920188402397641D+01
       C31= 0.1055588686048583D+01
       C32= 0.3351817267668938D+01
       C41= 0.3846869007049313D+01
       C42= 0.3427109241268180D+01
       C43=-0.2162408848753263D+01
       B1= 0.1845683240405840D+01
       B2= 0.1369796894360503D+00
       B3= 0.7129097783291559D+00
       B4= 0.6329113924050632D+00
       E1= 0.4831870177201765D-01
       E2=-0.6471108651049505D+00
       E3= 0.2186876660500240D+00
       E4=-0.6329113924050632D+00
       GAMMA= 0.3950000000000000D+00
       C2= 0.4380000000000000D+00
       C3= 0.8700000000000000D+00
       D1= 0.3950000000000000D+00
       D2=-0.3726723954840920D+00
       D3= 0.6629196544571492D-01
       D4= 0.4340946962568634D+00
      RETURN
      END
!
      SUBROUTINE GRK4T (A21,A31,A32,C21,C31,C32,C41,C42,C43,     &
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IMPLICIT REAL*8 (A-H,O-Z)
       A21= 0.2000000000000000D+01
       A31= 0.4524708207373116D+01
       A32= 0.4163528788597648D+01
       C21=-0.5071675338776316D+01
       C31= 0.6020152728650786D+01
       C32= 0.1597506846727117D+00
       C41=-0.1856343618686113D+01
       C42=-0.8505380858179826D+01
       C43=-0.2084075136023187D+01
       B1= 0.3957503746640777D+01
       B2= 0.4624892388363313D+01
       B3= 0.6174772638750108D+00
       B4= 0.1282612945269037D+01
       E1= 0.2302155402932996D+01
       E2= 0.3073634485392623D+01
       E3=-0.8732808018045032D+00
       E4=-0.1282612945269037D+01
       GAMMA= 0.2310000000000000D+00
       C2= 0.4620000000000000D+00
       C3= 0.8802083333333334D+00
       D1= 0.2310000000000000D+00
       D2=-0.3962966775244303D-01
       D3= 0.5507789395789127D+00
       D4=-0.5535098457052764D-01
      RETURN
      END
!
      SUBROUTINE VELDS (A21,A31,A32,C21,C31,C32,C41,C42,C43,     &
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
! --- METHOD GIVEN BY VAN VELDHUIZEN
      IMPLICIT REAL*8 (A-H,O-Z)
       A21= 0.2000000000000000D+01
       A31= 0.1750000000000000D+01
       A32= 0.2500000000000000D+00
       C21=-0.8000000000000000D+01
       C31=-0.8000000000000000D+01
       C32=-0.1000000000000000D+01
       C41= 0.5000000000000000D+00
       C42=-0.5000000000000000D+00
       C43= 0.2000000000000000D+01
       B1= 0.1333333333333333D+01
       B2= 0.6666666666666667D+00
       B3=-0.1333333333333333D+01
       B4= 0.1333333333333333D+01
       E1=-0.3333333333333333D+00
       E2=-0.3333333333333333D+00
       E3=-0.0000000000000000D+00
       E4=-0.1333333333333333D+01
       GAMMA= 0.5000000000000000D+00
       C2= 0.1000000000000000D+01
       C3= 0.5000000000000000D+00
       D1= 0.5000000000000000D+00
       D2=-0.1500000000000000D+01
       D3=-0.7500000000000000D+00
       D4= 0.2500000000000000D+00
      RETURN
      END
!
      SUBROUTINE VELDD (A21,A31,A32,C21,C31,C32,C41,C42,C43,     &
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
! --- METHOD GIVEN BY VAN VELDHUIZEN
      IMPLICIT REAL*8 (A-H,O-Z)
       A21= 0.2000000000000000D+01
       A31= 0.4812234362695436D+01
       A32= 0.4578146956747842D+01
       C21=-0.5333333333333331D+01
       C31= 0.6100529678848254D+01
       C32= 0.1804736797378427D+01
       C41=-0.2540515456634749D+01
       C42=-0.9443746328915205D+01
       C43=-0.1988471753215993D+01
       B1= 0.4289339254654537D+01
       B2= 0.5036098482851414D+01
       B3= 0.6085736420673917D+00
       B4= 0.1355958941201148D+01
       E1= 0.2175672787531755D+01
       E2= 0.2950911222575741D+01
       E3=-0.7859744544887430D+00
       E4=-0.1355958941201148D+01
       GAMMA= 0.2257081148225682D+00
       C2= 0.4514162296451364D+00
       C3= 0.8755928946018455D+00
       D1= 0.2257081148225682D+00
       D2=-0.4599403502680582D-01
       D3= 0.5177590504944076D+00
       D4=-0.3805623938054428D-01
      RETURN
      END
!
      SUBROUTINE LSTAB (A21,A31,A32,C21,C31,C32,C41,C42,C43,     &
                B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
! --- AN L-STABLE METHOD
      IMPLICIT REAL*8 (A-H,O-Z)
       A21= 0.2000000000000000D+01
       A31= 0.1867943637803922D+01
       A32= 0.2344449711399156D+00
       C21=-0.7137615036412310D+01
       C31= 0.2580708087951457D+01
       C32= 0.6515950076447975D+00
       C41=-0.2137148994382534D+01
       C42=-0.3214669691237626D+00
       C43=-0.6949742501781779D+00
       B1= 0.2255570073418735D+01
       B2= 0.2870493262186792D+00
       B3= 0.4353179431840180D+00
       B4= 0.1093502252409163D+01
       E1=-0.2815431932141155D+00
       E2=-0.7276199124938920D-01
       E3=-0.1082196201495311D+00
       E4=-0.1093502252409163D+01
       GAMMA= 0.5728200000000000D+00
       C2= 0.1145640000000000D+01
       C3= 0.6552168638155900D+00
       D1= 0.5728200000000000D+00
       D2=-0.1769193891319233D+01
       D3= 0.7592633437920482D+00
       D4=-0.1049021087100450D+00
      RETURN
      END

      SUBROUTINE DEC (N, NDIM, A, IP, IER)
! VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DOUBLE PRECISION A,T
      DIMENSION A(NDIM,N), IP(N)
!-----------------------------------------------------------------------
!  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
!  INPUT..
!     N = ORDER OF MATRIX.
!     NDIM = DECLARED DIMENSION OF ARRAY  A .
!     A = MATRIX TO BE TRIANGULARIZED.
!  OUTPUT..
!     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
!     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
!     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
!     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
!     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
!           SINGULAR AT STAGE K.
!  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
!  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
!  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
!
!  REFERENCE..
!     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
!     C.A.C.M. 15 (1972), P. 274.
!-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        DO 10 I = KP1,N
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        T = A(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(K,K)
        A(K,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = KP1,N
 30       A(I,K) = -A(I,K)*T
        DO 50 J = KP1,N
          T = A(M,J)
          A(M,J) = A(K,J)
          A(K,J) = T
          IF (T .EQ. 0.D0) GO TO 45
          DO 40 I = KP1,N
 40         A(I,J) = A(I,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(N,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
!----------------------- END OF SUBROUTINE DEC -------------------------
      END

      SUBROUTINE SOL (N, NDIM, A, B, IP)
! VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DOUBLE PRECISION A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
!-----------------------------------------------------------------------
!  SOLUTION OF LINEAR SYSTEM, A*X = B .
!  INPUT..
!    N = ORDER OF MATRIX.
!    NDIM = DECLARED DIMENSION OF ARRAY  A .
!    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
!    B = RIGHT HAND SIDE VECTOR.
!    IP = PIVOT VECTOR OBTAINED FROM DEC.
!  DO NOT USE IF DEC HAS SET IER <> 0.
!  OUTPUT..
!    B = SOLUTION VECTOR, X .
!-----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        DO 10 I = KP1,N
 10       B(I) = B(I) + A(I,K)*T
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        B(K) = B(K)/A(K,K)
        T = -B(K)
        DO 30 I = 1,KM1
 30       B(I) = B(I) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(1,1)
      RETURN
!----------------------- END OF SUBROUTINE SOL -------------------------
      END

      SUBROUTINE DECB (N, NDIM, A, ML, MU, IP, IER)
      REAL*8 A,T
      DIMENSION A(NDIM,N), IP(N)
!-----------------------------------------------------------------------
!  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED
!  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
!  INPUT..
!     N       ORDER OF THE ORIGINAL MATRIX A.
!     NDIM    DECLARED DIMENSION OF ARRAY  A.
!     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS
!                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND
!                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
!                ML+1 THROUGH 2*ML+MU+1 OF  A.
!     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!  OUTPUT..
!     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
!                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
!     IP      INDEX VECTOR OF PIVOT INDICES.
!     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
!     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
!                SINGULAR AT STAGE K.
!  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
!  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
!  IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO.
!
!  REFERENCE..
!     THIS IS A MODIFICATION OF
!     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
!     C.A.C.M. 15 (1972), P. 274.
!-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      MD = ML + MU + 1
      MD1 = MD + 1
      JU = 0
      IF (ML .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      IF (N .LT. MU+2) GO TO 7
      DO 5 J = MU+2,N
      DO 5 I = 1,ML
  5   A(I,J) = 0.D0
  7   NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = MD
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M + K - MD
        T = A(M,K)
        IF (M .EQ. MD) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(MD,K)
        A(MD,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = MD1,MDL
 30       A(I,K) = -A(I,K)*T
        JU = MIN0(MAX0(JU,MU+IP(K)),N)
        MM = MD
        IF (JU .LT. KP1) GO TO 55
        DO 50 J = KP1,JU
          M = M - 1
          MM = MM - 1
          T = A(M,J)
          IF (M .EQ. MM) GO TO 35
          A(M,J) = A(MM,J)
          A(MM,J) = T
 35       CONTINUE
          IF (T .EQ. 0.D0) GO TO 45
          JK = J - K
          DO 40 I = MD1,MDL
            IJK = I - JK
 40         A(IJK,J) = A(IJK,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 55     CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(MD,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
!----------------------- END OF SUBROUTINE DECB ------------------------
      END

      SUBROUTINE SOLB (N, NDIM, A, ML, MU, B, IP)
      REAL*8 A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
!-----------------------------------------------------------------------
!  SOLUTION OF LINEAR SYSTEM, A*X = B .
!  INPUT..
!    N      ORDER OF MATRIX A.
!    NDIM   DECLARED DIMENSION OF ARRAY  A .
!    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB.
!    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!    B      RIGHT HAND SIDE VECTOR.
!    IP     PIVOT VECTOR OBTAINED FROM DECB.
!  DO NOT USE IF DECB HAS SET IER  <> 0.
!  OUTPUT..
!    B      SOLUTION VECTOR, X .
!-----------------------------------------------------------------------
      MD = ML + MU + 1
      MD1 = MD + 1
      MDM = MD - 1
      NM1 = N - 1
      IF (ML .EQ. 0) GO TO 25
      IF (N .EQ. 1) GO TO 50
      DO 20 K = 1,NM1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IMD = I + K - MD
 10       B(IMD) = B(IMD) + A(I,K)*T
 20     CONTINUE
 25   CONTINUE
      DO 40 KB = 1,NM1
        K = N + 1 - KB
        B(K) = B(K)/A(MD,K)
        T = -B(K)
        KMD = MD - K
        LM = MAX0(1,KMD+1)
        DO 30 I = LM,MDM
          IMD = I - KMD
 30       B(IMD) = B(IMD) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(MD,1)
      RETURN
!----------------------- END OF SUBROUTINE SOLB ------------------------
      END

! end of file tros4.f90
