*DECK DPCHCS
      SUBROUTINE DPCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
C***BEGIN PROLOGUE  DPCHCS
C***SUBSIDIARY
C***PURPOSE  Adjusts derivative values for DPCHIC
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      DOUBLE PRECISION (PCHCS-S, DPCHCS-D)
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C         DPCHCS:  DPCHIC Monotonicity Switch Derivative Setter.
C
C     Called by  DPCHIC  to adjust the values of D in the vicinity of a
C     switch in direction of monotonicity, to produce a more "visually
C     pleasing" curve than that given by  DPCHIM .
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, IERR
C        DOUBLE PRECISION  SWITCH, H(N), SLOPE(N), D(INCFD,N)
C
C        CALL  DPCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
C
C   Parameters:
C
C     SWITCH -- (input) indicates the amount of control desired over
C           local excursions from data.
C
C     N -- (input) number of data points.  (assumes N.GT.2 .)
C
C     H -- (input) real*8 array of interval lengths.
C     SLOPE -- (input) real*8 array of data slopes.
C           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
C                  H(I) =  X(I+1)-X(I),
C              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
C
C     D -- (input) real*8 array of derivative values at the data points,
C           as determined by DPCHCI.
C          (output) derivatives in the vicinity of switches in direction
C           of monotonicity may be adjusted to produce a more "visually
C           pleasing" curve.
C           The value corresponding to X(I) is stored in
C                D(1+(I-1)*INCFD),  I=1(1)N.
C           No other entries in D are changed.
C
C     INCFD -- (input) increment between successive values in D.
C           This argument is provided primarily for 2-D applications.
C
C     IERR -- (output) error flag.  should be zero.
C           If negative, trouble in DPCHSW.  (should never happen.)
C
C    -------
C    WARNING:  This routine does no validity-checking of arguments.
C    -------
C
C  Fortran intrinsics used:  ABS, MAX, MIN.
C
C***SEE ALSO  DPCHIC
C***ROUTINES CALLED  DPCHST, DPCHSW
C***REVISION HISTORY  (YYMMDD)
C   820218  DATE WRITTEN
C   820617  Redesigned to (1) fix  problem with lack of continuity
C           approaching a flat-topped peak (2) be cleaner and
C           easier to verify.
C           Eliminated subroutines PCHSA and PCHSX in the process.
C   820622  1. Limited fact to not exceed one, so computed D is a
C             convex combination of DPCHCI value and DPCHSD value.
C           2. Changed fudge from 1 to 4 (based on experiments).
C   820623  Moved PCHSD to an inline function (eliminating MSWTYP).
C   820805  Converted to SLATEC library version.
C   870707  Corrected conversion to double precision.
C   870813  Minor cosmetic changes.
C   890411  Added SAVE statements (Vers. 3.2).
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891006  Modified spacing in computation of DFLOC.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated AUTHOR section in prologue.  (WRB)
C   930503  Improved purpose.  (FNF)
C***END PROLOGUE  DPCHCS
C
C  Programming notes:
C     1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
C        either argument is zero, +1 if they are of the same sign, and
C        -1 if they are of opposite sign.
C**End
C
C  DECLARE ARGUMENTS.
C
      INTEGER  N, INCFD, IERR
      DOUBLE PRECISION  SWITCH, H(*), SLOPE(*), D(INCFD,*)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, INDX, K, NLESS1
      DOUBLE PRECISION  DEL(3), DEXT, DFLOC, DFMX, FACT, FUDGE, ONE,
     *      SLMAX, WTAVE(2), ZERO
      SAVE ZERO, ONE, FUDGE
      DOUBLE PRECISION  DPCHST
C
C  DEFINE INLINE FUNCTION FOR WEIGHTED AVERAGE OF SLOPES.
C
      DOUBLE PRECISION  DPCHSD, S1, S2, H1, H2
      DPCHSD(S1,S2,H1,H2) = (H2/(H1+H2))*S1 + (H1/(H1+H2))*S2
C
C  INITIALIZE.
C
      DATA  ZERO /0.D0/,  ONE/1.D0/
      DATA  FUDGE /4.D0/
C***FIRST EXECUTABLE STATEMENT  DPCHCS
      IERR = 0
      NLESS1 = N - 1
C
C  LOOP OVER SEGMENTS.
C
      DO 900  I = 2, NLESS1
         IF ( DPCHST(SLOPE(I-1),SLOPE(I)) )  100, 300, 900
C             --------------------------
C
  100    CONTINUE
C
C....... SLOPE SWITCHES MONOTONICITY AT I-TH POINT .....................
C
C           DO NOT CHANGE D IF 'UP-DOWN-UP'.
            IF (I .GT. 2)  THEN
               IF ( DPCHST(SLOPE(I-2),SLOPE(I)) .GT. ZERO)  GO TO 900
C                   --------------------------
            ENDIF
            IF (I .LT. NLESS1)  THEN
               IF ( DPCHST(SLOPE(I+1),SLOPE(I-1)) .GT. ZERO)  GO TO 900
C                   ----------------------------
            ENDIF
C
C   ....... COMPUTE PROVISIONAL VALUE FOR D(1,I).
C
            DEXT = DPCHSD (SLOPE(I-1), SLOPE(I), H(I-1), H(I))
C
C   ....... DETERMINE WHICH INTERVAL CONTAINS THE EXTREMUM.
C
            IF ( DPCHST(DEXT, SLOPE(I-1)) )  200, 900, 250
C                -----------------------
C
  200       CONTINUE
C              DEXT AND SLOPE(I-1) HAVE OPPOSITE SIGNS --
C                        EXTREMUM IS IN (X(I-1),X(I)).
               K = I-1
C              SET UP TO COMPUTE NEW VALUES FOR D(1,I-1) AND D(1,I).
               WTAVE(2) = DEXT
               IF (K .GT. 1)
     *            WTAVE(1) = DPCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
               GO TO 400
C
  250       CONTINUE
C              DEXT AND SLOPE(I) HAVE OPPOSITE SIGNS --
C                        EXTREMUM IS IN (X(I),X(I+1)).
               K = I
C              SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
               WTAVE(1) = DEXT
               IF (K .LT. NLESS1)
     *            WTAVE(2) = DPCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))
               GO TO 400
C
  300    CONTINUE
C
C....... AT LEAST ONE OF SLOPE(I-1) AND SLOPE(I) IS ZERO --
C                     CHECK FOR FLAT-TOPPED PEAK .......................
C
            IF (I .EQ. NLESS1)  GO TO 900
            IF ( DPCHST(SLOPE(I-1), SLOPE(I+1)) .GE. ZERO)  GO TO 900
C                -----------------------------
C
C           WE HAVE FLAT-TOPPED PEAK ON (X(I),X(I+1)).
            K = I
C           SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
            WTAVE(1) = DPCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
            WTAVE(2) = DPCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))
C
  400    CONTINUE
C
C....... AT THIS POINT WE HAVE DETERMINED THAT THERE WILL BE AN EXTREMUM
C        ON (X(K),X(K+1)), WHERE K=I OR I-1, AND HAVE SET ARRAY WTAVE--
C           WTAVE(1) IS A WEIGHTED AVERAGE OF SLOPE(K-1) AND SLOPE(K),
C                    IF K.GT.1
C           WTAVE(2) IS A WEIGHTED AVERAGE OF SLOPE(K) AND SLOPE(K+1),
C                    IF K.LT.N-1
C
         SLMAX = ABS(SLOPE(K))
         IF (K .GT. 1)    SLMAX = MAX( SLMAX, ABS(SLOPE(K-1)) )
         IF (K.LT.NLESS1) SLMAX = MAX( SLMAX, ABS(SLOPE(K+1)) )
C
         IF (K .GT. 1)  DEL(1) = SLOPE(K-1) / SLMAX
         DEL(2) = SLOPE(K) / SLMAX
         IF (K.LT.NLESS1)  DEL(3) = SLOPE(K+1) / SLMAX
C
         IF ((K.GT.1) .AND. (K.LT.NLESS1))  THEN
C           NORMAL CASE -- EXTREMUM IS NOT IN A BOUNDARY INTERVAL.
            FACT = FUDGE* ABS(DEL(3)*(DEL(1)-DEL(2))*(WTAVE(2)/SLMAX))
            D(1,K) = D(1,K) + MIN(FACT,ONE)*(WTAVE(1) - D(1,K))
            FACT = FUDGE* ABS(DEL(1)*(DEL(3)-DEL(2))*(WTAVE(1)/SLMAX))
            D(1,K+1) = D(1,K+1) + MIN(FACT,ONE)*(WTAVE(2) - D(1,K+1))
         ELSE
C           SPECIAL CASE K=1 (WHICH CAN OCCUR ONLY IF I=2) OR
C                        K=NLESS1 (WHICH CAN OCCUR ONLY IF I=NLESS1).
            FACT = FUDGE* ABS(DEL(2))
            D(1,I) = MIN(FACT,ONE) * WTAVE(I-K+1)
C              NOTE THAT I-K+1 = 1 IF K=I  (=NLESS1),
C                        I-K+1 = 2 IF K=I-1(=1).
         ENDIF
C
C
C....... ADJUST IF NECESSARY TO LIMIT EXCURSIONS FROM DATA.
C
         IF (SWITCH .LE. ZERO)  GO TO 900
C
         DFLOC = H(K)*ABS(SLOPE(K))
         IF (K .GT. 1)    DFLOC = MAX( DFLOC, H(K-1)*ABS(SLOPE(K-1)) )
         IF (K.LT.NLESS1) DFLOC = MAX( DFLOC, H(K+1)*ABS(SLOPE(K+1)) )
         DFMX = SWITCH*DFLOC
         INDX = I-K+1
C        INDX = 1 IF K=I, 2 IF K=I-1.
C        ---------------------------------------------------------------
         CALL DPCHSW(DFMX, INDX, D(1,K), D(1,K+1), H(K), SLOPE(K), IERR)
C        ---------------------------------------------------------------
         IF (IERR .NE. 0)  RETURN
C
C....... END OF SEGMENT LOOP.
C
  900 CONTINUE
C
      RETURN
C------------- LAST LINE OF DPCHCS FOLLOWS -----------------------------
      END
