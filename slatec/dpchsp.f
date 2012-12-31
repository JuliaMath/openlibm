*DECK DPCHSP
      SUBROUTINE DPCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
C***BEGIN PROLOGUE  DPCHSP
C***PURPOSE  Set derivatives needed to determine the Hermite represen-
C            tation of the cubic spline interpolant to given data, with
C            specified boundary conditions.
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E1A
C***TYPE      DOUBLE PRECISION (PCHSP-S, DPCHSP-D)
C***KEYWORDS  CUBIC HERMITE INTERPOLATION, PCHIP,
C             PIECEWISE CUBIC INTERPOLATION, SPLINE INTERPOLATION
C***AUTHOR  Fritsch, F. N., (LLNL)
C             Lawrence Livermore National Laboratory
C             P.O. Box 808  (L-316)
C             Livermore, CA  94550
C             FTS 532-4275, (510) 422-4275
C***DESCRIPTION
C
C          DPCHSP:   Piecewise Cubic Hermite Spline
C
C     Computes the Hermite representation of the cubic spline inter-
C     polant to the data given in X and F satisfying the boundary
C     conditions specified by IC and VC.
C
C     To facilitate two-dimensional applications, includes an increment
C     between successive values of the F- and D-arrays.
C
C     The resulting piecewise cubic Hermite function may be evaluated
C     by DPCHFE or DPCHFD.
C
C     NOTE:  This is a modified version of C. de Boor's cubic spline
C            routine CUBSPL.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  IC(2), N, NWK, IERR
C        DOUBLE PRECISION  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
C
C        CALL  DPCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
C
C   Parameters:
C
C     IC -- (input) integer array of length 2 specifying desired
C           boundary conditions:
C           IC(1) = IBEG, desired condition at beginning of data.
C           IC(2) = IEND, desired condition at end of data.
C
C           IBEG = 0  to set D(1) so that the third derivative is con-
C              tinuous at X(2).  This is the "not a knot" condition
C              provided by de Boor's cubic spline routine CUBSPL.
C              < This is the default boundary condition. >
C           IBEG = 1  if first derivative at X(1) is given in VC(1).
C           IBEG = 2  if second derivative at X(1) is given in VC(1).
C           IBEG = 3  to use the 3-point difference formula for D(1).
C                     (Reverts to the default b.c. if N.LT.3 .)
C           IBEG = 4  to use the 4-point difference formula for D(1).
C                     (Reverts to the default b.c. if N.LT.4 .)
C          NOTES:
C           1. An error return is taken if IBEG is out of range.
C           2. For the "natural" boundary condition, use IBEG=2 and
C              VC(1)=0.
C
C           IEND may take on the same values as IBEG, but applied to
C           derivative at X(N).  In case IEND = 1 or 2, the value is
C           given in VC(2).
C
C          NOTES:
C           1. An error return is taken if IEND is out of range.
C           2. For the "natural" boundary condition, use IEND=2 and
C              VC(2)=0.
C
C     VC -- (input) real*8 array of length 2 specifying desired boundary
C           values, as indicated above.
C           VC(1) need be set only if IC(1) = 1 or 2 .
C           VC(2) need be set only if IC(2) = 1 or 2 .
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real*8 array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real*8 array of dependent variable values to be
C           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
C           X(I).
C
C     D -- (output) real*8 array of derivative values at the data
C           points.  These values will determine the cubic spline
C           interpolant with the requested boundary conditions.
C           The value corresponding to X(I) is stored in
C                D(1+(I-1)*INCFD),  I=1(1)N.
C           No other entries in D are changed.
C
C     INCFD -- (input) increment between successive values in F and D.
C           This argument is provided primarily for 2-D applications.
C           (Error return if  INCFD.LT.1 .)
C
C     WK -- (scratch) real*8 array of working storage.
C
C     NWK -- (input) length of work array.
C           (Error return if NWK.LT.2*N .)
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if IBEG.LT.0 or IBEG.GT.4 .
C              IERR = -5  if IEND.LT.0 of IEND.GT.4 .
C              IERR = -6  if both of the above are true.
C              IERR = -7  if NWK is too small.
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C             (The D-array has not been changed in any of these cases.)
C              IERR = -8  in case of trouble solving the linear system
C                         for the interior derivative values.
C             (The D-array may have been changed in this case.)
C             (             Do **NOT** use it!                )
C
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Springer-
C                 Verlag, New York, 1978, pp. 53-59.
C***ROUTINES CALLED  DPCHDF, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   820503  DATE WRITTEN
C   820804  Converted to SLATEC library version.
C   870707  Corrected XERROR calls for d.p. name(s).
C   890206  Corrected XERROR calls.
C   890411  Added SAVE statements (Vers. 3.2).
C   890703  Corrected category record.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920429  Revised format and order of references.  (WRB,FNF)
C***END PROLOGUE  DPCHSP
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DPCHSP to PCHSP wherever it occurs,
C        b. Change the double precision declarations to real, and
C        c. Change the constants ZERO, HALF, ... to single precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  IC(2), N, INCFD, NWK, IERR
      DOUBLE PRECISION  VC(2), X(*), F(INCFD,*), D(INCFD,*), WK(2,*)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  IBEG, IEND, INDEX, J, NM1
      DOUBLE PRECISION  G, HALF, ONE, STEMP(3), THREE, TWO, XTEMP(4),
     *  ZERO
      SAVE ZERO, HALF, ONE, TWO, THREE
      DOUBLE PRECISION  DPCHDF
C
      DATA  ZERO /0.D0/, HALF/.5D0/, ONE/1.D0/, TWO/2.D0/, THREE/3.D0/
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  DPCHSP
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  J = 2, N
         IF ( X(J).LE.X(J-1) )  GO TO 5003
    1 CONTINUE
C
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
      IF ( (IBEG.LT.0).OR.(IBEG.GT.4) )  IERR = IERR - 1
      IF ( (IEND.LT.0).OR.(IEND.GT.4) )  IERR = IERR - 2
      IF ( IERR.LT.0 )  GO TO 5004
C
C  FUNCTION DEFINITION IS OK -- GO ON.
C
      IF ( NWK .LT. 2*N )  GO TO 5007
C
C  COMPUTE FIRST DIFFERENCES OF X SEQUENCE AND STORE IN WK(1,.). ALSO,
C  COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN WK(2,.).
      DO 5  J=2,N
         WK(1,J) = X(J) - X(J-1)
         WK(2,J) = (F(1,J) - F(1,J-1))/WK(1,J)
    5 CONTINUE
C
C  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
C
      IF ( IBEG.GT.N )  IBEG = 0
      IF ( IEND.GT.N )  IEND = 0
C
C  SET UP FOR BOUNDARY CONDITIONS.
C
      IF ( (IBEG.EQ.1).OR.(IBEG.EQ.2) )  THEN
         D(1,1) = VC(1)
      ELSE IF (IBEG .GT. 2)  THEN
C        PICK UP FIRST IBEG POINTS, IN REVERSE ORDER.
         DO 10  J = 1, IBEG
            INDEX = IBEG-J+1
C           INDEX RUNS FROM IBEG DOWN TO 1.
            XTEMP(J) = X(INDEX)
            IF (J .LT. IBEG)  STEMP(J) = WK(2,INDEX)
   10    CONTINUE
C                 --------------------------------
         D(1,1) = DPCHDF (IBEG, XTEMP, STEMP, IERR)
C                 --------------------------------
         IF (IERR .NE. 0)  GO TO 5009
         IBEG = 1
      ENDIF
C
      IF ( (IEND.EQ.1).OR.(IEND.EQ.2) )  THEN
         D(1,N) = VC(2)
      ELSE IF (IEND .GT. 2)  THEN
C        PICK UP LAST IEND POINTS.
         DO 15  J = 1, IEND
            INDEX = N-IEND+J
C           INDEX RUNS FROM N+1-IEND UP TO N.
            XTEMP(J) = X(INDEX)
            IF (J .LT. IEND)  STEMP(J) = WK(2,INDEX+1)
   15    CONTINUE
C                 --------------------------------
         D(1,N) = DPCHDF (IEND, XTEMP, STEMP, IERR)
C                 --------------------------------
         IF (IERR .NE. 0)  GO TO 5009
         IEND = 1
      ENDIF
C
C --------------------( BEGIN CODING FROM CUBSPL )--------------------
C
C  **** A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(J) OF
C  F  AT X(J), J=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS ELIM-
C  INATION, WITH S(J) ENDING UP IN D(1,J), ALL J.
C     WK(1,.) AND WK(2,.) ARE USED FOR TEMPORARY STORAGE.
C
C  CONSTRUCT FIRST EQUATION FROM FIRST BOUNDARY CONDITION, OF THE FORM
C             WK(2,1)*S(1) + WK(1,1)*S(2) = D(1,1)
C
      IF (IBEG .EQ. 0)  THEN
         IF (N .EQ. 2)  THEN
C           NO CONDITION AT LEFT END AND N = 2.
            WK(2,1) = ONE
            WK(1,1) = ONE
            D(1,1) = TWO*WK(2,2)
         ELSE
C           NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
            WK(2,1) = WK(1,3)
            WK(1,1) = WK(1,2) + WK(1,3)
            D(1,1) =((WK(1,2) + TWO*WK(1,1))*WK(2,2)*WK(1,3)
     *                        + WK(1,2)**2*WK(2,3)) / WK(1,1)
         ENDIF
      ELSE IF (IBEG .EQ. 1)  THEN
C        SLOPE PRESCRIBED AT LEFT END.
         WK(2,1) = ONE
         WK(1,1) = ZERO
      ELSE
C        SECOND DERIVATIVE PRESCRIBED AT LEFT END.
         WK(2,1) = TWO
         WK(1,1) = ONE
         D(1,1) = THREE*WK(2,2) - HALF*WK(1,2)*D(1,1)
      ENDIF
C
C  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESPONDING EQUATIONS AND
C  CARRY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE J-TH
C  EQUATION READS    WK(2,J)*S(J) + WK(1,J)*S(J+1) = D(1,J).
C
      NM1 = N-1
      IF (NM1 .GT. 1)  THEN
         DO 20 J=2,NM1
            IF (WK(2,J-1) .EQ. ZERO)  GO TO 5008
            G = -WK(1,J+1)/WK(2,J-1)
            D(1,J) = G*D(1,J-1)
     *                  + THREE*(WK(1,J)*WK(2,J+1) + WK(1,J+1)*WK(2,J))
            WK(2,J) = G*WK(1,J-1) + TWO*(WK(1,J) + WK(1,J+1))
   20    CONTINUE
      ENDIF
C
C  CONSTRUCT LAST EQUATION FROM SECOND BOUNDARY CONDITION, OF THE FORM
C           (-G*WK(2,N-1))*S(N-1) + WK(2,N)*S(N) = D(1,N)
C
C     IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
C     SUBSTITUTION, SINCE ARRAYS HAPPEN TO BE SET UP JUST RIGHT FOR IT
C     AT THIS POINT.
      IF (IEND .EQ. 1)  GO TO 30
C
      IF (IEND .EQ. 0)  THEN
         IF (N.EQ.2 .AND. IBEG.EQ.0)  THEN
C           NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
            D(1,2) = WK(2,2)
            GO TO 30
         ELSE IF ((N.EQ.2) .OR. (N.EQ.3 .AND. IBEG.EQ.0))  THEN
C           EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND *NOT*
C           NOT-A-KNOT AT LEFT END POINT).
            D(1,N) = TWO*WK(2,N)
            WK(2,N) = ONE
            IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
            G = -ONE/WK(2,N-1)
         ELSE
C           NOT-A-KNOT AND N .GE. 3, AND EITHER N.GT.3 OR  ALSO NOT-A-
C           KNOT AT LEFT END POINT.
            G = WK(1,N-1) + WK(1,N)
C           DO NOT NEED TO CHECK FOLLOWING DENOMINATORS (X-DIFFERENCES).
            D(1,N) = ((WK(1,N)+TWO*G)*WK(2,N)*WK(1,N-1)
     *                  + WK(1,N)**2*(F(1,N-1)-F(1,N-2))/WK(1,N-1))/G
            IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
            G = -G/WK(2,N-1)
            WK(2,N) = WK(1,N-1)
         ENDIF
      ELSE
C        SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
         D(1,N) = THREE*WK(2,N) + HALF*WK(1,N)*D(1,N)
         WK(2,N) = TWO
         IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
         G = -ONE/WK(2,N-1)
      ENDIF
C
C  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
C
      WK(2,N) = G*WK(1,N-1) + WK(2,N)
      IF (WK(2,N) .EQ. ZERO)   GO TO 5008
      D(1,N) = (G*D(1,N-1) + D(1,N))/WK(2,N)
C
C  CARRY OUT BACK SUBSTITUTION
C
   30 CONTINUE
      DO 40 J=NM1,1,-1
         IF (WK(2,J) .EQ. ZERO)  GO TO 5008
         D(1,J) = (D(1,J) - WK(1,J)*D(1,J+1))/WK(2,J)
   40 CONTINUE
C --------------------(  END  CODING FROM CUBSPL )--------------------
C
C  NORMAL RETURN.
C
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERMSG ('SLATEC', 'DPCHSP',
     +   'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERMSG ('SLATEC', 'DPCHSP', 'INCREMENT LESS THAN ONE', IERR,
     +   1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERMSG ('SLATEC', 'DPCHSP',
     +   'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
      RETURN
C
 5004 CONTINUE
C     IC OUT OF RANGE RETURN.
      IERR = IERR - 3
      CALL XERMSG ('SLATEC', 'DPCHSP', 'IC OUT OF RANGE', IERR, 1)
      RETURN
C
 5007 CONTINUE
C     NWK TOO SMALL RETURN.
      IERR = -7
      CALL XERMSG ('SLATEC', 'DPCHSP', 'WORK ARRAY TOO SMALL', IERR, 1)
      RETURN
C
 5008 CONTINUE
C     SINGULAR SYSTEM.
C   *** THEORETICALLY, THIS CAN ONLY OCCUR IF SUCCESSIVE X-VALUES   ***
C   *** ARE EQUAL, WHICH SHOULD ALREADY HAVE BEEN CAUGHT (IERR=-3). ***
      IERR = -8
      CALL XERMSG ('SLATEC', 'DPCHSP', 'SINGULAR LINEAR SYSTEM', IERR,
     +   1)
      RETURN
C
 5009 CONTINUE
C     ERROR RETURN FROM DPCHDF.
C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -9
      CALL XERMSG ('SLATEC', 'DPCHSP', 'ERROR RETURN FROM DPCHDF',
     +   IERR, 1)
      RETURN
C------------- LAST LINE OF DPCHSP FOLLOWS -----------------------------
      END
