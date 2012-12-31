*DECK DPCHKT
      SUBROUTINE DPCHKT (N, X, KNOTYP, T)
C***BEGIN PROLOGUE  DPCHKT
C***SUBSIDIARY
C***PURPOSE  Compute B-spline knot sequence for DPCHBS.
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E3
C***TYPE      DOUBLE PRECISION (PCHKT-S, DPCHKT-D)
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C     Set a knot sequence for the B-spline representation of a PCH
C     function with breakpoints X.  All knots will be at least double.
C     Endknots are set as:
C        (1) quadruple knots at endpoints if KNOTYP=0;
C        (2) extrapolate the length of end interval if KNOTYP=1;
C        (3) periodic if KNOTYP=2.
C
C  Input arguments:  N, X, KNOTYP.
C  Output arguments:  T.
C
C  Restrictions/assumptions:
C     1. N.GE.2 .  (not checked)
C     2. X(i).LT.X(i+1), i=1,...,N .  (not checked)
C     3. 0.LE.KNOTYP.LE.2 .  (Acts like KNOTYP=0 for any other value.)
C
C***SEE ALSO  DPCHBS
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   870701  DATE WRITTEN
C   900405  Converted Fortran to upper case.
C   900410  Converted prologue to SLATEC 4.0 format.
C   900410  Minor cosmetic changes.
C   900430  Produced double precision version.
C   930514  Changed NKNOTS from an output to an input variable.  (FNF)
C   930604  Removed unused variable NKNOTS from argument list.  (FNF)
C***END PROLOGUE  DPCHKT
C
C*Internal Notes:
C
C  Since this is subsidiary to DPCHBS, which validates its input before
C  calling, it is unnecessary for such validation to be done here.
C
C**End
C
C  Declare arguments.
C
      INTEGER  N, KNOTYP
      DOUBLE PRECISION  X(*), T(*)
C
C  Declare local variables.
C
      INTEGER  J, K, NDIM
      DOUBLE PRECISION  HBEG, HEND
C***FIRST EXECUTABLE STATEMENT  DPCHKT
C
C  Initialize.
C
      NDIM = 2*N
C
C  Set interior knots.
C
      J = 1
      DO 20  K = 1, N
         J = J + 2
         T(J) = X(K)
         T(J+1) = T(J)
   20 CONTINUE
C     Assertion:  At this point T(3),...,T(NDIM+2) have been set and
C                 J=NDIM+1.
C
C  Set end knots according to KNOTYP.
C
      HBEG = X(2) - X(1)
      HEND = X(N) - X(N-1)
      IF (KNOTYP.EQ.1 )  THEN
C          Extrapolate.
         T(2) = X(1) - HBEG
         T(NDIM+3) = X(N) + HEND
      ELSE IF ( KNOTYP.EQ.2 )  THEN
C          Periodic.
         T(2) = X(1) - HEND
         T(NDIM+3) = X(N) + HBEG
      ELSE
C          Quadruple end knots.
         T(2) = X(1)
         T(NDIM+3) = X(N)
      ENDIF
      T(1) = T(2)
      T(NDIM+4) = T(NDIM+3)
C
C  Terminate.
C
      RETURN
C------------- LAST LINE OF DPCHKT FOLLOWS -----------------------------
      END
