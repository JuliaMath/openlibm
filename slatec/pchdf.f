*DECK PCHDF
      REAL FUNCTION PCHDF (K, X, S, IERR)
C***BEGIN PROLOGUE  PCHDF
C***SUBSIDIARY
C***PURPOSE  Computes divided differences for PCHCE and PCHSP
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      SINGLE PRECISION (PCHDF-S, DPCHDF-D)
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C          PCHDF:   PCHIP Finite Difference Formula
C
C     Uses a divided difference formulation to compute a K-point approx-
C     imation to the derivative at X(K) based on the data in X and S.
C
C     Called by  PCHCE  and  PCHSP  to compute 3- and 4-point boundary
C     derivative approximations.
C
C ----------------------------------------------------------------------
C
C     On input:
C        K      is the order of the desired derivative approximation.
C               K must be at least 3 (error return if not).
C        X      contains the K values of the independent variable.
C               X need not be ordered, but the values **MUST** be
C               distinct.  (Not checked here.)
C        S      contains the associated slope values:
C                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
C               (Note that S need only be of length K-1.)
C
C     On return:
C        S      will be destroyed.
C        IERR   will be set to -1 if K.LT.2 .
C        PCHDF  will be set to the desired derivative approximation if
C               IERR=0 or to zero if IERR=-1.
C
C ----------------------------------------------------------------------
C
C***SEE ALSO  PCHCE, PCHSP
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Springer-
C                 Verlag, New York, 1978, pp. 10-16.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   820503  DATE WRITTEN
C   820805  Converted to SLATEC library version.
C   870813  Minor cosmetic changes.
C   890411  Added SAVE statements (Vers. 3.2).
C   890411  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
C   920429  Revised format and order of references.  (WRB,FNF)
C   930503  Improved purpose.  (FNF)
C***END PROLOGUE  PCHDF
C
C**End
C
C  DECLARE ARGUMENTS.
C
      INTEGER  K, IERR
      REAL  X(K), S(K)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, J
      REAL  VALUE, ZERO
      SAVE ZERO
      DATA  ZERO /0./
C
C  CHECK FOR LEGAL VALUE OF K.
C
C***FIRST EXECUTABLE STATEMENT  PCHDF
      IF (K .LT. 3)  GO TO 5001
C
C  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
C
      DO 10  J = 2, K-1
         DO 9  I = 1, K-J
            S(I) = (S(I+1)-S(I))/(X(I+J)-X(I))
    9    CONTINUE
   10 CONTINUE
C
C  EVALUATE DERIVATIVE AT X(K).
C
      VALUE = S(1)
      DO 20  I = 2, K-1
         VALUE = S(I) + VALUE*(X(K)-X(I))
   20 CONTINUE
C
C  NORMAL RETURN.
C
      IERR = 0
      PCHDF = VALUE
      RETURN
C
C  ERROR RETURN.
C
 5001 CONTINUE
C     K.LT.3 RETURN.
      IERR = -1
      CALL XERMSG ('SLATEC', 'PCHDF', 'K LESS THAN THREE', IERR, 1)
      PCHDF = ZERO
      RETURN
C------------- LAST LINE OF PCHDF FOLLOWS ------------------------------
      END
