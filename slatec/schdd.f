*DECK SCHDD
      SUBROUTINE SCHDD (R, LDR, P, X, Z, LDZ, NZ, Y, RHO, C, S, INFO)
C***BEGIN PROLOGUE  SCHDD
C***PURPOSE  Downdate an augmented Cholesky decomposition or the
C            triangular factor of an augmented QR decomposition.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D7B
C***TYPE      SINGLE PRECISION (SCHDD-S, DCHDD-D, CCHDD-C)
C***KEYWORDS  CHOLESKY DECOMPOSITION, DOWNDATE, LINEAR ALGEBRA, LINPACK,
C             MATRIX
C***AUTHOR  Stewart, G. W., (U. of Maryland)
C***DESCRIPTION
C
C     SCHDD downdates an augmented Cholesky decomposition or the
C     triangular factor of an augmented QR decomposition.
C     Specifically, given an upper triangular matrix R of order P, a
C     row vector X, a column vector Z, and a scalar Y, SCHDD
C     determines an orthogonal matrix U and a scalar ZETA such that
C
C                        (R   Z )     (RR  ZZ)
C                    U * (      )  =  (      ) ,
C                        (0 ZETA)     ( X   Y)
C
C     where RR is upper triangular.  If R and Z have been obtained
C     from the factorization of a least squares problem, then
C     RR and ZZ are the factors corresponding to the problem
C     with the observation (X,Y) removed.  In this case, if RHO
C     is the norm of the residual vector, then the norm of
C     the residual vector of the downdated problem is
C     SQRT(RHO**2 - ZETA**2). SCHDD will simultaneously downdate
C     several triplets (Z,Y,RHO) along with R.
C     For a less terse description of what SCHDD does and how
C     it may be applied, see the LINPACK guide.
C
C     The matrix U is determined as the product U(1)*...*U(P)
C     where U(I) is a rotation in the (P+1,I)-plane of the
C     form
C
C                       ( C(I)     -S(I)     )
C                       (                    ) .
C                       ( S(I)       C(I)    )
C
C     The rotations are chosen so that C(I) is real.
C
C     The user is warned that a given downdating problem may
C     be impossible to accomplish or may produce
C     inaccurate results.  For example, this can happen
C     if X is near a vector whose removal will reduce the
C     rank of R.  Beware.
C
C     On Entry
C
C         R      REAL(LDR,P), where LDR .GE. P.
C                R contains the upper triangular matrix
C                that is to be downdated.  The part of  R
C                below the diagonal is not referenced.
C
C         LDR    INTEGER.
C                LDR is the leading dimension of the array R.
C
C         P      INTEGER.
C                P is the order of the matrix R.
C
C         X      REAL(P).
C                X contains the row vector that is to
C                be removed from R.  X is not altered by SCHDD.
C
C         Z      REAL(LDZ,NZ), where LDZ .GE. P.
C                Z is an array of NZ P-vectors which
C                are to be downdated along with R.
C
C         LDZ    INTEGER.
C                LDZ is the leading dimension of the array Z.
C
C         NZ     INTEGER.
C                NZ is the number of vectors to be downdated
C                NZ may be zero, in which case Z, Y, and RHO
C                are not referenced.
C
C         Y      REAL(NZ).
C                Y contains the scalars for the downdating
C                of the vectors Z.  Y is not altered by SCHDD.
C
C         RHO    REAL(NZ).
C                RHO contains the norms of the residual
C                vectors that are to be downdated.
C
C     On Return
C
C         R
C         Z      contain the downdated quantities.
C         RHO
C
C         C      REAL(P).
C                C contains the cosines of the transforming
C                rotations.
C
C         S      REAL(P).
C                S contains the sines of the transforming
C                rotations.
C
C         INFO   INTEGER.
C                INFO is set as follows.
C
C                   INFO = 0  if the entire downdating
C                             was successful.
C
C                   INFO =-1  if R could not be downdated.
C                             In this case, all quantities
C                             are left unaltered.
C
C                   INFO = 1  if some RHO could not be
C                             downdated.  The offending RHOs are
C                             set to -1.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  SDOT, SNRM2
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SCHDD
      INTEGER LDR,P,LDZ,NZ,INFO
      REAL R(LDR,*),X(*),Z(LDZ,*),Y(*),S(*)
      REAL RHO(*),C(*)
C
      INTEGER I,II,J
      REAL A,ALPHA,AZETA,NORM,SNRM2
      REAL SDOT,T,ZETA,B,XX
C
C     SOLVE THE SYSTEM TRANS(R)*A = X, PLACING THE RESULT
C     IN THE ARRAY S.
C
C***FIRST EXECUTABLE STATEMENT  SCHDD
      INFO = 0
      S(1) = X(1)/R(1,1)
      IF (P .LT. 2) GO TO 20
      DO 10 J = 2, P
         S(J) = X(J) - SDOT(J-1,R(1,J),1,S,1)
         S(J) = S(J)/R(J,J)
   10 CONTINUE
   20 CONTINUE
      NORM = SNRM2(P,S,1)
      IF (NORM .LT. 1.0E0) GO TO 30
         INFO = -1
      GO TO 120
   30 CONTINUE
         ALPHA = SQRT(1.0E0-NORM**2)
C
C        DETERMINE THE TRANSFORMATIONS.
C
         DO 40 II = 1, P
            I = P - II + 1
            SCALE = ALPHA + ABS(S(I))
            A = ALPHA/SCALE
            B = S(I)/SCALE
            NORM = SQRT(A**2+B**2)
            C(I) = A/NORM
            S(I) = B/NORM
            ALPHA = SCALE*NORM
   40    CONTINUE
C
C        APPLY THE TRANSFORMATIONS TO R.
C
         DO 60 J = 1, P
            XX = 0.0E0
            DO 50 II = 1, J
               I = J - II + 1
               T = C(I)*XX + S(I)*R(I,J)
               R(I,J) = C(I)*R(I,J) - S(I)*XX
               XX = T
   50       CONTINUE
   60    CONTINUE
C
C        IF REQUIRED, DOWNDATE Z AND RHO.
C
         IF (NZ .LT. 1) GO TO 110
         DO 100 J = 1, NZ
            ZETA = Y(J)
            DO 70 I = 1, P
               Z(I,J) = (Z(I,J) - S(I)*ZETA)/C(I)
               ZETA = C(I)*ZETA - S(I)*Z(I,J)
   70       CONTINUE
            AZETA = ABS(ZETA)
            IF (AZETA .LE. RHO(J)) GO TO 80
               INFO = 1
               RHO(J) = -1.0E0
            GO TO 90
   80       CONTINUE
               RHO(J) = RHO(J)*SQRT(1.0E0-(AZETA/RHO(J))**2)
   90       CONTINUE
  100    CONTINUE
  110    CONTINUE
  120 CONTINUE
      RETURN
      END
