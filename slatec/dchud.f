*DECK DCHUD
      SUBROUTINE DCHUD (R, LDR, P, X, Z, LDZ, NZ, Y, RHO, C, S)
C***BEGIN PROLOGUE  DCHUD
C***PURPOSE  Update an augmented Cholesky decomposition of the
C            triangular part of an augmented QR decomposition.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D7B
C***TYPE      DOUBLE PRECISION (SCHUD-S, DCHUD-D, CCHUD-C)
C***KEYWORDS  CHOLESKY DECOMPOSITION, LINEAR ALGEBRA, LINPACK, MATRIX,
C             UPDATE
C***AUTHOR  Stewart, G. W., (U. of Maryland)
C***DESCRIPTION
C
C     DCHUD updates an augmented Cholesky decomposition of the
C     triangular part of an augmented QR decomposition.  Specifically,
C     given an upper triangular matrix R of order P, a row vector
C     X, a column vector Z, and a scalar Y, DCHUD determines a
C     unitary matrix U and a scalar ZETA such that
C
C
C                              (R  Z)     (RR   ZZ )
C                         U  * (    )  =  (        ) ,
C                              (X  Y)     ( 0  ZETA)
C
C     where RR is upper triangular.  If R and Z have been
C     obtained from the factorization of a least squares
C     problem, then RR and ZZ are the factors corresponding to
C     the problem with the observation (X,Y) appended.  In this
C     case, if RHO is the norm of the residual vector, then the
C     norm of the residual vector of the updated problem is
C     SQRT(RHO**2 + ZETA**2).  DCHUD will simultaneously update
C     several triplets (Z,Y,RHO).
C     For a less terse description of what DCHUD does and how
C     it may be applied, see the LINPACK guide.
C
C     The matrix U is determined as the product U(P)*...*U(1),
C     where U(I) is a rotation in the (I,P+1) plane of the
C     form
C
C                       (     C(I)      S(I) )
C                       (                    ) .
C                       (    -S(I)      C(I) )
C
C     The rotations are chosen so that C(I) is double precision.
C
C     On Entry
C
C         R      DOUBLE PRECISION(LDR,P), where LDR .GE. P.
C                R contains the upper triangular matrix
C                that is to be updated.  The part of R
C                below the diagonal is not referenced.
C
C         LDR    INTEGER.
C                LDR is the leading dimension of the array R.
C
C         P      INTEGER.
C                P is the order of the matrix R.
C
C         X      DOUBLE PRECISION(P).
C                X contains the row to be added to R.  X is
C                not altered by DCHUD.
C
C         Z      DOUBLE PRECISION(LDZ,N)Z), where LDZ .GE. P.
C                Z is an array containing NZ P-vectors to
C                be updated with R.
C
C         LDZ    INTEGER.
C                LDZ is the leading dimension of the array Z.
C
C         NZ     INTEGER.
C                NZ is the number of vectors to be updated
C                NZ may be zero, in which case Z, Y, and RHO
C                are not referenced.
C
C         Y      DOUBLE PRECISION(NZ).
C                Y contains the scalars for updating the vectors
C                Z.  Y is not altered by DCHUD.
C
C         RHO    DOUBLE PRECISION(NZ).
C                RHO contains the norms of the residual
C                vectors that are to be updated.  If RHO(J)
C                is negative, it is left unaltered.
C
C     On Return
C
C         RC
C         RHO    contain the updated quantities.
C         Z
C
C         C      DOUBLE PRECISION(P).
C                C contains the cosines of the transforming
C                rotations.
C
C         S      DOUBLE PRECISION(P).
C                S contains the sines of the transforming
C                rotations.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DROTG
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DCHUD
      INTEGER LDR,P,LDZ,NZ
      DOUBLE PRECISION RHO(*),C(*)
      DOUBLE PRECISION R(LDR,*),X(*),Z(LDZ,*),Y(*),S(*)
C
      INTEGER I,J,JM1
      DOUBLE PRECISION AZETA,SCALE
      DOUBLE PRECISION T,XJ,ZETA
C
C     UPDATE R.
C
C***FIRST EXECUTABLE STATEMENT  DCHUD
      DO 30 J = 1, P
         XJ = X(J)
C
C        APPLY THE PREVIOUS ROTATIONS.
C
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            T = C(I)*R(I,J) + S(I)*XJ
            XJ = C(I)*XJ - S(I)*R(I,J)
            R(I,J) = T
   10    CONTINUE
   20    CONTINUE
C
C        COMPUTE THE NEXT ROTATION.
C
         CALL DROTG(R(J,J),XJ,C(J),S(J))
   30 CONTINUE
C
C     IF REQUIRED, UPDATE Z AND RHO.
C
      IF (NZ .LT. 1) GO TO 70
      DO 60 J = 1, NZ
         ZETA = Y(J)
         DO 40 I = 1, P
            T = C(I)*Z(I,J) + S(I)*ZETA
            ZETA = C(I)*ZETA - S(I)*Z(I,J)
            Z(I,J) = T
   40    CONTINUE
         AZETA = ABS(ZETA)
         IF (AZETA .EQ. 0.0D0 .OR. RHO(J) .LT. 0.0D0) GO TO 50
            SCALE = AZETA + RHO(J)
            RHO(J) = SCALE*SQRT((AZETA/SCALE)**2+(RHO(J)/SCALE)**2)
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      RETURN
      END
