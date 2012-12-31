*DECK COSGEN
      SUBROUTINE COSGEN (N, IJUMP, FNUM, FDEN, A)
C***BEGIN PROLOGUE  COSGEN
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (COSGEN-S, CMPCSG-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine computes required cosine values in ascending
C     order.  When IJUMP .GT. 1 the routine computes values
C
C        2*COS(J*PI/L) , J=1,2,...,L and J .NE. 0(MOD N/IJUMP+1)
C
C     where L = IJUMP*(N/IJUMP+1).
C
C
C     when IJUMP = 1 it computes
C
C            2*COS((J-FNUM)*PI/(N+FDEN)) ,  J=1, 2, ... ,N
C
C     where
C        FNUM = 0.5, FDEN = 0.0, for regular reduction values.
C        FNUM = 0.0, FDEN = 1.0, for B-R and C-R when ISTAG = 1
C        FNUM = 0.0, FDEN = 0.5, for B-R and C-R when ISTAG = 2
C        FNUM = 0.5, FDEN = 0.5, for B-R and C-R when ISTAG = 2
C                                in POISN2 only.
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  PIMACH
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  COSGEN
      DIMENSION       A(*)
C
C
C***FIRST EXECUTABLE STATEMENT  COSGEN
      PI = PIMACH(DUM)
      IF (N .EQ. 0) GO TO 105
      IF (IJUMP .EQ. 1) GO TO 103
      K3 = N/IJUMP+1
      K4 = K3-1
      PIBYN = PI/(N+IJUMP)
      DO 102 K=1,IJUMP
         K1 = (K-1)*K3
         K5 = (K-1)*K4
         DO 101 I=1,K4
            X = K1+I
            K2 = K5+I
            A(K2) = -2.*COS(X*PIBYN)
  101    CONTINUE
  102 CONTINUE
      GO TO 105
  103 CONTINUE
      NP1 = N+1
      Y = PI/(N+FDEN)
      DO 104 I=1,N
         X = NP1-I-FNUM
         A(I) = 2.*COS(X*Y)
  104 CONTINUE
  105 CONTINUE
      RETURN
      END
