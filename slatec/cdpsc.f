*DECK CDPSC
      SUBROUTINE CDPSC (KSGN, N, NQ, YH)
C***BEGIN PROLOGUE  CDPSC
C***SUBSIDIARY
C***PURPOSE  Subroutine CDPSC computes the predicted YH values by
C            effectively multiplying the YH array by the Pascal triangle
C            matrix when KSGN is +1, and performs the inverse function
C            when KSGN is -1.
C***LIBRARY   SLATEC (SDRIVE)
C***TYPE      COMPLEX (SDPSC-S, DDPSC-D, CDPSC-C)
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             Mail Stop D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  CDPSC
      INTEGER I, J, J1, J2, KSGN, N, NQ
      COMPLEX YH(N,*)
C***FIRST EXECUTABLE STATEMENT  CDPSC
      IF (KSGN .GT. 0) THEN
        DO 10 J1 = 1,NQ
          DO 10 J2 = J1,NQ
            J = NQ - J2 + J1
            DO 10 I = 1,N
 10           YH(I,J) = YH(I,J) + YH(I,J+1)
      ELSE
        DO 30 J1 = 1,NQ
          DO 30 J2 = J1,NQ
            J = NQ - J2 + J1
            DO 30 I = 1,N
 30           YH(I,J) = YH(I,J) - YH(I,J+1)
      END IF
      RETURN
      END
