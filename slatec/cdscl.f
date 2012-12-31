*DECK CDSCL
      SUBROUTINE CDSCL (HMAX, N, NQ, RMAX, H, RC, RH, YH)
C***BEGIN PROLOGUE  CDSCL
C***SUBSIDIARY
C***PURPOSE  Subroutine CDSCL rescales the YH array whenever the step
C            size is changed.
C***LIBRARY   SLATEC (SDRIVE)
C***TYPE      COMPLEX (SDSCL-S, DDSCL-D, CDSCL-C)
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
C***END PROLOGUE  CDSCL
      INTEGER I, J, N, NQ
      COMPLEX YH(N,*)
      REAL H, HMAX, RC, RH, RMAX, R1
C***FIRST EXECUTABLE STATEMENT  CDSCL
      IF (H .LT. 1.E0) THEN
        RH = MIN(ABS(H)*RH, ABS(H)*RMAX, HMAX)/ABS(H)
      ELSE
        RH = MIN(RH, RMAX, HMAX/ABS(H))
      END IF
      R1 = 1.E0
      DO 10 J = 1,NQ
        R1 = R1*RH
        DO 10 I = 1,N
 10       YH(I,J+1) = YH(I,J+1)*R1
      H = H*RH
      RC = RC*RH
      RETURN
      END
