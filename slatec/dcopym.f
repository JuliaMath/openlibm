*DECK DCOPYM
      SUBROUTINE DCOPYM (N, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DCOPYM
C***PURPOSE  Copy the negative of a vector to a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A5
C***TYPE      DOUBLE PRECISION (SCOPYM-S, DCOPYM-D)
C***KEYWORDS  BLAS, COPY, VECTOR
C***AUTHOR  Kahaner, D. K., (NBS)
C***DESCRIPTION
C
C       Description of Parameters
C           The * Flags Output Variables
C
C       N   Number of elements in vector(s)
C      DX   Double precision vector with N elements
C    INCX   Storage spacing between elements of DX
C      DY*  Double precision negative copy of DX
C    INCY   Storage spacing between elements of DY
C
C      ***  Note that DY = -DX  ***
C
C     Copy negative of d.p. DX to d.p. DY.  For I=0 to N-1,
C     copy  -DX(LX+I*INCX) to DY(LY+I*INCY), where LX=1 if
C     INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is defined
C     in a similar way using INCY.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C***END PROLOGUE  DCOPYM
      DOUBLE PRECISION DX(*), DY(*)
C***FIRST EXECUTABLE STATEMENT  DCOPYM
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
   5  IX=1
      IY=1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = -DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 7.
C
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = -DX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = -DX(I)
        DY(I+1) = -DX(I+1)
        DY(I+2) = -DX(I+2)
        DY(I+3) = -DX(I+3)
        DY(I+4) = -DX(I+4)
        DY(I+5) = -DX(I+5)
        DY(I+6) = -DX(I+6)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = -DX(I)
   70 CONTINUE
      RETURN
      END
