*DECK XADD
      SUBROUTINE XADD (X, IX, Y, IY, Z, IZ, IERROR)
C***BEGIN PROLOGUE  XADD
C***PURPOSE  To provide single-precision floating-point arithmetic
C            with an extended exponent range.
C***LIBRARY   SLATEC
C***CATEGORY  A3D
C***TYPE      SINGLE PRECISION (XADD-S, DXADD-D)
C***KEYWORDS  EXTENDED-RANGE SINGLE-PRECISION ARITHMETIC
C***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
C           Smith, John M., (NBS and George Mason University)
C***DESCRIPTION
C     REAL X, Y, Z
C     INTEGER IX, IY, IZ
C
C                  FORMS THE EXTENDED-RANGE SUM  (Z,IZ) =
C                  (X,IX) + (Y,IY).  (Z,IZ) IS ADJUSTED
C                  BEFORE RETURNING. THE INPUT OPERANDS
C                  NEED NOT BE IN ADJUSTED FORM, BUT THEIR
C                  PRINCIPAL PARTS MUST SATISFY
C                  RADIX**(-2L).LE.ABS(X).LE.RADIX**(2L),
C                  RADIX**(-2L).LE.ABS(Y).LE.RADIX**(2L).
C
C***SEE ALSO  XSET
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XADJ
C***COMMON BLOCKS    XBLK2
C***REVISION HISTORY  (YYMMDD)
C   820712  DATE WRITTEN
C   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Changed all specific intrinsics to generic.  (WRB)
C           Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  XADD
      REAL X, Y, Z
      INTEGER IX, IY, IZ
      REAL RADIX, RADIXL, RAD2L, DLG10R
      INTEGER L, L2, KMAX
      COMMON /XBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
      SAVE /XBLK2/
C
C
C   THE CONDITIONS IMPOSED ON L AND KMAX BY THIS SUBROUTINE
C ARE
C     (1) 1 .LT. L .LE. 0.5*LOGR(0.5*DZERO)
C
C     (2) NRADPL .LT. L .LE. KMAX/6
C
C     (3) KMAX .LE. (2**NBITS - 4*L - 1)/2
C
C THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING
C IN SUBROUTINE XSET.
C
C***FIRST EXECUTABLE STATEMENT  XADD
      IERROR=0
      IF (X.NE.0.0) GO TO 10
      Z = Y
      IZ = IY
      GO TO 220
   10 IF (Y.NE.0.0) GO TO 20
      Z = X
      IZ = IX
      GO TO 220
   20 CONTINUE
      IF (IX.GE.0 .AND. IY.GE.0) GO TO 40
      IF (IX.LT.0 .AND. IY.LT.0) GO TO 40
      IF (ABS(IX).LE.6*L .AND. ABS(IY).LE.6*L) GO TO 40
      IF (IX.GE.0) GO TO 30
      Z = Y
      IZ = IY
      GO TO 220
   30 CONTINUE
      Z = X
      IZ = IX
      GO TO 220
   40 I = IX - IY
      IF (I) 80, 50, 90
   50 IF (ABS(X).GT.1.0 .AND. ABS(Y).GT.1.0) GO TO 60
      IF (ABS(X).LT.1.0 .AND. ABS(Y).LT.1.0) GO TO 70
      Z = X + Y
      IZ = IX
      GO TO 220
   60 S = X/RADIXL
      T = Y/RADIXL
      Z = S + T
      IZ = IX + L
      GO TO 220
   70 S = X*RADIXL
      T = Y*RADIXL
      Z = S + T
      IZ = IX - L
      GO TO 220
   80 S = Y
      IS = IY
      T = X
      GO TO 100
   90 S = X
      IS = IX
      T = Y
  100 CONTINUE
C
C  AT THIS POINT, THE ONE OF (X,IX) OR (Y,IY) THAT HAS THE
C LARGER AUXILIARY INDEX IS STORED IN (S,IS). THE PRINCIPAL
C PART OF THE OTHER INPUT IS STORED IN T.
C
      I1 = ABS(I)/L
      I2 = MOD(ABS(I),L)
      IF (ABS(T).GE.RADIXL) GO TO 130
      IF (ABS(T).GE.1.0) GO TO 120
      IF (RADIXL*ABS(T).GE.1.0) GO TO 110
      J = I1 + 1
      T = T*RADIX**(L-I2)
      GO TO 140
  110 J = I1
      T = T*RADIX**(-I2)
      GO TO 140
  120 J = I1 - 1
      IF (J.LT.0) GO TO 110
      T = T*RADIX**(-I2)/RADIXL
      GO TO 140
  130 J = I1 - 2
      IF (J.LT.0) GO TO 120
      T = T*RADIX**(-I2)/RAD2L
  140 CONTINUE
C
C  AT THIS POINT, SOME OR ALL OF THE DIFFERENCE IN THE
C AUXILIARY INDICES HAS BEEN USED TO EFFECT A LEFT SHIFT
C OF T.  THE SHIFTED VALUE OF T SATISFIES
C
C       RADIX**(-2*L) .LE. ABS(T) .LE. 1.0
C
C AND, IF J=0, NO FURTHER SHIFTING REMAINS TO BE DONE.
C
      IF (J.EQ.0) GO TO 190
      IF (ABS(S).GE.RADIXL .OR. J.GT.3) GO TO 150
      IF (ABS(S).GE.1.0) GO TO (180, 150, 150), J
      IF (RADIXL*ABS(S).GE.1.0) GO TO (180, 170, 150), J
      GO TO (180, 170, 160), J
  150 Z = S
      IZ = IS
      GO TO 220
  160 S = S*RADIXL
  170 S = S*RADIXL
  180 S = S*RADIXL
  190 CONTINUE
C
C   AT THIS POINT, THE REMAINING DIFFERENCE IN THE
C AUXILIARY INDICES HAS BEEN USED TO EFFECT A RIGHT SHIFT
C OF S.  IF THE SHIFTED VALUE OF S WOULD HAVE EXCEEDED
C RADIX**L, THEN (S,IS) IS RETURNED AS THE VALUE OF THE
C SUM.
C
      IF (ABS(S).GT.1.0 .AND. ABS(T).GT.1.0) GO TO 200
      IF (ABS(S).LT.1.0 .AND. ABS(T).LT.1.0) GO TO 210
      Z = S + T
      IZ = IS - J*L
      GO TO 220
  200 S = S/RADIXL
      T = T/RADIXL
      Z = S + T
      IZ = IS - J*L + L
      GO TO 220
  210 S = S*RADIXL
      T = T*RADIXL
      Z = S + T
      IZ = IS - J*L - L
  220 CALL XADJ(Z, IZ,IERROR)
      RETURN
      END
