*DECK DXC210
      SUBROUTINE DXC210 (K, Z, J, IERROR)
C***BEGIN PROLOGUE  DXC210
C***PURPOSE  To provide double-precision floating-point arithmetic
C            with an extended exponent range.
C***LIBRARY   SLATEC
C***CATEGORY  A3D
C***TYPE      DOUBLE PRECISION (XC210-S, DXC210-D)
C***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
C***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
C           Smith, John M., (NBS and George Mason University)
C***DESCRIPTION
C     INTEGER K, J
C     DOUBLE PRECISION Z
C
C                  GIVEN K THIS SUBROUTINE COMPUTES J AND Z
C                  SUCH THAT  RADIX**K = Z*10**J, WHERE Z IS IN
C                  THE RANGE 1/10 .LE. Z .LT. 1.
C                  THE VALUE OF Z WILL BE ACCURATE TO FULL
C                  DOUBLE-PRECISION PROVIDED THE NUMBER
C                  OF DECIMAL PLACES IN THE LARGEST
C                  INTEGER PLUS THE NUMBER OF DECIMAL
C                  PLACES CARRIED IN DOUBLE-PRECISION DOES NOT
C                  EXCEED 60. DXC210 IS CALLED BY SUBROUTINE
C                  DXCON WHEN NECESSARY. THE USER SHOULD
C                  NEVER NEED TO CALL DXC210 DIRECTLY.
C
C***SEE ALSO  DXSET
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***COMMON BLOCKS    DXBLK3
C***REVISION HISTORY  (YYMMDD)
C   820712  DATE WRITTEN
C   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Changed all specific intrinsics to generic.  (WRB)
C           Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  DXC210
      DOUBLE PRECISION Z
      INTEGER K, J
      INTEGER NLG102, MLG102, LG102
      COMMON /DXBLK3/ NLG102, MLG102, LG102(21)
      SAVE /DXBLK3/
C
C   THE CONDITIONS IMPOSED ON NLG102, MLG102, AND LG102 BY
C THIS SUBROUTINE ARE
C
C     (1) NLG102 .GE. 2
C
C     (2) MLG102 .GE. 1
C
C     (3) 2*MLG102*(MLG102 - 1) .LE. 2**NBITS - 1
C
C THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING
C IN SUBROUTINE DXSET.
C
C***FIRST EXECUTABLE STATEMENT  DXC210
      IERROR=0
      IF (K.EQ.0) GO TO 70
      M = MLG102
      KA = ABS(K)
      KA1 = KA/M
      KA2 = MOD(KA,M)
      IF (KA1.GE.M) GO TO 60
      NM1 = NLG102 - 1
      NP1 = NLG102 + 1
      IT = KA2*LG102(NP1)
      IC = IT/M
      ID = MOD(IT,M)
      Z = ID
      IF (KA1.GT.0) GO TO 20
      DO 10 II=1,NM1
        I = NP1 - II
        IT = KA2*LG102(I) + IC
        IC = IT/M
        ID = MOD(IT,M)
        Z = Z/M + ID
   10 CONTINUE
      JA = KA*LG102(1) + IC
      GO TO 40
   20 CONTINUE
      DO 30 II=1,NM1
        I = NP1 - II
        IT = KA2*LG102(I) + KA1*LG102(I+1) + IC
        IC = IT/M
        ID = MOD(IT,M)
        Z = Z/M + ID
   30 CONTINUE
      JA = KA*LG102(1) + KA1*LG102(2) + IC
   40 CONTINUE
      Z = Z/M
      IF (K.GT.0) GO TO 50
      J = -JA
      Z = 10.0D0**(-Z)
      GO TO 80
   50 CONTINUE
      J = JA + 1
      Z = 10.0D0**(Z-1.0D0)
      GO TO 80
   60 CONTINUE
C   THIS ERROR OCCURS IF K EXCEEDS  MLG102**2 - 1  IN MAGNITUDE.
C
      CALL XERMSG ('SLATEC', 'DXC210', 'K too large', 208, 1)
      IERROR=208
      RETURN
   70 CONTINUE
      J = 0
      Z = 1.0D0
   80 RETURN
      END
