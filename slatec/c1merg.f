*DECK C1MERG
      SUBROUTINE C1MERG (TCOS, I1, M1, I2, M2, I3)
C***BEGIN PROLOGUE  C1MERG
C***SUBSIDIARY
C***PURPOSE  Merge two strings of complex numbers.  Each string is
C            ascending by the real part.
C***LIBRARY   SLATEC
C***TYPE      COMPLEX (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   This subroutine merges two ascending strings of numbers in the
C   array TCOS.  The first string is of length M1 and starts at
C   TCOS(I1+1).  The second string is of length M2 and starts at
C   TCOS(I2+1).  The merged string goes into TCOS(I3+1).  The ordering
C   is on the real part.
C
C***SEE ALSO  CMGNBN
C***ROUTINES CALLED  CCOPY
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   910408  Modified to use IF-THEN-ELSE.  Make it look like MERGE
C           which was modified earlier due to compiler problems on
C           the IBM RS6000.  (RWC)
C   920130  Code name changed from CMPMRG to C1MERG.  (WRB)
C***END PROLOGUE  C1MERG
      INTEGER I1, I2, I3, M1, M2
      COMPLEX TCOS(*)
C
      INTEGER J1, J2, J3
C
C***FIRST EXECUTABLE STATEMENT  C1MERG
      IF (M1.EQ.0 .AND. M2.EQ.0) RETURN
C
      IF (M1.EQ.0 .AND. M2.NE.0) THEN
         CALL CCOPY (M2, TCOS(I2+1), 1, TCOS(I3+1), 1)
         RETURN
      ENDIF
C
      IF (M1.NE.0 .AND. M2.EQ.0) THEN
         CALL CCOPY (M1, TCOS(I1+1), 1, TCOS(I3+1), 1)
         RETURN
      ENDIF
C
      J1 = 1
      J2 = 1
      J3 = 1
C
   10 IF (REAL(TCOS(J1+I1)) .LE. REAL(TCOS(I2+J2))) THEN
         TCOS(I3+J3) = TCOS(I1+J1)
         J1 = J1+1
         IF (J1 .GT. M1) THEN
            CALL CCOPY (M2-J2+1, TCOS(I2+J2), 1, TCOS(I3+J3+1), 1)
            RETURN
         ENDIF
      ELSE
         TCOS(I3+J3) = TCOS(I2+J2)
         J2 = J2+1
         IF (J2 .GT. M2) THEN
            CALL CCOPY (M1-J1+1, TCOS(I1+J1), 1, TCOS(I3+J3+1), 1)
            RETURN
         ENDIF
      ENDIF
      J3 = J3+1
      GO TO 10
      END
