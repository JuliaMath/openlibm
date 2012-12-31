*DECK D1MERG
      SUBROUTINE D1MERG (TCOS, I1, M1, I2, M2, I3)
C***BEGIN PROLOGUE  D1MERG
C***SUBSIDIARY
C***PURPOSE  Merge two strings of ascending double precision numbers.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (S1MERG-S, D1MERG-D, CMERGE-C, I1MERG-I)
C***AUTHOR  Boland, W. Robert, (LANL)
C           Clemens, Reginald, (PLK)
C***DESCRIPTION
C
C   This subroutine merges two ascending strings of numbers in the
C   array TCOS.  The first string is of length M1 and starts at
C   TCOS(I1+1).  The second string is of length M2 and starts at
C   TCOS(I2+1).  The merged string goes into TCOS(I3+1).
C
C   This routine is currently unused, but was added to complete
C   the set of routines S1MERG and C1MERG (both of which are used).
C
C***ROUTINES CALLED  DCOPY
C***REVISION HISTORY  (YYMMDD)
C   910819  DATE WRITTEN
C***END PROLOGUE  D1MERG
      INTEGER I1, I2, I3, M1, M2
      DOUBLE PRECISION TCOS(*)
C
      INTEGER J1, J2, J3
C
C***FIRST EXECUTABLE STATEMENT  D1MERG
      IF (M1.EQ.0 .AND. M2.EQ.0) RETURN
C
      IF (M1.EQ.0 .AND. M2.NE.0) THEN
         CALL DCOPY (M2, TCOS(I2+1), 1, TCOS(I3+1), 1)
         RETURN
      ENDIF
C
      IF (M1.NE.0 .AND. M2.EQ.0) THEN
         CALL DCOPY (M1, TCOS(I1+1), 1, TCOS(I3+1), 1)
         RETURN
      ENDIF
C
      J1 = 1
      J2 = 1
      J3 = 1
C
   10 IF (TCOS(I1+J1) .LE. TCOS(I2+J2)) THEN
         TCOS(I3+J3) = TCOS(I1+J1)
         J1 = J1+1
         IF (J1 .GT. M1) THEN
            CALL DCOPY (M2-J2+1, TCOS(I2+J2), 1, TCOS(I3+J3+1), 1)
            RETURN
         ENDIF
      ELSE
         TCOS(I3+J3) = TCOS(I2+J2)
         J2 = J2+1
         IF (J2 .GT. M2) THEN
            CALL DCOPY (M1-J1+1, TCOS(I1+J1), 1, TCOS(I3+J3+1), 1)
            RETURN
         ENDIF
      ENDIF
      J3 = J3+1
      GO TO 10
      END
