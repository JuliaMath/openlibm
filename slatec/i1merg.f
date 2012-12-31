*DECK I1MERG
      SUBROUTINE I1MERG (ICOS, I1, M1, I2, M2, I3)
C***BEGIN PROLOGUE  I1MERG
C***SUBSIDIARY
C***PURPOSE  Merge two strings of ascending integers.
C***LIBRARY   SLATEC
C***TYPE      INTEGER (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I)
C***AUTHOR  Boland, W. Robert, (LANL)
C           Clemens, Reginald, (PLK)
C***DESCRIPTION
C
C   This subroutine merges two ascending strings of integers in the
C   array ICOS.  The first string is of length M1 and starts at
C   ICOS(I1+1).  The second string is of length M2 and starts at
C   ICOS(I2+1).  The merged string goes into ICOS(I3+1).
C
C***ROUTINES CALLED  ICOPY
C***REVISION HISTORY  (YYMMDD)
C   920202  DATE WRITTEN
C***END PROLOGUE  I1MERG
      INTEGER I1, I2, I3, M1, M2
      REAL ICOS(*)
C
      INTEGER J1, J2, J3
C
C***FIRST EXECUTABLE STATEMENT  I1MERG
      IF (M1.EQ.0 .AND. M2.EQ.0) RETURN
C
      IF (M1.EQ.0 .AND. M2.NE.0) THEN
         CALL ICOPY (M2, ICOS(I2+1), 1, ICOS(I3+1), 1)
         RETURN
      ENDIF
C
      IF (M1.NE.0 .AND. M2.EQ.0) THEN
         CALL ICOPY (M1, ICOS(I1+1), 1, ICOS(I3+1), 1)
         RETURN
      ENDIF
C
      J1 = 1
      J2 = 1
      J3 = 1
C
   10 IF (ICOS(I1+J1) .LE. ICOS(I2+J2)) THEN
         ICOS(I3+J3) = ICOS(I1+J1)
         J1 = J1+1
         IF (J1 .GT. M1) THEN
            CALL ICOPY (M2-J2+1, ICOS(I2+J2), 1, ICOS(I3+J3+1), 1)
            RETURN
         ENDIF
      ELSE
         ICOS(I3+J3) = ICOS(I2+J2)
         J2 = J2+1
         IF (J2 .GT. M2) THEN
            CALL ICOPY (M1-J1+1, ICOS(I1+J1), 1, ICOS(I3+J3+1), 1)
            RETURN
         ENDIF
      ENDIF
      J3 = J3+1
      GO TO 10
      END
