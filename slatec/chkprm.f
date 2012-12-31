*DECK CHKPRM
      SUBROUTINE CHKPRM (INTL, IORDER, A, B, M, MBDCND, C, D, N, NBDCND,
     +   COFX, COFY, IDMN, IERROR)
C***BEGIN PROLOGUE  CHKPRM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPELI
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (CHKPRM-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This program checks the input parameters for errors.
C
C***SEE ALSO  SEPELI
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  CHKPRM
C
      EXTERNAL        COFX       ,COFY
C***FIRST EXECUTABLE STATEMENT  CHKPRM
      IERROR = 1
      IF (A.GE.B .OR. C.GE.D) RETURN
C
C     CHECK BOUNDARY SWITCHES
C
      IERROR = 2
      IF (MBDCND.LT.0 .OR. MBDCND.GT.4) RETURN
      IERROR = 3
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) RETURN
C
C     CHECK FIRST DIMENSION IN CALLING ROUTINE
C
      IERROR = 5
      IF (IDMN .LT. 7) RETURN
C
C     CHECK M
C
      IERROR = 6
      IF (M.GT.(IDMN-1) .OR. M.LT.6) RETURN
C
C     CHECK N
C
      IERROR = 7
      IF (N .LT. 5) RETURN
C
C     CHECK IORDER
C
      IERROR = 8
      IF (IORDER.NE.2 .AND. IORDER.NE.4) RETURN
C
C     CHECK INTL
C
      IERROR = 9
      IF (INTL.NE.0 .AND. INTL.NE.1) RETURN
C
C     CHECK THAT EQUATION IS ELLIPTIC
C
      DLX = (B-A)/M
      DLY = (D-C)/N
      DO  30 I=2,M
         XI = A+(I-1)*DLX
         CALL COFX (XI,AI,BI,CI)
         DO  20 J=2,N
            YJ = C+(J-1)*DLY
            CALL COFY (YJ,DJ,EJ,FJ)
            IF (AI*DJ .GT. 0.0) GO TO  10
            IERROR = 10
            RETURN
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
C
C     NO ERROR FOUND
C
      IERROR = 0
      RETURN
      END
