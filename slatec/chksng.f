*DECK CHKSNG
      SUBROUTINE CHKSNG (MBDCND, NBDCND, ALPHA, BETA, GAMA, XNU, COFX,
     +   COFY, SINGLR)
C***BEGIN PROLOGUE  CHKSNG
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPELI
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (CHKSNG-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine checks if the PDE SEPELI
C     must solve is a singular operator.
C
C***SEE ALSO  SEPELI
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    SPLPCM
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  CHKSNG
C
      COMMON /SPLPCM/ KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      LOGICAL         SINGLR
C***FIRST EXECUTABLE STATEMENT  CHKSNG
      SINGLR = .FALSE.
C
C     CHECK IF THE BOUNDARY CONDITIONS ARE
C     ENTIRELY PERIODIC AND/OR MIXED
C
      IF ((MBDCND.NE.0 .AND. MBDCND.NE.3) .OR.
     1    (NBDCND.NE.0 .AND. NBDCND.NE.3)) RETURN
C
C     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN
C
      IF (MBDCND .NE. 3) GO TO  10
      IF (ALPHA.NE.0.0 .OR. BETA.NE.0.0) RETURN
   10 IF (NBDCND .NE. 3) GO TO  20
      IF (GAMA.NE.0.0 .OR. XNU.NE.0.0) RETURN
   20 CONTINUE
C
C     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
C     ARE ZERO
C
      DO  30 I=IS,MS
         XI = AIT+(I-1)*DLX
         CALL COFX (XI,AI,BI,CI)
         IF (CI .NE. 0.0) RETURN
   30 CONTINUE
      DO  40 J=JS,NS
         YJ = CIT+(J-1)*DLY
         CALL COFY (YJ,DJ,EJ,FJ)
         IF (FJ .NE. 0.0) RETURN
   40 CONTINUE
C
C     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED
C
      SINGLR = .TRUE.
      RETURN
      END
