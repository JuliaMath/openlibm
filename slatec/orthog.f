*DECK ORTHOG
      SUBROUTINE ORTHOG (USOL, IDMN, ZN, ZM, PERTRB)
C***BEGIN PROLOGUE  ORTHOG
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPELI
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (ORTHOG-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine orthogonalizes the array USOL with respect to
C     the constant array in a weighted least squares norm.
C
C***SEE ALSO  SEPELI
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    SPLPCM
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  ORTHOG
C
      COMMON /SPLPCM/ KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      DIMENSION       USOL(IDMN,*)           ,ZN(*)      ,ZM(*)
C***FIRST EXECUTABLE STATEMENT  ORTHOG
      ISTR = IS
      IFNL = MS
      JSTR = JS
      JFNL = NS
C
C     COMPUTE WEIGHTED INNER PRODUCTS
C
      UTE = 0.0
      ETE = 0.0
      DO  20 I=IS,MS
         II = I-IS+1
         DO  10 J=JS,NS
            JJ = J-JS+1
            ETE = ETE+ZM(II)*ZN(JJ)
            UTE = UTE+USOL(I,J)*ZM(II)*ZN(JJ)
   10    CONTINUE
   20 CONTINUE
C
C     SET PERTURBATION PARAMETER
C
      PERTRB = UTE/ETE
C
C     SUBTRACT OFF CONSTANT PERTRB
C
      DO  40 I=ISTR,IFNL
         DO  30 J=JSTR,JFNL
            USOL(I,J) = USOL(I,J)-PERTRB
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
