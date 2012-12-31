*DECK DEFER
      SUBROUTINE DEFER (COFX, COFY, IDMN, USOL, GRHS)
C***BEGIN PROLOGUE  DEFER
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPELI
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (DEFER-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine first approximates the truncation error given by
C     TRUN1(X,Y)=DLX**2*TX+DLY**2*TY where
C     TX=AFUN(X)*UXXXX/12.0+BFUN(X)*UXXX/6.0 on the interior and
C     at the boundaries if periodic (here UXXX,UXXXX are the third
C     and fourth partial derivatives of U with respect to X).
C     TX is of the form AFUN(X)/3.0*(UXXXX/4.0+UXXX/DLX)
C     at X=A or X=B if the boundary condition there is mixed.
C     TX=0.0 along specified boundaries.  TY has symmetric form
C     in Y with X,AFUN(X),BFUN(X) replaced by Y,DFUN(Y),EFUN(Y).
C     The second order solution in USOL is used to approximate
C     (via second order finite differencing) the truncation error
C     and the result is added to the right hand side in GRHS
C     and then transferred to USOL to be used as a new right
C     hand side when calling BLKTRI for a fourth order solution.
C
C***SEE ALSO  SEPELI
C***ROUTINES CALLED  DX, DY
C***COMMON BLOCKS    SPLPCM
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  DEFER
C
      COMMON /SPLPCM/ KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      DIMENSION       GRHS(IDMN,*)           ,USOL(IDMN,*)
      EXTERNAL        COFX       ,COFY
C***FIRST EXECUTABLE STATEMENT  DEFER
      DO  40 J=JS,NS
         YJ = CIT+(J-1)*DLY
         CALL COFY (YJ,DJ,EJ,FJ)
         DO  30 I=IS,MS
            XI = AIT+(I-1)*DLX
            CALL COFX (XI,AI,BI,CI)
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ)
C
            CALL DX (USOL,IDMN,I,J,UXXX,UXXXX)
            CALL DY (USOL,IDMN,I,J,UYYY,UYYYY)
            TX = AI*UXXXX/12.0+BI*UXXX/6.0
            TY = DJ*UYYYY/12.0+EJ*UYYY/6.0
C
C     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC
C
            IF (KSWX.EQ.1 .OR. (I.GT.1 .AND. I.LT.K)) GO TO  10
            TX = AI/3.0*(UXXXX/4.0+UXXX/DLX)
   10       IF (KSWY.EQ.1 .OR. (J.GT.1 .AND. J.LT.L)) GO TO  20
            TY = DJ/3.0*(UYYYY/4.0+UYYY/DLY)
   20       GRHS(I,J) = GRHS(I,J)+DLX**2*TX+DLY**2*TY
   30    CONTINUE
   40 CONTINUE
C
C     RESET THE RIGHT HAND SIDE IN USOL
C
      DO  60 I=IS,MS
         DO  50 J=JS,NS
            USOL(I,J) = GRHS(I,J)
   50    CONTINUE
   60 CONTINUE
      RETURN
      END
