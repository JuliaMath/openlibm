*DECK DEFEHL
      SUBROUTINE DEFEHL (F, NEQ, T, Y, H, YP, F1, F2, F3, F4, F5, YS,
     +   RPAR, IPAR)
C***BEGIN PROLOGUE  DEFEHL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DERKF
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (DEFEHL-S, DFEHL-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     Fehlberg Fourth-Fifth order Runge-Kutta Method
C **********************************************************************
C
C    DEFEHL integrates a system of NEQ first order
C    ordinary differential equations of the form
C               dU/DX = F(X,U)
C    over one step when the vector Y(*) of initial values for U(*) and
C    the vector YP(*) of initial derivatives, satisfying  YP = F(T,Y),
C    are given at the starting point X=T.
C
C    DEFEHL advances the solution over the fixed step H and returns
C    the fifth order (sixth order accurate locally) solution
C    approximation at T+H in the array YS(*).
C    F1,---,F5 are arrays of dimension NEQ which are needed
C    for internal storage.
C    The formulas have been grouped to control loss of significance.
C    DEFEHL should be called with an H not smaller than 13 units of
C    roundoff in T so that the various independent arguments can be
C    distinguished.
C
C    This subroutine has been written with all variables and statement
C    numbers entirely compatible with DERKFS. For greater efficiency,
C    the call to DEFEHL can be replaced by the module beginning with
C    line 222 and extending to the last line just before the return
C    statement.
C
C **********************************************************************
C
C***SEE ALSO  DERKF
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800501  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced statement label.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DEFEHL
C
C
      DIMENSION Y(*),YP(*),F1(*),F2(*),F3(*),F4(*),F5(*),
     1          YS(*),RPAR(*),IPAR(*)
C
C***FIRST EXECUTABLE STATEMENT  DEFEHL
      CH=H/4.
      DO 230 K=1,NEQ
  230   YS(K)=Y(K)+CH*YP(K)
      CALL F(T+CH,YS,F1,RPAR,IPAR)
C
      CH=3.*H/32.
      DO 240 K=1,NEQ
  240   YS(K)=Y(K)+CH*(YP(K)+3.*F1(K))
      CALL F(T+3.*H/8.,YS,F2,RPAR,IPAR)
C
      CH=H/2197.
      DO 250 K=1,NEQ
  250   YS(K)=Y(K)+CH*(1932.*YP(K)+(7296.*F2(K)-7200.*F1(K)))
      CALL F(T+12.*H/13.,YS,F3,RPAR,IPAR)
C
      CH=H/4104.
      DO 260 K=1,NEQ
  260   YS(K)=Y(K)+CH*((8341.*YP(K)-845.*F3(K))+
     1                            (29440.*F2(K)-32832.*F1(K)))
      CALL F(T+H,YS,F4,RPAR,IPAR)
C
      CH=H/20520.
      DO 270 K=1,NEQ
  270   YS(K)=Y(K)+CH*((-6080.*YP(K)+(9295.*F3(K)-5643.*F4(K)))+
     1                             (41040.*F1(K)-28352.*F2(K)))
      CALL F(T+H/2.,YS,F5,RPAR,IPAR)
C
C     COMPUTE APPROXIMATE SOLUTION AT T+H
C
      CH=H/7618050.
      DO 290 K=1,NEQ
  290   YS(K)=Y(K)+CH*((902880.*YP(K)+(3855735.*F3(K)-1371249.*F4(K)))+
     1                (3953664.*F2(K)+277020.*F5(K)))
C
      RETURN
      END
