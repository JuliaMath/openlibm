*DECK DFEHL
      SUBROUTINE DFEHL (DF, NEQ, T, Y, H, YP, F1, F2, F3, F4, F5, YS,
     +   RPAR, IPAR)
C***BEGIN PROLOGUE  DFEHL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDERKF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (DEFEHL-S, DFEHL-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     Fehlberg Fourth-Fifth Order Runge-Kutta Method
C **********************************************************************
C
C    DFEHL integrates a system of NEQ first order
C    ordinary differential equations of the form
C               DU/DX = DF(X,U)
C    over one step when the vector Y(*) of initial values for U(*) and
C    the vector YP(*) of initial derivatives, satisfying  YP = DF(T,Y),
C    are given at the starting point X=T.
C
C    DFEHL advances the solution over the fixed step H and returns
C    the fifth order (sixth order accurate locally) solution
C    approximation at T+H in the array YS(*).
C    F1,---,F5 are arrays of dimension NEQ which are needed
C    for internal storage.
C    The formulas have been grouped to control loss of significance.
C    DFEHL should be called with an H not smaller than 13 units of
C    roundoff in T so that the various independent arguments can be
C    distinguished.
C
C    This subroutine has been written with all variables and statement
C    numbers entirely compatible with DRKFS. For greater efficiency,
C    the call to DFEHL can be replaced by the module beginning with
C    line 222 and extending to the last line just before the return
C    statement.
C
C **********************************************************************
C
C***SEE ALSO  DDERKF
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DFEHL
C
      INTEGER IPAR, K, NEQ
      DOUBLE PRECISION CH, F1, F2, F3, F4, F5, H, RPAR, T, Y, YP, YS
      DIMENSION Y(*),YP(*),F1(*),F2(*),F3(*),F4(*),F5(*),
     1          YS(*),RPAR(*),IPAR(*)
C
C***FIRST EXECUTABLE STATEMENT  DFEHL
      CH = H/4.0D0
      DO 10 K = 1, NEQ
         YS(K) = Y(K) + CH*YP(K)
   10 CONTINUE
      CALL DF(T+CH,YS,F1,RPAR,IPAR)
C
      CH = 3.0D0*H/32.0D0
      DO 20 K = 1, NEQ
         YS(K) = Y(K) + CH*(YP(K) + 3.0D0*F1(K))
   20 CONTINUE
      CALL DF(T+3.0D0*H/8.0D0,YS,F2,RPAR,IPAR)
C
      CH = H/2197.0D0
      DO 30 K = 1, NEQ
         YS(K) = Y(K)
     1           + CH
     2             *(1932.0D0*YP(K) + (7296.0D0*F2(K) - 7200.0D0*F1(K)))
   30 CONTINUE
      CALL DF(T+12.0D0*H/13.0D0,YS,F3,RPAR,IPAR)
C
      CH = H/4104.0D0
      DO 40 K = 1, NEQ
         YS(K) = Y(K)
     1           + CH
     2             *((8341.0D0*YP(K) - 845.0D0*F3(K))
     3               + (29440.0D0*F2(K) - 32832.0D0*F1(K)))
   40 CONTINUE
      CALL DF(T+H,YS,F4,RPAR,IPAR)
C
      CH = H/20520.0D0
      DO 50 K = 1, NEQ
         YS(K) = Y(K)
     1           + CH
     2             *((-6080.0D0*YP(K)
     3                + (9295.0D0*F3(K) - 5643.0D0*F4(K)))
     4               + (41040.0D0*F1(K) - 28352.0D0*F2(K)))
   50 CONTINUE
      CALL DF(T+H/2.0D0,YS,F5,RPAR,IPAR)
C
C     COMPUTE APPROXIMATE SOLUTION AT T+H
C
      CH = H/7618050.0D0
      DO 60 K = 1, NEQ
         YS(K) = Y(K)
     1           + CH
     2             *((902880.0D0*YP(K)
     3                + (3855735.0D0*F3(K) - 1371249.0D0*F4(K)))
     4               + (3953664.0D0*F2(K) + 277020.0D0*F5(K)))
   60 CONTINUE
C
      RETURN
      END
