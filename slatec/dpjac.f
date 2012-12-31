*DECK DPJAC
      SUBROUTINE DPJAC (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM, DF,
     +   DJAC, RPAR, IPAR)
C***BEGIN PROLOGUE  DPJAC
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (PJAC-S, DPJAC-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   DPJAC sets up the iteration matrix (involving the Jacobian) for the
C   integration package DDEBDF.
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  DGBFA, DGEFA, DVNRMS
C***COMMON BLOCKS    DDEBD1
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C   920422  Changed DIMENSION statement.  (WRB)
C***END PROLOGUE  DPJAC
C
      INTEGER I, I1, I2, IER, II, IOWND, IOWNS, IPAR, IWM, J, J1,
     1      JJ, JSTART, KFLAG, L, LENP, MAXORD, MBA, MBAND,
     2      MEB1, MEBAND, METH, MITER, ML, ML3, MU, N, NEQ,
     3      NFE, NJE, NQ, NQU, NST, NYH
      DOUBLE PRECISION CON, DI, DVNRMS, EL0, EWT,
     1      FAC, FTEM, H, HL0, HMIN, HMXI, HU, R, R0, ROWND, ROWNS,
     2      RPAR, SAVF, SRUR, TN, UROUND, WM, Y, YH, YI, YJ, YJJ
      EXTERNAL DF, DJAC
      DIMENSION Y(*),YH(NYH,*),EWT(*),FTEM(*),SAVF(*),WM(*),IWM(*),
     1          RPAR(*),IPAR(*)
      COMMON /DDEBD1/ ROWND,ROWNS(210),EL0,H,HMIN,HMXI,HU,TN,UROUND,
     1                IOWND(14),IOWNS(6),IER,JSTART,KFLAG,L,METH,MITER,
     2                MAXORD,N,NQ,NST,NFE,NJE,NQU
C     ------------------------------------------------------------------
C      DPJAC IS CALLED BY DSTOD  TO COMPUTE AND PROCESS THE MATRIX
C      P = I - H*EL(1)*J , WHERE J IS AN APPROXIMATION TO THE JACOBIAN.
C      HERE J IS COMPUTED BY THE USER-SUPPLIED ROUTINE DJAC IF
C      MITER = 1 OR 4, OR BY FINITE DIFFERENCING IF MITER = 2, 3, OR 5.
C      IF MITER = 3, A DIAGONAL APPROXIMATION TO J IS USED.
C      J IS STORED IN WM AND REPLACED BY P.  IF MITER .NE. 3, P IS THEN
C      SUBJECTED TO LU DECOMPOSITION IN PREPARATION FOR LATER SOLUTION
C      OF LINEAR SYSTEMS WITH P AS COEFFICIENT MATRIX. THIS IS DONE
C      BY DGEFA IF MITER = 1 OR 2, AND BY DGBFA IF MITER = 4 OR 5.
C
C      IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION
C      WITH DPJAC USES THE FOLLOWING..
C      Y    = ARRAY CONTAINING PREDICTED VALUES ON ENTRY.
C      FTEM = WORK ARRAY OF LENGTH N (ACOR IN DSTOD ).
C      SAVF = ARRAY CONTAINING DF EVALUATED AT PREDICTED Y.
C      WM   = DOUBLE PRECISION WORK SPACE FOR MATRICES.  ON OUTPUT IT
C      CONTAINS THE
C             INVERSE DIAGONAL MATRIX IF MITER = 3 AND THE LU
C             DECOMPOSITION OF P IF MITER IS 1, 2 , 4, OR 5.
C             STORAGE OF MATRIX ELEMENTS STARTS AT WM(3).
C             WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA..
C             WM(1) = SQRT(UROUND), USED IN NUMERICAL JACOBIAN
C             INCREMENTS.  WM(2) = H*EL0, SAVED FOR LATER USE IF MITER =
C             3.
C      IWM  = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING
C             AT IWM(21), IF MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS
C             THE BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) IF MITER
C             IS 4 OR 5.
C      EL0  = EL(1) (INPUT).
C      IER  = OUTPUT ERROR FLAG,  = 0 IF NO TROUBLE, .NE. 0 IF
C             P MATRIX FOUND TO BE SINGULAR.
C      THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, TN, UROUND,
C      MITER, N, NFE, AND NJE.
C-----------------------------------------------------------------------
C     BEGIN BLOCK PERMITTING ...EXITS TO 240
C        BEGIN BLOCK PERMITTING ...EXITS TO 220
C           BEGIN BLOCK PERMITTING ...EXITS TO 130
C              BEGIN BLOCK PERMITTING ...EXITS TO 70
C***FIRST EXECUTABLE STATEMENT  DPJAC
                  NJE = NJE + 1
                  HL0 = H*EL0
                  GO TO (10,40,90,140,170), MITER
C                 IF MITER = 1, CALL DJAC AND MULTIPLY BY SCALAR.
C                 -----------------------
   10             CONTINUE
                  LENP = N*N
                  DO 20 I = 1, LENP
                     WM(I+2) = 0.0D0
   20             CONTINUE
                  CALL DJAC(TN,Y,WM(3),N,RPAR,IPAR)
                  CON = -HL0
                  DO 30 I = 1, LENP
                     WM(I+2) = WM(I+2)*CON
   30             CONTINUE
C              ...EXIT
                  GO TO 70
C                 IF MITER = 2, MAKE N CALLS TO DF TO APPROXIMATE J.
C                 --------------------
   40             CONTINUE
                  FAC = DVNRMS(N,SAVF,EWT)
                  R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
                  IF (R0 .EQ. 0.0D0) R0 = 1.0D0
                  SRUR = WM(1)
                  J1 = 2
                  DO 60 J = 1, N
                     YJ = Y(J)
                     R = MAX(SRUR*ABS(YJ),R0*EWT(J))
                     Y(J) = Y(J) + R
                     FAC = -HL0/R
                     CALL DF(TN,Y,FTEM,RPAR,IPAR)
                     DO 50 I = 1, N
                        WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
   50                CONTINUE
                     Y(J) = YJ
                     J1 = J1 + N
   60             CONTINUE
                  NFE = NFE + N
   70          CONTINUE
C              ADD IDENTITY MATRIX.
C              -------------------------------------------------
               J = 3
               DO 80 I = 1, N
                  WM(J) = WM(J) + 1.0D0
                  J = J + (N + 1)
   80          CONTINUE
C              DO LU DECOMPOSITION ON P.
C              --------------------------------------------
               CALL DGEFA(WM(3),N,N,IWM(21),IER)
C     .........EXIT
               GO TO 240
C              IF MITER = 3, CONSTRUCT A DIAGONAL APPROXIMATION TO J AND
C              P. ---------
   90          CONTINUE
               WM(2) = HL0
               IER = 0
               R = EL0*0.1D0
               DO 100 I = 1, N
                  Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
  100          CONTINUE
               CALL DF(TN,Y,WM(3),RPAR,IPAR)
               NFE = NFE + 1
               DO 120 I = 1, N
                  R0 = H*SAVF(I) - YH(I,2)
                  DI = 0.1D0*R0 - H*(WM(I+2) - SAVF(I))
                  WM(I+2) = 1.0D0
                  IF (ABS(R0) .LT. UROUND*EWT(I)) GO TO 110
C           .........EXIT
                     IF (ABS(DI) .EQ. 0.0D0) GO TO 130
                     WM(I+2) = 0.1D0*R0/DI
  110             CONTINUE
  120          CONTINUE
C     .........EXIT
               GO TO 240
  130       CONTINUE
            IER = -1
C     ......EXIT
            GO TO 240
C           IF MITER = 4, CALL DJAC AND MULTIPLY BY SCALAR.
C           -----------------------
  140       CONTINUE
            ML = IWM(1)
            MU = IWM(2)
            ML3 = 3
            MBAND = ML + MU + 1
            MEBAND = MBAND + ML
            LENP = MEBAND*N
            DO 150 I = 1, LENP
               WM(I+2) = 0.0D0
  150       CONTINUE
            CALL DJAC(TN,Y,WM(ML3),MEBAND,RPAR,IPAR)
            CON = -HL0
            DO 160 I = 1, LENP
               WM(I+2) = WM(I+2)*CON
  160       CONTINUE
C        ...EXIT
            GO TO 220
C           IF MITER = 5, MAKE MBAND CALLS TO DF TO APPROXIMATE J.
C           ----------------
  170       CONTINUE
            ML = IWM(1)
            MU = IWM(2)
            MBAND = ML + MU + 1
            MBA = MIN(MBAND,N)
            MEBAND = MBAND + ML
            MEB1 = MEBAND - 1
            SRUR = WM(1)
            FAC = DVNRMS(N,SAVF,EWT)
            R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
            IF (R0 .EQ. 0.0D0) R0 = 1.0D0
            DO 210 J = 1, MBA
               DO 180 I = J, N, MBAND
                  YI = Y(I)
                  R = MAX(SRUR*ABS(YI),R0*EWT(I))
                  Y(I) = Y(I) + R
  180          CONTINUE
               CALL DF(TN,Y,FTEM,RPAR,IPAR)
               DO 200 JJ = J, N, MBAND
                  Y(JJ) = YH(JJ,1)
                  YJJ = Y(JJ)
                  R = MAX(SRUR*ABS(YJJ),R0*EWT(JJ))
                  FAC = -HL0/R
                  I1 = MAX(JJ-MU,1)
                  I2 = MIN(JJ+ML,N)
                  II = JJ*MEB1 - ML + 2
                  DO 190 I = I1, I2
                     WM(II+I) = (FTEM(I) - SAVF(I))*FAC
  190             CONTINUE
  200          CONTINUE
  210       CONTINUE
            NFE = NFE + MBA
  220    CONTINUE
C        ADD IDENTITY MATRIX.
C        -------------------------------------------------
         II = MBAND + 2
         DO 230 I = 1, N
            WM(II) = WM(II) + 1.0D0
            II = II + MEBAND
  230    CONTINUE
C        DO LU DECOMPOSITION OF P.
C        --------------------------------------------
         CALL DGBFA(WM(3),MEBAND,N,ML,MU,IWM(21),IER)
  240 CONTINUE
      RETURN
C     ----------------------- END OF SUBROUTINE DPJAC
C     -----------------------
      END
