*DECK DDPST
      SUBROUTINE DDPST (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML,
     8   MU, N, NDE, NQ, SAVE2, T, USERS, Y, YH, YWT, UROUND, NFE, NJE,
     8   A, DFDY, FAC, IER, IPVT, SAVE1, ISWFLG, BND, JSTATE)
C***BEGIN PROLOGUE  DDPST
C***SUBSIDIARY
C***PURPOSE  Subroutine DDPST evaluates the Jacobian matrix of the right
C            hand side of the differential equations.
C***LIBRARY   SLATEC (SDRIVE)
C***TYPE      DOUBLE PRECISION (SDPST-S, DDPST-D, CDPST-C)
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             Mail Stop D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***DESCRIPTION
C
C  If MITER is 1, 2, 4, or 5, the matrix
C  P = I - L(0)*H*Jacobian is stored in DFDY and subjected to LU
C  decomposition, with the results also stored in DFDY.
C
C***ROUTINES CALLED  DGBFA, DGEFA, DNRM2
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDPST
      INTEGER I, IFLAG, IMAX, IMPL, INFO, ISWFLG, J, J2, JSTATE, K,
     8        MATDIM, MITER, ML, MU, MW, N, NDE, NFE, NJE, NQ
      DOUBLE PRECISION A(MATDIM,*), BL, BND, BP, BR, BU, DFDY(MATDIM,*),
     8     DFDYMX, DIFF, DY, EL(13,12), FAC(*), FACMAX, FACMIN, FACTOR,
     8     H, SAVE1(*), SAVE2(*), SCALE, DNRM2, T, UROUND, Y(*),
     8     YH(N,*), YJ, YS, YWT(*)
      INTEGER IPVT(*)
      LOGICAL IER
      PARAMETER(FACMAX = .5D0, BU = 0.5D0)
C***FIRST EXECUTABLE STATEMENT  DDPST
      NJE = NJE + 1
      IER = .FALSE.
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (MITER .EQ. 1) THEN
          CALL JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
          IF (N .EQ. 0) THEN
            JSTATE = 8
            RETURN
          END IF
          IF (ISWFLG .EQ. 3) BND = DNRM2(N*N, DFDY, 1)
          FACTOR = -EL(1,NQ)*H
          DO 110 J = 1,N
            DO 110 I = 1,N
 110          DFDY(I,J) = FACTOR*DFDY(I,J)
        ELSE IF (MITER .EQ. 2) THEN
          BR = UROUND**(.875D0)
          BL = UROUND**(.75D0)
          BP = UROUND**(-.15D0)
          FACMIN = UROUND**(.78D0)
          DO 170 J = 1,N
            YS = MAX(ABS(YWT(J)), ABS(Y(J)))
 120        DY = FAC(J)*YS
            IF (DY .EQ. 0.D0) THEN
              IF (FAC(J) .LT. FACMAX) THEN
                FAC(J) = MIN(100.D0*FAC(J), FACMAX)
                GO TO 120
              ELSE
                DY = YS
              END IF
            END IF
            IF (NQ .EQ. 1) THEN
              DY = SIGN(DY, SAVE2(J))
            ELSE
              DY = SIGN(DY, YH(J,3))
            END IF
            DY = (Y(J) + DY) - Y(J)
            YJ = Y(J)
            Y(J) = Y(J) + DY
            CALL F (N, T, Y, SAVE1)
            IF (N .EQ. 0) THEN
              JSTATE = 6
              RETURN
            END IF
            Y(J) = YJ
            FACTOR = -EL(1,NQ)*H/DY
            DO 140 I = 1,N
 140          DFDY(I,J) = (SAVE1(I) - SAVE2(I))*FACTOR
C                                                                 Step 1
            DIFF = ABS(SAVE2(1) - SAVE1(1))
            IMAX = 1
            DO 150 I = 2,N
              IF (ABS(SAVE2(I) - SAVE1(I)) .GT. DIFF) THEN
                IMAX = I
                DIFF = ABS(SAVE2(I) - SAVE1(I))
              END IF
 150          CONTINUE
C                                                                 Step 2
            IF (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX))) .GT. 0.D0) THEN
              SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
C                                                                 Step 3
              IF (DIFF .GT. BU*SCALE) THEN
                FAC(J) = MAX(FACMIN, FAC(J)*.5D0)
              ELSE IF (BR*SCALE .LE. DIFF .AND. DIFF .LE. BL*SCALE) THEN
                FAC(J) = MIN(FAC(J)*2.D0, FACMAX)
C                                                                 Step 4
              ELSE IF (DIFF .LT. BR*SCALE) THEN
                FAC(J) = MIN(BP*FAC(J), FACMAX)
              END IF
            END IF
 170        CONTINUE
          IF (ISWFLG .EQ. 3) BND = DNRM2(N*N, DFDY, 1)/(-EL(1,NQ)*H)
          NFE = NFE + N
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 190 I = 1,N
 190        DFDY(I,I) = DFDY(I,I) + 1.D0
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 210 J = 1,N
            DO 210 I = 1,N
 210          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 230 I = 1,NDE
 230        DFDY(I,I) = DFDY(I,I) + A(I,1)
        ELSE IF (IMPL .EQ. 3) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 220 J = 1,NDE
            DO 220 I = 1,NDE
 220          DFDY(I,J) = DFDY(I,J) + A(I,J)
        END IF
        CALL DGEFA (DFDY, MATDIM, N, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (MITER .EQ. 4) THEN
          CALL JACOBN (N, T, Y, DFDY(ML+1,1), MATDIM, ML, MU)
          IF (N .EQ. 0) THEN
            JSTATE = 8
            RETURN
          END IF
          FACTOR = -EL(1,NQ)*H
          MW = ML + MU + 1
          DO 260 J = 1,N
            DO 260 I = MAX(ML+1, MW+1-J), MIN(MW+N-J, MW+ML)
 260          DFDY(I,J) = FACTOR*DFDY(I,J)
        ELSE IF (MITER .EQ. 5) THEN
          BR = UROUND**(.875D0)
          BL = UROUND**(.75D0)
          BP = UROUND**(-.15D0)
          FACMIN = UROUND**(.78D0)
          MW = ML + MU + 1
          J2 = MIN(MW, N)
          DO 340 J = 1,J2
            DO 290 K = J,N,MW
              YS = MAX(ABS(YWT(K)), ABS(Y(K)))
 280          DY = FAC(K)*YS
              IF (DY .EQ. 0.D0) THEN
                IF (FAC(K) .LT. FACMAX) THEN
                  FAC(K) = MIN(100.D0*FAC(K), FACMAX)
                  GO TO 280
                ELSE
                  DY = YS
                END IF
              END IF
              IF (NQ .EQ. 1) THEN
                DY = SIGN(DY, SAVE2(K))
              ELSE
                DY = SIGN(DY, YH(K,3))
              END IF
              DY = (Y(K) + DY) - Y(K)
              DFDY(MW,K) = Y(K)
 290          Y(K) = Y(K) + DY
            CALL F (N, T, Y, SAVE1)
            IF (N .EQ. 0) THEN
              JSTATE = 6
              RETURN
            END IF
            DO 330 K = J,N,MW
              Y(K) = DFDY(MW,K)
              YS = MAX(ABS(YWT(K)), ABS(Y(K)))
              DY = FAC(K)*YS
              IF (DY .EQ. 0.D0) DY = YS
              IF (NQ .EQ. 1) THEN
                DY = SIGN(DY, SAVE2(K))
              ELSE
                DY = SIGN(DY, YH(K,3))
              END IF
              DY = (Y(K) + DY) - Y(K)
              FACTOR = -EL(1,NQ)*H/DY
              DO 300 I = MAX(ML+1, MW+1-K), MIN(MW+N-K, MW+ML)
 300            DFDY(I,K) = FACTOR*(SAVE1(I+K-MW) - SAVE2(I+K-MW))
C                                                                 Step 1
              IMAX = MAX(1, K - MU)
              DIFF = ABS(SAVE2(IMAX) - SAVE1(IMAX))
              DO 310 I = MAX(1, K - MU)+1, MIN(K + ML, N)
                IF (ABS(SAVE2(I) - SAVE1(I)) .GT. DIFF) THEN
                  IMAX = I
                  DIFF = ABS(SAVE2(I) - SAVE1(I))
                END IF
 310            CONTINUE
C                                                                 Step 2
              IF (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX))) .GT.0.D0) THEN
                SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
C                                                                 Step 3
                IF (DIFF .GT. BU*SCALE) THEN
                  FAC(J) = MAX(FACMIN, FAC(J)*.5D0)
                ELSE IF (BR*SCALE .LE.DIFF .AND. DIFF .LE.BL*SCALE) THEN
                  FAC(J) = MIN(FAC(J)*2.D0, FACMAX)
C                                                                 Step 4
                ELSE IF (DIFF .LT. BR*SCALE) THEN
                  FAC(K) = MIN(BP*FAC(K), FACMAX)
                END IF
              END IF
 330          CONTINUE
 340        CONTINUE
          NFE = NFE + J2
        END IF
        IF (ISWFLG .EQ. 3) THEN
          DFDYMX = 0.D0
          DO 345 J = 1,N
            DO 345 I = MAX(ML+1, MW+1-J), MIN(MW+N-J, MW+ML)
 345          DFDYMX = MAX(DFDYMX, ABS(DFDY(I,J)))
          BND = 0.D0
          IF (DFDYMX .NE. 0.D0) THEN
            DO 350 J = 1,N
              DO 350 I = MAX(ML+1, MW+1-J), MIN(MW+N-J, MW+ML)
 350            BND = BND + (DFDY(I,J)/DFDYMX)**2
            BND = DFDYMX*SQRT(BND)/(-EL(1,NQ)*H)
          END IF
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 360 J = 1,N
 360        DFDY(MW,J) = DFDY(MW,J) + 1.D0
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 380 J = 1,N
            DO 380 I = MAX(ML+1, MW+1-J), MIN(MW+N-J, MW+ML)
 380          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 400 J = 1,NDE
 400        DFDY(MW,J) =  DFDY(MW,J) + A(J,1)
        ELSE IF (IMPL .EQ. 3) THEN
          CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 390 J = 1,NDE
            DO 390 I = MAX(ML+1, MW+1-J), MIN(MW+NDE-J, MW+ML)
 390          DFDY(I,J) = DFDY(I,J) + A(I,J)
        END IF
        CALL DGBFA (DFDY, MATDIM, N, ML, MU, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 1
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
        IF (IFLAG .EQ. -1) THEN
          IER = .TRUE.
          RETURN
        END IF
        IF (N .EQ. 0) THEN
          JSTATE = 10
          RETURN
        END IF
      END IF
      RETURN
      END
