*DECK CDPST
      SUBROUTINE CDPST (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML,
     8   MU, N, NDE, NQ, SAVE2, T, USERS, Y, YH, YWT, UROUND, NFE, NJE,
     8   A, DFDY, FAC, IER, IPVT, SAVE1, ISWFLG, BND, JSTATE)
C***BEGIN PROLOGUE  CDPST
C***SUBSIDIARY
C***PURPOSE  Subroutine CDPST evaluates the Jacobian matrix of the right
C            hand side of the differential equations.
C***LIBRARY   SLATEC (SDRIVE)
C***TYPE      COMPLEX (SDPST-S, DDPST-D, CDPST-C)
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
C***ROUTINES CALLED  CGBFA, CGEFA, SCNRM2
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  CDPST
      INTEGER I, IFLAG, IMAX, IMPL, INFO, ISWFLG, J, J2, JSTATE, K,
     8        MATDIM, MITER, ML, MU, MW, N, NDE, NFE, NJE, NQ
      COMPLEX A(MATDIM,*), CFCTR, DFDY(MATDIM,*), DY, FAC(*), SAVE1(*),
     8        SAVE2(*), Y(*), YH(N,*), YJ, YS, YWT(*)
      REAL BL, BND, BP, BR, BU, DFDYMX, DIFF, EL(13,12), FACMAX, FACMIN,
     8     FACTOR, H, SCALE, SCNRM2, T, UROUND, ZMAX, ZMIN
      INTEGER IPVT(*)
      LOGICAL IER
      PARAMETER(FACMAX = .5E0, BU = 0.5E0)
C***FIRST EXECUTABLE STATEMENT  CDPST
      NJE = NJE + 1
      IER = .FALSE.
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (MITER .EQ. 1) THEN
          CALL JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
          IF (N .EQ. 0) THEN
            JSTATE = 8
            RETURN
          END IF
          IF (ISWFLG .EQ. 3) BND = SCNRM2(N*N, DFDY, 1)
          FACTOR = -EL(1,NQ)*H
          DO 110 J = 1,N
            DO 110 I = 1,N
 110          DFDY(I,J) = FACTOR*DFDY(I,J)
        ELSE IF (MITER .EQ. 2) THEN
          BR = UROUND**(.875E0)
          BL = UROUND**(.75E0)
          BP = UROUND**(-.15E0)
          FACMIN = UROUND**(.78E0)
          DO 170 J = 1,N
            IF (ABS(Y(J)) .GT. ABS(YWT(J))) THEN
              YS = Y(J)
            ELSE
              YS = YWT(J)
            END IF
 120        DY = FAC(J)*YS
            IF (DY .EQ. 0.E0) THEN
              IF (REAL(FAC(J)) .LT. FACMAX) THEN
                FAC(J) = MIN(100.E0*REAL(FAC(J)), FACMAX)
                GO TO 120
              ELSE
                DY = YS
              END IF
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
            CFCTR = -EL(1,NQ)*H/DY
            DO 140 I = 1,N
 140          DFDY(I,J) = (SAVE1(I) - SAVE2(I))*CFCTR
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
            IF (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX))) .GT. 0.E0) THEN
              SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
C                                                                 Step 3
              IF (DIFF .GT. BU*SCALE) THEN
                FAC(J) = MAX(FACMIN, REAL(FAC(J))*.5E0)
              ELSE IF (BR*SCALE .LE. DIFF .AND. DIFF .LE. BL*SCALE) THEN
                FAC(J) = MIN(REAL(FAC(J))*2.E0, FACMAX)
C                                                                 Step 4
              ELSE IF (DIFF .LT. BR*SCALE) THEN
                FAC(J) = MIN(BP*REAL(FAC(J)), FACMAX)
              END IF
            END IF
 170        CONTINUE
          IF (ISWFLG .EQ. 3) BND = SCNRM2(N*N, DFDY, 1)/(-EL(1,NQ)*H)
          NFE = NFE + N
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 190 I = 1,N
 190        DFDY(I,I) = DFDY(I,I) + 1.E0
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
        CALL CGEFA (DFDY, MATDIM, N, IPVT, INFO)
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
          BR = UROUND**(.875E0)
          BL = UROUND**(.75E0)
          BP = UROUND**(-.15E0)
          FACMIN = UROUND**(.78E0)
          MW = ML + MU + 1
          J2 = MIN(MW, N)
          DO 340 J = 1,J2
            DO 290 K = J,N,MW
              IF (ABS(Y(K)) .GT. ABS(YWT(K))) THEN
                YS = Y(K)
              ELSE
                YS = YWT(K)
              END IF
 280          DY = FAC(K)*YS
              IF (DY .EQ. 0.E0) THEN
                IF (REAL(FAC(K)) .LT. FACMAX) THEN
                  FAC(K) = MIN(100.E0*REAL(FAC(K)), FACMAX)
                  GO TO 280
                ELSE
                  DY = YS
                END IF
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
              DY = Y(K) - DFDY(MW,K)
              Y(K) = DFDY(MW,K)
              CFCTR = -EL(1,NQ)*H/DY
              DO 300 I = MAX(ML+1, MW+1-K), MIN(MW+N-K, MW+ML)
 300            DFDY(I,K) = CFCTR*(SAVE1(I+K-MW) - SAVE2(I+K-MW))
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
              IF (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX))) .GT.0.E0) THEN
                SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
C                                                                 Step 3
                IF (DIFF .GT. BU*SCALE) THEN
                  FAC(J) = MAX(FACMIN, REAL(FAC(J))*.5E0)
                ELSE IF (BR*SCALE .LE.DIFF .AND. DIFF .LE.BL*SCALE) THEN
                  FAC(J) = MIN(REAL(FAC(J))*2.E0, FACMAX)
C                                                                 Step 4
                ELSE IF (DIFF .LT. BR*SCALE) THEN
                  FAC(K) = MIN(BP*REAL(FAC(K)), FACMAX)
                END IF
              END IF
 330          CONTINUE
 340        CONTINUE
          NFE = NFE + J2
        END IF
        IF (ISWFLG .EQ. 3) THEN
          DFDYMX = 0.E0
          DO 345 J = 1,N
            DO 345 I = MAX(ML+1, MW+1-J), MIN(MW+N-J, MW+ML)
              ZMAX = MAX(ABS(REAL(DFDY(I,J))), ABS(AIMAG(DFDY(I,J))))
              ZMIN = MIN(ABS(REAL(DFDY(I,J))), ABS(AIMAG(DFDY(I,J))))
              IF (ZMAX .NE. 0.E0)
     8        DFDYMX = MAX(DFDYMX, ZMAX*SQRT(1.E0+ (ZMIN/ZMAX)**2))
 345          CONTINUE
          BND = 0.E0
          IF (DFDYMX .NE. 0.E0) THEN
            DO 350 J = 1,N
              DO 350 I = MAX(ML+1, MW+1-J), MIN(MW+N-J, MW+ML)
                BND = BND + (REAL(DFDY(I,J))/DFDYMX)**2 +
     8          (AIMAG(DFDY(I,J))/DFDYMX)**2
 350            CONTINUE
            BND = DFDYMX*SQRT(BND)/(-EL(1,NQ)*H)
          END IF
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 360 J = 1,N
 360        DFDY(MW,J) = DFDY(MW,J) + 1.E0
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
        CALL CGBFA (DFDY, MATDIM, N, ML, MU, IPVT, INFO)
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
