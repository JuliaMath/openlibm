*DECK DDNTL
      SUBROUTINE DDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8   MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T, UROUND, USERS,
     8   Y, YWT, H, MNTOLD, MTROLD, NFE, RC, YH, A, CONVRG, EL, FAC,
     8   IER, IPVT, NQ, NWAIT, RH, RMAX, SAVE2, TQ, TREND, ISWFLG,
     8   JSTATE)
C***BEGIN PROLOGUE  DDNTL
C***SUBSIDIARY
C***PURPOSE  Subroutine DDNTL is called to set parameters on the first
C            call to DDSTP, on an internal restart, or when the user has
C            altered MINT, MITER, and/or H.
C***LIBRARY   SLATEC (SDRIVE)
C***TYPE      DOUBLE PRECISION (SDNTL-S, DDNTL-D, CDNTL-C)
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             Mail Stop D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***DESCRIPTION
C
C  On the first call, the order is set to 1 and the initial derivatives
C  are calculated.  RMAX is the maximum ratio by which H can be
C  increased in one step.  It is initially RMINIT to compensate
C  for the small initial H, but then is normally equal to RMNORM.
C  If a failure occurs (in corrector convergence or error test), RMAX
C  is set at RMFAIL for the next increase.
C  If the caller has changed MINT, or if JTASK = 0, DDCST is called
C  to set the coefficients of the method.  If the caller has changed H,
C  YH must be rescaled.  If H or MINT has been changed, NWAIT is
C  reset to NQ + 2 to prevent further increases in H for that many
C  steps.  Also, RC is reset.  RC is the ratio of new to old values of
C  the coefficient L(0)*H.  If the caller has changed MITER, RC is
C  set to 0 to force the partials to be updated, if partials are used.
C
C***ROUTINES CALLED  DDCST, DDSCL, DGBFA, DGBSL, DGEFA, DGESL, DNRM2
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDNTL
      INTEGER I, IFLAG, IMPL, INFO, ISWFLG, JSTATE, JTASK, MATDIM,
     8        MAXORD, MINT, MITER, ML, MNTOLD, MTROLD, MU, N, NDE, NFE,
     8        NQ, NWAIT
      DOUBLE PRECISION A(MATDIM,*), EL(13,12), EPS, FAC(*), H, HMAX,
     8     HOLD, OLDL0, RC, RH, RMAX, RMINIT, SAVE1(*), SAVE2(*), DNRM2,
     8     SUM, T, TQ(3,12), TREND, UROUND, Y(*), YH(N,*), YWT(*)
      INTEGER IPVT(*)
      LOGICAL CONVRG, IER
      PARAMETER(RMINIT = 10000.D0)
C***FIRST EXECUTABLE STATEMENT  DDNTL
      IER = .FALSE.
      IF (JTASK .GE. 0) THEN
        IF (JTASK .EQ. 0) THEN
          CALL DDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
          RMAX = RMINIT
        END IF
        RC = 0.D0
        CONVRG = .FALSE.
        TREND = 1.D0
        NQ = 1
        NWAIT = 3
        CALL F (N, T, Y, SAVE2)
        IF (N .EQ. 0) THEN
          JSTATE = 6
          RETURN
        END IF
        NFE = NFE + 1
        IF (IMPL .NE. 0) THEN
          IF (MITER .EQ. 3) THEN
            IFLAG = 0
            CALL USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL, IMPL, N,
     8                  NDE, IFLAG)
            IF (IFLAG .EQ. -1) THEN
              IER = .TRUE.
              RETURN
            END IF
            IF (N .EQ. 0) THEN
              JSTATE = 10
              RETURN
            END IF
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL DGEFA (A, MATDIM, N, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGESL (A, MATDIM, N, IPVT, SAVE2, 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL DGBFA (A, MATDIM, N, ML, MU, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGBSL (A, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
            DO 150 I = 1,NDE
              IF (A(I,1) .EQ. 0.D0) THEN
                IER = .TRUE.
                RETURN
              ELSE
                SAVE2(I) = SAVE2(I)/A(I,1)
              END IF
 150          CONTINUE
            DO 155 I = NDE+1,N
 155          A(I,1) = 0.D0
          ELSE IF (IMPL .EQ. 3) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL DGEFA (A, MATDIM, NDE, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGESL (A, MATDIM, NDE, IPVT, SAVE2, 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL DGBFA (A, MATDIM, NDE, ML, MU, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGBSL (A, MATDIM, NDE, ML, MU, IPVT, SAVE2, 0)
            END IF
          END IF
        END IF
        DO 170 I = 1,NDE
 170      SAVE1(I) = SAVE2(I)/MAX(1.D0, YWT(I))
        SUM = DNRM2(NDE, SAVE1, 1)/SQRT(DBLE(NDE))
        IF (SUM .GT. EPS/ABS(H)) H = SIGN(EPS/SUM, H)
        DO 180 I = 1,N
 180      YH(I,2) = H*SAVE2(I)
        IF (MITER .EQ. 2 .OR. MITER .EQ. 5 .OR. ISWFLG .EQ. 3) THEN
          DO 20 I = 1,N
 20         FAC(I) = SQRT(UROUND)
        END IF
      ELSE
        IF (MITER .NE. MTROLD) THEN
          MTROLD = MITER
          RC = 0.D0
          CONVRG = .FALSE.
        END IF
        IF (MINT .NE. MNTOLD) THEN
          MNTOLD = MINT
          OLDL0 = EL(1,NQ)
          CALL DDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
          RC = RC*EL(1,NQ)/OLDL0
          NWAIT = NQ + 2
        END IF
        IF (H .NE. HOLD) THEN
          NWAIT = NQ + 2
          RH = H/HOLD
          CALL DDSCL (HMAX, N, NQ, RMAX,  HOLD, RC, RH, YH)
        END IF
      END IF
      RETURN
      END
