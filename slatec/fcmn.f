*DECK FCMN
      SUBROUTINE FCMN (NDATA, XDATA, YDATA, SDDATA, NORD, NBKPT, BKPTIN,
     +   NCONST, XCONST, YCONST, NDERIV, MODE, COEFF, BF, XTEMP, PTEMP,
     +   BKPT, G, MDG, W, MDW, WORK, IWORK)
C***BEGIN PROLOGUE  FCMN
C***SUBSIDIARY
C***PURPOSE  Subsidiary to FC
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (FCMN-S, DFCMN-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This is a companion subprogram to FC( ).
C     The documentation for FC( ) has complete usage instructions.
C
C***SEE ALSO  FC
C***ROUTINES CALLED  BNDACC, BNDSOL, BSPLVD, BSPLVN, LSEI, SAXPY, SCOPY,
C                    SSCAL, SSORT, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   780801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and extensively revised (WRB & RWC)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C***END PROLOGUE  FCMN
      INTEGER IWORK(*), MDG, MDW, MODE, NBKPT, NCONST, NDATA, NDERIV(*),
     *   NORD
      REAL             BF(NORD,*), BKPT(*), BKPTIN(*), COEFF(*),
     *   G(MDG,*), PTEMP(*), SDDATA(*), W(MDW,*), WORK(*),
     *   XCONST(*), XDATA(*), XTEMP(*), YCONST(*), YDATA(*)
C
      EXTERNAL BNDACC, BNDSOL, BSPLVD, BSPLVN, LSEI, SAXPY, SCOPY,
     *    SSCAL, SSORT, XERMSG
C
      REAL             DUMMY, PRGOPT(10), RNORM, RNORME, RNORML, XMAX,
     *   XMIN, XVAL, YVAL
      INTEGER I, IDATA, IDERIV, ILEFT, INTRVL, INTW1, IP, IR, IROW,
     *   ITYPE, IW1, IW2, L, LW, MT, N, NB, NEQCON, NINCON, NORDM1,
     *   NORDP1, NP1
      LOGICAL BAND, NEW, VAR
      CHARACTER*8 XERN1
C
C***FIRST EXECUTABLE STATEMENT  FCMN
C
C     Analyze input.
C
      IF (NORD.LT.1 .OR. NORD.GT.20) THEN
         CALL XERMSG ('SLATEC', 'FCMN',
     +      'IN FC, THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.',
     +      2, 1)
         MODE = -1
         RETURN
C
      ELSEIF (NBKPT.LT.2*NORD) THEN
         CALL XERMSG ('SLATEC', 'FCMN',
     +      'IN FC, THE NUMBER OF KNOTS MUST BE AT LEAST TWICE ' //
     +      'THE B-SPLINE ORDER.', 2, 1)
         MODE = -1
         RETURN
      ENDIF
C
      IF (NDATA.LT.0) THEN
         CALL XERMSG ('SLATEC', 'FCMN',
     +      'IN FC, THE NUMBER OF DATA POINTS MUST BE NONNEGATIVE.',
     +      2, 1)
         MODE = -1
         RETURN
      ENDIF
C
C     Amount of storage allocated for W(*), IW(*).
C
      IW1 = IWORK(1)
      IW2 = IWORK(2)
      NB = (NBKPT-NORD+3)*(NORD+1) + 2*MAX(NDATA,NBKPT) + NBKPT +
     +     NORD**2
C
C     See if sufficient storage has been allocated.
C
      IF (IW1.LT.NB) THEN
         WRITE (XERN1, '(I8)') NB
         CALL XERMSG ('SLATEC', 'FCMN',
     *      'IN FC, INSUFFICIENT STORAGE FOR W(*).  CHECK NB = ' //
     *      XERN1, 2, 1)
         MODE = -1
         RETURN
      ENDIF
C
      IF (MODE.EQ.1) THEN
         BAND = .TRUE.
         VAR = .FALSE.
         NEW = .TRUE.
      ELSEIF (MODE.EQ.2) THEN
         BAND = .FALSE.
         VAR = .TRUE.
         NEW = .TRUE.
      ELSEIF (MODE.EQ.3) THEN
         BAND = .TRUE.
         VAR = .FALSE.
         NEW = .FALSE.
      ELSEIF (MODE.EQ.4) THEN
         BAND = .FALSE.
         VAR = .TRUE.
         NEW = .FALSE.
      ELSE
         CALL XERMSG ('SLATEC', 'FCMN',
     +      'IN FC, INPUT VALUE OF MODE MUST BE 1-4.', 2, 1)
         MODE = -1
         RETURN
      ENDIF
      MODE = 0
C
C     Sort the breakpoints.
C
      CALL SCOPY (NBKPT, BKPTIN, 1, BKPT, 1)
      CALL SSORT (BKPT, DUMMY, NBKPT, 1)
C
C     Initialize variables.
C
      NEQCON = 0
      NINCON = 0
      DO 100 I = 1,NCONST
         L = NDERIV(I)
         ITYPE = MOD(L,4)
         IF (ITYPE.LT.2) THEN
            NINCON = NINCON + 1
         ELSE
            NEQCON = NEQCON + 1
         ENDIF
  100 CONTINUE
C
C     Compute the number of variables.
C
      N = NBKPT - NORD
      NP1 = N + 1
      LW = NB + (NP1+NCONST)*NP1 + 2*(NEQCON+NP1) + (NINCON+NP1) +
     +     (NINCON+2)*(NP1+6)
      INTW1 = NINCON + 2*NP1
C
C     Save interval containing knots.
C
      XMIN = BKPT(NORD)
      XMAX = BKPT(NP1)
C
C     Find the smallest referenced independent variable value in any
C     constraint.
C
      DO 110 I = 1,NCONST
         XMIN = MIN(XMIN,XCONST(I))
         XMAX = MAX(XMAX,XCONST(I))
  110 CONTINUE
      NORDM1 = NORD - 1
      NORDP1 = NORD + 1
C
C     Define the option vector PRGOPT(1-10) for use in LSEI( ).
C
      PRGOPT(1) = 4
C
C     Set the covariance matrix computation flag.
C
      PRGOPT(2) = 1
      IF (VAR) THEN
         PRGOPT(3) = 1
      ELSE
         PRGOPT(3) = 0
      ENDIF
C
C     Increase the rank determination tolerances for both equality
C     constraint equations and least squares equations.
C
      PRGOPT(4) = 7
      PRGOPT(5) = 4
      PRGOPT(6) = 1.E-4
C
      PRGOPT(7) = 10
      PRGOPT(8) = 5
      PRGOPT(9) = 1.E-4
C
      PRGOPT(10) = 1
C
C     Turn off work array length checking in LSEI( ).
C
      IWORK(1) = 0
      IWORK(2) = 0
C
C     Initialize variables and analyze input.
C
      IF (NEW) THEN
C
C        To process least squares equations sort data and an array of
C        pointers.
C
         CALL SCOPY (NDATA, XDATA, 1, XTEMP, 1)
         DO 120 I = 1,NDATA
            PTEMP(I) = I
  120    CONTINUE
C
         IF (NDATA.GT.0) THEN
            CALL SSORT (XTEMP, PTEMP, NDATA, 2)
            XMIN = MIN(XMIN,XTEMP(1))
            XMAX = MAX(XMAX,XTEMP(NDATA))
         ENDIF
C
C        Fix breakpoint array if needed.
C
         DO 130 I = 1,NORD
            BKPT(I) = MIN(BKPT(I),XMIN)
  130    CONTINUE
C
         DO 140 I = NP1,NBKPT
            BKPT(I) = MAX(BKPT(I),XMAX)
  140    CONTINUE
C
C        Initialize parameters of banded matrix processor, BNDACC( ).
C
         MT = 0
         IP = 1
         IR = 1
         ILEFT = NORD
         DO 160 IDATA = 1,NDATA
C
C           Sorted indices are in PTEMP(*).
C
            L = PTEMP(IDATA)
            XVAL = XDATA(L)
C
C           When interval changes, process equations in the last block.
C
            IF (XVAL.GE.BKPT(ILEFT+1)) THEN
               CALL BNDACC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
               MT = 0
C
C              Move pointer up to have BKPT(ILEFT).LE.XVAL,
C                 ILEFT.LT.NP1.
C
  150          IF (XVAL.GE.BKPT(ILEFT+1) .AND. ILEFT.LT.N) THEN
                  ILEFT = ILEFT + 1
                  GO TO 150
               ENDIF
            ENDIF
C
C           Obtain B-spline function value.
C
            CALL BSPLVN (BKPT, NORD, 1, XVAL, ILEFT, BF)
C
C           Move row into place.
C
            IROW = IR + MT
            MT = MT + 1
            CALL SCOPY (NORD, BF, 1, G(IROW,1), MDG)
            G(IROW,NORDP1) = YDATA(L)
C
C           Scale data if uncertainty is nonzero.
C
            IF (SDDATA(L).NE.0.E0) CALL SSCAL (NORDP1, 1.E0/SDDATA(L),
     +                                  G(IROW,1), MDG)
C
C           When staging work area is exhausted, process rows.
C
            IF (IROW.EQ.MDG-1) THEN
               CALL BNDACC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
               MT = 0
            ENDIF
  160    CONTINUE
C
C        Process last block of equations.
C
         CALL BNDACC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
C
C        Last call to adjust block positioning.
C
         CALL SCOPY (NORDP1, 0.E0, 0, G(IR,1), MDG)
         CALL BNDACC (G, MDG, NORD, IP, IR, 1, NP1)
      ENDIF
C
      BAND = BAND .AND. NCONST.EQ.0
      DO 170 I = 1,N
         BAND = BAND .AND. G(I,1).NE.0.E0
  170 CONTINUE
C
C     Process banded least squares equations.
C
      IF (BAND) THEN
         CALL BNDSOL (1, G, MDG, NORD, IP, IR, COEFF, N, RNORM)
         RETURN
      ENDIF
C
C     Check further for sufficient storage in working arrays.
C
      IF (IW1.LT.LW) THEN
         WRITE (XERN1, '(I8)') LW
         CALL XERMSG ('SLATEC', 'FCMN',
     *      'IN FC, INSUFFICIENT STORAGE FOR W(*).  CHECK LW = ' //
     *      XERN1, 2, 1)
         MODE = -1
         RETURN
      ENDIF
C
      IF (IW2.LT.INTW1) THEN
         WRITE (XERN1, '(I8)') INTW1
         CALL XERMSG ('SLATEC', 'FCMN',
     *      'IN FC, INSUFFICIENT STORAGE FOR IW(*).  CHECK IW1 = ' //
     *      XERN1, 2, 1)
         MODE = -1
         RETURN
      ENDIF
C
C     Write equality constraints.
C     Analyze constraint indicators for an equality constraint.
C
      NEQCON = 0
      DO 220 IDATA = 1,NCONST
         L = NDERIV(IDATA)
         ITYPE = MOD(L,4)
         IF (ITYPE.GT.1) THEN
            IDERIV = L/4
            NEQCON = NEQCON + 1
            ILEFT = NORD
            XVAL = XCONST(IDATA)
C
  180       IF (XVAL.LT.BKPT(ILEFT+1) .OR. ILEFT.GE.N) GO TO 190
            ILEFT = ILEFT + 1
            GO TO 180
C
  190       CALL BSPLVD (BKPT, NORD, XVAL, ILEFT, BF, IDERIV+1)
            CALL SCOPY (NP1, 0.E0, 0, W(NEQCON,1), MDW)
            CALL SCOPY (NORD, BF(1,IDERIV+1), 1, W(NEQCON,ILEFT-NORDM1),
     +                  MDW)
C
            IF (ITYPE.EQ.2) THEN
               W(NEQCON,NP1) = YCONST(IDATA)
            ELSE
               ILEFT = NORD
               YVAL = YCONST(IDATA)
C
  200          IF (YVAL.LT.BKPT(ILEFT+1) .OR. ILEFT.GE.N) GO TO 210
               ILEFT = ILEFT + 1
               GO TO 200
C
  210          CALL BSPLVD (BKPT, NORD, YVAL, ILEFT, BF, IDERIV+1)
               CALL SAXPY (NORD, -1.E0, BF(1, IDERIV+1), 1,
     +                     W(NEQCON, ILEFT-NORDM1), MDW)
            ENDIF
         ENDIF
  220 CONTINUE
C
C     Transfer least squares data.
C
      DO 230 I = 1,NP1
         IROW = I + NEQCON
         CALL SCOPY (N, 0.E0, 0, W(IROW,1), MDW)
         CALL SCOPY (MIN(NP1-I, NORD), G(I,1), MDG, W(IROW,I), MDW)
         W(IROW,NP1) = G(I,NORDP1)
  230 CONTINUE
C
C     Write inequality constraints.
C     Analyze constraint indicators for inequality constraints.
C
      NINCON = 0
      DO 260 IDATA = 1,NCONST
         L = NDERIV(IDATA)
         ITYPE = MOD(L,4)
         IF (ITYPE.LT.2) THEN
            IDERIV = L/4
            NINCON = NINCON + 1
            ILEFT = NORD
            XVAL = XCONST(IDATA)
C
  240       IF (XVAL.LT.BKPT(ILEFT+1) .OR. ILEFT.GE.N) GO TO 250
            ILEFT = ILEFT + 1
            GO TO 240
C
  250       CALL BSPLVD (BKPT, NORD, XVAL, ILEFT, BF, IDERIV+1)
            IROW = NEQCON + NP1 + NINCON
            CALL SCOPY (N, 0.E0, 0, W(IROW,1), MDW)
            INTRVL = ILEFT - NORDM1
            CALL SCOPY (NORD, BF(1, IDERIV+1), 1, W(IROW, INTRVL), MDW)
C
            IF (ITYPE.EQ.1) THEN
               W(IROW,NP1) = YCONST(IDATA)
            ELSE
               W(IROW,NP1) = -YCONST(IDATA)
               CALL SSCAL (NORD, -1.E0, W(IROW, INTRVL), MDW)
            ENDIF
         ENDIF
  260 CONTINUE
C
C     Solve constrained least squares equations.
C
      CALL LSEI(W, MDW, NEQCON, NP1, NINCON, N, PRGOPT, COEFF, RNORME,
     +          RNORML, MODE, WORK, IWORK)
      RETURN
      END
