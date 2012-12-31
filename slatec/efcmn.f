*DECK EFCMN
      SUBROUTINE EFCMN (NDATA, XDATA, YDATA, SDDATA, NORD, NBKPT,
     +   BKPTIN, MDEIN, MDEOUT, COEFF, BF, XTEMP, PTEMP, BKPT, G, MDG,
     +   W, MDW, LW)
C***BEGIN PROLOGUE  EFCMN
C***SUBSIDIARY
C***PURPOSE  Subsidiary to EFC
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (EFCMN-S, DEFCMN-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C     This is a companion subprogram to EFC( ).
C     This subprogram does weighted least squares fitting of data by
C     B-spline curves.
C     The documentation for EFC( ) has complete usage instructions.
C
C***SEE ALSO  EFC
C***ROUTINES CALLED  BNDACC, BNDSOL, BSPLVN, SCOPY, SSCAL, SSORT, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and extensively revised (WRB & RWC)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C***END PROLOGUE  EFCMN
      INTEGER LW, MDEIN, MDEOUT, MDG, MDW, NBKPT, NDATA, NORD
      REAL             BF(NORD,*), BKPT(*), BKPTIN(*), COEFF(*),
     *   G(MDG,*), PTEMP(*), SDDATA(*), W(MDW,*), XDATA(*), XTEMP(*),
     *   YDATA(*)
C
      EXTERNAL BNDACC, BNDSOL, BSPLVN, SCOPY, SSCAL, SSORT, XERMSG
C
      REAL             DUMMY, RNORM, XMAX, XMIN, XVAL
      INTEGER I, IDATA, ILEFT, INTSEQ, IP, IR, IROW, L, MT, N, NB,
     *   NORDM1, NORDP1, NP1
      CHARACTER*8 XERN1, XERN2
C
C***FIRST EXECUTABLE STATEMENT  EFCMN
C
C     Initialize variables and analyze input.
C
      N = NBKPT - NORD
      NP1 = N + 1
C
C     Initially set all output coefficients to zero.
C
      CALL SCOPY (N, 0.E0, 0, COEFF, 1)
      MDEOUT = -1
      IF (NORD.LT.1 .OR. NORD.GT.20) THEN
         CALL XERMSG ('SLATEC', 'EFCMN',
     +      'IN EFC, THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.',
     +      3, 1)
         RETURN
      ENDIF
C
      IF (NBKPT.LT.2*NORD) THEN
         CALL XERMSG ('SLATEC', 'EFCMN',
     +      'IN EFC, THE NUMBER OF KNOTS MUST BE AT LEAST TWICE ' //
     +      'THE B-SPLINE ORDER.', 4, 1)
         RETURN
      ENDIF
C
      IF (NDATA.LT.0) THEN
         CALL XERMSG ('SLATEC', 'EFCMN',
     +      'IN EFC, THE NUMBER OF DATA POINTS MUST BE NONNEGATIVE.',
     +      5, 1)
         RETURN
      ENDIF
C
      NB = (NBKPT-NORD+3)*(NORD+1) + (NBKPT+1)*(NORD+1) +
     +     2*MAX(NBKPT,NDATA) + NBKPT + NORD**2
      IF (LW .LT. NB) THEN
         WRITE (XERN1, '(I8)') NB
         WRITE (XERN2, '(I8)') LW
         CALL XERMSG ('SLATEC', 'EFCMN',
     *      'IN EFC, INSUFFICIENT STORAGE FOR W(*).  CHECK FORMULA ' //
     *      'THAT READS LW.GE. ... .  NEED = ' // XERN1 //
     *      ' GIVEN = ' // XERN2, 6, 1)
         MDEOUT = -1
         RETURN
      ENDIF
C
      IF (MDEIN.NE.1 .AND. MDEIN.NE.2) THEN
         CALL XERMSG ('SLATEC', 'EFCMN',
     +      'IN EFC, INPUT VALUE OF MDEIN MUST BE 1-2.', 7, 1)
         RETURN
      ENDIF
C
C     Sort the breakpoints.
C
      CALL SCOPY (NBKPT, BKPTIN, 1, BKPT, 1)
      CALL SSORT (BKPT, DUMMY, NBKPT, 1)
C
C     Save interval containing knots.
C
      XMIN = BKPT(NORD)
      XMAX = BKPT(NP1)
      NORDM1 = NORD - 1
      NORDP1 = NORD + 1
C
C     Process least squares equations.
C
C     Sort data and an array of pointers.
C
      CALL SCOPY (NDATA, XDATA, 1, XTEMP, 1)
      DO 100 I = 1,NDATA
         PTEMP(I) = I
  100 CONTINUE
C
      IF (NDATA.GT.0) THEN
         CALL SSORT (XTEMP, PTEMP, NDATA, 2)
         XMIN = MIN(XMIN,XTEMP(1))
         XMAX = MAX(XMAX,XTEMP(NDATA))
      ENDIF
C
C     Fix breakpoint array if needed. This should only involve very
C     minor differences with the input array of breakpoints.
C
      DO 110 I = 1,NORD
         BKPT(I) = MIN(BKPT(I),XMIN)
  110 CONTINUE
C
      DO 120 I = NP1,NBKPT
         BKPT(I) = MAX(BKPT(I),XMAX)
  120 CONTINUE
C
C     Initialize parameters of banded matrix processor, BNDACC( ).
C
      MT = 0
      IP = 1
      IR = 1
      ILEFT = NORD
      INTSEQ = 1
      DO 150 IDATA = 1,NDATA
C
C        Sorted indices are in PTEMP(*).
C
         L = PTEMP(IDATA)
         XVAL = XDATA(L)
C
C        When interval changes, process equations in the last block.
C
         IF (XVAL.GE.BKPT(ILEFT+1)) THEN
            CALL BNDACC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
            MT = 0
C
C           Move pointer up to have BKPT(ILEFT).LE.XVAL, ILEFT.LE.N.
C
            DO 130 ILEFT = ILEFT,N
               IF (XVAL.LT.BKPT(ILEFT+1)) GO TO 140
               IF (MDEIN.EQ.2) THEN
C
C                 Data is being sequentially accumulated.
C                 Transfer previously accumulated rows from W(*,*) to
C                 G(*,*) and process them.
C
                  CALL SCOPY (NORDP1, W(INTSEQ,1), MDW, G(IR,1), MDG)
                  CALL BNDACC (G, MDG, NORD, IP, IR, 1, INTSEQ)
                  INTSEQ = INTSEQ + 1
               ENDIF
  130       CONTINUE
         ENDIF
C
C        Obtain B-spline function value.
C
  140    CALL BSPLVN (BKPT, NORD, 1, XVAL, ILEFT, BF)
C
C        Move row into place.
C
         IROW = IR + MT
         MT = MT + 1
         CALL SCOPY (NORD, BF, 1, G(IROW,1), MDG)
         G(IROW,NORDP1) = YDATA(L)
C
C        Scale data if uncertainty is nonzero.
C
         IF (SDDATA(L).NE.0.E0) CALL SSCAL (NORDP1, 1.E0/SDDATA(L),
     +                               G(IROW,1), MDG)
C
C        When staging work area is exhausted, process rows.
C
         IF (IROW.EQ.MDG-1) THEN
            CALL BNDACC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
            MT = 0
         ENDIF
  150 CONTINUE
C
C     Process last block of equations.
C
      CALL BNDACC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
C
C     Finish processing any previously accumulated rows from W(*,*)
C     to G(*,*).
C
      IF (MDEIN.EQ.2) THEN
         DO 160 I = INTSEQ,NP1
            CALL SCOPY (NORDP1, W(I,1), MDW, G(IR,1), MDG)
            CALL BNDACC (G, MDG, NORD, IP, IR, 1, MIN(N,I))
  160    CONTINUE
      ENDIF
C
C     Last call to adjust block positioning.
C
      CALL SCOPY (NORDP1, 0.E0, 0, G(IR,1), MDG)
      CALL BNDACC (G, MDG, NORD, IP, IR, 1, NP1)
C
C     Transfer accumulated rows from G(*,*) to W(*,*) for
C     possible later sequential accumulation.
C
      DO 170 I = 1,NP1
         CALL SCOPY (NORDP1, G(I,1), MDG, W(I,1), MDW)
  170 CONTINUE
C
C     Solve for coefficients when possible.
C
      DO 180 I = 1,N
         IF (G(I,1).EQ.0.E0) THEN
            MDEOUT = 2
            RETURN
         ENDIF
  180 CONTINUE
C
C     All the diagonal terms in the accumulated triangular
C     matrix are nonzero.  The solution can be computed but
C     it may be unsuitable for further use due to poor
C     conditioning or the lack of constraints.  No checking
C     for either of these is done here.
C
      CALL BNDSOL (1, G, MDG, NORD, IP, IR, COEFF, N, RNORM)
      MDEOUT = 1
      RETURN
      END
