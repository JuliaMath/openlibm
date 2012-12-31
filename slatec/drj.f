*DECK DRJ
      DOUBLE PRECISION FUNCTION DRJ (X, Y, Z, P, IER)
C***BEGIN PROLOGUE  DRJ
C***PURPOSE  Compute the incomplete or complete (X or Y or Z is zero)
C            elliptic integral of the 3rd kind.  For X, Y, and Z non-
C            negative, at most one of them zero, and P positive,
C             RJ(X,Y,Z,P) = Integral from zero to infinity of
C                              -1/2     -1/2     -1/2     -1
C                    (3/2)(t+X)    (t+Y)    (t+Z)    (t+P)  dt.
C***LIBRARY   SLATEC
C***CATEGORY  C14
C***TYPE      DOUBLE PRECISION (RJ-S, DRJ-D)
C***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
C             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE THIRD KIND,
C             TAYLOR SERIES
C***AUTHOR  Carlson, B. C.
C             Ames Laboratory-DOE
C             Iowa State University
C             Ames, IA  50011
C           Notis, E. M.
C             Ames Laboratory-DOE
C             Iowa State University
C             Ames, IA  50011
C           Pexton, R. L.
C             Lawrence Livermore National Laboratory
C             Livermore, CA  94550
C***DESCRIPTION
C
C   1.     DRJ
C          Standard FORTRAN function routine
C          Double precision version
C          The routine calculates an approximation result to
C          DRJ(X,Y,Z,P) = Integral from zero to infinity of
C
C                                -1/2     -1/2     -1/2     -1
C                      (3/2)(t+X)    (t+Y)    (t+Z)    (t+P)  dt,
C
C          where X, Y, and Z are nonnegative, at most one of them is
C          zero, and P is positive.  If X or Y or Z is zero, the
C          integral is COMPLETE.  The duplication theorem is iterated
C          until the variables are nearly equal, and the function is
C          then expanded in Taylor series to fifth order.
C
C
C   2.     Calling Sequence
C          DRJ( X, Y, Z, P, IER )
C
C          Parameters on Entry
C          Values assigned by the calling routine
C
C          X      - Double precision, nonnegative variable
C
C          Y      - Double precision, nonnegative variable
C
C          Z      - Double precision, nonnegative variable
C
C          P      - Double precision, positive variable
C
C
C          On  Return    (values assigned by the DRJ routine)
C
C          DRJ     - Double precision approximation to the integral
C
C          IER    - Integer
C
C                   IER = 0 Normal and reliable termination of the
C                           routine. It is assumed that the requested
C                           accuracy has been achieved.
C
C                   IER >  0 Abnormal termination of the routine
C
C
C          X, Y, Z, P are unaltered.
C
C
C   3.    Error Messages
C
C         Value of IER assigned by the DRJ routine
C
C              Value assigned         Error Message printed
C              IER = 1                MIN(X,Y,Z) .LT. 0.0D0
C                  = 2                MIN(X+Y,X+Z,Y+Z,P) .LT. LOLIM
C                  = 3                MAX(X,Y,Z,P) .GT. UPLIM
C
C
C
C   4.     Control Parameters
C
C                  Values of LOLIM, UPLIM, and ERRTOL are set by the
C                  routine.
C
C
C          LOLIM and UPLIM determine the valid range of X, Y, Z, and P
C
C          LOLIM is not less than the cube root of the value
C          of LOLIM used in the routine for DRC.
C
C          UPLIM is not greater than 0.3 times the cube root of
C          the value of UPLIM used in the routine for DRC.
C
C
C                     Acceptable values for:   LOLIM      UPLIM
C                     IBM 360/370 SERIES   :   2.0D-26     3.0D+24
C                     CDC 6000/7000 SERIES :   5.0D-98     3.0D+106
C                     UNIVAC 1100 SERIES   :   5.0D-103    6.0D+101
C                     CRAY                 :   1.32D-822   1.4D+821
C                     VAX 11 SERIES        :   2.5D-13     9.0D+11
C
C
C
C          ERRTOL determines the accuracy of the answer
C
C                 the value assigned by the routine will result
C                 in solution precision within 1-2 decimals of
C                 "machine precision".
C
C
C
C
C          Relative error due to truncation of the series for DRJ
C          is less than 3 * ERRTOL ** 6 / (1 - ERRTOL) ** 3/2.
C
C
C
C        The accuracy of the computed approximation to the integral
C        can be controlled by choosing the value of ERRTOL.
C        Truncation of a Taylor series after terms of fifth order
C        introduces an error less than the amount shown in the
C        second column of the following table for each value of
C        ERRTOL in the first column.  In addition to the truncation
C        error there will be round-off error, but in practice the
C        total error from both sources is usually less than the
C        amount given in the table.
C
C
C
C          Sample choices:  ERRTOL   Relative truncation
C                                    error less than
C                           1.0D-3    4.0D-18
C                           3.0D-3    3.0D-15
C                           1.0D-2    4.0D-12
C                           3.0D-2    3.0D-9
C                           1.0D-1    4.0D-6
C
C                    Decreasing ERRTOL by a factor of 10 yields six more
C                    decimal digits of accuracy at the expense of one or
C                    two more iterations of the duplication theorem.
C
C *Long Description:
C
C   DRJ Special Comments
C
C
C     Check by addition theorem: DRJ(X,X+Z,X+W,X+P)
C     + DRJ(Y,Y+Z,Y+W,Y+P) + (A-B) * DRJ(A,B,B,A) + 3.0D0 / SQRT(A)
C     = DRJ(0,Z,W,P), where X,Y,Z,W,P are positive and X * Y
C     = Z * W,  A = P * P * (X+Y+Z+W),  B = P * (P+X) * (P+Y),
C     and B - A = P * (P-Z) * (P-W).  The sum of the third and
C     fourth terms on the left side is 3.0D0 * DRC(A,B).
C
C
C          On Input:
C
C     X, Y, Z, and P are the variables in the integral DRJ(X,Y,Z,P).
C
C
C          On Output:
C
C
C          X, Y, Z, P are unaltered.
C
C          ********************************************************
C
C          WARNING: Changes in the program may improve speed at the
C                   expense of robustness.
C
C    -------------------------------------------------------------------
C
C
C   Special double precision functions via DRJ and DRF
C
C
C                  Legendre form of ELLIPTIC INTEGRAL of 3rd kind
C                  -----------------------------------------
C
C
C                          PHI         2         -1
C             P(PHI,K,N) = INT (1+N SIN (THETA) )   *
C                           0
C
C
C                                  2    2         -1/2
C                             *(1-K  SIN (THETA) )     D THETA
C
C
C                                           2          2   2
C                        = SIN (PHI) DRF(COS (PHI), 1-K SIN (PHI),1)
C
C                                   3             2         2   2
C                         -(N/3) SIN (PHI) DRJ(COS (PHI),1-K SIN (PHI),
C
C                                  2
C                         1,1+N SIN (PHI))
C
C
C
C                  Bulirsch form of ELLIPTIC INTEGRAL of 3rd kind
C                  -----------------------------------------
C
C
C                                            2 2    2
C                  EL3(X,KC,P) = X DRF(1,1+KC X ,1+X ) +
C
C                                            3           2 2    2     2
C                               +(1/3)(1-P) X  DRJ(1,1+KC X ,1+X ,1+PX )
C
C
C                                           2
C                  CEL(KC,P,A,B) = A RF(0,KC ,1) +
C
C
C                                                      2
C                                 +(1/3)(B-PA) DRJ(0,KC ,1,P)
C
C
C                  Heuman's LAMBDA function
C                  -----------------------------------------
C
C
C                                2                      2      2    1/2
C                  L(A,B,P) =(COS (A)SIN(B)COS(B)/(1-COS (A)SIN (B))   )
C
C                                            2         2       2
C                            *(SIN(P) DRF(COS (P),1-SIN (A) SIN (P),1)
C
C                                 2       3            2       2
C                            +(SIN (A) SIN (P)/(3(1-COS (A) SIN (B))))
C
C                                    2         2       2
C                            *DRJ(COS (P),1-SIN (A) SIN (P),1,1-
C
C                                2       2          2       2
C                            -SIN (A) SIN (P)/(1-COS (A) SIN (B))))
C
C
C
C                  (PI/2) LAMBDA0(A,B) =L(A,B,PI/2) =
C
C                   2                         2       2    -1/2
C              = COS (A)  SIN(B) COS(B) (1-COS (A) SIN (B))
C
C                           2                  2       2
C                 *DRF(0,COS (A),1) + (1/3) SIN (A) COS (A)
C
C                                      2       2    -3/2
C                 *SIN(B) COS(B) (1-COS (A) SIN (B))
C
C                           2         2       2          2       2
C                 *DRJ(0,COS (A),1,COS (A) COS (B)/(1-COS (A) SIN (B)))
C
C
C                  Jacobi ZETA function
C                  -----------------------------------------
C
C                        2                     2   2    1/2
C             Z(B,K) = (K/3) SIN(B) COS(B) (1-K SIN (B))
C
C
C                                  2      2   2                 2
C                        *DRJ(0,1-K ,1,1-K SIN (B)) / DRF (0,1-K ,1)
C
C
C  ---------------------------------------------------------------------
C
C***REFERENCES  B. C. Carlson and E. M. Notis, Algorithms for incomplete
C                 elliptic integrals, ACM Transactions on Mathematical
C                 Software 7, 3 (September 1981), pp. 398-403.
C               B. C. Carlson, Computing elliptic integrals by
C                 duplication, Numerische Mathematik 33, (1979),
C                 pp. 1-16.
C               B. C. Carlson, Elliptic integrals of the first kind,
C                 SIAM Journal of Mathematical Analysis 8, (1977),
C                 pp. 231-242.
C***ROUTINES CALLED  D1MACH, DRC, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891009  Removed unreferenced statement labels.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900510  Changed calls to XERMSG to standard form, and some
C           editorial changes.  (RWC)).
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DRJ
      INTEGER IER
      CHARACTER*16 XERN3, XERN4, XERN5, XERN6, XERN7
      DOUBLE PRECISION ALFA, BETA, C1, C2, C3, C4, EA, EB, EC, E2, E3
      DOUBLE PRECISION LOLIM, UPLIM, EPSLON, ERRTOL, D1MACH
      DOUBLE PRECISION LAMDA, MU, P, PN, PNDEV
      DOUBLE PRECISION POWER4, DRC, SIGMA, S1, S2, S3, X, XN, XNDEV
      DOUBLE PRECISION XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV,
     * ZNROOT
      LOGICAL FIRST
      SAVE ERRTOL,LOLIM,UPLIM,C1,C2,C3,C4,FIRST
      DATA FIRST /.TRUE./
C
C***FIRST EXECUTABLE STATEMENT  DRJ
      IF (FIRST) THEN
         ERRTOL = (D1MACH(3)/3.0D0)**(1.0D0/6.0D0)
         LOLIM  = (5.0D0 * D1MACH(1))**(1.0D0/3.0D0)
         UPLIM  = 0.30D0*( D1MACH(2) / 5.0D0)**(1.0D0/3.0D0)
C
         C1 = 3.0D0/14.0D0
         C2 = 1.0D0/3.0D0
         C3 = 3.0D0/22.0D0
         C4 = 3.0D0/26.0D0
      ENDIF
      FIRST = .FALSE.
C
C         CALL ERROR HANDLER IF NECESSARY.
C
      DRJ = 0.0D0
      IF (MIN(X,Y,Z).LT.0.0D0) THEN
         IER = 1
         WRITE (XERN3, '(1PE15.6)') X
         WRITE (XERN4, '(1PE15.6)') Y
         WRITE (XERN5, '(1PE15.6)') Z
         CALL XERMSG ('SLATEC', 'DRJ',
     *      'MIN(X,Y,Z).LT.0 WHERE X = ' // XERN3 // ' Y = ' // XERN4 //
     *      ' AND Z = ' // XERN5, 1, 1)
         RETURN
      ENDIF
C
      IF (MAX(X,Y,Z,P).GT.UPLIM) THEN
         IER = 3
         WRITE (XERN3, '(1PE15.6)') X
         WRITE (XERN4, '(1PE15.6)') Y
         WRITE (XERN5, '(1PE15.6)') Z
         WRITE (XERN6, '(1PE15.6)') P
         WRITE (XERN7, '(1PE15.6)') UPLIM
         CALL XERMSG ('SLATEC', 'DRJ',
     *      'MAX(X,Y,Z,P).GT.UPLIM WHERE X = ' // XERN3 // ' Y = ' //
     *      XERN4 // ' Z = ' // XERN5 // ' P = ' // XERN6 //
     *      ' AND UPLIM = ' // XERN7, 3, 1)
         RETURN
      ENDIF
C
      IF (MIN(X+Y,X+Z,Y+Z,P).LT.LOLIM) THEN
         IER = 2
         WRITE (XERN3, '(1PE15.6)') X
         WRITE (XERN4, '(1PE15.6)') Y
         WRITE (XERN5, '(1PE15.6)') Z
         WRITE (XERN6, '(1PE15.6)') P
         WRITE (XERN7, '(1PE15.6)') LOLIM
         CALL XERMSG ('SLATEC', 'RJ',
     *      'MIN(X+Y,X+Z,Y+Z,P).LT.LOLIM WHERE X = ' // XERN3 //
     *      ' Y = ' // XERN4 // ' Z = '  // XERN5 // ' P = ' // XERN6 //
     *      ' AND LOLIM = ', 2, 1)
         RETURN
      ENDIF
C
      IER = 0
      XN = X
      YN = Y
      ZN = Z
      PN = P
      SIGMA = 0.0D0
      POWER4 = 1.0D0
C
   30 MU = (XN+YN+ZN+PN+PN)*0.20D0
      XNDEV = (MU-XN)/MU
      YNDEV = (MU-YN)/MU
      ZNDEV = (MU-ZN)/MU
      PNDEV = (MU-PN)/MU
      EPSLON = MAX(ABS(XNDEV), ABS(YNDEV), ABS(ZNDEV), ABS(PNDEV))
      IF (EPSLON.LT.ERRTOL) GO TO 40
      XNROOT =  SQRT(XN)
      YNROOT =  SQRT(YN)
      ZNROOT =  SQRT(ZN)
      LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
      ALFA = PN*(XNROOT+YNROOT+ZNROOT) + XNROOT*YNROOT*ZNROOT
      ALFA = ALFA*ALFA
      BETA = PN*(PN+LAMDA)*(PN+LAMDA)
      SIGMA = SIGMA + POWER4*DRC(ALFA,BETA,IER)
      POWER4 = POWER4*0.250D0
      XN = (XN+LAMDA)*0.250D0
      YN = (YN+LAMDA)*0.250D0
      ZN = (ZN+LAMDA)*0.250D0
      PN = (PN+LAMDA)*0.250D0
      GO TO 30
C
   40 EA = XNDEV*(YNDEV+ZNDEV) + YNDEV*ZNDEV
      EB = XNDEV*YNDEV*ZNDEV
      EC = PNDEV*PNDEV
      E2 = EA - 3.0D0*EC
      E3 = EB + 2.0D0*PNDEV*(EA-EC)
      S1 = 1.0D0 + E2*(-C1+0.750D0*C3*E2-1.50D0*C4*E3)
      S2 = EB*(0.50D0*C2+PNDEV*(-C3-C3+PNDEV*C4))
      S3 = PNDEV*EA*(C2-PNDEV*C3) - C2*PNDEV*EC
      DRJ = 3.0D0*SIGMA + POWER4*(S1+S2+S3)/(MU* SQRT(MU))
      RETURN
      END
