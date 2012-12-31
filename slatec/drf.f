*DECK DRF
      DOUBLE PRECISION FUNCTION DRF (X, Y, Z, IER)
C***BEGIN PROLOGUE  DRF
C***PURPOSE  Compute the incomplete or complete elliptic integral of the
C            1st kind.  For X, Y, and Z non-negative and at most one of
C            them zero, RF(X,Y,Z) = Integral from zero to infinity of
C                                -1/2     -1/2     -1/2
C                      (1/2)(t+X)    (t+Y)    (t+Z)    dt.
C            If X, Y or Z is zero, the integral is complete.
C***LIBRARY   SLATEC
C***CATEGORY  C14
C***TYPE      DOUBLE PRECISION (RF-S, DRF-D)
C***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
C             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE FIRST KIND,
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
C   1.     DRF
C          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
C          of the first kind
C          Standard FORTRAN function routine
C          Double precision version
C          The routine calculates an approximation result to
C          DRF(X,Y,Z) = Integral from zero to infinity of
C
C                               -1/2     -1/2     -1/2
C                     (1/2)(t+X)    (t+Y)    (t+Z)    dt,
C
C          where X, Y, and Z are nonnegative and at most one of them
C          is zero.  If one of them  is zero, the integral is COMPLETE.
C          The duplication theorem is iterated until the variables are
C          nearly equal, and the function is then expanded in Taylor
C          series to fifth order.
C
C   2.     Calling sequence
C          DRF( X, Y, Z, IER )
C
C          Parameters On entry
C          Values assigned by the calling routine
C
C          X      - Double precision, nonnegative variable
C
C          Y      - Double precision, nonnegative variable
C
C          Z      - Double precision, nonnegative variable
C
C
C
C          On Return    (values assigned by the DRF routine)
C
C          DRF     - Double precision approximation to the integral
C
C          IER    - Integer
C
C                   IER = 0 Normal and reliable termination of the
C                           routine. It is assumed that the requested
C                           accuracy has been achieved.
C
C                   IER >  0 Abnormal termination of the routine
C
C          X, Y, Z are unaltered.
C
C
C   3.    Error Messages
C
C
C         Value of IER assigned by the DRF routine
C
C                  Value assigned         Error Message Printed
C                  IER = 1                MIN(X,Y,Z) .LT. 0.0D0
C                      = 2                MIN(X+Y,X+Z,Y+Z) .LT. LOLIM
C                      = 3                MAX(X,Y,Z) .GT. UPLIM
C
C
C
C   4.     Control Parameters
C
C                  Values of LOLIM, UPLIM, and ERRTOL are set by the
C                  routine.
C
C          LOLIM and UPLIM determine the valid range of X, Y and Z
C
C          LOLIM  - Lower limit of valid arguments
C
C                   Not less than 5 * (machine minimum).
C
C          UPLIM  - Upper limit of valid arguments
C
C                   Not greater than (machine maximum) / 5.
C
C
C                     Acceptable values for:   LOLIM      UPLIM
C                     IBM 360/370 SERIES   :   3.0D-78     1.0D+75
C                     CDC 6000/7000 SERIES :   1.0D-292    1.0D+321
C                     UNIVAC 1100 SERIES   :   1.0D-307    1.0D+307
C                     CRAY                 :   2.3D-2466   1.09D+2465
C                     VAX 11 SERIES        :   1.5D-38     3.0D+37
C
C
C
C          ERRTOL determines the accuracy of the answer
C
C                 The value assigned by the routine will result
C                 in solution precision within 1-2 decimals of
C                 "machine precision".
C
C
C
C          ERRTOL - Relative error due to truncation is less than
C                   ERRTOL ** 6 / (4 * (1-ERRTOL)  .
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
C
C
C          Sample choices:  ERRTOL   Relative Truncation
C                                    error less than
C                           1.0D-3    3.0D-19
C                           3.0D-3    2.0D-16
C                           1.0D-2    3.0D-13
C                           3.0D-2    2.0D-10
C                           1.0D-1    3.0D-7
C
C
C                    Decreasing ERRTOL by a factor of 10 yields six more
C                    decimal digits of accuracy at the expense of one or
C                    two more iterations of the duplication theorem.
C
C *Long Description:
C
C   DRF Special Comments
C
C
C
C          Check by addition theorem: DRF(X,X+Z,X+W) + DRF(Y,Y+Z,Y+W)
C          = DRF(0,Z,W), where X,Y,Z,W are positive and X * Y = Z * W.
C
C
C          On Input:
C
C          X, Y, and Z are the variables in the integral DRF(X,Y,Z).
C
C
C          On Output:
C
C
C          X, Y, Z are unaltered.
C
C
C
C          ********************************************************
C
C          WARNING: Changes in the program may improve speed at the
C                   expense of robustness.
C
C
C
C   Special double precision functions via DRF
C
C
C
C
C                  Legendre form of ELLIPTIC INTEGRAL of 1st kind
C
C                  -----------------------------------------
C
C
C
C                                             2         2   2
C                  F(PHI,K) = SIN(PHI) DRF(COS (PHI),1-K SIN (PHI),1)
C
C
C                                  2
C                  K(K) = DRF(0,1-K ,1)
C
C
C                         PI/2     2   2      -1/2
C                       = INT  (1-K SIN (PHI) )   D PHI
C                          0
C
C
C
C                  Bulirsch form of ELLIPTIC INTEGRAL of 1st kind
C
C                  -----------------------------------------
C
C
C                                          2 2    2
C                  EL1(X,KC) = X DRF(1,1+KC X ,1+X )
C
C
C                  Lemniscate constant A
C
C                  -----------------------------------------
C
C
C                       1      4 -1/2
C                  A = INT (1-S )    DS = DRF(0,1,2) = DRF(0,2,1)
C                       0
C
C
C
C    -------------------------------------------------------------------
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
C***ROUTINES CALLED  D1MACH, XERMSG
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
C           editorial changes.  (RWC))
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DRF
      CHARACTER*16 XERN3, XERN4, XERN5, XERN6
      INTEGER IER
      DOUBLE PRECISION LOLIM, UPLIM, EPSLON, ERRTOL, D1MACH
      DOUBLE PRECISION C1, C2, C3, E2, E3, LAMDA
      DOUBLE PRECISION MU, S, X, XN, XNDEV
      DOUBLE PRECISION XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV,
     * ZNROOT
      LOGICAL FIRST
      SAVE ERRTOL,LOLIM,UPLIM,C1,C2,C3,FIRST
      DATA FIRST /.TRUE./
C
C***FIRST EXECUTABLE STATEMENT  DRF
C
      IF (FIRST) THEN
         ERRTOL = (4.0D0*D1MACH(3))**(1.0D0/6.0D0)
         LOLIM  = 5.0D0 * D1MACH(1)
         UPLIM  = D1MACH(2)/5.0D0
C
         C1 = 1.0D0/24.0D0
         C2 = 3.0D0/44.0D0
         C3 = 1.0D0/14.0D0
      ENDIF
      FIRST = .FALSE.
C
C         CALL ERROR HANDLER IF NECESSARY.
C
      DRF = 0.0D0
      IF (MIN(X,Y,Z).LT.0.0D0) THEN
         IER = 1
         WRITE (XERN3, '(1PE15.6)') X
         WRITE (XERN4, '(1PE15.6)') Y
         WRITE (XERN5, '(1PE15.6)') Z
         CALL XERMSG ('SLATEC', 'DRF',
     *      'MIN(X,Y,Z).LT.0 WHERE X = ' // XERN3 // ' Y = ' // XERN4 //
     *      ' AND Z = ' // XERN5, 1, 1)
         RETURN
      ENDIF
C
      IF (MAX(X,Y,Z).GT.UPLIM) THEN
         IER = 3
         WRITE (XERN3, '(1PE15.6)') X
         WRITE (XERN4, '(1PE15.6)') Y
         WRITE (XERN5, '(1PE15.6)') Z
         WRITE (XERN6, '(1PE15.6)') UPLIM
         CALL XERMSG ('SLATEC', 'DRF',
     *      'MAX(X,Y,Z).GT.UPLIM WHERE X = '  // XERN3 // ' Y = ' //
     *      XERN4 // ' Z = ' // XERN5 // ' AND UPLIM = ' // XERN6, 3, 1)
         RETURN
      ENDIF
C
      IF (MIN(X+Y,X+Z,Y+Z).LT.LOLIM) THEN
         IER = 2
         WRITE (XERN3, '(1PE15.6)') X
         WRITE (XERN4, '(1PE15.6)') Y
         WRITE (XERN5, '(1PE15.6)') Z
         WRITE (XERN6, '(1PE15.6)') LOLIM
         CALL XERMSG ('SLATEC', 'DRF',
     *      'MIN(X+Y,X+Z,Y+Z).LT.LOLIM WHERE X = ' // XERN3 //
     *      ' Y = ' // XERN4 // ' Z = ' // XERN5 // ' AND LOLIM = ' //
     *      XERN6, 2, 1)
         RETURN
      ENDIF
C
      IER = 0
      XN = X
      YN = Y
      ZN = Z
C
   30 MU = (XN+YN+ZN)/3.0D0
      XNDEV = 2.0D0 - (MU+XN)/MU
      YNDEV = 2.0D0 - (MU+YN)/MU
      ZNDEV = 2.0D0 - (MU+ZN)/MU
      EPSLON = MAX(ABS(XNDEV),ABS(YNDEV),ABS(ZNDEV))
      IF (EPSLON.LT.ERRTOL) GO TO 40
      XNROOT = SQRT(XN)
      YNROOT = SQRT(YN)
      ZNROOT = SQRT(ZN)
      LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
      XN = (XN+LAMDA)*0.250D0
      YN = (YN+LAMDA)*0.250D0
      ZN = (ZN+LAMDA)*0.250D0
      GO TO 30
C
   40 E2 = XNDEV*YNDEV - ZNDEV*ZNDEV
      E3 = XNDEV*YNDEV*ZNDEV
      S  = 1.0D0 + (C1*E2-0.10D0-C2*E3)*E2 + C3*E3
      DRF = S/SQRT(MU)
C
      RETURN
      END
