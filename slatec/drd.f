*DECK DRD
      DOUBLE PRECISION FUNCTION DRD (X, Y, Z, IER)
C***BEGIN PROLOGUE  DRD
C***PURPOSE  Compute the incomplete or complete elliptic integral of
C            the 2nd kind. For X and Y nonnegative, X+Y and Z positive,
C            DRD(X,Y,Z) = Integral from zero to infinity of
C                                -1/2     -1/2     -3/2
C                      (3/2)(t+X)    (t+Y)    (t+Z)    dt.
C            If X or Y is zero, the integral is complete.
C***LIBRARY   SLATEC
C***CATEGORY  C14
C***TYPE      DOUBLE PRECISION (RD-S, DRD-D)
C***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
C             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE SECOND KIND,
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
C   1.     DRD
C          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
C          of the second kind
C          Standard FORTRAN function routine
C          Double precision version
C          The routine calculates an approximation result to
C          DRD(X,Y,Z) = Integral from zero to infinity of
C                              -1/2     -1/2     -3/2
C                    (3/2)(t+X)    (t+Y)    (t+Z)    dt,
C          where X and Y are nonnegative, X + Y is positive, and Z is
C          positive.  If X or Y is zero, the integral is COMPLETE.
C          The duplication theorem is iterated until the variables are
C          nearly equal, and the function is then expanded in Taylor
C          series to fifth order.
C
C   2.     Calling Sequence
C
C          DRD( X, Y, Z, IER )
C
C          Parameters On Entry
C          Values assigned by the calling routine
C
C          X      - Double precision, nonnegative variable
C
C          Y      - Double precision, nonnegative variable
C
C                   X + Y is positive
C
C          Z      - Double precision, positive variable
C
C
C
C          On Return    (values assigned by the DRD routine)
C
C          DRD     - Double precision approximation to the integral
C
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
C          X, Y, Z are unaltered.
C
C   3.    Error Messages
C
C         Value of IER assigned by the DRD routine
C
C                  Value assigned         Error message printed
C                  IER = 1                MIN(X,Y) .LT. 0.0D0
C                      = 2                MIN(X + Y, Z ) .LT. LOLIM
C                      = 3                MAX(X,Y,Z) .GT. UPLIM
C
C
C   4.     Control Parameters
C
C                  Values of LOLIM, UPLIM, and ERRTOL are set by the
C                  routine.
C
C          LOLIM and UPLIM determine the valid range of X, Y, and Z
C
C          LOLIM  - Lower limit of valid arguments
C
C                    Not less  than 2 / (machine maximum) ** (2/3).
C
C          UPLIM  - Upper limit of valid arguments
C
C                 Not greater than (0.1D0 * ERRTOL / machine
C                 minimum) ** (2/3), where ERRTOL is described below.
C                 In the following table it is assumed that ERRTOL will
C                 never be chosen smaller than 1.0D-5.
C
C
C                    Acceptable values for:   LOLIM      UPLIM
C                    IBM 360/370 SERIES   :   6.0D-51     1.0D+48
C                    CDC 6000/7000 SERIES :   5.0D-215    2.0D+191
C                    UNIVAC 1100 SERIES   :   1.0D-205    2.0D+201
C                    CRAY                 :   3.0D-1644   1.69D+1640
C                    VAX 11 SERIES        :   1.0D-25     4.5D+21
C
C
C          ERRTOL determines the accuracy of the answer
C
C                 The value assigned by the routine will result
C                 in solution precision within 1-2 decimals of
C                 "machine precision".
C
C          ERRTOL    Relative error due to truncation is less than
C                    3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2.
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
C          Sample choices:  ERRTOL   Relative truncation
C                                    error less than
C                           1.0D-3    4.0D-18
C                           3.0D-3    3.0D-15
C                           1.0D-2    4.0D-12
C                           3.0D-2    3.0D-9
C                           1.0D-1    4.0D-6
C
C
C                    Decreasing ERRTOL by a factor of 10 yields six more
C                    decimal digits of accuracy at the expense of one or
C                    two more iterations of the duplication theorem.
C
C *Long Description:
C
C   DRD Special Comments
C
C
C
C          Check: DRD(X,Y,Z) + DRD(Y,Z,X) + DRD(Z,X,Y)
C          = 3 / SQRT(X * Y * Z), where X, Y, and Z are positive.
C
C
C          On Input:
C
C          X, Y, and Z are the variables in the integral DRD(X,Y,Z).
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
C    -------------------------------------------------------------------
C
C
C   Special double precision functions via DRD and DRF
C
C
C                  Legendre form of ELLIPTIC INTEGRAL of 2nd kind
C
C                  -----------------------------------------
C
C
C                                             2         2   2
C                  E(PHI,K) = SIN(PHI) DRF(COS (PHI),1-K SIN (PHI),1) -
C
C                     2      3             2         2   2
C                  -(K/3) SIN (PHI) DRD(COS (PHI),1-K SIN (PHI),1)
C
C
C                                  2        2            2
C                  E(K) = DRF(0,1-K ,1) - (K/3) DRD(0,1-K ,1)
C
C                         PI/2     2   2      1/2
C                       = INT  (1-K SIN (PHI) )  D PHI
C                          0
C
C                  Bulirsch form of ELLIPTIC INTEGRAL of 2nd kind
C
C                  -----------------------------------------
C
C                                               2 2    2
C                  EL2(X,KC,A,B) = AX DRF(1,1+KC X ,1+X ) +
C
C                                              3          2 2    2
C                                 +(1/3)(B-A) X DRD(1,1+KC X ,1+X )
C
C
C
C
C                  Legendre form of alternative ELLIPTIC INTEGRAL
C                  of 2nd kind
C
C                  -----------------------------------------
C
C
C
C                            Q     2       2   2  -1/2
C                  D(Q,K) = INT SIN P  (1-K SIN P)     DP
C                            0
C
C
C
C                                     3          2     2   2
C                  D(Q,K) = (1/3) (SIN Q) DRD(COS Q,1-K SIN Q,1)
C
C
C
C
C                  Lemniscate constant  B
C
C                  -----------------------------------------
C
C
C
C
C                       1    2    4 -1/2
C                  B = INT  S (1-S )    DS
C                       0
C
C
C                  B = (1/3) DRD (0,2,1)
C
C
C                  Heuman's LAMBDA function
C
C                  -----------------------------------------
C
C
C
C                  (PI/2) LAMBDA0(A,B) =
C
C                                    2                2
C                 = SIN(B) (DRF(0,COS (A),1)-(1/3) SIN (A) *
C
C                            2               2         2       2
C                  *DRD(0,COS (A),1)) DRF(COS (B),1-COS (A) SIN (B),1)
C
C                            2       3             2
C                  -(1/3) COS (A) SIN (B) DRF(0,COS (A),1) *
C
C                           2         2       2
C                   *DRD(COS (B),1-COS (A) SIN (B),1)
C
C
C
C                  Jacobi ZETA function
C
C                  -----------------------------------------
C
C                             2                 2       2   2
C                  Z(B,K) = (K/3) SIN(B) DRF(COS (B),1-K SIN (B),1)
C
C
C                                       2             2
C                             *DRD(0,1-K ,1)/DRF(0,1-K ,1)
C
C                               2       3           2       2   2
C                            -(K /3) SIN (B) DRD(COS (B),1-K SIN (B),1)
C
C
C ---------------------------------------------------------------------
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
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900510  Modify calls to XERMSG to put in standard form.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DRD
      CHARACTER*16 XERN3, XERN4, XERN5, XERN6
      INTEGER IER
      DOUBLE PRECISION LOLIM, TUPLIM, UPLIM, EPSLON, ERRTOL, D1MACH
      DOUBLE PRECISION C1, C2, C3, C4, EA, EB, EC, ED, EF, LAMDA
      DOUBLE PRECISION MU, POWER4, SIGMA, S1, S2, X, XN, XNDEV
      DOUBLE PRECISION XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV,
     * ZNROOT
      LOGICAL FIRST
      SAVE ERRTOL, LOLIM, UPLIM, C1, C2, C3, C4, FIRST
      DATA FIRST /.TRUE./
C
C***FIRST EXECUTABLE STATEMENT  DRD
      IF (FIRST) THEN
         ERRTOL = (D1MACH(3)/3.0D0)**(1.0D0/6.0D0)
         LOLIM  = 2.0D0/(D1MACH(2))**(2.0D0/3.0D0)
         TUPLIM = D1MACH(1)**(1.0E0/3.0E0)
         TUPLIM = (0.10D0*ERRTOL)**(1.0E0/3.0E0)/TUPLIM
         UPLIM  = TUPLIM**2.0D0
C
         C1 = 3.0D0/14.0D0
         C2 = 1.0D0/6.0D0
         C3 = 9.0D0/22.0D0
         C4 = 3.0D0/26.0D0
      ENDIF
      FIRST = .FALSE.
C
C         CALL ERROR HANDLER IF NECESSARY.
C
      DRD = 0.0D0
      IF( MIN(X,Y).LT.0.0D0) THEN
         IER = 1
         WRITE (XERN3, '(1PE15.6)') X
         WRITE (XERN4, '(1PE15.6)') Y
         CALL XERMSG ('SLATEC', 'DRD',
     *      'MIN(X,Y).LT.0 WHERE X = ' // XERN3 // ' AND Y = ' //
     *      XERN4, 1, 1)
         RETURN
      ENDIF
C
      IF (MAX(X,Y,Z).GT.UPLIM) THEN
         IER = 3
         WRITE (XERN3, '(1PE15.6)') X
         WRITE (XERN4, '(1PE15.6)') Y
         WRITE (XERN5, '(1PE15.6)') Z
         WRITE (XERN6, '(1PE15.6)') UPLIM
         CALL XERMSG ('SLATEC', 'DRD',
     *      'MAX(X,Y,Z).GT.UPLIM WHERE X = ' // XERN3 // ' Y = ' //
     *      XERN4 // ' Z = ' // XERN5 // ' AND UPLIM = ' // XERN6,
     *      3, 1)
         RETURN
      ENDIF
C
      IF (MIN(X+Y,Z).LT.LOLIM) THEN
         IER = 2
         WRITE (XERN3, '(1PE15.6)') X
         WRITE (XERN4, '(1PE15.6)') Y
         WRITE (XERN5, '(1PE15.6)') Z
         WRITE (XERN6, '(1PE15.6)') LOLIM
         CALL XERMSG ('SLATEC', 'DRD',
     *      'MIN(X+Y,Z).LT.LOLIM WHERE X = ' // XERN3 // ' Y = ' //
     *      XERN4 // ' Z = ' // XERN5 // ' AND LOLIM = ' // XERN6,
     *      2, 1)
         RETURN
      ENDIF
C
      IER = 0
      XN = X
      YN = Y
      ZN = Z
      SIGMA = 0.0D0
      POWER4 = 1.0D0
C
   30 MU = (XN+YN+3.0D0*ZN)*0.20D0
      XNDEV = (MU-XN)/MU
      YNDEV = (MU-YN)/MU
      ZNDEV = (MU-ZN)/MU
      EPSLON = MAX(ABS(XNDEV), ABS(YNDEV), ABS(ZNDEV))
      IF (EPSLON.LT.ERRTOL) GO TO 40
      XNROOT = SQRT(XN)
      YNROOT = SQRT(YN)
      ZNROOT = SQRT(ZN)
      LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
      SIGMA = SIGMA + POWER4/(ZNROOT*(ZN+LAMDA))
      POWER4 = POWER4*0.250D0
      XN = (XN+LAMDA)*0.250D0
      YN = (YN+LAMDA)*0.250D0
      ZN = (ZN+LAMDA)*0.250D0
      GO TO 30
C
   40 EA = XNDEV*YNDEV
      EB = ZNDEV*ZNDEV
      EC = EA - EB
      ED = EA - 6.0D0*EB
      EF = ED + EC + EC
      S1 = ED*(-C1+0.250D0*C3*ED-1.50D0*C4*ZNDEV*EF)
      S2 = ZNDEV*(C2*EF+ZNDEV*(-C3*EC+ZNDEV*C4*EA))
      DRD = 3.0D0*SIGMA + POWER4*(1.0D0+S1+S2)/(MU*SQRT(MU))
C
      RETURN
      END
