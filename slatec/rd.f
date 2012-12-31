*DECK RD
      REAL FUNCTION RD (X, Y, Z, IER)
C***BEGIN PROLOGUE  RD
C***PURPOSE  Compute the incomplete or complete elliptic integral of the
C            2nd kind.  For X and Y nonnegative, X+Y and Z positive,
C             RD(X,Y,Z) = Integral from zero to infinity of
C                                -1/2     -1/2     -3/2
C                      (3/2)(t+X)    (t+Y)    (t+Z)    dt.
C            If X or Y is zero, the integral is complete.
C***LIBRARY   SLATEC
C***CATEGORY  C14
C***TYPE      SINGLE PRECISION (RD-S, DRD-D)
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
C   1.     RD
C          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
C          of the second kind
C          Standard FORTRAN function routine
C          Single precision version
C          The routine calculates an approximation result to
C          RD(X,Y,Z) = Integral from zero to infinity of
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
C          RD( X, Y, Z, IER )
C
C          Parameters on Entry
C          Values assigned by the calling routine
C
C          X      - Single precision, nonnegative variable
C
C          Y      - Single precision, nonnegative variable
C
C                   X + Y is positive
C
C          Z      - Real, positive variable
C
C
C
C          On Return     (values assigned by the RD routine)
C
C          RD     - Real approximation to the integral
C
C
C          IER    - Integer
C
C                   IER = 0 Normal and reliable termination of the
C                           routine.  It is assumed that the requested
C                           accuracy has been achieved.
C
C                   IER >  0 Abnormal termination of the routine
C
C
C          X, Y, Z are unaltered.
C
C   3.    Error Messages
C
C         Value of IER assigned by the RD routine
C
C                  Value Assigned         Error Message Printed
C                  IER = 1                MIN(X,Y) .LT. 0.0E0
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
C                    Not greater than (0.1E0 * ERRTOL / machine
C                    minimum) ** (2/3), where ERRTOL is described below.
C                    In the following table it is assumed that ERRTOL
C                    will never be chosen smaller than 1.0E-5.
C
C
C                    Acceptable Values For:   LOLIM      UPLIM
C                    IBM 360/370 SERIES   :   6.0E-51     1.0E+48
C                    CDC 6000/7000 SERIES :   5.0E-215    2.0E+191
C                    UNIVAC 1100 SERIES   :   1.0E-25     2.0E+21
C                    CRAY                 :   3.0E-1644   1.69E+1640
C                    VAX 11 SERIES        :   1.0E-25     4.5E+21
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
C              The accuracy of the computed approximation to the inte-
C              gral can be controlled by choosing the value of ERRTOL.
C              Truncation of a Taylor series after terms of fifth order
C              introduces an error less than the amount shown in the
C              second column of the following table for each value of
C              ERRTOL in the first column.  In addition to the trunca-
C              tion error there will be round-off error, but in prac-
C              tice the total error from both sources is usually less
C              than the amount given in the table.
C
C
C
C
C          Sample Choices:  ERRTOL   Relative Truncation
C                                    error less than
C                           1.0E-3    4.0E-18
C                           3.0E-3    3.0E-15
C                           1.0E-2    4.0E-12
C                           3.0E-2    3.0E-9
C                           1.0E-1    4.0E-6
C
C
C                    Decreasing ERRTOL by a factor of 10 yields six more
C                    decimal digits of accuracy at the expense of one or
C                    two more iterations of the duplication theorem.
C
C *Long Description:
C
C   RD Special Comments
C
C
C
C          Check: RD(X,Y,Z) + RD(Y,Z,X) + RD(Z,X,Y)
C          = 3 /  SQRT(X * Y * Z), where X, Y, and Z are positive.
C
C
C          On Input:
C
C          X, Y, and Z are the variables in the integral RD(X,Y,Z).
C
C
C          On Output:
C
C
C          X, Y, and Z are unaltered.
C
C
C
C          ********************************************************
C
C           WARNING: Changes in the program may improve speed at the
C                    expense of robustness.
C
C
C
C    -------------------------------------------------------------------
C
C
C   Special Functions via RD and RF
C
C
C                  Legendre form of ELLIPTIC INTEGRAL of 2nd kind
C                  ----------------------------------------------
C
C
C                                            2         2   2
C                  E(PHI,K) = SIN(PHI) RF(COS (PHI),1-K SIN (PHI),1) -
C
C                     2      3            2         2   2
C                  -(K/3) SIN (PHI) RD(COS (PHI),1-K SIN (PHI),1)
C
C
C                                 2        2           2
C                  E(K) = RF(0,1-K ,1) - (K/3) RD(0,1-K ,1)
C
C
C                         PI/2     2   2      1/2
C                       = INT  (1-K SIN (PHI) )  D PHI
C                          0
C
C
C
C                  Bulirsch form of ELLIPTIC INTEGRAL of 2nd kind
C                  ----------------------------------------------
C
C                                              2 2    2
C                  EL2(X,KC,A,B) = AX RF(1,1+KC X ,1+X ) +
C
C                                              3         2 2    2
C                                 +(1/3)(B-A) X RD(1,1+KC X ,1+X )
C
C
C
C                  Legendre form of alternative ELLIPTIC INTEGRAL of 2nd
C                  -----------------------------------------------------
C                        kind
C                        ----
C
C                            Q     2       2   2  -1/2
C                  D(Q,K) = INT SIN P  (1-K SIN P)     DP
C                            0
C
C
C
C                                   3          2     2   2
C                  D(Q,K) =(1/3)(SIN Q)  RD(COS Q,1-K SIN Q,1)
C
C
C
C
C
C                  Lemniscate constant B
C                  ---------------------
C
C
C
C                       1    2    4 -1/2
C                  B = INT  S (1-S )    DS
C                       0
C
C
C                  B =(1/3)RD (0,2,1)
C
C
C
C
C                  Heuman's LAMBDA function
C                  ------------------------
C
C
C
C                  (PI/2) LAMBDA0(A,B) =
C
C                                       2                2
C                     = SIN(B) (RF(0,COS (A),1)-(1/3) SIN (A) *
C
C                               2              2         2       2
C                      *RD(0,COS (A),1)) RF(COS (B),1-COS (A) SIN (B),1)
C
C                               2       3            2
C                     -(1/3) COS (A) SIN (B) RF(0,COS (A),1) *
C
C                             2         2       2
C                      *RD(COS (B),1-COS (A) SIN (B),1)
C
C
C
C                  Jacobi ZETA function
C                  --------------------
C
C
C                             2                2       2   2
C                  Z(B,K) = (K/3) SIN(B) RF(COS (B),1-K SIN (B),1)
C
C
C                                      2            2
C                             *RD(0,1-K ,1)/RF(0,1-K ,1)
C
C                               2       3          2       2   2
C                            -(K /3) SIN (B) RD(COS (B),1-K SIN (B),1)
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
C***ROUTINES CALLED  R1MACH, XERMSG
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
C***END PROLOGUE  RD
      CHARACTER*16 XERN3, XERN4, XERN5, XERN6
      INTEGER IER
      REAL LOLIM, UPLIM, EPSLON, ERRTOL
      REAL C1, C2, C3, C4, EA, EB, EC, ED, EF, LAMDA
      REAL MU, POWER4, SIGMA, S1, S2, X, XN, XNDEV
      REAL XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV, ZNROOT
      LOGICAL FIRST
      SAVE ERRTOL, LOLIM, UPLIM, C1, C2, C3, C4, FIRST
      DATA FIRST /.TRUE./
C
C***FIRST EXECUTABLE STATEMENT  RD
      IF (FIRST) THEN
         ERRTOL = (R1MACH(3)/3.0E0)**(1.0E0/6.0E0)
         LOLIM  = 2.0E0/(R1MACH(2))**(2.0E0/3.0E0)
         TUPLIM = R1MACH(1)**(1.0E0/3.0E0)
         TUPLIM = (0.10E0*ERRTOL)**(1.0E0/3.0E0)/TUPLIM
         UPLIM  = TUPLIM**2.0E0
C
         C1 = 3.0E0/14.0E0
         C2 = 1.0E0/6.0E0
         C3 = 9.0E0/22.0E0
         C4 = 3.0E0/26.0E0
      ENDIF
      FIRST = .FALSE.
C
C         CALL ERROR HANDLER IF NECESSARY.
C
      RD = 0.0E0
      IF( MIN(X,Y).LT.0.0E0) THEN
         IER = 1
         WRITE (XERN3, '(1PE15.6)') X
         WRITE (XERN4, '(1PE15.6)') Y
         CALL XERMSG ('SLATEC', 'RD',
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
         CALL XERMSG ('SLATEC', 'RD',
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
         CALL XERMSG ('SLATEC', 'RD',
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
      SIGMA = 0.0E0
      POWER4 = 1.0E0
C
   30 MU = (XN+YN+3.0E0*ZN)*0.20E0
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
      POWER4 = POWER4*0.250E0
      XN = (XN+LAMDA)*0.250E0
      YN = (YN+LAMDA)*0.250E0
      ZN = (ZN+LAMDA)*0.250E0
      GO TO 30
C
   40 EA = XNDEV*YNDEV
      EB = ZNDEV*ZNDEV
      EC = EA - EB
      ED = EA - 6.0E0*EB
      EF = ED + EC + EC
      S1 = ED*(-C1+0.250E0*C3*ED-1.50E0*C4*ZNDEV*EF)
      S2 = ZNDEV*(C2*EF+ZNDEV*(-C3*EC+ZNDEV*C4*EA))
      RD = 3.0E0*SIGMA + POWER4*(1.0E0+S1+S2)/(MU* SQRT(MU))
C
      RETURN
      END
