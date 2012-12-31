*DECK DRC
      DOUBLE PRECISION FUNCTION DRC (X, Y, IER)
C***BEGIN PROLOGUE  DRC
C***PURPOSE  Calculate a double precision approximation to
C             DRC(X,Y) = Integral from zero to infinity of
C                              -1/2     -1
C                    (1/2)(t+X)    (t+Y)  dt,
C            where X is nonnegative and Y is positive.
C***LIBRARY   SLATEC
C***CATEGORY  C14
C***TYPE      DOUBLE PRECISION (RC-S, DRC-D)
C***KEYWORDS  DUPLICATION THEOREM, ELEMENTARY FUNCTIONS,
C             ELLIPTIC INTEGRAL, TAYLOR SERIES
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
C   1.     DRC
C          Standard FORTRAN function routine
C          Double precision version
C          The routine calculates an approximation result to
C          DRC(X,Y) = integral from zero to infinity of
C
C                              -1/2     -1
C                    (1/2)(t+X)    (t+Y)  dt,
C
C          where X is nonnegative and Y is positive.  The duplication
C          theorem is iterated until the variables are nearly equal,
C          and the function is then expanded in Taylor series to fifth
C          order.  Logarithmic, inverse circular, and inverse hyper-
C          bolic functions can be expressed in terms of DRC.
C
C   2.     Calling Sequence
C          DRC( X, Y, IER )
C
C          Parameters On Entry
C          Values assigned by the calling routine
C
C          X      - Double precision, nonnegative variable
C
C          Y      - Double precision, positive variable
C
C
C
C          On Return  (values assigned by the DRC routine)
C
C          DRC    - Double precision approximation to the integral
C
C          IER    - Integer to indicate normal or abnormal termination.
C
C                     IER = 0 Normal and reliable termination of the
C                             routine.  It is assumed that the requested
C                             accuracy has been achieved.
C
C                     IER > 0 Abnormal termination of the routine
C
C          X and Y are unaltered.
C
C   3.    Error messages
C
C         Value of IER assigned by the DRC routine
C
C                  Value assigned         Error message printed
C                  IER = 1                X.LT.0.0D0.OR.Y.LE.0.0D0
C                      = 2                X+Y.LT.LOLIM
C                      = 3                MAX(X,Y) .GT. UPLIM
C
C   4.     Control parameters
C
C                  Values of LOLIM, UPLIM, and ERRTOL are set by the
C                  routine.
C
C          LOLIM and UPLIM determine the valid range of X and Y
C
C          LOLIM  - Lower limit of valid arguments
C
C                   Not less  than 5 * (machine minimum)  .
C
C          UPLIM  - Upper limit of valid arguments
C
C                   Not greater than (machine maximum) / 5 .
C
C
C                     Acceptable values for:   LOLIM       UPLIM
C                     IBM 360/370 SERIES   :   3.0D-78     1.0D+75
C                     CDC 6000/7000 SERIES :   1.0D-292    1.0D+321
C                     UNIVAC 1100 SERIES   :   1.0D-307    1.0D+307
C                     CRAY                 :   2.3D-2466   1.0D+2465
C                     VAX 11 SERIES        :   1.5D-38     3.0D+37
C
C          ERRTOL determines the accuracy of the answer
C
C                 The value assigned by the routine will result
C                 in solution precision within 1-2 decimals of
C                 "machine precision".
C
C
C          ERRTOL  - relative error due to truncation is less than
C                    16 * ERRTOL ** 6 / (1 - 2 * ERRTOL).
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
C          Sample choices:  ERRTOL   Relative truncation
C                                    error less than
C                           1.0D-3    2.0D-17
C                           3.0D-3    2.0D-14
C                           1.0D-2    2.0D-11
C                           3.0D-2    2.0D-8
C                           1.0D-1    2.0D-5
C
C
C                    Decreasing ERRTOL by a factor of 10 yields six more
C                    decimal digits of accuracy at the expense of one or
C                    two more iterations of the duplication theorem.
C
C *Long Description:
C
C   DRC special comments
C
C
C
C
C                  Check: DRC(X,X+Z) + DRC(Y,Y+Z) = DRC(0,Z)
C
C                  where X, Y, and Z are positive and X * Y = Z * Z
C
C
C          On Input:
C
C          X, and Y are the variables in the integral DRC(X,Y).
C
C          On Output:
C
C          X and Y are unaltered.
C
C
C
C                    DRC(0,1/4)=DRC(1/16,1/8)=PI=3.14159...
C
C                    DRC(9/4,2)=LN(2)
C
C
C
C          ********************************************************
C
C          WARNING: Changes in the program may improve speed at the
C                   expense of robustness.
C
C
C   --------------------------------------------------------------------
C
C   Special functions via DRC
C
C
C
C                  LN X                X .GT. 0
C
C                                             2
C                  LN(X) = (X-1) DRC(((1+X)/2)  , X )
C
C
C   --------------------------------------------------------------------
C
C                  ARCSIN X            -1 .LE. X .LE. 1
C
C                                       2
C                  ARCSIN X = X DRC (1-X  ,1 )
C
C   --------------------------------------------------------------------
C
C                  ARCCOS X            0 .LE. X .LE. 1
C
C
C                                     2       2
C                  ARCCOS X = SQRT(1-X ) DRC(X  ,1 )
C
C   --------------------------------------------------------------------
C
C                  ARCTAN X            -INF .LT. X .LT. +INF
C
C                                        2
C                  ARCTAN X = X DRC(1,1+X  )
C
C   --------------------------------------------------------------------
C
C                  ARCCOT X            0 .LE. X .LT. INF
C
C                                  2   2
C                  ARCCOT X = DRC(X  ,X +1 )
C
C   --------------------------------------------------------------------
C
C                  ARCSINH X           -INF .LT. X .LT. +INF
C
C                                       2
C                  ARCSINH X = X DRC(1+X  ,1 )
C
C   --------------------------------------------------------------------
C
C                  ARCCOSH X           X .GE. 1
C
C                                    2         2
C                  ARCCOSH X = SQRT(X -1) DRC(X  ,1 )
C
C   --------------------------------------------------------------------
C
C                  ARCTANH X           -1 .LT. X .LT. 1
C
C                                         2
C                  ARCTANH X = X DRC(1,1-X  )
C
C   --------------------------------------------------------------------
C
C                  ARCCOTH X           X .GT. 1
C
C                                   2   2
C                  ARCCOTH X = DRC(X  ,X -1 )
C
C   --------------------------------------------------------------------
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
C***END PROLOGUE  DRC
      CHARACTER*16 XERN3, XERN4, XERN5
      INTEGER IER
      DOUBLE PRECISION C1, C2, ERRTOL, LAMDA, LOLIM, D1MACH
      DOUBLE PRECISION MU, S, SN, UPLIM, X, XN, Y, YN
      LOGICAL FIRST
      SAVE ERRTOL,LOLIM,UPLIM,C1,C2,FIRST
      DATA FIRST /.TRUE./
C
C***FIRST EXECUTABLE STATEMENT  DRC
      IF (FIRST) THEN
         ERRTOL = (D1MACH(3)/16.0D0)**(1.0D0/6.0D0)
         LOLIM  = 5.0D0 * D1MACH(1)
         UPLIM  = D1MACH(2) / 5.0D0
C
         C1 = 1.0D0/7.0D0
         C2 = 9.0D0/22.0D0
      ENDIF
      FIRST = .FALSE.
C
C         CALL ERROR HANDLER IF NECESSARY.
C
      DRC = 0.0D0
      IF (X.LT.0.0D0.OR.Y.LE.0.0D0) THEN
         IER = 1
         WRITE (XERN3, '(1PE15.6)') X
         WRITE (XERN4, '(1PE15.6)') Y
         CALL XERMSG ('SLATEC', 'DRC',
     *      'X.LT.0 .OR. Y.LE.0 WHERE X = ' // XERN3 // ' AND Y = ' //
     *      XERN4, 1, 1)
         RETURN
      ENDIF
C
      IF (MAX(X,Y).GT.UPLIM) THEN
         IER = 3
         WRITE (XERN3, '(1PE15.6)') X
         WRITE (XERN4, '(1PE15.6)') Y
         WRITE (XERN5, '(1PE15.6)') UPLIM
         CALL XERMSG ('SLATEC', 'DRC',
     *      'MAX(X,Y).GT.UPLIM WHERE X = '  // XERN3 // ' Y = ' //
     *      XERN4 // ' AND UPLIM = ' // XERN5, 3, 1)
         RETURN
      ENDIF
C
      IF (X+Y.LT.LOLIM) THEN
         IER = 2
         WRITE (XERN3, '(1PE15.6)') X
         WRITE (XERN4, '(1PE15.6)') Y
         WRITE (XERN5, '(1PE15.6)') LOLIM
         CALL XERMSG ('SLATEC', 'DRC',
     *      'X+Y.LT.LOLIM WHERE X = ' // XERN3 // ' Y = ' // XERN4 //
     *      ' AND LOLIM = ' // XERN5, 2, 1)
         RETURN
      ENDIF
C
      IER = 0
      XN = X
      YN = Y
C
   30 MU = (XN+YN+YN)/3.0D0
      SN = (YN+MU)/MU - 2.0D0
      IF (ABS(SN).LT.ERRTOL) GO TO 40
      LAMDA = 2.0D0*SQRT(XN)*SQRT(YN) + YN
      XN = (XN+LAMDA)*0.250D0
      YN = (YN+LAMDA)*0.250D0
      GO TO 30
C
   40 S = SN*SN*(0.30D0+SN*(C1+SN*(0.3750D0+SN*C2)))
      DRC = (1.0D0+S)/SQRT(MU)
      RETURN
      END
