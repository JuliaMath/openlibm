*DECK CV
      REAL FUNCTION CV (XVAL, NDATA, NCONST, NORD, NBKPT, BKPT, W)
C***BEGIN PROLOGUE  CV
C***PURPOSE  Evaluate the variance function of the curve obtained
C            by the constrained B-spline fitting subprogram FC.
C***LIBRARY   SLATEC
C***CATEGORY  L7A3
C***TYPE      SINGLE PRECISION (CV-S, DCV-D)
C***KEYWORDS  ANALYSIS OF COVARIANCE, B-SPLINE,
C             CONSTRAINED LEAST SQUARES, CURVE FITTING
C***AUTHOR  Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C     CV( ) is a companion function subprogram for FC( ).  The
C     documentation for FC( ) has complete usage instructions.
C
C     CV( ) is used to evaluate the variance function of the curve
C     obtained by the constrained B-spline fitting subprogram, FC( ).
C     The variance function defines the square of the probable error
C     of the fitted curve at any point, XVAL.  One can use the square
C     root of this variance function to determine a probable error band
C     around the fitted curve.
C
C     CV( ) is used after a call to FC( ).  MODE, an input variable to
C     FC( ), is used to indicate if the variance function is desired.
C     In order to use CV( ), MODE must equal 2 or 4 on input to FC( ).
C     MODE is also used as an output flag from FC( ).  Check to make
C     sure that MODE = 0 after calling FC( ), indicating a successful
C     constrained curve fit.  The array SDDATA, as input to FC( ), must
C     also be defined with the standard deviation or uncertainty of the
C     Y values to use CV( ).
C
C     To evaluate the variance function after calling FC( ) as stated
C     above, use CV( ) as shown here
C
C          VAR=CV(XVAL,NDATA,NCONST,NORD,NBKPT,BKPT,W)
C
C     The variance function is given by
C
C          VAR=(transpose of B(XVAL))*C*B(XVAL)/MAX(NDATA-N,1)
C
C     where N = NBKPT - NORD.
C
C     The vector B(XVAL) is the B-spline basis function values at
C     X=XVAL.  The covariance matrix, C, of the solution coefficients
C     accounts only for the least squares equations and the explicitly
C     stated equality constraints.  This fact must be considered when
C     interpreting the variance function from a data fitting problem
C     that has inequality constraints on the fitted curve.
C
C     All the variables in the calling sequence for CV( ) are used in
C     FC( ) except the variable XVAL.  Do not change the values of these
C     variables between the call to FC( ) and the use of CV( ).
C
C     The following is a brief description of the variables
C
C     XVAL    The point where the variance is desired.
C
C     NDATA   The number of discrete (X,Y) pairs for which FC( )
C             calculated a piece-wise polynomial curve.
C
C     NCONST  The number of conditions that constrained the B-spline in
C             FC( ).
C
C     NORD    The order of the B-spline used in FC( ).
C             The value of NORD must satisfy 1 < NORD < 20 .
C
C             (The order of the spline is one more than the degree of
C             the piece-wise polynomial defined on each interval.  This
C             is consistent with the B-spline package convention.  For
C             example, NORD=4 when we are using piece-wise cubics.)
C
C     NBKPT   The number of knots in the array BKPT(*).
C             The value of NBKPT must satisfy NBKPT .GE. 2*NORD.
C
C     BKPT(*) The real array of knots.  Normally the problem data
C             interval will be included between the limits BKPT(NORD)
C             and BKPT(NBKPT-NORD+1).  The additional end knots
C             BKPT(I),I=1,...,NORD-1 and I=NBKPT-NORD+2,...,NBKPT, are
C             required by FC( ) to compute the functions used to fit
C             the data.
C
C     W(*)    Real work array as used in FC( ).  See FC( ) for the
C             required length of W(*).  The contents of W(*) must not
C             be modified by the user if the variance function is
C             desired.
C
C***REFERENCES  R. J. Hanson, Constrained least squares curve fitting
C                 to discrete data using B-splines, a users guide,
C                 Report SAND78-1291, Sandia Laboratories, December
C                 1978.
C***ROUTINES CALLED  BSPLVN, SDOT
C***REVISION HISTORY  (YYMMDD)
C   780801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CV
      DIMENSION BKPT(NBKPT), W(*), V(40)
C***FIRST EXECUTABLE STATEMENT  CV
      ZERO = 0.
      MDG = NBKPT - NORD + 3
      MDW = NBKPT - NORD + 1 + NCONST
      IS = MDG*(NORD+1) + 2*MAX(NDATA,NBKPT) + NBKPT + NORD**2
      LAST = NBKPT - NORD + 1
      ILEFT = NORD
   10 IF (.NOT.(XVAL.GE.BKPT(ILEFT+1) .AND. ILEFT.LT.LAST-1)) GO TO 20
      ILEFT = ILEFT + 1
      GO TO 10
   20 CALL BSPLVN(BKPT, NORD, 1, XVAL, ILEFT, V(NORD+1))
      ILEFT = ILEFT - NORD + 1
      IP = MDW*(ILEFT-1) + ILEFT + IS
      N = NBKPT - NORD
      DO 30 I=1,NORD
         V(I) = SDOT(NORD,W(IP),1,V(NORD+1),1)
         IP = IP + MDW
   30 CONTINUE
      CV = MAX(SDOT(NORD,V,1,V(NORD+1),1),ZERO)
C
C     SCALE THE VARIANCE SO IT IS AN UNBIASED ESTIMATE.
      CV = CV/MAX(NDATA-N,1)
      RETURN
      END
