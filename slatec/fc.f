*DECK FC
      SUBROUTINE FC (NDATA, XDATA, YDATA, SDDATA, NORD, NBKPT, BKPT,
     +   NCONST, XCONST, YCONST, NDERIV, MODE, COEFF, W, IW)
C***BEGIN PROLOGUE  FC
C***PURPOSE  Fit a piecewise polynomial curve to discrete data.
C            The piecewise polynomials are represented as B-splines.
C            The fitting is done in a weighted least squares sense.
C            Equality and inequality constraints can be imposed on the
C            fitted curve.
C***LIBRARY   SLATEC
C***CATEGORY  K1A1A1, K1A2A, L8A3
C***TYPE      SINGLE PRECISION (FC-S, DFC-D)
C***KEYWORDS  B-SPLINE, CONSTRAINED LEAST SQUARES, CURVE FITTING,
C             WEIGHTED LEAST SQUARES
C***AUTHOR  Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C      This subprogram fits a piecewise polynomial curve
C      to discrete data.  The piecewise polynomials are
C      represented as B-splines.
C      The fitting is done in a weighted least squares sense.
C      Equality and inequality constraints can be imposed on the
C      fitted curve.
C
C      For a description of the B-splines and usage instructions to
C      evaluate them, see
C
C      C. W. de Boor, Package for Calculating with B-Splines.
C                     SIAM J. Numer. Anal., p. 441, (June, 1977).
C
C      For further documentation and discussion of constrained
C      curve fitting using B-splines, see
C
C      R. J. Hanson, Constrained Least Squares Curve Fitting
C                   to Discrete Data Using B-Splines, a User's
C                   Guide. Sandia Labs. Tech. Rept. SAND-78-1291,
C                   December, (1978).
C
C  Input..
C      NDATA,XDATA(*),
C      YDATA(*),
C      SDDATA(*)
C                         The NDATA discrete (X,Y) pairs and the Y value
C                         standard deviation or uncertainty, SD, are in
C                         the respective arrays XDATA(*), YDATA(*), and
C                         SDDATA(*).  No sorting of XDATA(*) is
C                         required.  Any non-negative value of NDATA is
C                         allowed.  A negative value of NDATA is an
C                         error.  A zero value for any entry of
C                         SDDATA(*) will weight that data point as 1.
C                         Otherwise the weight of that data point is
C                         the reciprocal of this entry.
C
C      NORD,NBKPT,
C      BKPT(*)
C                         The NBKPT knots of the B-spline of order NORD
C                         are in the array BKPT(*).  Normally the
C                         problem data interval will be included between
C                         the limits BKPT(NORD) and BKPT(NBKPT-NORD+1).
C                         The additional end knots BKPT(I),I=1,...,
C                         NORD-1 and I=NBKPT-NORD+2,...,NBKPT, are
C                         required to compute the functions used to fit
C                         the data.  No sorting of BKPT(*) is required.
C                         Internal to  FC( ) the extreme end knots may
C                         be reduced and increased respectively to
C                         accommodate any data values that are exterior
C                         to the given knot values.  The contents of
C                         BKPT(*) is not changed.
C
C                         NORD must be in the range 1 .LE. NORD .LE. 20.
C                         The value of NBKPT must satisfy the condition
C                         NBKPT .GE. 2*NORD.
C                         Other values are considered errors.
C
C                         (The order of the spline is one more than the
C                         degree of the piecewise polynomial defined on
C                         each interval.  This is consistent with the
C                         B-spline package convention.  For example,
C                         NORD=4 when we are using piecewise cubics.)
C
C      NCONST,XCONST(*),
C      YCONST(*),NDERIV(*)
C                         The number of conditions that constrain the
C                         B-spline is NCONST.  A constraint is specified
C                         by an (X,Y) pair in the arrays XCONST(*) and
C                         YCONST(*), and by the type of constraint and
C                         derivative value encoded in the array
C                         NDERIV(*).  No sorting of XCONST(*) is
C                         required.  The value of NDERIV(*) is
C                         determined as follows.  Suppose the I-th
C                         constraint applies to the J-th derivative
C                         of the B-spline.  (Any non-negative value of
C                         J < NORD is permitted.  In particular the
C                         value J=0 refers to the B-spline itself.)
C                         For this I-th constraint, set
C                          XCONST(I)=X,
C                          YCONST(I)=Y, and
C                          NDERIV(I)=ITYPE+4*J, where
C
C                          ITYPE = 0,      if (J-th deriv. at X) .LE. Y.
C                                = 1,      if (J-th deriv. at X) .GE. Y.
C                                = 2,      if (J-th deriv. at X) .EQ. Y.
C                                = 3,      if (J-th deriv. at X) .EQ.
C                                             (J-th deriv. at Y).
C                          (A value of NDERIV(I)=-1 will cause this
C                          constraint to be ignored.  This subprogram
C                          feature is often useful when temporarily
C                          suppressing a constraint while still
C                          retaining the source code of the calling
C                          program.)
C
C        MODE
C                         An input flag that directs the least squares
C                         solution method used by FC( ).
C
C                         The variance function, referred to below,
C                         defines the square of the probable error of
C                         the fitted curve at any point, XVAL.
C                         This feature of  FC( ) allows one to use the
C                         square root of this variance function to
C                         determine a probable error band around the
C                         fitted curve.
C
C                         =1  a new problem.  No variance function.
C
C                         =2  a new problem.  Want variance function.
C
C                         =3  an old problem.  No variance function.
C
C                         =4  an old problem.  Want variance function.
C
C                         Any value of MODE other than 1-4 is an error.
C
C                         The user with a new problem can skip directly
C                         to the description of the input parameters
C                         IW(1), IW(2).
C
C                         If the user correctly specifies the new or old
C                         problem status, the subprogram FC( ) will
C                         perform more efficiently.
C                         By an old problem it is meant that subprogram
C                         FC( ) was last called with this same set of
C                         knots, data points and weights.
C
C                         Another often useful deployment of this old
C                         problem designation can occur when one has
C                         previously obtained a Q-R orthogonal
C                         decomposition of the matrix resulting from
C                         B-spline fitting of data (without constraints)
C                         at the breakpoints BKPT(I), I=1,...,NBKPT.
C                         For example, this matrix could be the result
C                         of sequential accumulation of the least
C                         squares equations for a very large data set.
C                         The user writes this code in a manner
C                         convenient for the application.  For the
C                         discussion here let
C
C                                      N=NBKPT-NORD, and K=N+3
C
C                         Let us assume that an equivalent least squares
C                         system
C
C                                      RC=D
C
C                         has been obtained.  Here R is an N+1 by N
C                         matrix and D is a vector with N+1 components.
C                         The last row of R is zero.  The matrix R is
C                         upper triangular and banded.  At most NORD of
C                         the diagonals are nonzero.
C                         The contents of R and D can be copied to the
C                         working array W(*) as follows.
C
C                         The I-th diagonal of R, which has N-I+1
C                         elements, is copied to W(*) starting at
C
C                                      W((I-1)*K+1),
C
C                         for I=1,...,NORD.
C                         The vector D is copied to W(*) starting at
C
C                                      W(NORD*K+1)
C
C                         The input value used for NDATA is arbitrary
C                         when an old problem is designated.  Because
C                         of the feature of FC( ) that checks the
C                         working storage array lengths, a value not
C                         exceeding NBKPT should be used.  For example,
C                         use NDATA=0.
C
C                         (The constraints or variance function request
C                         can change in each call to FC( ).)  A new
C                         problem is anything other than an old problem.
C
C      IW(1),IW(2)
C                         The amounts of working storage actually
C                         allocated for the working arrays W(*) and
C                         IW(*).  These quantities are compared with the
C                         actual amounts of storage needed in FC( ).
C                         Insufficient storage allocated for either
C                         W(*) or IW(*) is an error.  This feature was
C                         included in FC( ) because misreading the
C                         storage formulas for W(*) and IW(*) might very
C                         well lead to subtle and hard-to-find
C                         programming bugs.
C
C                         The length of W(*) must be at least
C
C                           NB=(NBKPT-NORD+3)*(NORD+1)+
C                               2*MAX(NDATA,NBKPT)+NBKPT+NORD**2
C
C                         Whenever possible the code uses banded matrix
C                         processors BNDACC( ) and BNDSOL( ).  These
C                         are utilized if there are no constraints,
C                         no variance function is required, and there
C                         is sufficient data to uniquely determine the
C                         B-spline coefficients.  If the band processors
C                         cannot be used to determine the solution,
C                         then the constrained least squares code LSEI
C                         is used.  In this case the subprogram requires
C                         an additional block of storage in W(*).  For
C                         the discussion here define the integers NEQCON
C                         and NINCON respectively as the number of
C                         equality (ITYPE=2,3) and inequality
C                         (ITYPE=0,1) constraints imposed on the fitted
C                         curve.  Define
C
C                           L=NBKPT-NORD+1
C
C                         and note that
C
C                           NCONST=NEQCON+NINCON.
C
C                         When the subprogram FC( ) uses LSEI( ) the
C                         length of the working array W(*) must be at
C                         least
C
C                           LW=NB+(L+NCONST)*L+
C                              2*(NEQCON+L)+(NINCON+L)+(NINCON+2)*(L+6)
C
C                         The length of the array IW(*) must be at least
C
C                           IW1=NINCON+2*L
C
C                         in any case.
C
C  Output..
C      MODE
C                         An output flag that indicates the status
C                         of the constrained curve fit.
C
C                         =-1  a usage error of FC( ) occurred.  The
C                              offending condition is noted with the
C                              SLATEC library error processor, XERMSG.
C                              In case the working arrays W(*) or IW(*)
C                              are not long enough, the minimal
C                              acceptable length is printed.
C
C                         = 0  successful constrained curve fit.
C
C                         = 1  the requested equality constraints
C                              are contradictory.
C
C                         = 2  the requested inequality constraints
C                              are contradictory.
C
C                         = 3  both equality and inequality constraints
C                              are contradictory.
C
C      COEFF(*)
C                         If the output value of MODE=0 or 1, this array
C                         contains the unknowns obtained from the least
C                         squares fitting process.  These N=NBKPT-NORD
C                         parameters are the B-spline coefficients.
C                         For MODE=1, the equality constraints are
C                         contradictory.  To make the fitting process
C                         more robust, the equality constraints are
C                         satisfied in a least squares sense.  In this
C                         case the array COEFF(*) contains B-spline
C                         coefficients for this extended concept of a
C                         solution.  If MODE=-1,2 or 3 on output, the
C                         array COEFF(*) is undefined.
C
C  Working Arrays..
C      W(*),IW(*)
C                         These arrays are respectively typed REAL and
C                         INTEGER.
C                         Their required lengths are specified as input
C                         parameters in IW(1), IW(2) noted above.  The
C                         contents of W(*) must not be modified by the
C                         user if the variance function is desired.
C
C  Evaluating the
C  Variance Function..
C                         To evaluate the variance function (assuming
C                         that the uncertainties of the Y values were
C                         provided to  FC( ) and an input value of
C                         MODE=2 or 4 was used), use the function
C                         subprogram  CV( )
C
C                           VAR=CV(XVAL,NDATA,NCONST,NORD,NBKPT,
C                                  BKPT,W)
C
C                         Here XVAL is the point where the variance is
C                         desired.  The other arguments have the same
C                         meaning as in the usage of FC( ).
C
C                         For those users employing the old problem
C                         designation, let MDATA be the number of data
C                         points in the problem.  (This may be different
C                         from NDATA if the old problem designation
C                         feature was used.)  The value, VAR, should be
C                         multiplied by the quantity
C
C                         REAL(MAX(NDATA-N,1))/MAX(MDATA-N,1)
C
C                         The output of this subprogram is not defined
C                         if an input value of MODE=1 or 3 was used in
C                         FC( ) or if an output value of MODE=-1, 2, or
C                         3 was obtained.  The variance function, except
C                         for the scaling factor noted above, is given
C                         by
C
C                           VAR=(transpose of B(XVAL))*C*B(XVAL)
C
C                         The vector B(XVAL) is the B-spline basis
C                         function values at X=XVAL.
C                         The covariance matrix, C, of the solution
C                         coefficients accounts only for the least
C                         squares equations and the explicitly stated
C                         equality constraints.  This fact must be
C                         considered when interpreting the variance
C                         function from a data fitting problem that has
C                         inequality constraints on the fitted curve.
C
C  Evaluating the
C  Fitted Curve..
C                         To evaluate derivative number IDER at XVAL,
C                         use the function subprogram BVALU( ).
C
C                           F = BVALU(BKPT,COEFF,NBKPT-NORD,NORD,IDER,
C                                      XVAL,INBV,WORKB)
C
C                         The output of this subprogram will not be
C                         defined unless an output value of MODE=0 or 1
C                         was obtained from  FC( ), XVAL is in the data
C                         interval, and IDER is nonnegative and .LT.
C                         NORD.
C
C                         The first time BVALU( ) is called, INBV=1
C                         must be specified.  This value of INBV is the
C                         overwritten by BVALU( ).  The array WORKB(*)
C                         must be of length at least 3*NORD, and must
C                         not be the same as the W(*) array used in
C                         the call to FC( ).
C
C                         BVALU( ) expects the breakpoint array BKPT(*)
C                         to be sorted.
C
C***REFERENCES  R. J. Hanson, Constrained least squares curve fitting
C                 to discrete data using B-splines, a users guide,
C                 Report SAND78-1291, Sandia Laboratories, December
C                 1978.
C***ROUTINES CALLED  FCMN
C***REVISION HISTORY  (YYMMDD)
C   780801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900510  Convert references to XERRWV to references to XERMSG.  (RWC)
C   900607  Editorial changes to Prologue to make Prologues for EFC,
C           DEFC, FC, and DFC look as much the same as possible.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  FC
      REAL             BKPT(*), COEFF(*), SDDATA(*), W(*), XCONST(*),
     *   XDATA(*), YCONST(*), YDATA(*)
      INTEGER IW(*), MODE, NBKPT, NCONST, NDATA, NDERIV(*), NORD
C
      EXTERNAL FCMN
C
      INTEGER I1, I2, I3, I4, I5, I6, I7, MDG, MDW
C
C***FIRST EXECUTABLE STATEMENT  FC
      MDG = NBKPT - NORD + 3
      MDW = NBKPT - NORD + 1 + NCONST
C                         USAGE IN FCMN( ) OF W(*)..
C     I1,...,I2-1      G(*,*)
C
C     I2,...,I3-1      XTEMP(*)
C
C     I3,...,I4-1      PTEMP(*)
C
C     I4,...,I5-1      BKPT(*) (LOCAL TO FCMN( ))
C
C     I5,...,I6-1      BF(*,*)
C
C     I6,...,I7-1      W(*,*)
C
C     I7,...           WORK(*) FOR LSEI( )
C
      I1 = 1
      I2 = I1 + MDG*(NORD+1)
      I3 = I2 + MAX(NDATA,NBKPT)
      I4 = I3 + MAX(NDATA,NBKPT)
      I5 = I4 + NBKPT
      I6 = I5 + NORD*NORD
      I7 = I6 + MDW*(NBKPT-NORD+1)
      CALL  FCMN(NDATA, XDATA, YDATA, SDDATA, NORD, NBKPT, BKPT, NCONST,
     1   XCONST, YCONST, NDERIV, MODE, COEFF, W(I5), W(I2), W(I3),
     2   W(I4), W(I1), MDG, W(I6), MDW, W(I7), IW)
      RETURN
      END
