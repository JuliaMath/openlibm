*DECK DEFC
      SUBROUTINE DEFC (NDATA, XDATA, YDATA, SDDATA, NORD, NBKPT, BKPT,
     +   MDEIN, MDEOUT, COEFF, LW, W)
C***BEGIN PROLOGUE  DEFC
C***PURPOSE  Fit a piecewise polynomial curve to discrete data.
C            The piecewise polynomials are represented as B-splines.
C            The fitting is done in a weighted least squares sense.
C***LIBRARY   SLATEC
C***CATEGORY  K1A1A1, K1A2A, L8A3
C***TYPE      DOUBLE PRECISION (EFC-S, DEFC-D)
C***KEYWORDS  B-SPLINE, CONSTRAINED LEAST SQUARES, CURVE FITTING
C***AUTHOR  Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C      This subprogram fits a piecewise polynomial curve
C      to discrete data.  The piecewise polynomials are
C      represented as B-splines.
C      The fitting is done in a weighted least squares sense.
C
C      The data can be processed in groups of modest size.
C      The size of the group is chosen by the user.  This feature
C      may be necessary for purposes of using constrained curve fitting
C      with subprogram DFC( ) on a very large data set.
C
C      For a description of the B-splines and usage instructions to
C      evaluate them, see
C
C      C. W. de Boor, Package for Calculating with B-Splines.
C                     SIAM J. Numer. Anal., p. 441, (June, 1977).
C
C      For further discussion of (constrained) curve fitting using
C      B-splines, see
C
C      R. J. Hanson, Constrained Least Squares Curve Fitting
C                   to Discrete Data Using B-Splines, a User's
C                   Guide. Sandia Labs. Tech. Rept. SAND-78-1291,
C                   December, (1978).
C
C  Input.. All TYPE REAL variables are DOUBLE PRECISION
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
C                         Internal to DEFC( ) the extreme end knots may
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
C        MDEIN
C                         An integer flag, with one of two possible
C                         values (1 or 2), that directs the subprogram
C                         action with regard to new data points provided
C                         by the user.
C
C                         =1  The first time that DEFC( ) has been
C                         entered.  There are NDATA points to process.
C
C                         =2  This is another entry to DEFC().  The sub-
C                         program DEFC( ) has been entered with MDEIN=1
C                         exactly once before for this problem.  There
C                         are NDATA new additional points to merge and
C                         process with any previous points.
C                         (When using DEFC( ) with MDEIN=2 it is import-
C                         ant that the set of knots remain fixed at the
C                         same values for all entries to DEFC( ).)
C       LW
C                         The amount of working storage actually
C                         allocated for the working array W(*).
C                         This quantity is compared with the
C                         actual amount of storage needed in DEFC( ).
C                         Insufficient storage allocated for W(*) is
C                         an error.  This feature was included in DEFC
C                         because misreading the storage formula
C                         for W(*) might very well lead to subtle
C                         and hard-to-find programming bugs.
C
C                         The length of the array W(*) must satisfy
C
C                         LW .GE. (NBKPT-NORD+3)*(NORD+1)+
C                                 (NBKPT+1)*(NORD+1)+
C                               2*MAX(NDATA,NBKPT)+NBKPT+NORD**2
C
C  Output.. All TYPE REAL variables are DOUBLE PRECISION
C      MDEOUT
C                         An output flag that indicates the status
C                         of the curve fit.
C
C                         =-1  A usage error of DEFC( ) occurred.  The
C                         offending condition is noted with the SLATEC
C                         library error processor, XERMSG( ).  In case
C                         the working array W(*) is not long enough, the
C                         minimal acceptable length is printed.
C
C                         =1  The B-spline coefficients for the fitted
C                         curve have been returned in array COEFF(*).
C
C                         =2  Not enough data has been processed to
C                         determine the B-spline coefficients.
C                         The user has one of two options.  Continue
C                         to process more data until a unique set
C                         of coefficients is obtained, or use the
C                         subprogram DFC( ) to obtain a specific
C                         set of coefficients.  The user should read
C                         the usage instructions for DFC( ) for further
C                         details if this second option is chosen.
C      COEFF(*)
C                         If the output value of MDEOUT=1, this array
C                         contains the unknowns obtained from the least
C                         squares fitting process.  These N=NBKPT-NORD
C                         parameters are the B-spline coefficients.
C                         For MDEOUT=2, not enough data was processed to
C                         uniquely determine the B-spline coefficients.
C                         In this case, and also when MDEOUT=-1, all
C                         values of COEFF(*) are set to zero.
C
C                         If the user is not satisfied with the fitted
C                         curve returned by DEFC( ), the constrained
C                         least squares curve fitting subprogram DFC( )
C                         may be required.  The work done within DEFC( )
C                         to accumulate the data can be utilized by
C                         the user, if so desired.  This involves
C                         saving the first (NBKPT-NORD+3)*(NORD+1)
C                         entries of W(*) and providing this data
C                         to DFC( ) with the "old problem" designation.
C                         The user should read the usage instructions
C                         for subprogram DFC( ) for further details.
C
C  Working Array.. All TYPE REAL variables are DOUBLE PRECISION
C      W(*)
C                         This array is typed DOUBLE PRECISION.
C                         Its length is  specified as an input parameter
C                         in LW as noted above.  The contents of W(*)
C                         must not be modified by the user between calls
C                         to DEFC( ) with values of MDEIN=1,2,2,... .
C                         The first (NBKPT-NORD+3)*(NORD+1) entries of
C                         W(*) are acceptable as direct input to DFC( )
C                         for an "old problem" only when MDEOUT=1 or 2.
C
C  Evaluating the
C  Fitted Curve..
C                         To evaluate derivative number IDER at XVAL,
C                         use the function subprogram DBVALU( ).
C
C                         F = DBVALU(BKPT,COEFF,NBKPT-NORD,NORD,IDER,
C                                      XVAL,INBV,WORKB)
C
C                         The output of this subprogram will not be
C                         defined unless an output value of MDEOUT=1
C                         was obtained from DEFC( ), XVAL is in the data
C                         interval, and IDER is nonnegative and .LT.
C                         NORD.
C
C                         The first time DBVALU( ) is called, INBV=1
C                         must be specified.  This value of INBV is the
C                         overwritten by DBVALU( ).  The array WORKB(*)
C                         must be of length at least 3*NORD, and must
C                         not be the same as the W(*) array used in the
C                         call to DEFC( ).
C
C                         DBVALU( ) expects the breakpoint array BKPT(*)
C                         to be sorted.
C
C***REFERENCES  R. J. Hanson, Constrained least squares curve fitting
C                 to discrete data using B-splines, a users guide,
C                 Report SAND78-1291, Sandia Laboratories, December
C                 1978.
C***ROUTINES CALLED  DEFCMN
C***REVISION HISTORY  (YYMMDD)
C   800801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900510  Change Prologue comments to refer to XERMSG.  (RWC)
C   900607  Editorial changes to Prologue to make Prologues for EFC,
C           DEFC, FC, and DFC look as much the same as possible.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DEFC
C
C      SUBROUTINE           FUNCTION/REMARKS
C
C      DBSPVN( )          COMPUTE FUNCTION VALUES OF B-SPLINES.  FROM
C                         THE B-SPLINE PACKAGE OF DE BOOR NOTED ABOVE.
C
C      DBNDAC( ),         BANDED LEAST SQUARES MATRIX PROCESSORS.
C      DBNDSL( )          FROM LAWSON-HANSON, SOLVING LEAST
C                         SQUARES PROBLEMS.
C
C      DSORT( )           DATA SORTING SUBROUTINE, FROM THE
C                         SANDIA MATH. LIBRARY, SAND77-1441.
C
C      XERMSG( )          ERROR HANDLING ROUTINE
C                         FOR THE SLATEC MATH. LIBRARY.
C                         SEE SAND78-1189, BY R. E. JONES.
C
C      DCOPY( ),DSCAL( )  SUBROUTINES FROM THE BLAS PACKAGE.
C
C                         WRITTEN BY R. HANSON, SANDIA NATL. LABS.,
C                         ALB., N. M., AUGUST-SEPTEMBER, 1980.
C
      DOUBLE PRECISION BKPT(*),COEFF(*),W(*),SDDATA(*),XDATA(*),YDATA(*)
      INTEGER LW, MDEIN, MDEOUT, NBKPT, NDATA, NORD
C
      EXTERNAL DEFCMN
C
      INTEGER LBF, LBKPT, LG, LPTEMP, LWW, LXTEMP, MDG, MDW
C
C***FIRST EXECUTABLE STATEMENT  DEFC
C     LWW=1               USAGE IN DEFCMN( ) OF W(*)..
C     LWW,...,LG-1        W(*,*)
C
C     LG,...,LXTEMP-1     G(*,*)
C
C     LXTEMP,...,LPTEMP-1 XTEMP(*)
C
C     LPTEMP,...,LBKPT-1  PTEMP(*)
C
C     LBKPT,...,LBF       BKPT(*) (LOCAL TO DEFCMN( ))
C
C     LBF,...,LBF+NORD**2 BF(*,*)
C
      MDG = NBKPT+1
      MDW = NBKPT-NORD+3
      LWW = 1
      LG = LWW + MDW*(NORD+1)
      LXTEMP = LG + MDG*(NORD+1)
      LPTEMP = LXTEMP + MAX(NDATA,NBKPT)
      LBKPT  = LPTEMP + MAX(NDATA,NBKPT)
      LBF = LBKPT + NBKPT
      CALL DEFCMN(NDATA,XDATA,YDATA,SDDATA,
     1         NORD,NBKPT,BKPT,
     2         MDEIN,MDEOUT,
     3         COEFF,
     4         W(LBF),W(LXTEMP),W(LPTEMP),W(LBKPT),
     5         W(LG),MDG,W(LWW),MDW,LW)
      RETURN
      END
