*DECK BSPDOC
      SUBROUTINE BSPDOC
C***BEGIN PROLOGUE  BSPDOC
C***PURPOSE  Documentation for BSPLINE, a package of subprograms for
C            working with piecewise polynomial functions
C            in B-representation.
C***LIBRARY   SLATEC
C***CATEGORY  E, E1A, K, Z
C***TYPE      ALL (BSPDOC-A)
C***KEYWORDS  B-SPLINE, DOCUMENTATION, SPLINES
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C         BSPDOC is a non-executable, B-spline documentary routine.
C         The narrative describes a B-spline and the routines
C         necessary to manipulate B-splines at a fairly high level.
C         The basic package described herein is that of reference
C         5 with names altered to prevent duplication and conflicts
C         with routines from reference 3.  The call lists used here
C         are also different.  Work vectors were added to ensure
C         portability and proper execution in an overlay environ-
C         ment.  These work arrays can be used for other purposes
C         except as noted in BSPVN.  While most of the original
C         routines in reference 5 were restricted to orders 20
C         or less, this restriction was removed from all routines
C         except the quadrature routine BSQAD.  (See the section
C         below on differentiation and integration for details.)
C
C         The subroutines referenced below are single precision
C         routines.  Corresponding double precision versions are also
C         part of the package, and these are referenced by prefixing
C         a D in front of the single precision name.  For example,
C         BVALU and DBVALU are the single and double precision
C         versions for evaluating a B-spline or any of its deriva-
C         tives in the B-representation.
C
C                ****Description of B-Splines****
C
C     A collection of polynomials of fixed degree K-1 defined on a
C     subdivision (X(I),X(I+1)), I=1,...,M-1 of (A,B) with X(1)=A,
C     X(M)=B is called a B-spline of order K.  If the spline has K-2
C     continuous derivatives on (A,B), then the B-spline is simply
C     called a spline of order K.  Each of the M-1 polynomial pieces
C     has K coefficients, making a total of K(M-1) parameters.  This
C     B-spline and its derivatives have M-2 jumps at the subdivision
C     points X(I), I=2,...,M-1.  Continuity requirements at these
C     subdivision points add constraints and reduce the number of free
C     parameters.  If a B-spline is continuous at each of the M-2 sub-
C     division points, there are K(M-1)-(M-2) free parameters; if in
C     addition the B-spline has continuous first derivatives, there
C     are K(M-1)-2(M-2) free parameters, etc., until we get to a
C     spline where we have K(M-1)-(K-1)(M-2) = M+K-2 free parameters.
C     Thus, the principle is that increasing the continuity of
C     derivatives decreases the number of free parameters and
C     conversely.
C
C     The points at which the polynomials are tied together by the
C     continuity conditions are called knots.  If two knots are
C     allowed to come together at some X(I), then we say that we
C     have a knot of multiplicity 2 there, and the knot values are
C     the X(I) value.  If we reverse the procedure of the first
C     paragraph, we find that adding a knot to increase multiplicity
C     increases the number of free parameters and, according to the
C     principle above, we thereby introduce a discontinuity in what
C     was the highest continuous derivative at that knot.  Thus, the
C     number of free parameters is N = NU+K-2 where NU is the sum
C     of multiplicities at the X(I) values with X(1) and X(M) of
C     multiplicity 1 (NU = M if all knots are simple, i.e., for a
C     spline, all knots have multiplicity 1.)  Each knot can have a
C     multiplicity of at most K.  A B-spline is commonly written in the
C     B-representation
C
C               Y(X) = sum( A(I)*B(I,X), I=1 , N)
C
C     to show the explicit dependence of the spline on the free
C     parameters or coefficients A(I)=BCOEF(I) and basis functions
C     B(I,X).  These basis functions are themselves special B-splines
C     which are zero except on (at most) K adjoining intervals where
C     each B(I,X) is positive and, in most cases, hat or bell-
C     shaped.  In order for the nonzero part of B(1,X) to be a spline
C     covering (X(1),X(2)), it is necessary to put K-1 knots to the
C     left of A and similarly for B(N,X) to the right of B.  Thus, the
C     total number of knots for this representation is NU+2K-2 = N+K.
C     These knots are carried in an array T(*) dimensioned by at least
C     N+K.  From the construction, A=T(K) and B=T(N+1) and the spline is
C     defined on T(K).LE.X.LE.T(N+1).  The nonzero part of each basis
C     function lies in the  Interval (T(I),T(I+K)).  In many problems
C     where extrapolation beyond A or B is not anticipated, it is common
C     practice to set T(1)=T(2)=...=T(K)=A and T(N+1)=T(N+2)=...=
C     T(N+K)=B.  In summary, since T(K) and T(N+1) as well as
C     interior knots can have multiplicity K, the number of free
C     parameters N = sum of multiplicities - K.  The fact that each
C     B(I,X) function is nonzero over at most K intervals means that
C     for a given X value, there are at most K nonzero terms of the
C     sum.  This leads to banded matrices in linear algebra problems,
C     and references 3 and 6 take advantage of this in con-
C     structing higher level routines to achieve speed and avoid
C     ill-conditioning.
C
C                     ****Basic Routines****
C
C     The basic routines which most casual users will need are those
C     concerned with direct evaluation of splines or B-splines.
C     Since the B-representation, denoted by (T,BCOEF,N,K), is
C     preferred because of numerical stability, the knots T(*), the
C     B-spline coefficients BCOEF(*), the number of coefficients N,
C     and the order K of the polynomial pieces (of degree K-1) are
C     usually given.  While the knot array runs from T(1) to T(N+K),
C     the B-spline is normally defined on the interval T(K).LE.X.LE.
C     T(N+1).  To evaluate the B-spline or any of its derivatives
C     on this interval, one can use
C
C                  Y = BVALU(T,BCOEF,N,K,ID,X,INBV,WORK)
C
C     where ID is an integer for the ID-th derivative, 0.LE.ID.LE.K-1.
C     ID=0 gives the zero-th derivative or B-spline value at X.
C     If X.LT.T(K) or X.GT.T(N+1), whether by mistake or the result
C     of round off accumulation in incrementing X, BVALU gives a
C     diagnostic.  INBV is an initialization parameter which is set
C     to 1 on the first call.  Distinct splines require distinct
C     INBV parameters.  WORK is a scratch vector of length at least
C     3*K.
C
C     When more conventional communication is needed for publication,
C     physical interpretation, etc., the B-spline coefficients can
C     be converted to piecewise polynomial (PP) coefficients.  Thus,
C     the breakpoints (distinct knots) XI(*), the number of
C     polynomial pieces LXI, and the (right) derivatives C(*,J) at
C     each breakpoint XI(J) are needed to define the Taylor
C     expansion to the right of XI(J) on each interval XI(J).LE.
C     X.LT.XI(J+1), J=1,LXI where XI(1)=A and XI(LXI+1)=B.
C     These are obtained from the (T,BCOEF,N,K) representation by
C
C                CALL BSPPP(T,BCOEF,N,K,LDC,C,XI,LXI,WORK)
C
C     where LDC.GE.K is the leading dimension of the matrix C and
C     WORK is a scratch vector of length at least K*(N+3).
C     Then the PP-representation (C,XI,LXI,K) of Y(X), denoted
C     by Y(J,X) on each interval XI(J).LE.X.LT.XI(J+1), is
C
C     Y(J,X) = sum( C(I,J)*((X-XI(J))**(I-1))/factorial(I-1), I=1,K)
C
C     for J=1,...,LXI.  One must view this conversion from the B-
C     to the PP-representation with some skepticism because the
C     conversion may lose significant digits when the B-spline
C     varies in an almost discontinuous fashion.  To evaluate
C     the B-spline or any of its derivatives using the PP-
C     representation, one uses
C
C                Y = PPVAL(LDC,C,XI,LXI,K,ID,X,INPPV)
C
C     where ID and INPPV have the same meaning and usage as ID and
C     INBV in BVALU.
C
C     To determine to what extent the conversion process loses
C     digits, compute the relative error ABS((Y1-Y2)/Y2) over
C     the X interval with Y1 from PPVAL and Y2 from BVALU.  A
C     major reason for considering PPVAL is that evaluation is
C     much faster than that from BVALU.
C
C     Recall that when multiple knots are encountered, jump type
C     discontinuities in the B-spline or its derivatives occur
C     at these knots, and we need to know that BVALU and PPVAL
C     return right limiting values at these knots except at
C     X=B where left limiting values are returned.  These values
C     are used for the Taylor expansions about left end points of
C     breakpoint intervals.  That is, the derivatives C(*,J) are
C     right derivatives.  Note also that a computed X value which,
C     mathematically, would be a knot value may differ from the knot
C     by a round off error.  When this happens in evaluating a dis-
C     continuous B-spline or some discontinuous derivative, the
C     value at the knot and the value at X can be radically
C     different.  In this case, setting X to a T or XI value makes
C     the computation precise.  For left limiting values at knots
C     other than X=B, see the prologues to BVALU and other
C     routines.
C
C                     ****Interpolation****
C
C     BINTK is used to generate B-spline parameters (T,BCOEF,N,K)
C     which will interpolate the data by calls to BVALU.  A similar
C     interpolation can also be done for cubic splines using BINT4
C     or the code in reference 7.  If the PP-representation is given,
C     one can evaluate this representation at an appropriate number of
C     abscissas to create data then use BINTK or BINT4 to generate
C     the B-representation.
C
C               ****Differentiation and Integration****
C
C     Derivatives of B-splines are obtained from BVALU or PPVAL.
C     Integrals are obtained from BSQAD using the B-representation
C     (T,BCOEF,N,K) and PPQAD using the PP-representation (C,XI,LXI,
C     K).  More complicated integrals involving the product of a
C     of a function F and some derivative of a B-spline can be
C     evaluated with BFQAD or PFQAD using the B- or PP- represen-
C     tations respectively.  All quadrature routines, except for PPQAD,
C     are limited in accuracy to 18 digits or working precision,
C     whichever is smaller.  PPQAD is limited to working precision
C     only.  In addition, the order K for BSQAD is limited to 20 or
C     less.  If orders greater than 20 are required, use BFQAD with
C     F(X) = 1.
C
C                      ****Extrapolation****
C
C     Extrapolation outside the interval (A,B) can be accomplished
C     easily by the PP-representation using PPVAL.  However,
C     caution should be exercised, especially when several knots
C     are located at A or B or when the extrapolation is carried
C     significantly beyond A or B.  On the other hand, direct
C     evaluation with BVALU outside A=T(K).LE.X.LE.T(N+1)=B
C     produces an error message, and some manipulation of the knots
C     and coefficients are needed to extrapolate with BVALU.  This
C     process is described in reference 6.
C
C                ****Curve Fitting and Smoothing****
C
C     Unless one has many accurate data points, direct inter-
C     polation is not recommended for summarizing data.  The
C     results are often not in accordance with intuition since the
C     fitted curve tends to oscillate through the set of points.
C     Monotone splines (reference 7) can help curb this undulating
C     tendency but constrained least squares is more likely to give an
C     acceptable fit with fewer parameters.  Subroutine FC, des-
C     cribed in reference 6, is recommended for this purpose.  The
C     output from this fitting process is the B-representation.
C
C              **** Routines in the B-Spline Package ****
C
C                      Single Precision Routines
C
C         The subroutines referenced below are SINGLE PRECISION
C         routines. Corresponding DOUBLE PRECISION versions are also
C         part of the package and these are referenced by prefixing
C         a D in front of the single precision name. For example,
C         BVALU and DBVALU are the SINGLE and DOUBLE PRECISION
C         versions for evaluating a B-spline or any of its deriva-
C         tives in the B-representation.
C
C     BINT4 - interpolates with splines of order 4
C     BINTK - interpolates with splines of order k
C     BSQAD - integrates the B-representation on subintervals
C     PPQAD - integrates the PP-representation
C     BFQAD - integrates the product of a function F and any spline
C             derivative in the B-representation
C     PFQAD - integrates the product of a function F and any spline
C             derivative in the PP-representation
C     BVALU - evaluates the B-representation or a derivative
C     PPVAL - evaluates the PP-representation or a derivative
C     INTRV - gets the largest index of the knot to the left of x
C     BSPPP - converts from B- to PP-representation
C     BSPVD - computes nonzero basis functions and derivatives at x
C     BSPDR - sets up difference array for BSPEV
C     BSPEV - evaluates the B-representation and derivatives
C     BSPVN - called by BSPEV, BSPVD, BSPPP and BINTK for function and
C             derivative evaluations
C                        Auxiliary Routines
C
C       BSGQ8,PPGQ8,BNSLV,BNFAC,XERMSG,DBSGQ8,DPPGQ8,DBNSLV,DBNFAC
C
C                    Machine Dependent Routines
C
C                      I1MACH, R1MACH, D1MACH
C
C***REFERENCES  1. D. E. Amos, Computation with splines and
C                 B-splines, Report SAND78-1968, Sandia
C                 Laboratories, March 1979.
C               2. D. E. Amos, Quadrature subroutines for splines and
C                 B-splines, Report SAND79-1825, Sandia Laboratories,
C                 December 1979.
C               3. Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C               4. Carl de Boor, On calculating with B-Splines, Journal
C                 of Approximation Theory 6, (1972), pp. 50-62.
C               5. Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C               6. R. J. Hanson, Constrained least squares curve fitting
C                 to discrete data using B-splines, a users guide,
C                 Report SAND78-1291, Sandia Laboratories, December
C                 1978.
C               7. F. N. Fritsch and R. E. Carlson, Monotone piecewise
C                 cubic interpolation, SIAM Journal on Numerical Ana-
C                 lysis 17, 2 (April 1980), pp. 238-246.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   810223  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900723  PURPOSE section revised.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BSPDOC
C***FIRST EXECUTABLE STATEMENT  BSPDOC
      RETURN
      END
