*DECK PCHDOC
      SUBROUTINE PCHDOC
C***BEGIN PROLOGUE  PCHDOC
C***PURPOSE  Documentation for PCHIP, a Fortran package for piecewise
C            cubic Hermite interpolation of data.
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E1A, Z
C***TYPE      ALL (PCHDOC-A)
C***KEYWORDS  CUBIC HERMITE INTERPOLATION, DOCUMENTATION,
C             MONOTONE INTERPOLATION, PCHIP,
C             PIECEWISE CUBIC INTERPOLATION
C***AUTHOR  Fritsch, F. N., (LLNL)
C             Lawrence Livermore National Laboratory
C             P.O. Box 808  (L-316)
C             Livermore, CA  94550
C             FTS 532-4275, (510) 422-4275
C***DESCRIPTION
C
C            PCHIP:  Piecewise Cubic Hermite Interpolation Package
C
C      This document describes the contents of PCHIP, which is a
C   Fortran package for piecewise cubic Hermite interpolation of data.
C   It features software to produce a monotone and "visually pleasing"
C   interpolant to monotone data.  As is demonstrated in Reference 4,
C   such an interpolant may be more reasonable than a cubic spline if
C   the data contains both "steep" and "flat" sections.  Interpola-
C   tion of cumulative probability distribution functions is another
C   application.  (See References 2-4 for examples.)
C
C
C      All piecewise cubic functions in PCHIP are represented in
C   cubic Hermite form; that is, f(x) is determined by its values
C   F(I) and derivatives D(I) at the breakpoints X(I), I=1(1)N.
C   Throughout the package a PCH function is represented by the
C   five variables  N, X, F, D, INCFD:
C     N     - number of data points;
C     X     - abscissa values for the data points;
C     F     - ordinates (function values) for the data points;
C     D     - slopes (derivative values) at the data points;
C     INCFD - increment between successive elements in the F- and
C             D-arrays (more on this later).
C   These appear together and in the same order in all calls.
C
C      The double precision equivalents of the PCHIP routines are
C   obtained from the single precision names by prefixing the
C   single precision names with a D.  For example, the double
C   precision equivalent of PCHIM is DPCHIM.
C
C      The contents of the package are as follows:
C
C   1. Determine Derivative Values.
C
C      NOTE:  These routines provide alternate ways of determining D
C             if these values are not already known.
C
C         PCHIM -- Piecewise Cubic Hermite Interpolation to Monotone
C               data.
C               Used if the data are monotonic or if the user wants
C               to guarantee that the interpolant stays within the
C               limits of the data.  (See Reference 3.)
C
C         PCHIC -- Piecewise Cubic Hermite Interpolation Coefficients.
C               Used if neither of the above conditions holds, or if
C               the user wishes control over boundary derivatives.
C               Will generally reproduce monotonicity on subintervals
C               over which the data are monotonic.
C
C         PCHSP -- Piecewise Cubic Hermite Spline.
C               Produces a cubic spline interpolator in cubic Hermite
C               form.  Provided primarily for easy comparison of the
C               spline with other piecewise cubic interpolants.  (A
C               modified version of de Boor's CUBSPL, Reference 1.)
C
C   2. Evaluate, Differentiate, or Integrate Resulting PCH Function.
C
C      NOTE:  If derivative values are available from some other
C             source, these routines can be used without calling
C             any of the previous routines.
C
C         CHFEV -- Cubic Hermite Function EValuator.
C               Evaluates a single cubic Hermite function at an array
C               of points.  Used when the interval is known, as in
C               graphing applications.  Called by PCHFE.
C
C         PCHFE -- Piecewise Cubic Hermite Function Evaluator.
C               Used when the interval is unknown or the evaluation
C               array spans more than one data interval.
C
C         CHFDV -- Cubic Hermite Function and Derivative Evaluator.
C               Evaluates a single cubic Hermite function and its
C               first derivative at an array of points.  Used when
C               the interval is known, as in graphing applications.
C               Called by PCHFD.
C
C         PCHFD -- Piecewise Cubic Hermite Function and Derivative
C               Evaluator.
C               Used when the interval is unknown or the evaluation
C               array spans more than one data interval.
C
C         PCHID -- Piecewise Cubic Hermite Integrator, Data Limits.
C               Computes the definite integral of a piecewise cubic
C               Hermite function when the integration limits are data
C               points.
C
C         PCHIA -- Piecewise Cubic Hermite Integrator, Arbitrary Limits.
C               Computes the definite integral of a piecewise cubic
C               Hermite function over an arbitrary finite interval.
C
C   3. Utility routines.
C
C         PCHBS -- Piecewise Cubic Hermite to B-Spline converter.
C               Converts a PCH function to B-representation, so that
C               it can be used with other elements of the B-spline
C               package (see BSPDOC).
C
C         PCHCM -- Piecewise Cubic Hermite, Check Monotonicity of.
C               Checks the monotonicity of an arbitrary PCH function.
C               Might be used with PCHSP to build a polyalgorithm for
C               piecewise C-2 interpolation.
C
C   4. Internal routines.
C
C         CHFIE -- Cubic Hermite Function Integral Evaluator.
C               (Real function called by PCHIA.)
C
C         CHFCM -- Cubic Hermite Function, Check Monotonicity of.
C               (Integer function called by PCHCM.)
C
C         PCHCE -- PCHIC End Derivative Setter.
C               (Called by PCHIC.)
C
C         PCHCI -- PCHIC Initial Derivative Setter.
C               (Called by PCHIC.)
C
C         PCHCS -- PCHIC Monotonicity Switch Derivative Setter.
C               (Called by PCHIC.)
C
C         PCHDF -- PCHIP Finite Difference Formula.
C               (Real function called by PCHCE and PCHSP.)
C
C         PCHST -- PCHIP Sign Testing Routine.
C               (Real function called by various PCHIP routines.)
C
C         PCHSW -- PCHCS Switch Excursion Adjuster.
C               (Called by PCHCS.)
C
C   The calling sequences for these routines are described in the
C   prologues of the respective routines.
C
C
C      INCFD, the increment between successive elements in the F-
C   and D-arrays is included in the representation of a PCH function
C   in this package to facilitate two-dimensional applications.  For
C   "normal" usage INCFD=1, and F and D are one-dimensional arrays.
C   one would call PCHxx (where "xx" is "IM", "IC", or "SP") with
C
C              N, X, F, D, 1  .
C
C   Suppose, however, that one has data on a rectangular mesh,
C
C         F2D(I,J) = value at (X(I), Y(J)),  I=1(1)NX,
C                                            J=1(1)NY.
C   Assume the following dimensions:
C
C         REAL  X(NXMAX), Y(NYMAX)
C         REAL  F2D(NXMAX,NYMAX), FX(NXMAX,NYMAX), FY(NXMAX,NYMAX)
C
C   where  2.LE.NX.LE.NXMAX AND 2.LE.NY.LE.NYMAX .  To interpolate
C   in X along the line  Y = Y(J), call PCHxx with
C
C              NX, X, F2D(1,J), FX(1,J), 1  .
C
C   To interpolate along the line X = X(I), call PCHxx with
C
C              NY, Y, F2D(I,1), FY(I,1), NXMAX  .
C
C   (This example assumes the usual columnwise storage of 2-D arrays
C    in Fortran.)
C
C***REFERENCES  1. Carl de Boor, A Practical Guide to Splines, Springer-
C                 Verlag, New York, 1978 (esp. Chapter IV, pp.49-62).
C               2. F. N. Fritsch, Piecewise Cubic Hermite Interpolation
C                 Package, Report UCRL-87285, Lawrence Livermore Natio-
C                 nal Laboratory, July 1982.  [Poster presented at the
C                 SIAM 30th Anniversary Meeting, 19-23 July 1982.]
C               3. F. N. Fritsch and J. Butland, A method for construc-
C                 ting local monotone piecewise cubic interpolants, SIAM
C                 Journal on Scientific and Statistical Computing 5, 2
C                 (June 1984), pp. 300-304.
C               4. F. N. Fritsch and R. E. Carlson, Monotone piecewise
C                 cubic interpolation, SIAM Journal on Numerical Ana-
C                 lysis 17, 2 (April 1980), pp. 238-246.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   811106  DATE WRITTEN
C   870930  Updated Reference 3.
C   890414  Changed PCHMC and CHFMC to PCHCM and CHFCM, respectively,
C           and augmented description of PCHCM.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   910826  1. Revised purpose, clarified role of argument INCFD,
C              corrected error in example, and removed redundant
C              reference list.
C           2. Added description of PCHBS.  (FNF)
C   920429  Revised format and order of references.  (WRB,FNF)
C   930505  Changed CHFIV to CHFIE.  (FNF)
C***END PROLOGUE  PCHDOC
C-----------------------------------------------------------------------
C     THIS IS A DUMMY SUBROUTINE, AND SHOULD NEVER BE CALLED.
C
C***FIRST EXECUTABLE STATEMENT  PCHDOC
      RETURN
C------------- LAST LINE OF PCHDOC FOLLOWS -----------------------------
      END
