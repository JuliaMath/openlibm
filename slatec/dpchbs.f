*DECK DPCHBS
      SUBROUTINE DPCHBS (N, X, F, D, INCFD, KNOTYP, NKNOTS, T, BCOEF,
     +   NDIM, KORD, IERR)
C***BEGIN PROLOGUE  DPCHBS
C***PURPOSE  Piecewise Cubic Hermite to B-Spline converter.
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E3
C***TYPE      DOUBLE PRECISION (PCHBS-S, DPCHBS-D)
C***KEYWORDS  B-SPLINES, CONVERSION, CUBIC HERMITE INTERPOLATION,
C             PIECEWISE CUBIC INTERPOLATION
C***AUTHOR  Fritsch, F. N., (LLNL)
C             Computing and Mathematics Research Division
C             Lawrence Livermore National Laboratory
C             P.O. Box 808  (L-316)
C             Livermore, CA  94550
C             FTS 532-4275, (510) 422-4275
C***DESCRIPTION
C
C *Usage:
C
C        INTEGER  N, INCFD, KNOTYP, NKNOTS, NDIM, KORD, IERR
C        PARAMETER  (INCFD = ...)
C        DOUBLE PRECISION  X(nmax), F(INCFD,nmax), D(INCFD,nmax),
C       *      T(2*nmax+4), BCOEF(2*nmax)
C
C        CALL DPCHBS (N, X, F, D, INCFD, KNOTYP, NKNOTS, T, BCOEF,
C       *             NDIM, KORD, IERR)
C
C *Arguments:
C
C     N:IN  is the number of data points, N.ge.2 .  (not checked)
C
C     X:IN  is the real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.   (not checked)
C           nmax, the dimension of X, must be .ge.N.
C
C     F:IN  is the real array of dependent variable values.
C           F(1+(I-1)*INCFD) is the value corresponding to X(I).
C           nmax, the second dimension of F, must be .ge.N.
C
C     D:IN  is the real array of derivative values at the data points.
C           D(1+(I-1)*INCFD) is the value corresponding to X(I).
C           nmax, the second dimension of D, must be .ge.N.
C
C     INCFD:IN  is the increment between successive values in F and D.
C           This argument is provided primarily for 2-D applications.
C           It may have the value 1 for one-dimensional applications,
C           in which case F and D may be singly-subscripted arrays.
C
C     KNOTYP:IN  is a flag to control the knot sequence.
C           The knot sequence T is normally computed from X by putting
C           a double knot at each X and setting the end knot pairs
C           according to the value of KNOTYP:
C              KNOTYP = 0:  Quadruple knots at X(1) and X(N).  (default)
C              KNOTYP = 1:  Replicate lengths of extreme subintervals:
C                           T( 1 ) = T( 2 ) = X(1) - (X(2)-X(1))  ;
C                           T(M+4) = T(M+3) = X(N) + (X(N)-X(N-1)).
C              KNOTYP = 2:  Periodic placement of boundary knots:
C                           T( 1 ) = T( 2 ) = X(1) - (X(N)-X(N-1));
C                           T(M+4) = T(M+3) = X(N) + (X(2)-X(1))  .
C              Here M=NDIM=2*N.
C           If the input value of KNOTYP is negative, however, it is
C           assumed that NKNOTS and T were set in a previous call.
C           This option is provided for improved efficiency when used
C           in a parametric setting.
C
C     NKNOTS:INOUT  is the number of knots.
C           If KNOTYP.GE.0, then NKNOTS will be set to NDIM+4.
C           If KNOTYP.LT.0, then NKNOTS is an input variable, and an
C              error return will be taken if it is not equal to NDIM+4.
C
C     T:INOUT  is the array of 2*N+4 knots for the B-representation.
C           If KNOTYP.GE.0, T will be returned by DPCHBS with the
C              interior double knots equal to the X-values and the
C              boundary knots set as indicated above.
C           If KNOTYP.LT.0, it is assumed that T was set by a
C              previous call to DPCHBS.  (This routine does **not**
C              verify that T forms a legitimate knot sequence.)
C
C     BCOEF:OUT  is the array of 2*N B-spline coefficients.
C
C     NDIM:OUT  is the dimension of the B-spline space.  (Set to 2*N.)
C
C     KORD:OUT  is the order of the B-spline.  (Set to 4.)
C
C     IERR:OUT  is an error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -4  if KNOTYP.GT.2 .
C              IERR = -5  if KNOTYP.LT.0 and NKNOTS.NE.(2*N+4).
C
C *Description:
C     DPCHBS computes the B-spline representation of the PCH function
C     determined by N,X,F,D.  To be compatible with the rest of PCHIP,
C     DPCHBS includes INCFD, the increment between successive values of
C     the F- and D-arrays.
C
C     The output is the B-representation for the function:  NKNOTS, T,
C     BCOEF, NDIM, KORD.
C
C *Caution:
C     Since it is assumed that the input PCH function has been
C     computed by one of the other routines in the package PCHIP,
C     input arguments N, X, INCFD are **not** checked for validity.
C
C *Restrictions/assumptions:
C     1. N.GE.2 .  (not checked)
C     2. X(i).LT.X(i+1), i=1,...,N .  (not checked)
C     3. INCFD.GT.0 .  (not checked)
C     4. KNOTYP.LE.2 .  (error return if not)
C    *5. NKNOTS = NDIM+4 = 2*N+4 .  (error return if not)
C    *6. T(2*k+1) = T(2*k) = X(k), k=1,...,N .  (not checked)
C
C       * Indicates this applies only if KNOTYP.LT.0 .
C
C *Portability:
C     Argument INCFD is used only to cause the compiler to generate
C     efficient code for the subscript expressions (1+(I-1)*INCFD) .
C     The normal usage, in which DPCHBS is called with one-dimensional
C     arrays F and D, is probably non-Fortran 77, in the strict sense,
C     but it works on all systems on which DPCHBS has been tested.
C
C *See Also:
C     PCHIC, PCHIM, or PCHSP can be used to determine an interpolating
C        PCH function from a set of data.
C     The B-spline routine DBVALU can be used to evaluate the
C        B-representation that is output by DPCHBS.
C        (See BSPDOC for more information.)
C
C***REFERENCES  F. N. Fritsch, "Representations for parametric cubic
C                 splines," Computer Aided Geometric Design 6 (1989),
C                 pp.79-82.
C***ROUTINES CALLED  DPCHKT, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   870701  DATE WRITTEN
C   900405  Converted Fortran to upper case.
C   900405  Removed requirement that X be dimensioned N+1.
C   900406  Modified to make PCHKT a subsidiary routine to simplify
C           usage.  In the process, added argument INCFD to be com-
C           patible with the rest of PCHIP.
C   900410  Converted prologue to SLATEC 4.0 format.
C   900410  Added calls to XERMSG and changed constant 3. to 3 to
C           reduce single/double differences.
C   900411  Added reference.
C   900430  Produced double precision version.
C   900501  Corrected declarations.
C   930317  Minor cosmetic changes.  (FNF)
C   930514  Corrected problems with dimensioning of arguments and
C           clarified DESCRIPTION.  (FNF)
C   930604  Removed  NKNOTS from DPCHKT call list.  (FNF)
C***END PROLOGUE  DPCHBS
C
C*Internal Notes:
C
C**End
C
C  Declare arguments.
C
      INTEGER  N, INCFD, KNOTYP, NKNOTS, NDIM, KORD, IERR
      DOUBLE PRECISION  X(*), F(INCFD,*), D(INCFD,*), T(*), BCOEF(*)
C
C  Declare local variables.
C
      INTEGER  K, KK
      DOUBLE PRECISION  DOV3, HNEW, HOLD
      CHARACTER*8  LIBNAM, SUBNAM
C***FIRST EXECUTABLE STATEMENT  DPCHBS
C
C  Initialize.
C
      NDIM = 2*N
      KORD = 4
      IERR = 0
      LIBNAM = 'SLATEC'
      SUBNAM = 'DPCHBS'
C
C  Check argument validity.  Set up knot sequence if OK.
C
      IF ( KNOTYP.GT.2 )  THEN
         IERR = -1
         CALL XERMSG (LIBNAM, SUBNAM, 'KNOTYP GREATER THAN 2', IERR, 1)
         RETURN
      ENDIF
      IF ( KNOTYP.LT.0 )  THEN
         IF ( NKNOTS.NE.NDIM+4 )  THEN
            IERR = -2
            CALL XERMSG (LIBNAM, SUBNAM,
     *                    'KNOTYP.LT.0 AND NKNOTS.NE.(2*N+4)', IERR, 1)
            RETURN
         ENDIF
      ELSE
C          Set up knot sequence.
         NKNOTS = NDIM + 4
         CALL DPCHKT (N, X, KNOTYP, T)
      ENDIF
C
C  Compute B-spline coefficients.
C
      HNEW = T(3) - T(1)
      DO 40  K = 1, N
         KK = 2*K
         HOLD = HNEW
C          The following requires mixed mode arithmetic.
         DOV3 = D(1,K)/3
         BCOEF(KK-1) = F(1,K) - HOLD*DOV3
C          The following assumes T(2*K+1) = X(K).
         HNEW = T(KK+3) - T(KK+1)
         BCOEF(KK) = F(1,K) + HNEW*DOV3
   40 CONTINUE
C
C  Terminate.
C
      RETURN
C------------- LAST LINE OF DPCHBS FOLLOWS -----------------------------
      END
