*DECK DDRIV3
      SUBROUTINE DDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS,
     8   EWT, IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8   LENW, IWORK, LENIW, JACOBN, FA, NDE, MXSTEP, G, USERS, IERFLG)
C***BEGIN PROLOGUE  DDRIV3
C***PURPOSE  The function of DDRIV3 is to solve N ordinary differential
C            equations of the form dY(I)/dT = F(Y(I),T), given the
C            initial conditions Y(I) = YI.  The program has options to
C            allow the solution of both stiff and non-stiff differential
C            equations.  Other important options are available.  DDRIV3
C            uses double precision arithmetic.
C***LIBRARY   SLATEC (SDRIVE)
C***CATEGORY  I1A2, I1A1B
C***TYPE      DOUBLE PRECISION (SDRIV3-S, DDRIV3-D, CDRIV3-C)
C***KEYWORDS  DOUBLE PRECISION, GEAR'S METHOD, INITIAL VALUE PROBLEMS,
C             ODE, ORDINARY DIFFERENTIAL EQUATIONS, SDRIVE, STIFF
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             Mail Stop D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***DESCRIPTION
C
C  I.  ABSTRACT  .......................................................
C
C    The primary function of DDRIV3 is to solve N ordinary differential
C    equations of the form dY(I)/dT = F(Y(I),T), given the initial
C    conditions Y(I) = YI.  The program has options to allow the
C    solution of both stiff and non-stiff differential equations.  In
C    addition, DDRIV3 may be used to solve:
C      1. The initial value problem, A*dY(I)/dT = F(Y(I),T), where A is
C         a non-singular matrix depending on Y and T.
C      2. The hybrid differential/algebraic initial value problem,
C         A*dY(I)/dT = F(Y(I),T), where A is a vector (whose values may
C         depend upon Y and T) some of whose components will be zero
C         corresponding to those equations which are algebraic rather
C         than differential.
C    DDRIV3 is to be called once for each output point of T.
C
C  II.  PARAMETERS  ....................................................
C       (REMEMBER--To run DDRIV3 correctly in double precision, ALL
C       non-integer arguments in the call sequence, including
C       arrays, MUST be declared double precision.)
C
C    The user should use parameter names in the call sequence of DDRIV3
C    for those quantities whose value may be altered by DDRIV3.  The
C    parameters in the call sequence are:
C
C    N      = (Input) The number of dependent functions whose solution
C             is desired.  N must not be altered during a problem.
C
C    T      = The independent variable.  On input for the first call, T
C             is the initial point.  On output, T is the point at which
C             the solution is given.
C
C    Y      = The vector of dependent variables.  Y is used as input on
C             the first call, to set the initial values.  On output, Y
C             is the computed solution vector.  This array Y is passed
C             in the call sequence of the user-provided routines F,
C             JACOBN, FA, USERS, and G.  Thus parameters required by
C             those routines can be stored in this array in components
C             N+1 and above.  (Note: Changes by the user to the first
C             N components of this array will take effect only after a
C             restart, i.e., after setting NSTATE to 1 .)
C
C    F      = A subroutine supplied by the user.  The name must be
C             declared EXTERNAL in the user's calling program.  This
C             subroutine is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   DOUBLE PRECISION Y(*), YDOT(*)
C                     .
C                     .
C                   YDOT(1) = ...
C                     .
C                     .
C                   YDOT(N) = ...
C                   END (Sample)
C             This computes YDOT = F(Y,T), the right hand side of the
C             differential equations.  Here Y is a vector of length at
C             least N.  The actual length of Y is determined by the
C             user's declaration in the program which calls DDRIV3.
C             Thus the dimensioning of Y in F, while required by FORTRAN
C             convention, does not actually allocate any storage.  When
C             this subroutine is called, the first N components of Y are
C             intermediate approximations to the solution components.
C             The user should not alter these values.  Here YDOT is a
C             vector of length N.  The user should only compute YDOT(I)
C             for I from 1 to N.  Normally a return from F passes
C             control back to  DDRIV3.  However, if the user would like
C             to abort the calculation, i.e., return control to the
C             program which calls DDRIV3, he should set N to zero.
C             DDRIV3 will signal this by returning a value of NSTATE
C             equal to 6 .  Altering the value of N in F has no effect
C             on the value of N in the call sequence of DDRIV3.
C
C    NSTATE = An integer describing the status of integration.  The
C             meaning of NSTATE is as follows:
C               1  (Input) Means the first call to the routine.  This
C                  value must be set by the user.  On all subsequent
C                  calls the value of NSTATE should be tested by the
C                  user, but must not be altered.  (As a convenience to
C                  the user who may wish to put out the initial
C                  conditions, DDRIV3 can be called with NSTATE=1, and
C                  TOUT=T.  In this case the program will return with
C                  NSTATE unchanged, i.e., NSTATE=1.)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  MXSTEP steps without reaching TOUT.  The user can
C                  continue the integration by simply calling DDRIV3
C                  again.
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling DDRIV3
C                  again.
C               5  (Output) A root was found at a point less than TOUT.
C                  The user can continue the integration toward TOUT by
C                  simply calling DDRIV3 again.
C               6  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE F.
C               7  (Output)(Unsuccessful) N has been set to zero in
C                  FUNCTION G.  See description of G below.
C               8  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE JACOBN.  See description of JACOBN below.
C               9  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE FA.  See description of FA below.
C              10  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE USERS.  See description of USERS below.
C              11  (Output)(Successful) For NTASK = 2 or 3, T is beyond
C                  TOUT.  The solution was obtained by interpolation.
C                  The user can continue the integration by simply
C                  advancing TOUT and calling DDRIV3 again.
C              12  (Output)(Unsuccessful) The solution could not be
C                  obtained.  The value of IERFLG (see description
C                  below) for a "Recoverable" situation indicates the
C                  type of difficulty encountered: either an illegal
C                  value for a parameter or an inability to continue the
C                  solution.  For this condition the user should take
C                  corrective action and reset NSTATE to 1 before
C                  calling DDRIV3 again.  Otherwise the program will
C                  terminate the run.
C
C    TOUT   = (Input) The point at which the solution is desired.  The
C             position of TOUT relative to T on the first call
C             determines the direction of integration.
C
C    NTASK  = (Input) An index specifying the manner of returning the
C             solution, according to the following:
C               NTASK = 1  Means DDRIV3 will integrate past TOUT and
C                          interpolate the solution.  This is the most
C                          efficient mode.
C               NTASK = 2  Means DDRIV3 will return the solution after
C                          each internal integration step, or at TOUT,
C                          whichever comes first.  In the latter case,
C                          the program integrates exactly to TOUT.
C               NTASK = 3  Means DDRIV3 will adjust its internal step to
C                          reach TOUT exactly (useful if a singularity
C                          exists beyond TOUT.)
C
C    NROOT  = (Input) The number of equations whose roots are desired.
C             If NROOT is zero, the root search is not active.  This
C             option is useful for obtaining output at points which are
C             not known in advance, but depend upon the solution, e.g.,
C             when some solution component takes on a specified value.
C             The root search is carried out using the user-written
C             function G (see description of G below.)  DDRIV3 attempts
C             to find the value of T at which one of the equations
C             changes sign.  DDRIV3 can find at most one root per
C             equation per internal integration step, and will then
C             return the solution either at TOUT or at a root, whichever
C             occurs first in the direction of integration.  The initial
C             point is never reported as a root.  The index of the
C             equation whose root is being reported is stored in the
C             sixth element of IWORK.
C             NOTE: NROOT is never altered by this program.
C
C    EPS    = On input, the requested relative accuracy in all solution
C             components.  EPS = 0 is allowed.  On output, the adjusted
C             relative accuracy if the input value was too small.  The
C             value of EPS should be set as large as is reasonable,
C             because the amount of work done by DDRIV3 increases as EPS
C             decreases.
C
C    EWT    = (Input) Problem zero, i.e., the smallest, nonzero,
C             physically meaningful value for the solution.  (Array,
C             possibly of length one.  See following description of
C             IERROR.)  Setting EWT smaller than necessary can adversely
C             affect the running time.
C
C    IERROR = (Input) Error control indicator.  A value of 3 is
C             suggested for most problems.  Other choices and detailed
C             explanations of EWT and IERROR are given below for those
C             who may need extra flexibility.
C
C             These last three input quantities EPS, EWT and IERROR
C             control the accuracy of the computed solution.  EWT and
C             IERROR are used internally to compute an array YWT.  One
C             step error estimates divided by YWT(I) are kept less than
C             EPS in root mean square norm.
C                 IERROR (Set by the user) =
C                 1  Means YWT(I) = 1. (Absolute error control)
C                                   EWT is ignored.
C                 2  Means YWT(I) = ABS(Y(I)),  (Relative error control)
C                                   EWT is ignored.
C                 3  Means YWT(I) = MAX(ABS(Y(I)), EWT(1)).
C                 4  Means YWT(I) = MAX(ABS(Y(I)), EWT(I)).
C                    This choice is useful when the solution components
C                    have differing scales.
C                 5  Means YWT(I) = EWT(I).
C             If IERROR is 3, EWT need only be dimensioned one.
C             If IERROR is 4 or 5, the user must dimension EWT at least
C             N, and set its values.
C
C    MINT   = (Input) The integration method indicator.
C               MINT = 1  Means the Adams methods, and is used for
C                         non-stiff problems.
C               MINT = 2  Means the stiff methods of Gear (i.e., the
C                         backward differentiation formulas), and is
C                         used for stiff problems.
C               MINT = 3  Means the program dynamically selects the
C                         Adams methods when the problem is non-stiff
C                         and the Gear methods when the problem is
C                         stiff.  When using the Adams methods, the
C                         program uses a value of MITER=0; when using
C                         the Gear methods, the program uses the value
C                         of MITER provided by the user.  Only a value
C                         of IMPL = 0 and a value of MITER = 1, 2, 4, or
C                         5 is allowed for this option.  The user may
C                         not alter the value of MINT or MITER without
C                         restarting, i.e., setting NSTATE to 1.
C
C    MITER  = (Input) The iteration method indicator.
C               MITER = 0  Means functional iteration.  This value is
C                          suggested for non-stiff problems.
C               MITER = 1  Means chord method with analytic Jacobian.
C                          In this case, the user supplies subroutine
C                          JACOBN (see description below).
C               MITER = 2  Means chord method with Jacobian calculated
C                          internally by finite differences.
C               MITER = 3  Means chord method with corrections computed
C                          by the user-written routine USERS (see
C                          description of USERS below.)  This option
C                          allows all matrix algebra and storage
C                          decisions to be made by the user.  When using
C                          a value of MITER = 3, the subroutine FA is
C                          not required, even if IMPL is not 0.  For
C                          further information on using this option, see
C                          Section IV-E below.
C               MITER = 4  Means the same as MITER = 1 but the A and
C                          Jacobian matrices are assumed to be banded.
C               MITER = 5  Means the same as MITER = 2 but the A and
C                          Jacobian matrices are assumed to be banded.
C
C    IMPL   = (Input) The implicit method indicator.
C               IMPL = 0    Means solving dY(I)/dT = F(Y(I),T).
C               IMPL = 1    Means solving A*dY(I)/dT = F(Y(I),T), non-
C                           singular A (see description of FA below.)
C                           Only MINT = 1 or 2, and MITER = 1, 2, 3, 4,
C                           or 5 are allowed for this option.
C               IMPL = 2,3  Means solving certain systems of hybrid
C                           differential/algebraic equations (see
C                           description of FA below.)  Only MINT = 2 and
C                           MITER = 1, 2, 3, 4, or 5, are allowed for
C                           this option.
C               The value of IMPL must not be changed during a problem.
C
C    ML     = (Input) The lower half-bandwidth in the case of a banded
C             A or Jacobian matrix.  (I.e., maximum(R-C) for nonzero
C             A(R,C).)
C
C    MU     = (Input) The upper half-bandwidth in the case of a banded
C             A or Jacobian matrix.  (I.e., maximum(C-R).)
C
C    MXORD  = (Input) The maximum order desired. This is .LE. 12 for
C             the Adams methods and .LE. 5 for the Gear methods.  Normal
C             value is 12 and 5, respectively.  If MINT is 3, the
C             maximum order used will be MIN(MXORD, 12) when using the
C             Adams methods, and MIN(MXORD, 5) when using the Gear
C             methods.  MXORD must not be altered during a problem.
C
C    HMAX   = (Input) The maximum magnitude of the step size that will
C             be used for the problem.  This is useful for ensuring that
C             important details are not missed.  If this is not the
C             case, a large value, such as the interval length, is
C             suggested.
C
C    WORK
C    LENW   = (Input)
C             WORK is an array of LENW double precision words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       DOUBLE PRECISION WORK(...)
C             The following table gives the required minimum value for
C             the length of WORK, depending on the value of IMPL and
C             MITER.  LENW should be set to the value used.  The
C             contents of WORK should not be disturbed between calls to
C             DDRIV3.
C
C      IMPL =   0            1               2             3
C              ---------------------------------------------------------
C MITER =  0   (MXORD+4)*N   Not allowed     Not allowed   Not allowed
C              + 2*NROOT
C              + 250
C
C         1,2  N*N +         2*N*N +         N*N +         N*(N + NDE)
C              (MXORD+5)*N   (MXORD+5)*N     (MXORD+6)*N   + (MXORD+5)*N
C              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT
C              + 250         + 250           + 250         + 250
C
C          3   (MXORD+4)*N   (MXORD+4)*N     (MXORD+4)*N   (MXORD+4)*N
C              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT
C              + 250         + 250           + 250         + 250
C
C         4,5  (2*ML+MU+1)   2*(2*ML+MU+1)   (2*ML+MU+1)   (2*ML+MU+1)*
C              *N +          *N +            *N +          (N+NDE) +
C              (MXORD+5)*N   (MXORD+5)*N     (MXORD+6)*N   + (MXORD+5)*N
C              + 2*NROOT     + 2*NROOT       + 2*NROOT     + 2*NROOT
C              + 250         + 250           + 250         + 250
C              ---------------------------------------------------------
C
C    IWORK
C    LENIW  = (Input)
C             IWORK is an integer array of length LENIW used internally
C             for temporary storage.  The user must allocate space for
C             this array in the calling program by a statement such as
C                       INTEGER IWORK(...)
C             The length of IWORK should be at least
C               50      if MITER is 0 or 3, or
C               N+50    if MITER is 1, 2, 4, or 5, or MINT is 3,
C             and LENIW should be set to the value used.  The contents
C             of IWORK should not be disturbed between calls to DDRIV3.
C
C    JACOBN = A subroutine supplied by the user, if MITER is 1 or 4.
C             If this is the case, the name must be declared EXTERNAL in
C             the user's calling program.  Given a system of N
C             differential equations, it is meaningful to speak about
C             the partial derivative of the I-th right hand side with
C             respect to the J-th dependent variable.  In general there
C             are N*N such quantities.  Often however the equations can
C             be ordered so that the I-th differential equation only
C             involves dependent variables with index near I, e.g., I+1,
C             I-2.  Such a system is called banded.  If, for all I, the
C             I-th equation depends on at most the variables
C               Y(I-ML), Y(I-ML+1), ... , Y(I), Y(I+1), ... , Y(I+MU)
C             then we call ML+MU+1 the bandwidth of the system.  In a
C             banded system many of the partial derivatives above are
C             automatically zero.  For the cases MITER = 1, 2, 4, and 5,
C             some of these partials are needed.  For the cases
C             MITER = 2 and 5 the necessary derivatives are
C             approximated numerically by DDRIV3, and we only ask the
C             user to tell DDRIV3 the value of ML and MU if the system
C             is banded.  For the cases MITER = 1 and 4 the user must
C             derive these partials algebraically and encode them in
C             subroutine JACOBN.  By computing these derivatives the
C             user can often save 20-30 per cent of the computing time.
C             Usually, however, the accuracy is not much affected and
C             most users will probably forego this option.  The optional
C             user-written subroutine JACOBN has the form:
C                   SUBROUTINE JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
C                   DOUBLE PRECISION Y(*), DFDY(MATDIM,*)
C                     .
C                     .
C                     Calculate values of DFDY
C                     .
C                     .
C                   END (Sample)
C             Here Y is a vector of length at least N.  The actual
C             length of Y is determined by the user's declaration in the
C             program which calls DDRIV3.  Thus the dimensioning of Y in
C             JACOBN, while required by FORTRAN convention, does not
C             actually allocate any storage.  When this subroutine is
C             called, the first N components of Y are intermediate
C             approximations to the solution components.  The user
C             should not alter these values.  If the system is not
C             banded (MITER=1), the partials of the I-th equation with
C             respect to the J-th dependent function are to be stored in
C             DFDY(I,J).  Thus partials of the I-th equation are stored
C             in the I-th row of DFDY.  If the system is banded
C             (MITER=4), then the partials of the I-th equation with
C             respect to Y(J) are to be stored in DFDY(K,J), where
C             K=I-J+MU+1 .  Normally a return from JACOBN passes control
C             back to DDRIV3.  However, if the user would like to abort
C             the calculation, i.e., return control to the program which
C             calls DDRIV3, he should set N to zero.  DDRIV3 will signal
C             this by returning a value of NSTATE equal to +8(-8).
C             Altering the value of N in JACOBN has no effect on the
C             value of N in the call sequence of DDRIV3.
C
C    FA     = A subroutine supplied by the user if IMPL is not zero, and
C             MITER is not 3.  If so, the name must be declared EXTERNAL
C             in the user's calling program.  This subroutine computes
C             the array A, where A*dY(I)/dT = F(Y(I),T).
C             There are three cases:
C
C               IMPL=1.
C               Subroutine FA is of the form:
C                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
C                   DOUBLE PRECISION Y(*), A(MATDIM,*)
C                     .
C                     .
C                     Calculate ALL values of A
C                     .
C                     .
C                   END (Sample)
C               In this case A is assumed to be a nonsingular matrix,
C               with the same structure as DFDY (see JACOBN description
C               above).  Programming considerations prevent complete
C               generality.  If MITER is 1 or 2, A is assumed to be full
C               and the user must compute and store all values of
C               A(I,J), I,J=1, ... ,N.  If MITER is 4 or 5, A is assumed
C               to be banded with lower and upper half bandwidth ML and
C               MU.  The left hand side of the I-th equation is a linear
C               combination of dY(I-ML)/dT, dY(I-ML+1)/dT, ... ,
C               dY(I)/dT, ... , dY(I+MU-1)/dT, dY(I+MU)/dT.  Thus in the
C               I-th equation, the coefficient of dY(J)/dT is to be
C               stored in A(K,J), where K=I-J+MU+1.
C               NOTE: The array A will be altered between calls to FA.
C
C               IMPL=2.
C               Subroutine FA is of the form:
C                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
C                   DOUBLE PRECISION Y(*), A(*)
C                     .
C                     .
C                     Calculate non-zero values of A(1),...,A(NDE)
C                     .
C                     .
C                   END (Sample)
C               In this case it is assumed that the system is ordered by
C               the user so that the differential equations appear
C               first, and the algebraic equations appear last.  The
C               algebraic equations must be written in the form:
C               0 = F(Y(I),T).  When using this option it is up to the
C               user to provide initial values for the Y(I) that satisfy
C               the algebraic equations as well as possible.  It is
C               further assumed that A is a vector of length NDE.  All
C               of the components of A, which may depend on T, Y(I),
C               etc., must be set by the user to non-zero values.
C
C               IMPL=3.
C               Subroutine FA is of the form:
C                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
C                   DOUBLE PRECISION Y(*), A(MATDIM,*)
C                     .
C                     .
C                     Calculate ALL values of A
C                     .
C                     .
C                   END (Sample)
C               In this case A is assumed to be a nonsingular NDE by NDE
C               matrix with the same structure as DFDY (see JACOBN
C               description above).  Programming considerations prevent
C               complete generality.  If MITER is 1 or 2, A is assumed
C               to be full and the user must compute and store all
C               values of A(I,J), I,J=1, ... ,NDE.  If MITER is 4 or 5,
C               A is assumed to be banded with lower and upper half
C               bandwidths ML and MU.  The left hand side of the I-th
C               equation is a linear combination of dY(I-ML)/dT,
C               dY(I-ML+1)/dT, ... , dY(I)/dT, ... , dY(I+MU-1)/dT,
C               dY(I+MU)/dT.  Thus in the I-th equation, the coefficient
C               of dY(J)/dT is to be stored in A(K,J), where K=I-J+MU+1.
C               It is assumed that the system is ordered by the user so
C               that the differential equations appear first, and the
C               algebraic equations appear last.  The algebraic
C               equations must be written in the form 0 = F(Y(I),T).
C               When using this option it is up to the user to provide
C               initial values for the Y(I) that satisfy the algebraic
C               equations as well as possible.
C               NOTE: For IMPL = 3, the array A will be altered between
C               calls to FA.
C             Here Y is a vector of length at least N.  The actual
C             length of Y is determined by the user's declaration in the
C             program which calls DDRIV3.  Thus the dimensioning of Y in
C             FA, while required by FORTRAN convention, does not
C             actually allocate any storage.  When this subroutine is
C             called, the first N components of Y are intermediate
C             approximations to the solution components.  The user
C             should not alter these values.  FA is always called
C             immediately after calling F, with the same values of T
C             and Y.  Normally a return from FA passes control back to
C             DDRIV3.  However, if the user would like to abort the
C             calculation, i.e., return control to the program which
C             calls DDRIV3, he should set N to zero.  DDRIV3 will signal
C             this by returning a value of NSTATE equal to +9(-9).
C             Altering the value of N in FA has no effect on the value
C             of N in the call sequence of DDRIV3.
C
C    NDE    = (Input) The number of differential equations.  This is
C             required only for IMPL = 2 or 3, with NDE .LT. N.
C
C    MXSTEP = (Input) The maximum number of internal steps allowed on
C             one call to DDRIV3.
C
C    G      = A double precision FORTRAN function supplied by the user
C             if NROOT is not 0.  In this case, the name must be
C             declared EXTERNAL in the user's calling program.  G is
C             repeatedly called with different values of IROOT to obtain
C             the value of each of the NROOT equations for which a root
C             is desired.  G is of the form:
C                   DOUBLE PRECISION FUNCTION G (N, T, Y, IROOT)
C                   DOUBLE PRECISION Y(*)
C                   GO TO (10, ...), IROOT
C              10   G = ...
C                     .
C                     .
C                   END (Sample)
C             Here, Y is a vector of length at least N, whose first N
C             components are the solution components at the point T.
C             The user should not alter these values.  The actual length
C             of Y is determined by the user's declaration in the
C             program which calls DDRIV3.  Thus the dimensioning of Y in
C             G, while required by FORTRAN convention, does not actually
C             allocate any storage.  Normally a return from G passes
C             control back to  DDRIV3.  However, if the user would like
C             to abort the calculation, i.e., return control to the
C             program which calls DDRIV3, he should set N to zero.
C             DDRIV3 will signal this by returning a value of NSTATE
C             equal to +7(-7).  In this case, the index of the equation
C             being evaluated is stored in the sixth element of IWORK.
C             Altering the value of N in G has no effect on the value of
C             N in the call sequence of DDRIV3.
C
C    USERS  = A subroutine supplied by the user, if MITER is 3.
C             If this is the case, the name must be declared EXTERNAL in
C             the user's calling program.  The routine USERS is called
C             by DDRIV3 when certain linear systems must be solved.  The
C             user may choose any method to form, store and solve these
C             systems in order to obtain the solution result that is
C             returned to DDRIV3.  In particular, this allows sparse
C             matrix methods to be used.  The call sequence for this
C             routine is:
C
C                SUBROUTINE USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL,
C               8                  IMPL, N, NDE, IFLAG)
C                DOUBLE PRECISION Y(*), YH(*), YWT(*), SAVE1(*),
C               8     SAVE2(*), T, H, EL
C
C             The input variable IFLAG indicates what action is to be
C             taken.  Subroutine USERS should perform the following
C             operations, depending on the value of IFLAG and IMPL.
C
C               IFLAG = 0
C                 IMPL = 0.  USERS is not called.
C                 IMPL = 1, 2 or 3.  Solve the system A*X = SAVE2,
C                   returning the result in SAVE2.  The array SAVE1 can
C                   be used as a work array.  For IMPL = 1, there are N
C                   components to the system, and for IMPL = 2 or 3,
C                   there are NDE components to the system.
C
C               IFLAG = 1
C                 IMPL = 0.  Compute, decompose and store the matrix
C                   (I - H*EL*J), where I is the identity matrix and J
C                   is the Jacobian matrix of the right hand side.  The
C                   array SAVE1 can be used as a work array.
C                 IMPL = 1, 2 or 3. Compute, decompose and store the
C                   matrix (A - H*EL*J).  The array SAVE1 can be used as
C                   a work array.
C
C               IFLAG = 2
C                 IMPL = 0.   Solve the system
C                     (I - H*EL*J)*X = H*SAVE2 - YH - SAVE1,
C                   returning the result in SAVE2.
C                 IMPL = 1, 2 or 3.  Solve the system
C                   (A - H*EL*J)*X = H*SAVE2 - A*(YH + SAVE1)
C                   returning the result in SAVE2.
C                 The array SAVE1 should not be altered.
C             If IFLAG is 0 and IMPL is 1 or 2 and the matrix A is
C             singular, or if IFLAG is 1 and one of the matrices
C             (I - H*EL*J), (A - H*EL*J) is singular, the INTEGER
C             variable IFLAG is to be set to -1 before RETURNing.
C             Normally a return from USERS passes control back to
C             DDRIV3.  However, if the user would like to abort the
C             calculation, i.e., return control to the program which
C             calls DDRIV3, he should set N to zero.  DDRIV3 will signal
C             this by returning a value of NSTATE equal to +10(-10).
C             Altering the value of N in USERS has no effect on the
C             value of N in the call sequence of DDRIV3.
C
C    IERFLG = An error flag.  The error number associated with a
C             diagnostic message (see Section III-A below) is the same
C             as the corresponding value of IERFLG.  The meaning of
C             IERFLG:
C               0  The routine completed successfully. (No message is
C                  issued.)
C               3  (Warning) The number of steps required to reach TOUT
C                  exceeds MXSTEP.
C               4  (Warning) The value of EPS is too small.
C              11  (Warning) For NTASK = 2 or 3, T is beyond TOUT.
C                  The solution was obtained by interpolation.
C              15  (Warning) The integration step size is below the
C                  roundoff level of T.  (The program issues this
C                  message as a warning but does not return control to
C                  the user.)
C              22  (Recoverable) N is not positive.
C              23  (Recoverable) MINT is less than 1 or greater than 3 .
C              24  (Recoverable) MITER is less than 0 or greater than
C                  5 .
C              25  (Recoverable) IMPL is less than 0 or greater than 3 .
C              26  (Recoverable) The value of NSTATE is less than 1 or
C                  greater than 12 .
C              27  (Recoverable) EPS is less than zero.
C              28  (Recoverable) MXORD is not positive.
C              29  (Recoverable) For MINT = 3, either MITER = 0 or 3, or
C                  IMPL = 0 .
C              30  (Recoverable) For MITER = 0, IMPL is not 0 .
C              31  (Recoverable) For MINT = 1, IMPL is 2 or 3 .
C              32  (Recoverable) Insufficient storage has been allocated
C                  for the WORK array.
C              33  (Recoverable) Insufficient storage has been allocated
C                  for the IWORK array.
C              41  (Recoverable) The integration step size has gone
C                  to zero.
C              42  (Recoverable) The integration step size has been
C                  reduced about 50 times without advancing the
C                  solution.  The problem setup may not be correct.
C              43  (Recoverable)  For IMPL greater than 0, the matrix A
C                  is singular.
C             999  (Fatal) The value of NSTATE is 12 .
C
C  III.  OTHER COMMUNICATION TO THE USER  ..............................
C
C    A. The solver communicates to the user through the parameters
C       above.  In addition it writes diagnostic messages through the
C       standard error handling program XERMSG.  A complete description
C       of XERMSG is given in "Guide to the SLATEC Common Mathematical
C       Library" by Kirby W. Fong et al..  At installations which do not
C       have this error handling package the short but serviceable
C       routine, XERMSG, available with this package, can be used.  That
C       program uses the file named OUTPUT to transmit messages.
C
C    B. The first three elements of WORK and the first five elements of
C       IWORK will contain the following statistical data:
C         AVGH     The average step size used.
C         HUSED    The step size last used (successfully).
C         AVGORD   The average order used.
C         IMXERR   The index of the element of the solution vector that
C                  contributed most to the last error test.
C         NQUSED   The order last used (successfully).
C         NSTEP    The number of steps taken since last initialization.
C         NFE      The number of evaluations of the right hand side.
C         NJE      The number of evaluations of the Jacobian matrix.
C
C  IV.  REMARKS  .......................................................
C
C    A. Other routines used:
C         DDNTP, DDZRO, DDSTP, DDNTL, DDPST, DDCOR, DDCST,
C         DDPSC, and DDSCL;
C         DGEFA, DGESL, DGBFA, DGBSL, and DNRM2 (from LINPACK)
C         D1MACH (from the Bell Laboratories Machine Constants Package)
C         XERMSG (from the SLATEC Common Math Library)
C       The last seven routines above, not having been written by the
C       present authors, are not explicitly part of this package.
C
C    B. On any return from DDRIV3 all information necessary to continue
C       the calculation is contained in the call sequence parameters,
C       including the work arrays.  Thus it is possible to suspend one
C       problem, integrate another, and then return to the first.
C
C    C. If this package is to be used in an overlay situation, the user
C       must declare in the primary overlay the variables in the call
C       sequence to DDRIV3.
C
C    D. Changing parameters during an integration.
C       The value of NROOT, EPS, EWT, IERROR, MINT, MITER, or HMAX may
C       be altered by the user between calls to DDRIV3.  For example, if
C       too much accuracy has been requested (the program returns with
C       NSTATE = 4 and an increased value of EPS) the user may wish to
C       increase EPS further.  In general, prudence is necessary when
C       making changes in parameters since such changes are not
C       implemented until the next integration step, which is not
C       necessarily the next call to DDRIV3.  This can happen if the
C       program has already integrated to a point which is beyond the
C       new point TOUT.
C
C    E. As the price for complete control of matrix algebra, the DDRIV3
C       USERS option puts all responsibility for Jacobian matrix
C       evaluation on the user.  It is often useful to approximate
C       numerically all or part of the Jacobian matrix.  However this
C       must be done carefully.  The FORTRAN sequence below illustrates
C       the method we recommend.  It can be inserted directly into
C       subroutine USERS to approximate Jacobian elements in rows I1
C       to I2 and columns J1 to J2.
C              DOUBLE PRECISION DFDY(N,N), EPSJ, H, R, D1MACH,
C             8     SAVE1(N), SAVE2(N), T, UROUND, Y(N), YJ, YWT(N)
C              UROUND = D1MACH(4)
C              EPSJ = SQRT(UROUND)
C              DO 30 J = J1,J2
C                R = EPSJ*MAX(ABS(YWT(J)), ABS(Y(J)))
C                IF (R .EQ. 0.D0) R = YWT(J)
C                YJ = Y(J)
C                Y(J) = Y(J) + R
C                CALL F (N, T, Y, SAVE1)
C                IF (N .EQ. 0) RETURN
C                Y(J) = YJ
C                DO 20 I = I1,I2
C         20       DFDY(I,J) = (SAVE1(I) - SAVE2(I))/R
C         30     CONTINUE
C       Many problems give rise to structured sparse Jacobians, e.g.,
C       block banded.  It is possible to approximate them with fewer
C       function evaluations than the above procedure uses; see Curtis,
C       Powell and Reid, J. Inst. Maths Applics, (1974), Vol. 13,
C       pp. 117-119.
C
C    F. When any of the routines JACOBN, FA, G, or USERS, is not
C       required, difficulties associated with unsatisfied externals can
C       be avoided by using the name of the routine which calculates the
C       right hand side of the differential equations in place of the
C       corresponding name in the call sequence of DDRIV3.
C
C***REFERENCES  C. W. Gear, Numerical Initial Value Problems in
C                 Ordinary Differential Equations, Prentice-Hall, 1971.
C***ROUTINES CALLED  D1MACH, DDNTP, DDSTP, DDZRO, DGBFA, DGBSL, DGEFA,
C                    DGESL, DNRM2, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  DDRIV3
      EXTERNAL F, JACOBN, FA, G, USERS
      DOUBLE PRECISION AE, BIG, EPS, EWT(*), G, GLAST, GNOW, H, HMAX,
     8     HSIGN, HUSED, NROUND, RE, D1MACH, SIZE, DNRM2, SUM, T, TLAST,
     8     TOUT, TROOT, UROUND, WORK(*), Y(*)
      INTEGER I, IA, IAVGH, IAVGRD, ICNVRG, IDFDY, IEL, IERFLG, IERROR,
     8        IFAC, IFLAG, IGNOW, IH, IHMAX, IHOLD, IHSIGN, IHUSED,
     8        IJROOT, IJSTPL, IJTASK, IMNT, IMNTLD, IMPL, IMTR, IMTRLD,
     8        IMTRSV, IMXERR, IMXORD, IMXRDS, INDMXR, INDPRT, INDPVT,
     8        INDTRT, INFE, INFO, INJE, INQ, INQUSE, INROOT, INRTLD,
     8        INSTEP, INWAIT, IRC, IRMAX, IROOT, IMACH1, IMACH4, ISAVE1,
     8        ISAVE2, IT, ITOUT, ITQ, ITREND, ITROOT, IWORK(*), IYH,
     8        IYWT, J, JSTATE, JTROOT, LENCHK, LENIW, LENW, LIWCHK,
     8        MATDIM, MAXORD, MINT, MITER, ML, MU, MXORD, MXSTEP, N,
     8        NDE, NDECOM, NPAR, NROOT, NSTATE, NSTEPL, NTASK
      LOGICAL CONVRG
      CHARACTER INTGR1*8, INTGR2*8, RL1*16, RL2*16
      PARAMETER(NROUND = 20.D0)
      PARAMETER(IAVGH = 1, IHUSED = 2, IAVGRD = 3,
     8          IEL = 4, IH = 160, IHMAX = 161, IHOLD = 162,
     8          IHSIGN = 163, IRC = 164, IRMAX = 165, IT = 166,
     8          ITOUT = 167, ITQ = 168, ITREND = 204, IMACH1 = 205,
     8          IMACH4 = 206, IYH = 251,
     8          INDMXR = 1, INQUSE = 2, INSTEP = 3, INFE = 4, INJE = 5,
     8          INROOT = 6, ICNVRG = 7, IJROOT = 8, IJTASK = 9,
     8          IMNTLD = 10, IMTRLD = 11, INQ = 12, INRTLD = 13,
     8          INDTRT = 14, INWAIT = 15, IMNT = 16, IMTRSV = 17,
     8          IMTR = 18, IMXRDS = 19, IMXORD = 20, INDPRT = 21,
     8          IJSTPL = 22, INDPVT = 51)
C***FIRST EXECUTABLE STATEMENT  DDRIV3
      IF (NSTATE .EQ. 12) THEN
        IERFLG = 999
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'Illegal input.  The value of NSTATE is 12 .', IERFLG, 2)
        RETURN
      ELSE IF (NSTATE .LT. 1 .OR. NSTATE .GT. 12) THEN
        WRITE(INTGR1, '(I8)') NSTATE
        IERFLG = 26
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'Illegal input.  Improper value for NSTATE(= '//INTGR1//').',
     8  IERFLG, 1)
        NSTATE = 12
        RETURN
      END IF
      NPAR = N
      IF (EPS .LT. 0.D0) THEN
        WRITE(RL1, '(D16.8)') EPS
        IERFLG = 27
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'Illegal input.  EPS, '//RL1//', is negative.', IERFLG, 1)
        NSTATE = 12
        RETURN
      END IF
      IF (N .LE. 0) THEN
        WRITE(INTGR1, '(I8)') N
        IERFLG = 22
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'Illegal input.  Number of equations, '//INTGR1//
     8  ', is not positive.', IERFLG, 1)
        NSTATE = 12
        RETURN
      END IF
      IF (MXORD .LE. 0) THEN
        WRITE(INTGR1, '(I8)') MXORD
        IERFLG = 28
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'Illegal input.  Maximum order, '//INTGR1//
     8  ', is not positive.', IERFLG, 1)
        NSTATE = 12
        RETURN
      END IF
      IF (MINT .LT. 1 .OR. MINT .GT. 3) THEN
        WRITE(INTGR1, '(I8)') MINT
        IERFLG = 23
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'Illegal input.  Improper value for the integration method '//
     8  'flag, '//INTGR1//' .', IERFLG, 1)
        NSTATE = 12
        RETURN
      ELSE IF (MITER .LT. 0 .OR. MITER .GT. 5) THEN
        WRITE(INTGR1, '(I8)') MITER
        IERFLG = 24
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'Illegal input.  Improper value for MITER(= '//INTGR1//').',
     8  IERFLG, 1)
        NSTATE = 12
        RETURN
      ELSE IF (IMPL .LT. 0 .OR. IMPL .GT. 3) THEN
        WRITE(INTGR1, '(I8)') IMPL
        IERFLG = 25
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'Illegal input.  Improper value for IMPL(= '//INTGR1//').',
     8  IERFLG, 1)
        NSTATE = 12
        RETURN
      ELSE IF (MINT .EQ. 3 .AND.
     8  (MITER .EQ. 0 .OR. MITER .EQ. 3 .OR. IMPL .NE. 0)) THEN
        WRITE(INTGR1, '(I8)') MITER
        WRITE(INTGR2, '(I8)') IMPL
        IERFLG = 29
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'Illegal input.  For MINT = 3, the value of MITER, '//INTGR1//
     8  ', and/or IMPL, '//INTGR2//', is not allowed.', IERFLG, 1)
        NSTATE = 12
        RETURN
      ELSE IF ((IMPL .GE. 1 .AND. IMPL .LE. 3) .AND. MITER .EQ. 0) THEN
        WRITE(INTGR1, '(I8)') IMPL
        IERFLG = 30
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'Illegal input.  For MITER = 0, the value of IMPL, '//INTGR1//
     8  ', is not allowed.', IERFLG, 1)
        NSTATE = 12
        RETURN
      ELSE IF ((IMPL .EQ. 2 .OR. IMPL .EQ. 3) .AND. MINT .EQ. 1) THEN
        WRITE(INTGR1, '(I8)') IMPL
        IERFLG = 31
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'Illegal input.  For MINT = 1, the value of IMPL, '//INTGR1//
     8  ', is not allowed.', IERFLG, 1)
        NSTATE = 12
        RETURN
      END IF
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        LIWCHK = INDPVT - 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2 .OR. MITER .EQ. 4 .OR.
     8  MITER .EQ. 5) THEN
        LIWCHK = INDPVT + N - 1
      END IF
      IF (LENIW .LT. LIWCHK) THEN
        WRITE(INTGR1, '(I8)') LIWCHK
        IERFLG = 33
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'Illegal input.  Insufficient storage allocated for the '//
     8  'IWORK array.  Based on the value of the input parameters '//
     8  'involved, the required storage is '//INTGR1//' .', IERFLG, 1)
        NSTATE = 12
        RETURN
      END IF
C                                                Allocate the WORK array
C                                         IYH is the index of YH in WORK
      IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN
        MAXORD = MIN(MXORD, 12)
      ELSE IF (MINT .EQ. 2) THEN
        MAXORD = MIN(MXORD, 5)
      END IF
      IDFDY = IYH + (MAXORD + 1)*N
C                                             IDFDY is the index of DFDY
C
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        IYWT = IDFDY
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IYWT = IDFDY + N*N
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IYWT = IDFDY + (2*ML + MU + 1)*N
      END IF
C                                               IYWT is the index of YWT
      ISAVE1 = IYWT + N
C                                           ISAVE1 is the index of SAVE1
      ISAVE2 = ISAVE1 + N
C                                           ISAVE2 is the index of SAVE2
      IGNOW = ISAVE2 + N
C                                             IGNOW is the index of GNOW
      ITROOT = IGNOW + NROOT
C                                           ITROOT is the index of TROOT
      IFAC = ITROOT + NROOT
C                                               IFAC is the index of FAC
      IF (MITER .EQ. 2 .OR. MITER .EQ. 5 .OR. MINT .EQ. 3) THEN
        IA = IFAC + N
      ELSE
        IA = IFAC
      END IF
C                                                   IA is the index of A
      IF (IMPL .EQ. 0 .OR. MITER .EQ. 3) THEN
        LENCHK = IA - 1
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 1 .OR. MITER .EQ. 2)) THEN
        LENCHK = IA - 1 + N*N
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 4 .OR. MITER .EQ. 5)) THEN
        LENCHK = IA - 1 + (2*ML + MU + 1)*N
      ELSE IF (IMPL .EQ. 2 .AND. MITER .NE. 3) THEN
        LENCHK = IA - 1 + N
      ELSE IF (IMPL .EQ. 3 .AND. (MITER .EQ. 1 .OR. MITER .EQ. 2)) THEN
        LENCHK = IA - 1 + N*NDE
      ELSE IF (IMPL .EQ. 3 .AND. (MITER .EQ. 4 .OR. MITER .EQ. 5)) THEN
        LENCHK = IA - 1 + (2*ML + MU + 1)*NDE
      END IF
      IF (LENW .LT. LENCHK) THEN
        WRITE(INTGR1, '(I8)') LENCHK
        IERFLG = 32
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'Illegal input.  Insufficient storage allocated for the '//
     8  'WORK array.  Based on the value of the input parameters '//
     8  'involved, the required storage is '//INTGR1//' .', IERFLG, 1)
        NSTATE = 12
        RETURN
      END IF
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        MATDIM = 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        MATDIM = N
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        MATDIM = 2*ML + MU + 1
      END IF
      IF (IMPL .EQ. 0 .OR. IMPL .EQ. 1) THEN
        NDECOM = N
      ELSE IF (IMPL .EQ. 2 .OR. IMPL .EQ. 3) THEN
        NDECOM = NDE
      END IF
      IF (NSTATE .EQ. 1) THEN
C                                                  Initialize parameters
        IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN
          IWORK(IMXORD) = MIN(MXORD, 12)
        ELSE IF (MINT .EQ. 2) THEN
          IWORK(IMXORD) = MIN(MXORD, 5)
        END IF
        IWORK(IMXRDS) = MXORD
        IF (MINT .EQ. 1 .OR. MINT .EQ. 2) THEN
          IWORK(IMNT) = MINT
          IWORK(IMTR) = MITER
          IWORK(IMNTLD) = MINT
          IWORK(IMTRLD) = MITER
        ELSE IF (MINT .EQ. 3) THEN
          IWORK(IMNT) = 1
          IWORK(IMTR) = 0
          IWORK(IMNTLD) = IWORK(IMNT)
          IWORK(IMTRLD) = IWORK(IMTR)
          IWORK(IMTRSV) = MITER
        END IF
        WORK(IHMAX) = HMAX
        UROUND = D1MACH (4)
        WORK(IMACH4) = UROUND
        WORK(IMACH1) = D1MACH (1)
        IF (NROOT .NE. 0) THEN
          RE = UROUND
          AE = WORK(IMACH1)
        END IF
        H = (TOUT - T)*(1.D0 - 4.D0*UROUND)
        H = SIGN(MIN(ABS(H), HMAX), H)
        WORK(IH) = H
        HSIGN = SIGN(1.D0, H)
        WORK(IHSIGN) = HSIGN
        IWORK(IJTASK) = 0
        WORK(IAVGH) = 0.D0
        WORK(IHUSED) = 0.D0
        WORK(IAVGRD) = 0.D0
        IWORK(INDMXR) = 0
        IWORK(INQUSE) = 0
        IWORK(INSTEP) = 0
        IWORK(IJSTPL) = 0
        IWORK(INFE) = 0
        IWORK(INJE) = 0
        IWORK(INROOT) = 0
        WORK(IT) = T
        IWORK(ICNVRG) = 0
        IWORK(INDPRT) = 0
C                                                 Set initial conditions
        DO 30 I = 1,N
 30       WORK(I+IYH-1) = Y(I)
        IF (T .EQ. TOUT) RETURN
        GO TO 180
      ELSE
        UROUND = WORK(IMACH4)
        IF (NROOT .NE. 0) THEN
          RE = UROUND
          AE = WORK(IMACH1)
        END IF
      END IF
C                                             On a continuation, check
C                                             that output points have
C                                             been or will be overtaken.
      IF (IWORK(ICNVRG) .EQ. 1) THEN
        CONVRG = .TRUE.
      ELSE
        CONVRG = .FALSE.
      END IF
      T = WORK(IT)
      H = WORK(IH)
      HSIGN = WORK(IHSIGN)
      IF (IWORK(IJTASK) .EQ. 0) GO TO 180
C
C                                   IWORK(IJROOT) flags unreported
C                                   roots, and is set to the value of
C                                   NTASK when a root was last selected.
C                                   It is set to zero when all roots
C                                   have been reported.  IWORK(INROOT)
C                                   contains the index and WORK(ITOUT)
C                                   contains the value of the root last
C                                   selected to be reported.
C                                   IWORK(INRTLD) contains the value of
C                                   NROOT and IWORK(INDTRT) contains
C                                   the value of ITROOT when the array
C                                   of roots was last calculated.
      IF (NROOT .NE. 0) THEN
        IF (IWORK(IJROOT) .GT. 0) THEN
C                                      TOUT has just been reported.
C                                      If TROOT .LE. TOUT, report TROOT.
          IF (NSTATE .NE. 5) THEN
            IF (TOUT*HSIGN .GE. WORK(ITOUT)*HSIGN) THEN
              TROOT = WORK(ITOUT)
              CALL DDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
              T = TROOT
              NSTATE = 5
              IERFLG = 0
              GO TO 580
            END IF
C                                         A root has just been reported.
C                                         Select the next root.
          ELSE
            TROOT = T
            IROOT = 0
            DO 50 I = 1,IWORK(INRTLD)
              JTROOT = I + IWORK(INDTRT) - 1
              IF (WORK(JTROOT)*HSIGN .LE. TROOT*HSIGN) THEN
C
C                                              Check for multiple roots.
C
                IF (WORK(JTROOT) .EQ. WORK(ITOUT) .AND.
     8          I .GT. IWORK(INROOT)) THEN
                  IROOT = I
                  TROOT = WORK(JTROOT)
                  GO TO 60
                END IF
                IF (WORK(JTROOT)*HSIGN .GT. WORK(ITOUT)*HSIGN) THEN
                  IROOT = I
                  TROOT = WORK(JTROOT)
                END IF
              END IF
 50           CONTINUE
 60         IWORK(INROOT) = IROOT
            WORK(ITOUT) = TROOT
            IWORK(IJROOT) = NTASK
            IF (NTASK .EQ. 1) THEN
              IF (IROOT .EQ. 0) THEN
                IWORK(IJROOT) = 0
              ELSE
                IF (TOUT*HSIGN .GE. TROOT*HSIGN) THEN
                  CALL DDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),
     8                        Y)
                  NSTATE = 5
                  T = TROOT
                  IERFLG = 0
                  GO TO 580
                END IF
              END IF
            ELSE IF (NTASK .EQ. 2 .OR. NTASK .EQ. 3) THEN
C
C                                     If there are no more roots, or the
C                                     user has altered TOUT to be less
C                                     than a root, set IJROOT to zero.
C
              IF (IROOT .EQ. 0 .OR. (TOUT*HSIGN .LT. TROOT*HSIGN)) THEN
                IWORK(IJROOT) = 0
              ELSE
                CALL DDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),
     8                      Y)
                NSTATE = 5
                IERFLG = 0
                T = TROOT
                GO TO 580
              END IF
            END IF
          END IF
        END IF
      END IF
C
      IF (NTASK .EQ. 1) THEN
        NSTATE = 2
        IF (T*HSIGN .GE. TOUT*HSIGN) THEN
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          IERFLG = 0
          GO TO 580
        END IF
      ELSE IF (NTASK .EQ. 2) THEN
C                                                      Check if TOUT has
C                                                      been reset .LT. T
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          WRITE(RL1, '(D16.8)') T
          WRITE(RL2, '(D16.8)') TOUT
          IERFLG = 11
          CALL XERMSG('SLATEC', 'DDRIV3',
     8    'While integrating exactly to TOUT, T, '//RL1//
     8    ', was beyond TOUT, '//RL2//' .  Solution obtained by '//
     8    'interpolation.', IERFLG, 0)
          NSTATE = 11
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
C                                   Determine if TOUT has been overtaken
C
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          NSTATE = 2
          IERFLG = 0
          GO TO 560
        END IF
C                                             If there are no more roots
C                                             to report, report T.
        IF (NSTATE .EQ. 5) THEN
          NSTATE = 2
          IERFLG = 0
          GO TO 560
        END IF
        NSTATE = 2
C                                                       See if TOUT will
C                                                       be overtaken.
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.D0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        NSTATE = 2
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          WRITE(RL1, '(D16.8)') T
          WRITE(RL2, '(D16.8)') TOUT
          IERFLG = 11
          CALL XERMSG('SLATEC', 'DDRIV3',
     8    'While integrating exactly to TOUT, T, '//RL1//
     8    ', was beyond TOUT, '//RL2//' .  Solution obtained by '//
     8    'interpolation.', IERFLG, 0)
          NSTATE = 11
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          IERFLG = 0
          GO TO 560
        END IF
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.D0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      END IF
C                         Implement changes in MINT, MITER, and/or HMAX.
C
      IF ((MINT .NE. IWORK(IMNTLD) .OR. MITER .NE. IWORK(IMTRLD)) .AND.
     8  MINT .NE. 3 .AND. IWORK(IMNTLD) .NE. 3) IWORK(IJTASK) = -1
      IF (HMAX .NE. WORK(IHMAX)) THEN
        H = SIGN(MIN(ABS(H), HMAX), H)
        IF (H .NE. WORK(IH)) THEN
          IWORK(IJTASK) = -1
          WORK(IH) = H
        END IF
        WORK(IHMAX) = HMAX
      END IF
C
 180  NSTEPL = IWORK(INSTEP)
      DO 190 I = 1,N
 190    Y(I) = WORK(I+IYH-1)
      IF (NROOT .NE. 0) THEN
        DO 200 I = 1,NROOT
          WORK(I+IGNOW-1) = G (NPAR, T, Y, I)
          IF (NPAR .EQ. 0) THEN
            IWORK(INROOT) = I
            NSTATE = 7
            RETURN
          END IF
 200     CONTINUE
      END IF
      IF (IERROR .EQ. 1) THEN
        DO 230 I = 1,N
 230      WORK(I+IYWT-1) = 1.D0
        GO TO 410
      ELSE IF (IERROR .EQ. 5) THEN
        DO 250 I = 1,N
 250      WORK(I+IYWT-1) = EWT(I)
        GO TO 410
      END IF
C                                       Reset YWT array.  Looping point.
 260  IF (IERROR .EQ. 2) THEN
        DO 280 I = 1,N
          IF (Y(I) .EQ. 0.D0) GO TO 290
 280      WORK(I+IYWT-1) = ABS(Y(I))
        GO TO 410
 290    IF (IWORK(IJTASK) .EQ. 0) THEN
          CALL F (NPAR, T, Y, WORK(ISAVE2))
          IF (NPAR .EQ. 0) THEN
            NSTATE = 6
            RETURN
          END IF
          IWORK(INFE) = IWORK(INFE) + 1
          IF (MITER .EQ. 3 .AND. IMPL .NE. 0) THEN
            IFLAG = 0
            CALL USERS (Y, WORK(IYH), WORK(IYWT), WORK(ISAVE1),
     8                  WORK(ISAVE2), T, H, WORK(IEL), IMPL, NPAR,
     8                  NDECOM, IFLAG)
            IF (IFLAG .EQ. -1) GO TO 690
            IF (NPAR .EQ. 0) THEN
              NSTATE = 10
              RETURN
            END IF
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL DGEFA (WORK(IA), MATDIM, N, IWORK(INDPVT), INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGESL (WORK(IA), MATDIM, N, IWORK(INDPVT),
     8                    WORK(ISAVE2), 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              CALL FA (NPAR, T, Y, WORK(IA+ML), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL DGBFA (WORK(IA), MATDIM, N, ML, MU, IWORK(INDPVT),
     8                    INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGBSL (WORK(IA), MATDIM, N, ML, MU, IWORK(INDPVT),
     8                    WORK(ISAVE2), 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
            IF (NPAR .EQ. 0) THEN
              NSTATE = 9
              RETURN
            END IF
            DO 340 I = 1,NDECOM
              IF (WORK(I+IA-1) .EQ. 0.D0) GO TO 690
 340          WORK(I+ISAVE2-1) = WORK(I+ISAVE2-1)/WORK(I+IA-1)
          ELSE IF (IMPL .EQ. 3) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL DGEFA (WORK(IA), MATDIM, NDE, IWORK(INDPVT), INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGESL (WORK(IA), MATDIM, NDE, IWORK(INDPVT),
     8                    WORK(ISAVE2), 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              CALL FA (NPAR, T, Y, WORK(IA+ML), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL DGBFA (WORK(IA), MATDIM, NDE, ML, MU, IWORK(INDPVT),
     8                    INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGBSL (WORK(IA), MATDIM, NDE, ML, MU, IWORK(INDPVT),
     8                    WORK(ISAVE2), 0)
            END IF
          END IF
        END IF
        DO 360 J = I,N
          IF (Y(J) .NE. 0.D0) THEN
            WORK(J+IYWT-1) = ABS(Y(J))
          ELSE
            IF (IWORK(IJTASK) .EQ. 0) THEN
              WORK(J+IYWT-1) = ABS(H*WORK(J+ISAVE2-1))
            ELSE
              WORK(J+IYWT-1) = ABS(WORK(J+IYH+N-1))
            END IF
          END IF
          IF (WORK(J+IYWT-1) .EQ. 0.D0) WORK(J+IYWT-1) = UROUND
 360      CONTINUE
      ELSE IF (IERROR .EQ. 3) THEN
        DO 380 I = 1,N
 380      WORK(I+IYWT-1) = MAX(EWT(1), ABS(Y(I)))
      ELSE IF (IERROR .EQ. 4) THEN
        DO 400 I = 1,N
 400      WORK(I+IYWT-1) = MAX(EWT(I), ABS(Y(I)))
      END IF
C
 410  DO 420 I = 1,N
 420    WORK(I+ISAVE2-1) = Y(I)/WORK(I+IYWT-1)
      SUM = DNRM2(N, WORK(ISAVE2), 1)/SQRT(DBLE(N))
      SUM = MAX(1.D0, SUM)
      IF (EPS .LT. SUM*UROUND) THEN
        EPS = SUM*UROUND*(1.D0 + 10.D0*UROUND)
        WRITE(RL1, '(D16.8)') T
        WRITE(RL2, '(D16.8)') EPS
        IERFLG = 4
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'At T, '//RL1//', the requested accuracy, EPS, was not '//
     8  'obtainable with the machine precision.  EPS has been '//
     8  'increased to '//RL2//' .', IERFLG, 0)
        NSTATE = 4
        GO TO 560
      END IF
      IF (ABS(H) .GE. UROUND*ABS(T)) THEN
        IWORK(INDPRT) = 0
      ELSE IF (IWORK(INDPRT) .EQ. 0) THEN
        WRITE(RL1, '(D16.8)') T
        WRITE(RL2, '(D16.8)') H
        IERFLG = 15
        CALL XERMSG('SLATEC', 'DDRIV3',
     8  'At T, '//RL1//', the step size, '//RL2//', is smaller '//
     8  'than the roundoff level of T.  This may occur if there is '//
     8  'an abrupt change in the right hand side of the '//
     8  'differential equations.', IERFLG, 0)
        IWORK(INDPRT) = 1
      END IF
      IF (NTASK.NE.2) THEN
        IF ((IWORK(INSTEP)-NSTEPL) .EQ. MXSTEP) THEN
          WRITE(RL1, '(D16.8)') T
          WRITE(INTGR1, '(I8)') MXSTEP
          WRITE(RL2, '(D16.8)') TOUT
          IERFLG = 3
          CALL XERMSG('SLATEC', 'DDRIV3',
     8    'At T, '//RL1//', '//INTGR1//' steps have been taken '//
     8    'without reaching TOUT, '//RL2//' .', IERFLG, 0)
          NSTATE = 3
          GO TO 560
        END IF
      END IF
C
C     CALL DDSTP (EPS, F, FA, HMAX, IMPL, IERROR, JACOBN, MATDIM,
C    8            MAXORD, MINT, MITER, ML, MU, N, NDE, YWT, UROUND,
C    8            USERS,  AVGH, AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD,
C    8            NFE, NJE, NQUSED, NSTEP, T, Y, YH,  A, CONVRG,
C    8            DFDY, EL, FAC, HOLD, IPVT, JSTATE, JSTEPL, NQ, NWAIT,
C    8            RC, RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG, MTRSV,
C    8            MXRDSV)
C
      CALL DDSTP (EPS, F, FA, WORK(IHMAX), IMPL, IERROR, JACOBN,
     8            MATDIM, IWORK(IMXORD), IWORK(IMNT), IWORK(IMTR), ML,
     8            MU, NPAR, NDECOM, WORK(IYWT), UROUND, USERS,
     8            WORK(IAVGH), WORK(IAVGRD), WORK(IH), HUSED,
     8            IWORK(IJTASK), IWORK(IMNTLD), IWORK(IMTRLD),
     8            IWORK(INFE), IWORK(INJE), IWORK(INQUSE),
     8            IWORK(INSTEP), WORK(IT), Y, WORK(IYH), WORK(IA),
     8            CONVRG, WORK(IDFDY), WORK(IEL), WORK(IFAC),
     8            WORK(IHOLD), IWORK(INDPVT), JSTATE, IWORK(IJSTPL),
     8            IWORK(INQ), IWORK(INWAIT), WORK(IRC), WORK(IRMAX),
     8            WORK(ISAVE1), WORK(ISAVE2), WORK(ITQ), WORK(ITREND),
     8            MINT, IWORK(IMTRSV), IWORK(IMXRDS))
      T = WORK(IT)
      H = WORK(IH)
      IF (CONVRG) THEN
        IWORK(ICNVRG) = 1
      ELSE
        IWORK(ICNVRG) = 0
      END IF
      GO TO (470, 670, 680, 690, 690, 660, 660, 660, 660, 660), JSTATE
 470  IWORK(IJTASK) = 1
C                                 Determine if a root has been overtaken
      IF (NROOT .NE. 0) THEN
        IROOT = 0
        DO 500 I = 1,NROOT
          GLAST = WORK(I+IGNOW-1)
          GNOW = G (NPAR, T, Y, I)
          IF (NPAR .EQ. 0) THEN
            IWORK(INROOT) = I
            NSTATE = 7
            RETURN
          END IF
          WORK(I+IGNOW-1) = GNOW
          IF (GLAST*GNOW .GT. 0.D0) THEN
            WORK(I+ITROOT-1) = T + H
          ELSE
            IF (GNOW .EQ. 0.D0) THEN
              WORK(I+ITROOT-1) = T
              IROOT = I
            ELSE
              IF (GLAST .EQ. 0.D0) THEN
                WORK(I+ITROOT-1) = T + H
              ELSE
                IF (ABS(HUSED) .GE. UROUND*ABS(T)) THEN
                  TLAST = T - HUSED
                  IROOT = I
                  TROOT = T
                  CALL DDZRO (AE, G, H, NPAR, IWORK(INQ), IROOT, RE, T,
     8                        WORK(IYH), UROUND,  TROOT, TLAST,
     8                        GNOW, GLAST,  Y)
                  DO 480 J = 1,N
 480                Y(J) = WORK(IYH+J-1)
                  IF (NPAR .EQ. 0) THEN
                    IWORK(INROOT) = I
                    NSTATE = 7
                    RETURN
                  END IF
                  WORK(I+ITROOT-1) = TROOT
                ELSE
                  WORK(I+ITROOT-1) = T
                  IROOT = I
                END IF
              END IF
            END IF
          END IF
 500      CONTINUE
        IF (IROOT .EQ. 0) THEN
          IWORK(IJROOT) = 0
C                                                  Select the first root
        ELSE
          IWORK(IJROOT) = NTASK
          IWORK(INRTLD) = NROOT
          IWORK(INDTRT) = ITROOT
          TROOT = T + H
          DO 510 I = 1,NROOT
            IF (WORK(I+ITROOT-1)*HSIGN .LT. TROOT*HSIGN) THEN
              TROOT = WORK(I+ITROOT-1)
              IROOT = I
            END IF
 510        CONTINUE
          IWORK(INROOT) = IROOT
          WORK(ITOUT) = TROOT
          IF (TROOT*HSIGN .LE. TOUT*HSIGN) THEN
            CALL DDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
            NSTATE = 5
            T = TROOT
            IERFLG = 0
            GO TO 580
          END IF
        END IF
      END IF
C                               Test for NTASK condition to be satisfied
      NSTATE = 2
      IF (NTASK .EQ. 1) THEN
        IF (T*HSIGN .LT. TOUT*HSIGN) GO TO 260
        CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
        T = TOUT
        IERFLG = 0
        GO TO 580
C                               TOUT is assumed to have been attained
C                               exactly if T is within twenty roundoff
C                               units of TOUT, relative to MAX(TOUT, T).
C
      ELSE IF (NTASK .EQ. 2) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.D0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.D0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
          GO TO 260
        END IF
      END IF
      IERFLG = 0
C                                      All returns are made through this
C                                      section.  IMXERR is determined.
 560  DO 570 I = 1,N
 570    Y(I) = WORK(I+IYH-1)
 580  IF (IWORK(IJTASK) .EQ. 0) RETURN
      BIG = 0.D0
      IMXERR = 1
      DO  590 I = 1,N
C                                            SIZE = ABS(ERROR(I)/YWT(I))
        SIZE = ABS(WORK(I+ISAVE1-1)/WORK(I+IYWT-1))
        IF (BIG .LT. SIZE) THEN
          BIG = SIZE
          IMXERR = I
        END IF
 590    CONTINUE
      IWORK(INDMXR) = IMXERR
      WORK(IHUSED) = HUSED
      RETURN
C
 660  NSTATE = JSTATE
      RETURN
C                                        Fatal errors are processed here
C
 670  WRITE(RL1, '(D16.8)') T
      IERFLG = 41
      CALL XERMSG('SLATEC', 'DDRIV3',
     8  'At T, '//RL1//', the attempted step size has gone to '//
     8  'zero.  Often this occurs if the problem setup is incorrect.',
     8  IERFLG, 1)
      NSTATE = 12
      RETURN
C
 680  WRITE(RL1, '(D16.8)') T
      IERFLG = 42
      CALL XERMSG('SLATEC', 'DDRIV3',
     8  'At T, '//RL1//', the step size has been reduced about 50 '//
     8  'times without advancing the solution.  Often this occurs '//
     8  'if the problem setup is incorrect.', IERFLG, 1)
      NSTATE = 12
      RETURN
C
 690  WRITE(RL1, '(D16.8)') T
      IERFLG = 43
      CALL XERMSG('SLATEC', 'DDRIV3',
     8  'At T, '//RL1//', while solving A*YDOT = F, A is singular.',
     8  IERFLG, 1)
      NSTATE = 12
      RETURN
      END
