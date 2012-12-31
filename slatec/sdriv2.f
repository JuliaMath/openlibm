*DECK SDRIV2
      SUBROUTINE SDRIV2 (N, T, Y, F, TOUT, MSTATE, NROOT, EPS, EWT,
     8   MINT, WORK, LENW, IWORK, LENIW, G, IERFLG)
C***BEGIN PROLOGUE  SDRIV2
C***PURPOSE  The function of SDRIV2 is to solve N ordinary differential
C            equations of the form dY(I)/dT = F(Y(I),T), given the
C            initial conditions Y(I) = YI.  The program has options to
C            allow the solution of both stiff and non-stiff differential
C            equations.  SDRIV2 uses single precision arithmetic.
C***LIBRARY   SLATEC (SDRIVE)
C***CATEGORY  I1A2, I1A1B
C***TYPE      SINGLE PRECISION (SDRIV2-S, DDRIV2-D, CDRIV2-C)
C***KEYWORDS  GEAR'S METHOD, INITIAL VALUE PROBLEMS, ODE,
C             ORDINARY DIFFERENTIAL EQUATIONS, SDRIVE, SINGLE PRECISION,
C             STIFF
C***AUTHOR  Kahaner, D. K., (NIST)
C             National Institute of Standards and Technology
C             Gaithersburg, MD  20899
C           Sutherland, C. D., (LANL)
C             Mail Stop D466
C             Los Alamos National Laboratory
C             Los Alamos, NM  87545
C***DESCRIPTION
C
C  I.  PARAMETERS  .....................................................
C
C    The user should use parameter names in the call sequence of SDRIV2
C    for those quantities whose value may be altered by SDRIV2.  The
C    parameters in the call sequence are:
C
C    N      = (Input) The number of differential equations.
C
C    T      = The independent variable.  On input for the first call, T
C             is the initial point.  On output, T is the point at which
C             the solution is given.
C
C    Y      = The vector of dependent variables.  Y is used as input on
C             the first call, to set the initial values.  On output, Y
C             is the computed solution vector.  This array Y is passed
C             in the call sequence of the user-provided routines F and
C             G.  Thus parameters required by F and G can be stored in
C             this array in components N+1 and above.  (Note: Changes
C             by the user to the first N components of this array will
C             take effect only after a restart, i.e., after setting
C             MSTATE to +1(-1).)
C
C    F      = A subroutine supplied by the user.  The name must be
C             declared EXTERNAL in the user's calling program.  This
C             subroutine is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   REAL Y(*), YDOT(*)
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
C             user's declaration in the program which calls SDRIV2.
C             Thus the dimensioning of Y in F, while required by FORTRAN
C             convention, does not actually allocate any storage.  When
C             this subroutine is called, the first N components of Y are
C             intermediate approximations to the solution components.
C             The user should not alter these values.  Here YDOT is a
C             vector of length N.  The user should only compute YDOT(I)
C             for I from 1 to N.  Normally a return from F passes
C             control back to  SDRIV2.  However, if the user would like
C             to abort the calculation, i.e., return control to the
C             program which calls SDRIV2, he should set N to zero.
C             SDRIV2 will signal this by returning a value of MSTATE
C             equal to +6(-6).  Altering the value of N in F has no
C             effect on the value of N in the call sequence of SDRIV2.
C
C    TOUT   = (Input) The point at which the solution is desired.
C
C    MSTATE = An integer describing the status of integration.  The user
C             must initialize MSTATE to +1 or -1.  If MSTATE is
C             positive, the routine will integrate past TOUT and
C             interpolate the solution.  This is the most efficient
C             mode.  If MSTATE is negative, the routine will adjust its
C             internal step to reach TOUT exactly (useful if a
C             singularity exists beyond TOUT.)  The meaning of the
C             magnitude of MSTATE:
C               1  (Input) Means the first call to the routine.  This
C                  value must be set by the user.  On all subsequent
C                  calls the value of MSTATE should be tested by the
C                  user.  Unless SDRIV2 is to be reinitialized, only the
C                  sign of MSTATE may be changed by the user.  (As a
C                  convenience to the user who may wish to put out the
C                  initial conditions, SDRIV2 can be called with
C                  MSTATE=+1(-1), and TOUT=T.  In this case the program
C                  will return with MSTATE unchanged, i.e.,
C                  MSTATE=+1(-1).)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  1000 steps without reaching TOUT.  The user can
C                  continue the integration by simply calling SDRIV2
C                  again.  Other than an error in problem setup, the
C                  most likely cause for this condition is trying to
C                  integrate a stiff set of equations with the non-stiff
C                  integrator option. (See description of MINT below.)
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling SDRIV2
C                  again.
C               5  (Output) A root was found at a point less than TOUT.
C                  The user can continue the integration toward TOUT by
C                  simply calling SDRIV2 again.
C               6  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE F.
C               7  (Output)(Unsuccessful) N has been set to zero in
C                  FUNCTION G.  See description of G below.
C               8  (Output)(Successful) For MSTATE negative, T is beyond
C                  TOUT.  The solution was obtained by interpolation.
C                  The user can continue the integration by simply
C                  advancing TOUT and calling SDRIV2 again.
C               9  (Output)(Unsuccessful) The solution could not be
C                  obtained.  The value of IERFLG (see description
C                  below) for a "Recoverable" situation indicates the
C                  type of difficulty encountered: either an illegal
C                  value for a parameter or an inability to continue the
C                  solution.  For this condition the user should take
C                  corrective action and reset MSTATE to +1(-1) before
C                  calling SDRIV2 again.  Otherwise the program will
C                  terminate the run.
C
C    NROOT  = (Input) The number of equations whose roots are desired.
C             If NROOT is zero, the root search is not active.  This
C             option is useful for obtaining output at points which are
C             not known in advance, but depend upon the solution, e.g.,
C             when some solution component takes on a specified value.
C             The root search is carried out using the user-written
C             function G (see description of G below.)  SDRIV2 attempts
C             to find the value of T at which one of the equations
C             changes sign.  SDRIV2 can find at most one root per
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
C             because the amount of work done by SDRIV2 increases as
C             EPS decreases.
C
C    EWT    = (Input) Problem zero, i.e., the smallest physically
C             meaningful value for the solution.  This is used inter-
C             nally to compute an array YWT(I) = MAX(ABS(Y(I)), EWT).
C             One step error estimates divided by YWT(I) are kept less
C             than EPS.  Setting EWT to zero provides pure relative
C             error control.  However, setting EWT smaller than
C             necessary can adversely affect the running time.
C
C    MINT   = (Input) The integration method flag.
C               MINT = 1  Means the Adams methods, and is used for
C                         non-stiff problems.
C               MINT = 2  Means the stiff methods of Gear (i.e., the
C                         backward differentiation formulas), and is
C                         used for stiff problems.
C               MINT = 3  Means the program dynamically selects the
C                         Adams methods when the problem is non-stiff
C                         and the Gear methods when the problem is
C                         stiff.
C             MINT may not be changed without restarting, i.e., setting
C             the magnitude of MSTATE to 1.
C
C    WORK
C    LENW   = (Input)
C             WORK is an array of LENW real words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       REAL WORK(...)
C             The length of WORK should be at least
C               16*N + 2*NROOT + 250         if MINT is 1, or
C               N*N + 10*N + 2*NROOT + 250   if MINT is 2, or
C               N*N + 17*N + 2*NROOT + 250   if MINT is 3,
C             and LENW should be set to the value used.  The contents of
C             WORK should not be disturbed between calls to SDRIV2.
C
C    IWORK
C    LENIW  = (Input)
C             IWORK is an integer array of length LENIW used internally
C             for temporary storage.  The user must allocate space for
C             this array in the calling program by a statement such as
C                       INTEGER IWORK(...)
C             The length of IWORK should be at least
C               50      if MINT is 1, or
C               N+50    if MINT is 2 or 3,
C             and LENIW should be set to the value used.  The contents
C             of IWORK should not be disturbed between calls to SDRIV2.
C
C    G      = A real FORTRAN function supplied by the user
C             if NROOT is not 0.  In this case, the name must be
C             declared EXTERNAL in the user's calling program.  G is
C             repeatedly called with different values of IROOT to
C             obtain the value of each of the NROOT equations for which
C             a root is desired.  G is of the form:
C                   REAL FUNCTION G (N, T, Y, IROOT)
C                   REAL Y(*)
C                   GO TO (10, ...), IROOT
C              10   G = ...
C                     .
C                     .
C                   END (Sample)
C             Here, Y is a vector of length at least N, whose first N
C             components are the solution components at the point T.
C             The user should not alter these values.  The actual length
C             of Y is determined by the user's declaration in the
C             program which calls SDRIV2.  Thus the dimensioning of Y in
C             G, while required by FORTRAN convention, does not actually
C             allocate any storage.  Normally a return from G passes
C             control back to  SDRIV2.  However, if the user would like
C             to abort the calculation, i.e., return control to the
C             program which calls SDRIV2, he should set N to zero.
C             SDRIV2 will signal this by returning a value of MSTATE
C             equal to +7(-7).  In this case, the index of the equation
C             being evaluated is stored in the sixth element of IWORK.
C             Altering the value of N in G has no effect on the value of
C             N in the call sequence of SDRIV2.
C
C    IERFLG = An error flag.  The error number associated with a
C             diagnostic message (see Section II-A below) is the same as
C             the corresponding value of IERFLG.  The meaning of IERFLG:
C               0  The routine completed successfully. (No message is
C                  issued.)
C               3  (Warning) The number of steps required to reach TOUT
C                  exceeds MXSTEP.
C               4  (Warning) The value of EPS is too small.
C              11  (Warning) For MSTATE negative, T is beyond TOUT.
C                  The solution was obtained by interpolation.
C              15  (Warning) The integration step size is below the
C                  roundoff level of T.  (The program issues this
C                  message as a warning but does not return control to
C                  the user.)
C              22  (Recoverable) N is not positive.
C              23  (Recoverable) MINT is less than 1 or greater than 3 .
C              26  (Recoverable) The magnitude of MSTATE is either 0 or
C                  greater than 9 .
C              27  (Recoverable) EPS is less than zero.
C              32  (Recoverable) Insufficient storage has been allocated
C                  for the WORK array.
C              33  (Recoverable) Insufficient storage has been allocated
C                  for the IWORK array.
C              41  (Recoverable) The integration step size has gone
C                  to zero.
C              42  (Recoverable) The integration step size has been
C                  reduced about 50 times without advancing the
C                  solution.  The problem setup may not be correct.
C             999  (Fatal) The magnitude of MSTATE is 9 .
C
C  II.  OTHER COMMUNICATION TO THE USER  ...............................
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
C  III.  REMARKS  ......................................................
C
C    A. On any return from SDRIV2 all information necessary to continue
C       the calculation is contained in the call sequence parameters,
C       including the work arrays.  Thus it is possible to suspend one
C       problem, integrate another, and then return to the first.
C
C    B. If this package is to be used in an overlay situation, the user
C       must declare in the primary overlay the variables in the call
C       sequence to SDRIV2.
C
C    C. When the routine G is not required, difficulties associated with
C       an unsatisfied external can be avoided by using the name of the
C       routine which calculates the right hand side of the differential
C       equations in place of G in the call sequence of SDRIV2.
C
C  IV.  USAGE  .........................................................
C
C               PROGRAM SAMPLE
C               EXTERNAL F
C               PARAMETER(MINT = 1, NROOT = 0, N = ...,
C              8          LENW = 16*N + 2*NROOT + 250, LENIW = 50)
C         C                                 N is the number of equations
C               REAL EPS, EWT, T, TOUT, WORK(LENW), Y(N)
C               INTEGER IWORK(LENIW)
C               OPEN(FILE='TAPE6', UNIT=6, STATUS='NEW')
C         C                                                Initial point
C               T = 0.
C         C                                       Set initial conditions
C               DO 10 I = 1,N
C          10     Y(I) = ...
C               TOUT = T
C               EWT = ...
C               MSTATE = 1
C               EPS = ...
C          20   CALL SDRIV2 (N, T, Y, F, TOUT, MSTATE, NROOT, EPS, EWT,
C              8             MINT, WORK, LENW, IWORK, LENIW, F, IERFLG)
C         C                                 Next to last argument is not
C         C                                    F if rootfinding is used.
C               IF (MSTATE .GT. 2) STOP
C               WRITE(6, 100) TOUT, (Y(I), I=1,N)
C               TOUT = TOUT + 1.
C               IF (TOUT .LE. 10.) GO TO 20
C          100  FORMAT(...)
C               END (Sample)
C
C***REFERENCES  C. W. Gear, Numerical Initial Value Problems in
C                 Ordinary Differential Equations, Prentice-Hall, 1971.
C***ROUTINES CALLED  SDRIV3, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  SDRIV2
      EXTERNAL F, G
      REAL EPS, EWT, EWTCOM(1), G, HMAX, T, TOUT,
     8     WORK(*), Y(*)
      INTEGER IWORK(*)
      INTEGER IERFLG, IERROR, IMPL, LENIW, LENW, MINT, MITER, ML,
     8        MSTATE, MU, MXORD, MXSTEP, N, NDE, NROOT, NSTATE, NTASK
      CHARACTER INTGR1*8
      PARAMETER(IMPL = 0, MXSTEP = 1000)
C***FIRST EXECUTABLE STATEMENT  SDRIV2
      IF (ABS(MSTATE) .EQ. 9) THEN
        IERFLG = 999
        CALL XERMSG('SLATEC', 'SDRIV2',
     8  'Illegal input.  The magnitude of MSTATE IS 9 .',
     8  IERFLG, 2)
        RETURN
      ELSE IF (ABS(MSTATE) .EQ. 0 .OR. ABS(MSTATE) .GT. 9) THEN
        WRITE(INTGR1, '(I8)') MSTATE
        IERFLG = 26
        CALL XERMSG('SLATEC', 'SDRIV2',
     8  'Illegal input.  The magnitude of MSTATE, '//INTGR1//
     8  ' is not in the range 1 to 8 .', IERFLG, 1)
        MSTATE = SIGN(9, MSTATE)
        RETURN
      END IF
      IF (MINT .LT. 1 .OR. MINT .GT. 3) THEN
        WRITE(INTGR1, '(I8)') MINT
        IERFLG = 23
        CALL XERMSG('SLATEC', 'SDRIV2',
     8  'Illegal input.  Improper value for the integration method '//
     8  'flag, '//INTGR1//' .', IERFLG, 1)
        MSTATE = SIGN(9, MSTATE)
        RETURN
      END IF
      IF (MSTATE .GE. 0) THEN
        NSTATE = MSTATE
        NTASK = 1
      ELSE
        NSTATE = - MSTATE
        NTASK = 3
      END IF
      EWTCOM(1) = EWT
      IF (EWT .NE. 0.E0) THEN
        IERROR = 3
      ELSE
        IERROR = 2
      END IF
      IF (MINT .EQ. 1) THEN
        MITER = 0
        MXORD = 12
      ELSE IF (MINT .EQ. 2) THEN
        MITER = 2
        MXORD = 5
      ELSE IF (MINT .EQ. 3) THEN
        MITER = 2
        MXORD = 12
      END IF
      HMAX = 2.E0*ABS(TOUT - T)
      CALL SDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, EWTCOM,
     8             IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8             LENW, IWORK, LENIW, F, F, NDE, MXSTEP, G, F, IERFLG)
      IF (NSTATE .LE. 7) THEN
        MSTATE = SIGN(NSTATE, MSTATE)
      ELSE IF (NSTATE .EQ. 11) THEN
        MSTATE = SIGN(8, MSTATE)
      ELSE IF (NSTATE .GT. 11) THEN
        MSTATE = SIGN(9, MSTATE)
      END IF
      RETURN
      END
