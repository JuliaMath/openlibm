*DECK SDRIV1
      SUBROUTINE SDRIV1 (N, T, Y, F, TOUT, MSTATE, EPS, WORK, LENW,
     8   IERFLG)
C***BEGIN PROLOGUE  SDRIV1
C***PURPOSE  The function of SDRIV1 is to solve N (200 or fewer)
C            ordinary differential equations of the form
C            dY(I)/dT = F(Y(I),T), given the initial conditions
C            Y(I) = YI.  SDRIV1 uses single precision arithmetic.
C***LIBRARY   SLATEC (SDRIVE)
C***CATEGORY  I1A2, I1A1B
C***TYPE      SINGLE PRECISION (SDRIV1-S, DDRIV1-D, CDRIV1-C)
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
C   Version 92.1
C
C  I.  CHOOSING THE CORRECT ROUTINE  ...................................
C
C     SDRIV
C     DDRIV
C     CDRIV
C           These are the generic names for three packages for solving
C           initial value problems for ordinary differential equations.
C           SDRIV uses single precision arithmetic.  DDRIV uses double
C           precision arithmetic.  CDRIV allows complex-valued
C           differential equations, integrated with respect to a single,
C           real, independent variable.
C
C    As an aid in selecting the proper program, the following is a
C    discussion of the important options or restrictions associated with
C    each program:
C
C      A. SDRIV1 should be tried first for those routine problems with
C         no more than 200 differential equations (SDRIV2 and SDRIV3
C         have no such restriction.)  Internally this routine has two
C         important technical defaults:
C           1. Numerical approximation of the Jacobian matrix of the
C              right hand side is used.
C           2. The stiff solver option is used.
C         Most users of SDRIV1 should not have to concern themselves
C         with these details.
C
C      B. SDRIV2 should be considered for those problems for which
C         SDRIV1 is inadequate.  For example, SDRIV1 may have difficulty
C         with problems having zero initial conditions and zero
C         derivatives.  In this case SDRIV2, with an appropriate value
C         of the parameter EWT, should perform more efficiently.  SDRIV2
C         provides three important additional options:
C           1. The nonstiff equation solver (as well as the stiff
C              solver) is available.
C           2. The root-finding option is available.
C           3. The program can dynamically select either the non-stiff
C              or the stiff methods.
C         Internally this routine also defaults to the numerical
C         approximation of the Jacobian matrix of the right hand side.
C
C      C. SDRIV3 is the most flexible, and hence the most complex, of
C         the programs.  Its important additional features include:
C           1. The ability to exploit band structure in the Jacobian
C              matrix.
C           2. The ability to solve some implicit differential
C              equations, i.e., those having the form:
C                   A(Y,T)*dY/dT = F(Y,T).
C           3. The option of integrating in the one step mode.
C           4. The option of allowing the user to provide a routine
C              which computes the analytic Jacobian matrix of the right
C              hand side.
C           5. The option of allowing the user to provide a routine
C              which does all the matrix algebra associated with
C              corrections to the solution components.
C
C  II.  PARAMETERS  ....................................................
C
C    The user should use parameter names in the call sequence of SDRIV1
C    for those quantities whose value may be altered by SDRIV1.  The
C    parameters in the call sequence are:
C
C    N      = (Input) The number of differential equations, N .LE. 200
C
C    T      = The independent variable.  On input for the first call, T
C             is the initial point.  On output, T is the point at which
C             the solution is given.
C
C    Y      = The vector of dependent variables.  Y is used as input on
C             the first call, to set the initial values.  On output, Y
C             is the computed solution vector.  This array Y is passed
C             in the call sequence of the user-provided routine F.  Thus
C             parameters required by F can be stored in this array in
C             components N+1 and above.  (Note: Changes by the user to
C             the first N components of this array will take effect only
C             after a restart, i.e., after setting MSTATE to +1(-1).)
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
C             user's declaration in the program which calls SDRIV1.
C             Thus the dimensioning of Y in F, while required by FORTRAN
C             convention, does not actually allocate any storage.  When
C             this subroutine is called, the first N components of Y are
C             intermediate approximations to the solution components.
C             The user should not alter these values.  Here YDOT is a
C             vector of length N.  The user should only compute YDOT(I)
C             for I from 1 to N.  Normally a return from F passes
C             control back to  SDRIV1.  However, if the user would like
C             to abort the calculation, i.e., return control to the
C             program which calls SDRIV1, he should set N to zero.
C             SDRIV1 will signal this by returning a value of MSTATE
C             equal to +5(-5).  Altering the value of N in F has no
C             effect on the value of N in the call sequence of SDRIV1.
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
C                  user.  Unless SDRIV1 is to be reinitialized, only the
C                  sign of MSTATE may be changed by the user.  (As a
C                  convenience to the user who may wish to put out the
C                  initial conditions, SDRIV1 can be called with
C                  MSTATE=+1(-1), and TOUT=T.  In this case the program
C                  will return with MSTATE unchanged, i.e.,
C                  MSTATE=+1(-1).)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  1000 steps without reaching TOUT.  The user can
C                  continue the integration by simply calling SDRIV1
C                  again.
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling SDRIV1
C                  again.
C               5  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE F.
C               6  (Output)(Successful) For MSTATE negative, T is beyond
C                  TOUT.  The solution was obtained by interpolation.
C                  The user can continue the integration by simply
C                  advancing TOUT and calling SDRIV1 again.
C               7  (Output)(Unsuccessful) The solution could not be
C                  obtained.  The value of IERFLG (see description
C                  below) for a "Recoverable" situation indicates the
C                  type of difficulty encountered: either an illegal
C                  value for a parameter or an inability to continue the
C                  solution.  For this condition the user should take
C                  corrective action and reset MSTATE to +1(-1) before
C                  calling SDRIV1 again.  Otherwise the program will
C                  terminate the run.
C
C    EPS    = On input, the requested relative accuracy in all solution
C             components.  On output, the adjusted relative accuracy if
C             the input value was too small.  The value of EPS should be
C             set as large as is reasonable, because the amount of work
C             done by SDRIV1 increases as EPS decreases.
C
C    WORK
C    LENW   = (Input)
C             WORK is an array of LENW real words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       REAL WORK(...)
C             The length of WORK should be at least N*N + 11*N + 300
C             and LENW should be set to the value used.  The contents of
C             WORK should not be disturbed between calls to SDRIV1.
C
C    IERFLG = An error flag.  The error number associated with a
C             diagnostic message (see Section IV-A below) is the same as
C             the corresponding value of IERFLG.  The meaning of IERFLG:
C               0  The routine completed successfully. (No message is
C                  issued.)
C               3  (Warning) The number of steps required to reach TOUT
C                  exceeds 1000 .
C               4  (Warning) The value of EPS is too small.
C              11  (Warning) For MSTATE negative, T is beyond TOUT.
C                  The solution was obtained by interpolation.
C              15  (Warning) The integration step size is below the
C                  roundoff level of T.  (The program issues this
C                  message as a warning but does not return control to
C                  the user.)
C              21  (Recoverable) N is greater than 200 .
C              22  (Recoverable) N is not positive.
C              26  (Recoverable) The magnitude of MSTATE is either 0 or
C                  greater than 7 .
C              27  (Recoverable) EPS is less than zero.
C              32  (Recoverable) Insufficient storage has been allocated
C                  for the WORK array.
C              41  (Recoverable) The integration step size has gone
C                  to zero.
C              42  (Recoverable) The integration step size has been
C                  reduced about 50 times without advancing the
C                  solution.  The problem setup may not be correct.
C             999  (Fatal) The magnitude of MSTATE is 7 .
C
C  III.  USAGE  ........................................................
C
C                PROGRAM SAMPLE
C                EXTERNAL F
C                REAL ALFA, EPS, T, TOUT
C          C                                N is the number of equations
C                PARAMETER(ALFA = 1.E0, N = 3, LENW = N*N + 11*N + 300)
C                REAL WORK(LENW), Y(N+1)
C          C                                               Initial point
C                T = 0.00001E0
C          C                                      Set initial conditions
C                Y(1) = 10.E0
C                Y(2) = 0.E0
C                Y(3) = 10.E0
C          C                                              Pass parameter
C                Y(4) = ALFA
C                TOUT = T
C                MSTATE = 1
C                EPS = .001E0
C           10   CALL SDRIV1 (N, T, Y, F, TOUT, MSTATE, EPS, WORK, LENW,
C               8             IERFLG)
C                IF (MSTATE .GT. 2) STOP
C                WRITE(*, '(4E12.3)') TOUT, (Y(I), I=1,3)
C                TOUT = 10.E0*TOUT
C                IF (TOUT .LT. 50.E0) GO TO 10
C                END
C
C                SUBROUTINE F (N, T, Y, YDOT)
C                REAL ALFA, T, Y(*), YDOT(*)
C                ALFA = Y(N+1)
C                YDOT(1) = 1.E0 + ALFA*(Y(2) - Y(1)) - Y(1)*Y(3)
C                YDOT(2) = ALFA*(Y(1) - Y(2)) - Y(2)*Y(3)
C                YDOT(3) = 1.E0 - Y(3)*(Y(1) + Y(2))
C                END
C
C  IV.  OTHER COMMUNICATION TO THE USER  ...............................
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
C    B. The number of evaluations of the right hand side can be found
C       in the WORK array in the location determined by:
C            LENW - (N + 50) + 4
C
C  V.  REMARKS  ........................................................
C
C    For other information, see Section IV of the writeup for SDRIV3.
C
C***REFERENCES  C. W. Gear, Numerical Initial Value Problems in
C                 Ordinary Differential Equations, Prentice-Hall, 1971.
C***ROUTINES CALLED  SDRIV3, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   900329  Initial submission to SLATEC.
C***END PROLOGUE  SDRIV1
      EXTERNAL F
      REAL EPS, EWTCOM(1), HMAX, T, TOUT, WORK(*), Y(*)
      INTEGER I, IDLIW, IERFLG, IERROR, IMPL, LENIW, LENW, LENWCM,
     8        LNWCHK, MINT, MITER, ML, MSTATE, MU, MXN, MXORD, MXSTEP,
     8        N, NDE, NROOT, NSTATE, NTASK
      PARAMETER(MXN = 200, IDLIW = 50)
      INTEGER IWORK(IDLIW+MXN)
      CHARACTER INTGR1*8
      PARAMETER(NROOT = 0, IERROR = 2, MINT = 2, MITER = 2, IMPL = 0,
     8          MXORD = 5, MXSTEP = 1000)
      DATA EWTCOM(1) /1.E0/
C***FIRST EXECUTABLE STATEMENT  SDRIV1
      IF (ABS(MSTATE) .EQ. 0 .OR. ABS(MSTATE) .GT. 7) THEN
        WRITE(INTGR1, '(I8)') MSTATE
        IERFLG = 26
        CALL XERMSG('SLATEC', 'SDRIV1',
     8  'Illegal input.  The magnitude of MSTATE, '//INTGR1//
     8  ', is not in the range 1 to 6 .', IERFLG, 1)
        MSTATE = SIGN(7, MSTATE)
        RETURN
      ELSE IF (ABS(MSTATE) .EQ. 7) THEN
        IERFLG = 999
        CALL XERMSG('SLATEC', 'SDRIV1',
     8  'Illegal input.  The magnitude of MSTATE is 7 .', IERFLG, 2)
        RETURN
      END IF
      IF (N .GT. MXN) THEN
        WRITE(INTGR1, '(I8)') N
        IERFLG = 21
        CALL XERMSG('SLATEC', 'SDRIV1',
     8  'Illegal input.  The number of equations, '//INTGR1//
     8  ', is greater than the maximum allowed: 200 .', IERFLG, 1)
        MSTATE = SIGN(7, MSTATE)
        RETURN
      END IF
      IF (MSTATE .GT. 0) THEN
        NSTATE = MSTATE
        NTASK = 1
      ELSE
        NSTATE = - MSTATE
        NTASK = 3
      END IF
      HMAX = 2.E0*ABS(TOUT - T)
      LENIW = N + IDLIW
      LENWCM = LENW - LENIW
      IF (LENWCM .LT. (N*N + 10*N + 250)) THEN
        LNWCHK = N*N + 10*N + 250 + LENIW
        WRITE(INTGR1, '(I8)') LNWCHK
        IERFLG = 32
        CALL XERMSG('SLATEC', 'SDRIV1',
     8  'Insufficient storage allocated for the work array.  '//
     8  'The required storage is at least '//INTGR1//' .', IERFLG, 1)
        MSTATE = SIGN(7, MSTATE)
        RETURN
      END IF
      IF (NSTATE .NE. 1) THEN
        DO 20 I = 1,LENIW
 20       IWORK(I) = WORK(I+LENWCM)
      END IF
      CALL SDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, EWTCOM,
     8             IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8             LENWCM, IWORK, LENIW, F, F, NDE, MXSTEP, F, F,
     8             IERFLG)
      DO 40 I = 1,LENIW
 40     WORK(I+LENWCM) = IWORK(I)
      IF (NSTATE .LE. 4) THEN
        MSTATE = SIGN(NSTATE, MSTATE)
      ELSE IF (NSTATE .EQ. 6) THEN
        MSTATE = SIGN(5, MSTATE)
      ELSE IF (IERFLG .EQ. 11) THEN
        MSTATE = SIGN(6, MSTATE)
      ELSE IF (IERFLG .GT. 11) THEN
        MSTATE = SIGN(7, MSTATE)
      END IF
      RETURN
      END
