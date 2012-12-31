*DECK DEABM
      SUBROUTINE DEABM (F, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID,
     +   RWORK, LRW, IWORK, LIW, RPAR, IPAR)
C***BEGIN PROLOGUE  DEABM
C***PURPOSE  Solve an initial value problem in ordinary differential
C            equations using an Adams-Bashforth method.
C***LIBRARY   SLATEC (DEPAC)
C***CATEGORY  I1A1B
C***TYPE      SINGLE PRECISION (DEABM-S, DDEABM-D)
C***KEYWORDS  ADAMS-BASHFORTH METHOD, DEPAC, INITIAL VALUE PROBLEMS,
C             ODE, ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR
C***AUTHOR  Shampine, L. F., (SNLA)
C           Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   This is the Adams code in the package of differential equation
C   solvers DEPAC, consisting of the codes DERKF, DEABM, and DEBDF.
C   Design of the package was by L. F. Shampine and H. A. Watts.
C   It is documented in
C        SAND79-2374 , DEPAC - Design of a User Oriented Package of ODE
C                              Solvers.
C   DEABM is a driver for a modification of the code ODE written by
C             L. F. Shampine and M. K. Gordon
C             Sandia Laboratories
C             Albuquerque, New Mexico 87185
C
C **********************************************************************
C **             DEPAC PACKAGE OVERVIEW           **
C **************************************************
C
C        You have a choice of three differential equation solvers from
C        DEPAC.  The following brief descriptions are meant to aid you
C        in choosing the most appropriate code for your problem.
C
C        DERKF is a fifth order Runge-Kutta code.  It is the simplest of
C        the three choices, both algorithmically and in the use of the
C        code.  DERKF is primarily designed to solve non-stiff and mild-
C        ly stiff differential equations when derivative evaluations are
C        not expensive.  It should generally not be used to get high
C        accuracy results nor answers at a great many specific points.
C        Because DERKF has very low overhead costs, it will usually
C        result in the least expensive integration when solving
C        problems requiring a modest amount of accuracy and having
C        equations that are not costly to evaluate.  DERKF attempts to
C        discover when it is not suitable for the task posed.
C
C        DEABM is a variable order (one through twelve) Adams code.
C        Its complexity lies somewhere between that of DERKF and DEBDF.
C        DEABM is primarily designed to solve non-stiff and mildly stiff
C        differential equations when derivative evaluations are
C        expensive, high accuracy results are needed or answers at
C        many specific points are required.  DEABM attempts to discover
C        when it is not suitable for the task posed.
C
C        DEBDF is a variable order (one through five) backward
C        differentiation formula code.  It is the most complicated of
C        the three choices.  DEBDF is primarily designed to solve stiff
C        differential equations at crude to moderate tolerances.
C        If the problem is very stiff at all, DERKF and DEABM will be
C        quite inefficient compared to DEBDF.  However, DEBDF will be
C        inefficient compared to DERKF and DEABM on non-stiff problems
C        because it uses much more storage, has a much larger overhead,
C        and the low order formulas will not give high accuracies
C        efficiently.
C
C        The concept of stiffness cannot be described in a few words.
C        If you do not know the problem to be stiff, try either DERKF
C        or DEABM.  Both of these codes will inform you of stiffness
C        when the cost of solving such problems becomes important.
C
C **********************************************************************
C ** ABSTRACT **
C **************
C
C   Subroutine DEABM uses the Adams-Bashforth-Moulton predictor-
C   corrector formulas of orders one through twelve to integrate a
C   system of NEQ first order ordinary differential equations of the
C   form
C                         DU/DX = F(X,U)
C   when the vector Y(*) of initial values for U(*) at X=T is given. The
C   subroutine integrates from T to TOUT.  It is easy to continue the
C   integration to get results at additional TOUT.  This is the interval
C   mode of operation.  It is also easy for the routine to return with
C   the solution at each intermediate step on the way to TOUT.  This is
C   the intermediate-output mode of operation.
C
C   DEABM uses subprograms DES, STEPS, SINTRP, HSTART, HVNRM, R1MACH and
C   the error handling routine XERMSG.  The only machine dependent
C   parameters to be assigned appear in R1MACH.
C
C **********************************************************************
C ** DESCRIPTION OF THE ARGUMENTS TO DEABM (AN OVERVIEW) **
C *********************************************************
C
C   The parameters are
C
C      F -- This is the name of a subroutine which you provide to
C             define the differential equations.
C
C      NEQ -- This is the number of (first order) differential
C             equations to be integrated.
C
C      T -- This is a value of the independent variable.
C
C      Y(*) -- This array contains the solution components at T.
C
C      TOUT -- This is a point at which a solution is desired.
C
C      INFO(*) -- The basic task of the code is to integrate the
C             differential equations from T to TOUT and return an
C             answer at TOUT.  INFO(*) is an integer array which is used
C             to communicate exactly how you want this task to be
C             carried out.
C
C      RTOL, ATOL -- These quantities represent relative and absolute
C             error tolerances which you provide to indicate how
C             accurately you wish the solution to be computed.  You may
C             choose them to be both scalars or else both vectors.
C
C      IDID -- This scalar quantity is an indicator reporting what
C             the code did.  You must monitor this integer variable to
C             decide what action to take next.
C
C      RWORK(*), LRW -- RWORK(*) is a real work array of length LRW
C             which provides the code with needed storage space.
C
C      IWORK(*), LIW -- IWORK(*) is an integer work array of length LIW
C             which provides the code with needed storage space and an
C             across call flag.
C
C      RPAR, IPAR -- These are real and integer parameter arrays which
C             you can use for communication between your calling
C             program and the F subroutine.
C
C  Quantities which are used as input items are
C             NEQ, T, Y(*), TOUT, INFO(*),
C             RTOL, ATOL, RWORK(1), LRW and LIW.
C
C  Quantities which may be altered by the code are
C             T, Y(*), INFO(1), RTOL, ATOL,
C             IDID, RWORK(*) and IWORK(*).
C
C **********************************************************************
C ** INPUT -- WHAT TO DO ON THE FIRST CALL TO DEABM **
C ****************************************************
C
C   The first call of the code is defined to be the start of each new
C   problem.  Read through the descriptions of all the following items,
C   provide sufficient storage space for designated arrays, set
C   appropriate variables for the initialization of the problem, and
C   give information about how you want the problem to be solved.
C
C
C      F -- Provide a subroutine of the form
C                               F(X,U,UPRIME,RPAR,IPAR)
C             to define the system of first order differential equations
C             which is to be solved.  For the given values of X and the
C             vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must
C             evaluate the NEQ components of the system of differential
C             equations  DU/DX = F(X,U)  and store the derivatives in
C             array UPRIME(*), that is,  UPRIME(I) = * DU(I)/DX *  for
C             equations I=1,...,NEQ.
C
C             Subroutine F must not alter X or U(*).  You must declare
C             the name F in an external statement in your program that
C             calls DEABM.  You must dimension U and UPRIME in F.
C
C             RPAR and IPAR are real and integer parameter arrays which
C             you can use for communication between your calling program
C             and subroutine F.  They are not used or altered by DEABM.
C             If you do not need RPAR or IPAR, ignore these parameters
C             by treating them as dummy arguments.  If you do choose to
C             use them, dimension them in your calling program and in F
C             as arrays of appropriate length.
C
C      NEQ -- Set it to the number of differential equations.
C             (NEQ .GE. 1)
C
C      T -- Set it to the initial point of the integration.
C             You must use a program variable for T because the code
C             changes its value.
C
C      Y(*) -- Set this vector to the initial values of the NEQ solution
C             components at the initial point.  You must dimension Y at
C             least NEQ in your calling program.
C
C      TOUT -- Set it to the first point at which a solution
C             is desired.  You can take TOUT = T, in which case the code
C             will evaluate the derivative of the solution at T and
C             return.  Integration either forward in T  (TOUT .GT.  T)
C             or backward in T  (TOUT .LT.  T)  is permitted.
C
C             The code advances the solution from T to TOUT using
C             step sizes which are automatically selected so as to
C             achieve the desired accuracy.  If you wish, the code will
C             return with the solution and its derivative following
C             each intermediate step (intermediate-output mode) so that
C             you can monitor them, but you still must provide TOUT in
C             accord with the basic aim of the code.
C
C             The first step taken by the code is a critical one
C             because it must reflect how fast the solution changes near
C             the initial point.  The code automatically selects an
C             initial step size which is practically always suitable for
C             the problem.  By using the fact that the code will not
C             step past TOUT in the first step, you could, if necessary,
C             restrict the length of the initial step size.
C
C             For some problems it may not be permissible to integrate
C             past a point TSTOP because a discontinuity occurs there
C             or the solution or its derivative is not defined beyond
C             TSTOP.  When you have declared a TSTOP point (see INFO(4)
C             and RWORK(1)), you have told the code not to integrate
C             past TSTOP.  In this case any TOUT beyond TSTOP is invalid
C             input.
C
C      INFO(*) -- Use the INFO array to give the code more details about
C             how you want your problem solved.  This array should be
C             dimensioned of length 15 to accommodate other members of
C             DEPAC or possible future extensions, though DEABM uses
C             only the first four entries.  You must respond to all of
C             the following items which are arranged as questions.  The
C             simplest use of the code corresponds to answering all
C             questions as YES ,i.e. setting all entries of INFO to 0.
C
C        INFO(1) -- This parameter enables the code to initialize
C             itself.  You must set it to indicate the start of every
C             new problem.
C
C            **** Is this the first call for this problem ...
C                  YES -- Set INFO(1) = 0
C                   NO -- Not applicable here.
C                         See below for continuation calls.  ****
C
C        INFO(2) -- How much accuracy you want of your solution
C               is specified by the error tolerances RTOL and ATOL.
C               The simplest use is to take them both to be scalars.
C               To obtain more flexibility, they can both be vectors.
C               The code must be told your choice.
C
C            **** Are both error tolerances RTOL, ATOL scalars ...
C                  YES -- Set INFO(2) = 0
C                         and input scalars for both RTOL and ATOL
C                   NO -- Set INFO(2) = 1
C                         and input arrays for both RTOL and ATOL ****
C
C        INFO(3) -- The code integrates from T in the direction
C               of TOUT by steps.  If you wish, it will return the
C               computed solution and derivative at the next
C               intermediate step (the intermediate-output mode) or
C               TOUT, whichever comes first.  This is a good way to
C               proceed if you want to see the behavior of the solution.
C               If you must have solutions at a great many specific
C               TOUT points, this code will compute them efficiently.
C
C            **** Do you want the solution only at
C                 TOUT (and not at the next intermediate step) ...
C                  YES -- Set INFO(3) = 0
C                   NO -- Set INFO(3) = 1 ****
C
C        INFO(4) -- To handle solutions at a great many specific
C               values TOUT efficiently, this code may integrate past
C               TOUT and interpolate to obtain the result at TOUT.
C               Sometimes it is not possible to integrate beyond some
C               point TSTOP because the equation changes there or it is
C               not defined past TSTOP.  Then you must tell the code
C               not to go past.
C
C            **** Can the integration be carried out without any
C                 restrictions on the independent variable T ...
C                  YES -- Set INFO(4)=0
C                   NO -- Set INFO(4)=1
C                         and define the stopping point TSTOP by
C                         setting RWORK(1)=TSTOP ****
C
C      RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL)
C             error tolerances to tell the code how accurately you want
C             the solution to be computed.  They must be defined as
C             program variables because the code may change them.  You
C             have two choices --
C                  both RTOL and ATOL are scalars. (INFO(2)=0)
C                  both RTOL and ATOL are vectors. (INFO(2)=1)
C             In either case all components must be non-negative.
C
C             The tolerances are used by the code in a local error test
C             at each step which requires roughly that
C                     ABS(LOCAL ERROR) .LE. RTOL*ABS(Y)+ATOL
C             for each vector component.
C             (More specifically, a Euclidean norm is used to measure
C             the size of vectors, and the error test uses the magnitude
C             of the solution at the beginning of the step.)
C
C             The true (global) error is the difference between the true
C             solution of the initial value problem and the computed
C             approximation.  Practically all present day codes,
C             including this one, control the local error at each step
C             and do not even attempt to control the global error
C             directly.  Roughly speaking, they produce a solution Y(T)
C             which satisfies the differential equations with a
C             residual R(T),    DY(T)/DT = F(T,Y(T)) + R(T)   ,
C             and, almost always, R(T) is bounded by the error
C             tolerances.  Usually, but not always, the true accuracy of
C             the computed Y is comparable to the error tolerances. This
C             code will usually, but not always, deliver a more accurate
C             solution if you reduce the tolerances and integrate again.
C             By comparing two such solutions you can get a fairly
C             reliable idea of the true error in the solution at the
C             bigger tolerances.
C
C             Setting ATOL=0.0 results in a pure relative error test on
C             that component.  Setting RTOL=0.0 results in a pure abso-
C             lute error test on that component.  A mixed test with non-
C             zero RTOL and ATOL corresponds roughly to a relative error
C             test when the solution component is much bigger than ATOL
C             and to an absolute error test when the solution component
C             is smaller than the threshold ATOL.
C
C             Proper selection of the absolute error control parameters
C             ATOL  requires you to have some idea of the scale of the
C             solution components.  To acquire this information may mean
C             that you will have to solve the problem more than once.
C             In the absence of scale information, you should ask for
C             some relative accuracy in all the components (by setting
C             RTOL values non-zero) and perhaps impose extremely small
C             absolute error tolerances to protect against the danger of
C             a solution component becoming zero.
C
C             The code will not attempt to compute a solution at an
C             accuracy unreasonable for the machine being used.  It will
C             advise you if you ask for too much accuracy and inform
C             you as to the maximum accuracy it believes possible.
C
C      RWORK(*) -- Dimension this real work array of length LRW in your
C             calling program.
C
C      RWORK(1) -- If you have set INFO(4)=0, you can ignore this
C             optional input parameter.  Otherwise you must define a
C             stopping point TSTOP by setting   RWORK(1) = TSTOP.
C             (for some problems it may not be permissible to integrate
C             past a point TSTOP because a discontinuity occurs there
C             or the solution or its derivative is not defined beyond
C             TSTOP.)
C
C      LRW -- Set it to the declared length of the RWORK array.
C             You must have  LRW .GE. 130+21*NEQ
C
C      IWORK(*) -- Dimension this integer work array of length LIW in
C             your calling program.
C
C      LIW -- Set it to the declared length of the IWORK array.
C             You must have  LIW .GE. 51
C
C      RPAR, IPAR -- These are parameter arrays, of real and integer
C             type, respectively.  You can use them for communication
C             between your program that calls DEABM and the  F
C             subroutine.  They are not used or altered by DEABM.  If
C             you do not need RPAR or IPAR, ignore these parameters by
C             treating them as dummy arguments.  If you do choose to use
C             them, dimension them in your calling program and in F as
C             arrays of appropriate length.
C
C **********************************************************************
C ** OUTPUT -- AFTER ANY RETURN FROM DEABM **
C *******************************************
C
C   The principal aim of the code is to return a computed solution at
C   TOUT, although it is also possible to obtain intermediate results
C   along the way.  To find out whether the code achieved its goal
C   or if the integration process was interrupted before the task was
C   completed, you must check the IDID parameter.
C
C
C      T -- The solution was successfully advanced to the
C             output value of T.
C
C      Y(*) -- Contains the computed solution approximation at T.
C             You may also be interested in the approximate derivative
C             of the solution at T.  It is contained in
C             RWORK(21),...,RWORK(20+NEQ).
C
C      IDID -- Reports what the code did
C
C                         *** Task Completed ***
C                   reported by positive values of IDID
C
C             IDID = 1 -- A step was successfully taken in the
C                       intermediate-output mode.  The code has not
C                       yet reached TOUT.
C
C             IDID = 2 -- The integration to TOUT was successfully
C                       completed (T=TOUT) by stepping exactly to TOUT.
C
C             IDID = 3 -- The integration to TOUT was successfully
C                       completed (T=TOUT) by stepping past TOUT.
C                       Y(*) is obtained by interpolation.
C
C                         *** Task Interrupted ***
C                   reported by negative values of IDID
C
C             IDID = -1 -- A large amount of work has been expended.
C                       (500 steps attempted)
C
C             IDID = -2 -- The error tolerances are too stringent.
C
C             IDID = -3 -- The local error test cannot be satisfied
C                       because you specified a zero component in ATOL
C                       and the corresponding computed solution
C                       component is zero.  Thus, a pure relative error
C                       test is impossible for this component.
C
C             IDID = -4 -- The problem appears to be stiff.
C
C             IDID = -5,-6,-7,..,-32  -- Not applicable for this code
C                       but used by other members of DEPAC or possible
C                       future extensions.
C
C                         *** Task Terminated ***
C                   reported by the value of IDID=-33
C
C             IDID = -33 -- The code has encountered trouble from which
C                       it cannot recover.  A message is printed
C                       explaining the trouble and control is returned
C                       to the calling program.  For example, this
C                       occurs when invalid input is detected.
C
C      RTOL, ATOL -- These quantities remain unchanged except when
C             IDID = -2.  In this case, the error tolerances have been
C             increased by the code to values which are estimated to be
C             appropriate for continuing the integration.  However, the
C             reported solution at T was obtained using the input values
C             of RTOL and ATOL.
C
C      RWORK, IWORK -- Contain information which is usually of no
C             interest to the user but necessary for subsequent calls.
C             However, you may find use for
C
C             RWORK(11)--Which contains the step size H to be
C                        attempted on the next step.
C
C             RWORK(12)--If the tolerances have been increased by the
C                        code (IDID = -2) , they were multiplied by the
C                        value in RWORK(12).
C
C             RWORK(13)--Which contains the current value of the
C                        independent variable, i.e.  the farthest point
C                        integration has reached.  This will be dif-
C                        ferent from T only when interpolation has been
C                        performed (IDID=3).
C
C             RWORK(20+I)--Which contains the approximate derivative of
C                        the solution component Y(I).  In DEABM, it is
C                        obtained by calling subroutine F to evaluate
C                        the differential equation using T and Y(*) when
C                        IDID=1 or 2, and by interpolation when IDID=3.
C
C **********************************************************************
C ** INPUT -- WHAT TO DO TO CONTINUE THE INTEGRATION **
C **             (CALLS AFTER THE FIRST)             **
C *****************************************************
C
C        This code is organized so that subsequent calls to continue the
C        integration involve little (if any) additional effort on your
C        part.  You must monitor the IDID parameter in order to
C        determine what to do next.
C
C        Recalling that the principal task of the code is to integrate
C        from T to TOUT (the interval mode), usually all you will need
C        to do is specify a new TOUT upon reaching the current TOUT.
C
C        Do not alter any quantity not specifically permitted below,
C        in particular do not alter NEQ, T, Y(*), RWORK(*), IWORK(*) or
C        the differential equation in subroutine F.  Any such alteration
C        constitutes a new problem and must be treated as such, i.e.
C        you must start afresh.
C
C        You cannot change from vector to scalar error control or vice
C        versa (INFO(2)) but you can change the size of the entries of
C        RTOL, ATOL.  Increasing a tolerance makes the equation easier
C        to integrate.  Decreasing a tolerance will make the equation
C        harder to integrate and should generally be avoided.
C
C        You can switch from the intermediate-output mode to the
C        interval mode (INFO(3)) or vice versa at any time.
C
C        If it has been necessary to prevent the integration from going
C        past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
C        code will not integrate to any TOUT beyond the currently
C        specified TSTOP.  Once TSTOP has been reached you must change
C        the value of TSTOP or set INFO(4)=0.  You may change INFO(4)
C        or TSTOP at any time but you must supply the value of TSTOP in
C        RWORK(1) whenever you set INFO(4)=1.
C
C        The parameter INFO(1) is used by the code to indicate the
C        beginning of a new problem and to indicate whether integration
C        is to be continued.  You must input the value  INFO(1) = 0
C        when starting a new problem.  You must input the value
C        INFO(1) = 1  if you wish to continue after an interrupted task.
C        Do not set  INFO(1) = 0  on a continuation call unless you
C        want the code to restart at the current T.
C
C                         *** Following a Completed Task ***
C         If
C             IDID = 1, call the code again to continue the integration
C                     another step in the direction of TOUT.
C
C             IDID = 2 or 3, define a new TOUT and call the code again.
C                     TOUT must be different from T.  You cannot change
C                     the direction of integration without restarting.
C
C                         *** Following an Interrupted Task ***
C                     To show the code that you realize the task was
C                     interrupted and that you want to continue, you
C                     must take appropriate action and reset INFO(1) = 1
C         If
C             IDID = -1, the code has attempted 500 steps.
C                     If you want to continue, set INFO(1) = 1 and
C                     call the code again.  An additional 500 steps
C                     will be allowed.
C
C             IDID = -2, the error tolerances RTOL, ATOL have been
C                     increased to values the code estimates appropriate
C                     for continuing.  You may want to change them
C                     yourself.  If you are sure you want to continue
C                     with relaxed error tolerances, set INFO(1)=1 and
C                     call the code again.
C
C             IDID = -3, a solution component is zero and you set the
C                     corresponding component of ATOL to zero.  If you
C                     are sure you want to continue, you must first
C                     alter the error criterion to use positive values
C                     for those components of ATOL corresponding to zero
C                     solution components, then set INFO(1)=1 and call
C                     the code again.
C
C             IDID = -4, the problem appears to be stiff.  It is very
C                     inefficient to solve such problems with DEABM. The
C                     code DEBDF in DEPAC handles this task efficiently.
C                     If you are absolutely sure you want to continue
C                     with DEABM, set INFO(1)=1 and call the code again.
C
C             IDID = -5,-6,-7,..,-32  --- cannot occur with this code
C                     but used by other members of DEPAC or possible
C                     future extensions.
C
C                         *** Following a Terminated Task ***
C         If
C             IDID = -33, you cannot continue the solution of this
C                     problem.  An attempt to do so will result in your
C                     run being terminated.
C
C **********************************************************************
C
C***REFERENCES  L. F. Shampine and H. A. Watts, DEPAC - design of a user
C                 oriented package of ODE solvers, Report SAND79-2374,
C                 Sandia Laboratories, 1979.
C***ROUTINES CALLED  DES, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800501  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891024  Changed references from VNORM to HVNRM.  (WRB)
C   891024  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DEABM
C
      LOGICAL START,PHASE1,NORND,STIFF,INTOUT
C
      DIMENSION Y(*),INFO(15),RTOL(*),ATOL(*),RWORK(*),IWORK(*),
     1          RPAR(*),IPAR(*)
C
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3
C
      EXTERNAL F
C
C     CHECK FOR AN APPARENT INFINITE LOOP
C
C***FIRST EXECUTABLE STATEMENT  DEABM
      IF ( INFO(1) .EQ. 0 ) IWORK(LIW) = 0
      IF (IWORK(LIW) .GE. 5) THEN
         IF (T .EQ. RWORK(21 + NEQ)) THEN
            WRITE (XERN3, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DEABM',
     *         'AN APPARENT INFINITE LOOP HAS BEEN DETECTED.$$' //
     *         'YOU HAVE MADE REPEATED CALLS AT T = ' // XERN3 //
     *         ' AND THE INTEGRATION HAS NOT ADVANCED.  CHECK THE ' //
     *         'WAY YOU HAVE SET PARAMETERS FOR THE CALL TO THE ' //
     *         'CODE, PARTICULARLY INFO(1).',  13, 2)
            RETURN
         ENDIF
      ENDIF
C
C     CHECK LRW AND LIW FOR SUFFICIENT STORAGE ALLOCATION
C
      IDID=0
      IF (LRW .LT. 130+21*NEQ) THEN
         WRITE (XERN1, '(I8)') LRW
         CALL XERMSG ('SLATEC', 'DEABM', 'THE LENGTH OF THE RWORK ' //
     *      'ARRAY MUST BE AT LEAST 130 + 21*NEQ.$$' //
     *      'YOU HAVE CALLED THE CODE WITH LRW = ' // XERN1, 1, 1)
         IDID=-33
      ENDIF
C
      IF (LIW .LT. 51) THEN
         WRITE (XERN1, '(I8)') LIW
         CALL XERMSG ('SLATEC', 'DEABM', 'THE LENGTH OF THE IWORK ' //
     *      'ARRAY MUST BE AT LEAST 51.$$YOU HAVE CALLED THE CODE ' //
     *      'WITH LIW = ' // XERN1, 2, 1)
         IDID=-33
      ENDIF
C
C     COMPUTE THE INDICES FOR THE ARRAYS TO BE STORED IN THE WORK ARRAY
C
      IYPOUT = 21
      ITSTAR = NEQ + 21
      IYP = 1 + ITSTAR
      IYY = NEQ + IYP
      IWT = NEQ + IYY
      IP = NEQ + IWT
      IPHI = NEQ + IP
      IALPHA = (NEQ*16) + IPHI
      IBETA = 12 + IALPHA
      IPSI = 12 + IBETA
      IV = 12 + IPSI
      IW = 12 + IV
      ISIG = 12 + IW
      IG = 13 + ISIG
      IGI = 13 + IG
      IXOLD = 11 + IGI
      IHOLD = 1 + IXOLD
      ITOLD = 1 + IHOLD
      IDELSN = 1 + ITOLD
      ITWOU = 1 + IDELSN
      IFOURU = 1 + ITWOU
C
      RWORK(ITSTAR) = T
      IF (INFO(1) .EQ. 0) GO TO 50
      START = IWORK(21) .NE. (-1)
      PHASE1 = IWORK(22) .NE. (-1)
      NORND = IWORK(23) .NE. (-1)
      STIFF = IWORK(24) .NE. (-1)
      INTOUT = IWORK(25) .NE. (-1)
C
 50   CALL DES(F,NEQ,T,Y,TOUT,INFO,RTOL,ATOL,IDID,RWORK(IYPOUT),
     1         RWORK(IYP),RWORK(IYY),RWORK(IWT),RWORK(IP),RWORK(IPHI),
     2         RWORK(IALPHA),RWORK(IBETA),RWORK(IPSI),RWORK(IV),
     3         RWORK(IW),RWORK(ISIG),RWORK(IG),RWORK(IGI),RWORK(11),
     4         RWORK(12),RWORK(13),RWORK(IXOLD),RWORK(IHOLD),
     5         RWORK(ITOLD),RWORK(IDELSN),RWORK(1),RWORK(ITWOU),
     5         RWORK(IFOURU),START,PHASE1,NORND,STIFF,INTOUT,IWORK(26),
     6         IWORK(27),IWORK(28),IWORK(29),IWORK(30),IWORK(31),
     7         IWORK(32),IWORK(33),IWORK(34),IWORK(35),IWORK(45),
     8         RPAR,IPAR)
C
      IWORK(21) = -1
      IF (START) IWORK(21) = 1
      IWORK(22) = -1
      IF (PHASE1) IWORK(22) = 1
      IWORK(23) = -1
      IF (NORND) IWORK(23) = 1
      IWORK(24) = -1
      IF (STIFF) IWORK(24) = 1
      IWORK(25) = -1
      IF (INTOUT) IWORK(25) = 1
C
      IF (IDID .NE. (-2)) IWORK(LIW) = IWORK(LIW) + 1
      IF (T .NE. RWORK(ITSTAR)) IWORK(LIW) = 0
C
      RETURN
      END
