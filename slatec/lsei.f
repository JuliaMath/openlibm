*DECK LSEI
      SUBROUTINE LSEI (W, MDW, ME, MA, MG, N, PRGOPT, X, RNORME, RNORML,
     +   MODE, WS, IP)
C***BEGIN PROLOGUE  LSEI
C***PURPOSE  Solve a linearly constrained least squares problem with
C            equality and inequality constraints, and optionally compute
C            a covariance matrix.
C***LIBRARY   SLATEC
C***CATEGORY  K1A2A, D9
C***TYPE      SINGLE PRECISION (LSEI-S, DLSEI-D)
C***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING,
C             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS,
C             QUADRATIC PROGRAMMING
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     Abstract
C
C     This subprogram solves a linearly constrained least squares
C     problem with both equality and inequality constraints, and, if the
C     user requests, obtains a covariance matrix of the solution
C     parameters.
C
C     Suppose there are given matrices E, A and G of respective
C     dimensions ME by N, MA by N and MG by N, and vectors F, B and H of
C     respective lengths ME, MA and MG.  This subroutine solves the
C     linearly constrained least squares problem
C
C                   EX = F, (E ME by N) (equations to be exactly
C                                       satisfied)
C                   AX = B, (A MA by N) (equations to be
C                                       approximately satisfied,
C                                       least squares sense)
C                   GX .GE. H,(G MG by N) (inequality constraints)
C
C     The inequalities GX .GE. H mean that every component of the
C     product GX must be .GE. the corresponding component of H.
C
C     In case the equality constraints cannot be satisfied, a
C     generalized inverse solution residual vector length is obtained
C     for F-EX.  This is the minimal length possible for F-EX.
C
C     Any values ME .GE. 0, MA .GE. 0, or MG .GE. 0 are permitted.  The
C     rank of the matrix E is estimated during the computation.  We call
C     this value KRANKE.  It is an output parameter in IP(1) defined
C     below.  Using a generalized inverse solution of EX=F, a reduced
C     least squares problem with inequality constraints is obtained.
C     The tolerances used in these tests for determining the rank
C     of E and the rank of the reduced least squares problem are
C     given in Sandia Tech. Rept. SAND-78-1290.  They can be
C     modified by the user if new values are provided in
C     the option list of the array PRGOPT(*).
C
C     The user must dimension all arrays appearing in the call list..
C     W(MDW,N+1),PRGOPT(*),X(N),WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2)
C     where K=MAX(MA+MG,N).  This allows for a solution of a range of
C     problems in the given working space.  The dimension of WS(*)
C     given is a necessary overestimate.  Once a particular problem
C     has been run, the output parameter IP(3) gives the actual
C     dimension required for that problem.
C
C     The parameters for LSEI( ) are
C
C     Input..
C
C     W(*,*),MDW,   The array W(*,*) is doubly subscripted with
C     ME,MA,MG,N    first dimensioning parameter equal to MDW.
C                   For this discussion let us call M = ME+MA+MG.  Then
C                   MDW must satisfy MDW .GE. M.  The condition
C                   MDW .LT. M is an error.
C
C                   The array W(*,*) contains the matrices and vectors
C
C                                  (E  F)
C                                  (A  B)
C                                  (G  H)
C
C                   in rows and columns 1,...,M and 1,...,N+1
C                   respectively.
C
C                   The integers ME, MA, and MG are the
C                   respective matrix row dimensions
C                   of E, A and G.  Each matrix has N columns.
C
C     PRGOPT(*)    This real-valued array is the option vector.
C                  If the user is satisfied with the nominal
C                  subprogram features set
C
C                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
C
C                  Otherwise PRGOPT(*) is a linked list consisting of
C                  groups of data of the following form
C
C                  LINK
C                  KEY
C                  DATA SET
C
C                  The parameters LINK and KEY are each one word.
C                  The DATA SET can be comprised of several words.
C                  The number of items depends on the value of KEY.
C                  The value of LINK points to the first
C                  entry of the next group of data within
C                  PRGOPT(*).  The exception is when there are
C                  no more options to change.  In that
C                  case, LINK=1 and the values KEY and DATA SET
C                  are not referenced.  The general layout of
C                  PRGOPT(*) is as follows.
C
C               ...PRGOPT(1) = LINK1 (link to first entry of next group)
C               .  PRGOPT(2) = KEY1 (key to the option change)
C               .  PRGOPT(3) = data value (data value for this change)
C               .       .
C               .       .
C               .       .
C               ...PRGOPT(LINK1)   = LINK2 (link to the first entry of
C               .                       next group)
C               .  PRGOPT(LINK1+1) = KEY2 (key to the option change)
C               .  PRGOPT(LINK1+2) = data value
C               ...     .
C               .       .
C               .       .
C               ...PRGOPT(LINK) = 1 (no more options to change)
C
C                  Values of LINK that are nonpositive are errors.
C                  A value of LINK .GT. NLINK=100000 is also an error.
C                  This helps prevent using invalid but positive
C                  values of LINK that will probably extend
C                  beyond the program limits of PRGOPT(*).
C                  Unrecognized values of KEY are ignored.  The
C                  order of the options is arbitrary and any number
C                  of options can be changed with the following
C                  restriction.  To prevent cycling in the
C                  processing of the option array, a count of the
C                  number of options changed is maintained.
C                  Whenever this count exceeds NOPT=1000, an error
C                  message is printed and the subprogram returns.
C
C                  Options..
C
C                  KEY=1
C                         Compute in W(*,*) the N by N
C                  covariance matrix of the solution variables
C                  as an output parameter.  Nominally the
C                  covariance matrix will not be computed.
C                  (This requires no user input.)
C                  The data set for this option is a single value.
C                  It must be nonzero when the covariance matrix
C                  is desired.  If it is zero, the covariance
C                  matrix is not computed.  When the covariance matrix
C                  is computed, the first dimensioning parameter
C                  of the array W(*,*) must satisfy MDW .GE. MAX(M,N).
C
C                  KEY=10
C                         Suppress scaling of the inverse of the
C                  normal matrix by the scale factor RNORM**2/
C                  MAX(1, no. of degrees of freedom).  This option
C                  only applies when the option for computing the
C                  covariance matrix (KEY=1) is used.  With KEY=1 and
C                  KEY=10 used as options the unscaled inverse of the
C                  normal matrix is returned in W(*,*).
C                  The data set for this option is a single value.
C                  When it is nonzero no scaling is done.  When it is
C                  zero scaling is done.  The nominal case is to do
C                  scaling so if option (KEY=1) is used alone, the
C                  matrix will be scaled on output.
C
C                  KEY=2
C                         Scale the nonzero columns of the
C                         entire data matrix.
C                  (E)
C                  (A)
C                  (G)
C
C                  to have length one.  The data set for this
C                  option is a single value.  It must be
C                  nonzero if unit length column scaling
C                  is desired.
C
C                  KEY=3
C                         Scale columns of the entire data matrix
C                  (E)
C                  (A)
C                  (G)
C
C                  with a user-provided diagonal matrix.
C                  The data set for this option consists
C                  of the N diagonal scaling factors, one for
C                  each matrix column.
C
C                  KEY=4
C                         Change the rank determination tolerance for
C                  the equality constraint equations from
C                  the nominal value of SQRT(SRELPR).  This quantity can
C                  be no smaller than SRELPR, the arithmetic-
C                  storage precision.  The quantity SRELPR is the
C                  largest positive number such that T=1.+SRELPR
C                  satisfies T .EQ. 1.  The quantity used
C                  here is internally restricted to be at
C                  least SRELPR.  The data set for this option
C                  is the new tolerance.
C
C                  KEY=5
C                         Change the rank determination tolerance for
C                  the reduced least squares equations from
C                  the nominal value of SQRT(SRELPR).  This quantity can
C                  be no smaller than SRELPR, the arithmetic-
C                  storage precision.  The quantity used
C                  here is internally restricted to be at
C                  least SRELPR.  The data set for this option
C                  is the new tolerance.
C
C                  For example, suppose we want to change
C                  the tolerance for the reduced least squares
C                  problem, compute the covariance matrix of
C                  the solution parameters, and provide
C                  column scaling for the data matrix.  For
C                  these options the dimension of PRGOPT(*)
C                  must be at least N+9.  The Fortran statements
C                  defining these options would be as follows:
C
C                  PRGOPT(1)=4 (link to entry 4 in PRGOPT(*))
C                  PRGOPT(2)=1 (covariance matrix key)
C                  PRGOPT(3)=1 (covariance matrix wanted)
C
C                  PRGOPT(4)=7 (link to entry 7 in PRGOPT(*))
C                  PRGOPT(5)=5 (least squares equas.  tolerance key)
C                  PRGOPT(6)=... (new value of the tolerance)
C
C                  PRGOPT(7)=N+9 (link to entry N+9 in PRGOPT(*))
C                  PRGOPT(8)=3 (user-provided column scaling key)
C
C                  CALL SCOPY (N, D, 1, PRGOPT(9), 1)  (Copy the N
C                    scaling factors from the user array D(*)
C                    to PRGOPT(9)-PRGOPT(N+8))
C
C                  PRGOPT(N+9)=1 (no more options to change)
C
C                  The contents of PRGOPT(*) are not modified
C                  by the subprogram.
C                  The options for WNNLS( ) can also be included
C                  in this array.  The values of KEY recognized
C                  by WNNLS( ) are 6, 7 and 8.  Their functions
C                  are documented in the usage instructions for
C                  subroutine WNNLS( ).  Normally these options
C                  do not need to be modified when using LSEI( ).
C
C     IP(1),       The amounts of working storage actually
C     IP(2)        allocated for the working arrays WS(*) and
C                  IP(*), respectively.  These quantities are
C                  compared with the actual amounts of storage
C                  needed by LSEI( ).  Insufficient storage
C                  allocated for either WS(*) or IP(*) is an
C                  error.  This feature was included in LSEI( )
C                  because miscalculating the storage formulas
C                  for WS(*) and IP(*) might very well lead to
C                  subtle and hard-to-find execution errors.
C
C                  The length of WS(*) must be at least
C
C                  LW = 2*(ME+N)+K+(MG+2)*(N+7)
C
C                  where K = max(MA+MG,N)
C                  This test will not be made if IP(1).LE.0.
C
C                  The length of IP(*) must be at least
C
C                  LIP = MG+2*N+2
C                  This test will not be made if IP(2).LE.0.
C
C     Output..
C
C     X(*),RNORME,  The array X(*) contains the solution parameters
C     RNORML        if the integer output flag MODE = 0 or 1.
C                   The definition of MODE is given directly below.
C                   When MODE = 0 or 1, RNORME and RNORML
C                   respectively contain the residual vector
C                   Euclidean lengths of F - EX and B - AX.  When
C                   MODE=1 the equality constraint equations EX=F
C                   are contradictory, so RNORME .NE. 0.  The residual
C                   vector F-EX has minimal Euclidean length.  For
C                   MODE .GE. 2, none of these parameters is defined.
C
C     MODE          Integer flag that indicates the subprogram
C                   status after completion.  If MODE .GE. 2, no
C                   solution has been computed.
C
C                   MODE =
C
C                   0  Both equality and inequality constraints
C                      are compatible and have been satisfied.
C
C                   1  Equality constraints are contradictory.
C                      A generalized inverse solution of EX=F was used
C                      to minimize the residual vector length F-EX.
C                      In this sense, the solution is still meaningful.
C
C                   2  Inequality constraints are contradictory.
C
C                   3  Both equality and inequality constraints
C                      are contradictory.
C
C                   The following interpretation of
C                   MODE=1,2 or 3 must be made.  The
C                   sets consisting of all solutions
C                   of the equality constraints EX=F
C                   and all vectors satisfying GX .GE. H
C                   have no points in common.  (In
C                   particular this does not say that
C                   each individual set has no points
C                   at all, although this could be the
C                   case.)
C
C                   4  Usage error occurred.  The value
C                      of MDW is .LT. ME+MA+MG, MDW is
C                      .LT. N and a covariance matrix is
C                      requested, or the option vector
C                      PRGOPT(*) is not properly defined,
C                      or the lengths of the working arrays
C                      WS(*) and IP(*), when specified in
C                      IP(1) and IP(2) respectively, are not
C                      long enough.
C
C     W(*,*)        The array W(*,*) contains the N by N symmetric
C                   covariance matrix of the solution parameters,
C                   provided this was requested on input with
C                   the option vector PRGOPT(*) and the output
C                   flag is returned with MODE = 0 or 1.
C
C     IP(*)         The integer working array has three entries
C                   that provide rank and working array length
C                   information after completion.
C
C                      IP(1) = rank of equality constraint
C                              matrix.  Define this quantity
C                              as KRANKE.
C
C                      IP(2) = rank of reduced least squares
C                              problem.
C
C                      IP(3) = the amount of storage in the
C                              working array WS(*) that was
C                              actually used by the subprogram.
C                              The formula given above for the length
C                              of WS(*) is a necessary overestimate.
C                              If exactly the same problem matrices
C                              are used in subsequent executions,
C                              the declared dimension of WS(*) can
C                              be reduced to this output value.
C     User Designated
C     Working Arrays..
C
C     WS(*),IP(*)              These are respectively type real
C                              and type integer working arrays.
C                              Their required minimal lengths are
C                              given above.
C
C***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for
C                 linear least squares problems with equality and
C                 nonnegativity constraints, Report SAND77-0552, Sandia
C                 Laboratories, June 1978.
C               K. H. Haskell and R. J. Hanson, Selected algorithms for
C                 the linearly constrained least squares problem - a
C                 users guide, Report SAND78-1290, Sandia Laboratories,
C                 August 1979.
C               K. H. Haskell and R. J. Hanson, An algorithm for
C                 linear least squares problems with equality and
C                 nonnegativity constraints, Mathematical Programming
C                 21 (1981), pp. 98-118.
C               R. J. Hanson and K. H. Haskell, Two algorithms for the
C                 linearly constrained least squares problem, ACM
C                 Transactions on Mathematical Software, September 1982.
C***ROUTINES CALLED  H12, LSI, R1MACH, SASUM, SAXPY, SCOPY, SDOT, SNRM2,
C                    SSCAL, SSWAP, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and extensively revised (WRB & RWC)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  LSEI
      INTEGER IP(3), MA, MDW, ME, MG, MODE, N
      REAL             PRGOPT(*), RNORME, RNORML, W(MDW,*), WS(*), X(*)
C
      EXTERNAL H12, LSI, R1MACH, SASUM, SAXPY, SCOPY, SDOT, SNRM2,
     *   SSCAL, SSWAP, XERMSG
      REAL             R1MACH, SASUM, SDOT, SNRM2
C
      REAL             ENORM, FNORM, GAM, RB, RN, RNMAX, SIZE, SN,
     *   SNMAX, SRELPR, T, TAU, UJ, UP, VJ, XNORM, XNRME
      INTEGER I, IMAX, J, JP1, K, KEY, KRANKE, LAST, LCHK, LINK, M,
     *   MAPKE1, MDEQC, MEND, MEP1, N1, N2, NEXT, NLINK, NOPT, NP1,
     *   NTIMES
      LOGICAL COV, FIRST
      CHARACTER*8 XERN1, XERN2, XERN3, XERN4
      SAVE FIRST, SRELPR
C
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  LSEI
C
C     Set the nominal tolerance used in the code for the equality
C     constraint equations.
C
      IF (FIRST) SRELPR = R1MACH(4)
      FIRST = .FALSE.
      TAU = SQRT(SRELPR)
C
C     Check that enough storage was allocated in WS(*) and IP(*).
C
      MODE = 4
      IF (MIN(N,ME,MA,MG) .LT. 0) THEN
         WRITE (XERN1, '(I8)') N
         WRITE (XERN2, '(I8)') ME
         WRITE (XERN3, '(I8)') MA
         WRITE (XERN4, '(I8)') MG
         CALL XERMSG ('SLATEC', 'LSEI', 'ALL OF THE VARIABLES N, ME,' //
     *      ' MA, MG MUST BE .GE. 0$$ENTERED ROUTINE WITH' //
     *      '$$N  = ' // XERN1 //
     *      '$$ME = ' // XERN2 //
     *      '$$MA = ' // XERN3 //
     *      '$$MG = ' // XERN4, 2, 1)
         RETURN
      ENDIF
C
      IF (IP(1).GT.0) THEN
         LCHK = 2*(ME+N) + MAX(MA+MG,N) + (MG+2)*(N+7)
         IF (IP(1).LT.LCHK) THEN
            WRITE (XERN1, '(I8)') LCHK
            CALL XERMSG ('SLATEC', 'LSEI', 'INSUFFICIENT STORAGE ' //
     *         'ALLOCATED FOR WS(*), NEED LW = ' // XERN1, 2, 1)
            RETURN
         ENDIF
      ENDIF
C
      IF (IP(2).GT.0) THEN
         LCHK = MG + 2*N + 2
         IF (IP(2).LT.LCHK) THEN
            WRITE (XERN1, '(I8)') LCHK
            CALL XERMSG ('SLATEC', 'LSEI', 'INSUFFICIENT STORAGE ' //
     *         'ALLOCATED FOR IP(*), NEED LIP = ' // XERN1, 2, 1)
            RETURN
         ENDIF
      ENDIF
C
C     Compute number of possible right multiplying Householder
C     transformations.
C
      M = ME + MA + MG
      IF (N.LE.0 .OR. M.LE.0) THEN
         MODE = 0
         RNORME = 0
         RNORML = 0
         RETURN
      ENDIF
C
      IF (MDW.LT.M) THEN
         CALL XERMSG ('SLATEC', 'LSEI', 'MDW.LT.ME+MA+MG IS AN ERROR',
     +      2, 1)
         RETURN
      ENDIF
C
      NP1 = N + 1
      KRANKE = MIN(ME,N)
      N1 = 2*KRANKE + 1
      N2 = N1 + N
C
C     Set nominal values.
C
C     The nominal column scaling used in the code is
C     the identity scaling.
C
      CALL SCOPY (N, 1.E0, 0, WS(N1), 1)
C
C     No covariance matrix is nominally computed.
C
      COV = .FALSE.
C
C     Process option vector.
C     Define bound for number of options to change.
C
      NOPT = 1000
      NTIMES = 0
C
C     Define bound for positive values of LINK.
C
      NLINK = 100000
      LAST = 1
      LINK = PRGOPT(1)
      IF (LINK.EQ.0 .OR. LINK.GT.NLINK) THEN
         CALL XERMSG ('SLATEC', 'LSEI',
     +      'THE OPTION VECTOR IS UNDEFINED', 2, 1)
         RETURN
      ENDIF
C
  100 IF (LINK.GT.1) THEN
         NTIMES = NTIMES + 1
         IF (NTIMES.GT.NOPT) THEN
            CALL XERMSG ('SLATEC', 'LSEI',
     +         'THE LINKS IN THE OPTION VECTOR ARE CYCLING.', 2, 1)
            RETURN
         ENDIF
C
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.1) THEN
            COV = PRGOPT(LAST+2) .NE. 0.E0
         ELSEIF (KEY.EQ.2 .AND. PRGOPT(LAST+2).NE.0.E0) THEN
            DO 110 J = 1,N
               T = SNRM2(M,W(1,J),1)
               IF (T.NE.0.E0) T = 1.E0/T
               WS(J+N1-1) = T
  110       CONTINUE
         ELSEIF (KEY.EQ.3) THEN
            CALL SCOPY (N, PRGOPT(LAST+2), 1, WS(N1), 1)
         ELSEIF (KEY.EQ.4) THEN
            TAU = MAX(SRELPR,PRGOPT(LAST+2))
         ENDIF
C
         NEXT = PRGOPT(LINK)
         IF (NEXT.LE.0 .OR. NEXT.GT.NLINK) THEN
         CALL XERMSG ('SLATEC', 'LSEI',
     +      'THE OPTION VECTOR IS UNDEFINED', 2, 1)
            RETURN
         ENDIF
C
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
C
      DO 120 J = 1,N
         CALL SSCAL (M, WS(N1+J-1), W(1,J), 1)
  120 CONTINUE
C
      IF (COV .AND. MDW.LT.N) THEN
         CALL XERMSG ('SLATEC', 'LSEI',
     +      'MDW .LT. N WHEN COV MATRIX NEEDED, IS AN ERROR', 2, 1)
         RETURN
      ENDIF
C
C     Problem definition and option vector OK.
C
      MODE = 0
C
C     Compute norm of equality constraint matrix and right side.
C
      ENORM = 0.E0
      DO 130 J = 1,N
         ENORM = MAX(ENORM,SASUM(ME,W(1,J),1))
  130 CONTINUE
C
      FNORM = SASUM(ME,W(1,NP1),1)
      SNMAX = 0.E0
      RNMAX = 0.E0
      DO 150 I = 1,KRANKE
C
C        Compute maximum ratio of vector lengths. Partition is at
C        column I.
C
         DO 140 K = I,ME
            SN = SDOT(N-I+1,W(K,I),MDW,W(K,I),MDW)
            RN = SDOT(I-1,W(K,1),MDW,W(K,1),MDW)
            IF (RN.EQ.0.E0 .AND. SN.GT.SNMAX) THEN
               SNMAX = SN
               IMAX = K
            ELSEIF (K.EQ.I .OR. SN*RNMAX.GT.RN*SNMAX) THEN
               SNMAX = SN
               RNMAX = RN
               IMAX = K
            ENDIF
  140    CONTINUE
C
C        Interchange rows if necessary.
C
         IF (I.NE.IMAX) CALL SSWAP (NP1, W(I,1), MDW, W(IMAX,1), MDW)
         IF (SNMAX.GT.RNMAX*TAU**2) THEN
C
C        Eliminate elements I+1,...,N in row I.
C
            CALL H12 (1, I, I+1, N, W(I,1), MDW, WS(I), W(I+1,1), MDW,
     +                1, M-I)
         ELSE
            KRANKE = I - 1
            GO TO 160
         ENDIF
  150 CONTINUE
C
C     Save diagonal terms of lower trapezoidal matrix.
C
  160 CALL SCOPY (KRANKE, W, MDW+1, WS(KRANKE+1), 1)
C
C     Use Householder transformation from left to achieve
C     KRANKE by KRANKE upper triangular form.
C
      IF (KRANKE.LT.ME) THEN
         DO 170 K = KRANKE,1,-1
C
C           Apply transformation to matrix cols. 1,...,K-1.
C
            CALL H12 (1, K, KRANKE+1, ME, W(1,K), 1, UP, W, 1, MDW, K-1)
C
C           Apply to rt side vector.
C
            CALL H12 (2, K, KRANKE+1, ME, W(1,K), 1, UP, W(1,NP1), 1, 1,
     +                1)
  170    CONTINUE
      ENDIF
C
C     Solve for variables 1,...,KRANKE in new coordinates.
C
      CALL SCOPY (KRANKE, W(1, NP1), 1, X, 1)
      DO 180 I = 1,KRANKE
         X(I) = (X(I)-SDOT(I-1,W(I,1),MDW,X,1))/W(I,I)
  180 CONTINUE
C
C     Compute residuals for reduced problem.
C
      MEP1 = ME + 1
      RNORML = 0.E0
      DO 190 I = MEP1,M
         W(I,NP1) = W(I,NP1) - SDOT(KRANKE,W(I,1),MDW,X,1)
         SN = SDOT(KRANKE,W(I,1),MDW,W(I,1),MDW)
         RN = SDOT(N-KRANKE,W(I,KRANKE+1),MDW,W(I,KRANKE+1),MDW)
         IF (RN.LE.SN*TAU**2 .AND. KRANKE.LT.N)
     *      CALL SCOPY (N-KRANKE, 0.E0, 0, W(I,KRANKE+1), MDW)
  190 CONTINUE
C
C     Compute equality constraint equations residual length.
C
      RNORME = SNRM2(ME-KRANKE,W(KRANKE+1,NP1),1)
C
C     Move reduced problem data upward if KRANKE.LT.ME.
C
      IF (KRANKE.LT.ME) THEN
         DO 200 J = 1,NP1
            CALL SCOPY (M-ME, W(ME+1,J), 1, W(KRANKE+1,J), 1)
  200    CONTINUE
      ENDIF
C
C     Compute solution of reduced problem.
C
      CALL LSI(W(KRANKE+1, KRANKE+1), MDW, MA, MG, N-KRANKE, PRGOPT,
     +         X(KRANKE+1), RNORML, MODE, WS(N2), IP(2))
C
C     Test for consistency of equality constraints.
C
      IF (ME.GT.0) THEN
         MDEQC = 0
         XNRME = SASUM(KRANKE,W(1,NP1),1)
         IF (RNORME.GT.TAU*(ENORM*XNRME+FNORM)) MDEQC = 1
         MODE = MODE + MDEQC
C
C        Check if solution to equality constraints satisfies inequality
C        constraints when there are no degrees of freedom left.
C
         IF (KRANKE.EQ.N .AND. MG.GT.0) THEN
            XNORM = SASUM(N,X,1)
            MAPKE1 = MA + KRANKE + 1
            MEND = MA + KRANKE + MG
            DO 210 I = MAPKE1,MEND
               SIZE = SASUM(N,W(I,1),MDW)*XNORM + ABS(W(I,NP1))
               IF (W(I,NP1).GT.TAU*SIZE) THEN
                  MODE = MODE + 2
                  GO TO 290
               ENDIF
  210       CONTINUE
         ENDIF
      ENDIF
C
C     Replace diagonal terms of lower trapezoidal matrix.
C
      IF (KRANKE.GT.0) THEN
         CALL SCOPY (KRANKE, WS(KRANKE+1), 1, W, MDW+1)
C
C        Reapply transformation to put solution in original coordinates.
C
         DO 220 I = KRANKE,1,-1
            CALL H12 (2, I, I+1, N, W(I,1), MDW, WS(I), X, 1, 1, 1)
  220    CONTINUE
C
C        Compute covariance matrix of equality constrained problem.
C
         IF (COV) THEN
            DO 270 J = MIN(KRANKE,N-1),1,-1
               RB = WS(J)*W(J,J)
               IF (RB.NE.0.E0) RB = 1.E0/RB
               JP1 = J + 1
               DO 230 I = JP1,N
                  W(I,J) = RB*SDOT(N-J,W(I,JP1),MDW,W(J,JP1),MDW)
  230          CONTINUE
C
               GAM = 0.5E0*RB*SDOT(N-J,W(JP1,J),1,W(J,JP1),MDW)
               CALL SAXPY (N-J, GAM, W(J,JP1), MDW, W(JP1,J), 1)
               DO 250 I = JP1,N
                  DO 240 K = I,N
                     W(I,K) = W(I,K) + W(J,I)*W(K,J) + W(I,J)*W(J,K)
                     W(K,I) = W(I,K)
  240             CONTINUE
  250          CONTINUE
               UJ = WS(J)
               VJ = GAM*UJ
               W(J,J) = UJ*VJ + UJ*VJ
               DO 260 I = JP1,N
                  W(J,I) = UJ*W(I,J) + VJ*W(J,I)
  260          CONTINUE
               CALL SCOPY (N-J, W(J, JP1), MDW, W(JP1,J), 1)
  270       CONTINUE
         ENDIF
      ENDIF
C
C     Apply the scaling to the covariance matrix.
C
      IF (COV) THEN
         DO 280 I = 1,N
            CALL SSCAL (N, WS(I+N1-1), W(I,1), MDW)
            CALL SSCAL (N, WS(I+N1-1), W(1,I), 1)
  280    CONTINUE
      ENDIF
C
C     Rescale solution vector.
C
  290 IF (MODE.LE.1) THEN
         DO 300 J = 1,N
            X(J) = X(J)*WS(N1+J-1)
  300    CONTINUE
      ENDIF
C
      IP(1) = KRANKE
      IP(3) = IP(3) + 2*KRANKE + N
      RETURN
      END
