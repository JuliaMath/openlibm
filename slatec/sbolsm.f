*DECK SBOLSM
      SUBROUTINE SBOLSM (W, MDW, MINPUT, NCOLS, BL, BU, IND, IOPT, X,
     +   RNORM, MODE, RW, WW, SCL, IBASIS, IBB)
C***BEGIN PROLOGUE  SBOLSM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SBOCLS and SBOLS
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SBOLSM-S, DBOLSM-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C          Solve E*X = F (least squares sense) with bounds on
C            selected X values.
C     The user must have DIMENSION statements of the form:
C
C       DIMENSION W(MDW,NCOLS+1), BL(NCOLS), BU(NCOLS),
C      * X(NCOLS+NX), RW(NCOLS), WW(NCOLS), SCL(NCOLS)
C       INTEGER IND(NCOLS), IOPT(1+NI), IBASIS(NCOLS), IBB(NCOLS)
C
C     (Here NX=number of extra locations required for options 1,...,7;
C     NX=0 for no options; here NI=number of extra locations possibly
C     required for options 1-7; NI=0 for no options; NI=14 if all the
C     options are simultaneously in use.)
C
C    INPUT
C    -----
C
C    --------------------
C    W(MDW,*),MINPUT,NCOLS
C    --------------------
C     The array W(*,*) contains the matrix [E:F] on entry. The matrix
C     [E:F] has MINPUT rows and NCOLS+1 columns. This data is placed in
C     the array W(*,*) with E occupying the first NCOLS columns and the
C     right side vector F in column NCOLS+1. The row dimension, MDW, of
C     the array W(*,*) must satisfy the inequality MDW .ge. MINPUT.
C     Other values of MDW are errors. The values of MINPUT and NCOLS
C     must be positive. Other values are errors.
C
C    ------------------
C    BL(*),BU(*),IND(*)
C    ------------------
C     These arrays contain the information about the bounds that the
C     solution values are to satisfy. The value of IND(J) tells the
C     type of bound and BL(J) and BU(J) give the explicit values for
C     the respective upper and lower bounds.
C
C    1.    For IND(J)=1, require X(J) .ge. BL(J).
C    2.    For IND(J)=2, require X(J) .le. BU(J).
C    3.    For IND(J)=3, require X(J) .ge. BL(J) and
C                                X(J) .le. BU(J).
C    4.    For IND(J)=4, no bounds on X(J) are required.
C     The values of BL(*),BL(*) are modified by the subprogram. Values
C     other than 1,2,3 or 4 for IND(J) are errors. In the case IND(J)=3
C     (upper and lower bounds) the condition BL(J) .gt. BU(J) is an
C     error.
C
C    -------
C    IOPT(*)
C    -------
C     This is the array where the user can specify nonstandard options
C     for SBOLSM. Most of the time this feature can be ignored by
C     setting the input value IOPT(1)=99. Occasionally users may have
C     needs that require use of the following subprogram options. For
C     details about how to use the options see below: IOPT(*) CONTENTS.
C
C     Option Number   Brief Statement of Purpose
C     ----- ------   ----- --------- -- -------
C           1         Move the IOPT(*) processing pointer.
C           2         Change rank determination tolerance.
C           3         Change blow-up factor that determines the
C                     size of variables being dropped from active
C                     status.
C           4         Reset the maximum number of iterations to use
C                     in solving the problem.
C           5         The data matrix is triangularized before the
C                     problem is solved whenever (NCOLS/MINPUT) .lt.
C                     FAC. Change the value of FAC.
C           6         Redefine the weighting matrix used for
C                     linear independence checking.
C           7         Debug output is desired.
C          99         No more options to change.
C
C    ----
C    X(*)
C    ----
C     This array is used to pass data associated with options 1,2,3 and
C     5. Ignore this input parameter if none of these options are used.
C     Otherwise see below: IOPT(*) CONTENTS.
C
C    ----------------
C    IBASIS(*),IBB(*)
C    ----------------
C     These arrays must be initialized by the user. The values
C         IBASIS(J)=J, J=1,...,NCOLS
C         IBB(J)   =1, J=1,...,NCOLS
C     are appropriate except when using nonstandard features.
C
C    ------
C    SCL(*)
C    ------
C     This is the array of scaling factors to use on the columns of the
C     matrix E. These values must be defined by the user. To suppress
C     any column scaling set SCL(J)=1.0, J=1,...,NCOLS.
C
C    OUTPUT
C    ------
C
C    ----------
C    X(*),RNORM
C    ----------
C     The array X(*) contains a solution (if MODE .ge. 0 or .eq. -22)
C     for the constrained least squares problem. The value RNORM is the
C     minimum residual vector length.
C
C    ----
C    MODE
C    ----
C     The sign of mode determines whether the subprogram has completed
C     normally, or encountered an error condition or abnormal status.
C     A value of MODE .ge. 0 signifies that the subprogram has completed
C     normally. The value of MODE (.ge. 0) is the number of variables
C     in an active status: not at a bound nor at the value ZERO, for
C     the case of free variables. A negative value of MODE will be one
C     of the 18 cases -38,-37,...,-22, or -1. Values .lt. -1 correspond
C     to an abnormal completion of the subprogram. To understand the
C     abnormal completion codes see below: ERROR MESSAGES for SBOLSM
C     An approximate solution will be returned to the user only when
C     maximum iterations is reached, MODE=-22.
C
C    -----------
C    RW(*),WW(*)
C    -----------
C     These are working arrays each with NCOLS entries. The array RW(*)
C     contains the working (scaled, nonactive) solution values. The
C     array WW(*) contains the working (scaled, active) gradient vector
C     values.
C
C    ----------------
C    IBASIS(*),IBB(*)
C    ----------------
C     These arrays contain information about the status of the solution
C     when MODE .ge. 0. The indices IBASIS(K), K=1,...,MODE, show the
C     nonactive variables; indices IBASIS(K), K=MODE+1,..., NCOLS are
C     the active variables. The value (IBB(J)-1) is the number of times
C     variable J was reflected from its upper bound. (Normally the user
C     can ignore these parameters.)
C
C    IOPT(*) CONTENTS
C    ------- --------
C     The option array allows a user to modify internal variables in
C     the subprogram without recompiling the source code. A central
C     goal of the initial software design was to do a good job for most
C     people. Thus the use of options will be restricted to a select
C     group of users. The processing of the option array proceeds as
C     follows: a pointer, here called LP, is initially set to the value
C     1. The value is updated as the options are processed.  At the
C     pointer position the option number is extracted and used for
C     locating other information that allows for options to be changed.
C     The portion of the array IOPT(*) that is used for each option is
C     fixed; the user and the subprogram both know how many locations
C     are needed for each option. A great deal of error checking is
C     done by the subprogram on the contents of the option array.
C     Nevertheless it is still possible to give the subprogram optional
C     input that is meaningless. For example, some of the options use
C     the location X(NCOLS+IOFF) for passing data. The user must manage
C     the allocation of these locations when more than one piece of
C     option data is being passed to the subprogram.
C
C   1
C   -
C     Move the processing pointer (either forward or backward) to the
C     location IOPT(LP+1). The processing pointer is moved to location
C     LP+2 of IOPT(*) in case IOPT(LP)=-1.  For example to skip over
C     locations 3,...,NCOLS+2 of IOPT(*),
C
C       IOPT(1)=1
C       IOPT(2)=NCOLS+3
C       (IOPT(I), I=3,...,NCOLS+2 are not defined here.)
C       IOPT(NCOLS+3)=99
C       CALL SBOLSM
C
C     CAUTION: Misuse of this option can yield some very hard-to-find
C     bugs.  Use it with care.
C
C   2
C   -
C     The algorithm that solves the bounded least squares problem
C     iteratively drops columns from the active set. This has the
C     effect of joining a new column vector to the QR factorization of
C     the rectangular matrix consisting of the partially triangularized
C     nonactive columns. After triangularizing this matrix a test is
C     made on the size of the pivot element. The column vector is
C     rejected as dependent if the magnitude of the pivot element is
C     .le. TOL* magnitude of the column in components strictly above
C     the pivot element. Nominally the value of this (rank) tolerance
C     is TOL = SQRT(R1MACH(4)). To change only the value of TOL, for
C     example,
C
C       X(NCOLS+1)=TOL
C       IOPT(1)=2
C       IOPT(2)=1
C       IOPT(3)=99
C       CALL SBOLSM
C
C     Generally, if LP is the processing pointer for IOPT(*),
C
C       X(NCOLS+IOFF)=TOL
C       IOPT(LP)=2
C       IOPT(LP+1)=IOFF
C        .
C       CALL SBOLSM
C
C     The required length of IOPT(*) is increased by 2 if option 2 is
C     used; The required length of X(*) is increased by 1. A value of
C     IOFF .le. 0 is an error. A value of TOL .le. R1MACH(4) gives a
C     warning message; it is not considered an error.
C
C   3
C   -
C     A solution component is left active (not used) if, roughly
C     speaking, it seems too large. Mathematically the new component is
C     left active if the magnitude is .ge.((vector norm of F)/(matrix
C     norm of E))/BLOWUP. Nominally the factor BLOWUP = SQRT(R1MACH(4)).
C     To change only the value of BLOWUP, for example,
C
C       X(NCOLS+2)=BLOWUP
C       IOPT(1)=3
C       IOPT(2)=2
C       IOPT(3)=99
C       CALL SBOLSM
C
C     Generally, if LP is the processing pointer for IOPT(*),
C
C       X(NCOLS+IOFF)=BLOWUP
C       IOPT(LP)=3
C       IOPT(LP+1)=IOFF
C        .
C       CALL SBOLSM
C
C     The required length of IOPT(*) is increased by 2 if option 3 is
C     used; the required length of X(*) is increased by 1. A value of
C     IOFF .le. 0 is an error. A value of BLOWUP .le. 0.0 is an error.
C
C   4
C   -
C     Normally the algorithm for solving the bounded least squares
C     problem requires between NCOLS/3 and NCOLS drop-add steps to
C     converge. (this remark is based on examining a small number of
C     test cases.) The amount of arithmetic for such problems is
C     typically about twice that required for linear least squares if
C     there are no bounds and if plane rotations are used in the
C     solution method. Convergence of the algorithm, while
C     mathematically certain, can be much slower than indicated. To
C     avoid this potential but unlikely event ITMAX drop-add steps are
C     permitted. Nominally ITMAX=5*(MAX(MINPUT,NCOLS)). To change the
C     value of ITMAX, for example,
C
C       IOPT(1)=4
C       IOPT(2)=ITMAX
C       IOPT(3)=99
C       CALL SBOLSM
C
C     Generally, if LP is the processing pointer for IOPT(*),
C
C       IOPT(LP)=4
C       IOPT(LP+1)=ITMAX
C        .
C       CALL SBOLSM
C
C     The value of ITMAX must be .gt. 0. Other values are errors. Use
C     of this option increases the required length of IOPT(*) by 2.
C
C   5
C   -
C     For purposes of increased efficiency the MINPUT by NCOLS+1 data
C     matrix [E:F] is triangularized as a first step whenever MINPUT
C     satisfies FAC*MINPUT .gt. NCOLS. Nominally FAC=0.75. To change the
C     value of FAC,
C
C       X(NCOLS+3)=FAC
C       IOPT(1)=5
C       IOPT(2)=3
C       IOPT(3)=99
C       CALL SBOLSM
C
C     Generally, if LP is the processing pointer for IOPT(*),
C
C       X(NCOLS+IOFF)=FAC
C       IOPT(LP)=5
C       IOPT(LP+1)=IOFF
C        .
C       CALL SBOLSM
C
C     The value of FAC must be nonnegative. Other values are errors.
C     Resetting FAC=0.0 suppresses the initial triangularization step.
C     Use of this option increases the required length of IOPT(*) by 2;
C     The required length of of X(*) is increased by 1.
C
C   6
C   -
C     The norm used in testing the magnitudes of the pivot element
C     compared to the mass of the column above the pivot line can be
C     changed. The type of change that this option allows is to weight
C     the components with an index larger than MVAL by the parameter
C     WT. Normally MVAL=0 and WT=1. To change both the values MVAL and
C     WT, where LP is the processing pointer for IOPT(*),
C
C       X(NCOLS+IOFF)=WT
C       IOPT(LP)=6
C       IOPT(LP+1)=IOFF
C       IOPT(LP+2)=MVAL
C
C     Use of this option increases the required length of IOPT(*) by 3.
C     The length of X(*) is increased by 1. Values of MVAL must be
C     nonnegative and not greater than MINPUT. Other values are errors.
C     The value of WT must be positive. Any other value is an error. If
C     either error condition is present a message will be printed.
C
C   7
C   -
C     Debug output, showing the detailed add-drop steps for the
C     constrained least squares problem, is desired. This option is
C     intended to be used to locate suspected bugs.
C
C   99
C   --
C     There are no more options to change.
C
C     The values for options are 1,...,7,99, and are the only ones
C     permitted. Other values are errors. Options -99,-1,...,-7 mean
C     that the repective options 99,1,...,7 are left at their default
C     values. An example is the option to modify the (rank) tolerance:
C
C       X(NCOLS+1)=TOL
C       IOPT(1)=-2
C       IOPT(2)=1
C       IOPT(3)=99
C
C    Error Messages for SBOLSM
C    ----- -------- --- ---------
C    -22    MORE THAN ITMAX = ... ITERATIONS SOLVING BOUNDED LEAST
C           SQUARES PROBLEM.
C
C    -23    THE OPTION NUMBER = ... IS NOT DEFINED.
C
C    -24    THE OFFSET = ... BEYOND POSTION NCOLS = ... MUST BE POSITIVE
C           FOR OPTION NUMBER 2.
C
C    -25    THE TOLERANCE FOR RANK DETERMINATION = ... IS LESS THAN
C           MACHINE PRECISION = ....
C
C    -26    THE OFFSET = ... BEYOND POSITION NCOLS = ... MUST BE POSTIVE
C           FOR OPTION NUMBER 3.
C
C    -27    THE RECIPROCAL OF THE BLOW-UP FACTOR FOR REJECTING VARIABLES
C           MUST BE POSITIVE. NOW = ....
C
C    -28    THE MAXIMUM NUMBER OF ITERATIONS = ... MUST BE POSITIVE.
C
C    -29    THE OFFSET = ... BEYOND POSITION NCOLS = ... MUST BE POSTIVE
C           FOR OPTION NUMBER 5.
C
C    -30    THE FACTOR (NCOLS/MINPUT) WHERE PRETRIANGULARIZING IS
C           PERFORMED MUST BE NONNEGATIVE. NOW = ....
C
C    -31    THE NUMBER OF ROWS = ... MUST BE POSITIVE.
C
C    -32    THE NUMBER OF COLUMNS = ... MUST BE POSTIVE.
C
C    -33    THE ROW DIMENSION OF W(,) = ... MUST BE .GE. THE NUMBER OF
C           ROWS = ....
C
C    -34    FOR J = ... THE CONSTRAINT INDICATOR MUST BE 1-4.
C
C    -35    FOR J = ... THE LOWER BOUND = ... IS .GT. THE UPPER BOUND =
C           ....
C
C    -36    THE INPUT ORDER OF COLUMNS = ... IS NOT BETWEEN 1 AND NCOLS
C           = ....
C
C    -37    THE BOUND POLARITY FLAG IN COMPONENT J = ... MUST BE
C           POSITIVE. NOW = ....
C
C    -38    THE ROW SEPARATOR TO APPLY WEIGHTING (...) MUST LIE BETWEEN
C           0 AND MINPUT = .... WEIGHT = ... MUST BE POSITIVE.
C
C***SEE ALSO  SBOCLS, SBOLS
C***ROUTINES CALLED  IVOUT, R1MACH, SAXPY, SCOPY, SDOT, SMOUT, SNRM2,
C                    SROT, SROTG, SSWAP, SVOUT, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   821220  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920422  Fixed usage of MINPUT.  (WRB)
C   901009  Editorial changes, code now reads from top to bottom.  (RWC)
C***END PROLOGUE  SBOLSM
C
C     PURPOSE
C     -------
C     THIS IS THE MAIN SUBPROGRAM THAT SOLVES THE BOUNDED
C     LEAST SQUARES PROBLEM.  THE PROBLEM SOLVED HERE IS:
C
C     SOLVE E*X =  F  (LEAST SQUARES SENSE)
C     WITH BOUNDS ON SELECTED X VALUES.
C
C     TO CHANGE THIS SUBPROGRAM FROM SINGLE TO DOUBLE PRECISION BEGIN
C     EDITING AT THE CARD 'C++'.
C     CHANGE THE SUBPROGRAM NAME TO DBOLSM AND THE STRINGS
C     /SAXPY/ TO /DAXPY/, /SCOPY/ TO /DCOPY/,
C     /SDOT/ TO /DDOT/, /SNRM2/ TO /DNRM2/,
C     /SROT/ TO /DROT/, /SROTG/ TO /DROTG/, /R1MACH/ TO /D1MACH/,
C     /SVOUT/ TO /DVOUT/, /SMOUT/ TO /DMOUT/,
C     /SSWAP/ TO /DSWAP/, /E0/ TO /D0/,
C     /REAL            / TO /DOUBLE PRECISION/.
C++
C
      REAL W(MDW,*),BL(*),BU(*)
      REAL X(*),RW(*),WW(*),SCL(*)
      REAL ALPHA,BETA,BOU,COLABV,COLBLO
      REAL CL1,CL2,CL3,ONE,BIG
      REAL FAC,RNORM,SC,SS,T,TOLIND,WT
      REAL TWO,T1,T2,WBIG,WLARGE,WMAG,XNEW
      REAL ZERO,SDOT,SNRM2
      REAL R1MACH,TOLSZE
      INTEGER IBASIS(*),IBB(*),IND(*),IOPT(*)
      LOGICAL FOUND,CONSTR
      CHARACTER*8 XERN1, XERN2
      CHARACTER*16 XERN3, XERN4
C
      PARAMETER (ZERO=0.0E0, ONE=1.0E0, TWO=2.0E0)
C
      INEXT(IDUM) = MIN(IDUM+1,MROWS)
C***FIRST EXECUTABLE STATEMENT  SBOLSM
C
C     Verify that the problem dimensions are defined properly.
C
      IF (MINPUT.LE.0) THEN
          WRITE (XERN1, '(I8)') MINPUT
          CALL XERMSG ('SLATEC', 'SBOLSM', 'THE NUMBER OF ROWS = ' //
     *       XERN1 // ' MUST BE POSITIVE.', 31, 1)
          MODE = -31
          RETURN
      ENDIF
C
      IF (NCOLS.LE.0) THEN
          WRITE (XERN1, '(I8)') NCOLS
          CALL XERMSG ('SLATEC', 'SBOLSM', 'THE NUMBER OF COLUMNS = ' //
     *       XERN1 // ' MUST BE POSITIVE.', 32, 1)
          MODE = -32
          RETURN
      ENDIF
C
      IF (MDW.LT.MINPUT) THEN
          WRITE (XERN1, '(I8)') MDW
          WRITE (XERN2, '(I8)') MINPUT
          CALL XERMSG ('SLATEC', 'SBOLSM',
     *       'THE ROW DIMENSION OF W(,) = ' // XERN1 //
     *       ' MUST BE .GE. THE NUMBER OF ROWS = ' // XERN2, 33, 1)
          MODE = -33
          RETURN
      ENDIF
C
C     Verify that bound information is correct.
C
      DO 10 J = 1,NCOLS
          IF (IND(J).LT.1 .OR. IND(J).GT.4) THEN
              WRITE (XERN1, '(I8)') J
              WRITE (XERN2, '(I8)') IND(J)
              CALL XERMSG ('SLATEC', 'SBOLSM', 'FOR J = ' // XERN1 //
     *           ' THE CONSTRAINT INDICATOR MUST BE 1-4', 34, 1)
              MODE = -34
              RETURN
          ENDIF
   10 CONTINUE
C
      DO 20 J = 1,NCOLS
          IF (IND(J).EQ.3) THEN
              IF (BU(J).LT.BL(J)) THEN
                  WRITE (XERN1, '(I8)') J
                  WRITE (XERN3, '(1PE15.6)') BL(J)
                  WRITE (XERN4, '(1PE15.6)') BU(J)
                  CALL XERMSG ('SLATEC', 'SBOLSM', 'FOR J = ' // XERN1
     *               // ' THE LOWER BOUND = ' // XERN3 //
     *               ' IS .GT. THE UPPER BOUND = ' // XERN4, 35, 1)
                  MODE = -35
                  RETURN
              ENDIF
          ENDIF
   20 CONTINUE
C
C     Check that permutation and polarity arrays have been set.
C
      DO 30 J = 1,NCOLS
          IF (IBASIS(J).LT.1 .OR. IBASIS(J).GT.NCOLS) THEN
              WRITE (XERN1, '(I8)') IBASIS(J)
              WRITE (XERN2, '(I8)') NCOLS
              CALL XERMSG ('SLATEC', 'SBOLSM',
     *           'THE INPUT ORDER OF COLUMNS = ' // XERN1 //
     *           ' IS NOT BETWEEN 1 AND NCOLS = ' // XERN2, 36, 1)
              MODE = -36
              RETURN
          ENDIF
C
          IF (IBB(J).LE.0) THEN
              WRITE (XERN1, '(I8)') J
              WRITE (XERN2, '(I8)') IBB(J)
              CALL XERMSG ('SLATEC', 'SBOLSM',
     *           'THE BOUND POLARITY FLAG IN COMPONENT J = ' // XERN1 //
     *           ' MUST BE POSITIVE.$$NOW = ' // XERN2, 37, 1)
              MODE = -37
              RETURN
          ENDIF
   30 CONTINUE
C
C     Process the option array.
C
      FAC = 0.75E0
      TOLIND = SQRT(R1MACH(4))
      TOLSZE = SQRT(R1MACH(4))
      ITMAX = 5*MAX(MINPUT,NCOLS)
      WT = ONE
      MVAL = 0
      IPRINT = 0
C
C     Changes to some parameters can occur through the option array,
C     IOPT(*).  Process this array looking carefully for input data
C     errors.
C
      LP = 0
      LDS = 0
C
C     Test for no more options.
C
  590 LP = LP + LDS
      IP = IOPT(LP+1)
      JP = ABS(IP)
      IF (IP.EQ.99) THEN
          GO TO 470
      ELSE IF (JP.EQ.99) THEN
          LDS = 1
      ELSE IF (JP.EQ.1) THEN
C
C         Move the IOPT(*) processing pointer.
C
          IF (IP.GT.0) THEN
              LP = IOPT(LP+2) - 1
              LDS = 0
          ELSE
              LDS = 2
          ENDIF
      ELSE IF (JP.EQ.2) THEN
C
C         Change tolerance for rank determination.
C
          IF (IP.GT.0) THEN
              IOFF = IOPT(LP+2)
              IF (IOFF.LE.0) THEN
                  WRITE (XERN1, '(I8)') IOFF
                  WRITE (XERN2, '(I8)') NCOLS
                  CALL XERMSG ('SLATEC', 'SBOLSM', 'THE OFFSET = ' //
     *               XERN1 // ' BEYOND POSITION NCOLS = ' // XERN2 //
     *               ' MUST BE POSITIVE FOR OPTION NUMBER 2.', 24, 1)
                  MODE = -24
                  RETURN
              ENDIF
C
              TOLIND = X(NCOLS+IOFF)
              IF (TOLIND.LT.R1MACH(4)) THEN
                  WRITE (XERN3, '(1PE15.6)') TOLIND
                  WRITE (XERN4, '(1PE15.6)') R1MACH(4)
                  CALL XERMSG ('SLATEC', 'SBOLSM',
     *               'THE TOLERANCE FOR RANK DETERMINATION = ' // XERN3
     *               // ' IS LESS THAN MACHINE PRECISION = ' // XERN4,
     *               25, 0)
                  MODE = -25
              ENDIF
          ENDIF
          LDS = 2
      ELSE IF (JP.EQ.3) THEN
C
C         Change blowup factor for allowing variables to become
C         inactive.
C
          IF (IP.GT.0) THEN
              IOFF = IOPT(LP+2)
              IF (IOFF.LE.0) THEN
                  WRITE (XERN1, '(I8)') IOFF
                  WRITE (XERN2, '(I8)') NCOLS
                  CALL XERMSG ('SLATEC', 'SBOLSM', 'THE OFFSET = ' //
     *               XERN1 // ' BEYOND POSITION NCOLS = ' // XERN2 //
     *               ' MUST BE POSITIVE FOR OPTION NUMBER 3.', 26, 1)
                  MODE = -26
                  RETURN
              ENDIF
C
              TOLSZE = X(NCOLS+IOFF)
              IF (TOLSZE.LE.ZERO) THEN
                  WRITE (XERN3, '(1PE15.6)') TOLSZE
                  CALL XERMSG ('SLATEC', 'SBOLSM', 'THE RECIPROCAL ' //
     *               'OF THE BLOW-UP FACTOR FOR REJECTING VARIABLES ' //
     *               'MUST BE POSITIVE.$$NOW = ' // XERN3, 27, 1)
                  MODE = -27
                  RETURN
              ENDIF
          ENDIF
          LDS = 2
      ELSE IF (JP.EQ.4) THEN
C
C         Change the maximum number of iterations allowed.
C
          IF (IP.GT.0) THEN
              ITMAX = IOPT(LP+2)
              IF (ITMAX.LE.0) THEN
                  WRITE (XERN1, '(I8)') ITMAX
                  CALL XERMSG ('SLATEC', 'SBOLSM',
     *               'THE MAXIMUM NUMBER OF ITERATIONS = ' // XERN1 //
     *               ' MUST BE POSITIVE.', 28, 1)
                  MODE = -28
                  RETURN
              ENDIF
          ENDIF
          LDS = 2
      ELSE IF (JP.EQ.5) THEN
C
C         Change the factor for pretriangularizing the data matrix.
C
          IF (IP.GT.0) THEN
              IOFF = IOPT(LP+2)
              IF (IOFF.LE.0) THEN
                  WRITE (XERN1, '(I8)') IOFF
                  WRITE (XERN2, '(I8)') NCOLS
                  CALL XERMSG ('SLATEC', 'SBOLSM', 'THE OFFSET = ' //
     *               XERN1 // ' BEYOND POSITION NCOLS = ' // XERN2 //
     *               ' MUST BE POSITIVE FOR OPTION NUMBER 5.', 29, 1)
                  MODE = -29
                  RETURN
              ENDIF
C
              FAC = X(NCOLS+IOFF)
              IF (FAC.LT.ZERO) THEN
                  WRITE (XERN3, '(1PE15.6)') FAC
                  CALL XERMSG ('SLATEC', 'SBOLSM',
     *               'THE FACTOR (NCOLS/MINPUT) WHERE PRE-' //
     *               'TRIANGULARIZING IS PERFORMED MUST BE NON-' //
     *               'NEGATIVE.$$NOW = ' // XERN3, 30, 0)
                  MODE = -30
                  RETURN
              ENDIF
          ENDIF
          LDS = 2
      ELSE IF (JP.EQ.6) THEN
C
C         Change the weighting factor (from 1.0) to apply to components
C         numbered .gt. MVAL (initially set to 1.)  This trick is needed
C         for applications of this subprogram to the heavily weighted
C         least squares problem that come from equality constraints.
C
          IF (IP.GT.0) THEN
              IOFF = IOPT(LP+2)
              MVAL = IOPT(LP+3)
              WT = X(NCOLS+IOFF)
          ENDIF
C
          IF (MVAL.LT.0 .OR. MVAL.GT.MINPUT .OR. WT.LE.ZERO) THEN
              WRITE (XERN1, '(I8)') MVAL
              WRITE (XERN2, '(I8)') MINPUT
              WRITE (XERN3, '(1PE15.6)') WT
              CALL XERMSG ('SLATEC', 'SBOLSM',
     *           'THE ROW SEPARATOR TO APPLY WEIGHTING (' // XERN1 //
     *           ') MUST LIE BETWEEN 0 AND MINPUT = ' // XERN2 //
     *           '.$$WEIGHT = ' // XERN3 // ' MUST BE POSITIVE.', 38, 0)
              MODE = -38
              RETURN
          ENDIF
          LDS = 3
      ELSE IF (JP.EQ.7) THEN
C
C         Turn on debug output.
C
          IF (IP.GT.0) IPRINT = 1
          LDS = 2
      ELSE
          WRITE (XERN1, '(I8)') IP
          CALL XERMSG ('SLATEC', 'SBOLSM', 'THE OPTION NUMBER = ' //
     *       XERN1 // ' IS NOT DEFINED.', 23, 1)
          MODE = -23
          RETURN
      ENDIF
      GO TO 590
C
C     Pretriangularize rectangular arrays of certain sizes for
C     increased efficiency.
C
  470 IF (FAC*MINPUT.GT.NCOLS) THEN
          DO 490 J = 1,NCOLS+1
              DO 480 I = MINPUT,J+MVAL+1,-1
                  CALL SROTG(W(I-1,J),W(I,J),SC,SS)
                  W(I,J) = ZERO
                  CALL SROT(NCOLS-J+1,W(I-1,J+1),MDW,W(I,J+1),MDW,SC,SS)
  480         CONTINUE
  490     CONTINUE
          MROWS = NCOLS + MVAL + 1
      ELSE
          MROWS = MINPUT
      ENDIF
C
C     Set the X(*) array to zero so all components are defined.
C
      CALL SCOPY(NCOLS,ZERO,0,X,1)
C
C     The arrays IBASIS(*) and IBB(*) are initialized by the calling
C     program and the column scaling is defined in the calling program.
C     'BIG' is plus infinity on this machine.
C
      BIG = R1MACH(2)
      DO 550 J = 1,NCOLS
          IF (IND(J).EQ.1) THEN
              BU(J) = BIG
          ELSE IF (IND(J).EQ.2) THEN
              BL(J) = -BIG
          ELSE IF (IND(J).EQ.4) THEN
              BL(J) = -BIG
              BU(J) = BIG
          ENDIF
  550 CONTINUE
C
      DO 570 J = 1,NCOLS
          IF ((BL(J).LE.ZERO.AND.ZERO.LE.BU(J).AND.ABS(BU(J)).LT.
     *        ABS(BL(J))) .OR. BU(J).LT.ZERO) THEN
              T = BU(J)
              BU(J) = -BL(J)
              BL(J) = -T
              SCL(J) = -SCL(J)
              DO 560 I = 1,MROWS
                  W(I,J) = -W(I,J)
  560         CONTINUE
          ENDIF
C
C         Indices in set T(=TIGHT) are denoted by negative values
C         of IBASIS(*).
C
          IF (BL(J).GE.ZERO) THEN
              IBASIS(J) = -IBASIS(J)
              T = -BL(J)
              BU(J) = BU(J) + T
              CALL SAXPY(MROWS,T,W(1,J),1,W(1,NCOLS+1),1)
          ENDIF
  570 CONTINUE
C
      NSETB = 0
      ITER = 0
C
      IF (IPRINT.GT.0) THEN
          CALL SMOUT(MROWS,NCOLS+1,MDW,W,'('' PRETRI. INPUT MATRIX'')',
     *               -4)
          CALL SVOUT(NCOLS,BL,'('' LOWER BOUNDS'')',-4)
          CALL SVOUT(NCOLS,BU,'('' UPPER BOUNDS'')',-4)
      ENDIF
C
  580 ITER = ITER + 1
      IF (ITER.GT.ITMAX) THEN
         WRITE (XERN1, '(I8)') ITMAX
         CALL XERMSG ('SLATEC', 'SBOLSM', 'MORE THAN ITMAX = ' // XERN1
     *      // ' ITERATIONS SOLVING BOUNDED LEAST SQUARES PROBLEM.',
     *      22, 1)
         MODE = -22
C
C        Rescale and translate variables.
C
         IGOPR = 1
         GO TO 130
      ENDIF
C
C     Find a variable to become non-active.
C                                                 T
C     Compute (negative) of gradient vector, W = E *(F-E*X).
C
      CALL SCOPY(NCOLS,ZERO,0,WW,1)
      DO 200 J = NSETB+1,NCOLS
          JCOL = ABS(IBASIS(J))
          WW(J) = SDOT(MROWS-NSETB,W(INEXT(NSETB),J),1,
     *            W(INEXT(NSETB),NCOLS+1),1)*ABS(SCL(JCOL))
  200 CONTINUE
C
      IF (IPRINT.GT.0) THEN
          CALL SVOUT(NCOLS,WW,'('' GRADIENT VALUES'')',-4)
          CALL IVOUT(NCOLS,IBASIS,'('' INTERNAL VARIABLE ORDER'')',-4)
          CALL IVOUT(NCOLS,IBB,'('' BOUND POLARITY'')',-4)
      ENDIF
C
C     If active set = number of total rows, quit.
C
  210 IF (NSETB.EQ.MROWS) THEN
          FOUND = .FALSE.
          GO TO 120
      ENDIF
C
C     Choose an extremal component of gradient vector for a candidate
C     to become non-active.
C
      WLARGE = -BIG
      WMAG = -BIG
      DO 220 J = NSETB+1,NCOLS
          T = WW(J)
          IF (T.EQ.BIG) GO TO 220
          ITEMP = IBASIS(J)
          JCOL = ABS(ITEMP)
          T1 = SNRM2(MVAL-NSETB,W(INEXT(NSETB),J),1)
          IF (ITEMP.LT.0) THEN
              IF (MOD(IBB(JCOL),2).EQ.0) T = -T
              IF (T.LT.ZERO) GO TO 220
              IF (MVAL.GT.NSETB) T = T1
              IF (T.GT.WLARGE) THEN
                  WLARGE = T
                  JLARGE = J
              ENDIF
          ELSE
              IF (MVAL.GT.NSETB) T = T1
              IF (ABS(T).GT.WMAG) THEN
                  WMAG = ABS(T)
                  JMAG = J
              ENDIF
          ENDIF
  220 CONTINUE
C
C     Choose magnitude of largest component of gradient for candidate.
C
      JBIG = 0
      WBIG = ZERO
      IF (WLARGE.GT.ZERO) THEN
          JBIG = JLARGE
          WBIG = WLARGE
      ENDIF
C
      IF (WMAG.GE.WBIG) THEN
          JBIG = JMAG
          WBIG = WMAG
      ENDIF
C
      IF (JBIG.EQ.0) THEN
          FOUND = .FALSE.
          IF (IPRINT.GT.0) THEN
              CALL IVOUT(0,I,'('' FOUND NO VARIABLE TO ENTER'')',-4)
          ENDIF
          GO TO 120
      ENDIF
C
C     See if the incoming column is sufficiently independent.  This
C     test is made before an elimination is performed.
C
      IF (IPRINT.GT.0)
     *    CALL IVOUT(1,JBIG,'('' TRY TO BRING IN THIS COL.'')',-4)
C
      IF (MVAL.LE.NSETB) THEN
          CL1 = SNRM2(MVAL,W(1,JBIG),1)
          CL2 = ABS(WT)*SNRM2(NSETB-MVAL,W(INEXT(MVAL),JBIG),1)
          CL3 = ABS(WT)*SNRM2(MROWS-NSETB,W(INEXT(NSETB),JBIG),1)
          CALL SROTG(CL1,CL2,SC,SS)
          COLABV = ABS(CL1)
          COLBLO = CL3
      ELSE
          CL1 = SNRM2(NSETB,W(1,JBIG),1)
          CL2 = SNRM2(MVAL-NSETB,W(INEXT(NSETB),JBIG),1)
          CL3 = ABS(WT)*SNRM2(MROWS-MVAL,W(INEXT(MVAL),JBIG),1)
          COLABV = CL1
          CALL SROTG(CL2,CL3,SC,SS)
          COLBLO = ABS(CL2)
      ENDIF
C
      IF (COLBLO.LE.TOLIND*COLABV) THEN
          WW(JBIG) = BIG
          IF (IPRINT.GT.0)
     *        CALL IVOUT(0,I,'('' VARIABLE IS DEPENDENT, NOT USED.'')',
     *           -4)
          GO TO 210
      ENDIF
C
C     Swap matrix columns NSETB+1 and JBIG, plus pointer information,
C     and gradient values.
C
      NSETB = NSETB + 1
      IF (NSETB.NE.JBIG) THEN
          CALL SSWAP(MROWS,W(1,NSETB),1,W(1,JBIG),1)
          CALL SSWAP(1,WW(NSETB),1,WW(JBIG),1)
          ITEMP = IBASIS(NSETB)
          IBASIS(NSETB) = IBASIS(JBIG)
          IBASIS(JBIG) = ITEMP
      ENDIF
C
C     Eliminate entries below the pivot line in column NSETB.
C
      IF (MROWS.GT.NSETB) THEN
          DO 230 I = MROWS,NSETB+1,-1
              IF (I.EQ.MVAL+1) GO TO 230
              CALL SROTG(W(I-1,NSETB),W(I,NSETB),SC,SS)
              W(I,NSETB) = ZERO
              CALL SROT(NCOLS-NSETB+1,W(I-1,NSETB+1),MDW,W(I,NSETB+1),
     *                  MDW,SC,SS)
  230     CONTINUE
C
          IF (MVAL.GE.NSETB .AND. MVAL.LT.MROWS) THEN
              CALL SROTG(W(NSETB,NSETB),W(MVAL+1,NSETB),SC,SS)
              W(MVAL+1,NSETB) = ZERO
              CALL SROT(NCOLS-NSETB+1,W(NSETB,NSETB+1),MDW,
     *                  W(MVAL+1,NSETB+1),MDW,SC,SS)
          ENDIF
      ENDIF
C
      IF (W(NSETB,NSETB).EQ.ZERO) THEN
          WW(NSETB) = BIG
          NSETB = NSETB - 1
          IF (IPRINT.GT.0) THEN
              CALL IVOUT(0,I,'('' PIVOT IS ZERO, NOT USED.'')',-4)
          ENDIF
          GO TO 210
      ENDIF
C
C     Check that new variable is moving in the right direction.
C
      ITEMP = IBASIS(NSETB)
      JCOL = ABS(ITEMP)
      XNEW = (W(NSETB,NCOLS+1)/W(NSETB,NSETB))/ABS(SCL(JCOL))
      IF (ITEMP.LT.0) THEN
C
C         IF(WW(NSETB).GE.ZERO.AND.XNEW.LE.ZERO) exit(quit)
C         IF(WW(NSETB).LE.ZERO.AND.XNEW.GE.ZERO) exit(quit)
C
          IF ((WW(NSETB).GE.ZERO.AND.XNEW.LE.ZERO) .OR.
     *        (WW(NSETB).LE.ZERO.AND.XNEW.GE.ZERO)) GO TO 240
      ENDIF
      FOUND = .TRUE.
      GO TO 120
C
  240 WW(NSETB) = BIG
      NSETB = NSETB - 1
      IF (IPRINT.GT.0)
     *    CALL IVOUT(0,I,'('' VARIABLE HAS BAD DIRECTION, NOT USED.'')',
     *       -4)
      GO TO 210
C
C     Solve the triangular system.
C
  270 CALL SCOPY(NSETB,W(1,NCOLS+1),1,RW,1)
      DO 280 J = NSETB,1,-1
          RW(J) = RW(J)/W(J,J)
          JCOL = ABS(IBASIS(J))
          T = RW(J)
          IF (MOD(IBB(JCOL),2).EQ.0) RW(J) = -RW(J)
          CALL SAXPY(J-1,-T,W(1,J),1,RW,1)
          RW(J) = RW(J)/ABS(SCL(JCOL))
  280 CONTINUE
C
      IF (IPRINT.GT.0) THEN
          CALL SVOUT(NSETB,RW,'('' SOLN. VALUES'')',-4)
          CALL IVOUT(NSETB,IBASIS,'('' COLS. USED'')',-4)
      ENDIF
C
      IF (LGOPR.EQ.2) THEN
          CALL SCOPY(NSETB,RW,1,X,1)
          DO 450 J = 1,NSETB
              ITEMP = IBASIS(J)
              JCOL = ABS(ITEMP)
              IF (ITEMP.LT.0) THEN
                  BOU = ZERO
              ELSE
                  BOU = BL(JCOL)
              ENDIF
C
              IF ((-BOU).NE.BIG) BOU = BOU/ABS(SCL(JCOL))
              IF (X(J).LE.BOU) THEN
                  JDROP1 = J
                  GO TO 340
              ENDIF
C
              BOU = BU(JCOL)
              IF (BOU.NE.BIG) BOU = BOU/ABS(SCL(JCOL))
              IF (X(J).GE.BOU) THEN
                  JDROP2 = J
                  GO TO 340
              ENDIF
  450     CONTINUE
          GO TO 340
      ENDIF
C
C     See if the unconstrained solution (obtained by solving the
C     triangular system) satisfies the problem bounds.
C
      ALPHA = TWO
      BETA = TWO
      X(NSETB) = ZERO
      DO 310 J = 1,NSETB
          ITEMP = IBASIS(J)
          JCOL = ABS(ITEMP)
          T1 = TWO
          T2 = TWO
          IF (ITEMP.LT.0) THEN
              BOU = ZERO
          ELSE
              BOU = BL(JCOL)
          ENDIF
          IF ((-BOU).NE.BIG) BOU = BOU/ABS(SCL(JCOL))
          IF (RW(J).LE.BOU) T1 = (X(J)-BOU)/ (X(J)-RW(J))
          BOU = BU(JCOL)
          IF (BOU.NE.BIG) BOU = BOU/ABS(SCL(JCOL))
          IF (RW(J).GE.BOU) T2 = (BOU-X(J))/ (RW(J)-X(J))
C
C     If not, then compute a step length so that the variables remain
C     feasible.
C
          IF (T1.LT.ALPHA) THEN
              ALPHA = T1
              JDROP1 = J
          ENDIF
C
          IF (T2.LT.BETA) THEN
              BETA = T2
              JDROP2 = J
          ENDIF
  310 CONTINUE
C
      CONSTR = ALPHA .LT. TWO .OR. BETA .LT. TWO
      IF (.NOT.CONSTR) THEN
C
C         Accept the candidate because it satisfies the stated bounds
C         on the variables.
C
          CALL SCOPY(NSETB,RW,1,X,1)
          GO TO 580
      ENDIF
C
C     Take a step that is as large as possible with all variables
C     remaining feasible.
C
      DO 330 J = 1,NSETB
          X(J) = X(J) + MIN(ALPHA,BETA)* (RW(J)-X(J))
  330 CONTINUE
C
      IF (ALPHA.LE.BETA) THEN
          JDROP2 = 0
      ELSE
          JDROP1 = 0
      ENDIF
C
  340 IF (JDROP1+JDROP2.LE.0 .OR. NSETB.LE.0) GO TO 580
  350 JDROP = JDROP1 + JDROP2
      ITEMP = IBASIS(JDROP)
      JCOL = ABS(ITEMP)
      IF (JDROP2.GT.0) THEN
C
C         Variable is at an upper bound.  Subtract multiple of this
C         column from right hand side.
C
          T = BU(JCOL)
          IF (ITEMP.GT.0) THEN
              BU(JCOL) = T - BL(JCOL)
              BL(JCOL) = -T
              ITEMP = -ITEMP
              SCL(JCOL) = -SCL(JCOL)
              DO 360 I = 1,JDROP
                  W(I,JDROP) = -W(I,JDROP)
  360         CONTINUE
          ELSE
              IBB(JCOL) = IBB(JCOL) + 1
              IF (MOD(IBB(JCOL),2).EQ.0) T = -T
          ENDIF
C
C     Variable is at a lower bound.
C
      ELSE
          IF (ITEMP.LT.ZERO) THEN
              T = ZERO
          ELSE
              T = -BL(JCOL)
              BU(JCOL) = BU(JCOL) + T
              ITEMP = -ITEMP
          ENDIF
      ENDIF
C
      CALL SAXPY(JDROP,T,W(1,JDROP),1,W(1,NCOLS+1),1)
C
C     Move certain columns left to achieve upper Hessenberg form.
C
      CALL SCOPY(JDROP,W(1,JDROP),1,RW,1)
      DO 370 J = JDROP+1,NSETB
          IBASIS(J-1) = IBASIS(J)
          X(J-1) = X(J)
          CALL SCOPY(J,W(1,J),1,W(1,J-1),1)
  370 CONTINUE
C
      IBASIS(NSETB) = ITEMP
      W(1,NSETB) = ZERO
      CALL SCOPY(MROWS-JDROP,W(1,NSETB),0,W(JDROP+1,NSETB),1)
      CALL SCOPY(JDROP,RW,1,W(1,NSETB),1)
C
C     Transform the matrix from upper Hessenberg form to upper
C     triangular form.
C
      NSETB = NSETB - 1
      DO 390 I = JDROP,NSETB
C
C         Look for small pivots and avoid mixing weighted and
C         nonweighted rows.
C
          IF (I.EQ.MVAL) THEN
              T = ZERO
              DO 380 J = I,NSETB
                  JCOL = ABS(IBASIS(J))
                  T1 = ABS(W(I,J)*SCL(JCOL))
                  IF (T1.GT.T) THEN
                      JBIG = J
                      T = T1
                  ENDIF
  380         CONTINUE
              GO TO 400
          ENDIF
          CALL SROTG(W(I,I),W(I+1,I),SC,SS)
          W(I+1,I) = ZERO
          CALL SROT(NCOLS-I+1,W(I,I+1),MDW,W(I+1,I+1),MDW,SC,SS)
  390 CONTINUE
      GO TO 430
C
C     The triangularization is completed by giving up the Hessenberg
C     form and triangularizing a rectangular matrix.
C
  400 CALL SSWAP(MROWS,W(1,I),1,W(1,JBIG),1)
      CALL SSWAP(1,WW(I),1,WW(JBIG),1)
      CALL SSWAP(1,X(I),1,X(JBIG),1)
      ITEMP = IBASIS(I)
      IBASIS(I) = IBASIS(JBIG)
      IBASIS(JBIG) = ITEMP
      JBIG = I
      DO 420 J = JBIG,NSETB
          DO 410 I = J+1,MROWS
              CALL SROTG(W(J,J),W(I,J),SC,SS)
              W(I,J) = ZERO
              CALL SROT(NCOLS-J+1,W(J,J+1),MDW,W(I,J+1),MDW,SC,SS)
  410     CONTINUE
  420 CONTINUE
C
C     See if the remaining coefficients are feasible.  They should be
C     because of the way MIN(ALPHA,BETA) was chosen.  Any that are not
C     feasible will be set to their bounds and appropriately translated.
C
  430 JDROP1 = 0
      JDROP2 = 0
      LGOPR = 2
      GO TO 270
C
C     Find a variable to become non-active.
C
  120 IF (FOUND) THEN
          LGOPR = 1
          GO TO 270
      ENDIF
C
C     Rescale and translate variables.
C
      IGOPR = 2
  130 CALL SCOPY(NSETB,X,1,RW,1)
      CALL SCOPY(NCOLS,ZERO,0,X,1)
      DO 140 J = 1,NSETB
          JCOL = ABS(IBASIS(J))
          X(JCOL) = RW(J)*ABS(SCL(JCOL))
  140 CONTINUE
C
      DO 150 J = 1,NCOLS
          IF (MOD(IBB(J),2).EQ.0) X(J) = BU(J) - X(J)
  150 CONTINUE
C
      DO 160 J = 1,NCOLS
          JCOL = IBASIS(J)
          IF (JCOL.LT.0) X(-JCOL) = BL(-JCOL) + X(-JCOL)
  160 CONTINUE
C
      DO 170 J = 1,NCOLS
          IF (SCL(J).LT.ZERO) X(J) = -X(J)
  170 CONTINUE
C
      I = MAX(NSETB,MVAL)
      RNORM = SNRM2(MROWS-I,W(INEXT(I),NCOLS+1),1)
C
      IF (IGOPR.EQ.2) MODE = NSETB
      RETURN
      END
