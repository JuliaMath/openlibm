*DECK DBOLS
      SUBROUTINE DBOLS (W, MDW, MROWS, NCOLS, BL, BU, IND, IOPT, X,
     +   RNORM, MODE, RW, IW)
C***BEGIN PROLOGUE  DBOLS
C***PURPOSE  Solve the problem
C                 E*X = F (in the least  squares  sense)
C            with bounds on selected X values.
C***LIBRARY   SLATEC
C***CATEGORY  K1A2A, G2E, G2H1, G2H2
C***TYPE      DOUBLE PRECISION (SBOLS-S, DBOLS-D)
C***KEYWORDS  BOUNDS, CONSTRAINTS, INEQUALITY, LEAST SQUARES, LINEAR
C***AUTHOR  Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C   **** All INPUT and OUTPUT real variables are DOUBLE PRECISION ****
C
C     The user must have dimension statements of the form:
C
C       DIMENSION W(MDW,NCOLS+1), BL(NCOLS), BU(NCOLS),
C      * X(NCOLS+NX), RW(5*NCOLS)
C       INTEGER IND(NCOLS), IOPT(1+NI), IW(2*NCOLS)
C
C     (Here NX=number of extra locations required for option 4; NX=0
C     for no options; NX=NCOLS if this option is in use. Here NI=number
C     of extra locations required for options 1-6; NI=0 for no
C     options.)
C
C   INPUT
C   -----
C
C    --------------------
C    W(MDW,*),MROWS,NCOLS
C    --------------------
C     The array W(*,*) contains the matrix [E:F] on entry. The matrix
C     [E:F] has MROWS rows and NCOLS+1 columns. This data is placed in
C     the array W(*,*) with E occupying the first NCOLS columns and the
C     right side vector F in column NCOLS+1. The row dimension, MDW, of
C     the array W(*,*) must satisfy the inequality MDW .ge. MROWS.
C     Other values of MDW are errors. The values of MROWS and NCOLS
C     must be positive. Other values are errors. There is an exception
C     to this when using option 1 for accumulation of blocks of
C     equations. In that case MROWS is an OUTPUT variable ONLY, and the
C     matrix data for [E:F] is placed in W(*,*), one block of rows at a
C     time.  MROWS contains the number of rows in the matrix after
C     triangularizing several blocks of equations. This is an OUTPUT
C     parameter ONLY when option 1 is used. See IOPT(*) CONTENTS
C     for details about option 1.
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
C          (the value of BU(J) is not used.)
C    2.    For IND(J)=2, require X(J) .le. BU(J).
C          (the value of BL(J) is not used.)
C    3.    For IND(J)=3, require X(J) .ge. BL(J) and
C                                X(J) .le. BU(J).
C    4.    For IND(J)=4, no bounds on X(J) are required.
C          (the values of BL(J) and BU(J) are not used.)
C
C     Values other than 1,2,3 or 4 for IND(J) are errors. In the case
C     IND(J)=3 (upper and lower bounds) the condition BL(J) .gt. BU(J)
C     is an error.
C
C    -------
C    IOPT(*)
C    -------
C     This is the array where the user can specify nonstandard options
C     for DBOLSM( ). Most of the time this feature can be ignored by
C     setting the input value IOPT(1)=99. Occasionally users may have
C     needs that require use of the following subprogram options. For
C     details about how to use the options see below: IOPT(*) CONTENTS.
C
C     Option Number   Brief Statement of Purpose
C     ------ ------   ----- --------- -- -------
C           1         Return to user for accumulation of blocks
C                     of least squares equations.
C           2         Check lengths of all arrays used in the
C                     subprogram.
C           3         Standard scaling of the data matrix, E.
C           4         User provides column scaling for matrix E.
C           5         Provide option array to the low-level
C                     subprogram DBOLSM( ).
C           6         Move the IOPT(*) processing pointer.
C          99         No more options to change.
C
C    ----
C    X(*)
C    ----
C     This array is used to pass data associated with option 4. Ignore
C     this parameter if this option is not used. Otherwise see below:
C     IOPT(*) CONTENTS.
C
C    OUTPUT
C    ------
C
C    ----------
C    X(*),RNORM
C    ----------
C     The array X(*) contains a solution (if MODE .ge.0 or .eq.-22) for
C     the constrained least squares problem. The value RNORM is the
C     minimum residual vector length.
C
C    ----
C    MODE
C    ----
C     The sign of MODE determines whether the subprogram has completed
C     normally, or encountered an error condition or abnormal status. A
C     value of MODE .ge. 0 signifies that the subprogram has completed
C     normally. The value of MODE (.GE. 0) is the number of variables
C     in an active status: not at a bound nor at the value ZERO, for
C     the case of free variables. A negative value of MODE will be one
C     of the cases -37,-36,...,-22, or -17,...,-2. Values .lt. -1
C     correspond to an abnormal completion of the subprogram. To
C     understand the abnormal completion codes see below: ERROR
C     MESSAGES for DBOLS( ). AN approximate solution will be returned
C     to the user only when max. iterations is reached, MODE=-22.
C     Values for MODE=-37,...,-22 come from the low-level subprogram
C     DBOLSM(). See the section ERROR MESSAGES for DBOLSM() in the
C     documentation for DBOLSM().
C
C    -----------
C    RW(*),IW(*)
C    -----------
C     These are working arrays with 5*NCOLS and 2*NCOLS entries.
C     (normally the user can ignore the contents of these arrays,
C     but they must be dimensioned properly.)
C
C    IOPT(*) CONTENTS
C    ------- --------
C     The option array allows a user to modify internal variables in
C     the subprogram without recompiling the source code. A central
C     goal of the initial software design was to do a good job for most
C     people. Thus the use of options will be restricted to a select
C     group of users. The processing of the option array proceeds as
C     follows: a pointer, here called LP, is initially set to the value
C     1. This value is updated as each option is processed. At the
C     pointer position the option number is extracted and used for
C     locating other information that allows for options to be changed.
C     The portion of the array IOPT(*) that is used for each option is
C     fixed; the user and the subprogram both know how many locations
C     are needed for each option. A great deal of error checking is
C     done by the subprogram on the contents of the option array.
C     Nevertheless it is still possible to give the subprogram optional
C     input that is meaningless. For example option 4 uses the
C     locations X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1) for passing
C     scaling data. The user must manage the allocation of these
C     locations.
C
C   1
C   -
C     This option allows the user to solve problems with a large number
C     of rows compared to the number of variables. The idea is that the
C     subprogram returns to the user (perhaps many times) and receives
C     new least squares equations from the calling program unit.
C     Eventually the user signals "that's all" and then computes the
C     solution with one final call to subprogram DBOLS( ). The value of
C     MROWS is an OUTPUT variable when this option is used. Its value
C     is always in the range 0 .le. MROWS .le. NCOLS+1. It is equal to
C     the number of rows after the triangularization of the entire set
C     of equations. If LP is the processing pointer for IOPT(*), the
C     usage for the sequential processing of blocks of equations is
C
C        IOPT(LP)=1
C        Move block of equations to W(*,*) starting at
C        the first row of W(*,*).
C        IOPT(LP+3)=# of rows in the block; user defined
C
C     The user now calls DBOLS( ) in a loop. The value of IOPT(LP+1)
C     directs the user's action. The value of IOPT(LP+2) points to
C     where the subsequent rows are to be placed in W(*,*).
C
C      .<LOOP
C      . CALL DBOLS()
C      . IF(IOPT(LP+1) .EQ. 1) THEN
C      .    IOPT(LP+3)=# OF ROWS IN THE NEW BLOCK; USER DEFINED
C      .    PLACE NEW BLOCK OF IOPT(LP+3) ROWS IN
C      .    W(*,*) STARTING AT ROW IOPT(LP+2).
C      .
C      .    IF( THIS IS THE LAST BLOCK OF EQUATIONS ) THEN
C      .       IOPT(LP+1)=2
C      .<------CYCLE LOOP
C      .    ELSE IF (IOPT(LP+1) .EQ. 2) THEN
C      <-------EXIT LOOP SOLUTION COMPUTED IF MODE .GE. 0
C      . ELSE
C      . ERROR CONDITION; SHOULD NOT HAPPEN.
C      .<END LOOP
C
C     Use of this option adds 4 to the required length of IOPT(*).
C
C
C   2
C   -
C     This option is useful for checking the lengths of all arrays used
C     by DBOLS() against their actual requirements for this problem.
C     The idea is simple: the user's program unit passes the declared
C     dimension information of the arrays. These values are compared
C     against the problem-dependent needs within the subprogram. If any
C     of the dimensions are too small an error message is printed and a
C     negative value of MODE is returned, -11 to -17. The printed error
C     message tells how long the dimension should be. If LP is the
C     processing pointer for IOPT(*),
C
C        IOPT(LP)=2
C        IOPT(LP+1)=Row dimension of W(*,*)
C        IOPT(LP+2)=Col. dimension of W(*,*)
C        IOPT(LP+3)=Dimensions of BL(*),BU(*),IND(*)
C        IOPT(LP+4)=Dimension of X(*)
C        IOPT(LP+5)=Dimension of RW(*)
C        IOPT(LP+6)=Dimension of IW(*)
C        IOPT(LP+7)=Dimension of IOPT(*)
C         .
C        CALL DBOLS()
C
C     Use of this option adds 8 to the required length of IOPT(*).
C
C   3
C   -
C     This option changes the type of scaling for the data matrix E.
C     Nominally each nonzero column of E is scaled so that the
C     magnitude of its largest entry is equal to the value ONE. If LP
C     is the processing pointer for IOPT(*),
C
C        IOPT(LP)=3
C        IOPT(LP+1)=1,2 or 3
C            1= Nominal scaling as noted;
C            2= Each nonzero column scaled to have length ONE;
C            3= Identity scaling; scaling effectively suppressed.
C         .
C        CALL DBOLS()
C
C     Use of this option adds 2 to the required length of IOPT(*).
C
C   4
C   -
C     This option allows the user to provide arbitrary (positive)
C     column scaling for the matrix E. If LP is the processing pointer
C     for IOPT(*),
C
C        IOPT(LP)=4
C        IOPT(LP+1)=IOFF
C        X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1)
C        = Positive scale factors for cols. of E.
C         .
C        CALL DBOLS()
C
C     Use of this option adds 2 to the required length of IOPT(*) and
C     NCOLS to the required length of X(*).
C
C   5
C   -
C     This option allows the user to provide an option array to the
C     low-level subprogram DBOLSM(). If LP is the processing pointer
C     for IOPT(*),
C
C        IOPT(LP)=5
C        IOPT(LP+1)= Position in IOPT(*) where option array
C                    data for DBOLSM() begins.
C         .
C        CALL DBOLS()
C
C     Use of this option adds 2 to the required length of IOPT(*).
C
C   6
C   -
C     Move the processing pointer (either forward or backward) to the
C     location IOPT(LP+1). The processing point is moved to entry
C     LP+2 of IOPT(*) if the option is left with -6 in IOPT(LP).  For
C     example to skip over locations 3,...,NCOLS+2 of IOPT(*),
C
C       IOPT(1)=6
C       IOPT(2)=NCOLS+3
C       (IOPT(I), I=3,...,NCOLS+2 are not defined here.)
C       IOPT(NCOLS+3)=99
C       CALL DBOLS()
C
C     CAUTION: Misuse of this option can yield some very hard
C     -to-find bugs.  Use it with care.
C
C   99
C   --
C     There are no more options to change.
C
C     Only option numbers -99, -6,-5,...,-1, 1,2,...,6, and 99 are
C     permitted. Other values are errors. Options -99,-1,...,-6 mean
C     that the respective options 99,1,...,6 are left at their default
C     values. An example is the option to modify the (rank) tolerance:
C
C       IOPT(1)=-3 Option is recognized but not changed
C       IOPT(2)=2  Scale nonzero cols. to have length ONE
C       IOPT(3)=99
C
C    ERROR MESSAGES for DBOLS()
C    ----- -------- --- -------
C
C WARNING IN...
C DBOLS(). MDW=(I1) MUST BE POSITIVE.
C           IN ABOVE MESSAGE, I1=         0
C ERROR NUMBER =         2
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). NCOLS=(I1) THE NO. OF VARIABLES MUST BE POSITIVE.
C           IN ABOVE MESSAGE, I1=         0
C ERROR NUMBER =         3
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). FOR J=(I1), IND(J)=(I2) MUST BE 1-4.
C           IN ABOVE MESSAGE, I1=         1
C           IN ABOVE MESSAGE, I2=         0
C ERROR NUMBER =         4
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). FOR J=(I1), BOUND BL(J)=(R1) IS .GT. BU(J)=(R2).
C           IN ABOVE MESSAGE, I1=         1
C           IN ABOVE MESSAGE, R1=    0.
C           IN ABOVE MESSAGE, R2=    ABOVE MESSAGE, I1=         0
C ERROR NUMBER =         6
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). ISCALE OPTION=(I1) MUST BE 1-3.
C           IN ABOVE MESSAGE, I1=         0
C ERROR NUMBER =         7
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). OFFSET PAST X(NCOLS) (I1) FOR USER-PROVIDED  COLUMN SCALING
C MUST BE POSITIVE.
C           IN ABOVE MESSAGE, I1=         0
C ERROR NUMBER =         8
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). EACH PROVIDED COL. SCALE FACTOR MUST BE POSITIVE.
C COMPONENT (I1) NOW = (R1).
C           IN ABOVE MESSAGE, I1=        ND. .LE. MDW=(I2).
C           IN ABOVE MESSAGE, I1=         1
C           IN ABOVE MESSAGE, I2=         0
C ERROR NUMBER =        10
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS().THE ROW DIMENSION OF W(,)=(I1) MUST BE .GE.THE NUMBER OF ROWS=
C (I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         1
C ERROR NUMBER =        11
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). THE COLUMN DIMENSION OF W(,)=(I1) MUST BE .GE. NCOLS+1=(I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         2
C ERROR NUMBER =        12
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS().THE DIMENSIONS OF THE ARRAYS BL(),BU(), AND IND()=(I1) MUST BE
C .GE. NCOLS=(I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         1
C ERROR NUMBER =        13
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). THE DIMENSION OF X()=(I1) MUST BE .GE. THE REQD. LENGTH=(I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         2
C ERROR NUMBER =        14
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). THE DIMENSION OF RW()=(I1) MUST BE .GE. 5*NCOLS=(I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         3
C ERROR NUMBER =        15
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS() THE DIMENSION OF IW()=(I1) MUST BE .GE. 2*NCOLS=(I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         2
C ERROR NUMBER =        16
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS() THE DIMENSION OF IOPT()=(I1) MUST BE .GE. THE REQD. LEN.=(I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         1
C ERROR NUMBER =        17
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C***REFERENCES  R. J. Hanson, Linear least squares with bounds and
C                 linear constraints, Report SAND82-1517, Sandia
C                 Laboratories, August 1982.
C***ROUTINES CALLED  DBOLSM, DCOPY, DNRM2, DROT, DROTG, IDAMAX, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   821220  DATE WRITTEN
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DBOLS
C
C     SOLVE LINEAR LEAST SQUARES SYSTEM WITH BOUNDS ON
C     SELECTED VARIABLES.
C     REVISED 850329-1400
C     REVISED YYMMDD-HHMM
C     TO CHANGE THIS SUBPROGRAM FROM SINGLE TO DOUBLE PRECISION BEGIN
C     EDITING AT THE CARD 'C++'.
C     CHANGE THIS SUBPROGRAM NAME TO DBOLS AND THE STRINGS
C     /SCOPY/ TO /DCOPY/, /SBOL/ TO /DBOL/,
C     /SNRM2/ TO /DNRM2/, /ISAMAX/ TO /IDAMAX/,
C     /SROTG/ TO /DROTG/, /SROT/ TO /DROT/, /E0/ TO /D0/,
C     /REAL            / TO /DOUBLE PRECISION/.
C ++
      DOUBLE PRECISION W(MDW,*),BL(*),BU(*),X(*),RW(*)
      DOUBLE PRECISION SC, SS, ONE, DNRM2, RNORM, ZERO
C
C     THIS VARIABLE SHOULD REMAIN TYPE REAL.
      INTEGER IND(*),IOPT(*),IW(*)
      LOGICAL CHECKL
      CHARACTER*8 XERN1, XERN2
      CHARACTER*16 XERN3, XERN4
      SAVE IGO,LOCACC,LOPT,ISCALE
      DATA IGO/0/
C***FIRST EXECUTABLE STATEMENT  DBOLS
      NERR = 0
      MODE = 0
      IF (IGO.EQ.0) THEN
C     DO(CHECK VALIDITY OF INPUT DATA)
C     PROCEDURE(CHECK VALIDITY OF INPUT DATA)
C
C     SEE THAT MDW IS .GT.0. GROSS CHECK ONLY.
          IF (MDW.LE.0) THEN
              WRITE (XERN1, '(I8)') MDW
              CALL XERMSG ('SLATEC', 'DBOLS', 'MDW = ' // XERN1 //
     *           ' MUST BE POSITIVE.', 2, 1)
C     DO(RETURN TO USER PROGRAM UNIT)
              GO TO 190
          ENDIF
C
C     SEE THAT NUMBER OF UNKNOWNS IS POSITIVE.
          IF (NCOLS.LE.0) THEN
              WRITE (XERN1, '(I8)') NCOLS
              CALL XERMSG ('SLATEC', 'DBOLS', 'NCOLS = ' // XERN1 //
     *           ' THE NO. OF VARIABLES MUST BE POSITIVE.', 3, 1)
C     DO(RETURN TO USER PROGRAM UNIT)
              GO TO 190
          ENDIF
C
C     SEE THAT CONSTRAINT INDICATORS ARE ALL WELL-DEFINED.
          DO 10 J = 1,NCOLS
              IF (IND(J).LT.1 .OR. IND(J).GT.4) THEN
                  WRITE (XERN1, '(I8)') J
                  WRITE (XERN2, '(I8)') IND(J)
                  CALL XERMSG ('SLATEC', 'DBOLS', 'IND(' // XERN1 //
     *               ') = ' // XERN2 // ' MUST BE 1-4.', 4, 1)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              ENDIF
   10     CONTINUE
C
C     SEE THAT BOUNDS ARE CONSISTENT.
          DO 20 J = 1,NCOLS
              IF (IND(J).EQ.3) THEN
                  IF (BL(J).GT.BU(J)) THEN
                      WRITE (XERN1, '(I8)') J
                      WRITE (XERN3, '(1PE15.6)') BL(J)
                      WRITE (XERN4, '(1PE15.6)') BU(J)
                      CALL XERMSG ('SLATEC', 'DBOLS', 'BOUND BL(' //
     *                   XERN1 // ') = ' // XERN3 // ' IS .GT. BU(' //
     *                   XERN1 // ') = ' // XERN4, 5, 1)
C     DO(RETURN TO USER PROGRAM UNIT)
                      GO TO 190
                  ENDIF
              ENDIF
   20     CONTINUE
C     END PROCEDURE
C     DO(PROCESS OPTION ARRAY)
C     PROCEDURE(PROCESS OPTION ARRAY)
          ZERO = 0.D0
          ONE = 1.D0
          CHECKL = .FALSE.
          LENX = NCOLS
          ISCALE = 1
          IGO = 2
          LOPT = 0
          LP = 0
          LDS = 0
   30     CONTINUE
          LP = LP + LDS
          IP = IOPT(LP+1)
          JP = ABS(IP)
C
C     TEST FOR NO MORE OPTIONS.
          IF (IP.EQ.99) THEN
              IF (LOPT.EQ.0) LOPT = LP + 1
              GO TO 50
          ELSE IF (JP.EQ.99) THEN
              LDS = 1
              GO TO 30
          ELSE IF (JP.EQ.1) THEN
              IF (IP.GT.0) THEN
C
C     SET UP DIRECTION FLAG, ROW STACKING POINTER
C     LOCATION, AND LOCATION FOR NUMBER OF NEW ROWS.
                  LOCACC = LP + 2
C
C                  IOPT(LOCACC-1)=OPTION NUMBER FOR SEQ. ACCUMULATION.
C     CONTENTS..   IOPT(LOCACC  )=USER DIRECTION FLAG, 1 OR 2.
C                  IOPT(LOCACC+1)=ROW STACKING POINTER.
C                  IOPT(LOCACC+2)=NUMBER OF NEW ROWS TO PROCESS.
C     USER ACTION WITH THIS OPTION..
C      (SET UP OPTION DATA FOR SEQ. ACCUMULATION IN IOPT(*).
C      MUST ALSO START PROCESS WITH IOPT(LOCACC)=1.)
C      (MOVE BLOCK OF EQUATIONS INTO W(*,*)  STARTING AT FIRST
C       ROW OF W(*,*).  SET IOPT(LOCACC+2)=NO. OF ROWS IN BLOCK.)
C              LOOP
C              CALL DBOLS()
C
C                  IF(IOPT(LOCACC) .EQ. 1) THEN
C                      STACK EQUAS., STARTING AT ROW IOPT(LOCACC+1),
C                       INTO W(*,*).
C                       SET IOPT(LOCACC+2)=NO. OF EQUAS.
C                      IF LAST BLOCK OF EQUAS., SET IOPT(LOCACC)=2.
C                  ELSE IF IOPT(LOCACC) .EQ. 2) THEN
C                      (PROCESS IS OVER. EXIT LOOP.)
C                  ELSE
C                      (ERROR CONDITION. SHOULD NOT HAPPEN.)
C                  END IF
C              END LOOP
C              SET IOPT(LOCACC-1)=-OPTION NUMBER FOR SEQ. ACCUMULATION.
C              CALL DBOLS( )
                  IOPT(LOCACC+1) = 1
                  IGO = 1
              ENDIF
              LDS = 4
              GO TO 30
          ELSE IF (JP.EQ.2) THEN
              IF (IP.GT.0) THEN
C
C     GET ACTUAL LENGTHS OF ARRAYS FOR CHECKING AGAINST NEEDS.
                  LOCDIM = LP + 2
C
C     LMDW.GE.MROWS
C     LNDW.GE.NCOLS+1
C     LLB .GE.NCOLS
C     LLX .GE.NCOLS+EXTRA REQD. IN OPTIONS.
C     LLRW.GE.5*NCOLS
C     LLIW.GE.2*NCOLS
C     LIOP.GE. AMOUNT REQD. FOR IOPTION ARRAY.
                  LMDW = IOPT(LOCDIM)
                  LNDW = IOPT(LOCDIM+1)
                  LLB = IOPT(LOCDIM+2)
                  LLX = IOPT(LOCDIM+3)
                  LLRW = IOPT(LOCDIM+4)
                  LLIW = IOPT(LOCDIM+5)
                  LIOPT = IOPT(LOCDIM+6)
                  CHECKL = .TRUE.
              ENDIF
              LDS = 8
              GO TO 30
C
C     OPTION TO MODIFY THE COLUMN SCALING.
          ELSE IF (JP.EQ.3) THEN
              IF (IP.GT.0) THEN
                  ISCALE = IOPT(LP+2)
C
C     SEE THAT ISCALE IS 1 THRU 3.
                  IF (ISCALE.LT.1 .OR. ISCALE.GT.3) THEN
                      WRITE (XERN1, '(I8)') ISCALE
                      CALL XERMSG ('SLATEC', 'DBOLS', 'ISCALE OPTION = '
     *                   // XERN1 // ' MUST BE 1-3', 7, 1)
C     DO(RETURN TO USER PROGRAM UNIT)
                      GO TO 190
                  ENDIF
              ENDIF
              LDS = 2
C     CYCLE FOREVER
              GO TO 30
C
C     IN THIS OPTION THE USER HAS PROVIDED SCALING.  THE
C     SCALE FACTORS FOR THE COLUMNS BEGIN IN X(NCOLS+IOPT(LP+2)).
          ELSE IF (JP.EQ.4) THEN
              IF (IP.GT.0) THEN
                  ISCALE = 4
                  IF (IOPT(LP+2).LE.0) THEN
                      WRITE (XERN1, '(I8)') IOPT(LP+2)
                      CALL XERMSG ('SLATEC', 'DBOLS',
     *                   'OFFSET PAST X(NCOLS) (' // XERN1 //
     *                   ') FOR USER-PROVIDED COLUMN SCALING MUST ' //
     *                   'BE POSITIVE.',  8, 1)
C     DO(RETURN TO USER PROGRAM UNIT)
                      GO TO 190
                  ENDIF
                  CALL DCOPY(NCOLS,X(NCOLS+IOPT(LP+2)),1,RW,1)
                  LENX = LENX + NCOLS
                  DO 40 J = 1,NCOLS
                      IF (RW(J).LE.ZERO) THEN
                          WRITE (XERN1, '(I8)') J
                          WRITE (XERN3, '(1PE15.6)') RW(J)
                          CALL XERMSG ('SLATEC', 'DBOLS',
     *                       'EACH PROVIDED COLUMN SCALE FACTOR ' //
     *                       'MUST BE POSITIVE.$$COMPONENT ' // XERN1 //
     *                       ' NOW = ' // XERN3, 9, 1)
                          GO TO 190
                      ENDIF
   40             CONTINUE
              ENDIF
              LDS = 2
C     CYCLE FOREVER
              GO TO 30
C
C     IN THIS OPTION AN OPTION ARRAY IS PROVIDED TO DBOLSM().
          ELSE IF (JP.EQ.5) THEN
              IF (IP.GT.0) THEN
                  LOPT = IOPT(LP+2)
              ENDIF
              LDS = 2
C     CYCLE FOREVER
              GO TO 30
C
C     THIS OPTION USES THE NEXT LOC OF IOPT(*) AS AN
C     INCREMENT TO SKIP.
          ELSE IF (JP.EQ.6) THEN
              IF (IP.GT.0) THEN
                  LP = IOPT(LP+2) - 1
                  LDS = 0
              ELSE
                  LDS = 2
              ENDIF
C     CYCLE FOREVER
              GO TO 30
C
C     NO VALID OPTION NUMBER WAS NOTED. THIS IS AN ERROR CONDITION.
          ELSE
              WRITE (XERN1, '(I8)') JP
              CALL XERMSG ('SLATEC', 'DBOLS', 'THE OPTION NUMBER = ' //
     *           XERN1 // ' IS NOT DEFINED.', 6, 1)
C     DO(RETURN TO USER PROGRAM UNIT)
              GO TO 190
          ENDIF
   50     CONTINUE
C     END PROCEDURE
          IF (CHECKL) THEN
C     DO(CHECK LENGTHS OF ARRAYS)
C     PROCEDURE(CHECK LENGTHS OF ARRAYS)
C
C     THIS FEATURE ALLOWS THE USER TO MAKE SURE THAT THE
C     ARRAYS ARE LONG ENOUGH FOR THE INTENDED PROBLEM SIZE AND USE.
              IF (LMDW.LT.MROWS) THEN
                  WRITE (XERN1, '(I8)') LMDW
                  WRITE (XERN2, '(I8)') MROWS
                  CALL XERMSG ('SLATEC', 'DBOLS',
     *               'THE ROW DIMENSION OF W(,) = ' // XERN1 //
     *               ' MUST BE .GE. THE NUMBER OF ROWS = ' // XERN2,
     *               11, 1)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              ENDIF
              IF (LNDW.LT.NCOLS+1) THEN
                  WRITE (XERN1, '(I8)') LNDW
                  WRITE (XERN2, '(I8)') NCOLS+1
                  CALL XERMSG ('SLATEC', 'DBOLS',
     *               'THE COLUMN DIMENSION OF W(,) = ' // XERN1 //
     *               ' MUST BE .GE. NCOLS+1 = ' // XERN2, 12, 1)
                  GO TO 190
              ENDIF
              IF (LLB.LT.NCOLS) THEN
                  WRITE (XERN1, '(I8)') LLB
                  WRITE (XERN2, '(I8)') NCOLS
                  CALL XERMSG ('SLATEC', 'DBOLS',
     *           'THE DIMENSIONS OF THE ARRAYS BL(), BU(), AND IND() = '
     *               // XERN1 // ' MUST BE .GE. NCOLS = ' // XERN2,
     *               13, 1)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              ENDIF
              IF (LLX.LT.LENX) THEN
                  WRITE (XERN1, '(I8)') LLX
                  WRITE (XERN2, '(I8)') LENX
                  CALL XERMSG ('SLATEC', 'DBOLS',
     *               'THE DIMENSION OF X() = ' // XERN1 //
     *               ' MUST BE .GE. THE REQUIRED LENGTH = ' // XERN2,
     *               14, 1)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              ENDIF
              IF (LLRW.LT.5*NCOLS) THEN
                  WRITE (XERN1, '(I8)') LLRW
                  WRITE (XERN2, '(I8)') 5*NCOLS
                  CALL XERMSG ('SLATEC', 'DBOLS',
     *               'THE DIMENSION OF RW() = ' // XERN1 //
     *               ' MUST BE .GE. 5*NCOLS = ' // XERN2, 15, 1)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              ENDIF
              IF (LLIW.LT.2*NCOLS) THEN
                  WRITE (XERN1, '(I8)') LLIW
                  WRITE (XERN2, '(I8)') 2*NCOLS
                  CALL XERMSG ('SLATEC', 'DBOLS',
     *               'THE DIMENSION OF IW() = ' // XERN1 //
     *               ' MUST BE .GE. 2*NCOLS = ' // XERN2, 16, 1)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              ENDIF
              IF (LIOPT.LT.LP+1) THEN
                  WRITE (XERN1, '(I8)') LIOPT
                  WRITE (XERN2, '(I8)') LP+1
                  CALL XERMSG ('SLATEC', 'DBOLS',
     *               'THE DIMENSION OF IOPT() = ' // XERN1 //
     *               ' MUST BE .GE. THE REQUIRED LEN = ' // XERN2, 17,1)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              ENDIF
C     END PROCEDURE
          ENDIF
      ENDIF
      GO TO (60,90),IGO
      GO TO 180
C
C     GO BACK TO THE USER FOR ACCUMULATION OF LEAST SQUARES
C     EQUATIONS AND DIRECTIONS TO QUIT PROCESSING.
C     CASE 1
   60 CONTINUE
C     DO(ACCUMULATE LEAST SQUARES EQUATIONS)
C     PROCEDURE(ACCUMULATE LEAST SQUARES EQUATIONS)
      MROWS = IOPT(LOCACC+1) - 1
      INROWS = IOPT(LOCACC+2)
      MNEW = MROWS + INROWS
      IF (MNEW.LT.0 .OR. MNEW.GT.MDW) THEN
          WRITE (XERN1, '(I8)') MNEW
          WRITE (XERN2, '(I8)') MDW
          CALL XERMSG ('SLATEC', 'DBOLS', 'NO. OF ROWS = ' // XERN1 //
     *       ' MUST BE .GE. 0 .AND. .LE. MDW = ' // XERN2, 10, 1)
C     DO(RETURN TO USER PROGRAM UNIT)
          GO TO 190
      ENDIF
      DO 80 J = 1,MIN(NCOLS+1,MNEW)
          DO 70 I = MNEW,MAX(MROWS,J) + 1,-1
              IBIG = IDAMAX(I-J,W(J,J),1) + J - 1
C
C     PIVOT FOR INCREASED STABILITY.
              CALL DROTG(W(IBIG,J),W(I,J),SC,SS)
              CALL DROT(NCOLS+1-J,W(IBIG,J+1),MDW,W(I,J+1),MDW,SC,SS)
              W(I,J) = ZERO
   70     CONTINUE
   80 CONTINUE
      MROWS = MIN(NCOLS+1,MNEW)
      IOPT(LOCACC+1) = MROWS + 1
      IGO = IOPT(LOCACC)
C     END PROCEDURE
      IF (IGO.EQ.2) THEN
          IGO = 0
      ENDIF
      GO TO 180
C     CASE 2
   90 CONTINUE
C     DO(INITIALIZE VARIABLES AND DATA VALUES)
C     PROCEDURE(INITIALIZE VARIABLES AND DATA VALUES)
      DO 150 J = 1,NCOLS
          GO TO (100,110,120,130),ISCALE
          GO TO 140
  100     CONTINUE
C     CASE 1
C
C     THIS IS THE NOMINAL SCALING. EACH NONZERO
C     COL. HAS MAX. NORM EQUAL TO ONE.
          IBIG = IDAMAX(MROWS,W(1,J),1)
          RW(J) = ABS(W(IBIG,J))
          IF (RW(J).EQ.ZERO) THEN
              RW(J) = ONE
          ELSE
              RW(J) = ONE/RW(J)
          ENDIF
          GO TO 140
  110     CONTINUE
C     CASE 2
C
C     THIS CHOICE OF SCALING MAKES EACH NONZERO COLUMN
C     HAVE EUCLIDEAN LENGTH EQUAL TO ONE.
          RW(J) = DNRM2(MROWS,W(1,J),1)
          IF (RW(J).EQ.ZERO) THEN
              RW(J) = ONE
          ELSE
              RW(J) = ONE/RW(J)
          ENDIF
          GO TO 140
  120     CONTINUE
C     CASE 3
C
C     THIS CASE EFFECTIVELY SUPPRESSES SCALING BY SETTING
C     THE SCALING MATRIX TO THE IDENTITY MATRIX.
          RW(1) = ONE
          CALL DCOPY(NCOLS,RW,0,RW,1)
          GO TO 160
  130     CONTINUE
C     CASE 4
          GO TO 160
  140     CONTINUE
  150 CONTINUE
  160 CONTINUE
C     END PROCEDURE
C     DO(SOLVE BOUNDED LEAST SQUARES PROBLEM)
C     PROCEDURE(SOLVE BOUNDED LEAST SQUARES PROBLEM)
C
C     INITIALIZE IBASIS(*), J=1,NCOLS, AND IBB(*), J=1,NCOLS,
C     TO =J,AND =1, FOR USE IN DBOLSM( ).
      DO 170 J = 1,NCOLS
          IW(J) = J
          IW(J+NCOLS) = 1
          RW(3*NCOLS+J) = BL(J)
          RW(4*NCOLS+J) = BU(J)
  170 CONTINUE
      CALL DBOLSM(W,MDW,MROWS,NCOLS,RW(3*NCOLS+1),RW(4*NCOLS+1),IND,
     .            IOPT(LOPT),X,RNORM,MODE,RW(NCOLS+1),RW(2*NCOLS+1),RW,
     .            IW,IW(NCOLS+1))
C     END PROCEDURE
      IGO = 0
  180 CONTINUE
      RETURN
C     PROCEDURE(RETURN TO USER PROGRAM UNIT)
  190 IF(MODE.GE.0)MODE = -NERR
      IGO = 0
      RETURN
C     END PROCEDURE
      END
