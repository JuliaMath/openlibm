*DECK SPLP
      SUBROUTINE SPLP (USRMAT, MRELAS, NVARS, COSTS, PRGOPT, DATTRV, BL,
     +   BU, IND, INFO, PRIMAL, DUALS, IBASIS, WORK, LW, IWORK, LIW)
C***BEGIN PROLOGUE  SPLP
C***PURPOSE  Solve linear programming problems involving at
C            most a few thousand constraints and variables.
C            Takes advantage of sparsity in the constraint matrix.
C***LIBRARY   SLATEC
C***CATEGORY  G2A2
C***TYPE      SINGLE PRECISION (SPLP-S, DSPLP-D)
C***KEYWORDS  LINEAR CONSTRAINTS, LINEAR OPTIMIZATION,
C             LINEAR PROGRAMMING, LP, SPARSE CONSTRAINTS
C***AUTHOR  Hanson, R. J., (SNLA)
C           Hiebert, K. L., (SNLA)
C***DESCRIPTION
C
C     These are the short usage instructions; for details about
C     other features, options and methods for defining the matrix
C     A, see the extended usage instructions which are contained in
C     the Long Description section below.
C
C   |------------|
C   |Introduction|
C   |------------|
C     The subprogram SPLP( ) solves a linear optimization problem.
C     The problem statement is as follows
C
C                         minimize (transpose of costs)*x
C                         subject to A*x=w.
C
C     The entries of the unknowns x and w may have simple lower or
C     upper bounds (or both), or be free to take on any value.  By
C     setting the bounds for x and w, the user is imposing the con-
C     straints of the problem.  The matrix A has MRELAS rows and
C     NVARS columns.  The vectors costs, x, and w respectively
C     have NVARS, NVARS, and MRELAS number of entries.
C
C     The input for the problem includes the problem dimensions,
C     MRELAS and NVARS, the array COSTS(*), data for the matrix
C     A, and the bound information for the unknowns x and w, BL(*),
C     BU(*), and IND(*).  Only the nonzero entries of the matrix A
C     are passed to SPLP( ).
C
C     The output from the problem (when output flag INFO=1) includes
C     optimal values for x and w in PRIMAL(*), optimal values for
C     dual variables of the equations A*x=w and the simple bounds
C     on x in  DUALS(*), and the indices of the basic columns,
C     IBASIS(*).
C
C  |------------------------------|
C  |Fortran Declarations Required:|
C  |------------------------------|
C
C     DIMENSION COSTS(NVARS),PRGOPT(*),DATTRV(*),
C    *BL(NVARS+MRELAS),BU(NVARS+MRELAS),IND(NVARS+MRELAS),
C    *PRIMAL(NVARS+MRELAS),DUALS(MRELAS+NVARS),IBASIS(NVARS+MRELAS),
C    *WORK(LW),IWORK(LIW)
C
C     EXTERNAL USRMAT
C
C     The dimensions of PRGOPT(*) and DATTRV(*) must be at least 1.
C     The exact lengths will be determined by user-required options and
C     data transferred to the subprogram USRMAT( ).
C
C     The values of LW and LIW, the lengths of the arrays WORK(*)
C     and IWORK(*), must satisfy the inequalities
C
C               LW .GE. 4*NVARS+ 8*MRELAS+LAMAT+  LBM
C               LIW.GE.   NVARS+11*MRELAS+LAMAT+2*LBM
C
C     It is an error if they do not both satisfy these inequalities.
C     (The subprogram will inform the user of the required lengths
C     if either LW or LIW is wrong.)  The values of LAMAT and LBM
C     nominally are
C
C               LAMAT=4*NVARS+7
C     and       LBM  =8*MRELAS
C
C     LAMAT determines the length of the sparse matrix storage area.
C     The value of LBM determines the amount of storage available
C     to decompose and update the active basis matrix.
C
C  |------|
C  |Input:|
C  |------|
C
C     MRELAS,NVARS
C     ------------
C     These parameters are respectively the number of constraints (the
C     linear relations A*x=w that the unknowns x and w are to satisfy)
C     and the number of entries in the vector x.  Both must be .GE. 1.
C     Other values are errors.
C
C     COSTS(*)
C     --------
C     The NVARS entries of this array are the coefficients of the
C     linear objective function.  The value COSTS(J) is the
C     multiplier for variable J of the unknown vector x.  Each
C     entry of this array must be defined.
C
C     USRMAT
C     ------
C     This is the name of a specific subprogram in the SPLP( ) package
C     used to define the matrix A.  In this usage mode of SPLP( )
C     the user places the nonzero entries of A in the
C     array DATTRV(*) as given in the description of that parameter.
C     The name USRMAT must appear in a Fortran EXTERNAL statement.
C
C     DATTRV(*)
C     ---------
C     The array DATTRV(*) contains data for the matrix A as follows:
C     Each column (numbered J) requires (floating point) data con-
C     sisting of the value (-J) followed by pairs of values.  Each pair
C     consists of the row index immediately followed by the value
C     of the matrix at that entry.  A value of J=0 signals that there
C     are no more columns.  The required length of
C     DATTRV(*) is 2*no. of nonzeros + NVARS + 1.
C
C     BL(*),BU(*),IND(*)
C     ------------------
C     The values of IND(*) are input parameters that define
C     the form of the bounds for the unknowns x and w.  The values for
C     the bounds are found in the arrays BL(*) and BU(*) as follows.
C
C     For values of J between 1 and NVARS,
C          if IND(J)=1, then X(J) .GE. BL(J); BU(J) is not used.
C          if IND(J)=2, then X(J) .LE. BU(J); BL(J) is not used.
C          if IND(J)=3, then BL(J) .LE. X(J) .LE. BU(J),(BL(J)=BU(J) ok)
C          if IND(J)=4, then X(J) is free to have any value,
C          and BL(J), BU(J) are not used.
C
C     For values of I between NVARS+1 and NVARS+MRELAS,
C          if IND(I)=1, then W(I-NVARS) .GE. BL(I); BU(I) is not used.
C          if IND(I)=2, then W(I-NVARS) .LE. BU(I); BL(I) is not used.
C          if IND(I)=3, then BL(I) .LE. W(I-NVARS) .LE. BU(I),
C          (BL(I)=BU(I) is ok).
C          if IND(I)=4, then W(I-NVARS) is free to have any value,
C          and BL(I), BU(I) are not used.
C
C     A value of IND(*) not equal to 1,2,3 or 4 is an error.  When
C     IND(I)=3, BL(I) must be .LE. BU(I).  The condition BL(I).GT.
C     BU(I) indicates infeasibility and is an error.
C
C     PRGOPT(*)
C     ---------
C     This array is used to redefine various parameters within SPLP( ).
C     Frequently, perhaps most of the time, a user will be satisfied
C     and obtain the solutions with no changes to any of these
C     parameters.  To try this, simply set PRGOPT(1)=1.E0.
C
C     For users with more sophisticated needs, SPLP( ) provides several
C     options that may be used to take advantage of more detailed
C     knowledge of the problem or satisfy other utilitarian needs.
C     The complete description of how to use this option array to
C     utilize additional subprogram features is found under the
C     heading  of SPLP( ) Subprogram Options in the Extended
C     Usage Instructions.
C
C     Briefly, the user should note the following value of the parameter
C     KEY and the corresponding task or feature desired before turning
C     to that document.
C
C     Value     Brief Statement of Purpose for Option
C     of KEY
C     ------    -------------------------------------
C     50        Change from a minimization problem to a
C               maximization problem.
C     51        Change the amount of printed output.
C               Normally, no printed output is obtained.
C     52        Redefine the line length and precision used
C               for the printed output.
C     53        Redefine the values of LAMAT and LBM that
C               were discussed above under the heading
C               Fortran Declarations Required.
C     54        Redefine the unit number where pages of the sparse
C               data matrix A are stored.  Normally, the unit
C               number is 1.
C     55        A computation, partially completed, is
C               being continued.  Read the up-to-date
C               partial results from unit number 2.
C     56        Redefine the unit number where the partial results
C               are stored.  Normally, the unit number is 2.
C     57        Save partial results on unit 2 either after
C               maximum iterations or at the optimum.
C     58        Redefine the value for the maximum number of
C               iterations.  Normally, the maximum number of
C               iterations is 3*(NVARS+MRELAS).
C     59        Provide SPLP( ) with a starting (feasible)
C               nonsingular basis.  Normally, SPLP( ) starts
C               with the identity matrix columns corresponding
C               to the vector w.
C     60        The user has provided scale factors for the
C               columns of A.  Normally, SPLP( ) computes scale
C               factors that are the reciprocals of the max. norm
C               of each column.
C     61        The user has provided a scale factor
C               for the vector costs.  Normally, SPLP( ) computes
C               a scale factor equal to the reciprocal of the
C               max. norm of the vector costs after the column
C               scaling for the data matrix has been applied.
C     62        Size parameters, namely the smallest and
C               largest magnitudes of nonzero entries in
C               the matrix A, are provided.  Values noted
C               outside this range are to be considered errors.
C     63        Redefine the tolerance required in
C               evaluating residuals for feasibility.
C               Normally, this value is set to RELPR,
C               where RELPR = relative precision of the arithmetic.
C     64        Change the criterion for bringing new variables
C               into the basis from the steepest edge (best
C               local move) to the minimum reduced cost.
C     65        Redefine the value for the number of iterations
C               between recalculating the error in the primal
C               solution.  Normally, this value is equal to ten.
C     66        Perform "partial pricing" on variable selection.
C               Redefine the value for the number of negative
C               reduced costs to compute (at most) when finding
C               a variable to enter the basis.  Normally this
C               value is set to NVARS.  This implies that no
C               "partial pricing" is used.
C     67        Adjust the tuning factor (normally one) to apply
C               to the primal and dual error estimates.
C     68        Pass  information to the  subprogram  FULMAT(),
C               provided with the SPLP() package, so that a Fortran
C               two-dimensional array can be used as the argument
C               DATTRV(*).
C     69        Pass an absolute tolerance to use for the feasibility
C               test when the usual relative error test indicates
C               infeasibility.  The nominal value of this tolerance,
C               TOLABS, is zero.
C
C
C  |---------------|
C  |Working Arrays:|
C  |---------------|
C
C     WORK(*),LW,
C     IWORK(*),LIW
C     ------------
C     The arrays WORK(*) and IWORK(*) are respectively floating point
C     and type INTEGER working arrays for SPLP( ) and its
C     subprograms.  The lengths of these arrays are respectively
C     LW and LIW.  These parameters must satisfy the inequalities
C     noted above under the heading "Fortran Declarations Required:"
C     It is an error if either value is too small.
C
C  |----------------------------|
C  |Input/Output files required:|
C  |----------------------------|
C
C     Fortran unit 1 is used by SPLP( ) to store the sparse matrix A
C     out of high-speed memory.  A crude
C     upper bound for the amount of information written on unit 1
C     is 6*nz, where nz is the number of nonzero entries in A.
C
C  |-------|
C  |Output:|
C  |-------|
C
C     INFO,PRIMAL(*),DUALS(*)
C     -----------------------
C     The integer flag INFO indicates why SPLP( ) has returned to the
C     user.  If INFO=1 the solution has been computed.  In this case
C     X(J)=PRIMAL(J) and W(I)=PRIMAL(I+NVARS).  The dual variables
C     for the equations A*x=w are in the array DUALS(I)=dual for
C     equation number I.  The dual value for the component X(J) that
C     has an upper or lower bound (or both) is returned in
C     DUALS(J+MRELAS).  The only other values for INFO are .LT. 0.
C     The meaning of these values can be found by reading
C     the diagnostic message in the output file, or by looking for
C     error number = (-INFO) in the Extended Usage Instructions
C     under the heading:
C
C          List of SPLP( ) Error and Diagnostic Messages.
C
C     BL(*),BU(*),IND(*)
C     ------------------
C     These arrays are output parameters only under the (unusual)
C     circumstances where the stated problem is infeasible, has an
C     unbounded optimum value, or both.  These respective conditions
C     correspond to INFO=-1,-2 or -3.    See the Extended
C     Usage Instructions for further details.
C
C     IBASIS(I),I=1,...,MRELAS
C     ------------------------
C     This array contains the indices of the variables that are
C     in the active basis set at the solution (INFO=1).  A value
C     of IBASIS(I) between 1 and NVARS corresponds to the variable
C     X(IBASIS(I)).  A value of IBASIS(I) between NVARS+1 and NVARS+
C     MRELAS corresponds to the variable W(IBASIS(I)-NVARS).
C
C *Long Description:
C
C     SUBROUTINE SPLP(USRMAT,MRELAS,NVARS,COSTS,PRGOPT,DATTRV,
C    *           BL,BU,IND,INFO,PRIMAL,DUALS,IBASIS,WORK,LW,IWORK,LIW)
C
C     |------------|
C     |Introduction|
C     |------------|
C     The subprogram SPLP( ) solves a linear optimization problem.
C     The problem statement is as follows
C
C                         minimize (transpose of costs)*x
C                         subject to A*x=w.
C
C     The entries of the unknowns x and w may have simple lower or
C     upper bounds (or both), or be free to take on any value.  By
C     setting the bounds for x and w, the user is imposing the con-
C     straints of the problem.
C
C     (The problem may also be stated as a maximization
C     problem.  This is done by means of input in the option array
C     PRGOPT(*).)  The matrix A has MRELAS rows and NVARS columns.  The
C     vectors costs, x, and w respectively have NVARS, NVARS, and
C     MRELAS number of entries.
C
C     The input for the problem includes the problem dimensions,
C     MRELAS and NVARS, the array COSTS(*), data for the matrix
C     A, and the bound information for the unknowns x and w, BL(*),
C     BU(*), and IND(*).
C
C     The output from the problem (when output flag INFO=1) includes
C     optimal values for x and w in PRIMAL(*), optimal values for
C     dual variables of the equations A*x=w and the simple bounds
C     on x in  DUALS(*), and the indices of the basic columns in
C     IBASIS(*).
C
C  |------------------------------|
C  |Fortran Declarations Required:|
C  |------------------------------|
C
C     DIMENSION COSTS(NVARS),PRGOPT(*),DATTRV(*),
C    *BL(NVARS+MRELAS),BU(NVARS+MRELAS),IND(NVARS+MRELAS),
C    *PRIMAL(NVARS+MRELAS),DUALS(MRELAS+NVARS),IBASIS(NVARS+MRELAS),
C    *WORK(LW),IWORK(LIW)
C
C     EXTERNAL USRMAT (or 'NAME', if user provides the subprogram)
C
C     The dimensions of PRGOPT(*) and DATTRV(*) must be at least 1.
C     The exact lengths will be determined by user-required options and
C     data transferred to the subprogram USRMAT( ) ( or 'NAME').
C
C     The values of LW and LIW, the lengths of the arrays WORK(*)
C     and IWORK(*), must satisfy the inequalities
C
C               LW .GE. 4*NVARS+ 8*MRELAS+LAMAT+  LBM
C               LIW.GE.   NVARS+11*MRELAS+LAMAT+2*LBM
C
C     It is an error if they do not both satisfy these inequalities.
C     (The subprogram will inform the user of the required lengths
C     if either LW or LIW is wrong.)  The values of LAMAT and LBM
C     nominally are
C
C               LAMAT=4*NVARS+7
C     and       LBM  =8*MRELAS
C
C     These values will be as shown unless the user changes them by
C     means of input in the option array PRGOPT(*).  The value of LAMAT
C     determines the length of the sparse matrix "staging" area.
C     For reasons of efficiency the user may want to increase the value
C     of LAMAT.  The value of LBM determines the amount of storage
C     available to decompose and update the active basis matrix.
C     Due to exhausting the working space because of fill-in,
C     it may be necessary for the user to increase the value of LBM.
C     (If this situation occurs an informative diagnostic is printed
C     and a value of INFO=-28 is obtained as an output parameter.)
C
C  |------|
C  |Input:|
C  |------|
C
C     MRELAS,NVARS
C     ------------
C     These parameters are respectively the number of constraints (the
C     linear relations A*x=w that the unknowns x and w are to satisfy)
C     and the number of entries in the vector x.  Both must be .GE. 1.
C     Other values are errors.
C
C     COSTS(*)
C     --------
C     The NVARS entries of this array are the coefficients of the
C     linear objective function.  The value COSTS(J) is the
C     multiplier for variable J of the unknown vector x.  Each
C     entry of this array must be defined.  This array can be changed
C     by the user between restarts.  See options with KEY=55,57 for
C     details of checkpointing and restarting.
C
C     USRMAT
C     ------
C     This is the name of a specific subprogram in the SPLP( ) package
C     that is used to define the matrix entries when this data is passed
C     to SPLP( ) as a linear array.  In this usage mode of SPLP( )
C     the user gives information about the nonzero entries of A
C     in DATTRV(*) as given under the description of that parameter.
C     The name USRMAT must appear in a Fortran EXTERNAL statement.
C     Users who are passing the matrix data with USRMAT( ) can skip
C     directly to the description of the input parameter DATTRV(*).
C     Also see option 68 for passing the constraint matrix data using
C     a standard Fortran two-dimensional array.
C
C     If the user chooses to provide a subprogram 'NAME'( ) to
C     define the matrix A, then DATTRV(*) may be used to pass floating
C     point data from the user's program unit to the subprogram
C     'NAME'( ). The content of DATTRV(*) is not changed in any way.
C
C     The subprogram 'NAME'( ) can be of the user's choice
C     but it must meet Fortran standards and it must appear in a
C     Fortran EXTERNAL statement.  The first statement of the subprogram
C     has the form
C
C          SUBROUTINE 'NAME'(I,J,AIJ, INDCAT, PRGOPT, DATTRV, IFLAG)
C
C     The variables I,J, INDCAT, IFLAG(10) are type INTEGER,
C          while  AIJ, PRGOPT(*),DATTRV(*) are type REAL.
C
C     The user interacts with the contents of IFLAG(*) to
C     direct the appropriate action.  The algorithmic steps are
C     as follows.
C
C          Test IFLAG(1).
C
C             IF(IFLAG(1).EQ.1) THEN
C
C               Initialize the necessary pointers and data
C               for defining the matrix A.  The contents
C               of IFLAG(K), K=2,...,10, may be used for
C               storage of the pointers.  This array remains intact
C               between calls to 'NAME'( ) by SPLP( ).
C               RETURN
C
C             END IF
C
C             IF(IFLAG(1).EQ.2) THEN
C
C               Define one set of values for I,J,AIJ, and INDCAT.
C               Each nonzero entry of A must be defined this way.
C               These values can be defined in any convenient order.
C               (It is most efficient to define the data by
C               columns in the order 1,...,NVARS; within each
C               column define the entries in the order 1,...,MRELAS.)
C               If this is the last matrix value to be
C               defined or updated, then set IFLAG(1)=3.
C               (When I and J are positive and respectively no larger
C               than MRELAS and NVARS, the value of AIJ is used to
C               define (or update) row I and column J of A.)
C               RETURN
C
C             END IF
C
C               END
C
C     Remarks:  The values of I and J are the row and column
C               indices for the nonzero entries of the matrix A.
C               The value of this entry is AIJ.
C               Set INDCAT=0 if this value defines that entry.
C               Set INDCAT=1 if this entry is to be updated,
C                            new entry=old entry+AIJ.
C               A value of I not between 1 and MRELAS, a value of J
C               not between 1 and NVARS, or a value of INDCAT
C               not equal to 0 or 1 are each errors.
C
C               The contents of IFLAG(K), K=2,...,10, can be used to
C               remember the status (of the process of defining the
C               matrix entries) between calls to 'NAME'( ) by SPLP( ).
C               On entry to 'NAME'( ), only the values 1 or 2 will be
C               in IFLAG(1).  More than 2*NVARS*MRELAS definitions of
C               the matrix elements is considered an error because
C               it suggests an infinite loop in the user-written
C               subprogram 'NAME'( ).  Any matrix element not
C               provided by 'NAME'( ) is defined to be zero.
C
C               The REAL arrays PRGOPT(*) and DATTRV(*) are passed as
C               arguments directly from SPLP( ) to 'NAME'( ).
C               The array PRGOPT(*) contains any user-defined program
C               options.  In this usage mode the array DATTRV(*) may
C               now contain any (type REAL) data that the user needs
C               to define the matrix A.  Both arrays PRGOPT(*) and
C               DATTRV(*) remain intact between calls to 'NAME'( )
C               by SPLP( ).
C     Here is a subprogram that communicates the matrix values for A,
C     as represented in DATTRV(*), to SPLP( ).  This subprogram,
C     called USRMAT( ), is included as part of the SPLP( ) package.
C     This subprogram 'decodes' the array DATTRV(*) and defines the
C     nonzero entries of the matrix  A for SPLP( ) to store.  This
C     listing is presented here as a guide and example
C     for the users who find it necessary to write their own subroutine
C     for this purpose.  The contents of DATTRV(*) are given below in
C     the description of that parameter.
C
C     SUBROUTINE USRMAT(I,J,AIJ, INDCAT,PRGOPT,DATTRV,IFLAG)
C     DIMENSION PRGOPT(*),DATTRV(*),IFLAG(10)
C
C     IF(IFLAG(1).EQ.1) THEN
C
C     THIS IS THE INITIALIZATION STEP.  THE VALUES OF IFLAG(K),K=2,3,4,
C     ARE RESPECTIVELY THE COLUMN INDEX, THE ROW INDEX (OR THE NEXT COL.
C     INDEX), AND THE POINTER TO THE MATRIX ENTRY'S VALUE WITHIN
C     DATTRV(*).  ALSO CHECK (DATTRV(1)=0.) SIGNIFYING NO DATA.
C          IF(DATTRV(1).EQ.0.) THEN
C          I = 0
C          J = 0
C          IFLAG(1) = 3
C          ELSE
C          IFLAG(2)=-DATTRV(1)
C          IFLAG(3)= DATTRV(2)
C          IFLAG(4)= 3
C          END IF
C
C          RETURN
C     ELSE
C          J=IFLAG(2)
C          I=IFLAG(3)
C          L=IFLAG(4)
C          IF(I.EQ.0) THEN
C
C     SIGNAL THAT ALL OF THE NONZERO ENTRIES HAVE BEEN DEFINED.
C               IFLAG(1)=3
C               RETURN
C          ELSE IF(I.LT.0) THEN
C
C     SIGNAL THAT A SWITCH IS MADE TO A NEW COLUMN.
C               J=-I
C               I=DATTRV(L)
C               L=L+1
C          END IF
C
C          AIJ=DATTRV(L)
C
C     UPDATE THE INDICES AND POINTERS FOR THE NEXT ENTRY.
C          IFLAG(2)=J
C          IFLAG(3)=DATTRV(L+1)
C          IFLAG(4)=L+2
C
C     INDCAT=0 DENOTES THAT ENTRIES OF THE MATRIX ARE ASSIGNED THE
C     VALUES FROM DATTRV(*).  NO ACCUMULATION IS PERFORMED.
C          INDCAT=0
C          RETURN
C     END IF
C     END
C
C     DATTRV(*)
C     ---------
C     If the user chooses to use the provided subprogram USRMAT( ) then
C     the array DATTRV(*) contains data for the matrix A as follows:
C     Each column (numbered J) requires (floating point) data con-
C     sisting of the value (-J) followed by pairs of values.  Each pair
C     consists of the row index immediately followed by the value
C     of the matrix at that entry.  A value of J=0 signals that there
C     are no more columns.  (See "Example of SPLP( ) Usage," below.)
C     The dimension of DATTRV(*) must be 2*no. of nonzeros
C     + NVARS + 1 in this usage.  No checking of the array
C     length is done by the subprogram package.
C
C     If the Save/Restore feature is in use (see options with
C     KEY=55,57 for details of checkpointing and restarting)
C     USRMAT( ) can be used to redefine entries of the matrix.
C     The matrix entries are redefined or overwritten.  No accum-
C     ulation is performed.
C     Any other nonzero entry of A, defined in a previous call to
C     SPLP( ), remain intact.
C
C     BL(*),BU(*),IND(*)
C     ------------------
C     The values of IND(*) are input parameters that define
C     the form of the bounds for the unknowns x and w.  The values for
C     the bounds are found in the arrays BL(*) and BU(*) as follows.
C
C     For values of J between 1 and NVARS,
C          if IND(J)=1, then X(J) .GE. BL(J); BU(J) is not used.
C          if IND(J)=2, then X(J) .LE. BU(J); BL(J) is not used.
C          if IND(J)=3, then BL(J) .LE. X(J) .LE. BU(J),(BL(J)=BU(J) ok)
C          if IND(J)=4, then X(J) is free to have any value,
C          and BL(J), BU(J) are not used.
C
C     For values of I between NVARS+1 and NVARS+MRELAS,
C          if IND(I)=1, then W(I-NVARS) .GE. BL(I); BU(I) is not used.
C          if IND(I)=2, then W(I-NVARS) .LE. BU(I); BL(I) is not used.
C          if IND(I)=3, then BL(I) .LE. W(I-NVARS) .LE. BU(I),
C          (BL(I)=BU(I) is ok).
C          if IND(I)=4, then W(I-NVARS) is free to have any value,
C          and BL(I), BU(I) are not used.
C
C     A value of IND(*) not equal to 1,2,3 or 4 is an error.  When
C     IND(I)=3, BL(I) must be .LE. BU(I).  The condition BL(I).GT.
C     BU(I) indicates infeasibility and is an error.  These
C     arrays can be changed by the user between restarts.  See
C     options with KEY=55,57 for details of checkpointing and
C     restarting.
C
C     PRGOPT(*)
C     ---------
C     This array is used to redefine various parameters within SPLP( ).
C     Frequently, perhaps most of the time, a user will be satisfied
C     and obtain the solutions with no changes to any of these
C     parameters.  To try this, simply set PRGOPT(1)=1.E0.
C
C     For users with more sophisticated needs, SPLP( ) provides several
C     options that may be used to take advantage of more detailed
C     knowledge of the problem or satisfy other utilitarian needs.
C     The complete description of how to use this option array to
C     utilize additional subprogram features is found under the
C     heading "Usage of SPLP( ) Subprogram Options."
C
C     Briefly, the user should note the following value of the parameter
C     KEY and the corresponding task or feature desired before turning
C     to that section.
C
C     Value     Brief Statement of Purpose for Option
C     of KEY
C     ------    -------------------------------------
C     50        Change from a minimization problem to a
C               maximization problem.
C     51        Change the amount of printed output.
C               Normally, no printed output is obtained.
C     52        Redefine the line length and precision used
C               for the printed output.
C     53        Redefine the values of LAMAT and LBM that
C               were discussed above under the heading
C               Fortran Declarations Required.
C     54        Redefine the unit number where pages of the sparse
C               data matrix A are stored.  Normally, the unit
C               number is 1.
C     55        A computation, partially completed, is
C               being continued.  Read the up-to-date
C               partial results from unit number 2.
C     56        Redefine the unit number where the partial results
C               are stored.  Normally, the unit number is 2.
C     57        Save partial results on unit 2 either after
C               maximum iterations or at the optimum.
C     58        Redefine the value for the maximum number of
C               iterations.  Normally, the maximum number of
C               iterations is 3*(NVARS+MRELAS).
C     59        Provide SPLP( ) with a starting (feasible)
C               nonsingular basis.  Normally, SPLP( ) starts
C               with the identity matrix columns corresponding
C               to the vector w.
C     60        The user has provided scale factors for the
C               columns of A.  Normally, SPLP( ) computes scale
C               factors that are the reciprocals of the max. norm
C               of each column.
C     61        The user has provided a scale factor
C               for the vector costs.  Normally, SPLP( ) computes
C               a scale factor equal to the reciprocal of the
C               max. norm of the vector costs after the column
C               scaling for the data matrix has been applied.
C     62        Size parameters, namely the smallest and
C               largest magnitudes of nonzero entries in
C               the matrix A, are provided.  Values noted
C               outside this range are to be considered errors.
C     63        Redefine the tolerance required in
C               evaluating residuals for feasibility.
C               Normally, this value is set to the value RELPR,
C               where RELPR = relative precision of the arithmetic.
C     64        Change the criterion for bringing new variables
C               into the basis from the steepest edge (best
C               local move) to the minimum reduced cost.
C     65        Redefine the value for the number of iterations
C               between recalculating the error in the primal
C               solution.  Normally, this value is equal to ten.
C     66        Perform "partial pricing" on variable selection.
C               Redefine the value for the number of negative
C               reduced costs to compute (at most) when finding
C               a variable to enter the basis.  Normally this
C               value is set to NVARS.  This implies that no
C               "partial pricing" is used.
C     67        Adjust the tuning factor (normally one) to apply
C               to the primal and dual error estimates.
C     68        Pass  information to the  subprogram  FULMAT(),
C               provided with the SPLP() package, so that a Fortran
C               two-dimensional array can be used as the argument
C               DATTRV(*).
C     69        Pass an absolute tolerance to use for the feasibility
C               test when the usual relative error test indicates
C               infeasibility.  The nominal value of this tolerance,
C               TOLABS, is zero.
C
C
C  |---------------|
C  |Working Arrays:|
C  |---------------|
C
C     WORK(*),LW,
C     IWORK(*),LIW
C     ------------
C     The arrays WORK(*) and IWORK(*) are respectively floating point
C     and type INTEGER working arrays for SPLP( ) and its
C     subprograms.  The lengths of these arrays are respectively
C     LW and LIW.  These parameters must satisfy the inequalities
C     noted above under the heading "Fortran Declarations Required."
C     It is an error if either value is too small.
C
C  |----------------------------|
C  |Input/Output files required:|
C  |----------------------------|
C
C     Fortran unit 1 is used by SPLP( ) to store the sparse matrix A
C     out of high-speed memory.  This direct access file is opened
C     within the package under the following two conditions.
C     1. When the Save/Restore feature is used.  2. When the
C     constraint matrix is so large that storage out of high-speed
C     memory is required.  The user may need to close unit 1
C     (with deletion from the job step) in the main program unit
C     when several calls are made to SPLP( ).  A crude
C     upper bound for the amount of information written on unit 1
C     is 6*nz, where nz is the number of nonzero entries in A.
C     The unit number may be redefined to any other positive value
C     by means of input in the option array PRGOPT(*).
C
C     Fortran unit 2 is used by SPLP( ) only when the Save/Restore
C     feature is desired.  Normally this feature is not used.  It is
C     activated by means of input in the option array PRGOPT(*).
C     On some computer systems the user may need to open unit
C     2 before executing a call to SPLP( ).  This file is type
C     sequential and is unformatted.
C
C     Fortran unit=I1MACH(2) (check local setting) is used by SPLP( )
C     when the printed output feature (KEY=51) is used.  Normally
C     this feature is not used.  It is activated by input in the
C     options array PRGOPT(*).  For many computer systems I1MACH(2)=6.
C
C  |-------|
C  |Output:|
C  |-------|
C
C     INFO,PRIMAL(*),DUALS(*)
C     -----------------------
C     The integer flag INFO indicates why SPLP( ) has returned to the
C     user.  If INFO=1 the solution has been computed.  In this case
C     X(J)=PRIMAL(J) and W(I)=PRIMAL(I+NVARS).  The dual variables
C     for the equations A*x=w are in the array DUALS(I)=dual for
C     equation number I.  The dual value for the component X(J) that
C     has an upper or lower bound (or both) is returned in
C     DUALS(J+MRELAS).  The only other values for INFO are .LT. 0.
C     The meaning of these values can be found by reading
C     the diagnostic message in the output file, or by looking for
C     error number = (-INFO) under the heading "List of SPLP( ) Error
C     and Diagnostic Messages."
C     The diagnostic messages are printed using the error processing
C     subprogram XERMSG( ) with error category LEVEL=1.
C     See the document "Brief Instr. for Using the Sandia Math.
C     Subroutine Library," SAND79-2382, Nov., 1980, for further inform-
C     ation about resetting the usual response to a diagnostic message.
C
C     BL(*),BU(*),IND(*)
C     ------------------
C     These arrays are output parameters only under the (unusual)
C     circumstances where the stated problem is infeasible, has an
C     unbounded optimum value, or both.  These respective conditions
C     correspond to INFO=-1,-2 or -3.  For INFO=-1 or -3 certain comp-
C     onents of the vectors x or w will not satisfy the input bounds.
C     If component J of X or component I of W does not satisfy its input
C     bound because of infeasibility, then IND(J)=-4 or IND(I+NVARS)=-4,
C     respectively.  For INFO=-2 or -3 certain
C     components of the vector x could not be used as basic variables
C     because the objective function would have become unbounded.
C     In particular if component J of x corresponds to such a variable,
C     then IND(J)=-3.  Further, if the input value of IND(J)
C                      =1, then BU(J)=BL(J);
C                      =2, then BL(J)=BU(J);
C                      =4, then BL(J)=0.,BU(J)=0.
C
C     (The J-th variable in x has been restricted to an appropriate
C     feasible value.)
C     The negative output value for IND(*) allows the user to identify
C     those constraints that are not satisfied or those variables that
C     would cause unbounded values of the objective function.  Note
C     that the absolute value of IND(*), together with BL(*) and BU(*),
C     are valid input to SPLP( ).  In the case of infeasibility the
C     sum of magnitudes of the infeasible values is minimized.  Thus
C     one could reenter SPLP( ) with these components of x or w now
C     fixed at their present values.  This involves setting
C     the appropriate components of IND(*) = 3, and BL(*) = BU(*).
C
C     IBASIS(I),I=1,...,MRELAS
C     ------------------------
C     This array contains the indices of the variables that are
C     in the active basis set at the solution (INFO=1).  A value
C     of IBASIS(I) between 1 and NVARS corresponds to the variable
C     X(IBASIS(I)).  A value of IBASIS(I) between NVARS+1 and NVARS+
C     MRELAS corresponds to the variable W(IBASIS(I)-NVARS).
C
C     Computing with the Matrix A after Calling SPLP( )
C     -------------------------------------------------
C     Following the return from SPLP( ), nonzero entries of the MRELAS
C     by NVARS matrix A are available for usage by the user.  The method
C     for obtaining the next nonzero in column J with a row index
C     strictly greater than I in value, is completed by executing
C
C         CALL PNNZRS(I,AIJ,IPLACE,WORK,IWORK,J)
C
C     The value of I is also an output parameter.  If I.LE.0 on output,
C     then there are no more nonzeroes in column J.  If I.GT.0, the
C     output value for component number I of column J is in AIJ.  The
C     parameters WORK(*) and IWORK(*) are the same arguments as in the
C     call to SPLP( ).  The parameter IPLACE is a single INTEGER
C     working variable.
C
C     The data structure used for storage of the matrix A within SPLP( )
C     corresponds to sequential storage by columns as defined in
C     SAND78-0785.  Note that the names of the subprograms LNNZRS(),
C     LCHNGS(),LINITM(),LLOC(),LRWPGE(), and LRWVIR() have been
C     changed to PNNZRS(),PCHNGS(),PINITM(),IPLOC(),PRWPGE(), and
C     PRWVIR() respectively.  The error processing subprogram LERROR()
C     is no longer used; XERMSG() is used instead.
C
C  |-------------------------------|
C  |Subprograms Required by SPLP( )|
C  |-------------------------------|
C     Called by SPLP() are SPLPMN(),SPLPUP(),SPINIT(),SPOPT(),
C                          SPLPDM(),SPLPCE(),SPINCW(),SPLPFL(),
C                          SPLPFE(),SPLPMU().
C
C     Error Processing Subprograms XERMSG(),I1MACH(),R1MACH()
C
C     Sparse Matrix Subprograms PNNZRS(),PCHNGS(),PRWPGE(),PRWVIR(),
C                               PINITM(),IPLOC()
C
C     Mass Storage File Subprograms SOPENM(),SCLOSM(),SREADP(),SWRITP()
C
C     Basic Linear Algebra Subprograms SCOPY(),SASUM(),SDOT()
C
C     Sparse Matrix Basis Handling Subprograms LA05AS(),LA05BS(),
C                                             LA05CS(),LA05ED(),MC20AS()
C
C     Vector Output Subprograms SVOUT(),IVOUT()
C
C     Machine-sensitive Subprograms I1MACH( ),R1MACH( ),
C                                SOPENM(),SCLOSM(),SREADP(),SWRITP().
C     COMMON Block Used
C     -----------------
C     /LA05DS/ SMALL,LP,LENL,LENU,NCP,LROW,LCOL
C     See the document AERE-R8269 for further details.
C    |------------------------|
C    |Example of SPLP( ) Usage|
C    |------------------------|
C     PROGRAM LPEX
C     THE OPTIMIZATION PROBLEM IS TO FIND X1, X2, X3 THAT
C     MINIMIZE X1 + X2 + X3, X1.GE.0, X2.GE.0, X3 UNCONSTRAINED.
C
C     THE UNKNOWNS X1,X2,X3 ARE TO SATISFY CONSTRAINTS
C
C        X1 -3*X2 +4*X3 = 5
C        X1 -2*X2     .LE.3
C            2*X2 - X3.GE.4
C
C     WE FIRST DEFINE THE DEPENDENT VARIABLES
C          W1=X1 -3*X2 +4*X3
C          W2=X1- 2*X2
C          W3=    2*X2 -X3
C
C     WE NOW SHOW HOW TO USE SPLP( ) TO SOLVE THIS LINEAR OPTIMIZATION
C     PROBLEM.  EACH REQUIRED STEP WILL BE SHOWN IN THIS EXAMPLE.
C     DIMENSION COSTS(03),PRGOPT(01),DATTRV(18),BL(06),BU(06),IND(06),
C    *PRIMAL(06),DUALS(06),IBASIS(06),WORK(079),IWORK(103)
C
C     EXTERNAL USRMAT
C     MRELAS=3
C     NVARS=3
C
C     DEFINE THE ARRAY COSTS(*) FOR THE OBJECTIVE FUNCTION.
C     COSTS(01)=1.
C     COSTS(02)=1.
C     COSTS(03)=1.
C
C     PLACE THE NONZERO INFORMATION ABOUT THE MATRIX IN DATTRV(*).
C     DEFINE COL. 1:
C     DATTRV(01)=-1
C     DATTRV(02)=1
C     DATTRV(03)=1.
C     DATTRV(04)=2
C     DATTRV(05)=1.
C
C     DEFINE COL. 2:
C     DATTRV(06)=-2
C     DATTRV(07)=1
C     DATTRV(08)=-3.
C     DATTRV(09)=2
C     DATTRV(10)=-2.
C     DATTRV(11)=3
C     DATTRV(12)=2.
C
C     DEFINE COL. 3:
C     DATTRV(13)=-3
C     DATTRV(14)=1
C     DATTRV(15)=4.
C     DATTRV(16)=3
C     DATTRV(17)=-1.
C
C     DATTRV(18)=0
C
C     CONSTRAIN X1,X2 TO BE NONNEGATIVE. LET X3 HAVE NO BOUNDS.
C     BL(1)=0.
C     IND(1)=1
C     BL(2)=0.
C     IND(2)=1
C     IND(3)=4
C
C     CONSTRAIN W1=5,W2.LE.3, AND W3.GE.4.
C     BL(4)=5.
C     BU(4)=5.
C     IND(4)=3
C     BU(5)=3.
C     IND(5)=2
C     BL(6)=4.
C     IND(6)=1
C
C     INDICATE THAT NO MODIFICATIONS TO OPTIONS ARE IN USE.
C     PRGOPT(01)=1
C
C     DEFINE THE WORKING ARRAY LENGTHS.
C     LW=079
C     LIW=103
C     CALL SPLP(USRMAT,MRELAS,NVARS,COSTS,PRGOPT,DATTRV,
C    *BL,BU,IND,INFO,PRIMAL,DUALS,IBASIS,WORK,LW,IWORK,LIW)
C
C     CALCULATE VAL, THE MINIMAL VALUE OF THE OBJECTIVE FUNCTION.
C     VAL=SDOT(NVARS,COSTS,1,PRIMAL,1)
C
C     STOP
C     END
C    |------------------------|
C    |End of Example of Usage |
C    |------------------------|
C
C    |------------------------------------|
C    |Usage of SPLP( ) Subprogram Options.|
C    |------------------------------------|
C
C     Users frequently have a large variety of requirements for linear
C     optimization software.  Allowing for these varied requirements
C     is at cross purposes with the desire to keep the usage of SPLP( )
C     as simple as possible. One solution to this dilemma is as follows.
C     (1) Provide a version of SPLP( ) that solves a wide class of
C     problems and is easy to use. (2) Identify parameters within SPLP()
C     that certain users may want to change.  (3) Provide a means
C     of changing any selected number of these parameters that does
C     not require changing all of them.
C
C     Changing selected parameters is done by requiring
C     that the user provide an option array, PRGOPT(*), to SPLP( ).
C     The contents of PRGOPT(*) inform SPLP( ) of just those options
C     that are going to be modified within the total set of possible
C     parameters that can be modified.  The array PRGOPT(*) is a linked
C     list consisting of groups of data of the following form
C
C          LINK
C          KEY
C          SWITCH
C          data set
C
C     that describe the desired options.  The parameters LINK, KEY and
C     switch are each one word and are always required.  The data set
C     can be comprised of several words or can be empty.  The number of
C     words in the data set for each option depends on the value of
C     the parameter KEY.
C
C     The value of LINK points to the first entry of the next group
C     of data within PRGOPT(*).  The exception is when there are no more
C     options to change.  In that case, LINK=1 and the values for KEY,
C     SWITCH and data set are not referenced.  The general layout of
C     PRGOPT(*) is as follows:
C          ...PRGOPT(1)=LINK1 (link to first entry of next group)
C          .  PRGOPT(2)=KEY1 (KEY to the option change)
C          .  PRGOPT(3)=SWITCH1 (on/off switch for the option)
C          .  PRGOPT(4)=data value
C          .       .
C          .       .
C          .       .
C          ...PRGOPT(LINK1)=LINK2 (link to first entry of next group)
C          .  PRGOPT(LINK1+1)=KEY2 (KEY to option change)
C          .  PRGOPT(LINK1+2)=SWITCH2 (on/off switch for the option)
C          .  PRGOPT(LINK1+3)=data value
C          ...     .
C          .       .
C          .       .
C          ...PRGOPT(LINK)=1 (no more options to change)
C
C     A value of LINK that is .LE.0 or .GT. 10000 is an error.
C     In this case SPLP( ) returns with an error message, INFO=-14.
C     This helps prevent using invalid but positive values of LINK that
C     will probably extend beyond the program limits of PRGOPT(*).
C     Unrecognized values of KEY are ignored.  If the value of SWITCH is
C     zero then the option is turned off.  For any other value of SWITCH
C     the option is turned on.  This is used to allow easy changing of
C     options without rewriting PRGOPT(*).  The order of the options is
C     arbitrary and any number of options can be changed with the
C     following restriction.  To prevent cycling in processing of the
C     option array PRGOPT(*), a count of the number of options changed
C     is maintained.  Whenever this count exceeds 1000 an error message
C     (INFO=-15) is printed and the subprogram returns.
C
C     In the following description of the options, the value of
C     LATP indicates the amount of additional storage that a particular
C     option requires.  The sum of all of these values (plus one) is
C     the minimum dimension for the array PRGOPT(*).
C
C     If a user is satisfied with the nominal form of SPLP( ),
C     set PRGOPT(1)=1 (or PRGOPT(1)=1.E0).
C
C     Options:
C
C -----KEY = 50.  Change from a minimization problem to a maximization
C     problem.
C     If SWITCH=0  option is off; solve minimization problem.
C              =1  option is on; solve maximization problem.
C     data set =empty
C     LATP=3
C
C -----KEY = 51.  Change the amount of printed output.  The nominal form
C     of SPLP( ) has no printed output.
C     The first level of output (SWITCH=1) includes
C
C     (1) Minimum dimensions for the arrays COSTS(*),BL(*),BU(*),IND(*),
C         PRIMAL(*),DUALS(*),IBASIS(*), and PRGOPT(*).
C     (2) Problem dimensions MRELAS,NVARS.
C     (3) The types of and values for the bounds on x and w,
C         and the values of the components of the vector costs.
C     (4) Whether optimization problem is minimization or
C         maximization.
C     (5) Whether steepest edge or smallest reduced cost criteria used
C         for exchanging variables in the revised simplex method.
C
C     Whenever a solution has been found, (INFO=1),
C
C     (6) the value of the objective function,
C     (7) the values of the vectors x and w,
C     (8) the dual variables for the constraints A*x=w and the
C         bounded components of x,
C     (9) the indices of the basic variables,
C    (10) the number of revised simplex method iterations,
C    (11) the number of full decompositions of the basis matrix.
C
C     The second level of output (SWITCH=2) includes all for SWITCH=1
C     plus
C
C    (12) the iteration number,
C    (13) the column number to enter the basis,
C    (14) the column number to leave the basis,
C    (15) the length of the step taken.
C
C     The third level of output (SWITCH=3) includes all for SWITCH=2
C     plus
C    (16) critical quantities required in the revised simplex method.
C          This output is rather voluminous.  It is intended to be used
C          as a diagnostic tool in case of a failure in SPLP( ).
C
C     If SWITCH=0  option is off; no printed output.
C              =1  summary output.
C              =2  lots of output.
C              =3  even more output.
C     data set =empty
C     LATP=3
C
C -----KEY = 52.  Redefine the parameter, IDIGIT, which determines the
C     format and precision used for the printed output.  In the printed
C     output, at least ABS(IDIGIT) decimal digits per number is printed.
C     If IDIGIT.LT.0, 72 printing columns are used.  IF IDIGIT.GT.0, 133
C     printing columns are used.
C     If SWITCH=0  option is off; IDIGIT=-4.
C              =1  option is on.
C     data set =IDIGIT
C     LATP=4
C
C -----KEY = 53.  Redefine LAMAT and LBM, the lengths of the portions of
C     WORK(*) and IWORK(*) that are allocated to the sparse matrix
C     storage and the sparse linear equation solver, respectively.
C     LAMAT must be .GE. NVARS+7 and LBM must be positive.
C     If SWITCH=0  option is off; LAMAT=4*NVARS+7
C                                 LBM  =8*MRELAS.
C              =1  option is on.
C     data set =LAMAT
C               LBM
C     LATP=5
C
C -----KEY = 54. Redefine IPAGEF, the file number where the pages of the
C     sparse data matrix are stored.  IPAGEF must be positive and
C     different from ISAVE (see option 56).
C     If SWITCH=0  option is off; IPAGEF=1.
C              =1  option is on.
C     data set =IPAGEF
C     LATP=4
C
C -----KEY = 55.  Partial results have been computed and stored on unit
C     number ISAVE (see option 56), during a previous run of
C     SPLP( ).  This is a continuation from these partial results.
C     The arrays COSTS(*),BL(*),BU(*),IND(*) do not have to have
C     the same values as they did when the checkpointing occurred.
C     This feature makes it possible for the user to do certain
C     types of parameter studies such as changing costs and varying
C     the constraints of the problem.  This file is rewound both be-
C     fore and after reading the partial results.
C     If SWITCH=0  option is off; start a new problem.
C              =1  option is on; continue from partial results
C                  that are stored in file ISAVE.
C     data set = empty
C     LATP=3
C
C -----KEY = 56.  Redefine ISAVE, the file number where the partial
C     results are stored (see option 57).  ISAVE must be positive and
C     different from IPAGEF (see option 54).
C     If SWITCH=0  option is off; ISAVE=2.
C              =1  option is on.
C     data set =ISAVE
C     LATP=4
C
C -----KEY = 57.  Save the partial results after maximum number of
C     iterations, MAXITR, or at the optimum.  When this option is on,
C     data essential to continuing the calculation is saved on a file
C     using a Fortran binary write operation.  The data saved includes
C     all the information about the sparse data matrix A.  Also saved
C     is information about the current basis.  Nominally the partial
C     results are saved on Fortran unit 2.  This unit number can be
C     redefined (see option 56).  If the save option is on,
C     this file must be opened (or declared) by the user prior to the
C     call to SPLP( ).  A crude upper bound for the number of words
C     written to this file is 6*nz.  Here nz= number of nonzeros in A.
C     If SWITCH=0  option is off; do not save partial results.
C              =1  option is on; save partial results.
C     data set = empty
C     LATP=3
C
C -----KEY = 58.  Redefine the maximum number of iterations, MAXITR, to
C     be taken before returning to the user.
C     If SWITCH=0  option is off; MAXITR=3*(NVARS+MRELAS).
C              =1  option is on.
C     data set =MAXITR
C     LATP=4
C
C -----KEY = 59.  Provide SPLP( ) with exactly MRELAS indices which
C     comprise a feasible, nonsingular basis.  The basis must define a
C     feasible point: values for x and w such that A*x=w and all the
C     stated bounds on x and w are satisfied.  The basis must also be
C     nonsingular.  The failure of either condition will cause an error
C     message (INFO=-23 or =-24, respectively).  Normally, SPLP( ) uses
C     identity matrix columns which correspond to the components of w.
C     This option would normally not be used when restarting from
C     a previously saved run (KEY=57).
C     In numbering the unknowns,
C     the components of x are numbered (1-NVARS) and the components
C     of w are numbered (NVARS+1)-(NVARS+MRELAS).  A value for an
C     index .LE. 0 or .GT. (NVARS+MRELAS) is an error (INFO=-16).
C     If SWITCH=0  option is off; SPLP( ) chooses the initial basis.
C              =1  option is on; user provides the initial basis.
C     data set =MRELAS indices of basis; order is arbitrary.
C     LATP=MRELAS+3
C
C -----KEY = 60.  Provide the scale factors for the columns of the data
C     matrix A.  Normally, SPLP( ) computes the scale factors as the
C     reciprocals of the max. norm of each column.
C     If SWITCH=0  option is off; SPLP( ) computes the scale factors.
C              =1  option is on; user provides the scale factors.
C     data set =scaling for column J, J=1,NVARS; order is sequential.
C     LATP=NVARS+3
C
C -----KEY = 61.  Provide a scale factor, COSTSC, for the vector of
C     costs.  Normally, SPLP( ) computes this scale factor to be the
C     reciprocal of the max. norm of the vector costs after the column
C     scaling has been applied.
C     If SWITCH=0  option is off; SPLP( ) computes COSTSC.
C              =1  option is on; user provides COSTSC.
C     data set =COSTSC
C     LATP=4
C
C -----KEY = 62.  Provide size parameters, ASMALL and ABIG, the smallest
C     and largest magnitudes of nonzero entries in the data matrix A,
C     respectively.  When this option is on, SPLP( ) will check the
C     nonzero entries of A to see if they are in the range of ASMALL and
C     ABIG.  If an entry of A is not within this range, SPLP( ) returns
C     an error message, INFO=-22. Both ASMALL and ABIG must be positive
C     with ASMALL .LE. ABIG.  Otherwise,  an error message is returned,
C     INFO=-17.
C     If SWITCH=0  option is off; no checking of the data matrix is done
C              =1  option is on; checking is done.
C     data set =ASMALL
C               ABIG
C     LATP=5
C
C -----KEY = 63.  Redefine the relative tolerance, TOLLS, used in
C     checking if the residuals are feasible.  Normally,
C     TOLLS=RELPR, where RELPR is the machine precision.
C     If SWITCH=0  option is off; TOLLS=RELPR.
C              =1  option is on.
C     data set =TOLLS
C     LATP=4
C
C -----KEY = 64. Use the minimum reduced cost pricing strategy to choose
C     columns to enter the basis.  Normally, SPLP( ) uses the steepest
C     edge pricing strategy which is the best local move.  The steepest
C     edge pricing strategy generally uses fewer iterations than the
C     minimum reduced cost pricing, but each iteration costs more in the
C     number of calculations done.  The steepest edge pricing is
C     considered to be more efficient.  However, this is very problem
C     dependent.  That is why SPLP( ) provides the option of either
C     pricing strategy.
C     If SWITCH=0  option is off; steepest option edge pricing is used.
C              =1  option is on; minimum reduced cost pricing is used.
C     data set =empty
C     LATP=3
C
C -----KEY = 65.  Redefine MXITBR, the number of iterations between
C     recalculating the error in the primal solution.  Normally, MXITBR
C     is set to 10.  The error in the primal solution is used to monitor
C     the error in solving the linear system.  This is an expensive
C     calculation and every tenth iteration is generally often enough.
C     If SWITCH=0  option is off; MXITBR=10.
C              =1  option is on.
C     data set =MXITBR
C     LATP=4
C
C -----KEY = 66.  Redefine NPP, the number of negative reduced costs
C     (at most) to be found at each iteration of choosing
C     a variable to enter the basis.  Normally NPP is set
C     to NVARS which implies that all of the reduced costs
C     are computed at each such step.  This "partial
C     pricing" may very well increase the total number
C     of iterations required.  However it decreases the
C     number of calculations at each iteration.
C     therefore the effect on overall efficiency is quite
C     problem-dependent.
C
C     if SWITCH=0 option is off; NPP=NVARS
C              =1 option is on.
C     data set =NPP
C     LATP=4
C
C -----KEY =  67.  Redefine the tuning factor (PHI) used to scale the
C     error estimates for the primal and dual linear algebraic systems
C     of equations.  Normally, PHI = 1.E0, but in some environments it
C     may be necessary to reset PHI to the range 0.001-0.01.  This is
C     particularly important for machines with short word lengths.
C
C     if SWITCH = 0 option is off; PHI=1.E0.
C               = 1 option is on.
C     Data Set  = PHI
C     LATP=4
C
C -----KEY = 68.  Used together with the subprogram FULMAT(), provided
C     with the SPLP() package, for passing a standard Fortran two-
C     dimensional array containing the constraint matrix.  Thus the sub-
C     program FULMAT must be declared in a Fortran EXTERNAL statement.
C     The two-dimensional array is passed as the argument DATTRV.
C     The information about the array and problem dimensions are passed
C     in the option array PRGOPT(*).  It is an error if FULMAT() is
C     used and this information is not passed in PRGOPT(*).
C
C     if SWITCH = 0 option is off; this is an error is FULMAT() is
C                                  used.
C               = 1 option is on.
C     Data Set  = IA = row dimension of two-dimensional array.
C                 MRELAS = number of constraint equations.
C                 NVARS  = number of dependent variables.
C     LATP = 6
C -----KEY = 69.  Normally a relative tolerance (TOLLS, see option 63)
C     is used to decide if the problem is feasible.  If this test fails
C     an absolute test will be applied using the value TOLABS.
C     Nominally TOLABS = zero.
C     If SWITCH = 0 option is off; TOLABS = zero.
C               = 1 option is on.
C     Data set  = TOLABS
C     LATP = 4
C
C    |-----------------------------|
C    |Example of Option array Usage|
C    |-----------------------------|
C     To illustrate the usage of the option array, let us suppose that
C     the user has the following nonstandard requirements:
C
C          a) Wants to change from minimization to maximization problem.
C          b) Wants to limit the number of simplex steps to 100.
C          c) Wants to save the partial results after 100 steps on
C             Fortran unit 2.
C
C     After these 100 steps are completed the user wants to continue the
C     problem (until completed) using the partial results saved on
C     Fortran unit 2.  Here are the entries of the array PRGOPT(*)
C     that accomplish these tasks.  (The definitions of the other
C     required input parameters are not shown.)
C
C     CHANGE TO A MAXIMIZATION PROBLEM; KEY=50.
C     PRGOPT(01)=4
C     PRGOPT(02)=50
C     PRGOPT(03)=1
C
C     LIMIT THE NUMBER OF SIMPLEX STEPS TO 100; KEY=58.
C     PRGOPT(04)=8
C     PRGOPT(05)=58
C     PRGOPT(06)=1
C     PRGOPT(07)=100
C
C     SAVE THE PARTIAL RESULTS, AFTER 100 STEPS, ON FORTRAN
C     UNIT 2; KEY=57.
C     PRGOPT(08)=11
C     PRGOPT(09)=57
C     PRGOPT(10)=1
C
C     NO MORE OPTIONS TO CHANGE.
C     PRGOPT(11)=1
C     The user makes the CALL statement for SPLP( ) at this point.
C     Now to restart, using the partial results after 100 steps, define
C     new values for the array PRGOPT(*):
C
C     AGAIN INFORM SPLP( ) THAT THIS IS A MAXIMIZATION PROBLEM.
C     PRGOPT(01)=4
C     PRGOPT(02)=50
C     PRGOPT(03)=1
C
C     RESTART, USING SAVED PARTIAL RESULTS; KEY=55.
C     PRGOPT(04)=7
C     PRGOPT(05)=55
C     PRGOPT(06)=1
C
C     NO MORE OPTIONS TO CHANGE.  THE SUBPROGRAM SPLP( ) IS NO LONGER
C     LIMITED TO 100 SIMPLEX STEPS BUT WILL RUN UNTIL COMPLETION OR
C     MAX.=3*(MRELAS+NVARS) ITERATIONS.
C     PRGOPT(07)=1
C     The user now makes a CALL to subprogram SPLP( ) to compute the
C     solution.
C    |-------------------------------------------|
C    |End of Usage of SPLP( ) Subprogram Options.|
C    |-------------------------------------------|
C
C     |----------------------------------------------|
C     |List of SPLP( ) Error and Diagnostic Messages.|
C     |----------------------------------------------|
C      This section may be required to understand the meanings of the
C     error flag =-INFO  that may be returned from SPLP( ).
C
C -----1. There is no set of values for x and w that satisfy A*x=w and
C     the stated bounds.  The problem can be made feasible by ident-
C     ifying components of w that are now infeasible and then rede-
C     signating them as free variables.  Subprogram SPLP( ) only
C     identifies an infeasible problem; it takes no other action to
C     change this condition.  Message:
C     SPLP( ). THE PROBLEM APPEARS TO BE INFEASIBLE.
C     ERROR NUMBER =         1
C
C     2. One of the variables in either the vector x or w was con-
C     strained at a bound.  Otherwise the objective function value,
C     (transpose of costs)*x, would not have a finite optimum.
C     Message:
C     SPLP( ). THE PROBLEM APPEARS TO HAVE NO FINITE SOLN.
C     ERROR NUMBER =         2
C
C     3.  Both of the conditions of 1. and 2. above have occurred.
C     Message:
C     SPLP( ). THE PROBLEM APPEARS TO BE INFEASIBLE AND TO
C     HAVE NO FINITE SOLN.
C     ERROR NUMBER =         3
C
C -----4.  The REAL and INTEGER working arrays, WORK(*) and IWORK(*),
C     are not long enough. The values (I1) and (I2) in the message
C     below will give you the minimum length required.  Also redefine
C     LW and LIW, the lengths of these arrays.  Message:
C     SPLP( ). WORK OR IWORK IS NOT LONG ENOUGH. LW MUST BE (I1)
C     AND LIW MUST BE (I2).
C               IN ABOVE MESSAGE, I1=         0
C               IN ABOVE MESSAGE, I2=         0
C     ERROR NUMBER =        4
C
C -----5. and 6.  These error messages often mean that one or more
C     arguments were left out of the call statement to SPLP( ) or
C     that the values of MRELAS and NVARS have been over-written
C     by garbage.  Messages:
C     SPLP( ). VALUE OF MRELAS MUST BE .GT.0. NOW=(I1).
C               IN ABOVE MESSAGE, I1=         0
C     ERROR NUMBER =         5
C
C     SPLP( ). VALUE OF NVARS MUST BE .GT.0. NOW=(I1).
C               IN ABOVE MESSAGE, I1=         0
C     ERROR NUMBER =         6
C
C -----7.,8., and 9.  These error messages can occur as the data matrix
C     is being defined by either USRMAT( ) or the user-supplied sub-
C     program, 'NAME'( ). They would indicate a mistake in the contents
C     of DATTRV(*), the user-written subprogram or that data has been
C     over-written.
C     Messages:
C     SPLP( ). MORE THAN 2*NVARS*MRELAS ITERS. DEFINING OR UPDATING
C     MATRIX DATA.
C     ERROR NUMBER =        7
C
C     SPLP( ). ROW INDEX (I1) OR COLUMN INDEX (I2) IS OUT OF RANGE.
C               IN ABOVE MESSAGE, I1=         1
C               IN ABOVE MESSAGE, I2=        12
C     ERROR NUMBER =        8
C
C     SPLP( ). INDICATION FLAG (I1) FOR MATRIX DATA MUST BE
C     EITHER 0 OR 1.
C               IN ABOVE MESSAGE, I1=        12
C     ERROR NUMBER =        9
C
C -----10. and 11.  The type of bound (even no bound) and the bounds
C     must be specified for each independent variable. If an independent
C     variable has both an upper and lower bound, the bounds must be
C     consistent.  The lower bound must be .LE. the upper bound.
C     Messages:
C     SPLP( ). INDEPENDENT VARIABLE (I1) IS NOT DEFINED.
C               IN ABOVE MESSAGE, I1=         1
C     ERROR NUMBER =        10
C
C     SPLP( ).  LOWER BOUND (R1) AND UPPER BOUND (R2) FOR INDEP.
C     VARIABLE (I1) ARE NOT CONSISTENT.
C               IN ABOVE MESSAGE, I1=         1
C               IN ABOVE MESSAGE, R1=    0.
C               IN ABOVE MESSAGE, R2=    -.1000000000E+01
C     ERROR NUMBER =        11
C
C -----12. and 13.  The type of bound (even no bound) and the bounds
C     must be specified for each dependent variable.  If a dependent
C     variable has both an upper and lower bound, the bounds must be
C     consistent. The lower bound must be .LE. the upper bound.
C     Messages:
C     SPLP( ). DEPENDENT VARIABLE (I1) IS NOT DEFINED.
C               IN ABOVE MESSAGE, I1=         1
C     ERROR NUMBER =        12
C
C     SPLP( ).  LOWER BOUND (R1) AND UPPER BOUND (R2) FOR DEP.
C      VARIABLE (I1) ARE NOT CONSISTENT.
C               IN ABOVE MESSAGE, I1=         1
C               IN ABOVE MESSAGE, R1=    0.
C               IN ABOVE MESSAGE, R2=    -.1000000000E+01
C     ERROR NUMBER =        13
C
C -----14. - 21.  These error messages can occur when processing the
C     option array, PRGOPT(*), supplied by the user.  They would
C     indicate a mistake in defining PRGOPT(*) or that data has been
C     over-written.  See heading Usage of SPLP( )
C     Subprogram Options, for details on how to define PRGOPT(*).
C     Messages:
C     SPLP( ). THE USER OPTION ARRAY HAS UNDEFINED DATA.
C     ERROR NUMBER =        14
C
C     SPLP( ). OPTION ARRAY PROCESSING IS CYCLING.
C     ERROR NUMBER =        15
C
C     SPLP( ). AN INDEX OF USER-SUPPLIED BASIS IS OUT OF RANGE.
C     ERROR NUMBER =        16
C
C     SPLP( ). SIZE PARAMETERS FOR MATRIX MUST BE SMALLEST AND LARGEST
C     MAGNITUDES OF NONZERO ENTRIES.
C     ERROR NUMBER =        17
C
C     SPLP( ). THE NUMBER OF REVISED SIMPLEX STEPS BETWEEN CHECK-POINTS
C     MUST BE POSITIVE.
C     ERROR NUMBER =        18
C
C     SPLP( ). FILE NUMBERS FOR SAVED DATA AND MATRIX PAGES MUST BE
C     POSITIVE AND NOT EQUAL.
C     ERROR NUMBER =        19
C
C     SPLP( ). USER-DEFINED VALUE OF LAMAT (I1)
C     MUST BE .GE. NVARS+7.
C               IN ABOVE MESSAGE, I1=         1
C     ERROR NUMBER =         20
C
C     SPLP( ). USER-DEFINED VALUE OF LBM MUST BE .GE. 0.
C     ERROR NUMBER =         21
C
C -----22.  The user-option, number 62, to check the size of the matrix
C     data has been used.  An element of the matrix does not lie within
C     the range of ASMALL and ABIG, parameters provided by the user.
C     (See the heading: Usage of SPLP( ) Subprogram Options,
C     for details about this feature.)  Message:
C     SPLP( ). A MATRIX ELEMENT'S SIZE IS OUT OF THE SPECIFIED RANGE.
C     ERROR NUMBER =        22
C
C -----23.  The user has provided an initial basis that is singular.
C     In this case, the user can remedy this problem by letting
C     subprogram SPLP( ) choose its own initial basis.  Message:
C     SPLP( ). A SINGULAR INITIAL BASIS WAS ENCOUNTERED.
C     ERROR NUMBER =         23
C
C -----24.  The user has provided an initial basis which is infeasible.
C     The x and w values it defines do not satisfy A*x=w and the stated
C     bounds.  In this case, the user can let subprogram SPLP( )
C     choose its own initial basis.  Message:
C     SPLP( ). AN INFEASIBLE INITIAL BASIS WAS ENCOUNTERED.
C     ERROR NUMBER =        24
C
C -----25. Subprogram SPLP( ) has completed the maximum specified number
C     of iterations.  (The nominal maximum number is 3*(MRELAS+NVARS).)
C     The results, necessary to continue on from
C     this point, can be saved on Fortran unit 2 by activating option
C     KEY=57.  If the user anticipates continuing the calculation, then
C     the contents of Fortran unit 2 must be retained intact.  This
C     is not done by subprogram SPLP( ), so the user needs to save unit
C     2 by using the appropriate system commands.  Message:
C     SPLP( ). MAX. ITERS. (I1) TAKEN. UP-TO-DATE RESULTS
C     SAVED ON FILE (I2). IF(I2)=0, NO SAVE.
C               IN ABOVE MESSAGE, I1=       500
C               IN ABOVE MESSAGE, I2=         2
C     ERROR NUMBER =        25
C
C -----26.  This error should never happen.  Message:
C     SPLP( ). MOVED TO A SINGULAR POINT. THIS SHOULD NOT HAPPEN.
C     ERROR NUMBER =        26
C
C -----27.  The subprogram LA05A( ), which decomposes the basis matrix,
C     has returned with an error flag (R1).  (See the document,
C     "Fortran subprograms for handling sparse linear programming
C     bases", AERE-R8269, J.K. Reid, Jan., 1976, H.M. Stationery Office,
C     for an explanation of this error.)  Message:
C     SPLP( ). LA05A( ) RETURNED ERROR FLAG (R1) BELOW.
C               IN ABOVE MESSAGE, R1=    -.5000000000E+01
C     ERROR NUMBER =        27
C
C -----28.  The sparse linear solver package, LA05*( ), requires more
C     space.  The value of LBM must be increased.  See the companion
C     document, Usage of SPLP( ) Subprogram Options, for details on how
C     to increase the value of LBM.  Message:
C     SPLP( ). SHORT ON STORAGE FOR LA05*( ) PACKAGE. USE PRGOPT(*)
C     TO GIVE MORE.
C     ERROR NUMBER =        28
C
C -----29.  The row dimension of the two-dimensional Fortran array,
C     the number of constraint equations (MRELAS), and the number
C     of variables (NVARS), were not passed to the subprogram
C     FULMAT().  See KEY = 68 for details.  Message:
C     FULMAT() OF SPLP() PACKAGE. ROW DIM., MRELAS, NVARS ARE
C     MISSING FROM PRGOPT(*).
C     ERROR NUMBER =        29
C
C     |------------------------------------------------------|
C     |End of List of SPLP( ) Error and Diagnostic Messages. |
C     |------------------------------------------------------|
C***REFERENCES  R. J. Hanson and K. L. Hiebert, A sparse linear
C                 programming subprogram, Report SAND81-0297, Sandia
C                 National Laboratories, 1981.
C***ROUTINES CALLED  SPLPMN, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   811215  DATE WRITTEN
C   890605  Corrected references to XERRWV.  (WRB)
C   890605  Removed unreferenced labels.  (WRB)
C   890605  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SPLP
      REAL             BL(*),BU(*),COSTS(*),DATTRV(*),DUALS(*),
     * PRGOPT(*),PRIMAL(*),WORK(*),ZERO
C
      INTEGER IBASIS(*),IND(*),IWORK(*)
      CHARACTER*8 XERN1, XERN2
C
      EXTERNAL USRMAT
C
C***FIRST EXECUTABLE STATEMENT  SPLP
      ZERO=0.E0
      IOPT=1
C
C     VERIFY THAT MRELAS, NVARS .GT. 0.
C
      IF (MRELAS.LE.0) THEN
         WRITE (XERN1, '(I8)') MRELAS
         CALL XERMSG ('SLATEC', 'SPLP', 'VALUE OF MRELAS MUST BE ' //
     *      '.GT. 0.  NOW = ' // XERN1, 5, 1)
         INFO = -5
         RETURN
      ENDIF
C
      IF (NVARS.LE.0) THEN
         WRITE (XERN1, '(I8)') NVARS
         CALL XERMSG ('SLATEC', 'SPLP', 'VALUE OF NVARS MUST BE ' //
     *      '.GT. 0.  NOW = ' // XERN1, 6, 1)
         INFO = -6
         RETURN
      ENDIF
C
      LMX=4*NVARS+7
      LBM=8*MRELAS
      LAST = 1
      IADBIG=10000
      ICTMAX=1000
      ICTOPT= 0
C
C     LOOK IN OPTION ARRAY FOR CHANGES TO WORK ARRAY LENGTHS.
20008 NEXT=PRGOPT(LAST)
      IF (.NOT.(NEXT.LE.0 .OR. NEXT.GT.IADBIG)) GO TO 20010
C
C     THE CHECKS FOR SMALL OR LARGE VALUES OF NEXT ARE TO PREVENT
C     WORKING WITH UNDEFINED DATA.
      NERR=14
      CALL XERMSG ('SLATEC', 'SPLP',
     +   'THE USER OPTION ARRAY HAS UNDEFINED DATA.', NERR, IOPT)
      INFO=-NERR
      RETURN
20010 IF (.NOT.(NEXT.EQ.1)) GO TO 10001
      GO TO 20009
10001 IF (.NOT.(ICTOPT.GT.ICTMAX)) GO TO 10002
      NERR=15
      CALL XERMSG ('SLATEC', 'SPLP',
     +   'OPTION ARRAY PROCESSING IS CYCLING.', NERR, IOPT)
      INFO=-NERR
      RETURN
10002 CONTINUE
      KEY = PRGOPT(LAST+1)
C
C     IF KEY = 53, USER MAY SPECIFY LENGTHS OF PORTIONS
C    OF WORK(*) AND IWORK(*) THAT ARE ALLOCATED TO THE
C     SPARSE MATRIX STORAGE AND SPARSE LINEAR EQUATION
C     SOLVING.
      IF (.NOT.(KEY.EQ.53)) GO TO 20013
      IF (.NOT.(PRGOPT(LAST+2).NE.ZERO)) GO TO 20016
      LMX=PRGOPT(LAST+3)
      LBM=PRGOPT(LAST+4)
20016 CONTINUE
20013 ICTOPT = ICTOPT+1
      LAST = NEXT
      GO TO 20008
C
C     CHECK LENGTH VALIDITY OF SPARSE MATRIX STAGING AREA.
C
20009 IF (LMX.LT.NVARS+7) THEN
         WRITE (XERN1, '(I8)') LMX
         CALL XERMSG ('SLATEC', 'SPLP', 'USER-DEFINED VALUE OF ' //
     *      'LAMAT = ' // XERN1 // ' MUST BE .GE. NVARS+7.', 20, 1)
         INFO = -20
         RETURN
      ENDIF
C
C     TRIVIAL CHECK ON LENGTH OF LA05*() MATRIX AREA.
      IF (.NOT.(LBM.LT.0)) GO TO 20022
      NERR=21
      CALL XERMSG ('SLATEC', 'SPLP',
     +   'USER-DEFINED VALUE OF LBM MUST BE .GE. 0.', NERR, IOPT)
      INFO=-NERR
      RETURN
20022 CONTINUE
C
C     DEFINE POINTERS FOR STARTS OF SUBARRAYS USED IN WORK(*)
C     AND IWORK(*) IN OTHER SUBPROGRAMS OF THE PACKAGE.
      LAMAT=1
      LCSC=LAMAT+LMX
      LCOLNR=LCSC+NVARS
      LERD=LCOLNR+NVARS
      LERP=LERD+MRELAS
      LBASMA=LERP+MRELAS
      LWR=LBASMA+LBM
      LRZ=LWR+MRELAS
      LRG=LRZ+NVARS+MRELAS
      LRPRIM=LRG+NVARS+MRELAS
      LRHS=LRPRIM+MRELAS
      LWW=LRHS+MRELAS
      LWORK=LWW+MRELAS-1
      LIMAT=1
      LIBB=LIMAT+LMX
      LIBRC=LIBB+NVARS+MRELAS
      LIPR=LIBRC+2*LBM
      LIWR=LIPR+2*MRELAS
      LIWORK=LIWR+8*MRELAS-1
C
C     CHECK ARRAY LENGTH VALIDITY OF WORK(*), IWORK(*).
C
      IF (LW.LT.LWORK .OR. LIW.LT.LIWORK) THEN
         WRITE (XERN1, '(I8)') LWORK
         WRITE (XERN2, '(I8)') LIWORK
         CALL XERMSG ('SLATEC', 'SPLP', 'WORK OR IWORK IS NOT LONG ' //
     *      'ENOUGH. LW MUST BE = ' // XERN1 // ' AND LIW MUST BE = ' //
     *      XERN2, 4, 1)
         INFO = -4
         RETURN
      ENDIF
C
      CALL SPLPMN(USRMAT,MRELAS,NVARS,COSTS,PRGOPT,DATTRV,
     * BL,BU,IND,INFO,PRIMAL,DUALS,WORK(LAMAT),
     * WORK(LCSC),WORK(LCOLNR),WORK(LERD),WORK(LERP),WORK(LBASMA),
     * WORK(LWR),WORK(LRZ),WORK(LRG),WORK(LRPRIM),WORK(LRHS),
     * WORK(LWW),LMX,LBM,IBASIS,IWORK(LIBB),IWORK(LIMAT),
     * IWORK(LIBRC),IWORK(LIPR),IWORK(LIWR))
C
C     CALL SPLPMN(USRMAT,MRELAS,NVARS,COSTS,PRGOPT,DATTRV,
C    1 BL,BU,IND,INFO,PRIMAL,DUALS,AMAT,
C    2 CSC,COLNRM,ERD,ERP,BASMAT,
C    3 WR,RZ,RG,RPRIM,RHS,
C    4 WW,LMX,LBM,IBASIS,IBB,IMAT,
C    5 IBRC,IPR,IWR)
C
      RETURN
      END
