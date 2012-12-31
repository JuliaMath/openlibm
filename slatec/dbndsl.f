*DECK DBNDSL
      SUBROUTINE DBNDSL (MODE, G, MDG, NB, IP, IR, X, N, RNORM)
C***BEGIN PROLOGUE  DBNDSL
C***PURPOSE  Solve the least squares problem for a banded matrix using
C            sequential accumulation of rows of the data matrix.
C            Exactly one right-hand side vector is permitted.
C***LIBRARY   SLATEC
C***CATEGORY  D9
C***TYPE      DOUBLE PRECISION (BNDSOL-S, DBNDSL-D)
C***KEYWORDS  BANDED MATRIX, CURVE FITTING, LEAST SQUARES
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C     These subroutines solve the least squares problem Ax = b for
C     banded matrices A using sequential accumulation of rows of the
C     data matrix.  Exactly one right-hand side vector is permitted.
C
C     These subroutines are intended for the type of least squares
C     systems that arise in applications such as curve or surface
C     fitting of data.  The least squares equations are accumulated and
C     processed using only part of the data.  This requires a certain
C     user interaction during the solution of Ax = b.
C
C     Specifically, suppose the data matrix (A B) is row partitioned
C     into Q submatrices.  Let (E F) be the T-th one of these
C     submatrices where E = (0 C 0).  Here the dimension of E is MT by N
C     and the dimension of C is MT by NB.  The value of NB is the
C     bandwidth of A.  The dimensions of the leading block of zeros in E
C     are MT by JT-1.
C
C     The user of the subroutine DBNDAC provides MT,JT,C and F for
C     T=1,...,Q.  Not all of this data must be supplied at once.
C
C     Following the processing of the various blocks (E F), the matrix
C     (A B) has been transformed to the form (R D) where R is upper
C     triangular and banded with bandwidth NB.  The least squares
C     system Rx = d is then easily solved using back substitution by
C     executing the statement CALL DBNDSL(1,...). The sequence of
C     values for JT must be nondecreasing.  This may require some
C     preliminary interchanges of rows and columns of the matrix A.
C
C     The primary reason for these subroutines is that the total
C     processing can take place in a working array of dimension MU by
C     NB+1.  An acceptable value for MU is
C
C                       MU = MAX(MT + N + 1),
C
C     where N is the number of unknowns.
C
C     Here the maximum is taken over all values of MT for T=1,...,Q.
C     Notice that MT can be taken to be a small as one, showing that
C     MU can be as small as N+2.  The subprogram DBNDAC processes the
C     rows more efficiently if MU is large enough so that each new
C     block (C F) has a distinct value of JT.
C
C     The four principle parts of these algorithms are obtained by the
C     following call statements
C
C     CALL DBNDAC(...)  Introduce new blocks of data.
C
C     CALL DBNDSL(1,...)Compute solution vector and length of
C                       residual vector.
C
C     CALL DBNDSL(2,...)Given any row vector H solve YR = H for the
C                       row vector Y.
C
C     CALL DBNDSL(3,...)Given any column vector W solve RZ = W for
C                       the column vector Z.
C
C     The dots in the above call statements indicate additional
C     arguments that will be specified in the following paragraphs.
C
C     The user must dimension the array appearing in the call list..
C     G(MDG,NB+1)
C
C     Description of calling sequence for DBNDAC..
C
C     The entire set of parameters for DBNDAC are
C
C     Input.. All Type REAL variables are DOUBLE PRECISION
C
C     G(*,*)            The working array into which the user will
C                       place the MT by NB+1 block (C F) in rows IR
C                       through IR+MT-1, columns 1 through NB+1.
C                       See descriptions of IR and MT below.
C
C     MDG               The number of rows in the working array
C                       G(*,*).  The value of MDG should be .GE. MU.
C                       The value of MU is defined in the abstract
C                       of these subprograms.
C
C     NB                The bandwidth of the data matrix A.
C
C     IP                Set by the user to the value 1 before the
C                       first call to DBNDAC.  Its subsequent value
C                       is controlled by DBNDAC to set up for the
C                       next call to DBNDAC.
C
C     IR                Index of the row of G(*,*) where the user is
C                       the user to the value 1 before the first call
C                       to DBNDAC.  Its subsequent value is controlled
C                       by DBNDAC. A value of IR .GT. MDG is considered
C                       an error.
C
C     MT,JT             Set by the user to indicate respectively the
C                       number of new rows of data in the block and
C                       the index of the first nonzero column in that
C                       set of rows (E F) = (0 C 0 F) being processed.
C     Output.. All Type REAL variables are DOUBLE PRECISION
C
C     G(*,*)            The working array which will contain the
C                       processed rows of that part of the data
C                       matrix which has been passed to DBNDAC.
C
C     IP,IR             The values of these arguments are advanced by
C                       DBNDAC to be ready for storing and processing
C                       a new block of data in G(*,*).
C
C     Description of calling sequence for DBNDSL..
C
C     The user must dimension the arrays appearing in the call list..
C
C     G(MDG,NB+1), X(N)
C
C     The entire set of parameters for DBNDSL are
C
C     Input..
C
C     MODE              Set by the user to one of the values 1, 2, or
C                       3.  These values respectively indicate that
C                       the solution of AX = B, YR = H or RZ = W is
C                       required.
C
C     G(*,*),MDG,       These arguments all have the same meaning and
C      NB,IP,IR         contents as following the last call to DBNDAC.
C
C     X(*)              With mode=2 or 3 this array contains,
C                       respectively, the right-side vectors H or W of
C                       the systems YR = H or RZ = W.
C
C     N                 The number of variables in the solution
C                       vector.  If any of the N diagonal terms are
C                       zero the subroutine DBNDSL prints an
C                       appropriate message.  This condition is
C                       considered an error.
C
C     Output..
C
C     X(*)              This array contains the solution vectors X,
C                       Y or Z of the systems AX = B, YR = H or
C                       RZ = W depending on the value of MODE=1,
C                       2 or 3.
C
C     RNORM             If MODE=1 RNORM is the Euclidean length of the
C                       residual vector AX-B.  When MODE=2 or 3 RNORM
C                       is set to zero.
C
C     Remarks..
C
C     To obtain the upper triangular matrix and transformed right-hand
C     side vector D so that the super diagonals of R form the columns
C     of G(*,*), execute the following Fortran statements.
C
C     NBP1=NB+1
C
C     DO 10 J=1, NBP1
C
C  10 G(IR,J) = 0.E0
C
C     MT=1
C
C     JT=N+1
C
C     CALL DBNDAC(G,MDG,NB,IP,IR,MT,JT)
C
C***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares
C                 Problems, Prentice-Hall, Inc., 1974, Chapter 27.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DBNDSL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(MDG,*),X(*)
C***FIRST EXECUTABLE STATEMENT  DBNDSL
      ZERO=0.D0
C
      RNORM=ZERO
      GO TO (10,90,50), MODE
C                                   ********************* MODE = 1
C                                   ALG. STEP 26
   10      DO 20 J=1,N
           X(J)=G(J,NB+1)
   20 CONTINUE
      RSQ=ZERO
      NP1=N+1
      IRM1=IR-1
      IF (NP1.GT.IRM1) GO TO 40
           DO 30 J=NP1,IRM1
           RSQ=RSQ+G(J,NB+1)**2
   30 CONTINUE
      RNORM=SQRT(RSQ)
   40 CONTINUE
C                                   ********************* MODE = 3
C                                   ALG. STEP 27
   50      DO 80 II=1,N
           I=N+1-II
C                                   ALG. STEP 28
           S=ZERO
           L=MAX(0,I-IP)
C                                   ALG. STEP 29
           IF (I.EQ.N) GO TO 70
C                                   ALG. STEP 30
           IE=MIN(N+1-I,NB)
                DO 60 J=2,IE
                JG=J+L
                IX=I-1+J
                S=S+G(I,JG)*X(IX)
   60 CONTINUE
C                                   ALG. STEP 31
   70      IF (G(I,L+1)) 80,130,80
   80      X(I)=(X(I)-S)/G(I,L+1)
C                                   ALG. STEP 32
      RETURN
C                                   ********************* MODE = 2
   90      DO 120 J=1,N
           S=ZERO
           IF (J.EQ.1) GO TO 110
           I1=MAX(1,J-NB+1)
           I2=J-1
                DO 100 I=I1,I2
                L=J-I+1+MAX(0,I-IP)
                S=S+X(I)*G(I,L)
  100 CONTINUE
  110      L=MAX(0,J-IP)
           IF (G(J,L+1)) 120,130,120
  120      X(J)=(X(J)-S)/G(J,L+1)
      RETURN
C
  130 CONTINUE
      NERR=1
      IOPT=2
      CALL XERMSG ('SLATEC', 'DBNDSL',
     +   'A ZERO DIAGONAL TERM IS IN THE N BY N UPPER TRIANGULAR ' //
     +   'MATRIX.', NERR, IOPT)
      RETURN
      END
