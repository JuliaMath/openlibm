*DECK SCGN
      SUBROUTINE SCGN (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC,
     +   MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, ATP,
     +   ATZ, DZ, ATDZ, RWORK, IWORK)
C***BEGIN PROLOGUE  SCGN
C***PURPOSE  Preconditioned CG Sparse Ax=b Solver for Normal Equations.
C            Routine to solve a general linear system  Ax = b  using the
C            Preconditioned Conjugate Gradient method applied to the
C            normal equations  AA'y = b, x=A'y.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      SINGLE PRECISION (SCGN-S, DCGN-D)
C***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM SOLVE,
C             NORMAL EQUATIONS., SLAP, SPARSE
C***AUTHOR  Greenbaum, Anne, (Courant Institute)
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C *Usage:
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
C     INTEGER  ITER, IERR, IUNIT, IWORK(USER DEFINED)
C     REAL     B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N)
C     REAL     P(N), ATP(N), ATZ(N), DZ(N), ATDZ(N)
C     REAL     RWORK(USER DEFINED)
C     EXTERNAL MATVEC, MTTVEC, MSOLVE
C
C     CALL SCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC,
C    $     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R,
C    $     Z, P, ATP, ATZ, DZ, ATDZ, RWORK, IWORK)
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right-hand side vector.
C X      :INOUT    Real X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Real A(NELT).
C         These arrays contain the matrix data structure for A.
C         It could take any form.  See "Description", below,
C         for more details.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C MATVEC :EXT      External.
C         Name of a routine which performs the matrix vector multiply
C         y = A*X given A and X.  The name of the MATVEC routine must
C         be declared external in the calling program.  The calling
C         sequence to MATVEC is:
C             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM )
C         Where N is the number of unknowns, Y is the product A*X
C         upon return X is an input vector, NELT is the number of
C         non-zeros in the SLAP-Column IA, JA, A storage for the matrix
C         A.  ISYM is a flag which, if non-zero, denotes that A is
C         symmetric and only the lower or upper triangle is stored.
C MTTVEC :EXT      External.
C         Name of a routine which performs the matrix transpose vector
C         multiply y = A'*X given A and X (where ' denotes transpose).
C         The name of the MTTVEC routine must be declared external in
C         the calling program.  The calling sequence to MTTVEC is the
C         same as that for MATVEC, viz.:
C             CALL MTTVEC( N, X, Y, NELT, IA, JA, A, ISYM )
C         Where N is the number of unknowns, Y is the product A'*X
C         upon return X is an input vector, NELT is the number of
C         non-zeros in the SLAP-Column IA, JA, A storage for the matrix
C         A.  ISYM is a flag which, if non-zero, denotes that A is
C         symmetric and only the lower or upper triangle is stored.
C MSOLVE :EXT      External.
C         Name of a routine which solves a linear system MZ = R for
C         Z given R with the preconditioning matrix M (M is supplied via
C         RWORK and IWORK arrays).  The name of the MSOLVE routine must
C         be declared external in the calling program.  The calling
C         sequence to MSOLVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         Where N is the number of unknowns, R is the right-hand side
C         vector and Z is the solution upon return.  NELT, IA, JA, A and
C         ISYM are defined as above.  RWORK is a real array that can
C         be used to pass necessary preconditioning information and/or
C         workspace to MSOLVE.  IWORK is an integer work array for
C         the same purpose as RWORK.
C ITOL   :IN       Integer.
C         Flag to indicate type of convergence criterion.
C         If ITOL=1, iteration stops when the 2-norm of the residual
C         divided by the 2-norm of the right-hand side is less than TOL.
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the
C         residual divided by the 2-norm of M-inv times the right hand
C         side is less than TOL, where M-inv is the inverse of the
C         diagonal of A.
C         ITOL=11 is often useful for checking and comparing different
C         routines.  For this case, the user must supply the "exact"
C         solution or a very accurate approximation (one with an error
C         much less than TOL) through a common block,
C             COMMON /SSLBLK/ SOLN( )
C         If ITOL=11, iteration stops when the 2-norm of the difference
C         between the iterative approximation and the user-supplied
C         solution divided by the 2-norm of the user-supplied solution
C         is less than TOL.  Note that this requires the user to set up
C         the "COMMON /SSLBLK/ SOLN(LENGTH)" in the calling routine.
C         The routine with this declaration should be loaded before the
C         stop test so that the correct length is used by the loader.
C         This procedure is not standard Fortran and may not work
C         correctly on your system (although it has worked on every
C         system the authors have tried).  If ITOL is not 11 then this
C         common block is indeed standard Fortran.
C TOL    :INOUT    Real.
C         Convergence criterion, as described above.  (Reset if IERR=4.)
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or
C         ITMAX+1 if convergence criterion could not be achieved in
C         ITMAX iterations.
C ERR    :OUT      Real.
C         Error estimate of error in final approximate solution, as
C         defined by ITOL.
C IERR   :OUT      Integer.
C         Return error flag.
C           IERR = 0 => All went well.
C           IERR = 1 => Insufficient space allocated for WORK or IWORK.
C           IERR = 2 => Method failed to converge in ITMAX steps.
C           IERR = 3 => Error in user input.
C                       Check input values of N, ITOL.
C           IERR = 4 => User error tolerance set too tight.
C                       Reset to 500*R1MACH(3).  Iteration proceeded.
C           IERR = 5 => Preconditioning matrix, M, is not positive
C                       definite.  (r,z) < 0.
C           IERR = 6 => Matrix A is not positive definite.  (p,Ap) < 0.
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration,
C         if this is desired for monitoring convergence.  If unit
C         number is 0, no writing will occur.
C R      :WORK     Real R(N).
C Z      :WORK     Real Z(N).
C P      :WORK     Real P(N).
C ATP    :WORK     Real ATP(N).
C ATZ    :WORK     Real ATZ(N).
C DZ     :WORK     Real DZ(N).
C ATDZ   :WORK     Real ATDZ(N).
C         Real arrays used for workspace.
C RWORK  :WORK     Real RWORK(USER DEFINED).
C         Real array that can be used by  MSOLVE.
C IWORK  :WORK     Integer IWORK(USER DEFINED).
C         Integer array that can be used by  MSOLVE.
C
C *Description:
C       This  routine applies the  preconditioned conjugate gradient
C       (PCG) method to a non-symmetric system of equations Ax=b. To
C       do this the normal equations are solved:
C               AA' y  = b, where  x  = A'y.
C       In PCG method the iteration count is determined by condition
C                               -1
C       number of the  matrix (M  A).   In the  situation where  the
C       normal equations are  used  to solve a  non-symmetric system
C       the condition number depends on  AA' and should therefore be
C       much worse than that of A.  This is the conventional wisdom.
C       When one has a good preconditioner for AA' this may not hold.
C       The latter is the situation when SCGN should be tried.
C
C       If one is trying to solve  a symmetric system, SCG should be
C       used instead.
C
C       This routine does  not care  what matrix data   structure is
C       used for A and M.  It simply calls MATVEC, MTTVEC and MSOLVE
C       routines, with arguments as described above.  The user could
C       write any type of structure, and  appropriate MATVEC, MTTVEC
C       and MSOLVE routines.  It is assumed  that A is stored in the
C       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is
C       stored  in  IWORK  and  RWORK)  in  some fashion.   The SLAP
C       routines SSDCGN and SSLUCN are examples of this procedure.
C
C       Two  examples  of  matrix  data structures  are the: 1) SLAP
C       Triad  format and 2) SLAP Column format.
C
C       =================== S L A P Triad format ===================
C
C       In  this   format only the  non-zeros are  stored.  They may
C       appear  in *ANY* order.   The user  supplies three arrays of
C       length NELT, where  NELT  is the number  of non-zeros in the
C       matrix:  (IA(NELT), JA(NELT),  A(NELT)).  For each  non-zero
C       the  user puts   the row  and  column index   of that matrix
C       element in the IA and JA arrays.  The  value of the non-zero
C       matrix  element is  placed in  the corresponding location of
C       the A  array.  This is  an extremely easy data  structure to
C       generate.  On  the other hand it  is  not too  efficient  on
C       vector  computers   for the  iterative  solution  of  linear
C       systems.  Hence, SLAP  changes this input  data structure to
C       the SLAP   Column  format for the  iteration (but   does not
C       change it back).
C
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       =================== S L A P Column format ==================
C
C       In  this format   the non-zeros are    stored counting  down
C       columns (except  for the diagonal  entry, which must  appear
C       first in each "column") and are  stored in the real array A.
C       In other words,  for  each column    in the matrix   put the
C       diagonal  entry  in A.   Then   put  in the  other  non-zero
C       elements going   down the  column (except  the  diagonal) in
C       order.  The IA array holds the row index  for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning   of   each  column.      That is,   IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1), A(JA(ICOL+1)-1)  points to the
C       end of the ICOL-th column.  Note that we always have JA(N+1)
C       = NELT+1, where N is the number of columns in the matrix and
C       NELT is the number of non-zeros in the matrix.
C
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C *Cautions:
C     This routine will attempt to write to the Fortran logical output
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
C     this logical unit is attached to a file or terminal before calling
C     this routine with a non-zero value for IUNIT.  This routine does
C     not check for the validity of a non-zero IUNIT unit number.
C
C***SEE ALSO  SSDCGN, SSLUCN, ISSCGN
C***REFERENCES  1. Mark K. Seager, A SLAP for the Masses, in
C                  G. F. Carey, Ed., Parallel Supercomputing: Methods,
C                  Algorithms and Applications, Wiley, 1989, pp.135-155.
C***ROUTINES CALLED  ISSCGN, R1MACH, SAXPY, SCOPY, SDOT
C***REVISION HISTORY  (YYMMDD)
C   871119  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890921  Removed TeX from comments.  (FNF)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   891004  Added new reference.
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910502  Removed MATVEC, MTTVEC and MSOLVE from ROUTINES CALLED
C           list.  (FNF)
C   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
C   920511  Added complete declaration section.  (WRB)
C   920929  Corrected format of reference.  (FNF)
C   921019  Changed 500.0 to 500 to reduce SP/DP differences.  (FNF)
C   921113  Corrected C***CATEGORY line.  (FNF)
C***END PROLOGUE  SCGN
C     .. Scalar Arguments ..
      REAL ERR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, N, NELT
C     .. Array Arguments ..
      REAL A(N), ATDZ(N), ATP(N), ATZ(N), B(N), DZ(N), P(N), R(N),
     +     RWORK(*), X(N), Z(N)
      INTEGER IA(NELT), IWORK(*), JA(NELT)
C     .. Subroutine Arguments ..
      EXTERNAL MATVEC, MSOLVE, MTTVEC
C     .. Local Scalars ..
      REAL AK, AKDEN, BK, BKDEN, BKNUM, BNRM, SOLNRM, TOLMIN
      INTEGER I, K
C     .. External Functions ..
      REAL R1MACH, SDOT
      INTEGER ISSCGN
      EXTERNAL R1MACH, SDOT, ISSCGN
C     .. External Subroutines ..
      EXTERNAL SAXPY, SCOPY
C***FIRST EXECUTABLE STATEMENT  SCGN
C
C         Check user input.
C
      ITER = 0
      IERR = 0
      IF( N.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      TOLMIN = 500*R1MACH(3)
      IF( TOL.LT.TOLMIN ) THEN
         TOL = TOLMIN
         IERR = 4
      ENDIF
C         Calculate initial residual and pseudo-residual, and check
C         stopping criterion.
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
      DO 10 I = 1, N
         R(I) = B(I) - R(I)
 10   CONTINUE
      CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
      CALL MTTVEC(N, Z, ATZ, NELT, IA, JA, A, ISYM)
C
      IF( ISSCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC, MSOLVE,
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, ATP, ATZ,
     $     DZ, ATDZ, RWORK, IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 )
     $     GO TO 200
      IF( IERR.NE.0 ) RETURN
C
C         ***** iteration loop *****
C
      DO 100 K=1,ITMAX
         ITER = K
C
C         Calculate coefficient BK and direction vector P.
         BKNUM = SDOT(N, Z, 1, R, 1)
         IF( BKNUM.LE.0.0E0 ) THEN
            IERR = 6
            RETURN
         ENDIF
         IF(ITER .EQ. 1) THEN
            CALL SCOPY(N, Z, 1, P, 1)
         ELSE
            BK = BKNUM/BKDEN
            DO 20 I = 1, N
               P(I) = Z(I) + BK*P(I)
 20         CONTINUE
         ENDIF
         BKDEN = BKNUM
C
C         Calculate coefficient AK, new iterate X, new residual R,
C         and new pseudo-residual ATZ.
         IF(ITER .NE. 1) CALL SAXPY(N, BK, ATP, 1, ATZ, 1)
         CALL SCOPY(N, ATZ, 1, ATP, 1)
         AKDEN = SDOT(N, ATP, 1, ATP, 1)
         IF( AKDEN.LE.0.0E0 ) THEN
            IERR = 6
            RETURN
         ENDIF
         AK = BKNUM/AKDEN
         CALL SAXPY(N, AK, ATP, 1, X, 1)
         CALL MATVEC(N, ATP, Z, NELT, IA, JA, A, ISYM)
         CALL SAXPY(N, -AK, Z, 1, R, 1)
         CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
         CALL MTTVEC(N, Z, ATZ, NELT, IA, JA, A, ISYM)
C
C         check stopping criterion.
         IF( ISSCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC,
     $        MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R,
     $        Z, P, ATP, ATZ, DZ, ATDZ, RWORK, IWORK, AK, BK, BNRM,
     $        SOLNRM) .NE. 0) GOTO 200
C
 100  CONTINUE
C
C         *****   end of loop  *****
C
C         stopping criterion not satisfied.
      ITER = ITMAX + 1
C
 200  RETURN
C------------- LAST LINE OF SCGN FOLLOWS ----------------------------
      END
