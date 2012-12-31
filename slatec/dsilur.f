*DECK DSILUR
      SUBROUTINE DSILUR (N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,
     +   ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)
C***BEGIN PROLOGUE  DSILUR
C***PURPOSE  Incomplete LU Iterative Refinement Sparse Ax = b Solver.
C            Routine to solve a general linear system  Ax = b  using
C            the incomplete LU decomposition with iterative refinement.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SSILUR-S, DSILUR-D)
C***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Greenbaum, Anne, (Courant Institute)
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
C     INTEGER ITER, IERR, IUNIT, LENW, IWORK(NL+NU+4*N+2), LENIW
C     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(NL+NU+4*N)
C
C     CALL DSILUR(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,
C    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right-hand side vector.
C X      :INOUT    Double Precision X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Double Precision A(NELT).
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See "Description",
C         below.  If the SLAP Triad format is chosen it is changed
C         internally to the SLAP Column format.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
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
C             COMMON /DSLBLK/ SOLN( )
C         If ITOL=11, iteration stops when the 2-norm of the difference
C         between the iterative approximation and the user-supplied
C         solution divided by the 2-norm of the user-supplied solution
C         is less than TOL.  Note that this requires the user to set up
C         the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.
C         The routine with this declaration should be loaded before the
C         stop test so that the correct length is used by the loader.
C         This procedure is not standard Fortran and may not work
C         correctly on your system (although it has worked on every
C         system the authors have tried).  If ITOL is not 11 then this
C         common block is indeed standard Fortran.
C TOL    :INOUT    Double Precision.
C         Convergence criterion, as described above.  (Reset if IERR=4.)
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or
C         ITMAX+1 if convergence criterion could not be achieved in
C         ITMAX iterations.
C ERR    :OUT      Double Precision.
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
C                       Reset to 500*D1MACH(3).  Iteration proceeded.
C           IERR = 5 => Preconditioning matrix, M, is not positive
C                       definite.  (r,z) < 0.
C           IERR = 6 => Matrix A is not positive definite.  (p,Ap) < 0.
C           IERR = 7 => Incomplete factorization broke down and was
C                       fudged.  Resulting preconditioning may be less
C                       than the best.
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration,
C         if this is desired for monitoring convergence.  If unit
C         number is 0, no writing will occur.
C RWORK  :WORK     Double Precision RWORK(LENW).
C         Double Precision array used for workspace.
C LENW   :IN       Integer.
C         Length of the double precision workspace, RWORK.
C         LENW >= NL+NU+4*N.
C         NL is the number of non-zeros in the lower triangle of the
C         matrix (including the diagonal).
C         NU is the number of non-zeros in the upper triangle of the
C         matrix (including the diagonal).
C IWORK  :WORK     Integer IWORK(LENIW).
C         Integer array used for workspace.
C         Upon return the following locations of IWORK hold information
C         which may be of use to the user:
C         IWORK(9)  Amount of Integer workspace actually used.
C         IWORK(10) Amount of Double Precision workspace actually used.
C LENIW  :IN       Integer.
C         Length of integer workspace, IWORK.  LENIW >= NL+NU+4*N+10.
C         NL is the number of non-zeros in the lower triangle of the
C         matrix (including the diagonal).
C         NU is the number of non-zeros in the upper triangle of the
C         matrix (including the diagonal).
C
C *Description
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
C       Column format.  The user can hand this routine either of the
C       of these data structures and SLAP  will figure out  which on
C       is being used and act accordingly.
C
C       =================== S L A P Triad format ===================
C
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero   matrix  element is  placed  in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
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
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except for  the diagonal entry, which
C       must appear first in each  "column")  and are stored  in the
C       double precision array A.   In other words,  for each column
C       in the matrix put the diagonal entry in  A.  Then put in the
C       other non-zero  elements going down  the column (except  the
C       diagonal) in order.   The  IA array holds the  row index for
C       each non-zero.  The JA array holds the offsets  into the IA,
C       A arrays  for  the  beginning  of each   column.   That  is,
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the
C       number of columns in  the matrix and NELT  is the number  of
C       non-zeros in the matrix.
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
C *Side Effects:
C       The SLAP Triad format (IA, JA, A) is modified internally to be
C       the SLAP Column format.  See above.
C
C *Cautions:
C     This routine will attempt to write to the Fortran logical output
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
C     this logical unit is attached to a file or terminal before calling
C     this routine with a non-zero value for IUNIT.  This routine does
C     not check for the validity of a non-zero IUNIT unit number.
C
C***SEE ALSO  DSJAC, DSGS, DIR
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DCHKW, DIR, DS2Y, DSILUS, DSLUI, DSMV
C***REVISION HISTORY  (YYMMDD)
C   871119  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890921  Removed TeX from comments.  (FNF)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)
C   920511  Added complete declaration section.  (WRB)
C   921019  Corrected NEL to NL.  (FNF)
C***END PROLOGUE  DSILUR
C     .. Parameters ..
      INTEGER LOCRB, LOCIB
      PARAMETER (LOCRB=1, LOCIB=11)
C     .. Scalar Arguments ..
      DOUBLE PRECISION ERR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LENIW, LENW, N, NELT
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), RWORK(LENW), X(N)
      INTEGER IA(NELT), IWORK(LENIW), JA(NELT)
C     .. Local Scalars ..
      INTEGER ICOL, J, JBGN, JEND, LOCDIN, LOCDZ, LOCIL, LOCIU, LOCIW,
     +        LOCJL, LOCJU, LOCL, LOCNC, LOCNR, LOCR, LOCU, LOCW, LOCZ,
     +        NL, NU
C     .. External Subroutines ..
      EXTERNAL DCHKW, DIR, DS2Y, DSILUS, DSLUI, DSMV
C***FIRST EXECUTABLE STATEMENT  DSILUR
C
      IERR = 0
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
C
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format.
      CALL DS2Y( N, NELT, IA, JA, A, ISYM )
C
C         Count number of Non-Zero elements in preconditioner ILU
C         matrix.  Then set up the work arrays.
      NL = 0
      NU = 0
      DO 20 ICOL = 1, N
C         Don't count diagonal.
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
CVD$ NOVECTOR
            DO 10 J = JBGN, JEND
               IF( IA(J).GT.ICOL ) THEN
                  NL = NL + 1
                  IF( ISYM.NE.0 ) NU = NU + 1
               ELSE
                  NU = NU + 1
               ENDIF
 10         CONTINUE
         ENDIF
 20   CONTINUE
C
      LOCIL = LOCIB
      LOCJL = LOCIL + N+1
      LOCIU = LOCJL + NL
      LOCJU = LOCIU + NU
      LOCNR = LOCJU + N+1
      LOCNC = LOCNR + N
      LOCIW = LOCNC + N
C
      LOCL = LOCRB
      LOCDIN = LOCL + NL
      LOCU = LOCDIN + N
      LOCR = LOCU + NU
      LOCZ = LOCR + N
      LOCDZ = LOCZ + N
      LOCW = LOCDZ + N
C
C         Check the workspace allocations.
      CALL DCHKW( 'DSILUR', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
C
      IWORK(1) = LOCIL
      IWORK(2) = LOCJL
      IWORK(3) = LOCIU
      IWORK(4) = LOCJU
      IWORK(5) = LOCL
      IWORK(6) = LOCDIN
      IWORK(7) = LOCU
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
C         Compute the Incomplete LU decomposition.
      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),
     $     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),
     $     IWORK(LOCJU), RWORK(LOCU), IWORK(LOCNR), IWORK(LOCNC) )
C
C         Do the Preconditioned Iterative Refinement iteration.
      CALL DIR(N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSLUI,
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK(LOCR),
     $     RWORK(LOCZ), RWORK(LOCDZ), RWORK, IWORK)
      RETURN
C------------- LAST LINE OF DSILUR FOLLOWS ----------------------------
      END
