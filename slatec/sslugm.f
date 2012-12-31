*DECK SSLUGM
      SUBROUTINE SSLUGM (N, B, X, NELT, IA, JA, A, ISYM, NSAVE, ITOL,
     +   TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)
C***BEGIN PROLOGUE  SSLUGM
C***PURPOSE  Incomplete LU GMRES Iterative Sparse Ax=b Solver.
C            This routine uses the generalized minimum residual
C            (GMRES) method with incomplete LU factorization for
C            preconditioning to solve possibly non-symmetric linear
C            systems of the form: Ax = b.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      SINGLE PRECISION (SSLUGM-S, DSLUGM-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C
C *Usage:
C      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL
C      INTEGER   ITMAX, ITER, IERR, IUNIT, LENW, IWORK(LENIW), LENIW
C      REAL      B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
C
C      CALL SSLUGM(N, B, X, NELT, IA, JA, A, ISYM, NSAVE,
C     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
C     $     RWORK, LENW, IWORK, LENIW)
C
C *Arguments:
C N      :IN       Integer.
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
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See "Description",
C         below.  If the SLAP Triad format is chosen it is changed
C         internally to the SLAP Column format.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C NSAVE  :IN       Integer.
C         Number of direction vectors to save and orthogonalize against.
C         Must be greater than 1.
C ITOL   :IN       Integer.
C         Flag to indicate the type of convergence criterion used.
C         ITOL=0  Means the  iteration stops when the test described
C                 below on  the  residual RL  is satisfied.  This is
C                 the  "Natural Stopping Criteria" for this routine.
C                 Other values  of   ITOL  cause  extra,   otherwise
C                 unnecessary, computation per iteration and     are
C                 therefore  much less  efficient.  See  ISSGMR (the
C                 stop test routine) for more information.
C         ITOL=1  Means   the  iteration stops   when the first test
C                 described below on  the residual RL  is satisfied,
C                 and there  is either right  or  no preconditioning
C                 being used.
C         ITOL=2  Implies     that   the  user    is   using    left
C                 preconditioning, and the second stopping criterion
C                 below is used.
C         ITOL=3  Means the  iteration stops   when  the  third test
C                 described below on Minv*Residual is satisfied, and
C                 there is either left  or no  preconditioning begin
C                 used.
C         ITOL=11 is    often  useful  for   checking  and comparing
C                 different routines.  For this case, the  user must
C                 supply  the  "exact" solution or  a  very accurate
C                 approximation (one with  an  error much less  than
C                 TOL) through a common block,
C                     COMMON /SSLBLK/ SOLN( )
C                 If ITOL=11, iteration stops when the 2-norm of the
C                 difference between the iterative approximation and
C                 the user-supplied solution  divided by the  2-norm
C                 of the  user-supplied solution  is  less than TOL.
C                 Note that this requires  the  user to  set up  the
C                 "COMMON     /SSLBLK/ SOLN(LENGTH)"  in the calling
C                 routine.  The routine with this declaration should
C                 be loaded before the stop test so that the correct
C                 length is used by  the loader.  This procedure  is
C                 not standard Fortran and may not work correctly on
C                 your   system (although  it  has  worked  on every
C                 system the authors have tried).  If ITOL is not 11
C                 then this common block is indeed standard Fortran.
C TOL    :INOUT    Real.
C         Convergence criterion, as described below.  If TOL is set
C         to zero on input, then a default value of 500*(the smallest
C         positive magnitude, machine epsilon) is used.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.  This routine uses the default
C         of NRMAX = ITMAX/NSAVE to determine the when each restart
C         should occur.  See the description of NRMAX and MAXL in
C         SGMRES for a full and frightfully interesting discussion of
C         this topic.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or
C         ITMAX+1 if convergence criterion could not be achieved in
C         ITMAX iterations.
C ERR    :OUT      Real.
C         Error estimate of error in final approximate solution, as
C         defined by ITOL.  Letting norm() denote the Euclidean
C         norm, ERR is defined as follows...
C         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               for right or no preconditioning, and
C                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               for left preconditioning.
C         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
C                               since right or no preconditioning
C                               being used.
C         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
C                                norm(SB*(M-inverse)*B),
C                               since left preconditioning is being
C                               used.
C         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
C                               i=1,n
C         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
C IERR   :OUT      Integer.
C         Return error flag.
C               IERR = 0 => All went well.
C               IERR = 1 => Insufficient storage allocated for
C                           RGWK or IGWK.
C               IERR = 2 => Routine SPIGMR failed to reduce the norm
C                           of the current residual on its last call,
C                           and so the iteration has stalled.  In
C                           this case, X equals the last computed
C                           approximation.  The user must either
C                           increase MAXL, or choose a different
C                           initial guess.
C               IERR =-1 => Insufficient length for RGWK array.
C                           IGWK(6) contains the required minimum
C                           length of the RGWK array.
C               IERR =-2 => Inconsistent ITOL and JPRE values.
C         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the
C         left-hand-side of the relevant stopping test defined
C         below associated with the residual for the current
C         approximation X(L).
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration,
C         if this is desired for monitoring convergence.  If unit
C         number is 0, no writing will occur.
C RWORK  :WORK    Real RWORK(LENW).
C         Real array of size LENW.
C LENW   :IN       Integer.
C         Length of the real workspace, RWORK.
C         LENW >= 1 + N*(NSAVE+7) +  NSAVE*(NSAVE+3)+NL+NU.
C         Here NL is the number of non-zeros in the lower triangle of
C         the matrix (including the diagonal) and NU is the number of
C         non-zeros in the upper triangle of the matrix (including the
C         diagonal).
C         For the recommended values,  RWORK  has size at least
C         131 + 17*N + NL + NU.
C IWORK  :INOUT    Integer IWORK(LENIW).
C         Used to hold pointers into the RWORK array.
C         Upon return the following locations of IWORK hold information
C         which may be of use to the user:
C         IWORK(9)  Amount of Integer workspace actually used.
C         IWORK(10) Amount of Real workspace actually used.
C LENIW  :IN       Integer.
C         Length of the integer workspace, IWORK.
C         LENIW >= NL+NU+4*N+32.
C
C *Description:
C       SSLUGM solves a linear system A*X = B rewritten in the form:
C
C        (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B,
C
C       with right preconditioning, or
C
C        (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B,
C
C       with left preconditioning, where A is an n-by-n real matrix,
C       X and  B are  N-vectors,  SB and  SX  are   diagonal scaling
C       matrices, and M is the Incomplete LU factorization of A.  It
C       uses preconditioned  Krylov subpace   methods  based on  the
C       generalized minimum residual  method (GMRES).   This routine
C       is a  driver  routine  which  assumes a SLAP   matrix   data
C       structure   and  sets  up  the  necessary  information to do
C       diagonal  preconditioning  and calls the main GMRES  routine
C       SGMRES for the   solution   of the linear   system.   SGMRES
C       optionally   performs  either  the full    orthogonalization
C       version of the  GMRES algorithm or  an incomplete variant of
C       it.  Both versions use restarting of the linear iteration by
C       default, although the user can disable this feature.
C
C       The GMRES  algorithm generates a sequence  of approximations
C       X(L) to the  true solution of the above  linear system.  The
C       convergence criteria for stopping the  iteration is based on
C       the size  of the  scaled norm of  the residual  R(L)  =  B -
C       A*X(L).  The actual stopping test is either:
C
C               norm(SB*(B-A*X(L))) .le. TOL*norm(SB*B),
C
C       for right preconditioning, or
C
C               norm(SB*(M-inverse)*(B-A*X(L))) .le.
C                       TOL*norm(SB*(M-inverse)*B),
C
C       for left preconditioning, where norm() denotes the Euclidean
C       norm, and TOL is  a positive scalar less  than one  input by
C       the user.  If TOL equals zero  when SSLUGM is called, then a
C       default  value  of 500*(the   smallest  positive  magnitude,
C       machine epsilon) is used.  If the  scaling arrays SB  and SX
C       are used, then  ideally they  should be chosen  so  that the
C       vectors SX*X(or SX*M*X) and  SB*B have all their  components
C       approximately equal  to  one in  magnitude.  If one wants to
C       use the same scaling in X  and B, then  SB and SX can be the
C       same array in the calling program.
C
C       The following is a list of the other routines and their
C       functions used by GMRES:
C       SGMRES  Contains the matrix structure independent driver
C               routine for GMRES.
C       SPIGMR  Contains the main iteration loop for GMRES.
C       SORTH   Orthogonalizes a new vector against older basis vectors.
C       SHEQR   Computes a QR decomposition of a Hessenberg matrix.
C       SHELS   Solves a Hessenberg least-squares system, using QR
C               factors.
C       RLCALC  Computes the scaled residual RL.
C       XLCALC  Computes the solution XL.
C       ISSGMR  User-replaceable stopping routine.
C
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
C       Column format.  The user can hand this routine either of the
C       of these data structures and SLAP  will figure out  which on
C       is being used and act accordingly.
C
C       =================== S L A P Triad format ===================
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
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
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
C***REFERENCES  1. Peter N. Brown and A. C. Hindmarsh, Reduced Storage
C                  Matrix Methods in Stiff ODE Systems, Lawrence Liver-
C                  more National Laboratory Report UCRL-95088, Rev. 1,
C                  Livermore, California, June 1987.
C***ROUTINES CALLED  SCHKW, SGMRES, SS2Y, SSILUS, SSLUI, SSMV
C***REVISION HISTORY  (YYMMDD)
C   880615  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
C   920511  Added complete declaration section.  (WRB)
C   920929  Corrected format of references.  (FNF)
C   921019  Corrected NEL to NL.  (FNF)
C***END PROLOGUE  SSLUGM
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Parameters ..
      INTEGER LOCRB, LOCIB
      PARAMETER (LOCRB=1, LOCIB=11)
C     .. Scalar Arguments ..
      REAL ERR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LENIW, LENW, N,
     +        NELT, NSAVE
C     .. Array Arguments ..
      REAL A(NELT), B(N), RWORK(LENW), X(N)
      INTEGER IA(NELT), IWORK(LENIW), JA(NELT)
C     .. Local Scalars ..
      INTEGER ICOL, J, JBGN, JEND, LOCDIN, LOCIGW, LOCIL, LOCIU, LOCIW,
     +        LOCJL, LOCJU, LOCL, LOCNC, LOCNR, LOCRGW, LOCU, LOCW,
     +        MYITOL, NL, NU
C     .. External Subroutines ..
      EXTERNAL SCHKW, SGMRES, SS2Y, SSILUS, SSLUI, SSMV
C***FIRST EXECUTABLE STATEMENT  SSLUGM
C
      IERR = 0
      ERR  = 0
      IF( NSAVE.LE.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
C
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format.
      CALL SS2Y( N, NELT, IA, JA, A, ISYM )
C
C         Count number of Non-Zero elements preconditioner ILU matrix.
C         Then set up the work arrays.  We assume MAXL=KMP=NSAVE.
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
      LOCIGW = LOCIB
      LOCIL = LOCIGW + 20
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
      LOCRGW = LOCU + NU
      LOCW = LOCRGW + 1+N*(NSAVE+6)+NSAVE*(NSAVE+3)
C
C         Check the workspace allocations.
      CALL SCHKW( 'SSLUGM', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
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
      CALL SSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),
     $     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),
     $     IWORK(LOCJU), RWORK(LOCU), IWORK(LOCNR), IWORK(LOCNC) )
C
C         Perform the Incomplete LU Preconditioned Generalized Minimum
C         Residual iteration algorithm.  The following SGMRES
C         defaults are used MAXL = KMP = NSAVE, JSCAL = 0,
C         JPRE = -1, NRMAX = ITMAX/NSAVE
      IWORK(LOCIGW  ) = NSAVE
      IWORK(LOCIGW+1) = NSAVE
      IWORK(LOCIGW+2) = 0
      IWORK(LOCIGW+3) = -1
      IWORK(LOCIGW+4) = ITMAX/NSAVE
      MYITOL = 0
C
      CALL SGMRES( N, B, X, NELT, IA, JA, A, ISYM, SSMV, SSLUI,
     $     MYITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, RWORK,
     $     RWORK(LOCRGW), LENW-LOCRGW, IWORK(LOCIGW), 20,
     $     RWORK, IWORK )
C
      IF( ITER.GT.ITMAX ) IERR = 2
      RETURN
C------------- LAST LINE OF SSLUGM FOLLOWS ----------------------------
      END
