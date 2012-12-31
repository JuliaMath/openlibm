*DECK SGMRES
      SUBROUTINE SGMRES (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
     +   ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, RGWK, LRGW,
     +   IGWK, LIGW, RWORK, IWORK)
C***BEGIN PROLOGUE  SGMRES
C***PURPOSE  Preconditioned GMRES Iterative Sparse Ax=b Solver.
C            This routine uses the generalized minimum residual
C            (GMRES) method with preconditioning to solve
C            non-symmetric linear systems of the form: Ax = b.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      SINGLE PRECISION (SGMRES-S, DGMRES-D)
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
C      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
C      INTEGER   ITER, IERR, IUNIT, LRGW, IGWK(LIGW), LIGW
C      INTEGER   IWORK(USER DEFINED)
C      REAL      B(N), X(N), A(NELT), TOL, ERR, SB(N), SX(N)
C      REAL      RGWK(LRGW), RWORK(USER DEFINED)
C      EXTERNAL  MATVEC, MSOLVE
C
C      CALL SGMRES(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
C     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX,
C     $     RGWK, LRGW, IGWK, LIGW, RWORK, IWORK)
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right-hand side vector.
C X      :INOUT    Real X(N).
C         On input X is your initial guess for the solution vector.
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
C         Y = A*X given A and X.  The name of the MATVEC routine must
C         be declared external in the calling program.  The calling
C         sequence to MATVEC is:
C             CALL MATVEC(N, X, Y, NELT, IA, JA, A, ISYM)
C         where N is the number of unknowns, Y is the product A*X
C         upon return, X is an input vector, and NELT is the number of
C         non-zeros in the SLAP IA, JA, A storage for the matrix A.
C         ISYM is a flag which, if non-zero, denotes that A is
C         symmetric and only the lower or upper triangle is stored.
C MSOLVE :EXT      External.
C         Name of the routine which solves a linear system Mz = r for
C         z given r with the preconditioning matrix M (M is supplied via
C         RWORK and IWORK arrays.  The name of the MSOLVE routine must
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
C                 there is either left  or no  preconditioning being
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
C ITMAX  :DUMMY    Integer.
C         Maximum number of iterations in most SLAP routines.  In
C         this routine this does not make sense.  The maximum number
C         of iterations here is given by ITMAX = MAXL*(NRMAX+1).
C         See IGWK for definitions of MAXL and NRMAX.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or
C         ITMAX if convergence criterion could not be achieved in
C         ITMAX iterations.
C ERR    :OUT      Real.
C         Error estimate of error in final approximate solution, as
C         defined by ITOL.  Letting norm() denote the Euclidean
C         norm, ERR is defined as follows..
C
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
C               IERR = 2 => Routine SGMRES failed to reduce the norm
C                           of the current residual on its last call,
C                           and so the iteration has stalled.  In
C                           this case, X equals the last computed
C                           approximation.  The user must either
C                           increase MAXL, or choose a different
C                           initial guess.
C               IERR =-1 => Insufficient length for RGWK array.
C                           IGWK(6) contains the required minimum
C                           length of the RGWK array.
C               IERR =-2 => Illegal value of ITOL, or ITOL and JPRE
C                           values are inconsistent.
C         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the
C         left-hand-side of the relevant stopping test defined
C         below associated with the residual for the current
C         approximation X(L).
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration,
C         if this is desired for monitoring convergence.  If unit
C         number is 0, no writing will occur.
C SB     :IN       Real SB(N).
C         Array of length N containing scale factors for the right
C         hand side vector B.  If JSCAL.eq.0 (see below), SB need
C         not be supplied.
C SX     :IN       Real SX(N).
C         Array of length N containing scale factors for the solution
C         vector X.  If JSCAL.eq.0 (see below), SX need not be
C         supplied.  SB and SX can be the same array in the calling
C         program if desired.
C RGWK   :INOUT    Real RGWK(LRGW).
C         Real array used for workspace by SGMRES.
C         On return, RGWK(1) = RHOL.  See IERR for definition of RHOL.
C LRGW   :IN       Integer.
C         Length of the real workspace, RGWK.
C         LRGW >= 1 + N*(MAXL+6) + MAXL*(MAXL+3).
C         See below for definition of MAXL.
C         For the default values, RGWK has size at least 131 + 16*N.
C IGWK   :INOUT    Integer IGWK(LIGW).
C         The following IGWK parameters should be set by the user
C         before calling this routine.
C         IGWK(1) = MAXL.  Maximum dimension of Krylov subspace in
C            which X - X0 is to be found (where, X0 is the initial
C            guess).  The default value of MAXL is 10.
C         IGWK(2) = KMP.  Maximum number of previous Krylov basis
C            vectors to which each new basis vector is made orthogonal.
C            The default value of KMP is MAXL.
C         IGWK(3) = JSCAL.  Flag indicating whether the scaling
C            arrays SB and SX are to be used.
C            JSCAL = 0 => SB and SX are not used and the algorithm
C               will perform as if all SB(I) = 1 and SX(I) = 1.
C            JSCAL = 1 =>  Only SX is used, and the algorithm
C               performs as if all SB(I) = 1.
C            JSCAL = 2 =>  Only SB is used, and the algorithm
C               performs as if all SX(I) = 1.
C            JSCAL = 3 =>  Both SB and SX are used.
C         IGWK(4) = JPRE.  Flag indicating whether preconditioning
C            is being used.
C            JPRE = 0  =>  There is no preconditioning.
C            JPRE > 0  =>  There is preconditioning on the right
C               only, and the solver will call routine MSOLVE.
C            JPRE < 0  =>  There is preconditioning on the left
C               only, and the solver will call routine MSOLVE.
C         IGWK(5) = NRMAX.  Maximum number of restarts of the
C            Krylov iteration.  The default value of NRMAX = 10.
C            if IWORK(5) = -1,  then no restarts are performed (in
C            this case, NRMAX is set to zero internally).
C         The following IWORK parameters are diagnostic information
C         made available to the user after this routine completes.
C         IGWK(6) = MLWK.  Required minimum length of RGWK array.
C         IGWK(7) = NMS.  The total number of calls to MSOLVE.
C LIGW   :IN       Integer.
C         Length of the integer workspace, IGWK.  LIGW >= 20.
C RWORK  :WORK     Real RWORK(USER DEFINED).
C         Real array that can be used for workspace in MSOLVE.
C IWORK  :WORK     Integer IWORK(USER DEFINED).
C         Integer array that can be used for workspace in MSOLVE.
C
C *Description:
C       SGMRES solves a linear system A*X = B rewritten in the form:
C
C        (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B,
C
C       with right preconditioning, or
C
C        (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B,
C
C       with left preconditioning, where A is an N-by-N real matrix,
C       X  and  B are N-vectors,   SB and SX   are  diagonal scaling
C       matrices,   and M is  a preconditioning    matrix.   It uses
C       preconditioned  Krylov   subpace  methods  based     on  the
C       generalized minimum residual  method (GMRES).   This routine
C       optionally performs  either  the  full     orthogonalization
C       version of the  GMRES  algorithm or an incomplete variant of
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
C       the user.  If TOL equals zero  when SGMRES is called, then a
C       default  value  of 500*(the   smallest  positive  magnitude,
C       machine epsilon) is used.  If the  scaling arrays SB  and SX
C       are used, then  ideally they  should be chosen  so  that the
C       vectors SX*X(or SX*M*X) and  SB*B have all their  components
C       approximately equal  to  one in  magnitude.  If one wants to
C       use the same scaling in X  and B, then  SB and SX can be the
C       same array in the calling program.
C
C       The following is a list of the other routines and their
C       functions used by SGMRES:
C       SPIGMR  Contains the main iteration loop for GMRES.
C       SORTH   Orthogonalizes a new vector against older basis vectors.
C       SHEQR   Computes a QR decomposition of a Hessenberg matrix.
C       SHELS   Solves a Hessenberg least-squares system, using QR
C               factors.
C       SRLCAL  Computes the scaled residual RL.
C       SXLCAL  Computes the solution XL.
C       ISSGMR  User-replaceable stopping routine.
C
C       This routine does  not care  what matrix data   structure is
C       used for  A and M.  It simply   calls  the MATVEC and MSOLVE
C       routines, with  the arguments as  described above.  The user
C       could write any type of structure and the appropriate MATVEC
C       and MSOLVE routines.  It is assumed  that A is stored in the
C       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is
C       stored  in  IWORK  and  RWORK   in  some fashion.   The SLAP
C       routines SSDCG and SSICCG are examples of this procedure.
C
C       Two  examples  of  matrix  data structures  are the: 1) SLAP
C       Triad  format and 2) SLAP Column format.
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
C               2. Mark K. Seager, A SLAP for the Masses, in
C                  G. F. Carey, Ed., Parallel Supercomputing: Methods,
C                  Algorithms and Applications, Wiley, 1989, pp.135-155.
C***ROUTINES CALLED  R1MACH, SCOPY, SNRM2, SPIGMR
C***REVISION HISTORY  (YYMMDD)
C   871001  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   891004  Added new reference.
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910506  Corrected errors in C***ROUTINES CALLED list.  (FNF)
C   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
C   920511  Added complete declaration section.  (WRB)
C   920929  Corrected format of references.  (FNF)
C   921019  Changed 500.0 to 500 to reduce SP/DP differences.  (FNF)
C   921026  Added check for valid value of ITOL.  (FNF)
C***END PROLOGUE  SGMRES
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      REAL ERR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LIGW, LRGW, N, NELT
C     .. Array Arguments ..
      REAL A(NELT), B(N), RGWK(LRGW), RWORK(*), SB(N), SX(N), X(N)
      INTEGER IA(NELT), IGWK(LIGW), IWORK(*), JA(NELT)
C     .. Subroutine Arguments ..
      EXTERNAL MATVEC, MSOLVE
C     .. Local Scalars ..
      REAL BNRM, RHOL, SUM
      INTEGER I, IFLAG, JPRE, JSCAL, KMP, LDL, LGMR, LHES, LQ, LR, LV,
     +        LW, LXL, LZ, LZM1, MAXL, MAXLP1, NMS, NMSL, NRMAX, NRSTS
C     .. External Functions ..
      REAL R1MACH, SNRM2
      EXTERNAL R1MACH, SNRM2
C     .. External Subroutines ..
      EXTERNAL SCOPY, SPIGMR
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C***FIRST EXECUTABLE STATEMENT  SGMRES
      IERR = 0
C   ------------------------------------------------------------------
C         Load method parameters with user values or defaults.
C   ------------------------------------------------------------------
      MAXL = IGWK(1)
      IF (MAXL .EQ. 0) MAXL = 10
      IF (MAXL .GT. N) MAXL = N
      KMP = IGWK(2)
      IF (KMP .EQ. 0) KMP = MAXL
      IF (KMP .GT. MAXL) KMP = MAXL
      JSCAL = IGWK(3)
      JPRE = IGWK(4)
C         Check for valid value of ITOL.
      IF( (ITOL.LT.0) .OR. ((ITOL.GT.3).AND.(ITOL.NE.11)) ) GOTO 650
C         Check for consistent values of ITOL and JPRE.
      IF( ITOL.EQ.1 .AND. JPRE.LT.0 ) GOTO 650
      IF( ITOL.EQ.2 .AND. JPRE.GE.0 ) GOTO 650
      NRMAX = IGWK(5)
      IF( NRMAX.EQ.0 ) NRMAX = 10
C         If NRMAX .eq. -1, then set NRMAX = 0 to turn off restarting.
      IF( NRMAX.EQ.-1 ) NRMAX = 0
C         If input value of TOL is zero, set it to its default value.
      IF( TOL.EQ.0.0E0 ) TOL = 500*R1MACH(3)
C
C         Initialize counters.
      ITER = 0
      NMS = 0
      NRSTS = 0
C   ------------------------------------------------------------------
C         Form work array segment pointers.
C   ------------------------------------------------------------------
      MAXLP1 = MAXL + 1
      LV = 1
      LR = LV + N*MAXLP1
      LHES = LR + N + 1
      LQ = LHES + MAXL*MAXLP1
      LDL = LQ + 2*MAXL
      LW = LDL + N
      LXL = LW + N
      LZ = LXL + N
C
C         Load IGWK(6) with required minimum length of the RGWK array.
      IGWK(6) = LZ + N - 1
      IF( LZ+N-1.GT.LRGW ) GOTO 640
C   ------------------------------------------------------------------
C         Calculate scaled-preconditioned norm of RHS vector b.
C   ------------------------------------------------------------------
      IF (JPRE .LT. 0) THEN
         CALL MSOLVE(N, B, RGWK(LR), NELT, IA, JA, A, ISYM,
     $        RWORK, IWORK)
         NMS = NMS + 1
      ELSE
         CALL SCOPY(N, B, 1, RGWK(LR), 1)
      ENDIF
      IF( JSCAL.EQ.2 .OR. JSCAL.EQ.3 ) THEN
         SUM = 0
         DO 10 I = 1,N
            SUM = SUM + (RGWK(LR-1+I)*SB(I))**2
 10      CONTINUE
         BNRM = SQRT(SUM)
      ELSE
         BNRM = SNRM2(N,RGWK(LR),1)
      ENDIF
C   ------------------------------------------------------------------
C         Calculate initial residual.
C   ------------------------------------------------------------------
      CALL MATVEC(N, X, RGWK(LR), NELT, IA, JA, A, ISYM)
      DO 50 I = 1,N
         RGWK(LR-1+I) = B(I) - RGWK(LR-1+I)
 50   CONTINUE
C   ------------------------------------------------------------------
C         If performing restarting, then load the residual into the
C         correct location in the RGWK array.
C   ------------------------------------------------------------------
 100  CONTINUE
      IF( NRSTS.GT.NRMAX ) GOTO 610
      IF( NRSTS.GT.0 ) THEN
C         Copy the current residual to a different location in the RGWK
C         array.
         CALL SCOPY(N, RGWK(LDL), 1, RGWK(LR), 1)
      ENDIF
C   ------------------------------------------------------------------
C         Use the SPIGMR algorithm to solve the linear system A*Z = R.
C   ------------------------------------------------------------------
      CALL SPIGMR(N, RGWK(LR), SB, SX, JSCAL, MAXL, MAXLP1, KMP,
     $       NRSTS, JPRE, MATVEC, MSOLVE, NMSL, RGWK(LZ), RGWK(LV),
     $       RGWK(LHES), RGWK(LQ), LGMR, RWORK, IWORK, RGWK(LW),
     $       RGWK(LDL), RHOL, NRMAX, B, BNRM, X, RGWK(LXL), ITOL,
     $       TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)
      ITER = ITER + LGMR
      NMS = NMS + NMSL
C
C         Increment X by the current approximate solution Z of A*Z = R.
C
      LZM1 = LZ - 1
      DO 110 I = 1,N
         X(I) = X(I) + RGWK(LZM1+I)
 110  CONTINUE
      IF( IFLAG.EQ.0 ) GOTO 600
      IF( IFLAG.EQ.1 ) THEN
         NRSTS = NRSTS + 1
         GOTO 100
      ENDIF
      IF( IFLAG.EQ.2 ) GOTO 620
C   ------------------------------------------------------------------
C         All returns are made through this section.
C   ------------------------------------------------------------------
C         The iteration has converged.
C
 600  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 0
      RETURN
C
C         Max number((NRMAX+1)*MAXL) of linear iterations performed.
 610  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 1
      RETURN
C
C         GMRES failed to reduce last residual in MAXL iterations.
C         The iteration has stalled.
 620  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 2
      RETURN
C         Error return.  Insufficient length for RGWK array.
 640  CONTINUE
      ERR = TOL
      IERR = -1
      RETURN
C         Error return.  Inconsistent ITOL and JPRE values.
 650  CONTINUE
      ERR = TOL
      IERR = -2
      RETURN
C------------- LAST LINE OF SGMRES FOLLOWS ----------------------------
      END
