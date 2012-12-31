*DECK DPIGMR
      SUBROUTINE DPIGMR (N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP, NRSTS,
     +   JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR, RPAR, IPAR, WK,
     +   DL, RHOL, NRMAX, B, BNRM, X, XL, ITOL, TOL, NELT, IA, JA, A,
     +   ISYM, IUNIT, IFLAG, ERR)
C***BEGIN PROLOGUE  DPIGMR
C***SUBSIDIARY
C***PURPOSE  Internal routine for DGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SPIGMR-S, DPIGMR-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C         This routine solves the linear system A * Z = R0 using a
C         scaled preconditioned version of the generalized minimum
C         residual method.  An initial guess of Z = 0 is assumed.
C
C *Usage:
C      INTEGER N, JSCAL, MAXL, MAXLP1, KMP, NRSTS, JPRE, NMSL, LGMR
C      INTEGER IPAR(USER DEFINED), NRMAX, ITOL, NELT, IA(NELT), JA(NELT)
C      INTEGER ISYM, IUNIT, IFLAG
C      DOUBLE PRECISION R0(N), SR(N), SZ(N), Z(N), V(N,MAXLP1),
C     $                 HES(MAXLP1,MAXL), Q(2*MAXL), RPAR(USER DEFINED),
C     $                 WK(N), DL(N), RHOL, B(N), BNRM, X(N), XL(N),
C     $                 TOL, A(NELT), ERR
C      EXTERNAL MATVEC, MSOLVE
C
C      CALL DPIGMR(N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP,
C     $     NRSTS, JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR,
C     $     RPAR, IPAR, WK, DL, RHOL, NRMAX, B, BNRM, X, XL,
C     $     ITOL, TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)
C
C *Arguments:
C N      :IN       Integer
C         The order of the matrix A, and the lengths
C         of the vectors SR, SZ, R0 and Z.
C R0     :IN       Double Precision R0(N)
C         R0 = the right hand side of the system A*Z = R0.
C         R0 is also used as workspace when computing
C         the final approximation.
C         (R0 is the same as V(*,MAXL+1) in the call to DPIGMR.)
C SR     :IN       Double Precision SR(N)
C         SR is a vector of length N containing the non-zero
C         elements of the diagonal scaling matrix for R0.
C SZ     :IN       Double Precision SZ(N)
C         SZ is a vector of length N containing the non-zero
C         elements of the diagonal scaling matrix for Z.
C JSCAL  :IN       Integer
C         A flag indicating whether arrays SR and SZ are used.
C         JSCAL=0 means SR and SZ are not used and the
C                 algorithm will perform as if all
C                 SR(i) = 1 and SZ(i) = 1.
C         JSCAL=1 means only SZ is used, and the algorithm
C                 performs as if all SR(i) = 1.
C         JSCAL=2 means only SR is used, and the algorithm
C                 performs as if all SZ(i) = 1.
C         JSCAL=3 means both SR and SZ are used.
C MAXL   :IN       Integer
C         The maximum allowable order of the matrix H.
C MAXLP1 :IN       Integer
C         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES.
C KMP    :IN       Integer
C         The number of previous vectors the new vector VNEW
C         must be made orthogonal to.  (KMP .le. MAXL)
C NRSTS  :IN       Integer
C         Counter for the number of restarts on the current
C         call to DGMRES.  If NRSTS .gt. 0, then the residual
C         R0 is already scaled, and so scaling of it is
C         not necessary.
C JPRE   :IN       Integer
C         Preconditioner type flag.
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
C         RPAR and IPAR arrays.  The name of the MSOLVE routine must
C         be declared external in the calling program.  The calling
C         sequence to MSOLVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
C         Where N is the number of unknowns, R is the right-hand side
C         vector and Z is the solution upon return.  NELT, IA, JA, A and
C         ISYM are defined as below.  RPAR is a double precision array
C         that can be used to pass necessary preconditioning information
C         and/or workspace to MSOLVE.  IPAR is an integer work array
C         for the same purpose as RPAR.
C NMSL   :OUT      Integer
C         The number of calls to MSOLVE.
C Z      :OUT      Double Precision Z(N)
C         The final computed approximation to the solution
C         of the system A*Z = R0.
C V      :OUT      Double Precision V(N,MAXLP1)
C         The N by (LGMR+1) array containing the LGMR
C         orthogonal vectors V(*,1) to V(*,LGMR).
C HES    :OUT      Double Precision HES(MAXLP1,MAXL)
C         The upper triangular factor of the QR decomposition
C         of the (LGMR+1) by LGMR upper Hessenberg matrix whose
C         entries are the scaled inner-products of A*V(*,I)
C         and V(*,K).
C Q      :OUT      Double Precision Q(2*MAXL)
C         A double precision array of length 2*MAXL containing the
C         components of the Givens rotations used in the QR
C         decomposition of HES.  It is loaded in DHEQR and used in
C         DHELS.
C LGMR   :OUT      Integer
C         The number of iterations performed and
C         the current order of the upper Hessenberg
C         matrix HES.
C RPAR   :IN       Double Precision RPAR(USER DEFINED)
C         Double Precision workspace passed directly to the MSOLVE
C         routine.
C IPAR   :IN       Integer IPAR(USER DEFINED)
C         Integer workspace passed directly to the MSOLVE routine.
C WK     :IN       Double Precision WK(N)
C         A double precision work array of length N used by routines
C         MATVEC and MSOLVE.
C DL     :INOUT    Double Precision DL(N)
C         On input, a double precision work array of length N used for
C         calculation of the residual norm RHO when the method is
C         incomplete (KMP.lt.MAXL), and/or when using restarting.
C         On output, the scaled residual vector RL.  It is only loaded
C         when performing restarts of the Krylov iteration.
C RHOL   :OUT      Double Precision
C         A double precision scalar containing the norm of the final
C         residual.
C NRMAX  :IN       Integer
C         The maximum number of restarts of the Krylov iteration.
C         NRMAX .gt. 0 means restarting is active, while
C         NRMAX = 0 means restarting is not being used.
C B      :IN       Double Precision B(N)
C         The right hand side of the linear system A*X = b.
C BNRM   :IN       Double Precision
C         The scaled norm of b.
C X      :IN       Double Precision X(N)
C         The current approximate solution as of the last
C         restart.
C XL     :IN       Double Precision XL(N)
C         An array of length N used to hold the approximate
C         solution X(L) when ITOL=11.
C ITOL   :IN       Integer
C         A flag to indicate the type of convergence criterion
C         used.  See the driver for its description.
C TOL    :IN       Double Precision
C         The tolerance on residuals R0-A*Z in scaled norm.
C NELT   :IN       Integer
C         The length of arrays IA, JA and A.
C IA     :IN       Integer IA(NELT)
C         An integer array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C JA     :IN       Integer JA(NELT)
C         An integer array of length NELT containing matrix data.
C         It is passed directly to the MATVEC and MSOLVE routines.
C A      :IN       Double Precision A(NELT)
C         A double precision array of length NELT containing matrix
C         data. It is passed directly to the MATVEC and MSOLVE routines.
C ISYM   :IN       Integer
C         A flag to indicate symmetric matrix storage.
C         If ISYM=0, all non-zero entries of the matrix are
C         stored.  If ISYM=1, the matrix is symmetric and
C         only the upper or lower triangular part is stored.
C IUNIT  :IN       Integer
C         The i/o unit number for writing intermediate residual
C         norm values.
C IFLAG  :OUT      Integer
C         An integer error flag..
C         0 means convergence in LGMR iterations, LGMR.le.MAXL.
C         1 means the convergence test did not pass in MAXL
C           iterations, but the residual norm is .lt. norm(R0),
C           and so Z is computed.
C         2 means the convergence test did not pass in MAXL
C           iterations, residual .ge. norm(R0), and Z = 0.
C ERR    :OUT      Double Precision.
C         Error estimate of error in final approximate solution, as
C         defined by ITOL.
C
C *Cautions:
C     This routine will attempt to write to the Fortran logical output
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
C     this logical unit is attached to a file or terminal before calling
C     this routine with a non-zero value for IUNIT.  This routine does
C     not check for the validity of a non-zero IUNIT unit number.
C
C***SEE ALSO  DGMRES
C***ROUTINES CALLED  DAXPY, DCOPY, DHELS, DHEQR, DNRM2, DORTH, DRLCAL,
C                    DSCAL, ISDGMR
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910502  Removed MATVEC and MSOLVE from ROUTINES CALLED list.  (FNF)
C   910506  Made subsidiary to DGMRES.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  DPIGMR
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      DOUBLE PRECISION BNRM, ERR, RHOL, TOL
      INTEGER IFLAG, ISYM, ITOL, IUNIT, JPRE, JSCAL, KMP, LGMR, MAXL,
     +        MAXLP1, N, NELT, NMSL, NRMAX, NRSTS
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(*), DL(*), HES(MAXLP1,*), Q(*), R0(*),
     +                 RPAR(*), SR(*), SZ(*), V(N,*), WK(*), X(*),
     +                 XL(*), Z(*)
      INTEGER IA(NELT), IPAR(*), JA(NELT)
C     .. Subroutine Arguments ..
      EXTERNAL MATVEC, MSOLVE
C     .. Local Scalars ..
      DOUBLE PRECISION C, DLNRM, PROD, R0NRM, RHO, S, SNORMW, TEM
      INTEGER I, I2, INFO, IP1, ITER, ITMAX, J, K, LL, LLP1
C     .. External Functions ..
      DOUBLE PRECISION DNRM2
      INTEGER ISDGMR
      EXTERNAL DNRM2, ISDGMR
C     .. External Subroutines ..
      EXTERNAL DAXPY, DCOPY, DHELS, DHEQR, DORTH, DRLCAL, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C***FIRST EXECUTABLE STATEMENT  DPIGMR
C
C         Zero out the Z array.
C
      DO 5 I = 1,N
         Z(I) = 0
 5    CONTINUE
C
      IFLAG = 0
      LGMR = 0
      NMSL = 0
C         Load ITMAX, the maximum number of iterations.
      ITMAX =(NRMAX+1)*MAXL
C   -------------------------------------------------------------------
C         The initial residual is the vector R0.
C         Apply left precon. if JPRE < 0 and this is not a restart.
C         Apply scaling to R0 if JSCAL = 2 or 3.
C   -------------------------------------------------------------------
      IF ((JPRE .LT. 0) .AND.(NRSTS .EQ. 0)) THEN
         CALL DCOPY(N, R0, 1, WK, 1)
         CALL MSOLVE(N, WK, R0, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
      IF (((JSCAL.EQ.2) .OR.(JSCAL.EQ.3)) .AND.(NRSTS.EQ.0)) THEN
         DO 10 I = 1,N
            V(I,1) = R0(I)*SR(I)
 10      CONTINUE
      ELSE
         DO 20 I = 1,N
            V(I,1) = R0(I)
 20      CONTINUE
      ENDIF
      R0NRM = DNRM2(N, V, 1)
      ITER = NRSTS*MAXL
C
C         Call stopping routine ISDGMR.
C
      IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,
     $    NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, V(1,1), Z, WK,
     $    RPAR, IPAR, R0NRM, BNRM, SR, SZ, JSCAL,
     $    KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,
     $    HES, JPRE) .NE. 0) RETURN
      TEM = 1.0D0/R0NRM
      CALL DSCAL(N, TEM, V(1,1), 1)
C
C         Zero out the HES array.
C
      DO 50 J = 1,MAXL
         DO 40 I = 1,MAXLP1
            HES(I,J) = 0
 40      CONTINUE
 50   CONTINUE
C   -------------------------------------------------------------------
C         Main loop to compute the vectors V(*,2) to V(*,MAXL).
C         The running product PROD is needed for the convergence test.
C   -------------------------------------------------------------------
      PROD = 1
      DO 90 LL = 1,MAXL
         LGMR = LL
C   -------------------------------------------------------------------
C        Unscale  the  current V(LL)  and store  in WK.  Call routine
C        MSOLVE    to   compute(M-inverse)*WK,   where    M   is  the
C        preconditioner matrix.  Save the answer in Z.   Call routine
C        MATVEC to compute  VNEW  = A*Z,  where  A is  the the system
C        matrix.  save the answer in  V(LL+1).  Scale V(LL+1).   Call
C        routine DORTH  to  orthogonalize the    new vector VNEW   =
C        V(*,LL+1).  Call routine DHEQR to update the factors of HES.
C   -------------------------------------------------------------------
        IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
           DO 60 I = 1,N
              WK(I) = V(I,LL)/SZ(I)
 60        CONTINUE
        ELSE
           CALL DCOPY(N, V(1,LL), 1, WK, 1)
        ENDIF
        IF (JPRE .GT. 0) THEN
           CALL MSOLVE(N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
           NMSL = NMSL + 1
           CALL MATVEC(N, Z, V(1,LL+1), NELT, IA, JA, A, ISYM)
        ELSE
           CALL MATVEC(N, WK, V(1,LL+1), NELT, IA, JA, A, ISYM)
        ENDIF
        IF (JPRE .LT. 0) THEN
           CALL DCOPY(N, V(1,LL+1), 1, WK, 1)
           CALL MSOLVE(N,WK,V(1,LL+1),NELT,IA,JA,A,ISYM,RPAR,IPAR)
           NMSL = NMSL + 1
        ENDIF
        IF ((JSCAL .EQ. 2) .OR.(JSCAL .EQ. 3)) THEN
           DO 65 I = 1,N
              V(I,LL+1) = V(I,LL+1)*SR(I)
 65        CONTINUE
        ENDIF
        CALL DORTH(V(1,LL+1), V, HES, N, LL, MAXLP1, KMP, SNORMW)
        HES(LL+1,LL) = SNORMW
        CALL DHEQR(HES, MAXLP1, LL, Q, INFO, LL)
        IF (INFO .EQ. LL) GO TO 120
C   -------------------------------------------------------------------
C         Update RHO, the estimate of the norm of the residual R0-A*ZL.
C         If KMP <  MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
C         necessarily orthogonal for LL > KMP.  The vector DL must then
C         be computed, and its norm used in the calculation of RHO.
C   -------------------------------------------------------------------
        PROD = PROD*Q(2*LL)
        RHO = ABS(PROD*R0NRM)
        IF ((LL.GT.KMP) .AND.(KMP.LT.MAXL)) THEN
           IF (LL .EQ. KMP+1) THEN
              CALL DCOPY(N, V(1,1), 1, DL, 1)
              DO 75 I = 1,KMP
                 IP1 = I + 1
                 I2 = I*2
                 S = Q(I2)
                 C = Q(I2-1)
                 DO 70 K = 1,N
                    DL(K) = S*DL(K) + C*V(K,IP1)
 70              CONTINUE
 75           CONTINUE
           ENDIF
           S = Q(2*LL)
           C = Q(2*LL-1)/SNORMW
           LLP1 = LL + 1
           DO 80 K = 1,N
              DL(K) = S*DL(K) + C*V(K,LLP1)
 80        CONTINUE
           DLNRM = DNRM2(N, DL, 1)
           RHO = RHO*DLNRM
        ENDIF
        RHOL = RHO
C   -------------------------------------------------------------------
C         Test for convergence.  If passed, compute approximation ZL.
C         If failed and LL < MAXL, then continue iterating.
C   -------------------------------------------------------------------
        ITER = NRSTS*MAXL + LGMR
        IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,
     $      NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, DL, Z, WK,
     $      RPAR, IPAR, RHOL, BNRM, SR, SZ, JSCAL,
     $      KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,
     $      HES, JPRE) .NE. 0) GO TO 200
        IF (LL .EQ. MAXL) GO TO 100
C   -------------------------------------------------------------------
C         Rescale so that the norm of V(1,LL+1) is one.
C   -------------------------------------------------------------------
        TEM = 1.0D0/SNORMW
        CALL DSCAL(N, TEM, V(1,LL+1), 1)
 90   CONTINUE
 100  CONTINUE
      IF (RHO .LT. R0NRM) GO TO 150
 120  CONTINUE
      IFLAG = 2
C
C         Load approximate solution with zero.
C
      DO 130 I = 1,N
         Z(I) = 0
 130  CONTINUE
      RETURN
 150  IFLAG = 1
C
C         Tolerance not met, but residual norm reduced.
C
      IF (NRMAX .GT. 0) THEN
C
C        If performing restarting (NRMAX > 0)  calculate the residual
C        vector RL and  store it in the DL  array.  If the incomplete
C        version is being used (KMP < MAXL) then DL has  already been
C        calculated up to a scaling factor.   Use DRLCAL to calculate
C        the scaled residual vector.
C
         CALL DRLCAL(N, KMP, MAXL, MAXL, V, Q, DL, SNORMW, PROD,
     $        R0NRM)
      ENDIF
C   -------------------------------------------------------------------
C         Compute the approximation ZL to the solution.  Since the
C         vector Z was used as workspace, and the initial guess
C         of the linear iteration is zero, Z must be reset to zero.
C   -------------------------------------------------------------------
 200  CONTINUE
      LL = LGMR
      LLP1 = LL + 1
      DO 210 K = 1,LLP1
         R0(K) = 0
 210  CONTINUE
      R0(1) = R0NRM
      CALL DHELS(HES, MAXLP1, LL, Q, R0)
      DO 220 K = 1,N
         Z(K) = 0
 220  CONTINUE
      DO 230 I = 1,LL
         CALL DAXPY(N, R0(I), V(1,I), 1, Z, 1)
 230  CONTINUE
      IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
         DO 240 I = 1,N
            Z(I) = Z(I)/SZ(I)
 240     CONTINUE
      ENDIF
      IF (JPRE .GT. 0) THEN
         CALL DCOPY(N, Z, 1, WK, 1)
         CALL MSOLVE(N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
      RETURN
C------------- LAST LINE OF DPIGMR FOLLOWS ----------------------------
      END
