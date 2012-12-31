*DECK ISSGMR
      INTEGER FUNCTION ISSGMR (N, B, X, XL, NELT, IA, JA, A, ISYM,
     +   MSOLVE, NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ,
     +   RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL, KMP, LGMR, MAXL,
     +   MAXLP1, V, Q, SNORMW, PROD, R0NRM, HES, JPRE)
C***BEGIN PROLOGUE  ISSGMR
C***SUBSIDIARY
C***PURPOSE  Generalized Minimum Residual Stop Test.
C            This routine calculates the stop test for the Generalized
C            Minimum RESidual (GMRES) iteration scheme.  It returns a
C            non-zero if the error estimate (the type of which is
C            determined by ITOL) is less than the user specified
C            tolerance TOL.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      SINGLE PRECISION (ISSGMR-S, ISDGMR-D)
C***KEYWORDS  GMRES, LINEAR SYSTEM, SLAP, SPARSE, STOP TEST
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C
C *Usage:
C      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, NMSL, ITOL
C      INTEGER ITMAX, ITER, IUNIT, IWORK(USER DEFINED), JSCAL
C      INTEGER KMP, LGMR, MAXL, MAXLP1, JPRE
C      REAL B(N), X(N), XL(MAXL), A(NELT), TOL, ERR, R(N), Z(N),
C     $     DZ(N), RWORK(USER DEFINED), RNRM, BNRM, SB(N), SX(N),
C     $     V(N,MAXLP1), Q(2*MAXL), SNORMW, PROD, R0NRM,
C     $     HES(MAXLP1,MAXL)
C      EXTERNAL MSOLVE
C
C      IF (ISSGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,
C     $     NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ,
C     $     RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL,
C     $     KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,
C     $     HES, JPRE) .NE. 0) THEN ITERATION DONE
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right-hand-side vector.
C X      :IN       Real X(N).
C         Approximate solution vector as of the last restart.
C XL     :OUT      Real XL(N)
C         An array of length N used to hold the approximate
C         solution as of the current iteration.  Only computed by
C         this routine when ITOL=11.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Real A(NELT).
C         These arrays contain the matrix data structure for A.
C         It could take any form.  See "Description", in the SGMRES,
C         SSLUGM and SSDGMR routines for more details.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C MSOLVE :EXT      External.
C         Name of a routine which solves a linear system Mz = r for  z
C         given r with the preconditioning matrix M (M is supplied via
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
C NMSL   :INOUT    Integer.
C         A counter for the number of calls to MSOLVE.
C ITOL   :IN       Integer.
C         Flag to indicate the type of convergence criterion used.
C         ITOL=0  Means the  iteration stops when the test described
C                 below on  the  residual RL  is satisfied.  This is
C                 the  "Natural Stopping Criteria" for this routine.
C                 Other values  of   ITOL  cause  extra,   otherwise
C                 unnecessary, computation per iteration and     are
C                 therefore much less efficient.
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
C TOL    :IN       Real.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :IN       Integer.
C         The iteration for which to check for convergence.
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
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration,
C         if this is desired for monitoring convergence.  If unit
C         number is 0, no writing will occur.
C R      :INOUT    Real R(N).
C         Work array used in calling routine.  It contains
C         information necessary to compute the residual RL = B-A*XL.
C Z      :WORK     Real Z(N).
C         Workspace used to hold the pseudo-residual M z = r.
C DZ     :WORK     Real DZ(N).
C         Workspace used to hold temporary vector(s).
C RWORK  :WORK     Real RWORK(USER DEFINED).
C         Real array that can be used by MSOLVE.
C IWORK  :WORK     Integer IWORK(USER DEFINED).
C         Integer array that can be used by MSOLVE.
C RNRM   :IN       Real.
C         Norm of the current residual.  Type of norm depends on ITOL.
C BNRM   :IN       Real.
C         Norm of the right hand side.  Type of norm depends on ITOL.
C SB     :IN       Real SB(N).
C         Scaling vector for B.
C SX     :IN       Real SX(N).
C         Scaling vector for X.
C JSCAL  :IN       Integer.
C         Flag indicating if scaling arrays SB and SX are being
C         used in the calling routine SPIGMR.
C         JSCAL=0 means SB and SX are not used and the
C                 algorithm will perform as if all
C                 SB(i) = 1 and SX(i) = 1.
C         JSCAL=1 means only SX is used, and the algorithm
C                 performs as if all SB(i) = 1.
C         JSCAL=2 means only SB is used, and the algorithm
C                 performs as if all SX(i) = 1.
C         JSCAL=3 means both SB and SX are used.
C KMP    :IN       Integer
C         The number of previous vectors the new vector VNEW
C         must be made orthogonal to.  (KMP .le. MAXL)
C LGMR   :IN       Integer
C         The number of GMRES iterations performed on the current call
C         to SPIGMR (i.e., # iterations since the last restart) and
C         the current order of the upper Hessenberg
C         matrix HES.
C MAXL   :IN       Integer
C         The maximum allowable order of the matrix H.
C MAXLP1 :IN       Integer
C         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES.
C V      :IN       Real V(N,MAXLP1)
C         The N by (LGMR+1) array containing the LGMR
C         orthogonal vectors V(*,1) to V(*,LGMR).
C Q      :IN       Real Q(2*MAXL)
C         A real array of length 2*MAXL containing the components
C         of the Givens rotations used in the QR decomposition
C         of HES.
C SNORMW :IN       Real
C         A scalar containing the scaled norm of VNEW before it
C         is renormalized in SPIGMR.
C PROD   :IN       Real
C         The product s1*s2*...*sl = the product of the sines of the
C         Givens rotations used in the QR factorization of the
C         Hessenberg matrix HES.
C R0NRM  :IN       Real
C         The scaled norm of initial residual R0.
C HES    :IN       Real HES(MAXLP1,MAXL)
C         The upper triangular factor of the QR decomposition
C         of the (LGMR+1) by LGMR upper Hessenberg matrix whose
C         entries are the scaled inner-products of A*V(*,I)
C         and V(*,K).
C JPRE   :IN       Integer
C         Preconditioner type flag.
C         (See description of IGWK(4) in SGMRES.)
C
C *Description
C       When using the GMRES solver,  the preferred value  for ITOL
C       is 0.  This is due to the fact that when ITOL=0 the norm of
C       the residual required in the stopping test is  obtained for
C       free, since this value is already  calculated  in the GMRES
C       algorithm.   The  variable  RNRM contains the   appropriate
C       norm, which is equal to norm(SB*(RL - A*XL))  when right or
C       no   preconditioning is  being  performed,   and equal   to
C       norm(SB*Minv*(RL - A*XL))  when using left preconditioning.
C       Here, norm() is the Euclidean norm.  Nonzero values of ITOL
C       require  additional work  to  calculate the  actual  scaled
C       residual  or its scaled/preconditioned  form,  and/or   the
C       approximate solution XL.  Hence, these values of  ITOL will
C       not be as efficient as ITOL=0.
C
C *Cautions:
C     This routine will attempt to write to the Fortran logical output
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
C     this logical unit is attached to a file or terminal before calling
C     this routine with a non-zero value for IUNIT.  This routine does
C     not check for the validity of a non-zero IUNIT unit number.
C
C     This routine does not verify that ITOL has a valid value.
C     The calling routine should make such a test before calling
C     ISSGMR, as is done in SGMRES.
C
C***SEE ALSO  SGMRES
C***ROUTINES CALLED  R1MACH, SCOPY, SNRM2, SRLCAL, SSCAL, SXLCAL
C***COMMON BLOCKS    SSLBLK
C***REVISION HISTORY  (YYMMDD)
C   871211  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910502  Corrected conversion errors, etc.  (FNF)
C   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
C   910506  Made subsidiary to SGMRES.  (FNF)
C   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
C   920511  Added complete declaration section.  (WRB)
C   921113  Corrected C***CATEGORY line.  (FNF)
C***END PROLOGUE  ISSGMR
C     .. Scalar Arguments ..
      REAL BNRM, ERR, PROD, R0NRM, RNRM, SNORMW, TOL
      INTEGER ISYM, ITER, ITMAX, ITOL, IUNIT, JPRE, JSCAL, KMP, LGMR,
     +        MAXL, MAXLP1, N, NELT, NMSL
C     .. Array Arguments ..
      REAL A(*), B(*), DZ(*), HES(MAXLP1, MAXL), Q(*), R(*), RWORK(*),
     +     SB(*), SX(*), V(N,*), X(*), XL(*), Z(*)
      INTEGER IA(*), IWORK(*), JA(*)
C     .. Subroutine Arguments ..
      EXTERNAL MSOLVE
C     .. Arrays in Common ..
      REAL SOLN(1)
C     .. Local Scalars ..
      REAL DXNRM, FUZZ, RAT, RATMAX, SOLNRM, TEM
      INTEGER I, IELMAX
C     .. External Functions ..
      REAL R1MACH, SNRM2
      EXTERNAL R1MACH, SNRM2
C     .. External Subroutines ..
      EXTERNAL SCOPY, SRLCAL, SSCAL, SXLCAL
C     .. Intrinsic Functions ..
      INTRINSIC ABS, MAX, SQRT
C     .. Common blocks ..
      COMMON /SSLBLK/ SOLN
C     .. Save statement ..
      SAVE SOLNRM
C***FIRST EXECUTABLE STATEMENT  ISSGMR
      ISSGMR = 0
      IF ( ITOL.EQ.0 ) THEN
C
C       Use input from SPIGMR to determine if stop conditions are met.
C
         ERR = RNRM/BNRM
      ENDIF
      IF ( (ITOL.GT.0) .AND. (ITOL.LE.3) ) THEN
C
C       Use SRLCAL to calculate the scaled residual vector.
C       Store answer in R.
C
         IF ( LGMR.NE.0 ) CALL SRLCAL(N, KMP, LGMR, MAXL, V, Q, R,
     $                                SNORMW, PROD, R0NRM)
         IF ( ITOL.LE.2 ) THEN
C         err = ||Residual||/||RightHandSide||(2-Norms).
            ERR = SNRM2(N, R, 1)/BNRM
C
C         Unscale R by R0NRM*PROD when KMP < MAXL.
C
            IF ( (KMP.LT.MAXL) .AND. (LGMR.NE.0) ) THEN
               TEM = 1.0E0/(R0NRM*PROD)
               CALL SSCAL(N, TEM, R, 1)
            ENDIF
         ELSEIF ( ITOL.EQ.3 ) THEN
C         err = Max |(Minv*Residual)(i)/x(i)|
C         When JPRE .lt. 0, R already contains Minv*Residual.
            IF ( JPRE.GT.0 ) THEN
               CALL MSOLVE(N, R, DZ, NELT, IA, JA, A, ISYM, RWORK,
     $              IWORK)
               NMSL = NMSL + 1
            ENDIF
C
C         Unscale R by R0NRM*PROD when KMP < MAXL.
C
            IF ( (KMP.LT.MAXL) .AND. (LGMR.NE.0) ) THEN
               TEM = 1.0E0/(R0NRM*PROD)
               CALL SSCAL(N, TEM, R, 1)
            ENDIF
C
            FUZZ = R1MACH(1)
            IELMAX = 1
            RATMAX = ABS(DZ(1))/MAX(ABS(X(1)),FUZZ)
            DO 25 I = 2, N
               RAT = ABS(DZ(I))/MAX(ABS(X(I)),FUZZ)
               IF( RAT.GT.RATMAX ) THEN
                  IELMAX = I
                  RATMAX = RAT
               ENDIF
 25         CONTINUE
            ERR = RATMAX
            IF( RATMAX.LE.TOL ) ISSGMR = 1
            IF( IUNIT.GT.0 ) WRITE(IUNIT,1020) ITER, IELMAX, RATMAX
            RETURN
         ENDIF
      ENDIF
      IF ( ITOL.EQ.11 ) THEN
C
C       Use SXLCAL to calculate the approximate solution XL.
C
         IF ( (LGMR.NE.0) .AND. (ITER.GT.0) ) THEN
            CALL SXLCAL(N, LGMR, X, XL, XL, HES, MAXLP1, Q, V, R0NRM,
     $           DZ, SX, JSCAL, JPRE, MSOLVE, NMSL, RWORK, IWORK,
     $           NELT, IA, JA, A, ISYM)
         ELSEIF ( ITER.EQ.0 ) THEN
C         Copy X to XL to check if initial guess is good enough.
            CALL SCOPY(N, X, 1, XL, 1)
         ELSE
C         Return since this is the first call to SPIGMR on a restart.
            RETURN
         ENDIF
C
         IF ((JSCAL .EQ. 0) .OR.(JSCAL .EQ. 2)) THEN
C         err = ||x-TrueSolution||/||TrueSolution||(2-Norms).
            IF ( ITER.EQ.0 ) SOLNRM = SNRM2(N, SOLN, 1)
            DO 30 I = 1, N
               DZ(I) = XL(I) - SOLN(I)
 30         CONTINUE
            ERR = SNRM2(N, DZ, 1)/SOLNRM
         ELSE
            IF (ITER .EQ. 0) THEN
               SOLNRM = 0
               DO 40 I = 1,N
                  SOLNRM = SOLNRM + (SX(I)*SOLN(I))**2
 40            CONTINUE
               SOLNRM = SQRT(SOLNRM)
            ENDIF
            DXNRM = 0
            DO 50 I = 1,N
               DXNRM = DXNRM + (SX(I)*(XL(I)-SOLN(I)))**2
 50         CONTINUE
            DXNRM = SQRT(DXNRM)
C         err = ||SX*(x-TrueSolution)||/||SX*TrueSolution|| (2-Norms).
            ERR = DXNRM/SOLNRM
         ENDIF
      ENDIF
C
      IF( IUNIT.NE.0 ) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL, MAXL, KMP
         ENDIF
         WRITE(IUNIT,1010) ITER, RNRM/BNRM, ERR
      ENDIF
      IF ( ERR.LE.TOL ) ISSGMR = 1
C
      RETURN
 1000 FORMAT(' Generalized Minimum Residual(',I3,I3,') for ',
     $     'N, ITOL = ',I5, I5,
     $     /' ITER','   Natural Err Est','   Error Estimate')
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7)
 1020 FORMAT(1X,' ITER = ',I5, ' IELMAX = ',I5,
     $     ' |R(IELMAX)/X(IELMAX)| = ',E12.5)
C------------- LAST LINE OF ISSGMR FOLLOWS ----------------------------
      END
