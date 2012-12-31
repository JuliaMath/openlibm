*DECK SHELS
      SUBROUTINE SHELS (A, LDA, N, Q, B)
C***BEGIN PROLOGUE  SHELS
C***SUBSIDIARY
C***PURPOSE  Internal routine for SGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      SINGLE PRECISION (SHELS-S, DHELS-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C        This routine is extracted from the LINPACK routine SGESL with
C        changes due to the fact that A is an upper Hessenberg matrix.
C
C        SHELS solves the least squares problem:
C
C                   MIN(B-A*X,B-A*X)
C
C        using the factors computed by SHEQR.
C
C *Usage:
C      INTEGER LDA, N
C      REAL A(LDA,N), Q(2*N), B(N+1)
C
C      CALL SHELS(A, LDA, N, Q, B)
C
C *Arguments:
C A       :IN       Real A(LDA,N)
C          The output from SHEQR which contains the upper
C          triangular factor R in the QR decomposition of A.
C LDA     :IN       Integer
C          The leading dimension of the array A.
C N       :IN       Integer
C          A is originally an (N+1) by N matrix.
C Q       :IN       Real Q(2*N)
C          The coefficients of the N Givens rotations
C          used in the QR factorization of A.
C B       :INOUT    Real B(N+1)
C          On input, B is the right hand side vector.
C          On output, B is the solution vector X.
C
C***SEE ALSO  SGMRES
C***ROUTINES CALLED  SAXPY
C***REVISION HISTORY  (YYMMDD)
C   871001  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910502  Added C***FIRST EXECUTABLE STATEMENT line.  (FNF)
C   910506  Made subsidiary to SGMRES.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  SHELS
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      INTEGER LDA, N
C     .. Array Arguments ..
      REAL A(LDA,*), B(*), Q(*)
C     .. Local Scalars ..
      REAL C, S, T, T1, T2
      INTEGER IQ, K, KB, KP1
C     .. External Subroutines ..
      EXTERNAL SAXPY
C***FIRST EXECUTABLE STATEMENT  SHELS
C
C         Minimize(B-A*X,B-A*X).  First form Q*B.
C
      DO 20 K = 1, N
         KP1 = K + 1
         IQ = 2*(K-1) + 1
         C = Q(IQ)
         S = Q(IQ+1)
         T1 = B(K)
         T2 = B(KP1)
         B(K) = C*T1 - S*T2
         B(KP1) = S*T1 + C*T2
 20   CONTINUE
C
C         Now solve  R*X = Q*B.
C
      DO 40 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T = -B(K)
         CALL SAXPY(K-1, T, A(1,K), 1, B(1), 1)
 40   CONTINUE
      RETURN
C------------- LAST LINE OF SHELS FOLLOWS ----------------------------
      END
