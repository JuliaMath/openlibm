*DECK DHELS
      SUBROUTINE DHELS (A, LDA, N, Q, B)
C***BEGIN PROLOGUE  DHELS
C***SUBSIDIARY
C***PURPOSE  Internal routine for DGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SHELS-S, DHELS-D)
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
C        DHELS solves the least squares problem:
C
C                   MIN(B-A*X,B-A*X)
C
C        using the factors computed by DHEQR.
C
C *Usage:
C      INTEGER LDA, N
C      DOUBLE PRECISION A(LDA,N), Q(2*N), B(N+1)
C
C      CALL DHELS(A, LDA, N, Q, B)
C
C *Arguments:
C A       :IN       Double Precision A(LDA,N)
C          The output from DHEQR which contains the upper
C          triangular factor R in the QR decomposition of A.
C LDA     :IN       Integer
C          The leading dimension of the array A.
C N       :IN       Integer
C          A is originally an (N+1) by N matrix.
C Q       :IN       Double Precision Q(2*N)
C          The coefficients of the N Givens rotations
C          used in the QR factorization of A.
C B       :INOUT    Double Precision B(N+1)
C          On input, B is the right hand side vector.
C          On output, B is the solution vector X.
C
C***SEE ALSO  DGMRES
C***ROUTINES CALLED  DAXPY
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910502  Added C***FIRST EXECUTABLE STATEMENT line.  (FNF)
C   910506  Made subsidiary to DGMRES.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  DHELS
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      INTEGER LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), B(*), Q(*)
C     .. Local Scalars ..
      DOUBLE PRECISION C, S, T, T1, T2
      INTEGER IQ, K, KB, KP1
C     .. External Subroutines ..
      EXTERNAL DAXPY
C***FIRST EXECUTABLE STATEMENT  DHELS
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
         CALL DAXPY(K-1, T, A(1,K), 1, B(1), 1)
 40   CONTINUE
      RETURN
C------------- LAST LINE OF DHELS FOLLOWS ----------------------------
      END
