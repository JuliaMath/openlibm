*DECK SRLCAL
      SUBROUTINE SRLCAL (N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD,
     +   R0NRM)
C***BEGIN PROLOGUE  SRLCAL
C***SUBSIDIARY
C***PURPOSE  Internal routine for SGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      SINGLE PRECISION (SRLCAL-S, DRLCAL-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C         This routine calculates the scaled residual RL from the
C         V(I)'s.
C *Usage:
C      INTEGER N, KMP, LL, MAXL
C      REAL V(N,LL), Q(2*MAXL), RL(N), SNORMW, PROD, R0NORM
C
C      CALL SRLCAL(N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD, R0NRM)
C
C *Arguments:
C N      :IN       Integer
C         The order of the matrix A, and the lengths
C         of the vectors SR, SZ, R0 and Z.
C KMP    :IN       Integer
C         The number of previous V vectors the new vector VNEW
C         must be made orthogonal to. (KMP .le. MAXL)
C LL     :IN       Integer
C         The current dimension of the Krylov subspace.
C MAXL   :IN       Integer
C         The maximum dimension of the Krylov subspace.
C V      :IN       Real V(N,LL)
C         The N x LL array containing the orthogonal vectors
C         V(*,1) to V(*,LL).
C Q      :IN       Real Q(2*MAXL)
C         A real array of length 2*MAXL containing the components
C         of the Givens rotations used in the QR decomposition
C         of HES.  It is loaded in SHEQR and used in SHELS.
C RL     :OUT      Real RL(N)
C         The residual vector RL.  This is either SB*(B-A*XL) if
C         not preconditioning or preconditioning on the right,
C         or SB*(M-inverse)*(B-A*XL) if preconditioning on the
C         left.
C SNORMW :IN       Real
C         Scale factor.
C PROD   :IN       Real
C         The product s1*s2*...*sl = the product of the sines of the
C         Givens rotations used in the QR factorization of
C         the Hessenberg matrix HES.
C R0NRM  :IN       Real
C         The scaled norm of initial residual R0.
C
C***SEE ALSO  SGMRES
C***ROUTINES CALLED  SCOPY, SSCAL
C***REVISION HISTORY  (YYMMDD)
C   871001  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910506  Made subsidiary to SGMRES.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  SRLCAL
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      REAL PROD, R0NRM, SNORMW
      INTEGER KMP, LL, MAXL, N
C     .. Array Arguments ..
      REAL Q(*), RL(N), V(N,*)
C     .. Local Scalars ..
      REAL C, S, TEM
      INTEGER I, I2, IP1, K, LLM1, LLP1
C     .. External Subroutines ..
      EXTERNAL SCOPY, SSCAL
C***FIRST EXECUTABLE STATEMENT  SRLCAL
      IF (KMP .EQ. MAXL) THEN
C
C         calculate RL.  Start by copying V(*,1) into RL.
C
         CALL SCOPY(N, V(1,1), 1, RL, 1)
         LLM1 = LL - 1
         DO 20 I = 1,LLM1
            IP1 = I + 1
            I2 = I*2
            S = Q(I2)
            C = Q(I2-1)
            DO 10 K = 1,N
               RL(K) = S*RL(K) + C*V(K,IP1)
 10         CONTINUE
 20      CONTINUE
         S = Q(2*LL)
         C = Q(2*LL-1)/SNORMW
         LLP1 = LL + 1
         DO 30 K = 1,N
            RL(K) = S*RL(K) + C*V(K,LLP1)
 30      CONTINUE
      ENDIF
C
C         When KMP < MAXL, RL vector already partially calculated.
C         Scale RL by R0NRM*PROD to obtain the residual RL.
C
      TEM = R0NRM*PROD
      CALL SSCAL(N, TEM, RL, 1)
      RETURN
C------------- LAST LINE OF SRLCAL FOLLOWS ----------------------------
      END
