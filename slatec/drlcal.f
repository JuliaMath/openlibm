*DECK DRLCAL
      SUBROUTINE DRLCAL (N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD,
     +   R0NRM)
C***BEGIN PROLOGUE  DRLCAL
C***SUBSIDIARY
C***PURPOSE  Internal routine for DGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SRLCAL-S, DRLCAL-D)
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
C      DOUBLE PRECISION V(N,LL), Q(2*MAXL), RL(N), SNORMW, PROD, R0NORM
C
C      CALL DRLCAL(N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD, R0NRM)
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
C V      :IN       Double Precision V(N,LL)
C         The N x LL array containing the orthogonal vectors
C         V(*,1) to V(*,LL).
C Q      :IN       Double Precision Q(2*MAXL)
C         A double precision array of length 2*MAXL containing the
C         components of the Givens rotations used in the QR
C         decomposition of HES.  It is loaded in DHEQR and used in
C         DHELS.
C RL     :OUT      Double Precision RL(N)
C         The residual vector RL.  This is either SB*(B-A*XL) if
C         not preconditioning or preconditioning on the right,
C         or SB*(M-inverse)*(B-A*XL) if preconditioning on the
C         left.
C SNORMW :IN       Double Precision
C         Scale factor.
C PROD   :IN       Double Precision
C         The product s1*s2*...*sl = the product of the sines of the
C         Givens rotations used in the QR factorization of
C         the Hessenberg matrix HES.
C R0NRM  :IN       Double Precision
C         The scaled norm of initial residual R0.
C
C***SEE ALSO  DGMRES
C***ROUTINES CALLED  DCOPY, DSCAL
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910506  Made subsidiary to DGMRES.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C***END PROLOGUE  DRLCAL
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      DOUBLE PRECISION PROD, R0NRM, SNORMW
      INTEGER KMP, LL, MAXL, N
C     .. Array Arguments ..
      DOUBLE PRECISION Q(*), RL(N), V(N,*)
C     .. Local Scalars ..
      DOUBLE PRECISION C, S, TEM
      INTEGER I, I2, IP1, K, LLM1, LLP1
C     .. External Subroutines ..
      EXTERNAL DCOPY, DSCAL
C***FIRST EXECUTABLE STATEMENT  DRLCAL
      IF (KMP .EQ. MAXL) THEN
C
C         calculate RL.  Start by copying V(*,1) into RL.
C
         CALL DCOPY(N, V(1,1), 1, RL, 1)
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
      CALL DSCAL(N, TEM, RL, 1)
      RETURN
C------------- LAST LINE OF DRLCAL FOLLOWS ----------------------------
      END
