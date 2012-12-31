*DECK DHEQR
      SUBROUTINE DHEQR (A, LDA, N, Q, INFO, IJOB)
C***BEGIN PROLOGUE  DHEQR
C***SUBSIDIARY
C***PURPOSE  Internal routine for DGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      DOUBLE PRECISION (SHEQR-S, DHEQR-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C        This   routine  performs  a QR   decomposition  of an  upper
C        Hessenberg matrix A using Givens  rotations.  There  are two
C        options  available: 1)  Performing  a fresh decomposition 2)
C        updating the QR factors by adding a row and  a column to the
C        matrix A.
C
C *Usage:
C      INTEGER LDA, N, INFO, IJOB
C      DOUBLE PRECISION A(LDA,N), Q(2*N)
C
C      CALL DHEQR(A, LDA, N, Q, INFO, IJOB)
C
C *Arguments:
C A      :INOUT    Double Precision A(LDA,N)
C         On input, the matrix to be decomposed.
C         On output, the upper triangular matrix R.
C         The factorization can be written Q*A = R, where
C         Q is a product of Givens rotations and R is upper
C         triangular.
C LDA    :IN       Integer
C         The leading dimension of the array A.
C N      :IN       Integer
C         A is an (N+1) by N Hessenberg matrix.
C Q      :OUT      Double Precision Q(2*N)
C         The factors c and s of each Givens rotation used
C         in decomposing A.
C INFO   :OUT      Integer
C         = 0  normal value.
C         = K  if  A(K,K) .eq. 0.0 .  This is not an error
C           condition for this subroutine, but it does
C           indicate that DHELS will divide by zero
C           if called.
C IJOB   :IN       Integer
C         = 1     means that a fresh decomposition of the
C                 matrix A is desired.
C         .ge. 2  means that the current decomposition of A
C                 will be updated by the addition of a row
C                 and a column.
C
C***SEE ALSO  DGMRES
C***ROUTINES CALLED  (NONE)
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
C***END PROLOGUE  DHEQR
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      INTEGER IJOB, INFO, LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), Q(*)
C     .. Local Scalars ..
      DOUBLE PRECISION C, S, T, T1, T2
      INTEGER I, IQ, J, K, KM1, KP1, NM1
C     .. Intrinsic Functions ..
      INTRINSIC ABS, SQRT
C***FIRST EXECUTABLE STATEMENT  DHEQR
      IF (IJOB .GT. 1) GO TO 70
C   -------------------------------------------------------------------
C         A new factorization is desired.
C   -------------------------------------------------------------------
C         QR decomposition without pivoting.
C
      INFO = 0
      DO 60 K = 1, N
         KM1 = K - 1
         KP1 = K + 1
C
C           Compute K-th column of R.
C           First, multiply the K-th column of A by the previous
C           K-1 Givens rotations.
C
         IF (KM1 .LT. 1) GO TO 20
         DO 10 J = 1, KM1
            I = 2*(J-1) + 1
            T1 = A(J,K)
            T2 = A(J+1,K)
            C = Q(I)
            S = Q(I+1)
            A(J,K) = C*T1 - S*T2
            A(J+1,K) = S*T1 + C*T2
 10      CONTINUE
C
C         Compute Givens components C and S.
C
 20      CONTINUE
         IQ = 2*KM1 + 1
         T1 = A(K,K)
         T2 = A(KP1,K)
         IF( T2.EQ.0.0D0 ) THEN
            C = 1
            S = 0
         ELSEIF( ABS(T2).GE.ABS(T1) ) THEN
            T = T1/T2
            S = -1.0D0/SQRT(1.0D0+T*T)
            C = -S*T
         ELSE
            T = T2/T1
            C = 1.0D0/SQRT(1.0D0+T*T)
            S = -C*T
         ENDIF
         Q(IQ) = C
         Q(IQ+1) = S
         A(K,K) = C*T1 - S*T2
         IF( A(K,K).EQ.0.0D0 ) INFO = K
 60   CONTINUE
      RETURN
C   -------------------------------------------------------------------
C         The old factorization of a will be updated.  A row and a
C         column has been added to the matrix A.  N by N-1 is now
C         the old size of the matrix.
C   -------------------------------------------------------------------
 70   CONTINUE
      NM1 = N - 1
C   -------------------------------------------------------------------
C         Multiply the new column by the N previous Givens rotations.
C   -------------------------------------------------------------------
      DO 100 K = 1,NM1
         I = 2*(K-1) + 1
         T1 = A(K,N)
         T2 = A(K+1,N)
         C = Q(I)
         S = Q(I+1)
         A(K,N) = C*T1 - S*T2
         A(K+1,N) = S*T1 + C*T2
 100  CONTINUE
C   -------------------------------------------------------------------
C         Complete update of decomposition by forming last Givens
C         rotation, and multiplying it times the column
C         vector(A(N,N),A(NP1,N)).
C   -------------------------------------------------------------------
      INFO = 0
      T1 = A(N,N)
      T2 = A(N+1,N)
      IF ( T2.EQ.0.0D0 ) THEN
         C = 1
         S = 0
      ELSEIF( ABS(T2).GE.ABS(T1) ) THEN
         T = T1/T2
         S = -1.0D0/SQRT(1.0D0+T*T)
         C = -S*T
      ELSE
         T = T2/T1
         C = 1.0D0/SQRT(1.0D0+T*T)
         S = -C*T
      ENDIF
      IQ = 2*N - 1
      Q(IQ) = C
      Q(IQ+1) = S
      A(N,N) = C*T1 - S*T2
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
C------------- LAST LINE OF DHEQR FOLLOWS ----------------------------
      END
