*DECK SORTH
      SUBROUTINE SORTH (VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
C***BEGIN PROLOGUE  SORTH
C***SUBSIDIARY
C***PURPOSE  Internal routine for SGMRES.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      SINGLE PRECISION (SORTH-S, DORTH-D)
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
C           Seager, Mark K., (LLNL), seager@llnl.gov
C             Lawrence Livermore National Laboratory
C             PO Box 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C***DESCRIPTION
C        This routine  orthogonalizes  the  vector  VNEW  against the
C        previous KMP  vectors in the   V array.  It uses  a modified
C        Gram-Schmidt   orthogonalization procedure with  conditional
C        reorthogonalization.
C
C *Usage:
C      INTEGER N, LL, LDHES, KMP
C      REAL VNEW(N), V(N,LL), HES(LDHES,LL), SNORMW
C
C      CALL SORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
C
C *Arguments:
C VNEW   :INOUT    Real VNEW(N)
C         On input, the vector of length N containing a scaled
C         product of the Jacobian and the vector V(*,LL).
C         On output, the new vector orthogonal to V(*,i0) to V(*,LL),
C         where i0 = max(1, LL-KMP+1).
C V      :IN       Real V(N,LL)
C         The N x LL array containing the previous LL
C         orthogonal vectors V(*,1) to V(*,LL).
C HES    :INOUT    Real HES(LDHES,LL)
C         On input, an LL x LL upper Hessenberg matrix containing,
C         in HES(I,K), K.lt.LL, the scaled inner products of
C         A*V(*,K) and V(*,i).
C         On return, column LL of HES is filled in with
C         the scaled inner products of A*V(*,LL) and V(*,i).
C N      :IN       Integer
C         The order of the matrix A, and the length of VNEW.
C LL     :IN       Integer
C         The current order of the matrix HES.
C LDHES  :IN       Integer
C         The leading dimension of the HES array.
C KMP    :IN       Integer
C         The number of previous vectors the new vector VNEW
C         must be made orthogonal to (KMP .le. MAXL).
C SNORMW :OUT      REAL
C         Scalar containing the l-2 norm of VNEW.
C
C***SEE ALSO  SGMRES
C***ROUTINES CALLED  SAXPY, SDOT, SNRM2
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
C***END PROLOGUE  SORTH
C         The following is for optimized compilation on LLNL/LTSS Crays.
CLLL. OPTIMIZE
C     .. Scalar Arguments ..
      REAL SNORMW
      INTEGER KMP, LDHES, LL, N
C     .. Array Arguments ..
      REAL HES(LDHES,*), V(N,*), VNEW(*)
C     .. Local Scalars ..
      REAL ARG, SUMDSQ, TEM, VNRM
      INTEGER I, I0
C     .. External Functions ..
      REAL SDOT, SNRM2
      EXTERNAL SDOT, SNRM2
C     .. External Subroutines ..
      EXTERNAL SAXPY
C     .. Intrinsic Functions ..
      INTRINSIC MAX, SQRT
C***FIRST EXECUTABLE STATEMENT  SORTH
C
C         Get norm of unaltered VNEW for later use.
C
      VNRM = SNRM2(N, VNEW, 1)
C   -------------------------------------------------------------------
C         Perform the modified Gram-Schmidt procedure on VNEW =A*V(LL).
C         Scaled inner products give new column of HES.
C         Projections of earlier vectors are subtracted from VNEW.
C   -------------------------------------------------------------------
      I0 = MAX(1,LL-KMP+1)
      DO 10 I = I0,LL
         HES(I,LL) = SDOT(N, V(1,I), 1, VNEW, 1)
         TEM = -HES(I,LL)
         CALL SAXPY(N, TEM, V(1,I), 1, VNEW, 1)
 10   CONTINUE
C   -------------------------------------------------------------------
C         Compute SNORMW = norm of VNEW.  If VNEW is small compared
C         to its input value (in norm), then reorthogonalize VNEW to
C         V(*,1) through V(*,LL).  Correct if relative correction
C         exceeds 1000*(unit roundoff).  Finally, correct SNORMW using
C         the dot products involved.
C   -------------------------------------------------------------------
      SNORMW = SNRM2(N, VNEW, 1)
      IF (VNRM + 0.001E0*SNORMW .NE. VNRM) RETURN
      SUMDSQ = 0
      DO 30 I = I0,LL
         TEM = -SDOT(N, V(1,I), 1, VNEW, 1)
         IF (HES(I,LL) + 0.001E0*TEM .EQ. HES(I,LL)) GO TO 30
         HES(I,LL) = HES(I,LL) - TEM
         CALL SAXPY(N, TEM, V(1,I), 1, VNEW, 1)
         SUMDSQ = SUMDSQ + TEM**2
 30   CONTINUE
      IF (SUMDSQ .EQ. 0.0E0) RETURN
      ARG = MAX(0.0E0,SNORMW**2 - SUMDSQ)
      SNORMW = SQRT(ARG)
C
      RETURN
C------------- LAST LINE OF SORTH FOLLOWS ----------------------------
      END
