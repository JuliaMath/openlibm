*DECK CGEEV
      SUBROUTINE CGEEV (A, LDA, N, E, V, LDV, WORK, JOB, INFO)
C***BEGIN PROLOGUE  CGEEV
C***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
C            of a complex general matrix.
C***LIBRARY   SLATEC
C***CATEGORY  D4A4
C***TYPE      COMPLEX (SGEEV-S, CGEEV-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, GENERAL MATRIX
C***AUTHOR  Kahaner, D. K., (NBS)
C           Moler, C. B., (U. of New Mexico)
C           Stewart, G. W., (U. of Maryland)
C***DESCRIPTION
C
C     Abstract
C      CGEEV computes the eigenvalues and, optionally,
C      the eigenvectors of a general complex matrix.
C
C     Call Sequence Parameters-
C       (The values of parameters marked with * (star) will be changed
C         by CGEEV.)
C
C        A*      COMPLEX(LDA,N)
C                complex nonsymmetric input matrix.
C
C        LDA     INTEGER
C                set by the user to
C                the leading dimension of the complex array A.
C
C        N       INTEGER
C                set by the user to
C                the order of the matrices A and V, and
C                the number of elements in E.
C
C        E*      COMPLEX(N)
C                on return from CGEEV E contains the eigenvalues of A.
C                See also INFO below.
C
C        V*      COMPLEX(LDV,N)
C                on return from CGEEV if the user has set JOB
C                = 0        V is not referenced.
C                = nonzero  the N eigenvectors of A are stored in the
C                first N columns of V.  See also INFO below.
C                (If the input matrix A is nearly degenerate, V
C                 will be badly conditioned, i.e. have nearly
C                 dependent columns.)
C
C        LDV     INTEGER
C                set by the user to
C                the leading dimension of the array V if JOB is also
C                set nonzero.  In that case N must be .LE. LDV.
C                If JOB is set to zero LDV is not referenced.
C
C        WORK*   REAL(3N)
C                temporary storage vector.  Contents changed by CGEEV.
C
C        JOB     INTEGER
C                set by the user to
C                = 0        eigenvalues only to be calculated by CGEEV.
C                           neither V nor LDV are referenced.
C                = nonzero  eigenvalues and vectors to be calculated.
C                           In this case A & V must be distinct arrays.
C                           Also,  if LDA > LDV,  CGEEV changes all the
C                           elements of A thru column N.  If LDA < LDV,
C                           CGEEV changes all the elements of V through
C                           column N.  If LDA = LDV only A(I,J) and V(I,
C                           J) for I,J = 1,...,N are changed by CGEEV.
C
C        INFO*   INTEGER
C                on return from CGEEV the value of INFO is
C                = 0  normal return, calculation successful.
C                = K  if the eigenvalue iteration fails to converge,
C                     eigenvalues K+1 through N are correct, but
C                     no eigenvectors were computed even if they were
C                     requested (JOB nonzero).
C
C      Error Messages
C           No. 1  recoverable  N is greater than LDA
C           No. 2  recoverable  N is less than one.
C           No. 3  recoverable  JOB is nonzero and N is greater than LDV
C           No. 4  warning      LDA > LDV,  elements of A other than the
C                               N by N input elements have been changed
C           No. 5  warning      LDA < LDV,  elements of V other than the
C                               N by N output elements have been changed
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CBABK2, CBAL, COMQR, COMQR2, CORTH, SCOPY, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800808  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  CGEEV
      INTEGER I,IHI,ILO,INFO,J,K,L,LDA,LDV,MDIM,N
      REAL A(*),E(*),WORK(*),V(*)
C***FIRST EXECUTABLE STATEMENT  CGEEV
      IF (N .GT. LDA) CALL XERMSG ('SLATEC', 'CGEEV', 'N .GT. LDA.', 1,
     +   1)
      IF(N .GT. LDA) RETURN
      IF (N .LT. 1) CALL XERMSG ('SLATEC', 'CGEEV', 'N .LT. 1', 2, 1)
      IF(N .LT. 1) RETURN
      IF(N .EQ. 1 .AND. JOB .EQ. 0) GO TO 35
      MDIM = 2 * LDA
      IF(JOB .EQ. 0) GO TO 5
      IF (N .GT. LDV) CALL XERMSG ('SLATEC', 'CGEEV',
     +   'JOB .NE. 0 AND N .GT. LDV.', 3, 1)
      IF(N .GT. LDV) RETURN
      IF(N .EQ. 1) GO TO 35
C
C       REARRANGE A IF NECESSARY WHEN LDA.GT.LDV AND JOB .NE.0
C
      MDIM = MIN(MDIM,2 * LDV)
      IF (LDA .LT. LDV) CALL XERMSG ('SLATEC', 'CGEEV',
     +   'LDA.LT.LDV,  ELEMENTS OF V OTHER THAN THE N BY N OUTPUT ' //
     +   'ELEMENTS HAVE BEEN CHANGED.', 5, 0)
      IF(LDA.LE.LDV) GO TO 5
      CALL XERMSG ('SLATEC', 'CGEEV',
     +   'LDA.GT.LDV, ELEMENTS OF A OTHER THAN THE N BY N INPUT ' //
     +   'ELEMENTS HAVE BEEN CHANGED.', 4, 0)
      L = N - 1
      DO 4 J=1,L
          I = 2 * N
         M = 1+J*2*LDV
         K = 1+J*2*LDA
         CALL SCOPY(I,A(K),1,A(M),1)
    4 CONTINUE
    5 CONTINUE
C
C     SEPARATE REAL AND IMAGINARY PARTS
C
      DO 6 J = 1, N
       K = (J-1) * MDIM +1
       L = K + N
       CALL SCOPY(N,A(K+1),2,WORK(1),1)
       CALL SCOPY(N,A(K),2,A(K),1)
       CALL SCOPY(N,WORK(1),1,A(L),1)
    6 CONTINUE
C
C     SCALE AND ORTHOGONAL REDUCTION TO HESSENBERG.
C
      CALL CBAL(MDIM,N,A(1),A(N+1),ILO,IHI,WORK(1))
      CALL CORTH(MDIM,N,ILO,IHI,A(1),A(N+1),WORK(N+1),WORK(2*N+1))
      IF(JOB .NE. 0) GO TO 10
C
C     EIGENVALUES ONLY
C
      CALL COMQR(MDIM,N,ILO,IHI,A(1),A(N+1),E(1),E(N+1),INFO)
      GO TO 30
C
C     EIGENVALUES AND EIGENVECTORS.
C
   10 CALL COMQR2(MDIM,N,ILO,IHI,WORK(N+1),WORK(2*N+1),A(1),A(N+1),
     1  E(1),E(N+1),V(1),V(N+1),INFO)
      IF (INFO .NE. 0) GO TO 30
      CALL CBABK2(MDIM,N,ILO,IHI,WORK(1),N,V(1),V(N+1))
C
C     CONVERT EIGENVECTORS TO COMPLEX STORAGE.
C
      DO 20 J = 1,N
       K = (J-1) * MDIM + 1
       I = (J-1) * 2 * LDV + 1
       L = K + N
       CALL SCOPY(N,V(K),1,WORK(1),1)
       CALL SCOPY(N,V(L),1,V(I+1),2)
       CALL SCOPY(N,WORK(1),1,V(I),2)
   20 CONTINUE
C
C     CONVERT EIGENVALUES TO COMPLEX STORAGE.
C
   30 CALL SCOPY(N,E(1),1,WORK(1),1)
      CALL SCOPY(N,E(N+1),1,E(2),2)
      CALL SCOPY(N,WORK(1),1,E(1),2)
      RETURN
C
C     TAKE CARE OF N=1 CASE
C
   35 E(1) = A(1)
      E(2) = A(2)
      INFO = 0
      IF(JOB .EQ. 0) RETURN
      V(1) = A(1)
      V(2) = A(2)
      RETURN
      END
