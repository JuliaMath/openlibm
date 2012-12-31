*DECK SSPEV
      SUBROUTINE SSPEV (A, N, E, V, LDV, WORK, JOB, INFO)
C***BEGIN PROLOGUE  SSPEV
C***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
C            of a real symmetric matrix stored in packed form.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A1
C***TYPE      SINGLE PRECISION (SSPEV-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK, PACKED, SYMMETRIC
C***AUTHOR  Kahaner, D. K., (NBS)
C           Moler, C. B., (U. of New Mexico)
C           Stewart, G. W., (U. of Maryland)
C***DESCRIPTION
C
C     Abstract
C      SSPEV computes the eigenvalues and, optionally, the eigenvectors
C      of a real symmetric matrix stored in packed form.
C
C     Call Sequence Parameters-
C       (The values of parameters marked with * (star) will be  changed
C         by SSPEV.)
C
C        A*      REAL(N*(N+1)/2)
C                real symmetric packed input matrix.  Contains upper
C                triangle and diagonal of A, by column (elements
C                11, 12, 22, 13, 23, 33, ...).
C
C        N       INTEGER
C                set by the user to
C                the order of the matrix A.
C
C        E*      REAL(N)
C                on return from SSPEV, E contains the eigenvalues of A.
C                See also INFO below.
C
C        V*      REAL(LDV,N)
C                on return from SSPEV, if the user has set JOB
C                = 0        V is not referenced.
C                = nonzero  the N eigenvectors of A are stored in the
C                first N columns of V.  See also INFO below.
C
C        LDV     INTEGER
C                set by the user to
C                the leading dimension of the array V if JOB is also
C                set nonzero.  In that case, N must be .LE. LDV.
C                If JOB is set to zero, LDV is not referenced.
C
C        WORK*   REAL(2N)
C                temporary storage vector.  Contents changed by SSPEV.
C
C        JOB     INTEGER
C                set by the user to
C                = 0        eigenvalues only to be calculated by SSPEV.
C                           Neither V nor LDV are referenced.
C                = nonzero  eigenvalues and vectors to be calculated.
C                           In this case, A & V must be distinct arrays.
C                           Also, if LDA .GT. LDV, SSPEV changes all the
C                           elements of A thru column N.  If LDA < LDV,
C                           SSPEV changes all the elements of V through
C                           column N.  If LDA=LDV, only A(I,J) and V(I,
C                           J) for I,J = 1,...,N are changed by SSPEV.
C
C       INFO*   INTEGER
C               on return from SSPEV, the value of INFO is
C               = 0 for normal return.
C               = K if the eigenvalue iteration fails to converge.
C                   Eigenvalues and vectors 1 through K-1 are correct.
C
C
C     Error Messages-
C          No. 1   recoverable  N is greater than LDV and JOB is nonzero
C          No. 2   recoverable  N is less than one
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IMTQL2, TQLRAT, TRBAK3, TRED3, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800808  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  SSPEV
      INTEGER I,INFO,J,LDV,M,N
      REAL A(*),E(*),V(LDV,*),WORK(*)
C***FIRST EXECUTABLE STATEMENT  SSPEV
       IF (N .GT. LDV) CALL XERMSG ('SLATEC', 'SSPEV', 'N .GT. LDV.',
     +    1, 1)
       IF(N .GT. LDV) RETURN
       IF (N .LT. 1) CALL XERMSG ('SLATEC', 'SSPEV', 'N .LT. 1', 2, 1)
       IF(N .LT. 1) RETURN
C
C       CHECK N=1 CASE
C
      E(1) = A(1)
      INFO = 0
      IF(N .EQ. 1) RETURN
C
      IF(JOB.NE.0) GO TO 20
C
C     EIGENVALUES ONLY
C
      CALL TRED3(N,1,A,E,WORK(1),WORK(N+1))
      CALL TQLRAT(N,E,WORK(N+1),INFO)
      RETURN
C
C     EIGENVALUES AND EIGENVECTORS
C
   20 CALL TRED3(N,1,A,E,WORK(1),WORK(1))
      DO 30 I = 1, N
        DO 25 J = 1, N
   25     V(I,J) = 0.
   30   V(I,I) = 1.
      CALL IMTQL2(LDV,N,E,WORK,V,INFO)
      M = N
      IF(INFO .NE. 0) M = INFO - 1
      CALL TRBAK3(LDV,N,1,A,M,V)
      RETURN
      END
