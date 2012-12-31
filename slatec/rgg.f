*DECK RGG
      SUBROUTINE RGG (NM, N, A, B, ALFR, ALFI, BETA, MATZ, Z, IERR)
C***BEGIN PROLOGUE  RGG
C***PURPOSE  Compute the eigenvalues and eigenvectors for a real
C            generalized eigenproblem.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4B2
C***TYPE      SINGLE PRECISION (RGG-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     to find the eigenvalues and eigenvectors (if desired)
C     for the REAL GENERAL GENERALIZED eigenproblem  Ax = (LAMBDA)Bx.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A, B, and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrices A and B.  N is an INTEGER
C          variable.  N must be less than or equal to NM.
C
C        A contains a real general matrix.  A is a two-dimensional
C          REAL array, dimensioned A(NM,N).
C
C        B contains a real general matrix.  B is a two-dimensional
C          REAL array, dimensioned B(NM,N).
C
C        MATZ is an INTEGER variable set equal to zero if only
C          eigenvalues are desired.  Otherwise, it is set to any
C          non-zero integer for both eigenvalues and eigenvectors.
C
C     On Output
C
C        A and B have been destroyed.
C
C        ALFR and ALFI contain the real and imaginary parts,
C          respectively, of the numerators of the eigenvalues.
C          ALFR and ALFI are one-dimensional REAL arrays,
C          dimensioned ALFR(N) and ALFI(N).
C
C        BETA contains the denominators of the eigenvalues,
C          which are thus given by the ratios  (ALFR+I*ALFI)/BETA.
C          Complex conjugate pairs of eigenvalues appear consecutively
C          with the eigenvalue having the positive imaginary part first.
C          BETA is a one-dimensional REAL array, dimensioned BETA(N).
C
C        Z contains the real and imaginary parts of the eigenvectors
C          if MATZ is not zero.  If the J-th eigenvalue is real, the
C          J-th column of  Z  contains its eigenvector.  If the J-th
C          eigenvalue is complex with positive imaginary part, the
C          J-th and (J+1)-th columns of  Z  contain the real and
C          imaginary parts of its eigenvector.  The conjugate of this
C          vector is the eigenvector for the conjugate eigenvalue.
C          Z is a two-dimensional REAL array, dimensioned Z(NM,N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          10*N       if N is greater than NM,
C          J          if the J-th eigenvalue has not been
C                     determined after a total of 30*N iterations.
C                     The eigenvalues should be correct for indices
C                     IERR+1, IERR+2, ..., N, but no eigenvectors are
C                     computed.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  QZHES, QZIT, QZVAL, QZVEC
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RGG
C
      INTEGER N,NM,IERR,MATZ
      REAL A(NM,*),B(NM,*),ALFR(*),ALFI(*),BETA(*),Z(NM,*)
      LOGICAL TF
C
C***FIRST EXECUTABLE STATEMENT  RGG
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      TF = .FALSE.
      CALL  QZHES(NM,N,A,B,TF,Z)
      CALL  QZIT(NM,N,A,B,0.0E0,TF,Z,IERR)
      CALL  QZVAL(NM,N,A,B,ALFR,ALFI,BETA,TF,Z)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 TF = .TRUE.
      CALL  QZHES(NM,N,A,B,TF,Z)
      CALL  QZIT(NM,N,A,B,0.0E0,TF,Z,IERR)
      CALL  QZVAL(NM,N,A,B,ALFR,ALFI,BETA,TF,Z)
      IF (IERR .NE. 0) GO TO 50
      CALL  QZVEC(NM,N,A,B,ALFR,ALFI,BETA,Z)
   50 RETURN
      END
