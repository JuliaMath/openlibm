*DECK RG
      SUBROUTINE RG (NM, N, A, WR, WI, MATZ, Z, IV1, FV1, IERR)
C***BEGIN PROLOGUE  RG
C***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
C            of a real general matrix.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A2
C***TYPE      SINGLE PRECISION (RG-S, CG-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     To find the eigenvalues and eigenvectors (if desired)
C     of a REAL GENERAL matrix.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix A.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        A contains the real general matrix.  A is a two-dimensional
C          REAL array, dimensioned A(NM,N).
C
C        MATZ is an INTEGER variable set equal to zero if only
C          eigenvalues are desired.  Otherwise, it is set to any
C          non-zero integer for both eigenvalues and eigenvectors.
C
C     On Output
C
C        A has been destroyed.
C
C        WR and WI contain the real and imaginary parts, respectively,
C          of the eigenvalues.  The eigenvalues are unordered except
C          that complex conjugate pairs of eigenvalues appear consecu-
C          tively with the eigenvalue having the positive imaginary part
C          first.  If an error exit is made, the eigenvalues should be
C          correct for indices IERR+1, IERR+2, ..., N.  WR and WI are
C          one-dimensional REAL arrays, dimensioned WR(N) and WI(N).
C
C        Z contains the real and imaginary parts of the eigenvectors
C          if MATZ is not zero.  If the J-th eigenvalue is real, the
C          J-th column of Z contains its eigenvector.  If the J-th
C          eigenvalue is complex with positive imaginary part, the
C          J-th and (J+1)-th columns of Z contain the real and
C          imaginary parts of its eigenvector.  The conjugate of this
C          vector is the eigenvector for the conjugate eigenvalue.
C          Z is a two-dimensional REAL array, dimensioned Z(NM,N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          10*N       if N is greater than NM,
C          J          if the J-th eigenvalue has not been
C                     determined after a total of 30 iterations.
C                     The eigenvalues should be correct for indices
C                     IERR+1, IERR+2, ..., N, but no eigenvectors are
C                     computed.
C
C        IV1 and FV1 are one-dimensional temporary storage arrays of
C          dimension N.  IV1 is of type INTEGER and FV1 of type REAL.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  BALANC, BALBAK, ELMHES, ELTRAN, HQR, HQR2
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   921103  Corrected description of IV1.  (DWL, FNF and WRB)
C***END PROLOGUE  RG
C
      INTEGER N,NM,IS1,IS2,IERR,MATZ
      REAL A(NM,*),WR(*),WI(*),Z(NM,*),FV1(*)
      INTEGER IV1(*)
C
C***FIRST EXECUTABLE STATEMENT  RG
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 CALL  BALANC(NM,N,A,IS1,IS2,FV1)
      CALL  ELMHES(NM,N,IS1,IS2,A,IV1)
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  HQR(NM,N,IS1,IS2,A,WR,WI,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  ELTRAN(NM,N,IS1,IS2,A,IV1,Z)
      CALL  HQR2(NM,N,IS1,IS2,A,WR,WI,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  BALBAK(NM,N,IS1,IS2,FV1,N,Z)
   50 RETURN
      END
