*DECK CH
      SUBROUTINE CH (NM, N, AR, AI, W, MATZ, ZR, ZI, FV1, FV2, FM1,
     +   IERR)
C***BEGIN PROLOGUE  CH
C***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
C            of a complex Hermitian matrix.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A3
C***TYPE      COMPLEX (RS-S, CH-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     to find the eigenvalues and eigenvectors (if desired)
C     of a COMPLEX HERMITIAN matrix.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, AR, AI, ZR and ZI, as declared in the
C          calling program dimension statement.  NM is an INTEGER
C          variable.
C
C        N is the order of the matrix A=(AR,AI).  N is an INTEGER
C          variable.  N must be less than or equal to NM.
C
C        AR and AI contain the real and imaginary parts, respectively,
C          of the complex Hermitian matrix.  AR and AI are
C          two-dimensional REAL arrays, dimensioned AR(NM,N)
C          and AI(NM,N).
C
C        MATZ is an INTEGER variable set equal to zero if only
C          eigenvalues are desired.  Otherwise, it is set to any
C          non-zero integer for both eigenvalues and eigenvectors.
C
C     On OUTPUT
C
C        W contains the eigenvalues in ascending order.
C          W is a one-dimensional REAL array, dimensioned W(N).
C
C        ZR and ZI contain the real and imaginary parts, respectively,
C          of the eigenvectors if MATZ is not zero.  ZR and ZI are
C          two-dimensional REAL arrays, dimensioned ZR(NM,N) and
C          ZI(NM,N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          10*N       if N is greater than NM,
C          J          if the J-th eigenvalue has not been
C                     determined after a total of 30 iterations.
C                     The eigenvalues should be correct for indices
C                     1, 2, ..., IERR-1, but no eigenvectors are
C                     computed.
C
C        FV1 and FV2 are one-dimensional REAL arrays used for
C          temporary storage, dimensioned FV1(N) and FV2(N).
C
C        FM1 is a two-dimensional REAL array used for temporary
C          storage, dimensioned FM1(2,N).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  HTRIBK, HTRIDI, TQL2, TQLRAT
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CH
C
      INTEGER I,J,N,NM,IERR,MATZ
      REAL AR(NM,*),AI(NM,*),W(*),ZR(NM,*),ZI(NM,*)
      REAL FV1(*),FV2(*),FM1(2,*)
C
C***FIRST EXECUTABLE STATEMENT  CH
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 CALL  HTRIDI(NM,N,AR,AI,W,FV1,FV2,FM1)
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 DO 40 I = 1, N
C
         DO 30 J = 1, N
            ZR(J,I) = 0.0E0
   30    CONTINUE
C
         ZR(I,I) = 1.0E0
   40 CONTINUE
C
      CALL  TQL2(NM,N,W,FV1,ZR,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  HTRIBK(NM,N,AR,AI,FM1,N,ZR,ZI)
   50 RETURN
      END
