*DECK CG
      SUBROUTINE CG (NM, N, AR, AI, WR, WI, MATZ, ZR, ZI, FV1, FV2, FV3,
     +   IERR)
C***BEGIN PROLOGUE  CG
C***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
C            of a complex general matrix.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A4
C***TYPE      COMPLEX (RG-S, CG-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     to find the eigenvalues and eigenvectors (if desired)
C     of a COMPLEX GENERAL matrix.
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
C          of the complex general matrix.  AR and AI are two-dimensional
C          REAL arrays, dimensioned AR(NM,N) and AI(NM,N).
C
C        MATZ is an INTEGER variable set equal to zero if only
C          eigenvalues are desired.  Otherwise, it is set to any
C          non-zero integer for both eigenvalues and eigenvectors.
C
C     On OUTPUT
C
C        WR and WI contain the real and imaginary parts, respectively,
C          of the eigenvalues.  WR and WI are one-dimensional REAL
C          arrays, dimensioned WR(N) and WI(N).
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
C                     IERR+1, IERR+2, ..., N, but no eigenvectors are
C                     computed.
C
C        FV1, FV2, and FV3 are one-dimensional REAL arrays used for
C          temporary storage, dimensioned FV1(N), FV2(N), and FV3(N).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  CBABK2, CBAL, COMQR, COMQR2, CORTH
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CG
C
      INTEGER N,NM,IS1,IS2,IERR,MATZ
      REAL AR(NM,*),AI(NM,*),WR(*),WI(*),ZR(NM,*),ZI(NM,*)
      REAL FV1(*),FV2(*),FV3(*)
C
C***FIRST EXECUTABLE STATEMENT  CG
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 CALL  CBAL(NM,N,AR,AI,IS1,IS2,FV1)
      CALL  CORTH(NM,N,IS1,IS2,AR,AI,FV2,FV3)
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  COMQR(NM,N,IS1,IS2,AR,AI,WR,WI,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  COMQR2(NM,N,IS1,IS2,FV2,FV3,AR,AI,WR,WI,ZR,ZI,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  CBABK2(NM,N,IS1,IS2,FV1,N,ZR,ZI)
   50 RETURN
      END
