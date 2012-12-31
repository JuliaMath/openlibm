*DECK EISDOC
      SUBROUTINE EISDOC
C***BEGIN PROLOGUE  EISDOC
C***PURPOSE  Documentation for EISPACK, a collection of subprograms for
C            solving matrix eigen-problems.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4, Z
C***TYPE      ALL (EISDOC-A)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Vandevender, W. H., (SNLA)
C***DESCRIPTION
C
C                 **********EISPACK Routines**********
C
C single double complx
C ------ ------ ------
C
C RS       -    CH     Computes eigenvalues and, optionally,
C                      eigenvectors of real symmetric
C                      (complex Hermitian) matrix.
C
C RSP      -      -    Compute eigenvalues and, optionally,
C                      eigenvectors of real symmetric matrix
C                      packed into a one dimensional array.
C
C RG       -    CG     Computes eigenvalues and, optionally,
C                      eigenvectors of a real (complex) general
C                      matrix.
C
C BISECT   -      -    Compute eigenvalues of symmetric tridiagonal
C                      matrix given interval using Sturm sequencing.
C
C IMTQL1   -      -    Computes eigenvalues of symmetric tridiagonal
C                      matrix implicit QL method.
C
C IMTQL2   -      -    Computes eigenvalues and eigenvectors of
C                      symmetric tridiagonal matrix using
C                      implicit QL method.
C
C IMTQLV   -      -    Computes eigenvalues of symmetric tridiagonal
C                      matrix by the implicit QL method.
C                      Eigenvectors may be computed later.
C
C RATQR    -      -    Computes largest or smallest eigenvalues
C                      of symmetric tridiagonal matrix using
C                      rational QR method with Newton correction.
C
C RST      -      -    Compute eigenvalues and, optionally,
C                      eigenvectors of real symmetric tridiagonal
C                      matrix.
C
C RT       -      -    Compute eigenvalues and eigenvectors of
C                      a special real tridiagonal matrix.
C
C TQL1     -      -    Compute eigenvalues of symmetric tridiagonal
C                      matrix by QL method.
C
C TQL2     -      -    Compute eigenvalues and eigenvectors
C                      of symmetric tridiagonal matrix.
C
C TQLRAT   -      -    Computes eigenvalues of symmetric
C                      tridiagonal matrix a rational variant
C                      of the QL method.
C
C TRIDIB   -      -    Computes eigenvalues of symmetric
C                      tridiagonal matrix given interval using
C                      Sturm sequencing.
C
C TSTURM   -      -    Computes eigenvalues of symmetric tridiagonal
C                      matrix given interval and eigenvectors
C                      by Sturm sequencing.  This subroutine
C                      is a translation of the ALGOL procedure
C                      TRISTURM by Peters and Wilkinson. HANDBOOK
C                      FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA,
C                      418-439(1971).
C
C BQR      -      -    Computes some of the eigenvalues of a real
C                      symmetric matrix using the QR method with
C                      shifts of origin.
C
C RSB      -      -    Computes eigenvalues and, optionally,
C                      eigenvectors of symmetric band matrix.
C
C RSG      -      -    Computes eigenvalues and, optionally,
C                      eigenvectors of symmetric generalized
C                      eigenproblem: A*X=(LAMBDA)*B*X
C
C RSGAB    -      -    Computes eigenvalues and, optionally,
C                      eigenvectors of symmetric generalized
C                      eigenproblem: A*B*X=(LAMBDA)*X
C
C RSGBA    -      -    Computes eigenvalues and, optionally,
C                      eigenvectors of symmetric generalized
C                      eigenproblem: B*A*X=(LAMBDA)*X
C
C RGG      -      -    Computes eigenvalues and eigenvectors
C                      for real generalized eigenproblem:
C                      A*X=(LAMBDA)*B*X.
C
C BALANC   -    CBAL   Balances a general real (complex)
C                      matrix and isolates eigenvalues whenever
C                      possible.
C
C BANDR    -      -    Reduces real symmetric band matrix
C                      to symmetric tridiagonal matrix and,
C                      optionally, accumulates orthogonal similarity
C                      transformations.
C
C HTRID3   -      -    Reduces complex Hermitian (packed) matrix
C                      to real symmetric tridiagonal matrix by unitary
C                      similarity transformations.
C
C HTRIDI   -      -    Reduces complex Hermitian matrix to real
C                      symmetric tridiagonal matrix using unitary
C                      similarity transformations.
C
C TRED1    -      -    Reduce real symmetric matrix to symmetric
C                      tridiagonal matrix using orthogonal
C                      similarity transformations.
C
C TRED2    -      -    Reduce real symmetric matrix to symmetric
C                      tridiagonal matrix using and accumulating
C                      orthogonal transformations.
C
C TRED3    -      -    Reduce  symmetric matrix stored in packed
C                      form to symmetric tridiagonal matrix using
C                      orthogonal transformations.
C
C ELMHES   -    COMHES Reduces real (complex) general matrix to
C                      upper Hessenberg form using stabilized
C                      elementary similarity transformations.
C
C ORTHES   -    CORTH  Reduces real (complex) general matrix to upper
C                      Hessenberg form orthogonal (unitary)
C                      similarity transformations.
C
C QZHES    -      -    The first step of the QZ algorithm for solving
C                      generalized matrix eigenproblems.  Accepts
C                      a pair of real general matrices and reduces
C                      one of them to upper Hessenberg and the other
C                      to upper triangular form using orthogonal
C                      transformations. Usually followed by QZIT,
C                      QZVAL, QZ
C
C QZIT     -      -    The second step of the QZ algorithm for
C                      generalized eigenproblems.  Accepts an upper
C                      Hessenberg and an upper triangular matrix
C                      and reduces the former to quasi-triangular
C                      form while preserving the form of the latter.
C                      Usually preceded by QZHES and followed by QZVAL
C                      and QZVEC.
C
C FIGI     -      -    Transforms certain real non-symmetric
C                      tridiagonal matrix to symmetric tridiagonal
C                      matrix.
C
C FIGI2    -      -    Transforms certain real non-symmetric
C                      tridiagonal matrix to symmetric tridiagonal
C                      matrix.
C
C REDUC    -      -    Reduces generalized symmetric eigenproblem
C                      A*X=(LAMBDA)*B*X, to standard symmetric
C                      eigenproblem using Cholesky factorization.
C
C REDUC2   -      -    Reduces certain generalized symmetric
C                      eigenproblems standard symmetric eigenproblem,
C                      using Cholesky factorization.
C
C   -      -    COMLR  Computes eigenvalues of a complex upper
C                      Hessenberg matrix using the modified LR method.
C
C   -      -    COMLR2 Computes eigenvalues and eigenvectors of
C                      complex upper Hessenberg matrix using
C                      modified LR method.
C
C HQR      -    COMQR  Computes eigenvalues of a real (complex)
C                      upper Hessenberg matrix using the QR method.
C
C HQR2     -    COMQR2 Computes eigenvalues and eigenvectors of
C                      real (complex) upper Hessenberg matrix
C                      using QR method.
C
C INVIT    -    CINVIT Computes eigenvectors of real (complex)
C                      Hessenberg matrix associated with specified
C                      eigenvalues by inverse iteration.
C
C QZVAL    -      -    The third step of the QZ algorithm for
C                      generalized eigenproblems.  Accepts a pair
C                      of real matrices, one quasi-triangular form
C                      and the other in upper triangular form and
C                      computes the eigenvalues of the associated
C                      eigenproblem.  Usually preceded by QZHES,
C                      QZIT, and followed by QZVEC.
C
C BANDV    -      -    Forms eigenvectors of real symmetric band
C                      matrix associated with a set of ordered
C                      approximate eigenvalue by inverse iteration.
C
C QZVEC    -      -    The optional fourth step of the QZ algorithm
C                      for generalized eigenproblems.  Accepts
C                      a matrix in quasi-triangular form and another
C                      in upper triangular and computes the
C                      eigenvectors of the triangular problem
C                      and transforms them back to the original
C                      coordinates Usually preceded by QZHES, QZIT,
C                      QZVAL.
C
C TINVIT   -      -    Eigenvectors of symmetric tridiagonal
C                      matrix corresponding to some specified
C                      eigenvalues, using inverse iteration.
C
C BAKVEC   -      -    Forms eigenvectors of certain real
C                      non-symmetric tridiagonal matrix from
C                      symmetric tridiagonal matrix output from FIGI.
C
C BALBAK   -    CBABK2 Forms eigenvectors of real (complex) general
C                      matrix from eigenvectors of matrix output
C                      from BALANC (CBAL).
C
C ELMBAK   -    COMBAK Forms eigenvectors of real (complex) general
C                      matrix from eigenvectors of upper Hessenberg
C                      matrix output from ELMHES (COMHES).
C
C ELTRAN   -      -    Accumulates the stabilized elementary
C                      similarity transformations used in the
C                      reduction of a real general matrix to upper
C                      Hessenberg form by ELMHES.
C
C HTRIB3   -      -    Computes eigenvectors of complex Hermitian
C                      matrix from eigenvectors of real symmetric
C                      tridiagonal matrix output from HTRID3.
C
C HTRIBK   -      -    Forms eigenvectors of complex Hermitian
C                      matrix from eigenvectors of real symmetric
C                      tridiagonal matrix output from HTRIDI.
C
C ORTBAK   -    CORTB  Forms eigenvectors of general real (complex)
C                      matrix from eigenvectors of upper Hessenberg
C                      matrix output from ORTHES (CORTH).
C
C ORTRAN   -      -    Accumulates orthogonal similarity
C                      transformations in reduction of real general
C                      matrix by ORTHES.
C
C REBAK    -      -    Forms eigenvectors of generalized symmetric
C                      eigensystem from eigenvectors of derived
C                      matrix output from REDUC or REDUC2.
C
C REBAKB   -      -    Forms eigenvectors of generalized symmetric
C                      eigensystem from eigenvectors of derived
C                      matrix output from REDUC2
C
C TRBAK1   -      -    Forms the eigenvectors of real symmetric
C                      matrix from eigenvectors of symmetric
C                      tridiagonal matrix formed by TRED1.
C
C TRBAK3   -      -    Forms eigenvectors of real symmetric matrix
C                      from the eigenvectors of symmetric tridiagonal
C                      matrix formed by TRED3.
C
C MINFIT   -      -    Compute Singular Value Decomposition
C                      of rectangular matrix and solve related
C                      Linear Least Squares problem.
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   811101  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900723  PURPOSE section revised.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  EISDOC
C***FIRST EXECUTABLE STATEMENT  EISDOC
      RETURN
      END
