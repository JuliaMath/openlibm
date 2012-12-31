*DECK SODS
      SUBROUTINE SODS (A, X, B, NEQ, NUK, NRDA, IFLAG, WORK, IWORK)
C***BEGIN PROLOGUE  SODS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SODS-S)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     SODS solves the overdetermined system of linear equations A X = B,
C     where A is NEQ by NUK and NEQ .GE. NUK. If rank A = NUK,
C     X is the UNIQUE least squares solution vector. That is,
C              R(1)**2 + ..... + R(NEQ)**2 = minimum
C     where R is the residual vector  R = B - A X.
C     If rank A .LT. NUK , the least squares solution of minimal
C     length can be provided.
C     SODS is an interfacing routine which calls subroutine LSSODS
C     for the solution. LSSODS in turn calls subroutine ORTHOL and
C     possibly subroutine OHTROR for the decomposition of A by
C     orthogonal transformations. In the process, ORTHOL calls upon
C     subroutine CSCALE for scaling.
C
C **********************************************************************
C   Input
C **********************************************************************
C
C     A -- Contains the matrix of NEQ equations in NUK unknowns and must
C          be dimensioned NRDA by NUK. The original A is destroyed
C     X -- Solution array of length at least NUK
C     B -- Given constant vector of length NEQ, B is destroyed
C     NEQ -- Number of equations, NEQ greater or equal to 1
C     NUK -- Number of columns in the matrix (which is also the number
C            of unknowns), NUK not larger than NEQ
C     NRDA -- Row dimension of A, NRDA greater or equal to NEQ
C     IFLAG -- Status indicator
C            =0 For the first call (and for each new problem defined by
C               a new matrix A) when the matrix data is treated as exact
C           =-K For the first call (and for each new problem defined by
C               a new matrix A) when the matrix data is assumed to be
C               accurate to about K digits
C            =1 For subsequent calls whenever the matrix A has already
C               been decomposed (problems with new vectors B but
C               same matrix a can be handled efficiently)
C     WORK(*),IWORK(*) -- Arrays for storage of internal information,
C                     WORK must be dimensioned at least  2 + 5*NUK
C                     IWORK must be dimensioned at least NUK+2
C     IWORK(2) -- Scaling indicator
C                 =-1 If the matrix A is to be pre-scaled by
C                 columns when appropriate
C                 If the scaling indicator is not equal to -1
C                 no scaling will be attempted
C              For most problems scaling will probably not be necessary
C
C **********************************************************************
C   OUTPUT
C **********************************************************************
C
C     IFLAG -- Status indicator
C            =1 If solution was obtained
C            =2 If improper input is detected
C            =3 If rank of matrix is less than NUK
C               If the minimal length least squares solution is
C               desired, simply reset IFLAG=1 and call the code again
C     X -- Least squares solution of  A X = B
C     A -- Contains the strictly upper triangular part of the reduced
C           matrix and the transformation information
C     WORK(*),IWORK(*) -- Contains information needed on subsequent
C                         Calls (IFLAG=1 case on input) which must not
C                         be altered
C                         WORK(1) contains the Euclidean norm of
C                         the residual vector
C                         WORK(2) contains the Euclidean norm of
C                         the solution vector
C                         IWORK(1) contains the numerically determined
C                         rank of the matrix A
C
C **********************************************************************
C
C***SEE ALSO  BVSUP
C***REFERENCES  G. Golub, Numerical methods for solving linear least
C                 squares problems, Numerische Mathematik 7, (1965),
C                 pp. 206-216.
C               P. Businger and G. Golub, Linear least squares
C                 solutions by Householder transformations, Numerische
C                 Mathematik  7, (1965), pp. 269-276.
C               H. A. Watts, Solving linear least squares problems
C                 using SODS/SUDS/CODS, Sandia Report SAND77-0683,
C                 Sandia Laboratories, 1977.
C***ROUTINES CALLED  LSSODS
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SODS
      DIMENSION A(NRDA,*),X(*),B(*),WORK(*),IWORK(*)
C
C***FIRST EXECUTABLE STATEMENT  SODS
      ITER=0
      IS=2
      IP=3
      KS=2
      KD=3
      KZ=KD+NUK
      KV=KZ+NUK
      KT=KV+NUK
      KC=KT+NUK
C
      CALL LSSODS(A,X,B,NEQ,NUK,NRDA,IFLAG,IWORK(1),IWORK(IS),A,
     1            WORK(KD),IWORK(IP),ITER,WORK(1),WORK(KS),
     2            WORK(KZ),B,WORK(KV),WORK(KT),WORK(KC))
C
      RETURN
      END
