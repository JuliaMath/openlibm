*DECK SUDS
      SUBROUTINE SUDS (A, X, B, NEQ, NUK, NRDA, IFLAG, MLSO, WORK,
     +   IWORK)
C***BEGIN PROLOGUE  SUDS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SUDS-S, DSUDS-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     SUDS solves the underdetermined system of linear equations A Z = B
C     where A is NEQ by NUK and NEQ .LE. NUK. In particular, if rank A
C     equals IRA, a vector X and a matrix U are determined such that
C     X is the UNIQUE solution of smallest length, satisfying A X = B,
C     and the columns of U form an orthonormal basis for the null
C     space of A, satisfying A U = 0 . Then all solutions Z are
C     given by
C              Z = X + C(1)*U(1) + ..... + C(NUK-IRA)*U(NUK-IRA)
C     where U(J) represents the J-th column of U and the C(J) are
C     arbitrary constants.
C     If the system of equations are not compatible, only the least
C     squares solution of minimal length is computed.
C     SUDS is an interfacing routine which calls subroutine LSSUDS
C     for the solution. LSSUDS in turn calls subroutine ORTHOR and
C     possibly subroutine OHTROL for the decomposition of A by
C     orthogonal transformations. In the process, ORTHOR calls upon
C     subroutine CSCALE for scaling.
C
C **********************************************************************
C   INPUT
C **********************************************************************
C
C     A -- Contains the matrix of NEQ equations in NUK unknowns and must
C          be dimensioned NRDA by NUK. The original A is destroyed.
C     X -- Solution array of length at least NUK
C     B -- Given constant vector of length NEQ, B is destroyed
C     NEQ -- Number of equations, NEQ greater or equal to 1
C     NUK -- Number of columns in the matrix (which is also the number
C            of unknowns), NUK not smaller than NEQ
C     NRDA -- Row dimension of A, NRDA greater or equal to NEQ
C     IFLAG -- Status indicator
C           =0  For the first call (and for each new problem defined by
C               a new matrix A) when the matrix data is treated as exact
C           =-K For the first call (and for each new problem defined by
C               a new matrix A) when the matrix data is assumed to be
C               accurate to about K digits
C           =1  For subsequent calls whenever the matrix A has already
C               been decomposed (problems with new vectors B but
C               same matrix A can be handled efficiently)
C     MLSO -- =0 If only the minimal length solution is wanted
C             =1 If the complete solution is wanted, includes the
C                linear space defined by the matrix U in the abstract
C     WORK(*),IWORK(*) -- Arrays for storage of internal information,
C                WORK must be dimensioned at least
C                       NUK + 3*NEQ + MLSO*NUK*(NUK-rank A)
C                where it is possible for   0 .LE. rank A .LE. NEQ
C                IWORK must be dimensioned at least   3 + NEQ
C     IWORK(2) -- Scaling indicator
C                 =-1 If the matrix is to be pre-scaled by
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
C            =3 If rank of matrix is less than NEQ
C               To continue simply reset IFLAG=1 and call SUDS again
C            =4 If the system of equations appears to be inconsistent.
C               However, the least squares solution of minimal length
C               was obtained.
C     X -- Minimal length least squares solution of  A X = B
C     A -- Contains the strictly upper triangular part of the reduced
C           matrix and transformation information
C     WORK(*),IWORK(*) -- Contains information needed on subsequent
C                         calls (IFLAG=1 case on input) which must not
C                         be altered.
C                         The matrix U described in the abstract is
C                         stored in the  NUK*(NUK-rank A) elements of
C                         the work array beginning at WORK(1+NUK+3*NEQ).
C                         However U is not defined when MLSO=0 or
C                         IFLAG=4.
C                         IWORK(1) Contains the numerically determined
C                         rank of the matrix A
C
C **********************************************************************
C
C***SEE ALSO  BVSUP
C***REFERENCES  H. A. Watts, Solving linear least squares problems
C                 using SODS/SUDS/CODS, Sandia Report SAND77-0683,
C                 Sandia Laboratories, 1977.
C***ROUTINES CALLED  LSSUDS
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SUDS
      DIMENSION A(NRDA,*),X(*),B(*),WORK(*),IWORK(*)
C
C***FIRST EXECUTABLE STATEMENT  SUDS
      IS=2
      IP=3
      IL=IP+NEQ
      KV=1+NEQ
      KT=KV+NEQ
      KS=KT+NEQ
      KU=KS+NUK
C
      CALL LSSUDS(A,X,B,NEQ,NUK,NRDA,WORK(KU),NUK,IFLAG,MLSO,IWORK(1),
     1            IWORK(IS),A,WORK(1),IWORK(IP),B,WORK(KV),WORK(KT),
     2            IWORK(IL),WORK(KS))
C
      RETURN
      END
