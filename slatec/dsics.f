*DECK DSICS
      SUBROUTINE DSICS (N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL, D,
     +   R, IWARN)
C***BEGIN PROLOGUE  DSICS
C***PURPOSE  Incompl. Cholesky Decomposition Preconditioner SLAP Set Up.
C            Routine to generate the Incomplete Cholesky decomposition,
C            L*D*L-trans, of a symmetric positive definite matrix, A,
C            which is stored in SLAP Column format.  The unit lower
C            triangular matrix L is stored by rows, and the inverse of
C            the diagonal matrix D is stored.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2E
C***TYPE      DOUBLE PRECISION (SSICS-S, DSICS-D)
C***KEYWORDS  INCOMPLETE CHOLESKY FACTORIZATION,
C             ITERATIVE PRECONDITION, LINEAR SYSTEM, SLAP SPARSE
C***AUTHOR  Greenbaum, Anne, (Courant Institute)
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     INTEGER NEL, IEL(NEL), JEL(NEL), IWARN
C     DOUBLE PRECISION A(NELT), EL(NEL), D(N), R(N)
C
C     CALL DSICS( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL, D, R,
C    $    IWARN )
C
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Double Precision A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C NEL    :OUT      Integer.
C         Number of non-zeros in the lower triangle of A.   Also
C         corresponds to the length of the IEL, JEL, EL arrays.
C IEL    :OUT      Integer IEL(NEL).
C JEL    :OUT      Integer JEL(NEL).
C EL     :OUT      Double Precision EL(NEL).
C         IEL, JEL, EL contain the unit lower triangular factor  of the
C         incomplete decomposition   of the A  matrix  stored  in  SLAP
C         Row format.   The Diagonal of   ones   *IS*   stored.     See
C         "Description", below for more details about the SLAP Row fmt.
C D      :OUT      Double Precision D(N)
C         Upon return this array holds D(I) = 1./DIAG(A).
C R      :WORK     Double Precision R(N).
C         Temporary double precision workspace needed for the
C         factorization.
C IWARN  :OUT      Integer.
C         This is a warning variable and is zero if the IC factoriza-
C         tion goes well.  It is set to the row index corresponding to
C         the last zero pivot found.  See "Description", below.
C
C *Description
C       =================== S L A P Column format ==================
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except for  the diagonal entry, which
C       must appear first in each  "column")  and are stored  in the
C       double precision array A.   In other words,  for each column
C       in the matrix put the diagonal entry in  A.  Then put in the
C       other non-zero  elements going down  the column (except  the
C       diagonal) in order.   The  IA array holds the  row index for
C       each non-zero.  The JA array holds the offsets  into the IA,
C       A arrays  for  the  beginning  of each   column.   That  is,
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the
C       number of columns in  the matrix and NELT  is the number  of
C       non-zeros in the matrix.
C
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       ==================== S L A P Row format ====================
C
C       This routine requires  that the matrix A  be  stored  in the
C       SLAP  Row format.   In this format  the non-zeros are stored
C       counting across  rows (except for the diagonal  entry, which
C       must  appear first  in each  "row")  and  are stored  in the
C       double precision  array A.  In other words, for each row  in
C       the matrix  put the diagonal  entry in A.   Then put in  the
C       other  non-zero elements  going across  the row  (except the
C       diagonal) in order.  The JA array holds the column index for
C       each non-zero.  The IA array holds the offsets  into the JA,
C       A  arrays  for  the   beginning  of  each  row.    That  is,
C       JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW-
C       th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
C       are  the last elements  of the  IROW-th row.   Note  that we
C       always have  IA(N+1) = NELT+1, where N is the number of rows
C       in the matrix  and  NELT is the  number of non-zeros  in the
C       matrix.
C
C       Here is an example of the SLAP Row storage format for a  5x5
C       Matrix (in the A and JA arrays '|' denotes the end of a row):
C
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       With the SLAP  format some  of  the   "inner  loops" of this
C       routine should vectorize  on  machines with hardware support
C       for vector   gather/scatter  operations.  Your compiler  may
C       require a compiler directive to  convince it that  there are
C       no  implicit  vector  dependencies.  Compiler directives for
C       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
C       supplied with the standard SLAP distribution.
C
C       The IC factorization does not always exist for SPD matrices.
C       In the event that a zero pivot is found it is set  to be 1.0
C       and the factorization proceeds.   The integer variable IWARN
C       is set to the last row where the Diagonal was fudged.  This
C       eventuality hardly ever occurs in practice.
C
C***SEE ALSO  DCG, DSICCG
C***REFERENCES  1. Gene Golub and Charles Van Loan, Matrix Computations,
C                  Johns Hopkins University Press, Baltimore, Maryland,
C                  1983.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   890404  DATE WRITTEN
C   890404  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   900805  Changed XERRWV calls to calls to XERMSG.  (RWC)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   920511  Added complete declaration section.  (WRB)
C   920929  Corrected format of reference.  (FNF)
C   930701  Updated CATEGORY section.  (FNF, WRB)
C***END PROLOGUE  DSICS
C     .. Scalar Arguments ..
      INTEGER ISYM, IWARN, N, NEL, NELT
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), D(N), EL(NEL), R(N)
      INTEGER IA(NELT), IEL(NEL), JA(NELT), JEL(NEL)
C     .. Local Scalars ..
      DOUBLE PRECISION ELTMP
      INTEGER I, IBGN, IC, ICBGN, ICEND, ICOL, IEND, IR, IRBGN, IREND,
     +        IROW, IRR, J, JBGN, JELTMP, JEND
      CHARACTER XERN1*8
C     .. External Subroutines ..
      EXTERNAL XERMSG
C***FIRST EXECUTABLE STATEMENT  DSICS
C
C         Set the lower triangle in IEL, JEL, EL
C
      IWARN = 0
C
C         All matrix elements stored in IA, JA, A.  Pick out the lower
C         triangle (making sure that the Diagonal of EL is one) and
C         store by rows.
C
      NEL = 1
      IEL(1) = 1
      JEL(1) = 1
      EL(1) = 1
      D(1) = A(1)
CVD$R NOCONCUR
      DO 30 IROW = 2, N
C         Put in the Diagonal.
         NEL = NEL + 1
         IEL(IROW) = NEL
         JEL(NEL) = IROW
         EL(NEL) = 1
         D(IROW) = A(JA(IROW))
C
C         Look in all the lower triangle columns for a matching row.
C         Since the matrix is symmetric, we can look across the
C         IROW-th row by looking down the IROW-th column (if it is
C         stored ISYM=0)...
         IF( ISYM.EQ.0 ) THEN
            ICBGN = JA(IROW)
            ICEND = JA(IROW+1)-1
         ELSE
            ICBGN = 1
            ICEND = IROW-1
         ENDIF
         DO 20 IC = ICBGN, ICEND
            IF( ISYM.EQ.0 ) THEN
               ICOL = IA(IC)
               IF( ICOL.GE.IROW ) GOTO 20
            ELSE
               ICOL = IC
            ENDIF
            JBGN = JA(ICOL)+1
            JEND = JA(ICOL+1)-1
            IF( JBGN.LE.JEND .AND. IA(JEND).GE.IROW ) THEN
CVD$ NOVECTOR
               DO 10 J = JBGN, JEND
                  IF( IA(J).EQ.IROW ) THEN
                     NEL = NEL + 1
                     JEL(NEL) = ICOL
                     EL(NEL)  = A(J)
                     GOTO 20
                  ENDIF
 10            CONTINUE
            ENDIF
 20      CONTINUE
 30   CONTINUE
      IEL(N+1) = NEL+1
C
C         Sort ROWS of lower triangle into descending order (count out
C         along rows out from Diagonal).
C
      DO 60 IROW = 2, N
         IBGN = IEL(IROW)+1
         IEND = IEL(IROW+1)-1
         IF( IBGN.LT.IEND ) THEN
            DO 50 I = IBGN, IEND-1
CVD$ NOVECTOR
               DO 40 J = I+1, IEND
                  IF( JEL(I).GT.JEL(J) ) THEN
                     JELTMP = JEL(J)
                     JEL(J) = JEL(I)
                     JEL(I) = JELTMP
                     ELTMP = EL(J)
                     EL(J) = EL(I)
                     EL(I) = ELTMP
                  ENDIF
 40            CONTINUE
 50         CONTINUE
         ENDIF
 60   CONTINUE
C
C         Perform the Incomplete Cholesky decomposition by looping
C         over the rows.
C         Scale the first column.  Use the structure of A to pick out
C         the rows with something in column 1.
C
      IRBGN = JA(1)+1
      IREND = JA(2)-1
      DO 65 IRR = IRBGN, IREND
         IR = IA(IRR)
C         Find the index into EL for EL(1,IR).
C         Hint: it's the second entry.
         I = IEL(IR)+1
         EL(I) = EL(I)/D(1)
 65   CONTINUE
C
      DO 110 IROW = 2, N
C
C         Update the IROW-th diagonal.
C
         DO 66 I = 1, IROW-1
            R(I) = 0
 66      CONTINUE
         IBGN = IEL(IROW)+1
         IEND = IEL(IROW+1)-1
         IF( IBGN.LE.IEND ) THEN
CLLL. OPTION ASSERT (NOHAZARD)
CDIR$ IVDEP
CVD$ NODEPCHK
            DO 70 I = IBGN, IEND
               R(JEL(I)) = EL(I)*D(JEL(I))
               D(IROW) = D(IROW) - EL(I)*R(JEL(I))
 70         CONTINUE
C
C         Check to see if we have a problem with the diagonal.
C
            IF( D(IROW).LE.0.0D0 ) THEN
               IF( IWARN.EQ.0 ) IWARN = IROW
               D(IROW) = 1
            ENDIF
         ENDIF
C
C         Update each EL(IROW+1:N,IROW), if there are any.
C         Use the structure of A to determine the Non-zero elements
C         of the IROW-th column of EL.
C
         IRBGN = JA(IROW)
         IREND = JA(IROW+1)-1
         DO 100 IRR = IRBGN, IREND
            IR = IA(IRR)
            IF( IR.LE.IROW ) GOTO 100
C         Find the index into EL for EL(IR,IROW)
            IBGN = IEL(IR)+1
            IEND = IEL(IR+1)-1
            IF( JEL(IBGN).GT.IROW ) GOTO 100
            DO 90 I = IBGN, IEND
               IF( JEL(I).EQ.IROW ) THEN
                  ICEND = IEND
 91               IF( JEL(ICEND).GE.IROW ) THEN
                     ICEND = ICEND - 1
                     GOTO 91
                  ENDIF
C         Sum up the EL(IR,1:IROW-1)*R(1:IROW-1) contributions.
CLLL. OPTION ASSERT (NOHAZARD)
CDIR$ IVDEP
CVD$ NODEPCHK
                  DO 80 IC = IBGN, ICEND
                     EL(I) = EL(I) - EL(IC)*R(JEL(IC))
 80               CONTINUE
                  EL(I) = EL(I)/D(IROW)
                  GOTO 100
               ENDIF
 90         CONTINUE
C
C         If we get here, we have real problems...
            WRITE (XERN1, '(I8)') IROW
            CALL XERMSG ('SLATEC', 'DSICS',
     $         'A and EL data structure mismatch in row '// XERN1, 1, 2)
 100     CONTINUE
 110  CONTINUE
C
C         Replace diagonals by their inverses.
C
CVD$ CONCUR
      DO 120 I =1, N
         D(I) = 1.0D0/D(I)
 120  CONTINUE
      RETURN
C------------- LAST LINE OF DSICS FOLLOWS ----------------------------
      END
