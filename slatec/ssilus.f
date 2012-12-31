*DECK SSILUS
      SUBROUTINE SSILUS (N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, DINV,
     +   NU, IU, JU, U, NROW, NCOL)
C***BEGIN PROLOGUE  SSILUS
C***PURPOSE  Incomplete LU Decomposition Preconditioner SLAP Set Up.
C            Routine to generate the incomplete LDU decomposition of a
C            matrix.  The unit lower triangular factor L is stored by
C            rows and the unit upper triangular factor U is stored by
C            columns.  The inverse of the diagonal matrix D is stored.
C            No fill in is allowed.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2E
C***TYPE      SINGLE PRECISION (SSILUS-S, DSILUS-D)
C***KEYWORDS  INCOMPLETE LU FACTORIZATION, ITERATIVE PRECONDITION,
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
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
C     INTEGER NL, IL(NL), JL(NL), NU, IU(NU), JU(NU)
C     INTEGER NROW(N), NCOL(N)
C     REAL    A(NELT), L(NL), DINV(N), U(NU)
C
C     CALL SSILUS( N, NELT, IA, JA, A, ISYM, NL, IL, JL, L,
C    $    DINV, NU, IU, JU, U, NROW, NCOL )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Real A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C NL     :OUT      Integer.
C         Number of non-zeros in the L array.
C IL     :OUT      Integer IL(NL).
C JL     :OUT      Integer JL(NL).
C L      :OUT      Real     L(NL).
C         IL, JL, L  contain the unit lower triangular factor of  the
C         incomplete decomposition  of some  matrix stored  in   SLAP
C         Row format.     The   Diagonal  of ones  *IS*  stored.  See
C         "DESCRIPTION", below for more details about the SLAP format.
C NU     :OUT      Integer.
C         Number of non-zeros in the U array.
C IU     :OUT      Integer IU(NU).
C JU     :OUT      Integer JU(NU).
C U      :OUT      Real     U(NU).
C         IU, JU, U contain   the unit upper triangular factor of the
C         incomplete  decomposition    of some matrix  stored in SLAP
C         Column  format.   The Diagonal of ones   *IS*  stored.  See
C         "Description", below  for  more  details  about  the   SLAP
C         format.
C NROW   :WORK     Integer NROW(N).
C         NROW(I) is the number of non-zero elements in the I-th row
C         of L.
C NCOL   :WORK     Integer NCOL(N).
C         NCOL(I) is the number of non-zero elements in the I-th
C         column of U.
C
C *Description
C       IL, JL, L should contain the unit  lower triangular factor of
C       the incomplete decomposition of the A matrix  stored in SLAP
C       Row format.  IU, JU, U should contain  the unit upper factor
C       of the  incomplete decomposition of  the A matrix  stored in
C       SLAP Column format This ILU factorization can be computed by
C       the SSILUS routine. The diagonals (which are all one's) are
C       stored.
C
C       =================== S L A P Column format ==================
C
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except for  the diagonal entry, which
C       must appear first in each  "column")  and are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
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
C       must appear first in each "row") and  are stored in the real
C       array A.  In other words, for each row in the matrix put the
C       diagonal entry in  A.   Then   put  in the   other  non-zero
C       elements   going  across the  row (except   the diagonal) in
C       order.   The  JA array  holds   the column   index for  each
C       non-zero.   The IA  array holds the  offsets into  the JA, A
C       arrays  for   the   beginning  of   each  row.   That    is,
C       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
C       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
C       points to the  end of the  IROW-th row.  Note that we always
C       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
C       the matrix  and NELT  is the  number   of  non-zeros in  the
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
C***SEE ALSO  SILUR
C***REFERENCES  1. Gene Golub and Charles Van Loan, Matrix Computations,
C                  Johns Hopkins University Press, Baltimore, Maryland,
C                  1983.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   871119  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   920511  Added complete declaration section.  (WRB)
C   920929  Corrected format of reference.  (FNF)
C   930701  Updated CATEGORY section.  (FNF, WRB)
C***END PROLOGUE  SSILUS
C     .. Scalar Arguments ..
      INTEGER ISYM, N, NELT, NL, NU
C     .. Array Arguments ..
      REAL A(NELT), DINV(N), L(NL), U(NU)
      INTEGER IA(NELT), IL(NL), IU(NU), JA(NELT), JL(NL), JU(NU),
     +        NCOL(N), NROW(N)
C     .. Local Scalars ..
      REAL TEMP
      INTEGER I, IBGN, ICOL, IEND, INDX, INDX1, INDX2, INDXC1, INDXC2,
     +        INDXR1, INDXR2, IROW, ITEMP, J, JBGN, JEND, JTEMP, K, KC,
     +        KR
C***FIRST EXECUTABLE STATEMENT  SSILUS
C
C         Count number of elements in each row of the lower triangle.
C
      DO 10 I=1,N
         NROW(I) = 0
         NCOL(I) = 0
 10   CONTINUE
CVD$R NOCONCUR
CVD$R NOVECTOR
      DO 30 ICOL = 1, N
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 20 J = JBGN, JEND
               IF( IA(J).LT.ICOL ) THEN
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
                  NROW(IA(J)) = NROW(IA(J)) + 1
                  IF( ISYM.NE.0 ) NCOL(IA(J)) = NCOL(IA(J)) + 1
               ENDIF
 20         CONTINUE
         ENDIF
 30   CONTINUE
      JU(1) = 1
      IL(1) = 1
      DO 40 ICOL = 1, N
         IL(ICOL+1) = IL(ICOL) + NROW(ICOL)
         JU(ICOL+1) = JU(ICOL) + NCOL(ICOL)
         NROW(ICOL) = IL(ICOL)
         NCOL(ICOL) = JU(ICOL)
 40   CONTINUE
C
C         Copy the matrix A into the L and U structures.
      DO 60 ICOL = 1, N
         DINV(ICOL) = A(JA(ICOL))
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 50 J = JBGN, JEND
               IROW = IA(J)
               IF( IROW.LT.ICOL ) THEN
C         Part of the upper triangle.
                  IU(NCOL(ICOL)) = IROW
                  U(NCOL(ICOL)) = A(J)
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
C         Part of the lower triangle (stored by row).
                  JL(NROW(IROW)) = ICOL
                  L(NROW(IROW)) = A(J)
                  NROW(IROW) = NROW(IROW) + 1
                  IF( ISYM.NE.0 ) THEN
C         Symmetric...Copy lower triangle into upper triangle as well.
                     IU(NCOL(IROW)) = ICOL
                     U(NCOL(IROW)) = A(J)
                     NCOL(IROW) = NCOL(IROW) + 1
                  ENDIF
               ENDIF
 50         CONTINUE
         ENDIF
 60   CONTINUE
C
C         Sort the rows of L and the columns of U.
      DO 110 K = 2, N
         JBGN = JU(K)
         JEND = JU(K+1)-1
         IF( JBGN.LT.JEND ) THEN
            DO 80 J = JBGN, JEND-1
               DO 70 I = J+1, JEND
                  IF( IU(J).GT.IU(I) ) THEN
                     ITEMP = IU(J)
                     IU(J) = IU(I)
                     IU(I) = ITEMP
                     TEMP = U(J)
                     U(J) = U(I)
                     U(I) = TEMP
                  ENDIF
 70            CONTINUE
 80         CONTINUE
         ENDIF
         IBGN = IL(K)
         IEND = IL(K+1)-1
         IF( IBGN.LT.IEND ) THEN
            DO 100 I = IBGN, IEND-1
               DO 90 J = I+1, IEND
                  IF( JL(I).GT.JL(J) ) THEN
                     JTEMP = JU(I)
                     JU(I) = JU(J)
                     JU(J) = JTEMP
                     TEMP = L(I)
                     L(I) = L(J)
                     L(J) = TEMP
                  ENDIF
 90            CONTINUE
 100        CONTINUE
         ENDIF
 110  CONTINUE
C
C         Perform the incomplete LDU decomposition.
      DO 300 I=2,N
C
C           I-th row of L
         INDX1 = IL(I)
         INDX2 = IL(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 200
         DO 190 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 180
            INDXR1 = INDX1
            INDXR2 = INDX - 1
            INDXC1 = JU(JL(INDX))
            INDXC2 = JU(JL(INDX)+1) - 1
            IF(INDXC1 .GT. INDXC2) GO TO 180
 160        KR = JL(INDXR1)
 170        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 170
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 160
            ELSEIF(KR .EQ. KC) THEN
               L(INDX) = L(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 160
            ENDIF
 180        L(INDX) = L(INDX)/DINV(JL(INDX))
 190     CONTINUE
C
C         I-th column of U
 200     INDX1 = JU(I)
         INDX2 = JU(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 260
         DO 250 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 240
            INDXC1 = INDX1
            INDXC2 = INDX - 1
            INDXR1 = IL(IU(INDX))
            INDXR2 = IL(IU(INDX)+1) - 1
            IF(INDXR1 .GT. INDXR2) GO TO 240
 210        KR = JL(INDXR1)
 220        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 220
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 210
            ELSEIF(KR .EQ. KC) THEN
               U(INDX) = U(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 210
            ENDIF
 240        U(INDX) = U(INDX)/DINV(IU(INDX))
 250     CONTINUE
C
C         I-th diagonal element
 260     INDXR1 = IL(I)
         INDXR2 = IL(I+1) - 1
         IF(INDXR1 .GT. INDXR2) GO TO 300
         INDXC1 = JU(I)
         INDXC2 = JU(I+1) - 1
         IF(INDXC1 .GT. INDXC2) GO TO 300
 270     KR = JL(INDXR1)
 280     KC = IU(INDXC1)
         IF(KR .GT. KC) THEN
            INDXC1 = INDXC1 + 1
            IF(INDXC1 .LE. INDXC2) GO TO 280
         ELSEIF(KR .LT. KC) THEN
            INDXR1 = INDXR1 + 1
            IF(INDXR1 .LE. INDXR2) GO TO 270
         ELSEIF(KR .EQ. KC) THEN
            DINV(I) = DINV(I) - L(INDXR1)*DINV(KC)*U(INDXC1)
            INDXR1 = INDXR1 + 1
            INDXC1 = INDXC1 + 1
            IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 270
         ENDIF
C
 300  CONTINUE
C
C         Replace diagonal elements by their inverses.
CVD$ VECTOR
      DO 430 I=1,N
         DINV(I) = 1.0E0/DINV(I)
 430  CONTINUE
C
      RETURN
C------------- LAST LINE OF SSILUS FOLLOWS ----------------------------
      END
