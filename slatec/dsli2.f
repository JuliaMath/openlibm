*DECK DSLI2
      SUBROUTINE DSLI2 (N, B, X, NEL, IEL, JEL, EL)
C***BEGIN PROLOGUE  DSLI2
C***PURPOSE  SLAP Lower Triangle Matrix Backsolve.
C            Routine to solve a system of the form  Lx = b , where L
C            is a lower triangular matrix.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2A3
C***TYPE      DOUBLE PRECISION (SSLI2-S, DSLI2-D)
C***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE
C***AUTHOR  Greenbaum, Anne, (Courant Institute)
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C *Usage:
C     INTEGER N, NEL, IEL(NEL), JEL(NEL)
C     DOUBLE PRECISION B(N), X(N), EL(NEL)
C
C     CALL DSLI2( N, B, X, NEL, IEL, JEL, EL )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right hand side vector.
C X      :OUT      Double Precision X(N).
C         Solution to Lx = b.
C NEL    :IN       Integer.
C         Number of non-zeros in the EL array.
C IEL    :IN       Integer IEL(NEL).
C JEL    :IN       Integer JEL(NEL).
C EL     :IN       Double Precision EL(NEL).
C         IEL, JEL, EL contain the unit lower triangular factor   of
C         the incomplete decomposition   of the A  matrix  stored in
C         SLAP Row format.  The diagonal of  ones *IS* stored.  This
C         structure can be  set up by the  DS2LT  routine.  See  the
C         "Description", below, for more details about the  SLAP Row
C         format.
C
C *Description:
C       This routine is supplied with the SLAP package  as a routine
C       to perform the MSOLVE operation in the DIR iteration routine
C       for the driver routine DSGS.  It must be called via the SLAP
C       MSOLVE calling sequence convention interface routine DSLI.
C         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****
C               **** SLAP MSOLVE CALLING CONVENTION ****
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
C       With  the SLAP  Row format  the "inner loop" of this routine
C       should vectorize   on machines with   hardware  support  for
C       vector gather/scatter operations.  Your compiler may require
C       a  compiler directive  to  convince   it that there  are  no
C       implicit vector  dependencies.  Compiler directives  for the
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
C       with the standard SLAP distribution.
C
C***SEE ALSO  DSLI
C***REFERENCES  (NONE)
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
C   921113  Corrected C***CATEGORY line.  (FNF)
C   930701  Updated CATEGORY section.  (FNF, WRB)
C***END PROLOGUE  DSLI2
C     .. Scalar Arguments ..
      INTEGER N, NEL
C     .. Array Arguments ..
      DOUBLE PRECISION B(N), EL(NEL), X(N)
      INTEGER IEL(NEL), JEL(NEL)
C     .. Local Scalars ..
      INTEGER I, ICOL, J, JBGN, JEND
C***FIRST EXECUTABLE STATEMENT  DSLI2
C
C         Initialize the solution by copying the right hands side
C         into it.
C
      DO 10 I=1,N
         X(I) = B(I)
 10   CONTINUE
C
CVD$ NOCONCUR
      DO 30 ICOL = 1, N
         X(ICOL) = X(ICOL)/EL(JEL(ICOL))
         JBGN = JEL(ICOL) + 1
         JEND = JEL(ICOL+1) - 1
         IF( JBGN.LE.JEND ) THEN
CLLL. OPTION ASSERT (NOHAZARD)
CDIR$ IVDEP
CVD$ NOCONCUR
CVD$ NODEPCHK
            DO 20 J = JBGN, JEND
               X(IEL(J)) = X(IEL(J)) - EL(J)*X(ICOL)
 20         CONTINUE
         ENDIF
 30   CONTINUE
C
      RETURN
C------------- LAST LINE OF DSLI2 FOLLOWS ----------------------------
      END
