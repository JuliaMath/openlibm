*DECK DLLTI2
      SUBROUTINE DLLTI2 (N, B, X, NEL, IEL, JEL, EL, DINV)
C***BEGIN PROLOGUE  DLLTI2
C***PURPOSE  SLAP Backsolve routine for LDL' Factorization.
C            Routine to solve a system of the form  L*D*L' X = B,
C            where L is a unit lower triangular matrix and D is a
C            diagonal matrix and ' means transpose.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2E
C***TYPE      DOUBLE PRECISION (SLLTI2-S, DLLTI2-D)
C***KEYWORDS  INCOMPLETE FACTORIZATION, ITERATIVE PRECONDITION, SLAP,
C             SPARSE, SYMMETRIC LINEAR SYSTEM SOLVE
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
C     DOUBLE PRECISION B(N), X(N), EL(NEL), DINV(N)
C
C     CALL DLLTI2( N, B, X, NEL, IEL, JEL, EL, DINV )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right hand side vector.
C X      :OUT      Double Precision X(N).
C         Solution to L*D*L' x = b.
C NEL    :IN       Integer.
C         Number of non-zeros in the EL array.
C IEL    :IN       Integer IEL(NEL).
C JEL    :IN       Integer JEL(NEL).
C EL     :IN       Double Precision     EL(NEL).
C         IEL, JEL, EL contain the unit lower triangular factor   of
C         the incomplete decomposition   of the A  matrix  stored in
C         SLAP Row format.   The diagonal of ones *IS* stored.  This
C         structure can be set  up  by  the DS2LT routine.  See  the
C         "Description", below for more details about the  SLAP  Row
C         format.
C DINV   :IN       Double Precision DINV(N).
C         Inverse of the diagonal matrix D.
C
C *Description:
C       This routine is supplied with  the SLAP package as a routine
C       to perform the MSOLVE operation in the SCG iteration routine
C       for  the driver  routine DSICCG.   It must be called via the
C       SLAP  MSOLVE calling sequence  convention  interface routine
C       DSLLI.
C         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****
C               **** SLAP MSOLVE CALLING CONVENTION ****
C
C       IEL, JEL, EL should contain the unit lower triangular factor
C       of  the incomplete decomposition of  the A matrix  stored in
C       SLAP Row format.   This IC factorization  can be computed by
C       the  DSICS routine.  The  diagonal  (which is all one's) is
C       stored.
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
C***SEE ALSO  DSICCG, DSICS
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
C***END PROLOGUE  DLLTI2
C     .. Scalar Arguments ..
      INTEGER N, NEL
C     .. Array Arguments ..
      DOUBLE PRECISION B(N), DINV(N), EL(NEL), X(N)
      INTEGER IEL(NEL), JEL(NEL)
C     .. Local Scalars ..
      INTEGER I, IBGN, IEND, IROW
C***FIRST EXECUTABLE STATEMENT  DLLTI2
C
C         Solve  L*y = b,  storing result in x.
C
      DO 10 I=1,N
         X(I) = B(I)
 10   CONTINUE
      DO 30 IROW = 1, N
         IBGN = IEL(IROW) + 1
         IEND = IEL(IROW+1) - 1
         IF( IBGN.LE.IEND ) THEN
CLLL. OPTION ASSERT (NOHAZARD)
CDIR$ IVDEP
CVD$ NOCONCUR
CVD$ NODEPCHK
            DO 20 I = IBGN, IEND
               X(IROW) = X(IROW) - EL(I)*X(JEL(I))
 20         CONTINUE
         ENDIF
 30   CONTINUE
C
C         Solve  D*Z = Y,  storing result in X.
C
      DO 40 I=1,N
         X(I) = X(I)*DINV(I)
 40   CONTINUE
C
C         Solve  L-trans*X = Z.
C
      DO 60 IROW = N, 2, -1
         IBGN = IEL(IROW) + 1
         IEND = IEL(IROW+1) - 1
         IF( IBGN.LE.IEND ) THEN
CLLL. OPTION ASSERT (NOHAZARD)
CDIR$ IVDEP
CVD$ NOCONCUR
CVD$ NODEPCHK
            DO 50 I = IBGN, IEND
               X(JEL(I)) = X(JEL(I)) - EL(I)*X(IROW)
 50         CONTINUE
         ENDIF
 60   CONTINUE
C
      RETURN
C------------- LAST LINE OF DLLTI2 FOLLOWS ----------------------------
      END
