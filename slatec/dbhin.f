*DECK DBHIN
      SUBROUTINE DBHIN (N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB)
C***BEGIN PROLOGUE  DBHIN
C***PURPOSE  Read a Sparse Linear System in the Boeing/Harwell Format.
C            The matrix is read in and if the right hand side is also
C            present in the input file then it too is read in.  The
C            matrix is then modified to be in the SLAP Column format.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  N1
C***TYPE      DOUBLE PRECISION (SBHIN-S, DBHIN-D)
C***KEYWORDS  LINEAR SYSTEM, MATRIX READ, SLAP SPARSE
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT, JOB
C     DOUBLE PRECISION A(NELT), SOLN(N), RHS(N)
C
C     CALL DBHIN( N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB )
C
C *Arguments:
C N      :OUT      Integer
C         Order of the Matrix.
C NELT   :INOUT    Integer.
C         On input NELT is the maximum number of non-zeros that
C         can be stored in the IA, JA, A arrays.
C         On output NELT is the number of non-zeros stored in A.
C IA     :OUT      Integer IA(NELT).
C JA     :OUT      Integer JA(NELT).
C A      :OUT      Double Precision A(NELT).
C         On output these arrays hold the matrix A in the SLAP
C         Triad format.  See "Description", below.
C ISYM   :OUT      Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C SOLN   :OUT      Double Precision SOLN(N).
C         The solution to the linear system, if present.  This array
C         is accessed if and only if JOB is set to read it in, see
C         below.  If the user requests that SOLN be read in, but it is
C         not in the file, then it is simply zeroed out.
C RHS    :OUT      Double Precision RHS(N).
C         The right hand side vector.  This array is accessed if and
C         only if JOB is set to read it in, see below.
C         If the user requests that RHS be read in, but it is not in
C         the file, then it is simply zeroed out.
C IUNIT  :IN       Integer.
C         Fortran logical I/O device unit number to read the matrix
C         from.  This unit must be connected in a system dependent
C         fashion to a file, or you will get a nasty message
C         from the Fortran I/O libraries.
C JOB    :INOUT    Integer.
C         Flag indicating what I/O operations to perform.
C         On input JOB indicates what Input operations to try to
C         perform.
C         JOB = 0 => Read only the matrix.
C         JOB = 1 => Read matrix and RHS (if present).
C         JOB = 2 => Read matrix and SOLN (if present).
C         JOB = 3 => Read matrix, RHS and SOLN (if present).
C         On output JOB indicates what operations were actually
C         performed.
C         JOB = -3 => Unable to parse matrix "CODE" from input file
C                     to determine if only the lower triangle of matrix
C                     is stored.
C         JOB = -2 => Number of non-zeros (NELT) too large.
C         JOB = -1 => System size (N) too large.
C         JOB =  0 => Read in only the matrix.
C         JOB =  1 => Read in the matrix and RHS.
C         JOB =  2 => Read in the matrix and SOLN.
C         JOB =  3 => Read in the matrix, RHS and SOLN.
C         JOB = 10 => Read in only the matrix *STRUCTURE*, but no
C                     non-zero entries.  Hence, A(*) is not referenced
C                     and has the return values the same as the input.
C         JOB = 11 => Read in the matrix *STRUCTURE* and RHS.
C         JOB = 12 => Read in the matrix *STRUCTURE* and SOLN.
C         JOB = 13 => Read in the matrix *STRUCTURE*, RHS and SOLN.
C
C *Description:
C       The format for the input is as follows.  The first line contains
C       a title to identify the data file.  On the second line (5I4) are
C       counters: NLINE, NPLS, NRILS, NNVLS, NRHSLS.
C        NLINE  Number of data lines (after the header) in the file.
C        NPLS   Number of lines for the Column Pointer data in the file.
C        NRILS  Number of lines for the Row indices in the file.
C        NNVLS  Number of lines for the Matrix elements in the file.
C        NRHSLS Number of lines for the RHS in the file.
C       The third line (A3,11X,4I4) contains a symmetry code and some
C       additional counters: CODE, NROW, NCOL, NIND, NELE.
C       On the fourth line (2A16,2A20) are formats to be used to read
C       the following data: PNTFNT, RINFMT, NVLFMT, RHSFMT.
C       Following that are the blocks of data in the order indicated.
C
C       =================== S L A P Triad format ===================
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero  matrix  element is  placed   in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C *Portability:
C         You must make sure that IUNIT is a valid Fortran logical
C         I/O device unit number and that the unit number has been
C         associated with a file or the console.  This is a system
C         dependent function.
C
C *Implementation note:
C         SOLN is not read by this version.  It will simply be
C         zeroed out if JOB = 2 or 3 and the returned value of
C         JOB will indicate SOLN has not been read.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   881107  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   911122  Added loop to zero out RHS if user wants to read RHS, but
C           it's not in the input file. (MKS)
C   911125  Minor improvements to prologue.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C   921007  Corrected description of input format.  (FNF)
C   921208  Added Implementation Note and code to zero out SOLN.  (FNF)
C   930701  Updated CATEGORY section.  (FNF, WRB)
C***END PROLOGUE  DBHIN
C     .. Scalar Arguments ..
      INTEGER ISYM, IUNIT, JOB, N, NELT
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), RHS(N), SOLN(N)
      INTEGER IA(NELT), JA(NELT)
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I, IBGN, ICOL, IEND, ITEMP, J, JOBRET, NCOL, NELE, NIND,
     +        NLINE, NNVLS, NPLS, NRHSLS, NRILS, NROW
      CHARACTER CODE*3, PNTFMT*16, RINFMT*16, NVLFMT*20, RHSFMT*20,
     +          TITLE*80
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C***FIRST EXECUTABLE STATEMENT  DBHIN
C
C         Read Matrices In BOEING-HARWELL format.
C
C TITLE  Header line to identify data file.
C NLINE  Number of data lines (after the header) in the file.
C NPLS   Number of lines for the Column Pointer data in the file.
C NRILS  Number of lines for the Row indices in the data file.
C NNVLS  Number of lines for the Matrix elements in the data file.
C NRHSLS Number of lines for the RHS in the data file.
C ---- Only those variables needed by SLAP are referenced. ----
C
      READ(IUNIT,9000) TITLE
      READ(IUNIT,9010) NLINE, NPLS, NRILS, NNVLS, NRHSLS
      READ(IUNIT,9020) CODE, NROW, NCOL, NIND, NELE
      READ(IUNIT,9030) PNTFMT, RINFMT, NVLFMT, RHSFMT
C
      IF( NROW.GT.N ) THEN
         N = NROW
         JOBRET = -1
         GOTO 999
      ENDIF
      IF( NIND.GT.NELT ) THEN
         NELT = NIND
         JOBRET = -2
         GOTO 999
      ENDIF
C
C         Set the parameters.
C
      N    = NROW
      NELT = NIND
      IF( CODE.EQ.'RUA' ) THEN
         ISYM = 0
      ELSE IF( CODE.EQ.'RSA' ) THEN
         ISYM = 1
      ELSE
         JOBRET = -3
         GOTO 999
      ENDIF
      READ(IUNIT,PNTFMT) (JA(I), I = 1, N+1)
      READ(IUNIT,RINFMT) (IA(I), I = 1, NELT)
      JOBRET = 10
      IF( NNVLS.GT.0 ) THEN
         READ(IUNIT,NVLFMT) (A(I),  I = 1, NELT)
         JOBRET = 0
      ENDIF
      IF( MOD(JOB,2).EQ.1 ) THEN
C
C         User requests that the RHS be read in.  If it is in the input
C         file, read it in; otherwise just zero it out.
C
         IF( NRHSLS.GT.0 ) THEN
            READ(5,RHSFMT) (RHS(I), I = 1, N)
            JOBRET = JOBRET + 1
         ELSE
            DO 10 I = 1, N
               RHS(I) = 0
 10         CONTINUE
         ENDIF
      ENDIF
      IF ( (JOB.EQ.2).OR.(JOB.EQ.3) ) THEN
C
C         User requests that the SOLN be read in.
C         Just zero out the array.
C
         DO 20 I = 1, N
            SOLN(I) = 0
 20      CONTINUE
      ENDIF
C
C         Now loop through the IA array making sure that the diagonal
C         matrix element appears first in the column.  Then sort the
C         rest of the column in ascending order.
C
CVD$R NOCONCUR
CVD$R NOVECTOR
      DO 70 ICOL = 1, N
         IBGN = JA(ICOL)
         IEND = JA(ICOL+1)-1
         DO 30 I = IBGN, IEND
            IF( IA(I).EQ.ICOL ) THEN
C
C              Swap the diagonal element with the first element in the
C              column.
C
               ITEMP = IA(I)
               IA(I) = IA(IBGN)
               IA(IBGN) = ITEMP
               TEMP = A(I)
               A(I) = A(IBGN)
               A(IBGN) = TEMP
               GOTO 40
            ENDIF
 30      CONTINUE
 40      IBGN = IBGN + 1
         IF( IBGN.LT.IEND ) THEN
            DO 60 I = IBGN, IEND
               DO 50 J = I+1, IEND
                  IF( IA(I).GT.IA(J) ) THEN
                     ITEMP = IA(I)
                     IA(I) = IA(J)
                     IA(J) = ITEMP
                     TEMP = A(I)
                     A(I) = A(J)
                     A(J) = TEMP
                  ENDIF
 50            CONTINUE
 60         CONTINUE
         ENDIF
 70   CONTINUE
C
C         Set return flag.
 999  JOB = JOBRET
      RETURN
 9000 FORMAT( A80 )
 9010 FORMAT( 5I14 )
 9020 FORMAT( A3, 11X, 4I14 )
 9030 FORMAT( 2A16, 2A20 )
C------------- LAST LINE OF DBHIN FOLLOWS ------------------------------
      END
