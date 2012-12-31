*DECK STIN
      SUBROUTINE STIN (N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB)
C***BEGIN PROLOGUE  STIN
C***PURPOSE  Read in SLAP Triad Format Linear System.
C            Routine to read in a SLAP Triad format matrix and right
C            hand side and solution to the system, if known.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  N1
C***TYPE      SINGLE PRECISION (STIN-S, DTIN-D)
C***KEYWORDS  DIAGNOSTICS, LINEAR SYSTEM, SLAP SPARSE
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT, JOB
C     REAL    A(NELT), SOLN(N), RHS(N)
C
C     CALL STIN( N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB )
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
C A      :OUT      Real A(NELT).
C         On output these arrays hold the matrix A in the SLAP
C         Triad format.  See "Description", below.
C ISYM   :OUT      Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C SOLN   :OUT      Real SOLN(N).
C         The solution to the linear system, if present.  This array
C         is accessed if and only if JOB to read it in, see below.
C         If the user requests that SOLN be read in, but it is not in
C         the file, then it is simply zeroed out.
C RHS    :OUT      Real RHS(N).
C         The right hand side vector.  This array is accessed if and
C         only if JOB is set to read it in, see below.
C         If the user requests that RHS be read in, but it is not in
C         the file, then it is simply zeroed out.
C IUNIT  :IN       Integer.
C         Fortran logical I/O device unit number to write the matrix
C         to.  This unit must be connected in a system dependent fashion
C         to a file or the console or you will get a nasty message
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
C         JOB = 0 => Read in only the matrix.
C         JOB = 1 => Read in the matrix and RHS.
C         JOB = 2 => Read in the matrix and SOLN.
C         JOB = 3 => Read in the matrix, RHS and SOLN.
C
C *Description:
C       The format for the  input is as follows.  On  the first line
C       are counters and flags: N, NELT, ISYM, IRHS, ISOLN.  N, NELT
C       and ISYM are described above.  IRHS is  a flag indicating if
C       the RHS was  written out (1 is  yes, 0 is  no).  ISOLN  is a
C       flag indicating if the SOLN was written out  (1 is yes, 0 is
C       no).  The format for the fist line is: 5i10.  Then comes the
C       NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format
C       for  these lines is   :  1X,I5,1X,I5,1X,E16.7.   Then  comes
C       RHS(I), I = 1, N, if IRHS = 1.  Then  comes SOLN(I), I  = 1,
C       N, if ISOLN = 1.  The format for these lines is: 1X,E16.7.
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
C *Cautions:
C     This routine will attempt to write to the Fortran logical output
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
C     this logical unit is attached to a file or terminal before calling
C     this routine with a non-zero value for IUNIT.  This routine does
C     not check for the validity of a non-zero IUNIT unit number.
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
C   930701  Updated CATEGORY section.  (FNF, WRB)
C***END PROLOGUE  STIN
C     .. Scalar Arguments ..
      INTEGER ISYM, IUNIT, JOB, N, NELT
C     .. Array Arguments ..
      REAL A(NELT), RHS(N), SOLN(N)
      INTEGER IA(NELT), JA(NELT)
C     .. Local Scalars ..
      INTEGER I, IRHS, ISOLN, JOBRET, NELTMX
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C***FIRST EXECUTABLE STATEMENT  STIN
C
C         Read in the information heading.
C
      NELTMX = NELT
      READ(IUNIT,1000) N, NELT, ISYM, IRHS, ISOLN
      NELT = MIN( NELT, NELTMX )
C
C         Read in the matrix non-zeros in Triad format.
      DO 10 I = 1, NELT
         READ(IUNIT,1010) IA(I), JA(I), A(I)
 10   CONTINUE
C
C         If requested, read in the rhs.
      JOBRET = 0
      IF( JOB.EQ.1 .OR. JOB.EQ.3 ) THEN
C
C         Check to see if rhs is in the file.
         IF( IRHS.EQ.1 ) THEN
            JOBRET = 1
            READ(IUNIT,1020) (RHS(I),I=1,N)
         ELSE
            DO 20 I = 1, N
               RHS(I) = 0
 20         CONTINUE
         ENDIF
      ENDIF
C
C         If requested, read in the solution.
      IF( JOB.GT.1 ) THEN
C
C         Check to see if solution is in the file.
         IF( ISOLN.EQ.1 ) THEN
            JOBRET = JOBRET + 2
            READ(IUNIT,1020) (SOLN(I),I=1,N)
         ELSE
            DO 30 I = 1, N
               SOLN(I) = 0
 30         CONTINUE
         ENDIF
      ENDIF
C
      JOB = JOBRET
      RETURN
 1000 FORMAT(5I10)
 1010 FORMAT(1X,I5,1X,I5,1X,E16.7)
 1020 FORMAT(1X,E16.7)
C------------- LAST LINE OF STIN FOLLOWS ----------------------------
      END
