*DECK DCPPLT
      SUBROUTINE DCPPLT (N, NELT, IA, JA, A, ISYM, IUNIT)
C***BEGIN PROLOGUE  DCPPLT
C***PURPOSE  Printer Plot of SLAP Column Format Matrix.
C            Routine to print out a SLAP Column format matrix in a
C            "printer plot" graphical representation.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  N1
C***TYPE      DOUBLE PRECISION (SCPPLT-S, DCPPLT-D)
C***KEYWORDS  DIAGNOSTICS, LINEAR SYSTEM, SLAP SPARSE
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT
C     DOUBLE PRECISION A(NELT)
C
C     CALL DCPPLT( N, NELT, IA, JA, A, ISYM, IUNIT )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C         If N.gt.MAXORD, only the leading MAXORD x MAXORD
C         submatrix will be printed.  (Currently MAXORD = 225.)
C NELT   :IN       Integer.
C         Number of non-zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays should hold the matrix A in the SLAP
C         Column format.  See "Description", below.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all non-zero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C IUNIT  :IN       Integer.
C         Fortran logical I/O device unit number to write the matrix
C         to.  This unit must be connected in a system dependent fashion
C         to a file or the console or you will get a nasty message
C         from the Fortran I/O libraries.
C
C *Description:
C       This routine prints out a SLAP  Column format matrix  to the
C       Fortran logical I/O unit   number  IUNIT.  The  numbers them
C       selves  are not printed  out, but   rather  a one  character
C       representation of the numbers.   Elements of the matrix that
C       are not represented in the (IA,JA,A)  arrays are  denoted by
C       ' ' character (a blank).  Elements of A that are *ZERO* (and
C       hence  should  really not be  stored) are  denoted  by a '0'
C       character.  Elements of A that are *POSITIVE* are denoted by
C       'D' if they are Diagonal elements  and '#' if  they are off
C       Diagonal  elements.  Elements of  A that are *NEGATIVE* are
C       denoted by 'N'  if they  are Diagonal  elements and  '*' if
C       they are off Diagonal elements.
C
C       =================== S L A P Column format ==================
C
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
C *Cautions:
C     This routine will attempt to write to the Fortran logical output
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
C     this logical unit is attached to a file or terminal before calling
C     this routine with a non-zero value for IUNIT.  This routine does
C     not check for the validity of a non-zero IUNIT unit number.
C
C *Portability:
C     This routine, as distributed, can generate lines up to 229
C     characters long.  Some Fortran systems have more restricted
C     line lengths.  Change parameter MAXORD and the large number
C     in FORMAT 1010 to reduce this line length.
C
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
C   921007  Replaced hard-wired 225 with parameter MAXORD.  (FNF)
C   921021  Corrected syntax of CHARACTER declaration.  (FNF)
C   921026  Corrected D to E in output format.  (FNF)
C   930701  Updated CATEGORY section.  (FNF, WRB)
C***END PROLOGUE  DCPPLT
C     .. Scalar Arguments ..
      INTEGER ISYM, IUNIT, N, NELT
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT)
      INTEGER IA(NELT), JA(NELT)
C     .. Parameters ..
      INTEGER  MAXORD
      PARAMETER (MAXORD=225)
C     .. Local Scalars ..
      INTEGER I, ICOL, IROW, J, JBGN, JEND, NMAX
C     .. Local Arrays ..
      CHARACTER CHMAT(MAXORD)*(MAXORD)
C     .. Intrinsic Functions ..
      INTRINSIC MIN, MOD, REAL
C***FIRST EXECUTABLE STATEMENT  DCPPLT
C
C         Set up the character matrix...
C
      NMAX = MIN( MAXORD, N )
      DO 10 I = 1, NMAX
         CHMAT(I)(1:NMAX) = ' '
 10   CONTINUE
      DO 30 ICOL = 1, NMAX
         JBGN = JA(ICOL)
         JEND = JA(ICOL+1)-1
         DO 20 J = JBGN, JEND
            IROW = IA(J)
            IF( IROW.LE.NMAX ) THEN
               IF( ISYM.NE.0 ) THEN
C         Put in non-sym part as well...
                  IF( A(J).EQ.0.0D0 ) THEN
                     CHMAT(IROW)(ICOL:ICOL) = '0'
                  ELSEIF( A(J).GT.0.0D0 ) THEN
                     CHMAT(IROW)(ICOL:ICOL) = '#'
                  ELSE
                     CHMAT(IROW)(ICOL:ICOL) = '*'
                  ENDIF
               ENDIF
               IF( IROW.EQ.ICOL ) THEN
C         Diagonal entry.
                  IF( A(J).EQ.0.0D0 ) THEN
                     CHMAT(IROW)(ICOL:ICOL) = '0'
                  ELSEIF( A(J).GT.0.0D0 ) THEN
                     CHMAT(IROW)(ICOL:ICOL) = 'D'
                  ELSE
                     CHMAT(IROW)(ICOL:ICOL) = 'N'
                  ENDIF
               ELSE
C         Off-Diagonal entry
                  IF( A(J).EQ.0.0D0 ) THEN
                     CHMAT(IROW)(ICOL:ICOL) = '0'
                  ELSEIF( A(J).GT.0.0D0 ) THEN
                     CHMAT(IROW)(ICOL:ICOL) = '#'
                  ELSE
                     CHMAT(IROW)(ICOL:ICOL) = '*'
                  ENDIF
               ENDIF
            ENDIF
 20      CONTINUE
 30   CONTINUE
C
C         Write out the heading.
      WRITE(IUNIT,1000) N, NELT, REAL(NELT)/(N*N)
      WRITE(IUNIT,1010) (MOD(I,10),I=1,NMAX)
C
C         Write out the character representations matrix elements.
      DO 40 IROW = 1, NMAX
         WRITE(IUNIT,1020) IROW, CHMAT(IROW)(1:NMAX)
 40   CONTINUE
      RETURN
C
 1000 FORMAT(/'**** Picture of Column SLAP matrix follows ****'/
     $     ' N, NELT and Density = ',2I10,D16.7)
C      The following assumes MAXORD.le.225.
 1010 FORMAT(4X,225(I1))
 1020 FORMAT(1X,I3,A)
C------------- LAST LINE OF DCPPLT FOLLOWS ----------------------------
      END
