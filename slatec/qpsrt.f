*DECK QPSRT
      SUBROUTINE QPSRT (LIMIT, LAST, MAXERR, ERMAX, ELIST, IORD, NRMAX)
C***BEGIN PROLOGUE  QPSRT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to QAGE, QAGIE, QAGPE, QAGSE, QAWCE, QAWOE and
C            QAWSE
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (QPSRT-S, DQPSRT-D)
C***KEYWORDS  SEQUENTIAL SORTING
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C 1.        QPSRT
C           Ordering Routine
C              Standard FORTRAN Subroutine
C              REAL Version
C
C 2.        PURPOSE
C              This routine maintains the descending ordering
C              in the list of the local error estimates resulting from
C              the interval subdivision process. At each call two error
C              estimates are inserted using the sequential search
C              method, top-down for the largest error estimate
C              and bottom-up for the smallest error estimate.
C
C 3.        CALLING SEQUENCE
C              CALL QPSRT(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)
C
C           PARAMETERS (MEANING AT OUTPUT)
C              LIMIT  - INTEGER
C                       Maximum number of error estimates the list
C                       can contain
C
C              LAST   - INTEGER
C                       Number of error estimates currently
C                       in the list
C
C              MAXERR - INTEGER
C                       MAXERR points to the NRMAX-th largest error
C                       estimate currently in the list
C
C              ERMAX  - REAL
C                       NRMAX-th largest error estimate
C                       ERMAX = ELIST(MAXERR)
C
C              ELIST  - REAL
C                       Vector of dimension LAST containing
C                       the error estimates
C
C              IORD   - INTEGER
C                       Vector of dimension LAST, the first K
C                       elements of which contain pointers
C                       to the error estimates, such that
C                       ELIST(IORD(1)),... , ELIST(IORD(K))
C                       form a decreasing sequence, with
C                       K = LAST if LAST.LE.(LIMIT/2+2), and
C                       K = LIMIT+1-LAST otherwise
C
C              NRMAX  - INTEGER
C                       MAXERR = IORD(NRMAX)
C
C***SEE ALSO  QAGE, QAGIE, QAGPE, QAGSE, QAWCE, QAWOE, QAWSE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  QPSRT
C
      REAL ELIST,ERMAX,ERRMAX,ERRMIN
      INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,
     1  NRMAX
      DIMENSION ELIST(*),IORD(*)
C
C           CHECK WHETHER THE LIST CONTAINS MORE THAN
C           TWO ERROR ESTIMATES.
C
C***FIRST EXECUTABLE STATEMENT  QPSRT
      IF(LAST.GT.2) GO TO 10
      IORD(1) = 1
      IORD(2) = 2
      GO TO 90
C
C           THIS PART OF THE ROUTINE IS ONLY EXECUTED
C           IF, DUE TO A DIFFICULT INTEGRAND, SUBDIVISION
C           INCREASED THE ERROR ESTIMATE. IN THE NORMAL CASE
C           THE INSERT PROCEDURE SHOULD START AFTER THE
C           NRMAX-TH LARGEST ERROR ESTIMATE.
C
   10 ERRMAX = ELIST(MAXERR)
      IF(NRMAX.EQ.1) GO TO 30
      IDO = NRMAX-1
      DO 20 I = 1,IDO
        ISUCC = IORD(NRMAX-1)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
        IORD(NRMAX) = ISUCC
        NRMAX = NRMAX-1
   20    CONTINUE
C
C           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO
C           BE MAINTAINED IN DESCENDING ORDER. THIS NUMBER
C           DEPENDS ON THE NUMBER OF SUBDIVISIONS STILL
C           ALLOWED.
C
   30 JUPBN = LAST
      IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
      ERRMIN = ELIST(LAST)
C
C           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
C           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
C
      JBND = JUPBN-1
      IBEG = NRMAX+1
      IF(IBEG.GT.JBND) GO TO 50
      DO 40 I=IBEG,JBND
        ISUCC = IORD(I)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
        IORD(I-1) = ISUCC
   40 CONTINUE
   50 IORD(JBND) = MAXERR
      IORD(JUPBN) = LAST
      GO TO 90
C
C           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
C
   60 IORD(I-1) = MAXERR
      K = JBND
      DO 70 J=I,JBND
        ISUCC = IORD(K)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
        IORD(K+1) = ISUCC
        K = K-1
   70 CONTINUE
      IORD(I) = LAST
      GO TO 90
   80 IORD(K+1) = LAST
C
C           SET MAXERR AND ERMAX.
C
   90 MAXERR = IORD(NRMAX)
      ERMAX = ELIST(MAXERR)
      RETURN
      END
