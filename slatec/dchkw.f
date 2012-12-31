*DECK DCHKW
      SUBROUTINE DCHKW (NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR)
C***BEGIN PROLOGUE  DCHKW
C***SUBSIDIARY
C***PURPOSE  SLAP WORK/IWORK Array Bounds Checker.
C            This routine checks the work array lengths and interfaces
C            to the SLATEC error handler if a problem is found.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  R2
C***TYPE      DOUBLE PRECISION (SCHKW-S, DCHKW-D)
C***KEYWORDS  ERROR CHECKING, SLAP, WORKSPACE CHECKING
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C
C *Usage:
C     CHARACTER*(*) NAME
C     INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER
C     DOUBLE PRECISION ERR
C
C     CALL DCHKW( NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
C
C *Arguments:
C NAME   :IN       Character*(*).
C         Name of the calling routine.  This is used in the output
C         message, if an error is detected.
C LOCIW  :IN       Integer.
C         Location of the first free element in the integer workspace
C         array.
C LENIW  :IN       Integer.
C         Length of the integer workspace array.
C LOCW   :IN       Integer.
C         Location of the first free element in the double precision
C         workspace array.
C LENRW  :IN       Integer.
C         Length of the double precision workspace array.
C IERR   :OUT      Integer.
C         Return error flag.
C               IERR = 0 => All went well.
C               IERR = 1 => Insufficient storage allocated for
C                           WORK or IWORK.
C ITER   :OUT      Integer.
C         Set to zero on return.
C ERR    :OUT      Double Precision.
C         Set to the smallest positive magnitude if all went well.
C         Set to a very large number if an error is detected.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   880225  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   900805  Changed XERRWV calls to calls to XERMSG.  (RWC)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910502  Corrected XERMSG calls to satisfy Section 6.2.2 of ANSI
C           X3.9-1978.  (FNF)
C   910506  Made subsidiary.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C   921015  Added code to initialize ITER and ERR when IERR=0.  (FNF)
C***END PROLOGUE  DCHKW
C     .. Scalar Arguments ..
      DOUBLE PRECISION ERR
      INTEGER IERR, ITER, LENIW, LENW, LOCIW, LOCW
      CHARACTER NAME*(*)
C     .. Local Scalars ..
      CHARACTER XERN1*8, XERN2*8, XERNAM*8
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
C     .. External Subroutines ..
      EXTERNAL XERMSG
C***FIRST EXECUTABLE STATEMENT  DCHKW
C
C         Check the Integer workspace situation.
C
      IERR = 0
      ITER = 0
      ERR = D1MACH(1)
      IF( LOCIW.GT.LENIW ) THEN
         IERR = 1
         ERR = D1MACH(2)
         XERNAM = NAME
         WRITE (XERN1, '(I8)') LOCIW
         WRITE (XERN2, '(I8)') LENIW
         CALL XERMSG ('SLATEC', 'DCHKW',
     $      'In ' // XERNAM // ', INTEGER work array too short.  ' //
     $      'IWORK needs ' // XERN1 // '; have allocated ' // XERN2,
     $      1, 1)
      ENDIF
C
C         Check the Double Precision workspace situation.
      IF( LOCW.GT.LENW ) THEN
         IERR = 1
         ERR = D1MACH(2)
         XERNAM = NAME
         WRITE (XERN1, '(I8)') LOCW
         WRITE (XERN2, '(I8)') LENW
         CALL XERMSG ('SLATEC', 'DCHKW',
     $      'In ' // XERNAM // ', DOUBLE PRECISION work array too ' //
     $      'short.  RWORK needs ' // XERN1 // '; have allocated ' //
     $      XERN2, 1, 1)
      ENDIF
      RETURN
C------------- LAST LINE OF DCHKW FOLLOWS ----------------------------
      END
