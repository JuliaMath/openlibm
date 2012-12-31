*DECK SPLPUP
      SUBROUTINE SPLPUP (USRMAT, MRELAS, NVARS, PRGOPT, DATTRV, BL, BU,
     +   IND, INFO, AMAT, IMAT, SIZEUP, ASMALL, ABIG)
C***BEGIN PROLOGUE  SPLPUP
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SPLP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SPLPUP-S, DPLPUP-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     /REAL (12 BLANKS)/DOUBLE PRECISION/.
C
C     REVISED 810613-1130
C     REVISED YYMMDD-HHMM
C
C     THIS SUBROUTINE COLLECTS INFORMATION ABOUT THE BOUNDS AND MATRIX
C     FROM THE USER.  IT IS PART OF THE SPLP( ) PACKAGE.
C
C***SEE ALSO  SPLP
C***ROUTINES CALLED  PCHNGS, PNNZRS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   811215  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890605  Corrected references to XERRWV.  (WRB)
C   890605  Removed unreferenced labels.  (WRB)
C   891009  Removed unreferenced variables.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls, changed do-it-yourself
C           DO loops to DO loops.  (RWC)
C   900602  Get rid of ASSIGNed GOTOs.  (RWC)
C***END PROLOGUE  SPLPUP
      REAL             ABIG,AIJ,AMAT(*),AMN,AMX,ASMALL,BL(*),
     * BU(*),DATTRV(*),PRGOPT(*),XVAL,ZERO
      INTEGER IFLAG(10),IMAT(*),IND(*)
      LOGICAL SIZEUP,FIRST
      CHARACTER*8 XERN1, XERN2
      CHARACTER*16 XERN3, XERN4
C
C***FIRST EXECUTABLE STATEMENT  SPLPUP
      ZERO = 0.E0
C
C     CHECK USER-SUPPLIED BOUNDS
C
C     CHECK THAT IND(*) VALUES ARE 1,2,3 OR 4.
C     ALSO CHECK CONSISTENCY OF UPPER AND LOWER BOUNDS.
C
      DO 10 J=1,NVARS
         IF (IND(J).LT.1 .OR. IND(J).GT.4) THEN
            WRITE (XERN1, '(I8)') J
            CALL XERMSG ('SLATEC', 'SPLPUP',
     *         'IN SPLP, INDEPENDENT VARIABLE = ' // XERN1 //
     *         ' IS NOT DEFINED.', 10, 1)
            INFO = -10
            RETURN
         ENDIF
C
         IF (IND(J).EQ.3) THEN
            IF (BL(J).GT.BU(J)) THEN
               WRITE (XERN1, '(I8)') J
               WRITE (XERN3, '(1PE15.6)') BL(J)
               WRITE (XERN4, '(1PE15.6)') BU(J)
               CALL XERMSG ('SLATEC', 'SPLPUP',
     *            'IN SPLP, LOWER BOUND = ' // XERN3 //
     *            ' AND UPPER BOUND = ' // XERN4 //
     *            ' FOR INDEPENDENT VARIABLE = ' // XERN1 //
     *            ' ARE NOT CONSISTENT.', 11, 1)
               RETURN
            ENDIF
         ENDIF
   10 CONTINUE
C
      DO 20 I=NVARS+1,NVARS+MRELAS
         IF (IND(I).LT.1 .OR. IND(I).GT.4) THEN
            WRITE (XERN1, '(I8)') I-NVARS
            CALL XERMSG ('SLATEC', 'SPLPUP',
     *         'IN SPLP, DEPENDENT VARIABLE = ' // XERN1 //
     *         ' IS NOT DEFINED.', 12, 1)
            INFO = -12
            RETURN
         ENDIF
C
         IF (IND(I).EQ.3) THEN
            IF (BL(I).GT.BU(I)) THEN
               WRITE (XERN1, '(I8)') I
               WRITE (XERN3, '(1PE15.6)') BL(I)
               WRITE (XERN4, '(1PE15.6)') BU(I)
               CALL XERMSG ('SLATEC', 'SPLPUP',
     *            'IN SPLP, LOWER BOUND = ' // XERN3 //
     *            ' AND UPPER BOUND = ' // XERN4 //
     *            ' FOR DEPENDANT VARIABLE = ' // XERN1 //
     *            ' ARE NOT CONSISTENT.',13,1)
               INFO = -13
               RETURN
            ENDIF
         ENDIF
   20 CONTINUE
C
C     GET UPDATES OR DATA FOR MATRIX FROM THE USER
C
C     GET THE ELEMENTS OF THE MATRIX FROM THE USER.  IT WILL BE STORED
C     BY COLUMNS USING THE SPARSE STORAGE CODES OF RJ HANSON AND
C     JA WISNIEWSKI.
C
      IFLAG(1) = 1
C
C     KEEP ACCEPTING ELEMENTS UNTIL THE USER IS FINISHED GIVING THEM.
C     LIMIT THIS LOOP TO 2*NVARS*MRELAS ITERATIONS.
C
      ITMAX = 2*NVARS*MRELAS+1
      ITCNT = 0
      FIRST = .TRUE.
C
C     CHECK ON THE ITERATION COUNT.
C
   30 ITCNT = ITCNT+1
      IF (ITCNT.GT.ITMAX) THEN
         CALL XERMSG ('SLATEC', 'SPLPUP',
     +      'IN SPLP, MORE THAN 2*NVARS*MRELAS ITERATIONS DEFINING ' //
     +      'OR UPDATING MATRIX DATA.', 7, 1)
         INFO = -7
         RETURN
      ENDIF
C
      AIJ = ZERO
      CALL USRMAT(I,J,AIJ,INDCAT,PRGOPT,DATTRV,IFLAG)
      IF (IFLAG(1).EQ.1) THEN
         IFLAG(1) = 2
         GO TO 30
      ENDIF
C
C     CHECK TO SEE THAT THE SUBSCRIPTS I AND J ARE VALID.
C
      IF (I.LT.1 .OR. I.GT.MRELAS .OR. J.LT.1 .OR. J.GT.NVARS) THEN
C
C        CHECK ON SIZE OF MATRIX DATA
C        RECORD THE LARGEST AND SMALLEST(IN MAGNITUDE) NONZERO ELEMENTS.
C
         IF (IFLAG(1).EQ.3) THEN
            IF (SIZEUP .AND. ABS(AIJ).NE.ZERO) THEN
               IF (FIRST) THEN
                  AMX = ABS(AIJ)
                  AMN = ABS(AIJ)
                  FIRST = .FALSE.
               ELSEIF (ABS(AIJ).GT.AMX) THEN
                  AMX = ABS(AIJ)
               ELSEIF (ABS(AIJ).LT.AMN) THEN
                  AMN = ABS(AIJ)
               ENDIF
            ENDIF
            GO TO 40
         ENDIF
C
         WRITE (XERN1, '(I8)') I
         WRITE (XERN2, '(I8)') J
         CALL XERMSG ('SLATEC', 'SPLPUP',
     *      'IN SPLP, ROW INDEX = ' // XERN1 // ' OR COLUMN INDEX = '
     *      // XERN2 // ' IS OUT OF RANGE.', 8, 1)
         INFO = -8
         RETURN
      ENDIF
C
C     IF INDCAT=0 THEN SET A(I,J)=AIJ.
C     IF INDCAT=1 THEN ACCUMULATE ELEMENT, A(I,J)=A(I,J)+AIJ.
C
      IF (INDCAT.EQ.0) THEN
         CALL PCHNGS(I,AIJ,IPLACE,AMAT,IMAT,J)
      ELSEIF (INDCAT.EQ.1) THEN
         INDEX = -(I-1)
         CALL PNNZRS(INDEX,XVAL,IPLACE,AMAT,IMAT,J)
         IF (INDEX.EQ.I) AIJ=AIJ+XVAL
         CALL PCHNGS(I,AIJ,IPLACE,AMAT,IMAT,J)
      ELSE
         WRITE (XERN1, '(I8)') INDCAT
         CALL XERMSG ('SLATEC', 'SPLPUP',
     *      'IN SPLP, INDICATION FLAG = ' // XERN1 //
     *      ' FOR MATRIX DATA MUST BE EITHER 0 OR 1.', 9, 1)
         INFO = -9
         RETURN
      ENDIF
C
C     CHECK ON SIZE OF MATRIX DATA
C     RECORD THE LARGEST AND SMALLEST(IN MAGNITUDE) NONZERO ELEMENTS.
C
      IF (SIZEUP .AND. ABS(AIJ).NE.ZERO) THEN
         IF (FIRST) THEN
            AMX = ABS(AIJ)
            AMN = ABS(AIJ)
            FIRST = .FALSE.
         ELSEIF (ABS(AIJ).GT.AMX) THEN
            AMX = ABS(AIJ)
         ELSEIF (ABS(AIJ).LT.AMN) THEN
            AMN = ABS(AIJ)
         ENDIF
      ENDIF
      IF (IFLAG(1).NE.3) GO TO 30
C
   40 IF (SIZEUP .AND. .NOT. FIRST) THEN
         IF (AMN.LT.ASMALL .OR. AMX.GT.ABIG) THEN
            CALL XERMSG ('SLATEC', 'SPLPUP',
     +         'IN SPLP, A MATRIX ELEMENT''S SIZE IS OUT OF THE ' //
     +         'SPECIFIED RANGE.', 22, 1)
            INFO = -22
            RETURN
         ENDIF
      ENDIF
      RETURN
      END
