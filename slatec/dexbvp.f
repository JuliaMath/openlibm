*DECK DEXBVP
      SUBROUTINE DEXBVP (Y, NROWY, XPTS, A, NROWA, ALPHA, B, NROWB,
     +   BETA, IFLAG, WORK, IWORK)
C***BEGIN PROLOGUE  DEXBVP
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (EXBVP-S, DEXBVP-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C  This subroutine is used to execute the basic technique for solving
C  the two-point boundary value problem.
C
C***SEE ALSO  DBVSUP
C***ROUTINES CALLED  DBVPOR, XERMSG
C***COMMON BLOCKS    DML15T, DML17B, DML18J, DML5MC, DML8SZ
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   890921  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DEXBVP
C
      INTEGER ICOCO, IEXP, IFLAG, IGOFX, INC, INDPVT, INFO, INHOMO,
     1     INTEG, ISTKOP, IVP, IWORK(*), K1, K10, K11, K2, K3,
     2     K4, K5, K6, K7, K8, K9, KKKINT, KKKZPW, KNSWOT, KOP, KOTC,
     3     L1, L2, LLLINT, LOTJP, LPAR, MNSWOT, MXNON, NCOMP, NDISK,
     4     NEEDIW, NEEDW, NEQ, NEQIVP, NFC, NFCC, NIC, NOPG,
     5     NPS, NROWA, NROWB, NROWY, NSAFIW, NSAFW, NSWOT, NTAPE, NTP,
     6     NUMORT, NXPTS
      DOUBLE PRECISION A(NROWA,*), AE, ALPHA(*), B(NROWB,*), BETA(*),
     1     C, EPS, FOURU, PWCND, PX, RE, SQOVFL, SRU, TND, TOL, TWOU,
     2     URO, WORK(*), X, XBEG, XEND, XL, XOP, XOT, XPTS(*), XSAV,
     3     Y(NROWY,*), ZQUIT
      CHARACTER*8 XERN1, XERN2
C
C     ******************************************************************
C
      COMMON /DML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMP,NFC
      COMMON /DML18J/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ,
     1                INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,NFCC,
     2                ICOCO
      COMMON /DML15T/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP,
     1                KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
      COMMON /DML17B/ KKKZPW,NEEDW,NEEDIW,K1,K2,K3,K4,K5,K6,K7,K8,K9,
     1                K10,K11,L1,L2,KKKINT,LLLINT
C
      COMMON /DML5MC/ URO,SRU,EPS,SQOVFL,TWOU,FOURU,LPAR
C
C***FIRST EXECUTABLE STATEMENT  DEXBVP
      KOTC = 1
      IEXP = 0
      IF (IWORK(7) .EQ. -1) IEXP = IWORK(8)
C
C     COMPUTE ORTHONORMALIZATION TOLERANCES.
C
   10 TOL = 10.0D0**((-LPAR - IEXP)*2)
C
      IWORK(8) = IEXP
      MXNON = IWORK(2)
C
C **********************************************************************
C **********************************************************************
C
      CALL DBVPOR(Y,NROWY,NCOMP,XPTS,NXPTS,A,NROWA,ALPHA,NIC,B,
     1            NROWB,BETA,NFC,IFLAG,WORK(1),MXNON,WORK(K1),NTP,
     2            IWORK(18),WORK(K2),IWORK(16),WORK(K3),WORK(K4),
     3            WORK(K5),WORK(K6),WORK(K7),WORK(K8),WORK(K9),
     4            WORK(K10),IWORK(L1),NFCC)
C
C **********************************************************************
C **********************************************************************
C     IF DMGSBV RETURNS WITH MESSAGE OF DEPENDENT VECTORS, WE REDUCE
C     ORTHONORMALIZATION TOLERANCE AND TRY AGAIN. THIS IS DONE
C     A MAXIMUM OF 2 TIMES.
C
      IF (IFLAG .NE. 30) GO TO 20
      IF (KOTC .EQ. 3  .OR.  NOPG .EQ. 1) GO TO 30
      KOTC = KOTC + 1
      IEXP = IEXP - 2
      GO TO 10
C
C **********************************************************************
C     IF DBVPOR RETURNS MESSAGE THAT THE MAXIMUM NUMBER OF
C     ORTHONORMALIZATIONS HAS BEEN ATTAINED AND WE CANNOT CONTINUE, THEN
C     WE ESTIMATE THE NEW STORAGE REQUIREMENTS IN ORDER TO SOLVE PROBLEM
C
   20 IF (IFLAG .NE. 13) GO TO 30
      XL = ABS(XEND-XBEG)
      ZQUIT = ABS(X-XBEG)
      INC = 1.5D0*XL/ZQUIT * (MXNON+1)
      IF (NDISK .NE. 1) THEN
         NSAFW = INC*KKKZPW + NEEDW
         NSAFIW = INC*NFCC + NEEDIW
      ELSE
         NSAFW = NEEDW + INC
         NSAFIW = NEEDIW
      ENDIF
C
      WRITE (XERN1, '(I8)') NSAFW
      WRITE (XERN2, '(I8)') NSAFIW
      CALL XERMSG ('SLATEC', 'DEXBVP',
     *   'IN DBVSUP, PREDICTED STORAGE ALLOCATION FOR WORK ARRAY IS ' //
     *   XERN1 // ', PREDICTED STORAGE ALLOCATION FOR IWORK ARRAY IS '
     *   // XERN2, 1, 0)
C
   30 IWORK(1) = MXNON
      RETURN
      END
