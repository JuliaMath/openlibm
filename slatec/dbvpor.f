*DECK DBVPOR
      SUBROUTINE DBVPOR (Y, NROWY, NCOMP, XPTS, NXPTS, A, NROWA, ALPHA,
     +   NIC, B, NROWB, BETA, NFC, IFLAG, Z, MXNON, P, NTP, IP, W, NIV,
     +   YHP, U, V, COEF, S, STOWA, G, WORK, IWORK, NFCC)
C***BEGIN PROLOGUE  DBVPOR
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (BVPOR-S, DBVPOR-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C **********************************************************************
C     INPUT to DBVPOR    (items not defined in DBVSUP comments)
C **********************************************************************
C
C     NOPG = 0 -- orthonormalization points not pre-assigned
C          = 1 -- orthonormalization points pre-assigned
C
C     MXNON = maximum number of orthogonalizations allowed.
C
C     NDISK = 0 -- in-core storage
C           = 1 -- disk storage.  Value of NTAPE in data statement
C                  is set to 13.  If another value is desired,
C                  the data statement must be changed.
C
C     INTEG = type of integrator and associated test to be used
C             to determine when to orthonormalize.
C
C             1 -- use GRAM-SCHMIDT test and DDERKF
C             2 -- use GRAM-SCHMIDT test and DDEABM
C
C     TOL = tolerance for allowable error in orthogonalization test.
C
C     NPS = 0 normalize particular solution to unit length at each
C             point of orthonormalization.
C         = 1 do not normalize particular solution.
C
C     NTP = must be .GE. NFC*(NFC+1)/2.
C
C     NFCC = 2*NFC for special treatment of a COMPLEX*16 valued problem
C
C     ICOCO = 0 skip final computations (superposition coefficients
C               and, hence, boundary problem solution)
C           = 1 calculate superposition coefficients and obtain
C               solution to the boundary value problem
C
C **********************************************************************
C     OUTPUT from DBVPOR
C **********************************************************************
C
C     Y(NROWY,NXPTS) = solution at specified output points.
C
C     MXNON = number of orthonormalizations performed by DBVPOR.
C
C     Z(MXNON+1) = locations of orthonormalizations performed by DBVPOR.
C
C     NIV = number of independent vectors returned from DMGSBV. Normally
C           this parameter will be meaningful only when DMGSBV returns
C           with MFLAG = 2.
C
C **********************************************************************
C
C     The following variables are in the argument list because of
C     variable dimensioning.  In general, they contain no information of
C     use to the user.  The amount of storage set aside by the user must
C     be greater than or equal to that indicated by the dimension
C     statements.  For the disk storage mode, NON = 0 and KPTS = 1,
C     while for the in-core storage mode, NON = MXNON and KPTS = NXPTS.
C
C     P(NTP,NON+1)
C     IP(NFCC,NON+1)
C     YHP(NCOMP,NFC+1)  plus an additional column of the length  NEQIVP
C     U(NCOMP,NFC,KPTS)
C     V(NCOMP,KPTS)
C     W(NFCC,NON+1)
C     COEF(NFCC)
C     S(NFC+1)
C     STOWA(NCOMP*(NFC+1)+NEQIVP+1)
C     G(NCOMP)
C     WORK(KKKWS)
C     IWORK(LLLIWS)
C
C **********************************************************************
C     SUBROUTINES used by DBVPOR
C         DLSSUD -- solves an underdetermined system of linear
C                   equations.  This routine is used to get a full
C                   set of initial conditions for integration.
C                   Called by DBVPOR.
C
C         DVECS -- obtains starting vectors for special treatment
C                   of COMPLEX*16 valued problems, called by DBVPOR.
C
C         DRKFAB -- routine which conducts integration using DDERKF or
C                   DDEABM.
C
C         DSTWAY -- storage for backup capability, called by
C                   DBVPOR and DREORT.
C
C         DSTOR1 -- storage at output points, called by DBVPOR,
C                   DRKFAB, DREORT and DSTWAY.
C
C         DDOT -- single precision vector inner product routine,
C                   called by DBVPOR, DCOEF, DLSSUD, DMGSBV,
C                   DBKSOL, DREORT and DPRVEC.
C         ** NOTE **
C         a considerable improvement in speed can be achieved if a
C         machine language version is used for DDOT.
C
C         DCOEF -- computes the superposition constants from the
C                   boundary conditions at XFINAL.
C
C         DBKSOL -- solves an upper triangular set of linear equations.
C
C **********************************************************************
C
C***SEE ALSO  DBVSUP
C***ROUTINES CALLED  DBKSOL, DCOEF, DDOT, DLSSUD, DRKFAB, DSTOR1,
C                    DSTWAY, DVECS
C***COMMON BLOCKS    DML15T, DML18J, DML8SZ
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DBVPOR
C
      DOUBLE PRECISION DDOT
      INTEGER I, I1, I2, IC, ICOCO, IFLAG, IGOFX, INDPVT, INFO, INHOMO,
     1     INTEG, IRA, ISFLG, ISTKOP, IVP, J,
     2     K, KNSWOT, KOD, KOP, KPTS, KWC, KWD, KWS, KWT, L, LOTJP, M,
     3     MNSWOT, MXNON, MXNOND, N, NCOMP, NCOMP2, NCOMPD, NDISK, NDW,
     4     NEQ, NEQIVP, NFC, NFCC, NFCCD, NFCD, NFCP1, NFCP2, NIC,
     5     NICD, NIV, NN, NON, NOPG, NPS, NROWA, NROWB, NROWY, NSWOT,
     6     NTAPE, NTP, NTPD, NUMORT, NXPTS, NXPTSD,
     7     IP(NFCC,*), IWORK(*)
      DOUBLE PRECISION A(NROWA,*), AE, ALPHA(*), B(NROWB,*),
     1     BETA(*), C, COEF(*), G(*), P(NTP,*), PWCND, PX,
     2     RE, S(*), STOWA(*), TND, TOL, U(NCOMP,NFC,*),
     3     V(NCOMP,*), W(NFCC,*), WORK(*), X, XBEG, XEND, XOP,
     4     XOT, XPTS(*), XSAV, Y(NROWY,*), YHP(NCOMP,*),
     5     Z(*)
C
C     ******************************************************************
C
      COMMON /DML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMPD,NFCD
      COMMON /DML15T/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP,
     1                KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
      COMMON /DML18J/ AE,RE,TOL,NXPTSD,NICD,NOPG,MXNOND,NDISK,NTAPE,
     1                NEQ,INDPVT,INTEG,NPS,NTPD,NEQIVP,NUMORT,NFCCD,
     2                ICOCO
C
C      *****************************************************************
C
C***FIRST EXECUTABLE STATEMENT  DBVPOR
      NFCP1 = NFC + 1
      NUMORT = 0
      C = 1.0D0
C
C     ******************************************************************
C         CALCULATE INITIAL CONDITIONS WHICH SATISFY
C                       A*YH(XINITIAL)=0  AND  A*YP(XINITIAL)=ALPHA.
C         WHEN NFC .NE. NFCC DLSSUD DEFINES VALUES YHP IN A MATRIX OF
C         SIZE (NFCC+1)*NCOMP AND ,HENCE, OVERFLOWS THE STORAGE
C         ALLOCATION INTO THE U ARRAY. HOWEVER, THIS IS OKAY SINCE
C         PLENTY OF SPACE IS AVAILABLE IN U AND IT HAS NOT YET BEEN
C         USED.
C
      NDW = NROWA*NCOMP
      KWS = NDW + NIC + 1
      KWD = KWS + NIC
      KWT = KWD + NIC
      KWC = KWT + NIC
      IFLAG = 0
      CALL DLSSUD(A,YHP(1,NFCC+1),ALPHA,NIC,NCOMP,NROWA,YHP,NCOMP,IFLAG,
     1            1,IRA,0,WORK(1),WORK(NDW+1),IWORK,WORK(KWS),WORK(KWD),
     2            WORK(KWT),ISFLG,WORK(KWC))
      IF (IFLAG .EQ. 1) GO TO 10
         IFLAG = -4
      GO TO 200
   10 CONTINUE
         IF (NFC .NE. NFCC)
     1      CALL DVECS(NCOMP,NFC,YHP,WORK,IWORK,INHOMO,IFLAG)
         IF (IFLAG .EQ. 1) GO TO 20
            IFLAG = -5
         GO TO 190
   20    CONTINUE
C
C           ************************************************************
C               DETERMINE THE NUMBER OF DIFFERENTIAL EQUATIONS TO BE
C               INTEGRATED, INITIALIZE VARIABLES FOR AUXILIARY INITIAL
C               VALUE PROBLEM AND STORE INITIAL CONDITIONS.
C
            NEQ = NCOMP*NFC
            IF (INHOMO .EQ. 1) NEQ = NEQ + NCOMP
            IVP = 0
            IF (NEQIVP .EQ. 0) GO TO 40
               IVP = NEQ
               NEQ = NEQ + NEQIVP
               NFCP2 = NFCP1
               IF (INHOMO .EQ. 1) NFCP2 = NFCP1 + 1
               DO 30 K = 1, NEQIVP
                  YHP(K,NFCP2) = ALPHA(NIC+K)
   30          CONTINUE
   40       CONTINUE
            CALL DSTOR1(U,YHP,V,YHP(1,NFCP1),0,NDISK,NTAPE)
C
C           ************************************************************
C               SET UP DATA FOR THE ORTHONORMALIZATION TESTING PROCEDURE
C               AND SAVE INITIAL CONDITIONS IN CASE A RESTART IS
C               NECESSARY.
C
            NSWOT = 1
            KNSWOT = 0
            LOTJP = 1
            TND = LOG10(10.0D0*TOL)
            PWCND = LOG10(SQRT(TOL))
            X = XBEG
            PX = X
            XOT = XEND
            XOP = X
            KOP = 1
            CALL DSTWAY(U,V,YHP,0,STOWA)
C
C           ************************************************************
C           ******** FORWARD INTEGRATION OF ALL INITIAL VALUE EQUATIONS
C           **********
C           ************************************************************
C
            CALL DRKFAB(NCOMP,XPTS,NXPTS,NFC,IFLAG,Z,MXNON,P,NTP,IP,YHP,
     1                  NIV,U,V,W,S,STOWA,G,WORK,IWORK,NFCC)
            IF (IFLAG .NE. 0 .OR. ICOCO .EQ. 0) GO TO 180
C
C              *********************************************************
C              **************** BACKWARD SWEEP TO OBTAIN SOLUTION
C              *******************
C              *********************************************************
C
C                  CALCULATE SUPERPOSITION COEFFICIENTS AT XFINAL.
C
C                FOR THE DISK STORAGE VERSION, IT IS NOT NECESSARY TO
C                READ  U  AND  V AT THE LAST OUTPUT POINT, SINCE THE
C                LOCAL COPY OF EACH STILL EXISTS.
C
               KOD = 1
               IF (NDISK .EQ. 0) KOD = NXPTS
               I1 = 1 + NFCC*NFCC
               I2 = I1 + NFCC
               CALL DCOEF(U(1,1,KOD),V(1,KOD),NCOMP,NROWB,NFC,NIC,B,
     1                     BETA,COEF,INHOMO,RE,AE,WORK,WORK(I1),
     2                     WORK(I2),IWORK,IFLAG,NFCC)
C
C              *********************************************************
C                  CALCULATE SOLUTION AT OUTPUT POINTS BY RECURRING
C                  BACKWARDS.  AS WE RECUR BACKWARDS FROM XFINAL TO
C                  XINITIAL WE MUST CALCULATE NEW SUPERPOSITION
C                  COEFFICIENTS EACH TIME WE CROSS A POINT OF
C                  ORTHONORMALIZATION.
C
               K = NUMORT
               NCOMP2 = NCOMP/2
               IC = 1
               IF (NFC .NE. NFCC) IC = 2
               DO 170 J = 1, NXPTS
                  KPTS = NXPTS - J + 1
                  KOD = KPTS
                  IF (NDISK .EQ. 1) KOD = 1
   50             CONTINUE
C                 ...EXIT
                     IF (K .EQ. 0) GO TO 120
C                 ...EXIT
                     IF (XEND .GT. XBEG .AND. XPTS(KPTS) .GE. Z(K))
     1                  GO TO 120
C                 ...EXIT
                     IF (XEND .LT. XBEG .AND. XPTS(KPTS) .LE. Z(K))
     1                  GO TO 120
                     NON = K
                     IF (NDISK .EQ. 0) GO TO 60
                        NON = 1
                        BACKSPACE NTAPE
                        READ (NTAPE)
     1                       (IP(I,1), I = 1, NFCC),(P(I,1), I = 1, NTP)
                        BACKSPACE NTAPE
   60                CONTINUE
                     IF (INHOMO .NE. 1) GO TO 90
                        IF (NDISK .EQ. 0) GO TO 70
                           BACKSPACE NTAPE
                           READ (NTAPE) (W(I,1), I = 1, NFCC)
                           BACKSPACE NTAPE
   70                   CONTINUE
                        DO 80 N = 1, NFCC
                           COEF(N) = COEF(N) - W(N,NON)
   80                   CONTINUE
   90                CONTINUE
                     CALL DBKSOL(NFCC,P(1,NON),COEF)
                     DO 100 M = 1, NFCC
                        WORK(M) = COEF(M)
  100                CONTINUE
                     DO 110 M = 1, NFCC
                        L = IP(M,NON)
                        COEF(L) = WORK(M)
  110                CONTINUE
                     K = K - 1
                  GO TO 50
  120             CONTINUE
                  IF (NDISK .EQ. 0) GO TO 130
                     BACKSPACE NTAPE
                     READ (NTAPE)
     1                    (V(I,1), I = 1, NCOMP),
     2                    ((U(I,M,1), I = 1, NCOMP), M = 1, NFC)
                     BACKSPACE NTAPE
  130             CONTINUE
                  DO 140 N = 1, NCOMP
                     Y(N,KPTS) = V(N,KOD)
     1                           + DDOT(NFC,U(N,1,KOD),NCOMP,COEF,IC)
  140             CONTINUE
                  IF (NFC .EQ. NFCC) GO TO 160
                     DO 150 N = 1, NCOMP2
                        NN = NCOMP2 + N
                        Y(N,KPTS) = Y(N,KPTS)
     1                              - DDOT(NFC,U(NN,1,KOD),NCOMP,
     2                                     COEF(2),2)
                        Y(NN,KPTS) = Y(NN,KPTS)
     1                               + DDOT(NFC,U(N,1,KOD),NCOMP,
     2                                      COEF(2),2)
  150                CONTINUE
  160             CONTINUE
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
C
C     ******************************************************************
C
      MXNON = NUMORT
      RETURN
      END
