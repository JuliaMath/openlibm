*DECK BVPOR
      SUBROUTINE BVPOR (Y, NROWY, NCOMP, XPTS, NXPTS, A, NROWA, ALPHA,
     +   NIC, B, NROWB, BETA, NFC, IFLAG, Z, MXNON, P, NTP, IP, W, NIV,
     +   YHP, U, V, COEF, S, STOWA, G, WORK, IWORK, NFCC)
C***BEGIN PROLOGUE  BVPOR
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BVPOR-S, DBVPOR-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C **********************************************************************
C     INPUT to BVPOR    (items not defined in BVSUP comments)
C **********************************************************************
C
C     NOPG = 0 -- Orthonormalization points not pre-assigned
C          = 1 -- Orthonormalization points pre-assigned
C
C     MXNON = Maximum number of orthogonalizations allowed.
C
C     NDISK = 0 -- IN-CORE storage
C           = 1 -- DISK storage.  Value of NTAPE in data statement
C                  is set to 13.  If another value is desired,
C                  the data statement must be changed.
C
C     INTEG = Type of integrator and associated test to be used
C             to determine when to orthonormalize.
C
C             1 -- Use GRAM-SCHMIDT test and DERKF
C             2 -- Use GRAM-SCHMIDT test and DEABM
C
C     TOL = Tolerance for allowable error in orthogonalization test.
C
C     NPS = 0 Normalize particular solution to unit length at each
C             point of orthonormalization.
C         = 1 Do not normalize particular solution.
C
C     NTP = Must be .GE. NFC*(NFC+1)/2.
C
C
C     NFCC = 2*NFC for special treatment of a complex valued problem
C
C     ICOCO = 0 Skip final computations (superposition coefficients
C               and ,hence, boundary problem solution)
C           = 1 Calculate superposition coefficients and obtain
C               solution to the boundary value problem
C
C **********************************************************************
C     OUTPUT from BVPOR
C **********************************************************************
C
C     Y(NROWY,NXPTS) = Solution at specified output points.
C
C     MXNON = Number of orthonormalizations performed by BVPOR.
C
C     Z(MXNON+1) = Locations of orthonormalizations performed by BVPOR.
C
C     NIV = Number of independent vectors returned from MGSBV. Normally
C        this parameter will be meaningful only when MGSBV returns with
C           MFLAG = 2.
C
C **********************************************************************
C
C     The following variables are in the argument list because of
C     variable dimensioning. In general, they contain no information of
C     use to the user.  The amount of storage set aside by the user must
C     be greater than or equal to that indicated by the dimension
C     statements.   For the DISK storage mode, NON = 0 and KPTS = 1,
C     while for the IN-CORE storage mode, NON = MXNON and KPTS = NXPTS.
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
C     Subroutines used by BVPOR
C         LSSUDS -- Solves an underdetermined system of linear
C                   equations.  This routine is used to get a full
C                   set of initial conditions for integration.
C                   Called by BVPOR
C
C         SVECS -- Obtains starting vectors for special treatment
C                  of complex valued problems , called by BVPOR
C
C         RKFAB -- Routine which conducts integration using DERKF or
C                   DEABM
C
C         STWAY -- Storage for backup capability, called by
C                   BVPOR and REORT
C
C         STOR1 -- Storage at output points, called by BVPOR,
C                  RKFAB, REORT and STWAY.
C
C         SDOT -- Single precision vector inner product routine,
C                   called by BVPOR, SCOEF, LSSUDS, MGSBV,
C                   BKSOL, REORT and PRVEC.
C         ** NOTE **
C         A considerable improvement in speed can be achieved if a
C         machine language version is used for SDOT.
C
C         SCOEF -- Computes the superposition constants from the
C                  boundary conditions at Xfinal.
C
C         BKSOL -- Solves an upper triangular set of linear equations.
C
C **********************************************************************
C
C***SEE ALSO  BVSUP
C***ROUTINES CALLED  BKSOL, LSSUDS, RKFAB, SCOEF, SDOT, STOR1, STWAY,
C                    SVECS
C***COMMON BLOCKS    ML15TO, ML18JR, ML8SZ
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  BVPOR
C
      DIMENSION Y(NROWY,*),A(NROWA,*),ALPHA(*),B(NROWB,*),
     1          BETA(*),P(NTP,*),IP(NFCC,*),
     2          U(NCOMP,NFC,*),V(NCOMP,*),W(NFCC,*),
     3          COEF(*),Z(*),YHP(NCOMP,*),XPTS(*),S(*),
     4          WORK(*),IWORK(*),STOWA(*),G(*)
C
C **********************************************************************
C
      COMMON /ML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMPD,NFCD
      COMMON /ML15TO/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP,
     1                KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
      COMMON /ML18JR/ AE,RE,TOL,NXPTSD,NICD,NOPG,MXNOND,NDISK,NTAPE,
     1                NEQ,INDPVT,INTEG,NPS,NTPD,NEQIVP,NUMORT,NFCCD,
     2                ICOCO
C
C **********************************************************************
C
C***FIRST EXECUTABLE STATEMENT  BVPOR
      NFCP1 = NFC + 1
      NUMORT = 0
      C = 1.0
C
C **********************************************************************
C     CALCULATE INITIAL CONDITIONS WHICH SATISFY
C                   A*YH(XINITIAL)=0  AND  A*YP(XINITIAL)=ALPHA.
C     WHEN NFC .NE. NFCC LSSUDS DEFINES VALUES YHP IN A MATRIX OF SIZE
C     (NFCC+1)*NCOMP AND ,HENCE, OVERFLOWS THE STORAGE ALLOCATION INTO
C     THE U ARRAY. HOWEVER, THIS IS OKAY SINCE PLENTY OF SPACE IS
C     AVAILABLE IN U AND IT HAS NOT YET BEEN USED.
C
      NDW = NROWA * NCOMP
      KWS = NDW + NIC + 1
      KWD = KWS + NIC
      KWT = KWD + NIC
      KWC = KWT + NIC
      IFLAG = 0
      CALL LSSUDS(A,YHP(1,NFCC+1),ALPHA,NIC,NCOMP,NROWA,YHP,NCOMP,
     1            IFLAG,1,IRA,0,WORK(1),WORK(NDW+1),IWORK,WORK(KWS),
     2            WORK(KWD),WORK(KWT),ISFLG,WORK(KWC))
      IF (IFLAG .EQ. 1) GO TO 3
      IFLAG=-4
      GO TO 250
    3 IF (NFC .NE. NFCC) CALL SVECS(NCOMP,NFC,YHP,WORK,IWORK,
     1                   INHOMO,IFLAG)
      IF (IFLAG .EQ. 1)  GO TO 5
      IFLAG=-5
      GO TO 250
C
C **********************************************************************
C     DETERMINE THE NUMBER OF DIFFERENTIAL EQUATIONS TO BE INTEGRATED,
C     INITIALIZE VARIABLES FOR AUXILIARY INITIAL VALUE PROBLEM AND
C     STORE INITIAL CONDITIONS.
C
    5 NEQ = NCOMP * NFC
      IF (INHOMO .EQ. 1)  NEQ = NEQ + NCOMP
      IVP = 0
      IF (NEQIVP .EQ. 0)  GO TO 10
      IVP = NEQ
      NEQ = NEQ + NEQIVP
      NFCP2 = NFCP1
      IF (INHOMO .EQ. 1)  NFCP2 = NFCP1 + 1
      DO 7 K = 1,NEQIVP
    7 YHP(K,NFCP2) = ALPHA(NIC+K)
   10 CALL STOR1(U,YHP,V,YHP(1,NFCP1),0,NDISK,NTAPE)
C
C **********************************************************************
C     SET UP DATA FOR THE ORTHONORMALIZATION TESTING PROCEDURE AND
C     SAVE INITIAL CONDITIONS IN CASE A RESTART IS NECESSARY.
C
      NSWOT=1
      KNSWOT=0
      LOTJP=1
      TND=LOG10(10.*TOL)
      PWCND=LOG10(SQRT(TOL))
      X=XBEG
      PX=X
      XOT=XEND
      XOP=X
      KOP=1
      CALL STWAY(U,V,YHP,0,STOWA)
C
C **********************************************************************
C ******** FORWARD INTEGRATION OF ALL INITIAL VALUE EQUATIONS **********
C **********************************************************************
C
      CALL RKFAB(NCOMP,XPTS,NXPTS,NFC,IFLAG,Z,MXNON,P,NTP,IP,
     1            YHP,NIV,U,V,W,S,STOWA,G,WORK,IWORK,NFCC)
      IF (IFLAG .NE. 0  .OR.  ICOCO .EQ. 0)  GO TO 250
C
C **********************************************************************
C **************** BACKWARD SWEEP TO OBTAIN SOLUTION *******************
C **********************************************************************
C
C     CALCULATE SUPERPOSITION COEFFICIENTS AT XFINAL.
C
C   FOR THE DISK STORAGE VERSION, IT IS NOT NECESSARY TO READ  U  AND  V
C   AT THE LAST OUTPUT POINT, SINCE THE LOCAL COPY OF EACH STILL EXISTS.
C
      KOD = 1
      IF (NDISK .EQ. 0)  KOD = NXPTS
      I1=1+NFCC*NFCC
      I2=I1+NFCC
      CALL SCOEF(U(1,1,KOD),V(1,KOD),NCOMP,NROWB,NFC,NIC,B,BETA,COEF,
     1           INHOMO,RE,AE,WORK,WORK(I1),WORK(I2),IWORK,IFLAG,NFCC)
C
C **********************************************************************
C     CALCULATE SOLUTION AT OUTPUT POINTS BY RECURRING BACKWARDS.
C     AS WE RECUR BACKWARDS FROM XFINAL TO XINITIAL WE MUST CALCULATE
C     NEW SUPERPOSITION COEFFICIENTS EACH TIME WE CROSS A POINT OF
C     ORTHONORMALIZATION.
C
      K = NUMORT
      NCOMP2=NCOMP/2
      IC=1
      IF (NFC .NE. NFCC) IC=2
      DO 200 J = 1,NXPTS
      KPTS = NXPTS - J + 1
      KOD = KPTS
      IF (NDISK .EQ. 1)  KOD = 1
  135 IF (K .EQ. 0)  GO TO 170
      IF (XEND.GT.XBEG .AND. XPTS(KPTS).GE.Z(K))  GO TO 170
      IF (XEND.LT.XBEG .AND. XPTS(KPTS).LE.Z(K))  GO TO 170
      NON = K
      IF (NDISK .EQ. 0)  GO TO 136
      NON = 1
      BACKSPACE NTAPE
      READ (NTAPE) (IP(I,1), I = 1,NFCC),(P(I,1), I = 1,NTP)
      BACKSPACE NTAPE
  136 IF (INHOMO .NE. 1)  GO TO 150
      IF (NDISK .EQ. 0)  GO TO 138
      BACKSPACE NTAPE
      READ (NTAPE) (W(I,1), I = 1,NFCC)
      BACKSPACE NTAPE
  138 DO 140 N = 1,NFCC
  140 COEF(N) = COEF(N) - W(N,NON)
  150 CALL BKSOL(NFCC,P(1,NON),COEF)
      DO 155 M = 1,NFCC
  155 WORK(M) = COEF(M)
      DO 160 M = 1,NFCC
      L = IP(M,NON)
  160 COEF(L) = WORK(M)
      K = K - 1
      GO TO 135
  170 IF (NDISK .EQ. 0)  GO TO 175
      BACKSPACE NTAPE
      READ (NTAPE) (V(I,1), I = 1,NCOMP),
     1             ((U(I,M,1), I = 1,NCOMP), M = 1,NFC)
      BACKSPACE NTAPE
  175 DO 180 N = 1,NCOMP
  180 Y(N,KPTS) = V(N,KOD) + SDOT(NFC,U(N,1,KOD),NCOMP,COEF,IC)
      IF (NFC .EQ. NFCC) GO TO 200
      DO 190 N=1,NCOMP2
      NN=NCOMP2+N
      Y(N,KPTS)=Y(N,KPTS) - SDOT(NFC,U(NN,1,KOD),NCOMP,COEF(2),2)
  190 Y(NN,KPTS)=Y(NN,KPTS) + SDOT(NFC,U(N,1,KOD),NCOMP,COEF(2),2)
  200 CONTINUE
C
C **********************************************************************
C
  250 MXNON = NUMORT
      RETURN
      END
