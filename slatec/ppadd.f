*DECK PPADD
      SUBROUTINE PPADD (N, IERROR, A, C, CBP, BP, BH)
C***BEGIN PROLOGUE  PPADD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BLKTRI
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (PPADD-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   PPADD computes the eigenvalues of the periodic tridiagonal matrix
C   with coefficients AN,BN,CN.
C
C   N    is the order of the BH and BP polynomials.
C   BP   contains the eigenvalues on output.
C   CBP  is the same as BP except type complex.
C   BH   is used to temporarily store the roots of the B HAT polynomial
C        which enters through BP.
C
C***SEE ALSO  BLKTRI
C***ROUTINES CALLED  BSRH, PPSGF, PPSPF, PSGF
C***COMMON BLOCKS    CBLKT
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  PPADD
C
      COMPLEX         CX         ,FSG        ,HSG        ,
     1                DD         ,F          ,FP         ,FPP        ,
     2                CDIS       ,R1         ,R2         ,R3         ,
     3                CBP
      DIMENSION       A(*)       ,C(*)       ,BP(*)      ,BH(*)      ,
     1                CBP(*)
      COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
      EXTERNAL        PSGF       ,PPSPF      ,PPSGF
C***FIRST EXECUTABLE STATEMENT  PPADD
      SCNV = SQRT(CNV)
      IZ = N
      IF (BP(N)-BP(1)) 101,142,103
  101 DO 102 J=1,N
         NT = N-J
         BH(J) = BP(NT+1)
  102 CONTINUE
      GO TO 105
  103 DO 104 J=1,N
         BH(J) = BP(J)
  104 CONTINUE
  105 NCMPLX = 0
      MODIZ = MOD(IZ,2)
      IS = 1
      IF (MODIZ) 106,107,106
  106 IF (A(1)) 110,142,107
  107 XL = BH(1)
      DB = BH(3)-BH(1)
  108 XL = XL-DB
      IF (PSGF(XL,IZ,C,A,BH)) 108,108,109
  109 SGN = -1.
      CBP(1) = CMPLX(BSRH(XL,BH(1),IZ,C,A,BH,PSGF,SGN),0.)
      IS = 2
  110 IF = IZ-1
      IF (MODIZ) 111,112,111
  111 IF (A(1)) 112,142,115
  112 XR = BH(IZ)
      DB = BH(IZ)-BH(IZ-2)
  113 XR = XR+DB
      IF (PSGF(XR,IZ,C,A,BH)) 113,114,114
  114 SGN = 1.
      CBP(IZ) = CMPLX(BSRH(BH(IZ),XR,IZ,C,A,BH,PSGF,SGN),0.)
      IF = IZ-2
  115 DO 136 IG=IS,IF,2
         XL = BH(IG)
         XR = BH(IG+1)
         SGN = -1.
         XM = BSRH(XL,XR,IZ,C,A,BH,PPSPF,SGN)
         PSG = PSGF(XM,IZ,C,A,BH)
         IF (ABS(PSG)-EPS) 118,118,116
  116    IF (PSG*PPSGF(XM,IZ,C,A,BH)) 117,118,119
C
C     CASE OF A REAL ZERO
C
  117    SGN = 1.
         CBP(IG) = CMPLX(BSRH(BH(IG),XM,IZ,C,A,BH,PSGF,SGN),0.)
         SGN = -1.
         CBP(IG+1) = CMPLX(BSRH(XM,BH(IG+1),IZ,C,A,BH,PSGF,SGN),0.)
         GO TO 136
C
C     CASE OF A MULTIPLE ZERO
C
  118    CBP(IG) = CMPLX(XM,0.)
         CBP(IG+1) = CMPLX(XM,0.)
         GO TO 136
C
C     CASE OF A COMPLEX ZERO
C
  119    IT = 0
         ICV = 0
         CX = CMPLX(XM,0.)
  120    FSG = (1.,0.)
         HSG = (1.,0.)
         FP = (0.,0.)
         FPP = (0.,0.)
         DO 121 J=1,IZ
            DD = 1./(CX-BH(J))
            FSG = FSG*A(J)*DD
            HSG = HSG*C(J)*DD
            FP = FP+DD
            FPP = FPP-DD*DD
  121    CONTINUE
         IF (MODIZ) 123,122,123
  122    F = (1.,0.)-FSG-HSG
         GO TO 124
  123    F = (1.,0.)+FSG+HSG
  124    I3 = 0
         IF (ABS(FP)) 126,126,125
  125    I3 = 1
         R3 = -F/FP
  126    IF (ABS(FPP)) 132,132,127
  127    CDIS = SQRT(FP**2-2.*F*FPP)
         R1 = CDIS-FP
         R2 = -FP-CDIS
         IF (ABS(R1)-ABS(R2)) 129,129,128
  128    R1 = R1/FPP
         GO TO 130
  129    R1 = R2/FPP
  130    R2 = 2.*F/FPP/R1
         IF (ABS(R2) .LT. ABS(R1)) R1 = R2
         IF (I3) 133,133,131
  131    IF (ABS(R3) .LT. ABS(R1)) R1 = R3
         GO TO 133
  132    R1 = R3
  133    CX = CX+R1
         IT = IT+1
         IF (IT .GT. 50) GO TO 142
         IF (ABS(R1) .GT. SCNV) GO TO 120
         IF (ICV) 134,134,135
  134    ICV = 1
         GO TO 120
  135    CBP(IG) = CX
         CBP(IG+1) = CONJG(CX)
  136 CONTINUE
      IF (ABS(CBP(N))-ABS(CBP(1))) 137,142,139
  137 NHALF = N/2
      DO 138 J=1,NHALF
         NT = N-J
         CX = CBP(J)
         CBP(J) = CBP(NT+1)
         CBP(NT+1) = CX
  138 CONTINUE
  139 NCMPLX = 1
      DO 140 J=2,IZ
         IF (AIMAG(CBP(J))) 143,140,143
  140 CONTINUE
      NCMPLX = 0
      DO 141 J=2,IZ
         BP(J) = REAL(CBP(J))
  141 CONTINUE
      GO TO 143
  142 IERROR = 4
  143 CONTINUE
      RETURN
      END
