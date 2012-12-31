*DECK XPSI
      REAL FUNCTION XPSI (A, IPSIK, IPSIX)
C***BEGIN PROLOGUE  XPSI
C***SUBSIDIARY
C***PURPOSE  To compute values of the Psi function for XLEGF.
C***LIBRARY   SLATEC
C***CATEGORY  C7C
C***TYPE      SINGLE PRECISION (XPSI-S, DXPSI-D)
C***KEYWORDS  PSI FUNCTION
C***AUTHOR  Smith, John M., (NBS and George Mason University)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820728  DATE WRITTEN
C   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Changed all specific intrinsics to generic.  (WRB)
C           Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  XPSI
      REAL A,B,C,CNUM,CDENOM
      DIMENSION CNUM(12),CDENOM(12)
      SAVE CNUM, CDENOM
C
C        CNUM(I) AND CDENOM(I) ARE THE ( REDUCED ) NUMERATOR
C        AND 2*I*DENOMINATOR RESPECTIVELY OF THE 2*I TH BERNOULLI
C        NUMBER.
C
      DATA CNUM(1),CNUM(2),CNUM(3),CNUM(4),CNUM(5),CNUM(6),CNUM(7),
     1CNUM(8),CNUM(9),CNUM(10),CNUM(11),CNUM(12)
     2    / 1.,     -1.,    1.,     -1., 1.,
     3   -691.,  1.,     -3617., 43867., -174611., 77683.,
     4   -236364091./
      DATA CDENOM(1),CDENOM(2),CDENOM(3),CDENOM(4),CDENOM(5),CDENOM(6),
     1 CDENOM(7),CDENOM(8),CDENOM(9),CDENOM(10),CDENOM(11),CDENOM(12)
     2/12.,120.,   252.,   240.,132.,
     3  32760., 12.,  8160., 14364., 6600., 276., 65520./
C***FIRST EXECUTABLE STATEMENT  XPSI
      N=MAX(0,IPSIX-INT(A))
      B=N+A
      K1=IPSIK-1
C
C        SERIES EXPANSION FOR A .GT. IPSIX USING IPSIK-1 TERMS.
C
      C=0.
      DO 12 I=1,K1
      K=IPSIK-I
   12 C=(C+CNUM(K)/CDENOM(K))/B**2
      XPSI=LOG(B)-(C+.5/B)
      IF(N.EQ.0) GO TO 20
      B=0.
C
C        RECURRENCE FOR A .LE. IPSIX.
C
      DO 15 M=1,N
   15 B=B+1./(N-M+A)
      XPSI=XPSI-B
   20 RETURN
      END
