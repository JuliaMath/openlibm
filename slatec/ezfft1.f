*DECK EZFFT1
      SUBROUTINE EZFFT1 (N, WA, IFAC)
C***BEGIN PROLOGUE  EZFFT1
C***SUBSIDIARY
C***PURPOSE  EZFFTI calls EZFFT1 with appropriate work array
C            partitioning.
C***LIBRARY   SLATEC (FFTPACK)
C***TYPE      SINGLE PRECISION (EZFFT1-S)
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           (a) changing dummy array size declarations (1) to (*),
C           (b) changing references to intrinsic function FLOAT
C               to REAL, and
C           (c) changing definition of variable TPI by using
C               FORTRAN intrinsic function ATAN instead of a DATA
C               statement.
C   881128  Modified by Dick Valent to meet prologue standards.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  EZFFT1
      DIMENSION WA(*), IFAC(*), NTRYH(4)
      SAVE NTRYH
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
C***FIRST EXECUTABLE STATEMENT  EZFFT1
      TPI = 8.*ATAN(1.)
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      ARGH = TPI/N
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 111 K1=1,NFM1
         IP = IFAC(K1+2)
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         ARG1 = L1*ARGH
         CH1 = 1.
         SH1 = 0.
         DCH1 = COS(ARG1)
         DSH1 = SIN(ARG1)
         DO 110 J=1,IPM
            CH1H = DCH1*CH1-DSH1*SH1
            SH1 = DCH1*SH1+DSH1*CH1
            CH1 = CH1H
            I = IS+2
            WA(I-1) = CH1
            WA(I) = SH1
            IF (IDO .LT. 5) GO TO 109
            DO 108 II=5,IDO,2
               I = I+2
               WA(I-1) = CH1*WA(I-3)-SH1*WA(I-2)
               WA(I) = CH1*WA(I-2)+SH1*WA(I-3)
  108       CONTINUE
  109       IS = IS+IDO
  110    CONTINUE
         L1 = L2
  111 CONTINUE
      RETURN
      END
