*DECK CFFTB1
      SUBROUTINE CFFTB1 (N, C, CH, WA, IFAC)
C***BEGIN PROLOGUE  CFFTB1
C***PURPOSE  Compute the unnormalized inverse of CFFTF1.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A2
C***TYPE      COMPLEX (RFFTB1-S, CFFTB1-C)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine CFFTB1 computes the backward complex discrete Fourier
C  transform (the Fourier synthesis).  Equivalently, CFFTB1 computes
C  a complex periodic sequence from its Fourier coefficients.
C  The transform is defined below at output parameter C.
C
C  A call of CFFTF1 followed by a call of CFFTB1 will multiply the
C  sequence by N.
C
C  The arrays WA and IFAC which are used by subroutine CFFTB1 must be
C  initialized by calling subroutine CFFTI1 (N, WA, IFAC).
C
C  Input Parameters
C
C  N       the length of the complex sequence C.  The method is
C          more efficient when N is the product of small primes.
C
C  C       a complex array of length N which contains the sequence
C
C  CH      a real work array of length at least 2*N
C
C  WA      a real work array which must be dimensioned at least 2*N.
C
C  IFAC    an integer work array which must be dimensioned at least 15.
C
C          The WA and IFAC arrays must be initialized by calling
C          subroutine CFFTI1 (N, WA, IFAC), and different WA and IFAC
C          arrays must be used for each different value of N.  This
C          initialization does not have to be repeated so long as N
C          remains unchanged.  Thus subsequent transforms can be
C          obtained faster than the first.  The same WA and IFAC arrays
C          can be used by CFFTF1 and CFFTB1.
C
C  Output Parameters
C
C  C       For J=1,...,N
C
C              C(J)=the sum from K=1,...,N of
C
C                 C(K)*EXP(I*(J-1)*(K-1)*2*PI/N)
C
C                         where I=SQRT(-1)
C
C  NOTE:   WA and IFAC contain initialization calculations which must
C          not be destroyed between calls of subroutine CFFTF1 or CFFTB1
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  PASSB, PASSB2, PASSB3, PASSB4, PASSB5
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*).
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900131  Routine changed from subsidiary to user-callable.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CFFTB1
      DIMENSION CH(*), C(*), WA(*), IFAC(*)
C***FIRST EXECUTABLE STATEMENT  CFFTB1
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL PASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL PASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL PASSB2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL PASSB2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL PASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL PASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL PASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL PASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL PASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL PASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
