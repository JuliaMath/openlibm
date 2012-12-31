*DECK CFFTF
      SUBROUTINE CFFTF (N, C, WSAVE)
C***BEGIN PROLOGUE  CFFTF
C***SUBSIDIARY
C***PURPOSE  Compute the forward transform of a complex, periodic
C            sequence.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A2
C***TYPE      COMPLEX (RFFTF-S, CFFTF-C)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  ********************************************************************
C  *   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   *
C  ********************************************************************
C  *                                                                  *
C  *   This routine uses non-standard Fortran 77 constructs and will  *
C  *   be removed from the library at a future date.  You are         *
C  *   requested to use CFFTF1.                                       *
C  *                                                                  *
C  ********************************************************************
C
C  Subroutine CFFTF computes the forward complex discrete Fourier
C  transform (the Fourier analysis).  Equivalently, CFFTF computes
C  the Fourier coefficients of a complex periodic sequence.
C  The transform is defined below at output parameter C.
C
C  The transform is not normalized.  To obtain a normalized transform
C  the output must be divided by N.  Otherwise a call of CFFTF
C  followed by a call of CFFTB will multiply the sequence by N.
C
C  The array WSAVE which is used by subroutine CFFTF must be
C  initialized by calling subroutine CFFTI(N,WSAVE).
C
C  Input Parameters
C
C  N       the length of the complex sequence C.  The method is
C          more efficient when N is the product of small primes.
C
C  C       a complex array of length N which contains the sequence
C
C  WSAVE   a real work array which must be dimensioned at least 4*N+15
C          in the program that calls CFFTF.  The WSAVE array must be
C          initialized by calling subroutine CFFTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C          The same WSAVE array can be used by CFFTF and CFFTB.
C
C  Output Parameters
C
C  C       For J=1,...,N
C
C              C(J)=the sum from K=1,...,N of
C
C                 C(K)*EXP(-I*(J-1)*(K-1)*2*PI/N)
C
C                         where I=SQRT(-1)
C
C  WSAVE   contains initialization calculations which must not be
C          destroyed between calls of subroutine CFFTF or CFFTB
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  CFFTF1
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*).
C   861211  REVISION DATE from Version 3.2
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900131  Routine changed from user-callable to subsidiary
C           because of non-standard Fortran 77 arguments in the
C           call to CFFTB1.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CFFTF
      COMPLEX C
      DIMENSION C(*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  CFFTF
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
