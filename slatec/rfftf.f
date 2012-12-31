*DECK RFFTF
      SUBROUTINE RFFTF (N, R, WSAVE)
C***BEGIN PROLOGUE  RFFTF
C***SUBSIDIARY
C***PURPOSE  Compute the forward transform of a real, periodic sequence.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A1
C***TYPE      SINGLE PRECISION (RFFTF-S, CFFTF-C)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C   ********************************************************************
C   *   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   *
C   ********************************************************************
C   *                                                                  *
C   *   This routine uses non-standard Fortran 77 constructs and will  *
C   *   be removed from the library at a future date.  You are         *
C   *   requested to use RFFTF1.                                       *
C   *                                                                  *
C   ********************************************************************
C
C   Subroutine RFFTF computes the Fourier coefficients of a real
C   periodic sequence (Fourier analysis).  The transform is defined
C   below at output parameter R.
C
C   Input Arguments
C
C   N       the length of the array R to be transformed.  The method
C           is most efficient when N is a product of small primes.
C           N may change so long as different work arrays are provided.
C
C   R       a real array of length N which contains the sequence
C           to be transformed.
C
C   WSAVE   a work array which must be dimensioned at least 2*N+15
C           in the program that calls RFFTF.  The WSAVE array must be
C           initialized by calling subroutine RFFTI, and a different
C           WSAVE array must be used for each different value of N.
C           This initialization does not have to be repeated so long as
C           remains unchanged.  Thus subsequent transforms can be
C           obtained faster than the first.  Moreover, the same WSAVE
C           array can be used by RFFTF and RFFTB as long as N remains
C           unchanged.
C
C   Output Argument
C
C   R       R(1) = the sum from I=1 to I=N of R(I)
C
C           If N is even set L = N/2; if N is odd set L = (N+1)/2
C
C             then for K = 2,...,L
C
C                R(2*K-2) = the sum from I = 1 to I = N of
C
C                     R(I)*COS((K-1)*(I-1)*2*PI/N)
C
C                R(2*K-1) = the sum from I = 1 to I = N of
C
C                    -R(I)*SIN((K-1)*(I-1)*2*PI/N)
C
C           If N is even
C
C                R(N) = the sum from I = 1 to I = N of
C
C                     (-1)**(I-1)*R(I)
C
C   Note:  This transform is unnormalized since a call of RFFTF
C          followed by a call of RFFTB will multiply the input
C          sequence by N.
C
C   WSAVE  contains results which must not be destroyed between
C          calls of RFFTF or RFFTB.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  RFFTF1
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
C***END PROLOGUE  RFFTF
      DIMENSION R(*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  RFFTF
      IF (N .EQ. 1) RETURN
      CALL RFFTF1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
