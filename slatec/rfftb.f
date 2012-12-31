*DECK RFFTB
      SUBROUTINE RFFTB (N, R, WSAVE)
C***BEGIN PROLOGUE  RFFTB
C***SUBSIDIARY
C***PURPOSE  Compute the backward fast Fourier transform of a real
C            coefficient array.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A1
C***TYPE      SINGLE PRECISION (RFFTB-S, CFFTB-C)
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
C   *   requested to use RFFTB1.                                       *
C   *                                                                  *
C   ********************************************************************
C
C   Subroutine RFFTB computes the real periodic sequence from its
C   Fourier coefficients (Fourier synthesis).  The transform is defined
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
C           in the program that calls RFFTB.  The WSAVE array must be
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
C   R       For N even and for I = 1,...,N
C
C                R(I) = R(1)+(-1)**(I-1)*R(N)
C
C                     plus the sum from K=2 to K=N/2 of
C
C                      2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
C
C                     -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
C
C           For N odd and for I = 1,...,N
C
C                R(I) = R(1) plus the sum from K=2 to K=(N+1)/2 of
C
C                     2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
C
C                    -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
C
C   Note:  This transform is unnormalized since a call of RFFTF
C          followed by a call of RFFTB will multiply the input
C          sequence by N.
C
C   WSAVE  contains results which must not be destroyed between
C          calls of RFFTB or RFFTF.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  RFFTB1
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
C***END PROLOGUE  RFFTB
      DIMENSION R(*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  RFFTB
      IF (N .EQ. 1) RETURN
      CALL RFFTB1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
