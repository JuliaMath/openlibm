*DECK CFFTI
      SUBROUTINE CFFTI (N, WSAVE)
C***BEGIN PROLOGUE  CFFTI
C***SUBSIDIARY
C***PURPOSE  Initialize a work array for CFFTF and CFFTB.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A2
C***TYPE      COMPLEX (RFFTI-S, CFFTI-C)
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
C  *   requested to use CFFTI1.                                       *
C  *                                                                  *
C  ********************************************************************
C
C  Subroutine CFFTI initializes the array WSAVE which is used in
C  both CFFTF and CFFTB.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least 4*N+15.
C          The same work array can be used for both CFFTF and CFFTB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C          WSAVE must not be changed between calls of CFFTF or CFFTB.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  CFFTI1
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
C***END PROLOGUE  CFFTI
      DIMENSION WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  CFFTI
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTI1 (N,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
