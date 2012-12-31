*DECK RFFTI
      SUBROUTINE RFFTI (N, WSAVE)
C***BEGIN PROLOGUE  RFFTI
C***SUBSIDIARY
C***PURPOSE  Initialize a work array for RFFTF and RFFTB.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A1
C***TYPE      SINGLE PRECISION (RFFTI-S, CFFTI-C)
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
C   *   requested to use RFFTI1.                                       *
C   *                                                                  *
C   ********************************************************************
C
C   Subroutine RFFTI initializes the array WSAVE which is used in
C   both RFFTF and RFFTB.  The prime factorization of N together with
C   a tabulation of the trigonometric functions are computed and
C   stored in WSAVE.
C
C   Input Argument
C
C   N       the length of the sequence to be transformed.
C
C   Output Argument
C
C   WSAVE   a work array which must be dimensioned at least 2*N+15.
C           The same work array can be used for both RFFTF and RFFTB
C           as long as N remains unchanged.  Different WSAVE arrays
C           are required for different values of N.  The contents of
C           WSAVE must not be changed between calls of RFFTF or RFFTB.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  RFFTI1
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
C***END PROLOGUE  RFFTI
      DIMENSION WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  RFFTI
      IF (N .EQ. 1) RETURN
      CALL RFFTI1 (N,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
