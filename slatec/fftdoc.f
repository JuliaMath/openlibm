*DECK FFTDOC
      SUBROUTINE FFTDOC
C***BEGIN PROLOGUE  FFTDOC
C***PURPOSE  Documentation for FFTPACK, a collection of Fast Fourier
C            Transform routines.
C***LIBRARY   SLATEC
C***CATEGORY  J1, Z
C***TYPE      ALL (FFTDOC-A)
C***KEYWORDS  DOCUMENTATION, FAST FOURIER TRANSFORM, FFT
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                       Version 3  June 1979
C
C          A Package of Fortran Subprograms for The Fast Fourier
C           Transform of Periodic and Other Symmetric Sequences
C                              By
C                       Paul N Swarztrauber
C
C    National Center For Atmospheric Research, Boulder, Colorado 80307
C        which is sponsored by the National Science Foundation
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     This package consists of programs which perform Fast Fourier
C     Transforms for both complex and real periodic sequences and
C     certain other symmetric sequences that are listed below.
C
C     1.   RFFTI     Initialize RFFTF and RFFTB
C     2.   RFFTF     Forward transform of a real periodic sequence
C     3.   RFFTB     Backward transform of a real coefficient array
C
C     4.   EZFFTI    Initialize EZFFTF and EZFFTB
C     5.   EZFFTF    A simplified real periodic forward transform
C     6.   EZFFTB    A simplified real periodic backward transform
C
C     7.   SINTI     Initialize SINT
C     8.   SINT      Sine transform of a real odd sequence
C
C     9.   COSTI     Initialize COST
C     10.  COST      Cosine transform of a real even sequence
C
C     11.  SINQI     Initialize SINQF and SINQB
C     12.  SINQF     Forward sine transform with odd wave numbers
C     13.  SINQB     Unnormalized inverse of SINQF
C
C     14.  COSQI     Initialize COSQF and COSQB
C     15.  COSQF     Forward cosine transform with odd wave numbers
C     16.  COSQB     Unnormalized inverse of COSQF
C
C     17.  CFFTI     Initialize CFFTF and CFFTB
C     18.  CFFTF     Forward transform of a complex periodic sequence
C     19.  CFFTB     Unnormalized inverse of CFFTF
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   780201  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900723  PURPOSE section revised.  (WRB)
C***END PROLOGUE  FFTDOC
C***FIRST EXECUTABLE STATEMENT  FFTDOC
       RETURN
      END
