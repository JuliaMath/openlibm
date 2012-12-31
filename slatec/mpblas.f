*DECK MPBLAS
      SUBROUTINE MPBLAS (I1)
C***BEGIN PROLOGUE  MPBLAS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPBLAS-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine is called to set up Brent's 'mp' package
C     for use by the extended precision inner products from the BLAS.
C
C     In the SLATEC library we require the Extended Precision MP number
C     to have a mantissa twice as long as Double Precision numbers.
C     The calculation of MPT (and MPMXR which is the actual array size)
C     in this routine will give 2x (or slightly more) on the machine
C     that we are running on.  The INTEGER array size of 30 was chosen
C     to be slightly longer than the longest INTEGER array needed on
C     any machine that we are currently aware of.
C
C***SEE ALSO  DQDOTA, DQDOTI
C***REFERENCES  R. P. Brent, A Fortran multiple-precision arithmetic
C                 package, ACM Transactions on Mathematical Software 4,
C                 1 (March 1978), pp. 57-70.
C               R. P. Brent, MP, a Fortran multiple-precision arithmetic
C                 package, Algorithm 524, ACM Transactions on Mathema-
C                 tical Software 4, 1 (March 1978), pp. 71-81.
C***ROUTINES CALLED  I1MACH, XERMSG
C***COMMON BLOCKS    MPCOM
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   930124  Increased Array size in MPCON for SUN -r8, and calculate
C               size for Quad Precision for 2x DP.  (RWC)
C***END PROLOGUE  MPBLAS
      COMMON /MPCOM/ MPB, MPT, MPM, MPLUN, MPMXR, MPR(30)
C***FIRST EXECUTABLE STATEMENT  MPBLAS
      I1 = 1
C
C     For full extended precision accuracy, MPB should be as large as
C     possible, subject to the restrictions in Brent's paper.
C
C     Statements below are for an integer wordlength of  48, 36, 32,
C     24, 18, and 16.  Pick one, or generate a new one.
C       48     MPB = 4194304
C       36     MPB =   65536
C       32     MPB =   16384
C       24     MPB =    1024
C       18     MPB =     128
C       16     MPB =      64
C
      MPBEXP = I1MACH(8)/2-2
      MPB = 2**MPBEXP
C
C     Set up remaining parameters
C                  UNIT FOR ERROR MESSAGES
      MPLUN = I1MACH(4)
C                  NUMBER OF MP DIGITS
      MPT = (2*I1MACH(14)+MPBEXP-1)/MPBEXP
C                  DIMENSION OF R
      MPMXR = MPT+4
C
      if (MPMXR.GT.30) THEN
         CALL XERMSG('SLATEC', 'MPBLAS',
     *      'Array space not sufficient for Quad Precision 2x ' //
     *      'Double Precision, Proceeding.', 1, 1)
         MPT = 26
         MPMXR = 30
      ENDIF
C                  EXPONENT RANGE
      MPM = MIN(32767,I1MACH(9)/4-1)
      RETURN
      END
