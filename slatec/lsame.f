*DECK LSAME
      LOGICAL FUNCTION LSAME (CA, CB)
C***BEGIN PROLOGUE  LSAME
C***SUBSIDIARY
C***PURPOSE  Test two characters to determine if they are the same
C            letter, except for case.
C***LIBRARY   SLATEC
C***CATEGORY  R, N3
C***TYPE      LOGICAL (LSAME-L)
C***KEYWORDS  CHARACTER COMPARISON, LEVEL 2 BLAS, LEVEL 3 BLAS
C***AUTHOR  Hanson, R., (SNLA)
C           Du Croz, J., (NAG)
C***DESCRIPTION
C
C  LSAME  tests if CA is the same letter as CB regardless of case.
C  CB is assumed to be an upper case letter. LSAME returns .TRUE. if
C  CA is either the same as CB or the equivalent lower case letter.
C
C  N.B. This version of the code is correct for both ASCII and EBCDIC
C       systems.  Installers must modify the routine for other
C       character-codes.
C
C       For CDC systems using 6-12 bit representations, the system-
C       specific code in comments must be activated.
C
C  Parameters
C  ==========
C
C  CA     - CHARACTER*1
C  CB     - CHARACTER*1
C           On entry, CA and CB specify characters to be compared.
C           Unchanged on exit.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   860720  DATE WRITTEN
C   910606  Modified to meet SLATEC prologue standards.  Only comment
C           lines were modified.  (BKS)
C   910607  Modified to handle ASCII and EBCDIC codes.  (WRB)
C   930201  Tests for equality and equivalence combined.  (RWC and WRB)
C***END PROLOGUE  LSAME
C     .. Scalar Arguments ..
      CHARACTER CA*1, CB*1
C     .. Local Scalars ..
      INTEGER IOFF
      LOGICAL FIRST
C     .. Intrinsic Functions ..
      INTRINSIC ICHAR
C     .. Save statement ..
      SAVE FIRST, IOFF
C     .. Data statements ..
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  LSAME
      IF (FIRST) IOFF = ICHAR('a') - ICHAR('A')
C
      FIRST = .FALSE.
C
C     Test if the characters are equal or equivalent.
C
      LSAME = (CA.EQ.CB) .OR. (ICHAR(CA)-IOFF.EQ.ICHAR(CB))
C
      RETURN
C
C  The following comments contain code for CDC systems using 6-12 bit
C  representations.
C
C     .. Parameters ..
C     INTEGER                ICIRFX
C     PARAMETER            ( ICIRFX=62 )
C     .. Scalar Arguments ..
C     CHARACTER*1            CB
C     .. Array Arguments ..
C     CHARACTER*1            CA(*)
C     .. Local Scalars ..
C     INTEGER                IVAL
C     .. Intrinsic Functions ..
C     INTRINSIC              ICHAR, CHAR
C     .. Executable Statements ..
C     INTRINSIC              ICHAR, CHAR
C
C     See if the first character in string CA equals string CB.
C
C     LSAME = CA(1) .EQ. CB .AND. CA(1) .NE. CHAR(ICIRFX)
C
C     IF (LSAME) RETURN
C
C     The characters are not identical. Now check them for equivalence.
C     Look for the 'escape' character, circumflex, followed by the
C     letter.
C
C     IVAL = ICHAR(CA(2))
C     IF (IVAL.GE.ICHAR('A') .AND. IVAL.LE.ICHAR('Z')) THEN
C        LSAME = CA(1) .EQ. CHAR(ICIRFX) .AND. CA(2) .EQ. CB
C     ENDIF
C
C     RETURN
C
C     End of LSAME.
C
      END
