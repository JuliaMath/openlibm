*DECK DPPERM
      SUBROUTINE DPPERM (DX, N, IPERM, IER)
C***BEGIN PROLOGUE  DPPERM
C***PURPOSE  Rearrange a given array according to a prescribed
C            permutation vector.
C***LIBRARY   SLATEC
C***CATEGORY  N8
C***TYPE      DOUBLE PRECISION (SPPERM-S, DPPERM-D, IPPERM-I, HPPERM-H)
C***KEYWORDS  PERMUTATION, REARRANGEMENT
C***AUTHOR  McClain, M. A., (NIST)
C           Rhoads, G. S., (NBS)
C***DESCRIPTION
C
C         DPPERM rearranges the data vector DX according to the
C         permutation IPERM: DX(I) <--- DX(IPERM(I)).  IPERM could come
C         from one of the sorting routines IPSORT, SPSORT, DPSORT or
C         HPSORT.
C
C     Description of Parameters
C         DX - input/output -- double precision array of values to be
C                   rearranged.
C         N - input -- number of values in double precision array DX.
C         IPERM - input -- permutation vector.
C         IER - output -- error indicator:
C             =  0  if no error,
C             =  1  if N is zero or negative,
C             =  2  if IPERM is not a valid permutation.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   901004  DATE WRITTEN
C   920507  Modified by M. McClain to revise prologue text.
C***END PROLOGUE  DPPERM
      INTEGER N, IPERM(*), I, IER, INDX, INDX0, ISTRT
      DOUBLE PRECISION DX(*), DTEMP
C***FIRST EXECUTABLE STATEMENT  DPPERM
      IER=0
      IF(N.LT.1)THEN
         IER=1
         CALL XERMSG ('SLATEC', 'DPPERM',
     +    'The number of values to be rearranged, N, is not positive.',
     +    IER, 1)
         RETURN
      ENDIF
C
C     CHECK WHETHER IPERM IS A VALID PERMUTATION
C
      DO 100 I=1,N
         INDX=ABS(IPERM(I))
         IF((INDX.GE.1).AND.(INDX.LE.N))THEN
            IF(IPERM(INDX).GT.0)THEN
               IPERM(INDX)=-IPERM(INDX)
               GOTO 100
            ENDIF
         ENDIF
         IER=2
         CALL XERMSG ('SLATEC', 'DPPERM',
     +    'The permutation vector, IPERM, is not valid.', IER, 1)
         RETURN
  100 CONTINUE
C
C     REARRANGE THE VALUES OF DX
C
C     USE THE IPERM VECTOR AS A FLAG.
C     IF IPERM(I) > 0, THEN THE I-TH VALUE IS IN CORRECT LOCATION
C
      DO 330 ISTRT = 1 , N
         IF (IPERM(ISTRT) .GT. 0) GOTO 330
         INDX = ISTRT
         INDX0 = INDX
         DTEMP = DX(ISTRT)
  320    CONTINUE
         IF (IPERM(INDX) .GE. 0) GOTO 325
            DX(INDX) = DX(-IPERM(INDX))
            INDX0 = INDX
            IPERM(INDX) = -IPERM(INDX)
            INDX = IPERM(INDX)
            GOTO 320
  325    CONTINUE
         DX(INDX0) = DTEMP
  330 CONTINUE
C
      RETURN
      END
