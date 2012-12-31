*DECK DPOLCF
      SUBROUTINE DPOLCF (XX, N, X, C, D, WORK)
C***BEGIN PROLOGUE  DPOLCF
C***PURPOSE  Compute the coefficients of the polynomial fit (including
C            Hermite polynomial fits) produced by a previous call to
C            POLINT.
C***LIBRARY   SLATEC
C***CATEGORY  E1B
C***TYPE      DOUBLE PRECISION (POLCOF-S, DPOLCF-D)
C***KEYWORDS  COEFFICIENTS, POLYNOMIAL
C***AUTHOR  Huddleston, R. E., (SNLL)
C***DESCRIPTION
C
C     Abstract
C        Subroutine DPOLCF computes the coefficients of the polynomial
C     fit (including Hermite polynomial fits ) produced by a previous
C     call to DPLINT.  The coefficients of the polynomial, expanded
C     about XX, are stored in the array D. The expansion is of the form
C     P(Z) = D(1) + D(2)*(Z-XX) +D(3)*((Z-XX)**2) + ... +
C                                                  D(N)*((Z-XX)**(N-1)).
C     Between the call to DPLINT and the call to DPOLCF the variable N
C     and the arrays X and C must not be altered.
C
C     *****  INPUT PARAMETERS
C      *** All TYPE REAL variables are DOUBLE PRECISION ***
C
C     XX   - The point about which the Taylor expansion is to be made.
C
C     N    - ****
C            *     N, X, and C must remain unchanged between the
C     X    - *     call to DPLINT and the call to DPOLCF.
C     C    - ****
C
C     *****  OUTPUT PARAMETER
C      *** All TYPE REAL variables are DOUBLE PRECISION ***
C
C     D    - The array of coefficients for the Taylor expansion as
C            explained in the abstract
C
C     *****  STORAGE PARAMETER
C
C     WORK - This is an array to provide internal working storage. It
C            must be dimensioned by at least 2*N in the calling program.
C
C
C     **** Note - There are two methods for evaluating the fit produced
C     by DPLINT. You may call DPOLVL to perform the task, or you may
C     call DPOLCF to obtain the coefficients of the Taylor expansion and
C     then write your own evaluation scheme. Due to the inherent errors
C     in the computations of the Taylor expansion from the Newton
C     coefficients produced by DPLINT, much more accuracy may be
C     expected by calling DPOLVL as opposed to writing your own scheme.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   890213  DATE WRITTEN
C   891006  Cosmetic changes to prologue.  (WRB)
C   891024  Corrected KEYWORD section.  (WRB)
C   891024  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DPOLCF
C
      INTEGER I,IM1,K,KM1,KM1PI,KM2N,KM2NPI,N,NM1,NMKP1,NPKM1
      DOUBLE PRECISION C(*),D(*),PONE,PTWO,X(*),XX,WORK(*)
C***FIRST EXECUTABLE STATEMENT  DPOLCF
      DO 10010 K=1,N
      D(K)=C(K)
10010 CONTINUE
      IF (N.EQ.1) RETURN
      WORK(1)=1.0D0
      PONE=C(1)
      NM1=N-1
      DO 10020 K=2,N
      KM1=K-1
      NPKM1=N+K-1
      WORK(NPKM1)=XX-X(KM1)
      WORK(K)=WORK(NPKM1)*WORK(KM1)
      PTWO=PONE+WORK(K)*C(K)
      PONE=PTWO
10020 CONTINUE
      D(1)=PTWO
      IF (N.EQ.2) RETURN
      DO 10030 K=2,NM1
      KM1=K-1
      KM2N=K-2+N
      NMKP1=N-K+1
      DO 10030 I=2,NMKP1
      KM2NPI=KM2N+I
      IM1=I-1
      KM1PI=KM1+I
      WORK(I)=WORK(KM2NPI)*WORK(IM1)+WORK(I)
      D(K)=D(K)+WORK(I)*D(KM1PI)
10030 CONTINUE
      RETURN
      END
