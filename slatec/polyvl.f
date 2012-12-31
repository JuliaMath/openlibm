*DECK POLYVL
      SUBROUTINE POLYVL (NDER, XX, YFIT, YP, N, X, C, WORK, IERR)
C***BEGIN PROLOGUE  POLYVL
C***PURPOSE  Calculate the value of a polynomial and its first NDER
C            derivatives where the polynomial was produced by a previous
C            call to POLINT.
C***LIBRARY   SLATEC
C***CATEGORY  E3
C***TYPE      SINGLE PRECISION (POLYVL-S, DPOLVL-D)
C***KEYWORDS  POLYNOMIAL EVALUATION
C***AUTHOR  Huddleston, R. E., (SNLL)
C***DESCRIPTION
C
C     Written by Robert E. Huddleston, Sandia Laboratories, Livermore
C
C     Abstract -
C        Subroutine POLYVL calculates the value of the polynomial and
C     its first NDER derivatives where the polynomial was produced by
C     a previous call to POLINT.
C        The variable N and the arrays X and C must not be altered
C     between the call to POLINT and the call to POLYVL.
C
C     ******  Dimensioning Information *******
C
C     YP   must be dimensioned by at least NDER
C     X    must be dimensioned by at least N (see the abstract )
C     C    must be dimensioned by at least N (see the abstract )
C     WORK must be dimensioned by at least 2*N if NDER is .GT. 0.
C
C     *** Note ***
C       If NDER=0, neither YP nor WORK need to be dimensioned variables.
C       If NDER=1, YP does not need to be a dimensioned variable.
C
C
C     *****  Input parameters
C
C     NDER - the number of derivatives to be evaluated
C
C     XX   - the argument at which the polynomial and its derivatives
C            are to be evaluated.
C
C     N    - *****
C            *       N, X, and C must not be altered between the call
C     X    - *       to POLINT and the call to POLYVL.
C     C    - *****
C
C
C     *****  Output Parameters
C
C     YFIT - the value of the polynomial at XX
C
C     YP   - the derivatives of the polynomial at XX.  The derivative of
C            order J at XX is stored in  YP(J) , J = 1,...,NDER.
C
C     IERR - Output error flag with the following possible values.
C          = 1  indicates normal execution
C
C     ***** Storage Parameters
C
C     WORK  = this is an array to provide internal working storage for
C             POLYVL.  It must be dimensioned by at least 2*N if NDER is
C             .GT. 0.  If NDER=0, WORK does not need to be a dimensioned
C             variable.
C
C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
C                 Curve fitting by polynomials in one variable, Report
C                 SLA-74-0270, Sandia Laboratories, June 1974.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   740601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  POLYVL
      DIMENSION  YP(*),X(*),C(*),WORK(*)
C***FIRST EXECUTABLE STATEMENT  POLYVL
      IERR=1
         IF (NDER.GT.0) GO TO 10020
C
C     *****   CODING FOR THE CASE NDER = 0
C
      PIONE=1.0
      PONE=C(1)
      YFIT=PONE
      IF (N.EQ.1) RETURN
      DO 10010 K=2,N
      PITWO=(XX-X(K-1))*PIONE
      PIONE=PITWO
      PTWO=PONE+PITWO*C(K)
      PONE=PTWO
10010 CONTINUE
      YFIT=PTWO
      RETURN
C
C     *****   END OF NDER = 0 CASE
C
10020 CONTINUE
         IF (N.GT.1) GO TO 10040
      YFIT=C(1)
C
C     *****  CODING FOR THE CASE  N=1 AND NDER .GT. 0
C
      DO 10030 K=1,NDER
      YP(K)=0.0
10030 CONTINUE
      RETURN
C
C     *****  END OF THE CASE  N = 1 AND  NDER .GT. 0
C
10040 CONTINUE
         IF (NDER.LT.N) GO TO 10050
C
C     *****  SET FLAGS FOR NUMBER OF DERIVATIVES AND FOR DERIVATIVES
C            IN EXCESS OF THE DEGREE (N-1) OF THE POLYNOMIAL.
C
      IZERO=1
      NDR=N-1
         GO TO 10060
10050 CONTINUE
      IZERO=0
      NDR=NDER
10060 CONTINUE
      M=NDR+1
      MM=M
C
C     *****  START OF THE CASE NDER .GT. 0  AND N .GT. 1
C     *****  THE POLYNOMIAL AND ITS DERIVATIVES WILL BE EVALUATED AT XX
C
      DO 10070 K=1,NDR
      YP(K)=C(K+1)
10070 CONTINUE
C
C     *****  THE FOLLOWING SECTION OF CODE IS EASIER TO READ IF ONE
C            BREAKS WORK INTO TWO ARRAYS W AND V. THE CODE WOULD THEN
C            READ
C                W(1) = 1.
C                PONE = C(1)
C               *DO   K = 2,N
C               *   V(K-1) =  XX - X(K-1)
C               *   W(K)   =  V(K-1)*W(K-1)
C               *   PTWO   =  PONE + W(K)*C(K)
C               *   PONE   =  PWO
C
C               YFIT = PTWO
C
      WORK(1)=1.0
      PONE=C(1)
      DO 10080 K=2,N
      KM1=K-1
      NPKM1=N+K-1
      WORK(NPKM1)=XX-X(KM1)
      WORK(K)=WORK(NPKM1)*WORK(KM1)
      PTWO=PONE+WORK(K)*C(K)
      PONE=PTWO
10080 CONTINUE
      YFIT=PTWO
C
C     ** AT THIS POINT THE POLYNOMIAL HAS BEEN EVALUATED AND INFORMATION
C        FOR THE DERIVATIVE EVALUATIONS HAVE BEEN STORED IN THE ARRAY
C        WORK
         IF (N.EQ.2) GO TO 10110
      IF (M.EQ.N) MM=NDR
C
C     ***** EVALUATE THE DERIVATIVES AT XX
C
C                  ******  DO K=2,MM   (FOR MOST CASES, MM = NDER + 1)
C                  *  ******  DO I=2,N-K+1
C                  *  *       W(I) = V(K-2+I)*W(I-1) + W(I)
C                  *  *       YP(K-1) = YP(K-1) + W(I)*C(K-1+I)
C                  ******  CONTINUE
C
      DO 10090 K=2,MM
      NMKP1=N-K+1
      KM1=K-1
      KM2PN=K-2+N
      DO 10090 I=2,NMKP1
      KM2PNI=KM2PN+I
      IM1=I-1
      KM1PI=KM1+I
      WORK(I)=WORK(KM2PNI)*WORK(IM1)+WORK(I)
      YP(KM1)=YP(KM1)+WORK(I)*C(KM1PI)
10090 CONTINUE
         IF (NDR.EQ.1) GO TO 10110
      FAC=1.0
      DO 10100 K=2,NDR
      XK=K
      FAC=XK*FAC
      YP(K)=FAC*YP(K)
10100 CONTINUE
C
C     ***** END OF DERIVATIVE EVALUATIONS
C
10110 CONTINUE
      IF (IZERO.EQ.0) RETURN
C
C     *****  SET EXCESS DERIVATIVES TO ZERO.
C
      DO 10120 K=N,NDER
      YP(K)=0.0
10120 CONTINUE
      RETURN
      END
