*DECK CPEVL
      SUBROUTINE CPEVL (N, M, A, Z, C, B, KBD)
C***BEGIN PROLOGUE  CPEVL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to CPZERO
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (CPEVL-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C        Evaluate a complex polynomial and its derivatives.
C        Optionally compute error bounds for these values.
C
C   INPUT...
C        N = Degree of the polynomial
C        M = Number of derivatives to be calculated,
C            M=0 evaluates only the function
C            M=1 evaluates the function and first derivative, etc.
C             if M .GT. N+1 function and all N derivatives will be
C                calculated.
C       A = Complex vector containing the N+1 coefficients of polynomial
C               A(I)= coefficient of Z**(N+1-I)
C        Z = Complex point at which the evaluation is to take place.
C        C = Array of 2(M+1) words into which values are placed.
C        B = Array of 2(M+1) words only needed if bounds are to be
C              calculated.  It is not used otherwise.
C        KBD = A logical variable, e.g. .TRUE. or .FALSE. which is
C              to be set .TRUE. if bounds are to be computed.
C
C  OUTPUT...
C        C =  C(I+1) contains the complex value of the I-th
C              derivative at Z, I=0,...,M
C        B =  B(I) contains the bounds on the real and imaginary parts
C              of C(I) if they were requested.
C
C***SEE ALSO  CPZERO
C***ROUTINES CALLED  I1MACH
C***REVISION HISTORY  (YYMMDD)
C   810223  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  CPEVL
C
      COMPLEX A(*),C(*),Z,CI,CIM1,B(*),BI,BIM1,T,ZA,Q
      LOGICAL KBD
      SAVE D1
      DATA D1 /0.0/
      ZA(Q)=CMPLX(ABS(REAL(Q)),ABS(AIMAG(Q)))
C***FIRST EXECUTABLE STATEMENT  CPEVL
      IF (D1 .EQ. 0.0) D1 = REAL(I1MACH(10))**(1-I1MACH(11))
      NP1=N+1
      DO 1 J=1,NP1
         CI=0.0
         CIM1=A(J)
         BI=0.0
         BIM1=0.0
         MINI=MIN(M+1,N+2-J)
            DO 1 I=1,MINI
               IF(J .NE. 1) CI=C(I)
               IF(I .NE. 1) CIM1=C(I-1)
               C(I)=CIM1+Z*CI
               IF(.NOT. KBD) GO TO 1
               IF(J .NE. 1) BI=B(I)
               IF(I .NE. 1) BIM1=B(I-1)
               T=BI+(3.*D1+4.*D1*D1)*ZA(CI)
               R=REAL(ZA(Z)*CMPLX(REAL(T),-AIMAG(T)))
               S=AIMAG(ZA(Z)*T)
               B(I)=(1.+8.*D1)*(BIM1+D1*ZA(CIM1)+CMPLX(R,S))
               IF(J .EQ. 1) B(I)=0.0
    1 CONTINUE
      RETURN
      END
