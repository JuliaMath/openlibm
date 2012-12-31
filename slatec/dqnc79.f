*DECK DQNC79
      SUBROUTINE DQNC79 (FUN, A, B, ERR, ANS, IERR, K)
C***BEGIN PROLOGUE  DQNC79
C***PURPOSE  Integrate a function using a 7-point adaptive Newton-Cotes
C            quadrature rule.
C***LIBRARY   SLATEC
C***CATEGORY  H2A1A1
C***TYPE      DOUBLE PRECISION (QNC79-S, DQNC79-D)
C***KEYWORDS  ADAPTIVE QUADRATURE, INTEGRATION, NEWTON-COTES
C***AUTHOR  Kahaner, D. K., (NBS)
C           Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract  *** a DOUBLE PRECISION routine ***
C       DQNC79 is a general purpose program for evaluation of
C       one dimensional integrals of user defined functions.
C       DQNC79 will pick its own points for evaluation of the
C       integrand and these will vary from problem to problem.
C       Thus, DQNC79 is not designed to integrate over data sets.
C       Moderately smooth integrands will be integrated efficiently
C       and reliably.  For problems with strong singularities,
C       oscillations etc., the user may wish to use more sophis-
C       ticated routines such as those in QUADPACK.  One measure
C       of the reliability of DQNC79 is the output parameter K,
C       giving the number of integrand evaluations that were needed.
C
C     Description of Arguments
C
C     --Input--* FUN, A, B, ERR are DOUBLE PRECISION *
C       FUN  - name of external function to be integrated.  This name
C              must be in an EXTERNAL statement in your calling
C              program.  You must write a Fortran function to evaluate
C              FUN.  This should be of the form
C                    DOUBLE PRECISION FUNCTION FUN (X)
C              C
C              C     X can vary from A to B
C              C     FUN(X) should be finite for all X on interval.
C              C
C                    FUN = ...
C                    RETURN
C                    END
C       A    - lower limit of integration
C       B    - upper limit of integration (may be less than A)
C       ERR  - is a requested error tolerance.  Normally, pick a value
C              0 .LT. ERR .LT. 1.0D-8.
C
C     --Output--
C       ANS  - computed value of the integral.  Hopefully, ANS is
C              accurate to within ERR * integral of ABS(FUN(X)).
C       IERR - a status code
C            - Normal codes
C               1  ANS most likely meets requested error tolerance.
C              -1  A equals B, or A and B are too nearly equal to
C                  allow normal integration.  ANS is set to zero.
C            - Abnormal code
C               2  ANS probably does not meet requested error tolerance.
C       K    - the number of function evaluations actually used to do
C              the integration.  A value of K .GT. 1000 indicates a
C              difficult problem; other programs may be more efficient.
C              DQNC79 will gracefully give up if K exceeds 2000.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, I1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920218  Code redone to parallel QNC79.  (WRB)
C   930120  Increase array size 80->99, and KMX 2000->5000 for SUN -r8
C           wordlength.  (RWC)
C***END PROLOGUE  DQNC79
C     .. Scalar Arguments ..
      DOUBLE PRECISION A, ANS, B, ERR
      INTEGER IERR, K
C     .. Function Arguments ..
      DOUBLE PRECISION FUN
      EXTERNAL FUN
C     .. Local Scalars ..
      DOUBLE PRECISION AE, AREA, BANK, BLOCAL, C, CE, EE, EF, EPS, Q13,
     +                 Q7, Q7L, SQ2, TEST, TOL, VR, W1, W2, W3, W4
      INTEGER I, KML, KMX, L, LMN, LMX, NBITS, NIB, NLMN, NLMX
      LOGICAL FIRST
C     .. Local Arrays ..
      DOUBLE PRECISION AA(99), F(13), F1(99), F2(99), F3(99), F4(99),
     +                 F5(99), F6(99), F7(99), HH(99), Q7R(99), VL(99)
      INTEGER LR(99)
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      INTEGER I1MACH
      EXTERNAL D1MACH, I1MACH
C     .. External Subroutines ..
      EXTERNAL XERMSG
C     .. Intrinsic Functions ..
      INTRINSIC ABS, LOG, MAX, MIN, SIGN, SQRT
C     .. Save statement ..
      SAVE NBITS, NLMX, FIRST, SQ2, W1, W2, W3, W4
C     .. Data statements ..
      DATA KML /7/, KMX /5000/, NLMN /2/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DQNC79
      IF (FIRST) THEN
        W1 = 41.0D0/140.0D0
        W2 = 216.0D0/140.0D0
        W3 = 27.0D0/140.0D0
        W4 = 272.0D0/140.0D0
        NBITS = D1MACH(5)*I1MACH(14)/0.30102000D0
        NLMX = MIN(99,(NBITS*4)/5)
        SQ2 = SQRT(2.0D0)
      ENDIF
      FIRST = .FALSE.
      ANS = 0.0D0
      IERR = 1
      CE = 0.0D0
      IF (A .EQ. B) GO TO 260
      LMX = NLMX
      LMN = NLMN
      IF (B .EQ. 0.0D0) GO TO 100
      IF (SIGN(1.0D0,B)*A .LE. 0.0D0) GO TO 100
      C = ABS(1.0D0-A/B)
      IF (C .GT. 0.1D0) GO TO 100
      IF (C .LE. 0.0D0) GO TO 260
      NIB = 0.5D0 - LOG(C)/LOG(2.0D0)
      LMX = MIN(NLMX,NBITS-NIB-4)
      IF (LMX .LT. 2) GO TO 260
      LMN = MIN(LMN,LMX)
  100 TOL = MAX(ABS(ERR),2.0D0**(5-NBITS))
      IF (ERR .EQ. 0.0D0) TOL = SQRT(D1MACH(4))
      EPS = TOL
      HH(1) = (B-A)/12.0D0
      AA(1) = A
      LR(1) = 1
      DO 110 I = 1,11,2
        F(I) = FUN(A+(I-1)*HH(1))
  110 CONTINUE
      BLOCAL = B
      F(13) = FUN(BLOCAL)
      K = 7
      L = 1
      AREA = 0.0D0
      Q7 = 0.0D0
      EF = 256.0D0/255.0D0
      BANK = 0.0D0
C
C     Compute refined estimates, estimate the error, etc.
C
  120 DO 130 I = 2,12,2
        F(I) = FUN(AA(L)+(I-1)*HH(L))
  130 CONTINUE
      K = K + 6
C
C     Compute left and right half estimates
C
      Q7L = HH(L)*((W1*(F(1)+F(7))+W2*(F(2)+F(6)))+
     +      (W3*(F(3)+F(5))+W4*F(4)))
      Q7R(L) = HH(L)*((W1*(F(7)+F(13))+W2*(F(8)+F(12)))+
     +         (W3*(F(9)+F(11))+W4*F(10)))
C
C     Update estimate of integral of absolute value
C
      AREA = AREA + (ABS(Q7L)+ABS(Q7R(L))-ABS(Q7))
C
C     Do not bother to test convergence before minimum refinement level
C
      IF (L .LT. LMN) GO TO 180
C
C     Estimate the error in new value for whole interval, Q13
C
      Q13 = Q7L + Q7R(L)
      EE = ABS(Q7-Q13)*EF
C
C     Compute nominal allowed error
C
      AE = EPS*AREA
C
C     Borrow from bank account, but not too much
C
      TEST = MIN(AE+0.8D0*BANK,10.0D0*AE)
C
C     Don't ask for excessive accuracy
C
      TEST = MAX(TEST,TOL*ABS(Q13),0.00003D0*TOL*AREA)
C
C     Now, did this interval pass or not?
C
      IF (EE-TEST) 150,150,170
C
C     Have hit maximum refinement level -- penalize the cumulative error
C
  140 CE = CE + (Q7-Q13)
      GO TO 160
C
C     On good intervals accumulate the theoretical estimate
C
  150 CE = CE + (Q7-Q13)/255.0D0
C
C     Update the bank account.  Don't go into debt.
C
  160 BANK = BANK + (AE-EE)
      IF (BANK .LT. 0.0D0) BANK = 0.0D0
C
C     Did we just finish a left half or a right half?
C
      IF (LR(L)) 190,190,210
C
C     Consider the left half of next deeper level
C
  170 IF (K .GT. KMX) LMX = MIN(KML,LMX)
      IF (L .GE. LMX) GO TO 140
  180 L = L + 1
      EPS = EPS*0.5D0
      IF (L .LE. 17) EF = EF/SQ2
      HH(L) = HH(L-1)*0.5D0
      LR(L) = -1
      AA(L) = AA(L-1)
      Q7 = Q7L
      F1(L) = F(7)
      F2(L) = F(8)
      F3(L) = F(9)
      F4(L) = F(10)
      F5(L) = F(11)
      F6(L) = F(12)
      F7(L) = F(13)
      F(13) = F(7)
      F(11) = F(6)
      F(9) = F(5)
      F(7) = F(4)
      F(5) = F(3)
      F(3) = F(2)
      GO TO 120
C
C     Proceed to right half at this level
C
  190 VL(L) = Q13
  200 Q7 = Q7R(L-1)
      LR(L) = 1
      AA(L) = AA(L) + 12.0D0*HH(L)
      F(1) = F1(L)
      F(3) = F2(L)
      F(5) = F3(L)
      F(7) = F4(L)
      F(9) = F5(L)
      F(11) = F6(L)
      F(13) = F7(L)
      GO TO 120
C
C     Left and right halves are done, so go back up a level
C
  210 VR = Q13
  220 IF (L .LE. 1) GO TO 250
      IF (L .LE. 17) EF = EF*SQ2
      EPS = EPS*2.0D0
      L = L - 1
      IF (LR(L)) 230,230,240
  230 VL(L) = VL(L+1) + VR
      GO TO 200
  240 VR = VL(L+1) + VR
      GO TO 220
C
C     Exit
C
  250 ANS = VR
      IF (ABS(CE) .LE. 2.0D0*TOL*AREA) GO TO 270
      IERR = 2
      CALL XERMSG ('SLATEC', 'DQNC79',
     +   'ANS is probably insufficiently accurate.', 2, 1)
      GO TO 270
  260 IERR = -1
      CALL XERMSG ('SLATEC', 'DQNC79',
     +   'A and B are too nearly equal to allow normal integration. $$'
     +   // 'ANS is set to zero and IERR to -1.', -1, -1)
  270 RETURN
      END
