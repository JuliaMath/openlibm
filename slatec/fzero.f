*DECK FZERO
      SUBROUTINE FZERO (F, B, C, R, RE, AE, IFLAG)
C***BEGIN PROLOGUE  FZERO
C***PURPOSE  Search for a zero of a function F(X) in a given interval
C            (B,C).  It is designed primarily for problems where F(B)
C            and F(C) have opposite signs.
C***LIBRARY   SLATEC
C***CATEGORY  F1B
C***TYPE      SINGLE PRECISION (FZERO-S, DFZERO-D)
C***KEYWORDS  BISECTION, NONLINEAR EQUATIONS, ROOTS, ZEROS
C***AUTHOR  Shampine, L. F., (SNLA)
C           Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     FZERO searches for a zero of a REAL function F(X) between the
C     given REAL values B and C until the width of the interval (B,C)
C     has collapsed to within a tolerance specified by the stopping
C     criterion,
C        ABS(B-C) .LE. 2.*(RW*ABS(B)+AE).
C     The method used is an efficient combination of bisection and the
C     secant rule and is due to T. J. Dekker.
C
C     Description Of Arguments
C
C   F     :EXT   - Name of the REAL external function.  This name must
C                  be in an EXTERNAL statement in the calling program.
C                  F must be a function of one REAL argument.
C
C   B     :INOUT - One end of the REAL interval (B,C).  The value
C                  returned for B usually is the better approximation
C                  to a zero of F.
C
C   C     :INOUT - The other end of the REAL interval (B,C)
C
C   R     :OUT   - A (better) REAL guess of a zero of F which could help
C                  in speeding up convergence.  If F(B) and F(R) have
C                  opposite signs, a root will be found in the interval
C                  (B,R); if not, but F(R) and F(C) have opposite signs,
C                  a root will be found in the interval (R,C);
C                  otherwise, the interval (B,C) will be searched for a
C                  possible root.  When no better guess is known, it is
C                  recommended that r be set to B or C, since if R is
C                  not interior to the interval (B,C), it will be
C                  ignored.
C
C   RE    :IN    - Relative error used for RW in the stopping criterion.
C                  If the requested RE is less than machine precision,
C                  then RW is set to approximately machine precision.
C
C   AE    :IN    - Absolute error used in the stopping criterion.  If
C                  the given interval (B,C) contains the origin, then a
C                  nonzero value should be chosen for AE.
C
C   IFLAG :OUT   - A status code.  User must check IFLAG after each
C                  call.  Control returns to the user from FZERO in all
C                  cases.
C
C                1  B is within the requested tolerance of a zero.
C                   The interval (B,C) collapsed to the requested
C                   tolerance, the function changes sign in (B,C), and
C                   F(X) decreased in magnitude as (B,C) collapsed.
C
C                2  F(B) = 0.  However, the interval (B,C) may not have
C                   collapsed to the requested tolerance.
C
C                3  B may be near a singular point of F(X).
C                   The interval (B,C) collapsed to the requested tol-
C                   erance and the function changes sign in (B,C), but
C                   F(X) increased in magnitude as (B,C) collapsed, i.e.
C                     ABS(F(B out)) .GT. MAX(ABS(F(B in)),ABS(F(C in)))
C
C                4  No change in sign of F(X) was found although the
C                   interval (B,C) collapsed to the requested tolerance.
C                   The user must examine this case and decide whether
C                   B is near a local minimum of F(X), or B is near a
C                   zero of even multiplicity, or neither of these.
C
C                5  Too many (.GT. 500) function evaluations used.
C
C***REFERENCES  L. F. Shampine and H. A. Watts, FZERO, a root-solving
C                 code, Report SC-TM-70-631, Sandia Laboratories,
C                 September 1970.
C               T. J. Dekker, Finding a zero by means of successive
C                 linear interpolation, Constructive Aspects of the
C                 Fundamental Theorem of Algebra, edited by B. Dejon
C                 and P. Henrici, Wiley-Interscience, 1969.
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   700901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  FZERO
      REAL A,ACBS,ACMB,AE,AW,B,C,CMB,ER,FA,FB,FC,FX,FZ,P,Q,R,
     +     RE,RW,T,TOL,Z
      INTEGER IC,IFLAG,KOUNT
C***FIRST EXECUTABLE STATEMENT  FZERO
C
C   ER is two times the computer unit roundoff value which is defined
C   here by the function R1MACH.
C
      ER = 2.0E0 * R1MACH(4)
C
C   Initialize.
C
      Z = R
      IF (R .LE. MIN(B,C)  .OR.  R .GE. MAX(B,C)) Z = C
      RW = MAX(RE,ER)
      AW = MAX(AE,0.E0)
      IC = 0
      T = Z
      FZ = F(T)
      FC = FZ
      T = B
      FB = F(T)
      KOUNT = 2
      IF (SIGN(1.0E0,FZ) .EQ. SIGN(1.0E0,FB)) GO TO 1
      C = Z
      GO TO 2
    1 IF (Z .EQ. C) GO TO 2
      T = C
      FC = F(T)
      KOUNT = 3
      IF (SIGN(1.0E0,FZ) .EQ. SIGN(1.0E0,FC)) GO TO 2
      B = Z
      FB = FZ
    2 A = C
      FA = FC
      ACBS = ABS(B-C)
      FX = MAX(ABS(FB),ABS(FC))
C
    3 IF (ABS(FC) .GE. ABS(FB)) GO TO 4
C
C   Perform interchange.
C
      A = B
      FA = FB
      B = C
      FB = FC
      C = A
      FC = FA
C
    4 CMB = 0.5E0*(C-B)
      ACMB = ABS(CMB)
      TOL = RW*ABS(B) + AW
C
C   Test stopping criterion and function count.
C
      IF (ACMB .LE. TOL) GO TO 10
      IF (FB .EQ. 0.E0) GO TO 11
      IF (KOUNT .GE. 500) GO TO 14
C
C   Calculate new iterate implicitly as B+P/Q, where we arrange
C   P .GE. 0.  The implicit form is used to prevent overflow.
C
      P = (B-A)*FB
      Q = FA - FB
      IF (P .GE. 0.E0) GO TO 5
      P = -P
      Q = -Q
C
C   Update A and check for satisfactory reduction in the size of the
C   bracketing interval.  If not, perform bisection.
C
    5 A = B
      FA = FB
      IC = IC + 1
      IF (IC .LT. 4) GO TO 6
      IF (8.0E0*ACMB .GE. ACBS) GO TO 8
      IC = 0
      ACBS = ACMB
C
C   Test for too small a change.
C
    6 IF (P .GT. ABS(Q)*TOL) GO TO 7
C
C   Increment by TOLerance.
C
      B = B + SIGN(TOL,CMB)
      GO TO 9
C
C   Root ought to be between B and (C+B)/2.
C
    7 IF (P .GE. CMB*Q) GO TO 8
C
C   Use secant rule.
C
      B = B + P/Q
      GO TO 9
C
C   Use bisection (C+B)/2.
C
    8 B = B + CMB
C
C   Have completed computation for new iterate B.
C
    9 T = B
      FB = F(T)
      KOUNT = KOUNT + 1
C
C   Decide whether next step is interpolation or extrapolation.
C
      IF (SIGN(1.0E0,FB) .NE. SIGN(1.0E0,FC)) GO TO 3
      C = A
      FC = FA
      GO TO 3
C
C   Finished.  Process results for proper setting of IFLAG.
C
   10 IF (SIGN(1.0E0,FB) .EQ. SIGN(1.0E0,FC)) GO TO 13
      IF (ABS(FB) .GT. FX) GO TO 12
      IFLAG = 1
      RETURN
   11 IFLAG = 2
      RETURN
   12 IFLAG = 3
      RETURN
   13 IFLAG = 4
      RETURN
   14 IFLAG = 5
      RETURN
      END
