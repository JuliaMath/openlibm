*DECK DCFOD
      SUBROUTINE DCFOD (METH, ELCO, TESCO)
C***BEGIN PROLOGUE  DCFOD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DDEBDF
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (CFOD-S, DCFOD-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   DCFOD defines coefficients needed in the integrator package DDEBDF
C
C***SEE ALSO  DDEBDF
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820301  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DCFOD
C
C
      INTEGER I, IB, METH, NQ, NQM1, NQP1
      DOUBLE PRECISION AGAMQ, ELCO, FNQ, FNQM1, PC, PINT, RAGQ,
     1      RQ1FAC, RQFAC, TESCO, TSIGN, XPIN
      DIMENSION ELCO(13,12),TESCO(3,12)
C     ------------------------------------------------------------------
C      DCFOD  IS CALLED BY THE INTEGRATOR ROUTINE TO SET COEFFICIENTS
C      NEEDED THERE.  THE COEFFICIENTS FOR THE CURRENT METHOD, AS
C      GIVEN BY THE VALUE OF METH, ARE SET FOR ALL ORDERS AND SAVED.
C      THE MAXIMUM ORDER ASSUMED HERE IS 12 IF METH = 1 AND 5 IF METH =
C      2.  (A SMALLER VALUE OF THE MAXIMUM ORDER IS ALSO ALLOWED.)
C      DCFOD  IS CALLED ONCE AT THE BEGINNING OF THE PROBLEM,
C      AND IS NOT CALLED AGAIN UNLESS AND UNTIL METH IS CHANGED.
C
C      THE ELCO ARRAY CONTAINS THE BASIC METHOD COEFFICIENTS.
C      THE COEFFICIENTS EL(I), 1 .LE. I .LE. NQ+1, FOR THE METHOD OF
C      ORDER NQ ARE STORED IN ELCO(I,NQ).  THEY ARE GIVEN BY A
C      GENERATING POLYNOMIAL, I.E.,
C          L(X) = EL(1) + EL(2)*X + ... + EL(NQ+1)*X**NQ.
C      FOR THE IMPLICIT ADAMS METHODS, L(X) IS GIVEN BY
C          DL/DX = (X+1)*(X+2)*...*(X+NQ-1)/FACTORIAL(NQ-1),    L(-1) =
C      0.  FOR THE BDF METHODS, L(X) IS GIVEN BY
C          L(X) = (X+1)*(X+2)* ... *(X+NQ)/K,
C      WHERE         K = FACTORIAL(NQ)*(1 + 1/2 + ... + 1/NQ).
C
C      THE TESCO ARRAY CONTAINS TEST CONSTANTS USED FOR THE
C      LOCAL ERROR TEST AND THE SELECTION OF STEP SIZE AND/OR ORDER.
C      AT ORDER NQ, TESCO(K,NQ) IS USED FOR THE SELECTION OF STEP
C      SIZE AT ORDER NQ - 1 IF K = 1, AT ORDER NQ IF K = 2, AND AT ORDER
C      NQ + 1 IF K = 3.
C     ------------------------------------------------------------------
      DIMENSION PC(12)
C
C***FIRST EXECUTABLE STATEMENT  DCFOD
      GO TO (10,60), METH
C
   10 CONTINUE
         ELCO(1,1) = 1.0D0
         ELCO(2,1) = 1.0D0
         TESCO(1,1) = 0.0D0
         TESCO(2,1) = 2.0D0
         TESCO(1,2) = 1.0D0
         TESCO(3,12) = 0.0D0
         PC(1) = 1.0D0
         RQFAC = 1.0D0
         DO 50 NQ = 2, 12
C           ------------------------------------------------------------
C            THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE
C                POLYNOMIAL P(X) = (X+1)*(X+2)*...*(X+NQ-1).
C            INITIALLY, P(X) = 1.
C           ------------------------------------------------------------
            RQ1FAC = RQFAC
            RQFAC = RQFAC/NQ
            NQM1 = NQ - 1
            FNQM1 = NQM1
            NQP1 = NQ + 1
C           FORM COEFFICIENTS OF P(X)*(X+NQ-1).
C           ----------------------------------
            PC(NQ) = 0.0D0
            DO 20 IB = 1, NQM1
               I = NQP1 - IB
               PC(I) = PC(I-1) + FNQM1*PC(I)
   20       CONTINUE
            PC(1) = FNQM1*PC(1)
C           COMPUTE INTEGRAL, -1 TO 0, OF P(X) AND X*P(X).
C           -----------------------
            PINT = PC(1)
            XPIN = PC(1)/2.0D0
            TSIGN = 1.0D0
            DO 30 I = 2, NQ
               TSIGN = -TSIGN
               PINT = PINT + TSIGN*PC(I)/I
               XPIN = XPIN + TSIGN*PC(I)/(I+1)
   30       CONTINUE
C           STORE COEFFICIENTS IN ELCO AND TESCO.
C           --------------------------------
            ELCO(1,NQ) = PINT*RQ1FAC
            ELCO(2,NQ) = 1.0D0
            DO 40 I = 2, NQ
               ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
   40       CONTINUE
            AGAMQ = RQFAC*XPIN
            RAGQ = 1.0D0/AGAMQ
            TESCO(2,NQ) = RAGQ
            IF (NQ .LT. 12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1
            TESCO(3,NQM1) = RAGQ
   50    CONTINUE
      GO TO 100
C
   60 CONTINUE
         PC(1) = 1.0D0
         RQ1FAC = 1.0D0
         DO 90 NQ = 1, 5
C           ------------------------------------------------------------
C            THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE
C                POLYNOMIAL P(X) = (X+1)*(X+2)*...*(X+NQ).
C            INITIALLY, P(X) = 1.
C           ------------------------------------------------------------
            FNQ = NQ
            NQP1 = NQ + 1
C           FORM COEFFICIENTS OF P(X)*(X+NQ).
C           ------------------------------------
            PC(NQP1) = 0.0D0
            DO 70 IB = 1, NQ
               I = NQ + 2 - IB
               PC(I) = PC(I-1) + FNQ*PC(I)
   70       CONTINUE
            PC(1) = FNQ*PC(1)
C           STORE COEFFICIENTS IN ELCO AND TESCO.
C           --------------------------------
            DO 80 I = 1, NQP1
               ELCO(I,NQ) = PC(I)/PC(2)
   80       CONTINUE
            ELCO(2,NQ) = 1.0D0
            TESCO(1,NQ) = RQ1FAC
            TESCO(2,NQ) = NQP1/ELCO(1,NQ)
            TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ)
            RQ1FAC = RQ1FAC/FNQ
   90    CONTINUE
  100 CONTINUE
      RETURN
C     ----------------------- END OF SUBROUTINE DCFOD
C     -----------------------
      END
