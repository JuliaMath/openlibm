*DECK DQELG
      SUBROUTINE DQELG (N, EPSTAB, RESULT, ABSERR, RES3LA, NRES)
C***BEGIN PROLOGUE  DQELG
C***SUBSIDIARY
C***PURPOSE  The routine determines the limit of a given sequence of
C            approximations, by means of the Epsilon algorithm of
C            P.Wynn. An estimate of the absolute error is also given.
C            The condensed Epsilon table is computed. Only those
C            elements needed for the computation of the next diagonal
C            are preserved.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (QELG-S, DQELG-D)
C***KEYWORDS  CONVERGENCE ACCELERATION, EPSILON ALGORITHM, EXTRAPOLATION
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Epsilon algorithm
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS
C              N      - Integer
C                       EPSTAB(N) contains the new element in the
C                       first column of the epsilon table.
C
C              EPSTAB - Double precision
C                       Vector of dimension 52 containing the elements
C                       of the two lower diagonals of the triangular
C                       epsilon table. The elements are numbered
C                       starting at the right-hand corner of the
C                       triangle.
C
C              RESULT - Double precision
C                       Resulting approximation to the integral
C
C              ABSERR - Double precision
C                       Estimate of the absolute error computed from
C                       RESULT and the 3 previous results
C
C              RES3LA - Double precision
C                       Vector of dimension 3 containing the last 3
C                       results
C
C              NRES   - Integer
C                       Number of calls to the routine
C                       (should be zero at first call)
C
C***SEE ALSO  DQAGIE, DQAGOE, DQAGPE, DQAGSE
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DQELG
C
      DOUBLE PRECISION ABSERR,DELTA1,DELTA2,DELTA3,D1MACH,
     1  EPMACH,EPSINF,EPSTAB,ERROR,ERR1,ERR2,ERR3,E0,E1,E1ABS,E2,E3,
     2  OFLOW,RES,RESULT,RES3LA,SS,TOL1,TOL2,TOL3
      INTEGER I,IB,IB2,IE,INDX,K1,K2,K3,LIMEXP,N,NEWELM,NRES,NUM
      DIMENSION EPSTAB(52),RES3LA(3)
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           E0     - THE 4 ELEMENTS ON WHICH THE COMPUTATION OF A NEW
C           E1       ELEMENT IN THE EPSILON TABLE IS BASED
C           E2
C           E3                 E0
C                        E3    E1    NEW
C                              E2
C           NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW
C                    DIAGONAL
C           ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2)
C           RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE
C                    OF ERROR
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
C           LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON
C           TABLE CAN CONTAIN. IF THIS NUMBER IS REACHED, THE UPPER
C           DIAGONAL OF THE EPSILON TABLE IS DELETED.
C
C***FIRST EXECUTABLE STATEMENT  DQELG
      EPMACH = D1MACH(4)
      OFLOW = D1MACH(2)
      NRES = NRES+1
      ABSERR = OFLOW
      RESULT = EPSTAB(N)
      IF(N.LT.3) GO TO 100
      LIMEXP = 50
      EPSTAB(N+2) = EPSTAB(N)
      NEWELM = (N-1)/2
      EPSTAB(N) = OFLOW
      NUM = N
      K1 = N
      DO 40 I = 1,NEWELM
        K2 = K1-1
        K3 = K1-2
        RES = EPSTAB(K1+2)
        E0 = EPSTAB(K3)
        E1 = EPSTAB(K2)
        E2 = RES
        E1ABS = ABS(E1)
        DELTA2 = E2-E1
        ERR2 = ABS(DELTA2)
        TOL2 = MAX(ABS(E2),E1ABS)*EPMACH
        DELTA3 = E1-E0
        ERR3 = ABS(DELTA3)
        TOL3 = MAX(E1ABS,ABS(E0))*EPMACH
        IF(ERR2.GT.TOL2.OR.ERR3.GT.TOL3) GO TO 10
C
C           IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE
C           ACCURACY, CONVERGENCE IS ASSUMED.
C           RESULT = E2
C           ABSERR = ABS(E1-E0)+ABS(E2-E1)
C
        RESULT = RES
        ABSERR = ERR2+ERR3
C ***JUMP OUT OF DO-LOOP
        GO TO 100
   10   E3 = EPSTAB(K1)
        EPSTAB(K1) = E1
        DELTA1 = E1-E3
        ERR1 = ABS(DELTA1)
        TOL1 = MAX(E1ABS,ABS(E3))*EPMACH
C
C           IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT
C           A PART OF THE TABLE BY ADJUSTING THE VALUE OF N
C
        IF(ERR1.LE.TOL1.OR.ERR2.LE.TOL2.OR.ERR3.LE.TOL3) GO TO 20
        SS = 0.1D+01/DELTA1+0.1D+01/DELTA2-0.1D+01/DELTA3
        EPSINF = ABS(SS*E1)
C
C           TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND
C           EVENTUALLY OMIT A PART OF THE TABLE ADJUSTING THE VALUE
C           OF N.
C
        IF(EPSINF.GT.0.1D-03) GO TO 30
   20   N = I+I-1
C ***JUMP OUT OF DO-LOOP
        GO TO 50
C
C           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST
C           THE VALUE OF RESULT.
C
   30   RES = E1+0.1D+01/SS
        EPSTAB(K1) = RES
        K1 = K1-2
        ERROR = ERR2+ABS(RES-E2)+ERR3
        IF(ERROR.GT.ABSERR) GO TO 40
        ABSERR = ERROR
        RESULT = RES
   40 CONTINUE
C
C           SHIFT THE TABLE.
C
   50 IF(N.EQ.LIMEXP) N = 2*(LIMEXP/2)-1
      IB = 1
      IF((NUM/2)*2.EQ.NUM) IB = 2
      IE = NEWELM+1
      DO 60 I=1,IE
        IB2 = IB+2
        EPSTAB(IB) = EPSTAB(IB2)
        IB = IB2
   60 CONTINUE
      IF(NUM.EQ.N) GO TO 80
      INDX = NUM-N+1
      DO 70 I = 1,N
        EPSTAB(I)= EPSTAB(INDX)
        INDX = INDX+1
   70 CONTINUE
   80 IF(NRES.GE.4) GO TO 90
      RES3LA(NRES) = RESULT
      ABSERR = OFLOW
      GO TO 100
C
C           COMPUTE ERROR ESTIMATE
C
   90 ABSERR = ABS(RESULT-RES3LA(3))+ABS(RESULT-RES3LA(2))
     1  +ABS(RESULT-RES3LA(1))
      RES3LA(1) = RES3LA(2)
      RES3LA(2) = RES3LA(3)
      RES3LA(3) = RESULT
  100 ABSERR = MAX(ABSERR,0.5D+01*EPMACH*ABS(RESULT))
      RETURN
      END
