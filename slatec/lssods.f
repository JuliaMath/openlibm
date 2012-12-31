*DECK LSSODS
      SUBROUTINE LSSODS (A, X, B, M, N, NRDA, IFLAG, IRANK, ISCALE, Q,
     +   DIAG, KPIVOT, ITER, RESNRM, XNORM, Z, R, DIV, TD, SCALES)
C***BEGIN PROLOGUE  LSSODS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (LSSODS-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     LSSODS solves the same problem as SODS (in fact, it is called by
C     SODS) but is somewhat more flexible in its use. In particular,
C     LSSODS allows for iterative refinement of the solution, makes the
C     transformation and triangular reduction information more
C     accessible, and enables the user to avoid destruction of the
C     original matrix A.
C
C     Modeled after the ALGOL codes in the articles in the REFERENCES
C     section.
C
C **********************************************************************
C   INPUT
C **********************************************************************
C
C     A -- Contains the matrix of M equations in N unknowns and must
C          be dimensioned NRDA by N. A remains unchanged
C     X -- Solution array of length at least N
C     B -- Given constant vector of length M, B remains unchanged
C     M -- Number of equations, M greater or equal to 1
C     N -- Number of unknowns, N not larger than M
C  NRDA -- Row dimension of A, NRDA greater or equal to M
C IFLAG -- Status indicator
C         = 0 for the first call (and for each new problem defined by
C             a new matrix A) when the matrix data is treated as exact
C         =-K for the first call (and for each new problem defined by
C             a new matrix A) when the matrix data is assumed to be
C             accurate to about K digits
C         = 1 for subsequent calls whenever the matrix A has already
C             been decomposed (problems with new vectors B but
C             same matrix a can be handled efficiently)
C ISCALE -- Scaling indicator
C         =-1 if the matrix A is to be pre-scaled by
C             columns when appropriate
C             If the scaling indicator is not equal to -1
C             no scaling will be attempted
C             For most problems scaling will probably not be necessary
C   ITER -- Maximum number of iterative improvement steps to be
C           performed,  0 .LE. ITER .LE. 10   (SODS uses ITER=0)
C      Q -- Matrix used for the transformation, must be dimensioned
C           NRDA by N  (SODS puts A in the Q location which conserves
C           storage but destroys A)
C           When iterative improvement of the solution is requested,
C           ITER .GT. 0, this additional storage for Q must be
C           made available
C DIAG,KPIVOT,Z,R, -- Arrays of length N (except for R which is M)
C   DIV,TD,SCALES     used for internal storage
C
C **********************************************************************
C   OUTPUT
C **********************************************************************
C
C  IFLAG -- Status indicator
C            =1 if solution was obtained
C            =2 if improper input is detected
C            =3 if rank of matrix is less than N
C               if the minimal length least squares solution is
C               desired, simply reset IFLAG=1 and call the code again
C
C       The next three IFLAG values can occur only when
C        the iterative improvement mode is being used.
C            =4 if the problem is ill-conditioned and maximal
C               machine accuracy is not achievable
C            =5 if the problem is very ill-conditioned and the solution
C               IS likely to have no correct digits
C            =6 if the allowable number of iterative improvement steps
C               has been completed without getting convergence
C      X -- Least squares solution of  A X = B
C  IRANK -- Contains the numerically determined matrix rank
C           the user must not alter this value on succeeding calls
C           with input values of IFLAG=1
C      Q -- Contains the strictly upper triangular part of the reduced
C           matrix and the transformation information in the lower
C           triangular part
C   DIAG -- Contains the diagonal elements of the triangular reduced
C           matrix
C KPIVOT -- Contains the pivotal information.  The column interchanges
C           performed on the original matrix are recorded here
C   ITER -- The actual number of iterative corrections used
C RESNRM -- The Euclidean norm of the residual vector  B - A X
C  XNORM -- The Euclidean norm of the solution vector
C DIV,TD -- Contains transformation information for rank
C           deficient problems
C SCALES -- Contains the column scaling parameters
C
C **********************************************************************
C
C***SEE ALSO  BVSUP
C***REFERENCES  G. Golub, Numerical methods for solving linear least
C                 squares problems, Numerische Mathematik 7, (1965),
C                 pp. 206-216.
C               P. Businger and G. Golub, Linear least squares
C                 solutions by Householder transformations, Numerische
C                 Mathematik  7, (1965), pp. 269-276.
C***ROUTINES CALLED  J4SAVE, OHTROR, ORTHOL, R1MACH, SDOT, SDSDOT,
C                    XERMAX, XERMSG, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900402  Added TYPE section.  (WRB)
C   910408  Updated the REFERENCES section.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  LSSODS
      DIMENSION A(NRDA,*),X(*),B(*),Q(NRDA,*),DIAG(*),
     1          Z(*),KPIVOT(*),R(*),DIV(*),TD(*),SCALES(*)
C
C **********************************************************************
C
C     MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
C     THE FUNCTION R1MACH.
C
C***FIRST EXECUTABLE STATEMENT  LSSODS
      URO = R1MACH(3)
C
C **********************************************************************
C
      IF (N .LT. 1  .OR.  M .LT. N  .OR.  NRDA .LT. M) GO TO 1
      IF (ITER .LT. 0) GO TO 1
      IF (IFLAG .LE. 0) GO TO 5
      IF (IFLAG .EQ. 1) GO TO 15
C
C     INVALID INPUT FOR LSSODS
    1 IFLAG=2
      CALL XERMSG ('SLATEC', 'LSSODS', 'INVALID INPUT PARAMETERS.', 2,
     +   1)
      RETURN
C
    5 CALL XGETF (NFATAL)
      MAXMES = J4SAVE (4,0,.FALSE.)
      IF (IFLAG .EQ. 0) GO TO 7
      NFAT = -1
      IF(NFATAL .EQ. 0) NFAT=0
      CALL XSETF (NFAT)
      CALL XERMAX (1)
C
C     COPY MATRIX A INTO MATRIX Q
C
    7 DO 10 J=1,N
         DO 10 K=1,M
   10       Q(K,J)=A(K,J)
C
C     USE ORTHOGONAL TRANSFORMATIONS TO REDUCE Q TO
C     UPPER TRIANGULAR FORM
C
      CALL ORTHOL(Q,M,N,NRDA,IFLAG,IRANK,ISCALE,DIAG,KPIVOT,SCALES,Z,TD)
C
      CALL XSETF (NFATAL)
      CALL XERMAX (MAXMES)
      IF (IRANK .EQ. N) GO TO 12
C
C     FOR RANK DEFICIENT PROBLEMS USE ADDITIONAL ORTHOGONAL
C     TRANSFORMATIONS TO FURTHER REDUCE Q
C
      IF (IRANK .NE. 0) CALL OHTROR(Q,N,NRDA,DIAG,IRANK,DIV,TD)
      RETURN
C
C     STORE DIVISORS FOR THE TRIANGULAR SOLUTION
C
   12 DO 13 K=1,N
   13    DIV(K)=DIAG(K)
C
   15 IRM=IRANK-1
      IRP=IRANK+1
      ITERP=MIN(ITER+1,11)
      ACC=10.*URO
C
C     ZERO OUT SOLUTION ARRAY
C
      DO 20 K=1,N
   20    X(K)=0.
C
      IF (IRANK .GT. 0) GO TO 25
C
C     SPECIAL CASE FOR THE NULL MATRIX
      ITER=0
      XNORM=0.
      RESNRM=SQRT(SDOT(M,B(1),1,B(1),1))
      RETURN
C
C     COPY CONSTANT VECTOR INTO R
C
   25 DO 30 K=1,M
   30    R(K)=B(K)
C
C **********************************************************************
C     SOLUTION SECTION
C     ITERATIVE REFINEMENT OF THE RESIDUAL VECTOR
C **********************************************************************
C
      DO 100 IT=1,ITERP
         ITER=IT-1
C
C        APPLY ORTHOGONAL TRANSFORMATION TO R
C
         DO 35 J=1,IRANK
            MJ=M-J+1
            GAMMA=SDOT(MJ,Q(J,J),1,R(J),1)/(DIAG(J)*Q(J,J))
            DO 35 K=J,M
   35          R(K)=R(K)+GAMMA*Q(K,J)
C
C        BACKWARD SUBSTITUTION FOR TRIANGULAR SYSTEM SOLUTION
C
         Z(IRANK)=R(IRANK)/DIV(IRANK)
         IF (IRM .EQ. 0) GO TO 45
         DO 40 L=1,IRM
            K=IRANK-L
            KP=K+1
   40       Z(K)=(R(K)-SDOT(L,Q(K,KP),NRDA,Z(KP),1))/DIV(K)
C
   45    IF (IRANK .EQ. N) GO TO 60
C
C        FOR RANK DEFICIENT PROBLEMS OBTAIN THE
C        MINIMAL LENGTH SOLUTION
C
         NMIR=N-IRANK
         DO 50 K=IRP,N
   50       Z(K)=0.
         DO 55 K=1,IRANK
            GAM=((TD(K)*Z(K))+SDOT(NMIR,Q(K,IRP),NRDA,Z(IRP),1))/
     1                (TD(K)*DIV(K))
            Z(K)=Z(K)+GAM*TD(K)
            DO 55 J=IRP,N
   55          Z(J)=Z(J)+GAM*Q(K,J)
C
C        REORDER SOLUTION COMPONENTS ACCORDING TO PIVOTAL POINTS
C        AND RESCALE ANSWERS AS DICTATED
C
   60    DO 65 K=1,N
            Z(K)=Z(K)*SCALES(K)
            L=KPIVOT(K)
   65       X(L)=X(L)+Z(K)
C
C        COMPUTE CORRECTION VECTOR NORM (SOLUTION NORM)
C
         ZNORM=SQRT(SDOT(N,Z(1),1,Z(1),1))
         IF (IT .EQ. 1) XNORM=ZNORM
         IF (ITERP .GT. 1) GO TO 80
C
C        NO ITERATIVE CORRECTIONS TO BE PERFORMED, SO COMPUTE
C        THE APPROXIMATE RESIDUAL NORM DEFINED BY THE EQUATIONS
C        WHICH ARE NOT SATISFIED BY THE SOLUTION
C        THEN WE ARE DONE
C
         MMIR=M-IRANK
         IF (MMIR .EQ. 0) GO TO 70
         RESNRM=SQRT(SDOT(MMIR,R(IRP),1,R(IRP),1))
         RETURN
   70    RESNRM=0.
         RETURN
C
C        COMPUTE RESIDUAL VECTOR FOR THE ITERATIVE IMPROVEMENT PROCESS
C
   80    DO 85 K=1,M
   85       R(K)=-SDSDOT(N,-B(K),A(K,1),NRDA,X(1),1)
         RESNRM=SQRT(SDOT(M,R(1),1,R(1),1))
         IF (IT .EQ. 1) GO TO 100
C
C        TEST FOR CONVERGENCE
C
         IF (ZNORM .LE. ACC*XNORM) RETURN
C
C        COMPARE SUCCESSIVE REFINEMENT VECTOR NORMS
C        FOR LOOP TERMINATION CRITERIA
C
         IF (ZNORM .LE. 0.25*ZNRM0) GO TO 100
         IF (IT .EQ. 2) GO TO 90
C
         IFLAG=4
         CALL XERMSG ('SLATEC', 'LSSODS',
     +   'PROBLEM MAY BE ILL-CONDITIONED.  MAXIMAL MACHINE ACCURACY ' //
     +   'IS NOT ACHIEVABLE.', 3, 1)
         RETURN
C
   90    IFLAG=5
         CALL XERMSG ('SLATEC', 'LSSODS',
     +      'PROBLEM IS VERY ILL-CONDITIONED.  ITERATIVE ' //
     +      'IMPROVEMENT IS INEFFECTIVE.', 8, 1)
         RETURN
C
  100    ZNRM0=ZNORM
C **********************************************************************
C
C **********************************************************************
      IFLAG=6
         CALL XERMSG ('SLATEC', 'LSSODS',
     +      'CONVERGENCE HAS NOT BEEN OBTAINED WITH ALLOWABLE ' //
     +      'NUMBER OF ITERATIVE IMPROVEMENT STEPS.', 8, 1)
C
      RETURN
      END
