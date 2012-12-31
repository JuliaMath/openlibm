*DECK LSSUDS
      SUBROUTINE LSSUDS (A, X, B, N, M, NRDA, U, NRDU, IFLAG, MLSO,
     +   IRANK, ISCALE, Q, DIAG, KPIVOT, S, DIV, TD, ISFLG, SCALES)
C***BEGIN PROLOGUE  LSSUDS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (LSSUDS-S, DLSSUD-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C    LSSUDS solves the underdetermined system of equations  A Z = B,
C    where A is N by M and N .LE. M.  In particular, if rank A equals
C    IRA, a vector X and a matrix U are determined such that X is the
C    UNIQUE solution of smallest length, satisfying A X = B, and the
C    columns of U form an orthonormal basis for the null space of A,
C    satisfying A U = 0 .  Then all solutions Z are given by
C              Z = X + C(1)*U(1) + ..... + C(M-IRA)*U(M-IRA)
C    where U(J) represents the J-th column of U and the C(J) are
C    arbitrary constants.
C    If the system of equations are not compatible, only the least
C    squares solution of minimal length is computed.
C
C *********************************************************************
C   INPUT
C *********************************************************************
C
C     A -- Contains the matrix of N equations in M unknowns, A remains
C          unchanged, must be dimensioned NRDA by M.
C     X -- Solution array of length at least M.
C     B -- Given constant vector of length N, B remains unchanged.
C     N -- Number of equations, N greater or equal to 1.
C     M -- Number of unknowns, M greater or equal to N.
C     NRDA -- Row dimension of A, NRDA greater or equal to N.
C     U -- Matrix used for solution, must be dimensioned NRDU by
C          (M - rank of A).
C          (storage for U may be ignored when only the minimal length
C           solution X is desired)
C     NRDU -- Row dimension of U, NRDU greater or equal to M.
C             (if only the minimal length solution is wanted,
C              NRDU=0 is acceptable)
C     IFLAG -- Status indicator
C           =0  for the first call (and for each new problem defined by
C               a new matrix A) when the matrix data is treated as exact
C           =-K for the first call (and for each new problem defined by
C               a new matrix A) when the matrix data is assumed to be
C               accurate to about K digits.
C           =1  for subsequent calls whenever the matrix A has already
C               been decomposed (problems with new vectors B but
C               same matrix A can be handled efficiently).
C     MLSO -- =0 if only the minimal length solution is wanted.
C             =1 if the complete solution is wanted, includes the
C                linear space defined by the matrix U.
C     IRANK -- Variable used for the rank of A, set by the code.
C     ISCALE -- Scaling indicator
C               =-1 if the matrix A is to be pre-scaled by
C               columns when appropriate.
C               If the scaling indicator is not equal to -1
C               no scaling will be attempted.
C            For most problems scaling will probably not be necessary.
C     Q -- Matrix used for the transformation, must be dimensioned
C            NRDA by M.
C     DIAG,KPIVOT,S, -- Arrays of length at least N used for internal
C      DIV,TD,SCALES    storage (except for SCALES which is M).
C     ISFLG -- Storage for an internal variable.
C
C *********************************************************************
C   OUTPUT
C *********************************************************************
C
C     IFLAG -- Status indicator
C            =1 if solution was obtained.
C            =2 if improper input is detected.
C            =3 if rank of matrix is less than N.
C               To continue, simply reset IFLAG=1 and call LSSUDS again.
C            =4 if the system of equations appears to be inconsistent.
C               However, the least squares solution of minimal length
C               was obtained.
C     X -- Minimal length least squares solution of A Z = B
C     IRANK -- Numerically determined rank of A, must not be altered
C              on succeeding calls with input values of IFLAG=1.
C     U -- Matrix whose M-IRANK columns are mutually orthogonal unit
C          vectors which span the null space of A. This is to be ignored
C          when MLSO was set to zero or IFLAG=4 on output.
C     Q -- Contains the strictly upper triangular part of the reduced
C           matrix and transformation information.
C     DIAG -- Contains the diagonal elements of the triangular reduced
C             matrix.
C     KPIVOT -- Contains the pivotal information.  The row interchanges
C               performed on the original matrix are recorded here.
C     S -- Contains the solution of the lower triangular system.
C     DIV,TD -- Contains transformation information for rank
C               deficient problems.
C     SCALES -- Contains the column scaling parameters.
C
C *********************************************************************
C
C***SEE ALSO  BVSUP
C***REFERENCES  H. A. Watts, Solving linear least squares problems
C                 using SODS/SUDS/CODS, Sandia Report SAND77-0683,
C                 Sandia Laboratories, 1977.
C***ROUTINES CALLED  J4SAVE, OHTROL, ORTHOR, R1MACH, SDOT, XERMAX,
C                    XERMSG, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C   900510  Fixed an error message.  (RWC)
C   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  LSSUDS
      DIMENSION A(NRDA,*),X(*),B(*),U(NRDU,*),Q(NRDA,*),
     1          DIAG(*),KPIVOT(*),S(*),DIV(*),TD(*),SCALES(*)
C
C **********************************************************************
C
C     MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
C     BY THE FUNCTION R1MACH.
C
C***FIRST EXECUTABLE STATEMENT  LSSUDS
      URO = R1MACH(4)
C
C **********************************************************************
C
      IF (N .LT. 1  .OR.  M .LT. N  .OR.  NRDA .LT. N) GO TO 1
      IF (NRDU .NE. 0  .AND.  NRDU .LT. M) GO TO 1
      IF (IFLAG .LE. 0) GO TO 5
      IF (IFLAG .EQ. 1) GO TO 25
C
C     INVALID INPUT FOR LSSUDS
    1 IFLAG=2
      CALL XERMSG ('SLATEC', 'LSSUDS', 'INVALID INPUT PARAMETERS.', 2,
     +   1)
      RETURN
C
    5 CALL XGETF(NFATAL)
      MAXMES = J4SAVE (4,0,.FALSE.)
      ISFLG=-15
      IF (IFLAG .EQ. 0) GO TO 7
      ISFLG=IFLAG
      NFAT = -1
      IF (NFATAL .EQ. 0) NFAT=0
      CALL XSETF(NFAT)
      CALL XERMAX(1)
C
C     COPY MATRIX A INTO MATRIX Q
C
    7 DO 10 K=1,M
         DO 10 J=1,N
   10       Q(J,K)=A(J,K)
C
C     USE ORTHOGONAL TRANSFORMATIONS TO REDUCE Q TO LOWER
C     TRIANGULAR FORM
C
      CALL ORTHOR(Q,N,M,NRDA,IFLAG,IRANK,ISCALE,DIAG,KPIVOT,SCALES,
     1            DIV,TD)
C
      CALL XSETF(NFATAL)
      CALL XERMAX(MAXMES)
      IF (IRANK .EQ. N) GO TO 15
C
C     FOR RANK DEFICIENT PROBLEMS USE ADDITIONAL ORTHOGONAL
C     TRANSFORMATIONS TO FURTHER REDUCE Q
C
      IF (IRANK .NE. 0) CALL OHTROL(Q,N,NRDA,DIAG,IRANK,DIV,TD)
      RETURN
C
C     STORE DIVISORS FOR THE TRIANGULAR SOLUTION
C
   15 DO 20 K=1,N
   20    DIV(K)=DIAG(K)
C
C
   25 IF (IRANK .GT. 0) GO TO 40
C
C     SPECIAL CASE FOR THE NULL MATRIX
      DO 35 K=1,M
         X(K)=0.
         IF (MLSO .EQ. 0) GO TO 35
         U(K,K)=1.
         DO 30 J=1,M
            IF (J .EQ. K) GO TO 30
            U(J,K)=0.
   30    CONTINUE
   35 CONTINUE
      DO 37 K=1,N
         IF (B(K) .GT. 0.) IFLAG=4
   37 CONTINUE
      RETURN
C
C     COPY CONSTANT VECTOR INTO S AFTER FIRST INTERCHANGING
C     THE ELEMENTS ACCORDING TO THE PIVOTAL SEQUENCE
C
   40 DO 45 K=1,N
         KP=KPIVOT(K)
   45    X(K)=B(KP)
      DO 50 K=1,N
   50    S(K)=X(K)
C
      IRP=IRANK+1
      NU=1
      IF (MLSO .EQ. 0) NU=0
      IF (IRANK .EQ. N) GO TO 60
C
C     FOR RANK DEFICIENT PROBLEMS WE MUST APPLY THE
C     ORTHOGONAL TRANSFORMATION TO S
C     WE ALSO CHECK TO SEE IF THE SYSTEM APPEARS TO BE INCONSISTENT
C
      NMIR=N-IRANK
      SS=SDOT(N,S(1),1,S(1),1)
      DO 55 L=1,IRANK
         K=IRP-L
         GAM=((TD(K)*S(K))+SDOT(NMIR,Q(IRP,K),1,S(IRP),1))/
     1             (TD(K)*DIV(K))
         S(K)=S(K)+GAM*TD(K)
         DO 55 J=IRP,N
   55       S(J)=S(J)+GAM*Q(J,K)
      RES=SDOT(NMIR,S(IRP),1,S(IRP),1)
      IF (RES .LE. SS*(10.*MAX(10.**ISFLG,10.*URO))**2) GO TO 60
C
C     INCONSISTENT SYSTEM
      IFLAG=4
      NU=0
C
C     APPLY FORWARD SUBSTITUTION TO SOLVE LOWER TRIANGULAR SYSTEM
C
   60 S(1)=S(1)/DIV(1)
      IF (IRANK .EQ. 1) GO TO 70
      DO 65 K=2,IRANK
   65    S(K)=(S(K)-SDOT(K-1,Q(K,1),NRDA,S(1),1))/DIV(K)
C
C     INITIALIZE X VECTOR AND THEN APPLY ORTHOGONAL TRANSFORMATION
C
   70 DO 75 K=1,M
         X(K)=0.
         IF (K .LE. IRANK) X(K)=S(K)
   75 CONTINUE
C
      DO 80 JR=1,IRANK
         J=IRP-JR
         MJ=M-J+1
         GAMMA=SDOT(MJ,Q(J,J),NRDA,X(J),1)/(DIAG(J)*Q(J,J))
         DO 80 K=J,M
   80       X(K)=X(K)+GAMMA*Q(J,K)
C
C     RESCALE ANSWERS AS DICTATED
C
      DO 85 K=1,M
   85    X(K)=X(K)*SCALES(K)
C
      IF ((NU .EQ. 0) .OR. (M .EQ. IRANK)) RETURN
C
C     INITIALIZE U MATRIX AND THEN APPLY ORTHOGONAL TRANSFORMATION
C
      L=M-IRANK
      DO 100 K=1,L
         DO 90 I=1,M
            U(I,K)=0.
            IF (I .EQ. IRANK+K) U(I,K)=1.
   90    CONTINUE
C
         DO 100 JR=1,IRANK
            J=IRP-JR
            MJ=M-J+1
            GAMMA=SDOT(MJ,Q(J,J),NRDA,U(J,K),1)/(DIAG(J)*Q(J,J))
            DO 100 I=J,M
  100          U(I,K)=U(I,K)+GAMMA*Q(J,I)
C
      RETURN
      END
