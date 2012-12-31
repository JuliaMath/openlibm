*DECK DLSSUD
      SUBROUTINE DLSSUD (A, X, B, N, M, NRDA, U, NRDU, IFLAG, MLSO,
     +   IRANK, ISCALE, Q, DIAG, KPIVOT, S, DIV, TD, ISFLG, SCALES)
C***BEGIN PROLOGUE  DLSSUD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP and DSUDS
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (LSSUDS-S, DLSSUD-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C    DLSSUD solves the underdetermined system of equations  A Z = B,
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
C               To continue, simply reset IFLAG=1 and call DLSSUD again.
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
C***SEE ALSO  DBVSUP, DSUDS
C***REFERENCES  H. A. Watts, Solving linear least squares problems
C                 using SODS/SUDS/CODS, Sandia Report SAND77-0683,
C                 Sandia Laboratories, 1977.
C***ROUTINES CALLED  D1MACH, DDOT, DOHTRL, DORTHR, J4SAVE, XERMAX,
C                    XERMSG, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DLSSUD
      INTEGER J4SAVE
      DOUBLE PRECISION DDOT, D1MACH
      INTEGER I, IFLAG, IRANK, IRP, ISCALE, ISFLG, J, JR, K, KP,
     1     KPIVOT(*), L, M, MAXMES, MJ, MLSO, N, NFAT, NFATAL, NMIR,
     2     NRDA, NRDU, NU
      DOUBLE PRECISION A(NRDA,*), B(*), DIAG(*), DIV(*), GAM, GAMMA,
     1     Q(NRDA,*), RES, S(*), SCALES(*), SS, TD(*), U(NRDU,*), URO,
     2     X(*)
C
C     ******************************************************************
C
C          MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
C          BY THE FUNCTION D1MACH.
C
C     ******************************************************************
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 310
C        BEGIN BLOCK PERMITTING ...EXITS TO 80
C***FIRST EXECUTABLE STATEMENT  DLSSUD
            URO = D1MACH(4)
C
            IF (N .LT. 1 .OR. M .LT. N .OR. NRDA .LT. N) GO TO 70
            IF (NRDU .NE. 0 .AND. NRDU .LT. M) GO TO 70
               IF (IFLAG .GT. 0) GO TO 60
C
                  CALL XGETF(NFATAL)
                  MAXMES = J4SAVE(4,0,.FALSE.)
                  ISFLG = -15
                  IF (IFLAG .EQ. 0) GO TO 10
                     ISFLG = IFLAG
                     NFAT = -1
                     IF (NFATAL .EQ. 0) NFAT = 0
                     CALL XSETF(NFAT)
                     CALL XERMAX(1)
   10             CONTINUE
C
C                 COPY MATRIX A INTO MATRIX Q
C
                  DO 30 K = 1, M
                     DO 20 J = 1, N
                        Q(J,K) = A(J,K)
   20                CONTINUE
   30             CONTINUE
C
C                 USE ORTHOGONAL TRANSFORMATIONS TO REDUCE Q TO LOWER
C                 TRIANGULAR FORM
C
                  CALL DORTHR(Q,N,M,NRDA,IFLAG,IRANK,ISCALE,DIAG,KPIVOT,
     1                        SCALES,DIV,TD)
C
                  CALL XSETF(NFATAL)
                  CALL XERMAX(MAXMES)
                  IF (IRANK .EQ. N) GO TO 40
C
C                    FOR RANK DEFICIENT PROBLEMS USE ADDITIONAL
C                    ORTHOGONAL TRANSFORMATIONS TO FURTHER REDUCE Q
C
                     IF (IRANK .NE. 0)
     1                  CALL DOHTRL(Q,N,NRDA,DIAG,IRANK,DIV,TD)
C     ...............EXIT
                     GO TO 310
   40             CONTINUE
C
C                 STORE DIVISORS FOR THE TRIANGULAR SOLUTION
C
                  DO 50 K = 1, N
                     DIV(K) = DIAG(K)
   50             CONTINUE
C        .........EXIT
                  GO TO 80
   60          CONTINUE
C        ......EXIT
               IF (IFLAG .EQ. 1) GO TO 80
   70       CONTINUE
C
C           INVALID INPUT FOR DLSSUD
            IFLAG = 2
            CALL XERMSG ('SLATEC', 'DLSSUD',
     +         'INVALID IMPUT PARAMETERS.', 2, 1)
C     ......EXIT
            GO TO 310
   80    CONTINUE
C
C
         IF (IRANK .GT. 0) GO TO 130
C
C           SPECIAL CASE FOR THE NULL MATRIX
            DO 110 K = 1, M
               X(K) = 0.0D0
               IF (MLSO .EQ. 0) GO TO 100
                  U(K,K) = 1.0D0
                  DO 90 J = 1, M
                     IF (J .NE. K) U(J,K) = 0.0D0
   90             CONTINUE
  100          CONTINUE
  110       CONTINUE
            DO 120 K = 1, N
               IF (B(K) .GT. 0.0D0) IFLAG = 4
  120       CONTINUE
         GO TO 300
  130    CONTINUE
C           BEGIN BLOCK PERMITTING ...EXITS TO 180
C
C              COPY CONSTANT VECTOR INTO S AFTER FIRST INTERCHANGING
C              THE ELEMENTS ACCORDING TO THE PIVOTAL SEQUENCE
C
               DO 140 K = 1, N
                  KP = KPIVOT(K)
                  X(K) = B(KP)
  140          CONTINUE
               DO 150 K = 1, N
                  S(K) = X(K)
  150          CONTINUE
C
               IRP = IRANK + 1
               NU = 1
               IF (MLSO .EQ. 0) NU = 0
C           ...EXIT
               IF (IRANK .EQ. N) GO TO 180
C
C              FOR RANK DEFICIENT PROBLEMS WE MUST APPLY THE
C              ORTHOGONAL TRANSFORMATION TO S
C              WE ALSO CHECK TO SEE IF THE SYSTEM APPEARS TO BE
C              INCONSISTENT
C
               NMIR = N - IRANK
               SS = DDOT(N,S(1),1,S(1),1)
               DO 170 L = 1, IRANK
                  K = IRP - L
                  GAM = ((TD(K)*S(K)) + DDOT(NMIR,Q(IRP,K),1,S(IRP),1))
     1                  /(TD(K)*DIV(K))
                  S(K) = S(K) + GAM*TD(K)
                  DO 160 J = IRP, N
                     S(J) = S(J) + GAM*Q(J,K)
  160             CONTINUE
  170          CONTINUE
               RES = DDOT(NMIR,S(IRP),1,S(IRP),1)
C           ...EXIT
               IF (RES
     1             .LE. SS*(10.0D0*MAX(10.0D0**ISFLG,10.0D0*URO))**2)
     2            GO TO 180
C
C              INCONSISTENT SYSTEM
               IFLAG = 4
               NU = 0
  180       CONTINUE
C
C           APPLY FORWARD SUBSTITUTION TO SOLVE LOWER TRIANGULAR SYSTEM
C
            S(1) = S(1)/DIV(1)
            IF (IRANK .LT. 2) GO TO 200
            DO 190 K = 2, IRANK
               S(K) = (S(K) - DDOT(K-1,Q(K,1),NRDA,S(1),1))/DIV(K)
  190       CONTINUE
  200       CONTINUE
C
C           INITIALIZE X VECTOR AND THEN APPLY ORTHOGONAL TRANSFORMATION
C
            DO 210 K = 1, M
               X(K) = 0.0D0
               IF (K .LE. IRANK) X(K) = S(K)
  210       CONTINUE
C
            DO 230 JR = 1, IRANK
               J = IRP - JR
               MJ = M - J + 1
               GAMMA = DDOT(MJ,Q(J,J),NRDA,X(J),1)/(DIAG(J)*Q(J,J))
               DO 220 K = J, M
                  X(K) = X(K) + GAMMA*Q(J,K)
  220          CONTINUE
  230       CONTINUE
C
C           RESCALE ANSWERS AS DICTATED
C
            DO 240 K = 1, M
               X(K) = X(K)*SCALES(K)
  240       CONTINUE
C
            IF (NU .EQ. 0 .OR. M .EQ. IRANK) GO TO 290
C
C              INITIALIZE U MATRIX AND THEN APPLY ORTHOGONAL
C              TRANSFORMATION
C
               L = M - IRANK
               DO 280 K = 1, L
                  DO 250 I = 1, M
                     U(I,K) = 0.0D0
                     IF (I .EQ. IRANK + K) U(I,K) = 1.0D0
  250             CONTINUE
C
                  DO 270 JR = 1, IRANK
                     J = IRP - JR
                     MJ = M - J + 1
                     GAMMA = DDOT(MJ,Q(J,J),NRDA,U(J,K),1)
     1                       /(DIAG(J)*Q(J,J))
                     DO 260 I = J, M
                        U(I,K) = U(I,K) + GAMMA*Q(J,I)
  260                CONTINUE
  270             CONTINUE
  280          CONTINUE
  290       CONTINUE
  300    CONTINUE
  310 CONTINUE
C
      RETURN
      END
