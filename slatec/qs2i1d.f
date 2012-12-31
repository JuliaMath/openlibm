*DECK QS2I1D
      SUBROUTINE QS2I1D (IA, JA, A, N, KFLAG)
C***BEGIN PROLOGUE  QS2I1D
C***SUBSIDIARY
C***PURPOSE  Sort an integer array, moving an integer and DP array.
C            This routine sorts the integer array IA and makes the same
C            interchanges in the integer array JA and the double pre-
C            cision array A.  The array IA may be sorted in increasing
C            order or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  N6A2A
C***TYPE      DOUBLE PRECISION (QS2I1R-S, QS2I1D-D)
C***KEYWORDS  SINGLETON QUICKSORT, SLAP, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Kahaner, D. K., (NBS)
C           Seager, M. K., (LLNL) seager@llnl.gov
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C     Written by Rondall E Jones
C     Modified by John A. Wisniewski to use the Singleton QUICKSORT
C     algorithm. date 18 November 1976.
C
C     Further modified by David K. Kahaner
C     National Bureau of Standards
C     August, 1981
C
C     Even further modification made to bring the code up to the
C     Fortran 77 level and make it more readable and to carry
C     along one integer array and one double precision array during
C     the sort by
C     Mark K. Seager
C     Lawrence Livermore National Laboratory
C     November, 1987
C     This routine was adapted from the ISORT routine.
C
C     ABSTRACT
C         This routine sorts an integer array IA and makes the same
C         interchanges in the integer array JA and the double precision
C         array A.
C         The array IA may be sorted in increasing order or decreasing
C         order.  A slightly modified quicksort algorithm is used.
C
C     DESCRIPTION OF PARAMETERS
C        IA - Integer array of values to be sorted.
C        JA - Integer array to be carried along.
C         A - Double Precision array to be carried along.
C         N - Number of values in integer array IA to be sorted.
C     KFLAG - Control parameter
C           = 1 means sort IA in INCREASING order.
C           =-1 means sort IA in DECREASING order.
C
C***SEE ALSO  DS2Y
C***REFERENCES  R. C. Singleton, Algorithm 347, An Efficient Algorithm
C                 for Sorting With Minimal Storage, Communications ACM
C                 12:3 (1969), pp.185-7.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   761118  DATE WRITTEN
C   890125  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   900805  Changed XERROR calls to calls to XERMSG.  (RWC)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   910506  Made subsidiary to DS2Y and corrected reference.  (FNF)
C   920511  Added complete declaration section.  (WRB)
C   920929  Corrected format of reference.  (FNF)
C   921012  Corrected all f.p. constants to double precision.  (FNF)
C***END PROLOGUE  QS2I1D
CVD$R NOVECTOR
CVD$R NOCONCUR
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      DOUBLE PRECISION A(N)
      INTEGER IA(N), JA(N)
C     .. Local Scalars ..
      DOUBLE PRECISION R, TA, TTA
      INTEGER I, IIT, IJ, IT, J, JJT, JT, K, KK, L, M, NN
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. External Subroutines ..
      EXTERNAL XERMSG
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***FIRST EXECUTABLE STATEMENT  QS2I1D
      NN = N
      IF (NN.LT.1) THEN
         CALL XERMSG ('SLATEC', 'QS2I1D',
     $      'The number of values to be sorted was not positive.', 1, 1)
         RETURN
      ENDIF
      IF( N.EQ.1 ) RETURN
      KK = ABS(KFLAG)
      IF ( KK.NE.1 ) THEN
         CALL XERMSG ('SLATEC', 'QS2I1D',
     $      'The sort control parameter, K, was not 1 or -1.', 2, 1)
         RETURN
      ENDIF
C
C     Alter array IA to get decreasing order if needed.
C
      IF( KFLAG.LT.1 ) THEN
         DO 20 I=1,NN
            IA(I) = -IA(I)
 20      CONTINUE
      ENDIF
C
C     Sort IA and carry JA and A along.
C     And now...Just a little black magic...
      M = 1
      I = 1
      J = NN
      R = .375D0
 210  IF( R.LE.0.5898437D0 ) THEN
         R = R + 3.90625D-2
      ELSE
         R = R-.21875D0
      ENDIF
 225  K = I
C
C     Select a central element of the array and save it in location
C     it, jt, at.
C
      IJ = I + INT ((J-I)*R)
      IT = IA(IJ)
      JT = JA(IJ)
      TA = A(IJ)
C
C     If first element of array is greater than it, interchange with it.
C
      IF( IA(I).GT.IT ) THEN
         IA(IJ) = IA(I)
         IA(I)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(I)
         JA(I)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(I)
         A(I)   = TA
         TA     = A(IJ)
      ENDIF
      L=J
C
C     If last element of array is less than it, swap with it.
C
      IF( IA(J).LT.IT ) THEN
         IA(IJ) = IA(J)
         IA(J)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(J)
         JA(J)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(J)
         A(J)   = TA
         TA     = A(IJ)
C
C     If first element of array is greater than it, swap with it.
C
         IF ( IA(I).GT.IT ) THEN
            IA(IJ) = IA(I)
            IA(I)  = IT
            IT     = IA(IJ)
            JA(IJ) = JA(I)
            JA(I)  = JT
            JT     = JA(IJ)
            A(IJ)  = A(I)
            A(I)   = TA
            TA     = A(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is
C     smaller than it.
C
  240 L=L-1
      IF( IA(L).GT.IT ) GO TO 240
C
C     Find an element in the first half of the array which is
C     greater than it.
C
  245 K=K+1
      IF( IA(K).LT.IT ) GO TO 245
C
C     Interchange these elements.
C
      IF( K.LE.L ) THEN
         IIT   = IA(L)
         IA(L) = IA(K)
         IA(K) = IIT
         JJT   = JA(L)
         JA(L) = JA(K)
         JA(K) = JJT
         TTA   = A(L)
         A(L)  = A(K)
         A(K)  = TTA
         GOTO 240
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted.
C
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
C
C     Begin again on another portion of the unsorted array.
C
  255 M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
  260 IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
  265 I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IA(I+1)
      JT = JA(I+1)
      TA =  A(I+1)
      IF( IA(I).LE.IT ) GO TO 265
      K=I
  270 IA(K+1) = IA(K)
      JA(K+1) = JA(K)
      A(K+1)  =  A(K)
      K = K-1
      IF( IT.LT.IA(K) ) GO TO 270
      IA(K+1) = IT
      JA(K+1) = JT
      A(K+1)  = TA
      GO TO 265
C
C     Clean up, if necessary.
C
  300 IF( KFLAG.LT.1 ) THEN
         DO 310 I=1,NN
            IA(I) = -IA(I)
 310     CONTINUE
      ENDIF
      RETURN
C------------- LAST LINE OF QS2I1D FOLLOWS ----------------------------
      END
