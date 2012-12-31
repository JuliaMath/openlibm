*DECK WNLSM
      SUBROUTINE WNLSM (W, MDW, MME, MA, N, L, PRGOPT, X, RNORM, MODE,
     +   IPIVOT, ITYPE, WD, H, SCALE, Z, TEMP, D)
C***BEGIN PROLOGUE  WNLSM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNNLS
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (WNLSM-S, DWNLSM-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     This is a companion subprogram to WNNLS.
C     The documentation for WNNLS has complete usage instructions.
C
C     In addition to the parameters discussed in the prologue to
C     subroutine WNNLS, the following work arrays are used in
C     subroutine WNLSM  (they are passed through the calling
C     sequence from WNNLS for purposes of variable dimensioning).
C     Their contents will in general be of no interest to the user.
C
C         IPIVOT(*)
C            An array of length N.  Upon completion it contains the
C         pivoting information for the cols of W(*,*).
C
C         ITYPE(*)
C            An array of length M which is used to keep track
C         of the classification of the equations.  ITYPE(I)=0
C         denotes equation I as an equality constraint.
C         ITYPE(I)=1 denotes equation I as a least squares
C         equation.
C
C         WD(*)
C            An array of length N.  Upon completion it contains the
C         dual solution vector.
C
C         H(*)
C            An array of length N.  Upon completion it contains the
C         pivot scalars of the Householder transformations performed
C         in the case KRANK.LT.L.
C
C         SCALE(*)
C            An array of length M which is used by the subroutine
C         to store the diagonal matrix of weights.
C         These are used to apply the modified Givens
C         transformations.
C
C         Z(*),TEMP(*)
C            Working arrays of length N.
C
C         D(*)
C            An array of length N that contains the
C         column scaling for the matrix (E).
C                                       (A)
C
C***SEE ALSO  WNNLS
C***ROUTINES CALLED  H12, ISAMAX, R1MACH, SASUM, SAXPY, SCOPY, SNRM2,
C                    SROTM, SROTMG, SSCAL, SSWAP, WNLIT, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and revised.  (WRB & RWC)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C   900510  Fixed an error message.  (RWC)
C***END PROLOGUE  WNLSM
      INTEGER IPIVOT(*), ITYPE(*), L, MA, MDW, MME, MODE, N
      REAL             D(*), H(*), PRGOPT(*), RNORM, SCALE(*), TEMP(*),
     *   W(MDW,*), WD(*), X(*), Z(*)
C
      EXTERNAL H12, ISAMAX, R1MACH, SASUM, SAXPY, SCOPY, SNRM2, SROTM,
     *   SROTMG, SSCAL, SSWAP, WNLIT, XERMSG
      REAL             R1MACH, SASUM, SNRM2
      INTEGER ISAMAX
C
      REAL             ALAMDA, ALPHA, ALSQ, AMAX, BLOWUP, BNORM,
     *   DOPE(3), EANORM, FAC, SM, SPARAM(5), SRELPR, T, TAU, WMAX, Z2,
     *   ZZ
      INTEGER I, IDOPE(3), IMAX, ISOL, ITEMP, ITER, ITMAX, IWMAX, J,
     *   JCON, JP, KEY, KRANK, L1, LAST, LINK, M, ME, NEXT, NIV, NLINK,
     *   NOPT, NSOLN, NTIMES
      LOGICAL DONE, FEASBL, FIRST, HITCON, POS
C
      SAVE SRELPR, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  WNLSM
C
C     Initialize variables.
C     SRELPR is the precision for the particular machine
C     being used.  This logic avoids resetting it every entry.
C
      IF (FIRST) SRELPR = R1MACH(4)
      FIRST = .FALSE.
C
C     Set the nominal tolerance used in the code.
C
      TAU = SQRT(SRELPR)
C
      M = MA + MME
      ME = MME
      MODE = 2
C
C     To process option vector
C
      FAC = 1.E-4
C
C     Set the nominal blow up factor used in the code.
C
      BLOWUP = TAU
C
C     The nominal column scaling used in the code is
C     the identity scaling.
C
      CALL SCOPY (N, 1.E0, 0, D, 1)
C
C     Define bound for number of options to change.
C
      NOPT = 1000
C
C     Define bound for positive value of LINK.
C
      NLINK = 100000
      NTIMES = 0
      LAST = 1
      LINK = PRGOPT(1)
      IF (LINK.LE.0 .OR. LINK.GT.NLINK) THEN
         CALL XERMSG ('SLATEC', 'WNLSM',
     +      'WNNLS, THE OPTION VECTOR IS UNDEFINED', 3, 1)
         RETURN
      ENDIF
C
  100 IF (LINK.GT.1) THEN
         NTIMES = NTIMES + 1
         IF (NTIMES.GT.NOPT) THEN
         CALL XERMSG ('SLATEC', 'WNLSM',
     +      'WNNLS, THE LINKS IN THE OPTION VECTOR ARE CYCLING.', 3, 1)
            RETURN
         ENDIF
C
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.6 .AND. PRGOPT(LAST+2).NE.0.E0) THEN
            DO 110 J = 1,N
               T = SNRM2(M,W(1,J),1)
               IF (T.NE.0.E0) T = 1.E0/T
               D(J) = T
  110       CONTINUE
         ENDIF
C
         IF (KEY.EQ.7) CALL SCOPY (N, PRGOPT(LAST+2), 1, D, 1)
         IF (KEY.EQ.8) TAU = MAX(SRELPR,PRGOPT(LAST+2))
         IF (KEY.EQ.9) BLOWUP = MAX(SRELPR,PRGOPT(LAST+2))
C
         NEXT = PRGOPT(LINK)
         IF (NEXT.LE.0 .OR. NEXT.GT.NLINK) THEN
            CALL XERMSG ('SLATEC', 'WNLSM',
     +         'WNNLS, THE OPTION VECTOR IS UNDEFINED', 3, 1)
            RETURN
         ENDIF
C
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
C
      DO 120 J = 1,N
         CALL SSCAL (M, D(J), W(1,J), 1)
  120 CONTINUE
C
C     Process option vector
C
      DONE = .FALSE.
      ITER = 0
      ITMAX = 3*(N-L)
      MODE = 0
      NSOLN = L
      L1 = MIN(M,L)
C
C     Compute scale factor to apply to equality constraint equations.
C
      DO 130 J = 1,N
         WD(J) = SASUM(M,W(1,J),1)
  130 CONTINUE
C
      IMAX = ISAMAX(N,WD,1)
      EANORM = WD(IMAX)
      BNORM = SASUM(M,W(1,N+1),1)
      ALAMDA = EANORM/(SRELPR*FAC)
C
C     Define scaling diagonal matrix for modified Givens usage and
C     classify equation types.
C
      ALSQ = ALAMDA**2
      DO 140 I = 1,M
C
C        When equation I is heavily weighted ITYPE(I)=0,
C        else ITYPE(I)=1.
C
         IF (I.LE.ME) THEN
            T = ALSQ
            ITEMP = 0
         ELSE
            T = 1.E0
            ITEMP = 1
         ENDIF
         SCALE(I) = T
         ITYPE(I) = ITEMP
  140 CONTINUE
C
C     Set the solution vector X(*) to zero and the column interchange
C     matrix to the identity.
C
      CALL SCOPY (N, 0.E0, 0, X, 1)
      DO 150 I = 1,N
         IPIVOT(I) = I
  150 CONTINUE
C
C     Perform initial triangularization in the submatrix
C     corresponding to the unconstrained variables.
C     Set first L components of dual vector to zero because
C     these correspond to the unconstrained variables.
C
      CALL SCOPY (L, 0.E0, 0, WD, 1)
C
C     The arrays IDOPE(*) and DOPE(*) are used to pass
C     information to WNLIT().  This was done to avoid
C     a long calling sequence or the use of COMMON.
C
      IDOPE(1) = ME
      IDOPE(2) = NSOLN
      IDOPE(3) = L1
C
      DOPE(1) = ALSQ
      DOPE(2) = EANORM
      DOPE(3) = TAU
      CALL WNLIT (W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE, RNORM,
     +            IDOPE, DOPE, DONE)
      ME    = IDOPE(1)
      KRANK = IDOPE(2)
      NIV   = IDOPE(3)
C
C     Perform WNNLS algorithm using the following steps.
C
C     Until(DONE)
C        compute search direction and feasible point
C        when (HITCON) add constraints
C        else perform multiplier test and drop a constraint
C        fin
C     Compute-Final-Solution
C
C     To compute search direction and feasible point,
C     solve the triangular system of currently non-active
C     variables and store the solution in Z(*).
C
C     To solve system
C     Copy right hand side into TEMP vector to use overwriting method.
C
  160 IF (DONE) GO TO 330
      ISOL = L + 1
      IF (NSOLN.GE.ISOL) THEN
         CALL SCOPY (NIV, W(1,N+1), 1, TEMP, 1)
         DO 170 J = NSOLN,ISOL,-1
            IF (J.GT.KRANK) THEN
               I = NIV - NSOLN + J
            ELSE
               I = J
            ENDIF
C
            IF (J.GT.KRANK .AND. J.LE.L) THEN
               Z(J) = 0.E0
            ELSE
               Z(J) = TEMP(I)/W(I,J)
               CALL SAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
            ENDIF
  170    CONTINUE
      ENDIF
C
C     Increment iteration counter and check against maximum number
C     of iterations.
C
      ITER = ITER + 1
      IF (ITER.GT.ITMAX) THEN
         MODE = 1
         DONE = .TRUE.
      ENDIF
C
C     Check to see if any constraints have become active.
C     If so, calculate an interpolation factor so that all
C     active constraints are removed from the basis.
C
      ALPHA = 2.E0
      HITCON = .FALSE.
      DO 180 J = L+1,NSOLN
         ZZ = Z(J)
         IF (ZZ.LE.0.E0) THEN
            T = X(J)/(X(J)-ZZ)
            IF (T.LT.ALPHA) THEN
               ALPHA = T
               JCON = J
            ENDIF
            HITCON = .TRUE.
         ENDIF
  180 CONTINUE
C
C     Compute search direction and feasible point
C
      IF (HITCON) THEN
C
C        To add constraints, use computed ALPHA to interpolate between
C        last feasible solution X(*) and current unconstrained (and
C        infeasible) solution Z(*).
C
         DO 190 J = L+1,NSOLN
            X(J) = X(J) + ALPHA*(Z(J)-X(J))
  190    CONTINUE
         FEASBL = .FALSE.
C
C        Remove column JCON and shift columns JCON+1 through N to the
C        left.  Swap column JCON into the N th position.  This achieves
C        upper Hessenberg form for the nonactive constraints and
C        leaves an upper Hessenberg matrix to retriangularize.
C
  200    DO 210 I = 1,M
            T = W(I,JCON)
            CALL SCOPY (N-JCON, W(I, JCON+1), MDW, W(I, JCON), MDW)
            W(I,N) = T
  210    CONTINUE
C
C        Update permuted index vector to reflect this shift and swap.
C
         ITEMP = IPIVOT(JCON)
         DO 220 I = JCON,N - 1
            IPIVOT(I) = IPIVOT(I+1)
  220    CONTINUE
         IPIVOT(N) = ITEMP
C
C        Similarly permute X(*) vector.
C
         CALL SCOPY (N-JCON, X(JCON+1), 1, X(JCON), 1)
         X(N) = 0.E0
         NSOLN = NSOLN - 1
         NIV = NIV - 1
C
C        Retriangularize upper Hessenberg matrix after adding
C        constraints.
C
         I = KRANK + JCON - L
         DO 230 J = JCON,NSOLN
            IF (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.0) THEN
C
C              Zero IP1 to I in column J
C
               IF (W(I+1,J).NE.0.E0) THEN
                  CALL SROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),
     +                         SPARAM)
                  W(I+1,J) = 0.E0
                  CALL SROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.1) THEN
C
C              Zero IP1 to I in column J
C
               IF (W(I+1,J).NE.0.E0) THEN
                  CALL SROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),
     +                         SPARAM)
                  W(I+1,J) = 0.E0
                  CALL SROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.0) THEN
               CALL SSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
               CALL SSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
               ITEMP = ITYPE(I+1)
               ITYPE(I+1) = ITYPE(I)
               ITYPE(I) = ITEMP
C
C              Swapped row was formerly a pivot element, so it will
C              be large enough to perform elimination.
C              Zero IP1 to I in column J.
C
               IF (W(I+1,J).NE.0.E0) THEN
                  CALL SROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),
     +                         SPARAM)
                  W(I+1,J) = 0.E0
                  CALL SROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.1) THEN
               IF (SCALE(I)*W(I,J)**2/ALSQ.GT.(TAU*EANORM)**2) THEN
C
C                 Zero IP1 to I in column J
C
                  IF (W(I+1,J).NE.0.E0) THEN
                     CALL SROTMG (SCALE(I), SCALE(I+1), W(I,J),
     +                            W(I+1,J), SPARAM)
                     W(I+1,J) = 0.E0
                     CALL SROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                           SPARAM)
                  ENDIF
               ELSE
                  CALL SSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
                  CALL SSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
                  ITEMP = ITYPE(I+1)
                  ITYPE(I+1) = ITYPE(I)
                  ITYPE(I) = ITEMP
                  W(I+1,J) = 0.E0
               ENDIF
            ENDIF
            I = I + 1
  230    CONTINUE
C
C        See if the remaining coefficients in the solution set are
C        feasible.  They should be because of the way ALPHA was
C        determined.  If any are infeasible, it is due to roundoff
C        error.  Any that are non-positive will be set to zero and
C        removed from the solution set.
C
         DO 240 JCON = L+1,NSOLN
            IF (X(JCON).LE.0.E0) GO TO 250
  240    CONTINUE
         FEASBL = .TRUE.
  250    IF (.NOT.FEASBL) GO TO 200
      ELSE
C
C        To perform multiplier test and drop a constraint.
C
         CALL SCOPY (NSOLN, Z, 1, X, 1)
         IF (NSOLN.LT.N) CALL SCOPY (N-NSOLN, 0.E0, 0, X(NSOLN+1), 1)
C
C        Reclassify least squares equations as equalities as necessary.
C
         I = NIV + 1
  260    IF (I.LE.ME) THEN
            IF (ITYPE(I).EQ.0) THEN
               I = I + 1
            ELSE
               CALL SSWAP (N+1, W(I,1), MDW, W(ME,1), MDW)
               CALL SSWAP (1, SCALE(I), 1, SCALE(ME), 1)
               ITEMP = ITYPE(I)
               ITYPE(I) = ITYPE(ME)
               ITYPE(ME) = ITEMP
               ME = ME - 1
            ENDIF
            GO TO 260
         ENDIF
C
C        Form inner product vector WD(*) of dual coefficients.
C
         DO 280 J = NSOLN+1,N
            SM = 0.E0
            DO 270 I = NSOLN+1,M
               SM = SM + SCALE(I)*W(I,J)*W(I,N+1)
  270       CONTINUE
            WD(J) = SM
  280    CONTINUE
C
C        Find J such that WD(J)=WMAX is maximum.  This determines
C        that the incoming column J will reduce the residual vector
C        and be positive.
C
  290    WMAX = 0.E0
         IWMAX = NSOLN + 1
         DO 300 J = NSOLN+1,N
            IF (WD(J).GT.WMAX) THEN
               WMAX = WD(J)
               IWMAX = J
            ENDIF
  300    CONTINUE
         IF (WMAX.LE.0.E0) GO TO 330
C
C        Set dual coefficients to zero for incoming column.
C
         WD(IWMAX) = 0.E0
C
C        WMAX .GT. 0.E0, so okay to move column IWMAX to solution set.
C        Perform transformation to retriangularize, and test for near
C        linear dependence.
C
C        Swap column IWMAX into NSOLN-th position to maintain upper
C        Hessenberg form of adjacent columns, and add new column to
C        triangular decomposition.
C
         NSOLN = NSOLN + 1
         NIV = NIV + 1
         IF (NSOLN.NE.IWMAX) THEN
            CALL SSWAP (M, W(1,NSOLN), 1, W(1,IWMAX), 1)
            WD(IWMAX) = WD(NSOLN)
            WD(NSOLN) = 0.E0
            ITEMP = IPIVOT(NSOLN)
            IPIVOT(NSOLN) = IPIVOT(IWMAX)
            IPIVOT(IWMAX) = ITEMP
         ENDIF
C
C        Reduce column NSOLN so that the matrix of nonactive constraints
C        variables is triangular.
C
         DO 320 J = M,NIV+1,-1
            JP = J - 1
C
C           When operating near the ME line, test to see if the pivot
C           element is near zero.  If so, use the largest element above
C           it as the pivot.  This is to maintain the sharp interface
C           between weighted and non-weighted rows in all cases.
C
            IF (J.EQ.ME+1) THEN
               IMAX = ME
               AMAX = SCALE(ME)*W(ME,NSOLN)**2
               DO 310 JP = J - 1,NIV,-1
                  T = SCALE(JP)*W(JP,NSOLN)**2
                  IF (T.GT.AMAX) THEN
                     IMAX = JP
                     AMAX = T
                  ENDIF
  310          CONTINUE
               JP = IMAX
            ENDIF
C
            IF (W(J,NSOLN).NE.0.E0) THEN
               CALL SROTMG (SCALE(JP), SCALE(J), W(JP,NSOLN),
     +                      W(J,NSOLN), SPARAM)
               W(J,NSOLN) = 0.E0
               CALL SROTM (N+1-NSOLN, W(JP,NSOLN+1), MDW, W(J,NSOLN+1),
     +                     MDW, SPARAM)
            ENDIF
  320    CONTINUE
C
C        Solve for Z(NSOLN)=proposed new value for X(NSOLN).  Test if
C        this is nonpositive or too large.  If this was true or if the
C        pivot term was zero, reject the column as dependent.
C
         IF (W(NIV,NSOLN).NE.0.E0) THEN
            ISOL = NIV
            Z2 = W(ISOL,N+1)/W(ISOL,NSOLN)
            Z(NSOLN) = Z2
            POS = Z2 .GT. 0.E0
            IF (Z2*EANORM.GE.BNORM .AND. POS) THEN
               POS = .NOT. (BLOWUP*Z2*EANORM.GE.BNORM)
            ENDIF
C
C           Try to add row ME+1 as an additional equality constraint.
C           Check size of proposed new solution component.
C           Reject it if it is too large.
C
         ELSEIF (NIV.LE.ME .AND. W(ME+1,NSOLN).NE.0.E0) THEN
            ISOL = ME + 1
            IF (POS) THEN
C
C              Swap rows ME+1 and NIV, and scale factors for these rows.
C
               CALL SSWAP (N+1, W(ME+1,1), MDW, W(NIV,1), MDW)
               CALL SSWAP (1, SCALE(ME+1), 1, SCALE(NIV), 1)
               ITEMP = ITYPE(ME+1)
               ITYPE(ME+1) = ITYPE(NIV)
               ITYPE(NIV) = ITEMP
               ME = ME + 1
            ENDIF
         ELSE
            POS = .FALSE.
         ENDIF
C
         IF (.NOT.POS) THEN
            NSOLN = NSOLN - 1
            NIV = NIV - 1
         ENDIF
         IF (.NOT.(POS.OR.DONE)) GO TO 290
      ENDIF
      GO TO 160
C
C     Else perform multiplier test and drop a constraint.  To compute
C     final solution.  Solve system, store results in X(*).
C
C     Copy right hand side into TEMP vector to use overwriting method.
C
  330 ISOL = 1
      IF (NSOLN.GE.ISOL) THEN
         CALL SCOPY (NIV, W(1,N+1), 1, TEMP, 1)
         DO 340 J = NSOLN,ISOL,-1
            IF (J.GT.KRANK) THEN
               I = NIV - NSOLN + J
            ELSE
               I = J
            ENDIF
C
            IF (J.GT.KRANK .AND. J.LE.L) THEN
               Z(J) = 0.E0
            ELSE
               Z(J) = TEMP(I)/W(I,J)
               CALL SAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
            ENDIF
  340    CONTINUE
      ENDIF
C
C     Solve system.
C
      CALL SCOPY (NSOLN, Z, 1, X, 1)
C
C     Apply Householder transformations to X(*) if KRANK.LT.L
C
      IF (KRANK.LT.L) THEN
         DO 350 I = 1,KRANK
            CALL H12 (2, I, KRANK+1, L, W(I,1), MDW, H(I), X, 1, 1, 1)
  350    CONTINUE
      ENDIF
C
C     Fill in trailing zeroes for constrained variables not in solution.
C
      IF (NSOLN.LT.N) CALL SCOPY (N-NSOLN, 0.E0, 0, X(NSOLN+1), 1)
C
C     Permute solution vector to natural order.
C
      DO 380 I = 1,N
         J = I
  360    IF (IPIVOT(J).EQ.I) GO TO 370
         J = J + 1
         GO TO 360
C
  370    IPIVOT(J) = IPIVOT(I)
         IPIVOT(I) = J
         CALL SSWAP (1, X(J), 1, X(I), 1)
  380 CONTINUE
C
C     Rescale the solution using the column scaling.
C
      DO 390 J = 1,N
         X(J) = X(J)*D(J)
  390 CONTINUE
C
      DO 400 I = NSOLN+1,M
         T = W(I,N+1)
         IF (I.LE.ME) T = T/ALAMDA
         T = (SCALE(I)*T)*T
         RNORM = RNORM + T
  400 CONTINUE
C
      RNORM = SQRT(RNORM)
      RETURN
      END
