*DECK DBVDER
      SUBROUTINE DBVDER (X, Y, YP, G, IPAR)
C***BEGIN PROLOGUE  DBVDER
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (BVDER-S, DBVDER-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C **********************************************************************
C     NFC = Number of base solution vectors
C
C     NCOMP = Number of components per solution vector
C
C              1 -- Nonzero particular solution
C     INHOMO =
C              2 or 3 -- Zero particular solution
C
C             0 -- Inhomogeneous vector term G(X) identically zero
C     IGOFX =
C             1 -- Inhomogeneous vector term G(X) not identically zero
C
C     G = Inhomogeneous vector term G(X)
C
C     XSAV = Previous value of X
C
C     C = Normalization factor for the particular solution
C
C           0   ( if  NEQIVP = 0 )
C     IVP =
C           Number of differential equations integrated due to
C           the original boundary value problem   ( if  NEQIVP .GT. 0 )
C
C     NOFST - For problems with auxiliary initial value equations,
C             NOFST communicates to the routine DFMAT how to access
C             the dependent variables corresponding to this initial
C             value problem.  For example, during any call to DFMAT,
C             the first dependent variable for the initial value
C             problem is in position  Y(NOFST + 1).
C             See example in SAND77-1328.
C **********************************************************************
C
C***SEE ALSO  DBVSUP
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    DML8SZ, DMLIVP
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910701  Corrected ROUTINES CALLED section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C   920618  Minor restructuring of code.  (RWC, WRB)
C***END PROLOGUE  DBVDER
      INTEGER IGOFX, INHOMO, IPAR, IVP, J, K, L, NA, NCOMP, NFC, NOFST
      DOUBLE PRECISION C, G(*), X, XSAV, Y(*), YP(*)
C
C **********************************************************************
C
      COMMON /DML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMP,NFC
C
C **********************************************************************
C     The COMMON block below is used to communicate with the user
C     supplied subroutine DFMAT.  The user should not alter this
C     COMMON block.
C
      COMMON /DMLIVP/ NOFST
C **********************************************************************
C
C***FIRST EXECUTABLE STATEMENT  DBVDER
      IF (IVP .GT. 0) CALL DUIVP(X,Y(IVP+1),YP(IVP+1))
      NOFST = IVP
      NA = 1
      DO 10 K=1,NFC
         CALL DFMAT(X,Y(NA),YP(NA))
         NOFST = NOFST - NCOMP
         NA = NA + NCOMP
   10 CONTINUE
C
      IF (INHOMO .NE. 1) RETURN
      CALL DFMAT(X,Y(NA),YP(NA))
C
      IF (IGOFX .EQ. 0) RETURN
      IF (X .NE. XSAV) THEN
         IF (IVP .EQ. 0) CALL DGVEC(X,G)
         IF (IVP .GT. 0) CALL DUVEC(X,Y(IVP+1),G)
         XSAV = X
      ENDIF
C
C     If the user has chosen not to normalize the particular
C     solution, then C is defined in DBVPOR to be 1.0
C
C     The following loop is just
C     CALL DAXPY (NCOMP, 1.0D0/C, G, 1, YP(NA), 1)
C
      DO 20 J=1,NCOMP
         L = NA + J - 1
         YP(L) = YP(L) + G(J)/C
   20 CONTINUE
      RETURN
      END
