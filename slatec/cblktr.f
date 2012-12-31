*DECK CBLKTR
      SUBROUTINE CBLKTR (IFLG, NP, N, AN, BN, CN, MP, M, AM, BM, CM,
     +   IDIMY, Y, IERROR, W)
C***BEGIN PROLOGUE  CBLKTR
C***PURPOSE  Solve a block tridiagonal system of linear equations
C            (usually resulting from the discretization of separable
C            two-dimensional elliptic equations).
C***LIBRARY   SLATEC (FISHPACK)
C***CATEGORY  I2B4B
C***TYPE      COMPLEX (BLKTRI-S, CBLKTR-C)
C***KEYWORDS  ELLIPTIC PDE, FISHPACK, TRIDIAGONAL LINEAR SYSTEM
C***AUTHOR  Adams, J., (NCAR)
C           Swarztrauber, P. N., (NCAR)
C           Sweet, R., (NCAR)
C***DESCRIPTION
C
C     Subroutine CBLKTR is a complex version of subroutine BLKTRI.
C     Both subroutines solve a system of linear equations of the form
C
C          AN(J)*X(I,J-1) + AM(I)*X(I-1,J) + (BN(J)+BM(I))*X(I,J)
C
C          + CN(J)*X(I,J+1) + CM(I)*X(I+1,J) = Y(I,J)
C
C               For I = 1,2,...,M  and  J = 1,2,...,N.
C
C     I+1 and I-1 are evaluated modulo M and J+1 and J-1 modulo N, i.e.,
C
C          X(I,0) = X(I,N),  X(I,N+1) = X(I,1),
C          X(0,J) = X(M,J),  X(M+1,J) = X(1,J).
C
C     These equations usually result from the discretization of
C     separable elliptic equations.  Boundary conditions may be
C     Dirichlet, Neumann, or periodic.
C
C
C     * * * * * * * * * *     On INPUT     * * * * * * * * * *
C
C     IFLG
C       = 0  Initialization only.  Certain quantities that depend on NP,
C            N, AN, BN, and CN are computed and stored in the work
C            array  W.
C       = 1  The quantities that were computed in the initialization are
C            used to obtain the solution X(I,J).
C
C       NOTE   A call with IFLG=0 takes approximately one half the time
C              time as a call with IFLG = 1.  However, the
C              initialization does not have to be repeated unless NP, N,
C              AN, BN, or CN change.
C
C     NP
C       = 0  If AN(1) and CN(N) are not zero, which corresponds to
C            periodic boundary conditions.
C       = 1  If AN(1) and CN(N) are zero.
C
C     N
C       The number of unknowns in the J-direction. N must be greater
C       than 4. The operation count is proportional to MNlog2(N), hence
C       N should be selected less than or equal to M.
C
C     AN,BN,CN
C       Real one-dimensional arrays of length N that specify the
C       coefficients in the linear equations given above.
C
C     MP
C       = 0  If AM(1) and CM(M) are not zero, which corresponds to
C            periodic boundary conditions.
C       = 1  If AM(1) = CM(M) = 0  .
C
C     M
C       The number of unknowns in the I-direction. M must be greater
C       than 4.
C
C     AM,BM,CM
C       Complex one-dimensional arrays of length M that specify the
C       coefficients in the linear equations given above.
C
C     IDIMY
C       The row (or first) dimension of the two-dimensional array Y as
C       it appears in the program calling BLKTRI.  This parameter is
C       used to specify the variable dimension of Y.  IDIMY must be at
C       least M.
C
C     Y
C       A complex two-dimensional array that specifies the values of
C       the right side of the linear system of equations given above.
C       Y must be dimensioned Y(IDIMY,N) with IDIMY .GE. M.
C
C     W
C       A one-dimensional array that must be provided by the user for
C       work space.
C             If NP=1 define K=INT(log2(N))+1 and set L=2**(K+1) then
C                     W must have dimension (K-2)*L+K+5+MAX(2N,12M)
C
C             If NP=0 define K=INT(log2(N-1))+1 and set L=2**(K+1) then
C                     W must have dimension (K-2)*L+K+5+2N+MAX(2N,12M)
C
C       **IMPORTANT** For purposes of checking, the required dimension
C                     of W is computed by BLKTRI and stored in W(1)
C                     in floating point format.
C
C     * * * * * * * * * *     On Output     * * * * * * * * * *
C
C     Y
C       Contains the solution X.
C
C     IERROR
C       An error flag that indicates invalid input parameters.  Except
C       for number zero, a solution is not attempted.
C
C       = 0  No error.
C       = 1  M is less than 5.
C       = 2  N is less than 5.
C       = 3  IDIMY is less than M.
C       = 4  BLKTRI failed while computing results that depend on the
C            coefficient arrays AN, BN, CN.  Check these arrays.
C       = 5  AN(J)*CN(J-1) is less than 0 for some J. Possible reasons
C            for this condition are
C            1. The arrays AN and CN are not correct.
C            2. Too large a grid spacing was used in the discretization
C               of the elliptic equation.
C            3. The linear equations resulted from a partial
C               differential equation which was not elliptic.
C
C     W
C       Contains intermediate values that must not be destroyed if
C       CBLKTR will be called again with IFLG=1.  W(1) contains the
C       number of locations required by W in floating point format.
C
C *Long Description:
C
C     * * * * * * *   Program Specifications    * * * * * * * * * * * *
C
C     Dimension of   AN(N),BN(N),CN(N),AM(M),BM(M),CM(M),Y(IDIMY,N)
C     Arguments      W(see argument list)
C
C     Latest         June 1979
C     Revision
C
C     Required       CBLKTR,CBLKT1,PROC,PROCP,CPROC,CPROCP,CCMPB,INXCA,
C     Subprograms    INXCB,INXCC,CPADD,PGSF,PPGSF,PPPSF,BCRH,TEVLC,
C                    R1MACH
C
C     Special        The algorithm may fail if ABS(BM(I)+BN(J)) is less
C     Conditions     than ABS(AM(I))+ABS(AN(J))+ABS(CM(I))+ABS(CN(J))
C                    for some I and J. The algorithm will also fail if
C                    AN(J)*CN(J-1) is less than zero for some J.
C                    See the description of the output parameter IERROR.
C
C     Common         CCBLK
C     Blocks
C
C     I/O            NONE
C
C     Precision      Single
C
C     Specialist     Paul Swarztrauber
C
C     Language       FORTRAN
C
C     History        CBLKTR is a complex version of BLKTRI (version 3)
C
C     Algorithm      Generalized Cyclic Reduction (see reference below)
C
C     Space
C     Required       CONTROL DATA 7600
C
C     Portability    American National Standards Institute FORTRAN.
C                    The machine accuracy is set using function R1MACH.
C
C     Required       NONE
C     Resident
C     Routines
C
C     References     Swarztrauber,P. and R. SWEET, 'Efficient Fortran
C                    Subprograms for the solution of elliptic equations'
C                    NCAR TN/IA-109, July, 1975, 138 PP.
C
C                    SWARZTRAUBER P. ,'A Direct Method for The Discrete
C                    Solution of Separable Elliptic Equations', SIAM
C                    J. Numer. Anal.,11(1974) PP. 1136-1150.
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran
C                 subprograms for the solution of elliptic equations,
C                 NCAR TN/IA-109, July 1975, 138 pp.
C               P. N. Swarztrauber, A direct method for the discrete
C                 solution of separable elliptic equations, SIAM Journal
C                 on Numerical Analysis 11, (1974), pp. 1136-1150.
C***ROUTINES CALLED  CBLKT1, CCMPB, CPROC, CPROCP, PROC, PROCP
C***COMMON BLOCKS    CCBLK
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CBLKTR
C
      DIMENSION       AN(*)      ,BN(*)      ,CN(*)      ,AM(*)      ,
     1                BM(*)      ,CM(*)      ,Y(IDIMY,*) ,W(*)
      EXTERNAL        PROC       ,PROCP      ,CPROC      ,CPROCP
      COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
      COMPLEX         AM         ,BM         ,CM         ,Y
C***FIRST EXECUTABLE STATEMENT  CBLKTR
      NM = N
      M2 = M+M
      IERROR = 0
      IF (M-5) 101,102,102
  101 IERROR = 1
      GO TO 119
  102 IF (NM-3) 103,104,104
  103 IERROR = 2
      GO TO 119
  104 IF (IDIMY-M) 105,106,106
  105 IERROR = 3
      GO TO 119
  106 NH = N
      NPP = NP
      IF (NPP) 107,108,107
  107 NH = NH+1
  108 IK = 2
      K = 1
  109 IK = IK+IK
      K = K+1
      IF (NH-IK) 110,110,109
  110 NL = IK
      IK = IK+IK
      NL = NL-1
      IWAH = (K-2)*IK+K+6
      IF (NPP) 111,112,111
C
C     DIVIDE W INTO WORKING SUB ARRAYS
C
  111 IW1 = IWAH
      IWBH = IW1+NM
      W(1) = IW1-1+MAX(2*NM,12*M)
      GO TO 113
  112 IWBH = IWAH+NM+NM
      IW1 = IWBH
      W(1) = IW1-1+MAX(2*NM,12*M)
      NM = NM-1
C
C SUBROUTINE CCMPB COMPUTES THE ROOTS OF THE B POLYNOMIALS
C
  113 IF (IERROR) 119,114,119
  114 IW2 = IW1+M2
      IW3 = IW2+M2
      IWD = IW3+M2
      IWW = IWD+M2
      IWU = IWW+M2
      IF (IFLG) 116,115,116
  115 CALL CCMPB (NL,IERROR,AN,BN,CN,W(2),W(IWAH),W(IWBH))
      GO TO 119
  116 IF (MP) 117,118,117
C
C SUBROUTINE CBLKT1 SOLVES THE LINEAR SYSTEM
C
  117 CALL CBLKT1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),W(IW1),W(IW2),
     1             W(IW3),W(IWD),W(IWW),W(IWU),PROC,CPROC)
      GO TO 119
  118 CALL CBLKT1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),W(IW1),W(IW2),
     1             W(IW3),W(IWD),W(IWW),W(IWU),PROCP,CPROCP)
  119 CONTINUE
      RETURN
      END
