*DECK BVSUP
      SUBROUTINE BVSUP (Y, NROWY, NCOMP, XPTS, NXPTS, A, NROWA, ALPHA,
     +   NIC, B, NROWB, BETA, NFC, IGOFX, RE, AE, IFLAG, WORK, NDW,
     +   IWORK, NDIW, NEQIVP)
C***BEGIN PROLOGUE  BVSUP
C***PURPOSE  Solve a linear two-point boundary value problem using
C            superposition coupled with an orthonormalization procedure
C            and a variable-step integration scheme.
C***LIBRARY   SLATEC
C***CATEGORY  I1B1
C***TYPE      SINGLE PRECISION (BVSUP-S, DBVSUP-D)
C***KEYWORDS  ORTHONORMALIZATION, SHOOTING,
C             TWO-POINT BOUNDARY VALUE PROBLEM
C***AUTHOR  Scott, M. R., (SNLA)
C           Watts, H. A., (SNLA)
C***DESCRIPTION
C
C **********************************************************************
C     Subroutine BVSUP solves a LINEAR two-point boundary-value problem
C     of the form
C                        dY/dX = MATRIX(X,U)*Y(X) + G(X,U)
C                A*Y(Xinitial) = ALPHA ,  B*Y(Xfinal) = BETA
C
C     Coupled with the solution of the initial value problem
C
C                        dU/dX = F(X,U)
C                      U(Xinitial) = ETA
C
C **********************************************************************
C     Abstract
C        The method of solution uses superposition coupled with an
C     orthonormalization procedure and a variable-step integration
C     scheme.  Each time the superposition solutions start to
C     lose their numerical linear independence, the vectors are
C     reorthonormalized before integration proceeds.  The underlying
C     principle of the algorithm is then to piece together the
C     intermediate (orthogonalized) solutions, defined on the various
C     subintervals, to obtain the desired solutions.
C
C **********************************************************************
C     INPUT to BVSUP
C **********************************************************************
C
C     NROWY = Actual row dimension of Y in calling program.
C             NROWY must be .GE. NCOMP
C
C     NCOMP = Number of components per solution vector.
C             NCOMP is equal to number of original differential
C             equations.  NCOMP = NIC + NFC.
C
C     XPTS = Desired output points for solution. They must be monotonic.
C            Xinitial = XPTS(1)
C            Xfinal = XPTS(NXPTS)
C
C     NXPTS = Number of output points
C
C     A(NROWA,NCOMP) = Boundary condition matrix at Xinitial,
C                      must be contained in (NIC,NCOMP) sub-matrix.
C
C     NROWA = Actual row dimension of A in calling program,
C             NROWA must be .GE. NIC.
C
C     ALPHA(NIC+NEQIVP) = Boundary conditions at Xinitial.
C                         If NEQIVP .GT. 0 (see below), the boundary
C                         conditions at Xinitial for the initial value
C                         equations must be stored starting in
C                         position (NIC + 1) of ALPHA.
C                         Thus,  ALPHA(NIC+K) = ETA(K).
C
C     NIC = Number of boundary conditions at Xinitial.
C
C     B(NROWB,NCOMP) = Boundary condition matrix at Xfinal,
C                      must be contained in (NFC,NCOMP) sub-matrix.
C
C     NROWB = Actual row dimension of B in calling program,
C             NROWB must be .GE. NFC.
C
C     BETA(NFC) = Boundary conditions at Xfinal.
C
C     NFC = Number of boundary conditions at Xfinal
C
C     IGOFX =0 -- The inhomogeneous term G(X) is identically zero.
C           =1 -- The inhomogeneous term G(X) is not identically zero.
C                 (if IGOFX=1, then subroutine GVEC (or UVEC) must be
C                  supplied).
C
C     RE = Relative error tolerance used by the integrator
C          (see one of the integrators)
C
C     AE = Absolute error tolerance used by the integrator
C          (see one of the integrators)
C **NOTE-  RE and AE should not both be zero.
C
C     IFLAG = A status parameter used principally for output.
C             However, for efficient solution of problems which
C             are originally defined as complex valued (but
C             converted to real systems to use this code), the
C             user must set IFLAG=13 on input. See the comment below
C             for more information on solving such problems.
C
C     WORK(NDW) = Floating point array used for internal storage.
C
C     NDW = Actual dimension of WORK array allocated by user.
C           An estimate for NDW can be computed from the following
C            NDW = 130 + NCOMP**2 * (6 + NXPTS/2 + expected number of
C                                                orthonormalizations/8)
C             For the DISK or TAPE storage mode,
C            NDW = 6 * NCOMP**2 + 10 * NCOMP + 130
C  However, when the ADAMS integrator is to be used, the estimates are
C            NDW = 130 + NCOMP**2 * (13 + NXPTS/2 + expected number of
C                                                orthonormalizations/8)
C    and     NDW = 13 * NCOMP**2 + 22 * NCOMP + 130   , respectively.
C
C     IWORK(NDIW) = Integer array used for internal storage.
C
C     NDIW = Actual dimension of IWORK array allocated by user.
C            An estimate for NDIW can be computed from the following
C            NDIW = 68 + NCOMP * (1 + expected number of
C                                        orthonormalizations)
C **NOTE --  The amount of storage required is problem dependent and may
C            be difficult to predict in advance. Experience has shown
C            that for most problems 20 or fewer orthonormalizations
C            should suffice. If the problem cannot be completed with the
C            allotted storage, then a message will be printed which
C            estimates the amount of storage necessary. In any case, the
C            user can examine the IWORK array for the actual storage
C            requirements, as described in the output information below.
C
C     NEQIVP = Number of auxiliary initial value equations being added
C              to the boundary value problem.
C **NOTE -- Occasionally the coefficients  MATRIX  and/or  G  may be
C           functions which depend on the independent variable  X  and
C           on  U, the solution of an auxiliary initial value problem.
C           In order to avoid the difficulties associated with
C           interpolation, the auxiliary equations may be solved
C           simultaneously with the given boundary value problem.
C           This initial value problem may be LINEAR or NONLINEAR.
C                 See SAND77-1328 for an example.
C
C
C     The user must supply subroutines FMAT, GVEC, UIVP and UVEC, when
C     needed (they MUST be so named), to evaluate the derivatives
C     as follows
C
C        A. FMAT must be supplied.
C
C              SUBROUTINE FMAT(X,Y,YP)
C              X = Independent variable (input to FMAT)
C              Y = Dependent variable vector (input to FMAT)
C              YP = dY/dX = Derivative vector (output from FMAT)
C
C            Compute the derivatives for the HOMOGENEOUS problem
C              YP(I) = dY(I)/dX = MATRIX(X) * Y(I)  , I = 1,...,NCOMP
C
C            When (NEQIVP .GT. 0) and  MATRIX  is dependent on  U  as
C            well as on  X, the following common statement must be
C            included in FMAT
C                    COMMON /MLIVP/ NOFST
C            For convenience, the  U  vector is stored at the bottom
C            of the  Y  array.  Thus, during any call to FMAT,
C            U(I) is referenced by  Y(NOFST + I).
C
C
C            Subroutine BVDER calls FMAT NFC times to evaluate the
C            homogeneous equations and, if necessary, it calls FMAT once
C            in evaluating the particular solution. Since X remains
C            unchanged in this sequence of calls it is possible to
C            realize considerable computational savings for complicated
C            and expensive evaluations of the MATRIX entries. To do this
C            the user merely passes a variable, say XS, via COMMON where
C            XS is defined in the main program to be any value except
C            the initial X. Then the non-constant elements of MATRIX(X)
C            appearing in the differential equations need only be
C            computed if X is unequal to XS, whereupon XS is reset to X.
C
C
C        B. If  NEQIVP .GT. 0 ,  UIVP must also be supplied.
C
C              SUBROUTINE UIVP(X,U,UP)
C              X = Independent variable (input to UIVP)
C              U = Dependent variable vector (input to UIVP)
C              UP = dU/dX = Derivative vector (output from UIVP)
C
C            Compute the derivatives for the auxiliary initial value eqs
C              UP(I) = dU(I)/dX, I = 1,...,NEQIVP.
C
C            Subroutine BVDER calls UIVP once to evaluate the
C            derivatives for the auxiliary initial value equations.
C
C
C        C. If  NEQIVP = 0  and  IGOFX = 1 ,  GVEC must be supplied.
C
C              SUBROUTINE GVEC(X,G)
C              X = Independent variable (input to GVEC)
C              G = Vector of inhomogeneous terms G(X) (output from GVEC)
C
C            Compute the inhomogeneous terms G(X)
C                G(I) = G(X) values for I = 1,...,NCOMP.
C
C            Subroutine BVDER calls GVEC in evaluating the particular
C            solution provided G(X) is NOT identically zero. Thus, when
C            IGOFX=0, the user need NOT write a GVEC subroutine. Also,
C            the user does not have to bother with the computational
C            savings scheme for GVEC as this is automatically achieved
C            via the BVDER subroutine.
C
C
C        D. If  NEQIVP .GT. 0  and  IGOFX = 1 ,  UVEC must be supplied.
C
C              SUBROUTINE UVEC(X,U,G)
C              X = Independent variable (input to UVEC)
C              U = Dependent variable vector from the auxiliary initial
C                  value problem    (input to UVEC)
C              G = Array of inhomogeneous terms G(X,U)(output from UVEC)
C
C            Compute the inhomogeneous terms G(X,U)
C                G(I) = G(X,U) values for I = 1,...,NCOMP.
C
C            Subroutine BVDER calls UVEC in evaluating the particular
C            solution provided G(X,U) is NOT identically zero.  Thus,
C            when IGOFX=0, the user need NOT write a UVEC subroutine.
C
C
C
C     The following is optional input to BVSUP to give the user more
C     flexibility in use of the code.  See SAND75-0198 , SAND77-1328 ,
C     SAND77-1690,SAND78-0522, and SAND78-1501 for more information.
C
C ****CAUTION -- The user MUST zero out IWORK(1),...,IWORK(15)
C                prior to calling BVSUP. These locations define optional
C                input and MUST be zero UNLESS set to special values by
C                the user as described below.
C
C     IWORK(1) -- Number of orthonormalization points.
C                 A value need be set only if IWORK(11) = 1
C
C     IWORK(9) -- Integrator and orthonormalization parameter
C                 (default value is 1)
C                 1 = RUNGE-KUTTA-FEHLBERG code using GRAM-SCHMIDT test.
C                 2 = ADAMS code using GRAM-SCHMIDT TEST.
C
C     IWORK(11) -- Orthonormalization points parameter
C                  (default value is 0)
C                  0 - Orthonormalization points not pre-assigned.
C                  1 - Orthonormalization points pre-assigned in
C                      the first IWORK(1) positions of WORK.
C
C     IWORK(12) -- Storage parameter
C                  (default value is 0)
C                  0 - All storage IN CORE
C                LUN - Homogeneous and inhomogeneous solutions at
C                     output points and orthonormalization information
C                     are stored on DISK.  The logical unit number to be
C                     used for DISK I/O (NTAPE) is set to IWORK(12).
C
C     WORK(1),... -- Pre-assigned orthonormalization points, stored
C                    monotonically, corresponding to the direction
C                    of integration.
C
C
C
C                 ******************************
C                 *** COMPLEX VALUED PROBLEM ***
C                 ******************************
C **NOTE***
C       Suppose the original boundary value problem is NC equations
C     of the form
C                   dW/dX = MAT(X,U)*W(X) + H(X,U)
C                 R*W(Xinitial)=GAMMA , S*W(Xfinal)=DELTA
C
C     where all variables are complex valued. The BVSUP code can be
C     used by converting to a real system of size 2*NC. To solve the
C     larger dimensioned problem efficiently,  the user must initialize
C     IFLAG=13 on input and order the vector components according to
C     Y(1)=real(W(1)),...,Y(NC)=real(W(NC)),Y(NC+1)=imag(W(1)),....,
C     Y(2*NC)=imag(W(NC)). Then define
C                        ...........................
C                        . real(MAT)    -imag(MAT) .
C            MATRIX  =   .                         .
C                        . imag(MAT)     real(MAT) .
C                        ...........................
C
C     The matrices A,B and vectors G,ALPHA,BETA must be defined
C     similarly. Further details can be found in SAND78-1501.
C
C
C **********************************************************************
C     OUTPUT from BVSUP
C **********************************************************************
C
C     Y(NROWY,NXPTS) = Solution at specified output points.
C
C     IFLAG output values
C            =-5 Algorithm ,for obtaining starting vectors for the
C                special complex problem structure, was unable to obtain
C                the initial vectors satisfying the necessary
C                independence criteria.
C            =-4 Rank of boundary condition matrix A is less than NIC,
C                as determined by LSSUDS.
C            =-2 Invalid input parameters.
C            =-1 Insufficient number of storage locations allocated for
C                WORK or IWORK.
C
C            =0 Indicates successful solution
C
C            =1 A computed solution is returned but UNIQUENESS of the
C               solution of the boundary-value problem is questionable.
C               For an eigenvalue problem, this should be treated as a
C               successful execution since this is the expected mode
C               of return.
C            =2 A computed solution is returned but the EXISTENCE of the
C               solution to the boundary-value problem is questionable.
C            =3 A nontrivial solution approximation is returned although
C               the boundary condition matrix B*Y(Xfinal) is found to be
C               nonsingular (to the desired accuracy level) while the
C               right hand side vector is zero. To eliminate this type
C               of return, the accuracy of the eigenvalue parameter
C               must be improved.
C           ***NOTE- We attempt to diagnose the correct problem behavior
C               and report possible difficulties by the appropriate
C               error flag.  However, the user should probably resolve
C               the problem using smaller error tolerances and/or
C               perturbations in the boundary conditions or other
C               parameters. This will often reveal the correct
C               interpretation for the problem posed.
C
C            =13 Maximum number of orthonormalizations attained before
C                reaching Xfinal.
C            =20-flag from integrator (DERKF or DEABM) values can range
C                from 21 to 25.
C            =30 Solution vectors form a dependent set.
C
C     WORK(1),...,WORK(IWORK(1)) = Orthonormalization points
C                                  determined by BVPOR.
C
C     IWORK(1) = Number of orthonormalizations performed by BVPOR.
C
C     IWORK(2) = Maximum number of orthonormalizations allowed as
C                calculated from storage allocated by user.
C
C     IWORK(3),IWORK(4),IWORK(5),IWORK(6)   Give information about
C                actual storage requirements for WORK and IWORK
C                arrays.  In particular,
C                       required storage for  WORK array is
C        IWORK(3) + IWORK(4)*(expected number of orthonormalizations)
C
C                       required storage for IWORK array is
C        IWORK(5) + IWORK(6)*(expected number of orthonormalizations)
C
C     IWORK(8) = Final value of exponent parameter used in tolerance
C                test for orthonormalization.
C
C     IWORK(16) = Number of independent vectors returned from MGSBV.
C                 It is only of interest when IFLAG=30 is obtained.
C
C     IWORK(17) = Numerically estimated rank of the boundary
C                 condition matrix defined from B*Y(Xfinal)
C
C **********************************************************************
C
C     Necessary machine constants are defined in the function
C     routine R1MACH. The user must make sure that the values
C     set in R1MACH are relevant to the computer being used.
C
C **********************************************************************
C
C***REFERENCES  M. R. Scott and H. A. Watts, SUPORT - a computer code
C                 for two-point boundary-value problems via
C                 orthonormalization, SIAM Journal of Numerical
C                 Analysis 14, (1977), pp. 40-70.
C               B. L. Darlow, M. R. Scott and H. A. Watts, Modifications
C                 of SUPORT, a linear boundary value problem solver
C                 Part I - pre-assigning orthonormalization points,
C                 auxiliary initial value problem, disk or tape storage,
C                 Report SAND77-1328, Sandia Laboratories, Albuquerque,
C                 New Mexico, 1977.
C               B. L. Darlow, M. R. Scott and H. A. Watts, Modifications
C                 of SUPORT, a linear boundary value problem solver
C                 Part II - inclusion of an Adams integrator, Report
C                 SAND77-1690, Sandia Laboratories, Albuquerque,
C                 New Mexico, 1977.
C               M. E. Lord and H. A. Watts, Modifications of SUPORT,
C                 a linear boundary value problem solver Part III -
C                 orthonormalization improvements, Report SAND78-0522,
C                 Sandia Laboratories, Albuquerque, New Mexico, 1978.
C               H. A. Watts, M. R. Scott and M. E. Lord, Computational
C                 solution of complex*16 valued boundary problems,
C                 Report SAND78-1501, Sandia Laboratories,
C                 Albuquerque, New Mexico, 1978.
C***ROUTINES CALLED  EXBVP, MACON, XERMSG
C***COMMON BLOCKS    ML15TO, ML17BW, ML18JR, ML5MCO, ML8SZ
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   890921  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BVSUP
C **********************************************************************
C
C
      DIMENSION Y(NROWY,*),A(NROWA,*),ALPHA(*),B(NROWB,*),
     1          BETA(*),WORK(*),IWORK(*),XPTS(*)
      CHARACTER*8 XERN1, XERN2, XERN3, XERN4
C
C **********************************************************************
C     THE COMMON BLOCK BELOW IS USED TO COMMUNICATE WITH SUBROUTINE
C     BVDER.  THE USER SHOULD NOT ALTER OR USE THIS COMMON BLOCK IN THE
C     CALLING PROGRAM.
C
      COMMON /ML8SZ/ C,XSAV,IGOFXD,INHOMO,IVP,NCOMPD,NFCD
C
C **********************************************************************
C     THESE COMMON BLOCKS AID IN REDUCING THE NUMBER OF SUBROUTINE
C     ARGUMENTS PREVALENT IN THIS MODULAR STRUCTURE
C
      COMMON /ML18JR/ AED,RED,TOL,NXPTSD,NICD,NOPG,MXNON,NDISK,NTAPE,
     1                NEQ,INDPVT,INTEG,NPS,NTP,NEQIVD,NUMORT,NFCC,
     2                ICOCO
      COMMON /ML17BW/ KKKZPW,NEEDW,NEEDIW,K1,K2,K3,K4,K5,K6,K7,K8,K9,
     1                K10,K11,L1,L2,KKKINT,LLLINT
C
C **********************************************************************
C     THIS COMMON BLOCK IS USED IN SUBROUTINES BVSUP,BVPOR,RKFAB,
C     REORT, AND STWAY. IT CONTAINS INFORMATION NECESSARY
C     FOR THE ORTHONORMALIZATION TESTING PROCEDURE AND A BACKUP
C     RESTARTING CAPABILITY.
C
      COMMON /ML15TO/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP,
     1                KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
C
C **********************************************************************
C     THIS COMMON BLOCK CONTAINS THE MACHINE DEPENDENT PARAMETERS
C     USED BY THE CODE
C
      COMMON /ML5MCO/ URO,SRU,EPS,SQOVFL,TWOU,FOURU,LPAR
C
C **********************************************************************
C     SET UP MACHINE DEPENDENT CONSTANTS.
C
C***FIRST EXECUTABLE STATEMENT  BVSUP
      CALL MACON
C
C **********************************************************************
C     TEST FOR INVALID INPUT
C
      IF (NROWY .LT. NCOMP)  GO TO 20
      IF (NCOMP .NE. NIC+NFC)  GO TO 20
      IF (NXPTS .LT. 2)  GO TO 20
      IF (NIC .LE. 0)  GO TO 20
      IF (NROWA .LT. NIC)  GO TO 20
      IF (NFC .LE. 0)  GO TO 20
      IF (NROWB .LT. NFC)  GO TO 20
      IF (IGOFX .LT. 0  .OR.  IGOFX .GT. 1) GO TO 20
      IF (RE .LT. 0.0)  GO TO 20
      IF (AE .LT. 0.0)  GO TO 20
      IF (RE .EQ. 0.0  .AND.  AE .EQ. 0.0)  GO TO 20
      IS = 1
      IF (XPTS(NXPTS) .LT. XPTS(1))  IS = 2
      NXPTSM = NXPTS - 1
      DO 13 K = 1,NXPTSM
      IF (IS .EQ. 2) GO TO 12
      IF (XPTS(K+1) .LE. XPTS(K))  GO TO 20
      GO TO 13
   12 IF (XPTS(K) .LE. XPTS(K+1))  GO TO 20
   13 CONTINUE
      GO TO 30
   20 IFLAG = -2
      RETURN
   30 CONTINUE
C
C **********************************************************************
C     CHECK FOR DISK STORAGE
C
      KPTS = NXPTS
      NDISK = 0
      IF (IWORK(12) .EQ. 0)  GO TO 35
      NTAPE = IWORK(12)
      KPTS = 1
      NDISK = 1
   35 CONTINUE
C
C **********************************************************************
C     SET INTEG PARAMETER ACCORDING TO CHOICE OF INTEGRATOR.
C
      INTEG = 1
      IF (IWORK(9) .EQ. 2)  INTEG = 2
C
C **********************************************************************
C     COMPUTE INHOMO
C
      IF (IGOFX .EQ. 1)  GO TO 43
      DO 40 J = 1,NIC
      IF (ALPHA(J) .NE. 0.0)  GO TO 43
   40 CONTINUE
      DO 41 J = 1,NFC
      IF (BETA(J) .NE. 0.0)  GO TO 42
   41 CONTINUE
      INHOMO = 3
      GO TO 45
   42 INHOMO = 2
      GO TO 45
   43 INHOMO = 1
   45 CONTINUE
C
C **********************************************************************
C     TO TAKE ADVANTAGE OF THE SPECIAL STRUCTURE WHEN SOLVING A
C     COMPLEX VALUED PROBLEM,WE INTRODUCE NFCC=NFC WHILE CHANGING
C     THE INTERNAL VALUE OF NFC
C
      NFCC=NFC
      IF (IFLAG .EQ. 13) NFC=NFC/2
C
C **********************************************************************
C     DETERMINE NECESSARY STORAGE REQUIREMENTS
C
C FOR BASIC ARRAYS IN BVPOR
      KKKYHP = NCOMP*(NFC+1) + NEQIVP
      KKKU   = NCOMP*NFC*KPTS
      KKKV   = NCOMP*KPTS
      KKKCOE = NFCC
      KKKS   = NFC+1
      KKKSTO = NCOMP*(NFC+1) + NEQIVP + 1
      KKKG   = NCOMP
C
C FOR ORTHONORMALIZATION RELATED MATTERS
      NTP = (NFCC*(NFCC+1))/2
      KKKZPW = 1 + NTP + NFCC
      LLLIP  = NFCC
C
C FOR ADDITIONAL REQUIRED WORK SPACE
C   (LSSUDS)
      KKKSUD = 4*NIC + (NROWA+1)*NCOMP
      LLLSUD = NIC
C   (SVECS)
      KKKSVC = 1 + 4*NFCC + 2*NFCC**2
      LLLSVC = 2*NFCC
C
      NDEQ=NCOMP*NFC+NEQIVP
      IF (INHOMO .EQ. 1) NDEQ=NDEQ+NCOMP
      GO TO (51,52),INTEG
C   (DERKF)
   51 KKKINT = 33 + 7*NDEQ
      LLLINT = 34
      GO TO 55
C   (DEABM)
   52 KKKINT = 130 + 21*NDEQ
      LLLINT = 51
C
C   (COEF)
   55 KKKCOF = 5*NFCC + NFCC**2
      LLLCOF = 3 + NFCC
C
      KKKWS  = MAX(KKKSUD,KKKSVC,KKKINT,KKKCOF)
      LLLIWS = MAX(LLLSUD,LLLSVC,LLLINT,LLLCOF)
C
      NEEDW  = KKKYHP + KKKU + KKKV + KKKCOE + KKKS + KKKSTO + KKKG +
     1         KKKZPW + KKKWS
      NEEDIW = 17 + LLLIP + LLLIWS
C **********************************************************************
C     COMPUTE THE NUMBER OF POSSIBLE ORTHONORMALIZATIONS WITH THE
C     ALLOTTED STORAGE
C
      IWORK(3) = NEEDW
      IWORK(4) = KKKZPW
      IWORK(5) = NEEDIW
      IWORK(6) = LLLIP
      NRTEMP = NDW - NEEDW
      NITEMP = NDIW - NEEDIW
      IF (NRTEMP .LT. 0)  GO TO 70
      IF (NITEMP .GE. 0)  GO TO 75
C
   70 IFLAG = -1
      IF (NDISK .NE. 1) THEN
         WRITE (XERN1, '(I8)') NEEDW
         WRITE (XERN2, '(I8)') KKKZPW
         WRITE (XERN3, '(I8)') NEEDIW
         WRITE (XERN4, '(I8)') LLLIP
         CALL XERMSG ('SLATEC', 'BVSUP',
     *      'REQUIRED STORAGE FOR WORK ARRAY IS '  // XERN1 // ' + ' //
     *      XERN2 // '*(EXPECTED NUMBER OF ORTHONORMALIZATIONS) $$'  //
     *      'REQUIRED STORAGE FOR IWORK ARRAY IS ' // XERN3 // ' + ' //
     *      XERN4 // '*(EXPECTED NUMBER OF ORTHONORMALIZATIONS)', 1, 0)
      ELSE
         WRITE (XERN1, '(I8)') NEEDW
         WRITE (XERN2, '(I8)') NEEDIW
         CALL XERMSG ('SLATEC', 'BVSUP',
     *      'REQUIRED STORAGE FOR WORK ARRAY IS '  // XERN1 //
     *      ' + NUMBER OF ORTHONOMALIZATIONS. $$'  //
     *      'REQUIRED STORAGE FOR IWORK ARRAY IS ' // XERN2, 1, 0)
      ENDIF
      RETURN
C
   75 IF (NDISK .EQ. 0)  GO TO 77
      NON = 0
      MXNON = NRTEMP
      GO TO 78
C
   77 MXNONR = NRTEMP / KKKZPW
      MXNONI = NITEMP / LLLIP
      MXNON = MIN(MXNONR,MXNONI)
      NON = MXNON
C
   78 IWORK(2) = MXNON
C
C **********************************************************************
C     CHECK FOR PRE-ASSIGNED ORTHONORMALIZATION POINTS
C
      NOPG = 0
      IF (IWORK(11) .NE. 1)  GO TO 85
      IF (MXNON .LT. IWORK(1))  GO TO 70
      NOPG = 1
      MXNON = IWORK(1)
      WORK(MXNON+1) = 2. * XPTS(NXPTS)  -  XPTS(1)
   85 CONTINUE
C
C **********************************************************************
C     ALLOCATE STORAGE FROM WORK AND IWORK ARRAYS
C
C  (Z)
      K1 = 1 + (MXNON+1)
C  (P)
      K2 = K1 + NTP*(NON+1)
C  (W)
      K3 = K2 + NFCC*(NON+1)
C  (YHP)
      K4 = K3 + KKKYHP
C  (U)
      K5 = K4 + KKKU
C  (V)
      K6 = K5 + KKKV
C  (COEF)
      K7 = K6 + KKKCOE
C  (S)
      K8 = K7 + KKKS
C  (STOWA)
      K9 = K8 + KKKSTO
C  (G)
      K10 = K9 + KKKG
      K11 = K10 + KKKWS
C            REQUIRED ADDITIONAL REAL WORK SPACE STARTS AT WORK(K10)
C            AND EXTENDS TO WORK(K11-1)
C
C     FIRST 17 LOCATIONS OF IWORK ARE USED FOR OPTIONAL
C     INPUT AND OUTPUT ITEMS
C  (IP)
      L1 = 18 + NFCC*(NON+1)
      L2 = L1 + LLLIWS
C            REQUIRED INTEGER WORK SPACE STARTS AT IWORK(L1)
C            AND EXTENDS TO IWORK(L2-1)
C
C **********************************************************************
C     SET INDICATOR FOR NORMALIZATION OF PARTICULAR SOLUTION
C
      NPS = 0
      IF (IWORK(10) .EQ. 1)  NPS = 1
C
C **********************************************************************
C     SET PIVOTING PARAMETER
C
      INDPVT=0
      IF (IWORK(15) .EQ. 1) INDPVT=1
C
C **********************************************************************
C     SET OTHER COMMON BLOCK PARAMETERS
C
      NFCD = NFC
      NCOMPD = NCOMP
      IGOFXD = IGOFX
      NXPTSD = NXPTS
      NICD = NIC
      RED = RE
      AED = AE
      NEQIVD = NEQIVP
      MNSWOT = 20
      IF (IWORK(13) .EQ. -1) MNSWOT=MAX(1,IWORK(14))
      XBEG=XPTS(1)
      XEND=XPTS(NXPTS)
      XSAV=XEND
      ICOCO=1
      IF (INHOMO .EQ. 3  .AND.  NOPG .EQ. 1) WORK(MXNON+1)=XEND
C
C **********************************************************************
C
      CALL EXBVP(Y,NROWY,XPTS,A,NROWA,ALPHA,B,NROWB,BETA,IFLAG,WORK,
     1           IWORK)
      NFC=NFCC
      IWORK(17)=IWORK(L1)
      RETURN
      END
