*DECK DY
      SUBROUTINE DY (U, IDMN, I, J, UYYY, UYYYY)
C***BEGIN PROLOGUE  DY
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SEPELI
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (DY-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This program computes second order finite difference
C     approximations to the third and fourth Y
C     partial derivatives of U at the (I,J) mesh point.
C
C***SEE ALSO  SEPELI
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    SPLPCM
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  DY
C
      COMMON /SPLPCM/ KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      DIMENSION       U(IDMN,*)
C***FIRST EXECUTABLE STATEMENT  DY
      IF (J.GT.2 .AND. J.LT.(L-1)) GO TO  50
      IF (J .EQ. 1) GO TO  10
      IF (J .EQ. 2) GO TO  30
      IF (J .EQ. L-1) GO TO  60
      IF (J .EQ. L) GO TO  80
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C
C
   10 IF (KSWY .EQ. 1) GO TO  20
      UYYY = (-5.0*U(I,1)+18.0*U(I,2)-24.0*U(I,3)+14.0*U(I,4)-
     1                                                 3.0*U(I,5))/TDLY3
      UYYYY = (3.0*U(I,1)-14.0*U(I,2)+26.0*U(I,3)-24.0*U(I,4)+
     1                                      11.0*U(I,5)-2.0*U(I,6))/DLY4
      RETURN
C
C     PERIODIC AT X=A
C
   20 UYYY = (-U(I,L-2)+2.0*U(I,L-1)-2.0*U(I,2)+U(I,3))/TDLY3
      UYYYY = (U(I,L-2)-4.0*U(I,L-1)+6.0*U(I,1)-4.0*U(I,2)+U(I,3))/DLY4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C+DLY
C
   30 IF (KSWY .EQ. 1) GO TO  40
      UYYY = (-3.0*U(I,1)+10.0*U(I,2)-12.0*U(I,3)+6.0*U(I,4)-U(I,5))/
     1       TDLY3
      UYYYY = (2.0*U(I,1)-9.0*U(I,2)+16.0*U(I,3)-14.0*U(I,4)+6.0*U(I,5)-
     1                                                      U(I,6))/DLY4
      RETURN
C
C     PERIODIC AT Y=C+DLY
C
   40 UYYY = (-U(I,L-1)+2.0*U(I,1)-2.0*U(I,3)+U(I,4))/TDLY3
      UYYYY = (U(I,L-1)-4.0*U(I,1)+6.0*U(I,2)-4.0*U(I,3)+U(I,4))/DLY4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
C
   50 CONTINUE
      UYYY = (-U(I,J-2)+2.0*U(I,J-1)-2.0*U(I,J+1)+U(I,J+2))/TDLY3
      UYYYY = (U(I,J-2)-4.0*U(I,J-1)+6.0*U(I,J)-4.0*U(I,J+1)+U(I,J+2))/
     1        DLY4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D-DLY
C
   60 IF (KSWY .EQ. 1) GO TO  70
      UYYY = (U(I,L-4)-6.0*U(I,L-3)+12.0*U(I,L-2)-10.0*U(I,L-1)+
     1                                                 3.0*U(I,L))/TDLY3
      UYYYY = (-U(I,L-5)+6.0*U(I,L-4)-14.0*U(I,L-3)+16.0*U(I,L-2)-
     1                                     9.0*U(I,L-1)+2.0*U(I,L))/DLY4
      RETURN
C
C     PERIODIC AT Y=D-DLY
C
   70 CONTINUE
      UYYY = (-U(I,L-3)+2.0*U(I,L-2)-2.0*U(I,1)+U(I,2))/TDLY3
      UYYYY = (U(I,L-3)-4.0*U(I,L-2)+6.0*U(I,L-1)-4.0*U(I,1)+U(I,2))/
     1        DLY4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D
C
   80 UYYY = -(3.0*U(I,L-4)-14.0*U(I,L-3)+24.0*U(I,L-2)-18.0*U(I,L-1)+
     1                                                 5.0*U(I,L))/TDLY3
      UYYYY = (-2.0*U(I,L-5)+11.0*U(I,L-4)-24.0*U(I,L-3)+26.0*U(I,L-2)-
     1                                    14.0*U(I,L-1)+3.0*U(I,L))/DLY4
      RETURN
      END
