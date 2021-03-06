C
C     file comf.f
C
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C PACKAGE COMF           THE ENTRIES IN THIS PACKAGE ARE LOWLEVEL
C                        ENTRIES, SUPPORTING FISHPACK ENTRIES BLKTRI
C                        AND CBLKTRI. THAT IS, THESE ROUTINES ARE
C                        NOT CALLED DIRECTLY BY USERS, BUT RATHER
C                        BY ENTRIES WITHIN BLKTRI AND CBLKTRI.
C                        DESCRIPTION OF ENTRIES EPMACH AND PIMACH
C                        FOLLOW BELOW.
C
C LATEST REVISION        JUNE 2004
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       NONE
C FILES
C
C LANGUAGE               FORTRAN 90
C ********************************************************************
C
C FUNCTION EPMACH (DUM)
C
C PURPOSE                TO COMPUTE AN APPROXIMATE MACHINE ACCURACY
C                        EPSILON ACCORDING TO THE FOLLOWING DEFINITION:
C                        EPSILON IS THE SMALLEST NUMBER SUCH THAT
C                        (1.+EPSILON).GT.1.)
C
C USAGE                  EPS = EPMACH (DUM)
C
C ARGUMENTS
C ON INPUT               DUM
C                          DUMMY VALUE
C
C ARGUMENTS
C ON OUTPUT              NONE
C
C HISTORY                THE ORIGINAL VERSION, WRITTEN WHEN THE
C                        BLKTRI PACKAGE WAS CONVERTED FROM THE
C                        CDC 7600 TO RUN ON THE CRAY-1, CALCULATED
C                        MACHINE ACCURACY BY SUCCESSIVE DIVISIONS
C                        BY 10.  USE OF THIS CONSTANT CAUSED BLKTRI
C                        TO COMPUTE SOLUTIONS ON THE CRAY-1 WITH FOUR
C                        FEWER PLACES OF ACCURACY THAN THE VERSION
C                        ON THE 7600.  IT WAS FOUND THAT COMPUTING
C                        MACHINE ACCURACY BY SUCCESSIVE DIVISIONS
C                        OF 2 PRODUCED A MACHINE ACCURACY 29% LESS
C                        THAN THE VALUE GENERATED BY SUCCESSIVE
C                        DIVISIONS BY 10, AND THAT USE OF THIS
C                        MACHINE CONSTANT IN THE BLKTRI PACKAGE
C                        RECOVERED THE ACCURACY THAT APPEARED TO
C                        BE LOST ON CONVERSION.
C
C ALGORITHM              COMPUTES MACHINE ACCURACY BY SUCCESSIVE
C                        DIVISIONS OF TWO.
C
C PORTABILITY            THIS CODE WILL EXECUTE ON MACHINES OTHER
C                        THAN THE CRAY1, BUT THE RETURNED VALUE MAY
C                        BE UNSATISFACTORY.  SEE HISTORY ABOVE.
C ********************************************************************
C
C FUNCTION PIMACH (DUM)
C
C PURPOSE                TO SUPPLY THE VALUE OF THE CONSTANT PI
C                        CORRECT TO MACHINE PRECISION WHERE
C                        PI=3.141592653589793238462643383279502884197
C                             1693993751058209749446
C
C USAGE                  PI = PIMACH (DUM)
C
C ARGUMENTS
C ON INPUT               DUM
C                          DUMMY VALUE
C
C ARGUMENTS
C ON OUTPUT              NONE
C
C ALGORITHM              THE VALUE OF PI IS SET TO 4.*ATAN(1.0)
C
C PORTABILITY            THIS ENTRY IS PORTABLE, BUT USERS SHOULD
C                        CHECK TO SEE WHETHER GREATER ACCURACY IS
C                        REQUIRED.
C
C***********************************************************************
      REAL FUNCTION EPMACH (DUM) 
      IMPLICIT NONE
      SAVE ALL
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL  :: DUM 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL :: ALL, EPS,V
      COMMON /VALUE/  V                                                 
      EPS = 1. 
      EPS = EPS/2. 
      CALL STRWRD (EPS + 1.) 
      DO WHILE(V - 1. > 0.) 
         EPS = EPS/2. 
         CALL STRWRD (EPS + 1.) 
      END DO 
      EPMACH = 100.*EPS 
      RETURN  
      END FUNCTION EPMACH 


      SUBROUTINE STRWRD(X) 
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL , INTENT(IN) :: X 
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /VALUE/ 
      COMMON /VALUE/ V 
      REAL   V 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL :: ALL 

      SAVE ALL 
C-----------------------------------------------
      V = X 
      RETURN  
      END SUBROUTINE STRWRD 


      REAL FUNCTION PIMACH (DUM) 
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL  :: DUM 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------
C     PI=3.1415926535897932384626433832795028841971693993751058209749446
C
      PIMACH = 4.*ATAN(1.0) 
      RETURN  
      END FUNCTION PIMACH 


      REAL FUNCTION PPSGF (X, IZ, C, A, BH) 
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IZ 
      REAL , INTENT(IN) :: X 
      REAL  :: C(*) 
      REAL  :: A(*) 
      REAL , INTENT(IN) :: BH(*) 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J 
      REAL :: ALL, SUM 

      SAVE ALL 
C-----------------------------------------------
      SUM = 0. 
      DO J = 1, IZ 
         SUM = SUM - 1./(X - BH(J))**2 
      END DO 
      PPSGF = SUM 
      RETURN  
      END FUNCTION PPSGF 


      REAL FUNCTION PPSPF (X, IZ, C, A, BH) 
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IZ 
      REAL , INTENT(IN) :: X 
      REAL  :: C(*) 
      REAL  :: A(*) 
      REAL , INTENT(IN) :: BH(*) 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J 
      REAL :: ALL, SUM 

      SAVE ALL 
C-----------------------------------------------
      SUM = 0. 
      DO J = 1, IZ 
         SUM = SUM + 1./(X - BH(J)) 
      END DO 
      PPSPF = SUM 
      RETURN  
      END FUNCTION PPSPF 


      REAL FUNCTION PSGF (X, IZ, C, A, BH) 
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IZ 
      REAL , INTENT(IN) :: X 
      REAL , INTENT(IN) :: C(*) 
      REAL , INTENT(IN) :: A(*) 
      REAL , INTENT(IN) :: BH(*) 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J 
      REAL :: ALL, FSG, HSG, DD 

      SAVE ALL 
C-----------------------------------------------
      FSG = 1. 
      HSG = 1. 
      DO J = 1, IZ 
         DD = 1./(X - BH(J)) 
         FSG = FSG*A(J)*DD 
         HSG = HSG*C(J)*DD 
      END DO 
      IF (MOD(IZ,2) == 0) THEN 
         PSGF = 1. - FSG - HSG 
         RETURN  
      ENDIF 
      PSGF = 1. + FSG + HSG 
      RETURN  
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
C June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END FUNCTION PSGF 
